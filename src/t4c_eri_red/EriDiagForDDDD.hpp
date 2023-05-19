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
#include "EriDiagVRRForSSSS.hpp"
#include "EriDiagVRRForSSSP.hpp"
#include "EriDiagVRRForSSSD.hpp"
#include "EriDiagVRRForSSSF.hpp"
#include "EriDiagVRRForSSSG.hpp"
#include "EriDiagVRRForSPSS.hpp"
#include "EriDiagVRRForSPSP.hpp"
#include "EriDiagVRRForSPSD.hpp"
#include "EriDiagVRRForSPSF.hpp"
#include "EriDiagVRRForSPSG.hpp"
#include "EriDiagVRRForSDSS.hpp"
#include "EriDiagVRRForSDSP.hpp"
#include "EriDiagVRRForSDSD.hpp"
#include "EriDiagVRRForSDSF.hpp"
#include "EriDiagVRRForSDSG.hpp"
#include "EriDiagVRRForSFSP.hpp"
#include "EriDiagVRRForSFSD.hpp"
#include "EriDiagVRRForSFSF.hpp"
#include "EriDiagVRRForSFSG.hpp"
#include "EriDiagVRRForSGSD.hpp"
#include "EriDiagVRRForSGSF.hpp"
#include "EriDiagVRRForSGSG.hpp"
#include "EriDiagHRRForSDPD.hpp"
#include "EriDiagHRRForSDPF.hpp"
#include "EriDiagHRRForSDDD.hpp"
#include "EriDiagHRRForSFPD.hpp"
#include "EriDiagHRRForSFPF.hpp"
#include "EriDiagHRRForSFDD.hpp"
#include "EriDiagHRRForSGPD.hpp"
#include "EriDiagHRRForSGPF.hpp"
#include "EriDiagHRRForSGDD.hpp"
#include "EriDiagHRRForPDDD.hpp"
#include "EriDiagHRRForPFDD.hpp"
#include "EriDiagHRRForDDDD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostDDDD(      T*                                 intsBuffer,
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

    BufferHostXY<T> bvals(9, ncpairs);

    CBoysFunc<T, 8> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(9, ncpairs);

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
    BufferHostXY<T> pbufSSSP7(3, ncpairs);

0_0_0_z_7 0_0_0_y_7 0_0_0_x_7 
    BufferHostXY<T> pbufSSSD0(6, ncpairs);

0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 
    BufferHostXY<T> pbufSSSD1(6, ncpairs);

0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 
    BufferHostXY<T> pbufSSSD2(6, ncpairs);

0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 
    BufferHostXY<T> pbufSSSD3(6, ncpairs);

0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 
    BufferHostXY<T> pbufSSSD4(6, ncpairs);

0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 
    BufferHostXY<T> pbufSSSD5(4, ncpairs);

0_0_0_zz_5 0_0_0_yz_5 0_0_0_yy_5 0_0_0_xx_5 
    BufferHostXY<T> pbufSSSD6(3, ncpairs);

0_0_0_zz_6 0_0_0_yy_6 0_0_0_xx_6 
    BufferHostXY<T> pbufSSSF0(10, ncpairs);

0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 
    BufferHostXY<T> pbufSSSF1(10, ncpairs);

0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 
    BufferHostXY<T> pbufSSSF2(10, ncpairs);

0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 
    BufferHostXY<T> pbufSSSF3(10, ncpairs);

0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 
    BufferHostXY<T> pbufSSSF4(10, ncpairs);

0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 
    BufferHostXY<T> pbufSSSF5(8, ncpairs);

0_0_0_zzz_5 0_0_0_yzz_5 0_0_0_yyz_5 0_0_0_yyy_5 0_0_0_xzz_5 0_0_0_xyy_5 0_0_0_xxz_5 0_0_0_xxx_5 
    BufferHostXY<T> pbufSSSG0(15, ncpairs);

0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 
    BufferHostXY<T> pbufSSSG1(15, ncpairs);

0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 
    BufferHostXY<T> pbufSSSG2(15, ncpairs);

0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 
    BufferHostXY<T> pbufSSSG3(15, ncpairs);

0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 
    BufferHostXY<T> pbufSSSG4(15, ncpairs);

0_0_0_zzzz_4 0_0_0_yzzz_4 0_0_0_yyzz_4 0_0_0_yyyz_4 0_0_0_yyyy_4 0_0_0_xzzz_4 0_0_0_xyzz_4 0_0_0_xyyz_4 0_0_0_xyyy_4 0_0_0_xxzz_4 0_0_0_xxyz_4 0_0_0_xxyy_4 0_0_0_xxxz_4 0_0_0_xxxy_4 0_0_0_xxxx_4 
    BufferHostXY<T> pbufSPSS2(3, ncpairs);

0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 
    BufferHostXY<T> pbufSPSS3(3, ncpairs);

0_z_0_0_3 0_y_0_0_3 0_x_0_0_3 
    BufferHostXY<T> pbufSPSP1(9, ncpairs);

0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSPSP2(9, ncpairs);

0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 
    BufferHostXY<T> pbufSPSP3(9, ncpairs);

0_z_0_z_3 0_z_0_y_3 0_z_0_x_3 0_y_0_z_3 0_y_0_y_3 0_y_0_x_3 0_x_0_z_3 0_x_0_y_3 0_x_0_x_3 
    BufferHostXY<T> pbufSPSD0(12, ncpairs);

0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 
    BufferHostXY<T> pbufSPSD1(18, ncpairs);

0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 
    BufferHostXY<T> pbufSPSD2(18, ncpairs);

0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 
    BufferHostXY<T> pbufSPSD3(15, ncpairs);

0_z_0_zz_3 0_z_0_yz_3 0_z_0_yy_3 0_z_0_xz_3 0_z_0_xy_3 0_z_0_xx_3 0_y_0_yz_3 0_y_0_yy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xx_3 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xx_3 
    BufferHostXY<T> pbufSPSF0(18, ncpairs);

0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 
    BufferHostXY<T> pbufSPSF1(18, ncpairs);

0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 
    BufferHostXY<T> pbufSPSF2(18, ncpairs);

0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 
    BufferHostXY<T> pbufSPSF3(15, ncpairs);

0_z_0_zzz_3 0_z_0_yzz_3 0_z_0_yyz_3 0_z_0_xzz_3 0_z_0_xyz_3 0_z_0_xxz_3 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xxy_3 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxx_3 
    BufferHostXY<T> pbufSPSG0(25, ncpairs);

0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 
    BufferHostXY<T> pbufSPSG1(25, ncpairs);

0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 
    BufferHostXY<T> pbufSPSG2(20, ncpairs);

0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 
    BufferHostXY<T> pbufSPSG3(15, ncpairs);

0_z_0_zzzz_3 0_z_0_yzzz_3 0_z_0_yyzz_3 0_z_0_xzzz_3 0_z_0_xyzz_3 0_z_0_xxzz_3 0_y_0_yyyz_3 0_y_0_yyyy_3 0_y_0_xyyz_3 0_y_0_xyyy_3 0_y_0_xxyy_3 0_x_0_xxyz_3 0_x_0_xxxz_3 0_x_0_xxxy_3 0_x_0_xxxx_3 
    BufferHostXY<T> pbufSDSS2(3, ncpairs);

0_zz_0_0_2 0_yy_0_0_2 0_xx_0_0_2 
    BufferHostXY<T> pbufSDSP1(9, ncpairs);

0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 
    BufferHostXY<T> pbufSDSP2(9, ncpairs);

0_zz_0_z_2 0_zz_0_y_2 0_zz_0_x_2 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_x_2 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_x_2 
    BufferHostXY<T> pbufSDSD0(15, ncpairs);

0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 
    BufferHostXY<T> pbufSDSD1(15, ncpairs);

0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 
    BufferHostXY<T> pbufSDSD2(12, ncpairs);

0_zz_0_zz_2 0_zz_0_yz_2 0_zz_0_xz_2 0_zz_0_xy_2 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_xz_2 0_yy_0_xy_2 0_xx_0_yz_2 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xx_2 
    BufferHostXY<T> pbufSDSF0(27, ncpairs);

0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 
    BufferHostXY<T> pbufSDSF1(25, ncpairs);

0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 
    BufferHostXY<T> pbufSDSF2(15, ncpairs);

0_zz_0_zzz_2 0_zz_0_yzz_2 0_zz_0_yyz_2 0_zz_0_xzz_2 0_zz_0_xyz_2 0_zz_0_xxz_2 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_xyz_2 0_yy_0_xyy_2 0_yy_0_xxy_2 0_xx_0_xyz_2 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxx_2 
    BufferHostXY<T> pbufSDSG0(36, ncpairs);

0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 
    BufferHostXY<T> pbufSDSG1(21, ncpairs);

0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 
    BufferHostXY<T> pbufSDSG2(15, ncpairs);

0_zz_0_zzzz_2 0_zz_0_yzzz_2 0_zz_0_yyzz_2 0_zz_0_xzzz_2 0_zz_0_xyzz_2 0_zz_0_xxzz_2 0_yy_0_yyyz_2 0_yy_0_yyyy_2 0_yy_0_xyyz_2 0_yy_0_xyyy_2 0_yy_0_xxyy_2 0_xx_0_xxyz_2 0_xx_0_xxxz_2 0_xx_0_xxxy_2 0_xx_0_xxxx_2 
    BufferHostXY<T> pbufSFSP1(9, ncpairs);

0_zzz_0_z_1 0_yzz_0_z_1 0_yzz_0_y_1 0_yyz_0_z_1 0_yyy_0_y_1 0_xzz_0_z_1 0_xyy_0_y_1 0_xxz_0_z_1 0_xxx_0_x_1 
    BufferHostXY<T> pbufSFSD0(27, ncpairs);

0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 
    BufferHostXY<T> pbufSFSD1(24, ncpairs);

0_zzz_0_zz_1 0_zzz_0_yz_1 0_zzz_0_xz_1 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_xy_1 0_yyz_0_yz_1 0_yyz_0_yy_1 0_yyz_0_xz_1 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_xy_1 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xx_1 0_xxy_0_xy_1 0_xxy_0_xx_1 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 
    BufferHostXY<T> pbufSFSF0(46, ncpairs);

0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 
    BufferHostXY<T> pbufSFSF1(21, ncpairs);

0_zzz_0_zzz_1 0_zzz_0_yzz_1 0_zzz_0_xzz_1 0_yzz_0_yzz_1 0_yzz_0_yyz_1 0_yzz_0_xyz_1 0_yyz_0_yyz_1 0_yyz_0_xyz_1 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_xyy_1 0_xzz_0_xzz_1 0_xzz_0_xxz_1 0_xyy_0_xyy_1 0_xyy_0_xxy_1 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxy_0_xxy_1 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 
    BufferHostXY<T> pbufSFSG0(30, ncpairs);

0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxyy_0 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxxz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 
    BufferHostXY<T> pbufSFSG1(15, ncpairs);

0_zzz_0_zzzz_1 0_zzz_0_yzzz_1 0_zzz_0_xzzz_1 0_yzz_0_yyzz_1 0_yzz_0_xyzz_1 0_yyz_0_xyyz_1 0_yyy_0_yyyz_1 0_yyy_0_yyyy_1 0_yyy_0_xyyy_1 0_xzz_0_xxzz_1 0_xyy_0_xxyy_1 0_xxz_0_xxyz_1 0_xxx_0_xxxz_1 0_xxx_0_xxxy_1 0_xxx_0_xxxx_1 
    BufferHostXY<T> pbufSGSD0(36, ncpairs);

0_zzzz_0_zz_0 0_yzzz_0_zz_0 0_yzzz_0_yz_0 0_yyzz_0_zz_0 0_yyzz_0_yz_0 0_yyzz_0_yy_0 0_yyyz_0_yz_0 0_yyyz_0_yy_0 0_yyyy_0_yy_0 0_xzzz_0_zz_0 0_xzzz_0_xz_0 0_xyzz_0_zz_0 0_xyzz_0_yz_0 0_xyzz_0_xz_0 0_xyzz_0_xy_0 0_xyyz_0_yz_0 0_xyyz_0_yy_0 0_xyyz_0_xz_0 0_xyyz_0_xy_0 0_xyyy_0_yy_0 0_xyyy_0_xy_0 0_xxzz_0_zz_0 0_xxzz_0_xz_0 0_xxzz_0_xx_0 0_xxyz_0_yz_0 0_xxyz_0_xz_0 0_xxyz_0_xy_0 0_xxyz_0_xx_0 0_xxyy_0_yy_0 0_xxyy_0_xy_0 0_xxyy_0_xx_0 0_xxxz_0_xz_0 0_xxxz_0_xx_0 0_xxxy_0_xy_0 0_xxxy_0_xx_0 0_xxxx_0_xx_0 
    BufferHostXY<T> pbufSGSF0(30, ncpairs);

0_zzzz_0_zzz_0 0_yzzz_0_zzz_0 0_yzzz_0_yzz_0 0_yyzz_0_yzz_0 0_yyzz_0_yyz_0 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyy_0_yyy_0 0_xzzz_0_zzz_0 0_xzzz_0_xzz_0 0_xyzz_0_yzz_0 0_xyzz_0_xzz_0 0_xyzz_0_xyz_0 0_xyyz_0_yyz_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyy_0_yyy_0 0_xyyy_0_xyy_0 0_xxzz_0_xzz_0 0_xxzz_0_xxz_0 0_xxyz_0_xyz_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyy_0_xyy_0 0_xxyy_0_xxy_0 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxx_0_xxx_0 
    BufferHostXY<T> pbufSGSG0(15, ncpairs);

0_zzzz_0_zzzz_0 0_yzzz_0_yzzz_0 0_yyzz_0_yyzz_0 0_yyyz_0_yyyz_0 0_yyyy_0_yyyy_0 0_xzzz_0_xzzz_0 0_xyzz_0_xyzz_0 0_xyyz_0_xyyz_0 0_xyyy_0_xyyy_0 0_xxzz_0_xxzz_0 0_xxyz_0_xxyz_0 0_xxyy_0_xxyy_0 0_xxxz_0_xxxz_0 0_xxxy_0_xxxy_0 0_xxxx_0_xxxx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSD(15, ncpairs);

0_zz_0_zz 0_zz_0_yz 0_zz_0_xz 0_zz_0_xy 0_yz_0_yz 0_yy_0_yz 0_yy_0_yy 0_yy_0_xz 0_yy_0_xy 0_xz_0_xz 0_xy_0_xy 0_xx_0_yz 0_xx_0_xz 0_xx_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSDSF(27, ncpairs);

0_zz_0_zzz 0_zz_0_yzz 0_zz_0_yyz 0_zz_0_xzz 0_zz_0_xyz 0_zz_0_xxz 0_yz_0_yzz 0_yz_0_yyz 0_yz_0_xyz 0_yy_0_yzz 0_yy_0_yyz 0_yy_0_yyy 0_yy_0_xyz 0_yy_0_xyy 0_yy_0_xxy 0_xz_0_xzz 0_xz_0_xyz 0_xz_0_xxz 0_xy_0_xyz 0_xy_0_xyy 0_xy_0_xxy 0_xx_0_xzz 0_xx_0_xyz 0_xx_0_xyy 0_xx_0_xxz 0_xx_0_xxy 0_xx_0_xxx 
    BufferHostXY<T> cbufSDSG(36, ncpairs);

0_zz_0_zzzz 0_zz_0_yzzz 0_zz_0_yyzz 0_zz_0_xzzz 0_zz_0_xyzz 0_zz_0_xxzz 0_yz_0_yzzz 0_yz_0_yyzz 0_yz_0_yyyz 0_yz_0_xyzz 0_yz_0_xyyz 0_yz_0_xxyz 0_yy_0_yyzz 0_yy_0_yyyz 0_yy_0_yyyy 0_yy_0_xyyz 0_yy_0_xyyy 0_yy_0_xxyy 0_xz_0_xzzz 0_xz_0_xyzz 0_xz_0_xyyz 0_xz_0_xxzz 0_xz_0_xxyz 0_xz_0_xxxz 0_xy_0_xyzz 0_xy_0_xyyz 0_xy_0_xyyy 0_xy_0_xxyz 0_xy_0_xxyy 0_xy_0_xxxy 0_xx_0_xxzz 0_xx_0_xxyz 0_xx_0_xxyy 0_xx_0_xxxz 0_xx_0_xxxy 0_xx_0_xxxx 
    BufferHostXY<T> cbufSFSD(27, ncpairs);

0_zzz_0_zz 0_zzz_0_yz 0_zzz_0_xz 0_yzz_0_zz 0_yzz_0_yz 0_yzz_0_xy 0_yyz_0_yz 0_yyz_0_yy 0_yyz_0_xz 0_yyy_0_yz 0_yyy_0_yy 0_yyy_0_xy 0_xzz_0_zz 0_xzz_0_xz 0_xyz_0_yz 0_xyz_0_xz 0_xyz_0_xy 0_xyy_0_yy 0_xyy_0_xy 0_xxz_0_yz 0_xxz_0_xz 0_xxz_0_xx 0_xxy_0_xy 0_xxy_0_xx 0_xxx_0_xz 0_xxx_0_xy 0_xxx_0_xx 
    BufferHostXY<T> cbufSFSF(46, ncpairs);

0_zzz_0_zzz 0_zzz_0_yzz 0_zzz_0_xzz 0_yzz_0_zzz 0_yzz_0_yzz 0_yzz_0_yyz 0_yzz_0_xzz 0_yzz_0_xyz 0_yyz_0_yzz 0_yyz_0_yyz 0_yyz_0_yyy 0_yyz_0_xyz 0_yyz_0_xyy 0_yyy_0_yyz 0_yyy_0_yyy 0_yyy_0_xyy 0_xzz_0_zzz 0_xzz_0_yzz 0_xzz_0_xzz 0_xzz_0_xyz 0_xzz_0_xxz 0_xyz_0_yzz 0_xyz_0_yyz 0_xyz_0_xzz 0_xyz_0_xyz 0_xyz_0_xyy 0_xyz_0_xxz 0_xyz_0_xxy 0_xyy_0_yyz 0_xyy_0_yyy 0_xyy_0_xyz 0_xyy_0_xyy 0_xyy_0_xxy 0_xxz_0_xzz 0_xxz_0_xyz 0_xxz_0_xxz 0_xxz_0_xxy 0_xxz_0_xxx 0_xxy_0_xyz 0_xxy_0_xyy 0_xxy_0_xxz 0_xxy_0_xxy 0_xxy_0_xxx 0_xxx_0_xxz 0_xxx_0_xxy 0_xxx_0_xxx 
    BufferHostXY<T> cbufSFSG(30, ncpairs);

0_zzz_0_zzzz 0_zzz_0_yzzz 0_zzz_0_xzzz 0_yzz_0_yzzz 0_yzz_0_yyzz 0_yzz_0_xyzz 0_yyz_0_yyzz 0_yyz_0_yyyz 0_yyz_0_xyyz 0_yyy_0_yyyz 0_yyy_0_yyyy 0_yyy_0_xyyy 0_xzz_0_xzzz 0_xzz_0_xyzz 0_xzz_0_xxzz 0_xyz_0_xyzz 0_xyz_0_xyyz 0_xyz_0_xxyz 0_xyy_0_xyyz 0_xyy_0_xyyy 0_xyy_0_xxyy 0_xxz_0_xxzz 0_xxz_0_xxyz 0_xxz_0_xxxz 0_xxy_0_xxyz 0_xxy_0_xxyy 0_xxy_0_xxxy 0_xxx_0_xxxz 0_xxx_0_xxxy 0_xxx_0_xxxx 
    BufferHostXY<T> cbufSGSD(36, ncpairs);

0_zzzz_0_zz 0_yzzz_0_zz 0_yzzz_0_yz 0_yyzz_0_zz 0_yyzz_0_yz 0_yyzz_0_yy 0_yyyz_0_yz 0_yyyz_0_yy 0_yyyy_0_yy 0_xzzz_0_zz 0_xzzz_0_xz 0_xyzz_0_zz 0_xyzz_0_yz 0_xyzz_0_xz 0_xyzz_0_xy 0_xyyz_0_yz 0_xyyz_0_yy 0_xyyz_0_xz 0_xyyz_0_xy 0_xyyy_0_yy 0_xyyy_0_xy 0_xxzz_0_zz 0_xxzz_0_xz 0_xxzz_0_xx 0_xxyz_0_yz 0_xxyz_0_xz 0_xxyz_0_xy 0_xxyz_0_xx 0_xxyy_0_yy 0_xxyy_0_xy 0_xxyy_0_xx 0_xxxz_0_xz 0_xxxz_0_xx 0_xxxy_0_xy 0_xxxy_0_xx 0_xxxx_0_xx 
    BufferHostXY<T> cbufSGSF(30, ncpairs);

0_zzzz_0_zzz 0_yzzz_0_zzz 0_yzzz_0_yzz 0_yyzz_0_yzz 0_yyzz_0_yyz 0_yyyz_0_yyz 0_yyyz_0_yyy 0_yyyy_0_yyy 0_xzzz_0_zzz 0_xzzz_0_xzz 0_xyzz_0_yzz 0_xyzz_0_xzz 0_xyzz_0_xyz 0_xyyz_0_yyz 0_xyyz_0_xyz 0_xyyz_0_xyy 0_xyyy_0_yyy 0_xyyy_0_xyy 0_xxzz_0_xzz 0_xxzz_0_xxz 0_xxyz_0_xyz 0_xxyz_0_xxz 0_xxyz_0_xxy 0_xxyy_0_xyy 0_xxyy_0_xxy 0_xxxz_0_xxz 0_xxxz_0_xxx 0_xxxy_0_xxy 0_xxxy_0_xxx 0_xxxx_0_xxx 
    BufferHostXY<T> cbufSGSG(15, ncpairs);

0_zzzz_0_zzzz 0_yzzz_0_yzzz 0_yyzz_0_yyzz 0_yyyz_0_yyyz 0_yyyy_0_yyyy 0_xzzz_0_xzzz 0_xyzz_0_xyzz 0_xyyz_0_xyyz 0_xyyy_0_xyyy 0_xxzz_0_xxzz 0_xxyz_0_xxyz 0_xxyy_0_xxyy 0_xxxz_0_xxxz 0_xxxy_0_xxxy 0_xxxx_0_xxxx 
    BufferHostXY<T> cbufSDPD(18, ncpairs);

0_zz_z_zz 0_zz_y_zz 0_zz_x_zz 0_yz_z_yz 0_yz_y_yz 0_yz_x_yz 0_yy_z_yy 0_yy_y_yy 0_yy_x_yy 0_xz_z_xz 0_xz_y_xz 0_xz_x_xz 0_xy_z_xy 0_xy_y_xy 0_xy_x_xy 0_xx_z_xx 0_xx_y_xx 0_xx_x_xx 
    BufferHostXY<T> cbufSDPF(36, ncpairs);

0_zz_z_zzz 0_zz_y_zzz 0_zz_y_yzz 0_zz_x_zzz 0_zz_x_yzz 0_zz_x_xzz 0_yz_z_yzz 0_yz_y_yzz 0_yz_y_yyz 0_yz_x_yzz 0_yz_x_yyz 0_yz_x_xyz 0_yy_z_yyz 0_yy_y_yyz 0_yy_y_yyy 0_yy_x_yyz 0_yy_x_yyy 0_yy_x_xyy 0_xz_z_xzz 0_xz_y_xzz 0_xz_y_xyz 0_xz_x_xzz 0_xz_x_xyz 0_xz_x_xxz 0_xy_z_xyz 0_xy_y_xyz 0_xy_y_xyy 0_xy_x_xyz 0_xy_x_xyy 0_xy_x_xxy 0_xx_z_xxz 0_xx_y_xxz 0_xx_y_xxy 0_xx_x_xxz 0_xx_x_xxy 0_xx_x_xxx 
    BufferHostXY<T> cbufSDDD(36, ncpairs);

0_zz_zz_zz 0_zz_yz_zz 0_zz_yy_zz 0_zz_xz_zz 0_zz_xy_zz 0_zz_xx_zz 0_yz_zz_yz 0_yz_yz_yz 0_yz_yy_yz 0_yz_xz_yz 0_yz_xy_yz 0_yz_xx_yz 0_yy_zz_yy 0_yy_yz_yy 0_yy_yy_yy 0_yy_xz_yy 0_yy_xy_yy 0_yy_xx_yy 0_xz_zz_xz 0_xz_yz_xz 0_xz_yy_xz 0_xz_xz_xz 0_xz_xy_xz 0_xz_xx_xz 0_xy_zz_xy 0_xy_yz_xy 0_xy_yy_xy 0_xy_xz_xy 0_xy_xy_xy 0_xy_xx_xy 0_xx_zz_xx 0_xx_yz_xx 0_xx_yy_xx 0_xx_xz_xx 0_xx_xy_xx 0_xx_xx_xx 
    BufferHostXY<T> cbufSFPD(36, ncpairs);

0_zzz_z_zz 0_zzz_y_zz 0_zzz_x_zz 0_yzz_z_yz 0_yzz_y_zz 0_yzz_y_yz 0_yzz_x_zz 0_yzz_x_yz 0_yyz_z_yy 0_yyz_y_yz 0_yyz_y_yy 0_yyz_x_yz 0_yyz_x_yy 0_yyy_y_yy 0_yyy_x_yy 0_xzz_z_xz 0_xzz_y_xz 0_xzz_x_zz 0_xzz_x_xz 0_xyz_z_xy 0_xyz_y_xz 0_xyz_y_xy 0_xyz_x_yz 0_xyz_x_xz 0_xyz_x_xy 0_xyy_y_xy 0_xyy_x_yy 0_xyy_x_xy 0_xxz_z_xx 0_xxz_y_xx 0_xxz_x_xz 0_xxz_x_xx 0_xxy_y_xx 0_xxy_x_xy 0_xxy_x_xx 0_xxx_x_xx 
    BufferHostXY<T> cbufSFPF(43, ncpairs);

0_zzz_z_zzz 0_zzz_y_zzz 0_zzz_x_zzz 0_yzz_z_yzz 0_yzz_y_zzz 0_yzz_y_yzz 0_yzz_x_yzz 0_yyz_z_yyz 0_yyz_y_yzz 0_yyz_y_yyz 0_yyz_x_yyz 0_yyy_y_yyz 0_yyy_y_yyy 0_yyy_x_yyy 0_xzz_z_xzz 0_xzz_y_xzz 0_xzz_x_zzz 0_xzz_x_yzz 0_xzz_x_xzz 0_xyz_z_xyz 0_xyz_y_xzz 0_xyz_y_xyz 0_xyz_x_yzz 0_xyz_x_yyz 0_xyz_x_xyz 0_xyy_y_xyz 0_xyy_y_xyy 0_xyy_x_yyz 0_xyy_x_yyy 0_xyy_x_xyy 0_xxz_z_xxz 0_xxz_y_xxz 0_xxz_x_xzz 0_xxz_x_xyz 0_xxz_x_xxz 0_xxy_y_xxz 0_xxy_y_xxy 0_xxy_x_xyz 0_xxy_x_xyy 0_xxy_x_xxy 0_xxx_x_xxz 0_xxx_x_xxy 0_xxx_x_xxx 
    BufferHostXY<T> cbufSFDD(54, ncpairs);

0_zzz_zz_zz 0_zzz_yz_zz 0_zzz_xz_zz 0_yzz_zz_yz 0_yzz_yz_zz 0_yzz_yz_yz 0_yzz_yy_zz 0_yzz_xz_yz 0_yzz_xy_zz 0_yyz_zz_yy 0_yyz_yz_yz 0_yyz_yz_yy 0_yyz_yy_yz 0_yyz_xz_yy 0_yyz_xy_yz 0_yyy_yz_yy 0_yyy_yy_yy 0_yyy_xy_yy 0_xzz_zz_xz 0_xzz_yz_xz 0_xzz_xz_zz 0_xzz_xz_xz 0_xzz_xy_zz 0_xzz_xx_zz 0_xyz_zz_xy 0_xyz_yz_xz 0_xyz_yz_xy 0_xyz_yy_xz 0_xyz_xz_yz 0_xyz_xz_xy 0_xyz_xy_yz 0_xyz_xy_xz 0_xyz_xx_yz 0_xyy_yz_xy 0_xyy_yy_xy 0_xyy_xz_yy 0_xyy_xy_yy 0_xyy_xy_xy 0_xyy_xx_yy 0_xxz_zz_xx 0_xxz_yz_xx 0_xxz_xz_xz 0_xxz_xz_xx 0_xxz_xy_xz 0_xxz_xx_xz 0_xxy_yz_xx 0_xxy_yy_xx 0_xxy_xz_xy 0_xxy_xy_xy 0_xxy_xy_xx 0_xxy_xx_xy 0_xxx_xz_xx 0_xxx_xy_xx 0_xxx_xx_xx 
    BufferHostXY<T> cbufSGPD(36, ncpairs);

0_zzzz_z_zz 0_yzzz_z_yz 0_yzzz_y_zz 0_yyzz_z_yy 0_yyzz_y_zz 0_yyzz_y_yz 0_yyyz_y_yz 0_yyyz_y_yy 0_yyyy_y_yy 0_xzzz_z_xz 0_xzzz_x_zz 0_xyzz_z_xy 0_xyzz_y_xz 0_xyzz_x_zz 0_xyzz_x_yz 0_xyyz_y_xz 0_xyyz_y_xy 0_xyyz_x_yz 0_xyyz_x_yy 0_xyyy_y_xy 0_xyyy_x_yy 0_xxzz_z_xx 0_xxzz_x_zz 0_xxzz_x_xz 0_xxyz_y_xx 0_xxyz_x_yz 0_xxyz_x_xz 0_xxyz_x_xy 0_xxyy_y_xx 0_xxyy_x_yy 0_xxyy_x_xy 0_xxxz_x_xz 0_xxxz_x_xx 0_xxxy_x_xy 0_xxxy_x_xx 0_xxxx_x_xx 
    BufferHostXY<T> cbufSGPF(25, ncpairs);

0_zzzz_z_zzz 0_yzzz_z_yzz 0_yzzz_y_zzz 0_yyzz_z_yyz 0_yyzz_y_yzz 0_yyyz_y_yyz 0_yyyy_y_yyy 0_xzzz_z_xzz 0_xzzz_x_zzz 0_xyzz_z_xyz 0_xyzz_y_xzz 0_xyzz_x_yzz 0_xyyz_y_xyz 0_xyyz_x_yyz 0_xyyy_y_xyy 0_xyyy_x_yyy 0_xxzz_z_xxz 0_xxzz_x_xzz 0_xxyz_y_xxz 0_xxyz_x_xyz 0_xxyy_y_xxy 0_xxyy_x_xyy 0_xxxz_x_xxz 0_xxxy_x_xxy 0_xxxx_x_xxx 
    BufferHostXY<T> cbufSGDD(36, ncpairs);

0_zzzz_zz_zz 0_yzzz_zz_yz 0_yzzz_yz_zz 0_yyzz_zz_yy 0_yyzz_yz_yz 0_yyzz_yy_zz 0_yyyz_yz_yy 0_yyyz_yy_yz 0_yyyy_yy_yy 0_xzzz_zz_xz 0_xzzz_xz_zz 0_xyzz_zz_xy 0_xyzz_yz_xz 0_xyzz_xz_yz 0_xyzz_xy_zz 0_xyyz_yz_xy 0_xyyz_yy_xz 0_xyyz_xz_yy 0_xyyz_xy_yz 0_xyyy_yy_xy 0_xyyy_xy_yy 0_xxzz_zz_xx 0_xxzz_xz_xz 0_xxzz_xx_zz 0_xxyz_yz_xx 0_xxyz_xz_xy 0_xxyz_xy_xz 0_xxyz_xx_yz 0_xxyy_yy_xx 0_xxyy_xy_xy 0_xxyy_xx_yy 0_xxxz_xz_xx 0_xxxz_xx_xz 0_xxxy_xy_xx 0_xxxy_xx_xy 0_xxxx_xx_xx 
    BufferHostXY<T> cbufPDDD(36, ncpairs);

z_zz_zz_zz z_yz_zz_yz z_yy_zz_yy z_xz_zz_xz z_xy_zz_xy z_xx_zz_xx y_zz_yz_zz y_zz_yy_zz y_yz_yz_yz y_yz_yy_yz y_yy_yz_yy y_yy_yy_yy y_xz_yz_xz y_xz_yy_xz y_xy_yz_xy y_xy_yy_xy y_xx_yz_xx y_xx_yy_xx x_zz_xz_zz x_zz_xy_zz x_zz_xx_zz x_yz_xz_yz x_yz_xy_yz x_yz_xx_yz x_yy_xz_yy x_yy_xy_yy x_yy_xx_yy x_xz_xz_xz x_xz_xy_xz x_xz_xx_xz x_xy_xz_xy x_xy_xy_xy x_xy_xx_xy x_xx_xz_xx x_xx_xy_xx x_xx_xx_xx 
    BufferHostXY<T> cbufPFDD(36, ncpairs);

z_zzz_zz_zz z_yzz_zz_yz z_yyz_zz_yy z_xzz_zz_xz z_xyz_zz_xy z_xxz_zz_xx y_zzz_yz_zz y_yzz_yz_yz y_yzz_yy_zz y_yyz_yz_yy y_yyz_yy_yz y_yyy_yy_yy y_xzz_yz_xz y_xyz_yz_xy y_xyz_yy_xz y_xyy_yy_xy y_xxz_yz_xx y_xxy_yy_xx x_zzz_xz_zz x_yzz_xz_yz x_yzz_xy_zz x_yyz_xz_yy x_yyz_xy_yz x_yyy_xy_yy x_xzz_xz_xz x_xzz_xx_zz x_xyz_xz_xy x_xyz_xy_xz x_xyz_xx_yz x_xyy_xy_xy x_xyy_xx_yy x_xxz_xz_xx x_xxz_xx_xz x_xxy_xy_xx x_xxy_xx_xy x_xxx_xx_xx 
    BufferHostXY<T> cbufDDDD(36, ncpairs);

zz_zz_zz_zz zz_yz_zz_yz zz_yy_zz_yy zz_xz_zz_xz zz_xy_zz_xy zz_xx_zz_xx yz_zz_yz_zz yz_yz_yz_yz yz_yy_yz_yy yz_xz_yz_xz yz_xy_yz_xy yz_xx_yz_xx yy_zz_yy_zz yy_yz_yy_yz yy_yy_yy_yy yy_xz_yy_xz yy_xy_yy_xy yy_xx_yy_xx xz_zz_xz_zz xz_yz_xz_yz xz_yy_xz_yy xz_xz_xz_xz xz_xy_xz_xy xz_xx_xz_xx xy_zz_xy_zz xy_yz_xy_yz xy_yy_xy_yy xy_xz_xy_xz xy_xy_xy_xy xy_xx_xy_xx xx_zz_xx_zz xx_yz_xx_yz xx_yy_xx_yy xx_xz_xx_xz xx_xy_xx_xy xx_xx_xx_xx 
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
z_zz_zz_zz_0 z_zzz_zz_zz_0 z_yz_zz_yz_0 z_yzz_zz_yz_0 z_yy_zz_yy_0 z_yyz_zz_yy_0 z_xz_zz_xz_0 z_xzz_zz_xz_0 z_xy_zz_xy_0 z_xyz_zz_xy_0 z_xx_zz_xx_0 z_xxz_zz_xx_0 zz_zz_zz_zz_0 zz_yz_zz_yz_0 zz_yy_zz_yy_0 zz_xz_zz_xz_0 zz_xy_zz_xy_0 zz_xx_zz_xx_0 y_zz_yz_zz_0 y_zz_yy_zz_0 y_zzz_yz_zz_0 y_yz_yz_yz_0 y_yz_yy_yz_0 y_yzz_yz_yz_0 y_yzz_yy_zz_0 y_yy_yz_yy_0 y_yy_yy_yy_0 y_yyz_yz_yy_0 y_yyz_yy_yz_0 y_yyy_yy_yy_0 y_xz_yz_xz_0 y_xz_yy_xz_0 y_xzz_yz_xz_0 y_xy_yz_xy_0 y_xy_yy_xy_0 y_xyz_yz_xy_0 y_xyz_yy_xz_0 y_xyy_yy_xy_0 y_xx_yz_xx_0 y_xx_yy_xx_0 y_xxz_yz_xx_0 y_xxy_yy_xx_0 yz_zz_yz_zz_0 yz_yz_yz_yz_0 yz_yy_yz_yy_0 yz_xz_yz_xz_0 yz_xy_yz_xy_0 yz_xx_yz_xx_0 yy_zz_yy_zz_0 yy_yz_yy_yz_0 yy_yy_yy_yy_0 yy_xz_yy_xz_0 yy_xy_yy_xy_0 yy_xx_yy_xx_0 x_zz_xz_zz_0 x_zz_xy_zz_0 x_zz_xx_zz_0 x_zzz_xz_zz_0 x_yz_xz_yz_0 x_yz_xy_yz_0 x_yz_xx_yz_0 x_yzz_xz_yz_0 x_yzz_xy_zz_0 x_yy_xz_yy_0 x_yy_xy_yy_0 x_yy_xx_yy_0 x_yyz_xz_yy_0 x_yyz_xy_yz_0 x_yyy_xy_yy_0 x_xz_xz_xz_0 x_xz_xy_xz_0 x_xz_xx_xz_0 x_xzz_xz_xz_0 x_xzz_xx_zz_0 x_xy_xz_xy_0 x_xy_xy_xy_0 x_xy_xx_xy_0 x_xyz_xz_xy_0 x_xyz_xy_xz_0 x_xyz_xx_yz_0 x_xyy_xy_xy_0 x_xyy_xx_yy_0 x_xx_xz_xx_0 x_xx_xy_xx_0 x_xx_xx_xx_0 x_xxz_xz_xx_0 x_xxz_xx_xz_0 x_xxy_xy_xx_0 x_xxy_xx_xy_0 x_xxx_xx_xx_0 xz_zz_xz_zz_0 xz_yz_xz_yz_0 xz_yy_xz_yy_0 xz_xz_xz_xz_0 xz_xy_xz_xy_0 xz_xx_xz_xx_0 xy_zz_xy_zz_0 xy_yz_xy_yz_0 xy_yy_xy_yy_0 xy_xz_xy_xz_0 xy_xy_xy_xy_0 xy_xx_xy_xx_0 xx_zz_xx_zz_0 xx_yz_xx_yz_0 xx_yy_xx_yy_0 xx_xz_xx_xz_0 xx_xy_xx_xy_0 xx_xx_xx_xx_0 

signature:
0_zzz_zz_zz_0 0_zzz_yz_zz_0 0_zzz_xz_zz_0 0_zzzz_zz_zz_0 0_yzz_zz_yz_0 0_yzz_yz_zz_0 0_yzz_yz_yz_0 0_yzz_yy_zz_0 0_yzz_xz_yz_0 0_yzz_xy_zz_0 0_yzzz_zz_yz_0 0_yzzz_yz_zz_0 0_yyz_zz_yy_0 0_yyz_yz_yz_0 0_yyz_yz_yy_0 0_yyz_yy_yz_0 0_yyz_xz_yy_0 0_yyz_xy_yz_0 0_yyzz_zz_yy_0 0_yyzz_yz_yz_0 0_yyzz_yy_zz_0 0_yyy_yz_yy_0 0_yyy_yy_yy_0 0_yyy_xy_yy_0 0_yyyz_yz_yy_0 0_yyyz_yy_yz_0 0_yyyy_yy_yy_0 0_xzz_zz_xz_0 0_xzz_yz_xz_0 0_xzz_xz_zz_0 0_xzz_xz_xz_0 0_xzz_xy_zz_0 0_xzz_xx_zz_0 0_xzzz_zz_xz_0 0_xzzz_xz_zz_0 0_xyz_zz_xy_0 0_xyz_yz_xz_0 0_xyz_yz_xy_0 0_xyz_yy_xz_0 0_xyz_xz_yz_0 0_xyz_xz_xy_0 0_xyz_xy_yz_0 0_xyz_xy_xz_0 0_xyz_xx_yz_0 0_xyzz_zz_xy_0 0_xyzz_yz_xz_0 0_xyzz_xz_yz_0 0_xyzz_xy_zz_0 0_xyy_yz_xy_0 0_xyy_yy_xy_0 0_xyy_xz_yy_0 0_xyy_xy_yy_0 0_xyy_xy_xy_0 0_xyy_xx_yy_0 0_xyyz_yz_xy_0 0_xyyz_yy_xz_0 0_xyyz_xz_yy_0 0_xyyz_xy_yz_0 0_xyyy_yy_xy_0 0_xyyy_xy_yy_0 0_xxz_zz_xx_0 0_xxz_yz_xx_0 0_xxz_xz_xz_0 0_xxz_xz_xx_0 0_xxz_xy_xz_0 0_xxz_xx_xz_0 0_xxzz_zz_xx_0 0_xxzz_xz_xz_0 0_xxzz_xx_zz_0 0_xxy_yz_xx_0 0_xxy_yy_xx_0 0_xxy_xz_xy_0 0_xxy_xy_xy_0 0_xxy_xy_xx_0 0_xxy_xx_xy_0 0_xxyz_yz_xx_0 0_xxyz_xz_xy_0 0_xxyz_xy_xz_0 0_xxyz_xx_yz_0 0_xxyy_yy_xx_0 0_xxyy_xy_xy_0 0_xxyy_xx_yy_0 0_xxx_xz_xx_0 0_xxx_xy_xx_0 0_xxx_xx_xx_0 0_xxxz_xz_xx_0 0_xxxz_xx_xz_0 0_xxxy_xy_xx_0 0_xxxy_xx_xy_0 0_xxxx_xx_xx_0 z_zzz_zz_zz_0 z_yzz_zz_yz_0 z_yyz_zz_yy_0 z_xzz_zz_xz_0 z_xyz_zz_xy_0 z_xxz_zz_xx_0 y_zzz_yz_zz_0 y_yzz_yz_yz_0 y_yzz_yy_zz_0 y_yyz_yz_yy_0 y_yyz_yy_yz_0 y_yyy_yy_yy_0 y_xzz_yz_xz_0 y_xyz_yz_xy_0 y_xyz_yy_xz_0 y_xyy_yy_xy_0 y_xxz_yz_xx_0 y_xxy_yy_xx_0 x_zzz_xz_zz_0 x_yzz_xz_yz_0 x_yzz_xy_zz_0 x_yyz_xz_yy_0 x_yyz_xy_yz_0 x_yyy_xy_yy_0 x_xzz_xz_xz_0 x_xzz_xx_zz_0 x_xyz_xz_xy_0 x_xyz_xy_xz_0 x_xyz_xx_yz_0 x_xyy_xy_xy_0 x_xyy_xx_yy_0 x_xxz_xz_xx_0 x_xxz_xx_xz_0 x_xxy_xy_xx_0 x_xxy_xx_xy_0 x_xxx_xx_xx_0 

signature:
0_zz_zz_zz_0 0_zz_yz_zz_0 0_zz_yy_zz_0 0_zz_xz_zz_0 0_zz_xy_zz_0 0_zz_xx_zz_0 0_zzz_zz_zz_0 0_zzz_yz_zz_0 0_zzz_xz_zz_0 0_yz_zz_yz_0 0_yz_yz_yz_0 0_yz_yy_yz_0 0_yz_xz_yz_0 0_yz_xy_yz_0 0_yz_xx_yz_0 0_yzz_zz_yz_0 0_yzz_yz_zz_0 0_yzz_yz_yz_0 0_yzz_yy_zz_0 0_yzz_xz_yz_0 0_yzz_xy_zz_0 0_yy_zz_yy_0 0_yy_yz_yy_0 0_yy_yy_yy_0 0_yy_xz_yy_0 0_yy_xy_yy_0 0_yy_xx_yy_0 0_yyz_zz_yy_0 0_yyz_yz_yz_0 0_yyz_yz_yy_0 0_yyz_yy_yz_0 0_yyz_xz_yy_0 0_yyz_xy_yz_0 0_yyy_yz_yy_0 0_yyy_yy_yy_0 0_yyy_xy_yy_0 0_xz_zz_xz_0 0_xz_yz_xz_0 0_xz_yy_xz_0 0_xz_xz_xz_0 0_xz_xy_xz_0 0_xz_xx_xz_0 0_xzz_zz_xz_0 0_xzz_yz_xz_0 0_xzz_xz_zz_0 0_xzz_xz_xz_0 0_xzz_xy_zz_0 0_xzz_xx_zz_0 0_xy_zz_xy_0 0_xy_yz_xy_0 0_xy_yy_xy_0 0_xy_xz_xy_0 0_xy_xy_xy_0 0_xy_xx_xy_0 0_xyz_zz_xy_0 0_xyz_yz_xz_0 0_xyz_yz_xy_0 0_xyz_yy_xz_0 0_xyz_xz_yz_0 0_xyz_xz_xy_0 0_xyz_xy_yz_0 0_xyz_xy_xz_0 0_xyz_xx_yz_0 0_xyy_yz_xy_0 0_xyy_yy_xy_0 0_xyy_xz_yy_0 0_xyy_xy_yy_0 0_xyy_xy_xy_0 0_xyy_xx_yy_0 0_xx_zz_xx_0 0_xx_yz_xx_0 0_xx_yy_xx_0 0_xx_xz_xx_0 0_xx_xy_xx_0 0_xx_xx_xx_0 0_xxz_zz_xx_0 0_xxz_yz_xx_0 0_xxz_xz_xz_0 0_xxz_xz_xx_0 0_xxz_xy_xz_0 0_xxz_xx_xz_0 0_xxy_yz_xx_0 0_xxy_yy_xx_0 0_xxy_xz_xy_0 0_xxy_xy_xy_0 0_xxy_xy_xx_0 0_xxy_xx_xy_0 0_xxx_xz_xx_0 0_xxx_xy_xx_0 0_xxx_xx_xx_0 z_zz_zz_zz_0 z_yz_zz_yz_0 z_yy_zz_yy_0 z_xz_zz_xz_0 z_xy_zz_xy_0 z_xx_zz_xx_0 y_zz_yz_zz_0 y_zz_yy_zz_0 y_yz_yz_yz_0 y_yz_yy_yz_0 y_yy_yz_yy_0 y_yy_yy_yy_0 y_xz_yz_xz_0 y_xz_yy_xz_0 y_xy_yz_xy_0 y_xy_yy_xy_0 y_xx_yz_xx_0 y_xx_yy_xx_0 x_zz_xz_zz_0 x_zz_xy_zz_0 x_zz_xx_zz_0 x_yz_xz_yz_0 x_yz_xy_yz_0 x_yz_xx_yz_0 x_yy_xz_yy_0 x_yy_xy_yy_0 x_yy_xx_yy_0 x_xz_xz_xz_0 x_xz_xy_xz_0 x_xz_xx_xz_0 x_xy_xz_xy_0 x_xy_xy_xy_0 x_xy_xx_xy_0 x_xx_xz_xx_0 x_xx_xy_xx_0 x_xx_xx_xx_0 

signature:
0_zzzz_z_zz_0 0_zzzz_z_zzz_0 0_zzzz_zz_zz_0 0_yzzz_z_yz_0 0_yzzz_z_yzz_0 0_yzzz_zz_yz_0 0_yzzz_y_zz_0 0_yzzz_y_zzz_0 0_yzzz_yz_zz_0 0_yyzz_z_yy_0 0_yyzz_z_yyz_0 0_yyzz_zz_yy_0 0_yyzz_y_zz_0 0_yyzz_y_yz_0 0_yyzz_y_yzz_0 0_yyzz_yz_yz_0 0_yyzz_yy_zz_0 0_yyyz_y_yz_0 0_yyyz_y_yy_0 0_yyyz_y_yyz_0 0_yyyz_yz_yy_0 0_yyyz_yy_yz_0 0_yyyy_y_yy_0 0_yyyy_y_yyy_0 0_yyyy_yy_yy_0 0_xzzz_z_xz_0 0_xzzz_z_xzz_0 0_xzzz_zz_xz_0 0_xzzz_x_zz_0 0_xzzz_x_zzz_0 0_xzzz_xz_zz_0 0_xyzz_z_xy_0 0_xyzz_z_xyz_0 0_xyzz_zz_xy_0 0_xyzz_y_xz_0 0_xyzz_y_xzz_0 0_xyzz_yz_xz_0 0_xyzz_x_zz_0 0_xyzz_x_yz_0 0_xyzz_x_yzz_0 0_xyzz_xz_yz_0 0_xyzz_xy_zz_0 0_xyyz_y_xz_0 0_xyyz_y_xy_0 0_xyyz_y_xyz_0 0_xyyz_yz_xy_0 0_xyyz_yy_xz_0 0_xyyz_x_yz_0 0_xyyz_x_yy_0 0_xyyz_x_yyz_0 0_xyyz_xz_yy_0 0_xyyz_xy_yz_0 0_xyyy_y_xy_0 0_xyyy_y_xyy_0 0_xyyy_yy_xy_0 0_xyyy_x_yy_0 0_xyyy_x_yyy_0 0_xyyy_xy_yy_0 0_xxzz_z_xx_0 0_xxzz_z_xxz_0 0_xxzz_zz_xx_0 0_xxzz_x_zz_0 0_xxzz_x_xz_0 0_xxzz_x_xzz_0 0_xxzz_xz_xz_0 0_xxzz_xx_zz_0 0_xxyz_y_xx_0 0_xxyz_y_xxz_0 0_xxyz_yz_xx_0 0_xxyz_x_yz_0 0_xxyz_x_xz_0 0_xxyz_x_xy_0 0_xxyz_x_xyz_0 0_xxyz_xz_xy_0 0_xxyz_xy_xz_0 0_xxyz_xx_yz_0 0_xxyy_y_xx_0 0_xxyy_y_xxy_0 0_xxyy_yy_xx_0 0_xxyy_x_yy_0 0_xxyy_x_xy_0 0_xxyy_x_xyy_0 0_xxyy_xy_xy_0 0_xxyy_xx_yy_0 0_xxxz_x_xz_0 0_xxxz_x_xx_0 0_xxxz_x_xxz_0 0_xxxz_xz_xx_0 0_xxxz_xx_xz_0 0_xxxy_x_xy_0 0_xxxy_x_xx_0 0_xxxy_x_xxy_0 0_xxxy_xy_xx_0 0_xxxy_xx_xy_0 0_xxxx_x_xx_0 0_xxxx_x_xxx_0 0_xxxx_xx_xx_0 

signature:
0_zzzz_0_zzz_0 0_zzzz_0_zzzz_0 0_zzzz_z_zzz_0 0_yzzz_0_zzz_0 0_yzzz_0_yzz_0 0_yzzz_0_yzzz_0 0_yzzz_z_yzz_0 0_yzzz_y_zzz_0 0_yyzz_0_yzz_0 0_yyzz_0_yyz_0 0_yyzz_0_yyzz_0 0_yyzz_z_yyz_0 0_yyzz_y_yzz_0 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyz_0_yyyz_0 0_yyyz_y_yyz_0 0_yyyy_0_yyy_0 0_yyyy_0_yyyy_0 0_yyyy_y_yyy_0 0_xzzz_0_zzz_0 0_xzzz_0_xzz_0 0_xzzz_0_xzzz_0 0_xzzz_z_xzz_0 0_xzzz_x_zzz_0 0_xyzz_0_yzz_0 0_xyzz_0_xzz_0 0_xyzz_0_xyz_0 0_xyzz_0_xyzz_0 0_xyzz_z_xyz_0 0_xyzz_y_xzz_0 0_xyzz_x_yzz_0 0_xyyz_0_yyz_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyz_0_xyyz_0 0_xyyz_y_xyz_0 0_xyyz_x_yyz_0 0_xyyy_0_yyy_0 0_xyyy_0_xyy_0 0_xyyy_0_xyyy_0 0_xyyy_y_xyy_0 0_xyyy_x_yyy_0 0_xxzz_0_xzz_0 0_xxzz_0_xxz_0 0_xxzz_0_xxzz_0 0_xxzz_z_xxz_0 0_xxzz_x_xzz_0 0_xxyz_0_xyz_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyz_0_xxyz_0 0_xxyz_y_xxz_0 0_xxyz_x_xyz_0 0_xxyy_0_xyy_0 0_xxyy_0_xxy_0 0_xxyy_0_xxyy_0 0_xxyy_y_xxy_0 0_xxyy_x_xyy_0 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxz_0_xxxz_0 0_xxxz_x_xxz_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxy_0_xxxy_0 0_xxxy_x_xxy_0 0_xxxx_0_xxx_0 0_xxxx_0_xxxx_0 0_xxxx_x_xxx_0 

signature:
0_zzzz_0_zz_0 0_zzzz_0_zzz_0 0_zzzz_z_zz_0 0_yzzz_0_zz_0 0_yzzz_0_zzz_0 0_yzzz_0_yz_0 0_yzzz_0_yzz_0 0_yzzz_z_yz_0 0_yzzz_y_zz_0 0_yyzz_0_zz_0 0_yyzz_0_yz_0 0_yyzz_0_yzz_0 0_yyzz_0_yy_0 0_yyzz_0_yyz_0 0_yyzz_z_yy_0 0_yyzz_y_zz_0 0_yyzz_y_yz_0 0_yyyz_0_yz_0 0_yyyz_0_yy_0 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyz_y_yz_0 0_yyyz_y_yy_0 0_yyyy_0_yy_0 0_yyyy_0_yyy_0 0_yyyy_y_yy_0 0_xzzz_0_zz_0 0_xzzz_0_zzz_0 0_xzzz_0_xz_0 0_xzzz_0_xzz_0 0_xzzz_z_xz_0 0_xzzz_x_zz_0 0_xyzz_0_zz_0 0_xyzz_0_yz_0 0_xyzz_0_yzz_0 0_xyzz_0_xz_0 0_xyzz_0_xzz_0 0_xyzz_0_xy_0 0_xyzz_0_xyz_0 0_xyzz_z_xy_0 0_xyzz_y_xz_0 0_xyzz_x_zz_0 0_xyzz_x_yz_0 0_xyyz_0_yz_0 0_xyyz_0_yy_0 0_xyyz_0_yyz_0 0_xyyz_0_xz_0 0_xyyz_0_xy_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyz_y_xz_0 0_xyyz_y_xy_0 0_xyyz_x_yz_0 0_xyyz_x_yy_0 0_xyyy_0_yy_0 0_xyyy_0_yyy_0 0_xyyy_0_xy_0 0_xyyy_0_xyy_0 0_xyyy_y_xy_0 0_xyyy_x_yy_0 0_xxzz_0_zz_0 0_xxzz_0_xz_0 0_xxzz_0_xzz_0 0_xxzz_0_xx_0 0_xxzz_0_xxz_0 0_xxzz_z_xx_0 0_xxzz_x_zz_0 0_xxzz_x_xz_0 0_xxyz_0_yz_0 0_xxyz_0_xz_0 0_xxyz_0_xy_0 0_xxyz_0_xyz_0 0_xxyz_0_xx_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyz_y_xx_0 0_xxyz_x_yz_0 0_xxyz_x_xz_0 0_xxyz_x_xy_0 0_xxyy_0_yy_0 0_xxyy_0_xy_0 0_xxyy_0_xyy_0 0_xxyy_0_xx_0 0_xxyy_0_xxy_0 0_xxyy_y_xx_0 0_xxyy_x_yy_0 0_xxyy_x_xy_0 0_xxxz_0_xz_0 0_xxxz_0_xx_0 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxz_x_xz_0 0_xxxz_x_xx_0 0_xxxy_0_xy_0 0_xxxy_0_xx_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxy_x_xy_0 0_xxxy_x_xx_0 0_xxxx_0_xx_0 0_xxxx_0_xxx_0 0_xxxx_x_xx_0 

signature:
0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yzzz_0 0_zz_0_yzzz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_xzzz_0 0_zz_0_xzzz_1 0_zz_0_xyzz_0 0_zz_0_xyzz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zzz_0_zzz_1 0_zzz_0_zzzz_0 0_zzz_0_zzzz_1 0_zzz_0_yzz_1 0_zzz_0_yzzz_0 0_zzz_0_yzzz_1 0_zzz_0_xzz_1 0_zzz_0_xzzz_0 0_zzz_0_xzzz_1 0_zzzz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyzz_1 0_yz_0_xyyz_0 0_yz_0_xyyz_1 0_yz_0_xxyz_0 0_yz_0_xxyz_1 0_yzz_0_yzz_1 0_yzz_0_yzzz_0 0_yzz_0_yyz_1 0_yzz_0_yyzz_0 0_yzz_0_yyzz_1 0_yzz_0_xyz_1 0_yzz_0_xyzz_0 0_yzz_0_xyzz_1 0_yzzz_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyzz_1 0_yy_0_yyyz_0 0_yy_0_yyyz_1 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xyyz_0 0_yy_0_xyyz_1 0_yy_0_xyyy_0 0_yy_0_xyyy_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yyz_0_yyz_1 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyz_1 0_yyz_0_xyyz_0 0_yyz_0_xyyz_1 0_yyzz_0_yyzz_0 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_yyyz_0 0_yyy_0_yyyz_1 0_yyy_0_yyyy_0 0_yyy_0_yyyy_1 0_yyy_0_xyy_1 0_yyy_0_xyyy_0 0_yyy_0_xyyy_1 0_yyyz_0_yyyz_0 0_yyyy_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xzz_0_xzz_1 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxz_1 0_xzz_0_xxzz_0 0_xzz_0_xxzz_1 0_xzzz_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyzz_0_xyzz_0 0_xyy_0_xyy_1 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxy_1 0_xyy_0_xxyy_0 0_xyy_0_xxyy_1 0_xyyz_0_xyyz_0 0_xyyy_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxzz_1 0_xx_0_xxyz_0 0_xx_0_xxyz_1 0_xx_0_xxyy_0 0_xx_0_xxyy_1 0_xx_0_xxxz_0 0_xx_0_xxxz_1 0_xx_0_xxxy_0 0_xx_0_xxxy_1 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxyz_1 0_xxz_0_xxxz_0 0_xxzz_0_xxzz_0 0_xxy_0_xxy_1 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxyz_0_xxyz_0 0_xxyy_0_xxyy_0 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 0_xxx_0_xxxz_0 0_xxx_0_xxxz_1 0_xxx_0_xxxy_0 0_xxx_0_xxxy_1 0_xxx_0_xxxx_0 0_xxx_0_xxxx_1 0_xxxz_0_xxxz_0 0_xxxy_0_xxxy_0 0_xxxx_0_xxxx_0 

signature:
0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zzz_0_zz_1 0_zzz_0_zzz_0 0_zzz_0_zzz_1 0_zzz_0_yz_1 0_zzz_0_yzz_0 0_zzz_0_yzz_1 0_zzz_0_xz_1 0_zzz_0_xzz_0 0_zzz_0_xzz_1 0_zzzz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yzz_0_zz_1 0_yzz_0_zzz_0 0_yzz_0_yz_1 0_yzz_0_yzz_0 0_yzz_0_yzz_1 0_yzz_0_yyz_0 0_yzz_0_yyz_1 0_yzz_0_xzz_0 0_yzz_0_xy_1 0_yzz_0_xyz_0 0_yzz_0_xyz_1 0_yzzz_0_zzz_0 0_yzzz_0_yzz_0 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yyz_0_yz_1 0_yyz_0_yzz_0 0_yyz_0_yy_1 0_yyz_0_yyz_0 0_yyz_0_yyz_1 0_yyz_0_yyy_0 0_yyz_0_xz_1 0_yyz_0_xyz_0 0_yyz_0_xyz_1 0_yyz_0_xyy_0 0_yyzz_0_yzz_0 0_yyzz_0_yyz_0 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_yyz_0 0_yyy_0_yyz_1 0_yyy_0_yyy_0 0_yyy_0_yyy_1 0_yyy_0_xy_1 0_yyy_0_xyy_0 0_yyy_0_xyy_1 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyy_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xzz_0_zz_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xz_1 0_xzz_0_xzz_0 0_xzz_0_xzz_1 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xzz_0_xxz_1 0_xzzz_0_zzz_0 0_xzzz_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyzz_0_yzz_0 0_xyzz_0_xzz_0 0_xyzz_0_xyz_0 0_xyy_0_yy_1 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xy_1 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xyy_1 0_xyy_0_xxy_0 0_xyy_0_xxy_1 0_xyyz_0_yyz_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyy_0_yyy_0 0_xyyy_0_xyy_0 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyz_1 0_xxz_0_xx_1 0_xxz_0_xxz_0 0_xxz_0_xxz_1 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxzz_0_xzz_0 0_xxzz_0_xxz_0 0_xxy_0_xy_1 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xx_1 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxy_1 0_xxy_0_xxx_0 0_xxyz_0_xyz_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyy_0_xyy_0 0_xxyy_0_xxy_0 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 0_xxx_0_xxz_0 0_xxx_0_xxz_1 0_xxx_0_xxy_0 0_xxx_0_xxy_1 0_xxx_0_xxx_0 0_xxx_0_xxx_1 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxx_0_xxx_0 

signature:
0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zzz_0_z_1 0_zzz_0_zz_0 0_zzz_0_zz_1 0_zzz_0_yz_0 0_zzz_0_yz_1 0_zzz_0_xz_0 0_zzz_0_xz_1 0_zzzz_0_zz_0 0_yz_0_yz_0 0_yz_0_yz_1 0_yzz_0_z_1 0_yzz_0_zz_0 0_yzz_0_zz_1 0_yzz_0_y_1 0_yzz_0_yz_0 0_yzz_0_yz_1 0_yzz_0_xy_0 0_yzz_0_xy_1 0_yzzz_0_zz_0 0_yzzz_0_yz_0 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yyz_0_z_1 0_yyz_0_yz_0 0_yyz_0_yz_1 0_yyz_0_yy_0 0_yyz_0_yy_1 0_yyz_0_xz_0 0_yyz_0_xz_1 0_yyzz_0_zz_0 0_yyzz_0_yz_0 0_yyzz_0_yy_0 0_yyy_0_y_1 0_yyy_0_yz_0 0_yyy_0_yz_1 0_yyy_0_yy_0 0_yyy_0_yy_1 0_yyy_0_xy_0 0_yyy_0_xy_1 0_yyyz_0_yz_0 0_yyyz_0_yy_0 0_yyyy_0_yy_0 0_xz_0_xz_0 0_xz_0_xz_1 0_xzz_0_z_1 0_xzz_0_zz_0 0_xzz_0_zz_1 0_xzz_0_xz_0 0_xzz_0_xz_1 0_xzzz_0_zz_0 0_xzzz_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyzz_0_zz_0 0_xyzz_0_yz_0 0_xyzz_0_xz_0 0_xyzz_0_xy_0 0_xyy_0_y_1 0_xyy_0_yy_0 0_xyy_0_yy_1 0_xyy_0_xy_0 0_xyy_0_xy_1 0_xyyz_0_yz_0 0_xyyz_0_yy_0 0_xyyz_0_xz_0 0_xyyz_0_xy_0 0_xyyy_0_yy_0 0_xyyy_0_xy_0 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_z_1 0_xxz_0_yz_0 0_xxz_0_yz_1 0_xxz_0_xz_0 0_xxz_0_xz_1 0_xxz_0_xx_0 0_xxz_0_xx_1 0_xxzz_0_zz_0 0_xxzz_0_xz_0 0_xxzz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xy_1 0_xxy_0_xx_0 0_xxy_0_xx_1 0_xxyz_0_yz_0 0_xxyz_0_xz_0 0_xxyz_0_xy_0 0_xxyz_0_xx_0 0_xxyy_0_yy_0 0_xxyy_0_xy_0 0_xxyy_0_xx_0 0_xxx_0_x_1 0_xxx_0_xz_0 0_xxx_0_xz_1 0_xxx_0_xy_0 0_xxx_0_xy_1 0_xxx_0_xx_0 0_xxx_0_xx_1 0_xxxz_0_xz_0 0_xxxz_0_xx_0 0_xxxy_0_xy_0 0_xxxy_0_xx_0 0_xxxx_0_xx_0 

signature:
0_zzz_z_zz_0 0_zzz_z_zzz_0 0_zzz_zz_zz_0 0_zzz_y_zz_0 0_zzz_y_zzz_0 0_zzz_yz_zz_0 0_zzz_x_zz_0 0_zzz_x_zzz_0 0_zzz_xz_zz_0 0_yzz_z_yz_0 0_yzz_z_yzz_0 0_yzz_zz_yz_0 0_yzz_y_zz_0 0_yzz_y_zzz_0 0_yzz_y_yz_0 0_yzz_y_yzz_0 0_yzz_yz_zz_0 0_yzz_yz_yz_0 0_yzz_yy_zz_0 0_yzz_x_zz_0 0_yzz_x_yz_0 0_yzz_x_yzz_0 0_yzz_xz_yz_0 0_yzz_xy_zz_0 0_yyz_z_yy_0 0_yyz_z_yyz_0 0_yyz_zz_yy_0 0_yyz_y_yz_0 0_yyz_y_yzz_0 0_yyz_y_yy_0 0_yyz_y_yyz_0 0_yyz_yz_yz_0 0_yyz_yz_yy_0 0_yyz_yy_yz_0 0_yyz_x_yz_0 0_yyz_x_yy_0 0_yyz_x_yyz_0 0_yyz_xz_yy_0 0_yyz_xy_yz_0 0_yyy_y_yy_0 0_yyy_y_yyz_0 0_yyy_y_yyy_0 0_yyy_yz_yy_0 0_yyy_yy_yy_0 0_yyy_x_yy_0 0_yyy_x_yyy_0 0_yyy_xy_yy_0 0_xzz_z_xz_0 0_xzz_z_xzz_0 0_xzz_zz_xz_0 0_xzz_y_xz_0 0_xzz_y_xzz_0 0_xzz_yz_xz_0 0_xzz_x_zz_0 0_xzz_x_zzz_0 0_xzz_x_yzz_0 0_xzz_x_xz_0 0_xzz_x_xzz_0 0_xzz_xz_zz_0 0_xzz_xz_xz_0 0_xzz_xy_zz_0 0_xzz_xx_zz_0 0_xyz_z_xy_0 0_xyz_z_xyz_0 0_xyz_zz_xy_0 0_xyz_y_xz_0 0_xyz_y_xzz_0 0_xyz_y_xy_0 0_xyz_y_xyz_0 0_xyz_yz_xz_0 0_xyz_yz_xy_0 0_xyz_yy_xz_0 0_xyz_x_yz_0 0_xyz_x_yzz_0 0_xyz_x_yyz_0 0_xyz_x_xz_0 0_xyz_x_xy_0 0_xyz_x_xyz_0 0_xyz_xz_yz_0 0_xyz_xz_xy_0 0_xyz_xy_yz_0 0_xyz_xy_xz_0 0_xyz_xx_yz_0 0_xyy_y_xy_0 0_xyy_y_xyz_0 0_xyy_y_xyy_0 0_xyy_yz_xy_0 0_xyy_yy_xy_0 0_xyy_x_yy_0 0_xyy_x_yyz_0 0_xyy_x_yyy_0 0_xyy_x_xy_0 0_xyy_x_xyy_0 0_xyy_xz_yy_0 0_xyy_xy_yy_0 0_xyy_xy_xy_0 0_xyy_xx_yy_0 0_xxz_z_xx_0 0_xxz_z_xxz_0 0_xxz_zz_xx_0 0_xxz_y_xx_0 0_xxz_y_xxz_0 0_xxz_yz_xx_0 0_xxz_x_xz_0 0_xxz_x_xzz_0 0_xxz_x_xyz_0 0_xxz_x_xx_0 0_xxz_x_xxz_0 0_xxz_xz_xz_0 0_xxz_xz_xx_0 0_xxz_xy_xz_0 0_xxz_xx_xz_0 0_xxy_y_xx_0 0_xxy_y_xxz_0 0_xxy_y_xxy_0 0_xxy_yz_xx_0 0_xxy_yy_xx_0 0_xxy_x_xy_0 0_xxy_x_xyz_0 0_xxy_x_xyy_0 0_xxy_x_xx_0 0_xxy_x_xxy_0 0_xxy_xz_xy_0 0_xxy_xy_xy_0 0_xxy_xy_xx_0 0_xxy_xx_xy_0 0_xxx_x_xx_0 0_xxx_x_xxz_0 0_xxx_x_xxy_0 0_xxx_x_xxx_0 0_xxx_xz_xx_0 0_xxx_xy_xx_0 0_xxx_xx_xx_0 

signature:
0_zzz_0_zzz_0 0_zzz_0_zzzz_0 0_zzz_0_yzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzz_0 0_zzz_0_xzzz_0 0_zzz_z_zzz_0 0_zzz_y_zzz_0 0_zzz_x_zzz_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyz_0 0_yzz_0_yyzz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyzz_0 0_yzz_z_yzz_0 0_yzz_y_zzz_0 0_yzz_y_yzz_0 0_yzz_x_yzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyzz_0 0_yyz_0_yyy_0 0_yyz_0_yyyz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xyyz_0 0_yyz_z_yyz_0 0_yyz_y_yzz_0 0_yyz_y_yyz_0 0_yyz_x_yyz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyy_0 0_yyy_0_xyyy_0 0_yyy_y_yyz_0 0_yyy_y_yyy_0 0_yyy_x_yyy_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xzzz_0 0_xzz_0_xyz_0 0_xzz_0_xyzz_0 0_xzz_0_xxz_0 0_xzz_0_xxzz_0 0_xzz_z_xzz_0 0_xzz_y_xzz_0 0_xzz_x_zzz_0 0_xzz_x_yzz_0 0_xzz_x_xzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyzz_0 0_xyz_0_xyy_0 0_xyz_0_xyyz_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxyz_0 0_xyz_z_xyz_0 0_xyz_y_xzz_0 0_xyz_y_xyz_0 0_xyz_x_yzz_0 0_xyz_x_yyz_0 0_xyz_x_xyz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxy_0 0_xyy_0_xxyy_0 0_xyy_y_xyz_0 0_xyy_y_xyy_0 0_xyy_x_yyz_0 0_xyy_x_yyy_0 0_xyy_x_xyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxzz_0 0_xxz_0_xxy_0 0_xxz_0_xxyz_0 0_xxz_0_xxx_0 0_xxz_0_xxxz_0 0_xxz_z_xxz_0 0_xxz_y_xxz_0 0_xxz_x_xzz_0 0_xxz_x_xyz_0 0_xxz_x_xxz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxx_0 0_xxy_0_xxxy_0 0_xxy_y_xxz_0 0_xxy_y_xxy_0 0_xxy_x_xyz_0 0_xxy_x_xyy_0 0_xxy_x_xxy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 0_xxx_x_xxz_0 0_xxx_x_xxy_0 0_xxx_x_xxx_0 

signature:
0_zzz_0_zz_0 0_zzz_0_zzz_0 0_zzz_0_yz_0 0_zzz_0_yzz_0 0_zzz_0_xz_0 0_zzz_0_xzz_0 0_zzz_z_zz_0 0_zzz_y_zz_0 0_zzz_x_zz_0 0_yzz_0_zz_0 0_yzz_0_zzz_0 0_yzz_0_yz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xy_0 0_yzz_0_xyz_0 0_yzz_z_yz_0 0_yzz_y_zz_0 0_yzz_y_yz_0 0_yzz_x_zz_0 0_yzz_x_yz_0 0_yyz_0_yz_0 0_yyz_0_yzz_0 0_yyz_0_yy_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_z_yy_0 0_yyz_y_yz_0 0_yyz_y_yy_0 0_yyz_x_yz_0 0_yyz_x_yy_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xy_0 0_yyy_0_xyy_0 0_yyy_y_yy_0 0_yyy_x_yy_0 0_xzz_0_zz_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xzz_z_xz_0 0_xzz_y_xz_0 0_xzz_x_zz_0 0_xzz_x_xz_0 0_xyz_0_yz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xz_0 0_xyz_0_xzz_0 0_xyz_0_xy_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_z_xy_0 0_xyz_y_xz_0 0_xyz_y_xy_0 0_xyz_x_yz_0 0_xyz_x_xz_0 0_xyz_x_xy_0 0_xyy_0_yy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xyy_y_xy_0 0_xyy_x_yy_0 0_xyy_x_xy_0 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xx_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxz_z_xx_0 0_xxz_y_xx_0 0_xxz_x_xz_0 0_xxz_x_xx_0 0_xxy_0_xy_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xx_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxy_y_xx_0 0_xxy_x_xy_0 0_xxy_x_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 0_xxx_x_xx_0 

signature:
0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxxz_0 0_zz_0_zzz_1 0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yzz_1 0_zz_0_yzzz_0 0_zz_0_yzzz_1 0_zz_0_yyz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_xzz_1 0_zz_0_xzzz_0 0_zz_0_xzzz_1 0_zz_0_xyz_1 0_zz_0_xyzz_0 0_zz_0_xyzz_1 0_zz_0_xxz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxy_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yy_0_yyz_1 0_yy_0_yyzz_0 0_yy_0_yyy_1 0_yy_0_yyyz_0 0_yy_0_yyyz_1 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xyyz_0 0_yy_0_xyyz_1 0_yy_0_xyyy_0 0_yy_0_xyyy_1 0_yy_0_xxy_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xzz_0_xxzz_0 0_xyy_0_xxyy_0 0_xx_0_xyz_1 0_xx_0_xxz_1 0_xx_0_xxzz_0 0_xx_0_xxy_1 0_xx_0_xxyz_0 0_xx_0_xxyz_1 0_xx_0_xxyy_0 0_xx_0_xxx_1 0_xx_0_xxxz_0 0_xx_0_xxxz_1 0_xx_0_xxxy_0 0_xx_0_xxxy_1 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_xxyz_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

signature:
0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_zz_0_zzz_1 0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yzz_1 0_zz_0_yzzz_0 0_zz_0_yzzz_1 0_zz_0_yyz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_xzz_1 0_zz_0_xzzz_0 0_zz_0_xzzz_1 0_zz_0_xyz_1 0_zz_0_xyzz_0 0_zz_0_xyzz_1 0_zz_0_xxz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_yz_0_yzz_1 0_yz_0_yzzz_0 0_yz_0_yyz_1 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyz_1 0_yz_0_xyzz_0 0_yz_0_xyzz_1 0_yz_0_xyyz_0 0_yz_0_xyyz_1 0_yz_0_xxyz_0 0_yz_0_xxyz_1 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyzz_0 0_yy_0_yyzz_1 0_yy_0_yyy_1 0_yy_0_yyyz_0 0_yy_0_yyyz_1 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xyyz_0 0_yy_0_xyyz_1 0_yy_0_xyyy_0 0_yy_0_xyyy_1 0_yy_0_xxy_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_xzz_1 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxz_1 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxzz_0 0_xy_0_xyzz_0 0_xy_0_xyy_1 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxy_1 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxyy_0 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxzz_0 0_xx_0_xxzz_1 0_xx_0_xxy_1 0_xx_0_xxyz_0 0_xx_0_xxyz_1 0_xx_0_xxyy_0 0_xx_0_xxyy_1 0_xx_0_xxx_1 0_xx_0_xxxz_0 0_xx_0_xxxz_1 0_xx_0_xxxy_0 0_xx_0_xxxy_1 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxxz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

signature:
0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xy_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xyz_0 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yyz_0_yyz_0 0_yyz_0_xyz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xxz_0 0_xzz_0_xzz_0 0_xzz_0_xxz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxy_0_xxy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

signature:
0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xy_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yz_0_yz_1 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xz_1 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xy_0_xy_1 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

signature:
0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_yz_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xy_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

signature:
0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xz_1 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

signature:
0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_0_1 0_zz_0_z_0 0_zz_0_z_1 0_zz_0_y_0 0_zz_0_y_1 0_zz_0_x_0 0_zz_0_x_1 0_zzz_0_z_0 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yzz_0_z_0 0_yzz_0_y_0 0_yy_0_0_1 0_yy_0_z_0 0_yy_0_z_1 0_yy_0_y_0 0_yy_0_y_1 0_yy_0_x_0 0_yy_0_x_1 0_yyz_0_z_0 0_yyy_0_y_0 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xzz_0_z_0 0_xyy_0_y_0 0_xx_0_0_1 0_xx_0_z_0 0_xx_0_z_1 0_xx_0_y_0 0_xx_0_y_1 0_xx_0_x_0 0_xx_0_x_1 0_xxz_0_z_0 0_xxx_0_x_0 

signature:
0_zz_z_zz_0 0_zz_z_zzz_0 0_zz_zz_zz_0 0_zz_y_zz_0 0_zz_y_zzz_0 0_zz_y_yzz_0 0_zz_yz_zz_0 0_zz_yy_zz_0 0_zz_x_zz_0 0_zz_x_zzz_0 0_zz_x_yzz_0 0_zz_x_xzz_0 0_zz_xz_zz_0 0_zz_xy_zz_0 0_zz_xx_zz_0 0_yz_z_yz_0 0_yz_z_yzz_0 0_yz_zz_yz_0 0_yz_y_yz_0 0_yz_y_yzz_0 0_yz_y_yyz_0 0_yz_yz_yz_0 0_yz_yy_yz_0 0_yz_x_yz_0 0_yz_x_yzz_0 0_yz_x_yyz_0 0_yz_x_xyz_0 0_yz_xz_yz_0 0_yz_xy_yz_0 0_yz_xx_yz_0 0_yy_z_yy_0 0_yy_z_yyz_0 0_yy_zz_yy_0 0_yy_y_yy_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_yz_yy_0 0_yy_yy_yy_0 0_yy_x_yy_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xyy_0 0_yy_xz_yy_0 0_yy_xy_yy_0 0_yy_xx_yy_0 0_xz_z_xz_0 0_xz_z_xzz_0 0_xz_zz_xz_0 0_xz_y_xz_0 0_xz_y_xzz_0 0_xz_y_xyz_0 0_xz_yz_xz_0 0_xz_yy_xz_0 0_xz_x_xz_0 0_xz_x_xzz_0 0_xz_x_xyz_0 0_xz_x_xxz_0 0_xz_xz_xz_0 0_xz_xy_xz_0 0_xz_xx_xz_0 0_xy_z_xy_0 0_xy_z_xyz_0 0_xy_zz_xy_0 0_xy_y_xy_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_yz_xy_0 0_xy_yy_xy_0 0_xy_x_xy_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xxy_0 0_xy_xz_xy_0 0_xy_xy_xy_0 0_xy_xx_xy_0 0_xx_z_xx_0 0_xx_z_xxz_0 0_xx_zz_xx_0 0_xx_y_xx_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_yz_xx_0 0_xx_yy_xx_0 0_xx_x_xx_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 0_xx_xz_xx_0 0_xx_xy_xx_0 0_xx_xx_xx_0 

signature:
0_zz_0_zzz_0 0_zz_0_zzzz_0 0_zz_0_yzz_0 0_zz_0_yzzz_0 0_zz_0_yyz_0 0_zz_0_yyzz_0 0_zz_0_xzz_0 0_zz_0_xzzz_0 0_zz_0_xyz_0 0_zz_0_xyzz_0 0_zz_0_xxz_0 0_zz_0_xxzz_0 0_zz_z_zzz_0 0_zz_y_zzz_0 0_zz_y_yzz_0 0_zz_x_zzz_0 0_zz_x_yzz_0 0_zz_x_xzz_0 0_yz_0_yzz_0 0_yz_0_yzzz_0 0_yz_0_yyz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yz_z_yzz_0 0_yz_y_yzz_0 0_yz_y_yyz_0 0_yz_x_yzz_0 0_yz_x_yyz_0 0_yz_x_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyzz_0 0_yy_0_yyy_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxy_0 0_yy_0_xxyy_0 0_yy_z_yyz_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xyy_0 0_xz_0_xzz_0 0_xz_0_xzzz_0 0_xz_0_xyz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xz_z_xzz_0 0_xz_y_xzz_0 0_xz_y_xyz_0 0_xz_x_xzz_0 0_xz_x_xyz_0 0_xz_x_xxz_0 0_xy_0_xyz_0 0_xy_0_xyzz_0 0_xy_0_xyy_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xy_z_xyz_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxzz_0 0_xx_0_xxy_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxx_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 0_xx_z_xxz_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 

signature:
0_zz_0_zz_0 0_zz_0_zzz_0 0_zz_0_yz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xz_0 0_zz_0_xzz_0 0_zz_0_xy_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_zz_z_zz_0 0_zz_y_zz_0 0_zz_x_zz_0 0_yz_0_yz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yz_z_yz_0 0_yz_y_yz_0 0_yz_x_yz_0 0_yy_0_yz_0 0_yy_0_yzz_0 0_yy_0_yy_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_yy_z_yy_0 0_yy_y_yy_0 0_yy_x_yy_0 0_xz_0_xz_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xz_z_xz_0 0_xz_y_xz_0 0_xz_x_xz_0 0_xy_0_xy_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xy_z_xy_0 0_xy_y_xy_0 0_xy_x_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xzz_0 0_xx_0_xy_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xx_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 0_xx_z_xx_0 0_xx_y_xx_0 0_xx_x_xx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxy_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_x_0_xyz_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xx_0_xxyz_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxxz_0 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxy_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_x_0_yz_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xx_0_xyz_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xxz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_y_0_z_1 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_1 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_x_0_z_1 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xx_0_0_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xxzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyy_0 0_x_0_xxyz_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xyz_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

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
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

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

INTEGRAL:0 : 1 : 7 : SSSS_7 Y SSSS_8 Y SSSP_7 Y 

0_0_0_0_7 0_0_0_0_8 0_0_0_z_7 0_0_0_y_7 0_0_0_x_7 

SSSS_7 SSSS_8 SSSP_7 

SSSS_7 : 0_0_0_0_7 

SSSS_8 : 0_0_0_0_8 

SSSP_7 : 0_0_0_z_7 0_0_0_y_7 0_0_0_x_7 

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

0_0_0_0_4 0_0_0_0_5 0_0_0_z_4 0_0_0_z_5 0_0_0_zz_4 0_0_0_y_4 0_0_0_y_5 0_0_0_yz_4 0_0_0_yy_4 0_0_0_x_4 0_0_0_x_5 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSS_4 SSSS_5 SSSP_4 SSSP_5 SSSD_4 

SSSS_4 : 0_0_0_0_4 

SSSS_5 : 0_0_0_0_5 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

INTEGRAL:0 : 2 : 5 : SSSS_5 Y SSSS_6 Y SSSP_5 Y SSSP_6 Y SSSD_5 Y 

0_0_0_0_5 0_0_0_0_6 0_0_0_z_5 0_0_0_z_6 0_0_0_zz_5 0_0_0_y_5 0_0_0_y_6 0_0_0_yz_5 0_0_0_yy_5 0_0_0_x_5 0_0_0_x_6 0_0_0_xx_5 

SSSS_5 SSSS_6 SSSP_5 SSSP_6 SSSD_5 

SSSS_5 : 0_0_0_0_5 

SSSS_6 : 0_0_0_0_6 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSP_6 : 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

SSSD_5 : 0_0_0_zz_5 0_0_0_yz_5 0_0_0_yy_5 0_0_0_xx_5 

INTEGRAL:0 : 2 : 6 : SSSS_6 Y SSSS_7 Y SSSP_6 Y SSSP_7 Y SSSD_6 Y 

0_0_0_0_6 0_0_0_0_7 0_0_0_z_6 0_0_0_z_7 0_0_0_zz_6 0_0_0_y_6 0_0_0_y_7 0_0_0_yy_6 0_0_0_x_6 0_0_0_x_7 0_0_0_xx_6 

SSSS_6 SSSS_7 SSSP_6 SSSP_7 SSSD_6 

SSSS_6 : 0_0_0_0_6 

SSSS_7 : 0_0_0_0_7 

SSSP_6 : 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

SSSP_7 : 0_0_0_z_7 0_0_0_y_7 0_0_0_x_7 

SSSD_6 : 0_0_0_zz_6 0_0_0_yy_6 0_0_0_xx_6 

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

INTEGRAL:0 : 3 : 3 : SSSP_3 Y SSSP_4 Y SSSD_3 N SSSD_4 N SSSF_3 Y 

0_0_0_z_3 0_0_0_z_4 0_0_0_zz_3 0_0_0_zz_4 0_0_0_zzz_3 0_0_0_y_3 0_0_0_y_4 0_0_0_yz_3 0_0_0_yz_4 0_0_0_yzz_3 0_0_0_yy_3 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_x_3 0_0_0_x_4 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xx_3 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSP_3 SSSP_4 SSSD_3 SSSD_4 SSSF_3 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

INTEGRAL:0 : 3 : 4 : SSSP_4 Y SSSP_5 Y SSSD_4 N SSSD_5 Y SSSF_4 Y 

0_0_0_z_4 0_0_0_z_5 0_0_0_zz_4 0_0_0_zz_5 0_0_0_zzz_4 0_0_0_y_4 0_0_0_y_5 0_0_0_yz_4 0_0_0_yz_5 0_0_0_yzz_4 0_0_0_yy_4 0_0_0_yy_5 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_x_4 0_0_0_x_5 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xx_4 0_0_0_xx_5 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

SSSP_4 SSSP_5 SSSD_4 SSSD_5 SSSF_4 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSD_5 : 0_0_0_zz_5 0_0_0_yz_5 0_0_0_yy_5 0_0_0_xx_5 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

INTEGRAL:0 : 3 : 5 : SSSP_5 Y SSSP_6 Y SSSD_5 N SSSD_6 Y SSSF_5 Y 

0_0_0_z_5 0_0_0_z_6 0_0_0_zz_5 0_0_0_zz_6 0_0_0_zzz_5 0_0_0_y_5 0_0_0_y_6 0_0_0_yzz_5 0_0_0_yy_5 0_0_0_yy_6 0_0_0_yyz_5 0_0_0_yyy_5 0_0_0_x_5 0_0_0_x_6 0_0_0_xzz_5 0_0_0_xyy_5 0_0_0_xx_5 0_0_0_xx_6 0_0_0_xxz_5 0_0_0_xxx_5 

SSSP_5 SSSP_6 SSSD_5 SSSD_6 SSSF_5 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSP_6 : 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

SSSD_5 : 0_0_0_zz_5 0_0_0_yz_5 0_0_0_yy_5 0_0_0_xx_5 

SSSD_6 : 0_0_0_zz_6 0_0_0_yy_6 0_0_0_xx_6 

SSSF_5 : 0_0_0_zzz_5 0_0_0_yzz_5 0_0_0_yyz_5 0_0_0_yyy_5 0_0_0_xzz_5 0_0_0_xyy_5 0_0_0_xxz_5 0_0_0_xxx_5 

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

INTEGRAL:0 : 4 : 3 : SSSD_3 N SSSD_4 N SSSF_3 N SSSF_4 N SSSG_3 Y 

0_0_0_zz_3 0_0_0_zz_4 0_0_0_zzz_3 0_0_0_zzz_4 0_0_0_zzzz_3 0_0_0_yzz_3 0_0_0_yzz_4 0_0_0_yzzz_3 0_0_0_yy_3 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyz_4 0_0_0_yyzz_3 0_0_0_yyy_3 0_0_0_yyy_4 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzz_3 0_0_0_xzz_4 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyy_3 0_0_0_xyy_4 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xx_3 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxz_4 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxx_3 0_0_0_xxx_4 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SSSD_3 SSSD_4 SSSF_3 SSSF_4 SSSG_3 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

INTEGRAL:0 : 4 : 4 : SSSD_4 N SSSD_5 N SSSF_4 N SSSF_5 Y SSSG_4 Y 

0_0_0_zz_4 0_0_0_zz_5 0_0_0_zzz_4 0_0_0_zzz_5 0_0_0_zzzz_4 0_0_0_yzz_4 0_0_0_yzz_5 0_0_0_yzzz_4 0_0_0_yy_4 0_0_0_yy_5 0_0_0_yyz_4 0_0_0_yyz_5 0_0_0_yyzz_4 0_0_0_yyy_4 0_0_0_yyy_5 0_0_0_yyyz_4 0_0_0_yyyy_4 0_0_0_xzz_4 0_0_0_xzz_5 0_0_0_xzzz_4 0_0_0_xyzz_4 0_0_0_xyy_4 0_0_0_xyy_5 0_0_0_xyyz_4 0_0_0_xyyy_4 0_0_0_xx_4 0_0_0_xx_5 0_0_0_xxz_4 0_0_0_xxz_5 0_0_0_xxzz_4 0_0_0_xxyz_4 0_0_0_xxyy_4 0_0_0_xxx_4 0_0_0_xxx_5 0_0_0_xxxz_4 0_0_0_xxxy_4 0_0_0_xxxx_4 

SSSD_4 SSSD_5 SSSF_4 SSSF_5 SSSG_4 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSD_5 : 0_0_0_zz_5 0_0_0_yz_5 0_0_0_yy_5 0_0_0_xx_5 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

SSSF_5 : 0_0_0_zzz_5 0_0_0_yzz_5 0_0_0_yyz_5 0_0_0_yyy_5 0_0_0_xzz_5 0_0_0_xyy_5 0_0_0_xxz_5 0_0_0_xxx_5 

SSSG_4 : 0_0_0_zzzz_4 0_0_0_yzzz_4 0_0_0_yyzz_4 0_0_0_yyyz_4 0_0_0_yyyy_4 0_0_0_xzzz_4 0_0_0_xyzz_4 0_0_0_xyyz_4 0_0_0_xyyy_4 0_0_0_xxzz_4 0_0_0_xxyz_4 0_0_0_xxyy_4 0_0_0_xxxz_4 0_0_0_xxxy_4 0_0_0_xxxx_4 

INTEGRAL:1 : 0 : 2 : SSSS_2 Y SSSS_3 Y SPSS_2 Y 

0_0_0_0_2 0_0_0_0_3 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SSSS_2 SSSS_3 SPSS_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

INTEGRAL:1 : 0 : 3 : SSSS_3 Y SSSS_4 Y SPSS_3 Y 

0_0_0_0_3 0_0_0_0_4 0_z_0_0_3 0_y_0_0_3 0_x_0_0_3 

SSSS_3 SSSS_4 SPSS_3 

SSSS_3 : 0_0_0_0_3 

SSSS_4 : 0_0_0_0_4 

SPSS_3 : 0_z_0_0_3 0_y_0_0_3 0_x_0_0_3 

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

INTEGRAL:1 : 1 : 3 : SSSS_4 Y SSSP_3 Y SSSP_4 Y SPSP_3 Y 

0_0_0_0_4 0_0_0_z_3 0_0_0_z_4 0_0_0_y_3 0_0_0_y_4 0_0_0_x_3 0_0_0_x_4 0_z_0_z_3 0_z_0_y_3 0_z_0_x_3 0_y_0_z_3 0_y_0_y_3 0_y_0_x_3 0_x_0_z_3 0_x_0_y_3 0_x_0_x_3 

SSSS_4 SSSP_3 SSSP_4 SPSP_3 

SSSS_4 : 0_0_0_0_4 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SPSP_3 : 0_z_0_z_3 0_z_0_y_3 0_z_0_x_3 0_y_0_z_3 0_y_0_y_3 0_y_0_x_3 0_x_0_z_3 0_x_0_y_3 0_x_0_x_3 

INTEGRAL:1 : 2 : 0 : SSSP_1 Y SSSD_0 Y SSSD_1 Y SPSD_0 Y 

0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SSSP_1 SSSD_0 SSSD_1 SPSD_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

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

INTEGRAL:1 : 2 : 3 : SSSP_4 Y SSSD_3 Y SSSD_4 Y SPSD_3 Y 

0_0_0_z_4 0_0_0_zz_3 0_0_0_zz_4 0_0_0_y_4 0_0_0_yz_3 0_0_0_yz_4 0_0_0_yy_3 0_0_0_yy_4 0_0_0_x_4 0_0_0_xz_3 0_0_0_xz_4 0_0_0_xy_3 0_0_0_xy_4 0_0_0_xx_3 0_0_0_xx_4 0_z_0_zz_3 0_z_0_yz_3 0_z_0_yy_3 0_z_0_xz_3 0_z_0_xy_3 0_z_0_xx_3 0_y_0_yz_3 0_y_0_yy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xx_3 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xx_3 

SSSP_4 SSSD_3 SSSD_4 SPSD_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SPSD_3 : 0_z_0_zz_3 0_z_0_yz_3 0_z_0_yy_3 0_z_0_xz_3 0_z_0_xy_3 0_z_0_xx_3 0_y_0_yz_3 0_y_0_yy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xx_3 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xx_3 

INTEGRAL:1 : 3 : 0 : SSSD_1 Y SSSF_0 Y SSSF_1 Y SPSF_0 Y 

0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SSSD_1 SSSF_0 SSSF_1 SPSF_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

INTEGRAL:1 : 3 : 1 : SSSD_2 Y SSSF_1 Y SSSF_2 Y SPSF_1 Y 

0_0_0_zz_2 0_0_0_zzz_1 0_0_0_zzz_2 0_0_0_yz_2 0_0_0_yzz_1 0_0_0_yzz_2 0_0_0_yy_2 0_0_0_yyz_1 0_0_0_yyz_2 0_0_0_yyy_1 0_0_0_yyy_2 0_0_0_xz_2 0_0_0_xzz_1 0_0_0_xzz_2 0_0_0_xy_2 0_0_0_xyz_1 0_0_0_xyz_2 0_0_0_xyy_1 0_0_0_xyy_2 0_0_0_xx_2 0_0_0_xxz_1 0_0_0_xxz_2 0_0_0_xxy_1 0_0_0_xxy_2 0_0_0_xxx_1 0_0_0_xxx_2 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SSSD_2 SSSF_1 SSSF_2 SPSF_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

INTEGRAL:1 : 3 : 2 : SSSD_3 Y SSSF_2 Y SSSF_3 Y SPSF_2 Y 

0_0_0_zz_3 0_0_0_zzz_2 0_0_0_zzz_3 0_0_0_yz_3 0_0_0_yzz_2 0_0_0_yzz_3 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyz_3 0_0_0_yyy_2 0_0_0_yyy_3 0_0_0_xz_3 0_0_0_xzz_2 0_0_0_xzz_3 0_0_0_xy_3 0_0_0_xyz_2 0_0_0_xyz_3 0_0_0_xyy_2 0_0_0_xyy_3 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxz_3 0_0_0_xxy_2 0_0_0_xxy_3 0_0_0_xxx_2 0_0_0_xxx_3 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SSSD_3 SSSF_2 SSSF_3 SPSF_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

INTEGRAL:1 : 3 : 3 : SSSD_4 Y SSSF_3 Y SSSF_4 Y SPSF_3 Y 

0_0_0_zz_4 0_0_0_zzz_3 0_0_0_zzz_4 0_0_0_yz_4 0_0_0_yzz_3 0_0_0_yzz_4 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyz_4 0_0_0_yyy_3 0_0_0_yyy_4 0_0_0_xz_4 0_0_0_xzz_3 0_0_0_xzz_4 0_0_0_xy_4 0_0_0_xyz_3 0_0_0_xyz_4 0_0_0_xyy_3 0_0_0_xyy_4 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxz_4 0_0_0_xxy_3 0_0_0_xxy_4 0_0_0_xxx_3 0_0_0_xxx_4 0_z_0_zzz_3 0_z_0_yzz_3 0_z_0_yyz_3 0_z_0_xzz_3 0_z_0_xyz_3 0_z_0_xxz_3 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xxy_3 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxx_3 

SSSD_4 SSSF_3 SSSF_4 SPSF_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xz_4 0_0_0_xy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

SPSF_3 : 0_z_0_zzz_3 0_z_0_yzz_3 0_z_0_yyz_3 0_z_0_xzz_3 0_z_0_xyz_3 0_z_0_xxz_3 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xxy_3 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxx_3 

INTEGRAL:1 : 4 : 0 : SSSF_1 Y SSSG_0 Y SSSG_1 Y SPSG_0 Y 

0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SSSF_1 SSSG_0 SSSG_1 SPSG_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

INTEGRAL:1 : 4 : 1 : SSSF_2 Y SSSG_1 Y SSSG_2 Y SPSG_1 Y 

0_0_0_zzz_2 0_0_0_zzzz_1 0_0_0_zzzz_2 0_0_0_yzz_2 0_0_0_yzzz_1 0_0_0_yzzz_2 0_0_0_yyz_2 0_0_0_yyzz_1 0_0_0_yyzz_2 0_0_0_yyy_2 0_0_0_yyyz_1 0_0_0_yyyz_2 0_0_0_yyyy_1 0_0_0_yyyy_2 0_0_0_xzz_2 0_0_0_xzzz_1 0_0_0_xzzz_2 0_0_0_xyz_2 0_0_0_xyzz_1 0_0_0_xyzz_2 0_0_0_xyy_2 0_0_0_xyyz_1 0_0_0_xyyz_2 0_0_0_xyyy_1 0_0_0_xyyy_2 0_0_0_xxz_2 0_0_0_xxzz_1 0_0_0_xxzz_2 0_0_0_xxy_2 0_0_0_xxyz_1 0_0_0_xxyz_2 0_0_0_xxyy_1 0_0_0_xxyy_2 0_0_0_xxx_2 0_0_0_xxxz_1 0_0_0_xxxz_2 0_0_0_xxxy_1 0_0_0_xxxy_2 0_0_0_xxxx_1 0_0_0_xxxx_2 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SSSF_2 SSSG_1 SSSG_2 SPSG_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

INTEGRAL:1 : 4 : 2 : SSSF_3 Y SSSG_2 Y SSSG_3 Y SPSG_2 Y 

0_0_0_zzz_3 0_0_0_zzzz_2 0_0_0_zzzz_3 0_0_0_yzz_3 0_0_0_yzzz_2 0_0_0_yzzz_3 0_0_0_yyz_3 0_0_0_yyzz_2 0_0_0_yyzz_3 0_0_0_yyy_3 0_0_0_yyyz_2 0_0_0_yyyz_3 0_0_0_yyyy_2 0_0_0_yyyy_3 0_0_0_xzz_3 0_0_0_xzzz_2 0_0_0_xzzz_3 0_0_0_xyz_3 0_0_0_xyzz_2 0_0_0_xyzz_3 0_0_0_xyy_3 0_0_0_xyyz_2 0_0_0_xyyz_3 0_0_0_xyyy_2 0_0_0_xyyy_3 0_0_0_xxz_3 0_0_0_xxzz_2 0_0_0_xxzz_3 0_0_0_xxy_3 0_0_0_xxyz_2 0_0_0_xxyz_3 0_0_0_xxyy_2 0_0_0_xxyy_3 0_0_0_xxx_3 0_0_0_xxxz_2 0_0_0_xxxz_3 0_0_0_xxxy_2 0_0_0_xxxy_3 0_0_0_xxxx_2 0_0_0_xxxx_3 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SSSF_3 SSSG_2 SSSG_3 SPSG_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

INTEGRAL:1 : 4 : 3 : SSSF_4 Y SSSG_3 Y SSSG_4 Y SPSG_3 Y 

0_0_0_zzz_4 0_0_0_zzzz_3 0_0_0_zzzz_4 0_0_0_yzz_4 0_0_0_yzzz_3 0_0_0_yzzz_4 0_0_0_yyz_4 0_0_0_yyzz_3 0_0_0_yyzz_4 0_0_0_yyy_4 0_0_0_yyyz_3 0_0_0_yyyz_4 0_0_0_yyyy_3 0_0_0_yyyy_4 0_0_0_xzz_4 0_0_0_xzzz_3 0_0_0_xzzz_4 0_0_0_xyz_4 0_0_0_xyzz_3 0_0_0_xyzz_4 0_0_0_xyy_4 0_0_0_xyyz_3 0_0_0_xyyz_4 0_0_0_xyyy_3 0_0_0_xyyy_4 0_0_0_xxz_4 0_0_0_xxzz_3 0_0_0_xxzz_4 0_0_0_xxy_4 0_0_0_xxyz_3 0_0_0_xxyz_4 0_0_0_xxyy_3 0_0_0_xxyy_4 0_0_0_xxx_4 0_0_0_xxxz_3 0_0_0_xxxz_4 0_0_0_xxxy_3 0_0_0_xxxy_4 0_0_0_xxxx_3 0_0_0_xxxx_4 0_z_0_zzzz_3 0_z_0_yzzz_3 0_z_0_yyzz_3 0_z_0_xzzz_3 0_z_0_xyzz_3 0_z_0_xxzz_3 0_y_0_yyyz_3 0_y_0_yyyy_3 0_y_0_xyyz_3 0_y_0_xyyy_3 0_y_0_xxyy_3 0_x_0_xxyz_3 0_x_0_xxxz_3 0_x_0_xxxy_3 0_x_0_xxxx_3 

SSSF_4 SSSG_3 SSSG_4 SPSG_3 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxy_4 0_0_0_xxx_4 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SSSG_4 : 0_0_0_zzzz_4 0_0_0_yzzz_4 0_0_0_yyzz_4 0_0_0_yyyz_4 0_0_0_yyyy_4 0_0_0_xzzz_4 0_0_0_xyzz_4 0_0_0_xyyz_4 0_0_0_xyyy_4 0_0_0_xxzz_4 0_0_0_xxyz_4 0_0_0_xxyy_4 0_0_0_xxxz_4 0_0_0_xxxy_4 0_0_0_xxxx_4 

SPSG_3 : 0_z_0_zzzz_3 0_z_0_yzzz_3 0_z_0_yyzz_3 0_z_0_xzzz_3 0_z_0_xyzz_3 0_z_0_xxzz_3 0_y_0_yyyz_3 0_y_0_yyyy_3 0_y_0_xyyz_3 0_y_0_xyyy_3 0_y_0_xxyy_3 0_x_0_xxyz_3 0_x_0_xxxz_3 0_x_0_xxxy_3 0_x_0_xxxx_3 

INTEGRAL:2 : 0 : 2 : SSSS_2 Y SSSS_3 Y SPSS_2 Y SPSS_3 Y SDSS_2 Y 

0_0_0_0_2 0_0_0_0_3 0_z_0_0_2 0_z_0_0_3 0_zz_0_0_2 0_y_0_0_2 0_y_0_0_3 0_yy_0_0_2 0_x_0_0_2 0_x_0_0_3 0_xx_0_0_2 

SSSS_2 SSSS_3 SPSS_2 SPSS_3 SDSS_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SPSS_3 : 0_z_0_0_3 0_y_0_0_3 0_x_0_0_3 

SDSS_2 : 0_zz_0_0_2 0_yy_0_0_2 0_xx_0_0_2 

INTEGRAL:2 : 1 : 1 : SSSP_1 Y SSSP_2 Y SPSS_2 Y SPSP_1 Y SPSP_2 Y SDSP_1 Y 

0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_0_2 0_z_0_z_1 0_z_0_z_2 0_z_0_y_1 0_z_0_y_2 0_z_0_x_1 0_z_0_x_2 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_y_0_0_2 0_y_0_z_1 0_y_0_z_2 0_y_0_y_1 0_y_0_y_2 0_y_0_x_1 0_y_0_x_2 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_x_0_0_2 0_x_0_z_1 0_x_0_z_2 0_x_0_y_1 0_x_0_y_2 0_x_0_x_1 0_x_0_x_2 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SSSP_1 SSSP_2 SPSS_2 SPSP_1 SPSP_2 SDSP_1 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

INTEGRAL:2 : 1 : 2 : SSSP_2 Y SSSP_3 Y SPSS_3 Y SPSP_2 Y SPSP_3 Y SDSP_2 Y 

0_0_0_z_2 0_0_0_z_3 0_0_0_y_2 0_0_0_y_3 0_0_0_x_2 0_0_0_x_3 0_z_0_0_3 0_z_0_z_2 0_z_0_z_3 0_z_0_y_2 0_z_0_y_3 0_z_0_x_2 0_z_0_x_3 0_zz_0_z_2 0_zz_0_y_2 0_zz_0_x_2 0_y_0_0_3 0_y_0_z_2 0_y_0_z_3 0_y_0_y_2 0_y_0_y_3 0_y_0_x_2 0_y_0_x_3 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_x_2 0_x_0_0_3 0_x_0_z_2 0_x_0_z_3 0_x_0_y_2 0_x_0_y_3 0_x_0_x_2 0_x_0_x_3 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_x_2 

SSSP_2 SSSP_3 SPSS_3 SPSP_2 SPSP_3 SDSP_2 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SPSS_3 : 0_z_0_0_3 0_y_0_0_3 0_x_0_0_3 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SPSP_3 : 0_z_0_z_3 0_z_0_y_3 0_z_0_x_3 0_y_0_z_3 0_y_0_y_3 0_y_0_x_3 0_x_0_z_3 0_x_0_y_3 0_x_0_x_3 

SDSP_2 : 0_zz_0_z_2 0_zz_0_y_2 0_zz_0_x_2 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_x_2 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_x_2 

INTEGRAL:2 : 2 : 0 : SSSD_0 Y SSSD_1 Y SPSP_1 Y SPSD_0 Y SPSD_1 N SDSD_0 Y 

0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_y_0_z_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_x_0_z_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SSSD_0 SSSD_1 SPSP_1 SPSD_0 SPSD_1 SDSD_0 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

INTEGRAL:2 : 2 : 1 : SSSD_1 Y SSSD_2 Y SPSP_2 Y SPSD_1 N SPSD_2 N SDSD_1 Y 

0_0_0_zz_1 0_0_0_zz_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yy_1 0_0_0_yy_2 0_0_0_xz_1 0_0_0_xz_2 0_0_0_xy_1 0_0_0_xy_2 0_0_0_xx_1 0_0_0_xx_2 0_z_0_z_2 0_z_0_zz_1 0_z_0_zz_2 0_z_0_y_2 0_z_0_yz_1 0_z_0_yz_2 0_z_0_x_2 0_z_0_xz_1 0_z_0_xz_2 0_z_0_xy_1 0_z_0_xy_2 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_y_0_z_2 0_y_0_y_2 0_y_0_yz_1 0_y_0_yz_2 0_y_0_yy_1 0_y_0_yy_2 0_y_0_x_2 0_y_0_xz_1 0_y_0_xz_2 0_y_0_xy_1 0_y_0_xy_2 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_x_0_z_2 0_x_0_y_2 0_x_0_yz_1 0_x_0_yz_2 0_x_0_x_2 0_x_0_xz_1 0_x_0_xz_2 0_x_0_xy_1 0_x_0_xy_2 0_x_0_xx_1 0_x_0_xx_2 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SSSD_1 SSSD_2 SPSP_2 SPSD_1 SPSD_2 SDSD_1 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

INTEGRAL:2 : 2 : 2 : SSSD_2 Y SSSD_3 Y SPSP_3 Y SPSD_2 N SPSD_3 N SDSD_2 Y 

0_0_0_zz_2 0_0_0_zz_3 0_0_0_yz_2 0_0_0_yz_3 0_0_0_yy_2 0_0_0_yy_3 0_0_0_xz_2 0_0_0_xz_3 0_0_0_xy_2 0_0_0_xy_3 0_0_0_xx_2 0_0_0_xx_3 0_z_0_z_3 0_z_0_zz_2 0_z_0_zz_3 0_z_0_y_3 0_z_0_yz_2 0_z_0_yz_3 0_z_0_x_3 0_z_0_xz_2 0_z_0_xz_3 0_z_0_xy_2 0_z_0_xy_3 0_zz_0_zz_2 0_zz_0_yz_2 0_zz_0_xz_2 0_zz_0_xy_2 0_y_0_z_3 0_y_0_y_3 0_y_0_yz_2 0_y_0_yz_3 0_y_0_yy_2 0_y_0_yy_3 0_y_0_x_3 0_y_0_xz_2 0_y_0_xz_3 0_y_0_xy_2 0_y_0_xy_3 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_xz_2 0_yy_0_xy_2 0_x_0_z_3 0_x_0_y_3 0_x_0_yz_2 0_x_0_yz_3 0_x_0_x_3 0_x_0_xz_2 0_x_0_xz_3 0_x_0_xy_2 0_x_0_xy_3 0_x_0_xx_2 0_x_0_xx_3 0_xx_0_yz_2 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xx_2 

SSSD_2 SSSD_3 SPSP_3 SPSD_2 SPSD_3 SDSD_2 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SPSP_3 : 0_z_0_z_3 0_z_0_y_3 0_z_0_x_3 0_y_0_z_3 0_y_0_y_3 0_y_0_x_3 0_x_0_z_3 0_x_0_y_3 0_x_0_x_3 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SPSD_3 : 0_z_0_zz_3 0_z_0_yz_3 0_z_0_yy_3 0_z_0_xz_3 0_z_0_xy_3 0_z_0_xx_3 0_y_0_yz_3 0_y_0_yy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xx_3 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xx_3 

SDSD_2 : 0_zz_0_zz_2 0_zz_0_yz_2 0_zz_0_xz_2 0_zz_0_xy_2 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_xz_2 0_yy_0_xy_2 0_xx_0_yz_2 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xx_2 

INTEGRAL:2 : 3 : 0 : SSSF_0 Y SSSF_1 Y SPSD_1 Y SPSF_0 Y SPSF_1 Y SDSF_0 Y 

0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxy_0 0_y_0_xxy_1 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SSSF_0 SSSF_1 SPSD_1 SPSF_0 SPSF_1 SDSF_0 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

INTEGRAL:2 : 3 : 1 : SSSF_1 Y SSSF_2 Y SPSD_2 Y SPSF_1 Y SPSF_2 Y SDSF_1 Y 

0_0_0_zzz_1 0_0_0_zzz_2 0_0_0_yzz_1 0_0_0_yzz_2 0_0_0_yyz_1 0_0_0_yyz_2 0_0_0_yyy_1 0_0_0_yyy_2 0_0_0_xzz_1 0_0_0_xzz_2 0_0_0_xyz_1 0_0_0_xyz_2 0_0_0_xyy_1 0_0_0_xyy_2 0_0_0_xxz_1 0_0_0_xxz_2 0_0_0_xxy_1 0_0_0_xxy_2 0_0_0_xxx_1 0_0_0_xxx_2 0_z_0_zz_2 0_z_0_zzz_1 0_z_0_zzz_2 0_z_0_yz_2 0_z_0_yzz_1 0_z_0_yzz_2 0_z_0_yy_2 0_z_0_yyz_1 0_z_0_yyz_2 0_z_0_xz_2 0_z_0_xzz_1 0_z_0_xzz_2 0_z_0_xy_2 0_z_0_xyz_1 0_z_0_xyz_2 0_z_0_xx_2 0_z_0_xxz_1 0_z_0_xxz_2 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yzz_1 0_y_0_yzz_2 0_y_0_yy_2 0_y_0_yyz_1 0_y_0_yyz_2 0_y_0_yyy_1 0_y_0_yyy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xyz_1 0_y_0_xyz_2 0_y_0_xyy_1 0_y_0_xyy_2 0_y_0_xx_2 0_y_0_xxy_1 0_y_0_xxy_2 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xzz_1 0_x_0_xzz_2 0_x_0_xy_2 0_x_0_xyz_1 0_x_0_xyz_2 0_x_0_xyy_1 0_x_0_xyy_2 0_x_0_xx_2 0_x_0_xxz_1 0_x_0_xxz_2 0_x_0_xxy_1 0_x_0_xxy_2 0_x_0_xxx_1 0_x_0_xxx_2 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SSSF_1 SSSF_2 SPSD_2 SPSF_1 SPSF_2 SDSF_1 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

INTEGRAL:2 : 3 : 2 : SSSF_2 Y SSSF_3 Y SPSD_3 Y SPSF_2 N SPSF_3 Y SDSF_2 Y 

0_0_0_zzz_2 0_0_0_zzz_3 0_0_0_yzz_2 0_0_0_yzz_3 0_0_0_yyz_2 0_0_0_yyz_3 0_0_0_yyy_2 0_0_0_yyy_3 0_0_0_xzz_2 0_0_0_xzz_3 0_0_0_xyz_2 0_0_0_xyz_3 0_0_0_xyy_2 0_0_0_xyy_3 0_0_0_xxz_2 0_0_0_xxz_3 0_0_0_xxy_2 0_0_0_xxy_3 0_0_0_xxx_2 0_0_0_xxx_3 0_z_0_zz_3 0_z_0_zzz_2 0_z_0_zzz_3 0_z_0_yz_3 0_z_0_yzz_2 0_z_0_yzz_3 0_z_0_yy_3 0_z_0_yyz_2 0_z_0_yyz_3 0_z_0_xz_3 0_z_0_xzz_2 0_z_0_xzz_3 0_z_0_xy_3 0_z_0_xyz_2 0_z_0_xyz_3 0_z_0_xx_3 0_z_0_xxz_2 0_z_0_xxz_3 0_zz_0_zzz_2 0_zz_0_yzz_2 0_zz_0_yyz_2 0_zz_0_xzz_2 0_zz_0_xyz_2 0_zz_0_xxz_2 0_y_0_yz_3 0_y_0_yy_3 0_y_0_yyz_2 0_y_0_yyz_3 0_y_0_yyy_2 0_y_0_yyy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xyz_2 0_y_0_xyz_3 0_y_0_xyy_2 0_y_0_xyy_3 0_y_0_xx_3 0_y_0_xxy_2 0_y_0_xxy_3 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_xyz_2 0_yy_0_xyy_2 0_yy_0_xxy_2 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xyz_2 0_x_0_xyz_3 0_x_0_xx_3 0_x_0_xxz_2 0_x_0_xxz_3 0_x_0_xxy_2 0_x_0_xxy_3 0_x_0_xxx_2 0_x_0_xxx_3 0_xx_0_xyz_2 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxx_2 

SSSF_2 SSSF_3 SPSD_3 SPSF_2 SPSF_3 SDSF_2 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SPSD_3 : 0_z_0_zz_3 0_z_0_yz_3 0_z_0_yy_3 0_z_0_xz_3 0_z_0_xy_3 0_z_0_xx_3 0_y_0_yz_3 0_y_0_yy_3 0_y_0_xz_3 0_y_0_xy_3 0_y_0_xx_3 0_x_0_yz_3 0_x_0_xz_3 0_x_0_xy_3 0_x_0_xx_3 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SPSF_3 : 0_z_0_zzz_3 0_z_0_yzz_3 0_z_0_yyz_3 0_z_0_xzz_3 0_z_0_xyz_3 0_z_0_xxz_3 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xxy_3 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxx_3 

SDSF_2 : 0_zz_0_zzz_2 0_zz_0_yzz_2 0_zz_0_yyz_2 0_zz_0_xzz_2 0_zz_0_xyz_2 0_zz_0_xxz_2 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_xyz_2 0_yy_0_xyy_2 0_yy_0_xxy_2 0_xx_0_xyz_2 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxx_2 

INTEGRAL:2 : 4 : 0 : SSSG_0 Y SSSG_1 Y SPSF_1 Y SPSG_0 Y SPSG_1 Y SDSG_0 Y 

0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

SSSG_0 SSSG_1 SPSF_1 SPSG_0 SPSG_1 SDSG_0 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SDSG_0 : 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

INTEGRAL:2 : 4 : 1 : SSSG_1 Y SSSG_2 Y SPSF_2 Y SPSG_1 N SPSG_2 Y SDSG_1 Y 

0_0_0_zzzz_1 0_0_0_zzzz_2 0_0_0_yzzz_1 0_0_0_yzzz_2 0_0_0_yyzz_1 0_0_0_yyzz_2 0_0_0_yyyz_1 0_0_0_yyyz_2 0_0_0_yyyy_1 0_0_0_yyyy_2 0_0_0_xzzz_1 0_0_0_xzzz_2 0_0_0_xyzz_1 0_0_0_xyzz_2 0_0_0_xyyz_1 0_0_0_xyyz_2 0_0_0_xyyy_1 0_0_0_xyyy_2 0_0_0_xxzz_1 0_0_0_xxzz_2 0_0_0_xxyz_1 0_0_0_xxyz_2 0_0_0_xxyy_1 0_0_0_xxyy_2 0_0_0_xxxz_1 0_0_0_xxxz_2 0_0_0_xxxy_1 0_0_0_xxxy_2 0_0_0_xxxx_1 0_0_0_xxxx_2 0_z_0_zzz_2 0_z_0_zzzz_1 0_z_0_zzzz_2 0_z_0_yzz_2 0_z_0_yzzz_1 0_z_0_yzzz_2 0_z_0_yyz_2 0_z_0_yyzz_1 0_z_0_yyzz_2 0_z_0_xzz_2 0_z_0_xzzz_1 0_z_0_xzzz_2 0_z_0_xyz_2 0_z_0_xyzz_1 0_z_0_xyzz_2 0_z_0_xyyz_1 0_z_0_xyyz_2 0_z_0_xxz_2 0_z_0_xxzz_1 0_z_0_xxzz_2 0_z_0_xxyz_1 0_z_0_xxyz_2 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyzz_1 0_y_0_yyzz_2 0_y_0_yyy_2 0_y_0_yyyz_1 0_y_0_yyyz_2 0_y_0_yyyy_1 0_y_0_yyyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xyyz_1 0_y_0_xyyz_2 0_y_0_xyyy_1 0_y_0_xyyy_2 0_y_0_xxy_2 0_y_0_xxyy_1 0_y_0_xxyy_2 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxzz_1 0_x_0_xxzz_2 0_x_0_xxy_2 0_x_0_xxyz_1 0_x_0_xxyz_2 0_x_0_xxyy_1 0_x_0_xxyy_2 0_x_0_xxx_2 0_x_0_xxxz_1 0_x_0_xxxz_2 0_x_0_xxxy_1 0_x_0_xxxy_2 0_x_0_xxxx_1 0_x_0_xxxx_2 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

SSSG_1 SSSG_2 SPSF_2 SPSG_1 SPSG_2 SDSG_1 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SDSG_1 : 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

INTEGRAL:2 : 4 : 2 : SSSG_2 Y SSSG_3 Y SPSF_3 Y SPSG_2 N SPSG_3 Y SDSG_2 Y 

0_0_0_zzzz_2 0_0_0_zzzz_3 0_0_0_yzzz_2 0_0_0_yzzz_3 0_0_0_yyzz_2 0_0_0_yyzz_3 0_0_0_yyyz_2 0_0_0_yyyz_3 0_0_0_yyyy_2 0_0_0_yyyy_3 0_0_0_xzzz_2 0_0_0_xzzz_3 0_0_0_xyzz_2 0_0_0_xyzz_3 0_0_0_xyyz_2 0_0_0_xyyz_3 0_0_0_xyyy_2 0_0_0_xyyy_3 0_0_0_xxzz_2 0_0_0_xxzz_3 0_0_0_xxyz_2 0_0_0_xxyz_3 0_0_0_xxyy_2 0_0_0_xxyy_3 0_0_0_xxxz_2 0_0_0_xxxz_3 0_0_0_xxxy_2 0_0_0_xxxy_3 0_0_0_xxxx_2 0_0_0_xxxx_3 0_z_0_zzz_3 0_z_0_zzzz_2 0_z_0_zzzz_3 0_z_0_yzz_3 0_z_0_yzzz_2 0_z_0_yzzz_3 0_z_0_yyz_3 0_z_0_yyzz_2 0_z_0_yyzz_3 0_z_0_xzz_3 0_z_0_xzzz_2 0_z_0_xzzz_3 0_z_0_xyz_3 0_z_0_xyzz_2 0_z_0_xyzz_3 0_z_0_xxz_3 0_z_0_xxzz_2 0_z_0_xxzz_3 0_zz_0_zzzz_2 0_zz_0_yzzz_2 0_zz_0_yyzz_2 0_zz_0_xzzz_2 0_zz_0_xyzz_2 0_zz_0_xxzz_2 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_yyyz_2 0_y_0_yyyz_3 0_y_0_yyyy_2 0_y_0_yyyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xyyz_2 0_y_0_xyyz_3 0_y_0_xyyy_2 0_y_0_xyyy_3 0_y_0_xxy_3 0_y_0_xxyy_2 0_y_0_xxyy_3 0_yy_0_yyyz_2 0_yy_0_yyyy_2 0_yy_0_xyyz_2 0_yy_0_xyyy_2 0_yy_0_xxyy_2 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxyz_2 0_x_0_xxyz_3 0_x_0_xxx_3 0_x_0_xxxz_2 0_x_0_xxxz_3 0_x_0_xxxy_2 0_x_0_xxxy_3 0_x_0_xxxx_2 0_x_0_xxxx_3 0_xx_0_xxyz_2 0_xx_0_xxxz_2 0_xx_0_xxxy_2 0_xx_0_xxxx_2 

SSSG_2 SSSG_3 SPSF_3 SPSG_2 SPSG_3 SDSG_2 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SPSF_3 : 0_z_0_zzz_3 0_z_0_yzz_3 0_z_0_yyz_3 0_z_0_xzz_3 0_z_0_xyz_3 0_z_0_xxz_3 0_y_0_yyz_3 0_y_0_yyy_3 0_y_0_xyz_3 0_y_0_xyy_3 0_y_0_xxy_3 0_x_0_xyz_3 0_x_0_xxz_3 0_x_0_xxy_3 0_x_0_xxx_3 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SPSG_3 : 0_z_0_zzzz_3 0_z_0_yzzz_3 0_z_0_yyzz_3 0_z_0_xzzz_3 0_z_0_xyzz_3 0_z_0_xxzz_3 0_y_0_yyyz_3 0_y_0_yyyy_3 0_y_0_xyyz_3 0_y_0_xyyy_3 0_y_0_xxyy_3 0_x_0_xxyz_3 0_x_0_xxxz_3 0_x_0_xxxy_3 0_x_0_xxxx_3 

SDSG_2 : 0_zz_0_zzzz_2 0_zz_0_yzzz_2 0_zz_0_yyzz_2 0_zz_0_xzzz_2 0_zz_0_xyzz_2 0_zz_0_xxzz_2 0_yy_0_yyyz_2 0_yy_0_yyyy_2 0_yy_0_xyyz_2 0_yy_0_xyyy_2 0_yy_0_xxyy_2 0_xx_0_xxyz_2 0_xx_0_xxxz_2 0_xx_0_xxxy_2 0_xx_0_xxxx_2 

INTEGRAL:3 : 1 : 1 : SPSP_1 N SPSP_2 N SDSS_2 Y SDSP_1 N SDSP_2 N SFSP_1 Y 

0_z_0_z_1 0_z_0_z_2 0_zz_0_0_2 0_zz_0_z_1 0_zz_0_z_2 0_zz_0_y_1 0_zz_0_y_2 0_zzz_0_z_1 0_y_0_y_1 0_y_0_y_2 0_yzz_0_z_1 0_yzz_0_y_1 0_yy_0_0_2 0_yy_0_z_1 0_yy_0_z_2 0_yy_0_y_1 0_yy_0_y_2 0_yyz_0_z_1 0_yyy_0_y_1 0_x_0_x_1 0_x_0_x_2 0_xzz_0_z_1 0_xyy_0_y_1 0_xx_0_0_2 0_xx_0_z_1 0_xx_0_z_2 0_xx_0_x_1 0_xx_0_x_2 0_xxz_0_z_1 0_xxx_0_x_1 

SPSP_1 SPSP_2 SDSS_2 SDSP_1 SDSP_2 SFSP_1 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SDSS_2 : 0_zz_0_0_2 0_yy_0_0_2 0_xx_0_0_2 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SDSP_2 : 0_zz_0_z_2 0_zz_0_y_2 0_zz_0_x_2 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_x_2 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_x_2 

SFSP_1 : 0_zzz_0_z_1 0_yzz_0_z_1 0_yzz_0_y_1 0_yyz_0_z_1 0_yyy_0_y_1 0_xzz_0_z_1 0_xyy_0_y_1 0_xxz_0_z_1 0_xxx_0_x_1 

INTEGRAL:3 : 2 : 0 : SPSD_0 N SPSD_1 N SDSP_1 Y SDSD_0 Y SDSD_1 Y SFSD_0 Y 

0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_xz_0 0_z_0_xz_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xy_0 0_y_0_xy_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xz_1 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

SPSD_0 SPSD_1 SDSP_1 SDSD_0 SDSD_1 SFSD_0 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_0_xy_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_x_0_yz_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SFSD_0 : 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

INTEGRAL:3 : 2 : 1 : SPSD_1 N SPSD_2 N SDSP_2 Y SDSD_1 N SDSD_2 Y SFSD_1 Y 

0_z_0_zz_1 0_z_0_zz_2 0_z_0_yz_1 0_z_0_yz_2 0_z_0_xz_1 0_z_0_xz_2 0_zz_0_z_2 0_zz_0_zz_1 0_zz_0_zz_2 0_zz_0_y_2 0_zz_0_yz_1 0_zz_0_yz_2 0_zz_0_x_2 0_zz_0_xz_1 0_zz_0_xz_2 0_zz_0_xy_1 0_zz_0_xy_2 0_zzz_0_zz_1 0_zzz_0_yz_1 0_zzz_0_xz_1 0_y_0_yz_1 0_y_0_yz_2 0_y_0_yy_1 0_y_0_yy_2 0_y_0_xy_1 0_y_0_xy_2 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_xy_1 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_yz_1 0_yy_0_yz_2 0_yy_0_yy_1 0_yy_0_yy_2 0_yy_0_x_2 0_yy_0_xz_1 0_yy_0_xz_2 0_yy_0_xy_1 0_yy_0_xy_2 0_yyz_0_yz_1 0_yyz_0_yy_1 0_yyz_0_xz_1 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_xy_1 0_x_0_xz_1 0_x_0_xz_2 0_x_0_xy_1 0_x_0_xy_2 0_x_0_xx_1 0_x_0_xx_2 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_yz_1 0_xx_0_yz_2 0_xx_0_x_2 0_xx_0_xz_1 0_xx_0_xz_2 0_xx_0_xy_1 0_xx_0_xy_2 0_xx_0_xx_1 0_xx_0_xx_2 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xx_1 0_xxy_0_xy_1 0_xxy_0_xx_1 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 

SPSD_1 SPSD_2 SDSP_2 SDSD_1 SDSD_2 SFSD_1 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SDSP_2 : 0_zz_0_z_2 0_zz_0_y_2 0_zz_0_x_2 0_yy_0_z_2 0_yy_0_y_2 0_yy_0_x_2 0_xx_0_z_2 0_xx_0_y_2 0_xx_0_x_2 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SDSD_2 : 0_zz_0_zz_2 0_zz_0_yz_2 0_zz_0_xz_2 0_zz_0_xy_2 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_xz_2 0_yy_0_xy_2 0_xx_0_yz_2 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xx_2 

SFSD_1 : 0_zzz_0_zz_1 0_zzz_0_yz_1 0_zzz_0_xz_1 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_xy_1 0_yyz_0_yz_1 0_yyz_0_yy_1 0_yyz_0_xz_1 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_xy_1 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xx_1 0_xxy_0_xy_1 0_xxy_0_xx_1 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 

INTEGRAL:3 : 3 : 0 : SPSF_0 N SPSF_1 N SDSD_1 N SDSF_0 N SDSF_1 Y SFSF_0 Y 

0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xyy_0 0_y_0_xyy_1 0_yz_0_yz_1 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

SPSF_0 SPSF_1 SDSD_1 SDSF_0 SDSF_1 SFSF_0 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xxz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SFSF_0 : 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

INTEGRAL:3 : 3 : 1 : SPSF_1 N SPSF_2 N SDSD_2 N SDSF_1 N SDSF_2 Y SFSF_1 Y 

0_z_0_zzz_1 0_z_0_zzz_2 0_z_0_yzz_1 0_z_0_yzz_2 0_z_0_xzz_1 0_z_0_xzz_2 0_zz_0_zz_2 0_zz_0_zzz_1 0_zz_0_zzz_2 0_zz_0_yz_2 0_zz_0_yzz_1 0_zz_0_yzz_2 0_zz_0_yyz_1 0_zz_0_yyz_2 0_zz_0_xz_2 0_zz_0_xzz_1 0_zz_0_xzz_2 0_zz_0_xyz_1 0_zz_0_xyz_2 0_zz_0_xxz_1 0_zz_0_xxz_2 0_zzz_0_zzz_1 0_zzz_0_yzz_1 0_zzz_0_xzz_1 0_y_0_yyz_1 0_y_0_yyz_2 0_y_0_yyy_1 0_y_0_yyy_2 0_y_0_xyy_1 0_y_0_xyy_2 0_yzz_0_yzz_1 0_yzz_0_yyz_1 0_yzz_0_xyz_1 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_yyz_1 0_yy_0_yyz_2 0_yy_0_yyy_1 0_yy_0_yyy_2 0_yy_0_xy_2 0_yy_0_xyz_1 0_yy_0_xyz_2 0_yy_0_xyy_1 0_yy_0_xyy_2 0_yy_0_xxy_1 0_yy_0_xxy_2 0_yyz_0_yyz_1 0_yyz_0_xyz_1 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_xyy_1 0_x_0_xxz_1 0_x_0_xxz_2 0_x_0_xxy_1 0_x_0_xxy_2 0_x_0_xxx_1 0_x_0_xxx_2 0_xzz_0_xzz_1 0_xzz_0_xxz_1 0_xyy_0_xyy_1 0_xyy_0_xxy_1 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xyz_1 0_xx_0_xyz_2 0_xx_0_xx_2 0_xx_0_xxz_1 0_xx_0_xxz_2 0_xx_0_xxy_1 0_xx_0_xxy_2 0_xx_0_xxx_1 0_xx_0_xxx_2 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxy_0_xxy_1 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 

SPSF_1 SPSF_2 SDSD_2 SDSF_1 SDSF_2 SFSF_1 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xxz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xxz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SDSD_2 : 0_zz_0_zz_2 0_zz_0_yz_2 0_zz_0_xz_2 0_zz_0_xy_2 0_yy_0_yz_2 0_yy_0_yy_2 0_yy_0_xz_2 0_yy_0_xy_2 0_xx_0_yz_2 0_xx_0_xz_2 0_xx_0_xy_2 0_xx_0_xx_2 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SDSF_2 : 0_zz_0_zzz_2 0_zz_0_yzz_2 0_zz_0_yyz_2 0_zz_0_xzz_2 0_zz_0_xyz_2 0_zz_0_xxz_2 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_xyz_2 0_yy_0_xyy_2 0_yy_0_xxy_2 0_xx_0_xyz_2 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxx_2 

SFSF_1 : 0_zzz_0_zzz_1 0_zzz_0_yzz_1 0_zzz_0_xzz_1 0_yzz_0_yzz_1 0_yzz_0_yyz_1 0_yzz_0_xyz_1 0_yyz_0_yyz_1 0_yyz_0_xyz_1 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_xyy_1 0_xzz_0_xzz_1 0_xzz_0_xxz_1 0_xyy_0_xyy_1 0_xyy_0_xxy_1 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxy_0_xxy_1 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 

INTEGRAL:3 : 4 : 0 : SPSG_0 N SPSG_1 N SDSF_1 N SDSG_0 N SDSG_1 Y SFSG_0 Y 

0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_zz_0_zzz_1 0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yzz_1 0_zz_0_yzzz_0 0_zz_0_yzzz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_xzz_1 0_zz_0_xzzz_0 0_zz_0_xzzz_1 0_zz_0_xyzz_0 0_zz_0_xyzz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yz_0_xyzz_0 0_yz_0_xyzz_1 0_yz_0_xyyz_0 0_yz_0_xyyz_1 0_yz_0_xxyz_0 0_yz_0_xxyz_1 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yy_0_yyz_1 0_yy_0_yyzz_0 0_yy_0_yyzz_1 0_yy_0_yyy_1 0_yy_0_yyyz_0 0_yy_0_yyyz_1 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xyy_1 0_yy_0_xyyz_0 0_yy_0_xyyz_1 0_yy_0_xyyy_0 0_yy_0_xyyy_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxyy_0 0_xx_0_xxz_1 0_xx_0_xxzz_0 0_xx_0_xxzz_1 0_xx_0_xxy_1 0_xx_0_xxyz_0 0_xx_0_xxyz_1 0_xx_0_xxyy_0 0_xx_0_xxyy_1 0_xx_0_xxx_1 0_xx_0_xxxz_0 0_xx_0_xxxz_1 0_xx_0_xxxy_0 0_xx_0_xxxy_1 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxxz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

SPSG_0 SPSG_1 SDSF_1 SDSG_0 SDSG_1 SFSG_0 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxxz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SDSG_0 : 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

SDSG_1 : 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

SFSG_0 : 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxyy_0 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxxz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

INTEGRAL:3 : 4 : 1 : SPSG_1 N SPSG_2 N SDSF_2 N SDSG_1 N SDSG_2 Y SFSG_1 Y 

0_z_0_zzzz_1 0_z_0_zzzz_2 0_z_0_yzzz_1 0_z_0_yzzz_2 0_z_0_xzzz_1 0_z_0_xzzz_2 0_zz_0_zzz_2 0_zz_0_zzzz_1 0_zz_0_zzzz_2 0_zz_0_yzz_2 0_zz_0_yzzz_1 0_zz_0_yzzz_2 0_zz_0_yyzz_1 0_zz_0_yyzz_2 0_zz_0_xzz_2 0_zz_0_xzzz_1 0_zz_0_xzzz_2 0_zz_0_xyzz_1 0_zz_0_xyzz_2 0_zz_0_xxzz_1 0_zz_0_xxzz_2 0_zzz_0_zzzz_1 0_zzz_0_yzzz_1 0_zzz_0_xzzz_1 0_y_0_yyyz_1 0_y_0_yyyz_2 0_y_0_yyyy_1 0_y_0_yyyy_2 0_y_0_xyyy_1 0_y_0_xyyy_2 0_yzz_0_yyzz_1 0_yzz_0_xyzz_1 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_yyyz_1 0_yy_0_yyyz_2 0_yy_0_yyyy_1 0_yy_0_yyyy_2 0_yy_0_xyy_2 0_yy_0_xyyz_1 0_yy_0_xyyz_2 0_yy_0_xyyy_1 0_yy_0_xyyy_2 0_yy_0_xxyy_1 0_yy_0_xxyy_2 0_yyz_0_xyyz_1 0_yyy_0_yyyz_1 0_yyy_0_yyyy_1 0_yyy_0_xyyy_1 0_x_0_xxxz_1 0_x_0_xxxz_2 0_x_0_xxxy_1 0_x_0_xxxy_2 0_x_0_xxxx_1 0_x_0_xxxx_2 0_xzz_0_xxzz_1 0_xyy_0_xxyy_1 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxyz_1 0_xx_0_xxyz_2 0_xx_0_xxx_2 0_xx_0_xxxz_1 0_xx_0_xxxz_2 0_xx_0_xxxy_1 0_xx_0_xxxy_2 0_xx_0_xxxx_1 0_xx_0_xxxx_2 0_xxz_0_xxyz_1 0_xxx_0_xxxz_1 0_xxx_0_xxxy_1 0_xxx_0_xxxx_1 

SPSG_1 SPSG_2 SDSF_2 SDSG_1 SDSG_2 SFSG_1 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxxz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SDSF_2 : 0_zz_0_zzz_2 0_zz_0_yzz_2 0_zz_0_yyz_2 0_zz_0_xzz_2 0_zz_0_xyz_2 0_zz_0_xxz_2 0_yy_0_yyz_2 0_yy_0_yyy_2 0_yy_0_xyz_2 0_yy_0_xyy_2 0_yy_0_xxy_2 0_xx_0_xyz_2 0_xx_0_xxz_2 0_xx_0_xxy_2 0_xx_0_xxx_2 

SDSG_1 : 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

SDSG_2 : 0_zz_0_zzzz_2 0_zz_0_yzzz_2 0_zz_0_yyzz_2 0_zz_0_xzzz_2 0_zz_0_xyzz_2 0_zz_0_xxzz_2 0_yy_0_yyyz_2 0_yy_0_yyyy_2 0_yy_0_xyyz_2 0_yy_0_xyyy_2 0_yy_0_xxyy_2 0_xx_0_xxyz_2 0_xx_0_xxxz_2 0_xx_0_xxxy_2 0_xx_0_xxxx_2 

SFSG_1 : 0_zzz_0_zzzz_1 0_zzz_0_yzzz_1 0_zzz_0_xzzz_1 0_yzz_0_yyzz_1 0_yzz_0_xyzz_1 0_yyz_0_xyyz_1 0_yyy_0_yyyz_1 0_yyy_0_yyyy_1 0_yyy_0_xyyy_1 0_xzz_0_xxzz_1 0_xyy_0_xxyy_1 0_xxz_0_xxyz_1 0_xxx_0_xxxz_1 0_xxx_0_xxxy_1 0_xxx_0_xxxx_1 

INTEGRAL:4 : 2 : 0 : SDSD_0 N SDSD_1 N SFSP_1 Y SFSD_0 N SFSD_1 Y SGSD_0 Y 

0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zzz_0_z_1 0_zzz_0_zz_0 0_zzz_0_zz_1 0_zzz_0_yz_0 0_zzz_0_yz_1 0_zzz_0_xz_0 0_zzz_0_xz_1 0_zzzz_0_zz_0 0_yzz_0_z_1 0_yzz_0_zz_0 0_yzz_0_zz_1 0_yzz_0_y_1 0_yzz_0_yz_0 0_yzz_0_yz_1 0_yzz_0_xy_0 0_yzz_0_xy_1 0_yzzz_0_zz_0 0_yzzz_0_yz_0 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yyz_0_z_1 0_yyz_0_yz_0 0_yyz_0_yz_1 0_yyz_0_yy_0 0_yyz_0_yy_1 0_yyz_0_xz_0 0_yyz_0_xz_1 0_yyzz_0_zz_0 0_yyzz_0_yz_0 0_yyzz_0_yy_0 0_yyy_0_y_1 0_yyy_0_yz_0 0_yyy_0_yz_1 0_yyy_0_yy_0 0_yyy_0_yy_1 0_yyy_0_xy_0 0_yyy_0_xy_1 0_yyyz_0_yz_0 0_yyyz_0_yy_0 0_yyyy_0_yy_0 0_xzz_0_z_1 0_xzz_0_zz_0 0_xzz_0_zz_1 0_xzz_0_xz_0 0_xzz_0_xz_1 0_xzzz_0_zz_0 0_xzzz_0_xz_0 0_xyzz_0_zz_0 0_xyzz_0_yz_0 0_xyzz_0_xz_0 0_xyzz_0_xy_0 0_xyy_0_y_1 0_xyy_0_yy_0 0_xyy_0_yy_1 0_xyy_0_xy_0 0_xyy_0_xy_1 0_xyyz_0_yz_0 0_xyyz_0_yy_0 0_xyyz_0_xz_0 0_xyyz_0_xy_0 0_xyyy_0_yy_0 0_xyyy_0_xy_0 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_z_1 0_xxz_0_yz_0 0_xxz_0_yz_1 0_xxz_0_xz_0 0_xxz_0_xz_1 0_xxz_0_xx_0 0_xxz_0_xx_1 0_xxzz_0_zz_0 0_xxzz_0_xz_0 0_xxzz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xy_1 0_xxy_0_xx_0 0_xxy_0_xx_1 0_xxyz_0_yz_0 0_xxyz_0_xz_0 0_xxyz_0_xy_0 0_xxyz_0_xx_0 0_xxyy_0_yy_0 0_xxyy_0_xy_0 0_xxyy_0_xx_0 0_xxx_0_x_1 0_xxx_0_xz_0 0_xxx_0_xz_1 0_xxx_0_xy_0 0_xxx_0_xy_1 0_xxx_0_xx_0 0_xxx_0_xx_1 0_xxxz_0_xz_0 0_xxxz_0_xx_0 0_xxxy_0_xy_0 0_xxxy_0_xx_0 0_xxxx_0_xx_0 

SDSD_0 SDSD_1 SFSP_1 SFSD_0 SFSD_1 SGSD_0 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_xz_0 0_zz_0_xy_0 0_yz_0_yz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_yz_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_xz_1 0_zz_0_xy_1 0_yz_0_yz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_xz_0_xz_1 0_xy_0_xy_1 0_xx_0_yz_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SFSP_1 : 0_zzz_0_z_1 0_yzz_0_z_1 0_yzz_0_y_1 0_yyz_0_z_1 0_yyy_0_y_1 0_xzz_0_z_1 0_xyy_0_y_1 0_xxz_0_z_1 0_xxx_0_x_1 

SFSD_0 : 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_xz_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_xy_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xy_0 0_xzz_0_zz_0 0_xzz_0_xz_0 0_xyz_0_yz_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyy_0_yy_0 0_xyy_0_xy_0 0_xxz_0_yz_0 0_xxz_0_xz_0 0_xxz_0_xx_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

SFSD_1 : 0_zzz_0_zz_1 0_zzz_0_yz_1 0_zzz_0_xz_1 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_xy_1 0_yyz_0_yz_1 0_yyz_0_yy_1 0_yyz_0_xz_1 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_xy_1 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xx_1 0_xxy_0_xy_1 0_xxy_0_xx_1 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 

SGSD_0 : 0_zzzz_0_zz_0 0_yzzz_0_zz_0 0_yzzz_0_yz_0 0_yyzz_0_zz_0 0_yyzz_0_yz_0 0_yyzz_0_yy_0 0_yyyz_0_yz_0 0_yyyz_0_yy_0 0_yyyy_0_yy_0 0_xzzz_0_zz_0 0_xzzz_0_xz_0 0_xyzz_0_zz_0 0_xyzz_0_yz_0 0_xyzz_0_xz_0 0_xyzz_0_xy_0 0_xyyz_0_yz_0 0_xyyz_0_yy_0 0_xyyz_0_xz_0 0_xyyz_0_xy_0 0_xyyy_0_yy_0 0_xyyy_0_xy_0 0_xxzz_0_zz_0 0_xxzz_0_xz_0 0_xxzz_0_xx_0 0_xxyz_0_yz_0 0_xxyz_0_xz_0 0_xxyz_0_xy_0 0_xxyz_0_xx_0 0_xxyy_0_yy_0 0_xxyy_0_xy_0 0_xxyy_0_xx_0 0_xxxz_0_xz_0 0_xxxz_0_xx_0 0_xxxy_0_xy_0 0_xxxy_0_xx_0 0_xxxx_0_xx_0 

INTEGRAL:4 : 3 : 0 : SDSF_0 N SDSF_1 N SFSD_1 N SFSF_0 N SFSF_1 Y SGSF_0 Y 

0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zzz_0_zz_1 0_zzz_0_zzz_0 0_zzz_0_zzz_1 0_zzz_0_yzz_0 0_zzz_0_yzz_1 0_zzz_0_xzz_0 0_zzz_0_xzz_1 0_zzzz_0_zzz_0 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_yzz_0 0_yzz_0_yzz_1 0_yzz_0_yyz_0 0_yzz_0_yyz_1 0_yzz_0_xyz_0 0_yzz_0_xyz_1 0_yzzz_0_zzz_0 0_yzzz_0_yzz_0 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yyz_0_yz_1 0_yyz_0_yyz_0 0_yyz_0_yyz_1 0_yyz_0_xyz_0 0_yyz_0_xyz_1 0_yyzz_0_yzz_0 0_yyzz_0_yyz_0 0_yyy_0_yy_1 0_yyy_0_yyz_0 0_yyy_0_yyz_1 0_yyy_0_yyy_0 0_yyy_0_yyy_1 0_yyy_0_xyy_0 0_yyy_0_xyy_1 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyy_0_yyy_0 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xzz_0_xzz_0 0_xzz_0_xzz_1 0_xzz_0_xxz_0 0_xzz_0_xxz_1 0_xzzz_0_zzz_0 0_xzzz_0_xzz_0 0_xyzz_0_yzz_0 0_xyzz_0_xzz_0 0_xyzz_0_xyz_0 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xyy_0_xyy_0 0_xyy_0_xyy_1 0_xyy_0_xxy_0 0_xyy_0_xxy_1 0_xyyz_0_yyz_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyy_0_yyy_0 0_xyyy_0_xyy_0 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_xz_1 0_xxz_0_xyz_0 0_xxz_0_xyz_1 0_xxz_0_xxz_0 0_xxz_0_xxz_1 0_xxzz_0_xzz_0 0_xxzz_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxy_1 0_xxyz_0_xyz_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyy_0_xyy_0 0_xxyy_0_xxy_0 0_xxx_0_xx_1 0_xxx_0_xxz_0 0_xxx_0_xxz_1 0_xxx_0_xxy_0 0_xxx_0_xxy_1 0_xxx_0_xxx_0 0_xxx_0_xxx_1 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxx_0_xxx_0 

SDSF_0 SDSF_1 SFSD_1 SFSF_0 SFSF_1 SGSF_0 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xxz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_xyz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xxz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xxz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_xyz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxy_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SFSD_1 : 0_zzz_0_zz_1 0_zzz_0_yz_1 0_zzz_0_xz_1 0_yzz_0_zz_1 0_yzz_0_yz_1 0_yzz_0_xy_1 0_yyz_0_yz_1 0_yyz_0_yy_1 0_yyz_0_xz_1 0_yyy_0_yz_1 0_yyy_0_yy_1 0_yyy_0_xy_1 0_xzz_0_zz_1 0_xzz_0_xz_1 0_xyy_0_yy_1 0_xyy_0_xy_1 0_xxz_0_yz_1 0_xxz_0_xz_1 0_xxz_0_xx_1 0_xxy_0_xy_1 0_xxy_0_xx_1 0_xxx_0_xz_1 0_xxx_0_xy_1 0_xxx_0_xx_1 

SFSF_0 : 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_xzz_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xyy_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xxz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

SFSF_1 : 0_zzz_0_zzz_1 0_zzz_0_yzz_1 0_zzz_0_xzz_1 0_yzz_0_yzz_1 0_yzz_0_yyz_1 0_yzz_0_xyz_1 0_yyz_0_yyz_1 0_yyz_0_xyz_1 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_xyy_1 0_xzz_0_xzz_1 0_xzz_0_xxz_1 0_xyy_0_xyy_1 0_xyy_0_xxy_1 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxy_0_xxy_1 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 

SGSF_0 : 0_zzzz_0_zzz_0 0_yzzz_0_zzz_0 0_yzzz_0_yzz_0 0_yyzz_0_yzz_0 0_yyzz_0_yyz_0 0_yyyz_0_yyz_0 0_yyyz_0_yyy_0 0_yyyy_0_yyy_0 0_xzzz_0_zzz_0 0_xzzz_0_xzz_0 0_xyzz_0_yzz_0 0_xyzz_0_xzz_0 0_xyzz_0_xyz_0 0_xyyz_0_yyz_0 0_xyyz_0_xyz_0 0_xyyz_0_xyy_0 0_xyyy_0_yyy_0 0_xyyy_0_xyy_0 0_xxzz_0_xzz_0 0_xxzz_0_xxz_0 0_xxyz_0_xyz_0 0_xxyz_0_xxz_0 0_xxyz_0_xxy_0 0_xxyy_0_xyy_0 0_xxyy_0_xxy_0 0_xxxz_0_xxz_0 0_xxxz_0_xxx_0 0_xxxy_0_xxy_0 0_xxxy_0_xxx_0 0_xxxx_0_xxx_0 

INTEGRAL:4 : 4 : 0 : SDSG_0 N SDSG_1 N SFSF_1 N SFSG_0 N SFSG_1 Y SGSG_0 Y 

0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zzz_0_zzz_1 0_zzz_0_zzzz_0 0_zzz_0_zzzz_1 0_zzz_0_yzzz_0 0_zzz_0_yzzz_1 0_zzz_0_xzzz_0 0_zzz_0_xzzz_1 0_zzzz_0_zzzz_0 0_yzz_0_yzz_1 0_yzz_0_yyzz_0 0_yzz_0_yyzz_1 0_yzz_0_xyzz_0 0_yzz_0_xyzz_1 0_yzzz_0_yzzz_0 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yyz_0_yyz_1 0_yyz_0_xyyz_0 0_yyz_0_xyyz_1 0_yyzz_0_yyzz_0 0_yyy_0_yyy_1 0_yyy_0_yyyz_0 0_yyy_0_yyyz_1 0_yyy_0_yyyy_0 0_yyy_0_yyyy_1 0_yyy_0_xyyy_0 0_yyy_0_xyyy_1 0_yyyz_0_yyyz_0 0_yyyy_0_yyyy_0 0_xzz_0_xzz_1 0_xzz_0_xxzz_0 0_xzz_0_xxzz_1 0_xzzz_0_xzzz_0 0_xyzz_0_xyzz_0 0_xyy_0_xyy_1 0_xyy_0_xxyy_0 0_xyy_0_xxyy_1 0_xyyz_0_xyyz_0 0_xyyy_0_xyyy_0 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_xxz_1 0_xxz_0_xxyz_0 0_xxz_0_xxyz_1 0_xxzz_0_xxzz_0 0_xxyz_0_xxyz_0 0_xxyy_0_xxyy_0 0_xxx_0_xxx_1 0_xxx_0_xxxz_0 0_xxx_0_xxxz_1 0_xxx_0_xxxy_0 0_xxx_0_xxxy_1 0_xxx_0_xxxx_0 0_xxx_0_xxxx_1 0_xxxz_0_xxxz_0 0_xxxy_0_xxxy_0 0_xxxx_0_xxxx_0 

SDSG_0 SDSG_1 SFSF_1 SFSG_0 SFSG_1 SGSG_0 

SDSG_0 : 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xxzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxxz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

SDSG_1 : 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xxzz_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

SFSF_1 : 0_zzz_0_zzz_1 0_zzz_0_yzz_1 0_zzz_0_xzz_1 0_yzz_0_yzz_1 0_yzz_0_yyz_1 0_yzz_0_xyz_1 0_yyz_0_yyz_1 0_yyz_0_xyz_1 0_yyy_0_yyz_1 0_yyy_0_yyy_1 0_yyy_0_xyy_1 0_xzz_0_xzz_1 0_xzz_0_xxz_1 0_xyy_0_xyy_1 0_xyy_0_xxy_1 0_xxz_0_xyz_1 0_xxz_0_xxz_1 0_xxy_0_xxy_1 0_xxx_0_xxz_1 0_xxx_0_xxy_1 0_xxx_0_xxx_1 

SFSG_0 : 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_xzzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_xyzz_0 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_xyyz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xyyy_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xxzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xxyz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxyy_0 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxxz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

SFSG_1 : 0_zzz_0_zzzz_1 0_zzz_0_yzzz_1 0_zzz_0_xzzz_1 0_yzz_0_yyzz_1 0_yzz_0_xyzz_1 0_yyz_0_xyyz_1 0_yyy_0_yyyz_1 0_yyy_0_yyyy_1 0_yyy_0_xyyy_1 0_xzz_0_xxzz_1 0_xyy_0_xxyy_1 0_xxz_0_xxyz_1 0_xxx_0_xxxz_1 0_xxx_0_xxxy_1 0_xxx_0_xxxx_1 

SGSG_0 : 0_zzzz_0_zzzz_0 0_yzzz_0_yzzz_0 0_yyzz_0_yyzz_0 0_yyyz_0_yyyz_0 0_yyyy_0_yyyy_0 0_xzzz_0_xzzz_0 0_xyzz_0_xyzz_0 0_xyyz_0_xyyz_0 0_xyyy_0_xyyy_0 0_xxzz_0_xxzz_0 0_xxyz_0_xxyz_0 0_xxyy_0_xxyy_0 0_xxxz_0_xxxz_0 0_xxxy_0_xxxy_0 0_xxxx_0_xxxx_0 

        }
    }
}


} // derirec namespace
