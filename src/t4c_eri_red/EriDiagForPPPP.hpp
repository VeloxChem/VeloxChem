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
#include "EriDiagVRRForSPSS.hpp"
#include "EriDiagVRRForSPSP.hpp"
#include "EriDiagVRRForSPSD.hpp"
#include "EriDiagVRRForSDSP.hpp"
#include "EriDiagVRRForSDSD.hpp"
#include "EriDiagHRRForSPPP.hpp"
#include "EriDiagHRRForSDPP.hpp"
#include "EriDiagHRRForPPPP.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPPPP(      T*                                 intsBuffer,
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
    BufferHostXY<T> pbufSSSD2(6, ncpairs);

0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 
    BufferHostXY<T> pbufSPSS1(3, ncpairs);

0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 
    BufferHostXY<T> pbufSPSP0(3, ncpairs);

0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 
    BufferHostXY<T> pbufSPSP1(3, ncpairs);

0_z_0_z_1 0_y_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSPSD0(9, ncpairs);

0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 
    BufferHostXY<T> pbufSPSD1(6, ncpairs);

0_z_0_zz_1 0_z_0_yz_1 0_z_0_xz_1 0_y_0_yy_1 0_y_0_xy_1 0_x_0_xx_1 
    BufferHostXY<T> pbufSDSP0(9, ncpairs);

0_zz_0_z_0 0_yz_0_z_0 0_yz_0_y_0 0_yy_0_y_0 0_xz_0_z_0 0_xz_0_x_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_x_0 
    BufferHostXY<T> pbufSDSD0(6, ncpairs);

0_zz_0_zz_0 0_yz_0_yz_0 0_yy_0_yy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_xx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSPSP(3, ncpairs);

0_z_0_z 0_y_0_y 0_x_0_x 
    BufferHostXY<T> cbufSPSD(9, ncpairs);

0_z_0_zz 0_z_0_yz 0_z_0_xz 0_y_0_yz 0_y_0_yy 0_y_0_xy 0_x_0_xz 0_x_0_xy 0_x_0_xx 
    BufferHostXY<T> cbufSDSP(9, ncpairs);

0_zz_0_z 0_yz_0_z 0_yz_0_y 0_yy_0_y 0_xz_0_z 0_xz_0_x 0_xy_0_y 0_xy_0_x 0_xx_0_x 
    BufferHostXY<T> cbufSDSD(6, ncpairs);

0_zz_0_zz 0_yz_0_yz 0_yy_0_yy 0_xz_0_xz 0_xy_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSPPP(9, ncpairs);

0_z_z_z 0_z_y_z 0_z_x_z 0_y_z_y 0_y_y_y 0_y_x_y 0_x_z_x 0_x_y_x 0_x_x_x 
    BufferHostXY<T> cbufSDPP(9, ncpairs);

0_zz_z_z 0_yz_z_y 0_yz_y_z 0_yy_y_y 0_xz_z_x 0_xz_x_z 0_xy_y_x 0_xy_x_y 0_xx_x_x 
    BufferHostXY<T> cbufPPPP(9, ncpairs);

z_z_z_z z_y_z_y z_x_z_x y_z_y_z y_y_y_y y_x_y_x x_z_x_z x_y_x_y x_x_x_x 
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
0_z_z_z_0 0_z_y_z_0 0_z_x_z_0 0_zz_z_z_0 0_y_z_y_0 0_y_y_y_0 0_y_x_y_0 0_yz_z_y_0 0_yz_y_z_0 0_yy_y_y_0 0_x_z_x_0 0_x_y_x_0 0_x_x_x_0 0_xz_z_x_0 0_xz_x_z_0 0_xy_y_x_0 0_xy_x_y_0 0_xx_x_x_0 z_z_z_z_0 z_y_z_y_0 z_x_z_x_0 y_z_y_z_0 y_y_y_y_0 y_x_y_x_0 x_z_x_z_0 x_y_x_y_0 x_x_x_x_0 

signature:
0_zz_0_z_0 0_zz_0_zz_0 0_zz_z_z_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_yz_0 0_yz_z_y_0 0_yz_y_z_0 0_yy_0_y_0 0_yy_0_yy_0 0_yy_y_y_0 0_xz_0_z_0 0_xz_0_x_0 0_xz_0_xz_0 0_xz_z_x_0 0_xz_x_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xy_0_xy_0 0_xy_y_x_0 0_xy_x_y_0 0_xx_0_x_0 0_xx_0_xx_0 0_xx_x_x_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_xz_0 0_z_0_xz_1 0_zz_0_zz_0 0_y_0_y_1 0_y_0_yz_0 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xy_0 0_y_0_xy_1 0_yz_0_yz_0 0_yy_0_yy_0 0_x_0_x_1 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_zz_0_z_0 0_y_0_0_1 0_y_0_y_0 0_y_0_y_1 0_yz_0_z_0 0_yz_0_y_0 0_yy_0_y_0 0_x_0_0_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_x_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_x_0 

signature:
0_z_0_z_0 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_z_z_z_0 0_z_y_z_0 0_z_x_z_0 0_y_0_y_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_y_z_y_0 0_y_y_y_0 0_y_x_y_0 0_x_0_x_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 0_x_z_x_0 0_x_y_x_0 0_x_x_x_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

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

0_0_0_0_2 0_0_0_0_3 0_0_0_z_2 0_0_0_z_3 0_0_0_zz_2 0_0_0_y_2 0_0_0_y_3 0_0_0_yz_2 0_0_0_yy_2 0_0_0_x_2 0_0_0_x_3 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSS_2 SSSS_3 SSSP_2 SSSP_3 SSSD_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

INTEGRAL:1 : 0 : 1 : SSSS_1 Y SSSS_2 Y SPSS_1 Y 

0_0_0_0_1 0_0_0_0_2 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SSSS_1 SSSS_2 SPSS_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

INTEGRAL:1 : 1 : 0 : SSSS_1 Y SSSP_0 Y SSSP_1 Y SPSP_0 Y 

0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 

SSSS_1 SSSP_0 SSSP_1 SPSP_0 

SSSS_1 : 0_0_0_0_1 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SPSP_0 : 0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 

INTEGRAL:1 : 1 : 1 : SSSS_2 Y SSSP_1 Y SSSP_2 Y SPSP_1 Y 

0_0_0_0_2 0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_z_1 0_y_0_y_1 0_x_0_x_1 

SSSS_2 SSSP_1 SSSP_2 SPSP_1 

SSSS_2 : 0_0_0_0_2 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSP_1 : 0_z_0_z_1 0_y_0_y_1 0_x_0_x_1 

INTEGRAL:1 : 2 : 0 : SSSP_1 Y SSSD_0 Y SSSD_1 Y SPSD_0 Y 

0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SSSP_1 SSSD_0 SSSD_1 SPSD_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

INTEGRAL:1 : 2 : 1 : SSSP_2 Y SSSD_1 Y SSSD_2 Y SPSD_1 Y 

0_0_0_z_2 0_0_0_zz_1 0_0_0_zz_2 0_0_0_y_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yy_1 0_0_0_yy_2 0_0_0_x_2 0_0_0_xz_1 0_0_0_xz_2 0_0_0_xy_1 0_0_0_xy_2 0_0_0_xx_1 0_0_0_xx_2 0_z_0_zz_1 0_z_0_yz_1 0_z_0_xz_1 0_y_0_yy_1 0_y_0_xy_1 0_x_0_xx_1 

SSSP_2 SSSD_1 SSSD_2 SPSD_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_xz_1 0_y_0_yy_1 0_y_0_xy_1 0_x_0_xx_1 

INTEGRAL:2 : 1 : 0 : SSSP_0 Y SSSP_1 Y SPSS_1 Y SPSP_0 Y SPSP_1 Y SDSP_0 Y 

0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_zz_0_z_0 0_y_0_0_1 0_y_0_y_0 0_y_0_y_1 0_yz_0_z_0 0_yz_0_y_0 0_yy_0_y_0 0_x_0_0_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_x_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_x_0 

SSSP_0 SSSP_1 SPSS_1 SPSP_0 SPSP_1 SDSP_0 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SPSP_0 : 0_z_0_z_0 0_y_0_y_0 0_x_0_x_0 

SPSP_1 : 0_z_0_z_1 0_y_0_y_1 0_x_0_x_1 

SDSP_0 : 0_zz_0_z_0 0_yz_0_z_0 0_yz_0_y_0 0_yy_0_y_0 0_xz_0_z_0 0_xz_0_x_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_x_0 

INTEGRAL:2 : 2 : 0 : SSSD_0 N SSSD_1 N SPSP_1 Y SPSD_0 N SPSD_1 Y SDSD_0 Y 

0_0_0_zz_0 0_0_0_zz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_xz_0 0_z_0_xz_1 0_zz_0_zz_0 0_y_0_y_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xy_0 0_y_0_xy_1 0_yz_0_yz_0 0_yy_0_yy_0 0_x_0_x_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_xx_0 

SSSD_0 SSSD_1 SPSP_1 SPSD_0 SPSD_1 SDSD_0 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSP_1 : 0_z_0_z_1 0_y_0_y_1 0_x_0_x_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_xz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_xz_1 0_y_0_yy_1 0_y_0_xy_1 0_x_0_xx_1 

SDSD_0 : 0_zz_0_zz_0 0_yz_0_yz_0 0_yy_0_yy_0 0_xz_0_xz_0 0_xy_0_xy_0 0_xx_0_xx_0 

        }
    }
}


} // derirec namespace
