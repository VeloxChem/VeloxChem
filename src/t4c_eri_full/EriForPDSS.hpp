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
#include "EriVRRForSPSS.hpp"
#include "EriVRRForSDSS.hpp"
#include "EriVRRForSFSS.hpp"
#include "EriHRRForPDSS.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPDSS(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions

    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();

    const auto ncpairs = ePosition - bPosition;

    // allocate Obara-Saika factors

    BufferHostMY<T, 3> osfacts(ncpairs);

    // allocate distances

    BufferHostMY<T, 3> rpq(ncpairs); 

    BufferHostMY<T, 3> rab(ncpairs);

    BufferHostMY<T, 3> rpb(ncpairs);

    BufferHostMY<T, 3> rwp(ncpairs);

    // allocate coordinates

    BufferHostMY<T, 3> rw(ncpairs); 

    // allocate Boys function data

    BufferHostX<T> bargs(ncpairs);

    BufferHostXY<T> bvals(4, ncpairs);

    CBoysFunc<T, 3> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(4, ncpairs);

    BufferHostXY<T> pbufSPSS0(3, ncpairs);

0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 
    BufferHostXY<T> pbufSPSS1(3, ncpairs);

0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 
    BufferHostXY<T> pbufSPSS2(3, ncpairs);

0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 
    BufferHostXY<T> pbufSDSS0(6, ncpairs);

0_zz_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 
    BufferHostXY<T> pbufSDSS1(4, ncpairs);

0_zz_0_0_1 0_yz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 
    BufferHostXY<T> pbufSFSS0(10, ncpairs);

0_zzz_0_0_0 0_yzz_0_0_0 0_yyz_0_0_0 0_yyy_0_0_0 0_xzz_0_0_0 0_xyz_0_0_0 0_xyy_0_0_0 0_xxz_0_0_0 0_xxy_0_0_0 0_xxx_0_0_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSS(6, ncpairs);

0_zz_0_0 0_yz_0_0 0_yy_0_0 0_xz_0_0 0_xy_0_0 0_xx_0_0 
    BufferHostXY<T> cbufSFSS(10, ncpairs);

0_zzz_0_0 0_yzz_0_0 0_yyz_0_0 0_yyy_0_0 0_xzz_0_0 0_xyz_0_0 0_xyy_0_0 0_xxz_0_0 0_xxy_0_0 0_xxx_0_0 
    BufferHostXY<T> cbufPDSS(18, ncpairs);

z_zz_0_0 z_yz_0_0 z_yy_0_0 z_xz_0_0 z_xy_0_0 z_xx_0_0 y_zz_0_0 y_yz_0_0 y_yy_0_0 y_xz_0_0 y_xy_0_0 y_xx_0_0 x_zz_0_0 x_yz_0_0 x_yy_0_0 x_xz_0_0 x_xy_0_0 x_xx_0_0 
    for (int32_t i = 0; i < nppairs; i++)
    {
        derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);

        derirec::compHostFactorsPartialZeta(fz, gtoPairBlock, bPosition, ePosition, i);

        for (int j = i; j < nppairs; j++)
        {
            derirec::compHostDistancesPQ(rpq, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsRho(frho, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsNorm(fnorm, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostCoordinatesW(rw, gtoPairBlock, bPosition, ePosition, i, j);

            if (i == j)
            {
                rwp.setZero();
            }
            else
            {
                derirec::compHostDistancesWT(rwp, rw, gtoPairBlock, bPosition, ePosition, i);
            }

            derirec::compHostBoysArguments(bargs, rpq, frho, ncpairs);

            bftable.compute(bvals, bargs);

            derirec::compHostVRRForSSSS(pbufSSSS, bvals, frho, fnorm, ncpairs);

signature:
0_zz_0_0_0 0_zzz_0_0_0 0_yz_0_0_0 0_yzz_0_0_0 0_yy_0_0_0 0_yyz_0_0_0 0_yyy_0_0_0 0_xz_0_0_0 0_xzz_0_0_0 0_xy_0_0_0 0_xyz_0_0_0 0_xyy_0_0_0 0_xx_0_0_0 0_xxz_0_0_0 0_xxy_0_0_0 0_xxx_0_0_0 z_zz_0_0_0 z_yz_0_0_0 z_yy_0_0_0 z_xz_0_0_0 z_xy_0_0_0 z_xx_0_0_0 y_zz_0_0_0 y_yz_0_0_0 y_yy_0_0_0 y_xz_0_0_0 y_xy_0_0_0 y_xx_0_0_0 x_zz_0_0_0 x_yz_0_0_0 x_yy_0_0_0 x_xz_0_0_0 x_xy_0_0_0 x_xx_0_0_0 

signature:
0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_zz_0_0_1 0_zzz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yz_0_0_1 0_yzz_0_0_0 0_yy_0_0_0 0_yy_0_0_1 0_yyz_0_0_0 0_yyy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xz_0_0_0 0_xzz_0_0_0 0_xy_0_0_0 0_xyz_0_0_0 0_xyy_0_0_0 0_xx_0_0_0 0_xx_0_0_1 0_xxz_0_0_0 0_xxy_0_0_0 0_xxx_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xx_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

INTEGRAL:1 : 0 : 0 : SSSS_0 Y SSSS_1 Y SPSS_0 Y 

0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

SSSS_0 SSSS_1 SPSS_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SPSS_0 : 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

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

INTEGRAL:2 : 0 : 0 : SSSS_0 Y SSSS_1 Y SPSS_0 Y SPSS_1 Y SDSS_0 Y 

0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

SSSS_0 SSSS_1 SPSS_0 SPSS_1 SDSS_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SPSS_0 : 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SDSS_0 : 0_zz_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

INTEGRAL:2 : 0 : 1 : SSSS_1 Y SSSS_2 Y SPSS_1 Y SPSS_2 Y SDSS_1 Y 

0_0_0_0_1 0_0_0_0_2 0_z_0_0_1 0_z_0_0_2 0_zz_0_0_1 0_y_0_0_1 0_y_0_0_2 0_yz_0_0_1 0_yy_0_0_1 0_x_0_0_1 0_x_0_0_2 0_xx_0_0_1 

SSSS_1 SSSS_2 SPSS_1 SPSS_2 SDSS_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SDSS_1 : 0_zz_0_0_1 0_yz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 

INTEGRAL:3 : 0 : 0 : SPSS_0 Y SPSS_1 Y SDSS_0 N SDSS_1 Y SFSS_0 Y 

0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_zz_0_0_1 0_zzz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yz_0_0_1 0_yzz_0_0_0 0_yy_0_0_0 0_yy_0_0_1 0_yyz_0_0_0 0_yyy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xzz_0_0_0 0_xyz_0_0_0 0_xyy_0_0_0 0_xx_0_0_0 0_xx_0_0_1 0_xxz_0_0_0 0_xxy_0_0_0 0_xxx_0_0_0 

SPSS_0 SPSS_1 SDSS_0 SDSS_1 SFSS_0 

SPSS_0 : 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SDSS_0 : 0_zz_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

SDSS_1 : 0_zz_0_0_1 0_yz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 

SFSS_0 : 0_zzz_0_0_0 0_yzz_0_0_0 0_yyz_0_0_0 0_yyy_0_0_0 0_xzz_0_0_0 0_xyz_0_0_0 0_xyy_0_0_0 0_xxz_0_0_0 0_xxy_0_0_0 0_xxx_0_0_0 

        }
    }
}


} // derirec namespace
