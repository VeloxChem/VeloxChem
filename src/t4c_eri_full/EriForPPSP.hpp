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
#include "EriVRRForSDSP.hpp"
#include "EriHRRForPPSP.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPPSP(      T*                                 intsBuffer,
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

    BufferHostXY<T> bvals(4, ncpairs);

    CBoysFunc<T, 3> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(4, ncpairs);

    BufferHostXY<T> pbufSSSP0(3, ncpairs);

0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 
    BufferHostXY<T> pbufSSSP1(3, ncpairs);

0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 
    BufferHostXY<T> pbufSSSP2(3, ncpairs);

0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 
    BufferHostXY<T> pbufSPSS1(3, ncpairs);

0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 
    BufferHostXY<T> pbufSPSP0(9, ncpairs);

0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 
    BufferHostXY<T> pbufSPSP1(9, ncpairs);

0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSDSP0(18, ncpairs);

0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSPSP(9, ncpairs);

0_z_0_z 0_z_0_y 0_z_0_x 0_y_0_z 0_y_0_y 0_y_0_x 0_x_0_z 0_x_0_y 0_x_0_x 
    BufferHostXY<T> cbufSDSP(18, ncpairs);

0_zz_0_z 0_zz_0_y 0_zz_0_x 0_yz_0_z 0_yz_0_y 0_yz_0_x 0_yy_0_z 0_yy_0_y 0_yy_0_x 0_xz_0_z 0_xz_0_y 0_xz_0_x 0_xy_0_z 0_xy_0_y 0_xy_0_x 0_xx_0_z 0_xx_0_y 0_xx_0_x 
    BufferHostXY<T> cbufPPSP(27, ncpairs);

z_z_0_z z_z_0_y z_z_0_x z_y_0_z z_y_0_y z_y_0_x z_x_0_z z_x_0_y z_x_0_x y_z_0_z y_z_0_y y_z_0_x y_y_0_z y_y_0_y y_y_0_x y_x_0_z y_x_0_y y_x_0_x x_z_0_z x_z_0_y x_z_0_x x_y_0_z x_y_0_y x_y_0_x x_x_0_z x_x_0_y x_x_0_x 
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
0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 z_z_0_z_0 z_z_0_y_0 z_z_0_x_0 z_y_0_z_0 z_y_0_y_0 z_y_0_x_0 z_x_0_z_0 z_x_0_y_0 z_x_0_x_0 y_z_0_z_0 y_z_0_y_0 y_z_0_x_0 y_y_0_z_0 y_y_0_y_0 y_y_0_x_0 y_x_0_z_0 y_x_0_y_0 y_x_0_x_0 x_z_0_z_0 x_z_0_y_0 x_z_0_x_0 x_y_0_z_0 x_y_0_y_0 x_y_0_x_0 x_x_0_z_0 x_x_0_y_0 x_x_0_x_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

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

INTEGRAL:1 : 0 : 1 : SSSS_1 Y SSSS_2 Y SPSS_1 Y 

0_0_0_0_1 0_0_0_0_2 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SSSS_1 SSSS_2 SPSS_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

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

INTEGRAL:2 : 1 : 0 : SSSP_0 Y SSSP_1 Y SPSS_1 Y SPSP_0 Y SPSP_1 Y SDSP_0 Y 

0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

SSSP_0 SSSP_1 SPSS_1 SPSP_0 SPSP_1 SDSP_0 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SPSP_0 : 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SDSP_0 : 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

        }
    }
}


} // derirec namespace
