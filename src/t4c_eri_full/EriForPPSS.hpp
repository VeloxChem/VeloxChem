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
#include "EriHRRForPPSS.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPPSS(      T*                                 intsBuffer,
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

    BufferHostXY<T> bvals(3, ncpairs);

    CBoysFunc<T, 2> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(3, ncpairs);

    BufferHostXY<T> pbufSPSS0(3, ncpairs);

0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 
    BufferHostXY<T> pbufSPSS1(3, ncpairs);

0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 
    BufferHostXY<T> pbufSDSS0(6, ncpairs);

0_zz_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSPSS(3, ncpairs);

0_z_0_0 0_y_0_0 0_x_0_0 
    BufferHostXY<T> cbufSDSS(6, ncpairs);

0_zz_0_0 0_yz_0_0 0_yy_0_0 0_xz_0_0 0_xy_0_0 0_xx_0_0 
    BufferHostXY<T> cbufPPSS(9, ncpairs);

z_z_0_0 z_y_0_0 z_x_0_0 y_z_0_0 y_y_0_0 y_x_0_0 x_z_0_0 x_y_0_0 x_x_0_0 
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
0_z_0_0_0 0_zz_0_0_0 0_y_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 z_z_0_0_0 z_y_0_0_0 z_x_0_0_0 y_z_0_0_0 y_y_0_0_0 y_x_0_0_0 x_z_0_0_0 x_y_0_0_0 x_x_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

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

INTEGRAL:2 : 0 : 0 : SSSS_0 Y SSSS_1 Y SPSS_0 Y SPSS_1 Y SDSS_0 Y 

0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yz_0_0_0 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

SSSS_0 SSSS_1 SPSS_0 SPSS_1 SDSS_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SPSS_0 : 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SDSS_0 : 0_zz_0_0_0 0_yz_0_0_0 0_yy_0_0_0 0_xz_0_0_0 0_xy_0_0_0 0_xx_0_0_0 

        }
    }
}


} // derirec namespace
