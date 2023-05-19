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
#include "EriHRRForSSPP.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostSSPP(      T*                                 intsBuffer,
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

    BufferHostMY<T, 3> rcd(ncpairs);

    BufferHostMY<T, 3> rqd(ncpairs);

    BufferHostMY<T, 3> rwq(ncpairs);

    // allocate coordinates

    BufferHostMY<T, 3> rw(ncpairs); 

    // allocate Boys function data

    BufferHostX<T> bargs(ncpairs);

    BufferHostXY<T> bvals(3, ncpairs);

    CBoysFunc<T, 2> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(3, ncpairs);

    BufferHostXY<T> pbufSSSP0(3, ncpairs);

0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 
    BufferHostXY<T> pbufSSSP1(3, ncpairs);

0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 
    BufferHostXY<T> pbufSSSD0(6, ncpairs);

0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSSSP(3, ncpairs);

0_0_0_z 0_0_0_y 0_0_0_x 
    BufferHostXY<T> cbufSSSD(6, ncpairs);

0_0_0_zz 0_0_0_yz 0_0_0_yy 0_0_0_xz 0_0_0_xy 0_0_0_xx 
    BufferHostXY<T> cbufSSPP(9, ncpairs);

0_0_z_z 0_0_z_y 0_0_z_x 0_0_y_z 0_0_y_y 0_0_y_x 0_0_x_z 0_0_x_y 0_0_x_x 
    for (int32_t i = 0; i < nppairs; i++)
    {
        for (int j = i; j < nppairs; j++)
        {
            derirec::compHostDistancesPQ(rpq, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsRho(frho, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsNorm(fnorm, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsPartialZeta(fe, gtoPairBlock, bPosition, ePosition, j);

            derirec::compHostDistancesPT(rqd, gtoPairBlock, bPosition, ePosition, j);

            derirec::compHostCoordinatesW(rw, gtoPairBlock, bPosition, ePosition, i, j);

            if (i == j)
            {
                rwq.setZero();
            }
            else
            {
                derirec::compHostDistancesWT(rwq, rw, gtoPairBlock, bPosition, ePosition, j);
            }

            derirec::compHostBoysArguments(bargs, rpq, frho, ncpairs);

            bftable.compute(bvals, bargs);

            derirec::compHostVRRForSSSS(pbufSSSS, bvals, frho, fnorm, ncpairs);

signature:
0_0_0_z_0 0_0_0_zz_0 0_0_0_y_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 0_0_z_z_0 0_0_z_y_0 0_0_z_x_0 0_0_y_z_0 0_0_y_y_0 0_0_y_x_0 0_0_x_z_0 0_0_x_y_0 0_0_x_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

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

INTEGRAL:0 : 2 : 0 : SSSS_0 Y SSSS_1 Y SSSP_0 Y SSSP_1 Y SSSD_0 Y 

0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSS_0 SSSS_1 SSSP_0 SSSP_1 SSSD_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

        }
    }
}


} // derirec namespace
