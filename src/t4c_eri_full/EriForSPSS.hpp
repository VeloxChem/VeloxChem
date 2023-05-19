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

namespace derirec { // derirec namespace

template <typename T>
auto
compHostSPSS(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions

    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();

    const auto ncpairs = ePosition - bPosition;

    // allocate Obara-Saika factors

    BufferHostMY<T, 1> osfacts(ncpairs);

    // allocate distances

    BufferHostMY<T, 3> rpq(ncpairs); 

    BufferHostMY<T, 3> rpb(ncpairs);

    BufferHostMY<T, 3> rwp(ncpairs);

    // allocate coordinates

    BufferHostMY<T, 3> rw(ncpairs); 

    // allocate Boys function data

    BufferHostX<T> bargs(ncpairs);

    BufferHostXY<T> bvals(2, ncpairs);

    CBoysFunc<T, 1> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(2, ncpairs);

    // Contracted integral buffers

    auto cbufSPSS = BufferHostXY<T>::Zero(3, ncpairs);

0_z_0_0 0_y_0_0 0_x_0_0 
    for (int32_t i = 0; i < nppairs; i++)
    {
        derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);

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
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

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

        }
    }
}


} // derirec namespace
