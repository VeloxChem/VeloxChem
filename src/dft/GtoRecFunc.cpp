//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoRecFunc.hpp"

namespace gtorec { // gtorec namespace
    
void
compGtoTypeSForLDA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CGtoContainer&        gtoContainer,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitives data
    
    auto pexp = gtoContainer.getExponents(iGtoBlock);
    
    auto pnorm = gtoContainer.getNormFactors(iGtoBlock);
    
    // primitive s-type GTOs
    
    auto g_0 = gtoValues.data(0);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g_0, pexp, pnorm, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        g_0[i] = pnorm[i] * std::exp(-pexp[i] * (rx[i] * rx[i] + ry[i] * ry[i]
                                     
                                     + rz[i] * rz[i]));
    }
}

void
compGtoTypePForLDA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive s-type GTOs
    
    auto g_0 = gtoValues.data(0);
    
    // primitive p-type GTOs
    
    auto g_x = gtoValues.data(1);
    
    auto g_y = gtoValues.data(2);
    
    auto g_z = gtoValues.data(3);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g_0, g_x, g_y, g_z, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        double fact = g_0[i];
        
        g_x[i] = rx[i] * fact;
        
        g_y[i] = ry[i] * fact;
        
        g_z[i] = rz[i] * fact;
    }
}
    
void
compGtoTypeDForLDA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive p-type GTOs
    
    auto g_x = gtoValues.data(1);
    
    auto g_y = gtoValues.data(2);
    
    auto g_z = gtoValues.data(3);
    
    // primitive d-type GTOs
    
    auto g_xx = gtoValues.data(4);
    
    auto g_xy = gtoValues.data(5);
    
    auto g_xz = gtoValues.data(6);
    
    auto g_yy = gtoValues.data(7);
    
    auto g_yz = gtoValues.data(8);
    
    auto g_zz = gtoValues.data(9);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g_x, g_y, g_z, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,\
                             rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // x leading index
        
        double fact = rx[i];
        
        g_xx[i] = fact * g_x[i];
        
        g_xy[i] = fact * g_y[i];
        
        g_xz[i] = fact * g_z[i];
        
        // y leading index
        
        fact = ry[i];
        
        g_yy[i] = fact * g_y[i];
        
        g_yz[i] = fact * g_z[i];
        
        // z leading index
        
        g_zz[i] = rz[i] * g_z[i];
    }
}
    
void
compGtoTypeFForLDA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive d-type GTOs
    
    auto g_xx = gtoValues.data(4);
    
    auto g_xy = gtoValues.data(5);
    
    auto g_xz = gtoValues.data(6);
    
    auto g_yy = gtoValues.data(7);
    
    auto g_yz = gtoValues.data(8);
    
    auto g_zz = gtoValues.data(9);
    
    // primitive f-type GTOs
    
    auto g_xxx = gtoValues.data(10);
    
    auto g_xxy = gtoValues.data(11);
    
    auto g_xxz = gtoValues.data(12);
    
    auto g_xyy = gtoValues.data(13);
    
    auto g_xyz = gtoValues.data(14);
    
    auto g_xzz = gtoValues.data(15);
    
    auto g_yyy = gtoValues.data(16);
    
    auto g_yyz = gtoValues.data(17);
    
    auto g_yzz = gtoValues.data(18);
    
    auto g_zzz = gtoValues.data(19);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, g_xxx, g_xxy,\
                             g_xxz, g_xyy, g_xyz, g_xzz, g_yyy, g_yyz, g_yzz,\
                             g_zzz, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // x leading index
        
        double fact = rx[i];
        
        g_xxx[i] = fact * g_xx[i];
        
        g_xxy[i] = fact * g_xy[i];
        
        g_xxz[i] = fact * g_xz[i];
        
        g_xyy[i] = fact * g_yy[i];
        
        g_xyz[i] = fact * g_yz[i];
        
        g_xzz[i] = fact * g_zz[i];
        
        // y leading index
        
        fact = ry[i];
        
        g_yyy[i] = fact * g_yy[i];
        
        g_yyz[i] = fact * g_yz[i];
        
        g_yzz[i] = fact * g_zz[i];
        
        // z leading index
        
        g_zzz[i] = rz[i] * g_zz[i];
    }
}
    
void
compGtoTypeGForLDA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive f-type GTOs
    
    auto g_xxx = gtoValues.data(10);
    
    auto g_xxy = gtoValues.data(11);
    
    auto g_xxz = gtoValues.data(12);
    
    auto g_xyy = gtoValues.data(13);
    
    auto g_xyz = gtoValues.data(14);
    
    auto g_xzz = gtoValues.data(15);
    
    auto g_yyy = gtoValues.data(16);
    
    auto g_yyz = gtoValues.data(17);
    
    auto g_yzz = gtoValues.data(18);
    
    auto g_zzz = gtoValues.data(19);
    
    // primitive g-type GTOs
    
    auto g_xxxx = gtoValues.data(20);
    
    auto g_xxxy = gtoValues.data(21);
    
    auto g_xxxz = gtoValues.data(22);
    
    auto g_xxyy = gtoValues.data(23);
    
    auto g_xxyz = gtoValues.data(24);
    
    auto g_xxzz = gtoValues.data(25);
    
    auto g_xyyy = gtoValues.data(26);
    
    auto g_xyyz = gtoValues.data(27);
    
    auto g_xyzz = gtoValues.data(28);
    
    auto g_xzzz = gtoValues.data(29);
    
    auto g_yyyy = gtoValues.data(30);
    
    auto g_yyyz = gtoValues.data(31);
    
    auto g_yyzz = gtoValues.data(32);
    
    auto g_yzzz = gtoValues.data(33);
    
    auto g_zzzz = gtoValues.data(34);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g_xxx, g_xxy, g_xxz, g_xyy, g_xyz, g_xzz, g_yyy,\
                             g_yyz, g_yzz, g_zzz, g_xxxx, g_xxxy, g_xxxz,\
                             g_xxyy, g_xxyz, g_xxzz, g_xyyy, g_xyyz, g_xyzz,\
                             g_xzzz, g_yyyy, g_yyyz, g_yyzz, g_yzzz, g_zzzz,\
                             rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // x leading index
        
        double fact = rx[i];
        
        g_xxxx[i] = fact * g_xxx[i];
        
        g_xxxy[i] = fact * g_xxy[i];
        
        g_xxxz[i] = fact * g_xxz[i];
        
        g_xxyy[i] = fact * g_xyy[i];
        
        g_xxyz[i] = fact * g_xyz[i];
        
        g_xxzz[i] = fact * g_xzz[i];
        
        g_xyyy[i] = fact * g_yyy[i];
        
        g_xyyz[i] = fact * g_yyz[i];
        
        g_xyzz[i] = fact * g_yzz[i];
        
        g_xzzz[i] = fact * g_zzz[i];
        
        // y leading index
        
        fact = ry[i];
        
        g_yyyy[i] = fact * g_yyy[i];
        
        g_yyyz[i] = fact * g_yyz[i];
        
        g_yyzz[i] = fact * g_yzz[i];
        
        g_yzzz[i] = fact * g_zzz[i];
        
        // z leading index
        
        g_zzzz[i] = rz[i] * g_zzz[i];
    }
}

void
compGtoTypeSForGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CGtoContainer&        gtoContainer,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitives data
    
    auto pexp = gtoContainer.getExponents(iGtoBlock);
    
    auto pnorm = gtoContainer.getNormFactors(iGtoBlock);
    
    // primitive s-type GTOs
    
    auto g0_0 = gtoValues.data(0);
    
    auto gx_0 = gtoValues.data(1);
    
    auto gy_0 = gtoValues.data(2);
    
    auto gz_0 = gtoValues.data(3);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_0, gx_0, gy_0, gz_0, pexp, pnorm, rx, ry,\
                             rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        g0_0[i] = pnorm[i] * std::exp(-pexp[i] * (rx[i] * rx[i] + ry[i] * ry[i]
                                     
                                      + rz[i] * rz[i]));
        
        double fact = -2.0 * pexp[i] * g0_0[i];
        
        gx_0[i] = fact * rx[i];
        
        gy_0[i] = fact * ry[i];
        
        gz_0[i] = fact * rz[i];
    }
}

void
compGtoTypePForGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive s-type GTOs
    
    auto g0_0 = gtoValues.data(0);
    
    auto gx_0 = gtoValues.data(1);
    
    auto gy_0 = gtoValues.data(2);
    
    auto gz_0 = gtoValues.data(3);
    
    // primitive p-type GTOs
    
    auto g0_x = gtoValues.data(4);
    
    auto gx_x = gtoValues.data(5);
    
    auto gy_x = gtoValues.data(6);
    
    auto gz_x = gtoValues.data(7);
    
    auto g0_y = gtoValues.data(8);
    
    auto gx_y = gtoValues.data(9);
    
    auto gy_y = gtoValues.data(10);
    
    auto gz_y = gtoValues.data(11);
    
    auto g0_z = gtoValues.data(12);
    
    auto gx_z = gtoValues.data(13);
    
    auto gy_z = gtoValues.data(14);
    
    auto gz_z = gtoValues.data(15);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_0, gx_0, gy_0, gz_0, g0_x, gx_x, gy_x, gz_x,\
                             g0_y, gx_y, gy_y, gz_y, g0_z, gx_z, gy_z, gz_z,\
                             rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // x component
        
        double fact = rx[i];
        
        g0_x[i] = fact * g0_0[i];
        
        gx_x[i] = g0_0[i] + fact * gx_0[i];
        
        gy_x[i] = fact * gy_0[i];
        
        gz_x[i] = fact * gz_0[i];
        
        // y component
        
        fact = ry[i];
        
        g0_y[i] = fact * g0_0[i];
        
        gx_y[i] = fact * gx_0[i];
        
        gy_y[i] = g0_0[i] + fact * gy_0[i];
        
        gz_y[i] = fact * gz_0[i];
        
        // z component
        
        fact = rz[i];
        
        g0_z[i] = fact * g0_0[i];
        
        gx_z[i] = fact * gx_0[i];
        
        gy_z[i] = fact * gy_0[i];
        
        gz_z[i] = g0_0[i] + fact * gz_0[i];
    }
}
    
void
compGtoTypeDForGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive p-type GTOs
    
    auto g0_x = gtoValues.data(4);
    
    auto gx_x = gtoValues.data(5);
    
    auto gy_x = gtoValues.data(6);
    
    auto gz_x = gtoValues.data(7);
    
    auto g0_y = gtoValues.data(8);
    
    auto gx_y = gtoValues.data(9);
    
    auto gy_y = gtoValues.data(10);
    
    auto gz_y = gtoValues.data(11);
    
    auto g0_z = gtoValues.data(12);
    
    auto gx_z = gtoValues.data(13);
    
    auto gy_z = gtoValues.data(14);
    
    auto gz_z = gtoValues.data(15);
    
    // primitive d-type GTOs
    
    auto g0_xx = gtoValues.data(16);
    
    auto gx_xx = gtoValues.data(17);
    
    auto gy_xx = gtoValues.data(18);
    
    auto gz_xx = gtoValues.data(19);
    
    auto g0_xy = gtoValues.data(20);
    
    auto gx_xy = gtoValues.data(21);
    
    auto gy_xy = gtoValues.data(22);
    
    auto gz_xy = gtoValues.data(23);
    
    auto g0_xz = gtoValues.data(24);
    
    auto gx_xz = gtoValues.data(25);
    
    auto gy_xz = gtoValues.data(26);
    
    auto gz_xz = gtoValues.data(27);
    
    auto g0_yy = gtoValues.data(28);
    
    auto gx_yy = gtoValues.data(29);
    
    auto gy_yy = gtoValues.data(30);
    
    auto gz_yy = gtoValues.data(31);
    
    auto g0_yz = gtoValues.data(32);
    
    auto gx_yz = gtoValues.data(33);
    
    auto gy_yz = gtoValues.data(34);
    
    auto gz_yz = gtoValues.data(35);
   
    auto g0_zz = gtoValues.data(36);
    
    auto gx_zz = gtoValues.data(37);
    
    auto gy_zz = gtoValues.data(38);
    
    auto gz_zz = gtoValues.data(39);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_x, gx_x, gy_x, gz_x, g0_y, gx_y, gy_y, gz_y,\
                             g0_z, gx_z, gy_z, gz_z, g0_xx, gx_xx, gy_xx,\
                             gz_xx, g0_xy, gx_xy, gy_xy, gz_xy, g0_xz, gx_xz,\
                             gy_xz, gz_xz, g0_yy, gx_yy, gy_yy, gz_yy, g0_yz,\
                             gx_yz, gy_yz, gz_yz, g0_zz, gx_zz, gy_zz, gz_zz,\
                             rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xx component
        
        double fact = rx[i];

        g0_xx[i] = fact * g0_x[i];
        
        gx_xx[i] = g0_x[i] + fact * gx_x[i];

        gy_xx[i] = fact * gy_x[i];

        gz_xx[i] = fact * gz_x[i];
        
        // xy index component
    
        fact = ry[i];
        
        g0_xy[i] = fact * g0_x[i];
        
        gx_xy[i] = fact * gx_x[i];
        
        gy_xy[i] = g0_x[i] + fact * gy_x[i];
        
        gz_xy[i] = fact * gz_x[i];
        
        // xz index component
        
        fact = rz[i];
        
        g0_xz[i] = fact * g0_x[i];
        
        gx_xz[i] = fact * gx_x[i];
        
        gy_xz[i] = fact * gy_x[i];
        
        gz_xz[i] = g0_x[i] + fact * gz_x[i];
        
        // yy component
        
        fact = ry[i];
        
        g0_yy[i] = fact * g0_y[i];
        
        gx_yy[i] = fact * gx_y[i];
        
        gy_yy[i] = g0_y[i] + fact * gy_y[i];
        
        gz_yy[i] = fact * gz_y[i];
        
        // yz component
        
        fact = rz[i];
        
        g0_yz[i] = fact * g0_y[i];
        
        gx_yz[i] = fact * gx_y[i];
        
        gy_yz[i] = fact * gy_y[i];
        
        gz_yz[i] = g0_y[i] + fact * gz_y[i];
        
        // zz component
        
        fact = rz[i];
        
        g0_zz[i] = fact * g0_z[i];
        
        gx_zz[i] = fact * gx_z[i];
        
        gy_zz[i] = fact * gy_z[i];
        
        gz_zz[i] = g0_z[i] + fact * gz_z[i];
    }
}
    
void
compGtoTypeFForGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive d-type GTOs
    
    auto g0_xx = gtoValues.data(16);
    
    auto gx_xx = gtoValues.data(17);
    
    auto gy_xx = gtoValues.data(18);
    
    auto gz_xx = gtoValues.data(19);
    
    auto g0_yy = gtoValues.data(28);
    
    auto gx_yy = gtoValues.data(29);
    
    auto gy_yy = gtoValues.data(30);
    
    auto gz_yy = gtoValues.data(31);
    
    auto g0_yz = gtoValues.data(32);
    
    auto gx_yz = gtoValues.data(33);
    
    auto gy_yz = gtoValues.data(34);
    
    auto gz_yz = gtoValues.data(35);
   
    auto g0_zz = gtoValues.data(36);
    
    auto gx_zz = gtoValues.data(37);
    
    auto gy_zz = gtoValues.data(38);
    
    auto gz_zz = gtoValues.data(39);
    
    // primitive f-type GTOs
    
    auto g0_xxx = gtoValues.data(40);
    
    auto gx_xxx = gtoValues.data(41);
    
    auto gy_xxx = gtoValues.data(42);
    
    auto gz_xxx = gtoValues.data(43);
    
    auto g0_xxy = gtoValues.data(44);
    
    auto gx_xxy = gtoValues.data(45);
    
    auto gy_xxy = gtoValues.data(46);
    
    auto gz_xxy = gtoValues.data(47);
    
    auto g0_xxz = gtoValues.data(48);
    
    auto gx_xxz = gtoValues.data(49);
    
    auto gy_xxz = gtoValues.data(50);
    
    auto gz_xxz = gtoValues.data(51);
    
    auto g0_xyy = gtoValues.data(52);
    
    auto gx_xyy = gtoValues.data(53);
    
    auto gy_xyy = gtoValues.data(54);
    
    auto gz_xyy = gtoValues.data(55);
    
    auto g0_xyz = gtoValues.data(56);
    
    auto gx_xyz = gtoValues.data(57);
    
    auto gy_xyz = gtoValues.data(58);
    
    auto gz_xyz = gtoValues.data(59);
    
    auto g0_xzz = gtoValues.data(60);
    
    auto gx_xzz = gtoValues.data(61);
    
    auto gy_xzz = gtoValues.data(62);
    
    auto gz_xzz = gtoValues.data(63);
    
    auto g0_yyy = gtoValues.data(64);
    
    auto gx_yyy = gtoValues.data(65);
    
    auto gy_yyy = gtoValues.data(66);
    
    auto gz_yyy = gtoValues.data(67);
    
    auto g0_yyz = gtoValues.data(68);
    
    auto gx_yyz = gtoValues.data(69);
    
    auto gy_yyz = gtoValues.data(70);
    
    auto gz_yyz = gtoValues.data(71);
    
    auto g0_yzz = gtoValues.data(72);
    
    auto gx_yzz = gtoValues.data(73);
    
    auto gy_yzz = gtoValues.data(74);
    
    auto gz_yzz = gtoValues.data(75);
    
    auto g0_zzz = gtoValues.data(76);
    
    auto gx_zzz = gtoValues.data(77);
    
    auto gy_zzz = gtoValues.data(78);
    
    auto gz_zzz = gtoValues.data(79);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_xx, gx_xx, gy_xx, gz_xx, g0_yy, gx_yy, gy_yy,\
                             gz_yy, g0_yz, gx_yz, gy_yz, gz_yz, g0_zz, gx_zz,\
                             gy_zz, gz_zz, g0_xxx, gx_xxx, gy_xxx, gz_xxx,\
                             g0_xxy, gx_xxy, gy_xxy, gz_xxy, g0_xxz, gx_xxz,\
                             gy_xxz, gz_xxz, g0_xyy, gx_xyy, gy_xyy, gz_xyy,\
                             g0_xyz, gx_xyz, gy_xyz, gz_xyz, g0_xzz, gx_xzz,\
                             gy_xzz, gz_xzz, g0_yyy, gx_yyy, gy_yyy, gz_yyy,\
                             g0_yyz, gx_yyz, gy_yyz, gz_yyz, g0_yzz, gx_yzz,\
                             gy_yzz, gz_yzz, g0_zzz, gx_zzz, gy_zzz, gz_zzz,\
                             rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xxx component
        
        double fact = rx[i];

        g0_xxx[i] = fact * g0_xx[i];
        
        gx_xxx[i] = g0_xx[i] + fact * gx_xx[i];

        gy_xxx[i] = fact * gy_xx[i];

        gz_xxx[i] = fact * gz_xx[i];
        
        // xxy index component
    
        fact = ry[i];
        
        g0_xxy[i] = fact * g0_xx[i];
        
        gx_xxy[i] = fact * gx_xx[i];
        
        gy_xxy[i] = g0_xx[i] + fact * gy_xx[i];
        
        gz_xxy[i] = fact * gz_xx[i];
        
        // xxz index component
        
        fact = rz[i];
        
        g0_xxz[i] = fact * g0_xx[i];
        
        gx_xxz[i] = fact * gx_xx[i];
        
        gy_xxz[i] = fact * gy_xx[i];
        
        gz_xxz[i] = g0_xx[i] + fact * gz_xx[i];
        
        // xyy component
        
        fact = rx[i];
        
        g0_xyy[i] = fact * g0_yy[i];
        
        gx_xyy[i] = g0_yy[i] + fact * gx_yy[i];
        
        gy_xyy[i] = fact * gy_yy[i];
        
        gz_xyy[i] = fact * gz_yy[i];
        
        // xyz component
        
        fact = rx[i];
        
        g0_xyz[i] = fact * g0_yz[i];
        
        gx_xyz[i] = g0_yz[i] + fact * gx_yz[i];
        
        gy_xyz[i] = fact * gy_yz[i];
        
        gz_xyz[i] = fact * gz_yz[i];
        
        // xzz component
        
        fact = rx[i];
        
        g0_xzz[i] = fact * g0_zz[i];
        
        gx_xzz[i] = g0_zz[i] + fact * gx_zz[i];
        
        gy_xzz[i] = fact * gy_zz[i];
        
        gz_xzz[i] = fact * gz_zz[i];
        
        // yyy component
        
        fact = ry[i];
        
        g0_yyy[i] = fact * g0_yy[i];
        
        gx_yyy[i] = fact * gx_yy[i];
        
        gy_yyy[i] = g0_yy[i] + fact * gy_yy[i];
        
        gz_yyy[i] = fact * gz_yy[i];
        
        // yyz component
        
        fact = rz[i];
        
        g0_yyz[i] = fact * g0_yy[i];
        
        gx_yyz[i] = fact * gx_yy[i];
        
        gy_yyz[i] = fact * gy_yy[i];
        
        gz_yyz[i] = g0_yy[i] + fact * gz_yy[i];
        
        // yzz component
        
        fact = ry[i];
        
        g0_yzz[i] = fact * g0_zz[i];
        
        gx_yzz[i] = fact * gx_zz[i];
        
        gy_yzz[i] = g0_zz[i] + fact * gy_zz[i];
        
        gz_yzz[i] = fact * gz_zz[i];
        
        // zzz component
        
        fact = rz[i];
        
        g0_zzz[i] = fact * g0_zz[i];
        
        gx_zzz[i] = fact * gx_zz[i];
        
        gy_zzz[i] = fact * gy_zz[i];
        
        gz_zzz[i] = g0_zz[i] + fact * gz_zz[i];
    }
}
    
void
compGtoTypeGForGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive f-type GTOs
    
    auto g0_xxx = gtoValues.data(40);
    
    auto gx_xxx = gtoValues.data(41);
    
    auto gy_xxx = gtoValues.data(42);
    
    auto gz_xxx = gtoValues.data(43);
    
    auto g0_xyy = gtoValues.data(52);
    
    auto gx_xyy = gtoValues.data(53);
    
    auto gy_xyy = gtoValues.data(54);
    
    auto gz_xyy = gtoValues.data(55);
    
    auto g0_xyz = gtoValues.data(56);
    
    auto gx_xyz = gtoValues.data(57);
    
    auto gy_xyz = gtoValues.data(58);
    
    auto gz_xyz = gtoValues.data(59);
    
    auto g0_xzz = gtoValues.data(60);
    
    auto gx_xzz = gtoValues.data(61);
    
    auto gy_xzz = gtoValues.data(62);
    
    auto gz_xzz = gtoValues.data(63);
    
    auto g0_yyy = gtoValues.data(64);
    
    auto gx_yyy = gtoValues.data(65);
    
    auto gy_yyy = gtoValues.data(66);
    
    auto gz_yyy = gtoValues.data(67);
    
    auto g0_yyz = gtoValues.data(68);
    
    auto gx_yyz = gtoValues.data(69);
    
    auto gy_yyz = gtoValues.data(70);
    
    auto gz_yyz = gtoValues.data(71);
    
    auto g0_yzz = gtoValues.data(72);
    
    auto gx_yzz = gtoValues.data(73);
    
    auto gy_yzz = gtoValues.data(74);
    
    auto gz_yzz = gtoValues.data(75);
    
    auto g0_zzz = gtoValues.data(76);
    
    auto gx_zzz = gtoValues.data(77);
    
    auto gy_zzz = gtoValues.data(78);
    
    auto gz_zzz = gtoValues.data(79);
    
    // primitive g-type GTOs
    
    auto g0_xxxx = gtoValues.data(80);
    
    auto gx_xxxx = gtoValues.data(81);
    
    auto gy_xxxx = gtoValues.data(82);
    
    auto gz_xxxx = gtoValues.data(83);
    
    auto g0_xxxy = gtoValues.data(84);
    
    auto gx_xxxy = gtoValues.data(85);
    
    auto gy_xxxy = gtoValues.data(86);
    
    auto gz_xxxy = gtoValues.data(87);
    
    auto g0_xxxz = gtoValues.data(88);
    
    auto gx_xxxz = gtoValues.data(89);
    
    auto gy_xxxz = gtoValues.data(90);
    
    auto gz_xxxz = gtoValues.data(91);
    
    auto g0_xxyy = gtoValues.data(92);
    
    auto gx_xxyy = gtoValues.data(93);
    
    auto gy_xxyy = gtoValues.data(94);
    
    auto gz_xxyy = gtoValues.data(95);
    
    auto g0_xxyz = gtoValues.data(96);
    
    auto gx_xxyz = gtoValues.data(97);
    
    auto gy_xxyz = gtoValues.data(98);
    
    auto gz_xxyz = gtoValues.data(99);
    
    auto g0_xxzz = gtoValues.data(100);
    
    auto gx_xxzz = gtoValues.data(101);
    
    auto gy_xxzz = gtoValues.data(102);
    
    auto gz_xxzz = gtoValues.data(103);
    
    auto g0_xyyy = gtoValues.data(104);
    
    auto gx_xyyy = gtoValues.data(105);
    
    auto gy_xyyy = gtoValues.data(106);
    
    auto gz_xyyy = gtoValues.data(107);
    
    auto g0_xyyz = gtoValues.data(108);
    
    auto gx_xyyz = gtoValues.data(109);
    
    auto gy_xyyz = gtoValues.data(110);
    
    auto gz_xyyz = gtoValues.data(111);
    
    auto g0_xyzz = gtoValues.data(112);
    
    auto gx_xyzz = gtoValues.data(113);
    
    auto gy_xyzz = gtoValues.data(114);
    
    auto gz_xyzz = gtoValues.data(115);
    
    auto g0_xzzz = gtoValues.data(116);
    
    auto gx_xzzz = gtoValues.data(117);
    
    auto gy_xzzz = gtoValues.data(118);
    
    auto gz_xzzz = gtoValues.data(119);
    
    auto g0_yyyy = gtoValues.data(120);
    
    auto gx_yyyy = gtoValues.data(121);
    
    auto gy_yyyy = gtoValues.data(122);
    
    auto gz_yyyy = gtoValues.data(123);
    
    auto g0_yyyz = gtoValues.data(124);
    
    auto gx_yyyz = gtoValues.data(125);
    
    auto gy_yyyz = gtoValues.data(126);
    
    auto gz_yyyz = gtoValues.data(127);
    
    auto g0_yyzz = gtoValues.data(128);
    
    auto gx_yyzz = gtoValues.data(129);
    
    auto gy_yyzz = gtoValues.data(130);
    
    auto gz_yyzz = gtoValues.data(131);
    
    auto g0_yzzz = gtoValues.data(132);
    
    auto gx_yzzz = gtoValues.data(133);
    
    auto gy_yzzz = gtoValues.data(134);
    
    auto gz_yzzz = gtoValues.data(135);
    
    auto g0_zzzz = gtoValues.data(136);
    
    auto gx_zzzz = gtoValues.data(137);
    
    auto gy_zzzz = gtoValues.data(138);
    
    auto gz_zzzz = gtoValues.data(139);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_xxx, gx_xxx, gy_xxx, gz_xxx, g0_xyy, gx_xyy,\
                             gy_xyy, gz_xyy, g0_xyz, gx_xyz, gy_xyz, gz_xyz,\
                             g0_xzz, gx_xzz, gy_xzz, gz_xzz, g0_yyy, gx_yyy,\
                             gy_yyy, gz_yyy, g0_yyz, gx_yyz, gy_yyz, gz_yyz,\
                             g0_yzz, gx_yzz, gy_yzz, gz_yzz, g0_zzz, gx_zzz,\
                             gy_zzz, gz_zzz, g0_xxxx, gx_xxxx, gy_xxxx,\
                             gz_xxxx, g0_xxxy, gx_xxxy, gy_xxxy, gz_xxxy,\
                             g0_xxxz, gx_xxxz, gy_xxxz, gz_xxxz, g0_xxyy,\
                             gx_xxyy, gy_xxyy, gz_xxyy, g0_xxyz, gx_xxyz,\
                             gy_xxyz, gz_xxyz, g0_xxzz, gx_xxzz, gy_xxzz,\
                             gz_xxzz, g0_xyyy, gx_xyyy, gy_xyyy, gz_xyyy,\
                             g0_xyyz, gx_xyyz, gy_xyyz, gz_xyyz, g0_xyzz,\
                             gx_xyzz, gy_xyzz, gz_xyzz, g0_xzzz, gx_xzzz,\
                             gy_xzzz, gz_xzzz, g0_yyyy, gx_yyyy, gy_yyyy,\
                             gz_yyyy, g0_yyyz, gx_yyyz, gy_yyyz, gz_yyyz,\
                             g0_yyzz, gx_yyzz, gy_yyzz, gz_yyzz, g0_yzzz,\
                             gx_yzzz, gy_yzzz, gz_yzzz, g0_zzzz, gx_zzzz,\
                             gy_zzzz, gz_zzzz, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xxxx component
        
        double fact = rx[i];

        g0_xxxx[i] = fact * g0_xxx[i];
        
        gx_xxxx[i] = g0_xxx[i] + fact * gx_xxx[i];

        gy_xxxx[i] = fact * gy_xxx[i];

        gz_xxxx[i] = fact * gz_xxx[i];
        
        // xxxy component
        
        fact = ry[i];
        
        g0_xxxy[i] = fact * g0_xxx[i];
        
        gx_xxxy[i] = fact * gx_xxx[i];
        
        gy_xxxy[i] = g0_xxx[i] + fact * gy_xxx[i];
        
        gz_xxxy[i] = fact * gz_xxx[i];
        
        // xxxz component
        
        fact = rz[i];
        
        g0_xxxz[i] = fact * g0_xxx[i];
        
        gx_xxxz[i] = fact * gx_xxx[i];
        
        gy_xxxz[i] = fact * gy_xxx[i];
        
        gz_xxxz[i] = g0_xxx[i] + fact * gz_xxx[i];
        
        // xxyy component
        
        fact = rx[i];
        
        g0_xxyy[i] = fact * g0_xyy[i];
        
        gx_xxyy[i] = g0_xyy[i] + fact * gx_xyy[i];
        
        gy_xxyy[i] = fact * gy_xyy[i];
        
        gz_xxyy[i] = fact * gz_xyy[i];
        
        // xxyz component
        
        fact = rx[i];
        
        g0_xxyz[i] = fact * g0_xyz[i];
        
        gx_xxyz[i] = g0_xyz[i] + fact * gx_xyz[i];
        
        gy_xxyz[i] = fact * gy_xyz[i];
        
        gz_xxyz[i] = fact * gz_xyz[i];
        
        // xxzz component
        
        fact = rx[i];
        
        g0_xxzz[i] = fact * g0_xzz[i];
        
        gx_xxzz[i] = g0_xzz[i] + fact * gx_xzz[i];
        
        gy_xxzz[i] = fact * gy_xzz[i];
        
        gz_xxzz[i] = fact * gz_xzz[i];
        
        // xyyy component
        
        fact = rx[i];
        
        g0_xyyy[i] = fact * g0_yyy[i];
        
        gx_xyyy[i] = g0_yyy[i] + fact * gx_yyy[i];
        
        gy_xyyy[i] = fact * gy_yyy[i];
        
        gz_xyyy[i] = fact * gz_yyy[i];
        
        // xyyz component
        
        fact = rx[i];
        
        g0_xyyz[i] = fact * g0_yyz[i];
        
        gx_xyyz[i] = g0_yyz[i] + fact * gx_yyz[i];
        
        gy_xyyz[i] = fact * gy_yyz[i];
        
        gz_xyyz[i] = fact * gz_yyz[i];
        
        // xyzz component
        
        fact = rx[i];
        
        g0_xyzz[i] = fact * g0_yzz[i];
        
        gx_xyzz[i] = g0_yzz[i] + fact * gx_yzz[i];
        
        gy_xyzz[i] = fact * gy_yzz[i];
        
        gz_xyzz[i] = fact * gz_yzz[i];
        
        // xzzz component
        
        fact = rx[i];
        
        g0_xzzz[i] = fact * g0_zzz[i];
        
        gx_xzzz[i] = g0_zzz[i] + fact * gx_zzz[i];
        
        gy_xzzz[i] = fact * gy_zzz[i];
        
        gz_xzzz[i] = fact * gz_zzz[i];
        
        // yyyy component
        
        fact = ry[i];
        
        g0_yyyy[i] = fact * g0_yyy[i];
        
        gx_yyyy[i] = fact * gx_yyy[i];
        
        gy_yyyy[i] = g0_yyy[i] + fact * gy_yyy[i];
        
        gz_yyyy[i] = fact * gz_yyy[i];
        
        // yyyz component
        
        fact = rz[i];
        
        g0_yyyz[i] = fact * g0_yyy[i];
        
        gx_yyyz[i] = fact * gx_yyy[i];
        
        gy_yyyz[i] = fact * gy_yyy[i];
        
        gz_yyyz[i] = g0_yyy[i] + fact * gz_yyy[i];
        
        // yyzz component
        
        fact = rz[i];
        
        g0_yyzz[i] = fact * g0_yyz[i];
        
        gx_yyzz[i] = fact * gx_yyz[i];
        
        gy_yyzz[i] = fact * gy_yyz[i];
        
        gz_yyzz[i] = g0_yyy[i] + fact * gz_yyz[i];
        
        // yzzz component
        
        fact = ry[i];
        
        g0_yzzz[i] = fact * g0_zzz[i];
        
        gx_yzzz[i] = fact * gx_zzz[i];
        
        gy_yzzz[i] = g0_zzz[i] + fact * gy_zzz[i];
        
        gz_yzzz[i] = fact * gz_zzz[i];
        
        // zzzz component
        
        fact = rz[i];
        
        g0_zzzz[i] = fact * g0_zzz[i];
        
        gx_zzzz[i] = fact * gx_zzz[i];
        
        gy_zzzz[i] = fact * gy_zzz[i];
        
        gz_zzzz[i] = g0_zzz[i] + fact * gz_zzz[i];
    }
}
    
void
compGtoTypeSForMGGA(      CMemBlock2D<double>&  gtoValues,
                    const CMemBlock2D<double>&  distances,
                    const CGtoContainer&        gtoContainer,
                    const CMemBlock2D<int32_t>& redDimensions,
                    const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitives data
    
    auto pexp = gtoContainer.getExponents(iGtoBlock);
    
    auto pnorm = gtoContainer.getNormFactors(iGtoBlock);
    
    // primitive s-type GTOs
    
    auto g0_0 = gtoValues.data(0);
    
    auto gx_0 = gtoValues.data(1);
    
    auto gy_0 = gtoValues.data(2);
    
    auto gz_0 = gtoValues.data(3);
    
    auto g2_0 = gtoValues.data(4);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_0, gx_0, gy_0, gz_0, g2_0, pexp, pnorm, rx, ry,\
                             rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        g0_0[i] = pnorm[i] * std::exp(-pexp[i] * (rx[i] * rx[i] + ry[i] * ry[i]
                                     
                                      + rz[i] * rz[i]));
        
        double fact = -2.0 * pexp[i] * g0_0[i];
        
        gx_0[i] = fact * rx[i];
        
        gy_0[i] = fact * ry[i];
        
        gz_0[i] = fact * rz[i];
        
        g2_0[i] = 3.0 * fact - 2.0 * pexp[i] * (rx[i] * gx_0[i] + ry[i] * gy_0[i]
                                            
                                                + rz[i] * gz_0[i]);
    }
}

void
compGtoTypePForMGGA(      CMemBlock2D<double>&  gtoValues,
                    const CMemBlock2D<double>&  distances,
                    const CMemBlock2D<int32_t>& redDimensions,
                    const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive s-type GTOs
    
    auto g0_0 = gtoValues.data(0);
    
    auto gx_0 = gtoValues.data(1);
    
    auto gy_0 = gtoValues.data(2);
    
    auto gz_0 = gtoValues.data(3);
    
    auto g2_0 = gtoValues.data(4);
    
    // primitive p-type GTOs
    
    auto g0_x = gtoValues.data(5);
    
    auto gx_x = gtoValues.data(6);
    
    auto gy_x = gtoValues.data(7);
    
    auto gz_x = gtoValues.data(8);
    
    auto g2_x = gtoValues.data(9);
    
    auto g0_y = gtoValues.data(10);
    
    auto gx_y = gtoValues.data(11);
    
    auto gy_y = gtoValues.data(12);
    
    auto gz_y = gtoValues.data(13);
    
    auto g2_y = gtoValues.data(14);
    
    auto g0_z = gtoValues.data(15);
    
    auto gx_z = gtoValues.data(16);
    
    auto gy_z = gtoValues.data(17);
    
    auto gz_z = gtoValues.data(18);
    
    auto g2_z = gtoValues.data(19);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_0, gx_0, gy_0, gz_0, g2_0, g0_x, gx_x, gy_x,\
                             gz_x, g2_x, g0_y, gx_y, gy_y, gz_y, g2_y, g0_z,\
                             gx_z, gy_z, gz_z, g2_z, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // x component
        
        double fact = rx[i];
        
        g0_x[i] = fact * g0_0[i];
        
        gx_x[i] = g0_0[i] + fact * gx_0[i];
        
        gy_x[i] = fact * gy_0[i];
        
        gz_x[i] = fact * gz_0[i];
        
        g2_x[i] = 2.0 * gx_0[i] + fact * g2_0[i];
        
        // y component
        
        fact = ry[i];
        
        g0_y[i] = fact * g0_0[i];
        
        gx_y[i] = fact * gx_0[i];
        
        gy_y[i] = g0_0[i] + fact * gy_0[i];
        
        gz_y[i] = fact * gz_0[i];
        
        g2_y[i] = 2.0 * gy_0[i] + fact * g2_0[i];
        
        // z component
        
        fact = rz[i];
        
        g0_z[i] = fact * g0_0[i];
        
        gx_z[i] = fact * gx_0[i];
        
        gy_z[i] = fact * gy_0[i];
        
        gz_z[i] = g0_0[i] + fact * gz_0[i];
        
        g2_z[i] = 2.0 * gz_0[i] + fact * g2_0[i];
    }
}

void
compGtoTypeDForMGGA(      CMemBlock2D<double>&  gtoValues,
                    const CMemBlock2D<double>&  distances,
                    const CMemBlock2D<int32_t>& redDimensions,
                    const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive p-type GTOs
    
    auto g0_x = gtoValues.data(5);
    
    auto gx_x = gtoValues.data(6);
    
    auto gy_x = gtoValues.data(7);
    
    auto gz_x = gtoValues.data(8);
    
    auto g2_x = gtoValues.data(9);
    
    auto g0_y = gtoValues.data(10);
    
    auto gx_y = gtoValues.data(11);
    
    auto gy_y = gtoValues.data(12);
    
    auto gz_y = gtoValues.data(13);
    
    auto g2_y = gtoValues.data(14);
    
    auto g0_z = gtoValues.data(15);
    
    auto gx_z = gtoValues.data(16);
    
    auto gy_z = gtoValues.data(17);
    
    auto gz_z = gtoValues.data(18);
    
    auto g2_z = gtoValues.data(19);
    
    // primitive d-type GTOs
    
    auto g0_xx = gtoValues.data(20);
    
    auto gx_xx = gtoValues.data(21);
    
    auto gy_xx = gtoValues.data(22);
    
    auto gz_xx = gtoValues.data(23);
    
    auto g2_xx = gtoValues.data(24);
    
    auto g0_xy = gtoValues.data(25);
    
    auto gx_xy = gtoValues.data(26);
    
    auto gy_xy = gtoValues.data(27);
    
    auto gz_xy = gtoValues.data(28);
    
    auto g2_xy = gtoValues.data(29);
    
    auto g0_xz = gtoValues.data(30);
    
    auto gx_xz = gtoValues.data(31);
    
    auto gy_xz = gtoValues.data(32);
    
    auto gz_xz = gtoValues.data(33);
    
    auto g2_xz = gtoValues.data(34);
    
    auto g0_yy = gtoValues.data(35);
    
    auto gx_yy = gtoValues.data(36);
    
    auto gy_yy = gtoValues.data(37);
    
    auto gz_yy = gtoValues.data(38);
    
    auto g2_yy = gtoValues.data(39);
    
    auto g0_yz = gtoValues.data(40);
    
    auto gx_yz = gtoValues.data(41);
    
    auto gy_yz = gtoValues.data(42);
    
    auto gz_yz = gtoValues.data(43);
    
    auto g2_yz = gtoValues.data(44);
   
    auto g0_zz = gtoValues.data(45);
    
    auto gx_zz = gtoValues.data(46);
    
    auto gy_zz = gtoValues.data(47);
    
    auto gz_zz = gtoValues.data(48);
    
    auto g2_zz = gtoValues.data(49);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_x, gx_x, gy_x, gz_x, g2_x, g0_y, gx_y, gy_y,\
                             gz_y, g2_y, g0_z, gx_z, gy_z, gz_z, g2_z, g0_xx,\
                             gx_xx, gy_xx, gz_xx, g2_xx, g0_xy, gx_xy, gy_xy,\
                             gz_xy, g2_xy, g0_xz, gx_xz, gy_xz, gz_xz, g2_xz,\
                             g0_yy, gx_yy, gy_yy, gz_yy, g2_yy, g0_yz, gx_yz,\
                             gy_yz, gz_yz, g2_yz, g0_zz, gx_zz, gy_zz, gz_zz,\
                             g2_zz, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xx component
        
        double fact = rx[i];

        g0_xx[i] = fact * g0_x[i];
        
        gx_xx[i] = g0_x[i] + fact * gx_x[i];

        gy_xx[i] = fact * gy_x[i];

        gz_xx[i] = fact * gz_x[i];
        
        g2_xx[i] = 2.0 * gx_x[i] + fact * g2_x[i];
        
        // xy index component
    
        fact = ry[i];
        
        g0_xy[i] = fact * g0_x[i];
        
        gx_xy[i] = fact * gx_x[i];
        
        gy_xy[i] = g0_x[i] + fact * gy_x[i];
        
        gz_xy[i] = fact * gz_x[i];
        
        g2_xy[i] = 2.0 * gy_x[i] + fact * g2_x[i];
        
        // xz index component
        
        fact = rz[i];
        
        g0_xz[i] = fact * g0_x[i];
        
        gx_xz[i] = fact * gx_x[i];
        
        gy_xz[i] = fact * gy_x[i];
        
        gz_xz[i] = g0_x[i] + fact * gz_x[i];
        
        g2_xz[i] = 2.0 * gz_x[i] + fact * g2_x[i];
        
        // yy component
        
        fact = ry[i];
        
        g0_yy[i] = fact * g0_y[i];
        
        gx_yy[i] = fact * gx_y[i];
        
        gy_yy[i] = g0_y[i] + fact * gy_y[i];
        
        gz_yy[i] = fact * gz_y[i];
        
        g2_yy[i] = 2.0 * gy_y[i] + fact * g2_y[i];
        
        // yz component
        
        fact = rz[i];
        
        g0_yz[i] = fact * g0_y[i];
        
        gx_yz[i] = fact * gx_y[i];
        
        gy_yz[i] = fact * gy_y[i];
        
        gz_yz[i] = g0_y[i] + fact * gz_y[i];
        
        g2_yz[i] = 2.0 * gz_y[i] + fact * g2_y[i];
        
        // zz component
        
        fact = rz[i];
        
        g0_zz[i] = fact * g0_z[i];
        
        gx_zz[i] = fact * gx_z[i];
        
        gy_zz[i] = fact * gy_z[i];
        
        gz_zz[i] = g0_z[i] + fact * gz_z[i];
        
        g2_zz[i] = 2.0 * gz_z[i] + fact * g2_z[i];
    }
}
    
void
compGtoTypeFForMGGA(      CMemBlock2D<double>&  gtoValues,
                   const CMemBlock2D<double>&  distances,
                   const CMemBlock2D<int32_t>& redDimensions,
                   const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive d-type GTOs
    
    auto g0_xx = gtoValues.data(20);
    
    auto gx_xx = gtoValues.data(21);
    
    auto gy_xx = gtoValues.data(22);
    
    auto gz_xx = gtoValues.data(23);
    
    auto g2_xx = gtoValues.data(24);
    
    auto g0_yy = gtoValues.data(35);
    
    auto gx_yy = gtoValues.data(36);
    
    auto gy_yy = gtoValues.data(37);
    
    auto gz_yy = gtoValues.data(38);
    
    auto g2_yy = gtoValues.data(39);
    
    auto g0_yz = gtoValues.data(40);
    
    auto gx_yz = gtoValues.data(41);
    
    auto gy_yz = gtoValues.data(42);
    
    auto gz_yz = gtoValues.data(43);
    
    auto g2_yz = gtoValues.data(44);
    
    auto g0_zz = gtoValues.data(45);
    
    auto gx_zz = gtoValues.data(46);
    
    auto gy_zz = gtoValues.data(47);
    
    auto gz_zz = gtoValues.data(48);
    
    auto g2_zz = gtoValues.data(49);
    
    // primitive f-type GTOs
    
    auto g0_xxx = gtoValues.data(50);
    
    auto gx_xxx = gtoValues.data(51);
    
    auto gy_xxx = gtoValues.data(52);
    
    auto gz_xxx = gtoValues.data(53);
    
    auto g2_xxx = gtoValues.data(54);
    
    auto g0_xxy = gtoValues.data(55);
    
    auto gx_xxy = gtoValues.data(56);
    
    auto gy_xxy = gtoValues.data(57);
    
    auto gz_xxy = gtoValues.data(58);
    
    auto g2_xxy = gtoValues.data(59);
    
    auto g0_xxz = gtoValues.data(60);
    
    auto gx_xxz = gtoValues.data(61);
    
    auto gy_xxz = gtoValues.data(62);
    
    auto gz_xxz = gtoValues.data(63);
    
    auto g2_xxz = gtoValues.data(64);
    
    auto g0_xyy = gtoValues.data(65);
    
    auto gx_xyy = gtoValues.data(66);
    
    auto gy_xyy = gtoValues.data(67);
    
    auto gz_xyy = gtoValues.data(68);
    
    auto g2_xyy = gtoValues.data(69);
    
    auto g0_xyz = gtoValues.data(70);
    
    auto gx_xyz = gtoValues.data(71);
    
    auto gy_xyz = gtoValues.data(72);
    
    auto gz_xyz = gtoValues.data(73);
    
    auto g2_xyz = gtoValues.data(74);
    
    auto g0_xzz = gtoValues.data(75);
    
    auto gx_xzz = gtoValues.data(76);
    
    auto gy_xzz = gtoValues.data(77);
    
    auto gz_xzz = gtoValues.data(78);
    
    auto g2_xzz = gtoValues.data(79);
    
    auto g0_yyy = gtoValues.data(80);
    
    auto gx_yyy = gtoValues.data(81);
    
    auto gy_yyy = gtoValues.data(82);
    
    auto gz_yyy = gtoValues.data(83);
    
    auto g2_yyy = gtoValues.data(84);
    
    auto g0_yyz = gtoValues.data(85);
    
    auto gx_yyz = gtoValues.data(86);
    
    auto gy_yyz = gtoValues.data(87);
    
    auto gz_yyz = gtoValues.data(88);
    
    auto g2_yyz = gtoValues.data(89);
    
    auto g0_yzz = gtoValues.data(90);
    
    auto gx_yzz = gtoValues.data(91);
    
    auto gy_yzz = gtoValues.data(92);
    
    auto gz_yzz = gtoValues.data(93);
    
    auto g2_yzz = gtoValues.data(94);
    
    auto g0_zzz = gtoValues.data(95);
    
    auto gx_zzz = gtoValues.data(96);
    
    auto gy_zzz = gtoValues.data(97);
    
    auto gz_zzz = gtoValues.data(98);
    
    auto g2_zzz = gtoValues.data(99);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_xx, gx_xx, gy_xx, gz_xx, g2_xx, g0_yy, gx_yy,\
                             gy_yy, gz_yy, g2_yy, g0_yz, gx_yz, gy_yz, gz_yz,\
                             g2_yz, g0_zz, gx_zz, gy_zz, gz_zz, g2_zz, g0_xxx,\
                             gx_xxx, gy_xxx, gz_xxx, g2_xxx, g0_xxy, gx_xxy,\
                             gy_xxy, gz_xxy, g2_xxy, g0_xxz, gx_xxz, gy_xxz,\
                             gz_xxz, g2_xxz, g0_xyy, gx_xyy, gy_xyy, gz_xyy,\
                             g2_xyy, g0_xyz, gx_xyz, gy_xyz, gz_xyz, g2_xyz,\
                             g0_xzz, gx_xzz, gy_xzz, gz_xzz, g2_xzz, g0_yyy,\
                             gx_yyy, gy_yyy, gz_yyy, g2_yyy, g0_yyz, gx_yyz,\
                             gy_yyz, gz_yyz, g2_yyz, g0_yzz, gx_yzz, gy_yzz,\
                             gz_yzz, g2_yzz, g0_zzz, gx_zzz, gy_zzz, gz_zzz,\
                             g2_zzz, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xxx component
        
        double fact = rx[i];

        g0_xxx[i] = fact * g0_xx[i];
        
        gx_xxx[i] = g0_xx[i] + fact * gx_xx[i];

        gy_xxx[i] = fact * gy_xx[i];

        gz_xxx[i] = fact * gz_xx[i];
        
        g2_xxx[i] = 2.0 * gx_xx[i] + fact * g2_xx[i];
        
        // xxy index component
    
        fact = ry[i];
        
        g0_xxy[i] = fact * g0_xx[i];
        
        gx_xxy[i] = fact * gx_xx[i];
        
        gy_xxy[i] = g0_xx[i] + fact * gy_xx[i];
        
        gz_xxy[i] = fact * gz_xx[i];
        
        g2_xxy[i] = 2.0 * gy_xx[i] + fact * g2_xx[i];
        
        // xxz index component
        
        fact = rz[i];
        
        g0_xxz[i] = fact * g0_xx[i];
        
        gx_xxz[i] = fact * gx_xx[i];
        
        gy_xxz[i] = fact * gy_xx[i];
        
        gz_xxz[i] = g0_xx[i] + fact * gz_xx[i];
        
        g2_xxz[i] = 2.0 * gz_xx[i] + fact * g2_xx[i];
        
        // xyy component
        
        fact = rx[i];
        
        g0_xyy[i] = fact * g0_yy[i];
        
        gx_xyy[i] = g0_yy[i] + fact * gx_yy[i];
        
        gy_xyy[i] = fact * gy_yy[i];
        
        gz_xyy[i] = fact * gz_yy[i];
        
        g2_xyy[i] = 2.0 * gx_yy[i] + fact * g2_yy[i];
        
        // xyz component
        
        fact = rx[i];
        
        g0_xyz[i] = fact * g0_yz[i];
        
        gx_xyz[i] = g0_yz[i] + fact * gx_yz[i];
        
        gy_xyz[i] = fact * gy_yz[i];
        
        gz_xyz[i] = fact * gz_yz[i];
        
        g2_xyz[i] = 2.0 * gx_yz[i] + fact * g2_yz[i];
        
        // xzz component
        
        fact = rx[i];
        
        g0_xzz[i] = fact * g0_zz[i];
        
        gx_xzz[i] = g0_zz[i] + fact * gx_zz[i];
        
        gy_xzz[i] = fact * gy_zz[i];
        
        gz_xzz[i] = fact * gz_zz[i];
        
        g2_xzz[i] = 2.0 * gx_zz[i] + fact * g2_zz[i];
        
        // yyy component
        
        fact = ry[i];
        
        g0_yyy[i] = fact * g0_yy[i];
        
        gx_yyy[i] = fact * gx_yy[i];
        
        gy_yyy[i] = g0_yy[i] + fact * gy_yy[i];
        
        gz_yyy[i] = fact * gz_yy[i];
        
        g2_yyy[i] = 2.0 * gy_yy[i] + fact * g2_yy[i];
        
        // yyz component
        
        fact = rz[i];
        
        g0_yyz[i] = fact * g0_yy[i];
        
        gx_yyz[i] = fact * gx_yy[i];
        
        gy_yyz[i] = fact * gy_yy[i];
        
        gz_yyz[i] = g0_yy[i] + fact * gz_yy[i];
        
        g2_yyz[i] = 2.0 * gz_yy[i] + fact * g2_yy[i];
        
        // yzz component
        
        fact = ry[i];
        
        g0_yzz[i] = fact * g0_zz[i];
        
        gx_yzz[i] = fact * gx_zz[i];
        
        gy_yzz[i] = g0_zz[i] + fact * gy_zz[i];
        
        gz_yzz[i] = fact * gz_zz[i];
        
        g2_yzz[i] = 2.0 * gy_zz[i] + fact * g2_zz[i];
        
        // zzz component
        
        fact = rz[i];
        
        g0_zzz[i] = fact * g0_zz[i];
        
        gx_zzz[i] = fact * gx_zz[i];
        
        gy_zzz[i] = fact * gy_zz[i];
        
        gz_zzz[i] = g0_zz[i] + fact * gz_zz[i];
        
        g2_zzz[i] = 2.0 * gz_zz[i] + fact * g2_zz[i];
    }
}

void
compGtoTypeGForMGGA(      CMemBlock2D<double>&  gtoValues,
                    const CMemBlock2D<double>&  distances,
                    const CMemBlock2D<int32_t>& redDimensions,
                    const int32_t               iGtoBlock)
{
    // set up number of primitive GTOs
    
    auto reddim = redDimensions.data(0);
    
    auto pnum = reddim[iGtoBlock];
    
    // set up distances
    
    auto rx = distances.data(0);
    
    auto ry = distances.data(1);
    
    auto rz = distances.data(2);
    
    // primitive f-type GTOs
    
    auto g0_xxx = gtoValues.data(50);
    
    auto gx_xxx = gtoValues.data(51);
    
    auto gy_xxx = gtoValues.data(52);
    
    auto gz_xxx = gtoValues.data(53);
    
    auto g2_xxx = gtoValues.data(54);
    
    auto g0_xyy = gtoValues.data(65);
    
    auto gx_xyy = gtoValues.data(66);
    
    auto gy_xyy = gtoValues.data(67);
    
    auto gz_xyy = gtoValues.data(68);
    
    auto g2_xyy = gtoValues.data(69);
    
    auto g0_xyz = gtoValues.data(70);
    
    auto gx_xyz = gtoValues.data(71);
    
    auto gy_xyz = gtoValues.data(72);
    
    auto gz_xyz = gtoValues.data(73);
    
    auto g2_xyz = gtoValues.data(74);
    
    auto g0_xzz = gtoValues.data(75);
    
    auto gx_xzz = gtoValues.data(76);
    
    auto gy_xzz = gtoValues.data(77);
    
    auto gz_xzz = gtoValues.data(78);
    
    auto g2_xzz = gtoValues.data(79);
    
    auto g0_yyy = gtoValues.data(80);
    
    auto gx_yyy = gtoValues.data(81);
    
    auto gy_yyy = gtoValues.data(82);
    
    auto gz_yyy = gtoValues.data(83);
    
    auto g2_yyy = gtoValues.data(84);
    
    auto g0_yyz = gtoValues.data(85);
    
    auto gx_yyz = gtoValues.data(86);
    
    auto gy_yyz = gtoValues.data(87);
    
    auto gz_yyz = gtoValues.data(88);
    
    auto g2_yyz = gtoValues.data(89);
    
    auto g0_yzz = gtoValues.data(90);
    
    auto gx_yzz = gtoValues.data(91);
    
    auto gy_yzz = gtoValues.data(92);
    
    auto gz_yzz = gtoValues.data(93);
    
    auto g2_yzz = gtoValues.data(94);
    
    auto g0_zzz = gtoValues.data(95);
    
    auto gx_zzz = gtoValues.data(96);
    
    auto gy_zzz = gtoValues.data(97);
    
    auto gz_zzz = gtoValues.data(98);
    
    auto g2_zzz = gtoValues.data(99);
    
    // primitive g-type GTOs
    
    auto g0_xxxx = gtoValues.data(100);
    
    auto gx_xxxx = gtoValues.data(101);
    
    auto gy_xxxx = gtoValues.data(102);
    
    auto gz_xxxx = gtoValues.data(103);
    
    auto g2_xxxx = gtoValues.data(104);
    
    auto g0_xxxy = gtoValues.data(105);
    
    auto gx_xxxy = gtoValues.data(106);
    
    auto gy_xxxy = gtoValues.data(107);
    
    auto gz_xxxy = gtoValues.data(108);
    
    auto g2_xxxy = gtoValues.data(109);
    
    auto g0_xxxz = gtoValues.data(110);
    
    auto gx_xxxz = gtoValues.data(111);
    
    auto gy_xxxz = gtoValues.data(112);
    
    auto gz_xxxz = gtoValues.data(113);
    
    auto g2_xxxz = gtoValues.data(114);
    
    auto g0_xxyy = gtoValues.data(115);
    
    auto gx_xxyy = gtoValues.data(116);
    
    auto gy_xxyy = gtoValues.data(117);
    
    auto gz_xxyy = gtoValues.data(118);
    
    auto g2_xxyy = gtoValues.data(119);
    
    auto g0_xxyz = gtoValues.data(120);
    
    auto gx_xxyz = gtoValues.data(121);
    
    auto gy_xxyz = gtoValues.data(122);
    
    auto gz_xxyz = gtoValues.data(123);
    
    auto g2_xxyz = gtoValues.data(124);
    
    auto g0_xxzz = gtoValues.data(125);
    
    auto gx_xxzz = gtoValues.data(126);
    
    auto gy_xxzz = gtoValues.data(127);
    
    auto gz_xxzz = gtoValues.data(128);
    
    auto g2_xxzz = gtoValues.data(129);
    
    auto g0_xyyy = gtoValues.data(130);
    
    auto gx_xyyy = gtoValues.data(131);
    
    auto gy_xyyy = gtoValues.data(132);
    
    auto gz_xyyy = gtoValues.data(133);
    
    auto g2_xyyy = gtoValues.data(134);
    
    auto g0_xyyz = gtoValues.data(135);
    
    auto gx_xyyz = gtoValues.data(136);
    
    auto gy_xyyz = gtoValues.data(137);
    
    auto gz_xyyz = gtoValues.data(138);
    
    auto g2_xyyz = gtoValues.data(139);
    
    auto g0_xyzz = gtoValues.data(140);
    
    auto gx_xyzz = gtoValues.data(141);
    
    auto gy_xyzz = gtoValues.data(142);
    
    auto gz_xyzz = gtoValues.data(143);
    
    auto g2_xyzz = gtoValues.data(144);
    
    auto g0_xzzz = gtoValues.data(145);
    
    auto gx_xzzz = gtoValues.data(146);
    
    auto gy_xzzz = gtoValues.data(147);
    
    auto gz_xzzz = gtoValues.data(148);
    
    auto g2_xzzz = gtoValues.data(149);
    
    auto g0_yyyy = gtoValues.data(150);
    
    auto gx_yyyy = gtoValues.data(151);
    
    auto gy_yyyy = gtoValues.data(152);
    
    auto gz_yyyy = gtoValues.data(153);
    
    auto g2_yyyy = gtoValues.data(154);
    
    auto g0_yyyz = gtoValues.data(155);
    
    auto gx_yyyz = gtoValues.data(156);
    
    auto gy_yyyz = gtoValues.data(157);
    
    auto gz_yyyz = gtoValues.data(158);
    
    auto g2_yyyz = gtoValues.data(159);
    
    auto g0_yyzz = gtoValues.data(160);
    
    auto gx_yyzz = gtoValues.data(161);
    
    auto gy_yyzz = gtoValues.data(162);
    
    auto gz_yyzz = gtoValues.data(163);
    
    auto g2_yyzz = gtoValues.data(164);
    
    auto g0_yzzz = gtoValues.data(165);
    
    auto gx_yzzz = gtoValues.data(166);
    
    auto gy_yzzz = gtoValues.data(167);
    
    auto gz_yzzz = gtoValues.data(168);
    
    auto g2_yzzz = gtoValues.data(169);
    
    auto g0_zzzz = gtoValues.data(170);
    
    auto gx_zzzz = gtoValues.data(171);
    
    auto gy_zzzz = gtoValues.data(172);
    
    auto gz_zzzz = gtoValues.data(173);
    
    auto g2_zzzz = gtoValues.data(174);
    
    // loop over primitive GTOs
    
    #pragma omp simd aligned(g0_xxx, gx_xxx, gy_xxx, gz_xxx, g2_xxx, g0_xyy,\
                             gx_xyy, gy_xyy, gz_xyy, g2_xyy, g0_xyz, gx_xyz,\
                             gy_xyz, gz_xyz, g2_xyz, g0_xzz, gx_xzz, gy_xzz,\
                             gz_xzz, g2_xzz, g0_yyy, gx_yyy, gy_yyy, gz_yyy,\
                             g2_yyy, g0_yyz, gx_yyz, gy_yyz, gz_yyz, g2_yyz,\
                             g0_yzz, gx_yzz, gy_yzz, gz_yzz, g2_yzz, g0_zzz,\
                             gx_zzz, gy_zzz, gz_zzz, g2_zzz, g0_xxxx, gx_xxxx,\
                             gy_xxxx, gz_xxxx, g2_xxxx, g0_xxxy, gx_xxxy,\
                             gy_xxxy, gz_xxxy, g2_xxxy, g0_xxxz, gx_xxxz,\
                             gy_xxxz, gz_xxxz, g2_xxxz, g0_xxyy, gx_xxyy,\
                             gy_xxyy, gz_xxyy, g2_xxyy, g0_xxyz, gx_xxyz,\
                             gy_xxyz, gz_xxyz, g2_xxyz, g0_xxzz, gx_xxzz,\
                             gy_xxzz, gz_xxzz, g2_xxzz, g0_xyyy, gx_xyyy,\
                             gy_xyyy, gz_xyyy, g2_xyyy, g0_xyyz, gx_xyyz,\
                             gy_xyyz, gz_xyyz, g2_xyyz, g0_xyzz, gx_xyzz,\
                             gy_xyzz, gz_xyzz, g2_xyzz, g0_xzzz, gx_xzzz,\
                             gy_xzzz, gz_xzzz, g2_xzzz, g0_yyyy, gx_yyyy,\
                             gy_yyyy, gz_yyyy, g2_yyyy, g0_yyyz, gx_yyyz,\
                             gy_yyyz, gz_yyyz, g2_yyyz, g0_yyzz, gx_yyzz,\
                             gy_yyzz, gz_yyzz, g2_yyzz, g0_yzzz, gx_yzzz,\
                             gy_yzzz, gz_yzzz, g2_yzzz, g0_zzzz, gx_zzzz,\
                             gy_zzzz, gz_zzzz, g2_zzzz, rx, ry, rz:VLX_ALIGN)
    for (int32_t i = 0; i < pnum; i++)
    {
        // xxxx component
        
        double fact = rx[i];

        g0_xxxx[i] = fact * g0_xxx[i];
        
        gx_xxxx[i] = g0_xxx[i] + fact * gx_xxx[i];

        gy_xxxx[i] = fact * gy_xxx[i];

        gz_xxxx[i] = fact * gz_xxx[i];
        
        g2_xxxx[i] = 2.0 * gx_xxx[i] + fact * g2_xxx[i];
        
        // xxxy component
        
        fact = ry[i];
        
        g0_xxxy[i] = fact * g0_xxx[i];
        
        gx_xxxy[i] = fact * gx_xxx[i];
        
        gy_xxxy[i] = g0_xxx[i] + fact * gy_xxx[i];
        
        gz_xxxy[i] = fact * gz_xxx[i];
        
        g2_xxxy[i] = 2.0 * gy_xxx[i] + fact * g2_xxx[i];
        
        // xxxz component
        
        fact = rz[i];
        
        g0_xxxz[i] = fact * g0_xxx[i];
        
        gx_xxxz[i] = fact * gx_xxx[i];
        
        gy_xxxz[i] = fact * gy_xxx[i];
        
        gz_xxxz[i] = g0_xxx[i] + fact * gz_xxx[i];
        
        g2_xxxz[i] = 2.0 * gz_xxx[i] + fact * g2_xxx[i];
        
        // xxyy component
        
        fact = rx[i];
        
        g0_xxyy[i] = fact * g0_xyy[i];
        
        gx_xxyy[i] = g0_xyy[i] + fact * gx_xyy[i];
        
        gy_xxyy[i] = fact * gy_xyy[i];
        
        gz_xxyy[i] = fact * gz_xyy[i];
        
        g2_xxyy[i] = 2.0 * gx_xyy[i] + fact * g2_xyy[i];
        
        // xxyz component
        
        fact = rx[i];
        
        g0_xxyz[i] = fact * g0_xyz[i];
        
        gx_xxyz[i] = g0_xyz[i] + fact * gx_xyz[i];
        
        gy_xxyz[i] = fact * gy_xyz[i];
        
        gz_xxyz[i] = fact * gz_xyz[i];
        
        g2_xxyz[i] = 2.0 * gx_xyz[i] + fact * g2_xyz[i];
        
        // xxzz component
        
        fact = rx[i];
        
        g0_xxzz[i] = fact * g0_xzz[i];
        
        gx_xxzz[i] = g0_xzz[i] + fact * gx_xzz[i];
        
        gy_xxzz[i] = fact * gy_xzz[i];
        
        gz_xxzz[i] = fact * gz_xzz[i];
        
        g2_xxzz[i] = 2.0 * gx_xzz[i] + fact * g2_xzz[i];
        
        // xyyy component
        
        fact = rx[i];
        
        g0_xyyy[i] = fact * g0_yyy[i];
        
        gx_xyyy[i] = g0_yyy[i] + fact * gx_yyy[i];
        
        gy_xyyy[i] = fact * gy_yyy[i];
        
        gz_xyyy[i] = fact * gz_yyy[i];
        
        g2_xyyy[i] = 2.0 * gx_yyy[i] + fact * g2_yyy[i];
        
        // xyyz component
        
        fact = rx[i];
        
        g0_xyyz[i] = fact * g0_yyz[i];
        
        gx_xyyz[i] = g0_yyz[i] + fact * gx_yyz[i];
        
        gy_xyyz[i] = fact * gy_yyz[i];
        
        gz_xyyz[i] = fact * gz_yyz[i];
        
        g2_xyyz[i] = 2.0 * gx_yyz[i] + fact * g2_yyz[i];
        
        // xyzz component
        
        fact = rx[i];
        
        g0_xyzz[i] = fact * g0_yzz[i];
        
        gx_xyzz[i] = g0_yzz[i] + fact * gx_yzz[i];
        
        gy_xyzz[i] = fact * gy_yzz[i];
        
        gz_xyzz[i] = fact * gz_yzz[i];
        
        g2_xyzz[i] = 2.0 * gx_yzz[i] + fact * g2_yzz[i];
        
        // xzzz component
        
        fact = rx[i];
        
        g0_xzzz[i] = fact * g0_zzz[i];
        
        gx_xzzz[i] = g0_zzz[i] + fact * gx_zzz[i];
        
        gy_xzzz[i] = fact * gy_zzz[i];
        
        gz_xzzz[i] = fact * gz_zzz[i];
        
        g2_xzzz[i] = 2.0 * gx_zzz[i] + fact * g2_zzz[i];
        
        // yyyy component
        
        fact = ry[i];
        
        g0_yyyy[i] = fact * g0_yyy[i];
        
        gx_yyyy[i] = fact * gx_yyy[i];
        
        gy_yyyy[i] = g0_yyy[i] + fact * gy_yyy[i];
        
        gz_yyyy[i] = fact * gz_yyy[i];
        
        g2_yyyy[i] = 2.0 * gy_yyy[i] + fact * g2_yyy[i];
        
        // yyyz component
        
        fact = rz[i];
        
        g0_yyyz[i] = fact * g0_yyy[i];
        
        gx_yyyz[i] = fact * gx_yyy[i];
        
        gy_yyyz[i] = fact * gy_yyy[i];
        
        gz_yyyz[i] = g0_yyy[i] + fact * gz_yyy[i];
        
        g2_yyyz[i] = 2.0 * gz_yyy[i] + fact * g2_yyy[i];
        
        // yyzz component
        
        fact = rz[i];
        
        g0_yyzz[i] = fact * g0_yyz[i];
        
        gx_yyzz[i] = fact * gx_yyz[i];
        
        gy_yyzz[i] = fact * gy_yyz[i];
        
        gz_yyzz[i] = g0_yyy[i] + fact * gz_yyz[i];
        
        g2_yyzz[i] = 2.0 * gy_yyz[i] + fact * g2_yyz[i];
        
        // yzzz component
        
        fact = ry[i];
        
        g0_yzzz[i] = fact * g0_zzz[i];
        
        gx_yzzz[i] = fact * gx_zzz[i];
        
        gy_yzzz[i] = g0_zzz[i] + fact * gy_zzz[i];
        
        gz_yzzz[i] = fact * gz_zzz[i];
        
        g2_yzzz[i] = 2.0 * gy_zzz[i] + fact * g2_zzz[i];
        
        // zzzz component
        
        fact = rz[i];
        
        g0_zzzz[i] = fact * g0_zzz[i];
        
        gx_zzzz[i] = fact * gx_zzz[i];
        
        gy_zzzz[i] = fact * gy_zzz[i];
        
        gz_zzzz[i] = g0_zzz[i] + fact * gz_zzz[i];
        
        g2_zzzz[i] = 2.0 * gz_zzz[i] + fact * g2_zzz[i];
    }
}
    
} // gtorec namespace
