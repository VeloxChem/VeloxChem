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

#include "GtoEvaluator.hpp"

#include <algorithm>
#include <cmath>

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

namespace gtoeval {  // gtoeval namespace

std::array<double, 6>
getGridBoxDimension(const int32_t gridBlockPosition,
                    const int32_t nGridPoints,
                    const double* xcoords,
                    const double* ycoords,
                    const double* zcoords)
{
    double xmin = xcoords[gridBlockPosition], ymin = ycoords[gridBlockPosition], zmin = zcoords[gridBlockPosition];

    double xmax = xcoords[gridBlockPosition], ymax = ycoords[gridBlockPosition], zmax = zcoords[gridBlockPosition];

    for (int32_t g = 0; g < nGridPoints; g++)
    {
        xmin = std::min(xmin, xcoords[gridBlockPosition + g]);

        ymin = std::min(ymin, ycoords[gridBlockPosition + g]);

        zmin = std::min(zmin, zcoords[gridBlockPosition + g]);

        xmax = std::max(xmax, xcoords[gridBlockPosition + g]);

        ymax = std::max(ymax, ycoords[gridBlockPosition + g]);

        zmax = std::max(zmax, zcoords[gridBlockPosition + g]);
    }

    return std::array<double, 6>({xmin, ymin, zmin, xmax, ymax, zmax});
}

void
preScreenGtos(CMemBlock<int32_t>&          skipCgtoIds,
              CMemBlock<int32_t>&          skipAOIds,
              const CGtoContainer*         gtoContainer,
              const int32_t                gtoDeriv,
              const double                 gtoThreshold,
              const std::array<double, 6>& boxDimension)
{
    skipCgtoIds.zero();

    skipAOIds.zero();

    double xmin = boxDimension[0], ymin = boxDimension[1], zmin = boxDimension[2];

    double xmax = boxDimension[3], ymax = boxDimension[4], zmax = boxDimension[5];

    for (int32_t i = 0, cgto_count = 0; i < gtoContainer->getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtoContainer->getGtoBlock(i);

        auto bang = bgtos.getAngularMomentum();

        auto bfnorms = bgtos.getNormFactors();

        auto bfexps = bgtos.getExponents();

        auto bfx = bgtos.getCoordinatesX();

        auto bfy = bgtos.getCoordinatesY();

        auto bfz = bgtos.getCoordinatesZ();

        auto spos = bgtos.getStartPositions();

        auto epos = bgtos.getEndPositions();

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++, cgto_count++)
        {
            // contracted GTO screening

            auto firstprim = spos[j];

            double rx = std::max({xmin - bfx[firstprim], bfx[firstprim] - xmax, 0.0});

            double ry = std::max({ymin - bfy[firstprim], bfy[firstprim] - ymax, 0.0});

            double rz = std::max({zmin - bfz[firstprim], bfz[firstprim] - zmax, 0.0});

            auto r2 = rx * rx + ry * ry + rz * rz;

            if (r2 > 1.0)
            {
                auto minexp = bfexps[firstprim];

                auto maxexp = bfexps[firstprim];

                auto maxcoef = std::fabs(bfnorms[firstprim]);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = std::fabs(bfnorms[iprim]);

                    minexp = std::min(minexp, bexp);

                    maxexp = std::max(maxexp, bexp);

                    maxcoef = std::max(maxcoef, bnorm);
                }

                // 0th-order derivative
                // gto:                    r^{ang}   |C| exp(-alpha r^2)

                // 1st-order derivative
                // gto_m:              ang r^{ang-1} |C| exp(-alpha r^2)
                // gto_p:           2alpha r^{ang+1} |C| exp(-alpha r^2)

                // 2nd-order derivative
                // gto_m2:     ang (ang-1) r^{ang-2} |C| exp(-alpha r^2)
                // gto   : 2alpha (2ang+1) r^{ang}   |C| exp(-alpha r^2)
                // gto_p2:        4alpha^2 r^{ang+2} |C| exp(-alpha r^2)

                auto gtolimit = maxcoef * std::exp(-minexp * r2);

                auto r = std::sqrt(r2);

                for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit *= r;

                if (gtoDeriv > 0)
                {
                    gtolimit = std::max(gtolimit, 2.0 * maxexp * r * gtolimit);

                    if (bang > 0) gtolimit = std::max(gtolimit, gtolimit / r * bang);
                }

                if (gtoDeriv > 1)
                {
                    gtolimit = std::max({gtolimit, 4.0 * maxexp * maxexp * r2 * gtolimit,

                                         2.0 * maxexp * (2 * bang + 1) * gtolimit});

                    if (bang > 1) gtolimit = std::max(gtolimit, gtolimit / r2 * bang * (bang - 1));
                }

                if (gtolimit < gtoThreshold)
                {
                    skipCgtoIds.data()[cgto_count] = 1;

                    auto bnspher = angmom::to_SphericalComponents(bang);

                    for (int32_t k = 0; k < bnspher; k++)
                    {
                        auto ao_idx = (bgtos.getIdentifiers(k))[j];

                        skipAOIds.data()[ao_idx] = 1;
                    }
                }
            }
        }
    }
}

void
computeGtosValuesForLDA(CMemBlock2D<double>&      gtoValues,
                        const CGtoContainer*      gtoContainer,
                        const double*             gridCoordinatesX,
                        const double*             gridCoordinatesY,
                        const double*             gridCoordinatesZ,
                        const int32_t             gridBlockPosition,
                        const int32_t             gridOffset,
                        const int32_t             nGridPoints,
                        const CMemBlock<int32_t>& skipCgtoIds)
{
    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);

    // loop over GTOs container data

    for (int32_t i = 0, cgto_count = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);

        // angular momentum data for bra and ket

        auto bang = bgtos.getAngularMomentum();

        // set up Cartesian GTOs buffer

        int32_t nvcomp = 1;  // xcfun::lda

        auto bncart = angmom::to_CartesianComponents(bang);

        auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();

        // set up spherical GTOs buffer

        auto bnspher = angmom::to_SphericalComponents(bang);

        CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);

        // set up pointers to primitives data

        auto bfnorms = bgtos.getNormFactors();

        auto bfexps = bgtos.getExponents();

        // set up pointers to primitives coordinates

        auto bfx = bgtos.getCoordinatesX();

        auto bfy = bgtos.getCoordinatesY();

        auto bfz = bgtos.getCoordinatesZ();

        // set up coordinates to primitives positions

        auto spos = bgtos.getStartPositions();

        auto epos = bgtos.getEndPositions();

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++, cgto_count++)
        {
            if (skipCgtoIds.data()[cgto_count]) continue;

            if (bang == 0)
            {
                // s-type GTOs on grid

                bspherbuff.zero();

                auto f0_0 = bspherbuff.data(0);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, f0_0 : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double g0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        f0_0[g] += g0;
                    }
                }
            }
            else if (bang == 1)
            {
                // p-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_x = bcartbuff.data(0);

                auto f0_y = bcartbuff.data(1);

                auto f0_z = bcartbuff.data(2);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                                             f0_x, f0_y, f0_z : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        f0_x[g] += f0_0 * dx;

                        f0_y[g] += f0_0 * dy;

                        f0_z[g] += f0_0 * dz;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(1), 0, 0, nGridPoints, 1);
            }
            else if (bang == 2)
            {
                // d-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_xx = bcartbuff.data(0);

                auto f0_xy = bcartbuff.data(1);

                auto f0_xz = bcartbuff.data(2);

                auto f0_yy = bcartbuff.data(3);

                auto f0_yz = bcartbuff.data(4);

                auto f0_zz = bcartbuff.data(5);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                                             f0_xx, f0_xy, f0_xz, f0_yy, f0_yz, f0_zz: VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double f0_x = dx * f0_0;

                        double f0_y = dy * f0_0;

                        double f0_z = dz * f0_0;

                        f0_xx[g] += f0_x * dx;

                        f0_xy[g] += f0_x * dy;

                        f0_xz[g] += f0_x * dz;

                        f0_yy[g] += f0_y * dy;

                        f0_yz[g] += f0_y * dz;

                        f0_zz[g] += f0_z * dz;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(2), 0, 0, nGridPoints, 1);
            }
            else
            {
                // FIX ME: implement l > 2 cases
            }

            // distribute j-th GTO values into grid values matrix

            for (int32_t k = 0; k < bnspher; k++)
            {
                auto bgaos = bspherbuff.data(k);

                auto idx = bgtos.getIdentifiers(k)[j];

                auto gvals = gtoValues.data(idx);

                #pragma omp simd aligned(bgaos, gvals : VLX_ALIGN)
                for (int32_t g = 0; g < nGridPoints; g++)
                {
                    gvals[gridOffset + g] = bgaos[g];
                }
            }
        }
    }
}

void
computeGtosValuesForGGA(CMemBlock2D<double>&      gtoValues,
                        CMemBlock2D<double>&      gtoValuesX,
                        CMemBlock2D<double>&      gtoValuesY,
                        CMemBlock2D<double>&      gtoValuesZ,
                        const CGtoContainer*      gtoContainer,
                        const double*             gridCoordinatesX,
                        const double*             gridCoordinatesY,
                        const double*             gridCoordinatesZ,
                        const int32_t             gridBlockPosition,
                        const int32_t             gridOffset,
                        const int32_t             nGridPoints,
                        const CMemBlock<int32_t>& skipCgtoIds)
{
    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);

    // loop over GTOs container data

    for (int32_t i = 0, cgto_count = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);

        // angular momentum data for bra and ket

        auto bang = bgtos.getAngularMomentum();

        // set up Cartesian GTOs buffer

        int32_t nvcomp = 4;  // xcfun::gga

        auto bncart = angmom::to_CartesianComponents(bang);

        auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();

        // set up spherical GTOs buffer

        auto bnspher = angmom::to_SphericalComponents(bang);

        CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);

        // set up pointers to primitives data

        auto bfnorms = bgtos.getNormFactors();

        auto bfexps = bgtos.getExponents();

        // set up pointers to primitives coordinates

        auto bfx = bgtos.getCoordinatesX();

        auto bfy = bgtos.getCoordinatesY();

        auto bfz = bgtos.getCoordinatesZ();

        // set up coordinates to primitives positions

        auto spos = bgtos.getStartPositions();

        auto epos = bgtos.getEndPositions();

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++, cgto_count++)
        {
            if (skipCgtoIds.data()[cgto_count]) continue;

            if (bang == 0)
            {
                // s-type GTOs on grid

                bspherbuff.zero();

                auto f0_0 = bspherbuff.data(0);

                auto fx_0 = bspherbuff.data(1);

                auto fy_0 = bspherbuff.data(2);

                auto fz_0 = bspherbuff.data(3);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, f0_0, fx_0, fy_0, fz_0 : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double g0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double g1 = -2.0 * bexp * g0;

                        f0_0[g] += g0;

                        fx_0[g] += dx * g1;

                        fy_0[g] += dy * g1;

                        fz_0[g] += dz * g1;
                    }
                }
            }
            else if (bang == 1)
            {
                // p-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_x = bcartbuff.data(0);

                auto fx_x = bcartbuff.data(1);

                auto fy_x = bcartbuff.data(2);

                auto fz_x = bcartbuff.data(3);

                auto f0_y = bcartbuff.data(4);

                auto fx_y = bcartbuff.data(5);

                auto fy_y = bcartbuff.data(6);

                auto fz_y = bcartbuff.data(7);

                auto f0_z = bcartbuff.data(8);

                auto fx_z = bcartbuff.data(9);

                auto fy_z = bcartbuff.data(10);

                auto fz_z = bcartbuff.data(11);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                                             f0_x, fx_x, fy_x, fz_x, \
                                             f0_y, fx_y, fy_y, fz_y, \
                                             f0_z, fx_z, fy_z, fz_z : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double fg_0 = -2.0 * bexp;

                        // leading p_x

                        f0_x[g] += f0_0 * dx;

                        fx_x[g] += f0_0 * (1.0 + dx * dx * fg_0);

                        fy_x[g] += f0_0 * dx * dy * fg_0;

                        fz_x[g] += f0_0 * dx * dz * fg_0;

                        // leading p_y

                        f0_y[g] += f0_0 * dy;

                        fx_y[g] += f0_0 * dy * dx * fg_0;

                        fy_y[g] += f0_0 * (1.0 + dy * dy * fg_0);

                        fz_y[g] += f0_0 * dy * dz * fg_0;

                        // leading p_z

                        f0_z[g] += f0_0 * dz;

                        fx_z[g] += f0_0 * dz * dx * fg_0;

                        fy_z[g] += f0_0 * dz * dy * fg_0;

                        fz_z[g] += f0_0 * (1.0 + dz * dz * fg_0);
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(1), 0, 0, nGridPoints, 4);
            }
            else if (bang == 2)
            {
                // d-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_xx = bcartbuff.data(0);

                auto fx_xx = bcartbuff.data(1);

                auto fy_xx = bcartbuff.data(2);

                auto fz_xx = bcartbuff.data(3);

                auto f0_xy = bcartbuff.data(4);

                auto fx_xy = bcartbuff.data(5);

                auto fy_xy = bcartbuff.data(6);

                auto fz_xy = bcartbuff.data(7);

                auto f0_xz = bcartbuff.data(8);

                auto fx_xz = bcartbuff.data(9);

                auto fy_xz = bcartbuff.data(10);

                auto fz_xz = bcartbuff.data(11);

                auto f0_yy = bcartbuff.data(12);

                auto fx_yy = bcartbuff.data(13);

                auto fy_yy = bcartbuff.data(14);

                auto fz_yy = bcartbuff.data(15);

                auto f0_yz = bcartbuff.data(16);

                auto fx_yz = bcartbuff.data(17);

                auto fy_yz = bcartbuff.data(18);

                auto fz_yz = bcartbuff.data(19);

                auto f0_zz = bcartbuff.data(20);

                auto fx_zz = bcartbuff.data(21);

                auto fy_zz = bcartbuff.data(22);

                auto fz_zz = bcartbuff.data(23);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                                             f0_xx, fx_xx, fy_xx, fz_xx, \
                                             f0_xy, fx_xy, fy_xy, fz_xy, \
                                             f0_xz, fx_xz, fy_xz, fz_xz, \
                                             f0_yy, fx_yy, fy_yy, fz_yy, \
                                             f0_yz, fx_yz, fy_yz, fz_yz, \
                                             f0_zz, fx_zz, fy_zz, fz_zz : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double fg_0 = -2.0 * bexp;

                        // leading xx component

                        double f0_x = dx * f0_0;

                        f0_xx[g] += f0_x * dx;

                        fx_xx[g] += f0_x * (2.0 + fg_0 * dx * dx);

                        fy_xx[g] += f0_x * fg_0 * dx * dy;

                        fz_xx[g] += f0_x * fg_0 * dx * dz;

                        // leading xy component

                        double f0_y = dy * f0_0;

                        f0_xy[g] += f0_x * dy;

                        fx_xy[g] += f0_y * (1.0 + fg_0 * dx * dx);

                        fy_xy[g] += f0_x * (1.0 + fg_0 * dy * dy);

                        fz_xy[g] += f0_x * fg_0 * dy * dz;

                        // leading xz component

                        double f0_z = dz * f0_0;

                        f0_xz[g] += f0_x * dz;

                        fx_xz[g] += f0_z * (1.0 + fg_0 * dx * dx);

                        fy_xz[g] += f0_x * fg_0 * dz * dy;

                        fz_xz[g] += f0_x * (1.0 + fg_0 * dz * dz);

                        // leading yy component

                        f0_yy[g] += f0_y * dy;

                        fx_yy[g] += f0_y * fg_0 * dy * dx;

                        fy_yy[g] += f0_y * (2.0 + fg_0 * dy * dy);

                        fz_yy[g] += f0_y * fg_0 * dy * dz;

                        // leading yz component

                        f0_yz[g] += f0_y * dz;

                        fx_yz[g] += f0_y * fg_0 * dz * dx;

                        fy_yz[g] += f0_z * (1.0 + fg_0 * dy * dy);

                        fz_yz[g] += f0_y * (1.0 + fg_0 * dz * dz);

                        // leading zz component

                        f0_zz[g] += f0_z * dz;

                        fx_zz[g] += f0_z * fg_0 * dz * dx;

                        fy_zz[g] += f0_z * fg_0 * dz * dy;

                        fz_zz[g] += f0_z * (2.0 + fg_0 * dz * dz);
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(2), 0, 0, nGridPoints, 4);
            }
            else
            {
                // FIX ME: implement l > 2 cases
            }

            // distribute j-th GTO values into grid values matrix

            for (int32_t k = 0; k < bnspher; k++)
            {
                auto bgaos = bspherbuff.data(4 * k);

                auto bgaox = bspherbuff.data(4 * k + 1);

                auto bgaoy = bspherbuff.data(4 * k + 2);

                auto bgaoz = bspherbuff.data(4 * k + 3);

                auto idx = bgtos.getIdentifiers(k)[j];

                auto gvals = gtoValues.data(idx);

                auto gvalx = gtoValuesX.data(idx);

                auto gvaly = gtoValuesY.data(idx);

                auto gvalz = gtoValuesZ.data(idx);

                #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz : VLX_ALIGN)
                for (int32_t g = 0; g < nGridPoints; g++)
                {
                    gvals[gridOffset + g] = bgaos[g];

                    gvalx[gridOffset + g] = bgaox[g];

                    gvaly[gridOffset + g] = bgaoy[g];

                    gvalz[gridOffset + g] = bgaoz[g];
                }
            }
        }
    }
}

void
computeGtosValuesForMetaGGA(CMemBlock2D<double>&      gtoValues,
                            CMemBlock2D<double>&      gtoValuesX,
                            CMemBlock2D<double>&      gtoValuesY,
                            CMemBlock2D<double>&      gtoValuesZ,
                            CMemBlock2D<double>&      gtoValuesXX,
                            CMemBlock2D<double>&      gtoValuesXY,
                            CMemBlock2D<double>&      gtoValuesXZ,
                            CMemBlock2D<double>&      gtoValuesYY,
                            CMemBlock2D<double>&      gtoValuesYZ,
                            CMemBlock2D<double>&      gtoValuesZZ,
                            const CGtoContainer*      gtoContainer,
                            const double*             gridCoordinatesX,
                            const double*             gridCoordinatesY,
                            const double*             gridCoordinatesZ,
                            const int32_t             gridBlockPosition,
                            const int32_t             gridOffset,
                            const int32_t             nGridPoints,
                            const CMemBlock<int32_t>& skipCgtoIds)
{
    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);

    // loop over GTOs container data

    for (int32_t i = 0, cgto_count = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);

        // angular momentum data for bra and ket

        auto bang = bgtos.getAngularMomentum();

        // set up Cartesian GTOs buffer

        int32_t nvcomp = 10;  // xcfun::mgga

        auto bncart = angmom::to_CartesianComponents(bang);

        auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();

        // set up spherical GTOs buffer

        auto bnspher = angmom::to_SphericalComponents(bang);

        CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);

        // set up pointers to primitives data

        auto bfnorms = bgtos.getNormFactors();

        auto bfexps = bgtos.getExponents();

        // set up pointers to primitives coordinates

        auto bfx = bgtos.getCoordinatesX();

        auto bfy = bgtos.getCoordinatesY();

        auto bfz = bgtos.getCoordinatesZ();

        // set up coordinates to primitives positions

        auto spos = bgtos.getStartPositions();

        auto epos = bgtos.getEndPositions();

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++, cgto_count++)
        {
            if (skipCgtoIds.data()[cgto_count]) continue;

            if (bang == 0)
            {
                // s-type GTOs on grid

                bspherbuff.zero();

                auto f0_0 = bspherbuff.data(0);

                auto fx_0 = bspherbuff.data(1);

                auto fy_0 = bspherbuff.data(2);

                auto fz_0 = bspherbuff.data(3);

                auto fxx_0 = bspherbuff.data(4);

                auto fxy_0 = bspherbuff.data(5);

                auto fxz_0 = bspherbuff.data(6);

                auto fyy_0 = bspherbuff.data(7);

                auto fyz_0 = bspherbuff.data(8);

                auto fzz_0 = bspherbuff.data(9);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                            f0_0, fx_0, fy_0, fz_0, fxx_0, fxy_0, fxz_0, fyy_0, fyz_0, fzz_0 : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double g0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double g1 = -2.0 * bexp * g0;

                        double g2 = -2.0 * bexp * g1;

                        f0_0[g] += g0;

                        fx_0[g] += dx * g1;

                        fy_0[g] += dy * g1;

                        fz_0[g] += dz * g1;

                        fxx_0[g] += (dx * dx) * g2 + g1;

                        fxy_0[g] += (dx * dy) * g2;

                        fxz_0[g] += (dx * dz) * g2;

                        fyy_0[g] += (dy * dy) * g2 + g1;

                        fyz_0[g] += (dy * dz) * g2;

                        fzz_0[g] += (dz * dz) * g2 + g1;
                    }
                }
            }
            else if (bang == 1)
            {
                // p-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_x = bcartbuff.data(0);

                auto fx_x = bcartbuff.data(1);

                auto fy_x = bcartbuff.data(2);

                auto fz_x = bcartbuff.data(3);

                auto fxx_x = bcartbuff.data(4);

                auto fxy_x = bcartbuff.data(5);

                auto fxz_x = bcartbuff.data(6);

                auto fyy_x = bcartbuff.data(7);

                auto fyz_x = bcartbuff.data(8);

                auto fzz_x = bcartbuff.data(9);

                auto f0_y = bcartbuff.data(10);

                auto fx_y = bcartbuff.data(11);

                auto fy_y = bcartbuff.data(12);

                auto fz_y = bcartbuff.data(13);

                auto fxx_y = bcartbuff.data(14);

                auto fxy_y = bcartbuff.data(15);

                auto fxz_y = bcartbuff.data(16);

                auto fyy_y = bcartbuff.data(17);

                auto fyz_y = bcartbuff.data(18);

                auto fzz_y = bcartbuff.data(19);

                auto f0_z = bcartbuff.data(20);

                auto fx_z = bcartbuff.data(21);

                auto fy_z = bcartbuff.data(22);

                auto fz_z = bcartbuff.data(23);

                auto fxx_z = bcartbuff.data(24);

                auto fxy_z = bcartbuff.data(25);

                auto fxz_z = bcartbuff.data(26);

                auto fyy_z = bcartbuff.data(27);

                auto fyz_z = bcartbuff.data(28);

                auto fzz_z = bcartbuff.data(29);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                            f0_x, fx_x, fy_x, fz_x, fxx_x, fxy_x, fxz_x, fyy_x, fyz_x, fzz_x, \
                            f0_y, fx_y, fy_y, fz_y, fxx_y, fxy_y, fxz_y, fyy_y, fyz_y, fzz_y, \
                            f0_z, fx_z, fy_z, fz_z, fxx_z, fxy_z, fxz_z, fyy_z, fyz_z, fzz_z : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double fg_0 = -2.0 * bexp;

                        double fg_1 = f0_0 * fg_0;

                        double fg_2 = fg_1 * fg_0;

                        // leading p_x

                        f0_x[g] += f0_0 * dx;

                        fx_x[g] += f0_0 * (1.0 + dx * dx * fg_0);

                        fy_x[g] += f0_0 * dx * dy * fg_0;

                        fz_x[g] += f0_0 * dx * dz * fg_0;

                        fxx_x[g] += fg_2 * dx * dx * dx + 3.0 * fg_1 * dx;

                        fxy_x[g] += fg_2 * dx * dx * dy + fg_1 * dy;

                        fxz_x[g] += fg_2 * dx * dx * dz + fg_1 * dz;

                        fyy_x[g] += fg_2 * dx * dy * dy + fg_1 * dx;

                        fyz_x[g] += fg_2 * dx * dy * dz;

                        fzz_x[g] += fg_2 * dx * dz * dz + fg_1 * dx;

                        // leading p_y

                        f0_y[g] += f0_0 * dy;

                        fx_y[g] += f0_0 * dy * dx * fg_0;

                        fy_y[g] += f0_0 * (1.0 + dy * dy * fg_0);

                        fz_y[g] += f0_0 * dy * dz * fg_0;

                        fxx_y[g] += fg_2 * dy * dx * dx + fg_1 * dy;

                        fxy_y[g] += fg_2 * dy * dy * dx + fg_1 * dx;

                        fxz_y[g] += fg_2 * dy * dx * dz;

                        fyy_y[g] += fg_2 * dy * dy * dy + 3.0 * fg_1 * dy;

                        fyz_y[g] += fg_2 * dy * dy * dz + fg_1 * dz;

                        fzz_y[g] += fg_2 * dy * dz * dz + fg_1 * dy;

                        // leading p_z

                        f0_z[g] += f0_0 * dz;

                        fx_z[g] += f0_0 * dz * dx * fg_0;

                        fy_z[g] += f0_0 * dz * dy * fg_0;

                        fz_z[g] += f0_0 * (1.0 + dz * dz * fg_0);

                        fxx_z[g] += fg_2 * dz * dx * dx + fg_1 * dz;

                        fxy_z[g] += fg_2 * dz * dy * dx;

                        fxz_z[g] += fg_2 * dz * dz * dx + fg_1 * dx;

                        fyy_z[g] += fg_2 * dz * dy * dy + fg_1 * dz;

                        fyz_z[g] += fg_2 * dz * dz * dy + fg_1 * dy;

                        fzz_z[g] += fg_2 * dz * dz * dz + 3.0 * fg_1 * dz;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(1), 0, 0, nGridPoints, 10);
            }
            else if (bang == 2)
            {
                // d-type GTOs on grid

                bspherbuff.zero();

                bcartbuff.zero();

                auto f0_xx = bcartbuff.data(0);

                auto fx_xx = bcartbuff.data(1);

                auto fy_xx = bcartbuff.data(2);

                auto fz_xx = bcartbuff.data(3);

                auto fxx_xx = bcartbuff.data(4);

                auto fxy_xx = bcartbuff.data(5);

                auto fxz_xx = bcartbuff.data(6);

                auto fyy_xx = bcartbuff.data(7);

                auto fyz_xx = bcartbuff.data(8);

                auto fzz_xx = bcartbuff.data(9);

                auto f0_xy = bcartbuff.data(10);

                auto fx_xy = bcartbuff.data(11);

                auto fy_xy = bcartbuff.data(12);

                auto fz_xy = bcartbuff.data(13);

                auto fxx_xy = bcartbuff.data(14);

                auto fxy_xy = bcartbuff.data(15);

                auto fxz_xy = bcartbuff.data(16);

                auto fyy_xy = bcartbuff.data(17);

                auto fyz_xy = bcartbuff.data(18);

                auto fzz_xy = bcartbuff.data(19);

                auto f0_xz = bcartbuff.data(20);

                auto fx_xz = bcartbuff.data(21);

                auto fy_xz = bcartbuff.data(22);

                auto fz_xz = bcartbuff.data(23);

                auto fxx_xz = bcartbuff.data(24);

                auto fxy_xz = bcartbuff.data(25);

                auto fxz_xz = bcartbuff.data(26);

                auto fyy_xz = bcartbuff.data(27);

                auto fyz_xz = bcartbuff.data(28);

                auto fzz_xz = bcartbuff.data(29);

                auto f0_yy = bcartbuff.data(30);

                auto fx_yy = bcartbuff.data(31);

                auto fy_yy = bcartbuff.data(32);

                auto fz_yy = bcartbuff.data(33);

                auto fxx_yy = bcartbuff.data(34);

                auto fxy_yy = bcartbuff.data(35);

                auto fxz_yy = bcartbuff.data(36);

                auto fyy_yy = bcartbuff.data(37);

                auto fyz_yy = bcartbuff.data(38);

                auto fzz_yy = bcartbuff.data(39);

                auto f0_yz = bcartbuff.data(40);

                auto fx_yz = bcartbuff.data(41);

                auto fy_yz = bcartbuff.data(42);

                auto fz_yz = bcartbuff.data(43);

                auto fxx_yz = bcartbuff.data(44);

                auto fxy_yz = bcartbuff.data(45);

                auto fxz_yz = bcartbuff.data(46);

                auto fyy_yz = bcartbuff.data(47);

                auto fyz_yz = bcartbuff.data(48);

                auto fzz_yz = bcartbuff.data(49);

                auto f0_zz = bcartbuff.data(50);

                auto fx_zz = bcartbuff.data(51);

                auto fy_zz = bcartbuff.data(52);

                auto fz_zz = bcartbuff.data(53);

                auto fxx_zz = bcartbuff.data(54);

                auto fxy_zz = bcartbuff.data(55);

                auto fxz_zz = bcartbuff.data(56);

                auto fyy_zz = bcartbuff.data(57);

                auto fyz_zz = bcartbuff.data(58);

                auto fzz_zz = bcartbuff.data(59);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];

                    auto ry = bfy[iprim];

                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                            f0_xx, fx_xx, fy_xx, fz_xx, fxx_xx, fxy_xx, fxz_xx, fyy_xx, fyz_xx, fzz_xx, \
                            f0_xy, fx_xy, fy_xy, fz_xy, fxx_xy, fxy_xy, fxz_xy, fyy_xy, fyz_xy, fzz_xy, \
                            f0_xz, fx_xz, fy_xz, fz_xz, fxx_xz, fxy_xz, fxz_xz, fyy_xz, fyz_xz, fzz_xz, \
                            f0_yy, fx_yy, fy_yy, fz_yy, fxx_yy, fxy_yy, fxz_yy, fyy_yy, fyz_yy, fzz_yy, \
                            f0_yz, fx_yz, fy_yz, fz_yz, fxx_yz, fxy_yz, fxz_yz, fyy_yz, fyz_yz, fzz_yz, \
                            f0_zz, fx_zz, fy_zz, fz_zz, fxx_zz, fxy_zz, fxz_zz, fyy_zz, fyz_zz, fzz_zz : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;

                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;

                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));

                        double fg_0 = -2.0 * bexp;

                        double fg_1 = f0_0 * fg_0;

                        double fg_2 = fg_1 * fg_0;

                        // leading xx component

                        double f0_x = dx * f0_0;

                        f0_xx[g] += f0_x * dx;

                        fx_xx[g] += f0_x * (2.0 + fg_0 * dx * dx);

                        fy_xx[g] += f0_x * fg_0 * dx * dy;

                        fz_xx[g] += f0_x * fg_0 * dx * dz;

                        fxx_xx[g] += fg_2 * dx * dx * dx * dx + 5.0 * fg_1 * dx * dx + 2.0 * f0_0;

                        fxy_xx[g] += fg_2 * dx * dx * dx * dy + 2.0 * fg_1 * dx * dy;

                        fxz_xx[g] += fg_2 * dx * dx * dx * dz + 2.0 * fg_1 * dx * dz;

                        fyy_xx[g] += fg_2 * dx * dx * dy * dy + 1.0 * fg_1 * dx * dx;

                        fyz_xx[g] += fg_2 * dx * dx * dy * dz;

                        fzz_xx[g] += fg_2 * dx * dx * dz * dz + 1.0 * fg_1 * dx * dx;

                        // leading xy component

                        double f0_y = dy * f0_0;

                        f0_xy[g] += f0_x * dy;

                        fx_xy[g] += f0_y * (1.0 + fg_0 * dx * dx);

                        fy_xy[g] += f0_x * (1.0 + fg_0 * dy * dy);

                        fz_xy[g] += f0_x * fg_0 * dy * dz;

                        fxx_xy[g] += fg_2 * dx * dx * dx * dy + 3.0 * fg_1 * dx * dy;

                        fxy_xy[g] += fg_2 * dx * dx * dy * dy + 1.0 * fg_1 * (dx * dx + dy * dy) + f0_0;

                        fxz_xy[g] += fg_2 * dx * dx * dy * dz + 1.0 * fg_1 * dy * dz;

                        fyy_xy[g] += fg_2 * dy * dy * dy * dx + 3.0 * fg_1 * dy * dx;

                        fyz_xy[g] += fg_2 * dy * dy * dx * dz + 1.0 * fg_1 * dx * dz;

                        fzz_xy[g] += fg_2 * dx * dy * dz * dz + 1.0 * fg_1 * dx * dy;

                        // leading xz component

                        double f0_z = dz * f0_0;

                        f0_xz[g] += f0_x * dz;

                        fx_xz[g] += f0_z * (1.0 + fg_0 * dx * dx);

                        fy_xz[g] += f0_x * fg_0 * dz * dy;

                        fz_xz[g] += f0_x * (1.0 + fg_0 * dz * dz);

                        fxx_xz[g] += fg_2 * dx * dx * dx * dz + 3.0 * fg_1 * dx * dz;

                        fxy_xz[g] += fg_2 * dx * dx * dz * dy + 1.0 * fg_1 * dz * dy;

                        fxz_xz[g] += fg_2 * dx * dx * dz * dz + 1.0 * fg_1 * (dx * dx + dz * dz) + f0_0;

                        fyy_xz[g] += fg_2 * dx * dz * dy * dy + 1.0 * fg_1 * dx * dz;

                        fyz_xz[g] += fg_2 * dz * dz * dx * dy + 1.0 * fg_1 * dx * dy;

                        fzz_xz[g] += fg_2 * dz * dz * dz * dx + 3.0 * fg_1 * dz * dx;

                        // leading yy component

                        f0_yy[g] += f0_y * dy;

                        fx_yy[g] += f0_y * fg_0 * dy * dx;

                        fy_yy[g] += f0_y * (2.0 + fg_0 * dy * dy);

                        fz_yy[g] += f0_y * fg_0 * dy * dz;

                        fxx_yy[g] += fg_2 * dy * dy * dx * dx + 1.0 * fg_1 * dy * dy;

                        fxy_yy[g] += fg_2 * dy * dy * dy * dx + 2.0 * fg_1 * dy * dx;

                        fxz_yy[g] += fg_2 * dy * dy * dx * dz;

                        fyy_yy[g] += fg_2 * dy * dy * dy * dy + 5.0 * fg_1 * dy * dy + 2.0 * f0_0;

                        fyz_yy[g] += fg_2 * dy * dy * dy * dz + 2.0 * fg_1 * dy * dz;

                        fzz_yy[g] += fg_2 * dy * dy * dz * dz + 1.0 * fg_1 * dy * dy;

                        // leading yz component

                        f0_yz[g] += f0_y * dz;

                        fx_yz[g] += f0_y * fg_0 * dz * dx;

                        fy_yz[g] += f0_z * (1.0 + fg_0 * dy * dy);

                        fz_yz[g] += f0_y * (1.0 + fg_0 * dz * dz);

                        fxx_yz[g] += fg_2 * dz * dy * dx * dx + 1.0 * fg_1 * dz * dy;

                        fxy_yz[g] += fg_2 * dy * dy * dz * dx + 1.0 * fg_1 * dz * dx;

                        fxz_yz[g] += fg_2 * dz * dz * dy * dx + 1.0 * fg_1 * dy * dx;

                        fyy_yz[g] += fg_2 * dy * dy * dy * dz + 3.0 * fg_1 * dy * dz;

                        fyz_yz[g] += fg_2 * dz * dz * dy * dy + 1.0 * fg_1 * (dz * dz + dy * dy) + f0_0;

                        fzz_yz[g] += fg_2 * dz * dz * dz * dy + 3.0 * fg_1 * dz * dy;

                        // leading zz component

                        f0_zz[g] += f0_z * dz;

                        fx_zz[g] += f0_z * fg_0 * dz * dx;

                        fy_zz[g] += f0_z * fg_0 * dz * dy;

                        fz_zz[g] += f0_z * (2.0 + fg_0 * dz * dz);

                        fxx_zz[g] += fg_2 * dz * dz * dx * dx + 1.0 * fg_1 * dz * dz;

                        fxy_zz[g] += fg_2 * dz * dz * dy * dx;

                        fxz_zz[g] += fg_2 * dz * dz * dz * dx + 2.0 * fg_1 * dz * dx;

                        fyy_zz[g] += fg_2 * dz * dz * dy * dy + 1.0 * fg_1 * dz * dz;

                        fyz_zz[g] += fg_2 * dz * dz * dz * dy + 2.0 * fg_1 * dz * dy;

                        fzz_zz[g] += fg_2 * dz * dz * dz * dz + 5.0 * fg_1 * dz * dz + 2.0 * f0_0;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(2), 0, 0, nGridPoints, 10);
            }
            else
            {
                // FIX ME: implement l > 2 cases
            }

            // distribute j-th GTO values into grid values matrix

            for (int32_t k = 0; k < bnspher; k++)
            {
                auto bgaos = bspherbuff.data(10 * k);

                auto bgaox = bspherbuff.data(10 * k + 1);

                auto bgaoy = bspherbuff.data(10 * k + 2);

                auto bgaoz = bspherbuff.data(10 * k + 3);

                auto bgaoxx = bspherbuff.data(10 * k + 4);

                auto bgaoxy = bspherbuff.data(10 * k + 5);

                auto bgaoxz = bspherbuff.data(10 * k + 6);

                auto bgaoyy = bspherbuff.data(10 * k + 7);

                auto bgaoyz = bspherbuff.data(10 * k + 8);

                auto bgaozz = bspherbuff.data(10 * k + 9);

                auto idx = bgtos.getIdentifiers(k)[j];

                auto gvals = gtoValues.data(idx);

                auto gvalx = gtoValuesX.data(idx);

                auto gvaly = gtoValuesY.data(idx);

                auto gvalz = gtoValuesZ.data(idx);

                auto gvalxx = gtoValuesXX.data(idx);

                auto gvalxy = gtoValuesXY.data(idx);

                auto gvalxz = gtoValuesXZ.data(idx);

                auto gvalyy = gtoValuesYY.data(idx);

                auto gvalyz = gtoValuesYZ.data(idx);

                auto gvalzz = gtoValuesZZ.data(idx);

                #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz : VLX_ALIGN)
                for (int32_t g = 0; g < nGridPoints; g++)
                {
                    gvals[gridOffset + g] = bgaos[g];

                    gvalx[gridOffset + g] = bgaox[g];

                    gvaly[gridOffset + g] = bgaoy[g];

                    gvalz[gridOffset + g] = bgaoz[g];

                    gvalxx[gridOffset + g] = bgaoxx[g];

                    gvalxy[gridOffset + g] = bgaoxy[g];

                    gvalxz[gridOffset + g] = bgaoxz[g];

                    gvalyy[gridOffset + g] = bgaoyy[g];

                    gvalyz[gridOffset + g] = bgaoyz[g];

                    gvalzz[gridOffset + g] = bgaozz[g];
                }
            }
        }
    }
}

}  // namespace gtoeval
