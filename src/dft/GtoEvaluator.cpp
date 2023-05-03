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

                // 3rd-order derivative
                // gto_m3: ang (ang-1) (ang-2) r^{ang-3} |C| exp(-alpha r^2)
                // gto_m1: 2alpha 3ang^2       r^{ang-1} |C| exp(-alpha r^2)
                // gto_p1: 4alpha^2 (3ang+3)   r^{ang+1} |C| exp(-alpha r^2)
                // gto_p3: 8alpha^3            r^{ang+3} |C| exp(-alpha r^2)

                auto gtolimit_base = maxcoef * std::exp(-minexp * r2);

                auto r = std::sqrt(r2);

                for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit_base *= r;

                auto gtolimit = gtolimit_base;

                if (gtoDeriv > 0)
                {
                    auto gtolimit_p = 2.0 * maxexp * r * gtolimit_base;  // gto_p

                    gtolimit = std::max(gtolimit, gtolimit_p);

                    if (bang > 0)
                    {
                        auto gtolimit_m = gtolimit_base / r * bang;  // gto_m

                        gtolimit = std::max(gtolimit, gtolimit_m);
                    }
                }

                if (gtoDeriv > 1)
                {
                    auto gtolimit_p2 = 4.0 * maxexp * maxexp * r2 * gtolimit_base;  // gto_p2

                    auto gtolimit_0 = 2.0 * maxexp * (2 * bang + 1) * gtolimit_base;  // gto

                    gtolimit = std::max({gtolimit, gtolimit_0, gtolimit_p2});

                    if (bang > 1)
                    {
                        auto gtolimit_m2 = gtolimit_base / r2 * bang * (bang - 1);  // gto_m2

                        gtolimit = std::max(gtolimit, gtolimit_m2);
                    }
                }

                if (gtoDeriv > 2)
                {
                    auto r3 = r * r2;

                    auto gtolimit_p3 = 8.0 * maxexp * maxexp * maxexp * r3 * gtolimit_base;  // gto_p3

                    auto gtolimit_p1 = 4.0 * maxexp * maxexp * (3 * bang + 3) * r * gtolimit_base;  // gto_p1

                    gtolimit = std::max({gtolimit, gtolimit_p1, gtolimit_p3});

                    if (bang > 0)
                    {
                        auto gtolimit_m1 = 2.0 * maxexp * gtolimit_base / r * (3 * bang * bang);  // gto_m1

                        gtolimit = std::max(gtolimit, gtolimit_m1);
                    }

                    if (bang > 2)
                    {
                        auto gtolimit_m3 = gtolimit_base / r3 * bang * (bang - 1) * (bang - 2);  // gto_m3

                        gtolimit = std::max(gtolimit, gtolimit_m3);
                    }
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

                std::string errangmom("gtoeval::computeGtosValuesForLDA: Only implemented up to d-orbitals");

                errors::assertMsgCritical(false, errangmom);
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

                std::string errangmom("gtoeval::computeGtosValuesForGGA: Only implemented up to d-orbitals");

                errors::assertMsgCritical(false, errangmom);
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

                std::string errangmom("gtoeval::computeGtosValuesForMetaGGA: Only implemented up to d-orbitals");

                errors::assertMsgCritical(false, errangmom);
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

void
computeGtosValuesForThirdOrder(CMemBlock2D<double>&      gtoValues,
                               CMemBlock2D<double>&      gtoValuesX,
                               CMemBlock2D<double>&      gtoValuesY,
                               CMemBlock2D<double>&      gtoValuesZ,
                               CMemBlock2D<double>&      gtoValuesXX,
                               CMemBlock2D<double>&      gtoValuesXY,
                               CMemBlock2D<double>&      gtoValuesXZ,
                               CMemBlock2D<double>&      gtoValuesYY,
                               CMemBlock2D<double>&      gtoValuesYZ,
                               CMemBlock2D<double>&      gtoValuesZZ,
                               CMemBlock2D<double>&      gtoValuesXXX,
                               CMemBlock2D<double>&      gtoValuesXXY,
                               CMemBlock2D<double>&      gtoValuesXXZ,
                               CMemBlock2D<double>&      gtoValuesXYY,
                               CMemBlock2D<double>&      gtoValuesXYZ,
                               CMemBlock2D<double>&      gtoValuesXZZ,
                               CMemBlock2D<double>&      gtoValuesYYY,
                               CMemBlock2D<double>&      gtoValuesYYZ,
                               CMemBlock2D<double>&      gtoValuesYZZ,
                               CMemBlock2D<double>&      gtoValuesZZZ,
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

        int32_t nvcomp = 20;  // up to third order derivatives (1+3+6+10)

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

                auto fxxx_0 = bspherbuff.data(10);
                auto fxxy_0 = bspherbuff.data(11);
                auto fxxz_0 = bspherbuff.data(12);
                auto fxyy_0 = bspherbuff.data(13);
                auto fxyz_0 = bspherbuff.data(14);
                auto fxzz_0 = bspherbuff.data(15);
                auto fyyy_0 = bspherbuff.data(16);
                auto fyyz_0 = bspherbuff.data(17);
                auto fyzz_0 = bspherbuff.data(18);
                auto fzzz_0 = bspherbuff.data(19);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];
                    auto bnorm = bfnorms[iprim];

                    auto rx = bfx[iprim];
                    auto ry = bfy[iprim];
                    auto rz = bfz[iprim];

                    #pragma omp simd aligned(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, \
                            f0_0, fx_0, fy_0, fz_0, fxx_0, fxy_0, fxz_0, fyy_0, fyz_0, fzz_0, \
                            fxxx_0, fxxy_0, fxxz_0, fxyy_0, fxyz_0, fxzz_0, fyyy_0, fyyz_0, fyzz_0, fzzz_0 : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;
                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;
                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double g0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));  // g
                        double g1 = -2.0 * bexp * g0;  // -2*a*g
                        double g2 = -2.0 * bexp * g1;  // 4*a*a*g
                        double g3 = -2.0 * bexp * g2;  // -8*a*a*a*g

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

                        fxxx_0[g] += g3 * dx * dx * dx + 3.0 * g2 * dx;
                        fxxy_0[g] += g3 * dx * dx * dy + g2 * dy;
                        fxxz_0[g] += g3 * dx * dx * dz + g2 * dz;
                        fxyy_0[g] += g3 * dx * dy * dy + g2 * dx;
                        fxyz_0[g] += g3 * dx * dy * dz;
                        fxzz_0[g] += g3 * dx * dz * dz + g2 * dx;
                        fyyy_0[g] += g3 * dy * dy * dy + 3.0 * g2 * dy;
                        fyyz_0[g] += g3 * dy * dy * dz + g2 * dz;
                        fyzz_0[g] += g3 * dy * dz * dz + g2 * dy;
                        fzzz_0[g] += g3 * dz * dz * dz + 3.0 * g2 * dz;
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

                auto fxxx_x = bcartbuff.data(10);
                auto fxxy_x = bcartbuff.data(11);
                auto fxxz_x = bcartbuff.data(12);
                auto fxyy_x = bcartbuff.data(13);
                auto fxyz_x = bcartbuff.data(14);
                auto fxzz_x = bcartbuff.data(15);
                auto fyyy_x = bcartbuff.data(16);
                auto fyyz_x = bcartbuff.data(17);
                auto fyzz_x = bcartbuff.data(18);
                auto fzzz_x = bcartbuff.data(19);

                auto f0_y = bcartbuff.data(20);

                auto fx_y = bcartbuff.data(21);
                auto fy_y = bcartbuff.data(22);
                auto fz_y = bcartbuff.data(23);

                auto fxx_y = bcartbuff.data(24);
                auto fxy_y = bcartbuff.data(25);
                auto fxz_y = bcartbuff.data(26);
                auto fyy_y = bcartbuff.data(27);
                auto fyz_y = bcartbuff.data(28);
                auto fzz_y = bcartbuff.data(29);

                auto fxxx_y = bcartbuff.data(30);
                auto fxxy_y = bcartbuff.data(31);
                auto fxxz_y = bcartbuff.data(32);
                auto fxyy_y = bcartbuff.data(33);
                auto fxyz_y = bcartbuff.data(34);
                auto fxzz_y = bcartbuff.data(35);
                auto fyyy_y = bcartbuff.data(36);
                auto fyyz_y = bcartbuff.data(37);
                auto fyzz_y = bcartbuff.data(38);
                auto fzzz_y = bcartbuff.data(39);

                auto f0_z = bcartbuff.data(40);

                auto fx_z = bcartbuff.data(41);
                auto fy_z = bcartbuff.data(42);
                auto fz_z = bcartbuff.data(43);

                auto fxx_z = bcartbuff.data(44);
                auto fxy_z = bcartbuff.data(45);
                auto fxz_z = bcartbuff.data(46);
                auto fyy_z = bcartbuff.data(47);
                auto fyz_z = bcartbuff.data(48);
                auto fzz_z = bcartbuff.data(49);

                auto fxxx_z = bcartbuff.data(50);
                auto fxxy_z = bcartbuff.data(51);
                auto fxxz_z = bcartbuff.data(52);
                auto fxyy_z = bcartbuff.data(53);
                auto fxyz_z = bcartbuff.data(54);
                auto fxzz_z = bcartbuff.data(55);
                auto fyyy_z = bcartbuff.data(56);
                auto fyyz_z = bcartbuff.data(57);
                auto fyzz_z = bcartbuff.data(58);
                auto fzzz_z = bcartbuff.data(59);

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
                            f0_z, fx_z, fy_z, fz_z, fxx_z, fxy_z, fxz_z, fyy_z, fyz_z, fzz_z, \
                            fxxx_x, fxxy_x, fxxz_x, fxyy_x, fxyz_x, fxzz_x, fyyy_x, fyyz_x, fyzz_x, fzzz_x, \
                            fxxx_y, fxxy_y, fxxz_y, fxyy_y, fxyz_y, fxzz_y, fyyy_y, fyyz_y, fyzz_y, fzzz_y, \
                            fxxx_z, fxxy_z, fxxz_z, fxyy_z, fxyz_z, fxzz_z, fyyy_z, fyyz_z, fyzz_z, fzzz_z : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;
                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;
                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                        double fg_0 = -2.0 * bexp;
                        double fg_1 = f0_0 * fg_0;  // -2*a*g
                        double fg_2 = fg_1 * fg_0;  // 4*a*a*g
                        double fg_3 = fg_2 * fg_0;  // -8*a*a*a*g

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

                        fxxx_x[g] += fg_3 * dx * dx * dx * dx + 6.0 * fg_2 * dx * dx + 3.0 * fg_1;
                        fxxy_x[g] += fg_3 * dx * dx * dx * dy + 3.0 * fg_2 * dx * dy;
                        fxxz_x[g] += fg_3 * dx * dx * dx * dz + 3.0 * fg_2 * dx * dz;
                        fxyy_x[g] += fg_3 * dx * dx * dy * dy + fg_2 * dx * dx + fg_2 * dy * dy + fg_1;
                        fxyz_x[g] += fg_3 * dx * dx * dy * dz + fg_2 * dy * dz;
                        fxzz_x[g] += fg_3 * dx * dx * dz * dz + fg_2 * dx * dx + fg_2 * dz * dz + fg_1;
                        fyyy_x[g] += fg_3 * dx * dy * dy * dy + 3.0 * fg_2 * dx * dy;
                        fyyz_x[g] += fg_3 * dx * dy * dy * dz + fg_2 * dx * dz;
                        fyzz_x[g] += fg_3 * dx * dy * dz * dz + fg_2 * dx * dy;
                        fzzz_x[g] += fg_3 * dx * dz * dz * dz + 3.0 * fg_2 * dx * dz;

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

                        fxxx_y[g] += fg_3 * dx * dx * dx * dy + 3.0 * fg_2 * dx * dy;
                        fxxy_y[g] += fg_3 * dx * dx * dy * dy + fg_2 * dx * dx + fg_2 * dy * dy + fg_1;
                        fxxz_y[g] += fg_3 * dx * dx * dy * dz + fg_2 * dy * dz;
                        fxyy_y[g] += fg_3 * dx * dy * dy * dy + 3.0 * fg_2 * dx * dy;
                        fxyz_y[g] += fg_3 * dx * dy * dy * dz + fg_2 * dx * dz;
                        fxzz_y[g] += fg_3 * dx * dy * dz * dz + fg_2 * dx * dy;
                        fyyy_y[g] += fg_3 * dy * dy * dy * dy + 6.0 * fg_2 * dy * dy + 3.0 * fg_1;
                        fyyz_y[g] += fg_3 * dy * dy * dy * dz + 3.0 * fg_2 * dy * dz;
                        fyzz_y[g] += fg_3 * dy * dy * dz * dz + fg_2 * dy * dy + fg_2 * dz * dz + fg_1;
                        fzzz_y[g] += fg_3 * dy * dz * dz * dz + 3.0 * fg_2 * dy * dz;

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

                        fxxx_z[g] += fg_3 * dx * dx * dx * dz + 3.0 * fg_2 * dx * dz;
                        fxxy_z[g] += fg_3 * dx * dx * dy * dz + fg_2 * dy * dz;
                        fxxz_z[g] += fg_3 * dx * dx * dz * dz + fg_2 * dx * dx + fg_2 * dz * dz + fg_1;
                        fxyy_z[g] += fg_3 * dx * dy * dy * dz + fg_2 * dx * dz;
                        fxyz_z[g] += fg_3 * dx * dy * dz * dz + fg_2 * dx * dy;
                        fxzz_z[g] += fg_3 * dx * dz * dz * dz + 3.0 * fg_2 * dx * dz;
                        fyyy_z[g] += fg_3 * dy * dy * dy * dz + 3.0 * fg_2 * dy * dz;
                        fyyz_z[g] += fg_3 * dy * dy * dz * dz + fg_2 * dy * dy + fg_2 * dz * dz + fg_1;
                        fyzz_z[g] += fg_3 * dy * dz * dz * dz + 3.0 * fg_2 * dy * dz;
                        fzzz_z[g] += fg_3 * dz * dz * dz * dz + 6.0 * fg_2 * dz * dz + 3.0 * fg_1;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(1), 0, 0, nGridPoints, 20);
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

                auto fxxx_xx = bcartbuff.data(10);
                auto fxxy_xx = bcartbuff.data(11);
                auto fxxz_xx = bcartbuff.data(12);
                auto fxyy_xx = bcartbuff.data(13);
                auto fxyz_xx = bcartbuff.data(14);
                auto fxzz_xx = bcartbuff.data(15);
                auto fyyy_xx = bcartbuff.data(16);
                auto fyyz_xx = bcartbuff.data(17);
                auto fyzz_xx = bcartbuff.data(18);
                auto fzzz_xx = bcartbuff.data(19);

                auto f0_xy = bcartbuff.data(20);

                auto fx_xy = bcartbuff.data(21);
                auto fy_xy = bcartbuff.data(22);
                auto fz_xy = bcartbuff.data(23);

                auto fxx_xy = bcartbuff.data(24);
                auto fxy_xy = bcartbuff.data(25);
                auto fxz_xy = bcartbuff.data(26);
                auto fyy_xy = bcartbuff.data(27);
                auto fyz_xy = bcartbuff.data(28);
                auto fzz_xy = bcartbuff.data(29);

                auto fxxx_xy = bcartbuff.data(30);
                auto fxxy_xy = bcartbuff.data(31);
                auto fxxz_xy = bcartbuff.data(32);
                auto fxyy_xy = bcartbuff.data(33);
                auto fxyz_xy = bcartbuff.data(34);
                auto fxzz_xy = bcartbuff.data(35);
                auto fyyy_xy = bcartbuff.data(36);
                auto fyyz_xy = bcartbuff.data(37);
                auto fyzz_xy = bcartbuff.data(38);
                auto fzzz_xy = bcartbuff.data(39);

                auto f0_xz = bcartbuff.data(40);

                auto fx_xz = bcartbuff.data(41);
                auto fy_xz = bcartbuff.data(42);
                auto fz_xz = bcartbuff.data(43);

                auto fxx_xz = bcartbuff.data(44);
                auto fxy_xz = bcartbuff.data(45);
                auto fxz_xz = bcartbuff.data(46);
                auto fyy_xz = bcartbuff.data(47);
                auto fyz_xz = bcartbuff.data(48);
                auto fzz_xz = bcartbuff.data(49);

                auto fxxx_xz = bcartbuff.data(50);
                auto fxxy_xz = bcartbuff.data(51);
                auto fxxz_xz = bcartbuff.data(52);
                auto fxyy_xz = bcartbuff.data(53);
                auto fxyz_xz = bcartbuff.data(54);
                auto fxzz_xz = bcartbuff.data(55);
                auto fyyy_xz = bcartbuff.data(56);
                auto fyyz_xz = bcartbuff.data(57);
                auto fyzz_xz = bcartbuff.data(58);
                auto fzzz_xz = bcartbuff.data(59);

                auto f0_yy = bcartbuff.data(60);

                auto fx_yy = bcartbuff.data(61);
                auto fy_yy = bcartbuff.data(62);
                auto fz_yy = bcartbuff.data(63);

                auto fxx_yy = bcartbuff.data(64);
                auto fxy_yy = bcartbuff.data(65);
                auto fxz_yy = bcartbuff.data(66);
                auto fyy_yy = bcartbuff.data(67);
                auto fyz_yy = bcartbuff.data(68);
                auto fzz_yy = bcartbuff.data(69);

                auto fxxx_yy = bcartbuff.data(70);
                auto fxxy_yy = bcartbuff.data(71);
                auto fxxz_yy = bcartbuff.data(72);
                auto fxyy_yy = bcartbuff.data(73);
                auto fxyz_yy = bcartbuff.data(74);
                auto fxzz_yy = bcartbuff.data(75);
                auto fyyy_yy = bcartbuff.data(76);
                auto fyyz_yy = bcartbuff.data(77);
                auto fyzz_yy = bcartbuff.data(78);
                auto fzzz_yy = bcartbuff.data(79);

                auto f0_yz = bcartbuff.data(80);

                auto fx_yz = bcartbuff.data(81);
                auto fy_yz = bcartbuff.data(82);
                auto fz_yz = bcartbuff.data(83);

                auto fxx_yz = bcartbuff.data(84);
                auto fxy_yz = bcartbuff.data(85);
                auto fxz_yz = bcartbuff.data(86);
                auto fyy_yz = bcartbuff.data(87);
                auto fyz_yz = bcartbuff.data(88);
                auto fzz_yz = bcartbuff.data(89);

                auto fxxx_yz = bcartbuff.data(90);
                auto fxxy_yz = bcartbuff.data(91);
                auto fxxz_yz = bcartbuff.data(92);
                auto fxyy_yz = bcartbuff.data(93);
                auto fxyz_yz = bcartbuff.data(94);
                auto fxzz_yz = bcartbuff.data(95);
                auto fyyy_yz = bcartbuff.data(96);
                auto fyyz_yz = bcartbuff.data(97);
                auto fyzz_yz = bcartbuff.data(98);
                auto fzzz_yz = bcartbuff.data(99);

                auto f0_zz = bcartbuff.data(100);

                auto fx_zz = bcartbuff.data(101);
                auto fy_zz = bcartbuff.data(102);
                auto fz_zz = bcartbuff.data(103);

                auto fxx_zz = bcartbuff.data(104);
                auto fxy_zz = bcartbuff.data(105);
                auto fxz_zz = bcartbuff.data(106);
                auto fyy_zz = bcartbuff.data(107);
                auto fyz_zz = bcartbuff.data(108);
                auto fzz_zz = bcartbuff.data(109);

                auto fxxx_zz = bcartbuff.data(110);
                auto fxxy_zz = bcartbuff.data(111);
                auto fxxz_zz = bcartbuff.data(112);
                auto fxyy_zz = bcartbuff.data(113);
                auto fxyz_zz = bcartbuff.data(114);
                auto fxzz_zz = bcartbuff.data(115);
                auto fyyy_zz = bcartbuff.data(116);
                auto fyyz_zz = bcartbuff.data(117);
                auto fyzz_zz = bcartbuff.data(118);
                auto fzzz_zz = bcartbuff.data(119);

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
                            f0_zz, fx_zz, fy_zz, fz_zz, fxx_zz, fxy_zz, fxz_zz, fyy_zz, fyz_zz, fzz_zz, \
                            fxxx_xx, fxxy_xx, fxxz_xx, fxyy_xx, fxyz_xx, fxzz_xx, fyyy_xx, fyyz_xx, fyzz_xx, fzzz_xx, \
                            fxxx_xy, fxxy_xy, fxxz_xy, fxyy_xy, fxyz_xy, fxzz_xy, fyyy_xy, fyyz_xy, fyzz_xy, fzzz_xy, \
                            fxxx_xz, fxxy_xz, fxxz_xz, fxyy_xz, fxyz_xz, fxzz_xz, fyyy_xz, fyyz_xz, fyzz_xz, fzzz_xz, \
                            fxxx_yy, fxxy_yy, fxxz_yy, fxyy_yy, fxyz_yy, fxzz_yy, fyyy_yy, fyyz_yy, fyzz_yy, fzzz_yy, \
                            fxxx_yz, fxxy_yz, fxxz_yz, fxyy_yz, fxyz_yz, fxzz_yz, fyyy_yz, fyyz_yz, fyzz_yz, fzzz_yz, \
                            fxxx_zz, fxxy_zz, fxxz_zz, fxyy_zz, fxyz_zz, fxzz_zz, fyyy_zz, fyyz_zz, fyzz_zz, fzzz_zz : VLX_ALIGN)
                    for (int32_t g = 0; g < nGridPoints; g++)
                    {
                        double dx = gridCoordinatesX[gridBlockPosition + gridOffset + g] - rx;
                        double dy = gridCoordinatesY[gridBlockPosition + gridOffset + g] - ry;
                        double dz = gridCoordinatesZ[gridBlockPosition + gridOffset + g] - rz;

                        double f0_0 = bnorm * std::exp(-bexp * (dx * dx + dy * dy + dz * dz));
                        double fg_0 = -2.0 * bexp;
                        double fg_1 = f0_0 * fg_0;
                        double fg_2 = fg_1 * fg_0;
                        double fg_3 = fg_2 * fg_0;  // -8*a*a*a*g

                        // leading xx component

                        double f0_x = dx * f0_0;

                        f0_xx[g] += f0_x * dx;

                        fx_xx[g] += fg_1 * dx * dx * dx + 2.0 * f0_0 * dx;
                        fy_xx[g] += fg_1 * dx * dx * dy;
                        fz_xx[g] += fg_1 * dx * dx * dz;

                        fxx_xx[g] += fg_2 * dx * dx * dx * dx + 5.0 * fg_1 * dx * dx + 2.0 * f0_0;
                        fxy_xx[g] += fg_2 * dx * dx * dx * dy + 2.0 * fg_1 * dx * dy;
                        fxz_xx[g] += fg_2 * dx * dx * dx * dz + 2.0 * fg_1 * dx * dz;
                        fyy_xx[g] += fg_2 * dx * dx * dy * dy + fg_1 * dx * dx;
                        fyz_xx[g] += fg_2 * dx * dx * dy * dz;
                        fzz_xx[g] += fg_2 * dx * dx * dz * dz + fg_1 * dx * dx;

                        fxxx_xx[g] += fg_3 * dx * dx * dx * dx * dx + 9.0 * fg_2 * dx * dx * dx + 12.0 * fg_1 * dx;
                        fxxy_xx[g] += fg_3 * dx * dx * dx * dx * dy + 5.0 * fg_2 * dx * dx * dy + 2.0 * fg_1 * dy;
                        fxxz_xx[g] += fg_3 * dx * dx * dx * dx * dz + 5.0 * fg_2 * dx * dx * dz + 2.0 * fg_1 * dz;
                        fxyy_xx[g] += fg_3 * dx * dx * dx * dy * dy + fg_2 * dx * dx * dx + 2.0 * fg_2 * dx * dy * dy + 2.0 * fg_1 * dx;
                        fxyz_xx[g] += fg_3 * dx * dx * dx * dy * dz + 2.0 * fg_2 * dx * dy * dz;
                        fxzz_xx[g] += fg_3 * dx * dx * dx * dz * dz + fg_2 * dx * dx * dx + 2.0 * fg_2 * dx * dz * dz + 2.0 * fg_1 * dx;
                        fyyy_xx[g] += fg_3 * dx * dx * dy * dy * dy + 3.0 * fg_2 * dx * dx * dy;
                        fyyz_xx[g] += fg_3 * dx * dx * dy * dy * dz + fg_2 * dx * dx * dz;
                        fyzz_xx[g] += fg_3 * dx * dx * dy * dz * dz + fg_2 * dx * dx * dy;
                        fzzz_xx[g] += fg_3 * dx * dx * dz * dz * dz + 3.0 * fg_2 * dx * dx * dz;

                        // leading xy component

                        double f0_y = dy * f0_0;

                        f0_xy[g] += f0_x * dy;

                        fx_xy[g] += fg_1 * dx * dx * dy + f0_0 * dy;
                        fy_xy[g] += fg_1 * dx * dy * dy + f0_0 * dx;
                        fz_xy[g] += fg_1 * dx * dy * dz;

                        fxx_xy[g] += fg_2 * dx * dx * dx * dy + 3.0 * fg_1 * dx * dy;
                        fxy_xy[g] += fg_2 * dx * dx * dy * dy + fg_1 * dx * dx + fg_1 * dy * dy + f0_0;
                        fxz_xy[g] += fg_2 * dx * dx * dy * dz + fg_1 * dy * dz;
                        fyy_xy[g] += fg_2 * dx * dy * dy * dy + 3.0 * fg_1 * dx * dy;
                        fyz_xy[g] += fg_2 * dx * dy * dy * dz + fg_1 * dx * dz;
                        fzz_xy[g] += fg_2 * dx * dy * dz * dz + fg_1 * dx * dy;

                        fxxx_xy[g] += fg_3 * dx * dx * dx * dx * dy + 6.0 * fg_2 * dx * dx * dy + 3.0 * fg_1 * dy;
                        fxxy_xy[g] += fg_3 * dx * dx * dx * dy * dy + fg_2 * dx * dx * dx + 3.0 * fg_2 * dx * dy * dy + 3.0 * fg_1 * dx;
                        fxxz_xy[g] += fg_3 * dx * dx * dx * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fxyy_xy[g] += fg_3 * dx * dx * dy * dy * dy + 3.0 * fg_2 * dx * dx * dy + fg_2 * dy * dy * dy + 3.0 * fg_1 * dy;
                        fxyz_xy[g] += fg_3 * dx * dx * dy * dy * dz + fg_2 * dx * dx * dz + fg_2 * dy * dy * dz + fg_1 * dz;
                        fxzz_xy[g] += fg_3 * dx * dx * dy * dz * dz + fg_2 * dx * dx * dy + fg_2 * dy * dz * dz + fg_1 * dy;
                        fyyy_xy[g] += fg_3 * dx * dy * dy * dy * dy + 6.0 * fg_2 * dx * dy * dy + 3.0 * fg_1 * dx;
                        fyyz_xy[g] += fg_3 * dx * dy * dy * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fyzz_xy[g] += fg_3 * dx * dy * dy * dz * dz + fg_2 * dx * dy * dy + fg_2 * dx * dz * dz + fg_1 * dx;
                        fzzz_xy[g] += fg_3 * dx * dy * dz * dz * dz + 3.0 * fg_2 * dx * dy * dz;

                        // leading xz component

                        double f0_z = dz * f0_0;

                        f0_xz[g] += f0_x * dz;

                        fx_xz[g] += fg_1 * dx * dx * dz + f0_0 * dz;
                        fy_xz[g] += fg_1 * dx * dy * dz;
                        fz_xz[g] += fg_1 * dx * dz * dz + f0_0 * dx;

                        fxx_xz[g] += fg_2 * dx * dx * dx * dz + 3.0 * fg_1 * dx * dz;
                        fxy_xz[g] += fg_2 * dx * dx * dy * dz + fg_1 * dy * dz;
                        fxz_xz[g] += fg_2 * dx * dx * dz * dz + fg_1 * dx * dx + fg_1 * dz * dz + f0_0;
                        fyy_xz[g] += fg_2 * dx * dy * dy * dz + fg_1 * dx * dz;
                        fyz_xz[g] += fg_2 * dx * dy * dz * dz + fg_1 * dx * dy;
                        fzz_xz[g] += fg_2 * dx * dz * dz * dz + 3.0 * fg_1 * dx * dz;

                        fxxx_xz[g] += fg_3 * dx * dx * dx * dx * dz + 6.0 * fg_2 * dx * dx * dz + 3.0 * fg_1 * dz;
                        fxxy_xz[g] += fg_3 * dx * dx * dx * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fxxz_xz[g] += fg_3 * dx * dx * dx * dz * dz + fg_2 * dx * dx * dx + 3.0 * fg_2 * dx * dz * dz + 3.0 * fg_1 * dx;
                        fxyy_xz[g] += fg_3 * dx * dx * dy * dy * dz + fg_2 * dx * dx * dz + fg_2 * dy * dy * dz + fg_1 * dz;
                        fxyz_xz[g] += fg_3 * dx * dx * dy * dz * dz + fg_2 * dx * dx * dy + fg_2 * dy * dz * dz + fg_1 * dy;
                        fxzz_xz[g] += fg_3 * dx * dx * dz * dz * dz + 3.0 * fg_2 * dx * dx * dz + fg_2 * dz * dz * dz + 3.0 * fg_1 * dz;
                        fyyy_xz[g] += fg_3 * dx * dy * dy * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fyyz_xz[g] += fg_3 * dx * dy * dy * dz * dz + fg_2 * dx * dy * dy + fg_2 * dx * dz * dz + fg_1 * dx;
                        fyzz_xz[g] += fg_3 * dx * dy * dz * dz * dz + 3.0 * fg_2 * dx * dy * dz;
                        fzzz_xz[g] += fg_3 * dx * dz * dz * dz * dz + 6.0 * fg_2 * dx * dz * dz + 3.0 * fg_1 * dx;

                        // leading yy component

                        f0_yy[g] += f0_y * dy;

                        fx_yy[g] += fg_1 * dx * dy * dy;
                        fy_yy[g] += fg_1 * dy * dy * dy + 2.0 * f0_0 * dy;
                        fz_yy[g] += fg_1 * dy * dy * dz;

                        fxx_yy[g] += fg_2 * dx * dx * dy * dy + fg_1 * dy * dy;
                        fxy_yy[g] += fg_2 * dx * dy * dy * dy + 2.0 * fg_1 * dx * dy;
                        fxz_yy[g] += fg_2 * dx * dy * dy * dz;
                        fyy_yy[g] += fg_2 * dy * dy * dy * dy + 5.0 * fg_1 * dy * dy + 2.0 * f0_0;
                        fyz_yy[g] += fg_2 * dy * dy * dy * dz + 2.0 * fg_1 * dy * dz;
                        fzz_yy[g] += fg_2 * dy * dy * dz * dz + fg_1 * dy * dy;

                        fxxx_yy[g] += fg_3 * dx * dx * dx * dy * dy + 3.0 * fg_2 * dx * dy * dy;
                        fxxy_yy[g] += fg_3 * dx * dx * dy * dy * dy + 2.0 * fg_2 * dx * dx * dy + fg_2 * dy * dy * dy + 2.0 * fg_1 * dy;
                        fxxz_yy[g] += fg_3 * dx * dx * dy * dy * dz + fg_2 * dy * dy * dz;
                        fxyy_yy[g] += fg_3 * dx * dy * dy * dy * dy + 5.0 * fg_2 * dx * dy * dy + 2.0 * fg_1 * dx;
                        fxyz_yy[g] += fg_3 * dx * dy * dy * dy * dz + 2.0 * fg_2 * dx * dy * dz;
                        fxzz_yy[g] += fg_3 * dx * dy * dy * dz * dz + fg_2 * dx * dy * dy;
                        fyyy_yy[g] += fg_3 * dy * dy * dy * dy * dy + 9.0 * fg_2 * dy * dy * dy + 12.0 * fg_1 * dy;
                        fyyz_yy[g] += fg_3 * dy * dy * dy * dy * dz + 5.0 * fg_2 * dy * dy * dz + 2.0 * fg_1 * dz;
                        fyzz_yy[g] += fg_3 * dy * dy * dy * dz * dz + fg_2 * dy * dy * dy + 2.0 * fg_2 * dy * dz * dz + 2.0 * fg_1 * dy;
                        fzzz_yy[g] += fg_3 * dy * dy * dz * dz * dz + 3.0 * fg_2 * dy * dy * dz;

                        // leading yz component

                        f0_yz[g] += f0_y * dz;

                        fx_yz[g] += fg_1 * dx * dy * dz;
                        fy_yz[g] += fg_1 * dy * dy * dz + f0_0 * dz;
                        fz_yz[g] += fg_1 * dy * dz * dz + f0_0 * dy;

                        fxx_yz[g] += fg_2 * dx * dx * dy * dz + fg_1 * dy * dz;
                        fxy_yz[g] += fg_2 * dx * dy * dy * dz + fg_1 * dx * dz;
                        fxz_yz[g] += fg_2 * dx * dy * dz * dz + fg_1 * dx * dy;
                        fyy_yz[g] += fg_2 * dy * dy * dy * dz + 3.0 * fg_1 * dy * dz;
                        fyz_yz[g] += fg_2 * dy * dy * dz * dz + fg_1 * dy * dy + fg_1 * dz * dz + f0_0;
                        fzz_yz[g] += fg_2 * dy * dz * dz * dz + 3.0 * fg_1 * dy * dz;

                        fxxx_yz[g] += fg_3 * dx * dx * dx * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fxxy_yz[g] += fg_3 * dx * dx * dy * dy * dz + fg_2 * dx * dx * dz + fg_2 * dy * dy * dz + fg_1 * dz;
                        fxxz_yz[g] += fg_3 * dx * dx * dy * dz * dz + fg_2 * dx * dx * dy + fg_2 * dy * dz * dz + fg_1 * dy;
                        fxyy_yz[g] += fg_3 * dx * dy * dy * dy * dz + 3.0 * fg_2 * dx * dy * dz;
                        fxyz_yz[g] += fg_3 * dx * dy * dy * dz * dz + fg_2 * dx * dy * dy + fg_2 * dx * dz * dz + fg_1 * dx;
                        fxzz_yz[g] += fg_3 * dx * dy * dz * dz * dz + 3.0 * fg_2 * dx * dy * dz;
                        fyyy_yz[g] += fg_3 * dy * dy * dy * dy * dz + 6.0 * fg_2 * dy * dy * dz + 3.0 * fg_1 * dz;
                        fyyz_yz[g] += fg_3 * dy * dy * dy * dz * dz + fg_2 * dy * dy * dy + 3.0 * fg_2 * dy * dz * dz + 3.0 * fg_1 * dy;
                        fyzz_yz[g] += fg_3 * dy * dy * dz * dz * dz + 3.0 * fg_2 * dy * dy * dz + fg_2 * dz * dz * dz + 3.0 * fg_1 * dz;
                        fzzz_yz[g] += fg_3 * dy * dz * dz * dz * dz + 6.0 * fg_2 * dy * dz * dz + 3.0 * fg_1 * dy;

                        // leading zz component

                        f0_zz[g] += f0_z * dz;

                        fx_zz[g] += fg_1 * dx * dz * dz;
                        fy_zz[g] += fg_1 * dy * dz * dz;
                        fz_zz[g] += fg_1 * dz * dz * dz + 2.0 * f0_0 * dz;

                        fxx_zz[g] += fg_2 * dx * dx * dz * dz + fg_1 * dz * dz;
                        fxy_zz[g] += fg_2 * dx * dy * dz * dz;
                        fxz_zz[g] += fg_2 * dx * dz * dz * dz + 2.0 * fg_1 * dx * dz;
                        fyy_zz[g] += fg_2 * dy * dy * dz * dz + fg_1 * dz * dz;
                        fyz_zz[g] += fg_2 * dy * dz * dz * dz + 2.0 * fg_1 * dy * dz;
                        fzz_zz[g] += fg_2 * dz * dz * dz * dz + 5.0 * fg_1 * dz * dz + 2.0 * f0_0;

                        fxxx_zz[g] += fg_3 * dx * dx * dx * dz * dz + 3.0 * fg_2 * dx * dz * dz;
                        fxxy_zz[g] += fg_3 * dx * dx * dy * dz * dz + fg_2 * dy * dz * dz;
                        fxxz_zz[g] += fg_3 * dx * dx * dz * dz * dz + 2.0 * fg_2 * dx * dx * dz + fg_2 * dz * dz * dz + 2.0 * fg_1 * dz;
                        fxyy_zz[g] += fg_3 * dx * dy * dy * dz * dz + fg_2 * dx * dz * dz;
                        fxyz_zz[g] += fg_3 * dx * dy * dz * dz * dz + 2.0 * fg_2 * dx * dy * dz;
                        fxzz_zz[g] += fg_3 * dx * dz * dz * dz * dz + 5.0 * fg_2 * dx * dz * dz + 2.0 * fg_1 * dx;
                        fyyy_zz[g] += fg_3 * dy * dy * dy * dz * dz + 3.0 * fg_2 * dy * dz * dz;
                        fyyz_zz[g] += fg_3 * dy * dy * dz * dz * dz + 2.0 * fg_2 * dy * dy * dz + fg_2 * dz * dz * dz + 2.0 * fg_1 * dz;
                        fyzz_zz[g] += fg_3 * dy * dz * dz * dz * dz + 5.0 * fg_2 * dy * dz * dz + 2.0 * fg_1 * dy;
                        fzzz_zz[g] += fg_3 * dz * dz * dz * dz * dz + 9.0 * fg_2 * dz * dz * dz + 12.0 * fg_1 * dz;
                    }
                }

                genfunc::transform(bspherbuff, bcartbuff, CSphericalMomentum(2), 0, 0, nGridPoints, 20);
            }
            else
            {
                // FIX ME: implement l > 2 cases

                std::string errangmom("gtoeval::computeGtosValuesForThirdOrder: Only implemented up to d-orbitals");

                errors::assertMsgCritical(false, errangmom);
            }

            // distribute j-th GTO values into grid values matrix

            for (int32_t k = 0; k < bnspher; k++)
            {
                auto bgaos = bspherbuff.data(20 * k);

                auto bgaox = bspherbuff.data(20 * k + 1);
                auto bgaoy = bspherbuff.data(20 * k + 2);
                auto bgaoz = bspherbuff.data(20 * k + 3);

                auto bgaoxx = bspherbuff.data(20 * k + 4);
                auto bgaoxy = bspherbuff.data(20 * k + 5);
                auto bgaoxz = bspherbuff.data(20 * k + 6);
                auto bgaoyy = bspherbuff.data(20 * k + 7);
                auto bgaoyz = bspherbuff.data(20 * k + 8);
                auto bgaozz = bspherbuff.data(20 * k + 9);

                auto bgaoxxx = bspherbuff.data(20 * k + 10);
                auto bgaoxxy = bspherbuff.data(20 * k + 11);
                auto bgaoxxz = bspherbuff.data(20 * k + 12);
                auto bgaoxyy = bspherbuff.data(20 * k + 13);
                auto bgaoxyz = bspherbuff.data(20 * k + 14);
                auto bgaoxzz = bspherbuff.data(20 * k + 15);
                auto bgaoyyy = bspherbuff.data(20 * k + 16);
                auto bgaoyyz = bspherbuff.data(20 * k + 17);
                auto bgaoyzz = bspherbuff.data(20 * k + 18);
                auto bgaozzz = bspherbuff.data(20 * k + 19);

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

                auto gvalxxx = gtoValuesXXX.data(idx);
                auto gvalxxy = gtoValuesXXY.data(idx);
                auto gvalxxz = gtoValuesXXZ.data(idx);
                auto gvalxyy = gtoValuesXYY.data(idx);
                auto gvalxyz = gtoValuesXYZ.data(idx);
                auto gvalxzz = gtoValuesXZZ.data(idx);
                auto gvalyyy = gtoValuesYYY.data(idx);
                auto gvalyyz = gtoValuesYYZ.data(idx);
                auto gvalyzz = gtoValuesYZZ.data(idx);
                auto gvalzzz = gtoValuesZZZ.data(idx);

                #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz, \
                    bgaoxx, bgaoxy, bgaoxz, bgaoyy, bgaoyz, bgaozz, \
                    gvalxx, gvalxy, gvalxz, gvalyy, gvalyz, gvalzz, \
                    bgaoxxx, bgaoxxy, bgaoxxz, bgaoxyy, bgaoxyz, bgaoxzz, bgaoyyy, bgaoyyz, bgaoyzz, bgaozzz, \
                    gvalxxx, gvalxxy, gvalxxz, gvalxyy, gvalxyz, gvalxzz, gvalyyy, gvalyyz, gvalyzz, gvalzzz : VLX_ALIGN)
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

                    gvalxxx[gridOffset + g] = bgaoxxx[g];
                    gvalxxy[gridOffset + g] = bgaoxxy[g];
                    gvalxxz[gridOffset + g] = bgaoxxz[g];
                    gvalxyy[gridOffset + g] = bgaoxyy[g];
                    gvalxyz[gridOffset + g] = bgaoxyz[g];
                    gvalxzz[gridOffset + g] = bgaoxzz[g];
                    gvalyyy[gridOffset + g] = bgaoyyy[g];
                    gvalyyz[gridOffset + g] = bgaoyyz[g];
                    gvalyzz[gridOffset + g] = bgaoyzz[g];
                    gvalzzz[gridOffset + g] = bgaozzz[g];
                }
            }
        }
    }
}

}  // namespace gtoeval
