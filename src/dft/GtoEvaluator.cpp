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

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

namespace gtoeval {  // gtoeval namespace

void
computeGtosValuesForGGA(CMemBlock2D<double>&         gtoValues,
                        CMemBlock2D<double>&         gtoValuesX,
                        CMemBlock2D<double>&         gtoValuesY,
                        CMemBlock2D<double>&         gtoValuesZ,
                        const CGtoContainer*         gtoContainer,
                        const double*                gridCoordinatesX,
                        const double*                gridCoordinatesY,
                        const double*                gridCoordinatesZ,
                        const int32_t                gridBlockPosition,
                        const int32_t                gridOffset,
                        const int32_t                nGridPoints,
                        const std::array<double, 6>& boxDimension,
                        const double                 gtoThreshold)
{
    // spatial extent of grid box

    auto xmin = boxDimension[0];

    auto ymin = boxDimension[1];

    auto zmin = boxDimension[2];

    auto xmax = boxDimension[3];

    auto ymax = boxDimension[4];

    auto zmax = boxDimension[5];

    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);

    // loop over GTOs container data

    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);

        // angular momentum data for bra and ket

        auto bang = bgtos.getAngularMomentum();

        // set up Cartesian GTOs buffer

        auto nvcomp = xcfun_components(xcfun::gga);

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

        // contracted GTO screening

        double rx = 0.0, ry = 0.0, rz = 0.0;

        if      (bfx[spos[0]] < xmin) rx = xmin - bfx[spos[0]];

        else if (bfx[spos[0]] > xmax) rx = bfx[spos[0]] - xmax;

        if      (bfy[spos[0]] < ymin) ry = ymin - bfy[spos[0]];

        else if (bfy[spos[0]] > ymax) ry = bfy[spos[0]] - ymax;

        if      (bfz[spos[0]] < zmin) rz = zmin - bfz[spos[0]];

        else if (bfz[spos[0]] > zmax) rz = bfz[spos[0]] - zmax;

        auto r2 = rx * rx + ry * ry + rz * rz;

        if (r2 > 1.0)
        {
            auto minexp = bfexps[spos[0]];

            auto maxexp = bfexps[spos[0]];

            auto maxcoef = std::fabs(bfnorms[spos[0]]);

            for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
            {
                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = std::fabs(bfnorms[iprim]);

                    if (minexp > bexp) minexp = bexp;

                    if (maxexp < bexp) maxexp = bexp;

                    if (maxcoef < bnorm) maxcoef = bnorm;
                }
            }

            // gto  :           r^{ang}   |C| exp(-alpha r^2)
            // gto_m:           r^{ang-1} |C| exp(-alpha r^2)
            // gto_p: (2 alpha) r^{ang+1} |C| exp(-alpha r^2)

            // Note that gto_m < gto (r > 1)

            auto r = std::sqrt(r2);

            auto gtolimit = maxcoef * std::exp(-minexp * r2);

            for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit *= r;

            auto gtolimit_p = 2.0 * maxexp * r * gtolimit;

            if ((gtolimit < gtoThreshold) && (gtolimit_p < gtoThreshold)) continue;
        }

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
        {
            auto rx = bfx[spos[j]];

            auto ry = bfy[spos[j]];

            auto rz = bfz[spos[j]];

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

                auto idx = (bgtos.getIdentifiers(k))[j];

                auto gvals = gtoValues.data(idx);

                auto gvalx = gtoValuesX.data(idx);

                auto gvaly = gtoValuesY.data(idx);

                auto gvalz = gtoValuesZ.data(idx);

                #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz : VLX_ALIGN)
                for (int32_t g = 0; g < nGridPoints; g++)
                {
                    gvals[g] = bgaos[g];

                    gvalx[g] = bgaox[g];

                    gvaly[g] = bgaoy[g];

                    gvalz[g] = bgaoz[g];
                }
            }
        }
    }
}

}  // namespace gtoeval
