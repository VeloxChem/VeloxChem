//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "Prescreener.hpp"

#include <algorithm>
#include <cmath>

#include "GtoFunc.hpp"
#include "GtoValuesRecD.hpp"
#include "GtoValuesRecP.hpp"
#include "GtoValuesRecS.hpp"

namespace prescr {  // prescr namespace

auto
getGridBoxDimension(const int gridBlockPosition, const int nGridPoints, const double* xcoords, const double* ycoords, const double* zcoords)
    -> std::array<double, 6>
{
    double xmin = xcoords[gridBlockPosition], ymin = ycoords[gridBlockPosition], zmin = zcoords[gridBlockPosition];

    double xmax = xcoords[gridBlockPosition], ymax = ycoords[gridBlockPosition], zmax = zcoords[gridBlockPosition];

    for (int g = 0; g < nGridPoints; g++)
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

auto
preScreenGtoBlock(const CGtoBlock& gtoBlock, const int gtoDeriv, const double gtoThreshold, const std::array<double, 6>& boxDimension)
    -> std::tuple<std::vector<int>, std::vector<int>>
{
    // grid box dimension

    double xmin = boxDimension[0], ymin = boxDimension[1], zmin = boxDimension[2];

    double xmax = boxDimension[3], ymax = boxDimension[4], zmax = boxDimension[5];

    // GTOs data

    const auto gto_exps = gtoBlock.exponents();

    const auto gto_norms = gtoBlock.normalization_factors();

    const auto gto_coords = gtoBlock.coordinates();

    // GTOs block dimensions

    const auto ncgtos = gtoBlock.number_of_basis_functions();

    const auto npgtos = gtoBlock.number_of_primitives();

    // pre-screen GTOs for grid box

    auto gto_ang = gtoBlock.angular_momentum();

    auto ncomps = gto_ang * 2 + 1;

    std::vector<int> cgto_mask(ncgtos, 1);

    std::vector<int> ao_mask(ncomps * ncgtos, 1);

    for (int i = 0; i < ncgtos; i++)
    {
        // GTO coordinates

        const auto gto_x = gto_coords[i].coordinates()[0];
        const auto gto_y = gto_coords[i].coordinates()[1];
        const auto gto_z = gto_coords[i].coordinates()[2];

        double rx = std::max({xmin - gto_x, gto_x - xmax, 0.0});
        double ry = std::max({ymin - gto_y, gto_y - ymax, 0.0});
        double rz = std::max({zmin - gto_z, gto_z - zmax, 0.0});

        auto r2 = rx * rx + ry * ry + rz * rz;

        if (r2 > 1.0)
        {
            int firstprim = 0 * ncgtos + i;

            auto minexp = gto_exps[firstprim];

            auto maxexp = gto_exps[firstprim];

            auto maxcoef = std::fabs(gto_norms[firstprim]);

            for (int j = 0; j < npgtos; j++)
            {
                int iprim = j * ncgtos + i;

                auto prim_exp = gto_exps[iprim];

                auto prim_norm = std::fabs(gto_norms[iprim]);

                minexp = std::min(minexp, prim_exp);

                maxexp = std::max(maxexp, prim_exp);

                maxcoef = std::max(maxcoef, prim_norm);
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

            for (int ipow = 0; ipow < gto_ang; ipow++)
            {
                gtolimit_base *= r;
            }

            auto gtolimit = gtolimit_base;

            if (gtoDeriv > 0)
            {
                auto gtolimit_p = 2.0 * maxexp * r * gtolimit_base;  // gto_p

                gtolimit = std::max(gtolimit, gtolimit_p);

                if (gto_ang > 0)
                {
                    auto gtolimit_m = gtolimit_base / r * gto_ang;  // gto_m

                    gtolimit = std::max(gtolimit, gtolimit_m);
                }
            }

            if (gtoDeriv > 1)
            {
                auto gtolimit_p2 = 4.0 * maxexp * maxexp * r2 * gtolimit_base;  // gto_p2

                auto gtolimit_0 = 2.0 * maxexp * (2 * gto_ang + 1) * gtolimit_base;  // gto

                gtolimit = std::max({gtolimit, gtolimit_0, gtolimit_p2});

                if (gto_ang > 1)
                {
                    auto gtolimit_m2 = gtolimit_base / r2 * gto_ang * (gto_ang - 1);  // gto_m2

                    gtolimit = std::max(gtolimit, gtolimit_m2);
                }
            }

            if (gtoDeriv > 2)
            {
                auto r3 = r * r2;

                auto gtolimit_p3 = 8.0 * maxexp * maxexp * maxexp * r3 * gtolimit_base;  // gto_p3

                auto gtolimit_p1 = 4.0 * maxexp * maxexp * (3 * gto_ang + 3) * r * gtolimit_base;  // gto_p1

                gtolimit = std::max({gtolimit, gtolimit_p1, gtolimit_p3});

                if (gto_ang > 0)
                {
                    auto gtolimit_m1 = 2.0 * maxexp * gtolimit_base / r * (3 * gto_ang * gto_ang);  // gto_m1

                    gtolimit = std::max(gtolimit, gtolimit_m1);
                }

                if (gto_ang > 2)
                {
                    auto gtolimit_m3 = gtolimit_base / r3 * gto_ang * (gto_ang - 1) * (gto_ang - 2);  // gto_m3

                    gtolimit = std::max(gtolimit, gtolimit_m3);
                }
            }

            if (gtolimit < gtoThreshold)
            {
                cgto_mask[i] = 0;

                for (int s = 0; s < ncomps; s++)
                {
                    ao_mask[s * ncgtos + i] = 0;
                }
            }
        }
    }

    // get pre-screened AO indices

    std::vector<int> pre_ao_inds;

    auto gto_ao_inds = gtoBlock.getAtomicOrbitalsIndexes();

    for (size_t i = 0; i < gto_ao_inds.size(); i++)
    {
        if (ao_mask[i] == 1) pre_ao_inds.push_back(gto_ao_inds[i]);
    }

    return {cgto_mask, pre_ao_inds};
}

}  // namespace prescr
