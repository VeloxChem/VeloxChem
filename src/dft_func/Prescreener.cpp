//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Prescreener.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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

        const auto gto_xyz = gto_coords[i].coordinates();

        double rx = std::max({xmin - gto_xyz[0], gto_xyz[0] - xmax, 0.0});
        double ry = std::max({ymin - gto_xyz[1], gto_xyz[1] - ymax, 0.0});
        double rz = std::max({zmin - gto_xyz[2], gto_xyz[2] - zmax, 0.0});

        auto r2 = rx * rx + ry * ry + rz * rz;

        if (r2 > 1.0)
        {
            // Sum the bounds of the primitives, each evaluated with its own
            // exponent and coefficient. By the triangle inequality the sum
            // bounds the contraction and its derivatives.

            auto r = std::sqrt(r2);

            auto rpow = 1.0;

            for (int ipow = 0; ipow < gto_ang; ipow++)
            {
                rpow *= r;
            }

            double gtolimit = 0.0;

            auto minexp = std::numeric_limits<double>::max();

            for (int j = 0; j < npgtos; j++)
            {
                int iprim = j * ncgtos + i;

                const auto prim_exp = gto_exps[iprim];

                const auto prim_norm = std::fabs(gto_norms[iprim]);

                minexp = std::min(minexp, prim_exp);

                const auto base = prim_norm * std::exp(-prim_exp * r2) * rpow;

                auto limit = base;

                if (gtoDeriv > 0)
                {
                    limit = std::max(limit, 2.0 * prim_exp * r * base);

                    if (gto_ang > 0) limit = std::max(limit, base / r * gto_ang);
                }

                if (gtoDeriv > 1)
                {
                    limit = std::max({limit,
                                      2.0 * prim_exp * (2 * gto_ang + 1) * base,
                                      4.0 * prim_exp * prim_exp * r2 * base});

                    if (gto_ang > 1) limit = std::max(limit, base / r2 * gto_ang * (gto_ang - 1));
                }

                if (gtoDeriv > 2)
                {
                    const auto r3 = r * r2;

                    limit = std::max({limit,
                                      4.0 * prim_exp * prim_exp * (3 * gto_ang + 3) * r * base,
                                      8.0 * prim_exp * prim_exp * prim_exp * r3 * base});

                    if (gto_ang > 0)
                    {
                        limit = std::max(limit, 2.0 * prim_exp * base / r * (3 * gto_ang * gto_ang));
                    }

                    if (gto_ang > 2)
                    {
                        limit = std::max(limit, base / r3 * gto_ang * (gto_ang - 1) * (gto_ang - 2));
                    }
                }

                gtolimit += limit;
            }

            // The bound is taken at the closest point, so screen only past the
            // radial peak sqrt((gto_ang + gtoDeriv) / (2 * minexp)).

            const auto past_peak = ((2.0 * minexp * r2) >= static_cast<double>(gto_ang + gtoDeriv));

            if ((gtolimit < gtoThreshold) && past_peak)
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
