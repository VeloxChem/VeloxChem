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

#include "PartitionFunc.hpp"

#include <algorithm>

#include "DenseMatrix.hpp"
#include "MathFunc.hpp"

namespace partfunc {  // partfunc namespace

auto
ssf(CDenseMatrix*   rawGridPoints,
    const double    minDistance,
    const int   gridOffset,
    const int   nGridPoints,
    const TPoint<double>*   atomCoordinates,
    const int   nAtoms,
    const int   idAtomic) -> void
{
    // partial weights

    std::vector<double> weights(nAtoms);

    auto pweights = weights.data();

    // set up pointers to grid data

    auto gridx = rawGridPoints->row(0) + gridOffset;

    auto gridy = rawGridPoints->row(1) + gridOffset;

    auto gridz = rawGridPoints->row(2) + gridOffset;

    auto gridw = rawGridPoints->row(3) + gridOffset;

    // loop over grid points

    for (int i = 0; i < nGridPoints; i++)
    {
        // grid coordinates

        double rgx = gridx[i];

        double rgy = gridy[i];

        double rgz = gridz[i];

        // weights screening

        auto atomxyz = atomCoordinates[idAtomic].coordinates();

        auto rig = mathfunc::distance(atomxyz[0], atomxyz[1], atomxyz[2], rgx, rgy, rgz);

        // min. distance scale 0.5 * (1 - a), SSF parameter a = 0.64

        if (rig < 0.18 * minDistance) continue;

        // initialize weights

        std::fill(pweights, pweights + nAtoms, 1.0);

        // outer loop over atoms

        for (int j = 0; j < nAtoms; j++)
        {
            // molecular coodinates

            auto raxyz = atomCoordinates[j].coordinates();

            // distance from grid point to j-th atom

            double rag = mathfunc::distance(raxyz[0], raxyz[1], raxyz[2], rgx, rgy, rgz);

            // loop over atoms

            for (int k = j + 1; k < nAtoms; k++)
            {
                // molecular coodinates

                auto rbxyz = atomCoordinates[k].coordinates();

                // distance from grid point to k-th atom

                double rbg = mathfunc::distance(rbxyz[0], rbxyz[1], rbxyz[2], rgx, rgy, rgz);

                // distance from j-th atom to k-th atom

                double rab = mathfunc::distance(raxyz[0], raxyz[1], raxyz[2], rbxyz[0], rbxyz[1], rbxyz[2]);

                // eliptical coordinate

                double mab = (rag - rbg) / rab;

                // scale partial weight

                pweights[j] *= 0.5 * (1.0 - partfunc::zeta(mab));

                pweights[k] *= 0.5 * (1.0 + partfunc::zeta(mab));
            }
        }

        //  adjust weight of i-th grid point

        double wsum = 0.0;

        for (int j = 0; j < nAtoms; j++)
        {
            wsum += pweights[j];
        }

        gridw[i] *= pweights[idAtomic] / wsum;
    }
}

inline auto
zeta(const double eRadius) -> double
{
    // SSF parameter a = 0.64

    // lower boundary

    if (eRadius <= -0.64) return -1.0;

    // upper boundary

    if (eRadius >= 0.64) return 1.0;

    // middle interval

    auto mab = 1.5625 * eRadius;

    auto mab2 = mab * mab;

    auto gab = 0.0625 * mab * (35.0 + mab2 * (-35.0 + mab2 * (21.0 - 5.0 * mab2)));

    return gab;
}

}  // namespace partfunc
