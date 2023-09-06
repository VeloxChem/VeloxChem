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

#include "PartitionFunc.hpp"

#include <algorithm>

#include "DenseMatrix.hpp"
#include "MathFunc.hpp"

namespace partfunc {  // partfunc namespace

auto
ssf(CDenseMatrix*   rawGridPoints,
    const double    minDistance,
    const int64_t   gridOffset,
    const int64_t   nGridPoints,
    const TPoint3D* atomCoordinates,
    const int64_t   nAtoms,
    const int64_t   idAtomic) -> void
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

    for (int64_t i = 0; i < nGridPoints; i++)
    {
        // grid coordinates

        double rgx = gridx[i];

        double rgy = gridy[i];

        double rgz = gridz[i];

        // weights screening

        auto rig = mathfunc::distance(atomCoordinates[idAtomic][0], atomCoordinates[idAtomic][1], atomCoordinates[idAtomic][2], rgx, rgy, rgz);

        // min. distance scale 0.5 * (1 - a), SSF parameter a = 0.64

        if (rig < 0.18 * minDistance) continue;

        // initialize weights

        std::fill(pweights, pweights + nAtoms, 1.0);

        // outer loop over atoms

        for (int64_t j = 0; j < nAtoms; j++)
        {
            // molecular coodinates

            double rax = atomCoordinates[j][0];

            double ray = atomCoordinates[j][1];

            double raz = atomCoordinates[j][2];

            // distance from grid point to j-th atom

            double rag = mathfunc::distance(rax, ray, raz, rgx, rgy, rgz);

            // loop over atoms

            for (int64_t k = j + 1; k < nAtoms; k++)
            {
                // molecular coodinates

                double rbx = atomCoordinates[k][0];

                double rby = atomCoordinates[k][1];

                double rbz = atomCoordinates[k][2];

                // distance from grid point to k-th atom

                double rbg = mathfunc::distance(rbx, rby, rbz, rgx, rgy, rgz);

                // distance from j-th atom to k-th atom

                double rab = mathfunc::distance(rax, ray, raz, rbx, rby, rbz);

                // eliptical coordinate

                double mab = (rag - rbg) / rab;

                // scale partial weight

                pweights[j] *= 0.5 * (1.0 - partfunc::zeta(mab));

                pweights[k] *= 0.5 * (1.0 + partfunc::zeta(mab));
            }
        }

        //  adjust weight of i-th grid point

        double wsum = 0.0;

        for (int64_t j = 0; j < nAtoms; j++)
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
