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
