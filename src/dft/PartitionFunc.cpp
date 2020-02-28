//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "PartitionFunc.hpp"

#include "MathFunc.hpp"

namespace partfunc {  // partfunc namespace

void
ssf(CMemBlock2D<double>* rawGridPoints,
    const double         minDistance,
    const int32_t        gridOffset,
    const int32_t        nGridPoints,
    const double*        atomCoordinatesX,
    const double*        atomCoordinatesY,
    const double*        atomCoordinatesZ,
    const int32_t        nAtoms,
    const int32_t        idAtomic)
{
    // partial weights

    CMemBlock<double> weights(nAtoms);

    auto pweights = weights.data();

    // set up pointers to grid data

    auto gridx = rawGridPoints->data(0, gridOffset);

    auto gridy = rawGridPoints->data(1, gridOffset);

    auto gridz = rawGridPoints->data(2, gridOffset);

    auto gridw = rawGridPoints->data(3, gridOffset);

    // loop over grid points

    for (int32_t i = 0; i < nGridPoints; i++)
    {
        // grid coordinates

        double rgx = gridx[i];

        double rgy = gridy[i];

        double rgz = gridz[i];

        // weights screening

        auto rig = mathfunc::distance(atomCoordinatesX[idAtomic], atomCoordinatesY[idAtomic], atomCoordinatesZ[idAtomic], rgx, rgy, rgz);

        // min. distance scale 0.5 * (1 - a), SSF parameter a = 0.64

        if (rig < 0.18 * minDistance) continue;

        // initialize weights

        mathfunc::set_to(pweights, 1.0, nAtoms);

        // outer loop over atoms

        for (int32_t j = 0; j < nAtoms; j++)
        {
            // molecular coodinates

            double rax = atomCoordinatesX[j];

            double ray = atomCoordinatesY[j];

            double raz = atomCoordinatesZ[j];

            // distance from grid point to j-th atom

            double rag = mathfunc::distance(rax, ray, raz, rgx, rgy, rgz);

            // loop over atoms

            for (int32_t k = j + 1; k < nAtoms; k++)
            {
                // molecular coodinates

                double rbx = atomCoordinatesX[k];

                double rby = atomCoordinatesY[k];

                double rbz = atomCoordinatesZ[k];

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

        gridw[i] *= pweights[idAtomic] / mathfunc::sum(pweights, nAtoms);
    }
}

inline double
zeta(const double eRadius)
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
