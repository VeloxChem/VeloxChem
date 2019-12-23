//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CoordinationNumber.hpp"

#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "Molecule.hpp"

namespace coordnum {  // coordnum namespace

std::vector<double>
getCovalentRadius()
{
    // Reference: dftd4 (v2.4.0)

    std::vector<double> cn({
        0.0,  0.32, 0.46,                                // H,He
        1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,  // Li-Ne
        1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,  // Na-Ar
        1.76, 1.54,                                      // K,Ca
        1.33, 1.22, 1.21, 1.10, 1.07,                    // Sc-
        1.04, 1.00, 0.99, 1.01, 1.09,                    // -Zn
        1.12, 1.09, 1.15, 1.10, 1.14, 1.17,              // Ga-Kr
        1.89, 1.67,                                      // Rb,Sr
        1.47, 1.39, 1.32, 1.24, 1.15,                    // Y-
        1.13, 1.13, 1.08, 1.15, 1.23,                    // -Cd
        1.28, 1.26, 1.26, 1.23, 1.32, 1.31,              // In-Xe
        2.09, 1.76,                                      // Cs,Ba
        1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,        // La-Eu
        1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,        // Gd-Yb
        1.46, 1.37, 1.31, 1.23, 1.18,                    // Lu-
        1.16, 1.11, 1.12, 1.13, 1.32,                    // -Hg
        1.30, 1.30, 1.36, 1.31, 1.38, 1.42,              // Tl-Rn
        2.01, 1.81,                                      // Fr,Ra
        1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,        // Ac-Am
        1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,        // Cm-No
        1.45, 1.41, 1.34, 1.29, 1.27,                    // Lr-
        1.21, 1.16, 1.15, 1.09, 1.22,                    // -Cn
        1.36, 1.43, 1.46, 1.58, 1.48, 1.57,              // Nh-Og
    });

    for (size_t i = 0; i < cn.size(); i++)
    {
        cn[i] /= units::getBohrValueInAngstroms();
    }

    return cn;
}

std::vector<double>
getCoordinationNumber(const CMolecule& molecule)
{
    CDenseMatrix dcndr;

    return getCoordinationNumber(molecule, dcndr);
}

std::vector<double>
getCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcndr)
{
    // Reference: dftd4 (v2.4.0)

    // prepare parameters

    const double k2 = 4.0 / 3.0;

    const double kn = 7.50;

    const double cn_thr = 1600.0;

    const double cnmax = 8.0;

    const double sqrtpi = std::sqrt(mathconst::getPiValue());

    auto covalent_radius = getCovalentRadius();

    // get molecular information

    auto natoms = molecule.getNumberOfAtoms();

    auto ids_elem = molecule.getIdsElemental();

    auto xcoord = molecule.getCoordinatesX();

    auto ycoord = molecule.getCoordinatesY();

    auto zcoord = molecule.getCoordinatesZ();

    // compute coordination numbers with error function

    std::vector<double> cn(natoms, 0.0);

    if (dcndr.getNumberOfElements() > 0)
    {
        std::string err_size("CoordinationNumber - Mismatch in dcndr matrix size");

        errors::assertMsgCritical(dcndr.getNumberOfRows() == 3 * natoms, err_size);

        errors::assertMsgCritical(dcndr.getNumberOfColumns() == natoms, err_size);
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < i; j++)
        {
            std::vector<double> rij(3);

            rij[0] = xcoord[j] - xcoord[i];

            rij[1] = ycoord[j] - ycoord[i];

            rij[2] = zcoord[j] - zcoord[i];

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            if (r2 > cn_thr) continue;

            double r = std::sqrt(r2);

            double rco = k2 * (covalent_radius[ids_elem[j]] + covalent_radius[ids_elem[i]]);

            double cn_val = 0.5 * (1.0 + std::erf(-kn * (r - rco) / rco));

            cn[i] += cn_val;

            cn[j] += cn_val;

            double arg = kn * (r - rco) / rco;

            double dcn_val = -kn / sqrtpi / rco * std::exp(-arg * arg);

            if (dcndr.getNumberOfElements() > 0)
            {
                for (int32_t d = 0; d < 3; d++)
                {
                    dcndr.values()[(d * natoms + i) * natoms + i] += dcn_val * rij[d] / r;

                    dcndr.values()[(d * natoms + j) * natoms + j] -= dcn_val * rij[d] / r;

                    dcndr.values()[(d * natoms + i) * natoms + j] += dcn_val * rij[d] / r;

                    dcndr.values()[(d * natoms + j) * natoms + i] -= dcn_val * rij[d] / r;
                }
            }
        }
    }

    // apply cutoff function for large coordination numbers

    for (int32_t i = 0; i < natoms; i++)
    {
        cn[i] = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn[i]));

        double dcnpdcn = std::exp(cnmax) / (std::exp(cnmax) + std::exp(cn[i]));

        if (dcndr.getNumberOfElements() > 0)
        {
            for (int32_t d = 0; d < 3; d++)
            {
                for (int32_t j = 0; j < natoms; j++)
                {
                    dcndr.values()[(d * natoms + j) * natoms + i] *= dcnpdcn;
                }
            }
        }
    }

    return cn;
}

}  // namespace coordnum
