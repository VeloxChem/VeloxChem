//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#include "CoordinationNumber.hpp"

#include <cmath>

#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "Molecule.hpp"

namespace coordnum {  // coordnum namespace

auto
getCovalentRadius() -> std::vector<double>
{
    // Reference: dftd4 (v2.4.0)

    std::vector<double> cn({
        0.0,                                             // dummy atom
        0.32, 0.46,                                      // H,He
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
        1.36, 1.43, 1.46, 1.58, 1.48, 1.57               // Nh-Og
    });

    for (size_t i = 0; i < cn.size(); i++)
    {
        // use dftd4 conversion factor
        cn[i] /= 0.52917726;
    }

    return cn;
}

auto
getCoordinationNumber(const CMolecule& molecule) -> std::vector<double>
{
    CDenseMatrix dcndr;

    return getCoordinationNumber(molecule, dcndr);
}

auto
getCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcndr) -> std::vector<double>
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

    auto xyzcoord = molecule.getCoordinates("bohr");

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
            std::vector<double> rij({xyzcoord[j][0] - xyzcoord[i][0], xyzcoord[j][1] - xyzcoord[i][1], xyzcoord[j][2] - xyzcoord[i][2]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            if (r2 > cn_thr) continue;

            double r = std::sqrt(r2);

            double rco = k2 * (covalent_radius[ids_elem[j]] + covalent_radius[ids_elem[i]]);

            double arg = kn * (r - rco) / rco;

            double cn_val = 0.5 * (1.0 + std::erf(-arg));

            cn[i] += cn_val;

            cn[j] += cn_val;

            double dcn_val = -kn / sqrtpi / rco * std::exp(-arg * arg);

            if (dcndr.getNumberOfElements() > 0)
            {
                for (int32_t d = 0; d < 3; d++)
                {
                    int32_t di = d * natoms + i;

                    int32_t dj = d * natoms + j;

                    dcndr.values()[di * natoms + i] += dcn_val * rij[d] / r;

                    dcndr.values()[dj * natoms + j] -= dcn_val * rij[d] / r;

                    dcndr.values()[di * natoms + j] += dcn_val * rij[d] / r;

                    dcndr.values()[dj * natoms + i] -= dcn_val * rij[d] / r;
                }
            }
        }
    }

    // apply cutoff function for large coordination numbers

    for (int32_t i = 0; i < natoms; i++)
    {
        double dcnpdcn = std::exp(cnmax) / (std::exp(cnmax) + std::exp(cn[i]));

        if (dcndr.getNumberOfElements() > 0)
        {
            for (int32_t dj = 0; dj < 3 * natoms; dj++)
            {
                dcndr.values()[dj * natoms + i] *= dcnpdcn;
            }
        }
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        cn[i] = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn[i]));
    }

    return cn;
}

auto
getPaulingElectronegativity() -> std::vector<double>
{
    // Reference: dftd4 (v2.4.0)

    return std::vector<double>({
        0.0,                                                         // dummy atom
        2.20, 3.00,                                                  // H,He
        0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50,              // Li-Ne
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,              // Na-Ar
        0.82, 1.00,                                                  // K,Ca
        1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  // Sc-Zn
        1.81, 2.01, 2.18, 2.55, 2.96, 3.00,                          // Ga-Kr
        0.82, 0.95,                                                  // Rb,Sr
        1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  // Y-Cd
        1.78, 1.96, 2.05, 2.10, 2.66, 2.60,                          // In-Xe
        0.79, 0.89,                                                  // Cs,Ba
        1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18,                    // La-Eu
        1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26,                    // Gd-Yb
        1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,  // Lu-Hg
        1.62, 2.33, 2.02, 2.00, 2.20, 2.20,                          // Tl-Rn
        1.50, 1.50,                                                  // Fr,Ra
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    // Ac-Am
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    // Cm-No
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,  // Rf-Cn
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50                           // Nh-Og
    });
}

auto
getCovalentCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcovcndr) -> std::vector<double>
{
    // Reference: dftd4 (v2.4.0)

    // prepare parameters

    const double k2 = 4.0 / 3.0;

    const double k4 = 4.10451;

    const double k5 = 19.08857;

    const double k6 = 2.0 * std::pow(11.28174, 2);

    const double kn = 7.50;

    const double cn_thr = 1600.0;

    const double sqrtpi = std::sqrt(mathconst::getPiValue());

    auto covalent_radius = getCovalentRadius();

    auto pauling_en = getPaulingElectronegativity();

    // get molecular information

    auto natoms = molecule.getNumberOfAtoms();

    auto ids_elem = molecule.getIdsElemental();

    auto xyzcoord = molecule.getCoordinates("bohr");

    // compute covalent coordination numbers

    std::vector<double> covcn(natoms, 0.0);

    if (dcovcndr.getNumberOfElements() > 0)
    {
        std::string err_size("CovalentCoordinationNumber - Mismatch in dcovcndr matrix size");

        errors::assertMsgCritical(dcovcndr.getNumberOfRows() == 3 * natoms, err_size);

        errors::assertMsgCritical(dcovcndr.getNumberOfColumns() == natoms, err_size);
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < i; j++)
        {
            std::vector<double> rij({xyzcoord[j][0] - xyzcoord[i][0], xyzcoord[j][1] - xyzcoord[i][1], xyzcoord[j][2] - xyzcoord[i][2]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            if (r2 > cn_thr) continue;

            double r = std::sqrt(r2);

            double rco = k2 * (covalent_radius[ids_elem[j]] + covalent_radius[ids_elem[i]]);

            double diff_en = std::fabs(pauling_en[ids_elem[i]] - pauling_en[ids_elem[j]]);

            double den = k4 * std::exp(-std::pow(diff_en + k5, 2) / k6);

            double arg = kn * (r - rco) / rco;

            double covcn_val = den * 0.5 * (1.0 + std::erf(-arg));

            covcn[i] += covcn_val;

            covcn[j] += covcn_val;

            double dcovcn_val = -den * kn / sqrtpi / rco * std::exp(-arg * arg);

            if (dcovcndr.getNumberOfElements() > 0)
            {
                for (int32_t d = 0; d < 3; d++)
                {
                    int32_t di = d * natoms + i;

                    int32_t dj = d * natoms + j;

                    dcovcndr.values()[di * natoms + i] -= dcovcn_val * rij[d] / r;

                    dcovcndr.values()[dj * natoms + j] += dcovcn_val * rij[d] / r;

                    dcovcndr.values()[di * natoms + j] -= dcovcn_val * rij[d] / r;

                    dcovcndr.values()[dj * natoms + i] += dcovcn_val * rij[d] / r;
                }
            }
        }
    }

    return covcn;
}

}  // namespace coordnum
