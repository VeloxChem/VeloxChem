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

#include "AtomicRadii.hpp"

#include <cstring>

#include "Codata.hpp"

namespace atomicradii {  // atomicradii namespace

std::vector<double>
buildVdwRadii()
{
    // VDW radii in Angstrom

    // clang-format off
    std::vector<double> radii({
        0.00,                          // dummy atom
        1.09, 1.40, 1.82, 2.00, 2.00,  // H-B
        1.70, 1.55, 1.52, 1.47, 1.54,  // C-Ne
        2.27, 1.73, 2.00, 2.10, 1.80,  // Na-P
        1.80, 1.75, 1.88, 2.75, 2.00,  // S-Ca
        2.00, 2.00, 2.00, 2.00, 2.00,  // Sc-Mn
        2.00, 2.00, 1.63, 1.40, 1.39,  // Fe-Zn
        1.87, 2.00, 1.85, 1.90, 1.85,  // Ga-Br
        2.02, 2.00, 2.00, 2.00, 2.00,  // Kr-Zr
        2.00, 2.00, 2.00, 2.00, 2.00,  // Nb-Rh
        1.63, 1.72, 1.58, 1.93, 2.17,  // Pd-Sn
        2.00, 2.06, 1.98, 2.16, 2.00,  // Sb-Cs
        2.00, 2.00, 2.00, 2.00, 2.00,  // Ba-Nd
        2.00, 2.00, 2.00, 2.00, 2.00,  // Pm-Tb
        2.00, 2.00, 2.00, 2.00, 2.00,  // Dy-Yb
        2.00, 2.00, 2.00, 2.00, 2.00,  // Lu-Re
        2.00, 2.00, 1.72, 1.66, 1.55,  // Os-Hg
        1.96, 2.02, 2.00, 2.00, 2.00,  // Tl-At
        2.00                           // Rn
    });
    // clang-format on

    // VDW radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::getBohrValueInAngstroms();
    }

    return radii;
}

std::vector<double>
buildMkRadii()
{
    // MK radii in Angstrom
    // (from: A. Alenazian, L. A. Burns, C. D. Sherrill, Int. J. Quantum Chem. 2020, 120, e26035)

    // clang-format off
    std::vector<double> radii({
        0.00,                          // dummy atom
        1.20, 1.20, 1.37, 1.45, 1.45,  // H-B
        1.50, 1.50, 1.40, 1.35, 1.30,  // C-Ne
        1.57, 1.36, 1.24, 1.17, 1.80,  // Na-P
        1.75, 1.70                     // S-Cl
    });
    // clang-format on

    // MK radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::getBohrValueInAngstroms();
    }

    return radii;
}

std::vector<double>
buildChelpgRadii()
{
    // CHELPG radii in Angstrom
    // C. M. Breneman, K. B. Wiberg, J. Comput. Chem. 1990, 11, 361-373

    // clang-format off
    std::vector<double> radii({
        0.00,                          // dummy atom
        1.45, 1.45, 1.50, 1.50, 1.50,  // H-B
        1.50, 1.70, 1.70, 1.70, 1.70,  // C-Ne
    });
    // clang-format on

    // CHELPG radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::getBohrValueInAngstroms();
    }

    return radii;
}

std::vector<double>
buildCovalentRadii()
{
    // covalent radii in Angstrom
    // (from: ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx)

    // clang-format off
    std::vector<double> radii({
        0.00,                          // dummy atom
        0.23, 1.50, 1.28, 0.96, 0.83,  // H-B
        0.68, 0.68, 0.68, 0.64, 1.50,  // C-Ne
        1.66, 1.41, 1.21, 1.20, 1.05,  // Na-P
        1.02, 0.99, 1.51, 2.03, 1.76,  // S-Ca
        1.70, 1.60, 1.53, 1.39, 1.61,  // Sc-Mn
        1.52, 1.26, 1.24, 1.32, 1.22,  // Fe-Zn
        1.22, 1.17, 1.21, 1.22, 1.21,  // Ga-Br
        1.50, 2.20, 1.95, 1.90, 1.75,  // Kr-Zr
        1.64, 1.54, 1.47, 1.46, 1.42,  // Nb-Rh
        1.39, 1.45, 1.54, 1.42, 1.39,  // Pd-Sn
        1.39, 1.47, 1.40, 1.50, 2.44,  // Sb-Cs
        2.15, 2.07, 2.04, 2.03, 2.01,  // Ba-Nd
        1.99, 1.98, 1.98, 1.96, 1.94,  // Pm-Tb
        1.92, 1.92, 1.89, 1.90, 1.87,  // Dy-Yb
        1.87, 1.75, 1.70, 1.62, 1.51,  // Lu-Re
        1.44, 1.41, 1.36, 1.36, 1.32,  // Os-Hg
        1.45, 1.46, 1.48, 1.40, 1.21,  // Tl-At
        1.50                           // Rn
    });
    // clang-format on

    // covalent radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::getBohrValueInAngstroms();
    }

    return radii;
}

}  // namespace atomicradii
