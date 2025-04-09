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

#include "AtomicRadii.hpp"

#include <cstring>

#include "Codata.hpp"

namespace atomicradii {  // atomicradii namespace

std::vector<double>
buildVdwRadii()
{
    // VDW radii in Angstrom

    // Hydrogen:       J. Phys. Chem.   1996, 100, 18, 7384-7391
    // Other elements: J. Phys. Chem.   1964,  68,  3,  441- 451
    //                 J. Phys. Chem. A 2009, 113, 19, 5806-5812

    // clang-format off
    std::vector<double> radii({
        0.00,                          // dummy atom
        1.10, 1.40, 1.81, 1.53, 1.92,  // H-B
        1.70, 1.55, 1.52, 1.47, 1.54,  // C-Ne
        2.27, 1.73, 1.84, 2.10, 1.80,  // Na-P
        1.80, 1.75, 1.88, 2.75, 2.31,  // S-Ca
        2.00, 2.00, 2.00, 2.00, 2.00,  // Sc-Mn
        2.00, 2.00, 1.63, 1.40, 1.39,  // Fe-Zn
        1.87, 2.11, 1.85, 1.90, 1.83,  // Ga-Br
        2.02, 3.03, 2.49, 2.00, 2.00,  // Kr-Zr
        2.00, 2.00, 2.00, 2.00, 2.00,  // Nb-Rh
        1.63, 1.72, 1.58, 1.93, 2.17,  // Pd-Sn
        2.06, 2.06, 1.98, 2.16, 3.43,  // Sb-Cs
        2.68, 2.00, 2.00, 2.00, 2.00,  // Ba-Nd
        2.00, 2.00, 2.00, 2.00, 2.00,  // Pm-Tb
        2.00, 2.00, 2.00, 2.00, 2.00,  // Dy-Yb
        2.00, 2.00, 2.00, 2.00, 2.00,  // Lu-Re
        2.00, 2.00, 1.72, 1.66, 1.55,  // Os-Hg
        1.96, 2.02, 2.07, 1.97, 2.02,  // Tl-At
        2.20                           // Rn
    });
    // clang-format on

    // VDW radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::bohr_in_angstrom();
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
        radii[i] /= units::bohr_in_angstrom();
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
        radii[i] /= units::bohr_in_angstrom();
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
        radii[i] /= units::bohr_in_angstrom();
    }

    return radii;
}

}  // namespace atomicradii
