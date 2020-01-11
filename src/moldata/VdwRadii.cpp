//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "VdwRadii.hpp"

#include <cstring>

#include "Codata.hpp"

namespace vdwradii {  // vdwradii namespace

std::vector<double>
buildVdwRadii()
{
    // VDW radii in Angstrom

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

    // VDW radii in Bohr

    for (size_t i = 0; i < radii.size(); i++)
    {
        radii[i] /= units::getBohrValueInAngstroms();
    }

    return radii;
}

}  // namespace vdwradii
