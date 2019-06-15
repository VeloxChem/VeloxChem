//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "Codata.hpp"

namespace units {  // units namespace

double
getBohrValueInAngstroms()
{
    return 0.52917721092;
}

double
getHartreeValueInElectronVolts()
{
    return 27.21138505;
}

double
getRotatoryStrengthInCGS()
{
    // Derived from Bohr radius (a0), elementary charge (e),
    // Bohr magneton (mu_B) and speed of light in vacuum (c).

    // 1 [a.u.] = 2 [e a0 mu_B] = 471.44364819 [10**(-40) cgs unit]

    return 471.44364819;
}

}  // namespace units
