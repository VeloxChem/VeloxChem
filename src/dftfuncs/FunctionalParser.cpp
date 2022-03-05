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

#include "FunctionalParser.hpp"

#include <algorithm>

#include "Becke88Functional.hpp"
#include "LYPFunctional.hpp"
#include "SlaterFunctional.hpp"
#include "StringFormat.hpp"
#include "VWN3Functional.hpp"
#include "PkzbFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace

std::vector<std::string>
getAvailableFunctionals()
{
    return std::vector<std::string>({"SLATER", "VWN3", "BECKE88", "LYP", "SLDA", "B88X", "BLYP", "B3LYP", "BHANDH", "BHANDHLYP","PKZB"});
}

CXCFunctional
getExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailableFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // pure spin-polarized Slater exchange functional

        if (fstr::upcase(xcLabel) == "SLATER") return vxcfuncs::setSlaterFunctional();

        // pure spin-polarized Vosko-Wilk-Nusair (Parameterization 3) correlation functional

        if (fstr::upcase(xcLabel) == "VWN3") return vxcfuncs::setVWN3Functional();

        // pure spin-polarized Becke (1988) exchange functional

        if (fstr::upcase(xcLabel) == "BECKE88") return vxcfuncs::setBecke88Functional();

        // pure spin-polarized Lee, Yang and Parr correlation functional

        if (fstr::upcase(xcLabel) == "LYP") return vxcfuncs::setLYPFunctional();

        // pure spin-polarized local density exchange-correlation functional

        if (fstr::upcase(xcLabel) == "SLDA")
        {
            return CXCFunctional({"SLDA"}, xcfun::lda, 0.0, {setPrimitiveSlaterFunctional(), setPrimitiveVWN3Functional()}, {1.0, 1.0});
        }

        // pure spin-polarized Becke/Slater exchange functional

        if (fstr::upcase(xcLabel) == "B88X")
        {
            return CXCFunctional({"B88X"}, xcfun::gga, 0.0, {setPrimitiveSlaterFunctional(), setPrimitiveBecke88Functional()}, {1.0, 1.0});
        }

        // pure spin-polarized BLYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BLYP")
        {
            return CXCFunctional({"BLYP"},
                                 xcfun::gga,
                                 0.0,
                                 {setPrimitiveSlaterFunctional(), setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional()},
                                 {1.0, 1.0, 1.0});
        }

        // pure spin-polarized hybrid B3LYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "B3LYP")
        {
            return CXCFunctional(
                {"B3LYP"},
                xcfun::gga,
                0.2,
                {setPrimitiveSlaterFunctional(), setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional(), setPrimitiveVWN3Functional()},
                {0.8, 0.72, 0.81, 0.19});
        }

        // pure spin-polarized hybrid BHANDH exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BHANDH")
        {
            return CXCFunctional({"BHANDH"}, xcfun::gga, 0.5, {setPrimitiveSlaterFunctional(), setPrimitiveLYPFunctional()}, {0.5, 1.0});
        }

        // pure spin-polarized hybrid B3LYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BHANDHLYP")
        {
            return CXCFunctional({"BHANDHLYP"},
                                 xcfun::gga,
                                 0.5,
                                 {setPrimitiveSlaterFunctional(), setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional()},
                                 {0.5, 0.5, 1.0});
        }

        // PKZB exchange functional

        if (fstr::upcase(xcLabel) == "PKZB") return vxcfuncs::setPkzbFunctional();

        // FIX ME: add other functionals here...
    }

    return CXCFunctional();
}

}  // namespace vxcfuncs
