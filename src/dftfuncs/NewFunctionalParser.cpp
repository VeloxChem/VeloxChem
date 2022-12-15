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

#include "NewFunctionalParser.hpp"

#include <algorithm>

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

namespace newvxcfuncs {  // newvxcfuncs namespace

std::vector<std::string>
getAvailableFunctionals()
{
    return std::vector<std::string>({"SLATER", "VWN_RPA", "BECKE88", "LYP", "SLDA", "B88X", "BLYP", "PBE", "B3LYP", "BHANDH", "BHANDHLYP"});
}

CXCNewFunctional
getExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailableFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // Slater exchange functional

        if (fstr::upcase(xcLabel) == "SLATER") return CXCNewFunctional("LDA_X");

        // Vosko-Wilk-Nusair RPA correlation functional

        if (fstr::upcase(xcLabel) == "VWN_RPA") return CXCNewFunctional("LDA_C_VWN_RPA");

        // Becke (1988) exchange functional

        if (fstr::upcase(xcLabel) == "BECKE88") return CXCNewFunctional({"GGA_X_B88", "LDA_X"}, {1.0, -1.0});

        // Lee, Yang and Parr correlation functional

        if (fstr::upcase(xcLabel) == "LYP") return CXCNewFunctional("GGA_C_LYP");

        // local density exchange-correlation functional

        if (fstr::upcase(xcLabel) == "SLDA") return CXCNewFunctional({"LDA_X", "LDA_C_VWN_RPA"}, {1.0, 1.0});

        // Becke/Slater exchange functional

        if (fstr::upcase(xcLabel) == "B88X") return CXCNewFunctional("GGA_X_B88");

        // BLYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BLYP") return CXCNewFunctional({"GGA_X_B88", "GGA_C_LYP"}, {1.0, 1.0});

        // PBE exchange-correlation functional

        if (fstr::upcase(xcLabel) == "PBE") return CXCNewFunctional({"GGA_X_PBE", "GGA_C_PBE"}, {1.0, 1.0});

        // hybrid B3LYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "B3LYP")
        {
            return CXCNewFunctional({"LDA_X", "GGA_X_B88", "LDA_C_VWN_RPA", "GGA_C_LYP"}, {0.08, 0.72, 0.19, 0.81}, 0.2);
        }

        // hybrid BHANDH exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BHANDH") return CXCNewFunctional({"LDA_X", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        // hybrid BHANDHLYP exchange-correlation functional

        if (fstr::upcase(xcLabel) == "BHANDHLYP") return CXCNewFunctional({"GGA_X_B88", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        // PKZB exchange functional
        // if (fstr::upcase(xcLabel) == "PKZB") return ...

        // FIX ME: add other functionals here...
    }

    std::string errmsg(std::string("getExchangeCorrelationFunctional: Cannot find functional ") + xcLabel);

    errors::assertMsgCritical(false, errmsg);

    return CXCNewFunctional({}, {});
}

std::vector<std::string>
getAvailablePairDensityFunctionals()
{
    return std::vector<std::string>({"PLDA", "PPBE"});
}

CXCPairDensityFunctional
getPairDensityExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailablePairDensityFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // Pair-density local density exchange-correlation functional

        if (fstr::upcase(xcLabel) == "PLDA") return CXCPairDensityFunctional("PLDA", {"PSLATER", "PVWN"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PPBE") return CXCPairDensityFunctional("PPBE", {"PPBE_X", "PPBE_C"}, {1.0, 1.0});

        // FIX ME: add other functionals here...
    }

    std::string errmsg(std::string("getPairDensityExchangeCorrelationFunctional: Cannot find functional ") + xcLabel);

    errors::assertMsgCritical(false, errmsg);

    return CXCPairDensityFunctional("Undefined", {}, {});
}

}  // namespace newvxcfuncs
