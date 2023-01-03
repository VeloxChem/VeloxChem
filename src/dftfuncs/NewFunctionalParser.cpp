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
    return std::vector<std::string>({
            "SLATER", "VWN_RPA", "SLDA",
            "B88X", "LYP", "BLYP", "B3LYP", "BHANDH", "BHANDHLYP", "PBE", "PBE0", "BP86", "PW91", "OLYP",
            "TPSS", "TPSSH", "PKZB", "SCAN", "M06", "M06L", "M11L"});
}

CXCNewFunctional
getExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailableFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // LDA

        if (fstr::upcase(xcLabel) == "SLATER") return CXCNewFunctional("SLATER", {"LDA_X"}, {1.0});

        if (fstr::upcase(xcLabel) == "VWN_RPA") return CXCNewFunctional("VWN_RPA", {"LDA_C_VWN_RPA"}, {1.0});

        if (fstr::upcase(xcLabel) == "SLDA") return CXCNewFunctional("SLDA", {"LDA_X", "LDA_C_VWN_RPA"}, {1.0, 1.0});

        // GGA

        if (fstr::upcase(xcLabel) == "B88X") return CXCNewFunctional("B88X", {"GGA_X_B88"}, {1.0});

        if (fstr::upcase(xcLabel) == "LYP") return CXCNewFunctional("LYP", {"GGA_C_LYP"}, {1.0});

        if (fstr::upcase(xcLabel) == "BLYP") return CXCNewFunctional("BLYP", {"GGA_X_B88", "GGA_C_LYP"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "B3LYP")
        {
            return CXCNewFunctional("B3LYP", {"LDA_X", "GGA_X_B88", "LDA_C_VWN_RPA", "GGA_C_LYP"}, {0.08, 0.72, 0.19, 0.81}, 0.2);
        }

        if (fstr::upcase(xcLabel) == "BHANDH") return CXCNewFunctional("BHANDH", {"LDA_X", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        if (fstr::upcase(xcLabel) == "BHANDHLYP") return CXCNewFunctional("BHANDHLYP", {"GGA_X_B88", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        if (fstr::upcase(xcLabel) == "PBE") return CXCNewFunctional("PBE", {"GGA_X_PBE", "GGA_C_PBE"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PBE0") return CXCNewFunctional("PBE0", {"GGA_X_PBE", "GGA_C_PBE"}, {0.75, 1.0}, 0.25);

        if (fstr::upcase(xcLabel) == "BP86") return CXCNewFunctional("BP86", {"GGA_X_B88", "GGA_C_P86"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PW91") return CXCNewFunctional("PW91", {"GGA_X_PW91", "GGA_C_PW91"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "OLYP") return CXCNewFunctional("OLYP", {"GGA_X_OPTX", "GGA_C_LYP"}, {1.0, 1.0});

        // meta-GGA

        if (fstr::upcase(xcLabel) == "TPSS") return CXCNewFunctional("TPSS", {"MGGA_X_TPSS", "MGGA_C_TPSS"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "TPSSH") return CXCNewFunctional("TPSSH", {"MGGA_X_TPSS", "MGGA_C_TPSS"}, {0.9, 1.0}, 0.1);

        if (fstr::upcase(xcLabel) == "SCAN") return CXCNewFunctional("SCAN", {"MGGA_X_SCAN", "MGGA_C_SCAN"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06") return CXCNewFunctional("M06", {"HYB_MGGA_X_M06", "MGGA_C_M06"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06L") return CXCNewFunctional("M06L", {"MGGA_X_M06_L", "MGGA_C_M06_L"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M11L") return CXCNewFunctional("M11L", {"MGGA_X_M11_L", "MGGA_C_M11_L"}, {1.0, 1.0});

        // TODO add more functionals here...
    }

    std::string errmsg(std::string("getExchangeCorrelationFunctional: Cannot find functional ") + xcLabel);

    errors::assertMsgCritical(false, errmsg);

    return CXCNewFunctional("Undefined", {}, {});
}

std::vector<std::string>
getAvailablePairDensityFunctionals()
{
    return std::vector<std::string>({"PLDA"});
}

CXCPairDensityFunctional
getPairDensityExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailablePairDensityFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // Pair-density local density exchange-correlation functional

        if (fstr::upcase(xcLabel) == "PLDA") return CXCPairDensityFunctional("PLDA", {"PSLATER", "PVWN"}, {1.0, 1.0});

        // FIX ME: add other functionals here...
    }

    std::string errmsg(std::string("getPairDensityExchangeCorrelationFunctional: Cannot find functional ") + xcLabel);

    errors::assertMsgCritical(false, errmsg);

    return CXCPairDensityFunctional("Undefined", {}, {});
}

}  // namespace newvxcfuncs
