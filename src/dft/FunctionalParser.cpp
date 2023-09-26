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

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

namespace vxcfuncs {  // vxcfuncs namespace

std::vector<std::string>
getAvailableFunctionals()
{
    return std::vector<std::string>({"SLATER", "SLDA",   "B88X",    "BLYP",  "B3LYP",   "BHANDH", "BHANDHLYP", "PBE",   "PBE0",   "REVPBE",
                                     "BP86",   "PW91",   "MPW1K",   "OLYP",  "O3LYP",   "X3LYP",  "B97",       "B97-1", "B97-2",  "B97-3",
                                     "TPSS",   "TPSSH",  "REVTPSS", "PKZB",  "SCAN",    "RSCAN",  "R2SCAN",    "M05",   "M05-2X", "M06",
                                     "M06-2X", "M06-HF", "M06-L",   "M11-L", "MPW1B95", "MPWB1K", "PW6B95",    "PWB6K"});
}

CXCFunctional
getExchangeCorrelationFunctional(const std::string &xcLabel)
{
    auto availFuncs = getAvailableFunctionals();

    if (std::find(availFuncs.begin(), availFuncs.end(), fstr::upcase(xcLabel)) != availFuncs.end())
    {
        // LDA

        if (fstr::upcase(xcLabel) == "SLATER") return CXCFunctional("SLATER", {"LDA_X"}, {1.0});

        if (fstr::upcase(xcLabel) == "SLDA") return CXCFunctional("SLDA", {"LDA_X", "LDA_C_VWN_RPA"}, {1.0, 1.0});

        // GGA

        if (fstr::upcase(xcLabel) == "B88X") return CXCFunctional("B88X", {"GGA_X_B88"}, {1.0});

        if (fstr::upcase(xcLabel) == "BLYP") return CXCFunctional("BLYP", {"GGA_X_B88", "GGA_C_LYP"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "B3LYP")
        {
            return CXCFunctional("B3LYP", {"LDA_X", "GGA_X_B88", "LDA_C_VWN_RPA", "GGA_C_LYP"}, {0.08, 0.72, 0.19, 0.81}, 0.2);
        }

        if (fstr::upcase(xcLabel) == "BHANDH") return CXCFunctional("BHANDH", {"LDA_X", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        if (fstr::upcase(xcLabel) == "BHANDHLYP") return CXCFunctional("BHANDHLYP", {"GGA_X_B88", "GGA_C_LYP"}, {0.5, 1.0}, 0.5);

        if (fstr::upcase(xcLabel) == "PBE") return CXCFunctional("PBE", {"GGA_X_PBE", "GGA_C_PBE"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "REVPBE") return CXCFunctional("REVPBE", {"GGA_X_PBE_R", "GGA_C_PBE"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PBE0") return CXCFunctional("PBE0", {"GGA_X_PBE", "GGA_C_PBE"}, {0.75, 1.0}, 0.25);

        if (fstr::upcase(xcLabel) == "BP86") return CXCFunctional("BP86", {"GGA_X_B88", "GGA_C_P86"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PW91") return CXCFunctional("PW91", {"GGA_X_PW91_MOD", "GGA_C_PW91"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "MPW1K") return CXCFunctional("MPW1K", {"HYB_GGA_XC_MPW1K"}, {1.0});

        if (fstr::upcase(xcLabel) == "OLYP") return CXCFunctional("OLYP", {"GGA_X_OPTX", "GGA_C_LYP"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "O3LYP") return CXCFunctional("O3LYP", {"HYB_GGA_XC_O3LYP"}, {1.0});

        if (fstr::upcase(xcLabel) == "X3LYP")
        {
            return CXCFunctional(
                "X3LYP", {"LDA_X", "GGA_X_B88", "LDA_C_VWN_RPA", "GGA_C_LYP", "GGA_X_PW91_MOD"}, {0.073, 0.542, 0.129, 0.871, 0.167}, 0.218);
        }

        if (fstr::upcase(xcLabel) == "B97") return CXCFunctional("B97", {"HYB_GGA_XC_B97"}, {1.0});

        if (fstr::upcase(xcLabel) == "B97-1") return CXCFunctional("B97-1", {"HYB_GGA_XC_B97_1"}, {1.0});

        if (fstr::upcase(xcLabel) == "B97-2") return CXCFunctional("B97-2", {"HYB_GGA_XC_B97_2"}, {1.0});

        if (fstr::upcase(xcLabel) == "B97-3") return CXCFunctional("B97-3", {"HYB_GGA_XC_B97_3"}, {1.0});

        // meta-GGA

        if (fstr::upcase(xcLabel) == "TPSS") return CXCFunctional("TPSS", {"MGGA_X_TPSS", "MGGA_C_TPSS"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "TPSSH") return CXCFunctional("TPSSH", {"MGGA_X_TPSS", "MGGA_C_TPSS"}, {0.9, 1.0}, 0.1);

        if (fstr::upcase(xcLabel) == "REVTPSS") return CXCFunctional("REVTPSS", {"MGGA_X_REVTPSS", "MGGA_C_REVTPSS"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "PKZB") return CXCFunctional("PKZB", {"MGGA_X_PKZB", "MGGA_C_PKZB"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "SCAN") return CXCFunctional("SCAN", {"MGGA_X_SCAN", "MGGA_C_SCAN"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "RSCAN") return CXCFunctional("RSCAN", {"MGGA_X_RSCAN", "MGGA_C_RSCAN"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "R2SCAN") return CXCFunctional("R2SCAN", {"MGGA_X_R2SCAN", "MGGA_C_R2SCAN"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M05") return CXCFunctional("M05", {"HYB_MGGA_X_M05", "MGGA_C_M05"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M05-2X") return CXCFunctional("M05-2X", {"HYB_MGGA_X_M05_2X", "MGGA_C_M05_2X"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06") return CXCFunctional("M06", {"HYB_MGGA_X_M06", "MGGA_C_M06"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06-2X") return CXCFunctional("M06-2X", {"HYB_MGGA_X_M06_2X", "MGGA_C_M06_2X"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06-HF") return CXCFunctional("M06-HF", {"HYB_MGGA_X_M06_HF", "MGGA_C_M06_HF"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M06-L") return CXCFunctional("M06-L", {"MGGA_X_M06_L", "MGGA_C_M06_L"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "M11-L") return CXCFunctional("M11-L", {"MGGA_X_M11_L", "MGGA_C_M11_L"}, {1.0, 1.0});

        if (fstr::upcase(xcLabel) == "MPW1B95") return CXCFunctional("MPW1B95", {"HYB_MGGA_XC_MPW1B95"}, {1.0});

        if (fstr::upcase(xcLabel) == "MPWB1K") return CXCFunctional("MPWB1K", {"HYB_MGGA_XC_MPWB1K"}, {1.0});

        if (fstr::upcase(xcLabel) == "PW6B95") return CXCFunctional("PW6B95", {"HYB_MGGA_XC_PW6B95"}, {1.0});

        if (fstr::upcase(xcLabel) == "PWB6K") return CXCFunctional("PWB6K", {"HYB_MGGA_XC_PWB6K"}, {1.0});

        // TODO add more functionals here...
    }

    std::string errmsg(std::string("getExchangeCorrelationFunctional: ") + xcLabel + std::string(" is not available"));

    errors::assertMsgCritical(false, errmsg);

    return CXCFunctional("Undefined", {}, {});
}

}  // namespace vxcfuncs
