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
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#include "DispersionParameters.hpp"

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CDispersionParameters::CDispersionParameters()
{
    _setDefaultParameters();
}

CDispersionParameters::CDispersionParameters(const std::string& xcLabel)
{
    _setDefaultParameters();

    _setFunctionalParameters(xcLabel);
}

CDispersionParameters::~CDispersionParameters()
{
}

void
CDispersionParameters::_setDefaultParameters()
{
    _s6 = 0.0;

    _s8 = 0.0;

    _s10 = 0.0;

    _a1 = 0.0;

    _a2 = 0.0;

    _s9 = 1.0;

    _alp = 16;

    _beta = 1.0;
}

void
CDispersionParameters::_setFunctionalParameters(const std::string& xcLabel)
{
    if (fstr::upcase(xcLabel) == "HF") _setFourParameters(1.0000, 1.61679827, 0.44959224, 3.35743605);

    else if (fstr::upcase(xcLabel) == "BLYP") _setFourParameters(1.0000, 2.34076671, 0.44488865, 4.09330090);

    else if (fstr::upcase(xcLabel) == "B3LYP") _setFourParameters(1.0000, 2.02929367, 0.40868035, 4.53807137);

    else if (fstr::upcase(xcLabel) == "PBE") _setFourParameters(1.0000, 0.95948085, 0.38574991, 4.80688534);

    else if (fstr::upcase(xcLabel) == "PBE0") _setFourParameters(1.0000, 1.20065498, 0.40085597, 5.02928789);

    else if (fstr::upcase(xcLabel) == "REVPBE") _setFourParameters(1.0000, 1.74676530, 0.53634900, 3.07261485);

    else if (fstr::upcase(xcLabel) == "PW91") _setFourParameters(1.0000, 0.77283111, 0.39581542, 4.93405761);

    else if (fstr::upcase(xcLabel) == "OLYP") _setFourParameters(1.0000, 2.74836820, 0.60184498, 2.53292167);

    else if (fstr::upcase(xcLabel) == "O3LYP") _setFourParameters(1.0000, 1.75762508, 0.10348980, 6.16233282);

    else if (fstr::upcase(xcLabel) == "X3LYP") _setFourParameters(1.0000, 1.54701429, 0.20318443, 5.61852648);

    else if (fstr::upcase(xcLabel) == "B97") _setFourParameters(1.0000, 0.87854260, 0.29319126, 4.51647719);

    else if (fstr::upcase(xcLabel) == "TPSS") _setFourParameters(1.0000, 1.76596355, 0.42822303, 4.54257102);

    else if (fstr::upcase(xcLabel) == "TPSSH") _setFourParameters(1.0000, 1.85897750, 0.44286966, 4.60230534);

    else if (fstr::upcase(xcLabel) == "REVTPSS") _setFourParameters(1.0000, 1.53089454, 0.44880597, 4.64042317);

    else if (fstr::upcase(xcLabel) == "SCAN") _setFourParameters(1.0000, 1.46126056, 0.62930855, 6.31284039);

    else if (fstr::upcase(xcLabel) == "M06") _setFourParameters(1.0000, 0.16366729, 0.53456413, 6.06192174);

    else if (fstr::upcase(xcLabel) == "M06-L") _setFourParameters(1.0000, 0.59493760, 0.71422359, 6.35314182);

    else if (fstr::upcase(xcLabel) == "MPW1B95") _setFourParameters(1.0000, 0.50093024, 0.41585097, 4.99154869);

    else if (fstr::upcase(xcLabel) == "MPWB1K") _setFourParameters(1.0000, 0.57338313, 0.44687975, 5.21266777);

    else if (fstr::upcase(xcLabel) == "PW6B95") _setFourParameters(1.0000, -0.31926054, 0.04142919, 5.84655608);

    else
    {
        std::string errmsg("DispersionParameters: D4 dispersion model not implemented for functional ");

        errors::assertMsgCritical(false, errmsg + xcLabel);
    }
}

void
CDispersionParameters::_setFourParameters(const double s6, const double s8, const double a1, const double a2)
{
    _s6 = s6;

    _s8 = s8;

    _a1 = a1;

    _a2 = a2;
}

double
CDispersionParameters::getS6() const
{
    return _s6;
}

double
CDispersionParameters::getS8() const
{
    return _s8;
}

double
CDispersionParameters::getS10() const
{
    return _s10;
}

double
CDispersionParameters::getA1() const
{
    return _a1;
}

double
CDispersionParameters::getA2() const
{
    return _a2;
}

int32_t
CDispersionParameters::getAlp() const
{
    return _alp;
}
