//
//                           VELOXCHEM 1.0-RC
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
    if (fstr::upcase(xcLabel) == "HF")
    {
        _setFourParameters(1.0000, 1.61679827, 0.44959224, 3.35743605);
    }
    else if (fstr::upcase(xcLabel) == "BLYP")
    {
        _setFourParameters(1.0000, 2.34076671, 0.44488865, 4.09330090);
    }
    else if (fstr::upcase(xcLabel) == "B3LYP")
    {
        _setFourParameters(1.0000, 2.02929367, 0.40868035, 4.53807137);
    }
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
