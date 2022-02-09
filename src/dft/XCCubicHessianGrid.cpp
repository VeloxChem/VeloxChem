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

#include "XCCubicHessianGrid.hpp"

#include <iostream>

CXCCubicHessianGrid::CXCCubicHessianGrid()

    : _densityGridType(dengrid::undefined)

    , _xcGridType(xcfun::undefined)

    , _xcValues(CMemBlock2D<double>())
{
}

CXCCubicHessianGrid::CXCCubicHessianGrid(const CMemBlock2D<double>& xcValues, const dengrid densityGridType, const xcfun xcGridType)

    : _densityGridType(densityGridType)

    , _xcGridType(xcGridType)

    , _xcValues(xcValues)
{
}

CXCCubicHessianGrid::CXCCubicHessianGrid(const int32_t nGridPoints, const dengrid densityGridType, const xcfun xcGridType)
{
    _densityGridType = densityGridType;

    _xcGridType = xcGridType;

    int32_t ncomp = 0;

    if (_xcGridType == xcfun::lda) ncomp = (_densityGridType == dengrid::ab) ? 4 : 1;

    if (_xcGridType == xcfun::gga) ncomp = (_densityGridType == dengrid::ab) ? 31 : 3;

    // NOTE: this needs to be checked with mgga functionals implementation
    if (_xcGridType == xcfun::mgga) ncomp = (_densityGridType == dengrid::ab) ? 28 : 6;

    _xcValues = CMemBlock2D<double>(nGridPoints, ncomp);
}

CXCCubicHessianGrid::CXCCubicHessianGrid(const CXCCubicHessianGrid& source)

    : _densityGridType(source._densityGridType)

    , _xcGridType(source._xcGridType)

    , _xcValues(source._xcValues)
{
}

CXCCubicHessianGrid::CXCCubicHessianGrid(CXCCubicHessianGrid&& source) noexcept

    : _densityGridType(std::move(source._densityGridType))

    , _xcGridType(std::move(source._xcGridType))

    , _xcValues(std::move(source._xcValues))
{
}

CXCCubicHessianGrid::~CXCCubicHessianGrid()
{
}

CXCCubicHessianGrid&
CXCCubicHessianGrid::operator=(const CXCCubicHessianGrid& source)
{
    if (this == &source) return *this;

    _densityGridType = source._densityGridType;

    _xcGridType = source._xcGridType;

    _xcValues = source._xcValues;

    return *this;
}

CXCCubicHessianGrid&
CXCCubicHessianGrid::operator=(CXCCubicHessianGrid&& source) noexcept
{
    if (this == &source) return *this;

    _densityGridType = std::move(source._densityGridType);

    _xcGridType = std::move(source._xcGridType);

    _xcValues = std::move(source._xcValues);

    return *this;
}

bool
CXCCubicHessianGrid::operator==(const CXCCubicHessianGrid& other) const
{
    if (_densityGridType != other._densityGridType) return false;

    if (_xcGridType != other._xcGridType) return false;

    if (_xcValues != other._xcValues) return false;

    return true;
}

bool
CXCCubicHessianGrid::operator!=(const CXCCubicHessianGrid& other) const
{
    return !(*this == other);
}

void
CXCCubicHessianGrid::zero()
{
    _xcValues.zero();
}

int32_t
CXCCubicHessianGrid::getNumberOfGridPoints() const
{
    return _xcValues.size(0);
}

const double*
CXCCubicHessianGrid::xcCubicHessianValues(const xcvars iVariableType, const xcvars jVariableType, const xcvars kVariableType) const
{
    if (_xcGridType == xcfun::lda)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(0);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return _xcValues.data(0);
                }

                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(3);

                    if (_densityGridType == dengrid::lima) return _xcValues.data(0);

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
    }
    if (_xcGridType == xcfun::gga)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(0);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(3);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(4);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(5);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(6);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(7);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(8);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(9);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(10);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(11);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(12);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(13);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(14);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(15);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(16);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(17);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(18);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(19);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(20);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(21);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(22);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::grada)
        {
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(23);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(24);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(25);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(26);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(27);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(28);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::gradb)
        {
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(29);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::gradab)
        {
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(30);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
    }
    return nullptr;
}

double*
CXCCubicHessianGrid::xcCubicHessianValues(const xcvars iVariableType, const xcvars jVariableType, const xcvars kVariableType)
{
    if (_xcGridType == xcfun::lda)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(0);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return _xcValues.data(0);
                }

                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(3);

                    if (_densityGridType == dengrid::lima) return _xcValues.data(0);

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
    }
    if (_xcGridType == xcfun::gga)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (kVariableType == xcvars::rhoa)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(0);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(1);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(2);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(3);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(4);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(5);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(6);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(7);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(8);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(9);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(10);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(11);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(12);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(13);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(14);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (kVariableType == xcvars::rhob)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(15);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(16);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }

                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(17);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(18);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(19);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(20);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(21);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(22);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::grada)
        {
            if (jVariableType == xcvars::grada)
            {
                if (kVariableType == xcvars::grada)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(23);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(24);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(25);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradb)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(26);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(27);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(28);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::gradb)
        {
            if (jVariableType == xcvars::gradb)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(29);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
        if (iVariableType == xcvars::gradab)
        {
            if (jVariableType == xcvars::gradab)
            {
                if (kVariableType == xcvars::gradab)
                {
                    if (_densityGridType == dengrid::ab) return _xcValues.data(30);

                    if (_densityGridType == dengrid::lima) return nullptr;

                    if (_densityGridType == dengrid::limb) return nullptr;
                }
            }
        }
    }

    return nullptr;
}

std::ostream&
operator<<(std::ostream& output, const CXCCubicHessianGrid& source)
{
    output << std::endl;

    output << "[CXCCubicHessianGrid (Object):" << &source << "]" << std::endl;

    output << "_densityGridType: " << to_string(source._densityGridType) << std::endl;

    output << "_xcGridType: " << to_string(source._xcGridType) << std::endl;

    output << "_xcValues: " << std::endl;

    output << source._xcValues << std::endl;

    return output;
}
