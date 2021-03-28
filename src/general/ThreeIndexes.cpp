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

#include "ThreeIndexes.hpp"

#include <utility>

CThreeIndexes::CThreeIndexes()

    : _iIndex(-1)

    , _jIndex(-1)

    , _kIndex(-1)
{
}

CThreeIndexes::CThreeIndexes(const int32_t iIndex, const int32_t jIndex, const int32_t kIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)
{
}

CThreeIndexes::CThreeIndexes(const CThreeIndexes& source)

    : _iIndex(source._iIndex)

    , _jIndex(source._jIndex)

    , _kIndex(source._kIndex)
{
}

CThreeIndexes::CThreeIndexes(CThreeIndexes&& source) noexcept

    : _iIndex(std::move(source._iIndex))

    , _jIndex(std::move(source._jIndex))

    , _kIndex(std::move(source._kIndex))
{
}

CThreeIndexes::~CThreeIndexes()
{
}

CThreeIndexes&
CThreeIndexes::operator=(const CThreeIndexes& source)
{
    if (this == &source) return *this;

    _iIndex = source._iIndex;

    _jIndex = source._jIndex;

    _kIndex = source._kIndex;

    return *this;
}

CThreeIndexes&
CThreeIndexes::operator=(CThreeIndexes&& source) noexcept
{
    if (this == &source) return *this;

    _iIndex = std::move(source._iIndex);

    _jIndex = std::move(source._jIndex);

    _kIndex = std::move(source._kIndex);

    return *this;
}

bool
CThreeIndexes::operator==(const CThreeIndexes& other) const
{
    if (this == &other) return true;

    if (_iIndex != other._iIndex) return false;

    if (_jIndex != other._jIndex) return false;

    if (_kIndex != other._kIndex) return false;

    return true;
}

bool
CThreeIndexes::operator!=(const CThreeIndexes& other) const
{
    return !((*this) == other);
}

int32_t
CThreeIndexes::first() const
{
    return _iIndex;
}

int32_t
CThreeIndexes::second() const
{
    return _jIndex;
}

int32_t
CThreeIndexes::third() const
{
    return _kIndex;
}

bool
CThreeIndexes::isValidTriple() const
{
    if (_iIndex < 0) return false;

    if (_jIndex < 0) return false;

    if (_kIndex < 0) return false;

    return true;
}

std::ostream&
operator<<(std::ostream& output, const CThreeIndexes& source)
{
    output << std::endl;

    output << "[CThreeIndexes (Object):" << &source << "]" << std::endl;

    output << "_iIndex: " << source._iIndex << std::endl;

    output << "_jIndex: " << source._jIndex << std::endl;

    output << "_kIndex: " << source._kIndex << std::endl;

    return output;
}
