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

#include "RecursionFunction.hpp"

CRecursionFunction::CRecursionFunction()

    : _label(std::string())

    , _funcAction(nullptr)
{
}

CRecursionFunction::CRecursionFunction(const std::string& label, const std::function<def_rec_func_typ>& funcAction)

    : _label(label)

    , _funcAction(funcAction)
{
}

CRecursionFunction::CRecursionFunction(const CRecursionFunction& source)

    : _label(source._label)

    , _funcAction(source._funcAction)
{
}

CRecursionFunction::CRecursionFunction(CRecursionFunction&& source) noexcept

    : _label(std::move(source._label))

    , _funcAction(std::move(source._funcAction))
{
}

CRecursionFunction::~CRecursionFunction()
{
}

CRecursionFunction&
CRecursionFunction::operator=(const CRecursionFunction& source)
{
    if (this == &source) return *this;

    _label = source._label;

    _funcAction = source._funcAction;

    return *this;
}

CRecursionFunction&
CRecursionFunction::operator=(CRecursionFunction&& source) noexcept
{
    if (this == &source) return *this;

    _label = std::move(source._label);

    _funcAction = std::move(source._funcAction);

    return *this;
}

bool
CRecursionFunction::operator==(const CRecursionFunction& other) const
{
    if (this == &other) return true;

    if (_label != other._label) return false;

    return true;
}

bool
CRecursionFunction::operator!=(const CRecursionFunction& other) const
{
    return !((*this) == other);
}

std::vector<CRecursionTerm>
CRecursionFunction::compute(const CRecursionTerm& recursionTerm) const
{
    return _funcAction(recursionTerm);
}

std::string
CRecursionFunction::getLabel() const
{
    return _label;
}

bool
CRecursionFunction::isMatch(const std::string label) const
{
    if (label != _label) return false;

    return true;
}

std::ostream&
operator<<(std::ostream& output, const CRecursionFunction& source)
{
    output << std::endl;

    output << "[CRecFunction (Object):" << &source << "]" << std::endl;

    output << "_label: " << source._label << std::endl;

    output << "_funcAction: " << &(source._funcAction) << std::endl;

    return output;
}
