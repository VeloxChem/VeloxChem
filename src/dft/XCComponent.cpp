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

#include "XCComponent.hpp"

#include <cstring>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CXCComponent::CXCComponent(const std::string& label, const double scalingFactor)

    : _label(fstr::upcase(label))

    , _scalingFactor(scalingFactor)
{
    _init_libxc_func();
}

CXCComponent::CXCComponent(const CXCComponent& source)

    : _label(source._label)

    , _scalingFactor(source._scalingFactor)
{
    _init_libxc_func();
}

CXCComponent::CXCComponent(CXCComponent&& source) noexcept

    : _label(std::move(source._label))

    , _scalingFactor(std::move(source._scalingFactor))
{
    _init_libxc_func();

    source._end_libxc_func();
}

CXCComponent::~CXCComponent()
{
    _end_libxc_func();
}

void
CXCComponent::_init_libxc_func()
{
    if (!_initialized)
    {
        auto funcID = xc_functional_get_number(_label.c_str());

        auto xc_err = xc_func_init(&_func, funcID, XC_POLARIZED);

        errors::assertMsgCritical(xc_err == 0, std::string("XCComponent: Invalid LibXC functional ") + _label);

        _initialized = true;
    }
}

void
CXCComponent::_end_libxc_func()
{
    if (_initialized)
    {
        xc_func_end(&_func);

        _initialized = false;
    }
}

void
CXCComponent::_reset_libxc_func()
{
    _end_libxc_func();

    _init_libxc_func();
}

CXCComponent&
CXCComponent::operator=(const CXCComponent& source)
{
    if (this == &source) return *this;

    _label = source._label;

    _scalingFactor = source._scalingFactor;

    _reset_libxc_func();

    return *this;
}

CXCComponent&
CXCComponent::operator=(CXCComponent&& source) noexcept
{
    if (this == &source) return *this;

    _label = std::move(source._label);

    _scalingFactor = std::move(source._scalingFactor);

    _reset_libxc_func();

    source._end_libxc_func();

    return *this;
}

bool
CXCComponent::operator==(const CXCComponent& other) const
{
    if (_label != other._label) return false;

    if (_scalingFactor != other._scalingFactor) return false;

    return true;
}

bool
CXCComponent::operator!=(const CXCComponent& other) const
{
    return !(*this == other);
}

std::string
CXCComponent::getLabel() const
{
    return _label;
}

double
CXCComponent::getScalingFactor() const
{
    return _scalingFactor;
}

const xc_func_type*
CXCComponent::getFunctionalPointer() const
{
    return &_func;
}

bool
CXCComponent::isLDA() const
{
    auto family = _func.info->family;

    return ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA));
}

bool
CXCComponent::isGGA() const
{
    auto family = _func.info->family;

    return ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA));
}

bool
CXCComponent::isMetaGGA() const
{
    auto family = _func.info->family;

    return ((family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA));
}
