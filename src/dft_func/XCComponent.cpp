//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "XCComponent.hpp"

#include <cstring>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CXCComponent::CXCComponent(const std::string& label, const double scalingFactor)

    : _label(format::upper_case(label))

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

xc_func_type*
CXCComponent::getFunctionalPointer()
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
