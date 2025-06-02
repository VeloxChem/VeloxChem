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

#include "BaseCorePotential.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "MathFunc.hpp"

CBaseCorePotential::CBaseCorePotential()

    : _exponents{}

    , _factors{}

    , _radial_orders{}
{
}

CBaseCorePotential::CBaseCorePotential(const std::vector<double> &exponents,
                                       const std::vector<double> &factors,
                                       const std::vector<int>    &radial_orders)

    : _exponents(exponents)

    , _factors(factors)

    , _radial_orders(radial_orders)
{
}

CBaseCorePotential::CBaseCorePotential(const CBaseCorePotential &other)

    : _exponents(other._exponents)

    , _factors(other._factors)

    , _radial_orders(other._radial_orders)
{
}

CBaseCorePotential::CBaseCorePotential(CBaseCorePotential &&other) noexcept

    : _exponents{}

    , _factors{}

    , _radial_orders{}
{
    std::swap(_exponents, other._exponents);

    std::swap(_factors, other._factors);

    std::swap(_radial_orders, other._radial_orders);
}

auto
CBaseCorePotential::operator=(const CBaseCorePotential &other) -> CBaseCorePotential &
{
    _exponents = other._exponents;

    _factors = other._factors;

    _radial_orders = other._radial_orders;

    return *this;
}

auto
CBaseCorePotential::operator=(CBaseCorePotential &&other) noexcept -> CBaseCorePotential &
{
    std::swap(_exponents, other._exponents);

    std::swap(_factors, other._factors);

    std::swap(_radial_orders, other._radial_orders);

    return *this;
}

auto
CBaseCorePotential::operator==(const CBaseCorePotential &other) const -> bool
{
    if (_radial_orders != other._radial_orders)
    {
        return false;
    }
    else if (!std::ranges::equal(
                 _exponents, other._exponents, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else
    {
        return std::ranges::equal(_factors, other._factors, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); });
    }
}

auto
CBaseCorePotential::operator!=(const CBaseCorePotential &other) const -> bool
{
    return !(*this == other);
}

auto
CBaseCorePotential::set_exponents(const std::vector<double> &exponents) -> void
{
    _exponents = exponents;
}

auto
CBaseCorePotential::set_factors(const std::vector<double> &factors) -> void
{
    _factors = factors;
}

auto
CBaseCorePotential::set_radial_orders(const std::vector<int> &radial_orders) -> void
{
    _radial_orders = radial_orders; 
}

auto
CBaseCorePotential::add(const double exponent, const double factor, const int radial_order) -> void
{
    _exponents.push_back(exponent);

    _factors.push_back(factor);
    
    _radial_orders.push_back(radial_order);
}

auto
CBaseCorePotential::get_exponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CBaseCorePotential::get_factors() const -> std::vector<double>
{
    return _factors;
}

auto
CBaseCorePotential::get_radial_orders() const -> std::vector<int>
{
    return _radial_orders;
}

auto
CBaseCorePotential::number_of_primitive_potentials() const -> size_t
{
    return _exponents.size();
}

auto
CBaseCorePotential::is_valid_radial_orders() const -> bool
{
    return (std::ranges::find_if(_radial_orders, [](int order) {return order != 2;}) == _radial_orders.end());
}
