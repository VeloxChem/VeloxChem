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

#include "BasisFunction.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "CustomViews.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"

CBasisFunction::CBasisFunction()

    : _exponents{}

    , _norms{}

    , _angular_momentum(-1)
{
}

CBasisFunction::CBasisFunction(const std::vector<double> &exponents, const std::vector<double> &norms, const int angular_momentum)

    : _exponents(exponents)

    , _norms(norms)

    , _angular_momentum(angular_momentum)
{
}

CBasisFunction::CBasisFunction(const CBasisFunction &other)

    : _exponents(other._exponents)

    , _norms(other._norms)

    , _angular_momentum(other._angular_momentum)
{
}

CBasisFunction::CBasisFunction(CBasisFunction &&other) noexcept

    : _exponents(std::move(other._exponents))

    , _norms(std::move(other._norms))

    , _angular_momentum(std::move(other._angular_momentum))
{
}

auto
CBasisFunction::operator=(const CBasisFunction &other) -> CBasisFunction &
{
    _exponents = other._exponents;

    _norms = other._norms;

    _angular_momentum = other._angular_momentum;

    return *this;
}

auto
CBasisFunction::operator=(CBasisFunction &&other) noexcept -> CBasisFunction &
{
    if (this != &other)
    {
        _exponents = std::move(other._exponents);

        _norms = std::move(other._norms);

        _angular_momentum = std::move(other._angular_momentum);
    }

    return *this;
}

auto
CBasisFunction::operator==(const CBasisFunction &other) const -> bool
{
    if (_angular_momentum != other._angular_momentum)
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
        return std::ranges::equal(_norms, other._norms, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); });
    }
}

auto
CBasisFunction::operator!=(const CBasisFunction &other) const -> bool
{
    return !(*this == other);
}

auto
CBasisFunction::set_exponents(const std::vector<double> &exponents) -> void
{
    _exponents = exponents;
}

auto
CBasisFunction::set_normalization_factors(const std::vector<double> &norms) -> void
{
    _norms = norms;
}

auto
CBasisFunction::set_angular_momentum(const int angular_momentum) -> void
{
    _angular_momentum = angular_momentum;
}

auto
CBasisFunction::add(const double exponent, const double norm) -> void
{
    _exponents.push_back(exponent);

    _norms.push_back(norm);
}

auto
CBasisFunction::normalize() -> void
{
    // TODO: Implemented for l > 6
    if (_angular_momentum > 6) return;

    if (_exponents.size() == 1) _norms[0] = 1.0;

    _rescale();

    double fact = 0.0;

    std::ranges::for_each(views::triangular(_exponents.size()), [&](const auto &index) { fact += _overlap(index); });

    fact = 1.0 / std::sqrt(fact);

    std::ranges::for_each(_norms, [=](double &norm) { norm *= fact; });
}

auto
CBasisFunction::get_exponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CBasisFunction::get_normalization_factors() const -> std::vector<double>
{
    return _norms;
}

auto
CBasisFunction::get_angular_momentum() const -> int
{
    return _angular_momentum;
}

auto
CBasisFunction::number_of_primitive_functions() const -> size_t
{
    return _exponents.size();
}

auto
CBasisFunction::_rescale() -> void
{
    // NOTE: Primitive Gaussians are normalized using standard solid harmonic
    // formula (Supporting Info Eq. A13):
    // N_eta = 1 / sqrt(integral |S^l_m(r)exp(-eta r^2)| dr)
    //       = (2 eta / pi)^3/4 sqrt[(4 eta)^l / (2l - 1)!!]
    // J. Chem. Theory Comput. 2020, https://doi.org/10.1021/acs.jctc.9b01296
    
    constexpr auto fpi = 2.0 / mathconst::pi_value();

    std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()), [&](const auto i) { _norms[i] *= std::pow(_exponents[i] * fpi, 0.75); });

    if (_angular_momentum == 1)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()), [&](const auto i) { _norms[i] *= 2.0 * std::sqrt(_exponents[i]); });
    }
    else if (_angular_momentum == 2)
    {
        const double fact = 4.0 / std::sqrt(3.0);

        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()), [&](const auto i) { _norms[i] *= fact * _exponents[i]; });
    }
    else if (_angular_momentum == 3)
    {
        const double fact = 8.0 / std::sqrt(15.0);

        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()),
                              [&](const auto i) { _norms[i] *= fact * _exponents[i] * std::sqrt(_exponents[i]); });
    }
    else if (_angular_momentum == 4)
    {
        const double fact = 16.0 / std::sqrt(105.0);

        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()),
                              [&](const auto i) { _norms[i] *= fact * _exponents[i] * _exponents[i]; });
    }
    else if (_angular_momentum == 5)
    {
        const double fact = 32.0 / std::sqrt(945.0);

        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()),
                              [&](const auto i) { _norms[i] *= fact * _exponents[i] * _exponents[i] * std::sqrt(_exponents[i]); });
    }
    else if (_angular_momentum == 6)
    {
        const double fact = 64.0 / std::sqrt(10395.0);

        std::ranges::for_each(std::views::iota(size_t{0}, _exponents.size()),
                              [&](const auto i) { _norms[i] *= fact * _exponents[i] * _exponents[i] * _exponents[i]; });
    }
    else
    {
        // TODO: implement l > 6
    }
}

auto
CBasisFunction::_overlap(const std::pair<size_t, size_t> &index) const -> double
{
    const auto [i, j] = index;

    const auto fab = 1.0 / (_exponents[i] + _exponents[j]);

    auto fovl = mathconst::pi_value() * fab;

    fovl = _norms[i] * _norms[j] * fovl * std::sqrt(fovl);

    if (i != j) fovl *= 2.0;

    if (_angular_momentum == 0)
    {
        return fovl;
    }
    else if (_angular_momentum == 1)
    {
        return 0.5 * fab * fovl;
    }
    else if (_angular_momentum == 2)
    {
        return 0.75 * fab * fab * fovl;
    }
    else if (_angular_momentum == 3)
    {
        return 1.875 * fab * fab * fab * fovl;
    }
    else if (_angular_momentum == 4)
    {
        return 6.5625 * fab * fab * fab * fab * fovl;
    }
    else if (_angular_momentum == 5)
    {
        return 29.53125 * fab * fab * fab * fab * fab * fovl;
    }
    else if (_angular_momentum == 6)
    {
        return 162.421875 * fab * fab * fab * fab * fab * fab * fovl;
    }
    else
    {
        // TODO: implement l > 6
        return 0.0;
    }
}
