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

#include <cmath>

#include "MathConst.hpp"

CBasisFunction::CBasisFunction(const std::vector<double>& exponents, const std::vector<double>& norms, const int64_t angmom)

    : _exponents(exponents)

    , _norms(norms)

    , _angmom(angmom)

{
}

auto
CBasisFunction::setExponents(const std::vector<double>& exponents) -> void
{
    _exponents = exponents;
}

auto
CBasisFunction::setNormalizationFactors(const std::vector<double>& norms) -> void
{
    _norms = norms;
}

auto
CBasisFunction::setAngularMomentum(const int64_t angmom) -> void
{
    _angmom = angmom;
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
    // NOTE: Currently implemented for l = 0..6

    if (_angmom > 6) return;

    if (const auto npgtos = _exponents.size(); npgtos > 0)
    {
        // uncontracted GTO

        if (npgtos == 1) _norms[0] = 1.0;

        // normalize primitive GTOs

        _rescale();

        // compute contracted GTO self-overlap

        double fovl = 0.0;

        for (size_t i = 0; i < npgtos; i++)
        {
            fovl += _overlap(i, i);

            for (size_t j = i + 1; j < npgtos; j++)
            {
                fovl += 2.0 * _overlap(i, j);
            }
        }

        // renormaliza primitive GTOs

        fovl = 1.0 / std::sqrt(fovl);

        for (size_t i = 0; i < npgtos; i++)
        {
            _norms[i] *= fovl;
        }
    }
}

auto
CBasisFunction::getExponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CBasisFunction::getNormalizationFactors() const -> std::vector<double>
{
    return _norms;
}

auto
CBasisFunction::getAngularMomentum() const -> int64_t
{
    return _angmom;
}

auto
CBasisFunction::getNumberOfPrimitiveFunctions() const -> int64_t
{
    return static_cast<int64_t>(_exponents.size());
}

auto
CBasisFunction::_rescale() -> void
{
    if (const auto npgtos = _exponents.size(); npgtos > 0)
    {
        const auto fpi = 2.0 / mathconst::getPiValue();

        for (size_t i = 0; i < npgtos; i++)
        {
            _norms[i] *= std::pow(_exponents[i] * fpi, 0.75);
        }

        if (_angmom == 1)
        {
            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= 2.0 * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 2)
        {
            const double fact = 2.0 / std::sqrt(3.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i];
            }
        }
        else if (_angmom == 3)
        {
            const double fact = 4.0 / std::sqrt(15.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 4)
        {
            const double fact = 2.0 / std::sqrt(105.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i];
            }
        }
        else if (_angmom == 5)
        {
            const double fact = 4.0 / std::sqrt(945.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i]

                             * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 6)
        {
            const double fact = 4.0 / std::sqrt(10395.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i]

                             * _exponents[i];
            }
        }
        else
        {
            // implement l > 6
        }
    }
}

auto
CBasisFunction::_overlap(const size_t i, const size_t j) const -> double
{
    const auto fab = 1.0 / (_exponents[i] + _exponents[j]);

    const auto fovl = _norms[i] * _norms[j]

                      * std::pow(mathconst::getPiValue() * fab, 1.5);

    if (_angmom == 0)
    {
        return fovl;
    }
    else if (_angmom == 1)
    {
        return 0.5 * fab * fovl;
    }
    else if (_angmom == 2)
    {
        return 3.0 * fab * fab * fovl;
    }
    else if (_angmom == 3)
    {
        return 7.5 * fab * fab * fab * fovl;
    }
    else if (_angmom == 4)
    {
        return 420.0 * fab * fab * fab * fab * fovl;
    }
    else if (_angmom == 5)
    {
        return 1890.0 * fab * fab * fab * fab * fab * fovl;
    }
    else if (_angmom == 6)
    {
        return 41580.0 * fab * fab * fab * fab * fab * fab * fovl;
    }
    else
    {
        return 0.0;
    }
}
