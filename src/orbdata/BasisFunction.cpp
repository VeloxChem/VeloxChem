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

#include "BasisFunction.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

#include <mpi.h>

#include "MathConst.hpp"
#include "MpiFunc.hpp"

CBasisFunction::CBasisFunction()

    : _angularMomentum(-1)
{
}

CBasisFunction::CBasisFunction(const std::vector<double>& exponents, const std::vector<double>& normFactors, const int32_t angularMomentum)

    : _exponents(exponents)

    , _normFactors(normFactors)

    , _angularMomentum(angularMomentum)

{
}

CBasisFunction::CBasisFunction(const CBasisFunction& source)

    : _exponents(source._exponents)

    , _normFactors(source._normFactors)

    , _angularMomentum(source._angularMomentum)
{
}

CBasisFunction::CBasisFunction(CBasisFunction&& source) noexcept

    : _exponents(std::move(source._exponents))

    , _normFactors(std::move(source._normFactors))

    , _angularMomentum(std::move(source._angularMomentum))
{
}

CBasisFunction::~CBasisFunction()
{
}

CBasisFunction&
CBasisFunction::operator=(const CBasisFunction& source)
{
    if (this == &source) return *this;

    _exponents = source._exponents;

    _normFactors = source._normFactors;

    _angularMomentum = source._angularMomentum;

    return *this;
}

CBasisFunction&
CBasisFunction::operator=(CBasisFunction&& source) noexcept
{
    if (this == &source) return *this;

    _exponents = std::move(source._exponents);

    _normFactors = std::move(source._normFactors);

    _angularMomentum = std::move(source._angularMomentum);

    return *this;
}

bool
CBasisFunction::operator==(const CBasisFunction& other) const
{
    if (_exponents.size() != other._exponents.size()) return false;

    for (size_t i = 0; i < _exponents.size(); i++)
    {
        if (std::fabs(_exponents[i] - other._exponents[i]) > 1.0e-13)
        {
            return false;
        }
    }

    if (_normFactors.size() != other._normFactors.size()) return false;

    for (size_t i = 0; i < _normFactors.size(); i++)
    {
        if (std::fabs(_normFactors[i] - other._normFactors[i]) > 1.0e-13)
        {
            return false;
        }
    }

    if (_angularMomentum != other._angularMomentum) return false;

    return true;
}

bool
CBasisFunction::operator!=(const CBasisFunction& other) const
{
    return !(*this == other);
}

void
CBasisFunction::setExponents(const std::vector<double>& exponents)
{
    _exponents = exponents;
}

void
CBasisFunction::setNormalizationFactors(const std::vector<double>& normFactors)
{
    _normFactors = normFactors;
}

void
CBasisFunction::setAngularMomentum(const int32_t angularMomentum)
{
    _angularMomentum = angularMomentum;
}

void
CBasisFunction::add(const double exponent, const double normFactor)
{
    _exponents.push_back(exponent);

    _normFactors.push_back(normFactor);
}

void
CBasisFunction::normalize()
{
    // NOTE: Currently implemented for l = 0..6

    if (_angularMomentum > 6) return;

    // uncontracted basis, set expansion coeficient to 1.0

    if (_normFactors.size() == 1) _normFactors[0] = 1.0;

    // normalize primitive GBFs

    _rescale();

    // compute overlap

    auto sz = _exponents.size();

    double ovl = 0.0;

    for (size_t i = 0; i < sz; i++)
    {
        ovl += _overlap(i, i);

        for (size_t j = i + 1; j < sz; j++)
        {
            ovl += 2.0 * _overlap(i, j);
        }
    }

    // renormaliza primitive BFs

    ovl = 1.0 / std::sqrt(ovl);

    for (size_t i = 0; i < sz; i++)
    {
        _normFactors[i] *= ovl;
    }
}

std::vector<double>
CBasisFunction::getExponents() const
{
    return _exponents;
}

std::vector<double>
CBasisFunction::getNormalizationFactors() const
{
    return _normFactors;
}

int32_t
CBasisFunction::getAngularMomentum() const
{
    return _angularMomentum;
}

int32_t
CBasisFunction::getNumberOfPrimitiveFunctions() const
{
    return static_cast<int32_t>(_exponents.size());
}

void
CBasisFunction::_rescale()
{
    auto fpi = 2.0 / mathconst::getPiValue();

    for (size_t i = 0; i < _exponents.size(); i++)
    {
        _normFactors[i] *= std::pow(_exponents[i] * fpi, 0.75);
    }

    if (_angularMomentum == 1)
    {
        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= 2.0 * std::sqrt(_exponents[i]);
        }

        return;
    }

    if (_angularMomentum == 2)
    {
        double f = 2.0 / std::sqrt(3.0);

        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= f * _exponents[i];
        }

        return;
    }

    if (_angularMomentum == 3)
    {
        double f = 4.0 / std::sqrt(15.0);

        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= f * _exponents[i] * std::sqrt(_exponents[i]);
        }

        return;
    }

    if (_angularMomentum == 4)
    {
        double f = 2.0 / std::sqrt(105.0);

        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= f * _exponents[i] * _exponents[i];
        }

        return;
    }

    if (_angularMomentum == 5)
    {
        double f = 4.0 / std::sqrt(945.0);

        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= f * _exponents[i] * _exponents[i]

                               * std::sqrt(_exponents[i]);
        }

        return;
    }

    if (_angularMomentum == 6)
    {
        double f = 4.0 / std::sqrt(10395.0);

        for (size_t i = 0; i < _exponents.size(); i++)
        {
            _normFactors[i] *= f * _exponents[i] * _exponents[i]

                               * _exponents[i];
        }

        return;
    }
}

double
CBasisFunction::_overlap(const size_t iComponent, const size_t jComponent) const
{
    auto fab = 1.0 / (_exponents[iComponent] + _exponents[jComponent]);

    auto ovl = _normFactors[iComponent] * _normFactors[jComponent]

               * std::pow(mathconst::getPiValue() * fab, 1.5);

    if (_angularMomentum == 0) return ovl;

    if (_angularMomentum == 1) return 0.5 * fab * ovl;

    auto fab2 = fab * fab;

    if (_angularMomentum == 2) return 3.0 * fab2 * ovl;

    if (_angularMomentum == 3) return 7.5 * fab2 * fab * ovl;

    if (_angularMomentum == 4) return 420.0 * fab2 * fab2 * ovl;

    if (_angularMomentum == 5) return 1890.0 * fab2 * fab2 * fab * ovl;

    if (_angularMomentum == 6) return 41580.0 * fab2 * fab2 * fab2 * ovl;

    return std::numeric_limits<double>::quiet_NaN();
}

void
CBasisFunction::broadcast(int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        mpi::bcast(_angularMomentum, comm);

        mpi::bcast(_exponents, rank, comm);

        mpi::bcast(_normFactors, rank, comm);
    }
}

std::string
CBasisFunction::repr() const
{
    std::ostringstream os;

    os << std::endl;

    os << "[CBasisFunction (Object):" << this << "]" << std::endl;

    os << "_angularMomentum: " << _angularMomentum << std::endl;

    os << "_exponents: " << std::endl;

    for (size_t i = 0; i < _exponents.size(); i++)
    {
        os << "_exponents[" << i << "]: " << std::endl;

        os << _exponents[i] << std::endl;
    }

    os << "_normFactors: " << std::endl;

    for (size_t i = 0; i < _normFactors.size(); i++)
    {
        os << "_normFactors[" << i << "]: " << std::endl;

        os << _normFactors[i] << std::endl;
    }

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CBasisFunction& source)
{
    return (output << source.repr());
}
