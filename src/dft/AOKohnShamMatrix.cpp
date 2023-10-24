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

#include "AOKohnShamMatrix.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

CAOKohnShamMatrix::CAOKohnShamMatrix()

    : _xcRestricted(true)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    _xcMatrices = std::vector<CDenseMatrix>({CDenseMatrix()});
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int64_t nRows, const int64_t nColumns, const bool xcRestricted)

    : _xcRestricted(xcRestricted)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    _xcMatrices = std::vector<CDenseMatrix>({CDenseMatrix(nRows, nColumns)});

    if (!_xcRestricted) _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const CAOKohnShamMatrix& source)

    : _xcMatrices(source._xcMatrices)

    , _xcRestricted(source._xcRestricted)

    , _xcElectrons(source._xcElectrons)

    , _xcEnergy(source._xcEnergy)
{
}

CAOKohnShamMatrix::CAOKohnShamMatrix(CAOKohnShamMatrix&& source) noexcept

    : _xcMatrices(std::move(source._xcMatrices))

    , _xcRestricted(std::move(source._xcRestricted))

    , _xcElectrons(std::move(source._xcElectrons))

    , _xcEnergy(std::move(source._xcEnergy))
{
}

CAOKohnShamMatrix::~CAOKohnShamMatrix()
{
}

auto
CAOKohnShamMatrix::operator=(const CAOKohnShamMatrix& source) -> CAOKohnShamMatrix&
{
    if (this == &source) return *this;

    _xcMatrices = source._xcMatrices;

    _xcRestricted = source._xcRestricted;

    _xcElectrons = source._xcElectrons;

    _xcEnergy = source._xcEnergy;

    return *this;
}

auto
CAOKohnShamMatrix::operator=(CAOKohnShamMatrix&& source) noexcept -> CAOKohnShamMatrix&
{
    if (this == &source) return *this;

    _xcMatrices = std::move(source._xcMatrices);

    _xcRestricted = std::move(source._xcRestricted);

    _xcElectrons = std::move(source._xcElectrons);

    _xcEnergy = std::move(source._xcEnergy);

    return *this;
}

auto
CAOKohnShamMatrix::operator==(const CAOKohnShamMatrix& other) const -> bool
{
    if (_xcMatrices.size() != other._xcMatrices.size()) return false;

    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        if (_xcMatrices[i] != other._xcMatrices[i]) return false;
    }

    if (_xcRestricted != other._xcRestricted) return false;

    if (std::fabs(_xcElectrons - other._xcElectrons) > 1.0e-13) return false;

    if (std::fabs(_xcEnergy - other._xcEnergy) > 1.0e-13) return false;

    return true;
}

auto
CAOKohnShamMatrix::operator!=(const CAOKohnShamMatrix& other) const -> bool
{
    return !(*this == other);
}

auto
CAOKohnShamMatrix::zero() -> void
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].zero();
    }
}

auto
CAOKohnShamMatrix::symmetrizeAndScale(const double factor) -> void
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].symmetrizeAndScale(factor);
    }
}

auto
CAOKohnShamMatrix::setNumberOfElectrons(const double xcElectrons) -> void
{
    _xcElectrons = xcElectrons;
}

auto
CAOKohnShamMatrix::setExchangeCorrelationEnergy(const double xcEnergy) -> void
{
    _xcEnergy = xcEnergy;
}

auto
CAOKohnShamMatrix::isRestricted() const -> bool
{
    return _xcRestricted;
}

auto
CAOKohnShamMatrix::getNumberOfElectrons() const -> double
{
    return _xcElectrons;
}

auto
CAOKohnShamMatrix::getExchangeCorrelationEnergy() const -> double
{
    return _xcEnergy;
}

auto
CAOKohnShamMatrix::getNumberOfRows() const -> int64_t
{
    return _xcMatrices[0].getNumberOfRows();
}

auto
CAOKohnShamMatrix::getNumberOfColumns() const -> int64_t
{
    return _xcMatrices[0].getNumberOfColumns();
}

auto
CAOKohnShamMatrix::getNumberOfElements() const -> int64_t
{
    return _xcMatrices[0].getNumberOfElements();
}

auto
CAOKohnShamMatrix::getReferenceToAlphaKohnSham() const -> const CDenseMatrix&
{
    return _xcMatrices[0];
}

auto
CAOKohnShamMatrix::getReferenceToBetaKohnSham() const -> const CDenseMatrix&
{
    if (_xcRestricted) return _xcMatrices[0];

    return _xcMatrices[1];
}

auto
CAOKohnShamMatrix::getPointerToAlphaValues() const -> const double*
{
    return _xcMatrices[0].values();
}

auto
CAOKohnShamMatrix::getPointerToAlphaValues() -> double*
{
    return _xcMatrices[0].values();
}

auto
CAOKohnShamMatrix::getPointerToBetaValues() const -> const double*
{
    if (_xcRestricted) return _xcMatrices[0].values();

    return _xcMatrices[1].values();
}

auto
CAOKohnShamMatrix::getPointerToBetaValues() -> double*
{
    if (_xcRestricted) return _xcMatrices[0].values();

    return _xcMatrices[1].values();
}
