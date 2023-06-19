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

#include "DensityGrid.hpp"

#include <cmath>
#include <iostream>

CDensityGrid::CDensityGrid()

    : _gridType(dengrid::undefined)

    , _nDensityMatrices(0)

    , _densityValues(CMemBlock2D<double>())
{
}

CDensityGrid::CDensityGrid(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _nDensityMatrices(1)

    , _densityValues(densityValues)
{
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;

    _nDensityMatrices = nDensityMatrices;

    int32_t ncomp = 0;

    if (_gridType == dengrid::ab)
    {
        if (xcFuncType == xcfun::lda) ncomp = 2;

        if (xcFuncType == xcfun::gga) ncomp = 11;

        if (xcFuncType == xcfun::mgga) ncomp = 15;
    }

    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints, const int32_t nDensityMatrices, const int32_t nComponents, const dengrid gridType)
{
    _gridType = gridType;

    _nDensityMatrices = nDensityMatrices;

    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * nComponents);
}

CDensityGrid::CDensityGrid(const CDensityGrid& source)

    : _gridType(source._gridType)

    , _nDensityMatrices(source._nDensityMatrices)

    , _densityValues(source._densityValues)
{
}

CDensityGrid::CDensityGrid(CDensityGrid&& source) noexcept

    : _gridType(std::move(source._gridType))

    , _nDensityMatrices(std::move(source._nDensityMatrices))

    , _densityValues(std::move(source._densityValues))
{
}

CDensityGrid::~CDensityGrid()
{
}

CDensityGrid&
CDensityGrid::operator=(const CDensityGrid& source)
{
    if (this == &source) return *this;

    _gridType = source._gridType;

    _nDensityMatrices = source._nDensityMatrices;

    _densityValues = source._densityValues;

    return *this;
}

CDensityGrid&
CDensityGrid::operator=(CDensityGrid&& source) noexcept
{
    if (this == &source) return *this;

    _gridType = std::move(source._gridType);

    _nDensityMatrices = std::move(source._nDensityMatrices);

    _densityValues = std::move(source._densityValues);

    return *this;
}

bool
CDensityGrid::operator==(const CDensityGrid& other) const
{
    if (_gridType != other._gridType) return false;

    if (_nDensityMatrices != other._nDensityMatrices) return false;

    if (_densityValues != other._densityValues) return false;

    return true;
}

bool
CDensityGrid::operator!=(const CDensityGrid& other) const
{
    return !(*this == other);
}

void
CDensityGrid::zero()
{
    _densityValues.zero();
}

void
CDensityGrid::slice(const int32_t nGridPoints)
{
    if (nGridPoints < getNumberOfGridPoints())
    {
        _densityValues = _densityValues.slice(0, nGridPoints);
    }
}

int32_t
CDensityGrid::getNumberOfGridPoints() const
{
    return _densityValues.size(0);
}

int32_t
CDensityGrid::getNumberOfDensityMatrices() const
{
    return _nDensityMatrices;
}

dengrid
CDensityGrid::getDensityGridType() const
{
    return _gridType;
}

const double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensitytau(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensitytau(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensitytau(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensitytau(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensitylapl(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensitylapl(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensitylapl(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensitylapl(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::getComponent(const int32_t iComponent) const
{
    return _densityValues.data(iComponent);
}

double*
CDensityGrid::getComponent(const int32_t iComponent)
{
    return _densityValues.data(iComponent);
}

std::ostream&
operator<<(std::ostream& output, const CDensityGrid& source)
{
    output << std::endl;

    output << "[CDensityGrid (Object):" << &source << "]" << std::endl;

    output << "_gridType: " << to_string(source._gridType) << std::endl;

    output << "_densityValues: " << std::endl;

    output << source._densityValues << std::endl;

    return output;
}
