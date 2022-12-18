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

#include "DensityGridCubic.hpp"

#include <cmath>
#include <iostream>

#include "StringFormat.hpp"

CDensityGridCubic::CDensityGridCubic()

    : _gridType(dengrid::undefined)

    , _nDensityMatrices(0)

    , _densityValues(CMemBlock2D<double>())
{
}

CDensityGridCubic::CDensityGridCubic(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _nDensityMatrices(1)

    , _densityValues(densityValues)
{
}

CDensityGridCubic::CDensityGridCubic(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;

    _nDensityMatrices = nDensityMatrices;

    int32_t ncomp = 0;

    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 2 : 2;

    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 53 : 53;

    // NOTE: this needs to be checked with mgga functionals implementation

    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 13 : 6;

    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
}

CDensityGridCubic::CDensityGridCubic(const CDensityGridCubic& source)

    : _gridType(source._gridType)

    , _nDensityMatrices(source._nDensityMatrices)

    , _densityValues(source._densityValues)
{
}

CDensityGridCubic::CDensityGridCubic(CDensityGridCubic&& source) noexcept

    : _gridType(std::move(source._gridType))

    , _nDensityMatrices(std::move(source._nDensityMatrices))

    , _densityValues(std::move(source._densityValues))
{
}

CDensityGridCubic::~CDensityGridCubic()
{
}

CDensityGridCubic&
CDensityGridCubic::operator=(const CDensityGridCubic& source)
{
    if (this == &source) return *this;

    _gridType = source._gridType;

    _nDensityMatrices = source._nDensityMatrices;

    _densityValues = source._densityValues;

    return *this;
}

CDensityGridCubic&
CDensityGridCubic::operator=(CDensityGridCubic&& source) noexcept
{
    if (this == &source) return *this;

    _gridType = std::move(source._gridType);

    _nDensityMatrices = std::move(source._nDensityMatrices);

    _densityValues = std::move(source._densityValues);

    return *this;
}

bool
CDensityGridCubic::operator==(const CDensityGridCubic& other) const
{
    if (_gridType != other._gridType) return false;

    if (_nDensityMatrices != other._nDensityMatrices) return false;

    if (_densityValues != other._densityValues) return false;

    return true;
}

bool
CDensityGridCubic::operator!=(const CDensityGridCubic& other) const
{
    return !(*this == other);
}

void
CDensityGridCubic::zero()
{
    _densityValues.zero();
}

int32_t
CDensityGridCubic::getNumberOfGridPoints() const
{
    return _densityValues.size(0);
}

int32_t
CDensityGridCubic::getNumberOfDensityMatrices() const
{
    return _nDensityMatrices;
}

dengrid
CDensityGridCubic::getDensityGridType() const
{
    return _gridType;
}


const double*
CDensityGridCubic::pi(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGridCubic::pi(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGridCubic::gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::gam(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}


const double*
CDensityGridCubic::gamXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::gamX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

void
CDensityGridCubic::DensityProd(const CDensityGrid& rwDensityGrid,
                              const CDensityGrid& rw2DensityGrid,
                              const xcfun         xcFuncType,
                              const int32_t             numdens,
                              const std::string&  CubeMode) 

{
    if (_gridType != dengrid::ab) return;

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    if (xcFuncType == xcfun::lda)
    {
        if (fstr::upcase(CubeMode) == "TPA")
            {
             for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto gamx_r = gam(6 * j);

                auto gamx_i = gam(6 * j + 1);

                auto gamy_r = gam(6 * j + 2);

                auto gamy_i = gam(6 * j + 3);

                auto gamz_r = gam(6 * j + 4);

                auto gamz_i = gam(6 * j + 5);

                auto pix_r = pi(6 * j);

                auto pix_i = pi(6 * j + 1);

                auto piy_r = pi(6 * j + 2 );

                auto piy_i = pi(6 * j + 3);

                auto piz_r = pi(6 * j + 4);

                auto piz_i = pi(6 * j + 5);


                // First-order densities

                auto rhoBx_r = rwDensityGrid.alphaDensity(12 * j);

                auto rhoBx_i = rwDensityGrid.alphaDensity(12 * j + 1);

                auto rhoBy_r = rwDensityGrid.alphaDensity(12 * j + 2);

                auto rhoBy_i = rwDensityGrid.alphaDensity(12 * j + 3);

                auto rhoBz_r = rwDensityGrid.alphaDensity(12 * j + 4);

                auto rhoBz_i = rwDensityGrid.alphaDensity(12 * j + 5);  

                auto rhoCx_r = rwDensityGrid.alphaDensity(12 * j + 6);

                auto rhoCx_i = rwDensityGrid.alphaDensity(12 * j + 7);

                auto rhoCy_r = rwDensityGrid.alphaDensity(12 * j + 8);

                auto rhoCy_i = rwDensityGrid.alphaDensity(12 * j + 9);

                auto rhoCz_r = rwDensityGrid.alphaDensity(12 * j + 10);

                auto rhoCz_i = rwDensityGrid.alphaDensity(12 * j + 11);  


                auto D_sig_xx_r = rw2DensityGrid.alphaDensity(24 * j);

                auto D_sig_xx_i = rw2DensityGrid.alphaDensity(24 * j + 1);

                auto D_sig_yy_r = rw2DensityGrid.alphaDensity(24 * j + 2);

                auto D_sig_yy_i = rw2DensityGrid.alphaDensity(24 * j + 3);

                auto D_sig_zz_r = rw2DensityGrid.alphaDensity(24 * j + 4);

                auto D_sig_zz_i = rw2DensityGrid.alphaDensity(24 * j + 5);

                auto D_sig_xy_r = rw2DensityGrid.alphaDensity(24 * j + 6);

                auto D_sig_xy_i = rw2DensityGrid.alphaDensity(24 * j + 7);

                auto D_sig_xz_r = rw2DensityGrid.alphaDensity(24 * j + 8);

                auto D_sig_xz_i = rw2DensityGrid.alphaDensity(24 * j + 9);

                auto D_sig_yz_r = rw2DensityGrid.alphaDensity(24 * j + 10);

                auto D_sig_yz_i = rw2DensityGrid.alphaDensity(24 * j + 11);

                auto D_lamtau_xx_r = rw2DensityGrid.alphaDensity(24 * j + 12);

                auto D_lamtau_xx_i = rw2DensityGrid.alphaDensity(24 * j + 13);

                auto D_lamtau_yy_r = rw2DensityGrid.alphaDensity(24 * j + 14);

                auto D_lamtau_yy_i = rw2DensityGrid.alphaDensity(24 * j + 15);

                auto D_lamtau_zz_r = rw2DensityGrid.alphaDensity(24 * j + 16);

                auto D_lamtau_zz_i = rw2DensityGrid.alphaDensity(24 * j + 17);

                auto D_lamtau_xy_r = rw2DensityGrid.alphaDensity(24 * j + 18);

                auto D_lamtau_xy_i = rw2DensityGrid.alphaDensity(24 * j + 19);

                auto D_lamtau_xz_r = rw2DensityGrid.alphaDensity(24 * j + 20);

                auto D_lamtau_xz_i = rw2DensityGrid.alphaDensity(24 * j + 21);

                auto D_lamtau_yz_r = rw2DensityGrid.alphaDensity(24 * j + 22);

                auto D_lamtau_yz_i = rw2DensityGrid.alphaDensity(24 * j + 23);


                for (int32_t i = 0; i < npoints; i++)
                {
                    double rhobx_r = rhoBx_r[i]; double rhobx_i = rhoBx_i[i]; double rhoby_r = rhoBy_r[i];

                    double rhoby_i = rhoBy_i[i]; double rhobz_r = rhoBz_r[i]; double rhobz_i = rhoBz_i[i];

                    double rhocx_r = rhoCx_r[i]; double rhocx_i = rhoCx_i[i]; double rhocy_r = rhoCy_r[i];

                    double rhocy_i = rhoCy_i[i]; double rhocz_r = rhoCz_r[i]; double rhocz_i = rhoCz_i[i];


                    double d_sig_xx_r = D_sig_xx_r[i]; double d_sig_xx_i = D_sig_xx_i[i]; double d_sig_yy_r = D_sig_yy_r[i];
                    
                    double d_sig_yy_i = D_sig_yy_i[i]; double d_sig_zz_r = D_sig_zz_r[i]; double d_sig_zz_i = D_sig_zz_i[i];

                    double d_sig_xy_r = D_sig_xy_r[i]; double d_sig_xy_i = D_sig_xy_i[i]; double d_sig_yx_r = D_sig_xy_r[i];

                    double d_sig_yx_i = D_sig_xy_i[i]; double d_sig_xz_r = D_sig_xz_r[i];  double d_sig_xz_i = D_sig_xz_i[i];

                    double d_sig_yz_r = D_sig_yz_r[i]; double d_sig_yz_i = D_sig_yz_i[i]; double d_sig_zx_r = D_sig_xz_r[i];

                    double d_sig_zx_i = D_sig_xz_i[i]; double d_sig_zy_r = D_sig_yz_r[i]; double d_sig_zy_i = D_sig_yz_i[i];

                    
                    double d_lamtau_xx_r = D_lamtau_xx_r[i]; double d_lamtau_xx_i = D_lamtau_xx_i[i];double d_lamtau_yy_r = D_lamtau_yy_r[i];
                    
                    double d_lamtau_yy_i = D_lamtau_yy_i[i]; double d_lamtau_zz_r = D_lamtau_zz_r[i]; double d_lamtau_zz_i = D_lamtau_zz_i[i];

                    double d_lamtau_xy_r = D_lamtau_xy_r[i]; double d_lamtau_xy_i = D_lamtau_xy_i[i]; double d_lamtau_xz_r = D_lamtau_xz_r[i];

                    double d_lamtau_xz_i = D_lamtau_xz_i[i]; double d_lamtau_yz_r = D_lamtau_yz_r[i]; double d_lamtau_yz_i = D_lamtau_yz_i[i];


                    double d_lamtau_yx_r = D_lamtau_xy_r[i]; double d_lamtau_yx_i = D_lamtau_xy_i[i]; double d_lamtau_zx_r = D_lamtau_xz_r[i];

                    double d_lamtau_zx_i = D_lamtau_xz_i[i]; double d_lamtau_zy_r = D_lamtau_yz_r[i]; double d_lamtau_zy_i = D_lamtau_yz_i[i];


                    pix_r[i] = 12.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 12.0 * (rhobx_r * rhocy_r * rhoby_r - rhobx_i * rhocy_r * rhoby_i - rhobx_r * rhocy_i * rhoby_i - rhobx_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocx_r * rhoby_r - rhoby_i * rhocx_r * rhoby_i - rhoby_r * rhocx_i * rhoby_i - rhoby_i * rhocx_i * rhoby_r)

                    + 12.0 * (rhobx_r * rhocz_r * rhobz_r - rhobx_i * rhocz_r * rhobz_i - rhobx_r * rhocz_i * rhobz_i - rhobx_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocx_r * rhobz_r - rhobz_i * rhocx_r * rhobz_i - rhobz_r * rhocx_i * rhobz_i - rhobz_i * rhocx_i * rhobz_r);


                    piy_r[i] = 12.0 * (rhoby_r * rhocx_r * rhobx_r - rhoby_i * rhocx_r * rhobx_i - rhoby_r * rhocx_i * rhobx_i - rhoby_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocy_r * rhobx_r - rhobx_i * rhocy_r * rhobx_i - rhobx_r * rhocy_i * rhobx_i - rhobx_i * rhocy_i * rhobx_r)

                    + 12.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 12.0 * (rhoby_r * rhocz_r * rhobz_r - rhoby_i * rhocz_r * rhobz_i - rhoby_r * rhocz_i * rhobz_i - rhoby_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocy_r * rhobz_r - rhobz_i * rhocy_r * rhobz_i - rhobz_r * rhocy_i * rhobz_i - rhobz_i * rhocy_i * rhobz_r);


                    piz_r[i] = 12.0 * (rhobz_r * rhocx_r * rhobx_r - rhobz_i * rhocx_r * rhobx_i - rhobz_r * rhocx_i * rhobx_i - rhobz_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocz_r * rhobx_r - rhobx_i * rhocz_r * rhobx_i - rhobx_r * rhocz_i * rhobx_i - rhobx_i * rhocz_i * rhobx_r)

                    + 12.0 * (rhobz_r * rhocy_r * rhoby_r - rhobz_i * rhocy_r * rhoby_i - rhobz_r * rhocy_i * rhoby_i - rhobz_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocz_r * rhoby_r - rhoby_i * rhocz_r * rhoby_i - rhoby_r * rhocz_i * rhoby_i - rhoby_i * rhocz_i * rhoby_r)

                    + 12.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r);



                    pix_i[i] = 12.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 12.0 * ( - rhobx_i * rhocy_i * rhoby_i + rhobx_i * rhocy_r * rhoby_r + rhobx_r * rhocy_i * rhoby_r + rhobx_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocx_i * rhoby_i + rhoby_i * rhocx_r * rhoby_r + rhoby_r * rhocx_i * rhoby_r + rhoby_r * rhocx_r * rhoby_i)

                    + 12.0 * ( - rhobx_i * rhocz_i * rhobz_i + rhobx_i * rhocz_r * rhobz_r + rhobx_r * rhocz_i * rhobz_r + rhobx_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocx_i * rhobz_i + rhobz_i * rhocx_r * rhobz_r + rhobz_r * rhocx_i * rhobz_r + rhobz_r * rhocx_r * rhobz_i);


                    piy_i[i] = 12.0 * ( - rhoby_i * rhocx_i * rhobx_i + rhoby_i * rhocx_r * rhobx_r + rhoby_r * rhocx_i * rhobx_r + rhoby_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocy_i * rhobx_i + rhobx_i * rhocy_r * rhobx_r + rhobx_r * rhocy_i * rhobx_r + rhobx_r * rhocy_r * rhobx_i)

                    + 12.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 12.0 * ( - rhoby_i * rhocz_i * rhobz_i + rhoby_i * rhocz_r * rhobz_r + rhoby_r * rhocz_i * rhobz_r + rhoby_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocy_i * rhobz_i + rhobz_i * rhocy_r * rhobz_r + rhobz_r * rhocy_i * rhobz_r + rhobz_r * rhocy_r * rhobz_i);


                    piz_i[i] = 12.0 * ( - rhobz_i * rhocx_i * rhobx_i + rhobz_i * rhocx_r * rhobx_r + rhobz_r * rhocx_i * rhobx_r + rhobz_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocz_i * rhobx_i + rhobx_i * rhocz_r * rhobx_r + rhobx_r * rhocz_i * rhobx_r + rhobx_r * rhocz_r * rhobx_i)

                    + 12.0 * ( - rhobz_i * rhocy_i * rhoby_i + rhobz_i * rhocy_r * rhoby_r + rhobz_r * rhocy_i * rhoby_r + rhobz_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocz_i * rhoby_i + rhoby_i * rhocz_r * rhoby_r + rhoby_r * rhocz_i * rhoby_r + rhoby_r * rhocz_r * rhoby_i)

                    + 12.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i);


gamx_r[i] = (d_sig_xx_r * rhobx_r- d_sig_xx_i * rhobx_i);

gamx_i[i] =  (d_sig_xx_i * rhobx_r + d_sig_xx_r * rhobx_i);

gamx_r[i] += (d_lamtau_xx_r * rhocx_r- d_lamtau_xx_i * rhocx_i);

gamx_i[i] +=  (d_lamtau_xx_i * rhocx_r + d_lamtau_xx_r * rhocx_i);

gamx_r[i] += (d_sig_xy_r * rhoby_r- d_sig_xy_i * rhoby_i);

gamx_i[i] +=  (d_sig_xy_i * rhoby_r + d_sig_xy_r * rhoby_i);

gamx_r[i] += (d_lamtau_xy_r * rhocy_r- d_lamtau_xy_i * rhocy_i);

gamx_i[i] +=  (d_lamtau_xy_i * rhocy_r + d_lamtau_xy_r * rhocy_i);

gamx_r[i] += (d_sig_xz_r * rhobz_r- d_sig_xz_i * rhobz_i);

gamx_i[i] +=  (d_sig_xz_i * rhobz_r + d_sig_xz_r * rhobz_i);

gamx_r[i] += (d_lamtau_xz_r * rhocz_r- d_lamtau_xz_i * rhocz_i);

gamx_i[i] +=  (d_lamtau_xz_i * rhocz_r + d_lamtau_xz_r * rhocz_i);

gamy_r[i] = (d_sig_yx_r * rhobx_r- d_sig_yx_i * rhobx_i);

gamy_i[i] =  (d_sig_yx_i * rhobx_r + d_sig_yx_r * rhobx_i);

gamy_r[i] += (d_lamtau_yx_r * rhocx_r- d_lamtau_yx_i * rhocx_i);

gamy_i[i] +=  (d_lamtau_yx_i * rhocx_r + d_lamtau_yx_r * rhocx_i);

gamy_r[i] += (d_sig_yy_r * rhoby_r- d_sig_yy_i * rhoby_i);

gamy_i[i] +=  (d_sig_yy_i * rhoby_r + d_sig_yy_r * rhoby_i);

gamy_r[i] += (d_lamtau_yy_r * rhocy_r- d_lamtau_yy_i * rhocy_i);

gamy_i[i] +=  (d_lamtau_yy_i * rhocy_r + d_lamtau_yy_r * rhocy_i);

gamy_r[i] += (d_sig_yz_r * rhobz_r- d_sig_yz_i * rhobz_i);

gamy_i[i] +=  (d_sig_yz_i * rhobz_r + d_sig_yz_r * rhobz_i);

gamy_r[i] += (d_lamtau_yz_r * rhocz_r- d_lamtau_yz_i * rhocz_i);

gamy_i[i] +=  (d_lamtau_yz_i * rhocz_r + d_lamtau_yz_r * rhocz_i);

gamz_r[i] = (d_sig_zx_r * rhobx_r- d_sig_zx_i * rhobx_i);

gamz_i[i] =  (d_sig_zx_i * rhobx_r + d_sig_zx_r * rhobx_i);

gamz_r[i] += (d_lamtau_zx_r * rhocx_r- d_lamtau_zx_i * rhocx_i);

gamz_i[i] +=  (d_lamtau_zx_i * rhocx_r + d_lamtau_zx_r * rhocx_i);

gamz_r[i] += (d_sig_zy_r * rhoby_r- d_sig_zy_i * rhoby_i);

gamz_i[i] +=  (d_sig_zy_i * rhoby_r + d_sig_zy_r * rhoby_i);

gamz_r[i] += (d_lamtau_zy_r * rhocy_r- d_lamtau_zy_i * rhocy_i);

gamz_i[i] +=  (d_lamtau_zy_i * rhocy_r + d_lamtau_zy_r * rhocy_i);

gamz_r[i] += (d_sig_zz_r * rhobz_r- d_sig_zz_i * rhobz_i);

gamz_i[i] +=  (d_sig_zz_i * rhobz_r + d_sig_zz_r * rhobz_i);

gamz_r[i] += (d_lamtau_zz_r * rhocz_r- d_lamtau_zz_i * rhocz_i);

gamz_i[i] +=  (d_lamtau_zz_i * rhocz_r + d_lamtau_zz_r * rhocz_i);
                }

            }

        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                
                auto gam_r = gam(2 * j);

                auto gam_i = gam(2 * j + 1);

                auto pi_r = pi(2 * j);

                auto pi_i = pi(2 * j + 1);

                
                // First-order densities

                auto rhoB_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhoB_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhoC_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhoC_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhoD_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhoD_i = rwDensityGrid.alphaDensity(6 * j + 5);


                auto rhoBC_r = rw2DensityGrid.alphaDensity(6 * j);

                auto rhoBC_i = rw2DensityGrid.alphaDensity(6 * j + 1);

                auto rhoBD_r = rw2DensityGrid.alphaDensity(6 * j + 2);

                auto rhoBD_i = rw2DensityGrid.alphaDensity(6 * j + 3);

                auto rhoCD_r = rw2DensityGrid.alphaDensity(6 * j + 4);

                auto rhoCD_i = rw2DensityGrid.alphaDensity(6 * j + 5);

                // pi = rho_b * rho_c * rho_d

                for (int32_t i = 0; i < npoints; i++)
                {
                    pi_r[i] = 6.0 * (rhoB_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * rhoD_r[i]);

                    pi_i[i] = 6.0 *  (- rhoB_i[i] * rhoC_i[i] * rhoD_i[i] 
                    
                                        + rhoB_i[i] * rhoC_r[i] * rhoD_r[i] 
                                        
                                        + rhoC_i[i] * rhoB_r[i] * rhoD_r[i] 
                                        
                                        + rhoD_i[i] * rhoB_r[i] * rhoC_r[i]);

                    gam_r[i] = 3.0 * (rhoD_r[i] * rhoBC_r[i]- rhoD_i[i] * rhoBC_i[i]);

                    gam_i[i] = 3.0 * (rhoD_i[i] * rhoBC_r[i] + rhoD_r[i] * rhoBC_i[i]);

                    gam_r[i] += 3.0 * (rhoB_r[i] * rhoCD_r[i]- rhoB_i[i] * rhoCD_i[i]);

                    gam_i[i] += 3.0 * (rhoB_i[i] * rhoCD_r[i] + rhoB_r[i] * rhoCD_i[i]);

                    gam_r[i] += 3.0 * (rhoC_r[i] * rhoBD_r[i]- rhoC_i[i] * rhoBD_i[i]);

                    gam_i[i] += 3.0 * (rhoC_i[i] * rhoBD_r[i] + rhoC_r[i] * rhoBD_i[i]);

                }
            }
        }

    }
    if (xcFuncType == xcfun::gga)
    {
        
        if (fstr::upcase(CubeMode) == "TPA")
        {
            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto gamx_r = gam(6 * j);
                auto gamx_X_r = gamX(6 * j);
                auto gamx_Y_r = gamY(6 * j);
                auto gamx_Z_r = gamZ(6 * j);
                auto gamx_XX_r = gamXX(6 * j);
                auto gamx_XY_r = gamXY(6 * j);
                auto gamx_XZ_r = gamXZ(6 * j);
                auto gamx_YX_r = gamYX(6 * j);
                auto gamx_YY_r = gamYY(6 * j);
                auto gamx_YZ_r = gamYZ(6 * j);
                auto gamx_ZX_r = gamZX(6 * j);
                auto gamx_ZY_r = gamZY(6 * j);
                auto gamx_ZZ_r = gamZZ(6 * j);

                auto gamx_i = gam(6 * j + 1);
                auto gamx_X_i = gamX(6 * j + 1);
                auto gamx_Y_i = gamY(6 * j + 1);
                auto gamx_Z_i = gamZ(6 * j + 1);
                auto gamx_XX_i = gamXX(6 * j + 1);
                auto gamx_XY_i = gamXY(6 * j + 1);
                auto gamx_XZ_i = gamXZ(6 * j + 1);
                auto gamx_YX_i = gamYX(6 * j + 1);
                auto gamx_YY_i = gamYY(6 * j + 1);
                auto gamx_YZ_i = gamYZ(6 * j + 1);
                auto gamx_ZX_i = gamZX(6 * j + 1);
                auto gamx_ZY_i = gamZY(6 * j + 1);
                auto gamx_ZZ_i = gamZZ(6 * j + 1);

                auto gamy_r = gam(6 * j + 2);
                auto gamy_X_r = gamX(6 * j + 2);
                auto gamy_Y_r = gamY(6 * j + 2);
                auto gamy_Z_r = gamZ(6 * j + 2);
                auto gamy_XX_r = gamXX(6 * j + 2);
                auto gamy_XY_r = gamXY(6 * j + 2);
                auto gamy_XZ_r = gamXZ(6 * j + 2);
                auto gamy_YX_r = gamYX(6 * j + 2);
                auto gamy_YY_r = gamYY(6 * j + 2);
                auto gamy_YZ_r = gamYZ(6 * j + 2);
                auto gamy_ZX_r = gamZX(6 * j + 2);
                auto gamy_ZY_r = gamZY(6 * j + 2);
                auto gamy_ZZ_r = gamZZ(6 * j + 2);

                auto gamy_i = gam(6 * j + 3);
                auto gamy_X_i = gamX(6 * j + 3);
                auto gamy_Y_i = gamY(6 * j + 3);
                auto gamy_Z_i = gamZ(6 * j + 3);
                auto gamy_XX_i = gamXX(6 * j + 3);
                auto gamy_XY_i = gamXY(6 * j + 3);
                auto gamy_XZ_i = gamXZ(6 * j + 3);
                auto gamy_YX_i = gamYX(6 * j + 3);
                auto gamy_YY_i = gamYY(6 * j + 3);
                auto gamy_YZ_i = gamYZ(6 * j + 3);
                auto gamy_ZX_i = gamZX(6 * j + 3);
                auto gamy_ZY_i = gamZY(6 * j + 3);
                auto gamy_ZZ_i = gamZZ(6 * j + 3);

                auto gamz_r = gam(6 * j + 4);
                auto gamz_X_r = gamX(6 * j + 4);
                auto gamz_Y_r = gamY(6 * j + 4);
                auto gamz_Z_r = gamZ(6 * j + 4);
                auto gamz_XX_r = gamXX(6 * j + 4);
                auto gamz_XY_r = gamXY(6 * j + 4);
                auto gamz_XZ_r = gamXZ(6 * j + 4);
                auto gamz_YX_r = gamYX(6 * j + 4);
                auto gamz_YY_r = gamYY(6 * j + 4);
                auto gamz_YZ_r = gamYZ(6 * j + 4);
                auto gamz_ZX_r = gamZX(6 * j + 4);
                auto gamz_ZY_r = gamZY(6 * j + 4);
                auto gamz_ZZ_r = gamZZ(6 * j + 4);

                auto gamz_i = gam(6 * j + 5);
                auto gamz_X_i = gamX(6 * j + 5);
                auto gamz_Y_i = gamY(6 * j + 5);
                auto gamz_Z_i = gamZ(6 * j + 5);
                auto gamz_XX_i = gamXX(6 * j + 5);
                auto gamz_XY_i = gamXY(6 * j + 5);
                auto gamz_XZ_i = gamXZ(6 * j + 5);
                auto gamz_YX_i = gamYX(6 * j + 5);
                auto gamz_YY_i = gamYY(6 * j + 5);
                auto gamz_YZ_i = gamYZ(6 * j + 5);
                auto gamz_ZX_i = gamZX(6 * j + 5);
                auto gamz_ZY_i = gamZY(6 * j + 5);
                auto gamz_ZZ_i = gamZZ(6 * j + 5);
                
                auto pix_r = pi(6 * j);
                auto pix_i = pi(6 * j + 1);
                auto piy_r = pi(6 * j + 2);
                auto piy_i = pi(6 * j + 3);
                auto piz_r = pi(6 * j + 4);
                auto piz_i = pi(6 * j + 5);

                auto pix_X_r = piX(6 * j);
                auto pix_Y_r = piY(6 * j);
                auto pix_Z_r = piZ(6 * j);

                auto pix_XX_r = piXX(6 * j);
                auto pix_XY_r = piXY(6 * j);
                auto pix_XZ_r = piXZ(6 * j);
                auto pix_YX_r = piYX(6 * j);
                auto pix_YY_r = piYY(6 * j);
                auto pix_YZ_r = piYZ(6 * j);
                auto pix_ZX_r = piZX(6 * j);
                auto pix_ZY_r = piZY(6 * j);
                auto pix_ZZ_r = piZZ(6 * j);

                auto pix_XXX_r = piXXX(6 * j);
                auto pix_XXY_r = piXXY(6 * j);
                auto pix_XXZ_r = piXXZ(6 * j);
                auto pix_XYX_r = piXYX(6 * j);
                auto pix_XYY_r = piXYY(6 * j);
                auto pix_XYZ_r = piXYZ(6 * j);
                auto pix_XZX_r = piXZX(6 * j);
                auto pix_XZY_r = piXZY(6 * j);
                auto pix_XZZ_r = piXZZ(6 * j);
                auto pix_YXX_r = piYXX(6 * j);
                auto pix_YXY_r = piYXY(6 * j);
                auto pix_YXZ_r = piYXZ(6 * j);
                auto pix_YYX_r = piYYX(6 * j);
                auto pix_YYY_r = piYYY(6 * j);
                auto pix_YYZ_r = piYYZ(6 * j);
                auto pix_YZX_r = piYZX(6 * j);
                auto pix_YZY_r = piYZY(6 * j);
                auto pix_YZZ_r = piYZZ(6 * j);
                auto pix_ZXX_r = piZXX(6 * j);
                auto pix_ZXY_r = piZXY(6 * j);
                auto pix_ZXZ_r = piZXZ(6 * j);
                auto pix_ZYX_r = piZYX(6 * j);
                auto pix_ZYY_r = piZYY(6 * j);
                auto pix_ZYZ_r = piZYZ(6 * j);
                auto pix_ZZX_r = piZZX(6 * j);
                auto pix_ZZY_r = piZZY(6 * j);
                auto pix_ZZZ_r = piZZZ(6 * j);

                auto pix_X_i = piX(6 * j + 1);
                auto pix_Y_i = piY(6 * j + 1);
                auto pix_Z_i = piZ(6 * j + 1);

                auto pix_XX_i = piXX(6 * j + 1);
                auto pix_XY_i = piXY(6 * j + 1);
                auto pix_XZ_i = piXZ(6 * j + 1);
                auto pix_YX_i = piYX(6 * j + 1);
                auto pix_YY_i = piYY(6 * j + 1);
                auto pix_YZ_i = piYZ(6 * j + 1);
                auto pix_ZX_i = piZX(6 * j + 1);
                auto pix_ZY_i = piZY(6 * j + 1);
                auto pix_ZZ_i = piZZ(6 * j + 1);

                auto pix_XXX_i = piXXX(6 * j + 1);
                auto pix_XXY_i = piXXY(6 * j + 1);
                auto pix_XXZ_i = piXXZ(6 * j + 1);
                auto pix_XYX_i = piXYX(6 * j + 1);
                auto pix_XYY_i = piXYY(6 * j + 1);
                auto pix_XYZ_i = piXYZ(6 * j + 1);
                auto pix_XZX_i = piXZX(6 * j + 1);
                auto pix_XZY_i = piXZY(6 * j + 1);
                auto pix_XZZ_i = piXZZ(6 * j + 1);
                auto pix_YXX_i = piYXX(6 * j + 1);
                auto pix_YXY_i = piYXY(6 * j + 1);
                auto pix_YXZ_i = piYXZ(6 * j + 1);
                auto pix_YYX_i = piYYX(6 * j + 1);
                auto pix_YYY_i = piYYY(6 * j + 1);
                auto pix_YYZ_i = piYYZ(6 * j + 1);
                auto pix_YZX_i = piYZX(6 * j + 1);
                auto pix_YZY_i = piYZY(6 * j + 1);
                auto pix_YZZ_i = piYZZ(6 * j + 1);
                auto pix_ZXX_i = piZXX(6 * j + 1);
                auto pix_ZXY_i = piZXY(6 * j + 1);
                auto pix_ZXZ_i = piZXZ(6 * j + 1);
                auto pix_ZYX_i = piZYX(6 * j + 1);
                auto pix_ZYY_i = piZYY(6 * j + 1);
                auto pix_ZYZ_i = piZYZ(6 * j + 1);
                auto pix_ZZX_i = piZZX(6 * j + 1);
                auto pix_ZZY_i = piZZY(6 * j + 1);
                auto pix_ZZZ_i = piZZZ(6 * j + 1);

                auto piy_X_r = piX(6 * j + 2);
                auto piy_Y_r = piY(6 * j + 2);
                auto piy_Z_r = piZ(6 * j + 2);

                auto piy_XX_r = piXX(6 * j + 2);
                auto piy_XY_r = piXY(6 * j + 2);
                auto piy_XZ_r = piXZ(6 * j + 2);
                auto piy_YX_r = piYX(6 * j + 2);
                auto piy_YY_r = piYY(6 * j + 2);
                auto piy_YZ_r = piYZ(6 * j + 2);
                auto piy_ZX_r = piZX(6 * j + 2);
                auto piy_ZY_r = piZY(6 * j + 2);
                auto piy_ZZ_r = piZZ(6 * j + 2);

                auto piy_XXX_r = piXXX(6 * j + 2);
                auto piy_XXY_r = piXXY(6 * j + 2);
                auto piy_XXZ_r = piXXZ(6 * j + 2);
                auto piy_XYX_r = piXYX(6 * j + 2);
                auto piy_XYY_r = piXYY(6 * j + 2);
                auto piy_XYZ_r = piXYZ(6 * j + 2);
                auto piy_XZX_r = piXZX(6 * j + 2);
                auto piy_XZY_r = piXZY(6 * j + 2);
                auto piy_XZZ_r = piXZZ(6 * j + 2);
                auto piy_YXX_r = piYXX(6 * j + 2);
                auto piy_YXY_r = piYXY(6 * j + 2);
                auto piy_YXZ_r = piYXZ(6 * j + 2);
                auto piy_YYX_r = piYYX(6 * j + 2);
                auto piy_YYY_r = piYYY(6 * j + 2);
                auto piy_YYZ_r = piYYZ(6 * j + 2);
                auto piy_YZX_r = piYZX(6 * j + 2);
                auto piy_YZY_r = piYZY(6 * j + 2);
                auto piy_YZZ_r = piYZZ(6 * j + 2);
                auto piy_ZXX_r = piZXX(6 * j + 2);
                auto piy_ZXY_r = piZXY(6 * j + 2);
                auto piy_ZXZ_r = piZXZ(6 * j + 2);
                auto piy_ZYX_r = piZYX(6 * j + 2);
                auto piy_ZYY_r = piZYY(6 * j + 2);
                auto piy_ZYZ_r = piZYZ(6 * j + 2);
                auto piy_ZZX_r = piZZX(6 * j + 2);
                auto piy_ZZY_r = piZZY(6 * j + 2);
                auto piy_ZZZ_r = piZZZ(6 * j + 2);

                auto piy_X_i = piX(6 * j + 3);
                auto piy_Y_i = piY(6 * j + 3);
                auto piy_Z_i = piZ(6 * j + 3);

                auto piy_XX_i = piXX(6 * j + 3);
                auto piy_XY_i = piXY(6 * j + 3);
                auto piy_XZ_i = piXZ(6 * j + 3);
                auto piy_YX_i = piYX(6 * j + 3);
                auto piy_YY_i = piYY(6 * j + 3);
                auto piy_YZ_i = piYZ(6 * j + 3);
                auto piy_ZX_i = piZX(6 * j + 3);
                auto piy_ZY_i = piZY(6 * j + 3);
                auto piy_ZZ_i = piZZ(6 * j + 3);

                auto piy_XXX_i = piXXX(6 * j + 3);
                auto piy_XXY_i = piXXY(6 * j + 3);
                auto piy_XXZ_i = piXXZ(6 * j + 3);
                auto piy_XYX_i = piXYX(6 * j + 3);
                auto piy_XYY_i = piXYY(6 * j + 3);
                auto piy_XYZ_i = piXYZ(6 * j + 3);
                auto piy_XZX_i = piXZX(6 * j + 3);
                auto piy_XZY_i = piXZY(6 * j + 3);
                auto piy_XZZ_i = piXZZ(6 * j + 3);
                auto piy_YXX_i = piYXX(6 * j + 3);
                auto piy_YXY_i = piYXY(6 * j + 3);
                auto piy_YXZ_i = piYXZ(6 * j + 3);
                auto piy_YYX_i = piYYX(6 * j + 3);
                auto piy_YYY_i = piYYY(6 * j + 3);
                auto piy_YYZ_i = piYYZ(6 * j + 3);
                auto piy_YZX_i = piYZX(6 * j + 3);
                auto piy_YZY_i = piYZY(6 * j + 3);
                auto piy_YZZ_i = piYZZ(6 * j + 3);
                auto piy_ZXX_i = piZXX(6 * j + 3);
                auto piy_ZXY_i = piZXY(6 * j + 3);
                auto piy_ZXZ_i = piZXZ(6 * j + 3);
                auto piy_ZYX_i = piZYX(6 * j + 3);
                auto piy_ZYY_i = piZYY(6 * j + 3);
                auto piy_ZYZ_i = piZYZ(6 * j + 3);
                auto piy_ZZX_i = piZZX(6 * j + 3);
                auto piy_ZZY_i = piZZY(6 * j + 3);
                auto piy_ZZZ_i = piZZZ(6 * j + 3);

                auto piz_X_r = piX(6 * j + 4);
                auto piz_Y_r = piY(6 * j + 4);
                auto piz_Z_r = piZ(6 * j + 4);

                auto piz_XX_r = piXX(6 * j + 4);
                auto piz_XY_r = piXY(6 * j + 4);
                auto piz_XZ_r = piXZ(6 * j + 4);
                auto piz_YX_r = piYX(6 * j + 4);
                auto piz_YY_r = piYY(6 * j + 4);
                auto piz_YZ_r = piYZ(6 * j + 4);
                auto piz_ZX_r = piZX(6 * j + 4);
                auto piz_ZY_r = piZY(6 * j + 4);
                auto piz_ZZ_r = piZZ(6 * j + 4);

                auto piz_XXX_r = piXXX(6 * j + 4);
                auto piz_XXY_r = piXXY(6 * j + 4);
                auto piz_XXZ_r = piXXZ(6 * j + 4);
                auto piz_XYX_r = piXYX(6 * j + 4);
                auto piz_XYY_r = piXYY(6 * j + 4);
                auto piz_XYZ_r = piXYZ(6 * j + 4);
                auto piz_XZX_r = piXZX(6 * j + 4);
                auto piz_XZY_r = piXZY(6 * j + 4);
                auto piz_XZZ_r = piXZZ(6 * j + 4);
                auto piz_YXX_r = piYXX(6 * j + 4);
                auto piz_YXY_r = piYXY(6 * j + 4);
                auto piz_YXZ_r = piYXZ(6 * j + 4);
                auto piz_YYX_r = piYYX(6 * j + 4);
                auto piz_YYY_r = piYYY(6 * j + 4);
                auto piz_YYZ_r = piYYZ(6 * j + 4);
                auto piz_YZX_r = piYZX(6 * j + 4);
                auto piz_YZY_r = piYZY(6 * j + 4);
                auto piz_YZZ_r = piYZZ(6 * j + 4);
                auto piz_ZXX_r = piZXX(6 * j + 4);
                auto piz_ZXY_r = piZXY(6 * j + 4);
                auto piz_ZXZ_r = piZXZ(6 * j + 4);
                auto piz_ZYX_r = piZYX(6 * j + 4);
                auto piz_ZYY_r = piZYY(6 * j + 4);
                auto piz_ZYZ_r = piZYZ(6 * j + 4);
                auto piz_ZZX_r = piZZX(6 * j + 4);
                auto piz_ZZY_r = piZZY(6 * j + 4);
                auto piz_ZZZ_r = piZZZ(6 * j + 4);


                auto piz_X_i = piX(6 * j + 5);
                auto piz_Y_i = piY(6 * j + 5);
                auto piz_Z_i = piZ(6 * j + 5);

                auto piz_XX_i = piXX(6 * j + 5);
                auto piz_XY_i = piXY(6 * j + 5);
                auto piz_XZ_i = piXZ(6 * j + 5);
                auto piz_YX_i = piYX(6 * j + 5);
                auto piz_YY_i = piYY(6 * j + 5);
                auto piz_YZ_i = piYZ(6 * j + 5);
                auto piz_ZX_i = piZX(6 * j + 5);
                auto piz_ZY_i = piZY(6 * j + 5);
                auto piz_ZZ_i = piZZ(6 * j + 5);

                auto piz_XXX_i = piXXX(6 * j + 5);
                auto piz_XXY_i = piXXY(6 * j + 5);
                auto piz_XXZ_i = piXXZ(6 * j + 5);
                auto piz_XYX_i = piXYX(6 * j + 5);
                auto piz_XYY_i = piXYY(6 * j + 5);
                auto piz_XYZ_i = piXYZ(6 * j + 5);
                auto piz_XZX_i = piXZX(6 * j + 5);
                auto piz_XZY_i = piXZY(6 * j + 5);
                auto piz_XZZ_i = piXZZ(6 * j + 5);
                auto piz_YXX_i = piYXX(6 * j + 5);
                auto piz_YXY_i = piYXY(6 * j + 5);
                auto piz_YXZ_i = piYXZ(6 * j + 5);
                auto piz_YYX_i = piYYX(6 * j + 5);
                auto piz_YYY_i = piYYY(6 * j + 5);
                auto piz_YYZ_i = piYYZ(6 * j + 5);
                auto piz_YZX_i = piYZX(6 * j + 5);
                auto piz_YZY_i = piYZY(6 * j + 5);
                auto piz_YZZ_i = piYZZ(6 * j + 5);
                auto piz_ZXX_i = piZXX(6 * j + 5);
                auto piz_ZXY_i = piZXY(6 * j + 5);
                auto piz_ZXZ_i = piZXZ(6 * j + 5);
                auto piz_ZYX_i = piZYX(6 * j + 5);
                auto piz_ZYY_i = piZYY(6 * j + 5);
                auto piz_ZYZ_i = piZYZ(6 * j + 5);
                auto piz_ZZX_i = piZZX(6 * j + 5);
                auto piz_ZZY_i = piZZY(6 * j + 5);
                auto piz_ZZZ_i = piZZZ(6 * j + 5);


                // First-order densities First variable is for the component of the dipole, the second variable is for the component of the gradient
                // B means positive frequency and C means negative frequency

                auto rhoBx_r = rwDensityGrid.alphaDensity(12 * j);
                auto gradBx_x_r = rwDensityGrid.alphaDensityGradientX(12 * j);
                auto gradBx_y_r = rwDensityGrid.alphaDensityGradientY(12 * j);
                auto gradBx_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j);

                
                auto rhoBx_i = rwDensityGrid.alphaDensity(12 * j + 1);
                auto gradBx_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 1);
                auto gradBx_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 1);
                auto gradBx_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 1);


                auto rhoBy_r = rwDensityGrid.alphaDensity(12 * j + 2);
                auto gradBy_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 2);
                auto gradBy_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 2);
                auto gradBy_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 2);


                auto rhoBy_i = rwDensityGrid.alphaDensity(12 * j + 3);
                auto gradBy_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 3);
                auto gradBy_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 3);
                auto gradBy_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 3);


                auto rhoBz_r = rwDensityGrid.alphaDensity(12 * j + 4);
                auto gradBz_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 4);
                auto gradBz_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 4);
                auto gradBz_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 4);


                auto rhoBz_i = rwDensityGrid.alphaDensity(12 * j + 5);  
                auto gradBz_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 5);
                auto gradBz_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 5);
                auto gradBz_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 5);


                auto rhoCx_r = rwDensityGrid.alphaDensity(12 * j + 6);
                auto gradCx_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 6);
                auto gradCx_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 6);
                auto gradCx_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 6);


                auto rhoCx_i = rwDensityGrid.alphaDensity(12 * j + 7);
                auto gradCx_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 7);
                auto gradCx_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 7);
                auto gradCx_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 7);


                auto rhoCy_r = rwDensityGrid.alphaDensity(12 * j + 8);
                auto gradCy_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 8);
                auto gradCy_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 8);
                auto gradCy_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 8);


                auto rhoCy_i = rwDensityGrid.alphaDensity(12 * j + 9);
                auto gradCy_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 9);
                auto gradCy_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 9);
                auto gradCy_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 9);


                auto rhoCz_r = rwDensityGrid.alphaDensity(12 * j + 10);
                auto gradCz_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 10);
                auto gradCz_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 10);
                auto gradCz_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 10);
                

                auto rhoCz_i = rwDensityGrid.alphaDensity(12 * j + 11);  
                auto gradCz_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 11);
                auto gradCz_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 11);
                auto gradCz_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 11);


                auto D_sig_xx_r = rw2DensityGrid.alphaDensity(24 * j);
                auto D_sig_xx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j);
                auto D_sig_xx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j);
                auto D_sig_xx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j);


                auto D_sig_xx_i = rw2DensityGrid.alphaDensity(24 * j + 1);
                auto D_sig_xx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 1);
                auto D_sig_xx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 1);
                auto D_sig_xx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 1);


                auto D_sig_yy_r = rw2DensityGrid.alphaDensity(24 * j + 2);
                auto D_sig_yy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 2);
                auto D_sig_yy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 2);
                auto D_sig_yy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 2);


                auto D_sig_yy_i = rw2DensityGrid.alphaDensity(24 * j + 3);
                auto D_sig_yy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 3);
                auto D_sig_yy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 3);
                auto D_sig_yy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 3);


                auto D_sig_zz_r = rw2DensityGrid.alphaDensity(24 * j + 4);
                auto D_sig_zz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 4);
                auto D_sig_zz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 4);
                auto D_sig_zz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 4);


                auto D_sig_zz_i = rw2DensityGrid.alphaDensity(24 * j + 5);
                auto D_sig_zz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 5);
                auto D_sig_zz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 5);
                auto D_sig_zz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 5);


                auto D_sig_xy_r = rw2DensityGrid.alphaDensity(24 * j + 6);
                auto D_sig_xy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 6);
                auto D_sig_xy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 6);
                auto D_sig_xy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 6);

                auto D_sig_xy_i = rw2DensityGrid.alphaDensity(24 * j + 7);
                auto D_sig_xy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 7);
                auto D_sig_xy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 7);
                auto D_sig_xy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 7);

                auto D_sig_yx_r = rw2DensityGrid.alphaDensity(24 * j + 6);
                auto D_sig_yx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 6);
                auto D_sig_yx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 6);
                auto D_sig_yx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 6);

                auto D_sig_yx_i = rw2DensityGrid.alphaDensity(24 * j + 7);
                auto D_sig_yx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 7);
                auto D_sig_yx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 7);
                auto D_sig_yx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 7);

                auto D_sig_xz_r = rw2DensityGrid.alphaDensity(24 * j + 8);
                auto D_sig_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 8);
                auto D_sig_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 8);
                auto D_sig_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 8);

                auto D_sig_xz_i = rw2DensityGrid.alphaDensity(24 * j + 9);
                auto D_sig_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 9);
                auto D_sig_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 9);
                auto D_sig_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 9);

                auto D_sig_zx_r = rw2DensityGrid.alphaDensity(24 * j + 8);
                auto D_sig_zx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 8);
                auto D_sig_zx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 8);
                auto D_sig_zx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 8);

                auto D_sig_zx_i = rw2DensityGrid.alphaDensity(24 * j + 9);
                auto D_sig_zx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 9);
                auto D_sig_zx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 9);
                auto D_sig_zx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 9);


                auto D_sig_yz_r = rw2DensityGrid.alphaDensity(24 * j + 10);
                auto D_sig_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 10);
                auto D_sig_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 10);
                auto D_sig_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 10);

                auto D_sig_yz_i = rw2DensityGrid.alphaDensity(24 * j + 11);
                auto D_sig_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 11);
                auto D_sig_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 11);
                auto D_sig_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 11);


                auto D_sig_zy_r = rw2DensityGrid.alphaDensity(24 * j + 10);
                auto D_sig_zy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 10);
                auto D_sig_zy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 10);
                auto D_sig_zy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 10);

                auto D_sig_zy_i = rw2DensityGrid.alphaDensity(24 * j + 11);
                auto D_sig_zy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 11);
                auto D_sig_zy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 11);
                auto D_sig_zy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 11);


                auto D_lamtau_xx_r = rw2DensityGrid.alphaDensity(24 * j + 12);
                auto D_lamtau_xx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 12);
                auto D_lamtau_xx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 12);
                auto D_lamtau_xx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 12);


                auto D_lamtau_xx_i = rw2DensityGrid.alphaDensity(24 * j + 13);
                auto D_lamtau_xx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 13);
                auto D_lamtau_xx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 13);
                auto D_lamtau_xx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 13);


                auto D_lamtau_yy_r = rw2DensityGrid.alphaDensity(24 * j + 14);
                auto D_lamtau_yy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 14);
                auto D_lamtau_yy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 14);
                auto D_lamtau_yy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 14);


                auto D_lamtau_yy_i = rw2DensityGrid.alphaDensity(24 * j + 15);
                auto D_lamtau_yy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 15);
                auto D_lamtau_yy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 15);
                auto D_lamtau_yy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 15);


                auto D_lamtau_zz_r = rw2DensityGrid.alphaDensity(24 * j + 16);
                auto D_lamtau_zz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 16);
                auto D_lamtau_zz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 16);
                auto D_lamtau_zz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 16);


                auto D_lamtau_zz_i = rw2DensityGrid.alphaDensity(24 * j + 17);
                auto D_lamtau_zz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 17);
                auto D_lamtau_zz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 17);
                auto D_lamtau_zz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 17);


                auto D_lamtau_xy_r = rw2DensityGrid.alphaDensity(24 * j + 18);
                auto D_lamtau_xy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 18);
                auto D_lamtau_xy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 18);
                auto D_lamtau_xy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 18);


                auto D_lamtau_xy_i = rw2DensityGrid.alphaDensity(24 * j + 19);
                auto D_lamtau_xy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 19);
                auto D_lamtau_xy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 19);
                auto D_lamtau_xy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 19);


                auto D_lamtau_yx_r = rw2DensityGrid.alphaDensity(24 * j + 18);
                auto D_lamtau_yx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 18);
                auto D_lamtau_yx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 18);
                auto D_lamtau_yx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 18);


                auto D_lamtau_yx_i = rw2DensityGrid.alphaDensity(24 * j + 19);
                auto D_lamtau_yx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 19);
                auto D_lamtau_yx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 19);
                auto D_lamtau_yx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 19);


                auto D_lamtau_xz_r = rw2DensityGrid.alphaDensity(24 * j + 20);
                auto D_lamtau_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 20);
                auto D_lamtau_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 20);
                auto D_lamtau_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 20);

                auto D_lamtau_xz_i = rw2DensityGrid.alphaDensity(24 * j + 21);
                auto D_lamtau_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 21);
                auto D_lamtau_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 21);
                auto D_lamtau_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 21);


                auto D_lamtau_zx_r = rw2DensityGrid.alphaDensity(24 * j + 20);
                auto D_lamtau_zx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 20);
                auto D_lamtau_zx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 20);
                auto D_lamtau_zx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 20);

                auto D_lamtau_zx_i = rw2DensityGrid.alphaDensity(24 * j + 21);
                auto D_lamtau_zx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 21);
                auto D_lamtau_zx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 21);
                auto D_lamtau_zx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 21);


                auto D_lamtau_yz_r = rw2DensityGrid.alphaDensity(24 * j + 22);
                auto D_lamtau_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 22);
                auto D_lamtau_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 22);
                auto D_lamtau_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 22);

                auto D_lamtau_zy_r = rw2DensityGrid.alphaDensity(24 * j + 22);
                auto D_lamtau_zy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 22);
                auto D_lamtau_zy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 22);
                auto D_lamtau_zy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 22);

                auto D_lamtau_yz_i = rw2DensityGrid.alphaDensity(24 * j + 23);
                auto D_lamtau_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 23);
                auto D_lamtau_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 23);
                auto D_lamtau_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 23);

                auto D_lamtau_zy_i = rw2DensityGrid.alphaDensity(24 * j + 23);
                auto D_lamtau_zy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 23);
                auto D_lamtau_zy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 23);
                auto D_lamtau_zy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 23);

                for (int32_t i = 0; i < npoints; i++)
                {
                   
                    pix_r[i] = 12.0 * (rhoBx_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                    -rhoBx_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                    -rhoBx_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                    -rhoBx_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_i[i] = 12.0 * ( -rhoBx_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_r[i] += 6.0 * (rhoBx_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_r[i] += 12.0 * (rhoBx_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    pix_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    pix_r[i] += 6.0 * (rhoBy_r[i] * rhoCx_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * rhoBy_r[i]);

                    pix_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCx_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * rhoBy_i[i]);

                    pix_r[i] += 12.0 * (rhoBx_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    pix_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    pix_r[i] += 6.0 * (rhoBz_r[i] * rhoCx_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * rhoBz_r[i]);

                    pix_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCx_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * rhoBz_i[i]);

                    piy_r[i] = 12.0 * (rhoBy_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piy_i[i] = 12.0 * ( -rhoBy_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piy_r[i] += 6.0 * (rhoBx_r[i] * rhoCy_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * rhoBx_r[i]);

                    piy_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCy_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * rhoBx_i[i]);

                    piy_r[i] += 12.0 * (rhoBy_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_r[i] += 6.0 * (rhoBy_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_r[i] += 12.0 * (rhoBy_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piy_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piy_r[i] += 6.0 * (rhoBz_r[i] * rhoCy_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * rhoBz_r[i]);

                    piy_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCy_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * rhoBz_i[i]);

                    piz_r[i] = 12.0 * (rhoBz_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piz_i[i] = 12.0 * ( -rhoBz_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piz_r[i] += 6.0 * (rhoBx_r[i] * rhoCz_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * rhoBx_r[i]);

                    piz_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCz_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * rhoBx_i[i]);

                    piz_r[i] += 12.0 * (rhoBz_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piz_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piz_r[i] += 6.0 * (rhoBy_r[i] * rhoCz_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * rhoBy_r[i]);

                    piz_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCz_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * rhoBy_i[i]);

                    piz_r[i] += 12.0 * (rhoBz_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_r[i] += 6.0 * (rhoBz_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    pix_X_i[i] = 12.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_X_r[i] += 6.0 * (gradBx_x_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_X_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBx_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_X_r[i] += 12.0 * (gradBx_x_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    pix_X_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    pix_X_r[i] += 6.0 * (gradBy_x_r[i] * rhoCx_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * rhoBy_r[i]);

                    pix_X_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * rhoBy_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBy_r[i] * gradCx_x_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * rhoBy_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * rhoBy_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBy_x_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBy_x_i[i]);

                    pix_X_r[i] += 12.0 * (gradBx_x_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    pix_X_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    pix_X_r[i] += 12.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    pix_X_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    pix_X_r[i] += 6.0 * (gradBz_x_r[i] * rhoCx_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * rhoBz_r[i]);

                    pix_X_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * rhoBz_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBz_r[i] * gradCx_x_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * rhoBz_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * rhoBz_i[i]);

                    pix_X_r[i] += 6.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBz_x_r[i]);

                    pix_X_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBz_x_i[i]);

                    pix_Y_r[i] = 12.0 * (gradBx_y_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_Y_i[i] = 12.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_Y_r[i] += 6.0 * (gradBx_y_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_Y_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBx_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_Y_r[i] += 12.0 * (gradBx_y_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    pix_Y_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    pix_Y_r[i] += 6.0 * (gradBy_y_r[i] * rhoCx_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * rhoBy_r[i]);

                    pix_Y_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * rhoBy_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBy_r[i] * gradCx_y_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * rhoBy_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * rhoBy_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBy_y_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBy_y_i[i]);

                    pix_Y_r[i] += 12.0 * (gradBx_y_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    pix_Y_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    pix_Y_r[i] += 12.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    pix_Y_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    pix_Y_r[i] += 6.0 * (gradBz_y_r[i] * rhoCx_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * rhoBz_r[i]);

                    pix_Y_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * rhoBz_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBz_r[i] * gradCx_y_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * rhoBz_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * rhoBz_i[i]);

                    pix_Y_r[i] += 6.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBz_y_r[i]);

                    pix_Y_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBz_y_i[i]);

                    pix_Z_r[i] = 12.0 * (gradBx_z_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_Z_i[i] = 12.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_Z_r[i] += 6.0 * (gradBx_z_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    pix_Z_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBx_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBx_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_Z_r[i] += 12.0 * (gradBx_z_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    pix_Z_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    pix_Z_r[i] += 6.0 * (gradBy_z_r[i] * rhoCx_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * rhoBy_r[i]);

                    pix_Z_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * rhoBy_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBy_r[i] * gradCx_z_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * rhoBy_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * rhoBy_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBy_z_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBy_z_i[i]);

                    pix_Z_r[i] += 12.0 * (gradBx_z_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    pix_Z_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    pix_Z_r[i] += 12.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    pix_Z_i[i] += 12.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_Z_r[i] += 6.0 * (gradBz_z_r[i] * rhoCx_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * rhoBz_r[i]);

                    pix_Z_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * rhoBz_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBz_r[i] * gradCx_z_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * rhoBz_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * rhoBz_i[i]);

                    pix_Z_r[i] += 6.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBz_z_r[i]);

                    pix_Z_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBz_z_i[i]);

                    piy_X_r[i] = 12.0 * (gradBy_x_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piy_X_i[i] = 12.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piy_X_r[i] += 6.0 * (gradBx_x_r[i] * rhoCy_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * rhoBx_r[i]);

                    piy_X_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * rhoBx_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBx_r[i] * gradCy_x_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * rhoBx_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * rhoBx_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBx_x_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBx_x_i[i]);

                    piy_X_r[i] += 12.0 * (gradBy_x_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_X_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_X_r[i] += 6.0 * (gradBy_x_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_X_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBy_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_X_r[i] += 12.0 * (gradBy_x_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piy_X_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piy_X_r[i] += 12.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piy_X_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piy_X_r[i] += 6.0 * (gradBz_x_r[i] * rhoCy_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * rhoBz_r[i]);

                    piy_X_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * rhoBz_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBz_r[i] * gradCy_x_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * rhoBz_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * rhoBz_i[i]);

                    piy_X_r[i] += 6.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBz_x_r[i]);

                    piy_X_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBz_x_i[i]);

                    piy_Y_r[i] = 12.0 * (gradBy_y_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piy_Y_i[i] = 12.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piy_Y_r[i] += 6.0 * (gradBx_y_r[i] * rhoCy_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * rhoBx_r[i]);

                    piy_Y_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * rhoBx_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBx_r[i] * gradCy_y_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * rhoBx_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * rhoBx_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBx_y_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBx_y_i[i]);

                    piy_Y_r[i] += 12.0 * (gradBy_y_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_Y_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_Y_r[i] += 6.0 * (gradBy_y_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_Y_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBy_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_Y_r[i] += 12.0 * (gradBy_y_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piy_Y_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piy_Y_r[i] += 12.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piy_Y_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piy_Y_r[i] += 6.0 * (gradBz_y_r[i] * rhoCy_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * rhoBz_r[i]);

                    piy_Y_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * rhoBz_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBz_r[i] * gradCy_y_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * rhoBz_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * rhoBz_i[i]);

                    piy_Y_r[i] += 6.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBz_y_r[i]);

                    piy_Y_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBz_y_i[i]);

                    piy_Z_r[i] = 12.0 * (gradBy_z_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piy_Z_i[i] = 12.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -rhoBy_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -rhoBy_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -rhoBy_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +rhoBy_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piy_Z_r[i] += 6.0 * (gradBx_z_r[i] * rhoCy_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * rhoBx_r[i]);

                    piy_Z_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * rhoBx_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBx_r[i] * gradCy_z_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * rhoBx_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * rhoBx_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBx_r[i] * rhoCy_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * rhoCy_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * rhoCy_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * rhoCy_i[i] * gradBx_z_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCy_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * rhoCy_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCy_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCy_r[i] * gradBx_z_i[i]);

                    piy_Z_r[i] += 12.0 * (gradBy_z_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_Z_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_Z_r[i] += 6.0 * (gradBy_z_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piy_Z_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBy_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBy_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_Z_r[i] += 12.0 * (gradBy_z_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piy_Z_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piy_Z_r[i] += 12.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piy_Z_i[i] += 12.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piy_Z_r[i] += 6.0 * (gradBz_z_r[i] * rhoCy_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * rhoBz_r[i]);

                    piy_Z_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * rhoBz_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBz_r[i] * gradCy_z_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * rhoBz_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * rhoBz_i[i]);

                    piy_Z_r[i] += 6.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBz_z_r[i]);

                    piy_Z_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBz_z_i[i]);

                    piz_X_r[i] = 12.0 * (gradBz_x_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piz_X_i[i] = 12.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piz_X_r[i] += 6.0 * (gradBx_x_r[i] * rhoCz_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * rhoBx_r[i]);

                    piz_X_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * rhoBx_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBx_r[i] * gradCz_x_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * rhoBx_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * rhoBx_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBx_x_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBx_x_i[i]);

                    piz_X_r[i] += 12.0 * (gradBz_x_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piz_X_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piz_X_r[i] += 6.0 * (gradBy_x_r[i] * rhoCz_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * rhoBy_r[i]);

                    piz_X_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * rhoBy_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBy_r[i] * gradCz_x_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * rhoBy_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * rhoBy_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBy_x_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBy_x_i[i]);

                    piz_X_r[i] += 12.0 * (gradBz_x_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_X_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_X_r[i] += 12.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_X_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_X_r[i] += 6.0 * (gradBz_x_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_X_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBz_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_X_r[i] += 6.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_X_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_Y_r[i] = 12.0 * (gradBz_y_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piz_Y_i[i] = 12.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piz_Y_r[i] += 6.0 * (gradBx_y_r[i] * rhoCz_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * rhoBx_r[i]);

                    piz_Y_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * rhoBx_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBx_r[i] * gradCz_y_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * rhoBx_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * rhoBx_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBx_y_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBx_y_i[i]);

                    piz_Y_r[i] += 12.0 * (gradBz_y_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piz_Y_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piz_Y_r[i] += 6.0 * (gradBy_y_r[i] * rhoCz_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * rhoBy_r[i]);

                    piz_Y_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * rhoBy_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBy_r[i] * gradCz_y_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * rhoBy_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * rhoBy_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBy_y_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBy_y_i[i]);

                    piz_Y_r[i] += 12.0 * (gradBz_y_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_Y_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_Y_r[i] += 12.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_Y_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_Y_r[i] += 6.0 * (gradBz_y_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_Y_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBz_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_Y_r[i] += 6.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_Y_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_Z_r[i] = 12.0 * (gradBz_z_r[i] * rhoCx_r[i] * rhoBx_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * rhoBx_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * rhoBx_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * rhoBx_r[i]);

                    piz_Z_i[i] = 12.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * rhoBx_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * rhoBx_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -rhoBz_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -rhoBz_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -rhoBz_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +rhoBz_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piz_Z_r[i] += 6.0 * (gradBx_z_r[i] * rhoCz_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * rhoBx_r[i]);

                    piz_Z_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * rhoBx_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBx_r[i] * gradCz_z_r[i] * rhoBx_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * rhoBx_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * rhoBx_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * rhoBx_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * rhoBx_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * rhoBx_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * rhoBx_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBx_r[i] * rhoCz_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * rhoCz_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * rhoCz_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * rhoCz_i[i] * gradBx_z_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBx_i[i] * rhoCz_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * rhoCz_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCz_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * rhoCz_r[i] * gradBx_z_i[i]);

                    piz_Z_r[i] += 12.0 * (gradBz_z_r[i] * rhoCy_r[i] * rhoBy_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * rhoBy_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * rhoBy_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * rhoBy_r[i]);

                    piz_Z_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * rhoBy_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * rhoBy_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -rhoBz_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -rhoBz_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -rhoBz_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +rhoBz_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piz_Z_r[i] += 6.0 * (gradBy_z_r[i] * rhoCz_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * rhoBy_r[i]);

                    piz_Z_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * rhoBy_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBy_r[i] * gradCz_z_r[i] * rhoBy_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * rhoBy_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * rhoBy_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * rhoBy_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * rhoBy_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * rhoBy_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * rhoBy_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBy_r[i] * rhoCz_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * rhoCz_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * rhoCz_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * rhoCz_i[i] * gradBy_z_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBy_i[i] * rhoCz_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * rhoCz_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCz_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * rhoCz_r[i] * gradBy_z_i[i]);

                    piz_Z_r[i] += 12.0 * (gradBz_z_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_Z_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_Z_r[i] += 12.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_Z_i[i] += 12.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_Z_r[i] += 6.0 * (gradBz_z_r[i] * rhoCz_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * rhoBz_r[i]);

                    piz_Z_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * rhoBz_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBz_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_Z_r[i] += 6.0 * (rhoBz_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_Z_i[i] += 6.0 * ( -rhoBz_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_XX_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_XX_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_XX_r[i] += 12.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_XX_r[i] += 12.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_XX_r[i] += 6.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_XX_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    pix_XX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    pix_XX_r[i] += 12.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    pix_XX_r[i] += 12.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_x_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * rhoBy_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * rhoBy_i[i]);

                    pix_XX_r[i] += 6.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBy_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBy_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBy_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBy_x_i[i]);

                    pix_XX_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    pix_XX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    pix_XX_r[i] += 12.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    pix_XX_r[i] += 12.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    pix_XX_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_x_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * rhoBz_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * rhoBz_i[i]);

                    pix_XX_r[i] += 6.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBz_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBz_x_i[i]);

                    pix_XX_r[i] += 6.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBz_x_r[i]);

                    pix_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBz_x_i[i]);

                    pix_XY_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_XY_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_XY_r[i] += 12.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_XY_r[i] += 12.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_XY_r[i] += 6.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_XY_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    pix_XY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    pix_XY_r[i] += 12.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    pix_XY_r[i] += 12.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_y_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * rhoBy_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * rhoBy_i[i]);

                    pix_XY_r[i] += 6.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBy_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBy_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBy_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBy_y_i[i]);

                    pix_XY_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    pix_XY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    pix_XY_r[i] += 12.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    pix_XY_r[i] += 12.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    pix_XY_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_y_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * rhoBz_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * rhoBz_i[i]);

                    pix_XY_r[i] += 6.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBz_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBz_y_i[i]);

                    pix_XY_r[i] += 6.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBz_y_r[i]);

                    pix_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBz_y_i[i]);

                    pix_XZ_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_XZ_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_XZ_r[i] += 12.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_XZ_r[i] += 12.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_XZ_r[i] += 6.0 * (rhoBx_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBx_x_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_XZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    pix_XZ_r[i] += 12.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    pix_XZ_r[i] += 12.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_z_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * rhoBy_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * rhoBy_i[i]);

                    pix_XZ_r[i] += 6.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBy_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBy_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBy_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBy_z_i[i]);

                    pix_XZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    pix_XZ_r[i] += 12.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    pix_XZ_r[i] += 12.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    pix_XZ_i[i] += 12.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_z_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * rhoBz_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * rhoBz_i[i]);

                    pix_XZ_r[i] += 6.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBz_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBz_z_i[i]);

                    pix_XZ_r[i] += 6.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBz_z_r[i]);

                    pix_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBz_z_i[i]);

                    pix_YX_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_YX_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_YX_r[i] += 12.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_YX_r[i] += 12.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_YX_r[i] += 6.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_YX_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    pix_YX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    pix_YX_r[i] += 12.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    pix_YX_r[i] += 12.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_x_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * rhoBy_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * rhoBy_i[i]);

                    pix_YX_r[i] += 6.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBy_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBy_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBy_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBy_x_i[i]);

                    pix_YX_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    pix_YX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    pix_YX_r[i] += 12.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    pix_YX_r[i] += 12.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    pix_YX_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_x_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * rhoBz_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * rhoBz_i[i]);

                    pix_YX_r[i] += 6.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBz_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBz_x_i[i]);

                    pix_YX_r[i] += 6.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBz_x_r[i]);

                    pix_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBz_x_i[i]);

                    pix_YY_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_YY_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_YY_r[i] += 12.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_YY_r[i] += 12.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_YY_r[i] += 6.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_YY_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    pix_YY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    pix_YY_r[i] += 12.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    pix_YY_r[i] += 12.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_y_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * rhoBy_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * rhoBy_i[i]);

                    pix_YY_r[i] += 6.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBy_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBy_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBy_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBy_y_i[i]);

                    pix_YY_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    pix_YY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    pix_YY_r[i] += 12.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    pix_YY_r[i] += 12.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    pix_YY_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_y_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * rhoBz_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * rhoBz_i[i]);

                    pix_YY_r[i] += 6.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBz_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBz_y_i[i]);

                    pix_YY_r[i] += 6.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBz_y_r[i]);

                    pix_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBz_y_i[i]);

                    pix_YZ_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_YZ_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_YZ_r[i] += 12.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_YZ_r[i] += 12.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_YZ_r[i] += 6.0 * (rhoBx_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBx_y_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_YZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    pix_YZ_r[i] += 12.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    pix_YZ_r[i] += 12.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_z_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * rhoBy_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * rhoBy_i[i]);

                    pix_YZ_r[i] += 6.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBy_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBy_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBy_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBy_z_i[i]);

                    pix_YZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    pix_YZ_r[i] += 12.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    pix_YZ_r[i] += 12.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    pix_YZ_i[i] += 12.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_z_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * rhoBz_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * rhoBz_i[i]);

                    pix_YZ_r[i] += 6.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBz_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBz_z_i[i]);

                    pix_YZ_r[i] += 6.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBz_z_r[i]);

                    pix_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBz_z_i[i]);

                    pix_ZX_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_ZX_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_ZX_r[i] += 12.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_ZX_r[i] += 12.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    pix_ZX_r[i] += 6.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    pix_ZX_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    pix_ZX_r[i] += 12.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    pix_ZX_r[i] += 12.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_x_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * rhoBy_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * rhoBy_i[i]);

                    pix_ZX_r[i] += 6.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBy_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBy_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBy_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBy_x_i[i]);

                    pix_ZX_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    pix_ZX_r[i] += 12.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    pix_ZX_r[i] += 12.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    pix_ZX_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_x_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * rhoBz_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * rhoBz_i[i]);

                    pix_ZX_r[i] += 6.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBz_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBz_x_i[i]);

                    pix_ZX_r[i] += 6.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBz_x_r[i]);

                    pix_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBz_x_i[i]);

                    pix_ZY_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_ZY_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_ZY_r[i] += 12.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_ZY_r[i] += 12.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    pix_ZY_r[i] += 6.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    pix_ZY_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    pix_ZY_r[i] += 12.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    pix_ZY_r[i] += 12.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_y_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * rhoBy_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * rhoBy_i[i]);

                    pix_ZY_r[i] += 6.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBy_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBy_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBy_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBy_y_i[i]);

                    pix_ZY_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    pix_ZY_r[i] += 12.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    pix_ZY_r[i] += 12.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    pix_ZY_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_y_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * rhoBz_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * rhoBz_i[i]);

                    pix_ZY_r[i] += 6.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBz_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBz_y_i[i]);

                    pix_ZY_r[i] += 6.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBz_y_r[i]);

                    pix_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBz_y_i[i]);

                    pix_ZZ_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_ZZ_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_ZZ_r[i] += 12.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_ZZ_r[i] += 12.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    pix_ZZ_r[i] += 6.0 * (rhoBx_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    pix_ZZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    pix_ZZ_r[i] += 12.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    pix_ZZ_r[i] += 12.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_z_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * rhoBy_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * rhoBy_i[i]);

                    pix_ZZ_r[i] += 6.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBy_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBy_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBy_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBy_z_i[i]);

                    pix_ZZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    pix_ZZ_r[i] += 12.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    pix_ZZ_r[i] += 12.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    pix_ZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_z_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * rhoBz_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * rhoBz_i[i]);

                    pix_ZZ_r[i] += 6.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBz_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBz_z_i[i]);

                    pix_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBz_z_r[i]);

                    pix_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBz_z_i[i]);

                    piy_XX_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piy_XX_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piy_XX_r[i] += 12.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piy_XX_r[i] += 12.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_x_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * rhoBx_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * rhoBx_i[i]);

                    piy_XX_r[i] += 6.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBx_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBx_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBx_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBx_x_i[i]);

                    piy_XX_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_XX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_XX_r[i] += 12.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_XX_r[i] += 12.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_XX_r[i] += 6.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_XX_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piy_XX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piy_XX_r[i] += 12.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piy_XX_r[i] += 12.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piy_XX_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_x_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * rhoBz_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * rhoBz_i[i]);

                    piy_XX_r[i] += 6.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBz_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBz_x_i[i]);

                    piy_XX_r[i] += 6.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBz_x_r[i]);

                    piy_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBz_x_i[i]);

                    piy_XY_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piy_XY_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piy_XY_r[i] += 12.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piy_XY_r[i] += 12.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_y_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * rhoBx_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * rhoBx_i[i]);

                    piy_XY_r[i] += 6.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBx_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBx_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBx_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBx_y_i[i]);

                    piy_XY_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_XY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_XY_r[i] += 12.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_XY_r[i] += 12.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_XY_r[i] += 6.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_XY_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piy_XY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piy_XY_r[i] += 12.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piy_XY_r[i] += 12.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piy_XY_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_y_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * rhoBz_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * rhoBz_i[i]);

                    piy_XY_r[i] += 6.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBz_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBz_y_i[i]);

                    piy_XY_r[i] += 6.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBz_y_r[i]);

                    piy_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBz_y_i[i]);

                    piy_XZ_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piy_XZ_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piy_XZ_r[i] += 12.0 * (rhoBy_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -rhoBy_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -rhoBy_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -rhoBy_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +rhoBy_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piy_XZ_r[i] += 12.0 * (gradBy_x_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBy_x_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBy_x_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBy_x_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBy_x_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_z_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * rhoBx_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * rhoBx_i[i]);

                    piy_XZ_r[i] += 6.0 * (rhoBx_r[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCy_x_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCy_x_i[i] * gradBx_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_x_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_x_r[i] * gradBx_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBx_x_r[i] * rhoCy_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * rhoCy_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * rhoCy_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * rhoCy_i[i] * gradBx_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCy_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * rhoCy_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCy_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCy_r[i] * gradBx_z_i[i]);

                    piy_XZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_XZ_r[i] += 12.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_XZ_r[i] += 12.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_XZ_r[i] += 6.0 * (rhoBy_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBy_x_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_XZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piy_XZ_r[i] += 12.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piy_XZ_r[i] += 12.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piy_XZ_i[i] += 12.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_z_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * rhoBz_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * rhoBz_i[i]);

                    piy_XZ_r[i] += 6.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBz_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBz_z_i[i]);

                    piy_XZ_r[i] += 6.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBz_z_r[i]);

                    piy_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBz_z_i[i]);

                    piy_YX_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piy_YX_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piy_YX_r[i] += 12.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piy_YX_r[i] += 12.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_x_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * rhoBx_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * rhoBx_i[i]);

                    piy_YX_r[i] += 6.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBx_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBx_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBx_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBx_x_i[i]);

                    piy_YX_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_YX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_YX_r[i] += 12.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_YX_r[i] += 12.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_YX_r[i] += 6.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_YX_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piy_YX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piy_YX_r[i] += 12.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piy_YX_r[i] += 12.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piy_YX_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_x_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * rhoBz_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * rhoBz_i[i]);

                    piy_YX_r[i] += 6.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBz_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBz_x_i[i]);

                    piy_YX_r[i] += 6.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBz_x_r[i]);

                    piy_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBz_x_i[i]);

                    piy_YY_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piy_YY_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piy_YY_r[i] += 12.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piy_YY_r[i] += 12.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_y_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * rhoBx_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * rhoBx_i[i]);

                    piy_YY_r[i] += 6.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBx_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBx_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBx_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBx_y_i[i]);

                    piy_YY_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_YY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_YY_r[i] += 12.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_YY_r[i] += 12.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_YY_r[i] += 6.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_YY_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piy_YY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piy_YY_r[i] += 12.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piy_YY_r[i] += 12.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piy_YY_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_y_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * rhoBz_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * rhoBz_i[i]);

                    piy_YY_r[i] += 6.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBz_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBz_y_i[i]);

                    piy_YY_r[i] += 6.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBz_y_r[i]);

                    piy_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBz_y_i[i]);

                    piy_YZ_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piy_YZ_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piy_YZ_r[i] += 12.0 * (rhoBy_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -rhoBy_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -rhoBy_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -rhoBy_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +rhoBy_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piy_YZ_r[i] += 12.0 * (gradBy_y_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBy_y_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBy_y_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBy_y_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBy_y_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_z_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * rhoBx_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * rhoBx_i[i]);

                    piy_YZ_r[i] += 6.0 * (rhoBx_r[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCy_y_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCy_y_i[i] * gradBx_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_y_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_y_r[i] * gradBx_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBx_y_r[i] * rhoCy_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * rhoCy_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * rhoCy_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * rhoCy_i[i] * gradBx_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCy_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * rhoCy_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCy_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCy_r[i] * gradBx_z_i[i]);

                    piy_YZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_YZ_r[i] += 12.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_YZ_r[i] += 12.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_YZ_r[i] += 6.0 * (rhoBy_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBy_y_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_YZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piy_YZ_r[i] += 12.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piy_YZ_r[i] += 12.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piy_YZ_i[i] += 12.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_z_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * rhoBz_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * rhoBz_i[i]);

                    piy_YZ_r[i] += 6.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBz_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBz_z_i[i]);

                    piy_YZ_r[i] += 6.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBz_z_r[i]);

                    piy_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBz_z_i[i]);

                    piy_ZX_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piy_ZX_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piy_ZX_r[i] += 12.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piy_ZX_r[i] += 12.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_x_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * rhoBx_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * rhoBx_i[i]);

                    piy_ZX_r[i] += 6.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBx_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBx_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBx_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBx_x_i[i]);

                    piy_ZX_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_ZX_r[i] += 12.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_ZX_r[i] += 12.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piy_ZX_r[i] += 6.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piy_ZX_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piy_ZX_r[i] += 12.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piy_ZX_r[i] += 12.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piy_ZX_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_x_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * rhoBz_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * rhoBz_i[i]);

                    piy_ZX_r[i] += 6.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBz_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBz_x_i[i]);

                    piy_ZX_r[i] += 6.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBz_x_r[i]);

                    piy_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBz_x_i[i]);

                    piy_ZY_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piy_ZY_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piy_ZY_r[i] += 12.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piy_ZY_r[i] += 12.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_y_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * rhoBx_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * rhoBx_i[i]);

                    piy_ZY_r[i] += 6.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBx_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBx_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBx_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBx_y_i[i]);

                    piy_ZY_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_ZY_r[i] += 12.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_ZY_r[i] += 12.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piy_ZY_r[i] += 6.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piy_ZY_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piy_ZY_r[i] += 12.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piy_ZY_r[i] += 12.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piy_ZY_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_y_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * rhoBz_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * rhoBz_i[i]);

                    piy_ZY_r[i] += 6.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBz_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBz_y_i[i]);

                    piy_ZY_r[i] += 6.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBz_y_r[i]);

                    piy_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBz_y_i[i]);

                    piy_ZZ_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piy_ZZ_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piy_ZZ_r[i] += 12.0 * (rhoBy_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -rhoBy_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -rhoBy_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -rhoBy_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +rhoBy_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +rhoBy_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piy_ZZ_r[i] += 12.0 * (gradBy_z_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBy_z_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBy_z_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBy_z_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBy_z_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_z_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * rhoBx_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * rhoBx_i[i]);

                    piy_ZZ_r[i] += 6.0 * (rhoBx_r[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCy_z_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCy_z_i[i] * gradBx_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_z_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCy_z_r[i] * gradBx_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * rhoCy_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * rhoCy_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * rhoCy_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * rhoCy_i[i] * gradBx_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCy_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * rhoCy_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCy_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCy_r[i] * gradBx_z_i[i]);

                    piy_ZZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_ZZ_r[i] += 12.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_ZZ_r[i] += 12.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piy_ZZ_r[i] += 6.0 * (rhoBy_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piy_ZZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piy_ZZ_r[i] += 12.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piy_ZZ_r[i] += 12.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piy_ZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_z_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * rhoBz_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * rhoBz_i[i]);

                    piy_ZZ_r[i] += 6.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBz_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBz_z_i[i]);

                    piy_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBz_z_r[i]);

                    piy_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBz_z_i[i]);

                    piz_XX_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piz_XX_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piz_XX_r[i] += 12.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piz_XX_r[i] += 12.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_x_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * rhoBx_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * rhoBx_i[i]);

                    piz_XX_r[i] += 6.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBx_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBx_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBx_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBx_x_i[i]);

                    piz_XX_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piz_XX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piz_XX_r[i] += 12.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piz_XX_r[i] += 12.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_x_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * rhoBy_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * rhoBy_i[i]);

                    piz_XX_r[i] += 6.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBy_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBy_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBy_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBy_x_i[i]);

                    piz_XX_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_XX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_XX_r[i] += 12.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_XX_r[i] += 12.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_XX_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_XX_r[i] += 6.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_XX_r[i] += 6.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_XX_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_XY_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piz_XY_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piz_XY_r[i] += 12.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piz_XY_r[i] += 12.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_y_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * rhoBx_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * rhoBx_i[i]);

                    piz_XY_r[i] += 6.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBx_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBx_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBx_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBx_y_i[i]);

                    piz_XY_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piz_XY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piz_XY_r[i] += 12.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piz_XY_r[i] += 12.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_y_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * rhoBy_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * rhoBy_i[i]);

                    piz_XY_r[i] += 6.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBy_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBy_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBy_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBy_y_i[i]);

                    piz_XY_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_XY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_XY_r[i] += 12.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_XY_r[i] += 12.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_XY_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_XY_r[i] += 6.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_XY_r[i] += 6.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_XY_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_XZ_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piz_XZ_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piz_XZ_r[i] += 12.0 * (rhoBz_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -rhoBz_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -rhoBz_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -rhoBz_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +rhoBz_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piz_XZ_r[i] += 12.0 * (gradBz_x_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBz_x_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBz_x_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBz_x_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBz_x_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_z_r[i] * rhoBx_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * rhoBx_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * rhoBx_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * rhoBx_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * rhoBx_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * rhoBx_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * rhoBx_i[i]);

                    piz_XZ_r[i] += 6.0 * (rhoBx_r[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCz_x_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCz_x_i[i] * gradBx_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_x_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_x_r[i] * gradBx_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBx_x_r[i] * rhoCz_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * rhoCz_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * rhoCz_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * rhoCz_i[i] * gradBx_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBx_x_i[i] * rhoCz_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * rhoCz_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCz_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * rhoCz_r[i] * gradBx_z_i[i]);

                    piz_XZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piz_XZ_r[i] += 12.0 * (rhoBz_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -rhoBz_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -rhoBz_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -rhoBz_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +rhoBz_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piz_XZ_r[i] += 12.0 * (gradBz_x_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBz_x_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBz_x_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBz_x_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBz_x_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_z_r[i] * rhoBy_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * rhoBy_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * rhoBy_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * rhoBy_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * rhoBy_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * rhoBy_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * rhoBy_i[i]);

                    piz_XZ_r[i] += 6.0 * (rhoBy_r[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCz_x_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCz_x_i[i] * gradBy_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_x_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_x_r[i] * gradBy_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBy_x_r[i] * rhoCz_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * rhoCz_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * rhoCz_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * rhoCz_i[i] * gradBy_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBy_x_i[i] * rhoCz_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * rhoCz_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCz_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * rhoCz_r[i] * gradBy_z_i[i]);

                    piz_XZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_XZ_r[i] += 12.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_XZ_r[i] += 12.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_XZ_i[i] += 12.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_XZ_r[i] += 6.0 * (rhoBz_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_XZ_r[i] += 6.0 * (gradBz_x_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_XZ_i[i] += 6.0 * ( -gradBz_x_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_YX_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piz_YX_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piz_YX_r[i] += 12.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piz_YX_r[i] += 12.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_x_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * rhoBx_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * rhoBx_i[i]);

                    piz_YX_r[i] += 6.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBx_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBx_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBx_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBx_x_i[i]);

                    piz_YX_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piz_YX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piz_YX_r[i] += 12.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piz_YX_r[i] += 12.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_x_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * rhoBy_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * rhoBy_i[i]);

                    piz_YX_r[i] += 6.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBy_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBy_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBy_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBy_x_i[i]);

                    piz_YX_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_YX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_YX_r[i] += 12.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_YX_r[i] += 12.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_YX_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_YX_r[i] += 6.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_YX_r[i] += 6.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_YX_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_YY_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piz_YY_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piz_YY_r[i] += 12.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piz_YY_r[i] += 12.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_y_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * rhoBx_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * rhoBx_i[i]);

                    piz_YY_r[i] += 6.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBx_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBx_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBx_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBx_y_i[i]);

                    piz_YY_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piz_YY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piz_YY_r[i] += 12.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piz_YY_r[i] += 12.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_y_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * rhoBy_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * rhoBy_i[i]);

                    piz_YY_r[i] += 6.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBy_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBy_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBy_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBy_y_i[i]);

                    piz_YY_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_YY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_YY_r[i] += 12.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_YY_r[i] += 12.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_YY_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_YY_r[i] += 6.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_YY_r[i] += 6.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_YY_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_YZ_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piz_YZ_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piz_YZ_r[i] += 12.0 * (rhoBz_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -rhoBz_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -rhoBz_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -rhoBz_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +rhoBz_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piz_YZ_r[i] += 12.0 * (gradBz_y_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBz_y_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBz_y_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBz_y_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBz_y_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_z_r[i] * rhoBx_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * rhoBx_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * rhoBx_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * rhoBx_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * rhoBx_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * rhoBx_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * rhoBx_i[i]);

                    piz_YZ_r[i] += 6.0 * (rhoBx_r[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCz_y_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCz_y_i[i] * gradBx_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_y_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_y_r[i] * gradBx_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBx_y_r[i] * rhoCz_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * rhoCz_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * rhoCz_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * rhoCz_i[i] * gradBx_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBx_y_i[i] * rhoCz_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * rhoCz_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCz_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * rhoCz_r[i] * gradBx_z_i[i]);

                    piz_YZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piz_YZ_r[i] += 12.0 * (rhoBz_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -rhoBz_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -rhoBz_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -rhoBz_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +rhoBz_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piz_YZ_r[i] += 12.0 * (gradBz_y_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBz_y_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBz_y_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBz_y_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBz_y_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_z_r[i] * rhoBy_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * rhoBy_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * rhoBy_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * rhoBy_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * rhoBy_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * rhoBy_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * rhoBy_i[i]);

                    piz_YZ_r[i] += 6.0 * (rhoBy_r[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCz_y_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCz_y_i[i] * gradBy_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_y_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_y_r[i] * gradBy_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBy_y_r[i] * rhoCz_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * rhoCz_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * rhoCz_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * rhoCz_i[i] * gradBy_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBy_y_i[i] * rhoCz_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * rhoCz_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCz_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * rhoCz_r[i] * gradBy_z_i[i]);

                    piz_YZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_YZ_r[i] += 12.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_YZ_r[i] += 12.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_YZ_i[i] += 12.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_YZ_r[i] += 6.0 * (rhoBz_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_YZ_r[i] += 6.0 * (gradBz_y_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_YZ_i[i] += 6.0 * ( -gradBz_y_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_ZX_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_x_r[i] * rhoBx_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * rhoBx_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * rhoBx_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * rhoBx_r[i]);

                    piz_ZX_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * rhoBx_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * rhoBx_i[i]);

                    piz_ZX_r[i] += 12.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piz_ZX_r[i] += 12.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBx_x_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBx_x_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBx_x_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBx_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBx_x_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBx_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_x_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * rhoBx_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * rhoBx_i[i]);

                    piz_ZX_r[i] += 6.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBx_x_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBx_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBx_x_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBx_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBx_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBx_x_i[i]);

                    piz_ZX_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_x_r[i] * rhoBy_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * rhoBy_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * rhoBy_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * rhoBy_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * rhoBy_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * rhoBy_i[i]);

                    piz_ZX_r[i] += 12.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piz_ZX_r[i] += 12.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBy_x_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBy_x_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBy_x_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBy_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBy_x_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBy_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_x_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * rhoBy_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * rhoBy_i[i]);

                    piz_ZX_r[i] += 6.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBy_x_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBy_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBy_x_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBy_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBy_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBy_x_i[i]);

                    piz_ZX_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_ZX_r[i] += 12.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_ZX_r[i] += 12.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_ZX_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_x_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * rhoBz_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * rhoBz_i[i]);

                    piz_ZX_r[i] += 6.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_ZX_r[i] += 6.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_x_r[i]);

                    piz_ZX_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_x_i[i]);

                    piz_ZY_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_y_r[i] * rhoBx_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * rhoBx_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * rhoBx_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * rhoBx_r[i]);

                    piz_ZY_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * rhoBx_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * rhoBx_i[i]);

                    piz_ZY_r[i] += 12.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piz_ZY_r[i] += 12.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBx_y_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBx_y_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBx_y_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBx_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBx_y_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBx_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_y_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * rhoBx_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * rhoBx_i[i]);

                    piz_ZY_r[i] += 6.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBx_y_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBx_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBx_y_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBx_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBx_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBx_y_i[i]);

                    piz_ZY_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_y_r[i] * rhoBy_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * rhoBy_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * rhoBy_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * rhoBy_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * rhoBy_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * rhoBy_i[i]);

                    piz_ZY_r[i] += 12.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piz_ZY_r[i] += 12.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBy_y_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBy_y_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBy_y_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBy_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBy_y_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBy_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_y_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * rhoBy_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * rhoBy_i[i]);

                    piz_ZY_r[i] += 6.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBy_y_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBy_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBy_y_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBy_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBy_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBy_y_i[i]);

                    piz_ZY_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_ZY_r[i] += 12.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_ZY_r[i] += 12.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_ZY_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_y_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * rhoBz_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * rhoBz_i[i]);

                    piz_ZY_r[i] += 6.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_ZY_r[i] += 6.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_y_r[i]);

                    piz_ZY_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_y_i[i]);

                    piz_ZZ_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_z_r[i] * rhoBx_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * rhoBx_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * rhoBx_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * rhoBx_r[i]);

                    piz_ZZ_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * rhoBx_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * rhoBx_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * rhoBx_i[i]);

                    piz_ZZ_r[i] += 12.0 * (rhoBz_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -rhoBz_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -rhoBz_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -rhoBz_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +rhoBz_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +rhoBz_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piz_ZZ_r[i] += 12.0 * (gradBz_z_r[i] * rhoCx_r[i] * gradBx_z_r[i]
                                -gradBz_z_i[i] * rhoCx_r[i] * gradBx_z_i[i]
                                -gradBz_z_r[i] * rhoCx_i[i] * gradBx_z_i[i]
                                -gradBz_z_i[i] * rhoCx_i[i] * gradBx_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCx_i[i] * gradBx_z_i[i]
                                +gradBz_z_i[i] * rhoCx_r[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * rhoCx_i[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * rhoCx_r[i] * gradBx_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_z_r[i] * rhoBx_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * rhoBx_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * rhoBx_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * rhoBx_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * rhoBx_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * rhoBx_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * rhoBx_i[i]);

                    piz_ZZ_r[i] += 6.0 * (rhoBx_r[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                -rhoBx_i[i] * gradCz_z_r[i] * gradBx_z_i[i]
                                -rhoBx_r[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                -rhoBx_i[i] * gradCz_z_i[i] * gradBx_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -rhoBx_i[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                +rhoBx_i[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_z_i[i] * gradBx_z_r[i]
                                +rhoBx_r[i] * gradCz_z_r[i] * gradBx_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBx_z_r[i] * rhoCz_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * rhoCz_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * rhoCz_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * rhoCz_i[i] * gradBx_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * rhoCz_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * rhoCz_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCz_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * rhoCz_r[i] * gradBx_z_i[i]);

                    piz_ZZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_z_r[i] * rhoBy_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * rhoBy_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * rhoBy_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * rhoBy_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * rhoBy_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * rhoBy_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * rhoBy_i[i]);

                    piz_ZZ_r[i] += 12.0 * (rhoBz_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -rhoBz_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -rhoBz_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -rhoBz_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +rhoBz_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +rhoBz_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piz_ZZ_r[i] += 12.0 * (gradBz_z_r[i] * rhoCy_r[i] * gradBy_z_r[i]
                                -gradBz_z_i[i] * rhoCy_r[i] * gradBy_z_i[i]
                                -gradBz_z_r[i] * rhoCy_i[i] * gradBy_z_i[i]
                                -gradBz_z_i[i] * rhoCy_i[i] * gradBy_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCy_i[i] * gradBy_z_i[i]
                                +gradBz_z_i[i] * rhoCy_r[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * rhoCy_i[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * rhoCy_r[i] * gradBy_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_z_r[i] * rhoBy_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * rhoBy_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * rhoBy_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * rhoBy_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * rhoBy_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * rhoBy_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * rhoBy_i[i]);

                    piz_ZZ_r[i] += 6.0 * (rhoBy_r[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                -rhoBy_i[i] * gradCz_z_r[i] * gradBy_z_i[i]
                                -rhoBy_r[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                -rhoBy_i[i] * gradCz_z_i[i] * gradBy_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -rhoBy_i[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                +rhoBy_i[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_z_i[i] * gradBy_z_r[i]
                                +rhoBy_r[i] * gradCz_z_r[i] * gradBy_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBy_z_r[i] * rhoCz_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * rhoCz_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * rhoCz_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * rhoCz_i[i] * gradBy_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * rhoCz_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * rhoCz_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCz_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * rhoCz_r[i] * gradBy_z_i[i]);

                    piz_ZZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_ZZ_r[i] += 12.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_ZZ_r[i] += 12.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_ZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_z_r[i] * rhoBz_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * rhoBz_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * rhoBz_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * rhoBz_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * rhoBz_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * rhoBz_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * rhoBz_i[i]);

                    piz_ZZ_r[i] += 6.0 * (rhoBz_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -rhoBz_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -rhoBz_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -rhoBz_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -rhoBz_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +rhoBz_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +rhoBz_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_ZZ_r[i] += 6.0 * (gradBz_z_r[i] * rhoCz_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * rhoCz_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * rhoCz_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * rhoCz_i[i] * gradBz_z_r[i]);

                    piz_ZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * rhoCz_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * rhoCz_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCz_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * rhoCz_r[i] * gradBz_z_i[i]);

                    pix_XXX_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_XXX_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_XXX_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_XXX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_XXX_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    pix_XXX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    pix_XXX_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_x_r[i]);

                    pix_XXX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBy_x_i[i]);

                    pix_XXX_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    pix_XXX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    pix_XXX_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_x_r[i]);

                    pix_XXX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBz_x_i[i]);

                    pix_XXY_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_XXY_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_XXY_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_XXY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_XXY_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    pix_XXY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    pix_XXY_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_y_r[i]);

                    pix_XXY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBy_y_i[i]);

                    pix_XXY_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    pix_XXY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    pix_XXY_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_y_r[i]);

                    pix_XXY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBz_y_i[i]);

                    pix_XXZ_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_XXZ_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_XXZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_XXZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_XXZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    pix_XXZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    pix_XXZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_z_r[i]);

                    pix_XXZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBy_z_i[i]);

                    pix_XXZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    pix_XXZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    pix_XXZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_z_r[i]);

                    pix_XXZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBz_z_i[i]);

                    pix_XYX_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_XYX_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_XYX_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_XYX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_XYX_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    pix_XYX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    pix_XYX_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_x_r[i]);

                    pix_XYX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBy_x_i[i]);

                    pix_XYX_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    pix_XYX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    pix_XYX_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_x_r[i]);

                    pix_XYX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBz_x_i[i]);

                    pix_XYY_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_XYY_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_XYY_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_XYY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_XYY_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    pix_XYY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    pix_XYY_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_y_r[i]);

                    pix_XYY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBy_y_i[i]);

                    pix_XYY_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    pix_XYY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    pix_XYY_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_y_r[i]);

                    pix_XYY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBz_y_i[i]);

                    pix_XYZ_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_XYZ_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_XYZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_XYZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_XYZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    pix_XYZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    pix_XYZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_z_r[i]);

                    pix_XYZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBy_z_i[i]);

                    pix_XYZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    pix_XYZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    pix_XYZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_z_r[i]);

                    pix_XYZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBz_z_i[i]);

                    pix_XZX_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_XZX_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_XZX_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_XZX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_XZX_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    pix_XZX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    pix_XZX_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_x_r[i]);

                    pix_XZX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBy_x_i[i]);

                    pix_XZX_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    pix_XZX_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    pix_XZX_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_x_r[i]);

                    pix_XZX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBz_x_i[i]);

                    pix_XZY_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_XZY_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_XZY_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_XZY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_XZY_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    pix_XZY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    pix_XZY_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_y_r[i]);

                    pix_XZY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBy_y_i[i]);

                    pix_XZY_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    pix_XZY_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    pix_XZY_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_y_r[i]);

                    pix_XZY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBz_y_i[i]);

                    pix_XZZ_r[i] = 12.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_XZZ_i[i] = 12.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_XZZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_XZZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_XZZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    pix_XZZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    pix_XZZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_z_r[i]);

                    pix_XZZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBy_z_i[i]);

                    pix_XZZ_r[i] += 12.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    pix_XZZ_i[i] += 12.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    pix_XZZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_z_r[i]);

                    pix_XZZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBz_z_i[i]);

                    pix_YXX_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_YXX_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_YXX_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_YXX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_YXX_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    pix_YXX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    pix_YXX_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_x_r[i]);

                    pix_YXX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBy_x_i[i]);

                    pix_YXX_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    pix_YXX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    pix_YXX_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_x_r[i]);

                    pix_YXX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBz_x_i[i]);

                    pix_YXY_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_YXY_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_YXY_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_YXY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_YXY_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    pix_YXY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    pix_YXY_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_y_r[i]);

                    pix_YXY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBy_y_i[i]);

                    pix_YXY_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    pix_YXY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    pix_YXY_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_y_r[i]);

                    pix_YXY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBz_y_i[i]);

                    pix_YXZ_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_YXZ_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_YXZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_YXZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_YXZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    pix_YXZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    pix_YXZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_z_r[i]);

                    pix_YXZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBy_z_i[i]);

                    pix_YXZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    pix_YXZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    pix_YXZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_z_r[i]);

                    pix_YXZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBz_z_i[i]);

                    pix_YYX_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_YYX_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_YYX_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_YYX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_YYX_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    pix_YYX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    pix_YYX_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_x_r[i]);

                    pix_YYX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBy_x_i[i]);

                    pix_YYX_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    pix_YYX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    pix_YYX_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_x_r[i]);

                    pix_YYX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBz_x_i[i]);

                    pix_YYY_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_YYY_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_YYY_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_YYY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_YYY_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    pix_YYY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    pix_YYY_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_y_r[i]);

                    pix_YYY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBy_y_i[i]);

                    pix_YYY_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    pix_YYY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    pix_YYY_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_y_r[i]);

                    pix_YYY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBz_y_i[i]);

                    pix_YYZ_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_YYZ_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_YYZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_YYZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_YYZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    pix_YYZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    pix_YYZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_z_r[i]);

                    pix_YYZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBy_z_i[i]);

                    pix_YYZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    pix_YYZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    pix_YYZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_z_r[i]);

                    pix_YYZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBz_z_i[i]);

                    pix_YZX_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_YZX_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_YZX_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_YZX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_YZX_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    pix_YZX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    pix_YZX_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_x_r[i]);

                    pix_YZX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBy_x_i[i]);

                    pix_YZX_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    pix_YZX_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    pix_YZX_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_x_r[i]);

                    pix_YZX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBz_x_i[i]);

                    pix_YZY_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_YZY_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_YZY_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_YZY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_YZY_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    pix_YZY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    pix_YZY_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_y_r[i]);

                    pix_YZY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBy_y_i[i]);

                    pix_YZY_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    pix_YZY_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    pix_YZY_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_y_r[i]);

                    pix_YZY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBz_y_i[i]);

                    pix_YZZ_r[i] = 12.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_YZZ_i[i] = 12.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_YZZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_YZZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_YZZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    pix_YZZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    pix_YZZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_z_r[i]);

                    pix_YZZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBy_z_i[i]);

                    pix_YZZ_r[i] += 12.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    pix_YZZ_i[i] += 12.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    pix_YZZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_z_r[i]);

                    pix_YZZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBz_z_i[i]);

                    pix_ZXX_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_ZXX_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_ZXX_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    pix_ZXX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    pix_ZXX_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    pix_ZXX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    pix_ZXX_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_x_r[i]);

                    pix_ZXX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBy_x_i[i]);

                    pix_ZXX_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    pix_ZXX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    pix_ZXX_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_x_r[i]);

                    pix_ZXX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBz_x_i[i]);

                    pix_ZXY_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_ZXY_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_ZXY_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    pix_ZXY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    pix_ZXY_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    pix_ZXY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    pix_ZXY_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_y_r[i]);

                    pix_ZXY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBy_y_i[i]);

                    pix_ZXY_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    pix_ZXY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    pix_ZXY_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_y_r[i]);

                    pix_ZXY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBz_y_i[i]);

                    pix_ZXZ_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_ZXZ_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_ZXZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    pix_ZXZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    pix_ZXZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    pix_ZXZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    pix_ZXZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_z_r[i]);

                    pix_ZXZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBy_z_i[i]);

                    pix_ZXZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    pix_ZXZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    pix_ZXZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_z_r[i]);

                    pix_ZXZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBz_z_i[i]);

                    pix_ZYX_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_ZYX_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_ZYX_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    pix_ZYX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    pix_ZYX_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    pix_ZYX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    pix_ZYX_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_x_r[i]);

                    pix_ZYX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBy_x_i[i]);

                    pix_ZYX_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    pix_ZYX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    pix_ZYX_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_x_r[i]);

                    pix_ZYX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBz_x_i[i]);

                    pix_ZYY_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_ZYY_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_ZYY_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    pix_ZYY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    pix_ZYY_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    pix_ZYY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    pix_ZYY_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_y_r[i]);

                    pix_ZYY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBy_y_i[i]);

                    pix_ZYY_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    pix_ZYY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    pix_ZYY_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_y_r[i]);

                    pix_ZYY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBz_y_i[i]);

                    pix_ZYZ_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_ZYZ_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_ZYZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    pix_ZYZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    pix_ZYZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    pix_ZYZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    pix_ZYZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_z_r[i]);

                    pix_ZYZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBy_z_i[i]);

                    pix_ZYZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    pix_ZYZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    pix_ZYZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_z_r[i]);

                    pix_ZYZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBz_z_i[i]);

                    pix_ZZX_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_ZZX_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_ZZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    pix_ZZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    pix_ZZX_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    pix_ZZX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    pix_ZZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_x_r[i]);

                    pix_ZZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBy_x_i[i]);

                    pix_ZZX_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    pix_ZZX_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    pix_ZZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_x_r[i]);

                    pix_ZZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBz_x_i[i]);

                    pix_ZZY_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_ZZY_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_ZZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    pix_ZZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    pix_ZZY_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    pix_ZZY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    pix_ZZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_y_r[i]);

                    pix_ZZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBy_y_i[i]);

                    pix_ZZY_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    pix_ZZY_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    pix_ZZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_y_r[i]);

                    pix_ZZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBz_y_i[i]);

                    pix_ZZZ_r[i] = 12.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_ZZZ_i[i] = 12.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_ZZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    pix_ZZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    pix_ZZZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    pix_ZZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    pix_ZZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_z_r[i]);

                    pix_ZZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBy_z_i[i]);

                    pix_ZZZ_r[i] += 12.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    pix_ZZZ_i[i] += 12.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    pix_ZZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_z_r[i]);

                    pix_ZZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBz_z_i[i]);

                    piy_XXX_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piy_XXX_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piy_XXX_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_x_r[i]);

                    piy_XXX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBx_x_i[i]);

                    piy_XXX_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_XXX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_XXX_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_XXX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_XXX_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piy_XXX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piy_XXX_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_x_r[i]);

                    piy_XXX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBz_x_i[i]);

                    piy_XXY_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piy_XXY_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piy_XXY_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_y_r[i]);

                    piy_XXY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBx_y_i[i]);

                    piy_XXY_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_XXY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_XXY_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_XXY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_XXY_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piy_XXY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piy_XXY_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_y_r[i]);

                    piy_XXY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBz_y_i[i]);

                    piy_XXZ_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBy_x_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBy_x_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piy_XXZ_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBy_x_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piy_XXZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCy_x_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_z_r[i]);

                    piy_XXZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_x_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_x_r[i] * gradBx_z_i[i]);

                    piy_XXZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_XXZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_XXZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_XXZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_XXZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piy_XXZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piy_XXZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_z_r[i]);

                    piy_XXZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBz_z_i[i]);

                    piy_XYX_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piy_XYX_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piy_XYX_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_x_r[i]);

                    piy_XYX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBx_x_i[i]);

                    piy_XYX_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_XYX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_XYX_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_XYX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_XYX_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piy_XYX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piy_XYX_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_x_r[i]);

                    piy_XYX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBz_x_i[i]);

                    piy_XYY_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piy_XYY_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piy_XYY_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_y_r[i]);

                    piy_XYY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBx_y_i[i]);

                    piy_XYY_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_XYY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_XYY_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_XYY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_XYY_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piy_XYY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piy_XYY_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_y_r[i]);

                    piy_XYY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBz_y_i[i]);

                    piy_XYZ_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBy_x_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBy_x_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piy_XYZ_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBy_x_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piy_XYZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCy_y_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_z_r[i]);

                    piy_XYZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_y_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_y_r[i] * gradBx_z_i[i]);

                    piy_XYZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_XYZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_XYZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_XYZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_XYZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piy_XYZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piy_XYZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_z_r[i]);

                    piy_XYZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBz_z_i[i]);

                    piy_XZX_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piy_XZX_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piy_XZX_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_x_r[i]);

                    piy_XZX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBx_x_i[i]);

                    piy_XZX_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_XZX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_XZX_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_XZX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_XZX_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piy_XZX_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piy_XZX_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_x_r[i]);

                    piy_XZX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBz_x_i[i]);

                    piy_XZY_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piy_XZY_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piy_XZY_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_y_r[i]);

                    piy_XZY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBx_y_i[i]);

                    piy_XZY_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_XZY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_XZY_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_XZY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_XZY_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piy_XZY_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piy_XZY_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_y_r[i]);

                    piy_XZY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBz_y_i[i]);

                    piy_XZZ_r[i] = 12.0 * (gradBy_x_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBy_x_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBy_x_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piy_XZZ_i[i] = 12.0 * ( -gradBy_x_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBy_x_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBy_x_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piy_XZZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCy_z_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_z_r[i]);

                    piy_XZZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_z_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCy_z_r[i] * gradBx_z_i[i]);

                    piy_XZZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_XZZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_XZZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_XZZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_XZZ_r[i] += 12.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piy_XZZ_i[i] += 12.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piy_XZZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_z_r[i]);

                    piy_XZZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBz_z_i[i]);

                    piy_YXX_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piy_YXX_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piy_YXX_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_x_r[i]);

                    piy_YXX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBx_x_i[i]);

                    piy_YXX_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_YXX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_YXX_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_YXX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_YXX_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piy_YXX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piy_YXX_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_x_r[i]);

                    piy_YXX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBz_x_i[i]);

                    piy_YXY_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piy_YXY_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piy_YXY_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_y_r[i]);

                    piy_YXY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBx_y_i[i]);

                    piy_YXY_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_YXY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_YXY_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_YXY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_YXY_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piy_YXY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piy_YXY_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_y_r[i]);

                    piy_YXY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBz_y_i[i]);

                    piy_YXZ_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBy_y_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBy_y_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piy_YXZ_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBy_y_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piy_YXZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCy_x_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_z_r[i]);

                    piy_YXZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_x_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_x_r[i] * gradBx_z_i[i]);

                    piy_YXZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_YXZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_YXZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_YXZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_YXZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piy_YXZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piy_YXZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_z_r[i]);

                    piy_YXZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBz_z_i[i]);

                    piy_YYX_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piy_YYX_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piy_YYX_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_x_r[i]);

                    piy_YYX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBx_x_i[i]);

                    piy_YYX_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_YYX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_YYX_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_YYX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_YYX_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piy_YYX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piy_YYX_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_x_r[i]);

                    piy_YYX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBz_x_i[i]);

                    piy_YYY_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piy_YYY_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piy_YYY_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_y_r[i]);

                    piy_YYY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBx_y_i[i]);

                    piy_YYY_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_YYY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_YYY_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_YYY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_YYY_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piy_YYY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piy_YYY_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_y_r[i]);

                    piy_YYY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBz_y_i[i]);

                    piy_YYZ_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBy_y_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBy_y_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piy_YYZ_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBy_y_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piy_YYZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCy_y_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_z_r[i]);

                    piy_YYZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_y_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_y_r[i] * gradBx_z_i[i]);

                    piy_YYZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_YYZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_YYZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_YYZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_YYZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piy_YYZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piy_YYZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_z_r[i]);

                    piy_YYZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBz_z_i[i]);

                    piy_YZX_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piy_YZX_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piy_YZX_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_x_r[i]);

                    piy_YZX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBx_x_i[i]);

                    piy_YZX_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_YZX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_YZX_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_YZX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_YZX_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piy_YZX_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piy_YZX_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_x_r[i]);

                    piy_YZX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBz_x_i[i]);

                    piy_YZY_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piy_YZY_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piy_YZY_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_y_r[i]);

                    piy_YZY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBx_y_i[i]);

                    piy_YZY_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_YZY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_YZY_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_YZY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_YZY_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piy_YZY_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piy_YZY_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_y_r[i]);

                    piy_YZY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBz_y_i[i]);

                    piy_YZZ_r[i] = 12.0 * (gradBy_y_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBy_y_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBy_y_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piy_YZZ_i[i] = 12.0 * ( -gradBy_y_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBy_y_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBy_y_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piy_YZZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCy_z_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_z_r[i]);

                    piy_YZZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_z_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCy_z_r[i] * gradBx_z_i[i]);

                    piy_YZZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_YZZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_YZZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_YZZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_YZZ_r[i] += 12.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piy_YZZ_i[i] += 12.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piy_YZZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_z_r[i]);

                    piy_YZZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBz_z_i[i]);

                    piy_ZXX_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piy_ZXX_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piy_ZXX_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_x_r[i]);

                    piy_ZXX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBx_x_i[i]);

                    piy_ZXX_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_ZXX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_ZXX_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piy_ZXX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piy_ZXX_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piy_ZXX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piy_ZXX_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_x_r[i]);

                    piy_ZXX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBz_x_i[i]);

                    piy_ZXY_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piy_ZXY_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piy_ZXY_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_y_r[i]);

                    piy_ZXY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBx_y_i[i]);

                    piy_ZXY_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_ZXY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_ZXY_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piy_ZXY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piy_ZXY_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piy_ZXY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piy_ZXY_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_y_r[i]);

                    piy_ZXY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBz_y_i[i]);

                    piy_ZXZ_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBy_z_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBy_z_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piy_ZXZ_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBy_z_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piy_ZXZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCy_x_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_z_r[i]);

                    piy_ZXZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_x_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCy_x_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_x_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_x_r[i] * gradBx_z_i[i]);

                    piy_ZXZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_ZXZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_ZXZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piy_ZXZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piy_ZXZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piy_ZXZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piy_ZXZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_z_r[i]);

                    piy_ZXZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBz_z_i[i]);

                    piy_ZYX_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piy_ZYX_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piy_ZYX_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_x_r[i]);

                    piy_ZYX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBx_x_i[i]);

                    piy_ZYX_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_ZYX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_ZYX_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piy_ZYX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piy_ZYX_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piy_ZYX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piy_ZYX_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_x_r[i]);

                    piy_ZYX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBz_x_i[i]);

                    piy_ZYY_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piy_ZYY_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piy_ZYY_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_y_r[i]);

                    piy_ZYY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBx_y_i[i]);

                    piy_ZYY_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_ZYY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_ZYY_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piy_ZYY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piy_ZYY_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piy_ZYY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piy_ZYY_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_y_r[i]);

                    piy_ZYY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBz_y_i[i]);

                    piy_ZYZ_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBy_z_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBy_z_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piy_ZYZ_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBy_z_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piy_ZYZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCy_y_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_z_r[i]);

                    piy_ZYZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_y_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCy_y_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_y_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_y_r[i] * gradBx_z_i[i]);

                    piy_ZYZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_ZYZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_ZYZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piy_ZYZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piy_ZYZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piy_ZYZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piy_ZYZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_z_r[i]);

                    piy_ZYZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBz_z_i[i]);

                    piy_ZZX_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piy_ZZX_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piy_ZZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_x_r[i]);

                    piy_ZZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBx_x_i[i]);

                    piy_ZZX_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_ZZX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_ZZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piy_ZZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piy_ZZX_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piy_ZZX_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piy_ZZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_x_r[i]);

                    piy_ZZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBz_x_i[i]);

                    piy_ZZY_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piy_ZZY_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piy_ZZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_y_r[i]);

                    piy_ZZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBx_y_i[i]);

                    piy_ZZY_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_ZZY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_ZZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piy_ZZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piy_ZZY_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piy_ZZY_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piy_ZZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_y_r[i]);

                    piy_ZZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBz_y_i[i]);

                    piy_ZZZ_r[i] = 12.0 * (gradBy_z_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBy_z_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBy_z_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piy_ZZZ_i[i] = 12.0 * ( -gradBy_z_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBy_z_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBy_z_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piy_ZZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCy_z_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_z_r[i]);

                    piy_ZZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCy_z_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCy_z_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_z_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCy_z_r[i] * gradBx_z_i[i]);

                    piy_ZZZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_ZZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_ZZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piy_ZZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piy_ZZZ_r[i] += 12.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piy_ZZZ_i[i] += 12.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piy_ZZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_z_r[i]);

                    piy_ZZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBz_z_i[i]);

                    piz_XXX_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piz_XXX_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piz_XXX_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_x_r[i]);

                    piz_XXX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBx_x_i[i]);

                    piz_XXX_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piz_XXX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piz_XXX_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_x_r[i]);

                    piz_XXX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBy_x_i[i]);

                    piz_XXX_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_XXX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_XXX_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_XXX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_XXY_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piz_XXY_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piz_XXY_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_y_r[i]);

                    piz_XXY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBx_y_i[i]);

                    piz_XXY_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piz_XXY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piz_XXY_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_y_r[i]);

                    piz_XXY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBy_y_i[i]);

                    piz_XXY_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_XXY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_XXY_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_XXY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_XXZ_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBz_x_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBz_x_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piz_XXZ_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBz_x_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piz_XXZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCz_x_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_z_r[i]);

                    piz_XXZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_x_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_x_r[i] * gradBx_z_i[i]);

                    piz_XXZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBz_x_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBz_x_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piz_XXZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBz_x_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piz_XXZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCz_x_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_z_r[i]);

                    piz_XXZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_x_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_x_r[i] * gradBy_z_i[i]);

                    piz_XXZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_XXZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_XXZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_XXZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_XYX_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piz_XYX_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piz_XYX_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_x_r[i]);

                    piz_XYX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBx_x_i[i]);

                    piz_XYX_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piz_XYX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piz_XYX_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_x_r[i]);

                    piz_XYX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBy_x_i[i]);

                    piz_XYX_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_XYX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_XYX_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_XYX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_XYY_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piz_XYY_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piz_XYY_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_y_r[i]);

                    piz_XYY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBx_y_i[i]);

                    piz_XYY_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piz_XYY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piz_XYY_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_y_r[i]);

                    piz_XYY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBy_y_i[i]);

                    piz_XYY_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_XYY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_XYY_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_XYY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_XYZ_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBz_x_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBz_x_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piz_XYZ_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBz_x_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piz_XYZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCz_y_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_z_r[i]);

                    piz_XYZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_y_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_y_r[i] * gradBx_z_i[i]);

                    piz_XYZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBz_x_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBz_x_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piz_XYZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBz_x_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piz_XYZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCz_y_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_z_r[i]);

                    piz_XYZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_y_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_y_r[i] * gradBy_z_i[i]);

                    piz_XYZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_XYZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_XYZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_XYZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_XZX_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piz_XZX_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piz_XZX_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBx_x_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_x_r[i]);

                    piz_XZX_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBx_x_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBx_x_i[i]);

                    piz_XZX_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piz_XZX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piz_XZX_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBy_x_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_x_r[i]);

                    piz_XZX_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBy_x_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBy_x_i[i]);

                    piz_XZX_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_XZX_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_XZX_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_XZX_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_XZY_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piz_XZY_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piz_XZY_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBx_y_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_y_r[i]);

                    piz_XZY_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBx_y_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBx_y_i[i]);

                    piz_XZY_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piz_XZY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piz_XZY_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBy_y_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_y_r[i]);

                    piz_XZY_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBy_y_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBy_y_i[i]);

                    piz_XZY_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_XZY_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_XZY_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_XZY_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_XZZ_r[i] = 12.0 * (gradBz_x_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBz_x_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBz_x_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piz_XZZ_i[i] = 12.0 * ( -gradBz_x_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBz_x_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBz_x_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piz_XZZ_r[i] += 6.0 * (gradBx_x_r[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                -gradBx_x_i[i] * gradCz_z_r[i] * gradBx_z_i[i]
                                -gradBx_x_r[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_z_r[i]);

                    piz_XZZ_i[i] += 6.0 * ( -gradBx_x_i[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                +gradBx_x_i[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_z_i[i] * gradBx_z_r[i]
                                +gradBx_x_r[i] * gradCz_z_r[i] * gradBx_z_i[i]);

                    piz_XZZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBz_x_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBz_x_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piz_XZZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBz_x_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBz_x_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piz_XZZ_r[i] += 6.0 * (gradBy_x_r[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                -gradBy_x_i[i] * gradCz_z_r[i] * gradBy_z_i[i]
                                -gradBy_x_r[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_z_r[i]);

                    piz_XZZ_i[i] += 6.0 * ( -gradBy_x_i[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                +gradBy_x_i[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_z_i[i] * gradBy_z_r[i]
                                +gradBy_x_r[i] * gradCz_z_r[i] * gradBy_z_i[i]);

                    piz_XZZ_r[i] += 12.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_XZZ_i[i] += 12.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_XZZ_r[i] += 6.0 * (gradBz_x_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_x_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_x_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_XZZ_i[i] += 6.0 * ( -gradBz_x_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_x_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_x_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_YXX_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piz_YXX_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piz_YXX_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_x_r[i]);

                    piz_YXX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBx_x_i[i]);

                    piz_YXX_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piz_YXX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piz_YXX_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_x_r[i]);

                    piz_YXX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBy_x_i[i]);

                    piz_YXX_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_YXX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_YXX_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_YXX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_YXY_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piz_YXY_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piz_YXY_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_y_r[i]);

                    piz_YXY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBx_y_i[i]);

                    piz_YXY_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piz_YXY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piz_YXY_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_y_r[i]);

                    piz_YXY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBy_y_i[i]);

                    piz_YXY_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_YXY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_YXY_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_YXY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_YXZ_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBz_y_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBz_y_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piz_YXZ_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBz_y_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piz_YXZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCz_x_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_z_r[i]);

                    piz_YXZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_x_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_x_r[i] * gradBx_z_i[i]);

                    piz_YXZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBz_y_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBz_y_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piz_YXZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBz_y_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piz_YXZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCz_x_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_z_r[i]);

                    piz_YXZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_x_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_x_r[i] * gradBy_z_i[i]);

                    piz_YXZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_YXZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_YXZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_YXZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_YYX_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piz_YYX_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piz_YYX_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_x_r[i]);

                    piz_YYX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBx_x_i[i]);

                    piz_YYX_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piz_YYX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piz_YYX_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_x_r[i]);

                    piz_YYX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBy_x_i[i]);

                    piz_YYX_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_YYX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_YYX_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_YYX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_YYY_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piz_YYY_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piz_YYY_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_y_r[i]);

                    piz_YYY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBx_y_i[i]);

                    piz_YYY_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piz_YYY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piz_YYY_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_y_r[i]);

                    piz_YYY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBy_y_i[i]);

                    piz_YYY_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_YYY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_YYY_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_YYY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_YYZ_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBz_y_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBz_y_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piz_YYZ_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBz_y_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piz_YYZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCz_y_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_z_r[i]);

                    piz_YYZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_y_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_y_r[i] * gradBx_z_i[i]);

                    piz_YYZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBz_y_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBz_y_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piz_YYZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBz_y_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piz_YYZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCz_y_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_z_r[i]);

                    piz_YYZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_y_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_y_r[i] * gradBy_z_i[i]);

                    piz_YYZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_YYZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_YYZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_YYZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_YZX_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piz_YZX_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piz_YZX_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBx_x_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_x_r[i]);

                    piz_YZX_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBx_x_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBx_x_i[i]);

                    piz_YZX_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piz_YZX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piz_YZX_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBy_x_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_x_r[i]);

                    piz_YZX_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBy_x_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBy_x_i[i]);

                    piz_YZX_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_YZX_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_YZX_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_YZX_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_YZY_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piz_YZY_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piz_YZY_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBx_y_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_y_r[i]);

                    piz_YZY_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBx_y_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBx_y_i[i]);

                    piz_YZY_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piz_YZY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piz_YZY_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBy_y_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_y_r[i]);

                    piz_YZY_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBy_y_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBy_y_i[i]);

                    piz_YZY_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_YZY_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_YZY_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_YZY_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_YZZ_r[i] = 12.0 * (gradBz_y_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBz_y_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBz_y_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piz_YZZ_i[i] = 12.0 * ( -gradBz_y_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBz_y_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBz_y_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piz_YZZ_r[i] += 6.0 * (gradBx_y_r[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                -gradBx_y_i[i] * gradCz_z_r[i] * gradBx_z_i[i]
                                -gradBx_y_r[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_z_r[i]);

                    piz_YZZ_i[i] += 6.0 * ( -gradBx_y_i[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                +gradBx_y_i[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_z_i[i] * gradBx_z_r[i]
                                +gradBx_y_r[i] * gradCz_z_r[i] * gradBx_z_i[i]);

                    piz_YZZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBz_y_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBz_y_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piz_YZZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBz_y_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBz_y_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piz_YZZ_r[i] += 6.0 * (gradBy_y_r[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                -gradBy_y_i[i] * gradCz_z_r[i] * gradBy_z_i[i]
                                -gradBy_y_r[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_z_r[i]);

                    piz_YZZ_i[i] += 6.0 * ( -gradBy_y_i[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                +gradBy_y_i[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_z_i[i] * gradBy_z_r[i]
                                +gradBy_y_r[i] * gradCz_z_r[i] * gradBy_z_i[i]);

                    piz_YZZ_r[i] += 12.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_YZZ_i[i] += 12.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_YZZ_r[i] += 6.0 * (gradBz_y_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_y_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_y_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_YZZ_i[i] += 6.0 * ( -gradBz_y_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_y_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_y_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_ZXX_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBx_x_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_x_r[i]);

                    piz_ZXX_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_x_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBx_x_i[i]);

                    piz_ZXX_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_x_r[i]);

                    piz_ZXX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBx_x_i[i]);

                    piz_ZXX_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBy_x_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_x_r[i]);

                    piz_ZXX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_x_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBy_x_i[i]);

                    piz_ZXX_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_x_r[i]);

                    piz_ZXX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBy_x_i[i]);

                    piz_ZXX_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_ZXX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_ZXX_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_x_r[i]);

                    piz_ZXX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_x_i[i]);

                    piz_ZXY_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBx_y_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_y_r[i]);

                    piz_ZXY_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_y_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBx_y_i[i]);

                    piz_ZXY_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_y_r[i]);

                    piz_ZXY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBx_y_i[i]);

                    piz_ZXY_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBy_y_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_y_r[i]);

                    piz_ZXY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_y_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBy_y_i[i]);

                    piz_ZXY_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_y_r[i]);

                    piz_ZXY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBy_y_i[i]);

                    piz_ZXY_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_ZXY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_ZXY_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_y_r[i]);

                    piz_ZXY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_y_i[i]);

                    piz_ZXZ_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                -gradBz_z_i[i] * gradCx_x_r[i] * gradBx_z_i[i]
                                -gradBz_z_r[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_z_r[i]);

                    piz_ZXZ_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_x_i[i] * gradBx_z_i[i]
                                +gradBz_z_i[i] * gradCx_x_r[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_x_i[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_x_r[i] * gradBx_z_i[i]);

                    piz_ZXZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCz_x_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_z_r[i]);

                    piz_ZXZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_x_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCz_x_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_x_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_x_r[i] * gradBx_z_i[i]);

                    piz_ZXZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                -gradBz_z_i[i] * gradCy_x_r[i] * gradBy_z_i[i]
                                -gradBz_z_r[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_z_r[i]);

                    piz_ZXZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_x_i[i] * gradBy_z_i[i]
                                +gradBz_z_i[i] * gradCy_x_r[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_x_i[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_x_r[i] * gradBy_z_i[i]);

                    piz_ZXZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCz_x_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_z_r[i]);

                    piz_ZXZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_x_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCz_x_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_x_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_x_r[i] * gradBy_z_i[i]);

                    piz_ZXZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_ZXZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_ZXZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_x_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_z_r[i]);

                    piz_ZXZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_x_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_x_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_x_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_x_r[i] * gradBz_z_i[i]);

                    piz_ZYX_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBx_x_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_x_r[i]);

                    piz_ZYX_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_x_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBx_x_i[i]);

                    piz_ZYX_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_x_r[i]);

                    piz_ZYX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBx_x_i[i]);

                    piz_ZYX_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBy_x_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_x_r[i]);

                    piz_ZYX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_x_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBy_x_i[i]);

                    piz_ZYX_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_x_r[i]);

                    piz_ZYX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBy_x_i[i]);

                    piz_ZYX_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_ZYX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_ZYX_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_x_r[i]);

                    piz_ZYX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_x_i[i]);

                    piz_ZYY_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBx_y_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_y_r[i]);

                    piz_ZYY_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_y_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBx_y_i[i]);

                    piz_ZYY_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_y_r[i]);

                    piz_ZYY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBx_y_i[i]);

                    piz_ZYY_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBy_y_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_y_r[i]);

                    piz_ZYY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_y_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBy_y_i[i]);

                    piz_ZYY_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_y_r[i]);

                    piz_ZYY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBy_y_i[i]);

                    piz_ZYY_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_ZYY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_ZYY_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_y_r[i]);

                    piz_ZYY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_y_i[i]);

                    piz_ZYZ_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                -gradBz_z_i[i] * gradCx_y_r[i] * gradBx_z_i[i]
                                -gradBz_z_r[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_z_r[i]);

                    piz_ZYZ_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_y_i[i] * gradBx_z_i[i]
                                +gradBz_z_i[i] * gradCx_y_r[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_y_i[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_y_r[i] * gradBx_z_i[i]);

                    piz_ZYZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCz_y_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_z_r[i]);

                    piz_ZYZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_y_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCz_y_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_y_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_y_r[i] * gradBx_z_i[i]);

                    piz_ZYZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                -gradBz_z_i[i] * gradCy_y_r[i] * gradBy_z_i[i]
                                -gradBz_z_r[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_z_r[i]);

                    piz_ZYZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_y_i[i] * gradBy_z_i[i]
                                +gradBz_z_i[i] * gradCy_y_r[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_y_i[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_y_r[i] * gradBy_z_i[i]);

                    piz_ZYZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCz_y_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_z_r[i]);

                    piz_ZYZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_y_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCz_y_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_y_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_y_r[i] * gradBy_z_i[i]);

                    piz_ZYZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_ZYZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_ZYZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_y_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_z_r[i]);

                    piz_ZYZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_y_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_y_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_y_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_y_r[i] * gradBz_z_i[i]);

                    piz_ZZX_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBx_x_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_x_r[i]);

                    piz_ZZX_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_x_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBx_x_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBx_x_i[i]);

                    piz_ZZX_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBx_x_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_x_r[i]);

                    piz_ZZX_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_x_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBx_x_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBx_x_i[i]);

                    piz_ZZX_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBy_x_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_x_r[i]);

                    piz_ZZX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_x_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBy_x_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBy_x_i[i]);

                    piz_ZZX_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBy_x_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_x_r[i]);

                    piz_ZZX_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_x_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBy_x_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBy_x_i[i]);

                    piz_ZZX_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_ZZX_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_ZZX_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_x_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_x_r[i]);

                    piz_ZZX_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_x_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_x_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_x_i[i]);

                    piz_ZZY_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBx_y_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_y_r[i]);

                    piz_ZZY_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_y_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBx_y_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBx_y_i[i]);

                    piz_ZZY_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBx_y_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_y_r[i]);

                    piz_ZZY_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_y_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBx_y_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBx_y_i[i]);

                    piz_ZZY_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBy_y_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_y_r[i]);

                    piz_ZZY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_y_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBy_y_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBy_y_i[i]);

                    piz_ZZY_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBy_y_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_y_r[i]);

                    piz_ZZY_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_y_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBy_y_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBy_y_i[i]);

                    piz_ZZY_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_ZZY_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_ZZY_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_y_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_y_r[i]);

                    piz_ZZY_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_y_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_y_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_y_i[i]);

                    piz_ZZZ_r[i] = 12.0 * (gradBz_z_r[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                -gradBz_z_i[i] * gradCx_z_r[i] * gradBx_z_i[i]
                                -gradBz_z_r[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_z_r[i]);

                    piz_ZZZ_i[i] = 12.0 * ( -gradBz_z_i[i] * gradCx_z_i[i] * gradBx_z_i[i]
                                +gradBz_z_i[i] * gradCx_z_r[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_z_i[i] * gradBx_z_r[i]
                                +gradBz_z_r[i] * gradCx_z_r[i] * gradBx_z_i[i]);

                    piz_ZZZ_r[i] += 6.0 * (gradBx_z_r[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                -gradBx_z_i[i] * gradCz_z_r[i] * gradBx_z_i[i]
                                -gradBx_z_r[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_z_r[i]);

                    piz_ZZZ_i[i] += 6.0 * ( -gradBx_z_i[i] * gradCz_z_i[i] * gradBx_z_i[i]
                                +gradBx_z_i[i] * gradCz_z_r[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_z_i[i] * gradBx_z_r[i]
                                +gradBx_z_r[i] * gradCz_z_r[i] * gradBx_z_i[i]);

                    piz_ZZZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                -gradBz_z_i[i] * gradCy_z_r[i] * gradBy_z_i[i]
                                -gradBz_z_r[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_z_r[i]);

                    piz_ZZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCy_z_i[i] * gradBy_z_i[i]
                                +gradBz_z_i[i] * gradCy_z_r[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_z_i[i] * gradBy_z_r[i]
                                +gradBz_z_r[i] * gradCy_z_r[i] * gradBy_z_i[i]);

                    piz_ZZZ_r[i] += 6.0 * (gradBy_z_r[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                -gradBy_z_i[i] * gradCz_z_r[i] * gradBy_z_i[i]
                                -gradBy_z_r[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_z_r[i]);

                    piz_ZZZ_i[i] += 6.0 * ( -gradBy_z_i[i] * gradCz_z_i[i] * gradBy_z_i[i]
                                +gradBy_z_i[i] * gradCz_z_r[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_z_i[i] * gradBy_z_r[i]
                                +gradBy_z_r[i] * gradCz_z_r[i] * gradBy_z_i[i]);

                    piz_ZZZ_r[i] += 12.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_ZZZ_i[i] += 12.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    piz_ZZZ_r[i] += 6.0 * (gradBz_z_r[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                -gradBz_z_i[i] * gradCz_z_r[i] * gradBz_z_i[i]
                                -gradBz_z_r[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_z_r[i]);

                    piz_ZZZ_i[i] += 6.0 * ( -gradBz_z_i[i] * gradCz_z_i[i] * gradBz_z_i[i]
                                +gradBz_z_i[i] * gradCz_z_r[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_z_i[i] * gradBz_z_r[i]
                                +gradBz_z_r[i] * gradCz_z_r[i] * gradBz_z_i[i]);

                    gamx_r[i] = (D_sig_xx_r[i] * rhoBx_r[i]- D_sig_xx_i[i] * rhoBx_i[i]);

                    gamx_i[i] =  (D_sig_xx_i[i] * rhoBx_r[i]+ D_sig_xx_r[i] * rhoBx_i[i]);

                    gamx_r[i] += (D_lamtau_xx_r[i] * rhoCx_r[i]- D_lamtau_xx_i[i] * rhoCx_i[i]);

                    gamx_i[i] +=  (D_lamtau_xx_i[i] * rhoCx_r[i]+ D_lamtau_xx_r[i] * rhoCx_i[i]);

                    gamx_r[i] += (D_sig_xy_r[i] * rhoBy_r[i]- D_sig_xy_i[i] * rhoBy_i[i]);

                    gamx_i[i] +=  (D_sig_xy_i[i] * rhoBy_r[i]+ D_sig_xy_r[i] * rhoBy_i[i]);

                    gamx_r[i] += (D_lamtau_xy_r[i] * rhoCy_r[i]- D_lamtau_xy_i[i] * rhoCy_i[i]);

                    gamx_i[i] +=  (D_lamtau_xy_i[i] * rhoCy_r[i]+ D_lamtau_xy_r[i] * rhoCy_i[i]);

                    gamx_r[i] += (D_sig_xz_r[i] * rhoBz_r[i]- D_sig_xz_i[i] * rhoBz_i[i]);

                    gamx_i[i] +=  (D_sig_xz_i[i] * rhoBz_r[i]+ D_sig_xz_r[i] * rhoBz_i[i]);

                    gamx_r[i] += (D_lamtau_xz_r[i] * rhoCz_r[i]- D_lamtau_xz_i[i] * rhoCz_i[i]);

                    gamx_i[i] +=  (D_lamtau_xz_i[i] * rhoCz_r[i]+ D_lamtau_xz_r[i] * rhoCz_i[i]);

                    gamy_r[i] = (D_sig_yx_r[i] * rhoBx_r[i]- D_sig_yx_i[i] * rhoBx_i[i]);

                    gamy_i[i] =  (D_sig_yx_i[i] * rhoBx_r[i]+ D_sig_yx_r[i] * rhoBx_i[i]);

                    gamy_r[i] += (D_lamtau_yx_r[i] * rhoCx_r[i]- D_lamtau_yx_i[i] * rhoCx_i[i]);

                    gamy_i[i] +=  (D_lamtau_yx_i[i] * rhoCx_r[i]+ D_lamtau_yx_r[i] * rhoCx_i[i]);

                    gamy_r[i] += (D_sig_yy_r[i] * rhoBy_r[i]- D_sig_yy_i[i] * rhoBy_i[i]);

                    gamy_i[i] +=  (D_sig_yy_i[i] * rhoBy_r[i]+ D_sig_yy_r[i] * rhoBy_i[i]);

                    gamy_r[i] += (D_lamtau_yy_r[i] * rhoCy_r[i]- D_lamtau_yy_i[i] * rhoCy_i[i]);

                    gamy_i[i] +=  (D_lamtau_yy_i[i] * rhoCy_r[i]+ D_lamtau_yy_r[i] * rhoCy_i[i]);

                    gamy_r[i] += (D_sig_yz_r[i] * rhoBz_r[i]- D_sig_yz_i[i] * rhoBz_i[i]);

                    gamy_i[i] +=  (D_sig_yz_i[i] * rhoBz_r[i]+ D_sig_yz_r[i] * rhoBz_i[i]);

                    gamy_r[i] += (D_lamtau_yz_r[i] * rhoCz_r[i]- D_lamtau_yz_i[i] * rhoCz_i[i]);

                    gamy_i[i] +=  (D_lamtau_yz_i[i] * rhoCz_r[i]+ D_lamtau_yz_r[i] * rhoCz_i[i]);

                    gamz_r[i] = (D_sig_zx_r[i] * rhoBx_r[i]- D_sig_zx_i[i] * rhoBx_i[i]);

                    gamz_i[i] =  (D_sig_zx_i[i] * rhoBx_r[i]+ D_sig_zx_r[i] * rhoBx_i[i]);

                    gamz_r[i] += (D_lamtau_zx_r[i] * rhoCx_r[i]- D_lamtau_zx_i[i] * rhoCx_i[i]);

                    gamz_i[i] +=  (D_lamtau_zx_i[i] * rhoCx_r[i]+ D_lamtau_zx_r[i] * rhoCx_i[i]);

                    gamz_r[i] += (D_sig_zy_r[i] * rhoBy_r[i]- D_sig_zy_i[i] * rhoBy_i[i]);

                    gamz_i[i] +=  (D_sig_zy_i[i] * rhoBy_r[i]+ D_sig_zy_r[i] * rhoBy_i[i]);

                    gamz_r[i] += (D_lamtau_zy_r[i] * rhoCy_r[i]- D_lamtau_zy_i[i] * rhoCy_i[i]);

                    gamz_i[i] +=  (D_lamtau_zy_i[i] * rhoCy_r[i]+ D_lamtau_zy_r[i] * rhoCy_i[i]);

                    gamz_r[i] += (D_sig_zz_r[i] * rhoBz_r[i]- D_sig_zz_i[i] * rhoBz_i[i]);

                    gamz_i[i] +=  (D_sig_zz_i[i] * rhoBz_r[i]+ D_sig_zz_r[i] * rhoBz_i[i]);

                    gamz_r[i] += (D_lamtau_zz_r[i] * rhoCz_r[i]- D_lamtau_zz_i[i] * rhoCz_i[i]);

                    gamz_i[i] +=  (D_lamtau_zz_i[i] * rhoCz_r[i]+ D_lamtau_zz_r[i] * rhoCz_i[i]);

                    gamx_X_r[i] = (D_sig_xx_x_r[i] * rhoBx_r[i]- D_sig_xx_x_i[i] * rhoBx_i[i]);

                    gamx_X_i[i] =  (D_sig_xx_x_i[i] * rhoBx_r[i]+ D_sig_xx_x_r[i] * rhoBx_i[i]);

                    gamx_X_r[i] += (D_sig_xx_r[i] * gradBx_x_r[i]- D_sig_xx_i[i] * gradBx_x_i[i]);

                    gamx_X_i[i] +=  (D_sig_xx_i[i] * gradBx_x_r[i]+ D_sig_xx_r[i] * gradBx_x_i[i]);

                    gamx_X_r[i] += (D_lamtau_xx_x_r[i] * rhoCx_r[i]- D_lamtau_xx_x_i[i] * rhoCx_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xx_x_i[i] * rhoCx_r[i]+ D_lamtau_xx_x_r[i] * rhoCx_i[i]);

                    gamx_X_r[i] += (D_lamtau_xx_r[i] * gradCx_x_r[i]- D_lamtau_xx_i[i] * gradCx_x_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xx_i[i] * gradCx_x_r[i]+ D_lamtau_xx_r[i] * gradCx_x_i[i]);

                    gamx_X_r[i] += (D_sig_xy_x_r[i] * rhoBy_r[i]- D_sig_xy_x_i[i] * rhoBy_i[i]);

                    gamx_X_i[i] +=  (D_sig_xy_x_i[i] * rhoBy_r[i]+ D_sig_xy_x_r[i] * rhoBy_i[i]);

                    gamx_X_r[i] += (D_sig_xy_r[i] * gradBy_x_r[i]- D_sig_xy_i[i] * gradBy_x_i[i]);

                    gamx_X_i[i] +=  (D_sig_xy_i[i] * gradBy_x_r[i]+ D_sig_xy_r[i] * gradBy_x_i[i]);

                    gamx_X_r[i] += (D_lamtau_xy_x_r[i] * rhoCy_r[i]- D_lamtau_xy_x_i[i] * rhoCy_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xy_x_i[i] * rhoCy_r[i]+ D_lamtau_xy_x_r[i] * rhoCy_i[i]);

                    gamx_X_r[i] += (D_lamtau_xy_r[i] * gradCy_x_r[i]- D_lamtau_xy_i[i] * gradCy_x_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xy_i[i] * gradCy_x_r[i]+ D_lamtau_xy_r[i] * gradCy_x_i[i]);

                    gamx_X_r[i] += (D_sig_xz_x_r[i] * rhoBz_r[i]- D_sig_xz_x_i[i] * rhoBz_i[i]);

                    gamx_X_i[i] +=  (D_sig_xz_x_i[i] * rhoBz_r[i]+ D_sig_xz_x_r[i] * rhoBz_i[i]);

                    gamx_X_r[i] += (D_sig_xz_r[i] * gradBz_x_r[i]- D_sig_xz_i[i] * gradBz_x_i[i]);

                    gamx_X_i[i] +=  (D_sig_xz_i[i] * gradBz_x_r[i]+ D_sig_xz_r[i] * gradBz_x_i[i]);

                    gamx_X_r[i] += (D_lamtau_xz_x_r[i] * rhoCz_r[i]- D_lamtau_xz_x_i[i] * rhoCz_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xz_x_i[i] * rhoCz_r[i]+ D_lamtau_xz_x_r[i] * rhoCz_i[i]);

                    gamx_X_r[i] += (D_lamtau_xz_r[i] * gradCz_x_r[i]- D_lamtau_xz_i[i] * gradCz_x_i[i]);

                    gamx_X_i[i] +=  (D_lamtau_xz_i[i] * gradCz_x_r[i]+ D_lamtau_xz_r[i] * gradCz_x_i[i]);

                    gamx_Y_r[i] = (D_sig_xx_y_r[i] * rhoBx_r[i]- D_sig_xx_y_i[i] * rhoBx_i[i]);

                    gamx_Y_i[i] =  (D_sig_xx_y_i[i] * rhoBx_r[i]+ D_sig_xx_y_r[i] * rhoBx_i[i]);

                    gamx_Y_r[i] += (D_sig_xx_r[i] * gradBx_y_r[i]- D_sig_xx_i[i] * gradBx_y_i[i]);

                    gamx_Y_i[i] +=  (D_sig_xx_i[i] * gradBx_y_r[i]+ D_sig_xx_r[i] * gradBx_y_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xx_y_r[i] * rhoCx_r[i]- D_lamtau_xx_y_i[i] * rhoCx_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xx_y_i[i] * rhoCx_r[i]+ D_lamtau_xx_y_r[i] * rhoCx_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xx_r[i] * gradCx_y_r[i]- D_lamtau_xx_i[i] * gradCx_y_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xx_i[i] * gradCx_y_r[i]+ D_lamtau_xx_r[i] * gradCx_y_i[i]);

                    gamx_Y_r[i] += (D_sig_xy_y_r[i] * rhoBy_r[i]- D_sig_xy_y_i[i] * rhoBy_i[i]);

                    gamx_Y_i[i] +=  (D_sig_xy_y_i[i] * rhoBy_r[i]+ D_sig_xy_y_r[i] * rhoBy_i[i]);

                    gamx_Y_r[i] += (D_sig_xy_r[i] * gradBy_y_r[i]- D_sig_xy_i[i] * gradBy_y_i[i]);

                    gamx_Y_i[i] +=  (D_sig_xy_i[i] * gradBy_y_r[i]+ D_sig_xy_r[i] * gradBy_y_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xy_y_r[i] * rhoCy_r[i]- D_lamtau_xy_y_i[i] * rhoCy_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xy_y_i[i] * rhoCy_r[i]+ D_lamtau_xy_y_r[i] * rhoCy_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xy_r[i] * gradCy_y_r[i]- D_lamtau_xy_i[i] * gradCy_y_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xy_i[i] * gradCy_y_r[i]+ D_lamtau_xy_r[i] * gradCy_y_i[i]);

                    gamx_Y_r[i] += (D_sig_xz_y_r[i] * rhoBz_r[i]- D_sig_xz_y_i[i] * rhoBz_i[i]);

                    gamx_Y_i[i] +=  (D_sig_xz_y_i[i] * rhoBz_r[i]+ D_sig_xz_y_r[i] * rhoBz_i[i]);

                    gamx_Y_r[i] += (D_sig_xz_r[i] * gradBz_y_r[i]- D_sig_xz_i[i] * gradBz_y_i[i]);

                    gamx_Y_i[i] +=  (D_sig_xz_i[i] * gradBz_y_r[i]+ D_sig_xz_r[i] * gradBz_y_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xz_y_r[i] * rhoCz_r[i]- D_lamtau_xz_y_i[i] * rhoCz_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xz_y_i[i] * rhoCz_r[i]+ D_lamtau_xz_y_r[i] * rhoCz_i[i]);

                    gamx_Y_r[i] += (D_lamtau_xz_r[i] * gradCz_y_r[i]- D_lamtau_xz_i[i] * gradCz_y_i[i]);

                    gamx_Y_i[i] +=  (D_lamtau_xz_i[i] * gradCz_y_r[i]+ D_lamtau_xz_r[i] * gradCz_y_i[i]);

                    gamx_Z_r[i] = (D_sig_xx_z_r[i] * rhoBx_r[i]- D_sig_xx_z_i[i] * rhoBx_i[i]);

                    gamx_Z_i[i] =  (D_sig_xx_z_i[i] * rhoBx_r[i]+ D_sig_xx_z_r[i] * rhoBx_i[i]);

                    gamx_Z_r[i] += (D_sig_xx_r[i] * gradBx_z_r[i]- D_sig_xx_i[i] * gradBx_z_i[i]);

                    gamx_Z_i[i] +=  (D_sig_xx_i[i] * gradBx_z_r[i]+ D_sig_xx_r[i] * gradBx_z_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xx_z_r[i] * rhoCx_r[i]- D_lamtau_xx_z_i[i] * rhoCx_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xx_z_i[i] * rhoCx_r[i]+ D_lamtau_xx_z_r[i] * rhoCx_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xx_r[i] * gradCx_z_r[i]- D_lamtau_xx_i[i] * gradCx_z_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xx_i[i] * gradCx_z_r[i]+ D_lamtau_xx_r[i] * gradCx_z_i[i]);

                    gamx_Z_r[i] += (D_sig_xy_z_r[i] * rhoBy_r[i]- D_sig_xy_z_i[i] * rhoBy_i[i]);

                    gamx_Z_i[i] +=  (D_sig_xy_z_i[i] * rhoBy_r[i]+ D_sig_xy_z_r[i] * rhoBy_i[i]);

                    gamx_Z_r[i] += (D_sig_xy_r[i] * gradBy_z_r[i]- D_sig_xy_i[i] * gradBy_z_i[i]);

                    gamx_Z_i[i] +=  (D_sig_xy_i[i] * gradBy_z_r[i]+ D_sig_xy_r[i] * gradBy_z_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xy_z_r[i] * rhoCy_r[i]- D_lamtau_xy_z_i[i] * rhoCy_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xy_z_i[i] * rhoCy_r[i]+ D_lamtau_xy_z_r[i] * rhoCy_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xy_r[i] * gradCy_z_r[i]- D_lamtau_xy_i[i] * gradCy_z_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xy_i[i] * gradCy_z_r[i]+ D_lamtau_xy_r[i] * gradCy_z_i[i]);

                    gamx_Z_r[i] += (D_sig_xz_z_r[i] * rhoBz_r[i]- D_sig_xz_z_i[i] * rhoBz_i[i]);

                    gamx_Z_i[i] +=  (D_sig_xz_z_i[i] * rhoBz_r[i]+ D_sig_xz_z_r[i] * rhoBz_i[i]);

                    gamx_Z_r[i] += (D_sig_xz_r[i] * gradBz_z_r[i]- D_sig_xz_i[i] * gradBz_z_i[i]);

                    gamx_Z_i[i] +=  (D_sig_xz_i[i] * gradBz_z_r[i]+ D_sig_xz_r[i] * gradBz_z_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xz_z_r[i] * rhoCz_r[i]- D_lamtau_xz_z_i[i] * rhoCz_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xz_z_i[i] * rhoCz_r[i]+ D_lamtau_xz_z_r[i] * rhoCz_i[i]);

                    gamx_Z_r[i] += (D_lamtau_xz_r[i] * gradCz_z_r[i]- D_lamtau_xz_i[i] * gradCz_z_i[i]);

                    gamx_Z_i[i] +=  (D_lamtau_xz_i[i] * gradCz_z_r[i]+ D_lamtau_xz_r[i] * gradCz_z_i[i]);

                    gamy_X_r[i] = (D_sig_yx_x_r[i] * rhoBx_r[i]- D_sig_yx_x_i[i] * rhoBx_i[i]);

                    gamy_X_i[i] =  (D_sig_yx_x_i[i] * rhoBx_r[i]+ D_sig_yx_x_r[i] * rhoBx_i[i]);

                    gamy_X_r[i] += (D_sig_yx_r[i] * gradBx_x_r[i]- D_sig_yx_i[i] * gradBx_x_i[i]);

                    gamy_X_i[i] +=  (D_sig_yx_i[i] * gradBx_x_r[i]+ D_sig_yx_r[i] * gradBx_x_i[i]);

                    gamy_X_r[i] += (D_lamtau_yx_x_r[i] * rhoCx_r[i]- D_lamtau_yx_x_i[i] * rhoCx_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yx_x_i[i] * rhoCx_r[i]+ D_lamtau_yx_x_r[i] * rhoCx_i[i]);

                    gamy_X_r[i] += (D_lamtau_yx_r[i] * gradCx_x_r[i]- D_lamtau_yx_i[i] * gradCx_x_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yx_i[i] * gradCx_x_r[i]+ D_lamtau_yx_r[i] * gradCx_x_i[i]);

                    gamy_X_r[i] += (D_sig_yy_x_r[i] * rhoBy_r[i]- D_sig_yy_x_i[i] * rhoBy_i[i]);

                    gamy_X_i[i] +=  (D_sig_yy_x_i[i] * rhoBy_r[i]+ D_sig_yy_x_r[i] * rhoBy_i[i]);

                    gamy_X_r[i] += (D_sig_yy_r[i] * gradBy_x_r[i]- D_sig_yy_i[i] * gradBy_x_i[i]);

                    gamy_X_i[i] +=  (D_sig_yy_i[i] * gradBy_x_r[i]+ D_sig_yy_r[i] * gradBy_x_i[i]);

                    gamy_X_r[i] += (D_lamtau_yy_x_r[i] * rhoCy_r[i]- D_lamtau_yy_x_i[i] * rhoCy_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yy_x_i[i] * rhoCy_r[i]+ D_lamtau_yy_x_r[i] * rhoCy_i[i]);

                    gamy_X_r[i] += (D_lamtau_yy_r[i] * gradCy_x_r[i]- D_lamtau_yy_i[i] * gradCy_x_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yy_i[i] * gradCy_x_r[i]+ D_lamtau_yy_r[i] * gradCy_x_i[i]);

                    gamy_X_r[i] += (D_sig_yz_x_r[i] * rhoBz_r[i]- D_sig_yz_x_i[i] * rhoBz_i[i]);

                    gamy_X_i[i] +=  (D_sig_yz_x_i[i] * rhoBz_r[i]+ D_sig_yz_x_r[i] * rhoBz_i[i]);

                    gamy_X_r[i] += (D_sig_yz_r[i] * gradBz_x_r[i]- D_sig_yz_i[i] * gradBz_x_i[i]);

                    gamy_X_i[i] +=  (D_sig_yz_i[i] * gradBz_x_r[i]+ D_sig_yz_r[i] * gradBz_x_i[i]);

                    gamy_X_r[i] += (D_lamtau_yz_x_r[i] * rhoCz_r[i]- D_lamtau_yz_x_i[i] * rhoCz_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yz_x_i[i] * rhoCz_r[i]+ D_lamtau_yz_x_r[i] * rhoCz_i[i]);

                    gamy_X_r[i] += (D_lamtau_yz_r[i] * gradCz_x_r[i]- D_lamtau_yz_i[i] * gradCz_x_i[i]);

                    gamy_X_i[i] +=  (D_lamtau_yz_i[i] * gradCz_x_r[i]+ D_lamtau_yz_r[i] * gradCz_x_i[i]);

                    gamy_Y_r[i] = (D_sig_yx_y_r[i] * rhoBx_r[i]- D_sig_yx_y_i[i] * rhoBx_i[i]);

                    gamy_Y_i[i] =  (D_sig_yx_y_i[i] * rhoBx_r[i]+ D_sig_yx_y_r[i] * rhoBx_i[i]);

                    gamy_Y_r[i] += (D_sig_yx_r[i] * gradBx_y_r[i]- D_sig_yx_i[i] * gradBx_y_i[i]);

                    gamy_Y_i[i] +=  (D_sig_yx_i[i] * gradBx_y_r[i]+ D_sig_yx_r[i] * gradBx_y_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yx_y_r[i] * rhoCx_r[i]- D_lamtau_yx_y_i[i] * rhoCx_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yx_y_i[i] * rhoCx_r[i]+ D_lamtau_yx_y_r[i] * rhoCx_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yx_r[i] * gradCx_y_r[i]- D_lamtau_yx_i[i] * gradCx_y_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yx_i[i] * gradCx_y_r[i]+ D_lamtau_yx_r[i] * gradCx_y_i[i]);

                    gamy_Y_r[i] += (D_sig_yy_y_r[i] * rhoBy_r[i]- D_sig_yy_y_i[i] * rhoBy_i[i]);

                    gamy_Y_i[i] +=  (D_sig_yy_y_i[i] * rhoBy_r[i]+ D_sig_yy_y_r[i] * rhoBy_i[i]);

                    gamy_Y_r[i] += (D_sig_yy_r[i] * gradBy_y_r[i]- D_sig_yy_i[i] * gradBy_y_i[i]);

                    gamy_Y_i[i] +=  (D_sig_yy_i[i] * gradBy_y_r[i]+ D_sig_yy_r[i] * gradBy_y_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yy_y_r[i] * rhoCy_r[i]- D_lamtau_yy_y_i[i] * rhoCy_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yy_y_i[i] * rhoCy_r[i]+ D_lamtau_yy_y_r[i] * rhoCy_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yy_r[i] * gradCy_y_r[i]- D_lamtau_yy_i[i] * gradCy_y_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yy_i[i] * gradCy_y_r[i]+ D_lamtau_yy_r[i] * gradCy_y_i[i]);

                    gamy_Y_r[i] += (D_sig_yz_y_r[i] * rhoBz_r[i]- D_sig_yz_y_i[i] * rhoBz_i[i]);

                    gamy_Y_i[i] +=  (D_sig_yz_y_i[i] * rhoBz_r[i]+ D_sig_yz_y_r[i] * rhoBz_i[i]);

                    gamy_Y_r[i] += (D_sig_yz_r[i] * gradBz_y_r[i]- D_sig_yz_i[i] * gradBz_y_i[i]);

                    gamy_Y_i[i] +=  (D_sig_yz_i[i] * gradBz_y_r[i]+ D_sig_yz_r[i] * gradBz_y_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yz_y_r[i] * rhoCz_r[i]- D_lamtau_yz_y_i[i] * rhoCz_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yz_y_i[i] * rhoCz_r[i]+ D_lamtau_yz_y_r[i] * rhoCz_i[i]);

                    gamy_Y_r[i] += (D_lamtau_yz_r[i] * gradCz_y_r[i]- D_lamtau_yz_i[i] * gradCz_y_i[i]);

                    gamy_Y_i[i] +=  (D_lamtau_yz_i[i] * gradCz_y_r[i]+ D_lamtau_yz_r[i] * gradCz_y_i[i]);

                    gamy_Z_r[i] = (D_sig_yx_z_r[i] * rhoBx_r[i]- D_sig_yx_z_i[i] * rhoBx_i[i]);

                    gamy_Z_i[i] =  (D_sig_yx_z_i[i] * rhoBx_r[i]+ D_sig_yx_z_r[i] * rhoBx_i[i]);

                    gamy_Z_r[i] += (D_sig_yx_r[i] * gradBx_z_r[i]- D_sig_yx_i[i] * gradBx_z_i[i]);

                    gamy_Z_i[i] +=  (D_sig_yx_i[i] * gradBx_z_r[i]+ D_sig_yx_r[i] * gradBx_z_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yx_z_r[i] * rhoCx_r[i]- D_lamtau_yx_z_i[i] * rhoCx_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yx_z_i[i] * rhoCx_r[i]+ D_lamtau_yx_z_r[i] * rhoCx_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yx_r[i] * gradCx_z_r[i]- D_lamtau_yx_i[i] * gradCx_z_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yx_i[i] * gradCx_z_r[i]+ D_lamtau_yx_r[i] * gradCx_z_i[i]);

                    gamy_Z_r[i] += (D_sig_yy_z_r[i] * rhoBy_r[i]- D_sig_yy_z_i[i] * rhoBy_i[i]);

                    gamy_Z_i[i] +=  (D_sig_yy_z_i[i] * rhoBy_r[i]+ D_sig_yy_z_r[i] * rhoBy_i[i]);

                    gamy_Z_r[i] += (D_sig_yy_r[i] * gradBy_z_r[i]- D_sig_yy_i[i] * gradBy_z_i[i]);

                    gamy_Z_i[i] +=  (D_sig_yy_i[i] * gradBy_z_r[i]+ D_sig_yy_r[i] * gradBy_z_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yy_z_r[i] * rhoCy_r[i]- D_lamtau_yy_z_i[i] * rhoCy_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yy_z_i[i] * rhoCy_r[i]+ D_lamtau_yy_z_r[i] * rhoCy_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yy_r[i] * gradCy_z_r[i]- D_lamtau_yy_i[i] * gradCy_z_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yy_i[i] * gradCy_z_r[i]+ D_lamtau_yy_r[i] * gradCy_z_i[i]);

                    gamy_Z_r[i] += (D_sig_yz_z_r[i] * rhoBz_r[i]- D_sig_yz_z_i[i] * rhoBz_i[i]);

                    gamy_Z_i[i] +=  (D_sig_yz_z_i[i] * rhoBz_r[i]+ D_sig_yz_z_r[i] * rhoBz_i[i]);

                    gamy_Z_r[i] += (D_sig_yz_r[i] * gradBz_z_r[i]- D_sig_yz_i[i] * gradBz_z_i[i]);

                    gamy_Z_i[i] +=  (D_sig_yz_i[i] * gradBz_z_r[i]+ D_sig_yz_r[i] * gradBz_z_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yz_z_r[i] * rhoCz_r[i]- D_lamtau_yz_z_i[i] * rhoCz_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yz_z_i[i] * rhoCz_r[i]+ D_lamtau_yz_z_r[i] * rhoCz_i[i]);

                    gamy_Z_r[i] += (D_lamtau_yz_r[i] * gradCz_z_r[i]- D_lamtau_yz_i[i] * gradCz_z_i[i]);

                    gamy_Z_i[i] +=  (D_lamtau_yz_i[i] * gradCz_z_r[i]+ D_lamtau_yz_r[i] * gradCz_z_i[i]);

                    gamz_X_r[i] = (D_sig_zx_x_r[i] * rhoBx_r[i]- D_sig_zx_x_i[i] * rhoBx_i[i]);

                    gamz_X_i[i] =  (D_sig_zx_x_i[i] * rhoBx_r[i]+ D_sig_zx_x_r[i] * rhoBx_i[i]);

                    gamz_X_r[i] += (D_sig_zx_r[i] * gradBx_x_r[i]- D_sig_zx_i[i] * gradBx_x_i[i]);

                    gamz_X_i[i] +=  (D_sig_zx_i[i] * gradBx_x_r[i]+ D_sig_zx_r[i] * gradBx_x_i[i]);

                    gamz_X_r[i] += (D_lamtau_zx_x_r[i] * rhoCx_r[i]- D_lamtau_zx_x_i[i] * rhoCx_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zx_x_i[i] * rhoCx_r[i]+ D_lamtau_zx_x_r[i] * rhoCx_i[i]);

                    gamz_X_r[i] += (D_lamtau_zx_r[i] * gradCx_x_r[i]- D_lamtau_zx_i[i] * gradCx_x_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zx_i[i] * gradCx_x_r[i]+ D_lamtau_zx_r[i] * gradCx_x_i[i]);

                    gamz_X_r[i] += (D_sig_zy_x_r[i] * rhoBy_r[i]- D_sig_zy_x_i[i] * rhoBy_i[i]);

                    gamz_X_i[i] +=  (D_sig_zy_x_i[i] * rhoBy_r[i]+ D_sig_zy_x_r[i] * rhoBy_i[i]);

                    gamz_X_r[i] += (D_sig_zy_r[i] * gradBy_x_r[i]- D_sig_zy_i[i] * gradBy_x_i[i]);

                    gamz_X_i[i] +=  (D_sig_zy_i[i] * gradBy_x_r[i]+ D_sig_zy_r[i] * gradBy_x_i[i]);

                    gamz_X_r[i] += (D_lamtau_zy_x_r[i] * rhoCy_r[i]- D_lamtau_zy_x_i[i] * rhoCy_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zy_x_i[i] * rhoCy_r[i]+ D_lamtau_zy_x_r[i] * rhoCy_i[i]);

                    gamz_X_r[i] += (D_lamtau_zy_r[i] * gradCy_x_r[i]- D_lamtau_zy_i[i] * gradCy_x_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zy_i[i] * gradCy_x_r[i]+ D_lamtau_zy_r[i] * gradCy_x_i[i]);

                    gamz_X_r[i] += (D_sig_zz_x_r[i] * rhoBz_r[i]- D_sig_zz_x_i[i] * rhoBz_i[i]);

                    gamz_X_i[i] +=  (D_sig_zz_x_i[i] * rhoBz_r[i]+ D_sig_zz_x_r[i] * rhoBz_i[i]);

                    gamz_X_r[i] += (D_sig_zz_r[i] * gradBz_x_r[i]- D_sig_zz_i[i] * gradBz_x_i[i]);

                    gamz_X_i[i] +=  (D_sig_zz_i[i] * gradBz_x_r[i]+ D_sig_zz_r[i] * gradBz_x_i[i]);

                    gamz_X_r[i] += (D_lamtau_zz_x_r[i] * rhoCz_r[i]- D_lamtau_zz_x_i[i] * rhoCz_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zz_x_i[i] * rhoCz_r[i]+ D_lamtau_zz_x_r[i] * rhoCz_i[i]);

                    gamz_X_r[i] += (D_lamtau_zz_r[i] * gradCz_x_r[i]- D_lamtau_zz_i[i] * gradCz_x_i[i]);

                    gamz_X_i[i] +=  (D_lamtau_zz_i[i] * gradCz_x_r[i]+ D_lamtau_zz_r[i] * gradCz_x_i[i]);

                    gamz_Y_r[i] = (D_sig_zx_y_r[i] * rhoBx_r[i]- D_sig_zx_y_i[i] * rhoBx_i[i]);

                    gamz_Y_i[i] =  (D_sig_zx_y_i[i] * rhoBx_r[i]+ D_sig_zx_y_r[i] * rhoBx_i[i]);

                    gamz_Y_r[i] += (D_sig_zx_r[i] * gradBx_y_r[i]- D_sig_zx_i[i] * gradBx_y_i[i]);

                    gamz_Y_i[i] +=  (D_sig_zx_i[i] * gradBx_y_r[i]+ D_sig_zx_r[i] * gradBx_y_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zx_y_r[i] * rhoCx_r[i]- D_lamtau_zx_y_i[i] * rhoCx_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zx_y_i[i] * rhoCx_r[i]+ D_lamtau_zx_y_r[i] * rhoCx_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zx_r[i] * gradCx_y_r[i]- D_lamtau_zx_i[i] * gradCx_y_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zx_i[i] * gradCx_y_r[i]+ D_lamtau_zx_r[i] * gradCx_y_i[i]);

                    gamz_Y_r[i] += (D_sig_zy_y_r[i] * rhoBy_r[i]- D_sig_zy_y_i[i] * rhoBy_i[i]);

                    gamz_Y_i[i] +=  (D_sig_zy_y_i[i] * rhoBy_r[i]+ D_sig_zy_y_r[i] * rhoBy_i[i]);

                    gamz_Y_r[i] += (D_sig_zy_r[i] * gradBy_y_r[i]- D_sig_zy_i[i] * gradBy_y_i[i]);

                    gamz_Y_i[i] +=  (D_sig_zy_i[i] * gradBy_y_r[i]+ D_sig_zy_r[i] * gradBy_y_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zy_y_r[i] * rhoCy_r[i]- D_lamtau_zy_y_i[i] * rhoCy_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zy_y_i[i] * rhoCy_r[i]+ D_lamtau_zy_y_r[i] * rhoCy_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zy_r[i] * gradCy_y_r[i]- D_lamtau_zy_i[i] * gradCy_y_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zy_i[i] * gradCy_y_r[i]+ D_lamtau_zy_r[i] * gradCy_y_i[i]);

                    gamz_Y_r[i] += (D_sig_zz_y_r[i] * rhoBz_r[i]- D_sig_zz_y_i[i] * rhoBz_i[i]);

                    gamz_Y_i[i] +=  (D_sig_zz_y_i[i] * rhoBz_r[i]+ D_sig_zz_y_r[i] * rhoBz_i[i]);

                    gamz_Y_r[i] += (D_sig_zz_r[i] * gradBz_y_r[i]- D_sig_zz_i[i] * gradBz_y_i[i]);

                    gamz_Y_i[i] +=  (D_sig_zz_i[i] * gradBz_y_r[i]+ D_sig_zz_r[i] * gradBz_y_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zz_y_r[i] * rhoCz_r[i]- D_lamtau_zz_y_i[i] * rhoCz_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zz_y_i[i] * rhoCz_r[i]+ D_lamtau_zz_y_r[i] * rhoCz_i[i]);

                    gamz_Y_r[i] += (D_lamtau_zz_r[i] * gradCz_y_r[i]- D_lamtau_zz_i[i] * gradCz_y_i[i]);

                    gamz_Y_i[i] +=  (D_lamtau_zz_i[i] * gradCz_y_r[i]+ D_lamtau_zz_r[i] * gradCz_y_i[i]);

                    gamz_Z_r[i] = (D_sig_zx_z_r[i] * rhoBx_r[i]- D_sig_zx_z_i[i] * rhoBx_i[i]);

                    gamz_Z_i[i] =  (D_sig_zx_z_i[i] * rhoBx_r[i]+ D_sig_zx_z_r[i] * rhoBx_i[i]);

                    gamz_Z_r[i] += (D_sig_zx_r[i] * gradBx_z_r[i]- D_sig_zx_i[i] * gradBx_z_i[i]);

                    gamz_Z_i[i] +=  (D_sig_zx_i[i] * gradBx_z_r[i]+ D_sig_zx_r[i] * gradBx_z_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zx_z_r[i] * rhoCx_r[i]- D_lamtau_zx_z_i[i] * rhoCx_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zx_z_i[i] * rhoCx_r[i]+ D_lamtau_zx_z_r[i] * rhoCx_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zx_r[i] * gradCx_z_r[i]- D_lamtau_zx_i[i] * gradCx_z_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zx_i[i] * gradCx_z_r[i]+ D_lamtau_zx_r[i] * gradCx_z_i[i]);

                    gamz_Z_r[i] += (D_sig_zy_z_r[i] * rhoBy_r[i]- D_sig_zy_z_i[i] * rhoBy_i[i]);

                    gamz_Z_i[i] +=  (D_sig_zy_z_i[i] * rhoBy_r[i]+ D_sig_zy_z_r[i] * rhoBy_i[i]);

                    gamz_Z_r[i] += (D_sig_zy_r[i] * gradBy_z_r[i]- D_sig_zy_i[i] * gradBy_z_i[i]);

                    gamz_Z_i[i] +=  (D_sig_zy_i[i] * gradBy_z_r[i]+ D_sig_zy_r[i] * gradBy_z_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zy_z_r[i] * rhoCy_r[i]- D_lamtau_zy_z_i[i] * rhoCy_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zy_z_i[i] * rhoCy_r[i]+ D_lamtau_zy_z_r[i] * rhoCy_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zy_r[i] * gradCy_z_r[i]- D_lamtau_zy_i[i] * gradCy_z_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zy_i[i] * gradCy_z_r[i]+ D_lamtau_zy_r[i] * gradCy_z_i[i]);

                    gamz_Z_r[i] += (D_sig_zz_z_r[i] * rhoBz_r[i]- D_sig_zz_z_i[i] * rhoBz_i[i]);

                    gamz_Z_i[i] +=  (D_sig_zz_z_i[i] * rhoBz_r[i]+ D_sig_zz_z_r[i] * rhoBz_i[i]);

                    gamz_Z_r[i] += (D_sig_zz_r[i] * gradBz_z_r[i]- D_sig_zz_i[i] * gradBz_z_i[i]);

                    gamz_Z_i[i] +=  (D_sig_zz_i[i] * gradBz_z_r[i]+ D_sig_zz_r[i] * gradBz_z_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zz_z_r[i] * rhoCz_r[i]- D_lamtau_zz_z_i[i] * rhoCz_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zz_z_i[i] * rhoCz_r[i]+ D_lamtau_zz_z_r[i] * rhoCz_i[i]);

                    gamz_Z_r[i] += (D_lamtau_zz_r[i] * gradCz_z_r[i]- D_lamtau_zz_i[i] * gradCz_z_i[i]);

                    gamz_Z_i[i] +=  (D_lamtau_zz_i[i] * gradCz_z_r[i]+ D_lamtau_zz_r[i] * gradCz_z_i[i]);

                    gamx_XX_r[i] = (D_sig_xx_x_r[i] * gradCx_x_r[i]- D_sig_xx_x_i[i] * gradCx_x_i[i]);

                    gamx_XX_i[i] =  (D_sig_xx_x_i[i] * gradCx_x_r[i]+ D_sig_xx_x_r[i] * gradCx_x_i[i]);

                    gamx_XX_r[i] += (D_lamtau_xx_x_r[i] * gradBx_x_r[i]- D_lamtau_xx_x_i[i] * gradBx_x_i[i]);

                    gamx_XX_i[i] +=  (D_lamtau_xx_x_i[i] * gradBx_x_r[i]+ D_lamtau_xx_x_r[i] * gradBx_x_i[i]);

                    gamx_XX_r[i] += (D_sig_xy_x_r[i] * gradCy_x_r[i]- D_sig_xy_x_i[i] * gradCy_x_i[i]);

                    gamx_XX_i[i] +=  (D_sig_xy_x_i[i] * gradCy_x_r[i]+ D_sig_xy_x_r[i] * gradCy_x_i[i]);

                    gamx_XX_r[i] += (D_lamtau_xy_x_r[i] * gradBy_x_r[i]- D_lamtau_xy_x_i[i] * gradBy_x_i[i]);

                    gamx_XX_i[i] +=  (D_lamtau_xy_x_i[i] * gradBy_x_r[i]+ D_lamtau_xy_x_r[i] * gradBy_x_i[i]);

                    gamx_XX_r[i] += (D_sig_xz_x_r[i] * gradCz_x_r[i]- D_sig_xz_x_i[i] * gradCz_x_i[i]);

                    gamx_XX_i[i] +=  (D_sig_xz_x_i[i] * gradCz_x_r[i]+ D_sig_xz_x_r[i] * gradCz_x_i[i]);

                    gamx_XX_r[i] += (D_lamtau_xz_x_r[i] * gradBz_x_r[i]- D_lamtau_xz_x_i[i] * gradBz_x_i[i]);

                    gamx_XX_i[i] +=  (D_lamtau_xz_x_i[i] * gradBz_x_r[i]+ D_lamtau_xz_x_r[i] * gradBz_x_i[i]);

                    gamx_XY_r[i] = (D_sig_xx_x_r[i] * gradCx_y_r[i]- D_sig_xx_x_i[i] * gradCx_y_i[i]);

                    gamx_XY_i[i] =  (D_sig_xx_x_i[i] * gradCx_y_r[i]+ D_sig_xx_x_r[i] * gradCx_y_i[i]);

                    gamx_XY_r[i] += (D_lamtau_xx_x_r[i] * gradBx_y_r[i]- D_lamtau_xx_x_i[i] * gradBx_y_i[i]);

                    gamx_XY_i[i] +=  (D_lamtau_xx_x_i[i] * gradBx_y_r[i]+ D_lamtau_xx_x_r[i] * gradBx_y_i[i]);

                    gamx_XY_r[i] += (D_sig_xy_x_r[i] * gradCy_y_r[i]- D_sig_xy_x_i[i] * gradCy_y_i[i]);

                    gamx_XY_i[i] +=  (D_sig_xy_x_i[i] * gradCy_y_r[i]+ D_sig_xy_x_r[i] * gradCy_y_i[i]);

                    gamx_XY_r[i] += (D_lamtau_xy_x_r[i] * gradBy_y_r[i]- D_lamtau_xy_x_i[i] * gradBy_y_i[i]);

                    gamx_XY_i[i] +=  (D_lamtau_xy_x_i[i] * gradBy_y_r[i]+ D_lamtau_xy_x_r[i] * gradBy_y_i[i]);

                    gamx_XY_r[i] += (D_sig_xz_x_r[i] * gradCz_y_r[i]- D_sig_xz_x_i[i] * gradCz_y_i[i]);

                    gamx_XY_i[i] +=  (D_sig_xz_x_i[i] * gradCz_y_r[i]+ D_sig_xz_x_r[i] * gradCz_y_i[i]);

                    gamx_XY_r[i] += (D_lamtau_xz_x_r[i] * gradBz_y_r[i]- D_lamtau_xz_x_i[i] * gradBz_y_i[i]);

                    gamx_XY_i[i] +=  (D_lamtau_xz_x_i[i] * gradBz_y_r[i]+ D_lamtau_xz_x_r[i] * gradBz_y_i[i]);

                    gamx_XZ_r[i] = (D_sig_xx_x_r[i] * gradCx_z_r[i]- D_sig_xx_x_i[i] * gradCx_z_i[i]);

                    gamx_XZ_i[i] =  (D_sig_xx_x_i[i] * gradCx_z_r[i]+ D_sig_xx_x_r[i] * gradCx_z_i[i]);

                    gamx_XZ_r[i] += (D_lamtau_xx_x_r[i] * gradBx_z_r[i]- D_lamtau_xx_x_i[i] * gradBx_z_i[i]);

                    gamx_XZ_i[i] +=  (D_lamtau_xx_x_i[i] * gradBx_z_r[i]+ D_lamtau_xx_x_r[i] * gradBx_z_i[i]);

                    gamx_XZ_r[i] += (D_sig_xy_x_r[i] * gradCy_z_r[i]- D_sig_xy_x_i[i] * gradCy_z_i[i]);

                    gamx_XZ_i[i] +=  (D_sig_xy_x_i[i] * gradCy_z_r[i]+ D_sig_xy_x_r[i] * gradCy_z_i[i]);

                    gamx_XZ_r[i] += (D_lamtau_xy_x_r[i] * gradBy_z_r[i]- D_lamtau_xy_x_i[i] * gradBy_z_i[i]);

                    gamx_XZ_i[i] +=  (D_lamtau_xy_x_i[i] * gradBy_z_r[i]+ D_lamtau_xy_x_r[i] * gradBy_z_i[i]);

                    gamx_XZ_r[i] += (D_sig_xz_x_r[i] * gradCz_z_r[i]- D_sig_xz_x_i[i] * gradCz_z_i[i]);

                    gamx_XZ_i[i] +=  (D_sig_xz_x_i[i] * gradCz_z_r[i]+ D_sig_xz_x_r[i] * gradCz_z_i[i]);

                    gamx_XZ_r[i] += (D_lamtau_xz_x_r[i] * gradBz_z_r[i]- D_lamtau_xz_x_i[i] * gradBz_z_i[i]);

                    gamx_XZ_i[i] +=  (D_lamtau_xz_x_i[i] * gradBz_z_r[i]+ D_lamtau_xz_x_r[i] * gradBz_z_i[i]);

                    gamx_YX_r[i] = (D_sig_xx_y_r[i] * gradCx_x_r[i]- D_sig_xx_y_i[i] * gradCx_x_i[i]);

                    gamx_YX_i[i] =  (D_sig_xx_y_i[i] * gradCx_x_r[i]+ D_sig_xx_y_r[i] * gradCx_x_i[i]);

                    gamx_YX_r[i] += (D_lamtau_xx_y_r[i] * gradBx_x_r[i]- D_lamtau_xx_y_i[i] * gradBx_x_i[i]);

                    gamx_YX_i[i] +=  (D_lamtau_xx_y_i[i] * gradBx_x_r[i]+ D_lamtau_xx_y_r[i] * gradBx_x_i[i]);

                    gamx_YX_r[i] += (D_sig_xy_y_r[i] * gradCy_x_r[i]- D_sig_xy_y_i[i] * gradCy_x_i[i]);

                    gamx_YX_i[i] +=  (D_sig_xy_y_i[i] * gradCy_x_r[i]+ D_sig_xy_y_r[i] * gradCy_x_i[i]);

                    gamx_YX_r[i] += (D_lamtau_xy_y_r[i] * gradBy_x_r[i]- D_lamtau_xy_y_i[i] * gradBy_x_i[i]);

                    gamx_YX_i[i] +=  (D_lamtau_xy_y_i[i] * gradBy_x_r[i]+ D_lamtau_xy_y_r[i] * gradBy_x_i[i]);

                    gamx_YX_r[i] += (D_sig_xz_y_r[i] * gradCz_x_r[i]- D_sig_xz_y_i[i] * gradCz_x_i[i]);

                    gamx_YX_i[i] +=  (D_sig_xz_y_i[i] * gradCz_x_r[i]+ D_sig_xz_y_r[i] * gradCz_x_i[i]);

                    gamx_YX_r[i] += (D_lamtau_xz_y_r[i] * gradBz_x_r[i]- D_lamtau_xz_y_i[i] * gradBz_x_i[i]);

                    gamx_YX_i[i] +=  (D_lamtau_xz_y_i[i] * gradBz_x_r[i]+ D_lamtau_xz_y_r[i] * gradBz_x_i[i]);

                    gamx_YY_r[i] = (D_sig_xx_y_r[i] * gradCx_y_r[i]- D_sig_xx_y_i[i] * gradCx_y_i[i]);

                    gamx_YY_i[i] =  (D_sig_xx_y_i[i] * gradCx_y_r[i]+ D_sig_xx_y_r[i] * gradCx_y_i[i]);

                    gamx_YY_r[i] += (D_lamtau_xx_y_r[i] * gradBx_y_r[i]- D_lamtau_xx_y_i[i] * gradBx_y_i[i]);

                    gamx_YY_i[i] +=  (D_lamtau_xx_y_i[i] * gradBx_y_r[i]+ D_lamtau_xx_y_r[i] * gradBx_y_i[i]);

                    gamx_YY_r[i] += (D_sig_xy_y_r[i] * gradCy_y_r[i]- D_sig_xy_y_i[i] * gradCy_y_i[i]);

                    gamx_YY_i[i] +=  (D_sig_xy_y_i[i] * gradCy_y_r[i]+ D_sig_xy_y_r[i] * gradCy_y_i[i]);

                    gamx_YY_r[i] += (D_lamtau_xy_y_r[i] * gradBy_y_r[i]- D_lamtau_xy_y_i[i] * gradBy_y_i[i]);

                    gamx_YY_i[i] +=  (D_lamtau_xy_y_i[i] * gradBy_y_r[i]+ D_lamtau_xy_y_r[i] * gradBy_y_i[i]);

                    gamx_YY_r[i] += (D_sig_xz_y_r[i] * gradCz_y_r[i]- D_sig_xz_y_i[i] * gradCz_y_i[i]);

                    gamx_YY_i[i] +=  (D_sig_xz_y_i[i] * gradCz_y_r[i]+ D_sig_xz_y_r[i] * gradCz_y_i[i]);

                    gamx_YY_r[i] += (D_lamtau_xz_y_r[i] * gradBz_y_r[i]- D_lamtau_xz_y_i[i] * gradBz_y_i[i]);

                    gamx_YY_i[i] +=  (D_lamtau_xz_y_i[i] * gradBz_y_r[i]+ D_lamtau_xz_y_r[i] * gradBz_y_i[i]);

                    gamx_YZ_r[i] = (D_sig_xx_y_r[i] * gradCx_z_r[i]- D_sig_xx_y_i[i] * gradCx_z_i[i]);

                    gamx_YZ_i[i] =  (D_sig_xx_y_i[i] * gradCx_z_r[i]+ D_sig_xx_y_r[i] * gradCx_z_i[i]);

                    gamx_YZ_r[i] += (D_lamtau_xx_y_r[i] * gradBx_z_r[i]- D_lamtau_xx_y_i[i] * gradBx_z_i[i]);

                    gamx_YZ_i[i] +=  (D_lamtau_xx_y_i[i] * gradBx_z_r[i]+ D_lamtau_xx_y_r[i] * gradBx_z_i[i]);

                    gamx_YZ_r[i] += (D_sig_xy_y_r[i] * gradCy_z_r[i]- D_sig_xy_y_i[i] * gradCy_z_i[i]);

                    gamx_YZ_i[i] +=  (D_sig_xy_y_i[i] * gradCy_z_r[i]+ D_sig_xy_y_r[i] * gradCy_z_i[i]);

                    gamx_YZ_r[i] += (D_lamtau_xy_y_r[i] * gradBy_z_r[i]- D_lamtau_xy_y_i[i] * gradBy_z_i[i]);

                    gamx_YZ_i[i] +=  (D_lamtau_xy_y_i[i] * gradBy_z_r[i]+ D_lamtau_xy_y_r[i] * gradBy_z_i[i]);

                    gamx_YZ_r[i] += (D_sig_xz_y_r[i] * gradCz_z_r[i]- D_sig_xz_y_i[i] * gradCz_z_i[i]);

                    gamx_YZ_i[i] +=  (D_sig_xz_y_i[i] * gradCz_z_r[i]+ D_sig_xz_y_r[i] * gradCz_z_i[i]);

                    gamx_YZ_r[i] += (D_lamtau_xz_y_r[i] * gradBz_z_r[i]- D_lamtau_xz_y_i[i] * gradBz_z_i[i]);

                    gamx_YZ_i[i] +=  (D_lamtau_xz_y_i[i] * gradBz_z_r[i]+ D_lamtau_xz_y_r[i] * gradBz_z_i[i]);

                    gamx_ZX_r[i] = (D_sig_xx_z_r[i] * gradCx_x_r[i]- D_sig_xx_z_i[i] * gradCx_x_i[i]);

                    gamx_ZX_i[i] =  (D_sig_xx_z_i[i] * gradCx_x_r[i]+ D_sig_xx_z_r[i] * gradCx_x_i[i]);

                    gamx_ZX_r[i] += (D_lamtau_xx_z_r[i] * gradBx_x_r[i]- D_lamtau_xx_z_i[i] * gradBx_x_i[i]);

                    gamx_ZX_i[i] +=  (D_lamtau_xx_z_i[i] * gradBx_x_r[i]+ D_lamtau_xx_z_r[i] * gradBx_x_i[i]);

                    gamx_ZX_r[i] += (D_sig_xy_z_r[i] * gradCy_x_r[i]- D_sig_xy_z_i[i] * gradCy_x_i[i]);

                    gamx_ZX_i[i] +=  (D_sig_xy_z_i[i] * gradCy_x_r[i]+ D_sig_xy_z_r[i] * gradCy_x_i[i]);

                    gamx_ZX_r[i] += (D_lamtau_xy_z_r[i] * gradBy_x_r[i]- D_lamtau_xy_z_i[i] * gradBy_x_i[i]);

                    gamx_ZX_i[i] +=  (D_lamtau_xy_z_i[i] * gradBy_x_r[i]+ D_lamtau_xy_z_r[i] * gradBy_x_i[i]);

                    gamx_ZX_r[i] += (D_sig_xz_z_r[i] * gradCz_x_r[i]- D_sig_xz_z_i[i] * gradCz_x_i[i]);

                    gamx_ZX_i[i] +=  (D_sig_xz_z_i[i] * gradCz_x_r[i]+ D_sig_xz_z_r[i] * gradCz_x_i[i]);

                    gamx_ZX_r[i] += (D_lamtau_xz_z_r[i] * gradBz_x_r[i]- D_lamtau_xz_z_i[i] * gradBz_x_i[i]);

                    gamx_ZX_i[i] +=  (D_lamtau_xz_z_i[i] * gradBz_x_r[i]+ D_lamtau_xz_z_r[i] * gradBz_x_i[i]);

                    gamx_ZY_r[i] = (D_sig_xx_z_r[i] * gradCx_y_r[i]- D_sig_xx_z_i[i] * gradCx_y_i[i]);

                    gamx_ZY_i[i] =  (D_sig_xx_z_i[i] * gradCx_y_r[i]+ D_sig_xx_z_r[i] * gradCx_y_i[i]);

                    gamx_ZY_r[i] += (D_lamtau_xx_z_r[i] * gradBx_y_r[i]- D_lamtau_xx_z_i[i] * gradBx_y_i[i]);

                    gamx_ZY_i[i] +=  (D_lamtau_xx_z_i[i] * gradBx_y_r[i]+ D_lamtau_xx_z_r[i] * gradBx_y_i[i]);

                    gamx_ZY_r[i] += (D_sig_xy_z_r[i] * gradCy_y_r[i]- D_sig_xy_z_i[i] * gradCy_y_i[i]);

                    gamx_ZY_i[i] +=  (D_sig_xy_z_i[i] * gradCy_y_r[i]+ D_sig_xy_z_r[i] * gradCy_y_i[i]);

                    gamx_ZY_r[i] += (D_lamtau_xy_z_r[i] * gradBy_y_r[i]- D_lamtau_xy_z_i[i] * gradBy_y_i[i]);

                    gamx_ZY_i[i] +=  (D_lamtau_xy_z_i[i] * gradBy_y_r[i]+ D_lamtau_xy_z_r[i] * gradBy_y_i[i]);

                    gamx_ZY_r[i] += (D_sig_xz_z_r[i] * gradCz_y_r[i]- D_sig_xz_z_i[i] * gradCz_y_i[i]);

                    gamx_ZY_i[i] +=  (D_sig_xz_z_i[i] * gradCz_y_r[i]+ D_sig_xz_z_r[i] * gradCz_y_i[i]);

                    gamx_ZY_r[i] += (D_lamtau_xz_z_r[i] * gradBz_y_r[i]- D_lamtau_xz_z_i[i] * gradBz_y_i[i]);

                    gamx_ZY_i[i] +=  (D_lamtau_xz_z_i[i] * gradBz_y_r[i]+ D_lamtau_xz_z_r[i] * gradBz_y_i[i]);

                    gamx_ZZ_r[i] = (D_sig_xx_z_r[i] * gradCx_z_r[i]- D_sig_xx_z_i[i] * gradCx_z_i[i]);

                    gamx_ZZ_i[i] =  (D_sig_xx_z_i[i] * gradCx_z_r[i]+ D_sig_xx_z_r[i] * gradCx_z_i[i]);

                    gamx_ZZ_r[i] += (D_lamtau_xx_z_r[i] * gradBx_z_r[i]- D_lamtau_xx_z_i[i] * gradBx_z_i[i]);

                    gamx_ZZ_i[i] +=  (D_lamtau_xx_z_i[i] * gradBx_z_r[i]+ D_lamtau_xx_z_r[i] * gradBx_z_i[i]);

                    gamx_ZZ_r[i] += (D_sig_xy_z_r[i] * gradCy_z_r[i]- D_sig_xy_z_i[i] * gradCy_z_i[i]);

                    gamx_ZZ_i[i] +=  (D_sig_xy_z_i[i] * gradCy_z_r[i]+ D_sig_xy_z_r[i] * gradCy_z_i[i]);

                    gamx_ZZ_r[i] += (D_lamtau_xy_z_r[i] * gradBy_z_r[i]- D_lamtau_xy_z_i[i] * gradBy_z_i[i]);

                    gamx_ZZ_i[i] +=  (D_lamtau_xy_z_i[i] * gradBy_z_r[i]+ D_lamtau_xy_z_r[i] * gradBy_z_i[i]);

                    gamx_ZZ_r[i] += (D_sig_xz_z_r[i] * gradCz_z_r[i]- D_sig_xz_z_i[i] * gradCz_z_i[i]);

                    gamx_ZZ_i[i] +=  (D_sig_xz_z_i[i] * gradCz_z_r[i]+ D_sig_xz_z_r[i] * gradCz_z_i[i]);

                    gamx_ZZ_r[i] += (D_lamtau_xz_z_r[i] * gradBz_z_r[i]- D_lamtau_xz_z_i[i] * gradBz_z_i[i]);

                    gamx_ZZ_i[i] +=  (D_lamtau_xz_z_i[i] * gradBz_z_r[i]+ D_lamtau_xz_z_r[i] * gradBz_z_i[i]);

                    gamy_XX_r[i] = (D_sig_yx_x_r[i] * gradCx_x_r[i]- D_sig_yx_x_i[i] * gradCx_x_i[i]);

                    gamy_XX_i[i] =  (D_sig_yx_x_i[i] * gradCx_x_r[i]+ D_sig_yx_x_r[i] * gradCx_x_i[i]);

                    gamy_XX_r[i] += (D_lamtau_yx_x_r[i] * gradBx_x_r[i]- D_lamtau_yx_x_i[i] * gradBx_x_i[i]);

                    gamy_XX_i[i] +=  (D_lamtau_yx_x_i[i] * gradBx_x_r[i]+ D_lamtau_yx_x_r[i] * gradBx_x_i[i]);

                    gamy_XX_r[i] += (D_sig_yy_x_r[i] * gradCy_x_r[i]- D_sig_yy_x_i[i] * gradCy_x_i[i]);

                    gamy_XX_i[i] +=  (D_sig_yy_x_i[i] * gradCy_x_r[i]+ D_sig_yy_x_r[i] * gradCy_x_i[i]);

                    gamy_XX_r[i] += (D_lamtau_yy_x_r[i] * gradBy_x_r[i]- D_lamtau_yy_x_i[i] * gradBy_x_i[i]);

                    gamy_XX_i[i] +=  (D_lamtau_yy_x_i[i] * gradBy_x_r[i]+ D_lamtau_yy_x_r[i] * gradBy_x_i[i]);

                    gamy_XX_r[i] += (D_sig_yz_x_r[i] * gradCz_x_r[i]- D_sig_yz_x_i[i] * gradCz_x_i[i]);

                    gamy_XX_i[i] +=  (D_sig_yz_x_i[i] * gradCz_x_r[i]+ D_sig_yz_x_r[i] * gradCz_x_i[i]);

                    gamy_XX_r[i] += (D_lamtau_yz_x_r[i] * gradBz_x_r[i]- D_lamtau_yz_x_i[i] * gradBz_x_i[i]);

                    gamy_XX_i[i] +=  (D_lamtau_yz_x_i[i] * gradBz_x_r[i]+ D_lamtau_yz_x_r[i] * gradBz_x_i[i]);

                    gamy_XY_r[i] = (D_sig_yx_x_r[i] * gradCx_y_r[i]- D_sig_yx_x_i[i] * gradCx_y_i[i]);

                    gamy_XY_i[i] =  (D_sig_yx_x_i[i] * gradCx_y_r[i]+ D_sig_yx_x_r[i] * gradCx_y_i[i]);

                    gamy_XY_r[i] += (D_lamtau_yx_x_r[i] * gradBx_y_r[i]- D_lamtau_yx_x_i[i] * gradBx_y_i[i]);

                    gamy_XY_i[i] +=  (D_lamtau_yx_x_i[i] * gradBx_y_r[i]+ D_lamtau_yx_x_r[i] * gradBx_y_i[i]);

                    gamy_XY_r[i] += (D_sig_yy_x_r[i] * gradCy_y_r[i]- D_sig_yy_x_i[i] * gradCy_y_i[i]);

                    gamy_XY_i[i] +=  (D_sig_yy_x_i[i] * gradCy_y_r[i]+ D_sig_yy_x_r[i] * gradCy_y_i[i]);

                    gamy_XY_r[i] += (D_lamtau_yy_x_r[i] * gradBy_y_r[i]- D_lamtau_yy_x_i[i] * gradBy_y_i[i]);

                    gamy_XY_i[i] +=  (D_lamtau_yy_x_i[i] * gradBy_y_r[i]+ D_lamtau_yy_x_r[i] * gradBy_y_i[i]);

                    gamy_XY_r[i] += (D_sig_yz_x_r[i] * gradCz_y_r[i]- D_sig_yz_x_i[i] * gradCz_y_i[i]);

                    gamy_XY_i[i] +=  (D_sig_yz_x_i[i] * gradCz_y_r[i]+ D_sig_yz_x_r[i] * gradCz_y_i[i]);

                    gamy_XY_r[i] += (D_lamtau_yz_x_r[i] * gradBz_y_r[i]- D_lamtau_yz_x_i[i] * gradBz_y_i[i]);

                    gamy_XY_i[i] +=  (D_lamtau_yz_x_i[i] * gradBz_y_r[i]+ D_lamtau_yz_x_r[i] * gradBz_y_i[i]);

                    gamy_XZ_r[i] = (D_sig_yx_x_r[i] * gradCx_z_r[i]- D_sig_yx_x_i[i] * gradCx_z_i[i]);

                    gamy_XZ_i[i] =  (D_sig_yx_x_i[i] * gradCx_z_r[i]+ D_sig_yx_x_r[i] * gradCx_z_i[i]);

                    gamy_XZ_r[i] += (D_lamtau_yx_x_r[i] * gradBx_z_r[i]- D_lamtau_yx_x_i[i] * gradBx_z_i[i]);

                    gamy_XZ_i[i] +=  (D_lamtau_yx_x_i[i] * gradBx_z_r[i]+ D_lamtau_yx_x_r[i] * gradBx_z_i[i]);

                    gamy_XZ_r[i] += (D_sig_yy_x_r[i] * gradCy_z_r[i]- D_sig_yy_x_i[i] * gradCy_z_i[i]);

                    gamy_XZ_i[i] +=  (D_sig_yy_x_i[i] * gradCy_z_r[i]+ D_sig_yy_x_r[i] * gradCy_z_i[i]);

                    gamy_XZ_r[i] += (D_lamtau_yy_x_r[i] * gradBy_z_r[i]- D_lamtau_yy_x_i[i] * gradBy_z_i[i]);

                    gamy_XZ_i[i] +=  (D_lamtau_yy_x_i[i] * gradBy_z_r[i]+ D_lamtau_yy_x_r[i] * gradBy_z_i[i]);

                    gamy_XZ_r[i] += (D_sig_yz_x_r[i] * gradCz_z_r[i]- D_sig_yz_x_i[i] * gradCz_z_i[i]);

                    gamy_XZ_i[i] +=  (D_sig_yz_x_i[i] * gradCz_z_r[i]+ D_sig_yz_x_r[i] * gradCz_z_i[i]);

                    gamy_XZ_r[i] += (D_lamtau_yz_x_r[i] * gradBz_z_r[i]- D_lamtau_yz_x_i[i] * gradBz_z_i[i]);

                    gamy_XZ_i[i] +=  (D_lamtau_yz_x_i[i] * gradBz_z_r[i]+ D_lamtau_yz_x_r[i] * gradBz_z_i[i]);

                    gamy_YX_r[i] = (D_sig_yx_y_r[i] * gradCx_x_r[i]- D_sig_yx_y_i[i] * gradCx_x_i[i]);

                    gamy_YX_i[i] =  (D_sig_yx_y_i[i] * gradCx_x_r[i]+ D_sig_yx_y_r[i] * gradCx_x_i[i]);

                    gamy_YX_r[i] += (D_lamtau_yx_y_r[i] * gradBx_x_r[i]- D_lamtau_yx_y_i[i] * gradBx_x_i[i]);

                    gamy_YX_i[i] +=  (D_lamtau_yx_y_i[i] * gradBx_x_r[i]+ D_lamtau_yx_y_r[i] * gradBx_x_i[i]);

                    gamy_YX_r[i] += (D_sig_yy_y_r[i] * gradCy_x_r[i]- D_sig_yy_y_i[i] * gradCy_x_i[i]);

                    gamy_YX_i[i] +=  (D_sig_yy_y_i[i] * gradCy_x_r[i]+ D_sig_yy_y_r[i] * gradCy_x_i[i]);

                    gamy_YX_r[i] += (D_lamtau_yy_y_r[i] * gradBy_x_r[i]- D_lamtau_yy_y_i[i] * gradBy_x_i[i]);

                    gamy_YX_i[i] +=  (D_lamtau_yy_y_i[i] * gradBy_x_r[i]+ D_lamtau_yy_y_r[i] * gradBy_x_i[i]);

                    gamy_YX_r[i] += (D_sig_yz_y_r[i] * gradCz_x_r[i]- D_sig_yz_y_i[i] * gradCz_x_i[i]);

                    gamy_YX_i[i] +=  (D_sig_yz_y_i[i] * gradCz_x_r[i]+ D_sig_yz_y_r[i] * gradCz_x_i[i]);

                    gamy_YX_r[i] += (D_lamtau_yz_y_r[i] * gradBz_x_r[i]- D_lamtau_yz_y_i[i] * gradBz_x_i[i]);

                    gamy_YX_i[i] +=  (D_lamtau_yz_y_i[i] * gradBz_x_r[i]+ D_lamtau_yz_y_r[i] * gradBz_x_i[i]);

                    gamy_YY_r[i] = (D_sig_yx_y_r[i] * gradCx_y_r[i]- D_sig_yx_y_i[i] * gradCx_y_i[i]);

                    gamy_YY_i[i] =  (D_sig_yx_y_i[i] * gradCx_y_r[i]+ D_sig_yx_y_r[i] * gradCx_y_i[i]);

                    gamy_YY_r[i] += (D_lamtau_yx_y_r[i] * gradBx_y_r[i]- D_lamtau_yx_y_i[i] * gradBx_y_i[i]);

                    gamy_YY_i[i] +=  (D_lamtau_yx_y_i[i] * gradBx_y_r[i]+ D_lamtau_yx_y_r[i] * gradBx_y_i[i]);

                    gamy_YY_r[i] += (D_sig_yy_y_r[i] * gradCy_y_r[i]- D_sig_yy_y_i[i] * gradCy_y_i[i]);

                    gamy_YY_i[i] +=  (D_sig_yy_y_i[i] * gradCy_y_r[i]+ D_sig_yy_y_r[i] * gradCy_y_i[i]);

                    gamy_YY_r[i] += (D_lamtau_yy_y_r[i] * gradBy_y_r[i]- D_lamtau_yy_y_i[i] * gradBy_y_i[i]);

                    gamy_YY_i[i] +=  (D_lamtau_yy_y_i[i] * gradBy_y_r[i]+ D_lamtau_yy_y_r[i] * gradBy_y_i[i]);

                    gamy_YY_r[i] += (D_sig_yz_y_r[i] * gradCz_y_r[i]- D_sig_yz_y_i[i] * gradCz_y_i[i]);

                    gamy_YY_i[i] +=  (D_sig_yz_y_i[i] * gradCz_y_r[i]+ D_sig_yz_y_r[i] * gradCz_y_i[i]);

                    gamy_YY_r[i] += (D_lamtau_yz_y_r[i] * gradBz_y_r[i]- D_lamtau_yz_y_i[i] * gradBz_y_i[i]);

                    gamy_YY_i[i] +=  (D_lamtau_yz_y_i[i] * gradBz_y_r[i]+ D_lamtau_yz_y_r[i] * gradBz_y_i[i]);

                    gamy_YZ_r[i] = (D_sig_yx_y_r[i] * gradCx_z_r[i]- D_sig_yx_y_i[i] * gradCx_z_i[i]);

                    gamy_YZ_i[i] =  (D_sig_yx_y_i[i] * gradCx_z_r[i]+ D_sig_yx_y_r[i] * gradCx_z_i[i]);

                    gamy_YZ_r[i] += (D_lamtau_yx_y_r[i] * gradBx_z_r[i]- D_lamtau_yx_y_i[i] * gradBx_z_i[i]);

                    gamy_YZ_i[i] +=  (D_lamtau_yx_y_i[i] * gradBx_z_r[i]+ D_lamtau_yx_y_r[i] * gradBx_z_i[i]);

                    gamy_YZ_r[i] += (D_sig_yy_y_r[i] * gradCy_z_r[i]- D_sig_yy_y_i[i] * gradCy_z_i[i]);

                    gamy_YZ_i[i] +=  (D_sig_yy_y_i[i] * gradCy_z_r[i]+ D_sig_yy_y_r[i] * gradCy_z_i[i]);

                    gamy_YZ_r[i] += (D_lamtau_yy_y_r[i] * gradBy_z_r[i]- D_lamtau_yy_y_i[i] * gradBy_z_i[i]);

                    gamy_YZ_i[i] +=  (D_lamtau_yy_y_i[i] * gradBy_z_r[i]+ D_lamtau_yy_y_r[i] * gradBy_z_i[i]);

                    gamy_YZ_r[i] += (D_sig_yz_y_r[i] * gradCz_z_r[i]- D_sig_yz_y_i[i] * gradCz_z_i[i]);

                    gamy_YZ_i[i] +=  (D_sig_yz_y_i[i] * gradCz_z_r[i]+ D_sig_yz_y_r[i] * gradCz_z_i[i]);

                    gamy_YZ_r[i] += (D_lamtau_yz_y_r[i] * gradBz_z_r[i]- D_lamtau_yz_y_i[i] * gradBz_z_i[i]);

                    gamy_YZ_i[i] +=  (D_lamtau_yz_y_i[i] * gradBz_z_r[i]+ D_lamtau_yz_y_r[i] * gradBz_z_i[i]);

                    gamy_ZX_r[i] = (D_sig_yx_z_r[i] * gradCx_x_r[i]- D_sig_yx_z_i[i] * gradCx_x_i[i]);

                    gamy_ZX_i[i] =  (D_sig_yx_z_i[i] * gradCx_x_r[i]+ D_sig_yx_z_r[i] * gradCx_x_i[i]);

                    gamy_ZX_r[i] += (D_lamtau_yx_z_r[i] * gradBx_x_r[i]- D_lamtau_yx_z_i[i] * gradBx_x_i[i]);

                    gamy_ZX_i[i] +=  (D_lamtau_yx_z_i[i] * gradBx_x_r[i]+ D_lamtau_yx_z_r[i] * gradBx_x_i[i]);

                    gamy_ZX_r[i] += (D_sig_yy_z_r[i] * gradCy_x_r[i]- D_sig_yy_z_i[i] * gradCy_x_i[i]);

                    gamy_ZX_i[i] +=  (D_sig_yy_z_i[i] * gradCy_x_r[i]+ D_sig_yy_z_r[i] * gradCy_x_i[i]);

                    gamy_ZX_r[i] += (D_lamtau_yy_z_r[i] * gradBy_x_r[i]- D_lamtau_yy_z_i[i] * gradBy_x_i[i]);

                    gamy_ZX_i[i] +=  (D_lamtau_yy_z_i[i] * gradBy_x_r[i]+ D_lamtau_yy_z_r[i] * gradBy_x_i[i]);

                    gamy_ZX_r[i] += (D_sig_yz_z_r[i] * gradCz_x_r[i]- D_sig_yz_z_i[i] * gradCz_x_i[i]);

                    gamy_ZX_i[i] +=  (D_sig_yz_z_i[i] * gradCz_x_r[i]+ D_sig_yz_z_r[i] * gradCz_x_i[i]);

                    gamy_ZX_r[i] += (D_lamtau_yz_z_r[i] * gradBz_x_r[i]- D_lamtau_yz_z_i[i] * gradBz_x_i[i]);

                    gamy_ZX_i[i] +=  (D_lamtau_yz_z_i[i] * gradBz_x_r[i]+ D_lamtau_yz_z_r[i] * gradBz_x_i[i]);

                    gamy_ZY_r[i] = (D_sig_yx_z_r[i] * gradCx_y_r[i]- D_sig_yx_z_i[i] * gradCx_y_i[i]);

                    gamy_ZY_i[i] =  (D_sig_yx_z_i[i] * gradCx_y_r[i]+ D_sig_yx_z_r[i] * gradCx_y_i[i]);

                    gamy_ZY_r[i] += (D_lamtau_yx_z_r[i] * gradBx_y_r[i]- D_lamtau_yx_z_i[i] * gradBx_y_i[i]);

                    gamy_ZY_i[i] +=  (D_lamtau_yx_z_i[i] * gradBx_y_r[i]+ D_lamtau_yx_z_r[i] * gradBx_y_i[i]);

                    gamy_ZY_r[i] += (D_sig_yy_z_r[i] * gradCy_y_r[i]- D_sig_yy_z_i[i] * gradCy_y_i[i]);

                    gamy_ZY_i[i] +=  (D_sig_yy_z_i[i] * gradCy_y_r[i]+ D_sig_yy_z_r[i] * gradCy_y_i[i]);

                    gamy_ZY_r[i] += (D_lamtau_yy_z_r[i] * gradBy_y_r[i]- D_lamtau_yy_z_i[i] * gradBy_y_i[i]);

                    gamy_ZY_i[i] +=  (D_lamtau_yy_z_i[i] * gradBy_y_r[i]+ D_lamtau_yy_z_r[i] * gradBy_y_i[i]);

                    gamy_ZY_r[i] += (D_sig_yz_z_r[i] * gradCz_y_r[i]- D_sig_yz_z_i[i] * gradCz_y_i[i]);

                    gamy_ZY_i[i] +=  (D_sig_yz_z_i[i] * gradCz_y_r[i]+ D_sig_yz_z_r[i] * gradCz_y_i[i]);

                    gamy_ZY_r[i] += (D_lamtau_yz_z_r[i] * gradBz_y_r[i]- D_lamtau_yz_z_i[i] * gradBz_y_i[i]);

                    gamy_ZY_i[i] +=  (D_lamtau_yz_z_i[i] * gradBz_y_r[i]+ D_lamtau_yz_z_r[i] * gradBz_y_i[i]);

                    gamy_ZZ_r[i] = (D_sig_yx_z_r[i] * gradCx_z_r[i]- D_sig_yx_z_i[i] * gradCx_z_i[i]);

                    gamy_ZZ_i[i] =  (D_sig_yx_z_i[i] * gradCx_z_r[i]+ D_sig_yx_z_r[i] * gradCx_z_i[i]);

                    gamy_ZZ_r[i] += (D_lamtau_yx_z_r[i] * gradBx_z_r[i]- D_lamtau_yx_z_i[i] * gradBx_z_i[i]);

                    gamy_ZZ_i[i] +=  (D_lamtau_yx_z_i[i] * gradBx_z_r[i]+ D_lamtau_yx_z_r[i] * gradBx_z_i[i]);

                    gamy_ZZ_r[i] += (D_sig_yy_z_r[i] * gradCy_z_r[i]- D_sig_yy_z_i[i] * gradCy_z_i[i]);

                    gamy_ZZ_i[i] +=  (D_sig_yy_z_i[i] * gradCy_z_r[i]+ D_sig_yy_z_r[i] * gradCy_z_i[i]);

                    gamy_ZZ_r[i] += (D_lamtau_yy_z_r[i] * gradBy_z_r[i]- D_lamtau_yy_z_i[i] * gradBy_z_i[i]);

                    gamy_ZZ_i[i] +=  (D_lamtau_yy_z_i[i] * gradBy_z_r[i]+ D_lamtau_yy_z_r[i] * gradBy_z_i[i]);

                    gamy_ZZ_r[i] += (D_sig_yz_z_r[i] * gradCz_z_r[i]- D_sig_yz_z_i[i] * gradCz_z_i[i]);

                    gamy_ZZ_i[i] +=  (D_sig_yz_z_i[i] * gradCz_z_r[i]+ D_sig_yz_z_r[i] * gradCz_z_i[i]);

                    gamy_ZZ_r[i] += (D_lamtau_yz_z_r[i] * gradBz_z_r[i]- D_lamtau_yz_z_i[i] * gradBz_z_i[i]);

                    gamy_ZZ_i[i] +=  (D_lamtau_yz_z_i[i] * gradBz_z_r[i]+ D_lamtau_yz_z_r[i] * gradBz_z_i[i]);

                    gamz_XX_r[i] = (D_sig_zx_x_r[i] * gradCx_x_r[i]- D_sig_zx_x_i[i] * gradCx_x_i[i]);

                    gamz_XX_i[i] =  (D_sig_zx_x_i[i] * gradCx_x_r[i]+ D_sig_zx_x_r[i] * gradCx_x_i[i]);

                    gamz_XX_r[i] += (D_lamtau_zx_x_r[i] * gradBx_x_r[i]- D_lamtau_zx_x_i[i] * gradBx_x_i[i]);

                    gamz_XX_i[i] +=  (D_lamtau_zx_x_i[i] * gradBx_x_r[i]+ D_lamtau_zx_x_r[i] * gradBx_x_i[i]);

                    gamz_XX_r[i] += (D_sig_zy_x_r[i] * gradCy_x_r[i]- D_sig_zy_x_i[i] * gradCy_x_i[i]);

                    gamz_XX_i[i] +=  (D_sig_zy_x_i[i] * gradCy_x_r[i]+ D_sig_zy_x_r[i] * gradCy_x_i[i]);

                    gamz_XX_r[i] += (D_lamtau_zy_x_r[i] * gradBy_x_r[i]- D_lamtau_zy_x_i[i] * gradBy_x_i[i]);

                    gamz_XX_i[i] +=  (D_lamtau_zy_x_i[i] * gradBy_x_r[i]+ D_lamtau_zy_x_r[i] * gradBy_x_i[i]);

                    gamz_XX_r[i] += (D_sig_zz_x_r[i] * gradCz_x_r[i]- D_sig_zz_x_i[i] * gradCz_x_i[i]);

                    gamz_XX_i[i] +=  (D_sig_zz_x_i[i] * gradCz_x_r[i]+ D_sig_zz_x_r[i] * gradCz_x_i[i]);

                    gamz_XX_r[i] += (D_lamtau_zz_x_r[i] * gradBz_x_r[i]- D_lamtau_zz_x_i[i] * gradBz_x_i[i]);

                    gamz_XX_i[i] +=  (D_lamtau_zz_x_i[i] * gradBz_x_r[i]+ D_lamtau_zz_x_r[i] * gradBz_x_i[i]);

                    gamz_XY_r[i] = (D_sig_zx_x_r[i] * gradCx_y_r[i]- D_sig_zx_x_i[i] * gradCx_y_i[i]);

                    gamz_XY_i[i] =  (D_sig_zx_x_i[i] * gradCx_y_r[i]+ D_sig_zx_x_r[i] * gradCx_y_i[i]);

                    gamz_XY_r[i] += (D_lamtau_zx_x_r[i] * gradBx_y_r[i]- D_lamtau_zx_x_i[i] * gradBx_y_i[i]);

                    gamz_XY_i[i] +=  (D_lamtau_zx_x_i[i] * gradBx_y_r[i]+ D_lamtau_zx_x_r[i] * gradBx_y_i[i]);

                    gamz_XY_r[i] += (D_sig_zy_x_r[i] * gradCy_y_r[i]- D_sig_zy_x_i[i] * gradCy_y_i[i]);

                    gamz_XY_i[i] +=  (D_sig_zy_x_i[i] * gradCy_y_r[i]+ D_sig_zy_x_r[i] * gradCy_y_i[i]);

                    gamz_XY_r[i] += (D_lamtau_zy_x_r[i] * gradBy_y_r[i]- D_lamtau_zy_x_i[i] * gradBy_y_i[i]);

                    gamz_XY_i[i] +=  (D_lamtau_zy_x_i[i] * gradBy_y_r[i]+ D_lamtau_zy_x_r[i] * gradBy_y_i[i]);

                    gamz_XY_r[i] += (D_sig_zz_x_r[i] * gradCz_y_r[i]- D_sig_zz_x_i[i] * gradCz_y_i[i]);

                    gamz_XY_i[i] +=  (D_sig_zz_x_i[i] * gradCz_y_r[i]+ D_sig_zz_x_r[i] * gradCz_y_i[i]);

                    gamz_XY_r[i] += (D_lamtau_zz_x_r[i] * gradBz_y_r[i]- D_lamtau_zz_x_i[i] * gradBz_y_i[i]);

                    gamz_XY_i[i] +=  (D_lamtau_zz_x_i[i] * gradBz_y_r[i]+ D_lamtau_zz_x_r[i] * gradBz_y_i[i]);

                    gamz_XZ_r[i] = (D_sig_zx_x_r[i] * gradCx_z_r[i]- D_sig_zx_x_i[i] * gradCx_z_i[i]);

                    gamz_XZ_i[i] =  (D_sig_zx_x_i[i] * gradCx_z_r[i]+ D_sig_zx_x_r[i] * gradCx_z_i[i]);

                    gamz_XZ_r[i] += (D_lamtau_zx_x_r[i] * gradBx_z_r[i]- D_lamtau_zx_x_i[i] * gradBx_z_i[i]);

                    gamz_XZ_i[i] +=  (D_lamtau_zx_x_i[i] * gradBx_z_r[i]+ D_lamtau_zx_x_r[i] * gradBx_z_i[i]);

                    gamz_XZ_r[i] += (D_sig_zy_x_r[i] * gradCy_z_r[i]- D_sig_zy_x_i[i] * gradCy_z_i[i]);

                    gamz_XZ_i[i] +=  (D_sig_zy_x_i[i] * gradCy_z_r[i]+ D_sig_zy_x_r[i] * gradCy_z_i[i]);

                    gamz_XZ_r[i] += (D_lamtau_zy_x_r[i] * gradBy_z_r[i]- D_lamtau_zy_x_i[i] * gradBy_z_i[i]);

                    gamz_XZ_i[i] +=  (D_lamtau_zy_x_i[i] * gradBy_z_r[i]+ D_lamtau_zy_x_r[i] * gradBy_z_i[i]);

                    gamz_XZ_r[i] += (D_sig_zz_x_r[i] * gradCz_z_r[i]- D_sig_zz_x_i[i] * gradCz_z_i[i]);

                    gamz_XZ_i[i] +=  (D_sig_zz_x_i[i] * gradCz_z_r[i]+ D_sig_zz_x_r[i] * gradCz_z_i[i]);

                    gamz_XZ_r[i] += (D_lamtau_zz_x_r[i] * gradBz_z_r[i]- D_lamtau_zz_x_i[i] * gradBz_z_i[i]);

                    gamz_XZ_i[i] +=  (D_lamtau_zz_x_i[i] * gradBz_z_r[i]+ D_lamtau_zz_x_r[i] * gradBz_z_i[i]);

                    gamz_YX_r[i] = (D_sig_zx_y_r[i] * gradCx_x_r[i]- D_sig_zx_y_i[i] * gradCx_x_i[i]);

                    gamz_YX_i[i] =  (D_sig_zx_y_i[i] * gradCx_x_r[i]+ D_sig_zx_y_r[i] * gradCx_x_i[i]);

                    gamz_YX_r[i] += (D_lamtau_zx_y_r[i] * gradBx_x_r[i]- D_lamtau_zx_y_i[i] * gradBx_x_i[i]);

                    gamz_YX_i[i] +=  (D_lamtau_zx_y_i[i] * gradBx_x_r[i]+ D_lamtau_zx_y_r[i] * gradBx_x_i[i]);

                    gamz_YX_r[i] += (D_sig_zy_y_r[i] * gradCy_x_r[i]- D_sig_zy_y_i[i] * gradCy_x_i[i]);

                    gamz_YX_i[i] +=  (D_sig_zy_y_i[i] * gradCy_x_r[i]+ D_sig_zy_y_r[i] * gradCy_x_i[i]);

                    gamz_YX_r[i] += (D_lamtau_zy_y_r[i] * gradBy_x_r[i]- D_lamtau_zy_y_i[i] * gradBy_x_i[i]);

                    gamz_YX_i[i] +=  (D_lamtau_zy_y_i[i] * gradBy_x_r[i]+ D_lamtau_zy_y_r[i] * gradBy_x_i[i]);

                    gamz_YX_r[i] += (D_sig_zz_y_r[i] * gradCz_x_r[i]- D_sig_zz_y_i[i] * gradCz_x_i[i]);

                    gamz_YX_i[i] +=  (D_sig_zz_y_i[i] * gradCz_x_r[i]+ D_sig_zz_y_r[i] * gradCz_x_i[i]);

                    gamz_YX_r[i] += (D_lamtau_zz_y_r[i] * gradBz_x_r[i]- D_lamtau_zz_y_i[i] * gradBz_x_i[i]);

                    gamz_YX_i[i] +=  (D_lamtau_zz_y_i[i] * gradBz_x_r[i]+ D_lamtau_zz_y_r[i] * gradBz_x_i[i]);

                    gamz_YY_r[i] = (D_sig_zx_y_r[i] * gradCx_y_r[i]- D_sig_zx_y_i[i] * gradCx_y_i[i]);

                    gamz_YY_i[i] =  (D_sig_zx_y_i[i] * gradCx_y_r[i]+ D_sig_zx_y_r[i] * gradCx_y_i[i]);

                    gamz_YY_r[i] += (D_lamtau_zx_y_r[i] * gradBx_y_r[i]- D_lamtau_zx_y_i[i] * gradBx_y_i[i]);

                    gamz_YY_i[i] +=  (D_lamtau_zx_y_i[i] * gradBx_y_r[i]+ D_lamtau_zx_y_r[i] * gradBx_y_i[i]);

                    gamz_YY_r[i] += (D_sig_zy_y_r[i] * gradCy_y_r[i]- D_sig_zy_y_i[i] * gradCy_y_i[i]);

                    gamz_YY_i[i] +=  (D_sig_zy_y_i[i] * gradCy_y_r[i]+ D_sig_zy_y_r[i] * gradCy_y_i[i]);

                    gamz_YY_r[i] += (D_lamtau_zy_y_r[i] * gradBy_y_r[i]- D_lamtau_zy_y_i[i] * gradBy_y_i[i]);

                    gamz_YY_i[i] +=  (D_lamtau_zy_y_i[i] * gradBy_y_r[i]+ D_lamtau_zy_y_r[i] * gradBy_y_i[i]);

                    gamz_YY_r[i] += (D_sig_zz_y_r[i] * gradCz_y_r[i]- D_sig_zz_y_i[i] * gradCz_y_i[i]);

                    gamz_YY_i[i] +=  (D_sig_zz_y_i[i] * gradCz_y_r[i]+ D_sig_zz_y_r[i] * gradCz_y_i[i]);

                    gamz_YY_r[i] += (D_lamtau_zz_y_r[i] * gradBz_y_r[i]- D_lamtau_zz_y_i[i] * gradBz_y_i[i]);

                    gamz_YY_i[i] +=  (D_lamtau_zz_y_i[i] * gradBz_y_r[i]+ D_lamtau_zz_y_r[i] * gradBz_y_i[i]);

                    gamz_YZ_r[i] = (D_sig_zx_y_r[i] * gradCx_z_r[i]- D_sig_zx_y_i[i] * gradCx_z_i[i]);

                    gamz_YZ_i[i] =  (D_sig_zx_y_i[i] * gradCx_z_r[i]+ D_sig_zx_y_r[i] * gradCx_z_i[i]);

                    gamz_YZ_r[i] += (D_lamtau_zx_y_r[i] * gradBx_z_r[i]- D_lamtau_zx_y_i[i] * gradBx_z_i[i]);

                    gamz_YZ_i[i] +=  (D_lamtau_zx_y_i[i] * gradBx_z_r[i]+ D_lamtau_zx_y_r[i] * gradBx_z_i[i]);

                    gamz_YZ_r[i] += (D_sig_zy_y_r[i] * gradCy_z_r[i]- D_sig_zy_y_i[i] * gradCy_z_i[i]);

                    gamz_YZ_i[i] +=  (D_sig_zy_y_i[i] * gradCy_z_r[i]+ D_sig_zy_y_r[i] * gradCy_z_i[i]);

                    gamz_YZ_r[i] += (D_lamtau_zy_y_r[i] * gradBy_z_r[i]- D_lamtau_zy_y_i[i] * gradBy_z_i[i]);

                    gamz_YZ_i[i] +=  (D_lamtau_zy_y_i[i] * gradBy_z_r[i]+ D_lamtau_zy_y_r[i] * gradBy_z_i[i]);

                    gamz_YZ_r[i] += (D_sig_zz_y_r[i] * gradCz_z_r[i]- D_sig_zz_y_i[i] * gradCz_z_i[i]);

                    gamz_YZ_i[i] +=  (D_sig_zz_y_i[i] * gradCz_z_r[i]+ D_sig_zz_y_r[i] * gradCz_z_i[i]);

                    gamz_YZ_r[i] += (D_lamtau_zz_y_r[i] * gradBz_z_r[i]- D_lamtau_zz_y_i[i] * gradBz_z_i[i]);

                    gamz_YZ_i[i] +=  (D_lamtau_zz_y_i[i] * gradBz_z_r[i]+ D_lamtau_zz_y_r[i] * gradBz_z_i[i]);

                    gamz_ZX_r[i] = (D_sig_zx_z_r[i] * gradCx_x_r[i]- D_sig_zx_z_i[i] * gradCx_x_i[i]);

                    gamz_ZX_i[i] =  (D_sig_zx_z_i[i] * gradCx_x_r[i]+ D_sig_zx_z_r[i] * gradCx_x_i[i]);

                    gamz_ZX_r[i] += (D_lamtau_zx_z_r[i] * gradBx_x_r[i]- D_lamtau_zx_z_i[i] * gradBx_x_i[i]);

                    gamz_ZX_i[i] +=  (D_lamtau_zx_z_i[i] * gradBx_x_r[i]+ D_lamtau_zx_z_r[i] * gradBx_x_i[i]);

                    gamz_ZX_r[i] += (D_sig_zy_z_r[i] * gradCy_x_r[i]- D_sig_zy_z_i[i] * gradCy_x_i[i]);

                    gamz_ZX_i[i] +=  (D_sig_zy_z_i[i] * gradCy_x_r[i]+ D_sig_zy_z_r[i] * gradCy_x_i[i]);

                    gamz_ZX_r[i] += (D_lamtau_zy_z_r[i] * gradBy_x_r[i]- D_lamtau_zy_z_i[i] * gradBy_x_i[i]);

                    gamz_ZX_i[i] +=  (D_lamtau_zy_z_i[i] * gradBy_x_r[i]+ D_lamtau_zy_z_r[i] * gradBy_x_i[i]);

                    gamz_ZX_r[i] += (D_sig_zz_z_r[i] * gradCz_x_r[i]- D_sig_zz_z_i[i] * gradCz_x_i[i]);

                    gamz_ZX_i[i] +=  (D_sig_zz_z_i[i] * gradCz_x_r[i]+ D_sig_zz_z_r[i] * gradCz_x_i[i]);

                    gamz_ZX_r[i] += (D_lamtau_zz_z_r[i] * gradBz_x_r[i]- D_lamtau_zz_z_i[i] * gradBz_x_i[i]);

                    gamz_ZX_i[i] +=  (D_lamtau_zz_z_i[i] * gradBz_x_r[i]+ D_lamtau_zz_z_r[i] * gradBz_x_i[i]);

                    gamz_ZY_r[i] = (D_sig_zx_z_r[i] * gradCx_y_r[i]- D_sig_zx_z_i[i] * gradCx_y_i[i]);

                    gamz_ZY_i[i] =  (D_sig_zx_z_i[i] * gradCx_y_r[i]+ D_sig_zx_z_r[i] * gradCx_y_i[i]);

                    gamz_ZY_r[i] += (D_lamtau_zx_z_r[i] * gradBx_y_r[i]- D_lamtau_zx_z_i[i] * gradBx_y_i[i]);

                    gamz_ZY_i[i] +=  (D_lamtau_zx_z_i[i] * gradBx_y_r[i]+ D_lamtau_zx_z_r[i] * gradBx_y_i[i]);

                    gamz_ZY_r[i] += (D_sig_zy_z_r[i] * gradCy_y_r[i]- D_sig_zy_z_i[i] * gradCy_y_i[i]);

                    gamz_ZY_i[i] +=  (D_sig_zy_z_i[i] * gradCy_y_r[i]+ D_sig_zy_z_r[i] * gradCy_y_i[i]);

                    gamz_ZY_r[i] += (D_lamtau_zy_z_r[i] * gradBy_y_r[i]- D_lamtau_zy_z_i[i] * gradBy_y_i[i]);

                    gamz_ZY_i[i] +=  (D_lamtau_zy_z_i[i] * gradBy_y_r[i]+ D_lamtau_zy_z_r[i] * gradBy_y_i[i]);

                    gamz_ZY_r[i] += (D_sig_zz_z_r[i] * gradCz_y_r[i]- D_sig_zz_z_i[i] * gradCz_y_i[i]);

                    gamz_ZY_i[i] +=  (D_sig_zz_z_i[i] * gradCz_y_r[i]+ D_sig_zz_z_r[i] * gradCz_y_i[i]);

                    gamz_ZY_r[i] += (D_lamtau_zz_z_r[i] * gradBz_y_r[i]- D_lamtau_zz_z_i[i] * gradBz_y_i[i]);

                    gamz_ZY_i[i] +=  (D_lamtau_zz_z_i[i] * gradBz_y_r[i]+ D_lamtau_zz_z_r[i] * gradBz_y_i[i]);

                    gamz_ZZ_r[i] = (D_sig_zx_z_r[i] * gradCx_z_r[i]- D_sig_zx_z_i[i] * gradCx_z_i[i]);

                    gamz_ZZ_i[i] =  (D_sig_zx_z_i[i] * gradCx_z_r[i]+ D_sig_zx_z_r[i] * gradCx_z_i[i]);

                    gamz_ZZ_r[i] += (D_lamtau_zx_z_r[i] * gradBx_z_r[i]- D_lamtau_zx_z_i[i] * gradBx_z_i[i]);

                    gamz_ZZ_i[i] +=  (D_lamtau_zx_z_i[i] * gradBx_z_r[i]+ D_lamtau_zx_z_r[i] * gradBx_z_i[i]);

                    gamz_ZZ_r[i] += (D_sig_zy_z_r[i] * gradCy_z_r[i]- D_sig_zy_z_i[i] * gradCy_z_i[i]);

                    gamz_ZZ_i[i] +=  (D_sig_zy_z_i[i] * gradCy_z_r[i]+ D_sig_zy_z_r[i] * gradCy_z_i[i]);

                    gamz_ZZ_r[i] += (D_lamtau_zy_z_r[i] * gradBy_z_r[i]- D_lamtau_zy_z_i[i] * gradBy_z_i[i]);

                    gamz_ZZ_i[i] +=  (D_lamtau_zy_z_i[i] * gradBy_z_r[i]+ D_lamtau_zy_z_r[i] * gradBy_z_i[i]);

                    gamz_ZZ_r[i] += (D_sig_zz_z_r[i] * gradCz_z_r[i]- D_sig_zz_z_i[i] * gradCz_z_i[i]);

                    gamz_ZZ_i[i] +=  (D_sig_zz_z_i[i] * gradCz_z_r[i]+ D_sig_zz_z_r[i] * gradCz_z_i[i]);

                    gamz_ZZ_r[i] += (D_lamtau_zz_z_r[i] * gradBz_z_r[i]- D_lamtau_zz_z_i[i] * gradBz_z_i[i]);

                    gamz_ZZ_i[i] +=  (D_lamtau_zz_z_i[i] * gradBz_z_r[i]+ D_lamtau_zz_z_r[i] * gradBz_z_i[i]);

                }
            }
        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                
                auto gam_r = gam(2 * j);
                auto gam_i = gam(2 * j + 1);
                auto gam_X_r = gamX(2 * j);
                auto gam_X_i = gamX(2 * j + 1);
                auto gam_Y_r = gamY(2 * j);
                auto gam_Y_i = gamY(2 * j + 1);
                auto gam_Z_r = gamZ(2 * j);
                auto gam_Z_i = gamZ(2 * j + 1);
                auto gam_XX_r = gamXX(2 * j);
                auto gam_XX_i = gamXX(2 * j + 1);
                auto gam_XY_r = gamXY(2 * j);
                auto gam_XY_i = gamXY(2 * j + 1);
                auto gam_XZ_r = gamXZ(2 * j);
                auto gam_XZ_i = gamXZ(2 * j + 1);
                auto gam_YX_r = gamYX(2 * j);
                auto gam_YX_i = gamYX(2 * j + 1);
                auto gam_YY_r = gamYY(2 * j);
                auto gam_YY_i = gamYY(2 * j + 1);
                auto gam_YZ_r = gamYZ(2 * j);
                auto gam_YZ_i = gamYZ(2 * j + 1);
                auto gam_ZX_r = gamZX(2 * j);
                auto gam_ZX_i = gamZX(2 * j + 1);
                auto gam_ZY_r = gamZY(2 * j);
                auto gam_ZY_i = gamZY(2 * j + 1);
                auto gam_ZZ_r = gamZZ(2 * j);
                auto gam_ZZ_i = gamZZ(2 * j + 1);


                auto pi_r = pi(2 * j);
                auto pi_i = pi(2 * j + 1);
                
                auto pi_X_r = piX(2 * j);
                auto pi_X_i = piX(2 * j + 1);
                auto pi_Y_r = piY(2 * j);
                auto pi_Y_i = piY(2 * j + 1);
                auto pi_Z_r = piZ(2 * j);
                auto pi_Z_i = piZ(2 * j + 1);
                auto pi_XX_r = piXX(2 * j);
                auto pi_XX_i = piXX(2 * j + 1);
                auto pi_XY_r = piXY(2 * j);
                auto pi_XY_i = piXY(2 * j + 1);
                auto pi_XZ_r = piXZ(2 * j);
                auto pi_XZ_i = piXZ(2 * j + 1);
                auto pi_YX_r = piYX(2 * j);
                auto pi_YX_i = piYX(2 * j + 1);
                auto pi_YY_r = piYY(2 * j);
                auto pi_YY_i = piYY(2 * j + 1);
                auto pi_YZ_r = piYZ(2 * j);
                auto pi_YZ_i = piYZ(2 * j + 1);
                auto pi_ZX_r = piZX(2 * j);
                auto pi_ZX_i = piZX(2 * j + 1);
                auto pi_ZY_r = piZY(2 * j);
                auto pi_ZY_i = piZY(2 * j + 1);
                auto pi_ZZ_r = piZZ(2 * j);
                auto pi_ZZ_i = piZZ(2 * j + 1);

                auto pi_XXX_r=piXXX(2 * j);
                auto pi_XXX_i=piXXX(2 * j + 1);
                auto pi_XXY_r=piXXY(2 * j);
                auto pi_XXY_i=piXXY(2 * j + 1);
                auto pi_XXZ_r=piXXZ(2 * j);
                auto pi_XXZ_i=piXXZ(2 * j + 1);
                auto pi_XYX_r=piXYX(2 * j);
                auto pi_XYX_i=piXYX(2 * j + 1);
                auto pi_XYY_r=piXYY(2 * j);
                auto pi_XYY_i=piXYY(2 * j + 1);
                auto pi_XYZ_r=piXYZ(2 * j);
                auto pi_XYZ_i=piXYZ(2 * j + 1);
                auto pi_XZX_r=piXZX(2 * j);
                auto pi_XZX_i=piXZX(2 * j + 1);
                auto pi_XZY_r=piXZY(2 * j);
                auto pi_XZY_i=piXZY(2 * j + 1);
                auto pi_XZZ_r=piXZZ(2 * j);
                auto pi_XZZ_i=piXZZ(2 * j + 1);
                auto pi_YXX_r=piYXX(2 * j);
                auto pi_YXX_i=piYXX(2 * j + 1);
                auto pi_YXY_r=piYXY(2 * j);
                auto pi_YXY_i=piYXY(2 * j + 1);
                auto pi_YXZ_r=piYXZ(2 * j);
                auto pi_YXZ_i=piYXZ(2 * j + 1);
                auto pi_YYX_r=piYYX(2 * j);
                auto pi_YYX_i=piYYX(2 * j + 1);
                auto pi_YYY_r=piYYY(2 * j);
                auto pi_YYY_i=piYYY(2 * j + 1);
                auto pi_YYZ_r=piYYZ(2 * j);
                auto pi_YYZ_i=piYYZ(2 * j + 1);
                auto pi_YZX_r=piYZX(2 * j);
                auto pi_YZX_i=piYZX(2 * j + 1);
                auto pi_YZY_r=piYZY(2 * j);
                auto pi_YZY_i=piYZY(2 * j + 1);
                auto pi_YZZ_r=piYZZ(2 * j);
                auto pi_YZZ_i=piYZZ(2 * j + 1);
                auto pi_ZXX_r=piZXX(2 * j);
                auto pi_ZXX_i=piZXX(2 * j + 1);
                auto pi_ZXY_r=piZXY(2 * j);
                auto pi_ZXY_i=piZXY(2 * j + 1);
                auto pi_ZXZ_r=piZXZ(2 * j);
                auto pi_ZXZ_i=piZXZ(2 * j + 1);
                auto pi_ZYX_r=piZYX(2 * j);
                auto pi_ZYX_i=piZYX(2 * j + 1);
                auto pi_ZYY_r=piZYY(2 * j);
                auto pi_ZYY_i=piZYY(2 * j + 1);
                auto pi_ZYZ_r=piZYZ(2 * j);
                auto pi_ZYZ_i=piZYZ(2 * j + 1);
                auto pi_ZZX_r=piZZX(2 * j);
                auto pi_ZZX_i=piZZX(2 * j + 1);
                auto pi_ZZY_r=piZZY(2 * j);
                auto pi_ZZY_i=piZZY(2 * j + 1);
                auto pi_ZZZ_r=piZZZ(2 * j);
                auto pi_ZZZ_i=piZZZ(2 * j + 1);

                // First-order densities

                auto rhoB_r = rwDensityGrid.alphaDensity(6 * j);

                auto gradB_x_r = rwDensityGrid.alphaDensityGradientX(6 * j);

                auto gradB_y_r = rwDensityGrid.alphaDensityGradientY(6 * j);

                auto gradB_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j);


                auto rhoB_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto gradB_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 1);

                auto gradB_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 1);

                auto gradB_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 1);


                auto rhoC_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto gradC_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 2);

                auto gradC_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 2);

                auto gradC_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 2);


                auto rhoC_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto gradC_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 3);

                auto gradC_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 3);

                auto gradC_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 3);

                
                auto rhoD_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto gradD_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 4);

                auto gradD_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 4);

                auto gradD_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 4);


                auto rhoD_i = rwDensityGrid.alphaDensity(6 * j + 5);

                auto gradD_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 5);

                auto gradD_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 5);

                auto gradD_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 5);


                auto rhoBC_r = rw2DensityGrid.alphaDensity(6 * j);

                auto gradBC_x_r = rw2DensityGrid.alphaDensityGradientX(6 * j);

                auto gradBC_y_r = rw2DensityGrid.alphaDensityGradientY(6 * j);

                auto gradBC_z_r = rw2DensityGrid.alphaDensityGradientZ(6 * j);


                auto rhoBC_i = rw2DensityGrid.alphaDensity(6 * j + 1);

                auto gradBC_x_i = rw2DensityGrid.alphaDensityGradientX(6 * j + 1);

                auto gradBC_y_i = rw2DensityGrid.alphaDensityGradientY(6 * j + 1);

                auto gradBC_z_i = rw2DensityGrid.alphaDensityGradientZ(6 * j + 1);


                auto rhoBD_r = rw2DensityGrid.alphaDensity(6 * j + 2);

                auto gradBD_x_r = rw2DensityGrid.alphaDensityGradientX(6 * j + 2);

                auto gradBD_y_r = rw2DensityGrid.alphaDensityGradientY(6 * j + 2);

                auto gradBD_z_r = rw2DensityGrid.alphaDensityGradientZ(6 * j + 2);
   

                auto rhoBD_i = rw2DensityGrid.alphaDensity(6 * j + 3);

                auto gradBD_x_i = rw2DensityGrid.alphaDensityGradientX(6 * j + 3);

                auto gradBD_y_i = rw2DensityGrid.alphaDensityGradientY(6 * j + 3);

                auto gradBD_z_i = rw2DensityGrid.alphaDensityGradientZ(6 * j + 3);          


                auto rhoCD_r = rw2DensityGrid.alphaDensity(6 * j + 4);

                auto gradCD_x_r = rw2DensityGrid.alphaDensityGradientX(6 * j + 4);

                auto gradCD_y_r = rw2DensityGrid.alphaDensityGradientY(6 * j + 4);

                auto gradCD_z_r = rw2DensityGrid.alphaDensityGradientZ(6 * j + 4);


                auto rhoCD_i = rw2DensityGrid.alphaDensity(6 * j + 5);

                auto gradCD_x_i = rw2DensityGrid.alphaDensityGradientX(6 * j + 5);

                auto gradCD_y_i = rw2DensityGrid.alphaDensityGradientY(6 * j + 5);

                auto gradCD_z_i = rw2DensityGrid.alphaDensityGradientZ(6 * j + 5);


                // pi = rho_b * rho_c * rho_d

                for (int32_t i = 0; i < npoints; i++)
                {
                    pi_r[i] = 6.0 * (rhoB_r[i] * rhoC_r[i] * rhoD_r[i]
                            -rhoB_i[i] * rhoC_r[i] * rhoD_i[i]
                            -rhoB_r[i] * rhoC_i[i] * rhoD_i[i]
                            -rhoB_i[i] * rhoC_i[i] * rhoD_r[i]);

                    pi_i[i] = 6.0 * ( -rhoB_i[i] * rhoC_i[i] * rhoD_i[i]
                                +rhoB_i[i] * rhoC_r[i] * rhoD_r[i]
                                +rhoB_r[i] * rhoC_i[i] * rhoD_r[i]
                                +rhoB_r[i] * rhoC_r[i] * rhoD_i[i]);

                    pi_X_r[i] = 6.0 * (gradB_x_r[i] * rhoC_r[i] * rhoD_r[i]
                                -gradB_x_i[i] * rhoC_r[i] * rhoD_i[i]
                                -gradB_x_r[i] * rhoC_i[i] * rhoD_i[i]
                                -gradB_x_i[i] * rhoC_i[i] * rhoD_r[i]);

                    pi_X_i[i] = 6.0 * ( -gradB_x_i[i] * rhoC_i[i] * rhoD_i[i]
                                +gradB_x_i[i] * rhoC_r[i] * rhoD_r[i]
                                +gradB_x_r[i] * rhoC_i[i] * rhoD_r[i]
                                +gradB_x_r[i] * rhoC_r[i] * rhoD_i[i]);

                    pi_X_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * rhoD_r[i]
                                -rhoB_i[i] * gradC_x_r[i] * rhoD_i[i]
                                -rhoB_r[i] * gradC_x_i[i] * rhoD_i[i]
                                -rhoB_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    pi_X_i[i] += 6.0 * ( -rhoB_i[i] * gradC_x_i[i] * rhoD_i[i]
                                +rhoB_i[i] * gradC_x_r[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_x_i[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_x_r[i] * rhoD_i[i]);

                    pi_X_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_x_r[i]
                                -rhoB_i[i] * rhoC_r[i] * gradD_x_i[i]
                                -rhoB_r[i] * rhoC_i[i] * gradD_x_i[i]
                                -rhoB_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    pi_X_i[i] += 6.0 * ( -rhoB_i[i] * rhoC_i[i] * gradD_x_i[i]
                                +rhoB_i[i] * rhoC_r[i] * gradD_x_r[i]
                                +rhoB_r[i] * rhoC_i[i] * gradD_x_r[i]
                                +rhoB_r[i] * rhoC_r[i] * gradD_x_i[i]);

                    pi_Y_r[i] = 6.0 * (gradB_y_r[i] * rhoC_r[i] * rhoD_r[i]
                                -gradB_y_i[i] * rhoC_r[i] * rhoD_i[i]
                                -gradB_y_r[i] * rhoC_i[i] * rhoD_i[i]
                                -gradB_y_i[i] * rhoC_i[i] * rhoD_r[i]);

                    pi_Y_i[i] = 6.0 * ( -gradB_y_i[i] * rhoC_i[i] * rhoD_i[i]
                                +gradB_y_i[i] * rhoC_r[i] * rhoD_r[i]
                                +gradB_y_r[i] * rhoC_i[i] * rhoD_r[i]
                                +gradB_y_r[i] * rhoC_r[i] * rhoD_i[i]);

                    pi_Y_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * rhoD_r[i]
                                -rhoB_i[i] * gradC_y_r[i] * rhoD_i[i]
                                -rhoB_r[i] * gradC_y_i[i] * rhoD_i[i]
                                -rhoB_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    pi_Y_i[i] += 6.0 * ( -rhoB_i[i] * gradC_y_i[i] * rhoD_i[i]
                                +rhoB_i[i] * gradC_y_r[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_y_i[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_y_r[i] * rhoD_i[i]);

                    pi_Y_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_y_r[i]
                                -rhoB_i[i] * rhoC_r[i] * gradD_y_i[i]
                                -rhoB_r[i] * rhoC_i[i] * gradD_y_i[i]
                                -rhoB_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    pi_Y_i[i] += 6.0 * ( -rhoB_i[i] * rhoC_i[i] * gradD_y_i[i]
                                +rhoB_i[i] * rhoC_r[i] * gradD_y_r[i]
                                +rhoB_r[i] * rhoC_i[i] * gradD_y_r[i]
                                +rhoB_r[i] * rhoC_r[i] * gradD_y_i[i]);

                    pi_Z_r[i] = 6.0 * (gradB_z_r[i] * rhoC_r[i] * rhoD_r[i]
                                -gradB_z_i[i] * rhoC_r[i] * rhoD_i[i]
                                -gradB_z_r[i] * rhoC_i[i] * rhoD_i[i]
                                -gradB_z_i[i] * rhoC_i[i] * rhoD_r[i]);

                    pi_Z_i[i] = 6.0 * ( -gradB_z_i[i] * rhoC_i[i] * rhoD_i[i]
                                +gradB_z_i[i] * rhoC_r[i] * rhoD_r[i]
                                +gradB_z_r[i] * rhoC_i[i] * rhoD_r[i]
                                +gradB_z_r[i] * rhoC_r[i] * rhoD_i[i]);

                    pi_Z_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * rhoD_r[i]
                                -rhoB_i[i] * gradC_z_r[i] * rhoD_i[i]
                                -rhoB_r[i] * gradC_z_i[i] * rhoD_i[i]
                                -rhoB_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    pi_Z_i[i] += 6.0 * ( -rhoB_i[i] * gradC_z_i[i] * rhoD_i[i]
                                +rhoB_i[i] * gradC_z_r[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_z_i[i] * rhoD_r[i]
                                +rhoB_r[i] * gradC_z_r[i] * rhoD_i[i]);

                    pi_Z_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_z_r[i]
                                -rhoB_i[i] * rhoC_r[i] * gradD_z_i[i]
                                -rhoB_r[i] * rhoC_i[i] * gradD_z_i[i]
                                -rhoB_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    pi_Z_i[i] += 6.0 * ( -rhoB_i[i] * rhoC_i[i] * gradD_z_i[i]
                                +rhoB_i[i] * rhoC_r[i] * gradD_z_r[i]
                                +rhoB_r[i] * rhoC_i[i] * gradD_z_r[i]
                                +rhoB_r[i] * rhoC_r[i] * gradD_z_i[i]);

                    pi_XX_r[i] = 6.0 * (gradB_x_r[i] * gradC_x_r[i] * rhoD_r[i]
                                -gradB_x_i[i] * gradC_x_r[i] * rhoD_i[i]
                                -gradB_x_r[i] * gradC_x_i[i] * rhoD_i[i]
                                -gradB_x_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    pi_XX_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_x_i[i] * rhoD_i[i]
                                +gradB_x_i[i] * gradC_x_r[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_x_i[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_x_r[i] * rhoD_i[i]);

                    pi_XY_r[i] = 6.0 * (gradB_x_r[i] * gradC_y_r[i] * rhoD_r[i]
                                -gradB_x_i[i] * gradC_y_r[i] * rhoD_i[i]
                                -gradB_x_r[i] * gradC_y_i[i] * rhoD_i[i]
                                -gradB_x_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    pi_XY_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_y_i[i] * rhoD_i[i]
                                +gradB_x_i[i] * gradC_y_r[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_y_i[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_y_r[i] * rhoD_i[i]);

                    pi_XZ_r[i] = 6.0 * (gradB_x_r[i] * gradC_z_r[i] * rhoD_r[i]
                                -gradB_x_i[i] * gradC_z_r[i] * rhoD_i[i]
                                -gradB_x_r[i] * gradC_z_i[i] * rhoD_i[i]
                                -gradB_x_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    pi_XZ_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_z_i[i] * rhoD_i[i]
                                +gradB_x_i[i] * gradC_z_r[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_z_i[i] * rhoD_r[i]
                                +gradB_x_r[i] * gradC_z_r[i] * rhoD_i[i]);

                    pi_YX_r[i] = 6.0 * (gradB_y_r[i] * gradC_x_r[i] * rhoD_r[i]
                                -gradB_y_i[i] * gradC_x_r[i] * rhoD_i[i]
                                -gradB_y_r[i] * gradC_x_i[i] * rhoD_i[i]
                                -gradB_y_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    pi_YX_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_x_i[i] * rhoD_i[i]
                                +gradB_y_i[i] * gradC_x_r[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_x_i[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_x_r[i] * rhoD_i[i]);

                    pi_YY_r[i] = 6.0 * (gradB_y_r[i] * gradC_y_r[i] * rhoD_r[i]
                                -gradB_y_i[i] * gradC_y_r[i] * rhoD_i[i]
                                -gradB_y_r[i] * gradC_y_i[i] * rhoD_i[i]
                                -gradB_y_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    pi_YY_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_y_i[i] * rhoD_i[i]
                                +gradB_y_i[i] * gradC_y_r[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_y_i[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_y_r[i] * rhoD_i[i]);

                    pi_YZ_r[i] = 6.0 * (gradB_y_r[i] * gradC_z_r[i] * rhoD_r[i]
                                -gradB_y_i[i] * gradC_z_r[i] * rhoD_i[i]
                                -gradB_y_r[i] * gradC_z_i[i] * rhoD_i[i]
                                -gradB_y_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    pi_YZ_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_z_i[i] * rhoD_i[i]
                                +gradB_y_i[i] * gradC_z_r[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_z_i[i] * rhoD_r[i]
                                +gradB_y_r[i] * gradC_z_r[i] * rhoD_i[i]);

                    pi_ZX_r[i] = 6.0 * (gradB_z_r[i] * gradC_x_r[i] * rhoD_r[i]
                                -gradB_z_i[i] * gradC_x_r[i] * rhoD_i[i]
                                -gradB_z_r[i] * gradC_x_i[i] * rhoD_i[i]
                                -gradB_z_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    pi_ZX_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_x_i[i] * rhoD_i[i]
                                +gradB_z_i[i] * gradC_x_r[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_x_i[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_x_r[i] * rhoD_i[i]);

                    pi_ZY_r[i] = 6.0 * (gradB_z_r[i] * gradC_y_r[i] * rhoD_r[i]
                                -gradB_z_i[i] * gradC_y_r[i] * rhoD_i[i]
                                -gradB_z_r[i] * gradC_y_i[i] * rhoD_i[i]
                                -gradB_z_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    pi_ZY_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_y_i[i] * rhoD_i[i]
                                +gradB_z_i[i] * gradC_y_r[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_y_i[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_y_r[i] * rhoD_i[i]);

                    pi_ZZ_r[i] = 6.0 * (gradB_z_r[i] * gradC_z_r[i] * rhoD_r[i]
                                -gradB_z_i[i] * gradC_z_r[i] * rhoD_i[i]
                                -gradB_z_r[i] * gradC_z_i[i] * rhoD_i[i]
                                -gradB_z_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    pi_ZZ_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_z_i[i] * rhoD_i[i]
                                +gradB_z_i[i] * gradC_z_r[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_z_i[i] * rhoD_r[i]
                                +gradB_z_r[i] * gradC_z_r[i] * rhoD_i[i]);

                    pi_XX_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                -rhoB_i[i] * gradC_x_r[i] * gradD_x_i[i]
                                -rhoB_r[i] * gradC_x_i[i] * gradD_x_i[i]
                                -rhoB_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    pi_XX_i[i] += 6.0 * ( -rhoB_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                +rhoB_i[i] * gradC_x_r[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_x_i[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_x_r[i] * gradD_x_i[i]);

                    pi_XY_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                -rhoB_i[i] * gradC_x_r[i] * gradD_y_i[i]
                                -rhoB_r[i] * gradC_x_i[i] * gradD_y_i[i]
                                -rhoB_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    pi_XY_i[i] += 6.0 * ( -rhoB_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                +rhoB_i[i] * gradC_x_r[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_x_i[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_x_r[i] * gradD_y_i[i]);

                    pi_XZ_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                -rhoB_i[i] * gradC_x_r[i] * gradD_z_i[i]
                                -rhoB_r[i] * gradC_x_i[i] * gradD_z_i[i]
                                -rhoB_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    pi_XZ_i[i] += 6.0 * ( -rhoB_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                +rhoB_i[i] * gradC_x_r[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_x_i[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_x_r[i] * gradD_z_i[i]);

                    pi_YX_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                -rhoB_i[i] * gradC_y_r[i] * gradD_x_i[i]
                                -rhoB_r[i] * gradC_y_i[i] * gradD_x_i[i]
                                -rhoB_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    pi_YX_i[i] += 6.0 * ( -rhoB_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                +rhoB_i[i] * gradC_y_r[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_y_i[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_y_r[i] * gradD_x_i[i]);

                    pi_YY_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                -rhoB_i[i] * gradC_y_r[i] * gradD_y_i[i]
                                -rhoB_r[i] * gradC_y_i[i] * gradD_y_i[i]
                                -rhoB_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    pi_YY_i[i] += 6.0 * ( -rhoB_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                +rhoB_i[i] * gradC_y_r[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_y_i[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_y_r[i] * gradD_y_i[i]);

                    pi_YZ_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                -rhoB_i[i] * gradC_y_r[i] * gradD_z_i[i]
                                -rhoB_r[i] * gradC_y_i[i] * gradD_z_i[i]
                                -rhoB_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    pi_YZ_i[i] += 6.0 * ( -rhoB_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                +rhoB_i[i] * gradC_y_r[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_y_i[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_y_r[i] * gradD_z_i[i]);

                    pi_ZX_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                -rhoB_i[i] * gradC_z_r[i] * gradD_x_i[i]
                                -rhoB_r[i] * gradC_z_i[i] * gradD_x_i[i]
                                -rhoB_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    pi_ZX_i[i] += 6.0 * ( -rhoB_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                +rhoB_i[i] * gradC_z_r[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_z_i[i] * gradD_x_r[i]
                                +rhoB_r[i] * gradC_z_r[i] * gradD_x_i[i]);

                    pi_ZY_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                -rhoB_i[i] * gradC_z_r[i] * gradD_y_i[i]
                                -rhoB_r[i] * gradC_z_i[i] * gradD_y_i[i]
                                -rhoB_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    pi_ZY_i[i] += 6.0 * ( -rhoB_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                +rhoB_i[i] * gradC_z_r[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_z_i[i] * gradD_y_r[i]
                                +rhoB_r[i] * gradC_z_r[i] * gradD_y_i[i]);

                    pi_ZZ_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                -rhoB_i[i] * gradC_z_r[i] * gradD_z_i[i]
                                -rhoB_r[i] * gradC_z_i[i] * gradD_z_i[i]
                                -rhoB_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    pi_ZZ_i[i] += 6.0 * ( -rhoB_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                +rhoB_i[i] * gradC_z_r[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_z_i[i] * gradD_z_r[i]
                                +rhoB_r[i] * gradC_z_r[i] * gradD_z_i[i]);

                    pi_XX_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_x_r[i]
                                -gradB_x_i[i] * rhoC_r[i] * gradD_x_i[i]
                                -gradB_x_r[i] * rhoC_i[i] * gradD_x_i[i]
                                -gradB_x_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    pi_XX_i[i] += 6.0 * ( -gradB_x_i[i] * rhoC_i[i] * gradD_x_i[i]
                                +gradB_x_i[i] * rhoC_r[i] * gradD_x_r[i]
                                +gradB_x_r[i] * rhoC_i[i] * gradD_x_r[i]
                                +gradB_x_r[i] * rhoC_r[i] * gradD_x_i[i]);

                    pi_XY_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_y_r[i]
                                -gradB_x_i[i] * rhoC_r[i] * gradD_y_i[i]
                                -gradB_x_r[i] * rhoC_i[i] * gradD_y_i[i]
                                -gradB_x_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    pi_XY_i[i] += 6.0 * ( -gradB_x_i[i] * rhoC_i[i] * gradD_y_i[i]
                                +gradB_x_i[i] * rhoC_r[i] * gradD_y_r[i]
                                +gradB_x_r[i] * rhoC_i[i] * gradD_y_r[i]
                                +gradB_x_r[i] * rhoC_r[i] * gradD_y_i[i]);

                    pi_XZ_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_z_r[i]
                                -gradB_x_i[i] * rhoC_r[i] * gradD_z_i[i]
                                -gradB_x_r[i] * rhoC_i[i] * gradD_z_i[i]
                                -gradB_x_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    pi_XZ_i[i] += 6.0 * ( -gradB_x_i[i] * rhoC_i[i] * gradD_z_i[i]
                                +gradB_x_i[i] * rhoC_r[i] * gradD_z_r[i]
                                +gradB_x_r[i] * rhoC_i[i] * gradD_z_r[i]
                                +gradB_x_r[i] * rhoC_r[i] * gradD_z_i[i]);

                    pi_YX_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_x_r[i]
                                -gradB_y_i[i] * rhoC_r[i] * gradD_x_i[i]
                                -gradB_y_r[i] * rhoC_i[i] * gradD_x_i[i]
                                -gradB_y_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    pi_YX_i[i] += 6.0 * ( -gradB_y_i[i] * rhoC_i[i] * gradD_x_i[i]
                                +gradB_y_i[i] * rhoC_r[i] * gradD_x_r[i]
                                +gradB_y_r[i] * rhoC_i[i] * gradD_x_r[i]
                                +gradB_y_r[i] * rhoC_r[i] * gradD_x_i[i]);

                    pi_YY_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_y_r[i]
                                -gradB_y_i[i] * rhoC_r[i] * gradD_y_i[i]
                                -gradB_y_r[i] * rhoC_i[i] * gradD_y_i[i]
                                -gradB_y_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    pi_YY_i[i] += 6.0 * ( -gradB_y_i[i] * rhoC_i[i] * gradD_y_i[i]
                                +gradB_y_i[i] * rhoC_r[i] * gradD_y_r[i]
                                +gradB_y_r[i] * rhoC_i[i] * gradD_y_r[i]
                                +gradB_y_r[i] * rhoC_r[i] * gradD_y_i[i]);

                    pi_YZ_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_z_r[i]
                                -gradB_y_i[i] * rhoC_r[i] * gradD_z_i[i]
                                -gradB_y_r[i] * rhoC_i[i] * gradD_z_i[i]
                                -gradB_y_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    pi_YZ_i[i] += 6.0 * ( -gradB_y_i[i] * rhoC_i[i] * gradD_z_i[i]
                                +gradB_y_i[i] * rhoC_r[i] * gradD_z_r[i]
                                +gradB_y_r[i] * rhoC_i[i] * gradD_z_r[i]
                                +gradB_y_r[i] * rhoC_r[i] * gradD_z_i[i]);

                    pi_ZX_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_x_r[i]
                                -gradB_z_i[i] * rhoC_r[i] * gradD_x_i[i]
                                -gradB_z_r[i] * rhoC_i[i] * gradD_x_i[i]
                                -gradB_z_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    pi_ZX_i[i] += 6.0 * ( -gradB_z_i[i] * rhoC_i[i] * gradD_x_i[i]
                                +gradB_z_i[i] * rhoC_r[i] * gradD_x_r[i]
                                +gradB_z_r[i] * rhoC_i[i] * gradD_x_r[i]
                                +gradB_z_r[i] * rhoC_r[i] * gradD_x_i[i]);

                    pi_ZY_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_y_r[i]
                                -gradB_z_i[i] * rhoC_r[i] * gradD_y_i[i]
                                -gradB_z_r[i] * rhoC_i[i] * gradD_y_i[i]
                                -gradB_z_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    pi_ZY_i[i] += 6.0 * ( -gradB_z_i[i] * rhoC_i[i] * gradD_y_i[i]
                                +gradB_z_i[i] * rhoC_r[i] * gradD_y_r[i]
                                +gradB_z_r[i] * rhoC_i[i] * gradD_y_r[i]
                                +gradB_z_r[i] * rhoC_r[i] * gradD_y_i[i]);

                    pi_ZZ_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_z_r[i]
                                -gradB_z_i[i] * rhoC_r[i] * gradD_z_i[i]
                                -gradB_z_r[i] * rhoC_i[i] * gradD_z_i[i]
                                -gradB_z_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    pi_ZZ_i[i] += 6.0 * ( -gradB_z_i[i] * rhoC_i[i] * gradD_z_i[i]
                                +gradB_z_i[i] * rhoC_r[i] * gradD_z_r[i]
                                +gradB_z_r[i] * rhoC_i[i] * gradD_z_r[i]
                                +gradB_z_r[i] * rhoC_r[i] * gradD_z_i[i]);

                    pi_XXX_r[i] = 6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                -gradB_x_i[i] * gradC_x_r[i] * gradD_x_i[i]
                                -gradB_x_r[i] * gradC_x_i[i] * gradD_x_i[i]
                                -gradB_x_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    pi_XXX_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                +gradB_x_i[i] * gradC_x_r[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_x_i[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_x_r[i] * gradD_x_i[i]);

                    pi_XXY_r[i] = 6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                -gradB_x_i[i] * gradC_x_r[i] * gradD_y_i[i]
                                -gradB_x_r[i] * gradC_x_i[i] * gradD_y_i[i]
                                -gradB_x_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    pi_XXY_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                +gradB_x_i[i] * gradC_x_r[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_x_i[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_x_r[i] * gradD_y_i[i]);

                    pi_XXZ_r[i] = 6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                -gradB_x_i[i] * gradC_x_r[i] * gradD_z_i[i]
                                -gradB_x_r[i] * gradC_x_i[i] * gradD_z_i[i]
                                -gradB_x_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    pi_XXZ_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                +gradB_x_i[i] * gradC_x_r[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_x_i[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_x_r[i] * gradD_z_i[i]);

                    pi_XYX_r[i] = 6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                -gradB_x_i[i] * gradC_y_r[i] * gradD_x_i[i]
                                -gradB_x_r[i] * gradC_y_i[i] * gradD_x_i[i]
                                -gradB_x_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    pi_XYX_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                +gradB_x_i[i] * gradC_y_r[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_y_i[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_y_r[i] * gradD_x_i[i]);

                    pi_XYY_r[i] = 6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                -gradB_x_i[i] * gradC_y_r[i] * gradD_y_i[i]
                                -gradB_x_r[i] * gradC_y_i[i] * gradD_y_i[i]
                                -gradB_x_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    pi_XYY_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                +gradB_x_i[i] * gradC_y_r[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_y_i[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_y_r[i] * gradD_y_i[i]);

                    pi_XYZ_r[i] = 6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                -gradB_x_i[i] * gradC_y_r[i] * gradD_z_i[i]
                                -gradB_x_r[i] * gradC_y_i[i] * gradD_z_i[i]
                                -gradB_x_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    pi_XYZ_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                +gradB_x_i[i] * gradC_y_r[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_y_i[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_y_r[i] * gradD_z_i[i]);

                    pi_XZX_r[i] = 6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                -gradB_x_i[i] * gradC_z_r[i] * gradD_x_i[i]
                                -gradB_x_r[i] * gradC_z_i[i] * gradD_x_i[i]
                                -gradB_x_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    pi_XZX_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                +gradB_x_i[i] * gradC_z_r[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_z_i[i] * gradD_x_r[i]
                                +gradB_x_r[i] * gradC_z_r[i] * gradD_x_i[i]);

                    pi_XZY_r[i] = 6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                -gradB_x_i[i] * gradC_z_r[i] * gradD_y_i[i]
                                -gradB_x_r[i] * gradC_z_i[i] * gradD_y_i[i]
                                -gradB_x_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    pi_XZY_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                +gradB_x_i[i] * gradC_z_r[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_z_i[i] * gradD_y_r[i]
                                +gradB_x_r[i] * gradC_z_r[i] * gradD_y_i[i]);

                    pi_XZZ_r[i] = 6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                -gradB_x_i[i] * gradC_z_r[i] * gradD_z_i[i]
                                -gradB_x_r[i] * gradC_z_i[i] * gradD_z_i[i]
                                -gradB_x_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    pi_XZZ_i[i] = 6.0 * ( -gradB_x_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                +gradB_x_i[i] * gradC_z_r[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_z_i[i] * gradD_z_r[i]
                                +gradB_x_r[i] * gradC_z_r[i] * gradD_z_i[i]);

                    pi_YXX_r[i] = 6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                -gradB_y_i[i] * gradC_x_r[i] * gradD_x_i[i]
                                -gradB_y_r[i] * gradC_x_i[i] * gradD_x_i[i]
                                -gradB_y_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    pi_YXX_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                +gradB_y_i[i] * gradC_x_r[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_x_i[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_x_r[i] * gradD_x_i[i]);

                    pi_YXY_r[i] = 6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                -gradB_y_i[i] * gradC_x_r[i] * gradD_y_i[i]
                                -gradB_y_r[i] * gradC_x_i[i] * gradD_y_i[i]
                                -gradB_y_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    pi_YXY_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                +gradB_y_i[i] * gradC_x_r[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_x_i[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_x_r[i] * gradD_y_i[i]);

                    pi_YXZ_r[i] = 6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                -gradB_y_i[i] * gradC_x_r[i] * gradD_z_i[i]
                                -gradB_y_r[i] * gradC_x_i[i] * gradD_z_i[i]
                                -gradB_y_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    pi_YXZ_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                +gradB_y_i[i] * gradC_x_r[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_x_i[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_x_r[i] * gradD_z_i[i]);

                    pi_YYX_r[i] = 6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                -gradB_y_i[i] * gradC_y_r[i] * gradD_x_i[i]
                                -gradB_y_r[i] * gradC_y_i[i] * gradD_x_i[i]
                                -gradB_y_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    pi_YYX_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                +gradB_y_i[i] * gradC_y_r[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_y_i[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_y_r[i] * gradD_x_i[i]);

                    pi_YYY_r[i] = 6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                -gradB_y_i[i] * gradC_y_r[i] * gradD_y_i[i]
                                -gradB_y_r[i] * gradC_y_i[i] * gradD_y_i[i]
                                -gradB_y_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    pi_YYY_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                +gradB_y_i[i] * gradC_y_r[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_y_i[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_y_r[i] * gradD_y_i[i]);

                    pi_YYZ_r[i] = 6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                -gradB_y_i[i] * gradC_y_r[i] * gradD_z_i[i]
                                -gradB_y_r[i] * gradC_y_i[i] * gradD_z_i[i]
                                -gradB_y_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    pi_YYZ_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                +gradB_y_i[i] * gradC_y_r[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_y_i[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_y_r[i] * gradD_z_i[i]);

                    pi_YZX_r[i] = 6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                -gradB_y_i[i] * gradC_z_r[i] * gradD_x_i[i]
                                -gradB_y_r[i] * gradC_z_i[i] * gradD_x_i[i]
                                -gradB_y_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    pi_YZX_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                +gradB_y_i[i] * gradC_z_r[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_z_i[i] * gradD_x_r[i]
                                +gradB_y_r[i] * gradC_z_r[i] * gradD_x_i[i]);

                    pi_YZY_r[i] = 6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                -gradB_y_i[i] * gradC_z_r[i] * gradD_y_i[i]
                                -gradB_y_r[i] * gradC_z_i[i] * gradD_y_i[i]
                                -gradB_y_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    pi_YZY_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                +gradB_y_i[i] * gradC_z_r[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_z_i[i] * gradD_y_r[i]
                                +gradB_y_r[i] * gradC_z_r[i] * gradD_y_i[i]);

                    pi_YZZ_r[i] = 6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                -gradB_y_i[i] * gradC_z_r[i] * gradD_z_i[i]
                                -gradB_y_r[i] * gradC_z_i[i] * gradD_z_i[i]
                                -gradB_y_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    pi_YZZ_i[i] = 6.0 * ( -gradB_y_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                +gradB_y_i[i] * gradC_z_r[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_z_i[i] * gradD_z_r[i]
                                +gradB_y_r[i] * gradC_z_r[i] * gradD_z_i[i]);

                    pi_ZXX_r[i] = 6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                -gradB_z_i[i] * gradC_x_r[i] * gradD_x_i[i]
                                -gradB_z_r[i] * gradC_x_i[i] * gradD_x_i[i]
                                -gradB_z_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    pi_ZXX_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                +gradB_z_i[i] * gradC_x_r[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_x_i[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_x_r[i] * gradD_x_i[i]);

                    pi_ZXY_r[i] = 6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                -gradB_z_i[i] * gradC_x_r[i] * gradD_y_i[i]
                                -gradB_z_r[i] * gradC_x_i[i] * gradD_y_i[i]
                                -gradB_z_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    pi_ZXY_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                +gradB_z_i[i] * gradC_x_r[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_x_i[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_x_r[i] * gradD_y_i[i]);

                    pi_ZXZ_r[i] = 6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                -gradB_z_i[i] * gradC_x_r[i] * gradD_z_i[i]
                                -gradB_z_r[i] * gradC_x_i[i] * gradD_z_i[i]
                                -gradB_z_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    pi_ZXZ_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                +gradB_z_i[i] * gradC_x_r[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_x_i[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_x_r[i] * gradD_z_i[i]);

                    pi_ZYX_r[i] = 6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                -gradB_z_i[i] * gradC_y_r[i] * gradD_x_i[i]
                                -gradB_z_r[i] * gradC_y_i[i] * gradD_x_i[i]
                                -gradB_z_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    pi_ZYX_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                +gradB_z_i[i] * gradC_y_r[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_y_i[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_y_r[i] * gradD_x_i[i]);

                    pi_ZYY_r[i] = 6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                -gradB_z_i[i] * gradC_y_r[i] * gradD_y_i[i]
                                -gradB_z_r[i] * gradC_y_i[i] * gradD_y_i[i]
                                -gradB_z_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    pi_ZYY_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                +gradB_z_i[i] * gradC_y_r[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_y_i[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_y_r[i] * gradD_y_i[i]);

                    pi_ZYZ_r[i] = 6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                -gradB_z_i[i] * gradC_y_r[i] * gradD_z_i[i]
                                -gradB_z_r[i] * gradC_y_i[i] * gradD_z_i[i]
                                -gradB_z_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    pi_ZYZ_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                +gradB_z_i[i] * gradC_y_r[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_y_i[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_y_r[i] * gradD_z_i[i]);

                    pi_ZZX_r[i] = 6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                -gradB_z_i[i] * gradC_z_r[i] * gradD_x_i[i]
                                -gradB_z_r[i] * gradC_z_i[i] * gradD_x_i[i]
                                -gradB_z_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    pi_ZZX_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                +gradB_z_i[i] * gradC_z_r[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_z_i[i] * gradD_x_r[i]
                                +gradB_z_r[i] * gradC_z_r[i] * gradD_x_i[i]);

                    pi_ZZY_r[i] = 6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                -gradB_z_i[i] * gradC_z_r[i] * gradD_y_i[i]
                                -gradB_z_r[i] * gradC_z_i[i] * gradD_y_i[i]
                                -gradB_z_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    pi_ZZY_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                +gradB_z_i[i] * gradC_z_r[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_z_i[i] * gradD_y_r[i]
                                +gradB_z_r[i] * gradC_z_r[i] * gradD_y_i[i]);

                    pi_ZZZ_r[i] = 6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                -gradB_z_i[i] * gradC_z_r[i] * gradD_z_i[i]
                                -gradB_z_r[i] * gradC_z_i[i] * gradD_z_i[i]
                                -gradB_z_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    pi_ZZZ_i[i] = 6.0 * ( -gradB_z_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                +gradB_z_i[i] * gradC_z_r[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_z_i[i] * gradD_z_r[i]
                                +gradB_z_r[i] * gradC_z_r[i] * gradD_z_i[i]);

                    gam_r[i] = 3.0 * (rhoD_r[i] * rhoBC_r[i]- rhoD_i[i] * rhoBC_i[i]);

                    gam_i[i] = 3.0 * (rhoD_i[i] * rhoBC_r[i] + rhoD_r[i] * rhoBC_i[i]);

                    gam_r[i] += 3.0 * (rhoB_r[i] * rhoCD_r[i]- rhoB_i[i] * rhoCD_i[i]);

                    gam_i[i] += 3.0 * (rhoB_i[i] * rhoCD_r[i] + rhoB_r[i] * rhoCD_i[i]);

                    gam_r[i] += 3.0 * (rhoC_r[i] * rhoBD_r[i]- rhoC_i[i] * rhoBD_i[i]);

                    gam_i[i] += 3.0 * (rhoC_i[i] * rhoBD_r[i] + rhoC_r[i] * rhoBD_i[i]);

                    gam_X_r[i] = 3.0 * (rhoD_r[i] * gradBC_x_r[i]- rhoD_i[i] * gradBC_x_i[i]);

                    gam_X_i[i] = 3.0 * (rhoD_i[i] * gradBC_x_r[i] + rhoD_r[i] * gradBC_x_i[i]);

                    gam_X_r[i] += 3.0 * (rhoB_r[i] * gradCD_x_r[i]- rhoB_i[i] * gradCD_x_i[i]);

                    gam_X_i[i] += 3.0 * (rhoB_i[i] * gradCD_x_r[i] + rhoB_r[i] * gradCD_x_i[i]);

                    gam_X_r[i] += 3.0 * (rhoC_r[i] * gradBD_x_r[i]- rhoC_i[i] * gradBD_x_i[i]);

                    gam_X_i[i] += 3.0 * (rhoC_i[i] * gradBD_x_r[i] + rhoC_r[i] * gradBD_x_i[i]);

                    gam_X_r[i] += 3.0 * (gradD_x_r[i] * rhoBC_r[i]- gradD_x_i[i] * rhoBC_i[i]);

                    gam_X_i[i] += 3.0 * (gradD_x_i[i] * rhoBC_r[i] + gradD_x_r[i] * rhoBC_i[i]);

                    gam_X_r[i] += 3.0 * (gradB_x_r[i] * rhoCD_r[i]- gradB_x_i[i] * rhoCD_i[i]);

                    gam_X_i[i] += 3.0 * (gradB_x_i[i] * rhoCD_r[i] + gradB_x_r[i] * rhoCD_i[i]);

                    gam_X_r[i] += 3.0 * (gradC_x_r[i] * rhoBD_r[i]- gradC_x_i[i] * rhoBD_i[i]);

                    gam_X_i[i] += 3.0 * (gradC_x_i[i] * rhoBD_r[i] + gradC_x_r[i] * rhoBD_i[i]);

                    gam_Y_r[i] = 3.0 * (rhoD_r[i] * gradBC_y_r[i]- rhoD_i[i] * gradBC_y_i[i]);

                    gam_Y_i[i] = 3.0 * (rhoD_i[i] * gradBC_y_r[i] + rhoD_r[i] * gradBC_y_i[i]);

                    gam_Y_r[i] += 3.0 * (rhoB_r[i] * gradCD_y_r[i]- rhoB_i[i] * gradCD_y_i[i]);

                    gam_Y_i[i] += 3.0 * (rhoB_i[i] * gradCD_y_r[i] + rhoB_r[i] * gradCD_y_i[i]);

                    gam_Y_r[i] += 3.0 * (rhoC_r[i] * gradBD_y_r[i]- rhoC_i[i] * gradBD_y_i[i]);

                    gam_Y_i[i] += 3.0 * (rhoC_i[i] * gradBD_y_r[i] + rhoC_r[i] * gradBD_y_i[i]);

                    gam_Y_r[i] += 3.0 * (gradD_y_r[i] * rhoBC_r[i]- gradD_y_i[i] * rhoBC_i[i]);

                    gam_Y_i[i] += 3.0 * (gradD_y_i[i] * rhoBC_r[i] + gradD_y_r[i] * rhoBC_i[i]);

                    gam_Y_r[i] += 3.0 * (gradB_y_r[i] * rhoCD_r[i]- gradB_y_i[i] * rhoCD_i[i]);

                    gam_Y_i[i] += 3.0 * (gradB_y_i[i] * rhoCD_r[i] + gradB_y_r[i] * rhoCD_i[i]);

                    gam_Y_r[i] += 3.0 * (gradC_y_r[i] * rhoBD_r[i]- gradC_y_i[i] * rhoBD_i[i]);

                    gam_Y_i[i] += 3.0 * (gradC_y_i[i] * rhoBD_r[i] + gradC_y_r[i] * rhoBD_i[i]);

                    gam_Z_r[i] = 3.0 * (rhoD_r[i] * gradBC_z_r[i]- rhoD_i[i] * gradBC_z_i[i]);

                    gam_Z_i[i] = 3.0 * (rhoD_i[i] * gradBC_z_r[i] + rhoD_r[i] * gradBC_z_i[i]);

                    gam_Z_r[i] += 3.0 * (rhoB_r[i] * gradCD_z_r[i]- rhoB_i[i] * gradCD_z_i[i]);

                    gam_Z_i[i] += 3.0 * (rhoB_i[i] * gradCD_z_r[i] + rhoB_r[i] * gradCD_z_i[i]);

                    gam_Z_r[i] += 3.0 * (rhoC_r[i] * gradBD_z_r[i]- rhoC_i[i] * gradBD_z_i[i]);

                    gam_Z_i[i] += 3.0 * (rhoC_i[i] * gradBD_z_r[i] + rhoC_r[i] * gradBD_z_i[i]);

                    gam_Z_r[i] += 3.0 * (gradD_z_r[i] * rhoBC_r[i]- gradD_z_i[i] * rhoBC_i[i]);

                    gam_Z_i[i] += 3.0 * (gradD_z_i[i] * rhoBC_r[i] + gradD_z_r[i] * rhoBC_i[i]);

                    gam_Z_r[i] += 3.0 * (gradB_z_r[i] * rhoCD_r[i]- gradB_z_i[i] * rhoCD_i[i]);

                    gam_Z_i[i] += 3.0 * (gradB_z_i[i] * rhoCD_r[i] + gradB_z_r[i] * rhoCD_i[i]);

                    gam_Z_r[i] += 3.0 * (gradC_z_r[i] * rhoBD_r[i]- gradC_z_i[i] * rhoBD_i[i]);

                    gam_Z_i[i] += 3.0 * (gradC_z_i[i] * rhoBD_r[i] + gradC_z_r[i] * rhoBD_i[i]);

                    gam_XX_r[i] = 3.0 * (gradD_x_r[i] * gradBC_x_r[i]- gradD_x_i[i] * gradBC_x_i[i]);

                    gam_XX_i[i] = 3.0 * (gradD_x_i[i] * gradBC_x_r[i] + gradD_x_r[i] * gradBC_x_i[i]);

                    gam_XX_r[i] += 3.0 * (gradB_x_r[i] * gradCD_x_r[i]- gradB_x_i[i] * gradCD_x_i[i]);

                    gam_XX_i[i] += 3.0 * (gradB_x_i[i] * gradCD_x_r[i] + gradB_x_r[i] * gradCD_x_i[i]);

                    gam_XX_r[i] += 3.0 * (gradC_x_r[i] * gradBD_x_r[i]- gradC_x_i[i] * gradBD_x_i[i]);

                    gam_XX_i[i] += 3.0 * (gradC_x_i[i] * gradBD_x_r[i] + gradC_x_r[i] * gradBD_x_i[i]);

                    gam_XY_r[i] = 3.0 * (gradD_x_r[i] * gradBC_y_r[i]- gradD_x_i[i] * gradBC_y_i[i]);

                    gam_XY_i[i] = 3.0 * (gradD_x_i[i] * gradBC_y_r[i] + gradD_x_r[i] * gradBC_y_i[i]);

                    gam_XY_r[i] += 3.0 * (gradB_x_r[i] * gradCD_y_r[i]- gradB_x_i[i] * gradCD_y_i[i]);

                    gam_XY_i[i] += 3.0 * (gradB_x_i[i] * gradCD_y_r[i] + gradB_x_r[i] * gradCD_y_i[i]);

                    gam_XY_r[i] += 3.0 * (gradC_x_r[i] * gradBD_y_r[i]- gradC_x_i[i] * gradBD_y_i[i]);

                    gam_XY_i[i] += 3.0 * (gradC_x_i[i] * gradBD_y_r[i] + gradC_x_r[i] * gradBD_y_i[i]);

                    gam_XZ_r[i] = 3.0 * (gradD_x_r[i] * gradBC_z_r[i]- gradD_x_i[i] * gradBC_z_i[i]);

                    gam_XZ_i[i] = 3.0 * (gradD_x_i[i] * gradBC_z_r[i] + gradD_x_r[i] * gradBC_z_i[i]);

                    gam_XZ_r[i] += 3.0 * (gradB_x_r[i] * gradCD_z_r[i]- gradB_x_i[i] * gradCD_z_i[i]);

                    gam_XZ_i[i] += 3.0 * (gradB_x_i[i] * gradCD_z_r[i] + gradB_x_r[i] * gradCD_z_i[i]);

                    gam_XZ_r[i] += 3.0 * (gradC_x_r[i] * gradBD_z_r[i]- gradC_x_i[i] * gradBD_z_i[i]);

                    gam_XZ_i[i] += 3.0 * (gradC_x_i[i] * gradBD_z_r[i] + gradC_x_r[i] * gradBD_z_i[i]);

                    gam_YX_r[i] = 3.0 * (gradD_y_r[i] * gradBC_x_r[i]- gradD_y_i[i] * gradBC_x_i[i]);

                    gam_YX_i[i] = 3.0 * (gradD_y_i[i] * gradBC_x_r[i] + gradD_y_r[i] * gradBC_x_i[i]);

                    gam_YX_r[i] += 3.0 * (gradB_y_r[i] * gradCD_x_r[i]- gradB_y_i[i] * gradCD_x_i[i]);

                    gam_YX_i[i] += 3.0 * (gradB_y_i[i] * gradCD_x_r[i] + gradB_y_r[i] * gradCD_x_i[i]);

                    gam_YX_r[i] += 3.0 * (gradC_y_r[i] * gradBD_x_r[i]- gradC_y_i[i] * gradBD_x_i[i]);

                    gam_YX_i[i] += 3.0 * (gradC_y_i[i] * gradBD_x_r[i] + gradC_y_r[i] * gradBD_x_i[i]);

                    gam_YY_r[i] = 3.0 * (gradD_y_r[i] * gradBC_y_r[i]- gradD_y_i[i] * gradBC_y_i[i]);

                    gam_YY_i[i] = 3.0 * (gradD_y_i[i] * gradBC_y_r[i] + gradD_y_r[i] * gradBC_y_i[i]);

                    gam_YY_r[i] += 3.0 * (gradB_y_r[i] * gradCD_y_r[i]- gradB_y_i[i] * gradCD_y_i[i]);

                    gam_YY_i[i] += 3.0 * (gradB_y_i[i] * gradCD_y_r[i] + gradB_y_r[i] * gradCD_y_i[i]);

                    gam_YY_r[i] += 3.0 * (gradC_y_r[i] * gradBD_y_r[i]- gradC_y_i[i] * gradBD_y_i[i]);

                    gam_YY_i[i] += 3.0 * (gradC_y_i[i] * gradBD_y_r[i] + gradC_y_r[i] * gradBD_y_i[i]);

                    gam_YZ_r[i] = 3.0 * (gradD_y_r[i] * gradBC_z_r[i]- gradD_y_i[i] * gradBC_z_i[i]);

                    gam_YZ_i[i] = 3.0 * (gradD_y_i[i] * gradBC_z_r[i] + gradD_y_r[i] * gradBC_z_i[i]);

                    gam_YZ_r[i] += 3.0 * (gradB_y_r[i] * gradCD_z_r[i]- gradB_y_i[i] * gradCD_z_i[i]);

                    gam_YZ_i[i] += 3.0 * (gradB_y_i[i] * gradCD_z_r[i] + gradB_y_r[i] * gradCD_z_i[i]);

                    gam_YZ_r[i] += 3.0 * (gradC_y_r[i] * gradBD_z_r[i]- gradC_y_i[i] * gradBD_z_i[i]);

                    gam_YZ_i[i] += 3.0 * (gradC_y_i[i] * gradBD_z_r[i] + gradC_y_r[i] * gradBD_z_i[i]);

                    gam_ZX_r[i] = 3.0 * (gradD_z_r[i] * gradBC_x_r[i]- gradD_z_i[i] * gradBC_x_i[i]);

                    gam_ZX_i[i] = 3.0 * (gradD_z_i[i] * gradBC_x_r[i] + gradD_z_r[i] * gradBC_x_i[i]);

                    gam_ZX_r[i] += 3.0 * (gradB_z_r[i] * gradCD_x_r[i]- gradB_z_i[i] * gradCD_x_i[i]);

                    gam_ZX_i[i] += 3.0 * (gradB_z_i[i] * gradCD_x_r[i] + gradB_z_r[i] * gradCD_x_i[i]);

                    gam_ZX_r[i] += 3.0 * (gradC_z_r[i] * gradBD_x_r[i]- gradC_z_i[i] * gradBD_x_i[i]);

                    gam_ZX_i[i] += 3.0 * (gradC_z_i[i] * gradBD_x_r[i] + gradC_z_r[i] * gradBD_x_i[i]);

                    gam_ZY_r[i] = 3.0 * (gradD_z_r[i] * gradBC_y_r[i]- gradD_z_i[i] * gradBC_y_i[i]);

                    gam_ZY_i[i] = 3.0 * (gradD_z_i[i] * gradBC_y_r[i] + gradD_z_r[i] * gradBC_y_i[i]);

                    gam_ZY_r[i] += 3.0 * (gradB_z_r[i] * gradCD_y_r[i]- gradB_z_i[i] * gradCD_y_i[i]);

                    gam_ZY_i[i] += 3.0 * (gradB_z_i[i] * gradCD_y_r[i] + gradB_z_r[i] * gradCD_y_i[i]);

                    gam_ZY_r[i] += 3.0 * (gradC_z_r[i] * gradBD_y_r[i]- gradC_z_i[i] * gradBD_y_i[i]);

                    gam_ZY_i[i] += 3.0 * (gradC_z_i[i] * gradBD_y_r[i] + gradC_z_r[i] * gradBD_y_i[i]);

                    gam_ZZ_r[i] = 3.0 * (gradD_z_r[i] * gradBC_z_r[i]- gradD_z_i[i] * gradBC_z_i[i]);

                    gam_ZZ_i[i] = 3.0 * (gradD_z_i[i] * gradBC_z_r[i] + gradD_z_r[i] * gradBC_z_i[i]);

                    gam_ZZ_r[i] += 3.0 * (gradB_z_r[i] * gradCD_z_r[i]- gradB_z_i[i] * gradCD_z_i[i]);

                    gam_ZZ_i[i] += 3.0 * (gradB_z_i[i] * gradCD_z_r[i] + gradB_z_r[i] * gradCD_z_i[i]);

                    gam_ZZ_r[i] += 3.0 * (gradC_z_r[i] * gradBD_z_r[i]- gradC_z_i[i] * gradBD_z_i[i]);

                    gam_ZZ_i[i] += 3.0 * (gradC_z_i[i] * gradBD_z_r[i] + gradC_z_r[i] * gradBD_z_i[i]);

                }
            }
        }
    }
}


std::ostream&
operator<<(std::ostream& output, const CDensityGridCubic& source)
{
    output << std::endl;

    output << "[CDensityGridCubic (Object):" << &source << "]" << std::endl;

    output << "_gridType: " << to_string(source._gridType) << std::endl;

    output << "_densityValues: " << std::endl;

    output << source._densityValues << std::endl;

    return output;
}
