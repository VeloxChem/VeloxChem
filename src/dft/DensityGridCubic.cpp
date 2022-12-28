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

    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 3 : 3;

    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 66 : 66;

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
CDensityGridCubic::gam2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::gam2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::piZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::piZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);
    return nullptr;
}

double*
CDensityGridCubic::piZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::piXXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piXZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piXZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piYZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piYZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::piZZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::piZZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}


const double*
CDensityGridCubic::gamXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::gamX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::gamZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(53 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gamZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(53 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2X(const int32_t iDensityMatrix) const
{
if (_gridType == dengrid::ab) return _densityValues.data(54 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}
double*
CDensityGridCubic::gam2X(const int32_t iDensityMatrix)
{
if (_gridType == dengrid::ab) return _densityValues.data(54 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}

const double*
CDensityGridCubic::gam2Y(const int32_t iDensityMatrix) const
{
if (_gridType == dengrid::ab) return _densityValues.data(55 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}
double*
CDensityGridCubic::gam2Y(const int32_t iDensityMatrix)
{
if (_gridType == dengrid::ab) return _densityValues.data(55 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}

const double*
CDensityGridCubic::gam2Z(const int32_t iDensityMatrix) const
{
if (_gridType == dengrid::ab) return _densityValues.data(56 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}
double*
CDensityGridCubic::gam2Z(const int32_t iDensityMatrix)
{
if (_gridType == dengrid::ab) return _densityValues.data(56 * _nDensityMatrices + iDensityMatrix);

return nullptr;
}

const double*
CDensityGridCubic::gam2XX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(57 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2XX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(57 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2XY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(58 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2XY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(58 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2XZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(59 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2XZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(59 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2YX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(60 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2YX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(60 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2YY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(61 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2YY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(61 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2YZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(62 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2YZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(62 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2ZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(63 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2ZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(63 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2ZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(64 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2ZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(64 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::gam2ZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(65 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::gam2ZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(65 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double
CDensityGridCubic::prod3_r(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i)
{
    double BCD  =  (B_r * C_r * D_r 
                  - B_i * D_i * C_r 
                  - C_i * D_i * B_r 
                  - B_i * C_i * D_r);
    return BCD;
}

double
CDensityGridCubic::prod3_i(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i)
{
    double BCD  =  (- B_i * C_i * D_i 
                    + B_i * C_r * D_r 
                    + C_i * B_r * D_r 
                    + D_i * B_r * C_r);
    return BCD;
}


double
CDensityGridCubic::prod2_r(double B_r, double B_i, double C_r, double C_i)
{
    double BC  =  (B_r * C_r- B_i * C_i);

    return BC;
}

double
CDensityGridCubic::prod2_i(double B_r, double B_i, double C_r, double C_i)
{
    double BC  =  (B_i * C_r + B_r * C_i);

    return BC;
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
             for (int32_t j = 0; j < numdens / 30; j++)
            {
                
                // Second-order 

                auto gam_sig_xx_r = gam2(24 * j);
                auto gam_sig_xx_i = gam2(24 * j + 1);
                auto gam_sig_yy_r = gam2(24 * j + 2);
                auto gam_sig_yy_i = gam2(24 * j + 3);
                auto gam_sig_zz_r = gam2(24 * j + 4);
                auto gam_sig_zz_i = gam2(24 * j + 5);
                auto gam_sig_xy_r = gam2(24 * j + 6);
                auto gam_sig_xy_i = gam2(24 * j + 7);
                auto gam_sig_xz_r = gam2(24 * j + 8);
                auto gam_sig_xz_i = gam2(24 * j + 9);
                auto gam_sig_yz_r = gam2(24 * j + 10);
                auto gam_sig_yz_i = gam2(24 * j + 11);
                auto gam_lamtau_xx_r = gam2(24 * j + 12);
                auto gam_lamtau_xx_i = gam2(24 * j + 13);
                auto gam_lamtau_yy_r = gam2(24 * j + 14);
                auto gam_lamtau_yy_i = gam2(24 * j + 15);
                auto gam_lamtau_zz_r = gam2(24 * j + 16);
                auto gam_lamtau_zz_i = gam2(24 * j + 17);
                auto gam_lamtau_xy_r = gam2(24 * j + 18);
                auto gam_lamtau_xy_i = gam2(24 * j + 19);
                auto gam_lamtau_xz_r = gam2(24 * j + 20);
                auto gam_lamtau_xz_i = gam2(24 * j + 21);
                auto gam_lamtau_yz_r = gam2(24 * j + 22);
                auto gam_lamtau_yz_i = gam2(24 * j + 23);

                // Third-order 
                auto gamx_r = gam(6 * j);
                auto gamx_i = gam(6 * j + 1);
                auto gamy_r = gam(6 * j + 2);
                auto gamy_i = gam(6 * j + 3);
                auto gamz_r = gam(6 * j + 4);
                auto gamz_i = gam(6 * j + 5);
                auto pix_r = pi(6 * j);
                auto pix_i = pi(6 * j + 1);
                auto piy_r = pi(6 * j + 2);
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
                auto rho_sig_xx_r = rw2DensityGrid.alphaDensity(24 * j);
                auto rho_sig_xx_i = rw2DensityGrid.alphaDensity(24 * j + 1);
                auto rho_sig_yy_r = rw2DensityGrid.alphaDensity(24 * j + 2);
                auto rho_sig_yy_i = rw2DensityGrid.alphaDensity(24 * j + 3);
                auto rho_sig_zz_r = rw2DensityGrid.alphaDensity(24 * j + 4);
                auto rho_sig_zz_i = rw2DensityGrid.alphaDensity(24 * j + 5);
                auto rho_sig_xy_r = rw2DensityGrid.alphaDensity(24 * j + 6);
                auto rho_sig_xy_i = rw2DensityGrid.alphaDensity(24 * j + 7);
                auto rho_sig_xz_r = rw2DensityGrid.alphaDensity(24 * j + 8);
                auto rho_sig_xz_i = rw2DensityGrid.alphaDensity(24 * j + 9);
                auto rho_sig_yz_r = rw2DensityGrid.alphaDensity(24 * j + 10);
                auto rho_sig_yz_i = rw2DensityGrid.alphaDensity(24 * j + 11);
                auto rho_lamtau_xx_r = rw2DensityGrid.alphaDensity(24 * j + 12);
                auto rho_lamtau_xx_i = rw2DensityGrid.alphaDensity(24 * j + 13);
                auto rho_lamtau_yy_r = rw2DensityGrid.alphaDensity(24 * j + 14);
                auto rho_lamtau_yy_i = rw2DensityGrid.alphaDensity(24 * j + 15);
                auto rho_lamtau_zz_r = rw2DensityGrid.alphaDensity(24 * j + 16);
                auto rho_lamtau_zz_i = rw2DensityGrid.alphaDensity(24 * j + 17);
                auto rho_lamtau_xy_r = rw2DensityGrid.alphaDensity(24 * j + 18);
                auto rho_lamtau_xy_i = rw2DensityGrid.alphaDensity(24 * j + 19);
                auto rho_lamtau_xz_r = rw2DensityGrid.alphaDensity(24 * j + 20);
                auto rho_lamtau_xz_i = rw2DensityGrid.alphaDensity(24 * j + 21);
                auto rho_lamtau_yz_r = rw2DensityGrid.alphaDensity(24 * j + 22);
                auto rho_lamtau_yz_i = rw2DensityGrid.alphaDensity(24 * j + 23);

                for (int32_t i = 0; i < npoints; i++)
                {

                  

                    double sig_term_r =  + 6.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);

                    double sig_term_i =  + 6.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_r =  + 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])
                                            + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])
                                            + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]);
                    double lamtau_term_i =  + 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])
                                            + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])
                                            + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]);
                    gam_sig_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_r;
                    gam_sig_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_r;
                    gam_sig_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_r;
                    gam_sig_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_sig_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_sig_xx_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_i;
                    gam_sig_yy_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_i;
                    gam_sig_zz_i[i] = 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_i;
                    gam_sig_xy_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_xz_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_sig_yz_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i]) 
                                       + 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])  + lamtau_term_r;
                    gam_lamtau_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])  + lamtau_term_r;
                    gam_lamtau_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i])  + lamtau_term_r;
                    gam_lamtau_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCy_r[i],rhoCy_i[i]);
                    gam_lamtau_xx_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i]) 
                                       + 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])  + lamtau_term_i;
                    gam_lamtau_yy_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])  + lamtau_term_i;
                    gam_lamtau_zz_i[i] = 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i])  + lamtau_term_i;
                    gam_lamtau_xy_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_xz_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_yz_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCy_r[i],rhoCy_i[i]);


                    gamx_r[i] = prod2_r(rho_sig_xx_r[i],rho_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gamy_r[i] = prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_yy_r[i],rho_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gamz_r[i] = prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_zz_r[i],rho_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gamx_i[i] = prod2_i(rho_sig_xx_r[i],rho_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gamy_i[i] = prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_yy_r[i],rho_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gamz_i[i] = prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_zz_r[i],rho_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);


                    pix_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_r[i],gam_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    piy_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yy_r[i],gam_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    piz_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_zz_r[i],gam_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    pix_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_r[i],gam_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    piy_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yy_r[i],gam_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);
                    piz_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_zz_r[i],gam_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);


                }

            }

        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 8; j++)
            {
                
                auto gam_bc_r = gam2(6 * j);

                auto gam_bc_i = gam2(6 * j + 1);

                auto gam_bd_r = gam2(6 * j + 2);

                auto gam_bd_i = gam2(6 * j + 3);

                auto gam_cd_r = gam2(6 * j + 4);

                auto gam_cd_i = gam2(6 * j + 5);


                auto gam_bcd_r = gam(2 * j);

                auto gam_bcd_i = gam(2 * j + 1);

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

                for (int32_t i = 0; i < npoints; i++)
                {

                    gam_bc_r[i] = 2.0 * prod2_r(rhoB_r[i],rhoB_i[i], rhoC_r[i],rhoC_i[i]);
                    gam_bc_i[i] = 2.0 * prod2_i(rhoB_r[i],rhoB_i[i], rhoC_r[i],rhoC_i[i]);

                    gam_bd_r[i] = 2.0 * prod2_r(rhoB_r[i],rhoB_i[i], rhoD_r[i],rhoD_i[i]);
                    gam_bd_i[i] = 2.0 * prod2_i(rhoB_r[i],rhoB_i[i], rhoD_r[i],rhoD_i[i]);

                    gam_cd_r[i] = 2.0 * prod2_r(rhoC_r[i],rhoC_i[i], rhoD_r[i],rhoD_i[i]);
                    gam_cd_i[i] = 2.0 * prod2_i(rhoC_r[i],rhoC_i[i], rhoD_r[i],rhoD_i[i]);

                    // D BC 
                    gam_bcd_r[i] =  3.0 * prod2_r(rhoD_r[i],rhoD_i[i], rhoBC_r[i],rhoBC_i[i]);
                    gam_bcd_i[i] =  3.0 * prod2_i(rhoD_r[i],rhoD_i[i], rhoBC_r[i],rhoBC_i[i]);

                    // B CD
                    gam_bcd_r[i] += 3.0 * prod2_r(rhoB_r[i],rhoB_i[i], rhoCD_r[i],rhoCD_i[i]);
                    gam_bcd_i[i] += 3.0 * prod2_i(rhoB_r[i],rhoB_i[i], rhoCD_r[i],rhoCD_i[i]);

                    // C BD
                    gam_bcd_r[i] += 3.0 * prod2_r(rhoC_r[i],rhoC_i[i], rhoBD_r[i],rhoBD_i[i]);
                    gam_bcd_i[i] += 3.0 * prod2_i(rhoC_r[i],rhoC_i[i], rhoBD_r[i],rhoBD_i[i]);

                    pi_r[i] = prod2_r(rhoB_r[i],rhoB_i[i], gam_cd_r[i],gam_cd_i[i])
                            + prod2_r(rhoC_r[i],rhoC_i[i], gam_bd_r[i],gam_bd_i[i])
                            + prod2_r(rhoD_r[i],rhoD_i[i], gam_bc_r[i],gam_bc_i[i]);

                    pi_i[i] = prod2_i(rhoB_r[i],rhoB_i[i], gam_cd_r[i],gam_cd_i[i])
                            + prod2_i(rhoC_r[i],rhoC_i[i], gam_bd_r[i],gam_bd_i[i])
                            + prod2_i(rhoD_r[i],rhoD_i[i], gam_bc_r[i],gam_bc_i[i]);

                }
            }
        }
    }
    if (xcFuncType == xcfun::gga)
    {
        
        if (fstr::upcase(CubeMode) == "TPA")
        {
            for (int32_t j = 0; j < numdens / 30; j++)
            {

                auto gam_sig_xx_r = gam2(24 * j);
                auto gam_sig_xx_i = gam2(24 * j + 1);
                auto gam_sig_yy_r = gam2(24 * j + 2);
                auto gam_sig_yy_i = gam2(24 * j + 3);
                auto gam_sig_zz_r = gam2(24 * j + 4);
                auto gam_sig_zz_i = gam2(24 * j + 5);
                auto gam_sig_xy_r = gam2(24 * j + 6);
                auto gam_sig_xy_i = gam2(24 * j + 7);
                auto gam_sig_xz_r = gam2(24 * j + 8);
                auto gam_sig_xz_i = gam2(24 * j + 9);
                auto gam_sig_yz_r = gam2(24 * j + 10);
                auto gam_sig_yz_i = gam2(24 * j + 11);
                auto gam_lamtau_xx_r = gam2(24 * j + 12);
                auto gam_lamtau_xx_i = gam2(24 * j + 13);
                auto gam_lamtau_yy_r = gam2(24 * j + 14);
                auto gam_lamtau_yy_i = gam2(24 * j + 15);
                auto gam_lamtau_zz_r = gam2(24 * j + 16);
                auto gam_lamtau_zz_i = gam2(24 * j + 17);
                auto gam_lamtau_xy_r = gam2(24 * j + 18);
                auto gam_lamtau_xy_i = gam2(24 * j + 19);
                auto gam_lamtau_xz_r = gam2(24 * j + 20);
                auto gam_lamtau_xz_i = gam2(24 * j + 21);
                auto gam_lamtau_yz_r = gam2(24 * j + 22);
                auto gam_lamtau_yz_i = gam2(24 * j + 23);
                
                auto gam_sig_xx_x_r = gam2X(24 * j);
                auto gam_sig_xx_x_i = gam2X(24 * j + 1);
                auto gam_sig_yy_x_r = gam2X(24 * j + 2);
                auto gam_sig_yy_x_i = gam2X(24 * j + 3);
                auto gam_sig_zz_x_r = gam2X(24 * j + 4);
                auto gam_sig_zz_x_i = gam2X(24 * j + 5);
                auto gam_sig_xy_x_r = gam2X(24 * j + 6);
                auto gam_sig_xy_x_i = gam2X(24 * j + 7);
                auto gam_sig_xz_x_r = gam2X(24 * j + 8);
                auto gam_sig_xz_x_i = gam2X(24 * j + 9);
                auto gam_sig_yz_x_r = gam2X(24 * j + 10);
                auto gam_sig_yz_x_i = gam2X(24 * j + 11);
                auto gam_lamtau_xx_x_r = gam2X(24 * j + 12);
                auto gam_lamtau_xx_x_i = gam2X(24 * j + 13);
                auto gam_lamtau_yy_x_r = gam2X(24 * j + 14);
                auto gam_lamtau_yy_x_i = gam2X(24 * j + 15);
                auto gam_lamtau_zz_x_r = gam2X(24 * j + 16);
                auto gam_lamtau_zz_x_i = gam2X(24 * j + 17);
                auto gam_lamtau_xy_x_r = gam2X(24 * j + 18);
                auto gam_lamtau_xy_x_i = gam2X(24 * j + 19);
                auto gam_lamtau_xz_x_r = gam2X(24 * j + 20);
                auto gam_lamtau_xz_x_i = gam2X(24 * j + 21);
                auto gam_lamtau_yz_x_r = gam2X(24 * j + 22);
                auto gam_lamtau_yz_x_i = gam2X(24 * j + 23);
                auto gam_sig_xx_y_r = gam2Y(24 * j);
                auto gam_sig_xx_y_i = gam2Y(24 * j + 1);
                auto gam_sig_yy_y_r = gam2Y(24 * j + 2);
                auto gam_sig_yy_y_i = gam2Y(24 * j + 3);
                auto gam_sig_zz_y_r = gam2Y(24 * j + 4);
                auto gam_sig_zz_y_i = gam2Y(24 * j + 5);
                auto gam_sig_xy_y_r = gam2Y(24 * j + 6);
                auto gam_sig_xy_y_i = gam2Y(24 * j + 7);
                auto gam_sig_xz_y_r = gam2Y(24 * j + 8);
                auto gam_sig_xz_y_i = gam2Y(24 * j + 9);
                auto gam_sig_yz_y_r = gam2Y(24 * j + 10);
                auto gam_sig_yz_y_i = gam2Y(24 * j + 11);
                auto gam_lamtau_xx_y_r = gam2Y(24 * j + 12);
                auto gam_lamtau_xx_y_i = gam2Y(24 * j + 13);
                auto gam_lamtau_yy_y_r = gam2Y(24 * j + 14);
                auto gam_lamtau_yy_y_i = gam2Y(24 * j + 15);
                auto gam_lamtau_zz_y_r = gam2Y(24 * j + 16);
                auto gam_lamtau_zz_y_i = gam2Y(24 * j + 17);
                auto gam_lamtau_xy_y_r = gam2Y(24 * j + 18);
                auto gam_lamtau_xy_y_i = gam2Y(24 * j + 19);
                auto gam_lamtau_xz_y_r = gam2Y(24 * j + 20);
                auto gam_lamtau_xz_y_i = gam2Y(24 * j + 21);
                auto gam_lamtau_yz_y_r = gam2Y(24 * j + 22);
                auto gam_lamtau_yz_y_i = gam2Y(24 * j + 23);
                auto gam_sig_xx_z_r = gam2Z(24 * j);
                auto gam_sig_xx_z_i = gam2Z(24 * j + 1);
                auto gam_sig_yy_z_r = gam2Z(24 * j + 2);
                auto gam_sig_yy_z_i = gam2Z(24 * j + 3);
                auto gam_sig_zz_z_r = gam2Z(24 * j + 4);
                auto gam_sig_zz_z_i = gam2Z(24 * j + 5);
                auto gam_sig_xy_z_r = gam2Z(24 * j + 6);
                auto gam_sig_xy_z_i = gam2Z(24 * j + 7);
                auto gam_sig_xz_z_r = gam2Z(24 * j + 8);
                auto gam_sig_xz_z_i = gam2Z(24 * j + 9);
                auto gam_sig_yz_z_r = gam2Z(24 * j + 10);
                auto gam_sig_yz_z_i = gam2Z(24 * j + 11);
                auto gam_lamtau_xx_z_r = gam2Z(24 * j + 12);
                auto gam_lamtau_xx_z_i = gam2Z(24 * j + 13);
                auto gam_lamtau_yy_z_r = gam2Z(24 * j + 14);
                auto gam_lamtau_yy_z_i = gam2Z(24 * j + 15);
                auto gam_lamtau_zz_z_r = gam2Z(24 * j + 16);
                auto gam_lamtau_zz_z_i = gam2Z(24 * j + 17);
                auto gam_lamtau_xy_z_r = gam2Z(24 * j + 18);
                auto gam_lamtau_xy_z_i = gam2Z(24 * j + 19);
                auto gam_lamtau_xz_z_r = gam2Z(24 * j + 20);
                auto gam_lamtau_xz_z_i = gam2Z(24 * j + 21);
                auto gam_lamtau_yz_z_r = gam2Z(24 * j + 22);
                auto gam_lamtau_yz_z_i = gam2Z(24 * j + 23);
                auto gam_sig_xx_xx_r = gam2XX(24 * j);
                auto gam_sig_xx_xx_i = gam2XX(24 * j + 1);
                auto gam_sig_yy_xx_r = gam2XX(24 * j + 2);
                auto gam_sig_yy_xx_i = gam2XX(24 * j + 3);
                auto gam_sig_zz_xx_r = gam2XX(24 * j + 4);
                auto gam_sig_zz_xx_i = gam2XX(24 * j + 5);
                auto gam_sig_xy_xx_r = gam2XX(24 * j + 6);
                auto gam_sig_xy_xx_i = gam2XX(24 * j + 7);
                auto gam_sig_xz_xx_r = gam2XX(24 * j + 8);
                auto gam_sig_xz_xx_i = gam2XX(24 * j + 9);
                auto gam_sig_yz_xx_r = gam2XX(24 * j + 10);
                auto gam_sig_yz_xx_i = gam2XX(24 * j + 11);
                auto gam_lamtau_xx_xx_r = gam2XX(24 * j + 12);
                auto gam_lamtau_xx_xx_i = gam2XX(24 * j + 13);
                auto gam_lamtau_yy_xx_r = gam2XX(24 * j + 14);
                auto gam_lamtau_yy_xx_i = gam2XX(24 * j + 15);
                auto gam_lamtau_zz_xx_r = gam2XX(24 * j + 16);
                auto gam_lamtau_zz_xx_i = gam2XX(24 * j + 17);
                auto gam_lamtau_xy_xx_r = gam2XX(24 * j + 18);
                auto gam_lamtau_xy_xx_i = gam2XX(24 * j + 19);
                auto gam_lamtau_xz_xx_r = gam2XX(24 * j + 20);
                auto gam_lamtau_xz_xx_i = gam2XX(24 * j + 21);
                auto gam_lamtau_yz_xx_r = gam2XX(24 * j + 22);
                auto gam_lamtau_yz_xx_i = gam2XX(24 * j + 23);
                auto gam_sig_xx_xy_r = gam2XY(24 * j);
                auto gam_sig_xx_xy_i = gam2XY(24 * j + 1);
                auto gam_sig_yy_xy_r = gam2XY(24 * j + 2);
                auto gam_sig_yy_xy_i = gam2XY(24 * j + 3);
                auto gam_sig_zz_xy_r = gam2XY(24 * j + 4);
                auto gam_sig_zz_xy_i = gam2XY(24 * j + 5);
                auto gam_sig_xy_xy_r = gam2XY(24 * j + 6);
                auto gam_sig_xy_xy_i = gam2XY(24 * j + 7);
                auto gam_sig_xz_xy_r = gam2XY(24 * j + 8);
                auto gam_sig_xz_xy_i = gam2XY(24 * j + 9);
                auto gam_sig_yz_xy_r = gam2XY(24 * j + 10);
                auto gam_sig_yz_xy_i = gam2XY(24 * j + 11);
                auto gam_lamtau_xx_xy_r = gam2XY(24 * j + 12);
                auto gam_lamtau_xx_xy_i = gam2XY(24 * j + 13);
                auto gam_lamtau_yy_xy_r = gam2XY(24 * j + 14);
                auto gam_lamtau_yy_xy_i = gam2XY(24 * j + 15);
                auto gam_lamtau_zz_xy_r = gam2XY(24 * j + 16);
                auto gam_lamtau_zz_xy_i = gam2XY(24 * j + 17);
                auto gam_lamtau_xy_xy_r = gam2XY(24 * j + 18);
                auto gam_lamtau_xy_xy_i = gam2XY(24 * j + 19);
                auto gam_lamtau_xz_xy_r = gam2XY(24 * j + 20);
                auto gam_lamtau_xz_xy_i = gam2XY(24 * j + 21);
                auto gam_lamtau_yz_xy_r = gam2XY(24 * j + 22);
                auto gam_lamtau_yz_xy_i = gam2XY(24 * j + 23);
                auto gam_sig_xx_xz_r = gam2XZ(24 * j);
                auto gam_sig_xx_xz_i = gam2XZ(24 * j + 1);
                auto gam_sig_yy_xz_r = gam2XZ(24 * j + 2);
                auto gam_sig_yy_xz_i = gam2XZ(24 * j + 3);
                auto gam_sig_zz_xz_r = gam2XZ(24 * j + 4);
                auto gam_sig_zz_xz_i = gam2XZ(24 * j + 5);
                auto gam_sig_xy_xz_r = gam2XZ(24 * j + 6);
                auto gam_sig_xy_xz_i = gam2XZ(24 * j + 7);
                auto gam_sig_xz_xz_r = gam2XZ(24 * j + 8);
                auto gam_sig_xz_xz_i = gam2XZ(24 * j + 9);
                auto gam_sig_yz_xz_r = gam2XZ(24 * j + 10);
                auto gam_sig_yz_xz_i = gam2XZ(24 * j + 11);
                auto gam_lamtau_xx_xz_r = gam2XZ(24 * j + 12);
                auto gam_lamtau_xx_xz_i = gam2XZ(24 * j + 13);
                auto gam_lamtau_yy_xz_r = gam2XZ(24 * j + 14);
                auto gam_lamtau_yy_xz_i = gam2XZ(24 * j + 15);
                auto gam_lamtau_zz_xz_r = gam2XZ(24 * j + 16);
                auto gam_lamtau_zz_xz_i = gam2XZ(24 * j + 17);
                auto gam_lamtau_xy_xz_r = gam2XZ(24 * j + 18);
                auto gam_lamtau_xy_xz_i = gam2XZ(24 * j + 19);
                auto gam_lamtau_xz_xz_r = gam2XZ(24 * j + 20);
                auto gam_lamtau_xz_xz_i = gam2XZ(24 * j + 21);
                auto gam_lamtau_yz_xz_r = gam2XZ(24 * j + 22);
                auto gam_lamtau_yz_xz_i = gam2XZ(24 * j + 23);
                auto gam_sig_xx_yx_r = gam2YX(24 * j);
                auto gam_sig_xx_yx_i = gam2YX(24 * j + 1);
                auto gam_sig_yy_yx_r = gam2YX(24 * j + 2);
                auto gam_sig_yy_yx_i = gam2YX(24 * j + 3);
                auto gam_sig_zz_yx_r = gam2YX(24 * j + 4);
                auto gam_sig_zz_yx_i = gam2YX(24 * j + 5);
                auto gam_sig_xy_yx_r = gam2YX(24 * j + 6);
                auto gam_sig_xy_yx_i = gam2YX(24 * j + 7);
                auto gam_sig_xz_yx_r = gam2YX(24 * j + 8);
                auto gam_sig_xz_yx_i = gam2YX(24 * j + 9);
                auto gam_sig_yz_yx_r = gam2YX(24 * j + 10);
                auto gam_sig_yz_yx_i = gam2YX(24 * j + 11);
                auto gam_lamtau_xx_yx_r = gam2YX(24 * j + 12);
                auto gam_lamtau_xx_yx_i = gam2YX(24 * j + 13);
                auto gam_lamtau_yy_yx_r = gam2YX(24 * j + 14);
                auto gam_lamtau_yy_yx_i = gam2YX(24 * j + 15);
                auto gam_lamtau_zz_yx_r = gam2YX(24 * j + 16);
                auto gam_lamtau_zz_yx_i = gam2YX(24 * j + 17);
                auto gam_lamtau_xy_yx_r = gam2YX(24 * j + 18);
                auto gam_lamtau_xy_yx_i = gam2YX(24 * j + 19);
                auto gam_lamtau_xz_yx_r = gam2YX(24 * j + 20);
                auto gam_lamtau_xz_yx_i = gam2YX(24 * j + 21);
                auto gam_lamtau_yz_yx_r = gam2YX(24 * j + 22);
                auto gam_lamtau_yz_yx_i = gam2YX(24 * j + 23);
                auto gam_sig_xx_yy_r = gam2YY(24 * j);
                auto gam_sig_xx_yy_i = gam2YY(24 * j + 1);
                auto gam_sig_yy_yy_r = gam2YY(24 * j + 2);
                auto gam_sig_yy_yy_i = gam2YY(24 * j + 3);
                auto gam_sig_zz_yy_r = gam2YY(24 * j + 4);
                auto gam_sig_zz_yy_i = gam2YY(24 * j + 5);
                auto gam_sig_xy_yy_r = gam2YY(24 * j + 6);
                auto gam_sig_xy_yy_i = gam2YY(24 * j + 7);
                auto gam_sig_xz_yy_r = gam2YY(24 * j + 8);
                auto gam_sig_xz_yy_i = gam2YY(24 * j + 9);
                auto gam_sig_yz_yy_r = gam2YY(24 * j + 10);
                auto gam_sig_yz_yy_i = gam2YY(24 * j + 11);
                auto gam_lamtau_xx_yy_r = gam2YY(24 * j + 12);
                auto gam_lamtau_xx_yy_i = gam2YY(24 * j + 13);
                auto gam_lamtau_yy_yy_r = gam2YY(24 * j + 14);
                auto gam_lamtau_yy_yy_i = gam2YY(24 * j + 15);
                auto gam_lamtau_zz_yy_r = gam2YY(24 * j + 16);
                auto gam_lamtau_zz_yy_i = gam2YY(24 * j + 17);
                auto gam_lamtau_xy_yy_r = gam2YY(24 * j + 18);
                auto gam_lamtau_xy_yy_i = gam2YY(24 * j + 19);
                auto gam_lamtau_xz_yy_r = gam2YY(24 * j + 20);
                auto gam_lamtau_xz_yy_i = gam2YY(24 * j + 21);
                auto gam_lamtau_yz_yy_r = gam2YY(24 * j + 22);
                auto gam_lamtau_yz_yy_i = gam2YY(24 * j + 23);
                auto gam_sig_xx_yz_r = gam2YZ(24 * j);
                auto gam_sig_xx_yz_i = gam2YZ(24 * j + 1);
                auto gam_sig_yy_yz_r = gam2YZ(24 * j + 2);
                auto gam_sig_yy_yz_i = gam2YZ(24 * j + 3);
                auto gam_sig_zz_yz_r = gam2YZ(24 * j + 4);
                auto gam_sig_zz_yz_i = gam2YZ(24 * j + 5);
                auto gam_sig_xy_yz_r = gam2YZ(24 * j + 6);
                auto gam_sig_xy_yz_i = gam2YZ(24 * j + 7);
                auto gam_sig_xz_yz_r = gam2YZ(24 * j + 8);
                auto gam_sig_xz_yz_i = gam2YZ(24 * j + 9);
                auto gam_sig_yz_yz_r = gam2YZ(24 * j + 10);
                auto gam_sig_yz_yz_i = gam2YZ(24 * j + 11);
                auto gam_lamtau_xx_yz_r = gam2YZ(24 * j + 12);
                auto gam_lamtau_xx_yz_i = gam2YZ(24 * j + 13);
                auto gam_lamtau_yy_yz_r = gam2YZ(24 * j + 14);
                auto gam_lamtau_yy_yz_i = gam2YZ(24 * j + 15);
                auto gam_lamtau_zz_yz_r = gam2YZ(24 * j + 16);
                auto gam_lamtau_zz_yz_i = gam2YZ(24 * j + 17);
                auto gam_lamtau_xy_yz_r = gam2YZ(24 * j + 18);
                auto gam_lamtau_xy_yz_i = gam2YZ(24 * j + 19);
                auto gam_lamtau_xz_yz_r = gam2YZ(24 * j + 20);
                auto gam_lamtau_xz_yz_i = gam2YZ(24 * j + 21);
                auto gam_lamtau_yz_yz_r = gam2YZ(24 * j + 22);
                auto gam_lamtau_yz_yz_i = gam2YZ(24 * j + 23);
                auto gam_sig_xx_zx_r = gam2ZX(24 * j);
                auto gam_sig_xx_zx_i = gam2ZX(24 * j + 1);
                auto gam_sig_yy_zx_r = gam2ZX(24 * j + 2);
                auto gam_sig_yy_zx_i = gam2ZX(24 * j + 3);
                auto gam_sig_zz_zx_r = gam2ZX(24 * j + 4);
                auto gam_sig_zz_zx_i = gam2ZX(24 * j + 5);
                auto gam_sig_xy_zx_r = gam2ZX(24 * j + 6);
                auto gam_sig_xy_zx_i = gam2ZX(24 * j + 7);
                auto gam_sig_xz_zx_r = gam2ZX(24 * j + 8);
                auto gam_sig_xz_zx_i = gam2ZX(24 * j + 9);
                auto gam_sig_yz_zx_r = gam2ZX(24 * j + 10);
                auto gam_sig_yz_zx_i = gam2ZX(24 * j + 11);
                auto gam_lamtau_xx_zx_r = gam2ZX(24 * j + 12);
                auto gam_lamtau_xx_zx_i = gam2ZX(24 * j + 13);
                auto gam_lamtau_yy_zx_r = gam2ZX(24 * j + 14);
                auto gam_lamtau_yy_zx_i = gam2ZX(24 * j + 15);
                auto gam_lamtau_zz_zx_r = gam2ZX(24 * j + 16);
                auto gam_lamtau_zz_zx_i = gam2ZX(24 * j + 17);
                auto gam_lamtau_xy_zx_r = gam2ZX(24 * j + 18);
                auto gam_lamtau_xy_zx_i = gam2ZX(24 * j + 19);
                auto gam_lamtau_xz_zx_r = gam2ZX(24 * j + 20);
                auto gam_lamtau_xz_zx_i = gam2ZX(24 * j + 21);
                auto gam_lamtau_yz_zx_r = gam2ZX(24 * j + 22);
                auto gam_lamtau_yz_zx_i = gam2ZX(24 * j + 23);
                auto gam_sig_xx_zy_r = gam2ZY(24 * j);
                auto gam_sig_xx_zy_i = gam2ZY(24 * j + 1);
                auto gam_sig_yy_zy_r = gam2ZY(24 * j + 2);
                auto gam_sig_yy_zy_i = gam2ZY(24 * j + 3);
                auto gam_sig_zz_zy_r = gam2ZY(24 * j + 4);
                auto gam_sig_zz_zy_i = gam2ZY(24 * j + 5);
                auto gam_sig_xy_zy_r = gam2ZY(24 * j + 6);
                auto gam_sig_xy_zy_i = gam2ZY(24 * j + 7);
                auto gam_sig_xz_zy_r = gam2ZY(24 * j + 8);
                auto gam_sig_xz_zy_i = gam2ZY(24 * j + 9);
                auto gam_sig_yz_zy_r = gam2ZY(24 * j + 10);
                auto gam_sig_yz_zy_i = gam2ZY(24 * j + 11);
                auto gam_lamtau_xx_zy_r = gam2ZY(24 * j + 12);
                auto gam_lamtau_xx_zy_i = gam2ZY(24 * j + 13);
                auto gam_lamtau_yy_zy_r = gam2ZY(24 * j + 14);
                auto gam_lamtau_yy_zy_i = gam2ZY(24 * j + 15);
                auto gam_lamtau_zz_zy_r = gam2ZY(24 * j + 16);
                auto gam_lamtau_zz_zy_i = gam2ZY(24 * j + 17);
                auto gam_lamtau_xy_zy_r = gam2ZY(24 * j + 18);
                auto gam_lamtau_xy_zy_i = gam2ZY(24 * j + 19);
                auto gam_lamtau_xz_zy_r = gam2ZY(24 * j + 20);
                auto gam_lamtau_xz_zy_i = gam2ZY(24 * j + 21);
                auto gam_lamtau_yz_zy_r = gam2ZY(24 * j + 22);
                auto gam_lamtau_yz_zy_i = gam2ZY(24 * j + 23);
                auto gam_sig_xx_zz_r = gam2ZZ(24 * j);
                auto gam_sig_xx_zz_i = gam2ZZ(24 * j + 1);
                auto gam_sig_yy_zz_r = gam2ZZ(24 * j + 2);
                auto gam_sig_yy_zz_i = gam2ZZ(24 * j + 3);
                auto gam_sig_zz_zz_r = gam2ZZ(24 * j + 4);
                auto gam_sig_zz_zz_i = gam2ZZ(24 * j + 5);
                auto gam_sig_xy_zz_r = gam2ZZ(24 * j + 6);
                auto gam_sig_xy_zz_i = gam2ZZ(24 * j + 7);
                auto gam_sig_xz_zz_r = gam2ZZ(24 * j + 8);
                auto gam_sig_xz_zz_i = gam2ZZ(24 * j + 9);
                auto gam_sig_yz_zz_r = gam2ZZ(24 * j + 10);
                auto gam_sig_yz_zz_i = gam2ZZ(24 * j + 11);
                auto gam_lamtau_xx_zz_r = gam2ZZ(24 * j + 12);
                auto gam_lamtau_xx_zz_i = gam2ZZ(24 * j + 13);
                auto gam_lamtau_yy_zz_r = gam2ZZ(24 * j + 14);
                auto gam_lamtau_yy_zz_i = gam2ZZ(24 * j + 15);
                auto gam_lamtau_zz_zz_r = gam2ZZ(24 * j + 16);
                auto gam_lamtau_zz_zz_i = gam2ZZ(24 * j + 17);
                auto gam_lamtau_xy_zz_r = gam2ZZ(24 * j + 18);
                auto gam_lamtau_xy_zz_i = gam2ZZ(24 * j + 19);
                auto gam_lamtau_xz_zz_r = gam2ZZ(24 * j + 20);
                auto gam_lamtau_xz_zz_i = gam2ZZ(24 * j + 21);
                auto gam_lamtau_yz_zz_r = gam2ZZ(24 * j + 22);
                auto gam_lamtau_yz_zz_i = gam2ZZ(24 * j + 23);


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

                auto rho_sig_xx_r = rw2DensityGrid.alphaDensity(24 * j  + 0 );
                auto rho_sig_xx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 0 );
                auto rho_sig_xx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 0 );
                auto rho_sig_xx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 0 );
                auto rho_sig_xx_i = rw2DensityGrid.alphaDensity(24 * j + 1 );
                auto rho_sig_xx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  1 );
                auto rho_sig_xx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  1 );
                auto rho_sig_xx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  1 );
                auto rho_sig_yy_r = rw2DensityGrid.alphaDensity(24 * j  + 2 );
                auto rho_sig_yy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 2 );
                auto rho_sig_yy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 2 );
                auto rho_sig_yy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 2 );
                auto rho_sig_yy_i = rw2DensityGrid.alphaDensity(24 * j + 3 );
                auto rho_sig_yy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  3 );
                auto rho_sig_yy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  3 );
                auto rho_sig_yy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  3 );
                auto rho_sig_zz_r = rw2DensityGrid.alphaDensity(24 * j  + 4 );
                auto rho_sig_zz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 4 );
                auto rho_sig_zz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 4 );
                auto rho_sig_zz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 4 );
                auto rho_sig_zz_i = rw2DensityGrid.alphaDensity(24 * j + 5 );
                auto rho_sig_zz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  5 );
                auto rho_sig_zz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  5 );
                auto rho_sig_zz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  5 );
                auto rho_sig_xy_r = rw2DensityGrid.alphaDensity(24 * j  + 6 );
                auto rho_sig_xy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 6 );
                auto rho_sig_xy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 6 );
                auto rho_sig_xy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 6 );
                auto rho_sig_xy_i = rw2DensityGrid.alphaDensity(24 * j + 7 );
                auto rho_sig_xy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  7 );
                auto rho_sig_xy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  7 );
                auto rho_sig_xy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  7 );
                auto rho_sig_xz_r = rw2DensityGrid.alphaDensity(24 * j  + 8 );
                auto rho_sig_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 8 );
                auto rho_sig_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 8 );
                auto rho_sig_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 8 );
                auto rho_sig_xz_i = rw2DensityGrid.alphaDensity(24 * j + 9 );
                auto rho_sig_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  9 );
                auto rho_sig_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  9 );
                auto rho_sig_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  9 );
                auto rho_sig_yz_r = rw2DensityGrid.alphaDensity(24 * j  + 10 );
                auto rho_sig_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 10 );
                auto rho_sig_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 10 );
                auto rho_sig_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 10 );
                auto rho_sig_yz_i = rw2DensityGrid.alphaDensity(24 * j + 11 );
                auto rho_sig_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  11 );
                auto rho_sig_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  11 );
                auto rho_sig_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  11 );
                auto rho_lamtau_xx_r =   rw2DensityGrid.alphaDensity(24 * j  + 12 );
                auto rho_lamtau_xx_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 12 );
                auto rho_lamtau_xx_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 12 );
                auto rho_lamtau_xx_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 12 );
                auto rho_lamtau_xx_i =   rw2DensityGrid.alphaDensity(24 * j + 13 );
                auto rho_lamtau_xx_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  13 );
                auto rho_lamtau_xx_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  13 );
                auto rho_lamtau_xx_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  13 );
                auto rho_lamtau_yy_r =   rw2DensityGrid.alphaDensity(24 * j  + 14 );
                auto rho_lamtau_yy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 14 );
                auto rho_lamtau_yy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 14 );
                auto rho_lamtau_yy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 14 );
                auto rho_lamtau_yy_i =   rw2DensityGrid.alphaDensity(24 * j + 15 );
                auto rho_lamtau_yy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  15 );
                auto rho_lamtau_yy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  15 );
                auto rho_lamtau_yy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  15 );
                auto rho_lamtau_zz_r =   rw2DensityGrid.alphaDensity(24 * j  + 16 );
                auto rho_lamtau_zz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 16 );
                auto rho_lamtau_zz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 16 );
                auto rho_lamtau_zz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 16 );
                auto rho_lamtau_zz_i =   rw2DensityGrid.alphaDensity(24 * j + 17 );
                auto rho_lamtau_zz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  17 );
                auto rho_lamtau_zz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  17 );
                auto rho_lamtau_zz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  17 );
                auto rho_lamtau_xy_r =   rw2DensityGrid.alphaDensity(24 * j  + 18 );
                auto rho_lamtau_xy_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 18 );
                auto rho_lamtau_xy_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 18 );
                auto rho_lamtau_xy_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 18 );
                auto rho_lamtau_xy_i =   rw2DensityGrid.alphaDensity(24 * j + 19 );
                auto rho_lamtau_xy_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  19 );
                auto rho_lamtau_xy_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  19 );
                auto rho_lamtau_xy_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  19 );
                auto rho_lamtau_xz_r =   rw2DensityGrid.alphaDensity(24 * j  + 20 );
                auto rho_lamtau_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 20 );
                auto rho_lamtau_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 20 );
                auto rho_lamtau_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 20 );
                auto rho_lamtau_xz_i =   rw2DensityGrid.alphaDensity(24 * j + 21 );
                auto rho_lamtau_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  21 );
                auto rho_lamtau_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  21 );
                auto rho_lamtau_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  21 );
                auto rho_lamtau_yz_r =   rw2DensityGrid.alphaDensity(24 * j  + 22 );
                auto rho_lamtau_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 22 );
                auto rho_lamtau_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 22 );
                auto rho_lamtau_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 22 );
                auto rho_lamtau_yz_i =   rw2DensityGrid.alphaDensity(24 * j + 23 );
                auto rho_lamtau_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j +  23 );
                auto rho_lamtau_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j +  23 );
                auto rho_lamtau_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j +  23 );

                for (int32_t i = 0; i < npoints; i++)

                {
                    double sig_term_r =  + 6.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);

                    double sig_term_i =  + 6.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    double lamtau_term_r =  + 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])
                                            + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])
                                            + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]);                  

                    double lamtau_term_i =  + 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])
                                            + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])
                                            + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]);                  

                    double sig_term_x_r =  + 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_y_r =  + 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_z_r =  + 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_x_i =  + 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_y_i =  + 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_z_i =  + 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_xx_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_xy_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_xz_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_yx_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_yy_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_yz_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_zx_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_zy_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_zz_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_xx_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_xy_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_xz_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_yx_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_yy_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_yz_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_zx_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_zy_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_zz_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double lamtau_term_x_r =  + 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_y_r =  + 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_z_r =  + 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_x_i =  + 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_y_i =  + 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_z_i =  + 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                           + 12.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                           + 12.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                           + 12.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    double lamtau_term_xx_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_xy_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_xz_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    double lamtau_term_yx_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_yy_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_yz_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    double lamtau_term_zx_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_zy_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_zz_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    double lamtau_term_xx_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_xy_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_xz_i =  + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    double lamtau_term_yx_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_yy_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_yz_i =  + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    double lamtau_term_zx_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                            + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                            + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                            + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    double lamtau_term_zy_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                            + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                            + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                            + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    double lamtau_term_zz_i =  + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                            + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                            + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                            + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_sig_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_r;                 

                    gam_sig_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_r;                 

                    gam_sig_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_r;                 

                    gam_sig_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);                  

                    gam_sig_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    gam_sig_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    gam_sig_xx_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_i;                 

                    gam_sig_yy_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_i;                 

                    gam_sig_zz_i[i] = 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_i;                 

                    gam_sig_xy_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_xz_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_sig_yz_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i]) 
                                       + 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])  + lamtau_term_r;                  

                    gam_lamtau_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])  + lamtau_term_r;                  

                    gam_lamtau_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i])  + lamtau_term_r;                  

                    gam_lamtau_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoCy_r[i],rhoCy_i[i]);
                    gam_lamtau_xx_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i]) 
                                       + 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCx_r[i],rhoCx_i[i])  + lamtau_term_i;                  

                    gam_lamtau_yy_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCy_r[i],rhoCy_i[i])  + lamtau_term_i;                  

                    gam_lamtau_zz_i[i] = 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCz_r[i],rhoCz_i[i])  + lamtau_term_i;                  

                    gam_lamtau_xy_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCy_r[i],rhoCy_i[i]) 
                                       + 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_xz_i[i] = 12.0 * prod2_i(rhoBx_r[i],rhoBx_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCx_r[i],rhoCx_i[i]);
                    gam_lamtau_yz_i[i] = 12.0 * prod2_i(rhoBy_r[i],rhoBy_i[i],rhoCz_r[i],rhoCz_i[i]) 
                                       + 12.0 * prod2_i(rhoBz_r[i],rhoBz_i[i],rhoCy_r[i],rhoCy_i[i]);
                    gam_sig_xx_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_x_r;                   

                    gam_sig_xx_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_y_r;                   

                    gam_sig_xx_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_z_r;                   

                    gam_sig_yy_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_x_r;                   

                    gam_sig_yy_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_y_r;                   

                    gam_sig_yy_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_z_r;                   

                    gam_sig_zz_x_r[i] = 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_x_r;                   

                    gam_sig_zz_y_r[i] = 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_y_r;                   

                    gam_sig_zz_z_r[i] = 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_z_r;                   

                    gam_sig_xy_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_yz_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_xx_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                       + 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_x_i;                  

                    gam_sig_xx_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                       + 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_y_i;                  

                    gam_sig_xx_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                       + 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_z_i;                  

                    gam_sig_yy_x_i[i] = 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                       + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_x_i;                  

                    gam_sig_yy_y_i[i] = 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                       + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_y_i;                  

                    gam_sig_yy_z_i[i] = 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                       + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_z_i;                  

                    gam_sig_zz_x_i[i] = 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                       + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_x_i;                  

                    gam_sig_zz_y_i[i] = 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                       + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_y_i;                  

                    gam_sig_zz_z_i[i] = 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                       + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_z_i;                  

                    gam_sig_xy_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_yz_x_i[i] = 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_y_i[i] = 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_z_i[i] = 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_xx_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_xx_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_xx_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_xx_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_xx_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_xx_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_xx_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_xx_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_xx_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_yy_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_yy_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_yy_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_yy_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_yy_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_yy_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_yy_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_yy_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_yy_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_zz_xx_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_zz_xy_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_zz_xz_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_zz_yx_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_zz_yy_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_zz_yz_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_zz_zx_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_zz_zy_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_zz_zz_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_xy_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_yz_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_xx_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_xx_i;                    

                    gam_sig_xx_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_xy_i;                    

                    gam_sig_xx_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_xz_i;                    

                    gam_sig_xx_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_yx_i;                    

                    gam_sig_xx_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_yy_i;                    

                    gam_sig_xx_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_yz_i;                    

                    gam_sig_xx_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_zx_i;                    

                    gam_sig_xx_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_zy_i;                    

                    gam_sig_xx_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_zz_i;                    

                    gam_sig_yy_xx_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_xx_i;                    

                    gam_sig_yy_xy_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_xy_i;                    

                    gam_sig_yy_xz_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_xz_i;                    

                    gam_sig_yy_yx_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_yx_i;                    

                    gam_sig_yy_yy_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_yy_i;                    

                    gam_sig_yy_yz_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_yz_i;                    

                    gam_sig_yy_zx_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_zx_i;                    

                    gam_sig_yy_zy_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_zy_i;                    

                    gam_sig_yy_zz_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_zz_i;                    

                    gam_sig_zz_xx_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_xx_i;                    

                    gam_sig_zz_xy_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_xy_i;                    

                    gam_sig_zz_xz_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_xz_i;                    

                    gam_sig_zz_yx_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_yx_i;                    

                    gam_sig_zz_yy_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_yy_i;                    

                    gam_sig_zz_yz_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_yz_i;                    

                    gam_sig_zz_zx_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_zx_i;                    

                    gam_sig_zz_zy_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_zy_i;                    

                    gam_sig_zz_zz_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_zz_i;                    

                    gam_sig_xy_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_yz_xx_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_xy_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_xz_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_yx_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_yy_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_yz_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_zx_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_zy_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_zz_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xx_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_x_r;                 

                    gam_lamtau_xx_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_y_r;                 

                    gam_lamtau_xx_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_z_r;                 

                    gam_lamtau_yy_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_x_r;                 

                    gam_lamtau_yy_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_y_r;                 

                    gam_lamtau_yy_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_z_r;                 

                    gam_lamtau_zz_x_r[i] = 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_x_r;                 

                    gam_lamtau_zz_y_r[i] = 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_y_r;                 

                    gam_lamtau_zz_z_r[i] = 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_z_r;                 

                    gam_lamtau_xy_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xy_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xy_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xz_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xz_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xz_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xx_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_x_i;                 

                    gam_lamtau_xx_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_y_i;                 

                    gam_lamtau_xx_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],rhoBx_r[i],rhoBx_i[i]) + lamtau_term_z_i;                 

                    gam_lamtau_yy_x_i[i] = 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_x_i;                 

                    gam_lamtau_yy_y_i[i] = 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_y_i;                 

                    gam_lamtau_yy_z_i[i] = 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],rhoBy_r[i],rhoBy_i[i]) + lamtau_term_z_i;                 

                    gam_lamtau_zz_x_i[i] = 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_x_i;                 

                    gam_lamtau_zz_y_i[i] = 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_y_i;                 

                    gam_lamtau_zz_z_i[i] = 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                         + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],rhoBz_r[i],rhoBz_i[i]) + lamtau_term_z_i;                 

                    gam_lamtau_xy_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xy_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xy_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_lamtau_xz_x_i[i] = 12.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xz_y_i[i] = 12.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xz_z_i[i] = 12.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoCx_r[i],rhoCx_i[i])
                                         + 12.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_x_i[i] = 12.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_y_i[i] = 12.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_yz_z_i[i] = 12.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],rhoCz_r[i],rhoCz_i[i])
                                         + 12.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 12.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],rhoCy_r[i],rhoCy_i[i])
                                         + 12.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],rhoBz_r[i],rhoBz_i[i]);
                    gam_lamtau_xx_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_xx_r;                  

                    gam_lamtau_xx_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_xy_r;                  

                    gam_lamtau_xx_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_xz_r;                  

                    gam_lamtau_xx_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_yx_r;                  

                    gam_lamtau_xx_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_yy_r;                  

                    gam_lamtau_xx_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_yz_r;                  

                    gam_lamtau_xx_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_zx_r;                  

                    gam_lamtau_xx_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_zy_r;                  

                    gam_lamtau_xx_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_zz_r;                  

                    gam_lamtau_yy_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_xx_r;                  

                    gam_lamtau_yy_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_xy_r;                  

                    gam_lamtau_yy_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_xz_r;                  

                    gam_lamtau_yy_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_yx_r;                  

                    gam_lamtau_yy_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_yy_r;                  

                    gam_lamtau_yy_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_yz_r;                  

                    gam_lamtau_yy_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_zx_r;                  

                    gam_lamtau_yy_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_zy_r;                  

                    gam_lamtau_yy_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_zz_r;                  

                    gam_lamtau_zz_xx_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_xx_r;                  

                    gam_lamtau_zz_xy_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_xy_r;                  

                    gam_lamtau_zz_xz_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_xz_r;                  

                    gam_lamtau_zz_yx_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_yx_r;                  

                    gam_lamtau_zz_yy_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_yy_r;                  

                    gam_lamtau_zz_yz_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_yz_r;                  

                    gam_lamtau_zz_zx_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_zx_r;                  

                    gam_lamtau_zz_zy_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_zy_r;                  

                    gam_lamtau_zz_zz_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_zz_r;                  

                    gam_lamtau_xy_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xy_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xy_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xz_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_x_r[i],gradCx_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_xz_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_y_r[i],gradCx_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_xz_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_r(gradCx_z_r[i],gradCx_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_x_r[i],gradCz_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_x_r[i],gradCy_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_y_r[i],gradCz_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_y_r[i],gradCy_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_r(gradCz_z_r[i],gradCz_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_r(gradCy_z_r[i],gradCy_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_xx_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_xx_i;                  

                    gam_lamtau_xx_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_xy_i;                  

                    gam_lamtau_xx_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_xz_i;                  

                    gam_lamtau_xx_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_yx_i;                  

                    gam_lamtau_xx_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_yy_i;                  

                    gam_lamtau_xx_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_yz_i;                  

                    gam_lamtau_xx_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_x_r[i],gradBx_x_i[i]) + lamtau_term_zx_i;                  

                    gam_lamtau_xx_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_y_r[i],gradBx_y_i[i]) + lamtau_term_zy_i;                  

                    gam_lamtau_xx_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBx_z_r[i],gradBx_z_i[i]) + lamtau_term_zz_i;                  

                    gam_lamtau_yy_xx_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_xx_i;                  

                    gam_lamtau_yy_xy_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_xy_i;                  

                    gam_lamtau_yy_xz_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_xz_i;                  

                    gam_lamtau_yy_yx_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_yx_i;                  

                    gam_lamtau_yy_yy_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_yy_i;                  

                    gam_lamtau_yy_yz_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_yz_i;                  

                    gam_lamtau_yy_zx_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_x_r[i],gradBy_x_i[i]) + lamtau_term_zx_i;                  

                    gam_lamtau_yy_zy_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_y_r[i],gradBy_y_i[i]) + lamtau_term_zy_i;                  

                    gam_lamtau_yy_zz_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBy_z_r[i],gradBy_z_i[i]) + lamtau_term_zz_i;                  

                    gam_lamtau_zz_xx_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_xx_i;                  

                    gam_lamtau_zz_xy_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_xy_i;                  

                    gam_lamtau_zz_xz_i[i] = 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_xz_i;                  

                    gam_lamtau_zz_yx_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_yx_i;                  

                    gam_lamtau_zz_yy_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_yy_i;                  

                    gam_lamtau_zz_yz_i[i] = 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_yz_i;                  

                    gam_lamtau_zz_zx_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]) + lamtau_term_zx_i;                  

                    gam_lamtau_zz_zy_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]) + lamtau_term_zy_i;                  

                    gam_lamtau_zz_zz_i[i] = 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]) + lamtau_term_zz_i;                  

                    gam_lamtau_xy_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xy_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xy_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_lamtau_xy_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_lamtau_xy_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_lamtau_xz_xx_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_xy_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_xz_i[i] = 6.0 * prod2_i(gradBx_x_r[i],gradBx_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_x_r[i],gradCx_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_xz_yx_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_yy_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_yz_i[i] = 6.0 * prod2_i(gradBx_y_r[i],gradBx_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_y_r[i],gradCx_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_xz_zx_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_xz_zy_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_xz_zz_i[i] = 6.0 * prod2_i(gradBx_z_r[i],gradBx_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                          + 6.0 * prod2_i(gradCx_z_r[i],gradCx_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_xx_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_xy_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_xz_i[i] = 6.0 * prod2_i(gradBy_x_r[i],gradBy_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_x_r[i],gradCz_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBz_x_r[i],gradBz_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_x_r[i],gradCy_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_yx_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_yy_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_yz_i[i] = 6.0 * prod2_i(gradBy_y_r[i],gradBy_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_y_r[i],gradCz_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBz_y_r[i],gradBz_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_y_r[i],gradCy_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                    gam_lamtau_yz_zx_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);
                    gam_lamtau_yz_zy_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);
                    gam_lamtau_yz_zz_i[i] = 6.0 * prod2_i(gradBy_z_r[i],gradBy_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                          + 6.0 * prod2_i(gradCz_z_r[i],gradCz_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                          + 6.0 * prod2_i(gradBz_z_r[i],gradBz_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                          + 6.0 * prod2_i(gradCy_z_r[i],gradCy_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_r[i] = prod2_r(rho_sig_xx_r[i],rho_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamy_r[i] = prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_yy_r[i],rho_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamz_r[i] = prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_zz_r[i],rho_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamx_i[i] = prod2_i(rho_sig_xx_r[i],rho_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamy_i[i] = prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_yy_r[i],rho_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamz_i[i] = prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_zz_r[i],rho_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                   

                    gamx_X_r[i] =  prod2_r(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_Y_r[i] =  prod2_r(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_Z_r[i] =  prod2_r(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_X_r[i] =  prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_Y_r[i] =  prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_Z_r[i] =  prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_X_r[i] =  prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_Y_r[i] =  prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_Z_r[i] =  prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_r(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_r(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_r(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_r(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_r(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_r(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_X_i[i] =  prod2_i(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_Y_i[i] =  prod2_i(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_Z_i[i] =  prod2_i(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xx_r[i],rho_sig_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xx_r[i],rho_lamtau_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_X_i[i] =  prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_Y_i[i] =  prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_Z_i[i] =  prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xy_r[i],rho_sig_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yy_r[i],rho_sig_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xy_r[i],rho_lamtau_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yy_r[i],rho_lamtau_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_X_i[i] =  prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_Y_i[i] =  prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_Z_i[i] =  prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],rhoCx_r[i],rhoCx_i[i])
                              + prod2_i(rho_sig_xz_r[i],rho_sig_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],rhoCy_r[i],rhoCy_i[i])
                              + prod2_i(rho_sig_yz_r[i],rho_sig_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],rhoCz_r[i],rhoCz_i[i])
                              + prod2_i(rho_sig_zz_r[i],rho_sig_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],rhoBx_r[i],rhoBx_i[i])
                              + prod2_i(rho_lamtau_xz_r[i],rho_lamtau_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],rhoBy_r[i],rhoBy_i[i])
                              + prod2_i(rho_lamtau_yz_r[i],rho_lamtau_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                              + prod2_i(rho_lamtau_zz_r[i],rho_lamtau_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_XX_r[i] = prod2_r(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_XY_r[i] = prod2_r(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_XZ_r[i] = prod2_r(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_YX_r[i] = prod2_r(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_YY_r[i] = prod2_r(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_YZ_r[i] = prod2_r(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_ZX_r[i] = prod2_r(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_ZY_r[i] = prod2_r(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_ZZ_r[i] = prod2_r(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_XX_r[i] = prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_XY_r[i] = prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_XZ_r[i] = prod2_r(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_YX_r[i] = prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_YY_r[i] = prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_YZ_r[i] = prod2_r(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_ZX_r[i] = prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_ZY_r[i] = prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_ZZ_r[i] = prod2_r(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_XX_r[i] = prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_XY_r[i] = prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_XZ_r[i] = prod2_r(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_YX_r[i] = prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_YY_r[i] = prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_YZ_r[i] = prod2_r(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_ZX_r[i] = prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_r(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_r(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_ZY_r[i] = prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_r(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_r(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_ZZ_r[i] = prod2_r(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_r(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_r(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_r(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_r(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_r(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_XX_i[i] = prod2_i(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_XY_i[i] = prod2_i(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_XZ_i[i] = prod2_i(rho_sig_xx_x_r[i],rho_sig_xx_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xx_x_r[i],rho_lamtau_xx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_YX_i[i] = prod2_i(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_YY_i[i] = prod2_i(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_YZ_i[i] = prod2_i(rho_sig_xx_y_r[i],rho_sig_xx_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xx_y_r[i],rho_lamtau_xx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamx_ZX_i[i] = prod2_i(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamx_ZY_i[i] = prod2_i(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamx_ZZ_i[i] = prod2_i(rho_sig_xx_z_r[i],rho_sig_xx_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xx_z_r[i],rho_lamtau_xx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_XX_i[i] = prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_XY_i[i] = prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_XZ_i[i] = prod2_i(rho_sig_xy_x_r[i],rho_sig_xy_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yy_x_r[i],rho_sig_yy_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xy_x_r[i],rho_lamtau_xy_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yy_x_r[i],rho_lamtau_yy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_YX_i[i] = prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_YY_i[i] = prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_YZ_i[i] = prod2_i(rho_sig_xy_y_r[i],rho_sig_xy_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yy_y_r[i],rho_sig_yy_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xy_y_r[i],rho_lamtau_xy_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yy_y_r[i],rho_lamtau_yy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamy_ZX_i[i] = prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamy_ZY_i[i] = prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamy_ZZ_i[i] = prod2_i(rho_sig_xy_z_r[i],rho_sig_xy_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yy_z_r[i],rho_sig_yy_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xy_z_r[i],rho_lamtau_xy_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yy_z_r[i],rho_lamtau_yy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_XX_i[i] = prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_XY_i[i] = prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_XZ_i[i] = prod2_i(rho_sig_xz_x_r[i],rho_sig_xz_x_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yz_x_r[i],rho_sig_yz_x_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_zz_x_r[i],rho_sig_zz_x_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xz_x_r[i],rho_lamtau_xz_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yz_x_r[i],rho_lamtau_yz_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_zz_x_r[i],rho_lamtau_zz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_YX_i[i] = prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_YY_i[i] = prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_YZ_i[i] = prod2_i(rho_sig_xz_y_r[i],rho_sig_xz_y_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yz_y_r[i],rho_sig_yz_y_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_zz_y_r[i],rho_sig_zz_y_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xz_y_r[i],rho_lamtau_xz_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yz_y_r[i],rho_lamtau_yz_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_zz_y_r[i],rho_lamtau_zz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    gamz_ZX_i[i] = prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_x_r[i],gradCx_x_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_x_r[i],gradCy_x_i[i])
                              + prod2_i(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_x_r[i],gradCz_x_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                              + prod2_i(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    gamz_ZY_i[i] = prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_y_r[i],gradCx_y_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_y_r[i],gradCy_y_i[i])
                              + prod2_i(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_y_r[i],gradCz_y_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                              + prod2_i(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    gamz_ZZ_i[i] = prod2_i(rho_sig_xz_z_r[i],rho_sig_xz_z_i[i],gradCx_z_r[i],gradCx_z_i[i])
                              + prod2_i(rho_sig_yz_z_r[i],rho_sig_yz_z_i[i],gradCy_z_r[i],gradCy_z_i[i])
                              + prod2_i(rho_sig_zz_z_r[i],rho_sig_zz_z_i[i],gradCz_z_r[i],gradCz_z_i[i])
                              + prod2_i(rho_lamtau_xz_z_r[i],rho_lamtau_xz_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                              + prod2_i(rho_lamtau_yz_z_r[i],rho_lamtau_yz_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                              + prod2_i(rho_lamtau_zz_z_r[i],rho_lamtau_zz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_r[i],gam_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                 


                    piy_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yy_r[i],gam_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                 


                    piz_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_r(gam_sig_zz_r[i],gam_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_r(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                 


                    pix_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_r[i],gam_sig_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    piy_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yy_r[i],gam_sig_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    piz_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                              + 1.0/3.0 * prod2_i(gam_sig_zz_r[i],gam_sig_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                              + 1.0/3.0 * prod2_i(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    pix_X_r[i] =  prod2_r(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_r(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    pix_Y_r[i] =  prod2_r(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_r(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    pix_Z_r[i] =  prod2_r(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_r(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    piy_X_r[i] =  prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_r(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_r(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    piy_Y_r[i] =  prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_r(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_r(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    piy_Z_r[i] =  prod2_r(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_r(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_r(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_r(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    piz_X_r[i] =  prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_r(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_r(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    piz_Y_r[i] =  prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_r(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_r(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    piz_Z_r[i] =  prod2_r(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_r(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_r(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_r(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_r(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_r(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    pix_X_i[i] =  prod2_i(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_i(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    pix_Y_i[i] =  prod2_i(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_i(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    pix_Z_i[i] =  prod2_i(gam_sig_xx_r[i],gam_sig_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_i(gam_lamtau_xx_r[i],gam_lamtau_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    piy_X_i[i] =  prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_i(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_i(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    piy_Y_i[i] =  prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_i(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_i(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    piy_Z_i[i] =  prod2_i(gam_sig_xy_r[i],gam_sig_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_i(gam_lamtau_xy_r[i],gam_lamtau_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_i(gam_sig_yy_r[i],gam_sig_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_i(gam_lamtau_yy_r[i],gam_lamtau_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    piz_X_i[i] =  prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                + prod2_i(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                + prod2_i(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    piz_Y_i[i] =  prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                + prod2_i(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                + prod2_i(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    piz_Z_i[i] =  prod2_i(gam_sig_xz_r[i],gam_sig_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                + prod2_i(gam_lamtau_xz_r[i],gam_lamtau_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                + prod2_i(gam_sig_yz_r[i],gam_sig_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                + prod2_i(gam_lamtau_yz_r[i],gam_lamtau_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                + prod2_i(gam_sig_zz_r[i],gam_sig_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                + prod2_i(gam_lamtau_zz_r[i],gam_lamtau_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    pix_XX_r[i] =  prod2_r(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XY_r[i] =  prod2_r(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XZ_r[i] =  prod2_r(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YX_r[i] =  prod2_r(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YY_r[i] =  prod2_r(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YZ_r[i] =  prod2_r(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZX_r[i] =  prod2_r(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZY_r[i] =  prod2_r(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZZ_r[i] =  prod2_r(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XX_r[i] =  prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XY_r[i] =  prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XZ_r[i] =  prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YX_r[i] =  prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YY_r[i] =  prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YZ_r[i] =  prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZX_r[i] =  prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZY_r[i] =  prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZZ_r[i] =  prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XX_r[i] =  prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XY_r[i] =  prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XZ_r[i] =  prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YX_r[i] =  prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YY_r[i] =  prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YZ_r[i] =  prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZX_r[i] =  prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZY_r[i] =  prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZZ_r[i] =  prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_r(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_r(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XX_i[i] =  prod2_i(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XY_i[i] =  prod2_i(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XZ_i[i] =  prod2_i(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YX_i[i] =  prod2_i(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YY_i[i] =  prod2_i(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_YZ_i[i] =  prod2_i(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZX_i[i] =  prod2_i(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZY_i[i] =  prod2_i(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_ZZ_i[i] =  prod2_i(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XX_i[i] =  prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XY_i[i] =  prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_XZ_i[i] =  prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YX_i[i] =  prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YY_i[i] =  prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_YZ_i[i] =  prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZX_i[i] =  prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZY_i[i] =  prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piy_ZZ_i[i] =  prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XX_i[i] =  prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XY_i[i] =  prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_XZ_i[i] =  prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YX_i[i] =  prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YY_i[i] =  prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_YZ_i[i] =  prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZX_i[i] =  prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZY_i[i] =  prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    piz_ZZ_i[i] =  prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],rhoCx_r[i],rhoCx_i[i])
                                 + prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],rhoBx_r[i],rhoBx_i[i])
                                 + prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],rhoCy_r[i],rhoCy_i[i])
                                 + prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],rhoBy_r[i],rhoBy_i[i])
                                 + prod2_i(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],rhoCz_r[i],rhoCz_i[i])
                                 + prod2_i(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    pix_XXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_XYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_XZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZXX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZXY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZXZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZYX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZYY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZYZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZZX_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZZY_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZZZ_r[i] =  1.0/3.0 * prod2_r(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_r(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_XXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xx_r[i],gam_sig_xx_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xx_r[i],gam_lamtau_xx_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_XYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xy_r[i],gam_sig_xx_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xy_r[i],gam_lamtau_xx_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_XZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_XZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_XZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_xz_r[i],gam_sig_xx_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_xz_r[i],gam_lamtau_xx_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yx_r[i],gam_sig_xx_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yx_r[i],gam_lamtau_xx_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yy_r[i],gam_sig_xx_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yy_r[i],gam_lamtau_xx_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_YZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_YZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_YZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_yz_r[i],gam_sig_xx_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_yz_r[i],gam_lamtau_xx_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zx_r[i],gam_sig_xx_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zx_r[i],gam_lamtau_xx_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zy_r[i],gam_sig_xx_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zy_r[i],gam_lamtau_xx_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    pix_ZZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    pix_ZZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    pix_ZZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xx_zz_r[i],gam_sig_xx_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xx_zz_r[i],gam_lamtau_xx_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xx_r[i],gam_sig_xy_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xx_r[i],gam_lamtau_xy_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xx_r[i],gam_sig_yy_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xx_r[i],gam_lamtau_yy_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xy_r[i],gam_sig_xy_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xy_r[i],gam_lamtau_xy_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xy_r[i],gam_sig_yy_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xy_r[i],gam_lamtau_yy_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_XZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_XZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_XZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_xz_r[i],gam_sig_xy_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_xz_r[i],gam_lamtau_xy_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_xz_r[i],gam_sig_yy_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_xz_r[i],gam_lamtau_yy_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yx_r[i],gam_sig_xy_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yx_r[i],gam_lamtau_xy_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yx_r[i],gam_sig_yy_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yx_r[i],gam_lamtau_yy_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yy_r[i],gam_sig_xy_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yy_r[i],gam_lamtau_xy_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yy_r[i],gam_sig_yy_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yy_r[i],gam_lamtau_yy_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_YZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_YZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_YZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_yz_r[i],gam_sig_xy_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_yz_r[i],gam_lamtau_xy_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_yz_r[i],gam_sig_yy_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_yz_r[i],gam_lamtau_yy_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zx_r[i],gam_sig_xy_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zx_r[i],gam_lamtau_xy_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zx_r[i],gam_sig_yy_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zx_r[i],gam_lamtau_yy_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zy_r[i],gam_sig_xy_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zy_r[i],gam_lamtau_xy_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zy_r[i],gam_sig_yy_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zy_r[i],gam_lamtau_yy_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piy_ZZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piy_ZZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piy_ZZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xy_zz_r[i],gam_sig_xy_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xy_zz_r[i],gam_lamtau_xy_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yy_zz_r[i],gam_sig_yy_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yy_zz_r[i],gam_lamtau_yy_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xx_r[i],gam_sig_xz_xx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xx_r[i],gam_lamtau_xz_xx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xx_r[i],gam_sig_yz_xx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xx_r[i],gam_lamtau_yz_xx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xx_r[i],gam_sig_zz_xx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xx_r[i],gam_lamtau_zz_xx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xy_r[i],gam_sig_xz_xy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xy_r[i],gam_lamtau_xz_xy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xy_r[i],gam_sig_yz_xy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xy_r[i],gam_lamtau_yz_xy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xy_r[i],gam_sig_zz_xy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xy_r[i],gam_lamtau_zz_xy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_XZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_XZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_XZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_xz_r[i],gam_sig_xz_xz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_xz_r[i],gam_lamtau_xz_xz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_xz_r[i],gam_sig_yz_xz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_xz_r[i],gam_lamtau_yz_xz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_xz_r[i],gam_sig_zz_xz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_xz_r[i],gam_lamtau_zz_xz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yx_r[i],gam_sig_xz_yx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yx_r[i],gam_lamtau_xz_yx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yx_r[i],gam_sig_yz_yx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yx_r[i],gam_lamtau_yz_yx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yx_r[i],gam_sig_zz_yx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yx_r[i],gam_lamtau_zz_yx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yy_r[i],gam_sig_xz_yy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yy_r[i],gam_lamtau_xz_yy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yy_r[i],gam_sig_yz_yy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yy_r[i],gam_lamtau_yz_yy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yy_r[i],gam_sig_zz_yy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yy_r[i],gam_lamtau_zz_yy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_YZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_YZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_YZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_yz_r[i],gam_sig_xz_yz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_yz_r[i],gam_lamtau_xz_yz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_yz_r[i],gam_sig_yz_yz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_yz_r[i],gam_lamtau_yz_yz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_yz_r[i],gam_sig_zz_yz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_yz_r[i],gam_lamtau_zz_yz_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZXX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZXY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZXZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zx_r[i],gam_sig_xz_zx_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zx_r[i],gam_lamtau_xz_zx_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zx_r[i],gam_sig_yz_zx_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zx_r[i],gam_lamtau_yz_zx_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zx_r[i],gam_sig_zz_zx_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zx_r[i],gam_lamtau_zz_zx_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZYX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZYY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZYZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zy_r[i],gam_sig_xz_zy_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zy_r[i],gam_lamtau_xz_zy_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zy_r[i],gam_sig_yz_zy_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zy_r[i],gam_lamtau_yz_zy_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zy_r[i],gam_sig_zz_zy_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zy_r[i],gam_lamtau_zz_zy_i[i],gradBz_z_r[i],gradBz_z_i[i]);                 

                    piz_ZZX_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_x_r[i],gradCx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_x_r[i],gradCy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_x_r[i],gradCz_x_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_x_r[i],gradBz_x_i[i]);                 

                    piz_ZZY_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_y_r[i],gradCx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_y_r[i],gradCy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_y_r[i],gradCz_y_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_y_r[i],gradBz_y_i[i]);                 

                    piz_ZZZ_i[i] =  1.0/3.0 * prod2_i(gam_sig_xz_zz_r[i],gam_sig_xz_zz_i[i],gradCx_z_r[i],gradCx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_xz_zz_r[i],gam_lamtau_xz_zz_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_yz_zz_r[i],gam_sig_yz_zz_i[i],gradCy_z_r[i],gradCy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_yz_zz_r[i],gam_lamtau_yz_zz_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_sig_zz_zz_r[i],gam_sig_zz_zz_i[i],gradCz_z_r[i],gradCz_z_i[i])
                                  + 1.0/3.0 * prod2_i(gam_lamtau_zz_zz_r[i],gam_lamtau_zz_zz_i[i],gradBz_z_r[i],gradBz_z_i[i]);
                }

            }
        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 8; j++)
            {
                
                auto gam_bc_r = gam2(6 * j);
                auto gam_bc_i = gam2(6 * j + 1);
                auto gam_bd_r = gam2(6 * j + 2);
                auto gam_bd_i = gam2(6 * j + 3);
                auto gam_cd_r = gam2(6 * j + 4);
                auto gam_cd_i = gam2(6 * j + 5);

                auto gam_bc_x_r = gam2X(6 * j);
                auto gam_bc_x_i = gam2X(6 * j + 1);
                auto gam_bd_x_r = gam2X(6 * j + 2);
                auto gam_bd_x_i = gam2X(6 * j + 3);
                auto gam_cd_x_r = gam2X(6 * j + 4);
                auto gam_cd_x_i = gam2X(6 * j + 5);
                auto gam_bc_y_r = gam2Y(6 * j);
                auto gam_bc_y_i = gam2Y(6 * j + 1);
                auto gam_bd_y_r = gam2Y(6 * j + 2);
                auto gam_bd_y_i = gam2Y(6 * j + 3);
                auto gam_cd_y_r = gam2Y(6 * j + 4);
                auto gam_cd_y_i = gam2Y(6 * j + 5);
                auto gam_bc_z_r = gam2Z(6 * j);
                auto gam_bc_z_i = gam2Z(6 * j + 1);
                auto gam_bd_z_r = gam2Z(6 * j + 2);
                auto gam_bd_z_i = gam2Z(6 * j + 3);
                auto gam_cd_z_r = gam2Z(6 * j + 4);
                auto gam_cd_z_i = gam2Z(6 * j + 5);

                auto gam_bc_xx_r = gam2XX(6 * j);
                auto gam_bc_xx_i = gam2XX(6 * j + 1);
                auto gam_bd_xx_r = gam2XX(6 * j + 2);
                auto gam_bd_xx_i = gam2XX(6 * j + 3);
                auto gam_cd_xx_r = gam2XX(6 * j + 4);
                auto gam_cd_xx_i = gam2XX(6 * j + 5);
                auto gam_bc_xy_r = gam2XY(6 * j);
                auto gam_bc_xy_i = gam2XY(6 * j + 1);
                auto gam_bd_xy_r = gam2XY(6 * j + 2);
                auto gam_bd_xy_i = gam2XY(6 * j + 3);
                auto gam_cd_xy_r = gam2XY(6 * j + 4);
                auto gam_cd_xy_i = gam2XY(6 * j + 5);
                auto gam_bc_xz_r = gam2XZ(6 * j);
                auto gam_bc_xz_i = gam2XZ(6 * j + 1);
                auto gam_bd_xz_r = gam2XZ(6 * j + 2);
                auto gam_bd_xz_i = gam2XZ(6 * j + 3);
                auto gam_cd_xz_r = gam2XZ(6 * j + 4);
                auto gam_cd_xz_i = gam2XZ(6 * j + 5);
                auto gam_bc_yx_r = gam2YX(6 * j);
                auto gam_bc_yx_i = gam2YX(6 * j + 1);
                auto gam_bd_yx_r = gam2YX(6 * j + 2);
                auto gam_bd_yx_i = gam2YX(6 * j + 3);
                auto gam_cd_yx_r = gam2YX(6 * j + 4);
                auto gam_cd_yx_i = gam2YX(6 * j + 5);
                auto gam_bc_yy_r = gam2YY(6 * j);
                auto gam_bc_yy_i = gam2YY(6 * j + 1);
                auto gam_bd_yy_r = gam2YY(6 * j + 2);
                auto gam_bd_yy_i = gam2YY(6 * j + 3);
                auto gam_cd_yy_r = gam2YY(6 * j + 4);
                auto gam_cd_yy_i = gam2YY(6 * j + 5);
                auto gam_bc_yz_r = gam2YZ(6 * j);
                auto gam_bc_yz_i = gam2YZ(6 * j + 1);
                auto gam_bd_yz_r = gam2YZ(6 * j + 2);
                auto gam_bd_yz_i = gam2YZ(6 * j + 3);
                auto gam_cd_yz_r = gam2YZ(6 * j + 4);
                auto gam_cd_yz_i = gam2YZ(6 * j + 5);
                auto gam_bc_zx_r = gam2ZX(6 * j);
                auto gam_bc_zx_i = gam2ZX(6 * j + 1);
                auto gam_bd_zx_r = gam2ZX(6 * j + 2);
                auto gam_bd_zx_i = gam2ZX(6 * j + 3);
                auto gam_cd_zx_r = gam2ZX(6 * j + 4);
                auto gam_cd_zx_i = gam2ZX(6 * j + 5);
                auto gam_bc_zy_r = gam2ZY(6 * j);
                auto gam_bc_zy_i = gam2ZY(6 * j + 1);
                auto gam_bd_zy_r = gam2ZY(6 * j + 2);
                auto gam_bd_zy_i = gam2ZY(6 * j + 3);
                auto gam_cd_zy_r = gam2ZY(6 * j + 4);
                auto gam_cd_zy_i = gam2ZY(6 * j + 5);
                auto gam_bc_zz_r = gam2ZZ(6 * j);
                auto gam_bc_zz_i = gam2ZZ(6 * j + 1);
                auto gam_bd_zz_r = gam2ZZ(6 * j + 2);
                auto gam_bd_zz_i = gam2ZZ(6 * j + 3);
                auto gam_cd_zz_r = gam2ZZ(6 * j + 4);
                auto gam_cd_zz_i = gam2ZZ(6 * j + 5);

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


                for (int32_t i = 0; i < npoints; i++)

                {
                    gam_bc_r[i] += prod2_r(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                                + prod2_r(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_r[i] += prod2_r(rhoB_r[i],rhoB_i[i],rhoD_r[i],rhoD_i[i])
                                + prod2_r(rhoD_r[i],rhoD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_r[i] += prod2_r(rhoC_r[i],rhoC_i[i],rhoD_r[i],rhoD_i[i])
                                + prod2_r(rhoD_r[i],rhoD_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_i[i] += prod2_i(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                                   + prod2_i(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_i[i] += prod2_i(rhoB_r[i],rhoB_i[i],rhoD_r[i],rhoD_i[i])
                                   + prod2_i(rhoD_r[i],rhoD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_i[i] += prod2_i(rhoC_r[i],rhoC_i[i],rhoD_r[i],rhoD_i[i])
                                   + prod2_i(rhoD_r[i],rhoD_i[i],rhoC_r[i],rhoC_i[i]);


                    gam_bc_x_r[i] += 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_x_r[i] += 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_x_r[i] += 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_y_r[i] += 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_y_r[i] += 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_y_r[i] += 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_z_r[i] += 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_z_r[i] += 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_z_r[i] += 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_x_i[i] += 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_x_i[i] += 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_x_i[i] += 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_y_i[i] += 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_y_i[i] += 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_y_i[i] += 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],rhoC_r[i],rhoC_i[i]);

                    gam_bc_z_i[i] += 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                   + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_bd_z_i[i] += 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_cd_z_i[i] += 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoD_r[i],rhoD_i[i])
                                   + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],rhoC_r[i],rhoC_i[i]);


                    gam_bc_xx_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_xx_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_xx_r[i] += prod2_r(gradC_x_r[i],gradC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_xy_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_xy_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_xy_r[i] += prod2_r(gradC_x_r[i],gradC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_xz_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_xz_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_xz_r[i] += prod2_r(gradC_x_r[i],gradC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_x_r[i],gradD_x_i[i],gradC_z_r[i],gradC_z_i[i]);

                    gam_bc_yx_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_yx_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_yx_r[i] += prod2_r(gradC_y_r[i],gradC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_yy_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_yy_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_yy_r[i] += prod2_r(gradC_y_r[i],gradC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_yz_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_yz_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_yz_r[i] += prod2_r(gradC_y_r[i],gradC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_y_r[i],gradD_y_i[i],gradC_z_r[i],gradC_z_i[i]);

                    gam_bc_zx_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_zx_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_zx_r[i] += prod2_r(gradC_z_r[i],gradC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_zy_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_zy_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_zy_r[i] += prod2_r(gradC_z_r[i],gradC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_zz_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_zz_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_zz_r[i] += prod2_r(gradC_z_r[i],gradC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_r(gradD_z_r[i],gradD_z_i[i],gradC_z_r[i],gradC_z_i[i]);

                    gam_bc_xx_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_xx_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_xx_i[i] += prod2_i(gradC_x_r[i],gradC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_xy_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_xy_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_xy_i[i] += prod2_i(gradC_x_r[i],gradC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_xz_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_xz_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_xz_i[i] += prod2_i(gradC_x_r[i],gradC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_x_r[i],gradD_x_i[i],gradC_z_r[i],gradC_z_i[i]);

                    gam_bc_yx_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_yx_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_yx_i[i] += prod2_i(gradC_y_r[i],gradC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_yy_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_yy_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_yy_i[i] += prod2_i(gradC_y_r[i],gradC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_yz_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_yz_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_yz_i[i] += prod2_i(gradC_y_r[i],gradC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_y_r[i],gradD_y_i[i],gradC_z_r[i],gradC_z_i[i]);

                    gam_bc_zx_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                   +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_bd_zx_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_cd_zx_i[i] += prod2_i(gradC_z_r[i],gradC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradC_x_r[i],gradC_x_i[i]);

                    gam_bc_zy_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                   +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_bd_zy_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_cd_zy_i[i] += prod2_i(gradC_z_r[i],gradC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradC_y_r[i],gradC_y_i[i]);

                    gam_bc_zz_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                   +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_bd_zz_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_cd_zz_i[i] += prod2_i(gradC_z_r[i],gradC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                   +  prod2_i(gradD_z_r[i],gradD_z_i[i],gradC_z_r[i],gradC_z_i[i]);


                    pi_r[i] = prod2_r(rhoB_r[i],rhoB_i[i],gam_cd_r[i],gam_cd_i[i])
                               + prod2_r(rhoC_r[i],rhoC_i[i],gam_bd_r[i],gam_bd_i[i])
                               + prod2_r(rhoD_r[i],rhoD_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_i[i] = prod2_i(rhoB_r[i],rhoB_i[i],gam_cd_r[i],gam_cd_i[i])
                               + prod2_i(rhoC_r[i],rhoC_i[i],gam_bd_r[i],gam_bd_i[i])
                               + prod2_i(rhoD_r[i],rhoD_i[i],gam_bc_r[i],gam_bc_i[i]);


                    pi_X_r[i] += 3.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_Y_r[i] += 3.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_Z_r[i] += 3.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_X_i[i] += 3.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_Y_i[i] += 3.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],gam_bc_r[i],gam_bc_i[i]);

                    pi_Z_i[i] += 3.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],gam_cd_r[i],gam_cd_i[i])
                               + 3.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],gam_bd_r[i],gam_bd_i[i])
                               + 3.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],gam_bc_r[i],gam_bc_i[i]);


                    pi_XX_r[i] += 3.0 * prod2_r(gam_cd_xx_r[i],gam_cd_xx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_xx_r[i],gam_bd_xx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_xx_r[i],gam_bc_xx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_XY_r[i] += 3.0 * prod2_r(gam_cd_xy_r[i],gam_cd_xy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_xy_r[i],gam_bd_xy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_xy_r[i],gam_bc_xy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_XZ_r[i] += 3.0 * prod2_r(gam_cd_xz_r[i],gam_cd_xz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_xz_r[i],gam_bd_xz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_xz_r[i],gam_bc_xz_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YX_r[i] += 3.0 * prod2_r(gam_cd_yx_r[i],gam_cd_yx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_yx_r[i],gam_bd_yx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_yx_r[i],gam_bc_yx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YY_r[i] += 3.0 * prod2_r(gam_cd_yy_r[i],gam_cd_yy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_yy_r[i],gam_bd_yy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_yy_r[i],gam_bc_yy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YZ_r[i] += 3.0 * prod2_r(gam_cd_yz_r[i],gam_cd_yz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_yz_r[i],gam_bd_yz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_yz_r[i],gam_bc_yz_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZX_r[i] += 3.0 * prod2_r(gam_cd_zx_r[i],gam_cd_zx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_zx_r[i],gam_bd_zx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_zx_r[i],gam_bc_zx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZY_r[i] += 3.0 * prod2_r(gam_cd_zy_r[i],gam_cd_zy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_zy_r[i],gam_bd_zy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_zy_r[i],gam_bc_zy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZZ_r[i] += 3.0 * prod2_r(gam_cd_zz_r[i],gam_cd_zz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_r(gam_bd_zz_r[i],gam_bd_zz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_r(gam_bc_zz_r[i],gam_bc_zz_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_XX_i[i] += 3.0 * prod2_i(gam_cd_xx_r[i],gam_cd_xx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_xx_r[i],gam_bd_xx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_xx_r[i],gam_bc_xx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_XY_i[i] += 3.0 * prod2_i(gam_cd_xy_r[i],gam_cd_xy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_xy_r[i],gam_bd_xy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_xy_r[i],gam_bc_xy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_XZ_i[i] += 3.0 * prod2_i(gam_cd_xz_r[i],gam_cd_xz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_xz_r[i],gam_bd_xz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_xz_r[i],gam_bc_xz_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YX_i[i] += 3.0 * prod2_i(gam_cd_yx_r[i],gam_cd_yx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_yx_r[i],gam_bd_yx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_yx_r[i],gam_bc_yx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YY_i[i] += 3.0 * prod2_i(gam_cd_yy_r[i],gam_cd_yy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_yy_r[i],gam_bd_yy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_yy_r[i],gam_bc_yy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_YZ_i[i] += 3.0 * prod2_i(gam_cd_yz_r[i],gam_cd_yz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_yz_r[i],gam_bd_yz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_yz_r[i],gam_bc_yz_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZX_i[i] += 3.0 * prod2_i(gam_cd_zx_r[i],gam_cd_zx_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_zx_r[i],gam_bd_zx_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_zx_r[i],gam_bc_zx_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZY_i[i] += 3.0 * prod2_i(gam_cd_zy_r[i],gam_cd_zy_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_zy_r[i],gam_bd_zy_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_zy_r[i],gam_bc_zy_i[i],rhoD_r[i],rhoD_i[i]);

                    pi_ZZ_i[i] += 3.0 * prod2_i(gam_cd_zz_r[i],gam_cd_zz_i[i],rhoB_r[i],rhoB_i[i])
                               +  3.0 * prod2_i(gam_bd_zz_r[i],gam_bd_zz_i[i],rhoC_r[i],rhoC_i[i])
                               +  3.0 * prod2_i(gam_bc_zz_r[i],gam_bc_zz_i[i],rhoD_r[i],rhoD_i[i]);


                    pi_XXX_r[i] += prod2_r(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XXY_r[i] += prod2_r(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XXZ_r[i] += prod2_r(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_XYX_r[i] += prod2_r(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XYY_r[i] += prod2_r(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XYZ_r[i] += prod2_r(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_XZX_r[i] += prod2_r(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XZY_r[i] += prod2_r(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XZZ_r[i] += prod2_r(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YXX_r[i] += prod2_r(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YXY_r[i] += prod2_r(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YXZ_r[i] += prod2_r(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YYX_r[i] += prod2_r(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YYY_r[i] += prod2_r(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YYZ_r[i] += prod2_r(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YZX_r[i] += prod2_r(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YZY_r[i] += prod2_r(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YZZ_r[i] += prod2_r(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZXX_r[i] += prod2_r(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZXY_r[i] += prod2_r(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZXZ_r[i] += prod2_r(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZYX_r[i] += prod2_r(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZYY_r[i] += prod2_r(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZYZ_r[i] += prod2_r(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZZX_r[i] += prod2_r(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZZY_r[i] += prod2_r(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZZZ_r[i] += prod2_r(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_XXX_i[i] += prod2_i(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XXY_i[i] += prod2_i(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XXZ_i[i] += prod2_i(gam_cd_xx_r[i],gam_cd_xx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_xx_r[i],gam_bd_xx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_xx_r[i],gam_bc_xx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_XYX_i[i] += prod2_i(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XYY_i[i] += prod2_i(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XYZ_i[i] += prod2_i(gam_cd_xy_r[i],gam_cd_xy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_xy_r[i],gam_bd_xy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_xy_r[i],gam_bc_xy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_XZX_i[i] += prod2_i(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_XZY_i[i] += prod2_i(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_XZZ_i[i] += prod2_i(gam_cd_xz_r[i],gam_cd_xz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_xz_r[i],gam_bd_xz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_xz_r[i],gam_bc_xz_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YXX_i[i] += prod2_i(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YXY_i[i] += prod2_i(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YXZ_i[i] += prod2_i(gam_cd_yx_r[i],gam_cd_yx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_yx_r[i],gam_bd_yx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_yx_r[i],gam_bc_yx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YYX_i[i] += prod2_i(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YYY_i[i] += prod2_i(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YYZ_i[i] += prod2_i(gam_cd_yy_r[i],gam_cd_yy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_yy_r[i],gam_bd_yy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_yy_r[i],gam_bc_yy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_YZX_i[i] += prod2_i(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_YZY_i[i] += prod2_i(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_YZZ_i[i] += prod2_i(gam_cd_yz_r[i],gam_cd_yz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_yz_r[i],gam_bd_yz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_yz_r[i],gam_bc_yz_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZXX_i[i] += prod2_i(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZXY_i[i] += prod2_i(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZXZ_i[i] += prod2_i(gam_cd_zx_r[i],gam_cd_zx_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_zx_r[i],gam_bd_zx_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_zx_r[i],gam_bc_zx_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZYX_i[i] += prod2_i(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZYY_i[i] += prod2_i(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZYZ_i[i] += prod2_i(gam_cd_zy_r[i],gam_cd_zy_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_zy_r[i],gam_bd_zy_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_zy_r[i],gam_bc_zy_i[i],gradD_z_r[i],gradD_z_i[i]);

                    pi_ZZX_i[i] += prod2_i(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_i(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_x_r[i],gradD_x_i[i]);

                    pi_ZZY_i[i] += prod2_i(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_i(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_y_r[i],gradD_y_i[i]);

                    pi_ZZZ_i[i] += prod2_i(gam_cd_zz_r[i],gam_cd_zz_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_i(gam_bd_zz_r[i],gam_bd_zz_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gam_bc_zz_r[i],gam_bc_zz_i[i],gradD_z_r[i],gradD_z_i[i]);


                    gam_r[i] += 3.0 * prod2_r(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 3.0 * prod2_r(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 3.0 * prod2_r(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_i[i] += 3.0 * prod2_i(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 3.0 * prod2_i(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 3.0 * prod2_i(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_X_r[i] += 3.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_Y_r[i] += 3.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_Z_r[i] += 3.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_X_i[i] += 3.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_Y_i[i] += 3.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_Z_i[i] += 3.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 3.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 3.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 3.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 3.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 3.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);

                    gam_XX_r[i] += 3.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XY_r[i] += 3.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_XZ_r[i] += 3.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_YX_r[i] += 3.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_YY_r[i] += 3.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YZ_r[i] += 3.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZX_r[i] += 3.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_ZY_r[i] += 3.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_ZZ_r[i] += 3.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_XX_i[i] += 3.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XY_i[i] += 3.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_XZ_i[i] += 3.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_YX_i[i] += 3.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_YY_i[i] += 3.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YZ_i[i] += 3.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZX_i[i] += 3.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + 3.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + 3.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_ZY_i[i] += 3.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + 3.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + 3.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_ZZ_i[i] += 3.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + 3.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + 3.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);
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
