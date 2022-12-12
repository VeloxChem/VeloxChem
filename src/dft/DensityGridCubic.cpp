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

    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 53 : 5;

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
CDensityGridCubic::rBrCrD(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGridCubic::rBrCrD(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGridCubic::RhoBCRhoD(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::RhoBCRhoD(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::rBrCrDXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::rBrCrDXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::rBrCrDYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridCubic::rBrCrDZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridCubic::rBrCrDZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::rBrCrDXXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDXZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDXZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(24 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(25 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(26 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(27 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(28 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(29 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(30 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDYZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDYZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(31 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(32 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(33 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(34 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(35 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(36 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(37 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(38 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(39 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
const double*
CDensityGridCubic::rBrCrDZZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}
double*
CDensityGridCubic::rBrCrDZZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(40 * _nDensityMatrices + iDensityMatrix); 

    return nullptr;  
}


const double*
CDensityGridCubic::RhoBCRhoDXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(41 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(42 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(43 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(44 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(45 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(46 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(47 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(48 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(49 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridCubic::RhoBCRhoDX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(50 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(51 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
const double*
CDensityGridCubic::RhoBCRhoDZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
double*
CDensityGridCubic::RhoBCRhoDZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(52 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

void
CDensityGridCubic::DensityProd(const CDensityGrid& rwDensityGrid,
                              const CDensityGrid& rw2DensityGrid,
                              const xcfun         xcFuncType,
                              const int32_t             numdens,
                              const std::string&  CubeMode) const

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
                
                auto RhoBCRhoDx_r = RhoBCRhoD(6 * j);

                auto RhoBCRhoDx_i = RhoBCRhoD(6 * j + 1);

                auto RhoBCRhoDy_r = RhoBCRhoD(6 * j + 2);

                auto RhoBCRhoDy_i = RhoBCRhoD(6 * j + 3);

                auto RhoBCRhoDz_r = RhoBCRhoD(6 * j + 4);

                auto RhoBCRhoDz_i = RhoBCRhoD(6 * j + 5);

                auto rBrCrDx_r = rBrCrD(6 * j);

                auto rBrCrDx_i = rBrCrD(6 * j + 1);

                auto rBrCrDy_r = rBrCrD(6 * j + 2 );

                auto rBrCrDy_i = rBrCrD(6 * j + 3);

                auto rBrCrDz_r = rBrCrD(6 * j + 4);

                auto rBrCrDz_i = rBrCrD(6 * j + 5);


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


                    rBrCrDx_r[i] = 12.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 12.0 * (rhobx_r * rhocy_r * rhoby_r - rhobx_i * rhocy_r * rhoby_i - rhobx_r * rhocy_i * rhoby_i - rhobx_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocx_r * rhoby_r - rhoby_i * rhocx_r * rhoby_i - rhoby_r * rhocx_i * rhoby_i - rhoby_i * rhocx_i * rhoby_r)

                    + 12.0 * (rhobx_r * rhocz_r * rhobz_r - rhobx_i * rhocz_r * rhobz_i - rhobx_r * rhocz_i * rhobz_i - rhobx_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocx_r * rhobz_r - rhobz_i * rhocx_r * rhobz_i - rhobz_r * rhocx_i * rhobz_i - rhobz_i * rhocx_i * rhobz_r);


                    rBrCrDy_r[i] = 12.0 * (rhoby_r * rhocx_r * rhobx_r - rhoby_i * rhocx_r * rhobx_i - rhoby_r * rhocx_i * rhobx_i - rhoby_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocy_r * rhobx_r - rhobx_i * rhocy_r * rhobx_i - rhobx_r * rhocy_i * rhobx_i - rhobx_i * rhocy_i * rhobx_r)

                    + 12.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 12.0 * (rhoby_r * rhocz_r * rhobz_r - rhoby_i * rhocz_r * rhobz_i - rhoby_r * rhocz_i * rhobz_i - rhoby_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocy_r * rhobz_r - rhobz_i * rhocy_r * rhobz_i - rhobz_r * rhocy_i * rhobz_i - rhobz_i * rhocy_i * rhobz_r);


                    rBrCrDz_r[i] = 12.0 * (rhobz_r * rhocx_r * rhobx_r - rhobz_i * rhocx_r * rhobx_i - rhobz_r * rhocx_i * rhobx_i - rhobz_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocz_r * rhobx_r - rhobx_i * rhocz_r * rhobx_i - rhobx_r * rhocz_i * rhobx_i - rhobx_i * rhocz_i * rhobx_r)

                    + 12.0 * (rhobz_r * rhocy_r * rhoby_r - rhobz_i * rhocy_r * rhoby_i - rhobz_r * rhocy_i * rhoby_i - rhobz_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocz_r * rhoby_r - rhoby_i * rhocz_r * rhoby_i - rhoby_r * rhocz_i * rhoby_i - rhoby_i * rhocz_i * rhoby_r)

                    + 12.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r);



                    rBrCrDx_i[i] = 12.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 12.0 * ( - rhobx_i * rhocy_i * rhoby_i + rhobx_i * rhocy_r * rhoby_r + rhobx_r * rhocy_i * rhoby_r + rhobx_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocx_i * rhoby_i + rhoby_i * rhocx_r * rhoby_r + rhoby_r * rhocx_i * rhoby_r + rhoby_r * rhocx_r * rhoby_i)

                    + 12.0 * ( - rhobx_i * rhocz_i * rhobz_i + rhobx_i * rhocz_r * rhobz_r + rhobx_r * rhocz_i * rhobz_r + rhobx_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocx_i * rhobz_i + rhobz_i * rhocx_r * rhobz_r + rhobz_r * rhocx_i * rhobz_r + rhobz_r * rhocx_r * rhobz_i);


                    rBrCrDy_i[i] = 12.0 * ( - rhoby_i * rhocx_i * rhobx_i + rhoby_i * rhocx_r * rhobx_r + rhoby_r * rhocx_i * rhobx_r + rhoby_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocy_i * rhobx_i + rhobx_i * rhocy_r * rhobx_r + rhobx_r * rhocy_i * rhobx_r + rhobx_r * rhocy_r * rhobx_i)

                    + 12.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 12.0 * ( - rhoby_i * rhocz_i * rhobz_i + rhoby_i * rhocz_r * rhobz_r + rhoby_r * rhocz_i * rhobz_r + rhoby_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocy_i * rhobz_i + rhobz_i * rhocy_r * rhobz_r + rhobz_r * rhocy_i * rhobz_r + rhobz_r * rhocy_r * rhobz_i);


                    rBrCrDz_i[i] = 12.0 * ( - rhobz_i * rhocx_i * rhobx_i + rhobz_i * rhocx_r * rhobx_r + rhobz_r * rhocx_i * rhobx_r + rhobz_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocz_i * rhobx_i + rhobx_i * rhocz_r * rhobx_r + rhobx_r * rhocz_i * rhobx_r + rhobx_r * rhocz_r * rhobx_i)

                    + 12.0 * ( - rhobz_i * rhocy_i * rhoby_i + rhobz_i * rhocy_r * rhoby_r + rhobz_r * rhocy_i * rhoby_r + rhobz_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocz_i * rhoby_i + rhoby_i * rhocz_r * rhoby_r + rhoby_r * rhocz_i * rhoby_r + rhoby_r * rhocz_r * rhoby_i)

                    + 12.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i);


                    RhoBCRhoDx_r[i] = (d_sig_xx_r* rhocx_r - d_sig_xx_i * rhocx_i 
                    
                                    + d_lamtau_xx_r* rhobx_r - d_lamtau_xx_i * rhobx_i 
                                    
                                    + d_sig_xy_r* rhocy_r - d_sig_xy_i * rhocy_i 

                                     + d_lamtau_xy_r* rhoby_r - d_lamtau_xy_i * rhoby_i 

                                    + d_sig_xz_r* rhocz_r - d_sig_xz_i * rhocz_i 
                                    
                                    + d_lamtau_xz_r* rhobz_r - d_lamtau_xz_i * rhobz_i);


                    RhoBCRhoDy_r[i] = (d_sig_yx_r* rhocx_r - d_sig_yx_i * rhocx_i + d_lamtau_yx_r* rhobx_r - d_lamtau_yx_i * rhobx_i +

                    d_sig_yy_r* rhocy_r - d_sig_yy_i * rhocy_i + d_lamtau_yy_r* rhoby_r - d_lamtau_yy_i * rhoby_i +

                    d_sig_yz_r* rhocz_r - d_sig_yz_i * rhocz_i + d_lamtau_yz_r* rhobz_r - d_lamtau_yz_i * rhobz_i);


                    RhoBCRhoDz_r[i] = (d_sig_zx_r* rhocx_r - d_sig_zx_i * rhocx_i + d_lamtau_zx_r* rhobx_r - d_lamtau_zx_i * rhobx_i +

                    d_sig_zy_r* rhocy_r - d_sig_zy_i * rhocy_i + d_lamtau_zy_r* rhoby_r - d_lamtau_zy_i * rhoby_i +

                    d_sig_zz_r* rhocz_r - d_sig_zz_i * rhocz_i + d_lamtau_zz_r* rhobz_r - d_lamtau_zz_i * rhobz_i);


                    RhoBCRhoDx_i[i] =  (d_sig_xx_i* rhocx_r + d_sig_xx_r * rhocx_i + d_lamtau_xx_i* rhobx_r + d_lamtau_xx_r * rhobx_i +

                    d_sig_xy_i* rhocy_r + d_sig_xy_r * rhocy_i + d_lamtau_xy_i* rhoby_r + d_lamtau_xy_r * rhoby_i +

                    d_sig_xz_i* rhocz_r + d_sig_xz_r * rhocz_i + d_lamtau_xz_i* rhobz_r + d_lamtau_xz_r * rhobz_i);


                    RhoBCRhoDy_i[i] =  (d_sig_yx_i* rhocx_r + d_sig_yx_r * rhocx_i + d_lamtau_yx_i* rhobx_r + d_lamtau_yx_r * rhobx_i +

                    d_sig_yy_i* rhocy_r + d_sig_yy_r * rhocy_i + d_lamtau_yy_i* rhoby_r + d_lamtau_yy_r * rhoby_i +

                    d_sig_yz_i* rhocz_r + d_sig_yz_r * rhocz_i + d_lamtau_yz_i* rhobz_r + d_lamtau_yz_r * rhobz_i);


                    RhoBCRhoDz_i[i] =  (d_sig_zx_i* rhocx_r + d_sig_zx_r * rhocx_i + d_lamtau_zx_i* rhobx_r + d_lamtau_zx_r * rhobx_i +

                    d_sig_zy_i* rhocy_r + d_sig_zy_r * rhocy_i + d_lamtau_zy_i* rhoby_r + d_lamtau_zy_r * rhoby_i +

                    d_sig_zz_i* rhocz_r + d_sig_zz_r * rhocz_i + d_lamtau_zz_i* rhobz_r + d_lamtau_zz_r * rhobz_i);
                    
                }

            }

        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                
                auto RhoBCRhoD_r = RhoBCRhoD(2 * j);

                auto RhoBCRhoD_i = RhoBCRhoD(2 * j + 1);

                auto rBrCrD_r = rBrCrD(2 * j);

                auto rBrCrD_i = rBrCrD(2 * j + 1);

                
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

                // rBrCrD = rho_b * rho_c * rho_d

                for (int32_t i = 0; i < npoints; i++)
                {
                    rBrCrD_r[i] = 6.0 * (rhoB_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * rhoD_r[i]);

                    rBrCrD_i[i] = 6.0 *  (- rhoB_i[i] * rhoC_i[i] * rhoD_i[i] 
                    
                                        + rhoB_i[i] * rhoC_r[i] * rhoD_r[i] 
                                        
                                        + rhoC_i[i] * rhoB_r[i] * rhoD_r[i] 
                                        
                                        + rhoD_i[i] * rhoB_r[i] * rhoC_r[i]);

                    RhoBCRhoD_r[i] = 3.0 * (rhoD_r[i] * rhoBC_r[i] - rhoD_i[i] * rhoBC_i[i] 
                    
                                         + rhoC_r[i] * rhoBD_r[i] - rhoC_i[i] * rhoBD_i[i] 
                                         
                                         + rhoB_r[i] * rhoCD_r[i] - rhoB_i[i] * rhoCD_i[i]) ;

                    RhoBCRhoD_i[i] = 3.0 * (rhoD_i[i] * rhoBC_r[i] + rhoD_r[i] * rhoBC_i[i]
                    
                                         + rhoC_i[i] * rhoBD_r[i] + rhoC_r[i] * rhoBD_i[i] 
                                         
                                        + rhoB_i[i] * rhoCD_r[i] + rhoB_r[i] * rhoCD_i[i]);

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
                
                auto RhoBCRhoDx_r = RhoBCRhoD(6 * j);

                auto RhoBCRhoDx_i = RhoBCRhoD(6 * j + 1);

                auto RhoBCRhoDy_r = RhoBCRhoD(6 * j + 2);

                auto RhoBCRhoDy_i = RhoBCRhoD(6 * j + 3);

                auto RhoBCRhoDz_r = RhoBCRhoD(6 * j + 4);

                auto RhoBCRhoDz_i = RhoBCRhoD(6 * j + 5);

                auto rBrCrDx_r = rBrCrD(6 * j);

                auto rBrCrDx_i = rBrCrD(6 * j + 1);

                auto rBrCrDy_r = rBrCrD(6 * j + 2 );

                auto rBrCrDy_i = rBrCrD(6 * j + 3);

                auto rBrCrDz_r = rBrCrD(6 * j + 4);

                auto rBrCrDz_i = rBrCrD(6 * j + 5);

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


                auto D_sig_xz_r = rw2DensityGrid.alphaDensity(24 * j + 8);

                auto D_sig_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 8);

                auto D_sig_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 8);

                auto D_sig_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 8);


                auto D_sig_xz_i = rw2DensityGrid.alphaDensity(24 * j + 9);

                auto D_sig_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 9);

                auto D_sig_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 9);

                auto D_sig_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 9);


                auto D_sig_yz_r = rw2DensityGrid.alphaDensity(24 * j + 10);

                auto D_sig_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 10);

                auto D_sig_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 10);

                auto D_sig_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 10);


                auto D_sig_yz_i = rw2DensityGrid.alphaDensity(24 * j + 11);

                auto D_sig_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 11);

                auto D_sig_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 11);

                auto D_sig_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 11);


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


                auto D_lamtau_xz_r = rw2DensityGrid.alphaDensity(24 * j + 20);

                auto D_lamtau_xz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 20);

                auto D_lamtau_xz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 20);

                auto D_lamtau_xz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 20);


                auto D_lamtau_xz_i = rw2DensityGrid.alphaDensity(24 * j + 21);

                auto D_lamtau_xz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 21);

                auto D_lamtau_xz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 21);

                auto D_lamtau_xz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 21);


                auto D_lamtau_yz_r = rw2DensityGrid.alphaDensity(24 * j + 22);

                auto D_lamtau_yz_x_r = rw2DensityGrid.alphaDensityGradientX(24 * j + 22);

                auto D_lamtau_yz_y_r = rw2DensityGrid.alphaDensityGradientY(24 * j + 22);

                auto D_lamtau_yz_z_r = rw2DensityGrid.alphaDensityGradientZ(24 * j + 22);


                auto D_lamtau_yz_i = rw2DensityGrid.alphaDensity(24 * j + 23);

                auto D_lamtau_yz_x_i = rw2DensityGrid.alphaDensityGradientX(24 * j + 23);

                auto D_lamtau_yz_y_i = rw2DensityGrid.alphaDensityGradientY(24 * j + 23);

                auto D_lamtau_yz_z_i = rw2DensityGrid.alphaDensityGradientZ(24 * j + 23);



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


                    rBrCrDx_r[i] = 12.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocx_r * rhobx_r - rhobx_i * rhocx_r * rhobx_i - rhobx_r * rhocx_i * rhobx_i - rhobx_i * rhocx_i * rhobx_r)

                    + 12.0 * (rhobx_r * rhocy_r * rhoby_r - rhobx_i * rhocy_r * rhoby_i - rhobx_r * rhocy_i * rhoby_i - rhobx_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocx_r * rhoby_r - rhoby_i * rhocx_r * rhoby_i - rhoby_r * rhocx_i * rhoby_i - rhoby_i * rhocx_i * rhoby_r)

                    + 12.0 * (rhobx_r * rhocz_r * rhobz_r - rhobx_i * rhocz_r * rhobz_i - rhobx_r * rhocz_i * rhobz_i - rhobx_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocx_r * rhobz_r - rhobz_i * rhocx_r * rhobz_i - rhobz_r * rhocx_i * rhobz_i - rhobz_i * rhocx_i * rhobz_r);


                    rBrCrDy_r[i] = 12.0 * (rhoby_r * rhocx_r * rhobx_r - rhoby_i * rhocx_r * rhobx_i - rhoby_r * rhocx_i * rhobx_i - rhoby_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocy_r * rhobx_r - rhobx_i * rhocy_r * rhobx_i - rhobx_r * rhocy_i * rhobx_i - rhobx_i * rhocy_i * rhobx_r)

                    + 12.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocy_r * rhoby_r - rhoby_i * rhocy_r * rhoby_i - rhoby_r * rhocy_i * rhoby_i - rhoby_i * rhocy_i * rhoby_r)

                    + 12.0 * (rhoby_r * rhocz_r * rhobz_r - rhoby_i * rhocz_r * rhobz_i - rhoby_r * rhocz_i * rhobz_i - rhoby_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocy_r * rhobz_r - rhobz_i * rhocy_r * rhobz_i - rhobz_r * rhocy_i * rhobz_i - rhobz_i * rhocy_i * rhobz_r);


                    rBrCrDz_r[i] = 12.0 * (rhobz_r * rhocx_r * rhobx_r - rhobz_i * rhocx_r * rhobx_i - rhobz_r * rhocx_i * rhobx_i - rhobz_i * rhocx_i * rhobx_r)

                    + 6.0 * (rhobx_r * rhocz_r * rhobx_r - rhobx_i * rhocz_r * rhobx_i - rhobx_r * rhocz_i * rhobx_i - rhobx_i * rhocz_i * rhobx_r)

                    + 12.0 * (rhobz_r * rhocy_r * rhoby_r - rhobz_i * rhocy_r * rhoby_i - rhobz_r * rhocy_i * rhoby_i - rhobz_i * rhocy_i * rhoby_r)

                    + 6.0 * (rhoby_r * rhocz_r * rhoby_r - rhoby_i * rhocz_r * rhoby_i - rhoby_r * rhocz_i * rhoby_i - rhoby_i * rhocz_i * rhoby_r)

                    + 12.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r)

                    + 6.0 * (rhobz_r * rhocz_r * rhobz_r - rhobz_i * rhocz_r * rhobz_i - rhobz_r * rhocz_i * rhobz_i - rhobz_i * rhocz_i * rhobz_r);


                    rBrCrDx_i[i] = 12.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocx_i * rhobx_i + rhobx_i * rhocx_r * rhobx_r + rhobx_r * rhocx_i * rhobx_r + rhobx_r * rhocx_r * rhobx_i)

                    + 12.0 * ( - rhobx_i * rhocy_i * rhoby_i + rhobx_i * rhocy_r * rhoby_r + rhobx_r * rhocy_i * rhoby_r + rhobx_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocx_i * rhoby_i + rhoby_i * rhocx_r * rhoby_r + rhoby_r * rhocx_i * rhoby_r + rhoby_r * rhocx_r * rhoby_i)

                    + 12.0 * ( - rhobx_i * rhocz_i * rhobz_i + rhobx_i * rhocz_r * rhobz_r + rhobx_r * rhocz_i * rhobz_r + rhobx_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocx_i * rhobz_i + rhobz_i * rhocx_r * rhobz_r + rhobz_r * rhocx_i * rhobz_r + rhobz_r * rhocx_r * rhobz_i);


                    rBrCrDy_i[i] = 12.0 * ( - rhoby_i * rhocx_i * rhobx_i + rhoby_i * rhocx_r * rhobx_r + rhoby_r * rhocx_i * rhobx_r + rhoby_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocy_i * rhobx_i + rhobx_i * rhocy_r * rhobx_r + rhobx_r * rhocy_i * rhobx_r + rhobx_r * rhocy_r * rhobx_i)

                    + 12.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocy_i * rhoby_i + rhoby_i * rhocy_r * rhoby_r + rhoby_r * rhocy_i * rhoby_r + rhoby_r * rhocy_r * rhoby_i)

                    + 12.0 * ( - rhoby_i * rhocz_i * rhobz_i + rhoby_i * rhocz_r * rhobz_r + rhoby_r * rhocz_i * rhobz_r + rhoby_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocy_i * rhobz_i + rhobz_i * rhocy_r * rhobz_r + rhobz_r * rhocy_i * rhobz_r + rhobz_r * rhocy_r * rhobz_i);


                    rBrCrDz_i[i] = 12.0 * ( - rhobz_i * rhocx_i * rhobx_i + rhobz_i * rhocx_r * rhobx_r + rhobz_r * rhocx_i * rhobx_r + rhobz_r * rhocx_r * rhobx_i)

                    + 6.0 * ( - rhobx_i * rhocz_i * rhobx_i + rhobx_i * rhocz_r * rhobx_r + rhobx_r * rhocz_i * rhobx_r + rhobx_r * rhocz_r * rhobx_i)

                    + 12.0 * ( - rhobz_i * rhocy_i * rhoby_i + rhobz_i * rhocy_r * rhoby_r + rhobz_r * rhocy_i * rhoby_r + rhobz_r * rhocy_r * rhoby_i)

                    + 6.0 * ( - rhoby_i * rhocz_i * rhoby_i + rhoby_i * rhocz_r * rhoby_r + rhoby_r * rhocz_i * rhoby_r + rhoby_r * rhocz_r * rhoby_i)

                    + 12.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i)

                    + 6.0 * ( - rhobz_i * rhocz_i * rhobz_i + rhobz_i * rhocz_r * rhobz_r + rhobz_r * rhocz_i * rhobz_r + rhobz_r * rhocz_r * rhobz_i);


                    RhoBCRhoDx_r[i] = (d_sig_xx_r* rhocx_r - d_sig_xx_i * rhocx_i 
                    
                                    + d_lamtau_xx_r* rhobx_r - d_lamtau_xx_i * rhobx_i 
                                    
                                    + d_sig_xy_r* rhocy_r - d_sig_xy_i * rhocy_i 

                                     + d_lamtau_xy_r* rhoby_r - d_lamtau_xy_i * rhoby_i 

                                    + d_sig_xz_r* rhocz_r - d_sig_xz_i * rhocz_i 
                                    
                                    + d_lamtau_xz_r* rhobz_r - d_lamtau_xz_i * rhobz_i);


                    RhoBCRhoDy_r[i] = (d_sig_yx_r* rhocx_r - d_sig_yx_i * rhocx_i + d_lamtau_yx_r* rhobx_r - d_lamtau_yx_i * rhobx_i +

                    d_sig_yy_r* rhocy_r - d_sig_yy_i * rhocy_i + d_lamtau_yy_r* rhoby_r - d_lamtau_yy_i * rhoby_i +

                    d_sig_yz_r* rhocz_r - d_sig_yz_i * rhocz_i + d_lamtau_yz_r* rhobz_r - d_lamtau_yz_i * rhobz_i);


                    RhoBCRhoDz_r[i] = (d_sig_zx_r* rhocx_r - d_sig_zx_i * rhocx_i + d_lamtau_zx_r* rhobx_r - d_lamtau_zx_i * rhobx_i +

                    d_sig_zy_r* rhocy_r - d_sig_zy_i * rhocy_i + d_lamtau_zy_r* rhoby_r - d_lamtau_zy_i * rhoby_i +

                    d_sig_zz_r* rhocz_r - d_sig_zz_i * rhocz_i + d_lamtau_zz_r* rhobz_r - d_lamtau_zz_i * rhobz_i);


                    RhoBCRhoDx_i[i] =  (d_sig_xx_i* rhocx_r + d_sig_xx_r * rhocx_i + d_lamtau_xx_i* rhobx_r + d_lamtau_xx_r * rhobx_i +

                    d_sig_xy_i* rhocy_r + d_sig_xy_r * rhocy_i + d_lamtau_xy_i* rhoby_r + d_lamtau_xy_r * rhoby_i +

                    d_sig_xz_i* rhocz_r + d_sig_xz_r * rhocz_i + d_lamtau_xz_i* rhobz_r + d_lamtau_xz_r * rhobz_i);


                    RhoBCRhoDy_i[i] =  (d_sig_yx_i* rhocx_r + d_sig_yx_r * rhocx_i + d_lamtau_yx_i* rhobx_r + d_lamtau_yx_r * rhobx_i +

                    d_sig_yy_i* rhocy_r + d_sig_yy_r * rhocy_i + d_lamtau_yy_i* rhoby_r + d_lamtau_yy_r * rhoby_i +

                    d_sig_yz_i* rhocz_r + d_sig_yz_r * rhocz_i + d_lamtau_yz_i* rhobz_r + d_lamtau_yz_r * rhobz_i);


                    RhoBCRhoDz_i[i] =  (d_sig_zx_i* rhocx_r + d_sig_zx_r * rhocx_i + d_lamtau_zx_i* rhobx_r + d_lamtau_zx_r * rhobx_i +

                    d_sig_zy_i* rhocy_r + d_sig_zy_r * rhocy_i + d_lamtau_zy_i* rhoby_r + d_lamtau_zy_r * rhoby_i +

                    d_sig_zz_i* rhocz_r + d_sig_zz_r * rhocz_i + d_lamtau_zz_i* rhobz_r + d_lamtau_zz_r * rhobz_i);
                    
                }
            }
        }
        if (fstr::upcase(CubeMode) == "CRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                
                auto RhoBCRhoD_r = RhoBCRhoD(2 * j);

                auto RhoBCRhoD_i = RhoBCRhoD(2 * j + 1);

                auto RhoBCRhoD_X_r = RhoBCRhoDX(2 * j);

                auto RhoBCRhoD_X_i = RhoBCRhoDX(2 * j + 1);

                auto RhoBCRhoD_Y_r = RhoBCRhoDY(2 * j);

                auto RhoBCRhoD_Y_i = RhoBCRhoDY(2 * j + 1);

                auto RhoBCRhoD_Z_r = RhoBCRhoDZ(2 * j);

                auto RhoBCRhoD_Z_i = RhoBCRhoDZ(2 * j + 1);

                auto RhoBCRhoD_XX_r = RhoBCRhoDXX(2 * j);

                auto RhoBCRhoD_XX_i = RhoBCRhoDXX(2 * j + 1);

                auto RhoBCRhoD_XY_r = RhoBCRhoDXY(2 * j);

                auto RhoBCRhoD_XY_i = RhoBCRhoDXY(2 * j + 1);

                auto RhoBCRhoD_XZ_r = RhoBCRhoDXZ(2 * j);

                auto RhoBCRhoD_XZ_i = RhoBCRhoDXZ(2 * j + 1);

                auto RhoBCRhoD_YX_r = RhoBCRhoDYX(2 * j);

                auto RhoBCRhoD_YX_i = RhoBCRhoDYX(2 * j + 1);

                auto RhoBCRhoD_YY_r = RhoBCRhoDYY(2 * j);

                auto RhoBCRhoD_YY_i = RhoBCRhoDYY(2 * j + 1);

                auto RhoBCRhoD_YZ_r = RhoBCRhoDYZ(2 * j);

                auto RhoBCRhoD_YZ_i = RhoBCRhoDYZ(2 * j + 1);

                auto RhoBCRhoD_ZX_r = RhoBCRhoDZX(2 * j);

                auto RhoBCRhoD_ZX_i = RhoBCRhoDZX(2 * j + 1);

                auto RhoBCRhoD_ZY_r = RhoBCRhoDZY(2 * j);

                auto RhoBCRhoD_ZY_i = RhoBCRhoDZY(2 * j + 1);

                auto RhoBCRhoD_ZZ_r = RhoBCRhoDZZ(2 * j);

                auto RhoBCRhoD_ZZ_i = RhoBCRhoDZZ(2 * j + 1);


                auto rBrCrD_r = rBrCrD(2 * j);

                auto rBrCrD_i = rBrCrD(2 * j + 1);
                

                auto rBrCrDX_r = rBrCrDX(2 * j);

                auto rBrCrDX_i = rBrCrDX(2 * j + 1);

                auto rBrCrDY_r = rBrCrDY(2 * j);

                auto rBrCrDY_i = rBrCrDY(2 * j + 1);
                
                auto rBrCrDZ_r = rBrCrDZ(2 * j);

                auto rBrCrDZ_i = rBrCrDZ(2 * j + 1);


                auto rBrCrDXX_r = rBrCrDXX(2 * j);

                auto rBrCrDXX_i = rBrCrDXX(2 * j + 1);

                auto rBrCrDXY_r = rBrCrDXY(2 * j);

                auto rBrCrDXY_i = rBrCrDXY(2 * j + 1);
                
                auto rBrCrDXZ_r = rBrCrDXZ(2 * j);

                auto rBrCrDXZ_i = rBrCrDXZ(2 * j + 1);


                auto rBrCrDYX_r = rBrCrDYX(2 * j);

                auto rBrCrDYX_i = rBrCrDYX(2 * j + 1);

                auto rBrCrDYY_r = rBrCrDYY(2 * j);

                auto rBrCrDYY_i = rBrCrDYY(2 * j + 1);
                
                auto rBrCrDYZ_r = rBrCrDYZ(2 * j);

                auto rBrCrDYZ_i = rBrCrDYZ(2 * j + 1);


                auto rBrCrDZX_r = rBrCrDZX(2 * j);

                auto rBrCrDZX_i = rBrCrDZX(2 * j + 1);

                auto rBrCrDZY_r = rBrCrDZY(2 * j);

                auto rBrCrDZY_i = rBrCrDZY(2 * j + 1);
                
                auto rBrCrDZZ_r = rBrCrDZZ(2 * j);

                auto rBrCrDZZ_i = rBrCrDZZ(2 * j + 1);

                auto rBrCrD_XXX_r=rBrCrDXXX(2 * j);

                auto rBrCrD_XXX_i=rBrCrDXXX(2 * j + 1);

                auto rBrCrD_XXY_r=rBrCrDXXY(2 * j);

                auto rBrCrD_XXY_i=rBrCrDXXY(2 * j + 1);

                auto rBrCrD_XXZ_r=rBrCrDXXZ(2 * j);

                auto rBrCrD_XXZ_i=rBrCrDXXZ(2 * j + 1);

                auto rBrCrD_XYX_r=rBrCrDXYX(2 * j);

                auto rBrCrD_XYX_i=rBrCrDXYX(2 * j + 1);

                auto rBrCrD_XYY_r=rBrCrDXYY(2 * j);

                auto rBrCrD_XYY_i=rBrCrDXYY(2 * j + 1);

                auto rBrCrD_XYZ_r=rBrCrDXYZ(2 * j);

                auto rBrCrD_XYZ_i=rBrCrDXYZ(2 * j + 1);

                auto rBrCrD_XZX_r=rBrCrDXZX(2 * j);

                auto rBrCrD_XZX_i=rBrCrDXZX(2 * j + 1);

                auto rBrCrD_XZY_r=rBrCrDXZY(2 * j);

                auto rBrCrD_XZY_i=rBrCrDXZY(2 * j + 1);

                auto rBrCrD_XZZ_r=rBrCrDXZZ(2 * j);

                auto rBrCrD_XZZ_i=rBrCrDXZZ(2 * j + 1);

                auto rBrCrD_YXX_r=rBrCrDYXX(2 * j);

                auto rBrCrD_YXX_i=rBrCrDYXX(2 * j + 1);

                auto rBrCrD_YXY_r=rBrCrDYXY(2 * j);

                auto rBrCrD_YXY_i=rBrCrDYXY(2 * j + 1);

                auto rBrCrD_YXZ_r=rBrCrDYXZ(2 * j);

                auto rBrCrD_YXZ_i=rBrCrDYXZ(2 * j + 1);

                auto rBrCrD_YYX_r=rBrCrDYYX(2 * j);

                auto rBrCrD_YYX_i=rBrCrDYYX(2 * j + 1);

                auto rBrCrD_YYY_r=rBrCrDYYY(2 * j);

                auto rBrCrD_YYY_i=rBrCrDYYY(2 * j + 1);

                auto rBrCrD_YYZ_r=rBrCrDYYZ(2 * j);

                auto rBrCrD_YYZ_i=rBrCrDYYZ(2 * j + 1);

                auto rBrCrD_YZX_r=rBrCrDYZX(2 * j);

                auto rBrCrD_YZX_i=rBrCrDYZX(2 * j + 1);

                auto rBrCrD_YZY_r=rBrCrDYZY(2 * j);

                auto rBrCrD_YZY_i=rBrCrDYZY(2 * j + 1);

                auto rBrCrD_YZZ_r=rBrCrDYZZ(2 * j);

                auto rBrCrD_YZZ_i=rBrCrDYZZ(2 * j + 1);

                auto rBrCrD_ZXX_r=rBrCrDZXX(2 * j);

                auto rBrCrD_ZXX_i=rBrCrDZXX(2 * j + 1);

                auto rBrCrD_ZXY_r=rBrCrDZXY(2 * j);

                auto rBrCrD_ZXY_i=rBrCrDZXY(2 * j + 1);

                auto rBrCrD_ZXZ_r=rBrCrDZXZ(2 * j);

                auto rBrCrD_ZXZ_i=rBrCrDZXZ(2 * j + 1);

                auto rBrCrD_ZYX_r=rBrCrDZYX(2 * j);

                auto rBrCrD_ZYX_i=rBrCrDZYX(2 * j + 1);

                auto rBrCrD_ZYY_r=rBrCrDZYY(2 * j);

                auto rBrCrD_ZYY_i=rBrCrDZYY(2 * j + 1);

                auto rBrCrD_ZYZ_r=rBrCrDZYZ(2 * j);

                auto rBrCrD_ZYZ_i=rBrCrDZYZ(2 * j + 1);

                auto rBrCrD_ZZX_r=rBrCrDZZX(2 * j);

                auto rBrCrD_ZZX_i=rBrCrDZZX(2 * j + 1);

                auto rBrCrD_ZZY_r=rBrCrDZZY(2 * j);

                auto rBrCrD_ZZY_i=rBrCrDZZY(2 * j + 1);

                auto rBrCrD_ZZZ_r=rBrCrDZZZ(2 * j);

                auto rBrCrD_ZZZ_i=rBrCrDZZZ(2 * j + 1);


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


                // rBrCrD = rho_b * rho_c * rho_d

                for (int32_t i = 0; i < npoints; i++)
                {

                    rBrCrD_r[i] = 6.0 * (rhoB_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * rhoD_r[i]);

                    rBrCrD_i[i] = 6.0 *  (- rhoB_i[i] * rhoC_i[i] * rhoD_i[i] 
                    
                                        + rhoB_i[i] * rhoC_r[i] * rhoD_r[i] 
                                        
                                        + rhoC_i[i] * rhoB_r[i] * rhoD_r[i] 
                                        
                                        + rhoD_i[i] * rhoB_r[i] * rhoC_r[i]);


                    rBrCrDX_r[i] = 6.0 * (gradB_x_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - gradB_x_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * rhoC_i[i] * rhoD_r[i]);


                    rBrCrDX_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_x_i[i] * rhoD_r[i]);


                    rBrCrDX_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_x_r[i] 
                    
                                    - rhoB_i[i] * gradD_x_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_x_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * gradD_x_r[i]);


                    rBrCrDY_r[i] = 6.0 * (gradB_y_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - gradB_y_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * rhoC_i[i] * rhoD_r[i]);


                    rBrCrDY_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_y_i[i] * rhoD_r[i]);


                    rBrCrDY_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_y_r[i] 
                    
                                    - rhoB_i[i] * gradD_y_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_y_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * gradD_y_r[i]);


                    rBrCrDZ_r[i] = 6.0 * (gradB_z_r[i] * rhoC_r[i] * rhoD_r[i] 
                    
                                    - gradB_z_i[i] * rhoD_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * rhoD_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * rhoC_i[i] * rhoD_r[i]);


                    rBrCrDZ_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * rhoD_r[i] 
                    
                                    - rhoB_i[i] * rhoD_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * rhoD_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_z_i[i] * rhoD_r[i]);


                    rBrCrDZ_r[i] += 6.0 * (rhoB_r[i] * rhoC_r[i] * gradD_z_r[i] 
                    
                                    - rhoB_i[i] * gradD_z_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_z_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * rhoC_i[i] * gradD_z_r[i]);


                    rBrCrDX_i[i] =  6.0 * ( - gradB_x_i[i] * rhoC_i[i] * rhoD_i[i] 
                                            + gradB_x_i[i] * rhoC_r[i] * rhoD_r[i]
                                            + rhoC_i[i] * gradB_x_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_x_r[i] * rhoC_r[i]);

                    rBrCrDX_i[i] +=  6.0 * (- rhoB_i[i] * gradC_x_i[i] * rhoD_i[i] 
                                            + rhoB_i[i] * gradC_x_r[i] * rhoD_r[i]
                                            + gradC_x_i[i] * rhoB_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * rhoB_r[i] * gradC_x_r[i]);

                    rBrCrDX_i[i] +=  6.0 *(- rhoB_i[i] * rhoC_i[i] * gradD_x_i[i]
                                            + rhoB_i[i] * rhoC_r[i] * gradD_x_r[i]
                                            + rhoC_i[i] * rhoB_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * rhoB_r[i] * rhoC_r[i]);

                    rBrCrDY_i[i] =  6.0 * ( - gradB_y_i[i] * rhoC_i[i] * rhoD_i[i] 
                                            + gradB_y_i[i] * rhoC_r[i] * rhoD_r[i]
                                            + rhoC_i[i] * gradB_y_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_y_r[i] * rhoC_r[i]);

                    rBrCrDY_i[i] +=  6.0 * (- rhoB_i[i] * gradC_y_i[i] * rhoD_i[i] 
                                            + rhoB_i[i] * gradC_y_r[i] * rhoD_r[i]
                                            + gradC_y_i[i] * rhoB_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * rhoB_r[i] * gradC_y_r[i]);

                    rBrCrDY_i[i] +=  6.0 *(- rhoB_i[i] * rhoC_i[i] * gradD_y_i[i]
                                            + rhoB_i[i] * rhoC_r[i] * gradD_y_r[i]
                                            + rhoC_i[i] * rhoB_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * rhoB_r[i] * rhoC_r[i]);

                    rBrCrDZ_i[i] =  6.0 * ( - gradB_z_i[i] * rhoC_i[i] * rhoD_i[i] 
                                            + gradB_z_i[i] * rhoC_r[i] * rhoD_r[i]
                                            + rhoC_i[i] * gradB_z_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_z_r[i] * rhoC_r[i]);

                    rBrCrDZ_i[i] +=  6.0 * (- rhoB_i[i] * gradC_z_i[i] * rhoD_i[i] 
                                            + rhoB_i[i] * gradC_z_r[i] * rhoD_r[i]
                                            + gradC_z_i[i] * rhoB_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * rhoB_r[i] * gradC_z_r[i]);

                    rBrCrDZ_i[i] +=  6.0 *(- rhoB_i[i] * rhoC_i[i] * gradD_z_i[i]
                                            + rhoB_i[i] * rhoC_r[i] * gradD_z_r[i]
                                            + rhoC_i[i] * rhoB_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * rhoB_r[i] * rhoC_r[i]);
                    // BC

                    rBrCrDXX_r[i] = 6.0 * (gradB_x_r[i] * gradC_x_r[i] * rhoD_r[i] 
                    
                                    - gradB_x_i[i] * rhoD_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * rhoD_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDXX_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_x_r[i] 
                    
                                    - gradB_x_i[i] * gradD_x_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_x_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    // CD
                    rBrCrDXX_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_x_r[i] 
                    
                                    - rhoB_i[i] * gradD_x_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * gradD_x_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_x_i[i] * gradD_x_r[i]);


                    // BC

                    rBrCrDXY_r[i] = 6.0 * (gradB_x_r[i] * gradC_y_r[i] * rhoD_r[i] 
                    
                                    - gradB_x_i[i] * rhoD_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * rhoD_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDXY_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_y_r[i] 
                    
                                    - gradB_x_i[i] * gradD_y_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_y_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    // CD
                    rBrCrDXY_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_y_r[i] 
                    
                                    - rhoB_i[i] * gradD_y_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * gradD_y_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_x_i[i] * gradD_y_r[i]);


                    // BC

                    rBrCrDXZ_r[i] = 6.0 * (gradB_x_r[i] * gradC_z_r[i] * rhoD_r[i] 
                    
                                    - gradB_x_i[i] * rhoD_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * rhoD_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDXZ_r[i] += 6.0 * (gradB_x_r[i] * rhoC_r[i] * gradD_z_r[i] 
                    
                                    - gradB_x_i[i] * gradD_z_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_z_i[i] * gradB_x_r[i] 
                                    
                                    - gradB_x_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    // CD

                    rBrCrDXZ_r[i] += 6.0 * (rhoB_r[i] * gradC_x_r[i] * gradD_z_r[i] 
                    
                                    - rhoB_i[i] * gradD_z_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * gradD_z_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    // BC

                    rBrCrDYX_r[i] = 6.0 * (gradB_y_r[i] * gradC_x_r[i] * rhoD_r[i] 
                    
                                    - gradB_y_i[i] * rhoD_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * rhoD_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDYX_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_x_r[i] 
                    
                                    - gradB_y_i[i] * gradD_x_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_x_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    // CD
                    rBrCrDYX_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_x_r[i] 
                    
                                    - rhoB_i[i] * gradD_x_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * gradD_x_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    // BC

                    rBrCrDYY_r[i] = 6.0 * (gradB_y_r[i] * gradC_y_r[i] * rhoD_r[i] 
                    
                                    - gradB_y_i[i] * rhoD_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * rhoD_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDYY_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_y_r[i] 
                    
                                    - gradB_y_i[i] * gradD_y_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_y_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    // CD
                    rBrCrDYY_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_y_r[i] 
                    
                                    - rhoB_i[i] * gradD_y_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * gradD_y_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_y_i[i] * gradD_y_r[i]);


                    // BC

                    rBrCrDYZ_r[i] = 6.0 * (gradB_y_r[i] * gradC_z_r[i] * rhoD_r[i] 
                    
                                    - gradB_y_i[i] * rhoD_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * rhoD_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDYZ_r[i] += 6.0 * (gradB_y_r[i] * rhoC_r[i] * gradD_z_r[i] 
                    
                                    - gradB_y_i[i] * gradD_z_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_z_i[i] * gradB_y_r[i] 
                                    
                                    - gradB_y_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    // CD
                    rBrCrDYZ_r[i] += 6.0 * (rhoB_r[i] * gradC_y_r[i] * gradD_z_r[i] 
                    
                                    - rhoB_i[i] * gradD_z_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * gradD_z_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_y_i[i] * gradD_z_r[i]);


                    // BC

                    rBrCrDZX_r[i] = 6.0 * (gradB_z_r[i] * gradC_x_r[i] * rhoD_r[i] 
                    
                                    - gradB_z_i[i] * rhoD_i[i] * gradC_x_r[i] 
                                    
                                    - gradC_x_i[i] * rhoD_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * gradC_x_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDZX_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_x_r[i] 
                    
                                    - gradB_z_i[i] * gradD_x_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_x_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * rhoC_i[i] * gradD_x_r[i]);

                    // CD
                    rBrCrDZX_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_x_r[i] 
                    
                                    - rhoB_i[i] * gradD_x_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * gradD_x_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    // BC

                    rBrCrDZY_r[i] = 6.0 * (gradB_z_r[i] * gradC_y_r[i] * rhoD_r[i] 
                    
                                    - gradB_z_i[i] * rhoD_i[i] * gradC_y_r[i] 
                                    
                                    - gradC_y_i[i] * rhoD_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * gradC_y_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDZY_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_y_r[i] 
                    
                                    - gradB_z_i[i] * gradD_y_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_y_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * rhoC_i[i] * gradD_y_r[i]);

                    // CD
                    rBrCrDZY_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_y_r[i] 
                    
                                    - rhoB_i[i] * gradD_y_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * gradD_y_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    // BC

                    rBrCrDZZ_r[i] = 6.0 * (gradB_z_r[i] * gradC_z_r[i] * rhoD_r[i] 
                    
                                    - gradB_z_i[i] * rhoD_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * rhoD_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * gradC_z_i[i] * rhoD_r[i]);

                    // BD

                    rBrCrDZZ_r[i] += 6.0 * (gradB_z_r[i] * rhoC_r[i] * gradD_z_r[i] 
                    
                                    - gradB_z_i[i] * gradD_z_i[i] * rhoC_r[i] 
                                    
                                    - rhoC_i[i] * gradD_z_i[i] * gradB_z_r[i] 
                                    
                                    - gradB_z_i[i] * rhoC_i[i] * gradD_z_r[i]);

                    // CD
                    rBrCrDZZ_r[i] += 6.0 * (rhoB_r[i] * gradC_z_r[i] * gradD_z_r[i] 
                    
                                    - rhoB_i[i] * gradD_z_i[i] * gradC_z_r[i] 
                                    
                                    - gradC_z_i[i] * gradD_z_i[i] * rhoB_r[i] 
                                    
                                    - rhoB_i[i] * gradC_z_i[i] * gradD_z_r[i]);


                    // BC
                    rBrCrDXX_i[i] =  6.0 * ( - gradB_x_i[i] * gradC_x_i[i] * rhoD_i[i] 
                                            + gradB_x_i[i] * gradC_x_r[i] * rhoD_r[i]
                                            + gradC_x_i[i] * gradB_x_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_x_r[i] * gradC_x_r[i]);

                    // BD
                    rBrCrDXX_i[i] +=  6.0 *(- gradB_x_i[i] * rhoC_i[i] * gradD_x_i[i]
                                            + gradB_x_i[i] * rhoC_r[i] * gradD_x_r[i]
                                            + rhoC_i[i] * gradB_x_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_x_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDXX_i[i] +=  6.0 * (- rhoB_i[i] * gradC_x_i[i] * gradD_x_i[i] 
                                            + rhoB_i[i] * gradC_x_r[i] * gradD_x_r[i]
                                            + gradC_x_i[i] * rhoB_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * rhoB_r[i] * gradC_x_r[i]);

                    // BC
                    rBrCrDXY_i[i] =  6.0 * ( - gradB_x_i[i] * gradC_y_i[i] * rhoD_i[i] 
                                            + gradB_x_i[i] * gradC_y_r[i] * rhoD_r[i]
                                            + gradC_y_i[i] * gradB_x_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_x_r[i] * gradC_y_r[i]);

                    // BD
                    rBrCrDXY_i[i] +=  6.0 *(- gradB_x_i[i] * rhoC_i[i] * gradD_y_i[i]
                                            + gradB_x_i[i] * rhoC_r[i] * gradD_y_r[i]
                                            + rhoC_i[i] * gradB_x_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_x_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDXY_i[i] +=  6.0 * (- rhoB_i[i] * gradC_x_i[i] * gradD_y_i[i] 
                                            + rhoB_i[i] * gradC_x_r[i] * gradD_y_r[i]
                                            + gradC_x_i[i] * rhoB_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * rhoB_r[i] * gradC_x_r[i]);

                    // BC
                    rBrCrDXZ_i[i] =  6.0 * ( - gradB_x_i[i] * gradC_z_i[i] * rhoD_i[i] 
                                            + gradB_x_i[i] * gradC_z_r[i] * rhoD_r[i]
                                            + gradC_z_i[i] * gradB_x_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_x_r[i] * gradC_z_r[i]);

                    // BD
                    rBrCrDXZ_i[i] +=  6.0 *(- gradB_x_i[i] * rhoC_i[i] * gradD_z_i[i]
                                            + gradB_x_i[i] * rhoC_r[i] * gradD_z_r[i]
                                            + rhoC_i[i] * gradB_x_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_x_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDXZ_i[i] +=  6.0 * (- rhoB_i[i] * gradC_x_i[i] * gradD_z_i[i] 
                                            + rhoB_i[i] * gradC_x_r[i] * gradD_z_r[i]
                                            + gradC_x_i[i] * rhoB_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * rhoB_r[i] * gradC_x_r[i]);

                    // BC
                    rBrCrDYX_i[i] =  6.0 * ( - gradB_y_i[i] * gradC_x_i[i] * rhoD_i[i] 
                                            + gradB_y_i[i] * gradC_x_r[i] * rhoD_r[i]
                                            + gradC_x_i[i] * gradB_y_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_y_r[i] * gradC_x_r[i]);

                    // BD
                    rBrCrDYX_i[i] +=  6.0 *(- gradB_y_i[i] * rhoC_i[i] * gradD_x_i[i]
                                            + gradB_y_i[i] * rhoC_r[i] * gradD_x_r[i]
                                            + rhoC_i[i] * gradB_y_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_y_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDYX_i[i] +=  6.0 * (- rhoB_i[i] * gradC_y_i[i] * gradD_x_i[i] 
                                            + rhoB_i[i] * gradC_y_r[i] * gradD_x_r[i]
                                            + gradC_y_i[i] * rhoB_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * rhoB_r[i] * gradC_y_r[i]);

                    // BC
                    rBrCrDYY_i[i] =  6.0 * ( - gradB_y_i[i] * gradC_y_i[i] * rhoD_i[i] 
                                            + gradB_y_i[i] * gradC_y_r[i] * rhoD_r[i]
                                            + gradC_y_i[i] * gradB_y_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_y_r[i] * gradC_y_r[i]);

                    // BD
                    rBrCrDYY_i[i] +=  6.0 *(- gradB_y_i[i] * rhoC_i[i] * gradD_y_i[i]
                                            + gradB_y_i[i] * rhoC_r[i] * gradD_y_r[i]
                                            + rhoC_i[i] * gradB_y_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_y_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDYY_i[i] +=  6.0 * (- rhoB_i[i] * gradC_y_i[i] * gradD_y_i[i] 
                                            + rhoB_i[i] * gradC_y_r[i] * gradD_y_r[i]
                                            + gradC_y_i[i] * rhoB_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * rhoB_r[i] * gradC_y_r[i]);

                    // BC
                    rBrCrDYZ_i[i] =  6.0 * ( - gradB_y_i[i] * gradC_z_i[i] * rhoD_i[i] 
                                            + gradB_y_i[i] * gradC_z_r[i] * rhoD_r[i]
                                            + gradC_z_i[i] * gradB_y_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_y_r[i] * gradC_z_r[i]);

                    // BD
                    rBrCrDYZ_i[i] +=  6.0 *(- gradB_y_i[i] * rhoC_i[i] * gradD_z_i[i]
                                            + gradB_y_i[i] * rhoC_r[i] * gradD_z_r[i]
                                            + rhoC_i[i] * gradB_y_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_y_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDYZ_i[i] +=  6.0 * (- rhoB_i[i] * gradC_y_i[i] * gradD_z_i[i] 
                                            + rhoB_i[i] * gradC_y_r[i] * gradD_z_r[i]
                                            + gradC_y_i[i] * rhoB_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * rhoB_r[i] * gradC_y_r[i]);

                    // BC
                    rBrCrDZX_i[i] =  6.0 * ( - gradB_z_i[i] * gradC_x_i[i] * rhoD_i[i] 
                                            + gradB_z_i[i] * gradC_x_r[i] * rhoD_r[i]
                                            + gradC_x_i[i] * gradB_z_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_z_r[i] * gradC_x_r[i]);

                    // BD
                    rBrCrDZX_i[i] +=  6.0 *(- gradB_z_i[i] * rhoC_i[i] * gradD_x_i[i]
                                            + gradB_z_i[i] * rhoC_r[i] * gradD_x_r[i]
                                            + rhoC_i[i] * gradB_z_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_z_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDZX_i[i] +=  6.0 * (- rhoB_i[i] * gradC_z_i[i] * gradD_x_i[i] 
                                            + rhoB_i[i] * gradC_z_r[i] * gradD_x_r[i]
                                            + gradC_z_i[i] * rhoB_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * rhoB_r[i] * gradC_z_r[i]);

                    // BC
                    rBrCrDZY_i[i] =  6.0 * ( - gradB_z_i[i] * gradC_y_i[i] * rhoD_i[i] 
                                            + gradB_z_i[i] * gradC_y_r[i] * rhoD_r[i]
                                            + gradC_y_i[i] * gradB_z_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_z_r[i] * gradC_y_r[i]);

                    // BD
                    rBrCrDZY_i[i] +=  6.0 *(- gradB_z_i[i] * rhoC_i[i] * gradD_y_i[i]
                                            + gradB_z_i[i] * rhoC_r[i] * gradD_y_r[i]
                                            + rhoC_i[i] * gradB_z_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_z_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDZY_i[i] +=  6.0 * (- rhoB_i[i] * gradC_z_i[i] * gradD_y_i[i] 
                                            + rhoB_i[i] * gradC_z_r[i] * gradD_y_r[i]
                                            + gradC_z_i[i] * rhoB_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * rhoB_r[i] * gradC_z_r[i]);

                    // BC
                    rBrCrDZZ_i[i] =  6.0 * ( - gradB_z_i[i] * gradC_z_i[i] * rhoD_i[i] 
                                            + gradB_z_i[i] * gradC_z_r[i] * rhoD_r[i]
                                            + gradC_z_i[i] * gradB_z_r[i] * rhoD_r[i] 
                                            + rhoD_i[i] * gradB_z_r[i] * gradC_z_r[i]);

                    // BD
                    rBrCrDZZ_i[i] +=  6.0 *(- gradB_z_i[i] * rhoC_i[i] * gradD_z_i[i]
                                            + gradB_z_i[i] * rhoC_r[i] * gradD_z_r[i]
                                            + rhoC_i[i] * gradB_z_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_z_r[i] * rhoC_r[i]);

                    // CD
                    rBrCrDZZ_i[i] +=  6.0 * (- rhoB_i[i] * gradC_z_i[i] * gradD_z_i[i] 
                                            + rhoB_i[i] * gradC_z_r[i] * gradD_z_r[i]
                                            + gradC_z_i[i] * rhoB_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * rhoB_r[i] * gradC_z_r[i]);


                    rBrCrD_XXX_r[i] =  6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                            - gradB_x_i[i] * gradD_x_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_x_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    rBrCrD_XXY_r[i] =  6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                            - gradB_x_i[i] * gradD_y_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_y_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    rBrCrD_XXZ_r[i] =  6.0 * (gradB_x_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                            - gradB_x_i[i] * gradD_z_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_z_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    rBrCrD_XYX_r[i] =  6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                            - gradB_x_i[i] * gradD_x_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_x_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    rBrCrD_XYY_r[i] =  6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                            - gradB_x_i[i] * gradD_y_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_y_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    rBrCrD_XYZ_r[i] =  6.0 * (gradB_x_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                            - gradB_x_i[i] * gradD_z_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_z_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    rBrCrD_XZX_r[i] =  6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                            - gradB_x_i[i] * gradD_x_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_x_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    rBrCrD_XZY_r[i] =  6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                            - gradB_x_i[i] * gradD_y_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_y_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    rBrCrD_XZZ_r[i] =  6.0 * (gradB_x_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                            - gradB_x_i[i] * gradD_z_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_z_i[i] * gradB_x_r[i] 
                                            - gradB_x_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    rBrCrD_YXX_r[i] =  6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                            - gradB_y_i[i] * gradD_x_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_x_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    rBrCrD_YXY_r[i] =  6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                            - gradB_y_i[i] * gradD_y_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_y_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    rBrCrD_YXZ_r[i] =  6.0 * (gradB_y_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                            - gradB_y_i[i] * gradD_z_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_z_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    rBrCrD_YYX_r[i] =  6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                            - gradB_y_i[i] * gradD_x_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_x_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    rBrCrD_YYY_r[i] =  6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                            - gradB_y_i[i] * gradD_y_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_y_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    rBrCrD_YYZ_r[i] =  6.0 * (gradB_y_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                            - gradB_y_i[i] * gradD_z_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_z_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    rBrCrD_YZX_r[i] =  6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                            - gradB_y_i[i] * gradD_x_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_x_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    rBrCrD_YZY_r[i] =  6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                            - gradB_y_i[i] * gradD_y_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_y_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    rBrCrD_YZZ_r[i] =  6.0 * (gradB_y_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                            - gradB_y_i[i] * gradD_z_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_z_i[i] * gradB_y_r[i] 
                                            - gradB_y_i[i] * gradC_z_i[i] * gradD_z_r[i]);

                    rBrCrD_ZXX_r[i] =  6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_x_r[i]
                                            - gradB_z_i[i] * gradD_x_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_x_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_x_i[i] * gradD_x_r[i]);

                    rBrCrD_ZXY_r[i] =  6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_y_r[i]
                                            - gradB_z_i[i] * gradD_y_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_y_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_x_i[i] * gradD_y_r[i]);

                    rBrCrD_ZXZ_r[i] =  6.0 * (gradB_z_r[i] * gradC_x_r[i] * gradD_z_r[i]
                                            - gradB_z_i[i] * gradD_z_i[i] * gradC_x_r[i]
                                            - gradC_x_i[i] * gradD_z_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_x_i[i] * gradD_z_r[i]);

                    rBrCrD_ZYX_r[i] =  6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_x_r[i]
                                            - gradB_z_i[i] * gradD_x_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_x_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_y_i[i] * gradD_x_r[i]);

                    rBrCrD_ZYY_r[i] =  6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_y_r[i]
                                            - gradB_z_i[i] * gradD_y_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_y_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_y_i[i] * gradD_y_r[i]);

                    rBrCrD_ZYZ_r[i] =  6.0 * (gradB_z_r[i] * gradC_y_r[i] * gradD_z_r[i]
                                            - gradB_z_i[i] * gradD_z_i[i] * gradC_y_r[i]
                                            - gradC_y_i[i] * gradD_z_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_y_i[i] * gradD_z_r[i]);

                    rBrCrD_ZZX_r[i] =  6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_x_r[i]
                                            - gradB_z_i[i] * gradD_x_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_x_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_z_i[i] * gradD_x_r[i]);

                    rBrCrD_ZZY_r[i] =  6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_y_r[i]
                                            - gradB_z_i[i] * gradD_y_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_y_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_z_i[i] * gradD_y_r[i]);

                    rBrCrD_ZZZ_r[i] =  6.0 * (gradB_z_r[i] * gradC_z_r[i] * gradD_z_r[i]
                                            - gradB_z_i[i] * gradD_z_i[i] * gradC_z_r[i]
                                            - gradC_z_i[i] * gradD_z_i[i] * gradB_z_r[i] 
                                            - gradB_z_i[i] * gradC_z_i[i] * gradD_z_r[i]);


                    rBrCrD_XXX_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                            + gradB_x_i[i] * gradC_x_r[i] * gradD_x_r[i] 
                                            + gradC_x_i[i] * gradB_x_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_x_r[i] * gradC_x_r[i]);

                    rBrCrD_XXY_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                            + gradB_x_i[i] * gradC_x_r[i] * gradD_y_r[i] 
                                            + gradC_x_i[i] * gradB_x_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_x_r[i] * gradC_x_r[i]);

                    rBrCrD_XXZ_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                            + gradB_x_i[i] * gradC_x_r[i] * gradD_z_r[i] 
                                            + gradC_x_i[i] * gradB_x_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_x_r[i] * gradC_x_r[i]);

                    rBrCrD_XYX_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                            + gradB_x_i[i] * gradC_y_r[i] * gradD_x_r[i] 
                                            + gradC_y_i[i] * gradB_x_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_x_r[i] * gradC_y_r[i]);

                    rBrCrD_XYY_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                            + gradB_x_i[i] * gradC_y_r[i] * gradD_y_r[i] 
                                            + gradC_y_i[i] * gradB_x_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_x_r[i] * gradC_y_r[i]);

                    rBrCrD_XYZ_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                            + gradB_x_i[i] * gradC_y_r[i] * gradD_z_r[i] 
                                            + gradC_y_i[i] * gradB_x_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_x_r[i] * gradC_y_r[i]);

                    rBrCrD_XZX_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                            + gradB_x_i[i] * gradC_z_r[i] * gradD_x_r[i] 
                                            + gradC_z_i[i] * gradB_x_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_x_r[i] * gradC_z_r[i]);

                    rBrCrD_XZY_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                            + gradB_x_i[i] * gradC_z_r[i] * gradD_y_r[i] 
                                            + gradC_z_i[i] * gradB_x_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_x_r[i] * gradC_z_r[i]);

                    rBrCrD_XZZ_i[i] = 6.0 *  (- gradB_x_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                            + gradB_x_i[i] * gradC_z_r[i] * gradD_z_r[i] 
                                            + gradC_z_i[i] * gradB_x_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_x_r[i] * gradC_z_r[i]);

                    rBrCrD_YXX_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                            + gradB_y_i[i] * gradC_x_r[i] * gradD_x_r[i] 
                                            + gradC_x_i[i] * gradB_y_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_y_r[i] * gradC_x_r[i]);

                    rBrCrD_YXY_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                            + gradB_y_i[i] * gradC_x_r[i] * gradD_y_r[i] 
                                            + gradC_x_i[i] * gradB_y_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_y_r[i] * gradC_x_r[i]);

                    rBrCrD_YXZ_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                            + gradB_y_i[i] * gradC_x_r[i] * gradD_z_r[i] 
                                            + gradC_x_i[i] * gradB_y_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_y_r[i] * gradC_x_r[i]);

                    rBrCrD_YYX_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                            + gradB_y_i[i] * gradC_y_r[i] * gradD_x_r[i] 
                                            + gradC_y_i[i] * gradB_y_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_y_r[i] * gradC_y_r[i]);

                    rBrCrD_YYY_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                            + gradB_y_i[i] * gradC_y_r[i] * gradD_y_r[i] 
                                            + gradC_y_i[i] * gradB_y_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_y_r[i] * gradC_y_r[i]);

                    rBrCrD_YYZ_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                            + gradB_y_i[i] * gradC_y_r[i] * gradD_z_r[i] 
                                            + gradC_y_i[i] * gradB_y_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_y_r[i] * gradC_y_r[i]);

                    rBrCrD_YZX_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                            + gradB_y_i[i] * gradC_z_r[i] * gradD_x_r[i] 
                                            + gradC_z_i[i] * gradB_y_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_y_r[i] * gradC_z_r[i]);

                    rBrCrD_YZY_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                            + gradB_y_i[i] * gradC_z_r[i] * gradD_y_r[i] 
                                            + gradC_z_i[i] * gradB_y_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_y_r[i] * gradC_z_r[i]);

                    rBrCrD_YZZ_i[i] = 6.0 *  (- gradB_y_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                            + gradB_y_i[i] * gradC_z_r[i] * gradD_z_r[i] 
                                            + gradC_z_i[i] * gradB_y_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_y_r[i] * gradC_z_r[i]);

                    rBrCrD_ZXX_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_x_i[i] * gradD_x_i[i]
                                            + gradB_z_i[i] * gradC_x_r[i] * gradD_x_r[i] 
                                            + gradC_x_i[i] * gradB_z_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_z_r[i] * gradC_x_r[i]);

                    rBrCrD_ZXY_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_x_i[i] * gradD_y_i[i]
                                            + gradB_z_i[i] * gradC_x_r[i] * gradD_y_r[i] 
                                            + gradC_x_i[i] * gradB_z_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_z_r[i] * gradC_x_r[i]);

                    rBrCrD_ZXZ_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_x_i[i] * gradD_z_i[i]
                                            + gradB_z_i[i] * gradC_x_r[i] * gradD_z_r[i] 
                                            + gradC_x_i[i] * gradB_z_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_z_r[i] * gradC_x_r[i]);

                    rBrCrD_ZYX_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_y_i[i] * gradD_x_i[i]
                                            + gradB_z_i[i] * gradC_y_r[i] * gradD_x_r[i] 
                                            + gradC_y_i[i] * gradB_z_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_z_r[i] * gradC_y_r[i]);

                    rBrCrD_ZYY_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_y_i[i] * gradD_y_i[i]
                                            + gradB_z_i[i] * gradC_y_r[i] * gradD_y_r[i] 
                                            + gradC_y_i[i] * gradB_z_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_z_r[i] * gradC_y_r[i]);

                    rBrCrD_ZYZ_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_y_i[i] * gradD_z_i[i]
                                            + gradB_z_i[i] * gradC_y_r[i] * gradD_z_r[i] 
                                            + gradC_y_i[i] * gradB_z_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_z_r[i] * gradC_y_r[i]);

                    rBrCrD_ZZX_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_z_i[i] * gradD_x_i[i]
                                            + gradB_z_i[i] * gradC_z_r[i] * gradD_x_r[i] 
                                            + gradC_z_i[i] * gradB_z_r[i] * gradD_x_r[i] 
                                            + gradD_x_i[i] * gradB_z_r[i] * gradC_z_r[i]);

                    rBrCrD_ZZY_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_z_i[i] * gradD_y_i[i]
                                            + gradB_z_i[i] * gradC_z_r[i] * gradD_y_r[i] 
                                            + gradC_z_i[i] * gradB_z_r[i] * gradD_y_r[i] 
                                            + gradD_y_i[i] * gradB_z_r[i] * gradC_z_r[i]);

                    rBrCrD_ZZZ_i[i] = 6.0 *  (- gradB_z_i[i] * gradC_z_i[i] * gradD_z_i[i]
                                            + gradB_z_i[i] * gradC_z_r[i] * gradD_z_r[i] 
                                            + gradC_z_i[i] * gradB_z_r[i] * gradD_z_r[i] 
                                            + gradD_z_i[i] * gradB_z_r[i] * gradC_z_r[i]);


                    RhoBCRhoD_r[i] = 3.0 * (rhoD_r[i] * rhoBC_r[i] - rhoD_i[i] * rhoBC_i[i] 
                    
                                         + rhoC_r[i] * rhoBD_r[i] - rhoC_i[i] * rhoBD_i[i] 
                                         
                                         + rhoB_r[i] * rhoCD_r[i] - rhoB_i[i] * rhoCD_i[i]) ;

                    RhoBCRhoD_i[i] = 3.0 * (rhoD_i[i] * rhoBC_r[i] + rhoD_r[i] * rhoBC_i[i]
                    
                                         + rhoC_i[i] * rhoBD_r[i] + rhoC_r[i] * rhoBD_i[i] 
                                         
                                        + rhoB_i[i] * rhoCD_r[i] + rhoB_r[i] * rhoCD_i[i]);


                    RhoBCRhoD_X_r[i] = 3.0 * ( gradD_x_r[i] * rhoBC_r[i] - gradD_x_i[i] * rhoBC_i[i]
                                            + gradC_x_r[i] * rhoBD_r[i] - gradC_x_i[i] * rhoBD_i[i]
                                            + gradB_x_r[i] * rhoCD_r[i] - gradB_x_i[i] * rhoCD_i[i]) ;

                    RhoBCRhoD_Y_r[i] = 3.0 * ( gradD_y_r[i] * rhoBC_r[i] - gradD_y_i[i] * rhoBC_i[i]
                                            + gradC_y_r[i] * rhoBD_r[i] - gradC_y_i[i] * rhoBD_i[i]
                                            + gradB_y_r[i] * rhoCD_r[i] - gradB_y_i[i] * rhoCD_i[i]) ;

                    RhoBCRhoD_Z_r[i] = 3.0 * ( gradD_z_r[i] * rhoBC_r[i] - gradD_z_i[i] * rhoBC_i[i]
                                            + gradC_z_r[i] * rhoBD_r[i] - gradC_z_i[i] * rhoBD_i[i]
                                            + gradB_z_r[i] * rhoCD_r[i] - gradB_z_i[i] * rhoCD_i[i]) ;

                    RhoBCRhoD_X_r[i] += 3.0 * ( rhoD_r[i] * gradBC_x_r[i] - rhoD_i[i] * gradBC_x_i[i]
                                            + rhoC_r[i] * gradBD_x_r[i] - rhoC_i[i] * gradBD_x_i[i]
                                            + rhoB_r[i] * gradCD_x_r[i] - rhoB_i[i] * gradCD_x_i[i]) ;

                    RhoBCRhoD_Y_r[i] += 3.0 * ( rhoD_r[i] * gradBC_y_r[i] - rhoD_i[i] * gradBC_y_i[i]
                                            + rhoC_r[i] * gradBD_y_r[i] - rhoC_i[i] * gradBD_y_i[i]
                                            + rhoB_r[i] * gradCD_y_r[i] - rhoB_i[i] * gradCD_y_i[i]) ;

                    RhoBCRhoD_Z_r[i] += 3.0 * ( rhoD_r[i] * gradBC_z_r[i] - rhoD_i[i] * gradBC_z_i[i]
                                            + rhoC_r[i] * gradBD_z_r[i] - rhoC_i[i] * gradBD_z_i[i]
                                            + rhoB_r[i] * gradCD_z_r[i] - rhoB_i[i] * gradCD_z_i[i]) ;

                    RhoBCRhoD_XX_r[i] = 3.0 * ( gradD_x_r[i] * gradBC_x_r[i] - gradD_x_i[i] * gradBC_x_i[i]
                                            + gradC_x_r[i] * gradBD_x_r[i] - gradC_x_i[i] * gradBD_x_i[i]
                                            + gradB_x_r[i] * gradCD_x_r[i] - gradB_x_i[i] * gradCD_x_i[i]) ;

                    RhoBCRhoD_XY_r[i] = 3.0 * ( gradD_x_r[i] * gradBC_y_r[i] - gradD_x_i[i] * gradBC_y_i[i]
                                            + gradC_x_r[i] * gradBD_y_r[i] - gradC_x_i[i] * gradBD_y_i[i]
                                            + gradB_x_r[i] * gradCD_y_r[i] - gradB_x_i[i] * gradCD_y_i[i]) ;

                    RhoBCRhoD_XZ_r[i] = 3.0 * ( gradD_x_r[i] * gradBC_z_r[i] - gradD_x_i[i] * gradBC_z_i[i]
                                            + gradC_x_r[i] * gradBD_z_r[i] - gradC_x_i[i] * gradBD_z_i[i]
                                            + gradB_x_r[i] * gradCD_z_r[i] - gradB_x_i[i] * gradCD_z_i[i]) ;

                    RhoBCRhoD_YX_r[i] = 3.0 * ( gradD_y_r[i] * gradBC_x_r[i] - gradD_y_i[i] * gradBC_x_i[i]
                                            + gradC_y_r[i] * gradBD_x_r[i] - gradC_y_i[i] * gradBD_x_i[i]
                                            + gradB_y_r[i] * gradCD_x_r[i] - gradB_y_i[i] * gradCD_x_i[i]) ;

                    RhoBCRhoD_YY_r[i] = 3.0 * ( gradD_y_r[i] * gradBC_y_r[i] - gradD_y_i[i] * gradBC_y_i[i]
                                            + gradC_y_r[i] * gradBD_y_r[i] - gradC_y_i[i] * gradBD_y_i[i]
                                            + gradB_y_r[i] * gradCD_y_r[i] - gradB_y_i[i] * gradCD_y_i[i]) ;

                    RhoBCRhoD_YZ_r[i] = 3.0 * ( gradD_y_r[i] * gradBC_z_r[i] - gradD_y_i[i] * gradBC_z_i[i]
                                            + gradC_y_r[i] * gradBD_z_r[i] - gradC_y_i[i] * gradBD_z_i[i]
                                            + gradB_y_r[i] * gradCD_z_r[i] - gradB_y_i[i] * gradCD_z_i[i]) ;

                    RhoBCRhoD_ZX_r[i] = 3.0 * ( gradD_z_r[i] * gradBC_x_r[i] - gradD_z_i[i] * gradBC_x_i[i]
                                            + gradC_z_r[i] * gradBD_x_r[i] - gradC_z_i[i] * gradBD_x_i[i]
                                            + gradB_z_r[i] * gradCD_x_r[i] - gradB_z_i[i] * gradCD_x_i[i]) ;

                    RhoBCRhoD_ZY_r[i] = 3.0 * ( gradD_z_r[i] * gradBC_y_r[i] - gradD_z_i[i] * gradBC_y_i[i]
                                            + gradC_z_r[i] * gradBD_y_r[i] - gradC_z_i[i] * gradBD_y_i[i]
                                            + gradB_z_r[i] * gradCD_y_r[i] - gradB_z_i[i] * gradCD_y_i[i]) ;

                    RhoBCRhoD_ZZ_r[i] = 3.0 * ( gradD_z_r[i] * gradBC_z_r[i] - gradD_z_i[i] * gradBC_z_i[i]
                                            + gradC_z_r[i] * gradBD_z_r[i] - gradC_z_i[i] * gradBD_z_i[i]
                                            + gradB_z_r[i] * gradCD_z_r[i] - gradB_z_i[i] * gradCD_z_i[i]) ;

                    RhoBCRhoD_X_i[i]= 3.0 * (gradD_x_i[i] * rhoBC_r[i] + gradD_x_r[i] * rhoBC_i[i]
                                            + gradC_x_i[i] * rhoBD_r[i] + gradC_x_r[i] * rhoBD_i[i]
                                            + gradB_x_i[i] * rhoCD_r[i] + gradB_x_r[i] * rhoCD_i[i]);

                    RhoBCRhoD_Y_i[i]= 3.0 * (gradD_y_i[i] * rhoBC_r[i] + gradD_y_r[i] * rhoBC_i[i]
                                            + gradC_y_i[i] * rhoBD_r[i] + gradC_y_r[i] * rhoBD_i[i]
                                            + gradB_y_i[i] * rhoCD_r[i] + gradB_y_r[i] * rhoCD_i[i]);

                    RhoBCRhoD_Z_i[i]= 3.0 * (gradD_z_i[i] * rhoBC_r[i] + gradD_z_r[i] * rhoBC_i[i]
                                            + gradC_z_i[i] * rhoBD_r[i] + gradC_z_r[i] * rhoBD_i[i]
                                            + gradB_z_i[i] * rhoCD_r[i] + gradB_z_r[i] * rhoCD_i[i]);

                    RhoBCRhoD_X_i[i]+= 3.0 * (rhoD_i[i] * gradBC_x_r[i] + rhoD_r[i] * gradBC_x_i[i]
                                            + rhoC_i[i] * gradBD_x_r[i] + rhoC_r[i] * gradBD_x_i[i]
                                            + rhoB_i[i] * gradCD_x_r[i] + rhoB_r[i] * gradCD_x_i[i]);

                    RhoBCRhoD_Y_i[i]+= 3.0 * (rhoD_i[i] * gradBC_y_r[i] + rhoD_r[i] * gradBC_y_i[i]
                                            + rhoC_i[i] * gradBD_y_r[i] + rhoC_r[i] * gradBD_y_i[i]
                                            + rhoB_i[i] * gradCD_y_r[i] + rhoB_r[i] * gradCD_y_i[i]);

                    RhoBCRhoD_Z_i[i]+= 3.0 * (rhoD_i[i] * gradBC_z_r[i] + rhoD_r[i] * gradBC_z_i[i]
                                            + rhoC_i[i] * gradBD_z_r[i] + rhoC_r[i] * gradBD_z_i[i]
                                            + rhoB_i[i] * gradCD_z_r[i] + rhoB_r[i] * gradCD_z_i[i]);

                    RhoBCRhoD_XX_i[i]= 3.0 * (gradD_x_i[i] * gradBC_x_r[i] + gradD_x_r[i] * gradBC_x_i[i]
                                            + gradC_x_i[i] * gradBD_x_r[i] + gradC_x_r[i] * gradBD_x_i[i]
                                            + gradB_x_i[i] * gradCD_x_r[i] + gradB_x_r[i] * gradCD_x_i[i]);

                    RhoBCRhoD_XY_i[i]= 3.0 * (gradD_x_i[i] * gradBC_y_r[i] + gradD_x_r[i] * gradBC_y_i[i]
                                            + gradC_x_i[i] * gradBD_y_r[i] + gradC_x_r[i] * gradBD_y_i[i]
                                            + gradB_x_i[i] * gradCD_y_r[i] + gradB_x_r[i] * gradCD_y_i[i]);

                    RhoBCRhoD_XZ_i[i]= 3.0 * (gradD_x_i[i] * gradBC_z_r[i] + gradD_x_r[i] * gradBC_z_i[i]
                                            + gradC_x_i[i] * gradBD_z_r[i] + gradC_x_r[i] * gradBD_z_i[i]
                                            + gradB_x_i[i] * gradCD_z_r[i] + gradB_x_r[i] * gradCD_z_i[i]);

                    RhoBCRhoD_YX_i[i]= 3.0 * (gradD_y_i[i] * gradBC_x_r[i] + gradD_y_r[i] * gradBC_x_i[i]
                                            + gradC_y_i[i] * gradBD_x_r[i] + gradC_y_r[i] * gradBD_x_i[i]
                                            + gradB_y_i[i] * gradCD_x_r[i] + gradB_y_r[i] * gradCD_x_i[i]);

                    RhoBCRhoD_YY_i[i]= 3.0 * (gradD_y_i[i] * gradBC_y_r[i] + gradD_y_r[i] * gradBC_y_i[i]
                                            + gradC_y_i[i] * gradBD_y_r[i] + gradC_y_r[i] * gradBD_y_i[i]
                                            + gradB_y_i[i] * gradCD_y_r[i] + gradB_y_r[i] * gradCD_y_i[i]);

                    RhoBCRhoD_YZ_i[i]= 3.0 * (gradD_y_i[i] * gradBC_z_r[i] + gradD_y_r[i] * gradBC_z_i[i]
                                            + gradC_y_i[i] * gradBD_z_r[i] + gradC_y_r[i] * gradBD_z_i[i]
                                            + gradB_y_i[i] * gradCD_z_r[i] + gradB_y_r[i] * gradCD_z_i[i]);

                    RhoBCRhoD_ZX_i[i]= 3.0 * (gradD_z_i[i] * gradBC_x_r[i] + gradD_z_r[i] * gradBC_x_i[i]
                                            + gradC_z_i[i] * gradBD_x_r[i] + gradC_z_r[i] * gradBD_x_i[i]
                                            + gradB_z_i[i] * gradCD_x_r[i] + gradB_z_r[i] * gradCD_x_i[i]);

                    RhoBCRhoD_ZY_i[i]= 3.0 * (gradD_z_i[i] * gradBC_y_r[i] + gradD_z_r[i] * gradBC_y_i[i]
                                            + gradC_z_i[i] * gradBD_y_r[i] + gradC_z_r[i] * gradBD_y_i[i]
                                            + gradB_z_i[i] * gradCD_y_r[i] + gradB_z_r[i] * gradCD_y_i[i]);

                    RhoBCRhoD_ZZ_i[i]= 3.0 * (gradD_z_i[i] * gradBC_z_r[i] + gradD_z_r[i] * gradBC_z_i[i]
                                            + gradC_z_i[i] * gradBD_z_r[i] + gradC_z_r[i] * gradBD_z_i[i]
                                            + gradB_z_i[i] * gradCD_z_r[i] + gradB_z_r[i] * gradCD_z_i[i]);


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
