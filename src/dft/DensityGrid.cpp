//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGrid.hpp"

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
    
    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 2 : 1;
    
    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 5 : 2;
    
    // NOTE: this needs to be checked with mgga functionals implementation
    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 7 : 3;
    
    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
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

const double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(1);
    
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
