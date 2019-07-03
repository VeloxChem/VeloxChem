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
{
}

CDensityGrid::CDensityGrid(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _densityValues(densityValues)
{
}

CDensityGrid::CDensityGrid(const CDensityGrid& source)

    : _gridType(source._gridType)

    , _densityValues(source._densityValues)
{
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;
    
    int32_t ncomp = 0;
    
    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 2 : 1;
    
    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 5 : 2;
    
    // NOTE: this needs to be checked with mgga functionals implementation
    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 7 : 3;

    _densityValues = CMemBlock2D<double>(nGridPoints, ncomp); 
}

CDensityGrid::CDensityGrid(CDensityGrid&& source) noexcept

    : _gridType(std::move(source._gridType))

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
    
    _densityValues = source._densityValues;
    
    return *this;
}

CDensityGrid&
CDensityGrid::operator=(CDensityGrid&& source) noexcept
{
    if (this == &source) return *this;
    
    _gridType = std::move(source._gridType);
    
    _densityValues = std::move(source._densityValues);
    
    return *this;
}

bool
CDensityGrid::operator==(const CDensityGrid& other) const
{
    if (_gridType != other._gridType) return false;
    
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

const double*
CDensityGrid::alphaDensity() const
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(0);
}

double*
CDensityGrid::alphaDensity()
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(0);
}

const double*
CDensityGrid::betaDensity() const
{
    if (_gridType == dengrid::ab) return _densityValues.data(1);
    
    if (_gridType == dengrid::lima) return _densityValues.data(0);
    
    return nullptr;
}

double*
CDensityGrid::betaDensity()
{
    if (_gridType == dengrid::ab) return _densityValues.data(1);
    
    if (_gridType == dengrid::lima) return _densityValues.data(0);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradient() const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2);
    
    if (_gridType == dengrid::limb) return _densityValues.data(1);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradient()
{
    if (_gridType == dengrid::ab) return _densityValues.data(2);
    
    if (_gridType == dengrid::limb) return _densityValues.data(1);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradient() const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3);
    
    if (_gridType == dengrid::lima) return _densityValues.data(1);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradient()
{
    if (_gridType == dengrid::ab) return _densityValues.data(3);
    
    if (_gridType == dengrid::lima) return _densityValues.data(1);
    
    return nullptr;
}

const double*
CDensityGrid::mixedDensityGradient() const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4);
    
    return nullptr;
}

double*
CDensityGrid::mixedDensityGradient()
{
    if (_gridType == dengrid::ab) return _densityValues.data(4);
    
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
