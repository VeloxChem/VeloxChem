//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "XCGradientGrid.hpp"

CXCGradientGrid::CXCGradientGrid()

    : _densityGridType(dengrid::undefined)

    , _xcGridType(xcfun::undefined)

    , _xcValues(CMemBlock2D<double>())
{
    
}

CXCGradientGrid::CXCGradientGrid(const CMemBlock2D<double>& xcValues,
                                 const dengrid              densityGridType,
                                 const xcfun                xcGridType)

    : _densityGridType(densityGridType)

    , _xcGridType(xcGridType)

    , _xcValues(xcValues)
{
    
}


CXCGradientGrid::CXCGradientGrid(const int32_t nGridPoints,
                                 const dengrid densityGridType,
                                 const xcfun   xcGridType)
{
    _densityGridType = densityGridType;
    
    _xcGridType = xcGridType;
    
    int32_t ncomp = 0;
    
    if (_xcGridType == xcfun::lda) ncomp = 3;
    
    if (_xcGridType == xcfun::gga) ncomp = 6;
    
    // NOTE: this needs to be checked with mgga functionals implementation
    if (_xcGridType == xcfun::mgga) ncomp = (_densityGridType == dengrid::ab) ? 8 : 4;
    
    _xcValues = CMemBlock2D<double>(nGridPoints, ncomp);
}

CXCGradientGrid::CXCGradientGrid(const CXCGradientGrid& source)

    : _densityGridType(source._densityGridType)

    , _xcGridType(source._xcGridType)

    , _xcValues(source._xcValues)
{
}

CXCGradientGrid::CXCGradientGrid(CXCGradientGrid&& source) noexcept

    : _densityGridType(std::move(source._densityGridType))

    , _xcGridType(std::move(source._xcGridType))

    , _xcValues(std::move(source._xcValues))
{
}

CXCGradientGrid::~CXCGradientGrid()
{
}

CXCGradientGrid&
CXCGradientGrid::operator=(const CXCGradientGrid& source)
{
    if (this == &source) return *this;
    
    _densityGridType = source._densityGridType;
    
    _xcGridType = source._xcGridType;
    
    _xcValues = source._xcValues;
    
    return *this;
}

CXCGradientGrid&
CXCGradientGrid::operator=(CXCGradientGrid&& source) noexcept
{
    if (this == &source) return *this;
    
    _densityGridType = std::move(source._densityGridType);
    
    _xcGridType = std::move(source._xcGridType);
    
    _xcValues = std::move(source._xcValues);
    
    return *this;
}

bool
CXCGradientGrid::operator==(const CXCGradientGrid& other) const
{
    if (_densityGridType != other._densityGridType) return false;
    
    if (_xcGridType != other._xcGridType) return false;
    
    if (_xcValues != other._xcValues) return false;
    
    return true;
}

bool
CXCGradientGrid::operator!=(const CXCGradientGrid& other) const
{
    return !(*this == other);
}

void
CXCGradientGrid::zero()
{
    _xcValues.zero();
}

int32_t
CXCGradientGrid::getNumberOfGridPoints() const
{
    return _xcValues.size(0);
}

const double*
CXCGradientGrid::xcFunctionalValues() const
{
    return _xcValues.data(0);
}

double*
CXCGradientGrid::xcFunctionalValues()
{
    return _xcValues.data(0);
}

const double*
CXCGradientGrid::xcGradientValues(const xcvars gradientType) const
{
    if (gradientType == xcvars::rhoa)
    {
        return _xcValues.data(1);
    }
    
    if (gradientType == xcvars::rhob)
    {
        return _xcValues.data(2);
    }
    
    if (gradientType == xcvars::grada)
    {
        return _xcValues.data(3);
    }
    
    if (gradientType == xcvars::gradb)
    {
        return _xcValues.data(4);

    }
    
    if (gradientType == xcvars::gradab)
    {
        return _xcValues.data(5);
    }
    
    return nullptr;
}

double*
CXCGradientGrid::xcGradientValues(const xcvars gradientType)
{
    if (gradientType == xcvars::rhoa)
    {
        return _xcValues.data(1);
    }
    
    if (gradientType == xcvars::rhob)
    {
         return _xcValues.data(2);
    }
    
    if (gradientType == xcvars::grada)
    {
         return _xcValues.data(3);
    }
    
    if (gradientType == xcvars::gradb)
    {
       return _xcValues.data(4);
    }
    
    if (gradientType == xcvars::gradab)
    {
        return _xcValues.data(5);

    }
    
    return nullptr;
}

dengrid
CXCGradientGrid::getDensityGridType() const
{
    return _densityGridType;
}

std::ostream&
operator<<(std::ostream& output, const CXCGradientGrid& source)
{
    output << std::endl;
    
    output << "[CXCGradientGrid (Object):" << &source << "]" << std::endl;
    
    output << "_densityGridType: " << to_string(source._densityGridType) << std::endl;
    
    output << "_xcGridType: " << to_string(source._xcGridType) << std::endl;
    
    output << "_xcValues: " << std::endl;
    
    output << source._xcValues << std::endl;
    
    return output;
}
