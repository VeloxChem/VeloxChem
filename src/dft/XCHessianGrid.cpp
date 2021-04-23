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

#include "XCHessianGrid.hpp"

CXCHessianGrid::CXCHessianGrid()

    : _densityGridType(dengrid::undefined)

    , _xcGridType(xcfun::undefined)

    , _xcValues(CMemBlock2D<double>())
{
    
}

CXCHessianGrid::CXCHessianGrid(const CMemBlock2D<double>& xcValues,
                                 const dengrid              densityGridType,
                                 const xcfun                xcGridType)

    : _densityGridType(densityGridType)

    , _xcGridType(xcGridType)

    , _xcValues(xcValues)
{
    
}

CXCHessianGrid::CXCHessianGrid(const int32_t nGridPoints,
                               const dengrid densityGridType,
                               const xcfun   xcGridType)
{
    _densityGridType = densityGridType;
    
    _xcGridType = xcGridType;
    
    int32_t ncomp = 0;
    
    if (_xcGridType == xcfun::lda) ncomp = (_densityGridType == dengrid::ab) ? 3 : 1;
    
    if (_xcGridType == xcfun::gga) ncomp = (_densityGridType == dengrid::ab) ? 15 : 3;
    
    // NOTE: this needs to be checked with mgga functionals implementation
    if (_xcGridType == xcfun::mgga) ncomp = (_densityGridType == dengrid::ab) ? 28 : 6;
    
    _xcValues = CMemBlock2D<double>(nGridPoints, ncomp);
}

CXCHessianGrid::CXCHessianGrid(const CXCHessianGrid& source)

    : _densityGridType(source._densityGridType)

    , _xcGridType(source._xcGridType)

    , _xcValues(source._xcValues)
{
}

CXCHessianGrid::CXCHessianGrid(CXCHessianGrid&& source) noexcept

    : _densityGridType(std::move(source._densityGridType))

    , _xcGridType(std::move(source._xcGridType))

    , _xcValues(std::move(source._xcValues))
{
}

CXCHessianGrid::~CXCHessianGrid()
{
}

CXCHessianGrid&
CXCHessianGrid::operator=(const CXCHessianGrid& source)
{
    if (this == &source) return *this;
    
    _densityGridType = source._densityGridType;
    
    _xcGridType = source._xcGridType;
    
    _xcValues = source._xcValues;
    
    return *this;
}

CXCHessianGrid&
CXCHessianGrid::operator=(CXCHessianGrid&& source) noexcept
{
    if (this == &source) return *this;
    
    _densityGridType = std::move(source._densityGridType);
    
    _xcGridType = std::move(source._xcGridType);
    
    _xcValues = std::move(source._xcValues);
    
    return *this;
}

bool
CXCHessianGrid::operator==(const CXCHessianGrid& other) const
{
    if (_densityGridType != other._densityGridType) return false;
    
    if (_xcGridType != other._xcGridType) return false;
    
    if (_xcValues != other._xcValues) return false;
    
    return true;
}

bool
CXCHessianGrid::operator!=(const CXCHessianGrid& other) const
{
    return !(*this == other);
}

void
CXCHessianGrid::zero()
{
    _xcValues.zero();
}

int32_t
CXCHessianGrid::getNumberOfGridPoints() const
{
    return _xcValues.size(0);
}

const double*
CXCHessianGrid::xcHessianValues(const xcvars iVariableType,
                                const xcvars jVariableType) const
{
    if (_xcGridType == xcfun::lda)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(0);
                
                if (_densityGridType == dengrid::lima) return nullptr;
                
                if (_densityGridType == dengrid::limb) return _xcValues.data(0);
            }
            
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(1);
                
                if (_densityGridType == dengrid::lima) return nullptr;
                
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
        
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(2);
                
                if (_densityGridType == dengrid::lima) return _xcValues.data(0);
                
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    }
    
    if (_xcGridType == xcfun::gga)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(0);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(0);
            }
        
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(1);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(2);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(1);
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(3);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(4);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(5);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(0);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(6);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(7);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(1);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(8);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::grada)
        {
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(9);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(2);
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(10);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(11);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::gradb)
        {
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(12);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(2);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(13);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::gradab)
        {
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(14);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    }
    
    return nullptr;
}

double*
CXCHessianGrid::xcHessianValues(const xcvars iVariableType,
                                const xcvars jVariableType)
{
    if (_xcGridType == xcfun::lda)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(0);
                
                if (_densityGridType == dengrid::lima) return nullptr;
                
                if (_densityGridType == dengrid::limb) return _xcValues.data(0);
            }
            
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(1);
                
                if (_densityGridType == dengrid::lima) return nullptr;
                
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
        
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(2);
                
                if (_densityGridType == dengrid::lima) return _xcValues.data(0);
                
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    }
    
    if (_xcGridType == xcfun::gga)
    {
        if (iVariableType == xcvars::rhoa)
        {
            if (jVariableType == xcvars::rhoa)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(0);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(0);
            }
        
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(1);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(2);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(1);
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(3);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(4);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::rhob)
        {
            if (jVariableType == xcvars::rhob)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(5);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(0);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(6);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(7);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(1);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(8);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::grada)
        {
            if (jVariableType == xcvars::grada)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(9);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return _xcValues.data(2);
            }
        
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(10);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(11);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::gradb)
        {
            if (jVariableType == xcvars::gradb)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(12);
            
                if (_densityGridType == dengrid::lima) return _xcValues.data(2);
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(13);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    
        if (iVariableType == xcvars::gradab)
        {
            if (jVariableType == xcvars::gradab)
            {
                if (_densityGridType == dengrid::ab) return _xcValues.data(14);
            
                if (_densityGridType == dengrid::lima) return nullptr;
            
                if (_densityGridType == dengrid::limb) return nullptr;
            }
        }
    }
    
    return nullptr;
}

std::ostream&
operator<<(std::ostream& output, const CXCHessianGrid& source)
{
    output << std::endl;
    
    output << "[CXCHessianGrid (Object):" << &source << "]" << std::endl;
    
    output << "_densityGridType: " << to_string(source._densityGridType) << std::endl;
    
    output << "_xcGridType: " << to_string(source._xcGridType) << std::endl;
    
    output << "_xcValues: " << std::endl;
    
    output << source._xcValues << std::endl;
    
    return output;
}
