//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "FourIndexes.hpp"

#include <utility>

CFourIndexes::CFourIndexes()

    : _iIndex(-1)

    , _jIndex(-1)

    , _kIndex(-1)

    , _lIndex(-1)
{
    
}

CFourIndexes::CFourIndexes(const int32_t iIndex,
                           const int32_t jIndex,
                           const int32_t kIndex,
                           const int32_t lIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)

    , _lIndex(lIndex)
{
    
}

CFourIndexes::CFourIndexes(const CFourIndexes& source)

    : _iIndex(source._iIndex)

    , _jIndex(source._jIndex)

    , _kIndex(source._kIndex)

    , _lIndex(source._lIndex)
{
    
}

CFourIndexes::CFourIndexes(CFourIndexes&& source) noexcept

    : _iIndex(std::move(source._iIndex))

    , _jIndex(std::move(source._jIndex))

    , _kIndex(std::move(source._kIndex))

    , _lIndex(std::move(source._lIndex))
{
    
}

CFourIndexes::~CFourIndexes()
{
    
}

CFourIndexes&
CFourIndexes::operator=(const CFourIndexes& source)
{
    if (this == &source) return *this;
    
    _iIndex = source._iIndex;
    
    _jIndex = source._jIndex;
    
    _kIndex = source._kIndex;
    
    _lIndex = source._lIndex;
    
    return *this;
}

CFourIndexes&
CFourIndexes::operator=(CFourIndexes&& source) noexcept
{
    if (this == &source) return *this;
    
    _iIndex = std::move(source._iIndex);
    
    _jIndex = std::move(source._jIndex);
    
    _kIndex = std::move(source._kIndex);
    
    _lIndex = std::move(source._lIndex);
    
    return *this;
}

bool
CFourIndexes::operator==(const CFourIndexes& other) const
{
    if (this == &other) return true;
    
    if (_iIndex != other._iIndex) return false;
    
    if (_jIndex != other._jIndex) return false;
    
    if (_kIndex != other._kIndex) return false;
    
    if (_lIndex != other._lIndex) return false;
    
    return true;
}

bool
CFourIndexes::operator!=(const CFourIndexes& other) const
{
    
    return !( (*this) == other);
}

int32_t
CFourIndexes::first() const
{
    return _iIndex;
}

int32_t
CFourIndexes::second() const
{
    return _jIndex;
}

int32_t
CFourIndexes::third() const
{
    return _kIndex;
}

int32_t
CFourIndexes::fourth() const
{
    return _lIndex;
}

bool
CFourIndexes::isValidQuadruple() const
{
    if (_iIndex < 0) return false;
    
    if (_jIndex < 0) return false;
    
    if (_kIndex < 0) return false;
    
    if (_lIndex < 0) return false;
    
    return true;
}

std::ostream&
operator<<(      std::ostream&  output,
           const CFourIndexes& source)
{
    output << std::endl;
    
    output << "[CFourIndexes (Object):" << &source << "]" << std::endl;
    
    output << "_iIndex: " << source._iIndex <<  std::endl;
    
    output << "_jIndex: " << source._jIndex <<  std::endl;
    
    output << "_kIndex: " << source._kIndex <<  std::endl;
    
    output << "_lIndex: " << source._lIndex <<  std::endl;
    
    return output;
}
