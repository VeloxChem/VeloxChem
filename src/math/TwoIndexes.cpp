//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoIndexes.hpp"

#include <utility>

CTwoIndexes::CTwoIndexes()

    : _iIndex(-1)

    , _jIndex(-1)
{

}

CTwoIndexes::CTwoIndexes(const int32_t iIndex,
                         const int32_t jIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)
{

}

CTwoIndexes::CTwoIndexes(const CTwoIndexes& source)

    : _iIndex(source._iIndex)

    , _jIndex(source._jIndex)
{

}

CTwoIndexes::CTwoIndexes(CTwoIndexes&& source) noexcept

    : _iIndex(std::move(source._iIndex))

    , _jIndex(std::move(source._jIndex))
{

}

CTwoIndexes::~CTwoIndexes()
{

}

CTwoIndexes&
CTwoIndexes::operator=(const CTwoIndexes& source)
{
    if (this == &source) return *this;

    _iIndex = source._iIndex;

    _jIndex = source._jIndex;

    return *this;
}

CTwoIndexes&
CTwoIndexes::operator=(CTwoIndexes&& source) noexcept
{
    if (this == &source) return *this;

    _iIndex = std::move(source._iIndex);

    _jIndex = std::move(source._jIndex);

    return *this;
}

bool
CTwoIndexes::operator==(const CTwoIndexes& other) const
{
    if (this == &other) return true;

    if (_iIndex != other._iIndex) return false;
    
    if (_jIndex != other._jIndex) return false;

    return true;
}

bool
CTwoIndexes::operator!=(const CTwoIndexes& other) const
{

    return !( (*this) == other);
}

int32_t
CTwoIndexes::first() const
{
    return _iIndex;
}

int32_t
CTwoIndexes::second() const
{
    return _jIndex;
}

bool
CTwoIndexes::isValidPair() const
{
    if (_iIndex < 0) return false;
    
    if (_jIndex < 0) return false; 
    
    return true;
}

std::ostream&
operator<<(      std::ostream& output,
           const CTwoIndexes&  source)
{
    output << std::endl;
    
    output << "[CTwoIndexes (Object):" << &source << "]" << std::endl;
    
    output << "_iIndex: " << source._iIndex <<  std::endl;
    
    output << "_jIndex: " << source._jIndex <<  std::endl;
    
    return output;
}
