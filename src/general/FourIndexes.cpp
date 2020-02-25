//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "FourIndexes.hpp"

#include <utility>

CFourIndexes::CFourIndexes()

    : _iIndex(-1)

    , _jIndex(-1)

    , _kIndex(-1)

    , _lIndex(-1)
{
}

CFourIndexes::CFourIndexes(const int32_t iIndex, const int32_t jIndex, const int32_t kIndex, const int32_t lIndex)

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
    return !((*this) == other);
}

void
CFourIndexes::shift(const int32_t shiftValue, const int32_t iComponent)
{
    if (iComponent == 0) _iIndex += shiftValue;

    if (iComponent == 1) _jIndex += shiftValue;

    if (iComponent == 2) _kIndex += shiftValue;

    if (iComponent == 3) _lIndex += shiftValue;
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

int32_t
CFourIndexes::value(const int32_t iComponent) const
{
    if (iComponent == 0) return _iIndex;

    if (iComponent == 1) return _jIndex;

    if (iComponent == 2) return _kIndex;

    if (iComponent == 3) return _lIndex;

    return -1;
}

std::ostream&
operator<<(std::ostream& output, const CFourIndexes& source)
{
    output << std::endl;

    output << "[CFourIndexes (Object):" << &source << "]" << std::endl;

    output << "_iIndex: " << source._iIndex << std::endl;

    output << "_jIndex: " << source._jIndex << std::endl;

    output << "_kIndex: " << source._kIndex << std::endl;

    output << "_lIndex: " << source._lIndex << std::endl;

    return output;
}
