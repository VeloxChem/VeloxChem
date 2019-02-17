//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ExcitationVector.hpp"

#include <utility>

CExcitationVector::CExcitationVector()

    : _excitationType(szblock::aa)
{
    
}

CExcitationVector::CExcitationVector(const szblock               excitationType,
                                     const std::vector<int32_t>& braIndexes,
                                     const std::vector<int32_t>& ketIndexes,
                                     const std::vector<double>&  coefficients)

    : _excitationType(excitationType)

    , _braIndexes(CMemBlock<int32_t>(braIndexes))

    , _ketIndexes(CMemBlock<int32_t>(ketIndexes))

    , _coefficents(CMemBlock<double>(coefficients))
{
    
}

CExcitationVector::CExcitationVector(const CExcitationVector& source)

    : _excitationType(source._excitationType)

    , _braIndexes(source._braIndexes)

    , _ketIndexes(source._ketIndexes)

    , _coefficents(source._coefficents)
{
    
}

CExcitationVector::CExcitationVector(CExcitationVector&& source) noexcept

    : _excitationType(std::move(source._excitationType))

    , _braIndexes(std::move(source._braIndexes))

    , _ketIndexes(std::move(source._ketIndexes))

    , _coefficents(std::move(source._coefficents))
{
    
}

CExcitationVector::~CExcitationVector()
{
    
}

CExcitationVector&
CExcitationVector::operator=(const CExcitationVector& source)
{
    if (this == &source) return *this;
    
    _excitationType = source._excitationType;
    
    _braIndexes = source._braIndexes;
    
    _ketIndexes = source._ketIndexes;
    
    _coefficents = source._coefficents;
    
    return *this;
}

CExcitationVector&
CExcitationVector::operator=(CExcitationVector&& source) noexcept
{
    if (this == &source) return *this;
    
    _excitationType = std::move(source._excitationType);
    
    _braIndexes = std::move(source._braIndexes);
    
    _ketIndexes = std::move(source._ketIndexes);
    
    _coefficents = std::move(source._coefficents);
    
    return *this;
}

bool
CExcitationVector::operator==(const CExcitationVector& other) const
{
    if (_excitationType != other._excitationType) return false;
    
    if (_braIndexes != other._braIndexes) return false;
    
    if (_ketIndexes != other._ketIndexes) return false;
    
    if (_coefficents != other._coefficents) return false;
    
    return true;
}

bool
CExcitationVector::operator!=(const CExcitationVector& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&  output,
           const CExcitationVector& source)
{
    output << std::endl;
    
    output << "[CExcitationVector (Object):" << &source << "]" << std::endl;
    
    output << "_excitationType: " << to_string(source._excitationType) << std::endl;
    
    output << "_braIndexes: " << source._braIndexes << std::endl;
    
    output << "_ketIndexes: " << source._ketIndexes <<  std::endl;
    
    output << "_coefficents: " << source._coefficents <<  std::endl;
    
    return output;
}
