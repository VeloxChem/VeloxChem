//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionTerm.hpp"

#include <utility>

CRecursionTerm::CRecursionTerm()

    : _labelOfOperator(std::string())

    , _orderOfOperator(-1)

    , _isReducedOperator(false)

    , _braAngularMomentum(CFourIndexes(-1, -1, -1, -1))

    , _ketAngularMomentum(CFourIndexes(-1, -1, -1, -1))

    , _braCenters(-1)

    , _ketCenters(-1)

    , _orderOfIntegral(-1)
{
    
}

CRecursionTerm::CRecursionTerm(const std::string&  labelOfOperator,
                               const int32_t       orderOfOperator,
                               const bool          isReducedOperator,
                               const CFourIndexes& braAngularMomentum,
                               const CFourIndexes& ketAngularMomentum,
                               const int32_t       braCenters,
                               const int32_t       ketCenters,
                               const int32_t       orderOfIntegral)

    : _labelOfOperator(labelOfOperator)

    , _orderOfOperator(orderOfOperator)

    , _isReducedOperator(isReducedOperator)

    , _braAngularMomentum(braAngularMomentum)

    , _ketAngularMomentum(ketAngularMomentum)

    , _braCenters(braCenters)

    , _ketCenters(ketCenters)

    , _orderOfIntegral(orderOfIntegral)
{
    
}

CRecursionTerm::CRecursionTerm(const CRecursionTerm& source)

    : _labelOfOperator(source._labelOfOperator)

    , _orderOfOperator(source._orderOfOperator)

    , _isReducedOperator(source._isReducedOperator)

    , _braAngularMomentum(source._braAngularMomentum)

    , _ketAngularMomentum(source._ketAngularMomentum)

    , _braCenters(source._braCenters)

    , _ketCenters(source._ketCenters)

    , _orderOfIntegral(source._orderOfIntegral)
{
    
}

CRecursionTerm::CRecursionTerm(CRecursionTerm&& source) noexcept

    : _labelOfOperator(std::move(source._labelOfOperator))

    , _orderOfOperator(std::move(source._orderOfOperator))

    , _isReducedOperator(std::move(source._isReducedOperator))

    , _braAngularMomentum(std::move(source._braAngularMomentum))

    , _ketAngularMomentum(std::move(source._ketAngularMomentum))

    , _braCenters(std::move(source._braCenters))

    , _ketCenters(std::move(source._ketCenters))

    , _orderOfIntegral(std::move(source._orderOfIntegral))
{
    
}

CRecursionTerm::~CRecursionTerm()
{
    
}

CRecursionTerm&
CRecursionTerm::operator=(const CRecursionTerm& source)
{
    if (this == &source) return *this;
    
    _labelOfOperator = source._labelOfOperator;
    
    _orderOfOperator = source._orderOfOperator;
    
    _isReducedOperator = source._isReducedOperator;
    
    _braAngularMomentum = source._braAngularMomentum;
    
    _ketAngularMomentum = source._ketAngularMomentum;
    
    _braCenters = source._braCenters;
    
    _ketCenters = source._ketCenters;
    
    _orderOfIntegral = source._orderOfIntegral;
    
    return *this;
}

CRecursionTerm&
CRecursionTerm::operator=(CRecursionTerm&& source) noexcept
{
    if (this == &source) return *this;
    
    _labelOfOperator = std::move(source._labelOfOperator);
    
    _orderOfOperator = std::move(source._orderOfOperator);
    
    _isReducedOperator = std::move(source._isReducedOperator);
    
    _braAngularMomentum = std::move(source._braAngularMomentum);
    
    _ketAngularMomentum = std::move(source._ketAngularMomentum);
    
    _braCenters = std::move(source._braCenters);
    
    _ketCenters = std::move(source._ketCenters);
    
    _orderOfIntegral = std::move(source._orderOfIntegral);
    
    return *this;
}

bool
CRecursionTerm::operator==(const CRecursionTerm& other) const
{
    if (_labelOfOperator != other._labelOfOperator) return false;
    
    if (_orderOfOperator != other._orderOfOperator) return false;
    
    if (_isReducedOperator != other._isReducedOperator) return false;
    
    if (_braAngularMomentum != other._braAngularMomentum) return false;
    
    if (_ketAngularMomentum != other._ketAngularMomentum) return false;
    
    if (_braCenters != other._braCenters) return false;
    
    if (_ketCenters != other._ketCenters) return false;
    
    if (_orderOfIntegral != other._orderOfIntegral) return false;
    
    return true;
}

bool
CRecursionTerm::operator!=(const CRecursionTerm& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&   output,
           const CRecursionTerm& source)
{
    output << std::endl;
    
    output << "[CRecursionTerm (Object):" << &source << "]" << std::endl;
    
    output << "_labelOfOperator: " << source._labelOfOperator << std::endl;
    
    output << "_orderOfOperator: " << source._orderOfOperator << std::endl;
    
    output << "_isReducedOperator: " << source._isReducedOperator <<  std::endl;
    
    output << "_braAngularMomentum: " << source._braAngularMomentum << std::endl;
    
    output << "_ketAngularMomentum: " << source._ketAngularMomentum << std::endl;
    
    output << "_braCenters: " << source._braCenters << std::endl;
    
    output << "_ketCenters: " << source._ketCenters << std::endl;
    
    output << "_orderOfIntegral: " << source._orderOfIntegral << std::endl;
    
    return output;
}

