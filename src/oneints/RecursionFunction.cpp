//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionFunction.hpp"

CRecursionFunction::CRecursionFunction()

    : _label(std::string())

    , _funcAction(nullptr)
{
    
}

CRecursionFunction::CRecursionFunction(const std::string&                     label,
                                       const std::function<def_rec_func_typ>& funcAction)

    : _label(label)

    , _funcAction(funcAction)
{
    
}

CRecursionFunction::CRecursionFunction(const CRecursionFunction& source)

    : _label(source._label)

    , _funcAction(source._funcAction)
{
    
}

CRecursionFunction::CRecursionFunction(CRecursionFunction&& source) noexcept

    : _label(std::move(source._label))

    , _funcAction(std::move(source._funcAction))
{
    
}

CRecursionFunction::~CRecursionFunction()
{
    
}

CRecursionFunction&
CRecursionFunction::operator=(const CRecursionFunction& source)
{
    if (this == &source) return *this;
    
    _label = source._label;
    
    _funcAction = source._funcAction;
    
    return *this;
}

CRecursionFunction&
CRecursionFunction::operator=(CRecursionFunction&& source) noexcept
{
    if (this == &source) return *this;
    
    _label = std::move(source._label);
    
    _funcAction = std::move(source._funcAction);
    
    return *this;
}

bool
CRecursionFunction::operator==(const CRecursionFunction& other) const
{
    if (this == &other) return true;
    
    if (_label != other._label) return false;
    
    return true;
}

bool
CRecursionFunction::operator!=(const CRecursionFunction& other) const
{
    return !( (*this) == other);
}

std::vector<CRecursionTerm>
CRecursionFunction::compute(const CRecursionTerm& recursionTerm) const
{
    return _funcAction(recursionTerm);
}

std::string
CRecursionFunction::getLabel() const
{
    return _label;
}

bool
CRecursionFunction::isMatch(const std::string label) const
{
    if (label != _label) return false;
    
    return true;
}

std::ostream&
operator<<(      std::ostream&       output,
           const CRecursionFunction& source)
{
    output << std::endl;
    
    output << "[CRecFunction (Object):" << &source << "]" << std::endl;
    
    output << "_label: " << source._label <<  std::endl;
    
    output << "_funcAction: " << &(source._funcAction) << std::endl;
    
    return output;
}
