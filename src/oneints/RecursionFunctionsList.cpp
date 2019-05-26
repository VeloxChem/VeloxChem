//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionFunctionsList.hpp"

CRecursionFunctionsList::CRecursionFunctionsList()

    : _recursionFunctions(std::vector<CRecursionFunction>())
{
    
}

CRecursionFunctionsList::CRecursionFunctionsList(const std::vector<CRecursionFunction>& recursionFunctions)

    : _recursionFunctions(recursionFunctions)
{
    
}

CRecursionFunctionsList::CRecursionFunctionsList(const CRecursionFunctionsList& source)

    : _recursionFunctions(source._recursionFunctions)
{
    
}

CRecursionFunctionsList::CRecursionFunctionsList(CRecursionFunctionsList&& source) noexcept

    : _recursionFunctions(std::move(source._recursionFunctions))
{
    
}

CRecursionFunctionsList::~CRecursionFunctionsList()
{
    
}

CRecursionFunctionsList&
CRecursionFunctionsList::operator=(const CRecursionFunctionsList& source)
{
    if (this == &source) return *this;
    
    _recursionFunctions = source._recursionFunctions;
    
    return *this;
}

CRecursionFunctionsList&
CRecursionFunctionsList::operator=(CRecursionFunctionsList&& source) noexcept
{
    if (this == &source) return *this;
    
    _recursionFunctions = std::move(source._recursionFunctions);
    
    return *this;
}

bool
CRecursionFunctionsList::operator==(const CRecursionFunctionsList& other) const
{
    if (this == &other) return true;
    
    if (_recursionFunctions.size() != other._recursionFunctions.size()) return false;
    
    for (size_t i = 0; i < _recursionFunctions.size(); i++)
    {
        if (_recursionFunctions[i] != other._recursionFunctions[i]) return false;
    }
    
    return true;
}

bool
CRecursionFunctionsList::operator!=(const CRecursionFunctionsList& other) const
{
    return !( (*this) == other);
}

void
CRecursionFunctionsList::add(const CRecursionFunction& recursionFunction)
{
    auto idx = find(recursionFunction.getLabel());
    
    if (idx == -1)
    {
        _recursionFunctions.push_back(recursionFunction);
    }
}

std::vector<CRecursionTerm>
CRecursionFunctionsList::compute(const CRecursionTerm& recursionTerm) const
{
    auto idx = find(recursionTerm.getLabel());
    
    if (idx != -1) return _recursionFunctions[idx].compute(recursionTerm);
    
    return std::vector<CRecursionTerm>();
}

int32_t
CRecursionFunctionsList::find(const std::string& label) const
{
    for (size_t i = 0; i < _recursionFunctions.size(); i++)
    {
        if (_recursionFunctions[i].isMatch(label)) return static_cast<int32_t>(i);
    }
    
    return -1;
}

std::ostream&
operator<<(      std::ostream& output,
           const CRecursionFunctionsList& source)
{
    output << std::endl;
    
    output << "[CRecursionFunctionsList (Object):" << &source << "]" << std::endl;
    
    output << "_recursionFunctions: " << std::endl;
    
    for (size_t i = 0; i < source._recursionFunctions.size(); i++)
    {
        output << "_recursionFunctions[" << i << "]: "<< std::endl;
        
        output << source._recursionFunctions[i] << std::endl;
    }
    
    return output;
}
