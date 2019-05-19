//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionMap.hpp"

CRecursionMap::CRecursionMap()
{
    _angularForm = recblock::cc;
}

CRecursionMap::CRecursionMap(const recblock angularForm)
{
    _angularForm = angularForm;
}

CRecursionMap::CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                             const std::vector<int32_t>&        recursionIndexes,
                             const recblock                     angularForm)

    : _recursionTerms(recursionTerms)

    , _recursionIndexes(recursionIndexes)

    , _angularForm(angularForm)
{
    
}

CRecursionMap::CRecursionMap(const CRecursionMap& source)

    : _recursionTerms(source._recursionTerms)

    , _recursionIndexes(source._recursionIndexes)

    , _angularForm(source._angularForm)
{
    
}

CRecursionMap::CRecursionMap(CRecursionMap&& source) noexcept

    : _recursionTerms(std::move(source._recursionTerms))

    , _recursionIndexes(std::move(source._recursionIndexes))

    , _angularForm(std::move(source._angularForm))
{
    
}

CRecursionMap::~CRecursionMap()
{
    
}

CRecursionMap&
CRecursionMap::operator=(const CRecursionMap& source)
{
    if (this == &source) return *this;
    
    _recursionTerms = source._recursionTerms;
    
    _recursionIndexes = source._recursionIndexes;
    
    _angularForm = source._angularForm;
    
    return *this;
}

CRecursionMap&
CRecursionMap::operator=(CRecursionMap&& source) noexcept
{
    if (this == &source) return *this;
    
    _recursionTerms = std::move(source._recursionTerms);
    
    _recursionIndexes = std::move(source._recursionIndexes);
    
    _angularForm = std::move(source._angularForm);
    
    return *this;
}

bool
CRecursionMap::operator==(const CRecursionMap& other) const
{
    if (_recursionTerms.size() != other._recursionTerms.size()) return false;
    
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (_recursionTerms[i] != other._recursionTerms[i]) return false;
    }
    
    if (_recursionIndexes.size() != other._recursionIndexes.size()) return false;
    
    for (size_t i = 0; i < _recursionIndexes.size(); i++)
    {
        if (_recursionIndexes[i] != other._recursionIndexes[i]) return false;
    }
    
    if (_angularForm != other._angularForm) return false;
    
    return true;
}

bool
CRecursionMap::operator!=(const CRecursionMap& other) const
{
    return !(*this == other);
}

void
CRecursionMap::add(const CRecursionTerm& recursionTerm)
{
    if (recursionTerm.isValid())
    {
        if (!find(recursionTerm))
        {
            auto ncomps = getNumberOfComponents();
        
            _recursionTerms.push_back(recursionTerm);
        
            _recursionIndexes.push_back(ncomps);
        }
    }
}

void
CRecursionMap::append(const CRecursionMap& source)
{
    if (_angularForm == source._angularForm)
    {
        for (size_t i = 0; i < source._recursionTerms.size(); i++)
        {
           add(source._recursionTerms[i]);
        }
    }
}

int32_t
CRecursionMap::getNumberOfComponents() const
{
    int32_t ncomps = 0;
    
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        ncomps += _recursionTerms[i].getNumberOfComponents(_angularForm);
    }
    
    return ncomps;
}

int32_t
CRecursionMap::getNumberOfTerms() const
{
    return static_cast<int32_t>(_recursionTerms.size()); 
}

CRecursionTerm
CRecursionMap::getTerm(const int32_t iRecursionTerm) const
{
    if (iRecursionTerm < getNumberOfComponents())
    {
        return _recursionTerms[iRecursionTerm];
    }
    
    return CRecursionTerm(); 
}

int32_t
CRecursionMap::getIndexOfTerm(const CRecursionTerm& recursionTerm) const
{
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (recursionTerm == _recursionTerms[i])
        {
            return _recursionIndexes[i]; 
        }
    }
    
    return -1;
}

int32_t
CRecursionMap::getMaxOrder(const std::string&  label,
                           const CFourIndexes& braAngularMomentum,
                           const CFourIndexes& ketAngularMomentum,
                           const int32_t       braCenters,
                           const int32_t       ketCenters) const
{
    int32_t mord = -1;
    
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (_recursionTerms[i].isIntegral(label, braAngularMomentum,
                                          ketAngularMomentum, braCenters,
                                          ketCenters))
        {
            auto cord = _recursionTerms[i].getOrder();
            
            if (cord > mord) mord = cord;
        }
    }
    
    return mord;
}

bool
CRecursionMap::find(const CRecursionTerm& recursionTerm) const
{
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (recursionTerm == _recursionTerms[i]) return true;
    }
    
    return false;
}

std::ostream&
operator<<(      std::ostream&  output,
           const CRecursionMap& source)
{
    output << std::endl;
    
    output << "[CRecursionMap (Object):" << &source << "]" << std::endl;
    
    output << "_recursionTerms: " << std::endl;
    
    for (size_t i = 0; i < source._recursionTerms.size(); i++)
    {
        output << "_recursionTerms[" << i << "]: "<< std::endl;
        
        output << source._recursionTerms[i] << std::endl;
    }

    output << "_recursionIndexes: " << std::endl;
    
    for (size_t i = 0; i < source._recursionIndexes.size(); i++)
    {
        output << "_recursionIndexes[" << i << "]: "<< std::endl;
        
        output << source._recursionIndexes[i] << std::endl;
    }
    
    output << "_angularForm: " << to_string(source._angularForm) << std::endl;
    
    return output;
}
