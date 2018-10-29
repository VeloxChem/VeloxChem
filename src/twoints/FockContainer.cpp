//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "FockContainer.hpp"

CFockContainer::CFockContainer()

    : _subFockMatrices({})
{
    
}

CFockContainer::CFockContainer(const std::vector<CFockSubMatrix>& subFockMatrices)

    : _subFockMatrices(subFockMatrices)
{
    
}

CFockContainer::CFockContainer(const CAOFockMatrix*  fockMatrix,
                               const CGtoPairsBlock& braGtoPairsBlock,
                               const CGtoPairsBlock& ketGtoPairsBlock)
{
    auto nfock = fockMatrix->getNumberOfFockMatrices();
    
    for (int32_t i = 0; i < nfock; i++)
    {
        _subFockMatrices.push_back(CFockSubMatrix(braGtoPairsBlock,
                                                  ketGtoPairsBlock,
                                                  fockMatrix->getFockType(i)));
    }
}

CFockContainer::CFockContainer(const CFockContainer& source)

    : _subFockMatrices(source._subFockMatrices)
{
    
}

CFockContainer::CFockContainer(CFockContainer&& source) noexcept

    : _subFockMatrices(std::move(source._subFockMatrices))
{
    
}

CFockContainer::~CFockContainer()
{
    
}

CFockContainer&
CFockContainer::operator=(const CFockContainer& source)
{
    if (this == &source) return *this;
    
    _subFockMatrices = source._subFockMatrices;
    
    return *this;
}

CFockContainer&
CFockContainer::operator=(CFockContainer&& source) noexcept
{
    if (this == &source) return *this;
    
    _subFockMatrices = std::move(source._subFockMatrices);
    
    return *this;
}

bool
CFockContainer::operator==(const CFockContainer& other) const
{
    if (_subFockMatrices.size() != other._subFockMatrices.size())
    {
        return false;
    }
    
    for (size_t i = 0; i < _subFockMatrices.size(); i++)
    {
        if (_subFockMatrices[i] != other._subFockMatrices[i]) return false;
    }
    
    return true;
}

bool
CFockContainer::operator!=(const CFockContainer& other) const
{
    return !(*this == other);
}

double*
CFockContainer::getSubMatrixData(const int32_t iFockMatrix,
                                 const int32_t iSubMatrix,
                                 const int32_t iComponent)
{
    return _subFockMatrices[iFockMatrix].getSubMatrixData(iSubMatrix,
                                                          iComponent);
}

const int32_t*
CFockContainer::getStartPositionsA(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getStartPositionsA();
}

const int32_t*
CFockContainer::getStartPositionsB(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getStartPositionsB();
}

const int32_t*
CFockContainer::getStartPositionsC(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getStartPositionsC();
}

const int32_t*
CFockContainer::getStartPositionsD(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getStartPositionsD();
}

int32_t
CFockContainer::getDimensionsA(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getDimensionsA();
}

int32_t
CFockContainer::getDimensionsB(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getDimensionsB();
}

int32_t
CFockContainer::getDimensionsC(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getDimensionsC();
}

int32_t
CFockContainer::getDimensionsD(const int32_t iFockMatrix) const
{
    return _subFockMatrices[iFockMatrix].getDimensionsD();
}

std::ostream&
operator<<(      std::ostream&   output,
           const CFockContainer& source)
{
    output << std::endl;
    
    output << "[CFockContainer (Object):" << &source << "]" << std::endl;
    
    output << "_subFockMatrices: " << std::endl;
    
    for (size_t i = 0; i < source._subFockMatrices.size(); i++)
    {
        output << "_subFockMatrices[" << i << "]: " << std::endl;
        
        output << source._subFockMatrices[i] << std::endl;
    }
    
    return output;
}
