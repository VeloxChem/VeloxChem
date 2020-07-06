//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ScreeningContainer.hpp"


CScreeningContainer::CScreeningContainer()
{
    
}

CScreeningContainer::CScreeningContainer(const CVecMemBlock<double>& braQValues,
                                         const CVecMemBlock<double>& ketQValues,
                                         const CGtoPairsContainer&   braGtoPairsContainer,
                                         const CGtoPairsContainer&   ketGtoPairsContainer,
                                         const ericut                screeningScheme,
                                         const double                threshold)
{
    // determine symmetry of GTOs containers
    
    bool symbk = (braGtoPairsContainer == ketGtoPairsContainer);
    
    // dimensions of bra and ket GTOs pairs containers
    
    auto bdim = braGtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    auto kdim = ketGtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    // create list of screened objects
    
    for (int32_t i = (bdim - 1); i >= 0; i--)
    {
        auto bpairs = braGtoPairsContainer.getGtoPairsBlock(i);
        
        auto joff = (symbk) ? i : 0;
        
        for (int32_t j = (kdim - 1); j >= joff; j--)
        {
            auto kpairs = ketGtoPairsContainer.getGtoPairsBlock(j);
            
            _screeners.push_back(CCauchySchwarzScreener(braQValues[i],
                                                        ketQValues[j],
                                                        bpairs, kpairs,
                                                        screeningScheme,
                                                        threshold));
        }
    }
}

CScreeningContainer::CScreeningContainer(const CScreeningContainer& source)

    : _screeners(source._screeners)
{
    
}

CScreeningContainer::CScreeningContainer(CScreeningContainer&& source) noexcept

    : _screeners(std::move(source._screeners))
{
    
}

CScreeningContainer::~CScreeningContainer()
{
    
}

CScreeningContainer&
CScreeningContainer::operator=(const CScreeningContainer& source)
{
    if (this == &source) return *this;
    
    _screeners = source._screeners;
    
    return *this;
}

CScreeningContainer&
CScreeningContainer::operator=(CScreeningContainer&& source) noexcept
{
    if (this == &source) return *this;
    
    _screeners = std::move(source._screeners);
    
    return *this;
}

bool
CScreeningContainer::operator==(const CScreeningContainer& other) const
{
    if (_screeners.size() != other._screeners.size()) return false;
    
    for (size_t i = 0; i < _screeners.size(); i++)
    {
        if (_screeners[i] != other._screeners[i]) return false;
    }
    
    return true;
}

bool
CScreeningContainer::operator!=(const CScreeningContainer& other) const
{
    return !(*this == other);
}

void
CScreeningContainer::setThreshold(const double threshold)
{
    for (size_t i = 0; i < _screeners.size(); i++)
    {
        _screeners[i].setThreshold(threshold); 
    }
}

bool
CScreeningContainer::isEmpty() const
{
    return _screeners.empty(); 
}

int32_t
CScreeningContainer::getNumberOfScreeners() const
{
    return static_cast<int32_t>(_screeners.size());
}

CCauchySchwarzScreener
CScreeningContainer::getScreener(const int32_t iScreener) const
{
    if (iScreener < getNumberOfScreeners())
    {
        return _screeners[iScreener];
    }
    
    return CCauchySchwarzScreener();
}

std::ostream&
operator<<(      std::ostream&        output,
           const CScreeningContainer& source)
{
    output << std::endl;
    
    output << "[CScreeningContainer (Object):" << &source << "]" << std::endl;
    
    output << "_screeners: " << std::endl;
    
    for (size_t i = 0; i < source._screeners.size(); i++)
    {
        output << "_screeners[" << i << "]: "<< std::endl;
        
        output << source._screeners[i] << std::endl;
    }
    
    return output;
}
