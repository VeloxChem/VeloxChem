//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GtoPairsContainer.hpp"

#include <sstream>

#include "GtoContainer.hpp"
#include "StringFormat.hpp"

CGtoPairsContainer::CGtoPairsContainer()
{
    
}

CGtoPairsContainer::CGtoPairsContainer(const std::vector<CGtoPairsBlock>& gtoPairsBlocks)

    : _gtoPairsBlocks(gtoPairsBlocks)
{
    
}

CGtoPairsContainer::CGtoPairsContainer(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           threshold)
{
    CGtoContainer gtoblks(molecule, basis);
    
    for (int32_t i = 0;  i < gtoblks.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtoblks.getGtoBlock(i);
        
        CGtoPairsBlock bpairs(bgtos, threshold);
        
        if (!bpairs.empty()) _gtoPairsBlocks.push_back(bpairs);
        
        for (int32_t j = i + 1; j < gtoblks.getNumberOfGtoBlocks(); j++)
        {
            auto kgtos = gtoblks.getGtoBlock(j);
            
            CGtoPairsBlock bkpairs(bgtos, kgtos, threshold);
            
            if (!bpairs.empty()) _gtoPairsBlocks.push_back(bkpairs);
        }
    }
}

CGtoPairsContainer::CGtoPairsContainer(const CGtoPairsContainer& source)

    : _gtoPairsBlocks(source._gtoPairsBlocks)
{
    
}

CGtoPairsContainer::CGtoPairsContainer(CGtoPairsContainer&& source) noexcept

    : _gtoPairsBlocks(std::move(source._gtoPairsBlocks))
{
    
}

CGtoPairsContainer::~CGtoPairsContainer()
{
    
}

CGtoPairsContainer&
CGtoPairsContainer::operator=(const CGtoPairsContainer& source)
{
    if (this == &source) return *this;
    
    _gtoPairsBlocks = source._gtoPairsBlocks;
    
    return *this;
}

CGtoPairsContainer&
CGtoPairsContainer::operator=(CGtoPairsContainer&& source) noexcept
{
    if (this == &source) return *this;
    
    _gtoPairsBlocks = std::move(source._gtoPairsBlocks);
    
    return *this;
}

bool
CGtoPairsContainer::operator==(const CGtoPairsContainer& other) const
{
    if (_gtoPairsBlocks.size() != other._gtoPairsBlocks.size()) return false;
    
    for (size_t i = 0; i < _gtoPairsBlocks.size(); i++)
    {
        if (_gtoPairsBlocks[i] != other._gtoPairsBlocks[i]) return false;
    }
    
    return true;
}

bool
CGtoPairsContainer::operator!=(const CGtoPairsContainer& other) const
{
    return !(*this == other);
}

CGtoPairsContainer
CGtoPairsContainer::split(const int32_t nodes) const
{
    std::vector<CGtoPairsBlock> ppvec;
    
    for (size_t i = 0; i < _gtoPairsBlocks.size(); i++)
    {
        auto cvec = _gtoPairsBlocks[i].split(nodes);
        
        for (size_t j = 0; j < cvec.size(); j++)
        {
            ppvec.push_back(cvec[j]);
        }
    }
    
   return CGtoPairsContainer(ppvec);
}

int32_t
CGtoPairsContainer::getNumberOfGtoPairsBlocks() const
{
    return static_cast<int32_t>(_gtoPairsBlocks.size());
}

CGtoPairsBlock
CGtoPairsContainer::getGtoPairsBlock(const int32_t iBlock) const
{
    return _gtoPairsBlocks[iBlock]; 
}

std::string
CGtoPairsContainer::printScreeningInfo() const
{
    std::string str("Contracted and Primitive GTO Pairs Screening: ");

    std::stringstream ss;
    
    ss << fstr::format(str, 80, fmt::left) << "\n\n";
    
    for (size_t i = 0; i < _gtoPairsBlocks.size(); i++)
    {
        std::string str(_gtoPairsBlocks[i].getPairType());
        
        str.append(" : ");
        
        str.append(_gtoPairsBlocks[i].getRawSizeString());
        
        str.append(" => ");
        
        str.append(_gtoPairsBlocks[i].getScreenedSizeString());
        
        ss << fstr::format(str, 80, fmt::left) << "\n";
    }

    return ss.str();
}

std::ostream&
operator<<(      std::ostream&       output,
           const CGtoPairsContainer& source)
{
    output << std::endl;
    
    output << "[CGtoPairsContainer (Object):" << &source << "]" << std::endl;
    
    output << "_gtoPairsBlocks: " << std::endl;
    
    for (size_t i = 0; i < source._gtoPairsBlocks.size(); i++)
    {
        output << "_gtoPairsBlocks[" << i << "]: "<< std::endl;
        
        output << source._gtoPairsBlocks[i] << std::endl;
    }
    
    return output;
}
