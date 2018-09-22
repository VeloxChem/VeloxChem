//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ScreeningContainer.hpp"


CScreeningContainer::CScreeningContainer()
{
    
}

CScreeningContainer::CScreeningContainer(const CGtoPairsContainer& braGtoPairsContainer,
                                         const CGtoPairsContainer& ketGtoPairsContainer,
                                         const ericut              screeningScheme,
                                         const double              threshold)
{
    // determine symmetry of GTOs containers
    
    bool symbk = (braGtoPairsContainer == ketGtoPairsContainer);
    
    // dimensions of bra and ket GTOs pairs containers
    
    auto bdim = braGtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    auto kdim = ketGtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    // create list of screened objects
    
    for (int32_t i = 0; i < bdim; i++)
    {
        auto bpairs = braGtoPairsContainer.getGtoPairsBlock(i);
        
        auto joff = (symbk) ? i : 0;
        
        for (int32_t j = joff; j < kdim; j++)
        {
            auto kpairs = ketGtoPairsContainer.getGtoPairsBlock(j);
            
            _screeners.push_back(CCauchySchwarzScreener(bpairs, kpairs,
                                                        screeningScheme,
                                                        threshold));
        }
    }
}

CScreeningContainer::~CScreeningContainer()
{
    
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
        output << "_braScreeners[" << i << "]: "<< std::endl;
        
        output << source._screeners[i] << std::endl;
    }
    
    return output;
}
