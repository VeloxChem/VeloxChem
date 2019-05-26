//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoCentersRecursionFunctions.hpp"

namespace t2crecfunc { // t2crecfunc namespace
    
    std::vector<CRecursionTerm>
    obRecursionForOverlap(const CRecursionTerm& recursionTerm)
    {
        std::vector<CRecursionTerm> obvec;
        
        auto rterm = recursionTerm;
        
        if (recursionTerm.isBraOfZeroOrder())
        {
            auto k1term = rterm.ketShift(-1, 0);
            
            if (k1term.isValid())
            {
                obvec.push_back(k1term);
            
                auto k2term = k1term.ketShift(-1, 0);
            
                if (k2term.isValid()) obvec.push_back(k2term);
            }
        }
        else
        {
            auto b1term = rterm.braShift(-1, 0);
            
            if (b1term.isValid())
            {
                obvec.push_back(b1term);
            
                auto b2term = b1term.braShift(-1, 0);
            
                if (b2term.isValid()) obvec.push_back(b2term);

                auto k1term = b1term.ketShift(-1, 0);
                
                if (k1term.isValid()) obvec.push_back(k1term);
            }
        }
        
        return obvec;
    }
    
    std::vector<CRecursionTerm>
    obRecursionForKineticEnergy(const CRecursionTerm& recursionTerm)
    {
        std::vector<CRecursionTerm> obvec;
        
        auto rterm = recursionTerm;
        
        if (recursionTerm.isBraOfZeroOrder())
        {
            auto k1term = rterm.ketShift(-1, 0);
            
            if (k1term.isValid())
            {
                obvec.push_back(k1term);
                
                auto k2term = k1term.ketShift(-1, 0);
                
                if (k2term.isValid())
                {
                    obvec.push_back(k2term);
                    
                    k2term.setLabel("Overlap");
                    
                    obvec.push_back(k2term); 
                }
            }
        }
        else
        {
            auto b1term = rterm.braShift(-1, 0);
            
            if (b1term.isValid())
            {
                obvec.push_back(b1term);
                
                auto b2term = b1term.braShift(-1, 0);
                
                if (b2term.isValid())
                {
                    obvec.push_back(b2term);
                    
                    b2term.setLabel("Overlap");
                    
                    obvec.push_back(b2term);
                }
                
                auto k1term = b1term.ketShift(-1, 0);
                
                if (k1term.isValid()) obvec.push_back(k1term);
            }
        }
        
        rterm.setLabel({"Overlap"});
        
        if (rterm.isValid()) obvec.push_back(rterm);
        
        return obvec;
    }
    
} // t2crecfunc namespace
