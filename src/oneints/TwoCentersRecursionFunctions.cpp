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
    
    std::vector<CRecursionTerm>
    obRecursionForNuclearPotential(const CRecursionTerm& recursionTerm)
    {
        std::vector<CRecursionTerm> obvec;
        
        auto rterm = recursionTerm;
        
        if (recursionTerm.isBraOfZeroOrder())
        {
            auto k1term = rterm.ketShift(-1, 0);
            
            if (k1term.isValid())
            {
                obvec.push_back(k1term);
                
                obvec.push_back(k1term.orderShift(1));
                
                auto k2term = k1term.ketShift(-1, 0);
                
                if (k2term.isValid())
                {
                    obvec.push_back(k2term);
                    
                    obvec.push_back(k2term.orderShift(1));
                }
            }
        }
        else
        {
            auto b1term = rterm.braShift(-1, 0);
            
            if (b1term.isValid())
            {
                obvec.push_back(b1term);
                
                obvec.push_back(b1term.orderShift(1));
                
                auto b2term = b1term.braShift(-1, 0);
                
                if (b2term.isValid())
                {
                    obvec.push_back(b2term);
                    
                    obvec.push_back(b2term.orderShift(1));
                }
                
                auto k1term = b1term.ketShift(-1, 0);
                
                if (k1term.isValid())
                {
                    obvec.push_back(k1term);
                    
                    obvec.push_back(k1term.orderShift(1));
                }
            }
        }
        
        // special case: (S|A|S) integral
        
        if (obvec.empty() && rterm.isValid())
        {
            rterm.setLabel("Overlap");
            
            obvec.push_back(rterm);
        }
        
        return obvec;
    }
    
    std::vector<CRecursionTerm>
    obRecursionForElectricField(const CRecursionTerm& recursionTerm)
    {
        std::vector<CRecursionTerm> obvec;
        
        auto rterm = recursionTerm;
        
        if (recursionTerm.isBraOfZeroOrder())
        {
            auto k1term = rterm.ketShift(-1, 0);
            
            if (k1term.isValid())
            {
                obvec.push_back(k1term);
                
                obvec.push_back(k1term.orderShift(1));
                
                auto k2term = k1term.ketShift(-1, 0);
                
                if (k2term.isValid())
                {
                    obvec.push_back(k2term);
                    
                    obvec.push_back(k2term.orderShift(1));
                }
                
                k1term.setLabel({"Nuclear Potential"});
                
                auto nterm = (k1term.orderShift(1)).operatorShift(-1);
                
                if (nterm.isValid()) obvec.push_back(nterm);
            }
        }
        else
        {
            auto b1term = rterm.braShift(-1, 0);
            
            if (b1term.isValid())
            {
                obvec.push_back(b1term);
                
                obvec.push_back(b1term.orderShift(1));
                
                auto b2term = b1term.braShift(-1, 0);
                
                if (b2term.isValid())
                {
                    obvec.push_back(b2term);
                    
                    obvec.push_back(b2term.orderShift(1));
                }
                
                auto k1term = b1term.ketShift(-1, 0);
                
                if (k1term.isValid())
                {
                    obvec.push_back(k1term);
                    
                    obvec.push_back(k1term.orderShift(1));
                }
                
                b1term.setLabel({"Nuclear Potential"});
                
                auto nterm = (b1term.orderShift(1)).operatorShift(-1);
                
                if (nterm.isValid()) obvec.push_back(nterm);
            }
        }
        
        // special case: (S|A(1)|S) integral
        
        if (obvec.empty() && rterm.isValid())
        {
            rterm.setLabel({"Nuclear Potential"});
            
            auto nterm = (rterm.orderShift(1)).operatorShift(-1);
            
            if (nterm.isValid()) obvec.push_back(nterm);
        }
        
        return obvec;
    }
    
    std::vector<CRecursionTerm>
    obRecursionForElectricFieldGradient(const CRecursionTerm& recursionTerm)
    {
        std::vector<CRecursionTerm> obvec;
        
        auto rterm = recursionTerm;
        
        if (recursionTerm.isBraOfZeroOrder())
        {
            auto k1term = rterm.ketShift(-1, 0);
            
            if (k1term.isValid())
            {
                obvec.push_back(k1term);
                
                obvec.push_back(k1term.orderShift(1));
                
                auto k2term = k1term.ketShift(-1, 0);
                
                if (k2term.isValid())
                {
                    obvec.push_back(k2term);
                    
                    obvec.push_back(k2term.orderShift(1));
                }
                
                k1term.setLabel({"Electric Field"});
                
                auto nterm = (k1term.orderShift(1)).operatorShift(-1);
                
                if (nterm.isValid()) obvec.push_back(nterm);
            }
        }
        else
        {
            auto b1term = rterm.braShift(-1, 0);
            
            if (b1term.isValid())
            {
                obvec.push_back(b1term);
                
                obvec.push_back(b1term.orderShift(1));
                
                auto b2term = b1term.braShift(-1, 0);
                
                if (b2term.isValid())
                {
                    obvec.push_back(b2term);
                    
                    obvec.push_back(b2term.orderShift(1));
                }
                
                auto k1term = b1term.ketShift(-1, 0);
                
                if (k1term.isValid())
                {
                    obvec.push_back(k1term);
                    
                    obvec.push_back(k1term.orderShift(1));
                }
                
                b1term.setLabel({"Electric Field"});
                
                auto nterm = (b1term.orderShift(1)).operatorShift(-1);
                
                if (nterm.isValid()) obvec.push_back(nterm);
            }
        }
        
        // special case: (S|A(2)|S) integral
        
        if (obvec.empty() && rterm.isValid())
        {
            rterm.setLabel({"Electric Field"});
            
            auto nterm = (rterm.orderShift(1)).operatorShift(-1);
            
            if (nterm.isValid()) obvec.push_back(nterm);
        }
        
        return obvec;
    }
    
} // t2crecfunc namespace
