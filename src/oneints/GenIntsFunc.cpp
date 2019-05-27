//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GenIntsFunc.hpp"

namespace gintsfunc { // gintsfunc namespace
    
    CRecursionMap
    genRecursionMap(const CRecursionTerm&          recursionTerm,
                    const recblock                 angularForm,
                    const CRecursionFunctionsList& recursionFunctionsList)
    {
        CRecursionMap recmap(angularForm);
        
        if (recursionTerm.isValid()) recmap.add(recursionTerm); 
        
        std::vector<CRecursionTerm> refterms({recursionTerm});
        
        while (!refterms.empty())
        {
            std::vector<CRecursionTerm> genterms;
            
            // generate new recursion terms
            
            for (size_t i = 0; i < refterms.size(); i++)
            {
                auto curterms = recursionFunctionsList.compute(refterms[i]);
                
                if (!curterms.empty())
                {
                    genterms.insert(genterms.cbegin(), curterms.cbegin(),
                                    curterms.cend());
                }
            }
            
            // add recursion terms to recursion map
            
            recmap.append(genterms);
            
            // reset reference terms
            
            refterms = genterms;
        }
        
        return recmap;
    }
    
    CRecursionTerm
    genIntegral(const std::string& labelOfOperator,
                const int32_t      braAngularMomentum,
                const int32_t      ketAngularMomentum,
                const int32_t      ordderOfOperator)
    {
        CFourIndexes bang(braAngularMomentum, -1, -1, -1);
        
        CFourIndexes kang(ketAngularMomentum, -1, -1, -1);
        
        // set number of operator component according to operator label
        
        int32_t opcomp = 0;
        
        if (labelOfOperator == std::string("Electric Field")) opcomp = 1;
        
        if (labelOfOperator == std::string("Electric Field Gradient")) opcomp = 2;
        
        return CRecursionTerm(labelOfOperator, opcomp, true, bang, kang, 1, 1, ordderOfOperator);;
    }
    
} // gintsfunc namespace
