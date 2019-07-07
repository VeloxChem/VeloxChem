//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef XCFunctional_hpp
#define XCFunctional_hpp

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>

#include "XCFuncType.hpp"
#include "PrimitiveFunctional.hpp"

/**
 Class CXCFunctional stores information about exchange-correlation functional and provides methods
 to perform actions with stored exchange-correlation functional.
 
 @author Z. Rinkevicius
 */
class CXCFunctional
{
    /**
     The label of exchange-correlation functional.
     */
    std::string _label;
    
    /**
     The type of exchange-correlation functional.
     */
    xcfun _xcFuncType;
    
    
    /**
     The fraction of exact Hatree-Fock exchange in functional.
     */
    double _fractionOfExactExchange;
    
    /**
     The vector of primitive exchange-correlation functionals.
     */
    std::vector<CPrimitiveFunctional> _primitiveFunctionals;
    
public:
    /**
     Creates an empty exchange-correlation functional object.
     */
    CXCFunctional();
    
    /**
     Creates a exchange-correlation functional object.
     
     @param label the label of primitive exchange-correlation functional.
     @param xcFuncType the type of primitive exchange-correlation functional.
     @param fractionOfExactExchange the fraction of exact Hatree-Fock exchange.
     
     */
    CXCFunctional(const std::string&                       label,
                  const xcfun                              xcFuncType,
                  const double                             fractionOfExactExchange,
                  const std::vector<CPrimitiveFunctional>& primitiveFunctionals);
    
    /**
     Creates a exchange-correlation functional object by copying other exchange-correlation functional object.
     
     @param source the recursion term object.
     */
    CXCFunctional(const CXCFunctional& source);
    
    /**
     Creates a exchange-correlation functional object by moving other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional(CXCFunctional&& source) noexcept;
    
    /**
     Destroys a exchange-correlation functional object.
     */
    ~CXCFunctional();
    
    /**
     Assigns a exchange-correlation functional object by copying other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional& operator=(const CXCFunctional& source);
    
    /**
     Assigns a exchange-correlation functional object by moving other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional& operator=(CXCFunctional&& source) noexcept;
    
    /**
     Compares exchange-correlation functional object with other exchange-correlation functional object.
     
     @param other the exchange-correlation functional object.
     @return true if exchange-correlation functional objects are equal, false otherwise.
     */
    bool operator==(const CXCFunctional& other) const;
    
    /**
     Compares exchange-correlation functional object with other exchange-correlation functional object.
     
     @param other the exchange-correlation functional object.
     @return true if exchange-correlation functional objects are not equal, false otherwise.
     */
    bool operator!=(const CXCFunctional& other) const;
    
    /**
     Gets label of exchange-correlation functional.
     
     @return the label.
     */
    std::string getLabel() const;
    
    /**
     Gets type of exchange-correlation functional.
     
     @return the type of exchange-correlation functional.
     */
    xcfun getFunctionalType() const;
    
    /**
     Gets fraction of exact Hatree-Fock exchange in exchange-correlation functional.

     @return the fraction of exact Hatree-Fock exchange.
     */
    double getFractionOfExactExchange() const;
    
    /**
     Determines if exchange-correlation functional is of hybrid type i.e. non-zero fraction of exact Hatree-Fock
     exchange.

     @return true if hybrid exchange-correlation functional, false otherwise.
     */
    bool isHybridFunctional() const;
    
    /**
     Computes first derivative of exchange-correlation functional for given density grid.
     
     @param xcGradientGrid the exchange-correlation gradient grid object.
     @param densityGrid the density grid object.
     */
    void compute(      CXCGradientGrid& xcGradientGrid,
                 const CDensityGrid&    densityGrid) const;
    
    /**
     Converts exchange-correlation functional object to text format and insert it into output text stream.
     
     @param output the output text stream.
     @param source the exchange-correlation functional object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CXCFunctional& source);
};

#endif /* XCFunctional_hpp */
