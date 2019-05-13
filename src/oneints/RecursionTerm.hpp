//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef RecursionTerm_hpp
#define RecursionTerm_hpp

#include <cstdlib>
#include <string>

#include "FourIndexes.hpp"

/**
 Class CRecursionTerm stores meta data for primitive or auxilary integral
 and provides set of methods for manipulating with meta data.
 
 @author Z. Rinkevicius
 */
class CRecursionTerm
{
    /**
     The label of integrand operator.
     */
    std::string _labelOfOperator;
    
    /**
     The tensorial order of integrand operator.
     */
    int32_t _orderOfOperator;
    
    /**
     The flag indicating form, full or reduced, of integrand operator.
     */
    int32_t _isReducedOperator;
    
    /**
     The angular momentum of bra side.
     */
    CFourIndexes _braAngularMomentum;
    
    /**
     The angular momentum of ket side.
     */
    CFourIndexes _ketAngularMomentum;
    
    /**
     The number of centers on bra side.
     */
    int32_t _braCenters;
    
    /**
     The number of centers on ket side.
     */
    int32_t _ketCenters;
    
    /**
     The order of integral.
     */
    int32_t _orderOfIntegral;
    
    /**
     Cheks if angular momentum is valid for given number of centers.

     @param angularMomentum the angular momentum.
     @param nCenters the number of centers.
     @return true if angular momentum is valids, false otherwise.
     */
    bool _isValidAngularMomentum(const CFourIndexes& angularMomentum,
                                 const int32_t       nCenters) const;
    
public:
    
    /**
     Creates an empty recursion term object.
     */
    CRecursionTerm();
    
    /**
     Creates a recursion term object from given meta data.
     
     @param labelOfOperator the label of integrand operator.
     @param orderOfOperator the tensorial order of integrand operator.
     @param isReducedOperator the form of integrand operator: reduced or full
            tensor.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @param braCenters the number of GTO centers on bra side.
     @param ketCenters the number of GTO centers on ket side.
     @param orderOfIntegral the order of integral.
     */
   CRecursionTerm(const std::string&  labelOfOperator,
                  const int32_t       orderOfOperator,
                  const bool          isReducedOperator,
                  const CFourIndexes& braAngularMomentum,
                  const CFourIndexes& ketAngularMomentum,
                  const int32_t       braCenters,
                  const int32_t       ketCenters,
                  const int32_t       orderOfIntegral);
    
    /**
     Creates a recursion term object by copying other recursion term object.
     
     @param source the recursion term object.
     */
    CRecursionTerm(const CRecursionTerm& source);
    
    /**
     Creates a recursion term object by moving other recursion term object.
     
     @param source the recursion term object.
     */
    CRecursionTerm(CRecursionTerm&& source) noexcept;
    
    /**
     Destroys a recursion term object.
     */
    ~CRecursionTerm();
    
    /**
     Assigns a recursion term object by copying other recursion term object.
     
     @param source the recursion term object.
     */
    CRecursionTerm& operator=(const CRecursionTerm& source);
    
    /**
     Assigns a recursion term object by moving other recursion term object.
     
     @param source the recursion term object.
     */
    CRecursionTerm& operator=(CRecursionTerm&& source) noexcept;
    
    /**
     Compares recursion term object with other recursion term object.
     
     @param other the recursion term object.
     @return true if recursion term objects are equal, false otherwise.
     */
    bool operator==(const CRecursionTerm& other) const;
    
    /**
     Compares recursion term object with other recursion term object.
     
     @param other the recursion term object.
     @return true if recursion term objects are not equal, false otherwise.
     */
    bool operator!=(const CRecursionTerm& other) const;
    
    /**
     Checks if recursion term object is valid recursion term.

     @return true if recursion term object is valid recursion term, false otherwise.
     */
    bool isValid() const;
    
    /**
     Converts recursion term object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the recursion term object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CRecursionTerm& source);
};

#endif /* RecursionTerm_hpp */
