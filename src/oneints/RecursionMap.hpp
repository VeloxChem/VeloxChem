//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef RecursionMap_hpp
#define RecursionMap_hpp

#include <cstdlib>
#include <vector>

#include "RecursionTerm.hpp"

/**
 Class CRecursionMap stores recursion map for primitiveintegral
 and provides set of methods for manipulating with recursion map data.
 
 @author Z. Rinkevicius
 */
class CRecursionMap
{
    /**
     The vector of recursion term objects.
     */
    std::vector<CRecursionTerm> _recursionTerms;
    
    /**
     The vector of recursion indexes.
     */
    std::vector<int32_t> _recursionIndexes;
    
public:
    
    /**
     Creates an empty recursion map object.
     */
    CRecursionMap();
    
    /**
     Creates a recursion map object from vector of recursion term objects.
     
     @param recursionTerms the vector of recursion term objects.
     @param recursionIndexes the vector of recursion indexes.
     */
    CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                  const std::vector<int32_t>&        recursionIndexes);
    
    /**
     Creates a recursion map object by copying other recursion map object.
     
     @param source the recursion map object.
     */
    CRecursionMap(const CRecursionMap& source);
    
    /**
     Creates a recursion map object by moving other recursion map object.
     
     @param source the recursion map object.
     */
    CRecursionMap(CRecursionMap&& source) noexcept;
    
    /**
     Destroys a recursion map object.
     */
    ~CRecursionMap();
    
    /**
     Assigns a recursion map object by copying other recursion map object.
     
     @param source the recursion map object.
     */
    CRecursionMap& operator=(const CRecursionMap& source);
    
    /**
     Assigns a recursion map object by moving other recursion map object.
     
     @param source the recursion map object.
     */
    CRecursionMap& operator=(CRecursionMap&& source) noexcept;
    
    /**
     Compares recursion map object with other recursion map object.
     
     @param other the recursion map object.
     @return true if recursion map objects are equal, false otherwise.
     */
    bool operator==(const CRecursionMap& other) const;
    
    /**
     Compares recursion map object with other recursion map object.
     
     @param other the recursion term object.
     @return true if recursion map objects are not equal, false otherwise.
     */
    bool operator!=(const CRecursionMap& other) const;
    
    /**
     Converts recursion map object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the recursion map object.
     */
    friend std::ostream& operator<<(      std::ostream&  output,
                                    const CRecursionMap& source);
};


#endif /* RecursionMap_hpp */
