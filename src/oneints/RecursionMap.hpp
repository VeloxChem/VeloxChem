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
#include "RecursionBlock.hpp"

/**
 Class CRecursionMap stores recursion map for primitiveintegral
 and provides set of methods for manipulating with recursion map data.
 
 @author Z. Rinkevicius
 */
class CRecursionMap
{
    /**
     The vector of unique recursion term objects.
     */
    std::vector<CRecursionTerm> _recursionTerms;
    
    /**
     The vector of recursion indexes.
     */
    std::vector<int32_t> _recursionIndexes;
    
    /**
     The angular form of recursion term objects.
     */
    recblock _angularForm;
    
public:
    
    /**
     Creates an empty recursion map object.
     */
    CRecursionMap();
    
    /**
     Creates an empty recursion map object with defined angular form of
     recursion term objects.
     
     @param angularForm the angular form of recursion term objects.
     */
    CRecursionMap(const recblock angularForm);
    
    /**
     Creates a recursion map object from vector of recursion term objects.
     
     @param recursionTerms the vector of recursion term objects.
     @param recursionIndexes the vector of recursion indexes.
     @param angularForm the angular form of recursion term objects. 
     */
    CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                  const std::vector<int32_t>&        recursionIndexes,
                  const recblock                     angularForm);
    
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
     Adds valid unique recursion term to recursion map object.

     @param recursionTerm the recursion term.
     */
    void add(const CRecursionTerm& recursionTerm);
    
    /**
     Appends valid unique recursion objects from recursion map object to other
     recursion map object.
     
     @param source the recursion map object.
     */
    void append(const CRecursionMap& source);
    
    /**
     Gets total number of integral components overl all recursion term objects
     included in recursion map object.
     
     @return the number of components.
     */
    int32_t getNumberOfComponents() const;
    
    /**
     Finds recursion term object in recursion map object.

     @param recursionTerm the recursion term object.
     @return true if recursion term object is in recursion map objects,
             false otherwise.
     */
    bool find(const CRecursionTerm& recursionTerm) const;
    
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
