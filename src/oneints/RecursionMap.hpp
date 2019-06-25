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

#include "RecursionBlock.hpp"
#include "RecursionTerm.hpp"
#include "MemBlock2D.hpp"

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

    /**
     The number of repeated units in recursion term subspace.
     */
    int32_t _repeatUnits;

   public:
    /**
     Creates an empty recursion map object.
     */
    CRecursionMap();

    /**
     Creates an empty recursion map object with defined angular form of
     recursion term objects.

     @param angularForm the angular form of recursion term objects.
     @param repeatUnits the maximum number of repeated units in recursion term
            subspace.
     */
    CRecursionMap(const recblock angularForm, const int32_t repeatUnits);

    /**
     Creates a recursion map object from vector of recursion term objects.

     @param recursionTerms the vector of recursion term objects.
     @param recursionIndexes the vector of recursion indexes.
     @param angularForm the angular form of recursion term objects.
     @param repeatUnits the maximum number of repeated units in recursion term
            subspace.
     */
    CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                  const std::vector<int32_t>&        recursionIndexes,
                  const recblock                     angularForm,
                  const int32_t                      repeatUnits);

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
     Appends valid unique recursion objects from vector of recursion term
     objects to recursion map object.

     @param recursionTerms the vector of recursion term objects.
     */
    void append(const std::vector<CRecursionTerm>& recursionTerms);
    
    /**
     Gets raw index of specific recursion term object.

     @param recursionTerm the recursion term object.
     @return the raw index of recursion term.
     */
    int32_t index(const CRecursionTerm& recursionTerm) const;

    /**
     Gets total number of integral components overl all recursion term objects
     included in recursion map object.

     @return the number of components.
     */
    int32_t getNumberOfComponents() const;

    /**
     Gets number of recursion term objects in recursion map object.

     @return the number of recursion term objects.
     */
    int32_t getNumberOfTerms() const;

    /**
     Gets recursion term from vector of recursion term objects.

     @param iRecursionTerm the index of recursion term in vector of recursion
            term  objects.
     @return the recursion term object.
     */
    CRecursionTerm getTerm(const int32_t iRecursionTerm) const;

    /**
     Gets index of recursion term object in space of recursion term objects
     components.

     @param recursionTerm the recursion term object.
     @return the index of recursion term object.
     */
    int32_t getIndexOfTerm(const CRecursionTerm& recursionTerm) const;

    /**
     Determines maximum order of recursion term object with specified properties
     in vector of recursion term objects.

     @param label the label of integrand operator.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @param braCenters the number of centers on bra side.
     @param ketCenters the number of centers on ket side.
     @return the maximum order of recursion term object.
     */
    int32_t getMaxOrder(const std::string&  label,
                        const CFourIndexes& braAngularMomentum,
                        const CFourIndexes& ketAngularMomentum,
                        const int32_t       braCenters,
                        const int32_t       ketCenters) const;

    /**
     Finds recursion term object in recursion map object.

     @param recursionTerm the recursion term object.
     @return true if recursion term object is in recursion map objects,
             false otherwise.
     */
    bool find(const CRecursionTerm& recursionTerm) const;
    
    /**
     Allocates array of 2D memory blocks for storing recursion data.
     
     @param length the length of buffer.
     @return the array of 2D memory blocks.
     */
    CMemBlock2D<double>* createBuffer(const int32_t length) const;
    
    /**
     Deallocates array of 2D memory blocks for storing recursion data.

     @param buffer the array of 2D memory blocks. 
     */
    void destroyBuffer(CMemBlock2D<double>* buffer) const;

    /**
     Converts recursion map object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the recursion map object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CRecursionMap& source);
};

#endif /* RecursionMap_hpp */
