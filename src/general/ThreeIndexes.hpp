//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ThreeIndexes_hpp
#define ThreeIndexes_hpp

#include <cstdint>
#include <ostream>

/**
 Class CThreeIndexes stores information about triple of indexes and provides
 functions to manipulate these indexes.

 @author Z. Rinkevicius
 */
class CThreeIndexes
{
    /**
     The first index from triple of indexes.
     */
    int32_t _iIndex;

    /**
     The second index from triple of indexes.
     */
    int32_t _jIndex;

    /**
     The third index from triple of indexes.
     */
    int32_t _kIndex;

   public:
    /**
     Creates an empty three indexes object.
     */
    CThreeIndexes();

    /**
     Creates a three indexes object.

     @param iIndex the first index from triple of indexes.
     @param jIndex the second index from triple of indexes.
     @param kIndex the third index from triple of indexes.
     */
    CThreeIndexes(const int32_t iIndex, const int32_t jIndex, const int32_t kIndex);

    /**
     Creates a three indexes object by copying other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes(const CThreeIndexes& source);

    /**
     Creates a three indexes object by moving other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes(CThreeIndexes&& source) noexcept;

    /**
     Destroys a three indexes object.
     */
    ~CThreeIndexes();

    /**
     Assigns a three indexes object by copying other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes& operator=(const CThreeIndexes& source);

    /**
     Assigns a three indexes object by moving other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes& operator=(CThreeIndexes&& source) noexcept;

    /**
     Compares three indexes object with other three indexes object.

     @param other the three indexes object.
     @return true if three indexes objects are equal, false otherwise.
     */
    bool operator==(const CThreeIndexes& other) const;

    /**
     Compares three indexes object with other three indexes object.

     @param other the three indexes object.
     @return true if three indexes objects are not equal, false otherwise.
     */
    bool operator!=(const CThreeIndexes& other) const;

    /**
     Gets first index from triple of indexes.

     @return the first index from triple of indexes.
     */
    int32_t first() const;

    /**
     Gets second index from triple of indexes.

     @return the second index from triple of indexes.
     */
    int32_t second() const;

    /**
     Gets third index from triple of indexes.

     @return the third index from triple of indexes.
     */
    int32_t third() const;

    /**
     Check if triple of indexes represents valid index triple (0..+Inf, +0..Inf,
     +0..Inf).

     @return true if triple of indexes is valid index triple, false - otherwise.
     */
    bool isValidTriple() const;

    /**
     Converts triple indexes object to text output and insert it into output text
     stream.

     @param output the output text stream.
     @param source the triple indexes object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CThreeIndexes& source);
};

#endif /* ThreeIndexes_hpp */
