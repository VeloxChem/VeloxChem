//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef TwoIndexes_hpp
#define TwoIndexes_hpp

#include <cstdint>
#include <ostream>

/**
 Class CTwoIndexes stores information about pair of indexes and provides
 functions to manipulate these indexes.
 
 @author Z. Rinkevicius
 */
class CTwoIndexes
{
    /**
     The first index from pair of indexes.
     */
    int32_t _iIndex;

    /**
     The second index from pair of indexes.
     */
    int32_t _jIndex;

public:

    /**
     Creates an empty two indexes object.
     */
    CTwoIndexes();

    /**
     Creates a two indexes object.
     
     @param iIndex the first index from pair of indexes.
     @param jIndex the second index from pair of indexes.
     */
    CTwoIndexes(const int32_t iIndex,
                const int32_t jIndex);

    /**
     Creates a two indexes object by copying other two indexes object.
     
     @param source the two indexes object.
     */
    CTwoIndexes(const CTwoIndexes& source);

    /**
     Creates a two indexes object by moving other two indexes object.
     
     @param source the two indexes object.
     */
    CTwoIndexes(CTwoIndexes&& source) noexcept;

    /**
     Destroys a two indexes object.
     */
    ~CTwoIndexes();

    /**
     Assigns a two indexes object by copying other two indexes object.
     
     @param source the two indexes object.
     */
    CTwoIndexes& operator=(const CTwoIndexes& source);

    /**
     Assigns a two indexes object by moving other two indexes object.
     
     @param source the two indexes object.
     */
    CTwoIndexes& operator=(CTwoIndexes&& source) noexcept;

    /**
     Compares two indexes object with other two indexes object.
     
     @param other the two indexes object.
     @return true if two indexes objects are equal, false otherwise.
     */
    bool operator==(const CTwoIndexes& other) const;

    /**
     Compares two indexes object with other two indexes object.
     
     @param other the two indexes object.
     @return true if two indexes objects are not equal, false otherwise.
     */
    bool operator!=(const CTwoIndexes& other) const;

    /**
     Gets first index from pair of indexes.

     @return the first index from pait of indexes.
     */
    int32_t first() const;

    /**
     Gets second index from pair of indexes.
     
     @return the second index from pait of indexes.
     */
    int32_t second() const;
    
    /**
     Check if pair of indexes represents valid index pair (0..+Inf, +0..Inf).

     @return true if pair of indexes is valid index pair, false - otherwise.
     */
    bool isValidPair() const;
    
    /**
     Converts two indexes object to text output and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the two indexes object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CTwoIndexes&  source);
};

#endif /* TwoIndexes_hpp */
