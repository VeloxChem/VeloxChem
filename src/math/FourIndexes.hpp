//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef FourIndexes_hpp
#define FourIndexes_hpp

#include <cstdint>
#include <ostream>

/**
 Class CFourIndexes stores information about triple of indexes and provides
 functions to manipulate these indexes.
 
 @author Z. Rinkevicius
 */
class CFourIndexes
{
    /**
     The first index from quadruple of indexes.
     */
    int32_t _iIndex;
    
    /**
     The second index from quadruple of indexes.
     */
    int32_t _jIndex;
    
    /**
     The third index from quadruple of indexes.
     */
    int32_t _kIndex;
    
    /**
     The fourth index from quadruple of indexes.
     */
    int32_t _lIndex;
    
public:
    
    /**
     Creates an empty four indexes object.
     */
    CFourIndexes();
    
    /**
     Creates a four indexes object.
     
     @param iIndex the first index from quadruple of indexes.
     @param jIndex the second index from quadruple of indexes.
     @param kIndex the third index from quadruple of indexes.
     @param lIndex the fourth index from quadruple of indexes.
     */
    CFourIndexes(const int32_t iIndex,
                 const int32_t jIndex,
                 const int32_t kIndex,
                 const int32_t lIndex);
    
    /**
     Creates a four indexes object by copying other four indexes object.
     
     @param source the four indexes object.
     */
    CFourIndexes(const CFourIndexes& source);
    
    /**
     Creates a four indexes object by moving other four indexes object.
     
     @param source the four indexes object.
     */
    CFourIndexes(CFourIndexes&& source) noexcept;
    
    /**
     Destroys a four indexes object.
     */
    ~CFourIndexes();
    
    /**
     Assigns a four indexes object by copying other four indexes object.
     
     @param source the four indexes object.
     */
    CFourIndexes& operator=(const CFourIndexes& source);
    
    /**
     Assigns a four indexes object by moving other four indexes object.
     
     @param source the four indexes object.
     */
    CFourIndexes& operator=(CFourIndexes&& source) noexcept;
    
    /**
     Compares four indexes object with other four indexes object.
     
     @param other the four indexes object.
     @return true if four indexes objects are equal, false otherwise.
     */
    bool operator==(const CFourIndexes& other) const;
    
    /**
     Compares four indexes object with other four indexes object.
     
     @param other the Four indexes object.
     @return true if four indexes objects are not equal, false otherwise.
     */
    bool operator!=(const CFourIndexes& other) const;
    
    /**
     Gets first index from quadruple of indexes.
     
     @return the first index from quadruple of indexes.
     */
    int32_t first() const;
    
    /**
     Gets second index from quadruple of indexes.
     
     @return the second index from quadruple of indexes.
     */
    int32_t second() const;
    
    /**
     Gets third index from quadruple of indexes.
     
     @return the third index from quadruple of indexes.
     */
    int32_t third() const;
    
    /**
     Gets fourth index from quadruple of indexes.
     
     @return the fourth index from quadruple of indexes.
     */
    int32_t fourth() const;
    
    /**
     Check if quadruple of indexes represents valid index quadruple (0..+Inf,
     +0..Inf, +0..Inf, +0..Inf).
     
     @return true if quadruple  of indexes is valid index quadruple,
             false - otherwise.
     */
    bool isValidQuadruple() const;
    
    /**
     Converts four indexes object to text output and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the triple indexes object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CFourIndexes&  source);
};


#endif /* FourIndexes_hpp */
