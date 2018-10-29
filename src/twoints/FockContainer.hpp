//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef FockContainer_hpp
#define FockContainer_hpp

#include <cstdint>
#include <vector>

#include "FockSubMatrix.hpp"
#include "AOFockMatrix.hpp"

/**
 Class CFockContainer stores vector of Fock submatrices and provides set of
 methods for handling of Fock submatrices data.
 
 @author Z. Rinkevicius
 */
class CFockContainer
{
    /**
     The vector of Fock submatrices.
     */
    std::vector<CFockSubMatrix> _subFockMatrices;
    
public:
    
    /**
     Creates an empty Fock container object.
     */
    CFockContainer();
    
    /**
     Creates a Fock container object.
     
     @param subFockMatrices the vector of Fock submatrices.
     */
    CFockContainer(const std::vector<CFockSubMatrix>& subFockMatrices);
    
    /**
     Creates a Fock container object.
     
     @param fockMatrix the pointer to AO Fock matrix.
     @param braGtoPairsBlock the GTOs pairsblock on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    CFockContainer(const CAOFockMatrix*  fockMatrix,
                   const CGtoPairsBlock& braGtoPairsBlock,
                   const CGtoPairsBlock& ketGtoPairsBlock);
    
    /**
     Creates a Fock container object by copying other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer(const CFockContainer& source);
    
    /**
     Creates a Fock container object by moving other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer(CFockContainer&& source) noexcept;
    
    /**
     Destroys a Fock container object.
     */
    ~CFockContainer();
    
    /**
     Assigns a Fock container object by copying other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer& operator=(const CFockContainer& source);
    
    /**
     Assigns a Fock container object by moving other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer& operator=(CFockContainer&& source) noexcept;
    
    /**
     Compares Fock container object with other Fock container object.
     
     @param other the Fock container object.
     @return true if Fock container objects are equal, false otherwise.
     */
    bool operator==(const CFockContainer& other) const;
    
    /**
     Compares Fock container object with other Fock container object.
     
     @param other the Fock container object.
     @return true if Fock container objects are not equal, false otherwise.
     */
    bool operator!=(const CFockContainer& other) const;
    
    /**
     Converts Fock container object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the Fock container object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CFockContainer& source);
};

#endif /* FockContainer_hpp */
