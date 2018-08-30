//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GtoPairsContainer_hpp
#define GtoPairsContainer_hpp

#include <vector>

#include "GtoPairsBlock.hpp"
#include "OutputStream.hpp"

/**
 Class CGtoPairsContainer stores vector of GTOs pairs block objects and provides
 set of methods for manipulating with basis function pairs of various angular
 momentum.
 
 @author Z. Rinkevicius
 */
class CGtoPairsContainer
{
    /**
     The vector of GTOs pairs block objects.
     */
    std::vector<CGtoPairsBlock> _gtoPairsBlocks;
    
public:
    
    /**
     Creates an empty GTOs pairs container object.
     */
    CGtoPairsContainer();
    
    /**
     Creates a GTOs pairs container object.
     
     @param gtoPairsBlocks the vector GTOs pairs block objects.
     */
    CGtoPairsContainer(const std::vector<CGtoPairsBlock>& gtoPairsBlocks);
    
    /**
     Creates a GTOs pairs container object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param threshold the primitive pairs cut-off threshold.
     */
    CGtoPairsContainer(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold);
    
    /**
     Creates a GTOs pairs container object by copying other GTOs pairs container
     object.
     
     @param source the GTOs pairs container object.
     */
    CGtoPairsContainer(const CGtoPairsContainer& source);
    
    /**
     Creates a GTOs pairs container object by moving other GTOs pairs container
     object.
     
     @param source the GTOs pairs container object.
     */
    CGtoPairsContainer(CGtoPairsContainer&& source) noexcept;
    
    /**
     Destroys a GTOs pairs container object.
     */
    ~CGtoPairsContainer();
    
    /**
     Assigns a GTOs pairs container object by copying other GTOs pairs container
     object.
     
     @param source the GTOs pairs container object.
     */
    CGtoPairsContainer& operator=(const CGtoPairsContainer& source);
    
    /**
     Assigns a GTOs pairs container object by moving other GTOs pairs container
     object.
     
     @param source the GTOs pairs container object.
     */
    CGtoPairsContainer& operator=(CGtoPairsContainer&& source) noexcept;
    
    /**
     Compares GTOs pairs container object with other GTOs pairs container object.
     
     @param other the GTOs pairs container object.
     @return true if GTOs pairs container objects are equal, false otherwise.
     */
    bool operator==(const CGtoPairsContainer& other) const;
    
    /**
     Compares GTOs pairs container object with other GTOs pairs container object.
     
     @param other the GTOs pairs container object.
     @return true if GTOs pairs container objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoPairsContainer& other) const;
    
    
    /**
     Creates a GTOs pairs container object by splitting the GTOs pairs container
     object.

     @param batchSize the size of batch.
     @return the GTOs pairs container object.
     */
    CGtoPairsContainer split(const int32_t batchSize) const;
    
    /**
     Gets numnber of GTOs pairs block objects in GTOs pairs container.

     @return the number of GTOs pairs block objects.
     */
    int32_t getNumberOfGtoPairsBlocks() const;
    
    /**
     Gets specific GTOs pairs block object from GTOs pairs container.

     @param iBlock the index of requested GTOs pairs block object.
     @return the GTOs pairs block object.
     */
    CGtoPairsBlock getGtoPairsBlock(const int32_t iBlock) const;
    
    /**
     Prints GTO pairs screening information to output stream.

     @param oStream the output stream.
     */
    void printScreeningInfo(COutputStream& oStream) const;
    
    /**
     Converts GTOs pairs container object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the GTOs pairs container object.
     */
    friend std::ostream& operator<<(      std::ostream&       output,
                                    const CGtoPairsContainer& source);
};

#endif /* GtoPairsContainer_hpp */
