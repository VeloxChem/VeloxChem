//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GtoContainer_hpp
#define GtoContainer_hpp

#include <vector>

/**
 Class CGtoBlock stores vector of GTOs block objects and provides set of methods
 for manipulating with basis functions of various angular momentum.
 
 @author Z. Rinkevicius
 */
#include "GtoBlock.hpp"

class CGtoContainer
{
    /**
     The maximum angular momentum.
     */
    int32_t _maxAngularMomentum;

    /**
     The vector of GTOs block objects.
     */
    std::vector<CGtoBlock> _gtoBlocks; 

public:

    /**
     Creates an empty GTOs container object.
     */
    CGtoContainer();

    /**
     Creates a GTOs container object.
     
     @param gtoBlocks the vector GTOs block objects.
     */
    CGtoContainer(const std::vector<CGtoBlock>& gtoBlocks);

    /**
     Creates a GTOs container object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     */
    CGtoContainer(const CMolecule&       molecule,
                  const CMolecularBasis& basis);

    /**
     Creates a GTOs container object by copying other GTOs container object.
     
     @param source the GTOs block object.
     */
    CGtoContainer(const CGtoContainer& source);

    /**
     Creates a GTOs container object by moving other GTOs container object.
     
     @param source the GTOs block object.
     */
    CGtoContainer(CGtoContainer&& source) noexcept;

    /**
     Destroys a GTOs container object.
     */
    ~CGtoContainer();

    /**
     Assigns a GTOs container object by copying other GTOs container object.
     
     @param source the GTOs block object.
     */
    CGtoContainer& operator=(const CGtoContainer& source);

    /**
     Assigns a GTOs container object by moving other GTOs container object.
     
     @param source the GTOs block object.
     */
    CGtoContainer& operator=(CGtoContainer&& source) noexcept;

    /**
     Compares GTOs container object with other GTOs container object.
     
     @param other the GTOs container object.
     @return true if GTOs container objects are equal, false otherwise.
     */
    bool operator==(const CGtoContainer& other) const;

    /**
     Compares GTOs container object with other GTOs container object.
     
     @param other the GTOs container object.
     @return true if GTOs container objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoContainer& other) const;
    
    /**
     Gets maximum angular momentum of GTOs block objects in GTOs container.

     @return the maxumum angular momentum.
     */
    int32_t getMaxAngularMomentum() const;
    
    /**
     Gets angular momentum of specific GTOs block object from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the angular momentum.
     */
    int32_t getAngularMomentum(const int32_t iBlock) const;
    
    /**
     Gets number of GTOs blocks in GTOs container.

     @return the number of GTOs blocks.
     */
    int32_t getNumberOfGtoBlocks() const;
    
    /**
     Gets maximum number of primitive Gaussian functions within GTOs container.

     @return the number of primitive Gaussian functions.
     */
    int32_t getMaxNumberOfPrimGtos() const;
    
    /**
     Gets maximum number of contracted basis functions within GTOs container.

     @return the number of contracted GTOs.
     */
    int32_t getMaxNumberOfContrGtos() const;
    
    /**
     Gets number of primitive Gaussian functions in specific GTOs block object
     from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimGtos(const int32_t iBlock) const;
    
    /**
     Gets number of contracted basis functions in specific GTOs block object
     from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the number of contracted basis functions.
     */
    int32_t getNumberOfContrGtos(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to basis function start positions in primitive
     Gaussian functions vector from specific GTOs block object in GTOs
     container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the start positions of basis fucntions.
     */
    const int32_t* getStartPositions(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to basis function end positions in primitive
     Gaussian functions vector from specific GTOs block object in GTOs
     container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the end positions of basis fucntions.
     */
    const int32_t* getEndPositions(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to basis function indexes in full AO basis for
     specific angular momentum component from specific GTOs block object in
     GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    const int32_t* getIdentifiers(const int32_t iBlock,
                                  const int32_t iComponent) const;
    
    /**
     Gets constant pointer to exponents of primitive Gaussian functions from
     specific GTOs block object in GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getExponents(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to normalization factors of primitive Gaussian
     functions from specific GTOs block object in GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the normalization factors of primitive Gaussian functions.
     */
    const double* getNormFactors(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to Cartesian X coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesX(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to Cartesian Y coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesY(const int32_t iBlock) const;
    
    /**
     Gets constant pointer to Cartesian Z coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.
     
     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesZ(const int32_t iBlock) const;
    
    /**
     Converts GTOs container object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the GTOs block object.
     */
    friend std::ostream& operator<<(      std::ostream&  output,
                                    const CGtoContainer& source);
};

#endif /* GtoContainer_hpp */
