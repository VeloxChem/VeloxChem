//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GtoBlock_hpp
#define GtoBlock_hpp

#include "MemBlock2D.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

/**
 Class CGtoBlock stores data about basis functions and provides set of methods
 for manipulating with basis functions.
 
 @author Z. Rinkevicius
 */
class CGtoBlock
{
    /**
     The angular momentum.
     */
    int32_t _angularMomentum;

    /**
     The contraction scheme (contracted basis functions start/end position in
     primitive Gaussian functions vector, contracted basis functions indexes).
     */
    CMemBlock2D<int32_t> _contrPattern;

    /**
     The primitives Gaussian functions data (exponents, norm. factors and
     coordinates).
     */
    CMemBlock2D<double> _gtoPrimitives; 

public:

    /**
     Creates an empty GTOs block object.
     */
    CGtoBlock();

    /**
     Creates a GTOs block object.
     
     @param gtoPrimitives the primitives Gaussian functions data (exponents,
                           norm. factors, coordinates).
     @param contrPattern the contraction pattern (contracted basis functions
                         start/end position in primitive Gaussian functions
                         vector, basis functions indexes).
     @param angularMomentum the angular momentum of contracted basis functions.
     */
    CGtoBlock(const CMemBlock2D<double>&  gtoPrimitives,
              const CMemBlock2D<int32_t>& contrPattern,
              const int32_t               angularMomentum);

    /**
     Creates a GTOs block object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param angularMomentum the angular momentum of contracted basis functions.
     */
    CGtoBlock(const CMolecule&       molecule,
              const CMolecularBasis& basis,
              const int32_t          angularMomentum);

    /**
     Creates a GTOs block object by copying other GTOs block object.
     
     @param source the GTOs block object.
     */
    CGtoBlock(const CGtoBlock& source);

    /**
     Creates a GTOs block object by moving other GTOs block object.
     
     @param source the GTOs block object.
     */
    CGtoBlock(CGtoBlock&& source) noexcept;

    /**
     Destroys a GTOs block object.
     */
    ~CGtoBlock();

    /**
     Assigns a GTOs block object by copying other GTOs block object.
     
     @param source the GTOs block object.
     */
    CGtoBlock& operator=(const CGtoBlock& source);

    /**
     Assigns a GTOs block object by moving other GTOs block object.
     
     @param source the GTOs block object.
     */
    CGtoBlock& operator=(CGtoBlock&& source) noexcept;

    /**
     Compares GTOs block object with other GTOs block object.
     
     @param other the GTOs block object.
     @return true if GTOs block objects are equal, false otherwise.
     */
    bool operator==(const CGtoBlock& other) const;

    /**
     Compares GTOs block object with other GTOs block object.
     
     @param other the GTOs block object.
     @return true if GTOs block objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoBlock& other) const;
    
    /**
     Gets angular momentum of GTOs block.

     @return the angular momentum.
     */
    int32_t getAngularMomentum() const;
    
    /**
     Checks if GTOs block is empty.

     @return true if GTOs block is empty, false otherwise.
     */
    bool empty() const;
    
    /**
     Gets number of primitive Gaussian functions in GTOs block.

     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimGtos() const;
    
    /**
     Gets number of contracted basis functions in GTOs block.

     @return the number of contracted basis functions.
     */
    int32_t getNumberOfContrGtos() const;
    
    /**
     Gets constant pointer to basis function start positions in primitive
     Gaussian functions vector.

     @return the start positions of basis fucntions.
     */
    const int32_t* getStartPositions() const;
    
    /**
     Gets constant pointer to basis function end positions in primitive
     Gaussian functions vector.

     @return the end positions of basis fucntions.
     */
    const int32_t* getEndPositions() const;
    
    /**
     Gets constant pointer to basis function indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    const int32_t* getIdentifiers(const int32_t iComponent) const; 
    
    /**
     Gets constant pointer to exponents of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    const double* getExponents() const;
    
    /**
     Gets constant pointer to normalization factors of primitive Gaussian
     functions.
     
     @return the normalization factors of primitive Gaussian functions.
     */
    const double* getNormFactors() const;
    
    /**
     Gets constant pointer to Cartesian X coordinates of primitive Gaussian
     functions.
     
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesX() const;
    
    /**
     Gets constant pointer to Cartesian Y coordinates of primitive Gaussian
     functions.
     
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesY() const;
    
    /**
     Gets constant pointer to Cartesian Z coordinates of primitive Gaussian
     functions.
     
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesZ() const;
    
    /**
     Converts GTOs block object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the GTOs block object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CGtoBlock&    source);
};

#endif /* GtoBlock_hpp */
