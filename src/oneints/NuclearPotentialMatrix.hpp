//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef NuclearPotentialMatrix_hpp
#define NuclearPotentialMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class CNuclearPotentialMatrix stores general nuclear potential matrix and
 provides set of methods for handling of nuclear potential matrix data.
 
 @author Z. Rinkevicius
 */
class CNuclearPotentialMatrix
{
    /**
     The generic dense nuclear potential matrix (rectangular or square).
     */
    CDenseMatrix _matrix;
    
public:
    
    /**
     Creates an empty nuclear potential matrix object.
     */
    CNuclearPotentialMatrix();
    
    /**
     Creates a nuclear potential matrix object.
     
     @param matrix the dense matrix with nuclear potential integrals.
     */
    CNuclearPotentialMatrix(const CDenseMatrix& matrix);
    
    /**
     Creates a nuclear potential matrix object by copying other nuclear potential
     matrix object.
     
     @param source the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix(const CNuclearPotentialMatrix& source);
    
    /**
     Creates a nuclear potential matrix object by moving other nuclear potential
     matrix object.
     
     @param source the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix(CNuclearPotentialMatrix&& source) noexcept;
    
    /**
     Destroys a nuclear potential matrix object.
     */
    ~CNuclearPotentialMatrix();
    
    /**
     Assigns a nuclear potential matrix object by copying other nuclear potential
     matrix object.
     
     @param source the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix& operator=(const CNuclearPotentialMatrix& source);
    
    /**
     Assigns a nuclear potential matrix object by moving other nuclear potential
     matrix object.
     
     @param source the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix& operator=(CNuclearPotentialMatrix&& source) noexcept;
    
    /**
     Compares nuclear potential matrix object with other nuclear potential
     matrix object.
     
     @param other the nuclear potential matrix object.
     @return true if nuclear potential matrix objects are equal, false otherwise.
     */
    bool operator==(const CNuclearPotentialMatrix& other) const;
    
    /**
     Compares nuclear potential matrix object with other nuclear potential
     matrix object.
     
     @param other the nuclear potential matrix object.
     @return true if nuclear potential matrix objects are not equal, false
             otherwise.
     */
    bool operator!=(const CNuclearPotentialMatrix& other) const;
    
    /**
     Gets string representation of nuclear potential matrix.
     
     @return a string for printing the nuclear potential matrix.
     */
    std::string getString() const;
    
    /**
     Gets number of rows in nuclear potential matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in nuclear potential matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in nuclear potential matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of nuclear potential matrix.
     
     @return the constant pointer to first element of nuclear potential matrix.
     */
    const double* values() const;
    
    /**
     Converts nuclear potential matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the nuclear potential matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&            output,
                                    const CNuclearPotentialMatrix& source);
};

#endif /* NuclearPotentialMatrix_hpp */
