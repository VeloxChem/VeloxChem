//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef NuclearPotentialMatrix_hpp
#define NuclearPotentialMatrix_hpp

#include "SparseMatrix.hpp"

/**
 Class CNuclearPotentialMatrix stores general nuclear potential matrix and
 provides set of methods for handling of nuclear potential matrix data.
 
 @author Z. Rinkevicius
 */
class CNuclearPotentialMatrix
{
    /**
     The generic sparse nuclear potential matrix (rectangular or square).
     */
    CSparseMatrix _matrix;
    
public:
    
    /**
     Creates an empty nuclear potential matrix object.
     */
    CNuclearPotentialMatrix();
    
    /**
     Creates a nuclear potential matrix object.
     
     @param matrix the sparse matrix with nuclear potential integrals.
     */
    CNuclearPotentialMatrix(const CSparseMatrix& matrix);
    
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
     Converts nuclear potential matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the nuclear potential matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&            output,
                                    const CNuclearPotentialMatrix& source);
};

#endif /* NuclearPotentialMatrix_hpp */
