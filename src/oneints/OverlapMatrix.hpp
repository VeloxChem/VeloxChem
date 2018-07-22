//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OverlapMatrix_hpp
#define OverlapMatrix_hpp

#include "SparseMatrix.hpp"

/**
 Class COverlapMatrix stores general one electron overlap matrix and provides
 set of methods for handling of overlap matrix data.
 
 @author Z. Rinkevicius
 */
class COverlapMatrix
{
    /**
     The generic sparse overlap matrix (rectangular or square).
     */
    CSparseMatrix _matrix;
    
public:
    
    /**
     Creates an empty overlap matrix object.
     */
    COverlapMatrix();
    
    /**
     Creates a overlap matrix object.
     
     @param matrix the sparse matrix with overlap integrals.
     */
    COverlapMatrix(const CSparseMatrix& matrix);
    
    /**
     Creates a overlap matrix object by copying other overlap matrix object.
     
     @param source the overlap matrix object.
     */
    COverlapMatrix(const COverlapMatrix& source);
    
    /**
     Creates a overlap matrix object by moving other overlap matrix object.
     
     @param source the overlap matrix object.
     */
    COverlapMatrix(COverlapMatrix&& source) noexcept;
    
    /**
     Destroys a overlap matrix object.
     */
    ~COverlapMatrix();
    
    /**
     Assigns a overlap matrix object by copying other overlap matrix object.
     
     @param source the overlap matrix object.
     */
    COverlapMatrix& operator=(const COverlapMatrix& source);
    
    /**
     Assigns a overlap matrix object by moving other overlap matrix object.
     
     @param source the basis function object.
     */
    COverlapMatrix& operator=(COverlapMatrix&& source) noexcept;
    
    /**
     Compares overlap matrix object with other overlap matrix object.
     
     @param other the overlap matrix object.
     @return true if overlap matrix objects are equal, false otherwise.
     */
    bool operator==(const COverlapMatrix& other) const;
    
    /**
     Compares overlap matrix object with other overlap matrix object.
     
     @param other the overlap matrix object.
     @return true if overlap matrix objects are not equal, false otherwise.
     */
    bool operator!=(const COverlapMatrix& other) const;
    
    /**
     Converts overlap matrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the overlap matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const COverlapMatrix& source);
};



#endif /* OverlapMatrix_hpp */
