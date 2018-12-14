//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef OverlapMatrix_hpp
#define OverlapMatrix_hpp

#include <string>

#include "DenseMatrix.hpp"
#include "SystemClock.hpp"

/**
 Class COverlapMatrix stores general one electron overlap matrix and provides
 set of methods for handling of overlap matrix data.
 
 @author Z. Rinkevicius
 */
class COverlapMatrix
{
    /**
     The generic dense overlap matrix (rectangular or square).
     */
    CDenseMatrix _matrix;
    
public:
    
    /**
     Creates an empty overlap matrix object.
     */
    COverlapMatrix();
    
    /**
     Creates a overlap matrix object.
     
     @param matrix the dense matrix with overlap integrals.
     */
    COverlapMatrix(const CDenseMatrix& matrix);
    
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
     
     @param source the overlap matrix object.
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
     Gets string representation of overlap matrix.
     
     @return a string for printing the overlap matrix.
     */
    std::string getString() const;
    
    /**
     Gets number of rows in overlap matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in overlap matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in overlap matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of overlap matrix.
     
     @return the constant pointer to first element of overlap matrix.
     */
    const double* values() const;
    
    /**
     Gets Löwdin symmetric orthogonalization matrix if overalp matrix
     eigenvalues is linearly independent, or canonical orthogonalization
     matrix if overlap matrix eigenvalues is linearly dependent.

     @param threshold the linear dependence threshold.
     @return the orthogonalization matrix.
     */
    CDenseMatrix getOrthogonalizationMatrix(const double threshold) const;
    
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
