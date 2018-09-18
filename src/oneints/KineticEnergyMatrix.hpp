//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef KineticEnergyMatrix_hpp
#define KineticEnergyMatrix_hpp

#include "DenseMatrix.hpp"

/**
 Class CKineticEnergyMatrix stores general kinetic energy matrix and provides
 set of methods for handling of kinetic energy matrix data.
 
 @author Z. Rinkevicius
 */
class CKineticEnergyMatrix
{
    /**
     The generic dense kinetic energy matrix (rectangular or square).
     */
    CDenseMatrix _matrix;
    
public:
    
    /**
     Creates an empty kinetic energy matrix object.
     */
    CKineticEnergyMatrix();
    
    /**
     Creates a kinetic energy matrix object.
     
     @param matrix the dense matrix with kinetic energy integrals.
     */
    CKineticEnergyMatrix(const CDenseMatrix& matrix);
    
    /**
     Creates a kinetic energy matrix object by copying other kinetic energy
     matrix object.
     
     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix(const CKineticEnergyMatrix& source);
    
    /**
     Creates a kinetic energy matrix object by moving other kinetic energy
     matrix object.
     
     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix(CKineticEnergyMatrix&& source) noexcept;
    
    /**
     Destroys a kinetic energy matrix object.
     */
    ~CKineticEnergyMatrix();
    
    /**
     Assigns a kinetic energy matrix object by copying other kinetic energy
     matrix object.
     
     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix& operator=(const CKineticEnergyMatrix& source);
    
    /**
     Assigns a kinetic energy matrix object by moving other kinetic energy
     matrix object.
     
     @param source the kinetic energy matrix object.
     */
    CKineticEnergyMatrix& operator=(CKineticEnergyMatrix&& source) noexcept;
    
    /**
     Compares kinetic energy matrix object with other kinetic energy matrix
     object.
     
     @param other the kinetic energy matrix object.
     @return true if kinetic energy matrix objects are equal, false otherwise.
     */
    bool operator==(const CKineticEnergyMatrix& other) const;
    
    /**
     Compares kinetic energy matrix object with other kinetic energy matrix
     object.
     
     @param other the kinetic energy matrix object.
     @return true if kinetic energy matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CKineticEnergyMatrix& other) const;
    
    /**
     Gets string representation of kinetic energy matrix.
     
     @return a string for printing the kinetic energy matrix.
     */
    std::string getString() const;
    
    /**
     Gets number of rows in kinetic energy matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in kinetic energy matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in kinetic energy matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of kinetic energy matrix.
     
     @return the constant pointer to first element of kinetic energy matrix.
     */
    const double* values() const;
    
    /**
     Converts kinetic energy matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the kinetic energy matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&         output,
                                    const CKineticEnergyMatrix& source);
};


#endif /* KineticEnergyMatrix_hpp */
