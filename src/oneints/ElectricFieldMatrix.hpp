//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricFieldMatrix_hpp
#define ElectricFieldMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class CElectricFieldMatrix stores electric field matrix and provides
 set of methods for handling of electric field matrix data.
 
 @author Z. Rinkevicius
 */
class CElectricFieldMatrix
{
    /**
     The generic dense x component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _xMatrix;
    
    /**
     The generic dense Y component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _yMatrix;
    
    /**
     The generic dense Z component of electric field matrix (rectangular or
     square).
     */
    CDenseMatrix _zMatrix;
    
public:
    
    /**
     Creates an empty electric field  matrix object.
     */
    CElectricFieldMatrix();
    
    /**
     Creates a electric field matrix object.
     
     @param xMatrix the dense matrix with X component of electric field
     integrals.
     @param yMatrix the dense matrix with Y component of electric field
     integrals.
     @param zMatrix the dense matrix with Z component of electric field
     integrals.
     */
    CElectricFieldMatrix(const CDenseMatrix& xMatrix,
                         const CDenseMatrix& yMatrix,
                         const CDenseMatrix& zMatrix);
    
    /**
     Creates a electric field matrix object by copying other electric field
     matrix object.
     
     @param source the electric field matrix object.
     */
    CElectricFieldMatrix(const CElectricFieldMatrix& source);
    
    /**
     Creates a electric field matrix object by moving other electric field
     matrix object.
     
     @param source the electric field matrix object.
     */
    CElectricFieldMatrix(CElectricFieldMatrix&& source) noexcept;
    
    /**
     Destroys a electric field matrix object.
     */
    ~CElectricFieldMatrix();
    
    /**
     Assigns a electric field matrix object by copying other electric field
     matrix object.
     
     @param source the electric field matrix object.
     */
    CElectricFieldMatrix& operator=(const CElectricFieldMatrix& source);
    
    /**
     Assigns a electric field matrix object by moving other electric field
     matrix object.
     
     @param source the electric field matrix object.
     */
    CElectricFieldMatrix& operator=(CElectricFieldMatrix&& source) noexcept;
    
    /**
     Compares electric field matrix object with other electric field matrix
     object.
     
     @param other the electric field matrix object.
     @return true if electric field matrix objects are equal, false otherwise.
     */
    bool operator==(const CElectricFieldMatrix& other) const;
    
    /**
     Compares electric field matrix object with other electric field matrix
     object.
     
     @param other the electric field matrix object.
     @return true if electric field matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CElectricFieldMatrix& other) const;
    
    /**
     Gets string representation of X component of electric field matrix.
     
     @return a string for printing the X component of electric field matrix.
     */
    std::string getStringForComponentX() const;
    
    /**
     Gets string representation of Y component of electric field matrix.
     
     @return a string for printing the Y component of electric field matrix.
     */
    std::string getStringForComponentY() const;
    
    /**
     Gets string representation of Z component of electric field matrix.
     
     @return a string for printing the Z component of electric field matrix.
     */
    std::string getStringForComponentZ() const;
    
    /**
     Gets number of rows in electric field matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in electric field matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in electric field matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of X component of electric field
     matrix.
     
     @return the constant pointer to first element of X component of electric
     field matrix.
     */
    const double* xvalues() const;
    
    /**
     Gets constant pointer to first element of Y component of electric field
     matrix.
     
     @return the constant pointer to first element of Y component of electric
     field matrix.
     */
    const double* yvalues() const;
    
    /**
     Gets constant pointer to first element of Z component of electric field
     matrix.
     
     @return the constant pointer to first element of Z component of electric
     field matrix.
     */
    const double* zvalues() const;
    
    /**
     Converts electric field matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the electric field matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&         output,
                                    const CElectricFieldMatrix& source);
};

#endif /* ElectricFieldMatrix_hpp */
