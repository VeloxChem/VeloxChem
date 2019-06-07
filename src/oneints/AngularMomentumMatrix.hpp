//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AngularMomentumMatrix_hpp
#define AngularMomentumMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class CAngularMomentumMatrix stores angular momentum matrix and provides
 set of methods for handling of angular momentum matrix data.
 
 @author Z. Rinkevicius
 */
class CAngularMomentumMatrix
{
    /**
     The X coordinate of angular momentum origin.
     */
    double _xOrigin;
    
    /**
     The Y coordinate of angular momentum origin.
     */
    double _yOrigin;
    
    /**
     The Z coordinate of angular momentum origin.
     */
    double _zOrigin;
    
    /**
     The generic dense x component of angular momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _xMatrix;
    
    /**
     The generic dense Y component of angular momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _yMatrix;
    
    /**
     The generic dense Z component of angular momentum matrix (rectangular or
     square).
     */
    CDenseMatrix _zMatrix;
    
public:
    
    /**
     Creates an empty angular momentum  matrix object.
     */
    CAngularMomentumMatrix();
    
    /**
     Creates a angular momentum matrix object.
     
     @param xMatrix the dense matrix with X component of angular momentum
            integrals.
     @param yMatrix the dense matrix with Y component of angular momentum
            integrals.
     @param zMatrix the dense matrix with Z component of angular momentum
            integrals.
     @param xOrigin the Cartesian X coordinate of angular momentum origin.
     @param yOrigin the Cartesian Y coordinate of angular momentum origin.
     @param zOrigin the Cartesian Z coordinate of angular momentum origin.
     */
    CAngularMomentumMatrix(const CDenseMatrix& xMatrix,
                          const CDenseMatrix& yMatrix,
                          const CDenseMatrix& zMatrix,
                          const double        xOrigin,
                          const double        yOrigin,
                          const double        zOrigin);
    
    /**
     Creates a angular momentum matrix object by copying other angular momentum
     matrix object.
     
     @param source the angular momentum matrix object.
     */
    CAngularMomentumMatrix(const CAngularMomentumMatrix& source);
    
    /**
     Creates a angular momentum  matrix object by moving other angular momentum
     matrix object.
     
     @param source the angular momentum matrix object.
     */
    CAngularMomentumMatrix(CAngularMomentumMatrix&& source) noexcept;
    
    /**
     Destroys a angular momentum matrix object.
     */
    ~CAngularMomentumMatrix();
    
    /**
     Assigns a angular momentum  matrix object by copying other angular momentum
     matrix object.
     
     @param source the angular momentum matrix object.
     */
    CAngularMomentumMatrix& operator=(const CAngularMomentumMatrix& source);
    
    /**
     Assigns a angular momentum matrix object by moving other angular momentum
     matrix object.
     
     @param source the angular momentum matrix object.
     */
    CAngularMomentumMatrix& operator=(CAngularMomentumMatrix&& source) noexcept;
    
    /**
     Compares angular momentum matrix object with other angular momentum matrix
     object.
     
     @param other the angular momentum matrix object.
     @return true if angular momentum matrix objects are equal, false otherwise.
     */
    bool operator==(const CAngularMomentumMatrix& other) const;
    
    /**
     Compares angular momentum  matrix object with other angular momentum matrix
     object.
     
     @param other the angular momentum  matrix object.
     @return true if angular momentum  matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CAngularMomentumMatrix& other) const;
    
    /**
     Sets coordinates of angular momentum origin.

     @param xOrigin the X coordinate of angular momentum origin.
     @param yOrigin the Y coordinate of angular momentum origin.
     @param zOrigin the Z coordinate of angular momentum origin.
     */
    void setOriginCoordinates(const double xOrigin,
                              const double yOrigin,
                              const double zOrigin);
    
    /**
     Gets string representation of X component of angular momentum matrix.
     
     @return a string for printing the X component of angular momentum matrix.
     */
    std::string getStringForComponentX() const;
    
    /**
     Gets string representation of Y component of angular momentum matrix.
     
     @return a string for printing the Y component of angular momentum matrix.
     */
    std::string getStringForComponentY() const;
    
    /**
     Gets string representation of Z component of angular momentum matrix.
     
     @return a string for printing the Z component of angular momentum matrix.
     */
    std::string getStringForComponentZ() const;
    
    /**
     Gets number of rows in angular momentum matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in angular momentum matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in angular momentum matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of X component of angular momentum
     matrix.
     
     @return the constant pointer to first element of X component of electric
     dipole matrix.
     */
    const double* xvalues() const;
    
    /**
     Gets constant pointer to first element of Y component of angular momentum
     matrix.
     
     @return the constant pointer to first element of Y component of electric
     dipole matrix.
     */
    const double* yvalues() const;
    
    /**
     Gets constant pointer to first element of Z component of angular momentum
     matrix.
     
     @return the constant pointer to first element of Z component of electric
     dipole matrix.
     */
    const double* zvalues() const;
    
    /**
     Gets X coordinate of angular momentum origin.

     @return the X coordinate of angular momentum origin.
     */
    double getOriginCoordinateX() const;
    
    /**
     Gets Y coordinate of angular momentum origin.
     
     @return the Y coordinate of angular momentum origin.
     */
    double getOriginCoordinateY() const;
    
    /**
     Gets Z coordinate of angular momentum origin.
     
     @return the Z coordinate of angular momentum origin.
     */
    double getOriginCoordinateZ() const;
    
    /**
     Converts angular momentum matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the angular momentum matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&          output,
                                    const CAngularMomentumMatrix& source);
};


#endif /* AngularMomentumMatrix_hpp */
