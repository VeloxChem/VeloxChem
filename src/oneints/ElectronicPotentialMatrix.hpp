//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectronicPotentialMatrix_hpp
#define ElectronicPotentialMatrix_hpp

#include <cstdint>
#include <string>

#include "DenseMatrix.hpp"

/**
 Class CElectronicPotenialMatrix stores general electronic potential matrix and
 provides set of methods for handling of electronic potential matrix data.
 
 @author Z. Rinkevicius
 */
class CElectronicPotentialMatrix
{
    /**
     The generic dense electronic potential matrix (rectangular or square).
     */
    CDenseMatrix _matrix;
    
public:
    
    /**
     Creates an empty electronic potential matrix object.
     */
    CElectronicPotentialMatrix();
    
    /**
     Creates a electronic potential matrix object.
     
     @param matrix the dense matrix with electronic potential integrals.
     */
    CElectronicPotentialMatrix(const CDenseMatrix& matrix);
    
    /**
     Creates a electronic potential matrix object by copying other electronic potential
     matrix object.
     
     @param source the electronic potential matrix object.
     */
    CElectronicPotentialMatrix(const CElectronicPotentialMatrix& source);
    
    /**
     Creates a electronic potential matrix object by moving other electronic potential
     matrix object.
     
     @param source the electronic potential matrix object.
     */
    CElectronicPotentialMatrix(CElectronicPotentialMatrix&& source) noexcept;
    
    /**
     Destroys a electronic potential matrix object.
     */
    ~CElectronicPotentialMatrix();
    
    /**
     Assigns a electronic potential matrix object by copying other electronic potential
     matrix object.
     
     @param source the electronic potential matrix object.
     */
    CElectronicPotentialMatrix& operator=(const CElectronicPotentialMatrix& source);
    
    /**
     Assigns a electronic potential matrix object by moving other electronic potential
     matrix object.
     
     @param source the electronic potential matrix object.
     */
    CElectronicPotentialMatrix& operator=(CElectronicPotentialMatrix&& source) noexcept;
    
    /**
     Compares electronic potential matrix object with other electronic potential
     matrix object.
     
     @param other the electronic potential matrix object.
     @return true if electronic potential matrix objects are equal, false otherwise.
     */
    bool operator==(const CElectronicPotentialMatrix& other) const;
    
    /**
     Compares electronic potential matrix object with other electronic potential
     matrix object.
     
     @param other the electronic potential matrix object.
     @return true if electronic potential matrix objects are not equal, false
     otherwise.
     */
    bool operator!=(const CElectronicPotentialMatrix& other) const;
    
    /**
     Converts electronic potential matrix object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the electronic potential matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&              output,
                                    const CElectronicPotentialMatrix& source);
};

#endif /* ElectronicPotentialMatrix_hpp */
