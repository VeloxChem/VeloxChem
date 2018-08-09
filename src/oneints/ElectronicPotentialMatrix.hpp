//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ElectronicPotentialMatrix_hpp
#define ElectronicPotentialMatrix_hpp

#include "SparseMatrix.hpp"

/**
 Class CElectronicPotenialMatrix stores general electronic potential matrix and
 provides set of methods for handling of electronic potential matrix data.
 
 @author Z. Rinkevicius
 */
class CElectronicPotentialMatrix
{
    /**
     The generic sparse electronic potential matrix (rectangular or square).
     */
    CSparseMatrix _matrix;
    
public:
    
    /**
     Creates an empty electronic potential matrix object.
     */
    CElectronicPotentialMatrix();
    
    /**
     Creates a electronic potential matrix object.
     
     @param matrix the sparse matrix with electronic potential integrals.
     */
    CElectronicPotentialMatrix(const CSparseMatrix& matrix);
    
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
