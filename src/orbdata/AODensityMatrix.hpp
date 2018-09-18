//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef AODensityMatrix_hpp
#define AODensityMatrix_hpp

#include <cstdint>
#include <vector>

#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"

/**
 Class CAODensityMatrix stores set of AO density matrices and provides
 set of methods for handling of AO density matrices data.
 
 @author Z. Rinkevicius
 */
class CAODensityMatrix
{
    /**
     The set of density matrices.
     */
    std::vector<CDenseMatrix> _denMatrices;
    
    /**
     The type of density matrices.
     */
    denmat _denType;
    
public:
    
    /**
     Creates an empty AO density matrix object.
     */
    CAODensityMatrix();
    
    /**
     Creates a AO density matrix object.
     
     @param denMatrices the set of density matrices.
     @param denType the type (restricted, unrestricted, etc) of density matrices.
     */
    CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices,
                     const denmat                     denType);
    
    /**
     Creates a AO density matrix object by copying other AO density matrix
     object.
     
     @param source the AO density matrix object.
     */
    CAODensityMatrix(const CAODensityMatrix& source);
    
    /**
     Creates a AO density matrix object by moving other AO density matrix
     object.
     
     @param source the AO density matrix object.
     */
    CAODensityMatrix(CAODensityMatrix&& source) noexcept;
    
    /**
     Destroys a AO density matrix object.
     */
    ~CAODensityMatrix();
    
    /**
     Assigns a AO density matrix object by copying other AO density matrix
     object.
     
     @param source the AO density matrix object.
     */
    CAODensityMatrix& operator=(const CAODensityMatrix& source);
    
    /**
     Assigns a AO density matrix object by moving other AO density matrix
     object.
     
     @param source the AO density matrix object.
     */
    CAODensityMatrix& operator=(CAODensityMatrix&& source) noexcept;
    
    /**
     Compares AO density matrix object with other AO density matrix object.
     
     @param other the AO density matrix object.
     @return true if AO density matrix objects are equal, false otherwise.
     */
    bool operator==(const CAODensityMatrix& other) const;
    
    /**
     Compares AO density matrix object with other AO density matrix object.
     
     @param other the AO density matrix object.
     @return true if AO density matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CAODensityMatrix& other) const;
    
    /**
     Gets number of rows in overlap matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the number of rows.
     */
    int32_t getNumberOfRows(const int32_t iDensityMatrix) const;
    
    /**
     Gets number of columns in specific matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the number of columns.
     */
    int32_t getNumberOfColumns(const int32_t iDensityMatrix) const;
    
    /**
     Gets number of elements in specific density matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the number of elements.
     */
    int32_t getNumberOfElements(const int32_t iDensityMatrix) const;
    
    /**
     Gets constant pointer to first element of specific restricted density
     matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of restricted density
     matrix.
     */
    const double* totalDensity(const int32_t iDensityMatrix) const;
    
    /**
     Gets constant pointer to first element of specific unrestricted alpha
     density matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of unrestricted alpha density
     matrix.
     */
    const double* alphaDensity(const int32_t iDensityMatrix) const;
    
    /**
     Gets constant pointer to first element of specific unrestricted beta
     density matrix.
     
     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of unrestricted beta density
     matrix.
     */
    const double* betaDensity(const int32_t iDensityMatrix) const;
    
    /**
     Converts AO density matrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the AO density matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&     output,
                                    const CAODensityMatrix& source);
};

#endif /* AODensityMatrix_hpp */
