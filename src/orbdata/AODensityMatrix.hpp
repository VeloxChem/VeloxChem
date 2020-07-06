//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef AODensityMatrix_hpp
#define AODensityMatrix_hpp

#include <cstdint>
#include <vector>

#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MpiFunc.hpp"

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
    CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices, const denmat denType);

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
     Sets AO density matrix type.

     @param denType the density matrix type.
     */
    void setDensityType(const denmat denType);

    /**
     Appends AO density matrix object to current AO density matrix object.

     @param other the AO density matrix object.
     */
    void append(const CAODensityMatrix& other);

    /**
     Creates difference AO density matrix between this AO density matrix and
     given AO density matrix.

     @param other the AO density matrix object.
     @return the AO density matrix object with difference between AO density
             matrices.
     */
    CAODensityMatrix sub(const CAODensityMatrix& other) const;

    /**
     Gets number of density matrices.

     @return the number of density matrices.
     */
    int32_t getNumberOfDensityMatrices() const;

    /**
     Gets total number of matrices stored in AO density matrix.

     @return the number of density matrices.
     */
    int32_t getNumberOfMatrices() const;

    /**
     Gets type of density matrix.

     @return the type of density matrix.
     */
    denmat getDensityType() const;

    /**
     Gets number of rows in specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the number of rows.
     */
    int32_t getNumberOfRows(const int32_t iDensityMatrix) const;

    /**
     Gets number of columns in specific density matrix.

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
     Gets constant pointer to first element of specific alpha density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of unrestricted alpha density
     matrix.
     */
    const double* alphaDensity(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to first element of specific beta density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of unrestricted beta density
     matrix.
     */
    const double* betaDensity(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to first element of specific density matrix.

     @param iDensityMatrix the index of density matrix.
     @return the constant pointer to first element of unrestricted beta density
     matrix.
     */
    const double* getDensity(const int32_t iDensityMatrix) const;

    /**
     Gets constant reference to density matrix as dense matrix object.

     @param iDensityMatrix the index of density matrix.
     @return the constant reference to density matrix.
     */
    const CDenseMatrix& getReferenceToDensity(const int32_t iDensityMatrix) const;

    /**
     Gets string representation of density matrix object.

     @return the string representation.
     */
    std::string getString() const;
    
    
    /**
     Checks if AO density matrix of spin restricted type.

     @return true if AO density of spin restricted type.
     */
    bool isRestricted() const;

    /**
     Broadcasts AO density matrix object within domain of MPI communicator.

     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);

    /**
     Converts AO density matrix object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the AO density matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CAODensityMatrix& source);
};

#endif /* AODensityMatrix_hpp */
