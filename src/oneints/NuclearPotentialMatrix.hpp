//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef NuclearPotentialMatrix_hpp
#define NuclearPotentialMatrix_hpp

#include <cstdint>
#include <string>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"

/**
 Class CNuclearPotentialMatrix stores general nuclear potential matrix and
 provides set of methods for handling of nuclear potential matrix data.

 @author Z. Rinkevicius
 */
class CNuclearPotentialMatrix
{
    /**
     The generic dense nuclear potential matrix (rectangular or square).
     */
    CDenseMatrix _matrix;

   public:
    /**
     Creates an empty nuclear potential matrix object.
     */
    CNuclearPotentialMatrix();

    /**
     Creates a nuclear potential matrix object.

     @param matrix the dense matrix with nuclear potential integrals.
     */
    CNuclearPotentialMatrix(const CDenseMatrix& matrix);

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
     Reduces nuclear potential matrix object from all MPI process within domain
     of MPI communicator into nuclear potential matrix object on master node by
     summing them.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void reduce_sum(int32_t  rank,
                    int32_t  nodes,
                    MPI_Comm comm);

    /**
     Gets string representation of nuclear potential matrix.

     @return a string for printing the nuclear potential matrix.
     */
    std::string getString() const;

    /**
     Gets number of rows in nuclear potential matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in nuclear potential matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in nuclear potential matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets constant pointer to first element of nuclear potential matrix.

     @return the constant pointer to first element of nuclear potential matrix.
     */
    const double* values() const;

    /**
     Computes nuclear potential energy for specific AO density matrix.

     @param aoDensityMatrix the AO density matrix object.
     @param iDensityMatrix the index of AO density matrix in AO density matrix object.
     @return the nuclear potential energy.
     */
    double getNuclearPotentialEnergy(const CAODensityMatrix& aoDensityMatrix, const int32_t iDensityMatrix) const;

    /**
     Converts nuclear potential matrix object to text output and insert it into
     output text stream.

     @param output the output text stream.
     @param source the nuclear potential matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CNuclearPotentialMatrix& source);
};

#endif /* NuclearPotentialMatrix_hpp */
