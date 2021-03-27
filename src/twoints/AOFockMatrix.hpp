//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef AOFockMatrix_hpp
#define AOFockMatrix_hpp

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "FockMatrixType.hpp"
#include "KineticEnergyMatrix.hpp"
#include "MpiFunc.hpp"
#include "NuclearPotentialMatrix.hpp"

/**
 Class CAOFockMatrix stores set of AO Fock matrices and provides set of methods
 for handling of AO Fock matrices data.

 @author Z. Rinkevicius
 */
class CAOFockMatrix
{
    /**
     The set of Fock matrices.
     */
    std::vector<CDenseMatrix> _fockMatrices;

    /**
     The set of Fock matrices types.
     */
    std::vector<fockmat> _fockTypes;

    /**
     The scaling factors for exchange contributions.
     */
    std::vector<double> _scaleFactors;

    /**
     The set of identifiers of density matrices used to construct Fock matrices.
     */
    std::vector<int32_t> _idDensityMatrices;

    /**
     Gets number of matrices per Fock. (1: closed-shell; 2: open-shell)

     @return the number of matrices per Fock.
     */
    int32_t _getNumberOfMatricesPerFock() const;

    /**
     Gets the actual index of a matrix in _fockMatrices.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     */
    int32_t _getMatrixID(const int32_t iFockMatrix, const std::string& spin) const;

   public:
    /**
     Creates an empty AO Fock matrix object.
     */
    CAOFockMatrix();

    /**
     Creates a AO Fock matrix object.

     @param fockMatrices the set of Fock matrices.
     @param fockTypes the set of Fock matrices types.
     @param scaleFactors the set of scaling factors for exchange contributions.
     @param idDensityMatrices the set of identifiers used to construct Fock
            matrices.
     */
    CAOFockMatrix(const std::vector<CDenseMatrix>& fockMatrices,
                  const std::vector<fockmat>&      fockTypes,
                  const std::vector<double>&       scaleFactors,
                  const std::vector<int32_t>&      idDensityMatrices);

    /**
     Creates a AO Fock matrix object.

     @param aoDensityMatrix the AO density matrix.
     */
    CAOFockMatrix(const CAODensityMatrix& aoDensityMatrix);

    /**
     Creates a AO Fock matrix object by copying other AO Fock matrix object.

     @param source the AO Fock matrix object.
     */
    CAOFockMatrix(const CAOFockMatrix& source);

    /**
     Creates a AO Fock matrix object by moving other AO Fock matrix object.

     @param source the AO Fock matrix object.
     */
    CAOFockMatrix(CAOFockMatrix&& source) noexcept;

    /**
     Destroys a AO Fock matrix object.
     */
    ~CAOFockMatrix();

    /**
     Assigns a AO Fock matrix object by copying other AO Fock matrix object.

     @param source the AO Fock matrix object.
     */
    CAOFockMatrix& operator=(const CAOFockMatrix& source);

    /**
     Assigns a AO Fock matrix object by moving other AO Fock matrix object.

     @param source the AO Fock matrix object.
     */
    CAOFockMatrix& operator=(CAOFockMatrix&& source) noexcept;

    /**
     Compares AO Fock matrix object with other AO Fock matrix object.

     @param other the AO Fock matrix object.
     @return true if AO Fock matrix objects are equal, false otherwise.
     */
    bool operator==(const CAOFockMatrix& other) const;

    /**
     Compares AO Fock matrix object with other AO Fock matrix object.

     @param other the AO Fock matrix object.
     @return true if AO Fock matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CAOFockMatrix& other) const;

    /**
     Checks if AO Fock matrix is of closed-shell type.

     @return true if AO Fock matrix is of closed-shell type.
     */
    bool isClosedShell() const;

    /**
     Sets type of specific Fock matrix.

     @param fockType the type of Fock matrix.
     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     */
    void setFockType(const fockmat& fockType, const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA"));

    /**
     Sets scaling factor of exact exchange of specific Fock matrix.

     @param factor the scaling factor of exact exchange.
     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     */
    void setFockScaleFactor(const double factor, const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA"));

    /**
     Resets all elements of AO Fock matrix to zero.
     */
    void zero();

    /**
     Symmetrizes AO Fock matrix according to it's type.
     */
    void symmetrize();

    /**
     Add AO Fock matrix to AO Fock matrix.

     @param source the AO Fock matrix.
     */
    void add(const CAOFockMatrix& source);

    /**
     Adds core Hamiltonian, kinetic energy and nuclear potential matrices, to
     specific Fock matrix.

     @param kineticEnergyMatrix the kinetic energy matrix.
     @param nuclearPotentialMatrix the nuclear potential.
     @param iFockMatrix the index of Fock matrix.
     */
    void addCoreHamiltonian(const CKineticEnergyMatrix&    kineticEnergyMatrix,
                            const CNuclearPotentialMatrix& nuclearPotentialMatrix,
                            const int32_t                  iFockMatrix);

    /**
     Adds one electron operator matrix to specific Fock matrix.

     @param oneElectronMatrix the nuclear potential.
     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
    */
    void addOneElectronMatrix(const CDenseMatrix& oneElectronMatrix, const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA"));

    /**
     Scales specific Fock matrix by factor.

     @param factor the scaling factor.
     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     */
    void scale(const double factor, const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA"));

    /**
     Reduces AO Fock matrix objects from all MPI process within domain of MPI
     communicator into AO Fock matrix object on master node by summing them.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm);

    /**
     Gets number of Fock matrices.

     @return the number of Fock matrices.
     */
    int32_t getNumberOfFockMatrices() const;

    /**
     Gets number of rows in specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @return the number of rows.
     */
    int32_t getNumberOfRows(const int32_t iFockMatrix) const;

    /**
     Gets number of columns in specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @return the number of columns.
     */
    int32_t getNumberOfColumns(const int32_t iFockMatrix) const;

    /**
     Gets number of elements in specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @return the number of elements.
     */
    int32_t getNumberOfElements(const int32_t iFockMatrix) const;

    /**
     Gets constant pointer to first element of specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the constant pointer to first element of Fock matrix.
     */
    const double* getFock(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA")) const;

    /**
     Gets pointer to first element of specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the pointer to first element of Fock matrix.
     */
    double* getFock(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA"));

    /**
     Gets constant reference to specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the constant reference to Fock matrix.
     */
    const CDenseMatrix& getReferenceToFock(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA")) const;

    /**
     Gets type of specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the type of Fock matrix.
     */
    fockmat getFockType(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA")) const;

    /**
     Gets scaling factor of exchange contribution for specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the scaling factor.
     */
    double getScaleFactor(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA")) const;

    /**
     Gets identifier of AO density matrix used to construct specific Fock matrix.

     @param iFockMatrix the index of Fock matrix.
     @param spin the spin of Fock matrix.
     @return the identifier of density matrix.
     */
    int32_t getDensityIdentifier(const int32_t iFockMatrix, const std::string& spin = std::string("ALPHA")) const;

    /**
     Checks if specific Fock matrix is symmetric.

     @param iFockMatrix the index of Fock matrix.
     @return true if Fock matrix is symmetric, false otherwise.
     */
    bool isSymmetric(const int32_t iFockMatrix) const;

    /**
     Computes electronic energy for specific AO density matrix.

     @param iFockMatrix the index of Fock matrix
     @param aoDensityMatrix the AO density matrix object.
     @param iDensityMatrix the index of AO density matrix in AO density matrix object.
     @return the electronic energy.
     */
    double getElectronicEnergy(const int32_t iFockMatrix, const CAODensityMatrix& aoDensityMatrix, const int32_t iDensityMatrix) const;

    /**
     Gets string representation of density matrix object.

     @return the string representation of density matrix.
     */
    std::string getString() const;

    /**
     Converts AO Fock matrix object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the AO Fock matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CAOFockMatrix& source);
};

#endif /* AOFockMatrix_hpp */
