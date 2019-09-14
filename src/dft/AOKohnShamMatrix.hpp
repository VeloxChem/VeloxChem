//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AOKohnShamMatrix_hpp
#define AOKohnShamMatrix_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"

/**
 Class CAOKohnShamMatrix stores set of AO Kohn-Sham matrices and provides set of methods
 for handling of AO Kohn-Sham matrices data.

 @author Z. Rinkevicius
 */
class CAOKohnShamMatrix
{
    /**
     The set of exchange-correlation matrices.
     */
    std::vector<CDenseMatrix> _xcMatrices;
    
    /**
     The flag indicating form of Kohn-Sham matrix: restricted or unrestricted.
     */
    bool _xcRestricted;
    
    /**
     The number of electrons obtained during integration of exchange-correlation matrices.
     */
    double _xcElectrons;
    
    /**
     The exchange-correlation energy associated with Kohn-Sham matrix.
     */
    double _xcEnergy;
    
public:
    
    /**
     Creates an empty AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix();
    
    /**
     Creates a AO Kohn-Sham matrix object.
     
     @param xcMatrices the set of exchange-correlation matrices.
     @param xcRestricted the flag indicating restricted on unrestricted form of Kohn-Sham matrix.
     @param xcElectrons the number of electrons obtained during integration of exchange-correlation matrices.
     @param xcEnergy the exchange-correlation energy associated with Kohn-Sham matrix.
     */
    CAOKohnShamMatrix(const std::vector<CDenseMatrix>& xcMatrices,
                      const bool                       xcRestricted,
                      const double                     xcElectrons,
                      const double                     xcEnergy);
    
    /**
    Creates a AO Kohn-Sham matrix object.
    
    @param nRows the number of rows in exchange-correlation matrices.
    @param nColumns the numner of columns in exchange-correlation matrices.
    @param xcRestricted the flag indicating restricted on unrestricted form of Kohn-Sham matrix.
    */
    CAOKohnShamMatrix(const int32_t nRows,
                      const int32_t nColumns,
                      const bool    xcRestricted);
    
    /**
     Creates a AO Kohn-Sham matrix object.
     
     @param aoDensityMatrix the AO density matrix.
     */
    CAOKohnShamMatrix(const CAODensityMatrix& aoDensityMatrix);
    
    /**
     Creates a AO Kohn-Sham matrix object by copying other AO Kohn-Sham matrix object.
     
     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix(const CAOKohnShamMatrix& source);
    
    /**
     Creates a AO Kohn-Sham matrix object by moving other AO Kohn-Sham matrix object.
     
     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix(CAOKohnShamMatrix&& source) noexcept;
    
    /**
     Destroys a AO Kohn-Sham matrix object.
     */
    ~CAOKohnShamMatrix();
    
    /**
     Assigns a AO Kohn-Sham matrix object by copying other AO Kohn-Sham matrix object.
     
     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix& operator=(const CAOKohnShamMatrix& source);
    
    /**
     Assigns a AO Kohn-Sham matrix object by moving other AO Kohn-Sham matrix object.
     
     @param source the AO Kohn-Sham matrix object.
     */
    CAOKohnShamMatrix& operator=(CAOKohnShamMatrix&& source) noexcept;
    
    /**
     Compares AO Kohn-Sham matrix object with other AO Kohn-Sham matrix object.
     
     @param other the AO Kohn-Sham matrix object.
     @return true if AO Kohn-Sham matrix objects are equal, false otherwise.
     */
    bool operator==(const CAOKohnShamMatrix& other) const;
    
    /**
     Compares AO Kohn-Sham matrix object with other AO Kohn-Sham matrix object.
     
     @param other the AO Kohn-Sham matrix object.
     @return true if AO Kohn-Sham matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CAOKohnShamMatrix& other) const;
    
    /**
     Sets number of electron obtained by integrating Kohn-Sham matrix.

     @param xcElectrons the number of electrons.
     */
    void setNumberOfElectrons(const double xcElectrons);
    
    /**
     Set exchange-correlation energy associated with Kohn-Sham matrix.

     @param xcEnergy the exchange-correlation energy.
     */
    void setExchangeCorrelationEnergy(const double xcEnergy);
    
    /**
     Resets all elements of AO Kohn-Sham matrix to zero.
     */
    void zero();
    
    /**
     Symmetrizes AO Kohn-Sham matrix according to it's type.
     */
    void symmetrize();
    
    /**
     Reduces AO Kohn-Sham matrix objects from all MPI process within domain of MPI
     communicator into AO Kohn-Sham matrix object on master node by summing them.
     
     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void reduce_sum(int32_t  rank,
                    int32_t  nodes,
                    MPI_Comm comm);

    /**
     Collect Kohn-Sham matrix object to master node of MPI communicator.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     @param source the rank of MPI process that provides the Kohn-Sham matrix.
     */
    void collect(int32_t rank, int32_t nodes, MPI_Comm comm, int32_t source);
    
    /**
     Checks if the Fock matrices are restricted.
     
     @return true if the Fock matrices are restricted.
     */
    bool isRestricted() const;
    
    
    /**
     Gets number of electrons obtaines by integrating Kohn-Sham matrix.

     @return the number of electrons.
     */
    double getNumberOfElectrons() const;
    
    /**
     Gets exchange-correlation energy associated with Kohn-Sham matrix.

     @return the exchange-correlation energy.
     */
    double getExchangeCorrelationEnergy() const;
    
    /**
     Gets number of rows in Kohn-Sham matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in Kohn-Sham matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in Kohn.Sham matrix.
     
     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of Kohn-Sham matrix.
     
     @param beta requires Kohn-Sham  matrix with beta spin.
     @return the constant pointer to first element of Kohn-Sham matrix.
     */
    const double* getKohnSham(const bool beta=false) const;
    
    /**
     Gets pointer to first element of Kohn-Sham matrix.
     
     @param beta requires Kohn-Sham matrix with beta spin.
     @return the pointer to first element of Kohn-Sham matrix.
     */
    double* getKohnSham(const bool beta=false);
    
    /**
     Gets constant reference to specific Kohn-Sham matrix.

     @param beta requires Kohn-Sham matrix with beta spin.
     @return the constant reference to Kohn-Sham matrix.
     */
    const CDenseMatrix& getReferenceToKohnSham(const bool beta=false) const;
    
    /**
     Gets string representation of Kohn-Sham matrix object.

     @return the string representation of Kohn-Shame matrix.
     */
    std::string getString() const;
    
    /**
     Converts AO Kohn-Sham matrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the AO Kohn-Sham matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&     output,
                                    const CAOKohnShamMatrix& source);
};

#endif /* AOKohnShamMatrix_hpp */
