//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricFieldIntegralsDriver_hpp
#define ElectricFieldIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "MemBlock2D.hpp"
#include "GtoContainer.hpp"
#include "ElectricFieldMatrix.hpp"
#include "OneIntsDistributor.hpp"
#include "VecIndexes.hpp"

/**
 Class CElectricFieldIntegralsDriver computes one-electron electric field
 integrals.
 
 @author Z. Rinkevicius
 */
class CElectricFieldIntegralsDriver
{
    /**
     The rank of associated global MPI process.
     */
    int32_t _globRank;
    
    /**
     The total number of global MPI processes.
     */
    int32_t _globNodes;
    
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;
    
    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;
    
    /**
     The flag for local execution mode.
     */
    bool _isLocalMode;
    
    /**
     Comutes electric field integrals for pair of GTOs containers.
     
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordines.
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix _compElectricFieldIntegrals(const CMemBlock2D<double>* dipoles,
                                                     const CMemBlock2D<double>* coordinates,
                                                     const CGtoContainer*       braGtoContainer,
                                                     const CGtoContainer*       ketGtoContainer) const;
    
    /**
     Computes electric field integrals for specific pair of GTOs blocks.
     
     @param distPatternX the pointer to X component of integrals distribution
            pattern.
     @param distPatternY the pointer to Y component of integrals distribution
            pattern.
     @param distPatternZ the pointer to Z component of integrals distribution
            pattern.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compElectricFieldForGtoBlocks(      COneIntsDistribution* distPatternX,
                                              COneIntsDistribution* distPatternY,
                                              COneIntsDistribution* distPatternZ,
                                        const CMemBlock2D<double>*  dipoles,
                                        const CMemBlock2D<double>*  coordinates,
                                        const CGtoBlock&            braGtoBlock,
                                        const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Gets Obara-Saika recursion pattern for specific combination of GTOs blocks
     on bra and ket sides.
     
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @return the vector of four indexes object with recursion pattern.
     */
    CVecFourIndexes _getRecursionPattern(const CGtoBlock& braGtoBlock,
                                         const CGtoBlock& ketGtoBlock) const;
    
public:
    
    /**
     Creates a electric field integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CElectricFieldIntegralsDriver(const int32_t  globRank,
                                  const int32_t  globNodes,
                                        MPI_Comm comm);
    
    /**
     Destroys a electric field integrals driver object.
     */
    ~CElectricFieldIntegralsDriver();
    
    /**
     Computes electric field integrals for molecule in specific basis set at
     specified position and stores results in electric field matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param coordinateX the Cartesian X coordinate of electric field position.
     @param coordinateY the Cartesian Y coordinate of electric field position.
     @param coordinateZ the Cartesian Z coordinate of electric field position.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ,
                                       MPI_Comm         comm) const;
    
    /**
     Computes electric field integrals for molecules in specific basis set for
     given set of point dipoles and stores results in nuclear potential matrix
     object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param comm the MPI communicator.
     @return the nuclear potential matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           molecule,
                                 const CMolecularBasis&     basis,
                                 const CMemBlock2D<double>* dipoles,
                                 const CMemBlock2D<double>* coordinates,
                                       MPI_Comm             comm) const;
    
    /**
     Computes electric field integrals for molecule in two basis sets at
     specified position and stores results in electric field matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param coordinateX the Cartesian X coordinate of electric field position.
     @param coordinateY the Cartesian Y coordinate of electric field position.
     @param coordinateZ the Cartesian Z coordinate of electric field position.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       molecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ,
                                       MPI_Comm         comm) const;
    
    /**
     Computes electric field integrals for molecule in two basis sets for given
     set of point dipoles and stores results in electric field matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           molecule,
                                 const CMolecularBasis&     braBasis,
                                 const CMolecularBasis&     ketBasis,
                                 const CMemBlock2D<double>* dipoles,
                                 const CMemBlock2D<double>* coordinates,
                                       MPI_Comm             comm) const;
    
    /**
     Computes electric field integrals for two molecules in specific basis
     set at specific position and stores results in electric field matrix
     object.
     
     @param braMolecule the molecule for bra side of electric field matrix.
     @param ketMolecule the molecule for ket side of electric field matrix.
     @param basis the molecular basis.
     @param coordinateX the Cartesian X coordinate of electric field position.
     @param coordinateY the Cartesian Y coordinate of electric field position.
     @param coordinateZ the Cartesian Z coordinate of electric field position.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& basis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ,
                                       MPI_Comm         comm) const;
    
    /**
     Computes electric field integrals for two molecules in specific basis
     set for set of point dipoles and stores results in electric field matrix
     object.
     
     @param braMolecule the molecule for bra side of electric field matrix.
     @param ketMolecule the molecule for ket side of electric field matrix.
     @param basis the molecular basis.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           braMolecule,
                                 const CMolecule&           ketMolecule,
                                 const CMolecularBasis&     basis,
                                 const CMemBlock2D<double>* dipoles,
                                 const CMemBlock2D<double>* coordinates,
                                       MPI_Comm             comm) const;
    
    /**
     Computes electric field integrals for two molecules in two basis sets
     at specified position and stores results in electric field matrix object.
     
     @param braMolecule the molecule for bra side of electric field matrix.
     @param ketMolecule the molecule for ket side of electric field matrix.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param coordinateX the Cartesian X coordinate of electric field position.
     @param coordinateY the Cartesian Y coordinate of electric field position.
     @param coordinateZ the Cartesian Z coordinate of electric field position.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ,
                                       MPI_Comm         comm) const;
    
    /**
     Computes electric field integrals for two molecules in two basis sets
     for set of point dipoles and stores results in electric field matrix
     object.
     
     @param braMolecule the molecule for bra side of electric field matrix.
     @param ketMolecule the molecule for ket side of electric field matrix.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param comm the MPI communicator.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           braMolecule,
                                 const CMolecule&           ketMolecule,
                                 const CMolecularBasis&     braBasis,
                                 const CMolecularBasis&     ketBasis,
                                 const CMemBlock2D<double>* dipoles,
                                 const CMemBlock2D<double>* coordinates,
                                       MPI_Comm             comm) const;
    
    /**
     Computes electric field integrals blocks for pair of GTOs blocks and
     stores them into integrals batch.
     
     @param intsBatchX the pointer to batch buffer of X component of electric
            field integrals.
     @param intsBatchY the pointer to batch buffer of Y component of electric
            field integrals.
     @param intsBatchZ the pointer to batch buffer of Z component of electric
            field integrals.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*              intsBatchX,
                       double*              intsBatchY,
                       double*              intsBatchZ,
                 const CMemBlock2D<double>* dipoles,
                 const CMemBlock2D<double>* coordinates,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock) const;
};

#endif /* ElectricFieldIntegralsDriver_hpp */
