//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricDipoleIntegralsDriver_hpp
#define ElectricDipoleIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "ExecMode.hpp"
#include "ElectricDipoleMatrix.hpp"
#include "GtoContainer.hpp"
#include "VecIndexes.hpp"
#include "OneIntsDistributor.hpp"

/**
 Class CElectricDipoleIntegralsDriver computes one-electron electric dipole
 integrals.
 
 @author Z. Rinkevicius
 */
class CElectricDipoleIntegralsDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;
    
    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;
    
    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;
    
    /**
     The Cartesian X coordinate of electric dipole origin.
     */
    double _xOrigin;
    
    /**
     The Cartesian Y coordinate of electric dipole origin.
     */
    double _yOrigin;
    
    /**
     The Cartesian Z coordinate of electric dipole origin.
     */
    double _zOrigin;
    
    /**
     Comutes electric dipole integrals for pair of GTOs containers.

     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix _compElectricDipoleIntegrals(const CGtoContainer* braGtoContainer,
                                                       const CGtoContainer* ketGtoContainer) const;
    
    /**
     Computes electric dipole integrals for specific pair of GTOs blocks and
     stores integrals in given distribution buffer.
     
     @param distPatternX the pointer to X component of integrals distribution
            pattern.
     @param distPatternY the pointer to Y component of integrals distribution
            pattern.
     @param distPatternZ the pointer to Z component of integrals distribution
            pattern.
     @param xOrigin the Cartesian X coordinate of electric dipole origin.
     @param yOrigin the Cartesian Y coordinate of electric dipole origin.
     @param zOrigin the Cartesian Z coordinate of electric dipole origin.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compElectricDipoleForGtoBlocks(      COneIntsDistribution* distPatternX,
                                               COneIntsDistribution* distPatternY,
                                               COneIntsDistribution* distPatternZ,
                                         const double                xOrigin,
                                         const double                yOrigin,
                                         const double                zOrigin,
                                         const CGtoBlock&            braGtoBlock,
                                         const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Computes batch of primitive electric dipole integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimElectricDipoleInts(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  abDistances,
                                     const CMemBlock2D<double>&  paDistances,
                                     const CMemBlock2D<double>&  pbDistances,
                                     const CMemBlock2D<double>&  pcDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoBlock&            ketGtoBlock,
                                     const int32_t               iContrGto) const;
    
    /**
     Gets Obara-Saika recursion pattern for specific combination of GTOs blocks
     on bra and ket sides.

     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @return the vector of two indexes object with recursion pattern.
     */
    CVecThreeIndexes _getRecursionPattern(const CGtoBlock& braGtoBlock,
                                          const CGtoBlock& ketGtoBlock) const;
    
    /**
     Gets vector of unified indexes of primitive GTOs buffer for specific
     Obara-Saika recursion pattern.

     @param recIndexes the vector of starting indexes of data blocks in recursion
            pattern.
     @param recPattern the recursion pattern.
     @param maxPrimGtos the maximum number of primitive GTOs in contracted
            GTO on bra side.
     @return the total number of blocks in recursion pattern.
     */
    int32_t _getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
                                           const CVecThreeIndexes&     recPattern,
                                           const int32_t               maxPrimGtos) const;
    
public:
    
    /**
     Creates a electric dipole integrals driver object using MPI info.
     
     @param comm the MPI communicator.
     */
    CElectricDipoleIntegralsDriver(MPI_Comm comm);
    
    /**
     Destroys a electric dipole integrals driver object.
     */
    ~CElectricDipoleIntegralsDriver();
    
    /**
     Sets origin of electric dipole.

     @param xOrigin the Cartesian X coordinate of electric dipole origin.
     @param yOrigin the Cartesian Y coordinate of electric dipole origin.
     @param zOrigin the Cartesian Z coordinate of electric dipole origin.
     */
    void setElectricDipoleOrigin(const double xOrigin,
                                 const double yOrigin,
                                 const double zOrigin);
    
    /**
     Computes electric dipole integrals for molecule in specific basis set and
     stores results in electric dipole matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule&       molecule,
                                  const CMolecularBasis& basis) const;
    
    /**
     Computes electric dipole integrals for molecule in two basis sets and stores
     results in electric dipole matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule&       molecule,
                                  const CMolecularBasis& braBasis,
                                  const CMolecularBasis& ketBasis) const;
    
    /**
     Computes electric dipole integrals for two molecules in basis set and
     stores results in electric dipole matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param basis the molecular basis.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule&       braMolecule,
                                  const CMolecule&       ketMolecule,
                                  const CMolecularBasis& basis) const;
    
    /**
     Computes electric dipole integrals for two molecules in different basis
     sets and stores results in electric dipole matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule&       braMolecule,
                                  const CMolecule&       ketMolecule,
                                  const CMolecularBasis& braBasis,
                                  const CMolecularBasis& ketBasis) const;
    
    /**
     Computes electric dipole integrals blocks for pair of GTOs blocks and
     stores them into integrals batch.
 
     @param intsBatchX the pointer to batch buffer of X component of electric
            dipole integrals.
     @param intsBatchY the pointer to batch buffer of Y component of electric
            dipole integrals.
     @param intsBatchZ the pointer to batch buffer of Z component of electric
            dipole integrals.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*    intsBatchX,
                       double*    intsBatchY,
                       double*    intsBatchZ,
                 const CGtoBlock& braGtoBlock,
                 const CGtoBlock& ketGtoBlock) const;
};

#endif /* ElectricDipoleIntegralsDriver_hpp */
