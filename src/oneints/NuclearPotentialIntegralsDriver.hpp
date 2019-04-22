//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef NuclearPotentialIntegralsDriver_hpp
#define NuclearPotentialIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "NuclearPotentialMatrix.hpp"
#include "GtoContainer.hpp"
#include "BoysFunction.hpp"
#include "OneIntsDistributor.hpp"

/**
 Class CNuclearPotentialIntegralsDriver computes one-electron nuclear potential
 integrals.
 
 @author Z. Rinkevicius
 */
class CNuclearPotentialIntegralsDriver
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
     Comutes nuclear potential integrals for pair of GTOs containers.
     
     @param charges the vector of point charges.
     @param coordinates the vector of point charges coordines.
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix _compNuclearPotentialIntegrals(const CMemBlock<double>*   charges,
                                                           const CMemBlock2D<double>* coordinates,
                                                           const CGtoContainer*       braGtoContainer,
                                                           const CGtoContainer*       ketGtoContainer) const;
    
    /**
     Computes nuclear potential integrals for specific pair of GTOs blocks.
     
     @param distPattern the pointer to integrals distribution pattern.
     @param charges the vector of point charges.
     @param coordinates the vector of point charges coordines.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compNuclearPotentialForGtoBlocks(      COneIntsDistribution* distPattern,
                                           const CMemBlock<double>*    charges,
                                           const CMemBlock2D<double>*  coordinates,
                                           const CGtoBlock&            braGtoBlock,
                                           const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Adds single point charge contribution from primitives recursion buffer to
     primitives accumulation buffer, which contains primitive nuclear potential
     integrals.

     @param accBuffer the primitive integrals accumulation buffer.
     @param primBuffer the primitives recursion buffer.
     @param charges the vector of point charges.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointCharge the index of point charge in vector point charges.
     */
    void _addPointChargeContribution(      CMemBlock2D<double>& accBuffer,
                                     const CMemBlock2D<double>& primBuffer,
                                     const CMemBlock<double>&   charges,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto,
                                     const int32_t              iPointCharge) const;
    
    /**
     Computes batch of primitive nuclear potential integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).
     
     @param primBuffer the primitives buffer.
     @param auxBuffer the auxilary integrals buffer.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param pcComponents the number of tensor components of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimNuclearPotentialInts(      CMemBlock2D<double>&  primBuffer,
                                             CMemBlock2D<double>&  auxBuffer,
                                       const CBoysFunction&        bfTable,
                                             CMemBlock<double>&    bfArguments,
                                             CMemBlock2D<double>&  bfValues,
                                       const int32_t               bfOrder,
                                       const CMemBlock2D<double>&  osFactors,
                                       const CMemBlock2D<double>&  abDistances,
                                       const CMemBlock2D<double>&  paDistances,
                                       const CMemBlock2D<double>&  pbDistances,
                                       const CMemBlock2D<double>&  pcDistances,
                                       const int32_t               pcComponents, 
                                       const CGtoBlock&            braGtoBlock,
                                       const CGtoBlock&            ketGtoBlock,
                                       const int32_t               iContrGto) const;
    
public:
    
    /**
     Creates a nuclear potential integrals driver object using MPI info.
     
     @param comm the MPI communicator.
     */
    CNuclearPotentialIntegralsDriver(MPI_Comm comm);
    
    /**
     Destroys a nuclear potential integrals driver object.
     */
    ~CNuclearPotentialIntegralsDriver();
    
    /**
     Computes nuclear potential integrals for molecule in specific basis set and
     stores results in nuclear potential matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       molecule,
                                    const CMolecularBasis& basis) const;
    
    /**
     Computes nuclear potential integrals for molecules in specific basis set
     and stores results in nuclear potential matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param pchgMolecule the molecule providing nuclear point charges.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       molecule,
                                    const CMolecularBasis& basis,
                                    const CMolecule&       pchgMolecule) const;

    /**
     Computes nuclear potential integrals for molecules in two basis sets and
     stores results in nuclear potential matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of nuclear potential matrix.
     @param ketBasis the molecular basis for ket side of nuclear potential matrix.
     @param pchgMolecule the molecule providing nuclear point charges.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       molecule,
                                    const CMolecularBasis& braBasis,
                                    const CMolecularBasis& ketBasis,
                                    const CMolecule&       pchgMolecule) const;

    /**
     Computes nuclear potential integrals for two molecules in specific basis
     set and stores results in nuclear potential matrix object.
     
     @param braMolecule the molecule for bra side of nuclear potential matrix.
     @param ketMolecule the molecule for ket side of nuclear potential matrix.
     @param basis the molecular basis.
     @param pchgMolecule the molecule providing nuclear point charges.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       braMolecule,
                                    const CMolecule&       ketMolecule,
                                    const CMolecularBasis& basis,
                                    const CMolecule&       pchgMolecule) const;

    /**
     Computes nuclear potential integrals for two molecules in two basis sets
     and stores results in nuclear potential matrix object.
     
     @param braMolecule the molecule for bra side of nuclear potential matrix.
     @param ketMolecule the molecule for ket side of nuclear potential matrix.
     @param braBasis the molecular basis for bra side of nuclear potential matrix.
     @param ketBasis the molecular basis for ket side of nuclear potential matrix.
     @param pchgMolecule the molecule providing nuclear point charges.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       braMolecule,
                                    const CMolecule&       ketMolecule,
                                    const CMolecularBasis& braBasis,
                                    const CMolecularBasis& ketBasis,
                                    const CMolecule&       pchgMolecule) const;

    /**
     Computes nuclear potential integrals blocks for pair of GTOs blocks and
     stores them into integrals batch.
     
     @param intsBatch the pointer to integrals batch buffer.
     @param charges the vector of point charges.
     @param coordinates the vector of point charges coordines.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*              intsBatch,
                 const CMemBlock<double>*   charges,
                 const CMemBlock2D<double>* coordinates,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock) const;
};

#endif /* NuclearPotentialIntegralsDriver_hpp */
