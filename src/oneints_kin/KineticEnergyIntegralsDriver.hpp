//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef KineticEnergyIntegralsDriver_hpp
#define KineticEnergyIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "ExecMode.hpp"
#include "GtoContainer.hpp"
#include "KineticEnergyMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionMap.hpp"
#include "SparseMatrix.hpp"
#include "VecIndexes.hpp"

/**
 Class CKineticEnergyIntegralsDriver computes one-electron kinetic energy
 integrals.

 @author Z. Rinkevicius
 */
class CKineticEnergyIntegralsDriver
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
     Comutes kinetic energy integrals for pair of GTOs containers.

     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the kinetic energy matrix object.
     */
    CKineticEnergyMatrix _compKineticEnergyIntegrals(const CGtoContainer* braGtoContainer, const CGtoContainer* ketGtoContainer) const;

    /**
     Computes kinetic energy integrals for specific pair of GTOs blocks and
     stores integrals into distribution buffer.

     @param distPattern the pointer to integrals distribution pattern.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compKineticEnergyForGtoBlocks(COneIntsDistribution* distPattern, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const;

    /**
     Computes batch of primitive kinetic energy integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).

     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).

     @param primBuffer the primitives buffer.
     @param recursionMap the recursion map for Obara-Saika recursion.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimKineticEnergyInts(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& abDistances,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto) const;

    /**
     Sets recursion map object for overlap integrals of specified angular
     momentum.

     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @return the recursion map for kinetic energy integrals.
     */
    CRecursionMap _setRecursionMap(const int32_t braAngularMomentum, const int32_t ketAngularMomentum, const int32_t maxNumberOfPrimitives) const;

   public:
    /**
     Creates a kinetic energy integrals driver object using MPI info.

     @param comm the MPI communicator.
     */
    CKineticEnergyIntegralsDriver(MPI_Comm comm);

    /**
     Destroys a kinetic energy integrals driver object.
     */
    ~CKineticEnergyIntegralsDriver();

    /**
     Computes kinetic energy integrals for molecule in specific basis set and
     stores results in kinetic energy matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @return the kinetic energy matrix object.
     */
    CKineticEnergyMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Computes kinetic energy integrals for molecule in two basis sets and stores
     results in kinetic matrix object.

     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of kinetic energy matrix.
     @param ketBasis the molecular basis for ket side of kinetic energy matrix.
     @return the kinetic energy matrix object.
     */
    CKineticEnergyMatrix compute(const CMolecule& molecule, const CMolecularBasis& braBasis, const CMolecularBasis& ketBasis) const;

    /**
     Computes kinetic energy integrals for two molecules in basis set and stores
     results in kinetic energy matrix object.

     @param braMolecule the molecule for bra side of kinetic energy matrix.
     @param ketMolecule the molecule for ket side of kinetic energy matrix.
     @param basis the molecular basis.
     @return the kinetic energy matrix object.
     */
    CKineticEnergyMatrix compute(const CMolecule& braMolecule, const CMolecule& ketMolecule, const CMolecularBasis& basis) const;

    /**
     Computes kinetic energy integrals for two molecules in different basis sets
     and stores results in kinetic energy matrix object.

     @param braMolecule the molecule for bra side of kinetic energy matrix.
     @param ketMolecule the molecule for ket side of kinetic energy matrix.
     @param braBasis the molecular basis for bra side of kinetic energy matrix.
     @param ketBasis the molecular basis for ket side of kinetic energy matrix.
     @return the kinetic energy matrix object.
     */
    CKineticEnergyMatrix compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis) const;

    /**
     Computes kinetic energy integrals blocks for pair of GTOs blocks and stores
     them into integrals batch.

     @param intsBatch the pointer to integrals batch buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(double* intsBatch, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const;
};

#endif /* KineticEnergyIntegralsDriver_hpp */
