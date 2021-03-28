//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef NuclearPotentialIntegralsDriver_hpp
#define NuclearPotentialIntegralsDriver_hpp

#include <cstdint>

#include <mpi.h>

#include "BoysFunction.hpp"
#include "GtoContainer.hpp"
#include "NuclearPotentialMatrix.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionMap.hpp"

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
    void _compNuclearPotentialForGtoBlocks(COneIntsDistribution*      distPattern,
                                           const CMemBlock<double>*   charges,
                                           const CMemBlock2D<double>* coordinates,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock) const;

    /**
     Adds single point charge contribution from primitives recursion buffer to
     primitives accumulation buffer, which contains primitive nuclear potential
     integrals.

     @param accBuffer the primitive integrals accumulation buffer.
     @param primBuffer the primitives recursion buffer.
     @param primIndex the index of primitive nuclear repulsion integrals in primitive integrals buffer.
     @param charges the vector of point charges.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointCharge the index of point charge in vector point charges.
     */
    void _addPointChargeContribution(CMemBlock2D<double>&       accBuffer,
                                     const CMemBlock2D<double>& primBuffer,
                                     const int32_t              primIndex,
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
     @param recursionMap the recursion map for Obara-Saika recursion.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimNuclearPotentialInts(CMemBlock2D<double>&       primBuffer,
                                       const CRecursionMap&       recursionMap,
                                       const CBoysFunction&       bfTable,
                                       CMemBlock<double>&         bfArguments,
                                       CMemBlock2D<double>&       bfValues,
                                       const int32_t              bfOrder,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& abDistances,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto) const;

    /**
     Sets recursion map object for nuclear potential integrals of specified angular
     momentum.

     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @return the recursion map for nuclear potential integrals.
     */
    CRecursionMap _setRecursionMap(const int32_t braAngularMomentum, const int32_t ketAngularMomentum, const int32_t maxNumberOfPrimitives) const;

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
    CNuclearPotentialMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Computes nuclear potential integrals for molecules in specific basis set
     and stores results in nuclear potential matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param pchgMolecule the molecule providing nuclear point charges.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecule& pchgMolecule) const;

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
     Computes nuclear potential integrals for a molecule in specific basis
     using pre-defined point charges and coordinates
     and stores results in nuclear potential matrix object.

     @param molecule the molecule for nuclear potential matrix.
     @param basis the molecular basis.
     @param charges the point charges
     @param coordinates the coordinates of the charges
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&           molecule,
                                    const CMolecularBasis&     basis,
                                    const CMemBlock<double>&   charges,
                                    const CMemBlock2D<double>& coordinates) const;

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
    void compute(double*                    intsBatch,
                 const CMemBlock<double>*   charges,
                 const CMemBlock2D<double>* coordinates,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock) const;
};

#endif /* NuclearPotentialIntegralsDriver_hpp */
