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

#ifndef ElectricFieldIntegralsDriver_hpp
#define ElectricFieldIntegralsDriver_hpp

#include <mpi.h>

#include <cstdint>

#include "BoysFunction.hpp"
#include "ElectricFieldMatrix.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionMap.hpp"

/**
 Class CElectricFieldIntegralsDriver computes one-electron electric field
 integrals.

 @author Z. Rinkevicius
 */
class CElectricFieldIntegralsDriver
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
     The local MPI communicator.
     */
    MPI_Comm _locComm;

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
    void _compElectricFieldForGtoBlocks(COneIntsDistribution*      distPatternX,
                                        COneIntsDistribution*      distPatternY,
                                        COneIntsDistribution*      distPatternZ,
                                        const CMemBlock2D<double>* dipoles,
                                        const CMemBlock2D<double>* coordinates,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock) const;

    /**
     Adds single point dipole contribution from primitives recursion buffer to
     primitives accumulation buffer, which contains primitive electric field
     integrals.

     @param accBuffer the primitive integrals accumulation buffer.
     @param primBuffer the primitives recursion buffer.
     @param primIndex the index of specific integrals in primitives recursion
     buffer.
     @param dipoles the vector of point dipoles.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointDipole the index of point dipole in vector of point dipoles.
     */
    void _addPointDipoleContribution(CMemBlock2D<double>&       accBuffer,
                                     const CMemBlock2D<double>& primBuffer,
                                     const int32_t              primIndex,
                                     const CMemBlock2D<double>& dipoles,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto,
                                     const int32_t              iPointDipole) const;

    /**
     Computes batch of primitive electric field integrals using Obara-Saika
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
    void _compPrimElectricFieldInts(CMemBlock2D<double>&       primBuffer,
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
     Sets recursion map object for electric field integrals of specified angular
     momentum.

     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @return the recursion map for electric dipole integrals.
     */
    CRecursionMap _setRecursionMap(const int32_t braAngularMomentum, const int32_t ketAngularMomentum, const int32_t maxNumberOfPrimitives) const;

   public:
    /**
     Creates a electric field integrals driver object using MPI info.

     @param comm the MPI communicator.
     */
    CElectricFieldIntegralsDriver(MPI_Comm comm);

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
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ) const;

    /**
     Computes electric field integrals for molecules in specific basis set for
     given set of point dipoles and stores results in electric field matrix
     object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           molecule,
                                 const CMolecularBasis&     basis,
                                 const CMemBlock2D<double>& dipoles,
                                 const CMemBlock2D<double>& coordinates) const;

    /**
     Computes electric field integrals for molecule in two basis sets at
     specified position and stores results in electric field matrix object.

     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param coordinateX the Cartesian X coordinate of electric field position.
     @param coordinateY the Cartesian Y coordinate of electric field position.
     @param coordinateZ the Cartesian Z coordinate of electric field position.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       molecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ) const;

    /**
     Computes electric field integrals for molecule in two basis sets for given
     set of point dipoles and stores results in electric field matrix object.

     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of electric field matrix.
     @param ketBasis the molecular basis for ket side of electric field matrix.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           molecule,
                                 const CMolecularBasis&     braBasis,
                                 const CMolecularBasis&     ketBasis,
                                 const CMemBlock2D<double>& dipoles,
                                 const CMemBlock2D<double>& coordinates) const;

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
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& basis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ) const;

    /**
     Computes electric field integrals for two molecules in specific basis
     set for set of point dipoles and stores results in electric field matrix
     object.

     @param braMolecule the molecule for bra side of electric field matrix.
     @param ketMolecule the molecule for ket side of electric field matrix.
     @param basis the molecular basis.
     @param dipoles the vector of point dipoles.
     @param coordinates the vector of point dipoles coordinates.
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           braMolecule,
                                 const CMolecule&           ketMolecule,
                                 const CMolecularBasis&     basis,
                                 const CMemBlock2D<double>& dipoles,
                                 const CMemBlock2D<double>& coordinates) const;

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
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                 const double           coordinateX,
                                 const double           coordinateY,
                                 const double           coordinateZ) const;

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
     @return the electric field matrix object.
     */
    CElectricFieldMatrix compute(const CMolecule&           braMolecule,
                                 const CMolecule&           ketMolecule,
                                 const CMolecularBasis&     braBasis,
                                 const CMolecularBasis&     ketBasis,
                                 const CMemBlock2D<double>& dipoles,
                                 const CMemBlock2D<double>& coordinates) const;

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
    void compute(double*                    intsBatchX,
                 double*                    intsBatchY,
                 double*                    intsBatchZ,
                 const CMemBlock2D<double>* dipoles,
                 const CMemBlock2D<double>* coordinates,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock) const;
};

#endif /* ElectricFieldIntegralsDriver_hpp */
