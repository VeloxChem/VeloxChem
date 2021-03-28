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

#ifndef ElectricDipoleIntegralsDriver_hpp
#define ElectricDipoleIntegralsDriver_hpp

#include <cstdint>

#include <mpi.h>

#include "ElectricDipoleMatrix.hpp"
#include "ExecMode.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionMap.hpp"

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
    CElectricDipoleMatrix _compElectricDipoleIntegrals(const CGtoContainer* braGtoContainer, const CGtoContainer* ketGtoContainer) const;

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
    void _compElectricDipoleForGtoBlocks(COneIntsDistribution* distPatternX,
                                         COneIntsDistribution* distPatternY,
                                         COneIntsDistribution* distPatternZ,
                                         const double          xOrigin,
                                         const double          yOrigin,
                                         const double          zOrigin,
                                         const CGtoBlock&      braGtoBlock,
                                         const CGtoBlock&      ketGtoBlock) const;

    /**
     Computes batch of primitive electric dipole integrals using Obara-Saika
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
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimElectricDipoleInts(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& abDistances,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pcDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto) const;

    /**
     Sets recursion map object for electric dipole integrals of specified angular
     momentum.

     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @return the recursion map for electric dipole integrals.
     */
    CRecursionMap _setRecursionMap(const int32_t braAngularMomentum, const int32_t ketAngularMomentum, const int32_t maxNumberOfPrimitives) const;

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
    void setElectricDipoleOrigin(const double xOrigin, const double yOrigin, const double zOrigin);

    /**
     Computes electric dipole integrals for molecule in specific basis set and
     stores results in electric dipole matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Computes electric dipole integrals for molecule in two basis sets and stores
     results in electric dipole matrix object.

     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule& molecule, const CMolecularBasis& braBasis, const CMolecularBasis& ketBasis) const;

    /**
     Computes electric dipole integrals for two molecules in basis set and
     stores results in electric dipole matrix object.

     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param basis the molecular basis.
     @return the electric dipole matrix object.
     */
    CElectricDipoleMatrix compute(const CMolecule& braMolecule, const CMolecule& ketMolecule, const CMolecularBasis& basis) const;

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
    void compute(double* intsBatchX, double* intsBatchY, double* intsBatchZ, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const;
};

#endif /* ElectricDipoleIntegralsDriver_hpp */
