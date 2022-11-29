//
//                           VELOXCHEM 1.0-RC2
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

#ifndef VisualizationDriver_hpp
#define VisualizationDriver_hpp

#include <mpi.h>

#include <array>
#include <vector>

#include "AODensityMatrix.hpp"
#include "Buffer.hpp"
#include "CubicGrid.hpp"
#include "MemBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularOrbitals.hpp"
#include "Molecule.hpp"

/**
 Class CVisualizationDriver computes wavefunction or density on grid points.

 @author X. Li
 */
class CVisualizationDriver
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
     Builds the components of Cartesian angular momentum for a shell.

     @param angl the angular momentum of the shell.
     @return a vector of vector, e.g. {{1,0,0}, {0,1,0}, {0,0,1}} for p-shell.
     */
    std::vector<std::vector<int32_t>> _buildCartesianAngularMomentum(int angl) const;

    /**
     Computes atomic orbitals at a given grid point.

     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param xp the X cooridnate of the grid point.
     @param yp the Y cooridnate of the grid point.
     @param zp the Z cooridnate of the grid point.
     @return phi values of the atomic orbitals at the grid point.
     */
    std::vector<double> _compPhiAtomicOrbitals(const CMolecule&       molecule,
                                               const CMolecularBasis& basis,
                                               const double           xp,
                                               const double           yp,
                                               const double           zp) const;

    /**
     Computes molecular orbital values at rank-local cubic grid points (OpenMP).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param molorb the molecular orbitals of the molecule.
     @param moidx the index of the molecular orbital (0-based).
     @param mospin the spin of the molecular orbital ('alpha' or 'beta').
     */
    void _computeLocalGrid(CCubicGrid&               grid,
                           const CMolecule&          molecule,
                           const CMolecularBasis&    basis,
                           const CMolecularOrbitals& molorb,
                           const int32_t             moidx,
                           const std::string&        mospin) const;

    /**
     Computes electronic densities at rank-local cubic grid points (OpenMP).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the density matrix of the molecule.
     @param densityIndex the index of the density matrix (0-based).
     @param densitySpin the spin of the density matrix ('alpha' or 'beta').
     */
    void _computeLocalGrid(CCubicGrid&             grid,
                           const CMolecule&        molecule,
                           const CMolecularBasis&  basis,
                           const CAODensityMatrix& density,
                           const int32_t           densityIndex,
                           const std::string&      densitySpin) const;

    /** Computes values of orbital-like quantities at rank-local cubic grid points (OpenMP).
     *
     * @param grid the cubic grid.
     * @param molecule the molecule.
     * @param basis the basis set for the molecule.
     * @param coeffs AO-basis expansion coefficients of the orbital-like
     * quantities to plot. This is of N_AO x N_orb dimensions.
     * @param idxs indices of the orbital-like quantities to plot (0-based).
     */
    void _computeLocalGrid(CCubicGrid&                 grid,
                           const CMolecule&            molecule,
                           const CMolecularBasis&      basis,
                           const BufferHostXYd&        coeffs,
                           const std::vector<int32_t>& idxs) const;

   public:
    /**
     Creates a visualization driver object.

     @param comm the MPI communicator.
     */
    CVisualizationDriver(MPI_Comm comm);

    /**
     Gets atomic orbital information.

     @param molecule the molecule.
     @param basis the molecular basis set.
     @return a vector of vector containing atomic number (nuclear charge),
             angular momentum, spherical harmonic index and basis function
             index for each atomic orbital.
     */
    std::vector<std::vector<int32_t>> getAtomicOrbitalInformation(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Maps atom to atomic orbitals.

     @param molecule the molecule.
     @param basis the molecular basis set.
     @return a vector of vector that maps atom index to atomic orbital indices.
     */
    std::vector<std::vector<int32_t>> mapAtomToAtomicOrbitals(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Computes atomic orbital (centered at origin) at cubic grid points (OpenMP).

     @param grid the cubic grid.
     @param basis the molecular basis set.
     @param aoinfo the atomic orbital information.
     */
    void computeAtomicOrbitalForGrid(CCubicGrid& grid, const CMolecularBasis& basis, const std::vector<int32_t>& aoinfo) const;

    /**
     Gets rank of the MPI process.

     @return rank of the MPI process.
     */
    int32_t getRank() const;

    /**
     Computes counts and displacements for distributing workloads to MPI processes.

     @param nx number of tasks.
     @return a vector of vector containing counts and displacements.
     */
    std::vector<std::vector<int32_t>> getCountsAndDisplacements(const int32_t nx) const;

    /**
     Computes molecular orbital values at cubic grid points (MPI).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param molorb the molecular orbitals of the molecule.
     @param moidx the index of the molecular orbital (0-based).
     @param mospin the spin of the molecular orbital ('alpha' or 'beta').
     */
    void compute(CCubicGrid&               grid,
                 const CMolecule&          molecule,
                 const CMolecularBasis&    basis,
                 const CMolecularOrbitals& molorb,
                 const int32_t             moidx,
                 const std::string&        mospin) const;

    /**
     Computes density values at cubic grid points (MPI).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the AO density matrix.
     @param denidx the index of the density matrix (0-based).
     @param denspin the spin of the density matrix ('alpha' or 'beta').
     */
    void compute(CCubicGrid&             grid,
                 const CMolecule&        molecule,
                 const CMolecularBasis&  basis,
                 const CAODensityMatrix& density,
                 const int32_t           denidx,
                 const std::string&      denspin) const;

    /** Computes values of orbital-like quantities at cubic grid points (MPI).
     *
     * @param grid the cubic grid.
     * @param molecule the molecule.
     * @param basis the basis set for the molecule.
     * @param coeffs AO-basis expansion coefficients of the orbital-like
     * quantities to plot. This is of N_AO x N_orb dimensions.
     * @param idxs indices of the orbital-like quantities to plot (0-based).
     */
    void compute(CCubicGrid&                 grid,
                 const CMolecule&            molecule,
                 const CMolecularBasis&      basis,
                 const BufferHostXYd&        coeffs,
                 const std::vector<int32_t>& idxs) const;

    /**
     Computes molecular orbital at given coordinates.

     @param coords the coordinates.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param mo the molecular orbitals of the molecule.
     @param moidx the index of the molecular orbital (0-based).
     @param mospin the spin of the molecular orbital ('alpha' or 'beta').
     */
    std::vector<double> getMO(const std::vector<std::vector<double>>& coords,
                              const CMolecule&                        molecule,
                              const CMolecularBasis&                  basis,
                              const CMolecularOrbitals&               mo,
                              const int32_t                           moidx,
                              const std::string&                      mospin) const;

    /**
     Computes densities at given coordinates.

     @param coords the coordinates.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the density matrix of the molecule.
     @param denspin the spin of the density matrix ('alpha' or 'beta').
     */
    std::vector<double> getDensity(const std::vector<std::vector<double>>& coords,
                                   const CMolecule&                        molecule,
                                   const CMolecularBasis&                  basis,
                                   const CAODensityMatrix&                 density,
                                   const std::string&                      denspin) const;

    /**
     Computes one-particle density gamma(x1;x2) at given coordinates.

     @param coords_1 the coordinates of x1.
     @param coords_2 the coordinates of x2.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the density matrix of the molecule.
     @param spin_1 the spin of x1 ('alpha' or 'beta').
     @param spin_2 the spin of x2 ('alpha' or 'beta').
     */
    std::vector<double> getOneParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                              const std::vector<std::vector<double>>& coords_2,
                                              const CMolecule&                        molecule,
                                              const CMolecularBasis&                  basis,
                                              const CAODensityMatrix&                 density,
                                              const std::string&                      spin_1,
                                              const std::string&                      spin_2) const;

    /**
     Computes two-particle density Gamma(x1,x2;x1,x2) at given coordinates.

     @param coords_1 the coordinate of x1.
     @param coords_2 the coordinate of x2.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the density matrix of the molecule.
     @param spin_1 the spin of x1 ('alpha' or 'beta').
     @param spin_2 the spin of x2 ('alpha' or 'beta').
     */
    std::vector<double> getTwoParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                              const std::vector<std::vector<double>>& coords_2,
                                              const CMolecule&                        molecule,
                                              const CMolecularBasis&                  basis,
                                              const CAODensityMatrix&                 density,
                                              const std::string&                      spin_1,
                                              const std::string&                      spin_2) const;
};

#endif /* VisualizationDriver_hpp */
