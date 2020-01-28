//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef VisualizationDriver_hpp
#define VisualizationDriver_hpp

#include <vector>

#include "AODensityMatrix.hpp"
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

   public:
    /**
     Creates a visualization driver object.

     @param comm the MPI communicator.
     */
    CVisualizationDriver(MPI_Comm comm);

    /**
     Destroys a visualization driver object.
     */
    ~CVisualizationDriver();

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
     @param mospin the spin of the molecular orbital ('a' or 'b').
     */
    void compute(CCubicGrid&               grid,
                 const CMolecule&          molecule,
                 const CMolecularBasis&    basis,
                 const CMolecularOrbitals& molorb,
                 const int32_t             moidx,
                 const std::string&        mospin) const;

    /**
     Computes molecular orbital values at cubic grid points (MPI).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the AO density matrix.
     @param denidx the index of the density matrix (0-based).
     @param denspin the spin of the density matrix ('a' or 'b').
     */
    void compute(CCubicGrid&             grid,
                 const CMolecule&        molecule,
                 const CMolecularBasis&  basis,
                 const CAODensityMatrix& density,
                 const int32_t           denidx,
                 const std::string&      denspin) const;

    /**
     Computes molecular orbital values at cubic grid points (OpenMP).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param molorb the molecular orbitals of the molecule.
     @param moidx the index of the molecular orbital (0-based).
     @param mospin the spin of the molecular orbital ('a' or 'b').
     */
    void compute_omp(CCubicGrid&               grid,
                     const CMolecule&          molecule,
                     const CMolecularBasis&    basis,
                     const CMolecularOrbitals& molorb,
                     const int32_t             moidx,
                     const std::string&        mospin) const;

    /**
     Computes electronic densities at cubic grid points (OpenMP).

     @param grid the cubic grid.
     @param molecule the molecule.
     @param basis the basis set for the molecule.
     @param density the density matrix of the molecule.
     @param densityIndex the index of the density matrix (0-based).
     @param densitySpin the spin of the density matrix ('a' or 'b').
     */
    void compute_omp(CCubicGrid&             grid,
                     const CMolecule&        molecule,
                     const CMolecularBasis&  basis,
                     const CAODensityMatrix& density,
                     const int32_t           densityIndex,
                     const std::string&      densitySpin) const;
};

#endif /* VisualizationDriver_hpp */
