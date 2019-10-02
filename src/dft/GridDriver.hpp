//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef GridDriver_hpp
#define GridDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "ExecMode.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"

/**
 Class CGridDriver generates grid points data for usage in numerical
 integration.

 @author Z. Rinkevicius
 */
class CGridDriver
{
    /**
     The accuracy level of grid (from 1 to 6).
     */
    int32_t _gridLevel;

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
     The threshold of weights screening.
     */
    double _thresholdOfWeight;

    /**
     The execution mode of grid driver object.
     */
    execmode _runMode;

    /**
     Determines number of radial grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of radial points.
     */
    int32_t _getNumberOfRadialPoints(const int32_t idElemental) const;

    /**
     Determines number of angular grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of angular points.
     */
    int32_t _getNumberOfAngularPoints(const int32_t idElemental) const;

    /**
     Prints start header with grid generation settings to output stream.

     @param molecule the molecule.
     @return the output string.
     */
    std::string _startHeader(const CMolecule& molecule) const;

    /**
     Prints finish header with grid generation settings to output stream.

     @param molecularGrid the molecular grid object.
     @return the output string.
     */
    std::string _finishHeader(const CMolecularGrid& molecularGrid) const;

    /**
     Creates molecular grid on master node by generating fraction of grid
     points on each MPI process within domain of MPI communicator. Grid points
     are generated using only CPUs.

     @param molecule the molecule.
     @return the molecular grid object.
     */
    CMolecularGrid _genGridPointsOnCPU(const CMolecule& molecule) const;

    /**
     Gets size of grid points batch.

     @param idsElemental the vector of chemical elements identifiers.
     @param offset the  in vector of chemical elements identifiers.
     @param nAtoms the number of atoms in batch.
     @return the number of grid points.
     */
    int32_t _getBatchSize(const int32_t* idsElemental, const int32_t offset, const int32_t nAtoms) const;

    /**
     Generates grid points for specific atom in molecule.

     @param rawGridPoints the raw grid points.
     @param minDistance the distance to closest neighbouring atom.
     @param gridOffset the atom grid points offset in raw grid points.
     @param atomCoordinatesX the vector of Cartesian X coordinates of atoms.
     @param atomCoordinatesY the vector of Cartesian Y coordinates of atoms.
     @param atomCoordinatesZ the vector of Cartesian Z coordinates of atoms.
     @param nAtoms the number of atoms.
     @param idElemental the chemical element identifier of atom.
     @param idAtomic the index of atom.
     */
    void _genAtomGridPoints(CMemBlock2D<double>* rawGridPoints,
                            const double         minDistance,
                            const int32_t        gridOffset,
                            const double*        atomCoordinatesX,
                            const double*        atomCoordinatesY,
                            const double*        atomCoordinatesZ,
                            const int32_t        nAtoms,
                            const int32_t        idElemental,
                            const int32_t        idAtomic) const;

    /**
     Prunes raw grid points by screening weights and discarding all grid point
     with weights bellow cutoff threshold.

     @param rawGridPoints the raw grid points.
     @return the number of pruned grid points.
     */
    int32_t _screenRawGridPoints(CMemBlock2D<double>* rawGridPoints) const;

   public:
    /**
     Creates a grid driver object using MPI info.

     @param comm the MPI communicator.
     */
    CGridDriver(MPI_Comm comm);

    /**
     Destroys a grid driver object.
     */
    ~CGridDriver();

    /**
     Sets accuracy level for grid generation. Level: 1-6, where 1 is coarse
     grid, 5 is ultrafine grid, 6 special benchmarking grid.

     @param gridLevel the accuracy level of generated grid.
     */
    void setLevel(const int32_t gridLevel);

    /**
     Generates molecular grid for molecule. Errors are printed to output stream.
     Grid generation is distributed within domain of MPI communicator.

     @param molecule the molecule.
     @return the molecular grid object.
     */
    CMolecularGrid generate(const CMolecule& molecule) const;
};

#endif /* GridDriver_hpp */
