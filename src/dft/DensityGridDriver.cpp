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

#include "DensityGridDriver.hpp"

#include <mpi.h>

#include <cmath>

#include "AODensityMatrix.hpp"
#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGrid.hpp"
#include "ExecMode.hpp"
#include "GenFunc.hpp"
#include "GtoContainer.hpp"
#include "GtoFunc.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MpiFunc.hpp"
#include "OMPTasks.hpp"
#include "SphericalMomentum.hpp"
#include "VecMemBlocks.hpp"
#include "XCFuncType.hpp"

CDensityGridDriver::CDensityGridDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;

    _thresholdOfDensity = 1.0e-13;

    _runMode = execmode::cpu;
}

CDensityGridDriver::~CDensityGridDriver()
{
}

CDensityGrid
CDensityGridDriver::generate(const CAODensityMatrix& aoDensityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const xcfun             xcFunctional)
{
    // initialize density grid

    CDensityGrid dgrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(), xcFunctional, dengrid::ab);

    dgrid.zero();

    // execution mode: CPU

    if (_runMode == execmode::cpu)
    {
        _genDensityGridOnCPU(dgrid, aoDensityMatrix, molecule, basis, molecularGrid, xcFunctional);
    }

    // execution mode: CPU/GPU

    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }

    return dgrid;
}

void
CDensityGridDriver::_genDensityGridOnCPU(CDensityGrid&           densityGrid,
                                         const CAODensityMatrix& aoDensityMatrix,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CMolecularGrid&   molecularGrid,
                                         const xcfun             xcFunctional)
{
    if (aoDensityMatrix.isClosedShell())
    {
        if (xcFunctional == xcfun::lda)
        {
            _genRestrictedDensityForLda(densityGrid, aoDensityMatrix, molecule, basis, molecularGrid);

            return;
        }

        if (xcFunctional == xcfun::gga)
        {
            _genRestrictedDensityForGga(densityGrid, aoDensityMatrix, molecule, basis, molecularGrid);

            return;
        }

        if (xcFunctional == xcfun::mgga)
        {
            _genRestrictedDensityForMgga(densityGrid, aoDensityMatrix, molecule, basis, molecularGrid);

            return;
        }
    }
    else
    {
        if (xcFunctional == xcfun::lda)
        {
            _genUnrestrictedDensityForLda(densityGrid, aoDensityMatrix, molecule, basis, molecularGrid);

            return;
        }

        if (xcFunctional == xcfun::gga)
        {
            _genUnrestrictedDensityForGga(densityGrid, aoDensityMatrix, molecule, basis, molecularGrid);

            return;
        }
    }
}

void
CDensityGridDriver::_genUnrestrictedDensityForLda(CDensityGrid&           densityGrid,
                                                  const CAODensityMatrix& aoDensityMatrix,
                                                  const CMolecule&        molecule,
                                                  const CMolecularBasis&  basis,
                                                  const CMolecularGrid&   molecularGrid) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfUnrestrictedDensityGridPointsForLda(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    delete gtovec;
}

void
CDensityGridDriver::_genBatchOfUnrestrictedDensityGridPointsForLda(CDensityGrid*           densityGrid,
                                                                   const CAODensityMatrix* aoDensityMatrix,
                                                                   const CGtoContainer*    gtoContainer,
                                                                   const double*           gridCoordinatesX,
                                                                   const double*           gridCoordinatesY,
                                                                   const double*           gridCoordinatesZ,
                                                                   const int32_t           gridOffset,
                                                                   const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distUnrestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distUnrestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_distUnrestrictedDensityValuesForLda(CDensityGrid*              densityGrid,
                                                         const CAODensityMatrix*    aoDensityMatrix,
                                                         const CMemBlock2D<double>& gtoValues,
                                                         const int32_t              gridOffset,
                                                         const int32_t              gridBlockPosition,
                                                         const int32_t              nGridPoints) const
{
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();

    for (int32_t i = 0; i < ndmat; i++)
    {
        // set up pointer to density grid data

        auto rhoa = densityGrid->alphaDensity(i);

        auto rhob = densityGrid->betaDensity(i);

        // set up poiinter to density matrix data

        auto denmat_a = aoDensityMatrix->alphaDensity(i);

        auto denmat_b = aoDensityMatrix->betaDensity(i);

        auto naos = aoDensityMatrix->getNumberOfRows(i);

        // loop over density matrix

        for (int32_t j = 0; j < naos; j++)
        {
            auto bgaos = gtoValues.data(j);

            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);

                const auto jkidx = j * naos + k;

                const auto kjidx = k * naos + j;

                auto dval_a = (j == k) ? denmat_a[jkidx] : denmat_a[kjidx] + denmat_a[kjidx];

                auto dval_b = (j == k) ? denmat_b[jkidx] : denmat_b[kjidx] + denmat_b[kjidx];

                if (std::fabs(dval_a + dval_b) > _thresholdOfDensity)
                {
                    #pragma omp simd
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        const double fact = bgaos[l] * kgaos[l];

                        rhoa[gridOffset + gridBlockPosition + l] += dval_a * fact;

                        rhob[gridOffset + gridBlockPosition + l] += dval_b * fact;
                    }
                }
            }
        }
    }
}

void
CDensityGridDriver::_genRestrictedDensityForLda(CDensityGrid&           densityGrid,
                                                const CAODensityMatrix& aoDensityMatrix,
                                                const CMolecule&        molecule,
                                                const CMolecularBasis&  basis,
                                                const CMolecularGrid&   molecularGrid) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfRestrictedDensityGridPointsForLda(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.updateBetaDensities();

    // destroy GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genRestrictedDensityForGga(CDensityGrid&           densityGrid,
                                                const CAODensityMatrix& aoDensityMatrix,
                                                const CMolecule&        molecule,
                                                const CMolecularBasis&  basis,
                                                const CMolecularGrid&   molecularGrid) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfRestrictedDensityGridPointsForGga(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.updateBetaDensities();

    densityGrid.computeDensityNorms();

    // destroy GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genRestrictedDensityForMgga(CDensityGrid&           densityGrid,
                                                 const CAODensityMatrix& aoDensityMatrix,
                                                 const CMolecule&        molecule,
                                                 const CMolecularBasis&  basis,
                                                 const CMolecularGrid&   molecularGrid) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfRestrictedDensityGridPointsForMgga(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.updateBetaDensities();

    densityGrid.computeDensityNorms();

    // destroy GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genUnrestrictedDensityForGga(CDensityGrid&           densityGrid,
                                                  const CAODensityMatrix& aoDensityMatrix,
                                                  const CMolecule&        molecule,
                                                  const CMolecularBasis&  basis,
                                                  const CMolecularGrid&   molecularGrid) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfUnrestrictedDensityGridPointsForGga(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.computeDensityNorms();

    // destroy GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genBatchOfRestrictedDensityGridPointsForLda(CDensityGrid*           densityGrid,
                                                                 const CAODensityMatrix* aoDensityMatrix,
                                                                 const CGtoContainer*    gtoContainer,
                                                                 const double*           gridCoordinatesX,
                                                                 const double*           gridCoordinatesY,
                                                                 const double*           gridCoordinatesZ,
                                                                 const int32_t           gridOffset,
                                                                 const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_genBatchOfRestrictedDensityGridPointsForGga(CDensityGrid*           densityGrid,
                                                                 const CAODensityMatrix* aoDensityMatrix,
                                                                 const CGtoContainer*    gtoContainer,
                                                                 const double*           gridCoordinatesX,
                                                                 const double*           gridCoordinatesY,
                                                                 const double*           gridCoordinatesZ,
                                                                 const int32_t           gridOffset,
                                                                 const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForGGA(
                gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        gtorec::computeGtosValuesForGGA(
            gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_genBatchOfRestrictedDensityGridPointsForMgga(CDensityGrid*           densityGrid,
                                                                  const CAODensityMatrix* aoDensityMatrix,
                                                                  const CGtoContainer*    gtoContainer,
                                                                  const double*           gridCoordinatesX,
                                                                  const double*           gridCoordinatesY,
                                                                  const double*           gridCoordinatesZ,
                                                                  const int32_t           gridOffset,
                                                                  const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForMGGA(
                gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForMgga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        gtorec::computeGtosValuesForMGGA(
            gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForMgga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_genBatchOfUnrestrictedDensityGridPointsForGga(CDensityGrid*           densityGrid,
                                                                   const CAODensityMatrix* aoDensityMatrix,
                                                                   const CGtoContainer*    gtoContainer,
                                                                   const double*           gridCoordinatesX,
                                                                   const double*           gridCoordinatesY,
                                                                   const double*           gridCoordinatesZ,
                                                                   const int32_t           gridOffset,
                                                                   const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForGGA(
                gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distUnrestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // compute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        gtorec::computeGtosValuesForGGA(
            gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distUnrestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_distRestrictedDensityValuesForLda(CDensityGrid*              densityGrid,
                                                       const CAODensityMatrix*    aoDensityMatrix,
                                                       const CMemBlock2D<double>& gtoValues,
                                                       const int32_t              gridOffset,
                                                       const int32_t              gridBlockPosition,
                                                       const int32_t              nGridPoints) const
{
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();

    for (int32_t i = 0; i < ndmat; i++)
    {
        // set up pointer to density grid data

        auto rhoa = densityGrid->alphaDensity(i);

        // set up poiinter to density matrix data

        auto denmat = aoDensityMatrix->alphaDensity(i);

        auto naos = aoDensityMatrix->getNumberOfRows(i);

        // loop over density matrix

        for (int32_t j = 0; j < naos; j++)
        {
            auto bgaos = gtoValues.data(j);

            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);

                auto dval = (j == k) ? denmat[j * naos + k] : denmat[j * naos + k] + denmat[k * naos + j];

                if (std::fabs(dval) > _thresholdOfDensity)
                {
                    #pragma omp simd
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        rhoa[gridOffset + gridBlockPosition + l] += dval * bgaos[l] * kgaos[l];
                    }
                }
            }
        }
    }
}

void
CDensityGridDriver::_distRestrictedDensityValuesForGga(CDensityGrid*              densityGrid,
                                                       const CAODensityMatrix*    aoDensityMatrix,
                                                       const CMemBlock2D<double>& gtoValues,
                                                       const CMemBlock2D<double>& gtoValuesX,
                                                       const CMemBlock2D<double>& gtoValuesY,
                                                       const CMemBlock2D<double>& gtoValuesZ,
                                                       const int32_t              gridOffset,
                                                       const int32_t              gridBlockPosition,
                                                       const int32_t              nGridPoints) const
{
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();

    for (int32_t i = 0; i < ndmat; i++)
    {
        // set up pointer to density grid data

        auto rhoa = densityGrid->alphaDensity(i);

        auto gradax = densityGrid->alphaDensityGradientX(i);

        auto graday = densityGrid->alphaDensityGradientY(i);

        auto gradaz = densityGrid->alphaDensityGradientZ(i);

        // set up poiinter to density matrix data

        auto denmat = aoDensityMatrix->alphaDensity(i);

        auto naos = aoDensityMatrix->getNumberOfRows(i);

        // loop over density matrix

        for (int32_t j = 0; j < naos; j++)
        {
            auto bgaos = gtoValues.data(j);

            auto bgaox = gtoValuesX.data(j);

            auto bgaoy = gtoValuesY.data(j);

            auto bgaoz = gtoValuesZ.data(j);

            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);

                auto kgaox = gtoValuesX.data(k);

                auto kgaoy = gtoValuesY.data(k);

                auto kgaoz = gtoValuesZ.data(k);

                auto dval = (j == k) ? denmat[j * naos + k] : denmat[j * naos + k] + denmat[k * naos + j];

                if (std::fabs(dval) > _thresholdOfDensity)
                {
                    #pragma omp simd
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        rhoa[gridOffset + gridBlockPosition + l] += dval * bgaos[l] * kgaos[l];

                        gradax[gridOffset + gridBlockPosition + l] += dval * (bgaox[l] * kgaos[l] + bgaos[l] * kgaox[l]);

                        graday[gridOffset + gridBlockPosition + l] += dval * (bgaoy[l] * kgaos[l] + bgaos[l] * kgaoy[l]);

                        gradaz[gridOffset + gridBlockPosition + l] += dval * (bgaoz[l] * kgaos[l] + bgaos[l] * kgaoz[l]);
                    }
                }
            }
        }
    }
}

void
CDensityGridDriver::_distRestrictedDensityValuesForMgga(CDensityGrid*              densityGrid,
                                                        const CAODensityMatrix*    aoDensityMatrix,
                                                        const CMemBlock2D<double>& gtoValues,
                                                        const CMemBlock2D<double>& gtoValuesX,
                                                        const CMemBlock2D<double>& gtoValuesY,
                                                        const CMemBlock2D<double>& gtoValuesZ,
                                                        const int32_t              gridOffset,
                                                        const int32_t              gridBlockPosition,
                                                        const int32_t              nGridPoints) const
{
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();

    for (int32_t i = 0; i < ndmat; i++)
    {
        // set up pointer to density grid data

        auto rhoa = densityGrid->alphaDensity(i);

        auto gradax = densityGrid->alphaDensityGradientX(i);

        auto graday = densityGrid->alphaDensityGradientY(i);

        auto gradaz = densityGrid->alphaDensityGradientZ(i);

        auto lapa = densityGrid->alphaDensitytau(i);

        // set up poiinter to density matrix data

        auto denmat = aoDensityMatrix->alphaDensity(i);

        auto naos = aoDensityMatrix->getNumberOfRows(i);

        // loop over density matrix

        for (int32_t j = 0; j < naos; j++)
        {
            auto bgaos = gtoValues.data(j);

            auto bgaox = gtoValuesX.data(j);

            auto bgaoy = gtoValuesY.data(j);

            auto bgaoz = gtoValuesZ.data(j);

            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);

                auto kgaox = gtoValuesX.data(k);

                auto kgaoy = gtoValuesY.data(k);

                auto kgaoz = gtoValuesZ.data(k);

                auto dval = (j == k) ? denmat[j * naos + k] : denmat[j * naos + k] + denmat[k * naos + j];

                if (std::fabs(dval) > _thresholdOfDensity)
                {
                    #pragma omp simd
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        rhoa[gridOffset + gridBlockPosition + l] += dval * bgaos[l] * kgaos[l];

                        gradax[gridOffset + gridBlockPosition + l] += dval * (bgaox[l] * kgaos[l] + bgaos[l] * kgaox[l]);

                        graday[gridOffset + gridBlockPosition + l] += dval * (bgaoy[l] * kgaos[l] + bgaos[l] * kgaoy[l]);

                        gradaz[gridOffset + gridBlockPosition + l] += dval * (bgaoz[l] * kgaos[l] + bgaos[l] * kgaoz[l]);

                        lapa[gridOffset + gridBlockPosition + l] += 0.5 * dval * (bgaox[l] * kgaox[l] + bgaoy[l] * kgaoy[l] + bgaoz[l] * kgaoz[l]);
                    }
                }
            }
        }
    }
}

void
CDensityGridDriver::_distUnrestrictedDensityValuesForGga(CDensityGrid*              densityGrid,
                                                         const CAODensityMatrix*    aoDensityMatrix,
                                                         const CMemBlock2D<double>& gtoValues,
                                                         const CMemBlock2D<double>& gtoValuesX,
                                                         const CMemBlock2D<double>& gtoValuesY,
                                                         const CMemBlock2D<double>& gtoValuesZ,
                                                         const int32_t              gridOffset,
                                                         const int32_t              gridBlockPosition,
                                                         const int32_t              nGridPoints) const
{
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();

    for (int32_t i = 0; i < ndmat; i++)
    {
        // set up pointer to density grid data

        auto rhoa = densityGrid->alphaDensity(i);

        auto gradax = densityGrid->alphaDensityGradientX(i);

        auto graday = densityGrid->alphaDensityGradientY(i);

        auto gradaz = densityGrid->alphaDensityGradientZ(i);

        auto rhob = densityGrid->betaDensity(i);

        auto gradbx = densityGrid->betaDensityGradientX(i);

        auto gradby = densityGrid->betaDensityGradientY(i);

        auto gradbz = densityGrid->betaDensityGradientZ(i);

        // set up poiinter to density matrix data

        auto denmata = aoDensityMatrix->alphaDensity(i);

        auto denmatb = aoDensityMatrix->betaDensity(i);

        auto naos = aoDensityMatrix->getNumberOfRows(i);

        // loop over density matrix

        for (int32_t j = 0; j < naos; j++)
        {
            auto bgaos = gtoValues.data(j);

            auto bgaox = gtoValuesX.data(j);

            auto bgaoy = gtoValuesY.data(j);

            auto bgaoz = gtoValuesZ.data(j);

            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);

                auto kgaox = gtoValuesX.data(k);

                auto kgaoy = gtoValuesY.data(k);

                auto kgaoz = gtoValuesZ.data(k);

                auto dvala = (j == k) ? denmata[j * naos + k] : denmata[j * naos + k] + denmata[k * naos + j];

                auto dvalb = (j == k) ? denmatb[j * naos + k] : denmatb[j * naos + k] + denmatb[k * naos + j];

                if (std::fabs(dvala + dvalb) > _thresholdOfDensity)
                {
                    #pragma omp simd
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        const double fgx = (bgaox[l] * kgaos[l] + bgaos[l] * kgaox[l]);

                        const double fgy = (bgaoy[l] * kgaos[l] + bgaos[l] * kgaoy[l]);

                        const double fgz = (bgaoz[l] * kgaos[l] + bgaos[l] * kgaoz[l]);

                        rhoa[gridOffset + gridBlockPosition + l] += dvala * bgaos[l] * kgaos[l];

                        gradax[gridOffset + gridBlockPosition + l] += dvala * fgx;

                        graday[gridOffset + gridBlockPosition + l] += dvala * fgy;

                        gradaz[gridOffset + gridBlockPosition + l] += dvala * fgz;

                        rhob[gridOffset + gridBlockPosition + l] += dvalb * bgaos[l] * kgaos[l];

                        gradbx[gridOffset + gridBlockPosition + l] += dvalb * fgx;

                        gradby[gridOffset + gridBlockPosition + l] += dvalb * fgy;

                        gradbz[gridOffset + gridBlockPosition + l] += dvalb * fgz;
                    }
                }
            }
        }
    }
}

int32_t
CDensityGridDriver::_getSizeOfBlock() const
{
    return 500;
}


CDensityGrid
CDensityGridDriver::pdft(const CAODensityMatrix& aoDensityMatrix,
                         double* twoDM,
                         double* activeMOs,
                         int nActive,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const xcfun             xcFunctional)
{
    // initialize density grid

    CDensityGrid densityGrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(),
                       xcFunctional, dengrid::ab);

    densityGrid.zero();

    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();
   
    auto mgz = molecularGrid.getCoordinatesZ();
   
    // set up pointer to density matrix
   
    auto denptr = &aoDensityMatrix;
   
    // set up poinet to density grid
   
    auto dgridptr = &densityGrid;
   
    // generate density on grid points
   
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                if (xcFunctional == xcfun::lda)
                {
                    #pragma omp task firstprivate(tbsize, tbposition)
                    {
                        _PDFT_Lda(dgridptr, denptr, twoDM, activeMOs, nActive, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                    }
                }
                else if (xcFunctional == xcfun::gga)
                {
                    #pragma omp task firstprivate(tbsize, tbposition)
                    {
                        _PDFT_Gga(dgridptr, denptr, twoDM, activeMOs, nActive, gtovec, mgx, mgy, mgz, tbposition, tbsize);
                    }
                }
                
            }
        }
    }
   
    //Why only GGA?
    if (xcFunctional == xcfun::gga)
    {
        // finalize density grid

        densityGrid.computeDensityNorms();
    }

    // destroy GTOs container

    delete gtovec;

    return densityGrid;
}

void
CDensityGridDriver::_PDFT_Lda(      CDensityGrid*     densityGrid,
                                    const CAODensityMatrix* aoDensityMatrix,
                                    double* twoDM,
                                    double* activeMOs,
                                    int nActive,
                                    const CGtoContainer*    gtoContainer,
                                    const double*           gridCoordinatesX,
                                    const double*           gridCoordinatesY,
                                    const double*           gridCoordinatesZ,
                                    const int32_t           gridOffset,
                                    const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);

            _distPDFT_LDA(densityGrid,twoDM, activeMOs, nActive, naos, gaos, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }
    // compute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);
        _distPDFT_LDA(densityGrid,twoDM, activeMOs, nActive, naos, gaos, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_distPDFT_LDA(      CDensityGrid*        densityGrid,
                                        double* twoDM,
                                        double* activeMOs,
                                        int nActive,
                                        int nAOs,
                                        const CMemBlock2D<double>& gtoValues,
                                        const int32_t              gridOffset,
                                        const int32_t              gridBlockPosition,
                                        const int32_t              nGridPoints) const
{
    int ndmat=1;

    for (int32_t i = 0; i < ndmat; i++)
    {   
        // set up pointer to density grid data
        
        auto rhoa = densityGrid->alphaDensity(i);

        auto rhob = densityGrid->betaDensity(i);

        // Compute MOs on grid

        double** ActMOs=new double*[nActive];

        for (int32_t iAct = 0; iAct < nActive; iAct++)
        {
            ActMOs[iAct] = new double[nGridPoints];
            #pragma omp simd 
            for (int32_t l = 0; l < nGridPoints; l++)
            {
                ActMOs[iAct][l] = 0.0;
            }

            auto MO = &activeMOs[iAct*nAOs];
            for (int32_t j = 0; j < nAOs; j++)
            {
                auto bgaos = gtoValues.data(j);
                #pragma omp simd 
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    ActMOs[iAct][l]+= bgaos[l] * MO[j];
                }
            }
        }

        // Compute on-top pair density on grid
        for (int32_t iAct = 0; iAct < nActive; iAct++)
        {
            auto iMO = ActMOs[iAct];
            for (int32_t jAct = 0; jAct < nActive; jAct++)
            {
                int ioff=iAct*nActive+jAct;
                auto jMO = ActMOs[jAct];
                for (int32_t kAct = 0; kAct < nActive; kAct++)
                {
                    int joff=ioff*nActive+kAct;
                    auto kMO = ActMOs[kAct];
                    for (int32_t lAct = 0; lAct < nActive; lAct++)
                    {
                        auto lMO = ActMOs[lAct];
                        double dval_b = twoDM[joff*nActive+lAct];
                        #pragma omp simd 
                        for (int32_t l = 0; l < nGridPoints; l++)
                        {
                            double fact=iMO[l] * jMO[l] * kMO[l] * lMO[l];
                            rhob[gridOffset + gridBlockPosition + l] += dval_b * fact;
                        }
                    }
                }
            }
        }
        // Compute the "effective" alpha and beta densities
        for (int32_t l = 0; l < nGridPoints; l++)
        {
            auto da = rhoa[gridOffset + gridBlockPosition + l];
            auto db = rhob[gridOffset + gridBlockPosition + l];

            auto delta=-db;
            if (delta>0)
            {
                delta=sqrt(2.0*delta);
            }
            else
            {
                delta=0.0;
            }
            rhoa[gridOffset + gridBlockPosition + l]=0.5*(da+delta);
            rhob[gridOffset + gridBlockPosition + l]=0.5*(da-delta);
        }
    }
}

void
CDensityGridDriver::_PDFT_Gga(      CDensityGrid*     densityGrid,
                                    const CAODensityMatrix* aoDensityMatrix,
                                    double* twoDM,
                                    double* activeMOs,
                                    int nActive,
                                    const CGtoContainer*    gtoContainer,
                                    const double*           gridCoordinatesX,
                                    const double*           gridCoordinatesY,
                                    const double*           gridCoordinatesZ,
                                    const int32_t           gridOffset,
                                    const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz,  gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                             gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
                                               gridOffset, igpnt, blockdim);
            _distPDFT_GGA(densityGrid,twoDM, activeMOs, nActive, naos, gaos, gaox, gaoy, gaoz, 
                                               gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }
    // compute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> gaox(blockdim, naos);

        CMemBlock2D<double> gaoy(blockdim, naos);

        CMemBlock2D<double> gaoz(blockdim, naos);

        gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                        gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
                                           gridOffset, igpnt, blockdim);

        _distPDFT_GGA(densityGrid,twoDM, activeMOs, nActive, naos, gaos, gaox, gaoy, gaoz,
                                           gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_distPDFT_GGA(      CDensityGrid*        densityGrid,
                                        double* twoDM,
                                        double* activeMOs,
                                        int nActive,
                                        int nAOs,
                                        const CMemBlock2D<double>& gtoValues,
                                        const CMemBlock2D<double>& gtoValuesX,
                                        const CMemBlock2D<double>& gtoValuesY,
                                        const CMemBlock2D<double>& gtoValuesZ,
                                        const int32_t              gridOffset,
                                        const int32_t              gridBlockPosition,
                                        const int32_t              nGridPoints) const
{
    int ndmat=1;
//#define FullTranslation

    for (int32_t i = 0; i < ndmat; i++)
    {   
        // set up pointer to density grid data
        
        auto rhoa = densityGrid->alphaDensity(i);

        auto gradax = densityGrid->alphaDensityGradientX(i);

        auto graday = densityGrid->alphaDensityGradientY(i);

        auto gradaz = densityGrid->alphaDensityGradientZ(i);

        auto rhob = densityGrid->betaDensity(i);

        auto gradbx = densityGrid->betaDensityGradientX(i);

        auto gradby = densityGrid->betaDensityGradientY(i);

        auto gradbz = densityGrid->betaDensityGradientZ(i);

        // Compute MOs on grid

        double** ActMOs=new double*[nActive];
#ifdef FullTranslation
        double** ActMOs_x=new double*[nActive];
        double** ActMOs_y=new double*[nActive];
        double** ActMOs_z=new double*[nActive];
#endif

        for (int32_t iAct = 0; iAct < nActive; iAct++)
        {
            ActMOs[iAct] = new double[nGridPoints];
#ifdef FullTranslation
            ActMOs_x[iAct] = new double[nGridPoints];
            ActMOs_y[iAct] = new double[nGridPoints];
            ActMOs_z[iAct] = new double[nGridPoints];
#endif
            #pragma omp simd 
            for (int32_t l = 0; l < nGridPoints; l++)
            {
                ActMOs[iAct][l] = 0.0;
#ifdef FullTranslation
                ActMOs_x[iAct][l] = 0.0;
                ActMOs_y[iAct][l] = 0.0;
                ActMOs_z[iAct][l] = 0.0;
#endif
            }

            auto MO = &activeMOs[iAct*nAOs];
            for (int32_t j = 0; j < nAOs; j++)
            {
                auto bgaos = gtoValues.data(j);
#ifdef FullTranslation
                auto bgaox = gtoValuesX.data(j);
                auto bgaoy = gtoValuesY.data(j);
                auto bgaoz = gtoValuesZ.data(j);
#endif
                #pragma omp simd 
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    ActMOs[iAct][l]+= bgaos[l] * MO[j];
#ifdef FullTranslation
                    ActMOs_x[iAct][l]+= bgaox[l] * MO[j];
                    ActMOs_y[iAct][l]+= bgaoy[l] * MO[j];
                    ActMOs_z[iAct][l]+= bgaoz[l] * MO[j];
#endif
                }
            }
        }

        // Compute on-top pair density on grid
        for (int32_t iAct = 0; iAct < nActive; iAct++)
        {
            auto iMO = ActMOs[iAct];
#ifdef FullTranslation
            auto iMO_x = ActMOs_x[iAct];
            auto iMO_y = ActMOs_y[iAct];
            auto iMO_z = ActMOs_z[iAct];
#endif
            for (int32_t jAct = 0; jAct < nActive; jAct++)
            {
                int ioff=iAct*nActive+jAct;
                auto jMO = ActMOs[jAct];
#ifdef FullTranslation
                auto jMO_x = ActMOs_x[jAct];
                auto jMO_y = ActMOs_y[jAct];
                auto jMO_z = ActMOs_z[jAct];
#endif
                for (int32_t kAct = 0; kAct < nActive; kAct++)
                {
                    int joff=ioff*nActive+kAct;
                    auto kMO = ActMOs[kAct];
#ifdef FullTranslation
                    auto kMO_x = ActMOs_x[kAct];
                    auto kMO_y = ActMOs_y[kAct];
                    auto kMO_z = ActMOs_z[kAct];
#endif
                    for (int32_t lAct = 0; lAct < nActive; lAct++)
                    {
                        auto lMO = ActMOs[lAct];
#ifdef FullTranslation
                        auto lMO_x = ActMOs_x[lAct];
                        auto lMO_y = ActMOs_y[lAct];
                        auto lMO_z = ActMOs_z[lAct];
#endif
                        double dval_b = twoDM[joff*nActive+lAct];
                        #pragma omp simd 
                        for (int32_t l = 0; l < nGridPoints; l++)
                        {
                            double fact=iMO[l] * jMO[l] * kMO[l] * lMO[l];
                            rhob[gridOffset + gridBlockPosition + l] += dval_b * fact;
#ifdef FullTranslation
                            double fgx= iMO_x[l] * jMO[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO_x[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO_x[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO[l] * lMO_x[l];
                            double fgy= iMO_y[l] * jMO[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO_y[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO_y[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO[l] * lMO_y[l];
                            double fgz= iMO_z[l] * jMO[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO_z[l] * kMO[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO_z[l] * lMO[l]
                                       +iMO[l] * jMO[l] * kMO[l] * lMO_z[l];
                            gradbx[gridOffset + gridBlockPosition + l] += dval_b * fgx;
                            gradby[gridOffset + gridBlockPosition + l] += dval_b * fgy;
                            gradbz[gridOffset + gridBlockPosition + l] += dval_b * fgz;
#endif
                        }
                    }
                }
            }
        }
        // Compute the "effective" alpha and beta densities
        for (int32_t l = 0; l < nGridPoints; l++)
        {
            auto da  = rhoa  [gridOffset + gridBlockPosition + l];
            auto dax = gradax[gridOffset + gridBlockPosition + l];
            auto day = graday[gridOffset + gridBlockPosition + l];
            auto daz = gradaz[gridOffset + gridBlockPosition + l];
            auto db = rhob   [gridOffset + gridBlockPosition + l];
#ifdef FullTranslation
            auto dbx = gradbx[gridOffset + gridBlockPosition + l];
            auto dby = gradby[gridOffset + gridBlockPosition + l];
            auto dbz = gradbz[gridOffset + gridBlockPosition + l];
#endif

            double delta=0.0;
            if (db<0)
            {
                delta=sqrt(-2.0*db);
            }
            rhoa[gridOffset + gridBlockPosition + l]=0.5*(da+delta);
            rhob[gridOffset + gridBlockPosition + l]=0.5*(da-delta);
#ifdef FullTranslation
//Correct formulas (except for potential bug)
            if (delta>1.0e-8)
            {
                gradax[gridOffset + gridBlockPosition + l]=0.5*(dax+dbx/delta);
                graday[gridOffset + gridBlockPosition + l]=0.5*(day+dby/delta);
                gradaz[gridOffset + gridBlockPosition + l]=0.5*(daz+dbz/delta);
                gradbx[gridOffset + gridBlockPosition + l]=0.5*(dax-dbx/delta);
                gradby[gridOffset + gridBlockPosition + l]=0.5*(day-dby/delta);
                gradbz[gridOffset + gridBlockPosition + l]=0.5*(daz-dbz/delta);
            }
#else
//"translated" formulas from Li Manni 2014
            if (da>1.0e-8)
            {
                gradax[gridOffset + gridBlockPosition + l]=0.5*(dax+delta*dax/da);
                graday[gridOffset + gridBlockPosition + l]=0.5*(day+delta*day/da);
                gradaz[gridOffset + gridBlockPosition + l]=0.5*(daz+delta*daz/da);
                gradbx[gridOffset + gridBlockPosition + l]=0.5*(dax-delta*dax/da);
                gradby[gridOffset + gridBlockPosition + l]=0.5*(day-delta*day/da);
                gradbz[gridOffset + gridBlockPosition + l]=0.5*(daz-delta*daz/da);
            }
#endif
            else
            {
                gradax[gridOffset + gridBlockPosition + l]=0.5*dax;
                graday[gridOffset + gridBlockPosition + l]=0.5*day;
                gradaz[gridOffset + gridBlockPosition + l]=0.5*daz;
                gradbx[gridOffset + gridBlockPosition + l]=0.5*dax;
                gradby[gridOffset + gridBlockPosition + l]=0.5*day;
                gradbz[gridOffset + gridBlockPosition + l]=0.5*daz;
            }
        }
    }
}
