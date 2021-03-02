//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "DensityGridDriver.hpp"

#include <cmath>

#include <mpi.h>

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
    
    CDensityGrid dgrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(),
                       xcFunctional, dengrid::ab);
    
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
CDensityGridDriver::_genDensityGridOnCPU(      CDensityGrid&     densityGrid,
                                         const CAODensityMatrix& aoDensityMatrix,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CMolecularGrid&   molecularGrid,
                                         const xcfun             xcFunctional)
{
    if (aoDensityMatrix.isRestricted())
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
    }
    
    if (aoDensityMatrix.isUnrestricted())
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
CDensityGridDriver::_genUnrestrictedDensityForLda(      CDensityGrid&     densityGrid,
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
CDensityGridDriver::_genBatchOfUnrestrictedDensityGridPointsForLda(      CDensityGrid*     densityGrid,
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
CDensityGridDriver::_distUnrestrictedDensityValuesForLda(      CDensityGrid*        densityGrid,
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
CDensityGridDriver::_genRestrictedDensityForLda(      CDensityGrid&     densityGrid,
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
                    _genBatchOfRestrictedDensityGridPointsForLda(dgridptr, denptr, gtovec, mgx, mgy, mgz,
                                                                 tbposition, tbsize);
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
CDensityGridDriver::_genRestrictedDensityForGga(      CDensityGrid&     densityGrid,
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
                    _genBatchOfRestrictedDensityGridPointsForGga(dgridptr, denptr, gtovec, mgx, mgy, mgz,
                                                                 tbposition, tbsize);
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
CDensityGridDriver::_genUnrestrictedDensityForGga(      CDensityGrid&     densityGrid,
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
                    _genBatchOfUnrestrictedDensityGridPointsForGga(dgridptr, denptr, gtovec, mgx, mgy, mgz,
                                                                   tbposition, tbsize);
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
CDensityGridDriver::_genBatchOfRestrictedDensityGridPointsForLda(      CDensityGrid*     densityGrid,
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
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                            gridOffset, igpnt, blockdim);
            
            _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim);
            
            igpnt += blockdim;
        }
    }
    
    // comopute remaining grid points block
    
    blockdim = nGridPoints % blockdim;
    
    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);
        
        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                        gridOffset, igpnt, blockdim);
        
        _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, gridOffset, igpnt, blockdim); 
    }
}

void
CDensityGridDriver::_genBatchOfRestrictedDensityGridPointsForGga(      CDensityGrid*     densityGrid,
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
            gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz,  gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
                                               gridOffset, igpnt, blockdim);
            
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
        
        gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                        gridCoordinatesZ, gridOffset, igpnt, blockdim);
        
        _distRestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
                                           gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_genBatchOfUnrestrictedDensityGridPointsForGga(      CDensityGrid*     densityGrid,
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
            gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz,  gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distUnrestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
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
        
        _distUnrestrictedDensityValuesForGga(densityGrid, aoDensityMatrix, gaos, gaox, gaoy, gaoz,
                                             gridOffset, igpnt, blockdim);
    }
}

void
CDensityGridDriver::_distRestrictedDensityValuesForLda(      CDensityGrid*        densityGrid,
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
CDensityGridDriver::_distRestrictedDensityValuesForGga(      CDensityGrid*        densityGrid,
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
CDensityGridDriver::_distUnrestrictedDensityValuesForGga(      CDensityGrid*        densityGrid,
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
                
                if (std::fabs(dvala+dvalb) > _thresholdOfDensity)
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
