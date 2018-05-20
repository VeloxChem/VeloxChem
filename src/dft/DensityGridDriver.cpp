//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DensityGridDriver.hpp"

#include "SystemClock.hpp"
#include "OMPTasks.hpp"

CDensityGridDriver::CDensityGridDriver(const int32_t  globRank,
                                       const int32_t  globNodes,
                                       const execmode runMode,
                                             MPI_Comm comm)
    :_globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)

    , _runMode(runMode)

    , _thresholdOfDensity(1.0e-13)

    , _thresholdOfPrimGTOs(1.0e-15)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CDensityGridDriver::~CDensityGridDriver()
{
    
}

void
CDensityGridDriver::generate(const CMolecule&       molecule,
                             const CMolecularBasis& basis,
                             const CMolecularGrid&  molGrid,
                             const xcfun            xcFunctional, 
                                   COutputStream&   oStream,
                                   MPI_Comm         comm)
{
    CSystemClock tim;
    
    // execution mode: CPU
    
    if (_runMode == execmode::cpu)
    {
        _genDensityGridOnCPU(molecule, basis, molGrid, xcFunctional,
                             oStream, comm);
    }
    
    // execution mode: CPU/GPU
    
    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }
    
    printf("Density Grid: rank %i points: %i time %lf s.\n",
           _locRank, molGrid.getNumberOfGridPoints(), tim.getElapsedTimeInSeconds());
    
    return;
}

void
CDensityGridDriver::_genDensityGridOnCPU(const CMolecule&       molecule,
                                         const CMolecularBasis& basis,
                                         const CMolecularGrid&  molGrid,
                                         const xcfun            xcFunctional,
                                               COutputStream&   oStream,
                                               MPI_Comm         comm)
{
    // set up OMP tasks
    
    COMPTasks omptaks(5);
    
    omptaks.set(molGrid.getNumberOfGridPoints());
    
    auto ntasks = omptaks.getNumberOfTasks();
    
    auto tbsizes = omptaks.getTaskSizes();
    
    auto tbpositions = omptaks.getTaskPositions();
    
    // set up molecular grid data
    
    auto mgx = molGrid.getCoordinatesX();
    
    auto mgy = molGrid.getCoordinatesY();
    
    auto mgz = molGrid.getCoordinatesZ();
    
    // create GTOs container
    
    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
    // generate density on grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz,\
                                gtovec)
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
                    _genBatchOfDensityGridPoints(gtovec, mgx, mgy, mgz,
                                                 tbposition, tbsize);
                }
            }
        }
    }
    
    // delete GTOs container
    
    delete gtovec; 
}

void
CDensityGridDriver::_genBatchOfDensityGridPoints(const CGtoContainer* gtoContainer,
                                                 const double*        gridCoordinatesX,
                                                 const double*        gridCoordinatesY,
                                                 const double*        gridCoordinatesZ,
                                                 const int32_t        gridOffset,
                                                 const int32_t        nGridPoints)
{
    // local copy of GTOs containers
    
    auto gtovec = CGtoContainer(*gtoContainer);
    
    auto cmpvec = gtovec;
    
    // allocate vector reduced indexes
    
    CMemBlock2D<int32_t> redidx(gtovec.getNumberOfGtoBlocks(), 2);
    
    // allocate screening data
    
    auto scrdata = gtovec.getPrimBuffer();
    
    // loop over batch of grid points
    
    for (int32_t i = 0; i < nGridPoints; i++)
    {
        // grid point coordinates
        
        auto gx = gridCoordinatesX[i];
        
        auto gy = gridCoordinatesY[i];

        auto gz = gridCoordinatesZ[i];
        
        // compute screening factors
        
        _compScreeningFactors(gtovec, gx, gy, gz, scrdata);
        
        // update screened GTOs container
        
        cmpvec.compress(gtovec, redidx, scrdata, _thresholdOfPrimGTOs);
    }
    
    printf("task: size: %i pos %i\n", nGridPoints, gridOffset);
}

void
CDensityGridDriver::_compScreeningFactors(const CGtoContainer&        gtoContainer,
                                          const double                gridCoordinateX,
                                          const double                gridCoordinateY,
                                          const double                gridCoordinateZ,
                                                CVecMemBlock<double>& screenFactors)
{
    for (int32_t i = 0; i < gtoContainer.getNumberOfGtoBlocks(); i++)
    {
        // set up primitives GTOs data
        
        auto bang = gtoContainer.getAngularMomentum(i);
        
        auto pexps = gtoContainer.getExponents(i);
        
        auto pfacts = gtoContainer.getNormFactors(i);
        
        // set up primitive GTOs coordinates
        
        auto coordsx = gtoContainer.getCoordinatesX(i);
        
        auto coordsy = gtoContainer.getCoordinatesY(i);
        
        auto coordsz = gtoContainer.getCoordinatesZ(i);
        
        // set up screening factors
        
        auto sfacts = screenFactors[i].data();
        
        // loop over GTOs in GTOs block
        
        for (int32_t j = 0; j < gtoContainer.getNumberOfPrimGtos(i); j++)
        {
            // compute distances
            
            auto rx = gridCoordinateX - coordsx[j];
            
            auto ry = gridCoordinateY - coordsy[j];
            
            auto rz = gridCoordinateZ - coordsz[j];
            
            // determine max distance
            
            auto rmax = (rx > ry) ? rx : ry;
            
            if (rz > rmax) rmax = rz;
            
            // compute zero order overlap
            
            sfacts[j] = pfacts[j] * std::exp(-pexps[j] * (rx * rx  + ry * ry +
                                                          
                                                          rz * rz));
            
            // compose screening factor
            
            sfacts[j] *= _getScaleFactor(bang, rmax);
        }
    }
}

double
CDensityGridDriver::_getScaleFactor(const int32_t angularMomentum,
                                    const double  maxRadius) const
{
    if (angularMomentum == 0) return 1.0;
    
    if (angularMomentum == 1) return maxRadius;
    
    auto r2max = maxRadius * maxRadius;
    
    if (angularMomentum == 2) return r2max;
    
    if (angularMomentum == 3) return r2max * maxRadius;
    
    auto r4max = r2max * r2max;
    
    if (angularMomentum == 4) return r4max;
    
    if (angularMomentum == 5) return r4max * maxRadius;
    
    if (angularMomentum == 6) return r4max * r2max;
    
    // TO DO: implement higher order if needed
    
    return 0.0;
}

