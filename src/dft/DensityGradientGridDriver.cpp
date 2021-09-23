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

#include "DensityGradientGridDriver.hpp"

#include "MpiFunc.hpp"

#include "OMPTasks.hpp"
#include "GtoFunc.hpp"

CDensityGradientGridDriver::CDensityGradientGridDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CDensityGradientGridDriver::~CDensityGradientGridDriver()
{
}

CDensityGrid
CDensityGradientGridDriver::generate(const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid,
                                     const int32_t           iAtom)
{
    // initialize density grid
    
    CDensityGrid dgrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(),
                       3, dengrid::ab);
    
    dgrid.zero();
    
    // currently only LDA implemented
    
    _genRestrictedDensityForLda(dgrid, aoDensityMatrix, molecule, basis, molecularGrid, iAtom);
    
    return dgrid;
}


void
CDensityGradientGridDriver::_genRestrictedDensityForLda(      CDensityGrid&     densityGrid,
                                                        const CAODensityMatrix& aoDensityMatrix,
                                                        const CMolecule&        molecule,
                                                        const CMolecularBasis&  basis,
                                                        const CMolecularGrid&   molecularGrid,
                                                        const int32_t           iAtom) const
{
    // set up OMP tasks
    
    COMPTasks omptaks(5);
    
    omptaks.set(molecularGrid.getNumberOfGridPoints());
    
    auto ntasks = omptaks.getNumberOfTasks();
    
    auto tbsizes = omptaks.getTaskSizes();
    
    auto tbpositions = omptaks.getTaskPositions();
    
    // create GTOs container
    
    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
    CGtoContainer* xgtovec = new CGtoContainer(molecule, basis, iAtom, 1);
    
    // set up molecular grid data
    
    auto mgx = molecularGrid.getCoordinatesX();
    
    auto mgy = molecularGrid.getCoordinatesY();
    
    auto mgz = molecularGrid.getCoordinatesZ();
    
    // set up pointer to density matrix
    
    auto denptr = &aoDensityMatrix;
    
    // set up poinet to density grid
    
    auto dgridptr = &densityGrid;
    
    // generate density on grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, xgtovec, denptr, dgridptr)
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
                    _genBatchOfRestrictedDensityGridPointsForLda(dgridptr, denptr, gtovec, xgtovec, mgx, mgy, mgz,
                                                                 tbposition, tbsize);
                }
            }
        }
    }
    
    // finalize density grid
    
    densityGrid.updateBetaDensities();
    
    // destroy GTOs container
    
    delete gtovec;
    
    delete xgtovec;
}

void
CDensityGradientGridDriver::_genBatchOfRestrictedDensityGridPointsForLda(      CDensityGrid*     densityGrid,
                                                                         const CAODensityMatrix* aoDensityMatrix,
                                                                         const CGtoContainer*    gtoContainer,
                                                                         const CGtoContainer*    atmGtoContainer,
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
        
        CMemBlock2D<double> xgaos(blockdim, naos);
        
        CMemBlock2D<double> xgaox(blockdim, naos);
        
        CMemBlock2D<double> xgaoy(blockdim, naos);
        
        
        CMemBlock2D<double> xgaoz(blockdim, naos);
        
        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                            gridOffset, igpnt, blockdim);
            
            xgaos.zero();
            
            xgaox.zero();
            
            xgaoy.zero();
            
            xgaoz.zero();
            
            gtorec::computeGtosValuesForGGA(xgaos, xgaox, xgaoy, xgaoz, atmGtoContainer, gridCoordinatesX,
                                            gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
             _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, xgaox, xgaoy, xgaoz, gridOffset, igpnt, blockdim);
            
            igpnt += blockdim;
        }
    }
    
    // comopute remaining grid points block
    
    blockdim = nGridPoints % blockdim;
    
    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);
        
        CMemBlock2D<double> xgaos(blockdim, naos);
        
        xgaos.zero();
        
        CMemBlock2D<double> xgaox(blockdim, naos);
        
        xgaox.zero();
        
        CMemBlock2D<double> xgaoy(blockdim, naos);
        
        xgaoy.zero();
        
        CMemBlock2D<double> xgaoz(blockdim, naos);
        
        xgaoz.zero();
        
        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                        gridOffset, igpnt, blockdim);
        
        gtorec::computeGtosValuesForGGA(xgaos, xgaox, xgaoy, xgaoz, atmGtoContainer, gridCoordinatesX,
                                        gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);
        
        _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, xgaox, xgaoy, xgaoz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGradientGridDriver::_distRestrictedDensityValuesForLda(      CDensityGrid*        densityGrid,
                                                               const CAODensityMatrix*    aoDensityMatrix,
                                                               const CMemBlock2D<double>& braGtoValues,
                                                               const CMemBlock2D<double>& ketGtoValuesX,
                                                               const CMemBlock2D<double>& ketGtoValuesY,
                                                               const CMemBlock2D<double>& ketGtoValuesZ,
                                                               const int32_t              gridOffset,
                                                               const int32_t              gridBlockPosition,
                                                               const int32_t              nGridPoints) const
{
    if (aoDensityMatrix->getNumberOfDensityMatrices() == 1)
    {
        // set up pointer to density grid data
        
        auto rhoax = densityGrid->getComponent(0);
        
        auto rhoay = densityGrid->getComponent(1);
        
        auto rhoaz = densityGrid->getComponent(2);
        
        // set up pointer to density matrix data
        
        auto denmat = aoDensityMatrix->alphaDensity(0);
        
        auto naos = aoDensityMatrix->getNumberOfRows(0);
        
        // loop over density matrix
        
        for (int32_t i = 0; i < naos; i++)
        {
            auto bgaos = braGtoValues.data(i);
            
            for (int32_t j = 0; j < naos; j++)
            {
                auto kgaox = ketGtoValuesX.data(j);
                
                auto kgaoy = ketGtoValuesY.data(j);
                
                auto kgaoz = ketGtoValuesZ.data(j);
                
                const auto dval = 2.0 * denmat[i * naos + j];
                
                #pragma omp simd
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    rhoax[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaox[k];
                    
                    rhoay[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoy[k];
                    
                    rhoaz[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoz[k];
                }
            }
        }
    }
}

int32_t
CDensityGradientGridDriver::_getSizeOfBlock() const
{
    return 500;
}

