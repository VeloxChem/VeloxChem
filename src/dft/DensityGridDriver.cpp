//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGridDriver.hpp"

#include <cmath>

#include "GenFunc.hpp"
#include "OMPTasks.hpp"
#include "MpiFunc.hpp"
#include "AngularMomentum.hpp"
#include "GtoFunc.hpp"


CDensityGridDriver::CDensityGridDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);

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
CDensityGridDriver::_genDensityGridOnCPU(      CDensityGrid&     densityGrid,
                                         const CAODensityMatrix& density,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CMolecularGrid&   molecularGrid,
                                         const xcfun             xcFunctional)
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
    // set up pointer to density matrix
    
    auto denptr = &density;
    
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
                    _genBatchOfDensityGridPoints(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize, xcFunctional);
                }
            }
        }
    }

    // delete GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genBatchOfDensityGridPoints(      CDensityGrid*     densityGrid, 
                                                 const CAODensityMatrix* aoDensityMatrix,
                                                 const CGtoContainer*    gtoContainer,
                                                 const double*           gridCoordinatesX,
                                                 const double*           gridCoordinatesY,
                                                 const double*           gridCoordinatesZ,
                                                 const int32_t           gridOffset,
                                                 const int32_t           nGridPoints,
                                                 const xcfun             xcFunctional)
{
    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        for (int32_t j = i; j < gtovec.getNumberOfGtoBlocks(); j++)
        {
            _compDensityForGtoBlocks(densityGrid, aoDensityMatrix, bgtos, gtovec.getGtoBlock(j),
                                     gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                     gridOffset, nGridPoints, xcFunctional);
        }
    }
}

void
CDensityGridDriver::_compDensityForGtoBlocks(      CDensityGrid*     densityGrid,
                                             const CAODensityMatrix* aoDensityMatrix,
                                             const CGtoBlock&        braGtoBlock,
                                             const CGtoBlock&        ketGtoBlock,
                                             const double*           gridCoordinatesX,
                                             const double*           gridCoordinatesY,
                                             const double*           gridCoordinatesZ,
                                             const int32_t           gridOffset,
                                             const int32_t           nGridPoints,
                                             const xcfun             xcFunctional)
{
    // determine symmetry of bra and ket sides
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
    // set up max density matrix elements for (i,j) pair of density indexes
    
    auto ndmat = aoDensityMatrix->getNumberOfMatrices();
    
    // angular momentum data for bra and ket
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // set up Cartesian GTOs buffers
    
    auto nvcomp = xcfun_components(xcFunctional);
    
    auto bncart = angmom::to_CartesianComponents(bang);
    
    auto kncart = angmom::to_CartesianComponents(kang);
    
    auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();
    
    auto kcartbuff = (kang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * kncart) : CMemBlock2D<double>();
    
    // set up spherical GTOs buffers
    
    auto bnspher = angmom::to_SphericalComponents(bang);
    
    auto knspher = angmom::to_SphericalComponents(kang);
    
    CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);
    
    CMemBlock2D<double> kspherbuff(nGridPoints, nvcomp * knspher);
    
    // density data for pair of contracted GTOs
    
    CMemBlock2D<double> denpair(ndmat, angmom::to_SphericalComponents(bang, kang));
    
    // density type: restricted or unrestricted
    
    bool isrest = aoDensityMatrix->isRestricted();
    
    for (int32_t i = 0; i < braGtoBlock.getNumberOfContrGtos(); i++)
    {
        gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                       braGtoBlock, i, xcFunctional);
        
        for (int32_t j = 0; j < ketGtoBlock.getNumberOfContrGtos(); j++)
        {
            if (_setDensityPair(denpair, aoDensityMatrix, braGtoBlock, ketGtoBlock, symbk, i, j))
            {
                gtorec::computeGtoValuesOnGrid(kspherbuff, kcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                               ketGtoBlock, j, xcFunctional);
                
                
                _addGtosPairContribution(densityGrid, denpair, isrest, bspherbuff, kspherbuff, bnspher, knspher, gridOffset, xcFunctional);
            }
        }
    }
}

bool
CDensityGridDriver::_setDensityPair(      CMemBlock2D<double>& densityPairs,
                                    const CAODensityMatrix*    aoDensityMatrix,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const bool                 isBraEqualKet,
                                    const int32_t              iBraContrGto,
                                    const int32_t              iKetContrGto) const
{
    bool islarge = false;
    
    // set up angular momentum data
    
    auto bcomp = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());
    
    auto kcomp = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());
    
    // loop over list of density matrices
    
    for (int32_t i = 0; i < aoDensityMatrix->getNumberOfMatrices(); i++)
    {
        auto denmat = aoDensityMatrix->getDensity(i);
        
        auto nrows = aoDensityMatrix->getNumberOfRows(i);
        
        for (int32_t j = 0; j < bcomp; j++)
        {
            auto bidx = (braGtoBlock.getIdentifiers(j))[iBraContrGto];
            
            for (int32_t k = 0; k < kcomp; k++)
            {
                auto kidx = (ketGtoBlock.getIdentifiers(k))[iKetContrGto];
                
                auto dval = (isBraEqualKet) ? denmat[bidx * nrows + kidx] : denmat[bidx * nrows + kidx] + denmat[kidx * nrows + bidx];
                
                (densityPairs.data(j * kcomp + k))[i] = dval;
                
                if (std::fabs(dval) > _thresholdOfDensity) islarge = true;
            }
        }
    }
    
    return islarge;
}

void
CDensityGridDriver::_addGtosPairContribution(      CDensityGrid*        densityGrid,
                                             const CMemBlock2D<double>& densityPairs,
                                             const bool                 isRestrictedDensity, 
                                             const CMemBlock2D<double>& braGtoValues,
                                             const CMemBlock2D<double>& ketGtoValues,
                                             const int32_t              braComponents,
                                             const int32_t              ketComponents,
                                             const int32_t              gridOffset,
                                             const xcfun                xcFunctional) const
{
    // determine number of density matrices
    
    auto ndmat = densityPairs.size(0);
    
    // set up number of grid points
    
    auto ngpoints = braGtoValues.size(0);
    
    // local density approximation
    
    if (xcFunctional == xcfun::lda)
    {
        if (isRestrictedDensity)
        {
            // restricted densities
            
            for (int32_t i = 0; i < ndmat; i++)
            {
                auto rhoa = densityGrid->alphaDensity(i);
                
                auto rhob = densityGrid->betaDensity(i);
                
                // loop over density pair components
                
                for (int32_t j = 0; j < braComponents; j++)
                {
                    auto bgto = braGtoValues.data(j);
                    
                    for (int32_t k = 0; k < ketComponents; k++)
                    {
                        auto kgto = ketGtoValues.data(k);
                        
                        auto fden = (densityPairs.data(j * ketComponents + k))[i];
                        
                        if (std::fabs(fden) > _thresholdOfDensity)
                        {
                            #pragma omp simd
                            for (int32_t l = 0; l < ngpoints; l++)
                            {
                                double fact = fden * bgto[l] * kgto[l];
                                
                                rhoa[gridOffset + l] += fact;
                                
                                rhob[gridOffset + l] += fact;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // unrestricted densities
            
            // FIX ME: implement unrestricted case
        }
    }
    
    // general gradient approximation
    
    if (xcFunctional == xcfun::gga)
    {
        if (isRestrictedDensity)
        {
            // restricted densities
            
            for (int32_t i = 0; i < ndmat; i++)
            {
                auto rhoa = densityGrid->alphaDensity(i);
                
                auto rhob = densityGrid->betaDensity(i);
                
                auto grada = densityGrid->alphaDensityGradient(i);
                
                auto gradb = densityGrid->betaDensityGradient(i);
                
                auto gradab = densityGrid->mixedDensityGradient(i);
                
                // loop over density pair components
                
                for (int32_t j = 0; j < braComponents; j++)
                {
                    auto bgto = braGtoValues.data(4 * j);
                    
                    auto bgto_x = braGtoValues.data(4 * j + 1);
                    
                    auto bgto_y = braGtoValues.data(4 * j + 2);
                    
                    auto bgto_z = braGtoValues.data(4 * j + 3);
                    
                    for (int32_t k = 0; k < ketComponents; k++)
                    {
                        auto kgto = ketGtoValues.data(4 * k);
                        
                        auto kgto_x = ketGtoValues.data(4 * k + 1);
                        
                        auto kgto_y = ketGtoValues.data(4 * k + 2);
                        
                        auto kgto_z = ketGtoValues.data(4 * k + 3);
                        
                        auto fden = (densityPairs.data(j * ketComponents + k))[i];
                        
                        if (std::fabs(fden) > _thresholdOfDensity)
                        {
                            #pragma omp simd
                            for (int32_t l = 0; l < ngpoints; l++)
                            {
                                double f0_0 = fden * bgto[l] * kgto[l];
                                
                                double fx_0 = fden * (bgto_x[l] * kgto[l] + bgto[l] * kgto_x[l]);
                                
                                double fy_0 = fden * (bgto_y[l] * kgto[l] + bgto[l] * kgto_y[l]);
                                
                                double fz_0 = fden * (bgto_z[l] * kgto[l] + bgto[l] * kgto_z[l]);
                                
                                double fg_0 = std::sqrt(fx_0 * fx_0 + fy_0 * fy_0 + fz_0 * fz_0);
                                
                                double fm_0 = fg_0 * fg_0;
                                
                                rhoa[gridOffset + l] += f0_0;
                                
                                rhob[gridOffset + l] += f0_0;
                                
                                grada[gridOffset + l] += fg_0;
                                
                                gradb[gridOffset + l] += fg_0;
                                
                                gradab[gridOffset + l] += fm_0;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // unrestricted densities
            
            // FIX ME: implement unrestricted case
        }
    }
    
    // FIX ME: meta general gradient approximation
    
}
