//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "XCIntegrator.hpp"

#include "MpiFunc.hpp"
#include "DensityGridDriver.hpp"
#include "FunctionalParser.hpp"
#include "AOKohnShamMatrix.hpp"
#include "OMPTasks.hpp"
#include "AngularMomentum.hpp"
#include "GtoFunc.hpp"


CXCIntegrator::CXCIntegrator(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    mpi::duplicate(comm, &_locComm);
    
    _thresholdOfDensity = 1.0e-13;
}

CXCIntegrator::~CXCIntegrator()
{
    
}

CAOKohnShamMatrix
CXCIntegrator::integrate(const CAODensityMatrix& aoDensityMatrix,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel) const
{
    // integrator handles single AO density only
    
    if (aoDensityMatrix.getNumberOfMatrices() == 1)
    {
        // parse exchange-correlation functional data
    
        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
    
        // create GTOs container
    
        CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
        // generate reference density grid
    
        CDensityGridDriver dgdrv(_locComm);
    
        auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());
    
        if (aoDensityMatrix.isRestricted())
        {
            // molecular and density grids
        
            std::vector<CMolecularGrid> mgrids(1, molecularGrid);
        
            std::vector<CDensityGrid> dgrids(1, CDensityGrid(molecularGrid.getNumberOfGridPoints(), 1, fvxc.getFunctionalType(), dengrid::ab));
        
            // generate screened density and molecular grids
        
            refdengrid.setScreenedGrids(dgrids, mgrids, _thresholdOfDensity, fvxc.getFunctionalType());
        
            // allocate XC gradient grid
            
            CXCGradientGrid vxcgrid(mgrids[0].getNumberOfGridPoints(), dgrids[0].getDensityGridType(), fvxc.getFunctionalType());
            
            // compute exchange-correlation functional first derrivatives
            
            fvxc.compute(vxcgrid, dgrids[0]);
            
            // compute Kohn-Sham matrix
            
            auto nrow = aoDensityMatrix.getNumberOfRows(0);
            
            auto ncol = aoDensityMatrix.getNumberOfColumns(0);
            
            CAOKohnShamMatrix ksmat(nrow, ncol, true); 
            
            _compRestrictedContribution(ksmat, gtovec, vxcgrid, dgrids[0], mgrids[0], fvxc.getFunctionalType());
            
            return ksmat;
        }
        else
        {
            // perform integration for spin unrestricted densities
        
            // FIX ME: implement partitioning and screening for density grid
        }
  
        // delete GTOs container
    
        delete gtovec;
    }
    
    return CAOKohnShamMatrix(); 
}

void
CXCIntegrator::integrate(      CAOFockMatrix&    aoFockMatrix,
                         const CAODensityMatrix& rwDensityMatrix,
                         const CAODensityMatrix& gsDensityMatrix,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel) const
{
    // parse exchange-correlation functional data
    
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
    
    // create GTOs container
    
    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
    // generate reference density grid
    
    CDensityGridDriver dgdrv(_locComm);
    
    auto refdengrid = dgdrv.generate(rwDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());

    // set up number of density matrices
    
    auto ndmat = refdengrid.getNumberOfDensityMatrices();
    
    if (rwDensityMatrix.isRestricted())
    {
        // molecular and density grids
        
        std::vector<CMolecularGrid> mgrids(ndmat, molecularGrid);
        
        std::vector<CDensityGrid> dgrids(ndmat, CDensityGrid(molecularGrid.getNumberOfGridPoints(), 1, fvxc.getFunctionalType(), dengrid::ab));
        
        // generate screened density and molecular grids
        
        refdengrid.setScreenedGrids(dgrids, mgrids, _thresholdOfDensity, fvxc.getFunctionalType());
        
        for (int32_t i = 0; i < ndmat; i++)
        {
            // rescale AO Fock matrix for non-hybrid functionals
            
            if (!fvxc.isHybridFunctional()) aoFockMatrix.scale(2.0, i);
            
            // compute ground state density for compressed grid
            
            auto curdengrid = dgdrv.generate(gsDensityMatrix, molecule, basis, mgrids[i], fvxc.getFunctionalType());
            
            // compute gradient and hessian of exchange-correlation functional
            
            CXCGradientGrid vxcgrid(mgrids[i].getNumberOfGridPoints(), curdengrid.getDensityGridType(), fvxc.getFunctionalType());
            
            CXCHessianGrid vxc2grid(mgrids[i].getNumberOfGridPoints(), curdengrid.getDensityGridType(), fvxc.getFunctionalType());
            
            fvxc.compute(vxcgrid, curdengrid);
            
            fvxc.compute(vxc2grid, curdengrid);
            
            // compute Kohn-Sham matrix
            
            auto nrow = gsDensityMatrix.getNumberOfRows(0);
            
            auto ncol = gsDensityMatrix.getNumberOfColumns(0);
            
            CAOKohnShamMatrix ksmat(nrow, ncol, true);
            
            // compute linear response contribution
            
            _compRestrictedContribution(ksmat, gtovec, vxcgrid, vxc2grid, dgrids[i], curdengrid, mgrids[i], fvxc.getFunctionalType());
            
            aoFockMatrix.addOneElectronMatrix(ksmat.getReferenceToKohnSham(), i); 
        }
    }
    else
    {
        // perform integration for spin unrestricted densities
        
        // FIX ME: implement partitioning and screening for density grid
    }
}

void
CXCIntegrator::_compRestrictedContribution(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCGradientGrid&   xcGradientGrid,
                                           const CDensityGrid&      densityGrid,
                                           const CMolecularGrid&    molecularGrid,
                                           const xcfun              xcFunctional) const
{
    // initialize Kohn-Sham matrix to zero
    
    aoKohnShamMatrix.zero();
    
    // set up OMP tasks
    
    COMPTasks omptaks(5);
    
    omptaks.set(molecularGrid.getNumberOfGridPoints());
    
    auto ntasks = omptaks.getNumberOfTasks();
    
    auto tbsizes = omptaks.getTaskSizes();
    
    auto tbpositions = omptaks.getTaskPositions();
    
    // set up pointer to molecular grid weigths
    
    auto mgx = molecularGrid.getCoordinatesX();
    
    auto mgy = molecularGrid.getCoordinatesY();
    
    auto mgz = molecularGrid.getCoordinatesZ();
    
    auto mgw = molecularGrid.getWeights();
    
    // set up pointer to density matrix
    
    auto xcgridptr = &xcGradientGrid;
    
    // set up poinet to density grid
    
    auto dgridptr = &densityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // number of electrons
    
    double xcele = 0.0;
    
    // exchange-correlation energy
    
    double xcene = 0.0;
    
    // generate density on grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, dgridptr, xcgridptr, ksmatprt, xcele, xcene)
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
                    _compRestrictedVXCForBatchOfGridPoints(ksmatprt, gtoContainer, xcgridptr, dgridptr, mgx, mgy, mgz, mgw,
                                                           tbposition, tbsize, xcFunctional);
                    
                    _compRestrictedEnergyForBatchOfGridPoints(xcele, xcene, xcgridptr, dgridptr, mgw, tbposition, tbsize);
                }
            }
        }
    }
    
    aoKohnShamMatrix.setNumberOfElectrons(xcele);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(xcene); 
}


void
CXCIntegrator::_compRestrictedContribution(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCGradientGrid&   xcGradientGrid,
                                           const CXCHessianGrid&    xcHessianGrid,
                                           const CDensityGrid&      rwDensityGrid,
                                           const CDensityGrid&      gsDensityGrid,
                                           const CMolecularGrid&    molecularGrid,
                                           const xcfun              xcFunctional) const
{
    // initialize Kohn-Sham matrix to zero
    
    aoKohnShamMatrix.zero();
    
    // set up OMP tasks
    
    COMPTasks omptaks(5);
    
    omptaks.set(molecularGrid.getNumberOfGridPoints());
    
    auto ntasks = omptaks.getNumberOfTasks();
    
    auto tbsizes = omptaks.getTaskSizes();
    
    auto tbpositions = omptaks.getTaskPositions();
    
    // set up pointer to molecular grid weigths
    
    auto mgx = molecularGrid.getCoordinatesX();
    
    auto mgy = molecularGrid.getCoordinatesY();
    
    auto mgz = molecularGrid.getCoordinatesZ();
    
    auto mgw = molecularGrid.getWeights();
    
    // set up pointer to gradient and hessian grids
    
    auto xcgridptr = &xcGradientGrid;
    
    auto xcgrid2ptr = &xcHessianGrid;
    
    // set up poinet to density grids
    
    auto drwptr = &rwDensityGrid;
    
    auto dgsptr = &gsDensityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // number of electrons
    
    double xcele = 0.0;
    
    // exchange-correlation energy
    
    double xcene = 0.0;
    
    // generate density on grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, drwptr, dgsptr, xcgridptr, xcgrid2ptr, ksmatprt, xcele, xcene)
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
                    _compRestrictedVXCForBatchOfGridPoints(ksmatprt, gtoContainer, xcgridptr, xcgrid2ptr, drwptr, dgsptr, mgx, mgy, mgz, mgw,
                                                           tbposition, tbsize, xcFunctional);
                    
                    _compRestrictedEnergyForBatchOfGridPoints(xcele, xcene, xcgridptr, drwptr, mgw, tbposition, tbsize);
                }
            }
        }
    }
    
    aoKohnShamMatrix.setNumberOfElectrons(xcele);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(xcene);
}

void
CXCIntegrator::_compRestrictedVXCForBatchOfGridPoints(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                                      const CGtoContainer*     gtoContainer,
                                                      const CXCGradientGrid*   xcGradientGrid,
                                                      const CDensityGrid*      densityGrid,
                                                      const double*            gridCoordinatesX,
                                                      const double*            gridCoordinatesY,
                                                      const double*            gridCoordinatesZ,
                                                      const double*            gridWeights,
                                                      const int32_t            gridOffset,
                                                      const int32_t            nGridPoints,
                                                      const xcfun              xcFunctional) const
{
    // local copy of GTOs containers
    
    auto gtovec = CGtoContainer(*gtoContainer);
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        for (int32_t j = i; j < gtovec.getNumberOfGtoBlocks(); j++)
        {
            _compRestrictedVXCForGtoBlocks(aoKohnShamMatrix, bgtos, gtovec.getGtoBlock(j), xcGradientGrid, densityGrid,
                                           gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridWeights,
                                           gridOffset, nGridPoints, xcFunctional);
        }
    }
}

void
CXCIntegrator::_compRestrictedVXCForBatchOfGridPoints(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                                      const CGtoContainer*     gtoContainer,
                                                      const CXCGradientGrid*   xcGradientGrid,
                                                      const CXCHessianGrid*    xcHessianGrid,
                                                      const CDensityGrid*      rwDensityGrid,
                                                      const CDensityGrid*      gsDensityGrid,
                                                      const double*            gridCoordinatesX,
                                                      const double*            gridCoordinatesY,
                                                      const double*            gridCoordinatesZ,
                                                      const double*            gridWeights,
                                                      const int32_t            gridOffset,
                                                      const int32_t            nGridPoints,
                                                      const xcfun              xcFunctional) const
{
    // local copy of GTOs containers
    
    auto gtovec = CGtoContainer(*gtoContainer);
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        for (int32_t j = i; j < gtovec.getNumberOfGtoBlocks(); j++)
        {
            _compRestrictedVXCForGtoBlocks(aoKohnShamMatrix, bgtos, gtovec.getGtoBlock(j), xcGradientGrid, xcHessianGrid,
                                           rwDensityGrid, gsDensityGrid, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridWeights,
                                           gridOffset, nGridPoints, xcFunctional);
        }
    }
}

void
CXCIntegrator::_compRestrictedEnergyForBatchOfGridPoints(      double&          xcElectrons,
                                                               double&          xcEnergy,
                                                         const CXCGradientGrid* xcGradientGrid,
                                                         const CDensityGrid*    densityGrid,
                                                         const double*          gridWeights,
                                                         const int32_t          gridOffset,
                                                         const int32_t          nGridPoints) const
{
    double nelec = 0.0, xcene = 0.0;
    
    // set up pointers to density grid
    
    auto rhoa = densityGrid->alphaDensity(0);
    
    auto rhob = densityGrid->betaDensity(0);
    
    // set up pointer to exchange-correlation energy grid
    
    auto efunc = xcGradientGrid->xcFunctionalValues();
    
    #pragma omp simd
    for (int32_t i = 0; i < nGridPoints; i++)
    {
        nelec += gridWeights[gridOffset + i] * (rhoa[gridOffset + i] + rhob[gridOffset + i]);
        
        xcene += gridWeights[gridOffset + i] * efunc[gridOffset + i];
    }
    
    #pragma omp critical (accelec)
    {
        xcElectrons += nelec;
        
        xcEnergy += xcene;
    }
}

void
CXCIntegrator::_compRestrictedVXCForGtoBlocks(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                              const CGtoBlock&         braGtoBlock,
                                              const CGtoBlock&         ketGtoBlock,
                                              const CXCGradientGrid*   xcGradientGrid,
                                              const CDensityGrid*      densityGrid,
                                              const double*            gridCoordinatesX,
                                              const double*            gridCoordinatesY,
                                              const double*            gridCoordinatesZ,
                                              const double*            gridWeights,
                                              const int32_t            gridOffset,
                                              const int32_t            nGridPoints,
                                              const xcfun              xcFunctional) const
{
    // determine symmetry of bra and ket sides
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
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
    
    // initialize GTOs pairs values
    
    CMemBlock<double> ppvalues(bnspher * knspher);
    
    for (int32_t i = 0; i < braGtoBlock.getNumberOfContrGtos(); i++)
    {
        gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                       braGtoBlock, i, xcFunctional);
        
        for (int32_t j = 0; j < ketGtoBlock.getNumberOfContrGtos(); j++)
        {
            gtorec::computeGtoValuesOnGrid(kspherbuff, kcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                           ketGtoBlock, j, xcFunctional);
            
            _compRestrictedVXCValueForGtosPair(ppvalues, bspherbuff, kspherbuff, bnspher, knspher, xcGradientGrid, densityGrid,
                                               gridWeights, gridOffset, xcFunctional);
            
            #pragma omp critical (ksmacc)
            _distRestrictedVXCValues(aoKohnShamMatrix, ppvalues, braGtoBlock, ketGtoBlock, symbk, i, j);
        }
    }
}

void
CXCIntegrator::_compRestrictedVXCForGtoBlocks(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                              const CGtoBlock&         braGtoBlock,
                                              const CGtoBlock&         ketGtoBlock,
                                              const CXCGradientGrid*   xcGradientGrid,
                                              const CXCHessianGrid*    xcHessianGrid,
                                              const CDensityGrid*      rwDensityGrid,
                                              const CDensityGrid*      gsDensityGrid,
                                              const double*            gridCoordinatesX,
                                              const double*            gridCoordinatesY,
                                              const double*            gridCoordinatesZ,
                                              const double*            gridWeights,
                                              const int32_t            gridOffset,
                                              const int32_t            nGridPoints,
                                              const xcfun              xcFunctional) const
{
    // determine symmetry of bra and ket sides
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
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
    
    // initialize GTOs pairs values
    
    CMemBlock<double> ppvalues(bnspher * knspher);
    
    for (int32_t i = 0; i < braGtoBlock.getNumberOfContrGtos(); i++)
    {
        gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                       braGtoBlock, i, xcFunctional);
        
        for (int32_t j = 0; j < ketGtoBlock.getNumberOfContrGtos(); j++)
        {
            gtorec::computeGtoValuesOnGrid(kspherbuff, kcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                           ketGtoBlock, j, xcFunctional);
            
            _compRestrictedVXCValueForGtosPair(ppvalues, bspherbuff, kspherbuff, bnspher, knspher, xcGradientGrid, xcHessianGrid,
                                               rwDensityGrid, gsDensityGrid, gridWeights, gridOffset, xcFunctional);
            
            #pragma omp critical (ksmacc)
            _distRestrictedVXCValues(aoKohnShamMatrix, ppvalues, braGtoBlock, ketGtoBlock, symbk, i, j);
        }
    }
}

void
CXCIntegrator::_compRestrictedVXCValueForGtosPair(      CMemBlock<double>&   pairValues,
                                                  const CMemBlock2D<double>& braGtoGridBuffer,
                                                  const CMemBlock2D<double>& ketGtoGridBuffer,
                                                  const int32_t              braAngularComponents,
                                                  const int32_t              ketAngularComponents,
                                                  const CXCGradientGrid*     xcGradientGrid,
                                                  const CDensityGrid*        densityGrid,
                                                  const double*              gridWeights, 
                                                  const int32_t              gridOffset,
                                                  const xcfun                xcFunctional) const
{
    // initialize pair values to zero
    
    pairValues.zero();
    
    // set up pointer to pair values
    
    auto ppvals = pairValues.data();
    
    // local density approximation
    
    if (xcFunctional == xcfun::lda)
    {
        auto ngpoints = braGtoGridBuffer.size(0);
        
        // set up pointers to gradient data
        
        auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);
        
        // NOTE: we compute F_a matrix, since F_a = F_b
        
        for (int32_t i = 0; i < braAngularComponents; i++)
        {
            auto bgto = braGtoGridBuffer.data(i);
            
            for (int32_t j = 0; j < ketAngularComponents; j++)
            {
                auto kgto = ketGtoGridBuffer.data(j); 
                
                double psum = 0.0;
                
                for (int32_t k = 0; k < ngpoints; k++)
                {
                    psum += gridWeights[gridOffset + k] * bgto[k] * kgto[k] * grhoa[gridOffset + k];
                }
                
                ppvals[i * ketAngularComponents + j] = psum;
            }
        }
        
        return;
    }
    
     // general gradient approximation
    
    if (xcFunctional == xcfun::gga)
    {
        auto ngpoints = braGtoGridBuffer.size(0);
        
        // set up pointers to gradient data
        
        auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);
        
        auto ggrada = xcGradientGrid->xcGradientValues(xcvars::grada);

        auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
        
        // set up pointers to density gradient norms
        
        auto ngrada = densityGrid->alphaDensityGradient(0);
        
        auto grada_x = densityGrid->alphaDensityGradientX(0);
        
        auto grada_y = densityGrid->alphaDensityGradientY(0);
        
        auto grada_z = densityGrid->alphaDensityGradientZ(0);
        
        // NOTE: we compute F_a matrix, since F_a = F_b
        
        for (int32_t i = 0; i < braAngularComponents; i++)
        {
            auto bgto = braGtoGridBuffer.data(4 * i);
            
            auto bgto_x = braGtoGridBuffer.data(4 * i + 1);
            
            auto bgto_y = braGtoGridBuffer.data(4 * i + 2);
            
            auto bgto_z = braGtoGridBuffer.data(4 * i + 3);
            
            for (int32_t j = 0; j < ketAngularComponents; j++)
            {
                auto kgto = ketGtoGridBuffer.data(4 * j);
                
                auto kgto_x = ketGtoGridBuffer.data(4 * j + 1);
                
                auto kgto_y = ketGtoGridBuffer.data(4 * j + 2);
                
                auto kgto_z = ketGtoGridBuffer.data(4 * j + 3);
                
                double psum = 0.0;
                
                for (int32_t k = 0; k < ngpoints; k++)
                {
                    double w = gridWeights[gridOffset + k];

                    double gx = grada_x[gridOffset + k];

                    double gy = grada_y[gridOffset + k];

                    double gz = grada_z[gridOffset + k];
		    
                    psum += w * bgto[k] * kgto[k] * grhoa[gridOffset + k];
                    
                    double fgrd = w * (ggrada[gridOffset + k] / ngrada[gridOffset + k] + ggradab[gridOffset + k]);
  
                    psum += fgrd * bgto[k] * (gx * kgto_x[k] + gy * kgto_y[k] + gz * kgto_z[k]);
                    
                    psum += fgrd * kgto[k] * (gx * bgto_x[k] + gy * bgto_y[k] + gz * bgto_z[k]);
                }
                
                ppvals[i * ketAngularComponents + j] = psum;
            }
        }
        
        return;
    }
    
    // FIX ME: impelemnt MGGA case
}

void
CXCIntegrator::_compRestrictedVXCValueForGtosPair(      CMemBlock<double>&   pairValues,
                                                  const CMemBlock2D<double>& braGtoGridBuffer,
                                                  const CMemBlock2D<double>& ketGtoGridBuffer,
                                                  const int32_t              braAngularComponents,
                                                  const int32_t              ketAngularComponents,
                                                  const CXCGradientGrid*     xcGradientGrid,
                                                  const CXCHessianGrid*      xcHessianGrid,
                                                  const CDensityGrid*        rwDensityGrid,
                                                  const CDensityGrid*        gsDensityGrid,
                                                  const double*              gridWeights,
                                                  const int32_t              gridOffset,
                                                  const xcfun                xcFunctional) const
{
    // initialize pair values to zero
    
    pairValues.zero();
    
    // set up pointer to pair values
    
    auto ppvals = pairValues.data();
    
    // local density approximation
    
    if (xcFunctional == xcfun::lda)
    {
        auto ngpoints = braGtoGridBuffer.size(0);
        
        // set up pointers to gradient data
        
        auto grho_aa = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        // set up pointer to perturbed density
        
        auto rhowa = rwDensityGrid->alphaDensity(0);
        
        // NOTE: we compute F_a matrix, since F_a = F_b
        
        for (int32_t i = 0; i < braAngularComponents; i++)
        {
            auto bgto = braGtoGridBuffer.data(i);
            
            for (int32_t j = 0; j < ketAngularComponents; j++)
            {
                auto kgto = ketGtoGridBuffer.data(j);
                
                double psum = 0.0;
                
                for (int32_t k = 0; k < ngpoints; k++)
                {
                    psum += gridWeights[gridOffset + k] * bgto[k] * kgto[k] * grho_aa[gridOffset + k] * rhowa[gridOffset + k];
                }
                
                ppvals[i * ketAngularComponents + j] = psum;
            }
        }
        
        return;
    }
    
    // general gradient approximation
    
    if (xcFunctional == xcfun::gga)
    {
      
        
        return;
    }
    
    // FIX ME: impelemnt MGGA case
}

void
CXCIntegrator::_distRestrictedVXCValues(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                        const CMemBlock<double>& pairValues,
                                        const CGtoBlock&         braGtoBlock,
                                        const CGtoBlock&         ketGtoBlock,
                                        const bool               isBraEqualKet,
                                        const int32_t            iBraContrGto,
                                        const int32_t            iKetContrGto) const
{
    // set up angular components
    
    auto bcomps = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());
    
    auto kcomps = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());
    
    // set up Kohn-Sham matrix data
    
    auto ksmat = aoKohnShamMatrix->getKohnSham();
    
    auto nrows = aoKohnShamMatrix->getNumberOfRows();
    
    // distribute VXC values for GTOs pair
    
    for (int32_t i = 0; i < bcomps; i++)
    {
        auto bidx = (braGtoBlock.getIdentifiers(i))[iBraContrGto];
        
        for (int32_t j = 0; j < kcomps; j++)
        {
            auto kidx = (ketGtoBlock.getIdentifiers(j))[iKetContrGto];
            
            auto fval = pairValues.at(i * kcomps + j);
            
            if (isBraEqualKet)
            {
                ksmat[bidx * nrows + kidx] += fval;
            }
            else
            {
                ksmat[bidx * nrows + kidx] += fval;
                
                ksmat[kidx * nrows + bidx] += fval;
            }
        }
    }
    
}
