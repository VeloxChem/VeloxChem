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
#include "DenseLinearAlgebra.hpp"


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
        
        // temporary test for matrix driven approach to computation Kohn-Sham matrix 
        
        //if (fvxc.getFunctionalType() == xcfun::lda)
        //{
        //    return integrate_m3(aoDensityMatrix, molecule, basis, molecularGrid, xcFuncLabel);
        //}
    
        // create GTOs container
    
        CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
        // generate reference density grid
    
        CDensityGridDriver dgdrv(_locComm);
    
        auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());
    
        if (aoDensityMatrix.isRestricted())
        {
            // generate screened molecular and density grids
        
            CMolecularGrid mgrid(molecularGrid);
        
            CDensityGrid dgrid;
        
            refdengrid.getScreenedGridsPair(dgrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());
        
            // allocate XC gradient grid
            
            CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dgrid.getDensityGridType(), fvxc.getFunctionalType());
            
            // compute exchange-correlation functional first derrivatives
            
            fvxc.compute(vxcgrid, dgrid);
            
            // compute Kohn-Sham matrix
            
            auto nrow = aoDensityMatrix.getNumberOfRows(0);
            
            auto ncol = aoDensityMatrix.getNumberOfColumns(0);
            
            CAOKohnShamMatrix ksmat(nrow, ncol, true); 
            
            _compRestrictedContribution(ksmat, gtovec, vxcgrid, dgrid, mgrid, fvxc.getFunctionalType());
            
            // delete GTOs container
            
            delete gtovec;
            
            // done with Kohn-Sham matrix
            
            return ksmat;
        }
        else
        {
            // perform integration for spin unrestricted densities
        
            // FIX ME: implement partitioning and screening for density grid
        }
    }
    
    return CAOKohnShamMatrix(); 
}

CAOKohnShamMatrix
CXCIntegrator::integrate_m3(const CAODensityMatrix& aoDensityMatrix,
                            const CMolecule&        molecule,
                            const CMolecularBasis&  basis,
                            const CMolecularGrid&   molecularGrid,
                            const std::string&      xcFuncLabel) const
{
    CAOKohnShamMatrix ksmat;
    
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
            // generate screened molecular and density grids
            
            CMolecularGrid mgrid(molecularGrid);
            
            CDensityGrid dgrid;
            
            refdengrid.getScreenedGridsPair(dgrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());
            
            // allocate XC gradient grid
            
            CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dgrid.getDensityGridType(), fvxc.getFunctionalType());
            
            // compute exchange-correlation functional first derrivatives
            
            fvxc.compute(vxcgrid, dgrid);
            
            // compute Kohn-Sham matrix
            
            ksmat = CAOKohnShamMatrix(aoDensityMatrix.getNumberOfRows(0), aoDensityMatrix.getNumberOfColumns(0), true);
            
            _compRestrictedContributionM3(ksmat, gtovec, vxcgrid, dgrid, mgrid, fvxc.getFunctionalType());
        }
        
        // delete GTOs container
        
        delete gtovec;
    }
    
    return ksmat;
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
    
    auto refdengrid = dgdrv.generate(gsDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());
    
    // generate screened density and molecular grids
    
    CMolecularGrid mgrid(molecularGrid);
    
    CDensityGrid gsdengrid;
    
    refdengrid.getScreenedGridsPair(gsdengrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());

    // set up number of density matrices
    
    auto ndmat = rwDensityMatrix.getNumberOfDensityMatrices();
    
    if (rwDensityMatrix.isRestricted())
    {
        for (int32_t i = 0; i < ndmat; i++)
        {
            // compute gradient and hessian of exchange-correlation functional
            
            CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());
            
            CXCHessianGrid vxc2grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());
            
            fvxc.compute(vxcgrid, gsdengrid);
            
            fvxc.compute(vxc2grid, gsdengrid);
            
            // compute perturbed density grids (we will need to refactor this...)
            
            CAODensityMatrix currden({rwDensityMatrix.getReferenceToDensity(i)}, rwDensityMatrix.getDensityType());
            
            auto rwdengrid = dgdrv.generate(currden, molecule, basis, mgrid, fvxc.getFunctionalType());
            
            // compute perturbed Kohn-Sham matrix
            
            auto nrow = gsDensityMatrix.getNumberOfRows(0);
            
            auto ncol = gsDensityMatrix.getNumberOfColumns(0);
            
            CAOKohnShamMatrix ksmat(nrow, ncol, true);
            
            // compute linear response contribution
            
            _compRestrictedContribution(ksmat, gtovec, vxcgrid, vxc2grid, rwdengrid, gsdengrid, mgrid,
                                        fvxc.getFunctionalType());
            
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
    
    // initialize number of electrons and XC energy
    
    double xcele = 0.0, xcene = 0.0;
    
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
    
    // set number of electrons and XC energy
    
    aoKohnShamMatrix.setNumberOfElectrons(xcele);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(xcene); 
}

void
CXCIntegrator::_compRestrictedContributionM3(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                             const CGtoContainer*     gtoContainer,
                                             const CXCGradientGrid&   xcGradientGrid,
                                             const CDensityGrid&      densityGrid,
                                             const CMolecularGrid&    molecularGrid,
                                             const xcfun              xcFunctional) const
{
    // initialize Kohn-Sham matrix to zero
    
    aoKohnShamMatrix.zero();
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmat = aoKohnShamMatrix.getKohnSham();
    
    // set number of rows in grid block matrix
    
    auto nrows = _getNumberOfGridRows();
    
    // determine number of AO grid blocks
    
    auto nblocks = molecularGrid.getNumberOfGridPoints() / nrows;

    // set up current grid point
    
    int32_t igpnt = 0;
    
    // testing overlap  value
    
    auto gw = molecularGrid.getWeights();
    
    double fovl = 0.0;
    
    // loop over grid points blocks
    
    if (nblocks > 0)
    {
        // allocate AOs grid matrices for bra and ket
        
        CDenseMatrix bmat(nrows, gtoContainer->getNumberOfAtomicOrbitals());
        
        CDenseMatrix kmat(nrows, gtoContainer->getNumberOfAtomicOrbitals());
        
        for (int32_t i = 0; i < nblocks; i++)
        {
            _compGtosMatrixForLDA(bmat, gtoContainer, molecularGrid, igpnt, nrows);
            
            auto gao = bmat.values();
            
            for (int32_t j = 0; j < nrows; j++)
            {
                fovl += gw[igpnt + j] * gao[j] * gao[j];
            }
            
            _compRestrictedVXCMatrixForLDA(kmat, bmat, xcGradientGrid, molecularGrid.getWeights(), igpnt, nrows);
            
            denblas::multAtB(ksmat, 1.0, bmat, kmat);
            
            igpnt += nrows;
        }
    }
    
    // comopute remaining grid points block
    
    nrows = molecularGrid.getNumberOfGridPoints() % nrows;
    
    if (nrows > 0)
    {
        // allocate AOs grid matrices for bra and ket
        
        CDenseMatrix bmat(nrows, gtoContainer->getNumberOfAtomicOrbitals());
        
        CDenseMatrix kmat(nrows, gtoContainer->getNumberOfAtomicOrbitals());
        
        _compGtosMatrixForLDA(bmat, gtoContainer, molecularGrid, igpnt, nrows);
        
        auto gao = bmat.values();
        
        for (int32_t j = 0; j < nrows; j++)
        {
            fovl += gw[igpnt + j] * gao[j] * gao[j];
        }
        
        _compRestrictedVXCMatrixForLDA(kmat, bmat, xcGradientGrid, molecularGrid.getWeights(), igpnt, nrows);
        
        denblas::multAtB(ksmat, 1.0, bmat, kmat);
    }
    
    printf("Overlap (0,0) = %lf\n", fovl); 
    
    // compute exchange-correlation energy and number of electrons
    
    auto xcdat = _compEnergyAndDensity(xcGradientGrid, densityGrid, molecularGrid);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(std::get<0>(xcdat));
    
    aoKohnShamMatrix.setNumberOfElectrons(std::get<1>(xcdat));
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
        
        auto grho_ab = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhob);
        
        // set up pointer to perturbed density
        
        auto rhowa = rwDensityGrid->alphaDensity(0);
        
        auto rhowb = rwDensityGrid->betaDensity(0);
        
        for (int32_t i = 0; i < braAngularComponents; i++)
        {
            auto bgto = braGtoGridBuffer.data(i);
            
            for (int32_t j = 0; j < ketAngularComponents; j++)
            {
                auto kgto = ketGtoGridBuffer.data(j);
                
                double psum = 0.0;
                
                for (int32_t k = 0; k < ngpoints; k++)
                {
                    psum += gridWeights[gridOffset + k] * bgto[k] * kgto[k] * (grho_aa[gridOffset + k] * rhowa[gridOffset + k] +
                                                                               
                                                                               grho_ab[gridOffset + k] * rhowb[gridOffset + k]);
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
        
        auto grho_aa = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        auto grho_ab = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhob);
        
        auto gmix_aa = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::grada);
        
        auto gmix_ab = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::gradb);
        
        auto gmix_ac = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::gradab);
        
        auto gmix_bc = xcHessianGrid->xcHessianValues(xcvars::rhob, xcvars::gradab);
        
        auto ggrad_aa = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::grada);
        
        auto ggrad_ab = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::gradb);
        
        auto ggrad_ac = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::gradab);
        
        auto ggrad_bc = xcHessianGrid->xcHessianValues(xcvars::gradb, xcvars::gradab);
        
        auto ggrad_cc = xcHessianGrid->xcHessianValues(xcvars::gradab, xcvars::gradab);
        
        auto ggrad_a = xcGradientGrid->xcGradientValues(xcvars::grada);
        
        auto ggrad_c = xcGradientGrid->xcGradientValues(xcvars::gradab);
        
        // set up pointers to ground state density gradient norms
        
        auto ngrada = gsDensityGrid->alphaDensityGradient(0);
        
        auto ngradb = gsDensityGrid->betaDensityGradient(0);
        
        auto grada_x = gsDensityGrid->alphaDensityGradientX(0);
        
        auto grada_y = gsDensityGrid->alphaDensityGradientY(0);
        
        auto grada_z = gsDensityGrid->alphaDensityGradientZ(0);
        
        auto gradb_x = gsDensityGrid->betaDensityGradientX(0);
        
        auto gradb_y = gsDensityGrid->betaDensityGradientY(0);
        
        auto gradb_z = gsDensityGrid->betaDensityGradientZ(0);
        
        // set up pointers to perturbed density gradient norms
        
        auto rhowa = rwDensityGrid->alphaDensity(0);
        
        auto rhowb = rwDensityGrid->betaDensity(0);
        
        auto gradwa_x = rwDensityGrid->alphaDensityGradientX(0);
        
        auto gradwa_y = rwDensityGrid->alphaDensityGradientY(0);
        
        auto gradwa_z = rwDensityGrid->alphaDensityGradientZ(0);
        
        auto gradwb_x = rwDensityGrid->betaDensityGradientX(0);
        
        auto gradwb_y = rwDensityGrid->betaDensityGradientY(0);
        
        auto gradwb_z = rwDensityGrid->betaDensityGradientZ(0);
        
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
                    
                    double znva = 1.0 / ngrada[gridOffset + k];
                    
                    double znvb = 1.0 / ngradb[gridOffset + k];
                    
                    double rxa = znva * grada_x[gridOffset + k];
                    
                    double rya = znva * grada_y[gridOffset + k];
                    
                    double rza = znva * grada_z[gridOffset + k];
                    
                    double rxb = znvb * gradb_x[gridOffset + k];
                    
                    double ryb = znvb * gradb_y[gridOffset + k];
                    
                    double rzb = znvb * gradb_z[gridOffset + k];
                    
                    double rxwa = gradwa_x[gridOffset + k];
                    
                    double rywa = gradwa_y[gridOffset + k];
                    
                    double rzwa = gradwa_z[gridOffset + k];
                    
                    double rxwb = gradwb_x[gridOffset + k];
                    
                    double rywb = gradwb_y[gridOffset + k];
                    
                    double rzwb = gradwb_z[gridOffset + k];
                    
                    // GTOs values
                    
                    auto a0 = bgto[k] * kgto[k];
                    
                    auto ax = bgto_x[k] * kgto[k] + bgto[k] * kgto_x[k];
                    
                    auto ay = bgto_y[k] * kgto[k] + bgto[k] * kgto_y[k];
                    
                    auto az = bgto_z[k] * kgto[k] + bgto[k] * kgto_z[k];
                    
                    //  variations of functionals variables
                    
                    double zetaa = rxwa * rxa + rywa * rya + rzwa * rza;
                    
                    double zetab = rxwb * rxb + rywb * ryb + rzwb * rzb;
                    
                    double zetac = grada_x[gridOffset + k] * rxwb + grada_y[gridOffset + k] * rywb
                    
                                 + grada_z[gridOffset + k] * rzwb + gradb_x[gridOffset + k] * rxwa
                    
                                 + gradb_y[gridOffset + k] * rywa + gradb_z[gridOffset + k] * rzwa;
                    
                    // first contribution
                    
                    double fac0 = gmix_aa[gridOffset + k] * zetaa + gmix_ab[gridOffset + k] * zetab
                    
                                + gmix_ac[gridOffset + k] * zetac + grho_aa[gridOffset + k] * rhowa[gridOffset + k]
                    
                                + grho_ab[gridOffset + k] * rhowb[gridOffset + k];
                    
                    psum += w * a0 * fac0;
                    
                    // second contribution
                    
                    double facr = gmix_aa[gridOffset + k] * rhowa[gridOffset + k] + gmix_ab[gridOffset + k] * rhowb[gridOffset + k]
                    
                                + ggrad_aa[gridOffset + k] * zetaa + ggrad_ab[gridOffset + k] * zetab + ggrad_ac[gridOffset + k] * zetac;
                    
                    double ar = ax * rxa + ay * rya + az * rza;
                    
                    psum += w * facr * ar;
                    
                    // third contribution
                    
                    double facz = gmix_ac[gridOffset + k] * rhowa[gridOffset + k] +  gmix_bc[gridOffset + k] * rhowb[gridOffset + k]
                    
                                + ggrad_ac[gridOffset + k] * zetaa + ggrad_bc[gridOffset + k] * zetab + ggrad_cc[gridOffset + k] * zetac;
                    
                    double arb = ax * grada_x[gridOffset + k] + ay * grada_y[gridOffset + k] + az * grada_z[gridOffset + k];
                    
                    psum += w * facz * arb;
                    
                    // fourth contribution
                    
                    double ab = ax * rxwa + ay * rywa + az * rzwa - ar * zetaa;
                    
                    psum += w * znva * ggrad_a[gridOffset + k] * ab;
                    
                    // fifth contribution
                    
                    double abw = ax * rxwa + ay * rywa + az * rzwa;
                    
                    psum += w * ggrad_c[gridOffset + k] * abw;
                }
                
                ppvals[i * ketAngularComponents + j] = psum;
            }
        }
        
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

int32_t
CXCIntegrator::_getNumberOfGridRows() const
{
    // FIX ME: add basis size dependence if needed
    
    return 10000;
}

void
CXCIntegrator::_compGtosMatrixForLDA(      CDenseMatrix&   gtoMatrix,
                                     const CGtoContainer*  gtoContainer,
                                     const CMolecularGrid& molecularGrid,
                                     const int32_t         gridOffset,
                                     const int32_t         nGridPoints) const
{
    // set up OMP tasks
    
    COMPTasks omptaks(3);
    
    omptaks.set(nGridPoints);
    
    auto ntasks = omptaks.getNumberOfTasks();
    
    auto tbsizes = omptaks.getTaskSizes();
    
    auto tbpositions = omptaks.getTaskPositions();
    
    // set up pointer to molecular grid weigths
    
    auto mgx = molecularGrid.getCoordinatesX();
    
    auto mgy = molecularGrid.getCoordinatesY();
    
    auto mgz = molecularGrid.getCoordinatesZ();
    
    // set up GTOs data
    
    auto pgaos = gtoMatrix.values();
    
    auto nrows = gtoMatrix.getNumberOfRows();
    
    auto ncols = gtoMatrix.getNumberOfColumns();
    
    // generate density on grid points
    
    #pragma omp parallel shared(pgaos, nrows, ncols, tbsizes, tbpositions, ntasks, mgx, mgy, mgz)
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
                    _compGtosValuesForLDA(pgaos, nrows, ncols, gtoContainer, mgx, mgy, mgz,
                                          tbposition + gridOffset, tbsize);
                }
            }
        }
    }
}

void
CXCIntegrator::_compGtosValuesForLDA(      double*        gtoMatrix,
                                     const int32_t        nRows,
                                     const int32_t        nColumns,
                                     const CGtoContainer* gtoContainer,
                                     const double*        gridCoordinatesX,
                                     const double*        gridCoordinatesY,
                                     const double*        gridCoordinatesZ,
                                     const int32_t        gridOffset,
                                     const int32_t        nGridPoints) const
{
    // local copy of GTOs containers
    
    auto gtovec = CGtoContainer(*gtoContainer);
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        // angular momentum data for bra and ket
        
        auto bang = bgtos.getAngularMomentum();
        
        // set up Cartesian GTOs buffer
        
        auto nvcomp = xcfun_components(xcfun::lda);
        
        auto bncart = angmom::to_CartesianComponents(bang);
        
        auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();
        
        // set up spherical GTOs buffer
        
        auto bnspher = angmom::to_SphericalComponents(bang);
        
        CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);
        
        // loop over contracted GTOs
        
        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
        {
            // compute j-th GTO values on batch of grid points
            
            gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                           gridOffset, bgtos, j, xcfun::lda);
            
            // distribute j-th GTO values into grid values matrix
            
            for (int32_t k = 0; k < bnspher; k++)
            {
                auto idx = (bgtos.getIdentifiers(k))[j];
                
                auto bgaos = bspherbuff.data(k);
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    gtoMatrix[l * nColumns + idx] = bgaos[l];
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedVXCMatrixForLDA(      CDenseMatrix&    ketGtoMatrix,
                                              const CDenseMatrix&    braGtoMatrix,
                                              const CXCGradientGrid& xcGradientGrid,
                                              const double*          gridWeights, 
                                              const int32_t          gridOffset,
                                              const int32_t          nGridPoints) const
{
    // set up pointers to GTOs matrices
    
    auto kvxc = ketGtoMatrix.values();
    
    auto bgao = braGtoMatrix.values();
    
    // set up dimensions of GTOs matrices
    
    auto ncols = ketGtoMatrix.getNumberOfColumns();
    
    // set up pointer to gradient of XC functional
    
    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
    
    #pragma omp parallel for
    for (int32_t i = 0; i < nGridPoints; i++)
    {
        auto fact = grhoa[gridOffset + i] * gridWeights[gridOffset + i];
        
        auto ioff = i * ncols;
        
        for (int32_t j = 0; j < ncols; j++)
        {
            kvxc[ioff + j] = fact * bgao[ioff + j];
        }
    }
}

std::tuple<double, double>
CXCIntegrator::_compEnergyAndDensity(const CXCGradientGrid& xcGradientGrid,
                                     const CDensityGrid&    densityGrid,
                                     const CMolecularGrid&  molecularGrid) const
{
    // set up pointers to density grid
    
    auto rhoa = densityGrid.alphaDensity(0);
    
    auto rhob = densityGrid.betaDensity(0);
    
    // set up pointer to exchange-correlation energy grid
    
    auto efunc = xcGradientGrid.xcFunctionalValues();
    
    // set up pointer to grid weights
    
    auto gw = molecularGrid.getWeights();
    
    // set up number of grid points
    
    auto gpoints = molecularGrid.getNumberOfGridPoints();
    
    // initialize exchange-correlation energy and gradient
    
    double xcene = 0.0, nele = 0.0;

    #pragma omp parallel for reduction(+:nele) reduction(+:xcene)
    for (int32_t i = 0; i < gpoints; i++)
    {
        nele += gw[i] * (rhoa[i] + rhob[i]);
        
        xcene += gw[i] * efunc[i];
    }
    
    return std::make_tuple(xcene, nele);
}


