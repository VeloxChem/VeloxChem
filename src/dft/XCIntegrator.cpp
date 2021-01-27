//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "XCIntegrator.hpp"

#include <cmath>

#include <mpi.h>

#include "MpiFunc.hpp"
#include "DensityGridDriver.hpp"
#include "FunctionalParser.hpp"
#include "AOKohnShamMatrix.hpp"
#include "OMPTasks.hpp"
#include "AngularMomentum.hpp"
#include "GtoFunc.hpp"
#include "DenseLinearAlgebra.hpp"
#include "MathConst.hpp"


CXCIntegrator::CXCIntegrator(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _locComm = comm;
    
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
            
            ksmat.zero(); 
            
            if (fvxc.getFunctionalType() == xcfun::lda)
            {
                _compRestrictedContributionForLDA(ksmat, gtovec, vxcgrid, dgrid, mgrid);
            }
            
            if (fvxc.getFunctionalType() == xcfun::gga)
            {
                _compRestrictedContributionForGGA(ksmat, gtovec, vxcgrid, dgrid, mgrid);
            }
            
            // symmetrize Kohn-Sham matrix
            
            ksmat.symmetrize();
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
    
    // compute gradient and hessian of exchange-correlation functional
    
    CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());
    
    CXCHessianGrid vxc2grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());
    
    fvxc.compute(vxcgrid, gsdengrid);
    
    fvxc.compute(vxc2grid, gsdengrid);
    
    // compute perturbed density grid
    
    auto rwdengrid = dgdrv.generate(rwDensityMatrix, molecule, basis, mgrid, fvxc.getFunctionalType());

    // set up number of perturbedf density matrix dimensions
    
    auto ndmat = rwDensityMatrix.getNumberOfDensityMatrices();
    
    auto nrow = gsDensityMatrix.getNumberOfRows(0);
    
    auto ncol = gsDensityMatrix.getNumberOfColumns(0);
    
    if (rwDensityMatrix.isRestricted())
    {
        // allocate perturbed Kohn-Sham matrix
        
        CAOKohnShamMatrix ksmat(nrow, ncol, ndmat, true);
        
        if (fvxc.getFunctionalType() == xcfun::lda)
        {
            _compRestrictedContributionForLDA(ksmat, gtovec, vxc2grid, rwdengrid, mgrid);
        }
        
        if (fvxc.getFunctionalType() == xcfun::gga)
        {
            _compRestrictedContributionForGGA(ksmat, gtovec, vxcgrid, vxc2grid, gsdengrid, rwdengrid, mgrid);
        }
        
        // symmetrize Kohn-Sham matrix
        
        ksmat.symmetrize();
        
        // add Kohn-Sham contribution to Fock matrix
        
        for (int32_t i = 0; i < ksmat.getNumberOfMatrices(); i++)
        {
            aoFockMatrix.addOneElectronMatrix(ksmat.getReferenceToMatrix(i), i); 
        }
    }
    else
    {
        // perform integration for spin unrestricted densities
        
        // FIX ME: implement partitioning and screening for density grid
    }
}

void
CXCIntegrator::_compRestrictedContributionForLDAWithNL(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                       const CGtoContainer*     gtoContainer,
                                                       const CXCGradientGrid&   xcGradientGrid,
                                                       const CDensityGrid&      densityGrid,
                                                       const CMolecularGrid&    molecularGrid) const
{
    // set up pointer to molecular grid weigths
    
    auto mgx = molecularGrid.getCoordinatesX();
    
    auto mgy = molecularGrid.getCoordinatesY();
    
    auto mgz = molecularGrid.getCoordinatesZ();
    
    auto mgw = molecularGrid.getWeights();
    
    // set up number of grid points
    
    auto gpoints = molecularGrid.getNumberOfGridPoints();
    
    // set up pointer to density matrix
    
    auto xcgridptr = &xcGradientGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    #pragma omp parallel shared(mgx, mgy, mgz, mgw, xcgridptr, ksmatprt, gpoints)
    {
        #pragma omp single nowait
        {
            auto nbra = gtoContainer->getNumberOfGtoBlocks();
            
            for (int32_t i = 0; i < nbra; i++)
            {
                auto bgtos = gtoContainer->getGtoBlock(i);
                
                for (int32_t j = i; j < nbra; j++)
                {
                    auto kgtos = gtoContainer->getGtoBlock(j);
                    
                    for (int32_t k = 0; k < bgtos.getNumberOfContrGtos(); k++)
                    {
                        #pragma omp task firstprivate(k)
                        {
                            _compRestrictedBatchForLDAWithNL(ksmatprt, bgtos, k, kgtos, xcgridptr,
                                                             mgx, mgy, mgz, mgw, gpoints);
                        }
                    }
                }
            }
        }
    }
    
    // compute exchange-correlation energy and number of electrons
    
    auto xcdat = _compEnergyAndDensity(xcGradientGrid, densityGrid, molecularGrid);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(std::get<0>(xcdat));
    
    aoKohnShamMatrix.setNumberOfElectrons(std::get<1>(xcdat));
}

void
CXCIntegrator::_compRestrictedBatchForLDAWithNL(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                                const CGtoBlock&         braGtoBlock,
                                                const int32_t            iBraContrGto,
                                                const CGtoBlock&         ketGtoBlock,
                                                const CXCGradientGrid*   xcGradientGrid,
                                                const double*            gridCoordinatesX,
                                                const double*            gridCoordinatesY,
                                                const double*            gridCoordinatesZ,
                                                const double*            gridWeights,
                                                const int32_t            nGridPoints) const
{
    // determine symmetry of bra and ket sides
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
    // angular momentum data for bra and ket sides
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    auto bncart = angmom::to_CartesianComponents(bang);
    
    auto kncart = angmom::to_CartesianComponents(kang);
    
    auto bnspher = angmom::to_SphericalComponents(bang);
    
    auto knspher = angmom::to_SphericalComponents(kang);
    
    // set up GTOs buffers for bra side
    
    auto nvcomp = xcfun_components(xcfun::lda);
    
    auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();
    
    CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);
    
    // compute GTO values on grid for bra side
    
    gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, 0,
                                   braGtoBlock, iBraContrGto, xcfun::lda);
    
    // determine GTO values screening pattern
    
    auto gpids = _getScreeningPattern(bspherbuff);
    
    // set up significant grid points
    
    auto redgrid = _getReducedGrid(gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridWeights, gpids);
    
    // compress GTO values on bra size
    
    auto bredbuff = _getReducedGtoValues(bspherbuff, gpids);
    
    // compress gradient
    
    auto redgrad = _getReducedRestrictedGradient(xcGradientGrid, gpids);
    
    // set up GTOs buffers for ket side
    
    auto ngdim = gpids.size();
    
    auto kcartbuff = (kang > 0) ? CMemBlock2D<double>(ngdim, nvcomp * kncart) : CMemBlock2D<double>();
    
    CMemBlock2D<double> kspherbuff(ngdim, nvcomp * knspher);
    
    // set up pointer to gradient data
    
    auto grhoa = redgrad.data();
    
    // set up pointer to grid weights
    
    auto gw = redgrid.data(3);
    
    // set up Kohn-Sham matrix data
    
    auto ncols = aoKohnShamMatrix->getNumberOfColumns();
    
    auto ksmat = aoKohnShamMatrix->getMatrix(0);
    
    // loop over ket side
    
    int32_t istart = (symbk) ? iBraContrGto : 0;
    
    for (int32_t i = istart; i < ketGtoBlock.getNumberOfContrGtos(); i++)
    {
        if (_isSignificantShellPair(braGtoBlock, iBraContrGto, ketGtoBlock, i))
        {
            // compute GTO values on grid for ket side
        
            gtorec::computeGtoValuesOnGrid(kspherbuff, kcartbuff, redgrid.data(0), redgrid.data(1), redgrid.data(2), 0,
                                           ketGtoBlock, i, xcfun::lda);
        
            for (int32_t j = 0; j < bnspher; j++)
            {
                auto bgaos = bredbuff.data(j);
            
                auto bidx = (braGtoBlock.getIdentifiers(j))[iBraContrGto];
            
                for (int32_t k = 0; k < knspher; k++)
                {
                    auto kgaos = kspherbuff.data(k);
                
                    auto kidx = (ketGtoBlock.getIdentifiers(k))[i];
                
                    double fval = 0.0;
                
                    #pragma omp simd reduction(+:fval) aligned(bgaos, kgaos, gw, grhoa: VLX_ALIGN)
                    for (int32_t l = 0; l < ngdim; l++)
                    {
                        fval += bgaos[l] * kgaos[l] * gw[l] * grhoa[l];
                    }
                
                    ksmat[bidx * ncols + kidx] = fval;
                    
                    ksmat[kidx * ncols + bidx] = fval;
                }
            }
        }
    }
}

bool
CXCIntegrator::_isSignificantShellPair(const CGtoBlock&      braGtoBlock,
                                       const int32_t         iBraContrGto,
                                       const CGtoBlock&      ketGtoBlock,
                                       const int32_t         iKetContrGto) const
{
    // get pi constant
    
    auto fpi = mathconst::getPiValue();
    
    // set up GTO data for bra
    
    auto bstart = (braGtoBlock.getStartPositions())[iBraContrGto];
    
    auto bend = (braGtoBlock.getEndPositions())[iBraContrGto];
    
    auto bexps = braGtoBlock.getExponents();
    
    auto bnorms = braGtoBlock.getNormFactors();
    
    auto brx = braGtoBlock.getCoordinatesX();
    
    auto bry = braGtoBlock.getCoordinatesY();
    
    auto brz = braGtoBlock.getCoordinatesZ();
    
    // compute zero order overlap for bra
    
    double bovl = 0.0;
    
    for (int32_t i = bstart; i < bend; i++)
    {
        bovl += bnorms[i] * bnorms[i] * std::pow(0.5 * fpi / bexps[i], 1.5);
        
        for (int32_t j = i + 1; j < bend; j++)
        {
            bovl += 2.0 * bnorms[i] * bnorms[j] * std::pow(fpi / (bexps[i] + bexps[j]), 1.5);
        }
    }
    
    // set up GTO data for ket
    
    auto kstart = (ketGtoBlock.getStartPositions())[iKetContrGto];
    
    auto kend = (ketGtoBlock.getEndPositions())[iKetContrGto];
    
    auto kexps = ketGtoBlock.getExponents();
    
    auto knorms = ketGtoBlock.getNormFactors();
    
    auto krx = ketGtoBlock.getCoordinatesX();
    
    auto kry = ketGtoBlock.getCoordinatesY();
    
    auto krz = ketGtoBlock.getCoordinatesZ();
    
    // compute zero order overlap for ket
    
    double kovl = 0.0;
    
    for (int32_t i = kstart; i < kend; i++)
    {
        kovl += knorms[i] * knorms[i] * std::pow(0.5 * fpi / kexps[i], 1.5);
        
        for (int32_t j = i + 1; j < kend; j++)
        {
            kovl += 2.0 * knorms[i] * knorms[j] * std::pow(fpi / (kexps[i] + kexps[j]), 1.5);
        }
    }
    
    // compute zero order bra/ket overlap
    
    double bkovl = 0.0;
    
    for (int32_t i = bstart; i < bend; i++)
    {
        for (int32_t j = kstart; j < kend; j++)
        {
            auto fab = 1.0 / (bexps[i] + kexps[j]);
            
            auto rabx = brx[i] - krx[j];
            
            auto raby = bry[i] - kry[j];
            
            auto rabz = brz[i] - krz[j];
            
            auto r2ab = rabx * rabx + raby * raby + rabz * rabz;
            
            auto fact = bnorms[i] * knorms[j] * std::pow(fpi * fab, 1.5);
            
            bkovl += fact * std::exp(-bexps[i] * kexps[j] * fab * r2ab);
        }
    }
    
    if (std::fabs(bkovl / std::sqrt(bovl * kovl)) > _thresholdOfDensity) return true;
    
    return false;
}

void
CXCIntegrator::_compRestrictedContributionForLDA(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                 const CGtoContainer*     gtoContainer,
                                                 const CXCGradientGrid&   xcGradientGrid,
                                                 const CDensityGrid&      densityGrid,
                                                 const CMolecularGrid&    molecularGrid) const
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
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // compute contributions to Kohn-Sham matrix from blocks of grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, xcgridptr, ksmatprt)
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
                    _compRestrictedBatchForLDA(ksmatprt, gtoContainer, xcgridptr, mgx, mgy, mgz, mgw,
                                               tbposition, tbsize);
                }
            }
        }
    }
    
    // compute exchange-correlation energy and number of electrons
    
    auto xcdat = _compEnergyAndDensity(xcGradientGrid, densityGrid, molecularGrid);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(std::get<0>(xcdat));
    
    aoKohnShamMatrix.setNumberOfElectrons(std::get<1>(xcdat));
}

void
CXCIntegrator::_compRestrictedContributionForLDA(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                 const CGtoContainer*     gtoContainer,
                                                 const CXCHessianGrid&    xcHessianGrid,
                                                 const CDensityGrid&      rwDensityGrid,
                                                 const CMolecularGrid&    molecularGrid) const
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
    
    // set up pointer to exchange-correlation derrivatives
    
    auto xchessptr = &xcHessianGrid;
    
    // set up pointer to perturbed density grid
    
    auto rwdenptr = &rwDensityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // compute contributions to Kohn-Sham matrix from blocks of grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, xchessptr, rwdenptr, ksmatprt)
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
                    _compRestrictedBatchForLDA(ksmatprt, gtoContainer, xchessptr, rwdenptr, mgx, mgy, mgz, mgw,
                                               tbposition, tbsize);
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedContributionForGGA(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                 const CGtoContainer*     gtoContainer,
                                                 const CXCGradientGrid&   xcGradientGrid,
                                                 const CDensityGrid&      densityGrid,
                                                 const CMolecularGrid&    molecularGrid) const
{
    // initialize Kohn-Sham matrix to zero
    
    aoKohnShamMatrix.zero();
    
    // set up OMP tasks
    
    COMPTasks omptaks(2);
    
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
    
    // set up number of GTOs blocks
    
    auto nblocks = gtoContainer->getNumberOfGtoBlocks();
    
    // loop over pairs of GTOs blocks
    
    for (int32_t i = 0; i < nblocks; i++)
    {
        auto bgtos = gtoContainer->getGtoBlock(i);
        
        for (int32_t j = i; j < nblocks; j++)
        {
            auto kgtos = gtoContainer->getGtoBlock(j);
            
            // compute contributions to Kohn-Sham matrix from blocks of grid points
            
            #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, dgridptr, xcgridptr, ksmatprt)
            {
                #pragma omp single nowait
                {
                    for (int32_t k = 0; k < ntasks; k++)
                    {
                        // set up task parameters
                        
                        auto tbsize = tbsizes[k];
                        
                        auto tbposition = tbpositions[k];
                        
                        // generate task
                        
                        #pragma omp task firstprivate(tbsize, tbposition)
                        {
                            _compRestrictedBatchForGGA(ksmatprt, bgtos, kgtos, xcgridptr, dgridptr,
                                                       mgx, mgy, mgz, mgw, tbposition, tbsize);
                        }
                    }
                }
            }
                
            
        }
    }
    
    // compute exchange-correlation energy and number of electrons
    
    auto xcdat = _compEnergyAndDensity(xcGradientGrid, densityGrid, molecularGrid);
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(std::get<0>(xcdat));
    
    aoKohnShamMatrix.setNumberOfElectrons(std::get<1>(xcdat));
}

void
CXCIntegrator::_compRestrictedContributionForGGA(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                 const CGtoContainer*     gtoContainer,
                                                 const CXCGradientGrid&   xcGradientGrid,
                                                 const CXCHessianGrid&    xcHessianGrid,
                                                 const CDensityGrid&      gsDensityGrid,
                                                 const CDensityGrid&      rwDensityGrid,
                                                 const CMolecularGrid&    molecularGrid) const
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
    
    // set up pointer to exchange-correlation derrivatives
    
    auto xcgradptr = &xcGradientGrid;
    
    auto xchessptr = &xcHessianGrid;
    
    // set up pointer to density grid
    
    auto gsdenptr = &gsDensityGrid;
    
    auto rwdenptr = &rwDensityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // compute contributions to Kohn-Sham matrix from blocks of grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, xcgradptr, xchessptr,\
                                gsdenptr, rwdenptr, ksmatprt)
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
                    _compRestrictedBatchForGGA(ksmatprt, gtoContainer, xcgradptr, xchessptr, gsdenptr, rwdenptr,
                                               mgx, mgy, mgz, mgw, tbposition, tbsize);
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedBatchForLDA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoContainer*     gtoContainer,
                                          const CXCGradientGrid*   xcGradientGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            gridOffset,
                                          const int32_t            nGridPoints) const
{
    // set up number of AOs
    
    auto naos = gtoContainer->getNumberOfAtomicOrbitals();
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // set up pointer to gradient data
    
    auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);
    
    // allocate XC contribution buffer
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    CMemBlock<double> vxcbuf(bufdim * naos);
    
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
            
            _distRestrictedBatchForLDA(aoKohnShamMatrix, vxcbuf, grhoa, gaos, gridWeights, gridOffset, igpnt, blockdim);
            
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
        
        _distRestrictedBatchForLDA(aoKohnShamMatrix, vxcbuf, grhoa, gaos, gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compRestrictedBatchForLDA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoContainer*     gtoContainer,
                                          const CXCHessianGrid*    xcHessianGrid,
                                          const CDensityGrid*      rwDensityGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            gridOffset,
                                          const int32_t            nGridPoints) const
{
    // set up number of AOs
    
    auto naos = gtoContainer->getNumberOfAtomicOrbitals();
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // allocate XC contribution buffer
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    CMemBlock<double> vxcbuf(bufdim * naos);
    
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
            
            _distRestrictedBatchForLDA(aoKohnShamMatrix, vxcbuf, xcHessianGrid, rwDensityGrid, gaos,
                                       gridWeights, gridOffset, igpnt, blockdim);
            
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
        
        _distRestrictedBatchForLDA(aoKohnShamMatrix, vxcbuf, xcHessianGrid, rwDensityGrid, gaos, gridWeights,
                                   gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoContainer*     gtoContainer,
                                          const CXCGradientGrid*   xcGradientGrid,
                                          const CDensityGrid*      densityGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            gridOffset,
                                          const int32_t            nGridPoints) const
{
    // set up local copy of density grid
    
    auto dengrid = *densityGrid;
    
    // set up local copy of exchange-correlation gradient
    
    auto xcgrad = *xcGradientGrid;
    
    // set up number of AOs
    
    auto naos = gtoContainer->getNumberOfAtomicOrbitals();
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // allocate XC contribution buffer
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    CMemBlock<double> vxcbuf(bufdim * naos);
    
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
            gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedBatchForGGA(aoKohnShamMatrix, vxcbuf, xcgrad, dengrid, gaos, gaox, gaoy, gaoz,
                                       gridWeights, gridOffset, igpnt, blockdim);
            
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
        
        gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                        gridOffset, igpnt, blockdim);
        
        _distRestrictedBatchForGGA(aoKohnShamMatrix, vxcbuf, xcgrad, dengrid, gaos, gaox, gaoy, gaoz,
                                   gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoBlock&         braGtoBlock,
                                          const CGtoBlock&         ketGtoBlock,
                                          const CXCGradientGrid*   xcGradientGrid,
                                          const CDensityGrid*      densityGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            gridOffset,
                                          const int32_t            nGridPoints) const
{
    // determine symmetry
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
    // set up angular momentum data
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // set up number of angular components
    
    auto bcomps = angmom::to_SphericalComponents(bang);
    
    auto kcomps = angmom::to_SphericalComponents(kang);
    
    // set up dimensions of bra anf ket sides
    
    auto bdim = bcomps * braGtoBlock.getNumberOfContrGtos();
    
    auto kdim = kcomps * ketGtoBlock.getNumberOfContrGtos();
    
    // initializing partial Kohn-Sham matrix
    
    CDenseMatrix submat(bdim, kdim);
    
    submat.zero(); 
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // set up current grid point
    
    int32_t igpnt = 0;
    
    // loop over grid points blocks
    
    if (nblocks > 0)
    {
        // bra GTOs values grid
        
        CMemBlock2D<double> baos(blockdim, bdim);
        
        CMemBlock2D<double> baox(blockdim, bdim);
        
        CMemBlock2D<double> baoy(blockdim, bdim);
        
        CMemBlock2D<double> baoz(blockdim, bdim);
        
        // ket GTOs values grid if needed
        
        auto kaos = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaox = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaoy = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaoz = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        for (int32_t i = 0; i < nblocks; i++)
        {
            gtorec::computeGtosValuesForGGA(baos, baox, baoy, baoz, braGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            if (symbk)
            {
                _distRestrictedBatchForGGA(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
                
                _distRestrictedBatchForGGA(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                           blockdim);
            }
            
            igpnt += blockdim;
        }
    }
    
    // comopute remaining grid points block
    
    blockdim = nGridPoints % blockdim;
    
    if (blockdim > 0)
    {
        // bra GTOs values grid
        
        CMemBlock2D<double> baos(blockdim, bdim);
        
        CMemBlock2D<double> baox(blockdim, bdim);
        
        CMemBlock2D<double> baoy(blockdim, bdim);
        
        CMemBlock2D<double> baoz(blockdim, bdim);
        
        // ket GTOs values grid if needed
        
        auto kaos = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaox = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaoy = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        auto kaoz = (symbk) ? CMemBlock2D<double>() : CMemBlock2D<double>(blockdim, kdim);
        
        gtorec::computeGtosValuesForGGA(baos, baox, baoy, baoz, braGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                        gridCoordinatesZ, gridOffset, igpnt, blockdim);
        
        if (symbk)
        {
            _distRestrictedBatchForGGA(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       gridWeights, gridOffset, igpnt, blockdim);
        }
        else
        {
            gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedBatchForGGA(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                       blockdim);
        }
    }
    
    // add submatrix to Kohn-Sham matrix
    
    #pragma omp critical
    _addSubMatrix(aoKohnShamMatrix, submat, braGtoBlock, ketGtoBlock);
}

void
CXCIntegrator::_compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoContainer*     gtoContainer,
                                          const CXCGradientGrid*   xcGradientGrid,
                                          const CXCHessianGrid*    xcHessianGrid,
                                          const CDensityGrid*      gsDensityGrid,
                                          const CDensityGrid*      rwDensityGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            gridOffset,
                                          const int32_t            nGridPoints) const
{
    // set up number of AOs
    
    auto naos = gtoContainer->getNumberOfAtomicOrbitals();
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // allocate XC contribution buffer
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    CMemBlock<double> vxcbuf(bufdim * naos);
    
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
            gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedBatchForGGA(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, gsDensityGrid, rwDensityGrid,
                                       gaos, gaox, gaoy, gaoz, gridWeights, gridOffset, igpnt, blockdim);
            
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
        
        gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                        gridOffset, igpnt, blockdim);
        
        _distRestrictedBatchForGGA(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, gsDensityGrid, rwDensityGrid,
                                   gaos, gaox, gaoy, gaoz, gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_distRestrictedBatchForLDA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const double*              xcGradient,
                                          const CMemBlock2D<double>& gtoValues,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up pointer to Kohn-Sham matrix
    
    auto ksmat = aoKohnShamMatrix->getKohnSham();
    
    // set up AOs blocks
    
    auto naos = aoKohnShamMatrix->getNumberOfRows();
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    auto nblock = naos / bufdim;
    
    // loop over AOs blocks
    
    int32_t curao = 0;
    
    for (int32_t i = 0; i < nblock; i++)
    {
        // compute block of Kohn-Sham matrix elements
        
        int32_t idx = 0;
        
        for (int32_t j = curao; j < (curao + bufdim); j++)
        {
            auto bgaos = gtoValues.data(j);
            
            for (int32_t k = j; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);
                
                double fvxc = 0.0;
                
                auto loff = gridOffset + gridBlockPosition;
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    fvxc += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient[loff + l];
                }
                
                xcBuffer.at(idx) = fvxc;
                
                idx++;
            }
        }
        
        // distribute block of Kohn-Sham matrix elements
        
        #pragma omp critical
        {
            idx = 0;
            
            for (int32_t j = curao; j < (curao + bufdim); j++)
            {
                for (int32_t k = j; k < naos; k++)
                {
                    auto fvxc = (j == k) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                
                    ksmat[j * naos + k] += fvxc;
                    
                    idx++;
                }
            }
        }
        
        // update AOs counter
        
        curao += bufdim;
    }
    
    // loop over last block
    
    if (naos % bufdim != 0)
    {
        // compute last block of Kohn-Sham matrix elements
        
        int32_t idx = 0;
        
        for (int32_t i = curao; i < naos; i++)
        {
            auto bgaos = gtoValues.data(i);
            
            for (int32_t j = i; j < naos; j++)
            {
                auto kgaos = gtoValues.data(j);
                
                double fvxc = 0.0;
                
                auto koff = gridOffset + gridBlockPosition;
                
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    fvxc += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient[koff + k];
                }
                
                xcBuffer.at(idx) = fvxc;
                
                idx++;
            }
        }
        
        // distribute last block of Kohn-Sham matrix elements
        
        #pragma omp critical
        {
            idx = 0;
            
            for (int32_t i = curao; i < naos; i++)
            {
                for (int32_t j = i; j < naos; j++)
                {
                    auto fvxc = (i == j) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                    
                    ksmat[i * naos + j] += fvxc;
                    
                    idx++;
                }
            }
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGGA(      CDenseMatrix&        subMatrix,
                                          const CXCGradientGrid*     xcGradientGrid,
                                          const CDensityGrid*        densityGrid,
                                          const CMemBlock2D<double>& gtoValues,
                                          const CMemBlock2D<double>& gtoValuesX,
                                          const CMemBlock2D<double>& gtoValuesY,
                                          const CMemBlock2D<double>& gtoValuesZ,
                                          const CGtoBlock&           gtoBlock,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up pointers to gradient data
    
    auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);
    
    auto ggrada = xcGradientGrid->xcGradientValues(xcvars::grada);
    
    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms
    
    auto ngrada = densityGrid->alphaDensityGradient(0);
    
    auto gradax = densityGrid->alphaDensityGradientX(0);
    
    auto graday = densityGrid->alphaDensityGradientY(0);
    
    auto gradaz = densityGrid->alphaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrix.getNumberOfRows();
    
    // set up pointer to submatrix values
    
    auto bmat = subMatrix.values();
    
    for (int32_t i = 0; i < bdim; i++)
    {
        auto bgaos = gtoValues.data(i);
        
        auto bgaox = gtoValuesX.data(i);
        
        auto bgaoy = gtoValuesY.data(i);
        
        auto bgaoz = gtoValuesZ.data(i);
        
        for (int32_t j = i; j < bdim; j++)
        {
            auto kgaos = gtoValues.data(j);
            
            auto kgaox = gtoValuesX.data(j);
            
            auto kgaoy = gtoValuesY.data(j);
            
            auto kgaoz = gtoValuesZ.data(j);
         
            // sum values
                
            double fvxc = 0.0;
                
            auto koff = gridOffset + gridBlockPosition;
                
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];
                    
                double gx = gradax[koff + k];
                    
                double gy = graday[koff + k];
                    
                double gz = gradaz[koff + k];
                    
                fvxc += w * bgaos[k] * kgaos[k] * grhoa[koff + k];
                    
                double fgrd = w * (ggrada[koff + k] / ngrada[koff + k] + ggradab[koff + k]);
                    
                fvxc += fgrd * bgaos[k] * (gx * kgaox[k] + gy * kgaoy[k] + gz * kgaoz[k]);
                    
                fvxc += fgrd * kgaos[k] * (gx * bgaox[k] + gy * bgaoy[k] + gz * bgaoz[k]);
            }
            
            if (i == j) fvxc *= 0.5;
            
            bmat[i * bdim + j] += fvxc;
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGGA(      CDenseMatrix&        subMatrix,
                                          const CXCGradientGrid*     xcGradientGrid,
                                          const CDensityGrid*        densityGrid,
                                          const CMemBlock2D<double>& braGtoValues,
                                          const CMemBlock2D<double>& braGtoValuesX,
                                          const CMemBlock2D<double>& braGtoValuesY,
                                          const CMemBlock2D<double>& braGtoValuesZ,
                                          const CGtoBlock&           braGtoBlock,
                                          const CMemBlock2D<double>& ketGtoValues,
                                          const CMemBlock2D<double>& ketGtoValuesX,
                                          const CMemBlock2D<double>& ketGtoValuesY,
                                          const CMemBlock2D<double>& ketGtoValuesZ,
                                          const CGtoBlock&           ketGtoBlock,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up pointers to gradient data
    
    auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);
    
    auto ggrada = xcGradientGrid->xcGradientValues(xcvars::grada);
    
    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms
    
    auto ngrada = densityGrid->alphaDensityGradient(0);
    
    auto gradax = densityGrid->alphaDensityGradientX(0);
    
    auto graday = densityGrid->alphaDensityGradientY(0);
    
    auto gradaz = densityGrid->alphaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrix.getNumberOfRows();
    
    auto kdim = subMatrix.getNumberOfColumns();
    
    // set up pointer to submatrix values
    
    auto bkmat = subMatrix.values();
    
    for (int32_t i = 0; i < bdim; i++)
    {
        auto bgaos = braGtoValues.data(i);
        
        auto bgaox = braGtoValuesX.data(i);
        
        auto bgaoy = braGtoValuesY.data(i);
        
        auto bgaoz = braGtoValuesZ.data(i);
        
        for (int32_t j = 0; j < kdim; j++)
        {
            auto kgaos = ketGtoValues.data(j);
            
            auto kgaox = ketGtoValuesX.data(j);
            
            auto kgaoy = ketGtoValuesY.data(j);
            
            auto kgaoz = ketGtoValuesZ.data(j);
            
            // sum values
            
            double fvxc = 0.0;
            
            auto koff = gridOffset + gridBlockPosition;
            
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];
                
                double gx = gradax[koff + k];
                
                double gy = graday[koff + k];
                
                double gz = gradaz[koff + k];
                
                fvxc += w * bgaos[k] * kgaos[k] * grhoa[koff + k];
                
                double fgrd = w * (ggrada[koff + k] / ngrada[koff + k] + ggradab[koff + k]);
                
                fvxc += fgrd * bgaos[k] * (gx * kgaox[k] + gy * kgaoy[k] + gz * kgaoz[k]);
                
                fvxc += fgrd * kgaos[k] * (gx * bgaox[k] + gy * bgaoy[k] + gz * bgaoz[k]);
            }
            
            bkmat[i * kdim + j] += fvxc;
        }
    }
    
}

void
CXCIntegrator::_distRestrictedBatchForLDA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CDensityGrid*        rwDensityGrid,
                                          const CMemBlock2D<double>& gtoValues,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up AOs blocks
    
    auto naos = aoKohnShamMatrix->getNumberOfRows();
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    auto nblock = naos / bufdim;
    
    // set up pointers to exchange-correlation functional derrivatives
    
    auto grho_aa = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhoa);
    
    auto grho_ab = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhob);
    
    // set up loop over number of matrices
    
    auto nmat = aoKohnShamMatrix->getNumberOfMatrices();
    
    for (int32_t i = 0; i < nmat; i++)
    {
        // set up pointer to Kohn-Sham matrix
        
        auto ksmat = aoKohnShamMatrix->getMatrix(i);
        
        // set up pointer to perturbed density
        
        auto rhowa = rwDensityGrid->alphaDensity(i);
        
        auto rhowb = rwDensityGrid->betaDensity(i);
        
        // loop over AOs blocks
        
        int32_t curao = 0;
        
        for (int32_t j = 0; j < nblock; j++)
        {
            // compute block of Kohn-Sham matrix elements
            
            int32_t idx = 0;
            
            for (int32_t k = curao; k < (curao + bufdim); k++)
            {
                auto bgaos = gtoValues.data(k);
                
                for (int32_t l = k; l < naos; l++)
                {
                    auto kgaos = gtoValues.data(l);
                    
                    double fvxc = 0.0;
                    
                    auto moff = gridOffset + gridBlockPosition;
                    
                    for (int32_t m = 0; m < nGridPoints; m++)
                    {
                        fvxc += bgaos[m] * kgaos[m] * gridWeights[moff + m]
                        
                             * (grho_aa[moff + m] * rhowa[moff + m] + grho_ab[moff + m] * rhowb[moff + m]);
                    }
                    
                    xcBuffer.at(idx) = fvxc;
                    
                    idx++;
                }
            }
            
            // distribute block of Kohn-Sham matrix elements
        
            #pragma omp critical
            {
                idx = 0;
            
                for (int32_t k = curao; k < (curao + bufdim); k++)
                {
                    for (int32_t l = k; l < naos; l++)
                    {
                        auto fvxc = (k == l) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                
                        ksmat[k * naos + l] += fvxc;
                    
                        idx++;
                    }
                }
            }
        
            // update AOs counter
        
            curao += bufdim;
        }
        
        // loop over last block
        
        if (naos % bufdim != 0)
        {
            // compute last block of Kohn-Sham matrix elements
            
            int32_t idx = 0;
            
            for (int32_t j = curao; j < naos; j++)
            {
                auto bgaos = gtoValues.data(j);
                
                for (int32_t k = j; k < naos; k++)
                {
                    auto kgaos = gtoValues.data(k);
                    
                    double fvxc = 0.0;
                    
                    auto loff = gridOffset + gridBlockPosition;
                    
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        fvxc += bgaos[l] * kgaos[l] * gridWeights[loff + l]
                        
                              * (grho_aa[loff + l] * rhowa[loff + l] + grho_ab[loff + l] * rhowb[loff + l]);
                    }
                    
                    xcBuffer.at(idx) = fvxc;
                    
                    idx++;
                }
            }
            
            // distribute last block of Kohn-Sham matrix elements
            
            #pragma omp critical
            {
                idx = 0;
                
                for (int32_t j = curao; j < naos; j++)
                {
                    for (int32_t k = j; k < naos; k++)
                    {
                        auto fvxc = (j == k) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                        
                        ksmat[j * naos + k] += fvxc;
                        
                        idx++;
                    }
                }
            }
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGGA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const CXCGradientGrid&     xcGradientGrid,
                                          const CDensityGrid&        densityGrid,
                                          const CMemBlock2D<double>& gtoValues,
                                          const CMemBlock2D<double>& gtoValuesX,
                                          const CMemBlock2D<double>& gtoValuesY,
                                          const CMemBlock2D<double>& gtoValuesZ,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up pointer to Kohn-Sham matrix
    
    auto ksmat = aoKohnShamMatrix->getKohnSham();
    
    // set up pointers to gradient data
    
    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
    
    auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);
    
    auto ggradab = xcGradientGrid.xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms
    
    auto ngrada = densityGrid.alphaDensityGradient(0);
    
    auto gradax = densityGrid.alphaDensityGradientX(0);
    
    auto graday = densityGrid.alphaDensityGradientY(0);
    
    auto gradaz = densityGrid.alphaDensityGradientZ(0);
    
    // set up AOs blocks
    
    auto naos = aoKohnShamMatrix->getNumberOfRows();
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    auto nblock = naos / bufdim;
    
    // loop over AOs blocks
    
    int32_t curao = 0;
    
    for (int32_t i = 0; i < nblock; i++)
    {
        // compute block of Kohn-Sham matrix elements
        
        int32_t idx = 0;
        
        for (int32_t j = curao; j < (curao + bufdim); j++)
        {
            auto bgaos = gtoValues.data(j);
            
            auto bgaox = gtoValuesX.data(j);
            
            auto bgaoy = gtoValuesY.data(j);
            
            auto bgaoz = gtoValuesZ.data(j);
            
            for (int32_t k = j ; k < naos; k++)
            {
                auto kgaos = gtoValues.data(k);
                
                auto kgaox = gtoValuesX.data(k);
                
                auto kgaoy = gtoValuesY.data(k);
                
                auto kgaoz = gtoValuesZ.data(k);
                
                double fvxc = 0.0;
                
                auto loff = gridOffset + gridBlockPosition;
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    double w = gridWeights[loff + l];
                    
                    double gx = gradax[loff + l];
                
                    double gy = graday[loff + l];
                    
                    double gz = gradaz[loff + l];
                    
                    fvxc += w * bgaos[l] * kgaos[l] * grhoa[loff + l];
                    
                    double fgrd = w * (ggrada[loff + l] / ngrada[loff + l] + ggradab[loff + l]);
                    
                    fvxc += fgrd * bgaos[l] * (gx * kgaox[l] + gy * kgaoy[l] + gz * kgaoz[l]);
                    
                    fvxc += fgrd * kgaos[l] * (gx * bgaox[l] + gy * bgaoy[l] + gz * bgaoz[l]);
                }
                
                xcBuffer.at(idx) = fvxc;
                
                idx++;
            }
        }
        
        // distribute block of Kohn-Sham matrix elements
        
        #pragma omp critical
        {
            idx = 0;
            
            for (int32_t j = curao; j < (curao + bufdim); j++)
            {
                for (int32_t k = j; k < naos; k++)
                {
                    auto fvxc = (j == k) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                    
                    ksmat[j * naos + k] += fvxc;
                    
                    idx++;
                }
            }
        }
        
        // update AOs counter
        
        curao += bufdim;
    }
    
    // loop over last block
    
    if (naos % bufdim != 0)
    {
        // compute last block of Kohn-Sham matrix elements
        
        int32_t idx = 0;
        
        for (int32_t i = curao; i < naos; i++)
        {
            auto bgaos = gtoValues.data(i);
            
            auto bgaox = gtoValuesX.data(i);
            
            auto bgaoy = gtoValuesY.data(i);
            
            auto bgaoz = gtoValuesZ.data(i);
            
            for (int32_t j = i; j < naos; j++)
            {
                auto kgaos = gtoValues.data(j);
                
                auto kgaox = gtoValuesX.data(j);
                
                auto kgaoy = gtoValuesY.data(j);
                
                auto kgaoz = gtoValuesZ.data(j);
                
                double fvxc = 0.0;
                
                auto koff = gridOffset + gridBlockPosition;
                
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    double w = gridWeights[koff + k];
                    
                    double gx = gradax[koff + k];
                    
                    double gy = graday[koff + k];
                    
                    double gz = gradaz[koff + k];
                    
                    fvxc += w * bgaos[k] * kgaos[k] * grhoa[koff + k];
                    
                    double fgrd = w * (ggrada[koff + k] / ngrada[koff + k] + ggradab[koff + k]);
                    
                    fvxc += fgrd * bgaos[k] * (gx * kgaox[k] + gy * kgaoy[k] + gz * kgaoz[k]);
                    
                    fvxc += fgrd * kgaos[k] * (gx * bgaox[k] + gy * bgaoy[k] + gz * bgaoz[k]);
                }
                
                xcBuffer.at(idx) = fvxc;
                
                idx++;
            }
        }
        
        // distribute last block of Kohn-Sham matrix elements
        
        #pragma omp critical
        {
            idx = 0;
            
            for (int32_t i = curao; i < naos; i++)
            {
                for (int32_t j = i ; j < naos; j++)
                {
                    auto fvxc = (i == j ) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                    
                    ksmat[i * naos + j] += fvxc;
                    
                    idx++;
                }
            }
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGGA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const CXCGradientGrid*     xcGradientGrid,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CDensityGrid*        gsDensityGrid,
                                          const CDensityGrid*        rwDensityGrid,
                                          const CMemBlock2D<double>& gtoValues,
                                          const CMemBlock2D<double>& gtoValuesX,
                                          const CMemBlock2D<double>& gtoValuesY,
                                          const CMemBlock2D<double>& gtoValuesZ,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{
    // set up AOs blocks
    
    auto naos = aoKohnShamMatrix->getNumberOfRows();
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    auto nblock = naos / bufdim;
    
    // set up pointers to exchange-correlation functional derrivatives
    
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
    
    // set up loop over number of matrices
    
    auto nmat = aoKohnShamMatrix->getNumberOfMatrices();
    
    for (int32_t i = 0; i < nmat; i++)
    {
        // set up pointer to Kohn-Sham matrix
        
        auto ksmat = aoKohnShamMatrix->getMatrix(i);
        
        // set up pointers to perturbed density gradient norms
        
        auto rhowa = rwDensityGrid->alphaDensity(i);
        
        auto rhowb = rwDensityGrid->betaDensity(i);
        
        auto gradwa_x = rwDensityGrid->alphaDensityGradientX(i);
        
        auto gradwa_y = rwDensityGrid->alphaDensityGradientY(i);
        
        auto gradwa_z = rwDensityGrid->alphaDensityGradientZ(i);
        
        auto gradwb_x = rwDensityGrid->betaDensityGradientX(i);
        
        auto gradwb_y = rwDensityGrid->betaDensityGradientY(i);
        
        auto gradwb_z = rwDensityGrid->betaDensityGradientZ(i);
        
        // loop over AOs blocks
        
        int32_t curao = 0;
        
        for (int32_t j = 0; j < nblock; j++)
        {
            // compute block of Kohn-Sham matrix elements
            
            int32_t idx = 0;
            
            for (int32_t k = curao; k < (curao + bufdim); k++)
            {
                auto bgaos = gtoValues.data(k);
                
                auto bgaox = gtoValuesX.data(k);
                
                auto bgaoy = gtoValuesY.data(k);
                
                auto bgaoz = gtoValuesZ.data(k);
                
                for (int32_t l = k ; l < naos; l++)
                {
                    auto kgaos = gtoValues.data(l);
                    
                    auto kgaox = gtoValuesX.data(l);
                    
                    auto kgaoy = gtoValuesY.data(l);
                    
                    auto kgaoz = gtoValuesZ.data(l);
                    
                    double fvxc = 0.0;
                    
                    auto moff = gridOffset + gridBlockPosition;
                    
                    for (int32_t m = 0; m < nGridPoints; m++)
                    {
                        double w = gridWeights[moff + m];
                        
                        double znva = 1.0 / ngrada[moff + m];
                        
                        double znvb = 1.0 / ngradb[moff + m];
                        
                        double rxa = znva * grada_x[moff + m];
                        
                        double rya = znva * grada_y[moff + m];
                        
                        double rza = znva * grada_z[moff + m];
                        
                        double rxb = znvb * gradb_x[moff + m];
                        
                        double ryb = znvb * gradb_y[moff + m];
                        
                        double rzb = znvb * gradb_z[moff + m];
                        
                        double rxwa = gradwa_x[moff + m];
                        
                        double rywa = gradwa_y[moff + m];
                        
                        double rzwa = gradwa_z[moff + m];
                        
                        double rxwb = gradwb_x[moff + m];
                        
                        double rywb = gradwb_y[moff + m];
                        
                        double rzwb = gradwb_z[moff + m];
                        
                        // GTOs values
                        
                        auto a0 = bgaos[m] * kgaos[m];
                        
                        auto ax = bgaox[m] * kgaos[m] + bgaos[m] * kgaox[m];
                        
                        auto ay = bgaoy[m] * kgaos[m] + bgaos[m] * kgaoy[m];
                        
                        auto az = bgaoz[m] * kgaos[m] + bgaos[m] * kgaoz[m];
                       
                        //  variations of functionals variables
                        
                        double zetaa = rxwa * rxa + rywa * rya + rzwa * rza;
                        
                        double zetab = rxwb * rxb + rywb * ryb + rzwb * rzb;
                        
                        double zetac = grada_x[moff + m] * rxwb + grada_y[moff + m] * rywb
                        
                                     + grada_z[moff + m] * rzwb + gradb_x[moff + m] * rxwa
                        
                                     + gradb_y[moff + m] * rywa + gradb_z[moff + m] * rzwa;
                        
                        // first contribution
                        
                        double fac0 = gmix_aa[moff + m] * zetaa + gmix_ab[moff + m] * zetab
                        
                                    + gmix_ac[moff + m] * zetac + grho_aa[moff + m] * rhowa[moff + m]
                        
                                    + grho_ab[moff + m] * rhowb[moff + m];
                        
                        fvxc += w * a0 * fac0;
                        
                        // second contribution
                        
                        double facr = gmix_aa[moff + m] * rhowa[moff + m] + gmix_ab[moff + m] * rhowb[moff + m]
                        
                                    + ggrad_aa[moff + m] * zetaa + ggrad_ab[moff + m] * zetab + ggrad_ac[moff + m] * zetac;
                        
                        double ar = ax * rxa + ay * rya + az * rza;
                        
                        fvxc += w * facr * ar;
                        
                        // third contribution
                        
                        double facz = gmix_ac[moff + m] * rhowa[moff + m] +  gmix_bc[moff + m] * rhowb[moff + m]
                        
                                    + ggrad_ac[moff + m] * zetaa + ggrad_bc[moff + m] * zetab + ggrad_cc[moff + m] * zetac;
                        
                        double arb = ax * grada_x[moff + m] + ay * grada_y[moff + m] + az * grada_z[moff + m];
                        
                        fvxc += w * facz * arb;
                        
                        // fourth contribution
                        
                        double ab = ax * rxwa + ay * rywa + az * rzwa - ar * zetaa;
                        
                        fvxc += w * znva * ggrad_a[moff + m] * ab;
                        
                        // fifth contribution
                        
                        double abw = ax * rxwa + ay * rywa + az * rzwa;
                        
                        fvxc += w * ggrad_c[moff + m] * abw;
                    }
                    
                    xcBuffer.at(idx) = fvxc;
                    
                    idx++;
                }
            }
            
            // distribute block of Kohn-Sham matrix elements
            
            #pragma omp critical
            {
                idx = 0;
                
                for (int32_t k = curao; k < (curao + bufdim); k++)
                {
                    for (int32_t l = k; l < naos; l++)
                    {
                        auto fvxc = (k == l) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                        
                        ksmat[k * naos + l] += fvxc;
                        
                        idx++;
                    }
                }
            }
            
            // update AOs counter
            
            curao += bufdim;
        }
        
        // loop over last block
        
        if (naos % bufdim != 0)
        {
            // compute last block of Kohn-Sham matrix elements
            
            int32_t idx = 0;
            
            for (int32_t j = curao; j < naos; j++)
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
                    
                    double fvxc = 0.0;
                    
                    auto loff = gridOffset + gridBlockPosition;
                    
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        double w = gridWeights[loff + l];
                        
                        double znva = 1.0 / ngrada[loff + l];
                        
                        double znvb = 1.0 / ngradb[loff + l];
                        
                        double rxa = znva * grada_x[loff + l];
                        
                        double rya = znva * grada_y[loff + l];
                        
                        double rza = znva * grada_z[loff + l];
                        
                        double rxb = znvb * gradb_x[loff + l];
                        
                        double ryb = znvb * gradb_y[loff + l];
                        
                        double rzb = znvb * gradb_z[loff + l];
                        
                        double rxwa = gradwa_x[loff + l];
                        
                        double rywa = gradwa_y[loff + l];
                        
                        double rzwa = gradwa_z[loff + l];
                        
                        double rxwb = gradwb_x[loff + l];
                        
                        double rywb = gradwb_y[loff + l];
                        
                        double rzwb = gradwb_z[loff + l];
                        
                        // GTOs values
                        
                        auto a0 = bgaos[l] * kgaos[l];
                        
                        auto ax = bgaox[l] * kgaos[l] + bgaos[l] * kgaox[l];
                        
                        auto ay = bgaoy[l] * kgaos[l] + bgaos[l] * kgaoy[l];
                        
                        auto az = bgaoz[l] * kgaos[l] + bgaos[l] * kgaoz[l];
                        
                        //  variations of functionals variables
                        
                        double zetaa = rxwa * rxa + rywa * rya + rzwa * rza;
                        
                        double zetab = rxwb * rxb + rywb * ryb + rzwb * rzb;
                        
                        double zetac = grada_x[loff + l] * rxwb + grada_y[loff + l] * rywb
                        
                                     + grada_z[loff + l] * rzwb + gradb_x[loff + l] * rxwa
                        
                                     + gradb_y[loff + l] * rywa + gradb_z[loff + l] * rzwa;
                        
                        // first contribution
                        
                        double fac0 = gmix_aa[loff + l] * zetaa + gmix_ab[loff + l] * zetab
                        
                                    + gmix_ac[loff + l] * zetac + grho_aa[loff + l] * rhowa[loff + l]
                        
                                    + grho_ab[loff + l] * rhowb[loff + l];
                        
                        fvxc += w * a0 * fac0;
                        
                        // second contribution
                        
                        double facr = gmix_aa[loff + l] * rhowa[loff + l] + gmix_ab[loff + l] * rhowb[loff + l]
                        
                                    + ggrad_aa[loff + l] * zetaa + ggrad_ab[loff + l] * zetab + ggrad_ac[loff + l] * zetac;
                        
                        double ar = ax * rxa + ay * rya + az * rza;
                        
                        fvxc += w * facr * ar;
                        
                        // third contribution
                        
                        double facz = gmix_ac[loff + l] * rhowa[loff + l] + gmix_bc[loff + l] * rhowb[loff + l]
                        
                                    + ggrad_ac[loff + l] * zetaa + ggrad_bc[loff + l] * zetab + ggrad_cc[loff + l] * zetac;
                        
                        double arb = ax * grada_x[loff + l] + ay * grada_y[loff + l] + az * grada_z[loff + l];
                        
                        fvxc += w * facz * arb;
                        
                        // fourth contribution
                        
                        double ab = ax * rxwa + ay * rywa + az * rzwa - ar * zetaa;
                        
                        fvxc += w * znva * ggrad_a[loff + l] * ab;
                        
                        // fifth contribution
                        
                        double abw = ax * rxwa + ay * rywa + az * rzwa;
                        
                        fvxc += w * ggrad_c[loff + l] * abw;
                    }
                    
                    xcBuffer.at(idx) = fvxc;
                    
                    idx++;
                }
            }
            
            // distribute last block of Kohn-Sham matrix elements
            
            #pragma omp critical
            {
                idx = 0;
                
                for (int32_t j = curao; j < naos; j++)
                {
                    for (int32_t k = j ; k < naos; k++)
                    {
                        auto fvxc = (j == k ) ? 0.5 * xcBuffer.at(idx) : xcBuffer.at(idx);
                        
                        ksmat[j * naos + k] += fvxc;
                        
                        idx++;
                    }
                }
            }
        }
    }
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
                    
                    //_compRestrictedEnergyForBatchOfGridPoints(xcele, xcene, xcgridptr, drwptr, mgw, tbposition, tbsize);
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

int32_t
CXCIntegrator::_getSizeOfBlock() const
{
    return 500;
}

int32_t
CXCIntegrator::_getNumberOfAOsInBuffer() const
{
    return 5;
}

void
CXCIntegrator::_addSubMatrix(      CAOKohnShamMatrix* aoKohnShamMatrix,
                             const CDenseMatrix&      subMatrix,
                             const CGtoBlock&         braGtoBlock,
                             const CGtoBlock&         ketGtoBlock) const
{
    // set up dimensions of Kohn-Sham matrix
    
    auto ncols = aoKohnShamMatrix->getNumberOfColumns();
    
    // set up submatrix dimensions

    auto bdim = subMatrix.getNumberOfRows();
    
    auto kdim = subMatrix.getNumberOfColumns();
    
    // set up submatrix offsets
    
    auto boff = (braGtoBlock.getIdentifiers(0))[0];
    
    auto koff = (ketGtoBlock.getIdentifiers(0))[0];
    
    // set up pointers to data
    
    auto ksmat = aoKohnShamMatrix->getMatrix(0);
    
    auto bkmat = subMatrix.values();
    
    // loop over submatrix
    
    for (int32_t i = 0; i < bdim; i++)
    {
        auto ioff = (boff + i) * ncols;
        
        for (int32_t j = 0; j < kdim; j++)
        {
            ksmat[ioff + koff + j] += bkmat[i * kdim + j];
        }
    }
}

CMemBlock<int32_t>
CXCIntegrator::_getScreeningPattern(const CMemBlock2D<double>& gtoValues) const
{
    auto ndim = gtoValues.size(0);
    
    // determine max values vector
    
    CMemBlock<double> maxvals(ndim);
    
    auto mvalsptr = maxvals.data();
    
    mathfunc::copy(mvalsptr, 0, gtoValues.data(0), 0, ndim);
    
    for (int32_t i = 0; i < ndim; i++) mvalsptr[i] = std::fabs(mvalsptr[i]);
    
    for (int32_t i = 1; i < gtoValues.blocks(); i++)
    {
        auto curvals = gtoValues.data(i);
        
        for (int32_t j = 0; j < ndim; j++)
        {
            auto fval = std::fabs(curvals[j]);
            
            if (fval > mvalsptr[j]) mvalsptr[j] = fval;
        }
    }
    
    // determine screening pattern
    
    CMemBlock<int32_t> patids(ndim);
    
    int32_t npoints = 0;
    
    for (int32_t i = 0; i < ndim; i++)
    {
        if (mvalsptr[i] > _thresholdOfDensity)
        {
            patids.at(npoints) = i;
            
            npoints++;
        }
    }
    
    if (npoints > 0) return patids.slice(0, npoints);
    
    return CMemBlock<int32_t>();
}

CMemBlock2D<double>
CXCIntegrator::_getReducedGrid(const double*             gridCoordinatesX,
                               const double*             gridCoordinatesY,
                               const double*             gridCoordinatesZ,
                               const double*             gridWeights,
                               const CMemBlock<int32_t>& screeningPattern) const
{
    auto ndim = screeningPattern.size();
    
    CMemBlock2D<double> mgrid(ndim, 4);
    
    // set up pointers to molecular grid
    
    auto mx = mgrid.data(0);
    
    auto my = mgrid.data(1);
    
    auto mz = mgrid.data(2);
    
    auto mw = mgrid.data(3);
    
    for (int32_t i = 0; i < ndim; i++)
    {
        auto idx = screeningPattern.at(i);
        
        mx[i] = gridCoordinatesX[idx];
        
        my[i] = gridCoordinatesY[idx];
        
        mz[i] = gridCoordinatesZ[idx];
        
        mw[i] = gridWeights[idx];
    }
    
    return mgrid;
}

CMemBlock2D<double>
CXCIntegrator::_getReducedGtoValues(const CMemBlock2D<double>& gtoValues,
                                    const CMemBlock<int32_t>&  screeningPattern) const
{
    auto ndim = screeningPattern.size();
    
    auto nblk = gtoValues.blocks();
    
    CMemBlock2D<double> gvals(ndim, nblk);
    
    auto pids = screeningPattern.data();
    
    for (int32_t i = 0; i < nblk; i++)
    {
        auto curvals = gtoValues.data(i);
        
        auto redvals = gvals.data(i);
        
        for (int32_t j = 0; j < ndim; j++)
        {
            redvals[j] = curvals[pids[j]];
        }
    }
    
    return gvals;
}

CMemBlock<double>
CXCIntegrator::_getReducedRestrictedGradient(const CXCGradientGrid*    xcGradientGrid,
                                             const CMemBlock<int32_t>& screeningPattern) const
{
    auto ndim = screeningPattern.size();
    
    CMemBlock<double> ggrid(ndim);
    
    auto redgrad = ggrid.data();
    
    auto curgrad = xcGradientGrid->xcGradientValues(xcvars::rhoa);
    
    auto pids = screeningPattern.data();
    
    for (int32_t i = 0; i < ndim; i++)
    {
        redgrad[i] = curgrad[pids[i]]; 
    }
    
    return ggrid;
}
