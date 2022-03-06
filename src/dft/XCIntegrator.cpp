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

#include "XCIntegrator.hpp"

#include <cmath>
#include <iostream> 

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

    if (aoDensityMatrix.getNumberOfDensityMatrices() == 1)
    {   
        // parse exchange-correlation functional data

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
        
        // create GTOs container
        
        CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
        
        // generate reference density grid
        
        CDensityGridDriver dgdrv(_locComm);
        
        auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());

        if (aoDensityMatrix.isClosedShell())
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
                _compRestrictedContributionForLda(ksmat, gtovec, vxcgrid, dgrid, mgrid);
            }
            
            if (fvxc.getFunctionalType() == xcfun::gga)
            {
                _compRestrictedContributionForGga(ksmat, gtovec, vxcgrid, dgrid, mgrid);
            }

            if (fvxc.getFunctionalType() == xcfun::mgga)
            {
                _compRestrictedContributionForMgga(ksmat, gtovec, vxcgrid, dgrid, mgrid);
            }
            
            // symmetrize Kohn-Sham matrix
            
            ksmat.symmetrize();
        }
        else
        {
            // generate screened molecular and density grids

            CMolecularGrid mgridab(molecularGrid);

            CMolecularGrid mgrida(molecularGrid);
            
            CMolecularGrid mgridb(molecularGrid);

            CDensityGrid dgridab;

            CDensityGrid dgrida;
            
            CDensityGrid dgridb;

            // screen the molecular and density grids 

            refdengrid.getScreenedGridPairUnrestricted(dgridab, dgrida, dgridb, mgridab, mgrida, mgridb, 
                                                       0, _thresholdOfDensity, fvxc.getFunctionalType());
            
            // allocate XC gradient grid

            CXCGradientGrid vxcgridab(mgridab.getNumberOfGridPoints(), dgridab.getDensityGridType(), fvxc.getFunctionalType());

            CXCGradientGrid vxcgrida(mgrida.getNumberOfGridPoints(), dgrida.getDensityGridType(), fvxc.getFunctionalType());

            CXCGradientGrid vxcgridb(mgridb.getNumberOfGridPoints(), dgridb.getDensityGridType(), fvxc.getFunctionalType());
            
            // compute exchange-correlation functional first derrivatives
            
            fvxc.compute(vxcgridab, dgridab);

            fvxc.compute(vxcgrida, dgrida);

            fvxc.compute(vxcgridb, dgridb);

            // compute Kohn-Sham matrix

            ksmat = CAOKohnShamMatrix(aoDensityMatrix.getNumberOfRows(0), aoDensityMatrix.getNumberOfColumns(0), false);
            
            ksmat.zero();
            
            if (fvxc.getFunctionalType() == xcfun::lda)
            {
               _compUnrestrictedContributionForLda(ksmat, gtovec, vxcgrida, dgrida, mgrida);

               _compUnrestrictedContributionForLda(ksmat, gtovec, vxcgridb, dgridb, mgridb);

               _compUnrestrictedContributionForLda(ksmat, gtovec, vxcgridab, dgridab, mgridab);

               auto xcdata = _compEnergyAndDensityUnrestricted(vxcgrida, dgrida, mgrida);

               auto xcdatb = _compEnergyAndDensityUnrestricted(vxcgridb, dgridb, mgridb);

               auto xcdatab = _compEnergyAndDensityUnrestricted(vxcgridab, dgridab, mgridab);
               
               ksmat.setExchangeCorrelationEnergy(std::get<0>(xcdata) + std::get<0>(xcdatb) + std::get<0>(xcdatab));
    
               ksmat.setNumberOfElectrons(std::get<1>(xcdata) + std::get<1>(xcdatb) + std::get<1>(xcdatab));
            }

            if (fvxc.getFunctionalType() == xcfun::gga)
            {
                _compUnrestrictedContributionForGga(ksmat, gtovec, vxcgridab, dgridab, mgridab);

                _compUnrestrictedContributionForGga(ksmat, gtovec, vxcgrida, dgrida, mgrida);

                _compUnrestrictedContributionForGga(ksmat, gtovec, vxcgridb, dgridb, mgridb);
            
               auto xcdata = _compEnergyAndDensityUnrestricted(vxcgrida, dgrida, mgrida);

               auto xcdatb = _compEnergyAndDensityUnrestricted(vxcgridb, dgridb, mgridb);

               auto xcdatab = _compEnergyAndDensityUnrestricted(vxcgridab, dgridab, mgridab);

               ksmat.setExchangeCorrelationEnergy(std::get<0>(xcdata) + std::get<0>(xcdatb) + std::get<0>(xcdatab));
    
               ksmat.setNumberOfElectrons(std::get<1>(xcdata) + std::get<1>(xcdatb) + std::get<1>(xcdatab));
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
    
    if (rwDensityMatrix.isClosedShell())
    {
        // allocate perturbed Kohn-Sham matrix
        
        CAOKohnShamMatrix ksmat(nrow, ncol, ndmat, true);
        
        if (fvxc.getFunctionalType() == xcfun::lda)
        {
            _compRestrictedContributionForLda(ksmat, gtovec, vxc2grid, rwdengrid, mgrid);
        }
        
        if (fvxc.getFunctionalType() == xcfun::gga)
        {
            _compRestrictedContributionForGga(ksmat, gtovec, vxcgrid, vxc2grid, gsdengrid, rwdengrid, mgrid);
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
CXCIntegrator::integrate(      CAOFockMatrix&    aoFockMatrix,
                         const CAODensityMatrix& rwDensityMatrix,
                         const CAODensityMatrix& rw2DensityMatrix,
                         const CAODensityMatrix& gsDensityMatrix,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel,
                         const std::string&      quadMode) const
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
    
    CXCCubicHessianGrid vxc3grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());  
    
    fvxc.compute(vxcgrid, gsdengrid);

    fvxc.compute(vxc2grid, gsdengrid);

    fvxc.compute(vxc3grid, gsdengrid);

    // compute perturbed density 
    
    auto rwdengrid = dgdrv.generate(rwDensityMatrix, molecule, basis, mgrid, fvxc.getFunctionalType());

    auto rw2dengrid = dgdrv.generate(rw2DensityMatrix, molecule, basis, mgrid, fvxc.getFunctionalType());

    auto rwdengridc = CDensityGridQuad(mgrid.getNumberOfGridPoints(), rw2DensityMatrix.getNumberOfDensityMatrices(), fvxc.getFunctionalType(), dengrid::ab);

    // Compute all and store all products of first-order transformed denisites 

    rwdengridc.DensityProd(rwdengridc,mgrid,rwdengrid,fvxc.getFunctionalType(),rw2DensityMatrix.getNumberOfDensityMatrices(),quadMode);
            
    // set up number of perturbed denstries and matrix dimensions
    
    auto ndmat = rw2DensityMatrix.getNumberOfDensityMatrices();
    
    auto nrow = gsDensityMatrix.getNumberOfRows(0);
    
    auto ncol = gsDensityMatrix.getNumberOfColumns(0);
    
    // allocate perturbed Kohn-Sham matrices
        
    CAOKohnShamMatrix ksmat(nrow, ncol, ndmat, true);

    if (fvxc.getFunctionalType() == xcfun::lda)
    {

        _compRestrictedContributionForLda(ksmat, gtovec, vxc2grid, vxc3grid, rwdengridc, rw2dengrid, mgrid,quadMode);
        
    }
        
    if (fvxc.getFunctionalType() == xcfun::gga)
    {
        _compRestrictedContributionForGga(ksmat, gtovec, vxcgrid, vxc2grid, vxc3grid, gsdengrid, rwdengridc, rw2dengrid, mgrid);
    }
    
    // symmetrize Kohn-Sham matrix
    
    ksmat.symmetrize(); 
    
    // add Kohn-Sham contribution to Fock matrix
    
    for (int32_t i = 0; i < ksmat.getNumberOfMatrices(); i++)
    {
        aoFockMatrix.addOneElectronMatrix(ksmat.getReferenceToMatrix(i), i); 

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
CXCIntegrator::_compRestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
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
                    _compRestrictedBatchForLda(ksmatprt, gtoContainer, xcgridptr, mgx, mgy, mgz, mgw,
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
CXCIntegrator::_compUnrestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                   const CGtoContainer*     gtoContainer,
                                                   const CXCGradientGrid&   xcGradientGrid,
                                                   const CDensityGrid&      densityGrid,
                                                   const CMolecularGrid&    molecularGrid) const
{      
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
                    _compUnrestrictedBatchForLda(ksmatprt, gtoContainer, xcgridptr, mgx, mgy, mgz, mgw,
                                               tbposition, tbsize);
                }
            }
        }
    }    
}

void
CXCIntegrator::_compRestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
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
                    _compRestrictedBatchForLda(ksmatprt, gtoContainer, xchessptr, rwdenptr, mgx, mgy, mgz, mgw,
                                               tbposition, tbsize);
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedContributionForLda(      CAOKohnShamMatrix&   aoKohnShamMatrix,
                                                 const CGtoContainer*       gtoContainer,
                                                 const CXCHessianGrid&      xcHessianGrid,
                                                 const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                 const CDensityGridQuad&    rwDensityGrid,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const CMolecularGrid&      molecularGrid,
                                                 const std::string&         quadMode) const
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

    auto xchess2ptr = &xcCubicHessianGrid;
    
    // set up pointer to perturbed density grid
    
    auto rwdenptr = &rwDensityGrid;

    auto rw2denptr = &rw2DensityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // compute contributions to Kohn-Sham matrix from blocks of grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, xchessptr, xchess2ptr, rwdenptr,rw2denptr, ksmatprt)
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
                    _compRestrictedBatchForLda(ksmatprt, gtoContainer, xchessptr, xchess2ptr, rwdenptr, rw2denptr, mgx, mgy, mgz, mgw,
                                               tbposition, tbsize,quadMode);
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
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
CXCIntegrator::_compRestrictedContributionForMgga(     CAOKohnShamMatrix& aoKohnShamMatrix,
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
                            _compRestrictedBatchForMGGA(ksmatprt, bgtos, kgtos, xcgridptr, dgridptr,
                                                       mgx, mgy, mgz, mgw, tbposition, tbsize);
                        }
                    }
                }
            }
        }
    }
    
    // compute exchange-correlation energy and number of electrons
    
    auto xcdat = _compEnergyAndDensity(xcGradientGrid, densityGrid, molecularGrid);

    std::cout << "xc energy " << std::get<0>(xcdat) << std::endl;
    
    aoKohnShamMatrix.setExchangeCorrelationEnergy(std::get<0>(xcdat));
    
    aoKohnShamMatrix.setNumberOfElectrons(std::get<1>(xcdat));
}

void
CXCIntegrator::_compUnrestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                   const CGtoContainer*     gtoContainer,
                                                   const CXCGradientGrid&   xcGradientGrid,
                                                   const CDensityGrid&      densityGrid,
                                                   const CMolecularGrid&    molecularGrid) const
{
        
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
                            _compUnrestrictedBatchForGga(ksmatprt, bgtos, kgtos, xcgridptr, dgridptr,
                                                         mgx, mgy, mgz, mgw, tbposition, tbsize);
                        }
                    }
                }
            }
        }
    }  
}

void
CXCIntegrator::_compRestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
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
CXCIntegrator::_compRestrictedContributionForGga(      CAOKohnShamMatrix&   aoKohnShamMatrix,
                                                 const CGtoContainer*       gtoContainer,
                                                 const CXCGradientGrid&     xcGradientGrid,
                                                 const CXCHessianGrid&      xcHessianGrid,
                                                 const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                 const CDensityGrid&        gsDensityGrid,
                                                 const CDensityGridQuad&    rwDensityGrid,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const CMolecularGrid&      molecularGrid) const
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

    auto xchess2ptr = &xcCubicHessianGrid;
    
    // set up pointer to density grid
    
    auto gsdenptr = &gsDensityGrid;
    
    auto rwdenptr = &rwDensityGrid;

    auto rw2denptr = &rw2DensityGrid;
    
    // set up pointer to Kohn-Sham matrix
    
    auto ksmatprt = &aoKohnShamMatrix;
    
    // compute contributions to Kohn-Sham matrix from blocks of grid points
    
    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, mgw, xcgradptr, xchessptr,xchess2ptr,\
                                gsdenptr, rwdenptr, rw2DensityGrid,ksmatprt)
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
                    _compRestrictedBatchForGGA(ksmatprt, gtoContainer, xcgradptr, xchessptr,xchess2ptr, gsdenptr, rwdenptr,rw2denptr,
                                               mgx, mgy, mgz, mgw, tbposition, tbsize);
                }
            }
        }
    }
}

void
CXCIntegrator::_compRestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
            
            _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf, grhoa, gaos, gridWeights, gridOffset, igpnt, blockdim);
            
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
        
        _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf, grhoa, gaos, gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compUnrestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    CMemBlock<double> vxcbuf_alpha(bufdim * naos);

    CMemBlock<double> vxcbuf_beta(bufdim * naos);
    
    int32_t igpnt = 0;
    
    if (xcGradientGrid->getDensityGridType() == dengrid::ab)
    {

        auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);

        auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

        if (nblocks > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            for (int32_t i = 0; i < nblocks; i++)
            {
                gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                                gridOffset, igpnt, blockdim);
                
                _distUnrestrictedBatchForLda(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa, grhob, gaos, gridWeights, gridOffset, igpnt, blockdim);
                
                igpnt += blockdim;
            }
        }
    
        // compute remaining grid points block
        
        blockdim = nGridPoints % blockdim;
        
        if (blockdim > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                            gridOffset, igpnt, blockdim);
            
            _distUnrestrictedBatchForLda(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa,grhob, gaos, gridWeights, gridOffset, igpnt, blockdim);
            
        }
    }

    if (xcGradientGrid->getDensityGridType() == dengrid::lima)
    {
        auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

        auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);

        if (nblocks > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            for (int32_t i = 0; i < nblocks; i++)
            {
                gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                                gridOffset, igpnt, blockdim);
                
                _distUnrestrictedBatchForLdaA(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa, grhob, gaos,
                                                 gridWeights, gridOffset, igpnt, blockdim);
                
                igpnt += blockdim;
            }
        }
        
        // compute remaining grid points block
        
        blockdim = nGridPoints % blockdim;
        
        if (blockdim > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                            gridOffset, igpnt, blockdim);
            
            _distUnrestrictedBatchForLdaA(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa,grhob, gaos,
                                             gridWeights, gridOffset, igpnt, blockdim);
            
        }
    }
    
    if (xcGradientGrid->getDensityGridType() == dengrid::limb)
    {
        auto grhoa = xcGradientGrid->xcGradientValues(xcvars::rhoa);

        auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);
        
        if (nblocks > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            for (int32_t i = 0; i < nblocks; i++)
            {
                gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                                gridOffset, igpnt, blockdim);
                
                _distUnrestrictedBatchForLdaB(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa, grhob, gaos, gridWeights, gridOffset, igpnt, blockdim);
                
                igpnt += blockdim;
            }
        }
        
        // compute remaining grid points block
        
        blockdim = nGridPoints % blockdim;
        
        if (blockdim > 0)
        {
            CMemBlock2D<double> gaos(blockdim, naos);
            
            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                            gridOffset, igpnt, blockdim);
            
            _distUnrestrictedBatchForLdaB(aoKohnShamMatrix, vxcbuf_alpha,vxcbuf_beta, grhoa,grhob, gaos, gridWeights, gridOffset, igpnt, blockdim);
            
        }
    }
}

void
CXCIntegrator::_compRestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
            
            _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf, xcHessianGrid, rwDensityGrid, gaos,
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
        
        _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf, xcHessianGrid, rwDensityGrid, gaos, gridWeights,
                                   gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                          const CGtoContainer*       gtoContainer,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CXCCubicHessianGrid* xcCubicHessianGrid,
                                          const CDensityGridQuad*        rwDensityGrid,
                                          const CDensityGrid*        rw2DensityGrid,
                                          const double*              gridCoordinatesX,
                                          const double*              gridCoordinatesY,
                                          const double*              gridCoordinatesZ,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              nGridPoints,
                                          const std::string&         quadMode) const
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
                                        

            _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf,xcHessianGrid, xcCubicHessianGrid, rwDensityGrid,rw2DensityGrid, gaos,
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
        
        
        _distRestrictedBatchForLda(aoKohnShamMatrix, vxcbuf,xcHessianGrid, xcCubicHessianGrid, rwDensityGrid,rw2DensityGrid, gaos,
                            gridWeights, gridOffset, igpnt, blockdim);
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
            
            _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcgrad, dengrid, gaos, gaox, gaoy, gaoz,
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
        
        _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcgrad, dengrid, gaos, gaox, gaoy, gaoz,
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
                _distRestrictedBatchForGga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
                
                _distRestrictedBatchForGga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                           blockdim);
            }
            
            igpnt += blockdim;
        }
    }
    
    // compute remaining grid points block
    
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
            _distRestrictedBatchForGga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       gridWeights, gridOffset, igpnt, blockdim);
        }
        else
        {
            gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedBatchForGga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                       blockdim);
        }
    }
    
    // add submatrix to Kohn-Sham matrix
    
    #pragma omp critical
    _addSubMatrix(aoKohnShamMatrix, submat, braGtoBlock, ketGtoBlock);
}

void
CXCIntegrator::_compRestrictedBatchForMGGA(     CAOKohnShamMatrix* aoKohnShamMatrix,
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
            gtorec::computeGtosValuesForMGGA(baos, baox, baoy, baoz, braGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            if (symbk)
            {
                _distRestrictedBatchForMgga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForMGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
                
                _distRestrictedBatchForMgga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                           kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                           blockdim);
            }
            
            igpnt += blockdim;
        }
    }
    
    // compute remaining grid points block
    
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
        
        gtorec::computeGtosValuesForMGGA(baos, baox, baoy, baoz, braGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                        gridCoordinatesZ, gridOffset, igpnt, blockdim);
        
        if (symbk)
        {
            _distRestrictedBatchForMgga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       gridWeights, gridOffset, igpnt, blockdim);
        }
        else
        {
            gtorec::computeGtosValuesForMGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                            gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
            _distRestrictedBatchForMgga(submat, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                       kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                       blockdim);
        }
    }
    
    // add submatrix to Kohn-Sham matrix
    
    #pragma omp critical
    _addSubMatrix(aoKohnShamMatrix, submat, braGtoBlock, ketGtoBlock);
}

void
CXCIntegrator::_compUnrestrictedBatchForGga(       CAOKohnShamMatrix* aoKohnShamMatrix,
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
    
    CDenseMatrix submata(bdim, kdim);

    CDenseMatrix submatb(bdim,kdim);

    submata.zero(); 

    submatb.zero();
    
    // determine number of grid blocks
    
    auto blockdim = _getSizeOfBlock();
    
    auto nblocks = nGridPoints / blockdim;
    
    // set up current grid point
    
    int32_t igpnt = 0;
    
    // loop over grid points blocks

    if (xcGradientGrid->getDensityGridType() == dengrid::ab)
    {
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
                    _distUnrestrictedBatchForGga(submata, submatb,xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                 gridWeights, gridOffset, igpnt, blockdim);

                }
                else
                {
                    gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                    gridCoordinatesZ, gridOffset, igpnt, blockdim);
                    
                    _distUnrestrictedBatchForGga(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                 kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                                 blockdim);

                }
                
                igpnt += blockdim;
            }
        }
    
        // compute remaining grid points block
    
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
                _distUnrestrictedBatchForGga(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                             gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
            
                _distUnrestrictedBatchForGga(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                             kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                             blockdim);
            }
        }
    }

    if (xcGradientGrid->getDensityGridType() == dengrid::lima)
    {
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
                    _distUnrestrictedBatchForGgaA(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                  gridWeights, gridOffset, igpnt, blockdim);
                }
                else
                {
                    gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                    gridCoordinatesZ, gridOffset, igpnt, blockdim);
                    
                    _distUnrestrictedBatchForGgaA(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                  kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                                  blockdim);
                }
                
                igpnt += blockdim;
            }
        }
    
        // compute remaining grid points block
    
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
                _distUnrestrictedBatchForGgaA(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                              gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
                
                _distUnrestrictedBatchForGgaA(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                              kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                              blockdim);
            }
        }
    }

    if (xcGradientGrid->getDensityGridType() == dengrid::limb)
    {
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
                    _distUnrestrictedBatchForGgaB(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                  gridWeights, gridOffset, igpnt, blockdim);
                }
                else
                {
                    gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                    gridCoordinatesZ, gridOffset, igpnt, blockdim);
                    
                    _distUnrestrictedBatchForGgaB(submata, submatb, xcGradientGrid, densityGrid, baos, baox, baoy, baoz, braGtoBlock,
                                                  kaos, kaox, kaoy, kaoz, ketGtoBlock, gridWeights, gridOffset, igpnt,
                                                  blockdim);
                }
                igpnt += blockdim;
            }
        }
        
        // compute remaining grid points block
        
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
                _distUnrestrictedBatchForGgaB(submata, submatb, xcGradientGrid, densityGrid, 
                                              baos, baox, baoy, baoz, braGtoBlock,gridWeights, gridOffset, igpnt, blockdim);
            }
            else
            {
                gtorec::computeGtosValuesForGGA(kaos, kaox, kaoy, kaoz, ketGtoBlock, gridCoordinatesX, gridCoordinatesY,
                                                gridCoordinatesZ, gridOffset, igpnt, blockdim);
                
                _distUnrestrictedBatchForGgaB(submata, submatb,xcGradientGrid, densityGrid, 
                                              baos, baox, baoy, baoz, braGtoBlock,kaos, kaox, kaoy, kaoz, ketGtoBlock, 
                                              gridWeights, gridOffset, igpnt,
                                              blockdim);
            }
        }
    }

    // add submatrix to Kohn-Sham matrix
    
    #pragma omp critical
    _addSubMatrixUnRest(aoKohnShamMatrix, submata,submatb, braGtoBlock, ketGtoBlock);
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
            
            _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, gsDensityGrid, rwDensityGrid,
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
        
        _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, gsDensityGrid, rwDensityGrid,
                                   gaos, gaox, gaoy, gaoz, gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_compRestrictedBatchForGGA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                          const CGtoContainer*       gtoContainer,
                                          const CXCGradientGrid*     xcGradientGrid,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CXCCubicHessianGrid* xcCubicHessianGrid,
                                          const CDensityGrid*        gsDensityGrid,
                                          const CDensityGridQuad*    rwDensityGrid,
                                          const CDensityGrid*        rw2DensityGrid,
                                          const double*              gridCoordinatesX,
                                          const double*              gridCoordinatesY,
                                          const double*              gridCoordinatesZ,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              nGridPoints) const
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
            
            _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, xcCubicHessianGrid, gsDensityGrid, rwDensityGrid,rw2DensityGrid,
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
        
        _distRestrictedBatchForGga(aoKohnShamMatrix, vxcbuf, xcGradientGrid, xcHessianGrid, xcCubicHessianGrid, gsDensityGrid, rwDensityGrid,rw2DensityGrid,
                                   gaos, gaox, gaoy, gaoz, gridWeights, gridOffset, igpnt, blockdim);
    }
}

void
CXCIntegrator::_distRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
CXCIntegrator::_distUnrestrictedBatchForLdaB(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                   CMemBlock<double>&   xcBuffera,
                                                   CMemBlock<double>&   xcBufferb,
                                             const double*              xcGradient_alpha,
                                             const double*              xcGradient_beta,
                                             const CMemBlock2D<double>& gtoValues,
                                             const double*              gridWeights,
                                             const int32_t              gridOffset,
                                             const int32_t              gridBlockPosition,
                                             const int32_t              nGridPoints) const
{
    // set up pointer to Kohn-Sham matrix
    
    auto ksmata = aoKohnShamMatrix->getKohnSham(false);

    auto ksmatb = aoKohnShamMatrix->getKohnSham(true);

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
                
                double fvxca = 0.0;

                double fvxcb = 0.0;
                
                auto loff = gridOffset + gridBlockPosition;
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    fvxca += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_alpha[loff + l];

                    fvxcb += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_beta[loff + l]; 
                }
                
                xcBuffera.at(idx) = fvxca;

                xcBufferb.at(idx) = fvxcb;
                
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
                    const auto fvxca = (j == k) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);

                    const auto fvxcb = (j == k) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);

                    ksmatb[j * naos + k] += fvxcb; 

                    ksmata[j * naos + k] += fvxca;
                    
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
                
                double fvxca = 0.0;

                double fvxcb = 0.0;
                
                auto koff = gridOffset + gridBlockPosition;
                
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    fvxca += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_alpha[koff + k];

                    fvxcb += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_beta[koff + k]; 
                }
                
                xcBuffera.at(idx) = fvxca;

                xcBufferb.at(idx) = fvxcb;
                
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
                    const auto fvxca = (i == j) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);

                    const auto fvxcb = (i == j) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);
                    
                    ksmata[i * naos + j] += fvxca;

                    ksmatb[i * naos + j] += fvxcb;
                    
                    idx++;
                }
            }
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForLdaA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                   CMemBlock<double>&   xcBuffera,
                                                   CMemBlock<double>&   xcBufferb,
                                             const double*              xcGradient_alpha,
                                             const double*              xcGradient_beta,
                                             const CMemBlock2D<double>& gtoValues,
                                             const double*              gridWeights,
                                             const int32_t              gridOffset,
                                             const int32_t              gridBlockPosition,
                                             const int32_t              nGridPoints) const
{
    // set up pointer to Kohn-Sham matrix
    
    auto ksmata = aoKohnShamMatrix->getKohnSham(false);

    auto ksmatb = aoKohnShamMatrix->getKohnSham(true);
    
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
                
                double fvxca = 0.0;

                double fvxcb = 0.0;
                
                auto loff = gridOffset + gridBlockPosition;
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    fvxcb += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_beta[loff + l];

                    fvxca += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_alpha[loff + l]; 
                }
                
                xcBuffera.at(idx) = fvxca;

                xcBufferb.at(idx) = fvxcb;
                
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
                    auto fvxcb = (j == k) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);

                    auto fvxca = (j == k) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);

                    ksmata[j * naos + k] += fvxca;

                    ksmatb[j * naos + k] += fvxcb;
                    
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
                
                double fvxcb = 0.0;

                double fvxca = 0.0;

                auto koff = gridOffset + gridBlockPosition;
                
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    fvxcb += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_beta[koff + k];
                    
                    fvxca += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_alpha[koff + k]; 
                }
                xcBufferb.at(idx) = fvxcb;

                xcBuffera.at(idx) = fvxca;
                
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
                    auto fvxcb = (i == j) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);
                    
                    auto fvxca = (i == j) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);
                    
                    ksmatb[i * naos + j] += fvxcb;

                    ksmata[i * naos + j] += fvxca;
                    
                    idx++;
                }
            }
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                  CMemBlock<double>&   xcBuffera,
                                                  CMemBlock<double>&   xcBufferb,
                                            const double*              xcGradient_alpha,
                                            const double*              xcGradient_beta,
                                            const CMemBlock2D<double>& gtoValues,
                                            const double*              gridWeights,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const
{
    // set up pointer to Kohn-Sham matrix
    
    auto ksmata = aoKohnShamMatrix->getKohnSham(false);

    auto ksmatb = aoKohnShamMatrix->getKohnSham(true);
    
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
                
                double fvxca = 0.0;

                double fvxcb = 0.0;
                
                auto loff = gridOffset + gridBlockPosition;
                
                for (int32_t l = 0; l < nGridPoints; l++)
                {
                    fvxca += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_alpha[loff + l];

                    fvxcb += bgaos[l] * kgaos[l] * gridWeights[loff + l] * xcGradient_beta[loff + l];
                }
                
                xcBuffera.at(idx) = fvxca;

                xcBufferb.at(idx) = fvxcb;
                
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
                    auto fvxca = (j == k) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);
                    
                    auto fvxcb = (j == k) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);
                    
                    ksmatb[j * naos + k] += fvxcb;
                
                    ksmata[j * naos + k] += fvxca;
                    
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
                
                double fvxca = 0.0;

                double fvxcb = 0.0;
                
                auto koff = gridOffset + gridBlockPosition;
                
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    fvxca += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_alpha[koff + k];

                    fvxcb += bgaos[k] * kgaos[k] * gridWeights[koff + k] * xcGradient_beta[koff + k];
                }
                
                xcBuffera.at(idx) = fvxca;

                xcBufferb.at(idx) = fvxcb;
                
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
                    auto fvxca = (i == j) ? 0.5 * xcBuffera.at(idx) : xcBuffera.at(idx);

                    auto fvxcb = (i == j) ? 0.5 * xcBufferb.at(idx) : xcBufferb.at(idx);
                    
                    ksmata[i * naos + j] += fvxca;

                    ksmatb[i * naos + j] += fvxcb;
                    
                    idx++;
                }
            }
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGga(      CDenseMatrix&        subMatrix,
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
CXCIntegrator::_distRestrictedBatchForMgga(      CDenseMatrix&        subMatrix,
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
    
    auto gtaua = xcGradientGrid->xcGradientValues(xcvars::taua);

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

                fvxc +=  w * 0.5 * gtaua[koff + k] * (bgaox[k] * kgaox[k] + bgaoy[k] * kgaoy[k] + bgaoz[k] * kgaoz[k]);
            }
            
            if (i == j) fvxc *= 0.5;
            
            bmat[i * bdim + j] += fvxc;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGga(      CDenseMatrix&        subMatrixa,
                                                  CDenseMatrix&        subMatrixb,
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

    auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

    auto ggrada = xcGradientGrid->xcGradientValues(xcvars::grada);

    auto ggradb = xcGradientGrid->xcGradientValues(xcvars::gradb);
    
    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms
    
    auto ngrada = densityGrid->alphaDensityGradient(0);
    
    auto gradax = densityGrid->alphaDensityGradientX(0);
    
    auto graday = densityGrid->alphaDensityGradientY(0);
    
    auto gradaz = densityGrid->alphaDensityGradientZ(0);

    auto ngradb = densityGrid->betaDensityGradient(0);
    
    auto gradbx = densityGrid->betaDensityGradientX(0);
    
    auto gradby = densityGrid->betaDensityGradientY(0);
    
    auto gradbz = densityGrid->betaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrixa.getNumberOfRows();
    
    // set up pointer to submatrix values
    
    auto bmat_alpha = subMatrixa.values();

    auto bmat_beta = subMatrixb.values();

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
                
            double fvxca = 0.0;

            double fvxcb = 0.0;
                
            auto koff = gridOffset + gridBlockPosition;
                
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];

                double gxa = gradax[koff + k];
                    
                double gya = graday[koff + k];
                    
                double gza = gradaz[koff + k];

                double gxb = gradbx[koff + k];
                    
                double gyb = gradby[koff + k];
                    
                double gzb = gradbz[koff + k];
                    
                fvxca += w * bgaos[k] * kgaos[k] * grhoa[koff + k];

                fvxcb += w * bgaos[k] * kgaos[k] * grhob[koff + k];
                    
                double gvala = bgaos[k] * (gxa * kgaox[k] + gya * kgaoy[k] + gza * kgaoz[k]) + kgaos[k] * (gxa * bgaox[k] + gya * bgaoy[k] + gza * bgaoz[k]);

                double gvalb = bgaos[k] * (gxb * kgaox[k] + gyb * kgaoy[k] + gzb * kgaoz[k]) + kgaos[k] * (gxb * bgaox[k] + gyb * bgaoy[k] + gzb * bgaoz[k]);

                fvxca += w * (ggrada[koff + k] * gvala/ ngrada[koff + k] + gvalb * ggradab[koff + k]);

                fvxcb += w * (ggradb[koff + k] *gvalb/ ngradb[koff + k] + gvala * ggradab[koff + k]);
             }
            
            if (i == j) 
            {
                fvxca *= 0.5;

                fvxcb *= 0.5;
            } 

            bmat_alpha[i * bdim + j] += fvxca;

            bmat_beta[i * bdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGgaA(      CDenseMatrix&        subMatrixa,
                                                   CDenseMatrix&        subMatrixb,
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

    auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

    auto ggradb = xcGradientGrid->xcGradientValues(xcvars::gradb);

    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
        
    // set up pointers to density gradient norms
    
    auto ngradb = densityGrid->betaDensityGradient(0);
    
    auto gradbx = densityGrid->betaDensityGradientX(0);
    
    auto gradby = densityGrid->betaDensityGradientY(0);
    
    auto gradbz = densityGrid->betaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrixb.getNumberOfRows();
    
    // set up pointer to submatrix values
    
    auto bmat_alpha = subMatrixa.values();

    auto bmat_beta = subMatrixb.values();

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
                
            double fvxca = 0.0;

            double fvxcb = 0.0;
                
            auto koff = gridOffset + gridBlockPosition;
                
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];

                double gxb = gradbx[koff + k];
                    
                double gyb = gradby[koff + k];
                    
                double gzb = gradbz[koff + k];
                
                fvxcb += w * bgaos[k] * kgaos[k] * grhob[koff + k];
                    
                double gvalb = bgaos[k] * (gxb * kgaox[k] + gyb * kgaoy[k] + gzb * kgaoz[k]) + kgaos[k] * (gxb * bgaox[k] + gyb * bgaoy[k] + gzb * bgaoz[k]);

                fvxca += w * (gvalb * ggradab[koff + k]);

                fvxcb += w * (ggradb[koff + k] *gvalb/ ngradb[koff + k]);
            }
            
            if (i == j) 
            {
                fvxca *= 0.5;

                fvxcb *= 0.5;
            } 

            bmat_alpha[i * bdim + j] += fvxca;

            bmat_beta[i * bdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGgaB(      CDenseMatrix&        subMatrixa,
                                                   CDenseMatrix&        subMatrixb,
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
    
    auto bdim = subMatrixa.getNumberOfRows();
    
    // set up pointer to submatrix values
    
    auto bmat_alpha = subMatrixa.values();

    auto bmat_beta = subMatrixb.values();

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
                
            double fvxca = 0.0;

            double fvxcb = 0.0;
                
            auto koff = gridOffset + gridBlockPosition;
                
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];

                double gxa = gradax[koff + k];
                    
                double gya = graday[koff + k];
                    
                double gza = gradaz[koff + k];

                fvxca += w * bgaos[k] * kgaos[k] * grhoa[koff + k];
                    
                double gvala = bgaos[k] * (gxa * kgaox[k] + gya * kgaoy[k] + gza * kgaoz[k]) + kgaos[k] * (gxa * bgaox[k] + gya * bgaoy[k] + gza * bgaoz[k]);

                fvxca += w * (ggrada[koff + k] * gvala/ ngrada[koff + k]);

                fvxcb += w * (gvala * ggradab[koff + k]);
            }
            
            if (i == j) 
            {
                fvxca *= 0.5;

                fvxcb *= 0.5;
            } 
            
            bmat_alpha[i * bdim + j] += fvxca;

            bmat_beta[i * bdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGga(      CDenseMatrix&        subMatrixa,
                                                  CDenseMatrix&        subMatrixb,
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

    auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

    auto ggrada = xcGradientGrid->xcGradientValues(xcvars::grada);

    auto ggradb = xcGradientGrid->xcGradientValues(xcvars::gradb);

    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms
    
    auto ngrada = densityGrid->alphaDensityGradient(0);
    
    auto gradax = densityGrid->alphaDensityGradientX(0);
    
    auto graday = densityGrid->alphaDensityGradientY(0);
    
    auto gradaz = densityGrid->alphaDensityGradientZ(0);

    auto ngradb = densityGrid->betaDensityGradient(0);
    
    auto gradbx = densityGrid->betaDensityGradientX(0);
    
    auto gradby = densityGrid->betaDensityGradientY(0);
    
    auto gradbz = densityGrid->betaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrixa.getNumberOfRows();
    
    auto kdim = subMatrixa.getNumberOfColumns();
    
    // set up pointer to submatrix values
    
    auto bkmata = subMatrixa.values();

    auto bkmatb = subMatrixb.values();
    
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
            
            double fvxca = 0.0;

            double fvxcb = 0.0;

            auto koff = gridOffset + gridBlockPosition;
            
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                double w = gridWeights[koff + k];

                double gxa = gradax[koff + k];
                
                double gya = graday[koff + k];
                
                double gza = gradaz[koff + k];

                double gxb = gradbx[koff + k];
                
                double gyb = gradby[koff + k];
                
                double gzb = gradbz[koff + k];

                fvxca += w * bgaos[k] * kgaos[k] * grhoa[koff + k];

                fvxcb += w * bgaos[k] * kgaos[k] * grhob[koff + k];
                    
                double gvala = bgaos[k] * (gxa * kgaox[k] + gya * kgaoy[k] + gza * kgaoz[k]) + kgaos[k] * (gxa * bgaox[k] + gya * bgaoy[k] + gza * bgaoz[k]);

                double gvalb = bgaos[k] * (gxb * kgaox[k] + gyb * kgaoy[k] + gzb * kgaoz[k]) + kgaos[k] * (gxb * bgaox[k] + gyb * bgaoy[k] + gzb * bgaoz[k]);

                fvxca += w * (ggrada[koff + k] * gvala/ ngrada[koff + k] + gvalb * ggradab[koff + k]);

                fvxcb += w * (ggradb[koff + k] *gvalb/ ngradb[koff + k] + gvala * ggradab[koff + k]);

            }
            
            bkmata[i * kdim + j] += fvxca;

            bkmatb[i * kdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGgaA(      CDenseMatrix&        subMatrixa,
                                                   CDenseMatrix&        subMatrixb,
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
    
    auto grhob = xcGradientGrid->xcGradientValues(xcvars::rhob);

    auto ggradb = xcGradientGrid->xcGradientValues(xcvars::gradb);

    auto ggradab = xcGradientGrid->xcGradientValues(xcvars::gradab);
    
    // set up pointers to density gradient norms

    auto ngradb = densityGrid->betaDensityGradient(0);
    
    auto gradbx = densityGrid->betaDensityGradientX(0);
    
    auto gradby = densityGrid->betaDensityGradientY(0);
    
    auto gradbz = densityGrid->betaDensityGradientZ(0);
    
    // set up dimensions of submatrix
    
    auto bdim = subMatrixa.getNumberOfRows();
    
    auto kdim = subMatrixa.getNumberOfColumns();
    
    // set up pointer to submatrix values
    
    auto bkmata = subMatrixa.values();

    auto bkmatb = subMatrixb.values();
    
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
            
            double fvxca = 0.0;

            double fvxcb = 0.0;

            auto koff = gridOffset + gridBlockPosition;
            
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                const double w = gridWeights[koff + k];
                
                const double gxb = gradbx[koff + k];
                
                const double gyb = gradby[koff + k];
                
                const double gzb = gradbz[koff + k];

                fvxcb += w * bgaos[k] * kgaos[k] * grhob[koff + k];
                    
                const double gvalb = bgaos[k] * (gxb * kgaox[k] + gyb * kgaoy[k] + gzb * kgaoz[k]) 
                
                                   + kgaos[k] * (gxb * bgaox[k] + gyb * bgaoy[k] + gzb * bgaoz[k]);

                fvxca += w * (gvalb * ggradab[koff + k]);

                fvxcb += w * (ggradb[koff + k] *gvalb/ ngradb[koff + k]);
            }
            
            bkmata[i * kdim + j] += fvxca;

            bkmatb[i * kdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distUnrestrictedBatchForGgaB(      CDenseMatrix&        subMatrixa,
                                                   CDenseMatrix&        subMatrixb,
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
    
    auto bdim = subMatrixa.getNumberOfRows();
    
    auto kdim = subMatrixa.getNumberOfColumns();
    
    // set up pointer to submatrix values
    
    auto bkmata = subMatrixa.values();

    auto bkmatb = subMatrixb.values();
    
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
            
            double fvxca = 0.0;

            double fvxcb = 0.0;

            auto koff = gridOffset + gridBlockPosition;
            
            for (int32_t k = 0; k < nGridPoints; k++)
            {
                const double w = gridWeights[koff + k];
                
                const double gxa = gradax[koff + k];
                
                const double gya = graday[koff + k];
                
                const double gza = gradaz[koff + k];

                fvxca += w * bgaos[k] * kgaos[k] * grhoa[koff + k];
                    
                const double gvala = bgaos[k] * (gxa * kgaox[k] + gya * kgaoy[k] + gza * kgaoz[k]) 
                          
                                   + kgaos[k] * (gxa * bgaox[k] + gya * bgaoy[k] + gza * bgaoz[k]);

                fvxca += w * (ggrada[koff + k] * gvala/ ngrada[koff + k]);

                fvxcb += w * (gvala * ggradab[koff + k]);
            }
            
            bkmata[i * kdim + j] += fvxca;

            bkmatb[i * kdim + j] += fvxcb;
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForGga(      CDenseMatrix&        subMatrix,
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
CXCIntegrator::_distRestrictedBatchForMgga(     CDenseMatrix&        subMatrix,
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

    auto gtaua = xcGradientGrid->xcGradientValues(xcvars::taua);
    
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

                fvxc +=  w * 0.5 * gtaua[koff + k] * (bgaox[k] * kgaox[k] + bgaoy[k] * kgaoy[k] + bgaoz[k] * kgaoz[k]);

            }
            
            bkmat[i * kdim + j] += fvxc;
        }
    }
}

void
CXCIntegrator::_distRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
CXCIntegrator::_distRestrictedBatchForLda(   CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CXCCubicHessianGrid* xcCubicHessianGrid,
                                          const CDensityGridQuad*    rwDensityGrid,
                                          const CDensityGrid*        rw2DensityGrid,
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
    
    auto grho_aaa = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);
    
    auto grho_aab = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

    auto grho_abb = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

    // set up loop over number of matrices
    
    auto nmat = aoKohnShamMatrix->getNumberOfMatrices();
    
    for (int32_t i = 0; i < nmat; i++)
    {
        // set up pointer to Kohn-Sham matrix

        auto ksmat = aoKohnShamMatrix->getMatrix(i);
        
        // set up pointer to perturbed density
        
        auto rhow1a = rwDensityGrid->rhow1rhow2(i);
        
        auto rhow12a = rw2DensityGrid->alphaDensity(i);
        
        auto rhow12b = rw2DensityGrid->betaDensity(i);

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
                       fvxc += bgaos[m] * kgaos[m] * gridWeights[moff + m] * 
                       
                            ( ( grho_aaa[moff + m] 
                             
                             +  grho_aab[moff + m] 

                             +  grho_aab[moff + m]
                             
                             +  grho_abb[moff + m] )* rhow1a[moff + m]
                             
                             +  grho_aa[moff + m] * rhow12a[moff + m] + grho_ab[moff + m] * rhow12b[moff + m] );                            
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
                        fvxc += bgaos[l] * kgaos[l] * gridWeights[loff + l] * 
                        
                            (  (grho_aaa[loff + l] 
                             
                             +  grho_aab[loff + l] 

                             +  grho_aab[loff + l] 
                             
                             +  grho_abb[loff + l]) * rhow1a[loff + l] 
                             
                             +  grho_aa[loff + l] * rhow12a[loff + l] + grho_ab[loff + l] * rhow12b[loff + l] ); 

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
CXCIntegrator::_distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
CXCIntegrator::_distRestrictedBatchForMgga(     CAOKohnShamMatrix*   aoKohnShamMatrix,
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

    auto gtaua = xcGradientGrid.xcGradientValues(xcvars::taua);

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

                    fvxc +=  w * 0.5 * gtaua[loff + l] * (bgaox[l] * kgaox[l] + bgaoy[l] * kgaoy[l] + bgaoz[l] * kgaoz[l]);

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

                    fvxc +=  w * 0.5 * gtaua[koff + k] * (bgaox[k] * kgaox[k] + bgaoy[k] * kgaoy[k] + bgaoz[k] * kgaoz[k]);

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
CXCIntegrator::_distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
CXCIntegrator::_distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer,
                                          const CXCGradientGrid*     xcGradientGrid,
                                          const CXCHessianGrid*      xcHessianGrid,
                                          const CXCCubicHessianGrid* xcCubicHessianGrid,
                                          const CDensityGrid*        gsDensityGrid,
                                          const CDensityGridQuad*    rwDensityGrid,
                                          const CDensityGrid*        rw2DensityGrid,
                                          const CMemBlock2D<double>& gtoValues,
                                          const CMemBlock2D<double>& gtoValuesX,
                                          const CMemBlock2D<double>& gtoValuesY,
                                          const CMemBlock2D<double>& gtoValuesZ,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const
{

    // set up pointers to exchange-correlation functional derrivatives

    auto df2000 = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhoa);
    
    auto df1100 = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::rhob);
    
    auto df1010 = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::grada);
    
    auto df1001 = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::gradb);
    
    auto df10001 = xcHessianGrid->xcHessianValues(xcvars::rhoa, xcvars::gradab);
        
    auto df0020 = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::grada);
    
    auto df0011 = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::gradb);
    
    auto df00101 = xcHessianGrid->xcHessianValues(xcvars::grada, xcvars::gradab);
    
    auto df00002 = xcHessianGrid->xcHessianValues(xcvars::gradab, xcvars::gradab);
    
    auto df0010 = xcGradientGrid->xcGradientValues(xcvars::grada);
    
    auto df00001 = xcGradientGrid->xcGradientValues(xcvars::gradab);

    auto df00011 = xcHessianGrid->xcHessianValues(xcvars::gradb ,xcvars::gradab);
    
    auto df01001 = xcHessianGrid->xcHessianValues(xcvars::rhob ,xcvars::gradab);

    auto df0110  = xcHessianGrid->xcHessianValues(xcvars::rhob ,xcvars::grada);

    auto df3000 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa );

    auto df2100 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa,xcvars::rhob );

    auto df1200 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhob,xcvars::rhob );

    auto df2010 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa,xcvars::grada );

    auto df0030 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::grada,xcvars::grada );

    auto df0021 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::grada,xcvars::gradb );

    auto df0012 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::gradb,xcvars::gradb );

    auto df00201 =  xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::grada,xcvars::gradab );

    auto df00111 =  xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::gradb,xcvars::gradab );

    auto df00102 =  xcCubicHessianGrid->xcCubicHessianValues(xcvars::grada, xcvars::gradab,xcvars::gradab );

    auto df00003 =  xcCubicHessianGrid->xcCubicHessianValues(xcvars::gradab, xcvars::gradab,xcvars::gradab );

    auto df2001 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa,xcvars::gradb );

    auto df1110 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhob,xcvars::grada );

    auto df1101 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhob,xcvars::gradb );

    auto df20001 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa,xcvars::gradab );

    auto df11001 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::rhob,xcvars::gradab );

    auto df1020 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::grada ,xcvars::grada );

    auto df1011 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::grada ,xcvars::gradb );

    auto df1002 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::gradb ,xcvars::gradb);

    auto df10101 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::grada ,xcvars::gradab);

    auto df10002 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::gradab ,xcvars::gradab);

    auto df01002 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::gradab ,xcvars::gradab);

    auto df0120 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::grada ,xcvars::grada);

    auto df0111 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::grada ,xcvars::gradb);

    auto df01101 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::grada ,xcvars::gradab);

    auto df10011 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhoa, xcvars::gradb ,xcvars::gradab);
    
    auto df01011 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::gradb ,xcvars::gradab);

    auto df0210 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::rhob ,xcvars::grada);

    auto df02001 = xcCubicHessianGrid->xcCubicHessianValues(xcvars::rhob, xcvars::rhob ,xcvars::gradab);

    auto df00021 =  xcCubicHessianGrid->xcCubicHessianValues(xcvars::gradb, xcvars::gradb ,xcvars::gradab);

    // set up pointers to ground state density gradient norms
    
    auto ngrada = gsDensityGrid->alphaDensityGradient(0);
        
    auto grada_x = gsDensityGrid->alphaDensityGradientX(0);
    
    auto grada_y = gsDensityGrid->alphaDensityGradientY(0);
    
    auto grada_z = gsDensityGrid->alphaDensityGradientZ(0);
    
    auto naos = aoKohnShamMatrix->getNumberOfRows();
    
    auto bufdim = _getNumberOfAOsInBuffer();
    
    auto nblock = naos / bufdim;
    
    // set up pointers to exchange-correlation functional derrivatives
    
    
    auto nmat = aoKohnShamMatrix->getNumberOfMatrices();
    
    for (int32_t i = 0; i < nmat; i++)
    {
        // set up pointer to Kohn-Sham matrix
        
        auto ksmat = aoKohnShamMatrix->getMatrix(i);
        
        // set up pointers to perturbed density gradient norms
        
        auto rhow1rhow2 = rwDensityGrid->rhow1rhow2(i);
        
        auto rxw1rhow2 = rwDensityGrid->rxw1rhow2(i);
        
        auto ryw1rhow2 = rwDensityGrid->ryw1rhow2(i);
        
        auto rzw1rhow2 = rwDensityGrid->rzw1rhow2(i);

        auto rxw1rxw2 = rwDensityGrid->rxw1rxw2(i);

        auto rxw1ryw2 = rwDensityGrid->rxw1ryw2(i); 

        auto rxw1rzw2 = rwDensityGrid->rxw1rzw2(i);
        
        auto ryw1rxw2 = rwDensityGrid->ryw1rxw2(i);
        
        auto ryw1ryw2 = rwDensityGrid->ryw1ryw2(i);
        
        auto ryw1rzw2 = rwDensityGrid->ryw1rzw2(i);
        
        auto rzw1rxw2 = rwDensityGrid->rzw1rxw2(i);
        
        auto rzw1ryw2 = rwDensityGrid->rzw1ryw2(i);

        auto rzw1rzw2 = rwDensityGrid->rzw1rzw2(i);

        auto rhow12a = rw2DensityGrid->alphaDensity(i);
        
        auto gradw12a_x = rw2DensityGrid->alphaDensityGradientX(i);
        
        auto gradw12a_y = rw2DensityGrid->alphaDensityGradientY(i);
        
        auto gradw12a_z = rw2DensityGrid->alphaDensityGradientZ(i);
        
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
                        
                        auto omega = bgaos[m] * kgaos[m];
                        
                        auto xomega = bgaox[m] * kgaos[m] + bgaos[m] * kgaox[m];
                        
                        auto yomega = bgaoy[m] * kgaos[m] + bgaos[m] * kgaoy[m];
                        
                        auto zomega = bgaoz[m] * kgaos[m] + bgaos[m] * kgaoz[m];

                        double znva =  1.0 / ngrada[moff + m];

                        double znva3 = 1.0 / std::pow(ngrada[moff + m],3.0);

                        double znva5 = 1.0 / std::pow(ngrada[moff + m],5.0);

                        double rxw12a = gradw12a_x[moff + m];
                        
                        double ryw12a = gradw12a_y[moff + m];
                        
                        double rzw12a = gradw12a_z[moff + m];

                        double xigrad_x = znva * grada_x[moff + m];
                        
                        double xigrad_y = znva * grada_y[moff + m];
                        
                        double xigrad_z = znva * grada_z[moff + m];

                        double xigrad_xx = (znva - grada_x[moff + m] * grada_x[moff + m] * znva3);

                        double xigrad_yy = (znva - grada_y[moff + m] * grada_y[moff + m] * znva3);

                        double xigrad_zz = (znva - grada_z[moff + m] * grada_z[moff + m] * znva3);

                        double xigrad_xy =  - grada_x[moff + m] * grada_y[moff + m] * znva3;

                        double xigrad_xz =  - grada_x[moff + m] * grada_z[moff + m] * znva3;

                        double xigrad_yz =  - grada_y[moff + m] * grada_z[moff + m] * znva3;

                        double xigrad_xxy =   3.0  * grada_x[moff + m] * grada_x[moff + m] * grada_y[moff + m] * znva5 - grada_y[moff + m] * znva3;

                        double xigrad_xxz =   3.0  * grada_x[moff + m] * grada_x[moff + m] * grada_z[moff + m] * znva5 - grada_z[moff + m] * znva3;

                        double xigrad_xyy =   3.0  * grada_x[moff + m] * grada_y[moff + m] * grada_y[moff + m] * znva5 - grada_x[moff + m] * znva3;

                        double xigrad_xzz =   3.0  * grada_x[moff + m] * grada_z[moff + m] * grada_z[moff + m] * znva5 - grada_x[moff + m] * znva3;
                        
                        double xigrad_yzz =   3.0  * grada_y[moff + m] * grada_z[moff + m] * grada_z[moff + m] * znva5 - grada_y[moff + m] * znva3;

                        double xigrad_yyz =   3.0  * grada_y[moff + m] * grada_y[moff + m] * grada_z[moff + m] * znva5 - grada_z[moff + m] * znva3;

                        double xigrad_xyz =   3.0  * grada_x[moff + m] * grada_y[moff + m] * grada_z[moff + m] * znva5;

                        double xigrad_xxx =  3.0 * grada_x[moff + m] * grada_x[moff + m] * grada_x[moff + m] * znva5 - 3.0  * grada_x[moff + m] * znva3;

                        double xigrad_yyy =  3.0 * grada_y[moff + m] * grada_y[moff + m] * grada_y[moff + m] * znva5 - 3.0  * grada_y[moff + m] * znva3;

                        double xigrad_zzz =  3.0 * grada_z[moff + m] * grada_z[moff + m] * grada_z[moff + m] * znva5 - 3.0  * grada_z[moff + m] * znva3;
                        
                        // Various required quantities

                        double xigrad_dot_rw12a = (xigrad_x * rxw12a + xigrad_y * ryw12a + xigrad_z * rzw12a);

                        double xigrad_dot_omega = (xigrad_x * xomega + xigrad_y * yomega + xigrad_z * zomega );

                        double xigrad_dot_rw1rw2 = xigrad_x * rxw1rhow2[moff + m] + xigrad_y * ryw1rhow2[moff + m] + xigrad_z * rzw1rhow2[moff + m];
                        
                        double rw1_dot_rw2 = rxw1rxw2[moff + m] + ryw1ryw2[moff + m] + rzw1rzw2[moff + m]; 

                        double xigrad_dot_rw1rhow2 = xigrad_x * rxw1rhow2[moff + m] + xigrad_y * ryw1rhow2[moff + m] + xigrad_z * rzw1rhow2[moff + m];

                        double grad_dot_rw12 = grada_x[moff + m] * rxw12a + grada_y[moff + m] * ryw12a + grada_z[moff + m] * rzw12a;

                        double grad_dot_omega = grada_x[moff + m] * xomega + grada_y[moff + m] * yomega + grada_z[moff + m] * zomega;

                        double grad_dot_rw1rw2 = grada_x[moff + m] * rxw1rhow2[moff + m] + grada_y[moff + m] * ryw1rhow2[moff + m] + grada_z[moff + m] * rzw1rhow2[moff + m];

                        double omega_dot_rw1rhow2 = xomega * rxw1rhow2[moff + m] + yomega * ryw1rhow2[moff + m] + zomega  * rzw1rhow2[moff + m];

                        double grad_dot_rw1rhow2 = grada_x[moff + m] * rxw1rhow2[moff + m] + grada_y[moff + m] * ryw1rhow2[moff + m] + grada_z[moff + m] * rzw1rhow2[moff + m];

                        double xigrad_dot_rw1_xigrad_dot_rw2  = xigrad_x * xigrad_x * rxw1rxw2[moff + m]
                                                              + xigrad_x * xigrad_y * rxw1ryw2[moff + m]
                                                              + xigrad_x * xigrad_z * rxw1rzw2[moff + m]
                                                              + xigrad_y * xigrad_x * ryw1rxw2[moff + m]
                                                              + xigrad_y * xigrad_y * ryw1ryw2[moff + m]
                                                              + xigrad_y * xigrad_z * ryw1rzw2[moff + m]
                                                              + xigrad_z * xigrad_x * rzw1rxw2[moff + m]
                                                              + xigrad_z * xigrad_y * rzw1ryw2[moff + m]
                                                              + xigrad_z * xigrad_z * rzw1rzw2[moff + m];

                        double first = w * (df2000[moff + m] * rhow12a[moff + m] + df1100[moff + m] * rhow12a[moff + m]) * omega;

                        double second = w * (df1010[moff + m] + df1001[moff + m])  * xigrad_dot_rw12a * omega;

                        second += w * df10001[moff + m] * 2.0 * grad_dot_rw12 * omega;  
                        
                        double third = w * (df1010[moff + m] + df0110[moff + m] ) * xigrad_dot_omega * rhow12a[moff + m];

                        third += w * (df10001[moff + m ] + df01001[moff + m]) * grad_dot_omega * rhow12a[moff + m];

                        double fourth = w * df0010[moff + m] * ( (xigrad_xx  * xomega + xigrad_xy  * yomega +  xigrad_xz  * zomega )* rxw12a

                                                     +  (xigrad_xy  * xomega + xigrad_yy  * yomega +  xigrad_yz  * zomega ) * ryw12a

                                                     +  (xigrad_xz  * xomega + xigrad_yz  * yomega +  xigrad_zz  * zomega ) * rzw12a 
                                                                                                                                
                                                                                                                                );

                        fourth += w * df00001[moff + m] * ( xomega * rxw12a + yomega * ryw12a + zomega * rzw12a );

                        fourth += w * df00101[moff + m ] * xigrad_dot_omega * (2.0 * grad_dot_rw12 ) ;

                        fourth += w * (df00101[moff + m ] * xigrad_dot_rw12a   + df00011[moff + m ] * xigrad_dot_rw12a) * grad_dot_omega; 

                        fourth += w * df00002[moff + m ] * grad_dot_omega * (2.0 * grad_dot_rw12);

                        fourth += w * (df0020[moff + m] + df0011[moff + m])* xigrad_dot_omega * xigrad_dot_rw12a;
                        
                        double fifth = w * (df3000[moff + m] + 2.0 * df2100[moff + m] + df1200[moff + m]) * rhow1rhow2[moff + m] * omega;

                        double seventh = w * (df2010[moff + m] + df2001[moff + m] ) * xigrad_dot_rw1rw2 * omega; 

                        seventh += w * (df1110[moff + m] + df1101[moff + m] ) * xigrad_dot_rw1rw2 * omega;     

                        seventh += w * 2.0 * (df20001[moff + m] + df11001[moff + m] ) * grad_dot_rw1rw2 * omega;

                        double eighth  = w * (df1020[moff + m] + 2.0 * df1011[moff + m] + df1002[moff + m] ) * xigrad_dot_rw1_xigrad_dot_rw2 * omega;

                        eighth  += w * (df1010[moff + m] + df1001[moff + m]  ) * (xigrad_xx *  rxw1rxw2[moff + m] +   xigrad_xy *  rxw1ryw2[moff + m] + xigrad_xz *  rxw1rzw2[moff + m] 

                                                   + xigrad_xy *  ryw1rxw2[moff + m] + xigrad_yy *  ryw1ryw2[moff + m] + xigrad_yz *  ryw1rzw2[moff + m]

                                                  + xigrad_xz *  rzw1rxw2[moff + m] + xigrad_yz *  rzw1ryw2[moff + m]  + xigrad_zz *  rzw1rzw2[moff + m]) * omega;


                        eighth += w * 2.0 * (df10101[moff + m] + df10101[moff + m]) * ngrada[moff + m] * xigrad_dot_rw1_xigrad_dot_rw2 *omega;
                        

                        eighth += w * 4.0 * df10002[moff + m]  * ngrada[moff + m] * ngrada[moff + m]*  xigrad_dot_rw1_xigrad_dot_rw2 *omega;
                        
                        eighth += w * 2.0 * df10001[moff + m] * rw1_dot_rw2 * omega;

                        double ninth = w * (df2010[moff + m] + 2.0 * df1110[moff + m] + df0210[moff + m] ) * rhow1rhow2[moff + m] * xigrad_dot_omega;

                        ninth += w * (df20001[moff + m] + 2.0 * df11001[moff + m] + df02001[moff + m]) *  grad_dot_omega * rhow1rhow2[moff + m];

                        double tenth = w * (df1010[moff + m] + df0110[moff + m] ) * 
                        
                                                ( (xigrad_xx  * rxw1rhow2[moff + m] + xigrad_xy  * ryw1rhow2[moff + m] +  xigrad_xz  * rzw1rhow2[moff + m] )* xomega

                                              +  (xigrad_xy  * rxw1rhow2[moff + m] + xigrad_yy  * ryw1rhow2[moff + m] +  xigrad_yz  * rzw1rhow2[moff + m] ) * yomega

                                              +  (xigrad_xz  * rxw1rhow2[moff + m] + xigrad_yz  * ryw1rhow2[moff + m] +  xigrad_zz  * rzw1rhow2[moff + m] ) * zomega 
                                                                                                                                
                                                                                                                                );


                        tenth += w * (df1020[moff + m ] + df1011[moff + m ] + df0120[moff + m ] + df0111[moff + m ] ) * xigrad_dot_rw1rhow2 * xigrad_dot_omega;

                        tenth += w * (df10001[moff + m] + df01001[moff + m]) * omega_dot_rw1rhow2;

                        tenth += w * (df10101[moff + m] + df10011[moff + m]+  df01101[moff + m] + df0111[moff + m] + df01011[moff + m] )* xigrad_dot_rw1rhow2 * grad_dot_omega;

                        tenth += w * (df10101[moff + m] + df10011[moff + m]+  df01101[moff + m] + df0111[moff + m] ) * grad_dot_rw1rhow2 * xigrad_dot_omega;

                        tenth += w * (df10002[moff + m] + df01002[moff + m] ) * grad_dot_rw1rhow2 * grad_dot_omega;
 
                        double twelthfirst =  xigrad_xxx * xomega * rxw1rxw2[moff + m]
                                            + xigrad_xxy * xomega * rxw1ryw2[moff + m]
                                            + xigrad_xxz * xomega * rxw1rzw2[moff + m]
                                            + xigrad_xxy * xomega * ryw1rxw2[moff + m]
                                            + xigrad_xyy * xomega * ryw1ryw2[moff + m]
                                            + xigrad_xyz * xomega * ryw1rzw2[moff + m]
                                            + xigrad_xxz * xomega * rzw1rxw2[moff + m]
                                            + xigrad_xyz * xomega * rzw1ryw2[moff + m]
                                            + xigrad_xzz * xomega * rzw1rzw2[moff + m]
                                            + xigrad_xxy * yomega * rxw1rxw2[moff + m]
                                            + xigrad_xyy * yomega * rxw1ryw2[moff + m]
                                            + xigrad_xyz * yomega * rxw1rzw2[moff + m]
                                            + xigrad_xyy * yomega * ryw1rxw2[moff + m]
                                            + xigrad_yyy * yomega * ryw1ryw2[moff + m]
                                            + xigrad_yyz * yomega * ryw1rzw2[moff + m]
                                            + xigrad_xyz * yomega * rzw1rxw2[moff + m]
                                            + xigrad_yyz * yomega * rzw1ryw2[moff + m]
                                            + xigrad_yzz * yomega * rzw1rzw2[moff + m]
                                            + xigrad_xxz * zomega * rxw1rxw2[moff + m]
                                            + xigrad_xyz * zomega * rxw1ryw2[moff + m]
                                            + xigrad_xzz * zomega * rxw1rzw2[moff + m]
                                            + xigrad_xyz * zomega * ryw1rxw2[moff + m]
                                            + xigrad_yyz * zomega * ryw1ryw2[moff + m]
                                            + xigrad_yzz * zomega * ryw1rzw2[moff + m]
                                            + xigrad_xzz * zomega * rzw1rxw2[moff + m]
                                            + xigrad_yzz * zomega * rzw1ryw2[moff + m]
                                            + xigrad_zzz * zomega * rzw1rzw2[moff + m];


                        double twelthsecond = xigrad_xx * xigrad_x * xomega * rxw1rxw2[moff + m]
                                            + xigrad_xx * xigrad_y * yomega * rxw1rxw2[moff + m]
                                            + xigrad_xx * xigrad_z * zomega * rxw1rxw2[moff + m]
                                            + xigrad_xy * xigrad_x * xomega * rxw1ryw2[moff + m]
                                            + xigrad_xy * xigrad_y * yomega * rxw1ryw2[moff + m]
                                            + xigrad_xy * xigrad_z * zomega * rxw1ryw2[moff + m]
                                            + xigrad_xz * xigrad_x * xomega * rxw1rzw2[moff + m]
                                            + xigrad_xz * xigrad_y * yomega * rxw1rzw2[moff + m]
                                            + xigrad_xz * xigrad_z * zomega * rxw1rzw2[moff + m]
                                            + xigrad_xy * xigrad_x * xomega * ryw1rxw2[moff + m]
                                            + xigrad_xy * xigrad_y * yomega * ryw1rxw2[moff + m]
                                            + xigrad_xy * xigrad_z * zomega * ryw1rxw2[moff + m]
                                            + xigrad_yy * xigrad_x * xomega * ryw1ryw2[moff + m]
                                            + xigrad_yy * xigrad_y * yomega * ryw1ryw2[moff + m]
                                            + xigrad_yy * xigrad_z * zomega * ryw1ryw2[moff + m]
                                            + xigrad_yz * xigrad_x * xomega * ryw1rzw2[moff + m]
                                            + xigrad_yz * xigrad_y * yomega * ryw1rzw2[moff + m]
                                            + xigrad_yz * xigrad_z * zomega * ryw1rzw2[moff + m]
                                            + xigrad_xz * xigrad_x * xomega * rzw1rxw2[moff + m]
                                            + xigrad_xz * xigrad_y * yomega * rzw1rxw2[moff + m]
                                            + xigrad_xz * xigrad_z * zomega * rzw1rxw2[moff + m]
                                            + xigrad_yz * xigrad_x * xomega * rzw1ryw2[moff + m]
                                            + xigrad_yz * xigrad_y * yomega * rzw1ryw2[moff + m]
                                            + xigrad_yz * xigrad_z * zomega * rzw1ryw2[moff + m]
                                            + xigrad_zz * xigrad_x * xomega * rzw1rzw2[moff + m]
                                            + xigrad_zz * xigrad_y * yomega * rzw1rzw2[moff + m]
                                            + xigrad_zz * xigrad_z * zomega * rzw1rzw2[moff + m];

                        double twelththird   =    xigrad_xx * xigrad_x * xomega * (rxw1rxw2[moff + m] + rxw1rxw2[moff + m])
                                                + xigrad_xx * xigrad_y * xomega * (ryw1rxw2[moff + m] + rxw1ryw2[moff + m])
                                                + xigrad_xx * xigrad_z * xomega * (rzw1rxw2[moff + m] + rxw1rzw2[moff + m])
                                                + xigrad_xy * xigrad_x * xomega * (rxw1ryw2[moff + m] + ryw1rxw2[moff + m])
                                                + xigrad_xy * xigrad_y * xomega * (ryw1ryw2[moff + m] + ryw1ryw2[moff + m])
                                                + xigrad_xy * xigrad_z * xomega * (rzw1ryw2[moff + m] + ryw1rzw2[moff + m])
                                                + xigrad_xz * xigrad_x * xomega * (rxw1rzw2[moff + m] + rzw1rxw2[moff + m])
                                                + xigrad_xz * xigrad_y * xomega * (ryw1rzw2[moff + m] + rzw1ryw2[moff + m])
                                                + xigrad_xz * xigrad_z * xomega * (rzw1rzw2[moff + m] + rzw1rzw2[moff + m])
                                                + xigrad_xy * xigrad_x * yomega * (rxw1rxw2[moff + m] + rxw1rxw2[moff + m])
                                                + xigrad_xy * xigrad_y * yomega * (ryw1rxw2[moff + m] + rxw1ryw2[moff + m])
                                                + xigrad_xy * xigrad_z * yomega * (rzw1rxw2[moff + m] + rxw1rzw2[moff + m])
                                                + xigrad_yy * xigrad_x * yomega * (rxw1ryw2[moff + m] + ryw1rxw2[moff + m])
                                                + xigrad_yy * xigrad_y * yomega * (ryw1ryw2[moff + m] + ryw1ryw2[moff + m])
                                                + xigrad_yy * xigrad_z * yomega * (rzw1ryw2[moff + m] + ryw1rzw2[moff + m])
                                                + xigrad_yz * xigrad_x * yomega * (rxw1rzw2[moff + m] + rzw1rxw2[moff + m])
                                                + xigrad_yz * xigrad_y * yomega * (ryw1rzw2[moff + m] + rzw1ryw2[moff + m])
                                                + xigrad_yz * xigrad_z * yomega * (rzw1rzw2[moff + m] + rzw1rzw2[moff + m])
                                                + xigrad_xz * xigrad_x * zomega * (rxw1rxw2[moff + m] + rxw1rxw2[moff + m])
                                                + xigrad_xz * xigrad_y * zomega * (ryw1rxw2[moff + m] + rxw1ryw2[moff + m])
                                                + xigrad_xz * xigrad_z * zomega * (rzw1rxw2[moff + m] + rxw1rzw2[moff + m])
                                                + xigrad_yz * xigrad_x * zomega * (rxw1ryw2[moff + m] + ryw1rxw2[moff + m])
                                                + xigrad_yz * xigrad_y * zomega * (ryw1ryw2[moff + m] + ryw1ryw2[moff + m])
                                                + xigrad_yz * xigrad_z * zomega * (rzw1ryw2[moff + m] + ryw1rzw2[moff + m])
                                                + xigrad_zz * xigrad_x * zomega * (rxw1rzw2[moff + m] + rzw1rxw2[moff + m])
                                                + xigrad_zz * xigrad_y * zomega * (ryw1rzw2[moff + m] + rzw1ryw2[moff + m])
                                                + xigrad_zz * xigrad_z * zomega * (rzw1rzw2[moff + m] + rzw1rzw2[moff + m]);

                        double twelthfourth_gam =       xigrad_x * xomega * rxw1rxw2[moff + m]
                                                    +   xigrad_x * yomega * rxw1ryw2[moff + m]
                                                    +   xigrad_x * zomega * rxw1rzw2[moff + m]
                                                    +   xigrad_y * xomega * ryw1rxw2[moff + m]
                                                    +   xigrad_y * yomega * ryw1ryw2[moff + m]
                                                    +   xigrad_y * zomega * ryw1rzw2[moff + m]
                                                    +   xigrad_z * xomega * rzw1rxw2[moff + m]
                                                    +   xigrad_z * yomega * rzw1ryw2[moff + m]
                                                    +   xigrad_z * zomega * rzw1rzw2[moff + m];

                        double twelthfifth_gam =    (xigrad_x * grada_x[moff + m] + grada_x[moff +m] * xigrad_x) * rxw1rxw2[moff +m]
                                                +   (xigrad_x * grada_y[moff + m] + grada_x[moff +m] * xigrad_y) * rxw1ryw2[moff +m]
                                                +   (xigrad_x * grada_z[moff + m] + grada_x[moff +m] * xigrad_z) * rxw1rzw2[moff +m]
                                                +   (xigrad_y * grada_x[moff + m] + grada_y[moff +m] * xigrad_x) * rxw1rxw2[moff +m]
                                                +   (xigrad_y * grada_y[moff + m] + grada_y[moff +m] * xigrad_y) * rxw1ryw2[moff +m]
                                                +   (xigrad_y * grada_z[moff + m] + grada_y[moff +m] * xigrad_z) * rxw1rzw2[moff +m]
                                                +   (xigrad_z * grada_x[moff + m] + grada_z[moff +m] * xigrad_x) * rxw1rxw2[moff +m]
                                                +   (xigrad_z * grada_y[moff + m] + grada_z[moff +m] * xigrad_y) * rxw1ryw2[moff +m]
                                                +   (xigrad_z * grada_z[moff + m] + grada_z[moff +m] * xigrad_z) * rxw1rzw2[moff +m];


                        double twelfth = w * df0010[moff + m] * twelthfirst;

                        twelfth += w * (df0020[moff + m] +  df0011[moff + m]) * (twelththird + twelthsecond);

                        twelfth += w * (df0030[moff + m] + 2.0 * df0021[moff + m] + df0012[moff + m]) * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;

                        // second-term mixed density gradient

                        twelfth += w * (df00101[moff + m] + df00011[moff +m]) * ngrada[moff + m] * twelthsecond;

                        twelfth += w * (df00101[moff + m] + df00011[moff +m]) * xigrad_dot_omega * rw1_dot_rw2;

                        twelfth += w * df00002[moff + m] * grad_dot_omega * rw1_dot_rw2;

                        // third-term mixed denisty gradient


                        twelfth += w *  (df00101[moff + m] + df00011[moff +m]) * twelthfourth_gam;


                        twelfth += w * df00101[moff + m] * ngrada[moff + m] * twelththird;

                        twelfth += w * df00002[moff +m] * ngrada[moff + m] * twelthfourth_gam;
                        
                        // Fourth terms 

                        twelfth += w * (df00201[moff + m] + df00111[moff + m]) * twelthfifth_gam * xigrad_dot_omega;

                        twelfth += w * df00102[moff + m] * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;
                        
                        twelfth += w * (df00201[moff + m] + 2 * df00111[moff + m] + df00021[moff + m]) * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;

                        twelfth += w * (df00102[moff + m] + df00011[moff + m]) * ngrada[moff + m] * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;

                        twelfth += w * df00003[moff + m] * ngrada[moff + m] * ngrada[moff +m] * grad_dot_omega;
                                     
                        fvxc += first;
                        fvxc += second;
                        fvxc += third;
                        fvxc += fourth;
                        fvxc += fifth;
                        fvxc += seventh;
                        fvxc += eighth ;
                        fvxc += ninth;
                        fvxc += tenth;
                        fvxc += twelfth;
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
                        
                        auto omega = bgaos[l] * kgaos[l];
                        
                        auto xomega = bgaox[l] * kgaos[l] + bgaos[l] * kgaox[l];
                        
                        auto yomega = bgaoy[l] * kgaos[l] + bgaos[l] * kgaoy[l];
                        
                        auto zomega = bgaoz[l] * kgaos[l] + bgaos[l] * kgaoz[l];

                        double znva =  1.0 / ngrada[loff + l];

                        double znva3 = 1.0 / std::pow(ngrada[loff + l],3.0);

                        double znva5 = 1.0 / std::pow(ngrada[loff + l],5.0);
                                                
                        double xigrad_x = znva * grada_x[loff + l];
                        
                        double xigrad_y = znva * grada_y[loff + l];
                        
                        double xigrad_z = znva * grada_z[loff + l];

                        double rxw12a = gradw12a_x[loff + l];
                        
                        double ryw12a = gradw12a_y[loff + l];
                        
                        double rzw12a = gradw12a_z[loff + l];

                        double xigrad_xx = (znva - grada_x[loff + l] * grada_x[loff + l] * znva3);

                        double xigrad_yy = (znva - grada_y[loff + l] * grada_y[loff + l] * znva3);

                        double xigrad_zz = (znva - grada_z[loff + l] * grada_z[loff + l] * znva3);

                        double xigrad_xy =  - grada_x[loff + l] * grada_y[loff + l] * znva3;

                        double xigrad_xz =  - grada_x[loff + l] * grada_z[loff + l] * znva3;

                        double xigrad_yz =  - grada_y[loff + l] * grada_z[loff + l] * znva3;
                        
                        double xigrad_dot_rw12a = (xigrad_x * rxw12a + xigrad_y * ryw12a + xigrad_z * rzw12a);

                        double xigrad_dot_omega = (xigrad_x * xomega + xigrad_y * yomega + xigrad_z * zomega );

                        double xigrad_dot_rw1rw2 = xigrad_x * rxw1rhow2[loff + l] + xigrad_y * ryw1rhow2[loff + l] + xigrad_z * rzw1rhow2[loff + l];
                        
                        double rw1_dot_rw2 = rxw1rxw2[loff + l] + ryw1ryw2[loff + l] + rzw1rzw2[loff + l]; 

                        double xigrad_dot_rw1rhow2 = xigrad_x * rxw1rhow2[loff + l] + xigrad_y * ryw1rhow2[loff + l] + xigrad_z * rzw1rhow2[loff + l];

                        double grad_dot_rw12 = grada_x[loff + l] * rxw12a + grada_y[loff + l] * ryw12a + grada_z[loff + l] * rzw12a;

                        double grad_dot_omega = grada_x[loff + l] * xomega + grada_y[loff + l] * yomega + grada_z[loff + l] * zomega;

                        double grad_dot_rw1rw2 = grada_x[loff + l] * rxw1rhow2[loff + l] + grada_y[loff + l] * ryw1rhow2[loff + l] + grada_z[loff + l] * rzw1rhow2[loff + l];

                        double omega_dot_rw1rhow2 = xomega * rxw1rhow2[loff + l] + yomega * ryw1rhow2[loff + l] + zomega  * rzw1rhow2[loff + l];

                        double grad_dot_rw1rhow2 = grada_x[loff + l] * rxw1rhow2[loff + l] + grada_y[loff + l] * ryw1rhow2[loff + l] + grada_z[loff + l] * rzw1rhow2[loff + l];

                        double xigrad_dot_rw1_xigrad_dot_rw2  = xigrad_x * xigrad_x * rxw1rxw2[loff + l]
                                                              + xigrad_x * xigrad_y * rxw1ryw2[loff + l]
                                                              + xigrad_x * xigrad_z * rxw1rzw2[loff + l]
                                                              + xigrad_y * xigrad_x * ryw1rxw2[loff + l]
                                                              + xigrad_y * xigrad_y * ryw1ryw2[loff + l]
                                                              + xigrad_y * xigrad_z * ryw1rzw2[loff + l]
                                                              + xigrad_z * xigrad_x * rzw1rxw2[loff + l]
                                                              + xigrad_z * xigrad_y * rzw1ryw2[loff + l]
                                                              + xigrad_z * xigrad_z * rzw1rzw2[loff + l];

                        double first = w * (df2000[loff + l] * rhow12a[loff + l] + df1100[loff + l] * rhow12a[loff + l]) * omega;

                        double second = w * (df1010[loff + l] + df1001[loff + l])  * xigrad_dot_rw12a * omega;

                        second += w * df10001[loff + l] * 2.0 * grad_dot_rw12 * omega;  
                        
                        double third = w * (df1010[loff + l] + df0110[loff + l] ) * xigrad_dot_omega * rhow12a[loff + l];

                        third += w * (df10001[loff + l ] + df01001[loff + l]) * grad_dot_omega * rhow12a[loff + l];

                        double fourth = w * df0010[loff + l] * ( (xigrad_xx  * xomega + xigrad_xy  * yomega +  xigrad_xz  * zomega )* rxw12a

                                                     +  (xigrad_xy  * xomega + xigrad_yy  * yomega +  xigrad_yz  * zomega ) * ryw12a

                                                     +  (xigrad_xz  * xomega + xigrad_yz  * yomega +  xigrad_zz  * zomega ) * rzw12a 
                                                                                                                                
                                                                                                                                );

                        fourth += w * df00001[loff + l] * ( xomega * rxw12a + yomega * ryw12a + zomega * rzw12a );

                        fourth += w * df00101[loff + l ] * xigrad_dot_omega * (2.0 * grad_dot_rw12 ) ;

                        fourth += w * (df00101[loff + l ] * xigrad_dot_rw12a   + df00011[loff + l ] * xigrad_dot_rw12a) * grad_dot_omega; 

                        fourth += w * df00002[loff + l ] * grad_dot_omega * (2.0 * grad_dot_rw12);

                        fourth += w * (df0020[loff + l] + df0011[loff + l])* xigrad_dot_omega * xigrad_dot_rw12a;
                        
                        double fifth = w * (df3000[loff + l] + 2.0 * df2100[loff + l] + df1200[loff + l]) * rhow1rhow2[loff + l] * omega;

                        double seventh = w * (df2010[loff + l] + df2001[loff + l] ) * xigrad_dot_rw1rw2 * omega; 

                        seventh += w * (df1110[loff + l] + df1101[loff + l] ) * xigrad_dot_rw1rw2 * omega;     

                        seventh += w * 2.0 * (df20001[loff + l] + df11001[loff + l] ) * grad_dot_rw1rw2 * omega;

                        double eighth  = w * (df1020[loff + l] + 2.0 * df1011[loff + l] + df1002[loff + l] ) * xigrad_dot_rw1_xigrad_dot_rw2 * omega;

                        eighth  += w * (df1010[loff + l] + df1001[loff + l]  ) * (xigrad_xx *  rxw1rxw2[loff + l] +   xigrad_xy *  rxw1ryw2[loff + l] + xigrad_xz *  rxw1rzw2[loff + l] 

                                                   + xigrad_xy *  ryw1rxw2[loff + l] + xigrad_yy *  ryw1ryw2[loff + l] + xigrad_yz *  ryw1rzw2[loff + l]

                                                  + xigrad_xz *  rzw1rxw2[loff + l] + xigrad_yz *  rzw1ryw2[loff + l]  + xigrad_zz *  rzw1rzw2[loff + l]) * omega;


                        eighth += w * 2.0 * (df10101[loff + l] + df10101[loff + l]) * ngrada[loff + l] * xigrad_dot_rw1_xigrad_dot_rw2 *omega;
                        

                        eighth += w * 4.0 * df10002[loff + l]  * ngrada[loff + l] * ngrada[loff + l]*  xigrad_dot_rw1_xigrad_dot_rw2 *omega;
                        
                        eighth += w * 2.0 * df10001[loff + l] * rw1_dot_rw2 * omega;

                        double ninth = w * (df2010[loff + l] + 2.0 * df1110[loff + l] + df0210[loff + l] ) * rhow1rhow2[loff + l] * xigrad_dot_omega;

                        ninth += w * (df20001[loff + l] + 2.0 * df11001[loff + l] + df02001[loff + l]) *  grad_dot_omega * rhow1rhow2[loff + l];

                        double tenth = w * (df1010[loff + l] + df0110[loff + l] ) * 
                        
                                                ( (xigrad_xx  * rxw1rhow2[loff + l] + xigrad_xy  * ryw1rhow2[loff + l] +  xigrad_xz  * rzw1rhow2[loff + l] )* xomega

                                              +  (xigrad_xy  * rxw1rhow2[loff + l] + xigrad_yy  * ryw1rhow2[loff + l] +  xigrad_yz  * rzw1rhow2[loff + l] ) * yomega

                                              +  (xigrad_xz  * rxw1rhow2[loff + l] + xigrad_yz  * ryw1rhow2[loff + l] +  xigrad_zz  * rzw1rhow2[loff + l] ) * zomega 
                                                                                                                                
                                                                                                                                );

                        tenth += w * (df1020[loff + l ] + df1011[loff + l ] + df0120[loff + l ] + df0111[loff + l ] ) * xigrad_dot_rw1rhow2 * xigrad_dot_omega;

                        tenth += w * (df10001[loff + l] + df01001[loff + l]) * omega_dot_rw1rhow2;

                        tenth += w * (df10101[loff + l] + df10011[loff + l]+  df01101[loff + l] + df0111[loff + l] +  df01011[loff + l] )* xigrad_dot_rw1rhow2 * grad_dot_omega;

                        tenth += w * (df10101[loff + l] + df10011[loff + l]+  df01101[loff + l] + df0111[loff + l] ) * grad_dot_rw1rhow2 * xigrad_dot_omega;

                        tenth += w * (df10002[loff + l] + df01002[loff + l] ) * grad_dot_rw1rhow2 * grad_dot_omega;

                        double xigrad_xxy =   3.0  * grada_x[loff + l] * grada_x[loff + l] * grada_y[loff + l] * znva5 - grada_y[loff + l] * znva3;

                        double xigrad_xxz =   3.0  * grada_x[loff + l] * grada_x[loff + l] * grada_z[loff + l] * znva5 - grada_z[loff + l] * znva3;

                        double xigrad_xyy =   3.0  * grada_x[loff + l] * grada_y[loff + l] * grada_y[loff + l] * znva5 - grada_x[loff + l] * znva3;

                        double xigrad_xzz =   3.0  * grada_x[loff + l] * grada_z[loff + l] * grada_z[loff + l] * znva5 - grada_x[loff + l] * znva3;
                        
                        double xigrad_yzz =   3.0  * grada_y[loff + l] * grada_z[loff + l] * grada_z[loff + l] * znva5 - grada_y[loff + l] * znva3;

                        double xigrad_yyz =   3.0  * grada_y[loff + l] * grada_y[loff + l] * grada_z[loff + l] * znva5 - grada_z[loff + l] * znva3;

                        double xigrad_xyz =   3.0  * grada_x[loff + l] * grada_y[loff + l] * grada_z[loff + l] * znva5;

                        double xigrad_xxx =  3.0 * grada_x[loff + l] * grada_x[loff + l] * grada_x[loff + l] * znva5 - 3.0  * grada_x[loff + l] * znva3;

                        double xigrad_yyy =  3.0 * grada_y[loff + l] * grada_y[loff + l] * grada_y[loff + l] * znva5 - 3.0  * grada_y[loff + l] * znva3;

                        double xigrad_zzz =  3.0 * grada_z[loff + l] * grada_z[loff + l] * grada_z[loff + l] * znva5 - 3.0  * grada_z[loff + l] * znva3;

                        double twelthfirst =  xigrad_xxx * xomega * rxw1rxw2[loff + l]
                                            + xigrad_xxy * xomega * rxw1ryw2[loff + l]
                                            + xigrad_xxz * xomega * rxw1rzw2[loff + l]
                                            + xigrad_xxy * xomega * ryw1rxw2[loff + l]
                                            + xigrad_xyy * xomega * ryw1ryw2[loff + l]
                                            + xigrad_xyz * xomega * ryw1rzw2[loff + l]
                                            + xigrad_xxz * xomega * rzw1rxw2[loff + l]
                                            + xigrad_xyz * xomega * rzw1ryw2[loff + l]
                                            + xigrad_xzz * xomega * rzw1rzw2[loff + l]
                                            + xigrad_xxy * yomega * rxw1rxw2[loff + l]
                                            + xigrad_xyy * yomega * rxw1ryw2[loff + l]
                                            + xigrad_xyz * yomega * rxw1rzw2[loff + l]
                                            + xigrad_xyy * yomega * ryw1rxw2[loff + l]
                                            + xigrad_yyy * yomega * ryw1ryw2[loff + l]
                                            + xigrad_yyz * yomega * ryw1rzw2[loff + l]
                                            + xigrad_xyz * yomega * rzw1rxw2[loff + l]
                                            + xigrad_yyz * yomega * rzw1ryw2[loff + l]
                                            + xigrad_yzz * yomega * rzw1rzw2[loff + l]
                                            + xigrad_xxz * zomega * rxw1rxw2[loff + l]
                                            + xigrad_xyz * zomega * rxw1ryw2[loff + l]
                                            + xigrad_xzz * zomega * rxw1rzw2[loff + l]
                                            + xigrad_xyz * zomega * ryw1rxw2[loff + l]
                                            + xigrad_yyz * zomega * ryw1ryw2[loff + l]
                                            + xigrad_yzz * zomega * ryw1rzw2[loff + l]
                                            + xigrad_xzz * zomega * rzw1rxw2[loff + l]
                                            + xigrad_yzz * zomega * rzw1ryw2[loff + l]
                                            + xigrad_zzz * zomega * rzw1rzw2[loff + l];

                        double twelthsecond  = xigrad_xx * xigrad_x * xomega * rxw1rxw2[loff + l]
                                            + xigrad_xx * xigrad_y * yomega * rxw1rxw2[loff + l]
                                            + xigrad_xx * xigrad_z * zomega * rxw1rxw2[loff + l]
                                            + xigrad_xy * xigrad_x * xomega * rxw1ryw2[loff + l]
                                            + xigrad_xy * xigrad_y * yomega * rxw1ryw2[loff + l]
                                            + xigrad_xy * xigrad_z * zomega * rxw1ryw2[loff + l]
                                            + xigrad_xz * xigrad_x * xomega * rxw1rzw2[loff + l]
                                            + xigrad_xz * xigrad_y * yomega * rxw1rzw2[loff + l]
                                            + xigrad_xz * xigrad_z * zomega * rxw1rzw2[loff + l]
                                            + xigrad_xy * xigrad_x * xomega * ryw1rxw2[loff + l]
                                            + xigrad_xy * xigrad_y * yomega * ryw1rxw2[loff + l]
                                            + xigrad_xy * xigrad_z * zomega * ryw1rxw2[loff + l]
                                            + xigrad_yy * xigrad_x * xomega * ryw1ryw2[loff + l]
                                            + xigrad_yy * xigrad_y * yomega * ryw1ryw2[loff + l]
                                            + xigrad_yy * xigrad_z * zomega * ryw1ryw2[loff + l]
                                            + xigrad_yz * xigrad_x * xomega * ryw1rzw2[loff + l]
                                            + xigrad_yz * xigrad_y * yomega * ryw1rzw2[loff + l]
                                            + xigrad_yz * xigrad_z * zomega * ryw1rzw2[loff + l]
                                            + xigrad_xz * xigrad_x * xomega * rzw1rxw2[loff + l]
                                            + xigrad_xz * xigrad_y * yomega * rzw1rxw2[loff + l]
                                            + xigrad_xz * xigrad_z * zomega * rzw1rxw2[loff + l]
                                            + xigrad_yz * xigrad_x * xomega * rzw1ryw2[loff + l]
                                            + xigrad_yz * xigrad_y * yomega * rzw1ryw2[loff + l]
                                            + xigrad_yz * xigrad_z * zomega * rzw1ryw2[loff + l]
                                            + xigrad_zz * xigrad_x * xomega * rzw1rzw2[loff + l]
                                            + xigrad_zz * xigrad_y * yomega * rzw1rzw2[loff + l]
                                            + xigrad_zz * xigrad_z * zomega * rzw1rzw2[loff + l];

                        double twelththird   =    xigrad_xx * xigrad_x * xomega * (rxw1rxw2[loff + l] + rxw1rxw2[loff + l])
                                                + xigrad_xx * xigrad_y * xomega * (ryw1rxw2[loff + l] + rxw1ryw2[loff + l])
                                                + xigrad_xx * xigrad_z * xomega * (rzw1rxw2[loff + l] + rxw1rzw2[loff + l])
                                                + xigrad_xy * xigrad_x * xomega * (rxw1ryw2[loff + l] + ryw1rxw2[loff + l])
                                                + xigrad_xy * xigrad_y * xomega * (ryw1ryw2[loff + l] + ryw1ryw2[loff + l])
                                                + xigrad_xy * xigrad_z * xomega * (rzw1ryw2[loff + l] + ryw1rzw2[loff + l])
                                                + xigrad_xz * xigrad_x * xomega * (rxw1rzw2[loff + l] + rzw1rxw2[loff + l])
                                                + xigrad_xz * xigrad_y * xomega * (ryw1rzw2[loff + l] + rzw1ryw2[loff + l])
                                                + xigrad_xz * xigrad_z * xomega * (rzw1rzw2[loff + l] + rzw1rzw2[loff + l])
                                                + xigrad_xy * xigrad_x * yomega * (rxw1rxw2[loff + l] + rxw1rxw2[loff + l])
                                                + xigrad_xy * xigrad_y * yomega * (ryw1rxw2[loff + l] + rxw1ryw2[loff + l])
                                                + xigrad_xy * xigrad_z * yomega * (rzw1rxw2[loff + l] + rxw1rzw2[loff + l])
                                                + xigrad_yy * xigrad_x * yomega * (rxw1ryw2[loff + l] + ryw1rxw2[loff + l])
                                                + xigrad_yy * xigrad_y * yomega * (ryw1ryw2[loff + l] + ryw1ryw2[loff + l])
                                                + xigrad_yy * xigrad_z * yomega * (rzw1ryw2[loff + l] + ryw1rzw2[loff + l])
                                                + xigrad_yz * xigrad_x * yomega * (rxw1rzw2[loff + l] + rzw1rxw2[loff + l])
                                                + xigrad_yz * xigrad_y * yomega * (ryw1rzw2[loff + l] + rzw1ryw2[loff + l])
                                                + xigrad_yz * xigrad_z * yomega * (rzw1rzw2[loff + l] + rzw1rzw2[loff + l])
                                                + xigrad_xz * xigrad_x * zomega * (rxw1rxw2[loff + l] + rxw1rxw2[loff + l])
                                                + xigrad_xz * xigrad_y * zomega * (ryw1rxw2[loff + l] + rxw1ryw2[loff + l])
                                                + xigrad_xz * xigrad_z * zomega * (rzw1rxw2[loff + l] + rxw1rzw2[loff + l])
                                                + xigrad_yz * xigrad_x * zomega * (rxw1ryw2[loff + l] + ryw1rxw2[loff + l])
                                                + xigrad_yz * xigrad_y * zomega * (ryw1ryw2[loff + l] + ryw1ryw2[loff + l])
                                                + xigrad_yz * xigrad_z * zomega * (rzw1ryw2[loff + l] + ryw1rzw2[loff + l])
                                                + xigrad_zz * xigrad_x * zomega * (rxw1rzw2[loff + l] + rzw1rxw2[loff + l])
                                                + xigrad_zz * xigrad_y * zomega * (ryw1rzw2[loff + l] + rzw1ryw2[loff + l])
                                                + xigrad_zz * xigrad_z * zomega * (rzw1rzw2[loff + l] + rzw1rzw2[loff + l]);
                        

                        double twelthfourth_gam =       xigrad_x * xomega * rxw1rxw2[loff + l]
                                                    +   xigrad_x * yomega * rxw1ryw2[loff + l]
                                                    +   xigrad_x * zomega * rxw1rzw2[loff + l]
                                                    +   xigrad_y * xomega * ryw1rxw2[loff + l]
                                                    +   xigrad_y * yomega * ryw1ryw2[loff + l]
                                                    +   xigrad_y * zomega * ryw1rzw2[loff + l]
                                                    +   xigrad_z * xomega * rzw1rxw2[loff + l]
                                                    +   xigrad_z * yomega * rzw1ryw2[loff + l]
                                                    +   xigrad_z * zomega * rzw1rzw2[loff + l];

                        double twelthfifth_gam =    (xigrad_x * grada_x[loff + l] + grada_x[loff + l] * xigrad_x) * rxw1rxw2[loff + l]
                                                +   (xigrad_x * grada_y[loff + l] + grada_x[loff + l] * xigrad_y) * rxw1ryw2[loff + l]
                                                +   (xigrad_x * grada_z[loff + l] + grada_x[loff + l] * xigrad_z) * rxw1rzw2[loff + l]
                                                +   (xigrad_y * grada_x[loff + l] + grada_y[loff + l] * xigrad_x) * rxw1rxw2[loff + l]
                                                +   (xigrad_y * grada_y[loff + l] + grada_y[loff + l] * xigrad_y) * rxw1ryw2[loff + l]
                                                +   (xigrad_y * grada_z[loff + l] + grada_y[loff + l] * xigrad_z) * rxw1rzw2[loff + l]
                                                +   (xigrad_z * grada_x[loff + l] + grada_z[loff + l] * xigrad_x) * rxw1rxw2[loff + l]
                                                +   (xigrad_z * grada_y[loff + l] + grada_z[loff + l] * xigrad_y) * rxw1ryw2[loff + l]
                                                +   (xigrad_z * grada_z[loff + l] + grada_z[loff + l] * xigrad_z) * rxw1rzw2[loff + l];


                        double twelfth = w * df0010[loff + l] * twelthfirst;

                        twelfth += w * (df0020[loff + l] +  df0011[loff + l]) * (twelththird + twelthsecond);

                        twelfth += w * (df0030[loff + l] + 2.0 * df0021[loff + l] + df0012[loff + l]) * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;

                        twelfth += w * (df00101[loff + l] + df00011[loff + l]) * ngrada[loff + l] * twelthsecond;

                        twelfth += w * (df00101[loff + l] + df00011[loff + l]) * xigrad_dot_omega * rw1_dot_rw2;

                        twelfth += w * df00002[loff + l] * grad_dot_omega * rw1_dot_rw2;

                        twelfth += w *  (df00101[loff + l] + df00011[loff + l]) * twelthfourth_gam;

                        twelfth += w * df00101[loff + l] * ngrada[loff + l] * twelththird;

                        twelfth += w * df00002[loff + l] * ngrada[loff + l] * twelthfourth_gam;
                        
                        twelfth += w * (df00201[loff + l] + df00111[loff + l]) * twelthfifth_gam * xigrad_dot_omega;

                        twelfth += w * df00102[loff + l] * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;
                        
                        twelfth += w * (df00201[loff + l] + 2 * df00111[loff + l] + df00021[loff + l]) * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;

                        twelfth += w * (df00102[loff + l] + df00011[loff + l]) * ngrada[loff + l] * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;

                        twelfth += w * df00003[loff + l] * ngrada[loff + l] * ngrada[loff + l] * grad_dot_omega;

                        fvxc += first;
                        fvxc += second;
                        fvxc += third;
                        fvxc += fourth;
                        fvxc += fifth;
                        fvxc += seventh;
                        fvxc += eighth ;
                        fvxc += ninth;
                        fvxc += tenth;
                        fvxc += twelfth;
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
CXCIntegrator:: _compEnergyAndDensityUnrestricted(const CXCGradientGrid& xcGradientGrid,
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

    if ((densityGrid.getDensityGridType() == dengrid::limb) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+:nele) reduction(+:xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * rhoa[i];
            
            xcene += gw[i] * efunc[i];
        }
    }

    if ((densityGrid.getDensityGridType() == dengrid::lima) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+:nele) reduction(+:xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * (rhob[i]);
            
            xcene += gw[i] * efunc[i];
        }
    }
    
    if ((densityGrid.getDensityGridType() == dengrid::ab) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+:nele) reduction(+:xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * (rhoa[i]+rhob[i]);
            
            xcene += gw[i] * efunc[i];
        }
    }
    
    
    return std::make_tuple(xcene, nele);
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

void
CXCIntegrator::_addSubMatrixUnRest(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                   const CDenseMatrix&      subMatrixa,
                                   const CDenseMatrix&      subMatrixb,
                                   const CGtoBlock&         braGtoBlock,
                                   const CGtoBlock&         ketGtoBlock) const
{
    // set up dimensions of Kohn-Sham matrix
    
    auto ncols = aoKohnShamMatrix->getNumberOfColumns();
    
    // set up submatrix dimensions

    auto bdim = subMatrixa.getNumberOfRows();

    auto kdim = subMatrixa.getNumberOfColumns();
    
    // set up submatrix offsets
    
    auto boff = (braGtoBlock.getIdentifiers(0))[0];

    auto koff = (ketGtoBlock.getIdentifiers(0))[0];
    
    // set up pointers to data
    
    auto ksmata = aoKohnShamMatrix->getKohnSham(false);

    auto ksmatb = aoKohnShamMatrix->getKohnSham(true);

    auto bkmata = subMatrixa.values();

    auto bkmatb = subMatrixb.values();
    
    // loop over submatrix
    
    for (int32_t i = 0; i < bdim; i++)
    {
        auto ioff = (boff + i) * ncols;
        
        for (int32_t j = 0; j < kdim; j++)
        {
            ksmata[ioff + koff + j] += bkmata[i * kdim + j];
            ksmatb[ioff + koff + j] += bkmatb[i * kdim + j];
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
