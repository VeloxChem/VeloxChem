//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GtoFunc.hpp"

#include "GtoFuncForLDA.hpp"
#include "GtoFuncForGGA.hpp"
#include "AngularMomentum.hpp"
#include "OMPTasks.hpp"

namespace gtorec {  // gtorec namespace
    
    void
    computeGtosValuesForLDA(      CMemBlock2D<double>& gtoValues,
                            const CGtoContainer*       gtoContainer,
                            const double*              gridCoordinatesX,
                            const double*              gridCoordinatesY,
                            const double*              gridCoordinatesZ,
                            const int32_t              gridOffset,
                            const int32_t              gridBlockPosition,
                            const int32_t              nGridPoints)
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
                                               gridBlockPosition + gridOffset, bgtos, j, xcfun::lda);
                
                // distribute j-th GTO values into grid values matrix
                
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto bgaos = bspherbuff.data(k);
                    
                    auto gvals = gtoValues.data((bgtos.getIdentifiers(k))[j]);
                    
                    #pragma omp simd aligned(bgaos, gvals: VLX_ALIGN)
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        gvals[l] = bgaos[l];
                    }
                }
            }
        }
    }
    
    void
    computeGtosValuesForGGA(      CMemBlock2D<double>& gtoValues,
                                  CMemBlock2D<double>& gtoValuesX,
                                  CMemBlock2D<double>& gtoValuesY,
                                  CMemBlock2D<double>& gtoValuesZ,
                            const CGtoContainer*       gtoContainer,
                            const double*              gridCoordinatesX,
                            const double*              gridCoordinatesY,
                            const double*              gridCoordinatesZ,
                            const int32_t              gridOffset,
                            const int32_t              gridBlockPosition,
                            const int32_t              nGridPoints)
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
            
            auto nvcomp = xcfun_components(xcfun::gga);
            
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
                                               gridBlockPosition + gridOffset, bgtos, j, xcfun::gga);
                
                // distribute j-th GTO values into grid values matrix
                
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto bgaos = bspherbuff.data(4 * k);
                    
                    auto bgaox = bspherbuff.data(4 * k + 1);
                    
                    auto bgaoy = bspherbuff.data(4 * k + 2);
                    
                    auto bgaoz = bspherbuff.data(4 * k + 3);
                    
                    auto idx = (bgtos.getIdentifiers(k))[j];
                    
                    auto gvals = gtoValues.data(idx);
                    
                    auto gvalx = gtoValuesX.data(idx);
                    
                    auto gvaly = gtoValuesY.data(idx);
                    
                    auto gvalz = gtoValuesZ.data(idx);
                    
                    #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz: VLX_ALIGN)
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        gvals[l] = bgaos[l];
                        
                        gvalx[l] = bgaox[l];
                        
                        gvaly[l] = bgaoy[l];
                        
                        gvalz[l] = bgaoz[l];
                    }
                }
            }
        }
    }
    
    void
    computeGtosValuesForGGA(      CMemBlock2D<double>& gtoValues,
                                  CMemBlock2D<double>& gtoValuesX,
                                  CMemBlock2D<double>& gtoValuesY,
                                  CMemBlock2D<double>& gtoValuesZ,
                            const CGtoBlock&           gtoBlock,
                            const double*              gridCoordinatesX,
                            const double*              gridCoordinatesY,
                            const double*              gridCoordinatesZ,
                            const int32_t              gridOffset,
                            const int32_t              gridBlockPosition,
                            const int32_t              nGridPoints)
    {
        // angular momentum data for bra and ket
        
        auto bang = gtoBlock.getAngularMomentum();
        
        // set up Cartesian GTOs buffer
        
        auto nvcomp = xcfun_components(xcfun::gga);
        
        auto bncart = angmom::to_CartesianComponents(bang);
        
        auto bcartbuff = (bang > 0) ? CMemBlock2D<double>(nGridPoints, nvcomp * bncart) : CMemBlock2D<double>();
        
        // set up spherical GTOs buffer
        
        auto bnspher = angmom::to_SphericalComponents(bang);
        
        CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);
        
        // set up offset of GTOs buffer
        
        auto boff = (gtoBlock.getIdentifiers(0))[0];
        
        // loop over contracted GTOs
        
        for (int32_t i = 0; i < gtoBlock.getNumberOfContrGtos(); i++)
        {
            // compute i-th GTO values on batch of grid points
            
            gtorec::computeGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                           gridBlockPosition + gridOffset, gtoBlock, i, xcfun::gga);
            
            // distribute i-th GTO values on batch of grid points
            
            for (int32_t j = 0; j < bnspher; j++)
            {
                auto bgaos = bspherbuff.data(4 * j);
                
                auto bgaox = bspherbuff.data(4 * j + 1);
                
                auto bgaoy = bspherbuff.data(4 * j + 2);
                
                auto bgaoz = bspherbuff.data(4 * j + 3);
                
                auto idx = (gtoBlock.getIdentifiers(j))[i] - boff;
                
                auto gvals = gtoValues.data(idx);
                
                auto gvalx = gtoValuesX.data(idx);
                
                auto gvaly = gtoValuesY.data(idx);
                
                auto gvalz = gtoValuesZ.data(idx);
                
                #pragma omp simd aligned(bgaos, bgaox, bgaoy, bgaoz, gvals, gvalx, gvaly, gvalz: VLX_ALIGN)
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    gvals[k] = bgaos[k];
                    
                    gvalx[k] = bgaox[k];
                    
                    gvaly[k] = bgaoy[k];
                    
                    gvalz[k] = bgaoz[k];
                }
            }
        }
    }
    
    void
    computeGtosMatrixForLDA(      CDenseMatrix&   gtoMatrix,
                            const CGtoContainer*  gtoContainer,
                            const CMolecularGrid& molecularGrid,
                            const int32_t         gridOffset,
                            const int32_t         nGridPoints)
    {
        // set up OMP tasks
        
        COMPTasks omptaks(5);
        
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
        
        // generate density on grid points
        
        #pragma omp parallel shared(pgaos, tbsizes, tbpositions, ntasks, mgx, mgy, mgz)
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
                        gtorec::computeGtosValuesForLDA(pgaos, gtoContainer, mgx, mgy, mgz,
                                                        gridOffset, tbposition, tbsize);
                    }
                }
            }
        }
    }
    
    void
    computeGtosMatrixForGGA(      CDenseMatrix&   gtoMatrix,
                                  CDenseMatrix&   gtoMatrixX,
                                  CDenseMatrix&   gtoMatrixY,
                                  CDenseMatrix&   gtoMatrixZ,
                            const CGtoContainer*  gtoContainer,
                            const CMolecularGrid& molecularGrid,
                            const int32_t         gridOffset,
                            const int32_t         nGridPoints)
    {
        // set up OMP tasks
        
        COMPTasks omptaks(5);
        
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
        
        auto pgaox = gtoMatrixX.values();
        
        auto pgaoy = gtoMatrixY.values();
        
        auto pgaoz = gtoMatrixZ.values();
        
        // generate density on grid points
        
        #pragma omp parallel shared(pgaos, pgaox, pgaoy, pgaoz, tbsizes, tbpositions, ntasks, mgx, mgy, mgz)
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
                        gtorec::computeGtosValuesForGGA(pgaos, pgaox, pgaoy, pgaoz, gtoContainer, mgx, mgy, mgz,
                                                        gridOffset, tbposition, tbsize);
                    }
                }
            }
        }
    }
    
    void
    computeGtosValuesForLDA(      double*        gtoMatrix,
                            const CGtoContainer* gtoContainer,
                            const double*        gridCoordinatesX,
                            const double*        gridCoordinatesY,
                            const double*        gridCoordinatesZ,
                            const int32_t        gridOffset,
                            const int32_t        gridBlockPosition,
                            const int32_t        nGridPoints)
    {
        // local copy of GTOs containers
        
        auto gtovec = CGtoContainer(*gtoContainer);
        
        // set up number of AOs
        
        auto naos = gtovec.getNumberOfAtomicOrbitals();
        
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
                                               gridBlockPosition + gridOffset, bgtos, j, xcfun::lda);
                
                // distribute j-th GTO values into grid values matrix
                
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto idx = (bgtos.getIdentifiers(k))[j];
                    
                    auto bgaos = bspherbuff.data(k);
                    
                    auto loff = gridBlockPosition * naos;
                    
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        gtoMatrix[loff + l * naos + idx] = bgaos[l];
                    }
                }
            }
        }
    }
    
    void
    computeGtosValuesForGGA(      double*        gtoMatrix,
                                  double*        gtoMatrixX,
                                  double*        gtoMatrixY,
                                  double*        gtoMatrixZ,
                            const CGtoContainer* gtoContainer,
                            const double*        gridCoordinatesX,
                            const double*        gridCoordinatesY,
                            const double*        gridCoordinatesZ,
                            const int32_t        gridOffset,
                            const int32_t        gridBlockPosition,
                            const int32_t        nGridPoints)
    {
        // local copy of GTOs containers
        
        auto gtovec = CGtoContainer(*gtoContainer);
        
        // set up number of AOs
        
        auto naos = gtovec.getNumberOfAtomicOrbitals();
        
        // loop over GTOs container data
        
        for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
        {
            auto bgtos = gtovec.getGtoBlock(i);
            
            // angular momentum data for bra and ket
            
            auto bang = bgtos.getAngularMomentum();
            
            // set up Cartesian GTOs buffer
            
            auto nvcomp = xcfun_components(xcfun::gga);
            
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
                                               gridBlockPosition + gridOffset, bgtos, j, xcfun::gga);
                
                // distribute j-th GTO values into grid values matrix
                
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto idx = (bgtos.getIdentifiers(k))[j];
                    
                    auto bgaos = bspherbuff.data(4 * k);
                    
                    auto bgaox = bspherbuff.data(4 * k + 1);
                    
                    auto bgaoy = bspherbuff.data(4 * k + 2);
                    
                    auto bgaoz = bspherbuff.data(4 * k + 3);
                    
                    auto loff = gridBlockPosition * naos;
                    
                    for (int32_t l = 0; l < nGridPoints; l++)
                    {
                        gtoMatrix[loff + l * naos + idx] = bgaos[l];
                        
                        gtoMatrixX[loff + l * naos + idx] = bgaox[l];
                        
                        gtoMatrixY[loff + l * naos + idx] = bgaoy[l];
                        
                        gtoMatrixZ[loff + l * naos + idx] = bgaoz[l]; 
                    }
                }
            }
        }
    }
    
    void
    computeGtoValuesOnGrid(      CMemBlock2D<double>& spherGtoGridBuffer,
                                 CMemBlock2D<double>& cartGtoGridBuffer,
                           const double*              gridCoordinatesX,
                           const double*              gridCoordinatesY,
                           const double*              gridCoordinatesZ,
                           const int32_t              gridOffset,
                           const CGtoBlock&           gtoBlock,
                           const int32_t              iContrGto,
                           const xcfun                xcFunctional)
    {
        // set up angular momentum
        
        auto mang = gtoBlock.getAngularMomentum();
        
        // local density approximation
        
        if (xcFunctional == xcfun::lda)
        {
            switch (mang)
            {
                    // s-type GTOs on grid
                    
                case 0:
                    ldarec::compGtoValuesForS(spherGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // p-type GTOs on grid
                    
                case 1:
                    ldarec::compGtoValuesForP(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // d-type GTOs on grid
                    
                case 2:
                    ldarec::compGtoValuesForD(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // FIX ME: implement l > 2 cases
                    
                default:
                    break;
            }
            
            return;
        }
        
        // generalized gradient approximation
        
        if (xcFunctional == xcfun::gga)
        {
            switch (mang)
            {
                    // s-type GTOs on grid
                    
                case 0:
                    ggarec::compGtoValuesForS(spherGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // p-type GTOs on grid
                    
                case 1:
                    ggarec::compGtoValuesForP(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // d-type GTOs on grid
                    
                case 2:
                    ggarec::compGtoValuesForD(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // FIX ME: implement l > 2 cases
                    
                default:
                    break;
            }
            
            return;
        }
        
        // meta generalized gradient approximation
        
        if (xcFunctional == xcfun::mgga)
        {
            // FIX ME: MGGA case
            
            return;
        }
    }
    
}  // namespace gtorec
