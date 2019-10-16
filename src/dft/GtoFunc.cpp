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

namespace gtorec {  // gtorec namespace
    
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
