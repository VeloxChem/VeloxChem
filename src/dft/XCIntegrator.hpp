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

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <cstdint>
#include <string>
#include <tuple>

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "XCCubicHessianGrid.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridQuad.hpp"
#include "GtoContainer.hpp"
#include "AOFockMatrix.hpp"
#include "MemBlock.hpp"
#include "OverlapMatrix.hpp"

/**
 Class CXCIntegrator implements exchange-correlation functional and it's derrivatives integraion.
 
 @author Z. Rinkevicius
 */
class CXCIntegrator
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;
    
    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;
    
    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;
    
    /**
     The threshold of density screening.
     */
    double _thresholdOfDensity;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix for spin-restricted LDA case.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForLDAWithNL(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                 const CGtoContainer*     gtoContainer,
                                                 const CXCGradientGrid&   xcGradientGrid,
                                                 const CDensityGrid&      densityGrid,
                                                 const CMolecularGrid&    molecularGrid) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param iBraContrGto the index of contracted GTO on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForLDAWithNL(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                          const CGtoBlock&         braGtoBlock,
                                          const int32_t            iBraContrGto,
                                          const CGtoBlock&         ketGtoBlock,
                                          const CXCGradientGrid*   xcGradientGrid,
                                          const double*            gridCoordinatesX,
                                          const double*            gridCoordinatesY,
                                          const double*            gridCoordinatesZ,
                                          const double*            gridWeights,
                                          const int32_t            nGridPoints) const;
    
    /**
     Checks if shell pair has significant contribution to Kohn-Sham matrix.
     
     @param braGtoBlock the GTOs block on bra side.
     @param iBraContrGto the index of contracted GTO on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iKetContrGto the index of contracted GTO on ket side.
     @return true if shell pair is significant, false otherwise.
     */
    bool _isSignificantShellPair(const CGtoBlock&      braGtoBlock,
                                 const int32_t         iBraContrGto,
                                 const CGtoBlock&      ketGtoBlock,
                                 const int32_t         iKetContrGto) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix for spin-restricted LDA case.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCGradientGrid&   xcGradientGrid,
                                           const CDensityGrid&      densityGrid,
                                           const CMolecularGrid&    molecularGrid) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix for spin-restricted LDA case.

     @param aoKohnShamMatrix the spin-restricted exchange-correlation contribution to perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed densities grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCHessianGrid&    xcHessianGrid,
                                           const CDensityGrid&      rwDensityGrid,
                                           const CMolecularGrid&    molecularGrid) const;
   /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix for spin-restricted LDA case.

     @param aoKohnShamMatrix the spin-restricted exchange-correlation contribution to perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed densities grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForLda(      CAOKohnShamMatrix&   aoKohnShamMatrix,
                                           const CGtoContainer*       gtoContainer,
                                           const CXCHessianGrid&      xcHessianGrid,
                                           const CXCCubicHessianGrid& xcCubicHessianGrid,
                                           const CDensityGridQuad&        rwDensityGrid,
                                           const CDensityGrid&        rw2DensityGrid,
                                           const CMolecularGrid&      molecularGrid,
                                           const std::string&         quadMode) const;

    /**
        Computes exchange-correlation contribution to Kohn-Sham matrix for spin-unrestricted LDA case.
        
        @param aoKohnShamMatrix the Kohn-Sham matrix.
        @param gtoContainer the container of GTOs blocks.
        @param xcGradientGrid the exchange-correlation functional gradient grid.
        @param densityGrid the density grid.
        @param molecularGrid the molecular grid.
        */

       void _compUnrestrictedContributionForLda(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                                const CGtoContainer*     gtoContainer,
                                                const CXCGradientGrid&   xcGradientGrid,
                                                const CDensityGrid&      densityGrid,
                                                const CMolecularGrid&    molecularGrid) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix for spin-restricted GGA case.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCGradientGrid&   xcGradientGrid,
                                           const CDensityGrid&      densityGrid,
                                           const CMolecularGrid&    molecularGrid) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix for spin-unrestricted GGA case.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     */
    void _compUnrestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                             const CGtoContainer*     gtoContainer,
                                             const CXCGradientGrid&   xcGradientGrid,
                                             const CDensityGrid&      densityGrid,
                                             const CMolecularGrid&    molecularGrid) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix for spin-restricted GGA case.
     
     @param aoKohnShamMatrix the spin-restricted exchange-correlation contribution to perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param gsDensityGrid the ground state density grid.
     @param rwDensityGrid the perturbed densities grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForGga(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                           const CGtoContainer*     gtoContainer,
                                           const CXCGradientGrid&   xcGradientGrid,
                                           const CXCHessianGrid&    xcHessianGrid,
                                           const CDensityGrid&      gsDensityGrid,
                                           const CDensityGrid&      rwDensityGrid,
                                           const CMolecularGrid&    molecularGrid) const;
    
   /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix for spin-restricted GGA case.
     
     @param aoKohnShamMatrix the spin-restricted exchange-correlation contribution to perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param gsDensityGrid the ground state density grid.
     @param rwDensityGrid the perturbed densities grid.
     @param molecularGrid the molecular grid.
     */
    void _compRestrictedContributionForGga(      CAOKohnShamMatrix&   aoKohnShamMatrix,
                                           const CGtoContainer*       gtoContainer,
                                           const CXCGradientGrid&     xcGradientGrid,
                                           const CXCHessianGrid&      xcHessianGrid,
                                           const CXCCubicHessianGrid& xcCubicHessianGrid,
                                           const CDensityGrid&        gsDensityGrid,
                                           const CDensityGridQuad&    rwDensityGrid,
                                           const CDensityGrid&        rw2DensityGrid,
                                           const CMolecularGrid&      molecularGrid) const;

    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                    const CGtoContainer*     gtoContainer,
                                    const CXCGradientGrid*   xcGradientGrid,
                                    const double*            gridCoordinatesX,
                                    const double*            gridCoordinatesY,
                                    const double*            gridCoordinatesZ,
                                    const double*            gridWeights,
                                    const int32_t            gridOffset,
                                    const int32_t            nGridPoints) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.
     
     @param aoKohnShamMatrix the perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                    const CGtoContainer*     gtoContainer,
                                    const CXCHessianGrid*    xcHessianGrid,
                                    const CDensityGrid*      rwDensityGrid,
                                    const double*            gridCoordinatesX,
                                    const double*            gridCoordinatesY,
                                    const double*            gridCoordinatesZ,
                                    const double*            gridWeights,
                                    const int32_t            gridOffset,
                                    const int32_t            nGridPoints) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.
     
     @param aoKohnShamMatrix the perturbed Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
                                    const std::string&         quadMode) const;

    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compUnrestrictedBatchForLda(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                      const CGtoContainer*     gtoContainer,
                                      const CXCGradientGrid*   xcGradientGrid,
                                      const double*            gridCoordinatesX,
                                      const double*            gridCoordinatesY,
                                      const double*            gridCoordinatesZ,
                                      const double*            gridWeights,
                                      const int32_t            gridOffset,
                                      const int32_t            nGridPoints) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the pointer to density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                    const CGtoContainer*     gtoContainer,
                                    const CXCGradientGrid*   xcGradientGrid,
                                    const CDensityGrid*      densityGrid,
                                    const double*            gridCoordinatesX,
                                    const double*            gridCoordinatesY,
                                    const double*            gridCoordinatesZ,
                                    const double*            gridWeights,
                                    const int32_t            gridOffset,
                                    const int32_t            nGridPoints) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the pointer to density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                    const CGtoBlock&         braGtoBlock,
                                    const CGtoBlock&         ketGtoBlock,
                                    const CXCGradientGrid*   xcGradientGrid,
                                    const CDensityGrid*      densityGrid,
                                    const double*            gridCoordinatesX,
                                    const double*            gridCoordinatesY,
                                    const double*            gridCoordinatesZ,
                                    const double*            gridWeights,
                                    const int32_t            gridOffset,
                                    const int32_t            nGridPoints) const;

    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the pointer to density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compUnrestrictedBatchForGga(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                      const CGtoBlock&         braGtoBlock,
                                      const CGtoBlock&         ketGtoBlock,
                                      const CXCGradientGrid*   xcGradientGrid,
                                      const CDensityGrid*      densityGrid,
                                      const double*            gridCoordinatesX,
                                      const double*            gridCoordinatesY,
                                      const double*            gridCoordinatesZ,
                                      const double*            gridWeights,
                                      const int32_t            gridOffset,
                                      const int32_t            nGridPoints) const;

    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param gsDensityGrid the pointer to ground state density grid.
     @param rwDensityGrid he pointer to perturbed density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForGGA(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
                                    const int32_t            nGridPoints) const;

    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param gsDensityGrid the pointer to ground state density grid.
     @param rwDensityGrid he pointer to perturbed density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points.
     */
    void _compRestrictedBatchForGGA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
                                    const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcGradient the pointer to exchange-correlation gradient grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                          CMemBlock<double>&   xcBuffer,
                                    const double*              xcGradient,
                                    const CMemBlock2D<double>& gtoValues,
                                    const double*              gridWeights,
                                    const int32_t              gridOffset,
                                    const int32_t              gridBlockPosition,
                                    const int32_t              nGridPoints) const;
    
    
    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer_alpha the exchange-correlation buffer.
     @param xcGradient_alpha the pointer to exchange-correlation gradient grid.
     @param xcBuffer_beta the exchange-correlation buffer.
     @param xcGradient_beta the pointer to the beta exchange-correlation gradient grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                            CMemBlock<double>&   xcBuffer_alpha,
                                            CMemBlock<double>&   xcBuffer_beta,
                                      const double*              xcGradient_alpha,
                                      const double*              xcGradient_beta,
                                      const CMemBlock2D<double>& gtoValues,
                                      const double*              gridWeights,
                                      const int32_t              gridOffset,
                                      const int32_t              gridBlockPosition,
                                      const int32_t              nGridPoints) const;
    
    /**
        Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
        for spin-unrestricted LDA case.

        @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
        @param xcBuffer_alpha the exchange-correlation buffer.
        @param xcGradient_alpha the pointer to exchange-correlation gradient grid.
        @param xcBuffer_beta the exchange-correlation buffer.
        @param xcGradient_beta the pointer to exchange-correlation gradient grid.
        @param gtoValues the pointer to GTOS values on grid.
        @param gridWeights the pointer to grid weights.
        @param gridOffset the offset of grids batch in density grid.
        @param gridBlockPosition the block position in grids batch.
        @param nGridPoints the number of grid points in grid block.
        */
       void _distUnrestrictedBatchForLdaA(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                                CMemBlock<double>&   xcBuffer_alpha,
                                                CMemBlock<double>&   xcBuffer_beta,
                                          const double*              xcGradient_alpha,
                                          const double*              xcGradient_beta,
                                          const CMemBlock2D<double>& gtoValues,
                                          const double*              gridWeights,
                                          const int32_t              gridOffset,
                                          const int32_t              gridBlockPosition,
                                          const int32_t              nGridPoints) const;
    
    /**
      Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
      for spin-unrestricted LDA case.

      @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
      @param xcBuffer_alpha the exchange-correlation buffer.
      @param xcGradient_alpha the pointer to exchange-correlation gradient grid.
      @param xcBuffer_beta the exchange-correlation buffer.
      @param xcGradient_beta the pointer to exchange-correlation gradient grid.
      @param gtoValues the pointer to GTOS values on grid.
      @param gridWeights the pointer to grid weights.
      @param gridOffset the offset of grids batch in density grid.
      @param gridBlockPosition the block position in grids batch.
      @param nGridPoints the number of grid points in grid block.
      */
      void _distUnrestrictedBatchForLdaB(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                               CMemBlock<double>&   xcBuffer_alpha,
                                               CMemBlock<double>&   xcBuffer_beta,
                                         const double*              xcGradient_alpha,
                                         const double*              xcGradient_beta,
                                         const CMemBlock2D<double>& gtoValues,
                                         const double*              gridWeights,
                                         const int32_t              gridOffset,
                                         const int32_t              gridBlockPosition,
                                         const int32_t              nGridPoints) const;
    
    /**
     Distributes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcHessianGrid the pointer to exchange-correlation hessian grid.
     @param rwDensityGrid the pointer to perturbed density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForLda(      CAOKohnShamMatrix*   aoKohnShamMatrix,
                                          CMemBlock<double>&   xcBuffer,
                                    const CXCHessianGrid*      xcHessianGrid,
                                    const CDensityGrid*        rwDensityGrid,
                                    const CMemBlock2D<double>& gtoValues,
                                    const double*              gridWeights,
                                    const int32_t              gridOffset,
                                    const int32_t              gridBlockPosition,
                                    const int32_t              nGridPoints) const;
    
    /**
     Distributes exchange-correlation contribution to perturbed Kohn-Sham matrix from batch of grid points
     for spin-restricted LDA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcHessianGrid the pointer to exchange-correlation hessian grid.
     @param rwDensityGrid the pointer to perturbed density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForLda(   CAOKohnShamMatrix*   aoKohnShamMatrix,
                                          CMemBlock<double>&   xcBuffer,
                                    const CXCHessianGrid*      xcHessianGrid,
                                    const CXCCubicHessianGrid* xcCubicHessianGrid,
                                    const CDensityGridQuad*        rwDensityGrid,
                                    const CDensityGrid*        rw2DensityGrid,
                                    const CMemBlock2D<double>& gtoValues,
                                    const double*              gridWeights,
                                    const int32_t              gridOffset,
                                    const int32_t              gridBlockPosition,
                                    const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcGradientGrid theexchange-correlation gradient grid.
     @param densityGrid the density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
                                    const int32_t              nGridPoints) const;
    
    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.
     
     @param subMatrix the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForGga(      CDenseMatrix&        subMatrix,
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
                                    const int32_t              nGridPoints) const;
    
     /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the partial Kohn-Sham matrix.
     @param subMatrix_beta the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGga(      CDenseMatrix&        subMatrix_alpha,
                                            CDenseMatrix&        subMatrix_beta,
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
                                      const int32_t              nGridPoints) const;

     /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the partial Kohn-Sham matrix.
     @param subMatrix_beta the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGgaA(      CDenseMatrix&        subMatrix_alpha,
                                             CDenseMatrix&        subMatrix_beta,
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
                                       const int32_t              nGridPoints) const;

     /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the partial Kohn-Sham matrix.
     @param subMatrix_beta the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGgaB( CDenseMatrix&        subMatrix_alpha,
                                            CDenseMatrix&        subMatrix_beta,
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
                                    const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.
     
     @param subMatrix the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param braGtoValues the pointer to GTOS values on grid on bra side.
     @param braGtoValuesX the GTOs gradient along X axis values buffer on bra side.
     @param braGtoValuesY the GTOs gradient along Y axis values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient along Z axis values buffer on bra side.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoValues the pointer to GTOS values on grid on ket side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param ketGtoBlock the GTOs block on ket side.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForGga(      CDenseMatrix&        subMatrix,
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
                                    const int32_t              nGridPoints) const;
    
    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the partial Kohn-Sham matrix.
     @param subMatrix_beta the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param braGtoValues the pointer to GTOS values on grid on bra side.
     @param braGtoValuesX the GTOs gradient along X axis values buffer on bra side.
     @param braGtoValuesY the GTOs gradient along Y axis values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient along Z axis values buffer on bra side.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoValues the pointer to GTOS values on grid on ket side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param ketGtoBlock the GTOs block on ket side.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGga(      CDenseMatrix&        subMatrix_alpha,
                                            CDenseMatrix&        subMatrix_beta,
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
                                      const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the partial Kohn-Sham matrix.
     @param subMatrix_beta the partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param braGtoValues the pointer to GTOS values on grid on bra side.
     @param braGtoValuesX the GTOs gradient along X axis values buffer on bra side.
     @param braGtoValuesY the GTOs gradient along Y axis values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient along Z axis values buffer on bra side.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoValues the pointer to GTOS values on grid on ket side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param ketGtoBlock the GTOs block on ket side.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGgaA(      CDenseMatrix&        subMatrix_alpha,
                                             CDenseMatrix&        subMatrix_beta,
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
                                       const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-unrestricted GGA case.
     
     @param subMatrix_alpha the alpha partial Kohn-Sham matrix.
     @param subMatrix_beta the beta partial Kohn-Sham matrix.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param densityGrid the pointer to density grid.
     @param braGtoValues the pointer to GTOS values on grid on bra side.
     @param braGtoValuesX the GTOs gradient along X axis values buffer on bra side.
     @param braGtoValuesY the GTOs gradient along Y axis values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient along Z axis values buffer on bra side.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoValues the pointer to GTOS values on grid on ket side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param ketGtoBlock the GTOs block on ket side.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distUnrestrictedBatchForGgaB(      CDenseMatrix&        subMatrix_alpha,
                                             CDenseMatrix&        subMatrix_beta,
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
                                       const int32_t              nGridPoints) const;

    /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param xcHessianGrid the pointer to exchange-correlation hessian grid.
     @param gsDensityGrid the pointer to ground state density grid.
     @param rwDensityGrid the pointer to perturbed density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
                                    const int32_t              nGridPoints) const;
   /**
     Distributes exchange-correlation contribution to Kohn-Sham matrix from batch of grid points
     for spin-restricted GGA case.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param xcBuffer the exchange-correlation buffer.
     @param xcGradientGrid the pointer to exchange-correlation gradient grid.
     @param xcHessianGrid the pointer to exchange-correlation hessian grid.
     @param gsDensityGrid the pointer to ground state density grid.
     @param rwDensityGrid the pointer to perturbed density grid.
     @param gtoValues the pointer to GTOS values on grid.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grids batch in density grid.
     @param gridBlockPosition the block position in grids batch.
     @param nGridPoints the number of grid points in grid block.
     */
    void _distRestrictedBatchForGga(      CAOKohnShamMatrix*   aoKohnShamMatrix,
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
                                    const int32_t              nGridPoints) const;


    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix from restricted density.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gsDensityGrid the ground state density grid.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compRestrictedContribution(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                     const CGtoContainer*     gtoContainer,
                                     const CXCGradientGrid&   xcGradientGrid,
                                     const CXCHessianGrid&    xcHessianGrid,
                                     const CDensityGrid&      rwDensityGrid,
                                     const CDensityGrid&      gsDensityGrid,
                                     const CMolecularGrid&    molecularGrid,
                                     const xcfun              xcFunctional) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham matrix for batches of GTOs blocks.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gsDensityGrid the ground state density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
            points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
           points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
            points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points,
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compRestrictedVXCForBatchOfGridPoints(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
                                                const xcfun              xcFunctional) const;
    
    /**
     Computes exchange-correlation contribution to perturbed Kohn-Sham for GTOs blocks.
     
     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gsDensityGrid the ground state density grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
            points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
            points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
            points.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points,
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compRestrictedVXCForGtoBlocks(      CAOKohnShamMatrix* aoKohnShamMatrix,
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
                                        const xcfun              xcFunctional) const;
    
    /**
     Computes exchange-correlation functional contribution from perturbed restricted density to pair of spherical contracted GTOs.
     
     @param pairValues the vector of partial Kohn-Sham elements for contracted GTOs pairs.
     @param braGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid for bra side.
     @param ketGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid for ket side.
     @param braAngularComponents the number of angular components on bra side.
     @param ketAngularComponents the number of angular components on ket side.
     @param xcGradientGrid the exchange-correlation functional grid.
     @param xcHessianGrid the exchange-correlation functional hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gsDensityGrid the ground state density grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the batch offset in vector grid points.
     @param xcFunctional the exchange-correlations functional type.
     */
    void _compRestrictedVXCValueForGtosPair(      CMemBlock<double>&   pairValues,
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
                                            const xcfun                xcFunctional) const;
    
    /**
     Distributes exchange-correlation functional contribution from pair of spherical contracted GTOs into Kohn-Sham matrix.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param pairValues the vector of partial Kohn-Sham elements for contracted GTOs pairs.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag indicating equality between bra and ket sides. 
     @param iBraContrGto the index of contracted GTO on bra side.
     @param iKetContrGto the index of contracted GTO on ket side.
     */
    void _distRestrictedVXCValues(      CAOKohnShamMatrix* aoKohnShamMatrix,
                                  const CMemBlock<double>& pairValues,
                                  const CGtoBlock&         braGtoBlock,
                                  const CGtoBlock&         ketGtoBlock,
                                  const bool               isBraEqualKet,
                                  const int32_t            iBraContrGto,
                                  const int32_t            iKetContrGto) const;
    
    /**
     Computes exchange-correlation energy and number of electrons for given density grid.

     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     @return the tuple (exchange-correlation energy, number of electrons). 
     */
    std::tuple<double, double> _compEnergyAndDensity(const CXCGradientGrid& xcGradientGrid,
                                                     const CDensityGrid&    densityGrid,
                                                     const CMolecularGrid&  molecularGrid) const;
    
    /**
     Computes exchange-correlation energy and number of electrons for given density grid.

     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     @return the tuple (exchange-correlation energy, number of electrons).
     */
    std::tuple<double, double>  _compEnergyAndDensityUnrestricted(const CXCGradientGrid& xcGradientGrid,
                                                                  const CDensityGrid&    densityGrid,
                                                                  const CMolecularGrid&  molecularGrid) const;
    
    /**
     Gets size of block in grid batch.
     
     @return the size of block in grid batch.
     */
    int32_t _getSizeOfBlock() const;
    
    /**
     Gets number of atomic orbitals included into accumulation buffer.

     @return the number of atomic orbitals.
     */
    int32_t _getNumberOfAOsInBuffer() const;
    
    /**
     Adds submatrix to Kohn-Sham matrix.

     @param aoKohnShamMatrix the AO Kohn-Sham matrix.
     @param subMatrix the submatrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _addSubMatrix(      CAOKohnShamMatrix* aoKohnShamMatrix,
                       const CDenseMatrix&      subMatrix,
                       const CGtoBlock&         braGtoBlock,
                       const CGtoBlock&         ketGtoBlock) const;

    /**
     Adds submatrix to Kohn-Sham matrix. Unrestricted

     @param aoKohnShamMatrix the AO Kohn-Sham matrix.
     @param subMatrix_alpha the alpha submatrix.
     @param subMatrix_beta the beta submatrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _addSubMatrixUnRest(      CAOKohnShamMatrix* aoKohnShamMatrix,
                             const CDenseMatrix&      subMatrix_alpha,
                             const CDenseMatrix&      subMatrix_beta,
                             const CGtoBlock&         braGtoBlock,
                             const CGtoBlock&         ketGtoBlock) const;

    /**
     Generates screening pattern for GTO values on grid.

     @param gtoValues the GTO values on grid.
     @return the screening pattern.
     */
    CMemBlock<int32_t> _getScreeningPattern(const CMemBlock2D<double>& gtoValues) const;
    
    /**
     Creates reduced integration grid with given screening pattern.

     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridWeights the pointer to grid weights.
     @param screeningPattern the screening pattern.
     @return the reduced integration grid.
     */
    CMemBlock2D<double> _getReducedGrid(const double*             gridCoordinatesX,
                                        const double*             gridCoordinatesY,
                                        const double*             gridCoordinatesZ,
                                        const double*             gridWeights,
                                        const CMemBlock<int32_t>& screeningPattern) const;
    
    /**
     Cretaes reduces GTO values on grid with given screening pattern.

     @param gtoValues the GTO values on grid.
     @param screeningPattern the screening pattern.
     @return the reduced GTO values on grid.
     */
    CMemBlock2D<double> _getReducedGtoValues(const CMemBlock2D<double>& gtoValues,
                                             const CMemBlock<int32_t>&  screeningPattern) const;
    
    /**
     Cretaes reduces spin restricted gradient grid with given screening pattern.

     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param screeningPattern the screening pattern.
     @return the redduced exchange-correlation functional gradient grid.
     */
    CMemBlock<double> _getReducedRestrictedGradient(const CXCGradientGrid*    xcGradientGrid,
                                                    const CMemBlock<int32_t>& screeningPattern) const;
   
public:
    
    /**
     Creates a XC integrator object using MPI info.
     
     @param comm the MPI communicator.
     */
    CXCIntegrator(MPI_Comm comm);
    
    /**
     Destroys a XC integrator object.
     */
    ~CXCIntegrator();
    
    /**
     Integrates exchnage-correlation functional contribution to zero order Kohn-Sham matrix.

     @param aoDensityMatrix the AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix integrate(const CAODensityMatrix& aoDensityMatrix,
                                const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const CMolecularGrid&   molecularGrid,
                                const std::string&      xcFuncLabel) const;
    
    /**
     Integrates exchnage-correlation functional contribution to first order Fock matrices and adds it to AO Fock matrix.
     
     @param aoFockMatrix the AO Fock matrix.
     @param rwDensityMatrix the perturbed AO density matrix object.
     @param gsDensityMatrix the ground state AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrate(      CAOFockMatrix&    aoFockMatrix,
                   const CAODensityMatrix& rwDensityMatrix,
                   const CAODensityMatrix& gsDensityMatrix,
                   const CMolecule&        molecule,
                   const CMolecularBasis&  basis,
                   const CMolecularGrid&   molecularGrid,
                   const std::string&      xcFuncLabel) const;
      
     /**
     Integrates exchnage-correlation functional contribution to second-order Fock matrices and adds it to AO Fock matrix.
     
     @param aoFockMatrix the AO Fock matrix.
     @param rwDensityMatrix the perturbed AO density matrix object.
     @param gsDensityMatrix the ground state AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrate(      CAOFockMatrix&    aoFockMatrix,
                   const CAODensityMatrix& rwDensityMatrix,
                   const CAODensityMatrix& rw2DensityMatrix,
                   const CAODensityMatrix& gsDensityMatrix,
                   const CMolecule&        molecule,
                   const CMolecularBasis&  basis,
                   const CMolecularGrid&   molecularGrid,
                   const std::string&      xcFuncLabel,
                   const std::string&      quadMode) const;
};

#endif /* XCIntegrator_hpp */
