//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <cstdint>
#include <string>
#include <tuple>

#include "mpi.h"

#include "AODensityMatrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DensityGrid.hpp"
#include "GtoContainer.hpp"
#include "AOFockMatrix.hpp"

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
     Computes exchange-correlation contribution to Kohn-Sham matrix from restricted density.

     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compRestrictedContribution(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                     const CGtoContainer*     gtoContainer,
                                     const CXCGradientGrid&   xcGradientGrid,
                                     const CDensityGrid&      densityGrid,
                                     const CMolecularGrid&    molecularGrid,
                                     const xcfun              xcFunctional) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham matrix from restricted density
     using matrix x matrix scheme.
     
     @param aoKohnShamMatrix the Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compRestrictedContributionM3(      CAOKohnShamMatrix& aoKohnShamMatrix,
                                       const CGtoContainer*     gtoContainer,
                                       const CXCGradientGrid&   xcGradientGrid,
                                       const CDensityGrid&      densityGrid,
                                       const CMolecularGrid&    molecularGrid,
                                       const xcfun              xcFunctional) const;
    
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
     Computes exchange-correlation contribution to Kohn-Sham matrix for batches of GTOs blocks.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param gtoContainer the container of GTOs blocks.
     @param xcGradientGrid the exchange-correlation functional grid.
     @param densityGrid the density grid.
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
                                                const CDensityGrid*      densityGrid,
                                                const double*            gridCoordinatesX,
                                                const double*            gridCoordinatesY,
                                                const double*            gridCoordinatesZ,
                                                const double*            gridWeights,
                                                const int32_t            gridOffset,
                                                const int32_t            nGridPoints,
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
     Computes number of electrons and exchange-correlation energy for batches of GTOs blocks.

     @param xcElectrons the number of electrons.
     @param xcEnergy the exchange-correlation energy.
     @param xcGradientGrid the exchange-correlation functional grid.
     @param densityGrid the density grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the grid offset.
     @param nGridPoints the number of grid points,
     */
    void _compRestrictedEnergyForBatchOfGridPoints(      double&          xcElectrons,
                                                         double&          xcEnergy,
                                                   const CXCGradientGrid* xcGradientGrid,
                                                   const CDensityGrid*    densityGrid,
                                                   const double*          gridWeights,
                                                   const int32_t          gridOffset,
                                                   const int32_t          nGridPoints) const;
    
    /**
     Computes exchange-correlation contribution to Kohn-Sham for GTOs blocks.

     @param aoKohnShamMatrix the pointer to Kohn-Sham matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param xcGradientGrid the exchange-correlation functional grid.
     @param densityGrid the density grid.
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
                                        const CDensityGrid*      densityGrid,
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
     Computes exchange-correlation functional contribution from restricted density to pair of spherical contracted GTOs.

     @param pairValues the vector of partial Kohn-Sham elements for contracted GTOs pairs.
     @param braGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid for bra side.
     @param ketGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid for ket side.
     @param braAngularComponents the number of angular components on bra side.
     @param ketAngularComponents the number of angular components on ket side.
     @param xcGradientGrid the exchange-correlation functional grid.
     @param densityGrid the density grid.
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
                                            const CDensityGrid*        densityGrid,
                                            const double*              gridWeights,
                                            const int32_t              gridOffset,
                                            const xcfun                xcFunctional) const;
    
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
     Gets number of rows in grid points matrix.

     @return the number of rows.
     */
    int32_t _getNumberOfGridRows() const;
    
    /**
     Computes GTOs values matrix for given block of grid points.

     @param gtoMatrix the GTOs values matrix.
     @param gtoContainer the container of GTOs blocks.
     @param molecularGrid the molecular grid.
     @param gridOffset the offset in molecular grid.
     @param nGridPoints the number of grid points. 
     */
    void _compGtosMatrixForLDA(      CDenseMatrix&   gtoMatrix,
                               const CGtoContainer*  gtoContainer,
                               const CMolecularGrid& molecularGrid,
                               const int32_t         gridOffset,
                               const int32_t         nGridPoints) const;
    
    
    /**
     Computes GTOs values for batch of grid points.

     @param gtoMatrix the pointer to GTOs values matrix.
     @param nRows the number of rows in GTOs values matrix.
     @param nColumns the number of columns in GTOs values matrix.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _compGtosValuesForLDA(      double*        gtoMatrix,
                               const int32_t        nRows,
                               const int32_t        nColumns,
                               const CGtoContainer* gtoContainer,
                               const double*        gridCoordinatesX,
                               const double*        gridCoordinatesY,
                               const double*        gridCoordinatesZ,
                               const int32_t        gridOffset,
                               const int32_t        nGridPoints) const;
    
    /**
     Computes XC contribution scaled GTOs matrix for spin-restricted LDA case.

     @param ketGtoMatrix the XC contribution scaled GTOs matrix.
     @param braGtoMatrix the GTOs matrix.
     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param gridWeights the pointer to grid weights.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _compRestrictedVXCMatrixForLDA(      CDenseMatrix&    ketGtoMatrix,
                                        const CDenseMatrix&    braGtoMatrix,
                                        const CXCGradientGrid& xcGradientGrid,
                                        const double*          gridWeights,
                                        const int32_t          gridOffset,
                                        const int32_t          nGridPoints) const;
    
    
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
     Integrates exchnage-correlation functional contribution to zero order Kohn-Sham matrix
     using matrix x matrix multiplication scheme.
     
     @param aoDensityMatrix the AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix integrate_m3(const CAODensityMatrix& aoDensityMatrix,
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
};

#endif /* XCIntegrator_hpp */
