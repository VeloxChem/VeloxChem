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

#include "mpi.h"

#include "AODensityMatrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "XCGradientGrid.hpp"
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
     @param xcGradientGrid the exchange-correlation functional grid.
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
};

#endif /* XCIntegrator_hpp */
