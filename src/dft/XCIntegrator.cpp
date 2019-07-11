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

void 
CXCIntegrator::integrate(const CAODensityMatrix& aoDensityMatrix,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel) const
{
    // parse exchange-correlation functional data
    
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
    
    // generate reference density grid
    
    CDensityGridDriver dgdrv(_locComm);
    
    auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());
    
    // set up number of density matrices
    
    auto ndmat = aoDensityMatrix.getNumberOfDensityMatrices();
    
    if (aoDensityMatrix.isRestricted())
    {
        // molecular and density grids
        
        std::vector<CMolecularGrid> mgrids(ndmat, molecularGrid);
        
        std::vector<CDensityGrid> dgrids(ndmat, CDensityGrid(molecularGrid.getNumberOfGridPoints(), 1,
                                                             fvxc.getFunctionalType(), dengrid::ab));
        
        // generate screened density and molecular grids
        
        refdengrid.setScreenedGrids(dgrids, mgrids, _thresholdOfDensity, fvxc.getFunctionalType());
        
        // perform integration for spin restricted densities
        
        for (int32_t i = 0; i < ndmat; i++)
        {
            // allocate XC gradient grid
            
            CXCGradientGrid vxcgrid (mgrids[i].getNumberOfGridPoints(), dgrids[i].getDensityGridType(), fvxc.getFunctionalType());
            
            // compute exchange-correlation functional first derrivatives
            
            fvxc.compute(vxcgrid, dgrids[i]);
            
            // compute Kohn-Sham matrix
            
            
        }
    }
    else
    {
        // perform integration for spin unrestricted densities
        
        // FIX ME: implement partitioning and screening for density grid
    }
}
