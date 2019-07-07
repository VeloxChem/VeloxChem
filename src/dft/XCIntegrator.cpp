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

std::tuple<double, double>
CXCIntegrator::integrate(const CAODensityMatrix& aoDensityMatrix,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel) const
{
    // initialize number of electrons
    
    double felec = 0.0;
    
    // initialize exchange-correlation functional
    
    double xcene = 0.0;
    
    // generate reference density grid
    
    CDensityGridDriver dgdrv(_locComm);
    
    auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, xcfun::lda); // FIX ME: generic functional
    
    // set up number of density matrices
    
    auto ndmat = aoDensityMatrix.getNumberOfDensityMatrices();
    
    if (aoDensityMatrix.isRestricted())
    {
        // molecular and density grids
        
        std::vector<CMolecularGrid> mgrids(ndmat, molecularGrid);
        
        std::vector<CDensityGrid> dgrids(ndmat, CDensityGrid(molecularGrid.getNumberOfGridPoints(), 1, xcfun::lda, dengrid::ab));
        
        // generate screened density and molecular grids
        
        refdengrid.setScreenedGrids(dgrids, mgrids, _thresholdOfDensity, xcfun::lda); 
        
        // perform integration for spin restricted densities
        
        for (int32_t i = 0; i < ndmat; i++)
        {
            printf("Number of grid points in grid %i: %i -> %i", i, molecularGrid.getNumberOfGridPoints(), mgrids[i].getNumberOfGridPoints()); 
        }
    }
    else
    {
        // perform integration for spin unrestricted densities
        
        // FIX ME: implement partitioning and screening for density grid
    }
    
    // computes functional and it's first derrivates for density grid
    
    return std::make_tuple(felec, xcene);
}
