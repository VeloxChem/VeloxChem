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
     @return the tuple <number of electrons, exchange-correlation energy>.
     */
    std::tuple<double, double> integrate(const CAODensityMatrix& aoDensityMatrix,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CMolecularGrid&   molecularGrid,
                                         const std::string&      xcFuncLabel) const;
};

#endif /* XCIntegrator_hpp */
