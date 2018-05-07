//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DensityGridDriver_hpp
#define DensityGridDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

class CDensityGridDriver
{
    int32_t _globRank;         // rank of associated global MPI process
    
    int32_t _globNodes;        // total number of global MPI processes
    
    int32_t _locRank;          // rank of associated local MPI process
    
    int32_t _locNodes;         // total number of local MPI processes
    
    bool _isLocalMode;         // flag for local execution mode
    
    double _thresholdOfDensity; // threshold of density screening
    
public:
    
    // CDensityGridDriver(const int32_t globRank, const int32_t globNodes,
    //                    MPI_Comm comm):
    //
    // Constructs a denisty grid griver using MPI info.
    //
    // Input:
    // globRank (int32_t)  - the rank of MPI process within domain of MPI
    //                       communicator.
    // globNodes (int32_t) - the number of MPI processes within domain of
    //                       MPI communicator.
    // comm (MPI_Comm)     - the MPI communicator.
    
    CDensityGridDriver(const int32_t  globRank,
                       const int32_t  globNodes,
                             MPI_Comm comm);
    
    // ~CDensityGridDriver():
    //
    // Destroys a density grid driver.
    
    ~CDensityGridDriver();
    
    // CDensityGrid generate(const CMolecule& molecule,
    //                       const CMolecularBasis& basis,
    //                       const CMolecularGrid& molGrid,
    //                       const xcfun            xcFunctional,
    //                       COutputStream& oStream, MPI_Comm comm):
    //
    // Generates partitioned density grid for given molecule and type of
    // exchange-correlation functional. Density grid generation is distributed
    // within domain of MPI communicator.
    //
    // Input/Output:
    // molecule (CMolecule&)     - the molecule.
    // basis (CMolecularBasis&)  - the molecular basis.
    // molGrid (CMolecularGrid&) - the molecular grid.
    // xcFunctional (xcfun)      - the type of exchange-correlation functional.
    // oStream (COutputStream&)  - the output stream.
    // comm (MPI_Comm)           - the MPI communicator.
    //
    // Output:
    // (CDensityGrid) - the density grid.
    
    void generate(const CMolecule&       molecule,
                  const CMolecularBasis& basis,
                  const CMolecularGrid&  molGrid,
                  const xcfun            xcFunctional,
                        COutputStream&   oStream,
                        MPI_Comm         comm);
};

#endif /* DensityGridDriver_hpp */
