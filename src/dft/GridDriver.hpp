//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GridDriver_hpp
#define GridDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "OutputStream.hpp"
#include "ExecMode.hpp"

/**
 Class CGridDriver generates grid points data for usage in numerical
 integration.
 
 @author Z. Rinkevicius
 */
class CGridDriver
{
    /**
     The accuracy level of grid (from 1 to 6).
     */
    int32_t _gridLevel;

    /**
     The rank of associated global MPI process.
     */
    int32_t _globRank;

    /**
     The total number of global MPI processes.
     */
    int32_t _globNodes;

    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The flag for local execution mode.
     */
    bool _isLocalMode;

    /**
     The threshold of weights screening.
     */
    double _thresholdOfWeight;
    
    /**
     The execution mode of grid driver object.
     */
    execmode _runMode;

    /**
     Determines number of radial grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of radial points.
     */
    int32_t _getNumberOfRadialPoints(const int32_t idElemental) const;

    /**
     Determines number of angular grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of angular points.
     */
    int32_t _getNumberOfAngularPoints(const int32_t idElemental) const;

public:

    /**
     Creates a grid driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode. 
     @param comm the MPI communicator.
     */
    CGridDriver(const int32_t globRank, const int32_t globNodes,
                execmode runMode, MPI_Comm comm);

    /**
     Destroys a grid driver object.
     */
    ~CGridDriver();

    /**
     Sets accuracy level for grid generation. Level: 1-6, where 1 is coarse
     grid, 5 is ultrafine grid, 6 special benchmarking grid.

     @param gridLevel the accuracy level of generated grid.
     @param comm the MPI communicator.
     */
    void setLevel(const int32_t gridLevel, MPI_Comm comm);
};

#endif /* GridDriver_hpp */
