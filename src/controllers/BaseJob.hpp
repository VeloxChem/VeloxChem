//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef BaseJob_hpp
#define BaseJob_hpp

#include <cstdint>

#include "mpi.h"

#include "InputData.hpp"
#include "OutputStream.hpp"
#include "ExecMode.hpp"

/**
 Class CBaseJob defines base class for description of single job.
 
 @author Z. Rinkevicius
 */
class CBaseJob
{
protected:

    /**
     The job state: true - no errors, false otherwise
     */
    bool _state;

    /**
     The rank of MPI process associated with base job object.
     */
    int32_t _globRank;

    /**
     The total number of MPI processes associated with base job object.
     */
    int32_t _globNodes;
    
    /**
     The execution mode of job.
     */
    execmode _runMode;
    
public:

    /**
     Creates a basic job object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode of job.
     */
    CBaseJob(const int32_t  globRank,
             const int32_t  globNodes,
             const execmode runMode);
    
    /**
     Destroys a basic job object.
     */
    virtual ~CBaseJob() { };

    /**
     Get a state of basic job object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Sets parameters of basic job.

     @param pathToBasisSets the path to basis sets library.
     @param pathToForceFields the path to force fields library.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    virtual void set(const std::string&   pathToBasisSets,
                     const std::string&   pathToForceFields,
                     const CInputData&    inputData,
                           COutputStream& oStream) = 0;

    /**
     Executes basic job.

     @param comm the MPI communicator.
     @param oStream the output stream.
     */
    virtual void run(COutputStream& oStream,
                     MPI_Comm       comm) = 0;
};

#endif /* BaseJob_hpp */
