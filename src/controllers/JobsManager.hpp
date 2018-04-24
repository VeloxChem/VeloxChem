//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef JobsManager_hpp
#define JobsManager_hpp

#include <cstdint>
#include <vector>

#include "InputData.hpp"
#include "OutputStream.hpp"
#include "BaseJob.hpp"
#include "ExecMode.hpp"

/**
 Class CJobsManager manages jobs list creation, ordering and execution workflow.
 
 @author Z. Rinkevicius
 */
class CJobsManager
{
    /**
     The state of jobs manager (true - no errors, false otherwise).
     */
    bool _state;

    /**
     The rank of MPI process associated with jobs manager object.
     */
    int32_t _globRank;

    /**
     The total number of MPI processes associated with jobs manager object.
     */
    int32_t _globNodes;
    
    /**
     The execution mode of jobs.
     */
    execmode _runMode;

    /**
     The vector of job objects.
     */
    std::vector<CBaseJob*> _listOfJobs;

    /**
     Creates vector of derrived job objects according to vector of given
     identifiers.

     @param listOfJobIds the vector of identifiers.
     */
    void _assignJobs(const std::vector<int32_t>& listOfJobIds);
    
    /**
     Assigns a run mode to jobs manager object.

     @param runMode the run mode.
     */
    void _assignRunMode(const execmode runMode); 

    /**
     Updates state of jobs manager object.

     @param state true to set normal state, false to abnormal state i.e.
            execution error is encountered.
     */
    void _updateState(const bool state);

public:

    /**
     Creates a jobs manager object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     */
    CJobsManager(const int32_t globRank,
                 const int32_t globNodes);

    /**
     Destroys a basic job object.
     */
    ~CJobsManager();

    /**
     Sets vector of job objects by reading input data object on master node and
     by creating up job objects on master and worker nodes according to encoded
     identifiers parsed from input data. Errors are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void setJobs(const CInputData&    inputData,
                       COutputStream& oStream);

    /**
     Executes job objects stored in vector of job objects using input data for
     additional parameters. Errors are printed to output stream.

     @param pathToBasisSets the path to basis set library.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void runJobs(const std::string&   pathToBasisSets,
                 const CInputData&    inputData,
                       COutputStream& oStream);

    /**
     Get a state of jobs manager.

     @return true if no execution errors, false otherwise.
     */
    bool getState() const;
};

#endif /* JobsManager_hpp */
