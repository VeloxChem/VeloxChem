//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "JobsManager.hpp"

#include "MpiFunc.hpp"
#include "JobsReader.hpp"
#include "JobType.hpp"
#include "SinglePointEnergy.hpp"
#include "OptimizationGeometry.hpp"

CJobsManager::CJobsManager(const int32_t globRank, const int32_t globNodes)

    : _state(true)

    , _globRank(globRank)

    , _globNodes(globNodes)

    , _runMode(execmode::cpu)
{

}

CJobsManager::~CJobsManager()
{
    for (size_t i = 0; i < _listOfJobs.size(); i++)
    {
        delete _listOfJobs[i];
    }
}

void CJobsManager::setJobs(const CInputData& inputData, COutputStream& oStream)
{
    std::vector<int32_t> idsjobs;

    if (_globRank == mpi::master())
    {
        CJobsReader rdrjobs;

        rdrjobs.parse(idsjobs, inputData, oStream);
        
        _runMode = rdrjobs.getRunMode();

        _updateState(rdrjobs.getState());
    }

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

    if (!_state) return;

    _assignRunMode(_runMode);
    
    mpi::bcast(idsjobs, _globRank, MPI_COMM_WORLD);
    
    _assignJobs(idsjobs);
}

void CJobsManager::runJobs(const std::string& pathToBasisSets,
                           const CInputData& inputData, COutputStream& oStream)
{
    for (size_t i = 0; i < _listOfJobs.size(); i++)
    {
        _listOfJobs[i]->set(pathToBasisSets, inputData, oStream);

        _updateState(_listOfJobs[i]->getState());

        mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

        if (!_state) return;

        _listOfJobs[i]->run(oStream, MPI_COMM_WORLD);

        _updateState(_listOfJobs[i]->getState());

        mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

        if (!_state) return;
    };
}

bool CJobsManager::getState() const
{
    return _state;
}

void CJobsManager::_assignJobs(const std::vector<int32_t>& listOfJobIds)
{
    for (size_t i = 0; i < listOfJobIds.size(); i++)
    {
        // add single point job

        if (listOfJobIds[i] == to_int(job::sp_energy))
        {
            _listOfJobs.push_back(new CSinglePointEnergy(_globRank,
                                                         _globNodes,
                                                         _runMode));
        }

        // add optimization job

        if (listOfJobIds[i] == to_int(job::opt_geometry))
        {
            _listOfJobs.push_back(new COptimizationGeometry(_globRank,
                                                            _globNodes,
                                                            _runMode));
        }

        // TODO: Add other types of jobs...
    }
}

void CJobsManager::_assignRunMode(const execmode runMode)
{
    int32_t keyval = to_int(runMode);
    
    mpi::bcast(keyval, MPI_COMM_WORLD);
    
    if (keyval == to_int(execmode::cpu)) _runMode = execmode::cpu;
    
    if (keyval == to_int(execmode::cpu_gpu)) _runMode = execmode::cpu_gpu;

    // TODO: add other execution mode...
}

void CJobsManager::_updateState(const bool state)
{
    if (_state) _state = state;
}
