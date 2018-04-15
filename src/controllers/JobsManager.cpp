//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

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
{

}

CJobsManager::~CJobsManager()
{
    for (size_t i = 0; i < _listOfJobs.size(); i++)
    {
        delete _listOfJobs[i];
    }
}

void
CJobsManager::setJobs(const CInputData&    inputData,
                            COutputStream& oStream)
{
    std::vector<int32_t> idsjobs;

    if (_globRank == mpi::master())
    {
        CJobsReader rdrjobs;

        rdrjobs.parse(idsjobs, inputData, oStream);

        _updateState(rdrjobs.getState());
    }

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

    if (!_state) return;

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

        _listOfJobs[i]->run(oStream);

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
                                                         _globNodes));
        }

        // add optimization job

        if (listOfJobIds[i] == to_int(job::opt_geometry))
        {
            _listOfJobs.push_back(new COptimizationGeometry(_globRank,
                                                            _globNodes));
        }

        // TODO: Add other types of jobs...
    }
}

void CJobsManager::_updateState(const bool state)
{
    if (_state) _state = state;
}
