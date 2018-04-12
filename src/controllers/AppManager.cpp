//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AppManager.hpp"

#include "MpiFunc.hpp"

CAppManager::CAppManager(int argc, char** argv)

    : _state(true)

    , _globRank(-1)

    , _globNodes(-1)
{
    // set up global MPI data

    _globRank  = mpi::rank(MPI_COMM_WORLD);

    _globNodes = mpi::nodes(MPI_COMM_WORLD);

    

    // update state of application manager across MPI processes

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
}

CAppManager::~CAppManager()
{

}

void CAppManager::execute()
{
   
}

void CAppManager::updateState(const bool state)
{
    if (_state) _state = state;
}

bool CAppManager::getState() const
{
    return _state;
}
