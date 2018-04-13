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

void
CAppManager::_printStartHeader(COutputStream& oStream)
{
    oStream << fmt::title << fmt::tsep << fmt::end;
    
    oStream << "VELOX CHEM MP (V.0.0 2018)" << fmt::end << fmt::end;
    
    oStream << "AN ELECTRONIC STRUCTURE CODE FOR NANOSCALE" << fmt::end;
    
    oStream << fmt::end << fmt::left;
    
    oStream << " Copyright (c) 2018 Velox Chem MP developers.";
    
    oStream << " All rights reserved.";
    
    oStream << fmt::end << fmt::tsep << fmt::center;
    
    oStream << " VeloxChem MP execution started";
    
    if (_globNodes > 1)
    {
        oStream <<  " on " <<std::to_string(_globNodes) << " compute nodes";
    }
    
    oStream << " at " << _sysClock.getStartDate() << ".";
    
    oStream << fmt::end << fmt::tsep << fmt::blank;
}

void
CAppManager::_printFinishHeader(COutputStream& oStream)
{
    oStream << fmt::blank << fmt::title << fmt::tsep << fmt::center;
    
    oStream << "VeloxChem MP execution ";
    
    if (_state)
    {
        oStream << "completed at ";
    }
    else
    {
        oStream << "abnormally terminated at ";
    }
    
    oStream << _sysClock.getCurrentDate() << "." << fmt::end << fmt::tsep;
    
    oStream << "Total execution time is " << _sysClock.getElapsedTime();
    
    oStream << fmt::end << fmt::tsep;
}
