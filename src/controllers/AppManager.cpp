//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AppManager.hpp"

#include "MpiFunc.hpp"
#include "InputStream.hpp"
#include "EnvironmentReader.hpp"
#include "JobsManager.hpp"
#include "DeviceProp.hpp"

CAppManager::CAppManager(const std::string& inputFilename,
                         const std::string& outputFilename)

    : _state(true)

    , _globRank(-1)

    , _globNodes(-1)
{
    // set up global MPI data

    _globRank  = mpi::rank(MPI_COMM_WORLD);

    _globNodes = mpi::nodes(MPI_COMM_WORLD);

    // read command line input
    
    if (_globRank == mpi::master())
    {
        _iFilename.assign(inputFilename);
            
        _oFilename.assign(outputFilename);
    }

    // update state of application manager across MPI processes

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
}

CAppManager::~CAppManager()
{

}

void
CAppManager::execute()
{
    COutputStream ostream(_oFilename);
    
    // print start header
    
    if (_globRank == mpi::master()) _printStartHeader(ostream);
    
    ostream.flush();
    
    CInputStream istream(_iFilename, ostream);
    
    updateState(istream.getState());

    // detect gpu

    if (_globRank == mpi::master())
    {
        if (_state)
        {
            #ifdef ENABLE_GPU

            gpu::get_device_prop(ostream);

            #endif
        }
    }
    
    // update state of application manager across MPI processes
    
    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
    
    if (_state)
    {
        // read input file
        
        CInputData inpdata;
        
        istream.read(inpdata, ostream);
        
        updateState(istream.getState());
        
        // update state of application manager across MPI processes
        
        mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
        
        if (_state)
        {
            CJobsManager jobsManager(_globRank, _globNodes);
            
            // set up list of jobs
            
            jobsManager.setJobs(inpdata, ostream);
            
            updateState(jobsManager.getState());
            
            // update state of application manager across MPI processes
            
            mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
            
            // read environment variables
            
            _setEnvironment(inpdata, ostream);
            
            // run jobs in list of jobs
            
            if (_state) jobsManager.runJobs(_pathToBasLib, inpdata, ostream);
            
            updateState(jobsManager.getState());
            
            // update state of application manager across MPI processes
            
            mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
        };
    }
    
    // print finish header
    
    if (_globRank == mpi::master()) _printFinishHeader(ostream);
    
    ostream.flush();
}

void
CAppManager::updateState(const bool state)
{
    if (_state) _state = state;
}

bool
CAppManager::getState() const
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

void
CAppManager::_setEnvironment(const CInputData&    inputData,
                                   COutputStream& oStream)
{
    if (_globRank == mpi::master())
    {
        CEnvironmentReader rdrenvironment;
        
        rdrenvironment.parse(inputData, oStream);
        
        updateState(rdrenvironment.getState());
        
        // assign environmental variables
        
        if (_state)
        {
            _pathToBasLib = rdrenvironment.getPathToBasisSets();
        
            // TODO: add other variables
        }
    }
    
    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
}
