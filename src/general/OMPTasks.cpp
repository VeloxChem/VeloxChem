//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "OMPTasks.hpp"

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "MathFunc.hpp"
#include "MpiFunc.hpp"

COMPTasks::COMPTasks(const int32_t nTasksPerThread)

    : _nOMPThreads(omp_get_max_threads())

    , _nTasksPerThread(nTasksPerThread)
{
    if (_nOMPThreads < 1) _nOMPThreads = 1;

    if (_nTasksPerThread < 1) _nTasksPerThread = 1;

    // set up tasks data

    _taskSizes = CMemBlock<int32_t>(_nTasksPerThread * _nOMPThreads);

    _taskPositions = CMemBlock<int32_t>(_nTasksPerThread * _nOMPThreads);
}

COMPTasks::~COMPTasks()
{
}

void
COMPTasks::set(const int32_t nElements)
{
    auto ndim = _nTasksPerThread * _nOMPThreads;

    // determine tasks sizes

    for (int32_t i = 0; i < ndim; i++)
    {
        _taskSizes.at(i) = mpi::batch_size(nElements, i, ndim);
    }

    // determining tasks start positions

    mathfunc::indexes(_taskPositions.data(), _taskSizes.data(), ndim);
}

int32_t
COMPTasks::getNumberOfTasks() const
{
    return _nTasksPerThread * _nOMPThreads;
}

const int32_t*
COMPTasks::getTaskSizes() const
{
    return _taskSizes.data();
}

const int32_t*
COMPTasks::getTaskPositions() const
{
    return _taskPositions.data();
}
