//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "OMPTasks.hpp"

#include "omp.h"

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
