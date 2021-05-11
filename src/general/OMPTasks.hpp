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

#ifndef OMPTasks_hpp
#define OMPTasks_hpp

#include <cstdint>

#include "MemBlock.hpp"

/**
 Class COMPTasks stores information about OMP tasks distribution for linear
 jobs lists.

 @author Z. Rinkevicius
 */
class COMPTasks
{
    /**
     The number of OMP threads.
     */
    int32_t _nOMPThreads;

    /**
     The number of OMP tasks assigned to single OMP thread.
     */
    int32_t _nTasksPerThread;

    /**
     The vector of task sizes in linear jobs list.
     */
    CMemBlock<int32_t> _taskSizes;

    /**
     The vector of task positions in linear jobs list.
     */
    CMemBlock<int32_t> _taskPositions;

   public:
    /**
     Creates a OMP tasks object by copying other atom basis object.

     @param nTasksPerThread the number of OMP tasks per a single OMP thread.
     */
    COMPTasks(const int32_t nTasksPerThread);

    /**
     Destroys a OMP tasks object.
     */
    ~COMPTasks();

    /**
     Sets the sizes and start positions of all tasks in OMP tasks object for
     specific linear jobs list.

     @param nElements the number of elements in linear jobs list.
     */
    void set(const int32_t nElements);

    /**
     Gets constant pointer to task sizes vector.

     @return the pointer to task sizes vector.
     */
    const int32_t* getTaskSizes() const;

    /**
     Gets constant pointer to task start positions vector.

     @return the pointer to task start positions vector.
     */
    const int32_t* getTaskPositions() const;

    /**
     Gets total number of tasks in OMP tasks object.

     @return the total number of tasks.
     */
    int32_t getNumberOfTasks() const;
};

#endif /* OMPTasks_hpp */
