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

#include "OMPTasksTest.hpp"

#include "omp.h"

#include "CheckFunctions.hpp"
#include "OMPTasks.hpp"

TEST_F(COMPTasksTest, GetNumberOfTasks)
{
    // set number of OMP threads

    omp_set_num_threads(4);

    COMPTasks tasks(5);

    ASSERT_EQ(20, tasks.getNumberOfTasks());
}

TEST_F(COMPTasksTest, GetTaskSizes)
{
    // set number of OMP threads

    omp_set_num_threads(2);

    COMPTasks tasks(3);

    tasks.set(26);

    vlxtest::compare({5, 5, 4, 4, 4, 4}, tasks.getTaskSizes());
}

TEST_F(COMPTasksTest, GetTaskPositions)
{
    // set number of OMP threads

    omp_set_num_threads(2);

    COMPTasks tasks(3);

    tasks.set(26);

    vlxtest::compare({0, 5, 10, 14, 18, 22}, tasks.getTaskPositions());
}
