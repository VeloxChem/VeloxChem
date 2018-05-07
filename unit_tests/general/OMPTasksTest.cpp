//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OMPTasksTest.hpp"

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "OMPTasks.hpp"
#include "CheckFunctions.hpp"

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
