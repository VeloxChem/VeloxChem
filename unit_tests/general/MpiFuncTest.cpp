//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MpiFuncTest.hpp"

#include "MpiFunc.hpp"

TEST_F(CMpiFuncTest, Master)
{
    ASSERT_EQ(0, mpi::master());
}

TEST_F(CMpiFuncTest, Rank)
{
    ASSERT_EQ(0, mpi::rank(MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, Nodes)
{
    ASSERT_EQ(1, mpi::nodes(MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, Compare)
{
    ASSERT_TRUE(mpi::compare(MPI_COMM_WORLD, MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, BcastInteger)
{
    int32_t mvalue = 3;
    
    mpi::bcast(mvalue, MPI_COMM_WORLD);
    
    ASSERT_EQ(3, mvalue);
}

TEST_F(CMpiFuncTest, BcastDouble)
{
    double mvalue = 3.4;
    
    mpi::bcast(mvalue, MPI_COMM_WORLD);
    
    ASSERT_NEAR(3.4, mvalue, 1.0e-13);
}

TEST_F(CMpiFuncTest, BcastBoolean)
{
    bool mvalue = false;
    
    mpi::bcast(mvalue, mpi::master(),  MPI_COMM_WORLD);
    
    ASSERT_FALSE(mvalue);
}

TEST_F(CMpiFuncTest, BcastIntegerVector)
{
    std::vector<int32_t> mvector({2, 3, 4, -7, 9, 18});
    
    mpi::bcast(mvector, mpi::master(),  MPI_COMM_WORLD);
    
    ASSERT_EQ(mvector[0], 2);
    
    ASSERT_EQ(mvector[1], 3);
    
    ASSERT_EQ(mvector[2], 4);
    
    ASSERT_EQ(mvector[3], -7);
    
    ASSERT_EQ(mvector[4], 9);
    
    ASSERT_EQ(mvector[5], 18);
}



