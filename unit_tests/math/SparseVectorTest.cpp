//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SparseVectorTest.hpp"

#include "SparseVector.hpp"

TEST_F(CSparseVectorTest, DefaultConstructor)
{
    CSparseVector veca;
    
    CSparseVector vecb(std::vector<double>{}, std::vector<int32_t>{}, 1.0e-13);
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 1.0e-13);
    
    CSparseVector vecb(veca);
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 1.0e-13);
    
    CSparseVector vecb(CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4},
                                     1.0e-13));
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 1.0e-13);
    
    CSparseVector vecb = veca;
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 1.0e-13);
    
    CSparseVector vecb = CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4},
                                        1.0e-13);
    
    ASSERT_EQ(veca, vecb);
}
