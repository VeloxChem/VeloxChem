//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "SparseVectorTest.hpp"

#include "SparseVector.hpp"
#include "CheckFunctions.hpp"

TEST_F(CSparseVectorTest, DefaultConstructor)
{
    CSparseVector veca;
    
    CSparseVector vecb(std::vector<double>{}, std::vector<int32_t>{}, 0, 1.0e-13);
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);
    
    CSparseVector vecb(veca);
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);
    
    CSparseVector vecb(CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4},
                                     10, 1.0e-13));
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);
    
    CSparseVector vecb = veca;
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10,  1.0e-13);
    
    CSparseVector vecb = CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4},
                                        10, 1.0e-13);
    
    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, ValuesConstant)
{
    const CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0}, veca.values());
}

TEST_F(CSparseVectorTest, Values)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0}, veca.values());
    
    auto vdat = veca.values();
    
    vdat[2] = 9.0;
    
    CSparseVector vecb({1.0, -1.0, 9.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);
    
    EXPECT_EQ(veca, vecb); 
}

TEST_F(CSparseVectorTest, IndexesConstant)
{
    const CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);
    
    vlxtest::compare({0, 1, 2, 4}, veca.indexes());
}

TEST_F(CSparseVectorTest, GetThreshold)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);
    
    ASSERT_NEAR(veca.getThreshold(), 0.28, 1.0e-13);
}

TEST_F(CSparseVectorTest, GetSparsity)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);
    
    ASSERT_NEAR(veca.getSparsity(), 0.40, 1.0e-13);
}
