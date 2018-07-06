//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SparseMatrixTest.hpp"

#include "SparseMatrix.hpp"

TEST_F(CSparseMatrixTest, DefaultConstructor)
{
    CSparseMatrix ma;
    
    CSparseMatrix mb(std::vector<double>{}, std::vector<int32_t>{},
                     std::vector<int32_t>{}, 0, 0, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, CopyConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                      {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb(ma);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, MoveConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb(CSparseMatrix({ 1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0,
                                    -4.0, 2.0, 7.0, 8.0, -5.0}, {0, 0, 0, 1, 1,
                                     2, 2, 2, 3, 3, 3, 4, 4}, {0, 1, 2, 0, 1, 2,
                                     3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13));
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, CopyAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb = ma;
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, MoveAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb = CSparseMatrix({ 1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0,
                                      -4.0, 2.0, 7.0, 8.0, -5.0}, {0, 0, 0, 1, 1,
                                       2, 2, 2, 3, 3, 3, 4, 4}, {0, 1, 2, 0, 1, 2,
                                       3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
}


