//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ExcitationVectorTest.hpp"

#include "ExcitationVector.hpp"

TEST_F(CExcitationVectorTest, DefaultConstructor)
{
    CExcitationVector ma;
    
    CExcitationVector mb(szblock::aa, std::vector<int32_t>{},
                         std::vector<int32_t>{},  std::vector<double>{});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, CopyConstructor)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2}, {1.0, -1.0, -3.0});
    
    CExcitationVector mb(ma);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, MoveConstructor)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2}, {1.0, -1.0, -3.0});
    
    CExcitationVector mb(CExcitationVector(szblock::aa, {0, 0, 1}, {2, 3, 2},
                                           {1.0, -1.0, -3.0}));
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, CopyAssignment)
{
    CExcitationVector ma(szblock::ab, {0, 0, 1}, {2, 3, 2}, {1.0, -1.0, -3.0});
    
    CExcitationVector mb = ma;
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, MoveAssignment)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2}, {1.0, -1.0, -3.0});
    
    CExcitationVector mb = CExcitationVector(szblock::bb, {0, 0, 1}, {2, 3, 2},
                                             {1.0, -1.0, -3.0});
    
    ASSERT_EQ(ma, mb);
}
