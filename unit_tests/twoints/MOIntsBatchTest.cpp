//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MOIntsBatchTest.hpp"

#include "MOIntsBatch.hpp"

TEST_F(CMOIntsBatchTest, DefaultConstructor)
{
    CMOIntsBatch mbatcha;
    
    CMOIntsBatch mbatchb({}, {}, {-1, -1}, moints::oooo);
    
    ASSERT_EQ(mbatcha, mbatchb);
}

TEST_F(CMOIntsBatchTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    CMOIntsBatch mbatchb(mbatcha);
    
    ASSERT_EQ(mbatcha, mbatcha);
}

TEST_F(CMOIntsBatchTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    CMOIntsBatch mbatchb(CMOIntsBatch({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1},
                                      {0, 1}}, {2, 3}, moints::vvvv));
    
    ASSERT_EQ(mbatcha, mbatcha);
}

TEST_F(CMOIntsBatchTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    CMOIntsBatch mbatchb = mbatcha;
    
    ASSERT_EQ(mbatcha, mbatcha);
}

TEST_F(CMOIntsBatchTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    CMOIntsBatch mbatchb = CMOIntsBatch({ma, ma, mb, mb}, {{0, 1}, {0, 1},
                                        {0, 1}, {0, 1}}, {2, 3}, moints::vvvv);
    
    ASSERT_EQ(mbatcha, mbatcha);
}

TEST_F(CMOIntsBatchTest, SetBatchType)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    mbatcha.setBatchType(moints::ovov);
    
    CMOIntsBatch mbatchb({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::ovov);
    
    ASSERT_EQ(mbatcha, mbatcha);
}

TEST_F(CMOIntsBatchTest, SetExternalIndexes)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMOIntsBatch mbatcha({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {2, 3}, moints::vvvv);
    
    mbatcha.setExternalIndexes({0, 2});
    
    CMOIntsBatch mbatchb({ma, ma, mb, mb}, {{0, 1}, {0, 1}, {0, 1}, {0, 1}},
                         {0,2}, moints::vvvv);
    
    ASSERT_EQ(mbatcha, mbatcha);
}

