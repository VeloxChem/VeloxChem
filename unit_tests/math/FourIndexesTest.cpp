//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "FourIndexesTest.hpp"

#include "FourIndexes.hpp"

TEST_F(CFourIndexesTest, DefaultConstructor)
{
    CFourIndexes ida(-1, -1, -1, -1);
    
    CFourIndexes idb;
    
    ASSERT_EQ(ida, idb);
}

TEST_F(CFourIndexesTest, CopyConstructor)
{
    CFourIndexes ida(2, 7, 3, 4);
    
    CFourIndexes idb(ida);
    
    ASSERT_EQ(ida, idb);
}

TEST_F(CFourIndexesTest, MoveConstructor)
{
    CFourIndexes ida(2, 7, 3, 4);
    
    CFourIndexes idb(CFourIndexes(2, 7, 3, 4));
    
    ASSERT_EQ(ida, idb);
}

TEST_F(CFourIndexesTest, CopyAssignment)
{
    CFourIndexes ida(2, 7, 3, 8);
    
    CFourIndexes idb = ida;
    
    ASSERT_EQ(ida, idb);
}

TEST_F(CFourIndexesTest, MoveAssignment)
{
    CFourIndexes ida(2, 7, 3, 1);
    
    CFourIndexes idb = CFourIndexes(2, 7, 3, 1);
    
    ASSERT_EQ(ida, idb);
}


TEST_F(CFourIndexesTest, First)
{
    CFourIndexes ida(2, 7, 1, 3);
    
    ASSERT_EQ(2, ida.first());
    
    CFourIndexes idb(3, 2, 4, 7);
    
    ASSERT_EQ(3, idb.first());
}

TEST_F(CFourIndexesTest, Second)
{
    CFourIndexes ida(2, 7, 1, 3);
    
    ASSERT_EQ(7, ida.second());
    
    CFourIndexes idb(3, 2, 4, 7);
    
    ASSERT_EQ(2, idb.second());
}

TEST_F(CFourIndexesTest, Third)
{
    CFourIndexes ida(2, 7, 1, 3);
    
    ASSERT_EQ(1, ida.third());
    
    CFourIndexes idb(3, 2, 4, 7);
    
    ASSERT_EQ(4, idb.third());
}

TEST_F(CFourIndexesTest, Fourth)
{
    CFourIndexes ida(2, 7, 1, 3);
    
    ASSERT_EQ(3, ida.fourth());
    
    CFourIndexes idb(3, 2, 4, 7);
    
    ASSERT_EQ(7, idb.fourth());
}

TEST_F(CFourIndexesTest, IsValidQuadruple)
{
    CFourIndexes ida(2, 7, 0, 3);
    
    ASSERT_TRUE(ida.isValidQuadruple());
    
    CFourIndexes idb(3, -3, 1, 0);
    
    ASSERT_FALSE(idb.isValidQuadruple());
    
    CFourIndexes idc(-3, 0, 1, 2);
    
    ASSERT_FALSE(idc.isValidQuadruple());
    
    CFourIndexes idd(3, 0, -1, 2);
    
    ASSERT_FALSE(idd.isValidQuadruple());
    
    CFourIndexes ide(3, 0, 1, -2);
    
    ASSERT_FALSE(ide.isValidQuadruple());
    
    CFourIndexes idf(0, 0, 0, 0);
    
    ASSERT_TRUE(idf.isValidQuadruple());
}
