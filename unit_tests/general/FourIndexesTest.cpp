//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

TEST_F(CFourIndexesTest, Shift)
{
    CFourIndexes ida(2, 7, 0, 3);

    ida.shift(2, 0);

    ASSERT_EQ(ida, CFourIndexes(4, 7, 0, 3));

    ida.shift(-2, 1);

    ASSERT_EQ(ida, CFourIndexes(4, 5, 0, 3));

    ida.shift(3, 2);

    ASSERT_EQ(ida, CFourIndexes(4, 5, 3, 3));

    ida.shift(2, 3);

    ASSERT_EQ(ida, CFourIndexes(4, 5, 3, 5));

    ida.shift(2, -1);

    ASSERT_EQ(ida, CFourIndexes(4, 5, 3, 5));

    ida.shift(2, 4);

    ASSERT_EQ(ida, CFourIndexes(4, 5, 3, 5));
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

TEST_F(CFourIndexesTest, Value)
{
    CFourIndexes ida(2, 7, 0, 3);

    ASSERT_EQ(2, ida.value(0));

    ASSERT_EQ(7, ida.value(1));

    ASSERT_EQ(0, ida.value(2));

    ASSERT_EQ(3, ida.value(3));

    ASSERT_EQ(-1, ida.value(-2));

    ASSERT_EQ(-1, ida.value(4));
}
