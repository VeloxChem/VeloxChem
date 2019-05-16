//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionMapTest.hpp"

#include "RecursionMap.hpp"

TEST_F(CRecursionMapTest, DefaultConstructor)
{
    CRecursionMap rma;
    
    CRecursionMap rmb({}, {}, recblock::cc);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, AlternativeConstructor)
{
    CRecursionMap rma(recblock::ss);
    
    CRecursionMap rmb({}, {}, recblock::ss);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc);
    
    CRecursionMap rmb(rma);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc);
    
    CRecursionMap rmb(CRecursionMap({rta, rtb}, {0, 4}, recblock::cc));
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc);
    
    CRecursionMap rmb = rma;
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc);
    
    CRecursionMap rmb = CRecursionMap({rta, rtb}, {0, 4}, recblock::cc);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, Add)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3},
                       1, 2, 0);
    
    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, -1, 7, 3},
                       1, 2, 0);
    
    CRecursionMap rma;
    
    rma.add(rta);
    
    rma.add(rtb);
    
    rma.add(rta);
    
    rma.add(rta);
    
    rma.add(rtb);
    
    rma.add(rtc);
    
    CRecursionMap rmb({rta, rtb}, {0, 18}, recblock::cc);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, Append)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3},
                       1, 2, 0);
    
    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, -1, 7, 3},
                       1, 2, 0);
    
    CRecursionTerm rtd({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3},
                       1, 2, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 18}, recblock::cc);
    
    CRecursionMap rmb({rtc, rtd, rta}, {0, 0, 0}, recblock::cc);
    
    rma.append(rmb);
    
    CRecursionMap rmc({rta, rtb, rtd}, {0, 18, 828}, recblock::cc);
    
    ASSERT_EQ(rma, rmc);
}

TEST_F(CRecursionMapTest, Find)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3},
                       1, 2, 0);
    
    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, 3, 7, 3},
                       1, 2, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 18}, recblock::cc);
    
    ASSERT_TRUE(rma.find(rta));
    
    ASSERT_TRUE(rma.find(rtb));
    
    ASSERT_FALSE(rma.find(rtc));
}

TEST_F(CRecursionMapTest, GetNumberOfComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3},
                       1, 2, 0);
    
    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3},
                       1, 2, 0);
    
    CRecursionMap rma({rta, rtb, rtc}, {0, 18, 828}, recblock::cc);
    
    ASSERT_EQ(1152, rma.getNumberOfComponents());
}

