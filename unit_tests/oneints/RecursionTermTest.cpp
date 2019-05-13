//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionTermTest.hpp"

#include "RecursionTerm.hpp"

TEST_F(CRecursionTermTest, DefaultConstructor)
{
    CRecursionTerm rta;
    
    CRecursionTerm rtb(std::string(), -1, false, {-1, -1, -1, -1}, {-1, -1, -1, -1},
                       -1, -1, -1);
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, CopyConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb(rta);
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, MoveConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb(CRecursionTerm({"Overlap"}, 0, true, {2, 3, 4, 5},
                                      {1, 7, 2, 3}, 1, 2, 5));
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, CopyAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb = rta;
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, MoveAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb = CRecursionTerm({"Overlap"}, 0, true, {2, 3, 4, 5},
                                        {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, IsValid)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_TRUE(rta.isValid());
    
    CRecursionTerm rtb({"Overlap"}, 0, false, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_TRUE(rtb.isValid());
    
    CRecursionTerm rtc({""}, 0, false, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_FALSE(rtc.isValid());
    
    CRecursionTerm rtd({"Overlap"}, -1, false, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_FALSE(rtd.isValid());
    
    CRecursionTerm rte({"Overlap"}, 0, false, {2, 3, 4, 5}, {1, 7, 2, 3}, -1, 2, 5);
    
    ASSERT_FALSE(rte.isValid());
    
    CRecursionTerm rtf({"Overlap"}, 0, false, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, -2, 5);
    
    ASSERT_FALSE(rtf.isValid());
    
    CRecursionTerm rtg({"Overlap"}, 0, false, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, -1);
    
    ASSERT_FALSE(rtg.isValid());
    
    CRecursionTerm rth({"Overlap"}, 0, false, {2, -1, -2, -3}, {1, 7, 2, 3}, 1, 2, 0);
    
    ASSERT_TRUE(rth.isValid());
    
    CRecursionTerm rti({"Overlap"}, 0, false, {2, -1, -2, -3}, {1, 7, 2, 3}, 2, 2, 0);
    
    ASSERT_FALSE(rti.isValid());
    
    CRecursionTerm rtk({"Overlap"}, 1, false, {2, -1, 3, 4}, {1, 7, -2, -3}, 1, 2, 5);
    
    ASSERT_TRUE(rtk.isValid());
    
    CRecursionTerm rtl({"Overlap"}, 1, false, {2, -1, 3, 4}, {-1, 7, -2, -3}, 1, 2, 5);
    
    ASSERT_FALSE(rtl.isValid());
    
    CRecursionTerm rtm({"Overlap"}, 1, false, {2, -1, 3, 4}, {1, -7, -2, -3}, 1, 2, 5);
    
    ASSERT_FALSE(rtm.isValid());
}

