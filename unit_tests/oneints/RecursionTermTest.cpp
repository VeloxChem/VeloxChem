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

TEST_F(CRecursionTermTest, SetLabel)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    rta.setLabel("Kinetic Energy");
    
    CRecursionTerm rtb({"Kinetic Energy"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta, rtb);
}

TEST_F(CRecursionTermTest, SetOrder)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 2);
    
    rta.setOrder(5);
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta, rtb);
}


TEST_F(CRecursionTermTest, BraShift)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {4, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta.braShift(2, 0), rtb);
    
    CRecursionTerm rtc;
    
    ASSERT_EQ(rta.braShift(2, -1), rtc);
    
    ASSERT_EQ(rta.braShift(2, 1), rtc);
}

TEST_F(CRecursionTermTest, KetShift)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {3, 7, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta.ketShift(2, 0), rtb);
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 4, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta.ketShift(-3, 1), rtc);
    
    CRecursionTerm rtd;
    
    ASSERT_EQ(rta.braShift(2, -1), rtd);
    
    ASSERT_EQ(rta.braShift(2, 2), rtd);
}

TEST_F(CRecursionTermTest, OrderShift)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 4);
    
    ASSERT_EQ(rta.orderShift(-1), rtb);
}

TEST_F(CRecursionTermTest, OperatorShift)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Overlap"}, 1, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(rta.operatorShift(1), rtb);
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

TEST_F(CRecursionTermTest, IsBraOfZeroOrder)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_FALSE(rta.isBraOfZeroOrder());
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {0, 0, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);
    
    ASSERT_TRUE(rtb.isBraOfZeroOrder());
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {0, 0, 4, 5}, {1, 7, 2, 3}, 2, 2, 5);
    
    ASSERT_TRUE(rtc.isBraOfZeroOrder());
    
    CRecursionTerm rtd({"Overlap"}, 0, true, {0, 0, 4, 5}, {1, 7, 2, 3}, 3, 2, 5);
    
    ASSERT_FALSE(rtd.isBraOfZeroOrder());
    
    CRecursionTerm rte({"Overlap"}, 0, true, {0, 0, 4, 5}, {1, 7, 2, 3}, -1, 2, 5);
    
    ASSERT_FALSE(rte.isBraOfZeroOrder());
}

TEST_F(CRecursionTermTest, GetLabel)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 2);
    
    ASSERT_EQ(std::string("Overlap"), rta.getLabel());
}

TEST_F(CRecursionTermTest, GetOrder)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 2);
    
    CRecursionTerm rtb({"Overlap"}, 3, false, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 5);
    
    ASSERT_EQ(2, rta.getOrder());
    
    ASSERT_EQ(5, rtb.getOrder());
}

TEST_F(CRecursionTermTest, GetNumberOfOperatorComponents)
{
    CRecursionTerm rta({"Overlap"}, 2, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(6, rta.getNumberOfOperatorComponents());
    
    CRecursionTerm rtb({"Overlap"}, 2, false, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 2, 5);
    
    ASSERT_EQ(9, rtb.getNumberOfOperatorComponents());
}

TEST_F(CRecursionTermTest, braCartesianComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(6, rta.braCartesianComponents());
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 2, 5);
    
    ASSERT_EQ(60, rtb.braCartesianComponents());
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       3, 2, 5);
    
    ASSERT_EQ(900, rtc.braCartesianComponents());
}

TEST_F(CRecursionTermTest, ketCartesianComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 5);
    
    ASSERT_EQ(3, rta.ketCartesianComponents());
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 2, 5);
    
    ASSERT_EQ(30, rtb.ketCartesianComponents());
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 3, 5);
    
    ASSERT_EQ(180, rtc.ketCartesianComponents());
}

TEST_F(CRecursionTermTest, braSphericalComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       1, 2, 5);
    
    ASSERT_EQ(5, rta.braSphericalComponents());
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 2, 5);
    
    ASSERT_EQ(35, rtb.braSphericalComponents());
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       3, 2, 5);
    
    ASSERT_EQ(315, rtc.braSphericalComponents());
}

TEST_F(CRecursionTermTest, ketSphericalComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 5);
    
    ASSERT_EQ(3, rta.ketSphericalComponents());
    
    CRecursionTerm rtb({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 2, 5);
    
    ASSERT_EQ(21, rtb.ketSphericalComponents());
    
    CRecursionTerm rtc({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 3, 5);
    
    ASSERT_EQ(105, rtc.ketSphericalComponents());
}

TEST_F(CRecursionTermTest, GetNumberOfComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 5);
    
    ASSERT_EQ(180, rta.getNumberOfComponents(recblock::cc));
    
    ASSERT_EQ(105, rta.getNumberOfComponents(recblock::sc));
    
    ASSERT_EQ(180, rta.getNumberOfComponents(recblock::cs));
    
    ASSERT_EQ(105, rta.getNumberOfComponents(recblock::ss));
}

TEST_F(CRecursionTermTest, IsIntegral)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 3, 2, 3},
                       2, 1, 5);
    
    ASSERT_TRUE(rta.isIntegral({"Overlap"}, {2, 3, 4, 5}, {1, 3, 2, 3}, 2, 1));
    
    ASSERT_FALSE(rta.isIntegral({"Overlap"}, {2, 3, 4, 5}, {1, 3, 2, 3}, 2, 2));
    
    ASSERT_FALSE(rta.isIntegral({"Overlap"}, {2, 3, 4, 5}, {1, 3, 2, 3}, 1, 2));
    
    ASSERT_FALSE(rta.isIntegral({"Overlap"}, {1, 3, 4, 5}, {1, 3, 2, 3}, 2, 1));
    
    ASSERT_FALSE(rta.isIntegral({"Overlap"}, {2, 1, 4, 5}, {1, 3, 2, 3}, 2, 1));
    
    ASSERT_FALSE(rta.isIntegral({"Overlap"}, {2, 3, 4, 5}, {3, 3, 2, 3}, 2, 1));
    
    ASSERT_FALSE(rta.isIntegral({"Kinetic Energy"}, {2, 3, 4, 5}, {1, 3, 2, 3}, 2, 1));
}


