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
