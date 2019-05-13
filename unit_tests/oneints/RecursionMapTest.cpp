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
    
    CRecursionMap rmb({}, {});
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4});
    
    CRecursionMap rmb(rma);
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4});
    
    CRecursionMap rmb(CRecursionMap({rta, rtb}, {0, 4}));
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4});
    
    CRecursionMap rmb = rma;
    
    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3},
                       1, 2, 5);
    
    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3},
                       2, 3, 0);
    
    CRecursionMap rma({rta, rtb}, {0, 4});
    
    CRecursionMap rmb = CRecursionMap({rta, rtb}, {0, 4});
    
    ASSERT_EQ(rma, rmb);
}
