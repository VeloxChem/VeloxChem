//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoCentersRecursionFunctionsTest.hpp"

#include "TwoCentersRecursionFunctions.hpp"

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForSS)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_TRUE(rvec.empty());
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForSP)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_EQ(1u, rvec.size());
    
    CRecursionTerm k1ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForSF)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k1ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForPS)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_EQ(1u, rvec.size());
    
    CRecursionTerm b1ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForDS)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm b1ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForOverlapForDF)
{
    CRecursionTerm rta({"Overlap"}, 0, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForOverlap(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm b1ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k1ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(k1ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForSS)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(1u, rvec.size());
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
     ASSERT_EQ(s0ta, rvec[0]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForSP)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k1ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(s0ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForSF)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm k1ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
    
    ASSERT_EQ(s2ta, rvec[2]);
    
    ASSERT_EQ(s0ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForPS)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm b1ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(s0ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForDS)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm b1ta({"Kinetic Energy"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(s2ta, rvec[2]);
    
    ASSERT_EQ(s0ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForKineticEnergyForDF)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForKineticEnergy(rta);
    
    ASSERT_EQ(5u, rvec.size());
    
    CRecursionTerm b1ta({"Kinetic Energy"}, 0, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Kinetic Energy"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k1ta({"Kinetic Energy"}, 0, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(s2ta, rvec[2]);
    
    ASSERT_EQ(k1ta, rvec[3]);
    
    ASSERT_EQ(s0ta, rvec[4]);
}


