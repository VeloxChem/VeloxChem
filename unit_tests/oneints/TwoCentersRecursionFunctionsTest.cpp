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

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForNuclearPotentialForSS)
{
    CRecursionTerm rta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForNuclearPotential(rta);
    
    ASSERT_EQ(1u, rvec.size());
    
    CRecursionTerm s0ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(s0ta, rvec[0]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForNuclearPotentialForSP)
{
    CRecursionTerm rta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForNuclearPotential(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k10ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                         1, 1, 0);
    
    CRecursionTerm k11ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                         1, 1, 1);
    
    ASSERT_EQ(k10ta, rvec[0]);
    
    ASSERT_EQ(k11ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForNuclearPotentialForSF)
{
    CRecursionTerm rta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 2);
    
    auto rvec = t2crecfunc::obRecursionForNuclearPotential(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm k10ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 2);
    
    CRecursionTerm k11ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 3);
    
    CRecursionTerm k20ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 2);
    
    CRecursionTerm k21ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 3);
    
    ASSERT_EQ(k10ta, rvec[0]);
    
    ASSERT_EQ(k11ta, rvec[1]);
    
    ASSERT_EQ(k20ta, rvec[2]);
    
    ASSERT_EQ(k21ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForNuclearPotentialForPP)
{
    CRecursionTerm rta({"Nuclear Potential"}, 0, false, {1, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForNuclearPotential(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm b10ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 0);
    
    CRecursionTerm b11ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 1);
    
    CRecursionTerm k10ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                         1, 1, 0);
    
    CRecursionTerm k11ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                         1, 1, 1);
    
    ASSERT_EQ(b10ta, rvec[0]);
    
    ASSERT_EQ(b11ta, rvec[1]);
    
    ASSERT_EQ(k10ta, rvec[2]);
    
    ASSERT_EQ(k11ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForNuclearPotentialForDD)
{
    CRecursionTerm rta({"Nuclear Potential"}, 0, false, {2, -1, -1, -1}, {2, -1, -1, -1},
                       1, 1, 3);
    
    auto rvec = t2crecfunc::obRecursionForNuclearPotential(rta);
    
    ASSERT_EQ(6u, rvec.size());
    
    CRecursionTerm b10ta({"Nuclear Potential"}, 0, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 3);
    
    CRecursionTerm b11ta({"Nuclear Potential"}, 0, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 4);
    
    CRecursionTerm b20ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 3);
    
    CRecursionTerm b21ta({"Nuclear Potential"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                         1, 1, 4);
    
    CRecursionTerm k10ta({"Nuclear Potential"}, 0, false, {1, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 3);
    
    CRecursionTerm k11ta({"Nuclear Potential"}, 0, false, {1, -1, -1, -1}, {1, -1, -1, -1},
                         1, 1, 4);
    
    ASSERT_EQ(b10ta, rvec[0]);
    
    ASSERT_EQ(b11ta, rvec[1]);
    
    ASSERT_EQ(b20ta, rvec[2]);
    
    ASSERT_EQ(b21ta, rvec[3]);
    
    ASSERT_EQ(k10ta, rvec[4]);
    
    ASSERT_EQ(k11ta, rvec[5]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForElectricDipoleForSP)
{
    CRecursionTerm rta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForElectricDipole(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k1ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForElectricDipoleForSF)
{
    CRecursionTerm rta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForElectricDipole(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm k1ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k3ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
    
    ASSERT_EQ(k3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForElectricDipoleForPS)
{
    CRecursionTerm rta({"Electric Dipole"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForElectricDipole(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm b1ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForElectricDipoleForDS)
{
    CRecursionTerm rta({"Electric Dipole"}, 1, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForElectricDipole(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm b1ta({"Electric Dipole"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(b3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForElectricDipoleForDF)
{
    CRecursionTerm rta({"Electric Dipole"}, 1, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForElectricDipole(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm b1ta({"Electric Dipole"}, 1, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Electric Dipole"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k1ta({"Electric Dipole"}, 1, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(k1ta, rvec[2]);
    
    ASSERT_EQ(b3ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForLinearMomentumForSP)
{
    CRecursionTerm rta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForLinearMomentum(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k1ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForLinearMomentumForSF)
{
    CRecursionTerm rta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForLinearMomentum(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm k1ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k3ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
    
    ASSERT_EQ(k3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForLinearMomentumForPS)
{
    CRecursionTerm rta({"Linear Momentum"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForLinearMomentum(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm b1ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForLinearMomentumForDS)
{
    CRecursionTerm rta({"Linear Momentum"}, 1, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForLinearMomentum(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm b1ta({"Linear Momentum"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(b3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForLinearMomentumForDF)
{
    CRecursionTerm rta({"Linear Momentum"}, 1, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForLinearMomentum(rta);
    
    ASSERT_EQ(4u, rvec.size());
    
    CRecursionTerm b1ta({"Linear Momentum"}, 1, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Linear Momentum"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k1ta({"Linear Momentum"}, 1, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(k1ta, rvec[2]);
    
    ASSERT_EQ(b3ta, rvec[3]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForAngularMomentumMomentumForSP)
{
    CRecursionTerm rta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForAngularMomentum(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm k1ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForAngularMomentumForSF)
{
    CRecursionTerm rta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForAngularMomentum(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm k1ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {1, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k3ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(k1ta, rvec[0]);
    
    ASSERT_EQ(k2ta, rvec[1]);
    
    ASSERT_EQ(k3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForAngularMomentumForPS)
{
    CRecursionTerm rta({"Angular Momentum"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForAngularMomentum(rta);
    
    ASSERT_EQ(2u, rvec.size());
    
    CRecursionTerm b1ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Overlap"}, 0, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForAngularMomentumForDS)
{
    CRecursionTerm rta({"Angular Momentum"}, 1, false, {2, -1, -1, -1}, {0, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForAngularMomentum(rta);
    
    ASSERT_EQ(3u, rvec.size());
    
    CRecursionTerm b1ta({"Angular Momentum"}, 1, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {0, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(b3ta, rvec[2]);
}

TEST_F(CTwoCentersRecursionFunctionsTest, ObRecursionForAngularMomentumForDF)
{
    CRecursionTerm rta({"Angular Momentum"}, 1, false, {2, -1, -1, -1}, {3, -1, -1, -1},
                       1, 1, 0);
    
    auto rvec = t2crecfunc::obRecursionForAngularMomentum(rta);
    
    ASSERT_EQ(5u, rvec.size());
    
    CRecursionTerm b1ta({"Angular Momentum"}, 1, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b2ta({"Angular Momentum"}, 1, false, {0, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k1ta({"Angular Momentum"}, 1, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm b3ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {3, -1, -1, -1},
                        1, 1, 0);
    
    CRecursionTerm k2ta({"Overlap"}, 0, false, {1, -1, -1, -1}, {2, -1, -1, -1},
                        1, 1, 0);
    
    ASSERT_EQ(b1ta, rvec[0]);
    
    ASSERT_EQ(b2ta, rvec[1]);
    
    ASSERT_EQ(k1ta, rvec[2]);
    
    ASSERT_EQ(b3ta, rvec[3]);
    
    ASSERT_EQ(k2ta, rvec[4]);
}
