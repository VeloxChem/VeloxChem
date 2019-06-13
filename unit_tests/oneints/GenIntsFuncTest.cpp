//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GenIntsFuncTest.hpp"

#include "GenIntsFunc.hpp"
#include "TwoCentersRecursionFunctions.hpp"

TEST_F(CGenIntsFuncTest, GenRecursionMap)
{
    CRecursionFunctionsList recfuncs;

    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));

    recfuncs.add(CRecursionFunction({"Kinetic Energy"}, &t2crecfunc::obRecursionForKineticEnergy));

    CRecursionTerm rta({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    auto recmap = gintsfunc::genRecursionMap(rta, recblock::cc, 1, recfuncs);

    CRecursionTerm t11({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm t10({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm t01({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm t00({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm s21({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm s11({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm s10({"Overlap"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm s01({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    CRecursionTerm s00({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);

    CRecursionMap refmap({rta, t11, t01, s01, t10, s21, s11, s10, t00, s00}, {0, 18, 27, 30, 33, 36, 54, 63, 66, 67}, recblock::cc, 1);

    ASSERT_EQ(recmap, refmap);
}

TEST_F(CGenIntsFuncTest, GenIntegral)
{
    CRecursionTerm rta({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    ASSERT_EQ(rta, gintsfunc::genIntegral({"Kinetic Energy"}, 2, 1, 0));

    CRecursionTerm rtb({"Electric Field Gradient"}, 2, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);

    ASSERT_EQ(rtb, gintsfunc::genIntegral({"Electric Field Gradient"}, 2, 1, 4));

    CRecursionTerm rtc({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);

    ASSERT_EQ(rtc, gintsfunc::genIntegral({"Linear Momentum"}, 2, 1, 4));

    CRecursionTerm rtd({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);

    ASSERT_EQ(rtd, gintsfunc::genIntegral({"Angular Momentum"}, 2, 1, 4));

    CRecursionTerm rte({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);

    ASSERT_EQ(rte, gintsfunc::genIntegral({"Electric Dipole"}, 2, 1, 4));

    CRecursionTerm rtf({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);

    ASSERT_EQ(rtf, gintsfunc::genIntegral({"Electric Field"}, 2, 1, 4));

    CRecursionTerm rtk({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    ASSERT_EQ(rtk, gintsfunc::genIntegral({"Overlap"}, 2, 1, 0));

    CRecursionTerm rtl({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);

    ASSERT_EQ(rtl, gintsfunc::genIntegral({"Nuclear Potential"}, 2, 1, 0));
}
