//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GenIntsFuncTest.hpp"

#include "TwoCentersRecursionFunctions.hpp"
#include "GenIntsFunc.hpp"

TEST_F(CGenIntsFuncTest, GenRecursionMap)
{
    CRecursionFunctionsList recfuncs;
    
    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));
    
    recfuncs.add(CRecursionFunction({"Kinetic Energy"}, &t2crecfunc::obRecursionForKineticEnergy));
    
    CRecursionTerm rta({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    auto recmap = gintsfunc::genRecursionMap(rta, recblock::cc, recfuncs);
    
    CRecursionTerm t11({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm t10({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm t01({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm t00({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm s21({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm s11({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm s10({"Overlap"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm s01({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    CRecursionTerm s00({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0);
    
    CRecursionMap refmap({rta, t11, t01, s01, t10, s21, s11, s10, t00, s00},
                         {  0,  18,  27,  30,  33,  36,  54,  63,  66,  67},
                         recblock::cc);
    
    ASSERT_EQ(recmap, refmap);
}

TEST_F(CGenIntsFuncTest, GenIntegral)
{
    CRecursionFunctionsList recfuncs;

    CRecursionTerm rta({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0);
    
    ASSERT_EQ(rta, gintsfunc::genIntegral({"Kinetic Energy"}, 2, 1, 0));
    
    CRecursionTerm rtb({"Electric Field Gradient"}, 2, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 4);
    
    ASSERT_EQ(rtb, gintsfunc::genIntegral({"Electric Field Gradient"}, 2, 1, 4));
}

