//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "PrimitiveFunctionalTest.hpp"

#include "PrimitiveFunctional.hpp"
#include "DummyFunctionals.hpp"

TEST_F(CPrimitiveFunctionalTest, DefaultConstructor)
{
    CPrimitiveFunctional rfa({}, xcfun::undefined, nullptr, nullptr, nullptr);
    
    CPrimitiveFunctional rfb;
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, CopyConstructor)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional rfb(rfa);
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, MoveConstructor)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional rfb(CPrimitiveFunctional({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b));
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, CopyAssignment)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional rfb = rfa;
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, MoveAssignment)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional rfb = CPrimitiveFunctional({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, GetLabel)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    ASSERT_EQ(std::string("Slater"), rfa.getLabel());
}

TEST_F(CPrimitiveFunctionalTest, GetFunctionalType)
{
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    ASSERT_EQ(xcfun::lda, rfa.getFunctionalType());
}

TEST_F(CPrimitiveFunctionalTest, Compute)
{
    CDensityGrid dgrid(2, 1, xcfun::lda, dengrid::ab);
   
    CXCGradientGrid xcgrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 2, 3), dengrid::ab, xcfun::lda);
    
    CPrimitiveFunctional rfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    rfa.compute(xcgrid, dgrid);
    
    ASSERT_EQ(xcgrid, CXCGradientGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 3), dengrid::ab, xcfun::lda));
}
