//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "XCFunctionalTest.hpp"

#include "XCFunctional.hpp"
#include "DummyFunctionals.hpp"

TEST_F(CXCFunctionalTest, DefaultConstructor)
{
    CXCFunctional rfa({}, xcfun::undefined, 0.0, {});
    
    CXCFunctional rfb;
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CXCFunctionalTest, CopyConstructor)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    CXCFunctional rfb(rfa);
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CXCFunctionalTest, MoveConstructor)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    CXCFunctional rfb(CXCFunctional({"LDA"}, xcfun::lda, 1.0, {pfa, pfb}));
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CXCFunctionalTest, CopyAssignment)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    CXCFunctional rfb = rfa;
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CXCFunctionalTest, MoveAssignment)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    CXCFunctional rfb = CXCFunctional({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    ASSERT_EQ(rfa, rfb);
}

TEST_F(CXCFunctionalTest, GetLabel)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    ASSERT_EQ(std::string("LDA"), rfa.getLabel());
}

TEST_F(CXCFunctionalTest, GetFunctionalType)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    ASSERT_EQ(xcfun::lda, rfa.getFunctionalType());
}

TEST_F(CXCFunctionalTest, GetFractionOfExactExchange)
{
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    ASSERT_NEAR(1.0, rfa.getFractionOfExactExchange(), 1.0e-13);
}

TEST_F(CXCFunctionalTest, Compute)
{
    CDensityGrid dgrid(2, 1, xcfun::lda, dengrid::ab);
    
    CXCGradientGrid xcgrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 2, 3), dengrid::ab, xcfun::lda);
    
    CPrimitiveFunctional pfa({"Slater"}, xcfun::lda,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CPrimitiveFunctional pfb({"Dirac"}, xcfun::gga,  &vlxtest::dummy_fvxc_ab, &vlxtest::dummy_fvxc_a, &vlxtest::dummy_fvxc_b);
    
    CXCFunctional rfa({"LDA"}, xcfun::lda, 1.0, {pfa, pfb});
    
    rfa.compute(xcgrid, dgrid);
    
    ASSERT_EQ(xcgrid, CXCGradientGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 3), dengrid::ab, xcfun::lda));
}
