//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "XCGradientGridTest.hpp"

#include "XCGradientGrid.hpp"
#include "CheckFunctions.hpp"

TEST_F(CXCGradientGridTest, DefaultConstructor)
{
    CXCGradientGrid dgrida;
    
    CMemBlock2D<double> mblock;
    
    CXCGradientGrid dgridb(0, dengrid::undefined, xcfun::undefined);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCGradientGridTest, ConstructorWithFunctionalType)
{
    CXCGradientGrid dgrida(5, dengrid::ab, xcfun::lda);
    
    ASSERT_EQ(dgrida, CXCGradientGrid(CMemBlock2D<double>(5, 3), dengrid::ab, xcfun::lda));
    
    CXCGradientGrid dgridb(5, dengrid::lima, xcfun::lda);
    
    ASSERT_EQ(dgridb, CXCGradientGrid(CMemBlock2D<double>(5, 3), dengrid::lima, xcfun::lda));
    
    CXCGradientGrid dgridc(5, dengrid::limb, xcfun::lda);
    
    ASSERT_EQ(dgridc, CXCGradientGrid(CMemBlock2D<double>(5, 3), dengrid::limb, xcfun::lda));
    
    CXCGradientGrid dgride(5, dengrid::ab, xcfun::gga);
    
    ASSERT_EQ(dgride, CXCGradientGrid(CMemBlock2D<double>(5, 6), dengrid::ab, xcfun::gga));
    
    CXCGradientGrid dgridf(5, dengrid::lima, xcfun::gga);
    
    ASSERT_EQ(dgridf, CXCGradientGrid(CMemBlock2D<double>(5, 6), dengrid::lima, xcfun::gga));
    
    CXCGradientGrid dgridg(5, dengrid::limb, xcfun::gga);
    
    ASSERT_EQ(dgridg, CXCGradientGrid(CMemBlock2D<double>(5, 6), dengrid::limb, xcfun::gga));
    
    CXCGradientGrid dgridk(5, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(dgridk, CXCGradientGrid(CMemBlock2D<double>(5, 8), dengrid::ab, xcfun::mgga));
    
    CXCGradientGrid dgridl(5, dengrid::lima, xcfun::mgga);
    
    ASSERT_EQ(dgridl, CXCGradientGrid(CMemBlock2D<double>(5, 4), dengrid::lima, xcfun::mgga));
    
    CXCGradientGrid dgridm(5, dengrid::limb, xcfun::mgga);
    
    ASSERT_EQ(dgridm, CXCGradientGrid(CMemBlock2D<double>(5, 4), dengrid::limb, xcfun::mgga));
}

TEST_F(CXCGradientGridTest, CopyConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::lima, xcfun::mgga);
    
    CXCGradientGrid dgridb(dgrida);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCGradientGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    CXCGradientGrid dgridb((CXCGradientGrid(mblock, dengrid::ab, xcfun::mgga)));
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCGradientGridTest, CopyAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::limb, xcfun::mgga);
    
    CXCGradientGrid dgridb = dgrida;
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCGradientGridTest, MoveAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    CXCGradientGrid dgridb = CXCGradientGrid(mblock, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCGradientGridTest, Zero)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    dgrida.zero();
    
    ASSERT_EQ(dgrida, CXCGradientGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 4), dengrid::ab, xcfun::mgga));
}

TEST_F(CXCGradientGridTest, GetNumberOfGridPoints)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCGradientGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(2, dgrida.getNumberOfGridPoints());
}

TEST_F(CXCGradientGridTest, XCFunctionalValuesConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 3.0, 5.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 4.0, 7.0}, 2, 2);
    
    const CXCGradientGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    const CXCGradientGrid dgridla(mblockb, dengrid::lima, xcfun::lda);
    
    const CXCGradientGrid dgridlb(mblocka, dengrid::limb, xcfun::lda);
    
    vlxtest::compare({1.0, 2.0}, dgridab.xcFunctionalValues());
    
    vlxtest::compare({1.0, 2.0}, dgridla.xcFunctionalValues());
    
    vlxtest::compare({3.0, 6.0}, dgridlb.xcFunctionalValues());
}

TEST_F(CXCGradientGridTest, XCFunctionalValues)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 3.0, 5.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 4.0, 7.0}, 2, 2);
    
    CXCGradientGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    CXCGradientGrid dgridla(mblockb, dengrid::lima, xcfun::lda);
    
    CXCGradientGrid dgridlb(mblocka, dengrid::limb, xcfun::lda);
    
    auto xcvab = dgridab.xcFunctionalValues();
    
    auto xcva = dgridla.xcFunctionalValues();
    
    auto xcvb = dgridlb.xcFunctionalValues();
    
    vlxtest::compare({1.0, 2.0}, xcvab);
    
    vlxtest::compare({1.0, 2.0}, xcva);
    
    vlxtest::compare({3.0, 6.0}, xcvb);
    
    xcvab[1] = 3.0;
    
    ASSERT_EQ(dgridab, CXCGradientGrid(CMemBlock2D<double>({1.0, 3.0, 3.0, 6.0, 1.0, 7.0}, 2, 3), dengrid::ab, xcfun::lda));
    
    xcva[0] = 0.0;
    
    ASSERT_EQ(dgridla, CXCGradientGrid(CMemBlock2D<double>({0.0, 2.0, 4.0, 7.0}, 2, 2), dengrid::lima, xcfun::lda));
    
    xcvb[1] = 0.0;
    
    ASSERT_EQ(dgridlb, CXCGradientGrid(CMemBlock2D<double>({3.0, 0.0, 3.0, 5.0}, 2, 2), dengrid::limb, xcfun::lda));
}

TEST_F(CXCGradientGridTest, XCGradientValuesConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);
    
    const CXCGradientGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    const CXCGradientGrid dgridla(mblockab, dengrid::lima, xcfun::lda);
    
    const CXCGradientGrid dgridlb(mblockab, dengrid::limb, xcfun::lda);
    
    // rho_ab case
    
    vlxtest::compare({3.0, 6.0}, dgridab.xcGradientValues(xcvars::rhoa));
    
    vlxtest::compare({1.0, 7.0}, dgridab.xcGradientValues(xcvars::rhob));
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::gradab) == nullptr);
    
    // lim_a case

    vlxtest::compare({3.0, 6.0}, dgridla.xcGradientValues(xcvars::rhoa));

    vlxtest::compare({1.0, 7.0}, dgridla.xcGradientValues(xcvars::rhob));
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradab) == nullptr);
    
    // lim_b case
    
    vlxtest::compare({3.0, 6.0}, dgridlb.xcGradientValues(xcvars::rhoa));

    vlxtest::compare({1.0, 7.0}, dgridlb.xcGradientValues(xcvars::rhob));
        
    ASSERT_TRUE(dgridlb.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcGradientValues(xcvars::gradab) == nullptr);
}

TEST_F(CXCGradientGridTest, XCGradientValues)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);

    
    CXCGradientGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    CXCGradientGrid dgridla(mblockab, dengrid::lima, xcfun::lda);
    
    CXCGradientGrid dgridlb(mblockab, dengrid::limb, xcfun::lda);
    
    // rho_ab case
    
    auto vxab_a = dgridab.xcGradientValues(xcvars::rhoa);
    
    auto vxab_b = dgridab.xcGradientValues(xcvars::rhob);
    
    vlxtest::compare({3.0, 6.0}, vxab_a);
    
    vlxtest::compare({1.0, 7.0}, vxab_b);
    
    vxab_a[1] = 9.0;
    
    vxab_b[0] = 1.5;
    
    ASSERT_EQ(dgridab, CXCGradientGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 9.0, 1.5, 7.0}, 2, 3), dengrid::ab, xcfun::lda));
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridab.xcGradientValues(xcvars::gradab) == nullptr);
    
    // lim_a case
    
    auto vxa_a = dgridla.xcGradientValues(xcvars::rhoa);
    
    auto vxa_b = dgridla.xcGradientValues(xcvars::rhob);
    
    vlxtest::compare({3.0, 6.0}, vxa_a);
    
    vlxtest::compare({1.0, 7.0}, vxa_b);
    
    vxa_a[1] = 9.0;
    
    vxa_b[0] = 1.5;
    
    ASSERT_EQ(dgridla, CXCGradientGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 9.0, 1.5, 7.0}, 2, 3), dengrid::lima, xcfun::lda));
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradab) == nullptr);
    
    // lim_b case
    
    auto vxb_a = dgridlb.xcGradientValues(xcvars::rhoa);
    
    auto vxb_b = dgridlb.xcGradientValues(xcvars::rhob);
    
    vlxtest::compare({3.0, 6.0}, vxb_a);
    
    vlxtest::compare({1.0, 7.0}, vxb_b);
    
    vxb_a[1] = 9.0;
    
    vxb_b[0] = 1.5;
    
    ASSERT_EQ(dgridlb, CXCGradientGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 9.0, 1.5, 7.0}, 2, 3), dengrid::limb, xcfun::lda));
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::grada) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradb) == nullptr);
    
    ASSERT_TRUE(dgridla.xcGradientValues(xcvars::gradab) == nullptr);
}
