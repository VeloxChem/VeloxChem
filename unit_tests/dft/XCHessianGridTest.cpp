//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "XCHessianGridTest.hpp"

#include "XCHessianGrid.hpp"
#include "CheckFunctions.hpp"

TEST_F(CXCHessianGridTest, DefaultConstructor)
{
    CXCHessianGrid dgrida;
    
    CMemBlock2D<double> mblock;
    
    CXCHessianGrid dgridb(0, dengrid::undefined, xcfun::undefined);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCHessianGridTest, ConstructorWithFunctionalType)
{
    CXCHessianGrid dgrida(5, dengrid::ab, xcfun::lda);
    
    ASSERT_EQ(dgrida, CXCHessianGrid(CMemBlock2D<double>(5, 3), dengrid::ab, xcfun::lda));
    
    CXCHessianGrid dgridb(5, dengrid::lima, xcfun::lda);
    
    ASSERT_EQ(dgridb, CXCHessianGrid(CMemBlock2D<double>(5, 1), dengrid::lima, xcfun::lda));
    
    CXCHessianGrid dgridc(5, dengrid::limb, xcfun::lda);
    
    ASSERT_EQ(dgridc, CXCHessianGrid(CMemBlock2D<double>(5, 1), dengrid::limb, xcfun::lda));
    
    CXCHessianGrid dgride(5, dengrid::ab, xcfun::gga);
    
    ASSERT_EQ(dgride, CXCHessianGrid(CMemBlock2D<double>(5, 15), dengrid::ab, xcfun::gga));
    
    CXCHessianGrid dgridf(5, dengrid::lima, xcfun::gga);
    
    ASSERT_EQ(dgridf, CXCHessianGrid(CMemBlock2D<double>(5, 3), dengrid::lima, xcfun::gga));
    
    CXCHessianGrid dgridg(5, dengrid::limb, xcfun::gga);
    
    ASSERT_EQ(dgridg, CXCHessianGrid(CMemBlock2D<double>(5, 3), dengrid::limb, xcfun::gga));
    
    CXCHessianGrid dgridk(5, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(dgridk, CXCHessianGrid(CMemBlock2D<double>(5, 28), dengrid::ab, xcfun::mgga));
    
    CXCHessianGrid dgridl(5, dengrid::lima, xcfun::mgga);
    
    ASSERT_EQ(dgridl, CXCHessianGrid(CMemBlock2D<double>(5, 6), dengrid::lima, xcfun::mgga));
    
    CXCHessianGrid dgridm(5, dengrid::limb, xcfun::mgga);
    
    ASSERT_EQ(dgridm, CXCHessianGrid(CMemBlock2D<double>(5, 6), dengrid::limb, xcfun::mgga));
}

TEST_F(CXCHessianGridTest, CopyConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::lima, xcfun::mgga);
    
    CXCHessianGrid dgridb(dgrida);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCHessianGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    CXCHessianGrid dgridb((CXCHessianGrid(mblock, dengrid::ab, xcfun::mgga)));
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCHessianGridTest, CopyAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::limb, xcfun::mgga);
    
    CXCHessianGrid dgridb = dgrida;
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCHessianGridTest, MoveAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    CXCHessianGrid dgridb = CXCHessianGrid(mblock, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CXCHessianGridTest, Zero)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    dgrida.zero();
    
    ASSERT_EQ(dgrida, CXCHessianGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 4), dengrid::ab, xcfun::mgga));
}

TEST_F(CXCHessianGridTest, GetNumberOfGridPoints)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CXCHessianGrid dgrida(mblock, dengrid::ab, xcfun::mgga);
    
    ASSERT_EQ(2, dgrida.getNumberOfGridPoints());
}

TEST_F(CXCHessianGridTest, XCHessianValuesConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);
    
    CMemBlock2D<double> mblocka({3.0, 5.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({4.0, 7.0}, 2, 1);
    
    const CXCHessianGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    const CXCHessianGrid dgridla(mblockb, dengrid::lima, xcfun::lda);
    
    const CXCHessianGrid dgridlb(mblocka, dengrid::limb, xcfun::lda);
    
    // rho_a, rho_a case
    
    vlxtest::compare({1.0, 2.0}, dgridab.xcHessianValues(xcvars::rhoa, xcvars::rhoa));
    
    vlxtest::compare({3.0, 6.0}, dgridab.xcHessianValues(xcvars::rhoa, xcvars::rhob));
    
    vlxtest::compare({1.0, 7.0}, dgridab.xcHessianValues(xcvars::rhob, xcvars::rhob));
    
    ASSERT_TRUE(dgridab.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    // lim_a case
    
    vlxtest::compare({4.0, 7.0}, dgridla.xcHessianValues(xcvars::rhob, xcvars::rhob));
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhoa, xcvars::rhoa) == nullptr);
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhoa, xcvars::rhob) == nullptr);
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    // lim_b case
    
    vlxtest::compare({3.0, 5.0}, dgridlb.xcHessianValues(xcvars::rhoa, xcvars::rhoa));
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhoa, xcvars::rhob) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhob, xcvars::rhob) == nullptr);
}

TEST_F(CXCHessianGridTest, XCHessianValues)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 1.0, 7.0}, 2, 3);
    
    CMemBlock2D<double> mblocka({3.0, 5.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({4.0, 7.0}, 2, 1);
    
    CXCHessianGrid dgridab(mblockab, dengrid::ab, xcfun::lda);
    
    CXCHessianGrid dgridla(mblockb, dengrid::lima, xcfun::lda);
    
    CXCHessianGrid dgridlb(mblocka, dengrid::limb, xcfun::lda);
    
    // rho_a, rho_a case
    
    auto rab_aa = dgridab.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
    
    auto rab_ab = dgridab.xcHessianValues(xcvars::rhoa, xcvars::rhob);
    
    auto rab_bb = dgridab.xcHessianValues(xcvars::rhob, xcvars::rhob);
    
    vlxtest::compare({1.0, 2.0}, rab_aa);
    
    vlxtest::compare({3.0, 6.0}, rab_ab);
    
    vlxtest::compare({1.0, 7.0}, rab_bb);
    
    rab_aa[0] = 3.0; rab_ab[1] = 4.0; rab_bb[1] = 1.0;
    
    CMemBlock2D<double> nblockab({3.0, 2.0, 3.0, 4.0, 1.0, 1.0}, 2, 3);
    
    ASSERT_EQ(dgridab, CXCHessianGrid(nblockab, dengrid::ab, xcfun::lda));
    
    ASSERT_TRUE(dgridab.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    // lim_a case
    
    auto rla_bb = dgridla.xcHessianValues(xcvars::rhob, xcvars::rhob);
    
    vlxtest::compare({4.0, 7.0}, rla_bb);
    
    rla_bb[1] = 3.0;
    
    CMemBlock2D<double> nblockb({4.0, 3.0}, 2, 1);
    
    ASSERT_EQ(dgridla, CXCHessianGrid(nblockb, dengrid::lima, xcfun::lda));
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhoa, xcvars::rhoa) == nullptr);
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhoa, xcvars::rhob) == nullptr);
    
    ASSERT_TRUE(dgridla.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    // lim_b case
    
    auto rlb_aa = dgridlb.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
    
    vlxtest::compare({3.0, 5.0}, rlb_aa);
    
    rlb_aa[0] = 4.0;
    
    CMemBlock2D<double> nblocka({4.0, 5.0}, 2, 1);
    
    ASSERT_EQ(dgridlb, CXCHessianGrid(nblocka, dengrid::limb, xcfun::lda));
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhoa, xcvars::rhob) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhob, xcvars::rhoa) == nullptr);
    
    ASSERT_TRUE(dgridlb.xcHessianValues(xcvars::rhob, xcvars::rhob) == nullptr);
}
