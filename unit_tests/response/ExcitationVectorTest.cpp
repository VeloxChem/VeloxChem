//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ExcitationVectorTest.hpp"

#include "ExcitationVector.hpp"
#include "CheckFunctions.hpp"

TEST_F(CExcitationVectorTest, DefaultConstructor)
{
    CExcitationVector ma;
    
    CExcitationVector mb(szblock::aa, std::vector<int32_t>{}, std::vector<int32_t>{},
                         std::vector<double>{}, std::vector<double>{});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, ConstructorWithPairs)
{
    CExcitationVector ma(szblock::aa, 0, 2, 2, 5);
    
    std::vector<double> coef({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    CExcitationVector mb(szblock::aa, {0, 0, 0, 1, 1, 1}, {2, 3, 4, 2, 3, 4},
                         coef, coef);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, CopyConstructor)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb(ma);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, MoveConstructor)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb(CExcitationVector(szblock::aa, {0, 0, 1}, {2, 3, 2},
                                           {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0}));
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, CopyAssignment)
{
    CExcitationVector ma(szblock::ab, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb = ma;
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, MoveAssignment)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb = CExcitationVector(szblock::bb, {0, 0, 1}, {2, 3, 2},
                                             {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, SetCoefficientsZY)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    
    ma.setCoefficientsZY(CMemBlock<double>({1.0, 2.0, 3.4}),
                         CMemBlock<double>({2.3, 1.0, 4.0}));
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, 2.0, 3.4}, {2.3, 1.0, 4.0});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, GetCoefficientsZ)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    auto zdat = ma.getCoefficientsZ();
    
    vlxtest::compare({1.0, -1.0, -3.0}, zdat);
    
    zdat[1] = 2.3;
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, 2.3, -3.0}, {2.0, 3.0, -2.0});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, GetCoefficientsZConstant)
{
    const CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                               {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({1.0, -1.0, -3.0}, ma.getCoefficientsZ());
}

TEST_F(CExcitationVectorTest, GetCoefficientsY)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    auto ydat = ma.getCoefficientsY();
    
    vlxtest::compare({2.0, 3.0, -2.0}, ydat);
    
    ydat[0] = 0.7;
    
    ydat[2] = 1.7;
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {0.7, 3.0, 1.7});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, GetCoefficientsYConstant)
{
    const CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                               {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({2.0, 3.0, -2.0}, ma.getCoefficientsY());
}

TEST_F(CExcitationVectorTest, GetNumberOfExcitations)
{
    CExcitationVector ma(szblock::aa, 0, 5, 5, 10);
    
    ASSERT_EQ(25, ma.getNumberOfExcitations());
}
