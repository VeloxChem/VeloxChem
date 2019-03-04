//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
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

TEST_F(CExcitationVectorTest, ConstructorWithPairsForTDA)
{
    CExcitationVector ma(szblock::aa, 0, 2, 2, 5, true);
    
    std::vector<double> coefz({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    std::vector<double> coefy;
    
    CExcitationVector mb(szblock::aa, {0, 0, 0, 1, 1, 1}, {2, 3, 4, 2, 3, 4},
                         coefz, coefy);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, ConstructorWithListForTDA)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {});
    
    CExcitationVector mb(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {2.0,  3.0, 4.0}, {});
    
    CExcitationVector mc(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {-2.0, 2.0, 2.0}, {});
    
    CExcitationVector md({2.0, 1.0, 3.0}, {ma, mb, mc});
    
    CExcitationVector mr(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {-2.0, 7.0, 4.0}, {});
    
    ASSERT_EQ(mr, md);
}

TEST_F(CExcitationVectorTest, ConstructorWithList)
{
    CExcitationVector ma(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, 5.0});
    
    CExcitationVector mb(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {2.0,  3.0, 4.0}, {1.0, -2.0, 2.0});
    
    CExcitationVector mc(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {-2.0, 2.0, 2.0}, {3.0, 4.0, 2.0});
    
    CExcitationVector md({2.0, 1.0, 3.0}, {ma, mb, mc});
    
    CExcitationVector mr(szblock::aa, {0, 0, 1}, {2, 3, 2},
                         {-2.0, 7.0, 4.0}, {14.0, 16.0, 18.0});
    
    ASSERT_EQ(mr, md);
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

TEST_F(CExcitationVectorTest, SetCoefficientZ)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    
    ma.setCoefficientZ(1.0, 0);
    
    ma.setCoefficientZ(3.4, 2);
    
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, 0.0, 3.4}, {0.0, 0.0, 0.0});
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CExcitationVectorTest, SetCoefficientY)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    
    ma.setCoefficientY(3.1, 1);
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {0.0, 0.0, 0.0}, {0.0, 3.1, 0.0});
    
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

TEST_F(CExcitationVectorTest, GetBraUniqueIndexes)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({0, 1}, ma.getBraUniqueIndexes());
}

TEST_F(CExcitationVectorTest, GetKetUniqueIndexes)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {4, 3, 4},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({3, 4}, ma.getKetUniqueIndexes());
}

TEST_F(CExcitationVectorTest, GetBraIndexes)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({0, 0, 1}, ma.getBraIndexes());
}

TEST_F(CExcitationVectorTest, GetKetIndexes)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    vlxtest::compare({2, 3, 2}, ma.getKetIndexes());
}

TEST_F(CExcitationVectorTest, GetMatrixZ)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 3},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CDenseMatrix refz({1.0, -1.0, 0.0, -3.0}, 2, 2);
    
    ASSERT_EQ(refz, ma.getMatrixZ());
}

TEST_F(CExcitationVectorTest, GetMatrixY)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 3},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CDenseMatrix refy({2.0, 0.0, 3.0, -2.0}, 2, 2);
    
    ASSERT_EQ(refy, ma.getMatrixY());
}

TEST_F(CExcitationVectorTest, GetDensityZ)
{
    CExcitationVector mvec(szblock::aa, {0, 0, 1}, {2, 3, 3},
                           {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);

    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    CDenseMatrix mden({  -7.0,  14.0,   8.0,   7.0,  13.0,
                         19.0, -38.0, -41.0,   8.0, -97.0,
                         10.0, -20.0,  -5.0, -19.0,   2.0,
                        -14.0,  28.0,  16.0,  14.0,  26.0,
                         13.0, -26.0, -32.0,  11.0, -79.0},
                       5, 5);
    
    CAODensityMatrix refden({mden}, denmat::rgen);
    
    ASSERT_EQ(refden, mvec.getDensityZ(mos));
}

TEST_F(CExcitationVectorTest, GetDensityY)
{
    CExcitationVector mvec(szblock::aa, {0, 0, 1}, {2, 3, 3},
                           {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);
    
    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    CDenseMatrix mden({ -16.0, -44.0,   52.0, -32.0, -48.0,
                         32.0,  88.0, -104.0,  64.0,  96.0,
                         19.0,  41.0,  -58.0,  38.0,  47.0,
                         15.0,  57.0,  -54.0,  30.0,  59.0,
                         32.0,  52.0,  -92.0,  64.0,  64.0},
                      5, 5);
    
    CAODensityMatrix refden({mden}, denmat::rgen);
    
    ASSERT_EQ(refden, mvec.getDensityY(mos));
}

TEST_F(CExcitationVectorTest, DotCoefficientsZ)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {2.0, 4.0, 1.0}, {2.0, 3.0, -2.0});
    
    ASSERT_NEAR(11.0, ma.dotCoefficientsZ(ma), 1.0e-13);
    
    ASSERT_NEAR(-5.0, ma.dotCoefficientsZ(mb), 1.0e-13);
    
    ASSERT_NEAR(-5.0, mb.dotCoefficientsZ(ma), 1.0e-13);
    
    ASSERT_NEAR(21.0, mb.dotCoefficientsZ(mb), 1.0e-13);
}

TEST_F(CExcitationVectorTest, DotCoefficientsY)
{
    CExcitationVector ma(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {1.0, -1.0, -3.0}, {2.0, 3.0, -2.0});
    
    CExcitationVector mb(szblock::bb, {0, 0, 1}, {2, 3, 2},
                         {2.0, 4.0, 1.0}, {1.0, 3.0, 1.0});
    
    ASSERT_NEAR(17.0, ma.dotCoefficientsY(ma), 1.0e-13);
    
    ASSERT_NEAR(9.0, ma.dotCoefficientsY(mb), 1.0e-13);
    
    ASSERT_NEAR(9.0, mb.dotCoefficientsY(ma), 1.0e-13);
    
    ASSERT_NEAR(11.0, mb.dotCoefficientsY(mb), 1.0e-13);
}

TEST_F(CExcitationVectorTest, GetBraEnergies)
{
    CExcitationVector mvec(szblock::ab, {0, 0, 1, 1}, {2, 3, 2, 3},
                           {1.0, -1.0, -3.0, 3.0}, {2.0, 3.0, -2.0, 1.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);
    
    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    vlxtest::compare(ea, mvec.getBraEnergies(mos));
}

TEST_F(CExcitationVectorTest, GetKetEnergies)
{
    CExcitationVector mvec(szblock::ba, {0, 0, 1, 1}, {2, 3, 2, 3},
                           {1.0, -1.0, -3.0, 3.0}, {2.0, 3.0, -2.0, 1.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);
    
    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    vlxtest::compare(ea, mvec.getKetEnergies(mos));
}

TEST_F(CExcitationVectorTest, GetSmallEnergyIdentifiers)
{
    CExcitationVector mvec(szblock::aa, {0, 0, 1, 1}, {2, 3, 2, 3},
                           {1.0, -1.0, -3.0, 3.0}, {2.0, 3.0, -2.0, 1.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);
    
    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    std::vector<int32_t> refidx4({0, 1, 2, 3});
    
    vlxtest::compare(refidx4, mvec.getSmallEnergyIdentifiers(mos, 4));
    
    std::vector<int32_t> refidx2({3, 1});
    
    vlxtest::compare(refidx2, mvec.getSmallEnergyIdentifiers(mos, 2));
}

TEST_F(CExcitationVectorTest, GetApproximateDiagonal)
{
    CExcitationVector mvec(szblock::aa, {0, 0, 1, 1}, {2, 3, 2, 3},
                           {1.0, -1.0, -3.0, 3.0}, {2.0, 3.0, -2.0, 1.0});
    
    CDenseMatrix ma({ 1.0, -1.0, -3.0, -2.0,
                      5.0,  4.0,  6.0,  4.0,
                     -4.0,  1.0,  2.0,  3.0,
                      2.0, -2.0,  5.0,  1.0,
                      5.0,  3.0,  1.0,  6.0},
                    5, 4);
    
    std::vector<double> ea({1.0, 2.0, 4.0, 2.0});
    
    const CMolecularOrbitals mos({ma}, {ea}, molorb::rest);
    
    CMemBlock<double> diagmat({3.0, 1.0, 2.0, 0.0});
    
    ASSERT_EQ(diagmat, mvec.getApproximateDiagonal(mos));
}
