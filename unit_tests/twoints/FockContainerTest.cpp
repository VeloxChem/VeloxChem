//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "FockContainerTest.hpp"

#include "FockContainer.hpp"
#include "AOFockMatrix.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CFockContainerTest, DefaultConstructor)
{
    CFockContainer fconta;
    
    CFockContainer fcontb(std::vector<CFockSubMatrix>({}));
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, ConstructorWithAOFockMatrix)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restjk}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, kpairs);
    
    CMemBlock<int32_t> sidx(1); sidx.zero();
    
    CMemBlock<int32_t> pidx({5, 8, 11});
    
    CMemBlock2D<double> mss(25, 1); mss.zero();
    
    CMemBlock2D<double> msp(15, 3); msp.zero();
    
    CFockSubMatrix submb({mss, msp, mss, msp, mss, msp},
                         sidx, sidx, sidx, pidx, 5, 5, 5, 3);
    
    CFockContainer fcontb({submb});
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, CopyConstructor)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restj}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, kpairs);
    
    CFockContainer fcontb(fconta);
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, MoveConstructor)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restj}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, kpairs);
    
    CFockContainer fcontb(CFockContainer(&fmata, bpairs, kpairs));
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, CopyAssignment)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restj}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, kpairs);
    
    CFockContainer fcontb = fconta;
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, MoveAssignment)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restj}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, kpairs);
    
    CFockContainer fcontb = CFockContainer(&fmata, bpairs, kpairs);
    
    ASSERT_EQ(fconta, fcontb);
}

TEST_F(CFockContainerTest, GetSubMatrixData)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CAOFockMatrix fmata({ma}, {fockmat::restjkx}, {1.0}, {0});
    
    CFockContainer fconta(&fmata, bpairs, bpairs);
    
    std::vector<double> vecsp(15, 0.0);
    
    for (int32_t i = 0; i < 3; i++)
    {
        vlxtest::compare(vecsp, fconta.getSubMatrixData(0, 0, i));
        
        vlxtest::compare(vecsp, fconta.getSubMatrixData(0, 1, i));
        
        vlxtest::compare(vecsp, fconta.getSubMatrixData(0, 3, i));
        
        vlxtest::compare(vecsp, fconta.getSubMatrixData(0, 4, i));
    }
    
    std::vector<double> vecss(25, 0.0);
    
    vlxtest::compare(vecsp, fconta.getSubMatrixData(0, 2, 0));
    
    std::vector<double> vecpp(9, 0.0);
    
    for (int32_t i = 0; i < 9; i++)
    {
        vlxtest::compare(vecpp, fconta.getSubMatrixData(0, 5, i));
    }
}

TEST_F(CFockContainerTest, GetStartPositionsA)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    vlxtest::compare({0}, fconta.getStartPositionsA(0));
    
    vlxtest::compare({5, 8, 11}, fcontb.getStartPositionsA(0));
}

TEST_F(CFockContainerTest, GetStartPositionsB)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    vlxtest::compare({5, 8, 11}, fconta.getStartPositionsB(0));
    
    vlxtest::compare({5, 8, 11}, fcontb.getStartPositionsB(0));
}

TEST_F(CFockContainerTest, GetStartPositionsC)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    vlxtest::compare({5, 8, 11}, fconta.getStartPositionsC(0));
    
    vlxtest::compare({0}, fcontb.getStartPositionsC(0));
}

TEST_F(CFockContainerTest, GetStartPositionsD)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    vlxtest::compare({5, 8, 11}, fconta.getStartPositionsD(0));
    
    vlxtest::compare({5, 8, 11}, fcontb.getStartPositionsD(0));
}

TEST_F(CFockContainerTest, GetDimensionsA)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    ASSERT_EQ(5, fconta.getDimensionsA(0));
    
    ASSERT_EQ(3, fcontb.getDimensionsA(0));
}

TEST_F(CFockContainerTest, GetDimensionsB)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    ASSERT_EQ(3, fconta.getDimensionsB(0));
    
    ASSERT_EQ(3, fcontb.getDimensionsB(0));
}

TEST_F(CFockContainerTest, GetDimensionsC)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    ASSERT_EQ(3, fconta.getDimensionsC(0));
    
    ASSERT_EQ(5, fcontb.getDimensionsC(0));
}

TEST_F(CFockContainerTest, GetDimensionsD)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    CFockContainer fconta({subma});
    
    CFockContainer fcontb({submb});
    
    ASSERT_EQ(3, fconta.getDimensionsD(0));
    
    ASSERT_EQ(3, fcontb.getDimensionsD(0));
}
