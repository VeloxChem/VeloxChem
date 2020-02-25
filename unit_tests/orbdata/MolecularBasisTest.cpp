//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "MolecularBasisTest.hpp"

#include "AtomBasisSetter.hpp"
#include "CheckFunctions.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CMolecularBasisTest, DefaultConstructor)
{
    CMolecularBasis ambas;

    CMolecularBasis bmbas = vlxbas::getMolecularBasisEmpty();

    ASSERT_EQ(ambas, bmbas);
}

TEST_F(CMolecularBasisTest, CopyConstructor)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    CMolecularBasis bmbas(ambas);

    ASSERT_EQ(ambas, bmbas);
}

TEST_F(CMolecularBasisTest, MoveConstructor)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    CMolecularBasis bmbas(vlxbas::getMolecularBasisForLiH());

    ASSERT_EQ(ambas, bmbas);
}

TEST_F(CMolecularBasisTest, CopyAssignment)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    CMolecularBasis bmbas = ambas;

    ASSERT_EQ(ambas, bmbas);
}

TEST_F(CMolecularBasisTest, MoveAssignment)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(ambas, vlxbas::getMolecularBasisForLiH());
}

TEST_F(CMolecularBasisTest, SetAndGetMaxAngularMomentum)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(1, ambas.getMaxAngularMomentum());

    ambas.setMaxAngularMomentum(3);

    ASSERT_EQ(3, ambas.getMaxAngularMomentum());
}

TEST_F(CMolecularBasisTest, AddAtomBasis)
{
    CMolecularBasis ambas;

    ambas.setLabel({"def2-SVP"});

    ambas.addAtomBasis(vlxbas::getAtomBasisForLi());

    ambas.addAtomBasis(vlxbas::getAtomBasisForH());

    ASSERT_EQ(ambas, vlxbas::getMolecularBasisForLiH());
}

TEST_F(CMolecularBasisTest, ReduceToValenceBasis)
{
    CMolecularBasis ambas;

    ambas.setLabel({"def2-SVP(VAL)"});

    auto libas = vlxbas::getAtomBasisForLi();

    auto hbas = vlxbas::getAtomBasisForH();

    ambas.addAtomBasis(libas.reduceToValenceBasis());

    ambas.addAtomBasis(hbas.reduceToValenceBasis());

    auto bmbas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(ambas, bmbas.reduceToValenceBasis());
}

TEST_F(CMolecularBasisTest, SetAndGetLabel)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(std::string("def2-SVP"), ambas.getLabel());

    ambas.setLabel("def2-TZVP");

    ASSERT_EQ(std::string("def2-TZVP"), ambas.getLabel());
}

TEST_F(CMolecularBasisTest, GetMaxAngularMomentumWithIdElemental)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(1, ambas.getMaxAngularMomentum(1));

    ASSERT_EQ(1, ambas.getMaxAngularMomentum(3));

    ASSERT_EQ(-1, ambas.getMaxAngularMomentum(2));
}

TEST_F(CMolecularBasisTest, GetMolecularMaxAngularMomentum)
{
    CMolecule h2o = vlxmol::getMoleculeH2O();

    CMolecularBasis ao_basis = vlxbas::getMolecularBasisForH2O();

    CMolecularBasis min_basis = vlxbas::getMinimalBasisForH2O();

    ASSERT_EQ(2, ao_basis.getMolecularMaxAngularMomentum(h2o));

    ASSERT_EQ(1, min_basis.getMolecularMaxAngularMomentum(h2o));
}

TEST_F(CMolecularBasisTest, GetNumberOfBasisFunctions)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(2, ambas.getNumberOfBasisFunctions(1, 0));

    ASSERT_EQ(1, ambas.getNumberOfBasisFunctions(1, 1));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(1, 2));

    ASSERT_EQ(3, ambas.getNumberOfBasisFunctions(3, 0));

    ASSERT_EQ(2, ambas.getNumberOfBasisFunctions(3, 1));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(3, 2));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(2, 0));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(2, 1));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(2, 2));
}

TEST_F(CMolecularBasisTest, GetNumberOfBasisFunctionsWithMolecule)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(5, ambas.getNumberOfBasisFunctions(lih, 0));

    ASSERT_EQ(3, ambas.getNumberOfBasisFunctions(lih, 1));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(lih, 2));
}

TEST_F(CMolecularBasisTest, GetNumberOfBasisFunctionsWithMoleculeAndAtomsList)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(3, ambas.getNumberOfBasisFunctions(lih, 0, 1, 0));

    ASSERT_EQ(2, ambas.getNumberOfBasisFunctions(lih, 1, 1, 0));

    ASSERT_EQ(5, ambas.getNumberOfBasisFunctions(lih, 0, 2, 0));

    ASSERT_EQ(2, ambas.getNumberOfBasisFunctions(lih, 0, 1, 1));

    ASSERT_EQ(1, ambas.getNumberOfBasisFunctions(lih, 1, 1, 1));

    ASSERT_EQ(3, ambas.getNumberOfBasisFunctions(lih, 0, 2, 1));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(lih, 0, 1, 2));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(lih, 1, 1, 2));

    ASSERT_EQ(0, ambas.getNumberOfBasisFunctions(lih, 0, 2, 2));
}

TEST_F(CMolecularBasisTest, GetNumberOfPrimitiveBasisFunctionsWithMolecule)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(11, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0));

    ASSERT_EQ(4, ambas.getNumberOfPrimitiveBasisFunctions(lih, 1));

    ASSERT_EQ(0, ambas.getNumberOfPrimitiveBasisFunctions(lih, 2));
}

TEST_F(CMolecularBasisTest, GetNumberOfPrimitiveBasisFunctionsWithMoleculeAndAtomsList)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(7, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 1, 0));

    ASSERT_EQ(4, ambas.getNumberOfPrimitiveBasisFunctions(lih, 1, 1, 0));

    ASSERT_EQ(11, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 2, 0));

    ASSERT_EQ(3, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 1, 1));

    ASSERT_EQ(1, ambas.getNumberOfPrimitiveBasisFunctions(lih, 1, 1, 1));

    ASSERT_EQ(4, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 2, 1));

    ASSERT_EQ(0, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 1, 2));

    ASSERT_EQ(0, ambas.getNumberOfPrimitiveBasisFunctions(lih, 1, 1, 2));

    ASSERT_EQ(0, ambas.getNumberOfPrimitiveBasisFunctions(lih, 0, 2, 2));
}

TEST_F(CMolecularBasisTest, GetDimensionsOfBasis)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(14, ambas.getDimensionsOfBasis(lih));
}

TEST_F(CMolecularBasisTest, GetPartialDimensionsOfBasis)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(0, ambas.getPartialDimensionsOfBasis(lih, 0));

    ASSERT_EQ(5, ambas.getPartialDimensionsOfBasis(lih, 1));

    ASSERT_EQ(14, ambas.getPartialDimensionsOfBasis(lih, 2));

    ASSERT_EQ(14, ambas.getPartialDimensionsOfBasis(lih, 3));
}

TEST_F(CMolecularBasisTest, GetDimensionsOfPrimitiveBasis)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    ASSERT_EQ(23, ambas.getDimensionsOfPrimitiveBasis(lih));
}

TEST_F(CMolecularBasisTest, GetAtomBasis)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    ASSERT_EQ(ambas.getAtomBasis(1), vlxbas::getAtomBasisForH());

    ASSERT_EQ(ambas.getAtomBasis(3), vlxbas::getAtomBasisForLi());

    ASSERT_EQ(ambas.getAtomBasis(5), vlxbas::getAtomBasisEmpty());
}

TEST_F(CMolecularBasisTest, GetBasisFunctions)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();

    // hydrogen atom

    auto bfa = mbas.getBasisFunctions(1, 0);

    ASSERT_EQ(2u, bfa.size());

    ASSERT_EQ(bfa[0],
              CBasisFunction(
                  {1.301070100000e+01, 1.962257200000e+00, 4.445379600000e-01}, {1.968215800000e-02, 1.379652400000e-01, 4.783193500000e-01}, 0));

    ASSERT_EQ(bfa[1], CBasisFunction({1.219496200000e-01}, {1.000000000000e+00}, 0));

    bfa = mbas.getBasisFunctions(1, 1);

    ASSERT_EQ(1u, bfa.size());

    ASSERT_EQ(bfa[0], CBasisFunction({8.000000000000e-01}, {1.000000000000e+00}, 1));

    bfa = mbas.getBasisFunctions(1, 2);

    ASSERT_EQ(0u, bfa.size());

    // boron atom

    bfa = mbas.getBasisFunctions(5, 0);

    ASSERT_EQ(0u, bfa.size());

    bfa = mbas.getBasisFunctions(5, 1);

    ASSERT_EQ(0u, bfa.size());

    bfa = mbas.getBasisFunctions(5, 2);

    ASSERT_EQ(0u, bfa.size());

    // lithium atom

    bfa = mbas.getBasisFunctions(3, 0);

    ASSERT_EQ(3u, bfa.size());

    ASSERT_EQ(bfa[0],
              CBasisFunction({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00, 2.450300905100e+00, 7.220957185500e-01},
                             {6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01, 4.860657481700e-01, 4.362697795500e-01},
                             0));

    ASSERT_EQ(bfa[1], CBasisFunction({5.281088472100e-02}, {1.000000000000e+00}, 0));

    ASSERT_EQ(bfa[2], CBasisFunction({2.096094879800e-02}, {1.000000000000e+00}, 0));

    bfa = mbas.getBasisFunctions(3, 1);

    ASSERT_EQ(2u, bfa.size());

    ASSERT_EQ(bfa[0], CBasisFunction({1.450000000000e+00, 3.000000000000e-01}, {2.586000000000e-01, 1.000000000000e+00}, 1));

    ASSERT_EQ(bfa[1], CBasisFunction({8.200000000000e-02}, {1.000000000000e+00}, 1));

    bfa = mbas.getBasisFunctions(3, 2);

    ASSERT_EQ(0u, bfa.size());
}

TEST_F(CMolecularBasisTest, GetAOBasisMap)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    auto strmap = ambas.getAOBasisMap(lih);

    ASSERT_EQ(14, strmap.size());

    ASSERT_EQ(std::string("   1 Li  1s  "), strmap[0]);

    ASSERT_EQ(std::string("   1 Li  2s  "), strmap[1]);

    ASSERT_EQ(std::string("   1 Li  3s  "), strmap[2]);

    ASSERT_EQ(std::string("   2 H   1s  "), strmap[3]);

    ASSERT_EQ(std::string("   2 H   2s  "), strmap[4]);

    ASSERT_EQ(std::string("   1 Li  1p-1"), strmap[5]);

    ASSERT_EQ(std::string("   1 Li  2p-1"), strmap[6]);

    ASSERT_EQ(std::string("   2 H   1p-1"), strmap[7]);

    ASSERT_EQ(std::string("   1 Li  1p0 "), strmap[8]);

    ASSERT_EQ(std::string("   1 Li  2p0 "), strmap[9]);

    ASSERT_EQ(std::string("   2 H   1p0 "), strmap[10]);

    ASSERT_EQ(std::string("   1 Li  1p+1"), strmap[11]);

    ASSERT_EQ(std::string("   1 Li  2p+1"), strmap[12]);

    ASSERT_EQ(std::string("   2 H   1p+1"), strmap[13]);
}

TEST_F(CMolecularBasisTest, GetIndexMapForDalton)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    auto idsmap = ambas.getIndexMapForDalton(lih);
    
    ASSERT_EQ(2, idsmap.blocks());
    
    ASSERT_EQ(14, idsmap.size(0));
    
    ASSERT_EQ(14, idsmap.size(1));
    
    vlxtest::compare({0, 1, 2, 11, 5, 8, 12, 6, 9, 3, 4, 13, 7, 10}, idsmap.data(0));
    
    vlxtest::compare({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, idsmap.data(1));
    
    std::cout << idsmap << std::endl; 
}

TEST_F(CMolecularBasisTest, GetPositionInAngularBlock)
{
    CMolecularBasis ambas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    ASSERT_EQ( 0, ambas.getPositionInAngularBlock(lih, 0, 0));
    
    ASSERT_EQ( 3, ambas.getPositionInAngularBlock(lih, 1, 0));
    
    ASSERT_EQ(-1, ambas.getPositionInAngularBlock(lih, 2, 0));
    
    ASSERT_EQ( 0, ambas.getPositionInAngularBlock(lih, 0, 1));
    
    ASSERT_EQ( 2, ambas.getPositionInAngularBlock(lih, 1, 1));
    
    ASSERT_EQ(-1, ambas.getPositionInAngularBlock(lih, 2, 1));
    
    ASSERT_EQ(-1, ambas.getPositionInAngularBlock(lih, 0, 2));
    
    ASSERT_EQ(-1, ambas.getPositionInAngularBlock(lih, 1, 2));
    
    ASSERT_EQ(-1, ambas.getPositionInAngularBlock(lih, 2, 2));
}
