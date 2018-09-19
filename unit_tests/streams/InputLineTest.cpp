//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "InputLineTest.hpp"

#include "InputLine.hpp"

TEST_F(CInputLineTest, DefaultConstructor)
{
    std::string str;

    CInputLine aln(str);

    CInputLine bln;

    ASSERT_EQ(aln, bln);
}

TEST_F(CInputLineTest, CopyConstructor)
{
    CInputLine aln(std::string(" Velox Chem "));

    CInputLine bln(aln);

    ASSERT_EQ(aln, bln);
}

TEST_F(CInputLineTest, MoveConstructor)
{
    CInputLine aln(std::string(" Velox Chem "));

    CInputLine bln(CInputLine(std::string(" Velox Chem ")));

    ASSERT_EQ(aln, bln);
}

TEST_F(CInputLineTest, CopyAssignment)
{
    CInputLine aln(std::string(" Velox Chem "));

    CInputLine bln = aln;

    ASSERT_EQ(aln, bln);
}

TEST_F(CInputLineTest, MoveAssignment)
{
    CInputLine aln(std::string(" Velox Chem "));

    CInputLine bln = CInputLine(std::string(" Velox Chem "));

    ASSERT_EQ(aln, bln);
}

TEST_F(CInputLineTest, GetKeyword)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_EQ(aln.getKeyword(0), std::string("Velox"));

    ASSERT_EQ(aln.getKeyword(1), std::string("Chem"));

    ASSERT_EQ(aln.getKeyword(2), std::string());
}

TEST_F(CInputLineTest, GetUpcasedKeyword)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_EQ(aln.getUpcasedKeyword(0), std::string("VELOX"));

    ASSERT_EQ(aln.getUpcasedKeyword(1), std::string("CHEM"));

    ASSERT_EQ(aln.getUpcasedKeyword(2), std::string());
}

TEST_F(CInputLineTest, GetNumberOfKeywords)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_EQ(aln.getNumberOfKeywords(), 2);
}

TEST_F(CInputLineTest, IsRealNumber)
{
    CInputLine aln(std::string(" Velox -2.0e-4 "));

    ASSERT_FALSE(aln.isRealNumber(0));

    ASSERT_TRUE(aln.isRealNumber(1));
}

TEST_F(CInputLineTest, GetRealNumber)
{
    CInputLine aln(std::string(" Velox -2.0e-4 "));

    ASSERT_NEAR(aln.getRealNumber(1), -2.0e-4, 1.0e-15);
}

TEST_F(CInputLineTest, IsIntegerNumber)
{
    CInputLine aln(std::string(" Velox 64 "));

    ASSERT_FALSE(aln.isIntegerNumber(0));

    ASSERT_TRUE(aln.isIntegerNumber(1));
}

TEST_F(CInputLineTest, GetIntegerNumber)
{
    CInputLine aln(std::string(" Velox 64 "));

    ASSERT_EQ(aln.getIntegerNumber(1), 64);
}

TEST_F(CInputLineTest, IsKeyword)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_TRUE(aln.isKeyword(0, std::string("VeloX")));

    ASSERT_FALSE(aln.isKeyword(0, std::string("chem")));

    ASSERT_FALSE(aln.isKeyword(1, std::string("VeloX")));

    ASSERT_TRUE(aln.isKeyword(1, std::string("chem")));
}

TEST_F(CInputLineTest, IsKeywordForCString)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_TRUE(aln.isKeyword(0, "VeloX"));

    ASSERT_FALSE(aln.isKeyword(0, "chem"));

    ASSERT_FALSE(aln.isKeyword(1, "VeloX"));

    ASSERT_TRUE(aln.isKeyword(1, "chem"));
}

TEST_F(CInputLineTest, IsControlKeyword)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_FALSE(aln.isControlKeyword(std::string("molxyz")));

    aln = CInputLine(" @Scf ! Scf input");

    ASSERT_FALSE(aln.isControlKeyword(std::string("molxyz")));

    ASSERT_TRUE(aln.isControlKeyword(std::string("scf")));
}

TEST_F(CInputLineTest, IsControlKeywordForCString)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_FALSE(aln.isControlKeyword("molxyz"));

    aln = CInputLine(" @Scf ! Scf input");

    ASSERT_FALSE(aln.isControlKeyword("molxyz"));

    ASSERT_TRUE(aln.isControlKeyword("scf"));
}

TEST_F(CInputLineTest, IsControlLine)
{
    CInputLine aln(std::string(" Velox Chem "));

    ASSERT_FALSE(aln.isControlLine());

    aln = CInputLine(" @Scf ! Scf input");

    ASSERT_TRUE(aln.isControlLine());
}

TEST_F(CInputLineTest, GetParsedString)
{
    CInputLine aln(" @Scf ! Scf input");

    ASSERT_EQ(aln.getParsedString(), std::string("@Scf"));
}

TEST_F(CInputLineTest, GetOriginalString)
{
    CInputLine aln(" @Scf ! Scf input");

    ASSERT_EQ(aln.getOriginalString(), std::string(" @Scf ! Scf input"));
}

TEST_F(CInputLineTest, ClearAndIsEmpty)
{
    CInputLine aln(" @Scf ! Scf input");

    aln.clear();

    CInputLine bln;

    ASSERT_EQ(aln, bln);

    ASSERT_TRUE(aln.isEmpty());

    ASSERT_TRUE(bln.isEmpty());
}
