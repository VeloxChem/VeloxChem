//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "StringFormatTest.hpp"

#include <string>

#include "FmtType.hpp"
#include "StringFormat.hpp"

TEST_F(CStringFormatTest, Upcase)
{
    auto ustr = fstr::upcase(std::string("Slater X"));

    ASSERT_EQ(ustr, std::string("SLATER X"));
}

TEST_F(CStringFormatTest, Format)
{
    std::string str("Slate");

    auto mstr = fstr::format(str, 9, fmt::center);

    ASSERT_EQ(mstr, std::string("  Slate  "));

    mstr = fstr::format(str, 9, fmt::left);

    ASSERT_EQ(mstr, std::string("Slate    "));

    mstr = fstr::format(str, 9, fmt::right);

    ASSERT_EQ(mstr, std::string("    Slate"));

    mstr = fstr::format(str, 2, fmt::center);

    ASSERT_EQ(mstr, std::string("Sl"));

    mstr = fstr::format(str, 3, fmt::left);

    ASSERT_EQ(mstr, std::string("Sla"));

    mstr = fstr::format(str, 4, fmt::right);

    ASSERT_EQ(mstr, std::string("Slat"));
}

TEST_F(CStringFormatTest, To_StringForCString)
{
    auto mstr = fstr::to_string("Slate", 9, fmt::center);

    ASSERT_EQ(mstr, std::string("  Slate  "));

    mstr = fstr::to_string("Slate", 9, fmt::left);

    ASSERT_EQ(mstr, std::string("Slate    "));

    mstr = fstr::to_string("Slate", 9, fmt::right);

    ASSERT_EQ(mstr, std::string("    Slate"));

    mstr = fstr::to_string("Slate", 2, fmt::center);

    ASSERT_EQ(mstr, std::string("Sl"));

    mstr = fstr::to_string("Slate", 3, fmt::left);

    ASSERT_EQ(mstr, std::string("Sla"));

    mstr = fstr::to_string("Slate", 4, fmt::right);

    ASSERT_EQ(mstr, std::string("Slat"));
}

TEST_F(CStringFormatTest, To_StringForCStringWithWidth)
{
    auto mstr = fstr::to_string("Slate", 9);

    ASSERT_EQ(mstr, std::string("  Slate  "));

    mstr = fstr::to_string("SlatE", 2, fmt::center);

    ASSERT_EQ(mstr, std::string("Sl"));
}

TEST_F(CStringFormatTest, To_StringForDouble)
{
    double f = 0.52917721067;

    auto mstr = fstr::to_string(f, 5, 12, fmt::center);

    ASSERT_EQ(mstr, std::string("   0.52918  "));

    mstr = fstr::to_string(-f, 5, 12, fmt::center);

    ASSERT_EQ(mstr, std::string("  -0.52918  "));

    mstr = fstr::to_string(f, 5, 12, fmt::left);

    ASSERT_EQ(mstr, std::string(" 0.52918    "));

    mstr = fstr::to_string(-f, 5, 12, fmt::left);

    ASSERT_EQ(mstr, std::string("-0.52918    "));

    mstr = fstr::to_string(f, 5, 12, fmt::right);

    ASSERT_EQ(mstr, std::string("     0.52918"));

    mstr = fstr::to_string(-f, 5, 12, fmt::right);

    ASSERT_EQ(mstr, std::string("    -0.52918"));
}

TEST_F(CStringFormatTest, To_StringForDoubleWithPrecision)
{
    double f = 0.52917721067;

    auto mstr = fstr::to_string(f, 3);

    ASSERT_EQ(mstr, std::string("0.529"));

    mstr = fstr::to_string(-f, 3);

    ASSERT_EQ(mstr, std::string("-0.529"));

    mstr = fstr::to_string(f, 2);

    ASSERT_EQ(mstr, std::string("0.53"));

    mstr = fstr::to_string(-f, 2);

    ASSERT_EQ(mstr, std::string("-0.53"));
}

TEST_F(CStringFormatTest, To_StringForInt32)
{
    int32_t f = 124;

    auto mstr = fstr::to_string(f, 6, fmt::center);

    ASSERT_EQ(mstr, std::string("  124 "));

    mstr = fstr::to_string(-f, 6, fmt::center);

    ASSERT_EQ(mstr, std::string(" -124 "));

    mstr = fstr::to_string(f, 6, fmt::left);

    ASSERT_EQ(mstr, std::string(" 124  "));

    mstr = fstr::to_string(-f, 6, fmt::left);

    ASSERT_EQ(mstr, std::string("-124  "));

    mstr = fstr::to_string(f, 6, fmt::right);

    ASSERT_EQ(mstr, std::string("   124"));

    mstr = fstr::to_string(-f, 6, fmt::right);

    ASSERT_EQ(mstr, std::string("  -124"));
}

TEST_F(CStringFormatTest, To_StringForSizeT)
{
    size_t f = 1024;

    auto mstr = fstr::to_string(f, 6, fmt::center);

    ASSERT_EQ(mstr, std::string(" 1024 "));

    mstr = fstr::to_string(f, 6, fmt::left);

    ASSERT_EQ(mstr, std::string("1024  "));

    mstr = fstr::to_string(f, 6, fmt::right);

    ASSERT_EQ(mstr, std::string("  1024"));
}

TEST_F(CStringFormatTest, To_String_Boolean)
{
    ASSERT_EQ(fstr::to_string(true), std::string("True"));

    ASSERT_EQ(fstr::to_string(false), std::string("False"));
}

TEST_F(CStringFormatTest, To_AngularMomentum)
{
    ASSERT_EQ(-1, fstr::to_AngularMomentum({"X"}));

    ASSERT_EQ(-1, fstr::to_AngularMomentum({"XZ"}));

    ASSERT_EQ(0, fstr::to_AngularMomentum({"S"}));

    ASSERT_EQ(0, fstr::to_AngularMomentum({"s"}));

    ASSERT_EQ(1, fstr::to_AngularMomentum({"P"}));

    ASSERT_EQ(1, fstr::to_AngularMomentum({"p"}));

    ASSERT_EQ(2, fstr::to_AngularMomentum({"D"}));

    ASSERT_EQ(2, fstr::to_AngularMomentum({"d"}));

    ASSERT_EQ(3, fstr::to_AngularMomentum({"F"}));

    ASSERT_EQ(3, fstr::to_AngularMomentum({"f"}));

    ASSERT_EQ(4, fstr::to_AngularMomentum({"G"}));

    ASSERT_EQ(4, fstr::to_AngularMomentum({"g"}));

    ASSERT_EQ(5, fstr::to_AngularMomentum({"H"}));

    ASSERT_EQ(5, fstr::to_AngularMomentum({"h"}));

    ASSERT_EQ(6, fstr::to_AngularMomentum({"I"}));

    ASSERT_EQ(6, fstr::to_AngularMomentum({"i"}));
}

TEST_F(CStringFormatTest, To_AngularMomentumFromNumber)
{
    ASSERT_EQ(std::string(), fstr::to_AngularMomentum(8));

    ASSERT_EQ(std::string("S"), fstr::to_AngularMomentum(0));

    ASSERT_EQ(std::string("P"), fstr::to_AngularMomentum(1));

    ASSERT_EQ(std::string("D"), fstr::to_AngularMomentum(2));

    ASSERT_EQ(std::string("F"), fstr::to_AngularMomentum(3));

    ASSERT_EQ(std::string("G"), fstr::to_AngularMomentum(4));

    ASSERT_EQ(std::string("H"), fstr::to_AngularMomentum(5));

    ASSERT_EQ(std::string("I"), fstr::to_AngularMomentum(6));
}
