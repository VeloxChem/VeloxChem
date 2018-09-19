//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OutputLineTest.hpp"

#include "OutputLine.hpp"

TEST_F(COutputLineTest, DefaultConstructor)
{
    COutputLine aln(std::string(), 0, ' ', ' ', ' ', 0);

    COutputLine bln;

    ASSERT_EQ(aln, bln);
}

TEST_F(COutputLineTest, CopyConstructor)
{
    COutputLine aln(std::string("Velox"), 10, '=', '*', '_', 120);

    COutputLine bln(aln);

    ASSERT_EQ(aln, bln);
}

TEST_F(COutputLineTest, MoveConstructor)
{
    COutputLine aln(std::string("Velox"), 10, '=', '*', '_', 120);

    COutputLine bln(COutputLine(std::string("Velox"), 10, '=', '*', '_', 120));

    ASSERT_EQ(aln, bln);
}

TEST_F(COutputLineTest, CopyAssignment)
{
    COutputLine aln(std::string("Velox"), 10, '=', '*', '_', 120);

    COutputLine bln = aln;

    ASSERT_EQ(aln, bln);
}

TEST_F(COutputLineTest, MoveAssignment)
{
    COutputLine aln(std::string("Velox"), 10, '=', '*', '_', 120);

    COutputLine bln = COutputLine(std::string("Velox"), 10, '=', '*', '_', 120);

    ASSERT_EQ(aln, bln);
}
