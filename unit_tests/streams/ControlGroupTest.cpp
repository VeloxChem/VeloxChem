//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ControlGroupTest.hpp"

#include "ControlGroup.hpp"
#include "ControlGroupSetter.hpp"

TEST_F(CControlGroupTest, DefaultConstructor)
{
    CControlGroup acg;

    ASSERT_TRUE(acg.isEmpty());
}

TEST_F(CControlGroupTest, CopyConstructor)
{
    CControlGroup acg = getMolXYZGroup();

    CControlGroup bcg(acg);

    ASSERT_EQ(acg, bcg);
}

TEST_F(CControlGroupTest, MoveConstructor)
{
    CControlGroup acg = getMolXYZGroup();

    CControlGroup bcg(getMolXYZGroup());

    ASSERT_EQ(acg, bcg);
}

TEST_F(CControlGroupTest, CopyAssignment)
{
    CControlGroup acg = getMolXYZGroup();

    CControlGroup bcg = acg;

    ASSERT_EQ(acg, bcg);
}

TEST_F(CControlGroupTest, MoveAssignment)
{
    CControlGroup acg = getMolXYZGroup();

    ASSERT_EQ(acg, getMolXYZGroup());
}

TEST_F(CControlGroupTest, SetHeader)
{
    CControlGroup acg;

    acg.setHeader(CInputLine("@delta"));

    ASSERT_TRUE(acg.isNameOfControlGroup("delta"));
}

TEST_F(CControlGroupTest, AddCommand)
{
    CControlGroup acg;

    ASSERT_TRUE(acg.isEmpty());

    acg.addCommand(CInputLine("0 1 ! singlet"));

    ASSERT_EQ(acg.getNumberOfCommands(), 1);

    acg.addCommand(CInputLine("0 3 ! triplet"));

    ASSERT_EQ(acg.getNumberOfCommands(), 2);
}

TEST_F(CControlGroupTest, ClearAndIsEmpty)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_FALSE(acg.isEmpty());

    acg.clear();

    ASSERT_TRUE(acg.isEmpty());
}

TEST_F(CControlGroupTest, IsNameOfControlGroup)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_FALSE(acg.isNameOfControlGroup(std::string("molxyz")));

    ASSERT_TRUE(acg.isNameOfControlGroup(std::string("scf")));
}

TEST_F(CControlGroupTest, IsNameOfControlGroupFotCString)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_FALSE(acg.isNameOfControlGroup("molxyz"));

    ASSERT_TRUE(acg.isNameOfControlGroup("scf"));
}

TEST_F(CControlGroupTest, GetNumberOfCommands)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_EQ(acg.getNumberOfCommands(), 2);

    acg = getMolXYZGroup();

    ASSERT_EQ(acg.getNumberOfCommands(), 3);
}

TEST_F(CControlGroupTest, GetNumberOfCommandsForStrig)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_EQ(acg.getNumberOfCommands(std::string("MaxIterations")), 1);

    ASSERT_EQ(acg.getNumberOfCommands(std::string("Accelerator")), 0);
}

TEST_F(CControlGroupTest, GetNumberOfCommandsForCString)
{
    CControlGroup acg = getSCFGroup();

    ASSERT_EQ(acg.getNumberOfCommands("MaxIterations"), 1);

    ASSERT_EQ(acg.getNumberOfCommands("Accelerator"), 0);
}

TEST_F(CControlGroupTest, GetCommand)
{
    CControlGroup acg = getMolXYZGroup();

    size_t idx = 0;

    ASSERT_EQ(acg.getCommand(idx), CInputLine(std::string("0 1 ! charge multiplicity")));

    idx = 1;

    ASSERT_EQ(acg.getCommand(idx), CInputLine(std::string("Li 0.0 0.0 0.0")));

    idx = 2;

    ASSERT_EQ(acg.getCommand(idx), CInputLine(std::string("H  0.0 0.0 1.2")));
}

TEST_F(CControlGroupTest, GetCommandForString)
{
    CControlGroup acg = getMolXYZGroup();

    ASSERT_EQ(acg.getCommand(std::string("LI")), CInputLine(std::string("Li 0.0 0.0 0.0")));

    ASSERT_EQ(acg.getCommand(std::string("h")), CInputLine(std::string("H  0.0 0.0 1.2")));
}

TEST_F(CControlGroupTest, GetCommandForCString)
{
    CControlGroup acg = getMolXYZGroup();

    ASSERT_EQ(acg.getCommand("Li"), CInputLine(std::string("Li 0.0 0.0 0.0")));

    ASSERT_EQ(acg.getCommand("H"), CInputLine(std::string("H  0.0 0.0 1.2")));
}
