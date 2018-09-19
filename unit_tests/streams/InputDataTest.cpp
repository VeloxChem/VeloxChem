//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.


#include "InputDataTest.hpp"

#include "InputData.hpp"
#include "ControlGroupSetter.hpp"

TEST_F(CInputDataTest, DefaultConstructor)
{
    CInputData idat;

    ASSERT_EQ(idat.getNumberOfControlGroups(), 0);
}

TEST_F(CInputDataTest, AddControlConstructor)
{
    CInputData idat;

    ASSERT_EQ(idat.getNumberOfControlGroups(), 0);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups(), 1);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups(), 2);
}

TEST_F(CInputDataTest, GetNumberOfControlGroupsForString)
{
    CInputData idat;

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("scf")), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("molxyz")), 0);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("scf")), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("molxyz")), 1);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("scf")), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups(std::string("molxyz")), 2);
}

TEST_F(CInputDataTest, GetNumberOfControlGroupsForCString)
{
    CInputData idat;

    ASSERT_EQ(idat.getNumberOfControlGroups("scf"), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups("molxyz"), 0);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups("scf"), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups("molxyz"), 1);

    idat.addControlGroup(getMolXYZGroup());

    ASSERT_EQ(idat.getNumberOfControlGroups("scf"), 0);

    ASSERT_EQ(idat.getNumberOfControlGroups("molxyz"), 2);
}

TEST_F(CInputDataTest, GetControlGroup)
{
    CInputData idat;

    idat.addControlGroup(getMolXYZGroup());

    idat.addControlGroup(getSCFGroup());

    size_t idx = 0;

    ASSERT_EQ(idat.getControlGroup(idx), getMolXYZGroup());

    idx = 1;

    ASSERT_EQ(idat.getControlGroup(idx), getSCFGroup());
}

TEST_F(CInputDataTest, GetControlGroupForString)
{
    CInputData idat;

    idat.addControlGroup(getMolXYZGroup());

    idat.addControlGroup(getSCFGroup());

    idat.addControlGroup(getMolXYZGroup());

    size_t idx = 0;

    ASSERT_EQ(idat.getControlGroup(idx, std::string("molxyz")), getMolXYZGroup());

    idx = 1;

    ASSERT_EQ(idat.getControlGroup(idx, std::string("molxyz")), getMolXYZGroup());
}

TEST_F(CInputDataTest, GetControlGroupForCString)
{
    CInputData idat;

    idat.addControlGroup(getMolXYZGroup());

    idat.addControlGroup(getSCFGroup());

    idat.addControlGroup(getMolXYZGroup());

    size_t idx = 0;

    ASSERT_EQ(idat.getControlGroup(idx, "molxyz"), getMolXYZGroup());

    idx = 1;

    ASSERT_EQ(idat.getControlGroup(idx, "molxyz"), getMolXYZGroup());
}
