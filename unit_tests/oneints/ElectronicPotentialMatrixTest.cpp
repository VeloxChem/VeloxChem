//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronicPotentialMatrixTest.hpp"

#include "ElectronicPotentialMatrix.hpp"

TEST_F(CElectronicPotentialMatrixTest, DefaultConstructor)
{
    CElectronicPotentialMatrix smata;

    CElectronicPotentialMatrix smatb(CDenseMatrix(0, 0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb = CElectronicPotentialMatrix({ma});

    ASSERT_EQ(smata, smatb);
}
