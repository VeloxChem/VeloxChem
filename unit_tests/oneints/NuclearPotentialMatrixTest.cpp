//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "NuclearPotentialMatrixTest.hpp"

#include "NuclearPotentialMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CNuclearPotentialMatrixTest, DefaultConstructor)
{
    CNuclearPotentialMatrix smata;
    
    CNuclearPotentialMatrix smatb(CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb(CNuclearPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb = CNuclearPotentialMatrix({ma});
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, GetString)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CNuclearPotentialMatrix smata(ma);
    
    ASSERT_EQ(ma.getString(), smata.getString());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CNuclearPotentialMatrix smata(ma);
    
    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CNuclearPotentialMatrix smata(ma);
    
    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CNuclearPotentialMatrix smata(ma);
    
    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CNuclearPotentialMatrixTest, Values)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CNuclearPotentialMatrix smata(ma);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, smata.values());
}

TEST_F(CNuclearPotentialMatrixTest, GetKineticEnergy)
{
    CDenseMatrix ma({1.0, -1.0, -3.0,
                    -2.0,  5.0,  4.0,
                     6.0,  4.0, -4.0},
                    3, 3);
    
    CNuclearPotentialMatrix smat(ma);
    
    CDenseMatrix da({ 1.0, -1.0, -3.0,
                     -2.0,  2.0,  3.0,
                      6.0,  3.0, -4.0},
                    3, 3);
    
    CDenseMatrix db({ 1.0, -1.0, -3.0,
                     -2.0,  1.0,  4.0,
                      3.0,  7.0,  3.0},
                    3, 3);
    
    CAODensityMatrix dmat({da, db}, denmat::rest);
    
    ASSERT_NEAR(19.0, smat.getNuclearPotentialEnergy(dmat, 0), 1.0e-13);
    
    ASSERT_NEAR(15.0, smat.getNuclearPotentialEnergy(dmat, 1), 1.0e-13);
    
    ASSERT_NEAR(0.0, smat.getNuclearPotentialEnergy(dmat, 2), 1.0e-13);
}
