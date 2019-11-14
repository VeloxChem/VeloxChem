//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularGridTest.hpp"

#include "MolecularGrid.hpp"
#include "CheckFunctions.hpp"

TEST_F(CMolecularGridTest, DefaultConstructor)
{
    CMolecularGrid mgrida;

    CMemBlock2D<double> mblock;

    CMolecularGrid mgridb(mblock);

    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, CopyConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);

    CMolecularGrid mgrida(mblock);

    CMolecularGrid mgridb(mgrida);

    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);

    CMolecularGrid mgrida(mblock);

    CMolecularGrid mgridb((CMolecularGrid(mblock)));

    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, CopyAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);

    CMolecularGrid mgrida(mblock);

    CMolecularGrid mgridb = mgrida;

    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, MoveAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);

    CMolecularGrid mgrida(mblock);

    CMolecularGrid mgridb = CMolecularGrid(mblock);

    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, Slice)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    mgrida.slice(4);
    
    ASSERT_EQ(mgrida, CMolecularGrid(mblock));
    
    mgrida.slice(1); 
    
    ASSERT_EQ(mgrida, CMolecularGrid(CMemBlock2D<double>({1.0, 3.0, 2.0, 5.0}, 1, 4)));
}

TEST_F(CMolecularGridTest, GetNumberOfGridPoints)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    ASSERT_EQ(2, mgrida.getNumberOfGridPoints());
}

TEST_F(CMolecularGridTest, GetCoordinatesXConstant)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    const CMolecularGrid mgrida(mblock);
    
    vlxtest::compare({1.0, 2.0}, mgrida.getCoordinatesX());
}

TEST_F(CMolecularGridTest, GetCoordinatesX)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    auto rx = mgrida.getCoordinatesX();
    
    vlxtest::compare({1.0, 2.0}, rx);
    
    rx[1] = 3.0;
    
    ASSERT_EQ(mgrida, CMolecularGrid(CMemBlock2D<double>({1.0, 3.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4)));
}

TEST_F(CMolecularGridTest, GetCoordinatesYConstant)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    const CMolecularGrid mgrida(mblock);
    
    vlxtest::compare({3.0, 6.0}, mgrida.getCoordinatesY());
}

TEST_F(CMolecularGridTest, GetCoordinatesY)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    auto ry = mgrida.getCoordinatesY();
    
    vlxtest::compare({3.0, 6.0}, ry);
    
    ry[0] = 5.0;
    
    ASSERT_EQ(mgrida, CMolecularGrid(CMemBlock2D<double>({1.0, 2.0, 5.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4)));
}

TEST_F(CMolecularGridTest, GetCoordinatesZConstant)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    const CMolecularGrid mgrida(mblock);
    
    vlxtest::compare({2.0, 4.0}, mgrida.getCoordinatesZ());
}

TEST_F(CMolecularGridTest, GetCoordinatesZ)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    auto rz = mgrida.getCoordinatesZ();
    
    vlxtest::compare({2.0, 4.0}, rz);
    
    rz[0] = 1.0; rz[1] = 3.0;
    
    ASSERT_EQ(mgrida, CMolecularGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 9.0}, 2, 4)));
}

TEST_F(CMolecularGridTest, GetWeightsConstant)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    const CMolecularGrid mgrida(mblock);
    
    vlxtest::compare({5.0, 9.0}, mgrida.getWeights());
}

TEST_F(CMolecularGridTest, GetWeights)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    auto w = mgrida.getWeights();
    
    vlxtest::compare({5.0, 9.0}, w);
    
    w[1] = 2.0;
    
    ASSERT_EQ(mgrida, CMolecularGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 2.0}, 2, 4)));
}

TEST_F(CMolecularGridTest, GetSpatialExtent)
{
    CMemBlock2D<double> mblock({ 1.0, 2.0, 1.0, 3.0,
                                -6.0, 2.0, 2.0, 4.0,
                                 5.0, 0.0, 9.0, 3.0,
                                 2.0, 4.0, 5.0, 2.0,}, 4, 4);
    
    CMolecularGrid mgrida(mblock);
    
    auto rext = mgrida.getSpatialExtent();
    
    ASSERT_NEAR( 1.0, rext[0], 1.0e-13);
    
    ASSERT_NEAR(-6.0, rext[1], 1.0e-13);
    
    ASSERT_NEAR( 0.0, rext[2], 1.0e-13);
    
    ASSERT_NEAR( 3.0, rext[3], 1.0e-13);
    
    ASSERT_NEAR( 4.0, rext[4], 1.0e-13);
    
    ASSERT_NEAR( 9.0, rext[5], 1.0e-13);
}
