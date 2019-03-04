//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentumTest.hpp"

#include "AngularMomentum.hpp"

TEST_F(CAngularMomentumTest, To_SphericalComponents)
{
    ASSERT_EQ(1, angmom::to_SphericalComponents(0));

    ASSERT_EQ(3, angmom::to_SphericalComponents(1));

    ASSERT_EQ(5, angmom::to_SphericalComponents(2));

    ASSERT_EQ(7, angmom::to_SphericalComponents(3));

    ASSERT_EQ(9, angmom::to_SphericalComponents(4));

    ASSERT_EQ(11, angmom::to_SphericalComponents(5));

    ASSERT_EQ(13, angmom::to_SphericalComponents(6));
}

TEST_F(CAngularMomentumTest, To_SphericalComponentsForPair)
{
    ASSERT_EQ(1, angmom::to_SphericalComponents(0, 0));
    
    ASSERT_EQ(3, angmom::to_SphericalComponents(1, 0));
    
    ASSERT_EQ(5, angmom::to_SphericalComponents(0, 2));
    
    ASSERT_EQ(9, angmom::to_SphericalComponents(1, 1));
    
    ASSERT_EQ(45, angmom::to_SphericalComponents(2, 4));
    
    ASSERT_EQ(45, angmom::to_SphericalComponents(4, 2));
}

TEST_F(CAngularMomentumTest, To_CartesianComponents)
{
    ASSERT_EQ(1, angmom::to_CartesianComponents(0));

    ASSERT_EQ(3, angmom::to_CartesianComponents(1));

    ASSERT_EQ(6, angmom::to_CartesianComponents(2));

    ASSERT_EQ(10, angmom::to_CartesianComponents(3));

    ASSERT_EQ(15, angmom::to_CartesianComponents(4));

    ASSERT_EQ(21, angmom::to_CartesianComponents(5));

    ASSERT_EQ(28, angmom::to_CartesianComponents(6));
}

TEST_F(CAngularMomentumTest, To_CartesianComponentsForPair)
{
    ASSERT_EQ(1, angmom::to_CartesianComponents(0, 0));
    
    ASSERT_EQ(3, angmom::to_CartesianComponents(1, 0));
    
    ASSERT_EQ(9, angmom::to_CartesianComponents(1, 1));
    
    ASSERT_EQ(36, angmom::to_CartesianComponents(2, 2));
    
    ASSERT_EQ(90, angmom::to_CartesianComponents(4, 2));
    
    ASSERT_EQ(90, angmom::to_CartesianComponents(2, 4));
}

TEST_F(CAngularMomentumTest, )
{
    ASSERT_EQ(std::string("s  "), angmom::getStringOfAngularMomentum(0, 0));
    
    ASSERT_EQ(std::string("p-1"), angmom::getStringOfAngularMomentum(1, 0));
    
    ASSERT_EQ(std::string("p0 "), angmom::getStringOfAngularMomentum(1, 1));
    
    ASSERT_EQ(std::string("p+1"), angmom::getStringOfAngularMomentum(1, 2));
    
    ASSERT_EQ(std::string("d-2"), angmom::getStringOfAngularMomentum(2, 0));
    
    ASSERT_EQ(std::string("d-1"), angmom::getStringOfAngularMomentum(2, 1));
    
    ASSERT_EQ(std::string("d0 "), angmom::getStringOfAngularMomentum(2, 2));
    
    ASSERT_EQ(std::string("d+1"), angmom::getStringOfAngularMomentum(2, 3));
    
    ASSERT_EQ(std::string("d+2"), angmom::getStringOfAngularMomentum(2, 4));
    
    ASSERT_EQ(std::string("f-3"), angmom::getStringOfAngularMomentum(3, 0));
    
    ASSERT_EQ(std::string("f-2"), angmom::getStringOfAngularMomentum(3, 1));
    
    ASSERT_EQ(std::string("f-1"), angmom::getStringOfAngularMomentum(3, 2));
    
    ASSERT_EQ(std::string("f0 "), angmom::getStringOfAngularMomentum(3, 3));
    
    ASSERT_EQ(std::string("f+1"), angmom::getStringOfAngularMomentum(3, 4));
    
    ASSERT_EQ(std::string("f+2"), angmom::getStringOfAngularMomentum(3, 5));
    
    ASSERT_EQ(std::string("f+3"), angmom::getStringOfAngularMomentum(3, 6));
    
    ASSERT_EQ(std::string("g-4"), angmom::getStringOfAngularMomentum(4, 0));
    
    ASSERT_EQ(std::string("g-3"), angmom::getStringOfAngularMomentum(4, 1));
    
    ASSERT_EQ(std::string("g-2"), angmom::getStringOfAngularMomentum(4, 2));
    
    ASSERT_EQ(std::string("g-1"), angmom::getStringOfAngularMomentum(4, 3));
    
    ASSERT_EQ(std::string("g0 "), angmom::getStringOfAngularMomentum(4, 4));
    
    ASSERT_EQ(std::string("g+1"), angmom::getStringOfAngularMomentum(4, 5));
    
    ASSERT_EQ(std::string("g+2"), angmom::getStringOfAngularMomentum(4, 6));
    
    ASSERT_EQ(std::string("g+3"), angmom::getStringOfAngularMomentum(4, 7));
    
    ASSERT_EQ(std::string("g+4"), angmom::getStringOfAngularMomentum(4, 8));
    
    ASSERT_EQ(std::string("h-5"), angmom::getStringOfAngularMomentum(5, 0));
    
    ASSERT_EQ(std::string("h-4"), angmom::getStringOfAngularMomentum(5, 1));
    
    ASSERT_EQ(std::string("h-3"), angmom::getStringOfAngularMomentum(5, 2));
    
    ASSERT_EQ(std::string("h-2"), angmom::getStringOfAngularMomentum(5, 3));
    
    ASSERT_EQ(std::string("h-1"), angmom::getStringOfAngularMomentum(5, 4));
    
    ASSERT_EQ(std::string("h0 "), angmom::getStringOfAngularMomentum(5, 5));
    
    ASSERT_EQ(std::string("h+1"), angmom::getStringOfAngularMomentum(5, 6));
    
    ASSERT_EQ(std::string("h+2"), angmom::getStringOfAngularMomentum(5, 7));
    
    ASSERT_EQ(std::string("h+3"), angmom::getStringOfAngularMomentum(5, 8));
    
    ASSERT_EQ(std::string("h+4"), angmom::getStringOfAngularMomentum(5, 9));
    
    ASSERT_EQ(std::string("h+5"), angmom::getStringOfAngularMomentum(5, 10));
    
    ASSERT_EQ(std::string("i-6"), angmom::getStringOfAngularMomentum(6, 0));
    
    ASSERT_EQ(std::string("i-5"), angmom::getStringOfAngularMomentum(6, 1));
    
    ASSERT_EQ(std::string("i-4"), angmom::getStringOfAngularMomentum(6, 2));
    
    ASSERT_EQ(std::string("i-3"), angmom::getStringOfAngularMomentum(6, 3));
    
    ASSERT_EQ(std::string("i-2"), angmom::getStringOfAngularMomentum(6, 4));
    
    ASSERT_EQ(std::string("i-1"), angmom::getStringOfAngularMomentum(6, 5));
    
    ASSERT_EQ(std::string("i0 "), angmom::getStringOfAngularMomentum(6, 6));
    
    ASSERT_EQ(std::string("i+1"), angmom::getStringOfAngularMomentum(6, 7));
    
    ASSERT_EQ(std::string("i+2"), angmom::getStringOfAngularMomentum(6, 8));
    
    ASSERT_EQ(std::string("i+3"), angmom::getStringOfAngularMomentum(6, 9));
    
    ASSERT_EQ(std::string("i+4"), angmom::getStringOfAngularMomentum(6, 10));
    
    ASSERT_EQ(std::string("i+5"), angmom::getStringOfAngularMomentum(6, 11));
    
    ASSERT_EQ(std::string("i+6"), angmom::getStringOfAngularMomentum(6, 12));
}
