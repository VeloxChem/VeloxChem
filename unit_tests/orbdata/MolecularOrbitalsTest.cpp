//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularOrbitalsTest.hpp"

#include "MolecularOrbitals.hpp"

TEST_F(CMolecularOrbitalsTest, DefaultConstructor)
{
    CMolecularOrbitals moa;
    
    CMolecularOrbitals mob({ }, { }, molorb::rest);
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    std::vector<double> eb({3.0, 5.0});
    
    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);
    
    CMolecularOrbitals mob(moa);
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    std::vector<double> eb({3.0, 5.0});
    
    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);
    
    CMolecularOrbitals mob(CMolecularOrbitals({ma, mb}, {ea, eb}, molorb::unrest));
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);
    
    CMolecularOrbitals mob = moa;
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);
    
    CMolecularOrbitals mob = CMolecularOrbitals({ma}, {ea}, molorb::rest);
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, GetAODensityForRestrictedCase)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);
    
    CDenseMatrix ref01mat({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0},
                          3, 3);
    
    ASSERT_EQ(moa.getAODensity(2), CAODensityMatrix({ref01mat}, denmat::rest));
    
    CDenseMatrix ref02mat({2.0, -7.0, 2.0, -7.0, 29.0, 8.0, 2.0, 8.0, 52.0},
                          3, 3);
    
    ASSERT_EQ(moa.getAODensity(4), CAODensityMatrix({ref02mat}, denmat::rest));
    
    CDenseMatrix ref03mat({11.0, -19.0, 14.0, -19.0, 45.0, -8.0, 14.0, -8.0, 68},
                          3, 3);
    
    ASSERT_EQ(moa.getAODensity(6), CAODensityMatrix({ref03mat}, denmat::rest));
}

TEST_F(CMolecularOrbitalsTest, GetAODensityForUnrestrictedCase)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    std::vector<double> ea({1.0, 2.0, 4.0});
    
    std::vector<double> eb({3.0, 5.0});
    
    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);
    
    CDenseMatrix ref01mata({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0},
                           3, 3);
    
    CDenseMatrix ref02mata({2.0, -7.0, 2.0, -7.0, 29.0, 8.0, 2.0, 8.0, 52.0},
                          3, 3);

    CDenseMatrix ref03mata({11.0, -19.0, 14.0, -19.0, 45.0, -8.0, 14.0, -8.0, 68},
                           3, 3);
    
    CDenseMatrix ref01matb({1.0, -3.0, 5.0, -3.0, 9.0, -15.0, 5.0, -15.0, 25.0},
                           3, 3);
    
    CDenseMatrix ref02matb({2.0, -1.0, 1.0, -1.0, 13.0, -23.0, 1.0, -23.0, 41.0},
                           3, 3);
    
    ASSERT_EQ(moa.getAODensity(2, 1), CAODensityMatrix({ref02mata, ref01matb}, denmat::unrest));
    
    ASSERT_EQ(moa.getAODensity(1, 1), CAODensityMatrix({ref01mata, ref01matb}, denmat::unrest));
    
    ASSERT_EQ(moa.getAODensity(2, 2), CAODensityMatrix({ref02mata, ref02matb}, denmat::unrest));
    
    ASSERT_EQ(moa.getAODensity(3, 2), CAODensityMatrix({ref03mata, ref02matb}, denmat::unrest));
}
