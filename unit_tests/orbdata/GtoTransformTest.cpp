//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "GtoTransformTest.hpp"

#include "GtoTransform.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CGtoTransformTest, To_VeloxChem)
{
    auto ambas = vlxbas::getMolecularBasisForHeAtom();
    
    auto he = vlxmol::getMoleculeHeAtom();
    
    CDenseMatrix dalmat({1.0, 2.0, 3.0, 4.0, 5.0,
                         1.1, 2.1, 3.1, 4.1, 5.1,
                         1.2, 2.2, 3.2, 4.2, 5.2,
                         1.3, 2.3, 3.3, 4.3, 5.3,
                         1.4, 2.4, 3.4, 4.4, 5.4},
                        5, 5);
    
    auto vlxmat = gtotra::to_veloxchem(dalmat, ambas, he);
    
    CDenseMatrix refmat({1.0, 2.0, 4.0, 5.0, 3.0,
                         1.1, 2.1, 4.1, 5.1, 3.1,
                         1.3, 2.3, 4.3, 5.3, 3.3,
                         1.4, 2.4, 4.4, 5.4, 3.4,
                         1.2, 2.2, 4.2, 5.2, 3.2},
                        5, 5);
    
    ASSERT_EQ(vlxmat, refmat);
}

TEST_F(CGtoTransformTest, To_Dalton)
{
    auto ambas = vlxbas::getMolecularBasisForHeAtom();
    
    auto he = vlxmol::getMoleculeHeAtom();
    
    CDenseMatrix vlxmat({1.0, 2.0, 4.0, 5.0, 3.0,
                         1.1, 2.1, 4.1, 5.1, 3.1,
                         1.3, 2.3, 4.3, 5.3, 3.3,
                         1.4, 2.4, 4.4, 5.4, 3.4,
                         1.2, 2.2, 4.2, 5.2, 3.2},
                        5, 5);
    
    auto dalmat = gtotra::to_dalton(vlxmat, ambas, he);
    
    CDenseMatrix refmat({1.0, 2.0, 3.0, 4.0, 5.0,
                         1.1, 2.1, 3.1, 4.1, 5.1,
                         1.2, 2.2, 3.2, 4.2, 5.2,
                         1.3, 2.3, 3.3, 4.3, 5.3,
                         1.4, 2.4, 3.4, 4.4, 5.4},
                        5, 5);
    
    ASSERT_EQ(dalmat, refmat); 
}
