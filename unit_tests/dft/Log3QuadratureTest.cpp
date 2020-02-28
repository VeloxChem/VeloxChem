//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "Log3QuadratureTest.hpp"

#include "CheckFunctions.hpp"
#include "Log3Quadrature.hpp"

TEST_F(CLog3QuadratureTest, ConstructorWith5Points)
{
    CLog3Quadrature rquad(5, 1);

    auto qpoints = rquad.generate();

    ASSERT_EQ(2, qpoints.blocks());

    ASSERT_EQ(5, qpoints.size(0));

    ASSERT_EQ(5, qpoints.size(1));

    std::vector<double> rad({4.536298573026920000, 2.040679201001270000, 0.800000000000003000, 0.219058105426335000, 0.023957400390054500});

    vlxtest::compare(rad, qpoints.data(0));

    std::vector<double> wgt({3.661009388224700000, 1.705129857584120000, 0.855642097864144000, 0.349387213052902000, 0.076565310395870600});

    vlxtest::compare(wgt, qpoints.data(1));
}
