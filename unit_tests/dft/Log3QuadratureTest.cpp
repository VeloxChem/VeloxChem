//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "Log3QuadratureTest.hpp"

#include "Log3Quadrature.hpp"
#include "CheckFunctions.hpp"

TEST_F(CLog3QuadratureTest, ConstructorWith5Points)
{
    CLog3Quadrature rquad(5, 1);

    auto qpoints = rquad.generate();

    ASSERT_EQ(2, qpoints.blocks());

    ASSERT_EQ(5, qpoints.size(0));
    
    ASSERT_EQ(5, qpoints.size(1));

    std::vector<double> rad({4.536298573026920000, 2.040679201001270000,
                             0.800000000000003000, 0.219058105426335000,
                             0.023957400390054500});

    vlxtest::compare(rad, qpoints.data(0));

    std::vector<double> wgt({3.661009388224700000, 1.705129857584120000,
                             0.855642097864144000, 0.349387213052902000,
                             0.076565310395870600});

    vlxtest::compare(wgt, qpoints.data(1));
}
