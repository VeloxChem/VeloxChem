//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OneIntsDistTypeTest.hpp"

#include "OneIntsDistType.hpp"

TEST_F(COneIntsDistTypeTest, To_String)
{
    ASSERT_EQ(to_string(dist1e::symsq), std::string("Symmetric Square Matrix"));
    
    ASSERT_EQ(to_string(dist1e::antisq), std::string("Anti-symmetric Square Matrix"));
    
    ASSERT_EQ(to_string(dist1e::rect), std::string("Rectangular Matrix"));
    
    ASSERT_EQ(to_string(dist1e::batch), std::string("Raw Integrals Batch"));
}
