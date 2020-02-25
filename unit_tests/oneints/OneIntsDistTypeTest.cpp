//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "OneIntsDistTypeTest.hpp"

#include "OneIntsDistType.hpp"

TEST_F(COneIntsDistTypeTest, To_String)
{
    ASSERT_EQ(to_string(dist1e::symsq), std::string("Symmetric Square Matrix"));

    ASSERT_EQ(to_string(dist1e::antisq), std::string("Anti-symmetric Square Matrix"));

    ASSERT_EQ(to_string(dist1e::rect), std::string("Rectangular Matrix"));

    ASSERT_EQ(to_string(dist1e::batch), std::string("Raw Integrals Batch"));
}
