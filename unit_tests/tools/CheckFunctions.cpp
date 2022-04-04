//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "CheckFunctions.hpp"

#include <gtest/gtest.h>

namespace vlxtest {  // namespace

void
compare(const std::vector<double>& aVector, const double* bVector)
{
    vlxtest::compare(aVector, bVector, 1.0e-13);
}

void
compare(const std::vector<double>& aVector, const double* bVector, const double threshod)
{
    for (size_t i = 0; i < aVector.size(); i++)
    {
        ASSERT_NEAR(aVector[i], bVector[i], threshod);
    }
}

void
compare(const std::vector<int32_t>& aVector, const int32_t* bVector)
{
    for (size_t i = 0; i < aVector.size(); i++)
    {
        ASSERT_EQ(aVector[i], bVector[i]);
    }
}

void
compare(const int32_t* aVector, const int32_t* bVector, const int32_t nElements)
{
    for (int32_t i = 0; i < nElements; i++)
    {
        ASSERT_EQ(aVector[i], bVector[i]);
    }
}

void
compare(const double* aVector, const double* bVector, const int32_t nElements)
{
    for (int32_t i = 0; i < nElements; i++)
    {
        ASSERT_NEAR(aVector[i], bVector[i], 1.0e-13);
    }
}

void
compare(const std::vector<double>& aVector, const std::vector<double>& bVector)
{
    ASSERT_EQ(aVector.size(), bVector.size());

    if (aVector.size() == bVector.size())
    {
        for (size_t i = 0; i < aVector.size(); i++)
        {
            ASSERT_NEAR(aVector[i], bVector[i], 1.0e-13);
        }
    }
}

void
compare(const std::vector<int32_t>& aVector, const std::vector<int32_t>& bVector)
{
    ASSERT_EQ(aVector.size(), bVector.size());

    if (aVector.size() == bVector.size())
    {
        for (size_t i = 0; i < aVector.size(); i++)
        {
            ASSERT_EQ(aVector[i], bVector[i]);
        }
    }
}

void
checkNorm(const double* aVector, const int32_t nElements)
{
    double fsum = 0.0;

    for (int32_t i = 0; i < nElements; i++)
        fsum += aVector[i];

    ASSERT_NEAR(fsum, 1.0, 1.0e-13);
}

}  // namespace vlxtest
