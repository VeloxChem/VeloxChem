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

#include "MpiFuncTest.hpp"

#include <mpi.h>

#include "MpiFunc.hpp"

TEST_F(CMpiFuncTest, Master)
{
    ASSERT_EQ(0, mpi::master());
}

TEST_F(CMpiFuncTest, Rank)
{
    ASSERT_EQ(0, mpi::rank(MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, Nodes)
{
    ASSERT_EQ(1, mpi::nodes(MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, Compare)
{
    ASSERT_TRUE(mpi::compare(MPI_COMM_WORLD, MPI_COMM_WORLD));
}

TEST_F(CMpiFuncTest, BcastInteger)
{
    int32_t mvalue = 3;

    mpi::bcast(mvalue, MPI_COMM_WORLD);

    ASSERT_EQ(3, mvalue);
}

TEST_F(CMpiFuncTest, BcastDouble)
{
    double mvalue = 3.4;

    mpi::bcast(mvalue, MPI_COMM_WORLD);

    ASSERT_NEAR(3.4, mvalue, 1.0e-13);
}

TEST_F(CMpiFuncTest, BcastBoolean)
{
    bool mvalue = false;

    mpi::bcast(mvalue, mpi::master(), MPI_COMM_WORLD);

    ASSERT_FALSE(mvalue);
}

TEST_F(CMpiFuncTest, BcastSymbol)
{
    char mvalue = 'x';

    mpi::bcast(mvalue, MPI_COMM_WORLD);

    ASSERT_EQ(mvalue, 'x');
}

TEST_F(CMpiFuncTest, BcastIntegersVector)
{
    std::vector<int32_t> mvector({2, 3, 4, -7, 9, 18});

    mpi::bcast(mvector, mpi::master(), MPI_COMM_WORLD);

    ASSERT_EQ(mvector[0], 2);

    ASSERT_EQ(mvector[1], 3);

    ASSERT_EQ(mvector[2], 4);

    ASSERT_EQ(mvector[3], -7);

    ASSERT_EQ(mvector[4], 9);

    ASSERT_EQ(mvector[5], 18);
}

TEST_F(CMpiFuncTest, BcastString)
{
    std::string str("Velox");

    mpi::bcast(str, mpi::master(), MPI_COMM_WORLD);

    ASSERT_EQ(str, std::string("Velox"));
}

TEST_F(CMpiFuncTest, Batch_size)
{
    ASSERT_EQ(6, mpi::batch_size(11, 0, 2));

    ASSERT_EQ(5, mpi::batch_size(11, 1, 2));
}

TEST_F(CMpiFuncTest, Batch_offset)
{
    ASSERT_EQ(0, mpi::batch_offset(11, 0, 2));

    ASSERT_EQ(6, mpi::batch_offset(11, 1, 2));
}

TEST_F(CMpiFuncTest, Batches_pattern)
{
    int32_t mpat[2];

    mpi::batches_pattern(mpat, 11, 2);

    ASSERT_EQ(mpat[0], 6);

    ASSERT_EQ(mpat[1], 5);
}

TEST_F(CMpiFuncTest, GatherInteger)
{
    int32_t mvalue = 0;

    mpi::gather(&mvalue, 9, mpi::master(), MPI_COMM_WORLD);

    ASSERT_EQ(mvalue, 9);
}
