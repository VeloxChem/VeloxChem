//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MemBlock2DTest.hpp"

#include "MatOrder.hpp"
#include "MemBlock2D.hpp"

TEST_F(CMemBlock2DTest, DefaultConstructor)
{
    CMemBlock2D<double> ma;

    CMemBlock2D<double> mb(0, 0);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ConstructorWithNumberOfElements)
{
    CMemBlock2D<double> ma(3, 2);

    ma.zero();

    CMemBlock2D<double> mb({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 3, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ConstructorWithVector)
{
    CMemBlock2D<double> ma(std::vector<int32_t>({3, 3}));

    ma.zero();

    CMemBlock2D<double> mb({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 3, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ConstructorWithDataVector)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, {2, 2});

    CMemBlock2D<double> mb({1.0, 2.0, 3.0, 6.0}, 2, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ConstructorWithDataVectorAndOrder)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, matorder::row_major, 2, 2);

    CMemBlock2D<double> mb({1.0, 2.0, 3.0, 6.0}, 2, 2);

    ASSERT_EQ(ma, mb);

    CMemBlock2D<double> mc({1.0, 2.0, 3.0, 6.0}, matorder::col_major, 2, 2);

    CMemBlock2D<double> md({1.0, 3.0, 2.0, 6.0}, 2, 2);

    ASSERT_EQ(mc, md);
}

TEST_F(CMemBlock2DTest, CopyConstructor)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    CMemBlock2D<double> mb(ma);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, MoveConstructor)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    CMemBlock2D<double> mb(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0}, 2, 2));

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, CopyAssignment)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    CMemBlock2D<double> mb = ma;

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, MoveAssignment)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    CMemBlock2D<double> mb = CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0}, 2, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, Slice)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    CMemBlock2D<double> mb({1.0, 3.0}, 1, 2);

    CMemBlock2D<double> mc({2.0, 6.0}, 1, 2);

    ASSERT_EQ(ma, ma.slice(0, 2));

    ASSERT_EQ(mb, ma.slice(0, 1));

    ASSERT_EQ(mc, ma.slice(1, 1));
}

TEST_F(CMemBlock2DTest, Zero)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    ma.zero();

    CMemBlock2D<double> mb = CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0}, 2, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, IsEmpty)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    ASSERT_FALSE(ma.isEmpty());

    CMemBlock2D<double> mb;

    ASSERT_TRUE(mb.isEmpty());
}

TEST_F(CMemBlock2DTest, Data)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    auto ar = ma.data(0);

    ASSERT_EQ(0u, ((size_t)ar) % VLX_ALIGN);

    ASSERT_NEAR(1.0, ar[0], 1.0e-13);

    ASSERT_NEAR(2.0, ar[1], 1.0e-13);

    ar[0] = 4.0;

    ar[1] = 5.0;

    auto br = ma.data(1);

    ASSERT_EQ(0u, ((size_t)br) % VLX_ALIGN);

    ASSERT_NEAR(3.0, br[0], 1.0e-13);

    ASSERT_NEAR(6.0, br[1], 1.0e-13);

    br[0] = 8.0;

    br[1] = 9.0;

    CMemBlock2D<double> mb = CMemBlock2D<double>({4.0, 5.0, 8.0, 9.0}, 2, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, DataConstant)
{
    const CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    auto ar = ma.data(0);

    ASSERT_EQ(0u, ((size_t)ar) % VLX_ALIGN);

    ASSERT_NEAR(1.0, ar[0], 1.0e-13);

    ASSERT_NEAR(2.0, ar[1], 1.0e-13);

    auto br = ma.data(1);

    ASSERT_EQ(0u, ((size_t)br) % VLX_ALIGN);

    ASSERT_NEAR(3.0, br[0], 1.0e-13);

    ASSERT_NEAR(6.0, br[1], 1.0e-13);
}

TEST_F(CMemBlock2DTest, DataAtElement)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    auto ar = ma.data(0, 1);

    ASSERT_NEAR(2.0, ar[0], 1.0e-13);

    ar[0] = 4.0;

    auto br = ma.data(1, 0);

    ASSERT_NEAR(3.0, br[0], 1.0e-13);

    ASSERT_NEAR(6.0, br[1], 1.0e-13);

    br[0] = 8.0;

    br[1] = 9.0;

    CMemBlock2D<double> mb = CMemBlock2D<double>({1.0, 4.0, 8.0, 9.0}, 2, 2);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, DataConstantAtElement)
{
    const CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0}, 2, 2);

    auto ar = ma.data(0, 0);

    ASSERT_NEAR(1.0, ar[0], 1.0e-13);

    ASSERT_NEAR(2.0, ar[1], 1.0e-13);

    auto br = ma.data(1, 1);

    ASSERT_NEAR(6.0, br[0], 1.0e-13);
}

TEST_F(CMemBlock2DTest, Size)
{
    CMemBlock2D<double> ma(std::vector<int32_t>({2, 4, 5}));

    ASSERT_EQ(ma.size(0), 2);

    ASSERT_EQ(ma.size(1), 4);

    ASSERT_EQ(ma.size(2), 5);
}

TEST_F(CMemBlock2DTest, pitched_size)
{
    CMemBlock2D<double> ma(std::vector<int32_t>({2, 9, 5}));
    
    ASSERT_EQ(ma.pitched_size(0), 8);
    
    ASSERT_EQ(ma.pitched_size(1), 16);
    
    ASSERT_EQ(ma.pitched_size(2), 8);
}

TEST_F(CMemBlock2DTest, Blocks)
{
    CMemBlock2D<double> ma(3, 2);

    ASSERT_EQ(ma.blocks(), 2);
}

TEST_F(CMemBlock2DTest, BroadcastIntegers)
{
    CMemBlock2D<int32_t> ma({1, 2, 3, 6, 5}, {2, 3});

    ma.broadcast(0, MPI_COMM_WORLD);

    CMemBlock2D<int32_t> mb({1, 2, 3, 6, 5}, {2, 3});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, BroadcastReals)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0, 5.0}, {2, 3});

    ma.broadcast(0, MPI_COMM_WORLD);

    CMemBlock2D<double> mb({1.0, 2.0, 3.0, 6.0, 5.0}, {2, 3});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, GatherIntegers)
{
    CMemBlock2D<int32_t> ma({1, 2, 3, 6, 5}, {2, 3});

    auto mb = ma.gather(0, 1, MPI_COMM_WORLD);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, GatherReals)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0, 5.0}, {2, 3});

    auto mb = ma.gather(0, 1, MPI_COMM_WORLD);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ScatterIntegers)
{
    CMemBlock2D<int32_t> ma({1, 2, 3, 6, 5}, {2, 3});

    ma.scatter(0, 1, MPI_COMM_WORLD);

    CMemBlock2D<int32_t> mb({1, 2, 3, 6, 5}, {2, 3});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlock2DTest, ScatterReals)
{
    CMemBlock2D<double> ma({1.0, 2.0, 3.0, 6.0, 5.0}, {2, 3});

    ma.scatter(0, 1, MPI_COMM_WORLD);

    CMemBlock2D<double> mb({1.0, 2.0, 3.0, 6.0, 5.0}, {2, 3});

    ASSERT_EQ(ma, mb);
}
