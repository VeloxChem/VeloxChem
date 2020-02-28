//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "MemBlockTest.hpp"

#include "MemBlock.hpp"

TEST_F(CMemBlockTest, DefaultConstructor)
{
    CMemBlock<double> ma;

    CMemBlock<double> mb(std::vector<double>{});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, ConstructorWithNumberOfElements)
{
    CMemBlock<double> ma(5);

    ma.zero();

    CMemBlock<double> mb({0.0, 0.0, 0.0, 0.0, 0.0});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, CopyConstructor)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    CMemBlock<double> mb(ma);

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, MoveConstructor)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    CMemBlock<double> mb(CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, CopyAssignment)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    CMemBlock<double> mb = ma;

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, MoveAssignment)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    CMemBlock<double> mb = CMemBlock<double>({1.0, 2.0, 3.0, 6.0});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, Data)
{
    CMemBlock<int32_t> ma({1, 2, 7});

    auto dptr = ma.data();

    ASSERT_EQ(0u, ((size_t)dptr % VLX_ALIGN));

    ASSERT_EQ(1, dptr[0]);

    ASSERT_EQ(2, dptr[1]);

    ASSERT_EQ(7, dptr[2]);

    dptr[1] = 8;

    ASSERT_EQ(ma, CMemBlock<int32_t>({1, 8, 7}));
}

TEST_F(CMemBlockTest, ConstData)
{
    const CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    auto dptr = ma.data();

    ASSERT_EQ(0u, ((size_t)dptr) % VLX_ALIGN);

    ASSERT_NEAR(1.0, dptr[0], 1.0e-13);

    ASSERT_NEAR(2.0, dptr[1], 1.0e-13);

    ASSERT_NEAR(3.0, dptr[2], 1.0e-13);

    ASSERT_NEAR(6.0, dptr[3], 1.0e-13);
}

TEST_F(CMemBlockTest, DataAtPosition)
{
    CMemBlock<int32_t> ma({1, 2, 7});

    auto dptr = ma.data(2);

    ASSERT_EQ(7, dptr[0]);

    dptr[0] = 8;

    ASSERT_EQ(ma, CMemBlock<int32_t>({1, 2, 8}));
}

TEST_F(CMemBlockTest, ConstDataAtPosition)
{
    const CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    auto dptr = ma.data(1);

    ASSERT_NEAR(2.0, dptr[0], 1.0e-13);

    ASSERT_NEAR(3.0, dptr[1], 1.0e-13);

    ASSERT_NEAR(6.0, dptr[2], 1.0e-13);
}

TEST_F(CMemBlockTest, Size)
{
    const CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ASSERT_EQ(4, ma.size());
}

TEST_F(CMemBlockTest, Zero)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ma.zero();

    ASSERT_EQ(ma, CMemBlock<double>({0.0, 0.0, 0.0, 0.0}));
}

TEST_F(CMemBlockTest, At)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ASSERT_NEAR(ma.at(0), 1.0, 1.0e-13);

    ASSERT_NEAR(ma.at(1), 2.0, 1.0e-13);

    ASSERT_NEAR(ma.at(2), 3.0, 1.0e-13);

    ASSERT_NEAR(ma.at(3), 6.0, 1.0e-13);

    ma.at(3) = 4.0;

    CMemBlock<double> mb({1.0, 2.0, 3.0, 4.0});

    ASSERT_EQ(ma, mb);
}

TEST_F(CMemBlockTest, AtConstant)
{
    const CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ASSERT_NEAR(ma.at(0), 1.0, 1.0e-13);

    ASSERT_NEAR(ma.at(1), 2.0, 1.0e-13);

    ASSERT_NEAR(ma.at(2), 3.0, 1.0e-13);

    ASSERT_NEAR(ma.at(3), 6.0, 1.0e-13);
}

TEST_F(CMemBlockTest, Pack)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    auto mb = ma.pack(2, 1);

    ASSERT_EQ(mb, CMemBlock<double>({1.0, 2.0}));

    auto mc = ma.pack(2, 2);

    ASSERT_EQ(mc, CMemBlock<double>({4.0, 8.0}));

    auto md = ma.pack(3, 1);

    ASSERT_EQ(md, CMemBlock<double>({1.0, 2.0, 3.0}));

    auto mf = ma.pack(1, 4);

    ASSERT_EQ(mf, CMemBlock<double>(std::vector<double>{12.0}));

    auto me = ma.pack(3, 2);

    ASSERT_EQ(me, CMemBlock<double>());
}

TEST_F(CMemBlockTest, Pick)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    auto mb = ma.pick(2, 2, 0);

    ASSERT_EQ(mb, CMemBlock<double>({1.0, 3.0}));

    auto mc = ma.pick(2, 2, 1);

    ASSERT_EQ(mc, CMemBlock<double>({2.0, 6.0}));

    auto md = ma.pick(3, 1, 0);

    ASSERT_EQ(md, CMemBlock<double>(std::vector<double>({1.0})));

    auto mf = ma.pick(4, 1, 1);

    ASSERT_EQ(mf, CMemBlock<double>(std::vector<double>({2.0})));

    auto me = ma.pick(3, 1, 2);

    ASSERT_EQ(me, CMemBlock<double>(std::vector<double>({3.0})));

    auto mg = ma.pick(4, 1, 3);

    ASSERT_EQ(mg, CMemBlock<double>(std::vector<double>({6.0})));

    auto mk = ma.pick(1, 4, 0);

    ASSERT_EQ(mk, CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));

    auto ml = ma.pick(5, 1, 0);

    ASSERT_EQ(ml, CMemBlock<double>());
}

TEST_F(CMemBlockTest, Shrink)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ma.shrink(4);

    ASSERT_EQ(ma, CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));

    ma.shrink(6);

    ASSERT_EQ(ma, CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));

    ma.shrink(3);

    ASSERT_EQ(ma, CMemBlock<double>({1.0, 2.0, 3.0}));

    ma.shrink(2);

    ASSERT_EQ(ma, CMemBlock<double>({1.0, 2.0}));
}

TEST_F(CMemBlockTest, Slice)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ASSERT_EQ(ma.slice(0, 5), CMemBlock<double>());

    ASSERT_EQ(ma.slice(2, 4), CMemBlock<double>());

    ASSERT_EQ(ma.slice(0, 4), CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));

    ASSERT_EQ(ma.slice(1, 3), CMemBlock<double>({2.0, 3.0, 6.0}));

    ASSERT_EQ(ma.slice(0, 3), CMemBlock<double>({1.0, 2.0, 3.0}));

    ASSERT_EQ(ma.slice(1, 2), CMemBlock<double>({2.0, 3.0}));

    ASSERT_EQ(ma.slice(2, 2), CMemBlock<double>({3.0, 6.0}));
}

TEST_F(CMemBlockTest, BroadcastIntegers)
{
    CMemBlock<int32_t> ma({1, 2, 3, 9});

    ma.broadcast(0, MPI_COMM_WORLD);

    EXPECT_EQ(ma, CMemBlock<int32_t>({1, 2, 3, 9}));
}

TEST_F(CMemBlockTest, BroadcastReals)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ma.broadcast(0, MPI_COMM_WORLD);

    EXPECT_EQ(ma, CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));
}

TEST_F(CMemBlockTest, GatherIntegers)
{
    CMemBlock<int32_t> ma({1, 2, 3, 9});

    auto mb = ma.gather(0, 1, MPI_COMM_WORLD);

    EXPECT_EQ(mb, ma);
}

TEST_F(CMemBlockTest, GatherReals)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    auto mb = ma.gather(0, 1, MPI_COMM_WORLD);

    EXPECT_EQ(mb, ma);
}

TEST_F(CMemBlockTest, ScatterIntegers)
{
    CMemBlock<int32_t> ma({1, 2, 3, 9});

    ma.scatter(0, 1, MPI_COMM_WORLD);

    EXPECT_EQ(ma, CMemBlock<int32_t>({1, 2, 3, 9}));
}

TEST_F(CMemBlockTest, ScatterReals)
{
    CMemBlock<double> ma({1.0, 2.0, 3.0, 6.0});

    ma.scatter(0, 1, MPI_COMM_WORLD);

    EXPECT_EQ(ma, CMemBlock<double>({1.0, 2.0, 3.0, 6.0}));
}
