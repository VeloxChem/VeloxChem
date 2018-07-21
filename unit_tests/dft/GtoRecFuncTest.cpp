//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoRecFuncTest.hpp"

#include "GtoRecFunc.hpp"
#include "MemBlock2D.hpp"

TEST_F(CGtoRecFuncTest, CompGtoTypePForLDA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 4);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 4);

    buffb.zero();

    auto pa = buffa.data(0);

    auto pb = buffb.data(0);

    pa[0] = 0.50; pb[0] = 0.50;

    pa[1] = 0.80; pb[1] = 0.80;

    pa[2] = 0.50; pb[2] = 0.50;

    pa[3] = 3.80; pb[3] = 3.80;

    pa[4] = 1.70; pb[4] = 1.70;

    gtorec::compGtoTypePForLDA(buffa, dist, ridx, 1);

    pb = buffb.data(1);

    pb[0] = 0.50;

    pb[1] = 0.96;

    pb[2] = 1.50;

    pb[3] = 15.20;

    pb = buffb.data(2);

    pb[0] = 1.00;

    pb[1] = 2.56;

    pb[2] = 0.40;

    pb[3] = 5.70;

    pb = buffb.data(3);

    pb[0] = 3.60;

    pb[1] = 6.24;

    pb[2] = 0.45;

    pb[3] = 7.22;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeDForLDA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 10);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 10);

    buffb.zero();

    auto pa = buffa.data(1);

    auto pb = buffb.data(1);

    pa[0] = 1.70; pb[0] = 1.70;

    pa[1] = 0.50; pb[1] = 0.50;

    pa[2] = 3.40; pb[2] = 3.40;

    pa[3] = 0.80; pb[3] = 0.80;

    pa[4] = 5.00; pb[4] = 5.00;

    pa = buffa.data(2);

    pb = buffb.data(2);

    pa[0] = 4.50; pb[0] = 4.50;

    pa[1] = 0.50; pb[1] = 0.50;

    pa[2] = 1.40; pb[2] = 1.40;

    pa[3] = 0.50; pb[3] = 0.50;

    pa[4] = 1.70; pb[4] = 1.70;

    pa = buffa.data(3);

    pb = buffb.data(3);

    pa[0] = 3.40; pb[0] = 3.40;

    pa[1] = 3.40; pb[1] = 3.40;

    pa[2] = 2.10; pb[2] = 2.10;

    pa[3] = 1.00; pb[3] = 1.00;

    pa[4] = 1.00; pb[4] = 1.00;

    gtorec::compGtoTypeDForLDA(buffa, dist, ridx, 1);

    pb = buffb.data(4);

    pb[0] = 1.70;

    pb[1] = 0.60;

    pb[2] = 10.20;

    pb[3] = 3.20;

    pb = buffb.data(5);

    pb[0] = 4.50;

    pb[1] = 0.60;

    pb[2] = 4.20;

    pb[3] = 2.00;

    pb = buffb.data(6);

    pb[0] = 3.40;

    pb[1] = 4.08;

    pb[2] = 6.30;

    pb[3] = 4.00;

    pb = buffb.data(7);

    pb[0] = 9.00;

    pb[1] = 1.60;

    pb[2] = 1.12;

    pb[3] = 0.75;

    pb = buffb.data(8);

    pb[0] = 6.80;

    pb[1] = 10.88;

    pb[2] = 1.68;

    pb[3] = 1.50;

    pb = buffb.data(9);

    pb[0] = 24.48;

    pb[1] = 26.52;

    pb[2] = 1.89;

    pb[3] = 1.90;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeFForLDA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 20);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 20);

    buffb.zero();

    auto pa = buffa.data(4);

    auto pb = buffb.data(4);

    pa[0] = 1.70; pb[0] = 1.70;

    pa[1] = 0.80; pb[1] = 0.80;

    pa[2] = 5.00; pb[2] = 5.00;

    pa[3] = 2.20; pb[3] = 2.20;

    pa[4] = 3.80; pb[4] = 3.80;

    pa = buffa.data(5);

    pb = buffb.data(5);

    pa[0] = 3.40; pb[0] = 3.40;

    pa[1] = 1.40; pb[1] = 1.40;

    pa[2] = 3.40; pb[2] = 3.40;

    pa[3] = 5.00; pb[3] = 5.00;

    pa[4] = 0.50; pb[4] = 0.50;

    pa = buffa.data(6);

    pb = buffb.data(6);

    pa[0] = 1.00; pb[0] = 1.00;

    pa[1] = 1.70; pb[1] = 1.70;

    pa[2] = 3.80; pb[2] = 3.80;

    pa[3] = 4.50; pb[3] = 4.50;

    pa[4] = 1.00; pb[4] = 1.00;

    pa = buffa.data(7);

    pb = buffb.data(7);

    pa[0] = 1.00; pb[0] = 1.00;

    pa[1] = 4.50; pb[1] = 4.50;

    pa[2] = 5.00; pb[2] = 5.00;

    pa[3] = 0.50; pb[3] = 0.50;

    pa[4] = 2.10; pb[4] = 2.10;

    pa = buffa.data(8);

    pb = buffb.data(8);

    pa[0] = 0.80; pb[0] = 0.80;

    pa[1] = 2.20; pb[1] = 2.20;

    pa[2] = 4.50; pb[2] = 4.50;

    pa[3] = 4.50; pb[3] = 4.50;

    pa[4] = 2.20; pb[4] = 2.20;

    pa = buffa.data(9);

    pb = buffb.data(9);

    pa[0] = 4.50; pb[0] = 4.50;

    pa[1] = 2.20; pb[1] = 2.20;

    pa[2] = 0.50; pb[2] = 0.50;

    pa[3] = 2.10; pb[3] = 2.10;

    pa[4] = 4.50; pb[4] = 4.50;

    gtorec::compGtoTypeFForLDA(buffa, dist, ridx, 1);

    pb = buffb.data(10);

    pb[0] = 1.70;

    pb[1] = 0.96;

    pb[2] = 15.00;

    pb[3] = 8.80;

    pb = buffb.data(11);

    pb[0] = 3.40;

    pb[1] = 1.68;

    pb[2] = 10.20;

    pb[3] = 20.00;

    pb = buffb.data(12);

    pb[0] = 1.00;

    pb[1] = 2.04;

    pb[2] = 11.40;

    pb[3] = 18.00;

    pb = buffb.data(13);

    pb[0] = 1.00;

    pb[1] = 5.40;

    pb[2] = 15.00;

    pb[3] = 2.00;

    pb = buffb.data(14);

    pb[0] = 0.80;

    pb[1] = 2.64;

    pb[2] = 13.50;

    pb[3] = 18.00;

    pb = buffb.data(15);

    pb[0] = 4.50;

    pb[1] = 2.64;

    pb[2] = 1.50;

    pb[3] = 8.40;

    pb = buffb.data(16);

    pb[0] = 2.00;

    pb[1] = 14.40;

    pb[2] = 4.00;

    pb[3] = 0.75;

    pb = buffb.data(17);

    pb[0] = 1.60;

    pb[1] = 7.04;

    pb[2] = 3.60;

    pb[3] = 6.75;

    pb = buffb.data(18);

    pb[0] = 9.00;

    pb[1] = 7.04;

    pb[2] = 0.40;

    pb[3] = 3.15;

    pb = buffb.data(19);

    pb[0] = 32.40;

    pb[1] = 17.16;

    pb[2] = 0.45;

    pb[3] = 3.99;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeGForLDA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 35);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 35);

    buffb.zero();

    auto pa = buffa.data(10);

    auto pb = buffb.data(10);

    pa[0] = 5.00; pb[0] = 5.00;

    pa[1] = 1.40; pb[1] = 1.40;

    pa[2] = 2.10; pb[2] = 2.10;

    pa[3] = 1.00; pb[3] = 1.00;

    pa[4] = 2.10; pb[4] = 2.10;

    pa = buffa.data(11);

    pb = buffb.data(11);

    pa[0] = 2.10; pb[0] = 2.10;

    pa[1] = 3.80; pb[1] = 3.80;

    pa[2] = 2.20; pb[2] = 2.20;

    pa[3] = 5.00; pb[3] = 5.00;

    pa[4] = 2.20; pb[4] = 2.20;

    pa = buffa.data(12);

    pb = buffb.data(12);

    pa[0] = 3.40; pb[0] = 3.40;

    pa[1] = 1.70; pb[1] = 1.70;

    pa[2] = 1.00; pb[2] = 1.00;

    pa[3] = 4.50; pb[3] = 4.50;

    pa[4] = 5.00; pb[4] = 5.00;

    pa = buffa.data(13);

    pb = buffb.data(13);

    pa[0] = 1.00; pb[0] = 1.00;

    pa[1] = 1.40; pb[1] = 1.40;

    pa[2] = 1.00; pb[2] = 1.00;

    pa[3] = 2.20; pb[3] = 2.20;

    pa[4] = 1.00; pb[4] = 1.00;

    pa = buffa.data(14);

    pb = buffb.data(14);

    pa[0] = 0.80; pb[0] = 0.80;

    pa[1] = 0.80; pb[1] = 0.80;

    pa[2] = 1.40; pb[2] = 1.40;

    pa[3] = 1.70; pb[3] = 1.70;

    pa[4] = 2.20; pb[4] = 2.20;

    pa = buffa.data(15);

    pb = buffb.data(15);

    pa[0] = 2.10; pb[0] = 2.10;

    pa[1] = 0.50; pb[1] = 0.50;

    pa[2] = 2.20; pb[2] = 2.20;

    pa[3] = 0.50; pb[3] = 0.50;

    pa[4] = 2.20; pb[4] = 2.20;

    pa = buffa.data(16);

    pb = buffb.data(16);

    pa[0] = 0.50; pb[0] = 0.50;

    pa[1] = 3.80; pb[1] = 3.80;

    pa[2] = 5.00; pb[2] = 5.00;

    pa[3] = 5.00; pb[3] = 5.00;

    pa[4] = 1.40; pb[4] = 1.40;

    pa = buffa.data(17);

    pb = buffb.data(17);

    pa[0] = 5.00; pb[0] = 5.00;

    pa[1] = 0.50; pb[1] = 0.50;

    pa[2] = 2.20; pb[2] = 2.20;

    pa[3] = 5.00; pb[3] = 5.00;

    pa[4] = 2.10; pb[4] = 2.10;

    pa = buffa.data(18);

    pb = buffb.data(18);

    pa[0] = 3.80; pb[0] = 3.80;

    pa[1] = 1.40; pb[1] = 1.40;

    pa[2] = 0.50; pb[2] = 0.50;

    pa[3] = 3.80; pb[3] = 3.80;

    pa[4] = 1.00; pb[4] = 1.00;

    pa = buffa.data(19);

    pb = buffb.data(19);

    pa[0] = 4.50; pb[0] = 4.50;

    pa[1] = 5.00; pb[1] = 5.00;

    pa[2] = 1.00; pb[2] = 1.00;

    pa[3] = 0.80; pb[3] = 0.80;

    pa[4] = 0.80; pb[4] = 0.80;

    gtorec::compGtoTypeGForLDA(buffa, dist, ridx, 1);

    pb = buffb.data(20);

    pb[0] = 5.00;

    pb[1] = 1.68;

    pb[2] = 6.30;

    pb[3] = 4.00;

    pb = buffb.data(21);

    pb[0] = 2.10;

    pb[1] = 4.56;

    pb[2] = 6.60;

    pb[3] = 20.00;

    pb = buffb.data(22);

    pb[0] = 3.40;

    pb[1] = 2.04;

    pb[2] = 3.00;

    pb[3] = 18.00;

    pb = buffb.data(23);

    pb[0] = 1.00;

    pb[1] = 1.68;

    pb[2] = 3.00;

    pb[3] = 8.80;

    pb = buffb.data(24);

    pb[0] = 0.80;

    pb[1] = 0.96;

    pb[2] = 4.20;

    pb[3] = 6.80;

    pb = buffb.data(25);

    pb[0] = 2.10;

    pb[1] = 0.60;

    pb[2] = 6.60;

    pb[3] = 2.00;

    pb = buffb.data(26);

    pb[0] = 0.50;

    pb[1] = 4.56;

    pb[2] = 15.00;

    pb[3] = 20.00;

    pb = buffb.data(27);

    pb[0] = 5.00;

    pb[1] = 0.60;

    pb[2] = 6.60;

    pb[3] = 20.00;

    pb = buffb.data(28);

    pb[0] = 3.80;

    pb[1] = 1.68;

    pb[2] = 1.50;

    pb[3] = 15.20;

    pb = buffb.data(29);

    pb[0] = 4.50;

    pb[1] = 6.00;

    pb[2] = 3.00;

    pb[3] = 3.20;

    pb = buffb.data(30);

    pb[0] = 1.00;

    pb[1] = 12.16;

    pb[2] = 4.00;

    pb[3] = 7.50;

    pb = buffb.data(31);

    pb[0] = 10.00;

    pb[1] = 1.60;

    pb[2] = 1.76;

    pb[3] = 7.50;

    pb = buffb.data(32);

    pb[0] = 7.60;

    pb[1] = 4.48;

    pb[2] = 0.40;

    pb[3] = 5.70;

    pb = buffb.data(33);

    pb[0] = 9.00;

    pb[1] = 16.00;

    pb[2] = 0.80;

    pb[3] = 1.20;

    pb = buffb.data(34);

    pb[0] = 32.40;

    pb[1] = 39.00;

    pb[2] = 0.90;

    pb[3] = 1.52;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypePForGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 16);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 16);

    buffb.zero();

    auto pa_0 = buffa.data(0);

    auto pb_0 = buffb.data(0);

    auto pa_x = buffa.data(1);

    auto pb_x = buffb.data(1);

    auto pa_y = buffa.data(2);

    auto pb_y = buffb.data(2);

    auto pa_z = buffa.data(3);

    auto pb_z = buffb.data(3);

    pa_0[0] = 4.50; pb_0[0] = 4.50;

    pa_x[0] = 1.70; pb_x[0] = 1.70;

    pa_y[0] = 0.80; pb_y[0] = 0.80;

    pa_z[0] = 2.10; pb_z[0] = 2.10;

    pa_0[1] = 1.70; pb_0[1] = 1.70;

    pa_x[1] = 0.80; pb_x[1] = 0.80;

    pa_y[1] = 2.10; pb_y[1] = 2.10;

    pa_z[1] = 1.70; pb_z[1] = 1.70;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 2.10; pb_z[2] = 2.10;

    pa_0[3] = 4.50; pb_0[3] = 4.50;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 1.40; pb_z[3] = 1.40;

    pa_0[4] = 4.50; pb_0[4] = 4.50;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 3.40; pb_y[4] = 3.40;

    pa_z[4] = 1.40; pb_z[4] = 1.40;

    gtorec::compGtoTypePForGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(4);

    pb_x = buffb.data(5);

    pb_y = buffb.data(6);

    pb_z = buffb.data(7);

    pb_0[0] = 4.50;

    pb_x[0] = 6.20;

    pb_y[0] = 0.80;

    pb_z[0] = 2.10;

    pb_0[1] = 2.04;

    pb_x[1] = 2.66;

    pb_y[1] = 2.52;

    pb_z[1] = 2.04;

    pb_0[2] = 11.40;

    pb_x[2] = 6.80;

    pb_y[2] = 10.20;

    pb_z[2] = 6.30;

    pb_0[3] = 18.00;

    pb_x[3] = 22.50;

    pb_y[3] = 2.00;

    pb_z[3] = 5.60;

    pb_0 = buffb.data(8);

    pb_x = buffb.data(9);

    pb_y = buffb.data(10);

    pb_z = buffb.data(11);

    pb_0[0] = 9.00;

    pb_x[0] = 3.40;

    pb_y[0] = 6.10;

    pb_z[0] = 4.20;

    pb_0[1] = 5.44;

    pb_x[1] = 2.56;

    pb_y[1] = 8.42;

    pb_z[1] = 5.44;

    pb_0[2] = 3.04;

    pb_x[2] = 0.80;

    pb_y[2] = 6.52;

    pb_z[2] = 1.68;

    pb_0[3] = 6.75;

    pb_x[3] = 6.75;

    pb_y[3] = 5.25;

    pb_z[3] = 2.10;

    pb_0 = buffb.data(12);

    pb_x = buffb.data(13);

    pb_y = buffb.data(14);

    pb_z = buffb.data(15);

    pb_0[0] = 32.40;

    pb_x[0] = 12.24;

    pb_y[0] = 5.76;

    pb_z[0] = 19.62;

    pb_0[1] = 13.26;

    pb_x[1] = 6.24;

    pb_y[1] = 16.38;

    pb_z[1] = 14.96;

    pb_0[2] = 3.42;

    pb_x[2] = 0.90;

    pb_y[2] = 3.06;

    pb_z[2] = 5.69;

    pb_0[3] = 8.55;

    pb_x[3] = 8.55;

    pb_y[3] = 0.95;

    pb_z[3] = 7.16;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeDForGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 40);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 40);

    buffb.zero();

    auto pa_0 = buffa.data(4);

    auto pb_0 = buffb.data(4);

    auto pa_x = buffa.data(5);

    auto pb_x = buffb.data(5);

    auto pa_y = buffa.data(6);

    auto pb_y = buffb.data(6);

    auto pa_z = buffa.data(7);

    auto pb_z = buffb.data(7);

    pa_0[0] = 4.50; pb_0[0] = 4.50;

    pa_x[0] = 3.80; pb_x[0] = 3.80;

    pa_y[0] = 2.20; pb_y[0] = 2.20;

    pa_z[0] = 3.40; pb_z[0] = 3.40;

    pa_0[1] = 3.80; pb_0[1] = 3.80;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 4.50; pb_z[1] = 4.50;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 1.00; pb_y[2] = 1.00;

    pa_z[2] = 1.00; pb_z[2] = 1.00;

    pa_0[3] = 3.80; pb_0[3] = 3.80;

    pa_x[3] = 2.10; pb_x[3] = 2.10;

    pa_y[3] = 5.00; pb_y[3] = 5.00;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 0.80; pb_x[4] = 0.80;

    pa_y[4] = 1.00; pb_y[4] = 1.00;

    pa_z[4] = 3.40; pb_z[4] = 3.40;

    pa_0 = buffa.data(8);

    pb_0 = buffb.data(8);

    pa_x = buffa.data(9);

    pb_x = buffb.data(9);

    pa_y = buffa.data(10);

    pb_y = buffb.data(10);

    pa_z = buffa.data(11);

    pb_z = buffb.data(11);

    pa_0[0] = 0.50; pb_0[0] = 0.50;

    pa_x[0] = 0.50; pb_x[0] = 0.50;

    pa_y[0] = 1.40; pb_y[0] = 1.40;

    pa_z[0] = 2.20; pb_z[0] = 2.20;

    pa_0[1] = 0.50; pb_0[1] = 0.50;

    pa_x[1] = 3.80; pb_x[1] = 3.80;

    pa_y[1] = 3.80; pb_y[1] = 3.80;

    pa_z[1] = 0.50; pb_z[1] = 0.50;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 2.20; pb_x[2] = 2.20;

    pa_y[2] = 3.80; pb_y[2] = 3.80;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 1.40; pb_x[3] = 1.40;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_0[4] = 3.80; pb_0[4] = 3.80;

    pa_x[4] = 0.80; pb_x[4] = 0.80;

    pa_y[4] = 5.00; pb_y[4] = 5.00;

    pa_z[4] = 3.40; pb_z[4] = 3.40;

    pa_0 = buffa.data(12);

    pb_0 = buffb.data(12);

    pa_x = buffa.data(13);

    pb_x = buffb.data(13);

    pa_y = buffa.data(14);

    pb_y = buffb.data(14);

    pa_z = buffa.data(15);

    pb_z = buffb.data(15);

    pa_0[0] = 0.80; pb_0[0] = 0.80;

    pa_x[0] = 0.80; pb_x[0] = 0.80;

    pa_y[0] = 4.50; pb_y[0] = 4.50;

    pa_z[0] = 4.50; pb_z[0] = 4.50;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 2.10; pb_x[1] = 2.10;

    pa_y[1] = 0.80; pb_y[1] = 0.80;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_0[2] = 4.50; pb_0[2] = 4.50;

    pa_x[2] = 2.20; pb_x[2] = 2.20;

    pa_y[2] = 4.50; pb_y[2] = 4.50;

    pa_z[2] = 0.80; pb_z[2] = 0.80;

    pa_0[3] = 0.80; pb_0[3] = 0.80;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 4.50; pb_z[3] = 4.50;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    gtorec::compGtoTypeDForGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(16);

    pb_x = buffb.data(17);

    pb_y = buffb.data(18);

    pb_z = buffb.data(19);

    pb_0[0] = 4.50;

    pb_x[0] = 8.30;

    pb_y[0] = 2.20;

    pb_z[0] = 3.40;

    pb_0[1] = 4.56;

    pb_x[1] = 5.84;

    pb_y[1] = 2.64;

    pb_z[1] = 5.40;

    pb_0[2] = 11.40;

    pb_x[2] = 6.80;

    pb_y[2] = 3.00;

    pb_z[2] = 3.00;

    pb_0[3] = 15.20;

    pb_x[3] = 12.20;

    pb_y[3] = 20.00;

    pb_z[3] = 2.00;

    pb_0 = buffb.data(20);

    pb_x = buffb.data(21);

    pb_y = buffb.data(22);

    pb_z = buffb.data(23);

    pb_0[0] = 9.00;

    pb_x[0] = 7.60;

    pb_y[0] = 8.90;

    pb_z[0] = 6.80;

    pb_0[1] = 12.16;

    pb_x[1] = 5.44;

    pb_y[1] = 10.84;

    pb_z[1] = 14.40;

    pb_0[2] = 3.04;

    pb_x[2] = 0.80;

    pb_y[2] = 4.60;

    pb_z[2] = 0.80;

    pb_0[3] = 5.70;

    pb_x[3] = 3.15;

    pb_y[3] = 11.30;

    pb_z[3] = 0.75;

    pb_0 = buffb.data(24);

    pb_x = buffb.data(25);

    pb_y = buffb.data(26);

    pb_z = buffb.data(27);

    pb_0[0] = 32.40;

    pb_x[0] = 27.36;

    pb_y[0] = 15.84;

    pb_z[0] = 28.98;

    pb_0[1] = 29.64;

    pb_x[1] = 13.26;

    pb_y[1] = 17.16;

    pb_z[1] = 38.90;

    pb_0[2] = 3.42;

    pb_x[2] = 0.90;

    pb_y[2] = 0.90;

    pb_z[2] = 4.70;

    pb_0[3] = 7.22;

    pb_x[3] = 3.99;

    pb_y[3] = 9.50;

    pb_z[3] = 4.75;

    pb_0 = buffb.data(28);

    pb_x = buffb.data(29);

    pb_y = buffb.data(30);

    pb_z = buffb.data(31);

    pb_0[0] = 1.00;

    pb_x[0] = 1.00;

    pb_y[0] = 3.30;

    pb_z[0] = 4.40;

    pb_0[1] = 1.60;

    pb_x[1] = 12.16;

    pb_y[1] = 12.66;

    pb_z[1] = 1.60;

    pb_0[2] = 3.04;

    pb_x[2] = 1.76;

    pb_y[2] = 6.84;

    pb_z[2] = 1.76;

    pb_0[3] = 7.50;

    pb_x[3] = 2.10;

    pb_y[3] = 10.10;

    pb_z[3] = 0.75;

    pb_0 = buffb.data(32);

    pb_x = buffb.data(33);

    pb_y = buffb.data(34);

    pb_z = buffb.data(35);

    pb_0[0] = 3.60;

    pb_x[0] = 3.60;

    pb_y[0] = 10.08;

    pb_z[0] = 16.34;

    pb_0[1] = 3.90;

    pb_x[1] = 29.64;

    pb_y[1] = 29.64;

    pb_z[1] = 4.40;

    pb_0[2] = 3.42;

    pb_x[2] = 1.98;

    pb_y[2] = 3.42;

    pb_z[2] = 5.78;

    pb_0[3] = 9.50;

    pb_x[3] = 2.66;

    pb_y[3] = 6.46;

    pb_z[3] = 5.95;

    pb_0 = buffb.data(36);

    pb_x = buffb.data(37);

    pb_y = buffb.data(38);

    pb_z = buffb.data(39);

    pb_0[0] = 5.76;

    pb_x[0] = 5.76;

    pb_y[0] = 32.40;

    pb_z[0] = 33.20;

    pb_0[1] = 16.38;

    pb_x[1] = 16.38;

    pb_y[1] = 6.24;

    pb_z[1] = 8.34;

    pb_0[2] = 4.05;

    pb_x[2] = 1.98;

    pb_y[2] = 4.05;

    pb_z[2] = 5.22;

    pb_0[3] = 1.52;

    pb_x[3] = 8.55;

    pb_y[3] = 6.46;

    pb_z[3] = 9.35;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeFForGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 80);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 80);

    buffb.zero();

    auto pa_0 = buffa.data(16);

    auto pb_0 = buffb.data(16);

    auto pa_x = buffa.data(17);

    auto pb_x = buffb.data(17);

    auto pa_y = buffa.data(18);

    auto pb_y = buffb.data(18);

    auto pa_z = buffa.data(19);

    auto pb_z = buffb.data(19);

    pa_0[0] = 1.70; pb_0[0] = 1.70;

    pa_x[0] = 5.00; pb_x[0] = 5.00;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 2.10; pb_z[0] = 2.10;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 4.50; pb_x[1] = 4.50;

    pa_y[1] = 0.80; pb_y[1] = 0.80;

    pa_z[1] = 0.50; pb_z[1] = 0.50;

    pa_0[2] = 1.00; pb_0[2] = 1.00;

    pa_x[2] = 2.10; pb_x[2] = 2.10;

    pa_y[2] = 5.00; pb_y[2] = 5.00;

    pa_z[2] = 3.80; pb_z[2] = 3.80;

    pa_0[3] = 1.40; pb_0[3] = 1.40;

    pa_x[3] = 3.80; pb_x[3] = 3.80;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_0[4] = 2.10; pb_0[4] = 2.10;

    pa_x[4] = 0.80; pb_x[4] = 0.80;

    pa_y[4] = 0.50; pb_y[4] = 0.50;

    pa_z[4] = 5.00; pb_z[4] = 5.00;

    pa_0 = buffa.data(20);

    pb_0 = buffb.data(20);

    pa_x = buffa.data(21);

    pb_x = buffb.data(21);

    pa_y = buffa.data(22);

    pb_y = buffb.data(22);

    pa_z = buffa.data(23);

    pb_z = buffb.data(23);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 0.80; pb_x[0] = 0.80;

    pa_y[0] = 1.40; pb_y[0] = 1.40;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 1.00; pb_z[1] = 1.00;

    pa_0[2] = 0.80; pb_0[2] = 0.80;

    pa_x[2] = 3.40; pb_x[2] = 3.40;

    pa_y[2] = 0.80; pb_y[2] = 0.80;

    pa_z[2] = 1.70; pb_z[2] = 1.70;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 0.80; pb_y[3] = 0.80;

    pa_z[3] = 1.00; pb_z[3] = 1.00;

    pa_0[4] = 1.40; pb_0[4] = 1.40;

    pa_x[4] = 0.50; pb_x[4] = 0.50;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    pa_0 = buffa.data(24);

    pb_0 = buffb.data(24);

    pa_x = buffa.data(25);

    pb_x = buffb.data(25);

    pa_y = buffa.data(26);

    pb_y = buffb.data(26);

    pa_z = buffa.data(27);

    pb_z = buffb.data(27);

    pa_0[0] = 2.20; pb_0[0] = 2.20;

    pa_x[0] = 1.40; pb_x[0] = 1.40;

    pa_y[0] = 0.50; pb_y[0] = 0.50;

    pa_z[0] = 1.00; pb_z[0] = 1.00;

    pa_0[1] = 4.50; pb_0[1] = 4.50;

    pa_x[1] = 1.40; pb_x[1] = 1.40;

    pa_y[1] = 0.80; pb_y[1] = 0.80;

    pa_z[1] = 1.00; pb_z[1] = 1.00;

    pa_0[2] = 3.40; pb_0[2] = 3.40;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 2.10; pb_y[2] = 2.10;

    pa_z[2] = 3.40; pb_z[2] = 3.40;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 1.70; pb_z[3] = 1.70;

    pa_0[4] = 1.40; pb_0[4] = 1.40;

    pa_x[4] = 0.50; pb_x[4] = 0.50;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 3.80; pb_z[4] = 3.80;

    pa_0 = buffa.data(28);

    pb_0 = buffb.data(28);

    pa_x = buffa.data(29);

    pb_x = buffb.data(29);

    pa_y = buffa.data(30);

    pb_y = buffb.data(30);

    pa_z = buffa.data(31);

    pb_z = buffb.data(31);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 2.10; pb_y[0] = 2.10;

    pa_z[0] = 1.70; pb_z[0] = 1.70;

    pa_0[1] = 0.80; pb_0[1] = 0.80;

    pa_x[1] = 3.40; pb_x[1] = 3.40;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 1.00; pb_z[1] = 1.00;

    pa_0[2] = 4.50; pb_0[2] = 4.50;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 5.00; pb_y[2] = 5.00;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 3.40; pb_0[3] = 3.40;

    pa_x[3] = 2.20; pb_x[3] = 2.20;

    pa_y[3] = 1.40; pb_y[3] = 1.40;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_0[4] = 4.50; pb_0[4] = 4.50;

    pa_x[4] = 1.70; pb_x[4] = 1.70;

    pa_y[4] = 4.50; pb_y[4] = 4.50;

    pa_z[4] = 4.50; pb_z[4] = 4.50;

    pa_0 = buffa.data(32);

    pb_0 = buffb.data(32);

    pa_x = buffa.data(33);

    pb_x = buffb.data(33);

    pa_y = buffa.data(34);

    pb_y = buffb.data(34);

    pa_z = buffa.data(35);

    pb_z = buffb.data(35);

    pa_0[0] = 0.80; pb_0[0] = 0.80;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 3.80; pb_z[0] = 3.80;

    pa_0[1] = 0.80; pb_0[1] = 0.80;

    pa_x[1] = 0.50; pb_x[1] = 0.50;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 3.40; pb_z[1] = 3.40;

    pa_0[2] = 4.50; pb_0[2] = 4.50;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 0.50; pb_y[2] = 0.50;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 0.50; pb_0[3] = 0.50;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 1.70; pb_y[3] = 1.70;

    pa_z[3] = 3.40; pb_z[3] = 3.40;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 1.70; pb_x[4] = 1.70;

    pa_y[4] = 2.20; pb_y[4] = 2.20;

    pa_z[4] = 1.00; pb_z[4] = 1.00;

    pa_0 = buffa.data(36);

    pb_0 = buffb.data(36);

    pa_x = buffa.data(37);

    pb_x = buffb.data(37);

    pa_y = buffa.data(38);

    pb_y = buffb.data(38);

    pa_z = buffa.data(39);

    pb_z = buffb.data(39);

    pa_0[0] = 1.70; pb_0[0] = 1.70;

    pa_x[0] = 2.10; pb_x[0] = 2.10;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 1.00; pb_z[0] = 1.00;

    pa_0[1] = 5.00; pb_0[1] = 5.00;

    pa_x[1] = 1.00; pb_x[1] = 1.00;

    pa_y[1] = 4.50; pb_y[1] = 4.50;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_0[2] = 1.00; pb_0[2] = 1.00;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 1.00; pb_y[2] = 1.00;

    pa_z[2] = 1.70; pb_z[2] = 1.70;

    pa_0[3] = 1.40; pb_0[3] = 1.40;

    pa_x[3] = 3.40; pb_x[3] = 3.40;

    pa_y[3] = 1.70; pb_y[3] = 1.70;

    pa_z[3] = 5.00; pb_z[3] = 5.00;

    pa_0[4] = 1.70; pb_0[4] = 1.70;

    pa_x[4] = 3.40; pb_x[4] = 3.40;

    pa_y[4] = 0.80; pb_y[4] = 0.80;

    pa_z[4] = 1.00; pb_z[4] = 1.00;

    gtorec::compGtoTypeFForGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(40);

    pb_x = buffb.data(41);

    pb_y = buffb.data(42);

    pb_z = buffb.data(43);

    pb_0[0] = 1.70;

    pb_x[0] = 6.70;

    pb_y[0] = 3.80;

    pb_z[0] = 2.10;

    pb_0[1] = 2.52;

    pb_x[1] = 7.50;

    pb_y[1] = 0.96;

    pb_z[1] = 0.60;

    pb_0[2] = 3.00;

    pb_x[2] = 7.30;

    pb_y[2] = 15.00;

    pb_z[2] = 11.40;

    pb_0[3] = 5.60;

    pb_x[3] = 16.60;

    pb_y[3] = 2.00;

    pb_z[3] = 8.80;

    pb_0 = buffb.data(44);

    pb_x = buffb.data(45);

    pb_y = buffb.data(46);

    pb_z = buffb.data(47);

    pb_0[0] = 3.40;

    pb_x[0] = 10.00;

    pb_y[0] = 9.30;

    pb_z[0] = 4.20;

    pb_0[1] = 6.72;

    pb_x[1] = 14.40;

    pb_y[1] = 4.66;

    pb_z[1] = 1.60;

    pb_0[2] = 0.80;

    pb_x[2] = 1.68;

    pb_y[2] = 5.00;

    pb_z[2] = 3.04;

    pb_0[3] = 2.10;

    pb_x[3] = 5.70;

    pb_y[3] = 2.15;

    pb_z[3] = 3.30;

    pb_0 = buffb.data(48);

    pb_x = buffb.data(49);

    pb_y = buffb.data(50);

    pb_z = buffb.data(51);

    pb_0[0] = 12.24;

    pb_x[0] = 36.00;

    pb_y[0] = 27.36;

    pb_z[0] = 16.82;

    pb_0[1] = 16.38;

    pb_x[1] = 35.10;

    pb_y[1] = 6.24;

    pb_z[1] = 6.00;

    pb_0[2] = 0.90;

    pb_x[2] = 1.89;

    pb_y[2] = 4.50;

    pb_z[2] = 4.42;

    pb_0[3] = 2.66;

    pb_x[3] = 7.22;

    pb_y[3] = 0.95;

    pb_z[3] = 5.58;

    pb_0 = buffb.data(52);

    pb_x = buffb.data(53);

    pb_y = buffb.data(54);

    pb_z = buffb.data(55);

    pb_0[0] = 3.40;

    pb_x[0] = 5.60;

    pb_y[0] = 2.10;

    pb_z[0] = 1.70;

    pb_0[1] = 0.96;

    pb_x[1] = 4.88;

    pb_y[1] = 2.64;

    pb_z[1] = 1.20;

    pb_0[2] = 13.50;

    pb_x[2] = 9.60;

    pb_y[2] = 15.00;

    pb_z[2] = 6.60;

    pb_0[3] = 13.60;

    pb_x[3] = 12.20;

    pb_y[3] = 5.60;

    pb_z[3] = 2.00;

    pb_0 = buffb.data(56);

    pb_x = buffb.data(57);

    pb_y = buffb.data(58);

    pb_z = buffb.data(59);

    pb_0[0] = 0.80;

    pb_x[0] = 3.00;

    pb_y[0] = 3.80;

    pb_z[0] = 3.80;

    pb_0[1] = 0.96;

    pb_x[1] = 1.40;

    pb_y[1] = 2.64;

    pb_z[1] = 4.08;

    pb_0[2] = 13.50;

    pb_x[2] = 7.50;

    pb_y[2] = 1.50;

    pb_z[2] = 6.60;

    pb_0[3] = 2.00;

    pb_x[3] = 4.50;

    pb_y[3] = 6.80;

    pb_z[3] = 13.60;

    pb_0 = buffb.data(60);

    pb_x = buffb.data(61);

    pb_y = buffb.data(62);

    pb_z = buffb.data(63);

    pb_0[0] = 1.70;

    pb_x[0] = 3.80;

    pb_y[0] = 3.80;

    pb_z[0] = 1.00;

    pb_0[1] = 6.00;

    pb_x[1] = 6.20;

    pb_y[1] = 5.40;

    pb_z[1] = 2.64;

    pb_0[2] = 3.00;

    pb_x[2] = 14.50;

    pb_y[2] = 3.00;

    pb_z[2] = 5.10;

    pb_0[3] = 5.60;

    pb_x[3] = 15.00;

    pb_y[3] = 6.80;

    pb_z[3] = 20.00;

    pb_0 = buffb.data(64);

    pb_x = buffb.data(65);

    pb_y = buffb.data(66);

    pb_z = buffb.data(67);

    pb_0[0] = 6.80;

    pb_x[0] = 4.40;

    pb_y[0] = 7.60;

    pb_z[0] = 3.40;

    pb_0[1] = 2.56;

    pb_x[1] = 10.88;

    pb_y[1] = 7.84;

    pb_z[1] = 3.20;

    pb_0[2] = 3.60;

    pb_x[2] = 1.36;

    pb_y[2] = 8.50;

    pb_z[2] = 1.76;

    pb_0[3] = 5.10;

    pb_x[3] = 3.30;

    pb_y[3] = 5.50;

    pb_z[3] = 0.75;

    pb_0 = buffb.data(68);

    pb_x = buffb.data(69);

    pb_y = buffb.data(70);

    pb_z = buffb.data(71);

    pb_0[0] = 24.48;

    pb_x[0] = 15.84;

    pb_y[0] = 15.12;

    pb_z[0] = 15.64;

    pb_0[1] = 6.24;

    pb_x[1] = 26.52;

    pb_y[1] = 17.16;

    pb_z[1] = 8.60;

    pb_0[2] = 4.05;

    pb_x[2] = 1.53;

    pb_y[2] = 4.50;

    pb_z[2] = 6.48;

    pb_0[3] = 6.46;

    pb_x[3] = 4.18;

    pb_y[3] = 2.66;

    pb_z[3] = 4.35;

    pb_0 = buffb.data(72);

    pb_x = buffb.data(73);

    pb_y = buffb.data(74);

    pb_z = buffb.data(75);

    pb_0[0] = 3.40;

    pb_x[0] = 4.20;

    pb_y[0] = 9.30;

    pb_z[0] = 2.00;

    pb_0[1] = 16.00;

    pb_x[1] = 3.20;

    pb_y[1] = 19.40;

    pb_z[1] = 7.04;

    pb_0[2] = 0.80;

    pb_x[2] = 3.60;

    pb_y[2] = 1.80;

    pb_z[2] = 1.36;

    pb_0[3] = 2.10;

    pb_x[3] = 5.10;

    pb_y[3] = 3.95;

    pb_z[3] = 7.50;

    pb_0 = buffb.data(76);

    pb_x = buffb.data(77);

    pb_y = buffb.data(78);

    pb_z = buffb.data(79);

    pb_0[0] = 12.24;

    pb_x[0] = 15.12;

    pb_y[0] = 27.36;

    pb_z[0] = 8.90;

    pb_0[1] = 39.00;

    pb_x[1] = 7.80;

    pb_y[1] = 35.10;

    pb_z[1] = 22.16;

    pb_0[2] = 0.90;

    pb_x[2] = 4.05;

    pb_y[2] = 0.90;

    pb_z[2] = 2.53;

    pb_0[3] = 2.66;

    pb_x[3] = 6.46;

    pb_y[3] = 3.23;

    pb_z[3] = 10.90;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeGForGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 140);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 140);

    buffb.zero();

    auto pa_0 = buffa.data(40);

    auto pb_0 = buffb.data(40);

    auto pa_x = buffa.data(41);

    auto pb_x = buffb.data(41);

    auto pa_y = buffa.data(42);

    auto pb_y = buffb.data(42);

    auto pa_z = buffa.data(43);

    auto pb_z = buffb.data(43);

    pa_0[0] = 1.00; pb_0[0] = 1.00;

    pa_x[0] = 2.10; pb_x[0] = 2.10;

    pa_y[0] = 0.50; pb_y[0] = 0.50;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 3.40; pb_x[1] = 3.40;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 3.40; pb_z[1] = 3.40;

    pa_0[2] = 2.20; pb_0[2] = 2.20;

    pa_x[2] = 0.80; pb_x[2] = 0.80;

    pa_y[2] = 1.70; pb_y[2] = 1.70;

    pa_z[2] = 0.50; pb_z[2] = 0.50;

    pa_0[3] = 1.00; pb_0[3] = 1.00;

    pa_x[3] = 3.80; pb_x[3] = 3.80;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 1.40; pb_z[3] = 1.40;

    pa_0[4] = 2.10; pb_0[4] = 2.10;

    pa_x[4] = 1.00; pb_x[4] = 1.00;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 0.50; pb_z[4] = 0.50;

    pa_0 = buffa.data(44);

    pb_0 = buffb.data(44);

    pa_x = buffa.data(45);

    pb_x = buffb.data(45);

    pa_y = buffa.data(46);

    pb_y = buffb.data(46);

    pa_z = buffa.data(47);

    pb_z = buffb.data(47);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 1.00; pb_x[0] = 1.00;

    pa_y[0] = 0.50; pb_y[0] = 0.50;

    pa_z[0] = 2.20; pb_z[0] = 2.20;

    pa_0[1] = 1.00; pb_0[1] = 1.00;

    pa_x[1] = 2.10; pb_x[1] = 2.10;

    pa_y[1] = 1.40; pb_y[1] = 1.40;

    pa_z[1] = 0.50; pb_z[1] = 0.50;

    pa_0[2] = 0.80; pb_0[2] = 0.80;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 4.50; pb_y[2] = 4.50;

    pa_z[2] = 3.40; pb_z[2] = 3.40;

    pa_0[3] = 1.40; pb_0[3] = 1.40;

    pa_x[3] = 1.70; pb_x[3] = 1.70;

    pa_y[3] = 1.00; pb_y[3] = 1.00;

    pa_z[3] = 1.70; pb_z[3] = 1.70;

    pa_0[4] = 1.70; pb_0[4] = 1.70;

    pa_x[4] = 5.00; pb_x[4] = 5.00;

    pa_y[4] = 3.40; pb_y[4] = 3.40;

    pa_z[4] = 5.00; pb_z[4] = 5.00;

    pa_0 = buffa.data(48);

    pb_0 = buffb.data(48);

    pa_x = buffa.data(49);

    pb_x = buffb.data(49);

    pa_y = buffa.data(50);

    pb_y = buffb.data(50);

    pa_z = buffa.data(51);

    pb_z = buffb.data(51);

    pa_0[0] = 5.00; pb_0[0] = 5.00;

    pa_x[0] = 1.70; pb_x[0] = 1.70;

    pa_y[0] = 0.50; pb_y[0] = 0.50;

    pa_z[0] = 5.00; pb_z[0] = 5.00;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 2.20; pb_x[1] = 2.20;

    pa_y[1] = 5.00; pb_y[1] = 5.00;

    pa_z[1] = 1.40; pb_z[1] = 1.40;

    pa_0[2] = 2.20; pb_0[2] = 2.20;

    pa_x[2] = 3.40; pb_x[2] = 3.40;

    pa_y[2] = 4.50; pb_y[2] = 4.50;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 1.40; pb_0[3] = 1.40;

    pa_x[3] = 2.20; pb_x[3] = 2.20;

    pa_y[3] = 5.00; pb_y[3] = 5.00;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_0[4] = 0.80; pb_0[4] = 0.80;

    pa_x[4] = 1.00; pb_x[4] = 1.00;

    pa_y[4] = 0.80; pb_y[4] = 0.80;

    pa_z[4] = 3.80; pb_z[4] = 3.80;

    pa_0 = buffa.data(52);

    pb_0 = buffb.data(52);

    pa_x = buffa.data(53);

    pb_x = buffb.data(53);

    pa_y = buffa.data(54);

    pb_y = buffb.data(54);

    pa_z = buffa.data(55);

    pb_z = buffb.data(55);

    pa_0[0] = 2.10; pb_0[0] = 2.10;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 2.20; pb_z[0] = 2.20;

    pa_0[1] = 1.40; pb_0[1] = 1.40;

    pa_x[1] = 2.20; pb_x[1] = 2.20;

    pa_y[1] = 3.80; pb_y[1] = 3.80;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_0[2] = 2.10; pb_0[2] = 2.10;

    pa_x[2] = 1.40; pb_x[2] = 1.40;

    pa_y[2] = 5.00; pb_y[2] = 5.00;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 1.70; pb_0[3] = 1.70;

    pa_x[3] = 2.10; pb_x[3] = 2.10;

    pa_y[3] = 1.40; pb_y[3] = 1.40;

    pa_z[3] = 4.50; pb_z[3] = 4.50;

    pa_0[4] = 1.70; pb_0[4] = 1.70;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 4.50; pb_y[4] = 4.50;

    pa_z[4] = 1.00; pb_z[4] = 1.00;

    pa_0 = buffa.data(56);

    pb_0 = buffb.data(56);

    pa_x = buffa.data(57);

    pb_x = buffb.data(57);

    pa_y = buffa.data(58);

    pb_y = buffb.data(58);

    pa_z = buffa.data(59);

    pb_z = buffb.data(59);

    pa_0[0] = 4.50; pb_0[0] = 4.50;

    pa_x[0] = 2.10; pb_x[0] = 2.10;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 2.10; pb_z[0] = 2.10;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 1.00; pb_x[1] = 1.00;

    pa_y[1] = 2.10; pb_y[1] = 2.10;

    pa_z[1] = 5.00; pb_z[1] = 5.00;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 1.70; pb_z[2] = 1.70;

    pa_0[3] = 1.00; pb_0[3] = 1.00;

    pa_x[3] = 1.40; pb_x[3] = 1.40;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_0[4] = 0.50; pb_0[4] = 0.50;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 0.80; pb_y[4] = 0.80;

    pa_z[4] = 2.10; pb_z[4] = 2.10;

    pa_0 = buffa.data(60);

    pb_0 = buffb.data(60);

    pa_x = buffa.data(61);

    pb_x = buffb.data(61);

    pa_y = buffa.data(62);

    pb_y = buffb.data(62);

    pa_z = buffa.data(63);

    pb_z = buffb.data(63);

    pa_0[0] = 0.50; pb_0[0] = 0.50;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 5.00; pb_y[0] = 5.00;

    pa_z[0] = 5.00; pb_z[0] = 5.00;

    pa_0[1] = 3.80; pb_0[1] = 3.80;

    pa_x[1] = 1.40; pb_x[1] = 1.40;

    pa_y[1] = 3.80; pb_y[1] = 3.80;

    pa_z[1] = 5.00; pb_z[1] = 5.00;

    pa_0[2] = 1.70; pb_0[2] = 1.70;

    pa_x[2] = 2.20; pb_x[2] = 2.20;

    pa_y[2] = 2.20; pb_y[2] = 2.20;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_0[3] = 2.10; pb_0[3] = 2.10;

    pa_x[3] = 0.50; pb_x[3] = 0.50;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_0[4] = 1.70; pb_0[4] = 1.70;

    pa_x[4] = 0.50; pb_x[4] = 0.50;

    pa_y[4] = 3.80; pb_y[4] = 3.80;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    pa_0 = buffa.data(64);

    pb_0 = buffb.data(64);

    pa_x = buffa.data(65);

    pb_x = buffb.data(65);

    pa_y = buffa.data(66);

    pb_y = buffb.data(66);

    pa_z = buffa.data(67);

    pb_z = buffb.data(67);

    pa_0[0] = 0.50; pb_0[0] = 0.50;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 4.50; pb_y[0] = 4.50;

    pa_z[0] = 0.50; pb_z[0] = 0.50;

    pa_0[1] = 4.50; pb_0[1] = 4.50;

    pa_x[1] = 2.10; pb_x[1] = 2.10;

    pa_y[1] = 2.10; pb_y[1] = 2.10;

    pa_z[1] = 0.50; pb_z[1] = 0.50;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 3.80; pb_y[2] = 3.80;

    pa_z[2] = 1.00; pb_z[2] = 1.00;

    pa_0[3] = 3.80; pb_0[3] = 3.80;

    pa_x[3] = 3.40; pb_x[3] = 3.40;

    pa_y[3] = 1.70; pb_y[3] = 1.70;

    pa_z[3] = 0.80; pb_z[3] = 0.80;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 1.00; pb_y[4] = 1.00;

    pa_z[4] = 3.80; pb_z[4] = 3.80;

    pa_0 = buffa.data(68);

    pb_0 = buffb.data(68);

    pa_x = buffa.data(69);

    pb_x = buffb.data(69);

    pa_y = buffa.data(70);

    pb_y = buffb.data(70);

    pa_z = buffa.data(71);

    pb_z = buffb.data(71);

    pa_0[0] = 1.40; pb_0[0] = 1.40;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 0.80; pb_y[0] = 0.80;

    pa_z[0] = 3.40; pb_z[0] = 3.40;

    pa_0[1] = 4.50; pb_0[1] = 4.50;

    pa_x[1] = 4.50; pb_x[1] = 4.50;

    pa_y[1] = 5.00; pb_y[1] = 5.00;

    pa_z[1] = 0.50; pb_z[1] = 0.50;

    pa_0[2] = 3.80; pb_0[2] = 3.80;

    pa_x[2] = 0.80; pb_x[2] = 0.80;

    pa_y[2] = 1.00; pb_y[2] = 1.00;

    pa_z[2] = 5.00; pb_z[2] = 5.00;

    pa_0[3] = 4.50; pb_0[3] = 4.50;

    pa_x[3] = 3.40; pb_x[3] = 3.40;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 1.70; pb_z[3] = 1.70;

    pa_0[4] = 0.80; pb_0[4] = 0.80;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 2.20; pb_y[4] = 2.20;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    pa_0 = buffa.data(72);

    pb_0 = buffb.data(72);

    pa_x = buffa.data(73);

    pb_x = buffb.data(73);

    pa_y = buffa.data(74);

    pb_y = buffb.data(74);

    pa_z = buffa.data(75);

    pb_z = buffb.data(75);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 3.40; pb_y[0] = 3.40;

    pa_z[0] = 3.40; pb_z[0] = 3.40;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 0.50; pb_x[1] = 0.50;

    pa_y[1] = 1.40; pb_y[1] = 1.40;

    pa_z[1] = 2.10; pb_z[1] = 2.10;

    pa_0[2] = 0.80; pb_0[2] = 0.80;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 0.50; pb_z[2] = 0.50;

    pa_0[3] = 0.50; pb_0[3] = 0.50;

    pa_x[3] = 1.70; pb_x[3] = 1.70;

    pa_y[3] = 2.10; pb_y[3] = 2.10;

    pa_z[3] = 1.40; pb_z[3] = 1.40;

    pa_0[4] = 4.50; pb_0[4] = 4.50;

    pa_x[4] = 5.00; pb_x[4] = 5.00;

    pa_y[4] = 3.80; pb_y[4] = 3.80;

    pa_z[4] = 5.00; pb_z[4] = 5.00;

    pa_0 = buffa.data(76);

    pb_0 = buffb.data(76);

    pa_x = buffa.data(77);

    pb_x = buffb.data(77);

    pa_y = buffa.data(78);

    pb_y = buffb.data(78);

    pa_z = buffa.data(79);

    pb_z = buffb.data(79);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 5.00; pb_x[0] = 5.00;

    pa_y[0] = 2.10; pb_y[0] = 2.10;

    pa_z[0] = 0.50; pb_z[0] = 0.50;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 1.40; pb_x[1] = 1.40;

    pa_y[1] = 4.50; pb_y[1] = 4.50;

    pa_z[1] = 2.10; pb_z[1] = 2.10;

    pa_0[2] = 2.10; pb_0[2] = 2.10;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 1.70; pb_z[2] = 1.70;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 1.00; pb_y[3] = 1.00;

    pa_z[3] = 4.50; pb_z[3] = 4.50;

    pa_0[4] = 0.50; pb_0[4] = 0.50;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    gtorec::compGtoTypeGForGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(80);

    pb_x = buffb.data(81);

    pb_y = buffb.data(82);

    pb_z = buffb.data(83);

    pb_0[0] = 1.00;

    pb_x[0] = 3.10;

    pb_y[0] = 0.50;

    pb_z[0] = 0.80;

    pb_0[1] = 2.52;

    pb_x[1] = 6.18;

    pb_y[1] = 2.64;

    pb_z[1] = 4.08;

    pb_0[2] = 6.60;

    pb_x[2] = 4.60;

    pb_y[2] = 5.10;

    pb_z[2] = 1.50;

    pb_0[3] = 4.00;

    pb_x[3] = 16.20;

    pb_y[3] = 2.00;

    pb_z[3] = 5.60;

    pb_0 = buffb.data(84);

    pb_x = buffb.data(85);

    pb_y = buffb.data(86);

    pb_z = buffb.data(87);

    pb_0[0] = 2.00;

    pb_x[0] = 4.20;

    pb_y[0] = 2.00;

    pb_z[0] = 1.60;

    pb_0[1] = 6.72;

    pb_x[1] = 10.88;

    pb_y[1] = 9.14;

    pb_z[1] = 10.88;

    pb_0[2] = 1.76;

    pb_x[2] = 0.64;

    pb_y[2] = 3.56;

    pb_z[2] = 0.40;

    pb_0[3] = 1.50;

    pb_x[3] = 5.70;

    pb_y[3] = 1.75;

    pb_z[3] = 2.10;

    pb_0 = buffb.data(88);

    pb_x = buffb.data(89);

    pb_y = buffb.data(90);

    pb_z = buffb.data(91);

    pb_0[0] = 7.20;

    pb_x[0] = 15.12;

    pb_y[0] = 3.60;

    pb_z[0] = 6.76;

    pb_0[1] = 16.38;

    pb_x[1] = 26.52;

    pb_y[1] = 17.16;

    pb_z[1] = 28.62;

    pb_0[2] = 1.98;

    pb_x[2] = 0.72;

    pb_y[2] = 1.53;

    pb_z[2] = 2.65;

    pb_0[3] = 1.90;

    pb_x[3] = 7.22;

    pb_y[3] = 0.95;

    pb_z[3] = 3.66;

    pb_0 = buffb.data(92);

    pb_x = buffb.data(93);

    pb_y = buffb.data(94);

    pb_z = buffb.data(95);

    pb_0[0] = 2.10;

    pb_x[0] = 4.30;

    pb_y[0] = 3.80;

    pb_z[0] = 2.20;

    pb_0[1] = 1.68;

    pb_x[1] = 4.04;

    pb_y[1] = 4.56;

    pb_z[1] = 0.96;

    pb_0[2] = 6.30;

    pb_x[2] = 6.30;

    pb_y[2] = 15.00;

    pb_z[2] = 6.60;

    pb_0[3] = 6.80;

    pb_x[3] = 10.10;

    pb_y[3] = 5.60;

    pb_z[3] = 18.00;

    pb_0 = buffb.data(96);

    pb_x = buffb.data(97);

    pb_y = buffb.data(98);

    pb_z = buffb.data(99);

    pb_0[0] = 4.50;

    pb_x[0] = 6.60;

    pb_y[0] = 3.80;

    pb_z[0] = 2.10;

    pb_0[1] = 4.08;

    pb_x[1] = 4.60;

    pb_y[1] = 2.52;

    pb_z[1] = 6.00;

    pb_0[2] = 11.40;

    pb_x[2] = 8.90;

    pb_y[2] = 10.20;

    pb_z[2] = 5.10;

    pb_0[3] = 4.00;

    pb_x[3] = 6.60;

    pb_y[3] = 2.00;

    pb_z[3] = 8.80;

    pb_0 = buffb.data(100);

    pb_x = buffb.data(101);

    pb_y = buffb.data(102);

    pb_z = buffb.data(103);

    pb_0[0] = 0.50;

    pb_x[0] = 5.00;

    pb_y[0] = 5.00;

    pb_z[0] = 5.00;

    pb_0[1] = 4.56;

    pb_x[1] = 5.48;

    pb_y[1] = 4.56;

    pb_z[1] = 6.00;

    pb_0[2] = 5.10;

    pb_x[2] = 8.30;

    pb_y[2] = 6.60;

    pb_z[2] = 6.60;

    pb_0[3] = 8.40;

    pb_x[3] = 4.10;

    pb_y[3] = 13.60;

    pb_z[3] = 8.80;

    pb_0 = buffb.data(104);

    pb_x = buffb.data(105);

    pb_y = buffb.data(106);

    pb_z = buffb.data(107);

    pb_0[0] = 0.50;

    pb_x[0] = 2.70;

    pb_y[0] = 4.50;

    pb_z[0] = 0.50;

    pb_0[1] = 5.40;

    pb_x[1] = 7.02;

    pb_y[1] = 2.52;

    pb_z[1] = 0.60;

    pb_0[2] = 11.40;

    pb_x[2] = 6.80;

    pb_y[2] = 11.40;

    pb_z[2] = 3.00;

    pb_0[3] = 15.20;

    pb_x[3] = 17.40;

    pb_y[3] = 6.80;

    pb_z[3] = 3.20;

    pb_0 = buffb.data(108);

    pb_x = buffb.data(109);

    pb_y = buffb.data(110);

    pb_z = buffb.data(111);

    pb_0[0] = 1.40;

    pb_x[0] = 3.60;

    pb_y[0] = 0.80;

    pb_z[0] = 3.40;

    pb_0[1] = 5.40;

    pb_x[1] = 9.90;

    pb_y[1] = 6.00;

    pb_z[1] = 0.60;

    pb_0[2] = 11.40;

    pb_x[2] = 6.20;

    pb_y[2] = 3.00;

    pb_z[2] = 15.00;

    pb_0[3] = 18.00;

    pb_x[3] = 18.10;

    pb_y[3] = 2.00;

    pb_z[3] = 6.80;

    pb_0 = buffb.data(112);

    pb_x = buffb.data(113);

    pb_y = buffb.data(114);

    pb_z = buffb.data(115);

    pb_0[0] = 3.80;

    pb_x[0] = 8.30;

    pb_y[0] = 3.40;

    pb_z[0] = 3.40;

    pb_0[1] = 2.52;

    pb_x[1] = 2.70;

    pb_y[1] = 1.68;

    pb_z[1] = 2.52;

    pb_0[2] = 2.40;

    pb_x[2] = 5.90;

    pb_y[2] = 10.20;

    pb_z[2] = 1.50;

    pb_0[3] = 2.00;

    pb_x[3] = 7.30;

    pb_y[3] = 8.40;

    pb_z[3] = 5.60;

    pb_0 = buffb.data(116);

    pb_x = buffb.data(117);

    pb_y = buffb.data(118);

    pb_z = buffb.data(119);

    pb_0[0] = 3.80;

    pb_x[0] = 8.80;

    pb_y[0] = 2.10;

    pb_z[0] = 0.50;

    pb_0[1] = 2.52;

    pb_x[1] = 3.78;

    pb_y[1] = 5.40;

    pb_z[1] = 2.52;

    pb_0[2] = 6.30;

    pb_x[2] = 7.20;

    pb_y[2] = 10.20;

    pb_z[2] = 5.10;

    pb_0[3] = 20.00;

    pb_x[3] = 23.00;

    pb_y[3] = 4.00;

    pb_z[3] = 18.00;

    pb_0 = buffb.data(120);

    pb_x = buffb.data(121);

    pb_y = buffb.data(122);

    pb_z = buffb.data(123);

    pb_0[0] = 1.00;

    pb_x[0] = 4.40;

    pb_y[0] = 9.50;

    pb_z[0] = 1.00;

    pb_0[1] = 14.40;

    pb_x[1] = 6.72;

    pb_y[1] = 11.22;

    pb_z[1] = 1.60;

    pb_0[2] = 3.04;

    pb_x[2] = 0.80;

    pb_y[2] = 6.84;

    pb_z[2] = 0.80;

    pb_0[3] = 5.70;

    pb_x[3] = 5.10;

    pb_y[3] = 6.35;

    pb_z[3] = 1.20;

    pb_0 = buffb.data(124);

    pb_x = buffb.data(125);

    pb_y = buffb.data(126);

    pb_z = buffb.data(127);

    pb_0[0] = 3.60;

    pb_x[0] = 15.84;

    pb_y[0] = 32.40;

    pb_z[0] = 4.10;

    pb_0[1] = 35.10;

    pb_x[1] = 16.38;

    pb_y[1] = 16.38;

    pb_z[1] = 8.40;

    pb_0[2] = 3.42;

    pb_x[2] = 0.90;

    pb_y[2] = 3.42;

    pb_z[2] = 4.70;

    pb_0[3] = 7.22;

    pb_x[3] = 6.46;

    pb_y[3] = 3.23;

    pb_z[3] = 5.32;

    pb_0 = buffb.data(128);

    pb_x = buffb.data(129);

    pb_y = buffb.data(130);

    pb_z = buffb.data(131);

    pb_0[0] = 10.08;

    pb_x[0] = 15.84;

    pb_y[0] = 5.76;

    pb_z[0] = 25.88;

    pb_0[1] = 35.10;

    pb_x[1] = 35.10;

    pb_y[1] = 39.00;

    pb_z[1] = 8.40;

    pb_0[2] = 3.42;

    pb_x[2] = 0.72;

    pb_y[2] = 0.90;

    pb_z[2] = 8.30;

    pb_0[3] = 8.55;

    pb_x[3] = 6.46;

    pb_y[3] = 0.95;

    pb_z[3] = 7.73;

    pb_0 = buffb.data(132);

    pb_x = buffb.data(133);

    pb_y = buffb.data(134);

    pb_z = buffb.data(135);

    pb_0[0] = 7.60;

    pb_x[0] = 10.00;

    pb_y[0] = 8.00;

    pb_z[0] = 1.00;

    pb_0[1] = 6.72;

    pb_x[1] = 4.48;

    pb_y[1] = 16.50;

    pb_z[1] = 6.72;

    pb_0[2] = 1.68;

    pb_x[2] = 1.36;

    pb_y[2] = 4.82;

    pb_z[2] = 1.36;

    pb_0[3] = 7.50;

    pb_x[3] = 6.75;

    pb_y[3] = 6.50;

    pb_z[3] = 6.75;

    pb_0 = buffb.data(136);

    pb_x = buffb.data(137);

    pb_y = buffb.data(138);

    pb_z = buffb.data(139);

    pb_0[0] = 27.36;

    pb_x[0] = 36.00;

    pb_y[0] = 15.12;

    pb_z[0] = 7.40;

    pb_0[1] = 16.38;

    pb_x[1] = 10.92;

    pb_y[1] = 35.10;

    pb_z[1] = 18.48;

    pb_0[2] = 1.89;

    pb_x[2] = 1.53;

    pb_y[2] = 3.06;

    pb_z[2] = 3.63;

    pb_0[3] = 9.50;

    pb_x[3] = 8.55;

    pb_y[3] = 1.90;

    pb_z[3] = 13.55;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypePForMGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 20);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 20);

    buffb.zero();

    auto pa_0 = buffa.data(0);

    auto pb_0 = buffb.data(0);

    auto pa_x = buffa.data(1);

    auto pb_x = buffb.data(1);

    auto pa_y = buffa.data(2);

    auto pb_y = buffb.data(2);

    auto pa_z = buffa.data(3);

    auto pb_z = buffb.data(3);

    auto pa_2 = buffa.data(4);

    auto pb_2 = buffb.data(4);

    pa_0[0] = 4.50; pb_0[0] = 4.50;

    pa_x[0] = 3.80; pb_x[0] = 3.80;

    pa_y[0] = 3.40; pb_y[0] = 3.40;

    pa_z[0] = 1.70; pb_z[0] = 1.70;

    pa_2[0] = 5.00; pb_2[0] = 5.00;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 2.10; pb_y[1] = 2.10;

    pa_z[1] = 1.40; pb_z[1] = 1.40;

    pa_2[1] = 1.00; pb_2[1] = 1.00;

    pa_0[2] = 2.10; pb_0[2] = 2.10;

    pa_x[2] = 3.80; pb_x[2] = 3.80;

    pa_y[2] = 2.20; pb_y[2] = 2.20;

    pa_z[2] = 5.00; pb_z[2] = 5.00;

    pa_2[2] = 1.40; pb_2[2] = 1.40;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 2.20; pb_x[3] = 2.20;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 3.40; pb_z[3] = 3.40;

    pa_2[3] = 2.20; pb_2[3] = 2.20;

    pa_0[4] = 0.80; pb_0[4] = 0.80;

    pa_x[4] = 1.70; pb_x[4] = 1.70;

    pa_y[4] = 4.50; pb_y[4] = 4.50;

    pa_z[4] = 2.20; pb_z[4] = 2.20;

    pa_2[4] = 1.00; pb_2[4] = 1.00;

    gtorec::compGtoTypePForMGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(5);

    pb_x = buffb.data(6);

    pb_y = buffb.data(7);

    pb_z = buffb.data(8);

    pb_2 = buffb.data(9);

    pb_0[0] = 4.50;

    pb_x[0] = 8.30;

    pb_y[0] = 3.40;

    pb_z[0] = 1.70;

    pb_2[0] = 12.60;

    pb_0[1] = 4.08;

    pb_x[1] = 9.40;

    pb_y[1] = 2.52;

    pb_z[1] = 1.68;

    pb_2[1] = 11.20;

    pb_0[2] = 6.30;

    pb_x[2] = 13.50;

    pb_y[2] = 6.60;

    pb_z[2] = 15.00;

    pb_2[2] = 11.80;

    pb_0[3] = 20.00;

    pb_x[3] = 13.80;

    pb_y[3] = 13.60;

    pb_z[3] = 13.60;

    pb_2[3] = 13.20;

    pb_0 = buffb.data(10);

    pb_x = buffb.data(11);

    pb_y = buffb.data(12);

    pb_z = buffb.data(13);

    pb_2 = buffb.data(14);

    pb_0[0] = 9.00;

    pb_x[0] = 7.60;

    pb_y[0] = 11.30;

    pb_z[0] = 3.40;

    pb_2[0] = 16.80;

    pb_0[1] = 10.88;

    pb_x[1] = 16.00;

    pb_y[1] = 10.12;

    pb_z[1] = 4.48;

    pb_2[1] = 7.40;

    pb_0[2] = 1.68;

    pb_x[2] = 3.04;

    pb_y[2] = 3.86;

    pb_z[2] = 4.00;

    pb_2[2] = 5.52;

    pb_0[3] = 7.50;

    pb_x[3] = 3.30;

    pb_y[3] = 10.10;

    pb_z[3] = 5.10;

    pb_2[3] = 10.10;

    pb_0 = buffb.data(15);

    pb_x = buffb.data(16);

    pb_y = buffb.data(17);

    pb_z = buffb.data(18);

    pb_2 = buffb.data(19);

    pb_0[0] = 32.40;

    pb_x[0] = 27.36;

    pb_y[0] = 24.48;

    pb_z[0] = 16.74;

    pb_2[0] = 39.40;

    pb_0[1] = 26.52;

    pb_x[1] = 39.00;

    pb_y[1] = 16.38;

    pb_z[1] = 14.32;

    pb_2[1] = 10.60;

    pb_0[2] = 1.89;

    pb_x[2] = 3.42;

    pb_y[2] = 1.98;

    pb_z[2] = 6.60;

    pb_2[2] = 11.26;

    pb_0[3] = 9.50;

    pb_x[3] = 4.18;

    pb_y[3] = 6.46;

    pb_z[3] = 11.46;

    pb_2[3] = 10.98;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeDForMGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 50);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 50);

    buffb.zero();

    auto pa_0 = buffa.data(5);

    auto pb_0 = buffb.data(5);

    auto pa_x = buffa.data(6);

    auto pb_x = buffb.data(6);

    auto pa_y = buffa.data(7);

    auto pb_y = buffb.data(7);

    auto pa_z = buffa.data(8);

    auto pb_z = buffb.data(8);

    auto pa_2 = buffa.data(9);

    auto pb_2 = buffb.data(9);

    pa_0[0] = 2.20; pb_0[0] = 2.20;

    pa_x[0] = 1.70; pb_x[0] = 1.70;

    pa_y[0] = 5.00; pb_y[0] = 5.00;

    pa_z[0] = 1.40; pb_z[0] = 1.40;

    pa_2[0] = 1.70; pb_2[0] = 1.70;

    pa_0[1] = 1.70; pb_0[1] = 1.70;

    pa_x[1] = 1.00; pb_x[1] = 1.00;

    pa_y[1] = 3.80; pb_y[1] = 3.80;

    pa_z[1] = 1.40; pb_z[1] = 1.40;

    pa_2[1] = 0.50; pb_2[1] = 0.50;

    pa_0[2] = 2.20; pb_0[2] = 2.20;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 1.40; pb_y[2] = 1.40;

    pa_z[2] = 5.00; pb_z[2] = 5.00;

    pa_2[2] = 0.80; pb_2[2] = 0.80;

    pa_0[3] = 0.80; pb_0[3] = 0.80;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 5.00; pb_z[3] = 5.00;

    pa_2[3] = 4.50; pb_2[3] = 4.50;

    pa_0[4] = 0.50; pb_0[4] = 0.50;

    pa_x[4] = 2.10; pb_x[4] = 2.10;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 0.50; pb_z[4] = 0.50;

    pa_2[4] = 2.20; pb_2[4] = 2.20;

    pa_0 = buffa.data(10);

    pb_0 = buffb.data(10);

    pa_x = buffa.data(11);

    pb_x = buffb.data(11);

    pa_y = buffa.data(12);

    pb_y = buffb.data(12);

    pa_z = buffa.data(13);

    pb_z = buffb.data(13);

    pa_2 = buffa.data(14);

    pb_2 = buffb.data(14);

    pa_0[0] = 2.10; pb_0[0] = 2.10;

    pa_x[0] = 3.80; pb_x[0] = 3.80;

    pa_y[0] = 0.80; pb_y[0] = 0.80;

    pa_z[0] = 2.20; pb_z[0] = 2.20;

    pa_2[0] = 2.20; pb_2[0] = 2.20;

    pa_0[1] = 1.00; pb_0[1] = 1.00;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 5.00; pb_y[1] = 5.00;

    pa_z[1] = 3.80; pb_z[1] = 3.80;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 1.40; pb_0[2] = 1.40;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 0.50; pb_y[2] = 0.50;

    pa_z[2] = 1.00; pb_z[2] = 1.00;

    pa_2[2] = 2.10; pb_2[2] = 2.10;

    pa_0[3] = 3.80; pb_0[3] = 3.80;

    pa_x[3] = 0.50; pb_x[3] = 0.50;

    pa_y[3] = 1.00; pb_y[3] = 1.00;

    pa_z[3] = 1.40; pb_z[3] = 1.40;

    pa_2[3] = 4.50; pb_2[3] = 4.50;

    pa_0[4] = 3.80; pb_0[4] = 3.80;

    pa_x[4] = 2.20; pb_x[4] = 2.20;

    pa_y[4] = 4.50; pb_y[4] = 4.50;

    pa_z[4] = 1.70; pb_z[4] = 1.70;

    pa_2[4] = 3.40; pb_2[4] = 3.40;

    pa_0 = buffa.data(15);

    pb_0 = buffb.data(15);

    pa_x = buffa.data(16);

    pb_x = buffb.data(16);

    pa_y = buffa.data(17);

    pb_y = buffb.data(17);

    pa_z = buffa.data(18);

    pb_z = buffb.data(18);

    pa_2 = buffa.data(19);

    pb_2 = buffb.data(19);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 1.40; pb_x[0] = 1.40;

    pa_y[0] = 2.10; pb_y[0] = 2.10;

    pa_z[0] = 2.20; pb_z[0] = 2.20;

    pa_2[0] = 2.10; pb_2[0] = 2.10;

    pa_0[1] = 0.50; pb_0[1] = 0.50;

    pa_x[1] = 2.10; pb_x[1] = 2.10;

    pa_y[1] = 0.80; pb_y[1] = 0.80;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_2[1] = 0.80; pb_2[1] = 0.80;

    pa_0[2] = 1.00; pb_0[2] = 1.00;

    pa_x[2] = 3.40; pb_x[2] = 3.40;

    pa_y[2] = 4.50; pb_y[2] = 4.50;

    pa_z[2] = 4.50; pb_z[2] = 4.50;

    pa_2[2] = 3.40; pb_2[2] = 3.40;

    pa_0[3] = 1.70; pb_0[3] = 1.70;

    pa_x[3] = 0.50; pb_x[3] = 0.50;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_2[3] = 2.10; pb_2[3] = 2.10;

    pa_0[4] = 1.70; pb_0[4] = 1.70;

    pa_x[4] = 1.70; pb_x[4] = 1.70;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 2.20; pb_z[4] = 2.20;

    pa_2[4] = 5.00; pb_2[4] = 5.00;

    gtorec::compGtoTypeDForMGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(20);

    pb_x = buffb.data(21);

    pb_y = buffb.data(22);

    pb_z = buffb.data(23);

    pb_2 = buffb.data(24);

    pb_0[0] = 2.20;

    pb_x[0] = 3.90;

    pb_y[0] = 5.00;

    pb_z[0] = 1.40;

    pb_2[0] = 5.10;

    pb_0[1] = 2.04;

    pb_x[1] = 2.90;

    pb_y[1] = 4.56;

    pb_z[1] = 1.68;

    pb_2[1] = 2.60;

    pb_0[2] = 6.60;

    pb_x[2] = 7.30;

    pb_y[2] = 4.20;

    pb_z[2] = 15.00;

    pb_2[2] = 5.80;

    pb_0[3] = 3.20;

    pb_x[3] = 4.80;

    pb_y[3] = 13.60;

    pb_z[3] = 20.00;

    pb_2[3] = 20.00;

    pb_0 = buffb.data(25);

    pb_x = buffb.data(26);

    pb_y = buffb.data(27);

    pb_z = buffb.data(28);

    pb_2 = buffb.data(29);

    pb_0[0] = 4.40;

    pb_x[0] = 3.40;

    pb_y[0] = 12.20;

    pb_z[0] = 2.80;

    pb_2[0] = 13.40;

    pb_0[1] = 5.44;

    pb_x[1] = 3.20;

    pb_y[1] = 13.86;

    pb_z[1] = 4.48;

    pb_2[1] = 9.20;

    pb_0[2] = 1.76;

    pb_x[2] = 1.36;

    pb_y[2] = 3.32;

    pb_z[2] = 4.00;

    pb_2[2] = 3.44;

    pb_0[3] = 1.20;

    pb_x[3] = 1.50;

    pb_y[3] = 5.90;

    pb_z[3] = 7.50;

    pb_2[3] = 13.55;

    pb_0 = buffb.data(30);

    pb_x = buffb.data(31);

    pb_y = buffb.data(32);

    pb_z = buffb.data(33);

    pb_2 = buffb.data(34);

    pb_0[0] = 15.84;

    pb_x[0] = 12.24;

    pb_y[0] = 36.00;

    pb_z[0] = 12.28;

    pb_2[0] = 15.04;

    pb_0[1] = 13.26;

    pb_x[1] = 7.80;

    pb_y[1] = 29.64;

    pb_z[1] = 12.62;

    pb_2[1] = 6.70;

    pb_0[2] = 1.98;

    pb_x[2] = 1.53;

    pb_y[2] = 1.26;

    pb_z[2] = 6.70;

    pb_2[2] = 10.72;

    pb_0[3] = 1.52;

    pb_x[3] = 1.90;

    pb_y[3] = 6.46;

    pb_z[3] = 10.30;

    pb_2[3] = 18.55;

    pb_0 = buffb.data(35);

    pb_x = buffb.data(36);

    pb_y = buffb.data(37);

    pb_z = buffb.data(38);

    pb_2 = buffb.data(39);

    pb_0[0] = 4.20;

    pb_x[0] = 7.60;

    pb_y[0] = 3.70;

    pb_z[0] = 4.40;

    pb_2[0] = 6.00;

    pb_0[1] = 3.20;

    pb_x[1] = 5.44;

    pb_y[1] = 17.00;

    pb_z[1] = 12.16;

    pb_2[1] = 26.00;

    pb_0[2] = 1.12;

    pb_x[2] = 0.80;

    pb_y[2] = 1.80;

    pb_z[2] = 0.80;

    pb_2[2] = 2.68;

    pb_0[3] = 5.70;

    pb_x[3] = 0.75;

    pb_y[3] = 5.30;

    pb_z[3] = 2.10;

    pb_2[3] = 8.75;

    pb_0 = buffb.data(40);

    pb_x = buffb.data(41);

    pb_y = buffb.data(42);

    pb_z = buffb.data(43);

    pb_2 = buffb.data(44);

    pb_0[0] = 15.12;

    pb_x[0] = 27.36;

    pb_y[0] = 5.76;

    pb_z[0] = 17.94;

    pb_2[0] = 20.24;

    pb_0[1] = 7.80;

    pb_x[1] = 13.26;

    pb_y[1] = 39.00;

    pb_z[1] = 30.64;

    pb_2[1] = 46.60;

    pb_0[2] = 1.26;

    pb_x[2] = 0.90;

    pb_y[2] = 0.45;

    pb_z[2] = 2.30;

    pb_2[2] = 3.89;

    pb_0[3] = 7.22;

    pb_x[3] = 0.95;

    pb_y[3] = 1.90;

    pb_z[3] = 6.46;

    pb_2[3] = 11.35;

    pb_0 = buffb.data(45);

    pb_x = buffb.data(46);

    pb_y = buffb.data(47);

    pb_z = buffb.data(48);

    pb_2 = buffb.data(49);

    pb_0[0] = 27.36;

    pb_x[0] = 10.08;

    pb_y[0] = 15.12;

    pb_z[0] = 19.64;

    pb_2[0] = 19.52;

    pb_0[1] = 3.90;

    pb_x[1] = 16.38;

    pb_y[1] = 6.24;

    pb_z[1] = 6.74;

    pb_2[1] = 7.84;

    pb_0[2] = 0.90;

    pb_x[2] = 3.06;

    pb_y[2] = 4.05;

    pb_z[2] = 5.05;

    pb_2[2] = 12.06;

    pb_0[3] = 3.23;

    pb_x[3] = 0.95;

    pb_y[3] = 6.46;

    pb_z[3] = 2.65;

    pb_2[3] = 4.99;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeFForMGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 100);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 100);

    buffb.zero();

    auto pa_0 = buffa.data(20);

    auto pb_0 = buffb.data(20);

    auto pa_x = buffa.data(21);

    auto pb_x = buffb.data(21);

    auto pa_y = buffa.data(22);

    auto pb_y = buffb.data(22);

    auto pa_z = buffa.data(23);

    auto pb_z = buffb.data(23);

    auto pa_2 = buffa.data(24);

    auto pb_2 = buffb.data(24);

    pa_0[0] = 4.50; pb_0[0] = 4.50;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 3.40; pb_y[0] = 3.40;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_2[0] = 3.40; pb_2[0] = 3.40;

    pa_0[1] = 1.70; pb_0[1] = 1.70;

    pa_x[1] = 4.50; pb_x[1] = 4.50;

    pa_y[1] = 0.50; pb_y[1] = 0.50;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_2[1] = 0.50; pb_2[1] = 0.50;

    pa_0[2] = 1.70; pb_0[2] = 1.70;

    pa_x[2] = 1.00; pb_x[2] = 1.00;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 3.80; pb_z[2] = 3.80;

    pa_2[2] = 1.40; pb_2[2] = 1.40;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 1.70; pb_x[3] = 1.70;

    pa_y[3] = 3.40; pb_y[3] = 3.40;

    pa_z[3] = 2.10; pb_z[3] = 2.10;

    pa_2[3] = 0.80; pb_2[3] = 0.80;

    pa_0[4] = 2.20; pb_0[4] = 2.20;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 3.80; pb_y[4] = 3.80;

    pa_z[4] = 1.00; pb_z[4] = 1.00;

    pa_2[4] = 5.00; pb_2[4] = 5.00;

    pa_0 = buffa.data(25);

    pb_0 = buffb.data(25);

    pa_x = buffa.data(26);

    pb_x = buffb.data(26);

    pa_y = buffa.data(27);

    pb_y = buffb.data(27);

    pa_z = buffa.data(28);

    pb_z = buffb.data(28);

    pa_2 = buffa.data(29);

    pb_2 = buffb.data(29);

    pa_0[0] = 1.00; pb_0[0] = 1.00;

    pa_x[0] = 1.00; pb_x[0] = 1.00;

    pa_y[0] = 1.70; pb_y[0] = 1.70;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_2[0] = 1.70; pb_2[0] = 1.70;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 0.80; pb_x[1] = 0.80;

    pa_y[1] = 1.40; pb_y[1] = 1.40;

    pa_z[1] = 1.00; pb_z[1] = 1.00;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 0.50; pb_0[2] = 0.50;

    pa_x[2] = 0.80; pb_x[2] = 0.80;

    pa_y[2] = 0.50; pb_y[2] = 0.50;

    pa_z[2] = 0.50; pb_z[2] = 0.50;

    pa_2[2] = 4.50; pb_2[2] = 4.50;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_2[3] = 3.80; pb_2[3] = 3.80;

    pa_0[4] = 2.20; pb_0[4] = 2.20;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    pa_2[4] = 0.80; pb_2[4] = 0.80;

    pa_0 = buffa.data(30);

    pb_0 = buffb.data(30);

    pa_x = buffa.data(31);

    pb_x = buffb.data(31);

    pa_y = buffa.data(32);

    pb_y = buffb.data(32);

    pa_z = buffa.data(33);

    pb_z = buffb.data(33);

    pa_2 = buffa.data(34);

    pb_2 = buffb.data(34);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 1.40; pb_y[0] = 1.40;

    pa_z[0] = 1.70; pb_z[0] = 1.70;

    pa_2[0] = 4.50; pb_2[0] = 4.50;

    pa_0[1] = 1.00; pb_0[1] = 1.00;

    pa_x[1] = 1.40; pb_x[1] = 1.40;

    pa_y[1] = 1.70; pb_y[1] = 1.70;

    pa_z[1] = 1.00; pb_z[1] = 1.00;

    pa_2[1] = 2.20; pb_2[1] = 2.20;

    pa_0[2] = 1.70; pb_0[2] = 1.70;

    pa_x[2] = 0.50; pb_x[2] = 0.50;

    pa_y[2] = 2.20; pb_y[2] = 2.20;

    pa_z[2] = 5.00; pb_z[2] = 5.00;

    pa_2[2] = 2.10; pb_2[2] = 2.10;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 3.80; pb_y[3] = 3.80;

    pa_z[3] = 4.50; pb_z[3] = 4.50;

    pa_2[3] = 2.20; pb_2[3] = 2.20;

    pa_0[4] = 4.50; pb_0[4] = 4.50;

    pa_x[4] = 5.00; pb_x[4] = 5.00;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 2.20; pb_z[4] = 2.20;

    pa_2[4] = 1.40; pb_2[4] = 1.40;

    pa_0 = buffa.data(35);

    pb_0 = buffb.data(35);

    pa_x = buffa.data(36);

    pb_x = buffb.data(36);

    pa_y = buffa.data(37);

    pb_y = buffb.data(37);

    pa_z = buffa.data(38);

    pb_z = buffb.data(38);

    pa_2 = buffa.data(39);

    pb_2 = buffb.data(39);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 1.00; pb_y[0] = 1.00;

    pa_z[0] = 4.50; pb_z[0] = 4.50;

    pa_2[0] = 1.00; pb_2[0] = 1.00;

    pa_0[1] = 2.10; pb_0[1] = 2.10;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 1.00; pb_y[1] = 1.00;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_2[1] = 2.20; pb_2[1] = 2.20;

    pa_0[2] = 2.20; pb_0[2] = 2.20;

    pa_x[2] = 2.20; pb_x[2] = 2.20;

    pa_y[2] = 1.70; pb_y[2] = 1.70;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_2[2] = 0.50; pb_2[2] = 0.50;

    pa_0[3] = 1.70; pb_0[3] = 1.70;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 0.80; pb_z[3] = 0.80;

    pa_2[3] = 3.80; pb_2[3] = 3.80;

    pa_0[4] = 2.20; pb_0[4] = 2.20;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 0.50; pb_y[4] = 0.50;

    pa_z[4] = 1.70; pb_z[4] = 1.70;

    pa_2[4] = 2.20; pb_2[4] = 2.20;

    pa_0 = buffa.data(40);

    pb_0 = buffb.data(40);

    pa_x = buffa.data(41);

    pb_x = buffb.data(41);

    pa_y = buffa.data(42);

    pb_y = buffb.data(42);

    pa_z = buffa.data(43);

    pb_z = buffb.data(43);

    pa_2 = buffa.data(44);

    pb_2 = buffb.data(44);

    pa_0[0] = 1.00; pb_0[0] = 1.00;

    pa_x[0] = 5.00; pb_x[0] = 5.00;

    pa_y[0] = 2.10; pb_y[0] = 2.10;

    pa_z[0] = 3.40; pb_z[0] = 3.40;

    pa_2[0] = 0.50; pb_2[0] = 0.50;

    pa_0[1] = 2.20; pb_0[1] = 2.20;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 3.80; pb_z[1] = 3.80;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 5.00; pb_0[2] = 5.00;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 3.80; pb_y[2] = 3.80;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_2[2] = 4.50; pb_2[2] = 4.50;

    pa_0[3] = 3.80; pb_0[3] = 3.80;

    pa_x[3] = 0.50; pb_x[3] = 0.50;

    pa_y[3] = 1.70; pb_y[3] = 1.70;

    pa_z[3] = 3.80; pb_z[3] = 3.80;

    pa_2[3] = 4.50; pb_2[3] = 4.50;

    pa_0[4] = 5.00; pb_0[4] = 5.00;

    pa_x[4] = 3.80; pb_x[4] = 3.80;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 4.50; pb_z[4] = 4.50;

    pa_2[4] = 1.00; pb_2[4] = 1.00;

    pa_0 = buffa.data(45);

    pb_0 = buffb.data(45);

    pa_x = buffa.data(46);

    pb_x = buffb.data(46);

    pa_y = buffa.data(47);

    pb_y = buffb.data(47);

    pa_z = buffa.data(48);

    pb_z = buffb.data(48);

    pa_2 = buffa.data(49);

    pb_2 = buffb.data(49);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 5.00; pb_y[0] = 5.00;

    pa_z[0] = 1.40; pb_z[0] = 1.40;

    pa_2[0] = 5.00; pb_2[0] = 5.00;

    pa_0[1] = 5.00; pb_0[1] = 5.00;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 4.50; pb_y[1] = 4.50;

    pa_z[1] = 4.50; pb_z[1] = 4.50;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 4.50; pb_0[2] = 4.50;

    pa_x[2] = 2.10; pb_x[2] = 2.10;

    pa_y[2] = 1.00; pb_y[2] = 1.00;

    pa_z[2] = 3.80; pb_z[2] = 3.80;

    pa_2[2] = 1.40; pb_2[2] = 1.40;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 0.80; pb_x[3] = 0.80;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 1.70; pb_z[3] = 1.70;

    pa_2[3] = 1.40; pb_2[3] = 1.40;

    pa_0[4] = 0.50; pb_0[4] = 0.50;

    pa_x[4] = 2.10; pb_x[4] = 2.10;

    pa_y[4] = 0.50; pb_y[4] = 0.50;

    pa_z[4] = 0.50; pb_z[4] = 0.50;

    pa_2[4] = 1.40; pb_2[4] = 1.40;

    gtorec::compGtoTypeFForMGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(50);

    pb_x = buffb.data(51);

    pb_y = buffb.data(52);

    pb_z = buffb.data(53);

    pb_2 = buffb.data(54);

    pb_0[0] = 4.50;

    pb_x[0] = 9.00;

    pb_y[0] = 3.40;

    pb_z[0] = 0.80;

    pb_2[0] = 12.40;

    pb_0[1] = 2.04;

    pb_x[1] = 7.10;

    pb_y[1] = 0.60;

    pb_z[1] = 2.64;

    pb_2[1] = 9.60;

    pb_0[2] = 5.10;

    pb_x[2] = 4.70;

    pb_y[2] = 10.20;

    pb_z[2] = 11.40;

    pb_2[2] = 6.20;

    pb_0[3] = 8.80;

    pb_x[3] = 9.00;

    pb_y[3] = 13.60;

    pb_z[3] = 8.40;

    pb_2[3] = 6.60;

    pb_0 = buffb.data(55);

    pb_x = buffb.data(56);

    pb_y = buffb.data(57);

    pb_z = buffb.data(58);

    pb_2 = buffb.data(59);

    pb_0[0] = 9.00;

    pb_x[0] = 9.00;

    pb_y[0] = 11.30;

    pb_z[0] = 1.60;

    pb_2[0] = 13.60;

    pb_0[1] = 5.44;

    pb_x[1] = 14.40;

    pb_y[1] = 3.30;

    pb_z[1] = 7.04;

    pb_2[1] = 2.60;

    pb_0[2] = 1.36;

    pb_x[2] = 0.80;

    pb_y[2] = 4.42;

    pb_z[2] = 3.04;

    pb_2[2] = 7.92;

    pb_0[3] = 3.30;

    pb_x[3] = 2.55;

    pb_y[3] = 7.30;

    pb_z[3] = 3.15;

    pb_2[3] = 8.00;

    pb_0 = buffb.data(60);

    pb_x = buffb.data(61);

    pb_y = buffb.data(62);

    pb_z = buffb.data(63);

    pb_2 = buffb.data(64);

    pb_0[0] = 32.40;

    pb_x[0] = 32.40;

    pb_y[0] = 24.48;

    pb_z[0] = 10.26;

    pb_2[0] = 26.08;

    pb_0[1] = 13.26;

    pb_x[1] = 35.10;

    pb_y[1] = 3.90;

    pb_z[1] = 18.86;

    pb_2[1] = 8.30;

    pb_0[2] = 1.53;

    pb_x[2] = 0.90;

    pb_y[2] = 3.06;

    pb_z[2] = 5.12;

    pb_2[2] = 8.86;

    pb_0[3] = 4.18;

    pb_x[3] = 3.23;

    pb_y[3] = 6.46;

    pb_z[3] = 6.19;

    pb_2[3] = 5.72;

    pb_0 = buffb.data(65);

    pb_x = buffb.data(66);

    pb_y = buffb.data(67);

    pb_z = buffb.data(68);

    pb_2 = buffb.data(69);

    pb_0[0] = 3.40;

    pb_x[0] = 7.90;

    pb_y[0] = 1.00;

    pb_z[0] = 4.50;

    pb_2[0] = 10.00;

    pb_0[1] = 2.52;

    pb_x[1] = 4.14;

    pb_y[1] = 1.20;

    pb_z[1] = 0.96;

    pb_2[1] = 6.04;

    pb_0[2] = 6.60;

    pb_x[2] = 8.80;

    pb_y[2] = 5.10;

    pb_z[2] = 6.60;

    pb_2[2] = 5.90;

    pb_0[3] = 6.80;

    pb_x[3] = 5.70;

    pb_y[3] = 2.00;

    pb_z[3] = 3.20;

    pb_2[3] = 17.20;

    pb_0 = buffb.data(70);

    pb_x = buffb.data(71);

    pb_y = buffb.data(72);

    pb_z = buffb.data(73);

    pb_2 = buffb.data(74);

    pb_0[0] = 1.00;

    pb_x[0] = 6.00;

    pb_y[0] = 2.10;

    pb_z[0] = 3.40;

    pb_2[0] = 10.50;

    pb_0[1] = 2.64;

    pb_x[1] = 8.20;

    pb_y[1] = 2.64;

    pb_z[1] = 4.56;

    pb_2[1] = 16.00;

    pb_0[2] = 15.00;

    pb_x[2] = 18.50;

    pb_y[2] = 11.40;

    pb_z[2] = 6.60;

    pb_2[2] = 22.50;

    pb_0[3] = 15.20;

    pb_x[3] = 5.80;

    pb_y[3] = 6.80;

    pb_z[3] = 15.20;

    pb_2[3] = 19.00;

    pb_0 = buffb.data(75);

    pb_x = buffb.data(76);

    pb_y = buffb.data(77);

    pb_z = buffb.data(78);

    pb_2 = buffb.data(79);

    pb_0[0] = 3.80;

    pb_x[0] = 6.00;

    pb_y[0] = 5.00;

    pb_z[0] = 1.40;

    pb_2[0] = 9.40;

    pb_0[1] = 6.00;

    pb_x[1] = 7.04;

    pb_y[1] = 5.40;

    pb_z[1] = 5.40;

    pb_2[1] = 9.40;

    pb_0[2] = 13.50;

    pb_x[2] = 10.80;

    pb_y[2] = 3.00;

    pb_z[2] = 11.40;

    pb_2[2] = 8.40;

    pb_0[3] = 8.80;

    pb_x[3] = 5.40;

    pb_y[3] = 2.00;

    pb_z[3] = 6.80;

    pb_2[3] = 7.20;

    pb_0 = buffb.data(80);

    pb_x = buffb.data(81);

    pb_y = buffb.data(82);

    pb_z = buffb.data(83);

    pb_2 = buffb.data(84);

    pb_0[0] = 6.80;

    pb_x[0] = 9.00;

    pb_y[0] = 5.40;

    pb_z[0] = 9.00;

    pb_2[0] = 4.00;

    pb_0[1] = 6.72;

    pb_x[1] = 5.44;

    pb_y[1] = 5.30;

    pb_z[1] = 2.56;

    pb_2[1] = 9.04;

    pb_0[2] = 1.76;

    pb_x[2] = 1.76;

    pb_y[2] = 3.56;

    pb_z[2] = 1.76;

    pb_2[2] = 3.80;

    pb_0[3] = 2.55;

    pb_x[3] = 1.50;

    pb_y[3] = 2.45;

    pb_z[3] = 1.20;

    pb_2[3] = 6.70;

    pb_0 = buffb.data(85);

    pb_x = buffb.data(86);

    pb_y = buffb.data(87);

    pb_z = buffb.data(88);

    pb_2 = buffb.data(89);

    pb_0[0] = 24.48;

    pb_x[0] = 32.40;

    pb_y[0] = 7.20;

    pb_z[0] = 35.80;

    pb_2[0] = 16.20;

    pb_0[1] = 16.38;

    pb_x[1] = 13.26;

    pb_y[1] = 7.80;

    pb_z[1] = 8.34;

    pb_2[1] = 18.76;

    pb_0[2] = 1.98;

    pb_x[2] = 1.98;

    pb_y[2] = 1.53;

    pb_z[2] = 4.18;

    pb_2[2] = 4.85;

    pb_0[3] = 3.23;

    pb_x[3] = 1.90;

    pb_y[3] = 0.95;

    pb_z[3] = 3.22;

    pb_2[3] = 8.82;

    pb_0 = buffb.data(90);

    pb_x = buffb.data(91);

    pb_y = buffb.data(92);

    pb_z = buffb.data(93);

    pb_2 = buffb.data(94);

    pb_0[0] = 7.60;

    pb_x[0] = 4.40;

    pb_y[0] = 13.80;

    pb_z[0] = 2.80;

    pb_2[0] = 20.00;

    pb_0[1] = 16.00;

    pb_x[1] = 5.44;

    pb_y[1] = 19.40;

    pb_z[1] = 14.40;

    pb_2[1] = 25.00;

    pb_0[2] = 3.60;

    pb_x[2] = 1.68;

    pb_y[2] = 5.30;

    pb_z[2] = 3.04;

    pb_2[2] = 3.12;

    pb_0[3] = 3.30;

    pb_x[3] = 1.20;

    pb_y[3] = 2.95;

    pb_z[3] = 2.55;

    pb_2[3] = 3.10;

    pb_0 = buffb.data(95);

    pb_x = buffb.data(96);

    pb_y = buffb.data(97);

    pb_z = buffb.data(98);

    pb_2 = buffb.data(99);

    pb_0[0] = 27.36;

    pb_x[0] = 15.84;

    pb_y[0] = 36.00;

    pb_z[0] = 13.88;

    pb_2[0] = 38.80;

    pb_0[1] = 39.00;

    pb_x[1] = 13.26;

    pb_y[1] = 35.10;

    pb_z[1] = 40.10;

    pb_2[1] = 48.00;

    pb_0[2] = 4.05;

    pb_x[2] = 1.89;

    pb_y[2] = 0.90;

    pb_z[2] = 7.92;

    pb_2[2] = 8.86;

    pb_0[3] = 4.18;

    pb_x[3] = 1.52;

    pb_y[3] = 0.95;

    pb_z[3] = 5.43;

    pb_2[3] = 6.06;

    ASSERT_EQ(buffa, buffb);
}

TEST_F(CGtoRecFuncTest, CompGtoTypeGForMGGA)
{
    CMemBlock2D<double> dist({1.00, 1.20, 3.00, 4.00, 5.60,
                              2.00, 3.20, 0.80, 1.50, 8.00,
                              7.20, 7.80, 0.90, 1.90, 2.10},
                              5, 3);

    CMemBlock2D<int32_t> ridx({2, 4, 1, 3}, 2, 2);

    CMemBlock2D<double> buffa(5, 175);

    buffa.zero();

    CMemBlock2D<double> buffb(5, 175);

    buffb.zero();

    auto pa_0 = buffa.data(50);

    auto pb_0 = buffb.data(50);

    auto pa_x = buffa.data(51);

    auto pb_x = buffb.data(51);

    auto pa_y = buffa.data(52);

    auto pb_y = buffb.data(52);

    auto pa_z = buffa.data(53);

    auto pb_z = buffb.data(53);

    auto pa_2 = buffa.data(54);

    auto pb_2 = buffb.data(54);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 2.10; pb_x[0] = 2.10;

    pa_y[0] = 2.20; pb_y[0] = 2.20;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_2[0] = 1.40; pb_2[0] = 1.40;

    pa_0[1] = 1.00; pb_0[1] = 1.00;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 2.10; pb_y[1] = 2.10;

    pa_z[1] = 2.10; pb_z[1] = 2.10;

    pa_2[1] = 3.80; pb_2[1] = 3.80;

    pa_0[2] = 1.70; pb_0[2] = 1.70;

    pa_x[2] = 2.20; pb_x[2] = 2.20;

    pa_y[2] = 1.70; pb_y[2] = 1.70;

    pa_z[2] = 5.00; pb_z[2] = 5.00;

    pa_2[2] = 3.80; pb_2[2] = 3.80;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 2.20; pb_y[3] = 2.20;

    pa_z[3] = 5.00; pb_z[3] = 5.00;

    pa_2[3] = 1.40; pb_2[3] = 1.40;

    pa_0[4] = 5.00; pb_0[4] = 5.00;

    pa_x[4] = 0.80; pb_x[4] = 0.80;

    pa_y[4] = 2.10; pb_y[4] = 2.10;

    pa_z[4] = 2.10; pb_z[4] = 2.10;

    pa_2[4] = 1.00; pb_2[4] = 1.00;

    pa_0 = buffa.data(55);

    pb_0 = buffb.data(55);

    pa_x = buffa.data(56);

    pb_x = buffb.data(56);

    pa_y = buffa.data(57);

    pb_y = buffb.data(57);

    pa_z = buffa.data(58);

    pb_z = buffb.data(58);

    pa_2 = buffa.data(59);

    pb_2 = buffb.data(59);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 3.40; pb_x[0] = 3.40;

    pa_y[0] = 1.40; pb_y[0] = 1.40;

    pa_z[0] = 2.10; pb_z[0] = 2.10;

    pa_2[0] = 4.50; pb_2[0] = 4.50;

    pa_0[1] = 2.20; pb_0[1] = 2.20;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 0.80; pb_y[1] = 0.80;

    pa_z[1] = 4.50; pb_z[1] = 4.50;

    pa_2[1] = 0.50; pb_2[1] = 0.50;

    pa_0[2] = 1.00; pb_0[2] = 1.00;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 4.50; pb_y[2] = 4.50;

    pa_z[2] = 4.50; pb_z[2] = 4.50;

    pa_2[2] = 1.70; pb_2[2] = 1.70;

    pa_0[3] = 1.00; pb_0[3] = 1.00;

    pa_x[3] = 0.80; pb_x[3] = 0.80;

    pa_y[3] = 5.00; pb_y[3] = 5.00;

    pa_z[3] = 4.50; pb_z[3] = 4.50;

    pa_2[3] = 0.80; pb_2[3] = 0.80;

    pa_0[4] = 4.50; pb_0[4] = 4.50;

    pa_x[4] = 0.80; pb_x[4] = 0.80;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 1.40; pb_z[4] = 1.40;

    pa_2[4] = 0.50; pb_2[4] = 0.50;

    pa_0 = buffa.data(60);

    pb_0 = buffb.data(60);

    pa_x = buffa.data(61);

    pb_x = buffb.data(61);

    pa_y = buffa.data(62);

    pb_y = buffb.data(62);

    pa_z = buffa.data(63);

    pb_z = buffb.data(63);

    pa_2 = buffa.data(64);

    pb_2 = buffb.data(64);

    pa_0[0] = 1.40; pb_0[0] = 1.40;

    pa_x[0] = 4.50; pb_x[0] = 4.50;

    pa_y[0] = 3.80; pb_y[0] = 3.80;

    pa_z[0] = 2.10; pb_z[0] = 2.10;

    pa_2[0] = 1.40; pb_2[0] = 1.40;

    pa_0[1] = 3.40; pb_0[1] = 3.40;

    pa_x[1] = 5.00; pb_x[1] = 5.00;

    pa_y[1] = 1.40; pb_y[1] = 1.40;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_2[1] = 2.10; pb_2[1] = 2.10;

    pa_0[2] = 3.40; pb_0[2] = 3.40;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 0.50; pb_y[2] = 0.50;

    pa_z[2] = 3.40; pb_z[2] = 3.40;

    pa_2[2] = 5.00; pb_2[2] = 5.00;

    pa_0[3] = 4.50; pb_0[3] = 4.50;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 3.80; pb_y[3] = 3.80;

    pa_z[3] = 2.10; pb_z[3] = 2.10;

    pa_2[3] = 3.80; pb_2[3] = 3.80;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 3.80; pb_x[4] = 3.80;

    pa_y[4] = 5.00; pb_y[4] = 5.00;

    pa_z[4] = 5.00; pb_z[4] = 5.00;

    pa_2[4] = 0.80; pb_2[4] = 0.80;

    pa_0 = buffa.data(65);

    pb_0 = buffb.data(65);

    pa_x = buffa.data(66);

    pb_x = buffb.data(66);

    pa_y = buffa.data(67);

    pb_y = buffb.data(67);

    pa_z = buffa.data(68);

    pb_z = buffb.data(68);

    pa_2 = buffa.data(69);

    pb_2 = buffb.data(69);

    pa_0[0] = 2.20; pb_0[0] = 2.20;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 2.10; pb_y[0] = 2.10;

    pa_z[0] = 0.50; pb_z[0] = 0.50;

    pa_2[0] = 1.40; pb_2[0] = 1.40;

    pa_0[1] = 1.40; pb_0[1] = 1.40;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 5.00; pb_y[1] = 5.00;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_2[1] = 0.50; pb_2[1] = 0.50;

    pa_0[2] = 2.10; pb_0[2] = 2.10;

    pa_x[2] = 0.50; pb_x[2] = 0.50;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 1.70; pb_z[2] = 1.70;

    pa_2[2] = 0.50; pb_2[2] = 0.50;

    pa_0[3] = 3.40; pb_0[3] = 3.40;

    pa_x[3] = 1.00; pb_x[3] = 1.00;

    pa_y[3] = 0.80; pb_y[3] = 0.80;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_2[3] = 2.10; pb_2[3] = 2.10;

    pa_0[4] = 5.00; pb_0[4] = 5.00;

    pa_x[4] = 2.20; pb_x[4] = 2.20;

    pa_y[4] = 3.80; pb_y[4] = 3.80;

    pa_z[4] = 4.50; pb_z[4] = 4.50;

    pa_2[4] = 0.50; pb_2[4] = 0.50;

    pa_0 = buffa.data(70);

    pb_0 = buffb.data(70);

    pa_x = buffa.data(71);

    pb_x = buffb.data(71);

    pa_y = buffa.data(72);

    pb_y = buffb.data(72);

    pa_z = buffa.data(73);

    pb_z = buffb.data(73);

    pa_2 = buffa.data(74);

    pb_2 = buffb.data(74);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 3.40; pb_x[0] = 3.40;

    pa_y[0] = 0.80; pb_y[0] = 0.80;

    pa_z[0] = 1.40; pb_z[0] = 1.40;

    pa_2[0] = 2.20; pb_2[0] = 2.20;

    pa_0[1] = 0.50; pb_0[1] = 0.50;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 2.20; pb_y[1] = 2.20;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_2[1] = 1.70; pb_2[1] = 1.70;

    pa_0[2] = 0.80; pb_0[2] = 0.80;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 1.70; pb_y[2] = 1.70;

    pa_z[2] = 1.40; pb_z[2] = 1.40;

    pa_2[2] = 0.80; pb_2[2] = 0.80;

    pa_0[3] = 1.00; pb_0[3] = 1.00;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 1.00; pb_y[3] = 1.00;

    pa_z[3] = 1.00; pb_z[3] = 1.00;

    pa_2[3] = 1.00; pb_2[3] = 1.00;

    pa_0[4] = 3.80; pb_0[4] = 3.80;

    pa_x[4] = 1.40; pb_x[4] = 1.40;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 0.80; pb_z[4] = 0.80;

    pa_2[4] = 3.40; pb_2[4] = 3.40;

    pa_0 = buffa.data(75);

    pb_0 = buffb.data(75);

    pa_x = buffa.data(76);

    pb_x = buffb.data(76);

    pa_y = buffa.data(77);

    pb_y = buffb.data(77);

    pa_z = buffa.data(78);

    pb_z = buffb.data(78);

    pa_2 = buffa.data(79);

    pb_2 = buffb.data(79);

    pa_0[0] = 1.70; pb_0[0] = 1.70;

    pa_x[0] = 0.80; pb_x[0] = 0.80;

    pa_y[0] = 2.20; pb_y[0] = 2.20;

    pa_z[0] = 1.00; pb_z[0] = 1.00;

    pa_2[0] = 1.40; pb_2[0] = 1.40;

    pa_0[1] = 1.70; pb_0[1] = 1.70;

    pa_x[1] = 2.10; pb_x[1] = 2.10;

    pa_y[1] = 1.70; pb_y[1] = 1.70;

    pa_z[1] = 3.80; pb_z[1] = 3.80;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 1.70; pb_0[2] = 1.70;

    pa_x[2] = 1.40; pb_x[2] = 1.40;

    pa_y[2] = 3.80; pb_y[2] = 3.80;

    pa_z[2] = 2.10; pb_z[2] = 2.10;

    pa_2[2] = 3.40; pb_2[2] = 3.40;

    pa_0[3] = 5.00; pb_0[3] = 5.00;

    pa_x[3] = 0.80; pb_x[3] = 0.80;

    pa_y[3] = 1.70; pb_y[3] = 1.70;

    pa_z[3] = 1.70; pb_z[3] = 1.70;

    pa_2[3] = 3.80; pb_2[3] = 3.80;

    pa_0[4] = 0.50; pb_0[4] = 0.50;

    pa_x[4] = 2.10; pb_x[4] = 2.10;

    pa_y[4] = 4.50; pb_y[4] = 4.50;

    pa_z[4] = 2.20; pb_z[4] = 2.20;

    pa_2[4] = 0.50; pb_2[4] = 0.50;

    pa_0 = buffa.data(80);

    pb_0 = buffb.data(80);

    pa_x = buffa.data(81);

    pb_x = buffb.data(81);

    pa_y = buffa.data(82);

    pb_y = buffb.data(82);

    pa_z = buffa.data(83);

    pb_z = buffb.data(83);

    pa_2 = buffa.data(84);

    pb_2 = buffb.data(84);

    pa_0[0] = 3.80; pb_0[0] = 3.80;

    pa_x[0] = 2.20; pb_x[0] = 2.20;

    pa_y[0] = 1.40; pb_y[0] = 1.40;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_2[0] = 2.10; pb_2[0] = 2.10;

    pa_0[1] = 1.00; pb_0[1] = 1.00;

    pa_x[1] = 1.40; pb_x[1] = 1.40;

    pa_y[1] = 5.00; pb_y[1] = 5.00;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_2[1] = 0.50; pb_2[1] = 0.50;

    pa_0[2] = 2.20; pb_0[2] = 2.20;

    pa_x[2] = 1.70; pb_x[2] = 1.70;

    pa_y[2] = 2.20; pb_y[2] = 2.20;

    pa_z[2] = 0.50; pb_z[2] = 0.50;

    pa_2[2] = 3.40; pb_2[2] = 3.40;

    pa_0[3] = 3.40; pb_0[3] = 3.40;

    pa_x[3] = 1.40; pb_x[3] = 1.40;

    pa_y[3] = 0.50; pb_y[3] = 0.50;

    pa_z[3] = 3.80; pb_z[3] = 3.80;

    pa_2[3] = 0.80; pb_2[3] = 0.80;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 1.70; pb_x[4] = 1.70;

    pa_y[4] = 5.00; pb_y[4] = 5.00;

    pa_z[4] = 4.50; pb_z[4] = 4.50;

    pa_2[4] = 5.00; pb_2[4] = 5.00;

    pa_0 = buffa.data(85);

    pb_0 = buffb.data(85);

    pa_x = buffa.data(86);

    pb_x = buffb.data(86);

    pa_y = buffa.data(87);

    pb_y = buffb.data(87);

    pa_z = buffa.data(88);

    pb_z = buffb.data(88);

    pa_2 = buffa.data(89);

    pb_2 = buffb.data(89);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 1.70; pb_x[0] = 1.70;

    pa_y[0] = 4.50; pb_y[0] = 4.50;

    pa_z[0] = 0.80; pb_z[0] = 0.80;

    pa_2[0] = 1.00; pb_2[0] = 1.00;

    pa_0[1] = 0.50; pb_0[1] = 0.50;

    pa_x[1] = 1.70; pb_x[1] = 1.70;

    pa_y[1] = 1.00; pb_y[1] = 1.00;

    pa_z[1] = 2.20; pb_z[1] = 2.20;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 2.10; pb_0[2] = 2.10;

    pa_x[2] = 4.50; pb_x[2] = 4.50;

    pa_y[2] = 1.70; pb_y[2] = 1.70;

    pa_z[2] = 1.40; pb_z[2] = 1.40;

    pa_2[2] = 1.40; pb_2[2] = 1.40;

    pa_0[3] = 2.20; pb_0[3] = 2.20;

    pa_x[3] = 0.50; pb_x[3] = 0.50;

    pa_y[3] = 5.00; pb_y[3] = 5.00;

    pa_z[3] = 1.00; pb_z[3] = 1.00;

    pa_2[3] = 0.80; pb_2[3] = 0.80;

    pa_0[4] = 1.00; pb_0[4] = 1.00;

    pa_x[4] = 3.80; pb_x[4] = 3.80;

    pa_y[4] = 1.70; pb_y[4] = 1.70;

    pa_z[4] = 1.40; pb_z[4] = 1.40;

    pa_2[4] = 3.80; pb_2[4] = 3.80;

    pa_0 = buffa.data(90);

    pb_0 = buffb.data(90);

    pa_x = buffa.data(91);

    pb_x = buffb.data(91);

    pa_y = buffa.data(92);

    pb_y = buffb.data(92);

    pa_z = buffa.data(93);

    pb_z = buffb.data(93);

    pa_2 = buffa.data(94);

    pb_2 = buffb.data(94);

    pa_0[0] = 5.00; pb_0[0] = 5.00;

    pa_x[0] = 3.80; pb_x[0] = 3.80;

    pa_y[0] = 1.00; pb_y[0] = 1.00;

    pa_z[0] = 1.40; pb_z[0] = 1.40;

    pa_2[0] = 0.80; pb_2[0] = 0.80;

    pa_0[1] = 4.50; pb_0[1] = 4.50;

    pa_x[1] = 3.40; pb_x[1] = 3.40;

    pa_y[1] = 1.00; pb_y[1] = 1.00;

    pa_z[1] = 4.50; pb_z[1] = 4.50;

    pa_2[1] = 5.00; pb_2[1] = 5.00;

    pa_0[2] = 0.50; pb_0[2] = 0.50;

    pa_x[2] = 1.40; pb_x[2] = 1.40;

    pa_y[2] = 0.80; pb_y[2] = 0.80;

    pa_z[2] = 0.80; pb_z[2] = 0.80;

    pa_2[2] = 0.50; pb_2[2] = 0.50;

    pa_0[3] = 4.50; pb_0[3] = 4.50;

    pa_x[3] = 1.40; pb_x[3] = 1.40;

    pa_y[3] = 4.50; pb_y[3] = 4.50;

    pa_z[3] = 0.50; pb_z[3] = 0.50;

    pa_2[3] = 1.40; pb_2[3] = 1.40;

    pa_0[4] = 3.40; pb_0[4] = 3.40;

    pa_x[4] = 4.50; pb_x[4] = 4.50;

    pa_y[4] = 1.00; pb_y[4] = 1.00;

    pa_z[4] = 2.10; pb_z[4] = 2.10;

    pa_2[4] = 5.00; pb_2[4] = 5.00;

    pa_0 = buffa.data(95);

    pb_0 = buffb.data(95);

    pa_x = buffa.data(96);

    pb_x = buffb.data(96);

    pa_y = buffa.data(97);

    pb_y = buffb.data(97);

    pa_z = buffa.data(98);

    pb_z = buffb.data(98);

    pa_2 = buffa.data(99);

    pb_2 = buffb.data(99);

    pa_0[0] = 3.40; pb_0[0] = 3.40;

    pa_x[0] = 3.40; pb_x[0] = 3.40;

    pa_y[0] = 0.50; pb_y[0] = 0.50;

    pa_z[0] = 5.00; pb_z[0] = 5.00;

    pa_2[0] = 0.50; pb_2[0] = 0.50;

    pa_0[1] = 0.50; pb_0[1] = 0.50;

    pa_x[1] = 3.80; pb_x[1] = 3.80;

    pa_y[1] = 1.40; pb_y[1] = 1.40;

    pa_z[1] = 0.80; pb_z[1] = 0.80;

    pa_2[1] = 3.80; pb_2[1] = 3.80;

    pa_0[2] = 1.40; pb_0[2] = 1.40;

    pa_x[2] = 3.80; pb_x[2] = 3.80;

    pa_y[2] = 3.40; pb_y[2] = 3.40;

    pa_z[2] = 2.20; pb_z[2] = 2.20;

    pa_2[2] = 5.00; pb_2[2] = 5.00;

    pa_0[3] = 3.80; pb_0[3] = 3.80;

    pa_x[3] = 4.50; pb_x[3] = 4.50;

    pa_y[3] = 4.50; pb_y[3] = 4.50;

    pa_z[3] = 2.20; pb_z[3] = 2.20;

    pa_2[3] = 1.70; pb_2[3] = 1.70;

    pa_0[4] = 3.40; pb_0[4] = 3.40;

    pa_x[4] = 3.40; pb_x[4] = 3.40;

    pa_y[4] = 1.40; pb_y[4] = 1.40;

    pa_z[4] = 3.80; pb_z[4] = 3.80;

    pa_2[4] = 2.10; pb_2[4] = 2.10;

    gtorec::compGtoTypeGForMGGA(buffa, dist, ridx, 1);

    pb_0 = buffb.data(100);

    pb_x = buffb.data(101);

    pb_y = buffb.data(102);

    pb_z = buffb.data(103);

    pb_2 = buffb.data(104);

    pb_0[0] = 3.40;

    pb_x[0] = 5.50;

    pb_y[0] = 2.20;

    pb_z[0] = 0.80;

    pb_2[0] = 5.60;

    pb_0[1] = 1.20;

    pb_x[1] = 7.00;

    pb_y[1] = 2.52;

    pb_z[1] = 2.52;

    pb_2[1] = 14.56;

    pb_0[2] = 5.10;

    pb_x[2] = 8.30;

    pb_y[2] = 5.10;

    pb_z[2] = 15.00;

    pb_2[2] = 15.80;

    pb_0[3] = 20.00;

    pb_x[3] = 23.00;

    pb_y[3] = 8.80;

    pb_z[3] = 20.00;

    pb_2[3] = 14.60;

    pb_0 = buffb.data(105);

    pb_x = buffb.data(106);

    pb_y = buffb.data(107);

    pb_z = buffb.data(108);

    pb_2 = buffb.data(109);

    pb_0[0] = 6.80;

    pb_x[0] = 4.20;

    pb_y[0] = 7.80;

    pb_z[0] = 1.60;

    pb_2[0] = 7.20;

    pb_0[1] = 3.20;

    pb_x[1] = 16.00;

    pb_y[1] = 7.72;

    pb_z[1] = 6.72;

    pb_2[1] = 16.36;

    pb_0[2] = 1.36;

    pb_x[2] = 1.76;

    pb_y[2] = 3.06;

    pb_z[2] = 4.00;

    pb_2[2] = 6.44;

    pb_0[3] = 7.50;

    pb_x[3] = 6.75;

    pb_y[3] = 8.30;

    pb_z[3] = 7.50;

    pb_2[3] = 6.50;

    pb_0 = buffb.data(110);

    pb_x = buffb.data(111);

    pb_y = buffb.data(112);

    pb_z = buffb.data(113);

    pb_2 = buffb.data(114);

    pb_0[0] = 24.48;

    pb_x[0] = 15.12;

    pb_y[0] = 15.84;

    pb_z[0] = 9.16;

    pb_2[0] = 11.68;

    pb_0[1] = 7.80;

    pb_x[1] = 39.00;

    pb_y[1] = 16.38;

    pb_z[1] = 17.38;

    pb_2[1] = 33.84;

    pb_0[2] = 1.53;

    pb_x[2] = 1.98;

    pb_y[2] = 1.53;

    pb_z[2] = 6.20;

    pb_2[2] = 13.42;

    pb_0[3] = 9.50;

    pb_x[3] = 8.55;

    pb_y[3] = 4.18;

    pb_z[3] = 14.50;

    pb_2[3] = 12.66;

    pb_0 = buffb.data(115);

    pb_x = buffb.data(116);

    pb_y = buffb.data(117);

    pb_z = buffb.data(118);

    pb_2 = buffb.data(119);

    pb_0[0] = 2.20;

    pb_x[0] = 4.40;

    pb_y[0] = 2.10;

    pb_z[0] = 0.50;

    pb_2[0] = 5.80;

    pb_0[1] = 1.68;

    pb_x[1] = 3.44;

    pb_y[1] = 6.00;

    pb_z[1] = 2.64;

    pb_2[1] = 4.00;

    pb_0[2] = 6.30;

    pb_x[2] = 3.60;

    pb_y[2] = 10.20;

    pb_z[2] = 5.10;

    pb_2[2] = 2.50;

    pb_0[3] = 13.60;

    pb_x[3] = 7.40;

    pb_y[3] = 3.20;

    pb_z[3] = 8.80;

    pb_2[3] = 10.40;

    pb_0 = buffb.data(120);

    pb_x = buffb.data(121);

    pb_y = buffb.data(122);

    pb_z = buffb.data(123);

    pb_2 = buffb.data(124);

    pb_0[0] = 3.40;

    pb_x[0] = 6.80;

    pb_y[0] = 0.80;

    pb_z[0] = 1.40;

    pb_2[0] = 9.00;

    pb_0[1] = 0.60;

    pb_x[1] = 2.54;

    pb_y[1] = 2.64;

    pb_z[1] = 0.96;

    pb_2[1] = 5.44;

    pb_0[2] = 2.40;

    pb_x[2] = 14.30;

    pb_y[2] = 5.10;

    pb_z[2] = 4.20;

    pb_2[2] = 11.40;

    pb_0[3] = 4.00;

    pb_x[3] = 19.00;

    pb_y[3] = 4.00;

    pb_z[3] = 4.00;

    pb_2[3] = 13.00;

    pb_0 = buffb.data(125);

    pb_x = buffb.data(126);

    pb_y = buffb.data(127);

    pb_z = buffb.data(128);

    pb_2 = buffb.data(129);

    pb_0[0] = 1.70;

    pb_x[0] = 2.50;

    pb_y[0] = 2.20;

    pb_z[0] = 1.00;

    pb_2[0] = 3.00;

    pb_0[1] = 2.04;

    pb_x[1] = 4.22;

    pb_y[1] = 2.04;

    pb_z[1] = 4.56;

    pb_2[1] = 10.20;

    pb_0[2] = 5.10;

    pb_x[2] = 5.90;

    pb_y[2] = 11.40;

    pb_z[2] = 6.30;

    pb_2[2] = 13.00;

    pb_0[3] = 20.00;

    pb_x[3] = 8.20;

    pb_y[3] = 6.80;

    pb_z[3] = 6.80;

    pb_2[3] = 16.80;

    pb_0 = buffb.data(130);

    pb_x = buffb.data(131);

    pb_y = buffb.data(132);

    pb_z = buffb.data(133);

    pb_2 = buffb.data(134);

    pb_0[0] = 3.80;

    pb_x[0] = 6.00;

    pb_y[0] = 1.40;

    pb_z[0] = 0.80;

    pb_2[0] = 6.50;

    pb_0[1] = 1.20;

    pb_x[1] = 2.68;

    pb_y[1] = 6.00;

    pb_z[1] = 2.64;

    pb_2[1] = 3.40;

    pb_0[2] = 6.60;

    pb_x[2] = 7.30;

    pb_y[2] = 6.60;

    pb_z[2] = 1.50;

    pb_2[2] = 13.60;

    pb_0[3] = 13.60;

    pb_x[3] = 9.00;

    pb_y[3] = 2.00;

    pb_z[3] = 15.20;

    pb_2[3] = 6.00;

    pb_0 = buffb.data(135);

    pb_x = buffb.data(136);

    pb_y = buffb.data(137);

    pb_z = buffb.data(138);

    pb_2 = buffb.data(139);

    pb_0[0] = 3.40;

    pb_x[0] = 5.10;

    pb_y[0] = 4.50;

    pb_z[0] = 0.80;

    pb_2[0] = 4.40;

    pb_0[1] = 0.60;

    pb_x[1] = 2.54;

    pb_y[1] = 1.20;

    pb_z[1] = 2.64;

    pb_2[1] = 9.40;

    pb_0[2] = 6.30;

    pb_x[2] = 15.60;

    pb_y[2] = 5.10;

    pb_z[2] = 4.20;

    pb_2[2] = 13.20;

    pb_0[3] = 8.80;

    pb_x[3] = 4.20;

    pb_y[3] = 20.00;

    pb_z[3] = 4.00;

    pb_2[3] = 4.20;

    pb_0 = buffb.data(140);

    pb_x = buffb.data(141);

    pb_y = buffb.data(142);

    pb_z = buffb.data(143);

    pb_2 = buffb.data(144);

    pb_0[0] = 5.00;

    pb_x[0] = 8.80;

    pb_y[0] = 1.00;

    pb_z[0] = 1.40;

    pb_2[0] = 8.40;

    pb_0[1] = 5.40;

    pb_x[1] = 8.58;

    pb_y[1] = 1.20;

    pb_z[1] = 5.40;

    pb_2[1] = 12.80;

    pb_0[2] = 1.50;

    pb_x[2] = 4.70;

    pb_y[2] = 2.40;

    pb_z[2] = 2.40;

    pb_2[2] = 4.30;

    pb_0[3] = 18.00;

    pb_x[3] = 10.10;

    pb_y[3] = 18.00;

    pb_z[3] = 2.00;

    pb_2[3] = 8.40;

    pb_0 = buffb.data(145);

    pb_x = buffb.data(146);

    pb_y = buffb.data(147);

    pb_z = buffb.data(148);

    pb_2 = buffb.data(149);

    pb_0[0] = 3.40;

    pb_x[0] = 6.80;

    pb_y[0] = 0.50;

    pb_z[0] = 5.00;

    pb_2[0] = 7.30;

    pb_0[1] = 0.60;

    pb_x[1] = 5.06;

    pb_y[1] = 1.68;

    pb_z[1] = 0.96;

    pb_2[1] = 12.16;

    pb_0[2] = 4.20;

    pb_x[2] = 12.80;

    pb_y[2] = 10.20;

    pb_z[2] = 6.60;

    pb_2[2] = 22.60;

    pb_0[3] = 15.20;

    pb_x[3] = 21.80;

    pb_y[3] = 18.00;

    pb_z[3] = 8.80;

    pb_2[3] = 15.80;

    pb_0 = buffb.data(150);

    pb_x = buffb.data(151);

    pb_y = buffb.data(152);

    pb_z = buffb.data(153);

    pb_2 = buffb.data(154);

    pb_0[0] = 7.60;

    pb_x[0] = 4.40;

    pb_y[0] = 6.60;

    pb_z[0] = 1.60;

    pb_2[0] = 7.00;

    pb_0[1] = 3.20;

    pb_x[1] = 4.48;

    pb_y[1] = 17.00;

    pb_z[1] = 7.04;

    pb_2[1] = 11.60;

    pb_0[2] = 1.76;

    pb_x[2] = 1.36;

    pb_y[2] = 3.96;

    pb_z[2] = 0.40;

    pb_2[2] = 7.12;

    pb_0[3] = 5.10;

    pb_x[3] = 2.10;

    pb_y[3] = 4.15;

    pb_z[3] = 5.70;

    pb_2[3] = 2.20;

    pb_0 = buffb.data(155);

    pb_x = buffb.data(156);

    pb_y = buffb.data(157);

    pb_z = buffb.data(158);

    pb_2 = buffb.data(159);

    pb_0[0] = 27.36;

    pb_x[0] = 15.84;

    pb_y[0] = 10.08;

    pb_z[0] = 9.56;

    pb_2[0] = 16.72;

    pb_0[1] = 7.80;

    pb_x[1] = 10.92;

    pb_y[1] = 39.00;

    pb_z[1] = 18.16;

    pb_2[1] = 8.30;

    pb_0[2] = 1.98;

    pb_x[2] = 1.53;

    pb_y[2] = 1.98;

    pb_z[2] = 2.65;

    pb_2[2] = 4.06;

    pb_0[3] = 6.46;

    pb_x[3] = 2.66;

    pb_y[3] = 0.95;

    pb_z[3] = 10.62;

    pb_2[3] = 9.12;

    pb_0 = buffb.data(160);

    pb_x = buffb.data(161);

    pb_y = buffb.data(162);

    pb_z = buffb.data(163);

    pb_2 = buffb.data(164);

    pb_0[0] = 24.48;

    pb_x[0] = 12.24;

    pb_y[0] = 32.40;

    pb_z[0] = 9.16;

    pb_2[0] = 8.80;

    pb_0[1] = 3.90;

    pb_x[1] = 13.26;

    pb_y[1] = 7.80;

    pb_z[1] = 17.66;

    pb_2[1] = 43.40;

    pb_0[2] = 1.89;

    pb_x[2] = 4.05;

    pb_y[2] = 1.53;

    pb_z[2] = 3.36;

    pb_2[2] = 4.06;

    pb_0[3] = 4.18;

    pb_x[3] = 0.95;

    pb_y[3] = 9.50;

    pb_z[3] = 4.10;

    pb_2[3] = 3.52;

    pb_0 = buffb.data(165);

    pb_x = buffb.data(166);

    pb_y = buffb.data(167);

    pb_z = buffb.data(168);

    pb_2 = buffb.data(169);

    pb_0[0] = 6.80;

    pb_x[0] = 6.80;

    pb_y[0] = 4.40;

    pb_z[0] = 10.00;

    pb_2[0] = 2.00;

    pb_0[1] = 1.60;

    pb_x[1] = 12.16;

    pb_y[1] = 4.98;

    pb_z[1] = 2.56;

    pb_2[1] = 14.96;

    pb_0[2] = 1.12;

    pb_x[2] = 3.04;

    pb_y[2] = 4.12;

    pb_z[2] = 1.76;

    pb_2[2] = 10.80;

    pb_0[3] = 5.70;

    pb_x[3] = 6.75;

    pb_y[3] = 10.55;

    pb_z[3] = 3.30;

    pb_2[3] = 11.55;

    pb_0 = buffb.data(170);

    pb_x = buffb.data(171);

    pb_y = buffb.data(172);

    pb_z = buffb.data(173);

    pb_2 = buffb.data(174);

    pb_0[0] = 24.48;

    pb_x[0] = 24.48;

    pb_y[0] = 3.60;

    pb_z[0] = 39.40;

    pb_2[0] = 13.60;

    pb_0[1] = 3.90;

    pb_x[1] = 29.64;

    pb_y[1] = 10.92;

    pb_z[1] = 6.74;

    pb_2[1] = 31.24;

    pb_0[2] = 1.26;

    pb_x[2] = 3.42;

    pb_y[2] = 3.06;

    pb_z[2] = 3.38;

    pb_2[2] = 8.90;

    pb_0[3] = 7.22;

    pb_x[3] = 8.55;

    pb_y[3] = 8.55;

    pb_z[3] = 7.98;

    pb_2[3] = 7.63;

    ASSERT_EQ(buffa, buffb);
}

