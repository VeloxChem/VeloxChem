#include "NuclearPotentialGeom010PrimRecFD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fd,
                                        const size_t              idx_npot_geom_010_0_pd,
                                        const size_t              idx_npot_geom_010_1_pd,
                                        const size_t              idx_npot_geom_010_0_dp,
                                        const size_t              idx_npot_geom_010_1_dp,
                                        const size_t              idx_npot_1_dd,
                                        const size_t              idx_npot_geom_010_0_dd,
                                        const size_t              idx_npot_geom_010_1_dd,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : PD

    auto ta1_x_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd);

    auto ta1_x_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 1);

    auto ta1_x_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 2);

    auto ta1_x_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 3);

    auto ta1_x_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 4);

    auto ta1_x_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 5);

    auto ta1_x_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 6);

    auto ta1_x_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 7);

    auto ta1_x_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 8);

    auto ta1_x_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 9);

    auto ta1_x_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 10);

    auto ta1_x_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 11);

    auto ta1_x_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 12);

    auto ta1_x_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 13);

    auto ta1_x_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 14);

    auto ta1_x_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 15);

    auto ta1_x_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 16);

    auto ta1_x_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 17);

    auto ta1_y_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 18);

    auto ta1_y_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 19);

    auto ta1_y_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 20);

    auto ta1_y_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 21);

    auto ta1_y_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 22);

    auto ta1_y_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 23);

    auto ta1_y_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 24);

    auto ta1_y_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 25);

    auto ta1_y_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 26);

    auto ta1_y_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 27);

    auto ta1_y_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 28);

    auto ta1_y_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 29);

    auto ta1_y_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 30);

    auto ta1_y_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 31);

    auto ta1_y_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 32);

    auto ta1_y_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 33);

    auto ta1_y_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 34);

    auto ta1_y_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 35);

    auto ta1_z_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 36);

    auto ta1_z_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 37);

    auto ta1_z_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 38);

    auto ta1_z_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 39);

    auto ta1_z_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 40);

    auto ta1_z_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 41);

    auto ta1_z_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 42);

    auto ta1_z_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 43);

    auto ta1_z_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 44);

    auto ta1_z_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 45);

    auto ta1_z_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 46);

    auto ta1_z_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 47);

    auto ta1_z_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 48);

    auto ta1_z_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 49);

    auto ta1_z_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 50);

    auto ta1_z_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 51);

    auto ta1_z_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 52);

    auto ta1_z_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 53);

    // Set up components of auxiliary buffer : PD

    auto ta1_x_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd);

    auto ta1_x_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 1);

    auto ta1_x_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 2);

    auto ta1_x_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 3);

    auto ta1_x_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 4);

    auto ta1_x_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 5);

    auto ta1_x_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 6);

    auto ta1_x_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 7);

    auto ta1_x_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 8);

    auto ta1_x_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 9);

    auto ta1_x_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 10);

    auto ta1_x_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 11);

    auto ta1_x_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 12);

    auto ta1_x_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 13);

    auto ta1_x_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 14);

    auto ta1_x_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 15);

    auto ta1_x_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 16);

    auto ta1_x_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 17);

    auto ta1_y_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 18);

    auto ta1_y_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 19);

    auto ta1_y_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 20);

    auto ta1_y_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 21);

    auto ta1_y_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 22);

    auto ta1_y_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 23);

    auto ta1_y_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 24);

    auto ta1_y_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 25);

    auto ta1_y_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 26);

    auto ta1_y_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 27);

    auto ta1_y_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 28);

    auto ta1_y_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 29);

    auto ta1_y_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 30);

    auto ta1_y_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 31);

    auto ta1_y_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 32);

    auto ta1_y_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 33);

    auto ta1_y_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 34);

    auto ta1_y_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 35);

    auto ta1_z_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 36);

    auto ta1_z_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 37);

    auto ta1_z_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 38);

    auto ta1_z_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 39);

    auto ta1_z_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 40);

    auto ta1_z_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 41);

    auto ta1_z_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 42);

    auto ta1_z_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 43);

    auto ta1_z_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 44);

    auto ta1_z_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 45);

    auto ta1_z_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 46);

    auto ta1_z_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 47);

    auto ta1_z_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 48);

    auto ta1_z_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 49);

    auto ta1_z_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 50);

    auto ta1_z_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 51);

    auto ta1_z_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 52);

    auto ta1_z_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 53);

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp);

    auto ta1_x_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 1);

    auto ta1_x_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 2);

    auto ta1_x_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 9);

    auto ta1_x_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 10);

    auto ta1_x_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 11);

    auto ta1_x_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 15);

    auto ta1_x_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 16);

    auto ta1_x_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 17);

    auto ta1_y_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 18);

    auto ta1_y_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 19);

    auto ta1_y_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 20);

    auto ta1_y_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 27);

    auto ta1_y_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 28);

    auto ta1_y_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 29);

    auto ta1_y_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 33);

    auto ta1_y_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 34);

    auto ta1_y_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 35);

    auto ta1_z_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 36);

    auto ta1_z_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 37);

    auto ta1_z_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 38);

    auto ta1_z_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 45);

    auto ta1_z_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 46);

    auto ta1_z_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 47);

    auto ta1_z_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 51);

    auto ta1_z_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 52);

    auto ta1_z_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 53);

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp);

    auto ta1_x_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 1);

    auto ta1_x_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 2);

    auto ta1_x_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 9);

    auto ta1_x_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 10);

    auto ta1_x_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 11);

    auto ta1_x_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 15);

    auto ta1_x_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 16);

    auto ta1_x_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 17);

    auto ta1_y_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 18);

    auto ta1_y_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 19);

    auto ta1_y_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 20);

    auto ta1_y_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 27);

    auto ta1_y_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 28);

    auto ta1_y_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 29);

    auto ta1_y_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 33);

    auto ta1_y_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 34);

    auto ta1_y_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 35);

    auto ta1_z_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 36);

    auto ta1_z_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 37);

    auto ta1_z_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 38);

    auto ta1_z_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 45);

    auto ta1_z_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 46);

    auto ta1_z_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 47);

    auto ta1_z_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 51);

    auto ta1_z_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 52);

    auto ta1_z_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 53);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_1 = pbuffer.data(idx_npot_1_dd);

    auto ta_xx_xy_1 = pbuffer.data(idx_npot_1_dd + 1);

    auto ta_xx_xz_1 = pbuffer.data(idx_npot_1_dd + 2);

    auto ta_xx_yy_1 = pbuffer.data(idx_npot_1_dd + 3);

    auto ta_xx_yz_1 = pbuffer.data(idx_npot_1_dd + 4);

    auto ta_xx_zz_1 = pbuffer.data(idx_npot_1_dd + 5);

    auto ta_xy_xy_1 = pbuffer.data(idx_npot_1_dd + 7);

    auto ta_xz_xz_1 = pbuffer.data(idx_npot_1_dd + 14);

    auto ta_yy_xx_1 = pbuffer.data(idx_npot_1_dd + 18);

    auto ta_yy_xy_1 = pbuffer.data(idx_npot_1_dd + 19);

    auto ta_yy_xz_1 = pbuffer.data(idx_npot_1_dd + 20);

    auto ta_yy_yy_1 = pbuffer.data(idx_npot_1_dd + 21);

    auto ta_yy_yz_1 = pbuffer.data(idx_npot_1_dd + 22);

    auto ta_yy_zz_1 = pbuffer.data(idx_npot_1_dd + 23);

    auto ta_yz_yz_1 = pbuffer.data(idx_npot_1_dd + 28);

    auto ta_zz_xx_1 = pbuffer.data(idx_npot_1_dd + 30);

    auto ta_zz_xy_1 = pbuffer.data(idx_npot_1_dd + 31);

    auto ta_zz_xz_1 = pbuffer.data(idx_npot_1_dd + 32);

    auto ta_zz_yy_1 = pbuffer.data(idx_npot_1_dd + 33);

    auto ta_zz_yz_1 = pbuffer.data(idx_npot_1_dd + 34);

    auto ta_zz_zz_1 = pbuffer.data(idx_npot_1_dd + 35);

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd);

    auto ta1_x_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 1);

    auto ta1_x_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 2);

    auto ta1_x_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 3);

    auto ta1_x_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 4);

    auto ta1_x_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 5);

    auto ta1_x_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 6);

    auto ta1_x_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 7);

    auto ta1_x_xy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 8);

    auto ta1_x_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 9);

    auto ta1_x_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 12);

    auto ta1_x_xz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 13);

    auto ta1_x_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 14);

    auto ta1_x_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 17);

    auto ta1_x_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 18);

    auto ta1_x_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 19);

    auto ta1_x_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 20);

    auto ta1_x_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 21);

    auto ta1_x_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 22);

    auto ta1_x_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 23);

    auto ta1_x_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 26);

    auto ta1_x_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 28);

    auto ta1_x_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 29);

    auto ta1_x_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 30);

    auto ta1_x_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 31);

    auto ta1_x_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 32);

    auto ta1_x_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 33);

    auto ta1_x_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 34);

    auto ta1_x_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 35);

    auto ta1_y_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 36);

    auto ta1_y_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 37);

    auto ta1_y_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 38);

    auto ta1_y_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 39);

    auto ta1_y_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 40);

    auto ta1_y_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 41);

    auto ta1_y_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 42);

    auto ta1_y_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 43);

    auto ta1_y_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 45);

    auto ta1_y_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 46);

    auto ta1_y_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 50);

    auto ta1_y_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 52);

    auto ta1_y_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 53);

    auto ta1_y_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 54);

    auto ta1_y_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 55);

    auto ta1_y_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 56);

    auto ta1_y_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 57);

    auto ta1_y_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 58);

    auto ta1_y_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 59);

    auto ta1_y_yz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 61);

    auto ta1_y_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 63);

    auto ta1_y_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 64);

    auto ta1_y_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 65);

    auto ta1_y_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 66);

    auto ta1_y_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 67);

    auto ta1_y_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 68);

    auto ta1_y_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 69);

    auto ta1_y_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 70);

    auto ta1_y_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 71);

    auto ta1_z_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 72);

    auto ta1_z_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 73);

    auto ta1_z_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 74);

    auto ta1_z_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 75);

    auto ta1_z_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 76);

    auto ta1_z_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 77);

    auto ta1_z_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 79);

    auto ta1_z_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 81);

    auto ta1_z_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 82);

    auto ta1_z_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 84);

    auto ta1_z_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 86);

    auto ta1_z_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 88);

    auto ta1_z_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 89);

    auto ta1_z_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 90);

    auto ta1_z_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 91);

    auto ta1_z_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 92);

    auto ta1_z_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 93);

    auto ta1_z_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 94);

    auto ta1_z_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 95);

    auto ta1_z_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 98);

    auto ta1_z_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 99);

    auto ta1_z_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 100);

    auto ta1_z_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 101);

    auto ta1_z_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 102);

    auto ta1_z_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 103);

    auto ta1_z_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 104);

    auto ta1_z_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 105);

    auto ta1_z_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 106);

    auto ta1_z_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 107);

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd);

    auto ta1_x_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 1);

    auto ta1_x_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 2);

    auto ta1_x_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 3);

    auto ta1_x_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 4);

    auto ta1_x_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 5);

    auto ta1_x_xy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 6);

    auto ta1_x_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 7);

    auto ta1_x_xy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 8);

    auto ta1_x_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 9);

    auto ta1_x_xz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 12);

    auto ta1_x_xz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 13);

    auto ta1_x_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 14);

    auto ta1_x_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 17);

    auto ta1_x_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 18);

    auto ta1_x_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 19);

    auto ta1_x_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 20);

    auto ta1_x_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 21);

    auto ta1_x_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 22);

    auto ta1_x_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 23);

    auto ta1_x_yz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 26);

    auto ta1_x_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 28);

    auto ta1_x_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 29);

    auto ta1_x_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 30);

    auto ta1_x_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 31);

    auto ta1_x_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 32);

    auto ta1_x_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 33);

    auto ta1_x_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 34);

    auto ta1_x_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 35);

    auto ta1_y_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 36);

    auto ta1_y_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 37);

    auto ta1_y_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 38);

    auto ta1_y_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 39);

    auto ta1_y_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 40);

    auto ta1_y_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 41);

    auto ta1_y_xy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 42);

    auto ta1_y_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 43);

    auto ta1_y_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 45);

    auto ta1_y_xy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 46);

    auto ta1_y_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 50);

    auto ta1_y_xz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 52);

    auto ta1_y_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 53);

    auto ta1_y_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 54);

    auto ta1_y_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 55);

    auto ta1_y_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 56);

    auto ta1_y_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 57);

    auto ta1_y_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 58);

    auto ta1_y_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 59);

    auto ta1_y_yz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 61);

    auto ta1_y_yz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 63);

    auto ta1_y_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 64);

    auto ta1_y_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 65);

    auto ta1_y_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 66);

    auto ta1_y_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 67);

    auto ta1_y_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 68);

    auto ta1_y_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 69);

    auto ta1_y_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 70);

    auto ta1_y_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 71);

    auto ta1_z_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 72);

    auto ta1_z_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 73);

    auto ta1_z_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 74);

    auto ta1_z_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 75);

    auto ta1_z_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 76);

    auto ta1_z_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 77);

    auto ta1_z_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 79);

    auto ta1_z_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 81);

    auto ta1_z_xy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 82);

    auto ta1_z_xz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 84);

    auto ta1_z_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 86);

    auto ta1_z_xz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 88);

    auto ta1_z_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 89);

    auto ta1_z_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 90);

    auto ta1_z_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 91);

    auto ta1_z_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 92);

    auto ta1_z_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 93);

    auto ta1_z_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 94);

    auto ta1_z_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 95);

    auto ta1_z_yz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 98);

    auto ta1_z_yz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 99);

    auto ta1_z_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 100);

    auto ta1_z_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 101);

    auto ta1_z_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 102);

    auto ta1_z_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 103);

    auto ta1_z_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 104);

    auto ta1_z_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 105);

    auto ta1_z_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 106);

    auto ta1_z_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 107);

    // Set up 0-6 components of targeted buffer : FD

    auto ta1_x_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd);

    auto ta1_x_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 1);

    auto ta1_x_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 2);

    auto ta1_x_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 3);

    auto ta1_x_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 4);

    auto ta1_x_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 5);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xy_0,   \
                             ta1_x_x_xy_1,   \
                             ta1_x_x_xz_0,   \
                             ta1_x_x_xz_1,   \
                             ta1_x_x_yy_0,   \
                             ta1_x_x_yy_1,   \
                             ta1_x_x_yz_0,   \
                             ta1_x_x_yz_1,   \
                             ta1_x_x_zz_0,   \
                             ta1_x_x_zz_1,   \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_xx_0,  \
                             ta1_x_xx_xx_1,  \
                             ta1_x_xx_xy_0,  \
                             ta1_x_xx_xy_1,  \
                             ta1_x_xx_xz_0,  \
                             ta1_x_xx_xz_1,  \
                             ta1_x_xx_y_0,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xx_yy_0,  \
                             ta1_x_xx_yy_1,  \
                             ta1_x_xx_yz_0,  \
                             ta1_x_xx_yz_1,  \
                             ta1_x_xx_z_0,   \
                             ta1_x_xx_z_1,   \
                             ta1_x_xx_zz_0,  \
                             ta1_x_xx_zz_1,  \
                             ta1_x_xxx_xx_0, \
                             ta1_x_xxx_xy_0, \
                             ta1_x_xxx_xz_0, \
                             ta1_x_xxx_yy_0, \
                             ta1_x_xxx_yz_0, \
                             ta1_x_xxx_zz_0, \
                             ta_xx_xx_1,     \
                             ta_xx_xy_1,     \
                             ta_xx_xz_1,     \
                             ta_xx_yy_1,     \
                             ta_xx_yz_1,     \
                             ta_xx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_xx_0[i] = 2.0 * ta1_x_x_xx_0[i] * fe_0 - 2.0 * ta1_x_x_xx_1[i] * fe_0 + 2.0 * ta1_x_xx_x_0[i] * fe_0 -
                            2.0 * ta1_x_xx_x_1[i] * fe_0 + ta_xx_xx_1[i] + ta1_x_xx_xx_0[i] * pa_x[i] - ta1_x_xx_xx_1[i] * pc_x[i];

        ta1_x_xxx_xy_0[i] = 2.0 * ta1_x_x_xy_0[i] * fe_0 - 2.0 * ta1_x_x_xy_1[i] * fe_0 + ta1_x_xx_y_0[i] * fe_0 - ta1_x_xx_y_1[i] * fe_0 +
                            ta_xx_xy_1[i] + ta1_x_xx_xy_0[i] * pa_x[i] - ta1_x_xx_xy_1[i] * pc_x[i];

        ta1_x_xxx_xz_0[i] = 2.0 * ta1_x_x_xz_0[i] * fe_0 - 2.0 * ta1_x_x_xz_1[i] * fe_0 + ta1_x_xx_z_0[i] * fe_0 - ta1_x_xx_z_1[i] * fe_0 +
                            ta_xx_xz_1[i] + ta1_x_xx_xz_0[i] * pa_x[i] - ta1_x_xx_xz_1[i] * pc_x[i];

        ta1_x_xxx_yy_0[i] =
            2.0 * ta1_x_x_yy_0[i] * fe_0 - 2.0 * ta1_x_x_yy_1[i] * fe_0 + ta_xx_yy_1[i] + ta1_x_xx_yy_0[i] * pa_x[i] - ta1_x_xx_yy_1[i] * pc_x[i];

        ta1_x_xxx_yz_0[i] =
            2.0 * ta1_x_x_yz_0[i] * fe_0 - 2.0 * ta1_x_x_yz_1[i] * fe_0 + ta_xx_yz_1[i] + ta1_x_xx_yz_0[i] * pa_x[i] - ta1_x_xx_yz_1[i] * pc_x[i];

        ta1_x_xxx_zz_0[i] =
            2.0 * ta1_x_x_zz_0[i] * fe_0 - 2.0 * ta1_x_x_zz_1[i] * fe_0 + ta_xx_zz_1[i] + ta1_x_xx_zz_0[i] * pa_x[i] - ta1_x_xx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto ta1_x_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 6);

    auto ta1_x_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 7);

    auto ta1_x_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 8);

    auto ta1_x_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 9);

    auto ta1_x_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 10);

    auto ta1_x_xxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 11);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_xx_0,  \
                             ta1_x_xx_xx_1,  \
                             ta1_x_xx_xy_0,  \
                             ta1_x_xx_xy_1,  \
                             ta1_x_xx_xz_0,  \
                             ta1_x_xx_xz_1,  \
                             ta1_x_xx_y_0,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xx_yy_0,  \
                             ta1_x_xx_yy_1,  \
                             ta1_x_xx_yz_0,  \
                             ta1_x_xx_yz_1,  \
                             ta1_x_xx_z_0,   \
                             ta1_x_xx_z_1,   \
                             ta1_x_xx_zz_0,  \
                             ta1_x_xx_zz_1,  \
                             ta1_x_xxy_xx_0, \
                             ta1_x_xxy_xy_0, \
                             ta1_x_xxy_xz_0, \
                             ta1_x_xxy_yy_0, \
                             ta1_x_xxy_yz_0, \
                             ta1_x_xxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_xx_0[i] = ta1_x_xx_xx_0[i] * pa_y[i] - ta1_x_xx_xx_1[i] * pc_y[i];

        ta1_x_xxy_xy_0[i] = ta1_x_xx_x_0[i] * fe_0 - ta1_x_xx_x_1[i] * fe_0 + ta1_x_xx_xy_0[i] * pa_y[i] - ta1_x_xx_xy_1[i] * pc_y[i];

        ta1_x_xxy_xz_0[i] = ta1_x_xx_xz_0[i] * pa_y[i] - ta1_x_xx_xz_1[i] * pc_y[i];

        ta1_x_xxy_yy_0[i] = 2.0 * ta1_x_xx_y_0[i] * fe_0 - 2.0 * ta1_x_xx_y_1[i] * fe_0 + ta1_x_xx_yy_0[i] * pa_y[i] - ta1_x_xx_yy_1[i] * pc_y[i];

        ta1_x_xxy_yz_0[i] = ta1_x_xx_z_0[i] * fe_0 - ta1_x_xx_z_1[i] * fe_0 + ta1_x_xx_yz_0[i] * pa_y[i] - ta1_x_xx_yz_1[i] * pc_y[i];

        ta1_x_xxy_zz_0[i] = ta1_x_xx_zz_0[i] * pa_y[i] - ta1_x_xx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto ta1_x_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 12);

    auto ta1_x_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 13);

    auto ta1_x_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 14);

    auto ta1_x_xxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 15);

    auto ta1_x_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 16);

    auto ta1_x_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 17);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_xx_0,  \
                             ta1_x_xx_xx_1,  \
                             ta1_x_xx_xy_0,  \
                             ta1_x_xx_xy_1,  \
                             ta1_x_xx_xz_0,  \
                             ta1_x_xx_xz_1,  \
                             ta1_x_xx_y_0,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xx_yy_0,  \
                             ta1_x_xx_yy_1,  \
                             ta1_x_xx_yz_0,  \
                             ta1_x_xx_yz_1,  \
                             ta1_x_xx_z_0,   \
                             ta1_x_xx_z_1,   \
                             ta1_x_xx_zz_0,  \
                             ta1_x_xx_zz_1,  \
                             ta1_x_xxz_xx_0, \
                             ta1_x_xxz_xy_0, \
                             ta1_x_xxz_xz_0, \
                             ta1_x_xxz_yy_0, \
                             ta1_x_xxz_yz_0, \
                             ta1_x_xxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_xx_0[i] = ta1_x_xx_xx_0[i] * pa_z[i] - ta1_x_xx_xx_1[i] * pc_z[i];

        ta1_x_xxz_xy_0[i] = ta1_x_xx_xy_0[i] * pa_z[i] - ta1_x_xx_xy_1[i] * pc_z[i];

        ta1_x_xxz_xz_0[i] = ta1_x_xx_x_0[i] * fe_0 - ta1_x_xx_x_1[i] * fe_0 + ta1_x_xx_xz_0[i] * pa_z[i] - ta1_x_xx_xz_1[i] * pc_z[i];

        ta1_x_xxz_yy_0[i] = ta1_x_xx_yy_0[i] * pa_z[i] - ta1_x_xx_yy_1[i] * pc_z[i];

        ta1_x_xxz_yz_0[i] = ta1_x_xx_y_0[i] * fe_0 - ta1_x_xx_y_1[i] * fe_0 + ta1_x_xx_yz_0[i] * pa_z[i] - ta1_x_xx_yz_1[i] * pc_z[i];

        ta1_x_xxz_zz_0[i] = 2.0 * ta1_x_xx_z_0[i] * fe_0 - 2.0 * ta1_x_xx_z_1[i] * fe_0 + ta1_x_xx_zz_0[i] * pa_z[i] - ta1_x_xx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto ta1_x_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 18);

    auto ta1_x_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 19);

    auto ta1_x_xyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 20);

    auto ta1_x_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 21);

    auto ta1_x_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 22);

    auto ta1_x_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 23);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xz_0,   \
                             ta1_x_x_xz_1,   \
                             ta1_x_xy_xx_0,  \
                             ta1_x_xy_xx_1,  \
                             ta1_x_xy_xz_0,  \
                             ta1_x_xy_xz_1,  \
                             ta1_x_xyy_xx_0, \
                             ta1_x_xyy_xy_0, \
                             ta1_x_xyy_xz_0, \
                             ta1_x_xyy_yy_0, \
                             ta1_x_xyy_yz_0, \
                             ta1_x_xyy_zz_0, \
                             ta1_x_yy_xy_0,  \
                             ta1_x_yy_xy_1,  \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_yy_0,  \
                             ta1_x_yy_yy_1,  \
                             ta1_x_yy_yz_0,  \
                             ta1_x_yy_yz_1,  \
                             ta1_x_yy_zz_0,  \
                             ta1_x_yy_zz_1,  \
                             ta_yy_xy_1,     \
                             ta_yy_yy_1,     \
                             ta_yy_yz_1,     \
                             ta_yy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_xx_0[i] = ta1_x_x_xx_0[i] * fe_0 - ta1_x_x_xx_1[i] * fe_0 + ta1_x_xy_xx_0[i] * pa_y[i] - ta1_x_xy_xx_1[i] * pc_y[i];

        ta1_x_xyy_xy_0[i] = ta1_x_yy_y_0[i] * fe_0 - ta1_x_yy_y_1[i] * fe_0 + ta_yy_xy_1[i] + ta1_x_yy_xy_0[i] * pa_x[i] - ta1_x_yy_xy_1[i] * pc_x[i];

        ta1_x_xyy_xz_0[i] = ta1_x_x_xz_0[i] * fe_0 - ta1_x_x_xz_1[i] * fe_0 + ta1_x_xy_xz_0[i] * pa_y[i] - ta1_x_xy_xz_1[i] * pc_y[i];

        ta1_x_xyy_yy_0[i] = ta_yy_yy_1[i] + ta1_x_yy_yy_0[i] * pa_x[i] - ta1_x_yy_yy_1[i] * pc_x[i];

        ta1_x_xyy_yz_0[i] = ta_yy_yz_1[i] + ta1_x_yy_yz_0[i] * pa_x[i] - ta1_x_yy_yz_1[i] * pc_x[i];

        ta1_x_xyy_zz_0[i] = ta_yy_zz_1[i] + ta1_x_yy_zz_0[i] * pa_x[i] - ta1_x_yy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto ta1_x_xyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 24);

    auto ta1_x_xyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 25);

    auto ta1_x_xyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 26);

    auto ta1_x_xyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 27);

    auto ta1_x_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 28);

    auto ta1_x_xyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_xy_xy_0,  \
                             ta1_x_xy_xy_1,  \
                             ta1_x_xy_yy_0,  \
                             ta1_x_xy_yy_1,  \
                             ta1_x_xyz_xx_0, \
                             ta1_x_xyz_xy_0, \
                             ta1_x_xyz_xz_0, \
                             ta1_x_xyz_yy_0, \
                             ta1_x_xyz_yz_0, \
                             ta1_x_xyz_zz_0, \
                             ta1_x_xz_xx_0,  \
                             ta1_x_xz_xx_1,  \
                             ta1_x_xz_xz_0,  \
                             ta1_x_xz_xz_1,  \
                             ta1_x_xz_zz_0,  \
                             ta1_x_xz_zz_1,  \
                             ta1_x_yz_yz_0,  \
                             ta1_x_yz_yz_1,  \
                             ta_yz_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyz_xx_0[i] = ta1_x_xz_xx_0[i] * pa_y[i] - ta1_x_xz_xx_1[i] * pc_y[i];

        ta1_x_xyz_xy_0[i] = ta1_x_xy_xy_0[i] * pa_z[i] - ta1_x_xy_xy_1[i] * pc_z[i];

        ta1_x_xyz_xz_0[i] = ta1_x_xz_xz_0[i] * pa_y[i] - ta1_x_xz_xz_1[i] * pc_y[i];

        ta1_x_xyz_yy_0[i] = ta1_x_xy_yy_0[i] * pa_z[i] - ta1_x_xy_yy_1[i] * pc_z[i];

        ta1_x_xyz_yz_0[i] = ta_yz_yz_1[i] + ta1_x_yz_yz_0[i] * pa_x[i] - ta1_x_yz_yz_1[i] * pc_x[i];

        ta1_x_xyz_zz_0[i] = ta1_x_xz_zz_0[i] * pa_y[i] - ta1_x_xz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto ta1_x_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 30);

    auto ta1_x_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 31);

    auto ta1_x_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 32);

    auto ta1_x_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 33);

    auto ta1_x_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 34);

    auto ta1_x_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 35);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xy_0,   \
                             ta1_x_x_xy_1,   \
                             ta1_x_xz_xx_0,  \
                             ta1_x_xz_xx_1,  \
                             ta1_x_xz_xy_0,  \
                             ta1_x_xz_xy_1,  \
                             ta1_x_xzz_xx_0, \
                             ta1_x_xzz_xy_0, \
                             ta1_x_xzz_xz_0, \
                             ta1_x_xzz_yy_0, \
                             ta1_x_xzz_yz_0, \
                             ta1_x_xzz_zz_0, \
                             ta1_x_zz_xz_0,  \
                             ta1_x_zz_xz_1,  \
                             ta1_x_zz_yy_0,  \
                             ta1_x_zz_yy_1,  \
                             ta1_x_zz_yz_0,  \
                             ta1_x_zz_yz_1,  \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             ta1_x_zz_zz_0,  \
                             ta1_x_zz_zz_1,  \
                             ta_zz_xz_1,     \
                             ta_zz_yy_1,     \
                             ta_zz_yz_1,     \
                             ta_zz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_xx_0[i] = ta1_x_x_xx_0[i] * fe_0 - ta1_x_x_xx_1[i] * fe_0 + ta1_x_xz_xx_0[i] * pa_z[i] - ta1_x_xz_xx_1[i] * pc_z[i];

        ta1_x_xzz_xy_0[i] = ta1_x_x_xy_0[i] * fe_0 - ta1_x_x_xy_1[i] * fe_0 + ta1_x_xz_xy_0[i] * pa_z[i] - ta1_x_xz_xy_1[i] * pc_z[i];

        ta1_x_xzz_xz_0[i] = ta1_x_zz_z_0[i] * fe_0 - ta1_x_zz_z_1[i] * fe_0 + ta_zz_xz_1[i] + ta1_x_zz_xz_0[i] * pa_x[i] - ta1_x_zz_xz_1[i] * pc_x[i];

        ta1_x_xzz_yy_0[i] = ta_zz_yy_1[i] + ta1_x_zz_yy_0[i] * pa_x[i] - ta1_x_zz_yy_1[i] * pc_x[i];

        ta1_x_xzz_yz_0[i] = ta_zz_yz_1[i] + ta1_x_zz_yz_0[i] * pa_x[i] - ta1_x_zz_yz_1[i] * pc_x[i];

        ta1_x_xzz_zz_0[i] = ta_zz_zz_1[i] + ta1_x_zz_zz_0[i] * pa_x[i] - ta1_x_zz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto ta1_x_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 36);

    auto ta1_x_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 37);

    auto ta1_x_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 38);

    auto ta1_x_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 39);

    auto ta1_x_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 40);

    auto ta1_x_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 41);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_y_xx_0,   \
                             ta1_x_y_xx_1,   \
                             ta1_x_y_xy_0,   \
                             ta1_x_y_xy_1,   \
                             ta1_x_y_xz_0,   \
                             ta1_x_y_xz_1,   \
                             ta1_x_y_yy_0,   \
                             ta1_x_y_yy_1,   \
                             ta1_x_y_yz_0,   \
                             ta1_x_y_yz_1,   \
                             ta1_x_y_zz_0,   \
                             ta1_x_y_zz_1,   \
                             ta1_x_yy_x_0,   \
                             ta1_x_yy_x_1,   \
                             ta1_x_yy_xx_0,  \
                             ta1_x_yy_xx_1,  \
                             ta1_x_yy_xy_0,  \
                             ta1_x_yy_xy_1,  \
                             ta1_x_yy_xz_0,  \
                             ta1_x_yy_xz_1,  \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_yy_0,  \
                             ta1_x_yy_yy_1,  \
                             ta1_x_yy_yz_0,  \
                             ta1_x_yy_yz_1,  \
                             ta1_x_yy_z_0,   \
                             ta1_x_yy_z_1,   \
                             ta1_x_yy_zz_0,  \
                             ta1_x_yy_zz_1,  \
                             ta1_x_yyy_xx_0, \
                             ta1_x_yyy_xy_0, \
                             ta1_x_yyy_xz_0, \
                             ta1_x_yyy_yy_0, \
                             ta1_x_yyy_yz_0, \
                             ta1_x_yyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_xx_0[i] = 2.0 * ta1_x_y_xx_0[i] * fe_0 - 2.0 * ta1_x_y_xx_1[i] * fe_0 + ta1_x_yy_xx_0[i] * pa_y[i] - ta1_x_yy_xx_1[i] * pc_y[i];

        ta1_x_yyy_xy_0[i] = 2.0 * ta1_x_y_xy_0[i] * fe_0 - 2.0 * ta1_x_y_xy_1[i] * fe_0 + ta1_x_yy_x_0[i] * fe_0 - ta1_x_yy_x_1[i] * fe_0 +
                            ta1_x_yy_xy_0[i] * pa_y[i] - ta1_x_yy_xy_1[i] * pc_y[i];

        ta1_x_yyy_xz_0[i] = 2.0 * ta1_x_y_xz_0[i] * fe_0 - 2.0 * ta1_x_y_xz_1[i] * fe_0 + ta1_x_yy_xz_0[i] * pa_y[i] - ta1_x_yy_xz_1[i] * pc_y[i];

        ta1_x_yyy_yy_0[i] = 2.0 * ta1_x_y_yy_0[i] * fe_0 - 2.0 * ta1_x_y_yy_1[i] * fe_0 + 2.0 * ta1_x_yy_y_0[i] * fe_0 -
                            2.0 * ta1_x_yy_y_1[i] * fe_0 + ta1_x_yy_yy_0[i] * pa_y[i] - ta1_x_yy_yy_1[i] * pc_y[i];

        ta1_x_yyy_yz_0[i] = 2.0 * ta1_x_y_yz_0[i] * fe_0 - 2.0 * ta1_x_y_yz_1[i] * fe_0 + ta1_x_yy_z_0[i] * fe_0 - ta1_x_yy_z_1[i] * fe_0 +
                            ta1_x_yy_yz_0[i] * pa_y[i] - ta1_x_yy_yz_1[i] * pc_y[i];

        ta1_x_yyy_zz_0[i] = 2.0 * ta1_x_y_zz_0[i] * fe_0 - 2.0 * ta1_x_y_zz_1[i] * fe_0 + ta1_x_yy_zz_0[i] * pa_y[i] - ta1_x_yy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto ta1_x_yyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 42);

    auto ta1_x_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 43);

    auto ta1_x_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 44);

    auto ta1_x_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 45);

    auto ta1_x_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 46);

    auto ta1_x_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 47);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_yy_xx_0,  \
                             ta1_x_yy_xx_1,  \
                             ta1_x_yy_xy_0,  \
                             ta1_x_yy_xy_1,  \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_yy_0,  \
                             ta1_x_yy_yy_1,  \
                             ta1_x_yy_yz_0,  \
                             ta1_x_yy_yz_1,  \
                             ta1_x_yyz_xx_0, \
                             ta1_x_yyz_xy_0, \
                             ta1_x_yyz_xz_0, \
                             ta1_x_yyz_yy_0, \
                             ta1_x_yyz_yz_0, \
                             ta1_x_yyz_zz_0, \
                             ta1_x_yz_xz_0,  \
                             ta1_x_yz_xz_1,  \
                             ta1_x_yz_zz_0,  \
                             ta1_x_yz_zz_1,  \
                             ta1_x_z_xz_0,   \
                             ta1_x_z_xz_1,   \
                             ta1_x_z_zz_0,   \
                             ta1_x_z_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_xx_0[i] = ta1_x_yy_xx_0[i] * pa_z[i] - ta1_x_yy_xx_1[i] * pc_z[i];

        ta1_x_yyz_xy_0[i] = ta1_x_yy_xy_0[i] * pa_z[i] - ta1_x_yy_xy_1[i] * pc_z[i];

        ta1_x_yyz_xz_0[i] = ta1_x_z_xz_0[i] * fe_0 - ta1_x_z_xz_1[i] * fe_0 + ta1_x_yz_xz_0[i] * pa_y[i] - ta1_x_yz_xz_1[i] * pc_y[i];

        ta1_x_yyz_yy_0[i] = ta1_x_yy_yy_0[i] * pa_z[i] - ta1_x_yy_yy_1[i] * pc_z[i];

        ta1_x_yyz_yz_0[i] = ta1_x_yy_y_0[i] * fe_0 - ta1_x_yy_y_1[i] * fe_0 + ta1_x_yy_yz_0[i] * pa_z[i] - ta1_x_yy_yz_1[i] * pc_z[i];

        ta1_x_yyz_zz_0[i] = ta1_x_z_zz_0[i] * fe_0 - ta1_x_z_zz_1[i] * fe_0 + ta1_x_yz_zz_0[i] * pa_y[i] - ta1_x_yz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto ta1_x_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 48);

    auto ta1_x_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 49);

    auto ta1_x_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 50);

    auto ta1_x_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 51);

    auto ta1_x_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 52);

    auto ta1_x_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 53);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_yzz_xx_0, \
                             ta1_x_yzz_xy_0, \
                             ta1_x_yzz_xz_0, \
                             ta1_x_yzz_yy_0, \
                             ta1_x_yzz_yz_0, \
                             ta1_x_yzz_zz_0, \
                             ta1_x_zz_x_0,   \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_xx_0,  \
                             ta1_x_zz_xx_1,  \
                             ta1_x_zz_xy_0,  \
                             ta1_x_zz_xy_1,  \
                             ta1_x_zz_xz_0,  \
                             ta1_x_zz_xz_1,  \
                             ta1_x_zz_y_0,   \
                             ta1_x_zz_y_1,   \
                             ta1_x_zz_yy_0,  \
                             ta1_x_zz_yy_1,  \
                             ta1_x_zz_yz_0,  \
                             ta1_x_zz_yz_1,  \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             ta1_x_zz_zz_0,  \
                             ta1_x_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_xx_0[i] = ta1_x_zz_xx_0[i] * pa_y[i] - ta1_x_zz_xx_1[i] * pc_y[i];

        ta1_x_yzz_xy_0[i] = ta1_x_zz_x_0[i] * fe_0 - ta1_x_zz_x_1[i] * fe_0 + ta1_x_zz_xy_0[i] * pa_y[i] - ta1_x_zz_xy_1[i] * pc_y[i];

        ta1_x_yzz_xz_0[i] = ta1_x_zz_xz_0[i] * pa_y[i] - ta1_x_zz_xz_1[i] * pc_y[i];

        ta1_x_yzz_yy_0[i] = 2.0 * ta1_x_zz_y_0[i] * fe_0 - 2.0 * ta1_x_zz_y_1[i] * fe_0 + ta1_x_zz_yy_0[i] * pa_y[i] - ta1_x_zz_yy_1[i] * pc_y[i];

        ta1_x_yzz_yz_0[i] = ta1_x_zz_z_0[i] * fe_0 - ta1_x_zz_z_1[i] * fe_0 + ta1_x_zz_yz_0[i] * pa_y[i] - ta1_x_zz_yz_1[i] * pc_y[i];

        ta1_x_yzz_zz_0[i] = ta1_x_zz_zz_0[i] * pa_y[i] - ta1_x_zz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto ta1_x_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 54);

    auto ta1_x_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 55);

    auto ta1_x_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 56);

    auto ta1_x_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 57);

    auto ta1_x_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 58);

    auto ta1_x_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 59);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_z_xx_0,   \
                             ta1_x_z_xx_1,   \
                             ta1_x_z_xy_0,   \
                             ta1_x_z_xy_1,   \
                             ta1_x_z_xz_0,   \
                             ta1_x_z_xz_1,   \
                             ta1_x_z_yy_0,   \
                             ta1_x_z_yy_1,   \
                             ta1_x_z_yz_0,   \
                             ta1_x_z_yz_1,   \
                             ta1_x_z_zz_0,   \
                             ta1_x_z_zz_1,   \
                             ta1_x_zz_x_0,   \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_xx_0,  \
                             ta1_x_zz_xx_1,  \
                             ta1_x_zz_xy_0,  \
                             ta1_x_zz_xy_1,  \
                             ta1_x_zz_xz_0,  \
                             ta1_x_zz_xz_1,  \
                             ta1_x_zz_y_0,   \
                             ta1_x_zz_y_1,   \
                             ta1_x_zz_yy_0,  \
                             ta1_x_zz_yy_1,  \
                             ta1_x_zz_yz_0,  \
                             ta1_x_zz_yz_1,  \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             ta1_x_zz_zz_0,  \
                             ta1_x_zz_zz_1,  \
                             ta1_x_zzz_xx_0, \
                             ta1_x_zzz_xy_0, \
                             ta1_x_zzz_xz_0, \
                             ta1_x_zzz_yy_0, \
                             ta1_x_zzz_yz_0, \
                             ta1_x_zzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_xx_0[i] = 2.0 * ta1_x_z_xx_0[i] * fe_0 - 2.0 * ta1_x_z_xx_1[i] * fe_0 + ta1_x_zz_xx_0[i] * pa_z[i] - ta1_x_zz_xx_1[i] * pc_z[i];

        ta1_x_zzz_xy_0[i] = 2.0 * ta1_x_z_xy_0[i] * fe_0 - 2.0 * ta1_x_z_xy_1[i] * fe_0 + ta1_x_zz_xy_0[i] * pa_z[i] - ta1_x_zz_xy_1[i] * pc_z[i];

        ta1_x_zzz_xz_0[i] = 2.0 * ta1_x_z_xz_0[i] * fe_0 - 2.0 * ta1_x_z_xz_1[i] * fe_0 + ta1_x_zz_x_0[i] * fe_0 - ta1_x_zz_x_1[i] * fe_0 +
                            ta1_x_zz_xz_0[i] * pa_z[i] - ta1_x_zz_xz_1[i] * pc_z[i];

        ta1_x_zzz_yy_0[i] = 2.0 * ta1_x_z_yy_0[i] * fe_0 - 2.0 * ta1_x_z_yy_1[i] * fe_0 + ta1_x_zz_yy_0[i] * pa_z[i] - ta1_x_zz_yy_1[i] * pc_z[i];

        ta1_x_zzz_yz_0[i] = 2.0 * ta1_x_z_yz_0[i] * fe_0 - 2.0 * ta1_x_z_yz_1[i] * fe_0 + ta1_x_zz_y_0[i] * fe_0 - ta1_x_zz_y_1[i] * fe_0 +
                            ta1_x_zz_yz_0[i] * pa_z[i] - ta1_x_zz_yz_1[i] * pc_z[i];

        ta1_x_zzz_zz_0[i] = 2.0 * ta1_x_z_zz_0[i] * fe_0 - 2.0 * ta1_x_z_zz_1[i] * fe_0 + 2.0 * ta1_x_zz_z_0[i] * fe_0 -
                            2.0 * ta1_x_zz_z_1[i] * fe_0 + ta1_x_zz_zz_0[i] * pa_z[i] - ta1_x_zz_zz_1[i] * pc_z[i];
    }

    // Set up 60-66 components of targeted buffer : FD

    auto ta1_y_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 60);

    auto ta1_y_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 61);

    auto ta1_y_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 62);

    auto ta1_y_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 63);

    auto ta1_y_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 64);

    auto ta1_y_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 65);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_x_xx_0,   \
                             ta1_y_x_xx_1,   \
                             ta1_y_x_xy_0,   \
                             ta1_y_x_xy_1,   \
                             ta1_y_x_xz_0,   \
                             ta1_y_x_xz_1,   \
                             ta1_y_x_yy_0,   \
                             ta1_y_x_yy_1,   \
                             ta1_y_x_yz_0,   \
                             ta1_y_x_yz_1,   \
                             ta1_y_x_zz_0,   \
                             ta1_y_x_zz_1,   \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_xx_0,  \
                             ta1_y_xx_xx_1,  \
                             ta1_y_xx_xy_0,  \
                             ta1_y_xx_xy_1,  \
                             ta1_y_xx_xz_0,  \
                             ta1_y_xx_xz_1,  \
                             ta1_y_xx_y_0,   \
                             ta1_y_xx_y_1,   \
                             ta1_y_xx_yy_0,  \
                             ta1_y_xx_yy_1,  \
                             ta1_y_xx_yz_0,  \
                             ta1_y_xx_yz_1,  \
                             ta1_y_xx_z_0,   \
                             ta1_y_xx_z_1,   \
                             ta1_y_xx_zz_0,  \
                             ta1_y_xx_zz_1,  \
                             ta1_y_xxx_xx_0, \
                             ta1_y_xxx_xy_0, \
                             ta1_y_xxx_xz_0, \
                             ta1_y_xxx_yy_0, \
                             ta1_y_xxx_yz_0, \
                             ta1_y_xxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_xx_0[i] = 2.0 * ta1_y_x_xx_0[i] * fe_0 - 2.0 * ta1_y_x_xx_1[i] * fe_0 + 2.0 * ta1_y_xx_x_0[i] * fe_0 -
                            2.0 * ta1_y_xx_x_1[i] * fe_0 + ta1_y_xx_xx_0[i] * pa_x[i] - ta1_y_xx_xx_1[i] * pc_x[i];

        ta1_y_xxx_xy_0[i] = 2.0 * ta1_y_x_xy_0[i] * fe_0 - 2.0 * ta1_y_x_xy_1[i] * fe_0 + ta1_y_xx_y_0[i] * fe_0 - ta1_y_xx_y_1[i] * fe_0 +
                            ta1_y_xx_xy_0[i] * pa_x[i] - ta1_y_xx_xy_1[i] * pc_x[i];

        ta1_y_xxx_xz_0[i] = 2.0 * ta1_y_x_xz_0[i] * fe_0 - 2.0 * ta1_y_x_xz_1[i] * fe_0 + ta1_y_xx_z_0[i] * fe_0 - ta1_y_xx_z_1[i] * fe_0 +
                            ta1_y_xx_xz_0[i] * pa_x[i] - ta1_y_xx_xz_1[i] * pc_x[i];

        ta1_y_xxx_yy_0[i] = 2.0 * ta1_y_x_yy_0[i] * fe_0 - 2.0 * ta1_y_x_yy_1[i] * fe_0 + ta1_y_xx_yy_0[i] * pa_x[i] - ta1_y_xx_yy_1[i] * pc_x[i];

        ta1_y_xxx_yz_0[i] = 2.0 * ta1_y_x_yz_0[i] * fe_0 - 2.0 * ta1_y_x_yz_1[i] * fe_0 + ta1_y_xx_yz_0[i] * pa_x[i] - ta1_y_xx_yz_1[i] * pc_x[i];

        ta1_y_xxx_zz_0[i] = 2.0 * ta1_y_x_zz_0[i] * fe_0 - 2.0 * ta1_y_x_zz_1[i] * fe_0 + ta1_y_xx_zz_0[i] * pa_x[i] - ta1_y_xx_zz_1[i] * pc_x[i];
    }

    // Set up 66-72 components of targeted buffer : FD

    auto ta1_y_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 66);

    auto ta1_y_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 67);

    auto ta1_y_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 68);

    auto ta1_y_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 69);

    auto ta1_y_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 70);

    auto ta1_y_xxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 71);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_xx_0,  \
                             ta1_y_xx_xx_1,  \
                             ta1_y_xx_xy_0,  \
                             ta1_y_xx_xy_1,  \
                             ta1_y_xx_xz_0,  \
                             ta1_y_xx_xz_1,  \
                             ta1_y_xx_zz_0,  \
                             ta1_y_xx_zz_1,  \
                             ta1_y_xxy_xx_0, \
                             ta1_y_xxy_xy_0, \
                             ta1_y_xxy_xz_0, \
                             ta1_y_xxy_yy_0, \
                             ta1_y_xxy_yz_0, \
                             ta1_y_xxy_zz_0, \
                             ta1_y_xy_yy_0,  \
                             ta1_y_xy_yy_1,  \
                             ta1_y_xy_yz_0,  \
                             ta1_y_xy_yz_1,  \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_y_yz_0,   \
                             ta1_y_y_yz_1,   \
                             ta_xx_xx_1,     \
                             ta_xx_xy_1,     \
                             ta_xx_xz_1,     \
                             ta_xx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_xx_0[i] = ta_xx_xx_1[i] + ta1_y_xx_xx_0[i] * pa_y[i] - ta1_y_xx_xx_1[i] * pc_y[i];

        ta1_y_xxy_xy_0[i] = ta1_y_xx_x_0[i] * fe_0 - ta1_y_xx_x_1[i] * fe_0 + ta_xx_xy_1[i] + ta1_y_xx_xy_0[i] * pa_y[i] - ta1_y_xx_xy_1[i] * pc_y[i];

        ta1_y_xxy_xz_0[i] = ta_xx_xz_1[i] + ta1_y_xx_xz_0[i] * pa_y[i] - ta1_y_xx_xz_1[i] * pc_y[i];

        ta1_y_xxy_yy_0[i] = ta1_y_y_yy_0[i] * fe_0 - ta1_y_y_yy_1[i] * fe_0 + ta1_y_xy_yy_0[i] * pa_x[i] - ta1_y_xy_yy_1[i] * pc_x[i];

        ta1_y_xxy_yz_0[i] = ta1_y_y_yz_0[i] * fe_0 - ta1_y_y_yz_1[i] * fe_0 + ta1_y_xy_yz_0[i] * pa_x[i] - ta1_y_xy_yz_1[i] * pc_x[i];

        ta1_y_xxy_zz_0[i] = ta_xx_zz_1[i] + ta1_y_xx_zz_0[i] * pa_y[i] - ta1_y_xx_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : FD

    auto ta1_y_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 72);

    auto ta1_y_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 73);

    auto ta1_y_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 74);

    auto ta1_y_xxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 75);

    auto ta1_y_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 76);

    auto ta1_y_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 77);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_xx_0,  \
                             ta1_y_xx_xx_1,  \
                             ta1_y_xx_xy_0,  \
                             ta1_y_xx_xy_1,  \
                             ta1_y_xx_xz_0,  \
                             ta1_y_xx_xz_1,  \
                             ta1_y_xx_yy_0,  \
                             ta1_y_xx_yy_1,  \
                             ta1_y_xxz_xx_0, \
                             ta1_y_xxz_xy_0, \
                             ta1_y_xxz_xz_0, \
                             ta1_y_xxz_yy_0, \
                             ta1_y_xxz_yz_0, \
                             ta1_y_xxz_zz_0, \
                             ta1_y_xz_yz_0,  \
                             ta1_y_xz_yz_1,  \
                             ta1_y_xz_zz_0,  \
                             ta1_y_xz_zz_1,  \
                             ta1_y_z_yz_0,   \
                             ta1_y_z_yz_1,   \
                             ta1_y_z_zz_0,   \
                             ta1_y_z_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_xx_0[i] = ta1_y_xx_xx_0[i] * pa_z[i] - ta1_y_xx_xx_1[i] * pc_z[i];

        ta1_y_xxz_xy_0[i] = ta1_y_xx_xy_0[i] * pa_z[i] - ta1_y_xx_xy_1[i] * pc_z[i];

        ta1_y_xxz_xz_0[i] = ta1_y_xx_x_0[i] * fe_0 - ta1_y_xx_x_1[i] * fe_0 + ta1_y_xx_xz_0[i] * pa_z[i] - ta1_y_xx_xz_1[i] * pc_z[i];

        ta1_y_xxz_yy_0[i] = ta1_y_xx_yy_0[i] * pa_z[i] - ta1_y_xx_yy_1[i] * pc_z[i];

        ta1_y_xxz_yz_0[i] = ta1_y_z_yz_0[i] * fe_0 - ta1_y_z_yz_1[i] * fe_0 + ta1_y_xz_yz_0[i] * pa_x[i] - ta1_y_xz_yz_1[i] * pc_x[i];

        ta1_y_xxz_zz_0[i] = ta1_y_z_zz_0[i] * fe_0 - ta1_y_z_zz_1[i] * fe_0 + ta1_y_xz_zz_0[i] * pa_x[i] - ta1_y_xz_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : FD

    auto ta1_y_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 78);

    auto ta1_y_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 79);

    auto ta1_y_xyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 80);

    auto ta1_y_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 81);

    auto ta1_y_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 82);

    auto ta1_y_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 83);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xyy_xx_0, \
                             ta1_y_xyy_xy_0, \
                             ta1_y_xyy_xz_0, \
                             ta1_y_xyy_yy_0, \
                             ta1_y_xyy_yz_0, \
                             ta1_y_xyy_zz_0, \
                             ta1_y_yy_x_0,   \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_xx_0,  \
                             ta1_y_yy_xx_1,  \
                             ta1_y_yy_xy_0,  \
                             ta1_y_yy_xy_1,  \
                             ta1_y_yy_xz_0,  \
                             ta1_y_yy_xz_1,  \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_yy_0,  \
                             ta1_y_yy_yy_1,  \
                             ta1_y_yy_yz_0,  \
                             ta1_y_yy_yz_1,  \
                             ta1_y_yy_z_0,   \
                             ta1_y_yy_z_1,   \
                             ta1_y_yy_zz_0,  \
                             ta1_y_yy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_xx_0[i] = 2.0 * ta1_y_yy_x_0[i] * fe_0 - 2.0 * ta1_y_yy_x_1[i] * fe_0 + ta1_y_yy_xx_0[i] * pa_x[i] - ta1_y_yy_xx_1[i] * pc_x[i];

        ta1_y_xyy_xy_0[i] = ta1_y_yy_y_0[i] * fe_0 - ta1_y_yy_y_1[i] * fe_0 + ta1_y_yy_xy_0[i] * pa_x[i] - ta1_y_yy_xy_1[i] * pc_x[i];

        ta1_y_xyy_xz_0[i] = ta1_y_yy_z_0[i] * fe_0 - ta1_y_yy_z_1[i] * fe_0 + ta1_y_yy_xz_0[i] * pa_x[i] - ta1_y_yy_xz_1[i] * pc_x[i];

        ta1_y_xyy_yy_0[i] = ta1_y_yy_yy_0[i] * pa_x[i] - ta1_y_yy_yy_1[i] * pc_x[i];

        ta1_y_xyy_yz_0[i] = ta1_y_yy_yz_0[i] * pa_x[i] - ta1_y_yy_yz_1[i] * pc_x[i];

        ta1_y_xyy_zz_0[i] = ta1_y_yy_zz_0[i] * pa_x[i] - ta1_y_yy_zz_1[i] * pc_x[i];
    }

    // Set up 84-90 components of targeted buffer : FD

    auto ta1_y_xyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 84);

    auto ta1_y_xyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 85);

    auto ta1_y_xyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 86);

    auto ta1_y_xyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 87);

    auto ta1_y_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 88);

    auto ta1_y_xyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_xy_xx_0,  \
                             ta1_y_xy_xx_1,  \
                             ta1_y_xy_xy_0,  \
                             ta1_y_xy_xy_1,  \
                             ta1_y_xyz_xx_0, \
                             ta1_y_xyz_xy_0, \
                             ta1_y_xyz_xz_0, \
                             ta1_y_xyz_yy_0, \
                             ta1_y_xyz_yz_0, \
                             ta1_y_xyz_zz_0, \
                             ta1_y_xz_xz_0,  \
                             ta1_y_xz_xz_1,  \
                             ta1_y_yz_yy_0,  \
                             ta1_y_yz_yy_1,  \
                             ta1_y_yz_yz_0,  \
                             ta1_y_yz_yz_1,  \
                             ta1_y_yz_zz_0,  \
                             ta1_y_yz_zz_1,  \
                             ta_xz_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyz_xx_0[i] = ta1_y_xy_xx_0[i] * pa_z[i] - ta1_y_xy_xx_1[i] * pc_z[i];

        ta1_y_xyz_xy_0[i] = ta1_y_xy_xy_0[i] * pa_z[i] - ta1_y_xy_xy_1[i] * pc_z[i];

        ta1_y_xyz_xz_0[i] = ta_xz_xz_1[i] + ta1_y_xz_xz_0[i] * pa_y[i] - ta1_y_xz_xz_1[i] * pc_y[i];

        ta1_y_xyz_yy_0[i] = ta1_y_yz_yy_0[i] * pa_x[i] - ta1_y_yz_yy_1[i] * pc_x[i];

        ta1_y_xyz_yz_0[i] = ta1_y_yz_yz_0[i] * pa_x[i] - ta1_y_yz_yz_1[i] * pc_x[i];

        ta1_y_xyz_zz_0[i] = ta1_y_yz_zz_0[i] * pa_x[i] - ta1_y_yz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : FD

    auto ta1_y_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 90);

    auto ta1_y_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 91);

    auto ta1_y_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 92);

    auto ta1_y_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 93);

    auto ta1_y_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 94);

    auto ta1_y_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 95);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xzz_xx_0, \
                             ta1_y_xzz_xy_0, \
                             ta1_y_xzz_xz_0, \
                             ta1_y_xzz_yy_0, \
                             ta1_y_xzz_yz_0, \
                             ta1_y_xzz_zz_0, \
                             ta1_y_zz_x_0,   \
                             ta1_y_zz_x_1,   \
                             ta1_y_zz_xx_0,  \
                             ta1_y_zz_xx_1,  \
                             ta1_y_zz_xy_0,  \
                             ta1_y_zz_xy_1,  \
                             ta1_y_zz_xz_0,  \
                             ta1_y_zz_xz_1,  \
                             ta1_y_zz_y_0,   \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_yy_0,  \
                             ta1_y_zz_yy_1,  \
                             ta1_y_zz_yz_0,  \
                             ta1_y_zz_yz_1,  \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             ta1_y_zz_zz_0,  \
                             ta1_y_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_xx_0[i] = 2.0 * ta1_y_zz_x_0[i] * fe_0 - 2.0 * ta1_y_zz_x_1[i] * fe_0 + ta1_y_zz_xx_0[i] * pa_x[i] - ta1_y_zz_xx_1[i] * pc_x[i];

        ta1_y_xzz_xy_0[i] = ta1_y_zz_y_0[i] * fe_0 - ta1_y_zz_y_1[i] * fe_0 + ta1_y_zz_xy_0[i] * pa_x[i] - ta1_y_zz_xy_1[i] * pc_x[i];

        ta1_y_xzz_xz_0[i] = ta1_y_zz_z_0[i] * fe_0 - ta1_y_zz_z_1[i] * fe_0 + ta1_y_zz_xz_0[i] * pa_x[i] - ta1_y_zz_xz_1[i] * pc_x[i];

        ta1_y_xzz_yy_0[i] = ta1_y_zz_yy_0[i] * pa_x[i] - ta1_y_zz_yy_1[i] * pc_x[i];

        ta1_y_xzz_yz_0[i] = ta1_y_zz_yz_0[i] * pa_x[i] - ta1_y_zz_yz_1[i] * pc_x[i];

        ta1_y_xzz_zz_0[i] = ta1_y_zz_zz_0[i] * pa_x[i] - ta1_y_zz_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : FD

    auto ta1_y_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 96);

    auto ta1_y_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 97);

    auto ta1_y_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 98);

    auto ta1_y_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 99);

    auto ta1_y_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 100);

    auto ta1_y_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 101);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_y_xx_0,   \
                             ta1_y_y_xx_1,   \
                             ta1_y_y_xy_0,   \
                             ta1_y_y_xy_1,   \
                             ta1_y_y_xz_0,   \
                             ta1_y_y_xz_1,   \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_y_yz_0,   \
                             ta1_y_y_yz_1,   \
                             ta1_y_y_zz_0,   \
                             ta1_y_y_zz_1,   \
                             ta1_y_yy_x_0,   \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_xx_0,  \
                             ta1_y_yy_xx_1,  \
                             ta1_y_yy_xy_0,  \
                             ta1_y_yy_xy_1,  \
                             ta1_y_yy_xz_0,  \
                             ta1_y_yy_xz_1,  \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_yy_0,  \
                             ta1_y_yy_yy_1,  \
                             ta1_y_yy_yz_0,  \
                             ta1_y_yy_yz_1,  \
                             ta1_y_yy_z_0,   \
                             ta1_y_yy_z_1,   \
                             ta1_y_yy_zz_0,  \
                             ta1_y_yy_zz_1,  \
                             ta1_y_yyy_xx_0, \
                             ta1_y_yyy_xy_0, \
                             ta1_y_yyy_xz_0, \
                             ta1_y_yyy_yy_0, \
                             ta1_y_yyy_yz_0, \
                             ta1_y_yyy_zz_0, \
                             ta_yy_xx_1,     \
                             ta_yy_xy_1,     \
                             ta_yy_xz_1,     \
                             ta_yy_yy_1,     \
                             ta_yy_yz_1,     \
                             ta_yy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_xx_0[i] =
            2.0 * ta1_y_y_xx_0[i] * fe_0 - 2.0 * ta1_y_y_xx_1[i] * fe_0 + ta_yy_xx_1[i] + ta1_y_yy_xx_0[i] * pa_y[i] - ta1_y_yy_xx_1[i] * pc_y[i];

        ta1_y_yyy_xy_0[i] = 2.0 * ta1_y_y_xy_0[i] * fe_0 - 2.0 * ta1_y_y_xy_1[i] * fe_0 + ta1_y_yy_x_0[i] * fe_0 - ta1_y_yy_x_1[i] * fe_0 +
                            ta_yy_xy_1[i] + ta1_y_yy_xy_0[i] * pa_y[i] - ta1_y_yy_xy_1[i] * pc_y[i];

        ta1_y_yyy_xz_0[i] =
            2.0 * ta1_y_y_xz_0[i] * fe_0 - 2.0 * ta1_y_y_xz_1[i] * fe_0 + ta_yy_xz_1[i] + ta1_y_yy_xz_0[i] * pa_y[i] - ta1_y_yy_xz_1[i] * pc_y[i];

        ta1_y_yyy_yy_0[i] = 2.0 * ta1_y_y_yy_0[i] * fe_0 - 2.0 * ta1_y_y_yy_1[i] * fe_0 + 2.0 * ta1_y_yy_y_0[i] * fe_0 -
                            2.0 * ta1_y_yy_y_1[i] * fe_0 + ta_yy_yy_1[i] + ta1_y_yy_yy_0[i] * pa_y[i] - ta1_y_yy_yy_1[i] * pc_y[i];

        ta1_y_yyy_yz_0[i] = 2.0 * ta1_y_y_yz_0[i] * fe_0 - 2.0 * ta1_y_y_yz_1[i] * fe_0 + ta1_y_yy_z_0[i] * fe_0 - ta1_y_yy_z_1[i] * fe_0 +
                            ta_yy_yz_1[i] + ta1_y_yy_yz_0[i] * pa_y[i] - ta1_y_yy_yz_1[i] * pc_y[i];

        ta1_y_yyy_zz_0[i] =
            2.0 * ta1_y_y_zz_0[i] * fe_0 - 2.0 * ta1_y_y_zz_1[i] * fe_0 + ta_yy_zz_1[i] + ta1_y_yy_zz_0[i] * pa_y[i] - ta1_y_yy_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : FD

    auto ta1_y_yyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 102);

    auto ta1_y_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 103);

    auto ta1_y_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 104);

    auto ta1_y_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 105);

    auto ta1_y_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 106);

    auto ta1_y_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 107);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_yy_x_0,   \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_xx_0,  \
                             ta1_y_yy_xx_1,  \
                             ta1_y_yy_xy_0,  \
                             ta1_y_yy_xy_1,  \
                             ta1_y_yy_xz_0,  \
                             ta1_y_yy_xz_1,  \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_yy_0,  \
                             ta1_y_yy_yy_1,  \
                             ta1_y_yy_yz_0,  \
                             ta1_y_yy_yz_1,  \
                             ta1_y_yy_z_0,   \
                             ta1_y_yy_z_1,   \
                             ta1_y_yy_zz_0,  \
                             ta1_y_yy_zz_1,  \
                             ta1_y_yyz_xx_0, \
                             ta1_y_yyz_xy_0, \
                             ta1_y_yyz_xz_0, \
                             ta1_y_yyz_yy_0, \
                             ta1_y_yyz_yz_0, \
                             ta1_y_yyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_xx_0[i] = ta1_y_yy_xx_0[i] * pa_z[i] - ta1_y_yy_xx_1[i] * pc_z[i];

        ta1_y_yyz_xy_0[i] = ta1_y_yy_xy_0[i] * pa_z[i] - ta1_y_yy_xy_1[i] * pc_z[i];

        ta1_y_yyz_xz_0[i] = ta1_y_yy_x_0[i] * fe_0 - ta1_y_yy_x_1[i] * fe_0 + ta1_y_yy_xz_0[i] * pa_z[i] - ta1_y_yy_xz_1[i] * pc_z[i];

        ta1_y_yyz_yy_0[i] = ta1_y_yy_yy_0[i] * pa_z[i] - ta1_y_yy_yy_1[i] * pc_z[i];

        ta1_y_yyz_yz_0[i] = ta1_y_yy_y_0[i] * fe_0 - ta1_y_yy_y_1[i] * fe_0 + ta1_y_yy_yz_0[i] * pa_z[i] - ta1_y_yy_yz_1[i] * pc_z[i];

        ta1_y_yyz_zz_0[i] = 2.0 * ta1_y_yy_z_0[i] * fe_0 - 2.0 * ta1_y_yy_z_1[i] * fe_0 + ta1_y_yy_zz_0[i] * pa_z[i] - ta1_y_yy_zz_1[i] * pc_z[i];
    }

    // Set up 108-114 components of targeted buffer : FD

    auto ta1_y_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 108);

    auto ta1_y_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 109);

    auto ta1_y_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 110);

    auto ta1_y_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 111);

    auto ta1_y_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 112);

    auto ta1_y_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 113);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_y_xy_0,   \
                             ta1_y_y_xy_1,   \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_yz_xy_0,  \
                             ta1_y_yz_xy_1,  \
                             ta1_y_yz_yy_0,  \
                             ta1_y_yz_yy_1,  \
                             ta1_y_yzz_xx_0, \
                             ta1_y_yzz_xy_0, \
                             ta1_y_yzz_xz_0, \
                             ta1_y_yzz_yy_0, \
                             ta1_y_yzz_yz_0, \
                             ta1_y_yzz_zz_0, \
                             ta1_y_zz_xx_0,  \
                             ta1_y_zz_xx_1,  \
                             ta1_y_zz_xz_0,  \
                             ta1_y_zz_xz_1,  \
                             ta1_y_zz_yz_0,  \
                             ta1_y_zz_yz_1,  \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             ta1_y_zz_zz_0,  \
                             ta1_y_zz_zz_1,  \
                             ta_zz_xx_1,     \
                             ta_zz_xz_1,     \
                             ta_zz_yz_1,     \
                             ta_zz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_xx_0[i] = ta_zz_xx_1[i] + ta1_y_zz_xx_0[i] * pa_y[i] - ta1_y_zz_xx_1[i] * pc_y[i];

        ta1_y_yzz_xy_0[i] = ta1_y_y_xy_0[i] * fe_0 - ta1_y_y_xy_1[i] * fe_0 + ta1_y_yz_xy_0[i] * pa_z[i] - ta1_y_yz_xy_1[i] * pc_z[i];

        ta1_y_yzz_xz_0[i] = ta_zz_xz_1[i] + ta1_y_zz_xz_0[i] * pa_y[i] - ta1_y_zz_xz_1[i] * pc_y[i];

        ta1_y_yzz_yy_0[i] = ta1_y_y_yy_0[i] * fe_0 - ta1_y_y_yy_1[i] * fe_0 + ta1_y_yz_yy_0[i] * pa_z[i] - ta1_y_yz_yy_1[i] * pc_z[i];

        ta1_y_yzz_yz_0[i] = ta1_y_zz_z_0[i] * fe_0 - ta1_y_zz_z_1[i] * fe_0 + ta_zz_yz_1[i] + ta1_y_zz_yz_0[i] * pa_y[i] - ta1_y_zz_yz_1[i] * pc_y[i];

        ta1_y_yzz_zz_0[i] = ta_zz_zz_1[i] + ta1_y_zz_zz_0[i] * pa_y[i] - ta1_y_zz_zz_1[i] * pc_y[i];
    }

    // Set up 114-120 components of targeted buffer : FD

    auto ta1_y_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 114);

    auto ta1_y_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 115);

    auto ta1_y_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 116);

    auto ta1_y_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 117);

    auto ta1_y_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 118);

    auto ta1_y_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 119);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_z_xx_0,   \
                             ta1_y_z_xx_1,   \
                             ta1_y_z_xy_0,   \
                             ta1_y_z_xy_1,   \
                             ta1_y_z_xz_0,   \
                             ta1_y_z_xz_1,   \
                             ta1_y_z_yy_0,   \
                             ta1_y_z_yy_1,   \
                             ta1_y_z_yz_0,   \
                             ta1_y_z_yz_1,   \
                             ta1_y_z_zz_0,   \
                             ta1_y_z_zz_1,   \
                             ta1_y_zz_x_0,   \
                             ta1_y_zz_x_1,   \
                             ta1_y_zz_xx_0,  \
                             ta1_y_zz_xx_1,  \
                             ta1_y_zz_xy_0,  \
                             ta1_y_zz_xy_1,  \
                             ta1_y_zz_xz_0,  \
                             ta1_y_zz_xz_1,  \
                             ta1_y_zz_y_0,   \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_yy_0,  \
                             ta1_y_zz_yy_1,  \
                             ta1_y_zz_yz_0,  \
                             ta1_y_zz_yz_1,  \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             ta1_y_zz_zz_0,  \
                             ta1_y_zz_zz_1,  \
                             ta1_y_zzz_xx_0, \
                             ta1_y_zzz_xy_0, \
                             ta1_y_zzz_xz_0, \
                             ta1_y_zzz_yy_0, \
                             ta1_y_zzz_yz_0, \
                             ta1_y_zzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_xx_0[i] = 2.0 * ta1_y_z_xx_0[i] * fe_0 - 2.0 * ta1_y_z_xx_1[i] * fe_0 + ta1_y_zz_xx_0[i] * pa_z[i] - ta1_y_zz_xx_1[i] * pc_z[i];

        ta1_y_zzz_xy_0[i] = 2.0 * ta1_y_z_xy_0[i] * fe_0 - 2.0 * ta1_y_z_xy_1[i] * fe_0 + ta1_y_zz_xy_0[i] * pa_z[i] - ta1_y_zz_xy_1[i] * pc_z[i];

        ta1_y_zzz_xz_0[i] = 2.0 * ta1_y_z_xz_0[i] * fe_0 - 2.0 * ta1_y_z_xz_1[i] * fe_0 + ta1_y_zz_x_0[i] * fe_0 - ta1_y_zz_x_1[i] * fe_0 +
                            ta1_y_zz_xz_0[i] * pa_z[i] - ta1_y_zz_xz_1[i] * pc_z[i];

        ta1_y_zzz_yy_0[i] = 2.0 * ta1_y_z_yy_0[i] * fe_0 - 2.0 * ta1_y_z_yy_1[i] * fe_0 + ta1_y_zz_yy_0[i] * pa_z[i] - ta1_y_zz_yy_1[i] * pc_z[i];

        ta1_y_zzz_yz_0[i] = 2.0 * ta1_y_z_yz_0[i] * fe_0 - 2.0 * ta1_y_z_yz_1[i] * fe_0 + ta1_y_zz_y_0[i] * fe_0 - ta1_y_zz_y_1[i] * fe_0 +
                            ta1_y_zz_yz_0[i] * pa_z[i] - ta1_y_zz_yz_1[i] * pc_z[i];

        ta1_y_zzz_zz_0[i] = 2.0 * ta1_y_z_zz_0[i] * fe_0 - 2.0 * ta1_y_z_zz_1[i] * fe_0 + 2.0 * ta1_y_zz_z_0[i] * fe_0 -
                            2.0 * ta1_y_zz_z_1[i] * fe_0 + ta1_y_zz_zz_0[i] * pa_z[i] - ta1_y_zz_zz_1[i] * pc_z[i];
    }

    // Set up 120-126 components of targeted buffer : FD

    auto ta1_z_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 120);

    auto ta1_z_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 121);

    auto ta1_z_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 122);

    auto ta1_z_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 123);

    auto ta1_z_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 124);

    auto ta1_z_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 125);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_x_xx_0,   \
                             ta1_z_x_xx_1,   \
                             ta1_z_x_xy_0,   \
                             ta1_z_x_xy_1,   \
                             ta1_z_x_xz_0,   \
                             ta1_z_x_xz_1,   \
                             ta1_z_x_yy_0,   \
                             ta1_z_x_yy_1,   \
                             ta1_z_x_yz_0,   \
                             ta1_z_x_yz_1,   \
                             ta1_z_x_zz_0,   \
                             ta1_z_x_zz_1,   \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_xx_0,  \
                             ta1_z_xx_xx_1,  \
                             ta1_z_xx_xy_0,  \
                             ta1_z_xx_xy_1,  \
                             ta1_z_xx_xz_0,  \
                             ta1_z_xx_xz_1,  \
                             ta1_z_xx_y_0,   \
                             ta1_z_xx_y_1,   \
                             ta1_z_xx_yy_0,  \
                             ta1_z_xx_yy_1,  \
                             ta1_z_xx_yz_0,  \
                             ta1_z_xx_yz_1,  \
                             ta1_z_xx_z_0,   \
                             ta1_z_xx_z_1,   \
                             ta1_z_xx_zz_0,  \
                             ta1_z_xx_zz_1,  \
                             ta1_z_xxx_xx_0, \
                             ta1_z_xxx_xy_0, \
                             ta1_z_xxx_xz_0, \
                             ta1_z_xxx_yy_0, \
                             ta1_z_xxx_yz_0, \
                             ta1_z_xxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_xx_0[i] = 2.0 * ta1_z_x_xx_0[i] * fe_0 - 2.0 * ta1_z_x_xx_1[i] * fe_0 + 2.0 * ta1_z_xx_x_0[i] * fe_0 -
                            2.0 * ta1_z_xx_x_1[i] * fe_0 + ta1_z_xx_xx_0[i] * pa_x[i] - ta1_z_xx_xx_1[i] * pc_x[i];

        ta1_z_xxx_xy_0[i] = 2.0 * ta1_z_x_xy_0[i] * fe_0 - 2.0 * ta1_z_x_xy_1[i] * fe_0 + ta1_z_xx_y_0[i] * fe_0 - ta1_z_xx_y_1[i] * fe_0 +
                            ta1_z_xx_xy_0[i] * pa_x[i] - ta1_z_xx_xy_1[i] * pc_x[i];

        ta1_z_xxx_xz_0[i] = 2.0 * ta1_z_x_xz_0[i] * fe_0 - 2.0 * ta1_z_x_xz_1[i] * fe_0 + ta1_z_xx_z_0[i] * fe_0 - ta1_z_xx_z_1[i] * fe_0 +
                            ta1_z_xx_xz_0[i] * pa_x[i] - ta1_z_xx_xz_1[i] * pc_x[i];

        ta1_z_xxx_yy_0[i] = 2.0 * ta1_z_x_yy_0[i] * fe_0 - 2.0 * ta1_z_x_yy_1[i] * fe_0 + ta1_z_xx_yy_0[i] * pa_x[i] - ta1_z_xx_yy_1[i] * pc_x[i];

        ta1_z_xxx_yz_0[i] = 2.0 * ta1_z_x_yz_0[i] * fe_0 - 2.0 * ta1_z_x_yz_1[i] * fe_0 + ta1_z_xx_yz_0[i] * pa_x[i] - ta1_z_xx_yz_1[i] * pc_x[i];

        ta1_z_xxx_zz_0[i] = 2.0 * ta1_z_x_zz_0[i] * fe_0 - 2.0 * ta1_z_x_zz_1[i] * fe_0 + ta1_z_xx_zz_0[i] * pa_x[i] - ta1_z_xx_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : FD

    auto ta1_z_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 126);

    auto ta1_z_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 127);

    auto ta1_z_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 128);

    auto ta1_z_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 129);

    auto ta1_z_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 130);

    auto ta1_z_xxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 131);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_xx_0,  \
                             ta1_z_xx_xx_1,  \
                             ta1_z_xx_xy_0,  \
                             ta1_z_xx_xy_1,  \
                             ta1_z_xx_xz_0,  \
                             ta1_z_xx_xz_1,  \
                             ta1_z_xx_zz_0,  \
                             ta1_z_xx_zz_1,  \
                             ta1_z_xxy_xx_0, \
                             ta1_z_xxy_xy_0, \
                             ta1_z_xxy_xz_0, \
                             ta1_z_xxy_yy_0, \
                             ta1_z_xxy_yz_0, \
                             ta1_z_xxy_zz_0, \
                             ta1_z_xy_yy_0,  \
                             ta1_z_xy_yy_1,  \
                             ta1_z_xy_yz_0,  \
                             ta1_z_xy_yz_1,  \
                             ta1_z_y_yy_0,   \
                             ta1_z_y_yy_1,   \
                             ta1_z_y_yz_0,   \
                             ta1_z_y_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_xx_0[i] = ta1_z_xx_xx_0[i] * pa_y[i] - ta1_z_xx_xx_1[i] * pc_y[i];

        ta1_z_xxy_xy_0[i] = ta1_z_xx_x_0[i] * fe_0 - ta1_z_xx_x_1[i] * fe_0 + ta1_z_xx_xy_0[i] * pa_y[i] - ta1_z_xx_xy_1[i] * pc_y[i];

        ta1_z_xxy_xz_0[i] = ta1_z_xx_xz_0[i] * pa_y[i] - ta1_z_xx_xz_1[i] * pc_y[i];

        ta1_z_xxy_yy_0[i] = ta1_z_y_yy_0[i] * fe_0 - ta1_z_y_yy_1[i] * fe_0 + ta1_z_xy_yy_0[i] * pa_x[i] - ta1_z_xy_yy_1[i] * pc_x[i];

        ta1_z_xxy_yz_0[i] = ta1_z_y_yz_0[i] * fe_0 - ta1_z_y_yz_1[i] * fe_0 + ta1_z_xy_yz_0[i] * pa_x[i] - ta1_z_xy_yz_1[i] * pc_x[i];

        ta1_z_xxy_zz_0[i] = ta1_z_xx_zz_0[i] * pa_y[i] - ta1_z_xx_zz_1[i] * pc_y[i];
    }

    // Set up 132-138 components of targeted buffer : FD

    auto ta1_z_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 132);

    auto ta1_z_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 133);

    auto ta1_z_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 134);

    auto ta1_z_xxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 135);

    auto ta1_z_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 136);

    auto ta1_z_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 137);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_xx_0,  \
                             ta1_z_xx_xx_1,  \
                             ta1_z_xx_xy_0,  \
                             ta1_z_xx_xy_1,  \
                             ta1_z_xx_xz_0,  \
                             ta1_z_xx_xz_1,  \
                             ta1_z_xx_yy_0,  \
                             ta1_z_xx_yy_1,  \
                             ta1_z_xxz_xx_0, \
                             ta1_z_xxz_xy_0, \
                             ta1_z_xxz_xz_0, \
                             ta1_z_xxz_yy_0, \
                             ta1_z_xxz_yz_0, \
                             ta1_z_xxz_zz_0, \
                             ta1_z_xz_yz_0,  \
                             ta1_z_xz_yz_1,  \
                             ta1_z_xz_zz_0,  \
                             ta1_z_xz_zz_1,  \
                             ta1_z_z_yz_0,   \
                             ta1_z_z_yz_1,   \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta_xx_xx_1,     \
                             ta_xx_xy_1,     \
                             ta_xx_xz_1,     \
                             ta_xx_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_xx_0[i] = ta_xx_xx_1[i] + ta1_z_xx_xx_0[i] * pa_z[i] - ta1_z_xx_xx_1[i] * pc_z[i];

        ta1_z_xxz_xy_0[i] = ta_xx_xy_1[i] + ta1_z_xx_xy_0[i] * pa_z[i] - ta1_z_xx_xy_1[i] * pc_z[i];

        ta1_z_xxz_xz_0[i] = ta1_z_xx_x_0[i] * fe_0 - ta1_z_xx_x_1[i] * fe_0 + ta_xx_xz_1[i] + ta1_z_xx_xz_0[i] * pa_z[i] - ta1_z_xx_xz_1[i] * pc_z[i];

        ta1_z_xxz_yy_0[i] = ta_xx_yy_1[i] + ta1_z_xx_yy_0[i] * pa_z[i] - ta1_z_xx_yy_1[i] * pc_z[i];

        ta1_z_xxz_yz_0[i] = ta1_z_z_yz_0[i] * fe_0 - ta1_z_z_yz_1[i] * fe_0 + ta1_z_xz_yz_0[i] * pa_x[i] - ta1_z_xz_yz_1[i] * pc_x[i];

        ta1_z_xxz_zz_0[i] = ta1_z_z_zz_0[i] * fe_0 - ta1_z_z_zz_1[i] * fe_0 + ta1_z_xz_zz_0[i] * pa_x[i] - ta1_z_xz_zz_1[i] * pc_x[i];
    }

    // Set up 138-144 components of targeted buffer : FD

    auto ta1_z_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 138);

    auto ta1_z_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 139);

    auto ta1_z_xyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 140);

    auto ta1_z_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 141);

    auto ta1_z_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 142);

    auto ta1_z_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 143);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xyy_xx_0, \
                             ta1_z_xyy_xy_0, \
                             ta1_z_xyy_xz_0, \
                             ta1_z_xyy_yy_0, \
                             ta1_z_xyy_yz_0, \
                             ta1_z_xyy_zz_0, \
                             ta1_z_yy_x_0,   \
                             ta1_z_yy_x_1,   \
                             ta1_z_yy_xx_0,  \
                             ta1_z_yy_xx_1,  \
                             ta1_z_yy_xy_0,  \
                             ta1_z_yy_xy_1,  \
                             ta1_z_yy_xz_0,  \
                             ta1_z_yy_xz_1,  \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_yy_0,  \
                             ta1_z_yy_yy_1,  \
                             ta1_z_yy_yz_0,  \
                             ta1_z_yy_yz_1,  \
                             ta1_z_yy_z_0,   \
                             ta1_z_yy_z_1,   \
                             ta1_z_yy_zz_0,  \
                             ta1_z_yy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_xx_0[i] = 2.0 * ta1_z_yy_x_0[i] * fe_0 - 2.0 * ta1_z_yy_x_1[i] * fe_0 + ta1_z_yy_xx_0[i] * pa_x[i] - ta1_z_yy_xx_1[i] * pc_x[i];

        ta1_z_xyy_xy_0[i] = ta1_z_yy_y_0[i] * fe_0 - ta1_z_yy_y_1[i] * fe_0 + ta1_z_yy_xy_0[i] * pa_x[i] - ta1_z_yy_xy_1[i] * pc_x[i];

        ta1_z_xyy_xz_0[i] = ta1_z_yy_z_0[i] * fe_0 - ta1_z_yy_z_1[i] * fe_0 + ta1_z_yy_xz_0[i] * pa_x[i] - ta1_z_yy_xz_1[i] * pc_x[i];

        ta1_z_xyy_yy_0[i] = ta1_z_yy_yy_0[i] * pa_x[i] - ta1_z_yy_yy_1[i] * pc_x[i];

        ta1_z_xyy_yz_0[i] = ta1_z_yy_yz_0[i] * pa_x[i] - ta1_z_yy_yz_1[i] * pc_x[i];

        ta1_z_xyy_zz_0[i] = ta1_z_yy_zz_0[i] * pa_x[i] - ta1_z_yy_zz_1[i] * pc_x[i];
    }

    // Set up 144-150 components of targeted buffer : FD

    auto ta1_z_xyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 144);

    auto ta1_z_xyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 145);

    auto ta1_z_xyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 146);

    auto ta1_z_xyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 147);

    auto ta1_z_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 148);

    auto ta1_z_xyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 149);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_xy_xy_0,  \
                             ta1_z_xy_xy_1,  \
                             ta1_z_xyz_xx_0, \
                             ta1_z_xyz_xy_0, \
                             ta1_z_xyz_xz_0, \
                             ta1_z_xyz_yy_0, \
                             ta1_z_xyz_yz_0, \
                             ta1_z_xyz_zz_0, \
                             ta1_z_xz_xx_0,  \
                             ta1_z_xz_xx_1,  \
                             ta1_z_xz_xz_0,  \
                             ta1_z_xz_xz_1,  \
                             ta1_z_yz_yy_0,  \
                             ta1_z_yz_yy_1,  \
                             ta1_z_yz_yz_0,  \
                             ta1_z_yz_yz_1,  \
                             ta1_z_yz_zz_0,  \
                             ta1_z_yz_zz_1,  \
                             ta_xy_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyz_xx_0[i] = ta1_z_xz_xx_0[i] * pa_y[i] - ta1_z_xz_xx_1[i] * pc_y[i];

        ta1_z_xyz_xy_0[i] = ta_xy_xy_1[i] + ta1_z_xy_xy_0[i] * pa_z[i] - ta1_z_xy_xy_1[i] * pc_z[i];

        ta1_z_xyz_xz_0[i] = ta1_z_xz_xz_0[i] * pa_y[i] - ta1_z_xz_xz_1[i] * pc_y[i];

        ta1_z_xyz_yy_0[i] = ta1_z_yz_yy_0[i] * pa_x[i] - ta1_z_yz_yy_1[i] * pc_x[i];

        ta1_z_xyz_yz_0[i] = ta1_z_yz_yz_0[i] * pa_x[i] - ta1_z_yz_yz_1[i] * pc_x[i];

        ta1_z_xyz_zz_0[i] = ta1_z_yz_zz_0[i] * pa_x[i] - ta1_z_yz_zz_1[i] * pc_x[i];
    }

    // Set up 150-156 components of targeted buffer : FD

    auto ta1_z_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 150);

    auto ta1_z_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 151);

    auto ta1_z_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 152);

    auto ta1_z_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 153);

    auto ta1_z_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 154);

    auto ta1_z_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 155);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xzz_xx_0, \
                             ta1_z_xzz_xy_0, \
                             ta1_z_xzz_xz_0, \
                             ta1_z_xzz_yy_0, \
                             ta1_z_xzz_yz_0, \
                             ta1_z_xzz_zz_0, \
                             ta1_z_zz_x_0,   \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_xx_0,  \
                             ta1_z_zz_xx_1,  \
                             ta1_z_zz_xy_0,  \
                             ta1_z_zz_xy_1,  \
                             ta1_z_zz_xz_0,  \
                             ta1_z_zz_xz_1,  \
                             ta1_z_zz_y_0,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_yy_0,  \
                             ta1_z_zz_yy_1,  \
                             ta1_z_zz_yz_0,  \
                             ta1_z_zz_yz_1,  \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta1_z_zz_zz_0,  \
                             ta1_z_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_xx_0[i] = 2.0 * ta1_z_zz_x_0[i] * fe_0 - 2.0 * ta1_z_zz_x_1[i] * fe_0 + ta1_z_zz_xx_0[i] * pa_x[i] - ta1_z_zz_xx_1[i] * pc_x[i];

        ta1_z_xzz_xy_0[i] = ta1_z_zz_y_0[i] * fe_0 - ta1_z_zz_y_1[i] * fe_0 + ta1_z_zz_xy_0[i] * pa_x[i] - ta1_z_zz_xy_1[i] * pc_x[i];

        ta1_z_xzz_xz_0[i] = ta1_z_zz_z_0[i] * fe_0 - ta1_z_zz_z_1[i] * fe_0 + ta1_z_zz_xz_0[i] * pa_x[i] - ta1_z_zz_xz_1[i] * pc_x[i];

        ta1_z_xzz_yy_0[i] = ta1_z_zz_yy_0[i] * pa_x[i] - ta1_z_zz_yy_1[i] * pc_x[i];

        ta1_z_xzz_yz_0[i] = ta1_z_zz_yz_0[i] * pa_x[i] - ta1_z_zz_yz_1[i] * pc_x[i];

        ta1_z_xzz_zz_0[i] = ta1_z_zz_zz_0[i] * pa_x[i] - ta1_z_zz_zz_1[i] * pc_x[i];
    }

    // Set up 156-162 components of targeted buffer : FD

    auto ta1_z_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 156);

    auto ta1_z_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 157);

    auto ta1_z_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 158);

    auto ta1_z_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 159);

    auto ta1_z_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 160);

    auto ta1_z_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 161);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_y_xx_0,   \
                             ta1_z_y_xx_1,   \
                             ta1_z_y_xy_0,   \
                             ta1_z_y_xy_1,   \
                             ta1_z_y_xz_0,   \
                             ta1_z_y_xz_1,   \
                             ta1_z_y_yy_0,   \
                             ta1_z_y_yy_1,   \
                             ta1_z_y_yz_0,   \
                             ta1_z_y_yz_1,   \
                             ta1_z_y_zz_0,   \
                             ta1_z_y_zz_1,   \
                             ta1_z_yy_x_0,   \
                             ta1_z_yy_x_1,   \
                             ta1_z_yy_xx_0,  \
                             ta1_z_yy_xx_1,  \
                             ta1_z_yy_xy_0,  \
                             ta1_z_yy_xy_1,  \
                             ta1_z_yy_xz_0,  \
                             ta1_z_yy_xz_1,  \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_yy_0,  \
                             ta1_z_yy_yy_1,  \
                             ta1_z_yy_yz_0,  \
                             ta1_z_yy_yz_1,  \
                             ta1_z_yy_z_0,   \
                             ta1_z_yy_z_1,   \
                             ta1_z_yy_zz_0,  \
                             ta1_z_yy_zz_1,  \
                             ta1_z_yyy_xx_0, \
                             ta1_z_yyy_xy_0, \
                             ta1_z_yyy_xz_0, \
                             ta1_z_yyy_yy_0, \
                             ta1_z_yyy_yz_0, \
                             ta1_z_yyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_xx_0[i] = 2.0 * ta1_z_y_xx_0[i] * fe_0 - 2.0 * ta1_z_y_xx_1[i] * fe_0 + ta1_z_yy_xx_0[i] * pa_y[i] - ta1_z_yy_xx_1[i] * pc_y[i];

        ta1_z_yyy_xy_0[i] = 2.0 * ta1_z_y_xy_0[i] * fe_0 - 2.0 * ta1_z_y_xy_1[i] * fe_0 + ta1_z_yy_x_0[i] * fe_0 - ta1_z_yy_x_1[i] * fe_0 +
                            ta1_z_yy_xy_0[i] * pa_y[i] - ta1_z_yy_xy_1[i] * pc_y[i];

        ta1_z_yyy_xz_0[i] = 2.0 * ta1_z_y_xz_0[i] * fe_0 - 2.0 * ta1_z_y_xz_1[i] * fe_0 + ta1_z_yy_xz_0[i] * pa_y[i] - ta1_z_yy_xz_1[i] * pc_y[i];

        ta1_z_yyy_yy_0[i] = 2.0 * ta1_z_y_yy_0[i] * fe_0 - 2.0 * ta1_z_y_yy_1[i] * fe_0 + 2.0 * ta1_z_yy_y_0[i] * fe_0 -
                            2.0 * ta1_z_yy_y_1[i] * fe_0 + ta1_z_yy_yy_0[i] * pa_y[i] - ta1_z_yy_yy_1[i] * pc_y[i];

        ta1_z_yyy_yz_0[i] = 2.0 * ta1_z_y_yz_0[i] * fe_0 - 2.0 * ta1_z_y_yz_1[i] * fe_0 + ta1_z_yy_z_0[i] * fe_0 - ta1_z_yy_z_1[i] * fe_0 +
                            ta1_z_yy_yz_0[i] * pa_y[i] - ta1_z_yy_yz_1[i] * pc_y[i];

        ta1_z_yyy_zz_0[i] = 2.0 * ta1_z_y_zz_0[i] * fe_0 - 2.0 * ta1_z_y_zz_1[i] * fe_0 + ta1_z_yy_zz_0[i] * pa_y[i] - ta1_z_yy_zz_1[i] * pc_y[i];
    }

    // Set up 162-168 components of targeted buffer : FD

    auto ta1_z_yyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 162);

    auto ta1_z_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 163);

    auto ta1_z_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 164);

    auto ta1_z_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 165);

    auto ta1_z_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 166);

    auto ta1_z_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 167);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_yy_xx_0,  \
                             ta1_z_yy_xx_1,  \
                             ta1_z_yy_xy_0,  \
                             ta1_z_yy_xy_1,  \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_yy_0,  \
                             ta1_z_yy_yy_1,  \
                             ta1_z_yy_yz_0,  \
                             ta1_z_yy_yz_1,  \
                             ta1_z_yyz_xx_0, \
                             ta1_z_yyz_xy_0, \
                             ta1_z_yyz_xz_0, \
                             ta1_z_yyz_yy_0, \
                             ta1_z_yyz_yz_0, \
                             ta1_z_yyz_zz_0, \
                             ta1_z_yz_xz_0,  \
                             ta1_z_yz_xz_1,  \
                             ta1_z_yz_zz_0,  \
                             ta1_z_yz_zz_1,  \
                             ta1_z_z_xz_0,   \
                             ta1_z_z_xz_1,   \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta_yy_xx_1,     \
                             ta_yy_xy_1,     \
                             ta_yy_yy_1,     \
                             ta_yy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_xx_0[i] = ta_yy_xx_1[i] + ta1_z_yy_xx_0[i] * pa_z[i] - ta1_z_yy_xx_1[i] * pc_z[i];

        ta1_z_yyz_xy_0[i] = ta_yy_xy_1[i] + ta1_z_yy_xy_0[i] * pa_z[i] - ta1_z_yy_xy_1[i] * pc_z[i];

        ta1_z_yyz_xz_0[i] = ta1_z_z_xz_0[i] * fe_0 - ta1_z_z_xz_1[i] * fe_0 + ta1_z_yz_xz_0[i] * pa_y[i] - ta1_z_yz_xz_1[i] * pc_y[i];

        ta1_z_yyz_yy_0[i] = ta_yy_yy_1[i] + ta1_z_yy_yy_0[i] * pa_z[i] - ta1_z_yy_yy_1[i] * pc_z[i];

        ta1_z_yyz_yz_0[i] = ta1_z_yy_y_0[i] * fe_0 - ta1_z_yy_y_1[i] * fe_0 + ta_yy_yz_1[i] + ta1_z_yy_yz_0[i] * pa_z[i] - ta1_z_yy_yz_1[i] * pc_z[i];

        ta1_z_yyz_zz_0[i] = ta1_z_z_zz_0[i] * fe_0 - ta1_z_z_zz_1[i] * fe_0 + ta1_z_yz_zz_0[i] * pa_y[i] - ta1_z_yz_zz_1[i] * pc_y[i];
    }

    // Set up 168-174 components of targeted buffer : FD

    auto ta1_z_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 168);

    auto ta1_z_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 169);

    auto ta1_z_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 170);

    auto ta1_z_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 171);

    auto ta1_z_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 172);

    auto ta1_z_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 173);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_yzz_xx_0, \
                             ta1_z_yzz_xy_0, \
                             ta1_z_yzz_xz_0, \
                             ta1_z_yzz_yy_0, \
                             ta1_z_yzz_yz_0, \
                             ta1_z_yzz_zz_0, \
                             ta1_z_zz_x_0,   \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_xx_0,  \
                             ta1_z_zz_xx_1,  \
                             ta1_z_zz_xy_0,  \
                             ta1_z_zz_xy_1,  \
                             ta1_z_zz_xz_0,  \
                             ta1_z_zz_xz_1,  \
                             ta1_z_zz_y_0,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_yy_0,  \
                             ta1_z_zz_yy_1,  \
                             ta1_z_zz_yz_0,  \
                             ta1_z_zz_yz_1,  \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta1_z_zz_zz_0,  \
                             ta1_z_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_xx_0[i] = ta1_z_zz_xx_0[i] * pa_y[i] - ta1_z_zz_xx_1[i] * pc_y[i];

        ta1_z_yzz_xy_0[i] = ta1_z_zz_x_0[i] * fe_0 - ta1_z_zz_x_1[i] * fe_0 + ta1_z_zz_xy_0[i] * pa_y[i] - ta1_z_zz_xy_1[i] * pc_y[i];

        ta1_z_yzz_xz_0[i] = ta1_z_zz_xz_0[i] * pa_y[i] - ta1_z_zz_xz_1[i] * pc_y[i];

        ta1_z_yzz_yy_0[i] = 2.0 * ta1_z_zz_y_0[i] * fe_0 - 2.0 * ta1_z_zz_y_1[i] * fe_0 + ta1_z_zz_yy_0[i] * pa_y[i] - ta1_z_zz_yy_1[i] * pc_y[i];

        ta1_z_yzz_yz_0[i] = ta1_z_zz_z_0[i] * fe_0 - ta1_z_zz_z_1[i] * fe_0 + ta1_z_zz_yz_0[i] * pa_y[i] - ta1_z_zz_yz_1[i] * pc_y[i];

        ta1_z_yzz_zz_0[i] = ta1_z_zz_zz_0[i] * pa_y[i] - ta1_z_zz_zz_1[i] * pc_y[i];
    }

    // Set up 174-180 components of targeted buffer : FD

    auto ta1_z_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 174);

    auto ta1_z_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 175);

    auto ta1_z_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 176);

    auto ta1_z_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 177);

    auto ta1_z_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 178);

    auto ta1_z_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 179);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_z_xx_0,   \
                             ta1_z_z_xx_1,   \
                             ta1_z_z_xy_0,   \
                             ta1_z_z_xy_1,   \
                             ta1_z_z_xz_0,   \
                             ta1_z_z_xz_1,   \
                             ta1_z_z_yy_0,   \
                             ta1_z_z_yy_1,   \
                             ta1_z_z_yz_0,   \
                             ta1_z_z_yz_1,   \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta1_z_zz_x_0,   \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_xx_0,  \
                             ta1_z_zz_xx_1,  \
                             ta1_z_zz_xy_0,  \
                             ta1_z_zz_xy_1,  \
                             ta1_z_zz_xz_0,  \
                             ta1_z_zz_xz_1,  \
                             ta1_z_zz_y_0,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_yy_0,  \
                             ta1_z_zz_yy_1,  \
                             ta1_z_zz_yz_0,  \
                             ta1_z_zz_yz_1,  \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta1_z_zz_zz_0,  \
                             ta1_z_zz_zz_1,  \
                             ta1_z_zzz_xx_0, \
                             ta1_z_zzz_xy_0, \
                             ta1_z_zzz_xz_0, \
                             ta1_z_zzz_yy_0, \
                             ta1_z_zzz_yz_0, \
                             ta1_z_zzz_zz_0, \
                             ta_zz_xx_1,     \
                             ta_zz_xy_1,     \
                             ta_zz_xz_1,     \
                             ta_zz_yy_1,     \
                             ta_zz_yz_1,     \
                             ta_zz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_xx_0[i] =
            2.0 * ta1_z_z_xx_0[i] * fe_0 - 2.0 * ta1_z_z_xx_1[i] * fe_0 + ta_zz_xx_1[i] + ta1_z_zz_xx_0[i] * pa_z[i] - ta1_z_zz_xx_1[i] * pc_z[i];

        ta1_z_zzz_xy_0[i] =
            2.0 * ta1_z_z_xy_0[i] * fe_0 - 2.0 * ta1_z_z_xy_1[i] * fe_0 + ta_zz_xy_1[i] + ta1_z_zz_xy_0[i] * pa_z[i] - ta1_z_zz_xy_1[i] * pc_z[i];

        ta1_z_zzz_xz_0[i] = 2.0 * ta1_z_z_xz_0[i] * fe_0 - 2.0 * ta1_z_z_xz_1[i] * fe_0 + ta1_z_zz_x_0[i] * fe_0 - ta1_z_zz_x_1[i] * fe_0 +
                            ta_zz_xz_1[i] + ta1_z_zz_xz_0[i] * pa_z[i] - ta1_z_zz_xz_1[i] * pc_z[i];

        ta1_z_zzz_yy_0[i] =
            2.0 * ta1_z_z_yy_0[i] * fe_0 - 2.0 * ta1_z_z_yy_1[i] * fe_0 + ta_zz_yy_1[i] + ta1_z_zz_yy_0[i] * pa_z[i] - ta1_z_zz_yy_1[i] * pc_z[i];

        ta1_z_zzz_yz_0[i] = 2.0 * ta1_z_z_yz_0[i] * fe_0 - 2.0 * ta1_z_z_yz_1[i] * fe_0 + ta1_z_zz_y_0[i] * fe_0 - ta1_z_zz_y_1[i] * fe_0 +
                            ta_zz_yz_1[i] + ta1_z_zz_yz_0[i] * pa_z[i] - ta1_z_zz_yz_1[i] * pc_z[i];

        ta1_z_zzz_zz_0[i] = 2.0 * ta1_z_z_zz_0[i] * fe_0 - 2.0 * ta1_z_z_zz_1[i] * fe_0 + 2.0 * ta1_z_zz_z_0[i] * fe_0 -
                            2.0 * ta1_z_zz_z_1[i] * fe_0 + ta_zz_zz_1[i] + ta1_z_zz_zz_0[i] * pa_z[i] - ta1_z_zz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
