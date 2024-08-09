#include "NuclearPotentialGeom010PrimRecDD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_dd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_dd,
                                        const size_t              idx_npot_geom_010_0_sd,
                                        const size_t              idx_npot_geom_010_1_sd,
                                        const size_t              idx_npot_geom_010_0_pp,
                                        const size_t              idx_npot_geom_010_1_pp,
                                        const size_t              idx_npot_1_pd,
                                        const size_t              idx_npot_geom_010_0_pd,
                                        const size_t              idx_npot_geom_010_1_pd,
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

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd);

    auto ta1_x_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 1);

    auto ta1_x_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 2);

    auto ta1_x_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 3);

    auto ta1_x_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 4);

    auto ta1_x_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 5);

    auto ta1_y_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 6);

    auto ta1_y_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 7);

    auto ta1_y_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 8);

    auto ta1_y_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 9);

    auto ta1_y_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 10);

    auto ta1_y_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 11);

    auto ta1_z_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 12);

    auto ta1_z_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 13);

    auto ta1_z_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 14);

    auto ta1_z_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 15);

    auto ta1_z_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 16);

    auto ta1_z_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 1);

    auto ta1_x_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 2);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 4);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 7);

    auto ta1_y_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 8);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 10);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 13);

    auto ta1_z_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 14);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 16);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

    // Set up components of auxiliary buffer : PP

    auto ta1_x_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp);

    auto ta1_x_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 1);

    auto ta1_x_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 2);

    auto ta1_x_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 3);

    auto ta1_x_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 4);

    auto ta1_x_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 5);

    auto ta1_x_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 6);

    auto ta1_x_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 7);

    auto ta1_x_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 8);

    auto ta1_y_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 9);

    auto ta1_y_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 10);

    auto ta1_y_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 11);

    auto ta1_y_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 12);

    auto ta1_y_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 13);

    auto ta1_y_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 14);

    auto ta1_y_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 15);

    auto ta1_y_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 16);

    auto ta1_y_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 17);

    auto ta1_z_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 18);

    auto ta1_z_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 19);

    auto ta1_z_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 20);

    auto ta1_z_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 21);

    auto ta1_z_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 22);

    auto ta1_z_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 23);

    auto ta1_z_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 24);

    auto ta1_z_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 25);

    auto ta1_z_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 26);

    // Set up components of auxiliary buffer : PP

    auto ta1_x_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp);

    auto ta1_x_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 1);

    auto ta1_x_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 2);

    auto ta1_x_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 3);

    auto ta1_x_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 4);

    auto ta1_x_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 5);

    auto ta1_x_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 6);

    auto ta1_x_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 7);

    auto ta1_x_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 8);

    auto ta1_y_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 9);

    auto ta1_y_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 10);

    auto ta1_y_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 11);

    auto ta1_y_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 12);

    auto ta1_y_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 13);

    auto ta1_y_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 14);

    auto ta1_y_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 15);

    auto ta1_y_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 16);

    auto ta1_y_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 17);

    auto ta1_z_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 18);

    auto ta1_z_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 19);

    auto ta1_z_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 20);

    auto ta1_z_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 21);

    auto ta1_z_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 22);

    auto ta1_z_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 23);

    auto ta1_z_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 24);

    auto ta1_z_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 25);

    auto ta1_z_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 26);

    // Set up components of auxiliary buffer : PD

    auto ta_x_xx_1 = pbuffer.data(idx_npot_1_pd);

    auto ta_x_xy_1 = pbuffer.data(idx_npot_1_pd + 1);

    auto ta_x_xz_1 = pbuffer.data(idx_npot_1_pd + 2);

    auto ta_x_yy_1 = pbuffer.data(idx_npot_1_pd + 3);

    auto ta_x_yz_1 = pbuffer.data(idx_npot_1_pd + 4);

    auto ta_x_zz_1 = pbuffer.data(idx_npot_1_pd + 5);

    auto ta_y_xx_1 = pbuffer.data(idx_npot_1_pd + 6);

    auto ta_y_xy_1 = pbuffer.data(idx_npot_1_pd + 7);

    auto ta_y_xz_1 = pbuffer.data(idx_npot_1_pd + 8);

    auto ta_y_yy_1 = pbuffer.data(idx_npot_1_pd + 9);

    auto ta_y_yz_1 = pbuffer.data(idx_npot_1_pd + 10);

    auto ta_y_zz_1 = pbuffer.data(idx_npot_1_pd + 11);

    auto ta_z_xx_1 = pbuffer.data(idx_npot_1_pd + 12);

    auto ta_z_xy_1 = pbuffer.data(idx_npot_1_pd + 13);

    auto ta_z_xz_1 = pbuffer.data(idx_npot_1_pd + 14);

    auto ta_z_yy_1 = pbuffer.data(idx_npot_1_pd + 15);

    auto ta_z_yz_1 = pbuffer.data(idx_npot_1_pd + 16);

    auto ta_z_zz_1 = pbuffer.data(idx_npot_1_pd + 17);

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

    // Set up 0-6 components of targeted buffer : DD

    auto ta1_x_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd);

    auto ta1_x_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 1);

    auto ta1_x_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 2);

    auto ta1_x_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 3);

    auto ta1_x_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 4);

    auto ta1_x_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 5);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_x_xx_0,  \
                             ta1_x_x_xx_1,  \
                             ta1_x_x_xy_0,  \
                             ta1_x_x_xy_1,  \
                             ta1_x_x_xz_0,  \
                             ta1_x_x_xz_1,  \
                             ta1_x_x_y_0,   \
                             ta1_x_x_y_1,   \
                             ta1_x_x_yy_0,  \
                             ta1_x_x_yy_1,  \
                             ta1_x_x_yz_0,  \
                             ta1_x_x_yz_1,  \
                             ta1_x_x_z_0,   \
                             ta1_x_x_z_1,   \
                             ta1_x_x_zz_0,  \
                             ta1_x_x_zz_1,  \
                             ta1_x_xx_xx_0, \
                             ta1_x_xx_xy_0, \
                             ta1_x_xx_xz_0, \
                             ta1_x_xx_yy_0, \
                             ta1_x_xx_yz_0, \
                             ta1_x_xx_zz_0, \
                             ta_x_xx_1,     \
                             ta_x_xy_1,     \
                             ta_x_xz_1,     \
                             ta_x_yy_1,     \
                             ta_x_yz_1,     \
                             ta_x_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_xx_0[i] = ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + 2.0 * ta1_x_x_x_0[i] * fe_0 -
                           2.0 * ta1_x_x_x_1[i] * fe_0 + ta_x_xx_1[i] + ta1_x_x_xx_0[i] * pa_x[i] -
                           ta1_x_x_xx_1[i] * pc_x[i];

        ta1_x_xx_xy_0[i] = ta1_x_0_xy_0[i] * fe_0 - ta1_x_0_xy_1[i] * fe_0 + ta1_x_x_y_0[i] * fe_0 -
                           ta1_x_x_y_1[i] * fe_0 + ta_x_xy_1[i] + ta1_x_x_xy_0[i] * pa_x[i] - ta1_x_x_xy_1[i] * pc_x[i];

        ta1_x_xx_xz_0[i] = ta1_x_0_xz_0[i] * fe_0 - ta1_x_0_xz_1[i] * fe_0 + ta1_x_x_z_0[i] * fe_0 -
                           ta1_x_x_z_1[i] * fe_0 + ta_x_xz_1[i] + ta1_x_x_xz_0[i] * pa_x[i] - ta1_x_x_xz_1[i] * pc_x[i];

        ta1_x_xx_yy_0[i] = ta1_x_0_yy_0[i] * fe_0 - ta1_x_0_yy_1[i] * fe_0 + ta_x_yy_1[i] + ta1_x_x_yy_0[i] * pa_x[i] -
                           ta1_x_x_yy_1[i] * pc_x[i];

        ta1_x_xx_yz_0[i] = ta1_x_0_yz_0[i] * fe_0 - ta1_x_0_yz_1[i] * fe_0 + ta_x_yz_1[i] + ta1_x_x_yz_0[i] * pa_x[i] -
                           ta1_x_x_yz_1[i] * pc_x[i];

        ta1_x_xx_zz_0[i] = ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + ta_x_zz_1[i] + ta1_x_x_zz_0[i] * pa_x[i] -
                           ta1_x_x_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto ta1_x_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 6);

    auto ta1_x_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 7);

    auto ta1_x_xy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 8);

    auto ta1_x_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 9);

    auto ta1_x_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 10);

    auto ta1_x_xy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 11);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_x_xx_0,  \
                             ta1_x_x_xx_1,  \
                             ta1_x_x_xy_0,  \
                             ta1_x_x_xy_1,  \
                             ta1_x_x_xz_0,  \
                             ta1_x_x_xz_1,  \
                             ta1_x_x_zz_0,  \
                             ta1_x_x_zz_1,  \
                             ta1_x_xy_xx_0, \
                             ta1_x_xy_xy_0, \
                             ta1_x_xy_xz_0, \
                             ta1_x_xy_yy_0, \
                             ta1_x_xy_yz_0, \
                             ta1_x_xy_zz_0, \
                             ta1_x_y_yy_0,  \
                             ta1_x_y_yy_1,  \
                             ta1_x_y_yz_0,  \
                             ta1_x_y_yz_1,  \
                             ta_y_yy_1,     \
                             ta_y_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xy_xx_0[i] = ta1_x_x_xx_0[i] * pa_y[i] - ta1_x_x_xx_1[i] * pc_y[i];

        ta1_x_xy_xy_0[i] =
            ta1_x_x_x_0[i] * fe_0 - ta1_x_x_x_1[i] * fe_0 + ta1_x_x_xy_0[i] * pa_y[i] - ta1_x_x_xy_1[i] * pc_y[i];

        ta1_x_xy_xz_0[i] = ta1_x_x_xz_0[i] * pa_y[i] - ta1_x_x_xz_1[i] * pc_y[i];

        ta1_x_xy_yy_0[i] = ta_y_yy_1[i] + ta1_x_y_yy_0[i] * pa_x[i] - ta1_x_y_yy_1[i] * pc_x[i];

        ta1_x_xy_yz_0[i] = ta_y_yz_1[i] + ta1_x_y_yz_0[i] * pa_x[i] - ta1_x_y_yz_1[i] * pc_x[i];

        ta1_x_xy_zz_0[i] = ta1_x_x_zz_0[i] * pa_y[i] - ta1_x_x_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto ta1_x_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 12);

    auto ta1_x_xz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 13);

    auto ta1_x_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 14);

    auto ta1_x_xz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 15);

    auto ta1_x_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 16);

    auto ta1_x_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 17);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_x_xx_0,  \
                             ta1_x_x_xx_1,  \
                             ta1_x_x_xy_0,  \
                             ta1_x_x_xy_1,  \
                             ta1_x_x_xz_0,  \
                             ta1_x_x_xz_1,  \
                             ta1_x_x_yy_0,  \
                             ta1_x_x_yy_1,  \
                             ta1_x_xz_xx_0, \
                             ta1_x_xz_xy_0, \
                             ta1_x_xz_xz_0, \
                             ta1_x_xz_yy_0, \
                             ta1_x_xz_yz_0, \
                             ta1_x_xz_zz_0, \
                             ta1_x_z_yz_0,  \
                             ta1_x_z_yz_1,  \
                             ta1_x_z_zz_0,  \
                             ta1_x_z_zz_1,  \
                             ta_z_yz_1,     \
                             ta_z_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xz_xx_0[i] = ta1_x_x_xx_0[i] * pa_z[i] - ta1_x_x_xx_1[i] * pc_z[i];

        ta1_x_xz_xy_0[i] = ta1_x_x_xy_0[i] * pa_z[i] - ta1_x_x_xy_1[i] * pc_z[i];

        ta1_x_xz_xz_0[i] =
            ta1_x_x_x_0[i] * fe_0 - ta1_x_x_x_1[i] * fe_0 + ta1_x_x_xz_0[i] * pa_z[i] - ta1_x_x_xz_1[i] * pc_z[i];

        ta1_x_xz_yy_0[i] = ta1_x_x_yy_0[i] * pa_z[i] - ta1_x_x_yy_1[i] * pc_z[i];

        ta1_x_xz_yz_0[i] = ta_z_yz_1[i] + ta1_x_z_yz_0[i] * pa_x[i] - ta1_x_z_yz_1[i] * pc_x[i];

        ta1_x_xz_zz_0[i] = ta_z_zz_1[i] + ta1_x_z_zz_0[i] * pa_x[i] - ta1_x_z_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto ta1_x_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 18);

    auto ta1_x_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 19);

    auto ta1_x_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 20);

    auto ta1_x_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 21);

    auto ta1_x_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 22);

    auto ta1_x_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 23);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_y_x_0,   \
                             ta1_x_y_x_1,   \
                             ta1_x_y_xx_0,  \
                             ta1_x_y_xx_1,  \
                             ta1_x_y_xy_0,  \
                             ta1_x_y_xy_1,  \
                             ta1_x_y_xz_0,  \
                             ta1_x_y_xz_1,  \
                             ta1_x_y_y_0,   \
                             ta1_x_y_y_1,   \
                             ta1_x_y_yy_0,  \
                             ta1_x_y_yy_1,  \
                             ta1_x_y_yz_0,  \
                             ta1_x_y_yz_1,  \
                             ta1_x_y_z_0,   \
                             ta1_x_y_z_1,   \
                             ta1_x_y_zz_0,  \
                             ta1_x_y_zz_1,  \
                             ta1_x_yy_xx_0, \
                             ta1_x_yy_xy_0, \
                             ta1_x_yy_xz_0, \
                             ta1_x_yy_yy_0, \
                             ta1_x_yy_yz_0, \
                             ta1_x_yy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_xx_0[i] =
            ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_y_xx_0[i] * pa_y[i] - ta1_x_y_xx_1[i] * pc_y[i];

        ta1_x_yy_xy_0[i] = ta1_x_0_xy_0[i] * fe_0 - ta1_x_0_xy_1[i] * fe_0 + ta1_x_y_x_0[i] * fe_0 -
                           ta1_x_y_x_1[i] * fe_0 + ta1_x_y_xy_0[i] * pa_y[i] - ta1_x_y_xy_1[i] * pc_y[i];

        ta1_x_yy_xz_0[i] =
            ta1_x_0_xz_0[i] * fe_0 - ta1_x_0_xz_1[i] * fe_0 + ta1_x_y_xz_0[i] * pa_y[i] - ta1_x_y_xz_1[i] * pc_y[i];

        ta1_x_yy_yy_0[i] = ta1_x_0_yy_0[i] * fe_0 - ta1_x_0_yy_1[i] * fe_0 + 2.0 * ta1_x_y_y_0[i] * fe_0 -
                           2.0 * ta1_x_y_y_1[i] * fe_0 + ta1_x_y_yy_0[i] * pa_y[i] - ta1_x_y_yy_1[i] * pc_y[i];

        ta1_x_yy_yz_0[i] = ta1_x_0_yz_0[i] * fe_0 - ta1_x_0_yz_1[i] * fe_0 + ta1_x_y_z_0[i] * fe_0 -
                           ta1_x_y_z_1[i] * fe_0 + ta1_x_y_yz_0[i] * pa_y[i] - ta1_x_y_yz_1[i] * pc_y[i];

        ta1_x_yy_zz_0[i] =
            ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + ta1_x_y_zz_0[i] * pa_y[i] - ta1_x_y_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto ta1_x_yz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 24);

    auto ta1_x_yz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 25);

    auto ta1_x_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 26);

    auto ta1_x_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 27);

    auto ta1_x_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 28);

    auto ta1_x_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 29);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_x_y_xy_0,  \
                             ta1_x_y_xy_1,  \
                             ta1_x_y_yy_0,  \
                             ta1_x_y_yy_1,  \
                             ta1_x_yz_xx_0, \
                             ta1_x_yz_xy_0, \
                             ta1_x_yz_xz_0, \
                             ta1_x_yz_yy_0, \
                             ta1_x_yz_yz_0, \
                             ta1_x_yz_zz_0, \
                             ta1_x_z_xx_0,  \
                             ta1_x_z_xx_1,  \
                             ta1_x_z_xz_0,  \
                             ta1_x_z_xz_1,  \
                             ta1_x_z_yz_0,  \
                             ta1_x_z_yz_1,  \
                             ta1_x_z_z_0,   \
                             ta1_x_z_z_1,   \
                             ta1_x_z_zz_0,  \
                             ta1_x_z_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yz_xx_0[i] = ta1_x_z_xx_0[i] * pa_y[i] - ta1_x_z_xx_1[i] * pc_y[i];

        ta1_x_yz_xy_0[i] = ta1_x_y_xy_0[i] * pa_z[i] - ta1_x_y_xy_1[i] * pc_z[i];

        ta1_x_yz_xz_0[i] = ta1_x_z_xz_0[i] * pa_y[i] - ta1_x_z_xz_1[i] * pc_y[i];

        ta1_x_yz_yy_0[i] = ta1_x_y_yy_0[i] * pa_z[i] - ta1_x_y_yy_1[i] * pc_z[i];

        ta1_x_yz_yz_0[i] =
            ta1_x_z_z_0[i] * fe_0 - ta1_x_z_z_1[i] * fe_0 + ta1_x_z_yz_0[i] * pa_y[i] - ta1_x_z_yz_1[i] * pc_y[i];

        ta1_x_yz_zz_0[i] = ta1_x_z_zz_0[i] * pa_y[i] - ta1_x_z_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto ta1_x_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 30);

    auto ta1_x_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 31);

    auto ta1_x_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 32);

    auto ta1_x_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 33);

    auto ta1_x_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 34);

    auto ta1_x_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 35);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_z_x_0,   \
                             ta1_x_z_x_1,   \
                             ta1_x_z_xx_0,  \
                             ta1_x_z_xx_1,  \
                             ta1_x_z_xy_0,  \
                             ta1_x_z_xy_1,  \
                             ta1_x_z_xz_0,  \
                             ta1_x_z_xz_1,  \
                             ta1_x_z_y_0,   \
                             ta1_x_z_y_1,   \
                             ta1_x_z_yy_0,  \
                             ta1_x_z_yy_1,  \
                             ta1_x_z_yz_0,  \
                             ta1_x_z_yz_1,  \
                             ta1_x_z_z_0,   \
                             ta1_x_z_z_1,   \
                             ta1_x_z_zz_0,  \
                             ta1_x_z_zz_1,  \
                             ta1_x_zz_xx_0, \
                             ta1_x_zz_xy_0, \
                             ta1_x_zz_xz_0, \
                             ta1_x_zz_yy_0, \
                             ta1_x_zz_yz_0, \
                             ta1_x_zz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_xx_0[i] =
            ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_z_xx_0[i] * pa_z[i] - ta1_x_z_xx_1[i] * pc_z[i];

        ta1_x_zz_xy_0[i] =
            ta1_x_0_xy_0[i] * fe_0 - ta1_x_0_xy_1[i] * fe_0 + ta1_x_z_xy_0[i] * pa_z[i] - ta1_x_z_xy_1[i] * pc_z[i];

        ta1_x_zz_xz_0[i] = ta1_x_0_xz_0[i] * fe_0 - ta1_x_0_xz_1[i] * fe_0 + ta1_x_z_x_0[i] * fe_0 -
                           ta1_x_z_x_1[i] * fe_0 + ta1_x_z_xz_0[i] * pa_z[i] - ta1_x_z_xz_1[i] * pc_z[i];

        ta1_x_zz_yy_0[i] =
            ta1_x_0_yy_0[i] * fe_0 - ta1_x_0_yy_1[i] * fe_0 + ta1_x_z_yy_0[i] * pa_z[i] - ta1_x_z_yy_1[i] * pc_z[i];

        ta1_x_zz_yz_0[i] = ta1_x_0_yz_0[i] * fe_0 - ta1_x_0_yz_1[i] * fe_0 + ta1_x_z_y_0[i] * fe_0 -
                           ta1_x_z_y_1[i] * fe_0 + ta1_x_z_yz_0[i] * pa_z[i] - ta1_x_z_yz_1[i] * pc_z[i];

        ta1_x_zz_zz_0[i] = ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + 2.0 * ta1_x_z_z_0[i] * fe_0 -
                           2.0 * ta1_x_z_z_1[i] * fe_0 + ta1_x_z_zz_0[i] * pa_z[i] - ta1_x_z_zz_1[i] * pc_z[i];
    }

    // Set up 36-42 components of targeted buffer : DD

    auto ta1_y_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 36);

    auto ta1_y_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 37);

    auto ta1_y_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 38);

    auto ta1_y_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 39);

    auto ta1_y_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 40);

    auto ta1_y_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 41);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_x_x_0,   \
                             ta1_y_x_x_1,   \
                             ta1_y_x_xx_0,  \
                             ta1_y_x_xx_1,  \
                             ta1_y_x_xy_0,  \
                             ta1_y_x_xy_1,  \
                             ta1_y_x_xz_0,  \
                             ta1_y_x_xz_1,  \
                             ta1_y_x_y_0,   \
                             ta1_y_x_y_1,   \
                             ta1_y_x_yy_0,  \
                             ta1_y_x_yy_1,  \
                             ta1_y_x_yz_0,  \
                             ta1_y_x_yz_1,  \
                             ta1_y_x_z_0,   \
                             ta1_y_x_z_1,   \
                             ta1_y_x_zz_0,  \
                             ta1_y_x_zz_1,  \
                             ta1_y_xx_xx_0, \
                             ta1_y_xx_xy_0, \
                             ta1_y_xx_xz_0, \
                             ta1_y_xx_yy_0, \
                             ta1_y_xx_yz_0, \
                             ta1_y_xx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_xx_0[i] = ta1_y_0_xx_0[i] * fe_0 - ta1_y_0_xx_1[i] * fe_0 + 2.0 * ta1_y_x_x_0[i] * fe_0 -
                           2.0 * ta1_y_x_x_1[i] * fe_0 + ta1_y_x_xx_0[i] * pa_x[i] - ta1_y_x_xx_1[i] * pc_x[i];

        ta1_y_xx_xy_0[i] = ta1_y_0_xy_0[i] * fe_0 - ta1_y_0_xy_1[i] * fe_0 + ta1_y_x_y_0[i] * fe_0 -
                           ta1_y_x_y_1[i] * fe_0 + ta1_y_x_xy_0[i] * pa_x[i] - ta1_y_x_xy_1[i] * pc_x[i];

        ta1_y_xx_xz_0[i] = ta1_y_0_xz_0[i] * fe_0 - ta1_y_0_xz_1[i] * fe_0 + ta1_y_x_z_0[i] * fe_0 -
                           ta1_y_x_z_1[i] * fe_0 + ta1_y_x_xz_0[i] * pa_x[i] - ta1_y_x_xz_1[i] * pc_x[i];

        ta1_y_xx_yy_0[i] =
            ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_x_yy_0[i] * pa_x[i] - ta1_y_x_yy_1[i] * pc_x[i];

        ta1_y_xx_yz_0[i] =
            ta1_y_0_yz_0[i] * fe_0 - ta1_y_0_yz_1[i] * fe_0 + ta1_y_x_yz_0[i] * pa_x[i] - ta1_y_x_yz_1[i] * pc_x[i];

        ta1_y_xx_zz_0[i] =
            ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + ta1_y_x_zz_0[i] * pa_x[i] - ta1_y_x_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : DD

    auto ta1_y_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 42);

    auto ta1_y_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 43);

    auto ta1_y_xy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 44);

    auto ta1_y_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 45);

    auto ta1_y_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 46);

    auto ta1_y_xy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 47);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_y_x_xx_0,  \
                             ta1_y_x_xx_1,  \
                             ta1_y_x_xz_0,  \
                             ta1_y_x_xz_1,  \
                             ta1_y_xy_xx_0, \
                             ta1_y_xy_xy_0, \
                             ta1_y_xy_xz_0, \
                             ta1_y_xy_yy_0, \
                             ta1_y_xy_yz_0, \
                             ta1_y_xy_zz_0, \
                             ta1_y_y_xy_0,  \
                             ta1_y_y_xy_1,  \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta1_y_y_yy_0,  \
                             ta1_y_y_yy_1,  \
                             ta1_y_y_yz_0,  \
                             ta1_y_y_yz_1,  \
                             ta1_y_y_zz_0,  \
                             ta1_y_y_zz_1,  \
                             ta_x_xx_1,     \
                             ta_x_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xy_xx_0[i] = ta_x_xx_1[i] + ta1_y_x_xx_0[i] * pa_y[i] - ta1_y_x_xx_1[i] * pc_y[i];

        ta1_y_xy_xy_0[i] =
            ta1_y_y_y_0[i] * fe_0 - ta1_y_y_y_1[i] * fe_0 + ta1_y_y_xy_0[i] * pa_x[i] - ta1_y_y_xy_1[i] * pc_x[i];

        ta1_y_xy_xz_0[i] = ta_x_xz_1[i] + ta1_y_x_xz_0[i] * pa_y[i] - ta1_y_x_xz_1[i] * pc_y[i];

        ta1_y_xy_yy_0[i] = ta1_y_y_yy_0[i] * pa_x[i] - ta1_y_y_yy_1[i] * pc_x[i];

        ta1_y_xy_yz_0[i] = ta1_y_y_yz_0[i] * pa_x[i] - ta1_y_y_yz_1[i] * pc_x[i];

        ta1_y_xy_zz_0[i] = ta1_y_y_zz_0[i] * pa_x[i] - ta1_y_y_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : DD

    auto ta1_y_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 48);

    auto ta1_y_xz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 49);

    auto ta1_y_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 50);

    auto ta1_y_xz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 51);

    auto ta1_y_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 52);

    auto ta1_y_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 53);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_y_x_xx_0,  \
                             ta1_y_x_xx_1,  \
                             ta1_y_x_xy_0,  \
                             ta1_y_x_xy_1,  \
                             ta1_y_xz_xx_0, \
                             ta1_y_xz_xy_0, \
                             ta1_y_xz_xz_0, \
                             ta1_y_xz_yy_0, \
                             ta1_y_xz_yz_0, \
                             ta1_y_xz_zz_0, \
                             ta1_y_z_xz_0,  \
                             ta1_y_z_xz_1,  \
                             ta1_y_z_yy_0,  \
                             ta1_y_z_yy_1,  \
                             ta1_y_z_yz_0,  \
                             ta1_y_z_yz_1,  \
                             ta1_y_z_z_0,   \
                             ta1_y_z_z_1,   \
                             ta1_y_z_zz_0,  \
                             ta1_y_z_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xz_xx_0[i] = ta1_y_x_xx_0[i] * pa_z[i] - ta1_y_x_xx_1[i] * pc_z[i];

        ta1_y_xz_xy_0[i] = ta1_y_x_xy_0[i] * pa_z[i] - ta1_y_x_xy_1[i] * pc_z[i];

        ta1_y_xz_xz_0[i] =
            ta1_y_z_z_0[i] * fe_0 - ta1_y_z_z_1[i] * fe_0 + ta1_y_z_xz_0[i] * pa_x[i] - ta1_y_z_xz_1[i] * pc_x[i];

        ta1_y_xz_yy_0[i] = ta1_y_z_yy_0[i] * pa_x[i] - ta1_y_z_yy_1[i] * pc_x[i];

        ta1_y_xz_yz_0[i] = ta1_y_z_yz_0[i] * pa_x[i] - ta1_y_z_yz_1[i] * pc_x[i];

        ta1_y_xz_zz_0[i] = ta1_y_z_zz_0[i] * pa_x[i] - ta1_y_z_zz_1[i] * pc_x[i];
    }

    // Set up 54-60 components of targeted buffer : DD

    auto ta1_y_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 54);

    auto ta1_y_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 55);

    auto ta1_y_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 56);

    auto ta1_y_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 57);

    auto ta1_y_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 58);

    auto ta1_y_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 59);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_y_x_0,   \
                             ta1_y_y_x_1,   \
                             ta1_y_y_xx_0,  \
                             ta1_y_y_xx_1,  \
                             ta1_y_y_xy_0,  \
                             ta1_y_y_xy_1,  \
                             ta1_y_y_xz_0,  \
                             ta1_y_y_xz_1,  \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta1_y_y_yy_0,  \
                             ta1_y_y_yy_1,  \
                             ta1_y_y_yz_0,  \
                             ta1_y_y_yz_1,  \
                             ta1_y_y_z_0,   \
                             ta1_y_y_z_1,   \
                             ta1_y_y_zz_0,  \
                             ta1_y_y_zz_1,  \
                             ta1_y_yy_xx_0, \
                             ta1_y_yy_xy_0, \
                             ta1_y_yy_xz_0, \
                             ta1_y_yy_yy_0, \
                             ta1_y_yy_yz_0, \
                             ta1_y_yy_zz_0, \
                             ta_y_xx_1,     \
                             ta_y_xy_1,     \
                             ta_y_xz_1,     \
                             ta_y_yy_1,     \
                             ta_y_yz_1,     \
                             ta_y_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_xx_0[i] = ta1_y_0_xx_0[i] * fe_0 - ta1_y_0_xx_1[i] * fe_0 + ta_y_xx_1[i] + ta1_y_y_xx_0[i] * pa_y[i] -
                           ta1_y_y_xx_1[i] * pc_y[i];

        ta1_y_yy_xy_0[i] = ta1_y_0_xy_0[i] * fe_0 - ta1_y_0_xy_1[i] * fe_0 + ta1_y_y_x_0[i] * fe_0 -
                           ta1_y_y_x_1[i] * fe_0 + ta_y_xy_1[i] + ta1_y_y_xy_0[i] * pa_y[i] - ta1_y_y_xy_1[i] * pc_y[i];

        ta1_y_yy_xz_0[i] = ta1_y_0_xz_0[i] * fe_0 - ta1_y_0_xz_1[i] * fe_0 + ta_y_xz_1[i] + ta1_y_y_xz_0[i] * pa_y[i] -
                           ta1_y_y_xz_1[i] * pc_y[i];

        ta1_y_yy_yy_0[i] = ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + 2.0 * ta1_y_y_y_0[i] * fe_0 -
                           2.0 * ta1_y_y_y_1[i] * fe_0 + ta_y_yy_1[i] + ta1_y_y_yy_0[i] * pa_y[i] -
                           ta1_y_y_yy_1[i] * pc_y[i];

        ta1_y_yy_yz_0[i] = ta1_y_0_yz_0[i] * fe_0 - ta1_y_0_yz_1[i] * fe_0 + ta1_y_y_z_0[i] * fe_0 -
                           ta1_y_y_z_1[i] * fe_0 + ta_y_yz_1[i] + ta1_y_y_yz_0[i] * pa_y[i] - ta1_y_y_yz_1[i] * pc_y[i];

        ta1_y_yy_zz_0[i] = ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + ta_y_zz_1[i] + ta1_y_y_zz_0[i] * pa_y[i] -
                           ta1_y_y_zz_1[i] * pc_y[i];
    }

    // Set up 60-66 components of targeted buffer : DD

    auto ta1_y_yz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 60);

    auto ta1_y_yz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 61);

    auto ta1_y_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 62);

    auto ta1_y_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 63);

    auto ta1_y_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 64);

    auto ta1_y_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 65);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_y_y_xx_0,  \
                             ta1_y_y_xx_1,  \
                             ta1_y_y_xy_0,  \
                             ta1_y_y_xy_1,  \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta1_y_y_yy_0,  \
                             ta1_y_y_yy_1,  \
                             ta1_y_y_yz_0,  \
                             ta1_y_y_yz_1,  \
                             ta1_y_yz_xx_0, \
                             ta1_y_yz_xy_0, \
                             ta1_y_yz_xz_0, \
                             ta1_y_yz_yy_0, \
                             ta1_y_yz_yz_0, \
                             ta1_y_yz_zz_0, \
                             ta1_y_z_xz_0,  \
                             ta1_y_z_xz_1,  \
                             ta1_y_z_zz_0,  \
                             ta1_y_z_zz_1,  \
                             ta_z_xz_1,     \
                             ta_z_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yz_xx_0[i] = ta1_y_y_xx_0[i] * pa_z[i] - ta1_y_y_xx_1[i] * pc_z[i];

        ta1_y_yz_xy_0[i] = ta1_y_y_xy_0[i] * pa_z[i] - ta1_y_y_xy_1[i] * pc_z[i];

        ta1_y_yz_xz_0[i] = ta_z_xz_1[i] + ta1_y_z_xz_0[i] * pa_y[i] - ta1_y_z_xz_1[i] * pc_y[i];

        ta1_y_yz_yy_0[i] = ta1_y_y_yy_0[i] * pa_z[i] - ta1_y_y_yy_1[i] * pc_z[i];

        ta1_y_yz_yz_0[i] =
            ta1_y_y_y_0[i] * fe_0 - ta1_y_y_y_1[i] * fe_0 + ta1_y_y_yz_0[i] * pa_z[i] - ta1_y_y_yz_1[i] * pc_z[i];

        ta1_y_yz_zz_0[i] = ta_z_zz_1[i] + ta1_y_z_zz_0[i] * pa_y[i] - ta1_y_z_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : DD

    auto ta1_y_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 66);

    auto ta1_y_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 67);

    auto ta1_y_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 68);

    auto ta1_y_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 69);

    auto ta1_y_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 70);

    auto ta1_y_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 71);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_z_x_0,   \
                             ta1_y_z_x_1,   \
                             ta1_y_z_xx_0,  \
                             ta1_y_z_xx_1,  \
                             ta1_y_z_xy_0,  \
                             ta1_y_z_xy_1,  \
                             ta1_y_z_xz_0,  \
                             ta1_y_z_xz_1,  \
                             ta1_y_z_y_0,   \
                             ta1_y_z_y_1,   \
                             ta1_y_z_yy_0,  \
                             ta1_y_z_yy_1,  \
                             ta1_y_z_yz_0,  \
                             ta1_y_z_yz_1,  \
                             ta1_y_z_z_0,   \
                             ta1_y_z_z_1,   \
                             ta1_y_z_zz_0,  \
                             ta1_y_z_zz_1,  \
                             ta1_y_zz_xx_0, \
                             ta1_y_zz_xy_0, \
                             ta1_y_zz_xz_0, \
                             ta1_y_zz_yy_0, \
                             ta1_y_zz_yz_0, \
                             ta1_y_zz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_xx_0[i] =
            ta1_y_0_xx_0[i] * fe_0 - ta1_y_0_xx_1[i] * fe_0 + ta1_y_z_xx_0[i] * pa_z[i] - ta1_y_z_xx_1[i] * pc_z[i];

        ta1_y_zz_xy_0[i] =
            ta1_y_0_xy_0[i] * fe_0 - ta1_y_0_xy_1[i] * fe_0 + ta1_y_z_xy_0[i] * pa_z[i] - ta1_y_z_xy_1[i] * pc_z[i];

        ta1_y_zz_xz_0[i] = ta1_y_0_xz_0[i] * fe_0 - ta1_y_0_xz_1[i] * fe_0 + ta1_y_z_x_0[i] * fe_0 -
                           ta1_y_z_x_1[i] * fe_0 + ta1_y_z_xz_0[i] * pa_z[i] - ta1_y_z_xz_1[i] * pc_z[i];

        ta1_y_zz_yy_0[i] =
            ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_z_yy_0[i] * pa_z[i] - ta1_y_z_yy_1[i] * pc_z[i];

        ta1_y_zz_yz_0[i] = ta1_y_0_yz_0[i] * fe_0 - ta1_y_0_yz_1[i] * fe_0 + ta1_y_z_y_0[i] * fe_0 -
                           ta1_y_z_y_1[i] * fe_0 + ta1_y_z_yz_0[i] * pa_z[i] - ta1_y_z_yz_1[i] * pc_z[i];

        ta1_y_zz_zz_0[i] = ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + 2.0 * ta1_y_z_z_0[i] * fe_0 -
                           2.0 * ta1_y_z_z_1[i] * fe_0 + ta1_y_z_zz_0[i] * pa_z[i] - ta1_y_z_zz_1[i] * pc_z[i];
    }

    // Set up 72-78 components of targeted buffer : DD

    auto ta1_z_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 72);

    auto ta1_z_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 73);

    auto ta1_z_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 74);

    auto ta1_z_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 75);

    auto ta1_z_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 76);

    auto ta1_z_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 77);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_x_x_0,   \
                             ta1_z_x_x_1,   \
                             ta1_z_x_xx_0,  \
                             ta1_z_x_xx_1,  \
                             ta1_z_x_xy_0,  \
                             ta1_z_x_xy_1,  \
                             ta1_z_x_xz_0,  \
                             ta1_z_x_xz_1,  \
                             ta1_z_x_y_0,   \
                             ta1_z_x_y_1,   \
                             ta1_z_x_yy_0,  \
                             ta1_z_x_yy_1,  \
                             ta1_z_x_yz_0,  \
                             ta1_z_x_yz_1,  \
                             ta1_z_x_z_0,   \
                             ta1_z_x_z_1,   \
                             ta1_z_x_zz_0,  \
                             ta1_z_x_zz_1,  \
                             ta1_z_xx_xx_0, \
                             ta1_z_xx_xy_0, \
                             ta1_z_xx_xz_0, \
                             ta1_z_xx_yy_0, \
                             ta1_z_xx_yz_0, \
                             ta1_z_xx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_xx_0[i] = ta1_z_0_xx_0[i] * fe_0 - ta1_z_0_xx_1[i] * fe_0 + 2.0 * ta1_z_x_x_0[i] * fe_0 -
                           2.0 * ta1_z_x_x_1[i] * fe_0 + ta1_z_x_xx_0[i] * pa_x[i] - ta1_z_x_xx_1[i] * pc_x[i];

        ta1_z_xx_xy_0[i] = ta1_z_0_xy_0[i] * fe_0 - ta1_z_0_xy_1[i] * fe_0 + ta1_z_x_y_0[i] * fe_0 -
                           ta1_z_x_y_1[i] * fe_0 + ta1_z_x_xy_0[i] * pa_x[i] - ta1_z_x_xy_1[i] * pc_x[i];

        ta1_z_xx_xz_0[i] = ta1_z_0_xz_0[i] * fe_0 - ta1_z_0_xz_1[i] * fe_0 + ta1_z_x_z_0[i] * fe_0 -
                           ta1_z_x_z_1[i] * fe_0 + ta1_z_x_xz_0[i] * pa_x[i] - ta1_z_x_xz_1[i] * pc_x[i];

        ta1_z_xx_yy_0[i] =
            ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + ta1_z_x_yy_0[i] * pa_x[i] - ta1_z_x_yy_1[i] * pc_x[i];

        ta1_z_xx_yz_0[i] =
            ta1_z_0_yz_0[i] * fe_0 - ta1_z_0_yz_1[i] * fe_0 + ta1_z_x_yz_0[i] * pa_x[i] - ta1_z_x_yz_1[i] * pc_x[i];

        ta1_z_xx_zz_0[i] =
            ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_x_zz_0[i] * pa_x[i] - ta1_z_x_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : DD

    auto ta1_z_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 78);

    auto ta1_z_xy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 79);

    auto ta1_z_xy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 80);

    auto ta1_z_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 81);

    auto ta1_z_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 82);

    auto ta1_z_xy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 83);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_z_x_xx_0,  \
                             ta1_z_x_xx_1,  \
                             ta1_z_x_xz_0,  \
                             ta1_z_x_xz_1,  \
                             ta1_z_xy_xx_0, \
                             ta1_z_xy_xy_0, \
                             ta1_z_xy_xz_0, \
                             ta1_z_xy_yy_0, \
                             ta1_z_xy_yz_0, \
                             ta1_z_xy_zz_0, \
                             ta1_z_y_xy_0,  \
                             ta1_z_y_xy_1,  \
                             ta1_z_y_y_0,   \
                             ta1_z_y_y_1,   \
                             ta1_z_y_yy_0,  \
                             ta1_z_y_yy_1,  \
                             ta1_z_y_yz_0,  \
                             ta1_z_y_yz_1,  \
                             ta1_z_y_zz_0,  \
                             ta1_z_y_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xy_xx_0[i] = ta1_z_x_xx_0[i] * pa_y[i] - ta1_z_x_xx_1[i] * pc_y[i];

        ta1_z_xy_xy_0[i] =
            ta1_z_y_y_0[i] * fe_0 - ta1_z_y_y_1[i] * fe_0 + ta1_z_y_xy_0[i] * pa_x[i] - ta1_z_y_xy_1[i] * pc_x[i];

        ta1_z_xy_xz_0[i] = ta1_z_x_xz_0[i] * pa_y[i] - ta1_z_x_xz_1[i] * pc_y[i];

        ta1_z_xy_yy_0[i] = ta1_z_y_yy_0[i] * pa_x[i] - ta1_z_y_yy_1[i] * pc_x[i];

        ta1_z_xy_yz_0[i] = ta1_z_y_yz_0[i] * pa_x[i] - ta1_z_y_yz_1[i] * pc_x[i];

        ta1_z_xy_zz_0[i] = ta1_z_y_zz_0[i] * pa_x[i] - ta1_z_y_zz_1[i] * pc_x[i];
    }

    // Set up 84-90 components of targeted buffer : DD

    auto ta1_z_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 84);

    auto ta1_z_xz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 85);

    auto ta1_z_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 86);

    auto ta1_z_xz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 87);

    auto ta1_z_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 88);

    auto ta1_z_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_z_x_xx_0,  \
                             ta1_z_x_xx_1,  \
                             ta1_z_x_xy_0,  \
                             ta1_z_x_xy_1,  \
                             ta1_z_xz_xx_0, \
                             ta1_z_xz_xy_0, \
                             ta1_z_xz_xz_0, \
                             ta1_z_xz_yy_0, \
                             ta1_z_xz_yz_0, \
                             ta1_z_xz_zz_0, \
                             ta1_z_z_xz_0,  \
                             ta1_z_z_xz_1,  \
                             ta1_z_z_yy_0,  \
                             ta1_z_z_yy_1,  \
                             ta1_z_z_yz_0,  \
                             ta1_z_z_yz_1,  \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta1_z_z_zz_0,  \
                             ta1_z_z_zz_1,  \
                             ta_x_xx_1,     \
                             ta_x_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xz_xx_0[i] = ta_x_xx_1[i] + ta1_z_x_xx_0[i] * pa_z[i] - ta1_z_x_xx_1[i] * pc_z[i];

        ta1_z_xz_xy_0[i] = ta_x_xy_1[i] + ta1_z_x_xy_0[i] * pa_z[i] - ta1_z_x_xy_1[i] * pc_z[i];

        ta1_z_xz_xz_0[i] =
            ta1_z_z_z_0[i] * fe_0 - ta1_z_z_z_1[i] * fe_0 + ta1_z_z_xz_0[i] * pa_x[i] - ta1_z_z_xz_1[i] * pc_x[i];

        ta1_z_xz_yy_0[i] = ta1_z_z_yy_0[i] * pa_x[i] - ta1_z_z_yy_1[i] * pc_x[i];

        ta1_z_xz_yz_0[i] = ta1_z_z_yz_0[i] * pa_x[i] - ta1_z_z_yz_1[i] * pc_x[i];

        ta1_z_xz_zz_0[i] = ta1_z_z_zz_0[i] * pa_x[i] - ta1_z_z_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : DD

    auto ta1_z_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 90);

    auto ta1_z_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 91);

    auto ta1_z_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 92);

    auto ta1_z_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 93);

    auto ta1_z_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 94);

    auto ta1_z_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 95);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_y_x_0,   \
                             ta1_z_y_x_1,   \
                             ta1_z_y_xx_0,  \
                             ta1_z_y_xx_1,  \
                             ta1_z_y_xy_0,  \
                             ta1_z_y_xy_1,  \
                             ta1_z_y_xz_0,  \
                             ta1_z_y_xz_1,  \
                             ta1_z_y_y_0,   \
                             ta1_z_y_y_1,   \
                             ta1_z_y_yy_0,  \
                             ta1_z_y_yy_1,  \
                             ta1_z_y_yz_0,  \
                             ta1_z_y_yz_1,  \
                             ta1_z_y_z_0,   \
                             ta1_z_y_z_1,   \
                             ta1_z_y_zz_0,  \
                             ta1_z_y_zz_1,  \
                             ta1_z_yy_xx_0, \
                             ta1_z_yy_xy_0, \
                             ta1_z_yy_xz_0, \
                             ta1_z_yy_yy_0, \
                             ta1_z_yy_yz_0, \
                             ta1_z_yy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_xx_0[i] =
            ta1_z_0_xx_0[i] * fe_0 - ta1_z_0_xx_1[i] * fe_0 + ta1_z_y_xx_0[i] * pa_y[i] - ta1_z_y_xx_1[i] * pc_y[i];

        ta1_z_yy_xy_0[i] = ta1_z_0_xy_0[i] * fe_0 - ta1_z_0_xy_1[i] * fe_0 + ta1_z_y_x_0[i] * fe_0 -
                           ta1_z_y_x_1[i] * fe_0 + ta1_z_y_xy_0[i] * pa_y[i] - ta1_z_y_xy_1[i] * pc_y[i];

        ta1_z_yy_xz_0[i] =
            ta1_z_0_xz_0[i] * fe_0 - ta1_z_0_xz_1[i] * fe_0 + ta1_z_y_xz_0[i] * pa_y[i] - ta1_z_y_xz_1[i] * pc_y[i];

        ta1_z_yy_yy_0[i] = ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + 2.0 * ta1_z_y_y_0[i] * fe_0 -
                           2.0 * ta1_z_y_y_1[i] * fe_0 + ta1_z_y_yy_0[i] * pa_y[i] - ta1_z_y_yy_1[i] * pc_y[i];

        ta1_z_yy_yz_0[i] = ta1_z_0_yz_0[i] * fe_0 - ta1_z_0_yz_1[i] * fe_0 + ta1_z_y_z_0[i] * fe_0 -
                           ta1_z_y_z_1[i] * fe_0 + ta1_z_y_yz_0[i] * pa_y[i] - ta1_z_y_yz_1[i] * pc_y[i];

        ta1_z_yy_zz_0[i] =
            ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_y_zz_0[i] * pa_y[i] - ta1_z_y_zz_1[i] * pc_y[i];
    }

    // Set up 96-102 components of targeted buffer : DD

    auto ta1_z_yz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 96);

    auto ta1_z_yz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 97);

    auto ta1_z_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 98);

    auto ta1_z_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 99);

    auto ta1_z_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 100);

    auto ta1_z_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 101);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_z_y_xy_0,  \
                             ta1_z_y_xy_1,  \
                             ta1_z_y_yy_0,  \
                             ta1_z_y_yy_1,  \
                             ta1_z_yz_xx_0, \
                             ta1_z_yz_xy_0, \
                             ta1_z_yz_xz_0, \
                             ta1_z_yz_yy_0, \
                             ta1_z_yz_yz_0, \
                             ta1_z_yz_zz_0, \
                             ta1_z_z_xx_0,  \
                             ta1_z_z_xx_1,  \
                             ta1_z_z_xz_0,  \
                             ta1_z_z_xz_1,  \
                             ta1_z_z_yz_0,  \
                             ta1_z_z_yz_1,  \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta1_z_z_zz_0,  \
                             ta1_z_z_zz_1,  \
                             ta_y_xy_1,     \
                             ta_y_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yz_xx_0[i] = ta1_z_z_xx_0[i] * pa_y[i] - ta1_z_z_xx_1[i] * pc_y[i];

        ta1_z_yz_xy_0[i] = ta_y_xy_1[i] + ta1_z_y_xy_0[i] * pa_z[i] - ta1_z_y_xy_1[i] * pc_z[i];

        ta1_z_yz_xz_0[i] = ta1_z_z_xz_0[i] * pa_y[i] - ta1_z_z_xz_1[i] * pc_y[i];

        ta1_z_yz_yy_0[i] = ta_y_yy_1[i] + ta1_z_y_yy_0[i] * pa_z[i] - ta1_z_y_yy_1[i] * pc_z[i];

        ta1_z_yz_yz_0[i] =
            ta1_z_z_z_0[i] * fe_0 - ta1_z_z_z_1[i] * fe_0 + ta1_z_z_yz_0[i] * pa_y[i] - ta1_z_z_yz_1[i] * pc_y[i];

        ta1_z_yz_zz_0[i] = ta1_z_z_zz_0[i] * pa_y[i] - ta1_z_z_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : DD

    auto ta1_z_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 102);

    auto ta1_z_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 103);

    auto ta1_z_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 104);

    auto ta1_z_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 105);

    auto ta1_z_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 106);

    auto ta1_z_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 107);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_z_x_0,   \
                             ta1_z_z_x_1,   \
                             ta1_z_z_xx_0,  \
                             ta1_z_z_xx_1,  \
                             ta1_z_z_xy_0,  \
                             ta1_z_z_xy_1,  \
                             ta1_z_z_xz_0,  \
                             ta1_z_z_xz_1,  \
                             ta1_z_z_y_0,   \
                             ta1_z_z_y_1,   \
                             ta1_z_z_yy_0,  \
                             ta1_z_z_yy_1,  \
                             ta1_z_z_yz_0,  \
                             ta1_z_z_yz_1,  \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta1_z_z_zz_0,  \
                             ta1_z_z_zz_1,  \
                             ta1_z_zz_xx_0, \
                             ta1_z_zz_xy_0, \
                             ta1_z_zz_xz_0, \
                             ta1_z_zz_yy_0, \
                             ta1_z_zz_yz_0, \
                             ta1_z_zz_zz_0, \
                             ta_z_xx_1,     \
                             ta_z_xy_1,     \
                             ta_z_xz_1,     \
                             ta_z_yy_1,     \
                             ta_z_yz_1,     \
                             ta_z_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_xx_0[i] = ta1_z_0_xx_0[i] * fe_0 - ta1_z_0_xx_1[i] * fe_0 + ta_z_xx_1[i] + ta1_z_z_xx_0[i] * pa_z[i] -
                           ta1_z_z_xx_1[i] * pc_z[i];

        ta1_z_zz_xy_0[i] = ta1_z_0_xy_0[i] * fe_0 - ta1_z_0_xy_1[i] * fe_0 + ta_z_xy_1[i] + ta1_z_z_xy_0[i] * pa_z[i] -
                           ta1_z_z_xy_1[i] * pc_z[i];

        ta1_z_zz_xz_0[i] = ta1_z_0_xz_0[i] * fe_0 - ta1_z_0_xz_1[i] * fe_0 + ta1_z_z_x_0[i] * fe_0 -
                           ta1_z_z_x_1[i] * fe_0 + ta_z_xz_1[i] + ta1_z_z_xz_0[i] * pa_z[i] - ta1_z_z_xz_1[i] * pc_z[i];

        ta1_z_zz_yy_0[i] = ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + ta_z_yy_1[i] + ta1_z_z_yy_0[i] * pa_z[i] -
                           ta1_z_z_yy_1[i] * pc_z[i];

        ta1_z_zz_yz_0[i] = ta1_z_0_yz_0[i] * fe_0 - ta1_z_0_yz_1[i] * fe_0 + ta1_z_z_y_0[i] * fe_0 -
                           ta1_z_z_y_1[i] * fe_0 + ta_z_yz_1[i] + ta1_z_z_yz_0[i] * pa_z[i] - ta1_z_z_yz_1[i] * pc_z[i];

        ta1_z_zz_zz_0[i] = ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + 2.0 * ta1_z_z_z_0[i] * fe_0 -
                           2.0 * ta1_z_z_z_1[i] * fe_0 + ta_z_zz_1[i] + ta1_z_z_zz_0[i] * pa_z[i] -
                           ta1_z_z_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
