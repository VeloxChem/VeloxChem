#include "NuclearPotentialGeom020PrimRecDD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_dd(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_dd,
                                        const size_t idx_npot_geom_020_0_sd,
                                        const size_t idx_npot_geom_020_1_sd,
                                        const size_t idx_npot_geom_020_0_pp,
                                        const size_t idx_npot_geom_020_1_pp,
                                        const size_t idx_npot_geom_010_1_pd,
                                        const size_t idx_npot_geom_020_0_pd,
                                        const size_t idx_npot_geom_020_1_pd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
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

    auto ta2_xx_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd);

    auto ta2_xx_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 1);

    auto ta2_xx_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 2);

    auto ta2_xx_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 3);

    auto ta2_xx_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 4);

    auto ta2_xx_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 5);

    auto ta2_xy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 6);

    auto ta2_xy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 7);

    auto ta2_xy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 8);

    auto ta2_xy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 9);

    auto ta2_xy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 10);

    auto ta2_xy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 11);

    auto ta2_xz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 12);

    auto ta2_xz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 13);

    auto ta2_xz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 14);

    auto ta2_xz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 15);

    auto ta2_xz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 16);

    auto ta2_xz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 17);

    auto ta2_yy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 18);

    auto ta2_yy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 19);

    auto ta2_yy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 20);

    auto ta2_yy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 21);

    auto ta2_yy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 22);

    auto ta2_yy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 23);

    auto ta2_yz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 24);

    auto ta2_yz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 25);

    auto ta2_yz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 26);

    auto ta2_yz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 27);

    auto ta2_yz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 28);

    auto ta2_yz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 29);

    auto ta2_zz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 30);

    auto ta2_zz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 31);

    auto ta2_zz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 32);

    auto ta2_zz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 33);

    auto ta2_zz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 34);

    auto ta2_zz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 35);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd);

    auto ta2_xx_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 1);

    auto ta2_xx_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 2);

    auto ta2_xx_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 3);

    auto ta2_xx_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 4);

    auto ta2_xx_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 5);

    auto ta2_xy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 6);

    auto ta2_xy_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 7);

    auto ta2_xy_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 8);

    auto ta2_xy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 9);

    auto ta2_xy_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 10);

    auto ta2_xy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 11);

    auto ta2_xz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 12);

    auto ta2_xz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 13);

    auto ta2_xz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 14);

    auto ta2_xz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 15);

    auto ta2_xz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 16);

    auto ta2_xz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 17);

    auto ta2_yy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 18);

    auto ta2_yy_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 19);

    auto ta2_yy_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 20);

    auto ta2_yy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 21);

    auto ta2_yy_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 22);

    auto ta2_yy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 23);

    auto ta2_yz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 24);

    auto ta2_yz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 25);

    auto ta2_yz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 26);

    auto ta2_yz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 27);

    auto ta2_yz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 28);

    auto ta2_yz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 29);

    auto ta2_zz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 30);

    auto ta2_zz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 31);

    auto ta2_zz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 32);

    auto ta2_zz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 33);

    auto ta2_zz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 34);

    auto ta2_zz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 35);

    // Set up components of auxiliary buffer : PP

    auto ta2_xx_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp);

    auto ta2_xx_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 1);

    auto ta2_xx_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 2);

    auto ta2_xx_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 3);

    auto ta2_xx_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 4);

    auto ta2_xx_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 5);

    auto ta2_xx_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 6);

    auto ta2_xx_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 7);

    auto ta2_xx_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 8);

    auto ta2_xy_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 9);

    auto ta2_xy_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 10);

    auto ta2_xy_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 11);

    auto ta2_xy_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 12);

    auto ta2_xy_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 13);

    auto ta2_xy_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 14);

    auto ta2_xy_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 15);

    auto ta2_xy_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 16);

    auto ta2_xy_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 17);

    auto ta2_xz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 18);

    auto ta2_xz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 19);

    auto ta2_xz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 20);

    auto ta2_xz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 21);

    auto ta2_xz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 22);

    auto ta2_xz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 23);

    auto ta2_xz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 24);

    auto ta2_xz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 25);

    auto ta2_xz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 26);

    auto ta2_yy_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 27);

    auto ta2_yy_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 28);

    auto ta2_yy_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 29);

    auto ta2_yy_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 30);

    auto ta2_yy_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 31);

    auto ta2_yy_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 32);

    auto ta2_yy_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 33);

    auto ta2_yy_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 34);

    auto ta2_yy_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 35);

    auto ta2_yz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 36);

    auto ta2_yz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 37);

    auto ta2_yz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 38);

    auto ta2_yz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 39);

    auto ta2_yz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 40);

    auto ta2_yz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 41);

    auto ta2_yz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 42);

    auto ta2_yz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 43);

    auto ta2_yz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 44);

    auto ta2_zz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 45);

    auto ta2_zz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 46);

    auto ta2_zz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 47);

    auto ta2_zz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 48);

    auto ta2_zz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 49);

    auto ta2_zz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 50);

    auto ta2_zz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 51);

    auto ta2_zz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 52);

    auto ta2_zz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 53);

    // Set up components of auxiliary buffer : PP

    auto ta2_xx_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp);

    auto ta2_xx_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 1);

    auto ta2_xx_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 2);

    auto ta2_xx_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 3);

    auto ta2_xx_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 4);

    auto ta2_xx_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 5);

    auto ta2_xx_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 6);

    auto ta2_xx_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 7);

    auto ta2_xx_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 8);

    auto ta2_xy_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 9);

    auto ta2_xy_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 10);

    auto ta2_xy_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 11);

    auto ta2_xy_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 12);

    auto ta2_xy_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 13);

    auto ta2_xy_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 14);

    auto ta2_xy_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 15);

    auto ta2_xy_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 16);

    auto ta2_xy_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 17);

    auto ta2_xz_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 18);

    auto ta2_xz_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 19);

    auto ta2_xz_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 20);

    auto ta2_xz_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 21);

    auto ta2_xz_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 22);

    auto ta2_xz_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 23);

    auto ta2_xz_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 24);

    auto ta2_xz_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 25);

    auto ta2_xz_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 26);

    auto ta2_yy_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 27);

    auto ta2_yy_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 28);

    auto ta2_yy_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 29);

    auto ta2_yy_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 30);

    auto ta2_yy_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 31);

    auto ta2_yy_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 32);

    auto ta2_yy_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 33);

    auto ta2_yy_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 34);

    auto ta2_yy_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 35);

    auto ta2_yz_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 36);

    auto ta2_yz_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 37);

    auto ta2_yz_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 38);

    auto ta2_yz_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 39);

    auto ta2_yz_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 40);

    auto ta2_yz_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 41);

    auto ta2_yz_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 42);

    auto ta2_yz_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 43);

    auto ta2_yz_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 44);

    auto ta2_zz_x_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 45);

    auto ta2_zz_x_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 46);

    auto ta2_zz_x_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 47);

    auto ta2_zz_y_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 48);

    auto ta2_zz_y_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 49);

    auto ta2_zz_y_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 50);

    auto ta2_zz_z_x_1 = pbuffer.data(idx_npot_geom_020_1_pp + 51);

    auto ta2_zz_z_y_1 = pbuffer.data(idx_npot_geom_020_1_pp + 52);

    auto ta2_zz_z_z_1 = pbuffer.data(idx_npot_geom_020_1_pp + 53);

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

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd);

    auto ta2_xx_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 1);

    auto ta2_xx_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 2);

    auto ta2_xx_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 3);

    auto ta2_xx_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 4);

    auto ta2_xx_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 5);

    auto ta2_xx_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 6);

    auto ta2_xx_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 7);

    auto ta2_xx_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 8);

    auto ta2_xx_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 9);

    auto ta2_xx_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 10);

    auto ta2_xx_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 11);

    auto ta2_xx_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 12);

    auto ta2_xx_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 13);

    auto ta2_xx_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 14);

    auto ta2_xx_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 15);

    auto ta2_xx_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 16);

    auto ta2_xx_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 17);

    auto ta2_xy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 18);

    auto ta2_xy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 19);

    auto ta2_xy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 20);

    auto ta2_xy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 21);

    auto ta2_xy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 22);

    auto ta2_xy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 23);

    auto ta2_xy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 24);

    auto ta2_xy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 25);

    auto ta2_xy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 26);

    auto ta2_xy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 27);

    auto ta2_xy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 28);

    auto ta2_xy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 29);

    auto ta2_xy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 30);

    auto ta2_xy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 31);

    auto ta2_xy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 32);

    auto ta2_xy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 33);

    auto ta2_xy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 34);

    auto ta2_xy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 35);

    auto ta2_xz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 36);

    auto ta2_xz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 37);

    auto ta2_xz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 38);

    auto ta2_xz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 39);

    auto ta2_xz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 40);

    auto ta2_xz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 41);

    auto ta2_xz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 42);

    auto ta2_xz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 43);

    auto ta2_xz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 44);

    auto ta2_xz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 45);

    auto ta2_xz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 46);

    auto ta2_xz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 47);

    auto ta2_xz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 48);

    auto ta2_xz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 49);

    auto ta2_xz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 50);

    auto ta2_xz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 51);

    auto ta2_xz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 52);

    auto ta2_xz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 53);

    auto ta2_yy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 54);

    auto ta2_yy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 55);

    auto ta2_yy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 56);

    auto ta2_yy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 57);

    auto ta2_yy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 58);

    auto ta2_yy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 59);

    auto ta2_yy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 60);

    auto ta2_yy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 61);

    auto ta2_yy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 62);

    auto ta2_yy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 63);

    auto ta2_yy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 64);

    auto ta2_yy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 65);

    auto ta2_yy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 66);

    auto ta2_yy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 67);

    auto ta2_yy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 68);

    auto ta2_yy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 69);

    auto ta2_yy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 70);

    auto ta2_yy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 71);

    auto ta2_yz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 72);

    auto ta2_yz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 73);

    auto ta2_yz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 74);

    auto ta2_yz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 75);

    auto ta2_yz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 76);

    auto ta2_yz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 77);

    auto ta2_yz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 78);

    auto ta2_yz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 79);

    auto ta2_yz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 80);

    auto ta2_yz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 81);

    auto ta2_yz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 82);

    auto ta2_yz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 83);

    auto ta2_yz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 84);

    auto ta2_yz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 85);

    auto ta2_yz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 86);

    auto ta2_yz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 87);

    auto ta2_yz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 88);

    auto ta2_yz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 89);

    auto ta2_zz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 90);

    auto ta2_zz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 91);

    auto ta2_zz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 92);

    auto ta2_zz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 93);

    auto ta2_zz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 94);

    auto ta2_zz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 95);

    auto ta2_zz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 96);

    auto ta2_zz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 97);

    auto ta2_zz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 98);

    auto ta2_zz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 99);

    auto ta2_zz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 100);

    auto ta2_zz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 101);

    auto ta2_zz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 102);

    auto ta2_zz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 103);

    auto ta2_zz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 104);

    auto ta2_zz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 105);

    auto ta2_zz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 106);

    auto ta2_zz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 107);

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd);

    auto ta2_xx_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 1);

    auto ta2_xx_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 2);

    auto ta2_xx_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 3);

    auto ta2_xx_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 4);

    auto ta2_xx_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 5);

    auto ta2_xx_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 6);

    auto ta2_xx_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 7);

    auto ta2_xx_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 8);

    auto ta2_xx_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 9);

    auto ta2_xx_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 10);

    auto ta2_xx_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 11);

    auto ta2_xx_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 12);

    auto ta2_xx_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 13);

    auto ta2_xx_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 14);

    auto ta2_xx_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 15);

    auto ta2_xx_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 16);

    auto ta2_xx_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 17);

    auto ta2_xy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 18);

    auto ta2_xy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 19);

    auto ta2_xy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 20);

    auto ta2_xy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 21);

    auto ta2_xy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 22);

    auto ta2_xy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 23);

    auto ta2_xy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 24);

    auto ta2_xy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 25);

    auto ta2_xy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 26);

    auto ta2_xy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 27);

    auto ta2_xy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 28);

    auto ta2_xy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 29);

    auto ta2_xy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 30);

    auto ta2_xy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 31);

    auto ta2_xy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 32);

    auto ta2_xy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 33);

    auto ta2_xy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 34);

    auto ta2_xy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 35);

    auto ta2_xz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 36);

    auto ta2_xz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 37);

    auto ta2_xz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 38);

    auto ta2_xz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 39);

    auto ta2_xz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 40);

    auto ta2_xz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 41);

    auto ta2_xz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 42);

    auto ta2_xz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 43);

    auto ta2_xz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 44);

    auto ta2_xz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 45);

    auto ta2_xz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 46);

    auto ta2_xz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 47);

    auto ta2_xz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 48);

    auto ta2_xz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 49);

    auto ta2_xz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 50);

    auto ta2_xz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 51);

    auto ta2_xz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 52);

    auto ta2_xz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 53);

    auto ta2_yy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 54);

    auto ta2_yy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 55);

    auto ta2_yy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 56);

    auto ta2_yy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 57);

    auto ta2_yy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 58);

    auto ta2_yy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 59);

    auto ta2_yy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 60);

    auto ta2_yy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 61);

    auto ta2_yy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 62);

    auto ta2_yy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 63);

    auto ta2_yy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 64);

    auto ta2_yy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 65);

    auto ta2_yy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 66);

    auto ta2_yy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 67);

    auto ta2_yy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 68);

    auto ta2_yy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 69);

    auto ta2_yy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 70);

    auto ta2_yy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 71);

    auto ta2_yz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 72);

    auto ta2_yz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 73);

    auto ta2_yz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 74);

    auto ta2_yz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 75);

    auto ta2_yz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 76);

    auto ta2_yz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 77);

    auto ta2_yz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 78);

    auto ta2_yz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 79);

    auto ta2_yz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 80);

    auto ta2_yz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 81);

    auto ta2_yz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 82);

    auto ta2_yz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 83);

    auto ta2_yz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 84);

    auto ta2_yz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 85);

    auto ta2_yz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 86);

    auto ta2_yz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 87);

    auto ta2_yz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 88);

    auto ta2_yz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 89);

    auto ta2_zz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 90);

    auto ta2_zz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 91);

    auto ta2_zz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 92);

    auto ta2_zz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 93);

    auto ta2_zz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 94);

    auto ta2_zz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 95);

    auto ta2_zz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 96);

    auto ta2_zz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 97);

    auto ta2_zz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 98);

    auto ta2_zz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 99);

    auto ta2_zz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 100);

    auto ta2_zz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 101);

    auto ta2_zz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 102);

    auto ta2_zz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 103);

    auto ta2_zz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 104);

    auto ta2_zz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 105);

    auto ta2_zz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 106);

    auto ta2_zz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 107);

    // Set up 0-6 components of targeted buffer : DD

    auto ta2_xx_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd);

    auto ta2_xx_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 1);

    auto ta2_xx_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 2);

    auto ta2_xx_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 3);

    auto ta2_xx_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 4);

    auto ta2_xx_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 5);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_x_xx_1, ta1_x_x_xy_1, ta1_x_x_xz_1, ta1_x_x_yy_1, ta1_x_x_yz_1, ta1_x_x_zz_1, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xy_0, ta2_xx_x_xy_1, ta2_xx_x_xz_0, ta2_xx_x_xz_1, ta2_xx_x_y_0, ta2_xx_x_y_1, ta2_xx_x_yy_0, ta2_xx_x_yy_1, ta2_xx_x_yz_0, ta2_xx_x_yz_1, ta2_xx_x_z_0, ta2_xx_x_z_1, ta2_xx_x_zz_0, ta2_xx_x_zz_1, ta2_xx_xx_xx_0, ta2_xx_xx_xy_0, ta2_xx_xx_xz_0, ta2_xx_xx_yy_0, ta2_xx_xx_yz_0, ta2_xx_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xx_xx_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + 2.0 * ta2_xx_x_x_0[i] * fe_0 - 2.0 * ta2_xx_x_x_1[i] * fe_0 + 2.0 * ta1_x_x_xx_1[i] + ta2_xx_x_xx_0[i] * pa_x[i] - ta2_xx_x_xx_1[i] * pc_x[i];

        ta2_xx_xx_xy_0[i] = ta2_xx_0_xy_0[i] * fe_0 - ta2_xx_0_xy_1[i] * fe_0 + ta2_xx_x_y_0[i] * fe_0 - ta2_xx_x_y_1[i] * fe_0 + 2.0 * ta1_x_x_xy_1[i] + ta2_xx_x_xy_0[i] * pa_x[i] - ta2_xx_x_xy_1[i] * pc_x[i];

        ta2_xx_xx_xz_0[i] = ta2_xx_0_xz_0[i] * fe_0 - ta2_xx_0_xz_1[i] * fe_0 + ta2_xx_x_z_0[i] * fe_0 - ta2_xx_x_z_1[i] * fe_0 + 2.0 * ta1_x_x_xz_1[i] + ta2_xx_x_xz_0[i] * pa_x[i] - ta2_xx_x_xz_1[i] * pc_x[i];

        ta2_xx_xx_yy_0[i] = ta2_xx_0_yy_0[i] * fe_0 - ta2_xx_0_yy_1[i] * fe_0 + 2.0 * ta1_x_x_yy_1[i] + ta2_xx_x_yy_0[i] * pa_x[i] - ta2_xx_x_yy_1[i] * pc_x[i];

        ta2_xx_xx_yz_0[i] = ta2_xx_0_yz_0[i] * fe_0 - ta2_xx_0_yz_1[i] * fe_0 + 2.0 * ta1_x_x_yz_1[i] + ta2_xx_x_yz_0[i] * pa_x[i] - ta2_xx_x_yz_1[i] * pc_x[i];

        ta2_xx_xx_zz_0[i] = ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + 2.0 * ta1_x_x_zz_1[i] + ta2_xx_x_zz_0[i] * pa_x[i] - ta2_xx_x_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto ta2_xx_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 6);

    auto ta2_xx_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 7);

    auto ta2_xx_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 8);

    auto ta2_xx_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 9);

    auto ta2_xx_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 10);

    auto ta2_xx_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 11);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_y_yy_1, ta1_x_y_yz_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xy_0, ta2_xx_x_xy_1, ta2_xx_x_xz_0, ta2_xx_x_xz_1, ta2_xx_x_zz_0, ta2_xx_x_zz_1, ta2_xx_xy_xx_0, ta2_xx_xy_xy_0, ta2_xx_xy_xz_0, ta2_xx_xy_yy_0, ta2_xx_xy_yz_0, ta2_xx_xy_zz_0, ta2_xx_y_yy_0, ta2_xx_y_yy_1, ta2_xx_y_yz_0, ta2_xx_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xy_xx_0[i] = ta2_xx_x_xx_0[i] * pa_y[i] - ta2_xx_x_xx_1[i] * pc_y[i];

        ta2_xx_xy_xy_0[i] = ta2_xx_x_x_0[i] * fe_0 - ta2_xx_x_x_1[i] * fe_0 + ta2_xx_x_xy_0[i] * pa_y[i] - ta2_xx_x_xy_1[i] * pc_y[i];

        ta2_xx_xy_xz_0[i] = ta2_xx_x_xz_0[i] * pa_y[i] - ta2_xx_x_xz_1[i] * pc_y[i];

        ta2_xx_xy_yy_0[i] = 2.0 * ta1_x_y_yy_1[i] + ta2_xx_y_yy_0[i] * pa_x[i] - ta2_xx_y_yy_1[i] * pc_x[i];

        ta2_xx_xy_yz_0[i] = 2.0 * ta1_x_y_yz_1[i] + ta2_xx_y_yz_0[i] * pa_x[i] - ta2_xx_y_yz_1[i] * pc_x[i];

        ta2_xx_xy_zz_0[i] = ta2_xx_x_zz_0[i] * pa_y[i] - ta2_xx_x_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto ta2_xx_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 12);

    auto ta2_xx_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 13);

    auto ta2_xx_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 14);

    auto ta2_xx_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 15);

    auto ta2_xx_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 16);

    auto ta2_xx_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 17);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_z_yz_1, ta1_x_z_zz_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xy_0, ta2_xx_x_xy_1, ta2_xx_x_xz_0, ta2_xx_x_xz_1, ta2_xx_x_yy_0, ta2_xx_x_yy_1, ta2_xx_xz_xx_0, ta2_xx_xz_xy_0, ta2_xx_xz_xz_0, ta2_xx_xz_yy_0, ta2_xx_xz_yz_0, ta2_xx_xz_zz_0, ta2_xx_z_yz_0, ta2_xx_z_yz_1, ta2_xx_z_zz_0, ta2_xx_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xz_xx_0[i] = ta2_xx_x_xx_0[i] * pa_z[i] - ta2_xx_x_xx_1[i] * pc_z[i];

        ta2_xx_xz_xy_0[i] = ta2_xx_x_xy_0[i] * pa_z[i] - ta2_xx_x_xy_1[i] * pc_z[i];

        ta2_xx_xz_xz_0[i] = ta2_xx_x_x_0[i] * fe_0 - ta2_xx_x_x_1[i] * fe_0 + ta2_xx_x_xz_0[i] * pa_z[i] - ta2_xx_x_xz_1[i] * pc_z[i];

        ta2_xx_xz_yy_0[i] = ta2_xx_x_yy_0[i] * pa_z[i] - ta2_xx_x_yy_1[i] * pc_z[i];

        ta2_xx_xz_yz_0[i] = 2.0 * ta1_x_z_yz_1[i] + ta2_xx_z_yz_0[i] * pa_x[i] - ta2_xx_z_yz_1[i] * pc_x[i];

        ta2_xx_xz_zz_0[i] = 2.0 * ta1_x_z_zz_1[i] + ta2_xx_z_zz_0[i] * pa_x[i] - ta2_xx_z_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto ta2_xx_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 18);

    auto ta2_xx_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 19);

    auto ta2_xx_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 20);

    auto ta2_xx_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 21);

    auto ta2_xx_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 22);

    auto ta2_xx_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 23);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_y_x_0, ta2_xx_y_x_1, ta2_xx_y_xx_0, ta2_xx_y_xx_1, ta2_xx_y_xy_0, ta2_xx_y_xy_1, ta2_xx_y_xz_0, ta2_xx_y_xz_1, ta2_xx_y_y_0, ta2_xx_y_y_1, ta2_xx_y_yy_0, ta2_xx_y_yy_1, ta2_xx_y_yz_0, ta2_xx_y_yz_1, ta2_xx_y_z_0, ta2_xx_y_z_1, ta2_xx_y_zz_0, ta2_xx_y_zz_1, ta2_xx_yy_xx_0, ta2_xx_yy_xy_0, ta2_xx_yy_xz_0, ta2_xx_yy_yy_0, ta2_xx_yy_yz_0, ta2_xx_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yy_xx_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_y_xx_0[i] * pa_y[i] - ta2_xx_y_xx_1[i] * pc_y[i];

        ta2_xx_yy_xy_0[i] = ta2_xx_0_xy_0[i] * fe_0 - ta2_xx_0_xy_1[i] * fe_0 + ta2_xx_y_x_0[i] * fe_0 - ta2_xx_y_x_1[i] * fe_0 + ta2_xx_y_xy_0[i] * pa_y[i] - ta2_xx_y_xy_1[i] * pc_y[i];

        ta2_xx_yy_xz_0[i] = ta2_xx_0_xz_0[i] * fe_0 - ta2_xx_0_xz_1[i] * fe_0 + ta2_xx_y_xz_0[i] * pa_y[i] - ta2_xx_y_xz_1[i] * pc_y[i];

        ta2_xx_yy_yy_0[i] = ta2_xx_0_yy_0[i] * fe_0 - ta2_xx_0_yy_1[i] * fe_0 + 2.0 * ta2_xx_y_y_0[i] * fe_0 - 2.0 * ta2_xx_y_y_1[i] * fe_0 + ta2_xx_y_yy_0[i] * pa_y[i] - ta2_xx_y_yy_1[i] * pc_y[i];

        ta2_xx_yy_yz_0[i] = ta2_xx_0_yz_0[i] * fe_0 - ta2_xx_0_yz_1[i] * fe_0 + ta2_xx_y_z_0[i] * fe_0 - ta2_xx_y_z_1[i] * fe_0 + ta2_xx_y_yz_0[i] * pa_y[i] - ta2_xx_y_yz_1[i] * pc_y[i];

        ta2_xx_yy_zz_0[i] = ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + ta2_xx_y_zz_0[i] * pa_y[i] - ta2_xx_y_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto ta2_xx_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 24);

    auto ta2_xx_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 25);

    auto ta2_xx_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 26);

    auto ta2_xx_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 27);

    auto ta2_xx_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 28);

    auto ta2_xx_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 29);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_y_xy_0, ta2_xx_y_xy_1, ta2_xx_y_yy_0, ta2_xx_y_yy_1, ta2_xx_yz_xx_0, ta2_xx_yz_xy_0, ta2_xx_yz_xz_0, ta2_xx_yz_yy_0, ta2_xx_yz_yz_0, ta2_xx_yz_zz_0, ta2_xx_z_xx_0, ta2_xx_z_xx_1, ta2_xx_z_xz_0, ta2_xx_z_xz_1, ta2_xx_z_yz_0, ta2_xx_z_yz_1, ta2_xx_z_z_0, ta2_xx_z_z_1, ta2_xx_z_zz_0, ta2_xx_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yz_xx_0[i] = ta2_xx_z_xx_0[i] * pa_y[i] - ta2_xx_z_xx_1[i] * pc_y[i];

        ta2_xx_yz_xy_0[i] = ta2_xx_y_xy_0[i] * pa_z[i] - ta2_xx_y_xy_1[i] * pc_z[i];

        ta2_xx_yz_xz_0[i] = ta2_xx_z_xz_0[i] * pa_y[i] - ta2_xx_z_xz_1[i] * pc_y[i];

        ta2_xx_yz_yy_0[i] = ta2_xx_y_yy_0[i] * pa_z[i] - ta2_xx_y_yy_1[i] * pc_z[i];

        ta2_xx_yz_yz_0[i] = ta2_xx_z_z_0[i] * fe_0 - ta2_xx_z_z_1[i] * fe_0 + ta2_xx_z_yz_0[i] * pa_y[i] - ta2_xx_z_yz_1[i] * pc_y[i];

        ta2_xx_yz_zz_0[i] = ta2_xx_z_zz_0[i] * pa_y[i] - ta2_xx_z_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto ta2_xx_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 30);

    auto ta2_xx_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 31);

    auto ta2_xx_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 32);

    auto ta2_xx_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 33);

    auto ta2_xx_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 34);

    auto ta2_xx_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 35);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_z_x_0, ta2_xx_z_x_1, ta2_xx_z_xx_0, ta2_xx_z_xx_1, ta2_xx_z_xy_0, ta2_xx_z_xy_1, ta2_xx_z_xz_0, ta2_xx_z_xz_1, ta2_xx_z_y_0, ta2_xx_z_y_1, ta2_xx_z_yy_0, ta2_xx_z_yy_1, ta2_xx_z_yz_0, ta2_xx_z_yz_1, ta2_xx_z_z_0, ta2_xx_z_z_1, ta2_xx_z_zz_0, ta2_xx_z_zz_1, ta2_xx_zz_xx_0, ta2_xx_zz_xy_0, ta2_xx_zz_xz_0, ta2_xx_zz_yy_0, ta2_xx_zz_yz_0, ta2_xx_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zz_xx_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_z_xx_0[i] * pa_z[i] - ta2_xx_z_xx_1[i] * pc_z[i];

        ta2_xx_zz_xy_0[i] = ta2_xx_0_xy_0[i] * fe_0 - ta2_xx_0_xy_1[i] * fe_0 + ta2_xx_z_xy_0[i] * pa_z[i] - ta2_xx_z_xy_1[i] * pc_z[i];

        ta2_xx_zz_xz_0[i] = ta2_xx_0_xz_0[i] * fe_0 - ta2_xx_0_xz_1[i] * fe_0 + ta2_xx_z_x_0[i] * fe_0 - ta2_xx_z_x_1[i] * fe_0 + ta2_xx_z_xz_0[i] * pa_z[i] - ta2_xx_z_xz_1[i] * pc_z[i];

        ta2_xx_zz_yy_0[i] = ta2_xx_0_yy_0[i] * fe_0 - ta2_xx_0_yy_1[i] * fe_0 + ta2_xx_z_yy_0[i] * pa_z[i] - ta2_xx_z_yy_1[i] * pc_z[i];

        ta2_xx_zz_yz_0[i] = ta2_xx_0_yz_0[i] * fe_0 - ta2_xx_0_yz_1[i] * fe_0 + ta2_xx_z_y_0[i] * fe_0 - ta2_xx_z_y_1[i] * fe_0 + ta2_xx_z_yz_0[i] * pa_z[i] - ta2_xx_z_yz_1[i] * pc_z[i];

        ta2_xx_zz_zz_0[i] = ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + 2.0 * ta2_xx_z_z_0[i] * fe_0 - 2.0 * ta2_xx_z_z_1[i] * fe_0 + ta2_xx_z_zz_0[i] * pa_z[i] - ta2_xx_z_zz_1[i] * pc_z[i];
    }

    // Set up 36-42 components of targeted buffer : DD

    auto ta2_xy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 36);

    auto ta2_xy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 37);

    auto ta2_xy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 38);

    auto ta2_xy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 39);

    auto ta2_xy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 40);

    auto ta2_xy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 41);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_x_xx_1, ta1_y_x_xy_1, ta1_y_x_xz_1, ta1_y_x_yy_1, ta1_y_x_yz_1, ta1_y_x_zz_1, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_x_x_0, ta2_xy_x_x_1, ta2_xy_x_xx_0, ta2_xy_x_xx_1, ta2_xy_x_xy_0, ta2_xy_x_xy_1, ta2_xy_x_xz_0, ta2_xy_x_xz_1, ta2_xy_x_y_0, ta2_xy_x_y_1, ta2_xy_x_yy_0, ta2_xy_x_yy_1, ta2_xy_x_yz_0, ta2_xy_x_yz_1, ta2_xy_x_z_0, ta2_xy_x_z_1, ta2_xy_x_zz_0, ta2_xy_x_zz_1, ta2_xy_xx_xx_0, ta2_xy_xx_xy_0, ta2_xy_xx_xz_0, ta2_xy_xx_yy_0, ta2_xy_xx_yz_0, ta2_xy_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xx_xx_0[i] = ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + 2.0 * ta2_xy_x_x_0[i] * fe_0 - 2.0 * ta2_xy_x_x_1[i] * fe_0 + ta1_y_x_xx_1[i] + ta2_xy_x_xx_0[i] * pa_x[i] - ta2_xy_x_xx_1[i] * pc_x[i];

        ta2_xy_xx_xy_0[i] = ta2_xy_0_xy_0[i] * fe_0 - ta2_xy_0_xy_1[i] * fe_0 + ta2_xy_x_y_0[i] * fe_0 - ta2_xy_x_y_1[i] * fe_0 + ta1_y_x_xy_1[i] + ta2_xy_x_xy_0[i] * pa_x[i] - ta2_xy_x_xy_1[i] * pc_x[i];

        ta2_xy_xx_xz_0[i] = ta2_xy_0_xz_0[i] * fe_0 - ta2_xy_0_xz_1[i] * fe_0 + ta2_xy_x_z_0[i] * fe_0 - ta2_xy_x_z_1[i] * fe_0 + ta1_y_x_xz_1[i] + ta2_xy_x_xz_0[i] * pa_x[i] - ta2_xy_x_xz_1[i] * pc_x[i];

        ta2_xy_xx_yy_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta1_y_x_yy_1[i] + ta2_xy_x_yy_0[i] * pa_x[i] - ta2_xy_x_yy_1[i] * pc_x[i];

        ta2_xy_xx_yz_0[i] = ta2_xy_0_yz_0[i] * fe_0 - ta2_xy_0_yz_1[i] * fe_0 + ta1_y_x_yz_1[i] + ta2_xy_x_yz_0[i] * pa_x[i] - ta2_xy_x_yz_1[i] * pc_x[i];

        ta2_xy_xx_zz_0[i] = ta2_xy_0_zz_0[i] * fe_0 - ta2_xy_0_zz_1[i] * fe_0 + ta1_y_x_zz_1[i] + ta2_xy_x_zz_0[i] * pa_x[i] - ta2_xy_x_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : DD

    auto ta2_xy_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 42);

    auto ta2_xy_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 43);

    auto ta2_xy_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 44);

    auto ta2_xy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 45);

    auto ta2_xy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 46);

    auto ta2_xy_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 47);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_x_xx_1, ta1_x_x_xz_1, ta1_y_y_xy_1, ta1_y_y_yy_1, ta1_y_y_yz_1, ta1_y_y_zz_1, ta2_xy_x_xx_0, ta2_xy_x_xx_1, ta2_xy_x_xz_0, ta2_xy_x_xz_1, ta2_xy_xy_xx_0, ta2_xy_xy_xy_0, ta2_xy_xy_xz_0, ta2_xy_xy_yy_0, ta2_xy_xy_yz_0, ta2_xy_xy_zz_0, ta2_xy_y_xy_0, ta2_xy_y_xy_1, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_y_yz_0, ta2_xy_y_yz_1, ta2_xy_y_zz_0, ta2_xy_y_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xy_xx_0[i] = ta1_x_x_xx_1[i] + ta2_xy_x_xx_0[i] * pa_y[i] - ta2_xy_x_xx_1[i] * pc_y[i];

        ta2_xy_xy_xy_0[i] = ta2_xy_y_y_0[i] * fe_0 - ta2_xy_y_y_1[i] * fe_0 + ta1_y_y_xy_1[i] + ta2_xy_y_xy_0[i] * pa_x[i] - ta2_xy_y_xy_1[i] * pc_x[i];

        ta2_xy_xy_xz_0[i] = ta1_x_x_xz_1[i] + ta2_xy_x_xz_0[i] * pa_y[i] - ta2_xy_x_xz_1[i] * pc_y[i];

        ta2_xy_xy_yy_0[i] = ta1_y_y_yy_1[i] + ta2_xy_y_yy_0[i] * pa_x[i] - ta2_xy_y_yy_1[i] * pc_x[i];

        ta2_xy_xy_yz_0[i] = ta1_y_y_yz_1[i] + ta2_xy_y_yz_0[i] * pa_x[i] - ta2_xy_y_yz_1[i] * pc_x[i];

        ta2_xy_xy_zz_0[i] = ta1_y_y_zz_1[i] + ta2_xy_y_zz_0[i] * pa_x[i] - ta2_xy_y_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : DD

    auto ta2_xy_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 48);

    auto ta2_xy_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 49);

    auto ta2_xy_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 50);

    auto ta2_xy_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 51);

    auto ta2_xy_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 52);

    auto ta2_xy_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 53);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_z_yz_1, ta1_y_z_zz_1, ta2_xy_x_x_0, ta2_xy_x_x_1, ta2_xy_x_xx_0, ta2_xy_x_xx_1, ta2_xy_x_xy_0, ta2_xy_x_xy_1, ta2_xy_x_xz_0, ta2_xy_x_xz_1, ta2_xy_x_yy_0, ta2_xy_x_yy_1, ta2_xy_xz_xx_0, ta2_xy_xz_xy_0, ta2_xy_xz_xz_0, ta2_xy_xz_yy_0, ta2_xy_xz_yz_0, ta2_xy_xz_zz_0, ta2_xy_z_yz_0, ta2_xy_z_yz_1, ta2_xy_z_zz_0, ta2_xy_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xz_xx_0[i] = ta2_xy_x_xx_0[i] * pa_z[i] - ta2_xy_x_xx_1[i] * pc_z[i];

        ta2_xy_xz_xy_0[i] = ta2_xy_x_xy_0[i] * pa_z[i] - ta2_xy_x_xy_1[i] * pc_z[i];

        ta2_xy_xz_xz_0[i] = ta2_xy_x_x_0[i] * fe_0 - ta2_xy_x_x_1[i] * fe_0 + ta2_xy_x_xz_0[i] * pa_z[i] - ta2_xy_x_xz_1[i] * pc_z[i];

        ta2_xy_xz_yy_0[i] = ta2_xy_x_yy_0[i] * pa_z[i] - ta2_xy_x_yy_1[i] * pc_z[i];

        ta2_xy_xz_yz_0[i] = ta1_y_z_yz_1[i] + ta2_xy_z_yz_0[i] * pa_x[i] - ta2_xy_z_yz_1[i] * pc_x[i];

        ta2_xy_xz_zz_0[i] = ta1_y_z_zz_1[i] + ta2_xy_z_zz_0[i] * pa_x[i] - ta2_xy_z_zz_1[i] * pc_x[i];
    }

    // Set up 54-60 components of targeted buffer : DD

    auto ta2_xy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 54);

    auto ta2_xy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 55);

    auto ta2_xy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 56);

    auto ta2_xy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 57);

    auto ta2_xy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 58);

    auto ta2_xy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 59);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_y_xx_1, ta1_x_y_xy_1, ta1_x_y_xz_1, ta1_x_y_yy_1, ta1_x_y_yz_1, ta1_x_y_zz_1, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_y_x_0, ta2_xy_y_x_1, ta2_xy_y_xx_0, ta2_xy_y_xx_1, ta2_xy_y_xy_0, ta2_xy_y_xy_1, ta2_xy_y_xz_0, ta2_xy_y_xz_1, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_y_yz_0, ta2_xy_y_yz_1, ta2_xy_y_z_0, ta2_xy_y_z_1, ta2_xy_y_zz_0, ta2_xy_y_zz_1, ta2_xy_yy_xx_0, ta2_xy_yy_xy_0, ta2_xy_yy_xz_0, ta2_xy_yy_yy_0, ta2_xy_yy_yz_0, ta2_xy_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yy_xx_0[i] = ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + ta1_x_y_xx_1[i] + ta2_xy_y_xx_0[i] * pa_y[i] - ta2_xy_y_xx_1[i] * pc_y[i];

        ta2_xy_yy_xy_0[i] = ta2_xy_0_xy_0[i] * fe_0 - ta2_xy_0_xy_1[i] * fe_0 + ta2_xy_y_x_0[i] * fe_0 - ta2_xy_y_x_1[i] * fe_0 + ta1_x_y_xy_1[i] + ta2_xy_y_xy_0[i] * pa_y[i] - ta2_xy_y_xy_1[i] * pc_y[i];

        ta2_xy_yy_xz_0[i] = ta2_xy_0_xz_0[i] * fe_0 - ta2_xy_0_xz_1[i] * fe_0 + ta1_x_y_xz_1[i] + ta2_xy_y_xz_0[i] * pa_y[i] - ta2_xy_y_xz_1[i] * pc_y[i];

        ta2_xy_yy_yy_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + 2.0 * ta2_xy_y_y_0[i] * fe_0 - 2.0 * ta2_xy_y_y_1[i] * fe_0 + ta1_x_y_yy_1[i] + ta2_xy_y_yy_0[i] * pa_y[i] - ta2_xy_y_yy_1[i] * pc_y[i];

        ta2_xy_yy_yz_0[i] = ta2_xy_0_yz_0[i] * fe_0 - ta2_xy_0_yz_1[i] * fe_0 + ta2_xy_y_z_0[i] * fe_0 - ta2_xy_y_z_1[i] * fe_0 + ta1_x_y_yz_1[i] + ta2_xy_y_yz_0[i] * pa_y[i] - ta2_xy_y_yz_1[i] * pc_y[i];

        ta2_xy_yy_zz_0[i] = ta2_xy_0_zz_0[i] * fe_0 - ta2_xy_0_zz_1[i] * fe_0 + ta1_x_y_zz_1[i] + ta2_xy_y_zz_0[i] * pa_y[i] - ta2_xy_y_zz_1[i] * pc_y[i];
    }

    // Set up 60-66 components of targeted buffer : DD

    auto ta2_xy_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 60);

    auto ta2_xy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 61);

    auto ta2_xy_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 62);

    auto ta2_xy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 63);

    auto ta2_xy_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 64);

    auto ta2_xy_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 65);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_z_xz_1, ta1_x_z_zz_1, ta2_xy_y_xx_0, ta2_xy_y_xx_1, ta2_xy_y_xy_0, ta2_xy_y_xy_1, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_y_yz_0, ta2_xy_y_yz_1, ta2_xy_yz_xx_0, ta2_xy_yz_xy_0, ta2_xy_yz_xz_0, ta2_xy_yz_yy_0, ta2_xy_yz_yz_0, ta2_xy_yz_zz_0, ta2_xy_z_xz_0, ta2_xy_z_xz_1, ta2_xy_z_zz_0, ta2_xy_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yz_xx_0[i] = ta2_xy_y_xx_0[i] * pa_z[i] - ta2_xy_y_xx_1[i] * pc_z[i];

        ta2_xy_yz_xy_0[i] = ta2_xy_y_xy_0[i] * pa_z[i] - ta2_xy_y_xy_1[i] * pc_z[i];

        ta2_xy_yz_xz_0[i] = ta1_x_z_xz_1[i] + ta2_xy_z_xz_0[i] * pa_y[i] - ta2_xy_z_xz_1[i] * pc_y[i];

        ta2_xy_yz_yy_0[i] = ta2_xy_y_yy_0[i] * pa_z[i] - ta2_xy_y_yy_1[i] * pc_z[i];

        ta2_xy_yz_yz_0[i] = ta2_xy_y_y_0[i] * fe_0 - ta2_xy_y_y_1[i] * fe_0 + ta2_xy_y_yz_0[i] * pa_z[i] - ta2_xy_y_yz_1[i] * pc_z[i];

        ta2_xy_yz_zz_0[i] = ta1_x_z_zz_1[i] + ta2_xy_z_zz_0[i] * pa_y[i] - ta2_xy_z_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : DD

    auto ta2_xy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 66);

    auto ta2_xy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 67);

    auto ta2_xy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 68);

    auto ta2_xy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 69);

    auto ta2_xy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 70);

    auto ta2_xy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 71);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_z_x_0, ta2_xy_z_x_1, ta2_xy_z_xx_0, ta2_xy_z_xx_1, ta2_xy_z_xy_0, ta2_xy_z_xy_1, ta2_xy_z_xz_0, ta2_xy_z_xz_1, ta2_xy_z_y_0, ta2_xy_z_y_1, ta2_xy_z_yy_0, ta2_xy_z_yy_1, ta2_xy_z_yz_0, ta2_xy_z_yz_1, ta2_xy_z_z_0, ta2_xy_z_z_1, ta2_xy_z_zz_0, ta2_xy_z_zz_1, ta2_xy_zz_xx_0, ta2_xy_zz_xy_0, ta2_xy_zz_xz_0, ta2_xy_zz_yy_0, ta2_xy_zz_yz_0, ta2_xy_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zz_xx_0[i] = ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + ta2_xy_z_xx_0[i] * pa_z[i] - ta2_xy_z_xx_1[i] * pc_z[i];

        ta2_xy_zz_xy_0[i] = ta2_xy_0_xy_0[i] * fe_0 - ta2_xy_0_xy_1[i] * fe_0 + ta2_xy_z_xy_0[i] * pa_z[i] - ta2_xy_z_xy_1[i] * pc_z[i];

        ta2_xy_zz_xz_0[i] = ta2_xy_0_xz_0[i] * fe_0 - ta2_xy_0_xz_1[i] * fe_0 + ta2_xy_z_x_0[i] * fe_0 - ta2_xy_z_x_1[i] * fe_0 + ta2_xy_z_xz_0[i] * pa_z[i] - ta2_xy_z_xz_1[i] * pc_z[i];

        ta2_xy_zz_yy_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta2_xy_z_yy_0[i] * pa_z[i] - ta2_xy_z_yy_1[i] * pc_z[i];

        ta2_xy_zz_yz_0[i] = ta2_xy_0_yz_0[i] * fe_0 - ta2_xy_0_yz_1[i] * fe_0 + ta2_xy_z_y_0[i] * fe_0 - ta2_xy_z_y_1[i] * fe_0 + ta2_xy_z_yz_0[i] * pa_z[i] - ta2_xy_z_yz_1[i] * pc_z[i];

        ta2_xy_zz_zz_0[i] = ta2_xy_0_zz_0[i] * fe_0 - ta2_xy_0_zz_1[i] * fe_0 + 2.0 * ta2_xy_z_z_0[i] * fe_0 - 2.0 * ta2_xy_z_z_1[i] * fe_0 + ta2_xy_z_zz_0[i] * pa_z[i] - ta2_xy_z_zz_1[i] * pc_z[i];
    }

    // Set up 72-78 components of targeted buffer : DD

    auto ta2_xz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 72);

    auto ta2_xz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 73);

    auto ta2_xz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 74);

    auto ta2_xz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 75);

    auto ta2_xz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 76);

    auto ta2_xz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 77);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_x_xx_1, ta1_z_x_xy_1, ta1_z_x_xz_1, ta1_z_x_yy_1, ta1_z_x_yz_1, ta1_z_x_zz_1, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_x_x_0, ta2_xz_x_x_1, ta2_xz_x_xx_0, ta2_xz_x_xx_1, ta2_xz_x_xy_0, ta2_xz_x_xy_1, ta2_xz_x_xz_0, ta2_xz_x_xz_1, ta2_xz_x_y_0, ta2_xz_x_y_1, ta2_xz_x_yy_0, ta2_xz_x_yy_1, ta2_xz_x_yz_0, ta2_xz_x_yz_1, ta2_xz_x_z_0, ta2_xz_x_z_1, ta2_xz_x_zz_0, ta2_xz_x_zz_1, ta2_xz_xx_xx_0, ta2_xz_xx_xy_0, ta2_xz_xx_xz_0, ta2_xz_xx_yy_0, ta2_xz_xx_yz_0, ta2_xz_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xx_xx_0[i] = ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + 2.0 * ta2_xz_x_x_0[i] * fe_0 - 2.0 * ta2_xz_x_x_1[i] * fe_0 + ta1_z_x_xx_1[i] + ta2_xz_x_xx_0[i] * pa_x[i] - ta2_xz_x_xx_1[i] * pc_x[i];

        ta2_xz_xx_xy_0[i] = ta2_xz_0_xy_0[i] * fe_0 - ta2_xz_0_xy_1[i] * fe_0 + ta2_xz_x_y_0[i] * fe_0 - ta2_xz_x_y_1[i] * fe_0 + ta1_z_x_xy_1[i] + ta2_xz_x_xy_0[i] * pa_x[i] - ta2_xz_x_xy_1[i] * pc_x[i];

        ta2_xz_xx_xz_0[i] = ta2_xz_0_xz_0[i] * fe_0 - ta2_xz_0_xz_1[i] * fe_0 + ta2_xz_x_z_0[i] * fe_0 - ta2_xz_x_z_1[i] * fe_0 + ta1_z_x_xz_1[i] + ta2_xz_x_xz_0[i] * pa_x[i] - ta2_xz_x_xz_1[i] * pc_x[i];

        ta2_xz_xx_yy_0[i] = ta2_xz_0_yy_0[i] * fe_0 - ta2_xz_0_yy_1[i] * fe_0 + ta1_z_x_yy_1[i] + ta2_xz_x_yy_0[i] * pa_x[i] - ta2_xz_x_yy_1[i] * pc_x[i];

        ta2_xz_xx_yz_0[i] = ta2_xz_0_yz_0[i] * fe_0 - ta2_xz_0_yz_1[i] * fe_0 + ta1_z_x_yz_1[i] + ta2_xz_x_yz_0[i] * pa_x[i] - ta2_xz_x_yz_1[i] * pc_x[i];

        ta2_xz_xx_zz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta1_z_x_zz_1[i] + ta2_xz_x_zz_0[i] * pa_x[i] - ta2_xz_x_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : DD

    auto ta2_xz_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 78);

    auto ta2_xz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 79);

    auto ta2_xz_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 80);

    auto ta2_xz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 81);

    auto ta2_xz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 82);

    auto ta2_xz_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 83);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_y_yy_1, ta1_z_y_yz_1, ta2_xz_x_x_0, ta2_xz_x_x_1, ta2_xz_x_xx_0, ta2_xz_x_xx_1, ta2_xz_x_xy_0, ta2_xz_x_xy_1, ta2_xz_x_xz_0, ta2_xz_x_xz_1, ta2_xz_x_zz_0, ta2_xz_x_zz_1, ta2_xz_xy_xx_0, ta2_xz_xy_xy_0, ta2_xz_xy_xz_0, ta2_xz_xy_yy_0, ta2_xz_xy_yz_0, ta2_xz_xy_zz_0, ta2_xz_y_yy_0, ta2_xz_y_yy_1, ta2_xz_y_yz_0, ta2_xz_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xy_xx_0[i] = ta2_xz_x_xx_0[i] * pa_y[i] - ta2_xz_x_xx_1[i] * pc_y[i];

        ta2_xz_xy_xy_0[i] = ta2_xz_x_x_0[i] * fe_0 - ta2_xz_x_x_1[i] * fe_0 + ta2_xz_x_xy_0[i] * pa_y[i] - ta2_xz_x_xy_1[i] * pc_y[i];

        ta2_xz_xy_xz_0[i] = ta2_xz_x_xz_0[i] * pa_y[i] - ta2_xz_x_xz_1[i] * pc_y[i];

        ta2_xz_xy_yy_0[i] = ta1_z_y_yy_1[i] + ta2_xz_y_yy_0[i] * pa_x[i] - ta2_xz_y_yy_1[i] * pc_x[i];

        ta2_xz_xy_yz_0[i] = ta1_z_y_yz_1[i] + ta2_xz_y_yz_0[i] * pa_x[i] - ta2_xz_y_yz_1[i] * pc_x[i];

        ta2_xz_xy_zz_0[i] = ta2_xz_x_zz_0[i] * pa_y[i] - ta2_xz_x_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : DD

    auto ta2_xz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 84);

    auto ta2_xz_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 85);

    auto ta2_xz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 86);

    auto ta2_xz_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 87);

    auto ta2_xz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 88);

    auto ta2_xz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 89);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_x_xx_1, ta1_x_x_xy_1, ta1_z_z_xz_1, ta1_z_z_yy_1, ta1_z_z_yz_1, ta1_z_z_zz_1, ta2_xz_x_xx_0, ta2_xz_x_xx_1, ta2_xz_x_xy_0, ta2_xz_x_xy_1, ta2_xz_xz_xx_0, ta2_xz_xz_xy_0, ta2_xz_xz_xz_0, ta2_xz_xz_yy_0, ta2_xz_xz_yz_0, ta2_xz_xz_zz_0, ta2_xz_z_xz_0, ta2_xz_z_xz_1, ta2_xz_z_yy_0, ta2_xz_z_yy_1, ta2_xz_z_yz_0, ta2_xz_z_yz_1, ta2_xz_z_z_0, ta2_xz_z_z_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xz_xx_0[i] = ta1_x_x_xx_1[i] + ta2_xz_x_xx_0[i] * pa_z[i] - ta2_xz_x_xx_1[i] * pc_z[i];

        ta2_xz_xz_xy_0[i] = ta1_x_x_xy_1[i] + ta2_xz_x_xy_0[i] * pa_z[i] - ta2_xz_x_xy_1[i] * pc_z[i];

        ta2_xz_xz_xz_0[i] = ta2_xz_z_z_0[i] * fe_0 - ta2_xz_z_z_1[i] * fe_0 + ta1_z_z_xz_1[i] + ta2_xz_z_xz_0[i] * pa_x[i] - ta2_xz_z_xz_1[i] * pc_x[i];

        ta2_xz_xz_yy_0[i] = ta1_z_z_yy_1[i] + ta2_xz_z_yy_0[i] * pa_x[i] - ta2_xz_z_yy_1[i] * pc_x[i];

        ta2_xz_xz_yz_0[i] = ta1_z_z_yz_1[i] + ta2_xz_z_yz_0[i] * pa_x[i] - ta2_xz_z_yz_1[i] * pc_x[i];

        ta2_xz_xz_zz_0[i] = ta1_z_z_zz_1[i] + ta2_xz_z_zz_0[i] * pa_x[i] - ta2_xz_z_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : DD

    auto ta2_xz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 90);

    auto ta2_xz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 91);

    auto ta2_xz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 92);

    auto ta2_xz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 93);

    auto ta2_xz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 94);

    auto ta2_xz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 95);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_y_x_0, ta2_xz_y_x_1, ta2_xz_y_xx_0, ta2_xz_y_xx_1, ta2_xz_y_xy_0, ta2_xz_y_xy_1, ta2_xz_y_xz_0, ta2_xz_y_xz_1, ta2_xz_y_y_0, ta2_xz_y_y_1, ta2_xz_y_yy_0, ta2_xz_y_yy_1, ta2_xz_y_yz_0, ta2_xz_y_yz_1, ta2_xz_y_z_0, ta2_xz_y_z_1, ta2_xz_y_zz_0, ta2_xz_y_zz_1, ta2_xz_yy_xx_0, ta2_xz_yy_xy_0, ta2_xz_yy_xz_0, ta2_xz_yy_yy_0, ta2_xz_yy_yz_0, ta2_xz_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yy_xx_0[i] = ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + ta2_xz_y_xx_0[i] * pa_y[i] - ta2_xz_y_xx_1[i] * pc_y[i];

        ta2_xz_yy_xy_0[i] = ta2_xz_0_xy_0[i] * fe_0 - ta2_xz_0_xy_1[i] * fe_0 + ta2_xz_y_x_0[i] * fe_0 - ta2_xz_y_x_1[i] * fe_0 + ta2_xz_y_xy_0[i] * pa_y[i] - ta2_xz_y_xy_1[i] * pc_y[i];

        ta2_xz_yy_xz_0[i] = ta2_xz_0_xz_0[i] * fe_0 - ta2_xz_0_xz_1[i] * fe_0 + ta2_xz_y_xz_0[i] * pa_y[i] - ta2_xz_y_xz_1[i] * pc_y[i];

        ta2_xz_yy_yy_0[i] = ta2_xz_0_yy_0[i] * fe_0 - ta2_xz_0_yy_1[i] * fe_0 + 2.0 * ta2_xz_y_y_0[i] * fe_0 - 2.0 * ta2_xz_y_y_1[i] * fe_0 + ta2_xz_y_yy_0[i] * pa_y[i] - ta2_xz_y_yy_1[i] * pc_y[i];

        ta2_xz_yy_yz_0[i] = ta2_xz_0_yz_0[i] * fe_0 - ta2_xz_0_yz_1[i] * fe_0 + ta2_xz_y_z_0[i] * fe_0 - ta2_xz_y_z_1[i] * fe_0 + ta2_xz_y_yz_0[i] * pa_y[i] - ta2_xz_y_yz_1[i] * pc_y[i];

        ta2_xz_yy_zz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta2_xz_y_zz_0[i] * pa_y[i] - ta2_xz_y_zz_1[i] * pc_y[i];
    }

    // Set up 96-102 components of targeted buffer : DD

    auto ta2_xz_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 96);

    auto ta2_xz_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 97);

    auto ta2_xz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 98);

    auto ta2_xz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 99);

    auto ta2_xz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 100);

    auto ta2_xz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 101);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_y_xy_1, ta1_x_y_yy_1, ta2_xz_y_xy_0, ta2_xz_y_xy_1, ta2_xz_y_yy_0, ta2_xz_y_yy_1, ta2_xz_yz_xx_0, ta2_xz_yz_xy_0, ta2_xz_yz_xz_0, ta2_xz_yz_yy_0, ta2_xz_yz_yz_0, ta2_xz_yz_zz_0, ta2_xz_z_xx_0, ta2_xz_z_xx_1, ta2_xz_z_xz_0, ta2_xz_z_xz_1, ta2_xz_z_yz_0, ta2_xz_z_yz_1, ta2_xz_z_z_0, ta2_xz_z_z_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yz_xx_0[i] = ta2_xz_z_xx_0[i] * pa_y[i] - ta2_xz_z_xx_1[i] * pc_y[i];

        ta2_xz_yz_xy_0[i] = ta1_x_y_xy_1[i] + ta2_xz_y_xy_0[i] * pa_z[i] - ta2_xz_y_xy_1[i] * pc_z[i];

        ta2_xz_yz_xz_0[i] = ta2_xz_z_xz_0[i] * pa_y[i] - ta2_xz_z_xz_1[i] * pc_y[i];

        ta2_xz_yz_yy_0[i] = ta1_x_y_yy_1[i] + ta2_xz_y_yy_0[i] * pa_z[i] - ta2_xz_y_yy_1[i] * pc_z[i];

        ta2_xz_yz_yz_0[i] = ta2_xz_z_z_0[i] * fe_0 - ta2_xz_z_z_1[i] * fe_0 + ta2_xz_z_yz_0[i] * pa_y[i] - ta2_xz_z_yz_1[i] * pc_y[i];

        ta2_xz_yz_zz_0[i] = ta2_xz_z_zz_0[i] * pa_y[i] - ta2_xz_z_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : DD

    auto ta2_xz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 102);

    auto ta2_xz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 103);

    auto ta2_xz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 104);

    auto ta2_xz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 105);

    auto ta2_xz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 106);

    auto ta2_xz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 107);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_z_xx_1, ta1_x_z_xy_1, ta1_x_z_xz_1, ta1_x_z_yy_1, ta1_x_z_yz_1, ta1_x_z_zz_1, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_z_x_0, ta2_xz_z_x_1, ta2_xz_z_xx_0, ta2_xz_z_xx_1, ta2_xz_z_xy_0, ta2_xz_z_xy_1, ta2_xz_z_xz_0, ta2_xz_z_xz_1, ta2_xz_z_y_0, ta2_xz_z_y_1, ta2_xz_z_yy_0, ta2_xz_z_yy_1, ta2_xz_z_yz_0, ta2_xz_z_yz_1, ta2_xz_z_z_0, ta2_xz_z_z_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, ta2_xz_zz_xx_0, ta2_xz_zz_xy_0, ta2_xz_zz_xz_0, ta2_xz_zz_yy_0, ta2_xz_zz_yz_0, ta2_xz_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zz_xx_0[i] = ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + ta1_x_z_xx_1[i] + ta2_xz_z_xx_0[i] * pa_z[i] - ta2_xz_z_xx_1[i] * pc_z[i];

        ta2_xz_zz_xy_0[i] = ta2_xz_0_xy_0[i] * fe_0 - ta2_xz_0_xy_1[i] * fe_0 + ta1_x_z_xy_1[i] + ta2_xz_z_xy_0[i] * pa_z[i] - ta2_xz_z_xy_1[i] * pc_z[i];

        ta2_xz_zz_xz_0[i] = ta2_xz_0_xz_0[i] * fe_0 - ta2_xz_0_xz_1[i] * fe_0 + ta2_xz_z_x_0[i] * fe_0 - ta2_xz_z_x_1[i] * fe_0 + ta1_x_z_xz_1[i] + ta2_xz_z_xz_0[i] * pa_z[i] - ta2_xz_z_xz_1[i] * pc_z[i];

        ta2_xz_zz_yy_0[i] = ta2_xz_0_yy_0[i] * fe_0 - ta2_xz_0_yy_1[i] * fe_0 + ta1_x_z_yy_1[i] + ta2_xz_z_yy_0[i] * pa_z[i] - ta2_xz_z_yy_1[i] * pc_z[i];

        ta2_xz_zz_yz_0[i] = ta2_xz_0_yz_0[i] * fe_0 - ta2_xz_0_yz_1[i] * fe_0 + ta2_xz_z_y_0[i] * fe_0 - ta2_xz_z_y_1[i] * fe_0 + ta1_x_z_yz_1[i] + ta2_xz_z_yz_0[i] * pa_z[i] - ta2_xz_z_yz_1[i] * pc_z[i];

        ta2_xz_zz_zz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + 2.0 * ta2_xz_z_z_0[i] * fe_0 - 2.0 * ta2_xz_z_z_1[i] * fe_0 + ta1_x_z_zz_1[i] + ta2_xz_z_zz_0[i] * pa_z[i] - ta2_xz_z_zz_1[i] * pc_z[i];
    }

    // Set up 108-114 components of targeted buffer : DD

    auto ta2_yy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 108);

    auto ta2_yy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 109);

    auto ta2_yy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 110);

    auto ta2_yy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 111);

    auto ta2_yy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 112);

    auto ta2_yy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 113);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_x_x_0, ta2_yy_x_x_1, ta2_yy_x_xx_0, ta2_yy_x_xx_1, ta2_yy_x_xy_0, ta2_yy_x_xy_1, ta2_yy_x_xz_0, ta2_yy_x_xz_1, ta2_yy_x_y_0, ta2_yy_x_y_1, ta2_yy_x_yy_0, ta2_yy_x_yy_1, ta2_yy_x_yz_0, ta2_yy_x_yz_1, ta2_yy_x_z_0, ta2_yy_x_z_1, ta2_yy_x_zz_0, ta2_yy_x_zz_1, ta2_yy_xx_xx_0, ta2_yy_xx_xy_0, ta2_yy_xx_xz_0, ta2_yy_xx_yy_0, ta2_yy_xx_yz_0, ta2_yy_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xx_xx_0[i] = ta2_yy_0_xx_0[i] * fe_0 - ta2_yy_0_xx_1[i] * fe_0 + 2.0 * ta2_yy_x_x_0[i] * fe_0 - 2.0 * ta2_yy_x_x_1[i] * fe_0 + ta2_yy_x_xx_0[i] * pa_x[i] - ta2_yy_x_xx_1[i] * pc_x[i];

        ta2_yy_xx_xy_0[i] = ta2_yy_0_xy_0[i] * fe_0 - ta2_yy_0_xy_1[i] * fe_0 + ta2_yy_x_y_0[i] * fe_0 - ta2_yy_x_y_1[i] * fe_0 + ta2_yy_x_xy_0[i] * pa_x[i] - ta2_yy_x_xy_1[i] * pc_x[i];

        ta2_yy_xx_xz_0[i] = ta2_yy_0_xz_0[i] * fe_0 - ta2_yy_0_xz_1[i] * fe_0 + ta2_yy_x_z_0[i] * fe_0 - ta2_yy_x_z_1[i] * fe_0 + ta2_yy_x_xz_0[i] * pa_x[i] - ta2_yy_x_xz_1[i] * pc_x[i];

        ta2_yy_xx_yy_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_x_yy_0[i] * pa_x[i] - ta2_yy_x_yy_1[i] * pc_x[i];

        ta2_yy_xx_yz_0[i] = ta2_yy_0_yz_0[i] * fe_0 - ta2_yy_0_yz_1[i] * fe_0 + ta2_yy_x_yz_0[i] * pa_x[i] - ta2_yy_x_yz_1[i] * pc_x[i];

        ta2_yy_xx_zz_0[i] = ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + ta2_yy_x_zz_0[i] * pa_x[i] - ta2_yy_x_zz_1[i] * pc_x[i];
    }

    // Set up 114-120 components of targeted buffer : DD

    auto ta2_yy_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 114);

    auto ta2_yy_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 115);

    auto ta2_yy_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 116);

    auto ta2_yy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 117);

    auto ta2_yy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 118);

    auto ta2_yy_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 119);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_x_xx_1, ta1_y_x_xz_1, ta2_yy_x_xx_0, ta2_yy_x_xx_1, ta2_yy_x_xz_0, ta2_yy_x_xz_1, ta2_yy_xy_xx_0, ta2_yy_xy_xy_0, ta2_yy_xy_xz_0, ta2_yy_xy_yy_0, ta2_yy_xy_yz_0, ta2_yy_xy_zz_0, ta2_yy_y_xy_0, ta2_yy_y_xy_1, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_y_yz_0, ta2_yy_y_yz_1, ta2_yy_y_zz_0, ta2_yy_y_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xy_xx_0[i] = 2.0 * ta1_y_x_xx_1[i] + ta2_yy_x_xx_0[i] * pa_y[i] - ta2_yy_x_xx_1[i] * pc_y[i];

        ta2_yy_xy_xy_0[i] = ta2_yy_y_y_0[i] * fe_0 - ta2_yy_y_y_1[i] * fe_0 + ta2_yy_y_xy_0[i] * pa_x[i] - ta2_yy_y_xy_1[i] * pc_x[i];

        ta2_yy_xy_xz_0[i] = 2.0 * ta1_y_x_xz_1[i] + ta2_yy_x_xz_0[i] * pa_y[i] - ta2_yy_x_xz_1[i] * pc_y[i];

        ta2_yy_xy_yy_0[i] = ta2_yy_y_yy_0[i] * pa_x[i] - ta2_yy_y_yy_1[i] * pc_x[i];

        ta2_yy_xy_yz_0[i] = ta2_yy_y_yz_0[i] * pa_x[i] - ta2_yy_y_yz_1[i] * pc_x[i];

        ta2_yy_xy_zz_0[i] = ta2_yy_y_zz_0[i] * pa_x[i] - ta2_yy_y_zz_1[i] * pc_x[i];
    }

    // Set up 120-126 components of targeted buffer : DD

    auto ta2_yy_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 120);

    auto ta2_yy_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 121);

    auto ta2_yy_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 122);

    auto ta2_yy_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 123);

    auto ta2_yy_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 124);

    auto ta2_yy_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 125);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_x_xx_0, ta2_yy_x_xx_1, ta2_yy_x_xy_0, ta2_yy_x_xy_1, ta2_yy_xz_xx_0, ta2_yy_xz_xy_0, ta2_yy_xz_xz_0, ta2_yy_xz_yy_0, ta2_yy_xz_yz_0, ta2_yy_xz_zz_0, ta2_yy_z_xz_0, ta2_yy_z_xz_1, ta2_yy_z_yy_0, ta2_yy_z_yy_1, ta2_yy_z_yz_0, ta2_yy_z_yz_1, ta2_yy_z_z_0, ta2_yy_z_z_1, ta2_yy_z_zz_0, ta2_yy_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xz_xx_0[i] = ta2_yy_x_xx_0[i] * pa_z[i] - ta2_yy_x_xx_1[i] * pc_z[i];

        ta2_yy_xz_xy_0[i] = ta2_yy_x_xy_0[i] * pa_z[i] - ta2_yy_x_xy_1[i] * pc_z[i];

        ta2_yy_xz_xz_0[i] = ta2_yy_z_z_0[i] * fe_0 - ta2_yy_z_z_1[i] * fe_0 + ta2_yy_z_xz_0[i] * pa_x[i] - ta2_yy_z_xz_1[i] * pc_x[i];

        ta2_yy_xz_yy_0[i] = ta2_yy_z_yy_0[i] * pa_x[i] - ta2_yy_z_yy_1[i] * pc_x[i];

        ta2_yy_xz_yz_0[i] = ta2_yy_z_yz_0[i] * pa_x[i] - ta2_yy_z_yz_1[i] * pc_x[i];

        ta2_yy_xz_zz_0[i] = ta2_yy_z_zz_0[i] * pa_x[i] - ta2_yy_z_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : DD

    auto ta2_yy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 126);

    auto ta2_yy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 127);

    auto ta2_yy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 128);

    auto ta2_yy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 129);

    auto ta2_yy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 130);

    auto ta2_yy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 131);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_y_xx_1, ta1_y_y_xy_1, ta1_y_y_xz_1, ta1_y_y_yy_1, ta1_y_y_yz_1, ta1_y_y_zz_1, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_y_x_0, ta2_yy_y_x_1, ta2_yy_y_xx_0, ta2_yy_y_xx_1, ta2_yy_y_xy_0, ta2_yy_y_xy_1, ta2_yy_y_xz_0, ta2_yy_y_xz_1, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_y_yz_0, ta2_yy_y_yz_1, ta2_yy_y_z_0, ta2_yy_y_z_1, ta2_yy_y_zz_0, ta2_yy_y_zz_1, ta2_yy_yy_xx_0, ta2_yy_yy_xy_0, ta2_yy_yy_xz_0, ta2_yy_yy_yy_0, ta2_yy_yy_yz_0, ta2_yy_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yy_xx_0[i] = ta2_yy_0_xx_0[i] * fe_0 - ta2_yy_0_xx_1[i] * fe_0 + 2.0 * ta1_y_y_xx_1[i] + ta2_yy_y_xx_0[i] * pa_y[i] - ta2_yy_y_xx_1[i] * pc_y[i];

        ta2_yy_yy_xy_0[i] = ta2_yy_0_xy_0[i] * fe_0 - ta2_yy_0_xy_1[i] * fe_0 + ta2_yy_y_x_0[i] * fe_0 - ta2_yy_y_x_1[i] * fe_0 + 2.0 * ta1_y_y_xy_1[i] + ta2_yy_y_xy_0[i] * pa_y[i] - ta2_yy_y_xy_1[i] * pc_y[i];

        ta2_yy_yy_xz_0[i] = ta2_yy_0_xz_0[i] * fe_0 - ta2_yy_0_xz_1[i] * fe_0 + 2.0 * ta1_y_y_xz_1[i] + ta2_yy_y_xz_0[i] * pa_y[i] - ta2_yy_y_xz_1[i] * pc_y[i];

        ta2_yy_yy_yy_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + 2.0 * ta2_yy_y_y_0[i] * fe_0 - 2.0 * ta2_yy_y_y_1[i] * fe_0 + 2.0 * ta1_y_y_yy_1[i] + ta2_yy_y_yy_0[i] * pa_y[i] - ta2_yy_y_yy_1[i] * pc_y[i];

        ta2_yy_yy_yz_0[i] = ta2_yy_0_yz_0[i] * fe_0 - ta2_yy_0_yz_1[i] * fe_0 + ta2_yy_y_z_0[i] * fe_0 - ta2_yy_y_z_1[i] * fe_0 + 2.0 * ta1_y_y_yz_1[i] + ta2_yy_y_yz_0[i] * pa_y[i] - ta2_yy_y_yz_1[i] * pc_y[i];

        ta2_yy_yy_zz_0[i] = ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + 2.0 * ta1_y_y_zz_1[i] + ta2_yy_y_zz_0[i] * pa_y[i] - ta2_yy_y_zz_1[i] * pc_y[i];
    }

    // Set up 132-138 components of targeted buffer : DD

    auto ta2_yy_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 132);

    auto ta2_yy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 133);

    auto ta2_yy_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 134);

    auto ta2_yy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 135);

    auto ta2_yy_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 136);

    auto ta2_yy_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 137);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_z_xz_1, ta1_y_z_zz_1, ta2_yy_y_xx_0, ta2_yy_y_xx_1, ta2_yy_y_xy_0, ta2_yy_y_xy_1, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_y_yz_0, ta2_yy_y_yz_1, ta2_yy_yz_xx_0, ta2_yy_yz_xy_0, ta2_yy_yz_xz_0, ta2_yy_yz_yy_0, ta2_yy_yz_yz_0, ta2_yy_yz_zz_0, ta2_yy_z_xz_0, ta2_yy_z_xz_1, ta2_yy_z_zz_0, ta2_yy_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yz_xx_0[i] = ta2_yy_y_xx_0[i] * pa_z[i] - ta2_yy_y_xx_1[i] * pc_z[i];

        ta2_yy_yz_xy_0[i] = ta2_yy_y_xy_0[i] * pa_z[i] - ta2_yy_y_xy_1[i] * pc_z[i];

        ta2_yy_yz_xz_0[i] = 2.0 * ta1_y_z_xz_1[i] + ta2_yy_z_xz_0[i] * pa_y[i] - ta2_yy_z_xz_1[i] * pc_y[i];

        ta2_yy_yz_yy_0[i] = ta2_yy_y_yy_0[i] * pa_z[i] - ta2_yy_y_yy_1[i] * pc_z[i];

        ta2_yy_yz_yz_0[i] = ta2_yy_y_y_0[i] * fe_0 - ta2_yy_y_y_1[i] * fe_0 + ta2_yy_y_yz_0[i] * pa_z[i] - ta2_yy_y_yz_1[i] * pc_z[i];

        ta2_yy_yz_zz_0[i] = 2.0 * ta1_y_z_zz_1[i] + ta2_yy_z_zz_0[i] * pa_y[i] - ta2_yy_z_zz_1[i] * pc_y[i];
    }

    // Set up 138-144 components of targeted buffer : DD

    auto ta2_yy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 138);

    auto ta2_yy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 139);

    auto ta2_yy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 140);

    auto ta2_yy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 141);

    auto ta2_yy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 142);

    auto ta2_yy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 143);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_z_x_0, ta2_yy_z_x_1, ta2_yy_z_xx_0, ta2_yy_z_xx_1, ta2_yy_z_xy_0, ta2_yy_z_xy_1, ta2_yy_z_xz_0, ta2_yy_z_xz_1, ta2_yy_z_y_0, ta2_yy_z_y_1, ta2_yy_z_yy_0, ta2_yy_z_yy_1, ta2_yy_z_yz_0, ta2_yy_z_yz_1, ta2_yy_z_z_0, ta2_yy_z_z_1, ta2_yy_z_zz_0, ta2_yy_z_zz_1, ta2_yy_zz_xx_0, ta2_yy_zz_xy_0, ta2_yy_zz_xz_0, ta2_yy_zz_yy_0, ta2_yy_zz_yz_0, ta2_yy_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zz_xx_0[i] = ta2_yy_0_xx_0[i] * fe_0 - ta2_yy_0_xx_1[i] * fe_0 + ta2_yy_z_xx_0[i] * pa_z[i] - ta2_yy_z_xx_1[i] * pc_z[i];

        ta2_yy_zz_xy_0[i] = ta2_yy_0_xy_0[i] * fe_0 - ta2_yy_0_xy_1[i] * fe_0 + ta2_yy_z_xy_0[i] * pa_z[i] - ta2_yy_z_xy_1[i] * pc_z[i];

        ta2_yy_zz_xz_0[i] = ta2_yy_0_xz_0[i] * fe_0 - ta2_yy_0_xz_1[i] * fe_0 + ta2_yy_z_x_0[i] * fe_0 - ta2_yy_z_x_1[i] * fe_0 + ta2_yy_z_xz_0[i] * pa_z[i] - ta2_yy_z_xz_1[i] * pc_z[i];

        ta2_yy_zz_yy_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_z_yy_0[i] * pa_z[i] - ta2_yy_z_yy_1[i] * pc_z[i];

        ta2_yy_zz_yz_0[i] = ta2_yy_0_yz_0[i] * fe_0 - ta2_yy_0_yz_1[i] * fe_0 + ta2_yy_z_y_0[i] * fe_0 - ta2_yy_z_y_1[i] * fe_0 + ta2_yy_z_yz_0[i] * pa_z[i] - ta2_yy_z_yz_1[i] * pc_z[i];

        ta2_yy_zz_zz_0[i] = ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + 2.0 * ta2_yy_z_z_0[i] * fe_0 - 2.0 * ta2_yy_z_z_1[i] * fe_0 + ta2_yy_z_zz_0[i] * pa_z[i] - ta2_yy_z_zz_1[i] * pc_z[i];
    }

    // Set up 144-150 components of targeted buffer : DD

    auto ta2_yz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 144);

    auto ta2_yz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 145);

    auto ta2_yz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 146);

    auto ta2_yz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 147);

    auto ta2_yz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 148);

    auto ta2_yz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 149);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_x_x_0, ta2_yz_x_x_1, ta2_yz_x_xx_0, ta2_yz_x_xx_1, ta2_yz_x_xy_0, ta2_yz_x_xy_1, ta2_yz_x_xz_0, ta2_yz_x_xz_1, ta2_yz_x_y_0, ta2_yz_x_y_1, ta2_yz_x_yy_0, ta2_yz_x_yy_1, ta2_yz_x_yz_0, ta2_yz_x_yz_1, ta2_yz_x_z_0, ta2_yz_x_z_1, ta2_yz_x_zz_0, ta2_yz_x_zz_1, ta2_yz_xx_xx_0, ta2_yz_xx_xy_0, ta2_yz_xx_xz_0, ta2_yz_xx_yy_0, ta2_yz_xx_yz_0, ta2_yz_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xx_xx_0[i] = ta2_yz_0_xx_0[i] * fe_0 - ta2_yz_0_xx_1[i] * fe_0 + 2.0 * ta2_yz_x_x_0[i] * fe_0 - 2.0 * ta2_yz_x_x_1[i] * fe_0 + ta2_yz_x_xx_0[i] * pa_x[i] - ta2_yz_x_xx_1[i] * pc_x[i];

        ta2_yz_xx_xy_0[i] = ta2_yz_0_xy_0[i] * fe_0 - ta2_yz_0_xy_1[i] * fe_0 + ta2_yz_x_y_0[i] * fe_0 - ta2_yz_x_y_1[i] * fe_0 + ta2_yz_x_xy_0[i] * pa_x[i] - ta2_yz_x_xy_1[i] * pc_x[i];

        ta2_yz_xx_xz_0[i] = ta2_yz_0_xz_0[i] * fe_0 - ta2_yz_0_xz_1[i] * fe_0 + ta2_yz_x_z_0[i] * fe_0 - ta2_yz_x_z_1[i] * fe_0 + ta2_yz_x_xz_0[i] * pa_x[i] - ta2_yz_x_xz_1[i] * pc_x[i];

        ta2_yz_xx_yy_0[i] = ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + ta2_yz_x_yy_0[i] * pa_x[i] - ta2_yz_x_yy_1[i] * pc_x[i];

        ta2_yz_xx_yz_0[i] = ta2_yz_0_yz_0[i] * fe_0 - ta2_yz_0_yz_1[i] * fe_0 + ta2_yz_x_yz_0[i] * pa_x[i] - ta2_yz_x_yz_1[i] * pc_x[i];

        ta2_yz_xx_zz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta2_yz_x_zz_0[i] * pa_x[i] - ta2_yz_x_zz_1[i] * pc_x[i];
    }

    // Set up 150-156 components of targeted buffer : DD

    auto ta2_yz_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 150);

    auto ta2_yz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 151);

    auto ta2_yz_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 152);

    auto ta2_yz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 153);

    auto ta2_yz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 154);

    auto ta2_yz_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 155);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_x_xx_1, ta1_z_x_xz_1, ta2_yz_x_xx_0, ta2_yz_x_xx_1, ta2_yz_x_xz_0, ta2_yz_x_xz_1, ta2_yz_xy_xx_0, ta2_yz_xy_xy_0, ta2_yz_xy_xz_0, ta2_yz_xy_yy_0, ta2_yz_xy_yz_0, ta2_yz_xy_zz_0, ta2_yz_y_xy_0, ta2_yz_y_xy_1, ta2_yz_y_y_0, ta2_yz_y_y_1, ta2_yz_y_yy_0, ta2_yz_y_yy_1, ta2_yz_y_yz_0, ta2_yz_y_yz_1, ta2_yz_y_zz_0, ta2_yz_y_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xy_xx_0[i] = ta1_z_x_xx_1[i] + ta2_yz_x_xx_0[i] * pa_y[i] - ta2_yz_x_xx_1[i] * pc_y[i];

        ta2_yz_xy_xy_0[i] = ta2_yz_y_y_0[i] * fe_0 - ta2_yz_y_y_1[i] * fe_0 + ta2_yz_y_xy_0[i] * pa_x[i] - ta2_yz_y_xy_1[i] * pc_x[i];

        ta2_yz_xy_xz_0[i] = ta1_z_x_xz_1[i] + ta2_yz_x_xz_0[i] * pa_y[i] - ta2_yz_x_xz_1[i] * pc_y[i];

        ta2_yz_xy_yy_0[i] = ta2_yz_y_yy_0[i] * pa_x[i] - ta2_yz_y_yy_1[i] * pc_x[i];

        ta2_yz_xy_yz_0[i] = ta2_yz_y_yz_0[i] * pa_x[i] - ta2_yz_y_yz_1[i] * pc_x[i];

        ta2_yz_xy_zz_0[i] = ta2_yz_y_zz_0[i] * pa_x[i] - ta2_yz_y_zz_1[i] * pc_x[i];
    }

    // Set up 156-162 components of targeted buffer : DD

    auto ta2_yz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 156);

    auto ta2_yz_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 157);

    auto ta2_yz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 158);

    auto ta2_yz_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 159);

    auto ta2_yz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 160);

    auto ta2_yz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 161);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_x_xx_1, ta1_y_x_xy_1, ta2_yz_x_xx_0, ta2_yz_x_xx_1, ta2_yz_x_xy_0, ta2_yz_x_xy_1, ta2_yz_xz_xx_0, ta2_yz_xz_xy_0, ta2_yz_xz_xz_0, ta2_yz_xz_yy_0, ta2_yz_xz_yz_0, ta2_yz_xz_zz_0, ta2_yz_z_xz_0, ta2_yz_z_xz_1, ta2_yz_z_yy_0, ta2_yz_z_yy_1, ta2_yz_z_yz_0, ta2_yz_z_yz_1, ta2_yz_z_z_0, ta2_yz_z_z_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xz_xx_0[i] = ta1_y_x_xx_1[i] + ta2_yz_x_xx_0[i] * pa_z[i] - ta2_yz_x_xx_1[i] * pc_z[i];

        ta2_yz_xz_xy_0[i] = ta1_y_x_xy_1[i] + ta2_yz_x_xy_0[i] * pa_z[i] - ta2_yz_x_xy_1[i] * pc_z[i];

        ta2_yz_xz_xz_0[i] = ta2_yz_z_z_0[i] * fe_0 - ta2_yz_z_z_1[i] * fe_0 + ta2_yz_z_xz_0[i] * pa_x[i] - ta2_yz_z_xz_1[i] * pc_x[i];

        ta2_yz_xz_yy_0[i] = ta2_yz_z_yy_0[i] * pa_x[i] - ta2_yz_z_yy_1[i] * pc_x[i];

        ta2_yz_xz_yz_0[i] = ta2_yz_z_yz_0[i] * pa_x[i] - ta2_yz_z_yz_1[i] * pc_x[i];

        ta2_yz_xz_zz_0[i] = ta2_yz_z_zz_0[i] * pa_x[i] - ta2_yz_z_zz_1[i] * pc_x[i];
    }

    // Set up 162-168 components of targeted buffer : DD

    auto ta2_yz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 162);

    auto ta2_yz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 163);

    auto ta2_yz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 164);

    auto ta2_yz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 165);

    auto ta2_yz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 166);

    auto ta2_yz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 167);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_y_xx_1, ta1_z_y_xy_1, ta1_z_y_xz_1, ta1_z_y_yy_1, ta1_z_y_yz_1, ta1_z_y_zz_1, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_y_x_0, ta2_yz_y_x_1, ta2_yz_y_xx_0, ta2_yz_y_xx_1, ta2_yz_y_xy_0, ta2_yz_y_xy_1, ta2_yz_y_xz_0, ta2_yz_y_xz_1, ta2_yz_y_y_0, ta2_yz_y_y_1, ta2_yz_y_yy_0, ta2_yz_y_yy_1, ta2_yz_y_yz_0, ta2_yz_y_yz_1, ta2_yz_y_z_0, ta2_yz_y_z_1, ta2_yz_y_zz_0, ta2_yz_y_zz_1, ta2_yz_yy_xx_0, ta2_yz_yy_xy_0, ta2_yz_yy_xz_0, ta2_yz_yy_yy_0, ta2_yz_yy_yz_0, ta2_yz_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yy_xx_0[i] = ta2_yz_0_xx_0[i] * fe_0 - ta2_yz_0_xx_1[i] * fe_0 + ta1_z_y_xx_1[i] + ta2_yz_y_xx_0[i] * pa_y[i] - ta2_yz_y_xx_1[i] * pc_y[i];

        ta2_yz_yy_xy_0[i] = ta2_yz_0_xy_0[i] * fe_0 - ta2_yz_0_xy_1[i] * fe_0 + ta2_yz_y_x_0[i] * fe_0 - ta2_yz_y_x_1[i] * fe_0 + ta1_z_y_xy_1[i] + ta2_yz_y_xy_0[i] * pa_y[i] - ta2_yz_y_xy_1[i] * pc_y[i];

        ta2_yz_yy_xz_0[i] = ta2_yz_0_xz_0[i] * fe_0 - ta2_yz_0_xz_1[i] * fe_0 + ta1_z_y_xz_1[i] + ta2_yz_y_xz_0[i] * pa_y[i] - ta2_yz_y_xz_1[i] * pc_y[i];

        ta2_yz_yy_yy_0[i] = ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + 2.0 * ta2_yz_y_y_0[i] * fe_0 - 2.0 * ta2_yz_y_y_1[i] * fe_0 + ta1_z_y_yy_1[i] + ta2_yz_y_yy_0[i] * pa_y[i] - ta2_yz_y_yy_1[i] * pc_y[i];

        ta2_yz_yy_yz_0[i] = ta2_yz_0_yz_0[i] * fe_0 - ta2_yz_0_yz_1[i] * fe_0 + ta2_yz_y_z_0[i] * fe_0 - ta2_yz_y_z_1[i] * fe_0 + ta1_z_y_yz_1[i] + ta2_yz_y_yz_0[i] * pa_y[i] - ta2_yz_y_yz_1[i] * pc_y[i];

        ta2_yz_yy_zz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta1_z_y_zz_1[i] + ta2_yz_y_zz_0[i] * pa_y[i] - ta2_yz_y_zz_1[i] * pc_y[i];
    }

    // Set up 168-174 components of targeted buffer : DD

    auto ta2_yz_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 168);

    auto ta2_yz_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 169);

    auto ta2_yz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 170);

    auto ta2_yz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 171);

    auto ta2_yz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 172);

    auto ta2_yz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 173);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_y_xy_1, ta1_y_y_yy_1, ta1_z_z_xx_1, ta1_z_z_xz_1, ta1_z_z_yz_1, ta1_z_z_zz_1, ta2_yz_y_xy_0, ta2_yz_y_xy_1, ta2_yz_y_yy_0, ta2_yz_y_yy_1, ta2_yz_yz_xx_0, ta2_yz_yz_xy_0, ta2_yz_yz_xz_0, ta2_yz_yz_yy_0, ta2_yz_yz_yz_0, ta2_yz_yz_zz_0, ta2_yz_z_xx_0, ta2_yz_z_xx_1, ta2_yz_z_xz_0, ta2_yz_z_xz_1, ta2_yz_z_yz_0, ta2_yz_z_yz_1, ta2_yz_z_z_0, ta2_yz_z_z_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yz_xx_0[i] = ta1_z_z_xx_1[i] + ta2_yz_z_xx_0[i] * pa_y[i] - ta2_yz_z_xx_1[i] * pc_y[i];

        ta2_yz_yz_xy_0[i] = ta1_y_y_xy_1[i] + ta2_yz_y_xy_0[i] * pa_z[i] - ta2_yz_y_xy_1[i] * pc_z[i];

        ta2_yz_yz_xz_0[i] = ta1_z_z_xz_1[i] + ta2_yz_z_xz_0[i] * pa_y[i] - ta2_yz_z_xz_1[i] * pc_y[i];

        ta2_yz_yz_yy_0[i] = ta1_y_y_yy_1[i] + ta2_yz_y_yy_0[i] * pa_z[i] - ta2_yz_y_yy_1[i] * pc_z[i];

        ta2_yz_yz_yz_0[i] = ta2_yz_z_z_0[i] * fe_0 - ta2_yz_z_z_1[i] * fe_0 + ta1_z_z_yz_1[i] + ta2_yz_z_yz_0[i] * pa_y[i] - ta2_yz_z_yz_1[i] * pc_y[i];

        ta2_yz_yz_zz_0[i] = ta1_z_z_zz_1[i] + ta2_yz_z_zz_0[i] * pa_y[i] - ta2_yz_z_zz_1[i] * pc_y[i];
    }

    // Set up 174-180 components of targeted buffer : DD

    auto ta2_yz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 174);

    auto ta2_yz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 175);

    auto ta2_yz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 176);

    auto ta2_yz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 177);

    auto ta2_yz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 178);

    auto ta2_yz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 179);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_z_xx_1, ta1_y_z_xy_1, ta1_y_z_xz_1, ta1_y_z_yy_1, ta1_y_z_yz_1, ta1_y_z_zz_1, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_z_x_0, ta2_yz_z_x_1, ta2_yz_z_xx_0, ta2_yz_z_xx_1, ta2_yz_z_xy_0, ta2_yz_z_xy_1, ta2_yz_z_xz_0, ta2_yz_z_xz_1, ta2_yz_z_y_0, ta2_yz_z_y_1, ta2_yz_z_yy_0, ta2_yz_z_yy_1, ta2_yz_z_yz_0, ta2_yz_z_yz_1, ta2_yz_z_z_0, ta2_yz_z_z_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, ta2_yz_zz_xx_0, ta2_yz_zz_xy_0, ta2_yz_zz_xz_0, ta2_yz_zz_yy_0, ta2_yz_zz_yz_0, ta2_yz_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zz_xx_0[i] = ta2_yz_0_xx_0[i] * fe_0 - ta2_yz_0_xx_1[i] * fe_0 + ta1_y_z_xx_1[i] + ta2_yz_z_xx_0[i] * pa_z[i] - ta2_yz_z_xx_1[i] * pc_z[i];

        ta2_yz_zz_xy_0[i] = ta2_yz_0_xy_0[i] * fe_0 - ta2_yz_0_xy_1[i] * fe_0 + ta1_y_z_xy_1[i] + ta2_yz_z_xy_0[i] * pa_z[i] - ta2_yz_z_xy_1[i] * pc_z[i];

        ta2_yz_zz_xz_0[i] = ta2_yz_0_xz_0[i] * fe_0 - ta2_yz_0_xz_1[i] * fe_0 + ta2_yz_z_x_0[i] * fe_0 - ta2_yz_z_x_1[i] * fe_0 + ta1_y_z_xz_1[i] + ta2_yz_z_xz_0[i] * pa_z[i] - ta2_yz_z_xz_1[i] * pc_z[i];

        ta2_yz_zz_yy_0[i] = ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + ta1_y_z_yy_1[i] + ta2_yz_z_yy_0[i] * pa_z[i] - ta2_yz_z_yy_1[i] * pc_z[i];

        ta2_yz_zz_yz_0[i] = ta2_yz_0_yz_0[i] * fe_0 - ta2_yz_0_yz_1[i] * fe_0 + ta2_yz_z_y_0[i] * fe_0 - ta2_yz_z_y_1[i] * fe_0 + ta1_y_z_yz_1[i] + ta2_yz_z_yz_0[i] * pa_z[i] - ta2_yz_z_yz_1[i] * pc_z[i];

        ta2_yz_zz_zz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + 2.0 * ta2_yz_z_z_0[i] * fe_0 - 2.0 * ta2_yz_z_z_1[i] * fe_0 + ta1_y_z_zz_1[i] + ta2_yz_z_zz_0[i] * pa_z[i] - ta2_yz_z_zz_1[i] * pc_z[i];
    }

    // Set up 180-186 components of targeted buffer : DD

    auto ta2_zz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 180);

    auto ta2_zz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 181);

    auto ta2_zz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 182);

    auto ta2_zz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 183);

    auto ta2_zz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 184);

    auto ta2_zz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 185);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_x_x_0, ta2_zz_x_x_1, ta2_zz_x_xx_0, ta2_zz_x_xx_1, ta2_zz_x_xy_0, ta2_zz_x_xy_1, ta2_zz_x_xz_0, ta2_zz_x_xz_1, ta2_zz_x_y_0, ta2_zz_x_y_1, ta2_zz_x_yy_0, ta2_zz_x_yy_1, ta2_zz_x_yz_0, ta2_zz_x_yz_1, ta2_zz_x_z_0, ta2_zz_x_z_1, ta2_zz_x_zz_0, ta2_zz_x_zz_1, ta2_zz_xx_xx_0, ta2_zz_xx_xy_0, ta2_zz_xx_xz_0, ta2_zz_xx_yy_0, ta2_zz_xx_yz_0, ta2_zz_xx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xx_xx_0[i] = ta2_zz_0_xx_0[i] * fe_0 - ta2_zz_0_xx_1[i] * fe_0 + 2.0 * ta2_zz_x_x_0[i] * fe_0 - 2.0 * ta2_zz_x_x_1[i] * fe_0 + ta2_zz_x_xx_0[i] * pa_x[i] - ta2_zz_x_xx_1[i] * pc_x[i];

        ta2_zz_xx_xy_0[i] = ta2_zz_0_xy_0[i] * fe_0 - ta2_zz_0_xy_1[i] * fe_0 + ta2_zz_x_y_0[i] * fe_0 - ta2_zz_x_y_1[i] * fe_0 + ta2_zz_x_xy_0[i] * pa_x[i] - ta2_zz_x_xy_1[i] * pc_x[i];

        ta2_zz_xx_xz_0[i] = ta2_zz_0_xz_0[i] * fe_0 - ta2_zz_0_xz_1[i] * fe_0 + ta2_zz_x_z_0[i] * fe_0 - ta2_zz_x_z_1[i] * fe_0 + ta2_zz_x_xz_0[i] * pa_x[i] - ta2_zz_x_xz_1[i] * pc_x[i];

        ta2_zz_xx_yy_0[i] = ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + ta2_zz_x_yy_0[i] * pa_x[i] - ta2_zz_x_yy_1[i] * pc_x[i];

        ta2_zz_xx_yz_0[i] = ta2_zz_0_yz_0[i] * fe_0 - ta2_zz_0_yz_1[i] * fe_0 + ta2_zz_x_yz_0[i] * pa_x[i] - ta2_zz_x_yz_1[i] * pc_x[i];

        ta2_zz_xx_zz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_x_zz_0[i] * pa_x[i] - ta2_zz_x_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : DD

    auto ta2_zz_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 186);

    auto ta2_zz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 187);

    auto ta2_zz_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 188);

    auto ta2_zz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 189);

    auto ta2_zz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 190);

    auto ta2_zz_xy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 191);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_x_xx_0, ta2_zz_x_xx_1, ta2_zz_x_xz_0, ta2_zz_x_xz_1, ta2_zz_xy_xx_0, ta2_zz_xy_xy_0, ta2_zz_xy_xz_0, ta2_zz_xy_yy_0, ta2_zz_xy_yz_0, ta2_zz_xy_zz_0, ta2_zz_y_xy_0, ta2_zz_y_xy_1, ta2_zz_y_y_0, ta2_zz_y_y_1, ta2_zz_y_yy_0, ta2_zz_y_yy_1, ta2_zz_y_yz_0, ta2_zz_y_yz_1, ta2_zz_y_zz_0, ta2_zz_y_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xy_xx_0[i] = ta2_zz_x_xx_0[i] * pa_y[i] - ta2_zz_x_xx_1[i] * pc_y[i];

        ta2_zz_xy_xy_0[i] = ta2_zz_y_y_0[i] * fe_0 - ta2_zz_y_y_1[i] * fe_0 + ta2_zz_y_xy_0[i] * pa_x[i] - ta2_zz_y_xy_1[i] * pc_x[i];

        ta2_zz_xy_xz_0[i] = ta2_zz_x_xz_0[i] * pa_y[i] - ta2_zz_x_xz_1[i] * pc_y[i];

        ta2_zz_xy_yy_0[i] = ta2_zz_y_yy_0[i] * pa_x[i] - ta2_zz_y_yy_1[i] * pc_x[i];

        ta2_zz_xy_yz_0[i] = ta2_zz_y_yz_0[i] * pa_x[i] - ta2_zz_y_yz_1[i] * pc_x[i];

        ta2_zz_xy_zz_0[i] = ta2_zz_y_zz_0[i] * pa_x[i] - ta2_zz_y_zz_1[i] * pc_x[i];
    }

    // Set up 192-198 components of targeted buffer : DD

    auto ta2_zz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 192);

    auto ta2_zz_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 193);

    auto ta2_zz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 194);

    auto ta2_zz_xz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 195);

    auto ta2_zz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 196);

    auto ta2_zz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 197);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_x_xx_1, ta1_z_x_xy_1, ta2_zz_x_xx_0, ta2_zz_x_xx_1, ta2_zz_x_xy_0, ta2_zz_x_xy_1, ta2_zz_xz_xx_0, ta2_zz_xz_xy_0, ta2_zz_xz_xz_0, ta2_zz_xz_yy_0, ta2_zz_xz_yz_0, ta2_zz_xz_zz_0, ta2_zz_z_xz_0, ta2_zz_z_xz_1, ta2_zz_z_yy_0, ta2_zz_z_yy_1, ta2_zz_z_yz_0, ta2_zz_z_yz_1, ta2_zz_z_z_0, ta2_zz_z_z_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xz_xx_0[i] = 2.0 * ta1_z_x_xx_1[i] + ta2_zz_x_xx_0[i] * pa_z[i] - ta2_zz_x_xx_1[i] * pc_z[i];

        ta2_zz_xz_xy_0[i] = 2.0 * ta1_z_x_xy_1[i] + ta2_zz_x_xy_0[i] * pa_z[i] - ta2_zz_x_xy_1[i] * pc_z[i];

        ta2_zz_xz_xz_0[i] = ta2_zz_z_z_0[i] * fe_0 - ta2_zz_z_z_1[i] * fe_0 + ta2_zz_z_xz_0[i] * pa_x[i] - ta2_zz_z_xz_1[i] * pc_x[i];

        ta2_zz_xz_yy_0[i] = ta2_zz_z_yy_0[i] * pa_x[i] - ta2_zz_z_yy_1[i] * pc_x[i];

        ta2_zz_xz_yz_0[i] = ta2_zz_z_yz_0[i] * pa_x[i] - ta2_zz_z_yz_1[i] * pc_x[i];

        ta2_zz_xz_zz_0[i] = ta2_zz_z_zz_0[i] * pa_x[i] - ta2_zz_z_zz_1[i] * pc_x[i];
    }

    // Set up 198-204 components of targeted buffer : DD

    auto ta2_zz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 198);

    auto ta2_zz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 199);

    auto ta2_zz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 200);

    auto ta2_zz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 201);

    auto ta2_zz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 202);

    auto ta2_zz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 203);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_y_x_0, ta2_zz_y_x_1, ta2_zz_y_xx_0, ta2_zz_y_xx_1, ta2_zz_y_xy_0, ta2_zz_y_xy_1, ta2_zz_y_xz_0, ta2_zz_y_xz_1, ta2_zz_y_y_0, ta2_zz_y_y_1, ta2_zz_y_yy_0, ta2_zz_y_yy_1, ta2_zz_y_yz_0, ta2_zz_y_yz_1, ta2_zz_y_z_0, ta2_zz_y_z_1, ta2_zz_y_zz_0, ta2_zz_y_zz_1, ta2_zz_yy_xx_0, ta2_zz_yy_xy_0, ta2_zz_yy_xz_0, ta2_zz_yy_yy_0, ta2_zz_yy_yz_0, ta2_zz_yy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yy_xx_0[i] = ta2_zz_0_xx_0[i] * fe_0 - ta2_zz_0_xx_1[i] * fe_0 + ta2_zz_y_xx_0[i] * pa_y[i] - ta2_zz_y_xx_1[i] * pc_y[i];

        ta2_zz_yy_xy_0[i] = ta2_zz_0_xy_0[i] * fe_0 - ta2_zz_0_xy_1[i] * fe_0 + ta2_zz_y_x_0[i] * fe_0 - ta2_zz_y_x_1[i] * fe_0 + ta2_zz_y_xy_0[i] * pa_y[i] - ta2_zz_y_xy_1[i] * pc_y[i];

        ta2_zz_yy_xz_0[i] = ta2_zz_0_xz_0[i] * fe_0 - ta2_zz_0_xz_1[i] * fe_0 + ta2_zz_y_xz_0[i] * pa_y[i] - ta2_zz_y_xz_1[i] * pc_y[i];

        ta2_zz_yy_yy_0[i] = ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + 2.0 * ta2_zz_y_y_0[i] * fe_0 - 2.0 * ta2_zz_y_y_1[i] * fe_0 + ta2_zz_y_yy_0[i] * pa_y[i] - ta2_zz_y_yy_1[i] * pc_y[i];

        ta2_zz_yy_yz_0[i] = ta2_zz_0_yz_0[i] * fe_0 - ta2_zz_0_yz_1[i] * fe_0 + ta2_zz_y_z_0[i] * fe_0 - ta2_zz_y_z_1[i] * fe_0 + ta2_zz_y_yz_0[i] * pa_y[i] - ta2_zz_y_yz_1[i] * pc_y[i];

        ta2_zz_yy_zz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_y_zz_0[i] * pa_y[i] - ta2_zz_y_zz_1[i] * pc_y[i];
    }

    // Set up 204-210 components of targeted buffer : DD

    auto ta2_zz_yz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 204);

    auto ta2_zz_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 205);

    auto ta2_zz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 206);

    auto ta2_zz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 207);

    auto ta2_zz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 208);

    auto ta2_zz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 209);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_y_xy_1, ta1_z_y_yy_1, ta2_zz_y_xy_0, ta2_zz_y_xy_1, ta2_zz_y_yy_0, ta2_zz_y_yy_1, ta2_zz_yz_xx_0, ta2_zz_yz_xy_0, ta2_zz_yz_xz_0, ta2_zz_yz_yy_0, ta2_zz_yz_yz_0, ta2_zz_yz_zz_0, ta2_zz_z_xx_0, ta2_zz_z_xx_1, ta2_zz_z_xz_0, ta2_zz_z_xz_1, ta2_zz_z_yz_0, ta2_zz_z_yz_1, ta2_zz_z_z_0, ta2_zz_z_z_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yz_xx_0[i] = ta2_zz_z_xx_0[i] * pa_y[i] - ta2_zz_z_xx_1[i] * pc_y[i];

        ta2_zz_yz_xy_0[i] = 2.0 * ta1_z_y_xy_1[i] + ta2_zz_y_xy_0[i] * pa_z[i] - ta2_zz_y_xy_1[i] * pc_z[i];

        ta2_zz_yz_xz_0[i] = ta2_zz_z_xz_0[i] * pa_y[i] - ta2_zz_z_xz_1[i] * pc_y[i];

        ta2_zz_yz_yy_0[i] = 2.0 * ta1_z_y_yy_1[i] + ta2_zz_y_yy_0[i] * pa_z[i] - ta2_zz_y_yy_1[i] * pc_z[i];

        ta2_zz_yz_yz_0[i] = ta2_zz_z_z_0[i] * fe_0 - ta2_zz_z_z_1[i] * fe_0 + ta2_zz_z_yz_0[i] * pa_y[i] - ta2_zz_z_yz_1[i] * pc_y[i];

        ta2_zz_yz_zz_0[i] = ta2_zz_z_zz_0[i] * pa_y[i] - ta2_zz_z_zz_1[i] * pc_y[i];
    }

    // Set up 210-216 components of targeted buffer : DD

    auto ta2_zz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 210);

    auto ta2_zz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 211);

    auto ta2_zz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 212);

    auto ta2_zz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 213);

    auto ta2_zz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 214);

    auto ta2_zz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 215);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_z_xx_1, ta1_z_z_xy_1, ta1_z_z_xz_1, ta1_z_z_yy_1, ta1_z_z_yz_1, ta1_z_z_zz_1, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_z_x_0, ta2_zz_z_x_1, ta2_zz_z_xx_0, ta2_zz_z_xx_1, ta2_zz_z_xy_0, ta2_zz_z_xy_1, ta2_zz_z_xz_0, ta2_zz_z_xz_1, ta2_zz_z_y_0, ta2_zz_z_y_1, ta2_zz_z_yy_0, ta2_zz_z_yy_1, ta2_zz_z_yz_0, ta2_zz_z_yz_1, ta2_zz_z_z_0, ta2_zz_z_z_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, ta2_zz_zz_xx_0, ta2_zz_zz_xy_0, ta2_zz_zz_xz_0, ta2_zz_zz_yy_0, ta2_zz_zz_yz_0, ta2_zz_zz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zz_xx_0[i] = ta2_zz_0_xx_0[i] * fe_0 - ta2_zz_0_xx_1[i] * fe_0 + 2.0 * ta1_z_z_xx_1[i] + ta2_zz_z_xx_0[i] * pa_z[i] - ta2_zz_z_xx_1[i] * pc_z[i];

        ta2_zz_zz_xy_0[i] = ta2_zz_0_xy_0[i] * fe_0 - ta2_zz_0_xy_1[i] * fe_0 + 2.0 * ta1_z_z_xy_1[i] + ta2_zz_z_xy_0[i] * pa_z[i] - ta2_zz_z_xy_1[i] * pc_z[i];

        ta2_zz_zz_xz_0[i] = ta2_zz_0_xz_0[i] * fe_0 - ta2_zz_0_xz_1[i] * fe_0 + ta2_zz_z_x_0[i] * fe_0 - ta2_zz_z_x_1[i] * fe_0 + 2.0 * ta1_z_z_xz_1[i] + ta2_zz_z_xz_0[i] * pa_z[i] - ta2_zz_z_xz_1[i] * pc_z[i];

        ta2_zz_zz_yy_0[i] = ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + 2.0 * ta1_z_z_yy_1[i] + ta2_zz_z_yy_0[i] * pa_z[i] - ta2_zz_z_yy_1[i] * pc_z[i];

        ta2_zz_zz_yz_0[i] = ta2_zz_0_yz_0[i] * fe_0 - ta2_zz_0_yz_1[i] * fe_0 + ta2_zz_z_y_0[i] * fe_0 - ta2_zz_z_y_1[i] * fe_0 + 2.0 * ta1_z_z_yz_1[i] + ta2_zz_z_yz_0[i] * pa_z[i] - ta2_zz_z_yz_1[i] * pc_z[i];

        ta2_zz_zz_zz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + 2.0 * ta2_zz_z_z_0[i] * fe_0 - 2.0 * ta2_zz_z_z_1[i] * fe_0 + 2.0 * ta1_z_z_zz_1[i] + ta2_zz_z_zz_0[i] * pa_z[i] - ta2_zz_z_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

