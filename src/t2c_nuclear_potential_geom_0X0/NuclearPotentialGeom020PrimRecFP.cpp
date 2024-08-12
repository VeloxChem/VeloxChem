#include "NuclearPotentialGeom020PrimRecFP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_fp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_fp,
                                        const size_t              idx_npot_geom_020_0_pp,
                                        const size_t              idx_npot_geom_020_1_pp,
                                        const size_t              idx_npot_geom_020_0_ds,
                                        const size_t              idx_npot_geom_020_1_ds,
                                        const size_t              idx_npot_geom_010_1_dp,
                                        const size_t              idx_npot_geom_020_0_dp,
                                        const size_t              idx_npot_geom_020_1_dp,
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

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds);

    auto ta2_xx_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 3);

    auto ta2_xx_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 5);

    auto ta2_xy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 6);

    auto ta2_xy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 9);

    auto ta2_xy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 11);

    auto ta2_xz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 12);

    auto ta2_xz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 15);

    auto ta2_xz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 17);

    auto ta2_yy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 18);

    auto ta2_yy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 21);

    auto ta2_yy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 23);

    auto ta2_yz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 24);

    auto ta2_yz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 27);

    auto ta2_yz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 29);

    auto ta2_zz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 30);

    auto ta2_zz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 33);

    auto ta2_zz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 35);

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds);

    auto ta2_xx_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 3);

    auto ta2_xx_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 5);

    auto ta2_xy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 6);

    auto ta2_xy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 9);

    auto ta2_xy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 11);

    auto ta2_xz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 12);

    auto ta2_xz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 15);

    auto ta2_xz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 17);

    auto ta2_yy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 18);

    auto ta2_yy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 21);

    auto ta2_yy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 23);

    auto ta2_yz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 24);

    auto ta2_yz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 27);

    auto ta2_yz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 29);

    auto ta2_zz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 30);

    auto ta2_zz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 33);

    auto ta2_zz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 35);

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

    auto ta1_y_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 22);

    auto ta1_y_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 27);

    auto ta1_y_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 28);

    auto ta1_y_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 29);

    auto ta1_y_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 32);

    auto ta1_y_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 33);

    auto ta1_y_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 34);

    auto ta1_y_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 35);

    auto ta1_z_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 36);

    auto ta1_z_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 37);

    auto ta1_z_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 38);

    auto ta1_z_xz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 42);

    auto ta1_z_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 44);

    auto ta1_z_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 45);

    auto ta1_z_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 46);

    auto ta1_z_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 47);

    auto ta1_z_yz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 49);

    auto ta1_z_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 50);

    auto ta1_z_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 51);

    auto ta1_z_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 52);

    auto ta1_z_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 53);

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp);

    auto ta2_xx_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 1);

    auto ta2_xx_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 2);

    auto ta2_xx_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 3);

    auto ta2_xx_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 4);

    auto ta2_xx_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 6);

    auto ta2_xx_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 8);

    auto ta2_xx_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 9);

    auto ta2_xx_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 10);

    auto ta2_xx_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 11);

    auto ta2_xx_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 14);

    auto ta2_xx_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 15);

    auto ta2_xx_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 16);

    auto ta2_xx_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 17);

    auto ta2_xy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 18);

    auto ta2_xy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 19);

    auto ta2_xy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 20);

    auto ta2_xy_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 21);

    auto ta2_xy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 22);

    auto ta2_xy_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 24);

    auto ta2_xy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 27);

    auto ta2_xy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 28);

    auto ta2_xy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 29);

    auto ta2_xy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 31);

    auto ta2_xy_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 32);

    auto ta2_xy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 33);

    auto ta2_xy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 34);

    auto ta2_xy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 35);

    auto ta2_xz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 36);

    auto ta2_xz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 37);

    auto ta2_xz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 38);

    auto ta2_xz_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 39);

    auto ta2_xz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 42);

    auto ta2_xz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 44);

    auto ta2_xz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 45);

    auto ta2_xz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 46);

    auto ta2_xz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 47);

    auto ta2_xz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 49);

    auto ta2_xz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 50);

    auto ta2_xz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 51);

    auto ta2_xz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 52);

    auto ta2_xz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 53);

    auto ta2_yy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 54);

    auto ta2_yy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 55);

    auto ta2_yy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 56);

    auto ta2_yy_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 57);

    auto ta2_yy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 58);

    auto ta2_yy_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 62);

    auto ta2_yy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 63);

    auto ta2_yy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 64);

    auto ta2_yy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 65);

    auto ta2_yy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 67);

    auto ta2_yy_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 68);

    auto ta2_yy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 69);

    auto ta2_yy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 70);

    auto ta2_yy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 71);

    auto ta2_yz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 72);

    auto ta2_yz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 73);

    auto ta2_yz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 74);

    auto ta2_yz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 76);

    auto ta2_yz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 78);

    auto ta2_yz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 80);

    auto ta2_yz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 81);

    auto ta2_yz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 82);

    auto ta2_yz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 83);

    auto ta2_yz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 85);

    auto ta2_yz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 86);

    auto ta2_yz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 87);

    auto ta2_yz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 88);

    auto ta2_yz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 89);

    auto ta2_zz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 90);

    auto ta2_zz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 91);

    auto ta2_zz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 92);

    auto ta2_zz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 94);

    auto ta2_zz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 96);

    auto ta2_zz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 98);

    auto ta2_zz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 99);

    auto ta2_zz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 100);

    auto ta2_zz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 101);

    auto ta2_zz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 103);

    auto ta2_zz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 104);

    auto ta2_zz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 105);

    auto ta2_zz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 106);

    auto ta2_zz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 107);

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp);

    auto ta2_xx_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 1);

    auto ta2_xx_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 2);

    auto ta2_xx_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 3);

    auto ta2_xx_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 4);

    auto ta2_xx_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 6);

    auto ta2_xx_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 8);

    auto ta2_xx_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 9);

    auto ta2_xx_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 10);

    auto ta2_xx_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 11);

    auto ta2_xx_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 14);

    auto ta2_xx_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 15);

    auto ta2_xx_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 16);

    auto ta2_xx_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 17);

    auto ta2_xy_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 18);

    auto ta2_xy_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 19);

    auto ta2_xy_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 20);

    auto ta2_xy_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 21);

    auto ta2_xy_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 22);

    auto ta2_xy_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 24);

    auto ta2_xy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 27);

    auto ta2_xy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 28);

    auto ta2_xy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 29);

    auto ta2_xy_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 31);

    auto ta2_xy_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 32);

    auto ta2_xy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 33);

    auto ta2_xy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 34);

    auto ta2_xy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 35);

    auto ta2_xz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 36);

    auto ta2_xz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 37);

    auto ta2_xz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 38);

    auto ta2_xz_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 39);

    auto ta2_xz_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 42);

    auto ta2_xz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 44);

    auto ta2_xz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 45);

    auto ta2_xz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 46);

    auto ta2_xz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 47);

    auto ta2_xz_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 49);

    auto ta2_xz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 50);

    auto ta2_xz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 51);

    auto ta2_xz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 52);

    auto ta2_xz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 53);

    auto ta2_yy_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 54);

    auto ta2_yy_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 55);

    auto ta2_yy_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 56);

    auto ta2_yy_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 57);

    auto ta2_yy_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 58);

    auto ta2_yy_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 62);

    auto ta2_yy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 63);

    auto ta2_yy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 64);

    auto ta2_yy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 65);

    auto ta2_yy_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 67);

    auto ta2_yy_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 68);

    auto ta2_yy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 69);

    auto ta2_yy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 70);

    auto ta2_yy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 71);

    auto ta2_yz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 72);

    auto ta2_yz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 73);

    auto ta2_yz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 74);

    auto ta2_yz_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 76);

    auto ta2_yz_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 78);

    auto ta2_yz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 80);

    auto ta2_yz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 81);

    auto ta2_yz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 82);

    auto ta2_yz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 83);

    auto ta2_yz_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 85);

    auto ta2_yz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 86);

    auto ta2_yz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 87);

    auto ta2_yz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 88);

    auto ta2_yz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 89);

    auto ta2_zz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 90);

    auto ta2_zz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 91);

    auto ta2_zz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 92);

    auto ta2_zz_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 94);

    auto ta2_zz_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 96);

    auto ta2_zz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 98);

    auto ta2_zz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 99);

    auto ta2_zz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 100);

    auto ta2_zz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 101);

    auto ta2_zz_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 103);

    auto ta2_zz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 104);

    auto ta2_zz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 105);

    auto ta2_zz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 106);

    auto ta2_zz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 107);

    // Set up 0-3 components of targeted buffer : FP

    auto ta2_xx_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp);

    auto ta2_xx_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 1);

    auto ta2_xx_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 2);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xx_z_1,   \
                             ta2_xx_x_x_0,   \
                             ta2_xx_x_x_1,   \
                             ta2_xx_x_y_0,   \
                             ta2_xx_x_y_1,   \
                             ta2_xx_x_z_0,   \
                             ta2_xx_x_z_1,   \
                             ta2_xx_xx_0_0,  \
                             ta2_xx_xx_0_1,  \
                             ta2_xx_xx_x_0,  \
                             ta2_xx_xx_x_1,  \
                             ta2_xx_xx_y_0,  \
                             ta2_xx_xx_y_1,  \
                             ta2_xx_xx_z_0,  \
                             ta2_xx_xx_z_1,  \
                             ta2_xx_xxx_x_0, \
                             ta2_xx_xxx_y_0, \
                             ta2_xx_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxx_x_0[i] = 2.0 * ta2_xx_x_x_0[i] * fe_0 - 2.0 * ta2_xx_x_x_1[i] * fe_0 + ta2_xx_xx_0_0[i] * fe_0 - ta2_xx_xx_0_1[i] * fe_0 +
                            2.0 * ta1_x_xx_x_1[i] + ta2_xx_xx_x_0[i] * pa_x[i] - ta2_xx_xx_x_1[i] * pc_x[i];

        ta2_xx_xxx_y_0[i] = 2.0 * ta2_xx_x_y_0[i] * fe_0 - 2.0 * ta2_xx_x_y_1[i] * fe_0 + 2.0 * ta1_x_xx_y_1[i] + ta2_xx_xx_y_0[i] * pa_x[i] -
                            ta2_xx_xx_y_1[i] * pc_x[i];

        ta2_xx_xxx_z_0[i] = 2.0 * ta2_xx_x_z_0[i] * fe_0 - 2.0 * ta2_xx_x_z_1[i] * fe_0 + 2.0 * ta1_x_xx_z_1[i] + ta2_xx_xx_z_0[i] * pa_x[i] -
                            ta2_xx_xx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto ta2_xx_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 3);

    auto ta2_xx_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 4);

    auto ta2_xx_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 5);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xx_xx_0_0,  \
                             ta2_xx_xx_0_1,  \
                             ta2_xx_xx_x_0,  \
                             ta2_xx_xx_x_1,  \
                             ta2_xx_xx_y_0,  \
                             ta2_xx_xx_y_1,  \
                             ta2_xx_xx_z_0,  \
                             ta2_xx_xx_z_1,  \
                             ta2_xx_xxy_x_0, \
                             ta2_xx_xxy_y_0, \
                             ta2_xx_xxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxy_x_0[i] = ta2_xx_xx_x_0[i] * pa_y[i] - ta2_xx_xx_x_1[i] * pc_y[i];

        ta2_xx_xxy_y_0[i] = ta2_xx_xx_0_0[i] * fe_0 - ta2_xx_xx_0_1[i] * fe_0 + ta2_xx_xx_y_0[i] * pa_y[i] - ta2_xx_xx_y_1[i] * pc_y[i];

        ta2_xx_xxy_z_0[i] = ta2_xx_xx_z_0[i] * pa_y[i] - ta2_xx_xx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto ta2_xx_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 6);

    auto ta2_xx_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 7);

    auto ta2_xx_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 8);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xx_xx_0_0,  \
                             ta2_xx_xx_0_1,  \
                             ta2_xx_xx_x_0,  \
                             ta2_xx_xx_x_1,  \
                             ta2_xx_xx_y_0,  \
                             ta2_xx_xx_y_1,  \
                             ta2_xx_xx_z_0,  \
                             ta2_xx_xx_z_1,  \
                             ta2_xx_xxz_x_0, \
                             ta2_xx_xxz_y_0, \
                             ta2_xx_xxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxz_x_0[i] = ta2_xx_xx_x_0[i] * pa_z[i] - ta2_xx_xx_x_1[i] * pc_z[i];

        ta2_xx_xxz_y_0[i] = ta2_xx_xx_y_0[i] * pa_z[i] - ta2_xx_xx_y_1[i] * pc_z[i];

        ta2_xx_xxz_z_0[i] = ta2_xx_xx_0_0[i] * fe_0 - ta2_xx_xx_0_1[i] * fe_0 + ta2_xx_xx_z_0[i] * pa_z[i] - ta2_xx_xx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto ta2_xx_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 9);

    auto ta2_xx_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 10);

    auto ta2_xx_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 11);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_z_1,   \
                             ta2_xx_x_x_0,   \
                             ta2_xx_x_x_1,   \
                             ta2_xx_xy_x_0,  \
                             ta2_xx_xy_x_1,  \
                             ta2_xx_xyy_x_0, \
                             ta2_xx_xyy_y_0, \
                             ta2_xx_xyy_z_0, \
                             ta2_xx_yy_y_0,  \
                             ta2_xx_yy_y_1,  \
                             ta2_xx_yy_z_0,  \
                             ta2_xx_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyy_x_0[i] = ta2_xx_x_x_0[i] * fe_0 - ta2_xx_x_x_1[i] * fe_0 + ta2_xx_xy_x_0[i] * pa_y[i] - ta2_xx_xy_x_1[i] * pc_y[i];

        ta2_xx_xyy_y_0[i] = 2.0 * ta1_x_yy_y_1[i] + ta2_xx_yy_y_0[i] * pa_x[i] - ta2_xx_yy_y_1[i] * pc_x[i];

        ta2_xx_xyy_z_0[i] = 2.0 * ta1_x_yy_z_1[i] + ta2_xx_yy_z_0[i] * pa_x[i] - ta2_xx_yy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto ta2_xx_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 12);

    auto ta2_xx_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 13);

    auto ta2_xx_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 14);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta2_xx_xy_y_0,  \
                             ta2_xx_xy_y_1,  \
                             ta2_xx_xyz_x_0, \
                             ta2_xx_xyz_y_0, \
                             ta2_xx_xyz_z_0, \
                             ta2_xx_xz_x_0,  \
                             ta2_xx_xz_x_1,  \
                             ta2_xx_xz_z_0,  \
                             ta2_xx_xz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xyz_x_0[i] = ta2_xx_xz_x_0[i] * pa_y[i] - ta2_xx_xz_x_1[i] * pc_y[i];

        ta2_xx_xyz_y_0[i] = ta2_xx_xy_y_0[i] * pa_z[i] - ta2_xx_xy_y_1[i] * pc_z[i];

        ta2_xx_xyz_z_0[i] = ta2_xx_xz_z_0[i] * pa_y[i] - ta2_xx_xz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto ta2_xx_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 15);

    auto ta2_xx_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 16);

    auto ta2_xx_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 17);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_zz_y_1,   \
                             ta1_x_zz_z_1,   \
                             ta2_xx_x_x_0,   \
                             ta2_xx_x_x_1,   \
                             ta2_xx_xz_x_0,  \
                             ta2_xx_xz_x_1,  \
                             ta2_xx_xzz_x_0, \
                             ta2_xx_xzz_y_0, \
                             ta2_xx_xzz_z_0, \
                             ta2_xx_zz_y_0,  \
                             ta2_xx_zz_y_1,  \
                             ta2_xx_zz_z_0,  \
                             ta2_xx_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzz_x_0[i] = ta2_xx_x_x_0[i] * fe_0 - ta2_xx_x_x_1[i] * fe_0 + ta2_xx_xz_x_0[i] * pa_z[i] - ta2_xx_xz_x_1[i] * pc_z[i];

        ta2_xx_xzz_y_0[i] = 2.0 * ta1_x_zz_y_1[i] + ta2_xx_zz_y_0[i] * pa_x[i] - ta2_xx_zz_y_1[i] * pc_x[i];

        ta2_xx_xzz_z_0[i] = 2.0 * ta1_x_zz_z_1[i] + ta2_xx_zz_z_0[i] * pa_x[i] - ta2_xx_zz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto ta2_xx_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 18);

    auto ta2_xx_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 19);

    auto ta2_xx_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 20);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xx_y_x_0,   \
                             ta2_xx_y_x_1,   \
                             ta2_xx_y_y_0,   \
                             ta2_xx_y_y_1,   \
                             ta2_xx_y_z_0,   \
                             ta2_xx_y_z_1,   \
                             ta2_xx_yy_0_0,  \
                             ta2_xx_yy_0_1,  \
                             ta2_xx_yy_x_0,  \
                             ta2_xx_yy_x_1,  \
                             ta2_xx_yy_y_0,  \
                             ta2_xx_yy_y_1,  \
                             ta2_xx_yy_z_0,  \
                             ta2_xx_yy_z_1,  \
                             ta2_xx_yyy_x_0, \
                             ta2_xx_yyy_y_0, \
                             ta2_xx_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyy_x_0[i] = 2.0 * ta2_xx_y_x_0[i] * fe_0 - 2.0 * ta2_xx_y_x_1[i] * fe_0 + ta2_xx_yy_x_0[i] * pa_y[i] - ta2_xx_yy_x_1[i] * pc_y[i];

        ta2_xx_yyy_y_0[i] = 2.0 * ta2_xx_y_y_0[i] * fe_0 - 2.0 * ta2_xx_y_y_1[i] * fe_0 + ta2_xx_yy_0_0[i] * fe_0 - ta2_xx_yy_0_1[i] * fe_0 +
                            ta2_xx_yy_y_0[i] * pa_y[i] - ta2_xx_yy_y_1[i] * pc_y[i];

        ta2_xx_yyy_z_0[i] = 2.0 * ta2_xx_y_z_0[i] * fe_0 - 2.0 * ta2_xx_y_z_1[i] * fe_0 + ta2_xx_yy_z_0[i] * pa_y[i] - ta2_xx_yy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto ta2_xx_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 21);

    auto ta2_xx_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 22);

    auto ta2_xx_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 23);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta2_xx_yy_x_0,  \
                             ta2_xx_yy_x_1,  \
                             ta2_xx_yy_y_0,  \
                             ta2_xx_yy_y_1,  \
                             ta2_xx_yyz_x_0, \
                             ta2_xx_yyz_y_0, \
                             ta2_xx_yyz_z_0, \
                             ta2_xx_yz_z_0,  \
                             ta2_xx_yz_z_1,  \
                             ta2_xx_z_z_0,   \
                             ta2_xx_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyz_x_0[i] = ta2_xx_yy_x_0[i] * pa_z[i] - ta2_xx_yy_x_1[i] * pc_z[i];

        ta2_xx_yyz_y_0[i] = ta2_xx_yy_y_0[i] * pa_z[i] - ta2_xx_yy_y_1[i] * pc_z[i];

        ta2_xx_yyz_z_0[i] = ta2_xx_z_z_0[i] * fe_0 - ta2_xx_z_z_1[i] * fe_0 + ta2_xx_yz_z_0[i] * pa_y[i] - ta2_xx_yz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto ta2_xx_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 24);

    auto ta2_xx_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 25);

    auto ta2_xx_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 26);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xx_yzz_x_0, \
                             ta2_xx_yzz_y_0, \
                             ta2_xx_yzz_z_0, \
                             ta2_xx_zz_0_0,  \
                             ta2_xx_zz_0_1,  \
                             ta2_xx_zz_x_0,  \
                             ta2_xx_zz_x_1,  \
                             ta2_xx_zz_y_0,  \
                             ta2_xx_zz_y_1,  \
                             ta2_xx_zz_z_0,  \
                             ta2_xx_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzz_x_0[i] = ta2_xx_zz_x_0[i] * pa_y[i] - ta2_xx_zz_x_1[i] * pc_y[i];

        ta2_xx_yzz_y_0[i] = ta2_xx_zz_0_0[i] * fe_0 - ta2_xx_zz_0_1[i] * fe_0 + ta2_xx_zz_y_0[i] * pa_y[i] - ta2_xx_zz_y_1[i] * pc_y[i];

        ta2_xx_yzz_z_0[i] = ta2_xx_zz_z_0[i] * pa_y[i] - ta2_xx_zz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto ta2_xx_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 27);

    auto ta2_xx_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 28);

    auto ta2_xx_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 29);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xx_z_x_0,   \
                             ta2_xx_z_x_1,   \
                             ta2_xx_z_y_0,   \
                             ta2_xx_z_y_1,   \
                             ta2_xx_z_z_0,   \
                             ta2_xx_z_z_1,   \
                             ta2_xx_zz_0_0,  \
                             ta2_xx_zz_0_1,  \
                             ta2_xx_zz_x_0,  \
                             ta2_xx_zz_x_1,  \
                             ta2_xx_zz_y_0,  \
                             ta2_xx_zz_y_1,  \
                             ta2_xx_zz_z_0,  \
                             ta2_xx_zz_z_1,  \
                             ta2_xx_zzz_x_0, \
                             ta2_xx_zzz_y_0, \
                             ta2_xx_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzz_x_0[i] = 2.0 * ta2_xx_z_x_0[i] * fe_0 - 2.0 * ta2_xx_z_x_1[i] * fe_0 + ta2_xx_zz_x_0[i] * pa_z[i] - ta2_xx_zz_x_1[i] * pc_z[i];

        ta2_xx_zzz_y_0[i] = 2.0 * ta2_xx_z_y_0[i] * fe_0 - 2.0 * ta2_xx_z_y_1[i] * fe_0 + ta2_xx_zz_y_0[i] * pa_z[i] - ta2_xx_zz_y_1[i] * pc_z[i];

        ta2_xx_zzz_z_0[i] = 2.0 * ta2_xx_z_z_0[i] * fe_0 - 2.0 * ta2_xx_z_z_1[i] * fe_0 + ta2_xx_zz_0_0[i] * fe_0 - ta2_xx_zz_0_1[i] * fe_0 +
                            ta2_xx_zz_z_0[i] * pa_z[i] - ta2_xx_zz_z_1[i] * pc_z[i];
    }

    // Set up 30-33 components of targeted buffer : FP

    auto ta2_xy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 30);

    auto ta2_xy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 31);

    auto ta2_xy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 32);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_y_1,   \
                             ta1_y_xx_z_1,   \
                             ta2_xy_x_x_0,   \
                             ta2_xy_x_x_1,   \
                             ta2_xy_x_y_0,   \
                             ta2_xy_x_y_1,   \
                             ta2_xy_x_z_0,   \
                             ta2_xy_x_z_1,   \
                             ta2_xy_xx_0_0,  \
                             ta2_xy_xx_0_1,  \
                             ta2_xy_xx_x_0,  \
                             ta2_xy_xx_x_1,  \
                             ta2_xy_xx_y_0,  \
                             ta2_xy_xx_y_1,  \
                             ta2_xy_xx_z_0,  \
                             ta2_xy_xx_z_1,  \
                             ta2_xy_xxx_x_0, \
                             ta2_xy_xxx_y_0, \
                             ta2_xy_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxx_x_0[i] = 2.0 * ta2_xy_x_x_0[i] * fe_0 - 2.0 * ta2_xy_x_x_1[i] * fe_0 + ta2_xy_xx_0_0[i] * fe_0 - ta2_xy_xx_0_1[i] * fe_0 +
                            ta1_y_xx_x_1[i] + ta2_xy_xx_x_0[i] * pa_x[i] - ta2_xy_xx_x_1[i] * pc_x[i];

        ta2_xy_xxx_y_0[i] =
            2.0 * ta2_xy_x_y_0[i] * fe_0 - 2.0 * ta2_xy_x_y_1[i] * fe_0 + ta1_y_xx_y_1[i] + ta2_xy_xx_y_0[i] * pa_x[i] - ta2_xy_xx_y_1[i] * pc_x[i];

        ta2_xy_xxx_z_0[i] =
            2.0 * ta2_xy_x_z_0[i] * fe_0 - 2.0 * ta2_xy_x_z_1[i] * fe_0 + ta1_y_xx_z_1[i] + ta2_xy_xx_z_0[i] * pa_x[i] - ta2_xy_xx_z_1[i] * pc_x[i];
    }

    // Set up 33-36 components of targeted buffer : FP

    auto ta2_xy_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 33);

    auto ta2_xy_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 34);

    auto ta2_xy_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 35);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_z_1,   \
                             ta1_y_xy_y_1,   \
                             ta2_xy_xx_x_0,  \
                             ta2_xy_xx_x_1,  \
                             ta2_xy_xx_z_0,  \
                             ta2_xy_xx_z_1,  \
                             ta2_xy_xxy_x_0, \
                             ta2_xy_xxy_y_0, \
                             ta2_xy_xxy_z_0, \
                             ta2_xy_xy_y_0,  \
                             ta2_xy_xy_y_1,  \
                             ta2_xy_y_y_0,   \
                             ta2_xy_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxy_x_0[i] = ta1_x_xx_x_1[i] + ta2_xy_xx_x_0[i] * pa_y[i] - ta2_xy_xx_x_1[i] * pc_y[i];

        ta2_xy_xxy_y_0[i] =
            ta2_xy_y_y_0[i] * fe_0 - ta2_xy_y_y_1[i] * fe_0 + ta1_y_xy_y_1[i] + ta2_xy_xy_y_0[i] * pa_x[i] - ta2_xy_xy_y_1[i] * pc_x[i];

        ta2_xy_xxy_z_0[i] = ta1_x_xx_z_1[i] + ta2_xy_xx_z_0[i] * pa_y[i] - ta2_xy_xx_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : FP

    auto ta2_xy_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 36);

    auto ta2_xy_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 37);

    auto ta2_xy_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 38);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xy_xx_0_0,  \
                             ta2_xy_xx_0_1,  \
                             ta2_xy_xx_x_0,  \
                             ta2_xy_xx_x_1,  \
                             ta2_xy_xx_y_0,  \
                             ta2_xy_xx_y_1,  \
                             ta2_xy_xx_z_0,  \
                             ta2_xy_xx_z_1,  \
                             ta2_xy_xxz_x_0, \
                             ta2_xy_xxz_y_0, \
                             ta2_xy_xxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxz_x_0[i] = ta2_xy_xx_x_0[i] * pa_z[i] - ta2_xy_xx_x_1[i] * pc_z[i];

        ta2_xy_xxz_y_0[i] = ta2_xy_xx_y_0[i] * pa_z[i] - ta2_xy_xx_y_1[i] * pc_z[i];

        ta2_xy_xxz_z_0[i] = ta2_xy_xx_0_0[i] * fe_0 - ta2_xy_xx_0_1[i] * fe_0 + ta2_xy_xx_z_0[i] * pa_z[i] - ta2_xy_xx_z_1[i] * pc_z[i];
    }

    // Set up 39-42 components of targeted buffer : FP

    auto ta2_xy_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 39);

    auto ta2_xy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 40);

    auto ta2_xy_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 41);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_z_1,   \
                             ta2_xy_xyy_x_0, \
                             ta2_xy_xyy_y_0, \
                             ta2_xy_xyy_z_0, \
                             ta2_xy_yy_0_0,  \
                             ta2_xy_yy_0_1,  \
                             ta2_xy_yy_x_0,  \
                             ta2_xy_yy_x_1,  \
                             ta2_xy_yy_y_0,  \
                             ta2_xy_yy_y_1,  \
                             ta2_xy_yy_z_0,  \
                             ta2_xy_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyy_x_0[i] =
            ta2_xy_yy_0_0[i] * fe_0 - ta2_xy_yy_0_1[i] * fe_0 + ta1_y_yy_x_1[i] + ta2_xy_yy_x_0[i] * pa_x[i] - ta2_xy_yy_x_1[i] * pc_x[i];

        ta2_xy_xyy_y_0[i] = ta1_y_yy_y_1[i] + ta2_xy_yy_y_0[i] * pa_x[i] - ta2_xy_yy_y_1[i] * pc_x[i];

        ta2_xy_xyy_z_0[i] = ta1_y_yy_z_1[i] + ta2_xy_yy_z_0[i] * pa_x[i] - ta2_xy_yy_z_1[i] * pc_x[i];
    }

    // Set up 42-45 components of targeted buffer : FP

    auto ta2_xy_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 42);

    auto ta2_xy_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 43);

    auto ta2_xy_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 44);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_yz_z_1,   \
                             ta2_xy_xy_x_0,  \
                             ta2_xy_xy_x_1,  \
                             ta2_xy_xy_y_0,  \
                             ta2_xy_xy_y_1,  \
                             ta2_xy_xyz_x_0, \
                             ta2_xy_xyz_y_0, \
                             ta2_xy_xyz_z_0, \
                             ta2_xy_yz_z_0,  \
                             ta2_xy_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xyz_x_0[i] = ta2_xy_xy_x_0[i] * pa_z[i] - ta2_xy_xy_x_1[i] * pc_z[i];

        ta2_xy_xyz_y_0[i] = ta2_xy_xy_y_0[i] * pa_z[i] - ta2_xy_xy_y_1[i] * pc_z[i];

        ta2_xy_xyz_z_0[i] = ta1_y_yz_z_1[i] + ta2_xy_yz_z_0[i] * pa_x[i] - ta2_xy_yz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : FP

    auto ta2_xy_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 45);

    auto ta2_xy_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 46);

    auto ta2_xy_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 47);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_z_1,   \
                             ta2_xy_x_x_0,   \
                             ta2_xy_x_x_1,   \
                             ta2_xy_xz_x_0,  \
                             ta2_xy_xz_x_1,  \
                             ta2_xy_xzz_x_0, \
                             ta2_xy_xzz_y_0, \
                             ta2_xy_xzz_z_0, \
                             ta2_xy_zz_y_0,  \
                             ta2_xy_zz_y_1,  \
                             ta2_xy_zz_z_0,  \
                             ta2_xy_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzz_x_0[i] = ta2_xy_x_x_0[i] * fe_0 - ta2_xy_x_x_1[i] * fe_0 + ta2_xy_xz_x_0[i] * pa_z[i] - ta2_xy_xz_x_1[i] * pc_z[i];

        ta2_xy_xzz_y_0[i] = ta1_y_zz_y_1[i] + ta2_xy_zz_y_0[i] * pa_x[i] - ta2_xy_zz_y_1[i] * pc_x[i];

        ta2_xy_xzz_z_0[i] = ta1_y_zz_z_1[i] + ta2_xy_zz_z_0[i] * pa_x[i] - ta2_xy_zz_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : FP

    auto ta2_xy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 48);

    auto ta2_xy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 49);

    auto ta2_xy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 50);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_yy_x_1,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_z_1,   \
                             ta2_xy_y_x_0,   \
                             ta2_xy_y_x_1,   \
                             ta2_xy_y_y_0,   \
                             ta2_xy_y_y_1,   \
                             ta2_xy_y_z_0,   \
                             ta2_xy_y_z_1,   \
                             ta2_xy_yy_0_0,  \
                             ta2_xy_yy_0_1,  \
                             ta2_xy_yy_x_0,  \
                             ta2_xy_yy_x_1,  \
                             ta2_xy_yy_y_0,  \
                             ta2_xy_yy_y_1,  \
                             ta2_xy_yy_z_0,  \
                             ta2_xy_yy_z_1,  \
                             ta2_xy_yyy_x_0, \
                             ta2_xy_yyy_y_0, \
                             ta2_xy_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyy_x_0[i] =
            2.0 * ta2_xy_y_x_0[i] * fe_0 - 2.0 * ta2_xy_y_x_1[i] * fe_0 + ta1_x_yy_x_1[i] + ta2_xy_yy_x_0[i] * pa_y[i] - ta2_xy_yy_x_1[i] * pc_y[i];

        ta2_xy_yyy_y_0[i] = 2.0 * ta2_xy_y_y_0[i] * fe_0 - 2.0 * ta2_xy_y_y_1[i] * fe_0 + ta2_xy_yy_0_0[i] * fe_0 - ta2_xy_yy_0_1[i] * fe_0 +
                            ta1_x_yy_y_1[i] + ta2_xy_yy_y_0[i] * pa_y[i] - ta2_xy_yy_y_1[i] * pc_y[i];

        ta2_xy_yyy_z_0[i] =
            2.0 * ta2_xy_y_z_0[i] * fe_0 - 2.0 * ta2_xy_y_z_1[i] * fe_0 + ta1_x_yy_z_1[i] + ta2_xy_yy_z_0[i] * pa_y[i] - ta2_xy_yy_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : FP

    auto ta2_xy_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 51);

    auto ta2_xy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 52);

    auto ta2_xy_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 53);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xy_yy_0_0,  \
                             ta2_xy_yy_0_1,  \
                             ta2_xy_yy_x_0,  \
                             ta2_xy_yy_x_1,  \
                             ta2_xy_yy_y_0,  \
                             ta2_xy_yy_y_1,  \
                             ta2_xy_yy_z_0,  \
                             ta2_xy_yy_z_1,  \
                             ta2_xy_yyz_x_0, \
                             ta2_xy_yyz_y_0, \
                             ta2_xy_yyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyz_x_0[i] = ta2_xy_yy_x_0[i] * pa_z[i] - ta2_xy_yy_x_1[i] * pc_z[i];

        ta2_xy_yyz_y_0[i] = ta2_xy_yy_y_0[i] * pa_z[i] - ta2_xy_yy_y_1[i] * pc_z[i];

        ta2_xy_yyz_z_0[i] = ta2_xy_yy_0_0[i] * fe_0 - ta2_xy_yy_0_1[i] * fe_0 + ta2_xy_yy_z_0[i] * pa_z[i] - ta2_xy_yy_z_1[i] * pc_z[i];
    }

    // Set up 54-57 components of targeted buffer : FP

    auto ta2_xy_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 54);

    auto ta2_xy_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 55);

    auto ta2_xy_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 56);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_z_1,   \
                             ta2_xy_y_y_0,   \
                             ta2_xy_y_y_1,   \
                             ta2_xy_yz_y_0,  \
                             ta2_xy_yz_y_1,  \
                             ta2_xy_yzz_x_0, \
                             ta2_xy_yzz_y_0, \
                             ta2_xy_yzz_z_0, \
                             ta2_xy_zz_x_0,  \
                             ta2_xy_zz_x_1,  \
                             ta2_xy_zz_z_0,  \
                             ta2_xy_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzz_x_0[i] = ta1_x_zz_x_1[i] + ta2_xy_zz_x_0[i] * pa_y[i] - ta2_xy_zz_x_1[i] * pc_y[i];

        ta2_xy_yzz_y_0[i] = ta2_xy_y_y_0[i] * fe_0 - ta2_xy_y_y_1[i] * fe_0 + ta2_xy_yz_y_0[i] * pa_z[i] - ta2_xy_yz_y_1[i] * pc_z[i];

        ta2_xy_yzz_z_0[i] = ta1_x_zz_z_1[i] + ta2_xy_zz_z_0[i] * pa_y[i] - ta2_xy_zz_z_1[i] * pc_y[i];
    }

    // Set up 57-60 components of targeted buffer : FP

    auto ta2_xy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 57);

    auto ta2_xy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 58);

    auto ta2_xy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 59);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xy_z_x_0,   \
                             ta2_xy_z_x_1,   \
                             ta2_xy_z_y_0,   \
                             ta2_xy_z_y_1,   \
                             ta2_xy_z_z_0,   \
                             ta2_xy_z_z_1,   \
                             ta2_xy_zz_0_0,  \
                             ta2_xy_zz_0_1,  \
                             ta2_xy_zz_x_0,  \
                             ta2_xy_zz_x_1,  \
                             ta2_xy_zz_y_0,  \
                             ta2_xy_zz_y_1,  \
                             ta2_xy_zz_z_0,  \
                             ta2_xy_zz_z_1,  \
                             ta2_xy_zzz_x_0, \
                             ta2_xy_zzz_y_0, \
                             ta2_xy_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzz_x_0[i] = 2.0 * ta2_xy_z_x_0[i] * fe_0 - 2.0 * ta2_xy_z_x_1[i] * fe_0 + ta2_xy_zz_x_0[i] * pa_z[i] - ta2_xy_zz_x_1[i] * pc_z[i];

        ta2_xy_zzz_y_0[i] = 2.0 * ta2_xy_z_y_0[i] * fe_0 - 2.0 * ta2_xy_z_y_1[i] * fe_0 + ta2_xy_zz_y_0[i] * pa_z[i] - ta2_xy_zz_y_1[i] * pc_z[i];

        ta2_xy_zzz_z_0[i] = 2.0 * ta2_xy_z_z_0[i] * fe_0 - 2.0 * ta2_xy_z_z_1[i] * fe_0 + ta2_xy_zz_0_0[i] * fe_0 - ta2_xy_zz_0_1[i] * fe_0 +
                            ta2_xy_zz_z_0[i] * pa_z[i] - ta2_xy_zz_z_1[i] * pc_z[i];
    }

    // Set up 60-63 components of targeted buffer : FP

    auto ta2_xz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 60);

    auto ta2_xz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 61);

    auto ta2_xz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 62);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_y_1,   \
                             ta1_z_xx_z_1,   \
                             ta2_xz_x_x_0,   \
                             ta2_xz_x_x_1,   \
                             ta2_xz_x_y_0,   \
                             ta2_xz_x_y_1,   \
                             ta2_xz_x_z_0,   \
                             ta2_xz_x_z_1,   \
                             ta2_xz_xx_0_0,  \
                             ta2_xz_xx_0_1,  \
                             ta2_xz_xx_x_0,  \
                             ta2_xz_xx_x_1,  \
                             ta2_xz_xx_y_0,  \
                             ta2_xz_xx_y_1,  \
                             ta2_xz_xx_z_0,  \
                             ta2_xz_xx_z_1,  \
                             ta2_xz_xxx_x_0, \
                             ta2_xz_xxx_y_0, \
                             ta2_xz_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxx_x_0[i] = 2.0 * ta2_xz_x_x_0[i] * fe_0 - 2.0 * ta2_xz_x_x_1[i] * fe_0 + ta2_xz_xx_0_0[i] * fe_0 - ta2_xz_xx_0_1[i] * fe_0 +
                            ta1_z_xx_x_1[i] + ta2_xz_xx_x_0[i] * pa_x[i] - ta2_xz_xx_x_1[i] * pc_x[i];

        ta2_xz_xxx_y_0[i] =
            2.0 * ta2_xz_x_y_0[i] * fe_0 - 2.0 * ta2_xz_x_y_1[i] * fe_0 + ta1_z_xx_y_1[i] + ta2_xz_xx_y_0[i] * pa_x[i] - ta2_xz_xx_y_1[i] * pc_x[i];

        ta2_xz_xxx_z_0[i] =
            2.0 * ta2_xz_x_z_0[i] * fe_0 - 2.0 * ta2_xz_x_z_1[i] * fe_0 + ta1_z_xx_z_1[i] + ta2_xz_xx_z_0[i] * pa_x[i] - ta2_xz_xx_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : FP

    auto ta2_xz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 63);

    auto ta2_xz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 64);

    auto ta2_xz_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 65);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xz_xx_0_0,  \
                             ta2_xz_xx_0_1,  \
                             ta2_xz_xx_x_0,  \
                             ta2_xz_xx_x_1,  \
                             ta2_xz_xx_y_0,  \
                             ta2_xz_xx_y_1,  \
                             ta2_xz_xx_z_0,  \
                             ta2_xz_xx_z_1,  \
                             ta2_xz_xxy_x_0, \
                             ta2_xz_xxy_y_0, \
                             ta2_xz_xxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxy_x_0[i] = ta2_xz_xx_x_0[i] * pa_y[i] - ta2_xz_xx_x_1[i] * pc_y[i];

        ta2_xz_xxy_y_0[i] = ta2_xz_xx_0_0[i] * fe_0 - ta2_xz_xx_0_1[i] * fe_0 + ta2_xz_xx_y_0[i] * pa_y[i] - ta2_xz_xx_y_1[i] * pc_y[i];

        ta2_xz_xxy_z_0[i] = ta2_xz_xx_z_0[i] * pa_y[i] - ta2_xz_xx_z_1[i] * pc_y[i];
    }

    // Set up 66-69 components of targeted buffer : FP

    auto ta2_xz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 66);

    auto ta2_xz_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 67);

    auto ta2_xz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 68);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_y_1,   \
                             ta1_z_xz_z_1,   \
                             ta2_xz_xx_x_0,  \
                             ta2_xz_xx_x_1,  \
                             ta2_xz_xx_y_0,  \
                             ta2_xz_xx_y_1,  \
                             ta2_xz_xxz_x_0, \
                             ta2_xz_xxz_y_0, \
                             ta2_xz_xxz_z_0, \
                             ta2_xz_xz_z_0,  \
                             ta2_xz_xz_z_1,  \
                             ta2_xz_z_z_0,   \
                             ta2_xz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxz_x_0[i] = ta1_x_xx_x_1[i] + ta2_xz_xx_x_0[i] * pa_z[i] - ta2_xz_xx_x_1[i] * pc_z[i];

        ta2_xz_xxz_y_0[i] = ta1_x_xx_y_1[i] + ta2_xz_xx_y_0[i] * pa_z[i] - ta2_xz_xx_y_1[i] * pc_z[i];

        ta2_xz_xxz_z_0[i] =
            ta2_xz_z_z_0[i] * fe_0 - ta2_xz_z_z_1[i] * fe_0 + ta1_z_xz_z_1[i] + ta2_xz_xz_z_0[i] * pa_x[i] - ta2_xz_xz_z_1[i] * pc_x[i];
    }

    // Set up 69-72 components of targeted buffer : FP

    auto ta2_xz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 69);

    auto ta2_xz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 70);

    auto ta2_xz_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 71);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_z_1,   \
                             ta2_xz_x_x_0,   \
                             ta2_xz_x_x_1,   \
                             ta2_xz_xy_x_0,  \
                             ta2_xz_xy_x_1,  \
                             ta2_xz_xyy_x_0, \
                             ta2_xz_xyy_y_0, \
                             ta2_xz_xyy_z_0, \
                             ta2_xz_yy_y_0,  \
                             ta2_xz_yy_y_1,  \
                             ta2_xz_yy_z_0,  \
                             ta2_xz_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyy_x_0[i] = ta2_xz_x_x_0[i] * fe_0 - ta2_xz_x_x_1[i] * fe_0 + ta2_xz_xy_x_0[i] * pa_y[i] - ta2_xz_xy_x_1[i] * pc_y[i];

        ta2_xz_xyy_y_0[i] = ta1_z_yy_y_1[i] + ta2_xz_yy_y_0[i] * pa_x[i] - ta2_xz_yy_y_1[i] * pc_x[i];

        ta2_xz_xyy_z_0[i] = ta1_z_yy_z_1[i] + ta2_xz_yy_z_0[i] * pa_x[i] - ta2_xz_yy_z_1[i] * pc_x[i];
    }

    // Set up 72-75 components of targeted buffer : FP

    auto ta2_xz_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 72);

    auto ta2_xz_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 73);

    auto ta2_xz_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 74);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_yz_y_1,   \
                             ta2_xz_xyz_x_0, \
                             ta2_xz_xyz_y_0, \
                             ta2_xz_xyz_z_0, \
                             ta2_xz_xz_x_0,  \
                             ta2_xz_xz_x_1,  \
                             ta2_xz_xz_z_0,  \
                             ta2_xz_xz_z_1,  \
                             ta2_xz_yz_y_0,  \
                             ta2_xz_yz_y_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xyz_x_0[i] = ta2_xz_xz_x_0[i] * pa_y[i] - ta2_xz_xz_x_1[i] * pc_y[i];

        ta2_xz_xyz_y_0[i] = ta1_z_yz_y_1[i] + ta2_xz_yz_y_0[i] * pa_x[i] - ta2_xz_yz_y_1[i] * pc_x[i];

        ta2_xz_xyz_z_0[i] = ta2_xz_xz_z_0[i] * pa_y[i] - ta2_xz_xz_z_1[i] * pc_y[i];
    }

    // Set up 75-78 components of targeted buffer : FP

    auto ta2_xz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 75);

    auto ta2_xz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 76);

    auto ta2_xz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 77);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_z_1,   \
                             ta2_xz_xzz_x_0, \
                             ta2_xz_xzz_y_0, \
                             ta2_xz_xzz_z_0, \
                             ta2_xz_zz_0_0,  \
                             ta2_xz_zz_0_1,  \
                             ta2_xz_zz_x_0,  \
                             ta2_xz_zz_x_1,  \
                             ta2_xz_zz_y_0,  \
                             ta2_xz_zz_y_1,  \
                             ta2_xz_zz_z_0,  \
                             ta2_xz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzz_x_0[i] =
            ta2_xz_zz_0_0[i] * fe_0 - ta2_xz_zz_0_1[i] * fe_0 + ta1_z_zz_x_1[i] + ta2_xz_zz_x_0[i] * pa_x[i] - ta2_xz_zz_x_1[i] * pc_x[i];

        ta2_xz_xzz_y_0[i] = ta1_z_zz_y_1[i] + ta2_xz_zz_y_0[i] * pa_x[i] - ta2_xz_zz_y_1[i] * pc_x[i];

        ta2_xz_xzz_z_0[i] = ta1_z_zz_z_1[i] + ta2_xz_zz_z_0[i] * pa_x[i] - ta2_xz_zz_z_1[i] * pc_x[i];
    }

    // Set up 78-81 components of targeted buffer : FP

    auto ta2_xz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 78);

    auto ta2_xz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 79);

    auto ta2_xz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 80);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xz_y_x_0,   \
                             ta2_xz_y_x_1,   \
                             ta2_xz_y_y_0,   \
                             ta2_xz_y_y_1,   \
                             ta2_xz_y_z_0,   \
                             ta2_xz_y_z_1,   \
                             ta2_xz_yy_0_0,  \
                             ta2_xz_yy_0_1,  \
                             ta2_xz_yy_x_0,  \
                             ta2_xz_yy_x_1,  \
                             ta2_xz_yy_y_0,  \
                             ta2_xz_yy_y_1,  \
                             ta2_xz_yy_z_0,  \
                             ta2_xz_yy_z_1,  \
                             ta2_xz_yyy_x_0, \
                             ta2_xz_yyy_y_0, \
                             ta2_xz_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyy_x_0[i] = 2.0 * ta2_xz_y_x_0[i] * fe_0 - 2.0 * ta2_xz_y_x_1[i] * fe_0 + ta2_xz_yy_x_0[i] * pa_y[i] - ta2_xz_yy_x_1[i] * pc_y[i];

        ta2_xz_yyy_y_0[i] = 2.0 * ta2_xz_y_y_0[i] * fe_0 - 2.0 * ta2_xz_y_y_1[i] * fe_0 + ta2_xz_yy_0_0[i] * fe_0 - ta2_xz_yy_0_1[i] * fe_0 +
                            ta2_xz_yy_y_0[i] * pa_y[i] - ta2_xz_yy_y_1[i] * pc_y[i];

        ta2_xz_yyy_z_0[i] = 2.0 * ta2_xz_y_z_0[i] * fe_0 - 2.0 * ta2_xz_y_z_1[i] * fe_0 + ta2_xz_yy_z_0[i] * pa_y[i] - ta2_xz_yy_z_1[i] * pc_y[i];
    }

    // Set up 81-84 components of targeted buffer : FP

    auto ta2_xz_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 81);

    auto ta2_xz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 82);

    auto ta2_xz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 83);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_yy_x_1,   \
                             ta1_x_yy_y_1,   \
                             ta2_xz_yy_x_0,  \
                             ta2_xz_yy_x_1,  \
                             ta2_xz_yy_y_0,  \
                             ta2_xz_yy_y_1,  \
                             ta2_xz_yyz_x_0, \
                             ta2_xz_yyz_y_0, \
                             ta2_xz_yyz_z_0, \
                             ta2_xz_yz_z_0,  \
                             ta2_xz_yz_z_1,  \
                             ta2_xz_z_z_0,   \
                             ta2_xz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyz_x_0[i] = ta1_x_yy_x_1[i] + ta2_xz_yy_x_0[i] * pa_z[i] - ta2_xz_yy_x_1[i] * pc_z[i];

        ta2_xz_yyz_y_0[i] = ta1_x_yy_y_1[i] + ta2_xz_yy_y_0[i] * pa_z[i] - ta2_xz_yy_y_1[i] * pc_z[i];

        ta2_xz_yyz_z_0[i] = ta2_xz_z_z_0[i] * fe_0 - ta2_xz_z_z_1[i] * fe_0 + ta2_xz_yz_z_0[i] * pa_y[i] - ta2_xz_yz_z_1[i] * pc_y[i];
    }

    // Set up 84-87 components of targeted buffer : FP

    auto ta2_xz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 84);

    auto ta2_xz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 85);

    auto ta2_xz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 86);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xz_yzz_x_0, \
                             ta2_xz_yzz_y_0, \
                             ta2_xz_yzz_z_0, \
                             ta2_xz_zz_0_0,  \
                             ta2_xz_zz_0_1,  \
                             ta2_xz_zz_x_0,  \
                             ta2_xz_zz_x_1,  \
                             ta2_xz_zz_y_0,  \
                             ta2_xz_zz_y_1,  \
                             ta2_xz_zz_z_0,  \
                             ta2_xz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzz_x_0[i] = ta2_xz_zz_x_0[i] * pa_y[i] - ta2_xz_zz_x_1[i] * pc_y[i];

        ta2_xz_yzz_y_0[i] = ta2_xz_zz_0_0[i] * fe_0 - ta2_xz_zz_0_1[i] * fe_0 + ta2_xz_zz_y_0[i] * pa_y[i] - ta2_xz_zz_y_1[i] * pc_y[i];

        ta2_xz_yzz_z_0[i] = ta2_xz_zz_z_0[i] * pa_y[i] - ta2_xz_zz_z_1[i] * pc_y[i];
    }

    // Set up 87-90 components of targeted buffer : FP

    auto ta2_xz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 87);

    auto ta2_xz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 88);

    auto ta2_xz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 89);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_y_1,   \
                             ta1_x_zz_z_1,   \
                             ta2_xz_z_x_0,   \
                             ta2_xz_z_x_1,   \
                             ta2_xz_z_y_0,   \
                             ta2_xz_z_y_1,   \
                             ta2_xz_z_z_0,   \
                             ta2_xz_z_z_1,   \
                             ta2_xz_zz_0_0,  \
                             ta2_xz_zz_0_1,  \
                             ta2_xz_zz_x_0,  \
                             ta2_xz_zz_x_1,  \
                             ta2_xz_zz_y_0,  \
                             ta2_xz_zz_y_1,  \
                             ta2_xz_zz_z_0,  \
                             ta2_xz_zz_z_1,  \
                             ta2_xz_zzz_x_0, \
                             ta2_xz_zzz_y_0, \
                             ta2_xz_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzz_x_0[i] =
            2.0 * ta2_xz_z_x_0[i] * fe_0 - 2.0 * ta2_xz_z_x_1[i] * fe_0 + ta1_x_zz_x_1[i] + ta2_xz_zz_x_0[i] * pa_z[i] - ta2_xz_zz_x_1[i] * pc_z[i];

        ta2_xz_zzz_y_0[i] =
            2.0 * ta2_xz_z_y_0[i] * fe_0 - 2.0 * ta2_xz_z_y_1[i] * fe_0 + ta1_x_zz_y_1[i] + ta2_xz_zz_y_0[i] * pa_z[i] - ta2_xz_zz_y_1[i] * pc_z[i];

        ta2_xz_zzz_z_0[i] = 2.0 * ta2_xz_z_z_0[i] * fe_0 - 2.0 * ta2_xz_z_z_1[i] * fe_0 + ta2_xz_zz_0_0[i] * fe_0 - ta2_xz_zz_0_1[i] * fe_0 +
                            ta1_x_zz_z_1[i] + ta2_xz_zz_z_0[i] * pa_z[i] - ta2_xz_zz_z_1[i] * pc_z[i];
    }

    // Set up 90-93 components of targeted buffer : FP

    auto ta2_yy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 90);

    auto ta2_yy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 91);

    auto ta2_yy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 92);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yy_x_x_0,   \
                             ta2_yy_x_x_1,   \
                             ta2_yy_x_y_0,   \
                             ta2_yy_x_y_1,   \
                             ta2_yy_x_z_0,   \
                             ta2_yy_x_z_1,   \
                             ta2_yy_xx_0_0,  \
                             ta2_yy_xx_0_1,  \
                             ta2_yy_xx_x_0,  \
                             ta2_yy_xx_x_1,  \
                             ta2_yy_xx_y_0,  \
                             ta2_yy_xx_y_1,  \
                             ta2_yy_xx_z_0,  \
                             ta2_yy_xx_z_1,  \
                             ta2_yy_xxx_x_0, \
                             ta2_yy_xxx_y_0, \
                             ta2_yy_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxx_x_0[i] = 2.0 * ta2_yy_x_x_0[i] * fe_0 - 2.0 * ta2_yy_x_x_1[i] * fe_0 + ta2_yy_xx_0_0[i] * fe_0 - ta2_yy_xx_0_1[i] * fe_0 +
                            ta2_yy_xx_x_0[i] * pa_x[i] - ta2_yy_xx_x_1[i] * pc_x[i];

        ta2_yy_xxx_y_0[i] = 2.0 * ta2_yy_x_y_0[i] * fe_0 - 2.0 * ta2_yy_x_y_1[i] * fe_0 + ta2_yy_xx_y_0[i] * pa_x[i] - ta2_yy_xx_y_1[i] * pc_x[i];

        ta2_yy_xxx_z_0[i] = 2.0 * ta2_yy_x_z_0[i] * fe_0 - 2.0 * ta2_yy_x_z_1[i] * fe_0 + ta2_yy_xx_z_0[i] * pa_x[i] - ta2_yy_xx_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : FP

    auto ta2_yy_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 93);

    auto ta2_yy_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 94);

    auto ta2_yy_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 95);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_z_1,   \
                             ta2_yy_xx_x_0,  \
                             ta2_yy_xx_x_1,  \
                             ta2_yy_xx_z_0,  \
                             ta2_yy_xx_z_1,  \
                             ta2_yy_xxy_x_0, \
                             ta2_yy_xxy_y_0, \
                             ta2_yy_xxy_z_0, \
                             ta2_yy_xy_y_0,  \
                             ta2_yy_xy_y_1,  \
                             ta2_yy_y_y_0,   \
                             ta2_yy_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxy_x_0[i] = 2.0 * ta1_y_xx_x_1[i] + ta2_yy_xx_x_0[i] * pa_y[i] - ta2_yy_xx_x_1[i] * pc_y[i];

        ta2_yy_xxy_y_0[i] = ta2_yy_y_y_0[i] * fe_0 - ta2_yy_y_y_1[i] * fe_0 + ta2_yy_xy_y_0[i] * pa_x[i] - ta2_yy_xy_y_1[i] * pc_x[i];

        ta2_yy_xxy_z_0[i] = 2.0 * ta1_y_xx_z_1[i] + ta2_yy_xx_z_0[i] * pa_y[i] - ta2_yy_xx_z_1[i] * pc_y[i];
    }

    // Set up 96-99 components of targeted buffer : FP

    auto ta2_yy_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 96);

    auto ta2_yy_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 97);

    auto ta2_yy_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 98);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta2_yy_xx_x_0,  \
                             ta2_yy_xx_x_1,  \
                             ta2_yy_xx_y_0,  \
                             ta2_yy_xx_y_1,  \
                             ta2_yy_xxz_x_0, \
                             ta2_yy_xxz_y_0, \
                             ta2_yy_xxz_z_0, \
                             ta2_yy_xz_z_0,  \
                             ta2_yy_xz_z_1,  \
                             ta2_yy_z_z_0,   \
                             ta2_yy_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxz_x_0[i] = ta2_yy_xx_x_0[i] * pa_z[i] - ta2_yy_xx_x_1[i] * pc_z[i];

        ta2_yy_xxz_y_0[i] = ta2_yy_xx_y_0[i] * pa_z[i] - ta2_yy_xx_y_1[i] * pc_z[i];

        ta2_yy_xxz_z_0[i] = ta2_yy_z_z_0[i] * fe_0 - ta2_yy_z_z_1[i] * fe_0 + ta2_yy_xz_z_0[i] * pa_x[i] - ta2_yy_xz_z_1[i] * pc_x[i];
    }

    // Set up 99-102 components of targeted buffer : FP

    auto ta2_yy_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 99);

    auto ta2_yy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 100);

    auto ta2_yy_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 101);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yy_xyy_x_0, \
                             ta2_yy_xyy_y_0, \
                             ta2_yy_xyy_z_0, \
                             ta2_yy_yy_0_0,  \
                             ta2_yy_yy_0_1,  \
                             ta2_yy_yy_x_0,  \
                             ta2_yy_yy_x_1,  \
                             ta2_yy_yy_y_0,  \
                             ta2_yy_yy_y_1,  \
                             ta2_yy_yy_z_0,  \
                             ta2_yy_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyy_x_0[i] = ta2_yy_yy_0_0[i] * fe_0 - ta2_yy_yy_0_1[i] * fe_0 + ta2_yy_yy_x_0[i] * pa_x[i] - ta2_yy_yy_x_1[i] * pc_x[i];

        ta2_yy_xyy_y_0[i] = ta2_yy_yy_y_0[i] * pa_x[i] - ta2_yy_yy_y_1[i] * pc_x[i];

        ta2_yy_xyy_z_0[i] = ta2_yy_yy_z_0[i] * pa_x[i] - ta2_yy_yy_z_1[i] * pc_x[i];
    }

    // Set up 102-105 components of targeted buffer : FP

    auto ta2_yy_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 102);

    auto ta2_yy_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 103);

    auto ta2_yy_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 104);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta2_yy_xy_x_0,  \
                             ta2_yy_xy_x_1,  \
                             ta2_yy_xyz_x_0, \
                             ta2_yy_xyz_y_0, \
                             ta2_yy_xyz_z_0, \
                             ta2_yy_yz_y_0,  \
                             ta2_yy_yz_y_1,  \
                             ta2_yy_yz_z_0,  \
                             ta2_yy_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xyz_x_0[i] = ta2_yy_xy_x_0[i] * pa_z[i] - ta2_yy_xy_x_1[i] * pc_z[i];

        ta2_yy_xyz_y_0[i] = ta2_yy_yz_y_0[i] * pa_x[i] - ta2_yy_yz_y_1[i] * pc_x[i];

        ta2_yy_xyz_z_0[i] = ta2_yy_yz_z_0[i] * pa_x[i] - ta2_yy_yz_z_1[i] * pc_x[i];
    }

    // Set up 105-108 components of targeted buffer : FP

    auto ta2_yy_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 105);

    auto ta2_yy_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 106);

    auto ta2_yy_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 107);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yy_xzz_x_0, \
                             ta2_yy_xzz_y_0, \
                             ta2_yy_xzz_z_0, \
                             ta2_yy_zz_0_0,  \
                             ta2_yy_zz_0_1,  \
                             ta2_yy_zz_x_0,  \
                             ta2_yy_zz_x_1,  \
                             ta2_yy_zz_y_0,  \
                             ta2_yy_zz_y_1,  \
                             ta2_yy_zz_z_0,  \
                             ta2_yy_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzz_x_0[i] = ta2_yy_zz_0_0[i] * fe_0 - ta2_yy_zz_0_1[i] * fe_0 + ta2_yy_zz_x_0[i] * pa_x[i] - ta2_yy_zz_x_1[i] * pc_x[i];

        ta2_yy_xzz_y_0[i] = ta2_yy_zz_y_0[i] * pa_x[i] - ta2_yy_zz_y_1[i] * pc_x[i];

        ta2_yy_xzz_z_0[i] = ta2_yy_zz_z_0[i] * pa_x[i] - ta2_yy_zz_z_1[i] * pc_x[i];
    }

    // Set up 108-111 components of targeted buffer : FP

    auto ta2_yy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 108);

    auto ta2_yy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 109);

    auto ta2_yy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 110);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_z_1,   \
                             ta2_yy_y_x_0,   \
                             ta2_yy_y_x_1,   \
                             ta2_yy_y_y_0,   \
                             ta2_yy_y_y_1,   \
                             ta2_yy_y_z_0,   \
                             ta2_yy_y_z_1,   \
                             ta2_yy_yy_0_0,  \
                             ta2_yy_yy_0_1,  \
                             ta2_yy_yy_x_0,  \
                             ta2_yy_yy_x_1,  \
                             ta2_yy_yy_y_0,  \
                             ta2_yy_yy_y_1,  \
                             ta2_yy_yy_z_0,  \
                             ta2_yy_yy_z_1,  \
                             ta2_yy_yyy_x_0, \
                             ta2_yy_yyy_y_0, \
                             ta2_yy_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyy_x_0[i] = 2.0 * ta2_yy_y_x_0[i] * fe_0 - 2.0 * ta2_yy_y_x_1[i] * fe_0 + 2.0 * ta1_y_yy_x_1[i] + ta2_yy_yy_x_0[i] * pa_y[i] -
                            ta2_yy_yy_x_1[i] * pc_y[i];

        ta2_yy_yyy_y_0[i] = 2.0 * ta2_yy_y_y_0[i] * fe_0 - 2.0 * ta2_yy_y_y_1[i] * fe_0 + ta2_yy_yy_0_0[i] * fe_0 - ta2_yy_yy_0_1[i] * fe_0 +
                            2.0 * ta1_y_yy_y_1[i] + ta2_yy_yy_y_0[i] * pa_y[i] - ta2_yy_yy_y_1[i] * pc_y[i];

        ta2_yy_yyy_z_0[i] = 2.0 * ta2_yy_y_z_0[i] * fe_0 - 2.0 * ta2_yy_y_z_1[i] * fe_0 + 2.0 * ta1_y_yy_z_1[i] + ta2_yy_yy_z_0[i] * pa_y[i] -
                            ta2_yy_yy_z_1[i] * pc_y[i];
    }

    // Set up 111-114 components of targeted buffer : FP

    auto ta2_yy_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 111);

    auto ta2_yy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 112);

    auto ta2_yy_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 113);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_yy_yy_0_0,  \
                             ta2_yy_yy_0_1,  \
                             ta2_yy_yy_x_0,  \
                             ta2_yy_yy_x_1,  \
                             ta2_yy_yy_y_0,  \
                             ta2_yy_yy_y_1,  \
                             ta2_yy_yy_z_0,  \
                             ta2_yy_yy_z_1,  \
                             ta2_yy_yyz_x_0, \
                             ta2_yy_yyz_y_0, \
                             ta2_yy_yyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyz_x_0[i] = ta2_yy_yy_x_0[i] * pa_z[i] - ta2_yy_yy_x_1[i] * pc_z[i];

        ta2_yy_yyz_y_0[i] = ta2_yy_yy_y_0[i] * pa_z[i] - ta2_yy_yy_y_1[i] * pc_z[i];

        ta2_yy_yyz_z_0[i] = ta2_yy_yy_0_0[i] * fe_0 - ta2_yy_yy_0_1[i] * fe_0 + ta2_yy_yy_z_0[i] * pa_z[i] - ta2_yy_yy_z_1[i] * pc_z[i];
    }

    // Set up 114-117 components of targeted buffer : FP

    auto ta2_yy_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 114);

    auto ta2_yy_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 115);

    auto ta2_yy_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 116);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_zz_x_1,   \
                             ta1_y_zz_z_1,   \
                             ta2_yy_y_y_0,   \
                             ta2_yy_y_y_1,   \
                             ta2_yy_yz_y_0,  \
                             ta2_yy_yz_y_1,  \
                             ta2_yy_yzz_x_0, \
                             ta2_yy_yzz_y_0, \
                             ta2_yy_yzz_z_0, \
                             ta2_yy_zz_x_0,  \
                             ta2_yy_zz_x_1,  \
                             ta2_yy_zz_z_0,  \
                             ta2_yy_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzz_x_0[i] = 2.0 * ta1_y_zz_x_1[i] + ta2_yy_zz_x_0[i] * pa_y[i] - ta2_yy_zz_x_1[i] * pc_y[i];

        ta2_yy_yzz_y_0[i] = ta2_yy_y_y_0[i] * fe_0 - ta2_yy_y_y_1[i] * fe_0 + ta2_yy_yz_y_0[i] * pa_z[i] - ta2_yy_yz_y_1[i] * pc_z[i];

        ta2_yy_yzz_z_0[i] = 2.0 * ta1_y_zz_z_1[i] + ta2_yy_zz_z_0[i] * pa_y[i] - ta2_yy_zz_z_1[i] * pc_y[i];
    }

    // Set up 117-120 components of targeted buffer : FP

    auto ta2_yy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 117);

    auto ta2_yy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 118);

    auto ta2_yy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 119);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_yy_z_x_0,   \
                             ta2_yy_z_x_1,   \
                             ta2_yy_z_y_0,   \
                             ta2_yy_z_y_1,   \
                             ta2_yy_z_z_0,   \
                             ta2_yy_z_z_1,   \
                             ta2_yy_zz_0_0,  \
                             ta2_yy_zz_0_1,  \
                             ta2_yy_zz_x_0,  \
                             ta2_yy_zz_x_1,  \
                             ta2_yy_zz_y_0,  \
                             ta2_yy_zz_y_1,  \
                             ta2_yy_zz_z_0,  \
                             ta2_yy_zz_z_1,  \
                             ta2_yy_zzz_x_0, \
                             ta2_yy_zzz_y_0, \
                             ta2_yy_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzz_x_0[i] = 2.0 * ta2_yy_z_x_0[i] * fe_0 - 2.0 * ta2_yy_z_x_1[i] * fe_0 + ta2_yy_zz_x_0[i] * pa_z[i] - ta2_yy_zz_x_1[i] * pc_z[i];

        ta2_yy_zzz_y_0[i] = 2.0 * ta2_yy_z_y_0[i] * fe_0 - 2.0 * ta2_yy_z_y_1[i] * fe_0 + ta2_yy_zz_y_0[i] * pa_z[i] - ta2_yy_zz_y_1[i] * pc_z[i];

        ta2_yy_zzz_z_0[i] = 2.0 * ta2_yy_z_z_0[i] * fe_0 - 2.0 * ta2_yy_z_z_1[i] * fe_0 + ta2_yy_zz_0_0[i] * fe_0 - ta2_yy_zz_0_1[i] * fe_0 +
                            ta2_yy_zz_z_0[i] * pa_z[i] - ta2_yy_zz_z_1[i] * pc_z[i];
    }

    // Set up 120-123 components of targeted buffer : FP

    auto ta2_yz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 120);

    auto ta2_yz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 121);

    auto ta2_yz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 122);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yz_x_x_0,   \
                             ta2_yz_x_x_1,   \
                             ta2_yz_x_y_0,   \
                             ta2_yz_x_y_1,   \
                             ta2_yz_x_z_0,   \
                             ta2_yz_x_z_1,   \
                             ta2_yz_xx_0_0,  \
                             ta2_yz_xx_0_1,  \
                             ta2_yz_xx_x_0,  \
                             ta2_yz_xx_x_1,  \
                             ta2_yz_xx_y_0,  \
                             ta2_yz_xx_y_1,  \
                             ta2_yz_xx_z_0,  \
                             ta2_yz_xx_z_1,  \
                             ta2_yz_xxx_x_0, \
                             ta2_yz_xxx_y_0, \
                             ta2_yz_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxx_x_0[i] = 2.0 * ta2_yz_x_x_0[i] * fe_0 - 2.0 * ta2_yz_x_x_1[i] * fe_0 + ta2_yz_xx_0_0[i] * fe_0 - ta2_yz_xx_0_1[i] * fe_0 +
                            ta2_yz_xx_x_0[i] * pa_x[i] - ta2_yz_xx_x_1[i] * pc_x[i];

        ta2_yz_xxx_y_0[i] = 2.0 * ta2_yz_x_y_0[i] * fe_0 - 2.0 * ta2_yz_x_y_1[i] * fe_0 + ta2_yz_xx_y_0[i] * pa_x[i] - ta2_yz_xx_y_1[i] * pc_x[i];

        ta2_yz_xxx_z_0[i] = 2.0 * ta2_yz_x_z_0[i] * fe_0 - 2.0 * ta2_yz_x_z_1[i] * fe_0 + ta2_yz_xx_z_0[i] * pa_x[i] - ta2_yz_xx_z_1[i] * pc_x[i];
    }

    // Set up 123-126 components of targeted buffer : FP

    auto ta2_yz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 123);

    auto ta2_yz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 124);

    auto ta2_yz_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 125);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_z_1,   \
                             ta2_yz_xx_x_0,  \
                             ta2_yz_xx_x_1,  \
                             ta2_yz_xx_z_0,  \
                             ta2_yz_xx_z_1,  \
                             ta2_yz_xxy_x_0, \
                             ta2_yz_xxy_y_0, \
                             ta2_yz_xxy_z_0, \
                             ta2_yz_xy_y_0,  \
                             ta2_yz_xy_y_1,  \
                             ta2_yz_y_y_0,   \
                             ta2_yz_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxy_x_0[i] = ta1_z_xx_x_1[i] + ta2_yz_xx_x_0[i] * pa_y[i] - ta2_yz_xx_x_1[i] * pc_y[i];

        ta2_yz_xxy_y_0[i] = ta2_yz_y_y_0[i] * fe_0 - ta2_yz_y_y_1[i] * fe_0 + ta2_yz_xy_y_0[i] * pa_x[i] - ta2_yz_xy_y_1[i] * pc_x[i];

        ta2_yz_xxy_z_0[i] = ta1_z_xx_z_1[i] + ta2_yz_xx_z_0[i] * pa_y[i] - ta2_yz_xx_z_1[i] * pc_y[i];
    }

    // Set up 126-129 components of targeted buffer : FP

    auto ta2_yz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 126);

    auto ta2_yz_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 127);

    auto ta2_yz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 128);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_y_1,   \
                             ta2_yz_xx_x_0,  \
                             ta2_yz_xx_x_1,  \
                             ta2_yz_xx_y_0,  \
                             ta2_yz_xx_y_1,  \
                             ta2_yz_xxz_x_0, \
                             ta2_yz_xxz_y_0, \
                             ta2_yz_xxz_z_0, \
                             ta2_yz_xz_z_0,  \
                             ta2_yz_xz_z_1,  \
                             ta2_yz_z_z_0,   \
                             ta2_yz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxz_x_0[i] = ta1_y_xx_x_1[i] + ta2_yz_xx_x_0[i] * pa_z[i] - ta2_yz_xx_x_1[i] * pc_z[i];

        ta2_yz_xxz_y_0[i] = ta1_y_xx_y_1[i] + ta2_yz_xx_y_0[i] * pa_z[i] - ta2_yz_xx_y_1[i] * pc_z[i];

        ta2_yz_xxz_z_0[i] = ta2_yz_z_z_0[i] * fe_0 - ta2_yz_z_z_1[i] * fe_0 + ta2_yz_xz_z_0[i] * pa_x[i] - ta2_yz_xz_z_1[i] * pc_x[i];
    }

    // Set up 129-132 components of targeted buffer : FP

    auto ta2_yz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 129);

    auto ta2_yz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 130);

    auto ta2_yz_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 131);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yz_xyy_x_0, \
                             ta2_yz_xyy_y_0, \
                             ta2_yz_xyy_z_0, \
                             ta2_yz_yy_0_0,  \
                             ta2_yz_yy_0_1,  \
                             ta2_yz_yy_x_0,  \
                             ta2_yz_yy_x_1,  \
                             ta2_yz_yy_y_0,  \
                             ta2_yz_yy_y_1,  \
                             ta2_yz_yy_z_0,  \
                             ta2_yz_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyy_x_0[i] = ta2_yz_yy_0_0[i] * fe_0 - ta2_yz_yy_0_1[i] * fe_0 + ta2_yz_yy_x_0[i] * pa_x[i] - ta2_yz_yy_x_1[i] * pc_x[i];

        ta2_yz_xyy_y_0[i] = ta2_yz_yy_y_0[i] * pa_x[i] - ta2_yz_yy_y_1[i] * pc_x[i];

        ta2_yz_xyy_z_0[i] = ta2_yz_yy_z_0[i] * pa_x[i] - ta2_yz_yy_z_1[i] * pc_x[i];
    }

    // Set up 132-135 components of targeted buffer : FP

    auto ta2_yz_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 132);

    auto ta2_yz_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 133);

    auto ta2_yz_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 134);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xz_x_1,   \
                             ta2_yz_xyz_x_0, \
                             ta2_yz_xyz_y_0, \
                             ta2_yz_xyz_z_0, \
                             ta2_yz_xz_x_0,  \
                             ta2_yz_xz_x_1,  \
                             ta2_yz_yz_y_0,  \
                             ta2_yz_yz_y_1,  \
                             ta2_yz_yz_z_0,  \
                             ta2_yz_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xyz_x_0[i] = ta1_z_xz_x_1[i] + ta2_yz_xz_x_0[i] * pa_y[i] - ta2_yz_xz_x_1[i] * pc_y[i];

        ta2_yz_xyz_y_0[i] = ta2_yz_yz_y_0[i] * pa_x[i] - ta2_yz_yz_y_1[i] * pc_x[i];

        ta2_yz_xyz_z_0[i] = ta2_yz_yz_z_0[i] * pa_x[i] - ta2_yz_yz_z_1[i] * pc_x[i];
    }

    // Set up 135-138 components of targeted buffer : FP

    auto ta2_yz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 135);

    auto ta2_yz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 136);

    auto ta2_yz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 137);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yz_xzz_x_0, \
                             ta2_yz_xzz_y_0, \
                             ta2_yz_xzz_z_0, \
                             ta2_yz_zz_0_0,  \
                             ta2_yz_zz_0_1,  \
                             ta2_yz_zz_x_0,  \
                             ta2_yz_zz_x_1,  \
                             ta2_yz_zz_y_0,  \
                             ta2_yz_zz_y_1,  \
                             ta2_yz_zz_z_0,  \
                             ta2_yz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzz_x_0[i] = ta2_yz_zz_0_0[i] * fe_0 - ta2_yz_zz_0_1[i] * fe_0 + ta2_yz_zz_x_0[i] * pa_x[i] - ta2_yz_zz_x_1[i] * pc_x[i];

        ta2_yz_xzz_y_0[i] = ta2_yz_zz_y_0[i] * pa_x[i] - ta2_yz_zz_y_1[i] * pc_x[i];

        ta2_yz_xzz_z_0[i] = ta2_yz_zz_z_0[i] * pa_x[i] - ta2_yz_zz_z_1[i] * pc_x[i];
    }

    // Set up 138-141 components of targeted buffer : FP

    auto ta2_yz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 138);

    auto ta2_yz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 139);

    auto ta2_yz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 140);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_yy_x_1,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_z_1,   \
                             ta2_yz_y_x_0,   \
                             ta2_yz_y_x_1,   \
                             ta2_yz_y_y_0,   \
                             ta2_yz_y_y_1,   \
                             ta2_yz_y_z_0,   \
                             ta2_yz_y_z_1,   \
                             ta2_yz_yy_0_0,  \
                             ta2_yz_yy_0_1,  \
                             ta2_yz_yy_x_0,  \
                             ta2_yz_yy_x_1,  \
                             ta2_yz_yy_y_0,  \
                             ta2_yz_yy_y_1,  \
                             ta2_yz_yy_z_0,  \
                             ta2_yz_yy_z_1,  \
                             ta2_yz_yyy_x_0, \
                             ta2_yz_yyy_y_0, \
                             ta2_yz_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyy_x_0[i] =
            2.0 * ta2_yz_y_x_0[i] * fe_0 - 2.0 * ta2_yz_y_x_1[i] * fe_0 + ta1_z_yy_x_1[i] + ta2_yz_yy_x_0[i] * pa_y[i] - ta2_yz_yy_x_1[i] * pc_y[i];

        ta2_yz_yyy_y_0[i] = 2.0 * ta2_yz_y_y_0[i] * fe_0 - 2.0 * ta2_yz_y_y_1[i] * fe_0 + ta2_yz_yy_0_0[i] * fe_0 - ta2_yz_yy_0_1[i] * fe_0 +
                            ta1_z_yy_y_1[i] + ta2_yz_yy_y_0[i] * pa_y[i] - ta2_yz_yy_y_1[i] * pc_y[i];

        ta2_yz_yyy_z_0[i] =
            2.0 * ta2_yz_y_z_0[i] * fe_0 - 2.0 * ta2_yz_y_z_1[i] * fe_0 + ta1_z_yy_z_1[i] + ta2_yz_yy_z_0[i] * pa_y[i] - ta2_yz_yy_z_1[i] * pc_y[i];
    }

    // Set up 141-144 components of targeted buffer : FP

    auto ta2_yz_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 141);

    auto ta2_yz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 142);

    auto ta2_yz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 143);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_y_1,   \
                             ta1_z_yz_z_1,   \
                             ta2_yz_yy_x_0,  \
                             ta2_yz_yy_x_1,  \
                             ta2_yz_yy_y_0,  \
                             ta2_yz_yy_y_1,  \
                             ta2_yz_yyz_x_0, \
                             ta2_yz_yyz_y_0, \
                             ta2_yz_yyz_z_0, \
                             ta2_yz_yz_z_0,  \
                             ta2_yz_yz_z_1,  \
                             ta2_yz_z_z_0,   \
                             ta2_yz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyz_x_0[i] = ta1_y_yy_x_1[i] + ta2_yz_yy_x_0[i] * pa_z[i] - ta2_yz_yy_x_1[i] * pc_z[i];

        ta2_yz_yyz_y_0[i] = ta1_y_yy_y_1[i] + ta2_yz_yy_y_0[i] * pa_z[i] - ta2_yz_yy_y_1[i] * pc_z[i];

        ta2_yz_yyz_z_0[i] =
            ta2_yz_z_z_0[i] * fe_0 - ta2_yz_z_z_1[i] * fe_0 + ta1_z_yz_z_1[i] + ta2_yz_yz_z_0[i] * pa_y[i] - ta2_yz_yz_z_1[i] * pc_y[i];
    }

    // Set up 144-147 components of targeted buffer : FP

    auto ta2_yz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 144);

    auto ta2_yz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 145);

    auto ta2_yz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 146);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_z_1,   \
                             ta2_yz_yzz_x_0, \
                             ta2_yz_yzz_y_0, \
                             ta2_yz_yzz_z_0, \
                             ta2_yz_zz_0_0,  \
                             ta2_yz_zz_0_1,  \
                             ta2_yz_zz_x_0,  \
                             ta2_yz_zz_x_1,  \
                             ta2_yz_zz_y_0,  \
                             ta2_yz_zz_y_1,  \
                             ta2_yz_zz_z_0,  \
                             ta2_yz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzz_x_0[i] = ta1_z_zz_x_1[i] + ta2_yz_zz_x_0[i] * pa_y[i] - ta2_yz_zz_x_1[i] * pc_y[i];

        ta2_yz_yzz_y_0[i] =
            ta2_yz_zz_0_0[i] * fe_0 - ta2_yz_zz_0_1[i] * fe_0 + ta1_z_zz_y_1[i] + ta2_yz_zz_y_0[i] * pa_y[i] - ta2_yz_zz_y_1[i] * pc_y[i];

        ta2_yz_yzz_z_0[i] = ta1_z_zz_z_1[i] + ta2_yz_zz_z_0[i] * pa_y[i] - ta2_yz_zz_z_1[i] * pc_y[i];
    }

    // Set up 147-150 components of targeted buffer : FP

    auto ta2_yz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 147);

    auto ta2_yz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 148);

    auto ta2_yz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 149);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_zz_x_1,   \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_z_1,   \
                             ta2_yz_z_x_0,   \
                             ta2_yz_z_x_1,   \
                             ta2_yz_z_y_0,   \
                             ta2_yz_z_y_1,   \
                             ta2_yz_z_z_0,   \
                             ta2_yz_z_z_1,   \
                             ta2_yz_zz_0_0,  \
                             ta2_yz_zz_0_1,  \
                             ta2_yz_zz_x_0,  \
                             ta2_yz_zz_x_1,  \
                             ta2_yz_zz_y_0,  \
                             ta2_yz_zz_y_1,  \
                             ta2_yz_zz_z_0,  \
                             ta2_yz_zz_z_1,  \
                             ta2_yz_zzz_x_0, \
                             ta2_yz_zzz_y_0, \
                             ta2_yz_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzz_x_0[i] =
            2.0 * ta2_yz_z_x_0[i] * fe_0 - 2.0 * ta2_yz_z_x_1[i] * fe_0 + ta1_y_zz_x_1[i] + ta2_yz_zz_x_0[i] * pa_z[i] - ta2_yz_zz_x_1[i] * pc_z[i];

        ta2_yz_zzz_y_0[i] =
            2.0 * ta2_yz_z_y_0[i] * fe_0 - 2.0 * ta2_yz_z_y_1[i] * fe_0 + ta1_y_zz_y_1[i] + ta2_yz_zz_y_0[i] * pa_z[i] - ta2_yz_zz_y_1[i] * pc_z[i];

        ta2_yz_zzz_z_0[i] = 2.0 * ta2_yz_z_z_0[i] * fe_0 - 2.0 * ta2_yz_z_z_1[i] * fe_0 + ta2_yz_zz_0_0[i] * fe_0 - ta2_yz_zz_0_1[i] * fe_0 +
                            ta1_y_zz_z_1[i] + ta2_yz_zz_z_0[i] * pa_z[i] - ta2_yz_zz_z_1[i] * pc_z[i];
    }

    // Set up 150-153 components of targeted buffer : FP

    auto ta2_zz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 150);

    auto ta2_zz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 151);

    auto ta2_zz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 152);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_zz_x_x_0,   \
                             ta2_zz_x_x_1,   \
                             ta2_zz_x_y_0,   \
                             ta2_zz_x_y_1,   \
                             ta2_zz_x_z_0,   \
                             ta2_zz_x_z_1,   \
                             ta2_zz_xx_0_0,  \
                             ta2_zz_xx_0_1,  \
                             ta2_zz_xx_x_0,  \
                             ta2_zz_xx_x_1,  \
                             ta2_zz_xx_y_0,  \
                             ta2_zz_xx_y_1,  \
                             ta2_zz_xx_z_0,  \
                             ta2_zz_xx_z_1,  \
                             ta2_zz_xxx_x_0, \
                             ta2_zz_xxx_y_0, \
                             ta2_zz_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxx_x_0[i] = 2.0 * ta2_zz_x_x_0[i] * fe_0 - 2.0 * ta2_zz_x_x_1[i] * fe_0 + ta2_zz_xx_0_0[i] * fe_0 - ta2_zz_xx_0_1[i] * fe_0 +
                            ta2_zz_xx_x_0[i] * pa_x[i] - ta2_zz_xx_x_1[i] * pc_x[i];

        ta2_zz_xxx_y_0[i] = 2.0 * ta2_zz_x_y_0[i] * fe_0 - 2.0 * ta2_zz_x_y_1[i] * fe_0 + ta2_zz_xx_y_0[i] * pa_x[i] - ta2_zz_xx_y_1[i] * pc_x[i];

        ta2_zz_xxx_z_0[i] = 2.0 * ta2_zz_x_z_0[i] * fe_0 - 2.0 * ta2_zz_x_z_1[i] * fe_0 + ta2_zz_xx_z_0[i] * pa_x[i] - ta2_zz_xx_z_1[i] * pc_x[i];
    }

    // Set up 153-156 components of targeted buffer : FP

    auto ta2_zz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 153);

    auto ta2_zz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 154);

    auto ta2_zz_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 155);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta2_zz_xx_x_0,  \
                             ta2_zz_xx_x_1,  \
                             ta2_zz_xx_z_0,  \
                             ta2_zz_xx_z_1,  \
                             ta2_zz_xxy_x_0, \
                             ta2_zz_xxy_y_0, \
                             ta2_zz_xxy_z_0, \
                             ta2_zz_xy_y_0,  \
                             ta2_zz_xy_y_1,  \
                             ta2_zz_y_y_0,   \
                             ta2_zz_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxy_x_0[i] = ta2_zz_xx_x_0[i] * pa_y[i] - ta2_zz_xx_x_1[i] * pc_y[i];

        ta2_zz_xxy_y_0[i] = ta2_zz_y_y_0[i] * fe_0 - ta2_zz_y_y_1[i] * fe_0 + ta2_zz_xy_y_0[i] * pa_x[i] - ta2_zz_xy_y_1[i] * pc_x[i];

        ta2_zz_xxy_z_0[i] = ta2_zz_xx_z_0[i] * pa_y[i] - ta2_zz_xx_z_1[i] * pc_y[i];
    }

    // Set up 156-159 components of targeted buffer : FP

    auto ta2_zz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 156);

    auto ta2_zz_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 157);

    auto ta2_zz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 158);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_y_1,   \
                             ta2_zz_xx_x_0,  \
                             ta2_zz_xx_x_1,  \
                             ta2_zz_xx_y_0,  \
                             ta2_zz_xx_y_1,  \
                             ta2_zz_xxz_x_0, \
                             ta2_zz_xxz_y_0, \
                             ta2_zz_xxz_z_0, \
                             ta2_zz_xz_z_0,  \
                             ta2_zz_xz_z_1,  \
                             ta2_zz_z_z_0,   \
                             ta2_zz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxz_x_0[i] = 2.0 * ta1_z_xx_x_1[i] + ta2_zz_xx_x_0[i] * pa_z[i] - ta2_zz_xx_x_1[i] * pc_z[i];

        ta2_zz_xxz_y_0[i] = 2.0 * ta1_z_xx_y_1[i] + ta2_zz_xx_y_0[i] * pa_z[i] - ta2_zz_xx_y_1[i] * pc_z[i];

        ta2_zz_xxz_z_0[i] = ta2_zz_z_z_0[i] * fe_0 - ta2_zz_z_z_1[i] * fe_0 + ta2_zz_xz_z_0[i] * pa_x[i] - ta2_zz_xz_z_1[i] * pc_x[i];
    }

    // Set up 159-162 components of targeted buffer : FP

    auto ta2_zz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 159);

    auto ta2_zz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 160);

    auto ta2_zz_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 161);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_zz_xyy_x_0, \
                             ta2_zz_xyy_y_0, \
                             ta2_zz_xyy_z_0, \
                             ta2_zz_yy_0_0,  \
                             ta2_zz_yy_0_1,  \
                             ta2_zz_yy_x_0,  \
                             ta2_zz_yy_x_1,  \
                             ta2_zz_yy_y_0,  \
                             ta2_zz_yy_y_1,  \
                             ta2_zz_yy_z_0,  \
                             ta2_zz_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyy_x_0[i] = ta2_zz_yy_0_0[i] * fe_0 - ta2_zz_yy_0_1[i] * fe_0 + ta2_zz_yy_x_0[i] * pa_x[i] - ta2_zz_yy_x_1[i] * pc_x[i];

        ta2_zz_xyy_y_0[i] = ta2_zz_yy_y_0[i] * pa_x[i] - ta2_zz_yy_y_1[i] * pc_x[i];

        ta2_zz_xyy_z_0[i] = ta2_zz_yy_z_0[i] * pa_x[i] - ta2_zz_yy_z_1[i] * pc_x[i];
    }

    // Set up 162-165 components of targeted buffer : FP

    auto ta2_zz_xyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 162);

    auto ta2_zz_xyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 163);

    auto ta2_zz_xyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 164);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta2_zz_xyz_x_0, \
                             ta2_zz_xyz_y_0, \
                             ta2_zz_xyz_z_0, \
                             ta2_zz_xz_x_0,  \
                             ta2_zz_xz_x_1,  \
                             ta2_zz_yz_y_0,  \
                             ta2_zz_yz_y_1,  \
                             ta2_zz_yz_z_0,  \
                             ta2_zz_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xyz_x_0[i] = ta2_zz_xz_x_0[i] * pa_y[i] - ta2_zz_xz_x_1[i] * pc_y[i];

        ta2_zz_xyz_y_0[i] = ta2_zz_yz_y_0[i] * pa_x[i] - ta2_zz_yz_y_1[i] * pc_x[i];

        ta2_zz_xyz_z_0[i] = ta2_zz_yz_z_0[i] * pa_x[i] - ta2_zz_yz_z_1[i] * pc_x[i];
    }

    // Set up 165-168 components of targeted buffer : FP

    auto ta2_zz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 165);

    auto ta2_zz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 166);

    auto ta2_zz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 167);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_zz_xzz_x_0, \
                             ta2_zz_xzz_y_0, \
                             ta2_zz_xzz_z_0, \
                             ta2_zz_zz_0_0,  \
                             ta2_zz_zz_0_1,  \
                             ta2_zz_zz_x_0,  \
                             ta2_zz_zz_x_1,  \
                             ta2_zz_zz_y_0,  \
                             ta2_zz_zz_y_1,  \
                             ta2_zz_zz_z_0,  \
                             ta2_zz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzz_x_0[i] = ta2_zz_zz_0_0[i] * fe_0 - ta2_zz_zz_0_1[i] * fe_0 + ta2_zz_zz_x_0[i] * pa_x[i] - ta2_zz_zz_x_1[i] * pc_x[i];

        ta2_zz_xzz_y_0[i] = ta2_zz_zz_y_0[i] * pa_x[i] - ta2_zz_zz_y_1[i] * pc_x[i];

        ta2_zz_xzz_z_0[i] = ta2_zz_zz_z_0[i] * pa_x[i] - ta2_zz_zz_z_1[i] * pc_x[i];
    }

    // Set up 168-171 components of targeted buffer : FP

    auto ta2_zz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 168);

    auto ta2_zz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 169);

    auto ta2_zz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 170);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_zz_y_x_0,   \
                             ta2_zz_y_x_1,   \
                             ta2_zz_y_y_0,   \
                             ta2_zz_y_y_1,   \
                             ta2_zz_y_z_0,   \
                             ta2_zz_y_z_1,   \
                             ta2_zz_yy_0_0,  \
                             ta2_zz_yy_0_1,  \
                             ta2_zz_yy_x_0,  \
                             ta2_zz_yy_x_1,  \
                             ta2_zz_yy_y_0,  \
                             ta2_zz_yy_y_1,  \
                             ta2_zz_yy_z_0,  \
                             ta2_zz_yy_z_1,  \
                             ta2_zz_yyy_x_0, \
                             ta2_zz_yyy_y_0, \
                             ta2_zz_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyy_x_0[i] = 2.0 * ta2_zz_y_x_0[i] * fe_0 - 2.0 * ta2_zz_y_x_1[i] * fe_0 + ta2_zz_yy_x_0[i] * pa_y[i] - ta2_zz_yy_x_1[i] * pc_y[i];

        ta2_zz_yyy_y_0[i] = 2.0 * ta2_zz_y_y_0[i] * fe_0 - 2.0 * ta2_zz_y_y_1[i] * fe_0 + ta2_zz_yy_0_0[i] * fe_0 - ta2_zz_yy_0_1[i] * fe_0 +
                            ta2_zz_yy_y_0[i] * pa_y[i] - ta2_zz_yy_y_1[i] * pc_y[i];

        ta2_zz_yyy_z_0[i] = 2.0 * ta2_zz_y_z_0[i] * fe_0 - 2.0 * ta2_zz_y_z_1[i] * fe_0 + ta2_zz_yy_z_0[i] * pa_y[i] - ta2_zz_yy_z_1[i] * pc_y[i];
    }

    // Set up 171-174 components of targeted buffer : FP

    auto ta2_zz_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 171);

    auto ta2_zz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 172);

    auto ta2_zz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 173);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_yy_x_1,   \
                             ta1_z_yy_y_1,   \
                             ta2_zz_yy_x_0,  \
                             ta2_zz_yy_x_1,  \
                             ta2_zz_yy_y_0,  \
                             ta2_zz_yy_y_1,  \
                             ta2_zz_yyz_x_0, \
                             ta2_zz_yyz_y_0, \
                             ta2_zz_yyz_z_0, \
                             ta2_zz_yz_z_0,  \
                             ta2_zz_yz_z_1,  \
                             ta2_zz_z_z_0,   \
                             ta2_zz_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyz_x_0[i] = 2.0 * ta1_z_yy_x_1[i] + ta2_zz_yy_x_0[i] * pa_z[i] - ta2_zz_yy_x_1[i] * pc_z[i];

        ta2_zz_yyz_y_0[i] = 2.0 * ta1_z_yy_y_1[i] + ta2_zz_yy_y_0[i] * pa_z[i] - ta2_zz_yy_y_1[i] * pc_z[i];

        ta2_zz_yyz_z_0[i] = ta2_zz_z_z_0[i] * fe_0 - ta2_zz_z_z_1[i] * fe_0 + ta2_zz_yz_z_0[i] * pa_y[i] - ta2_zz_yz_z_1[i] * pc_y[i];
    }

    // Set up 174-177 components of targeted buffer : FP

    auto ta2_zz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 174);

    auto ta2_zz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 175);

    auto ta2_zz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 176);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_zz_yzz_x_0, \
                             ta2_zz_yzz_y_0, \
                             ta2_zz_yzz_z_0, \
                             ta2_zz_zz_0_0,  \
                             ta2_zz_zz_0_1,  \
                             ta2_zz_zz_x_0,  \
                             ta2_zz_zz_x_1,  \
                             ta2_zz_zz_y_0,  \
                             ta2_zz_zz_y_1,  \
                             ta2_zz_zz_z_0,  \
                             ta2_zz_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzz_x_0[i] = ta2_zz_zz_x_0[i] * pa_y[i] - ta2_zz_zz_x_1[i] * pc_y[i];

        ta2_zz_yzz_y_0[i] = ta2_zz_zz_0_0[i] * fe_0 - ta2_zz_zz_0_1[i] * fe_0 + ta2_zz_zz_y_0[i] * pa_y[i] - ta2_zz_zz_y_1[i] * pc_y[i];

        ta2_zz_yzz_z_0[i] = ta2_zz_zz_z_0[i] * pa_y[i] - ta2_zz_zz_z_1[i] * pc_y[i];
    }

    // Set up 177-180 components of targeted buffer : FP

    auto ta2_zz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 177);

    auto ta2_zz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 178);

    auto ta2_zz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 179);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_z_1,   \
                             ta2_zz_z_x_0,   \
                             ta2_zz_z_x_1,   \
                             ta2_zz_z_y_0,   \
                             ta2_zz_z_y_1,   \
                             ta2_zz_z_z_0,   \
                             ta2_zz_z_z_1,   \
                             ta2_zz_zz_0_0,  \
                             ta2_zz_zz_0_1,  \
                             ta2_zz_zz_x_0,  \
                             ta2_zz_zz_x_1,  \
                             ta2_zz_zz_y_0,  \
                             ta2_zz_zz_y_1,  \
                             ta2_zz_zz_z_0,  \
                             ta2_zz_zz_z_1,  \
                             ta2_zz_zzz_x_0, \
                             ta2_zz_zzz_y_0, \
                             ta2_zz_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzz_x_0[i] = 2.0 * ta2_zz_z_x_0[i] * fe_0 - 2.0 * ta2_zz_z_x_1[i] * fe_0 + 2.0 * ta1_z_zz_x_1[i] + ta2_zz_zz_x_0[i] * pa_z[i] -
                            ta2_zz_zz_x_1[i] * pc_z[i];

        ta2_zz_zzz_y_0[i] = 2.0 * ta2_zz_z_y_0[i] * fe_0 - 2.0 * ta2_zz_z_y_1[i] * fe_0 + 2.0 * ta1_z_zz_y_1[i] + ta2_zz_zz_y_0[i] * pa_z[i] -
                            ta2_zz_zz_y_1[i] * pc_z[i];

        ta2_zz_zzz_z_0[i] = 2.0 * ta2_zz_z_z_0[i] * fe_0 - 2.0 * ta2_zz_z_z_1[i] * fe_0 + ta2_zz_zz_0_0[i] * fe_0 - ta2_zz_zz_0_1[i] * fe_0 +
                            2.0 * ta1_z_zz_z_1[i] + ta2_zz_zz_z_0[i] * pa_z[i] - ta2_zz_zz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
