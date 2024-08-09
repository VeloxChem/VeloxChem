#include "NuclearPotentialGeom020PrimRecDP.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_dp(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_dp,
                                        const size_t idx_npot_geom_020_0_sp,
                                        const size_t idx_npot_geom_020_1_sp,
                                        const size_t idx_npot_geom_020_0_ps,
                                        const size_t idx_npot_geom_020_1_ps,
                                        const size_t idx_npot_geom_010_1_pp,
                                        const size_t idx_npot_geom_020_0_pp,
                                        const size_t idx_npot_geom_020_1_pp,
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

    // Set up components of auxiliary buffer : SP

    auto ta2_xx_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp);

    auto ta2_xx_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 1);

    auto ta2_xx_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 2);

    auto ta2_xy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 3);

    auto ta2_xy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 4);

    auto ta2_xy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 5);

    auto ta2_xz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 6);

    auto ta2_xz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 7);

    auto ta2_xz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 8);

    auto ta2_yy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 9);

    auto ta2_yy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 10);

    auto ta2_yy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 11);

    auto ta2_yz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 12);

    auto ta2_yz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 13);

    auto ta2_yz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 14);

    auto ta2_zz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 15);

    auto ta2_zz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 16);

    auto ta2_zz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 17);

    // Set up components of auxiliary buffer : SP

    auto ta2_xx_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp);

    auto ta2_xx_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 1);

    auto ta2_xx_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 2);

    auto ta2_xy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 3);

    auto ta2_xy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 4);

    auto ta2_xy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 5);

    auto ta2_xz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 6);

    auto ta2_xz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 7);

    auto ta2_xz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 8);

    auto ta2_yy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 9);

    auto ta2_yy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 10);

    auto ta2_yy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 11);

    auto ta2_yz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 12);

    auto ta2_yz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 13);

    auto ta2_yz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 14);

    auto ta2_zz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 15);

    auto ta2_zz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 16);

    auto ta2_zz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 17);

    // Set up components of auxiliary buffer : PS

    auto ta2_xx_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps);

    auto ta2_xx_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 1);

    auto ta2_xx_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 2);

    auto ta2_xy_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 3);

    auto ta2_xy_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 4);

    auto ta2_xy_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 5);

    auto ta2_xz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 6);

    auto ta2_xz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 7);

    auto ta2_xz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 8);

    auto ta2_yy_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 9);

    auto ta2_yy_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 10);

    auto ta2_yy_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 11);

    auto ta2_yz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 12);

    auto ta2_yz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 13);

    auto ta2_yz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 14);

    auto ta2_zz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 15);

    auto ta2_zz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 16);

    auto ta2_zz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 17);

    // Set up components of auxiliary buffer : PS

    auto ta2_xx_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps);

    auto ta2_xx_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 1);

    auto ta2_xx_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 2);

    auto ta2_xy_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 3);

    auto ta2_xy_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 4);

    auto ta2_xy_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 5);

    auto ta2_xz_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 6);

    auto ta2_xz_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 7);

    auto ta2_xz_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 8);

    auto ta2_yy_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 9);

    auto ta2_yy_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 10);

    auto ta2_yy_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 11);

    auto ta2_yz_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 12);

    auto ta2_yz_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 13);

    auto ta2_yz_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 14);

    auto ta2_zz_x_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 15);

    auto ta2_zz_y_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 16);

    auto ta2_zz_z_0_1 = pbuffer.data(idx_npot_geom_020_1_ps + 17);

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

    // Set up 0-3 components of targeted buffer : DP

    auto ta2_xx_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp);

    auto ta2_xx_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 1);

    auto ta2_xx_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 2);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_x_x_1, ta1_x_x_y_1, ta1_x_x_z_1, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_x_0_0, ta2_xx_x_0_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_y_0, ta2_xx_x_y_1, ta2_xx_x_z_0, ta2_xx_x_z_1, ta2_xx_xx_x_0, ta2_xx_xx_y_0, ta2_xx_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xx_x_0[i] = ta2_xx_0_x_0[i] * fe_0 - ta2_xx_0_x_1[i] * fe_0 + ta2_xx_x_0_0[i] * fe_0 - ta2_xx_x_0_1[i] * fe_0 + 2.0 * ta1_x_x_x_1[i] + ta2_xx_x_x_0[i] * pa_x[i] - ta2_xx_x_x_1[i] * pc_x[i];

        ta2_xx_xx_y_0[i] = ta2_xx_0_y_0[i] * fe_0 - ta2_xx_0_y_1[i] * fe_0 + 2.0 * ta1_x_x_y_1[i] + ta2_xx_x_y_0[i] * pa_x[i] - ta2_xx_x_y_1[i] * pc_x[i];

        ta2_xx_xx_z_0[i] = ta2_xx_0_z_0[i] * fe_0 - ta2_xx_0_z_1[i] * fe_0 + 2.0 * ta1_x_x_z_1[i] + ta2_xx_x_z_0[i] * pa_x[i] - ta2_xx_x_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto ta2_xx_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 3);

    auto ta2_xx_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 4);

    auto ta2_xx_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 5);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_y_y_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_z_0, ta2_xx_x_z_1, ta2_xx_xy_x_0, ta2_xx_xy_y_0, ta2_xx_xy_z_0, ta2_xx_y_y_0, ta2_xx_y_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xy_x_0[i] = ta2_xx_x_x_0[i] * pa_y[i] - ta2_xx_x_x_1[i] * pc_y[i];

        ta2_xx_xy_y_0[i] = 2.0 * ta1_x_y_y_1[i] + ta2_xx_y_y_0[i] * pa_x[i] - ta2_xx_y_y_1[i] * pc_x[i];

        ta2_xx_xy_z_0[i] = ta2_xx_x_z_0[i] * pa_y[i] - ta2_xx_x_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto ta2_xx_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 6);

    auto ta2_xx_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 7);

    auto ta2_xx_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 8);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_z_z_1, ta2_xx_x_x_0, ta2_xx_x_x_1, ta2_xx_x_y_0, ta2_xx_x_y_1, ta2_xx_xz_x_0, ta2_xx_xz_y_0, ta2_xx_xz_z_0, ta2_xx_z_z_0, ta2_xx_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xz_x_0[i] = ta2_xx_x_x_0[i] * pa_z[i] - ta2_xx_x_x_1[i] * pc_z[i];

        ta2_xx_xz_y_0[i] = ta2_xx_x_y_0[i] * pa_z[i] - ta2_xx_x_y_1[i] * pc_z[i];

        ta2_xx_xz_z_0[i] = 2.0 * ta1_x_z_z_1[i] + ta2_xx_z_z_0[i] * pa_x[i] - ta2_xx_z_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto ta2_xx_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 9);

    auto ta2_xx_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 10);

    auto ta2_xx_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 11);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_y_0_0, ta2_xx_y_0_1, ta2_xx_y_x_0, ta2_xx_y_x_1, ta2_xx_y_y_0, ta2_xx_y_y_1, ta2_xx_y_z_0, ta2_xx_y_z_1, ta2_xx_yy_x_0, ta2_xx_yy_y_0, ta2_xx_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yy_x_0[i] = ta2_xx_0_x_0[i] * fe_0 - ta2_xx_0_x_1[i] * fe_0 + ta2_xx_y_x_0[i] * pa_y[i] - ta2_xx_y_x_1[i] * pc_y[i];

        ta2_xx_yy_y_0[i] = ta2_xx_0_y_0[i] * fe_0 - ta2_xx_0_y_1[i] * fe_0 + ta2_xx_y_0_0[i] * fe_0 - ta2_xx_y_0_1[i] * fe_0 + ta2_xx_y_y_0[i] * pa_y[i] - ta2_xx_y_y_1[i] * pc_y[i];

        ta2_xx_yy_z_0[i] = ta2_xx_0_z_0[i] * fe_0 - ta2_xx_0_z_1[i] * fe_0 + ta2_xx_y_z_0[i] * pa_y[i] - ta2_xx_y_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto ta2_xx_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 12);

    auto ta2_xx_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 13);

    auto ta2_xx_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 14);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_y_y_0, ta2_xx_y_y_1, ta2_xx_yz_x_0, ta2_xx_yz_y_0, ta2_xx_yz_z_0, ta2_xx_z_x_0, ta2_xx_z_x_1, ta2_xx_z_z_0, ta2_xx_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_yz_x_0[i] = ta2_xx_z_x_0[i] * pa_y[i] - ta2_xx_z_x_1[i] * pc_y[i];

        ta2_xx_yz_y_0[i] = ta2_xx_y_y_0[i] * pa_z[i] - ta2_xx_y_y_1[i] * pc_z[i];

        ta2_xx_yz_z_0[i] = ta2_xx_z_z_0[i] * pa_y[i] - ta2_xx_z_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto ta2_xx_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 15);

    auto ta2_xx_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 16);

    auto ta2_xx_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 17);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_z_0_0, ta2_xx_z_0_1, ta2_xx_z_x_0, ta2_xx_z_x_1, ta2_xx_z_y_0, ta2_xx_z_y_1, ta2_xx_z_z_0, ta2_xx_z_z_1, ta2_xx_zz_x_0, ta2_xx_zz_y_0, ta2_xx_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zz_x_0[i] = ta2_xx_0_x_0[i] * fe_0 - ta2_xx_0_x_1[i] * fe_0 + ta2_xx_z_x_0[i] * pa_z[i] - ta2_xx_z_x_1[i] * pc_z[i];

        ta2_xx_zz_y_0[i] = ta2_xx_0_y_0[i] * fe_0 - ta2_xx_0_y_1[i] * fe_0 + ta2_xx_z_y_0[i] * pa_z[i] - ta2_xx_z_y_1[i] * pc_z[i];

        ta2_xx_zz_z_0[i] = ta2_xx_0_z_0[i] * fe_0 - ta2_xx_0_z_1[i] * fe_0 + ta2_xx_z_0_0[i] * fe_0 - ta2_xx_z_0_1[i] * fe_0 + ta2_xx_z_z_0[i] * pa_z[i] - ta2_xx_z_z_1[i] * pc_z[i];
    }

    // Set up 18-21 components of targeted buffer : DP

    auto ta2_xy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 18);

    auto ta2_xy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 19);

    auto ta2_xy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 20);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_x_x_1, ta1_y_x_y_1, ta1_y_x_z_1, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_x_0_0, ta2_xy_x_0_1, ta2_xy_x_x_0, ta2_xy_x_x_1, ta2_xy_x_y_0, ta2_xy_x_y_1, ta2_xy_x_z_0, ta2_xy_x_z_1, ta2_xy_xx_x_0, ta2_xy_xx_y_0, ta2_xy_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xx_x_0[i] = ta2_xy_0_x_0[i] * fe_0 - ta2_xy_0_x_1[i] * fe_0 + ta2_xy_x_0_0[i] * fe_0 - ta2_xy_x_0_1[i] * fe_0 + ta1_y_x_x_1[i] + ta2_xy_x_x_0[i] * pa_x[i] - ta2_xy_x_x_1[i] * pc_x[i];

        ta2_xy_xx_y_0[i] = ta2_xy_0_y_0[i] * fe_0 - ta2_xy_0_y_1[i] * fe_0 + ta1_y_x_y_1[i] + ta2_xy_x_y_0[i] * pa_x[i] - ta2_xy_x_y_1[i] * pc_x[i];

        ta2_xy_xx_z_0[i] = ta2_xy_0_z_0[i] * fe_0 - ta2_xy_0_z_1[i] * fe_0 + ta1_y_x_z_1[i] + ta2_xy_x_z_0[i] * pa_x[i] - ta2_xy_x_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : DP

    auto ta2_xy_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 21);

    auto ta2_xy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 22);

    auto ta2_xy_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 23);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_x_x_1, ta1_y_y_y_1, ta1_y_y_z_1, ta2_xy_x_x_0, ta2_xy_x_x_1, ta2_xy_xy_x_0, ta2_xy_xy_y_0, ta2_xy_xy_z_0, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_y_z_0, ta2_xy_y_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xy_x_0[i] = ta1_x_x_x_1[i] + ta2_xy_x_x_0[i] * pa_y[i] - ta2_xy_x_x_1[i] * pc_y[i];

        ta2_xy_xy_y_0[i] = ta1_y_y_y_1[i] + ta2_xy_y_y_0[i] * pa_x[i] - ta2_xy_y_y_1[i] * pc_x[i];

        ta2_xy_xy_z_0[i] = ta1_y_y_z_1[i] + ta2_xy_y_z_0[i] * pa_x[i] - ta2_xy_y_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : DP

    auto ta2_xy_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 24);

    auto ta2_xy_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 25);

    auto ta2_xy_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 26);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_z_z_1, ta2_xy_x_x_0, ta2_xy_x_x_1, ta2_xy_x_y_0, ta2_xy_x_y_1, ta2_xy_xz_x_0, ta2_xy_xz_y_0, ta2_xy_xz_z_0, ta2_xy_z_z_0, ta2_xy_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xz_x_0[i] = ta2_xy_x_x_0[i] * pa_z[i] - ta2_xy_x_x_1[i] * pc_z[i];

        ta2_xy_xz_y_0[i] = ta2_xy_x_y_0[i] * pa_z[i] - ta2_xy_x_y_1[i] * pc_z[i];

        ta2_xy_xz_z_0[i] = ta1_y_z_z_1[i] + ta2_xy_z_z_0[i] * pa_x[i] - ta2_xy_z_z_1[i] * pc_x[i];
    }

    // Set up 27-30 components of targeted buffer : DP

    auto ta2_xy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 27);

    auto ta2_xy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 28);

    auto ta2_xy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 29);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_y_x_1, ta1_x_y_y_1, ta1_x_y_z_1, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_y_0_0, ta2_xy_y_0_1, ta2_xy_y_x_0, ta2_xy_y_x_1, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_y_z_0, ta2_xy_y_z_1, ta2_xy_yy_x_0, ta2_xy_yy_y_0, ta2_xy_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yy_x_0[i] = ta2_xy_0_x_0[i] * fe_0 - ta2_xy_0_x_1[i] * fe_0 + ta1_x_y_x_1[i] + ta2_xy_y_x_0[i] * pa_y[i] - ta2_xy_y_x_1[i] * pc_y[i];

        ta2_xy_yy_y_0[i] = ta2_xy_0_y_0[i] * fe_0 - ta2_xy_0_y_1[i] * fe_0 + ta2_xy_y_0_0[i] * fe_0 - ta2_xy_y_0_1[i] * fe_0 + ta1_x_y_y_1[i] + ta2_xy_y_y_0[i] * pa_y[i] - ta2_xy_y_y_1[i] * pc_y[i];

        ta2_xy_yy_z_0[i] = ta2_xy_0_z_0[i] * fe_0 - ta2_xy_0_z_1[i] * fe_0 + ta1_x_y_z_1[i] + ta2_xy_y_z_0[i] * pa_y[i] - ta2_xy_y_z_1[i] * pc_y[i];
    }

    // Set up 30-33 components of targeted buffer : DP

    auto ta2_xy_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 30);

    auto ta2_xy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 31);

    auto ta2_xy_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 32);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_z_z_1, ta2_xy_y_x_0, ta2_xy_y_x_1, ta2_xy_y_y_0, ta2_xy_y_y_1, ta2_xy_yz_x_0, ta2_xy_yz_y_0, ta2_xy_yz_z_0, ta2_xy_z_z_0, ta2_xy_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_yz_x_0[i] = ta2_xy_y_x_0[i] * pa_z[i] - ta2_xy_y_x_1[i] * pc_z[i];

        ta2_xy_yz_y_0[i] = ta2_xy_y_y_0[i] * pa_z[i] - ta2_xy_y_y_1[i] * pc_z[i];

        ta2_xy_yz_z_0[i] = ta1_x_z_z_1[i] + ta2_xy_z_z_0[i] * pa_y[i] - ta2_xy_z_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : DP

    auto ta2_xy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 33);

    auto ta2_xy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 34);

    auto ta2_xy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 35);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_z_0_0, ta2_xy_z_0_1, ta2_xy_z_x_0, ta2_xy_z_x_1, ta2_xy_z_y_0, ta2_xy_z_y_1, ta2_xy_z_z_0, ta2_xy_z_z_1, ta2_xy_zz_x_0, ta2_xy_zz_y_0, ta2_xy_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zz_x_0[i] = ta2_xy_0_x_0[i] * fe_0 - ta2_xy_0_x_1[i] * fe_0 + ta2_xy_z_x_0[i] * pa_z[i] - ta2_xy_z_x_1[i] * pc_z[i];

        ta2_xy_zz_y_0[i] = ta2_xy_0_y_0[i] * fe_0 - ta2_xy_0_y_1[i] * fe_0 + ta2_xy_z_y_0[i] * pa_z[i] - ta2_xy_z_y_1[i] * pc_z[i];

        ta2_xy_zz_z_0[i] = ta2_xy_0_z_0[i] * fe_0 - ta2_xy_0_z_1[i] * fe_0 + ta2_xy_z_0_0[i] * fe_0 - ta2_xy_z_0_1[i] * fe_0 + ta2_xy_z_z_0[i] * pa_z[i] - ta2_xy_z_z_1[i] * pc_z[i];
    }

    // Set up 36-39 components of targeted buffer : DP

    auto ta2_xz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 36);

    auto ta2_xz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 37);

    auto ta2_xz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 38);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_x_x_1, ta1_z_x_y_1, ta1_z_x_z_1, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_x_0_0, ta2_xz_x_0_1, ta2_xz_x_x_0, ta2_xz_x_x_1, ta2_xz_x_y_0, ta2_xz_x_y_1, ta2_xz_x_z_0, ta2_xz_x_z_1, ta2_xz_xx_x_0, ta2_xz_xx_y_0, ta2_xz_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xx_x_0[i] = ta2_xz_0_x_0[i] * fe_0 - ta2_xz_0_x_1[i] * fe_0 + ta2_xz_x_0_0[i] * fe_0 - ta2_xz_x_0_1[i] * fe_0 + ta1_z_x_x_1[i] + ta2_xz_x_x_0[i] * pa_x[i] - ta2_xz_x_x_1[i] * pc_x[i];

        ta2_xz_xx_y_0[i] = ta2_xz_0_y_0[i] * fe_0 - ta2_xz_0_y_1[i] * fe_0 + ta1_z_x_y_1[i] + ta2_xz_x_y_0[i] * pa_x[i] - ta2_xz_x_y_1[i] * pc_x[i];

        ta2_xz_xx_z_0[i] = ta2_xz_0_z_0[i] * fe_0 - ta2_xz_0_z_1[i] * fe_0 + ta1_z_x_z_1[i] + ta2_xz_x_z_0[i] * pa_x[i] - ta2_xz_x_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : DP

    auto ta2_xz_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 39);

    auto ta2_xz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 40);

    auto ta2_xz_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 41);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_y_y_1, ta2_xz_x_x_0, ta2_xz_x_x_1, ta2_xz_x_z_0, ta2_xz_x_z_1, ta2_xz_xy_x_0, ta2_xz_xy_y_0, ta2_xz_xy_z_0, ta2_xz_y_y_0, ta2_xz_y_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xy_x_0[i] = ta2_xz_x_x_0[i] * pa_y[i] - ta2_xz_x_x_1[i] * pc_y[i];

        ta2_xz_xy_y_0[i] = ta1_z_y_y_1[i] + ta2_xz_y_y_0[i] * pa_x[i] - ta2_xz_y_y_1[i] * pc_x[i];

        ta2_xz_xy_z_0[i] = ta2_xz_x_z_0[i] * pa_y[i] - ta2_xz_x_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : DP

    auto ta2_xz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 42);

    auto ta2_xz_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 43);

    auto ta2_xz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 44);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_x_x_1, ta1_z_z_y_1, ta1_z_z_z_1, ta2_xz_x_x_0, ta2_xz_x_x_1, ta2_xz_xz_x_0, ta2_xz_xz_y_0, ta2_xz_xz_z_0, ta2_xz_z_y_0, ta2_xz_z_y_1, ta2_xz_z_z_0, ta2_xz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xz_x_0[i] = ta1_x_x_x_1[i] + ta2_xz_x_x_0[i] * pa_z[i] - ta2_xz_x_x_1[i] * pc_z[i];

        ta2_xz_xz_y_0[i] = ta1_z_z_y_1[i] + ta2_xz_z_y_0[i] * pa_x[i] - ta2_xz_z_y_1[i] * pc_x[i];

        ta2_xz_xz_z_0[i] = ta1_z_z_z_1[i] + ta2_xz_z_z_0[i] * pa_x[i] - ta2_xz_z_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : DP

    auto ta2_xz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 45);

    auto ta2_xz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 46);

    auto ta2_xz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 47);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_y_0_0, ta2_xz_y_0_1, ta2_xz_y_x_0, ta2_xz_y_x_1, ta2_xz_y_y_0, ta2_xz_y_y_1, ta2_xz_y_z_0, ta2_xz_y_z_1, ta2_xz_yy_x_0, ta2_xz_yy_y_0, ta2_xz_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yy_x_0[i] = ta2_xz_0_x_0[i] * fe_0 - ta2_xz_0_x_1[i] * fe_0 + ta2_xz_y_x_0[i] * pa_y[i] - ta2_xz_y_x_1[i] * pc_y[i];

        ta2_xz_yy_y_0[i] = ta2_xz_0_y_0[i] * fe_0 - ta2_xz_0_y_1[i] * fe_0 + ta2_xz_y_0_0[i] * fe_0 - ta2_xz_y_0_1[i] * fe_0 + ta2_xz_y_y_0[i] * pa_y[i] - ta2_xz_y_y_1[i] * pc_y[i];

        ta2_xz_yy_z_0[i] = ta2_xz_0_z_0[i] * fe_0 - ta2_xz_0_z_1[i] * fe_0 + ta2_xz_y_z_0[i] * pa_y[i] - ta2_xz_y_z_1[i] * pc_y[i];
    }

    // Set up 48-51 components of targeted buffer : DP

    auto ta2_xz_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 48);

    auto ta2_xz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 49);

    auto ta2_xz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 50);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_y_y_1, ta2_xz_y_y_0, ta2_xz_y_y_1, ta2_xz_yz_x_0, ta2_xz_yz_y_0, ta2_xz_yz_z_0, ta2_xz_z_x_0, ta2_xz_z_x_1, ta2_xz_z_z_0, ta2_xz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_yz_x_0[i] = ta2_xz_z_x_0[i] * pa_y[i] - ta2_xz_z_x_1[i] * pc_y[i];

        ta2_xz_yz_y_0[i] = ta1_x_y_y_1[i] + ta2_xz_y_y_0[i] * pa_z[i] - ta2_xz_y_y_1[i] * pc_z[i];

        ta2_xz_yz_z_0[i] = ta2_xz_z_z_0[i] * pa_y[i] - ta2_xz_z_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : DP

    auto ta2_xz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 51);

    auto ta2_xz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 52);

    auto ta2_xz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 53);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_z_x_1, ta1_x_z_y_1, ta1_x_z_z_1, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_z_0_0, ta2_xz_z_0_1, ta2_xz_z_x_0, ta2_xz_z_x_1, ta2_xz_z_y_0, ta2_xz_z_y_1, ta2_xz_z_z_0, ta2_xz_z_z_1, ta2_xz_zz_x_0, ta2_xz_zz_y_0, ta2_xz_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zz_x_0[i] = ta2_xz_0_x_0[i] * fe_0 - ta2_xz_0_x_1[i] * fe_0 + ta1_x_z_x_1[i] + ta2_xz_z_x_0[i] * pa_z[i] - ta2_xz_z_x_1[i] * pc_z[i];

        ta2_xz_zz_y_0[i] = ta2_xz_0_y_0[i] * fe_0 - ta2_xz_0_y_1[i] * fe_0 + ta1_x_z_y_1[i] + ta2_xz_z_y_0[i] * pa_z[i] - ta2_xz_z_y_1[i] * pc_z[i];

        ta2_xz_zz_z_0[i] = ta2_xz_0_z_0[i] * fe_0 - ta2_xz_0_z_1[i] * fe_0 + ta2_xz_z_0_0[i] * fe_0 - ta2_xz_z_0_1[i] * fe_0 + ta1_x_z_z_1[i] + ta2_xz_z_z_0[i] * pa_z[i] - ta2_xz_z_z_1[i] * pc_z[i];
    }

    // Set up 54-57 components of targeted buffer : DP

    auto ta2_yy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 54);

    auto ta2_yy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 55);

    auto ta2_yy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 56);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_x_0_0, ta2_yy_x_0_1, ta2_yy_x_x_0, ta2_yy_x_x_1, ta2_yy_x_y_0, ta2_yy_x_y_1, ta2_yy_x_z_0, ta2_yy_x_z_1, ta2_yy_xx_x_0, ta2_yy_xx_y_0, ta2_yy_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xx_x_0[i] = ta2_yy_0_x_0[i] * fe_0 - ta2_yy_0_x_1[i] * fe_0 + ta2_yy_x_0_0[i] * fe_0 - ta2_yy_x_0_1[i] * fe_0 + ta2_yy_x_x_0[i] * pa_x[i] - ta2_yy_x_x_1[i] * pc_x[i];

        ta2_yy_xx_y_0[i] = ta2_yy_0_y_0[i] * fe_0 - ta2_yy_0_y_1[i] * fe_0 + ta2_yy_x_y_0[i] * pa_x[i] - ta2_yy_x_y_1[i] * pc_x[i];

        ta2_yy_xx_z_0[i] = ta2_yy_0_z_0[i] * fe_0 - ta2_yy_0_z_1[i] * fe_0 + ta2_yy_x_z_0[i] * pa_x[i] - ta2_yy_x_z_1[i] * pc_x[i];
    }

    // Set up 57-60 components of targeted buffer : DP

    auto ta2_yy_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 57);

    auto ta2_yy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 58);

    auto ta2_yy_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 59);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_x_x_1, ta2_yy_x_x_0, ta2_yy_x_x_1, ta2_yy_xy_x_0, ta2_yy_xy_y_0, ta2_yy_xy_z_0, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_y_z_0, ta2_yy_y_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xy_x_0[i] = 2.0 * ta1_y_x_x_1[i] + ta2_yy_x_x_0[i] * pa_y[i] - ta2_yy_x_x_1[i] * pc_y[i];

        ta2_yy_xy_y_0[i] = ta2_yy_y_y_0[i] * pa_x[i] - ta2_yy_y_y_1[i] * pc_x[i];

        ta2_yy_xy_z_0[i] = ta2_yy_y_z_0[i] * pa_x[i] - ta2_yy_y_z_1[i] * pc_x[i];
    }

    // Set up 60-63 components of targeted buffer : DP

    auto ta2_yy_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 60);

    auto ta2_yy_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 61);

    auto ta2_yy_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 62);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_x_x_0, ta2_yy_x_x_1, ta2_yy_xz_x_0, ta2_yy_xz_y_0, ta2_yy_xz_z_0, ta2_yy_z_y_0, ta2_yy_z_y_1, ta2_yy_z_z_0, ta2_yy_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xz_x_0[i] = ta2_yy_x_x_0[i] * pa_z[i] - ta2_yy_x_x_1[i] * pc_z[i];

        ta2_yy_xz_y_0[i] = ta2_yy_z_y_0[i] * pa_x[i] - ta2_yy_z_y_1[i] * pc_x[i];

        ta2_yy_xz_z_0[i] = ta2_yy_z_z_0[i] * pa_x[i] - ta2_yy_z_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : DP

    auto ta2_yy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 63);

    auto ta2_yy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 64);

    auto ta2_yy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 65);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_y_x_1, ta1_y_y_y_1, ta1_y_y_z_1, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_y_0_0, ta2_yy_y_0_1, ta2_yy_y_x_0, ta2_yy_y_x_1, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_y_z_0, ta2_yy_y_z_1, ta2_yy_yy_x_0, ta2_yy_yy_y_0, ta2_yy_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yy_x_0[i] = ta2_yy_0_x_0[i] * fe_0 - ta2_yy_0_x_1[i] * fe_0 + 2.0 * ta1_y_y_x_1[i] + ta2_yy_y_x_0[i] * pa_y[i] - ta2_yy_y_x_1[i] * pc_y[i];

        ta2_yy_yy_y_0[i] = ta2_yy_0_y_0[i] * fe_0 - ta2_yy_0_y_1[i] * fe_0 + ta2_yy_y_0_0[i] * fe_0 - ta2_yy_y_0_1[i] * fe_0 + 2.0 * ta1_y_y_y_1[i] + ta2_yy_y_y_0[i] * pa_y[i] - ta2_yy_y_y_1[i] * pc_y[i];

        ta2_yy_yy_z_0[i] = ta2_yy_0_z_0[i] * fe_0 - ta2_yy_0_z_1[i] * fe_0 + 2.0 * ta1_y_y_z_1[i] + ta2_yy_y_z_0[i] * pa_y[i] - ta2_yy_y_z_1[i] * pc_y[i];
    }

    // Set up 66-69 components of targeted buffer : DP

    auto ta2_yy_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 66);

    auto ta2_yy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 67);

    auto ta2_yy_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 68);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_z_z_1, ta2_yy_y_x_0, ta2_yy_y_x_1, ta2_yy_y_y_0, ta2_yy_y_y_1, ta2_yy_yz_x_0, ta2_yy_yz_y_0, ta2_yy_yz_z_0, ta2_yy_z_z_0, ta2_yy_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_yz_x_0[i] = ta2_yy_y_x_0[i] * pa_z[i] - ta2_yy_y_x_1[i] * pc_z[i];

        ta2_yy_yz_y_0[i] = ta2_yy_y_y_0[i] * pa_z[i] - ta2_yy_y_y_1[i] * pc_z[i];

        ta2_yy_yz_z_0[i] = 2.0 * ta1_y_z_z_1[i] + ta2_yy_z_z_0[i] * pa_y[i] - ta2_yy_z_z_1[i] * pc_y[i];
    }

    // Set up 69-72 components of targeted buffer : DP

    auto ta2_yy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 69);

    auto ta2_yy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 70);

    auto ta2_yy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 71);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_z_0_0, ta2_yy_z_0_1, ta2_yy_z_x_0, ta2_yy_z_x_1, ta2_yy_z_y_0, ta2_yy_z_y_1, ta2_yy_z_z_0, ta2_yy_z_z_1, ta2_yy_zz_x_0, ta2_yy_zz_y_0, ta2_yy_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zz_x_0[i] = ta2_yy_0_x_0[i] * fe_0 - ta2_yy_0_x_1[i] * fe_0 + ta2_yy_z_x_0[i] * pa_z[i] - ta2_yy_z_x_1[i] * pc_z[i];

        ta2_yy_zz_y_0[i] = ta2_yy_0_y_0[i] * fe_0 - ta2_yy_0_y_1[i] * fe_0 + ta2_yy_z_y_0[i] * pa_z[i] - ta2_yy_z_y_1[i] * pc_z[i];

        ta2_yy_zz_z_0[i] = ta2_yy_0_z_0[i] * fe_0 - ta2_yy_0_z_1[i] * fe_0 + ta2_yy_z_0_0[i] * fe_0 - ta2_yy_z_0_1[i] * fe_0 + ta2_yy_z_z_0[i] * pa_z[i] - ta2_yy_z_z_1[i] * pc_z[i];
    }

    // Set up 72-75 components of targeted buffer : DP

    auto ta2_yz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 72);

    auto ta2_yz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 73);

    auto ta2_yz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 74);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_x_0_0, ta2_yz_x_0_1, ta2_yz_x_x_0, ta2_yz_x_x_1, ta2_yz_x_y_0, ta2_yz_x_y_1, ta2_yz_x_z_0, ta2_yz_x_z_1, ta2_yz_xx_x_0, ta2_yz_xx_y_0, ta2_yz_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xx_x_0[i] = ta2_yz_0_x_0[i] * fe_0 - ta2_yz_0_x_1[i] * fe_0 + ta2_yz_x_0_0[i] * fe_0 - ta2_yz_x_0_1[i] * fe_0 + ta2_yz_x_x_0[i] * pa_x[i] - ta2_yz_x_x_1[i] * pc_x[i];

        ta2_yz_xx_y_0[i] = ta2_yz_0_y_0[i] * fe_0 - ta2_yz_0_y_1[i] * fe_0 + ta2_yz_x_y_0[i] * pa_x[i] - ta2_yz_x_y_1[i] * pc_x[i];

        ta2_yz_xx_z_0[i] = ta2_yz_0_z_0[i] * fe_0 - ta2_yz_0_z_1[i] * fe_0 + ta2_yz_x_z_0[i] * pa_x[i] - ta2_yz_x_z_1[i] * pc_x[i];
    }

    // Set up 75-78 components of targeted buffer : DP

    auto ta2_yz_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 75);

    auto ta2_yz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 76);

    auto ta2_yz_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 77);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_x_x_1, ta2_yz_x_x_0, ta2_yz_x_x_1, ta2_yz_xy_x_0, ta2_yz_xy_y_0, ta2_yz_xy_z_0, ta2_yz_y_y_0, ta2_yz_y_y_1, ta2_yz_y_z_0, ta2_yz_y_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xy_x_0[i] = ta1_z_x_x_1[i] + ta2_yz_x_x_0[i] * pa_y[i] - ta2_yz_x_x_1[i] * pc_y[i];

        ta2_yz_xy_y_0[i] = ta2_yz_y_y_0[i] * pa_x[i] - ta2_yz_y_y_1[i] * pc_x[i];

        ta2_yz_xy_z_0[i] = ta2_yz_y_z_0[i] * pa_x[i] - ta2_yz_y_z_1[i] * pc_x[i];
    }

    // Set up 78-81 components of targeted buffer : DP

    auto ta2_yz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 78);

    auto ta2_yz_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 79);

    auto ta2_yz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 80);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_x_x_1, ta2_yz_x_x_0, ta2_yz_x_x_1, ta2_yz_xz_x_0, ta2_yz_xz_y_0, ta2_yz_xz_z_0, ta2_yz_z_y_0, ta2_yz_z_y_1, ta2_yz_z_z_0, ta2_yz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xz_x_0[i] = ta1_y_x_x_1[i] + ta2_yz_x_x_0[i] * pa_z[i] - ta2_yz_x_x_1[i] * pc_z[i];

        ta2_yz_xz_y_0[i] = ta2_yz_z_y_0[i] * pa_x[i] - ta2_yz_z_y_1[i] * pc_x[i];

        ta2_yz_xz_z_0[i] = ta2_yz_z_z_0[i] * pa_x[i] - ta2_yz_z_z_1[i] * pc_x[i];
    }

    // Set up 81-84 components of targeted buffer : DP

    auto ta2_yz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 81);

    auto ta2_yz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 82);

    auto ta2_yz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 83);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_y_x_1, ta1_z_y_y_1, ta1_z_y_z_1, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_y_0_0, ta2_yz_y_0_1, ta2_yz_y_x_0, ta2_yz_y_x_1, ta2_yz_y_y_0, ta2_yz_y_y_1, ta2_yz_y_z_0, ta2_yz_y_z_1, ta2_yz_yy_x_0, ta2_yz_yy_y_0, ta2_yz_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yy_x_0[i] = ta2_yz_0_x_0[i] * fe_0 - ta2_yz_0_x_1[i] * fe_0 + ta1_z_y_x_1[i] + ta2_yz_y_x_0[i] * pa_y[i] - ta2_yz_y_x_1[i] * pc_y[i];

        ta2_yz_yy_y_0[i] = ta2_yz_0_y_0[i] * fe_0 - ta2_yz_0_y_1[i] * fe_0 + ta2_yz_y_0_0[i] * fe_0 - ta2_yz_y_0_1[i] * fe_0 + ta1_z_y_y_1[i] + ta2_yz_y_y_0[i] * pa_y[i] - ta2_yz_y_y_1[i] * pc_y[i];

        ta2_yz_yy_z_0[i] = ta2_yz_0_z_0[i] * fe_0 - ta2_yz_0_z_1[i] * fe_0 + ta1_z_y_z_1[i] + ta2_yz_y_z_0[i] * pa_y[i] - ta2_yz_y_z_1[i] * pc_y[i];
    }

    // Set up 84-87 components of targeted buffer : DP

    auto ta2_yz_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 84);

    auto ta2_yz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 85);

    auto ta2_yz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 86);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_y_y_1, ta1_z_z_x_1, ta1_z_z_z_1, ta2_yz_y_y_0, ta2_yz_y_y_1, ta2_yz_yz_x_0, ta2_yz_yz_y_0, ta2_yz_yz_z_0, ta2_yz_z_x_0, ta2_yz_z_x_1, ta2_yz_z_z_0, ta2_yz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_yz_x_0[i] = ta1_z_z_x_1[i] + ta2_yz_z_x_0[i] * pa_y[i] - ta2_yz_z_x_1[i] * pc_y[i];

        ta2_yz_yz_y_0[i] = ta1_y_y_y_1[i] + ta2_yz_y_y_0[i] * pa_z[i] - ta2_yz_y_y_1[i] * pc_z[i];

        ta2_yz_yz_z_0[i] = ta1_z_z_z_1[i] + ta2_yz_z_z_0[i] * pa_y[i] - ta2_yz_z_z_1[i] * pc_y[i];
    }

    // Set up 87-90 components of targeted buffer : DP

    auto ta2_yz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 87);

    auto ta2_yz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 88);

    auto ta2_yz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 89);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_z_x_1, ta1_y_z_y_1, ta1_y_z_z_1, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_z_0_0, ta2_yz_z_0_1, ta2_yz_z_x_0, ta2_yz_z_x_1, ta2_yz_z_y_0, ta2_yz_z_y_1, ta2_yz_z_z_0, ta2_yz_z_z_1, ta2_yz_zz_x_0, ta2_yz_zz_y_0, ta2_yz_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zz_x_0[i] = ta2_yz_0_x_0[i] * fe_0 - ta2_yz_0_x_1[i] * fe_0 + ta1_y_z_x_1[i] + ta2_yz_z_x_0[i] * pa_z[i] - ta2_yz_z_x_1[i] * pc_z[i];

        ta2_yz_zz_y_0[i] = ta2_yz_0_y_0[i] * fe_0 - ta2_yz_0_y_1[i] * fe_0 + ta1_y_z_y_1[i] + ta2_yz_z_y_0[i] * pa_z[i] - ta2_yz_z_y_1[i] * pc_z[i];

        ta2_yz_zz_z_0[i] = ta2_yz_0_z_0[i] * fe_0 - ta2_yz_0_z_1[i] * fe_0 + ta2_yz_z_0_0[i] * fe_0 - ta2_yz_z_0_1[i] * fe_0 + ta1_y_z_z_1[i] + ta2_yz_z_z_0[i] * pa_z[i] - ta2_yz_z_z_1[i] * pc_z[i];
    }

    // Set up 90-93 components of targeted buffer : DP

    auto ta2_zz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 90);

    auto ta2_zz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 91);

    auto ta2_zz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 92);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_x_0_0, ta2_zz_x_0_1, ta2_zz_x_x_0, ta2_zz_x_x_1, ta2_zz_x_y_0, ta2_zz_x_y_1, ta2_zz_x_z_0, ta2_zz_x_z_1, ta2_zz_xx_x_0, ta2_zz_xx_y_0, ta2_zz_xx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xx_x_0[i] = ta2_zz_0_x_0[i] * fe_0 - ta2_zz_0_x_1[i] * fe_0 + ta2_zz_x_0_0[i] * fe_0 - ta2_zz_x_0_1[i] * fe_0 + ta2_zz_x_x_0[i] * pa_x[i] - ta2_zz_x_x_1[i] * pc_x[i];

        ta2_zz_xx_y_0[i] = ta2_zz_0_y_0[i] * fe_0 - ta2_zz_0_y_1[i] * fe_0 + ta2_zz_x_y_0[i] * pa_x[i] - ta2_zz_x_y_1[i] * pc_x[i];

        ta2_zz_xx_z_0[i] = ta2_zz_0_z_0[i] * fe_0 - ta2_zz_0_z_1[i] * fe_0 + ta2_zz_x_z_0[i] * pa_x[i] - ta2_zz_x_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : DP

    auto ta2_zz_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 93);

    auto ta2_zz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 94);

    auto ta2_zz_xy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 95);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_x_x_0, ta2_zz_x_x_1, ta2_zz_xy_x_0, ta2_zz_xy_y_0, ta2_zz_xy_z_0, ta2_zz_y_y_0, ta2_zz_y_y_1, ta2_zz_y_z_0, ta2_zz_y_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xy_x_0[i] = ta2_zz_x_x_0[i] * pa_y[i] - ta2_zz_x_x_1[i] * pc_y[i];

        ta2_zz_xy_y_0[i] = ta2_zz_y_y_0[i] * pa_x[i] - ta2_zz_y_y_1[i] * pc_x[i];

        ta2_zz_xy_z_0[i] = ta2_zz_y_z_0[i] * pa_x[i] - ta2_zz_y_z_1[i] * pc_x[i];
    }

    // Set up 96-99 components of targeted buffer : DP

    auto ta2_zz_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 96);

    auto ta2_zz_xz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 97);

    auto ta2_zz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 98);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_x_x_1, ta2_zz_x_x_0, ta2_zz_x_x_1, ta2_zz_xz_x_0, ta2_zz_xz_y_0, ta2_zz_xz_z_0, ta2_zz_z_y_0, ta2_zz_z_y_1, ta2_zz_z_z_0, ta2_zz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xz_x_0[i] = 2.0 * ta1_z_x_x_1[i] + ta2_zz_x_x_0[i] * pa_z[i] - ta2_zz_x_x_1[i] * pc_z[i];

        ta2_zz_xz_y_0[i] = ta2_zz_z_y_0[i] * pa_x[i] - ta2_zz_z_y_1[i] * pc_x[i];

        ta2_zz_xz_z_0[i] = ta2_zz_z_z_0[i] * pa_x[i] - ta2_zz_z_z_1[i] * pc_x[i];
    }

    // Set up 99-102 components of targeted buffer : DP

    auto ta2_zz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 99);

    auto ta2_zz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 100);

    auto ta2_zz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 101);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_y_0_0, ta2_zz_y_0_1, ta2_zz_y_x_0, ta2_zz_y_x_1, ta2_zz_y_y_0, ta2_zz_y_y_1, ta2_zz_y_z_0, ta2_zz_y_z_1, ta2_zz_yy_x_0, ta2_zz_yy_y_0, ta2_zz_yy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yy_x_0[i] = ta2_zz_0_x_0[i] * fe_0 - ta2_zz_0_x_1[i] * fe_0 + ta2_zz_y_x_0[i] * pa_y[i] - ta2_zz_y_x_1[i] * pc_y[i];

        ta2_zz_yy_y_0[i] = ta2_zz_0_y_0[i] * fe_0 - ta2_zz_0_y_1[i] * fe_0 + ta2_zz_y_0_0[i] * fe_0 - ta2_zz_y_0_1[i] * fe_0 + ta2_zz_y_y_0[i] * pa_y[i] - ta2_zz_y_y_1[i] * pc_y[i];

        ta2_zz_yy_z_0[i] = ta2_zz_0_z_0[i] * fe_0 - ta2_zz_0_z_1[i] * fe_0 + ta2_zz_y_z_0[i] * pa_y[i] - ta2_zz_y_z_1[i] * pc_y[i];
    }

    // Set up 102-105 components of targeted buffer : DP

    auto ta2_zz_yz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 102);

    auto ta2_zz_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 103);

    auto ta2_zz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 104);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_y_y_1, ta2_zz_y_y_0, ta2_zz_y_y_1, ta2_zz_yz_x_0, ta2_zz_yz_y_0, ta2_zz_yz_z_0, ta2_zz_z_x_0, ta2_zz_z_x_1, ta2_zz_z_z_0, ta2_zz_z_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_yz_x_0[i] = ta2_zz_z_x_0[i] * pa_y[i] - ta2_zz_z_x_1[i] * pc_y[i];

        ta2_zz_yz_y_0[i] = 2.0 * ta1_z_y_y_1[i] + ta2_zz_y_y_0[i] * pa_z[i] - ta2_zz_y_y_1[i] * pc_z[i];

        ta2_zz_yz_z_0[i] = ta2_zz_z_z_0[i] * pa_y[i] - ta2_zz_z_z_1[i] * pc_y[i];
    }

    // Set up 105-108 components of targeted buffer : DP

    auto ta2_zz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 105);

    auto ta2_zz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 106);

    auto ta2_zz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 107);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_z_x_1, ta1_z_z_y_1, ta1_z_z_z_1, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_z_0_0, ta2_zz_z_0_1, ta2_zz_z_x_0, ta2_zz_z_x_1, ta2_zz_z_y_0, ta2_zz_z_y_1, ta2_zz_z_z_0, ta2_zz_z_z_1, ta2_zz_zz_x_0, ta2_zz_zz_y_0, ta2_zz_zz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zz_x_0[i] = ta2_zz_0_x_0[i] * fe_0 - ta2_zz_0_x_1[i] * fe_0 + 2.0 * ta1_z_z_x_1[i] + ta2_zz_z_x_0[i] * pa_z[i] - ta2_zz_z_x_1[i] * pc_z[i];

        ta2_zz_zz_y_0[i] = ta2_zz_0_y_0[i] * fe_0 - ta2_zz_0_y_1[i] * fe_0 + 2.0 * ta1_z_z_y_1[i] + ta2_zz_z_y_0[i] * pa_z[i] - ta2_zz_z_y_1[i] * pc_z[i];

        ta2_zz_zz_z_0[i] = ta2_zz_0_z_0[i] * fe_0 - ta2_zz_0_z_1[i] * fe_0 + ta2_zz_z_0_0[i] * fe_0 - ta2_zz_z_0_1[i] * fe_0 + 2.0 * ta1_z_z_z_1[i] + ta2_zz_z_z_0[i] * pa_z[i] - ta2_zz_z_z_1[i] * pc_z[i];
    }

}

} // npotrec namespace

