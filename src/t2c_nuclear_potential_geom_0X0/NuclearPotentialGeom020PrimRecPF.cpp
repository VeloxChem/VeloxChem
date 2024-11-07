#include "NuclearPotentialGeom020PrimRecPF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_pf(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_pf,
                                        const size_t              idx_npot_geom_020_0_sd,
                                        const size_t              idx_npot_geom_020_1_sd,
                                        const size_t              idx_npot_geom_010_1_sf,
                                        const size_t              idx_npot_geom_020_0_sf,
                                        const size_t              idx_npot_geom_020_1_sf,
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

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 1);

    auto ta1_x_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 2);

    auto ta1_x_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 3);

    auto ta1_x_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 4);

    auto ta1_x_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 5);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 7);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 11);

    auto ta1_y_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 12);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 14);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 21);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 24);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 27);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up components of auxiliary buffer : SF

    auto ta2_xx_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf);

    auto ta2_xx_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 1);

    auto ta2_xx_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 2);

    auto ta2_xx_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 3);

    auto ta2_xx_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 4);

    auto ta2_xx_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 5);

    auto ta2_xx_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 6);

    auto ta2_xx_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 7);

    auto ta2_xx_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 8);

    auto ta2_xx_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 9);

    auto ta2_xy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 10);

    auto ta2_xy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 11);

    auto ta2_xy_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 12);

    auto ta2_xy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 13);

    auto ta2_xy_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 14);

    auto ta2_xy_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 15);

    auto ta2_xy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 16);

    auto ta2_xy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 17);

    auto ta2_xy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 18);

    auto ta2_xy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 19);

    auto ta2_xz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 20);

    auto ta2_xz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 21);

    auto ta2_xz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 22);

    auto ta2_xz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 23);

    auto ta2_xz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 24);

    auto ta2_xz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 25);

    auto ta2_xz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 26);

    auto ta2_xz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 27);

    auto ta2_xz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 28);

    auto ta2_xz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 29);

    auto ta2_yy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 30);

    auto ta2_yy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 31);

    auto ta2_yy_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 32);

    auto ta2_yy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 33);

    auto ta2_yy_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 34);

    auto ta2_yy_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 35);

    auto ta2_yy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 36);

    auto ta2_yy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 37);

    auto ta2_yy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 38);

    auto ta2_yy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 39);

    auto ta2_yz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 40);

    auto ta2_yz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 41);

    auto ta2_yz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 42);

    auto ta2_yz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 43);

    auto ta2_yz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 44);

    auto ta2_yz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 45);

    auto ta2_yz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 46);

    auto ta2_yz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 47);

    auto ta2_yz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 48);

    auto ta2_yz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 49);

    auto ta2_zz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 50);

    auto ta2_zz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 51);

    auto ta2_zz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 52);

    auto ta2_zz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 53);

    auto ta2_zz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 54);

    auto ta2_zz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 55);

    auto ta2_zz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 56);

    auto ta2_zz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 57);

    auto ta2_zz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 58);

    auto ta2_zz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 59);

    // Set up components of auxiliary buffer : SF

    auto ta2_xx_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf);

    auto ta2_xx_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 1);

    auto ta2_xx_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 2);

    auto ta2_xx_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 3);

    auto ta2_xx_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 4);

    auto ta2_xx_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 5);

    auto ta2_xx_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 6);

    auto ta2_xx_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 7);

    auto ta2_xx_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 8);

    auto ta2_xx_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 9);

    auto ta2_xy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 10);

    auto ta2_xy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 11);

    auto ta2_xy_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 12);

    auto ta2_xy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 13);

    auto ta2_xy_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 14);

    auto ta2_xy_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 15);

    auto ta2_xy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 16);

    auto ta2_xy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 17);

    auto ta2_xy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 18);

    auto ta2_xy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 19);

    auto ta2_xz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 20);

    auto ta2_xz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 21);

    auto ta2_xz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 22);

    auto ta2_xz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 23);

    auto ta2_xz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 24);

    auto ta2_xz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 25);

    auto ta2_xz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 26);

    auto ta2_xz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 27);

    auto ta2_xz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 28);

    auto ta2_xz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 29);

    auto ta2_yy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 30);

    auto ta2_yy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 31);

    auto ta2_yy_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 32);

    auto ta2_yy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 33);

    auto ta2_yy_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 34);

    auto ta2_yy_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 35);

    auto ta2_yy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 36);

    auto ta2_yy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 37);

    auto ta2_yy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 38);

    auto ta2_yy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 39);

    auto ta2_yz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 40);

    auto ta2_yz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 41);

    auto ta2_yz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 42);

    auto ta2_yz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 43);

    auto ta2_yz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 44);

    auto ta2_yz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 45);

    auto ta2_yz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 46);

    auto ta2_yz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 47);

    auto ta2_yz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 48);

    auto ta2_yz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 49);

    auto ta2_zz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 50);

    auto ta2_zz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 51);

    auto ta2_zz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 52);

    auto ta2_zz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 53);

    auto ta2_zz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 54);

    auto ta2_zz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 55);

    auto ta2_zz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 56);

    auto ta2_zz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 57);

    auto ta2_zz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 58);

    auto ta2_zz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 59);

    // Set up 0-10 components of targeted buffer : PF

    auto ta2_xx_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf);

    auto ta2_xx_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 1);

    auto ta2_xx_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 2);

    auto ta2_xx_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 3);

    auto ta2_xx_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 4);

    auto ta2_xx_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 5);

    auto ta2_xx_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 6);

    auto ta2_xx_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 7);

    auto ta2_xx_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 8);

    auto ta2_xx_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 9);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_1,  \
                             ta2_xx_0_xx_0,  \
                             ta2_xx_0_xx_1,  \
                             ta2_xx_0_xxx_0, \
                             ta2_xx_0_xxx_1, \
                             ta2_xx_0_xxy_0, \
                             ta2_xx_0_xxy_1, \
                             ta2_xx_0_xxz_0, \
                             ta2_xx_0_xxz_1, \
                             ta2_xx_0_xy_0,  \
                             ta2_xx_0_xy_1,  \
                             ta2_xx_0_xyy_0, \
                             ta2_xx_0_xyy_1, \
                             ta2_xx_0_xyz_0, \
                             ta2_xx_0_xyz_1, \
                             ta2_xx_0_xz_0,  \
                             ta2_xx_0_xz_1,  \
                             ta2_xx_0_xzz_0, \
                             ta2_xx_0_xzz_1, \
                             ta2_xx_0_yy_0,  \
                             ta2_xx_0_yy_1,  \
                             ta2_xx_0_yyy_0, \
                             ta2_xx_0_yyy_1, \
                             ta2_xx_0_yyz_0, \
                             ta2_xx_0_yyz_1, \
                             ta2_xx_0_yz_0,  \
                             ta2_xx_0_yz_1,  \
                             ta2_xx_0_yzz_0, \
                             ta2_xx_0_yzz_1, \
                             ta2_xx_0_zz_0,  \
                             ta2_xx_0_zz_1,  \
                             ta2_xx_0_zzz_0, \
                             ta2_xx_0_zzz_1, \
                             ta2_xx_x_xxx_0, \
                             ta2_xx_x_xxy_0, \
                             ta2_xx_x_xxz_0, \
                             ta2_xx_x_xyy_0, \
                             ta2_xx_x_xyz_0, \
                             ta2_xx_x_xzz_0, \
                             ta2_xx_x_yyy_0, \
                             ta2_xx_x_yyz_0, \
                             ta2_xx_x_yzz_0, \
                             ta2_xx_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_x_xxx_0[i] = 3.0 * ta2_xx_0_xx_0[i] * fe_0 - 3.0 * ta2_xx_0_xx_1[i] * fe_0 + 2.0 * ta1_x_0_xxx_1[i] + ta2_xx_0_xxx_0[i] * pa_x[i] -
                            ta2_xx_0_xxx_1[i] * pc_x[i];

        ta2_xx_x_xxy_0[i] = 2.0 * ta2_xx_0_xy_0[i] * fe_0 - 2.0 * ta2_xx_0_xy_1[i] * fe_0 + 2.0 * ta1_x_0_xxy_1[i] + ta2_xx_0_xxy_0[i] * pa_x[i] -
                            ta2_xx_0_xxy_1[i] * pc_x[i];

        ta2_xx_x_xxz_0[i] = 2.0 * ta2_xx_0_xz_0[i] * fe_0 - 2.0 * ta2_xx_0_xz_1[i] * fe_0 + 2.0 * ta1_x_0_xxz_1[i] + ta2_xx_0_xxz_0[i] * pa_x[i] -
                            ta2_xx_0_xxz_1[i] * pc_x[i];

        ta2_xx_x_xyy_0[i] =
            ta2_xx_0_yy_0[i] * fe_0 - ta2_xx_0_yy_1[i] * fe_0 + 2.0 * ta1_x_0_xyy_1[i] + ta2_xx_0_xyy_0[i] * pa_x[i] - ta2_xx_0_xyy_1[i] * pc_x[i];

        ta2_xx_x_xyz_0[i] =
            ta2_xx_0_yz_0[i] * fe_0 - ta2_xx_0_yz_1[i] * fe_0 + 2.0 * ta1_x_0_xyz_1[i] + ta2_xx_0_xyz_0[i] * pa_x[i] - ta2_xx_0_xyz_1[i] * pc_x[i];

        ta2_xx_x_xzz_0[i] =
            ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + 2.0 * ta1_x_0_xzz_1[i] + ta2_xx_0_xzz_0[i] * pa_x[i] - ta2_xx_0_xzz_1[i] * pc_x[i];

        ta2_xx_x_yyy_0[i] = 2.0 * ta1_x_0_yyy_1[i] + ta2_xx_0_yyy_0[i] * pa_x[i] - ta2_xx_0_yyy_1[i] * pc_x[i];

        ta2_xx_x_yyz_0[i] = 2.0 * ta1_x_0_yyz_1[i] + ta2_xx_0_yyz_0[i] * pa_x[i] - ta2_xx_0_yyz_1[i] * pc_x[i];

        ta2_xx_x_yzz_0[i] = 2.0 * ta1_x_0_yzz_1[i] + ta2_xx_0_yzz_0[i] * pa_x[i] - ta2_xx_0_yzz_1[i] * pc_x[i];

        ta2_xx_x_zzz_0[i] = 2.0 * ta1_x_0_zzz_1[i] + ta2_xx_0_zzz_0[i] * pa_x[i] - ta2_xx_0_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto ta2_xx_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 10);

    auto ta2_xx_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 11);

    auto ta2_xx_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 12);

    auto ta2_xx_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 13);

    auto ta2_xx_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 14);

    auto ta2_xx_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 15);

    auto ta2_xx_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 16);

    auto ta2_xx_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 17);

    auto ta2_xx_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 18);

    auto ta2_xx_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 19);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xx_0_xx_0,  \
                             ta2_xx_0_xx_1,  \
                             ta2_xx_0_xxx_0, \
                             ta2_xx_0_xxx_1, \
                             ta2_xx_0_xxy_0, \
                             ta2_xx_0_xxy_1, \
                             ta2_xx_0_xxz_0, \
                             ta2_xx_0_xxz_1, \
                             ta2_xx_0_xy_0,  \
                             ta2_xx_0_xy_1,  \
                             ta2_xx_0_xyy_0, \
                             ta2_xx_0_xyy_1, \
                             ta2_xx_0_xyz_0, \
                             ta2_xx_0_xyz_1, \
                             ta2_xx_0_xz_0,  \
                             ta2_xx_0_xz_1,  \
                             ta2_xx_0_xzz_0, \
                             ta2_xx_0_xzz_1, \
                             ta2_xx_0_yy_0,  \
                             ta2_xx_0_yy_1,  \
                             ta2_xx_0_yyy_0, \
                             ta2_xx_0_yyy_1, \
                             ta2_xx_0_yyz_0, \
                             ta2_xx_0_yyz_1, \
                             ta2_xx_0_yz_0,  \
                             ta2_xx_0_yz_1,  \
                             ta2_xx_0_yzz_0, \
                             ta2_xx_0_yzz_1, \
                             ta2_xx_0_zz_0,  \
                             ta2_xx_0_zz_1,  \
                             ta2_xx_0_zzz_0, \
                             ta2_xx_0_zzz_1, \
                             ta2_xx_y_xxx_0, \
                             ta2_xx_y_xxy_0, \
                             ta2_xx_y_xxz_0, \
                             ta2_xx_y_xyy_0, \
                             ta2_xx_y_xyz_0, \
                             ta2_xx_y_xzz_0, \
                             ta2_xx_y_yyy_0, \
                             ta2_xx_y_yyz_0, \
                             ta2_xx_y_yzz_0, \
                             ta2_xx_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_y_xxx_0[i] = ta2_xx_0_xxx_0[i] * pa_y[i] - ta2_xx_0_xxx_1[i] * pc_y[i];

        ta2_xx_y_xxy_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_0_xxy_0[i] * pa_y[i] - ta2_xx_0_xxy_1[i] * pc_y[i];

        ta2_xx_y_xxz_0[i] = ta2_xx_0_xxz_0[i] * pa_y[i] - ta2_xx_0_xxz_1[i] * pc_y[i];

        ta2_xx_y_xyy_0[i] = 2.0 * ta2_xx_0_xy_0[i] * fe_0 - 2.0 * ta2_xx_0_xy_1[i] * fe_0 + ta2_xx_0_xyy_0[i] * pa_y[i] - ta2_xx_0_xyy_1[i] * pc_y[i];

        ta2_xx_y_xyz_0[i] = ta2_xx_0_xz_0[i] * fe_0 - ta2_xx_0_xz_1[i] * fe_0 + ta2_xx_0_xyz_0[i] * pa_y[i] - ta2_xx_0_xyz_1[i] * pc_y[i];

        ta2_xx_y_xzz_0[i] = ta2_xx_0_xzz_0[i] * pa_y[i] - ta2_xx_0_xzz_1[i] * pc_y[i];

        ta2_xx_y_yyy_0[i] = 3.0 * ta2_xx_0_yy_0[i] * fe_0 - 3.0 * ta2_xx_0_yy_1[i] * fe_0 + ta2_xx_0_yyy_0[i] * pa_y[i] - ta2_xx_0_yyy_1[i] * pc_y[i];

        ta2_xx_y_yyz_0[i] = 2.0 * ta2_xx_0_yz_0[i] * fe_0 - 2.0 * ta2_xx_0_yz_1[i] * fe_0 + ta2_xx_0_yyz_0[i] * pa_y[i] - ta2_xx_0_yyz_1[i] * pc_y[i];

        ta2_xx_y_yzz_0[i] = ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + ta2_xx_0_yzz_0[i] * pa_y[i] - ta2_xx_0_yzz_1[i] * pc_y[i];

        ta2_xx_y_zzz_0[i] = ta2_xx_0_zzz_0[i] * pa_y[i] - ta2_xx_0_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto ta2_xx_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 20);

    auto ta2_xx_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 21);

    auto ta2_xx_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 22);

    auto ta2_xx_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 23);

    auto ta2_xx_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 24);

    auto ta2_xx_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 25);

    auto ta2_xx_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 26);

    auto ta2_xx_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 27);

    auto ta2_xx_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 28);

    auto ta2_xx_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 29);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xx_0_xx_0,  \
                             ta2_xx_0_xx_1,  \
                             ta2_xx_0_xxx_0, \
                             ta2_xx_0_xxx_1, \
                             ta2_xx_0_xxy_0, \
                             ta2_xx_0_xxy_1, \
                             ta2_xx_0_xxz_0, \
                             ta2_xx_0_xxz_1, \
                             ta2_xx_0_xy_0,  \
                             ta2_xx_0_xy_1,  \
                             ta2_xx_0_xyy_0, \
                             ta2_xx_0_xyy_1, \
                             ta2_xx_0_xyz_0, \
                             ta2_xx_0_xyz_1, \
                             ta2_xx_0_xz_0,  \
                             ta2_xx_0_xz_1,  \
                             ta2_xx_0_xzz_0, \
                             ta2_xx_0_xzz_1, \
                             ta2_xx_0_yy_0,  \
                             ta2_xx_0_yy_1,  \
                             ta2_xx_0_yyy_0, \
                             ta2_xx_0_yyy_1, \
                             ta2_xx_0_yyz_0, \
                             ta2_xx_0_yyz_1, \
                             ta2_xx_0_yz_0,  \
                             ta2_xx_0_yz_1,  \
                             ta2_xx_0_yzz_0, \
                             ta2_xx_0_yzz_1, \
                             ta2_xx_0_zz_0,  \
                             ta2_xx_0_zz_1,  \
                             ta2_xx_0_zzz_0, \
                             ta2_xx_0_zzz_1, \
                             ta2_xx_z_xxx_0, \
                             ta2_xx_z_xxy_0, \
                             ta2_xx_z_xxz_0, \
                             ta2_xx_z_xyy_0, \
                             ta2_xx_z_xyz_0, \
                             ta2_xx_z_xzz_0, \
                             ta2_xx_z_yyy_0, \
                             ta2_xx_z_yyz_0, \
                             ta2_xx_z_yzz_0, \
                             ta2_xx_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_z_xxx_0[i] = ta2_xx_0_xxx_0[i] * pa_z[i] - ta2_xx_0_xxx_1[i] * pc_z[i];

        ta2_xx_z_xxy_0[i] = ta2_xx_0_xxy_0[i] * pa_z[i] - ta2_xx_0_xxy_1[i] * pc_z[i];

        ta2_xx_z_xxz_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_0_xxz_0[i] * pa_z[i] - ta2_xx_0_xxz_1[i] * pc_z[i];

        ta2_xx_z_xyy_0[i] = ta2_xx_0_xyy_0[i] * pa_z[i] - ta2_xx_0_xyy_1[i] * pc_z[i];

        ta2_xx_z_xyz_0[i] = ta2_xx_0_xy_0[i] * fe_0 - ta2_xx_0_xy_1[i] * fe_0 + ta2_xx_0_xyz_0[i] * pa_z[i] - ta2_xx_0_xyz_1[i] * pc_z[i];

        ta2_xx_z_xzz_0[i] = 2.0 * ta2_xx_0_xz_0[i] * fe_0 - 2.0 * ta2_xx_0_xz_1[i] * fe_0 + ta2_xx_0_xzz_0[i] * pa_z[i] - ta2_xx_0_xzz_1[i] * pc_z[i];

        ta2_xx_z_yyy_0[i] = ta2_xx_0_yyy_0[i] * pa_z[i] - ta2_xx_0_yyy_1[i] * pc_z[i];

        ta2_xx_z_yyz_0[i] = ta2_xx_0_yy_0[i] * fe_0 - ta2_xx_0_yy_1[i] * fe_0 + ta2_xx_0_yyz_0[i] * pa_z[i] - ta2_xx_0_yyz_1[i] * pc_z[i];

        ta2_xx_z_yzz_0[i] = 2.0 * ta2_xx_0_yz_0[i] * fe_0 - 2.0 * ta2_xx_0_yz_1[i] * fe_0 + ta2_xx_0_yzz_0[i] * pa_z[i] - ta2_xx_0_yzz_1[i] * pc_z[i];

        ta2_xx_z_zzz_0[i] = 3.0 * ta2_xx_0_zz_0[i] * fe_0 - 3.0 * ta2_xx_0_zz_1[i] * fe_0 + ta2_xx_0_zzz_0[i] * pa_z[i] - ta2_xx_0_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : PF

    auto ta2_xy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 30);

    auto ta2_xy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 31);

    auto ta2_xy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 32);

    auto ta2_xy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 33);

    auto ta2_xy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 34);

    auto ta2_xy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 35);

    auto ta2_xy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 36);

    auto ta2_xy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 37);

    auto ta2_xy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 38);

    auto ta2_xy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 39);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_1,  \
                             ta2_xy_0_xx_0,  \
                             ta2_xy_0_xx_1,  \
                             ta2_xy_0_xxx_0, \
                             ta2_xy_0_xxx_1, \
                             ta2_xy_0_xxy_0, \
                             ta2_xy_0_xxy_1, \
                             ta2_xy_0_xxz_0, \
                             ta2_xy_0_xxz_1, \
                             ta2_xy_0_xy_0,  \
                             ta2_xy_0_xy_1,  \
                             ta2_xy_0_xyy_0, \
                             ta2_xy_0_xyy_1, \
                             ta2_xy_0_xyz_0, \
                             ta2_xy_0_xyz_1, \
                             ta2_xy_0_xz_0,  \
                             ta2_xy_0_xz_1,  \
                             ta2_xy_0_xzz_0, \
                             ta2_xy_0_xzz_1, \
                             ta2_xy_0_yy_0,  \
                             ta2_xy_0_yy_1,  \
                             ta2_xy_0_yyy_0, \
                             ta2_xy_0_yyy_1, \
                             ta2_xy_0_yyz_0, \
                             ta2_xy_0_yyz_1, \
                             ta2_xy_0_yz_0,  \
                             ta2_xy_0_yz_1,  \
                             ta2_xy_0_yzz_0, \
                             ta2_xy_0_yzz_1, \
                             ta2_xy_0_zz_0,  \
                             ta2_xy_0_zz_1,  \
                             ta2_xy_0_zzz_0, \
                             ta2_xy_0_zzz_1, \
                             ta2_xy_x_xxx_0, \
                             ta2_xy_x_xxy_0, \
                             ta2_xy_x_xxz_0, \
                             ta2_xy_x_xyy_0, \
                             ta2_xy_x_xyz_0, \
                             ta2_xy_x_xzz_0, \
                             ta2_xy_x_yyy_0, \
                             ta2_xy_x_yyz_0, \
                             ta2_xy_x_yzz_0, \
                             ta2_xy_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_x_xxx_0[i] = 3.0 * ta2_xy_0_xx_0[i] * fe_0 - 3.0 * ta2_xy_0_xx_1[i] * fe_0 + ta1_y_0_xxx_1[i] + ta2_xy_0_xxx_0[i] * pa_x[i] -
                            ta2_xy_0_xxx_1[i] * pc_x[i];

        ta2_xy_x_xxy_0[i] = 2.0 * ta2_xy_0_xy_0[i] * fe_0 - 2.0 * ta2_xy_0_xy_1[i] * fe_0 + ta1_y_0_xxy_1[i] + ta2_xy_0_xxy_0[i] * pa_x[i] -
                            ta2_xy_0_xxy_1[i] * pc_x[i];

        ta2_xy_x_xxz_0[i] = 2.0 * ta2_xy_0_xz_0[i] * fe_0 - 2.0 * ta2_xy_0_xz_1[i] * fe_0 + ta1_y_0_xxz_1[i] + ta2_xy_0_xxz_0[i] * pa_x[i] -
                            ta2_xy_0_xxz_1[i] * pc_x[i];

        ta2_xy_x_xyy_0[i] =
            ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta1_y_0_xyy_1[i] + ta2_xy_0_xyy_0[i] * pa_x[i] - ta2_xy_0_xyy_1[i] * pc_x[i];

        ta2_xy_x_xyz_0[i] =
            ta2_xy_0_yz_0[i] * fe_0 - ta2_xy_0_yz_1[i] * fe_0 + ta1_y_0_xyz_1[i] + ta2_xy_0_xyz_0[i] * pa_x[i] - ta2_xy_0_xyz_1[i] * pc_x[i];

        ta2_xy_x_xzz_0[i] =
            ta2_xy_0_zz_0[i] * fe_0 - ta2_xy_0_zz_1[i] * fe_0 + ta1_y_0_xzz_1[i] + ta2_xy_0_xzz_0[i] * pa_x[i] - ta2_xy_0_xzz_1[i] * pc_x[i];

        ta2_xy_x_yyy_0[i] = ta1_y_0_yyy_1[i] + ta2_xy_0_yyy_0[i] * pa_x[i] - ta2_xy_0_yyy_1[i] * pc_x[i];

        ta2_xy_x_yyz_0[i] = ta1_y_0_yyz_1[i] + ta2_xy_0_yyz_0[i] * pa_x[i] - ta2_xy_0_yyz_1[i] * pc_x[i];

        ta2_xy_x_yzz_0[i] = ta1_y_0_yzz_1[i] + ta2_xy_0_yzz_0[i] * pa_x[i] - ta2_xy_0_yzz_1[i] * pc_x[i];

        ta2_xy_x_zzz_0[i] = ta1_y_0_zzz_1[i] + ta2_xy_0_zzz_0[i] * pa_x[i] - ta2_xy_0_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : PF

    auto ta2_xy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 40);

    auto ta2_xy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 41);

    auto ta2_xy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 42);

    auto ta2_xy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 43);

    auto ta2_xy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 44);

    auto ta2_xy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 45);

    auto ta2_xy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 46);

    auto ta2_xy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 47);

    auto ta2_xy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 48);

    auto ta2_xy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 49);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_1,  \
                             ta2_xy_0_xx_0,  \
                             ta2_xy_0_xx_1,  \
                             ta2_xy_0_xxx_0, \
                             ta2_xy_0_xxx_1, \
                             ta2_xy_0_xxy_0, \
                             ta2_xy_0_xxy_1, \
                             ta2_xy_0_xxz_0, \
                             ta2_xy_0_xxz_1, \
                             ta2_xy_0_xy_0,  \
                             ta2_xy_0_xy_1,  \
                             ta2_xy_0_xyy_0, \
                             ta2_xy_0_xyy_1, \
                             ta2_xy_0_xyz_0, \
                             ta2_xy_0_xyz_1, \
                             ta2_xy_0_xz_0,  \
                             ta2_xy_0_xz_1,  \
                             ta2_xy_0_xzz_0, \
                             ta2_xy_0_xzz_1, \
                             ta2_xy_0_yy_0,  \
                             ta2_xy_0_yy_1,  \
                             ta2_xy_0_yyy_0, \
                             ta2_xy_0_yyy_1, \
                             ta2_xy_0_yyz_0, \
                             ta2_xy_0_yyz_1, \
                             ta2_xy_0_yz_0,  \
                             ta2_xy_0_yz_1,  \
                             ta2_xy_0_yzz_0, \
                             ta2_xy_0_yzz_1, \
                             ta2_xy_0_zz_0,  \
                             ta2_xy_0_zz_1,  \
                             ta2_xy_0_zzz_0, \
                             ta2_xy_0_zzz_1, \
                             ta2_xy_y_xxx_0, \
                             ta2_xy_y_xxy_0, \
                             ta2_xy_y_xxz_0, \
                             ta2_xy_y_xyy_0, \
                             ta2_xy_y_xyz_0, \
                             ta2_xy_y_xzz_0, \
                             ta2_xy_y_yyy_0, \
                             ta2_xy_y_yyz_0, \
                             ta2_xy_y_yzz_0, \
                             ta2_xy_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_y_xxx_0[i] = ta1_x_0_xxx_1[i] + ta2_xy_0_xxx_0[i] * pa_y[i] - ta2_xy_0_xxx_1[i] * pc_y[i];

        ta2_xy_y_xxy_0[i] =
            ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + ta1_x_0_xxy_1[i] + ta2_xy_0_xxy_0[i] * pa_y[i] - ta2_xy_0_xxy_1[i] * pc_y[i];

        ta2_xy_y_xxz_0[i] = ta1_x_0_xxz_1[i] + ta2_xy_0_xxz_0[i] * pa_y[i] - ta2_xy_0_xxz_1[i] * pc_y[i];

        ta2_xy_y_xyy_0[i] = 2.0 * ta2_xy_0_xy_0[i] * fe_0 - 2.0 * ta2_xy_0_xy_1[i] * fe_0 + ta1_x_0_xyy_1[i] + ta2_xy_0_xyy_0[i] * pa_y[i] -
                            ta2_xy_0_xyy_1[i] * pc_y[i];

        ta2_xy_y_xyz_0[i] =
            ta2_xy_0_xz_0[i] * fe_0 - ta2_xy_0_xz_1[i] * fe_0 + ta1_x_0_xyz_1[i] + ta2_xy_0_xyz_0[i] * pa_y[i] - ta2_xy_0_xyz_1[i] * pc_y[i];

        ta2_xy_y_xzz_0[i] = ta1_x_0_xzz_1[i] + ta2_xy_0_xzz_0[i] * pa_y[i] - ta2_xy_0_xzz_1[i] * pc_y[i];

        ta2_xy_y_yyy_0[i] = 3.0 * ta2_xy_0_yy_0[i] * fe_0 - 3.0 * ta2_xy_0_yy_1[i] * fe_0 + ta1_x_0_yyy_1[i] + ta2_xy_0_yyy_0[i] * pa_y[i] -
                            ta2_xy_0_yyy_1[i] * pc_y[i];

        ta2_xy_y_yyz_0[i] = 2.0 * ta2_xy_0_yz_0[i] * fe_0 - 2.0 * ta2_xy_0_yz_1[i] * fe_0 + ta1_x_0_yyz_1[i] + ta2_xy_0_yyz_0[i] * pa_y[i] -
                            ta2_xy_0_yyz_1[i] * pc_y[i];

        ta2_xy_y_yzz_0[i] =
            ta2_xy_0_zz_0[i] * fe_0 - ta2_xy_0_zz_1[i] * fe_0 + ta1_x_0_yzz_1[i] + ta2_xy_0_yzz_0[i] * pa_y[i] - ta2_xy_0_yzz_1[i] * pc_y[i];

        ta2_xy_y_zzz_0[i] = ta1_x_0_zzz_1[i] + ta2_xy_0_zzz_0[i] * pa_y[i] - ta2_xy_0_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : PF

    auto ta2_xy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 50);

    auto ta2_xy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 51);

    auto ta2_xy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 52);

    auto ta2_xy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 53);

    auto ta2_xy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 54);

    auto ta2_xy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 55);

    auto ta2_xy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 56);

    auto ta2_xy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 57);

    auto ta2_xy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 58);

    auto ta2_xy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 59);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_xy_0_xx_0,  \
                             ta2_xy_0_xx_1,  \
                             ta2_xy_0_xxx_0, \
                             ta2_xy_0_xxx_1, \
                             ta2_xy_0_xxy_0, \
                             ta2_xy_0_xxy_1, \
                             ta2_xy_0_xxz_0, \
                             ta2_xy_0_xxz_1, \
                             ta2_xy_0_xy_0,  \
                             ta2_xy_0_xy_1,  \
                             ta2_xy_0_xyy_0, \
                             ta2_xy_0_xyy_1, \
                             ta2_xy_0_xyz_0, \
                             ta2_xy_0_xyz_1, \
                             ta2_xy_0_xz_0,  \
                             ta2_xy_0_xz_1,  \
                             ta2_xy_0_xzz_0, \
                             ta2_xy_0_xzz_1, \
                             ta2_xy_0_yy_0,  \
                             ta2_xy_0_yy_1,  \
                             ta2_xy_0_yyy_0, \
                             ta2_xy_0_yyy_1, \
                             ta2_xy_0_yyz_0, \
                             ta2_xy_0_yyz_1, \
                             ta2_xy_0_yz_0,  \
                             ta2_xy_0_yz_1,  \
                             ta2_xy_0_yzz_0, \
                             ta2_xy_0_yzz_1, \
                             ta2_xy_0_zz_0,  \
                             ta2_xy_0_zz_1,  \
                             ta2_xy_0_zzz_0, \
                             ta2_xy_0_zzz_1, \
                             ta2_xy_z_xxx_0, \
                             ta2_xy_z_xxy_0, \
                             ta2_xy_z_xxz_0, \
                             ta2_xy_z_xyy_0, \
                             ta2_xy_z_xyz_0, \
                             ta2_xy_z_xzz_0, \
                             ta2_xy_z_yyy_0, \
                             ta2_xy_z_yyz_0, \
                             ta2_xy_z_yzz_0, \
                             ta2_xy_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_z_xxx_0[i] = ta2_xy_0_xxx_0[i] * pa_z[i] - ta2_xy_0_xxx_1[i] * pc_z[i];

        ta2_xy_z_xxy_0[i] = ta2_xy_0_xxy_0[i] * pa_z[i] - ta2_xy_0_xxy_1[i] * pc_z[i];

        ta2_xy_z_xxz_0[i] = ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + ta2_xy_0_xxz_0[i] * pa_z[i] - ta2_xy_0_xxz_1[i] * pc_z[i];

        ta2_xy_z_xyy_0[i] = ta2_xy_0_xyy_0[i] * pa_z[i] - ta2_xy_0_xyy_1[i] * pc_z[i];

        ta2_xy_z_xyz_0[i] = ta2_xy_0_xy_0[i] * fe_0 - ta2_xy_0_xy_1[i] * fe_0 + ta2_xy_0_xyz_0[i] * pa_z[i] - ta2_xy_0_xyz_1[i] * pc_z[i];

        ta2_xy_z_xzz_0[i] = 2.0 * ta2_xy_0_xz_0[i] * fe_0 - 2.0 * ta2_xy_0_xz_1[i] * fe_0 + ta2_xy_0_xzz_0[i] * pa_z[i] - ta2_xy_0_xzz_1[i] * pc_z[i];

        ta2_xy_z_yyy_0[i] = ta2_xy_0_yyy_0[i] * pa_z[i] - ta2_xy_0_yyy_1[i] * pc_z[i];

        ta2_xy_z_yyz_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta2_xy_0_yyz_0[i] * pa_z[i] - ta2_xy_0_yyz_1[i] * pc_z[i];

        ta2_xy_z_yzz_0[i] = 2.0 * ta2_xy_0_yz_0[i] * fe_0 - 2.0 * ta2_xy_0_yz_1[i] * fe_0 + ta2_xy_0_yzz_0[i] * pa_z[i] - ta2_xy_0_yzz_1[i] * pc_z[i];

        ta2_xy_z_zzz_0[i] = 3.0 * ta2_xy_0_zz_0[i] * fe_0 - 3.0 * ta2_xy_0_zz_1[i] * fe_0 + ta2_xy_0_zzz_0[i] * pa_z[i] - ta2_xy_0_zzz_1[i] * pc_z[i];
    }

    // Set up 60-70 components of targeted buffer : PF

    auto ta2_xz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 60);

    auto ta2_xz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 61);

    auto ta2_xz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 62);

    auto ta2_xz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 63);

    auto ta2_xz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 64);

    auto ta2_xz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 65);

    auto ta2_xz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 66);

    auto ta2_xz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 67);

    auto ta2_xz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 68);

    auto ta2_xz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 69);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_1,  \
                             ta2_xz_0_xx_0,  \
                             ta2_xz_0_xx_1,  \
                             ta2_xz_0_xxx_0, \
                             ta2_xz_0_xxx_1, \
                             ta2_xz_0_xxy_0, \
                             ta2_xz_0_xxy_1, \
                             ta2_xz_0_xxz_0, \
                             ta2_xz_0_xxz_1, \
                             ta2_xz_0_xy_0,  \
                             ta2_xz_0_xy_1,  \
                             ta2_xz_0_xyy_0, \
                             ta2_xz_0_xyy_1, \
                             ta2_xz_0_xyz_0, \
                             ta2_xz_0_xyz_1, \
                             ta2_xz_0_xz_0,  \
                             ta2_xz_0_xz_1,  \
                             ta2_xz_0_xzz_0, \
                             ta2_xz_0_xzz_1, \
                             ta2_xz_0_yy_0,  \
                             ta2_xz_0_yy_1,  \
                             ta2_xz_0_yyy_0, \
                             ta2_xz_0_yyy_1, \
                             ta2_xz_0_yyz_0, \
                             ta2_xz_0_yyz_1, \
                             ta2_xz_0_yz_0,  \
                             ta2_xz_0_yz_1,  \
                             ta2_xz_0_yzz_0, \
                             ta2_xz_0_yzz_1, \
                             ta2_xz_0_zz_0,  \
                             ta2_xz_0_zz_1,  \
                             ta2_xz_0_zzz_0, \
                             ta2_xz_0_zzz_1, \
                             ta2_xz_x_xxx_0, \
                             ta2_xz_x_xxy_0, \
                             ta2_xz_x_xxz_0, \
                             ta2_xz_x_xyy_0, \
                             ta2_xz_x_xyz_0, \
                             ta2_xz_x_xzz_0, \
                             ta2_xz_x_yyy_0, \
                             ta2_xz_x_yyz_0, \
                             ta2_xz_x_yzz_0, \
                             ta2_xz_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_x_xxx_0[i] = 3.0 * ta2_xz_0_xx_0[i] * fe_0 - 3.0 * ta2_xz_0_xx_1[i] * fe_0 + ta1_z_0_xxx_1[i] + ta2_xz_0_xxx_0[i] * pa_x[i] -
                            ta2_xz_0_xxx_1[i] * pc_x[i];

        ta2_xz_x_xxy_0[i] = 2.0 * ta2_xz_0_xy_0[i] * fe_0 - 2.0 * ta2_xz_0_xy_1[i] * fe_0 + ta1_z_0_xxy_1[i] + ta2_xz_0_xxy_0[i] * pa_x[i] -
                            ta2_xz_0_xxy_1[i] * pc_x[i];

        ta2_xz_x_xxz_0[i] = 2.0 * ta2_xz_0_xz_0[i] * fe_0 - 2.0 * ta2_xz_0_xz_1[i] * fe_0 + ta1_z_0_xxz_1[i] + ta2_xz_0_xxz_0[i] * pa_x[i] -
                            ta2_xz_0_xxz_1[i] * pc_x[i];

        ta2_xz_x_xyy_0[i] =
            ta2_xz_0_yy_0[i] * fe_0 - ta2_xz_0_yy_1[i] * fe_0 + ta1_z_0_xyy_1[i] + ta2_xz_0_xyy_0[i] * pa_x[i] - ta2_xz_0_xyy_1[i] * pc_x[i];

        ta2_xz_x_xyz_0[i] =
            ta2_xz_0_yz_0[i] * fe_0 - ta2_xz_0_yz_1[i] * fe_0 + ta1_z_0_xyz_1[i] + ta2_xz_0_xyz_0[i] * pa_x[i] - ta2_xz_0_xyz_1[i] * pc_x[i];

        ta2_xz_x_xzz_0[i] =
            ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta1_z_0_xzz_1[i] + ta2_xz_0_xzz_0[i] * pa_x[i] - ta2_xz_0_xzz_1[i] * pc_x[i];

        ta2_xz_x_yyy_0[i] = ta1_z_0_yyy_1[i] + ta2_xz_0_yyy_0[i] * pa_x[i] - ta2_xz_0_yyy_1[i] * pc_x[i];

        ta2_xz_x_yyz_0[i] = ta1_z_0_yyz_1[i] + ta2_xz_0_yyz_0[i] * pa_x[i] - ta2_xz_0_yyz_1[i] * pc_x[i];

        ta2_xz_x_yzz_0[i] = ta1_z_0_yzz_1[i] + ta2_xz_0_yzz_0[i] * pa_x[i] - ta2_xz_0_yzz_1[i] * pc_x[i];

        ta2_xz_x_zzz_0[i] = ta1_z_0_zzz_1[i] + ta2_xz_0_zzz_0[i] * pa_x[i] - ta2_xz_0_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : PF

    auto ta2_xz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 70);

    auto ta2_xz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 71);

    auto ta2_xz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 72);

    auto ta2_xz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 73);

    auto ta2_xz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 74);

    auto ta2_xz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 75);

    auto ta2_xz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 76);

    auto ta2_xz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 77);

    auto ta2_xz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 78);

    auto ta2_xz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 79);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_xz_0_xx_0,  \
                             ta2_xz_0_xx_1,  \
                             ta2_xz_0_xxx_0, \
                             ta2_xz_0_xxx_1, \
                             ta2_xz_0_xxy_0, \
                             ta2_xz_0_xxy_1, \
                             ta2_xz_0_xxz_0, \
                             ta2_xz_0_xxz_1, \
                             ta2_xz_0_xy_0,  \
                             ta2_xz_0_xy_1,  \
                             ta2_xz_0_xyy_0, \
                             ta2_xz_0_xyy_1, \
                             ta2_xz_0_xyz_0, \
                             ta2_xz_0_xyz_1, \
                             ta2_xz_0_xz_0,  \
                             ta2_xz_0_xz_1,  \
                             ta2_xz_0_xzz_0, \
                             ta2_xz_0_xzz_1, \
                             ta2_xz_0_yy_0,  \
                             ta2_xz_0_yy_1,  \
                             ta2_xz_0_yyy_0, \
                             ta2_xz_0_yyy_1, \
                             ta2_xz_0_yyz_0, \
                             ta2_xz_0_yyz_1, \
                             ta2_xz_0_yz_0,  \
                             ta2_xz_0_yz_1,  \
                             ta2_xz_0_yzz_0, \
                             ta2_xz_0_yzz_1, \
                             ta2_xz_0_zz_0,  \
                             ta2_xz_0_zz_1,  \
                             ta2_xz_0_zzz_0, \
                             ta2_xz_0_zzz_1, \
                             ta2_xz_y_xxx_0, \
                             ta2_xz_y_xxy_0, \
                             ta2_xz_y_xxz_0, \
                             ta2_xz_y_xyy_0, \
                             ta2_xz_y_xyz_0, \
                             ta2_xz_y_xzz_0, \
                             ta2_xz_y_yyy_0, \
                             ta2_xz_y_yyz_0, \
                             ta2_xz_y_yzz_0, \
                             ta2_xz_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_y_xxx_0[i] = ta2_xz_0_xxx_0[i] * pa_y[i] - ta2_xz_0_xxx_1[i] * pc_y[i];

        ta2_xz_y_xxy_0[i] = ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + ta2_xz_0_xxy_0[i] * pa_y[i] - ta2_xz_0_xxy_1[i] * pc_y[i];

        ta2_xz_y_xxz_0[i] = ta2_xz_0_xxz_0[i] * pa_y[i] - ta2_xz_0_xxz_1[i] * pc_y[i];

        ta2_xz_y_xyy_0[i] = 2.0 * ta2_xz_0_xy_0[i] * fe_0 - 2.0 * ta2_xz_0_xy_1[i] * fe_0 + ta2_xz_0_xyy_0[i] * pa_y[i] - ta2_xz_0_xyy_1[i] * pc_y[i];

        ta2_xz_y_xyz_0[i] = ta2_xz_0_xz_0[i] * fe_0 - ta2_xz_0_xz_1[i] * fe_0 + ta2_xz_0_xyz_0[i] * pa_y[i] - ta2_xz_0_xyz_1[i] * pc_y[i];

        ta2_xz_y_xzz_0[i] = ta2_xz_0_xzz_0[i] * pa_y[i] - ta2_xz_0_xzz_1[i] * pc_y[i];

        ta2_xz_y_yyy_0[i] = 3.0 * ta2_xz_0_yy_0[i] * fe_0 - 3.0 * ta2_xz_0_yy_1[i] * fe_0 + ta2_xz_0_yyy_0[i] * pa_y[i] - ta2_xz_0_yyy_1[i] * pc_y[i];

        ta2_xz_y_yyz_0[i] = 2.0 * ta2_xz_0_yz_0[i] * fe_0 - 2.0 * ta2_xz_0_yz_1[i] * fe_0 + ta2_xz_0_yyz_0[i] * pa_y[i] - ta2_xz_0_yyz_1[i] * pc_y[i];

        ta2_xz_y_yzz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta2_xz_0_yzz_0[i] * pa_y[i] - ta2_xz_0_yzz_1[i] * pc_y[i];

        ta2_xz_y_zzz_0[i] = ta2_xz_0_zzz_0[i] * pa_y[i] - ta2_xz_0_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : PF

    auto ta2_xz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 80);

    auto ta2_xz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 81);

    auto ta2_xz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 82);

    auto ta2_xz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 83);

    auto ta2_xz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 84);

    auto ta2_xz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 85);

    auto ta2_xz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 86);

    auto ta2_xz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 87);

    auto ta2_xz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 88);

    auto ta2_xz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 89);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_1,  \
                             ta2_xz_0_xx_0,  \
                             ta2_xz_0_xx_1,  \
                             ta2_xz_0_xxx_0, \
                             ta2_xz_0_xxx_1, \
                             ta2_xz_0_xxy_0, \
                             ta2_xz_0_xxy_1, \
                             ta2_xz_0_xxz_0, \
                             ta2_xz_0_xxz_1, \
                             ta2_xz_0_xy_0,  \
                             ta2_xz_0_xy_1,  \
                             ta2_xz_0_xyy_0, \
                             ta2_xz_0_xyy_1, \
                             ta2_xz_0_xyz_0, \
                             ta2_xz_0_xyz_1, \
                             ta2_xz_0_xz_0,  \
                             ta2_xz_0_xz_1,  \
                             ta2_xz_0_xzz_0, \
                             ta2_xz_0_xzz_1, \
                             ta2_xz_0_yy_0,  \
                             ta2_xz_0_yy_1,  \
                             ta2_xz_0_yyy_0, \
                             ta2_xz_0_yyy_1, \
                             ta2_xz_0_yyz_0, \
                             ta2_xz_0_yyz_1, \
                             ta2_xz_0_yz_0,  \
                             ta2_xz_0_yz_1,  \
                             ta2_xz_0_yzz_0, \
                             ta2_xz_0_yzz_1, \
                             ta2_xz_0_zz_0,  \
                             ta2_xz_0_zz_1,  \
                             ta2_xz_0_zzz_0, \
                             ta2_xz_0_zzz_1, \
                             ta2_xz_z_xxx_0, \
                             ta2_xz_z_xxy_0, \
                             ta2_xz_z_xxz_0, \
                             ta2_xz_z_xyy_0, \
                             ta2_xz_z_xyz_0, \
                             ta2_xz_z_xzz_0, \
                             ta2_xz_z_yyy_0, \
                             ta2_xz_z_yyz_0, \
                             ta2_xz_z_yzz_0, \
                             ta2_xz_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_z_xxx_0[i] = ta1_x_0_xxx_1[i] + ta2_xz_0_xxx_0[i] * pa_z[i] - ta2_xz_0_xxx_1[i] * pc_z[i];

        ta2_xz_z_xxy_0[i] = ta1_x_0_xxy_1[i] + ta2_xz_0_xxy_0[i] * pa_z[i] - ta2_xz_0_xxy_1[i] * pc_z[i];

        ta2_xz_z_xxz_0[i] =
            ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + ta1_x_0_xxz_1[i] + ta2_xz_0_xxz_0[i] * pa_z[i] - ta2_xz_0_xxz_1[i] * pc_z[i];

        ta2_xz_z_xyy_0[i] = ta1_x_0_xyy_1[i] + ta2_xz_0_xyy_0[i] * pa_z[i] - ta2_xz_0_xyy_1[i] * pc_z[i];

        ta2_xz_z_xyz_0[i] =
            ta2_xz_0_xy_0[i] * fe_0 - ta2_xz_0_xy_1[i] * fe_0 + ta1_x_0_xyz_1[i] + ta2_xz_0_xyz_0[i] * pa_z[i] - ta2_xz_0_xyz_1[i] * pc_z[i];

        ta2_xz_z_xzz_0[i] = 2.0 * ta2_xz_0_xz_0[i] * fe_0 - 2.0 * ta2_xz_0_xz_1[i] * fe_0 + ta1_x_0_xzz_1[i] + ta2_xz_0_xzz_0[i] * pa_z[i] -
                            ta2_xz_0_xzz_1[i] * pc_z[i];

        ta2_xz_z_yyy_0[i] = ta1_x_0_yyy_1[i] + ta2_xz_0_yyy_0[i] * pa_z[i] - ta2_xz_0_yyy_1[i] * pc_z[i];

        ta2_xz_z_yyz_0[i] =
            ta2_xz_0_yy_0[i] * fe_0 - ta2_xz_0_yy_1[i] * fe_0 + ta1_x_0_yyz_1[i] + ta2_xz_0_yyz_0[i] * pa_z[i] - ta2_xz_0_yyz_1[i] * pc_z[i];

        ta2_xz_z_yzz_0[i] = 2.0 * ta2_xz_0_yz_0[i] * fe_0 - 2.0 * ta2_xz_0_yz_1[i] * fe_0 + ta1_x_0_yzz_1[i] + ta2_xz_0_yzz_0[i] * pa_z[i] -
                            ta2_xz_0_yzz_1[i] * pc_z[i];

        ta2_xz_z_zzz_0[i] = 3.0 * ta2_xz_0_zz_0[i] * fe_0 - 3.0 * ta2_xz_0_zz_1[i] * fe_0 + ta1_x_0_zzz_1[i] + ta2_xz_0_zzz_0[i] * pa_z[i] -
                            ta2_xz_0_zzz_1[i] * pc_z[i];
    }

    // Set up 90-100 components of targeted buffer : PF

    auto ta2_yy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 90);

    auto ta2_yy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 91);

    auto ta2_yy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 92);

    auto ta2_yy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 93);

    auto ta2_yy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 94);

    auto ta2_yy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 95);

    auto ta2_yy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 96);

    auto ta2_yy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 97);

    auto ta2_yy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 98);

    auto ta2_yy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 99);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yy_0_xx_0,  \
                             ta2_yy_0_xx_1,  \
                             ta2_yy_0_xxx_0, \
                             ta2_yy_0_xxx_1, \
                             ta2_yy_0_xxy_0, \
                             ta2_yy_0_xxy_1, \
                             ta2_yy_0_xxz_0, \
                             ta2_yy_0_xxz_1, \
                             ta2_yy_0_xy_0,  \
                             ta2_yy_0_xy_1,  \
                             ta2_yy_0_xyy_0, \
                             ta2_yy_0_xyy_1, \
                             ta2_yy_0_xyz_0, \
                             ta2_yy_0_xyz_1, \
                             ta2_yy_0_xz_0,  \
                             ta2_yy_0_xz_1,  \
                             ta2_yy_0_xzz_0, \
                             ta2_yy_0_xzz_1, \
                             ta2_yy_0_yy_0,  \
                             ta2_yy_0_yy_1,  \
                             ta2_yy_0_yyy_0, \
                             ta2_yy_0_yyy_1, \
                             ta2_yy_0_yyz_0, \
                             ta2_yy_0_yyz_1, \
                             ta2_yy_0_yz_0,  \
                             ta2_yy_0_yz_1,  \
                             ta2_yy_0_yzz_0, \
                             ta2_yy_0_yzz_1, \
                             ta2_yy_0_zz_0,  \
                             ta2_yy_0_zz_1,  \
                             ta2_yy_0_zzz_0, \
                             ta2_yy_0_zzz_1, \
                             ta2_yy_x_xxx_0, \
                             ta2_yy_x_xxy_0, \
                             ta2_yy_x_xxz_0, \
                             ta2_yy_x_xyy_0, \
                             ta2_yy_x_xyz_0, \
                             ta2_yy_x_xzz_0, \
                             ta2_yy_x_yyy_0, \
                             ta2_yy_x_yyz_0, \
                             ta2_yy_x_yzz_0, \
                             ta2_yy_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_x_xxx_0[i] = 3.0 * ta2_yy_0_xx_0[i] * fe_0 - 3.0 * ta2_yy_0_xx_1[i] * fe_0 + ta2_yy_0_xxx_0[i] * pa_x[i] - ta2_yy_0_xxx_1[i] * pc_x[i];

        ta2_yy_x_xxy_0[i] = 2.0 * ta2_yy_0_xy_0[i] * fe_0 - 2.0 * ta2_yy_0_xy_1[i] * fe_0 + ta2_yy_0_xxy_0[i] * pa_x[i] - ta2_yy_0_xxy_1[i] * pc_x[i];

        ta2_yy_x_xxz_0[i] = 2.0 * ta2_yy_0_xz_0[i] * fe_0 - 2.0 * ta2_yy_0_xz_1[i] * fe_0 + ta2_yy_0_xxz_0[i] * pa_x[i] - ta2_yy_0_xxz_1[i] * pc_x[i];

        ta2_yy_x_xyy_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_0_xyy_0[i] * pa_x[i] - ta2_yy_0_xyy_1[i] * pc_x[i];

        ta2_yy_x_xyz_0[i] = ta2_yy_0_yz_0[i] * fe_0 - ta2_yy_0_yz_1[i] * fe_0 + ta2_yy_0_xyz_0[i] * pa_x[i] - ta2_yy_0_xyz_1[i] * pc_x[i];

        ta2_yy_x_xzz_0[i] = ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + ta2_yy_0_xzz_0[i] * pa_x[i] - ta2_yy_0_xzz_1[i] * pc_x[i];

        ta2_yy_x_yyy_0[i] = ta2_yy_0_yyy_0[i] * pa_x[i] - ta2_yy_0_yyy_1[i] * pc_x[i];

        ta2_yy_x_yyz_0[i] = ta2_yy_0_yyz_0[i] * pa_x[i] - ta2_yy_0_yyz_1[i] * pc_x[i];

        ta2_yy_x_yzz_0[i] = ta2_yy_0_yzz_0[i] * pa_x[i] - ta2_yy_0_yzz_1[i] * pc_x[i];

        ta2_yy_x_zzz_0[i] = ta2_yy_0_zzz_0[i] * pa_x[i] - ta2_yy_0_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : PF

    auto ta2_yy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 100);

    auto ta2_yy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 101);

    auto ta2_yy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 102);

    auto ta2_yy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 103);

    auto ta2_yy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 104);

    auto ta2_yy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 105);

    auto ta2_yy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 106);

    auto ta2_yy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 107);

    auto ta2_yy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 108);

    auto ta2_yy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 109);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_1,  \
                             ta2_yy_0_xx_0,  \
                             ta2_yy_0_xx_1,  \
                             ta2_yy_0_xxx_0, \
                             ta2_yy_0_xxx_1, \
                             ta2_yy_0_xxy_0, \
                             ta2_yy_0_xxy_1, \
                             ta2_yy_0_xxz_0, \
                             ta2_yy_0_xxz_1, \
                             ta2_yy_0_xy_0,  \
                             ta2_yy_0_xy_1,  \
                             ta2_yy_0_xyy_0, \
                             ta2_yy_0_xyy_1, \
                             ta2_yy_0_xyz_0, \
                             ta2_yy_0_xyz_1, \
                             ta2_yy_0_xz_0,  \
                             ta2_yy_0_xz_1,  \
                             ta2_yy_0_xzz_0, \
                             ta2_yy_0_xzz_1, \
                             ta2_yy_0_yy_0,  \
                             ta2_yy_0_yy_1,  \
                             ta2_yy_0_yyy_0, \
                             ta2_yy_0_yyy_1, \
                             ta2_yy_0_yyz_0, \
                             ta2_yy_0_yyz_1, \
                             ta2_yy_0_yz_0,  \
                             ta2_yy_0_yz_1,  \
                             ta2_yy_0_yzz_0, \
                             ta2_yy_0_yzz_1, \
                             ta2_yy_0_zz_0,  \
                             ta2_yy_0_zz_1,  \
                             ta2_yy_0_zzz_0, \
                             ta2_yy_0_zzz_1, \
                             ta2_yy_y_xxx_0, \
                             ta2_yy_y_xxy_0, \
                             ta2_yy_y_xxz_0, \
                             ta2_yy_y_xyy_0, \
                             ta2_yy_y_xyz_0, \
                             ta2_yy_y_xzz_0, \
                             ta2_yy_y_yyy_0, \
                             ta2_yy_y_yyz_0, \
                             ta2_yy_y_yzz_0, \
                             ta2_yy_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_y_xxx_0[i] = 2.0 * ta1_y_0_xxx_1[i] + ta2_yy_0_xxx_0[i] * pa_y[i] - ta2_yy_0_xxx_1[i] * pc_y[i];

        ta2_yy_y_xxy_0[i] =
            ta2_yy_0_xx_0[i] * fe_0 - ta2_yy_0_xx_1[i] * fe_0 + 2.0 * ta1_y_0_xxy_1[i] + ta2_yy_0_xxy_0[i] * pa_y[i] - ta2_yy_0_xxy_1[i] * pc_y[i];

        ta2_yy_y_xxz_0[i] = 2.0 * ta1_y_0_xxz_1[i] + ta2_yy_0_xxz_0[i] * pa_y[i] - ta2_yy_0_xxz_1[i] * pc_y[i];

        ta2_yy_y_xyy_0[i] = 2.0 * ta2_yy_0_xy_0[i] * fe_0 - 2.0 * ta2_yy_0_xy_1[i] * fe_0 + 2.0 * ta1_y_0_xyy_1[i] + ta2_yy_0_xyy_0[i] * pa_y[i] -
                            ta2_yy_0_xyy_1[i] * pc_y[i];

        ta2_yy_y_xyz_0[i] =
            ta2_yy_0_xz_0[i] * fe_0 - ta2_yy_0_xz_1[i] * fe_0 + 2.0 * ta1_y_0_xyz_1[i] + ta2_yy_0_xyz_0[i] * pa_y[i] - ta2_yy_0_xyz_1[i] * pc_y[i];

        ta2_yy_y_xzz_0[i] = 2.0 * ta1_y_0_xzz_1[i] + ta2_yy_0_xzz_0[i] * pa_y[i] - ta2_yy_0_xzz_1[i] * pc_y[i];

        ta2_yy_y_yyy_0[i] = 3.0 * ta2_yy_0_yy_0[i] * fe_0 - 3.0 * ta2_yy_0_yy_1[i] * fe_0 + 2.0 * ta1_y_0_yyy_1[i] + ta2_yy_0_yyy_0[i] * pa_y[i] -
                            ta2_yy_0_yyy_1[i] * pc_y[i];

        ta2_yy_y_yyz_0[i] = 2.0 * ta2_yy_0_yz_0[i] * fe_0 - 2.0 * ta2_yy_0_yz_1[i] * fe_0 + 2.0 * ta1_y_0_yyz_1[i] + ta2_yy_0_yyz_0[i] * pa_y[i] -
                            ta2_yy_0_yyz_1[i] * pc_y[i];

        ta2_yy_y_yzz_0[i] =
            ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + 2.0 * ta1_y_0_yzz_1[i] + ta2_yy_0_yzz_0[i] * pa_y[i] - ta2_yy_0_yzz_1[i] * pc_y[i];

        ta2_yy_y_zzz_0[i] = 2.0 * ta1_y_0_zzz_1[i] + ta2_yy_0_zzz_0[i] * pa_y[i] - ta2_yy_0_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : PF

    auto ta2_yy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 110);

    auto ta2_yy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 111);

    auto ta2_yy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 112);

    auto ta2_yy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 113);

    auto ta2_yy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 114);

    auto ta2_yy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 115);

    auto ta2_yy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 116);

    auto ta2_yy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 117);

    auto ta2_yy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 118);

    auto ta2_yy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 119);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta2_yy_0_xx_0,  \
                             ta2_yy_0_xx_1,  \
                             ta2_yy_0_xxx_0, \
                             ta2_yy_0_xxx_1, \
                             ta2_yy_0_xxy_0, \
                             ta2_yy_0_xxy_1, \
                             ta2_yy_0_xxz_0, \
                             ta2_yy_0_xxz_1, \
                             ta2_yy_0_xy_0,  \
                             ta2_yy_0_xy_1,  \
                             ta2_yy_0_xyy_0, \
                             ta2_yy_0_xyy_1, \
                             ta2_yy_0_xyz_0, \
                             ta2_yy_0_xyz_1, \
                             ta2_yy_0_xz_0,  \
                             ta2_yy_0_xz_1,  \
                             ta2_yy_0_xzz_0, \
                             ta2_yy_0_xzz_1, \
                             ta2_yy_0_yy_0,  \
                             ta2_yy_0_yy_1,  \
                             ta2_yy_0_yyy_0, \
                             ta2_yy_0_yyy_1, \
                             ta2_yy_0_yyz_0, \
                             ta2_yy_0_yyz_1, \
                             ta2_yy_0_yz_0,  \
                             ta2_yy_0_yz_1,  \
                             ta2_yy_0_yzz_0, \
                             ta2_yy_0_yzz_1, \
                             ta2_yy_0_zz_0,  \
                             ta2_yy_0_zz_1,  \
                             ta2_yy_0_zzz_0, \
                             ta2_yy_0_zzz_1, \
                             ta2_yy_z_xxx_0, \
                             ta2_yy_z_xxy_0, \
                             ta2_yy_z_xxz_0, \
                             ta2_yy_z_xyy_0, \
                             ta2_yy_z_xyz_0, \
                             ta2_yy_z_xzz_0, \
                             ta2_yy_z_yyy_0, \
                             ta2_yy_z_yyz_0, \
                             ta2_yy_z_yzz_0, \
                             ta2_yy_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_z_xxx_0[i] = ta2_yy_0_xxx_0[i] * pa_z[i] - ta2_yy_0_xxx_1[i] * pc_z[i];

        ta2_yy_z_xxy_0[i] = ta2_yy_0_xxy_0[i] * pa_z[i] - ta2_yy_0_xxy_1[i] * pc_z[i];

        ta2_yy_z_xxz_0[i] = ta2_yy_0_xx_0[i] * fe_0 - ta2_yy_0_xx_1[i] * fe_0 + ta2_yy_0_xxz_0[i] * pa_z[i] - ta2_yy_0_xxz_1[i] * pc_z[i];

        ta2_yy_z_xyy_0[i] = ta2_yy_0_xyy_0[i] * pa_z[i] - ta2_yy_0_xyy_1[i] * pc_z[i];

        ta2_yy_z_xyz_0[i] = ta2_yy_0_xy_0[i] * fe_0 - ta2_yy_0_xy_1[i] * fe_0 + ta2_yy_0_xyz_0[i] * pa_z[i] - ta2_yy_0_xyz_1[i] * pc_z[i];

        ta2_yy_z_xzz_0[i] = 2.0 * ta2_yy_0_xz_0[i] * fe_0 - 2.0 * ta2_yy_0_xz_1[i] * fe_0 + ta2_yy_0_xzz_0[i] * pa_z[i] - ta2_yy_0_xzz_1[i] * pc_z[i];

        ta2_yy_z_yyy_0[i] = ta2_yy_0_yyy_0[i] * pa_z[i] - ta2_yy_0_yyy_1[i] * pc_z[i];

        ta2_yy_z_yyz_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_0_yyz_0[i] * pa_z[i] - ta2_yy_0_yyz_1[i] * pc_z[i];

        ta2_yy_z_yzz_0[i] = 2.0 * ta2_yy_0_yz_0[i] * fe_0 - 2.0 * ta2_yy_0_yz_1[i] * fe_0 + ta2_yy_0_yzz_0[i] * pa_z[i] - ta2_yy_0_yzz_1[i] * pc_z[i];

        ta2_yy_z_zzz_0[i] = 3.0 * ta2_yy_0_zz_0[i] * fe_0 - 3.0 * ta2_yy_0_zz_1[i] * fe_0 + ta2_yy_0_zzz_0[i] * pa_z[i] - ta2_yy_0_zzz_1[i] * pc_z[i];
    }

    // Set up 120-130 components of targeted buffer : PF

    auto ta2_yz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 120);

    auto ta2_yz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 121);

    auto ta2_yz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 122);

    auto ta2_yz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 123);

    auto ta2_yz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 124);

    auto ta2_yz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 125);

    auto ta2_yz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 126);

    auto ta2_yz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 127);

    auto ta2_yz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 128);

    auto ta2_yz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 129);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_yz_0_xx_0,  \
                             ta2_yz_0_xx_1,  \
                             ta2_yz_0_xxx_0, \
                             ta2_yz_0_xxx_1, \
                             ta2_yz_0_xxy_0, \
                             ta2_yz_0_xxy_1, \
                             ta2_yz_0_xxz_0, \
                             ta2_yz_0_xxz_1, \
                             ta2_yz_0_xy_0,  \
                             ta2_yz_0_xy_1,  \
                             ta2_yz_0_xyy_0, \
                             ta2_yz_0_xyy_1, \
                             ta2_yz_0_xyz_0, \
                             ta2_yz_0_xyz_1, \
                             ta2_yz_0_xz_0,  \
                             ta2_yz_0_xz_1,  \
                             ta2_yz_0_xzz_0, \
                             ta2_yz_0_xzz_1, \
                             ta2_yz_0_yy_0,  \
                             ta2_yz_0_yy_1,  \
                             ta2_yz_0_yyy_0, \
                             ta2_yz_0_yyy_1, \
                             ta2_yz_0_yyz_0, \
                             ta2_yz_0_yyz_1, \
                             ta2_yz_0_yz_0,  \
                             ta2_yz_0_yz_1,  \
                             ta2_yz_0_yzz_0, \
                             ta2_yz_0_yzz_1, \
                             ta2_yz_0_zz_0,  \
                             ta2_yz_0_zz_1,  \
                             ta2_yz_0_zzz_0, \
                             ta2_yz_0_zzz_1, \
                             ta2_yz_x_xxx_0, \
                             ta2_yz_x_xxy_0, \
                             ta2_yz_x_xxz_0, \
                             ta2_yz_x_xyy_0, \
                             ta2_yz_x_xyz_0, \
                             ta2_yz_x_xzz_0, \
                             ta2_yz_x_yyy_0, \
                             ta2_yz_x_yyz_0, \
                             ta2_yz_x_yzz_0, \
                             ta2_yz_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_x_xxx_0[i] = 3.0 * ta2_yz_0_xx_0[i] * fe_0 - 3.0 * ta2_yz_0_xx_1[i] * fe_0 + ta2_yz_0_xxx_0[i] * pa_x[i] - ta2_yz_0_xxx_1[i] * pc_x[i];

        ta2_yz_x_xxy_0[i] = 2.0 * ta2_yz_0_xy_0[i] * fe_0 - 2.0 * ta2_yz_0_xy_1[i] * fe_0 + ta2_yz_0_xxy_0[i] * pa_x[i] - ta2_yz_0_xxy_1[i] * pc_x[i];

        ta2_yz_x_xxz_0[i] = 2.0 * ta2_yz_0_xz_0[i] * fe_0 - 2.0 * ta2_yz_0_xz_1[i] * fe_0 + ta2_yz_0_xxz_0[i] * pa_x[i] - ta2_yz_0_xxz_1[i] * pc_x[i];

        ta2_yz_x_xyy_0[i] = ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + ta2_yz_0_xyy_0[i] * pa_x[i] - ta2_yz_0_xyy_1[i] * pc_x[i];

        ta2_yz_x_xyz_0[i] = ta2_yz_0_yz_0[i] * fe_0 - ta2_yz_0_yz_1[i] * fe_0 + ta2_yz_0_xyz_0[i] * pa_x[i] - ta2_yz_0_xyz_1[i] * pc_x[i];

        ta2_yz_x_xzz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta2_yz_0_xzz_0[i] * pa_x[i] - ta2_yz_0_xzz_1[i] * pc_x[i];

        ta2_yz_x_yyy_0[i] = ta2_yz_0_yyy_0[i] * pa_x[i] - ta2_yz_0_yyy_1[i] * pc_x[i];

        ta2_yz_x_yyz_0[i] = ta2_yz_0_yyz_0[i] * pa_x[i] - ta2_yz_0_yyz_1[i] * pc_x[i];

        ta2_yz_x_yzz_0[i] = ta2_yz_0_yzz_0[i] * pa_x[i] - ta2_yz_0_yzz_1[i] * pc_x[i];

        ta2_yz_x_zzz_0[i] = ta2_yz_0_zzz_0[i] * pa_x[i] - ta2_yz_0_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : PF

    auto ta2_yz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 130);

    auto ta2_yz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 131);

    auto ta2_yz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 132);

    auto ta2_yz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 133);

    auto ta2_yz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 134);

    auto ta2_yz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 135);

    auto ta2_yz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 136);

    auto ta2_yz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 137);

    auto ta2_yz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 138);

    auto ta2_yz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 139);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_1,  \
                             ta2_yz_0_xx_0,  \
                             ta2_yz_0_xx_1,  \
                             ta2_yz_0_xxx_0, \
                             ta2_yz_0_xxx_1, \
                             ta2_yz_0_xxy_0, \
                             ta2_yz_0_xxy_1, \
                             ta2_yz_0_xxz_0, \
                             ta2_yz_0_xxz_1, \
                             ta2_yz_0_xy_0,  \
                             ta2_yz_0_xy_1,  \
                             ta2_yz_0_xyy_0, \
                             ta2_yz_0_xyy_1, \
                             ta2_yz_0_xyz_0, \
                             ta2_yz_0_xyz_1, \
                             ta2_yz_0_xz_0,  \
                             ta2_yz_0_xz_1,  \
                             ta2_yz_0_xzz_0, \
                             ta2_yz_0_xzz_1, \
                             ta2_yz_0_yy_0,  \
                             ta2_yz_0_yy_1,  \
                             ta2_yz_0_yyy_0, \
                             ta2_yz_0_yyy_1, \
                             ta2_yz_0_yyz_0, \
                             ta2_yz_0_yyz_1, \
                             ta2_yz_0_yz_0,  \
                             ta2_yz_0_yz_1,  \
                             ta2_yz_0_yzz_0, \
                             ta2_yz_0_yzz_1, \
                             ta2_yz_0_zz_0,  \
                             ta2_yz_0_zz_1,  \
                             ta2_yz_0_zzz_0, \
                             ta2_yz_0_zzz_1, \
                             ta2_yz_y_xxx_0, \
                             ta2_yz_y_xxy_0, \
                             ta2_yz_y_xxz_0, \
                             ta2_yz_y_xyy_0, \
                             ta2_yz_y_xyz_0, \
                             ta2_yz_y_xzz_0, \
                             ta2_yz_y_yyy_0, \
                             ta2_yz_y_yyz_0, \
                             ta2_yz_y_yzz_0, \
                             ta2_yz_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_y_xxx_0[i] = ta1_z_0_xxx_1[i] + ta2_yz_0_xxx_0[i] * pa_y[i] - ta2_yz_0_xxx_1[i] * pc_y[i];

        ta2_yz_y_xxy_0[i] =
            ta2_yz_0_xx_0[i] * fe_0 - ta2_yz_0_xx_1[i] * fe_0 + ta1_z_0_xxy_1[i] + ta2_yz_0_xxy_0[i] * pa_y[i] - ta2_yz_0_xxy_1[i] * pc_y[i];

        ta2_yz_y_xxz_0[i] = ta1_z_0_xxz_1[i] + ta2_yz_0_xxz_0[i] * pa_y[i] - ta2_yz_0_xxz_1[i] * pc_y[i];

        ta2_yz_y_xyy_0[i] = 2.0 * ta2_yz_0_xy_0[i] * fe_0 - 2.0 * ta2_yz_0_xy_1[i] * fe_0 + ta1_z_0_xyy_1[i] + ta2_yz_0_xyy_0[i] * pa_y[i] -
                            ta2_yz_0_xyy_1[i] * pc_y[i];

        ta2_yz_y_xyz_0[i] =
            ta2_yz_0_xz_0[i] * fe_0 - ta2_yz_0_xz_1[i] * fe_0 + ta1_z_0_xyz_1[i] + ta2_yz_0_xyz_0[i] * pa_y[i] - ta2_yz_0_xyz_1[i] * pc_y[i];

        ta2_yz_y_xzz_0[i] = ta1_z_0_xzz_1[i] + ta2_yz_0_xzz_0[i] * pa_y[i] - ta2_yz_0_xzz_1[i] * pc_y[i];

        ta2_yz_y_yyy_0[i] = 3.0 * ta2_yz_0_yy_0[i] * fe_0 - 3.0 * ta2_yz_0_yy_1[i] * fe_0 + ta1_z_0_yyy_1[i] + ta2_yz_0_yyy_0[i] * pa_y[i] -
                            ta2_yz_0_yyy_1[i] * pc_y[i];

        ta2_yz_y_yyz_0[i] = 2.0 * ta2_yz_0_yz_0[i] * fe_0 - 2.0 * ta2_yz_0_yz_1[i] * fe_0 + ta1_z_0_yyz_1[i] + ta2_yz_0_yyz_0[i] * pa_y[i] -
                            ta2_yz_0_yyz_1[i] * pc_y[i];

        ta2_yz_y_yzz_0[i] =
            ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta1_z_0_yzz_1[i] + ta2_yz_0_yzz_0[i] * pa_y[i] - ta2_yz_0_yzz_1[i] * pc_y[i];

        ta2_yz_y_zzz_0[i] = ta1_z_0_zzz_1[i] + ta2_yz_0_zzz_0[i] * pa_y[i] - ta2_yz_0_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : PF

    auto ta2_yz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 140);

    auto ta2_yz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 141);

    auto ta2_yz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 142);

    auto ta2_yz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 143);

    auto ta2_yz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 144);

    auto ta2_yz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 145);

    auto ta2_yz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 146);

    auto ta2_yz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 147);

    auto ta2_yz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 148);

    auto ta2_yz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 149);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_1,  \
                             ta2_yz_0_xx_0,  \
                             ta2_yz_0_xx_1,  \
                             ta2_yz_0_xxx_0, \
                             ta2_yz_0_xxx_1, \
                             ta2_yz_0_xxy_0, \
                             ta2_yz_0_xxy_1, \
                             ta2_yz_0_xxz_0, \
                             ta2_yz_0_xxz_1, \
                             ta2_yz_0_xy_0,  \
                             ta2_yz_0_xy_1,  \
                             ta2_yz_0_xyy_0, \
                             ta2_yz_0_xyy_1, \
                             ta2_yz_0_xyz_0, \
                             ta2_yz_0_xyz_1, \
                             ta2_yz_0_xz_0,  \
                             ta2_yz_0_xz_1,  \
                             ta2_yz_0_xzz_0, \
                             ta2_yz_0_xzz_1, \
                             ta2_yz_0_yy_0,  \
                             ta2_yz_0_yy_1,  \
                             ta2_yz_0_yyy_0, \
                             ta2_yz_0_yyy_1, \
                             ta2_yz_0_yyz_0, \
                             ta2_yz_0_yyz_1, \
                             ta2_yz_0_yz_0,  \
                             ta2_yz_0_yz_1,  \
                             ta2_yz_0_yzz_0, \
                             ta2_yz_0_yzz_1, \
                             ta2_yz_0_zz_0,  \
                             ta2_yz_0_zz_1,  \
                             ta2_yz_0_zzz_0, \
                             ta2_yz_0_zzz_1, \
                             ta2_yz_z_xxx_0, \
                             ta2_yz_z_xxy_0, \
                             ta2_yz_z_xxz_0, \
                             ta2_yz_z_xyy_0, \
                             ta2_yz_z_xyz_0, \
                             ta2_yz_z_xzz_0, \
                             ta2_yz_z_yyy_0, \
                             ta2_yz_z_yyz_0, \
                             ta2_yz_z_yzz_0, \
                             ta2_yz_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_z_xxx_0[i] = ta1_y_0_xxx_1[i] + ta2_yz_0_xxx_0[i] * pa_z[i] - ta2_yz_0_xxx_1[i] * pc_z[i];

        ta2_yz_z_xxy_0[i] = ta1_y_0_xxy_1[i] + ta2_yz_0_xxy_0[i] * pa_z[i] - ta2_yz_0_xxy_1[i] * pc_z[i];

        ta2_yz_z_xxz_0[i] =
            ta2_yz_0_xx_0[i] * fe_0 - ta2_yz_0_xx_1[i] * fe_0 + ta1_y_0_xxz_1[i] + ta2_yz_0_xxz_0[i] * pa_z[i] - ta2_yz_0_xxz_1[i] * pc_z[i];

        ta2_yz_z_xyy_0[i] = ta1_y_0_xyy_1[i] + ta2_yz_0_xyy_0[i] * pa_z[i] - ta2_yz_0_xyy_1[i] * pc_z[i];

        ta2_yz_z_xyz_0[i] =
            ta2_yz_0_xy_0[i] * fe_0 - ta2_yz_0_xy_1[i] * fe_0 + ta1_y_0_xyz_1[i] + ta2_yz_0_xyz_0[i] * pa_z[i] - ta2_yz_0_xyz_1[i] * pc_z[i];

        ta2_yz_z_xzz_0[i] = 2.0 * ta2_yz_0_xz_0[i] * fe_0 - 2.0 * ta2_yz_0_xz_1[i] * fe_0 + ta1_y_0_xzz_1[i] + ta2_yz_0_xzz_0[i] * pa_z[i] -
                            ta2_yz_0_xzz_1[i] * pc_z[i];

        ta2_yz_z_yyy_0[i] = ta1_y_0_yyy_1[i] + ta2_yz_0_yyy_0[i] * pa_z[i] - ta2_yz_0_yyy_1[i] * pc_z[i];

        ta2_yz_z_yyz_0[i] =
            ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + ta1_y_0_yyz_1[i] + ta2_yz_0_yyz_0[i] * pa_z[i] - ta2_yz_0_yyz_1[i] * pc_z[i];

        ta2_yz_z_yzz_0[i] = 2.0 * ta2_yz_0_yz_0[i] * fe_0 - 2.0 * ta2_yz_0_yz_1[i] * fe_0 + ta1_y_0_yzz_1[i] + ta2_yz_0_yzz_0[i] * pa_z[i] -
                            ta2_yz_0_yzz_1[i] * pc_z[i];

        ta2_yz_z_zzz_0[i] = 3.0 * ta2_yz_0_zz_0[i] * fe_0 - 3.0 * ta2_yz_0_zz_1[i] * fe_0 + ta1_y_0_zzz_1[i] + ta2_yz_0_zzz_0[i] * pa_z[i] -
                            ta2_yz_0_zzz_1[i] * pc_z[i];
    }

    // Set up 150-160 components of targeted buffer : PF

    auto ta2_zz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 150);

    auto ta2_zz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 151);

    auto ta2_zz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 152);

    auto ta2_zz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 153);

    auto ta2_zz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 154);

    auto ta2_zz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 155);

    auto ta2_zz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 156);

    auto ta2_zz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 157);

    auto ta2_zz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 158);

    auto ta2_zz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 159);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta2_zz_0_xx_0,  \
                             ta2_zz_0_xx_1,  \
                             ta2_zz_0_xxx_0, \
                             ta2_zz_0_xxx_1, \
                             ta2_zz_0_xxy_0, \
                             ta2_zz_0_xxy_1, \
                             ta2_zz_0_xxz_0, \
                             ta2_zz_0_xxz_1, \
                             ta2_zz_0_xy_0,  \
                             ta2_zz_0_xy_1,  \
                             ta2_zz_0_xyy_0, \
                             ta2_zz_0_xyy_1, \
                             ta2_zz_0_xyz_0, \
                             ta2_zz_0_xyz_1, \
                             ta2_zz_0_xz_0,  \
                             ta2_zz_0_xz_1,  \
                             ta2_zz_0_xzz_0, \
                             ta2_zz_0_xzz_1, \
                             ta2_zz_0_yy_0,  \
                             ta2_zz_0_yy_1,  \
                             ta2_zz_0_yyy_0, \
                             ta2_zz_0_yyy_1, \
                             ta2_zz_0_yyz_0, \
                             ta2_zz_0_yyz_1, \
                             ta2_zz_0_yz_0,  \
                             ta2_zz_0_yz_1,  \
                             ta2_zz_0_yzz_0, \
                             ta2_zz_0_yzz_1, \
                             ta2_zz_0_zz_0,  \
                             ta2_zz_0_zz_1,  \
                             ta2_zz_0_zzz_0, \
                             ta2_zz_0_zzz_1, \
                             ta2_zz_x_xxx_0, \
                             ta2_zz_x_xxy_0, \
                             ta2_zz_x_xxz_0, \
                             ta2_zz_x_xyy_0, \
                             ta2_zz_x_xyz_0, \
                             ta2_zz_x_xzz_0, \
                             ta2_zz_x_yyy_0, \
                             ta2_zz_x_yyz_0, \
                             ta2_zz_x_yzz_0, \
                             ta2_zz_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_x_xxx_0[i] = 3.0 * ta2_zz_0_xx_0[i] * fe_0 - 3.0 * ta2_zz_0_xx_1[i] * fe_0 + ta2_zz_0_xxx_0[i] * pa_x[i] - ta2_zz_0_xxx_1[i] * pc_x[i];

        ta2_zz_x_xxy_0[i] = 2.0 * ta2_zz_0_xy_0[i] * fe_0 - 2.0 * ta2_zz_0_xy_1[i] * fe_0 + ta2_zz_0_xxy_0[i] * pa_x[i] - ta2_zz_0_xxy_1[i] * pc_x[i];

        ta2_zz_x_xxz_0[i] = 2.0 * ta2_zz_0_xz_0[i] * fe_0 - 2.0 * ta2_zz_0_xz_1[i] * fe_0 + ta2_zz_0_xxz_0[i] * pa_x[i] - ta2_zz_0_xxz_1[i] * pc_x[i];

        ta2_zz_x_xyy_0[i] = ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + ta2_zz_0_xyy_0[i] * pa_x[i] - ta2_zz_0_xyy_1[i] * pc_x[i];

        ta2_zz_x_xyz_0[i] = ta2_zz_0_yz_0[i] * fe_0 - ta2_zz_0_yz_1[i] * fe_0 + ta2_zz_0_xyz_0[i] * pa_x[i] - ta2_zz_0_xyz_1[i] * pc_x[i];

        ta2_zz_x_xzz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_0_xzz_0[i] * pa_x[i] - ta2_zz_0_xzz_1[i] * pc_x[i];

        ta2_zz_x_yyy_0[i] = ta2_zz_0_yyy_0[i] * pa_x[i] - ta2_zz_0_yyy_1[i] * pc_x[i];

        ta2_zz_x_yyz_0[i] = ta2_zz_0_yyz_0[i] * pa_x[i] - ta2_zz_0_yyz_1[i] * pc_x[i];

        ta2_zz_x_yzz_0[i] = ta2_zz_0_yzz_0[i] * pa_x[i] - ta2_zz_0_yzz_1[i] * pc_x[i];

        ta2_zz_x_zzz_0[i] = ta2_zz_0_zzz_0[i] * pa_x[i] - ta2_zz_0_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : PF

    auto ta2_zz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 160);

    auto ta2_zz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 161);

    auto ta2_zz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 162);

    auto ta2_zz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 163);

    auto ta2_zz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 164);

    auto ta2_zz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 165);

    auto ta2_zz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 166);

    auto ta2_zz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 167);

    auto ta2_zz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 168);

    auto ta2_zz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 169);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta2_zz_0_xx_0,  \
                             ta2_zz_0_xx_1,  \
                             ta2_zz_0_xxx_0, \
                             ta2_zz_0_xxx_1, \
                             ta2_zz_0_xxy_0, \
                             ta2_zz_0_xxy_1, \
                             ta2_zz_0_xxz_0, \
                             ta2_zz_0_xxz_1, \
                             ta2_zz_0_xy_0,  \
                             ta2_zz_0_xy_1,  \
                             ta2_zz_0_xyy_0, \
                             ta2_zz_0_xyy_1, \
                             ta2_zz_0_xyz_0, \
                             ta2_zz_0_xyz_1, \
                             ta2_zz_0_xz_0,  \
                             ta2_zz_0_xz_1,  \
                             ta2_zz_0_xzz_0, \
                             ta2_zz_0_xzz_1, \
                             ta2_zz_0_yy_0,  \
                             ta2_zz_0_yy_1,  \
                             ta2_zz_0_yyy_0, \
                             ta2_zz_0_yyy_1, \
                             ta2_zz_0_yyz_0, \
                             ta2_zz_0_yyz_1, \
                             ta2_zz_0_yz_0,  \
                             ta2_zz_0_yz_1,  \
                             ta2_zz_0_yzz_0, \
                             ta2_zz_0_yzz_1, \
                             ta2_zz_0_zz_0,  \
                             ta2_zz_0_zz_1,  \
                             ta2_zz_0_zzz_0, \
                             ta2_zz_0_zzz_1, \
                             ta2_zz_y_xxx_0, \
                             ta2_zz_y_xxy_0, \
                             ta2_zz_y_xxz_0, \
                             ta2_zz_y_xyy_0, \
                             ta2_zz_y_xyz_0, \
                             ta2_zz_y_xzz_0, \
                             ta2_zz_y_yyy_0, \
                             ta2_zz_y_yyz_0, \
                             ta2_zz_y_yzz_0, \
                             ta2_zz_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_y_xxx_0[i] = ta2_zz_0_xxx_0[i] * pa_y[i] - ta2_zz_0_xxx_1[i] * pc_y[i];

        ta2_zz_y_xxy_0[i] = ta2_zz_0_xx_0[i] * fe_0 - ta2_zz_0_xx_1[i] * fe_0 + ta2_zz_0_xxy_0[i] * pa_y[i] - ta2_zz_0_xxy_1[i] * pc_y[i];

        ta2_zz_y_xxz_0[i] = ta2_zz_0_xxz_0[i] * pa_y[i] - ta2_zz_0_xxz_1[i] * pc_y[i];

        ta2_zz_y_xyy_0[i] = 2.0 * ta2_zz_0_xy_0[i] * fe_0 - 2.0 * ta2_zz_0_xy_1[i] * fe_0 + ta2_zz_0_xyy_0[i] * pa_y[i] - ta2_zz_0_xyy_1[i] * pc_y[i];

        ta2_zz_y_xyz_0[i] = ta2_zz_0_xz_0[i] * fe_0 - ta2_zz_0_xz_1[i] * fe_0 + ta2_zz_0_xyz_0[i] * pa_y[i] - ta2_zz_0_xyz_1[i] * pc_y[i];

        ta2_zz_y_xzz_0[i] = ta2_zz_0_xzz_0[i] * pa_y[i] - ta2_zz_0_xzz_1[i] * pc_y[i];

        ta2_zz_y_yyy_0[i] = 3.0 * ta2_zz_0_yy_0[i] * fe_0 - 3.0 * ta2_zz_0_yy_1[i] * fe_0 + ta2_zz_0_yyy_0[i] * pa_y[i] - ta2_zz_0_yyy_1[i] * pc_y[i];

        ta2_zz_y_yyz_0[i] = 2.0 * ta2_zz_0_yz_0[i] * fe_0 - 2.0 * ta2_zz_0_yz_1[i] * fe_0 + ta2_zz_0_yyz_0[i] * pa_y[i] - ta2_zz_0_yyz_1[i] * pc_y[i];

        ta2_zz_y_yzz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_0_yzz_0[i] * pa_y[i] - ta2_zz_0_yzz_1[i] * pc_y[i];

        ta2_zz_y_zzz_0[i] = ta2_zz_0_zzz_0[i] * pa_y[i] - ta2_zz_0_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : PF

    auto ta2_zz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 170);

    auto ta2_zz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 171);

    auto ta2_zz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 172);

    auto ta2_zz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 173);

    auto ta2_zz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 174);

    auto ta2_zz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 175);

    auto ta2_zz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 176);

    auto ta2_zz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 177);

    auto ta2_zz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 178);

    auto ta2_zz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 179);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_1,  \
                             ta2_zz_0_xx_0,  \
                             ta2_zz_0_xx_1,  \
                             ta2_zz_0_xxx_0, \
                             ta2_zz_0_xxx_1, \
                             ta2_zz_0_xxy_0, \
                             ta2_zz_0_xxy_1, \
                             ta2_zz_0_xxz_0, \
                             ta2_zz_0_xxz_1, \
                             ta2_zz_0_xy_0,  \
                             ta2_zz_0_xy_1,  \
                             ta2_zz_0_xyy_0, \
                             ta2_zz_0_xyy_1, \
                             ta2_zz_0_xyz_0, \
                             ta2_zz_0_xyz_1, \
                             ta2_zz_0_xz_0,  \
                             ta2_zz_0_xz_1,  \
                             ta2_zz_0_xzz_0, \
                             ta2_zz_0_xzz_1, \
                             ta2_zz_0_yy_0,  \
                             ta2_zz_0_yy_1,  \
                             ta2_zz_0_yyy_0, \
                             ta2_zz_0_yyy_1, \
                             ta2_zz_0_yyz_0, \
                             ta2_zz_0_yyz_1, \
                             ta2_zz_0_yz_0,  \
                             ta2_zz_0_yz_1,  \
                             ta2_zz_0_yzz_0, \
                             ta2_zz_0_yzz_1, \
                             ta2_zz_0_zz_0,  \
                             ta2_zz_0_zz_1,  \
                             ta2_zz_0_zzz_0, \
                             ta2_zz_0_zzz_1, \
                             ta2_zz_z_xxx_0, \
                             ta2_zz_z_xxy_0, \
                             ta2_zz_z_xxz_0, \
                             ta2_zz_z_xyy_0, \
                             ta2_zz_z_xyz_0, \
                             ta2_zz_z_xzz_0, \
                             ta2_zz_z_yyy_0, \
                             ta2_zz_z_yyz_0, \
                             ta2_zz_z_yzz_0, \
                             ta2_zz_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_z_xxx_0[i] = 2.0 * ta1_z_0_xxx_1[i] + ta2_zz_0_xxx_0[i] * pa_z[i] - ta2_zz_0_xxx_1[i] * pc_z[i];

        ta2_zz_z_xxy_0[i] = 2.0 * ta1_z_0_xxy_1[i] + ta2_zz_0_xxy_0[i] * pa_z[i] - ta2_zz_0_xxy_1[i] * pc_z[i];

        ta2_zz_z_xxz_0[i] =
            ta2_zz_0_xx_0[i] * fe_0 - ta2_zz_0_xx_1[i] * fe_0 + 2.0 * ta1_z_0_xxz_1[i] + ta2_zz_0_xxz_0[i] * pa_z[i] - ta2_zz_0_xxz_1[i] * pc_z[i];

        ta2_zz_z_xyy_0[i] = 2.0 * ta1_z_0_xyy_1[i] + ta2_zz_0_xyy_0[i] * pa_z[i] - ta2_zz_0_xyy_1[i] * pc_z[i];

        ta2_zz_z_xyz_0[i] =
            ta2_zz_0_xy_0[i] * fe_0 - ta2_zz_0_xy_1[i] * fe_0 + 2.0 * ta1_z_0_xyz_1[i] + ta2_zz_0_xyz_0[i] * pa_z[i] - ta2_zz_0_xyz_1[i] * pc_z[i];

        ta2_zz_z_xzz_0[i] = 2.0 * ta2_zz_0_xz_0[i] * fe_0 - 2.0 * ta2_zz_0_xz_1[i] * fe_0 + 2.0 * ta1_z_0_xzz_1[i] + ta2_zz_0_xzz_0[i] * pa_z[i] -
                            ta2_zz_0_xzz_1[i] * pc_z[i];

        ta2_zz_z_yyy_0[i] = 2.0 * ta1_z_0_yyy_1[i] + ta2_zz_0_yyy_0[i] * pa_z[i] - ta2_zz_0_yyy_1[i] * pc_z[i];

        ta2_zz_z_yyz_0[i] =
            ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + 2.0 * ta1_z_0_yyz_1[i] + ta2_zz_0_yyz_0[i] * pa_z[i] - ta2_zz_0_yyz_1[i] * pc_z[i];

        ta2_zz_z_yzz_0[i] = 2.0 * ta2_zz_0_yz_0[i] * fe_0 - 2.0 * ta2_zz_0_yz_1[i] * fe_0 + 2.0 * ta1_z_0_yzz_1[i] + ta2_zz_0_yzz_0[i] * pa_z[i] -
                            ta2_zz_0_yzz_1[i] * pc_z[i];

        ta2_zz_z_zzz_0[i] = 3.0 * ta2_zz_0_zz_0[i] * fe_0 - 3.0 * ta2_zz_0_zz_1[i] * fe_0 + 2.0 * ta1_z_0_zzz_1[i] + ta2_zz_0_zzz_0[i] * pa_z[i] -
                            ta2_zz_0_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
