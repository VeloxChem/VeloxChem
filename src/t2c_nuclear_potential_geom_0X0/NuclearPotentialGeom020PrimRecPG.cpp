#include "NuclearPotentialGeom020PrimRecPG.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_pg(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_pg,
                                        const size_t idx_npot_geom_020_0_sf,
                                        const size_t idx_npot_geom_020_1_sf,
                                        const size_t idx_npot_geom_010_1_sg,
                                        const size_t idx_npot_geom_020_0_sg,
                                        const size_t idx_npot_geom_020_1_sg,
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

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg);

    auto ta1_x_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 1);

    auto ta1_x_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 2);

    auto ta1_x_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 3);

    auto ta1_x_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 4);

    auto ta1_x_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 5);

    auto ta1_x_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 6);

    auto ta1_x_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 7);

    auto ta1_x_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 8);

    auto ta1_x_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 9);

    auto ta1_x_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 10);

    auto ta1_x_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 11);

    auto ta1_x_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 12);

    auto ta1_x_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 13);

    auto ta1_x_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 14);

    auto ta1_y_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 15);

    auto ta1_y_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 16);

    auto ta1_y_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 17);

    auto ta1_y_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 18);

    auto ta1_y_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 19);

    auto ta1_y_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 20);

    auto ta1_y_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 21);

    auto ta1_y_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 22);

    auto ta1_y_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 23);

    auto ta1_y_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 24);

    auto ta1_y_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 25);

    auto ta1_y_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 26);

    auto ta1_y_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 27);

    auto ta1_y_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 28);

    auto ta1_y_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 29);

    auto ta1_z_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 30);

    auto ta1_z_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 31);

    auto ta1_z_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 32);

    auto ta1_z_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 33);

    auto ta1_z_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 34);

    auto ta1_z_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 35);

    auto ta1_z_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 36);

    auto ta1_z_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 37);

    auto ta1_z_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 38);

    auto ta1_z_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 39);

    auto ta1_z_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 40);

    auto ta1_z_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 41);

    auto ta1_z_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 42);

    auto ta1_z_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 43);

    auto ta1_z_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 44);

    // Set up components of auxiliary buffer : SG

    auto ta2_xx_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg);

    auto ta2_xx_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 1);

    auto ta2_xx_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 2);

    auto ta2_xx_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 3);

    auto ta2_xx_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 4);

    auto ta2_xx_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 5);

    auto ta2_xx_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 6);

    auto ta2_xx_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 7);

    auto ta2_xx_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 8);

    auto ta2_xx_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 9);

    auto ta2_xx_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 10);

    auto ta2_xx_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 11);

    auto ta2_xx_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 12);

    auto ta2_xx_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 13);

    auto ta2_xx_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 14);

    auto ta2_xy_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 15);

    auto ta2_xy_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 16);

    auto ta2_xy_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 17);

    auto ta2_xy_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 18);

    auto ta2_xy_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 19);

    auto ta2_xy_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 20);

    auto ta2_xy_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 21);

    auto ta2_xy_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 22);

    auto ta2_xy_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 23);

    auto ta2_xy_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 24);

    auto ta2_xy_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 25);

    auto ta2_xy_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 26);

    auto ta2_xy_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 27);

    auto ta2_xy_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 28);

    auto ta2_xy_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 29);

    auto ta2_xz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 30);

    auto ta2_xz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 31);

    auto ta2_xz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 32);

    auto ta2_xz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 33);

    auto ta2_xz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 34);

    auto ta2_xz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 35);

    auto ta2_xz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 36);

    auto ta2_xz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 37);

    auto ta2_xz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 38);

    auto ta2_xz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 39);

    auto ta2_xz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 40);

    auto ta2_xz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 41);

    auto ta2_xz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 42);

    auto ta2_xz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 43);

    auto ta2_xz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 44);

    auto ta2_yy_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 45);

    auto ta2_yy_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 46);

    auto ta2_yy_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 47);

    auto ta2_yy_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 48);

    auto ta2_yy_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 49);

    auto ta2_yy_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 50);

    auto ta2_yy_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 51);

    auto ta2_yy_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 52);

    auto ta2_yy_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 53);

    auto ta2_yy_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 54);

    auto ta2_yy_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 55);

    auto ta2_yy_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 56);

    auto ta2_yy_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 57);

    auto ta2_yy_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 58);

    auto ta2_yy_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 59);

    auto ta2_yz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 60);

    auto ta2_yz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 61);

    auto ta2_yz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 62);

    auto ta2_yz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 63);

    auto ta2_yz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 64);

    auto ta2_yz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 65);

    auto ta2_yz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 66);

    auto ta2_yz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 67);

    auto ta2_yz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 68);

    auto ta2_yz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 69);

    auto ta2_yz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 70);

    auto ta2_yz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 71);

    auto ta2_yz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 72);

    auto ta2_yz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 73);

    auto ta2_yz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 74);

    auto ta2_zz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 75);

    auto ta2_zz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 76);

    auto ta2_zz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 77);

    auto ta2_zz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 78);

    auto ta2_zz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 79);

    auto ta2_zz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 80);

    auto ta2_zz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 81);

    auto ta2_zz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 82);

    auto ta2_zz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 83);

    auto ta2_zz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 84);

    auto ta2_zz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 85);

    auto ta2_zz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 86);

    auto ta2_zz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 87);

    auto ta2_zz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 88);

    auto ta2_zz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 89);

    // Set up components of auxiliary buffer : SG

    auto ta2_xx_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg);

    auto ta2_xx_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 1);

    auto ta2_xx_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 2);

    auto ta2_xx_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 3);

    auto ta2_xx_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 4);

    auto ta2_xx_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 5);

    auto ta2_xx_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 6);

    auto ta2_xx_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 7);

    auto ta2_xx_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 8);

    auto ta2_xx_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 9);

    auto ta2_xx_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 10);

    auto ta2_xx_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 11);

    auto ta2_xx_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 12);

    auto ta2_xx_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 13);

    auto ta2_xx_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 14);

    auto ta2_xy_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 15);

    auto ta2_xy_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 16);

    auto ta2_xy_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 17);

    auto ta2_xy_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 18);

    auto ta2_xy_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 19);

    auto ta2_xy_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 20);

    auto ta2_xy_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 21);

    auto ta2_xy_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 22);

    auto ta2_xy_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 23);

    auto ta2_xy_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 24);

    auto ta2_xy_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 25);

    auto ta2_xy_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 26);

    auto ta2_xy_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 27);

    auto ta2_xy_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 28);

    auto ta2_xy_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 29);

    auto ta2_xz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 30);

    auto ta2_xz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 31);

    auto ta2_xz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 32);

    auto ta2_xz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 33);

    auto ta2_xz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 34);

    auto ta2_xz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 35);

    auto ta2_xz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 36);

    auto ta2_xz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 37);

    auto ta2_xz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 38);

    auto ta2_xz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 39);

    auto ta2_xz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 40);

    auto ta2_xz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 41);

    auto ta2_xz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 42);

    auto ta2_xz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 43);

    auto ta2_xz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 44);

    auto ta2_yy_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 45);

    auto ta2_yy_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 46);

    auto ta2_yy_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 47);

    auto ta2_yy_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 48);

    auto ta2_yy_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 49);

    auto ta2_yy_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 50);

    auto ta2_yy_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 51);

    auto ta2_yy_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 52);

    auto ta2_yy_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 53);

    auto ta2_yy_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 54);

    auto ta2_yy_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 55);

    auto ta2_yy_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 56);

    auto ta2_yy_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 57);

    auto ta2_yy_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 58);

    auto ta2_yy_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 59);

    auto ta2_yz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 60);

    auto ta2_yz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 61);

    auto ta2_yz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 62);

    auto ta2_yz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 63);

    auto ta2_yz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 64);

    auto ta2_yz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 65);

    auto ta2_yz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 66);

    auto ta2_yz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 67);

    auto ta2_yz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 68);

    auto ta2_yz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 69);

    auto ta2_yz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 70);

    auto ta2_yz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 71);

    auto ta2_yz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 72);

    auto ta2_yz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 73);

    auto ta2_yz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 74);

    auto ta2_zz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 75);

    auto ta2_zz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 76);

    auto ta2_zz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 77);

    auto ta2_zz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 78);

    auto ta2_zz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 79);

    auto ta2_zz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 80);

    auto ta2_zz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 81);

    auto ta2_zz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 82);

    auto ta2_zz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 83);

    auto ta2_zz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 84);

    auto ta2_zz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 85);

    auto ta2_zz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 86);

    auto ta2_zz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 87);

    auto ta2_zz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 88);

    auto ta2_zz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 89);

    // Set up 0-15 components of targeted buffer : PG

    auto ta2_xx_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg);

    auto ta2_xx_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 1);

    auto ta2_xx_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 2);

    auto ta2_xx_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 3);

    auto ta2_xx_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 4);

    auto ta2_xx_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 5);

    auto ta2_xx_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 6);

    auto ta2_xx_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 7);

    auto ta2_xx_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 8);

    auto ta2_xx_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 9);

    auto ta2_xx_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 10);

    auto ta2_xx_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 11);

    auto ta2_xx_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 12);

    auto ta2_xx_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 13);

    auto ta2_xx_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 14);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_0_xxxx_1, ta1_x_0_xxxy_1, ta1_x_0_xxxz_1, ta1_x_0_xxyy_1, ta1_x_0_xxyz_1, ta1_x_0_xxzz_1, ta1_x_0_xyyy_1, ta1_x_0_xyyz_1, ta1_x_0_xyzz_1, ta1_x_0_xzzz_1, ta1_x_0_yyyy_1, ta1_x_0_yyyz_1, ta1_x_0_yyzz_1, ta1_x_0_yzzz_1, ta1_x_0_zzzz_1, ta2_xx_0_xxx_0, ta2_xx_0_xxx_1, ta2_xx_0_xxxx_0, ta2_xx_0_xxxx_1, ta2_xx_0_xxxy_0, ta2_xx_0_xxxy_1, ta2_xx_0_xxxz_0, ta2_xx_0_xxxz_1, ta2_xx_0_xxy_0, ta2_xx_0_xxy_1, ta2_xx_0_xxyy_0, ta2_xx_0_xxyy_1, ta2_xx_0_xxyz_0, ta2_xx_0_xxyz_1, ta2_xx_0_xxz_0, ta2_xx_0_xxz_1, ta2_xx_0_xxzz_0, ta2_xx_0_xxzz_1, ta2_xx_0_xyy_0, ta2_xx_0_xyy_1, ta2_xx_0_xyyy_0, ta2_xx_0_xyyy_1, ta2_xx_0_xyyz_0, ta2_xx_0_xyyz_1, ta2_xx_0_xyz_0, ta2_xx_0_xyz_1, ta2_xx_0_xyzz_0, ta2_xx_0_xyzz_1, ta2_xx_0_xzz_0, ta2_xx_0_xzz_1, ta2_xx_0_xzzz_0, ta2_xx_0_xzzz_1, ta2_xx_0_yyy_0, ta2_xx_0_yyy_1, ta2_xx_0_yyyy_0, ta2_xx_0_yyyy_1, ta2_xx_0_yyyz_0, ta2_xx_0_yyyz_1, ta2_xx_0_yyz_0, ta2_xx_0_yyz_1, ta2_xx_0_yyzz_0, ta2_xx_0_yyzz_1, ta2_xx_0_yzz_0, ta2_xx_0_yzz_1, ta2_xx_0_yzzz_0, ta2_xx_0_yzzz_1, ta2_xx_0_zzz_0, ta2_xx_0_zzz_1, ta2_xx_0_zzzz_0, ta2_xx_0_zzzz_1, ta2_xx_x_xxxx_0, ta2_xx_x_xxxy_0, ta2_xx_x_xxxz_0, ta2_xx_x_xxyy_0, ta2_xx_x_xxyz_0, ta2_xx_x_xxzz_0, ta2_xx_x_xyyy_0, ta2_xx_x_xyyz_0, ta2_xx_x_xyzz_0, ta2_xx_x_xzzz_0, ta2_xx_x_yyyy_0, ta2_xx_x_yyyz_0, ta2_xx_x_yyzz_0, ta2_xx_x_yzzz_0, ta2_xx_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_x_xxxx_0[i] = 4.0 * ta2_xx_0_xxx_0[i] * fe_0 - 4.0 * ta2_xx_0_xxx_1[i] * fe_0 + 2.0 * ta1_x_0_xxxx_1[i] + ta2_xx_0_xxxx_0[i] * pa_x[i] - ta2_xx_0_xxxx_1[i] * pc_x[i];

        ta2_xx_x_xxxy_0[i] = 3.0 * ta2_xx_0_xxy_0[i] * fe_0 - 3.0 * ta2_xx_0_xxy_1[i] * fe_0 + 2.0 * ta1_x_0_xxxy_1[i] + ta2_xx_0_xxxy_0[i] * pa_x[i] - ta2_xx_0_xxxy_1[i] * pc_x[i];

        ta2_xx_x_xxxz_0[i] = 3.0 * ta2_xx_0_xxz_0[i] * fe_0 - 3.0 * ta2_xx_0_xxz_1[i] * fe_0 + 2.0 * ta1_x_0_xxxz_1[i] + ta2_xx_0_xxxz_0[i] * pa_x[i] - ta2_xx_0_xxxz_1[i] * pc_x[i];

        ta2_xx_x_xxyy_0[i] = 2.0 * ta2_xx_0_xyy_0[i] * fe_0 - 2.0 * ta2_xx_0_xyy_1[i] * fe_0 + 2.0 * ta1_x_0_xxyy_1[i] + ta2_xx_0_xxyy_0[i] * pa_x[i] - ta2_xx_0_xxyy_1[i] * pc_x[i];

        ta2_xx_x_xxyz_0[i] = 2.0 * ta2_xx_0_xyz_0[i] * fe_0 - 2.0 * ta2_xx_0_xyz_1[i] * fe_0 + 2.0 * ta1_x_0_xxyz_1[i] + ta2_xx_0_xxyz_0[i] * pa_x[i] - ta2_xx_0_xxyz_1[i] * pc_x[i];

        ta2_xx_x_xxzz_0[i] = 2.0 * ta2_xx_0_xzz_0[i] * fe_0 - 2.0 * ta2_xx_0_xzz_1[i] * fe_0 + 2.0 * ta1_x_0_xxzz_1[i] + ta2_xx_0_xxzz_0[i] * pa_x[i] - ta2_xx_0_xxzz_1[i] * pc_x[i];

        ta2_xx_x_xyyy_0[i] = ta2_xx_0_yyy_0[i] * fe_0 - ta2_xx_0_yyy_1[i] * fe_0 + 2.0 * ta1_x_0_xyyy_1[i] + ta2_xx_0_xyyy_0[i] * pa_x[i] - ta2_xx_0_xyyy_1[i] * pc_x[i];

        ta2_xx_x_xyyz_0[i] = ta2_xx_0_yyz_0[i] * fe_0 - ta2_xx_0_yyz_1[i] * fe_0 + 2.0 * ta1_x_0_xyyz_1[i] + ta2_xx_0_xyyz_0[i] * pa_x[i] - ta2_xx_0_xyyz_1[i] * pc_x[i];

        ta2_xx_x_xyzz_0[i] = ta2_xx_0_yzz_0[i] * fe_0 - ta2_xx_0_yzz_1[i] * fe_0 + 2.0 * ta1_x_0_xyzz_1[i] + ta2_xx_0_xyzz_0[i] * pa_x[i] - ta2_xx_0_xyzz_1[i] * pc_x[i];

        ta2_xx_x_xzzz_0[i] = ta2_xx_0_zzz_0[i] * fe_0 - ta2_xx_0_zzz_1[i] * fe_0 + 2.0 * ta1_x_0_xzzz_1[i] + ta2_xx_0_xzzz_0[i] * pa_x[i] - ta2_xx_0_xzzz_1[i] * pc_x[i];

        ta2_xx_x_yyyy_0[i] = 2.0 * ta1_x_0_yyyy_1[i] + ta2_xx_0_yyyy_0[i] * pa_x[i] - ta2_xx_0_yyyy_1[i] * pc_x[i];

        ta2_xx_x_yyyz_0[i] = 2.0 * ta1_x_0_yyyz_1[i] + ta2_xx_0_yyyz_0[i] * pa_x[i] - ta2_xx_0_yyyz_1[i] * pc_x[i];

        ta2_xx_x_yyzz_0[i] = 2.0 * ta1_x_0_yyzz_1[i] + ta2_xx_0_yyzz_0[i] * pa_x[i] - ta2_xx_0_yyzz_1[i] * pc_x[i];

        ta2_xx_x_yzzz_0[i] = 2.0 * ta1_x_0_yzzz_1[i] + ta2_xx_0_yzzz_0[i] * pa_x[i] - ta2_xx_0_yzzz_1[i] * pc_x[i];

        ta2_xx_x_zzzz_0[i] = 2.0 * ta1_x_0_zzzz_1[i] + ta2_xx_0_zzzz_0[i] * pa_x[i] - ta2_xx_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto ta2_xx_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 15);

    auto ta2_xx_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 16);

    auto ta2_xx_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 17);

    auto ta2_xx_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 18);

    auto ta2_xx_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 19);

    auto ta2_xx_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 20);

    auto ta2_xx_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 21);

    auto ta2_xx_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 22);

    auto ta2_xx_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 23);

    auto ta2_xx_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 24);

    auto ta2_xx_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 25);

    auto ta2_xx_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 26);

    auto ta2_xx_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 27);

    auto ta2_xx_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 28);

    auto ta2_xx_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 29);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_0_xxx_0, ta2_xx_0_xxx_1, ta2_xx_0_xxxx_0, ta2_xx_0_xxxx_1, ta2_xx_0_xxxy_0, ta2_xx_0_xxxy_1, ta2_xx_0_xxxz_0, ta2_xx_0_xxxz_1, ta2_xx_0_xxy_0, ta2_xx_0_xxy_1, ta2_xx_0_xxyy_0, ta2_xx_0_xxyy_1, ta2_xx_0_xxyz_0, ta2_xx_0_xxyz_1, ta2_xx_0_xxz_0, ta2_xx_0_xxz_1, ta2_xx_0_xxzz_0, ta2_xx_0_xxzz_1, ta2_xx_0_xyy_0, ta2_xx_0_xyy_1, ta2_xx_0_xyyy_0, ta2_xx_0_xyyy_1, ta2_xx_0_xyyz_0, ta2_xx_0_xyyz_1, ta2_xx_0_xyz_0, ta2_xx_0_xyz_1, ta2_xx_0_xyzz_0, ta2_xx_0_xyzz_1, ta2_xx_0_xzz_0, ta2_xx_0_xzz_1, ta2_xx_0_xzzz_0, ta2_xx_0_xzzz_1, ta2_xx_0_yyy_0, ta2_xx_0_yyy_1, ta2_xx_0_yyyy_0, ta2_xx_0_yyyy_1, ta2_xx_0_yyyz_0, ta2_xx_0_yyyz_1, ta2_xx_0_yyz_0, ta2_xx_0_yyz_1, ta2_xx_0_yyzz_0, ta2_xx_0_yyzz_1, ta2_xx_0_yzz_0, ta2_xx_0_yzz_1, ta2_xx_0_yzzz_0, ta2_xx_0_yzzz_1, ta2_xx_0_zzz_0, ta2_xx_0_zzz_1, ta2_xx_0_zzzz_0, ta2_xx_0_zzzz_1, ta2_xx_y_xxxx_0, ta2_xx_y_xxxy_0, ta2_xx_y_xxxz_0, ta2_xx_y_xxyy_0, ta2_xx_y_xxyz_0, ta2_xx_y_xxzz_0, ta2_xx_y_xyyy_0, ta2_xx_y_xyyz_0, ta2_xx_y_xyzz_0, ta2_xx_y_xzzz_0, ta2_xx_y_yyyy_0, ta2_xx_y_yyyz_0, ta2_xx_y_yyzz_0, ta2_xx_y_yzzz_0, ta2_xx_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_y_xxxx_0[i] = ta2_xx_0_xxxx_0[i] * pa_y[i] - ta2_xx_0_xxxx_1[i] * pc_y[i];

        ta2_xx_y_xxxy_0[i] = ta2_xx_0_xxx_0[i] * fe_0 - ta2_xx_0_xxx_1[i] * fe_0 + ta2_xx_0_xxxy_0[i] * pa_y[i] - ta2_xx_0_xxxy_1[i] * pc_y[i];

        ta2_xx_y_xxxz_0[i] = ta2_xx_0_xxxz_0[i] * pa_y[i] - ta2_xx_0_xxxz_1[i] * pc_y[i];

        ta2_xx_y_xxyy_0[i] = 2.0 * ta2_xx_0_xxy_0[i] * fe_0 - 2.0 * ta2_xx_0_xxy_1[i] * fe_0 + ta2_xx_0_xxyy_0[i] * pa_y[i] - ta2_xx_0_xxyy_1[i] * pc_y[i];

        ta2_xx_y_xxyz_0[i] = ta2_xx_0_xxz_0[i] * fe_0 - ta2_xx_0_xxz_1[i] * fe_0 + ta2_xx_0_xxyz_0[i] * pa_y[i] - ta2_xx_0_xxyz_1[i] * pc_y[i];

        ta2_xx_y_xxzz_0[i] = ta2_xx_0_xxzz_0[i] * pa_y[i] - ta2_xx_0_xxzz_1[i] * pc_y[i];

        ta2_xx_y_xyyy_0[i] = 3.0 * ta2_xx_0_xyy_0[i] * fe_0 - 3.0 * ta2_xx_0_xyy_1[i] * fe_0 + ta2_xx_0_xyyy_0[i] * pa_y[i] - ta2_xx_0_xyyy_1[i] * pc_y[i];

        ta2_xx_y_xyyz_0[i] = 2.0 * ta2_xx_0_xyz_0[i] * fe_0 - 2.0 * ta2_xx_0_xyz_1[i] * fe_0 + ta2_xx_0_xyyz_0[i] * pa_y[i] - ta2_xx_0_xyyz_1[i] * pc_y[i];

        ta2_xx_y_xyzz_0[i] = ta2_xx_0_xzz_0[i] * fe_0 - ta2_xx_0_xzz_1[i] * fe_0 + ta2_xx_0_xyzz_0[i] * pa_y[i] - ta2_xx_0_xyzz_1[i] * pc_y[i];

        ta2_xx_y_xzzz_0[i] = ta2_xx_0_xzzz_0[i] * pa_y[i] - ta2_xx_0_xzzz_1[i] * pc_y[i];

        ta2_xx_y_yyyy_0[i] = 4.0 * ta2_xx_0_yyy_0[i] * fe_0 - 4.0 * ta2_xx_0_yyy_1[i] * fe_0 + ta2_xx_0_yyyy_0[i] * pa_y[i] - ta2_xx_0_yyyy_1[i] * pc_y[i];

        ta2_xx_y_yyyz_0[i] = 3.0 * ta2_xx_0_yyz_0[i] * fe_0 - 3.0 * ta2_xx_0_yyz_1[i] * fe_0 + ta2_xx_0_yyyz_0[i] * pa_y[i] - ta2_xx_0_yyyz_1[i] * pc_y[i];

        ta2_xx_y_yyzz_0[i] = 2.0 * ta2_xx_0_yzz_0[i] * fe_0 - 2.0 * ta2_xx_0_yzz_1[i] * fe_0 + ta2_xx_0_yyzz_0[i] * pa_y[i] - ta2_xx_0_yyzz_1[i] * pc_y[i];

        ta2_xx_y_yzzz_0[i] = ta2_xx_0_zzz_0[i] * fe_0 - ta2_xx_0_zzz_1[i] * fe_0 + ta2_xx_0_yzzz_0[i] * pa_y[i] - ta2_xx_0_yzzz_1[i] * pc_y[i];

        ta2_xx_y_zzzz_0[i] = ta2_xx_0_zzzz_0[i] * pa_y[i] - ta2_xx_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

    auto ta2_xx_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 30);

    auto ta2_xx_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 31);

    auto ta2_xx_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 32);

    auto ta2_xx_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 33);

    auto ta2_xx_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 34);

    auto ta2_xx_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 35);

    auto ta2_xx_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 36);

    auto ta2_xx_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 37);

    auto ta2_xx_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 38);

    auto ta2_xx_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 39);

    auto ta2_xx_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 40);

    auto ta2_xx_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 41);

    auto ta2_xx_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 42);

    auto ta2_xx_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 43);

    auto ta2_xx_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 44);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_0_xxx_0, ta2_xx_0_xxx_1, ta2_xx_0_xxxx_0, ta2_xx_0_xxxx_1, ta2_xx_0_xxxy_0, ta2_xx_0_xxxy_1, ta2_xx_0_xxxz_0, ta2_xx_0_xxxz_1, ta2_xx_0_xxy_0, ta2_xx_0_xxy_1, ta2_xx_0_xxyy_0, ta2_xx_0_xxyy_1, ta2_xx_0_xxyz_0, ta2_xx_0_xxyz_1, ta2_xx_0_xxz_0, ta2_xx_0_xxz_1, ta2_xx_0_xxzz_0, ta2_xx_0_xxzz_1, ta2_xx_0_xyy_0, ta2_xx_0_xyy_1, ta2_xx_0_xyyy_0, ta2_xx_0_xyyy_1, ta2_xx_0_xyyz_0, ta2_xx_0_xyyz_1, ta2_xx_0_xyz_0, ta2_xx_0_xyz_1, ta2_xx_0_xyzz_0, ta2_xx_0_xyzz_1, ta2_xx_0_xzz_0, ta2_xx_0_xzz_1, ta2_xx_0_xzzz_0, ta2_xx_0_xzzz_1, ta2_xx_0_yyy_0, ta2_xx_0_yyy_1, ta2_xx_0_yyyy_0, ta2_xx_0_yyyy_1, ta2_xx_0_yyyz_0, ta2_xx_0_yyyz_1, ta2_xx_0_yyz_0, ta2_xx_0_yyz_1, ta2_xx_0_yyzz_0, ta2_xx_0_yyzz_1, ta2_xx_0_yzz_0, ta2_xx_0_yzz_1, ta2_xx_0_yzzz_0, ta2_xx_0_yzzz_1, ta2_xx_0_zzz_0, ta2_xx_0_zzz_1, ta2_xx_0_zzzz_0, ta2_xx_0_zzzz_1, ta2_xx_z_xxxx_0, ta2_xx_z_xxxy_0, ta2_xx_z_xxxz_0, ta2_xx_z_xxyy_0, ta2_xx_z_xxyz_0, ta2_xx_z_xxzz_0, ta2_xx_z_xyyy_0, ta2_xx_z_xyyz_0, ta2_xx_z_xyzz_0, ta2_xx_z_xzzz_0, ta2_xx_z_yyyy_0, ta2_xx_z_yyyz_0, ta2_xx_z_yyzz_0, ta2_xx_z_yzzz_0, ta2_xx_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_z_xxxx_0[i] = ta2_xx_0_xxxx_0[i] * pa_z[i] - ta2_xx_0_xxxx_1[i] * pc_z[i];

        ta2_xx_z_xxxy_0[i] = ta2_xx_0_xxxy_0[i] * pa_z[i] - ta2_xx_0_xxxy_1[i] * pc_z[i];

        ta2_xx_z_xxxz_0[i] = ta2_xx_0_xxx_0[i] * fe_0 - ta2_xx_0_xxx_1[i] * fe_0 + ta2_xx_0_xxxz_0[i] * pa_z[i] - ta2_xx_0_xxxz_1[i] * pc_z[i];

        ta2_xx_z_xxyy_0[i] = ta2_xx_0_xxyy_0[i] * pa_z[i] - ta2_xx_0_xxyy_1[i] * pc_z[i];

        ta2_xx_z_xxyz_0[i] = ta2_xx_0_xxy_0[i] * fe_0 - ta2_xx_0_xxy_1[i] * fe_0 + ta2_xx_0_xxyz_0[i] * pa_z[i] - ta2_xx_0_xxyz_1[i] * pc_z[i];

        ta2_xx_z_xxzz_0[i] = 2.0 * ta2_xx_0_xxz_0[i] * fe_0 - 2.0 * ta2_xx_0_xxz_1[i] * fe_0 + ta2_xx_0_xxzz_0[i] * pa_z[i] - ta2_xx_0_xxzz_1[i] * pc_z[i];

        ta2_xx_z_xyyy_0[i] = ta2_xx_0_xyyy_0[i] * pa_z[i] - ta2_xx_0_xyyy_1[i] * pc_z[i];

        ta2_xx_z_xyyz_0[i] = ta2_xx_0_xyy_0[i] * fe_0 - ta2_xx_0_xyy_1[i] * fe_0 + ta2_xx_0_xyyz_0[i] * pa_z[i] - ta2_xx_0_xyyz_1[i] * pc_z[i];

        ta2_xx_z_xyzz_0[i] = 2.0 * ta2_xx_0_xyz_0[i] * fe_0 - 2.0 * ta2_xx_0_xyz_1[i] * fe_0 + ta2_xx_0_xyzz_0[i] * pa_z[i] - ta2_xx_0_xyzz_1[i] * pc_z[i];

        ta2_xx_z_xzzz_0[i] = 3.0 * ta2_xx_0_xzz_0[i] * fe_0 - 3.0 * ta2_xx_0_xzz_1[i] * fe_0 + ta2_xx_0_xzzz_0[i] * pa_z[i] - ta2_xx_0_xzzz_1[i] * pc_z[i];

        ta2_xx_z_yyyy_0[i] = ta2_xx_0_yyyy_0[i] * pa_z[i] - ta2_xx_0_yyyy_1[i] * pc_z[i];

        ta2_xx_z_yyyz_0[i] = ta2_xx_0_yyy_0[i] * fe_0 - ta2_xx_0_yyy_1[i] * fe_0 + ta2_xx_0_yyyz_0[i] * pa_z[i] - ta2_xx_0_yyyz_1[i] * pc_z[i];

        ta2_xx_z_yyzz_0[i] = 2.0 * ta2_xx_0_yyz_0[i] * fe_0 - 2.0 * ta2_xx_0_yyz_1[i] * fe_0 + ta2_xx_0_yyzz_0[i] * pa_z[i] - ta2_xx_0_yyzz_1[i] * pc_z[i];

        ta2_xx_z_yzzz_0[i] = 3.0 * ta2_xx_0_yzz_0[i] * fe_0 - 3.0 * ta2_xx_0_yzz_1[i] * fe_0 + ta2_xx_0_yzzz_0[i] * pa_z[i] - ta2_xx_0_yzzz_1[i] * pc_z[i];

        ta2_xx_z_zzzz_0[i] = 4.0 * ta2_xx_0_zzz_0[i] * fe_0 - 4.0 * ta2_xx_0_zzz_1[i] * fe_0 + ta2_xx_0_zzzz_0[i] * pa_z[i] - ta2_xx_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : PG

    auto ta2_xy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 45);

    auto ta2_xy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 46);

    auto ta2_xy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 47);

    auto ta2_xy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 48);

    auto ta2_xy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 49);

    auto ta2_xy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 50);

    auto ta2_xy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 51);

    auto ta2_xy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 52);

    auto ta2_xy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 53);

    auto ta2_xy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 54);

    auto ta2_xy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 55);

    auto ta2_xy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 56);

    auto ta2_xy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 57);

    auto ta2_xy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 58);

    auto ta2_xy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 59);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_0_xxxx_1, ta1_y_0_xxxy_1, ta1_y_0_xxxz_1, ta1_y_0_xxyy_1, ta1_y_0_xxyz_1, ta1_y_0_xxzz_1, ta1_y_0_xyyy_1, ta1_y_0_xyyz_1, ta1_y_0_xyzz_1, ta1_y_0_xzzz_1, ta1_y_0_yyyy_1, ta1_y_0_yyyz_1, ta1_y_0_yyzz_1, ta1_y_0_yzzz_1, ta1_y_0_zzzz_1, ta2_xy_0_xxx_0, ta2_xy_0_xxx_1, ta2_xy_0_xxxx_0, ta2_xy_0_xxxx_1, ta2_xy_0_xxxy_0, ta2_xy_0_xxxy_1, ta2_xy_0_xxxz_0, ta2_xy_0_xxxz_1, ta2_xy_0_xxy_0, ta2_xy_0_xxy_1, ta2_xy_0_xxyy_0, ta2_xy_0_xxyy_1, ta2_xy_0_xxyz_0, ta2_xy_0_xxyz_1, ta2_xy_0_xxz_0, ta2_xy_0_xxz_1, ta2_xy_0_xxzz_0, ta2_xy_0_xxzz_1, ta2_xy_0_xyy_0, ta2_xy_0_xyy_1, ta2_xy_0_xyyy_0, ta2_xy_0_xyyy_1, ta2_xy_0_xyyz_0, ta2_xy_0_xyyz_1, ta2_xy_0_xyz_0, ta2_xy_0_xyz_1, ta2_xy_0_xyzz_0, ta2_xy_0_xyzz_1, ta2_xy_0_xzz_0, ta2_xy_0_xzz_1, ta2_xy_0_xzzz_0, ta2_xy_0_xzzz_1, ta2_xy_0_yyy_0, ta2_xy_0_yyy_1, ta2_xy_0_yyyy_0, ta2_xy_0_yyyy_1, ta2_xy_0_yyyz_0, ta2_xy_0_yyyz_1, ta2_xy_0_yyz_0, ta2_xy_0_yyz_1, ta2_xy_0_yyzz_0, ta2_xy_0_yyzz_1, ta2_xy_0_yzz_0, ta2_xy_0_yzz_1, ta2_xy_0_yzzz_0, ta2_xy_0_yzzz_1, ta2_xy_0_zzz_0, ta2_xy_0_zzz_1, ta2_xy_0_zzzz_0, ta2_xy_0_zzzz_1, ta2_xy_x_xxxx_0, ta2_xy_x_xxxy_0, ta2_xy_x_xxxz_0, ta2_xy_x_xxyy_0, ta2_xy_x_xxyz_0, ta2_xy_x_xxzz_0, ta2_xy_x_xyyy_0, ta2_xy_x_xyyz_0, ta2_xy_x_xyzz_0, ta2_xy_x_xzzz_0, ta2_xy_x_yyyy_0, ta2_xy_x_yyyz_0, ta2_xy_x_yyzz_0, ta2_xy_x_yzzz_0, ta2_xy_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_x_xxxx_0[i] = 4.0 * ta2_xy_0_xxx_0[i] * fe_0 - 4.0 * ta2_xy_0_xxx_1[i] * fe_0 + ta1_y_0_xxxx_1[i] + ta2_xy_0_xxxx_0[i] * pa_x[i] - ta2_xy_0_xxxx_1[i] * pc_x[i];

        ta2_xy_x_xxxy_0[i] = 3.0 * ta2_xy_0_xxy_0[i] * fe_0 - 3.0 * ta2_xy_0_xxy_1[i] * fe_0 + ta1_y_0_xxxy_1[i] + ta2_xy_0_xxxy_0[i] * pa_x[i] - ta2_xy_0_xxxy_1[i] * pc_x[i];

        ta2_xy_x_xxxz_0[i] = 3.0 * ta2_xy_0_xxz_0[i] * fe_0 - 3.0 * ta2_xy_0_xxz_1[i] * fe_0 + ta1_y_0_xxxz_1[i] + ta2_xy_0_xxxz_0[i] * pa_x[i] - ta2_xy_0_xxxz_1[i] * pc_x[i];

        ta2_xy_x_xxyy_0[i] = 2.0 * ta2_xy_0_xyy_0[i] * fe_0 - 2.0 * ta2_xy_0_xyy_1[i] * fe_0 + ta1_y_0_xxyy_1[i] + ta2_xy_0_xxyy_0[i] * pa_x[i] - ta2_xy_0_xxyy_1[i] * pc_x[i];

        ta2_xy_x_xxyz_0[i] = 2.0 * ta2_xy_0_xyz_0[i] * fe_0 - 2.0 * ta2_xy_0_xyz_1[i] * fe_0 + ta1_y_0_xxyz_1[i] + ta2_xy_0_xxyz_0[i] * pa_x[i] - ta2_xy_0_xxyz_1[i] * pc_x[i];

        ta2_xy_x_xxzz_0[i] = 2.0 * ta2_xy_0_xzz_0[i] * fe_0 - 2.0 * ta2_xy_0_xzz_1[i] * fe_0 + ta1_y_0_xxzz_1[i] + ta2_xy_0_xxzz_0[i] * pa_x[i] - ta2_xy_0_xxzz_1[i] * pc_x[i];

        ta2_xy_x_xyyy_0[i] = ta2_xy_0_yyy_0[i] * fe_0 - ta2_xy_0_yyy_1[i] * fe_0 + ta1_y_0_xyyy_1[i] + ta2_xy_0_xyyy_0[i] * pa_x[i] - ta2_xy_0_xyyy_1[i] * pc_x[i];

        ta2_xy_x_xyyz_0[i] = ta2_xy_0_yyz_0[i] * fe_0 - ta2_xy_0_yyz_1[i] * fe_0 + ta1_y_0_xyyz_1[i] + ta2_xy_0_xyyz_0[i] * pa_x[i] - ta2_xy_0_xyyz_1[i] * pc_x[i];

        ta2_xy_x_xyzz_0[i] = ta2_xy_0_yzz_0[i] * fe_0 - ta2_xy_0_yzz_1[i] * fe_0 + ta1_y_0_xyzz_1[i] + ta2_xy_0_xyzz_0[i] * pa_x[i] - ta2_xy_0_xyzz_1[i] * pc_x[i];

        ta2_xy_x_xzzz_0[i] = ta2_xy_0_zzz_0[i] * fe_0 - ta2_xy_0_zzz_1[i] * fe_0 + ta1_y_0_xzzz_1[i] + ta2_xy_0_xzzz_0[i] * pa_x[i] - ta2_xy_0_xzzz_1[i] * pc_x[i];

        ta2_xy_x_yyyy_0[i] = ta1_y_0_yyyy_1[i] + ta2_xy_0_yyyy_0[i] * pa_x[i] - ta2_xy_0_yyyy_1[i] * pc_x[i];

        ta2_xy_x_yyyz_0[i] = ta1_y_0_yyyz_1[i] + ta2_xy_0_yyyz_0[i] * pa_x[i] - ta2_xy_0_yyyz_1[i] * pc_x[i];

        ta2_xy_x_yyzz_0[i] = ta1_y_0_yyzz_1[i] + ta2_xy_0_yyzz_0[i] * pa_x[i] - ta2_xy_0_yyzz_1[i] * pc_x[i];

        ta2_xy_x_yzzz_0[i] = ta1_y_0_yzzz_1[i] + ta2_xy_0_yzzz_0[i] * pa_x[i] - ta2_xy_0_yzzz_1[i] * pc_x[i];

        ta2_xy_x_zzzz_0[i] = ta1_y_0_zzzz_1[i] + ta2_xy_0_zzzz_0[i] * pa_x[i] - ta2_xy_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : PG

    auto ta2_xy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 60);

    auto ta2_xy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 61);

    auto ta2_xy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 62);

    auto ta2_xy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 63);

    auto ta2_xy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 64);

    auto ta2_xy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 65);

    auto ta2_xy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 66);

    auto ta2_xy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 67);

    auto ta2_xy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 68);

    auto ta2_xy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 69);

    auto ta2_xy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 70);

    auto ta2_xy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 71);

    auto ta2_xy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 72);

    auto ta2_xy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 73);

    auto ta2_xy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 74);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_0_xxxx_1, ta1_x_0_xxxy_1, ta1_x_0_xxxz_1, ta1_x_0_xxyy_1, ta1_x_0_xxyz_1, ta1_x_0_xxzz_1, ta1_x_0_xyyy_1, ta1_x_0_xyyz_1, ta1_x_0_xyzz_1, ta1_x_0_xzzz_1, ta1_x_0_yyyy_1, ta1_x_0_yyyz_1, ta1_x_0_yyzz_1, ta1_x_0_yzzz_1, ta1_x_0_zzzz_1, ta2_xy_0_xxx_0, ta2_xy_0_xxx_1, ta2_xy_0_xxxx_0, ta2_xy_0_xxxx_1, ta2_xy_0_xxxy_0, ta2_xy_0_xxxy_1, ta2_xy_0_xxxz_0, ta2_xy_0_xxxz_1, ta2_xy_0_xxy_0, ta2_xy_0_xxy_1, ta2_xy_0_xxyy_0, ta2_xy_0_xxyy_1, ta2_xy_0_xxyz_0, ta2_xy_0_xxyz_1, ta2_xy_0_xxz_0, ta2_xy_0_xxz_1, ta2_xy_0_xxzz_0, ta2_xy_0_xxzz_1, ta2_xy_0_xyy_0, ta2_xy_0_xyy_1, ta2_xy_0_xyyy_0, ta2_xy_0_xyyy_1, ta2_xy_0_xyyz_0, ta2_xy_0_xyyz_1, ta2_xy_0_xyz_0, ta2_xy_0_xyz_1, ta2_xy_0_xyzz_0, ta2_xy_0_xyzz_1, ta2_xy_0_xzz_0, ta2_xy_0_xzz_1, ta2_xy_0_xzzz_0, ta2_xy_0_xzzz_1, ta2_xy_0_yyy_0, ta2_xy_0_yyy_1, ta2_xy_0_yyyy_0, ta2_xy_0_yyyy_1, ta2_xy_0_yyyz_0, ta2_xy_0_yyyz_1, ta2_xy_0_yyz_0, ta2_xy_0_yyz_1, ta2_xy_0_yyzz_0, ta2_xy_0_yyzz_1, ta2_xy_0_yzz_0, ta2_xy_0_yzz_1, ta2_xy_0_yzzz_0, ta2_xy_0_yzzz_1, ta2_xy_0_zzz_0, ta2_xy_0_zzz_1, ta2_xy_0_zzzz_0, ta2_xy_0_zzzz_1, ta2_xy_y_xxxx_0, ta2_xy_y_xxxy_0, ta2_xy_y_xxxz_0, ta2_xy_y_xxyy_0, ta2_xy_y_xxyz_0, ta2_xy_y_xxzz_0, ta2_xy_y_xyyy_0, ta2_xy_y_xyyz_0, ta2_xy_y_xyzz_0, ta2_xy_y_xzzz_0, ta2_xy_y_yyyy_0, ta2_xy_y_yyyz_0, ta2_xy_y_yyzz_0, ta2_xy_y_yzzz_0, ta2_xy_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_y_xxxx_0[i] = ta1_x_0_xxxx_1[i] + ta2_xy_0_xxxx_0[i] * pa_y[i] - ta2_xy_0_xxxx_1[i] * pc_y[i];

        ta2_xy_y_xxxy_0[i] = ta2_xy_0_xxx_0[i] * fe_0 - ta2_xy_0_xxx_1[i] * fe_0 + ta1_x_0_xxxy_1[i] + ta2_xy_0_xxxy_0[i] * pa_y[i] - ta2_xy_0_xxxy_1[i] * pc_y[i];

        ta2_xy_y_xxxz_0[i] = ta1_x_0_xxxz_1[i] + ta2_xy_0_xxxz_0[i] * pa_y[i] - ta2_xy_0_xxxz_1[i] * pc_y[i];

        ta2_xy_y_xxyy_0[i] = 2.0 * ta2_xy_0_xxy_0[i] * fe_0 - 2.0 * ta2_xy_0_xxy_1[i] * fe_0 + ta1_x_0_xxyy_1[i] + ta2_xy_0_xxyy_0[i] * pa_y[i] - ta2_xy_0_xxyy_1[i] * pc_y[i];

        ta2_xy_y_xxyz_0[i] = ta2_xy_0_xxz_0[i] * fe_0 - ta2_xy_0_xxz_1[i] * fe_0 + ta1_x_0_xxyz_1[i] + ta2_xy_0_xxyz_0[i] * pa_y[i] - ta2_xy_0_xxyz_1[i] * pc_y[i];

        ta2_xy_y_xxzz_0[i] = ta1_x_0_xxzz_1[i] + ta2_xy_0_xxzz_0[i] * pa_y[i] - ta2_xy_0_xxzz_1[i] * pc_y[i];

        ta2_xy_y_xyyy_0[i] = 3.0 * ta2_xy_0_xyy_0[i] * fe_0 - 3.0 * ta2_xy_0_xyy_1[i] * fe_0 + ta1_x_0_xyyy_1[i] + ta2_xy_0_xyyy_0[i] * pa_y[i] - ta2_xy_0_xyyy_1[i] * pc_y[i];

        ta2_xy_y_xyyz_0[i] = 2.0 * ta2_xy_0_xyz_0[i] * fe_0 - 2.0 * ta2_xy_0_xyz_1[i] * fe_0 + ta1_x_0_xyyz_1[i] + ta2_xy_0_xyyz_0[i] * pa_y[i] - ta2_xy_0_xyyz_1[i] * pc_y[i];

        ta2_xy_y_xyzz_0[i] = ta2_xy_0_xzz_0[i] * fe_0 - ta2_xy_0_xzz_1[i] * fe_0 + ta1_x_0_xyzz_1[i] + ta2_xy_0_xyzz_0[i] * pa_y[i] - ta2_xy_0_xyzz_1[i] * pc_y[i];

        ta2_xy_y_xzzz_0[i] = ta1_x_0_xzzz_1[i] + ta2_xy_0_xzzz_0[i] * pa_y[i] - ta2_xy_0_xzzz_1[i] * pc_y[i];

        ta2_xy_y_yyyy_0[i] = 4.0 * ta2_xy_0_yyy_0[i] * fe_0 - 4.0 * ta2_xy_0_yyy_1[i] * fe_0 + ta1_x_0_yyyy_1[i] + ta2_xy_0_yyyy_0[i] * pa_y[i] - ta2_xy_0_yyyy_1[i] * pc_y[i];

        ta2_xy_y_yyyz_0[i] = 3.0 * ta2_xy_0_yyz_0[i] * fe_0 - 3.0 * ta2_xy_0_yyz_1[i] * fe_0 + ta1_x_0_yyyz_1[i] + ta2_xy_0_yyyz_0[i] * pa_y[i] - ta2_xy_0_yyyz_1[i] * pc_y[i];

        ta2_xy_y_yyzz_0[i] = 2.0 * ta2_xy_0_yzz_0[i] * fe_0 - 2.0 * ta2_xy_0_yzz_1[i] * fe_0 + ta1_x_0_yyzz_1[i] + ta2_xy_0_yyzz_0[i] * pa_y[i] - ta2_xy_0_yyzz_1[i] * pc_y[i];

        ta2_xy_y_yzzz_0[i] = ta2_xy_0_zzz_0[i] * fe_0 - ta2_xy_0_zzz_1[i] * fe_0 + ta1_x_0_yzzz_1[i] + ta2_xy_0_yzzz_0[i] * pa_y[i] - ta2_xy_0_yzzz_1[i] * pc_y[i];

        ta2_xy_y_zzzz_0[i] = ta1_x_0_zzzz_1[i] + ta2_xy_0_zzzz_0[i] * pa_y[i] - ta2_xy_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : PG

    auto ta2_xy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 75);

    auto ta2_xy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 76);

    auto ta2_xy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 77);

    auto ta2_xy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 78);

    auto ta2_xy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 79);

    auto ta2_xy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 80);

    auto ta2_xy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 81);

    auto ta2_xy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 82);

    auto ta2_xy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 83);

    auto ta2_xy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 84);

    auto ta2_xy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 85);

    auto ta2_xy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 86);

    auto ta2_xy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 87);

    auto ta2_xy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 88);

    auto ta2_xy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 89);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_0_xxx_0, ta2_xy_0_xxx_1, ta2_xy_0_xxxx_0, ta2_xy_0_xxxx_1, ta2_xy_0_xxxy_0, ta2_xy_0_xxxy_1, ta2_xy_0_xxxz_0, ta2_xy_0_xxxz_1, ta2_xy_0_xxy_0, ta2_xy_0_xxy_1, ta2_xy_0_xxyy_0, ta2_xy_0_xxyy_1, ta2_xy_0_xxyz_0, ta2_xy_0_xxyz_1, ta2_xy_0_xxz_0, ta2_xy_0_xxz_1, ta2_xy_0_xxzz_0, ta2_xy_0_xxzz_1, ta2_xy_0_xyy_0, ta2_xy_0_xyy_1, ta2_xy_0_xyyy_0, ta2_xy_0_xyyy_1, ta2_xy_0_xyyz_0, ta2_xy_0_xyyz_1, ta2_xy_0_xyz_0, ta2_xy_0_xyz_1, ta2_xy_0_xyzz_0, ta2_xy_0_xyzz_1, ta2_xy_0_xzz_0, ta2_xy_0_xzz_1, ta2_xy_0_xzzz_0, ta2_xy_0_xzzz_1, ta2_xy_0_yyy_0, ta2_xy_0_yyy_1, ta2_xy_0_yyyy_0, ta2_xy_0_yyyy_1, ta2_xy_0_yyyz_0, ta2_xy_0_yyyz_1, ta2_xy_0_yyz_0, ta2_xy_0_yyz_1, ta2_xy_0_yyzz_0, ta2_xy_0_yyzz_1, ta2_xy_0_yzz_0, ta2_xy_0_yzz_1, ta2_xy_0_yzzz_0, ta2_xy_0_yzzz_1, ta2_xy_0_zzz_0, ta2_xy_0_zzz_1, ta2_xy_0_zzzz_0, ta2_xy_0_zzzz_1, ta2_xy_z_xxxx_0, ta2_xy_z_xxxy_0, ta2_xy_z_xxxz_0, ta2_xy_z_xxyy_0, ta2_xy_z_xxyz_0, ta2_xy_z_xxzz_0, ta2_xy_z_xyyy_0, ta2_xy_z_xyyz_0, ta2_xy_z_xyzz_0, ta2_xy_z_xzzz_0, ta2_xy_z_yyyy_0, ta2_xy_z_yyyz_0, ta2_xy_z_yyzz_0, ta2_xy_z_yzzz_0, ta2_xy_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_z_xxxx_0[i] = ta2_xy_0_xxxx_0[i] * pa_z[i] - ta2_xy_0_xxxx_1[i] * pc_z[i];

        ta2_xy_z_xxxy_0[i] = ta2_xy_0_xxxy_0[i] * pa_z[i] - ta2_xy_0_xxxy_1[i] * pc_z[i];

        ta2_xy_z_xxxz_0[i] = ta2_xy_0_xxx_0[i] * fe_0 - ta2_xy_0_xxx_1[i] * fe_0 + ta2_xy_0_xxxz_0[i] * pa_z[i] - ta2_xy_0_xxxz_1[i] * pc_z[i];

        ta2_xy_z_xxyy_0[i] = ta2_xy_0_xxyy_0[i] * pa_z[i] - ta2_xy_0_xxyy_1[i] * pc_z[i];

        ta2_xy_z_xxyz_0[i] = ta2_xy_0_xxy_0[i] * fe_0 - ta2_xy_0_xxy_1[i] * fe_0 + ta2_xy_0_xxyz_0[i] * pa_z[i] - ta2_xy_0_xxyz_1[i] * pc_z[i];

        ta2_xy_z_xxzz_0[i] = 2.0 * ta2_xy_0_xxz_0[i] * fe_0 - 2.0 * ta2_xy_0_xxz_1[i] * fe_0 + ta2_xy_0_xxzz_0[i] * pa_z[i] - ta2_xy_0_xxzz_1[i] * pc_z[i];

        ta2_xy_z_xyyy_0[i] = ta2_xy_0_xyyy_0[i] * pa_z[i] - ta2_xy_0_xyyy_1[i] * pc_z[i];

        ta2_xy_z_xyyz_0[i] = ta2_xy_0_xyy_0[i] * fe_0 - ta2_xy_0_xyy_1[i] * fe_0 + ta2_xy_0_xyyz_0[i] * pa_z[i] - ta2_xy_0_xyyz_1[i] * pc_z[i];

        ta2_xy_z_xyzz_0[i] = 2.0 * ta2_xy_0_xyz_0[i] * fe_0 - 2.0 * ta2_xy_0_xyz_1[i] * fe_0 + ta2_xy_0_xyzz_0[i] * pa_z[i] - ta2_xy_0_xyzz_1[i] * pc_z[i];

        ta2_xy_z_xzzz_0[i] = 3.0 * ta2_xy_0_xzz_0[i] * fe_0 - 3.0 * ta2_xy_0_xzz_1[i] * fe_0 + ta2_xy_0_xzzz_0[i] * pa_z[i] - ta2_xy_0_xzzz_1[i] * pc_z[i];

        ta2_xy_z_yyyy_0[i] = ta2_xy_0_yyyy_0[i] * pa_z[i] - ta2_xy_0_yyyy_1[i] * pc_z[i];

        ta2_xy_z_yyyz_0[i] = ta2_xy_0_yyy_0[i] * fe_0 - ta2_xy_0_yyy_1[i] * fe_0 + ta2_xy_0_yyyz_0[i] * pa_z[i] - ta2_xy_0_yyyz_1[i] * pc_z[i];

        ta2_xy_z_yyzz_0[i] = 2.0 * ta2_xy_0_yyz_0[i] * fe_0 - 2.0 * ta2_xy_0_yyz_1[i] * fe_0 + ta2_xy_0_yyzz_0[i] * pa_z[i] - ta2_xy_0_yyzz_1[i] * pc_z[i];

        ta2_xy_z_yzzz_0[i] = 3.0 * ta2_xy_0_yzz_0[i] * fe_0 - 3.0 * ta2_xy_0_yzz_1[i] * fe_0 + ta2_xy_0_yzzz_0[i] * pa_z[i] - ta2_xy_0_yzzz_1[i] * pc_z[i];

        ta2_xy_z_zzzz_0[i] = 4.0 * ta2_xy_0_zzz_0[i] * fe_0 - 4.0 * ta2_xy_0_zzz_1[i] * fe_0 + ta2_xy_0_zzzz_0[i] * pa_z[i] - ta2_xy_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 90-105 components of targeted buffer : PG

    auto ta2_xz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 90);

    auto ta2_xz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 91);

    auto ta2_xz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 92);

    auto ta2_xz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 93);

    auto ta2_xz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 94);

    auto ta2_xz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 95);

    auto ta2_xz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 96);

    auto ta2_xz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 97);

    auto ta2_xz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 98);

    auto ta2_xz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 99);

    auto ta2_xz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 100);

    auto ta2_xz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 101);

    auto ta2_xz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 102);

    auto ta2_xz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 103);

    auto ta2_xz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 104);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_0_xxxx_1, ta1_z_0_xxxy_1, ta1_z_0_xxxz_1, ta1_z_0_xxyy_1, ta1_z_0_xxyz_1, ta1_z_0_xxzz_1, ta1_z_0_xyyy_1, ta1_z_0_xyyz_1, ta1_z_0_xyzz_1, ta1_z_0_xzzz_1, ta1_z_0_yyyy_1, ta1_z_0_yyyz_1, ta1_z_0_yyzz_1, ta1_z_0_yzzz_1, ta1_z_0_zzzz_1, ta2_xz_0_xxx_0, ta2_xz_0_xxx_1, ta2_xz_0_xxxx_0, ta2_xz_0_xxxx_1, ta2_xz_0_xxxy_0, ta2_xz_0_xxxy_1, ta2_xz_0_xxxz_0, ta2_xz_0_xxxz_1, ta2_xz_0_xxy_0, ta2_xz_0_xxy_1, ta2_xz_0_xxyy_0, ta2_xz_0_xxyy_1, ta2_xz_0_xxyz_0, ta2_xz_0_xxyz_1, ta2_xz_0_xxz_0, ta2_xz_0_xxz_1, ta2_xz_0_xxzz_0, ta2_xz_0_xxzz_1, ta2_xz_0_xyy_0, ta2_xz_0_xyy_1, ta2_xz_0_xyyy_0, ta2_xz_0_xyyy_1, ta2_xz_0_xyyz_0, ta2_xz_0_xyyz_1, ta2_xz_0_xyz_0, ta2_xz_0_xyz_1, ta2_xz_0_xyzz_0, ta2_xz_0_xyzz_1, ta2_xz_0_xzz_0, ta2_xz_0_xzz_1, ta2_xz_0_xzzz_0, ta2_xz_0_xzzz_1, ta2_xz_0_yyy_0, ta2_xz_0_yyy_1, ta2_xz_0_yyyy_0, ta2_xz_0_yyyy_1, ta2_xz_0_yyyz_0, ta2_xz_0_yyyz_1, ta2_xz_0_yyz_0, ta2_xz_0_yyz_1, ta2_xz_0_yyzz_0, ta2_xz_0_yyzz_1, ta2_xz_0_yzz_0, ta2_xz_0_yzz_1, ta2_xz_0_yzzz_0, ta2_xz_0_yzzz_1, ta2_xz_0_zzz_0, ta2_xz_0_zzz_1, ta2_xz_0_zzzz_0, ta2_xz_0_zzzz_1, ta2_xz_x_xxxx_0, ta2_xz_x_xxxy_0, ta2_xz_x_xxxz_0, ta2_xz_x_xxyy_0, ta2_xz_x_xxyz_0, ta2_xz_x_xxzz_0, ta2_xz_x_xyyy_0, ta2_xz_x_xyyz_0, ta2_xz_x_xyzz_0, ta2_xz_x_xzzz_0, ta2_xz_x_yyyy_0, ta2_xz_x_yyyz_0, ta2_xz_x_yyzz_0, ta2_xz_x_yzzz_0, ta2_xz_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_x_xxxx_0[i] = 4.0 * ta2_xz_0_xxx_0[i] * fe_0 - 4.0 * ta2_xz_0_xxx_1[i] * fe_0 + ta1_z_0_xxxx_1[i] + ta2_xz_0_xxxx_0[i] * pa_x[i] - ta2_xz_0_xxxx_1[i] * pc_x[i];

        ta2_xz_x_xxxy_0[i] = 3.0 * ta2_xz_0_xxy_0[i] * fe_0 - 3.0 * ta2_xz_0_xxy_1[i] * fe_0 + ta1_z_0_xxxy_1[i] + ta2_xz_0_xxxy_0[i] * pa_x[i] - ta2_xz_0_xxxy_1[i] * pc_x[i];

        ta2_xz_x_xxxz_0[i] = 3.0 * ta2_xz_0_xxz_0[i] * fe_0 - 3.0 * ta2_xz_0_xxz_1[i] * fe_0 + ta1_z_0_xxxz_1[i] + ta2_xz_0_xxxz_0[i] * pa_x[i] - ta2_xz_0_xxxz_1[i] * pc_x[i];

        ta2_xz_x_xxyy_0[i] = 2.0 * ta2_xz_0_xyy_0[i] * fe_0 - 2.0 * ta2_xz_0_xyy_1[i] * fe_0 + ta1_z_0_xxyy_1[i] + ta2_xz_0_xxyy_0[i] * pa_x[i] - ta2_xz_0_xxyy_1[i] * pc_x[i];

        ta2_xz_x_xxyz_0[i] = 2.0 * ta2_xz_0_xyz_0[i] * fe_0 - 2.0 * ta2_xz_0_xyz_1[i] * fe_0 + ta1_z_0_xxyz_1[i] + ta2_xz_0_xxyz_0[i] * pa_x[i] - ta2_xz_0_xxyz_1[i] * pc_x[i];

        ta2_xz_x_xxzz_0[i] = 2.0 * ta2_xz_0_xzz_0[i] * fe_0 - 2.0 * ta2_xz_0_xzz_1[i] * fe_0 + ta1_z_0_xxzz_1[i] + ta2_xz_0_xxzz_0[i] * pa_x[i] - ta2_xz_0_xxzz_1[i] * pc_x[i];

        ta2_xz_x_xyyy_0[i] = ta2_xz_0_yyy_0[i] * fe_0 - ta2_xz_0_yyy_1[i] * fe_0 + ta1_z_0_xyyy_1[i] + ta2_xz_0_xyyy_0[i] * pa_x[i] - ta2_xz_0_xyyy_1[i] * pc_x[i];

        ta2_xz_x_xyyz_0[i] = ta2_xz_0_yyz_0[i] * fe_0 - ta2_xz_0_yyz_1[i] * fe_0 + ta1_z_0_xyyz_1[i] + ta2_xz_0_xyyz_0[i] * pa_x[i] - ta2_xz_0_xyyz_1[i] * pc_x[i];

        ta2_xz_x_xyzz_0[i] = ta2_xz_0_yzz_0[i] * fe_0 - ta2_xz_0_yzz_1[i] * fe_0 + ta1_z_0_xyzz_1[i] + ta2_xz_0_xyzz_0[i] * pa_x[i] - ta2_xz_0_xyzz_1[i] * pc_x[i];

        ta2_xz_x_xzzz_0[i] = ta2_xz_0_zzz_0[i] * fe_0 - ta2_xz_0_zzz_1[i] * fe_0 + ta1_z_0_xzzz_1[i] + ta2_xz_0_xzzz_0[i] * pa_x[i] - ta2_xz_0_xzzz_1[i] * pc_x[i];

        ta2_xz_x_yyyy_0[i] = ta1_z_0_yyyy_1[i] + ta2_xz_0_yyyy_0[i] * pa_x[i] - ta2_xz_0_yyyy_1[i] * pc_x[i];

        ta2_xz_x_yyyz_0[i] = ta1_z_0_yyyz_1[i] + ta2_xz_0_yyyz_0[i] * pa_x[i] - ta2_xz_0_yyyz_1[i] * pc_x[i];

        ta2_xz_x_yyzz_0[i] = ta1_z_0_yyzz_1[i] + ta2_xz_0_yyzz_0[i] * pa_x[i] - ta2_xz_0_yyzz_1[i] * pc_x[i];

        ta2_xz_x_yzzz_0[i] = ta1_z_0_yzzz_1[i] + ta2_xz_0_yzzz_0[i] * pa_x[i] - ta2_xz_0_yzzz_1[i] * pc_x[i];

        ta2_xz_x_zzzz_0[i] = ta1_z_0_zzzz_1[i] + ta2_xz_0_zzzz_0[i] * pa_x[i] - ta2_xz_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : PG

    auto ta2_xz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 105);

    auto ta2_xz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 106);

    auto ta2_xz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 107);

    auto ta2_xz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 108);

    auto ta2_xz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 109);

    auto ta2_xz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 110);

    auto ta2_xz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 111);

    auto ta2_xz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 112);

    auto ta2_xz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 113);

    auto ta2_xz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 114);

    auto ta2_xz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 115);

    auto ta2_xz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 116);

    auto ta2_xz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 117);

    auto ta2_xz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 118);

    auto ta2_xz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 119);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_0_xxx_0, ta2_xz_0_xxx_1, ta2_xz_0_xxxx_0, ta2_xz_0_xxxx_1, ta2_xz_0_xxxy_0, ta2_xz_0_xxxy_1, ta2_xz_0_xxxz_0, ta2_xz_0_xxxz_1, ta2_xz_0_xxy_0, ta2_xz_0_xxy_1, ta2_xz_0_xxyy_0, ta2_xz_0_xxyy_1, ta2_xz_0_xxyz_0, ta2_xz_0_xxyz_1, ta2_xz_0_xxz_0, ta2_xz_0_xxz_1, ta2_xz_0_xxzz_0, ta2_xz_0_xxzz_1, ta2_xz_0_xyy_0, ta2_xz_0_xyy_1, ta2_xz_0_xyyy_0, ta2_xz_0_xyyy_1, ta2_xz_0_xyyz_0, ta2_xz_0_xyyz_1, ta2_xz_0_xyz_0, ta2_xz_0_xyz_1, ta2_xz_0_xyzz_0, ta2_xz_0_xyzz_1, ta2_xz_0_xzz_0, ta2_xz_0_xzz_1, ta2_xz_0_xzzz_0, ta2_xz_0_xzzz_1, ta2_xz_0_yyy_0, ta2_xz_0_yyy_1, ta2_xz_0_yyyy_0, ta2_xz_0_yyyy_1, ta2_xz_0_yyyz_0, ta2_xz_0_yyyz_1, ta2_xz_0_yyz_0, ta2_xz_0_yyz_1, ta2_xz_0_yyzz_0, ta2_xz_0_yyzz_1, ta2_xz_0_yzz_0, ta2_xz_0_yzz_1, ta2_xz_0_yzzz_0, ta2_xz_0_yzzz_1, ta2_xz_0_zzz_0, ta2_xz_0_zzz_1, ta2_xz_0_zzzz_0, ta2_xz_0_zzzz_1, ta2_xz_y_xxxx_0, ta2_xz_y_xxxy_0, ta2_xz_y_xxxz_0, ta2_xz_y_xxyy_0, ta2_xz_y_xxyz_0, ta2_xz_y_xxzz_0, ta2_xz_y_xyyy_0, ta2_xz_y_xyyz_0, ta2_xz_y_xyzz_0, ta2_xz_y_xzzz_0, ta2_xz_y_yyyy_0, ta2_xz_y_yyyz_0, ta2_xz_y_yyzz_0, ta2_xz_y_yzzz_0, ta2_xz_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_y_xxxx_0[i] = ta2_xz_0_xxxx_0[i] * pa_y[i] - ta2_xz_0_xxxx_1[i] * pc_y[i];

        ta2_xz_y_xxxy_0[i] = ta2_xz_0_xxx_0[i] * fe_0 - ta2_xz_0_xxx_1[i] * fe_0 + ta2_xz_0_xxxy_0[i] * pa_y[i] - ta2_xz_0_xxxy_1[i] * pc_y[i];

        ta2_xz_y_xxxz_0[i] = ta2_xz_0_xxxz_0[i] * pa_y[i] - ta2_xz_0_xxxz_1[i] * pc_y[i];

        ta2_xz_y_xxyy_0[i] = 2.0 * ta2_xz_0_xxy_0[i] * fe_0 - 2.0 * ta2_xz_0_xxy_1[i] * fe_0 + ta2_xz_0_xxyy_0[i] * pa_y[i] - ta2_xz_0_xxyy_1[i] * pc_y[i];

        ta2_xz_y_xxyz_0[i] = ta2_xz_0_xxz_0[i] * fe_0 - ta2_xz_0_xxz_1[i] * fe_0 + ta2_xz_0_xxyz_0[i] * pa_y[i] - ta2_xz_0_xxyz_1[i] * pc_y[i];

        ta2_xz_y_xxzz_0[i] = ta2_xz_0_xxzz_0[i] * pa_y[i] - ta2_xz_0_xxzz_1[i] * pc_y[i];

        ta2_xz_y_xyyy_0[i] = 3.0 * ta2_xz_0_xyy_0[i] * fe_0 - 3.0 * ta2_xz_0_xyy_1[i] * fe_0 + ta2_xz_0_xyyy_0[i] * pa_y[i] - ta2_xz_0_xyyy_1[i] * pc_y[i];

        ta2_xz_y_xyyz_0[i] = 2.0 * ta2_xz_0_xyz_0[i] * fe_0 - 2.0 * ta2_xz_0_xyz_1[i] * fe_0 + ta2_xz_0_xyyz_0[i] * pa_y[i] - ta2_xz_0_xyyz_1[i] * pc_y[i];

        ta2_xz_y_xyzz_0[i] = ta2_xz_0_xzz_0[i] * fe_0 - ta2_xz_0_xzz_1[i] * fe_0 + ta2_xz_0_xyzz_0[i] * pa_y[i] - ta2_xz_0_xyzz_1[i] * pc_y[i];

        ta2_xz_y_xzzz_0[i] = ta2_xz_0_xzzz_0[i] * pa_y[i] - ta2_xz_0_xzzz_1[i] * pc_y[i];

        ta2_xz_y_yyyy_0[i] = 4.0 * ta2_xz_0_yyy_0[i] * fe_0 - 4.0 * ta2_xz_0_yyy_1[i] * fe_0 + ta2_xz_0_yyyy_0[i] * pa_y[i] - ta2_xz_0_yyyy_1[i] * pc_y[i];

        ta2_xz_y_yyyz_0[i] = 3.0 * ta2_xz_0_yyz_0[i] * fe_0 - 3.0 * ta2_xz_0_yyz_1[i] * fe_0 + ta2_xz_0_yyyz_0[i] * pa_y[i] - ta2_xz_0_yyyz_1[i] * pc_y[i];

        ta2_xz_y_yyzz_0[i] = 2.0 * ta2_xz_0_yzz_0[i] * fe_0 - 2.0 * ta2_xz_0_yzz_1[i] * fe_0 + ta2_xz_0_yyzz_0[i] * pa_y[i] - ta2_xz_0_yyzz_1[i] * pc_y[i];

        ta2_xz_y_yzzz_0[i] = ta2_xz_0_zzz_0[i] * fe_0 - ta2_xz_0_zzz_1[i] * fe_0 + ta2_xz_0_yzzz_0[i] * pa_y[i] - ta2_xz_0_yzzz_1[i] * pc_y[i];

        ta2_xz_y_zzzz_0[i] = ta2_xz_0_zzzz_0[i] * pa_y[i] - ta2_xz_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : PG

    auto ta2_xz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 120);

    auto ta2_xz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 121);

    auto ta2_xz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 122);

    auto ta2_xz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 123);

    auto ta2_xz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 124);

    auto ta2_xz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 125);

    auto ta2_xz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 126);

    auto ta2_xz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 127);

    auto ta2_xz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 128);

    auto ta2_xz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 129);

    auto ta2_xz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 130);

    auto ta2_xz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 131);

    auto ta2_xz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 132);

    auto ta2_xz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 133);

    auto ta2_xz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 134);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_0_xxxx_1, ta1_x_0_xxxy_1, ta1_x_0_xxxz_1, ta1_x_0_xxyy_1, ta1_x_0_xxyz_1, ta1_x_0_xxzz_1, ta1_x_0_xyyy_1, ta1_x_0_xyyz_1, ta1_x_0_xyzz_1, ta1_x_0_xzzz_1, ta1_x_0_yyyy_1, ta1_x_0_yyyz_1, ta1_x_0_yyzz_1, ta1_x_0_yzzz_1, ta1_x_0_zzzz_1, ta2_xz_0_xxx_0, ta2_xz_0_xxx_1, ta2_xz_0_xxxx_0, ta2_xz_0_xxxx_1, ta2_xz_0_xxxy_0, ta2_xz_0_xxxy_1, ta2_xz_0_xxxz_0, ta2_xz_0_xxxz_1, ta2_xz_0_xxy_0, ta2_xz_0_xxy_1, ta2_xz_0_xxyy_0, ta2_xz_0_xxyy_1, ta2_xz_0_xxyz_0, ta2_xz_0_xxyz_1, ta2_xz_0_xxz_0, ta2_xz_0_xxz_1, ta2_xz_0_xxzz_0, ta2_xz_0_xxzz_1, ta2_xz_0_xyy_0, ta2_xz_0_xyy_1, ta2_xz_0_xyyy_0, ta2_xz_0_xyyy_1, ta2_xz_0_xyyz_0, ta2_xz_0_xyyz_1, ta2_xz_0_xyz_0, ta2_xz_0_xyz_1, ta2_xz_0_xyzz_0, ta2_xz_0_xyzz_1, ta2_xz_0_xzz_0, ta2_xz_0_xzz_1, ta2_xz_0_xzzz_0, ta2_xz_0_xzzz_1, ta2_xz_0_yyy_0, ta2_xz_0_yyy_1, ta2_xz_0_yyyy_0, ta2_xz_0_yyyy_1, ta2_xz_0_yyyz_0, ta2_xz_0_yyyz_1, ta2_xz_0_yyz_0, ta2_xz_0_yyz_1, ta2_xz_0_yyzz_0, ta2_xz_0_yyzz_1, ta2_xz_0_yzz_0, ta2_xz_0_yzz_1, ta2_xz_0_yzzz_0, ta2_xz_0_yzzz_1, ta2_xz_0_zzz_0, ta2_xz_0_zzz_1, ta2_xz_0_zzzz_0, ta2_xz_0_zzzz_1, ta2_xz_z_xxxx_0, ta2_xz_z_xxxy_0, ta2_xz_z_xxxz_0, ta2_xz_z_xxyy_0, ta2_xz_z_xxyz_0, ta2_xz_z_xxzz_0, ta2_xz_z_xyyy_0, ta2_xz_z_xyyz_0, ta2_xz_z_xyzz_0, ta2_xz_z_xzzz_0, ta2_xz_z_yyyy_0, ta2_xz_z_yyyz_0, ta2_xz_z_yyzz_0, ta2_xz_z_yzzz_0, ta2_xz_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_z_xxxx_0[i] = ta1_x_0_xxxx_1[i] + ta2_xz_0_xxxx_0[i] * pa_z[i] - ta2_xz_0_xxxx_1[i] * pc_z[i];

        ta2_xz_z_xxxy_0[i] = ta1_x_0_xxxy_1[i] + ta2_xz_0_xxxy_0[i] * pa_z[i] - ta2_xz_0_xxxy_1[i] * pc_z[i];

        ta2_xz_z_xxxz_0[i] = ta2_xz_0_xxx_0[i] * fe_0 - ta2_xz_0_xxx_1[i] * fe_0 + ta1_x_0_xxxz_1[i] + ta2_xz_0_xxxz_0[i] * pa_z[i] - ta2_xz_0_xxxz_1[i] * pc_z[i];

        ta2_xz_z_xxyy_0[i] = ta1_x_0_xxyy_1[i] + ta2_xz_0_xxyy_0[i] * pa_z[i] - ta2_xz_0_xxyy_1[i] * pc_z[i];

        ta2_xz_z_xxyz_0[i] = ta2_xz_0_xxy_0[i] * fe_0 - ta2_xz_0_xxy_1[i] * fe_0 + ta1_x_0_xxyz_1[i] + ta2_xz_0_xxyz_0[i] * pa_z[i] - ta2_xz_0_xxyz_1[i] * pc_z[i];

        ta2_xz_z_xxzz_0[i] = 2.0 * ta2_xz_0_xxz_0[i] * fe_0 - 2.0 * ta2_xz_0_xxz_1[i] * fe_0 + ta1_x_0_xxzz_1[i] + ta2_xz_0_xxzz_0[i] * pa_z[i] - ta2_xz_0_xxzz_1[i] * pc_z[i];

        ta2_xz_z_xyyy_0[i] = ta1_x_0_xyyy_1[i] + ta2_xz_0_xyyy_0[i] * pa_z[i] - ta2_xz_0_xyyy_1[i] * pc_z[i];

        ta2_xz_z_xyyz_0[i] = ta2_xz_0_xyy_0[i] * fe_0 - ta2_xz_0_xyy_1[i] * fe_0 + ta1_x_0_xyyz_1[i] + ta2_xz_0_xyyz_0[i] * pa_z[i] - ta2_xz_0_xyyz_1[i] * pc_z[i];

        ta2_xz_z_xyzz_0[i] = 2.0 * ta2_xz_0_xyz_0[i] * fe_0 - 2.0 * ta2_xz_0_xyz_1[i] * fe_0 + ta1_x_0_xyzz_1[i] + ta2_xz_0_xyzz_0[i] * pa_z[i] - ta2_xz_0_xyzz_1[i] * pc_z[i];

        ta2_xz_z_xzzz_0[i] = 3.0 * ta2_xz_0_xzz_0[i] * fe_0 - 3.0 * ta2_xz_0_xzz_1[i] * fe_0 + ta1_x_0_xzzz_1[i] + ta2_xz_0_xzzz_0[i] * pa_z[i] - ta2_xz_0_xzzz_1[i] * pc_z[i];

        ta2_xz_z_yyyy_0[i] = ta1_x_0_yyyy_1[i] + ta2_xz_0_yyyy_0[i] * pa_z[i] - ta2_xz_0_yyyy_1[i] * pc_z[i];

        ta2_xz_z_yyyz_0[i] = ta2_xz_0_yyy_0[i] * fe_0 - ta2_xz_0_yyy_1[i] * fe_0 + ta1_x_0_yyyz_1[i] + ta2_xz_0_yyyz_0[i] * pa_z[i] - ta2_xz_0_yyyz_1[i] * pc_z[i];

        ta2_xz_z_yyzz_0[i] = 2.0 * ta2_xz_0_yyz_0[i] * fe_0 - 2.0 * ta2_xz_0_yyz_1[i] * fe_0 + ta1_x_0_yyzz_1[i] + ta2_xz_0_yyzz_0[i] * pa_z[i] - ta2_xz_0_yyzz_1[i] * pc_z[i];

        ta2_xz_z_yzzz_0[i] = 3.0 * ta2_xz_0_yzz_0[i] * fe_0 - 3.0 * ta2_xz_0_yzz_1[i] * fe_0 + ta1_x_0_yzzz_1[i] + ta2_xz_0_yzzz_0[i] * pa_z[i] - ta2_xz_0_yzzz_1[i] * pc_z[i];

        ta2_xz_z_zzzz_0[i] = 4.0 * ta2_xz_0_zzz_0[i] * fe_0 - 4.0 * ta2_xz_0_zzz_1[i] * fe_0 + ta1_x_0_zzzz_1[i] + ta2_xz_0_zzzz_0[i] * pa_z[i] - ta2_xz_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 135-150 components of targeted buffer : PG

    auto ta2_yy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 135);

    auto ta2_yy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 136);

    auto ta2_yy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 137);

    auto ta2_yy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 138);

    auto ta2_yy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 139);

    auto ta2_yy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 140);

    auto ta2_yy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 141);

    auto ta2_yy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 142);

    auto ta2_yy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 143);

    auto ta2_yy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 144);

    auto ta2_yy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 145);

    auto ta2_yy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 146);

    auto ta2_yy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 147);

    auto ta2_yy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 148);

    auto ta2_yy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 149);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_0_xxx_0, ta2_yy_0_xxx_1, ta2_yy_0_xxxx_0, ta2_yy_0_xxxx_1, ta2_yy_0_xxxy_0, ta2_yy_0_xxxy_1, ta2_yy_0_xxxz_0, ta2_yy_0_xxxz_1, ta2_yy_0_xxy_0, ta2_yy_0_xxy_1, ta2_yy_0_xxyy_0, ta2_yy_0_xxyy_1, ta2_yy_0_xxyz_0, ta2_yy_0_xxyz_1, ta2_yy_0_xxz_0, ta2_yy_0_xxz_1, ta2_yy_0_xxzz_0, ta2_yy_0_xxzz_1, ta2_yy_0_xyy_0, ta2_yy_0_xyy_1, ta2_yy_0_xyyy_0, ta2_yy_0_xyyy_1, ta2_yy_0_xyyz_0, ta2_yy_0_xyyz_1, ta2_yy_0_xyz_0, ta2_yy_0_xyz_1, ta2_yy_0_xyzz_0, ta2_yy_0_xyzz_1, ta2_yy_0_xzz_0, ta2_yy_0_xzz_1, ta2_yy_0_xzzz_0, ta2_yy_0_xzzz_1, ta2_yy_0_yyy_0, ta2_yy_0_yyy_1, ta2_yy_0_yyyy_0, ta2_yy_0_yyyy_1, ta2_yy_0_yyyz_0, ta2_yy_0_yyyz_1, ta2_yy_0_yyz_0, ta2_yy_0_yyz_1, ta2_yy_0_yyzz_0, ta2_yy_0_yyzz_1, ta2_yy_0_yzz_0, ta2_yy_0_yzz_1, ta2_yy_0_yzzz_0, ta2_yy_0_yzzz_1, ta2_yy_0_zzz_0, ta2_yy_0_zzz_1, ta2_yy_0_zzzz_0, ta2_yy_0_zzzz_1, ta2_yy_x_xxxx_0, ta2_yy_x_xxxy_0, ta2_yy_x_xxxz_0, ta2_yy_x_xxyy_0, ta2_yy_x_xxyz_0, ta2_yy_x_xxzz_0, ta2_yy_x_xyyy_0, ta2_yy_x_xyyz_0, ta2_yy_x_xyzz_0, ta2_yy_x_xzzz_0, ta2_yy_x_yyyy_0, ta2_yy_x_yyyz_0, ta2_yy_x_yyzz_0, ta2_yy_x_yzzz_0, ta2_yy_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_x_xxxx_0[i] = 4.0 * ta2_yy_0_xxx_0[i] * fe_0 - 4.0 * ta2_yy_0_xxx_1[i] * fe_0 + ta2_yy_0_xxxx_0[i] * pa_x[i] - ta2_yy_0_xxxx_1[i] * pc_x[i];

        ta2_yy_x_xxxy_0[i] = 3.0 * ta2_yy_0_xxy_0[i] * fe_0 - 3.0 * ta2_yy_0_xxy_1[i] * fe_0 + ta2_yy_0_xxxy_0[i] * pa_x[i] - ta2_yy_0_xxxy_1[i] * pc_x[i];

        ta2_yy_x_xxxz_0[i] = 3.0 * ta2_yy_0_xxz_0[i] * fe_0 - 3.0 * ta2_yy_0_xxz_1[i] * fe_0 + ta2_yy_0_xxxz_0[i] * pa_x[i] - ta2_yy_0_xxxz_1[i] * pc_x[i];

        ta2_yy_x_xxyy_0[i] = 2.0 * ta2_yy_0_xyy_0[i] * fe_0 - 2.0 * ta2_yy_0_xyy_1[i] * fe_0 + ta2_yy_0_xxyy_0[i] * pa_x[i] - ta2_yy_0_xxyy_1[i] * pc_x[i];

        ta2_yy_x_xxyz_0[i] = 2.0 * ta2_yy_0_xyz_0[i] * fe_0 - 2.0 * ta2_yy_0_xyz_1[i] * fe_0 + ta2_yy_0_xxyz_0[i] * pa_x[i] - ta2_yy_0_xxyz_1[i] * pc_x[i];

        ta2_yy_x_xxzz_0[i] = 2.0 * ta2_yy_0_xzz_0[i] * fe_0 - 2.0 * ta2_yy_0_xzz_1[i] * fe_0 + ta2_yy_0_xxzz_0[i] * pa_x[i] - ta2_yy_0_xxzz_1[i] * pc_x[i];

        ta2_yy_x_xyyy_0[i] = ta2_yy_0_yyy_0[i] * fe_0 - ta2_yy_0_yyy_1[i] * fe_0 + ta2_yy_0_xyyy_0[i] * pa_x[i] - ta2_yy_0_xyyy_1[i] * pc_x[i];

        ta2_yy_x_xyyz_0[i] = ta2_yy_0_yyz_0[i] * fe_0 - ta2_yy_0_yyz_1[i] * fe_0 + ta2_yy_0_xyyz_0[i] * pa_x[i] - ta2_yy_0_xyyz_1[i] * pc_x[i];

        ta2_yy_x_xyzz_0[i] = ta2_yy_0_yzz_0[i] * fe_0 - ta2_yy_0_yzz_1[i] * fe_0 + ta2_yy_0_xyzz_0[i] * pa_x[i] - ta2_yy_0_xyzz_1[i] * pc_x[i];

        ta2_yy_x_xzzz_0[i] = ta2_yy_0_zzz_0[i] * fe_0 - ta2_yy_0_zzz_1[i] * fe_0 + ta2_yy_0_xzzz_0[i] * pa_x[i] - ta2_yy_0_xzzz_1[i] * pc_x[i];

        ta2_yy_x_yyyy_0[i] = ta2_yy_0_yyyy_0[i] * pa_x[i] - ta2_yy_0_yyyy_1[i] * pc_x[i];

        ta2_yy_x_yyyz_0[i] = ta2_yy_0_yyyz_0[i] * pa_x[i] - ta2_yy_0_yyyz_1[i] * pc_x[i];

        ta2_yy_x_yyzz_0[i] = ta2_yy_0_yyzz_0[i] * pa_x[i] - ta2_yy_0_yyzz_1[i] * pc_x[i];

        ta2_yy_x_yzzz_0[i] = ta2_yy_0_yzzz_0[i] * pa_x[i] - ta2_yy_0_yzzz_1[i] * pc_x[i];

        ta2_yy_x_zzzz_0[i] = ta2_yy_0_zzzz_0[i] * pa_x[i] - ta2_yy_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : PG

    auto ta2_yy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 150);

    auto ta2_yy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 151);

    auto ta2_yy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 152);

    auto ta2_yy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 153);

    auto ta2_yy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 154);

    auto ta2_yy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 155);

    auto ta2_yy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 156);

    auto ta2_yy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 157);

    auto ta2_yy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 158);

    auto ta2_yy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 159);

    auto ta2_yy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 160);

    auto ta2_yy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 161);

    auto ta2_yy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 162);

    auto ta2_yy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 163);

    auto ta2_yy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 164);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_0_xxxx_1, ta1_y_0_xxxy_1, ta1_y_0_xxxz_1, ta1_y_0_xxyy_1, ta1_y_0_xxyz_1, ta1_y_0_xxzz_1, ta1_y_0_xyyy_1, ta1_y_0_xyyz_1, ta1_y_0_xyzz_1, ta1_y_0_xzzz_1, ta1_y_0_yyyy_1, ta1_y_0_yyyz_1, ta1_y_0_yyzz_1, ta1_y_0_yzzz_1, ta1_y_0_zzzz_1, ta2_yy_0_xxx_0, ta2_yy_0_xxx_1, ta2_yy_0_xxxx_0, ta2_yy_0_xxxx_1, ta2_yy_0_xxxy_0, ta2_yy_0_xxxy_1, ta2_yy_0_xxxz_0, ta2_yy_0_xxxz_1, ta2_yy_0_xxy_0, ta2_yy_0_xxy_1, ta2_yy_0_xxyy_0, ta2_yy_0_xxyy_1, ta2_yy_0_xxyz_0, ta2_yy_0_xxyz_1, ta2_yy_0_xxz_0, ta2_yy_0_xxz_1, ta2_yy_0_xxzz_0, ta2_yy_0_xxzz_1, ta2_yy_0_xyy_0, ta2_yy_0_xyy_1, ta2_yy_0_xyyy_0, ta2_yy_0_xyyy_1, ta2_yy_0_xyyz_0, ta2_yy_0_xyyz_1, ta2_yy_0_xyz_0, ta2_yy_0_xyz_1, ta2_yy_0_xyzz_0, ta2_yy_0_xyzz_1, ta2_yy_0_xzz_0, ta2_yy_0_xzz_1, ta2_yy_0_xzzz_0, ta2_yy_0_xzzz_1, ta2_yy_0_yyy_0, ta2_yy_0_yyy_1, ta2_yy_0_yyyy_0, ta2_yy_0_yyyy_1, ta2_yy_0_yyyz_0, ta2_yy_0_yyyz_1, ta2_yy_0_yyz_0, ta2_yy_0_yyz_1, ta2_yy_0_yyzz_0, ta2_yy_0_yyzz_1, ta2_yy_0_yzz_0, ta2_yy_0_yzz_1, ta2_yy_0_yzzz_0, ta2_yy_0_yzzz_1, ta2_yy_0_zzz_0, ta2_yy_0_zzz_1, ta2_yy_0_zzzz_0, ta2_yy_0_zzzz_1, ta2_yy_y_xxxx_0, ta2_yy_y_xxxy_0, ta2_yy_y_xxxz_0, ta2_yy_y_xxyy_0, ta2_yy_y_xxyz_0, ta2_yy_y_xxzz_0, ta2_yy_y_xyyy_0, ta2_yy_y_xyyz_0, ta2_yy_y_xyzz_0, ta2_yy_y_xzzz_0, ta2_yy_y_yyyy_0, ta2_yy_y_yyyz_0, ta2_yy_y_yyzz_0, ta2_yy_y_yzzz_0, ta2_yy_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_y_xxxx_0[i] = 2.0 * ta1_y_0_xxxx_1[i] + ta2_yy_0_xxxx_0[i] * pa_y[i] - ta2_yy_0_xxxx_1[i] * pc_y[i];

        ta2_yy_y_xxxy_0[i] = ta2_yy_0_xxx_0[i] * fe_0 - ta2_yy_0_xxx_1[i] * fe_0 + 2.0 * ta1_y_0_xxxy_1[i] + ta2_yy_0_xxxy_0[i] * pa_y[i] - ta2_yy_0_xxxy_1[i] * pc_y[i];

        ta2_yy_y_xxxz_0[i] = 2.0 * ta1_y_0_xxxz_1[i] + ta2_yy_0_xxxz_0[i] * pa_y[i] - ta2_yy_0_xxxz_1[i] * pc_y[i];

        ta2_yy_y_xxyy_0[i] = 2.0 * ta2_yy_0_xxy_0[i] * fe_0 - 2.0 * ta2_yy_0_xxy_1[i] * fe_0 + 2.0 * ta1_y_0_xxyy_1[i] + ta2_yy_0_xxyy_0[i] * pa_y[i] - ta2_yy_0_xxyy_1[i] * pc_y[i];

        ta2_yy_y_xxyz_0[i] = ta2_yy_0_xxz_0[i] * fe_0 - ta2_yy_0_xxz_1[i] * fe_0 + 2.0 * ta1_y_0_xxyz_1[i] + ta2_yy_0_xxyz_0[i] * pa_y[i] - ta2_yy_0_xxyz_1[i] * pc_y[i];

        ta2_yy_y_xxzz_0[i] = 2.0 * ta1_y_0_xxzz_1[i] + ta2_yy_0_xxzz_0[i] * pa_y[i] - ta2_yy_0_xxzz_1[i] * pc_y[i];

        ta2_yy_y_xyyy_0[i] = 3.0 * ta2_yy_0_xyy_0[i] * fe_0 - 3.0 * ta2_yy_0_xyy_1[i] * fe_0 + 2.0 * ta1_y_0_xyyy_1[i] + ta2_yy_0_xyyy_0[i] * pa_y[i] - ta2_yy_0_xyyy_1[i] * pc_y[i];

        ta2_yy_y_xyyz_0[i] = 2.0 * ta2_yy_0_xyz_0[i] * fe_0 - 2.0 * ta2_yy_0_xyz_1[i] * fe_0 + 2.0 * ta1_y_0_xyyz_1[i] + ta2_yy_0_xyyz_0[i] * pa_y[i] - ta2_yy_0_xyyz_1[i] * pc_y[i];

        ta2_yy_y_xyzz_0[i] = ta2_yy_0_xzz_0[i] * fe_0 - ta2_yy_0_xzz_1[i] * fe_0 + 2.0 * ta1_y_0_xyzz_1[i] + ta2_yy_0_xyzz_0[i] * pa_y[i] - ta2_yy_0_xyzz_1[i] * pc_y[i];

        ta2_yy_y_xzzz_0[i] = 2.0 * ta1_y_0_xzzz_1[i] + ta2_yy_0_xzzz_0[i] * pa_y[i] - ta2_yy_0_xzzz_1[i] * pc_y[i];

        ta2_yy_y_yyyy_0[i] = 4.0 * ta2_yy_0_yyy_0[i] * fe_0 - 4.0 * ta2_yy_0_yyy_1[i] * fe_0 + 2.0 * ta1_y_0_yyyy_1[i] + ta2_yy_0_yyyy_0[i] * pa_y[i] - ta2_yy_0_yyyy_1[i] * pc_y[i];

        ta2_yy_y_yyyz_0[i] = 3.0 * ta2_yy_0_yyz_0[i] * fe_0 - 3.0 * ta2_yy_0_yyz_1[i] * fe_0 + 2.0 * ta1_y_0_yyyz_1[i] + ta2_yy_0_yyyz_0[i] * pa_y[i] - ta2_yy_0_yyyz_1[i] * pc_y[i];

        ta2_yy_y_yyzz_0[i] = 2.0 * ta2_yy_0_yzz_0[i] * fe_0 - 2.0 * ta2_yy_0_yzz_1[i] * fe_0 + 2.0 * ta1_y_0_yyzz_1[i] + ta2_yy_0_yyzz_0[i] * pa_y[i] - ta2_yy_0_yyzz_1[i] * pc_y[i];

        ta2_yy_y_yzzz_0[i] = ta2_yy_0_zzz_0[i] * fe_0 - ta2_yy_0_zzz_1[i] * fe_0 + 2.0 * ta1_y_0_yzzz_1[i] + ta2_yy_0_yzzz_0[i] * pa_y[i] - ta2_yy_0_yzzz_1[i] * pc_y[i];

        ta2_yy_y_zzzz_0[i] = 2.0 * ta1_y_0_zzzz_1[i] + ta2_yy_0_zzzz_0[i] * pa_y[i] - ta2_yy_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : PG

    auto ta2_yy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 165);

    auto ta2_yy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 166);

    auto ta2_yy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 167);

    auto ta2_yy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 168);

    auto ta2_yy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 169);

    auto ta2_yy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 170);

    auto ta2_yy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 171);

    auto ta2_yy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 172);

    auto ta2_yy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 173);

    auto ta2_yy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 174);

    auto ta2_yy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 175);

    auto ta2_yy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 176);

    auto ta2_yy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 177);

    auto ta2_yy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 178);

    auto ta2_yy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 179);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_0_xxx_0, ta2_yy_0_xxx_1, ta2_yy_0_xxxx_0, ta2_yy_0_xxxx_1, ta2_yy_0_xxxy_0, ta2_yy_0_xxxy_1, ta2_yy_0_xxxz_0, ta2_yy_0_xxxz_1, ta2_yy_0_xxy_0, ta2_yy_0_xxy_1, ta2_yy_0_xxyy_0, ta2_yy_0_xxyy_1, ta2_yy_0_xxyz_0, ta2_yy_0_xxyz_1, ta2_yy_0_xxz_0, ta2_yy_0_xxz_1, ta2_yy_0_xxzz_0, ta2_yy_0_xxzz_1, ta2_yy_0_xyy_0, ta2_yy_0_xyy_1, ta2_yy_0_xyyy_0, ta2_yy_0_xyyy_1, ta2_yy_0_xyyz_0, ta2_yy_0_xyyz_1, ta2_yy_0_xyz_0, ta2_yy_0_xyz_1, ta2_yy_0_xyzz_0, ta2_yy_0_xyzz_1, ta2_yy_0_xzz_0, ta2_yy_0_xzz_1, ta2_yy_0_xzzz_0, ta2_yy_0_xzzz_1, ta2_yy_0_yyy_0, ta2_yy_0_yyy_1, ta2_yy_0_yyyy_0, ta2_yy_0_yyyy_1, ta2_yy_0_yyyz_0, ta2_yy_0_yyyz_1, ta2_yy_0_yyz_0, ta2_yy_0_yyz_1, ta2_yy_0_yyzz_0, ta2_yy_0_yyzz_1, ta2_yy_0_yzz_0, ta2_yy_0_yzz_1, ta2_yy_0_yzzz_0, ta2_yy_0_yzzz_1, ta2_yy_0_zzz_0, ta2_yy_0_zzz_1, ta2_yy_0_zzzz_0, ta2_yy_0_zzzz_1, ta2_yy_z_xxxx_0, ta2_yy_z_xxxy_0, ta2_yy_z_xxxz_0, ta2_yy_z_xxyy_0, ta2_yy_z_xxyz_0, ta2_yy_z_xxzz_0, ta2_yy_z_xyyy_0, ta2_yy_z_xyyz_0, ta2_yy_z_xyzz_0, ta2_yy_z_xzzz_0, ta2_yy_z_yyyy_0, ta2_yy_z_yyyz_0, ta2_yy_z_yyzz_0, ta2_yy_z_yzzz_0, ta2_yy_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_z_xxxx_0[i] = ta2_yy_0_xxxx_0[i] * pa_z[i] - ta2_yy_0_xxxx_1[i] * pc_z[i];

        ta2_yy_z_xxxy_0[i] = ta2_yy_0_xxxy_0[i] * pa_z[i] - ta2_yy_0_xxxy_1[i] * pc_z[i];

        ta2_yy_z_xxxz_0[i] = ta2_yy_0_xxx_0[i] * fe_0 - ta2_yy_0_xxx_1[i] * fe_0 + ta2_yy_0_xxxz_0[i] * pa_z[i] - ta2_yy_0_xxxz_1[i] * pc_z[i];

        ta2_yy_z_xxyy_0[i] = ta2_yy_0_xxyy_0[i] * pa_z[i] - ta2_yy_0_xxyy_1[i] * pc_z[i];

        ta2_yy_z_xxyz_0[i] = ta2_yy_0_xxy_0[i] * fe_0 - ta2_yy_0_xxy_1[i] * fe_0 + ta2_yy_0_xxyz_0[i] * pa_z[i] - ta2_yy_0_xxyz_1[i] * pc_z[i];

        ta2_yy_z_xxzz_0[i] = 2.0 * ta2_yy_0_xxz_0[i] * fe_0 - 2.0 * ta2_yy_0_xxz_1[i] * fe_0 + ta2_yy_0_xxzz_0[i] * pa_z[i] - ta2_yy_0_xxzz_1[i] * pc_z[i];

        ta2_yy_z_xyyy_0[i] = ta2_yy_0_xyyy_0[i] * pa_z[i] - ta2_yy_0_xyyy_1[i] * pc_z[i];

        ta2_yy_z_xyyz_0[i] = ta2_yy_0_xyy_0[i] * fe_0 - ta2_yy_0_xyy_1[i] * fe_0 + ta2_yy_0_xyyz_0[i] * pa_z[i] - ta2_yy_0_xyyz_1[i] * pc_z[i];

        ta2_yy_z_xyzz_0[i] = 2.0 * ta2_yy_0_xyz_0[i] * fe_0 - 2.0 * ta2_yy_0_xyz_1[i] * fe_0 + ta2_yy_0_xyzz_0[i] * pa_z[i] - ta2_yy_0_xyzz_1[i] * pc_z[i];

        ta2_yy_z_xzzz_0[i] = 3.0 * ta2_yy_0_xzz_0[i] * fe_0 - 3.0 * ta2_yy_0_xzz_1[i] * fe_0 + ta2_yy_0_xzzz_0[i] * pa_z[i] - ta2_yy_0_xzzz_1[i] * pc_z[i];

        ta2_yy_z_yyyy_0[i] = ta2_yy_0_yyyy_0[i] * pa_z[i] - ta2_yy_0_yyyy_1[i] * pc_z[i];

        ta2_yy_z_yyyz_0[i] = ta2_yy_0_yyy_0[i] * fe_0 - ta2_yy_0_yyy_1[i] * fe_0 + ta2_yy_0_yyyz_0[i] * pa_z[i] - ta2_yy_0_yyyz_1[i] * pc_z[i];

        ta2_yy_z_yyzz_0[i] = 2.0 * ta2_yy_0_yyz_0[i] * fe_0 - 2.0 * ta2_yy_0_yyz_1[i] * fe_0 + ta2_yy_0_yyzz_0[i] * pa_z[i] - ta2_yy_0_yyzz_1[i] * pc_z[i];

        ta2_yy_z_yzzz_0[i] = 3.0 * ta2_yy_0_yzz_0[i] * fe_0 - 3.0 * ta2_yy_0_yzz_1[i] * fe_0 + ta2_yy_0_yzzz_0[i] * pa_z[i] - ta2_yy_0_yzzz_1[i] * pc_z[i];

        ta2_yy_z_zzzz_0[i] = 4.0 * ta2_yy_0_zzz_0[i] * fe_0 - 4.0 * ta2_yy_0_zzz_1[i] * fe_0 + ta2_yy_0_zzzz_0[i] * pa_z[i] - ta2_yy_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 180-195 components of targeted buffer : PG

    auto ta2_yz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 180);

    auto ta2_yz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 181);

    auto ta2_yz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 182);

    auto ta2_yz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 183);

    auto ta2_yz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 184);

    auto ta2_yz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 185);

    auto ta2_yz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 186);

    auto ta2_yz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 187);

    auto ta2_yz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 188);

    auto ta2_yz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 189);

    auto ta2_yz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 190);

    auto ta2_yz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 191);

    auto ta2_yz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 192);

    auto ta2_yz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 193);

    auto ta2_yz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 194);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_0_xxx_0, ta2_yz_0_xxx_1, ta2_yz_0_xxxx_0, ta2_yz_0_xxxx_1, ta2_yz_0_xxxy_0, ta2_yz_0_xxxy_1, ta2_yz_0_xxxz_0, ta2_yz_0_xxxz_1, ta2_yz_0_xxy_0, ta2_yz_0_xxy_1, ta2_yz_0_xxyy_0, ta2_yz_0_xxyy_1, ta2_yz_0_xxyz_0, ta2_yz_0_xxyz_1, ta2_yz_0_xxz_0, ta2_yz_0_xxz_1, ta2_yz_0_xxzz_0, ta2_yz_0_xxzz_1, ta2_yz_0_xyy_0, ta2_yz_0_xyy_1, ta2_yz_0_xyyy_0, ta2_yz_0_xyyy_1, ta2_yz_0_xyyz_0, ta2_yz_0_xyyz_1, ta2_yz_0_xyz_0, ta2_yz_0_xyz_1, ta2_yz_0_xyzz_0, ta2_yz_0_xyzz_1, ta2_yz_0_xzz_0, ta2_yz_0_xzz_1, ta2_yz_0_xzzz_0, ta2_yz_0_xzzz_1, ta2_yz_0_yyy_0, ta2_yz_0_yyy_1, ta2_yz_0_yyyy_0, ta2_yz_0_yyyy_1, ta2_yz_0_yyyz_0, ta2_yz_0_yyyz_1, ta2_yz_0_yyz_0, ta2_yz_0_yyz_1, ta2_yz_0_yyzz_0, ta2_yz_0_yyzz_1, ta2_yz_0_yzz_0, ta2_yz_0_yzz_1, ta2_yz_0_yzzz_0, ta2_yz_0_yzzz_1, ta2_yz_0_zzz_0, ta2_yz_0_zzz_1, ta2_yz_0_zzzz_0, ta2_yz_0_zzzz_1, ta2_yz_x_xxxx_0, ta2_yz_x_xxxy_0, ta2_yz_x_xxxz_0, ta2_yz_x_xxyy_0, ta2_yz_x_xxyz_0, ta2_yz_x_xxzz_0, ta2_yz_x_xyyy_0, ta2_yz_x_xyyz_0, ta2_yz_x_xyzz_0, ta2_yz_x_xzzz_0, ta2_yz_x_yyyy_0, ta2_yz_x_yyyz_0, ta2_yz_x_yyzz_0, ta2_yz_x_yzzz_0, ta2_yz_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_x_xxxx_0[i] = 4.0 * ta2_yz_0_xxx_0[i] * fe_0 - 4.0 * ta2_yz_0_xxx_1[i] * fe_0 + ta2_yz_0_xxxx_0[i] * pa_x[i] - ta2_yz_0_xxxx_1[i] * pc_x[i];

        ta2_yz_x_xxxy_0[i] = 3.0 * ta2_yz_0_xxy_0[i] * fe_0 - 3.0 * ta2_yz_0_xxy_1[i] * fe_0 + ta2_yz_0_xxxy_0[i] * pa_x[i] - ta2_yz_0_xxxy_1[i] * pc_x[i];

        ta2_yz_x_xxxz_0[i] = 3.0 * ta2_yz_0_xxz_0[i] * fe_0 - 3.0 * ta2_yz_0_xxz_1[i] * fe_0 + ta2_yz_0_xxxz_0[i] * pa_x[i] - ta2_yz_0_xxxz_1[i] * pc_x[i];

        ta2_yz_x_xxyy_0[i] = 2.0 * ta2_yz_0_xyy_0[i] * fe_0 - 2.0 * ta2_yz_0_xyy_1[i] * fe_0 + ta2_yz_0_xxyy_0[i] * pa_x[i] - ta2_yz_0_xxyy_1[i] * pc_x[i];

        ta2_yz_x_xxyz_0[i] = 2.0 * ta2_yz_0_xyz_0[i] * fe_0 - 2.0 * ta2_yz_0_xyz_1[i] * fe_0 + ta2_yz_0_xxyz_0[i] * pa_x[i] - ta2_yz_0_xxyz_1[i] * pc_x[i];

        ta2_yz_x_xxzz_0[i] = 2.0 * ta2_yz_0_xzz_0[i] * fe_0 - 2.0 * ta2_yz_0_xzz_1[i] * fe_0 + ta2_yz_0_xxzz_0[i] * pa_x[i] - ta2_yz_0_xxzz_1[i] * pc_x[i];

        ta2_yz_x_xyyy_0[i] = ta2_yz_0_yyy_0[i] * fe_0 - ta2_yz_0_yyy_1[i] * fe_0 + ta2_yz_0_xyyy_0[i] * pa_x[i] - ta2_yz_0_xyyy_1[i] * pc_x[i];

        ta2_yz_x_xyyz_0[i] = ta2_yz_0_yyz_0[i] * fe_0 - ta2_yz_0_yyz_1[i] * fe_0 + ta2_yz_0_xyyz_0[i] * pa_x[i] - ta2_yz_0_xyyz_1[i] * pc_x[i];

        ta2_yz_x_xyzz_0[i] = ta2_yz_0_yzz_0[i] * fe_0 - ta2_yz_0_yzz_1[i] * fe_0 + ta2_yz_0_xyzz_0[i] * pa_x[i] - ta2_yz_0_xyzz_1[i] * pc_x[i];

        ta2_yz_x_xzzz_0[i] = ta2_yz_0_zzz_0[i] * fe_0 - ta2_yz_0_zzz_1[i] * fe_0 + ta2_yz_0_xzzz_0[i] * pa_x[i] - ta2_yz_0_xzzz_1[i] * pc_x[i];

        ta2_yz_x_yyyy_0[i] = ta2_yz_0_yyyy_0[i] * pa_x[i] - ta2_yz_0_yyyy_1[i] * pc_x[i];

        ta2_yz_x_yyyz_0[i] = ta2_yz_0_yyyz_0[i] * pa_x[i] - ta2_yz_0_yyyz_1[i] * pc_x[i];

        ta2_yz_x_yyzz_0[i] = ta2_yz_0_yyzz_0[i] * pa_x[i] - ta2_yz_0_yyzz_1[i] * pc_x[i];

        ta2_yz_x_yzzz_0[i] = ta2_yz_0_yzzz_0[i] * pa_x[i] - ta2_yz_0_yzzz_1[i] * pc_x[i];

        ta2_yz_x_zzzz_0[i] = ta2_yz_0_zzzz_0[i] * pa_x[i] - ta2_yz_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : PG

    auto ta2_yz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 195);

    auto ta2_yz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 196);

    auto ta2_yz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 197);

    auto ta2_yz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 198);

    auto ta2_yz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 199);

    auto ta2_yz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 200);

    auto ta2_yz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 201);

    auto ta2_yz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 202);

    auto ta2_yz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 203);

    auto ta2_yz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 204);

    auto ta2_yz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 205);

    auto ta2_yz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 206);

    auto ta2_yz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 207);

    auto ta2_yz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 208);

    auto ta2_yz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 209);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_0_xxxx_1, ta1_z_0_xxxy_1, ta1_z_0_xxxz_1, ta1_z_0_xxyy_1, ta1_z_0_xxyz_1, ta1_z_0_xxzz_1, ta1_z_0_xyyy_1, ta1_z_0_xyyz_1, ta1_z_0_xyzz_1, ta1_z_0_xzzz_1, ta1_z_0_yyyy_1, ta1_z_0_yyyz_1, ta1_z_0_yyzz_1, ta1_z_0_yzzz_1, ta1_z_0_zzzz_1, ta2_yz_0_xxx_0, ta2_yz_0_xxx_1, ta2_yz_0_xxxx_0, ta2_yz_0_xxxx_1, ta2_yz_0_xxxy_0, ta2_yz_0_xxxy_1, ta2_yz_0_xxxz_0, ta2_yz_0_xxxz_1, ta2_yz_0_xxy_0, ta2_yz_0_xxy_1, ta2_yz_0_xxyy_0, ta2_yz_0_xxyy_1, ta2_yz_0_xxyz_0, ta2_yz_0_xxyz_1, ta2_yz_0_xxz_0, ta2_yz_0_xxz_1, ta2_yz_0_xxzz_0, ta2_yz_0_xxzz_1, ta2_yz_0_xyy_0, ta2_yz_0_xyy_1, ta2_yz_0_xyyy_0, ta2_yz_0_xyyy_1, ta2_yz_0_xyyz_0, ta2_yz_0_xyyz_1, ta2_yz_0_xyz_0, ta2_yz_0_xyz_1, ta2_yz_0_xyzz_0, ta2_yz_0_xyzz_1, ta2_yz_0_xzz_0, ta2_yz_0_xzz_1, ta2_yz_0_xzzz_0, ta2_yz_0_xzzz_1, ta2_yz_0_yyy_0, ta2_yz_0_yyy_1, ta2_yz_0_yyyy_0, ta2_yz_0_yyyy_1, ta2_yz_0_yyyz_0, ta2_yz_0_yyyz_1, ta2_yz_0_yyz_0, ta2_yz_0_yyz_1, ta2_yz_0_yyzz_0, ta2_yz_0_yyzz_1, ta2_yz_0_yzz_0, ta2_yz_0_yzz_1, ta2_yz_0_yzzz_0, ta2_yz_0_yzzz_1, ta2_yz_0_zzz_0, ta2_yz_0_zzz_1, ta2_yz_0_zzzz_0, ta2_yz_0_zzzz_1, ta2_yz_y_xxxx_0, ta2_yz_y_xxxy_0, ta2_yz_y_xxxz_0, ta2_yz_y_xxyy_0, ta2_yz_y_xxyz_0, ta2_yz_y_xxzz_0, ta2_yz_y_xyyy_0, ta2_yz_y_xyyz_0, ta2_yz_y_xyzz_0, ta2_yz_y_xzzz_0, ta2_yz_y_yyyy_0, ta2_yz_y_yyyz_0, ta2_yz_y_yyzz_0, ta2_yz_y_yzzz_0, ta2_yz_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_y_xxxx_0[i] = ta1_z_0_xxxx_1[i] + ta2_yz_0_xxxx_0[i] * pa_y[i] - ta2_yz_0_xxxx_1[i] * pc_y[i];

        ta2_yz_y_xxxy_0[i] = ta2_yz_0_xxx_0[i] * fe_0 - ta2_yz_0_xxx_1[i] * fe_0 + ta1_z_0_xxxy_1[i] + ta2_yz_0_xxxy_0[i] * pa_y[i] - ta2_yz_0_xxxy_1[i] * pc_y[i];

        ta2_yz_y_xxxz_0[i] = ta1_z_0_xxxz_1[i] + ta2_yz_0_xxxz_0[i] * pa_y[i] - ta2_yz_0_xxxz_1[i] * pc_y[i];

        ta2_yz_y_xxyy_0[i] = 2.0 * ta2_yz_0_xxy_0[i] * fe_0 - 2.0 * ta2_yz_0_xxy_1[i] * fe_0 + ta1_z_0_xxyy_1[i] + ta2_yz_0_xxyy_0[i] * pa_y[i] - ta2_yz_0_xxyy_1[i] * pc_y[i];

        ta2_yz_y_xxyz_0[i] = ta2_yz_0_xxz_0[i] * fe_0 - ta2_yz_0_xxz_1[i] * fe_0 + ta1_z_0_xxyz_1[i] + ta2_yz_0_xxyz_0[i] * pa_y[i] - ta2_yz_0_xxyz_1[i] * pc_y[i];

        ta2_yz_y_xxzz_0[i] = ta1_z_0_xxzz_1[i] + ta2_yz_0_xxzz_0[i] * pa_y[i] - ta2_yz_0_xxzz_1[i] * pc_y[i];

        ta2_yz_y_xyyy_0[i] = 3.0 * ta2_yz_0_xyy_0[i] * fe_0 - 3.0 * ta2_yz_0_xyy_1[i] * fe_0 + ta1_z_0_xyyy_1[i] + ta2_yz_0_xyyy_0[i] * pa_y[i] - ta2_yz_0_xyyy_1[i] * pc_y[i];

        ta2_yz_y_xyyz_0[i] = 2.0 * ta2_yz_0_xyz_0[i] * fe_0 - 2.0 * ta2_yz_0_xyz_1[i] * fe_0 + ta1_z_0_xyyz_1[i] + ta2_yz_0_xyyz_0[i] * pa_y[i] - ta2_yz_0_xyyz_1[i] * pc_y[i];

        ta2_yz_y_xyzz_0[i] = ta2_yz_0_xzz_0[i] * fe_0 - ta2_yz_0_xzz_1[i] * fe_0 + ta1_z_0_xyzz_1[i] + ta2_yz_0_xyzz_0[i] * pa_y[i] - ta2_yz_0_xyzz_1[i] * pc_y[i];

        ta2_yz_y_xzzz_0[i] = ta1_z_0_xzzz_1[i] + ta2_yz_0_xzzz_0[i] * pa_y[i] - ta2_yz_0_xzzz_1[i] * pc_y[i];

        ta2_yz_y_yyyy_0[i] = 4.0 * ta2_yz_0_yyy_0[i] * fe_0 - 4.0 * ta2_yz_0_yyy_1[i] * fe_0 + ta1_z_0_yyyy_1[i] + ta2_yz_0_yyyy_0[i] * pa_y[i] - ta2_yz_0_yyyy_1[i] * pc_y[i];

        ta2_yz_y_yyyz_0[i] = 3.0 * ta2_yz_0_yyz_0[i] * fe_0 - 3.0 * ta2_yz_0_yyz_1[i] * fe_0 + ta1_z_0_yyyz_1[i] + ta2_yz_0_yyyz_0[i] * pa_y[i] - ta2_yz_0_yyyz_1[i] * pc_y[i];

        ta2_yz_y_yyzz_0[i] = 2.0 * ta2_yz_0_yzz_0[i] * fe_0 - 2.0 * ta2_yz_0_yzz_1[i] * fe_0 + ta1_z_0_yyzz_1[i] + ta2_yz_0_yyzz_0[i] * pa_y[i] - ta2_yz_0_yyzz_1[i] * pc_y[i];

        ta2_yz_y_yzzz_0[i] = ta2_yz_0_zzz_0[i] * fe_0 - ta2_yz_0_zzz_1[i] * fe_0 + ta1_z_0_yzzz_1[i] + ta2_yz_0_yzzz_0[i] * pa_y[i] - ta2_yz_0_yzzz_1[i] * pc_y[i];

        ta2_yz_y_zzzz_0[i] = ta1_z_0_zzzz_1[i] + ta2_yz_0_zzzz_0[i] * pa_y[i] - ta2_yz_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : PG

    auto ta2_yz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 210);

    auto ta2_yz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 211);

    auto ta2_yz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 212);

    auto ta2_yz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 213);

    auto ta2_yz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 214);

    auto ta2_yz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 215);

    auto ta2_yz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 216);

    auto ta2_yz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 217);

    auto ta2_yz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 218);

    auto ta2_yz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 219);

    auto ta2_yz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 220);

    auto ta2_yz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 221);

    auto ta2_yz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 222);

    auto ta2_yz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 223);

    auto ta2_yz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 224);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_0_xxxx_1, ta1_y_0_xxxy_1, ta1_y_0_xxxz_1, ta1_y_0_xxyy_1, ta1_y_0_xxyz_1, ta1_y_0_xxzz_1, ta1_y_0_xyyy_1, ta1_y_0_xyyz_1, ta1_y_0_xyzz_1, ta1_y_0_xzzz_1, ta1_y_0_yyyy_1, ta1_y_0_yyyz_1, ta1_y_0_yyzz_1, ta1_y_0_yzzz_1, ta1_y_0_zzzz_1, ta2_yz_0_xxx_0, ta2_yz_0_xxx_1, ta2_yz_0_xxxx_0, ta2_yz_0_xxxx_1, ta2_yz_0_xxxy_0, ta2_yz_0_xxxy_1, ta2_yz_0_xxxz_0, ta2_yz_0_xxxz_1, ta2_yz_0_xxy_0, ta2_yz_0_xxy_1, ta2_yz_0_xxyy_0, ta2_yz_0_xxyy_1, ta2_yz_0_xxyz_0, ta2_yz_0_xxyz_1, ta2_yz_0_xxz_0, ta2_yz_0_xxz_1, ta2_yz_0_xxzz_0, ta2_yz_0_xxzz_1, ta2_yz_0_xyy_0, ta2_yz_0_xyy_1, ta2_yz_0_xyyy_0, ta2_yz_0_xyyy_1, ta2_yz_0_xyyz_0, ta2_yz_0_xyyz_1, ta2_yz_0_xyz_0, ta2_yz_0_xyz_1, ta2_yz_0_xyzz_0, ta2_yz_0_xyzz_1, ta2_yz_0_xzz_0, ta2_yz_0_xzz_1, ta2_yz_0_xzzz_0, ta2_yz_0_xzzz_1, ta2_yz_0_yyy_0, ta2_yz_0_yyy_1, ta2_yz_0_yyyy_0, ta2_yz_0_yyyy_1, ta2_yz_0_yyyz_0, ta2_yz_0_yyyz_1, ta2_yz_0_yyz_0, ta2_yz_0_yyz_1, ta2_yz_0_yyzz_0, ta2_yz_0_yyzz_1, ta2_yz_0_yzz_0, ta2_yz_0_yzz_1, ta2_yz_0_yzzz_0, ta2_yz_0_yzzz_1, ta2_yz_0_zzz_0, ta2_yz_0_zzz_1, ta2_yz_0_zzzz_0, ta2_yz_0_zzzz_1, ta2_yz_z_xxxx_0, ta2_yz_z_xxxy_0, ta2_yz_z_xxxz_0, ta2_yz_z_xxyy_0, ta2_yz_z_xxyz_0, ta2_yz_z_xxzz_0, ta2_yz_z_xyyy_0, ta2_yz_z_xyyz_0, ta2_yz_z_xyzz_0, ta2_yz_z_xzzz_0, ta2_yz_z_yyyy_0, ta2_yz_z_yyyz_0, ta2_yz_z_yyzz_0, ta2_yz_z_yzzz_0, ta2_yz_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_z_xxxx_0[i] = ta1_y_0_xxxx_1[i] + ta2_yz_0_xxxx_0[i] * pa_z[i] - ta2_yz_0_xxxx_1[i] * pc_z[i];

        ta2_yz_z_xxxy_0[i] = ta1_y_0_xxxy_1[i] + ta2_yz_0_xxxy_0[i] * pa_z[i] - ta2_yz_0_xxxy_1[i] * pc_z[i];

        ta2_yz_z_xxxz_0[i] = ta2_yz_0_xxx_0[i] * fe_0 - ta2_yz_0_xxx_1[i] * fe_0 + ta1_y_0_xxxz_1[i] + ta2_yz_0_xxxz_0[i] * pa_z[i] - ta2_yz_0_xxxz_1[i] * pc_z[i];

        ta2_yz_z_xxyy_0[i] = ta1_y_0_xxyy_1[i] + ta2_yz_0_xxyy_0[i] * pa_z[i] - ta2_yz_0_xxyy_1[i] * pc_z[i];

        ta2_yz_z_xxyz_0[i] = ta2_yz_0_xxy_0[i] * fe_0 - ta2_yz_0_xxy_1[i] * fe_0 + ta1_y_0_xxyz_1[i] + ta2_yz_0_xxyz_0[i] * pa_z[i] - ta2_yz_0_xxyz_1[i] * pc_z[i];

        ta2_yz_z_xxzz_0[i] = 2.0 * ta2_yz_0_xxz_0[i] * fe_0 - 2.0 * ta2_yz_0_xxz_1[i] * fe_0 + ta1_y_0_xxzz_1[i] + ta2_yz_0_xxzz_0[i] * pa_z[i] - ta2_yz_0_xxzz_1[i] * pc_z[i];

        ta2_yz_z_xyyy_0[i] = ta1_y_0_xyyy_1[i] + ta2_yz_0_xyyy_0[i] * pa_z[i] - ta2_yz_0_xyyy_1[i] * pc_z[i];

        ta2_yz_z_xyyz_0[i] = ta2_yz_0_xyy_0[i] * fe_0 - ta2_yz_0_xyy_1[i] * fe_0 + ta1_y_0_xyyz_1[i] + ta2_yz_0_xyyz_0[i] * pa_z[i] - ta2_yz_0_xyyz_1[i] * pc_z[i];

        ta2_yz_z_xyzz_0[i] = 2.0 * ta2_yz_0_xyz_0[i] * fe_0 - 2.0 * ta2_yz_0_xyz_1[i] * fe_0 + ta1_y_0_xyzz_1[i] + ta2_yz_0_xyzz_0[i] * pa_z[i] - ta2_yz_0_xyzz_1[i] * pc_z[i];

        ta2_yz_z_xzzz_0[i] = 3.0 * ta2_yz_0_xzz_0[i] * fe_0 - 3.0 * ta2_yz_0_xzz_1[i] * fe_0 + ta1_y_0_xzzz_1[i] + ta2_yz_0_xzzz_0[i] * pa_z[i] - ta2_yz_0_xzzz_1[i] * pc_z[i];

        ta2_yz_z_yyyy_0[i] = ta1_y_0_yyyy_1[i] + ta2_yz_0_yyyy_0[i] * pa_z[i] - ta2_yz_0_yyyy_1[i] * pc_z[i];

        ta2_yz_z_yyyz_0[i] = ta2_yz_0_yyy_0[i] * fe_0 - ta2_yz_0_yyy_1[i] * fe_0 + ta1_y_0_yyyz_1[i] + ta2_yz_0_yyyz_0[i] * pa_z[i] - ta2_yz_0_yyyz_1[i] * pc_z[i];

        ta2_yz_z_yyzz_0[i] = 2.0 * ta2_yz_0_yyz_0[i] * fe_0 - 2.0 * ta2_yz_0_yyz_1[i] * fe_0 + ta1_y_0_yyzz_1[i] + ta2_yz_0_yyzz_0[i] * pa_z[i] - ta2_yz_0_yyzz_1[i] * pc_z[i];

        ta2_yz_z_yzzz_0[i] = 3.0 * ta2_yz_0_yzz_0[i] * fe_0 - 3.0 * ta2_yz_0_yzz_1[i] * fe_0 + ta1_y_0_yzzz_1[i] + ta2_yz_0_yzzz_0[i] * pa_z[i] - ta2_yz_0_yzzz_1[i] * pc_z[i];

        ta2_yz_z_zzzz_0[i] = 4.0 * ta2_yz_0_zzz_0[i] * fe_0 - 4.0 * ta2_yz_0_zzz_1[i] * fe_0 + ta1_y_0_zzzz_1[i] + ta2_yz_0_zzzz_0[i] * pa_z[i] - ta2_yz_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 225-240 components of targeted buffer : PG

    auto ta2_zz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 225);

    auto ta2_zz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 226);

    auto ta2_zz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 227);

    auto ta2_zz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 228);

    auto ta2_zz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 229);

    auto ta2_zz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 230);

    auto ta2_zz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 231);

    auto ta2_zz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 232);

    auto ta2_zz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 233);

    auto ta2_zz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 234);

    auto ta2_zz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 235);

    auto ta2_zz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 236);

    auto ta2_zz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 237);

    auto ta2_zz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 238);

    auto ta2_zz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 239);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_0_xxx_0, ta2_zz_0_xxx_1, ta2_zz_0_xxxx_0, ta2_zz_0_xxxx_1, ta2_zz_0_xxxy_0, ta2_zz_0_xxxy_1, ta2_zz_0_xxxz_0, ta2_zz_0_xxxz_1, ta2_zz_0_xxy_0, ta2_zz_0_xxy_1, ta2_zz_0_xxyy_0, ta2_zz_0_xxyy_1, ta2_zz_0_xxyz_0, ta2_zz_0_xxyz_1, ta2_zz_0_xxz_0, ta2_zz_0_xxz_1, ta2_zz_0_xxzz_0, ta2_zz_0_xxzz_1, ta2_zz_0_xyy_0, ta2_zz_0_xyy_1, ta2_zz_0_xyyy_0, ta2_zz_0_xyyy_1, ta2_zz_0_xyyz_0, ta2_zz_0_xyyz_1, ta2_zz_0_xyz_0, ta2_zz_0_xyz_1, ta2_zz_0_xyzz_0, ta2_zz_0_xyzz_1, ta2_zz_0_xzz_0, ta2_zz_0_xzz_1, ta2_zz_0_xzzz_0, ta2_zz_0_xzzz_1, ta2_zz_0_yyy_0, ta2_zz_0_yyy_1, ta2_zz_0_yyyy_0, ta2_zz_0_yyyy_1, ta2_zz_0_yyyz_0, ta2_zz_0_yyyz_1, ta2_zz_0_yyz_0, ta2_zz_0_yyz_1, ta2_zz_0_yyzz_0, ta2_zz_0_yyzz_1, ta2_zz_0_yzz_0, ta2_zz_0_yzz_1, ta2_zz_0_yzzz_0, ta2_zz_0_yzzz_1, ta2_zz_0_zzz_0, ta2_zz_0_zzz_1, ta2_zz_0_zzzz_0, ta2_zz_0_zzzz_1, ta2_zz_x_xxxx_0, ta2_zz_x_xxxy_0, ta2_zz_x_xxxz_0, ta2_zz_x_xxyy_0, ta2_zz_x_xxyz_0, ta2_zz_x_xxzz_0, ta2_zz_x_xyyy_0, ta2_zz_x_xyyz_0, ta2_zz_x_xyzz_0, ta2_zz_x_xzzz_0, ta2_zz_x_yyyy_0, ta2_zz_x_yyyz_0, ta2_zz_x_yyzz_0, ta2_zz_x_yzzz_0, ta2_zz_x_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_x_xxxx_0[i] = 4.0 * ta2_zz_0_xxx_0[i] * fe_0 - 4.0 * ta2_zz_0_xxx_1[i] * fe_0 + ta2_zz_0_xxxx_0[i] * pa_x[i] - ta2_zz_0_xxxx_1[i] * pc_x[i];

        ta2_zz_x_xxxy_0[i] = 3.0 * ta2_zz_0_xxy_0[i] * fe_0 - 3.0 * ta2_zz_0_xxy_1[i] * fe_0 + ta2_zz_0_xxxy_0[i] * pa_x[i] - ta2_zz_0_xxxy_1[i] * pc_x[i];

        ta2_zz_x_xxxz_0[i] = 3.0 * ta2_zz_0_xxz_0[i] * fe_0 - 3.0 * ta2_zz_0_xxz_1[i] * fe_0 + ta2_zz_0_xxxz_0[i] * pa_x[i] - ta2_zz_0_xxxz_1[i] * pc_x[i];

        ta2_zz_x_xxyy_0[i] = 2.0 * ta2_zz_0_xyy_0[i] * fe_0 - 2.0 * ta2_zz_0_xyy_1[i] * fe_0 + ta2_zz_0_xxyy_0[i] * pa_x[i] - ta2_zz_0_xxyy_1[i] * pc_x[i];

        ta2_zz_x_xxyz_0[i] = 2.0 * ta2_zz_0_xyz_0[i] * fe_0 - 2.0 * ta2_zz_0_xyz_1[i] * fe_0 + ta2_zz_0_xxyz_0[i] * pa_x[i] - ta2_zz_0_xxyz_1[i] * pc_x[i];

        ta2_zz_x_xxzz_0[i] = 2.0 * ta2_zz_0_xzz_0[i] * fe_0 - 2.0 * ta2_zz_0_xzz_1[i] * fe_0 + ta2_zz_0_xxzz_0[i] * pa_x[i] - ta2_zz_0_xxzz_1[i] * pc_x[i];

        ta2_zz_x_xyyy_0[i] = ta2_zz_0_yyy_0[i] * fe_0 - ta2_zz_0_yyy_1[i] * fe_0 + ta2_zz_0_xyyy_0[i] * pa_x[i] - ta2_zz_0_xyyy_1[i] * pc_x[i];

        ta2_zz_x_xyyz_0[i] = ta2_zz_0_yyz_0[i] * fe_0 - ta2_zz_0_yyz_1[i] * fe_0 + ta2_zz_0_xyyz_0[i] * pa_x[i] - ta2_zz_0_xyyz_1[i] * pc_x[i];

        ta2_zz_x_xyzz_0[i] = ta2_zz_0_yzz_0[i] * fe_0 - ta2_zz_0_yzz_1[i] * fe_0 + ta2_zz_0_xyzz_0[i] * pa_x[i] - ta2_zz_0_xyzz_1[i] * pc_x[i];

        ta2_zz_x_xzzz_0[i] = ta2_zz_0_zzz_0[i] * fe_0 - ta2_zz_0_zzz_1[i] * fe_0 + ta2_zz_0_xzzz_0[i] * pa_x[i] - ta2_zz_0_xzzz_1[i] * pc_x[i];

        ta2_zz_x_yyyy_0[i] = ta2_zz_0_yyyy_0[i] * pa_x[i] - ta2_zz_0_yyyy_1[i] * pc_x[i];

        ta2_zz_x_yyyz_0[i] = ta2_zz_0_yyyz_0[i] * pa_x[i] - ta2_zz_0_yyyz_1[i] * pc_x[i];

        ta2_zz_x_yyzz_0[i] = ta2_zz_0_yyzz_0[i] * pa_x[i] - ta2_zz_0_yyzz_1[i] * pc_x[i];

        ta2_zz_x_yzzz_0[i] = ta2_zz_0_yzzz_0[i] * pa_x[i] - ta2_zz_0_yzzz_1[i] * pc_x[i];

        ta2_zz_x_zzzz_0[i] = ta2_zz_0_zzzz_0[i] * pa_x[i] - ta2_zz_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : PG

    auto ta2_zz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 240);

    auto ta2_zz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 241);

    auto ta2_zz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 242);

    auto ta2_zz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 243);

    auto ta2_zz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 244);

    auto ta2_zz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 245);

    auto ta2_zz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 246);

    auto ta2_zz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 247);

    auto ta2_zz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 248);

    auto ta2_zz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 249);

    auto ta2_zz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 250);

    auto ta2_zz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 251);

    auto ta2_zz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 252);

    auto ta2_zz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 253);

    auto ta2_zz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 254);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_0_xxx_0, ta2_zz_0_xxx_1, ta2_zz_0_xxxx_0, ta2_zz_0_xxxx_1, ta2_zz_0_xxxy_0, ta2_zz_0_xxxy_1, ta2_zz_0_xxxz_0, ta2_zz_0_xxxz_1, ta2_zz_0_xxy_0, ta2_zz_0_xxy_1, ta2_zz_0_xxyy_0, ta2_zz_0_xxyy_1, ta2_zz_0_xxyz_0, ta2_zz_0_xxyz_1, ta2_zz_0_xxz_0, ta2_zz_0_xxz_1, ta2_zz_0_xxzz_0, ta2_zz_0_xxzz_1, ta2_zz_0_xyy_0, ta2_zz_0_xyy_1, ta2_zz_0_xyyy_0, ta2_zz_0_xyyy_1, ta2_zz_0_xyyz_0, ta2_zz_0_xyyz_1, ta2_zz_0_xyz_0, ta2_zz_0_xyz_1, ta2_zz_0_xyzz_0, ta2_zz_0_xyzz_1, ta2_zz_0_xzz_0, ta2_zz_0_xzz_1, ta2_zz_0_xzzz_0, ta2_zz_0_xzzz_1, ta2_zz_0_yyy_0, ta2_zz_0_yyy_1, ta2_zz_0_yyyy_0, ta2_zz_0_yyyy_1, ta2_zz_0_yyyz_0, ta2_zz_0_yyyz_1, ta2_zz_0_yyz_0, ta2_zz_0_yyz_1, ta2_zz_0_yyzz_0, ta2_zz_0_yyzz_1, ta2_zz_0_yzz_0, ta2_zz_0_yzz_1, ta2_zz_0_yzzz_0, ta2_zz_0_yzzz_1, ta2_zz_0_zzz_0, ta2_zz_0_zzz_1, ta2_zz_0_zzzz_0, ta2_zz_0_zzzz_1, ta2_zz_y_xxxx_0, ta2_zz_y_xxxy_0, ta2_zz_y_xxxz_0, ta2_zz_y_xxyy_0, ta2_zz_y_xxyz_0, ta2_zz_y_xxzz_0, ta2_zz_y_xyyy_0, ta2_zz_y_xyyz_0, ta2_zz_y_xyzz_0, ta2_zz_y_xzzz_0, ta2_zz_y_yyyy_0, ta2_zz_y_yyyz_0, ta2_zz_y_yyzz_0, ta2_zz_y_yzzz_0, ta2_zz_y_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_y_xxxx_0[i] = ta2_zz_0_xxxx_0[i] * pa_y[i] - ta2_zz_0_xxxx_1[i] * pc_y[i];

        ta2_zz_y_xxxy_0[i] = ta2_zz_0_xxx_0[i] * fe_0 - ta2_zz_0_xxx_1[i] * fe_0 + ta2_zz_0_xxxy_0[i] * pa_y[i] - ta2_zz_0_xxxy_1[i] * pc_y[i];

        ta2_zz_y_xxxz_0[i] = ta2_zz_0_xxxz_0[i] * pa_y[i] - ta2_zz_0_xxxz_1[i] * pc_y[i];

        ta2_zz_y_xxyy_0[i] = 2.0 * ta2_zz_0_xxy_0[i] * fe_0 - 2.0 * ta2_zz_0_xxy_1[i] * fe_0 + ta2_zz_0_xxyy_0[i] * pa_y[i] - ta2_zz_0_xxyy_1[i] * pc_y[i];

        ta2_zz_y_xxyz_0[i] = ta2_zz_0_xxz_0[i] * fe_0 - ta2_zz_0_xxz_1[i] * fe_0 + ta2_zz_0_xxyz_0[i] * pa_y[i] - ta2_zz_0_xxyz_1[i] * pc_y[i];

        ta2_zz_y_xxzz_0[i] = ta2_zz_0_xxzz_0[i] * pa_y[i] - ta2_zz_0_xxzz_1[i] * pc_y[i];

        ta2_zz_y_xyyy_0[i] = 3.0 * ta2_zz_0_xyy_0[i] * fe_0 - 3.0 * ta2_zz_0_xyy_1[i] * fe_0 + ta2_zz_0_xyyy_0[i] * pa_y[i] - ta2_zz_0_xyyy_1[i] * pc_y[i];

        ta2_zz_y_xyyz_0[i] = 2.0 * ta2_zz_0_xyz_0[i] * fe_0 - 2.0 * ta2_zz_0_xyz_1[i] * fe_0 + ta2_zz_0_xyyz_0[i] * pa_y[i] - ta2_zz_0_xyyz_1[i] * pc_y[i];

        ta2_zz_y_xyzz_0[i] = ta2_zz_0_xzz_0[i] * fe_0 - ta2_zz_0_xzz_1[i] * fe_0 + ta2_zz_0_xyzz_0[i] * pa_y[i] - ta2_zz_0_xyzz_1[i] * pc_y[i];

        ta2_zz_y_xzzz_0[i] = ta2_zz_0_xzzz_0[i] * pa_y[i] - ta2_zz_0_xzzz_1[i] * pc_y[i];

        ta2_zz_y_yyyy_0[i] = 4.0 * ta2_zz_0_yyy_0[i] * fe_0 - 4.0 * ta2_zz_0_yyy_1[i] * fe_0 + ta2_zz_0_yyyy_0[i] * pa_y[i] - ta2_zz_0_yyyy_1[i] * pc_y[i];

        ta2_zz_y_yyyz_0[i] = 3.0 * ta2_zz_0_yyz_0[i] * fe_0 - 3.0 * ta2_zz_0_yyz_1[i] * fe_0 + ta2_zz_0_yyyz_0[i] * pa_y[i] - ta2_zz_0_yyyz_1[i] * pc_y[i];

        ta2_zz_y_yyzz_0[i] = 2.0 * ta2_zz_0_yzz_0[i] * fe_0 - 2.0 * ta2_zz_0_yzz_1[i] * fe_0 + ta2_zz_0_yyzz_0[i] * pa_y[i] - ta2_zz_0_yyzz_1[i] * pc_y[i];

        ta2_zz_y_yzzz_0[i] = ta2_zz_0_zzz_0[i] * fe_0 - ta2_zz_0_zzz_1[i] * fe_0 + ta2_zz_0_yzzz_0[i] * pa_y[i] - ta2_zz_0_yzzz_1[i] * pc_y[i];

        ta2_zz_y_zzzz_0[i] = ta2_zz_0_zzzz_0[i] * pa_y[i] - ta2_zz_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : PG

    auto ta2_zz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 255);

    auto ta2_zz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 256);

    auto ta2_zz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 257);

    auto ta2_zz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 258);

    auto ta2_zz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 259);

    auto ta2_zz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 260);

    auto ta2_zz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 261);

    auto ta2_zz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 262);

    auto ta2_zz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 263);

    auto ta2_zz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 264);

    auto ta2_zz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 265);

    auto ta2_zz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 266);

    auto ta2_zz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 267);

    auto ta2_zz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 268);

    auto ta2_zz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 269);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_0_xxxx_1, ta1_z_0_xxxy_1, ta1_z_0_xxxz_1, ta1_z_0_xxyy_1, ta1_z_0_xxyz_1, ta1_z_0_xxzz_1, ta1_z_0_xyyy_1, ta1_z_0_xyyz_1, ta1_z_0_xyzz_1, ta1_z_0_xzzz_1, ta1_z_0_yyyy_1, ta1_z_0_yyyz_1, ta1_z_0_yyzz_1, ta1_z_0_yzzz_1, ta1_z_0_zzzz_1, ta2_zz_0_xxx_0, ta2_zz_0_xxx_1, ta2_zz_0_xxxx_0, ta2_zz_0_xxxx_1, ta2_zz_0_xxxy_0, ta2_zz_0_xxxy_1, ta2_zz_0_xxxz_0, ta2_zz_0_xxxz_1, ta2_zz_0_xxy_0, ta2_zz_0_xxy_1, ta2_zz_0_xxyy_0, ta2_zz_0_xxyy_1, ta2_zz_0_xxyz_0, ta2_zz_0_xxyz_1, ta2_zz_0_xxz_0, ta2_zz_0_xxz_1, ta2_zz_0_xxzz_0, ta2_zz_0_xxzz_1, ta2_zz_0_xyy_0, ta2_zz_0_xyy_1, ta2_zz_0_xyyy_0, ta2_zz_0_xyyy_1, ta2_zz_0_xyyz_0, ta2_zz_0_xyyz_1, ta2_zz_0_xyz_0, ta2_zz_0_xyz_1, ta2_zz_0_xyzz_0, ta2_zz_0_xyzz_1, ta2_zz_0_xzz_0, ta2_zz_0_xzz_1, ta2_zz_0_xzzz_0, ta2_zz_0_xzzz_1, ta2_zz_0_yyy_0, ta2_zz_0_yyy_1, ta2_zz_0_yyyy_0, ta2_zz_0_yyyy_1, ta2_zz_0_yyyz_0, ta2_zz_0_yyyz_1, ta2_zz_0_yyz_0, ta2_zz_0_yyz_1, ta2_zz_0_yyzz_0, ta2_zz_0_yyzz_1, ta2_zz_0_yzz_0, ta2_zz_0_yzz_1, ta2_zz_0_yzzz_0, ta2_zz_0_yzzz_1, ta2_zz_0_zzz_0, ta2_zz_0_zzz_1, ta2_zz_0_zzzz_0, ta2_zz_0_zzzz_1, ta2_zz_z_xxxx_0, ta2_zz_z_xxxy_0, ta2_zz_z_xxxz_0, ta2_zz_z_xxyy_0, ta2_zz_z_xxyz_0, ta2_zz_z_xxzz_0, ta2_zz_z_xyyy_0, ta2_zz_z_xyyz_0, ta2_zz_z_xyzz_0, ta2_zz_z_xzzz_0, ta2_zz_z_yyyy_0, ta2_zz_z_yyyz_0, ta2_zz_z_yyzz_0, ta2_zz_z_yzzz_0, ta2_zz_z_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_z_xxxx_0[i] = 2.0 * ta1_z_0_xxxx_1[i] + ta2_zz_0_xxxx_0[i] * pa_z[i] - ta2_zz_0_xxxx_1[i] * pc_z[i];

        ta2_zz_z_xxxy_0[i] = 2.0 * ta1_z_0_xxxy_1[i] + ta2_zz_0_xxxy_0[i] * pa_z[i] - ta2_zz_0_xxxy_1[i] * pc_z[i];

        ta2_zz_z_xxxz_0[i] = ta2_zz_0_xxx_0[i] * fe_0 - ta2_zz_0_xxx_1[i] * fe_0 + 2.0 * ta1_z_0_xxxz_1[i] + ta2_zz_0_xxxz_0[i] * pa_z[i] - ta2_zz_0_xxxz_1[i] * pc_z[i];

        ta2_zz_z_xxyy_0[i] = 2.0 * ta1_z_0_xxyy_1[i] + ta2_zz_0_xxyy_0[i] * pa_z[i] - ta2_zz_0_xxyy_1[i] * pc_z[i];

        ta2_zz_z_xxyz_0[i] = ta2_zz_0_xxy_0[i] * fe_0 - ta2_zz_0_xxy_1[i] * fe_0 + 2.0 * ta1_z_0_xxyz_1[i] + ta2_zz_0_xxyz_0[i] * pa_z[i] - ta2_zz_0_xxyz_1[i] * pc_z[i];

        ta2_zz_z_xxzz_0[i] = 2.0 * ta2_zz_0_xxz_0[i] * fe_0 - 2.0 * ta2_zz_0_xxz_1[i] * fe_0 + 2.0 * ta1_z_0_xxzz_1[i] + ta2_zz_0_xxzz_0[i] * pa_z[i] - ta2_zz_0_xxzz_1[i] * pc_z[i];

        ta2_zz_z_xyyy_0[i] = 2.0 * ta1_z_0_xyyy_1[i] + ta2_zz_0_xyyy_0[i] * pa_z[i] - ta2_zz_0_xyyy_1[i] * pc_z[i];

        ta2_zz_z_xyyz_0[i] = ta2_zz_0_xyy_0[i] * fe_0 - ta2_zz_0_xyy_1[i] * fe_0 + 2.0 * ta1_z_0_xyyz_1[i] + ta2_zz_0_xyyz_0[i] * pa_z[i] - ta2_zz_0_xyyz_1[i] * pc_z[i];

        ta2_zz_z_xyzz_0[i] = 2.0 * ta2_zz_0_xyz_0[i] * fe_0 - 2.0 * ta2_zz_0_xyz_1[i] * fe_0 + 2.0 * ta1_z_0_xyzz_1[i] + ta2_zz_0_xyzz_0[i] * pa_z[i] - ta2_zz_0_xyzz_1[i] * pc_z[i];

        ta2_zz_z_xzzz_0[i] = 3.0 * ta2_zz_0_xzz_0[i] * fe_0 - 3.0 * ta2_zz_0_xzz_1[i] * fe_0 + 2.0 * ta1_z_0_xzzz_1[i] + ta2_zz_0_xzzz_0[i] * pa_z[i] - ta2_zz_0_xzzz_1[i] * pc_z[i];

        ta2_zz_z_yyyy_0[i] = 2.0 * ta1_z_0_yyyy_1[i] + ta2_zz_0_yyyy_0[i] * pa_z[i] - ta2_zz_0_yyyy_1[i] * pc_z[i];

        ta2_zz_z_yyyz_0[i] = ta2_zz_0_yyy_0[i] * fe_0 - ta2_zz_0_yyy_1[i] * fe_0 + 2.0 * ta1_z_0_yyyz_1[i] + ta2_zz_0_yyyz_0[i] * pa_z[i] - ta2_zz_0_yyyz_1[i] * pc_z[i];

        ta2_zz_z_yyzz_0[i] = 2.0 * ta2_zz_0_yyz_0[i] * fe_0 - 2.0 * ta2_zz_0_yyz_1[i] * fe_0 + 2.0 * ta1_z_0_yyzz_1[i] + ta2_zz_0_yyzz_0[i] * pa_z[i] - ta2_zz_0_yyzz_1[i] * pc_z[i];

        ta2_zz_z_yzzz_0[i] = 3.0 * ta2_zz_0_yzz_0[i] * fe_0 - 3.0 * ta2_zz_0_yzz_1[i] * fe_0 + 2.0 * ta1_z_0_yzzz_1[i] + ta2_zz_0_yzzz_0[i] * pa_z[i] - ta2_zz_0_yzzz_1[i] * pc_z[i];

        ta2_zz_z_zzzz_0[i] = 4.0 * ta2_zz_0_zzz_0[i] * fe_0 - 4.0 * ta2_zz_0_zzz_1[i] * fe_0 + 2.0 * ta1_z_0_zzzz_1[i] + ta2_zz_0_zzzz_0[i] * pa_z[i] - ta2_zz_0_zzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

