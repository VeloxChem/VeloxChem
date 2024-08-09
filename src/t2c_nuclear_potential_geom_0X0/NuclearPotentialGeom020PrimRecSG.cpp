#include "NuclearPotentialGeom020PrimRecSG.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_sg(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_sg,
                                        const size_t idx_npot_geom_020_0_sd,
                                        const size_t idx_npot_geom_020_1_sd,
                                        const size_t idx_npot_geom_010_1_sf,
                                        const size_t idx_npot_geom_020_0_sf,
                                        const size_t idx_npot_geom_020_1_sf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpb,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd);

    auto ta2_xx_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 3);

    auto ta2_xx_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 5);

    auto ta2_xy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 6);

    auto ta2_xy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 9);

    auto ta2_xy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 11);

    auto ta2_xz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 12);

    auto ta2_xz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 15);

    auto ta2_xz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 17);

    auto ta2_yy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 18);

    auto ta2_yy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 21);

    auto ta2_yy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 23);

    auto ta2_yz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 24);

    auto ta2_yz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 27);

    auto ta2_yz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 29);

    auto ta2_zz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 30);

    auto ta2_zz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 33);

    auto ta2_zz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 35);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd);

    auto ta2_xx_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 3);

    auto ta2_xx_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 5);

    auto ta2_xy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 6);

    auto ta2_xy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 9);

    auto ta2_xy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 11);

    auto ta2_xz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 12);

    auto ta2_xz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 15);

    auto ta2_xz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 17);

    auto ta2_yy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 18);

    auto ta2_yy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 21);

    auto ta2_yy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 23);

    auto ta2_yz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 24);

    auto ta2_yz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 27);

    auto ta2_yz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 29);

    auto ta2_zz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 30);

    auto ta2_zz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 33);

    auto ta2_zz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 35);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

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

    auto ta2_xx_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 5);

    auto ta2_xx_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 6);

    auto ta2_xx_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 8);

    auto ta2_xx_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 9);

    auto ta2_xy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 10);

    auto ta2_xy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 11);

    auto ta2_xy_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 12);

    auto ta2_xy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 13);

    auto ta2_xy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 16);

    auto ta2_xy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 17);

    auto ta2_xy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 18);

    auto ta2_xy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 19);

    auto ta2_xz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 20);

    auto ta2_xz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 21);

    auto ta2_xz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 22);

    auto ta2_xz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 25);

    auto ta2_xz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 26);

    auto ta2_xz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 27);

    auto ta2_xz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 28);

    auto ta2_xz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 29);

    auto ta2_yy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 30);

    auto ta2_yy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 31);

    auto ta2_yy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 33);

    auto ta2_yy_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 35);

    auto ta2_yy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 36);

    auto ta2_yy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 37);

    auto ta2_yy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 38);

    auto ta2_yy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 39);

    auto ta2_yz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 40);

    auto ta2_yz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 42);

    auto ta2_yz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 43);

    auto ta2_yz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 45);

    auto ta2_yz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 46);

    auto ta2_yz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 47);

    auto ta2_yz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 48);

    auto ta2_yz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 49);

    auto ta2_zz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 50);

    auto ta2_zz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 52);

    auto ta2_zz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 53);

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

    auto ta2_xx_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 5);

    auto ta2_xx_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 6);

    auto ta2_xx_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 8);

    auto ta2_xx_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 9);

    auto ta2_xy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 10);

    auto ta2_xy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 11);

    auto ta2_xy_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 12);

    auto ta2_xy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 13);

    auto ta2_xy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 16);

    auto ta2_xy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 17);

    auto ta2_xy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 18);

    auto ta2_xy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 19);

    auto ta2_xz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 20);

    auto ta2_xz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 21);

    auto ta2_xz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 22);

    auto ta2_xz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 25);

    auto ta2_xz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 26);

    auto ta2_xz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 27);

    auto ta2_xz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 28);

    auto ta2_xz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 29);

    auto ta2_yy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 30);

    auto ta2_yy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 31);

    auto ta2_yy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 33);

    auto ta2_yy_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 35);

    auto ta2_yy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 36);

    auto ta2_yy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 37);

    auto ta2_yy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 38);

    auto ta2_yy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 39);

    auto ta2_yz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 40);

    auto ta2_yz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 42);

    auto ta2_yz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 43);

    auto ta2_yz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 45);

    auto ta2_yz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 46);

    auto ta2_yz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 47);

    auto ta2_yz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 48);

    auto ta2_yz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 49);

    auto ta2_zz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 50);

    auto ta2_zz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 52);

    auto ta2_zz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 53);

    auto ta2_zz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 55);

    auto ta2_zz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 56);

    auto ta2_zz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 57);

    auto ta2_zz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 58);

    auto ta2_zz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 59);

    // Set up components of targeted buffer : SG

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

    #pragma omp simd aligned(pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta1_x_0_xxx_1, ta1_x_0_yyy_1, ta1_x_0_zzz_1, ta1_y_0_xxx_1, ta1_y_0_xyy_1, ta1_y_0_yyy_1, ta1_y_0_yzz_1, ta1_y_0_zzz_1, ta1_z_0_xxx_1, ta1_z_0_xxz_1, ta1_z_0_xzz_1, ta1_z_0_yyy_1, ta1_z_0_yyz_1, ta1_z_0_yzz_1, ta1_z_0_zzz_1, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xxx_0, ta2_xx_0_xxx_1, ta2_xx_0_xxxx_0, ta2_xx_0_xxxy_0, ta2_xx_0_xxxz_0, ta2_xx_0_xxy_0, ta2_xx_0_xxy_1, ta2_xx_0_xxyy_0, ta2_xx_0_xxyz_0, ta2_xx_0_xxz_0, ta2_xx_0_xxz_1, ta2_xx_0_xxzz_0, ta2_xx_0_xyy_0, ta2_xx_0_xyy_1, ta2_xx_0_xyyy_0, ta2_xx_0_xyyz_0, ta2_xx_0_xyzz_0, ta2_xx_0_xzz_0, ta2_xx_0_xzz_1, ta2_xx_0_xzzz_0, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yyy_0, ta2_xx_0_yyy_1, ta2_xx_0_yyyy_0, ta2_xx_0_yyyz_0, ta2_xx_0_yyzz_0, ta2_xx_0_yzz_0, ta2_xx_0_yzz_1, ta2_xx_0_yzzz_0, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_0_zzz_0, ta2_xx_0_zzz_1, ta2_xx_0_zzzz_0, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xxx_0, ta2_xy_0_xxx_1, ta2_xy_0_xxxx_0, ta2_xy_0_xxxy_0, ta2_xy_0_xxxz_0, ta2_xy_0_xxy_0, ta2_xy_0_xxy_1, ta2_xy_0_xxyy_0, ta2_xy_0_xxyz_0, ta2_xy_0_xxz_0, ta2_xy_0_xxz_1, ta2_xy_0_xxzz_0, ta2_xy_0_xyy_0, ta2_xy_0_xyy_1, ta2_xy_0_xyyy_0, ta2_xy_0_xyyz_0, ta2_xy_0_xyzz_0, ta2_xy_0_xzzz_0, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yyy_0, ta2_xy_0_yyy_1, ta2_xy_0_yyyy_0, ta2_xy_0_yyyz_0, ta2_xy_0_yyz_0, ta2_xy_0_yyz_1, ta2_xy_0_yyzz_0, ta2_xy_0_yzz_0, ta2_xy_0_yzz_1, ta2_xy_0_yzzz_0, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_0_zzz_0, ta2_xy_0_zzz_1, ta2_xy_0_zzzz_0, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xxx_0, ta2_xz_0_xxx_1, ta2_xz_0_xxxx_0, ta2_xz_0_xxxy_0, ta2_xz_0_xxxz_0, ta2_xz_0_xxy_0, ta2_xz_0_xxy_1, ta2_xz_0_xxyy_0, ta2_xz_0_xxyz_0, ta2_xz_0_xxz_0, ta2_xz_0_xxz_1, ta2_xz_0_xxzz_0, ta2_xz_0_xyyy_0, ta2_xz_0_xyyz_0, ta2_xz_0_xyzz_0, ta2_xz_0_xzz_0, ta2_xz_0_xzz_1, ta2_xz_0_xzzz_0, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yyy_0, ta2_xz_0_yyy_1, ta2_xz_0_yyyy_0, ta2_xz_0_yyyz_0, ta2_xz_0_yyz_0, ta2_xz_0_yyz_1, ta2_xz_0_yyzz_0, ta2_xz_0_yzz_0, ta2_xz_0_yzz_1, ta2_xz_0_yzzz_0, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_0_zzz_0, ta2_xz_0_zzz_1, ta2_xz_0_zzzz_0, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xxx_0, ta2_yy_0_xxx_1, ta2_yy_0_xxxx_0, ta2_yy_0_xxxy_0, ta2_yy_0_xxxz_0, ta2_yy_0_xxy_0, ta2_yy_0_xxy_1, ta2_yy_0_xxyy_0, ta2_yy_0_xxyz_0, ta2_yy_0_xxzz_0, ta2_yy_0_xyy_0, ta2_yy_0_xyy_1, ta2_yy_0_xyyy_0, ta2_yy_0_xyyz_0, ta2_yy_0_xyzz_0, ta2_yy_0_xzz_0, ta2_yy_0_xzz_1, ta2_yy_0_xzzz_0, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yyy_0, ta2_yy_0_yyy_1, ta2_yy_0_yyyy_0, ta2_yy_0_yyyz_0, ta2_yy_0_yyz_0, ta2_yy_0_yyz_1, ta2_yy_0_yyzz_0, ta2_yy_0_yzz_0, ta2_yy_0_yzz_1, ta2_yy_0_yzzz_0, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_0_zzz_0, ta2_yy_0_zzz_1, ta2_yy_0_zzzz_0, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xxx_0, ta2_yz_0_xxx_1, ta2_yz_0_xxxx_0, ta2_yz_0_xxxy_0, ta2_yz_0_xxxz_0, ta2_yz_0_xxyy_0, ta2_yz_0_xxyz_0, ta2_yz_0_xxz_0, ta2_yz_0_xxz_1, ta2_yz_0_xxzz_0, ta2_yz_0_xyy_0, ta2_yz_0_xyy_1, ta2_yz_0_xyyy_0, ta2_yz_0_xyyz_0, ta2_yz_0_xyzz_0, ta2_yz_0_xzz_0, ta2_yz_0_xzz_1, ta2_yz_0_xzzz_0, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yyy_0, ta2_yz_0_yyy_1, ta2_yz_0_yyyy_0, ta2_yz_0_yyyz_0, ta2_yz_0_yyz_0, ta2_yz_0_yyz_1, ta2_yz_0_yyzz_0, ta2_yz_0_yzz_0, ta2_yz_0_yzz_1, ta2_yz_0_yzzz_0, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_0_zzz_0, ta2_yz_0_zzz_1, ta2_yz_0_zzzz_0, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xxx_0, ta2_zz_0_xxx_1, ta2_zz_0_xxxx_0, ta2_zz_0_xxxy_0, ta2_zz_0_xxxz_0, ta2_zz_0_xxyy_0, ta2_zz_0_xxyz_0, ta2_zz_0_xxz_0, ta2_zz_0_xxz_1, ta2_zz_0_xxzz_0, ta2_zz_0_xyy_0, ta2_zz_0_xyy_1, ta2_zz_0_xyyy_0, ta2_zz_0_xyyz_0, ta2_zz_0_xyzz_0, ta2_zz_0_xzz_0, ta2_zz_0_xzz_1, ta2_zz_0_xzzz_0, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yyy_0, ta2_zz_0_yyy_1, ta2_zz_0_yyyy_0, ta2_zz_0_yyyz_0, ta2_zz_0_yyz_0, ta2_zz_0_yyz_1, ta2_zz_0_yyzz_0, ta2_zz_0_yzz_0, ta2_zz_0_yzz_1, ta2_zz_0_yzzz_0, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_0_zzz_0, ta2_zz_0_zzz_1, ta2_zz_0_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_0_xxxx_0[i] = 3.0 * ta2_xx_0_xx_0[i] * fe_0 - 3.0 * ta2_xx_0_xx_1[i] * fe_0 + 2.0 * ta1_x_0_xxx_1[i] + ta2_xx_0_xxx_0[i] * pb_x[i] - ta2_xx_0_xxx_1[i] * pc_x[i];

        ta2_xx_0_xxxy_0[i] = ta2_xx_0_xxx_0[i] * pb_y[i] - ta2_xx_0_xxx_1[i] * pc_y[i];

        ta2_xx_0_xxxz_0[i] = ta2_xx_0_xxx_0[i] * pb_z[i] - ta2_xx_0_xxx_1[i] * pc_z[i];

        ta2_xx_0_xxyy_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_0_xxy_0[i] * pb_y[i] - ta2_xx_0_xxy_1[i] * pc_y[i];

        ta2_xx_0_xxyz_0[i] = ta2_xx_0_xxz_0[i] * pb_y[i] - ta2_xx_0_xxz_1[i] * pc_y[i];

        ta2_xx_0_xxzz_0[i] = ta2_xx_0_xx_0[i] * fe_0 - ta2_xx_0_xx_1[i] * fe_0 + ta2_xx_0_xxz_0[i] * pb_z[i] - ta2_xx_0_xxz_1[i] * pc_z[i];

        ta2_xx_0_xyyy_0[i] = 2.0 * ta1_x_0_yyy_1[i] + ta2_xx_0_yyy_0[i] * pb_x[i] - ta2_xx_0_yyy_1[i] * pc_x[i];

        ta2_xx_0_xyyz_0[i] = ta2_xx_0_xyy_0[i] * pb_z[i] - ta2_xx_0_xyy_1[i] * pc_z[i];

        ta2_xx_0_xyzz_0[i] = ta2_xx_0_xzz_0[i] * pb_y[i] - ta2_xx_0_xzz_1[i] * pc_y[i];

        ta2_xx_0_xzzz_0[i] = 2.0 * ta1_x_0_zzz_1[i] + ta2_xx_0_zzz_0[i] * pb_x[i] - ta2_xx_0_zzz_1[i] * pc_x[i];

        ta2_xx_0_yyyy_0[i] = 3.0 * ta2_xx_0_yy_0[i] * fe_0 - 3.0 * ta2_xx_0_yy_1[i] * fe_0 + ta2_xx_0_yyy_0[i] * pb_y[i] - ta2_xx_0_yyy_1[i] * pc_y[i];

        ta2_xx_0_yyyz_0[i] = ta2_xx_0_yyy_0[i] * pb_z[i] - ta2_xx_0_yyy_1[i] * pc_z[i];

        ta2_xx_0_yyzz_0[i] = ta2_xx_0_zz_0[i] * fe_0 - ta2_xx_0_zz_1[i] * fe_0 + ta2_xx_0_yzz_0[i] * pb_y[i] - ta2_xx_0_yzz_1[i] * pc_y[i];

        ta2_xx_0_yzzz_0[i] = ta2_xx_0_zzz_0[i] * pb_y[i] - ta2_xx_0_zzz_1[i] * pc_y[i];

        ta2_xx_0_zzzz_0[i] = 3.0 * ta2_xx_0_zz_0[i] * fe_0 - 3.0 * ta2_xx_0_zz_1[i] * fe_0 + ta2_xx_0_zzz_0[i] * pb_z[i] - ta2_xx_0_zzz_1[i] * pc_z[i];

        ta2_xy_0_xxxx_0[i] = 3.0 * ta2_xy_0_xx_0[i] * fe_0 - 3.0 * ta2_xy_0_xx_1[i] * fe_0 + ta1_y_0_xxx_1[i] + ta2_xy_0_xxx_0[i] * pb_x[i] - ta2_xy_0_xxx_1[i] * pc_x[i];

        ta2_xy_0_xxxy_0[i] = ta1_x_0_xxx_1[i] + ta2_xy_0_xxx_0[i] * pb_y[i] - ta2_xy_0_xxx_1[i] * pc_y[i];

        ta2_xy_0_xxxz_0[i] = ta2_xy_0_xxx_0[i] * pb_z[i] - ta2_xy_0_xxx_1[i] * pc_z[i];

        ta2_xy_0_xxyy_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta1_y_0_xyy_1[i] + ta2_xy_0_xyy_0[i] * pb_x[i] - ta2_xy_0_xyy_1[i] * pc_x[i];

        ta2_xy_0_xxyz_0[i] = ta2_xy_0_xxy_0[i] * pb_z[i] - ta2_xy_0_xxy_1[i] * pc_z[i];

        ta2_xy_0_xxzz_0[i] = ta2_xy_0_xx_0[i] * fe_0 - ta2_xy_0_xx_1[i] * fe_0 + ta2_xy_0_xxz_0[i] * pb_z[i] - ta2_xy_0_xxz_1[i] * pc_z[i];

        ta2_xy_0_xyyy_0[i] = ta1_y_0_yyy_1[i] + ta2_xy_0_yyy_0[i] * pb_x[i] - ta2_xy_0_yyy_1[i] * pc_x[i];

        ta2_xy_0_xyyz_0[i] = ta2_xy_0_xyy_0[i] * pb_z[i] - ta2_xy_0_xyy_1[i] * pc_z[i];

        ta2_xy_0_xyzz_0[i] = ta1_y_0_yzz_1[i] + ta2_xy_0_yzz_0[i] * pb_x[i] - ta2_xy_0_yzz_1[i] * pc_x[i];

        ta2_xy_0_xzzz_0[i] = ta1_y_0_zzz_1[i] + ta2_xy_0_zzz_0[i] * pb_x[i] - ta2_xy_0_zzz_1[i] * pc_x[i];

        ta2_xy_0_yyyy_0[i] = 3.0 * ta2_xy_0_yy_0[i] * fe_0 - 3.0 * ta2_xy_0_yy_1[i] * fe_0 + ta1_x_0_yyy_1[i] + ta2_xy_0_yyy_0[i] * pb_y[i] - ta2_xy_0_yyy_1[i] * pc_y[i];

        ta2_xy_0_yyyz_0[i] = ta2_xy_0_yyy_0[i] * pb_z[i] - ta2_xy_0_yyy_1[i] * pc_z[i];

        ta2_xy_0_yyzz_0[i] = ta2_xy_0_yy_0[i] * fe_0 - ta2_xy_0_yy_1[i] * fe_0 + ta2_xy_0_yyz_0[i] * pb_z[i] - ta2_xy_0_yyz_1[i] * pc_z[i];

        ta2_xy_0_yzzz_0[i] = ta1_x_0_zzz_1[i] + ta2_xy_0_zzz_0[i] * pb_y[i] - ta2_xy_0_zzz_1[i] * pc_y[i];

        ta2_xy_0_zzzz_0[i] = 3.0 * ta2_xy_0_zz_0[i] * fe_0 - 3.0 * ta2_xy_0_zz_1[i] * fe_0 + ta2_xy_0_zzz_0[i] * pb_z[i] - ta2_xy_0_zzz_1[i] * pc_z[i];

        ta2_xz_0_xxxx_0[i] = 3.0 * ta2_xz_0_xx_0[i] * fe_0 - 3.0 * ta2_xz_0_xx_1[i] * fe_0 + ta1_z_0_xxx_1[i] + ta2_xz_0_xxx_0[i] * pb_x[i] - ta2_xz_0_xxx_1[i] * pc_x[i];

        ta2_xz_0_xxxy_0[i] = ta2_xz_0_xxx_0[i] * pb_y[i] - ta2_xz_0_xxx_1[i] * pc_y[i];

        ta2_xz_0_xxxz_0[i] = ta1_x_0_xxx_1[i] + ta2_xz_0_xxx_0[i] * pb_z[i] - ta2_xz_0_xxx_1[i] * pc_z[i];

        ta2_xz_0_xxyy_0[i] = ta2_xz_0_xx_0[i] * fe_0 - ta2_xz_0_xx_1[i] * fe_0 + ta2_xz_0_xxy_0[i] * pb_y[i] - ta2_xz_0_xxy_1[i] * pc_y[i];

        ta2_xz_0_xxyz_0[i] = ta2_xz_0_xxz_0[i] * pb_y[i] - ta2_xz_0_xxz_1[i] * pc_y[i];

        ta2_xz_0_xxzz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta1_z_0_xzz_1[i] + ta2_xz_0_xzz_0[i] * pb_x[i] - ta2_xz_0_xzz_1[i] * pc_x[i];

        ta2_xz_0_xyyy_0[i] = ta1_z_0_yyy_1[i] + ta2_xz_0_yyy_0[i] * pb_x[i] - ta2_xz_0_yyy_1[i] * pc_x[i];

        ta2_xz_0_xyyz_0[i] = ta1_z_0_yyz_1[i] + ta2_xz_0_yyz_0[i] * pb_x[i] - ta2_xz_0_yyz_1[i] * pc_x[i];

        ta2_xz_0_xyzz_0[i] = ta2_xz_0_xzz_0[i] * pb_y[i] - ta2_xz_0_xzz_1[i] * pc_y[i];

        ta2_xz_0_xzzz_0[i] = ta1_z_0_zzz_1[i] + ta2_xz_0_zzz_0[i] * pb_x[i] - ta2_xz_0_zzz_1[i] * pc_x[i];

        ta2_xz_0_yyyy_0[i] = 3.0 * ta2_xz_0_yy_0[i] * fe_0 - 3.0 * ta2_xz_0_yy_1[i] * fe_0 + ta2_xz_0_yyy_0[i] * pb_y[i] - ta2_xz_0_yyy_1[i] * pc_y[i];

        ta2_xz_0_yyyz_0[i] = ta1_x_0_yyy_1[i] + ta2_xz_0_yyy_0[i] * pb_z[i] - ta2_xz_0_yyy_1[i] * pc_z[i];

        ta2_xz_0_yyzz_0[i] = ta2_xz_0_zz_0[i] * fe_0 - ta2_xz_0_zz_1[i] * fe_0 + ta2_xz_0_yzz_0[i] * pb_y[i] - ta2_xz_0_yzz_1[i] * pc_y[i];

        ta2_xz_0_yzzz_0[i] = ta2_xz_0_zzz_0[i] * pb_y[i] - ta2_xz_0_zzz_1[i] * pc_y[i];

        ta2_xz_0_zzzz_0[i] = 3.0 * ta2_xz_0_zz_0[i] * fe_0 - 3.0 * ta2_xz_0_zz_1[i] * fe_0 + ta1_x_0_zzz_1[i] + ta2_xz_0_zzz_0[i] * pb_z[i] - ta2_xz_0_zzz_1[i] * pc_z[i];

        ta2_yy_0_xxxx_0[i] = 3.0 * ta2_yy_0_xx_0[i] * fe_0 - 3.0 * ta2_yy_0_xx_1[i] * fe_0 + ta2_yy_0_xxx_0[i] * pb_x[i] - ta2_yy_0_xxx_1[i] * pc_x[i];

        ta2_yy_0_xxxy_0[i] = 2.0 * ta1_y_0_xxx_1[i] + ta2_yy_0_xxx_0[i] * pb_y[i] - ta2_yy_0_xxx_1[i] * pc_y[i];

        ta2_yy_0_xxxz_0[i] = ta2_yy_0_xxx_0[i] * pb_z[i] - ta2_yy_0_xxx_1[i] * pc_z[i];

        ta2_yy_0_xxyy_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_0_xyy_0[i] * pb_x[i] - ta2_yy_0_xyy_1[i] * pc_x[i];

        ta2_yy_0_xxyz_0[i] = ta2_yy_0_xxy_0[i] * pb_z[i] - ta2_yy_0_xxy_1[i] * pc_z[i];

        ta2_yy_0_xxzz_0[i] = ta2_yy_0_zz_0[i] * fe_0 - ta2_yy_0_zz_1[i] * fe_0 + ta2_yy_0_xzz_0[i] * pb_x[i] - ta2_yy_0_xzz_1[i] * pc_x[i];

        ta2_yy_0_xyyy_0[i] = ta2_yy_0_yyy_0[i] * pb_x[i] - ta2_yy_0_yyy_1[i] * pc_x[i];

        ta2_yy_0_xyyz_0[i] = ta2_yy_0_yyz_0[i] * pb_x[i] - ta2_yy_0_yyz_1[i] * pc_x[i];

        ta2_yy_0_xyzz_0[i] = ta2_yy_0_yzz_0[i] * pb_x[i] - ta2_yy_0_yzz_1[i] * pc_x[i];

        ta2_yy_0_xzzz_0[i] = ta2_yy_0_zzz_0[i] * pb_x[i] - ta2_yy_0_zzz_1[i] * pc_x[i];

        ta2_yy_0_yyyy_0[i] = 3.0 * ta2_yy_0_yy_0[i] * fe_0 - 3.0 * ta2_yy_0_yy_1[i] * fe_0 + 2.0 * ta1_y_0_yyy_1[i] + ta2_yy_0_yyy_0[i] * pb_y[i] - ta2_yy_0_yyy_1[i] * pc_y[i];

        ta2_yy_0_yyyz_0[i] = ta2_yy_0_yyy_0[i] * pb_z[i] - ta2_yy_0_yyy_1[i] * pc_z[i];

        ta2_yy_0_yyzz_0[i] = ta2_yy_0_yy_0[i] * fe_0 - ta2_yy_0_yy_1[i] * fe_0 + ta2_yy_0_yyz_0[i] * pb_z[i] - ta2_yy_0_yyz_1[i] * pc_z[i];

        ta2_yy_0_yzzz_0[i] = 2.0 * ta1_y_0_zzz_1[i] + ta2_yy_0_zzz_0[i] * pb_y[i] - ta2_yy_0_zzz_1[i] * pc_y[i];

        ta2_yy_0_zzzz_0[i] = 3.0 * ta2_yy_0_zz_0[i] * fe_0 - 3.0 * ta2_yy_0_zz_1[i] * fe_0 + ta2_yy_0_zzz_0[i] * pb_z[i] - ta2_yy_0_zzz_1[i] * pc_z[i];

        ta2_yz_0_xxxx_0[i] = 3.0 * ta2_yz_0_xx_0[i] * fe_0 - 3.0 * ta2_yz_0_xx_1[i] * fe_0 + ta2_yz_0_xxx_0[i] * pb_x[i] - ta2_yz_0_xxx_1[i] * pc_x[i];

        ta2_yz_0_xxxy_0[i] = ta1_z_0_xxx_1[i] + ta2_yz_0_xxx_0[i] * pb_y[i] - ta2_yz_0_xxx_1[i] * pc_y[i];

        ta2_yz_0_xxxz_0[i] = ta1_y_0_xxx_1[i] + ta2_yz_0_xxx_0[i] * pb_z[i] - ta2_yz_0_xxx_1[i] * pc_z[i];

        ta2_yz_0_xxyy_0[i] = ta2_yz_0_yy_0[i] * fe_0 - ta2_yz_0_yy_1[i] * fe_0 + ta2_yz_0_xyy_0[i] * pb_x[i] - ta2_yz_0_xyy_1[i] * pc_x[i];

        ta2_yz_0_xxyz_0[i] = ta1_z_0_xxz_1[i] + ta2_yz_0_xxz_0[i] * pb_y[i] - ta2_yz_0_xxz_1[i] * pc_y[i];

        ta2_yz_0_xxzz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta2_yz_0_xzz_0[i] * pb_x[i] - ta2_yz_0_xzz_1[i] * pc_x[i];

        ta2_yz_0_xyyy_0[i] = ta2_yz_0_yyy_0[i] * pb_x[i] - ta2_yz_0_yyy_1[i] * pc_x[i];

        ta2_yz_0_xyyz_0[i] = ta2_yz_0_yyz_0[i] * pb_x[i] - ta2_yz_0_yyz_1[i] * pc_x[i];

        ta2_yz_0_xyzz_0[i] = ta2_yz_0_yzz_0[i] * pb_x[i] - ta2_yz_0_yzz_1[i] * pc_x[i];

        ta2_yz_0_xzzz_0[i] = ta2_yz_0_zzz_0[i] * pb_x[i] - ta2_yz_0_zzz_1[i] * pc_x[i];

        ta2_yz_0_yyyy_0[i] = 3.0 * ta2_yz_0_yy_0[i] * fe_0 - 3.0 * ta2_yz_0_yy_1[i] * fe_0 + ta1_z_0_yyy_1[i] + ta2_yz_0_yyy_0[i] * pb_y[i] - ta2_yz_0_yyy_1[i] * pc_y[i];

        ta2_yz_0_yyyz_0[i] = ta1_y_0_yyy_1[i] + ta2_yz_0_yyy_0[i] * pb_z[i] - ta2_yz_0_yyy_1[i] * pc_z[i];

        ta2_yz_0_yyzz_0[i] = ta2_yz_0_zz_0[i] * fe_0 - ta2_yz_0_zz_1[i] * fe_0 + ta1_z_0_yzz_1[i] + ta2_yz_0_yzz_0[i] * pb_y[i] - ta2_yz_0_yzz_1[i] * pc_y[i];

        ta2_yz_0_yzzz_0[i] = ta1_z_0_zzz_1[i] + ta2_yz_0_zzz_0[i] * pb_y[i] - ta2_yz_0_zzz_1[i] * pc_y[i];

        ta2_yz_0_zzzz_0[i] = 3.0 * ta2_yz_0_zz_0[i] * fe_0 - 3.0 * ta2_yz_0_zz_1[i] * fe_0 + ta1_y_0_zzz_1[i] + ta2_yz_0_zzz_0[i] * pb_z[i] - ta2_yz_0_zzz_1[i] * pc_z[i];

        ta2_zz_0_xxxx_0[i] = 3.0 * ta2_zz_0_xx_0[i] * fe_0 - 3.0 * ta2_zz_0_xx_1[i] * fe_0 + ta2_zz_0_xxx_0[i] * pb_x[i] - ta2_zz_0_xxx_1[i] * pc_x[i];

        ta2_zz_0_xxxy_0[i] = ta2_zz_0_xxx_0[i] * pb_y[i] - ta2_zz_0_xxx_1[i] * pc_y[i];

        ta2_zz_0_xxxz_0[i] = 2.0 * ta1_z_0_xxx_1[i] + ta2_zz_0_xxx_0[i] * pb_z[i] - ta2_zz_0_xxx_1[i] * pc_z[i];

        ta2_zz_0_xxyy_0[i] = ta2_zz_0_yy_0[i] * fe_0 - ta2_zz_0_yy_1[i] * fe_0 + ta2_zz_0_xyy_0[i] * pb_x[i] - ta2_zz_0_xyy_1[i] * pc_x[i];

        ta2_zz_0_xxyz_0[i] = ta2_zz_0_xxz_0[i] * pb_y[i] - ta2_zz_0_xxz_1[i] * pc_y[i];

        ta2_zz_0_xxzz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_0_xzz_0[i] * pb_x[i] - ta2_zz_0_xzz_1[i] * pc_x[i];

        ta2_zz_0_xyyy_0[i] = ta2_zz_0_yyy_0[i] * pb_x[i] - ta2_zz_0_yyy_1[i] * pc_x[i];

        ta2_zz_0_xyyz_0[i] = ta2_zz_0_yyz_0[i] * pb_x[i] - ta2_zz_0_yyz_1[i] * pc_x[i];

        ta2_zz_0_xyzz_0[i] = ta2_zz_0_yzz_0[i] * pb_x[i] - ta2_zz_0_yzz_1[i] * pc_x[i];

        ta2_zz_0_xzzz_0[i] = ta2_zz_0_zzz_0[i] * pb_x[i] - ta2_zz_0_zzz_1[i] * pc_x[i];

        ta2_zz_0_yyyy_0[i] = 3.0 * ta2_zz_0_yy_0[i] * fe_0 - 3.0 * ta2_zz_0_yy_1[i] * fe_0 + ta2_zz_0_yyy_0[i] * pb_y[i] - ta2_zz_0_yyy_1[i] * pc_y[i];

        ta2_zz_0_yyyz_0[i] = 2.0 * ta1_z_0_yyy_1[i] + ta2_zz_0_yyy_0[i] * pb_z[i] - ta2_zz_0_yyy_1[i] * pc_z[i];

        ta2_zz_0_yyzz_0[i] = ta2_zz_0_zz_0[i] * fe_0 - ta2_zz_0_zz_1[i] * fe_0 + ta2_zz_0_yzz_0[i] * pb_y[i] - ta2_zz_0_yzz_1[i] * pc_y[i];

        ta2_zz_0_yzzz_0[i] = ta2_zz_0_zzz_0[i] * pb_y[i] - ta2_zz_0_zzz_1[i] * pc_y[i];

        ta2_zz_0_zzzz_0[i] = 3.0 * ta2_zz_0_zz_0[i] * fe_0 - 3.0 * ta2_zz_0_zz_1[i] * fe_0 + 2.0 * ta1_z_0_zzz_1[i] + ta2_zz_0_zzz_0[i] * pb_z[i] - ta2_zz_0_zzz_1[i] * pc_z[i];
    }
}

} // npotrec namespace

