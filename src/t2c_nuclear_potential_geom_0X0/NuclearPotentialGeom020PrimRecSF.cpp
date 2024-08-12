#include "NuclearPotentialGeom020PrimRecSF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_sf(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_sf,
                                        const size_t              idx_npot_geom_020_0_sp,
                                        const size_t              idx_npot_geom_020_1_sp,
                                        const size_t              idx_npot_geom_010_1_sd,
                                        const size_t              idx_npot_geom_020_0_sd,
                                        const size_t              idx_npot_geom_020_1_sd,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpb,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
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

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd);

    auto ta2_xx_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 2);

    auto ta2_xx_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 3);

    auto ta2_xx_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 5);

    auto ta2_xy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 6);

    auto ta2_xy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 7);

    auto ta2_xy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 9);

    auto ta2_xy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 11);

    auto ta2_xz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 12);

    auto ta2_xz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 14);

    auto ta2_xz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 15);

    auto ta2_xz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 17);

    auto ta2_yy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 18);

    auto ta2_yy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 21);

    auto ta2_yy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 22);

    auto ta2_yy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 23);

    auto ta2_yz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 24);

    auto ta2_yz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 27);

    auto ta2_yz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 28);

    auto ta2_yz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 29);

    auto ta2_zz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 30);

    auto ta2_zz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 33);

    auto ta2_zz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 34);

    auto ta2_zz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 35);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd);

    auto ta2_xx_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 2);

    auto ta2_xx_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 3);

    auto ta2_xx_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 5);

    auto ta2_xy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 6);

    auto ta2_xy_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 7);

    auto ta2_xy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 9);

    auto ta2_xy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 11);

    auto ta2_xz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 12);

    auto ta2_xz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 14);

    auto ta2_xz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 15);

    auto ta2_xz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 17);

    auto ta2_yy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 18);

    auto ta2_yy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 21);

    auto ta2_yy_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 22);

    auto ta2_yy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 23);

    auto ta2_yz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 24);

    auto ta2_yz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 27);

    auto ta2_yz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 28);

    auto ta2_yz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 29);

    auto ta2_zz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 30);

    auto ta2_zz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 33);

    auto ta2_zz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 34);

    auto ta2_zz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 35);

    // Set up components of targeted buffer : SF

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

#pragma omp simd aligned(pb_x,               \
                             pb_y,           \
                             pb_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_0_xx_1,   \
                             ta1_x_0_yy_1,   \
                             ta1_x_0_zz_1,   \
                             ta1_y_0_xx_1,   \
                             ta1_y_0_yy_1,   \
                             ta1_y_0_zz_1,   \
                             ta1_z_0_xx_1,   \
                             ta1_z_0_yy_1,   \
                             ta1_z_0_zz_1,   \
                             ta2_xx_0_x_0,   \
                             ta2_xx_0_x_1,   \
                             ta2_xx_0_xx_0,  \
                             ta2_xx_0_xx_1,  \
                             ta2_xx_0_xxx_0, \
                             ta2_xx_0_xxy_0, \
                             ta2_xx_0_xxz_0, \
                             ta2_xx_0_xyy_0, \
                             ta2_xx_0_xyz_0, \
                             ta2_xx_0_xz_0,  \
                             ta2_xx_0_xz_1,  \
                             ta2_xx_0_xzz_0, \
                             ta2_xx_0_y_0,   \
                             ta2_xx_0_y_1,   \
                             ta2_xx_0_yy_0,  \
                             ta2_xx_0_yy_1,  \
                             ta2_xx_0_yyy_0, \
                             ta2_xx_0_yyz_0, \
                             ta2_xx_0_yzz_0, \
                             ta2_xx_0_z_0,   \
                             ta2_xx_0_z_1,   \
                             ta2_xx_0_zz_0,  \
                             ta2_xx_0_zz_1,  \
                             ta2_xx_0_zzz_0, \
                             ta2_xy_0_x_0,   \
                             ta2_xy_0_x_1,   \
                             ta2_xy_0_xx_0,  \
                             ta2_xy_0_xx_1,  \
                             ta2_xy_0_xxx_0, \
                             ta2_xy_0_xxy_0, \
                             ta2_xy_0_xxz_0, \
                             ta2_xy_0_xy_0,  \
                             ta2_xy_0_xy_1,  \
                             ta2_xy_0_xyy_0, \
                             ta2_xy_0_xyz_0, \
                             ta2_xy_0_xzz_0, \
                             ta2_xy_0_y_0,   \
                             ta2_xy_0_y_1,   \
                             ta2_xy_0_yy_0,  \
                             ta2_xy_0_yy_1,  \
                             ta2_xy_0_yyy_0, \
                             ta2_xy_0_yyz_0, \
                             ta2_xy_0_yzz_0, \
                             ta2_xy_0_z_0,   \
                             ta2_xy_0_z_1,   \
                             ta2_xy_0_zz_0,  \
                             ta2_xy_0_zz_1,  \
                             ta2_xy_0_zzz_0, \
                             ta2_xz_0_x_0,   \
                             ta2_xz_0_x_1,   \
                             ta2_xz_0_xx_0,  \
                             ta2_xz_0_xx_1,  \
                             ta2_xz_0_xxx_0, \
                             ta2_xz_0_xxy_0, \
                             ta2_xz_0_xxz_0, \
                             ta2_xz_0_xyy_0, \
                             ta2_xz_0_xyz_0, \
                             ta2_xz_0_xz_0,  \
                             ta2_xz_0_xz_1,  \
                             ta2_xz_0_xzz_0, \
                             ta2_xz_0_y_0,   \
                             ta2_xz_0_y_1,   \
                             ta2_xz_0_yy_0,  \
                             ta2_xz_0_yy_1,  \
                             ta2_xz_0_yyy_0, \
                             ta2_xz_0_yyz_0, \
                             ta2_xz_0_yzz_0, \
                             ta2_xz_0_z_0,   \
                             ta2_xz_0_z_1,   \
                             ta2_xz_0_zz_0,  \
                             ta2_xz_0_zz_1,  \
                             ta2_xz_0_zzz_0, \
                             ta2_yy_0_x_0,   \
                             ta2_yy_0_x_1,   \
                             ta2_yy_0_xx_0,  \
                             ta2_yy_0_xx_1,  \
                             ta2_yy_0_xxx_0, \
                             ta2_yy_0_xxy_0, \
                             ta2_yy_0_xxz_0, \
                             ta2_yy_0_xyy_0, \
                             ta2_yy_0_xyz_0, \
                             ta2_yy_0_xzz_0, \
                             ta2_yy_0_y_0,   \
                             ta2_yy_0_y_1,   \
                             ta2_yy_0_yy_0,  \
                             ta2_yy_0_yy_1,  \
                             ta2_yy_0_yyy_0, \
                             ta2_yy_0_yyz_0, \
                             ta2_yy_0_yz_0,  \
                             ta2_yy_0_yz_1,  \
                             ta2_yy_0_yzz_0, \
                             ta2_yy_0_z_0,   \
                             ta2_yy_0_z_1,   \
                             ta2_yy_0_zz_0,  \
                             ta2_yy_0_zz_1,  \
                             ta2_yy_0_zzz_0, \
                             ta2_yz_0_x_0,   \
                             ta2_yz_0_x_1,   \
                             ta2_yz_0_xx_0,  \
                             ta2_yz_0_xx_1,  \
                             ta2_yz_0_xxx_0, \
                             ta2_yz_0_xxy_0, \
                             ta2_yz_0_xxz_0, \
                             ta2_yz_0_xyy_0, \
                             ta2_yz_0_xyz_0, \
                             ta2_yz_0_xzz_0, \
                             ta2_yz_0_y_0,   \
                             ta2_yz_0_y_1,   \
                             ta2_yz_0_yy_0,  \
                             ta2_yz_0_yy_1,  \
                             ta2_yz_0_yyy_0, \
                             ta2_yz_0_yyz_0, \
                             ta2_yz_0_yz_0,  \
                             ta2_yz_0_yz_1,  \
                             ta2_yz_0_yzz_0, \
                             ta2_yz_0_z_0,   \
                             ta2_yz_0_z_1,   \
                             ta2_yz_0_zz_0,  \
                             ta2_yz_0_zz_1,  \
                             ta2_yz_0_zzz_0, \
                             ta2_zz_0_x_0,   \
                             ta2_zz_0_x_1,   \
                             ta2_zz_0_xx_0,  \
                             ta2_zz_0_xx_1,  \
                             ta2_zz_0_xxx_0, \
                             ta2_zz_0_xxy_0, \
                             ta2_zz_0_xxz_0, \
                             ta2_zz_0_xyy_0, \
                             ta2_zz_0_xyz_0, \
                             ta2_zz_0_xzz_0, \
                             ta2_zz_0_y_0,   \
                             ta2_zz_0_y_1,   \
                             ta2_zz_0_yy_0,  \
                             ta2_zz_0_yy_1,  \
                             ta2_zz_0_yyy_0, \
                             ta2_zz_0_yyz_0, \
                             ta2_zz_0_yz_0,  \
                             ta2_zz_0_yz_1,  \
                             ta2_zz_0_yzz_0, \
                             ta2_zz_0_z_0,   \
                             ta2_zz_0_z_1,   \
                             ta2_zz_0_zz_0,  \
                             ta2_zz_0_zz_1,  \
                             ta2_zz_0_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_0_xxx_0[i] = 2.0 * ta2_xx_0_x_0[i] * fe_0 - 2.0 * ta2_xx_0_x_1[i] * fe_0 + 2.0 * ta1_x_0_xx_1[i] + ta2_xx_0_xx_0[i] * pb_x[i] -
                            ta2_xx_0_xx_1[i] * pc_x[i];

        ta2_xx_0_xxy_0[i] = ta2_xx_0_xx_0[i] * pb_y[i] - ta2_xx_0_xx_1[i] * pc_y[i];

        ta2_xx_0_xxz_0[i] = ta2_xx_0_xx_0[i] * pb_z[i] - ta2_xx_0_xx_1[i] * pc_z[i];

        ta2_xx_0_xyy_0[i] = 2.0 * ta1_x_0_yy_1[i] + ta2_xx_0_yy_0[i] * pb_x[i] - ta2_xx_0_yy_1[i] * pc_x[i];

        ta2_xx_0_xyz_0[i] = ta2_xx_0_xz_0[i] * pb_y[i] - ta2_xx_0_xz_1[i] * pc_y[i];

        ta2_xx_0_xzz_0[i] = 2.0 * ta1_x_0_zz_1[i] + ta2_xx_0_zz_0[i] * pb_x[i] - ta2_xx_0_zz_1[i] * pc_x[i];

        ta2_xx_0_yyy_0[i] = 2.0 * ta2_xx_0_y_0[i] * fe_0 - 2.0 * ta2_xx_0_y_1[i] * fe_0 + ta2_xx_0_yy_0[i] * pb_y[i] - ta2_xx_0_yy_1[i] * pc_y[i];

        ta2_xx_0_yyz_0[i] = ta2_xx_0_yy_0[i] * pb_z[i] - ta2_xx_0_yy_1[i] * pc_z[i];

        ta2_xx_0_yzz_0[i] = ta2_xx_0_zz_0[i] * pb_y[i] - ta2_xx_0_zz_1[i] * pc_y[i];

        ta2_xx_0_zzz_0[i] = 2.0 * ta2_xx_0_z_0[i] * fe_0 - 2.0 * ta2_xx_0_z_1[i] * fe_0 + ta2_xx_0_zz_0[i] * pb_z[i] - ta2_xx_0_zz_1[i] * pc_z[i];

        ta2_xy_0_xxx_0[i] =
            2.0 * ta2_xy_0_x_0[i] * fe_0 - 2.0 * ta2_xy_0_x_1[i] * fe_0 + ta1_y_0_xx_1[i] + ta2_xy_0_xx_0[i] * pb_x[i] - ta2_xy_0_xx_1[i] * pc_x[i];

        ta2_xy_0_xxy_0[i] = ta1_x_0_xx_1[i] + ta2_xy_0_xx_0[i] * pb_y[i] - ta2_xy_0_xx_1[i] * pc_y[i];

        ta2_xy_0_xxz_0[i] = ta2_xy_0_xx_0[i] * pb_z[i] - ta2_xy_0_xx_1[i] * pc_z[i];

        ta2_xy_0_xyy_0[i] = ta1_y_0_yy_1[i] + ta2_xy_0_yy_0[i] * pb_x[i] - ta2_xy_0_yy_1[i] * pc_x[i];

        ta2_xy_0_xyz_0[i] = ta2_xy_0_xy_0[i] * pb_z[i] - ta2_xy_0_xy_1[i] * pc_z[i];

        ta2_xy_0_xzz_0[i] = ta1_y_0_zz_1[i] + ta2_xy_0_zz_0[i] * pb_x[i] - ta2_xy_0_zz_1[i] * pc_x[i];

        ta2_xy_0_yyy_0[i] =
            2.0 * ta2_xy_0_y_0[i] * fe_0 - 2.0 * ta2_xy_0_y_1[i] * fe_0 + ta1_x_0_yy_1[i] + ta2_xy_0_yy_0[i] * pb_y[i] - ta2_xy_0_yy_1[i] * pc_y[i];

        ta2_xy_0_yyz_0[i] = ta2_xy_0_yy_0[i] * pb_z[i] - ta2_xy_0_yy_1[i] * pc_z[i];

        ta2_xy_0_yzz_0[i] = ta1_x_0_zz_1[i] + ta2_xy_0_zz_0[i] * pb_y[i] - ta2_xy_0_zz_1[i] * pc_y[i];

        ta2_xy_0_zzz_0[i] = 2.0 * ta2_xy_0_z_0[i] * fe_0 - 2.0 * ta2_xy_0_z_1[i] * fe_0 + ta2_xy_0_zz_0[i] * pb_z[i] - ta2_xy_0_zz_1[i] * pc_z[i];

        ta2_xz_0_xxx_0[i] =
            2.0 * ta2_xz_0_x_0[i] * fe_0 - 2.0 * ta2_xz_0_x_1[i] * fe_0 + ta1_z_0_xx_1[i] + ta2_xz_0_xx_0[i] * pb_x[i] - ta2_xz_0_xx_1[i] * pc_x[i];

        ta2_xz_0_xxy_0[i] = ta2_xz_0_xx_0[i] * pb_y[i] - ta2_xz_0_xx_1[i] * pc_y[i];

        ta2_xz_0_xxz_0[i] = ta1_x_0_xx_1[i] + ta2_xz_0_xx_0[i] * pb_z[i] - ta2_xz_0_xx_1[i] * pc_z[i];

        ta2_xz_0_xyy_0[i] = ta1_z_0_yy_1[i] + ta2_xz_0_yy_0[i] * pb_x[i] - ta2_xz_0_yy_1[i] * pc_x[i];

        ta2_xz_0_xyz_0[i] = ta2_xz_0_xz_0[i] * pb_y[i] - ta2_xz_0_xz_1[i] * pc_y[i];

        ta2_xz_0_xzz_0[i] = ta1_z_0_zz_1[i] + ta2_xz_0_zz_0[i] * pb_x[i] - ta2_xz_0_zz_1[i] * pc_x[i];

        ta2_xz_0_yyy_0[i] = 2.0 * ta2_xz_0_y_0[i] * fe_0 - 2.0 * ta2_xz_0_y_1[i] * fe_0 + ta2_xz_0_yy_0[i] * pb_y[i] - ta2_xz_0_yy_1[i] * pc_y[i];

        ta2_xz_0_yyz_0[i] = ta1_x_0_yy_1[i] + ta2_xz_0_yy_0[i] * pb_z[i] - ta2_xz_0_yy_1[i] * pc_z[i];

        ta2_xz_0_yzz_0[i] = ta2_xz_0_zz_0[i] * pb_y[i] - ta2_xz_0_zz_1[i] * pc_y[i];

        ta2_xz_0_zzz_0[i] =
            2.0 * ta2_xz_0_z_0[i] * fe_0 - 2.0 * ta2_xz_0_z_1[i] * fe_0 + ta1_x_0_zz_1[i] + ta2_xz_0_zz_0[i] * pb_z[i] - ta2_xz_0_zz_1[i] * pc_z[i];

        ta2_yy_0_xxx_0[i] = 2.0 * ta2_yy_0_x_0[i] * fe_0 - 2.0 * ta2_yy_0_x_1[i] * fe_0 + ta2_yy_0_xx_0[i] * pb_x[i] - ta2_yy_0_xx_1[i] * pc_x[i];

        ta2_yy_0_xxy_0[i] = 2.0 * ta1_y_0_xx_1[i] + ta2_yy_0_xx_0[i] * pb_y[i] - ta2_yy_0_xx_1[i] * pc_y[i];

        ta2_yy_0_xxz_0[i] = ta2_yy_0_xx_0[i] * pb_z[i] - ta2_yy_0_xx_1[i] * pc_z[i];

        ta2_yy_0_xyy_0[i] = ta2_yy_0_yy_0[i] * pb_x[i] - ta2_yy_0_yy_1[i] * pc_x[i];

        ta2_yy_0_xyz_0[i] = ta2_yy_0_yz_0[i] * pb_x[i] - ta2_yy_0_yz_1[i] * pc_x[i];

        ta2_yy_0_xzz_0[i] = ta2_yy_0_zz_0[i] * pb_x[i] - ta2_yy_0_zz_1[i] * pc_x[i];

        ta2_yy_0_yyy_0[i] = 2.0 * ta2_yy_0_y_0[i] * fe_0 - 2.0 * ta2_yy_0_y_1[i] * fe_0 + 2.0 * ta1_y_0_yy_1[i] + ta2_yy_0_yy_0[i] * pb_y[i] -
                            ta2_yy_0_yy_1[i] * pc_y[i];

        ta2_yy_0_yyz_0[i] = ta2_yy_0_yy_0[i] * pb_z[i] - ta2_yy_0_yy_1[i] * pc_z[i];

        ta2_yy_0_yzz_0[i] = 2.0 * ta1_y_0_zz_1[i] + ta2_yy_0_zz_0[i] * pb_y[i] - ta2_yy_0_zz_1[i] * pc_y[i];

        ta2_yy_0_zzz_0[i] = 2.0 * ta2_yy_0_z_0[i] * fe_0 - 2.0 * ta2_yy_0_z_1[i] * fe_0 + ta2_yy_0_zz_0[i] * pb_z[i] - ta2_yy_0_zz_1[i] * pc_z[i];

        ta2_yz_0_xxx_0[i] = 2.0 * ta2_yz_0_x_0[i] * fe_0 - 2.0 * ta2_yz_0_x_1[i] * fe_0 + ta2_yz_0_xx_0[i] * pb_x[i] - ta2_yz_0_xx_1[i] * pc_x[i];

        ta2_yz_0_xxy_0[i] = ta1_z_0_xx_1[i] + ta2_yz_0_xx_0[i] * pb_y[i] - ta2_yz_0_xx_1[i] * pc_y[i];

        ta2_yz_0_xxz_0[i] = ta1_y_0_xx_1[i] + ta2_yz_0_xx_0[i] * pb_z[i] - ta2_yz_0_xx_1[i] * pc_z[i];

        ta2_yz_0_xyy_0[i] = ta2_yz_0_yy_0[i] * pb_x[i] - ta2_yz_0_yy_1[i] * pc_x[i];

        ta2_yz_0_xyz_0[i] = ta2_yz_0_yz_0[i] * pb_x[i] - ta2_yz_0_yz_1[i] * pc_x[i];

        ta2_yz_0_xzz_0[i] = ta2_yz_0_zz_0[i] * pb_x[i] - ta2_yz_0_zz_1[i] * pc_x[i];

        ta2_yz_0_yyy_0[i] =
            2.0 * ta2_yz_0_y_0[i] * fe_0 - 2.0 * ta2_yz_0_y_1[i] * fe_0 + ta1_z_0_yy_1[i] + ta2_yz_0_yy_0[i] * pb_y[i] - ta2_yz_0_yy_1[i] * pc_y[i];

        ta2_yz_0_yyz_0[i] = ta1_y_0_yy_1[i] + ta2_yz_0_yy_0[i] * pb_z[i] - ta2_yz_0_yy_1[i] * pc_z[i];

        ta2_yz_0_yzz_0[i] = ta1_z_0_zz_1[i] + ta2_yz_0_zz_0[i] * pb_y[i] - ta2_yz_0_zz_1[i] * pc_y[i];

        ta2_yz_0_zzz_0[i] =
            2.0 * ta2_yz_0_z_0[i] * fe_0 - 2.0 * ta2_yz_0_z_1[i] * fe_0 + ta1_y_0_zz_1[i] + ta2_yz_0_zz_0[i] * pb_z[i] - ta2_yz_0_zz_1[i] * pc_z[i];

        ta2_zz_0_xxx_0[i] = 2.0 * ta2_zz_0_x_0[i] * fe_0 - 2.0 * ta2_zz_0_x_1[i] * fe_0 + ta2_zz_0_xx_0[i] * pb_x[i] - ta2_zz_0_xx_1[i] * pc_x[i];

        ta2_zz_0_xxy_0[i] = ta2_zz_0_xx_0[i] * pb_y[i] - ta2_zz_0_xx_1[i] * pc_y[i];

        ta2_zz_0_xxz_0[i] = 2.0 * ta1_z_0_xx_1[i] + ta2_zz_0_xx_0[i] * pb_z[i] - ta2_zz_0_xx_1[i] * pc_z[i];

        ta2_zz_0_xyy_0[i] = ta2_zz_0_yy_0[i] * pb_x[i] - ta2_zz_0_yy_1[i] * pc_x[i];

        ta2_zz_0_xyz_0[i] = ta2_zz_0_yz_0[i] * pb_x[i] - ta2_zz_0_yz_1[i] * pc_x[i];

        ta2_zz_0_xzz_0[i] = ta2_zz_0_zz_0[i] * pb_x[i] - ta2_zz_0_zz_1[i] * pc_x[i];

        ta2_zz_0_yyy_0[i] = 2.0 * ta2_zz_0_y_0[i] * fe_0 - 2.0 * ta2_zz_0_y_1[i] * fe_0 + ta2_zz_0_yy_0[i] * pb_y[i] - ta2_zz_0_yy_1[i] * pc_y[i];

        ta2_zz_0_yyz_0[i] = 2.0 * ta1_z_0_yy_1[i] + ta2_zz_0_yy_0[i] * pb_z[i] - ta2_zz_0_yy_1[i] * pc_z[i];

        ta2_zz_0_yzz_0[i] = ta2_zz_0_zz_0[i] * pb_y[i] - ta2_zz_0_zz_1[i] * pc_y[i];

        ta2_zz_0_zzz_0[i] = 2.0 * ta2_zz_0_z_0[i] * fe_0 - 2.0 * ta2_zz_0_z_1[i] * fe_0 + 2.0 * ta1_z_0_zz_1[i] + ta2_zz_0_zz_0[i] * pb_z[i] -
                            ta2_zz_0_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
