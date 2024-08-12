#include "NuclearPotentialGeom020PrimRecFS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_fs(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_fs,
                                        const size_t              idx_npot_geom_020_0_ps,
                                        const size_t              idx_npot_geom_020_1_ps,
                                        const size_t              idx_npot_geom_010_1_ds,
                                        const size_t              idx_npot_geom_020_0_ds,
                                        const size_t              idx_npot_geom_020_1_ds,
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

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds);

    auto ta1_x_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 3);

    auto ta1_x_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 5);

    auto ta1_y_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 6);

    auto ta1_y_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 9);

    auto ta1_y_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 11);

    auto ta1_z_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 12);

    auto ta1_z_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 15);

    auto ta1_z_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 17);

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds);

    auto ta2_xx_xz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 2);

    auto ta2_xx_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 3);

    auto ta2_xx_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 5);

    auto ta2_xy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 6);

    auto ta2_xy_xy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 7);

    auto ta2_xy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 9);

    auto ta2_xy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 11);

    auto ta2_xz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 12);

    auto ta2_xz_xz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 14);

    auto ta2_xz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 15);

    auto ta2_xz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 17);

    auto ta2_yy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 18);

    auto ta2_yy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 21);

    auto ta2_yy_yz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 22);

    auto ta2_yy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 23);

    auto ta2_yz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 24);

    auto ta2_yz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 27);

    auto ta2_yz_yz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 28);

    auto ta2_yz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 29);

    auto ta2_zz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 30);

    auto ta2_zz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 33);

    auto ta2_zz_yz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 34);

    auto ta2_zz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 35);

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds);

    auto ta2_xx_xz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 2);

    auto ta2_xx_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 3);

    auto ta2_xx_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 5);

    auto ta2_xy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 6);

    auto ta2_xy_xy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 7);

    auto ta2_xy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 9);

    auto ta2_xy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 11);

    auto ta2_xz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 12);

    auto ta2_xz_xz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 14);

    auto ta2_xz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 15);

    auto ta2_xz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 17);

    auto ta2_yy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 18);

    auto ta2_yy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 21);

    auto ta2_yy_yz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 22);

    auto ta2_yy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 23);

    auto ta2_yz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 24);

    auto ta2_yz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 27);

    auto ta2_yz_yz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 28);

    auto ta2_yz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 29);

    auto ta2_zz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 30);

    auto ta2_zz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 33);

    auto ta2_zz_yz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 34);

    auto ta2_zz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 35);

    // Set up components of targeted buffer : FS

    auto ta2_xx_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs);

    auto ta2_xx_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 1);

    auto ta2_xx_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 2);

    auto ta2_xx_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 3);

    auto ta2_xx_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 4);

    auto ta2_xx_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 5);

    auto ta2_xx_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 6);

    auto ta2_xx_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 7);

    auto ta2_xx_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 8);

    auto ta2_xx_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 9);

    auto ta2_xy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 10);

    auto ta2_xy_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 11);

    auto ta2_xy_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 12);

    auto ta2_xy_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 13);

    auto ta2_xy_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 14);

    auto ta2_xy_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 15);

    auto ta2_xy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 16);

    auto ta2_xy_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 17);

    auto ta2_xy_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 18);

    auto ta2_xy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 19);

    auto ta2_xz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 20);

    auto ta2_xz_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 21);

    auto ta2_xz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 22);

    auto ta2_xz_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 23);

    auto ta2_xz_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 24);

    auto ta2_xz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 25);

    auto ta2_xz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 26);

    auto ta2_xz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 27);

    auto ta2_xz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 28);

    auto ta2_xz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 29);

    auto ta2_yy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 30);

    auto ta2_yy_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 31);

    auto ta2_yy_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 32);

    auto ta2_yy_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 33);

    auto ta2_yy_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 34);

    auto ta2_yy_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 35);

    auto ta2_yy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 36);

    auto ta2_yy_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 37);

    auto ta2_yy_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 38);

    auto ta2_yy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 39);

    auto ta2_yz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 40);

    auto ta2_yz_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 41);

    auto ta2_yz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 42);

    auto ta2_yz_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 43);

    auto ta2_yz_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 44);

    auto ta2_yz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 45);

    auto ta2_yz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 46);

    auto ta2_yz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 47);

    auto ta2_yz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 48);

    auto ta2_yz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 49);

    auto ta2_zz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 50);

    auto ta2_zz_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 51);

    auto ta2_zz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 52);

    auto ta2_zz_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 53);

    auto ta2_zz_xyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 54);

    auto ta2_zz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 55);

    auto ta2_zz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 56);

    auto ta2_zz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 57);

    auto ta2_zz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 58);

    auto ta2_zz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_xx_0_1,   \
                             ta1_x_yy_0_1,   \
                             ta1_x_zz_0_1,   \
                             ta1_y_xx_0_1,   \
                             ta1_y_yy_0_1,   \
                             ta1_y_zz_0_1,   \
                             ta1_z_xx_0_1,   \
                             ta1_z_yy_0_1,   \
                             ta1_z_zz_0_1,   \
                             ta2_xx_x_0_0,   \
                             ta2_xx_x_0_1,   \
                             ta2_xx_xx_0_0,  \
                             ta2_xx_xx_0_1,  \
                             ta2_xx_xxx_0_0, \
                             ta2_xx_xxy_0_0, \
                             ta2_xx_xxz_0_0, \
                             ta2_xx_xyy_0_0, \
                             ta2_xx_xyz_0_0, \
                             ta2_xx_xz_0_0,  \
                             ta2_xx_xz_0_1,  \
                             ta2_xx_xzz_0_0, \
                             ta2_xx_y_0_0,   \
                             ta2_xx_y_0_1,   \
                             ta2_xx_yy_0_0,  \
                             ta2_xx_yy_0_1,  \
                             ta2_xx_yyy_0_0, \
                             ta2_xx_yyz_0_0, \
                             ta2_xx_yzz_0_0, \
                             ta2_xx_z_0_0,   \
                             ta2_xx_z_0_1,   \
                             ta2_xx_zz_0_0,  \
                             ta2_xx_zz_0_1,  \
                             ta2_xx_zzz_0_0, \
                             ta2_xy_x_0_0,   \
                             ta2_xy_x_0_1,   \
                             ta2_xy_xx_0_0,  \
                             ta2_xy_xx_0_1,  \
                             ta2_xy_xxx_0_0, \
                             ta2_xy_xxy_0_0, \
                             ta2_xy_xxz_0_0, \
                             ta2_xy_xy_0_0,  \
                             ta2_xy_xy_0_1,  \
                             ta2_xy_xyy_0_0, \
                             ta2_xy_xyz_0_0, \
                             ta2_xy_xzz_0_0, \
                             ta2_xy_y_0_0,   \
                             ta2_xy_y_0_1,   \
                             ta2_xy_yy_0_0,  \
                             ta2_xy_yy_0_1,  \
                             ta2_xy_yyy_0_0, \
                             ta2_xy_yyz_0_0, \
                             ta2_xy_yzz_0_0, \
                             ta2_xy_z_0_0,   \
                             ta2_xy_z_0_1,   \
                             ta2_xy_zz_0_0,  \
                             ta2_xy_zz_0_1,  \
                             ta2_xy_zzz_0_0, \
                             ta2_xz_x_0_0,   \
                             ta2_xz_x_0_1,   \
                             ta2_xz_xx_0_0,  \
                             ta2_xz_xx_0_1,  \
                             ta2_xz_xxx_0_0, \
                             ta2_xz_xxy_0_0, \
                             ta2_xz_xxz_0_0, \
                             ta2_xz_xyy_0_0, \
                             ta2_xz_xyz_0_0, \
                             ta2_xz_xz_0_0,  \
                             ta2_xz_xz_0_1,  \
                             ta2_xz_xzz_0_0, \
                             ta2_xz_y_0_0,   \
                             ta2_xz_y_0_1,   \
                             ta2_xz_yy_0_0,  \
                             ta2_xz_yy_0_1,  \
                             ta2_xz_yyy_0_0, \
                             ta2_xz_yyz_0_0, \
                             ta2_xz_yzz_0_0, \
                             ta2_xz_z_0_0,   \
                             ta2_xz_z_0_1,   \
                             ta2_xz_zz_0_0,  \
                             ta2_xz_zz_0_1,  \
                             ta2_xz_zzz_0_0, \
                             ta2_yy_x_0_0,   \
                             ta2_yy_x_0_1,   \
                             ta2_yy_xx_0_0,  \
                             ta2_yy_xx_0_1,  \
                             ta2_yy_xxx_0_0, \
                             ta2_yy_xxy_0_0, \
                             ta2_yy_xxz_0_0, \
                             ta2_yy_xyy_0_0, \
                             ta2_yy_xyz_0_0, \
                             ta2_yy_xzz_0_0, \
                             ta2_yy_y_0_0,   \
                             ta2_yy_y_0_1,   \
                             ta2_yy_yy_0_0,  \
                             ta2_yy_yy_0_1,  \
                             ta2_yy_yyy_0_0, \
                             ta2_yy_yyz_0_0, \
                             ta2_yy_yz_0_0,  \
                             ta2_yy_yz_0_1,  \
                             ta2_yy_yzz_0_0, \
                             ta2_yy_z_0_0,   \
                             ta2_yy_z_0_1,   \
                             ta2_yy_zz_0_0,  \
                             ta2_yy_zz_0_1,  \
                             ta2_yy_zzz_0_0, \
                             ta2_yz_x_0_0,   \
                             ta2_yz_x_0_1,   \
                             ta2_yz_xx_0_0,  \
                             ta2_yz_xx_0_1,  \
                             ta2_yz_xxx_0_0, \
                             ta2_yz_xxy_0_0, \
                             ta2_yz_xxz_0_0, \
                             ta2_yz_xyy_0_0, \
                             ta2_yz_xyz_0_0, \
                             ta2_yz_xzz_0_0, \
                             ta2_yz_y_0_0,   \
                             ta2_yz_y_0_1,   \
                             ta2_yz_yy_0_0,  \
                             ta2_yz_yy_0_1,  \
                             ta2_yz_yyy_0_0, \
                             ta2_yz_yyz_0_0, \
                             ta2_yz_yz_0_0,  \
                             ta2_yz_yz_0_1,  \
                             ta2_yz_yzz_0_0, \
                             ta2_yz_z_0_0,   \
                             ta2_yz_z_0_1,   \
                             ta2_yz_zz_0_0,  \
                             ta2_yz_zz_0_1,  \
                             ta2_yz_zzz_0_0, \
                             ta2_zz_x_0_0,   \
                             ta2_zz_x_0_1,   \
                             ta2_zz_xx_0_0,  \
                             ta2_zz_xx_0_1,  \
                             ta2_zz_xxx_0_0, \
                             ta2_zz_xxy_0_0, \
                             ta2_zz_xxz_0_0, \
                             ta2_zz_xyy_0_0, \
                             ta2_zz_xyz_0_0, \
                             ta2_zz_xzz_0_0, \
                             ta2_zz_y_0_0,   \
                             ta2_zz_y_0_1,   \
                             ta2_zz_yy_0_0,  \
                             ta2_zz_yy_0_1,  \
                             ta2_zz_yyy_0_0, \
                             ta2_zz_yyz_0_0, \
                             ta2_zz_yz_0_0,  \
                             ta2_zz_yz_0_1,  \
                             ta2_zz_yzz_0_0, \
                             ta2_zz_z_0_0,   \
                             ta2_zz_z_0_1,   \
                             ta2_zz_zz_0_0,  \
                             ta2_zz_zz_0_1,  \
                             ta2_zz_zzz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxx_0_0[i] = 2.0 * ta2_xx_x_0_0[i] * fe_0 - 2.0 * ta2_xx_x_0_1[i] * fe_0 + 2.0 * ta1_x_xx_0_1[i] + ta2_xx_xx_0_0[i] * pa_x[i] -
                            ta2_xx_xx_0_1[i] * pc_x[i];

        ta2_xx_xxy_0_0[i] = ta2_xx_xx_0_0[i] * pa_y[i] - ta2_xx_xx_0_1[i] * pc_y[i];

        ta2_xx_xxz_0_0[i] = ta2_xx_xx_0_0[i] * pa_z[i] - ta2_xx_xx_0_1[i] * pc_z[i];

        ta2_xx_xyy_0_0[i] = 2.0 * ta1_x_yy_0_1[i] + ta2_xx_yy_0_0[i] * pa_x[i] - ta2_xx_yy_0_1[i] * pc_x[i];

        ta2_xx_xyz_0_0[i] = ta2_xx_xz_0_0[i] * pa_y[i] - ta2_xx_xz_0_1[i] * pc_y[i];

        ta2_xx_xzz_0_0[i] = 2.0 * ta1_x_zz_0_1[i] + ta2_xx_zz_0_0[i] * pa_x[i] - ta2_xx_zz_0_1[i] * pc_x[i];

        ta2_xx_yyy_0_0[i] = 2.0 * ta2_xx_y_0_0[i] * fe_0 - 2.0 * ta2_xx_y_0_1[i] * fe_0 + ta2_xx_yy_0_0[i] * pa_y[i] - ta2_xx_yy_0_1[i] * pc_y[i];

        ta2_xx_yyz_0_0[i] = ta2_xx_yy_0_0[i] * pa_z[i] - ta2_xx_yy_0_1[i] * pc_z[i];

        ta2_xx_yzz_0_0[i] = ta2_xx_zz_0_0[i] * pa_y[i] - ta2_xx_zz_0_1[i] * pc_y[i];

        ta2_xx_zzz_0_0[i] = 2.0 * ta2_xx_z_0_0[i] * fe_0 - 2.0 * ta2_xx_z_0_1[i] * fe_0 + ta2_xx_zz_0_0[i] * pa_z[i] - ta2_xx_zz_0_1[i] * pc_z[i];

        ta2_xy_xxx_0_0[i] =
            2.0 * ta2_xy_x_0_0[i] * fe_0 - 2.0 * ta2_xy_x_0_1[i] * fe_0 + ta1_y_xx_0_1[i] + ta2_xy_xx_0_0[i] * pa_x[i] - ta2_xy_xx_0_1[i] * pc_x[i];

        ta2_xy_xxy_0_0[i] = ta1_x_xx_0_1[i] + ta2_xy_xx_0_0[i] * pa_y[i] - ta2_xy_xx_0_1[i] * pc_y[i];

        ta2_xy_xxz_0_0[i] = ta2_xy_xx_0_0[i] * pa_z[i] - ta2_xy_xx_0_1[i] * pc_z[i];

        ta2_xy_xyy_0_0[i] = ta1_y_yy_0_1[i] + ta2_xy_yy_0_0[i] * pa_x[i] - ta2_xy_yy_0_1[i] * pc_x[i];

        ta2_xy_xyz_0_0[i] = ta2_xy_xy_0_0[i] * pa_z[i] - ta2_xy_xy_0_1[i] * pc_z[i];

        ta2_xy_xzz_0_0[i] = ta1_y_zz_0_1[i] + ta2_xy_zz_0_0[i] * pa_x[i] - ta2_xy_zz_0_1[i] * pc_x[i];

        ta2_xy_yyy_0_0[i] =
            2.0 * ta2_xy_y_0_0[i] * fe_0 - 2.0 * ta2_xy_y_0_1[i] * fe_0 + ta1_x_yy_0_1[i] + ta2_xy_yy_0_0[i] * pa_y[i] - ta2_xy_yy_0_1[i] * pc_y[i];

        ta2_xy_yyz_0_0[i] = ta2_xy_yy_0_0[i] * pa_z[i] - ta2_xy_yy_0_1[i] * pc_z[i];

        ta2_xy_yzz_0_0[i] = ta1_x_zz_0_1[i] + ta2_xy_zz_0_0[i] * pa_y[i] - ta2_xy_zz_0_1[i] * pc_y[i];

        ta2_xy_zzz_0_0[i] = 2.0 * ta2_xy_z_0_0[i] * fe_0 - 2.0 * ta2_xy_z_0_1[i] * fe_0 + ta2_xy_zz_0_0[i] * pa_z[i] - ta2_xy_zz_0_1[i] * pc_z[i];

        ta2_xz_xxx_0_0[i] =
            2.0 * ta2_xz_x_0_0[i] * fe_0 - 2.0 * ta2_xz_x_0_1[i] * fe_0 + ta1_z_xx_0_1[i] + ta2_xz_xx_0_0[i] * pa_x[i] - ta2_xz_xx_0_1[i] * pc_x[i];

        ta2_xz_xxy_0_0[i] = ta2_xz_xx_0_0[i] * pa_y[i] - ta2_xz_xx_0_1[i] * pc_y[i];

        ta2_xz_xxz_0_0[i] = ta1_x_xx_0_1[i] + ta2_xz_xx_0_0[i] * pa_z[i] - ta2_xz_xx_0_1[i] * pc_z[i];

        ta2_xz_xyy_0_0[i] = ta1_z_yy_0_1[i] + ta2_xz_yy_0_0[i] * pa_x[i] - ta2_xz_yy_0_1[i] * pc_x[i];

        ta2_xz_xyz_0_0[i] = ta2_xz_xz_0_0[i] * pa_y[i] - ta2_xz_xz_0_1[i] * pc_y[i];

        ta2_xz_xzz_0_0[i] = ta1_z_zz_0_1[i] + ta2_xz_zz_0_0[i] * pa_x[i] - ta2_xz_zz_0_1[i] * pc_x[i];

        ta2_xz_yyy_0_0[i] = 2.0 * ta2_xz_y_0_0[i] * fe_0 - 2.0 * ta2_xz_y_0_1[i] * fe_0 + ta2_xz_yy_0_0[i] * pa_y[i] - ta2_xz_yy_0_1[i] * pc_y[i];

        ta2_xz_yyz_0_0[i] = ta1_x_yy_0_1[i] + ta2_xz_yy_0_0[i] * pa_z[i] - ta2_xz_yy_0_1[i] * pc_z[i];

        ta2_xz_yzz_0_0[i] = ta2_xz_zz_0_0[i] * pa_y[i] - ta2_xz_zz_0_1[i] * pc_y[i];

        ta2_xz_zzz_0_0[i] =
            2.0 * ta2_xz_z_0_0[i] * fe_0 - 2.0 * ta2_xz_z_0_1[i] * fe_0 + ta1_x_zz_0_1[i] + ta2_xz_zz_0_0[i] * pa_z[i] - ta2_xz_zz_0_1[i] * pc_z[i];

        ta2_yy_xxx_0_0[i] = 2.0 * ta2_yy_x_0_0[i] * fe_0 - 2.0 * ta2_yy_x_0_1[i] * fe_0 + ta2_yy_xx_0_0[i] * pa_x[i] - ta2_yy_xx_0_1[i] * pc_x[i];

        ta2_yy_xxy_0_0[i] = 2.0 * ta1_y_xx_0_1[i] + ta2_yy_xx_0_0[i] * pa_y[i] - ta2_yy_xx_0_1[i] * pc_y[i];

        ta2_yy_xxz_0_0[i] = ta2_yy_xx_0_0[i] * pa_z[i] - ta2_yy_xx_0_1[i] * pc_z[i];

        ta2_yy_xyy_0_0[i] = ta2_yy_yy_0_0[i] * pa_x[i] - ta2_yy_yy_0_1[i] * pc_x[i];

        ta2_yy_xyz_0_0[i] = ta2_yy_yz_0_0[i] * pa_x[i] - ta2_yy_yz_0_1[i] * pc_x[i];

        ta2_yy_xzz_0_0[i] = ta2_yy_zz_0_0[i] * pa_x[i] - ta2_yy_zz_0_1[i] * pc_x[i];

        ta2_yy_yyy_0_0[i] = 2.0 * ta2_yy_y_0_0[i] * fe_0 - 2.0 * ta2_yy_y_0_1[i] * fe_0 + 2.0 * ta1_y_yy_0_1[i] + ta2_yy_yy_0_0[i] * pa_y[i] -
                            ta2_yy_yy_0_1[i] * pc_y[i];

        ta2_yy_yyz_0_0[i] = ta2_yy_yy_0_0[i] * pa_z[i] - ta2_yy_yy_0_1[i] * pc_z[i];

        ta2_yy_yzz_0_0[i] = 2.0 * ta1_y_zz_0_1[i] + ta2_yy_zz_0_0[i] * pa_y[i] - ta2_yy_zz_0_1[i] * pc_y[i];

        ta2_yy_zzz_0_0[i] = 2.0 * ta2_yy_z_0_0[i] * fe_0 - 2.0 * ta2_yy_z_0_1[i] * fe_0 + ta2_yy_zz_0_0[i] * pa_z[i] - ta2_yy_zz_0_1[i] * pc_z[i];

        ta2_yz_xxx_0_0[i] = 2.0 * ta2_yz_x_0_0[i] * fe_0 - 2.0 * ta2_yz_x_0_1[i] * fe_0 + ta2_yz_xx_0_0[i] * pa_x[i] - ta2_yz_xx_0_1[i] * pc_x[i];

        ta2_yz_xxy_0_0[i] = ta1_z_xx_0_1[i] + ta2_yz_xx_0_0[i] * pa_y[i] - ta2_yz_xx_0_1[i] * pc_y[i];

        ta2_yz_xxz_0_0[i] = ta1_y_xx_0_1[i] + ta2_yz_xx_0_0[i] * pa_z[i] - ta2_yz_xx_0_1[i] * pc_z[i];

        ta2_yz_xyy_0_0[i] = ta2_yz_yy_0_0[i] * pa_x[i] - ta2_yz_yy_0_1[i] * pc_x[i];

        ta2_yz_xyz_0_0[i] = ta2_yz_yz_0_0[i] * pa_x[i] - ta2_yz_yz_0_1[i] * pc_x[i];

        ta2_yz_xzz_0_0[i] = ta2_yz_zz_0_0[i] * pa_x[i] - ta2_yz_zz_0_1[i] * pc_x[i];

        ta2_yz_yyy_0_0[i] =
            2.0 * ta2_yz_y_0_0[i] * fe_0 - 2.0 * ta2_yz_y_0_1[i] * fe_0 + ta1_z_yy_0_1[i] + ta2_yz_yy_0_0[i] * pa_y[i] - ta2_yz_yy_0_1[i] * pc_y[i];

        ta2_yz_yyz_0_0[i] = ta1_y_yy_0_1[i] + ta2_yz_yy_0_0[i] * pa_z[i] - ta2_yz_yy_0_1[i] * pc_z[i];

        ta2_yz_yzz_0_0[i] = ta1_z_zz_0_1[i] + ta2_yz_zz_0_0[i] * pa_y[i] - ta2_yz_zz_0_1[i] * pc_y[i];

        ta2_yz_zzz_0_0[i] =
            2.0 * ta2_yz_z_0_0[i] * fe_0 - 2.0 * ta2_yz_z_0_1[i] * fe_0 + ta1_y_zz_0_1[i] + ta2_yz_zz_0_0[i] * pa_z[i] - ta2_yz_zz_0_1[i] * pc_z[i];

        ta2_zz_xxx_0_0[i] = 2.0 * ta2_zz_x_0_0[i] * fe_0 - 2.0 * ta2_zz_x_0_1[i] * fe_0 + ta2_zz_xx_0_0[i] * pa_x[i] - ta2_zz_xx_0_1[i] * pc_x[i];

        ta2_zz_xxy_0_0[i] = ta2_zz_xx_0_0[i] * pa_y[i] - ta2_zz_xx_0_1[i] * pc_y[i];

        ta2_zz_xxz_0_0[i] = 2.0 * ta1_z_xx_0_1[i] + ta2_zz_xx_0_0[i] * pa_z[i] - ta2_zz_xx_0_1[i] * pc_z[i];

        ta2_zz_xyy_0_0[i] = ta2_zz_yy_0_0[i] * pa_x[i] - ta2_zz_yy_0_1[i] * pc_x[i];

        ta2_zz_xyz_0_0[i] = ta2_zz_yz_0_0[i] * pa_x[i] - ta2_zz_yz_0_1[i] * pc_x[i];

        ta2_zz_xzz_0_0[i] = ta2_zz_zz_0_0[i] * pa_x[i] - ta2_zz_zz_0_1[i] * pc_x[i];

        ta2_zz_yyy_0_0[i] = 2.0 * ta2_zz_y_0_0[i] * fe_0 - 2.0 * ta2_zz_y_0_1[i] * fe_0 + ta2_zz_yy_0_0[i] * pa_y[i] - ta2_zz_yy_0_1[i] * pc_y[i];

        ta2_zz_yyz_0_0[i] = 2.0 * ta1_z_yy_0_1[i] + ta2_zz_yy_0_0[i] * pa_z[i] - ta2_zz_yy_0_1[i] * pc_z[i];

        ta2_zz_yzz_0_0[i] = ta2_zz_zz_0_0[i] * pa_y[i] - ta2_zz_zz_0_1[i] * pc_y[i];

        ta2_zz_zzz_0_0[i] = 2.0 * ta2_zz_z_0_0[i] * fe_0 - 2.0 * ta2_zz_z_0_1[i] * fe_0 + 2.0 * ta1_z_zz_0_1[i] + ta2_zz_zz_0_0[i] * pa_z[i] -
                            ta2_zz_zz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
