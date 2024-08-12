#include "NuclearPotentialGeom020PrimRecGP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_gp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_gp,
                                        const size_t              idx_npot_geom_020_0_dp,
                                        const size_t              idx_npot_geom_020_1_dp,
                                        const size_t              idx_npot_geom_020_0_fs,
                                        const size_t              idx_npot_geom_020_1_fs,
                                        const size_t              idx_npot_geom_010_1_fp,
                                        const size_t              idx_npot_geom_020_0_fp,
                                        const size_t              idx_npot_geom_020_1_fp,
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

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp);

    auto ta2_xx_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 1);

    auto ta2_xx_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 2);

    auto ta2_xx_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 3);

    auto ta2_xx_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 6);

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

    auto ta2_xy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 22);

    auto ta2_xy_xz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 24);

    auto ta2_xy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 27);

    auto ta2_xy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 28);

    auto ta2_xy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 29);

    auto ta2_xy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 31);

    auto ta2_xy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 33);

    auto ta2_xy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 34);

    auto ta2_xy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 35);

    auto ta2_xz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 36);

    auto ta2_xz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 37);

    auto ta2_xz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 38);

    auto ta2_xz_xy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 39);

    auto ta2_xz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 44);

    auto ta2_xz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 45);

    auto ta2_xz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 46);

    auto ta2_xz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 47);

    auto ta2_xz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 50);

    auto ta2_xz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 51);

    auto ta2_xz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 52);

    auto ta2_xz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 53);

    auto ta2_yy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 54);

    auto ta2_yy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 55);

    auto ta2_yy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 56);

    auto ta2_yy_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 58);

    auto ta2_yy_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 62);

    auto ta2_yy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 63);

    auto ta2_yy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 64);

    auto ta2_yy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 65);

    auto ta2_yy_yz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 67);

    auto ta2_yy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 69);

    auto ta2_yy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 70);

    auto ta2_yy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 71);

    auto ta2_yz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 72);

    auto ta2_yz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 73);

    auto ta2_yz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 74);

    auto ta2_yz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 76);

    auto ta2_yz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 80);

    auto ta2_yz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 81);

    auto ta2_yz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 82);

    auto ta2_yz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 83);

    auto ta2_yz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 86);

    auto ta2_yz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 87);

    auto ta2_yz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 88);

    auto ta2_yz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 89);

    auto ta2_zz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 90);

    auto ta2_zz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 91);

    auto ta2_zz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 92);

    auto ta2_zz_xy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 94);

    auto ta2_zz_xz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 98);

    auto ta2_zz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 99);

    auto ta2_zz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 100);

    auto ta2_zz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 101);

    auto ta2_zz_yz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 104);

    auto ta2_zz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 105);

    auto ta2_zz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 106);

    auto ta2_zz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 107);

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp);

    auto ta2_xx_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 1);

    auto ta2_xx_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 2);

    auto ta2_xx_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 3);

    auto ta2_xx_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 6);

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

    auto ta2_xy_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 22);

    auto ta2_xy_xz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 24);

    auto ta2_xy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 27);

    auto ta2_xy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 28);

    auto ta2_xy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 29);

    auto ta2_xy_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 31);

    auto ta2_xy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 33);

    auto ta2_xy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 34);

    auto ta2_xy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 35);

    auto ta2_xz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 36);

    auto ta2_xz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 37);

    auto ta2_xz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 38);

    auto ta2_xz_xy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 39);

    auto ta2_xz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 44);

    auto ta2_xz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 45);

    auto ta2_xz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 46);

    auto ta2_xz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 47);

    auto ta2_xz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 50);

    auto ta2_xz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 51);

    auto ta2_xz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 52);

    auto ta2_xz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 53);

    auto ta2_yy_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 54);

    auto ta2_yy_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 55);

    auto ta2_yy_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 56);

    auto ta2_yy_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 58);

    auto ta2_yy_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 62);

    auto ta2_yy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 63);

    auto ta2_yy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 64);

    auto ta2_yy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 65);

    auto ta2_yy_yz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 67);

    auto ta2_yy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 69);

    auto ta2_yy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 70);

    auto ta2_yy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 71);

    auto ta2_yz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 72);

    auto ta2_yz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 73);

    auto ta2_yz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 74);

    auto ta2_yz_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 76);

    auto ta2_yz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 80);

    auto ta2_yz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 81);

    auto ta2_yz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 82);

    auto ta2_yz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 83);

    auto ta2_yz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 86);

    auto ta2_yz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 87);

    auto ta2_yz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 88);

    auto ta2_yz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 89);

    auto ta2_zz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 90);

    auto ta2_zz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 91);

    auto ta2_zz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 92);

    auto ta2_zz_xy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 94);

    auto ta2_zz_xz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 98);

    auto ta2_zz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 99);

    auto ta2_zz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 100);

    auto ta2_zz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 101);

    auto ta2_zz_yz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 104);

    auto ta2_zz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 105);

    auto ta2_zz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 106);

    auto ta2_zz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 107);

    // Set up components of auxiliary buffer : FS

    auto ta2_xx_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs);

    auto ta2_xx_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 6);

    auto ta2_xx_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 9);

    auto ta2_xy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 10);

    auto ta2_xy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 16);

    auto ta2_xy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 19);

    auto ta2_xz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 20);

    auto ta2_xz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 26);

    auto ta2_xz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 29);

    auto ta2_yy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 30);

    auto ta2_yy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 36);

    auto ta2_yy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 39);

    auto ta2_yz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 40);

    auto ta2_yz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 46);

    auto ta2_yz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 49);

    auto ta2_zz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 50);

    auto ta2_zz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 56);

    auto ta2_zz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 59);

    // Set up components of auxiliary buffer : FS

    auto ta2_xx_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs);

    auto ta2_xx_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 6);

    auto ta2_xx_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 9);

    auto ta2_xy_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 10);

    auto ta2_xy_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 16);

    auto ta2_xy_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 19);

    auto ta2_xz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 20);

    auto ta2_xz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 26);

    auto ta2_xz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 29);

    auto ta2_yy_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 30);

    auto ta2_yy_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 36);

    auto ta2_yy_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 39);

    auto ta2_yz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 40);

    auto ta2_yz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 46);

    auto ta2_yz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 49);

    auto ta2_zz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 50);

    auto ta2_zz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 56);

    auto ta2_zz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 59);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp);

    auto ta1_x_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 1);

    auto ta1_x_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 2);

    auto ta1_x_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 3);

    auto ta1_x_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 4);

    auto ta1_x_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 6);

    auto ta1_x_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 8);

    auto ta1_x_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 9);

    auto ta1_x_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 10);

    auto ta1_x_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 15);

    auto ta1_x_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 17);

    auto ta1_x_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 18);

    auto ta1_x_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 19);

    auto ta1_x_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 20);

    auto ta1_x_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 22);

    auto ta1_x_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 23);

    auto ta1_x_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 25);

    auto ta1_x_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 26);

    auto ta1_x_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 27);

    auto ta1_x_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 28);

    auto ta1_x_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 29);

    auto ta1_y_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 30);

    auto ta1_y_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 31);

    auto ta1_y_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 32);

    auto ta1_y_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 33);

    auto ta1_y_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 34);

    auto ta1_y_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 36);

    auto ta1_y_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 38);

    auto ta1_y_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 39);

    auto ta1_y_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 40);

    auto ta1_y_xyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 41);

    auto ta1_y_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 45);

    auto ta1_y_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 47);

    auto ta1_y_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 48);

    auto ta1_y_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 49);

    auto ta1_y_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 50);

    auto ta1_y_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 52);

    auto ta1_y_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 53);

    auto ta1_y_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 55);

    auto ta1_y_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 56);

    auto ta1_y_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 57);

    auto ta1_y_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 58);

    auto ta1_y_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 59);

    auto ta1_z_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 60);

    auto ta1_z_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 61);

    auto ta1_z_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 62);

    auto ta1_z_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 63);

    auto ta1_z_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 64);

    auto ta1_z_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 66);

    auto ta1_z_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 68);

    auto ta1_z_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 69);

    auto ta1_z_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 70);

    auto ta1_z_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 75);

    auto ta1_z_xzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 76);

    auto ta1_z_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 77);

    auto ta1_z_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 78);

    auto ta1_z_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 79);

    auto ta1_z_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 80);

    auto ta1_z_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 82);

    auto ta1_z_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 83);

    auto ta1_z_yzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 84);

    auto ta1_z_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 85);

    auto ta1_z_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 86);

    auto ta1_z_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 87);

    auto ta1_z_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 88);

    auto ta1_z_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 89);

    // Set up components of auxiliary buffer : FP

    auto ta2_xx_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp);

    auto ta2_xx_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 1);

    auto ta2_xx_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 2);

    auto ta2_xx_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 3);

    auto ta2_xx_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 4);

    auto ta2_xx_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 5);

    auto ta2_xx_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 6);

    auto ta2_xx_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 7);

    auto ta2_xx_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 8);

    auto ta2_xx_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 9);

    auto ta2_xx_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 10);

    auto ta2_xx_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 15);

    auto ta2_xx_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 17);

    auto ta2_xx_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 18);

    auto ta2_xx_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 19);

    auto ta2_xx_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 20);

    auto ta2_xx_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 22);

    auto ta2_xx_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 23);

    auto ta2_xx_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 24);

    auto ta2_xx_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 25);

    auto ta2_xx_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 26);

    auto ta2_xx_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 27);

    auto ta2_xx_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 28);

    auto ta2_xx_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 29);

    auto ta2_xy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 30);

    auto ta2_xy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 31);

    auto ta2_xy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 32);

    auto ta2_xy_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 33);

    auto ta2_xy_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 34);

    auto ta2_xy_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 36);

    auto ta2_xy_xxz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 37);

    auto ta2_xy_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 38);

    auto ta2_xy_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 39);

    auto ta2_xy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 40);

    auto ta2_xy_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 41);

    auto ta2_xy_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 45);

    auto ta2_xy_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 47);

    auto ta2_xy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 48);

    auto ta2_xy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 49);

    auto ta2_xy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 50);

    auto ta2_xy_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 51);

    auto ta2_xy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 52);

    auto ta2_xy_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 53);

    auto ta2_xy_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 55);

    auto ta2_xy_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 56);

    auto ta2_xy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 57);

    auto ta2_xy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 58);

    auto ta2_xy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 59);

    auto ta2_xz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 60);

    auto ta2_xz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 61);

    auto ta2_xz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 62);

    auto ta2_xz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 63);

    auto ta2_xz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 64);

    auto ta2_xz_xxy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 65);

    auto ta2_xz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 66);

    auto ta2_xz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 68);

    auto ta2_xz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 69);

    auto ta2_xz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 70);

    auto ta2_xz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 75);

    auto ta2_xz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 76);

    auto ta2_xz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 77);

    auto ta2_xz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 78);

    auto ta2_xz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 79);

    auto ta2_xz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 80);

    auto ta2_xz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 82);

    auto ta2_xz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 83);

    auto ta2_xz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 84);

    auto ta2_xz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 85);

    auto ta2_xz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 86);

    auto ta2_xz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 87);

    auto ta2_xz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 88);

    auto ta2_xz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 89);

    auto ta2_yy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 90);

    auto ta2_yy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 91);

    auto ta2_yy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 92);

    auto ta2_yy_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 93);

    auto ta2_yy_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 94);

    auto ta2_yy_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 96);

    auto ta2_yy_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 98);

    auto ta2_yy_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 99);

    auto ta2_yy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 100);

    auto ta2_yy_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 101);

    auto ta2_yy_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 105);

    auto ta2_yy_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 106);

    auto ta2_yy_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 107);

    auto ta2_yy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 108);

    auto ta2_yy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 109);

    auto ta2_yy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 110);

    auto ta2_yy_yyz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 111);

    auto ta2_yy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 112);

    auto ta2_yy_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 113);

    auto ta2_yy_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 115);

    auto ta2_yy_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 116);

    auto ta2_yy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 117);

    auto ta2_yy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 118);

    auto ta2_yy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 119);

    auto ta2_yz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 120);

    auto ta2_yz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 121);

    auto ta2_yz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 122);

    auto ta2_yz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 123);

    auto ta2_yz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 124);

    auto ta2_yz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 126);

    auto ta2_yz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 128);

    auto ta2_yz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 129);

    auto ta2_yz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 130);

    auto ta2_yz_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 131);

    auto ta2_yz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 135);

    auto ta2_yz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 136);

    auto ta2_yz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 137);

    auto ta2_yz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 138);

    auto ta2_yz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 139);

    auto ta2_yz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 140);

    auto ta2_yz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 142);

    auto ta2_yz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 143);

    auto ta2_yz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 144);

    auto ta2_yz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 145);

    auto ta2_yz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 146);

    auto ta2_yz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 147);

    auto ta2_yz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 148);

    auto ta2_yz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 149);

    auto ta2_zz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 150);

    auto ta2_zz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 151);

    auto ta2_zz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 152);

    auto ta2_zz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 153);

    auto ta2_zz_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 154);

    auto ta2_zz_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 156);

    auto ta2_zz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 158);

    auto ta2_zz_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 159);

    auto ta2_zz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 160);

    auto ta2_zz_xyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 161);

    auto ta2_zz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 165);

    auto ta2_zz_xzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 166);

    auto ta2_zz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 167);

    auto ta2_zz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 168);

    auto ta2_zz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 169);

    auto ta2_zz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 170);

    auto ta2_zz_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 172);

    auto ta2_zz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 173);

    auto ta2_zz_yzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 174);

    auto ta2_zz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 175);

    auto ta2_zz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 176);

    auto ta2_zz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 177);

    auto ta2_zz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 178);

    auto ta2_zz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 179);

    // Set up components of auxiliary buffer : FP

    auto ta2_xx_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp);

    auto ta2_xx_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 1);

    auto ta2_xx_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 2);

    auto ta2_xx_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 3);

    auto ta2_xx_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 4);

    auto ta2_xx_xxy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 5);

    auto ta2_xx_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 6);

    auto ta2_xx_xxz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 7);

    auto ta2_xx_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 8);

    auto ta2_xx_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 9);

    auto ta2_xx_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 10);

    auto ta2_xx_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 15);

    auto ta2_xx_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 17);

    auto ta2_xx_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 18);

    auto ta2_xx_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 19);

    auto ta2_xx_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 20);

    auto ta2_xx_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 22);

    auto ta2_xx_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 23);

    auto ta2_xx_yzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 24);

    auto ta2_xx_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 25);

    auto ta2_xx_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 26);

    auto ta2_xx_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 27);

    auto ta2_xx_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 28);

    auto ta2_xx_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 29);

    auto ta2_xy_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 30);

    auto ta2_xy_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 31);

    auto ta2_xy_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 32);

    auto ta2_xy_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 33);

    auto ta2_xy_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 34);

    auto ta2_xy_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 36);

    auto ta2_xy_xxz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 37);

    auto ta2_xy_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 38);

    auto ta2_xy_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 39);

    auto ta2_xy_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 40);

    auto ta2_xy_xyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 41);

    auto ta2_xy_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 45);

    auto ta2_xy_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 47);

    auto ta2_xy_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 48);

    auto ta2_xy_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 49);

    auto ta2_xy_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 50);

    auto ta2_xy_yyz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 51);

    auto ta2_xy_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 52);

    auto ta2_xy_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 53);

    auto ta2_xy_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 55);

    auto ta2_xy_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 56);

    auto ta2_xy_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 57);

    auto ta2_xy_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 58);

    auto ta2_xy_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 59);

    auto ta2_xz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 60);

    auto ta2_xz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 61);

    auto ta2_xz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 62);

    auto ta2_xz_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 63);

    auto ta2_xz_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 64);

    auto ta2_xz_xxy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 65);

    auto ta2_xz_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 66);

    auto ta2_xz_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 68);

    auto ta2_xz_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 69);

    auto ta2_xz_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 70);

    auto ta2_xz_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 75);

    auto ta2_xz_xzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 76);

    auto ta2_xz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 77);

    auto ta2_xz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 78);

    auto ta2_xz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 79);

    auto ta2_xz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 80);

    auto ta2_xz_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 82);

    auto ta2_xz_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 83);

    auto ta2_xz_yzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 84);

    auto ta2_xz_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 85);

    auto ta2_xz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 86);

    auto ta2_xz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 87);

    auto ta2_xz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 88);

    auto ta2_xz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 89);

    auto ta2_yy_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 90);

    auto ta2_yy_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 91);

    auto ta2_yy_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 92);

    auto ta2_yy_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 93);

    auto ta2_yy_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 94);

    auto ta2_yy_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 96);

    auto ta2_yy_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 98);

    auto ta2_yy_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 99);

    auto ta2_yy_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 100);

    auto ta2_yy_xyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 101);

    auto ta2_yy_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 105);

    auto ta2_yy_xzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 106);

    auto ta2_yy_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 107);

    auto ta2_yy_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 108);

    auto ta2_yy_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 109);

    auto ta2_yy_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 110);

    auto ta2_yy_yyz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 111);

    auto ta2_yy_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 112);

    auto ta2_yy_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 113);

    auto ta2_yy_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 115);

    auto ta2_yy_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 116);

    auto ta2_yy_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 117);

    auto ta2_yy_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 118);

    auto ta2_yy_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 119);

    auto ta2_yz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 120);

    auto ta2_yz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 121);

    auto ta2_yz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 122);

    auto ta2_yz_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 123);

    auto ta2_yz_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 124);

    auto ta2_yz_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 126);

    auto ta2_yz_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 128);

    auto ta2_yz_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 129);

    auto ta2_yz_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 130);

    auto ta2_yz_xyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 131);

    auto ta2_yz_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 135);

    auto ta2_yz_xzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 136);

    auto ta2_yz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 137);

    auto ta2_yz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 138);

    auto ta2_yz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 139);

    auto ta2_yz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 140);

    auto ta2_yz_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 142);

    auto ta2_yz_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 143);

    auto ta2_yz_yzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 144);

    auto ta2_yz_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 145);

    auto ta2_yz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 146);

    auto ta2_yz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 147);

    auto ta2_yz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 148);

    auto ta2_yz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 149);

    auto ta2_zz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 150);

    auto ta2_zz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 151);

    auto ta2_zz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 152);

    auto ta2_zz_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 153);

    auto ta2_zz_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 154);

    auto ta2_zz_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 156);

    auto ta2_zz_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 158);

    auto ta2_zz_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 159);

    auto ta2_zz_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 160);

    auto ta2_zz_xyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 161);

    auto ta2_zz_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 165);

    auto ta2_zz_xzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 166);

    auto ta2_zz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 167);

    auto ta2_zz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 168);

    auto ta2_zz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 169);

    auto ta2_zz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 170);

    auto ta2_zz_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 172);

    auto ta2_zz_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 173);

    auto ta2_zz_yzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 174);

    auto ta2_zz_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 175);

    auto ta2_zz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 176);

    auto ta2_zz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 177);

    auto ta2_zz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 178);

    auto ta2_zz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 179);

    // Set up 0-3 components of targeted buffer : GP

    auto ta2_xx_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp);

    auto ta2_xx_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 1);

    auto ta2_xx_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 2);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxx_z_1,   \
                             ta2_xx_xx_x_0,   \
                             ta2_xx_xx_x_1,   \
                             ta2_xx_xx_y_0,   \
                             ta2_xx_xx_y_1,   \
                             ta2_xx_xx_z_0,   \
                             ta2_xx_xx_z_1,   \
                             ta2_xx_xxx_0_0,  \
                             ta2_xx_xxx_0_1,  \
                             ta2_xx_xxx_x_0,  \
                             ta2_xx_xxx_x_1,  \
                             ta2_xx_xxx_y_0,  \
                             ta2_xx_xxx_y_1,  \
                             ta2_xx_xxx_z_0,  \
                             ta2_xx_xxx_z_1,  \
                             ta2_xx_xxxx_x_0, \
                             ta2_xx_xxxx_y_0, \
                             ta2_xx_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxx_x_0[i] = 3.0 * ta2_xx_xx_x_0[i] * fe_0 - 3.0 * ta2_xx_xx_x_1[i] * fe_0 + ta2_xx_xxx_0_0[i] * fe_0 - ta2_xx_xxx_0_1[i] * fe_0 +
                             2.0 * ta1_x_xxx_x_1[i] + ta2_xx_xxx_x_0[i] * pa_x[i] - ta2_xx_xxx_x_1[i] * pc_x[i];

        ta2_xx_xxxx_y_0[i] = 3.0 * ta2_xx_xx_y_0[i] * fe_0 - 3.0 * ta2_xx_xx_y_1[i] * fe_0 + 2.0 * ta1_x_xxx_y_1[i] + ta2_xx_xxx_y_0[i] * pa_x[i] -
                             ta2_xx_xxx_y_1[i] * pc_x[i];

        ta2_xx_xxxx_z_0[i] = 3.0 * ta2_xx_xx_z_0[i] * fe_0 - 3.0 * ta2_xx_xx_z_1[i] * fe_0 + 2.0 * ta1_x_xxx_z_1[i] + ta2_xx_xxx_z_0[i] * pa_x[i] -
                             ta2_xx_xxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto ta2_xx_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 3);

    auto ta2_xx_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 4);

    auto ta2_xx_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 5);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xx_xxx_0_0,  \
                             ta2_xx_xxx_0_1,  \
                             ta2_xx_xxx_x_0,  \
                             ta2_xx_xxx_x_1,  \
                             ta2_xx_xxx_y_0,  \
                             ta2_xx_xxx_y_1,  \
                             ta2_xx_xxx_z_0,  \
                             ta2_xx_xxx_z_1,  \
                             ta2_xx_xxxy_x_0, \
                             ta2_xx_xxxy_y_0, \
                             ta2_xx_xxxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxy_x_0[i] = ta2_xx_xxx_x_0[i] * pa_y[i] - ta2_xx_xxx_x_1[i] * pc_y[i];

        ta2_xx_xxxy_y_0[i] = ta2_xx_xxx_0_0[i] * fe_0 - ta2_xx_xxx_0_1[i] * fe_0 + ta2_xx_xxx_y_0[i] * pa_y[i] - ta2_xx_xxx_y_1[i] * pc_y[i];

        ta2_xx_xxxy_z_0[i] = ta2_xx_xxx_z_0[i] * pa_y[i] - ta2_xx_xxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto ta2_xx_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 6);

    auto ta2_xx_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 7);

    auto ta2_xx_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 8);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xx_xxx_0_0,  \
                             ta2_xx_xxx_0_1,  \
                             ta2_xx_xxx_x_0,  \
                             ta2_xx_xxx_x_1,  \
                             ta2_xx_xxx_y_0,  \
                             ta2_xx_xxx_y_1,  \
                             ta2_xx_xxx_z_0,  \
                             ta2_xx_xxx_z_1,  \
                             ta2_xx_xxxz_x_0, \
                             ta2_xx_xxxz_y_0, \
                             ta2_xx_xxxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxz_x_0[i] = ta2_xx_xxx_x_0[i] * pa_z[i] - ta2_xx_xxx_x_1[i] * pc_z[i];

        ta2_xx_xxxz_y_0[i] = ta2_xx_xxx_y_0[i] * pa_z[i] - ta2_xx_xxx_y_1[i] * pc_z[i];

        ta2_xx_xxxz_z_0[i] = ta2_xx_xxx_0_0[i] * fe_0 - ta2_xx_xxx_0_1[i] * fe_0 + ta2_xx_xxx_z_0[i] * pa_z[i] - ta2_xx_xxx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto ta2_xx_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 9);

    auto ta2_xx_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 10);

    auto ta2_xx_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 11);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xyy_y_1,   \
                             ta2_xx_xx_x_0,   \
                             ta2_xx_xx_x_1,   \
                             ta2_xx_xx_z_0,   \
                             ta2_xx_xx_z_1,   \
                             ta2_xx_xxy_x_0,  \
                             ta2_xx_xxy_x_1,  \
                             ta2_xx_xxy_z_0,  \
                             ta2_xx_xxy_z_1,  \
                             ta2_xx_xxyy_x_0, \
                             ta2_xx_xxyy_y_0, \
                             ta2_xx_xxyy_z_0, \
                             ta2_xx_xyy_y_0,  \
                             ta2_xx_xyy_y_1,  \
                             ta2_xx_yy_y_0,   \
                             ta2_xx_yy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyy_x_0[i] = ta2_xx_xx_x_0[i] * fe_0 - ta2_xx_xx_x_1[i] * fe_0 + ta2_xx_xxy_x_0[i] * pa_y[i] - ta2_xx_xxy_x_1[i] * pc_y[i];

        ta2_xx_xxyy_y_0[i] =
            ta2_xx_yy_y_0[i] * fe_0 - ta2_xx_yy_y_1[i] * fe_0 + 2.0 * ta1_x_xyy_y_1[i] + ta2_xx_xyy_y_0[i] * pa_x[i] - ta2_xx_xyy_y_1[i] * pc_x[i];

        ta2_xx_xxyy_z_0[i] = ta2_xx_xx_z_0[i] * fe_0 - ta2_xx_xx_z_1[i] * fe_0 + ta2_xx_xxy_z_0[i] * pa_y[i] - ta2_xx_xxy_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto ta2_xx_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 12);

    auto ta2_xx_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 13);

    auto ta2_xx_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 14);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta2_xx_xxy_y_0,  \
                             ta2_xx_xxy_y_1,  \
                             ta2_xx_xxyz_x_0, \
                             ta2_xx_xxyz_y_0, \
                             ta2_xx_xxyz_z_0, \
                             ta2_xx_xxz_x_0,  \
                             ta2_xx_xxz_x_1,  \
                             ta2_xx_xxz_z_0,  \
                             ta2_xx_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xxyz_x_0[i] = ta2_xx_xxz_x_0[i] * pa_y[i] - ta2_xx_xxz_x_1[i] * pc_y[i];

        ta2_xx_xxyz_y_0[i] = ta2_xx_xxy_y_0[i] * pa_z[i] - ta2_xx_xxy_y_1[i] * pc_z[i];

        ta2_xx_xxyz_z_0[i] = ta2_xx_xxz_z_0[i] * pa_y[i] - ta2_xx_xxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto ta2_xx_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 15);

    auto ta2_xx_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 16);

    auto ta2_xx_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 17);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xzz_z_1,   \
                             ta2_xx_xx_x_0,   \
                             ta2_xx_xx_x_1,   \
                             ta2_xx_xx_y_0,   \
                             ta2_xx_xx_y_1,   \
                             ta2_xx_xxz_x_0,  \
                             ta2_xx_xxz_x_1,  \
                             ta2_xx_xxz_y_0,  \
                             ta2_xx_xxz_y_1,  \
                             ta2_xx_xxzz_x_0, \
                             ta2_xx_xxzz_y_0, \
                             ta2_xx_xxzz_z_0, \
                             ta2_xx_xzz_z_0,  \
                             ta2_xx_xzz_z_1,  \
                             ta2_xx_zz_z_0,   \
                             ta2_xx_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxzz_x_0[i] = ta2_xx_xx_x_0[i] * fe_0 - ta2_xx_xx_x_1[i] * fe_0 + ta2_xx_xxz_x_0[i] * pa_z[i] - ta2_xx_xxz_x_1[i] * pc_z[i];

        ta2_xx_xxzz_y_0[i] = ta2_xx_xx_y_0[i] * fe_0 - ta2_xx_xx_y_1[i] * fe_0 + ta2_xx_xxz_y_0[i] * pa_z[i] - ta2_xx_xxz_y_1[i] * pc_z[i];

        ta2_xx_xxzz_z_0[i] =
            ta2_xx_zz_z_0[i] * fe_0 - ta2_xx_zz_z_1[i] * fe_0 + 2.0 * ta1_x_xzz_z_1[i] + ta2_xx_xzz_z_0[i] * pa_x[i] - ta2_xx_xzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto ta2_xx_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 18);

    auto ta2_xx_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 19);

    auto ta2_xx_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 20);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_z_1,   \
                             ta2_xx_xy_x_0,   \
                             ta2_xx_xy_x_1,   \
                             ta2_xx_xyy_x_0,  \
                             ta2_xx_xyy_x_1,  \
                             ta2_xx_xyyy_x_0, \
                             ta2_xx_xyyy_y_0, \
                             ta2_xx_xyyy_z_0, \
                             ta2_xx_yyy_y_0,  \
                             ta2_xx_yyy_y_1,  \
                             ta2_xx_yyy_z_0,  \
                             ta2_xx_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyy_x_0[i] =
            2.0 * ta2_xx_xy_x_0[i] * fe_0 - 2.0 * ta2_xx_xy_x_1[i] * fe_0 + ta2_xx_xyy_x_0[i] * pa_y[i] - ta2_xx_xyy_x_1[i] * pc_y[i];

        ta2_xx_xyyy_y_0[i] = 2.0 * ta1_x_yyy_y_1[i] + ta2_xx_yyy_y_0[i] * pa_x[i] - ta2_xx_yyy_y_1[i] * pc_x[i];

        ta2_xx_xyyy_z_0[i] = 2.0 * ta1_x_yyy_z_1[i] + ta2_xx_yyy_z_0[i] * pa_x[i] - ta2_xx_yyy_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto ta2_xx_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 21);

    auto ta2_xx_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 22);

    auto ta2_xx_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 23);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_yyz_z_1,   \
                             ta2_xx_xyy_x_0,  \
                             ta2_xx_xyy_x_1,  \
                             ta2_xx_xyy_y_0,  \
                             ta2_xx_xyy_y_1,  \
                             ta2_xx_xyyz_x_0, \
                             ta2_xx_xyyz_y_0, \
                             ta2_xx_xyyz_z_0, \
                             ta2_xx_yyz_z_0,  \
                             ta2_xx_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xyyz_x_0[i] = ta2_xx_xyy_x_0[i] * pa_z[i] - ta2_xx_xyy_x_1[i] * pc_z[i];

        ta2_xx_xyyz_y_0[i] = ta2_xx_xyy_y_0[i] * pa_z[i] - ta2_xx_xyy_y_1[i] * pc_z[i];

        ta2_xx_xyyz_z_0[i] = 2.0 * ta1_x_yyz_z_1[i] + ta2_xx_yyz_z_0[i] * pa_x[i] - ta2_xx_yyz_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto ta2_xx_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 24);

    auto ta2_xx_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 25);

    auto ta2_xx_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 26);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_yzz_y_1,   \
                             ta2_xx_xyzz_x_0, \
                             ta2_xx_xyzz_y_0, \
                             ta2_xx_xyzz_z_0, \
                             ta2_xx_xzz_x_0,  \
                             ta2_xx_xzz_x_1,  \
                             ta2_xx_xzz_z_0,  \
                             ta2_xx_xzz_z_1,  \
                             ta2_xx_yzz_y_0,  \
                             ta2_xx_yzz_y_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xyzz_x_0[i] = ta2_xx_xzz_x_0[i] * pa_y[i] - ta2_xx_xzz_x_1[i] * pc_y[i];

        ta2_xx_xyzz_y_0[i] = 2.0 * ta1_x_yzz_y_1[i] + ta2_xx_yzz_y_0[i] * pa_x[i] - ta2_xx_yzz_y_1[i] * pc_x[i];

        ta2_xx_xyzz_z_0[i] = ta2_xx_xzz_z_0[i] * pa_y[i] - ta2_xx_xzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto ta2_xx_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 27);

    auto ta2_xx_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 28);

    auto ta2_xx_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_zzz_y_1,   \
                             ta1_x_zzz_z_1,   \
                             ta2_xx_xz_x_0,   \
                             ta2_xx_xz_x_1,   \
                             ta2_xx_xzz_x_0,  \
                             ta2_xx_xzz_x_1,  \
                             ta2_xx_xzzz_x_0, \
                             ta2_xx_xzzz_y_0, \
                             ta2_xx_xzzz_z_0, \
                             ta2_xx_zzz_y_0,  \
                             ta2_xx_zzz_y_1,  \
                             ta2_xx_zzz_z_0,  \
                             ta2_xx_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzzz_x_0[i] =
            2.0 * ta2_xx_xz_x_0[i] * fe_0 - 2.0 * ta2_xx_xz_x_1[i] * fe_0 + ta2_xx_xzz_x_0[i] * pa_z[i] - ta2_xx_xzz_x_1[i] * pc_z[i];

        ta2_xx_xzzz_y_0[i] = 2.0 * ta1_x_zzz_y_1[i] + ta2_xx_zzz_y_0[i] * pa_x[i] - ta2_xx_zzz_y_1[i] * pc_x[i];

        ta2_xx_xzzz_z_0[i] = 2.0 * ta1_x_zzz_z_1[i] + ta2_xx_zzz_z_0[i] * pa_x[i] - ta2_xx_zzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto ta2_xx_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 30);

    auto ta2_xx_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 31);

    auto ta2_xx_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 32);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xx_yy_x_0,   \
                             ta2_xx_yy_x_1,   \
                             ta2_xx_yy_y_0,   \
                             ta2_xx_yy_y_1,   \
                             ta2_xx_yy_z_0,   \
                             ta2_xx_yy_z_1,   \
                             ta2_xx_yyy_0_0,  \
                             ta2_xx_yyy_0_1,  \
                             ta2_xx_yyy_x_0,  \
                             ta2_xx_yyy_x_1,  \
                             ta2_xx_yyy_y_0,  \
                             ta2_xx_yyy_y_1,  \
                             ta2_xx_yyy_z_0,  \
                             ta2_xx_yyy_z_1,  \
                             ta2_xx_yyyy_x_0, \
                             ta2_xx_yyyy_y_0, \
                             ta2_xx_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyy_x_0[i] =
            3.0 * ta2_xx_yy_x_0[i] * fe_0 - 3.0 * ta2_xx_yy_x_1[i] * fe_0 + ta2_xx_yyy_x_0[i] * pa_y[i] - ta2_xx_yyy_x_1[i] * pc_y[i];

        ta2_xx_yyyy_y_0[i] = 3.0 * ta2_xx_yy_y_0[i] * fe_0 - 3.0 * ta2_xx_yy_y_1[i] * fe_0 + ta2_xx_yyy_0_0[i] * fe_0 - ta2_xx_yyy_0_1[i] * fe_0 +
                             ta2_xx_yyy_y_0[i] * pa_y[i] - ta2_xx_yyy_y_1[i] * pc_y[i];

        ta2_xx_yyyy_z_0[i] =
            3.0 * ta2_xx_yy_z_0[i] * fe_0 - 3.0 * ta2_xx_yy_z_1[i] * fe_0 + ta2_xx_yyy_z_0[i] * pa_y[i] - ta2_xx_yyy_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto ta2_xx_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 33);

    auto ta2_xx_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 34);

    auto ta2_xx_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 35);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta2_xx_yyy_x_0,  \
                             ta2_xx_yyy_x_1,  \
                             ta2_xx_yyy_y_0,  \
                             ta2_xx_yyy_y_1,  \
                             ta2_xx_yyyz_x_0, \
                             ta2_xx_yyyz_y_0, \
                             ta2_xx_yyyz_z_0, \
                             ta2_xx_yyz_z_0,  \
                             ta2_xx_yyz_z_1,  \
                             ta2_xx_yz_z_0,   \
                             ta2_xx_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyz_x_0[i] = ta2_xx_yyy_x_0[i] * pa_z[i] - ta2_xx_yyy_x_1[i] * pc_z[i];

        ta2_xx_yyyz_y_0[i] = ta2_xx_yyy_y_0[i] * pa_z[i] - ta2_xx_yyy_y_1[i] * pc_z[i];

        ta2_xx_yyyz_z_0[i] =
            2.0 * ta2_xx_yz_z_0[i] * fe_0 - 2.0 * ta2_xx_yz_z_1[i] * fe_0 + ta2_xx_yyz_z_0[i] * pa_y[i] - ta2_xx_yyz_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto ta2_xx_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 36);

    auto ta2_xx_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 37);

    auto ta2_xx_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 38);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta2_xx_yy_y_0,   \
                             ta2_xx_yy_y_1,   \
                             ta2_xx_yyz_y_0,  \
                             ta2_xx_yyz_y_1,  \
                             ta2_xx_yyzz_x_0, \
                             ta2_xx_yyzz_y_0, \
                             ta2_xx_yyzz_z_0, \
                             ta2_xx_yzz_x_0,  \
                             ta2_xx_yzz_x_1,  \
                             ta2_xx_yzz_z_0,  \
                             ta2_xx_yzz_z_1,  \
                             ta2_xx_zz_x_0,   \
                             ta2_xx_zz_x_1,   \
                             ta2_xx_zz_z_0,   \
                             ta2_xx_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyzz_x_0[i] = ta2_xx_zz_x_0[i] * fe_0 - ta2_xx_zz_x_1[i] * fe_0 + ta2_xx_yzz_x_0[i] * pa_y[i] - ta2_xx_yzz_x_1[i] * pc_y[i];

        ta2_xx_yyzz_y_0[i] = ta2_xx_yy_y_0[i] * fe_0 - ta2_xx_yy_y_1[i] * fe_0 + ta2_xx_yyz_y_0[i] * pa_z[i] - ta2_xx_yyz_y_1[i] * pc_z[i];

        ta2_xx_yyzz_z_0[i] = ta2_xx_zz_z_0[i] * fe_0 - ta2_xx_zz_z_1[i] * fe_0 + ta2_xx_yzz_z_0[i] * pa_y[i] - ta2_xx_yzz_z_1[i] * pc_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto ta2_xx_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 39);

    auto ta2_xx_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 40);

    auto ta2_xx_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 41);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xx_yzzz_x_0, \
                             ta2_xx_yzzz_y_0, \
                             ta2_xx_yzzz_z_0, \
                             ta2_xx_zzz_0_0,  \
                             ta2_xx_zzz_0_1,  \
                             ta2_xx_zzz_x_0,  \
                             ta2_xx_zzz_x_1,  \
                             ta2_xx_zzz_y_0,  \
                             ta2_xx_zzz_y_1,  \
                             ta2_xx_zzz_z_0,  \
                             ta2_xx_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzzz_x_0[i] = ta2_xx_zzz_x_0[i] * pa_y[i] - ta2_xx_zzz_x_1[i] * pc_y[i];

        ta2_xx_yzzz_y_0[i] = ta2_xx_zzz_0_0[i] * fe_0 - ta2_xx_zzz_0_1[i] * fe_0 + ta2_xx_zzz_y_0[i] * pa_y[i] - ta2_xx_zzz_y_1[i] * pc_y[i];

        ta2_xx_yzzz_z_0[i] = ta2_xx_zzz_z_0[i] * pa_y[i] - ta2_xx_zzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto ta2_xx_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 42);

    auto ta2_xx_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 43);

    auto ta2_xx_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 44);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xx_zz_x_0,   \
                             ta2_xx_zz_x_1,   \
                             ta2_xx_zz_y_0,   \
                             ta2_xx_zz_y_1,   \
                             ta2_xx_zz_z_0,   \
                             ta2_xx_zz_z_1,   \
                             ta2_xx_zzz_0_0,  \
                             ta2_xx_zzz_0_1,  \
                             ta2_xx_zzz_x_0,  \
                             ta2_xx_zzz_x_1,  \
                             ta2_xx_zzz_y_0,  \
                             ta2_xx_zzz_y_1,  \
                             ta2_xx_zzz_z_0,  \
                             ta2_xx_zzz_z_1,  \
                             ta2_xx_zzzz_x_0, \
                             ta2_xx_zzzz_y_0, \
                             ta2_xx_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzzz_x_0[i] =
            3.0 * ta2_xx_zz_x_0[i] * fe_0 - 3.0 * ta2_xx_zz_x_1[i] * fe_0 + ta2_xx_zzz_x_0[i] * pa_z[i] - ta2_xx_zzz_x_1[i] * pc_z[i];

        ta2_xx_zzzz_y_0[i] =
            3.0 * ta2_xx_zz_y_0[i] * fe_0 - 3.0 * ta2_xx_zz_y_1[i] * fe_0 + ta2_xx_zzz_y_0[i] * pa_z[i] - ta2_xx_zzz_y_1[i] * pc_z[i];

        ta2_xx_zzzz_z_0[i] = 3.0 * ta2_xx_zz_z_0[i] * fe_0 - 3.0 * ta2_xx_zz_z_1[i] * fe_0 + ta2_xx_zzz_0_0[i] * fe_0 - ta2_xx_zzz_0_1[i] * fe_0 +
                             ta2_xx_zzz_z_0[i] * pa_z[i] - ta2_xx_zzz_z_1[i] * pc_z[i];
    }

    // Set up 45-48 components of targeted buffer : GP

    auto ta2_xy_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 45);

    auto ta2_xy_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 46);

    auto ta2_xy_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 47);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_y_1,   \
                             ta1_y_xxx_z_1,   \
                             ta2_xy_xx_x_0,   \
                             ta2_xy_xx_x_1,   \
                             ta2_xy_xx_y_0,   \
                             ta2_xy_xx_y_1,   \
                             ta2_xy_xx_z_0,   \
                             ta2_xy_xx_z_1,   \
                             ta2_xy_xxx_0_0,  \
                             ta2_xy_xxx_0_1,  \
                             ta2_xy_xxx_x_0,  \
                             ta2_xy_xxx_x_1,  \
                             ta2_xy_xxx_y_0,  \
                             ta2_xy_xxx_y_1,  \
                             ta2_xy_xxx_z_0,  \
                             ta2_xy_xxx_z_1,  \
                             ta2_xy_xxxx_x_0, \
                             ta2_xy_xxxx_y_0, \
                             ta2_xy_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxx_x_0[i] = 3.0 * ta2_xy_xx_x_0[i] * fe_0 - 3.0 * ta2_xy_xx_x_1[i] * fe_0 + ta2_xy_xxx_0_0[i] * fe_0 - ta2_xy_xxx_0_1[i] * fe_0 +
                             ta1_y_xxx_x_1[i] + ta2_xy_xxx_x_0[i] * pa_x[i] - ta2_xy_xxx_x_1[i] * pc_x[i];

        ta2_xy_xxxx_y_0[i] = 3.0 * ta2_xy_xx_y_0[i] * fe_0 - 3.0 * ta2_xy_xx_y_1[i] * fe_0 + ta1_y_xxx_y_1[i] + ta2_xy_xxx_y_0[i] * pa_x[i] -
                             ta2_xy_xxx_y_1[i] * pc_x[i];

        ta2_xy_xxxx_z_0[i] = 3.0 * ta2_xy_xx_z_0[i] * fe_0 - 3.0 * ta2_xy_xx_z_1[i] * fe_0 + ta1_y_xxx_z_1[i] + ta2_xy_xxx_z_0[i] * pa_x[i] -
                             ta2_xy_xxx_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : GP

    auto ta2_xy_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 48);

    auto ta2_xy_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 49);

    auto ta2_xy_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 50);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_z_1,   \
                             ta1_y_xxy_y_1,   \
                             ta2_xy_xxx_x_0,  \
                             ta2_xy_xxx_x_1,  \
                             ta2_xy_xxx_z_0,  \
                             ta2_xy_xxx_z_1,  \
                             ta2_xy_xxxy_x_0, \
                             ta2_xy_xxxy_y_0, \
                             ta2_xy_xxxy_z_0, \
                             ta2_xy_xxy_y_0,  \
                             ta2_xy_xxy_y_1,  \
                             ta2_xy_xy_y_0,   \
                             ta2_xy_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxy_x_0[i] = ta1_x_xxx_x_1[i] + ta2_xy_xxx_x_0[i] * pa_y[i] - ta2_xy_xxx_x_1[i] * pc_y[i];

        ta2_xy_xxxy_y_0[i] = 2.0 * ta2_xy_xy_y_0[i] * fe_0 - 2.0 * ta2_xy_xy_y_1[i] * fe_0 + ta1_y_xxy_y_1[i] + ta2_xy_xxy_y_0[i] * pa_x[i] -
                             ta2_xy_xxy_y_1[i] * pc_x[i];

        ta2_xy_xxxy_z_0[i] = ta1_x_xxx_z_1[i] + ta2_xy_xxx_z_0[i] * pa_y[i] - ta2_xy_xxx_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : GP

    auto ta2_xy_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 51);

    auto ta2_xy_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 52);

    auto ta2_xy_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 53);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xy_xxx_0_0,  \
                             ta2_xy_xxx_0_1,  \
                             ta2_xy_xxx_x_0,  \
                             ta2_xy_xxx_x_1,  \
                             ta2_xy_xxx_y_0,  \
                             ta2_xy_xxx_y_1,  \
                             ta2_xy_xxx_z_0,  \
                             ta2_xy_xxx_z_1,  \
                             ta2_xy_xxxz_x_0, \
                             ta2_xy_xxxz_y_0, \
                             ta2_xy_xxxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxz_x_0[i] = ta2_xy_xxx_x_0[i] * pa_z[i] - ta2_xy_xxx_x_1[i] * pc_z[i];

        ta2_xy_xxxz_y_0[i] = ta2_xy_xxx_y_0[i] * pa_z[i] - ta2_xy_xxx_y_1[i] * pc_z[i];

        ta2_xy_xxxz_z_0[i] = ta2_xy_xxx_0_0[i] * fe_0 - ta2_xy_xxx_0_1[i] * fe_0 + ta2_xy_xxx_z_0[i] * pa_z[i] - ta2_xy_xxx_z_1[i] * pc_z[i];
    }

    // Set up 54-57 components of targeted buffer : GP

    auto ta2_xy_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 54);

    auto ta2_xy_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 55);

    auto ta2_xy_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 56);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xxy_x_1,   \
                             ta1_y_xyy_y_1,   \
                             ta1_y_xyy_z_1,   \
                             ta2_xy_xx_x_0,   \
                             ta2_xy_xx_x_1,   \
                             ta2_xy_xxy_x_0,  \
                             ta2_xy_xxy_x_1,  \
                             ta2_xy_xxyy_x_0, \
                             ta2_xy_xxyy_y_0, \
                             ta2_xy_xxyy_z_0, \
                             ta2_xy_xyy_y_0,  \
                             ta2_xy_xyy_y_1,  \
                             ta2_xy_xyy_z_0,  \
                             ta2_xy_xyy_z_1,  \
                             ta2_xy_yy_y_0,   \
                             ta2_xy_yy_y_1,   \
                             ta2_xy_yy_z_0,   \
                             ta2_xy_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyy_x_0[i] =
            ta2_xy_xx_x_0[i] * fe_0 - ta2_xy_xx_x_1[i] * fe_0 + ta1_x_xxy_x_1[i] + ta2_xy_xxy_x_0[i] * pa_y[i] - ta2_xy_xxy_x_1[i] * pc_y[i];

        ta2_xy_xxyy_y_0[i] =
            ta2_xy_yy_y_0[i] * fe_0 - ta2_xy_yy_y_1[i] * fe_0 + ta1_y_xyy_y_1[i] + ta2_xy_xyy_y_0[i] * pa_x[i] - ta2_xy_xyy_y_1[i] * pc_x[i];

        ta2_xy_xxyy_z_0[i] =
            ta2_xy_yy_z_0[i] * fe_0 - ta2_xy_yy_z_1[i] * fe_0 + ta1_y_xyy_z_1[i] + ta2_xy_xyy_z_0[i] * pa_x[i] - ta2_xy_xyy_z_1[i] * pc_x[i];
    }

    // Set up 57-60 components of targeted buffer : GP

    auto ta2_xy_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 57);

    auto ta2_xy_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 58);

    auto ta2_xy_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 59);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxz_z_1,   \
                             ta2_xy_xxy_x_0,  \
                             ta2_xy_xxy_x_1,  \
                             ta2_xy_xxy_y_0,  \
                             ta2_xy_xxy_y_1,  \
                             ta2_xy_xxyz_x_0, \
                             ta2_xy_xxyz_y_0, \
                             ta2_xy_xxyz_z_0, \
                             ta2_xy_xxz_z_0,  \
                             ta2_xy_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xxyz_x_0[i] = ta2_xy_xxy_x_0[i] * pa_z[i] - ta2_xy_xxy_x_1[i] * pc_z[i];

        ta2_xy_xxyz_y_0[i] = ta2_xy_xxy_y_0[i] * pa_z[i] - ta2_xy_xxy_y_1[i] * pc_z[i];

        ta2_xy_xxyz_z_0[i] = ta1_x_xxz_z_1[i] + ta2_xy_xxz_z_0[i] * pa_y[i] - ta2_xy_xxz_z_1[i] * pc_y[i];
    }

    // Set up 60-63 components of targeted buffer : GP

    auto ta2_xy_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 60);

    auto ta2_xy_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 61);

    auto ta2_xy_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 62);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xzz_z_1,   \
                             ta2_xy_xx_x_0,   \
                             ta2_xy_xx_x_1,   \
                             ta2_xy_xx_y_0,   \
                             ta2_xy_xx_y_1,   \
                             ta2_xy_xxz_x_0,  \
                             ta2_xy_xxz_x_1,  \
                             ta2_xy_xxz_y_0,  \
                             ta2_xy_xxz_y_1,  \
                             ta2_xy_xxzz_x_0, \
                             ta2_xy_xxzz_y_0, \
                             ta2_xy_xxzz_z_0, \
                             ta2_xy_xzz_z_0,  \
                             ta2_xy_xzz_z_1,  \
                             ta2_xy_zz_z_0,   \
                             ta2_xy_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxzz_x_0[i] = ta2_xy_xx_x_0[i] * fe_0 - ta2_xy_xx_x_1[i] * fe_0 + ta2_xy_xxz_x_0[i] * pa_z[i] - ta2_xy_xxz_x_1[i] * pc_z[i];

        ta2_xy_xxzz_y_0[i] = ta2_xy_xx_y_0[i] * fe_0 - ta2_xy_xx_y_1[i] * fe_0 + ta2_xy_xxz_y_0[i] * pa_z[i] - ta2_xy_xxz_y_1[i] * pc_z[i];

        ta2_xy_xxzz_z_0[i] =
            ta2_xy_zz_z_0[i] * fe_0 - ta2_xy_zz_z_1[i] * fe_0 + ta1_y_xzz_z_1[i] + ta2_xy_xzz_z_0[i] * pa_x[i] - ta2_xy_xzz_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : GP

    auto ta2_xy_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 63);

    auto ta2_xy_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 64);

    auto ta2_xy_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 65);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_z_1,   \
                             ta2_xy_xyyy_x_0, \
                             ta2_xy_xyyy_y_0, \
                             ta2_xy_xyyy_z_0, \
                             ta2_xy_yyy_0_0,  \
                             ta2_xy_yyy_0_1,  \
                             ta2_xy_yyy_x_0,  \
                             ta2_xy_yyy_x_1,  \
                             ta2_xy_yyy_y_0,  \
                             ta2_xy_yyy_y_1,  \
                             ta2_xy_yyy_z_0,  \
                             ta2_xy_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyy_x_0[i] =
            ta2_xy_yyy_0_0[i] * fe_0 - ta2_xy_yyy_0_1[i] * fe_0 + ta1_y_yyy_x_1[i] + ta2_xy_yyy_x_0[i] * pa_x[i] - ta2_xy_yyy_x_1[i] * pc_x[i];

        ta2_xy_xyyy_y_0[i] = ta1_y_yyy_y_1[i] + ta2_xy_yyy_y_0[i] * pa_x[i] - ta2_xy_yyy_y_1[i] * pc_x[i];

        ta2_xy_xyyy_z_0[i] = ta1_y_yyy_z_1[i] + ta2_xy_yyy_z_0[i] * pa_x[i] - ta2_xy_yyy_z_1[i] * pc_x[i];
    }

    // Set up 66-69 components of targeted buffer : GP

    auto ta2_xy_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 66);

    auto ta2_xy_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 67);

    auto ta2_xy_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 68);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_yyz_z_1,   \
                             ta2_xy_xyy_x_0,  \
                             ta2_xy_xyy_x_1,  \
                             ta2_xy_xyy_y_0,  \
                             ta2_xy_xyy_y_1,  \
                             ta2_xy_xyyz_x_0, \
                             ta2_xy_xyyz_y_0, \
                             ta2_xy_xyyz_z_0, \
                             ta2_xy_yyz_z_0,  \
                             ta2_xy_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xyyz_x_0[i] = ta2_xy_xyy_x_0[i] * pa_z[i] - ta2_xy_xyy_x_1[i] * pc_z[i];

        ta2_xy_xyyz_y_0[i] = ta2_xy_xyy_y_0[i] * pa_z[i] - ta2_xy_xyy_y_1[i] * pc_z[i];

        ta2_xy_xyyz_z_0[i] = ta1_y_yyz_z_1[i] + ta2_xy_yyz_z_0[i] * pa_x[i] - ta2_xy_yyz_z_1[i] * pc_x[i];
    }

    // Set up 69-72 components of targeted buffer : GP

    auto ta2_xy_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 69);

    auto ta2_xy_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 70);

    auto ta2_xy_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 71);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xzz_x_1,   \
                             ta1_y_yzz_y_1,   \
                             ta1_y_yzz_z_1,   \
                             ta2_xy_xyzz_x_0, \
                             ta2_xy_xyzz_y_0, \
                             ta2_xy_xyzz_z_0, \
                             ta2_xy_xzz_x_0,  \
                             ta2_xy_xzz_x_1,  \
                             ta2_xy_yzz_y_0,  \
                             ta2_xy_yzz_y_1,  \
                             ta2_xy_yzz_z_0,  \
                             ta2_xy_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xyzz_x_0[i] = ta1_x_xzz_x_1[i] + ta2_xy_xzz_x_0[i] * pa_y[i] - ta2_xy_xzz_x_1[i] * pc_y[i];

        ta2_xy_xyzz_y_0[i] = ta1_y_yzz_y_1[i] + ta2_xy_yzz_y_0[i] * pa_x[i] - ta2_xy_yzz_y_1[i] * pc_x[i];

        ta2_xy_xyzz_z_0[i] = ta1_y_yzz_z_1[i] + ta2_xy_yzz_z_0[i] * pa_x[i] - ta2_xy_yzz_z_1[i] * pc_x[i];
    }

    // Set up 72-75 components of targeted buffer : GP

    auto ta2_xy_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 72);

    auto ta2_xy_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 73);

    auto ta2_xy_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 74);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_z_1,   \
                             ta2_xy_xz_x_0,   \
                             ta2_xy_xz_x_1,   \
                             ta2_xy_xzz_x_0,  \
                             ta2_xy_xzz_x_1,  \
                             ta2_xy_xzzz_x_0, \
                             ta2_xy_xzzz_y_0, \
                             ta2_xy_xzzz_z_0, \
                             ta2_xy_zzz_y_0,  \
                             ta2_xy_zzz_y_1,  \
                             ta2_xy_zzz_z_0,  \
                             ta2_xy_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzzz_x_0[i] =
            2.0 * ta2_xy_xz_x_0[i] * fe_0 - 2.0 * ta2_xy_xz_x_1[i] * fe_0 + ta2_xy_xzz_x_0[i] * pa_z[i] - ta2_xy_xzz_x_1[i] * pc_z[i];

        ta2_xy_xzzz_y_0[i] = ta1_y_zzz_y_1[i] + ta2_xy_zzz_y_0[i] * pa_x[i] - ta2_xy_zzz_y_1[i] * pc_x[i];

        ta2_xy_xzzz_z_0[i] = ta1_y_zzz_z_1[i] + ta2_xy_zzz_z_0[i] * pa_x[i] - ta2_xy_zzz_z_1[i] * pc_x[i];
    }

    // Set up 75-78 components of targeted buffer : GP

    auto ta2_xy_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 75);

    auto ta2_xy_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 76);

    auto ta2_xy_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 77);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yyy_x_1,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_z_1,   \
                             ta2_xy_yy_x_0,   \
                             ta2_xy_yy_x_1,   \
                             ta2_xy_yy_y_0,   \
                             ta2_xy_yy_y_1,   \
                             ta2_xy_yy_z_0,   \
                             ta2_xy_yy_z_1,   \
                             ta2_xy_yyy_0_0,  \
                             ta2_xy_yyy_0_1,  \
                             ta2_xy_yyy_x_0,  \
                             ta2_xy_yyy_x_1,  \
                             ta2_xy_yyy_y_0,  \
                             ta2_xy_yyy_y_1,  \
                             ta2_xy_yyy_z_0,  \
                             ta2_xy_yyy_z_1,  \
                             ta2_xy_yyyy_x_0, \
                             ta2_xy_yyyy_y_0, \
                             ta2_xy_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyy_x_0[i] = 3.0 * ta2_xy_yy_x_0[i] * fe_0 - 3.0 * ta2_xy_yy_x_1[i] * fe_0 + ta1_x_yyy_x_1[i] + ta2_xy_yyy_x_0[i] * pa_y[i] -
                             ta2_xy_yyy_x_1[i] * pc_y[i];

        ta2_xy_yyyy_y_0[i] = 3.0 * ta2_xy_yy_y_0[i] * fe_0 - 3.0 * ta2_xy_yy_y_1[i] * fe_0 + ta2_xy_yyy_0_0[i] * fe_0 - ta2_xy_yyy_0_1[i] * fe_0 +
                             ta1_x_yyy_y_1[i] + ta2_xy_yyy_y_0[i] * pa_y[i] - ta2_xy_yyy_y_1[i] * pc_y[i];

        ta2_xy_yyyy_z_0[i] = 3.0 * ta2_xy_yy_z_0[i] * fe_0 - 3.0 * ta2_xy_yy_z_1[i] * fe_0 + ta1_x_yyy_z_1[i] + ta2_xy_yyy_z_0[i] * pa_y[i] -
                             ta2_xy_yyy_z_1[i] * pc_y[i];
    }

    // Set up 78-81 components of targeted buffer : GP

    auto ta2_xy_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 78);

    auto ta2_xy_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 79);

    auto ta2_xy_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 80);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xy_yyy_0_0,  \
                             ta2_xy_yyy_0_1,  \
                             ta2_xy_yyy_x_0,  \
                             ta2_xy_yyy_x_1,  \
                             ta2_xy_yyy_y_0,  \
                             ta2_xy_yyy_y_1,  \
                             ta2_xy_yyy_z_0,  \
                             ta2_xy_yyy_z_1,  \
                             ta2_xy_yyyz_x_0, \
                             ta2_xy_yyyz_y_0, \
                             ta2_xy_yyyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyz_x_0[i] = ta2_xy_yyy_x_0[i] * pa_z[i] - ta2_xy_yyy_x_1[i] * pc_z[i];

        ta2_xy_yyyz_y_0[i] = ta2_xy_yyy_y_0[i] * pa_z[i] - ta2_xy_yyy_y_1[i] * pc_z[i];

        ta2_xy_yyyz_z_0[i] = ta2_xy_yyy_0_0[i] * fe_0 - ta2_xy_yyy_0_1[i] * fe_0 + ta2_xy_yyy_z_0[i] * pa_z[i] - ta2_xy_yyy_z_1[i] * pc_z[i];
    }

    // Set up 81-84 components of targeted buffer : GP

    auto ta2_xy_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 81);

    auto ta2_xy_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 82);

    auto ta2_xy_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 83);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yzz_z_1,   \
                             ta2_xy_yy_x_0,   \
                             ta2_xy_yy_x_1,   \
                             ta2_xy_yy_y_0,   \
                             ta2_xy_yy_y_1,   \
                             ta2_xy_yyz_x_0,  \
                             ta2_xy_yyz_x_1,  \
                             ta2_xy_yyz_y_0,  \
                             ta2_xy_yyz_y_1,  \
                             ta2_xy_yyzz_x_0, \
                             ta2_xy_yyzz_y_0, \
                             ta2_xy_yyzz_z_0, \
                             ta2_xy_yzz_z_0,  \
                             ta2_xy_yzz_z_1,  \
                             ta2_xy_zz_z_0,   \
                             ta2_xy_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyzz_x_0[i] = ta2_xy_yy_x_0[i] * fe_0 - ta2_xy_yy_x_1[i] * fe_0 + ta2_xy_yyz_x_0[i] * pa_z[i] - ta2_xy_yyz_x_1[i] * pc_z[i];

        ta2_xy_yyzz_y_0[i] = ta2_xy_yy_y_0[i] * fe_0 - ta2_xy_yy_y_1[i] * fe_0 + ta2_xy_yyz_y_0[i] * pa_z[i] - ta2_xy_yyz_y_1[i] * pc_z[i];

        ta2_xy_yyzz_z_0[i] =
            ta2_xy_zz_z_0[i] * fe_0 - ta2_xy_zz_z_1[i] * fe_0 + ta1_x_yzz_z_1[i] + ta2_xy_yzz_z_0[i] * pa_y[i] - ta2_xy_yzz_z_1[i] * pc_y[i];
    }

    // Set up 84-87 components of targeted buffer : GP

    auto ta2_xy_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 84);

    auto ta2_xy_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 85);

    auto ta2_xy_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 86);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_z_1,   \
                             ta2_xy_yz_y_0,   \
                             ta2_xy_yz_y_1,   \
                             ta2_xy_yzz_y_0,  \
                             ta2_xy_yzz_y_1,  \
                             ta2_xy_yzzz_x_0, \
                             ta2_xy_yzzz_y_0, \
                             ta2_xy_yzzz_z_0, \
                             ta2_xy_zzz_x_0,  \
                             ta2_xy_zzz_x_1,  \
                             ta2_xy_zzz_z_0,  \
                             ta2_xy_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzzz_x_0[i] = ta1_x_zzz_x_1[i] + ta2_xy_zzz_x_0[i] * pa_y[i] - ta2_xy_zzz_x_1[i] * pc_y[i];

        ta2_xy_yzzz_y_0[i] =
            2.0 * ta2_xy_yz_y_0[i] * fe_0 - 2.0 * ta2_xy_yz_y_1[i] * fe_0 + ta2_xy_yzz_y_0[i] * pa_z[i] - ta2_xy_yzz_y_1[i] * pc_z[i];

        ta2_xy_yzzz_z_0[i] = ta1_x_zzz_z_1[i] + ta2_xy_zzz_z_0[i] * pa_y[i] - ta2_xy_zzz_z_1[i] * pc_y[i];
    }

    // Set up 87-90 components of targeted buffer : GP

    auto ta2_xy_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 87);

    auto ta2_xy_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 88);

    auto ta2_xy_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 89);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xy_zz_x_0,   \
                             ta2_xy_zz_x_1,   \
                             ta2_xy_zz_y_0,   \
                             ta2_xy_zz_y_1,   \
                             ta2_xy_zz_z_0,   \
                             ta2_xy_zz_z_1,   \
                             ta2_xy_zzz_0_0,  \
                             ta2_xy_zzz_0_1,  \
                             ta2_xy_zzz_x_0,  \
                             ta2_xy_zzz_x_1,  \
                             ta2_xy_zzz_y_0,  \
                             ta2_xy_zzz_y_1,  \
                             ta2_xy_zzz_z_0,  \
                             ta2_xy_zzz_z_1,  \
                             ta2_xy_zzzz_x_0, \
                             ta2_xy_zzzz_y_0, \
                             ta2_xy_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzzz_x_0[i] =
            3.0 * ta2_xy_zz_x_0[i] * fe_0 - 3.0 * ta2_xy_zz_x_1[i] * fe_0 + ta2_xy_zzz_x_0[i] * pa_z[i] - ta2_xy_zzz_x_1[i] * pc_z[i];

        ta2_xy_zzzz_y_0[i] =
            3.0 * ta2_xy_zz_y_0[i] * fe_0 - 3.0 * ta2_xy_zz_y_1[i] * fe_0 + ta2_xy_zzz_y_0[i] * pa_z[i] - ta2_xy_zzz_y_1[i] * pc_z[i];

        ta2_xy_zzzz_z_0[i] = 3.0 * ta2_xy_zz_z_0[i] * fe_0 - 3.0 * ta2_xy_zz_z_1[i] * fe_0 + ta2_xy_zzz_0_0[i] * fe_0 - ta2_xy_zzz_0_1[i] * fe_0 +
                             ta2_xy_zzz_z_0[i] * pa_z[i] - ta2_xy_zzz_z_1[i] * pc_z[i];
    }

    // Set up 90-93 components of targeted buffer : GP

    auto ta2_xz_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 90);

    auto ta2_xz_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 91);

    auto ta2_xz_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 92);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_y_1,   \
                             ta1_z_xxx_z_1,   \
                             ta2_xz_xx_x_0,   \
                             ta2_xz_xx_x_1,   \
                             ta2_xz_xx_y_0,   \
                             ta2_xz_xx_y_1,   \
                             ta2_xz_xx_z_0,   \
                             ta2_xz_xx_z_1,   \
                             ta2_xz_xxx_0_0,  \
                             ta2_xz_xxx_0_1,  \
                             ta2_xz_xxx_x_0,  \
                             ta2_xz_xxx_x_1,  \
                             ta2_xz_xxx_y_0,  \
                             ta2_xz_xxx_y_1,  \
                             ta2_xz_xxx_z_0,  \
                             ta2_xz_xxx_z_1,  \
                             ta2_xz_xxxx_x_0, \
                             ta2_xz_xxxx_y_0, \
                             ta2_xz_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxx_x_0[i] = 3.0 * ta2_xz_xx_x_0[i] * fe_0 - 3.0 * ta2_xz_xx_x_1[i] * fe_0 + ta2_xz_xxx_0_0[i] * fe_0 - ta2_xz_xxx_0_1[i] * fe_0 +
                             ta1_z_xxx_x_1[i] + ta2_xz_xxx_x_0[i] * pa_x[i] - ta2_xz_xxx_x_1[i] * pc_x[i];

        ta2_xz_xxxx_y_0[i] = 3.0 * ta2_xz_xx_y_0[i] * fe_0 - 3.0 * ta2_xz_xx_y_1[i] * fe_0 + ta1_z_xxx_y_1[i] + ta2_xz_xxx_y_0[i] * pa_x[i] -
                             ta2_xz_xxx_y_1[i] * pc_x[i];

        ta2_xz_xxxx_z_0[i] = 3.0 * ta2_xz_xx_z_0[i] * fe_0 - 3.0 * ta2_xz_xx_z_1[i] * fe_0 + ta1_z_xxx_z_1[i] + ta2_xz_xxx_z_0[i] * pa_x[i] -
                             ta2_xz_xxx_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : GP

    auto ta2_xz_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 93);

    auto ta2_xz_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 94);

    auto ta2_xz_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 95);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xz_xxx_0_0,  \
                             ta2_xz_xxx_0_1,  \
                             ta2_xz_xxx_x_0,  \
                             ta2_xz_xxx_x_1,  \
                             ta2_xz_xxx_y_0,  \
                             ta2_xz_xxx_y_1,  \
                             ta2_xz_xxx_z_0,  \
                             ta2_xz_xxx_z_1,  \
                             ta2_xz_xxxy_x_0, \
                             ta2_xz_xxxy_y_0, \
                             ta2_xz_xxxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxy_x_0[i] = ta2_xz_xxx_x_0[i] * pa_y[i] - ta2_xz_xxx_x_1[i] * pc_y[i];

        ta2_xz_xxxy_y_0[i] = ta2_xz_xxx_0_0[i] * fe_0 - ta2_xz_xxx_0_1[i] * fe_0 + ta2_xz_xxx_y_0[i] * pa_y[i] - ta2_xz_xxx_y_1[i] * pc_y[i];

        ta2_xz_xxxy_z_0[i] = ta2_xz_xxx_z_0[i] * pa_y[i] - ta2_xz_xxx_z_1[i] * pc_y[i];
    }

    // Set up 96-99 components of targeted buffer : GP

    auto ta2_xz_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 96);

    auto ta2_xz_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 97);

    auto ta2_xz_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 98);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_y_1,   \
                             ta1_z_xxz_z_1,   \
                             ta2_xz_xxx_x_0,  \
                             ta2_xz_xxx_x_1,  \
                             ta2_xz_xxx_y_0,  \
                             ta2_xz_xxx_y_1,  \
                             ta2_xz_xxxz_x_0, \
                             ta2_xz_xxxz_y_0, \
                             ta2_xz_xxxz_z_0, \
                             ta2_xz_xxz_z_0,  \
                             ta2_xz_xxz_z_1,  \
                             ta2_xz_xz_z_0,   \
                             ta2_xz_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxz_x_0[i] = ta1_x_xxx_x_1[i] + ta2_xz_xxx_x_0[i] * pa_z[i] - ta2_xz_xxx_x_1[i] * pc_z[i];

        ta2_xz_xxxz_y_0[i] = ta1_x_xxx_y_1[i] + ta2_xz_xxx_y_0[i] * pa_z[i] - ta2_xz_xxx_y_1[i] * pc_z[i];

        ta2_xz_xxxz_z_0[i] = 2.0 * ta2_xz_xz_z_0[i] * fe_0 - 2.0 * ta2_xz_xz_z_1[i] * fe_0 + ta1_z_xxz_z_1[i] + ta2_xz_xxz_z_0[i] * pa_x[i] -
                             ta2_xz_xxz_z_1[i] * pc_x[i];
    }

    // Set up 99-102 components of targeted buffer : GP

    auto ta2_xz_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 99);

    auto ta2_xz_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 100);

    auto ta2_xz_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 101);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xyy_y_1,   \
                             ta2_xz_xx_x_0,   \
                             ta2_xz_xx_x_1,   \
                             ta2_xz_xx_z_0,   \
                             ta2_xz_xx_z_1,   \
                             ta2_xz_xxy_x_0,  \
                             ta2_xz_xxy_x_1,  \
                             ta2_xz_xxy_z_0,  \
                             ta2_xz_xxy_z_1,  \
                             ta2_xz_xxyy_x_0, \
                             ta2_xz_xxyy_y_0, \
                             ta2_xz_xxyy_z_0, \
                             ta2_xz_xyy_y_0,  \
                             ta2_xz_xyy_y_1,  \
                             ta2_xz_yy_y_0,   \
                             ta2_xz_yy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyy_x_0[i] = ta2_xz_xx_x_0[i] * fe_0 - ta2_xz_xx_x_1[i] * fe_0 + ta2_xz_xxy_x_0[i] * pa_y[i] - ta2_xz_xxy_x_1[i] * pc_y[i];

        ta2_xz_xxyy_y_0[i] =
            ta2_xz_yy_y_0[i] * fe_0 - ta2_xz_yy_y_1[i] * fe_0 + ta1_z_xyy_y_1[i] + ta2_xz_xyy_y_0[i] * pa_x[i] - ta2_xz_xyy_y_1[i] * pc_x[i];

        ta2_xz_xxyy_z_0[i] = ta2_xz_xx_z_0[i] * fe_0 - ta2_xz_xx_z_1[i] * fe_0 + ta2_xz_xxy_z_0[i] * pa_y[i] - ta2_xz_xxy_z_1[i] * pc_y[i];
    }

    // Set up 102-105 components of targeted buffer : GP

    auto ta2_xz_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 102);

    auto ta2_xz_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 103);

    auto ta2_xz_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 104);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxy_y_1,   \
                             ta2_xz_xxy_y_0,  \
                             ta2_xz_xxy_y_1,  \
                             ta2_xz_xxyz_x_0, \
                             ta2_xz_xxyz_y_0, \
                             ta2_xz_xxyz_z_0, \
                             ta2_xz_xxz_x_0,  \
                             ta2_xz_xxz_x_1,  \
                             ta2_xz_xxz_z_0,  \
                             ta2_xz_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xxyz_x_0[i] = ta2_xz_xxz_x_0[i] * pa_y[i] - ta2_xz_xxz_x_1[i] * pc_y[i];

        ta2_xz_xxyz_y_0[i] = ta1_x_xxy_y_1[i] + ta2_xz_xxy_y_0[i] * pa_z[i] - ta2_xz_xxy_y_1[i] * pc_z[i];

        ta2_xz_xxyz_z_0[i] = ta2_xz_xxz_z_0[i] * pa_y[i] - ta2_xz_xxz_z_1[i] * pc_y[i];
    }

    // Set up 105-108 components of targeted buffer : GP

    auto ta2_xz_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 105);

    auto ta2_xz_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 106);

    auto ta2_xz_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 107);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xxz_x_1,   \
                             ta1_z_xzz_y_1,   \
                             ta1_z_xzz_z_1,   \
                             ta2_xz_xx_x_0,   \
                             ta2_xz_xx_x_1,   \
                             ta2_xz_xxz_x_0,  \
                             ta2_xz_xxz_x_1,  \
                             ta2_xz_xxzz_x_0, \
                             ta2_xz_xxzz_y_0, \
                             ta2_xz_xxzz_z_0, \
                             ta2_xz_xzz_y_0,  \
                             ta2_xz_xzz_y_1,  \
                             ta2_xz_xzz_z_0,  \
                             ta2_xz_xzz_z_1,  \
                             ta2_xz_zz_y_0,   \
                             ta2_xz_zz_y_1,   \
                             ta2_xz_zz_z_0,   \
                             ta2_xz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxzz_x_0[i] =
            ta2_xz_xx_x_0[i] * fe_0 - ta2_xz_xx_x_1[i] * fe_0 + ta1_x_xxz_x_1[i] + ta2_xz_xxz_x_0[i] * pa_z[i] - ta2_xz_xxz_x_1[i] * pc_z[i];

        ta2_xz_xxzz_y_0[i] =
            ta2_xz_zz_y_0[i] * fe_0 - ta2_xz_zz_y_1[i] * fe_0 + ta1_z_xzz_y_1[i] + ta2_xz_xzz_y_0[i] * pa_x[i] - ta2_xz_xzz_y_1[i] * pc_x[i];

        ta2_xz_xxzz_z_0[i] =
            ta2_xz_zz_z_0[i] * fe_0 - ta2_xz_zz_z_1[i] * fe_0 + ta1_z_xzz_z_1[i] + ta2_xz_xzz_z_0[i] * pa_x[i] - ta2_xz_xzz_z_1[i] * pc_x[i];
    }

    // Set up 108-111 components of targeted buffer : GP

    auto ta2_xz_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 108);

    auto ta2_xz_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 109);

    auto ta2_xz_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 110);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_z_1,   \
                             ta2_xz_xy_x_0,   \
                             ta2_xz_xy_x_1,   \
                             ta2_xz_xyy_x_0,  \
                             ta2_xz_xyy_x_1,  \
                             ta2_xz_xyyy_x_0, \
                             ta2_xz_xyyy_y_0, \
                             ta2_xz_xyyy_z_0, \
                             ta2_xz_yyy_y_0,  \
                             ta2_xz_yyy_y_1,  \
                             ta2_xz_yyy_z_0,  \
                             ta2_xz_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyy_x_0[i] =
            2.0 * ta2_xz_xy_x_0[i] * fe_0 - 2.0 * ta2_xz_xy_x_1[i] * fe_0 + ta2_xz_xyy_x_0[i] * pa_y[i] - ta2_xz_xyy_x_1[i] * pc_y[i];

        ta2_xz_xyyy_y_0[i] = ta1_z_yyy_y_1[i] + ta2_xz_yyy_y_0[i] * pa_x[i] - ta2_xz_yyy_y_1[i] * pc_x[i];

        ta2_xz_xyyy_z_0[i] = ta1_z_yyy_z_1[i] + ta2_xz_yyy_z_0[i] * pa_x[i] - ta2_xz_yyy_z_1[i] * pc_x[i];
    }

    // Set up 111-114 components of targeted buffer : GP

    auto ta2_xz_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 111);

    auto ta2_xz_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 112);

    auto ta2_xz_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 113);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xyy_x_1,   \
                             ta1_z_yyz_y_1,   \
                             ta1_z_yyz_z_1,   \
                             ta2_xz_xyy_x_0,  \
                             ta2_xz_xyy_x_1,  \
                             ta2_xz_xyyz_x_0, \
                             ta2_xz_xyyz_y_0, \
                             ta2_xz_xyyz_z_0, \
                             ta2_xz_yyz_y_0,  \
                             ta2_xz_yyz_y_1,  \
                             ta2_xz_yyz_z_0,  \
                             ta2_xz_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xyyz_x_0[i] = ta1_x_xyy_x_1[i] + ta2_xz_xyy_x_0[i] * pa_z[i] - ta2_xz_xyy_x_1[i] * pc_z[i];

        ta2_xz_xyyz_y_0[i] = ta1_z_yyz_y_1[i] + ta2_xz_yyz_y_0[i] * pa_x[i] - ta2_xz_yyz_y_1[i] * pc_x[i];

        ta2_xz_xyyz_z_0[i] = ta1_z_yyz_z_1[i] + ta2_xz_yyz_z_0[i] * pa_x[i] - ta2_xz_yyz_z_1[i] * pc_x[i];
    }

    // Set up 114-117 components of targeted buffer : GP

    auto ta2_xz_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 114);

    auto ta2_xz_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 115);

    auto ta2_xz_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 116);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_yzz_y_1,   \
                             ta2_xz_xyzz_x_0, \
                             ta2_xz_xyzz_y_0, \
                             ta2_xz_xyzz_z_0, \
                             ta2_xz_xzz_x_0,  \
                             ta2_xz_xzz_x_1,  \
                             ta2_xz_xzz_z_0,  \
                             ta2_xz_xzz_z_1,  \
                             ta2_xz_yzz_y_0,  \
                             ta2_xz_yzz_y_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xyzz_x_0[i] = ta2_xz_xzz_x_0[i] * pa_y[i] - ta2_xz_xzz_x_1[i] * pc_y[i];

        ta2_xz_xyzz_y_0[i] = ta1_z_yzz_y_1[i] + ta2_xz_yzz_y_0[i] * pa_x[i] - ta2_xz_yzz_y_1[i] * pc_x[i];

        ta2_xz_xyzz_z_0[i] = ta2_xz_xzz_z_0[i] * pa_y[i] - ta2_xz_xzz_z_1[i] * pc_y[i];
    }

    // Set up 117-120 components of targeted buffer : GP

    auto ta2_xz_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 117);

    auto ta2_xz_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 118);

    auto ta2_xz_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 119);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_z_1,   \
                             ta2_xz_xzzz_x_0, \
                             ta2_xz_xzzz_y_0, \
                             ta2_xz_xzzz_z_0, \
                             ta2_xz_zzz_0_0,  \
                             ta2_xz_zzz_0_1,  \
                             ta2_xz_zzz_x_0,  \
                             ta2_xz_zzz_x_1,  \
                             ta2_xz_zzz_y_0,  \
                             ta2_xz_zzz_y_1,  \
                             ta2_xz_zzz_z_0,  \
                             ta2_xz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzzz_x_0[i] =
            ta2_xz_zzz_0_0[i] * fe_0 - ta2_xz_zzz_0_1[i] * fe_0 + ta1_z_zzz_x_1[i] + ta2_xz_zzz_x_0[i] * pa_x[i] - ta2_xz_zzz_x_1[i] * pc_x[i];

        ta2_xz_xzzz_y_0[i] = ta1_z_zzz_y_1[i] + ta2_xz_zzz_y_0[i] * pa_x[i] - ta2_xz_zzz_y_1[i] * pc_x[i];

        ta2_xz_xzzz_z_0[i] = ta1_z_zzz_z_1[i] + ta2_xz_zzz_z_0[i] * pa_x[i] - ta2_xz_zzz_z_1[i] * pc_x[i];
    }

    // Set up 120-123 components of targeted buffer : GP

    auto ta2_xz_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 120);

    auto ta2_xz_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 121);

    auto ta2_xz_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 122);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xz_yy_x_0,   \
                             ta2_xz_yy_x_1,   \
                             ta2_xz_yy_y_0,   \
                             ta2_xz_yy_y_1,   \
                             ta2_xz_yy_z_0,   \
                             ta2_xz_yy_z_1,   \
                             ta2_xz_yyy_0_0,  \
                             ta2_xz_yyy_0_1,  \
                             ta2_xz_yyy_x_0,  \
                             ta2_xz_yyy_x_1,  \
                             ta2_xz_yyy_y_0,  \
                             ta2_xz_yyy_y_1,  \
                             ta2_xz_yyy_z_0,  \
                             ta2_xz_yyy_z_1,  \
                             ta2_xz_yyyy_x_0, \
                             ta2_xz_yyyy_y_0, \
                             ta2_xz_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyy_x_0[i] =
            3.0 * ta2_xz_yy_x_0[i] * fe_0 - 3.0 * ta2_xz_yy_x_1[i] * fe_0 + ta2_xz_yyy_x_0[i] * pa_y[i] - ta2_xz_yyy_x_1[i] * pc_y[i];

        ta2_xz_yyyy_y_0[i] = 3.0 * ta2_xz_yy_y_0[i] * fe_0 - 3.0 * ta2_xz_yy_y_1[i] * fe_0 + ta2_xz_yyy_0_0[i] * fe_0 - ta2_xz_yyy_0_1[i] * fe_0 +
                             ta2_xz_yyy_y_0[i] * pa_y[i] - ta2_xz_yyy_y_1[i] * pc_y[i];

        ta2_xz_yyyy_z_0[i] =
            3.0 * ta2_xz_yy_z_0[i] * fe_0 - 3.0 * ta2_xz_yy_z_1[i] * fe_0 + ta2_xz_yyy_z_0[i] * pa_y[i] - ta2_xz_yyy_z_1[i] * pc_y[i];
    }

    // Set up 123-126 components of targeted buffer : GP

    auto ta2_xz_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 123);

    auto ta2_xz_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 124);

    auto ta2_xz_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 125);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyy_x_1,   \
                             ta1_x_yyy_y_1,   \
                             ta2_xz_yyy_x_0,  \
                             ta2_xz_yyy_x_1,  \
                             ta2_xz_yyy_y_0,  \
                             ta2_xz_yyy_y_1,  \
                             ta2_xz_yyyz_x_0, \
                             ta2_xz_yyyz_y_0, \
                             ta2_xz_yyyz_z_0, \
                             ta2_xz_yyz_z_0,  \
                             ta2_xz_yyz_z_1,  \
                             ta2_xz_yz_z_0,   \
                             ta2_xz_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyz_x_0[i] = ta1_x_yyy_x_1[i] + ta2_xz_yyy_x_0[i] * pa_z[i] - ta2_xz_yyy_x_1[i] * pc_z[i];

        ta2_xz_yyyz_y_0[i] = ta1_x_yyy_y_1[i] + ta2_xz_yyy_y_0[i] * pa_z[i] - ta2_xz_yyy_y_1[i] * pc_z[i];

        ta2_xz_yyyz_z_0[i] =
            2.0 * ta2_xz_yz_z_0[i] * fe_0 - 2.0 * ta2_xz_yz_z_1[i] * fe_0 + ta2_xz_yyz_z_0[i] * pa_y[i] - ta2_xz_yyz_z_1[i] * pc_y[i];
    }

    // Set up 126-129 components of targeted buffer : GP

    auto ta2_xz_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 126);

    auto ta2_xz_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 127);

    auto ta2_xz_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 128);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyz_y_1,   \
                             ta2_xz_yy_y_0,   \
                             ta2_xz_yy_y_1,   \
                             ta2_xz_yyz_y_0,  \
                             ta2_xz_yyz_y_1,  \
                             ta2_xz_yyzz_x_0, \
                             ta2_xz_yyzz_y_0, \
                             ta2_xz_yyzz_z_0, \
                             ta2_xz_yzz_x_0,  \
                             ta2_xz_yzz_x_1,  \
                             ta2_xz_yzz_z_0,  \
                             ta2_xz_yzz_z_1,  \
                             ta2_xz_zz_x_0,   \
                             ta2_xz_zz_x_1,   \
                             ta2_xz_zz_z_0,   \
                             ta2_xz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyzz_x_0[i] = ta2_xz_zz_x_0[i] * fe_0 - ta2_xz_zz_x_1[i] * fe_0 + ta2_xz_yzz_x_0[i] * pa_y[i] - ta2_xz_yzz_x_1[i] * pc_y[i];

        ta2_xz_yyzz_y_0[i] =
            ta2_xz_yy_y_0[i] * fe_0 - ta2_xz_yy_y_1[i] * fe_0 + ta1_x_yyz_y_1[i] + ta2_xz_yyz_y_0[i] * pa_z[i] - ta2_xz_yyz_y_1[i] * pc_z[i];

        ta2_xz_yyzz_z_0[i] = ta2_xz_zz_z_0[i] * fe_0 - ta2_xz_zz_z_1[i] * fe_0 + ta2_xz_yzz_z_0[i] * pa_y[i] - ta2_xz_yzz_z_1[i] * pc_y[i];
    }

    // Set up 129-132 components of targeted buffer : GP

    auto ta2_xz_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 129);

    auto ta2_xz_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 130);

    auto ta2_xz_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 131);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xz_yzzz_x_0, \
                             ta2_xz_yzzz_y_0, \
                             ta2_xz_yzzz_z_0, \
                             ta2_xz_zzz_0_0,  \
                             ta2_xz_zzz_0_1,  \
                             ta2_xz_zzz_x_0,  \
                             ta2_xz_zzz_x_1,  \
                             ta2_xz_zzz_y_0,  \
                             ta2_xz_zzz_y_1,  \
                             ta2_xz_zzz_z_0,  \
                             ta2_xz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzzz_x_0[i] = ta2_xz_zzz_x_0[i] * pa_y[i] - ta2_xz_zzz_x_1[i] * pc_y[i];

        ta2_xz_yzzz_y_0[i] = ta2_xz_zzz_0_0[i] * fe_0 - ta2_xz_zzz_0_1[i] * fe_0 + ta2_xz_zzz_y_0[i] * pa_y[i] - ta2_xz_zzz_y_1[i] * pc_y[i];

        ta2_xz_yzzz_z_0[i] = ta2_xz_zzz_z_0[i] * pa_y[i] - ta2_xz_zzz_z_1[i] * pc_y[i];
    }

    // Set up 132-135 components of targeted buffer : GP

    auto ta2_xz_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 132);

    auto ta2_xz_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 133);

    auto ta2_xz_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 134);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_y_1,   \
                             ta1_x_zzz_z_1,   \
                             ta2_xz_zz_x_0,   \
                             ta2_xz_zz_x_1,   \
                             ta2_xz_zz_y_0,   \
                             ta2_xz_zz_y_1,   \
                             ta2_xz_zz_z_0,   \
                             ta2_xz_zz_z_1,   \
                             ta2_xz_zzz_0_0,  \
                             ta2_xz_zzz_0_1,  \
                             ta2_xz_zzz_x_0,  \
                             ta2_xz_zzz_x_1,  \
                             ta2_xz_zzz_y_0,  \
                             ta2_xz_zzz_y_1,  \
                             ta2_xz_zzz_z_0,  \
                             ta2_xz_zzz_z_1,  \
                             ta2_xz_zzzz_x_0, \
                             ta2_xz_zzzz_y_0, \
                             ta2_xz_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzzz_x_0[i] = 3.0 * ta2_xz_zz_x_0[i] * fe_0 - 3.0 * ta2_xz_zz_x_1[i] * fe_0 + ta1_x_zzz_x_1[i] + ta2_xz_zzz_x_0[i] * pa_z[i] -
                             ta2_xz_zzz_x_1[i] * pc_z[i];

        ta2_xz_zzzz_y_0[i] = 3.0 * ta2_xz_zz_y_0[i] * fe_0 - 3.0 * ta2_xz_zz_y_1[i] * fe_0 + ta1_x_zzz_y_1[i] + ta2_xz_zzz_y_0[i] * pa_z[i] -
                             ta2_xz_zzz_y_1[i] * pc_z[i];

        ta2_xz_zzzz_z_0[i] = 3.0 * ta2_xz_zz_z_0[i] * fe_0 - 3.0 * ta2_xz_zz_z_1[i] * fe_0 + ta2_xz_zzz_0_0[i] * fe_0 - ta2_xz_zzz_0_1[i] * fe_0 +
                             ta1_x_zzz_z_1[i] + ta2_xz_zzz_z_0[i] * pa_z[i] - ta2_xz_zzz_z_1[i] * pc_z[i];
    }

    // Set up 135-138 components of targeted buffer : GP

    auto ta2_yy_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 135);

    auto ta2_yy_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 136);

    auto ta2_yy_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 137);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yy_xx_x_0,   \
                             ta2_yy_xx_x_1,   \
                             ta2_yy_xx_y_0,   \
                             ta2_yy_xx_y_1,   \
                             ta2_yy_xx_z_0,   \
                             ta2_yy_xx_z_1,   \
                             ta2_yy_xxx_0_0,  \
                             ta2_yy_xxx_0_1,  \
                             ta2_yy_xxx_x_0,  \
                             ta2_yy_xxx_x_1,  \
                             ta2_yy_xxx_y_0,  \
                             ta2_yy_xxx_y_1,  \
                             ta2_yy_xxx_z_0,  \
                             ta2_yy_xxx_z_1,  \
                             ta2_yy_xxxx_x_0, \
                             ta2_yy_xxxx_y_0, \
                             ta2_yy_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxx_x_0[i] = 3.0 * ta2_yy_xx_x_0[i] * fe_0 - 3.0 * ta2_yy_xx_x_1[i] * fe_0 + ta2_yy_xxx_0_0[i] * fe_0 - ta2_yy_xxx_0_1[i] * fe_0 +
                             ta2_yy_xxx_x_0[i] * pa_x[i] - ta2_yy_xxx_x_1[i] * pc_x[i];

        ta2_yy_xxxx_y_0[i] =
            3.0 * ta2_yy_xx_y_0[i] * fe_0 - 3.0 * ta2_yy_xx_y_1[i] * fe_0 + ta2_yy_xxx_y_0[i] * pa_x[i] - ta2_yy_xxx_y_1[i] * pc_x[i];

        ta2_yy_xxxx_z_0[i] =
            3.0 * ta2_yy_xx_z_0[i] * fe_0 - 3.0 * ta2_yy_xx_z_1[i] * fe_0 + ta2_yy_xxx_z_0[i] * pa_x[i] - ta2_yy_xxx_z_1[i] * pc_x[i];
    }

    // Set up 138-141 components of targeted buffer : GP

    auto ta2_yy_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 138);

    auto ta2_yy_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 139);

    auto ta2_yy_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 140);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_z_1,   \
                             ta2_yy_xxx_x_0,  \
                             ta2_yy_xxx_x_1,  \
                             ta2_yy_xxx_z_0,  \
                             ta2_yy_xxx_z_1,  \
                             ta2_yy_xxxy_x_0, \
                             ta2_yy_xxxy_y_0, \
                             ta2_yy_xxxy_z_0, \
                             ta2_yy_xxy_y_0,  \
                             ta2_yy_xxy_y_1,  \
                             ta2_yy_xy_y_0,   \
                             ta2_yy_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxy_x_0[i] = 2.0 * ta1_y_xxx_x_1[i] + ta2_yy_xxx_x_0[i] * pa_y[i] - ta2_yy_xxx_x_1[i] * pc_y[i];

        ta2_yy_xxxy_y_0[i] =
            2.0 * ta2_yy_xy_y_0[i] * fe_0 - 2.0 * ta2_yy_xy_y_1[i] * fe_0 + ta2_yy_xxy_y_0[i] * pa_x[i] - ta2_yy_xxy_y_1[i] * pc_x[i];

        ta2_yy_xxxy_z_0[i] = 2.0 * ta1_y_xxx_z_1[i] + ta2_yy_xxx_z_0[i] * pa_y[i] - ta2_yy_xxx_z_1[i] * pc_y[i];
    }

    // Set up 141-144 components of targeted buffer : GP

    auto ta2_yy_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 141);

    auto ta2_yy_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 142);

    auto ta2_yy_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 143);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta2_yy_xxx_x_0,  \
                             ta2_yy_xxx_x_1,  \
                             ta2_yy_xxx_y_0,  \
                             ta2_yy_xxx_y_1,  \
                             ta2_yy_xxxz_x_0, \
                             ta2_yy_xxxz_y_0, \
                             ta2_yy_xxxz_z_0, \
                             ta2_yy_xxz_z_0,  \
                             ta2_yy_xxz_z_1,  \
                             ta2_yy_xz_z_0,   \
                             ta2_yy_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxz_x_0[i] = ta2_yy_xxx_x_0[i] * pa_z[i] - ta2_yy_xxx_x_1[i] * pc_z[i];

        ta2_yy_xxxz_y_0[i] = ta2_yy_xxx_y_0[i] * pa_z[i] - ta2_yy_xxx_y_1[i] * pc_z[i];

        ta2_yy_xxxz_z_0[i] =
            2.0 * ta2_yy_xz_z_0[i] * fe_0 - 2.0 * ta2_yy_xz_z_1[i] * fe_0 + ta2_yy_xxz_z_0[i] * pa_x[i] - ta2_yy_xxz_z_1[i] * pc_x[i];
    }

    // Set up 144-147 components of targeted buffer : GP

    auto ta2_yy_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 144);

    auto ta2_yy_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 145);

    auto ta2_yy_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 146);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxy_x_1,   \
                             ta2_yy_xx_x_0,   \
                             ta2_yy_xx_x_1,   \
                             ta2_yy_xxy_x_0,  \
                             ta2_yy_xxy_x_1,  \
                             ta2_yy_xxyy_x_0, \
                             ta2_yy_xxyy_y_0, \
                             ta2_yy_xxyy_z_0, \
                             ta2_yy_xyy_y_0,  \
                             ta2_yy_xyy_y_1,  \
                             ta2_yy_xyy_z_0,  \
                             ta2_yy_xyy_z_1,  \
                             ta2_yy_yy_y_0,   \
                             ta2_yy_yy_y_1,   \
                             ta2_yy_yy_z_0,   \
                             ta2_yy_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyy_x_0[i] =
            ta2_yy_xx_x_0[i] * fe_0 - ta2_yy_xx_x_1[i] * fe_0 + 2.0 * ta1_y_xxy_x_1[i] + ta2_yy_xxy_x_0[i] * pa_y[i] - ta2_yy_xxy_x_1[i] * pc_y[i];

        ta2_yy_xxyy_y_0[i] = ta2_yy_yy_y_0[i] * fe_0 - ta2_yy_yy_y_1[i] * fe_0 + ta2_yy_xyy_y_0[i] * pa_x[i] - ta2_yy_xyy_y_1[i] * pc_x[i];

        ta2_yy_xxyy_z_0[i] = ta2_yy_yy_z_0[i] * fe_0 - ta2_yy_yy_z_1[i] * fe_0 + ta2_yy_xyy_z_0[i] * pa_x[i] - ta2_yy_xyy_z_1[i] * pc_x[i];
    }

    // Set up 147-150 components of targeted buffer : GP

    auto ta2_yy_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 147);

    auto ta2_yy_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 148);

    auto ta2_yy_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 149);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_xxz_z_1,   \
                             ta2_yy_xxy_x_0,  \
                             ta2_yy_xxy_x_1,  \
                             ta2_yy_xxy_y_0,  \
                             ta2_yy_xxy_y_1,  \
                             ta2_yy_xxyz_x_0, \
                             ta2_yy_xxyz_y_0, \
                             ta2_yy_xxyz_z_0, \
                             ta2_yy_xxz_z_0,  \
                             ta2_yy_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xxyz_x_0[i] = ta2_yy_xxy_x_0[i] * pa_z[i] - ta2_yy_xxy_x_1[i] * pc_z[i];

        ta2_yy_xxyz_y_0[i] = ta2_yy_xxy_y_0[i] * pa_z[i] - ta2_yy_xxy_y_1[i] * pc_z[i];

        ta2_yy_xxyz_z_0[i] = 2.0 * ta1_y_xxz_z_1[i] + ta2_yy_xxz_z_0[i] * pa_y[i] - ta2_yy_xxz_z_1[i] * pc_y[i];
    }

    // Set up 150-153 components of targeted buffer : GP

    auto ta2_yy_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 150);

    auto ta2_yy_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 151);

    auto ta2_yy_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 152);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta2_yy_xx_x_0,   \
                             ta2_yy_xx_x_1,   \
                             ta2_yy_xxz_x_0,  \
                             ta2_yy_xxz_x_1,  \
                             ta2_yy_xxzz_x_0, \
                             ta2_yy_xxzz_y_0, \
                             ta2_yy_xxzz_z_0, \
                             ta2_yy_xzz_y_0,  \
                             ta2_yy_xzz_y_1,  \
                             ta2_yy_xzz_z_0,  \
                             ta2_yy_xzz_z_1,  \
                             ta2_yy_zz_y_0,   \
                             ta2_yy_zz_y_1,   \
                             ta2_yy_zz_z_0,   \
                             ta2_yy_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxzz_x_0[i] = ta2_yy_xx_x_0[i] * fe_0 - ta2_yy_xx_x_1[i] * fe_0 + ta2_yy_xxz_x_0[i] * pa_z[i] - ta2_yy_xxz_x_1[i] * pc_z[i];

        ta2_yy_xxzz_y_0[i] = ta2_yy_zz_y_0[i] * fe_0 - ta2_yy_zz_y_1[i] * fe_0 + ta2_yy_xzz_y_0[i] * pa_x[i] - ta2_yy_xzz_y_1[i] * pc_x[i];

        ta2_yy_xxzz_z_0[i] = ta2_yy_zz_z_0[i] * fe_0 - ta2_yy_zz_z_1[i] * fe_0 + ta2_yy_xzz_z_0[i] * pa_x[i] - ta2_yy_xzz_z_1[i] * pc_x[i];
    }

    // Set up 153-156 components of targeted buffer : GP

    auto ta2_yy_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 153);

    auto ta2_yy_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 154);

    auto ta2_yy_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 155);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yy_xyyy_x_0, \
                             ta2_yy_xyyy_y_0, \
                             ta2_yy_xyyy_z_0, \
                             ta2_yy_yyy_0_0,  \
                             ta2_yy_yyy_0_1,  \
                             ta2_yy_yyy_x_0,  \
                             ta2_yy_yyy_x_1,  \
                             ta2_yy_yyy_y_0,  \
                             ta2_yy_yyy_y_1,  \
                             ta2_yy_yyy_z_0,  \
                             ta2_yy_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyy_x_0[i] = ta2_yy_yyy_0_0[i] * fe_0 - ta2_yy_yyy_0_1[i] * fe_0 + ta2_yy_yyy_x_0[i] * pa_x[i] - ta2_yy_yyy_x_1[i] * pc_x[i];

        ta2_yy_xyyy_y_0[i] = ta2_yy_yyy_y_0[i] * pa_x[i] - ta2_yy_yyy_y_1[i] * pc_x[i];

        ta2_yy_xyyy_z_0[i] = ta2_yy_yyy_z_0[i] * pa_x[i] - ta2_yy_yyy_z_1[i] * pc_x[i];
    }

    // Set up 156-159 components of targeted buffer : GP

    auto ta2_yy_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 156);

    auto ta2_yy_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 157);

    auto ta2_yy_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 158);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta2_yy_xyy_x_0,  \
                             ta2_yy_xyy_x_1,  \
                             ta2_yy_xyyz_x_0, \
                             ta2_yy_xyyz_y_0, \
                             ta2_yy_xyyz_z_0, \
                             ta2_yy_yyz_y_0,  \
                             ta2_yy_yyz_y_1,  \
                             ta2_yy_yyz_z_0,  \
                             ta2_yy_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xyyz_x_0[i] = ta2_yy_xyy_x_0[i] * pa_z[i] - ta2_yy_xyy_x_1[i] * pc_z[i];

        ta2_yy_xyyz_y_0[i] = ta2_yy_yyz_y_0[i] * pa_x[i] - ta2_yy_yyz_y_1[i] * pc_x[i];

        ta2_yy_xyyz_z_0[i] = ta2_yy_yyz_z_0[i] * pa_x[i] - ta2_yy_yyz_z_1[i] * pc_x[i];
    }

    // Set up 159-162 components of targeted buffer : GP

    auto ta2_yy_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 159);

    auto ta2_yy_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 160);

    auto ta2_yy_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 161);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xzz_x_1,   \
                             ta2_yy_xyzz_x_0, \
                             ta2_yy_xyzz_y_0, \
                             ta2_yy_xyzz_z_0, \
                             ta2_yy_xzz_x_0,  \
                             ta2_yy_xzz_x_1,  \
                             ta2_yy_yzz_y_0,  \
                             ta2_yy_yzz_y_1,  \
                             ta2_yy_yzz_z_0,  \
                             ta2_yy_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xyzz_x_0[i] = 2.0 * ta1_y_xzz_x_1[i] + ta2_yy_xzz_x_0[i] * pa_y[i] - ta2_yy_xzz_x_1[i] * pc_y[i];

        ta2_yy_xyzz_y_0[i] = ta2_yy_yzz_y_0[i] * pa_x[i] - ta2_yy_yzz_y_1[i] * pc_x[i];

        ta2_yy_xyzz_z_0[i] = ta2_yy_yzz_z_0[i] * pa_x[i] - ta2_yy_yzz_z_1[i] * pc_x[i];
    }

    // Set up 162-165 components of targeted buffer : GP

    auto ta2_yy_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 162);

    auto ta2_yy_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 163);

    auto ta2_yy_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 164);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yy_xzzz_x_0, \
                             ta2_yy_xzzz_y_0, \
                             ta2_yy_xzzz_z_0, \
                             ta2_yy_zzz_0_0,  \
                             ta2_yy_zzz_0_1,  \
                             ta2_yy_zzz_x_0,  \
                             ta2_yy_zzz_x_1,  \
                             ta2_yy_zzz_y_0,  \
                             ta2_yy_zzz_y_1,  \
                             ta2_yy_zzz_z_0,  \
                             ta2_yy_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzzz_x_0[i] = ta2_yy_zzz_0_0[i] * fe_0 - ta2_yy_zzz_0_1[i] * fe_0 + ta2_yy_zzz_x_0[i] * pa_x[i] - ta2_yy_zzz_x_1[i] * pc_x[i];

        ta2_yy_xzzz_y_0[i] = ta2_yy_zzz_y_0[i] * pa_x[i] - ta2_yy_zzz_y_1[i] * pc_x[i];

        ta2_yy_xzzz_z_0[i] = ta2_yy_zzz_z_0[i] * pa_x[i] - ta2_yy_zzz_z_1[i] * pc_x[i];
    }

    // Set up 165-168 components of targeted buffer : GP

    auto ta2_yy_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 165);

    auto ta2_yy_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 166);

    auto ta2_yy_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 167);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_z_1,   \
                             ta2_yy_yy_x_0,   \
                             ta2_yy_yy_x_1,   \
                             ta2_yy_yy_y_0,   \
                             ta2_yy_yy_y_1,   \
                             ta2_yy_yy_z_0,   \
                             ta2_yy_yy_z_1,   \
                             ta2_yy_yyy_0_0,  \
                             ta2_yy_yyy_0_1,  \
                             ta2_yy_yyy_x_0,  \
                             ta2_yy_yyy_x_1,  \
                             ta2_yy_yyy_y_0,  \
                             ta2_yy_yyy_y_1,  \
                             ta2_yy_yyy_z_0,  \
                             ta2_yy_yyy_z_1,  \
                             ta2_yy_yyyy_x_0, \
                             ta2_yy_yyyy_y_0, \
                             ta2_yy_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyy_x_0[i] = 3.0 * ta2_yy_yy_x_0[i] * fe_0 - 3.0 * ta2_yy_yy_x_1[i] * fe_0 + 2.0 * ta1_y_yyy_x_1[i] + ta2_yy_yyy_x_0[i] * pa_y[i] -
                             ta2_yy_yyy_x_1[i] * pc_y[i];

        ta2_yy_yyyy_y_0[i] = 3.0 * ta2_yy_yy_y_0[i] * fe_0 - 3.0 * ta2_yy_yy_y_1[i] * fe_0 + ta2_yy_yyy_0_0[i] * fe_0 - ta2_yy_yyy_0_1[i] * fe_0 +
                             2.0 * ta1_y_yyy_y_1[i] + ta2_yy_yyy_y_0[i] * pa_y[i] - ta2_yy_yyy_y_1[i] * pc_y[i];

        ta2_yy_yyyy_z_0[i] = 3.0 * ta2_yy_yy_z_0[i] * fe_0 - 3.0 * ta2_yy_yy_z_1[i] * fe_0 + 2.0 * ta1_y_yyy_z_1[i] + ta2_yy_yyy_z_0[i] * pa_y[i] -
                             ta2_yy_yyy_z_1[i] * pc_y[i];
    }

    // Set up 168-171 components of targeted buffer : GP

    auto ta2_yy_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 168);

    auto ta2_yy_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 169);

    auto ta2_yy_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 170);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_yy_yyy_0_0,  \
                             ta2_yy_yyy_0_1,  \
                             ta2_yy_yyy_x_0,  \
                             ta2_yy_yyy_x_1,  \
                             ta2_yy_yyy_y_0,  \
                             ta2_yy_yyy_y_1,  \
                             ta2_yy_yyy_z_0,  \
                             ta2_yy_yyy_z_1,  \
                             ta2_yy_yyyz_x_0, \
                             ta2_yy_yyyz_y_0, \
                             ta2_yy_yyyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyz_x_0[i] = ta2_yy_yyy_x_0[i] * pa_z[i] - ta2_yy_yyy_x_1[i] * pc_z[i];

        ta2_yy_yyyz_y_0[i] = ta2_yy_yyy_y_0[i] * pa_z[i] - ta2_yy_yyy_y_1[i] * pc_z[i];

        ta2_yy_yyyz_z_0[i] = ta2_yy_yyy_0_0[i] * fe_0 - ta2_yy_yyy_0_1[i] * fe_0 + ta2_yy_yyy_z_0[i] * pa_z[i] - ta2_yy_yyy_z_1[i] * pc_z[i];
    }

    // Set up 171-174 components of targeted buffer : GP

    auto ta2_yy_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 171);

    auto ta2_yy_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 172);

    auto ta2_yy_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 173);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yzz_z_1,   \
                             ta2_yy_yy_x_0,   \
                             ta2_yy_yy_x_1,   \
                             ta2_yy_yy_y_0,   \
                             ta2_yy_yy_y_1,   \
                             ta2_yy_yyz_x_0,  \
                             ta2_yy_yyz_x_1,  \
                             ta2_yy_yyz_y_0,  \
                             ta2_yy_yyz_y_1,  \
                             ta2_yy_yyzz_x_0, \
                             ta2_yy_yyzz_y_0, \
                             ta2_yy_yyzz_z_0, \
                             ta2_yy_yzz_z_0,  \
                             ta2_yy_yzz_z_1,  \
                             ta2_yy_zz_z_0,   \
                             ta2_yy_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyzz_x_0[i] = ta2_yy_yy_x_0[i] * fe_0 - ta2_yy_yy_x_1[i] * fe_0 + ta2_yy_yyz_x_0[i] * pa_z[i] - ta2_yy_yyz_x_1[i] * pc_z[i];

        ta2_yy_yyzz_y_0[i] = ta2_yy_yy_y_0[i] * fe_0 - ta2_yy_yy_y_1[i] * fe_0 + ta2_yy_yyz_y_0[i] * pa_z[i] - ta2_yy_yyz_y_1[i] * pc_z[i];

        ta2_yy_yyzz_z_0[i] =
            ta2_yy_zz_z_0[i] * fe_0 - ta2_yy_zz_z_1[i] * fe_0 + 2.0 * ta1_y_yzz_z_1[i] + ta2_yy_yzz_z_0[i] * pa_y[i] - ta2_yy_yzz_z_1[i] * pc_y[i];
    }

    // Set up 174-177 components of targeted buffer : GP

    auto ta2_yy_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 174);

    auto ta2_yy_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 175);

    auto ta2_yy_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 176);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_zzz_x_1,   \
                             ta1_y_zzz_z_1,   \
                             ta2_yy_yz_y_0,   \
                             ta2_yy_yz_y_1,   \
                             ta2_yy_yzz_y_0,  \
                             ta2_yy_yzz_y_1,  \
                             ta2_yy_yzzz_x_0, \
                             ta2_yy_yzzz_y_0, \
                             ta2_yy_yzzz_z_0, \
                             ta2_yy_zzz_x_0,  \
                             ta2_yy_zzz_x_1,  \
                             ta2_yy_zzz_z_0,  \
                             ta2_yy_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzzz_x_0[i] = 2.0 * ta1_y_zzz_x_1[i] + ta2_yy_zzz_x_0[i] * pa_y[i] - ta2_yy_zzz_x_1[i] * pc_y[i];

        ta2_yy_yzzz_y_0[i] =
            2.0 * ta2_yy_yz_y_0[i] * fe_0 - 2.0 * ta2_yy_yz_y_1[i] * fe_0 + ta2_yy_yzz_y_0[i] * pa_z[i] - ta2_yy_yzz_y_1[i] * pc_z[i];

        ta2_yy_yzzz_z_0[i] = 2.0 * ta1_y_zzz_z_1[i] + ta2_yy_zzz_z_0[i] * pa_y[i] - ta2_yy_zzz_z_1[i] * pc_y[i];
    }

    // Set up 177-180 components of targeted buffer : GP

    auto ta2_yy_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 177);

    auto ta2_yy_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 178);

    auto ta2_yy_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 179);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_yy_zz_x_0,   \
                             ta2_yy_zz_x_1,   \
                             ta2_yy_zz_y_0,   \
                             ta2_yy_zz_y_1,   \
                             ta2_yy_zz_z_0,   \
                             ta2_yy_zz_z_1,   \
                             ta2_yy_zzz_0_0,  \
                             ta2_yy_zzz_0_1,  \
                             ta2_yy_zzz_x_0,  \
                             ta2_yy_zzz_x_1,  \
                             ta2_yy_zzz_y_0,  \
                             ta2_yy_zzz_y_1,  \
                             ta2_yy_zzz_z_0,  \
                             ta2_yy_zzz_z_1,  \
                             ta2_yy_zzzz_x_0, \
                             ta2_yy_zzzz_y_0, \
                             ta2_yy_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzzz_x_0[i] =
            3.0 * ta2_yy_zz_x_0[i] * fe_0 - 3.0 * ta2_yy_zz_x_1[i] * fe_0 + ta2_yy_zzz_x_0[i] * pa_z[i] - ta2_yy_zzz_x_1[i] * pc_z[i];

        ta2_yy_zzzz_y_0[i] =
            3.0 * ta2_yy_zz_y_0[i] * fe_0 - 3.0 * ta2_yy_zz_y_1[i] * fe_0 + ta2_yy_zzz_y_0[i] * pa_z[i] - ta2_yy_zzz_y_1[i] * pc_z[i];

        ta2_yy_zzzz_z_0[i] = 3.0 * ta2_yy_zz_z_0[i] * fe_0 - 3.0 * ta2_yy_zz_z_1[i] * fe_0 + ta2_yy_zzz_0_0[i] * fe_0 - ta2_yy_zzz_0_1[i] * fe_0 +
                             ta2_yy_zzz_z_0[i] * pa_z[i] - ta2_yy_zzz_z_1[i] * pc_z[i];
    }

    // Set up 180-183 components of targeted buffer : GP

    auto ta2_yz_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 180);

    auto ta2_yz_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 181);

    auto ta2_yz_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 182);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yz_xx_x_0,   \
                             ta2_yz_xx_x_1,   \
                             ta2_yz_xx_y_0,   \
                             ta2_yz_xx_y_1,   \
                             ta2_yz_xx_z_0,   \
                             ta2_yz_xx_z_1,   \
                             ta2_yz_xxx_0_0,  \
                             ta2_yz_xxx_0_1,  \
                             ta2_yz_xxx_x_0,  \
                             ta2_yz_xxx_x_1,  \
                             ta2_yz_xxx_y_0,  \
                             ta2_yz_xxx_y_1,  \
                             ta2_yz_xxx_z_0,  \
                             ta2_yz_xxx_z_1,  \
                             ta2_yz_xxxx_x_0, \
                             ta2_yz_xxxx_y_0, \
                             ta2_yz_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxx_x_0[i] = 3.0 * ta2_yz_xx_x_0[i] * fe_0 - 3.0 * ta2_yz_xx_x_1[i] * fe_0 + ta2_yz_xxx_0_0[i] * fe_0 - ta2_yz_xxx_0_1[i] * fe_0 +
                             ta2_yz_xxx_x_0[i] * pa_x[i] - ta2_yz_xxx_x_1[i] * pc_x[i];

        ta2_yz_xxxx_y_0[i] =
            3.0 * ta2_yz_xx_y_0[i] * fe_0 - 3.0 * ta2_yz_xx_y_1[i] * fe_0 + ta2_yz_xxx_y_0[i] * pa_x[i] - ta2_yz_xxx_y_1[i] * pc_x[i];

        ta2_yz_xxxx_z_0[i] =
            3.0 * ta2_yz_xx_z_0[i] * fe_0 - 3.0 * ta2_yz_xx_z_1[i] * fe_0 + ta2_yz_xxx_z_0[i] * pa_x[i] - ta2_yz_xxx_z_1[i] * pc_x[i];
    }

    // Set up 183-186 components of targeted buffer : GP

    auto ta2_yz_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 183);

    auto ta2_yz_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 184);

    auto ta2_yz_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 185);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_z_1,   \
                             ta2_yz_xxx_x_0,  \
                             ta2_yz_xxx_x_1,  \
                             ta2_yz_xxx_z_0,  \
                             ta2_yz_xxx_z_1,  \
                             ta2_yz_xxxy_x_0, \
                             ta2_yz_xxxy_y_0, \
                             ta2_yz_xxxy_z_0, \
                             ta2_yz_xxy_y_0,  \
                             ta2_yz_xxy_y_1,  \
                             ta2_yz_xy_y_0,   \
                             ta2_yz_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxy_x_0[i] = ta1_z_xxx_x_1[i] + ta2_yz_xxx_x_0[i] * pa_y[i] - ta2_yz_xxx_x_1[i] * pc_y[i];

        ta2_yz_xxxy_y_0[i] =
            2.0 * ta2_yz_xy_y_0[i] * fe_0 - 2.0 * ta2_yz_xy_y_1[i] * fe_0 + ta2_yz_xxy_y_0[i] * pa_x[i] - ta2_yz_xxy_y_1[i] * pc_x[i];

        ta2_yz_xxxy_z_0[i] = ta1_z_xxx_z_1[i] + ta2_yz_xxx_z_0[i] * pa_y[i] - ta2_yz_xxx_z_1[i] * pc_y[i];
    }

    // Set up 186-189 components of targeted buffer : GP

    auto ta2_yz_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 186);

    auto ta2_yz_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 187);

    auto ta2_yz_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 188);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_y_1,   \
                             ta2_yz_xxx_x_0,  \
                             ta2_yz_xxx_x_1,  \
                             ta2_yz_xxx_y_0,  \
                             ta2_yz_xxx_y_1,  \
                             ta2_yz_xxxz_x_0, \
                             ta2_yz_xxxz_y_0, \
                             ta2_yz_xxxz_z_0, \
                             ta2_yz_xxz_z_0,  \
                             ta2_yz_xxz_z_1,  \
                             ta2_yz_xz_z_0,   \
                             ta2_yz_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxz_x_0[i] = ta1_y_xxx_x_1[i] + ta2_yz_xxx_x_0[i] * pa_z[i] - ta2_yz_xxx_x_1[i] * pc_z[i];

        ta2_yz_xxxz_y_0[i] = ta1_y_xxx_y_1[i] + ta2_yz_xxx_y_0[i] * pa_z[i] - ta2_yz_xxx_y_1[i] * pc_z[i];

        ta2_yz_xxxz_z_0[i] =
            2.0 * ta2_yz_xz_z_0[i] * fe_0 - 2.0 * ta2_yz_xz_z_1[i] * fe_0 + ta2_yz_xxz_z_0[i] * pa_x[i] - ta2_yz_xxz_z_1[i] * pc_x[i];
    }

    // Set up 189-192 components of targeted buffer : GP

    auto ta2_yz_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 189);

    auto ta2_yz_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 190);

    auto ta2_yz_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 191);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxy_x_1,   \
                             ta2_yz_xx_x_0,   \
                             ta2_yz_xx_x_1,   \
                             ta2_yz_xxy_x_0,  \
                             ta2_yz_xxy_x_1,  \
                             ta2_yz_xxyy_x_0, \
                             ta2_yz_xxyy_y_0, \
                             ta2_yz_xxyy_z_0, \
                             ta2_yz_xyy_y_0,  \
                             ta2_yz_xyy_y_1,  \
                             ta2_yz_xyy_z_0,  \
                             ta2_yz_xyy_z_1,  \
                             ta2_yz_yy_y_0,   \
                             ta2_yz_yy_y_1,   \
                             ta2_yz_yy_z_0,   \
                             ta2_yz_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyy_x_0[i] =
            ta2_yz_xx_x_0[i] * fe_0 - ta2_yz_xx_x_1[i] * fe_0 + ta1_z_xxy_x_1[i] + ta2_yz_xxy_x_0[i] * pa_y[i] - ta2_yz_xxy_x_1[i] * pc_y[i];

        ta2_yz_xxyy_y_0[i] = ta2_yz_yy_y_0[i] * fe_0 - ta2_yz_yy_y_1[i] * fe_0 + ta2_yz_xyy_y_0[i] * pa_x[i] - ta2_yz_xyy_y_1[i] * pc_x[i];

        ta2_yz_xxyy_z_0[i] = ta2_yz_yy_z_0[i] * fe_0 - ta2_yz_yy_z_1[i] * fe_0 + ta2_yz_xyy_z_0[i] * pa_x[i] - ta2_yz_xyy_z_1[i] * pc_x[i];
    }

    // Set up 192-195 components of targeted buffer : GP

    auto ta2_yz_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 192);

    auto ta2_yz_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 193);

    auto ta2_yz_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 194);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_xxy_y_1,   \
                             ta1_z_xxz_x_1,   \
                             ta1_z_xxz_z_1,   \
                             ta2_yz_xxy_y_0,  \
                             ta2_yz_xxy_y_1,  \
                             ta2_yz_xxyz_x_0, \
                             ta2_yz_xxyz_y_0, \
                             ta2_yz_xxyz_z_0, \
                             ta2_yz_xxz_x_0,  \
                             ta2_yz_xxz_x_1,  \
                             ta2_yz_xxz_z_0,  \
                             ta2_yz_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xxyz_x_0[i] = ta1_z_xxz_x_1[i] + ta2_yz_xxz_x_0[i] * pa_y[i] - ta2_yz_xxz_x_1[i] * pc_y[i];

        ta2_yz_xxyz_y_0[i] = ta1_y_xxy_y_1[i] + ta2_yz_xxy_y_0[i] * pa_z[i] - ta2_yz_xxy_y_1[i] * pc_z[i];

        ta2_yz_xxyz_z_0[i] = ta1_z_xxz_z_1[i] + ta2_yz_xxz_z_0[i] * pa_y[i] - ta2_yz_xxz_z_1[i] * pc_y[i];
    }

    // Set up 195-198 components of targeted buffer : GP

    auto ta2_yz_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 195);

    auto ta2_yz_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 196);

    auto ta2_yz_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 197);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxz_x_1,   \
                             ta2_yz_xx_x_0,   \
                             ta2_yz_xx_x_1,   \
                             ta2_yz_xxz_x_0,  \
                             ta2_yz_xxz_x_1,  \
                             ta2_yz_xxzz_x_0, \
                             ta2_yz_xxzz_y_0, \
                             ta2_yz_xxzz_z_0, \
                             ta2_yz_xzz_y_0,  \
                             ta2_yz_xzz_y_1,  \
                             ta2_yz_xzz_z_0,  \
                             ta2_yz_xzz_z_1,  \
                             ta2_yz_zz_y_0,   \
                             ta2_yz_zz_y_1,   \
                             ta2_yz_zz_z_0,   \
                             ta2_yz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxzz_x_0[i] =
            ta2_yz_xx_x_0[i] * fe_0 - ta2_yz_xx_x_1[i] * fe_0 + ta1_y_xxz_x_1[i] + ta2_yz_xxz_x_0[i] * pa_z[i] - ta2_yz_xxz_x_1[i] * pc_z[i];

        ta2_yz_xxzz_y_0[i] = ta2_yz_zz_y_0[i] * fe_0 - ta2_yz_zz_y_1[i] * fe_0 + ta2_yz_xzz_y_0[i] * pa_x[i] - ta2_yz_xzz_y_1[i] * pc_x[i];

        ta2_yz_xxzz_z_0[i] = ta2_yz_zz_z_0[i] * fe_0 - ta2_yz_zz_z_1[i] * fe_0 + ta2_yz_xzz_z_0[i] * pa_x[i] - ta2_yz_xzz_z_1[i] * pc_x[i];
    }

    // Set up 198-201 components of targeted buffer : GP

    auto ta2_yz_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 198);

    auto ta2_yz_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 199);

    auto ta2_yz_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 200);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yz_xyyy_x_0, \
                             ta2_yz_xyyy_y_0, \
                             ta2_yz_xyyy_z_0, \
                             ta2_yz_yyy_0_0,  \
                             ta2_yz_yyy_0_1,  \
                             ta2_yz_yyy_x_0,  \
                             ta2_yz_yyy_x_1,  \
                             ta2_yz_yyy_y_0,  \
                             ta2_yz_yyy_y_1,  \
                             ta2_yz_yyy_z_0,  \
                             ta2_yz_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyy_x_0[i] = ta2_yz_yyy_0_0[i] * fe_0 - ta2_yz_yyy_0_1[i] * fe_0 + ta2_yz_yyy_x_0[i] * pa_x[i] - ta2_yz_yyy_x_1[i] * pc_x[i];

        ta2_yz_xyyy_y_0[i] = ta2_yz_yyy_y_0[i] * pa_x[i] - ta2_yz_yyy_y_1[i] * pc_x[i];

        ta2_yz_xyyy_z_0[i] = ta2_yz_yyy_z_0[i] * pa_x[i] - ta2_yz_yyy_z_1[i] * pc_x[i];
    }

    // Set up 201-204 components of targeted buffer : GP

    auto ta2_yz_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 201);

    auto ta2_yz_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 202);

    auto ta2_yz_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 203);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xyy_x_1,   \
                             ta2_yz_xyy_x_0,  \
                             ta2_yz_xyy_x_1,  \
                             ta2_yz_xyyz_x_0, \
                             ta2_yz_xyyz_y_0, \
                             ta2_yz_xyyz_z_0, \
                             ta2_yz_yyz_y_0,  \
                             ta2_yz_yyz_y_1,  \
                             ta2_yz_yyz_z_0,  \
                             ta2_yz_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xyyz_x_0[i] = ta1_y_xyy_x_1[i] + ta2_yz_xyy_x_0[i] * pa_z[i] - ta2_yz_xyy_x_1[i] * pc_z[i];

        ta2_yz_xyyz_y_0[i] = ta2_yz_yyz_y_0[i] * pa_x[i] - ta2_yz_yyz_y_1[i] * pc_x[i];

        ta2_yz_xyyz_z_0[i] = ta2_yz_yyz_z_0[i] * pa_x[i] - ta2_yz_yyz_z_1[i] * pc_x[i];
    }

    // Set up 204-207 components of targeted buffer : GP

    auto ta2_yz_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 204);

    auto ta2_yz_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 205);

    auto ta2_yz_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 206);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xzz_x_1,   \
                             ta2_yz_xyzz_x_0, \
                             ta2_yz_xyzz_y_0, \
                             ta2_yz_xyzz_z_0, \
                             ta2_yz_xzz_x_0,  \
                             ta2_yz_xzz_x_1,  \
                             ta2_yz_yzz_y_0,  \
                             ta2_yz_yzz_y_1,  \
                             ta2_yz_yzz_z_0,  \
                             ta2_yz_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xyzz_x_0[i] = ta1_z_xzz_x_1[i] + ta2_yz_xzz_x_0[i] * pa_y[i] - ta2_yz_xzz_x_1[i] * pc_y[i];

        ta2_yz_xyzz_y_0[i] = ta2_yz_yzz_y_0[i] * pa_x[i] - ta2_yz_yzz_y_1[i] * pc_x[i];

        ta2_yz_xyzz_z_0[i] = ta2_yz_yzz_z_0[i] * pa_x[i] - ta2_yz_yzz_z_1[i] * pc_x[i];
    }

    // Set up 207-210 components of targeted buffer : GP

    auto ta2_yz_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 207);

    auto ta2_yz_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 208);

    auto ta2_yz_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 209);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yz_xzzz_x_0, \
                             ta2_yz_xzzz_y_0, \
                             ta2_yz_xzzz_z_0, \
                             ta2_yz_zzz_0_0,  \
                             ta2_yz_zzz_0_1,  \
                             ta2_yz_zzz_x_0,  \
                             ta2_yz_zzz_x_1,  \
                             ta2_yz_zzz_y_0,  \
                             ta2_yz_zzz_y_1,  \
                             ta2_yz_zzz_z_0,  \
                             ta2_yz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzzz_x_0[i] = ta2_yz_zzz_0_0[i] * fe_0 - ta2_yz_zzz_0_1[i] * fe_0 + ta2_yz_zzz_x_0[i] * pa_x[i] - ta2_yz_zzz_x_1[i] * pc_x[i];

        ta2_yz_xzzz_y_0[i] = ta2_yz_zzz_y_0[i] * pa_x[i] - ta2_yz_zzz_y_1[i] * pc_x[i];

        ta2_yz_xzzz_z_0[i] = ta2_yz_zzz_z_0[i] * pa_x[i] - ta2_yz_zzz_z_1[i] * pc_x[i];
    }

    // Set up 210-213 components of targeted buffer : GP

    auto ta2_yz_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 210);

    auto ta2_yz_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 211);

    auto ta2_yz_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 212);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yyy_x_1,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_z_1,   \
                             ta2_yz_yy_x_0,   \
                             ta2_yz_yy_x_1,   \
                             ta2_yz_yy_y_0,   \
                             ta2_yz_yy_y_1,   \
                             ta2_yz_yy_z_0,   \
                             ta2_yz_yy_z_1,   \
                             ta2_yz_yyy_0_0,  \
                             ta2_yz_yyy_0_1,  \
                             ta2_yz_yyy_x_0,  \
                             ta2_yz_yyy_x_1,  \
                             ta2_yz_yyy_y_0,  \
                             ta2_yz_yyy_y_1,  \
                             ta2_yz_yyy_z_0,  \
                             ta2_yz_yyy_z_1,  \
                             ta2_yz_yyyy_x_0, \
                             ta2_yz_yyyy_y_0, \
                             ta2_yz_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyy_x_0[i] = 3.0 * ta2_yz_yy_x_0[i] * fe_0 - 3.0 * ta2_yz_yy_x_1[i] * fe_0 + ta1_z_yyy_x_1[i] + ta2_yz_yyy_x_0[i] * pa_y[i] -
                             ta2_yz_yyy_x_1[i] * pc_y[i];

        ta2_yz_yyyy_y_0[i] = 3.0 * ta2_yz_yy_y_0[i] * fe_0 - 3.0 * ta2_yz_yy_y_1[i] * fe_0 + ta2_yz_yyy_0_0[i] * fe_0 - ta2_yz_yyy_0_1[i] * fe_0 +
                             ta1_z_yyy_y_1[i] + ta2_yz_yyy_y_0[i] * pa_y[i] - ta2_yz_yyy_y_1[i] * pc_y[i];

        ta2_yz_yyyy_z_0[i] = 3.0 * ta2_yz_yy_z_0[i] * fe_0 - 3.0 * ta2_yz_yy_z_1[i] * fe_0 + ta1_z_yyy_z_1[i] + ta2_yz_yyy_z_0[i] * pa_y[i] -
                             ta2_yz_yyy_z_1[i] * pc_y[i];
    }

    // Set up 213-216 components of targeted buffer : GP

    auto ta2_yz_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 213);

    auto ta2_yz_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 214);

    auto ta2_yz_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 215);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_y_1,   \
                             ta1_z_yyz_z_1,   \
                             ta2_yz_yyy_x_0,  \
                             ta2_yz_yyy_x_1,  \
                             ta2_yz_yyy_y_0,  \
                             ta2_yz_yyy_y_1,  \
                             ta2_yz_yyyz_x_0, \
                             ta2_yz_yyyz_y_0, \
                             ta2_yz_yyyz_z_0, \
                             ta2_yz_yyz_z_0,  \
                             ta2_yz_yyz_z_1,  \
                             ta2_yz_yz_z_0,   \
                             ta2_yz_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyz_x_0[i] = ta1_y_yyy_x_1[i] + ta2_yz_yyy_x_0[i] * pa_z[i] - ta2_yz_yyy_x_1[i] * pc_z[i];

        ta2_yz_yyyz_y_0[i] = ta1_y_yyy_y_1[i] + ta2_yz_yyy_y_0[i] * pa_z[i] - ta2_yz_yyy_y_1[i] * pc_z[i];

        ta2_yz_yyyz_z_0[i] = 2.0 * ta2_yz_yz_z_0[i] * fe_0 - 2.0 * ta2_yz_yz_z_1[i] * fe_0 + ta1_z_yyz_z_1[i] + ta2_yz_yyz_z_0[i] * pa_y[i] -
                             ta2_yz_yyz_z_1[i] * pc_y[i];
    }

    // Set up 216-219 components of targeted buffer : GP

    auto ta2_yz_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 216);

    auto ta2_yz_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 217);

    auto ta2_yz_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 218);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yyz_y_1,   \
                             ta1_z_yzz_x_1,   \
                             ta1_z_yzz_z_1,   \
                             ta2_yz_yy_y_0,   \
                             ta2_yz_yy_y_1,   \
                             ta2_yz_yyz_y_0,  \
                             ta2_yz_yyz_y_1,  \
                             ta2_yz_yyzz_x_0, \
                             ta2_yz_yyzz_y_0, \
                             ta2_yz_yyzz_z_0, \
                             ta2_yz_yzz_x_0,  \
                             ta2_yz_yzz_x_1,  \
                             ta2_yz_yzz_z_0,  \
                             ta2_yz_yzz_z_1,  \
                             ta2_yz_zz_x_0,   \
                             ta2_yz_zz_x_1,   \
                             ta2_yz_zz_z_0,   \
                             ta2_yz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyzz_x_0[i] =
            ta2_yz_zz_x_0[i] * fe_0 - ta2_yz_zz_x_1[i] * fe_0 + ta1_z_yzz_x_1[i] + ta2_yz_yzz_x_0[i] * pa_y[i] - ta2_yz_yzz_x_1[i] * pc_y[i];

        ta2_yz_yyzz_y_0[i] =
            ta2_yz_yy_y_0[i] * fe_0 - ta2_yz_yy_y_1[i] * fe_0 + ta1_y_yyz_y_1[i] + ta2_yz_yyz_y_0[i] * pa_z[i] - ta2_yz_yyz_y_1[i] * pc_z[i];

        ta2_yz_yyzz_z_0[i] =
            ta2_yz_zz_z_0[i] * fe_0 - ta2_yz_zz_z_1[i] * fe_0 + ta1_z_yzz_z_1[i] + ta2_yz_yzz_z_0[i] * pa_y[i] - ta2_yz_yzz_z_1[i] * pc_y[i];
    }

    // Set up 219-222 components of targeted buffer : GP

    auto ta2_yz_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 219);

    auto ta2_yz_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 220);

    auto ta2_yz_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 221);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_z_1,   \
                             ta2_yz_yzzz_x_0, \
                             ta2_yz_yzzz_y_0, \
                             ta2_yz_yzzz_z_0, \
                             ta2_yz_zzz_0_0,  \
                             ta2_yz_zzz_0_1,  \
                             ta2_yz_zzz_x_0,  \
                             ta2_yz_zzz_x_1,  \
                             ta2_yz_zzz_y_0,  \
                             ta2_yz_zzz_y_1,  \
                             ta2_yz_zzz_z_0,  \
                             ta2_yz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzzz_x_0[i] = ta1_z_zzz_x_1[i] + ta2_yz_zzz_x_0[i] * pa_y[i] - ta2_yz_zzz_x_1[i] * pc_y[i];

        ta2_yz_yzzz_y_0[i] =
            ta2_yz_zzz_0_0[i] * fe_0 - ta2_yz_zzz_0_1[i] * fe_0 + ta1_z_zzz_y_1[i] + ta2_yz_zzz_y_0[i] * pa_y[i] - ta2_yz_zzz_y_1[i] * pc_y[i];

        ta2_yz_yzzz_z_0[i] = ta1_z_zzz_z_1[i] + ta2_yz_zzz_z_0[i] * pa_y[i] - ta2_yz_zzz_z_1[i] * pc_y[i];
    }

    // Set up 222-225 components of targeted buffer : GP

    auto ta2_yz_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 222);

    auto ta2_yz_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 223);

    auto ta2_yz_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 224);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_zzz_x_1,   \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_z_1,   \
                             ta2_yz_zz_x_0,   \
                             ta2_yz_zz_x_1,   \
                             ta2_yz_zz_y_0,   \
                             ta2_yz_zz_y_1,   \
                             ta2_yz_zz_z_0,   \
                             ta2_yz_zz_z_1,   \
                             ta2_yz_zzz_0_0,  \
                             ta2_yz_zzz_0_1,  \
                             ta2_yz_zzz_x_0,  \
                             ta2_yz_zzz_x_1,  \
                             ta2_yz_zzz_y_0,  \
                             ta2_yz_zzz_y_1,  \
                             ta2_yz_zzz_z_0,  \
                             ta2_yz_zzz_z_1,  \
                             ta2_yz_zzzz_x_0, \
                             ta2_yz_zzzz_y_0, \
                             ta2_yz_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzzz_x_0[i] = 3.0 * ta2_yz_zz_x_0[i] * fe_0 - 3.0 * ta2_yz_zz_x_1[i] * fe_0 + ta1_y_zzz_x_1[i] + ta2_yz_zzz_x_0[i] * pa_z[i] -
                             ta2_yz_zzz_x_1[i] * pc_z[i];

        ta2_yz_zzzz_y_0[i] = 3.0 * ta2_yz_zz_y_0[i] * fe_0 - 3.0 * ta2_yz_zz_y_1[i] * fe_0 + ta1_y_zzz_y_1[i] + ta2_yz_zzz_y_0[i] * pa_z[i] -
                             ta2_yz_zzz_y_1[i] * pc_z[i];

        ta2_yz_zzzz_z_0[i] = 3.0 * ta2_yz_zz_z_0[i] * fe_0 - 3.0 * ta2_yz_zz_z_1[i] * fe_0 + ta2_yz_zzz_0_0[i] * fe_0 - ta2_yz_zzz_0_1[i] * fe_0 +
                             ta1_y_zzz_z_1[i] + ta2_yz_zzz_z_0[i] * pa_z[i] - ta2_yz_zzz_z_1[i] * pc_z[i];
    }

    // Set up 225-228 components of targeted buffer : GP

    auto ta2_zz_xxxx_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 225);

    auto ta2_zz_xxxx_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 226);

    auto ta2_zz_xxxx_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 227);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_zz_xx_x_0,   \
                             ta2_zz_xx_x_1,   \
                             ta2_zz_xx_y_0,   \
                             ta2_zz_xx_y_1,   \
                             ta2_zz_xx_z_0,   \
                             ta2_zz_xx_z_1,   \
                             ta2_zz_xxx_0_0,  \
                             ta2_zz_xxx_0_1,  \
                             ta2_zz_xxx_x_0,  \
                             ta2_zz_xxx_x_1,  \
                             ta2_zz_xxx_y_0,  \
                             ta2_zz_xxx_y_1,  \
                             ta2_zz_xxx_z_0,  \
                             ta2_zz_xxx_z_1,  \
                             ta2_zz_xxxx_x_0, \
                             ta2_zz_xxxx_y_0, \
                             ta2_zz_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxx_x_0[i] = 3.0 * ta2_zz_xx_x_0[i] * fe_0 - 3.0 * ta2_zz_xx_x_1[i] * fe_0 + ta2_zz_xxx_0_0[i] * fe_0 - ta2_zz_xxx_0_1[i] * fe_0 +
                             ta2_zz_xxx_x_0[i] * pa_x[i] - ta2_zz_xxx_x_1[i] * pc_x[i];

        ta2_zz_xxxx_y_0[i] =
            3.0 * ta2_zz_xx_y_0[i] * fe_0 - 3.0 * ta2_zz_xx_y_1[i] * fe_0 + ta2_zz_xxx_y_0[i] * pa_x[i] - ta2_zz_xxx_y_1[i] * pc_x[i];

        ta2_zz_xxxx_z_0[i] =
            3.0 * ta2_zz_xx_z_0[i] * fe_0 - 3.0 * ta2_zz_xx_z_1[i] * fe_0 + ta2_zz_xxx_z_0[i] * pa_x[i] - ta2_zz_xxx_z_1[i] * pc_x[i];
    }

    // Set up 228-231 components of targeted buffer : GP

    auto ta2_zz_xxxy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 228);

    auto ta2_zz_xxxy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 229);

    auto ta2_zz_xxxy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 230);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta2_zz_xxx_x_0,  \
                             ta2_zz_xxx_x_1,  \
                             ta2_zz_xxx_z_0,  \
                             ta2_zz_xxx_z_1,  \
                             ta2_zz_xxxy_x_0, \
                             ta2_zz_xxxy_y_0, \
                             ta2_zz_xxxy_z_0, \
                             ta2_zz_xxy_y_0,  \
                             ta2_zz_xxy_y_1,  \
                             ta2_zz_xy_y_0,   \
                             ta2_zz_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxy_x_0[i] = ta2_zz_xxx_x_0[i] * pa_y[i] - ta2_zz_xxx_x_1[i] * pc_y[i];

        ta2_zz_xxxy_y_0[i] =
            2.0 * ta2_zz_xy_y_0[i] * fe_0 - 2.0 * ta2_zz_xy_y_1[i] * fe_0 + ta2_zz_xxy_y_0[i] * pa_x[i] - ta2_zz_xxy_y_1[i] * pc_x[i];

        ta2_zz_xxxy_z_0[i] = ta2_zz_xxx_z_0[i] * pa_y[i] - ta2_zz_xxx_z_1[i] * pc_y[i];
    }

    // Set up 231-234 components of targeted buffer : GP

    auto ta2_zz_xxxz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 231);

    auto ta2_zz_xxxz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 232);

    auto ta2_zz_xxxz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 233);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_y_1,   \
                             ta2_zz_xxx_x_0,  \
                             ta2_zz_xxx_x_1,  \
                             ta2_zz_xxx_y_0,  \
                             ta2_zz_xxx_y_1,  \
                             ta2_zz_xxxz_x_0, \
                             ta2_zz_xxxz_y_0, \
                             ta2_zz_xxxz_z_0, \
                             ta2_zz_xxz_z_0,  \
                             ta2_zz_xxz_z_1,  \
                             ta2_zz_xz_z_0,   \
                             ta2_zz_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxz_x_0[i] = 2.0 * ta1_z_xxx_x_1[i] + ta2_zz_xxx_x_0[i] * pa_z[i] - ta2_zz_xxx_x_1[i] * pc_z[i];

        ta2_zz_xxxz_y_0[i] = 2.0 * ta1_z_xxx_y_1[i] + ta2_zz_xxx_y_0[i] * pa_z[i] - ta2_zz_xxx_y_1[i] * pc_z[i];

        ta2_zz_xxxz_z_0[i] =
            2.0 * ta2_zz_xz_z_0[i] * fe_0 - 2.0 * ta2_zz_xz_z_1[i] * fe_0 + ta2_zz_xxz_z_0[i] * pa_x[i] - ta2_zz_xxz_z_1[i] * pc_x[i];
    }

    // Set up 234-237 components of targeted buffer : GP

    auto ta2_zz_xxyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 234);

    auto ta2_zz_xxyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 235);

    auto ta2_zz_xxyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 236);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta2_zz_xx_x_0,   \
                             ta2_zz_xx_x_1,   \
                             ta2_zz_xxy_x_0,  \
                             ta2_zz_xxy_x_1,  \
                             ta2_zz_xxyy_x_0, \
                             ta2_zz_xxyy_y_0, \
                             ta2_zz_xxyy_z_0, \
                             ta2_zz_xyy_y_0,  \
                             ta2_zz_xyy_y_1,  \
                             ta2_zz_xyy_z_0,  \
                             ta2_zz_xyy_z_1,  \
                             ta2_zz_yy_y_0,   \
                             ta2_zz_yy_y_1,   \
                             ta2_zz_yy_z_0,   \
                             ta2_zz_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyy_x_0[i] = ta2_zz_xx_x_0[i] * fe_0 - ta2_zz_xx_x_1[i] * fe_0 + ta2_zz_xxy_x_0[i] * pa_y[i] - ta2_zz_xxy_x_1[i] * pc_y[i];

        ta2_zz_xxyy_y_0[i] = ta2_zz_yy_y_0[i] * fe_0 - ta2_zz_yy_y_1[i] * fe_0 + ta2_zz_xyy_y_0[i] * pa_x[i] - ta2_zz_xyy_y_1[i] * pc_x[i];

        ta2_zz_xxyy_z_0[i] = ta2_zz_yy_z_0[i] * fe_0 - ta2_zz_yy_z_1[i] * fe_0 + ta2_zz_xyy_z_0[i] * pa_x[i] - ta2_zz_xyy_z_1[i] * pc_x[i];
    }

    // Set up 237-240 components of targeted buffer : GP

    auto ta2_zz_xxyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 237);

    auto ta2_zz_xxyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 238);

    auto ta2_zz_xxyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 239);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_xxy_y_1,   \
                             ta2_zz_xxy_y_0,  \
                             ta2_zz_xxy_y_1,  \
                             ta2_zz_xxyz_x_0, \
                             ta2_zz_xxyz_y_0, \
                             ta2_zz_xxyz_z_0, \
                             ta2_zz_xxz_x_0,  \
                             ta2_zz_xxz_x_1,  \
                             ta2_zz_xxz_z_0,  \
                             ta2_zz_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xxyz_x_0[i] = ta2_zz_xxz_x_0[i] * pa_y[i] - ta2_zz_xxz_x_1[i] * pc_y[i];

        ta2_zz_xxyz_y_0[i] = 2.0 * ta1_z_xxy_y_1[i] + ta2_zz_xxy_y_0[i] * pa_z[i] - ta2_zz_xxy_y_1[i] * pc_z[i];

        ta2_zz_xxyz_z_0[i] = ta2_zz_xxz_z_0[i] * pa_y[i] - ta2_zz_xxz_z_1[i] * pc_y[i];
    }

    // Set up 240-243 components of targeted buffer : GP

    auto ta2_zz_xxzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 240);

    auto ta2_zz_xxzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 241);

    auto ta2_zz_xxzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 242);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxz_x_1,   \
                             ta2_zz_xx_x_0,   \
                             ta2_zz_xx_x_1,   \
                             ta2_zz_xxz_x_0,  \
                             ta2_zz_xxz_x_1,  \
                             ta2_zz_xxzz_x_0, \
                             ta2_zz_xxzz_y_0, \
                             ta2_zz_xxzz_z_0, \
                             ta2_zz_xzz_y_0,  \
                             ta2_zz_xzz_y_1,  \
                             ta2_zz_xzz_z_0,  \
                             ta2_zz_xzz_z_1,  \
                             ta2_zz_zz_y_0,   \
                             ta2_zz_zz_y_1,   \
                             ta2_zz_zz_z_0,   \
                             ta2_zz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxzz_x_0[i] =
            ta2_zz_xx_x_0[i] * fe_0 - ta2_zz_xx_x_1[i] * fe_0 + 2.0 * ta1_z_xxz_x_1[i] + ta2_zz_xxz_x_0[i] * pa_z[i] - ta2_zz_xxz_x_1[i] * pc_z[i];

        ta2_zz_xxzz_y_0[i] = ta2_zz_zz_y_0[i] * fe_0 - ta2_zz_zz_y_1[i] * fe_0 + ta2_zz_xzz_y_0[i] * pa_x[i] - ta2_zz_xzz_y_1[i] * pc_x[i];

        ta2_zz_xxzz_z_0[i] = ta2_zz_zz_z_0[i] * fe_0 - ta2_zz_zz_z_1[i] * fe_0 + ta2_zz_xzz_z_0[i] * pa_x[i] - ta2_zz_xzz_z_1[i] * pc_x[i];
    }

    // Set up 243-246 components of targeted buffer : GP

    auto ta2_zz_xyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 243);

    auto ta2_zz_xyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 244);

    auto ta2_zz_xyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 245);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_zz_xyyy_x_0, \
                             ta2_zz_xyyy_y_0, \
                             ta2_zz_xyyy_z_0, \
                             ta2_zz_yyy_0_0,  \
                             ta2_zz_yyy_0_1,  \
                             ta2_zz_yyy_x_0,  \
                             ta2_zz_yyy_x_1,  \
                             ta2_zz_yyy_y_0,  \
                             ta2_zz_yyy_y_1,  \
                             ta2_zz_yyy_z_0,  \
                             ta2_zz_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyy_x_0[i] = ta2_zz_yyy_0_0[i] * fe_0 - ta2_zz_yyy_0_1[i] * fe_0 + ta2_zz_yyy_x_0[i] * pa_x[i] - ta2_zz_yyy_x_1[i] * pc_x[i];

        ta2_zz_xyyy_y_0[i] = ta2_zz_yyy_y_0[i] * pa_x[i] - ta2_zz_yyy_y_1[i] * pc_x[i];

        ta2_zz_xyyy_z_0[i] = ta2_zz_yyy_z_0[i] * pa_x[i] - ta2_zz_yyy_z_1[i] * pc_x[i];
    }

    // Set up 246-249 components of targeted buffer : GP

    auto ta2_zz_xyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 246);

    auto ta2_zz_xyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 247);

    auto ta2_zz_xyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 248);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xyy_x_1,   \
                             ta2_zz_xyy_x_0,  \
                             ta2_zz_xyy_x_1,  \
                             ta2_zz_xyyz_x_0, \
                             ta2_zz_xyyz_y_0, \
                             ta2_zz_xyyz_z_0, \
                             ta2_zz_yyz_y_0,  \
                             ta2_zz_yyz_y_1,  \
                             ta2_zz_yyz_z_0,  \
                             ta2_zz_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xyyz_x_0[i] = 2.0 * ta1_z_xyy_x_1[i] + ta2_zz_xyy_x_0[i] * pa_z[i] - ta2_zz_xyy_x_1[i] * pc_z[i];

        ta2_zz_xyyz_y_0[i] = ta2_zz_yyz_y_0[i] * pa_x[i] - ta2_zz_yyz_y_1[i] * pc_x[i];

        ta2_zz_xyyz_z_0[i] = ta2_zz_yyz_z_0[i] * pa_x[i] - ta2_zz_yyz_z_1[i] * pc_x[i];
    }

    // Set up 249-252 components of targeted buffer : GP

    auto ta2_zz_xyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 249);

    auto ta2_zz_xyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 250);

    auto ta2_zz_xyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 251);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta2_zz_xyzz_x_0, \
                             ta2_zz_xyzz_y_0, \
                             ta2_zz_xyzz_z_0, \
                             ta2_zz_xzz_x_0,  \
                             ta2_zz_xzz_x_1,  \
                             ta2_zz_yzz_y_0,  \
                             ta2_zz_yzz_y_1,  \
                             ta2_zz_yzz_z_0,  \
                             ta2_zz_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xyzz_x_0[i] = ta2_zz_xzz_x_0[i] * pa_y[i] - ta2_zz_xzz_x_1[i] * pc_y[i];

        ta2_zz_xyzz_y_0[i] = ta2_zz_yzz_y_0[i] * pa_x[i] - ta2_zz_yzz_y_1[i] * pc_x[i];

        ta2_zz_xyzz_z_0[i] = ta2_zz_yzz_z_0[i] * pa_x[i] - ta2_zz_yzz_z_1[i] * pc_x[i];
    }

    // Set up 252-255 components of targeted buffer : GP

    auto ta2_zz_xzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 252);

    auto ta2_zz_xzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 253);

    auto ta2_zz_xzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 254);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_zz_xzzz_x_0, \
                             ta2_zz_xzzz_y_0, \
                             ta2_zz_xzzz_z_0, \
                             ta2_zz_zzz_0_0,  \
                             ta2_zz_zzz_0_1,  \
                             ta2_zz_zzz_x_0,  \
                             ta2_zz_zzz_x_1,  \
                             ta2_zz_zzz_y_0,  \
                             ta2_zz_zzz_y_1,  \
                             ta2_zz_zzz_z_0,  \
                             ta2_zz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzzz_x_0[i] = ta2_zz_zzz_0_0[i] * fe_0 - ta2_zz_zzz_0_1[i] * fe_0 + ta2_zz_zzz_x_0[i] * pa_x[i] - ta2_zz_zzz_x_1[i] * pc_x[i];

        ta2_zz_xzzz_y_0[i] = ta2_zz_zzz_y_0[i] * pa_x[i] - ta2_zz_zzz_y_1[i] * pc_x[i];

        ta2_zz_xzzz_z_0[i] = ta2_zz_zzz_z_0[i] * pa_x[i] - ta2_zz_zzz_z_1[i] * pc_x[i];
    }

    // Set up 255-258 components of targeted buffer : GP

    auto ta2_zz_yyyy_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 255);

    auto ta2_zz_yyyy_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 256);

    auto ta2_zz_yyyy_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 257);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_zz_yy_x_0,   \
                             ta2_zz_yy_x_1,   \
                             ta2_zz_yy_y_0,   \
                             ta2_zz_yy_y_1,   \
                             ta2_zz_yy_z_0,   \
                             ta2_zz_yy_z_1,   \
                             ta2_zz_yyy_0_0,  \
                             ta2_zz_yyy_0_1,  \
                             ta2_zz_yyy_x_0,  \
                             ta2_zz_yyy_x_1,  \
                             ta2_zz_yyy_y_0,  \
                             ta2_zz_yyy_y_1,  \
                             ta2_zz_yyy_z_0,  \
                             ta2_zz_yyy_z_1,  \
                             ta2_zz_yyyy_x_0, \
                             ta2_zz_yyyy_y_0, \
                             ta2_zz_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyy_x_0[i] =
            3.0 * ta2_zz_yy_x_0[i] * fe_0 - 3.0 * ta2_zz_yy_x_1[i] * fe_0 + ta2_zz_yyy_x_0[i] * pa_y[i] - ta2_zz_yyy_x_1[i] * pc_y[i];

        ta2_zz_yyyy_y_0[i] = 3.0 * ta2_zz_yy_y_0[i] * fe_0 - 3.0 * ta2_zz_yy_y_1[i] * fe_0 + ta2_zz_yyy_0_0[i] * fe_0 - ta2_zz_yyy_0_1[i] * fe_0 +
                             ta2_zz_yyy_y_0[i] * pa_y[i] - ta2_zz_yyy_y_1[i] * pc_y[i];

        ta2_zz_yyyy_z_0[i] =
            3.0 * ta2_zz_yy_z_0[i] * fe_0 - 3.0 * ta2_zz_yy_z_1[i] * fe_0 + ta2_zz_yyy_z_0[i] * pa_y[i] - ta2_zz_yyy_z_1[i] * pc_y[i];
    }

    // Set up 258-261 components of targeted buffer : GP

    auto ta2_zz_yyyz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 258);

    auto ta2_zz_yyyz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 259);

    auto ta2_zz_yyyz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 260);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyy_x_1,   \
                             ta1_z_yyy_y_1,   \
                             ta2_zz_yyy_x_0,  \
                             ta2_zz_yyy_x_1,  \
                             ta2_zz_yyy_y_0,  \
                             ta2_zz_yyy_y_1,  \
                             ta2_zz_yyyz_x_0, \
                             ta2_zz_yyyz_y_0, \
                             ta2_zz_yyyz_z_0, \
                             ta2_zz_yyz_z_0,  \
                             ta2_zz_yyz_z_1,  \
                             ta2_zz_yz_z_0,   \
                             ta2_zz_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyz_x_0[i] = 2.0 * ta1_z_yyy_x_1[i] + ta2_zz_yyy_x_0[i] * pa_z[i] - ta2_zz_yyy_x_1[i] * pc_z[i];

        ta2_zz_yyyz_y_0[i] = 2.0 * ta1_z_yyy_y_1[i] + ta2_zz_yyy_y_0[i] * pa_z[i] - ta2_zz_yyy_y_1[i] * pc_z[i];

        ta2_zz_yyyz_z_0[i] =
            2.0 * ta2_zz_yz_z_0[i] * fe_0 - 2.0 * ta2_zz_yz_z_1[i] * fe_0 + ta2_zz_yyz_z_0[i] * pa_y[i] - ta2_zz_yyz_z_1[i] * pc_y[i];
    }

    // Set up 261-264 components of targeted buffer : GP

    auto ta2_zz_yyzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 261);

    auto ta2_zz_yyzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 262);

    auto ta2_zz_yyzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 263);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyz_y_1,   \
                             ta2_zz_yy_y_0,   \
                             ta2_zz_yy_y_1,   \
                             ta2_zz_yyz_y_0,  \
                             ta2_zz_yyz_y_1,  \
                             ta2_zz_yyzz_x_0, \
                             ta2_zz_yyzz_y_0, \
                             ta2_zz_yyzz_z_0, \
                             ta2_zz_yzz_x_0,  \
                             ta2_zz_yzz_x_1,  \
                             ta2_zz_yzz_z_0,  \
                             ta2_zz_yzz_z_1,  \
                             ta2_zz_zz_x_0,   \
                             ta2_zz_zz_x_1,   \
                             ta2_zz_zz_z_0,   \
                             ta2_zz_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyzz_x_0[i] = ta2_zz_zz_x_0[i] * fe_0 - ta2_zz_zz_x_1[i] * fe_0 + ta2_zz_yzz_x_0[i] * pa_y[i] - ta2_zz_yzz_x_1[i] * pc_y[i];

        ta2_zz_yyzz_y_0[i] =
            ta2_zz_yy_y_0[i] * fe_0 - ta2_zz_yy_y_1[i] * fe_0 + 2.0 * ta1_z_yyz_y_1[i] + ta2_zz_yyz_y_0[i] * pa_z[i] - ta2_zz_yyz_y_1[i] * pc_z[i];

        ta2_zz_yyzz_z_0[i] = ta2_zz_zz_z_0[i] * fe_0 - ta2_zz_zz_z_1[i] * fe_0 + ta2_zz_yzz_z_0[i] * pa_y[i] - ta2_zz_yzz_z_1[i] * pc_y[i];
    }

    // Set up 264-267 components of targeted buffer : GP

    auto ta2_zz_yzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 264);

    auto ta2_zz_yzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 265);

    auto ta2_zz_yzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 266);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_zz_yzzz_x_0, \
                             ta2_zz_yzzz_y_0, \
                             ta2_zz_yzzz_z_0, \
                             ta2_zz_zzz_0_0,  \
                             ta2_zz_zzz_0_1,  \
                             ta2_zz_zzz_x_0,  \
                             ta2_zz_zzz_x_1,  \
                             ta2_zz_zzz_y_0,  \
                             ta2_zz_zzz_y_1,  \
                             ta2_zz_zzz_z_0,  \
                             ta2_zz_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzzz_x_0[i] = ta2_zz_zzz_x_0[i] * pa_y[i] - ta2_zz_zzz_x_1[i] * pc_y[i];

        ta2_zz_yzzz_y_0[i] = ta2_zz_zzz_0_0[i] * fe_0 - ta2_zz_zzz_0_1[i] * fe_0 + ta2_zz_zzz_y_0[i] * pa_y[i] - ta2_zz_zzz_y_1[i] * pc_y[i];

        ta2_zz_yzzz_z_0[i] = ta2_zz_zzz_z_0[i] * pa_y[i] - ta2_zz_zzz_z_1[i] * pc_y[i];
    }

    // Set up 267-270 components of targeted buffer : GP

    auto ta2_zz_zzzz_x_0 = pbuffer.data(idx_npot_geom_020_0_gp + 267);

    auto ta2_zz_zzzz_y_0 = pbuffer.data(idx_npot_geom_020_0_gp + 268);

    auto ta2_zz_zzzz_z_0 = pbuffer.data(idx_npot_geom_020_0_gp + 269);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_z_1,   \
                             ta2_zz_zz_x_0,   \
                             ta2_zz_zz_x_1,   \
                             ta2_zz_zz_y_0,   \
                             ta2_zz_zz_y_1,   \
                             ta2_zz_zz_z_0,   \
                             ta2_zz_zz_z_1,   \
                             ta2_zz_zzz_0_0,  \
                             ta2_zz_zzz_0_1,  \
                             ta2_zz_zzz_x_0,  \
                             ta2_zz_zzz_x_1,  \
                             ta2_zz_zzz_y_0,  \
                             ta2_zz_zzz_y_1,  \
                             ta2_zz_zzz_z_0,  \
                             ta2_zz_zzz_z_1,  \
                             ta2_zz_zzzz_x_0, \
                             ta2_zz_zzzz_y_0, \
                             ta2_zz_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzzz_x_0[i] = 3.0 * ta2_zz_zz_x_0[i] * fe_0 - 3.0 * ta2_zz_zz_x_1[i] * fe_0 + 2.0 * ta1_z_zzz_x_1[i] + ta2_zz_zzz_x_0[i] * pa_z[i] -
                             ta2_zz_zzz_x_1[i] * pc_z[i];

        ta2_zz_zzzz_y_0[i] = 3.0 * ta2_zz_zz_y_0[i] * fe_0 - 3.0 * ta2_zz_zz_y_1[i] * fe_0 + 2.0 * ta1_z_zzz_y_1[i] + ta2_zz_zzz_y_0[i] * pa_z[i] -
                             ta2_zz_zzz_y_1[i] * pc_z[i];

        ta2_zz_zzzz_z_0[i] = 3.0 * ta2_zz_zz_z_0[i] * fe_0 - 3.0 * ta2_zz_zz_z_1[i] * fe_0 + ta2_zz_zzz_0_0[i] * fe_0 - ta2_zz_zzz_0_1[i] * fe_0 +
                             2.0 * ta1_z_zzz_z_1[i] + ta2_zz_zzz_z_0[i] * pa_z[i] - ta2_zz_zzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
