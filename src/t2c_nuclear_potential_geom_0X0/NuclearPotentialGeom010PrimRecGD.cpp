#include "NuclearPotentialGeom010PrimRecGD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gd,
                                        const size_t              idx_npot_geom_010_0_dd,
                                        const size_t              idx_npot_geom_010_1_dd,
                                        const size_t              idx_npot_geom_010_0_fp,
                                        const size_t              idx_npot_geom_010_1_fp,
                                        const size_t              idx_npot_1_fd,
                                        const size_t              idx_npot_geom_010_0_fd,
                                        const size_t              idx_npot_geom_010_1_fd,
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

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd);

    auto ta1_x_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 1);

    auto ta1_x_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 2);

    auto ta1_x_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 3);

    auto ta1_x_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 4);

    auto ta1_x_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 5);

    auto ta1_x_xy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 6);

    auto ta1_x_xy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 8);

    auto ta1_x_xz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 12);

    auto ta1_x_xz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 13);

    auto ta1_x_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 14);

    auto ta1_x_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 18);

    auto ta1_x_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 19);

    auto ta1_x_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 20);

    auto ta1_x_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 21);

    auto ta1_x_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 22);

    auto ta1_x_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 23);

    auto ta1_x_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 26);

    auto ta1_x_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 29);

    auto ta1_x_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 30);

    auto ta1_x_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 31);

    auto ta1_x_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 32);

    auto ta1_x_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 33);

    auto ta1_x_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 34);

    auto ta1_x_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 35);

    auto ta1_y_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 36);

    auto ta1_y_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 37);

    auto ta1_y_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 38);

    auto ta1_y_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 39);

    auto ta1_y_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 40);

    auto ta1_y_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 41);

    auto ta1_y_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 45);

    auto ta1_y_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 46);

    auto ta1_y_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 52);

    auto ta1_y_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 53);

    auto ta1_y_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 54);

    auto ta1_y_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 55);

    auto ta1_y_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 56);

    auto ta1_y_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 57);

    auto ta1_y_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 58);

    auto ta1_y_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 59);

    auto ta1_y_yz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 61);

    auto ta1_y_yz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 63);

    auto ta1_y_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 64);

    auto ta1_y_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 66);

    auto ta1_y_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 67);

    auto ta1_y_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 68);

    auto ta1_y_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 69);

    auto ta1_y_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 70);

    auto ta1_y_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 71);

    auto ta1_z_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 72);

    auto ta1_z_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 73);

    auto ta1_z_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 74);

    auto ta1_z_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 75);

    auto ta1_z_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 76);

    auto ta1_z_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 77);

    auto ta1_z_xy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 81);

    auto ta1_z_xy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 82);

    auto ta1_z_xz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 88);

    auto ta1_z_xz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 89);

    auto ta1_z_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 90);

    auto ta1_z_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 91);

    auto ta1_z_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 92);

    auto ta1_z_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 93);

    auto ta1_z_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 94);

    auto ta1_z_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 95);

    auto ta1_z_yz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 98);

    auto ta1_z_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 100);

    auto ta1_z_yz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 101);

    auto ta1_z_zz_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 102);

    auto ta1_z_zz_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 103);

    auto ta1_z_zz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 104);

    auto ta1_z_zz_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 105);

    auto ta1_z_zz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 106);

    auto ta1_z_zz_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 107);

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd);

    auto ta1_x_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 1);

    auto ta1_x_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 2);

    auto ta1_x_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 3);

    auto ta1_x_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 4);

    auto ta1_x_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 5);

    auto ta1_x_xy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 6);

    auto ta1_x_xy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 8);

    auto ta1_x_xz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 12);

    auto ta1_x_xz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 13);

    auto ta1_x_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 14);

    auto ta1_x_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 18);

    auto ta1_x_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 19);

    auto ta1_x_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 20);

    auto ta1_x_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 21);

    auto ta1_x_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 22);

    auto ta1_x_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 23);

    auto ta1_x_yz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 26);

    auto ta1_x_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 29);

    auto ta1_x_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 30);

    auto ta1_x_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 31);

    auto ta1_x_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 32);

    auto ta1_x_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 33);

    auto ta1_x_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 34);

    auto ta1_x_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 35);

    auto ta1_y_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 36);

    auto ta1_y_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 37);

    auto ta1_y_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 38);

    auto ta1_y_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 39);

    auto ta1_y_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 40);

    auto ta1_y_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 41);

    auto ta1_y_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 45);

    auto ta1_y_xy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 46);

    auto ta1_y_xz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 52);

    auto ta1_y_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 53);

    auto ta1_y_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 54);

    auto ta1_y_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 55);

    auto ta1_y_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 56);

    auto ta1_y_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 57);

    auto ta1_y_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 58);

    auto ta1_y_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 59);

    auto ta1_y_yz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 61);

    auto ta1_y_yz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 63);

    auto ta1_y_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 64);

    auto ta1_y_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 66);

    auto ta1_y_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 67);

    auto ta1_y_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 68);

    auto ta1_y_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 69);

    auto ta1_y_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 70);

    auto ta1_y_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 71);

    auto ta1_z_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 72);

    auto ta1_z_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 73);

    auto ta1_z_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 74);

    auto ta1_z_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 75);

    auto ta1_z_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 76);

    auto ta1_z_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 77);

    auto ta1_z_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 81);

    auto ta1_z_xy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 82);

    auto ta1_z_xz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 88);

    auto ta1_z_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 89);

    auto ta1_z_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 90);

    auto ta1_z_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 91);

    auto ta1_z_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 92);

    auto ta1_z_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 93);

    auto ta1_z_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 94);

    auto ta1_z_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 95);

    auto ta1_z_yz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 98);

    auto ta1_z_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 100);

    auto ta1_z_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 101);

    auto ta1_z_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 102);

    auto ta1_z_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 103);

    auto ta1_z_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 104);

    auto ta1_z_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 105);

    auto ta1_z_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 106);

    auto ta1_z_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 107);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp);

    auto ta1_x_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 1);

    auto ta1_x_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 2);

    auto ta1_x_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 3);

    auto ta1_x_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 6);

    auto ta1_x_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 8);

    auto ta1_x_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 15);

    auto ta1_x_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 18);

    auto ta1_x_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 19);

    auto ta1_x_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 20);

    auto ta1_x_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 26);

    auto ta1_x_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 27);

    auto ta1_x_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 28);

    auto ta1_x_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 29);

    auto ta1_y_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 30);

    auto ta1_y_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 31);

    auto ta1_y_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 32);

    auto ta1_y_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 40);

    auto ta1_y_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 47);

    auto ta1_y_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 48);

    auto ta1_y_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 49);

    auto ta1_y_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 50);

    auto ta1_y_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 52);

    auto ta1_y_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 53);

    auto ta1_y_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 55);

    auto ta1_y_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 57);

    auto ta1_y_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 58);

    auto ta1_y_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 59);

    auto ta1_z_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 60);

    auto ta1_z_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 61);

    auto ta1_z_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 62);

    auto ta1_z_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 70);

    auto ta1_z_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 77);

    auto ta1_z_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 78);

    auto ta1_z_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 79);

    auto ta1_z_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 80);

    auto ta1_z_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 83);

    auto ta1_z_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 85);

    auto ta1_z_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 86);

    auto ta1_z_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 87);

    auto ta1_z_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 88);

    auto ta1_z_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 89);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp);

    auto ta1_x_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 1);

    auto ta1_x_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 2);

    auto ta1_x_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 3);

    auto ta1_x_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 6);

    auto ta1_x_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 8);

    auto ta1_x_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 15);

    auto ta1_x_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 18);

    auto ta1_x_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 19);

    auto ta1_x_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 20);

    auto ta1_x_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 26);

    auto ta1_x_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 27);

    auto ta1_x_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 28);

    auto ta1_x_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 29);

    auto ta1_y_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 30);

    auto ta1_y_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 31);

    auto ta1_y_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 32);

    auto ta1_y_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 40);

    auto ta1_y_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 47);

    auto ta1_y_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 48);

    auto ta1_y_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 49);

    auto ta1_y_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 50);

    auto ta1_y_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 52);

    auto ta1_y_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 53);

    auto ta1_y_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 55);

    auto ta1_y_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 57);

    auto ta1_y_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 58);

    auto ta1_y_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 59);

    auto ta1_z_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 60);

    auto ta1_z_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 61);

    auto ta1_z_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 62);

    auto ta1_z_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 70);

    auto ta1_z_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 77);

    auto ta1_z_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 78);

    auto ta1_z_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 79);

    auto ta1_z_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 80);

    auto ta1_z_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 83);

    auto ta1_z_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 85);

    auto ta1_z_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 86);

    auto ta1_z_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 87);

    auto ta1_z_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 88);

    auto ta1_z_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 89);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_1 = pbuffer.data(idx_npot_1_fd);

    auto ta_xxx_xy_1 = pbuffer.data(idx_npot_1_fd + 1);

    auto ta_xxx_xz_1 = pbuffer.data(idx_npot_1_fd + 2);

    auto ta_xxx_yy_1 = pbuffer.data(idx_npot_1_fd + 3);

    auto ta_xxx_yz_1 = pbuffer.data(idx_npot_1_fd + 4);

    auto ta_xxx_zz_1 = pbuffer.data(idx_npot_1_fd + 5);

    auto ta_xxy_xx_1 = pbuffer.data(idx_npot_1_fd + 6);

    auto ta_xxy_xy_1 = pbuffer.data(idx_npot_1_fd + 7);

    auto ta_xxy_xz_1 = pbuffer.data(idx_npot_1_fd + 8);

    auto ta_xxy_yy_1 = pbuffer.data(idx_npot_1_fd + 9);

    auto ta_xxz_xx_1 = pbuffer.data(idx_npot_1_fd + 12);

    auto ta_xxz_xy_1 = pbuffer.data(idx_npot_1_fd + 13);

    auto ta_xxz_xz_1 = pbuffer.data(idx_npot_1_fd + 14);

    auto ta_xxz_zz_1 = pbuffer.data(idx_npot_1_fd + 17);

    auto ta_xyy_xx_1 = pbuffer.data(idx_npot_1_fd + 18);

    auto ta_xyy_xy_1 = pbuffer.data(idx_npot_1_fd + 19);

    auto ta_xyy_yy_1 = pbuffer.data(idx_npot_1_fd + 21);

    auto ta_xyy_yz_1 = pbuffer.data(idx_npot_1_fd + 22);

    auto ta_xzz_xx_1 = pbuffer.data(idx_npot_1_fd + 30);

    auto ta_xzz_xz_1 = pbuffer.data(idx_npot_1_fd + 32);

    auto ta_xzz_yz_1 = pbuffer.data(idx_npot_1_fd + 34);

    auto ta_xzz_zz_1 = pbuffer.data(idx_npot_1_fd + 35);

    auto ta_yyy_xx_1 = pbuffer.data(idx_npot_1_fd + 36);

    auto ta_yyy_xy_1 = pbuffer.data(idx_npot_1_fd + 37);

    auto ta_yyy_xz_1 = pbuffer.data(idx_npot_1_fd + 38);

    auto ta_yyy_yy_1 = pbuffer.data(idx_npot_1_fd + 39);

    auto ta_yyy_yz_1 = pbuffer.data(idx_npot_1_fd + 40);

    auto ta_yyy_zz_1 = pbuffer.data(idx_npot_1_fd + 41);

    auto ta_yyz_xy_1 = pbuffer.data(idx_npot_1_fd + 43);

    auto ta_yyz_yy_1 = pbuffer.data(idx_npot_1_fd + 45);

    auto ta_yyz_yz_1 = pbuffer.data(idx_npot_1_fd + 46);

    auto ta_yyz_zz_1 = pbuffer.data(idx_npot_1_fd + 47);

    auto ta_yzz_xz_1 = pbuffer.data(idx_npot_1_fd + 50);

    auto ta_yzz_yy_1 = pbuffer.data(idx_npot_1_fd + 51);

    auto ta_yzz_yz_1 = pbuffer.data(idx_npot_1_fd + 52);

    auto ta_yzz_zz_1 = pbuffer.data(idx_npot_1_fd + 53);

    auto ta_zzz_xx_1 = pbuffer.data(idx_npot_1_fd + 54);

    auto ta_zzz_xy_1 = pbuffer.data(idx_npot_1_fd + 55);

    auto ta_zzz_xz_1 = pbuffer.data(idx_npot_1_fd + 56);

    auto ta_zzz_yy_1 = pbuffer.data(idx_npot_1_fd + 57);

    auto ta_zzz_yz_1 = pbuffer.data(idx_npot_1_fd + 58);

    auto ta_zzz_zz_1 = pbuffer.data(idx_npot_1_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto ta1_x_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd);

    auto ta1_x_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 1);

    auto ta1_x_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 2);

    auto ta1_x_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 3);

    auto ta1_x_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 4);

    auto ta1_x_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 5);

    auto ta1_x_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 6);

    auto ta1_x_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 7);

    auto ta1_x_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 8);

    auto ta1_x_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 9);

    auto ta1_x_xxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 11);

    auto ta1_x_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 12);

    auto ta1_x_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 13);

    auto ta1_x_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 14);

    auto ta1_x_xxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 15);

    auto ta1_x_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 16);

    auto ta1_x_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 17);

    auto ta1_x_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 18);

    auto ta1_x_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 19);

    auto ta1_x_xyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 20);

    auto ta1_x_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 21);

    auto ta1_x_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 22);

    auto ta1_x_xyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 26);

    auto ta1_x_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 30);

    auto ta1_x_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 31);

    auto ta1_x_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 32);

    auto ta1_x_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 34);

    auto ta1_x_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 35);

    auto ta1_x_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 36);

    auto ta1_x_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 37);

    auto ta1_x_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 38);

    auto ta1_x_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 39);

    auto ta1_x_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 40);

    auto ta1_x_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 41);

    auto ta1_x_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 43);

    auto ta1_x_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 44);

    auto ta1_x_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 45);

    auto ta1_x_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 46);

    auto ta1_x_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 47);

    auto ta1_x_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 48);

    auto ta1_x_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 50);

    auto ta1_x_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 51);

    auto ta1_x_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 52);

    auto ta1_x_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 53);

    auto ta1_x_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 54);

    auto ta1_x_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 55);

    auto ta1_x_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 56);

    auto ta1_x_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 57);

    auto ta1_x_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 58);

    auto ta1_x_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 59);

    auto ta1_y_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 60);

    auto ta1_y_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 61);

    auto ta1_y_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 62);

    auto ta1_y_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 63);

    auto ta1_y_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 64);

    auto ta1_y_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 65);

    auto ta1_y_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 66);

    auto ta1_y_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 67);

    auto ta1_y_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 68);

    auto ta1_y_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 69);

    auto ta1_y_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 70);

    auto ta1_y_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 72);

    auto ta1_y_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 73);

    auto ta1_y_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 74);

    auto ta1_y_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 76);

    auto ta1_y_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 77);

    auto ta1_y_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 78);

    auto ta1_y_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 79);

    auto ta1_y_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 81);

    auto ta1_y_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 82);

    auto ta1_y_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 83);

    auto ta1_y_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 88);

    auto ta1_y_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 90);

    auto ta1_y_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 92);

    auto ta1_y_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 93);

    auto ta1_y_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 94);

    auto ta1_y_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 95);

    auto ta1_y_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 96);

    auto ta1_y_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 97);

    auto ta1_y_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 98);

    auto ta1_y_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 99);

    auto ta1_y_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 100);

    auto ta1_y_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 101);

    auto ta1_y_yyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 102);

    auto ta1_y_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 103);

    auto ta1_y_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 104);

    auto ta1_y_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 105);

    auto ta1_y_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 106);

    auto ta1_y_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 107);

    auto ta1_y_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 109);

    auto ta1_y_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 110);

    auto ta1_y_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 111);

    auto ta1_y_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 112);

    auto ta1_y_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 113);

    auto ta1_y_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 114);

    auto ta1_y_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 115);

    auto ta1_y_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 116);

    auto ta1_y_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 117);

    auto ta1_y_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 118);

    auto ta1_y_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 119);

    auto ta1_z_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 120);

    auto ta1_z_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 121);

    auto ta1_z_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 122);

    auto ta1_z_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 123);

    auto ta1_z_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 124);

    auto ta1_z_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 125);

    auto ta1_z_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 126);

    auto ta1_z_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 127);

    auto ta1_z_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 128);

    auto ta1_z_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 129);

    auto ta1_z_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 130);

    auto ta1_z_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 132);

    auto ta1_z_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 133);

    auto ta1_z_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 134);

    auto ta1_z_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 136);

    auto ta1_z_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 137);

    auto ta1_z_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 138);

    auto ta1_z_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 139);

    auto ta1_z_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 141);

    auto ta1_z_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 142);

    auto ta1_z_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 143);

    auto ta1_z_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 148);

    auto ta1_z_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 150);

    auto ta1_z_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 152);

    auto ta1_z_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 153);

    auto ta1_z_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 154);

    auto ta1_z_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 155);

    auto ta1_z_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 156);

    auto ta1_z_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 157);

    auto ta1_z_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 158);

    auto ta1_z_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 159);

    auto ta1_z_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 160);

    auto ta1_z_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 161);

    auto ta1_z_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 163);

    auto ta1_z_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 164);

    auto ta1_z_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 165);

    auto ta1_z_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 166);

    auto ta1_z_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 167);

    auto ta1_z_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 168);

    auto ta1_z_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 169);

    auto ta1_z_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 170);

    auto ta1_z_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 171);

    auto ta1_z_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 172);

    auto ta1_z_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 173);

    auto ta1_z_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 174);

    auto ta1_z_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 175);

    auto ta1_z_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 176);

    auto ta1_z_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 177);

    auto ta1_z_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 178);

    auto ta1_z_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 179);

    // Set up components of auxiliary buffer : FD

    auto ta1_x_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd);

    auto ta1_x_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 1);

    auto ta1_x_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 2);

    auto ta1_x_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 3);

    auto ta1_x_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 4);

    auto ta1_x_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 5);

    auto ta1_x_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 6);

    auto ta1_x_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 7);

    auto ta1_x_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 8);

    auto ta1_x_xxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 9);

    auto ta1_x_xxy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 11);

    auto ta1_x_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 12);

    auto ta1_x_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 13);

    auto ta1_x_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 14);

    auto ta1_x_xxz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 15);

    auto ta1_x_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 16);

    auto ta1_x_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 17);

    auto ta1_x_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 18);

    auto ta1_x_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 19);

    auto ta1_x_xyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 20);

    auto ta1_x_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 21);

    auto ta1_x_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 22);

    auto ta1_x_xyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 26);

    auto ta1_x_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 30);

    auto ta1_x_xzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 31);

    auto ta1_x_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 32);

    auto ta1_x_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 34);

    auto ta1_x_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 35);

    auto ta1_x_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 36);

    auto ta1_x_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 37);

    auto ta1_x_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 38);

    auto ta1_x_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 39);

    auto ta1_x_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 40);

    auto ta1_x_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 41);

    auto ta1_x_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 43);

    auto ta1_x_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 44);

    auto ta1_x_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 45);

    auto ta1_x_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 46);

    auto ta1_x_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 47);

    auto ta1_x_yzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 48);

    auto ta1_x_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 50);

    auto ta1_x_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 51);

    auto ta1_x_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 52);

    auto ta1_x_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 53);

    auto ta1_x_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 54);

    auto ta1_x_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 55);

    auto ta1_x_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 56);

    auto ta1_x_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 57);

    auto ta1_x_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 58);

    auto ta1_x_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 59);

    auto ta1_y_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 60);

    auto ta1_y_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 61);

    auto ta1_y_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 62);

    auto ta1_y_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 63);

    auto ta1_y_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 64);

    auto ta1_y_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 65);

    auto ta1_y_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 66);

    auto ta1_y_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 67);

    auto ta1_y_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 68);

    auto ta1_y_xxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 69);

    auto ta1_y_xxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 70);

    auto ta1_y_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 72);

    auto ta1_y_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 73);

    auto ta1_y_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 74);

    auto ta1_y_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 76);

    auto ta1_y_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 77);

    auto ta1_y_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 78);

    auto ta1_y_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 79);

    auto ta1_y_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 81);

    auto ta1_y_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 82);

    auto ta1_y_xyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 83);

    auto ta1_y_xyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 88);

    auto ta1_y_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 90);

    auto ta1_y_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 92);

    auto ta1_y_xzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 93);

    auto ta1_y_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 94);

    auto ta1_y_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 95);

    auto ta1_y_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 96);

    auto ta1_y_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 97);

    auto ta1_y_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 98);

    auto ta1_y_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 99);

    auto ta1_y_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 100);

    auto ta1_y_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 101);

    auto ta1_y_yyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 102);

    auto ta1_y_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 103);

    auto ta1_y_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 104);

    auto ta1_y_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 105);

    auto ta1_y_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 106);

    auto ta1_y_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 107);

    auto ta1_y_yzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 109);

    auto ta1_y_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 110);

    auto ta1_y_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 111);

    auto ta1_y_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 112);

    auto ta1_y_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 113);

    auto ta1_y_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 114);

    auto ta1_y_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 115);

    auto ta1_y_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 116);

    auto ta1_y_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 117);

    auto ta1_y_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 118);

    auto ta1_y_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 119);

    auto ta1_z_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 120);

    auto ta1_z_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 121);

    auto ta1_z_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 122);

    auto ta1_z_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 123);

    auto ta1_z_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 124);

    auto ta1_z_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 125);

    auto ta1_z_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 126);

    auto ta1_z_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 127);

    auto ta1_z_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 128);

    auto ta1_z_xxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 129);

    auto ta1_z_xxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 130);

    auto ta1_z_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 132);

    auto ta1_z_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 133);

    auto ta1_z_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 134);

    auto ta1_z_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 136);

    auto ta1_z_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 137);

    auto ta1_z_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 138);

    auto ta1_z_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 139);

    auto ta1_z_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 141);

    auto ta1_z_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 142);

    auto ta1_z_xyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 143);

    auto ta1_z_xyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 148);

    auto ta1_z_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 150);

    auto ta1_z_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 152);

    auto ta1_z_xzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 153);

    auto ta1_z_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 154);

    auto ta1_z_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 155);

    auto ta1_z_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 156);

    auto ta1_z_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 157);

    auto ta1_z_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 158);

    auto ta1_z_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 159);

    auto ta1_z_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 160);

    auto ta1_z_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 161);

    auto ta1_z_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 163);

    auto ta1_z_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 164);

    auto ta1_z_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 165);

    auto ta1_z_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 166);

    auto ta1_z_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 167);

    auto ta1_z_yzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 168);

    auto ta1_z_yzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 169);

    auto ta1_z_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 170);

    auto ta1_z_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 171);

    auto ta1_z_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 172);

    auto ta1_z_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 173);

    auto ta1_z_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 174);

    auto ta1_z_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 175);

    auto ta1_z_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 176);

    auto ta1_z_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 177);

    auto ta1_z_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 178);

    auto ta1_z_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 179);

    // Set up 0-6 components of targeted buffer : GD

    auto ta1_x_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd);

    auto ta1_x_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 1);

    auto ta1_x_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 2);

    auto ta1_x_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 3);

    auto ta1_x_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 4);

    auto ta1_x_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 5);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_yy_0,   \
                             ta1_x_xx_yy_1,   \
                             ta1_x_xx_yz_0,   \
                             ta1_x_xx_yz_1,   \
                             ta1_x_xx_zz_0,   \
                             ta1_x_xx_zz_1,   \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_xx_0,  \
                             ta1_x_xxx_xx_1,  \
                             ta1_x_xxx_xy_0,  \
                             ta1_x_xxx_xy_1,  \
                             ta1_x_xxx_xz_0,  \
                             ta1_x_xxx_xz_1,  \
                             ta1_x_xxx_y_0,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxx_yy_0,  \
                             ta1_x_xxx_yy_1,  \
                             ta1_x_xxx_yz_0,  \
                             ta1_x_xxx_yz_1,  \
                             ta1_x_xxx_z_0,   \
                             ta1_x_xxx_z_1,   \
                             ta1_x_xxx_zz_0,  \
                             ta1_x_xxx_zz_1,  \
                             ta1_x_xxxx_xx_0, \
                             ta1_x_xxxx_xy_0, \
                             ta1_x_xxxx_xz_0, \
                             ta1_x_xxxx_yy_0, \
                             ta1_x_xxxx_yz_0, \
                             ta1_x_xxxx_zz_0, \
                             ta_xxx_xx_1,     \
                             ta_xxx_xy_1,     \
                             ta_xxx_xz_1,     \
                             ta_xxx_yy_1,     \
                             ta_xxx_yz_1,     \
                             ta_xxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_xx_0[i] = 3.0 * ta1_x_xx_xx_0[i] * fe_0 - 3.0 * ta1_x_xx_xx_1[i] * fe_0 +
                             2.0 * ta1_x_xxx_x_0[i] * fe_0 - 2.0 * ta1_x_xxx_x_1[i] * fe_0 + ta_xxx_xx_1[i] +
                             ta1_x_xxx_xx_0[i] * pa_x[i] - ta1_x_xxx_xx_1[i] * pc_x[i];

        ta1_x_xxxx_xy_0[i] = 3.0 * ta1_x_xx_xy_0[i] * fe_0 - 3.0 * ta1_x_xx_xy_1[i] * fe_0 + ta1_x_xxx_y_0[i] * fe_0 -
                             ta1_x_xxx_y_1[i] * fe_0 + ta_xxx_xy_1[i] + ta1_x_xxx_xy_0[i] * pa_x[i] -
                             ta1_x_xxx_xy_1[i] * pc_x[i];

        ta1_x_xxxx_xz_0[i] = 3.0 * ta1_x_xx_xz_0[i] * fe_0 - 3.0 * ta1_x_xx_xz_1[i] * fe_0 + ta1_x_xxx_z_0[i] * fe_0 -
                             ta1_x_xxx_z_1[i] * fe_0 + ta_xxx_xz_1[i] + ta1_x_xxx_xz_0[i] * pa_x[i] -
                             ta1_x_xxx_xz_1[i] * pc_x[i];

        ta1_x_xxxx_yy_0[i] = 3.0 * ta1_x_xx_yy_0[i] * fe_0 - 3.0 * ta1_x_xx_yy_1[i] * fe_0 + ta_xxx_yy_1[i] +
                             ta1_x_xxx_yy_0[i] * pa_x[i] - ta1_x_xxx_yy_1[i] * pc_x[i];

        ta1_x_xxxx_yz_0[i] = 3.0 * ta1_x_xx_yz_0[i] * fe_0 - 3.0 * ta1_x_xx_yz_1[i] * fe_0 + ta_xxx_yz_1[i] +
                             ta1_x_xxx_yz_0[i] * pa_x[i] - ta1_x_xxx_yz_1[i] * pc_x[i];

        ta1_x_xxxx_zz_0[i] = 3.0 * ta1_x_xx_zz_0[i] * fe_0 - 3.0 * ta1_x_xx_zz_1[i] * fe_0 + ta_xxx_zz_1[i] +
                             ta1_x_xxx_zz_0[i] * pa_x[i] - ta1_x_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto ta1_x_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 6);

    auto ta1_x_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 7);

    auto ta1_x_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 8);

    auto ta1_x_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 9);

    auto ta1_x_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 10);

    auto ta1_x_xxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 11);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_xx_0,  \
                             ta1_x_xxx_xx_1,  \
                             ta1_x_xxx_xy_0,  \
                             ta1_x_xxx_xy_1,  \
                             ta1_x_xxx_xz_0,  \
                             ta1_x_xxx_xz_1,  \
                             ta1_x_xxx_y_0,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxx_yy_0,  \
                             ta1_x_xxx_yy_1,  \
                             ta1_x_xxx_yz_0,  \
                             ta1_x_xxx_yz_1,  \
                             ta1_x_xxx_z_0,   \
                             ta1_x_xxx_z_1,   \
                             ta1_x_xxx_zz_0,  \
                             ta1_x_xxx_zz_1,  \
                             ta1_x_xxxy_xx_0, \
                             ta1_x_xxxy_xy_0, \
                             ta1_x_xxxy_xz_0, \
                             ta1_x_xxxy_yy_0, \
                             ta1_x_xxxy_yz_0, \
                             ta1_x_xxxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_xx_0[i] = ta1_x_xxx_xx_0[i] * pa_y[i] - ta1_x_xxx_xx_1[i] * pc_y[i];

        ta1_x_xxxy_xy_0[i] = ta1_x_xxx_x_0[i] * fe_0 - ta1_x_xxx_x_1[i] * fe_0 + ta1_x_xxx_xy_0[i] * pa_y[i] -
                             ta1_x_xxx_xy_1[i] * pc_y[i];

        ta1_x_xxxy_xz_0[i] = ta1_x_xxx_xz_0[i] * pa_y[i] - ta1_x_xxx_xz_1[i] * pc_y[i];

        ta1_x_xxxy_yy_0[i] = 2.0 * ta1_x_xxx_y_0[i] * fe_0 - 2.0 * ta1_x_xxx_y_1[i] * fe_0 +
                             ta1_x_xxx_yy_0[i] * pa_y[i] - ta1_x_xxx_yy_1[i] * pc_y[i];

        ta1_x_xxxy_yz_0[i] = ta1_x_xxx_z_0[i] * fe_0 - ta1_x_xxx_z_1[i] * fe_0 + ta1_x_xxx_yz_0[i] * pa_y[i] -
                             ta1_x_xxx_yz_1[i] * pc_y[i];

        ta1_x_xxxy_zz_0[i] = ta1_x_xxx_zz_0[i] * pa_y[i] - ta1_x_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto ta1_x_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 12);

    auto ta1_x_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 13);

    auto ta1_x_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 14);

    auto ta1_x_xxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 15);

    auto ta1_x_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 16);

    auto ta1_x_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 17);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_xx_0,  \
                             ta1_x_xxx_xx_1,  \
                             ta1_x_xxx_xy_0,  \
                             ta1_x_xxx_xy_1,  \
                             ta1_x_xxx_xz_0,  \
                             ta1_x_xxx_xz_1,  \
                             ta1_x_xxx_y_0,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxx_yy_0,  \
                             ta1_x_xxx_yy_1,  \
                             ta1_x_xxx_yz_0,  \
                             ta1_x_xxx_yz_1,  \
                             ta1_x_xxx_z_0,   \
                             ta1_x_xxx_z_1,   \
                             ta1_x_xxx_zz_0,  \
                             ta1_x_xxx_zz_1,  \
                             ta1_x_xxxz_xx_0, \
                             ta1_x_xxxz_xy_0, \
                             ta1_x_xxxz_xz_0, \
                             ta1_x_xxxz_yy_0, \
                             ta1_x_xxxz_yz_0, \
                             ta1_x_xxxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_xx_0[i] = ta1_x_xxx_xx_0[i] * pa_z[i] - ta1_x_xxx_xx_1[i] * pc_z[i];

        ta1_x_xxxz_xy_0[i] = ta1_x_xxx_xy_0[i] * pa_z[i] - ta1_x_xxx_xy_1[i] * pc_z[i];

        ta1_x_xxxz_xz_0[i] = ta1_x_xxx_x_0[i] * fe_0 - ta1_x_xxx_x_1[i] * fe_0 + ta1_x_xxx_xz_0[i] * pa_z[i] -
                             ta1_x_xxx_xz_1[i] * pc_z[i];

        ta1_x_xxxz_yy_0[i] = ta1_x_xxx_yy_0[i] * pa_z[i] - ta1_x_xxx_yy_1[i] * pc_z[i];

        ta1_x_xxxz_yz_0[i] = ta1_x_xxx_y_0[i] * fe_0 - ta1_x_xxx_y_1[i] * fe_0 + ta1_x_xxx_yz_0[i] * pa_z[i] -
                             ta1_x_xxx_yz_1[i] * pc_z[i];

        ta1_x_xxxz_zz_0[i] = 2.0 * ta1_x_xxx_z_0[i] * fe_0 - 2.0 * ta1_x_xxx_z_1[i] * fe_0 +
                             ta1_x_xxx_zz_0[i] * pa_z[i] - ta1_x_xxx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto ta1_x_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 18);

    auto ta1_x_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 19);

    auto ta1_x_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 20);

    auto ta1_x_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 21);

    auto ta1_x_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 22);

    auto ta1_x_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 23);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_zz_0,   \
                             ta1_x_xx_zz_1,   \
                             ta1_x_xxy_x_0,   \
                             ta1_x_xxy_x_1,   \
                             ta1_x_xxy_xx_0,  \
                             ta1_x_xxy_xx_1,  \
                             ta1_x_xxy_xy_0,  \
                             ta1_x_xxy_xy_1,  \
                             ta1_x_xxy_xz_0,  \
                             ta1_x_xxy_xz_1,  \
                             ta1_x_xxy_zz_0,  \
                             ta1_x_xxy_zz_1,  \
                             ta1_x_xxyy_xx_0, \
                             ta1_x_xxyy_xy_0, \
                             ta1_x_xxyy_xz_0, \
                             ta1_x_xxyy_yy_0, \
                             ta1_x_xxyy_yz_0, \
                             ta1_x_xxyy_zz_0, \
                             ta1_x_xyy_yy_0,  \
                             ta1_x_xyy_yy_1,  \
                             ta1_x_xyy_yz_0,  \
                             ta1_x_xyy_yz_1,  \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yy_yz_0,   \
                             ta1_x_yy_yz_1,   \
                             ta_xyy_yy_1,     \
                             ta_xyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_xx_0[i] = ta1_x_xx_xx_0[i] * fe_0 - ta1_x_xx_xx_1[i] * fe_0 + ta1_x_xxy_xx_0[i] * pa_y[i] -
                             ta1_x_xxy_xx_1[i] * pc_y[i];

        ta1_x_xxyy_xy_0[i] = ta1_x_xx_xy_0[i] * fe_0 - ta1_x_xx_xy_1[i] * fe_0 + ta1_x_xxy_x_0[i] * fe_0 -
                             ta1_x_xxy_x_1[i] * fe_0 + ta1_x_xxy_xy_0[i] * pa_y[i] - ta1_x_xxy_xy_1[i] * pc_y[i];

        ta1_x_xxyy_xz_0[i] = ta1_x_xx_xz_0[i] * fe_0 - ta1_x_xx_xz_1[i] * fe_0 + ta1_x_xxy_xz_0[i] * pa_y[i] -
                             ta1_x_xxy_xz_1[i] * pc_y[i];

        ta1_x_xxyy_yy_0[i] = ta1_x_yy_yy_0[i] * fe_0 - ta1_x_yy_yy_1[i] * fe_0 + ta_xyy_yy_1[i] +
                             ta1_x_xyy_yy_0[i] * pa_x[i] - ta1_x_xyy_yy_1[i] * pc_x[i];

        ta1_x_xxyy_yz_0[i] = ta1_x_yy_yz_0[i] * fe_0 - ta1_x_yy_yz_1[i] * fe_0 + ta_xyy_yz_1[i] +
                             ta1_x_xyy_yz_0[i] * pa_x[i] - ta1_x_xyy_yz_1[i] * pc_x[i];

        ta1_x_xxyy_zz_0[i] = ta1_x_xx_zz_0[i] * fe_0 - ta1_x_xx_zz_1[i] * fe_0 + ta1_x_xxy_zz_0[i] * pa_y[i] -
                             ta1_x_xxy_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto ta1_x_xxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 24);

    auto ta1_x_xxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 25);

    auto ta1_x_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 26);

    auto ta1_x_xxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 27);

    auto ta1_x_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 28);

    auto ta1_x_xxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 29);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxy_xy_0,  \
                             ta1_x_xxy_xy_1,  \
                             ta1_x_xxy_yy_0,  \
                             ta1_x_xxy_yy_1,  \
                             ta1_x_xxyz_xx_0, \
                             ta1_x_xxyz_xy_0, \
                             ta1_x_xxyz_xz_0, \
                             ta1_x_xxyz_yy_0, \
                             ta1_x_xxyz_yz_0, \
                             ta1_x_xxyz_zz_0, \
                             ta1_x_xxz_xx_0,  \
                             ta1_x_xxz_xx_1,  \
                             ta1_x_xxz_xz_0,  \
                             ta1_x_xxz_xz_1,  \
                             ta1_x_xxz_yz_0,  \
                             ta1_x_xxz_yz_1,  \
                             ta1_x_xxz_z_0,   \
                             ta1_x_xxz_z_1,   \
                             ta1_x_xxz_zz_0,  \
                             ta1_x_xxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyz_xx_0[i] = ta1_x_xxz_xx_0[i] * pa_y[i] - ta1_x_xxz_xx_1[i] * pc_y[i];

        ta1_x_xxyz_xy_0[i] = ta1_x_xxy_xy_0[i] * pa_z[i] - ta1_x_xxy_xy_1[i] * pc_z[i];

        ta1_x_xxyz_xz_0[i] = ta1_x_xxz_xz_0[i] * pa_y[i] - ta1_x_xxz_xz_1[i] * pc_y[i];

        ta1_x_xxyz_yy_0[i] = ta1_x_xxy_yy_0[i] * pa_z[i] - ta1_x_xxy_yy_1[i] * pc_z[i];

        ta1_x_xxyz_yz_0[i] = ta1_x_xxz_z_0[i] * fe_0 - ta1_x_xxz_z_1[i] * fe_0 + ta1_x_xxz_yz_0[i] * pa_y[i] -
                             ta1_x_xxz_yz_1[i] * pc_y[i];

        ta1_x_xxyz_zz_0[i] = ta1_x_xxz_zz_0[i] * pa_y[i] - ta1_x_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto ta1_x_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 30);

    auto ta1_x_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 31);

    auto ta1_x_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 32);

    auto ta1_x_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 33);

    auto ta1_x_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 34);

    auto ta1_x_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 35);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_yy_0,   \
                             ta1_x_xx_yy_1,   \
                             ta1_x_xxz_x_0,   \
                             ta1_x_xxz_x_1,   \
                             ta1_x_xxz_xx_0,  \
                             ta1_x_xxz_xx_1,  \
                             ta1_x_xxz_xy_0,  \
                             ta1_x_xxz_xy_1,  \
                             ta1_x_xxz_xz_0,  \
                             ta1_x_xxz_xz_1,  \
                             ta1_x_xxz_yy_0,  \
                             ta1_x_xxz_yy_1,  \
                             ta1_x_xxzz_xx_0, \
                             ta1_x_xxzz_xy_0, \
                             ta1_x_xxzz_xz_0, \
                             ta1_x_xxzz_yy_0, \
                             ta1_x_xxzz_yz_0, \
                             ta1_x_xxzz_zz_0, \
                             ta1_x_xzz_yz_0,  \
                             ta1_x_xzz_yz_1,  \
                             ta1_x_xzz_zz_0,  \
                             ta1_x_xzz_zz_1,  \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             ta_xzz_yz_1,     \
                             ta_xzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_xx_0[i] = ta1_x_xx_xx_0[i] * fe_0 - ta1_x_xx_xx_1[i] * fe_0 + ta1_x_xxz_xx_0[i] * pa_z[i] -
                             ta1_x_xxz_xx_1[i] * pc_z[i];

        ta1_x_xxzz_xy_0[i] = ta1_x_xx_xy_0[i] * fe_0 - ta1_x_xx_xy_1[i] * fe_0 + ta1_x_xxz_xy_0[i] * pa_z[i] -
                             ta1_x_xxz_xy_1[i] * pc_z[i];

        ta1_x_xxzz_xz_0[i] = ta1_x_xx_xz_0[i] * fe_0 - ta1_x_xx_xz_1[i] * fe_0 + ta1_x_xxz_x_0[i] * fe_0 -
                             ta1_x_xxz_x_1[i] * fe_0 + ta1_x_xxz_xz_0[i] * pa_z[i] - ta1_x_xxz_xz_1[i] * pc_z[i];

        ta1_x_xxzz_yy_0[i] = ta1_x_xx_yy_0[i] * fe_0 - ta1_x_xx_yy_1[i] * fe_0 + ta1_x_xxz_yy_0[i] * pa_z[i] -
                             ta1_x_xxz_yy_1[i] * pc_z[i];

        ta1_x_xxzz_yz_0[i] = ta1_x_zz_yz_0[i] * fe_0 - ta1_x_zz_yz_1[i] * fe_0 + ta_xzz_yz_1[i] +
                             ta1_x_xzz_yz_0[i] * pa_x[i] - ta1_x_xzz_yz_1[i] * pc_x[i];

        ta1_x_xxzz_zz_0[i] = ta1_x_zz_zz_0[i] * fe_0 - ta1_x_zz_zz_1[i] * fe_0 + ta_xzz_zz_1[i] +
                             ta1_x_xzz_zz_0[i] * pa_x[i] - ta1_x_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto ta1_x_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 36);

    auto ta1_x_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 37);

    auto ta1_x_xyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 38);

    auto ta1_x_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 39);

    auto ta1_x_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 40);

    auto ta1_x_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 41);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xy_xx_0,   \
                             ta1_x_xy_xx_1,   \
                             ta1_x_xy_xz_0,   \
                             ta1_x_xy_xz_1,   \
                             ta1_x_xyy_xx_0,  \
                             ta1_x_xyy_xx_1,  \
                             ta1_x_xyy_xz_0,  \
                             ta1_x_xyy_xz_1,  \
                             ta1_x_xyyy_xx_0, \
                             ta1_x_xyyy_xy_0, \
                             ta1_x_xyyy_xz_0, \
                             ta1_x_xyyy_yy_0, \
                             ta1_x_xyyy_yz_0, \
                             ta1_x_xyyy_zz_0, \
                             ta1_x_yyy_xy_0,  \
                             ta1_x_yyy_xy_1,  \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_yy_0,  \
                             ta1_x_yyy_yy_1,  \
                             ta1_x_yyy_yz_0,  \
                             ta1_x_yyy_yz_1,  \
                             ta1_x_yyy_zz_0,  \
                             ta1_x_yyy_zz_1,  \
                             ta_yyy_xy_1,     \
                             ta_yyy_yy_1,     \
                             ta_yyy_yz_1,     \
                             ta_yyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_xx_0[i] = 2.0 * ta1_x_xy_xx_0[i] * fe_0 - 2.0 * ta1_x_xy_xx_1[i] * fe_0 +
                             ta1_x_xyy_xx_0[i] * pa_y[i] - ta1_x_xyy_xx_1[i] * pc_y[i];

        ta1_x_xyyy_xy_0[i] = ta1_x_yyy_y_0[i] * fe_0 - ta1_x_yyy_y_1[i] * fe_0 + ta_yyy_xy_1[i] +
                             ta1_x_yyy_xy_0[i] * pa_x[i] - ta1_x_yyy_xy_1[i] * pc_x[i];

        ta1_x_xyyy_xz_0[i] = 2.0 * ta1_x_xy_xz_0[i] * fe_0 - 2.0 * ta1_x_xy_xz_1[i] * fe_0 +
                             ta1_x_xyy_xz_0[i] * pa_y[i] - ta1_x_xyy_xz_1[i] * pc_y[i];

        ta1_x_xyyy_yy_0[i] = ta_yyy_yy_1[i] + ta1_x_yyy_yy_0[i] * pa_x[i] - ta1_x_yyy_yy_1[i] * pc_x[i];

        ta1_x_xyyy_yz_0[i] = ta_yyy_yz_1[i] + ta1_x_yyy_yz_0[i] * pa_x[i] - ta1_x_yyy_yz_1[i] * pc_x[i];

        ta1_x_xyyy_zz_0[i] = ta_yyy_zz_1[i] + ta1_x_yyy_zz_0[i] * pa_x[i] - ta1_x_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto ta1_x_xyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 42);

    auto ta1_x_xyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 43);

    auto ta1_x_xyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 44);

    auto ta1_x_xyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 45);

    auto ta1_x_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 46);

    auto ta1_x_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 47);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xyy_xx_0,  \
                             ta1_x_xyy_xx_1,  \
                             ta1_x_xyy_xy_0,  \
                             ta1_x_xyy_xy_1,  \
                             ta1_x_xyy_yy_0,  \
                             ta1_x_xyy_yy_1,  \
                             ta1_x_xyyz_xx_0, \
                             ta1_x_xyyz_xy_0, \
                             ta1_x_xyyz_xz_0, \
                             ta1_x_xyyz_yy_0, \
                             ta1_x_xyyz_yz_0, \
                             ta1_x_xyyz_zz_0, \
                             ta1_x_xyz_xz_0,  \
                             ta1_x_xyz_xz_1,  \
                             ta1_x_xz_xz_0,   \
                             ta1_x_xz_xz_1,   \
                             ta1_x_yyz_yz_0,  \
                             ta1_x_yyz_yz_1,  \
                             ta1_x_yyz_zz_0,  \
                             ta1_x_yyz_zz_1,  \
                             ta_yyz_yz_1,     \
                             ta_yyz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyz_xx_0[i] = ta1_x_xyy_xx_0[i] * pa_z[i] - ta1_x_xyy_xx_1[i] * pc_z[i];

        ta1_x_xyyz_xy_0[i] = ta1_x_xyy_xy_0[i] * pa_z[i] - ta1_x_xyy_xy_1[i] * pc_z[i];

        ta1_x_xyyz_xz_0[i] = ta1_x_xz_xz_0[i] * fe_0 - ta1_x_xz_xz_1[i] * fe_0 + ta1_x_xyz_xz_0[i] * pa_y[i] -
                             ta1_x_xyz_xz_1[i] * pc_y[i];

        ta1_x_xyyz_yy_0[i] = ta1_x_xyy_yy_0[i] * pa_z[i] - ta1_x_xyy_yy_1[i] * pc_z[i];

        ta1_x_xyyz_yz_0[i] = ta_yyz_yz_1[i] + ta1_x_yyz_yz_0[i] * pa_x[i] - ta1_x_yyz_yz_1[i] * pc_x[i];

        ta1_x_xyyz_zz_0[i] = ta_yyz_zz_1[i] + ta1_x_yyz_zz_0[i] * pa_x[i] - ta1_x_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto ta1_x_xyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 48);

    auto ta1_x_xyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 49);

    auto ta1_x_xyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 50);

    auto ta1_x_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 51);

    auto ta1_x_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 52);

    auto ta1_x_xyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 53);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xyzz_xx_0, \
                             ta1_x_xyzz_xy_0, \
                             ta1_x_xyzz_xz_0, \
                             ta1_x_xyzz_yy_0, \
                             ta1_x_xyzz_yz_0, \
                             ta1_x_xyzz_zz_0, \
                             ta1_x_xzz_x_0,   \
                             ta1_x_xzz_x_1,   \
                             ta1_x_xzz_xx_0,  \
                             ta1_x_xzz_xx_1,  \
                             ta1_x_xzz_xy_0,  \
                             ta1_x_xzz_xy_1,  \
                             ta1_x_xzz_xz_0,  \
                             ta1_x_xzz_xz_1,  \
                             ta1_x_xzz_zz_0,  \
                             ta1_x_xzz_zz_1,  \
                             ta1_x_yzz_yy_0,  \
                             ta1_x_yzz_yy_1,  \
                             ta1_x_yzz_yz_0,  \
                             ta1_x_yzz_yz_1,  \
                             ta_yzz_yy_1,     \
                             ta_yzz_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzz_xx_0[i] = ta1_x_xzz_xx_0[i] * pa_y[i] - ta1_x_xzz_xx_1[i] * pc_y[i];

        ta1_x_xyzz_xy_0[i] = ta1_x_xzz_x_0[i] * fe_0 - ta1_x_xzz_x_1[i] * fe_0 + ta1_x_xzz_xy_0[i] * pa_y[i] -
                             ta1_x_xzz_xy_1[i] * pc_y[i];

        ta1_x_xyzz_xz_0[i] = ta1_x_xzz_xz_0[i] * pa_y[i] - ta1_x_xzz_xz_1[i] * pc_y[i];

        ta1_x_xyzz_yy_0[i] = ta_yzz_yy_1[i] + ta1_x_yzz_yy_0[i] * pa_x[i] - ta1_x_yzz_yy_1[i] * pc_x[i];

        ta1_x_xyzz_yz_0[i] = ta_yzz_yz_1[i] + ta1_x_yzz_yz_0[i] * pa_x[i] - ta1_x_yzz_yz_1[i] * pc_x[i];

        ta1_x_xyzz_zz_0[i] = ta1_x_xzz_zz_0[i] * pa_y[i] - ta1_x_xzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto ta1_x_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 54);

    auto ta1_x_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 55);

    auto ta1_x_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 56);

    auto ta1_x_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 57);

    auto ta1_x_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 58);

    auto ta1_x_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 59);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xz_xx_0,   \
                             ta1_x_xz_xx_1,   \
                             ta1_x_xz_xy_0,   \
                             ta1_x_xz_xy_1,   \
                             ta1_x_xzz_xx_0,  \
                             ta1_x_xzz_xx_1,  \
                             ta1_x_xzz_xy_0,  \
                             ta1_x_xzz_xy_1,  \
                             ta1_x_xzzz_xx_0, \
                             ta1_x_xzzz_xy_0, \
                             ta1_x_xzzz_xz_0, \
                             ta1_x_xzzz_yy_0, \
                             ta1_x_xzzz_yz_0, \
                             ta1_x_xzzz_zz_0, \
                             ta1_x_zzz_xz_0,  \
                             ta1_x_zzz_xz_1,  \
                             ta1_x_zzz_yy_0,  \
                             ta1_x_zzz_yy_1,  \
                             ta1_x_zzz_yz_0,  \
                             ta1_x_zzz_yz_1,  \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             ta1_x_zzz_zz_0,  \
                             ta1_x_zzz_zz_1,  \
                             ta_zzz_xz_1,     \
                             ta_zzz_yy_1,     \
                             ta_zzz_yz_1,     \
                             ta_zzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_xx_0[i] = 2.0 * ta1_x_xz_xx_0[i] * fe_0 - 2.0 * ta1_x_xz_xx_1[i] * fe_0 +
                             ta1_x_xzz_xx_0[i] * pa_z[i] - ta1_x_xzz_xx_1[i] * pc_z[i];

        ta1_x_xzzz_xy_0[i] = 2.0 * ta1_x_xz_xy_0[i] * fe_0 - 2.0 * ta1_x_xz_xy_1[i] * fe_0 +
                             ta1_x_xzz_xy_0[i] * pa_z[i] - ta1_x_xzz_xy_1[i] * pc_z[i];

        ta1_x_xzzz_xz_0[i] = ta1_x_zzz_z_0[i] * fe_0 - ta1_x_zzz_z_1[i] * fe_0 + ta_zzz_xz_1[i] +
                             ta1_x_zzz_xz_0[i] * pa_x[i] - ta1_x_zzz_xz_1[i] * pc_x[i];

        ta1_x_xzzz_yy_0[i] = ta_zzz_yy_1[i] + ta1_x_zzz_yy_0[i] * pa_x[i] - ta1_x_zzz_yy_1[i] * pc_x[i];

        ta1_x_xzzz_yz_0[i] = ta_zzz_yz_1[i] + ta1_x_zzz_yz_0[i] * pa_x[i] - ta1_x_zzz_yz_1[i] * pc_x[i];

        ta1_x_xzzz_zz_0[i] = ta_zzz_zz_1[i] + ta1_x_zzz_zz_0[i] * pa_x[i] - ta1_x_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto ta1_x_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 60);

    auto ta1_x_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 61);

    auto ta1_x_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 62);

    auto ta1_x_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 63);

    auto ta1_x_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 64);

    auto ta1_x_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 65);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yy_xx_0,   \
                             ta1_x_yy_xx_1,   \
                             ta1_x_yy_xy_0,   \
                             ta1_x_yy_xy_1,   \
                             ta1_x_yy_xz_0,   \
                             ta1_x_yy_xz_1,   \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yy_yz_0,   \
                             ta1_x_yy_yz_1,   \
                             ta1_x_yy_zz_0,   \
                             ta1_x_yy_zz_1,   \
                             ta1_x_yyy_x_0,   \
                             ta1_x_yyy_x_1,   \
                             ta1_x_yyy_xx_0,  \
                             ta1_x_yyy_xx_1,  \
                             ta1_x_yyy_xy_0,  \
                             ta1_x_yyy_xy_1,  \
                             ta1_x_yyy_xz_0,  \
                             ta1_x_yyy_xz_1,  \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_yy_0,  \
                             ta1_x_yyy_yy_1,  \
                             ta1_x_yyy_yz_0,  \
                             ta1_x_yyy_yz_1,  \
                             ta1_x_yyy_z_0,   \
                             ta1_x_yyy_z_1,   \
                             ta1_x_yyy_zz_0,  \
                             ta1_x_yyy_zz_1,  \
                             ta1_x_yyyy_xx_0, \
                             ta1_x_yyyy_xy_0, \
                             ta1_x_yyyy_xz_0, \
                             ta1_x_yyyy_yy_0, \
                             ta1_x_yyyy_yz_0, \
                             ta1_x_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_xx_0[i] = 3.0 * ta1_x_yy_xx_0[i] * fe_0 - 3.0 * ta1_x_yy_xx_1[i] * fe_0 +
                             ta1_x_yyy_xx_0[i] * pa_y[i] - ta1_x_yyy_xx_1[i] * pc_y[i];

        ta1_x_yyyy_xy_0[i] = 3.0 * ta1_x_yy_xy_0[i] * fe_0 - 3.0 * ta1_x_yy_xy_1[i] * fe_0 + ta1_x_yyy_x_0[i] * fe_0 -
                             ta1_x_yyy_x_1[i] * fe_0 + ta1_x_yyy_xy_0[i] * pa_y[i] - ta1_x_yyy_xy_1[i] * pc_y[i];

        ta1_x_yyyy_xz_0[i] = 3.0 * ta1_x_yy_xz_0[i] * fe_0 - 3.0 * ta1_x_yy_xz_1[i] * fe_0 +
                             ta1_x_yyy_xz_0[i] * pa_y[i] - ta1_x_yyy_xz_1[i] * pc_y[i];

        ta1_x_yyyy_yy_0[i] = 3.0 * ta1_x_yy_yy_0[i] * fe_0 - 3.0 * ta1_x_yy_yy_1[i] * fe_0 +
                             2.0 * ta1_x_yyy_y_0[i] * fe_0 - 2.0 * ta1_x_yyy_y_1[i] * fe_0 +
                             ta1_x_yyy_yy_0[i] * pa_y[i] - ta1_x_yyy_yy_1[i] * pc_y[i];

        ta1_x_yyyy_yz_0[i] = 3.0 * ta1_x_yy_yz_0[i] * fe_0 - 3.0 * ta1_x_yy_yz_1[i] * fe_0 + ta1_x_yyy_z_0[i] * fe_0 -
                             ta1_x_yyy_z_1[i] * fe_0 + ta1_x_yyy_yz_0[i] * pa_y[i] - ta1_x_yyy_yz_1[i] * pc_y[i];

        ta1_x_yyyy_zz_0[i] = 3.0 * ta1_x_yy_zz_0[i] * fe_0 - 3.0 * ta1_x_yy_zz_1[i] * fe_0 +
                             ta1_x_yyy_zz_0[i] * pa_y[i] - ta1_x_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto ta1_x_yyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 66);

    auto ta1_x_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 67);

    auto ta1_x_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 68);

    auto ta1_x_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 69);

    auto ta1_x_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 70);

    auto ta1_x_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 71);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyy_xx_0,  \
                             ta1_x_yyy_xx_1,  \
                             ta1_x_yyy_xy_0,  \
                             ta1_x_yyy_xy_1,  \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_yy_0,  \
                             ta1_x_yyy_yy_1,  \
                             ta1_x_yyy_yz_0,  \
                             ta1_x_yyy_yz_1,  \
                             ta1_x_yyyz_xx_0, \
                             ta1_x_yyyz_xy_0, \
                             ta1_x_yyyz_xz_0, \
                             ta1_x_yyyz_yy_0, \
                             ta1_x_yyyz_yz_0, \
                             ta1_x_yyyz_zz_0, \
                             ta1_x_yyz_xz_0,  \
                             ta1_x_yyz_xz_1,  \
                             ta1_x_yyz_zz_0,  \
                             ta1_x_yyz_zz_1,  \
                             ta1_x_yz_xz_0,   \
                             ta1_x_yz_xz_1,   \
                             ta1_x_yz_zz_0,   \
                             ta1_x_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_xx_0[i] = ta1_x_yyy_xx_0[i] * pa_z[i] - ta1_x_yyy_xx_1[i] * pc_z[i];

        ta1_x_yyyz_xy_0[i] = ta1_x_yyy_xy_0[i] * pa_z[i] - ta1_x_yyy_xy_1[i] * pc_z[i];

        ta1_x_yyyz_xz_0[i] = 2.0 * ta1_x_yz_xz_0[i] * fe_0 - 2.0 * ta1_x_yz_xz_1[i] * fe_0 +
                             ta1_x_yyz_xz_0[i] * pa_y[i] - ta1_x_yyz_xz_1[i] * pc_y[i];

        ta1_x_yyyz_yy_0[i] = ta1_x_yyy_yy_0[i] * pa_z[i] - ta1_x_yyy_yy_1[i] * pc_z[i];

        ta1_x_yyyz_yz_0[i] = ta1_x_yyy_y_0[i] * fe_0 - ta1_x_yyy_y_1[i] * fe_0 + ta1_x_yyy_yz_0[i] * pa_z[i] -
                             ta1_x_yyy_yz_1[i] * pc_z[i];

        ta1_x_yyyz_zz_0[i] = 2.0 * ta1_x_yz_zz_0[i] * fe_0 - 2.0 * ta1_x_yz_zz_1[i] * fe_0 +
                             ta1_x_yyz_zz_0[i] * pa_y[i] - ta1_x_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto ta1_x_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 72);

    auto ta1_x_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 73);

    auto ta1_x_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 74);

    auto ta1_x_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 75);

    auto ta1_x_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 76);

    auto ta1_x_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 77);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yy_xy_0,   \
                             ta1_x_yy_xy_1,   \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yyz_xy_0,  \
                             ta1_x_yyz_xy_1,  \
                             ta1_x_yyz_yy_0,  \
                             ta1_x_yyz_yy_1,  \
                             ta1_x_yyzz_xx_0, \
                             ta1_x_yyzz_xy_0, \
                             ta1_x_yyzz_xz_0, \
                             ta1_x_yyzz_yy_0, \
                             ta1_x_yyzz_yz_0, \
                             ta1_x_yyzz_zz_0, \
                             ta1_x_yzz_xx_0,  \
                             ta1_x_yzz_xx_1,  \
                             ta1_x_yzz_xz_0,  \
                             ta1_x_yzz_xz_1,  \
                             ta1_x_yzz_yz_0,  \
                             ta1_x_yzz_yz_1,  \
                             ta1_x_yzz_z_0,   \
                             ta1_x_yzz_z_1,   \
                             ta1_x_yzz_zz_0,  \
                             ta1_x_yzz_zz_1,  \
                             ta1_x_zz_xx_0,   \
                             ta1_x_zz_xx_1,   \
                             ta1_x_zz_xz_0,   \
                             ta1_x_zz_xz_1,   \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_xx_0[i] = ta1_x_zz_xx_0[i] * fe_0 - ta1_x_zz_xx_1[i] * fe_0 + ta1_x_yzz_xx_0[i] * pa_y[i] -
                             ta1_x_yzz_xx_1[i] * pc_y[i];

        ta1_x_yyzz_xy_0[i] = ta1_x_yy_xy_0[i] * fe_0 - ta1_x_yy_xy_1[i] * fe_0 + ta1_x_yyz_xy_0[i] * pa_z[i] -
                             ta1_x_yyz_xy_1[i] * pc_z[i];

        ta1_x_yyzz_xz_0[i] = ta1_x_zz_xz_0[i] * fe_0 - ta1_x_zz_xz_1[i] * fe_0 + ta1_x_yzz_xz_0[i] * pa_y[i] -
                             ta1_x_yzz_xz_1[i] * pc_y[i];

        ta1_x_yyzz_yy_0[i] = ta1_x_yy_yy_0[i] * fe_0 - ta1_x_yy_yy_1[i] * fe_0 + ta1_x_yyz_yy_0[i] * pa_z[i] -
                             ta1_x_yyz_yy_1[i] * pc_z[i];

        ta1_x_yyzz_yz_0[i] = ta1_x_zz_yz_0[i] * fe_0 - ta1_x_zz_yz_1[i] * fe_0 + ta1_x_yzz_z_0[i] * fe_0 -
                             ta1_x_yzz_z_1[i] * fe_0 + ta1_x_yzz_yz_0[i] * pa_y[i] - ta1_x_yzz_yz_1[i] * pc_y[i];

        ta1_x_yyzz_zz_0[i] = ta1_x_zz_zz_0[i] * fe_0 - ta1_x_zz_zz_1[i] * fe_0 + ta1_x_yzz_zz_0[i] * pa_y[i] -
                             ta1_x_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto ta1_x_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 78);

    auto ta1_x_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 79);

    auto ta1_x_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 80);

    auto ta1_x_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 81);

    auto ta1_x_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 82);

    auto ta1_x_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 83);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yzzz_xx_0, \
                             ta1_x_yzzz_xy_0, \
                             ta1_x_yzzz_xz_0, \
                             ta1_x_yzzz_yy_0, \
                             ta1_x_yzzz_yz_0, \
                             ta1_x_yzzz_zz_0, \
                             ta1_x_zzz_x_0,   \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_xx_0,  \
                             ta1_x_zzz_xx_1,  \
                             ta1_x_zzz_xy_0,  \
                             ta1_x_zzz_xy_1,  \
                             ta1_x_zzz_xz_0,  \
                             ta1_x_zzz_xz_1,  \
                             ta1_x_zzz_y_0,   \
                             ta1_x_zzz_y_1,   \
                             ta1_x_zzz_yy_0,  \
                             ta1_x_zzz_yy_1,  \
                             ta1_x_zzz_yz_0,  \
                             ta1_x_zzz_yz_1,  \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             ta1_x_zzz_zz_0,  \
                             ta1_x_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_xx_0[i] = ta1_x_zzz_xx_0[i] * pa_y[i] - ta1_x_zzz_xx_1[i] * pc_y[i];

        ta1_x_yzzz_xy_0[i] = ta1_x_zzz_x_0[i] * fe_0 - ta1_x_zzz_x_1[i] * fe_0 + ta1_x_zzz_xy_0[i] * pa_y[i] -
                             ta1_x_zzz_xy_1[i] * pc_y[i];

        ta1_x_yzzz_xz_0[i] = ta1_x_zzz_xz_0[i] * pa_y[i] - ta1_x_zzz_xz_1[i] * pc_y[i];

        ta1_x_yzzz_yy_0[i] = 2.0 * ta1_x_zzz_y_0[i] * fe_0 - 2.0 * ta1_x_zzz_y_1[i] * fe_0 +
                             ta1_x_zzz_yy_0[i] * pa_y[i] - ta1_x_zzz_yy_1[i] * pc_y[i];

        ta1_x_yzzz_yz_0[i] = ta1_x_zzz_z_0[i] * fe_0 - ta1_x_zzz_z_1[i] * fe_0 + ta1_x_zzz_yz_0[i] * pa_y[i] -
                             ta1_x_zzz_yz_1[i] * pc_y[i];

        ta1_x_yzzz_zz_0[i] = ta1_x_zzz_zz_0[i] * pa_y[i] - ta1_x_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto ta1_x_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 84);

    auto ta1_x_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 85);

    auto ta1_x_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 86);

    auto ta1_x_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 87);

    auto ta1_x_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 88);

    auto ta1_x_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 89);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_zz_xx_0,   \
                             ta1_x_zz_xx_1,   \
                             ta1_x_zz_xy_0,   \
                             ta1_x_zz_xy_1,   \
                             ta1_x_zz_xz_0,   \
                             ta1_x_zz_xz_1,   \
                             ta1_x_zz_yy_0,   \
                             ta1_x_zz_yy_1,   \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             ta1_x_zzz_x_0,   \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_xx_0,  \
                             ta1_x_zzz_xx_1,  \
                             ta1_x_zzz_xy_0,  \
                             ta1_x_zzz_xy_1,  \
                             ta1_x_zzz_xz_0,  \
                             ta1_x_zzz_xz_1,  \
                             ta1_x_zzz_y_0,   \
                             ta1_x_zzz_y_1,   \
                             ta1_x_zzz_yy_0,  \
                             ta1_x_zzz_yy_1,  \
                             ta1_x_zzz_yz_0,  \
                             ta1_x_zzz_yz_1,  \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             ta1_x_zzz_zz_0,  \
                             ta1_x_zzz_zz_1,  \
                             ta1_x_zzzz_xx_0, \
                             ta1_x_zzzz_xy_0, \
                             ta1_x_zzzz_xz_0, \
                             ta1_x_zzzz_yy_0, \
                             ta1_x_zzzz_yz_0, \
                             ta1_x_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_xx_0[i] = 3.0 * ta1_x_zz_xx_0[i] * fe_0 - 3.0 * ta1_x_zz_xx_1[i] * fe_0 +
                             ta1_x_zzz_xx_0[i] * pa_z[i] - ta1_x_zzz_xx_1[i] * pc_z[i];

        ta1_x_zzzz_xy_0[i] = 3.0 * ta1_x_zz_xy_0[i] * fe_0 - 3.0 * ta1_x_zz_xy_1[i] * fe_0 +
                             ta1_x_zzz_xy_0[i] * pa_z[i] - ta1_x_zzz_xy_1[i] * pc_z[i];

        ta1_x_zzzz_xz_0[i] = 3.0 * ta1_x_zz_xz_0[i] * fe_0 - 3.0 * ta1_x_zz_xz_1[i] * fe_0 + ta1_x_zzz_x_0[i] * fe_0 -
                             ta1_x_zzz_x_1[i] * fe_0 + ta1_x_zzz_xz_0[i] * pa_z[i] - ta1_x_zzz_xz_1[i] * pc_z[i];

        ta1_x_zzzz_yy_0[i] = 3.0 * ta1_x_zz_yy_0[i] * fe_0 - 3.0 * ta1_x_zz_yy_1[i] * fe_0 +
                             ta1_x_zzz_yy_0[i] * pa_z[i] - ta1_x_zzz_yy_1[i] * pc_z[i];

        ta1_x_zzzz_yz_0[i] = 3.0 * ta1_x_zz_yz_0[i] * fe_0 - 3.0 * ta1_x_zz_yz_1[i] * fe_0 + ta1_x_zzz_y_0[i] * fe_0 -
                             ta1_x_zzz_y_1[i] * fe_0 + ta1_x_zzz_yz_0[i] * pa_z[i] - ta1_x_zzz_yz_1[i] * pc_z[i];

        ta1_x_zzzz_zz_0[i] = 3.0 * ta1_x_zz_zz_0[i] * fe_0 - 3.0 * ta1_x_zz_zz_1[i] * fe_0 +
                             2.0 * ta1_x_zzz_z_0[i] * fe_0 - 2.0 * ta1_x_zzz_z_1[i] * fe_0 +
                             ta1_x_zzz_zz_0[i] * pa_z[i] - ta1_x_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 90-96 components of targeted buffer : GD

    auto ta1_y_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 90);

    auto ta1_y_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 91);

    auto ta1_y_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 92);

    auto ta1_y_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 93);

    auto ta1_y_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 94);

    auto ta1_y_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 95);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xy_0,   \
                             ta1_y_xx_xy_1,   \
                             ta1_y_xx_xz_0,   \
                             ta1_y_xx_xz_1,   \
                             ta1_y_xx_yy_0,   \
                             ta1_y_xx_yy_1,   \
                             ta1_y_xx_yz_0,   \
                             ta1_y_xx_yz_1,   \
                             ta1_y_xx_zz_0,   \
                             ta1_y_xx_zz_1,   \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_xx_0,  \
                             ta1_y_xxx_xx_1,  \
                             ta1_y_xxx_xy_0,  \
                             ta1_y_xxx_xy_1,  \
                             ta1_y_xxx_xz_0,  \
                             ta1_y_xxx_xz_1,  \
                             ta1_y_xxx_y_0,   \
                             ta1_y_xxx_y_1,   \
                             ta1_y_xxx_yy_0,  \
                             ta1_y_xxx_yy_1,  \
                             ta1_y_xxx_yz_0,  \
                             ta1_y_xxx_yz_1,  \
                             ta1_y_xxx_z_0,   \
                             ta1_y_xxx_z_1,   \
                             ta1_y_xxx_zz_0,  \
                             ta1_y_xxx_zz_1,  \
                             ta1_y_xxxx_xx_0, \
                             ta1_y_xxxx_xy_0, \
                             ta1_y_xxxx_xz_0, \
                             ta1_y_xxxx_yy_0, \
                             ta1_y_xxxx_yz_0, \
                             ta1_y_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_xx_0[i] = 3.0 * ta1_y_xx_xx_0[i] * fe_0 - 3.0 * ta1_y_xx_xx_1[i] * fe_0 +
                             2.0 * ta1_y_xxx_x_0[i] * fe_0 - 2.0 * ta1_y_xxx_x_1[i] * fe_0 +
                             ta1_y_xxx_xx_0[i] * pa_x[i] - ta1_y_xxx_xx_1[i] * pc_x[i];

        ta1_y_xxxx_xy_0[i] = 3.0 * ta1_y_xx_xy_0[i] * fe_0 - 3.0 * ta1_y_xx_xy_1[i] * fe_0 + ta1_y_xxx_y_0[i] * fe_0 -
                             ta1_y_xxx_y_1[i] * fe_0 + ta1_y_xxx_xy_0[i] * pa_x[i] - ta1_y_xxx_xy_1[i] * pc_x[i];

        ta1_y_xxxx_xz_0[i] = 3.0 * ta1_y_xx_xz_0[i] * fe_0 - 3.0 * ta1_y_xx_xz_1[i] * fe_0 + ta1_y_xxx_z_0[i] * fe_0 -
                             ta1_y_xxx_z_1[i] * fe_0 + ta1_y_xxx_xz_0[i] * pa_x[i] - ta1_y_xxx_xz_1[i] * pc_x[i];

        ta1_y_xxxx_yy_0[i] = 3.0 * ta1_y_xx_yy_0[i] * fe_0 - 3.0 * ta1_y_xx_yy_1[i] * fe_0 +
                             ta1_y_xxx_yy_0[i] * pa_x[i] - ta1_y_xxx_yy_1[i] * pc_x[i];

        ta1_y_xxxx_yz_0[i] = 3.0 * ta1_y_xx_yz_0[i] * fe_0 - 3.0 * ta1_y_xx_yz_1[i] * fe_0 +
                             ta1_y_xxx_yz_0[i] * pa_x[i] - ta1_y_xxx_yz_1[i] * pc_x[i];

        ta1_y_xxxx_zz_0[i] = 3.0 * ta1_y_xx_zz_0[i] * fe_0 - 3.0 * ta1_y_xx_zz_1[i] * fe_0 +
                             ta1_y_xxx_zz_0[i] * pa_x[i] - ta1_y_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : GD

    auto ta1_y_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 96);

    auto ta1_y_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 97);

    auto ta1_y_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 98);

    auto ta1_y_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 99);

    auto ta1_y_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 100);

    auto ta1_y_xxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 101);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_xx_0,  \
                             ta1_y_xxx_xx_1,  \
                             ta1_y_xxx_xy_0,  \
                             ta1_y_xxx_xy_1,  \
                             ta1_y_xxx_xz_0,  \
                             ta1_y_xxx_xz_1,  \
                             ta1_y_xxx_zz_0,  \
                             ta1_y_xxx_zz_1,  \
                             ta1_y_xxxy_xx_0, \
                             ta1_y_xxxy_xy_0, \
                             ta1_y_xxxy_xz_0, \
                             ta1_y_xxxy_yy_0, \
                             ta1_y_xxxy_yz_0, \
                             ta1_y_xxxy_zz_0, \
                             ta1_y_xxy_yy_0,  \
                             ta1_y_xxy_yy_1,  \
                             ta1_y_xxy_yz_0,  \
                             ta1_y_xxy_yz_1,  \
                             ta1_y_xy_yy_0,   \
                             ta1_y_xy_yy_1,   \
                             ta1_y_xy_yz_0,   \
                             ta1_y_xy_yz_1,   \
                             ta_xxx_xx_1,     \
                             ta_xxx_xy_1,     \
                             ta_xxx_xz_1,     \
                             ta_xxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_xx_0[i] = ta_xxx_xx_1[i] + ta1_y_xxx_xx_0[i] * pa_y[i] - ta1_y_xxx_xx_1[i] * pc_y[i];

        ta1_y_xxxy_xy_0[i] = ta1_y_xxx_x_0[i] * fe_0 - ta1_y_xxx_x_1[i] * fe_0 + ta_xxx_xy_1[i] +
                             ta1_y_xxx_xy_0[i] * pa_y[i] - ta1_y_xxx_xy_1[i] * pc_y[i];

        ta1_y_xxxy_xz_0[i] = ta_xxx_xz_1[i] + ta1_y_xxx_xz_0[i] * pa_y[i] - ta1_y_xxx_xz_1[i] * pc_y[i];

        ta1_y_xxxy_yy_0[i] = 2.0 * ta1_y_xy_yy_0[i] * fe_0 - 2.0 * ta1_y_xy_yy_1[i] * fe_0 +
                             ta1_y_xxy_yy_0[i] * pa_x[i] - ta1_y_xxy_yy_1[i] * pc_x[i];

        ta1_y_xxxy_yz_0[i] = 2.0 * ta1_y_xy_yz_0[i] * fe_0 - 2.0 * ta1_y_xy_yz_1[i] * fe_0 +
                             ta1_y_xxy_yz_0[i] * pa_x[i] - ta1_y_xxy_yz_1[i] * pc_x[i];

        ta1_y_xxxy_zz_0[i] = ta_xxx_zz_1[i] + ta1_y_xxx_zz_0[i] * pa_y[i] - ta1_y_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : GD

    auto ta1_y_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 102);

    auto ta1_y_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 103);

    auto ta1_y_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 104);

    auto ta1_y_xxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 105);

    auto ta1_y_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 106);

    auto ta1_y_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 107);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_xx_0,  \
                             ta1_y_xxx_xx_1,  \
                             ta1_y_xxx_xy_0,  \
                             ta1_y_xxx_xy_1,  \
                             ta1_y_xxx_xz_0,  \
                             ta1_y_xxx_xz_1,  \
                             ta1_y_xxx_yy_0,  \
                             ta1_y_xxx_yy_1,  \
                             ta1_y_xxxz_xx_0, \
                             ta1_y_xxxz_xy_0, \
                             ta1_y_xxxz_xz_0, \
                             ta1_y_xxxz_yy_0, \
                             ta1_y_xxxz_yz_0, \
                             ta1_y_xxxz_zz_0, \
                             ta1_y_xxz_yz_0,  \
                             ta1_y_xxz_yz_1,  \
                             ta1_y_xxz_zz_0,  \
                             ta1_y_xxz_zz_1,  \
                             ta1_y_xz_yz_0,   \
                             ta1_y_xz_yz_1,   \
                             ta1_y_xz_zz_0,   \
                             ta1_y_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_xx_0[i] = ta1_y_xxx_xx_0[i] * pa_z[i] - ta1_y_xxx_xx_1[i] * pc_z[i];

        ta1_y_xxxz_xy_0[i] = ta1_y_xxx_xy_0[i] * pa_z[i] - ta1_y_xxx_xy_1[i] * pc_z[i];

        ta1_y_xxxz_xz_0[i] = ta1_y_xxx_x_0[i] * fe_0 - ta1_y_xxx_x_1[i] * fe_0 + ta1_y_xxx_xz_0[i] * pa_z[i] -
                             ta1_y_xxx_xz_1[i] * pc_z[i];

        ta1_y_xxxz_yy_0[i] = ta1_y_xxx_yy_0[i] * pa_z[i] - ta1_y_xxx_yy_1[i] * pc_z[i];

        ta1_y_xxxz_yz_0[i] = 2.0 * ta1_y_xz_yz_0[i] * fe_0 - 2.0 * ta1_y_xz_yz_1[i] * fe_0 +
                             ta1_y_xxz_yz_0[i] * pa_x[i] - ta1_y_xxz_yz_1[i] * pc_x[i];

        ta1_y_xxxz_zz_0[i] = 2.0 * ta1_y_xz_zz_0[i] * fe_0 - 2.0 * ta1_y_xz_zz_1[i] * fe_0 +
                             ta1_y_xxz_zz_0[i] * pa_x[i] - ta1_y_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 108-114 components of targeted buffer : GD

    auto ta1_y_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 108);

    auto ta1_y_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 109);

    auto ta1_y_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 110);

    auto ta1_y_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 111);

    auto ta1_y_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 112);

    auto ta1_y_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 113);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xz_0,   \
                             ta1_y_xx_xz_1,   \
                             ta1_y_xxy_xx_0,  \
                             ta1_y_xxy_xx_1,  \
                             ta1_y_xxy_xz_0,  \
                             ta1_y_xxy_xz_1,  \
                             ta1_y_xxyy_xx_0, \
                             ta1_y_xxyy_xy_0, \
                             ta1_y_xxyy_xz_0, \
                             ta1_y_xxyy_yy_0, \
                             ta1_y_xxyy_yz_0, \
                             ta1_y_xxyy_zz_0, \
                             ta1_y_xyy_xy_0,  \
                             ta1_y_xyy_xy_1,  \
                             ta1_y_xyy_y_0,   \
                             ta1_y_xyy_y_1,   \
                             ta1_y_xyy_yy_0,  \
                             ta1_y_xyy_yy_1,  \
                             ta1_y_xyy_yz_0,  \
                             ta1_y_xyy_yz_1,  \
                             ta1_y_xyy_zz_0,  \
                             ta1_y_xyy_zz_1,  \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yy_zz_0,   \
                             ta1_y_yy_zz_1,   \
                             ta_xxy_xx_1,     \
                             ta_xxy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_xx_0[i] = ta1_y_xx_xx_0[i] * fe_0 - ta1_y_xx_xx_1[i] * fe_0 + ta_xxy_xx_1[i] +
                             ta1_y_xxy_xx_0[i] * pa_y[i] - ta1_y_xxy_xx_1[i] * pc_y[i];

        ta1_y_xxyy_xy_0[i] = ta1_y_yy_xy_0[i] * fe_0 - ta1_y_yy_xy_1[i] * fe_0 + ta1_y_xyy_y_0[i] * fe_0 -
                             ta1_y_xyy_y_1[i] * fe_0 + ta1_y_xyy_xy_0[i] * pa_x[i] - ta1_y_xyy_xy_1[i] * pc_x[i];

        ta1_y_xxyy_xz_0[i] = ta1_y_xx_xz_0[i] * fe_0 - ta1_y_xx_xz_1[i] * fe_0 + ta_xxy_xz_1[i] +
                             ta1_y_xxy_xz_0[i] * pa_y[i] - ta1_y_xxy_xz_1[i] * pc_y[i];

        ta1_y_xxyy_yy_0[i] = ta1_y_yy_yy_0[i] * fe_0 - ta1_y_yy_yy_1[i] * fe_0 + ta1_y_xyy_yy_0[i] * pa_x[i] -
                             ta1_y_xyy_yy_1[i] * pc_x[i];

        ta1_y_xxyy_yz_0[i] = ta1_y_yy_yz_0[i] * fe_0 - ta1_y_yy_yz_1[i] * fe_0 + ta1_y_xyy_yz_0[i] * pa_x[i] -
                             ta1_y_xyy_yz_1[i] * pc_x[i];

        ta1_y_xxyy_zz_0[i] = ta1_y_yy_zz_0[i] * fe_0 - ta1_y_yy_zz_1[i] * fe_0 + ta1_y_xyy_zz_0[i] * pa_x[i] -
                             ta1_y_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 114-120 components of targeted buffer : GD

    auto ta1_y_xxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 114);

    auto ta1_y_xxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 115);

    auto ta1_y_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 116);

    auto ta1_y_xxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 117);

    auto ta1_y_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 118);

    auto ta1_y_xxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 119);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_xxy_xx_0,  \
                             ta1_y_xxy_xx_1,  \
                             ta1_y_xxy_xy_0,  \
                             ta1_y_xxy_xy_1,  \
                             ta1_y_xxy_yy_0,  \
                             ta1_y_xxy_yy_1,  \
                             ta1_y_xxyz_xx_0, \
                             ta1_y_xxyz_xy_0, \
                             ta1_y_xxyz_xz_0, \
                             ta1_y_xxyz_yy_0, \
                             ta1_y_xxyz_yz_0, \
                             ta1_y_xxyz_zz_0, \
                             ta1_y_xxz_xz_0,  \
                             ta1_y_xxz_xz_1,  \
                             ta1_y_xxz_zz_0,  \
                             ta1_y_xxz_zz_1,  \
                             ta1_y_xyz_yz_0,  \
                             ta1_y_xyz_yz_1,  \
                             ta1_y_yz_yz_0,   \
                             ta1_y_yz_yz_1,   \
                             ta_xxz_xz_1,     \
                             ta_xxz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyz_xx_0[i] = ta1_y_xxy_xx_0[i] * pa_z[i] - ta1_y_xxy_xx_1[i] * pc_z[i];

        ta1_y_xxyz_xy_0[i] = ta1_y_xxy_xy_0[i] * pa_z[i] - ta1_y_xxy_xy_1[i] * pc_z[i];

        ta1_y_xxyz_xz_0[i] = ta_xxz_xz_1[i] + ta1_y_xxz_xz_0[i] * pa_y[i] - ta1_y_xxz_xz_1[i] * pc_y[i];

        ta1_y_xxyz_yy_0[i] = ta1_y_xxy_yy_0[i] * pa_z[i] - ta1_y_xxy_yy_1[i] * pc_z[i];

        ta1_y_xxyz_yz_0[i] = ta1_y_yz_yz_0[i] * fe_0 - ta1_y_yz_yz_1[i] * fe_0 + ta1_y_xyz_yz_0[i] * pa_x[i] -
                             ta1_y_xyz_yz_1[i] * pc_x[i];

        ta1_y_xxyz_zz_0[i] = ta_xxz_zz_1[i] + ta1_y_xxz_zz_0[i] * pa_y[i] - ta1_y_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 120-126 components of targeted buffer : GD

    auto ta1_y_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 120);

    auto ta1_y_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 121);

    auto ta1_y_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 122);

    auto ta1_y_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 123);

    auto ta1_y_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 124);

    auto ta1_y_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 125);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xy_0,   \
                             ta1_y_xx_xy_1,   \
                             ta1_y_xxz_xx_0,  \
                             ta1_y_xxz_xx_1,  \
                             ta1_y_xxz_xy_0,  \
                             ta1_y_xxz_xy_1,  \
                             ta1_y_xxzz_xx_0, \
                             ta1_y_xxzz_xy_0, \
                             ta1_y_xxzz_xz_0, \
                             ta1_y_xxzz_yy_0, \
                             ta1_y_xxzz_yz_0, \
                             ta1_y_xxzz_zz_0, \
                             ta1_y_xzz_xz_0,  \
                             ta1_y_xzz_xz_1,  \
                             ta1_y_xzz_yy_0,  \
                             ta1_y_xzz_yy_1,  \
                             ta1_y_xzz_yz_0,  \
                             ta1_y_xzz_yz_1,  \
                             ta1_y_xzz_z_0,   \
                             ta1_y_xzz_z_1,   \
                             ta1_y_xzz_zz_0,  \
                             ta1_y_xzz_zz_1,  \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_yy_0,   \
                             ta1_y_zz_yy_1,   \
                             ta1_y_zz_yz_0,   \
                             ta1_y_zz_yz_1,   \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_xx_0[i] = ta1_y_xx_xx_0[i] * fe_0 - ta1_y_xx_xx_1[i] * fe_0 + ta1_y_xxz_xx_0[i] * pa_z[i] -
                             ta1_y_xxz_xx_1[i] * pc_z[i];

        ta1_y_xxzz_xy_0[i] = ta1_y_xx_xy_0[i] * fe_0 - ta1_y_xx_xy_1[i] * fe_0 + ta1_y_xxz_xy_0[i] * pa_z[i] -
                             ta1_y_xxz_xy_1[i] * pc_z[i];

        ta1_y_xxzz_xz_0[i] = ta1_y_zz_xz_0[i] * fe_0 - ta1_y_zz_xz_1[i] * fe_0 + ta1_y_xzz_z_0[i] * fe_0 -
                             ta1_y_xzz_z_1[i] * fe_0 + ta1_y_xzz_xz_0[i] * pa_x[i] - ta1_y_xzz_xz_1[i] * pc_x[i];

        ta1_y_xxzz_yy_0[i] = ta1_y_zz_yy_0[i] * fe_0 - ta1_y_zz_yy_1[i] * fe_0 + ta1_y_xzz_yy_0[i] * pa_x[i] -
                             ta1_y_xzz_yy_1[i] * pc_x[i];

        ta1_y_xxzz_yz_0[i] = ta1_y_zz_yz_0[i] * fe_0 - ta1_y_zz_yz_1[i] * fe_0 + ta1_y_xzz_yz_0[i] * pa_x[i] -
                             ta1_y_xzz_yz_1[i] * pc_x[i];

        ta1_y_xxzz_zz_0[i] = ta1_y_zz_zz_0[i] * fe_0 - ta1_y_zz_zz_1[i] * fe_0 + ta1_y_xzz_zz_0[i] * pa_x[i] -
                             ta1_y_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : GD

    auto ta1_y_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 126);

    auto ta1_y_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 127);

    auto ta1_y_xyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 128);

    auto ta1_y_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 129);

    auto ta1_y_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 130);

    auto ta1_y_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 131);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xyyy_xx_0, \
                             ta1_y_xyyy_xy_0, \
                             ta1_y_xyyy_xz_0, \
                             ta1_y_xyyy_yy_0, \
                             ta1_y_xyyy_yz_0, \
                             ta1_y_xyyy_zz_0, \
                             ta1_y_yyy_x_0,   \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_xx_0,  \
                             ta1_y_yyy_xx_1,  \
                             ta1_y_yyy_xy_0,  \
                             ta1_y_yyy_xy_1,  \
                             ta1_y_yyy_xz_0,  \
                             ta1_y_yyy_xz_1,  \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_yy_0,  \
                             ta1_y_yyy_yy_1,  \
                             ta1_y_yyy_yz_0,  \
                             ta1_y_yyy_yz_1,  \
                             ta1_y_yyy_z_0,   \
                             ta1_y_yyy_z_1,   \
                             ta1_y_yyy_zz_0,  \
                             ta1_y_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_xx_0[i] = 2.0 * ta1_y_yyy_x_0[i] * fe_0 - 2.0 * ta1_y_yyy_x_1[i] * fe_0 +
                             ta1_y_yyy_xx_0[i] * pa_x[i] - ta1_y_yyy_xx_1[i] * pc_x[i];

        ta1_y_xyyy_xy_0[i] = ta1_y_yyy_y_0[i] * fe_0 - ta1_y_yyy_y_1[i] * fe_0 + ta1_y_yyy_xy_0[i] * pa_x[i] -
                             ta1_y_yyy_xy_1[i] * pc_x[i];

        ta1_y_xyyy_xz_0[i] = ta1_y_yyy_z_0[i] * fe_0 - ta1_y_yyy_z_1[i] * fe_0 + ta1_y_yyy_xz_0[i] * pa_x[i] -
                             ta1_y_yyy_xz_1[i] * pc_x[i];

        ta1_y_xyyy_yy_0[i] = ta1_y_yyy_yy_0[i] * pa_x[i] - ta1_y_yyy_yy_1[i] * pc_x[i];

        ta1_y_xyyy_yz_0[i] = ta1_y_yyy_yz_0[i] * pa_x[i] - ta1_y_yyy_yz_1[i] * pc_x[i];

        ta1_y_xyyy_zz_0[i] = ta1_y_yyy_zz_0[i] * pa_x[i] - ta1_y_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 132-138 components of targeted buffer : GD

    auto ta1_y_xyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 132);

    auto ta1_y_xyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 133);

    auto ta1_y_xyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 134);

    auto ta1_y_xyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 135);

    auto ta1_y_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 136);

    auto ta1_y_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 137);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xyy_xx_0,  \
                             ta1_y_xyy_xx_1,  \
                             ta1_y_xyy_xy_0,  \
                             ta1_y_xyy_xy_1,  \
                             ta1_y_xyyz_xx_0, \
                             ta1_y_xyyz_xy_0, \
                             ta1_y_xyyz_xz_0, \
                             ta1_y_xyyz_yy_0, \
                             ta1_y_xyyz_yz_0, \
                             ta1_y_xyyz_zz_0, \
                             ta1_y_yyz_xz_0,  \
                             ta1_y_yyz_xz_1,  \
                             ta1_y_yyz_yy_0,  \
                             ta1_y_yyz_yy_1,  \
                             ta1_y_yyz_yz_0,  \
                             ta1_y_yyz_yz_1,  \
                             ta1_y_yyz_z_0,   \
                             ta1_y_yyz_z_1,   \
                             ta1_y_yyz_zz_0,  \
                             ta1_y_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyz_xx_0[i] = ta1_y_xyy_xx_0[i] * pa_z[i] - ta1_y_xyy_xx_1[i] * pc_z[i];

        ta1_y_xyyz_xy_0[i] = ta1_y_xyy_xy_0[i] * pa_z[i] - ta1_y_xyy_xy_1[i] * pc_z[i];

        ta1_y_xyyz_xz_0[i] = ta1_y_yyz_z_0[i] * fe_0 - ta1_y_yyz_z_1[i] * fe_0 + ta1_y_yyz_xz_0[i] * pa_x[i] -
                             ta1_y_yyz_xz_1[i] * pc_x[i];

        ta1_y_xyyz_yy_0[i] = ta1_y_yyz_yy_0[i] * pa_x[i] - ta1_y_yyz_yy_1[i] * pc_x[i];

        ta1_y_xyyz_yz_0[i] = ta1_y_yyz_yz_0[i] * pa_x[i] - ta1_y_yyz_yz_1[i] * pc_x[i];

        ta1_y_xyyz_zz_0[i] = ta1_y_yyz_zz_0[i] * pa_x[i] - ta1_y_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 138-144 components of targeted buffer : GD

    auto ta1_y_xyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 138);

    auto ta1_y_xyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 139);

    auto ta1_y_xyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 140);

    auto ta1_y_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 141);

    auto ta1_y_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 142);

    auto ta1_y_xyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 143);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xyzz_xx_0, \
                             ta1_y_xyzz_xy_0, \
                             ta1_y_xyzz_xz_0, \
                             ta1_y_xyzz_yy_0, \
                             ta1_y_xyzz_yz_0, \
                             ta1_y_xyzz_zz_0, \
                             ta1_y_xzz_xx_0,  \
                             ta1_y_xzz_xx_1,  \
                             ta1_y_xzz_xz_0,  \
                             ta1_y_xzz_xz_1,  \
                             ta1_y_yzz_xy_0,  \
                             ta1_y_yzz_xy_1,  \
                             ta1_y_yzz_y_0,   \
                             ta1_y_yzz_y_1,   \
                             ta1_y_yzz_yy_0,  \
                             ta1_y_yzz_yy_1,  \
                             ta1_y_yzz_yz_0,  \
                             ta1_y_yzz_yz_1,  \
                             ta1_y_yzz_zz_0,  \
                             ta1_y_yzz_zz_1,  \
                             ta_xzz_xx_1,     \
                             ta_xzz_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzz_xx_0[i] = ta_xzz_xx_1[i] + ta1_y_xzz_xx_0[i] * pa_y[i] - ta1_y_xzz_xx_1[i] * pc_y[i];

        ta1_y_xyzz_xy_0[i] = ta1_y_yzz_y_0[i] * fe_0 - ta1_y_yzz_y_1[i] * fe_0 + ta1_y_yzz_xy_0[i] * pa_x[i] -
                             ta1_y_yzz_xy_1[i] * pc_x[i];

        ta1_y_xyzz_xz_0[i] = ta_xzz_xz_1[i] + ta1_y_xzz_xz_0[i] * pa_y[i] - ta1_y_xzz_xz_1[i] * pc_y[i];

        ta1_y_xyzz_yy_0[i] = ta1_y_yzz_yy_0[i] * pa_x[i] - ta1_y_yzz_yy_1[i] * pc_x[i];

        ta1_y_xyzz_yz_0[i] = ta1_y_yzz_yz_0[i] * pa_x[i] - ta1_y_yzz_yz_1[i] * pc_x[i];

        ta1_y_xyzz_zz_0[i] = ta1_y_yzz_zz_0[i] * pa_x[i] - ta1_y_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 144-150 components of targeted buffer : GD

    auto ta1_y_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 144);

    auto ta1_y_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 145);

    auto ta1_y_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 146);

    auto ta1_y_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 147);

    auto ta1_y_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 148);

    auto ta1_y_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 149);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xzzz_xx_0, \
                             ta1_y_xzzz_xy_0, \
                             ta1_y_xzzz_xz_0, \
                             ta1_y_xzzz_yy_0, \
                             ta1_y_xzzz_yz_0, \
                             ta1_y_xzzz_zz_0, \
                             ta1_y_zzz_x_0,   \
                             ta1_y_zzz_x_1,   \
                             ta1_y_zzz_xx_0,  \
                             ta1_y_zzz_xx_1,  \
                             ta1_y_zzz_xy_0,  \
                             ta1_y_zzz_xy_1,  \
                             ta1_y_zzz_xz_0,  \
                             ta1_y_zzz_xz_1,  \
                             ta1_y_zzz_y_0,   \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_yy_0,  \
                             ta1_y_zzz_yy_1,  \
                             ta1_y_zzz_yz_0,  \
                             ta1_y_zzz_yz_1,  \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             ta1_y_zzz_zz_0,  \
                             ta1_y_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_xx_0[i] = 2.0 * ta1_y_zzz_x_0[i] * fe_0 - 2.0 * ta1_y_zzz_x_1[i] * fe_0 +
                             ta1_y_zzz_xx_0[i] * pa_x[i] - ta1_y_zzz_xx_1[i] * pc_x[i];

        ta1_y_xzzz_xy_0[i] = ta1_y_zzz_y_0[i] * fe_0 - ta1_y_zzz_y_1[i] * fe_0 + ta1_y_zzz_xy_0[i] * pa_x[i] -
                             ta1_y_zzz_xy_1[i] * pc_x[i];

        ta1_y_xzzz_xz_0[i] = ta1_y_zzz_z_0[i] * fe_0 - ta1_y_zzz_z_1[i] * fe_0 + ta1_y_zzz_xz_0[i] * pa_x[i] -
                             ta1_y_zzz_xz_1[i] * pc_x[i];

        ta1_y_xzzz_yy_0[i] = ta1_y_zzz_yy_0[i] * pa_x[i] - ta1_y_zzz_yy_1[i] * pc_x[i];

        ta1_y_xzzz_yz_0[i] = ta1_y_zzz_yz_0[i] * pa_x[i] - ta1_y_zzz_yz_1[i] * pc_x[i];

        ta1_y_xzzz_zz_0[i] = ta1_y_zzz_zz_0[i] * pa_x[i] - ta1_y_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 150-156 components of targeted buffer : GD

    auto ta1_y_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 150);

    auto ta1_y_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 151);

    auto ta1_y_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 152);

    auto ta1_y_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 153);

    auto ta1_y_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 154);

    auto ta1_y_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 155);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_yy_xx_0,   \
                             ta1_y_yy_xx_1,   \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_xz_0,   \
                             ta1_y_yy_xz_1,   \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yy_zz_0,   \
                             ta1_y_yy_zz_1,   \
                             ta1_y_yyy_x_0,   \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_xx_0,  \
                             ta1_y_yyy_xx_1,  \
                             ta1_y_yyy_xy_0,  \
                             ta1_y_yyy_xy_1,  \
                             ta1_y_yyy_xz_0,  \
                             ta1_y_yyy_xz_1,  \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_yy_0,  \
                             ta1_y_yyy_yy_1,  \
                             ta1_y_yyy_yz_0,  \
                             ta1_y_yyy_yz_1,  \
                             ta1_y_yyy_z_0,   \
                             ta1_y_yyy_z_1,   \
                             ta1_y_yyy_zz_0,  \
                             ta1_y_yyy_zz_1,  \
                             ta1_y_yyyy_xx_0, \
                             ta1_y_yyyy_xy_0, \
                             ta1_y_yyyy_xz_0, \
                             ta1_y_yyyy_yy_0, \
                             ta1_y_yyyy_yz_0, \
                             ta1_y_yyyy_zz_0, \
                             ta_yyy_xx_1,     \
                             ta_yyy_xy_1,     \
                             ta_yyy_xz_1,     \
                             ta_yyy_yy_1,     \
                             ta_yyy_yz_1,     \
                             ta_yyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_xx_0[i] = 3.0 * ta1_y_yy_xx_0[i] * fe_0 - 3.0 * ta1_y_yy_xx_1[i] * fe_0 + ta_yyy_xx_1[i] +
                             ta1_y_yyy_xx_0[i] * pa_y[i] - ta1_y_yyy_xx_1[i] * pc_y[i];

        ta1_y_yyyy_xy_0[i] = 3.0 * ta1_y_yy_xy_0[i] * fe_0 - 3.0 * ta1_y_yy_xy_1[i] * fe_0 + ta1_y_yyy_x_0[i] * fe_0 -
                             ta1_y_yyy_x_1[i] * fe_0 + ta_yyy_xy_1[i] + ta1_y_yyy_xy_0[i] * pa_y[i] -
                             ta1_y_yyy_xy_1[i] * pc_y[i];

        ta1_y_yyyy_xz_0[i] = 3.0 * ta1_y_yy_xz_0[i] * fe_0 - 3.0 * ta1_y_yy_xz_1[i] * fe_0 + ta_yyy_xz_1[i] +
                             ta1_y_yyy_xz_0[i] * pa_y[i] - ta1_y_yyy_xz_1[i] * pc_y[i];

        ta1_y_yyyy_yy_0[i] = 3.0 * ta1_y_yy_yy_0[i] * fe_0 - 3.0 * ta1_y_yy_yy_1[i] * fe_0 +
                             2.0 * ta1_y_yyy_y_0[i] * fe_0 - 2.0 * ta1_y_yyy_y_1[i] * fe_0 + ta_yyy_yy_1[i] +
                             ta1_y_yyy_yy_0[i] * pa_y[i] - ta1_y_yyy_yy_1[i] * pc_y[i];

        ta1_y_yyyy_yz_0[i] = 3.0 * ta1_y_yy_yz_0[i] * fe_0 - 3.0 * ta1_y_yy_yz_1[i] * fe_0 + ta1_y_yyy_z_0[i] * fe_0 -
                             ta1_y_yyy_z_1[i] * fe_0 + ta_yyy_yz_1[i] + ta1_y_yyy_yz_0[i] * pa_y[i] -
                             ta1_y_yyy_yz_1[i] * pc_y[i];

        ta1_y_yyyy_zz_0[i] = 3.0 * ta1_y_yy_zz_0[i] * fe_0 - 3.0 * ta1_y_yy_zz_1[i] * fe_0 + ta_yyy_zz_1[i] +
                             ta1_y_yyy_zz_0[i] * pa_y[i] - ta1_y_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 156-162 components of targeted buffer : GD

    auto ta1_y_yyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 156);

    auto ta1_y_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 157);

    auto ta1_y_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 158);

    auto ta1_y_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 159);

    auto ta1_y_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 160);

    auto ta1_y_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 161);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_yyy_x_0,   \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_xx_0,  \
                             ta1_y_yyy_xx_1,  \
                             ta1_y_yyy_xy_0,  \
                             ta1_y_yyy_xy_1,  \
                             ta1_y_yyy_xz_0,  \
                             ta1_y_yyy_xz_1,  \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_yy_0,  \
                             ta1_y_yyy_yy_1,  \
                             ta1_y_yyy_yz_0,  \
                             ta1_y_yyy_yz_1,  \
                             ta1_y_yyy_z_0,   \
                             ta1_y_yyy_z_1,   \
                             ta1_y_yyy_zz_0,  \
                             ta1_y_yyy_zz_1,  \
                             ta1_y_yyyz_xx_0, \
                             ta1_y_yyyz_xy_0, \
                             ta1_y_yyyz_xz_0, \
                             ta1_y_yyyz_yy_0, \
                             ta1_y_yyyz_yz_0, \
                             ta1_y_yyyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_xx_0[i] = ta1_y_yyy_xx_0[i] * pa_z[i] - ta1_y_yyy_xx_1[i] * pc_z[i];

        ta1_y_yyyz_xy_0[i] = ta1_y_yyy_xy_0[i] * pa_z[i] - ta1_y_yyy_xy_1[i] * pc_z[i];

        ta1_y_yyyz_xz_0[i] = ta1_y_yyy_x_0[i] * fe_0 - ta1_y_yyy_x_1[i] * fe_0 + ta1_y_yyy_xz_0[i] * pa_z[i] -
                             ta1_y_yyy_xz_1[i] * pc_z[i];

        ta1_y_yyyz_yy_0[i] = ta1_y_yyy_yy_0[i] * pa_z[i] - ta1_y_yyy_yy_1[i] * pc_z[i];

        ta1_y_yyyz_yz_0[i] = ta1_y_yyy_y_0[i] * fe_0 - ta1_y_yyy_y_1[i] * fe_0 + ta1_y_yyy_yz_0[i] * pa_z[i] -
                             ta1_y_yyy_yz_1[i] * pc_z[i];

        ta1_y_yyyz_zz_0[i] = 2.0 * ta1_y_yyy_z_0[i] * fe_0 - 2.0 * ta1_y_yyy_z_1[i] * fe_0 +
                             ta1_y_yyy_zz_0[i] * pa_z[i] - ta1_y_yyy_zz_1[i] * pc_z[i];
    }

    // Set up 162-168 components of targeted buffer : GD

    auto ta1_y_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 162);

    auto ta1_y_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 163);

    auto ta1_y_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 164);

    auto ta1_y_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 165);

    auto ta1_y_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 166);

    auto ta1_y_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 167);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yy_xx_0,   \
                             ta1_y_yy_xx_1,   \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yyz_xx_0,  \
                             ta1_y_yyz_xx_1,  \
                             ta1_y_yyz_xy_0,  \
                             ta1_y_yyz_xy_1,  \
                             ta1_y_yyz_y_0,   \
                             ta1_y_yyz_y_1,   \
                             ta1_y_yyz_yy_0,  \
                             ta1_y_yyz_yy_1,  \
                             ta1_y_yyz_yz_0,  \
                             ta1_y_yyz_yz_1,  \
                             ta1_y_yyzz_xx_0, \
                             ta1_y_yyzz_xy_0, \
                             ta1_y_yyzz_xz_0, \
                             ta1_y_yyzz_yy_0, \
                             ta1_y_yyzz_yz_0, \
                             ta1_y_yyzz_zz_0, \
                             ta1_y_yzz_xz_0,  \
                             ta1_y_yzz_xz_1,  \
                             ta1_y_yzz_zz_0,  \
                             ta1_y_yzz_zz_1,  \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             ta_yzz_xz_1,     \
                             ta_yzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_xx_0[i] = ta1_y_yy_xx_0[i] * fe_0 - ta1_y_yy_xx_1[i] * fe_0 + ta1_y_yyz_xx_0[i] * pa_z[i] -
                             ta1_y_yyz_xx_1[i] * pc_z[i];

        ta1_y_yyzz_xy_0[i] = ta1_y_yy_xy_0[i] * fe_0 - ta1_y_yy_xy_1[i] * fe_0 + ta1_y_yyz_xy_0[i] * pa_z[i] -
                             ta1_y_yyz_xy_1[i] * pc_z[i];

        ta1_y_yyzz_xz_0[i] = ta1_y_zz_xz_0[i] * fe_0 - ta1_y_zz_xz_1[i] * fe_0 + ta_yzz_xz_1[i] +
                             ta1_y_yzz_xz_0[i] * pa_y[i] - ta1_y_yzz_xz_1[i] * pc_y[i];

        ta1_y_yyzz_yy_0[i] = ta1_y_yy_yy_0[i] * fe_0 - ta1_y_yy_yy_1[i] * fe_0 + ta1_y_yyz_yy_0[i] * pa_z[i] -
                             ta1_y_yyz_yy_1[i] * pc_z[i];

        ta1_y_yyzz_yz_0[i] = ta1_y_yy_yz_0[i] * fe_0 - ta1_y_yy_yz_1[i] * fe_0 + ta1_y_yyz_y_0[i] * fe_0 -
                             ta1_y_yyz_y_1[i] * fe_0 + ta1_y_yyz_yz_0[i] * pa_z[i] - ta1_y_yyz_yz_1[i] * pc_z[i];

        ta1_y_yyzz_zz_0[i] = ta1_y_zz_zz_0[i] * fe_0 - ta1_y_zz_zz_1[i] * fe_0 + ta_yzz_zz_1[i] +
                             ta1_y_yzz_zz_0[i] * pa_y[i] - ta1_y_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 168-174 components of targeted buffer : GD

    auto ta1_y_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 168);

    auto ta1_y_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 169);

    auto ta1_y_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 170);

    auto ta1_y_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 171);

    auto ta1_y_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 172);

    auto ta1_y_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 173);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yz_xy_0,   \
                             ta1_y_yz_xy_1,   \
                             ta1_y_yz_yy_0,   \
                             ta1_y_yz_yy_1,   \
                             ta1_y_yzz_xy_0,  \
                             ta1_y_yzz_xy_1,  \
                             ta1_y_yzz_yy_0,  \
                             ta1_y_yzz_yy_1,  \
                             ta1_y_yzzz_xx_0, \
                             ta1_y_yzzz_xy_0, \
                             ta1_y_yzzz_xz_0, \
                             ta1_y_yzzz_yy_0, \
                             ta1_y_yzzz_yz_0, \
                             ta1_y_yzzz_zz_0, \
                             ta1_y_zzz_xx_0,  \
                             ta1_y_zzz_xx_1,  \
                             ta1_y_zzz_xz_0,  \
                             ta1_y_zzz_xz_1,  \
                             ta1_y_zzz_yz_0,  \
                             ta1_y_zzz_yz_1,  \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             ta1_y_zzz_zz_0,  \
                             ta1_y_zzz_zz_1,  \
                             ta_zzz_xx_1,     \
                             ta_zzz_xz_1,     \
                             ta_zzz_yz_1,     \
                             ta_zzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_xx_0[i] = ta_zzz_xx_1[i] + ta1_y_zzz_xx_0[i] * pa_y[i] - ta1_y_zzz_xx_1[i] * pc_y[i];

        ta1_y_yzzz_xy_0[i] = 2.0 * ta1_y_yz_xy_0[i] * fe_0 - 2.0 * ta1_y_yz_xy_1[i] * fe_0 +
                             ta1_y_yzz_xy_0[i] * pa_z[i] - ta1_y_yzz_xy_1[i] * pc_z[i];

        ta1_y_yzzz_xz_0[i] = ta_zzz_xz_1[i] + ta1_y_zzz_xz_0[i] * pa_y[i] - ta1_y_zzz_xz_1[i] * pc_y[i];

        ta1_y_yzzz_yy_0[i] = 2.0 * ta1_y_yz_yy_0[i] * fe_0 - 2.0 * ta1_y_yz_yy_1[i] * fe_0 +
                             ta1_y_yzz_yy_0[i] * pa_z[i] - ta1_y_yzz_yy_1[i] * pc_z[i];

        ta1_y_yzzz_yz_0[i] = ta1_y_zzz_z_0[i] * fe_0 - ta1_y_zzz_z_1[i] * fe_0 + ta_zzz_yz_1[i] +
                             ta1_y_zzz_yz_0[i] * pa_y[i] - ta1_y_zzz_yz_1[i] * pc_y[i];

        ta1_y_yzzz_zz_0[i] = ta_zzz_zz_1[i] + ta1_y_zzz_zz_0[i] * pa_y[i] - ta1_y_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 174-180 components of targeted buffer : GD

    auto ta1_y_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 174);

    auto ta1_y_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 175);

    auto ta1_y_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 176);

    auto ta1_y_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 177);

    auto ta1_y_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 178);

    auto ta1_y_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 179);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_zz_xx_0,   \
                             ta1_y_zz_xx_1,   \
                             ta1_y_zz_xy_0,   \
                             ta1_y_zz_xy_1,   \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_yy_0,   \
                             ta1_y_zz_yy_1,   \
                             ta1_y_zz_yz_0,   \
                             ta1_y_zz_yz_1,   \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             ta1_y_zzz_x_0,   \
                             ta1_y_zzz_x_1,   \
                             ta1_y_zzz_xx_0,  \
                             ta1_y_zzz_xx_1,  \
                             ta1_y_zzz_xy_0,  \
                             ta1_y_zzz_xy_1,  \
                             ta1_y_zzz_xz_0,  \
                             ta1_y_zzz_xz_1,  \
                             ta1_y_zzz_y_0,   \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_yy_0,  \
                             ta1_y_zzz_yy_1,  \
                             ta1_y_zzz_yz_0,  \
                             ta1_y_zzz_yz_1,  \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             ta1_y_zzz_zz_0,  \
                             ta1_y_zzz_zz_1,  \
                             ta1_y_zzzz_xx_0, \
                             ta1_y_zzzz_xy_0, \
                             ta1_y_zzzz_xz_0, \
                             ta1_y_zzzz_yy_0, \
                             ta1_y_zzzz_yz_0, \
                             ta1_y_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_xx_0[i] = 3.0 * ta1_y_zz_xx_0[i] * fe_0 - 3.0 * ta1_y_zz_xx_1[i] * fe_0 +
                             ta1_y_zzz_xx_0[i] * pa_z[i] - ta1_y_zzz_xx_1[i] * pc_z[i];

        ta1_y_zzzz_xy_0[i] = 3.0 * ta1_y_zz_xy_0[i] * fe_0 - 3.0 * ta1_y_zz_xy_1[i] * fe_0 +
                             ta1_y_zzz_xy_0[i] * pa_z[i] - ta1_y_zzz_xy_1[i] * pc_z[i];

        ta1_y_zzzz_xz_0[i] = 3.0 * ta1_y_zz_xz_0[i] * fe_0 - 3.0 * ta1_y_zz_xz_1[i] * fe_0 + ta1_y_zzz_x_0[i] * fe_0 -
                             ta1_y_zzz_x_1[i] * fe_0 + ta1_y_zzz_xz_0[i] * pa_z[i] - ta1_y_zzz_xz_1[i] * pc_z[i];

        ta1_y_zzzz_yy_0[i] = 3.0 * ta1_y_zz_yy_0[i] * fe_0 - 3.0 * ta1_y_zz_yy_1[i] * fe_0 +
                             ta1_y_zzz_yy_0[i] * pa_z[i] - ta1_y_zzz_yy_1[i] * pc_z[i];

        ta1_y_zzzz_yz_0[i] = 3.0 * ta1_y_zz_yz_0[i] * fe_0 - 3.0 * ta1_y_zz_yz_1[i] * fe_0 + ta1_y_zzz_y_0[i] * fe_0 -
                             ta1_y_zzz_y_1[i] * fe_0 + ta1_y_zzz_yz_0[i] * pa_z[i] - ta1_y_zzz_yz_1[i] * pc_z[i];

        ta1_y_zzzz_zz_0[i] = 3.0 * ta1_y_zz_zz_0[i] * fe_0 - 3.0 * ta1_y_zz_zz_1[i] * fe_0 +
                             2.0 * ta1_y_zzz_z_0[i] * fe_0 - 2.0 * ta1_y_zzz_z_1[i] * fe_0 +
                             ta1_y_zzz_zz_0[i] * pa_z[i] - ta1_y_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 180-186 components of targeted buffer : GD

    auto ta1_z_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 180);

    auto ta1_z_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 181);

    auto ta1_z_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 182);

    auto ta1_z_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 183);

    auto ta1_z_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 184);

    auto ta1_z_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 185);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xy_0,   \
                             ta1_z_xx_xy_1,   \
                             ta1_z_xx_xz_0,   \
                             ta1_z_xx_xz_1,   \
                             ta1_z_xx_yy_0,   \
                             ta1_z_xx_yy_1,   \
                             ta1_z_xx_yz_0,   \
                             ta1_z_xx_yz_1,   \
                             ta1_z_xx_zz_0,   \
                             ta1_z_xx_zz_1,   \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_xx_0,  \
                             ta1_z_xxx_xx_1,  \
                             ta1_z_xxx_xy_0,  \
                             ta1_z_xxx_xy_1,  \
                             ta1_z_xxx_xz_0,  \
                             ta1_z_xxx_xz_1,  \
                             ta1_z_xxx_y_0,   \
                             ta1_z_xxx_y_1,   \
                             ta1_z_xxx_yy_0,  \
                             ta1_z_xxx_yy_1,  \
                             ta1_z_xxx_yz_0,  \
                             ta1_z_xxx_yz_1,  \
                             ta1_z_xxx_z_0,   \
                             ta1_z_xxx_z_1,   \
                             ta1_z_xxx_zz_0,  \
                             ta1_z_xxx_zz_1,  \
                             ta1_z_xxxx_xx_0, \
                             ta1_z_xxxx_xy_0, \
                             ta1_z_xxxx_xz_0, \
                             ta1_z_xxxx_yy_0, \
                             ta1_z_xxxx_yz_0, \
                             ta1_z_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_xx_0[i] = 3.0 * ta1_z_xx_xx_0[i] * fe_0 - 3.0 * ta1_z_xx_xx_1[i] * fe_0 +
                             2.0 * ta1_z_xxx_x_0[i] * fe_0 - 2.0 * ta1_z_xxx_x_1[i] * fe_0 +
                             ta1_z_xxx_xx_0[i] * pa_x[i] - ta1_z_xxx_xx_1[i] * pc_x[i];

        ta1_z_xxxx_xy_0[i] = 3.0 * ta1_z_xx_xy_0[i] * fe_0 - 3.0 * ta1_z_xx_xy_1[i] * fe_0 + ta1_z_xxx_y_0[i] * fe_0 -
                             ta1_z_xxx_y_1[i] * fe_0 + ta1_z_xxx_xy_0[i] * pa_x[i] - ta1_z_xxx_xy_1[i] * pc_x[i];

        ta1_z_xxxx_xz_0[i] = 3.0 * ta1_z_xx_xz_0[i] * fe_0 - 3.0 * ta1_z_xx_xz_1[i] * fe_0 + ta1_z_xxx_z_0[i] * fe_0 -
                             ta1_z_xxx_z_1[i] * fe_0 + ta1_z_xxx_xz_0[i] * pa_x[i] - ta1_z_xxx_xz_1[i] * pc_x[i];

        ta1_z_xxxx_yy_0[i] = 3.0 * ta1_z_xx_yy_0[i] * fe_0 - 3.0 * ta1_z_xx_yy_1[i] * fe_0 +
                             ta1_z_xxx_yy_0[i] * pa_x[i] - ta1_z_xxx_yy_1[i] * pc_x[i];

        ta1_z_xxxx_yz_0[i] = 3.0 * ta1_z_xx_yz_0[i] * fe_0 - 3.0 * ta1_z_xx_yz_1[i] * fe_0 +
                             ta1_z_xxx_yz_0[i] * pa_x[i] - ta1_z_xxx_yz_1[i] * pc_x[i];

        ta1_z_xxxx_zz_0[i] = 3.0 * ta1_z_xx_zz_0[i] * fe_0 - 3.0 * ta1_z_xx_zz_1[i] * fe_0 +
                             ta1_z_xxx_zz_0[i] * pa_x[i] - ta1_z_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : GD

    auto ta1_z_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 186);

    auto ta1_z_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 187);

    auto ta1_z_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 188);

    auto ta1_z_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 189);

    auto ta1_z_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 190);

    auto ta1_z_xxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 191);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_xx_0,  \
                             ta1_z_xxx_xx_1,  \
                             ta1_z_xxx_xy_0,  \
                             ta1_z_xxx_xy_1,  \
                             ta1_z_xxx_xz_0,  \
                             ta1_z_xxx_xz_1,  \
                             ta1_z_xxx_zz_0,  \
                             ta1_z_xxx_zz_1,  \
                             ta1_z_xxxy_xx_0, \
                             ta1_z_xxxy_xy_0, \
                             ta1_z_xxxy_xz_0, \
                             ta1_z_xxxy_yy_0, \
                             ta1_z_xxxy_yz_0, \
                             ta1_z_xxxy_zz_0, \
                             ta1_z_xxy_yy_0,  \
                             ta1_z_xxy_yy_1,  \
                             ta1_z_xxy_yz_0,  \
                             ta1_z_xxy_yz_1,  \
                             ta1_z_xy_yy_0,   \
                             ta1_z_xy_yy_1,   \
                             ta1_z_xy_yz_0,   \
                             ta1_z_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_xx_0[i] = ta1_z_xxx_xx_0[i] * pa_y[i] - ta1_z_xxx_xx_1[i] * pc_y[i];

        ta1_z_xxxy_xy_0[i] = ta1_z_xxx_x_0[i] * fe_0 - ta1_z_xxx_x_1[i] * fe_0 + ta1_z_xxx_xy_0[i] * pa_y[i] -
                             ta1_z_xxx_xy_1[i] * pc_y[i];

        ta1_z_xxxy_xz_0[i] = ta1_z_xxx_xz_0[i] * pa_y[i] - ta1_z_xxx_xz_1[i] * pc_y[i];

        ta1_z_xxxy_yy_0[i] = 2.0 * ta1_z_xy_yy_0[i] * fe_0 - 2.0 * ta1_z_xy_yy_1[i] * fe_0 +
                             ta1_z_xxy_yy_0[i] * pa_x[i] - ta1_z_xxy_yy_1[i] * pc_x[i];

        ta1_z_xxxy_yz_0[i] = 2.0 * ta1_z_xy_yz_0[i] * fe_0 - 2.0 * ta1_z_xy_yz_1[i] * fe_0 +
                             ta1_z_xxy_yz_0[i] * pa_x[i] - ta1_z_xxy_yz_1[i] * pc_x[i];

        ta1_z_xxxy_zz_0[i] = ta1_z_xxx_zz_0[i] * pa_y[i] - ta1_z_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 192-198 components of targeted buffer : GD

    auto ta1_z_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 192);

    auto ta1_z_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 193);

    auto ta1_z_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 194);

    auto ta1_z_xxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 195);

    auto ta1_z_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 196);

    auto ta1_z_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 197);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_xx_0,  \
                             ta1_z_xxx_xx_1,  \
                             ta1_z_xxx_xy_0,  \
                             ta1_z_xxx_xy_1,  \
                             ta1_z_xxx_xz_0,  \
                             ta1_z_xxx_xz_1,  \
                             ta1_z_xxx_yy_0,  \
                             ta1_z_xxx_yy_1,  \
                             ta1_z_xxxz_xx_0, \
                             ta1_z_xxxz_xy_0, \
                             ta1_z_xxxz_xz_0, \
                             ta1_z_xxxz_yy_0, \
                             ta1_z_xxxz_yz_0, \
                             ta1_z_xxxz_zz_0, \
                             ta1_z_xxz_yz_0,  \
                             ta1_z_xxz_yz_1,  \
                             ta1_z_xxz_zz_0,  \
                             ta1_z_xxz_zz_1,  \
                             ta1_z_xz_yz_0,   \
                             ta1_z_xz_yz_1,   \
                             ta1_z_xz_zz_0,   \
                             ta1_z_xz_zz_1,   \
                             ta_xxx_xx_1,     \
                             ta_xxx_xy_1,     \
                             ta_xxx_xz_1,     \
                             ta_xxx_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_xx_0[i] = ta_xxx_xx_1[i] + ta1_z_xxx_xx_0[i] * pa_z[i] - ta1_z_xxx_xx_1[i] * pc_z[i];

        ta1_z_xxxz_xy_0[i] = ta_xxx_xy_1[i] + ta1_z_xxx_xy_0[i] * pa_z[i] - ta1_z_xxx_xy_1[i] * pc_z[i];

        ta1_z_xxxz_xz_0[i] = ta1_z_xxx_x_0[i] * fe_0 - ta1_z_xxx_x_1[i] * fe_0 + ta_xxx_xz_1[i] +
                             ta1_z_xxx_xz_0[i] * pa_z[i] - ta1_z_xxx_xz_1[i] * pc_z[i];

        ta1_z_xxxz_yy_0[i] = ta_xxx_yy_1[i] + ta1_z_xxx_yy_0[i] * pa_z[i] - ta1_z_xxx_yy_1[i] * pc_z[i];

        ta1_z_xxxz_yz_0[i] = 2.0 * ta1_z_xz_yz_0[i] * fe_0 - 2.0 * ta1_z_xz_yz_1[i] * fe_0 +
                             ta1_z_xxz_yz_0[i] * pa_x[i] - ta1_z_xxz_yz_1[i] * pc_x[i];

        ta1_z_xxxz_zz_0[i] = 2.0 * ta1_z_xz_zz_0[i] * fe_0 - 2.0 * ta1_z_xz_zz_1[i] * fe_0 +
                             ta1_z_xxz_zz_0[i] * pa_x[i] - ta1_z_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 198-204 components of targeted buffer : GD

    auto ta1_z_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 198);

    auto ta1_z_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 199);

    auto ta1_z_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 200);

    auto ta1_z_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 201);

    auto ta1_z_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 202);

    auto ta1_z_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 203);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xz_0,   \
                             ta1_z_xx_xz_1,   \
                             ta1_z_xxy_xx_0,  \
                             ta1_z_xxy_xx_1,  \
                             ta1_z_xxy_xz_0,  \
                             ta1_z_xxy_xz_1,  \
                             ta1_z_xxyy_xx_0, \
                             ta1_z_xxyy_xy_0, \
                             ta1_z_xxyy_xz_0, \
                             ta1_z_xxyy_yy_0, \
                             ta1_z_xxyy_yz_0, \
                             ta1_z_xxyy_zz_0, \
                             ta1_z_xyy_xy_0,  \
                             ta1_z_xyy_xy_1,  \
                             ta1_z_xyy_y_0,   \
                             ta1_z_xyy_y_1,   \
                             ta1_z_xyy_yy_0,  \
                             ta1_z_xyy_yy_1,  \
                             ta1_z_xyy_yz_0,  \
                             ta1_z_xyy_yz_1,  \
                             ta1_z_xyy_zz_0,  \
                             ta1_z_xyy_zz_1,  \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yy_yz_0,   \
                             ta1_z_yy_yz_1,   \
                             ta1_z_yy_zz_0,   \
                             ta1_z_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_xx_0[i] = ta1_z_xx_xx_0[i] * fe_0 - ta1_z_xx_xx_1[i] * fe_0 + ta1_z_xxy_xx_0[i] * pa_y[i] -
                             ta1_z_xxy_xx_1[i] * pc_y[i];

        ta1_z_xxyy_xy_0[i] = ta1_z_yy_xy_0[i] * fe_0 - ta1_z_yy_xy_1[i] * fe_0 + ta1_z_xyy_y_0[i] * fe_0 -
                             ta1_z_xyy_y_1[i] * fe_0 + ta1_z_xyy_xy_0[i] * pa_x[i] - ta1_z_xyy_xy_1[i] * pc_x[i];

        ta1_z_xxyy_xz_0[i] = ta1_z_xx_xz_0[i] * fe_0 - ta1_z_xx_xz_1[i] * fe_0 + ta1_z_xxy_xz_0[i] * pa_y[i] -
                             ta1_z_xxy_xz_1[i] * pc_y[i];

        ta1_z_xxyy_yy_0[i] = ta1_z_yy_yy_0[i] * fe_0 - ta1_z_yy_yy_1[i] * fe_0 + ta1_z_xyy_yy_0[i] * pa_x[i] -
                             ta1_z_xyy_yy_1[i] * pc_x[i];

        ta1_z_xxyy_yz_0[i] = ta1_z_yy_yz_0[i] * fe_0 - ta1_z_yy_yz_1[i] * fe_0 + ta1_z_xyy_yz_0[i] * pa_x[i] -
                             ta1_z_xyy_yz_1[i] * pc_x[i];

        ta1_z_xxyy_zz_0[i] = ta1_z_yy_zz_0[i] * fe_0 - ta1_z_yy_zz_1[i] * fe_0 + ta1_z_xyy_zz_0[i] * pa_x[i] -
                             ta1_z_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 204-210 components of targeted buffer : GD

    auto ta1_z_xxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 204);

    auto ta1_z_xxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 205);

    auto ta1_z_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 206);

    auto ta1_z_xxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 207);

    auto ta1_z_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 208);

    auto ta1_z_xxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 209);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_xxy_xy_0,  \
                             ta1_z_xxy_xy_1,  \
                             ta1_z_xxy_yy_0,  \
                             ta1_z_xxy_yy_1,  \
                             ta1_z_xxyz_xx_0, \
                             ta1_z_xxyz_xy_0, \
                             ta1_z_xxyz_xz_0, \
                             ta1_z_xxyz_yy_0, \
                             ta1_z_xxyz_yz_0, \
                             ta1_z_xxyz_zz_0, \
                             ta1_z_xxz_xx_0,  \
                             ta1_z_xxz_xx_1,  \
                             ta1_z_xxz_xz_0,  \
                             ta1_z_xxz_xz_1,  \
                             ta1_z_xxz_zz_0,  \
                             ta1_z_xxz_zz_1,  \
                             ta1_z_xyz_yz_0,  \
                             ta1_z_xyz_yz_1,  \
                             ta1_z_yz_yz_0,   \
                             ta1_z_yz_yz_1,   \
                             ta_xxy_xy_1,     \
                             ta_xxy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyz_xx_0[i] = ta1_z_xxz_xx_0[i] * pa_y[i] - ta1_z_xxz_xx_1[i] * pc_y[i];

        ta1_z_xxyz_xy_0[i] = ta_xxy_xy_1[i] + ta1_z_xxy_xy_0[i] * pa_z[i] - ta1_z_xxy_xy_1[i] * pc_z[i];

        ta1_z_xxyz_xz_0[i] = ta1_z_xxz_xz_0[i] * pa_y[i] - ta1_z_xxz_xz_1[i] * pc_y[i];

        ta1_z_xxyz_yy_0[i] = ta_xxy_yy_1[i] + ta1_z_xxy_yy_0[i] * pa_z[i] - ta1_z_xxy_yy_1[i] * pc_z[i];

        ta1_z_xxyz_yz_0[i] = ta1_z_yz_yz_0[i] * fe_0 - ta1_z_yz_yz_1[i] * fe_0 + ta1_z_xyz_yz_0[i] * pa_x[i] -
                             ta1_z_xyz_yz_1[i] * pc_x[i];

        ta1_z_xxyz_zz_0[i] = ta1_z_xxz_zz_0[i] * pa_y[i] - ta1_z_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 210-216 components of targeted buffer : GD

    auto ta1_z_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 210);

    auto ta1_z_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 211);

    auto ta1_z_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 212);

    auto ta1_z_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 213);

    auto ta1_z_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 214);

    auto ta1_z_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 215);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xy_0,   \
                             ta1_z_xx_xy_1,   \
                             ta1_z_xxz_xx_0,  \
                             ta1_z_xxz_xx_1,  \
                             ta1_z_xxz_xy_0,  \
                             ta1_z_xxz_xy_1,  \
                             ta1_z_xxzz_xx_0, \
                             ta1_z_xxzz_xy_0, \
                             ta1_z_xxzz_xz_0, \
                             ta1_z_xxzz_yy_0, \
                             ta1_z_xxzz_yz_0, \
                             ta1_z_xxzz_zz_0, \
                             ta1_z_xzz_xz_0,  \
                             ta1_z_xzz_xz_1,  \
                             ta1_z_xzz_yy_0,  \
                             ta1_z_xzz_yy_1,  \
                             ta1_z_xzz_yz_0,  \
                             ta1_z_xzz_yz_1,  \
                             ta1_z_xzz_z_0,   \
                             ta1_z_xzz_z_1,   \
                             ta1_z_xzz_zz_0,  \
                             ta1_z_xzz_zz_1,  \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_yy_0,   \
                             ta1_z_zz_yy_1,   \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta_xxz_xx_1,     \
                             ta_xxz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_xx_0[i] = ta1_z_xx_xx_0[i] * fe_0 - ta1_z_xx_xx_1[i] * fe_0 + ta_xxz_xx_1[i] +
                             ta1_z_xxz_xx_0[i] * pa_z[i] - ta1_z_xxz_xx_1[i] * pc_z[i];

        ta1_z_xxzz_xy_0[i] = ta1_z_xx_xy_0[i] * fe_0 - ta1_z_xx_xy_1[i] * fe_0 + ta_xxz_xy_1[i] +
                             ta1_z_xxz_xy_0[i] * pa_z[i] - ta1_z_xxz_xy_1[i] * pc_z[i];

        ta1_z_xxzz_xz_0[i] = ta1_z_zz_xz_0[i] * fe_0 - ta1_z_zz_xz_1[i] * fe_0 + ta1_z_xzz_z_0[i] * fe_0 -
                             ta1_z_xzz_z_1[i] * fe_0 + ta1_z_xzz_xz_0[i] * pa_x[i] - ta1_z_xzz_xz_1[i] * pc_x[i];

        ta1_z_xxzz_yy_0[i] = ta1_z_zz_yy_0[i] * fe_0 - ta1_z_zz_yy_1[i] * fe_0 + ta1_z_xzz_yy_0[i] * pa_x[i] -
                             ta1_z_xzz_yy_1[i] * pc_x[i];

        ta1_z_xxzz_yz_0[i] = ta1_z_zz_yz_0[i] * fe_0 - ta1_z_zz_yz_1[i] * fe_0 + ta1_z_xzz_yz_0[i] * pa_x[i] -
                             ta1_z_xzz_yz_1[i] * pc_x[i];

        ta1_z_xxzz_zz_0[i] = ta1_z_zz_zz_0[i] * fe_0 - ta1_z_zz_zz_1[i] * fe_0 + ta1_z_xzz_zz_0[i] * pa_x[i] -
                             ta1_z_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 216-222 components of targeted buffer : GD

    auto ta1_z_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 216);

    auto ta1_z_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 217);

    auto ta1_z_xyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 218);

    auto ta1_z_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 219);

    auto ta1_z_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 220);

    auto ta1_z_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 221);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xyyy_xx_0, \
                             ta1_z_xyyy_xy_0, \
                             ta1_z_xyyy_xz_0, \
                             ta1_z_xyyy_yy_0, \
                             ta1_z_xyyy_yz_0, \
                             ta1_z_xyyy_zz_0, \
                             ta1_z_yyy_x_0,   \
                             ta1_z_yyy_x_1,   \
                             ta1_z_yyy_xx_0,  \
                             ta1_z_yyy_xx_1,  \
                             ta1_z_yyy_xy_0,  \
                             ta1_z_yyy_xy_1,  \
                             ta1_z_yyy_xz_0,  \
                             ta1_z_yyy_xz_1,  \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_yy_0,  \
                             ta1_z_yyy_yy_1,  \
                             ta1_z_yyy_yz_0,  \
                             ta1_z_yyy_yz_1,  \
                             ta1_z_yyy_z_0,   \
                             ta1_z_yyy_z_1,   \
                             ta1_z_yyy_zz_0,  \
                             ta1_z_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_xx_0[i] = 2.0 * ta1_z_yyy_x_0[i] * fe_0 - 2.0 * ta1_z_yyy_x_1[i] * fe_0 +
                             ta1_z_yyy_xx_0[i] * pa_x[i] - ta1_z_yyy_xx_1[i] * pc_x[i];

        ta1_z_xyyy_xy_0[i] = ta1_z_yyy_y_0[i] * fe_0 - ta1_z_yyy_y_1[i] * fe_0 + ta1_z_yyy_xy_0[i] * pa_x[i] -
                             ta1_z_yyy_xy_1[i] * pc_x[i];

        ta1_z_xyyy_xz_0[i] = ta1_z_yyy_z_0[i] * fe_0 - ta1_z_yyy_z_1[i] * fe_0 + ta1_z_yyy_xz_0[i] * pa_x[i] -
                             ta1_z_yyy_xz_1[i] * pc_x[i];

        ta1_z_xyyy_yy_0[i] = ta1_z_yyy_yy_0[i] * pa_x[i] - ta1_z_yyy_yy_1[i] * pc_x[i];

        ta1_z_xyyy_yz_0[i] = ta1_z_yyy_yz_0[i] * pa_x[i] - ta1_z_yyy_yz_1[i] * pc_x[i];

        ta1_z_xyyy_zz_0[i] = ta1_z_yyy_zz_0[i] * pa_x[i] - ta1_z_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 222-228 components of targeted buffer : GD

    auto ta1_z_xyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 222);

    auto ta1_z_xyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 223);

    auto ta1_z_xyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 224);

    auto ta1_z_xyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 225);

    auto ta1_z_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 226);

    auto ta1_z_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 227);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xyy_xx_0,  \
                             ta1_z_xyy_xx_1,  \
                             ta1_z_xyy_xy_0,  \
                             ta1_z_xyy_xy_1,  \
                             ta1_z_xyyz_xx_0, \
                             ta1_z_xyyz_xy_0, \
                             ta1_z_xyyz_xz_0, \
                             ta1_z_xyyz_yy_0, \
                             ta1_z_xyyz_yz_0, \
                             ta1_z_xyyz_zz_0, \
                             ta1_z_yyz_xz_0,  \
                             ta1_z_yyz_xz_1,  \
                             ta1_z_yyz_yy_0,  \
                             ta1_z_yyz_yy_1,  \
                             ta1_z_yyz_yz_0,  \
                             ta1_z_yyz_yz_1,  \
                             ta1_z_yyz_z_0,   \
                             ta1_z_yyz_z_1,   \
                             ta1_z_yyz_zz_0,  \
                             ta1_z_yyz_zz_1,  \
                             ta_xyy_xx_1,     \
                             ta_xyy_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyz_xx_0[i] = ta_xyy_xx_1[i] + ta1_z_xyy_xx_0[i] * pa_z[i] - ta1_z_xyy_xx_1[i] * pc_z[i];

        ta1_z_xyyz_xy_0[i] = ta_xyy_xy_1[i] + ta1_z_xyy_xy_0[i] * pa_z[i] - ta1_z_xyy_xy_1[i] * pc_z[i];

        ta1_z_xyyz_xz_0[i] = ta1_z_yyz_z_0[i] * fe_0 - ta1_z_yyz_z_1[i] * fe_0 + ta1_z_yyz_xz_0[i] * pa_x[i] -
                             ta1_z_yyz_xz_1[i] * pc_x[i];

        ta1_z_xyyz_yy_0[i] = ta1_z_yyz_yy_0[i] * pa_x[i] - ta1_z_yyz_yy_1[i] * pc_x[i];

        ta1_z_xyyz_yz_0[i] = ta1_z_yyz_yz_0[i] * pa_x[i] - ta1_z_yyz_yz_1[i] * pc_x[i];

        ta1_z_xyyz_zz_0[i] = ta1_z_yyz_zz_0[i] * pa_x[i] - ta1_z_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 228-234 components of targeted buffer : GD

    auto ta1_z_xyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 228);

    auto ta1_z_xyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 229);

    auto ta1_z_xyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 230);

    auto ta1_z_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 231);

    auto ta1_z_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 232);

    auto ta1_z_xyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 233);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xyzz_xx_0, \
                             ta1_z_xyzz_xy_0, \
                             ta1_z_xyzz_xz_0, \
                             ta1_z_xyzz_yy_0, \
                             ta1_z_xyzz_yz_0, \
                             ta1_z_xyzz_zz_0, \
                             ta1_z_xzz_xx_0,  \
                             ta1_z_xzz_xx_1,  \
                             ta1_z_xzz_xz_0,  \
                             ta1_z_xzz_xz_1,  \
                             ta1_z_yzz_xy_0,  \
                             ta1_z_yzz_xy_1,  \
                             ta1_z_yzz_y_0,   \
                             ta1_z_yzz_y_1,   \
                             ta1_z_yzz_yy_0,  \
                             ta1_z_yzz_yy_1,  \
                             ta1_z_yzz_yz_0,  \
                             ta1_z_yzz_yz_1,  \
                             ta1_z_yzz_zz_0,  \
                             ta1_z_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzz_xx_0[i] = ta1_z_xzz_xx_0[i] * pa_y[i] - ta1_z_xzz_xx_1[i] * pc_y[i];

        ta1_z_xyzz_xy_0[i] = ta1_z_yzz_y_0[i] * fe_0 - ta1_z_yzz_y_1[i] * fe_0 + ta1_z_yzz_xy_0[i] * pa_x[i] -
                             ta1_z_yzz_xy_1[i] * pc_x[i];

        ta1_z_xyzz_xz_0[i] = ta1_z_xzz_xz_0[i] * pa_y[i] - ta1_z_xzz_xz_1[i] * pc_y[i];

        ta1_z_xyzz_yy_0[i] = ta1_z_yzz_yy_0[i] * pa_x[i] - ta1_z_yzz_yy_1[i] * pc_x[i];

        ta1_z_xyzz_yz_0[i] = ta1_z_yzz_yz_0[i] * pa_x[i] - ta1_z_yzz_yz_1[i] * pc_x[i];

        ta1_z_xyzz_zz_0[i] = ta1_z_yzz_zz_0[i] * pa_x[i] - ta1_z_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 234-240 components of targeted buffer : GD

    auto ta1_z_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 234);

    auto ta1_z_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 235);

    auto ta1_z_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 236);

    auto ta1_z_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 237);

    auto ta1_z_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 238);

    auto ta1_z_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 239);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xzzz_xx_0, \
                             ta1_z_xzzz_xy_0, \
                             ta1_z_xzzz_xz_0, \
                             ta1_z_xzzz_yy_0, \
                             ta1_z_xzzz_yz_0, \
                             ta1_z_xzzz_zz_0, \
                             ta1_z_zzz_x_0,   \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_xx_0,  \
                             ta1_z_zzz_xx_1,  \
                             ta1_z_zzz_xy_0,  \
                             ta1_z_zzz_xy_1,  \
                             ta1_z_zzz_xz_0,  \
                             ta1_z_zzz_xz_1,  \
                             ta1_z_zzz_y_0,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_yy_0,  \
                             ta1_z_zzz_yy_1,  \
                             ta1_z_zzz_yz_0,  \
                             ta1_z_zzz_yz_1,  \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta1_z_zzz_zz_0,  \
                             ta1_z_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_xx_0[i] = 2.0 * ta1_z_zzz_x_0[i] * fe_0 - 2.0 * ta1_z_zzz_x_1[i] * fe_0 +
                             ta1_z_zzz_xx_0[i] * pa_x[i] - ta1_z_zzz_xx_1[i] * pc_x[i];

        ta1_z_xzzz_xy_0[i] = ta1_z_zzz_y_0[i] * fe_0 - ta1_z_zzz_y_1[i] * fe_0 + ta1_z_zzz_xy_0[i] * pa_x[i] -
                             ta1_z_zzz_xy_1[i] * pc_x[i];

        ta1_z_xzzz_xz_0[i] = ta1_z_zzz_z_0[i] * fe_0 - ta1_z_zzz_z_1[i] * fe_0 + ta1_z_zzz_xz_0[i] * pa_x[i] -
                             ta1_z_zzz_xz_1[i] * pc_x[i];

        ta1_z_xzzz_yy_0[i] = ta1_z_zzz_yy_0[i] * pa_x[i] - ta1_z_zzz_yy_1[i] * pc_x[i];

        ta1_z_xzzz_yz_0[i] = ta1_z_zzz_yz_0[i] * pa_x[i] - ta1_z_zzz_yz_1[i] * pc_x[i];

        ta1_z_xzzz_zz_0[i] = ta1_z_zzz_zz_0[i] * pa_x[i] - ta1_z_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 240-246 components of targeted buffer : GD

    auto ta1_z_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 240);

    auto ta1_z_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 241);

    auto ta1_z_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 242);

    auto ta1_z_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 243);

    auto ta1_z_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 244);

    auto ta1_z_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 245);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yy_xx_0,   \
                             ta1_z_yy_xx_1,   \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_xz_0,   \
                             ta1_z_yy_xz_1,   \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yy_yz_0,   \
                             ta1_z_yy_yz_1,   \
                             ta1_z_yy_zz_0,   \
                             ta1_z_yy_zz_1,   \
                             ta1_z_yyy_x_0,   \
                             ta1_z_yyy_x_1,   \
                             ta1_z_yyy_xx_0,  \
                             ta1_z_yyy_xx_1,  \
                             ta1_z_yyy_xy_0,  \
                             ta1_z_yyy_xy_1,  \
                             ta1_z_yyy_xz_0,  \
                             ta1_z_yyy_xz_1,  \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_yy_0,  \
                             ta1_z_yyy_yy_1,  \
                             ta1_z_yyy_yz_0,  \
                             ta1_z_yyy_yz_1,  \
                             ta1_z_yyy_z_0,   \
                             ta1_z_yyy_z_1,   \
                             ta1_z_yyy_zz_0,  \
                             ta1_z_yyy_zz_1,  \
                             ta1_z_yyyy_xx_0, \
                             ta1_z_yyyy_xy_0, \
                             ta1_z_yyyy_xz_0, \
                             ta1_z_yyyy_yy_0, \
                             ta1_z_yyyy_yz_0, \
                             ta1_z_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_xx_0[i] = 3.0 * ta1_z_yy_xx_0[i] * fe_0 - 3.0 * ta1_z_yy_xx_1[i] * fe_0 +
                             ta1_z_yyy_xx_0[i] * pa_y[i] - ta1_z_yyy_xx_1[i] * pc_y[i];

        ta1_z_yyyy_xy_0[i] = 3.0 * ta1_z_yy_xy_0[i] * fe_0 - 3.0 * ta1_z_yy_xy_1[i] * fe_0 + ta1_z_yyy_x_0[i] * fe_0 -
                             ta1_z_yyy_x_1[i] * fe_0 + ta1_z_yyy_xy_0[i] * pa_y[i] - ta1_z_yyy_xy_1[i] * pc_y[i];

        ta1_z_yyyy_xz_0[i] = 3.0 * ta1_z_yy_xz_0[i] * fe_0 - 3.0 * ta1_z_yy_xz_1[i] * fe_0 +
                             ta1_z_yyy_xz_0[i] * pa_y[i] - ta1_z_yyy_xz_1[i] * pc_y[i];

        ta1_z_yyyy_yy_0[i] = 3.0 * ta1_z_yy_yy_0[i] * fe_0 - 3.0 * ta1_z_yy_yy_1[i] * fe_0 +
                             2.0 * ta1_z_yyy_y_0[i] * fe_0 - 2.0 * ta1_z_yyy_y_1[i] * fe_0 +
                             ta1_z_yyy_yy_0[i] * pa_y[i] - ta1_z_yyy_yy_1[i] * pc_y[i];

        ta1_z_yyyy_yz_0[i] = 3.0 * ta1_z_yy_yz_0[i] * fe_0 - 3.0 * ta1_z_yy_yz_1[i] * fe_0 + ta1_z_yyy_z_0[i] * fe_0 -
                             ta1_z_yyy_z_1[i] * fe_0 + ta1_z_yyy_yz_0[i] * pa_y[i] - ta1_z_yyy_yz_1[i] * pc_y[i];

        ta1_z_yyyy_zz_0[i] = 3.0 * ta1_z_yy_zz_0[i] * fe_0 - 3.0 * ta1_z_yy_zz_1[i] * fe_0 +
                             ta1_z_yyy_zz_0[i] * pa_y[i] - ta1_z_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 246-252 components of targeted buffer : GD

    auto ta1_z_yyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 246);

    auto ta1_z_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 247);

    auto ta1_z_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 248);

    auto ta1_z_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 249);

    auto ta1_z_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 250);

    auto ta1_z_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 251);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyy_xx_0,  \
                             ta1_z_yyy_xx_1,  \
                             ta1_z_yyy_xy_0,  \
                             ta1_z_yyy_xy_1,  \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_yy_0,  \
                             ta1_z_yyy_yy_1,  \
                             ta1_z_yyy_yz_0,  \
                             ta1_z_yyy_yz_1,  \
                             ta1_z_yyyz_xx_0, \
                             ta1_z_yyyz_xy_0, \
                             ta1_z_yyyz_xz_0, \
                             ta1_z_yyyz_yy_0, \
                             ta1_z_yyyz_yz_0, \
                             ta1_z_yyyz_zz_0, \
                             ta1_z_yyz_xz_0,  \
                             ta1_z_yyz_xz_1,  \
                             ta1_z_yyz_zz_0,  \
                             ta1_z_yyz_zz_1,  \
                             ta1_z_yz_xz_0,   \
                             ta1_z_yz_xz_1,   \
                             ta1_z_yz_zz_0,   \
                             ta1_z_yz_zz_1,   \
                             ta_yyy_xx_1,     \
                             ta_yyy_xy_1,     \
                             ta_yyy_yy_1,     \
                             ta_yyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_xx_0[i] = ta_yyy_xx_1[i] + ta1_z_yyy_xx_0[i] * pa_z[i] - ta1_z_yyy_xx_1[i] * pc_z[i];

        ta1_z_yyyz_xy_0[i] = ta_yyy_xy_1[i] + ta1_z_yyy_xy_0[i] * pa_z[i] - ta1_z_yyy_xy_1[i] * pc_z[i];

        ta1_z_yyyz_xz_0[i] = 2.0 * ta1_z_yz_xz_0[i] * fe_0 - 2.0 * ta1_z_yz_xz_1[i] * fe_0 +
                             ta1_z_yyz_xz_0[i] * pa_y[i] - ta1_z_yyz_xz_1[i] * pc_y[i];

        ta1_z_yyyz_yy_0[i] = ta_yyy_yy_1[i] + ta1_z_yyy_yy_0[i] * pa_z[i] - ta1_z_yyy_yy_1[i] * pc_z[i];

        ta1_z_yyyz_yz_0[i] = ta1_z_yyy_y_0[i] * fe_0 - ta1_z_yyy_y_1[i] * fe_0 + ta_yyy_yz_1[i] +
                             ta1_z_yyy_yz_0[i] * pa_z[i] - ta1_z_yyy_yz_1[i] * pc_z[i];

        ta1_z_yyyz_zz_0[i] = 2.0 * ta1_z_yz_zz_0[i] * fe_0 - 2.0 * ta1_z_yz_zz_1[i] * fe_0 +
                             ta1_z_yyz_zz_0[i] * pa_y[i] - ta1_z_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 252-258 components of targeted buffer : GD

    auto ta1_z_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 252);

    auto ta1_z_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 253);

    auto ta1_z_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 254);

    auto ta1_z_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 255);

    auto ta1_z_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 256);

    auto ta1_z_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 257);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yyz_xy_0,  \
                             ta1_z_yyz_xy_1,  \
                             ta1_z_yyz_yy_0,  \
                             ta1_z_yyz_yy_1,  \
                             ta1_z_yyzz_xx_0, \
                             ta1_z_yyzz_xy_0, \
                             ta1_z_yyzz_xz_0, \
                             ta1_z_yyzz_yy_0, \
                             ta1_z_yyzz_yz_0, \
                             ta1_z_yyzz_zz_0, \
                             ta1_z_yzz_xx_0,  \
                             ta1_z_yzz_xx_1,  \
                             ta1_z_yzz_xz_0,  \
                             ta1_z_yzz_xz_1,  \
                             ta1_z_yzz_yz_0,  \
                             ta1_z_yzz_yz_1,  \
                             ta1_z_yzz_z_0,   \
                             ta1_z_yzz_z_1,   \
                             ta1_z_yzz_zz_0,  \
                             ta1_z_yzz_zz_1,  \
                             ta1_z_zz_xx_0,   \
                             ta1_z_zz_xx_1,   \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta_yyz_xy_1,     \
                             ta_yyz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_xx_0[i] = ta1_z_zz_xx_0[i] * fe_0 - ta1_z_zz_xx_1[i] * fe_0 + ta1_z_yzz_xx_0[i] * pa_y[i] -
                             ta1_z_yzz_xx_1[i] * pc_y[i];

        ta1_z_yyzz_xy_0[i] = ta1_z_yy_xy_0[i] * fe_0 - ta1_z_yy_xy_1[i] * fe_0 + ta_yyz_xy_1[i] +
                             ta1_z_yyz_xy_0[i] * pa_z[i] - ta1_z_yyz_xy_1[i] * pc_z[i];

        ta1_z_yyzz_xz_0[i] = ta1_z_zz_xz_0[i] * fe_0 - ta1_z_zz_xz_1[i] * fe_0 + ta1_z_yzz_xz_0[i] * pa_y[i] -
                             ta1_z_yzz_xz_1[i] * pc_y[i];

        ta1_z_yyzz_yy_0[i] = ta1_z_yy_yy_0[i] * fe_0 - ta1_z_yy_yy_1[i] * fe_0 + ta_yyz_yy_1[i] +
                             ta1_z_yyz_yy_0[i] * pa_z[i] - ta1_z_yyz_yy_1[i] * pc_z[i];

        ta1_z_yyzz_yz_0[i] = ta1_z_zz_yz_0[i] * fe_0 - ta1_z_zz_yz_1[i] * fe_0 + ta1_z_yzz_z_0[i] * fe_0 -
                             ta1_z_yzz_z_1[i] * fe_0 + ta1_z_yzz_yz_0[i] * pa_y[i] - ta1_z_yzz_yz_1[i] * pc_y[i];

        ta1_z_yyzz_zz_0[i] = ta1_z_zz_zz_0[i] * fe_0 - ta1_z_zz_zz_1[i] * fe_0 + ta1_z_yzz_zz_0[i] * pa_y[i] -
                             ta1_z_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 258-264 components of targeted buffer : GD

    auto ta1_z_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 258);

    auto ta1_z_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 259);

    auto ta1_z_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 260);

    auto ta1_z_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 261);

    auto ta1_z_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 262);

    auto ta1_z_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 263);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yzzz_xx_0, \
                             ta1_z_yzzz_xy_0, \
                             ta1_z_yzzz_xz_0, \
                             ta1_z_yzzz_yy_0, \
                             ta1_z_yzzz_yz_0, \
                             ta1_z_yzzz_zz_0, \
                             ta1_z_zzz_x_0,   \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_xx_0,  \
                             ta1_z_zzz_xx_1,  \
                             ta1_z_zzz_xy_0,  \
                             ta1_z_zzz_xy_1,  \
                             ta1_z_zzz_xz_0,  \
                             ta1_z_zzz_xz_1,  \
                             ta1_z_zzz_y_0,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_yy_0,  \
                             ta1_z_zzz_yy_1,  \
                             ta1_z_zzz_yz_0,  \
                             ta1_z_zzz_yz_1,  \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta1_z_zzz_zz_0,  \
                             ta1_z_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_xx_0[i] = ta1_z_zzz_xx_0[i] * pa_y[i] - ta1_z_zzz_xx_1[i] * pc_y[i];

        ta1_z_yzzz_xy_0[i] = ta1_z_zzz_x_0[i] * fe_0 - ta1_z_zzz_x_1[i] * fe_0 + ta1_z_zzz_xy_0[i] * pa_y[i] -
                             ta1_z_zzz_xy_1[i] * pc_y[i];

        ta1_z_yzzz_xz_0[i] = ta1_z_zzz_xz_0[i] * pa_y[i] - ta1_z_zzz_xz_1[i] * pc_y[i];

        ta1_z_yzzz_yy_0[i] = 2.0 * ta1_z_zzz_y_0[i] * fe_0 - 2.0 * ta1_z_zzz_y_1[i] * fe_0 +
                             ta1_z_zzz_yy_0[i] * pa_y[i] - ta1_z_zzz_yy_1[i] * pc_y[i];

        ta1_z_yzzz_yz_0[i] = ta1_z_zzz_z_0[i] * fe_0 - ta1_z_zzz_z_1[i] * fe_0 + ta1_z_zzz_yz_0[i] * pa_y[i] -
                             ta1_z_zzz_yz_1[i] * pc_y[i];

        ta1_z_yzzz_zz_0[i] = ta1_z_zzz_zz_0[i] * pa_y[i] - ta1_z_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 264-270 components of targeted buffer : GD

    auto ta1_z_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 264);

    auto ta1_z_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 265);

    auto ta1_z_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 266);

    auto ta1_z_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 267);

    auto ta1_z_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 268);

    auto ta1_z_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 269);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_zz_xx_0,   \
                             ta1_z_zz_xx_1,   \
                             ta1_z_zz_xy_0,   \
                             ta1_z_zz_xy_1,   \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_yy_0,   \
                             ta1_z_zz_yy_1,   \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta1_z_zzz_x_0,   \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_xx_0,  \
                             ta1_z_zzz_xx_1,  \
                             ta1_z_zzz_xy_0,  \
                             ta1_z_zzz_xy_1,  \
                             ta1_z_zzz_xz_0,  \
                             ta1_z_zzz_xz_1,  \
                             ta1_z_zzz_y_0,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_yy_0,  \
                             ta1_z_zzz_yy_1,  \
                             ta1_z_zzz_yz_0,  \
                             ta1_z_zzz_yz_1,  \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta1_z_zzz_zz_0,  \
                             ta1_z_zzz_zz_1,  \
                             ta1_z_zzzz_xx_0, \
                             ta1_z_zzzz_xy_0, \
                             ta1_z_zzzz_xz_0, \
                             ta1_z_zzzz_yy_0, \
                             ta1_z_zzzz_yz_0, \
                             ta1_z_zzzz_zz_0, \
                             ta_zzz_xx_1,     \
                             ta_zzz_xy_1,     \
                             ta_zzz_xz_1,     \
                             ta_zzz_yy_1,     \
                             ta_zzz_yz_1,     \
                             ta_zzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_xx_0[i] = 3.0 * ta1_z_zz_xx_0[i] * fe_0 - 3.0 * ta1_z_zz_xx_1[i] * fe_0 + ta_zzz_xx_1[i] +
                             ta1_z_zzz_xx_0[i] * pa_z[i] - ta1_z_zzz_xx_1[i] * pc_z[i];

        ta1_z_zzzz_xy_0[i] = 3.0 * ta1_z_zz_xy_0[i] * fe_0 - 3.0 * ta1_z_zz_xy_1[i] * fe_0 + ta_zzz_xy_1[i] +
                             ta1_z_zzz_xy_0[i] * pa_z[i] - ta1_z_zzz_xy_1[i] * pc_z[i];

        ta1_z_zzzz_xz_0[i] = 3.0 * ta1_z_zz_xz_0[i] * fe_0 - 3.0 * ta1_z_zz_xz_1[i] * fe_0 + ta1_z_zzz_x_0[i] * fe_0 -
                             ta1_z_zzz_x_1[i] * fe_0 + ta_zzz_xz_1[i] + ta1_z_zzz_xz_0[i] * pa_z[i] -
                             ta1_z_zzz_xz_1[i] * pc_z[i];

        ta1_z_zzzz_yy_0[i] = 3.0 * ta1_z_zz_yy_0[i] * fe_0 - 3.0 * ta1_z_zz_yy_1[i] * fe_0 + ta_zzz_yy_1[i] +
                             ta1_z_zzz_yy_0[i] * pa_z[i] - ta1_z_zzz_yy_1[i] * pc_z[i];

        ta1_z_zzzz_yz_0[i] = 3.0 * ta1_z_zz_yz_0[i] * fe_0 - 3.0 * ta1_z_zz_yz_1[i] * fe_0 + ta1_z_zzz_y_0[i] * fe_0 -
                             ta1_z_zzz_y_1[i] * fe_0 + ta_zzz_yz_1[i] + ta1_z_zzz_yz_0[i] * pa_z[i] -
                             ta1_z_zzz_yz_1[i] * pc_z[i];

        ta1_z_zzzz_zz_0[i] = 3.0 * ta1_z_zz_zz_0[i] * fe_0 - 3.0 * ta1_z_zz_zz_1[i] * fe_0 +
                             2.0 * ta1_z_zzz_z_0[i] * fe_0 - 2.0 * ta1_z_zzz_z_1[i] * fe_0 + ta_zzz_zz_1[i] +
                             ta1_z_zzz_zz_0[i] * pa_z[i] - ta1_z_zzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
