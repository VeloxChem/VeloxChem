#include "NuclearPotentialGeom010PrimRecFF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_ff(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_ff,
                                        const size_t              idx_npot_geom_010_0_pf,
                                        const size_t              idx_npot_geom_010_1_pf,
                                        const size_t              idx_npot_geom_010_0_dd,
                                        const size_t              idx_npot_geom_010_1_dd,
                                        const size_t              idx_npot_1_df,
                                        const size_t              idx_npot_geom_010_0_df,
                                        const size_t              idx_npot_geom_010_1_df,
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

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf);

    auto ta1_x_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 1);

    auto ta1_x_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 2);

    auto ta1_x_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 3);

    auto ta1_x_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 4);

    auto ta1_x_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 5);

    auto ta1_x_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 6);

    auto ta1_x_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 7);

    auto ta1_x_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 8);

    auto ta1_x_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 9);

    auto ta1_x_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 10);

    auto ta1_x_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 11);

    auto ta1_x_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 12);

    auto ta1_x_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 13);

    auto ta1_x_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 14);

    auto ta1_x_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 15);

    auto ta1_x_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 16);

    auto ta1_x_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 17);

    auto ta1_x_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 18);

    auto ta1_x_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 19);

    auto ta1_x_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 20);

    auto ta1_x_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 21);

    auto ta1_x_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 22);

    auto ta1_x_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 23);

    auto ta1_x_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 24);

    auto ta1_x_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 25);

    auto ta1_x_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 26);

    auto ta1_x_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 27);

    auto ta1_x_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 28);

    auto ta1_x_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 29);

    auto ta1_y_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 30);

    auto ta1_y_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 31);

    auto ta1_y_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 32);

    auto ta1_y_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 33);

    auto ta1_y_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 34);

    auto ta1_y_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 35);

    auto ta1_y_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 36);

    auto ta1_y_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 37);

    auto ta1_y_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 38);

    auto ta1_y_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 39);

    auto ta1_y_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 40);

    auto ta1_y_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 41);

    auto ta1_y_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 42);

    auto ta1_y_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 43);

    auto ta1_y_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 44);

    auto ta1_y_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 45);

    auto ta1_y_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 46);

    auto ta1_y_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 47);

    auto ta1_y_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 48);

    auto ta1_y_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 49);

    auto ta1_y_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 50);

    auto ta1_y_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 51);

    auto ta1_y_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 52);

    auto ta1_y_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 53);

    auto ta1_y_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 54);

    auto ta1_y_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 55);

    auto ta1_y_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 56);

    auto ta1_y_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 57);

    auto ta1_y_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 58);

    auto ta1_y_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 59);

    auto ta1_z_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 60);

    auto ta1_z_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 61);

    auto ta1_z_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 62);

    auto ta1_z_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 63);

    auto ta1_z_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 64);

    auto ta1_z_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 65);

    auto ta1_z_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 66);

    auto ta1_z_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 67);

    auto ta1_z_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 68);

    auto ta1_z_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 69);

    auto ta1_z_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 70);

    auto ta1_z_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 71);

    auto ta1_z_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 72);

    auto ta1_z_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 73);

    auto ta1_z_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 74);

    auto ta1_z_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 75);

    auto ta1_z_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 76);

    auto ta1_z_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 77);

    auto ta1_z_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 78);

    auto ta1_z_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 79);

    auto ta1_z_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 80);

    auto ta1_z_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 81);

    auto ta1_z_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 82);

    auto ta1_z_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 83);

    auto ta1_z_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 84);

    auto ta1_z_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 85);

    auto ta1_z_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 86);

    auto ta1_z_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 87);

    auto ta1_z_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 88);

    auto ta1_z_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 89);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf);

    auto ta1_x_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 1);

    auto ta1_x_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 2);

    auto ta1_x_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 3);

    auto ta1_x_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 4);

    auto ta1_x_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 5);

    auto ta1_x_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 6);

    auto ta1_x_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 7);

    auto ta1_x_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 8);

    auto ta1_x_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 9);

    auto ta1_x_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 10);

    auto ta1_x_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 11);

    auto ta1_x_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 12);

    auto ta1_x_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 13);

    auto ta1_x_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 14);

    auto ta1_x_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 15);

    auto ta1_x_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 16);

    auto ta1_x_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 17);

    auto ta1_x_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 18);

    auto ta1_x_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 19);

    auto ta1_x_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 20);

    auto ta1_x_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 21);

    auto ta1_x_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 22);

    auto ta1_x_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 23);

    auto ta1_x_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 24);

    auto ta1_x_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 25);

    auto ta1_x_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 26);

    auto ta1_x_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 27);

    auto ta1_x_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 28);

    auto ta1_x_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 29);

    auto ta1_y_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 30);

    auto ta1_y_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 31);

    auto ta1_y_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 32);

    auto ta1_y_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 33);

    auto ta1_y_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 34);

    auto ta1_y_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 35);

    auto ta1_y_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 36);

    auto ta1_y_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 37);

    auto ta1_y_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 38);

    auto ta1_y_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 39);

    auto ta1_y_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 40);

    auto ta1_y_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 41);

    auto ta1_y_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 42);

    auto ta1_y_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 43);

    auto ta1_y_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 44);

    auto ta1_y_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 45);

    auto ta1_y_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 46);

    auto ta1_y_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 47);

    auto ta1_y_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 48);

    auto ta1_y_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 49);

    auto ta1_y_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 50);

    auto ta1_y_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 51);

    auto ta1_y_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 52);

    auto ta1_y_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 53);

    auto ta1_y_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 54);

    auto ta1_y_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 55);

    auto ta1_y_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 56);

    auto ta1_y_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 57);

    auto ta1_y_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 58);

    auto ta1_y_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 59);

    auto ta1_z_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 60);

    auto ta1_z_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 61);

    auto ta1_z_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 62);

    auto ta1_z_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 63);

    auto ta1_z_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 64);

    auto ta1_z_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 65);

    auto ta1_z_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 66);

    auto ta1_z_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 67);

    auto ta1_z_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 68);

    auto ta1_z_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 69);

    auto ta1_z_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 70);

    auto ta1_z_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 71);

    auto ta1_z_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 72);

    auto ta1_z_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 73);

    auto ta1_z_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 74);

    auto ta1_z_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 75);

    auto ta1_z_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 76);

    auto ta1_z_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 77);

    auto ta1_z_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 78);

    auto ta1_z_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 79);

    auto ta1_z_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 80);

    auto ta1_z_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 81);

    auto ta1_z_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 82);

    auto ta1_z_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 83);

    auto ta1_z_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 84);

    auto ta1_z_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 85);

    auto ta1_z_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 86);

    auto ta1_z_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 87);

    auto ta1_z_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 88);

    auto ta1_z_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 89);

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd);

    auto ta1_x_xx_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 1);

    auto ta1_x_xx_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 2);

    auto ta1_x_xx_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 3);

    auto ta1_x_xx_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 4);

    auto ta1_x_xx_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 5);

    auto ta1_x_xz_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 14);

    auto ta1_x_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 18);

    auto ta1_x_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 19);

    auto ta1_x_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 20);

    auto ta1_x_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 21);

    auto ta1_x_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 22);

    auto ta1_x_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 23);

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

    auto ta1_y_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 54);

    auto ta1_y_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 55);

    auto ta1_y_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 56);

    auto ta1_y_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 57);

    auto ta1_y_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 58);

    auto ta1_y_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 59);

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

    auto ta1_z_yy_xx_0 = pbuffer.data(idx_npot_geom_010_0_dd + 90);

    auto ta1_z_yy_xy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 91);

    auto ta1_z_yy_xz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 92);

    auto ta1_z_yy_yy_0 = pbuffer.data(idx_npot_geom_010_0_dd + 93);

    auto ta1_z_yy_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 94);

    auto ta1_z_yy_zz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 95);

    auto ta1_z_yz_yz_0 = pbuffer.data(idx_npot_geom_010_0_dd + 100);

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

    auto ta1_x_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 14);

    auto ta1_x_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 18);

    auto ta1_x_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 19);

    auto ta1_x_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 20);

    auto ta1_x_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 21);

    auto ta1_x_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 22);

    auto ta1_x_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 23);

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

    auto ta1_y_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 54);

    auto ta1_y_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 55);

    auto ta1_y_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 56);

    auto ta1_y_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 57);

    auto ta1_y_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 58);

    auto ta1_y_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 59);

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

    auto ta1_z_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 90);

    auto ta1_z_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 91);

    auto ta1_z_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 92);

    auto ta1_z_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 93);

    auto ta1_z_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 94);

    auto ta1_z_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 95);

    auto ta1_z_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 100);

    auto ta1_z_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 102);

    auto ta1_z_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 103);

    auto ta1_z_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 104);

    auto ta1_z_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 105);

    auto ta1_z_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 106);

    auto ta1_z_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 107);

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_1 = pbuffer.data(idx_npot_1_df);

    auto ta_xx_xxy_1 = pbuffer.data(idx_npot_1_df + 1);

    auto ta_xx_xxz_1 = pbuffer.data(idx_npot_1_df + 2);

    auto ta_xx_xyy_1 = pbuffer.data(idx_npot_1_df + 3);

    auto ta_xx_xyz_1 = pbuffer.data(idx_npot_1_df + 4);

    auto ta_xx_xzz_1 = pbuffer.data(idx_npot_1_df + 5);

    auto ta_xx_yyy_1 = pbuffer.data(idx_npot_1_df + 6);

    auto ta_xx_yyz_1 = pbuffer.data(idx_npot_1_df + 7);

    auto ta_xx_yzz_1 = pbuffer.data(idx_npot_1_df + 8);

    auto ta_xx_zzz_1 = pbuffer.data(idx_npot_1_df + 9);

    auto ta_xy_xxy_1 = pbuffer.data(idx_npot_1_df + 11);

    auto ta_xy_xyy_1 = pbuffer.data(idx_npot_1_df + 13);

    auto ta_xz_xxz_1 = pbuffer.data(idx_npot_1_df + 22);

    auto ta_xz_xzz_1 = pbuffer.data(idx_npot_1_df + 25);

    auto ta_yy_xxx_1 = pbuffer.data(idx_npot_1_df + 30);

    auto ta_yy_xxy_1 = pbuffer.data(idx_npot_1_df + 31);

    auto ta_yy_xxz_1 = pbuffer.data(idx_npot_1_df + 32);

    auto ta_yy_xyy_1 = pbuffer.data(idx_npot_1_df + 33);

    auto ta_yy_xyz_1 = pbuffer.data(idx_npot_1_df + 34);

    auto ta_yy_xzz_1 = pbuffer.data(idx_npot_1_df + 35);

    auto ta_yy_yyy_1 = pbuffer.data(idx_npot_1_df + 36);

    auto ta_yy_yyz_1 = pbuffer.data(idx_npot_1_df + 37);

    auto ta_yy_yzz_1 = pbuffer.data(idx_npot_1_df + 38);

    auto ta_yy_zzz_1 = pbuffer.data(idx_npot_1_df + 39);

    auto ta_yz_yyz_1 = pbuffer.data(idx_npot_1_df + 47);

    auto ta_yz_yzz_1 = pbuffer.data(idx_npot_1_df + 48);

    auto ta_zz_xxx_1 = pbuffer.data(idx_npot_1_df + 50);

    auto ta_zz_xxy_1 = pbuffer.data(idx_npot_1_df + 51);

    auto ta_zz_xxz_1 = pbuffer.data(idx_npot_1_df + 52);

    auto ta_zz_xyy_1 = pbuffer.data(idx_npot_1_df + 53);

    auto ta_zz_xyz_1 = pbuffer.data(idx_npot_1_df + 54);

    auto ta_zz_xzz_1 = pbuffer.data(idx_npot_1_df + 55);

    auto ta_zz_yyy_1 = pbuffer.data(idx_npot_1_df + 56);

    auto ta_zz_yyz_1 = pbuffer.data(idx_npot_1_df + 57);

    auto ta_zz_yzz_1 = pbuffer.data(idx_npot_1_df + 58);

    auto ta_zz_zzz_1 = pbuffer.data(idx_npot_1_df + 59);

    // Set up components of auxiliary buffer : DF

    auto ta1_x_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df);

    auto ta1_x_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 1);

    auto ta1_x_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 2);

    auto ta1_x_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 3);

    auto ta1_x_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 4);

    auto ta1_x_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 5);

    auto ta1_x_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 6);

    auto ta1_x_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 7);

    auto ta1_x_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 8);

    auto ta1_x_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 9);

    auto ta1_x_xy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 10);

    auto ta1_x_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 11);

    auto ta1_x_xy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 12);

    auto ta1_x_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 13);

    auto ta1_x_xy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 15);

    auto ta1_x_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 16);

    auto ta1_x_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 20);

    auto ta1_x_xz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 21);

    auto ta1_x_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 22);

    auto ta1_x_xz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 23);

    auto ta1_x_xz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 24);

    auto ta1_x_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 25);

    auto ta1_x_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 29);

    auto ta1_x_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 30);

    auto ta1_x_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 31);

    auto ta1_x_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 32);

    auto ta1_x_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 33);

    auto ta1_x_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 34);

    auto ta1_x_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 35);

    auto ta1_x_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 36);

    auto ta1_x_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 37);

    auto ta1_x_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 38);

    auto ta1_x_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 39);

    auto ta1_x_yz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 42);

    auto ta1_x_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 45);

    auto ta1_x_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 47);

    auto ta1_x_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 48);

    auto ta1_x_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 49);

    auto ta1_x_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 50);

    auto ta1_x_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 51);

    auto ta1_x_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 52);

    auto ta1_x_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 53);

    auto ta1_x_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 54);

    auto ta1_x_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 55);

    auto ta1_x_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 56);

    auto ta1_x_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 57);

    auto ta1_x_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 58);

    auto ta1_x_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 59);

    auto ta1_y_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 60);

    auto ta1_y_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 61);

    auto ta1_y_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 62);

    auto ta1_y_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 63);

    auto ta1_y_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 64);

    auto ta1_y_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 65);

    auto ta1_y_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 66);

    auto ta1_y_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 67);

    auto ta1_y_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 68);

    auto ta1_y_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 69);

    auto ta1_y_xy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 70);

    auto ta1_y_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 71);

    auto ta1_y_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 73);

    auto ta1_y_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 76);

    auto ta1_y_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 77);

    auto ta1_y_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 78);

    auto ta1_y_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 82);

    auto ta1_y_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 85);

    auto ta1_y_xz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 87);

    auto ta1_y_xz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 88);

    auto ta1_y_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 89);

    auto ta1_y_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 90);

    auto ta1_y_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 91);

    auto ta1_y_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 92);

    auto ta1_y_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 93);

    auto ta1_y_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 94);

    auto ta1_y_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 95);

    auto ta1_y_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 96);

    auto ta1_y_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 97);

    auto ta1_y_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 98);

    auto ta1_y_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 99);

    auto ta1_y_yz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 101);

    auto ta1_y_yz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 103);

    auto ta1_y_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 104);

    auto ta1_y_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 106);

    auto ta1_y_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 107);

    auto ta1_y_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 108);

    auto ta1_y_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 109);

    auto ta1_y_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 110);

    auto ta1_y_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 111);

    auto ta1_y_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 112);

    auto ta1_y_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 113);

    auto ta1_y_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 114);

    auto ta1_y_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 115);

    auto ta1_y_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 116);

    auto ta1_y_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 117);

    auto ta1_y_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 118);

    auto ta1_y_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 119);

    auto ta1_z_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 120);

    auto ta1_z_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 121);

    auto ta1_z_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 122);

    auto ta1_z_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 123);

    auto ta1_z_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 124);

    auto ta1_z_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 125);

    auto ta1_z_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 126);

    auto ta1_z_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 127);

    auto ta1_z_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 128);

    auto ta1_z_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 129);

    auto ta1_z_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 131);

    auto ta1_z_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 133);

    auto ta1_z_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 136);

    auto ta1_z_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 137);

    auto ta1_z_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 138);

    auto ta1_z_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 140);

    auto ta1_z_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 142);

    auto ta1_z_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 145);

    auto ta1_z_xz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 147);

    auto ta1_z_xz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 148);

    auto ta1_z_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 149);

    auto ta1_z_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 150);

    auto ta1_z_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 151);

    auto ta1_z_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 152);

    auto ta1_z_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 153);

    auto ta1_z_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 154);

    auto ta1_z_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 155);

    auto ta1_z_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 156);

    auto ta1_z_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 157);

    auto ta1_z_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 158);

    auto ta1_z_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 159);

    auto ta1_z_yz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 162);

    auto ta1_z_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 164);

    auto ta1_z_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 165);

    auto ta1_z_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 166);

    auto ta1_z_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 167);

    auto ta1_z_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 168);

    auto ta1_z_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 169);

    auto ta1_z_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 170);

    auto ta1_z_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 171);

    auto ta1_z_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 172);

    auto ta1_z_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 173);

    auto ta1_z_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 174);

    auto ta1_z_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 175);

    auto ta1_z_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 176);

    auto ta1_z_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 177);

    auto ta1_z_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 178);

    auto ta1_z_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 179);

    // Set up components of auxiliary buffer : DF

    auto ta1_x_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df);

    auto ta1_x_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 1);

    auto ta1_x_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 2);

    auto ta1_x_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 3);

    auto ta1_x_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 4);

    auto ta1_x_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 5);

    auto ta1_x_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 6);

    auto ta1_x_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 7);

    auto ta1_x_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 8);

    auto ta1_x_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 9);

    auto ta1_x_xy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 10);

    auto ta1_x_xy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 11);

    auto ta1_x_xy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 12);

    auto ta1_x_xy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 13);

    auto ta1_x_xy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 15);

    auto ta1_x_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 16);

    auto ta1_x_xz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 20);

    auto ta1_x_xz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 21);

    auto ta1_x_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 22);

    auto ta1_x_xz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 23);

    auto ta1_x_xz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 24);

    auto ta1_x_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 25);

    auto ta1_x_xz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 29);

    auto ta1_x_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 30);

    auto ta1_x_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 31);

    auto ta1_x_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 32);

    auto ta1_x_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 33);

    auto ta1_x_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 34);

    auto ta1_x_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 35);

    auto ta1_x_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 36);

    auto ta1_x_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 37);

    auto ta1_x_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 38);

    auto ta1_x_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 39);

    auto ta1_x_yz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 42);

    auto ta1_x_yz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 45);

    auto ta1_x_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 47);

    auto ta1_x_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 48);

    auto ta1_x_yz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 49);

    auto ta1_x_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 50);

    auto ta1_x_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 51);

    auto ta1_x_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 52);

    auto ta1_x_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 53);

    auto ta1_x_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 54);

    auto ta1_x_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 55);

    auto ta1_x_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 56);

    auto ta1_x_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 57);

    auto ta1_x_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 58);

    auto ta1_x_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 59);

    auto ta1_y_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 60);

    auto ta1_y_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 61);

    auto ta1_y_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 62);

    auto ta1_y_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 63);

    auto ta1_y_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 64);

    auto ta1_y_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 65);

    auto ta1_y_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 66);

    auto ta1_y_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 67);

    auto ta1_y_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 68);

    auto ta1_y_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 69);

    auto ta1_y_xy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 70);

    auto ta1_y_xy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 71);

    auto ta1_y_xy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 73);

    auto ta1_y_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 76);

    auto ta1_y_xy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 77);

    auto ta1_y_xy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 78);

    auto ta1_y_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 82);

    auto ta1_y_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 85);

    auto ta1_y_xz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 87);

    auto ta1_y_xz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 88);

    auto ta1_y_xz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 89);

    auto ta1_y_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 90);

    auto ta1_y_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 91);

    auto ta1_y_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 92);

    auto ta1_y_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 93);

    auto ta1_y_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 94);

    auto ta1_y_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 95);

    auto ta1_y_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 96);

    auto ta1_y_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 97);

    auto ta1_y_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 98);

    auto ta1_y_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 99);

    auto ta1_y_yz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 101);

    auto ta1_y_yz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 103);

    auto ta1_y_yz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 104);

    auto ta1_y_yz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 106);

    auto ta1_y_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 107);

    auto ta1_y_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 108);

    auto ta1_y_yz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 109);

    auto ta1_y_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 110);

    auto ta1_y_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 111);

    auto ta1_y_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 112);

    auto ta1_y_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 113);

    auto ta1_y_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 114);

    auto ta1_y_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 115);

    auto ta1_y_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 116);

    auto ta1_y_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 117);

    auto ta1_y_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 118);

    auto ta1_y_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 119);

    auto ta1_z_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 120);

    auto ta1_z_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 121);

    auto ta1_z_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 122);

    auto ta1_z_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 123);

    auto ta1_z_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 124);

    auto ta1_z_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 125);

    auto ta1_z_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 126);

    auto ta1_z_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 127);

    auto ta1_z_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 128);

    auto ta1_z_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 129);

    auto ta1_z_xy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 131);

    auto ta1_z_xy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 133);

    auto ta1_z_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 136);

    auto ta1_z_xy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 137);

    auto ta1_z_xy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 138);

    auto ta1_z_xz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 140);

    auto ta1_z_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 142);

    auto ta1_z_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 145);

    auto ta1_z_xz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 147);

    auto ta1_z_xz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 148);

    auto ta1_z_xz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 149);

    auto ta1_z_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 150);

    auto ta1_z_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 151);

    auto ta1_z_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 152);

    auto ta1_z_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 153);

    auto ta1_z_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 154);

    auto ta1_z_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 155);

    auto ta1_z_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 156);

    auto ta1_z_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 157);

    auto ta1_z_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 158);

    auto ta1_z_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 159);

    auto ta1_z_yz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 162);

    auto ta1_z_yz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 164);

    auto ta1_z_yz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 165);

    auto ta1_z_yz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 166);

    auto ta1_z_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 167);

    auto ta1_z_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 168);

    auto ta1_z_yz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 169);

    auto ta1_z_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 170);

    auto ta1_z_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 171);

    auto ta1_z_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 172);

    auto ta1_z_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 173);

    auto ta1_z_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 174);

    auto ta1_z_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 175);

    auto ta1_z_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 176);

    auto ta1_z_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 177);

    auto ta1_z_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 178);

    auto ta1_z_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 179);

    // Set up 0-10 components of targeted buffer : FF

    auto ta1_x_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff);

    auto ta1_x_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 1);

    auto ta1_x_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 2);

    auto ta1_x_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 3);

    auto ta1_x_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 4);

    auto ta1_x_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 5);

    auto ta1_x_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 6);

    auto ta1_x_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 7);

    auto ta1_x_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 8);

    auto ta1_x_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 9);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxy_0,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xxz_0,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xyy_0,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_x_xyz_0,   \
                             ta1_x_x_xyz_1,   \
                             ta1_x_x_xzz_0,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_x_yyy_0,   \
                             ta1_x_x_yyy_1,   \
                             ta1_x_x_yyz_0,   \
                             ta1_x_x_yyz_1,   \
                             ta1_x_x_yzz_0,   \
                             ta1_x_x_yzz_1,   \
                             ta1_x_x_zzz_0,   \
                             ta1_x_x_zzz_1,   \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xxx_0,  \
                             ta1_x_xx_xxx_1,  \
                             ta1_x_xx_xxy_0,  \
                             ta1_x_xx_xxy_1,  \
                             ta1_x_xx_xxz_0,  \
                             ta1_x_xx_xxz_1,  \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xyy_0,  \
                             ta1_x_xx_xyy_1,  \
                             ta1_x_xx_xyz_0,  \
                             ta1_x_xx_xyz_1,  \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_xzz_0,  \
                             ta1_x_xx_xzz_1,  \
                             ta1_x_xx_yy_0,   \
                             ta1_x_xx_yy_1,   \
                             ta1_x_xx_yyy_0,  \
                             ta1_x_xx_yyy_1,  \
                             ta1_x_xx_yyz_0,  \
                             ta1_x_xx_yyz_1,  \
                             ta1_x_xx_yz_0,   \
                             ta1_x_xx_yz_1,   \
                             ta1_x_xx_yzz_0,  \
                             ta1_x_xx_yzz_1,  \
                             ta1_x_xx_zz_0,   \
                             ta1_x_xx_zz_1,   \
                             ta1_x_xx_zzz_0,  \
                             ta1_x_xx_zzz_1,  \
                             ta1_x_xxx_xxx_0, \
                             ta1_x_xxx_xxy_0, \
                             ta1_x_xxx_xxz_0, \
                             ta1_x_xxx_xyy_0, \
                             ta1_x_xxx_xyz_0, \
                             ta1_x_xxx_xzz_0, \
                             ta1_x_xxx_yyy_0, \
                             ta1_x_xxx_yyz_0, \
                             ta1_x_xxx_yzz_0, \
                             ta1_x_xxx_zzz_0, \
                             ta_xx_xxx_1,     \
                             ta_xx_xxy_1,     \
                             ta_xx_xxz_1,     \
                             ta_xx_xyy_1,     \
                             ta_xx_xyz_1,     \
                             ta_xx_xzz_1,     \
                             ta_xx_yyy_1,     \
                             ta_xx_yyz_1,     \
                             ta_xx_yzz_1,     \
                             ta_xx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_xxx_0[i] = 2.0 * ta1_x_x_xxx_0[i] * fe_0 - 2.0 * ta1_x_x_xxx_1[i] * fe_0 + 3.0 * ta1_x_xx_xx_0[i] * fe_0 -
                             3.0 * ta1_x_xx_xx_1[i] * fe_0 + ta_xx_xxx_1[i] + ta1_x_xx_xxx_0[i] * pa_x[i] - ta1_x_xx_xxx_1[i] * pc_x[i];

        ta1_x_xxx_xxy_0[i] = 2.0 * ta1_x_x_xxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxy_1[i] * fe_0 + 2.0 * ta1_x_xx_xy_0[i] * fe_0 -
                             2.0 * ta1_x_xx_xy_1[i] * fe_0 + ta_xx_xxy_1[i] + ta1_x_xx_xxy_0[i] * pa_x[i] - ta1_x_xx_xxy_1[i] * pc_x[i];

        ta1_x_xxx_xxz_0[i] = 2.0 * ta1_x_x_xxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxz_1[i] * fe_0 + 2.0 * ta1_x_xx_xz_0[i] * fe_0 -
                             2.0 * ta1_x_xx_xz_1[i] * fe_0 + ta_xx_xxz_1[i] + ta1_x_xx_xxz_0[i] * pa_x[i] - ta1_x_xx_xxz_1[i] * pc_x[i];

        ta1_x_xxx_xyy_0[i] = 2.0 * ta1_x_x_xyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyy_1[i] * fe_0 + ta1_x_xx_yy_0[i] * fe_0 - ta1_x_xx_yy_1[i] * fe_0 +
                             ta_xx_xyy_1[i] + ta1_x_xx_xyy_0[i] * pa_x[i] - ta1_x_xx_xyy_1[i] * pc_x[i];

        ta1_x_xxx_xyz_0[i] = 2.0 * ta1_x_x_xyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyz_1[i] * fe_0 + ta1_x_xx_yz_0[i] * fe_0 - ta1_x_xx_yz_1[i] * fe_0 +
                             ta_xx_xyz_1[i] + ta1_x_xx_xyz_0[i] * pa_x[i] - ta1_x_xx_xyz_1[i] * pc_x[i];

        ta1_x_xxx_xzz_0[i] = 2.0 * ta1_x_x_xzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzz_1[i] * fe_0 + ta1_x_xx_zz_0[i] * fe_0 - ta1_x_xx_zz_1[i] * fe_0 +
                             ta_xx_xzz_1[i] + ta1_x_xx_xzz_0[i] * pa_x[i] - ta1_x_xx_xzz_1[i] * pc_x[i];

        ta1_x_xxx_yyy_0[i] = 2.0 * ta1_x_x_yyy_0[i] * fe_0 - 2.0 * ta1_x_x_yyy_1[i] * fe_0 + ta_xx_yyy_1[i] + ta1_x_xx_yyy_0[i] * pa_x[i] -
                             ta1_x_xx_yyy_1[i] * pc_x[i];

        ta1_x_xxx_yyz_0[i] = 2.0 * ta1_x_x_yyz_0[i] * fe_0 - 2.0 * ta1_x_x_yyz_1[i] * fe_0 + ta_xx_yyz_1[i] + ta1_x_xx_yyz_0[i] * pa_x[i] -
                             ta1_x_xx_yyz_1[i] * pc_x[i];

        ta1_x_xxx_yzz_0[i] = 2.0 * ta1_x_x_yzz_0[i] * fe_0 - 2.0 * ta1_x_x_yzz_1[i] * fe_0 + ta_xx_yzz_1[i] + ta1_x_xx_yzz_0[i] * pa_x[i] -
                             ta1_x_xx_yzz_1[i] * pc_x[i];

        ta1_x_xxx_zzz_0[i] = 2.0 * ta1_x_x_zzz_0[i] * fe_0 - 2.0 * ta1_x_x_zzz_1[i] * fe_0 + ta_xx_zzz_1[i] + ta1_x_xx_zzz_0[i] * pa_x[i] -
                             ta1_x_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto ta1_x_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 10);

    auto ta1_x_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 11);

    auto ta1_x_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 12);

    auto ta1_x_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 13);

    auto ta1_x_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 14);

    auto ta1_x_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 15);

    auto ta1_x_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 16);

    auto ta1_x_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 17);

    auto ta1_x_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 18);

    auto ta1_x_xxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 19);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xxx_0,  \
                             ta1_x_xx_xxx_1,  \
                             ta1_x_xx_xxy_0,  \
                             ta1_x_xx_xxy_1,  \
                             ta1_x_xx_xxz_0,  \
                             ta1_x_xx_xxz_1,  \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xyy_0,  \
                             ta1_x_xx_xyy_1,  \
                             ta1_x_xx_xyz_0,  \
                             ta1_x_xx_xyz_1,  \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_xzz_0,  \
                             ta1_x_xx_xzz_1,  \
                             ta1_x_xx_yy_0,   \
                             ta1_x_xx_yy_1,   \
                             ta1_x_xx_yyy_0,  \
                             ta1_x_xx_yyy_1,  \
                             ta1_x_xx_yyz_0,  \
                             ta1_x_xx_yyz_1,  \
                             ta1_x_xx_yz_0,   \
                             ta1_x_xx_yz_1,   \
                             ta1_x_xx_yzz_0,  \
                             ta1_x_xx_yzz_1,  \
                             ta1_x_xx_zz_0,   \
                             ta1_x_xx_zz_1,   \
                             ta1_x_xx_zzz_0,  \
                             ta1_x_xx_zzz_1,  \
                             ta1_x_xxy_xxx_0, \
                             ta1_x_xxy_xxy_0, \
                             ta1_x_xxy_xxz_0, \
                             ta1_x_xxy_xyy_0, \
                             ta1_x_xxy_xyz_0, \
                             ta1_x_xxy_xzz_0, \
                             ta1_x_xxy_yyy_0, \
                             ta1_x_xxy_yyz_0, \
                             ta1_x_xxy_yzz_0, \
                             ta1_x_xxy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_xxx_0[i] = ta1_x_xx_xxx_0[i] * pa_y[i] - ta1_x_xx_xxx_1[i] * pc_y[i];

        ta1_x_xxy_xxy_0[i] = ta1_x_xx_xx_0[i] * fe_0 - ta1_x_xx_xx_1[i] * fe_0 + ta1_x_xx_xxy_0[i] * pa_y[i] - ta1_x_xx_xxy_1[i] * pc_y[i];

        ta1_x_xxy_xxz_0[i] = ta1_x_xx_xxz_0[i] * pa_y[i] - ta1_x_xx_xxz_1[i] * pc_y[i];

        ta1_x_xxy_xyy_0[i] =
            2.0 * ta1_x_xx_xy_0[i] * fe_0 - 2.0 * ta1_x_xx_xy_1[i] * fe_0 + ta1_x_xx_xyy_0[i] * pa_y[i] - ta1_x_xx_xyy_1[i] * pc_y[i];

        ta1_x_xxy_xyz_0[i] = ta1_x_xx_xz_0[i] * fe_0 - ta1_x_xx_xz_1[i] * fe_0 + ta1_x_xx_xyz_0[i] * pa_y[i] - ta1_x_xx_xyz_1[i] * pc_y[i];

        ta1_x_xxy_xzz_0[i] = ta1_x_xx_xzz_0[i] * pa_y[i] - ta1_x_xx_xzz_1[i] * pc_y[i];

        ta1_x_xxy_yyy_0[i] =
            3.0 * ta1_x_xx_yy_0[i] * fe_0 - 3.0 * ta1_x_xx_yy_1[i] * fe_0 + ta1_x_xx_yyy_0[i] * pa_y[i] - ta1_x_xx_yyy_1[i] * pc_y[i];

        ta1_x_xxy_yyz_0[i] =
            2.0 * ta1_x_xx_yz_0[i] * fe_0 - 2.0 * ta1_x_xx_yz_1[i] * fe_0 + ta1_x_xx_yyz_0[i] * pa_y[i] - ta1_x_xx_yyz_1[i] * pc_y[i];

        ta1_x_xxy_yzz_0[i] = ta1_x_xx_zz_0[i] * fe_0 - ta1_x_xx_zz_1[i] * fe_0 + ta1_x_xx_yzz_0[i] * pa_y[i] - ta1_x_xx_yzz_1[i] * pc_y[i];

        ta1_x_xxy_zzz_0[i] = ta1_x_xx_zzz_0[i] * pa_y[i] - ta1_x_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto ta1_x_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 20);

    auto ta1_x_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 21);

    auto ta1_x_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 22);

    auto ta1_x_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 23);

    auto ta1_x_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 24);

    auto ta1_x_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 25);

    auto ta1_x_xxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 26);

    auto ta1_x_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 27);

    auto ta1_x_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 28);

    auto ta1_x_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 29);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_xx_xx_0,   \
                             ta1_x_xx_xx_1,   \
                             ta1_x_xx_xxx_0,  \
                             ta1_x_xx_xxx_1,  \
                             ta1_x_xx_xxy_0,  \
                             ta1_x_xx_xxy_1,  \
                             ta1_x_xx_xxz_0,  \
                             ta1_x_xx_xxz_1,  \
                             ta1_x_xx_xy_0,   \
                             ta1_x_xx_xy_1,   \
                             ta1_x_xx_xyy_0,  \
                             ta1_x_xx_xyy_1,  \
                             ta1_x_xx_xyz_0,  \
                             ta1_x_xx_xyz_1,  \
                             ta1_x_xx_xz_0,   \
                             ta1_x_xx_xz_1,   \
                             ta1_x_xx_xzz_0,  \
                             ta1_x_xx_xzz_1,  \
                             ta1_x_xx_yy_0,   \
                             ta1_x_xx_yy_1,   \
                             ta1_x_xx_yyy_0,  \
                             ta1_x_xx_yyy_1,  \
                             ta1_x_xx_yyz_0,  \
                             ta1_x_xx_yyz_1,  \
                             ta1_x_xx_yz_0,   \
                             ta1_x_xx_yz_1,   \
                             ta1_x_xx_yzz_0,  \
                             ta1_x_xx_yzz_1,  \
                             ta1_x_xx_zz_0,   \
                             ta1_x_xx_zz_1,   \
                             ta1_x_xx_zzz_0,  \
                             ta1_x_xx_zzz_1,  \
                             ta1_x_xxz_xxx_0, \
                             ta1_x_xxz_xxy_0, \
                             ta1_x_xxz_xxz_0, \
                             ta1_x_xxz_xyy_0, \
                             ta1_x_xxz_xyz_0, \
                             ta1_x_xxz_xzz_0, \
                             ta1_x_xxz_yyy_0, \
                             ta1_x_xxz_yyz_0, \
                             ta1_x_xxz_yzz_0, \
                             ta1_x_xxz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_xxx_0[i] = ta1_x_xx_xxx_0[i] * pa_z[i] - ta1_x_xx_xxx_1[i] * pc_z[i];

        ta1_x_xxz_xxy_0[i] = ta1_x_xx_xxy_0[i] * pa_z[i] - ta1_x_xx_xxy_1[i] * pc_z[i];

        ta1_x_xxz_xxz_0[i] = ta1_x_xx_xx_0[i] * fe_0 - ta1_x_xx_xx_1[i] * fe_0 + ta1_x_xx_xxz_0[i] * pa_z[i] - ta1_x_xx_xxz_1[i] * pc_z[i];

        ta1_x_xxz_xyy_0[i] = ta1_x_xx_xyy_0[i] * pa_z[i] - ta1_x_xx_xyy_1[i] * pc_z[i];

        ta1_x_xxz_xyz_0[i] = ta1_x_xx_xy_0[i] * fe_0 - ta1_x_xx_xy_1[i] * fe_0 + ta1_x_xx_xyz_0[i] * pa_z[i] - ta1_x_xx_xyz_1[i] * pc_z[i];

        ta1_x_xxz_xzz_0[i] =
            2.0 * ta1_x_xx_xz_0[i] * fe_0 - 2.0 * ta1_x_xx_xz_1[i] * fe_0 + ta1_x_xx_xzz_0[i] * pa_z[i] - ta1_x_xx_xzz_1[i] * pc_z[i];

        ta1_x_xxz_yyy_0[i] = ta1_x_xx_yyy_0[i] * pa_z[i] - ta1_x_xx_yyy_1[i] * pc_z[i];

        ta1_x_xxz_yyz_0[i] = ta1_x_xx_yy_0[i] * fe_0 - ta1_x_xx_yy_1[i] * fe_0 + ta1_x_xx_yyz_0[i] * pa_z[i] - ta1_x_xx_yyz_1[i] * pc_z[i];

        ta1_x_xxz_yzz_0[i] =
            2.0 * ta1_x_xx_yz_0[i] * fe_0 - 2.0 * ta1_x_xx_yz_1[i] * fe_0 + ta1_x_xx_yzz_0[i] * pa_z[i] - ta1_x_xx_yzz_1[i] * pc_z[i];

        ta1_x_xxz_zzz_0[i] =
            3.0 * ta1_x_xx_zz_0[i] * fe_0 - 3.0 * ta1_x_xx_zz_1[i] * fe_0 + ta1_x_xx_zzz_0[i] * pa_z[i] - ta1_x_xx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto ta1_x_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 30);

    auto ta1_x_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 31);

    auto ta1_x_xyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 32);

    auto ta1_x_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 33);

    auto ta1_x_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 34);

    auto ta1_x_xyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 35);

    auto ta1_x_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 36);

    auto ta1_x_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 37);

    auto ta1_x_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 38);

    auto ta1_x_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 39);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxz_0,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xzz_0,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_xy_xxx_0,  \
                             ta1_x_xy_xxx_1,  \
                             ta1_x_xy_xxz_0,  \
                             ta1_x_xy_xxz_1,  \
                             ta1_x_xy_xzz_0,  \
                             ta1_x_xy_xzz_1,  \
                             ta1_x_xyy_xxx_0, \
                             ta1_x_xyy_xxy_0, \
                             ta1_x_xyy_xxz_0, \
                             ta1_x_xyy_xyy_0, \
                             ta1_x_xyy_xyz_0, \
                             ta1_x_xyy_xzz_0, \
                             ta1_x_xyy_yyy_0, \
                             ta1_x_xyy_yyz_0, \
                             ta1_x_xyy_yzz_0, \
                             ta1_x_xyy_zzz_0, \
                             ta1_x_yy_xxy_0,  \
                             ta1_x_yy_xxy_1,  \
                             ta1_x_yy_xy_0,   \
                             ta1_x_yy_xy_1,   \
                             ta1_x_yy_xyy_0,  \
                             ta1_x_yy_xyy_1,  \
                             ta1_x_yy_xyz_0,  \
                             ta1_x_yy_xyz_1,  \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yy_yyy_0,  \
                             ta1_x_yy_yyy_1,  \
                             ta1_x_yy_yyz_0,  \
                             ta1_x_yy_yyz_1,  \
                             ta1_x_yy_yz_0,   \
                             ta1_x_yy_yz_1,   \
                             ta1_x_yy_yzz_0,  \
                             ta1_x_yy_yzz_1,  \
                             ta1_x_yy_zzz_0,  \
                             ta1_x_yy_zzz_1,  \
                             ta_yy_xxy_1,     \
                             ta_yy_xyy_1,     \
                             ta_yy_xyz_1,     \
                             ta_yy_yyy_1,     \
                             ta_yy_yyz_1,     \
                             ta_yy_yzz_1,     \
                             ta_yy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_xxx_0[i] = ta1_x_x_xxx_0[i] * fe_0 - ta1_x_x_xxx_1[i] * fe_0 + ta1_x_xy_xxx_0[i] * pa_y[i] - ta1_x_xy_xxx_1[i] * pc_y[i];

        ta1_x_xyy_xxy_0[i] = 2.0 * ta1_x_yy_xy_0[i] * fe_0 - 2.0 * ta1_x_yy_xy_1[i] * fe_0 + ta_yy_xxy_1[i] + ta1_x_yy_xxy_0[i] * pa_x[i] -
                             ta1_x_yy_xxy_1[i] * pc_x[i];

        ta1_x_xyy_xxz_0[i] = ta1_x_x_xxz_0[i] * fe_0 - ta1_x_x_xxz_1[i] * fe_0 + ta1_x_xy_xxz_0[i] * pa_y[i] - ta1_x_xy_xxz_1[i] * pc_y[i];

        ta1_x_xyy_xyy_0[i] =
            ta1_x_yy_yy_0[i] * fe_0 - ta1_x_yy_yy_1[i] * fe_0 + ta_yy_xyy_1[i] + ta1_x_yy_xyy_0[i] * pa_x[i] - ta1_x_yy_xyy_1[i] * pc_x[i];

        ta1_x_xyy_xyz_0[i] =
            ta1_x_yy_yz_0[i] * fe_0 - ta1_x_yy_yz_1[i] * fe_0 + ta_yy_xyz_1[i] + ta1_x_yy_xyz_0[i] * pa_x[i] - ta1_x_yy_xyz_1[i] * pc_x[i];

        ta1_x_xyy_xzz_0[i] = ta1_x_x_xzz_0[i] * fe_0 - ta1_x_x_xzz_1[i] * fe_0 + ta1_x_xy_xzz_0[i] * pa_y[i] - ta1_x_xy_xzz_1[i] * pc_y[i];

        ta1_x_xyy_yyy_0[i] = ta_yy_yyy_1[i] + ta1_x_yy_yyy_0[i] * pa_x[i] - ta1_x_yy_yyy_1[i] * pc_x[i];

        ta1_x_xyy_yyz_0[i] = ta_yy_yyz_1[i] + ta1_x_yy_yyz_0[i] * pa_x[i] - ta1_x_yy_yyz_1[i] * pc_x[i];

        ta1_x_xyy_yzz_0[i] = ta_yy_yzz_1[i] + ta1_x_yy_yzz_0[i] * pa_x[i] - ta1_x_yy_yzz_1[i] * pc_x[i];

        ta1_x_xyy_zzz_0[i] = ta_yy_zzz_1[i] + ta1_x_yy_zzz_0[i] * pa_x[i] - ta1_x_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto ta1_x_xyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 40);

    auto ta1_x_xyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 41);

    auto ta1_x_xyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 42);

    auto ta1_x_xyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 43);

    auto ta1_x_xyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 44);

    auto ta1_x_xyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 45);

    auto ta1_x_xyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 46);

    auto ta1_x_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 47);

    auto ta1_x_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 48);

    auto ta1_x_xyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 49);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xy_xxy_0,  \
                             ta1_x_xy_xxy_1,  \
                             ta1_x_xy_xyy_0,  \
                             ta1_x_xy_xyy_1,  \
                             ta1_x_xy_yyy_0,  \
                             ta1_x_xy_yyy_1,  \
                             ta1_x_xyz_xxx_0, \
                             ta1_x_xyz_xxy_0, \
                             ta1_x_xyz_xxz_0, \
                             ta1_x_xyz_xyy_0, \
                             ta1_x_xyz_xyz_0, \
                             ta1_x_xyz_xzz_0, \
                             ta1_x_xyz_yyy_0, \
                             ta1_x_xyz_yyz_0, \
                             ta1_x_xyz_yzz_0, \
                             ta1_x_xyz_zzz_0, \
                             ta1_x_xz_xxx_0,  \
                             ta1_x_xz_xxx_1,  \
                             ta1_x_xz_xxz_0,  \
                             ta1_x_xz_xxz_1,  \
                             ta1_x_xz_xyz_0,  \
                             ta1_x_xz_xyz_1,  \
                             ta1_x_xz_xz_0,   \
                             ta1_x_xz_xz_1,   \
                             ta1_x_xz_xzz_0,  \
                             ta1_x_xz_xzz_1,  \
                             ta1_x_xz_zzz_0,  \
                             ta1_x_xz_zzz_1,  \
                             ta1_x_yz_yyz_0,  \
                             ta1_x_yz_yyz_1,  \
                             ta1_x_yz_yzz_0,  \
                             ta1_x_yz_yzz_1,  \
                             ta_yz_yyz_1,     \
                             ta_yz_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyz_xxx_0[i] = ta1_x_xz_xxx_0[i] * pa_y[i] - ta1_x_xz_xxx_1[i] * pc_y[i];

        ta1_x_xyz_xxy_0[i] = ta1_x_xy_xxy_0[i] * pa_z[i] - ta1_x_xy_xxy_1[i] * pc_z[i];

        ta1_x_xyz_xxz_0[i] = ta1_x_xz_xxz_0[i] * pa_y[i] - ta1_x_xz_xxz_1[i] * pc_y[i];

        ta1_x_xyz_xyy_0[i] = ta1_x_xy_xyy_0[i] * pa_z[i] - ta1_x_xy_xyy_1[i] * pc_z[i];

        ta1_x_xyz_xyz_0[i] = ta1_x_xz_xz_0[i] * fe_0 - ta1_x_xz_xz_1[i] * fe_0 + ta1_x_xz_xyz_0[i] * pa_y[i] - ta1_x_xz_xyz_1[i] * pc_y[i];

        ta1_x_xyz_xzz_0[i] = ta1_x_xz_xzz_0[i] * pa_y[i] - ta1_x_xz_xzz_1[i] * pc_y[i];

        ta1_x_xyz_yyy_0[i] = ta1_x_xy_yyy_0[i] * pa_z[i] - ta1_x_xy_yyy_1[i] * pc_z[i];

        ta1_x_xyz_yyz_0[i] = ta_yz_yyz_1[i] + ta1_x_yz_yyz_0[i] * pa_x[i] - ta1_x_yz_yyz_1[i] * pc_x[i];

        ta1_x_xyz_yzz_0[i] = ta_yz_yzz_1[i] + ta1_x_yz_yzz_0[i] * pa_x[i] - ta1_x_yz_yzz_1[i] * pc_x[i];

        ta1_x_xyz_zzz_0[i] = ta1_x_xz_zzz_0[i] * pa_y[i] - ta1_x_xz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto ta1_x_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 50);

    auto ta1_x_xzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 51);

    auto ta1_x_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 52);

    auto ta1_x_xzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 53);

    auto ta1_x_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 54);

    auto ta1_x_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 55);

    auto ta1_x_xzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 56);

    auto ta1_x_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 57);

    auto ta1_x_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 58);

    auto ta1_x_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 59);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxy_0,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xyy_0,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_xz_xxx_0,  \
                             ta1_x_xz_xxx_1,  \
                             ta1_x_xz_xxy_0,  \
                             ta1_x_xz_xxy_1,  \
                             ta1_x_xz_xyy_0,  \
                             ta1_x_xz_xyy_1,  \
                             ta1_x_xzz_xxx_0, \
                             ta1_x_xzz_xxy_0, \
                             ta1_x_xzz_xxz_0, \
                             ta1_x_xzz_xyy_0, \
                             ta1_x_xzz_xyz_0, \
                             ta1_x_xzz_xzz_0, \
                             ta1_x_xzz_yyy_0, \
                             ta1_x_xzz_yyz_0, \
                             ta1_x_xzz_yzz_0, \
                             ta1_x_xzz_zzz_0, \
                             ta1_x_zz_xxz_0,  \
                             ta1_x_zz_xxz_1,  \
                             ta1_x_zz_xyz_0,  \
                             ta1_x_zz_xyz_1,  \
                             ta1_x_zz_xz_0,   \
                             ta1_x_zz_xz_1,   \
                             ta1_x_zz_xzz_0,  \
                             ta1_x_zz_xzz_1,  \
                             ta1_x_zz_yyy_0,  \
                             ta1_x_zz_yyy_1,  \
                             ta1_x_zz_yyz_0,  \
                             ta1_x_zz_yyz_1,  \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_yzz_0,  \
                             ta1_x_zz_yzz_1,  \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             ta1_x_zz_zzz_0,  \
                             ta1_x_zz_zzz_1,  \
                             ta_zz_xxz_1,     \
                             ta_zz_xyz_1,     \
                             ta_zz_xzz_1,     \
                             ta_zz_yyy_1,     \
                             ta_zz_yyz_1,     \
                             ta_zz_yzz_1,     \
                             ta_zz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_xxx_0[i] = ta1_x_x_xxx_0[i] * fe_0 - ta1_x_x_xxx_1[i] * fe_0 + ta1_x_xz_xxx_0[i] * pa_z[i] - ta1_x_xz_xxx_1[i] * pc_z[i];

        ta1_x_xzz_xxy_0[i] = ta1_x_x_xxy_0[i] * fe_0 - ta1_x_x_xxy_1[i] * fe_0 + ta1_x_xz_xxy_0[i] * pa_z[i] - ta1_x_xz_xxy_1[i] * pc_z[i];

        ta1_x_xzz_xxz_0[i] = 2.0 * ta1_x_zz_xz_0[i] * fe_0 - 2.0 * ta1_x_zz_xz_1[i] * fe_0 + ta_zz_xxz_1[i] + ta1_x_zz_xxz_0[i] * pa_x[i] -
                             ta1_x_zz_xxz_1[i] * pc_x[i];

        ta1_x_xzz_xyy_0[i] = ta1_x_x_xyy_0[i] * fe_0 - ta1_x_x_xyy_1[i] * fe_0 + ta1_x_xz_xyy_0[i] * pa_z[i] - ta1_x_xz_xyy_1[i] * pc_z[i];

        ta1_x_xzz_xyz_0[i] =
            ta1_x_zz_yz_0[i] * fe_0 - ta1_x_zz_yz_1[i] * fe_0 + ta_zz_xyz_1[i] + ta1_x_zz_xyz_0[i] * pa_x[i] - ta1_x_zz_xyz_1[i] * pc_x[i];

        ta1_x_xzz_xzz_0[i] =
            ta1_x_zz_zz_0[i] * fe_0 - ta1_x_zz_zz_1[i] * fe_0 + ta_zz_xzz_1[i] + ta1_x_zz_xzz_0[i] * pa_x[i] - ta1_x_zz_xzz_1[i] * pc_x[i];

        ta1_x_xzz_yyy_0[i] = ta_zz_yyy_1[i] + ta1_x_zz_yyy_0[i] * pa_x[i] - ta1_x_zz_yyy_1[i] * pc_x[i];

        ta1_x_xzz_yyz_0[i] = ta_zz_yyz_1[i] + ta1_x_zz_yyz_0[i] * pa_x[i] - ta1_x_zz_yyz_1[i] * pc_x[i];

        ta1_x_xzz_yzz_0[i] = ta_zz_yzz_1[i] + ta1_x_zz_yzz_0[i] * pa_x[i] - ta1_x_zz_yzz_1[i] * pc_x[i];

        ta1_x_xzz_zzz_0[i] = ta_zz_zzz_1[i] + ta1_x_zz_zzz_0[i] * pa_x[i] - ta1_x_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto ta1_x_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 60);

    auto ta1_x_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 61);

    auto ta1_x_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 62);

    auto ta1_x_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 63);

    auto ta1_x_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 64);

    auto ta1_x_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 65);

    auto ta1_x_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 66);

    auto ta1_x_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 67);

    auto ta1_x_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 68);

    auto ta1_x_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 69);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_y_xxx_0,   \
                             ta1_x_y_xxx_1,   \
                             ta1_x_y_xxy_0,   \
                             ta1_x_y_xxy_1,   \
                             ta1_x_y_xxz_0,   \
                             ta1_x_y_xxz_1,   \
                             ta1_x_y_xyy_0,   \
                             ta1_x_y_xyy_1,   \
                             ta1_x_y_xyz_0,   \
                             ta1_x_y_xyz_1,   \
                             ta1_x_y_xzz_0,   \
                             ta1_x_y_xzz_1,   \
                             ta1_x_y_yyy_0,   \
                             ta1_x_y_yyy_1,   \
                             ta1_x_y_yyz_0,   \
                             ta1_x_y_yyz_1,   \
                             ta1_x_y_yzz_0,   \
                             ta1_x_y_yzz_1,   \
                             ta1_x_y_zzz_0,   \
                             ta1_x_y_zzz_1,   \
                             ta1_x_yy_xx_0,   \
                             ta1_x_yy_xx_1,   \
                             ta1_x_yy_xxx_0,  \
                             ta1_x_yy_xxx_1,  \
                             ta1_x_yy_xxy_0,  \
                             ta1_x_yy_xxy_1,  \
                             ta1_x_yy_xxz_0,  \
                             ta1_x_yy_xxz_1,  \
                             ta1_x_yy_xy_0,   \
                             ta1_x_yy_xy_1,   \
                             ta1_x_yy_xyy_0,  \
                             ta1_x_yy_xyy_1,  \
                             ta1_x_yy_xyz_0,  \
                             ta1_x_yy_xyz_1,  \
                             ta1_x_yy_xz_0,   \
                             ta1_x_yy_xz_1,   \
                             ta1_x_yy_xzz_0,  \
                             ta1_x_yy_xzz_1,  \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yy_yyy_0,  \
                             ta1_x_yy_yyy_1,  \
                             ta1_x_yy_yyz_0,  \
                             ta1_x_yy_yyz_1,  \
                             ta1_x_yy_yz_0,   \
                             ta1_x_yy_yz_1,   \
                             ta1_x_yy_yzz_0,  \
                             ta1_x_yy_yzz_1,  \
                             ta1_x_yy_zz_0,   \
                             ta1_x_yy_zz_1,   \
                             ta1_x_yy_zzz_0,  \
                             ta1_x_yy_zzz_1,  \
                             ta1_x_yyy_xxx_0, \
                             ta1_x_yyy_xxy_0, \
                             ta1_x_yyy_xxz_0, \
                             ta1_x_yyy_xyy_0, \
                             ta1_x_yyy_xyz_0, \
                             ta1_x_yyy_xzz_0, \
                             ta1_x_yyy_yyy_0, \
                             ta1_x_yyy_yyz_0, \
                             ta1_x_yyy_yzz_0, \
                             ta1_x_yyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_xxx_0[i] =
            2.0 * ta1_x_y_xxx_0[i] * fe_0 - 2.0 * ta1_x_y_xxx_1[i] * fe_0 + ta1_x_yy_xxx_0[i] * pa_y[i] - ta1_x_yy_xxx_1[i] * pc_y[i];

        ta1_x_yyy_xxy_0[i] = 2.0 * ta1_x_y_xxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxy_1[i] * fe_0 + ta1_x_yy_xx_0[i] * fe_0 - ta1_x_yy_xx_1[i] * fe_0 +
                             ta1_x_yy_xxy_0[i] * pa_y[i] - ta1_x_yy_xxy_1[i] * pc_y[i];

        ta1_x_yyy_xxz_0[i] =
            2.0 * ta1_x_y_xxz_0[i] * fe_0 - 2.0 * ta1_x_y_xxz_1[i] * fe_0 + ta1_x_yy_xxz_0[i] * pa_y[i] - ta1_x_yy_xxz_1[i] * pc_y[i];

        ta1_x_yyy_xyy_0[i] = 2.0 * ta1_x_y_xyy_0[i] * fe_0 - 2.0 * ta1_x_y_xyy_1[i] * fe_0 + 2.0 * ta1_x_yy_xy_0[i] * fe_0 -
                             2.0 * ta1_x_yy_xy_1[i] * fe_0 + ta1_x_yy_xyy_0[i] * pa_y[i] - ta1_x_yy_xyy_1[i] * pc_y[i];

        ta1_x_yyy_xyz_0[i] = 2.0 * ta1_x_y_xyz_0[i] * fe_0 - 2.0 * ta1_x_y_xyz_1[i] * fe_0 + ta1_x_yy_xz_0[i] * fe_0 - ta1_x_yy_xz_1[i] * fe_0 +
                             ta1_x_yy_xyz_0[i] * pa_y[i] - ta1_x_yy_xyz_1[i] * pc_y[i];

        ta1_x_yyy_xzz_0[i] =
            2.0 * ta1_x_y_xzz_0[i] * fe_0 - 2.0 * ta1_x_y_xzz_1[i] * fe_0 + ta1_x_yy_xzz_0[i] * pa_y[i] - ta1_x_yy_xzz_1[i] * pc_y[i];

        ta1_x_yyy_yyy_0[i] = 2.0 * ta1_x_y_yyy_0[i] * fe_0 - 2.0 * ta1_x_y_yyy_1[i] * fe_0 + 3.0 * ta1_x_yy_yy_0[i] * fe_0 -
                             3.0 * ta1_x_yy_yy_1[i] * fe_0 + ta1_x_yy_yyy_0[i] * pa_y[i] - ta1_x_yy_yyy_1[i] * pc_y[i];

        ta1_x_yyy_yyz_0[i] = 2.0 * ta1_x_y_yyz_0[i] * fe_0 - 2.0 * ta1_x_y_yyz_1[i] * fe_0 + 2.0 * ta1_x_yy_yz_0[i] * fe_0 -
                             2.0 * ta1_x_yy_yz_1[i] * fe_0 + ta1_x_yy_yyz_0[i] * pa_y[i] - ta1_x_yy_yyz_1[i] * pc_y[i];

        ta1_x_yyy_yzz_0[i] = 2.0 * ta1_x_y_yzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzz_1[i] * fe_0 + ta1_x_yy_zz_0[i] * fe_0 - ta1_x_yy_zz_1[i] * fe_0 +
                             ta1_x_yy_yzz_0[i] * pa_y[i] - ta1_x_yy_yzz_1[i] * pc_y[i];

        ta1_x_yyy_zzz_0[i] =
            2.0 * ta1_x_y_zzz_0[i] * fe_0 - 2.0 * ta1_x_y_zzz_1[i] * fe_0 + ta1_x_yy_zzz_0[i] * pa_y[i] - ta1_x_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto ta1_x_yyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 70);

    auto ta1_x_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 71);

    auto ta1_x_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 72);

    auto ta1_x_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 73);

    auto ta1_x_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 74);

    auto ta1_x_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 75);

    auto ta1_x_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 76);

    auto ta1_x_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 77);

    auto ta1_x_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 78);

    auto ta1_x_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 79);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yy_xxx_0,  \
                             ta1_x_yy_xxx_1,  \
                             ta1_x_yy_xxy_0,  \
                             ta1_x_yy_xxy_1,  \
                             ta1_x_yy_xy_0,   \
                             ta1_x_yy_xy_1,   \
                             ta1_x_yy_xyy_0,  \
                             ta1_x_yy_xyy_1,  \
                             ta1_x_yy_xyz_0,  \
                             ta1_x_yy_xyz_1,  \
                             ta1_x_yy_yy_0,   \
                             ta1_x_yy_yy_1,   \
                             ta1_x_yy_yyy_0,  \
                             ta1_x_yy_yyy_1,  \
                             ta1_x_yy_yyz_0,  \
                             ta1_x_yy_yyz_1,  \
                             ta1_x_yy_yz_0,   \
                             ta1_x_yy_yz_1,   \
                             ta1_x_yy_yzz_0,  \
                             ta1_x_yy_yzz_1,  \
                             ta1_x_yyz_xxx_0, \
                             ta1_x_yyz_xxy_0, \
                             ta1_x_yyz_xxz_0, \
                             ta1_x_yyz_xyy_0, \
                             ta1_x_yyz_xyz_0, \
                             ta1_x_yyz_xzz_0, \
                             ta1_x_yyz_yyy_0, \
                             ta1_x_yyz_yyz_0, \
                             ta1_x_yyz_yzz_0, \
                             ta1_x_yyz_zzz_0, \
                             ta1_x_yz_xxz_0,  \
                             ta1_x_yz_xxz_1,  \
                             ta1_x_yz_xzz_0,  \
                             ta1_x_yz_xzz_1,  \
                             ta1_x_yz_zzz_0,  \
                             ta1_x_yz_zzz_1,  \
                             ta1_x_z_xxz_0,   \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xzz_0,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_zzz_0,   \
                             ta1_x_z_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_xxx_0[i] = ta1_x_yy_xxx_0[i] * pa_z[i] - ta1_x_yy_xxx_1[i] * pc_z[i];

        ta1_x_yyz_xxy_0[i] = ta1_x_yy_xxy_0[i] * pa_z[i] - ta1_x_yy_xxy_1[i] * pc_z[i];

        ta1_x_yyz_xxz_0[i] = ta1_x_z_xxz_0[i] * fe_0 - ta1_x_z_xxz_1[i] * fe_0 + ta1_x_yz_xxz_0[i] * pa_y[i] - ta1_x_yz_xxz_1[i] * pc_y[i];

        ta1_x_yyz_xyy_0[i] = ta1_x_yy_xyy_0[i] * pa_z[i] - ta1_x_yy_xyy_1[i] * pc_z[i];

        ta1_x_yyz_xyz_0[i] = ta1_x_yy_xy_0[i] * fe_0 - ta1_x_yy_xy_1[i] * fe_0 + ta1_x_yy_xyz_0[i] * pa_z[i] - ta1_x_yy_xyz_1[i] * pc_z[i];

        ta1_x_yyz_xzz_0[i] = ta1_x_z_xzz_0[i] * fe_0 - ta1_x_z_xzz_1[i] * fe_0 + ta1_x_yz_xzz_0[i] * pa_y[i] - ta1_x_yz_xzz_1[i] * pc_y[i];

        ta1_x_yyz_yyy_0[i] = ta1_x_yy_yyy_0[i] * pa_z[i] - ta1_x_yy_yyy_1[i] * pc_z[i];

        ta1_x_yyz_yyz_0[i] = ta1_x_yy_yy_0[i] * fe_0 - ta1_x_yy_yy_1[i] * fe_0 + ta1_x_yy_yyz_0[i] * pa_z[i] - ta1_x_yy_yyz_1[i] * pc_z[i];

        ta1_x_yyz_yzz_0[i] =
            2.0 * ta1_x_yy_yz_0[i] * fe_0 - 2.0 * ta1_x_yy_yz_1[i] * fe_0 + ta1_x_yy_yzz_0[i] * pa_z[i] - ta1_x_yy_yzz_1[i] * pc_z[i];

        ta1_x_yyz_zzz_0[i] = ta1_x_z_zzz_0[i] * fe_0 - ta1_x_z_zzz_1[i] * fe_0 + ta1_x_yz_zzz_0[i] * pa_y[i] - ta1_x_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto ta1_x_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 80);

    auto ta1_x_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 81);

    auto ta1_x_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 82);

    auto ta1_x_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 83);

    auto ta1_x_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 84);

    auto ta1_x_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 85);

    auto ta1_x_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 86);

    auto ta1_x_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 87);

    auto ta1_x_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 88);

    auto ta1_x_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 89);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yzz_xxx_0, \
                             ta1_x_yzz_xxy_0, \
                             ta1_x_yzz_xxz_0, \
                             ta1_x_yzz_xyy_0, \
                             ta1_x_yzz_xyz_0, \
                             ta1_x_yzz_xzz_0, \
                             ta1_x_yzz_yyy_0, \
                             ta1_x_yzz_yyz_0, \
                             ta1_x_yzz_yzz_0, \
                             ta1_x_yzz_zzz_0, \
                             ta1_x_zz_xx_0,   \
                             ta1_x_zz_xx_1,   \
                             ta1_x_zz_xxx_0,  \
                             ta1_x_zz_xxx_1,  \
                             ta1_x_zz_xxy_0,  \
                             ta1_x_zz_xxy_1,  \
                             ta1_x_zz_xxz_0,  \
                             ta1_x_zz_xxz_1,  \
                             ta1_x_zz_xy_0,   \
                             ta1_x_zz_xy_1,   \
                             ta1_x_zz_xyy_0,  \
                             ta1_x_zz_xyy_1,  \
                             ta1_x_zz_xyz_0,  \
                             ta1_x_zz_xyz_1,  \
                             ta1_x_zz_xz_0,   \
                             ta1_x_zz_xz_1,   \
                             ta1_x_zz_xzz_0,  \
                             ta1_x_zz_xzz_1,  \
                             ta1_x_zz_yy_0,   \
                             ta1_x_zz_yy_1,   \
                             ta1_x_zz_yyy_0,  \
                             ta1_x_zz_yyy_1,  \
                             ta1_x_zz_yyz_0,  \
                             ta1_x_zz_yyz_1,  \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_yzz_0,  \
                             ta1_x_zz_yzz_1,  \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             ta1_x_zz_zzz_0,  \
                             ta1_x_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_xxx_0[i] = ta1_x_zz_xxx_0[i] * pa_y[i] - ta1_x_zz_xxx_1[i] * pc_y[i];

        ta1_x_yzz_xxy_0[i] = ta1_x_zz_xx_0[i] * fe_0 - ta1_x_zz_xx_1[i] * fe_0 + ta1_x_zz_xxy_0[i] * pa_y[i] - ta1_x_zz_xxy_1[i] * pc_y[i];

        ta1_x_yzz_xxz_0[i] = ta1_x_zz_xxz_0[i] * pa_y[i] - ta1_x_zz_xxz_1[i] * pc_y[i];

        ta1_x_yzz_xyy_0[i] =
            2.0 * ta1_x_zz_xy_0[i] * fe_0 - 2.0 * ta1_x_zz_xy_1[i] * fe_0 + ta1_x_zz_xyy_0[i] * pa_y[i] - ta1_x_zz_xyy_1[i] * pc_y[i];

        ta1_x_yzz_xyz_0[i] = ta1_x_zz_xz_0[i] * fe_0 - ta1_x_zz_xz_1[i] * fe_0 + ta1_x_zz_xyz_0[i] * pa_y[i] - ta1_x_zz_xyz_1[i] * pc_y[i];

        ta1_x_yzz_xzz_0[i] = ta1_x_zz_xzz_0[i] * pa_y[i] - ta1_x_zz_xzz_1[i] * pc_y[i];

        ta1_x_yzz_yyy_0[i] =
            3.0 * ta1_x_zz_yy_0[i] * fe_0 - 3.0 * ta1_x_zz_yy_1[i] * fe_0 + ta1_x_zz_yyy_0[i] * pa_y[i] - ta1_x_zz_yyy_1[i] * pc_y[i];

        ta1_x_yzz_yyz_0[i] =
            2.0 * ta1_x_zz_yz_0[i] * fe_0 - 2.0 * ta1_x_zz_yz_1[i] * fe_0 + ta1_x_zz_yyz_0[i] * pa_y[i] - ta1_x_zz_yyz_1[i] * pc_y[i];

        ta1_x_yzz_yzz_0[i] = ta1_x_zz_zz_0[i] * fe_0 - ta1_x_zz_zz_1[i] * fe_0 + ta1_x_zz_yzz_0[i] * pa_y[i] - ta1_x_zz_yzz_1[i] * pc_y[i];

        ta1_x_yzz_zzz_0[i] = ta1_x_zz_zzz_0[i] * pa_y[i] - ta1_x_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto ta1_x_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 90);

    auto ta1_x_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 91);

    auto ta1_x_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 92);

    auto ta1_x_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 93);

    auto ta1_x_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 94);

    auto ta1_x_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 95);

    auto ta1_x_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 96);

    auto ta1_x_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 97);

    auto ta1_x_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 98);

    auto ta1_x_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 99);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_z_xxx_0,   \
                             ta1_x_z_xxx_1,   \
                             ta1_x_z_xxy_0,   \
                             ta1_x_z_xxy_1,   \
                             ta1_x_z_xxz_0,   \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xyy_0,   \
                             ta1_x_z_xyy_1,   \
                             ta1_x_z_xyz_0,   \
                             ta1_x_z_xyz_1,   \
                             ta1_x_z_xzz_0,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_yyy_0,   \
                             ta1_x_z_yyy_1,   \
                             ta1_x_z_yyz_0,   \
                             ta1_x_z_yyz_1,   \
                             ta1_x_z_yzz_0,   \
                             ta1_x_z_yzz_1,   \
                             ta1_x_z_zzz_0,   \
                             ta1_x_z_zzz_1,   \
                             ta1_x_zz_xx_0,   \
                             ta1_x_zz_xx_1,   \
                             ta1_x_zz_xxx_0,  \
                             ta1_x_zz_xxx_1,  \
                             ta1_x_zz_xxy_0,  \
                             ta1_x_zz_xxy_1,  \
                             ta1_x_zz_xxz_0,  \
                             ta1_x_zz_xxz_1,  \
                             ta1_x_zz_xy_0,   \
                             ta1_x_zz_xy_1,   \
                             ta1_x_zz_xyy_0,  \
                             ta1_x_zz_xyy_1,  \
                             ta1_x_zz_xyz_0,  \
                             ta1_x_zz_xyz_1,  \
                             ta1_x_zz_xz_0,   \
                             ta1_x_zz_xz_1,   \
                             ta1_x_zz_xzz_0,  \
                             ta1_x_zz_xzz_1,  \
                             ta1_x_zz_yy_0,   \
                             ta1_x_zz_yy_1,   \
                             ta1_x_zz_yyy_0,  \
                             ta1_x_zz_yyy_1,  \
                             ta1_x_zz_yyz_0,  \
                             ta1_x_zz_yyz_1,  \
                             ta1_x_zz_yz_0,   \
                             ta1_x_zz_yz_1,   \
                             ta1_x_zz_yzz_0,  \
                             ta1_x_zz_yzz_1,  \
                             ta1_x_zz_zz_0,   \
                             ta1_x_zz_zz_1,   \
                             ta1_x_zz_zzz_0,  \
                             ta1_x_zz_zzz_1,  \
                             ta1_x_zzz_xxx_0, \
                             ta1_x_zzz_xxy_0, \
                             ta1_x_zzz_xxz_0, \
                             ta1_x_zzz_xyy_0, \
                             ta1_x_zzz_xyz_0, \
                             ta1_x_zzz_xzz_0, \
                             ta1_x_zzz_yyy_0, \
                             ta1_x_zzz_yyz_0, \
                             ta1_x_zzz_yzz_0, \
                             ta1_x_zzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_xxx_0[i] =
            2.0 * ta1_x_z_xxx_0[i] * fe_0 - 2.0 * ta1_x_z_xxx_1[i] * fe_0 + ta1_x_zz_xxx_0[i] * pa_z[i] - ta1_x_zz_xxx_1[i] * pc_z[i];

        ta1_x_zzz_xxy_0[i] =
            2.0 * ta1_x_z_xxy_0[i] * fe_0 - 2.0 * ta1_x_z_xxy_1[i] * fe_0 + ta1_x_zz_xxy_0[i] * pa_z[i] - ta1_x_zz_xxy_1[i] * pc_z[i];

        ta1_x_zzz_xxz_0[i] = 2.0 * ta1_x_z_xxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxz_1[i] * fe_0 + ta1_x_zz_xx_0[i] * fe_0 - ta1_x_zz_xx_1[i] * fe_0 +
                             ta1_x_zz_xxz_0[i] * pa_z[i] - ta1_x_zz_xxz_1[i] * pc_z[i];

        ta1_x_zzz_xyy_0[i] =
            2.0 * ta1_x_z_xyy_0[i] * fe_0 - 2.0 * ta1_x_z_xyy_1[i] * fe_0 + ta1_x_zz_xyy_0[i] * pa_z[i] - ta1_x_zz_xyy_1[i] * pc_z[i];

        ta1_x_zzz_xyz_0[i] = 2.0 * ta1_x_z_xyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyz_1[i] * fe_0 + ta1_x_zz_xy_0[i] * fe_0 - ta1_x_zz_xy_1[i] * fe_0 +
                             ta1_x_zz_xyz_0[i] * pa_z[i] - ta1_x_zz_xyz_1[i] * pc_z[i];

        ta1_x_zzz_xzz_0[i] = 2.0 * ta1_x_z_xzz_0[i] * fe_0 - 2.0 * ta1_x_z_xzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xz_0[i] * fe_0 -
                             2.0 * ta1_x_zz_xz_1[i] * fe_0 + ta1_x_zz_xzz_0[i] * pa_z[i] - ta1_x_zz_xzz_1[i] * pc_z[i];

        ta1_x_zzz_yyy_0[i] =
            2.0 * ta1_x_z_yyy_0[i] * fe_0 - 2.0 * ta1_x_z_yyy_1[i] * fe_0 + ta1_x_zz_yyy_0[i] * pa_z[i] - ta1_x_zz_yyy_1[i] * pc_z[i];

        ta1_x_zzz_yyz_0[i] = 2.0 * ta1_x_z_yyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyz_1[i] * fe_0 + ta1_x_zz_yy_0[i] * fe_0 - ta1_x_zz_yy_1[i] * fe_0 +
                             ta1_x_zz_yyz_0[i] * pa_z[i] - ta1_x_zz_yyz_1[i] * pc_z[i];

        ta1_x_zzz_yzz_0[i] = 2.0 * ta1_x_z_yzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzz_1[i] * fe_0 + 2.0 * ta1_x_zz_yz_0[i] * fe_0 -
                             2.0 * ta1_x_zz_yz_1[i] * fe_0 + ta1_x_zz_yzz_0[i] * pa_z[i] - ta1_x_zz_yzz_1[i] * pc_z[i];

        ta1_x_zzz_zzz_0[i] = 2.0 * ta1_x_z_zzz_0[i] * fe_0 - 2.0 * ta1_x_z_zzz_1[i] * fe_0 + 3.0 * ta1_x_zz_zz_0[i] * fe_0 -
                             3.0 * ta1_x_zz_zz_1[i] * fe_0 + ta1_x_zz_zzz_0[i] * pa_z[i] - ta1_x_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 100-110 components of targeted buffer : FF

    auto ta1_y_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 100);

    auto ta1_y_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 101);

    auto ta1_y_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 102);

    auto ta1_y_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 103);

    auto ta1_y_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 104);

    auto ta1_y_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 105);

    auto ta1_y_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 106);

    auto ta1_y_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 107);

    auto ta1_y_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 108);

    auto ta1_y_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 109);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_x_xxx_0,   \
                             ta1_y_x_xxx_1,   \
                             ta1_y_x_xxy_0,   \
                             ta1_y_x_xxy_1,   \
                             ta1_y_x_xxz_0,   \
                             ta1_y_x_xxz_1,   \
                             ta1_y_x_xyy_0,   \
                             ta1_y_x_xyy_1,   \
                             ta1_y_x_xyz_0,   \
                             ta1_y_x_xyz_1,   \
                             ta1_y_x_xzz_0,   \
                             ta1_y_x_xzz_1,   \
                             ta1_y_x_yyy_0,   \
                             ta1_y_x_yyy_1,   \
                             ta1_y_x_yyz_0,   \
                             ta1_y_x_yyz_1,   \
                             ta1_y_x_yzz_0,   \
                             ta1_y_x_yzz_1,   \
                             ta1_y_x_zzz_0,   \
                             ta1_y_x_zzz_1,   \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xxx_0,  \
                             ta1_y_xx_xxx_1,  \
                             ta1_y_xx_xxy_0,  \
                             ta1_y_xx_xxy_1,  \
                             ta1_y_xx_xxz_0,  \
                             ta1_y_xx_xxz_1,  \
                             ta1_y_xx_xy_0,   \
                             ta1_y_xx_xy_1,   \
                             ta1_y_xx_xyy_0,  \
                             ta1_y_xx_xyy_1,  \
                             ta1_y_xx_xyz_0,  \
                             ta1_y_xx_xyz_1,  \
                             ta1_y_xx_xz_0,   \
                             ta1_y_xx_xz_1,   \
                             ta1_y_xx_xzz_0,  \
                             ta1_y_xx_xzz_1,  \
                             ta1_y_xx_yy_0,   \
                             ta1_y_xx_yy_1,   \
                             ta1_y_xx_yyy_0,  \
                             ta1_y_xx_yyy_1,  \
                             ta1_y_xx_yyz_0,  \
                             ta1_y_xx_yyz_1,  \
                             ta1_y_xx_yz_0,   \
                             ta1_y_xx_yz_1,   \
                             ta1_y_xx_yzz_0,  \
                             ta1_y_xx_yzz_1,  \
                             ta1_y_xx_zz_0,   \
                             ta1_y_xx_zz_1,   \
                             ta1_y_xx_zzz_0,  \
                             ta1_y_xx_zzz_1,  \
                             ta1_y_xxx_xxx_0, \
                             ta1_y_xxx_xxy_0, \
                             ta1_y_xxx_xxz_0, \
                             ta1_y_xxx_xyy_0, \
                             ta1_y_xxx_xyz_0, \
                             ta1_y_xxx_xzz_0, \
                             ta1_y_xxx_yyy_0, \
                             ta1_y_xxx_yyz_0, \
                             ta1_y_xxx_yzz_0, \
                             ta1_y_xxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_xxx_0[i] = 2.0 * ta1_y_x_xxx_0[i] * fe_0 - 2.0 * ta1_y_x_xxx_1[i] * fe_0 + 3.0 * ta1_y_xx_xx_0[i] * fe_0 -
                             3.0 * ta1_y_xx_xx_1[i] * fe_0 + ta1_y_xx_xxx_0[i] * pa_x[i] - ta1_y_xx_xxx_1[i] * pc_x[i];

        ta1_y_xxx_xxy_0[i] = 2.0 * ta1_y_x_xxy_0[i] * fe_0 - 2.0 * ta1_y_x_xxy_1[i] * fe_0 + 2.0 * ta1_y_xx_xy_0[i] * fe_0 -
                             2.0 * ta1_y_xx_xy_1[i] * fe_0 + ta1_y_xx_xxy_0[i] * pa_x[i] - ta1_y_xx_xxy_1[i] * pc_x[i];

        ta1_y_xxx_xxz_0[i] = 2.0 * ta1_y_x_xxz_0[i] * fe_0 - 2.0 * ta1_y_x_xxz_1[i] * fe_0 + 2.0 * ta1_y_xx_xz_0[i] * fe_0 -
                             2.0 * ta1_y_xx_xz_1[i] * fe_0 + ta1_y_xx_xxz_0[i] * pa_x[i] - ta1_y_xx_xxz_1[i] * pc_x[i];

        ta1_y_xxx_xyy_0[i] = 2.0 * ta1_y_x_xyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyy_1[i] * fe_0 + ta1_y_xx_yy_0[i] * fe_0 - ta1_y_xx_yy_1[i] * fe_0 +
                             ta1_y_xx_xyy_0[i] * pa_x[i] - ta1_y_xx_xyy_1[i] * pc_x[i];

        ta1_y_xxx_xyz_0[i] = 2.0 * ta1_y_x_xyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyz_1[i] * fe_0 + ta1_y_xx_yz_0[i] * fe_0 - ta1_y_xx_yz_1[i] * fe_0 +
                             ta1_y_xx_xyz_0[i] * pa_x[i] - ta1_y_xx_xyz_1[i] * pc_x[i];

        ta1_y_xxx_xzz_0[i] = 2.0 * ta1_y_x_xzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzz_1[i] * fe_0 + ta1_y_xx_zz_0[i] * fe_0 - ta1_y_xx_zz_1[i] * fe_0 +
                             ta1_y_xx_xzz_0[i] * pa_x[i] - ta1_y_xx_xzz_1[i] * pc_x[i];

        ta1_y_xxx_yyy_0[i] =
            2.0 * ta1_y_x_yyy_0[i] * fe_0 - 2.0 * ta1_y_x_yyy_1[i] * fe_0 + ta1_y_xx_yyy_0[i] * pa_x[i] - ta1_y_xx_yyy_1[i] * pc_x[i];

        ta1_y_xxx_yyz_0[i] =
            2.0 * ta1_y_x_yyz_0[i] * fe_0 - 2.0 * ta1_y_x_yyz_1[i] * fe_0 + ta1_y_xx_yyz_0[i] * pa_x[i] - ta1_y_xx_yyz_1[i] * pc_x[i];

        ta1_y_xxx_yzz_0[i] =
            2.0 * ta1_y_x_yzz_0[i] * fe_0 - 2.0 * ta1_y_x_yzz_1[i] * fe_0 + ta1_y_xx_yzz_0[i] * pa_x[i] - ta1_y_xx_yzz_1[i] * pc_x[i];

        ta1_y_xxx_zzz_0[i] =
            2.0 * ta1_y_x_zzz_0[i] * fe_0 - 2.0 * ta1_y_x_zzz_1[i] * fe_0 + ta1_y_xx_zzz_0[i] * pa_x[i] - ta1_y_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 110-120 components of targeted buffer : FF

    auto ta1_y_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 110);

    auto ta1_y_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 111);

    auto ta1_y_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 112);

    auto ta1_y_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 113);

    auto ta1_y_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 114);

    auto ta1_y_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 115);

    auto ta1_y_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 116);

    auto ta1_y_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 117);

    auto ta1_y_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 118);

    auto ta1_y_xxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 119);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xxx_0,  \
                             ta1_y_xx_xxx_1,  \
                             ta1_y_xx_xxy_0,  \
                             ta1_y_xx_xxy_1,  \
                             ta1_y_xx_xxz_0,  \
                             ta1_y_xx_xxz_1,  \
                             ta1_y_xx_xy_0,   \
                             ta1_y_xx_xy_1,   \
                             ta1_y_xx_xyy_0,  \
                             ta1_y_xx_xyy_1,  \
                             ta1_y_xx_xyz_0,  \
                             ta1_y_xx_xyz_1,  \
                             ta1_y_xx_xz_0,   \
                             ta1_y_xx_xz_1,   \
                             ta1_y_xx_xzz_0,  \
                             ta1_y_xx_xzz_1,  \
                             ta1_y_xx_zzz_0,  \
                             ta1_y_xx_zzz_1,  \
                             ta1_y_xxy_xxx_0, \
                             ta1_y_xxy_xxy_0, \
                             ta1_y_xxy_xxz_0, \
                             ta1_y_xxy_xyy_0, \
                             ta1_y_xxy_xyz_0, \
                             ta1_y_xxy_xzz_0, \
                             ta1_y_xxy_yyy_0, \
                             ta1_y_xxy_yyz_0, \
                             ta1_y_xxy_yzz_0, \
                             ta1_y_xxy_zzz_0, \
                             ta1_y_xy_yyy_0,  \
                             ta1_y_xy_yyy_1,  \
                             ta1_y_xy_yyz_0,  \
                             ta1_y_xy_yyz_1,  \
                             ta1_y_xy_yzz_0,  \
                             ta1_y_xy_yzz_1,  \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyz_0,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yzz_0,   \
                             ta1_y_y_yzz_1,   \
                             ta_xx_xxx_1,     \
                             ta_xx_xxy_1,     \
                             ta_xx_xxz_1,     \
                             ta_xx_xyy_1,     \
                             ta_xx_xyz_1,     \
                             ta_xx_xzz_1,     \
                             ta_xx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_xxx_0[i] = ta_xx_xxx_1[i] + ta1_y_xx_xxx_0[i] * pa_y[i] - ta1_y_xx_xxx_1[i] * pc_y[i];

        ta1_y_xxy_xxy_0[i] =
            ta1_y_xx_xx_0[i] * fe_0 - ta1_y_xx_xx_1[i] * fe_0 + ta_xx_xxy_1[i] + ta1_y_xx_xxy_0[i] * pa_y[i] - ta1_y_xx_xxy_1[i] * pc_y[i];

        ta1_y_xxy_xxz_0[i] = ta_xx_xxz_1[i] + ta1_y_xx_xxz_0[i] * pa_y[i] - ta1_y_xx_xxz_1[i] * pc_y[i];

        ta1_y_xxy_xyy_0[i] = 2.0 * ta1_y_xx_xy_0[i] * fe_0 - 2.0 * ta1_y_xx_xy_1[i] * fe_0 + ta_xx_xyy_1[i] + ta1_y_xx_xyy_0[i] * pa_y[i] -
                             ta1_y_xx_xyy_1[i] * pc_y[i];

        ta1_y_xxy_xyz_0[i] =
            ta1_y_xx_xz_0[i] * fe_0 - ta1_y_xx_xz_1[i] * fe_0 + ta_xx_xyz_1[i] + ta1_y_xx_xyz_0[i] * pa_y[i] - ta1_y_xx_xyz_1[i] * pc_y[i];

        ta1_y_xxy_xzz_0[i] = ta_xx_xzz_1[i] + ta1_y_xx_xzz_0[i] * pa_y[i] - ta1_y_xx_xzz_1[i] * pc_y[i];

        ta1_y_xxy_yyy_0[i] = ta1_y_y_yyy_0[i] * fe_0 - ta1_y_y_yyy_1[i] * fe_0 + ta1_y_xy_yyy_0[i] * pa_x[i] - ta1_y_xy_yyy_1[i] * pc_x[i];

        ta1_y_xxy_yyz_0[i] = ta1_y_y_yyz_0[i] * fe_0 - ta1_y_y_yyz_1[i] * fe_0 + ta1_y_xy_yyz_0[i] * pa_x[i] - ta1_y_xy_yyz_1[i] * pc_x[i];

        ta1_y_xxy_yzz_0[i] = ta1_y_y_yzz_0[i] * fe_0 - ta1_y_y_yzz_1[i] * fe_0 + ta1_y_xy_yzz_0[i] * pa_x[i] - ta1_y_xy_yzz_1[i] * pc_x[i];

        ta1_y_xxy_zzz_0[i] = ta_xx_zzz_1[i] + ta1_y_xx_zzz_0[i] * pa_y[i] - ta1_y_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : FF

    auto ta1_y_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 120);

    auto ta1_y_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 121);

    auto ta1_y_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 122);

    auto ta1_y_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 123);

    auto ta1_y_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 124);

    auto ta1_y_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 125);

    auto ta1_y_xxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 126);

    auto ta1_y_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 127);

    auto ta1_y_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 128);

    auto ta1_y_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 129);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xx_xx_0,   \
                             ta1_y_xx_xx_1,   \
                             ta1_y_xx_xxx_0,  \
                             ta1_y_xx_xxx_1,  \
                             ta1_y_xx_xxy_0,  \
                             ta1_y_xx_xxy_1,  \
                             ta1_y_xx_xxz_0,  \
                             ta1_y_xx_xxz_1,  \
                             ta1_y_xx_xy_0,   \
                             ta1_y_xx_xy_1,   \
                             ta1_y_xx_xyy_0,  \
                             ta1_y_xx_xyy_1,  \
                             ta1_y_xx_xyz_0,  \
                             ta1_y_xx_xyz_1,  \
                             ta1_y_xx_xz_0,   \
                             ta1_y_xx_xz_1,   \
                             ta1_y_xx_xzz_0,  \
                             ta1_y_xx_xzz_1,  \
                             ta1_y_xx_yyy_0,  \
                             ta1_y_xx_yyy_1,  \
                             ta1_y_xxz_xxx_0, \
                             ta1_y_xxz_xxy_0, \
                             ta1_y_xxz_xxz_0, \
                             ta1_y_xxz_xyy_0, \
                             ta1_y_xxz_xyz_0, \
                             ta1_y_xxz_xzz_0, \
                             ta1_y_xxz_yyy_0, \
                             ta1_y_xxz_yyz_0, \
                             ta1_y_xxz_yzz_0, \
                             ta1_y_xxz_zzz_0, \
                             ta1_y_xz_yyz_0,  \
                             ta1_y_xz_yyz_1,  \
                             ta1_y_xz_yzz_0,  \
                             ta1_y_xz_yzz_1,  \
                             ta1_y_xz_zzz_0,  \
                             ta1_y_xz_zzz_1,  \
                             ta1_y_z_yyz_0,   \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yzz_0,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_zzz_0,   \
                             ta1_y_z_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_xxx_0[i] = ta1_y_xx_xxx_0[i] * pa_z[i] - ta1_y_xx_xxx_1[i] * pc_z[i];

        ta1_y_xxz_xxy_0[i] = ta1_y_xx_xxy_0[i] * pa_z[i] - ta1_y_xx_xxy_1[i] * pc_z[i];

        ta1_y_xxz_xxz_0[i] = ta1_y_xx_xx_0[i] * fe_0 - ta1_y_xx_xx_1[i] * fe_0 + ta1_y_xx_xxz_0[i] * pa_z[i] - ta1_y_xx_xxz_1[i] * pc_z[i];

        ta1_y_xxz_xyy_0[i] = ta1_y_xx_xyy_0[i] * pa_z[i] - ta1_y_xx_xyy_1[i] * pc_z[i];

        ta1_y_xxz_xyz_0[i] = ta1_y_xx_xy_0[i] * fe_0 - ta1_y_xx_xy_1[i] * fe_0 + ta1_y_xx_xyz_0[i] * pa_z[i] - ta1_y_xx_xyz_1[i] * pc_z[i];

        ta1_y_xxz_xzz_0[i] =
            2.0 * ta1_y_xx_xz_0[i] * fe_0 - 2.0 * ta1_y_xx_xz_1[i] * fe_0 + ta1_y_xx_xzz_0[i] * pa_z[i] - ta1_y_xx_xzz_1[i] * pc_z[i];

        ta1_y_xxz_yyy_0[i] = ta1_y_xx_yyy_0[i] * pa_z[i] - ta1_y_xx_yyy_1[i] * pc_z[i];

        ta1_y_xxz_yyz_0[i] = ta1_y_z_yyz_0[i] * fe_0 - ta1_y_z_yyz_1[i] * fe_0 + ta1_y_xz_yyz_0[i] * pa_x[i] - ta1_y_xz_yyz_1[i] * pc_x[i];

        ta1_y_xxz_yzz_0[i] = ta1_y_z_yzz_0[i] * fe_0 - ta1_y_z_yzz_1[i] * fe_0 + ta1_y_xz_yzz_0[i] * pa_x[i] - ta1_y_xz_yzz_1[i] * pc_x[i];

        ta1_y_xxz_zzz_0[i] = ta1_y_z_zzz_0[i] * fe_0 - ta1_y_z_zzz_1[i] * fe_0 + ta1_y_xz_zzz_0[i] * pa_x[i] - ta1_y_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : FF

    auto ta1_y_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 130);

    auto ta1_y_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 131);

    auto ta1_y_xyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 132);

    auto ta1_y_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 133);

    auto ta1_y_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 134);

    auto ta1_y_xyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 135);

    auto ta1_y_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 136);

    auto ta1_y_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 137);

    auto ta1_y_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 138);

    auto ta1_y_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 139);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xyy_xxx_0, \
                             ta1_y_xyy_xxy_0, \
                             ta1_y_xyy_xxz_0, \
                             ta1_y_xyy_xyy_0, \
                             ta1_y_xyy_xyz_0, \
                             ta1_y_xyy_xzz_0, \
                             ta1_y_xyy_yyy_0, \
                             ta1_y_xyy_yyz_0, \
                             ta1_y_xyy_yzz_0, \
                             ta1_y_xyy_zzz_0, \
                             ta1_y_yy_xx_0,   \
                             ta1_y_yy_xx_1,   \
                             ta1_y_yy_xxx_0,  \
                             ta1_y_yy_xxx_1,  \
                             ta1_y_yy_xxy_0,  \
                             ta1_y_yy_xxy_1,  \
                             ta1_y_yy_xxz_0,  \
                             ta1_y_yy_xxz_1,  \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_xyy_0,  \
                             ta1_y_yy_xyy_1,  \
                             ta1_y_yy_xyz_0,  \
                             ta1_y_yy_xyz_1,  \
                             ta1_y_yy_xz_0,   \
                             ta1_y_yy_xz_1,   \
                             ta1_y_yy_xzz_0,  \
                             ta1_y_yy_xzz_1,  \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yyy_0,  \
                             ta1_y_yy_yyy_1,  \
                             ta1_y_yy_yyz_0,  \
                             ta1_y_yy_yyz_1,  \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yy_yzz_0,  \
                             ta1_y_yy_yzz_1,  \
                             ta1_y_yy_zz_0,   \
                             ta1_y_yy_zz_1,   \
                             ta1_y_yy_zzz_0,  \
                             ta1_y_yy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_xxx_0[i] =
            3.0 * ta1_y_yy_xx_0[i] * fe_0 - 3.0 * ta1_y_yy_xx_1[i] * fe_0 + ta1_y_yy_xxx_0[i] * pa_x[i] - ta1_y_yy_xxx_1[i] * pc_x[i];

        ta1_y_xyy_xxy_0[i] =
            2.0 * ta1_y_yy_xy_0[i] * fe_0 - 2.0 * ta1_y_yy_xy_1[i] * fe_0 + ta1_y_yy_xxy_0[i] * pa_x[i] - ta1_y_yy_xxy_1[i] * pc_x[i];

        ta1_y_xyy_xxz_0[i] =
            2.0 * ta1_y_yy_xz_0[i] * fe_0 - 2.0 * ta1_y_yy_xz_1[i] * fe_0 + ta1_y_yy_xxz_0[i] * pa_x[i] - ta1_y_yy_xxz_1[i] * pc_x[i];

        ta1_y_xyy_xyy_0[i] = ta1_y_yy_yy_0[i] * fe_0 - ta1_y_yy_yy_1[i] * fe_0 + ta1_y_yy_xyy_0[i] * pa_x[i] - ta1_y_yy_xyy_1[i] * pc_x[i];

        ta1_y_xyy_xyz_0[i] = ta1_y_yy_yz_0[i] * fe_0 - ta1_y_yy_yz_1[i] * fe_0 + ta1_y_yy_xyz_0[i] * pa_x[i] - ta1_y_yy_xyz_1[i] * pc_x[i];

        ta1_y_xyy_xzz_0[i] = ta1_y_yy_zz_0[i] * fe_0 - ta1_y_yy_zz_1[i] * fe_0 + ta1_y_yy_xzz_0[i] * pa_x[i] - ta1_y_yy_xzz_1[i] * pc_x[i];

        ta1_y_xyy_yyy_0[i] = ta1_y_yy_yyy_0[i] * pa_x[i] - ta1_y_yy_yyy_1[i] * pc_x[i];

        ta1_y_xyy_yyz_0[i] = ta1_y_yy_yyz_0[i] * pa_x[i] - ta1_y_yy_yyz_1[i] * pc_x[i];

        ta1_y_xyy_yzz_0[i] = ta1_y_yy_yzz_0[i] * pa_x[i] - ta1_y_yy_yzz_1[i] * pc_x[i];

        ta1_y_xyy_zzz_0[i] = ta1_y_yy_zzz_0[i] * pa_x[i] - ta1_y_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 140-150 components of targeted buffer : FF

    auto ta1_y_xyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 140);

    auto ta1_y_xyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 141);

    auto ta1_y_xyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 142);

    auto ta1_y_xyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 143);

    auto ta1_y_xyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 144);

    auto ta1_y_xyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 145);

    auto ta1_y_xyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 146);

    auto ta1_y_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 147);

    auto ta1_y_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 148);

    auto ta1_y_xyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_xy_xxx_0,  \
                             ta1_y_xy_xxx_1,  \
                             ta1_y_xy_xxy_0,  \
                             ta1_y_xy_xxy_1,  \
                             ta1_y_xy_xyy_0,  \
                             ta1_y_xy_xyy_1,  \
                             ta1_y_xyz_xxx_0, \
                             ta1_y_xyz_xxy_0, \
                             ta1_y_xyz_xxz_0, \
                             ta1_y_xyz_xyy_0, \
                             ta1_y_xyz_xyz_0, \
                             ta1_y_xyz_xzz_0, \
                             ta1_y_xyz_yyy_0, \
                             ta1_y_xyz_yyz_0, \
                             ta1_y_xyz_yzz_0, \
                             ta1_y_xyz_zzz_0, \
                             ta1_y_xz_xxz_0,  \
                             ta1_y_xz_xxz_1,  \
                             ta1_y_xz_xzz_0,  \
                             ta1_y_xz_xzz_1,  \
                             ta1_y_yz_xyz_0,  \
                             ta1_y_yz_xyz_1,  \
                             ta1_y_yz_yyy_0,  \
                             ta1_y_yz_yyy_1,  \
                             ta1_y_yz_yyz_0,  \
                             ta1_y_yz_yyz_1,  \
                             ta1_y_yz_yz_0,   \
                             ta1_y_yz_yz_1,   \
                             ta1_y_yz_yzz_0,  \
                             ta1_y_yz_yzz_1,  \
                             ta1_y_yz_zzz_0,  \
                             ta1_y_yz_zzz_1,  \
                             ta_xz_xxz_1,     \
                             ta_xz_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyz_xxx_0[i] = ta1_y_xy_xxx_0[i] * pa_z[i] - ta1_y_xy_xxx_1[i] * pc_z[i];

        ta1_y_xyz_xxy_0[i] = ta1_y_xy_xxy_0[i] * pa_z[i] - ta1_y_xy_xxy_1[i] * pc_z[i];

        ta1_y_xyz_xxz_0[i] = ta_xz_xxz_1[i] + ta1_y_xz_xxz_0[i] * pa_y[i] - ta1_y_xz_xxz_1[i] * pc_y[i];

        ta1_y_xyz_xyy_0[i] = ta1_y_xy_xyy_0[i] * pa_z[i] - ta1_y_xy_xyy_1[i] * pc_z[i];

        ta1_y_xyz_xyz_0[i] = ta1_y_yz_yz_0[i] * fe_0 - ta1_y_yz_yz_1[i] * fe_0 + ta1_y_yz_xyz_0[i] * pa_x[i] - ta1_y_yz_xyz_1[i] * pc_x[i];

        ta1_y_xyz_xzz_0[i] = ta_xz_xzz_1[i] + ta1_y_xz_xzz_0[i] * pa_y[i] - ta1_y_xz_xzz_1[i] * pc_y[i];

        ta1_y_xyz_yyy_0[i] = ta1_y_yz_yyy_0[i] * pa_x[i] - ta1_y_yz_yyy_1[i] * pc_x[i];

        ta1_y_xyz_yyz_0[i] = ta1_y_yz_yyz_0[i] * pa_x[i] - ta1_y_yz_yyz_1[i] * pc_x[i];

        ta1_y_xyz_yzz_0[i] = ta1_y_yz_yzz_0[i] * pa_x[i] - ta1_y_yz_yzz_1[i] * pc_x[i];

        ta1_y_xyz_zzz_0[i] = ta1_y_yz_zzz_0[i] * pa_x[i] - ta1_y_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : FF

    auto ta1_y_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 150);

    auto ta1_y_xzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 151);

    auto ta1_y_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 152);

    auto ta1_y_xzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 153);

    auto ta1_y_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 154);

    auto ta1_y_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 155);

    auto ta1_y_xzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 156);

    auto ta1_y_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 157);

    auto ta1_y_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 158);

    auto ta1_y_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 159);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xzz_xxx_0, \
                             ta1_y_xzz_xxy_0, \
                             ta1_y_xzz_xxz_0, \
                             ta1_y_xzz_xyy_0, \
                             ta1_y_xzz_xyz_0, \
                             ta1_y_xzz_xzz_0, \
                             ta1_y_xzz_yyy_0, \
                             ta1_y_xzz_yyz_0, \
                             ta1_y_xzz_yzz_0, \
                             ta1_y_xzz_zzz_0, \
                             ta1_y_zz_xx_0,   \
                             ta1_y_zz_xx_1,   \
                             ta1_y_zz_xxx_0,  \
                             ta1_y_zz_xxx_1,  \
                             ta1_y_zz_xxy_0,  \
                             ta1_y_zz_xxy_1,  \
                             ta1_y_zz_xxz_0,  \
                             ta1_y_zz_xxz_1,  \
                             ta1_y_zz_xy_0,   \
                             ta1_y_zz_xy_1,   \
                             ta1_y_zz_xyy_0,  \
                             ta1_y_zz_xyy_1,  \
                             ta1_y_zz_xyz_0,  \
                             ta1_y_zz_xyz_1,  \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_xzz_0,  \
                             ta1_y_zz_xzz_1,  \
                             ta1_y_zz_yy_0,   \
                             ta1_y_zz_yy_1,   \
                             ta1_y_zz_yyy_0,  \
                             ta1_y_zz_yyy_1,  \
                             ta1_y_zz_yyz_0,  \
                             ta1_y_zz_yyz_1,  \
                             ta1_y_zz_yz_0,   \
                             ta1_y_zz_yz_1,   \
                             ta1_y_zz_yzz_0,  \
                             ta1_y_zz_yzz_1,  \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             ta1_y_zz_zzz_0,  \
                             ta1_y_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_xxx_0[i] =
            3.0 * ta1_y_zz_xx_0[i] * fe_0 - 3.0 * ta1_y_zz_xx_1[i] * fe_0 + ta1_y_zz_xxx_0[i] * pa_x[i] - ta1_y_zz_xxx_1[i] * pc_x[i];

        ta1_y_xzz_xxy_0[i] =
            2.0 * ta1_y_zz_xy_0[i] * fe_0 - 2.0 * ta1_y_zz_xy_1[i] * fe_0 + ta1_y_zz_xxy_0[i] * pa_x[i] - ta1_y_zz_xxy_1[i] * pc_x[i];

        ta1_y_xzz_xxz_0[i] =
            2.0 * ta1_y_zz_xz_0[i] * fe_0 - 2.0 * ta1_y_zz_xz_1[i] * fe_0 + ta1_y_zz_xxz_0[i] * pa_x[i] - ta1_y_zz_xxz_1[i] * pc_x[i];

        ta1_y_xzz_xyy_0[i] = ta1_y_zz_yy_0[i] * fe_0 - ta1_y_zz_yy_1[i] * fe_0 + ta1_y_zz_xyy_0[i] * pa_x[i] - ta1_y_zz_xyy_1[i] * pc_x[i];

        ta1_y_xzz_xyz_0[i] = ta1_y_zz_yz_0[i] * fe_0 - ta1_y_zz_yz_1[i] * fe_0 + ta1_y_zz_xyz_0[i] * pa_x[i] - ta1_y_zz_xyz_1[i] * pc_x[i];

        ta1_y_xzz_xzz_0[i] = ta1_y_zz_zz_0[i] * fe_0 - ta1_y_zz_zz_1[i] * fe_0 + ta1_y_zz_xzz_0[i] * pa_x[i] - ta1_y_zz_xzz_1[i] * pc_x[i];

        ta1_y_xzz_yyy_0[i] = ta1_y_zz_yyy_0[i] * pa_x[i] - ta1_y_zz_yyy_1[i] * pc_x[i];

        ta1_y_xzz_yyz_0[i] = ta1_y_zz_yyz_0[i] * pa_x[i] - ta1_y_zz_yyz_1[i] * pc_x[i];

        ta1_y_xzz_yzz_0[i] = ta1_y_zz_yzz_0[i] * pa_x[i] - ta1_y_zz_yzz_1[i] * pc_x[i];

        ta1_y_xzz_zzz_0[i] = ta1_y_zz_zzz_0[i] * pa_x[i] - ta1_y_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : FF

    auto ta1_y_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 160);

    auto ta1_y_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 161);

    auto ta1_y_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 162);

    auto ta1_y_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 163);

    auto ta1_y_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 164);

    auto ta1_y_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 165);

    auto ta1_y_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 166);

    auto ta1_y_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 167);

    auto ta1_y_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 168);

    auto ta1_y_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 169);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_y_xxx_0,   \
                             ta1_y_y_xxx_1,   \
                             ta1_y_y_xxy_0,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xxz_0,   \
                             ta1_y_y_xxz_1,   \
                             ta1_y_y_xyy_0,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyz_0,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_xzz_0,   \
                             ta1_y_y_xzz_1,   \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyz_0,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yzz_0,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_zzz_0,   \
                             ta1_y_y_zzz_1,   \
                             ta1_y_yy_xx_0,   \
                             ta1_y_yy_xx_1,   \
                             ta1_y_yy_xxx_0,  \
                             ta1_y_yy_xxx_1,  \
                             ta1_y_yy_xxy_0,  \
                             ta1_y_yy_xxy_1,  \
                             ta1_y_yy_xxz_0,  \
                             ta1_y_yy_xxz_1,  \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_xyy_0,  \
                             ta1_y_yy_xyy_1,  \
                             ta1_y_yy_xyz_0,  \
                             ta1_y_yy_xyz_1,  \
                             ta1_y_yy_xz_0,   \
                             ta1_y_yy_xz_1,   \
                             ta1_y_yy_xzz_0,  \
                             ta1_y_yy_xzz_1,  \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yyy_0,  \
                             ta1_y_yy_yyy_1,  \
                             ta1_y_yy_yyz_0,  \
                             ta1_y_yy_yyz_1,  \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yy_yzz_0,  \
                             ta1_y_yy_yzz_1,  \
                             ta1_y_yy_zz_0,   \
                             ta1_y_yy_zz_1,   \
                             ta1_y_yy_zzz_0,  \
                             ta1_y_yy_zzz_1,  \
                             ta1_y_yyy_xxx_0, \
                             ta1_y_yyy_xxy_0, \
                             ta1_y_yyy_xxz_0, \
                             ta1_y_yyy_xyy_0, \
                             ta1_y_yyy_xyz_0, \
                             ta1_y_yyy_xzz_0, \
                             ta1_y_yyy_yyy_0, \
                             ta1_y_yyy_yyz_0, \
                             ta1_y_yyy_yzz_0, \
                             ta1_y_yyy_zzz_0, \
                             ta_yy_xxx_1,     \
                             ta_yy_xxy_1,     \
                             ta_yy_xxz_1,     \
                             ta_yy_xyy_1,     \
                             ta_yy_xyz_1,     \
                             ta_yy_xzz_1,     \
                             ta_yy_yyy_1,     \
                             ta_yy_yyz_1,     \
                             ta_yy_yzz_1,     \
                             ta_yy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_xxx_0[i] = 2.0 * ta1_y_y_xxx_0[i] * fe_0 - 2.0 * ta1_y_y_xxx_1[i] * fe_0 + ta_yy_xxx_1[i] + ta1_y_yy_xxx_0[i] * pa_y[i] -
                             ta1_y_yy_xxx_1[i] * pc_y[i];

        ta1_y_yyy_xxy_0[i] = 2.0 * ta1_y_y_xxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxy_1[i] * fe_0 + ta1_y_yy_xx_0[i] * fe_0 - ta1_y_yy_xx_1[i] * fe_0 +
                             ta_yy_xxy_1[i] + ta1_y_yy_xxy_0[i] * pa_y[i] - ta1_y_yy_xxy_1[i] * pc_y[i];

        ta1_y_yyy_xxz_0[i] = 2.0 * ta1_y_y_xxz_0[i] * fe_0 - 2.0 * ta1_y_y_xxz_1[i] * fe_0 + ta_yy_xxz_1[i] + ta1_y_yy_xxz_0[i] * pa_y[i] -
                             ta1_y_yy_xxz_1[i] * pc_y[i];

        ta1_y_yyy_xyy_0[i] = 2.0 * ta1_y_y_xyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyy_1[i] * fe_0 + 2.0 * ta1_y_yy_xy_0[i] * fe_0 -
                             2.0 * ta1_y_yy_xy_1[i] * fe_0 + ta_yy_xyy_1[i] + ta1_y_yy_xyy_0[i] * pa_y[i] - ta1_y_yy_xyy_1[i] * pc_y[i];

        ta1_y_yyy_xyz_0[i] = 2.0 * ta1_y_y_xyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyz_1[i] * fe_0 + ta1_y_yy_xz_0[i] * fe_0 - ta1_y_yy_xz_1[i] * fe_0 +
                             ta_yy_xyz_1[i] + ta1_y_yy_xyz_0[i] * pa_y[i] - ta1_y_yy_xyz_1[i] * pc_y[i];

        ta1_y_yyy_xzz_0[i] = 2.0 * ta1_y_y_xzz_0[i] * fe_0 - 2.0 * ta1_y_y_xzz_1[i] * fe_0 + ta_yy_xzz_1[i] + ta1_y_yy_xzz_0[i] * pa_y[i] -
                             ta1_y_yy_xzz_1[i] * pc_y[i];

        ta1_y_yyy_yyy_0[i] = 2.0 * ta1_y_y_yyy_0[i] * fe_0 - 2.0 * ta1_y_y_yyy_1[i] * fe_0 + 3.0 * ta1_y_yy_yy_0[i] * fe_0 -
                             3.0 * ta1_y_yy_yy_1[i] * fe_0 + ta_yy_yyy_1[i] + ta1_y_yy_yyy_0[i] * pa_y[i] - ta1_y_yy_yyy_1[i] * pc_y[i];

        ta1_y_yyy_yyz_0[i] = 2.0 * ta1_y_y_yyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyz_1[i] * fe_0 + 2.0 * ta1_y_yy_yz_0[i] * fe_0 -
                             2.0 * ta1_y_yy_yz_1[i] * fe_0 + ta_yy_yyz_1[i] + ta1_y_yy_yyz_0[i] * pa_y[i] - ta1_y_yy_yyz_1[i] * pc_y[i];

        ta1_y_yyy_yzz_0[i] = 2.0 * ta1_y_y_yzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzz_1[i] * fe_0 + ta1_y_yy_zz_0[i] * fe_0 - ta1_y_yy_zz_1[i] * fe_0 +
                             ta_yy_yzz_1[i] + ta1_y_yy_yzz_0[i] * pa_y[i] - ta1_y_yy_yzz_1[i] * pc_y[i];

        ta1_y_yyy_zzz_0[i] = 2.0 * ta1_y_y_zzz_0[i] * fe_0 - 2.0 * ta1_y_y_zzz_1[i] * fe_0 + ta_yy_zzz_1[i] + ta1_y_yy_zzz_0[i] * pa_y[i] -
                             ta1_y_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : FF

    auto ta1_y_yyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 170);

    auto ta1_y_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 171);

    auto ta1_y_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 172);

    auto ta1_y_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 173);

    auto ta1_y_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 174);

    auto ta1_y_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 175);

    auto ta1_y_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 176);

    auto ta1_y_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 177);

    auto ta1_y_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 178);

    auto ta1_y_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 179);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_yy_xx_0,   \
                             ta1_y_yy_xx_1,   \
                             ta1_y_yy_xxx_0,  \
                             ta1_y_yy_xxx_1,  \
                             ta1_y_yy_xxy_0,  \
                             ta1_y_yy_xxy_1,  \
                             ta1_y_yy_xxz_0,  \
                             ta1_y_yy_xxz_1,  \
                             ta1_y_yy_xy_0,   \
                             ta1_y_yy_xy_1,   \
                             ta1_y_yy_xyy_0,  \
                             ta1_y_yy_xyy_1,  \
                             ta1_y_yy_xyz_0,  \
                             ta1_y_yy_xyz_1,  \
                             ta1_y_yy_xz_0,   \
                             ta1_y_yy_xz_1,   \
                             ta1_y_yy_xzz_0,  \
                             ta1_y_yy_xzz_1,  \
                             ta1_y_yy_yy_0,   \
                             ta1_y_yy_yy_1,   \
                             ta1_y_yy_yyy_0,  \
                             ta1_y_yy_yyy_1,  \
                             ta1_y_yy_yyz_0,  \
                             ta1_y_yy_yyz_1,  \
                             ta1_y_yy_yz_0,   \
                             ta1_y_yy_yz_1,   \
                             ta1_y_yy_yzz_0,  \
                             ta1_y_yy_yzz_1,  \
                             ta1_y_yy_zz_0,   \
                             ta1_y_yy_zz_1,   \
                             ta1_y_yy_zzz_0,  \
                             ta1_y_yy_zzz_1,  \
                             ta1_y_yyz_xxx_0, \
                             ta1_y_yyz_xxy_0, \
                             ta1_y_yyz_xxz_0, \
                             ta1_y_yyz_xyy_0, \
                             ta1_y_yyz_xyz_0, \
                             ta1_y_yyz_xzz_0, \
                             ta1_y_yyz_yyy_0, \
                             ta1_y_yyz_yyz_0, \
                             ta1_y_yyz_yzz_0, \
                             ta1_y_yyz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_xxx_0[i] = ta1_y_yy_xxx_0[i] * pa_z[i] - ta1_y_yy_xxx_1[i] * pc_z[i];

        ta1_y_yyz_xxy_0[i] = ta1_y_yy_xxy_0[i] * pa_z[i] - ta1_y_yy_xxy_1[i] * pc_z[i];

        ta1_y_yyz_xxz_0[i] = ta1_y_yy_xx_0[i] * fe_0 - ta1_y_yy_xx_1[i] * fe_0 + ta1_y_yy_xxz_0[i] * pa_z[i] - ta1_y_yy_xxz_1[i] * pc_z[i];

        ta1_y_yyz_xyy_0[i] = ta1_y_yy_xyy_0[i] * pa_z[i] - ta1_y_yy_xyy_1[i] * pc_z[i];

        ta1_y_yyz_xyz_0[i] = ta1_y_yy_xy_0[i] * fe_0 - ta1_y_yy_xy_1[i] * fe_0 + ta1_y_yy_xyz_0[i] * pa_z[i] - ta1_y_yy_xyz_1[i] * pc_z[i];

        ta1_y_yyz_xzz_0[i] =
            2.0 * ta1_y_yy_xz_0[i] * fe_0 - 2.0 * ta1_y_yy_xz_1[i] * fe_0 + ta1_y_yy_xzz_0[i] * pa_z[i] - ta1_y_yy_xzz_1[i] * pc_z[i];

        ta1_y_yyz_yyy_0[i] = ta1_y_yy_yyy_0[i] * pa_z[i] - ta1_y_yy_yyy_1[i] * pc_z[i];

        ta1_y_yyz_yyz_0[i] = ta1_y_yy_yy_0[i] * fe_0 - ta1_y_yy_yy_1[i] * fe_0 + ta1_y_yy_yyz_0[i] * pa_z[i] - ta1_y_yy_yyz_1[i] * pc_z[i];

        ta1_y_yyz_yzz_0[i] =
            2.0 * ta1_y_yy_yz_0[i] * fe_0 - 2.0 * ta1_y_yy_yz_1[i] * fe_0 + ta1_y_yy_yzz_0[i] * pa_z[i] - ta1_y_yy_yzz_1[i] * pc_z[i];

        ta1_y_yyz_zzz_0[i] =
            3.0 * ta1_y_yy_zz_0[i] * fe_0 - 3.0 * ta1_y_yy_zz_1[i] * fe_0 + ta1_y_yy_zzz_0[i] * pa_z[i] - ta1_y_yy_zzz_1[i] * pc_z[i];
    }

    // Set up 180-190 components of targeted buffer : FF

    auto ta1_y_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 180);

    auto ta1_y_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 181);

    auto ta1_y_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 182);

    auto ta1_y_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 183);

    auto ta1_y_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 184);

    auto ta1_y_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 185);

    auto ta1_y_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 186);

    auto ta1_y_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 187);

    auto ta1_y_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 188);

    auto ta1_y_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 189);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_y_xxy_0,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xyy_0,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_yz_xxy_0,  \
                             ta1_y_yz_xxy_1,  \
                             ta1_y_yz_xyy_0,  \
                             ta1_y_yz_xyy_1,  \
                             ta1_y_yz_yyy_0,  \
                             ta1_y_yz_yyy_1,  \
                             ta1_y_yzz_xxx_0, \
                             ta1_y_yzz_xxy_0, \
                             ta1_y_yzz_xxz_0, \
                             ta1_y_yzz_xyy_0, \
                             ta1_y_yzz_xyz_0, \
                             ta1_y_yzz_xzz_0, \
                             ta1_y_yzz_yyy_0, \
                             ta1_y_yzz_yyz_0, \
                             ta1_y_yzz_yzz_0, \
                             ta1_y_yzz_zzz_0, \
                             ta1_y_zz_xxx_0,  \
                             ta1_y_zz_xxx_1,  \
                             ta1_y_zz_xxz_0,  \
                             ta1_y_zz_xxz_1,  \
                             ta1_y_zz_xyz_0,  \
                             ta1_y_zz_xyz_1,  \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_xzz_0,  \
                             ta1_y_zz_xzz_1,  \
                             ta1_y_zz_yyz_0,  \
                             ta1_y_zz_yyz_1,  \
                             ta1_y_zz_yz_0,   \
                             ta1_y_zz_yz_1,   \
                             ta1_y_zz_yzz_0,  \
                             ta1_y_zz_yzz_1,  \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             ta1_y_zz_zzz_0,  \
                             ta1_y_zz_zzz_1,  \
                             ta_zz_xxx_1,     \
                             ta_zz_xxz_1,     \
                             ta_zz_xyz_1,     \
                             ta_zz_xzz_1,     \
                             ta_zz_yyz_1,     \
                             ta_zz_yzz_1,     \
                             ta_zz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_xxx_0[i] = ta_zz_xxx_1[i] + ta1_y_zz_xxx_0[i] * pa_y[i] - ta1_y_zz_xxx_1[i] * pc_y[i];

        ta1_y_yzz_xxy_0[i] = ta1_y_y_xxy_0[i] * fe_0 - ta1_y_y_xxy_1[i] * fe_0 + ta1_y_yz_xxy_0[i] * pa_z[i] - ta1_y_yz_xxy_1[i] * pc_z[i];

        ta1_y_yzz_xxz_0[i] = ta_zz_xxz_1[i] + ta1_y_zz_xxz_0[i] * pa_y[i] - ta1_y_zz_xxz_1[i] * pc_y[i];

        ta1_y_yzz_xyy_0[i] = ta1_y_y_xyy_0[i] * fe_0 - ta1_y_y_xyy_1[i] * fe_0 + ta1_y_yz_xyy_0[i] * pa_z[i] - ta1_y_yz_xyy_1[i] * pc_z[i];

        ta1_y_yzz_xyz_0[i] =
            ta1_y_zz_xz_0[i] * fe_0 - ta1_y_zz_xz_1[i] * fe_0 + ta_zz_xyz_1[i] + ta1_y_zz_xyz_0[i] * pa_y[i] - ta1_y_zz_xyz_1[i] * pc_y[i];

        ta1_y_yzz_xzz_0[i] = ta_zz_xzz_1[i] + ta1_y_zz_xzz_0[i] * pa_y[i] - ta1_y_zz_xzz_1[i] * pc_y[i];

        ta1_y_yzz_yyy_0[i] = ta1_y_y_yyy_0[i] * fe_0 - ta1_y_y_yyy_1[i] * fe_0 + ta1_y_yz_yyy_0[i] * pa_z[i] - ta1_y_yz_yyy_1[i] * pc_z[i];

        ta1_y_yzz_yyz_0[i] = 2.0 * ta1_y_zz_yz_0[i] * fe_0 - 2.0 * ta1_y_zz_yz_1[i] * fe_0 + ta_zz_yyz_1[i] + ta1_y_zz_yyz_0[i] * pa_y[i] -
                             ta1_y_zz_yyz_1[i] * pc_y[i];

        ta1_y_yzz_yzz_0[i] =
            ta1_y_zz_zz_0[i] * fe_0 - ta1_y_zz_zz_1[i] * fe_0 + ta_zz_yzz_1[i] + ta1_y_zz_yzz_0[i] * pa_y[i] - ta1_y_zz_yzz_1[i] * pc_y[i];

        ta1_y_yzz_zzz_0[i] = ta_zz_zzz_1[i] + ta1_y_zz_zzz_0[i] * pa_y[i] - ta1_y_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 190-200 components of targeted buffer : FF

    auto ta1_y_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 190);

    auto ta1_y_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 191);

    auto ta1_y_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 192);

    auto ta1_y_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 193);

    auto ta1_y_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 194);

    auto ta1_y_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 195);

    auto ta1_y_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 196);

    auto ta1_y_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 197);

    auto ta1_y_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 198);

    auto ta1_y_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 199);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_z_xxx_0,   \
                             ta1_y_z_xxx_1,   \
                             ta1_y_z_xxy_0,   \
                             ta1_y_z_xxy_1,   \
                             ta1_y_z_xxz_0,   \
                             ta1_y_z_xxz_1,   \
                             ta1_y_z_xyy_0,   \
                             ta1_y_z_xyy_1,   \
                             ta1_y_z_xyz_0,   \
                             ta1_y_z_xyz_1,   \
                             ta1_y_z_xzz_0,   \
                             ta1_y_z_xzz_1,   \
                             ta1_y_z_yyy_0,   \
                             ta1_y_z_yyy_1,   \
                             ta1_y_z_yyz_0,   \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yzz_0,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_zzz_0,   \
                             ta1_y_z_zzz_1,   \
                             ta1_y_zz_xx_0,   \
                             ta1_y_zz_xx_1,   \
                             ta1_y_zz_xxx_0,  \
                             ta1_y_zz_xxx_1,  \
                             ta1_y_zz_xxy_0,  \
                             ta1_y_zz_xxy_1,  \
                             ta1_y_zz_xxz_0,  \
                             ta1_y_zz_xxz_1,  \
                             ta1_y_zz_xy_0,   \
                             ta1_y_zz_xy_1,   \
                             ta1_y_zz_xyy_0,  \
                             ta1_y_zz_xyy_1,  \
                             ta1_y_zz_xyz_0,  \
                             ta1_y_zz_xyz_1,  \
                             ta1_y_zz_xz_0,   \
                             ta1_y_zz_xz_1,   \
                             ta1_y_zz_xzz_0,  \
                             ta1_y_zz_xzz_1,  \
                             ta1_y_zz_yy_0,   \
                             ta1_y_zz_yy_1,   \
                             ta1_y_zz_yyy_0,  \
                             ta1_y_zz_yyy_1,  \
                             ta1_y_zz_yyz_0,  \
                             ta1_y_zz_yyz_1,  \
                             ta1_y_zz_yz_0,   \
                             ta1_y_zz_yz_1,   \
                             ta1_y_zz_yzz_0,  \
                             ta1_y_zz_yzz_1,  \
                             ta1_y_zz_zz_0,   \
                             ta1_y_zz_zz_1,   \
                             ta1_y_zz_zzz_0,  \
                             ta1_y_zz_zzz_1,  \
                             ta1_y_zzz_xxx_0, \
                             ta1_y_zzz_xxy_0, \
                             ta1_y_zzz_xxz_0, \
                             ta1_y_zzz_xyy_0, \
                             ta1_y_zzz_xyz_0, \
                             ta1_y_zzz_xzz_0, \
                             ta1_y_zzz_yyy_0, \
                             ta1_y_zzz_yyz_0, \
                             ta1_y_zzz_yzz_0, \
                             ta1_y_zzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_xxx_0[i] =
            2.0 * ta1_y_z_xxx_0[i] * fe_0 - 2.0 * ta1_y_z_xxx_1[i] * fe_0 + ta1_y_zz_xxx_0[i] * pa_z[i] - ta1_y_zz_xxx_1[i] * pc_z[i];

        ta1_y_zzz_xxy_0[i] =
            2.0 * ta1_y_z_xxy_0[i] * fe_0 - 2.0 * ta1_y_z_xxy_1[i] * fe_0 + ta1_y_zz_xxy_0[i] * pa_z[i] - ta1_y_zz_xxy_1[i] * pc_z[i];

        ta1_y_zzz_xxz_0[i] = 2.0 * ta1_y_z_xxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxz_1[i] * fe_0 + ta1_y_zz_xx_0[i] * fe_0 - ta1_y_zz_xx_1[i] * fe_0 +
                             ta1_y_zz_xxz_0[i] * pa_z[i] - ta1_y_zz_xxz_1[i] * pc_z[i];

        ta1_y_zzz_xyy_0[i] =
            2.0 * ta1_y_z_xyy_0[i] * fe_0 - 2.0 * ta1_y_z_xyy_1[i] * fe_0 + ta1_y_zz_xyy_0[i] * pa_z[i] - ta1_y_zz_xyy_1[i] * pc_z[i];

        ta1_y_zzz_xyz_0[i] = 2.0 * ta1_y_z_xyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyz_1[i] * fe_0 + ta1_y_zz_xy_0[i] * fe_0 - ta1_y_zz_xy_1[i] * fe_0 +
                             ta1_y_zz_xyz_0[i] * pa_z[i] - ta1_y_zz_xyz_1[i] * pc_z[i];

        ta1_y_zzz_xzz_0[i] = 2.0 * ta1_y_z_xzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xz_0[i] * fe_0 -
                             2.0 * ta1_y_zz_xz_1[i] * fe_0 + ta1_y_zz_xzz_0[i] * pa_z[i] - ta1_y_zz_xzz_1[i] * pc_z[i];

        ta1_y_zzz_yyy_0[i] =
            2.0 * ta1_y_z_yyy_0[i] * fe_0 - 2.0 * ta1_y_z_yyy_1[i] * fe_0 + ta1_y_zz_yyy_0[i] * pa_z[i] - ta1_y_zz_yyy_1[i] * pc_z[i];

        ta1_y_zzz_yyz_0[i] = 2.0 * ta1_y_z_yyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyz_1[i] * fe_0 + ta1_y_zz_yy_0[i] * fe_0 - ta1_y_zz_yy_1[i] * fe_0 +
                             ta1_y_zz_yyz_0[i] * pa_z[i] - ta1_y_zz_yyz_1[i] * pc_z[i];

        ta1_y_zzz_yzz_0[i] = 2.0 * ta1_y_z_yzz_0[i] * fe_0 - 2.0 * ta1_y_z_yzz_1[i] * fe_0 + 2.0 * ta1_y_zz_yz_0[i] * fe_0 -
                             2.0 * ta1_y_zz_yz_1[i] * fe_0 + ta1_y_zz_yzz_0[i] * pa_z[i] - ta1_y_zz_yzz_1[i] * pc_z[i];

        ta1_y_zzz_zzz_0[i] = 2.0 * ta1_y_z_zzz_0[i] * fe_0 - 2.0 * ta1_y_z_zzz_1[i] * fe_0 + 3.0 * ta1_y_zz_zz_0[i] * fe_0 -
                             3.0 * ta1_y_zz_zz_1[i] * fe_0 + ta1_y_zz_zzz_0[i] * pa_z[i] - ta1_y_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 200-210 components of targeted buffer : FF

    auto ta1_z_xxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 200);

    auto ta1_z_xxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 201);

    auto ta1_z_xxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 202);

    auto ta1_z_xxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 203);

    auto ta1_z_xxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 204);

    auto ta1_z_xxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 205);

    auto ta1_z_xxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 206);

    auto ta1_z_xxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 207);

    auto ta1_z_xxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 208);

    auto ta1_z_xxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 209);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_x_xxx_0,   \
                             ta1_z_x_xxx_1,   \
                             ta1_z_x_xxy_0,   \
                             ta1_z_x_xxy_1,   \
                             ta1_z_x_xxz_0,   \
                             ta1_z_x_xxz_1,   \
                             ta1_z_x_xyy_0,   \
                             ta1_z_x_xyy_1,   \
                             ta1_z_x_xyz_0,   \
                             ta1_z_x_xyz_1,   \
                             ta1_z_x_xzz_0,   \
                             ta1_z_x_xzz_1,   \
                             ta1_z_x_yyy_0,   \
                             ta1_z_x_yyy_1,   \
                             ta1_z_x_yyz_0,   \
                             ta1_z_x_yyz_1,   \
                             ta1_z_x_yzz_0,   \
                             ta1_z_x_yzz_1,   \
                             ta1_z_x_zzz_0,   \
                             ta1_z_x_zzz_1,   \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xxx_0,  \
                             ta1_z_xx_xxx_1,  \
                             ta1_z_xx_xxy_0,  \
                             ta1_z_xx_xxy_1,  \
                             ta1_z_xx_xxz_0,  \
                             ta1_z_xx_xxz_1,  \
                             ta1_z_xx_xy_0,   \
                             ta1_z_xx_xy_1,   \
                             ta1_z_xx_xyy_0,  \
                             ta1_z_xx_xyy_1,  \
                             ta1_z_xx_xyz_0,  \
                             ta1_z_xx_xyz_1,  \
                             ta1_z_xx_xz_0,   \
                             ta1_z_xx_xz_1,   \
                             ta1_z_xx_xzz_0,  \
                             ta1_z_xx_xzz_1,  \
                             ta1_z_xx_yy_0,   \
                             ta1_z_xx_yy_1,   \
                             ta1_z_xx_yyy_0,  \
                             ta1_z_xx_yyy_1,  \
                             ta1_z_xx_yyz_0,  \
                             ta1_z_xx_yyz_1,  \
                             ta1_z_xx_yz_0,   \
                             ta1_z_xx_yz_1,   \
                             ta1_z_xx_yzz_0,  \
                             ta1_z_xx_yzz_1,  \
                             ta1_z_xx_zz_0,   \
                             ta1_z_xx_zz_1,   \
                             ta1_z_xx_zzz_0,  \
                             ta1_z_xx_zzz_1,  \
                             ta1_z_xxx_xxx_0, \
                             ta1_z_xxx_xxy_0, \
                             ta1_z_xxx_xxz_0, \
                             ta1_z_xxx_xyy_0, \
                             ta1_z_xxx_xyz_0, \
                             ta1_z_xxx_xzz_0, \
                             ta1_z_xxx_yyy_0, \
                             ta1_z_xxx_yyz_0, \
                             ta1_z_xxx_yzz_0, \
                             ta1_z_xxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_xxx_0[i] = 2.0 * ta1_z_x_xxx_0[i] * fe_0 - 2.0 * ta1_z_x_xxx_1[i] * fe_0 + 3.0 * ta1_z_xx_xx_0[i] * fe_0 -
                             3.0 * ta1_z_xx_xx_1[i] * fe_0 + ta1_z_xx_xxx_0[i] * pa_x[i] - ta1_z_xx_xxx_1[i] * pc_x[i];

        ta1_z_xxx_xxy_0[i] = 2.0 * ta1_z_x_xxy_0[i] * fe_0 - 2.0 * ta1_z_x_xxy_1[i] * fe_0 + 2.0 * ta1_z_xx_xy_0[i] * fe_0 -
                             2.0 * ta1_z_xx_xy_1[i] * fe_0 + ta1_z_xx_xxy_0[i] * pa_x[i] - ta1_z_xx_xxy_1[i] * pc_x[i];

        ta1_z_xxx_xxz_0[i] = 2.0 * ta1_z_x_xxz_0[i] * fe_0 - 2.0 * ta1_z_x_xxz_1[i] * fe_0 + 2.0 * ta1_z_xx_xz_0[i] * fe_0 -
                             2.0 * ta1_z_xx_xz_1[i] * fe_0 + ta1_z_xx_xxz_0[i] * pa_x[i] - ta1_z_xx_xxz_1[i] * pc_x[i];

        ta1_z_xxx_xyy_0[i] = 2.0 * ta1_z_x_xyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyy_1[i] * fe_0 + ta1_z_xx_yy_0[i] * fe_0 - ta1_z_xx_yy_1[i] * fe_0 +
                             ta1_z_xx_xyy_0[i] * pa_x[i] - ta1_z_xx_xyy_1[i] * pc_x[i];

        ta1_z_xxx_xyz_0[i] = 2.0 * ta1_z_x_xyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyz_1[i] * fe_0 + ta1_z_xx_yz_0[i] * fe_0 - ta1_z_xx_yz_1[i] * fe_0 +
                             ta1_z_xx_xyz_0[i] * pa_x[i] - ta1_z_xx_xyz_1[i] * pc_x[i];

        ta1_z_xxx_xzz_0[i] = 2.0 * ta1_z_x_xzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzz_1[i] * fe_0 + ta1_z_xx_zz_0[i] * fe_0 - ta1_z_xx_zz_1[i] * fe_0 +
                             ta1_z_xx_xzz_0[i] * pa_x[i] - ta1_z_xx_xzz_1[i] * pc_x[i];

        ta1_z_xxx_yyy_0[i] =
            2.0 * ta1_z_x_yyy_0[i] * fe_0 - 2.0 * ta1_z_x_yyy_1[i] * fe_0 + ta1_z_xx_yyy_0[i] * pa_x[i] - ta1_z_xx_yyy_1[i] * pc_x[i];

        ta1_z_xxx_yyz_0[i] =
            2.0 * ta1_z_x_yyz_0[i] * fe_0 - 2.0 * ta1_z_x_yyz_1[i] * fe_0 + ta1_z_xx_yyz_0[i] * pa_x[i] - ta1_z_xx_yyz_1[i] * pc_x[i];

        ta1_z_xxx_yzz_0[i] =
            2.0 * ta1_z_x_yzz_0[i] * fe_0 - 2.0 * ta1_z_x_yzz_1[i] * fe_0 + ta1_z_xx_yzz_0[i] * pa_x[i] - ta1_z_xx_yzz_1[i] * pc_x[i];

        ta1_z_xxx_zzz_0[i] =
            2.0 * ta1_z_x_zzz_0[i] * fe_0 - 2.0 * ta1_z_x_zzz_1[i] * fe_0 + ta1_z_xx_zzz_0[i] * pa_x[i] - ta1_z_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : FF

    auto ta1_z_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 210);

    auto ta1_z_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 211);

    auto ta1_z_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 212);

    auto ta1_z_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 213);

    auto ta1_z_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 214);

    auto ta1_z_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 215);

    auto ta1_z_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 216);

    auto ta1_z_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 217);

    auto ta1_z_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 218);

    auto ta1_z_xxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 219);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xxx_0,  \
                             ta1_z_xx_xxx_1,  \
                             ta1_z_xx_xxy_0,  \
                             ta1_z_xx_xxy_1,  \
                             ta1_z_xx_xxz_0,  \
                             ta1_z_xx_xxz_1,  \
                             ta1_z_xx_xy_0,   \
                             ta1_z_xx_xy_1,   \
                             ta1_z_xx_xyy_0,  \
                             ta1_z_xx_xyy_1,  \
                             ta1_z_xx_xyz_0,  \
                             ta1_z_xx_xyz_1,  \
                             ta1_z_xx_xz_0,   \
                             ta1_z_xx_xz_1,   \
                             ta1_z_xx_xzz_0,  \
                             ta1_z_xx_xzz_1,  \
                             ta1_z_xx_zzz_0,  \
                             ta1_z_xx_zzz_1,  \
                             ta1_z_xxy_xxx_0, \
                             ta1_z_xxy_xxy_0, \
                             ta1_z_xxy_xxz_0, \
                             ta1_z_xxy_xyy_0, \
                             ta1_z_xxy_xyz_0, \
                             ta1_z_xxy_xzz_0, \
                             ta1_z_xxy_yyy_0, \
                             ta1_z_xxy_yyz_0, \
                             ta1_z_xxy_yzz_0, \
                             ta1_z_xxy_zzz_0, \
                             ta1_z_xy_yyy_0,  \
                             ta1_z_xy_yyy_1,  \
                             ta1_z_xy_yyz_0,  \
                             ta1_z_xy_yyz_1,  \
                             ta1_z_xy_yzz_0,  \
                             ta1_z_xy_yzz_1,  \
                             ta1_z_y_yyy_0,   \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyz_0,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yzz_0,   \
                             ta1_z_y_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_xxx_0[i] = ta1_z_xx_xxx_0[i] * pa_y[i] - ta1_z_xx_xxx_1[i] * pc_y[i];

        ta1_z_xxy_xxy_0[i] = ta1_z_xx_xx_0[i] * fe_0 - ta1_z_xx_xx_1[i] * fe_0 + ta1_z_xx_xxy_0[i] * pa_y[i] - ta1_z_xx_xxy_1[i] * pc_y[i];

        ta1_z_xxy_xxz_0[i] = ta1_z_xx_xxz_0[i] * pa_y[i] - ta1_z_xx_xxz_1[i] * pc_y[i];

        ta1_z_xxy_xyy_0[i] =
            2.0 * ta1_z_xx_xy_0[i] * fe_0 - 2.0 * ta1_z_xx_xy_1[i] * fe_0 + ta1_z_xx_xyy_0[i] * pa_y[i] - ta1_z_xx_xyy_1[i] * pc_y[i];

        ta1_z_xxy_xyz_0[i] = ta1_z_xx_xz_0[i] * fe_0 - ta1_z_xx_xz_1[i] * fe_0 + ta1_z_xx_xyz_0[i] * pa_y[i] - ta1_z_xx_xyz_1[i] * pc_y[i];

        ta1_z_xxy_xzz_0[i] = ta1_z_xx_xzz_0[i] * pa_y[i] - ta1_z_xx_xzz_1[i] * pc_y[i];

        ta1_z_xxy_yyy_0[i] = ta1_z_y_yyy_0[i] * fe_0 - ta1_z_y_yyy_1[i] * fe_0 + ta1_z_xy_yyy_0[i] * pa_x[i] - ta1_z_xy_yyy_1[i] * pc_x[i];

        ta1_z_xxy_yyz_0[i] = ta1_z_y_yyz_0[i] * fe_0 - ta1_z_y_yyz_1[i] * fe_0 + ta1_z_xy_yyz_0[i] * pa_x[i] - ta1_z_xy_yyz_1[i] * pc_x[i];

        ta1_z_xxy_yzz_0[i] = ta1_z_y_yzz_0[i] * fe_0 - ta1_z_y_yzz_1[i] * fe_0 + ta1_z_xy_yzz_0[i] * pa_x[i] - ta1_z_xy_yzz_1[i] * pc_x[i];

        ta1_z_xxy_zzz_0[i] = ta1_z_xx_zzz_0[i] * pa_y[i] - ta1_z_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 220-230 components of targeted buffer : FF

    auto ta1_z_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 220);

    auto ta1_z_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 221);

    auto ta1_z_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 222);

    auto ta1_z_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 223);

    auto ta1_z_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 224);

    auto ta1_z_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 225);

    auto ta1_z_xxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 226);

    auto ta1_z_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 227);

    auto ta1_z_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 228);

    auto ta1_z_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 229);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xx_xx_0,   \
                             ta1_z_xx_xx_1,   \
                             ta1_z_xx_xxx_0,  \
                             ta1_z_xx_xxx_1,  \
                             ta1_z_xx_xxy_0,  \
                             ta1_z_xx_xxy_1,  \
                             ta1_z_xx_xxz_0,  \
                             ta1_z_xx_xxz_1,  \
                             ta1_z_xx_xy_0,   \
                             ta1_z_xx_xy_1,   \
                             ta1_z_xx_xyy_0,  \
                             ta1_z_xx_xyy_1,  \
                             ta1_z_xx_xyz_0,  \
                             ta1_z_xx_xyz_1,  \
                             ta1_z_xx_xz_0,   \
                             ta1_z_xx_xz_1,   \
                             ta1_z_xx_xzz_0,  \
                             ta1_z_xx_xzz_1,  \
                             ta1_z_xx_yyy_0,  \
                             ta1_z_xx_yyy_1,  \
                             ta1_z_xxz_xxx_0, \
                             ta1_z_xxz_xxy_0, \
                             ta1_z_xxz_xxz_0, \
                             ta1_z_xxz_xyy_0, \
                             ta1_z_xxz_xyz_0, \
                             ta1_z_xxz_xzz_0, \
                             ta1_z_xxz_yyy_0, \
                             ta1_z_xxz_yyz_0, \
                             ta1_z_xxz_yzz_0, \
                             ta1_z_xxz_zzz_0, \
                             ta1_z_xz_yyz_0,  \
                             ta1_z_xz_yyz_1,  \
                             ta1_z_xz_yzz_0,  \
                             ta1_z_xz_yzz_1,  \
                             ta1_z_xz_zzz_0,  \
                             ta1_z_xz_zzz_1,  \
                             ta1_z_z_yyz_0,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yzz_0,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta_xx_xxx_1,     \
                             ta_xx_xxy_1,     \
                             ta_xx_xxz_1,     \
                             ta_xx_xyy_1,     \
                             ta_xx_xyz_1,     \
                             ta_xx_xzz_1,     \
                             ta_xx_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_xxx_0[i] = ta_xx_xxx_1[i] + ta1_z_xx_xxx_0[i] * pa_z[i] - ta1_z_xx_xxx_1[i] * pc_z[i];

        ta1_z_xxz_xxy_0[i] = ta_xx_xxy_1[i] + ta1_z_xx_xxy_0[i] * pa_z[i] - ta1_z_xx_xxy_1[i] * pc_z[i];

        ta1_z_xxz_xxz_0[i] =
            ta1_z_xx_xx_0[i] * fe_0 - ta1_z_xx_xx_1[i] * fe_0 + ta_xx_xxz_1[i] + ta1_z_xx_xxz_0[i] * pa_z[i] - ta1_z_xx_xxz_1[i] * pc_z[i];

        ta1_z_xxz_xyy_0[i] = ta_xx_xyy_1[i] + ta1_z_xx_xyy_0[i] * pa_z[i] - ta1_z_xx_xyy_1[i] * pc_z[i];

        ta1_z_xxz_xyz_0[i] =
            ta1_z_xx_xy_0[i] * fe_0 - ta1_z_xx_xy_1[i] * fe_0 + ta_xx_xyz_1[i] + ta1_z_xx_xyz_0[i] * pa_z[i] - ta1_z_xx_xyz_1[i] * pc_z[i];

        ta1_z_xxz_xzz_0[i] = 2.0 * ta1_z_xx_xz_0[i] * fe_0 - 2.0 * ta1_z_xx_xz_1[i] * fe_0 + ta_xx_xzz_1[i] + ta1_z_xx_xzz_0[i] * pa_z[i] -
                             ta1_z_xx_xzz_1[i] * pc_z[i];

        ta1_z_xxz_yyy_0[i] = ta_xx_yyy_1[i] + ta1_z_xx_yyy_0[i] * pa_z[i] - ta1_z_xx_yyy_1[i] * pc_z[i];

        ta1_z_xxz_yyz_0[i] = ta1_z_z_yyz_0[i] * fe_0 - ta1_z_z_yyz_1[i] * fe_0 + ta1_z_xz_yyz_0[i] * pa_x[i] - ta1_z_xz_yyz_1[i] * pc_x[i];

        ta1_z_xxz_yzz_0[i] = ta1_z_z_yzz_0[i] * fe_0 - ta1_z_z_yzz_1[i] * fe_0 + ta1_z_xz_yzz_0[i] * pa_x[i] - ta1_z_xz_yzz_1[i] * pc_x[i];

        ta1_z_xxz_zzz_0[i] = ta1_z_z_zzz_0[i] * fe_0 - ta1_z_z_zzz_1[i] * fe_0 + ta1_z_xz_zzz_0[i] * pa_x[i] - ta1_z_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 230-240 components of targeted buffer : FF

    auto ta1_z_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 230);

    auto ta1_z_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 231);

    auto ta1_z_xyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 232);

    auto ta1_z_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 233);

    auto ta1_z_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 234);

    auto ta1_z_xyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 235);

    auto ta1_z_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 236);

    auto ta1_z_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 237);

    auto ta1_z_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 238);

    auto ta1_z_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 239);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xyy_xxx_0, \
                             ta1_z_xyy_xxy_0, \
                             ta1_z_xyy_xxz_0, \
                             ta1_z_xyy_xyy_0, \
                             ta1_z_xyy_xyz_0, \
                             ta1_z_xyy_xzz_0, \
                             ta1_z_xyy_yyy_0, \
                             ta1_z_xyy_yyz_0, \
                             ta1_z_xyy_yzz_0, \
                             ta1_z_xyy_zzz_0, \
                             ta1_z_yy_xx_0,   \
                             ta1_z_yy_xx_1,   \
                             ta1_z_yy_xxx_0,  \
                             ta1_z_yy_xxx_1,  \
                             ta1_z_yy_xxy_0,  \
                             ta1_z_yy_xxy_1,  \
                             ta1_z_yy_xxz_0,  \
                             ta1_z_yy_xxz_1,  \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_xyy_0,  \
                             ta1_z_yy_xyy_1,  \
                             ta1_z_yy_xyz_0,  \
                             ta1_z_yy_xyz_1,  \
                             ta1_z_yy_xz_0,   \
                             ta1_z_yy_xz_1,   \
                             ta1_z_yy_xzz_0,  \
                             ta1_z_yy_xzz_1,  \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yy_yyy_0,  \
                             ta1_z_yy_yyy_1,  \
                             ta1_z_yy_yyz_0,  \
                             ta1_z_yy_yyz_1,  \
                             ta1_z_yy_yz_0,   \
                             ta1_z_yy_yz_1,   \
                             ta1_z_yy_yzz_0,  \
                             ta1_z_yy_yzz_1,  \
                             ta1_z_yy_zz_0,   \
                             ta1_z_yy_zz_1,   \
                             ta1_z_yy_zzz_0,  \
                             ta1_z_yy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_xxx_0[i] =
            3.0 * ta1_z_yy_xx_0[i] * fe_0 - 3.0 * ta1_z_yy_xx_1[i] * fe_0 + ta1_z_yy_xxx_0[i] * pa_x[i] - ta1_z_yy_xxx_1[i] * pc_x[i];

        ta1_z_xyy_xxy_0[i] =
            2.0 * ta1_z_yy_xy_0[i] * fe_0 - 2.0 * ta1_z_yy_xy_1[i] * fe_0 + ta1_z_yy_xxy_0[i] * pa_x[i] - ta1_z_yy_xxy_1[i] * pc_x[i];

        ta1_z_xyy_xxz_0[i] =
            2.0 * ta1_z_yy_xz_0[i] * fe_0 - 2.0 * ta1_z_yy_xz_1[i] * fe_0 + ta1_z_yy_xxz_0[i] * pa_x[i] - ta1_z_yy_xxz_1[i] * pc_x[i];

        ta1_z_xyy_xyy_0[i] = ta1_z_yy_yy_0[i] * fe_0 - ta1_z_yy_yy_1[i] * fe_0 + ta1_z_yy_xyy_0[i] * pa_x[i] - ta1_z_yy_xyy_1[i] * pc_x[i];

        ta1_z_xyy_xyz_0[i] = ta1_z_yy_yz_0[i] * fe_0 - ta1_z_yy_yz_1[i] * fe_0 + ta1_z_yy_xyz_0[i] * pa_x[i] - ta1_z_yy_xyz_1[i] * pc_x[i];

        ta1_z_xyy_xzz_0[i] = ta1_z_yy_zz_0[i] * fe_0 - ta1_z_yy_zz_1[i] * fe_0 + ta1_z_yy_xzz_0[i] * pa_x[i] - ta1_z_yy_xzz_1[i] * pc_x[i];

        ta1_z_xyy_yyy_0[i] = ta1_z_yy_yyy_0[i] * pa_x[i] - ta1_z_yy_yyy_1[i] * pc_x[i];

        ta1_z_xyy_yyz_0[i] = ta1_z_yy_yyz_0[i] * pa_x[i] - ta1_z_yy_yyz_1[i] * pc_x[i];

        ta1_z_xyy_yzz_0[i] = ta1_z_yy_yzz_0[i] * pa_x[i] - ta1_z_yy_yzz_1[i] * pc_x[i];

        ta1_z_xyy_zzz_0[i] = ta1_z_yy_zzz_0[i] * pa_x[i] - ta1_z_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 240-250 components of targeted buffer : FF

    auto ta1_z_xyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 240);

    auto ta1_z_xyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 241);

    auto ta1_z_xyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 242);

    auto ta1_z_xyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 243);

    auto ta1_z_xyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 244);

    auto ta1_z_xyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 245);

    auto ta1_z_xyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 246);

    auto ta1_z_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 247);

    auto ta1_z_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 248);

    auto ta1_z_xyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 249);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_xy_xxy_0,  \
                             ta1_z_xy_xxy_1,  \
                             ta1_z_xy_xyy_0,  \
                             ta1_z_xy_xyy_1,  \
                             ta1_z_xyz_xxx_0, \
                             ta1_z_xyz_xxy_0, \
                             ta1_z_xyz_xxz_0, \
                             ta1_z_xyz_xyy_0, \
                             ta1_z_xyz_xyz_0, \
                             ta1_z_xyz_xzz_0, \
                             ta1_z_xyz_yyy_0, \
                             ta1_z_xyz_yyz_0, \
                             ta1_z_xyz_yzz_0, \
                             ta1_z_xyz_zzz_0, \
                             ta1_z_xz_xxx_0,  \
                             ta1_z_xz_xxx_1,  \
                             ta1_z_xz_xxz_0,  \
                             ta1_z_xz_xxz_1,  \
                             ta1_z_xz_xzz_0,  \
                             ta1_z_xz_xzz_1,  \
                             ta1_z_yz_xyz_0,  \
                             ta1_z_yz_xyz_1,  \
                             ta1_z_yz_yyy_0,  \
                             ta1_z_yz_yyy_1,  \
                             ta1_z_yz_yyz_0,  \
                             ta1_z_yz_yyz_1,  \
                             ta1_z_yz_yz_0,   \
                             ta1_z_yz_yz_1,   \
                             ta1_z_yz_yzz_0,  \
                             ta1_z_yz_yzz_1,  \
                             ta1_z_yz_zzz_0,  \
                             ta1_z_yz_zzz_1,  \
                             ta_xy_xxy_1,     \
                             ta_xy_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyz_xxx_0[i] = ta1_z_xz_xxx_0[i] * pa_y[i] - ta1_z_xz_xxx_1[i] * pc_y[i];

        ta1_z_xyz_xxy_0[i] = ta_xy_xxy_1[i] + ta1_z_xy_xxy_0[i] * pa_z[i] - ta1_z_xy_xxy_1[i] * pc_z[i];

        ta1_z_xyz_xxz_0[i] = ta1_z_xz_xxz_0[i] * pa_y[i] - ta1_z_xz_xxz_1[i] * pc_y[i];

        ta1_z_xyz_xyy_0[i] = ta_xy_xyy_1[i] + ta1_z_xy_xyy_0[i] * pa_z[i] - ta1_z_xy_xyy_1[i] * pc_z[i];

        ta1_z_xyz_xyz_0[i] = ta1_z_yz_yz_0[i] * fe_0 - ta1_z_yz_yz_1[i] * fe_0 + ta1_z_yz_xyz_0[i] * pa_x[i] - ta1_z_yz_xyz_1[i] * pc_x[i];

        ta1_z_xyz_xzz_0[i] = ta1_z_xz_xzz_0[i] * pa_y[i] - ta1_z_xz_xzz_1[i] * pc_y[i];

        ta1_z_xyz_yyy_0[i] = ta1_z_yz_yyy_0[i] * pa_x[i] - ta1_z_yz_yyy_1[i] * pc_x[i];

        ta1_z_xyz_yyz_0[i] = ta1_z_yz_yyz_0[i] * pa_x[i] - ta1_z_yz_yyz_1[i] * pc_x[i];

        ta1_z_xyz_yzz_0[i] = ta1_z_yz_yzz_0[i] * pa_x[i] - ta1_z_yz_yzz_1[i] * pc_x[i];

        ta1_z_xyz_zzz_0[i] = ta1_z_yz_zzz_0[i] * pa_x[i] - ta1_z_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 250-260 components of targeted buffer : FF

    auto ta1_z_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 250);

    auto ta1_z_xzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 251);

    auto ta1_z_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 252);

    auto ta1_z_xzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 253);

    auto ta1_z_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 254);

    auto ta1_z_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 255);

    auto ta1_z_xzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 256);

    auto ta1_z_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 257);

    auto ta1_z_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 258);

    auto ta1_z_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 259);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xzz_xxx_0, \
                             ta1_z_xzz_xxy_0, \
                             ta1_z_xzz_xxz_0, \
                             ta1_z_xzz_xyy_0, \
                             ta1_z_xzz_xyz_0, \
                             ta1_z_xzz_xzz_0, \
                             ta1_z_xzz_yyy_0, \
                             ta1_z_xzz_yyz_0, \
                             ta1_z_xzz_yzz_0, \
                             ta1_z_xzz_zzz_0, \
                             ta1_z_zz_xx_0,   \
                             ta1_z_zz_xx_1,   \
                             ta1_z_zz_xxx_0,  \
                             ta1_z_zz_xxx_1,  \
                             ta1_z_zz_xxy_0,  \
                             ta1_z_zz_xxy_1,  \
                             ta1_z_zz_xxz_0,  \
                             ta1_z_zz_xxz_1,  \
                             ta1_z_zz_xy_0,   \
                             ta1_z_zz_xy_1,   \
                             ta1_z_zz_xyy_0,  \
                             ta1_z_zz_xyy_1,  \
                             ta1_z_zz_xyz_0,  \
                             ta1_z_zz_xyz_1,  \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_xzz_0,  \
                             ta1_z_zz_xzz_1,  \
                             ta1_z_zz_yy_0,   \
                             ta1_z_zz_yy_1,   \
                             ta1_z_zz_yyy_0,  \
                             ta1_z_zz_yyy_1,  \
                             ta1_z_zz_yyz_0,  \
                             ta1_z_zz_yyz_1,  \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_yzz_0,  \
                             ta1_z_zz_yzz_1,  \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta1_z_zz_zzz_0,  \
                             ta1_z_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_xxx_0[i] =
            3.0 * ta1_z_zz_xx_0[i] * fe_0 - 3.0 * ta1_z_zz_xx_1[i] * fe_0 + ta1_z_zz_xxx_0[i] * pa_x[i] - ta1_z_zz_xxx_1[i] * pc_x[i];

        ta1_z_xzz_xxy_0[i] =
            2.0 * ta1_z_zz_xy_0[i] * fe_0 - 2.0 * ta1_z_zz_xy_1[i] * fe_0 + ta1_z_zz_xxy_0[i] * pa_x[i] - ta1_z_zz_xxy_1[i] * pc_x[i];

        ta1_z_xzz_xxz_0[i] =
            2.0 * ta1_z_zz_xz_0[i] * fe_0 - 2.0 * ta1_z_zz_xz_1[i] * fe_0 + ta1_z_zz_xxz_0[i] * pa_x[i] - ta1_z_zz_xxz_1[i] * pc_x[i];

        ta1_z_xzz_xyy_0[i] = ta1_z_zz_yy_0[i] * fe_0 - ta1_z_zz_yy_1[i] * fe_0 + ta1_z_zz_xyy_0[i] * pa_x[i] - ta1_z_zz_xyy_1[i] * pc_x[i];

        ta1_z_xzz_xyz_0[i] = ta1_z_zz_yz_0[i] * fe_0 - ta1_z_zz_yz_1[i] * fe_0 + ta1_z_zz_xyz_0[i] * pa_x[i] - ta1_z_zz_xyz_1[i] * pc_x[i];

        ta1_z_xzz_xzz_0[i] = ta1_z_zz_zz_0[i] * fe_0 - ta1_z_zz_zz_1[i] * fe_0 + ta1_z_zz_xzz_0[i] * pa_x[i] - ta1_z_zz_xzz_1[i] * pc_x[i];

        ta1_z_xzz_yyy_0[i] = ta1_z_zz_yyy_0[i] * pa_x[i] - ta1_z_zz_yyy_1[i] * pc_x[i];

        ta1_z_xzz_yyz_0[i] = ta1_z_zz_yyz_0[i] * pa_x[i] - ta1_z_zz_yyz_1[i] * pc_x[i];

        ta1_z_xzz_yzz_0[i] = ta1_z_zz_yzz_0[i] * pa_x[i] - ta1_z_zz_yzz_1[i] * pc_x[i];

        ta1_z_xzz_zzz_0[i] = ta1_z_zz_zzz_0[i] * pa_x[i] - ta1_z_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 260-270 components of targeted buffer : FF

    auto ta1_z_yyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 260);

    auto ta1_z_yyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 261);

    auto ta1_z_yyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 262);

    auto ta1_z_yyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 263);

    auto ta1_z_yyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 264);

    auto ta1_z_yyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 265);

    auto ta1_z_yyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 266);

    auto ta1_z_yyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 267);

    auto ta1_z_yyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 268);

    auto ta1_z_yyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 269);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_y_xxx_0,   \
                             ta1_z_y_xxx_1,   \
                             ta1_z_y_xxy_0,   \
                             ta1_z_y_xxy_1,   \
                             ta1_z_y_xxz_0,   \
                             ta1_z_y_xxz_1,   \
                             ta1_z_y_xyy_0,   \
                             ta1_z_y_xyy_1,   \
                             ta1_z_y_xyz_0,   \
                             ta1_z_y_xyz_1,   \
                             ta1_z_y_xzz_0,   \
                             ta1_z_y_xzz_1,   \
                             ta1_z_y_yyy_0,   \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyz_0,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yzz_0,   \
                             ta1_z_y_yzz_1,   \
                             ta1_z_y_zzz_0,   \
                             ta1_z_y_zzz_1,   \
                             ta1_z_yy_xx_0,   \
                             ta1_z_yy_xx_1,   \
                             ta1_z_yy_xxx_0,  \
                             ta1_z_yy_xxx_1,  \
                             ta1_z_yy_xxy_0,  \
                             ta1_z_yy_xxy_1,  \
                             ta1_z_yy_xxz_0,  \
                             ta1_z_yy_xxz_1,  \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_xyy_0,  \
                             ta1_z_yy_xyy_1,  \
                             ta1_z_yy_xyz_0,  \
                             ta1_z_yy_xyz_1,  \
                             ta1_z_yy_xz_0,   \
                             ta1_z_yy_xz_1,   \
                             ta1_z_yy_xzz_0,  \
                             ta1_z_yy_xzz_1,  \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yy_yyy_0,  \
                             ta1_z_yy_yyy_1,  \
                             ta1_z_yy_yyz_0,  \
                             ta1_z_yy_yyz_1,  \
                             ta1_z_yy_yz_0,   \
                             ta1_z_yy_yz_1,   \
                             ta1_z_yy_yzz_0,  \
                             ta1_z_yy_yzz_1,  \
                             ta1_z_yy_zz_0,   \
                             ta1_z_yy_zz_1,   \
                             ta1_z_yy_zzz_0,  \
                             ta1_z_yy_zzz_1,  \
                             ta1_z_yyy_xxx_0, \
                             ta1_z_yyy_xxy_0, \
                             ta1_z_yyy_xxz_0, \
                             ta1_z_yyy_xyy_0, \
                             ta1_z_yyy_xyz_0, \
                             ta1_z_yyy_xzz_0, \
                             ta1_z_yyy_yyy_0, \
                             ta1_z_yyy_yyz_0, \
                             ta1_z_yyy_yzz_0, \
                             ta1_z_yyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_xxx_0[i] =
            2.0 * ta1_z_y_xxx_0[i] * fe_0 - 2.0 * ta1_z_y_xxx_1[i] * fe_0 + ta1_z_yy_xxx_0[i] * pa_y[i] - ta1_z_yy_xxx_1[i] * pc_y[i];

        ta1_z_yyy_xxy_0[i] = 2.0 * ta1_z_y_xxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxy_1[i] * fe_0 + ta1_z_yy_xx_0[i] * fe_0 - ta1_z_yy_xx_1[i] * fe_0 +
                             ta1_z_yy_xxy_0[i] * pa_y[i] - ta1_z_yy_xxy_1[i] * pc_y[i];

        ta1_z_yyy_xxz_0[i] =
            2.0 * ta1_z_y_xxz_0[i] * fe_0 - 2.0 * ta1_z_y_xxz_1[i] * fe_0 + ta1_z_yy_xxz_0[i] * pa_y[i] - ta1_z_yy_xxz_1[i] * pc_y[i];

        ta1_z_yyy_xyy_0[i] = 2.0 * ta1_z_y_xyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyy_1[i] * fe_0 + 2.0 * ta1_z_yy_xy_0[i] * fe_0 -
                             2.0 * ta1_z_yy_xy_1[i] * fe_0 + ta1_z_yy_xyy_0[i] * pa_y[i] - ta1_z_yy_xyy_1[i] * pc_y[i];

        ta1_z_yyy_xyz_0[i] = 2.0 * ta1_z_y_xyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyz_1[i] * fe_0 + ta1_z_yy_xz_0[i] * fe_0 - ta1_z_yy_xz_1[i] * fe_0 +
                             ta1_z_yy_xyz_0[i] * pa_y[i] - ta1_z_yy_xyz_1[i] * pc_y[i];

        ta1_z_yyy_xzz_0[i] =
            2.0 * ta1_z_y_xzz_0[i] * fe_0 - 2.0 * ta1_z_y_xzz_1[i] * fe_0 + ta1_z_yy_xzz_0[i] * pa_y[i] - ta1_z_yy_xzz_1[i] * pc_y[i];

        ta1_z_yyy_yyy_0[i] = 2.0 * ta1_z_y_yyy_0[i] * fe_0 - 2.0 * ta1_z_y_yyy_1[i] * fe_0 + 3.0 * ta1_z_yy_yy_0[i] * fe_0 -
                             3.0 * ta1_z_yy_yy_1[i] * fe_0 + ta1_z_yy_yyy_0[i] * pa_y[i] - ta1_z_yy_yyy_1[i] * pc_y[i];

        ta1_z_yyy_yyz_0[i] = 2.0 * ta1_z_y_yyz_0[i] * fe_0 - 2.0 * ta1_z_y_yyz_1[i] * fe_0 + 2.0 * ta1_z_yy_yz_0[i] * fe_0 -
                             2.0 * ta1_z_yy_yz_1[i] * fe_0 + ta1_z_yy_yyz_0[i] * pa_y[i] - ta1_z_yy_yyz_1[i] * pc_y[i];

        ta1_z_yyy_yzz_0[i] = 2.0 * ta1_z_y_yzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzz_1[i] * fe_0 + ta1_z_yy_zz_0[i] * fe_0 - ta1_z_yy_zz_1[i] * fe_0 +
                             ta1_z_yy_yzz_0[i] * pa_y[i] - ta1_z_yy_yzz_1[i] * pc_y[i];

        ta1_z_yyy_zzz_0[i] =
            2.0 * ta1_z_y_zzz_0[i] * fe_0 - 2.0 * ta1_z_y_zzz_1[i] * fe_0 + ta1_z_yy_zzz_0[i] * pa_y[i] - ta1_z_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 270-280 components of targeted buffer : FF

    auto ta1_z_yyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 270);

    auto ta1_z_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 271);

    auto ta1_z_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 272);

    auto ta1_z_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 273);

    auto ta1_z_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 274);

    auto ta1_z_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 275);

    auto ta1_z_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 276);

    auto ta1_z_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 277);

    auto ta1_z_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 278);

    auto ta1_z_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 279);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yy_xxx_0,  \
                             ta1_z_yy_xxx_1,  \
                             ta1_z_yy_xxy_0,  \
                             ta1_z_yy_xxy_1,  \
                             ta1_z_yy_xy_0,   \
                             ta1_z_yy_xy_1,   \
                             ta1_z_yy_xyy_0,  \
                             ta1_z_yy_xyy_1,  \
                             ta1_z_yy_xyz_0,  \
                             ta1_z_yy_xyz_1,  \
                             ta1_z_yy_yy_0,   \
                             ta1_z_yy_yy_1,   \
                             ta1_z_yy_yyy_0,  \
                             ta1_z_yy_yyy_1,  \
                             ta1_z_yy_yyz_0,  \
                             ta1_z_yy_yyz_1,  \
                             ta1_z_yy_yz_0,   \
                             ta1_z_yy_yz_1,   \
                             ta1_z_yy_yzz_0,  \
                             ta1_z_yy_yzz_1,  \
                             ta1_z_yyz_xxx_0, \
                             ta1_z_yyz_xxy_0, \
                             ta1_z_yyz_xxz_0, \
                             ta1_z_yyz_xyy_0, \
                             ta1_z_yyz_xyz_0, \
                             ta1_z_yyz_xzz_0, \
                             ta1_z_yyz_yyy_0, \
                             ta1_z_yyz_yyz_0, \
                             ta1_z_yyz_yzz_0, \
                             ta1_z_yyz_zzz_0, \
                             ta1_z_yz_xxz_0,  \
                             ta1_z_yz_xxz_1,  \
                             ta1_z_yz_xzz_0,  \
                             ta1_z_yz_xzz_1,  \
                             ta1_z_yz_zzz_0,  \
                             ta1_z_yz_zzz_1,  \
                             ta1_z_z_xxz_0,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xzz_0,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta_yy_xxx_1,     \
                             ta_yy_xxy_1,     \
                             ta_yy_xyy_1,     \
                             ta_yy_xyz_1,     \
                             ta_yy_yyy_1,     \
                             ta_yy_yyz_1,     \
                             ta_yy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_xxx_0[i] = ta_yy_xxx_1[i] + ta1_z_yy_xxx_0[i] * pa_z[i] - ta1_z_yy_xxx_1[i] * pc_z[i];

        ta1_z_yyz_xxy_0[i] = ta_yy_xxy_1[i] + ta1_z_yy_xxy_0[i] * pa_z[i] - ta1_z_yy_xxy_1[i] * pc_z[i];

        ta1_z_yyz_xxz_0[i] = ta1_z_z_xxz_0[i] * fe_0 - ta1_z_z_xxz_1[i] * fe_0 + ta1_z_yz_xxz_0[i] * pa_y[i] - ta1_z_yz_xxz_1[i] * pc_y[i];

        ta1_z_yyz_xyy_0[i] = ta_yy_xyy_1[i] + ta1_z_yy_xyy_0[i] * pa_z[i] - ta1_z_yy_xyy_1[i] * pc_z[i];

        ta1_z_yyz_xyz_0[i] =
            ta1_z_yy_xy_0[i] * fe_0 - ta1_z_yy_xy_1[i] * fe_0 + ta_yy_xyz_1[i] + ta1_z_yy_xyz_0[i] * pa_z[i] - ta1_z_yy_xyz_1[i] * pc_z[i];

        ta1_z_yyz_xzz_0[i] = ta1_z_z_xzz_0[i] * fe_0 - ta1_z_z_xzz_1[i] * fe_0 + ta1_z_yz_xzz_0[i] * pa_y[i] - ta1_z_yz_xzz_1[i] * pc_y[i];

        ta1_z_yyz_yyy_0[i] = ta_yy_yyy_1[i] + ta1_z_yy_yyy_0[i] * pa_z[i] - ta1_z_yy_yyy_1[i] * pc_z[i];

        ta1_z_yyz_yyz_0[i] =
            ta1_z_yy_yy_0[i] * fe_0 - ta1_z_yy_yy_1[i] * fe_0 + ta_yy_yyz_1[i] + ta1_z_yy_yyz_0[i] * pa_z[i] - ta1_z_yy_yyz_1[i] * pc_z[i];

        ta1_z_yyz_yzz_0[i] = 2.0 * ta1_z_yy_yz_0[i] * fe_0 - 2.0 * ta1_z_yy_yz_1[i] * fe_0 + ta_yy_yzz_1[i] + ta1_z_yy_yzz_0[i] * pa_z[i] -
                             ta1_z_yy_yzz_1[i] * pc_z[i];

        ta1_z_yyz_zzz_0[i] = ta1_z_z_zzz_0[i] * fe_0 - ta1_z_z_zzz_1[i] * fe_0 + ta1_z_yz_zzz_0[i] * pa_y[i] - ta1_z_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 280-290 components of targeted buffer : FF

    auto ta1_z_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 280);

    auto ta1_z_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 281);

    auto ta1_z_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 282);

    auto ta1_z_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 283);

    auto ta1_z_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 284);

    auto ta1_z_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 285);

    auto ta1_z_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 286);

    auto ta1_z_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 287);

    auto ta1_z_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 288);

    auto ta1_z_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 289);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yzz_xxx_0, \
                             ta1_z_yzz_xxy_0, \
                             ta1_z_yzz_xxz_0, \
                             ta1_z_yzz_xyy_0, \
                             ta1_z_yzz_xyz_0, \
                             ta1_z_yzz_xzz_0, \
                             ta1_z_yzz_yyy_0, \
                             ta1_z_yzz_yyz_0, \
                             ta1_z_yzz_yzz_0, \
                             ta1_z_yzz_zzz_0, \
                             ta1_z_zz_xx_0,   \
                             ta1_z_zz_xx_1,   \
                             ta1_z_zz_xxx_0,  \
                             ta1_z_zz_xxx_1,  \
                             ta1_z_zz_xxy_0,  \
                             ta1_z_zz_xxy_1,  \
                             ta1_z_zz_xxz_0,  \
                             ta1_z_zz_xxz_1,  \
                             ta1_z_zz_xy_0,   \
                             ta1_z_zz_xy_1,   \
                             ta1_z_zz_xyy_0,  \
                             ta1_z_zz_xyy_1,  \
                             ta1_z_zz_xyz_0,  \
                             ta1_z_zz_xyz_1,  \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_xzz_0,  \
                             ta1_z_zz_xzz_1,  \
                             ta1_z_zz_yy_0,   \
                             ta1_z_zz_yy_1,   \
                             ta1_z_zz_yyy_0,  \
                             ta1_z_zz_yyy_1,  \
                             ta1_z_zz_yyz_0,  \
                             ta1_z_zz_yyz_1,  \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_yzz_0,  \
                             ta1_z_zz_yzz_1,  \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta1_z_zz_zzz_0,  \
                             ta1_z_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_xxx_0[i] = ta1_z_zz_xxx_0[i] * pa_y[i] - ta1_z_zz_xxx_1[i] * pc_y[i];

        ta1_z_yzz_xxy_0[i] = ta1_z_zz_xx_0[i] * fe_0 - ta1_z_zz_xx_1[i] * fe_0 + ta1_z_zz_xxy_0[i] * pa_y[i] - ta1_z_zz_xxy_1[i] * pc_y[i];

        ta1_z_yzz_xxz_0[i] = ta1_z_zz_xxz_0[i] * pa_y[i] - ta1_z_zz_xxz_1[i] * pc_y[i];

        ta1_z_yzz_xyy_0[i] =
            2.0 * ta1_z_zz_xy_0[i] * fe_0 - 2.0 * ta1_z_zz_xy_1[i] * fe_0 + ta1_z_zz_xyy_0[i] * pa_y[i] - ta1_z_zz_xyy_1[i] * pc_y[i];

        ta1_z_yzz_xyz_0[i] = ta1_z_zz_xz_0[i] * fe_0 - ta1_z_zz_xz_1[i] * fe_0 + ta1_z_zz_xyz_0[i] * pa_y[i] - ta1_z_zz_xyz_1[i] * pc_y[i];

        ta1_z_yzz_xzz_0[i] = ta1_z_zz_xzz_0[i] * pa_y[i] - ta1_z_zz_xzz_1[i] * pc_y[i];

        ta1_z_yzz_yyy_0[i] =
            3.0 * ta1_z_zz_yy_0[i] * fe_0 - 3.0 * ta1_z_zz_yy_1[i] * fe_0 + ta1_z_zz_yyy_0[i] * pa_y[i] - ta1_z_zz_yyy_1[i] * pc_y[i];

        ta1_z_yzz_yyz_0[i] =
            2.0 * ta1_z_zz_yz_0[i] * fe_0 - 2.0 * ta1_z_zz_yz_1[i] * fe_0 + ta1_z_zz_yyz_0[i] * pa_y[i] - ta1_z_zz_yyz_1[i] * pc_y[i];

        ta1_z_yzz_yzz_0[i] = ta1_z_zz_zz_0[i] * fe_0 - ta1_z_zz_zz_1[i] * fe_0 + ta1_z_zz_yzz_0[i] * pa_y[i] - ta1_z_zz_yzz_1[i] * pc_y[i];

        ta1_z_yzz_zzz_0[i] = ta1_z_zz_zzz_0[i] * pa_y[i] - ta1_z_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 290-300 components of targeted buffer : FF

    auto ta1_z_zzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 290);

    auto ta1_z_zzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 291);

    auto ta1_z_zzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 292);

    auto ta1_z_zzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 293);

    auto ta1_z_zzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 294);

    auto ta1_z_zzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 295);

    auto ta1_z_zzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 296);

    auto ta1_z_zzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 297);

    auto ta1_z_zzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 298);

    auto ta1_z_zzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 299);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_z_xxx_0,   \
                             ta1_z_z_xxx_1,   \
                             ta1_z_z_xxy_0,   \
                             ta1_z_z_xxy_1,   \
                             ta1_z_z_xxz_0,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xyy_0,   \
                             ta1_z_z_xyy_1,   \
                             ta1_z_z_xyz_0,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xzz_0,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_yyy_0,   \
                             ta1_z_z_yyy_1,   \
                             ta1_z_z_yyz_0,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yzz_0,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta1_z_zz_xx_0,   \
                             ta1_z_zz_xx_1,   \
                             ta1_z_zz_xxx_0,  \
                             ta1_z_zz_xxx_1,  \
                             ta1_z_zz_xxy_0,  \
                             ta1_z_zz_xxy_1,  \
                             ta1_z_zz_xxz_0,  \
                             ta1_z_zz_xxz_1,  \
                             ta1_z_zz_xy_0,   \
                             ta1_z_zz_xy_1,   \
                             ta1_z_zz_xyy_0,  \
                             ta1_z_zz_xyy_1,  \
                             ta1_z_zz_xyz_0,  \
                             ta1_z_zz_xyz_1,  \
                             ta1_z_zz_xz_0,   \
                             ta1_z_zz_xz_1,   \
                             ta1_z_zz_xzz_0,  \
                             ta1_z_zz_xzz_1,  \
                             ta1_z_zz_yy_0,   \
                             ta1_z_zz_yy_1,   \
                             ta1_z_zz_yyy_0,  \
                             ta1_z_zz_yyy_1,  \
                             ta1_z_zz_yyz_0,  \
                             ta1_z_zz_yyz_1,  \
                             ta1_z_zz_yz_0,   \
                             ta1_z_zz_yz_1,   \
                             ta1_z_zz_yzz_0,  \
                             ta1_z_zz_yzz_1,  \
                             ta1_z_zz_zz_0,   \
                             ta1_z_zz_zz_1,   \
                             ta1_z_zz_zzz_0,  \
                             ta1_z_zz_zzz_1,  \
                             ta1_z_zzz_xxx_0, \
                             ta1_z_zzz_xxy_0, \
                             ta1_z_zzz_xxz_0, \
                             ta1_z_zzz_xyy_0, \
                             ta1_z_zzz_xyz_0, \
                             ta1_z_zzz_xzz_0, \
                             ta1_z_zzz_yyy_0, \
                             ta1_z_zzz_yyz_0, \
                             ta1_z_zzz_yzz_0, \
                             ta1_z_zzz_zzz_0, \
                             ta_zz_xxx_1,     \
                             ta_zz_xxy_1,     \
                             ta_zz_xxz_1,     \
                             ta_zz_xyy_1,     \
                             ta_zz_xyz_1,     \
                             ta_zz_xzz_1,     \
                             ta_zz_yyy_1,     \
                             ta_zz_yyz_1,     \
                             ta_zz_yzz_1,     \
                             ta_zz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_xxx_0[i] = 2.0 * ta1_z_z_xxx_0[i] * fe_0 - 2.0 * ta1_z_z_xxx_1[i] * fe_0 + ta_zz_xxx_1[i] + ta1_z_zz_xxx_0[i] * pa_z[i] -
                             ta1_z_zz_xxx_1[i] * pc_z[i];

        ta1_z_zzz_xxy_0[i] = 2.0 * ta1_z_z_xxy_0[i] * fe_0 - 2.0 * ta1_z_z_xxy_1[i] * fe_0 + ta_zz_xxy_1[i] + ta1_z_zz_xxy_0[i] * pa_z[i] -
                             ta1_z_zz_xxy_1[i] * pc_z[i];

        ta1_z_zzz_xxz_0[i] = 2.0 * ta1_z_z_xxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxz_1[i] * fe_0 + ta1_z_zz_xx_0[i] * fe_0 - ta1_z_zz_xx_1[i] * fe_0 +
                             ta_zz_xxz_1[i] + ta1_z_zz_xxz_0[i] * pa_z[i] - ta1_z_zz_xxz_1[i] * pc_z[i];

        ta1_z_zzz_xyy_0[i] = 2.0 * ta1_z_z_xyy_0[i] * fe_0 - 2.0 * ta1_z_z_xyy_1[i] * fe_0 + ta_zz_xyy_1[i] + ta1_z_zz_xyy_0[i] * pa_z[i] -
                             ta1_z_zz_xyy_1[i] * pc_z[i];

        ta1_z_zzz_xyz_0[i] = 2.0 * ta1_z_z_xyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyz_1[i] * fe_0 + ta1_z_zz_xy_0[i] * fe_0 - ta1_z_zz_xy_1[i] * fe_0 +
                             ta_zz_xyz_1[i] + ta1_z_zz_xyz_0[i] * pa_z[i] - ta1_z_zz_xyz_1[i] * pc_z[i];

        ta1_z_zzz_xzz_0[i] = 2.0 * ta1_z_z_xzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xz_0[i] * fe_0 -
                             2.0 * ta1_z_zz_xz_1[i] * fe_0 + ta_zz_xzz_1[i] + ta1_z_zz_xzz_0[i] * pa_z[i] - ta1_z_zz_xzz_1[i] * pc_z[i];

        ta1_z_zzz_yyy_0[i] = 2.0 * ta1_z_z_yyy_0[i] * fe_0 - 2.0 * ta1_z_z_yyy_1[i] * fe_0 + ta_zz_yyy_1[i] + ta1_z_zz_yyy_0[i] * pa_z[i] -
                             ta1_z_zz_yyy_1[i] * pc_z[i];

        ta1_z_zzz_yyz_0[i] = 2.0 * ta1_z_z_yyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyz_1[i] * fe_0 + ta1_z_zz_yy_0[i] * fe_0 - ta1_z_zz_yy_1[i] * fe_0 +
                             ta_zz_yyz_1[i] + ta1_z_zz_yyz_0[i] * pa_z[i] - ta1_z_zz_yyz_1[i] * pc_z[i];

        ta1_z_zzz_yzz_0[i] = 2.0 * ta1_z_z_yzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzz_1[i] * fe_0 + 2.0 * ta1_z_zz_yz_0[i] * fe_0 -
                             2.0 * ta1_z_zz_yz_1[i] * fe_0 + ta_zz_yzz_1[i] + ta1_z_zz_yzz_0[i] * pa_z[i] - ta1_z_zz_yzz_1[i] * pc_z[i];

        ta1_z_zzz_zzz_0[i] = 2.0 * ta1_z_z_zzz_0[i] * fe_0 - 2.0 * ta1_z_z_zzz_1[i] * fe_0 + 3.0 * ta1_z_zz_zz_0[i] * fe_0 -
                             3.0 * ta1_z_zz_zz_1[i] * fe_0 + ta_zz_zzz_1[i] + ta1_z_zz_zzz_0[i] * pa_z[i] - ta1_z_zz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
