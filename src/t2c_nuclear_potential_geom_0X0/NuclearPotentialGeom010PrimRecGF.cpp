#include "NuclearPotentialGeom010PrimRecGF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gf(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gf,
                                        const size_t              idx_npot_geom_010_0_df,
                                        const size_t              idx_npot_geom_010_1_df,
                                        const size_t              idx_npot_geom_010_0_fd,
                                        const size_t              idx_npot_geom_010_1_fd,
                                        const size_t              idx_npot_1_ff,
                                        const size_t              idx_npot_geom_010_0_ff,
                                        const size_t              idx_npot_geom_010_1_ff,
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

    auto ta1_x_xy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 12);

    auto ta1_x_xy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 15);

    auto ta1_x_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 20);

    auto ta1_x_xz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 21);

    auto ta1_x_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 22);

    auto ta1_x_xz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 23);

    auto ta1_x_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 25);

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

    auto ta1_y_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 76);

    auto ta1_y_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 77);

    auto ta1_y_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 78);

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

    auto ta1_y_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 106);

    auto ta1_y_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 107);

    auto ta1_y_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 108);

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

    auto ta1_z_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 136);

    auto ta1_z_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 137);

    auto ta1_z_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 138);

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

    auto ta1_z_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 165);

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

    auto ta1_x_xy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 12);

    auto ta1_x_xy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 15);

    auto ta1_x_xz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 20);

    auto ta1_x_xz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 21);

    auto ta1_x_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 22);

    auto ta1_x_xz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 23);

    auto ta1_x_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 25);

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

    auto ta1_y_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 76);

    auto ta1_y_xy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 77);

    auto ta1_y_xy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 78);

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

    auto ta1_y_yz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 106);

    auto ta1_y_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 107);

    auto ta1_y_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 108);

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

    auto ta1_z_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 136);

    auto ta1_z_xy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 137);

    auto ta1_z_xy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 138);

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

    auto ta1_z_yz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 165);

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

    auto ta1_x_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 12);

    auto ta1_x_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 13);

    auto ta1_x_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 14);

    auto ta1_x_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 16);

    auto ta1_x_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 17);

    auto ta1_x_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 19);

    auto ta1_x_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 30);

    auto ta1_x_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 31);

    auto ta1_x_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 32);

    auto ta1_x_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 36);

    auto ta1_x_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 37);

    auto ta1_x_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 38);

    auto ta1_x_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 39);

    auto ta1_x_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 40);

    auto ta1_x_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 41);

    auto ta1_x_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 50);

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

    auto ta1_y_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 67);

    auto ta1_y_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 79);

    auto ta1_y_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 81);

    auto ta1_y_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 82);

    auto ta1_y_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 92);

    auto ta1_y_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 94);

    auto ta1_y_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 95);

    auto ta1_y_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 96);

    auto ta1_y_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 97);

    auto ta1_y_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 98);

    auto ta1_y_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 99);

    auto ta1_y_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 100);

    auto ta1_y_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 101);

    auto ta1_y_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 103);

    auto ta1_y_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 104);

    auto ta1_y_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 105);

    auto ta1_y_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 106);

    auto ta1_y_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 107);

    auto ta1_y_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 109);

    auto ta1_y_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 111);

    auto ta1_y_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 112);

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

    auto ta1_z_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 134);

    auto ta1_z_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 139);

    auto ta1_z_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 141);

    auto ta1_z_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 142);

    auto ta1_z_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 152);

    auto ta1_z_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 154);

    auto ta1_z_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 155);

    auto ta1_z_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 156);

    auto ta1_z_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 157);

    auto ta1_z_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 158);

    auto ta1_z_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 159);

    auto ta1_z_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 160);

    auto ta1_z_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 161);

    auto ta1_z_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 164);

    auto ta1_z_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 166);

    auto ta1_z_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 167);

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

    auto ta1_x_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 12);

    auto ta1_x_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 13);

    auto ta1_x_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 14);

    auto ta1_x_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 16);

    auto ta1_x_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 17);

    auto ta1_x_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 19);

    auto ta1_x_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 30);

    auto ta1_x_xzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 31);

    auto ta1_x_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 32);

    auto ta1_x_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 36);

    auto ta1_x_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 37);

    auto ta1_x_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 38);

    auto ta1_x_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 39);

    auto ta1_x_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 40);

    auto ta1_x_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 41);

    auto ta1_x_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 50);

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

    auto ta1_y_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 67);

    auto ta1_y_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 79);

    auto ta1_y_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 81);

    auto ta1_y_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 82);

    auto ta1_y_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 92);

    auto ta1_y_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 94);

    auto ta1_y_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 95);

    auto ta1_y_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 96);

    auto ta1_y_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 97);

    auto ta1_y_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 98);

    auto ta1_y_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 99);

    auto ta1_y_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 100);

    auto ta1_y_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 101);

    auto ta1_y_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 103);

    auto ta1_y_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 104);

    auto ta1_y_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 105);

    auto ta1_y_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 106);

    auto ta1_y_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 107);

    auto ta1_y_yzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 109);

    auto ta1_y_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 111);

    auto ta1_y_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 112);

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

    auto ta1_z_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 134);

    auto ta1_z_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 139);

    auto ta1_z_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 141);

    auto ta1_z_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 142);

    auto ta1_z_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 152);

    auto ta1_z_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 154);

    auto ta1_z_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 155);

    auto ta1_z_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 156);

    auto ta1_z_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 157);

    auto ta1_z_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 158);

    auto ta1_z_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 159);

    auto ta1_z_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 160);

    auto ta1_z_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 161);

    auto ta1_z_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 164);

    auto ta1_z_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 166);

    auto ta1_z_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 167);

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

    // Set up components of auxiliary buffer : FF

    auto ta_xxx_xxx_1 = pbuffer.data(idx_npot_1_ff);

    auto ta_xxx_xxy_1 = pbuffer.data(idx_npot_1_ff + 1);

    auto ta_xxx_xxz_1 = pbuffer.data(idx_npot_1_ff + 2);

    auto ta_xxx_xyy_1 = pbuffer.data(idx_npot_1_ff + 3);

    auto ta_xxx_xyz_1 = pbuffer.data(idx_npot_1_ff + 4);

    auto ta_xxx_xzz_1 = pbuffer.data(idx_npot_1_ff + 5);

    auto ta_xxx_yyy_1 = pbuffer.data(idx_npot_1_ff + 6);

    auto ta_xxx_yyz_1 = pbuffer.data(idx_npot_1_ff + 7);

    auto ta_xxx_yzz_1 = pbuffer.data(idx_npot_1_ff + 8);

    auto ta_xxx_zzz_1 = pbuffer.data(idx_npot_1_ff + 9);

    auto ta_xxy_xxx_1 = pbuffer.data(idx_npot_1_ff + 10);

    auto ta_xxy_xxy_1 = pbuffer.data(idx_npot_1_ff + 11);

    auto ta_xxy_xxz_1 = pbuffer.data(idx_npot_1_ff + 12);

    auto ta_xxy_xyy_1 = pbuffer.data(idx_npot_1_ff + 13);

    auto ta_xxy_xzz_1 = pbuffer.data(idx_npot_1_ff + 15);

    auto ta_xxy_yyy_1 = pbuffer.data(idx_npot_1_ff + 16);

    auto ta_xxz_xxx_1 = pbuffer.data(idx_npot_1_ff + 20);

    auto ta_xxz_xxy_1 = pbuffer.data(idx_npot_1_ff + 21);

    auto ta_xxz_xxz_1 = pbuffer.data(idx_npot_1_ff + 22);

    auto ta_xxz_xyy_1 = pbuffer.data(idx_npot_1_ff + 23);

    auto ta_xxz_xzz_1 = pbuffer.data(idx_npot_1_ff + 25);

    auto ta_xxz_zzz_1 = pbuffer.data(idx_npot_1_ff + 29);

    auto ta_xyy_xxx_1 = pbuffer.data(idx_npot_1_ff + 30);

    auto ta_xyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 31);

    auto ta_xyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 33);

    auto ta_xyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 36);

    auto ta_xyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 37);

    auto ta_xyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 38);

    auto ta_xzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 50);

    auto ta_xzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 52);

    auto ta_xzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 55);

    auto ta_xzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 57);

    auto ta_xzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 58);

    auto ta_xzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 59);

    auto ta_yyy_xxx_1 = pbuffer.data(idx_npot_1_ff + 60);

    auto ta_yyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 61);

    auto ta_yyy_xxz_1 = pbuffer.data(idx_npot_1_ff + 62);

    auto ta_yyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 63);

    auto ta_yyy_xyz_1 = pbuffer.data(idx_npot_1_ff + 64);

    auto ta_yyy_xzz_1 = pbuffer.data(idx_npot_1_ff + 65);

    auto ta_yyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 66);

    auto ta_yyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 67);

    auto ta_yyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 68);

    auto ta_yyy_zzz_1 = pbuffer.data(idx_npot_1_ff + 69);

    auto ta_yyz_xxy_1 = pbuffer.data(idx_npot_1_ff + 71);

    auto ta_yyz_xyy_1 = pbuffer.data(idx_npot_1_ff + 73);

    auto ta_yyz_yyy_1 = pbuffer.data(idx_npot_1_ff + 76);

    auto ta_yyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 77);

    auto ta_yyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 78);

    auto ta_yyz_zzz_1 = pbuffer.data(idx_npot_1_ff + 79);

    auto ta_yzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 82);

    auto ta_yzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 85);

    auto ta_yzz_yyy_1 = pbuffer.data(idx_npot_1_ff + 86);

    auto ta_yzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 87);

    auto ta_yzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 88);

    auto ta_yzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 89);

    auto ta_zzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 90);

    auto ta_zzz_xxy_1 = pbuffer.data(idx_npot_1_ff + 91);

    auto ta_zzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 92);

    auto ta_zzz_xyy_1 = pbuffer.data(idx_npot_1_ff + 93);

    auto ta_zzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 94);

    auto ta_zzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 95);

    auto ta_zzz_yyy_1 = pbuffer.data(idx_npot_1_ff + 96);

    auto ta_zzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 97);

    auto ta_zzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 98);

    auto ta_zzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 99);

    // Set up components of auxiliary buffer : FF

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

    auto ta1_x_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 10);

    auto ta1_x_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 11);

    auto ta1_x_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 12);

    auto ta1_x_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 13);

    auto ta1_x_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 14);

    auto ta1_x_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 15);

    auto ta1_x_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 16);

    auto ta1_x_xxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 19);

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

    auto ta1_x_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 30);

    auto ta1_x_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 31);

    auto ta1_x_xyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 32);

    auto ta1_x_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 33);

    auto ta1_x_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 34);

    auto ta1_x_xyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 35);

    auto ta1_x_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 36);

    auto ta1_x_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 37);

    auto ta1_x_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 38);

    auto ta1_x_xyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 42);

    auto ta1_x_xyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 45);

    auto ta1_x_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 50);

    auto ta1_x_xzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 51);

    auto ta1_x_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 52);

    auto ta1_x_xzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 53);

    auto ta1_x_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 54);

    auto ta1_x_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 55);

    auto ta1_x_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 57);

    auto ta1_x_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 58);

    auto ta1_x_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 59);

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

    auto ta1_x_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 71);

    auto ta1_x_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 72);

    auto ta1_x_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 73);

    auto ta1_x_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 75);

    auto ta1_x_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 76);

    auto ta1_x_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 77);

    auto ta1_x_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 78);

    auto ta1_x_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 79);

    auto ta1_x_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 80);

    auto ta1_x_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 82);

    auto ta1_x_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 84);

    auto ta1_x_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 85);

    auto ta1_x_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 86);

    auto ta1_x_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 87);

    auto ta1_x_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 88);

    auto ta1_x_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 89);

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

    auto ta1_y_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 110);

    auto ta1_y_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 111);

    auto ta1_y_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 112);

    auto ta1_y_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 113);

    auto ta1_y_xxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 114);

    auto ta1_y_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 115);

    auto ta1_y_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 116);

    auto ta1_y_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 117);

    auto ta1_y_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 118);

    auto ta1_y_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 120);

    auto ta1_y_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 121);

    auto ta1_y_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 122);

    auto ta1_y_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 123);

    auto ta1_y_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 125);

    auto ta1_y_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 127);

    auto ta1_y_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 128);

    auto ta1_y_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 129);

    auto ta1_y_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 130);

    auto ta1_y_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 131);

    auto ta1_y_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 133);

    auto ta1_y_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 134);

    auto ta1_y_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 136);

    auto ta1_y_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 137);

    auto ta1_y_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 138);

    auto ta1_y_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 139);

    auto ta1_y_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 147);

    auto ta1_y_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 148);

    auto ta1_y_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 150);

    auto ta1_y_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 152);

    auto ta1_y_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 154);

    auto ta1_y_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 155);

    auto ta1_y_xzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 156);

    auto ta1_y_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 157);

    auto ta1_y_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 158);

    auto ta1_y_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 159);

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

    auto ta1_y_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 181);

    auto ta1_y_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 182);

    auto ta1_y_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 183);

    auto ta1_y_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 184);

    auto ta1_y_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 185);

    auto ta1_y_yzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 186);

    auto ta1_y_yzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 187);

    auto ta1_y_yzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 188);

    auto ta1_y_yzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 189);

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

    auto ta1_z_xxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 210);

    auto ta1_z_xxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 211);

    auto ta1_z_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 212);

    auto ta1_z_xxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 213);

    auto ta1_z_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 215);

    auto ta1_z_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 216);

    auto ta1_z_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 217);

    auto ta1_z_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 218);

    auto ta1_z_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 220);

    auto ta1_z_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 221);

    auto ta1_z_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 222);

    auto ta1_z_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 223);

    auto ta1_z_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 224);

    auto ta1_z_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 225);

    auto ta1_z_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 227);

    auto ta1_z_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 228);

    auto ta1_z_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 229);

    auto ta1_z_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 230);

    auto ta1_z_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 231);

    auto ta1_z_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 233);

    auto ta1_z_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 234);

    auto ta1_z_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 236);

    auto ta1_z_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 237);

    auto ta1_z_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 238);

    auto ta1_z_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 239);

    auto ta1_z_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 247);

    auto ta1_z_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 248);

    auto ta1_z_xzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 250);

    auto ta1_z_xzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 252);

    auto ta1_z_xzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 254);

    auto ta1_z_xzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 255);

    auto ta1_z_xzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 256);

    auto ta1_z_xzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 257);

    auto ta1_z_xzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 258);

    auto ta1_z_xzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 259);

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

    auto ta1_z_yyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 271);

    auto ta1_z_yyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 272);

    auto ta1_z_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 273);

    auto ta1_z_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 274);

    auto ta1_z_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 275);

    auto ta1_z_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 276);

    auto ta1_z_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 277);

    auto ta1_z_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 278);

    auto ta1_z_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 279);

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

    // Set up components of auxiliary buffer : FF

    auto ta1_x_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff);

    auto ta1_x_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 1);

    auto ta1_x_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 2);

    auto ta1_x_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 3);

    auto ta1_x_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 4);

    auto ta1_x_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 5);

    auto ta1_x_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 6);

    auto ta1_x_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 7);

    auto ta1_x_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 8);

    auto ta1_x_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 9);

    auto ta1_x_xxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 10);

    auto ta1_x_xxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 11);

    auto ta1_x_xxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 12);

    auto ta1_x_xxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 13);

    auto ta1_x_xxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 14);

    auto ta1_x_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 15);

    auto ta1_x_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 16);

    auto ta1_x_xxy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 19);

    auto ta1_x_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 20);

    auto ta1_x_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 21);

    auto ta1_x_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 22);

    auto ta1_x_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 23);

    auto ta1_x_xxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 24);

    auto ta1_x_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 25);

    auto ta1_x_xxz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 26);

    auto ta1_x_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 27);

    auto ta1_x_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 28);

    auto ta1_x_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 29);

    auto ta1_x_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 30);

    auto ta1_x_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 31);

    auto ta1_x_xyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 32);

    auto ta1_x_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 33);

    auto ta1_x_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 34);

    auto ta1_x_xyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 35);

    auto ta1_x_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 36);

    auto ta1_x_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 37);

    auto ta1_x_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 38);

    auto ta1_x_xyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 42);

    auto ta1_x_xyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 45);

    auto ta1_x_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 50);

    auto ta1_x_xzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 51);

    auto ta1_x_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 52);

    auto ta1_x_xzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 53);

    auto ta1_x_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 54);

    auto ta1_x_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 55);

    auto ta1_x_xzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 57);

    auto ta1_x_xzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 58);

    auto ta1_x_xzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 59);

    auto ta1_x_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 60);

    auto ta1_x_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 61);

    auto ta1_x_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 62);

    auto ta1_x_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 63);

    auto ta1_x_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 64);

    auto ta1_x_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 65);

    auto ta1_x_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 66);

    auto ta1_x_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 67);

    auto ta1_x_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 68);

    auto ta1_x_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 69);

    auto ta1_x_yyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 71);

    auto ta1_x_yyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 72);

    auto ta1_x_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 73);

    auto ta1_x_yyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 75);

    auto ta1_x_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 76);

    auto ta1_x_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 77);

    auto ta1_x_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 78);

    auto ta1_x_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 79);

    auto ta1_x_yzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 80);

    auto ta1_x_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 82);

    auto ta1_x_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 84);

    auto ta1_x_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 85);

    auto ta1_x_yzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 86);

    auto ta1_x_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 87);

    auto ta1_x_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 88);

    auto ta1_x_yzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 89);

    auto ta1_x_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 90);

    auto ta1_x_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 91);

    auto ta1_x_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 92);

    auto ta1_x_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 93);

    auto ta1_x_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 94);

    auto ta1_x_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 95);

    auto ta1_x_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 96);

    auto ta1_x_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 97);

    auto ta1_x_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 98);

    auto ta1_x_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 99);

    auto ta1_y_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 100);

    auto ta1_y_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 101);

    auto ta1_y_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 102);

    auto ta1_y_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 103);

    auto ta1_y_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 104);

    auto ta1_y_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 105);

    auto ta1_y_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 106);

    auto ta1_y_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 107);

    auto ta1_y_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 108);

    auto ta1_y_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 109);

    auto ta1_y_xxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 110);

    auto ta1_y_xxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 111);

    auto ta1_y_xxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 112);

    auto ta1_y_xxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 113);

    auto ta1_y_xxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 114);

    auto ta1_y_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 115);

    auto ta1_y_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 116);

    auto ta1_y_xxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 117);

    auto ta1_y_xxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 118);

    auto ta1_y_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 120);

    auto ta1_y_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 121);

    auto ta1_y_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 122);

    auto ta1_y_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 123);

    auto ta1_y_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 125);

    auto ta1_y_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 127);

    auto ta1_y_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 128);

    auto ta1_y_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 129);

    auto ta1_y_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 130);

    auto ta1_y_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 131);

    auto ta1_y_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 133);

    auto ta1_y_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 134);

    auto ta1_y_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 136);

    auto ta1_y_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 137);

    auto ta1_y_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 138);

    auto ta1_y_xyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 139);

    auto ta1_y_xyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 147);

    auto ta1_y_xyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 148);

    auto ta1_y_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 150);

    auto ta1_y_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 152);

    auto ta1_y_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 154);

    auto ta1_y_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 155);

    auto ta1_y_xzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 156);

    auto ta1_y_xzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 157);

    auto ta1_y_xzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 158);

    auto ta1_y_xzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 159);

    auto ta1_y_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 160);

    auto ta1_y_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 161);

    auto ta1_y_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 162);

    auto ta1_y_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 163);

    auto ta1_y_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 164);

    auto ta1_y_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 165);

    auto ta1_y_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 166);

    auto ta1_y_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 167);

    auto ta1_y_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 168);

    auto ta1_y_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 169);

    auto ta1_y_yyz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 170);

    auto ta1_y_yyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 171);

    auto ta1_y_yyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 172);

    auto ta1_y_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 173);

    auto ta1_y_yyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 174);

    auto ta1_y_yyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 175);

    auto ta1_y_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 176);

    auto ta1_y_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 177);

    auto ta1_y_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 178);

    auto ta1_y_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 179);

    auto ta1_y_yzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 181);

    auto ta1_y_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 182);

    auto ta1_y_yzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 183);

    auto ta1_y_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 184);

    auto ta1_y_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 185);

    auto ta1_y_yzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 186);

    auto ta1_y_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 187);

    auto ta1_y_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 188);

    auto ta1_y_yzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 189);

    auto ta1_y_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 190);

    auto ta1_y_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 191);

    auto ta1_y_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 192);

    auto ta1_y_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 193);

    auto ta1_y_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 194);

    auto ta1_y_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 195);

    auto ta1_y_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 196);

    auto ta1_y_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 197);

    auto ta1_y_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 198);

    auto ta1_y_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 199);

    auto ta1_z_xxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 200);

    auto ta1_z_xxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 201);

    auto ta1_z_xxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 202);

    auto ta1_z_xxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 203);

    auto ta1_z_xxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 204);

    auto ta1_z_xxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 205);

    auto ta1_z_xxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 206);

    auto ta1_z_xxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 207);

    auto ta1_z_xxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 208);

    auto ta1_z_xxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 209);

    auto ta1_z_xxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 210);

    auto ta1_z_xxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 211);

    auto ta1_z_xxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 212);

    auto ta1_z_xxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 213);

    auto ta1_z_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 215);

    auto ta1_z_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 216);

    auto ta1_z_xxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 217);

    auto ta1_z_xxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 218);

    auto ta1_z_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 220);

    auto ta1_z_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 221);

    auto ta1_z_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 222);

    auto ta1_z_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 223);

    auto ta1_z_xxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 224);

    auto ta1_z_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 225);

    auto ta1_z_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 227);

    auto ta1_z_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 228);

    auto ta1_z_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 229);

    auto ta1_z_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 230);

    auto ta1_z_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 231);

    auto ta1_z_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 233);

    auto ta1_z_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 234);

    auto ta1_z_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 236);

    auto ta1_z_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 237);

    auto ta1_z_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 238);

    auto ta1_z_xyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 239);

    auto ta1_z_xyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 247);

    auto ta1_z_xyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 248);

    auto ta1_z_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 250);

    auto ta1_z_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 252);

    auto ta1_z_xzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 254);

    auto ta1_z_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 255);

    auto ta1_z_xzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 256);

    auto ta1_z_xzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 257);

    auto ta1_z_xzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 258);

    auto ta1_z_xzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 259);

    auto ta1_z_yyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 260);

    auto ta1_z_yyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 261);

    auto ta1_z_yyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 262);

    auto ta1_z_yyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 263);

    auto ta1_z_yyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 264);

    auto ta1_z_yyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 265);

    auto ta1_z_yyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 266);

    auto ta1_z_yyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 267);

    auto ta1_z_yyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 268);

    auto ta1_z_yyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 269);

    auto ta1_z_yyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 271);

    auto ta1_z_yyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 272);

    auto ta1_z_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 273);

    auto ta1_z_yyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 274);

    auto ta1_z_yyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 275);

    auto ta1_z_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 276);

    auto ta1_z_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 277);

    auto ta1_z_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 278);

    auto ta1_z_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 279);

    auto ta1_z_yzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 280);

    auto ta1_z_yzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 281);

    auto ta1_z_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 282);

    auto ta1_z_yzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 283);

    auto ta1_z_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 284);

    auto ta1_z_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 285);

    auto ta1_z_yzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 286);

    auto ta1_z_yzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 287);

    auto ta1_z_yzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 288);

    auto ta1_z_yzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 289);

    auto ta1_z_zzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 290);

    auto ta1_z_zzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 291);

    auto ta1_z_zzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 292);

    auto ta1_z_zzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 293);

    auto ta1_z_zzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 294);

    auto ta1_z_zzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 295);

    auto ta1_z_zzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 296);

    auto ta1_z_zzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 297);

    auto ta1_z_zzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 298);

    auto ta1_z_zzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 299);

    // Set up 0-10 components of targeted buffer : GF

    auto ta1_x_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf);

    auto ta1_x_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 1);

    auto ta1_x_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 2);

    auto ta1_x_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 3);

    auto ta1_x_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 4);

    auto ta1_x_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 5);

    auto ta1_x_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 6);

    auto ta1_x_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 7);

    auto ta1_x_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 8);

    auto ta1_x_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 9);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_yyy_0,   \
                             ta1_x_xx_yyy_1,   \
                             ta1_x_xx_yyz_0,   \
                             ta1_x_xx_yyz_1,   \
                             ta1_x_xx_yzz_0,   \
                             ta1_x_xx_yzz_1,   \
                             ta1_x_xx_zzz_0,   \
                             ta1_x_xx_zzz_1,   \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xxx_0,  \
                             ta1_x_xxx_xxx_1,  \
                             ta1_x_xxx_xxy_0,  \
                             ta1_x_xxx_xxy_1,  \
                             ta1_x_xxx_xxz_0,  \
                             ta1_x_xxx_xxz_1,  \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xyy_0,  \
                             ta1_x_xxx_xyy_1,  \
                             ta1_x_xxx_xyz_0,  \
                             ta1_x_xxx_xyz_1,  \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_xzz_0,  \
                             ta1_x_xxx_xzz_1,  \
                             ta1_x_xxx_yy_0,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxx_yyy_0,  \
                             ta1_x_xxx_yyy_1,  \
                             ta1_x_xxx_yyz_0,  \
                             ta1_x_xxx_yyz_1,  \
                             ta1_x_xxx_yz_0,   \
                             ta1_x_xxx_yz_1,   \
                             ta1_x_xxx_yzz_0,  \
                             ta1_x_xxx_yzz_1,  \
                             ta1_x_xxx_zz_0,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_x_xxx_zzz_0,  \
                             ta1_x_xxx_zzz_1,  \
                             ta1_x_xxxx_xxx_0, \
                             ta1_x_xxxx_xxy_0, \
                             ta1_x_xxxx_xxz_0, \
                             ta1_x_xxxx_xyy_0, \
                             ta1_x_xxxx_xyz_0, \
                             ta1_x_xxxx_xzz_0, \
                             ta1_x_xxxx_yyy_0, \
                             ta1_x_xxxx_yyz_0, \
                             ta1_x_xxxx_yzz_0, \
                             ta1_x_xxxx_zzz_0, \
                             ta_xxx_xxx_1,     \
                             ta_xxx_xxy_1,     \
                             ta_xxx_xxz_1,     \
                             ta_xxx_xyy_1,     \
                             ta_xxx_xyz_1,     \
                             ta_xxx_xzz_1,     \
                             ta_xxx_yyy_1,     \
                             ta_xxx_yyz_1,     \
                             ta_xxx_yzz_1,     \
                             ta_xxx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_xxx_0[i] = 3.0 * ta1_x_xx_xxx_0[i] * fe_0 - 3.0 * ta1_x_xx_xxx_1[i] * fe_0 + 3.0 * ta1_x_xxx_xx_0[i] * fe_0 -
                              3.0 * ta1_x_xxx_xx_1[i] * fe_0 + ta_xxx_xxx_1[i] + ta1_x_xxx_xxx_0[i] * pa_x[i] - ta1_x_xxx_xxx_1[i] * pc_x[i];

        ta1_x_xxxx_xxy_0[i] = 3.0 * ta1_x_xx_xxy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xy_0[i] * fe_0 -
                              2.0 * ta1_x_xxx_xy_1[i] * fe_0 + ta_xxx_xxy_1[i] + ta1_x_xxx_xxy_0[i] * pa_x[i] - ta1_x_xxx_xxy_1[i] * pc_x[i];

        ta1_x_xxxx_xxz_0[i] = 3.0 * ta1_x_xx_xxz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xz_0[i] * fe_0 -
                              2.0 * ta1_x_xxx_xz_1[i] * fe_0 + ta_xxx_xxz_1[i] + ta1_x_xxx_xxz_0[i] * pa_x[i] - ta1_x_xxx_xxz_1[i] * pc_x[i];

        ta1_x_xxxx_xyy_0[i] = 3.0 * ta1_x_xx_xyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xyy_1[i] * fe_0 + ta1_x_xxx_yy_0[i] * fe_0 - ta1_x_xxx_yy_1[i] * fe_0 +
                              ta_xxx_xyy_1[i] + ta1_x_xxx_xyy_0[i] * pa_x[i] - ta1_x_xxx_xyy_1[i] * pc_x[i];

        ta1_x_xxxx_xyz_0[i] = 3.0 * ta1_x_xx_xyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyz_1[i] * fe_0 + ta1_x_xxx_yz_0[i] * fe_0 - ta1_x_xxx_yz_1[i] * fe_0 +
                              ta_xxx_xyz_1[i] + ta1_x_xxx_xyz_0[i] * pa_x[i] - ta1_x_xxx_xyz_1[i] * pc_x[i];

        ta1_x_xxxx_xzz_0[i] = 3.0 * ta1_x_xx_xzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xzz_1[i] * fe_0 + ta1_x_xxx_zz_0[i] * fe_0 - ta1_x_xxx_zz_1[i] * fe_0 +
                              ta_xxx_xzz_1[i] + ta1_x_xxx_xzz_0[i] * pa_x[i] - ta1_x_xxx_xzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyy_0[i] = 3.0 * ta1_x_xx_yyy_0[i] * fe_0 - 3.0 * ta1_x_xx_yyy_1[i] * fe_0 + ta_xxx_yyy_1[i] + ta1_x_xxx_yyy_0[i] * pa_x[i] -
                              ta1_x_xxx_yyy_1[i] * pc_x[i];

        ta1_x_xxxx_yyz_0[i] = 3.0 * ta1_x_xx_yyz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyz_1[i] * fe_0 + ta_xxx_yyz_1[i] + ta1_x_xxx_yyz_0[i] * pa_x[i] -
                              ta1_x_xxx_yyz_1[i] * pc_x[i];

        ta1_x_xxxx_yzz_0[i] = 3.0 * ta1_x_xx_yzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yzz_1[i] * fe_0 + ta_xxx_yzz_1[i] + ta1_x_xxx_yzz_0[i] * pa_x[i] -
                              ta1_x_xxx_yzz_1[i] * pc_x[i];

        ta1_x_xxxx_zzz_0[i] = 3.0 * ta1_x_xx_zzz_0[i] * fe_0 - 3.0 * ta1_x_xx_zzz_1[i] * fe_0 + ta_xxx_zzz_1[i] + ta1_x_xxx_zzz_0[i] * pa_x[i] -
                              ta1_x_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto ta1_x_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 10);

    auto ta1_x_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 11);

    auto ta1_x_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 12);

    auto ta1_x_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 13);

    auto ta1_x_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 14);

    auto ta1_x_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 15);

    auto ta1_x_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 16);

    auto ta1_x_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 17);

    auto ta1_x_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 18);

    auto ta1_x_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 19);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xxx_0,  \
                             ta1_x_xxx_xxx_1,  \
                             ta1_x_xxx_xxy_0,  \
                             ta1_x_xxx_xxy_1,  \
                             ta1_x_xxx_xxz_0,  \
                             ta1_x_xxx_xxz_1,  \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xyy_0,  \
                             ta1_x_xxx_xyy_1,  \
                             ta1_x_xxx_xyz_0,  \
                             ta1_x_xxx_xyz_1,  \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_xzz_0,  \
                             ta1_x_xxx_xzz_1,  \
                             ta1_x_xxx_yy_0,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxx_yyy_0,  \
                             ta1_x_xxx_yyy_1,  \
                             ta1_x_xxx_yyz_0,  \
                             ta1_x_xxx_yyz_1,  \
                             ta1_x_xxx_yz_0,   \
                             ta1_x_xxx_yz_1,   \
                             ta1_x_xxx_yzz_0,  \
                             ta1_x_xxx_yzz_1,  \
                             ta1_x_xxx_zz_0,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_x_xxx_zzz_0,  \
                             ta1_x_xxx_zzz_1,  \
                             ta1_x_xxxy_xxx_0, \
                             ta1_x_xxxy_xxy_0, \
                             ta1_x_xxxy_xxz_0, \
                             ta1_x_xxxy_xyy_0, \
                             ta1_x_xxxy_xyz_0, \
                             ta1_x_xxxy_xzz_0, \
                             ta1_x_xxxy_yyy_0, \
                             ta1_x_xxxy_yyz_0, \
                             ta1_x_xxxy_yzz_0, \
                             ta1_x_xxxy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_xxx_0[i] = ta1_x_xxx_xxx_0[i] * pa_y[i] - ta1_x_xxx_xxx_1[i] * pc_y[i];

        ta1_x_xxxy_xxy_0[i] = ta1_x_xxx_xx_0[i] * fe_0 - ta1_x_xxx_xx_1[i] * fe_0 + ta1_x_xxx_xxy_0[i] * pa_y[i] - ta1_x_xxx_xxy_1[i] * pc_y[i];

        ta1_x_xxxy_xxz_0[i] = ta1_x_xxx_xxz_0[i] * pa_y[i] - ta1_x_xxx_xxz_1[i] * pc_y[i];

        ta1_x_xxxy_xyy_0[i] =
            2.0 * ta1_x_xxx_xy_0[i] * fe_0 - 2.0 * ta1_x_xxx_xy_1[i] * fe_0 + ta1_x_xxx_xyy_0[i] * pa_y[i] - ta1_x_xxx_xyy_1[i] * pc_y[i];

        ta1_x_xxxy_xyz_0[i] = ta1_x_xxx_xz_0[i] * fe_0 - ta1_x_xxx_xz_1[i] * fe_0 + ta1_x_xxx_xyz_0[i] * pa_y[i] - ta1_x_xxx_xyz_1[i] * pc_y[i];

        ta1_x_xxxy_xzz_0[i] = ta1_x_xxx_xzz_0[i] * pa_y[i] - ta1_x_xxx_xzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyy_0[i] =
            3.0 * ta1_x_xxx_yy_0[i] * fe_0 - 3.0 * ta1_x_xxx_yy_1[i] * fe_0 + ta1_x_xxx_yyy_0[i] * pa_y[i] - ta1_x_xxx_yyy_1[i] * pc_y[i];

        ta1_x_xxxy_yyz_0[i] =
            2.0 * ta1_x_xxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yz_1[i] * fe_0 + ta1_x_xxx_yyz_0[i] * pa_y[i] - ta1_x_xxx_yyz_1[i] * pc_y[i];

        ta1_x_xxxy_yzz_0[i] = ta1_x_xxx_zz_0[i] * fe_0 - ta1_x_xxx_zz_1[i] * fe_0 + ta1_x_xxx_yzz_0[i] * pa_y[i] - ta1_x_xxx_yzz_1[i] * pc_y[i];

        ta1_x_xxxy_zzz_0[i] = ta1_x_xxx_zzz_0[i] * pa_y[i] - ta1_x_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto ta1_x_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 20);

    auto ta1_x_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 21);

    auto ta1_x_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 22);

    auto ta1_x_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 23);

    auto ta1_x_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 24);

    auto ta1_x_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 25);

    auto ta1_x_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 26);

    auto ta1_x_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 27);

    auto ta1_x_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 28);

    auto ta1_x_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 29);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xxx_0,  \
                             ta1_x_xxx_xxx_1,  \
                             ta1_x_xxx_xxy_0,  \
                             ta1_x_xxx_xxy_1,  \
                             ta1_x_xxx_xxz_0,  \
                             ta1_x_xxx_xxz_1,  \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xyy_0,  \
                             ta1_x_xxx_xyy_1,  \
                             ta1_x_xxx_xyz_0,  \
                             ta1_x_xxx_xyz_1,  \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_xzz_0,  \
                             ta1_x_xxx_xzz_1,  \
                             ta1_x_xxx_yy_0,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxx_yyy_0,  \
                             ta1_x_xxx_yyy_1,  \
                             ta1_x_xxx_yyz_0,  \
                             ta1_x_xxx_yyz_1,  \
                             ta1_x_xxx_yz_0,   \
                             ta1_x_xxx_yz_1,   \
                             ta1_x_xxx_yzz_0,  \
                             ta1_x_xxx_yzz_1,  \
                             ta1_x_xxx_zz_0,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_x_xxx_zzz_0,  \
                             ta1_x_xxx_zzz_1,  \
                             ta1_x_xxxz_xxx_0, \
                             ta1_x_xxxz_xxy_0, \
                             ta1_x_xxxz_xxz_0, \
                             ta1_x_xxxz_xyy_0, \
                             ta1_x_xxxz_xyz_0, \
                             ta1_x_xxxz_xzz_0, \
                             ta1_x_xxxz_yyy_0, \
                             ta1_x_xxxz_yyz_0, \
                             ta1_x_xxxz_yzz_0, \
                             ta1_x_xxxz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_xxx_0[i] = ta1_x_xxx_xxx_0[i] * pa_z[i] - ta1_x_xxx_xxx_1[i] * pc_z[i];

        ta1_x_xxxz_xxy_0[i] = ta1_x_xxx_xxy_0[i] * pa_z[i] - ta1_x_xxx_xxy_1[i] * pc_z[i];

        ta1_x_xxxz_xxz_0[i] = ta1_x_xxx_xx_0[i] * fe_0 - ta1_x_xxx_xx_1[i] * fe_0 + ta1_x_xxx_xxz_0[i] * pa_z[i] - ta1_x_xxx_xxz_1[i] * pc_z[i];

        ta1_x_xxxz_xyy_0[i] = ta1_x_xxx_xyy_0[i] * pa_z[i] - ta1_x_xxx_xyy_1[i] * pc_z[i];

        ta1_x_xxxz_xyz_0[i] = ta1_x_xxx_xy_0[i] * fe_0 - ta1_x_xxx_xy_1[i] * fe_0 + ta1_x_xxx_xyz_0[i] * pa_z[i] - ta1_x_xxx_xyz_1[i] * pc_z[i];

        ta1_x_xxxz_xzz_0[i] =
            2.0 * ta1_x_xxx_xz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xz_1[i] * fe_0 + ta1_x_xxx_xzz_0[i] * pa_z[i] - ta1_x_xxx_xzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyy_0[i] = ta1_x_xxx_yyy_0[i] * pa_z[i] - ta1_x_xxx_yyy_1[i] * pc_z[i];

        ta1_x_xxxz_yyz_0[i] = ta1_x_xxx_yy_0[i] * fe_0 - ta1_x_xxx_yy_1[i] * fe_0 + ta1_x_xxx_yyz_0[i] * pa_z[i] - ta1_x_xxx_yyz_1[i] * pc_z[i];

        ta1_x_xxxz_yzz_0[i] =
            2.0 * ta1_x_xxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yz_1[i] * fe_0 + ta1_x_xxx_yzz_0[i] * pa_z[i] - ta1_x_xxx_yzz_1[i] * pc_z[i];

        ta1_x_xxxz_zzz_0[i] =
            3.0 * ta1_x_xxx_zz_0[i] * fe_0 - 3.0 * ta1_x_xxx_zz_1[i] * fe_0 + ta1_x_xxx_zzz_0[i] * pa_z[i] - ta1_x_xxx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto ta1_x_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 30);

    auto ta1_x_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 31);

    auto ta1_x_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 32);

    auto ta1_x_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 33);

    auto ta1_x_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 34);

    auto ta1_x_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 35);

    auto ta1_x_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 36);

    auto ta1_x_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 37);

    auto ta1_x_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 38);

    auto ta1_x_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 39);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_zzz_0,   \
                             ta1_x_xx_zzz_1,   \
                             ta1_x_xxy_xx_0,   \
                             ta1_x_xxy_xx_1,   \
                             ta1_x_xxy_xxx_0,  \
                             ta1_x_xxy_xxx_1,  \
                             ta1_x_xxy_xxy_0,  \
                             ta1_x_xxy_xxy_1,  \
                             ta1_x_xxy_xxz_0,  \
                             ta1_x_xxy_xxz_1,  \
                             ta1_x_xxy_xy_0,   \
                             ta1_x_xxy_xy_1,   \
                             ta1_x_xxy_xyy_0,  \
                             ta1_x_xxy_xyy_1,  \
                             ta1_x_xxy_xyz_0,  \
                             ta1_x_xxy_xyz_1,  \
                             ta1_x_xxy_xz_0,   \
                             ta1_x_xxy_xz_1,   \
                             ta1_x_xxy_xzz_0,  \
                             ta1_x_xxy_xzz_1,  \
                             ta1_x_xxy_zzz_0,  \
                             ta1_x_xxy_zzz_1,  \
                             ta1_x_xxyy_xxx_0, \
                             ta1_x_xxyy_xxy_0, \
                             ta1_x_xxyy_xxz_0, \
                             ta1_x_xxyy_xyy_0, \
                             ta1_x_xxyy_xyz_0, \
                             ta1_x_xxyy_xzz_0, \
                             ta1_x_xxyy_yyy_0, \
                             ta1_x_xxyy_yyz_0, \
                             ta1_x_xxyy_yzz_0, \
                             ta1_x_xxyy_zzz_0, \
                             ta1_x_xyy_yyy_0,  \
                             ta1_x_xyy_yyy_1,  \
                             ta1_x_xyy_yyz_0,  \
                             ta1_x_xyy_yyz_1,  \
                             ta1_x_xyy_yzz_0,  \
                             ta1_x_xyy_yzz_1,  \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yy_yyz_0,   \
                             ta1_x_yy_yyz_1,   \
                             ta1_x_yy_yzz_0,   \
                             ta1_x_yy_yzz_1,   \
                             ta_xyy_yyy_1,     \
                             ta_xyy_yyz_1,     \
                             ta_xyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_xxx_0[i] = ta1_x_xx_xxx_0[i] * fe_0 - ta1_x_xx_xxx_1[i] * fe_0 + ta1_x_xxy_xxx_0[i] * pa_y[i] - ta1_x_xxy_xxx_1[i] * pc_y[i];

        ta1_x_xxyy_xxy_0[i] = ta1_x_xx_xxy_0[i] * fe_0 - ta1_x_xx_xxy_1[i] * fe_0 + ta1_x_xxy_xx_0[i] * fe_0 - ta1_x_xxy_xx_1[i] * fe_0 +
                              ta1_x_xxy_xxy_0[i] * pa_y[i] - ta1_x_xxy_xxy_1[i] * pc_y[i];

        ta1_x_xxyy_xxz_0[i] = ta1_x_xx_xxz_0[i] * fe_0 - ta1_x_xx_xxz_1[i] * fe_0 + ta1_x_xxy_xxz_0[i] * pa_y[i] - ta1_x_xxy_xxz_1[i] * pc_y[i];

        ta1_x_xxyy_xyy_0[i] = ta1_x_xx_xyy_0[i] * fe_0 - ta1_x_xx_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxy_xy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xy_1[i] * fe_0 +
                              ta1_x_xxy_xyy_0[i] * pa_y[i] - ta1_x_xxy_xyy_1[i] * pc_y[i];

        ta1_x_xxyy_xyz_0[i] = ta1_x_xx_xyz_0[i] * fe_0 - ta1_x_xx_xyz_1[i] * fe_0 + ta1_x_xxy_xz_0[i] * fe_0 - ta1_x_xxy_xz_1[i] * fe_0 +
                              ta1_x_xxy_xyz_0[i] * pa_y[i] - ta1_x_xxy_xyz_1[i] * pc_y[i];

        ta1_x_xxyy_xzz_0[i] = ta1_x_xx_xzz_0[i] * fe_0 - ta1_x_xx_xzz_1[i] * fe_0 + ta1_x_xxy_xzz_0[i] * pa_y[i] - ta1_x_xxy_xzz_1[i] * pc_y[i];

        ta1_x_xxyy_yyy_0[i] =
            ta1_x_yy_yyy_0[i] * fe_0 - ta1_x_yy_yyy_1[i] * fe_0 + ta_xyy_yyy_1[i] + ta1_x_xyy_yyy_0[i] * pa_x[i] - ta1_x_xyy_yyy_1[i] * pc_x[i];

        ta1_x_xxyy_yyz_0[i] =
            ta1_x_yy_yyz_0[i] * fe_0 - ta1_x_yy_yyz_1[i] * fe_0 + ta_xyy_yyz_1[i] + ta1_x_xyy_yyz_0[i] * pa_x[i] - ta1_x_xyy_yyz_1[i] * pc_x[i];

        ta1_x_xxyy_yzz_0[i] =
            ta1_x_yy_yzz_0[i] * fe_0 - ta1_x_yy_yzz_1[i] * fe_0 + ta_xyy_yzz_1[i] + ta1_x_xyy_yzz_0[i] * pa_x[i] - ta1_x_xyy_yzz_1[i] * pc_x[i];

        ta1_x_xxyy_zzz_0[i] = ta1_x_xx_zzz_0[i] * fe_0 - ta1_x_xx_zzz_1[i] * fe_0 + ta1_x_xxy_zzz_0[i] * pa_y[i] - ta1_x_xxy_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto ta1_x_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 40);

    auto ta1_x_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 41);

    auto ta1_x_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 42);

    auto ta1_x_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 43);

    auto ta1_x_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 44);

    auto ta1_x_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 45);

    auto ta1_x_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 46);

    auto ta1_x_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 47);

    auto ta1_x_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 48);

    auto ta1_x_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 49);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxy_xxy_0,  \
                             ta1_x_xxy_xxy_1,  \
                             ta1_x_xxy_xyy_0,  \
                             ta1_x_xxy_xyy_1,  \
                             ta1_x_xxy_yyy_0,  \
                             ta1_x_xxy_yyy_1,  \
                             ta1_x_xxyz_xxx_0, \
                             ta1_x_xxyz_xxy_0, \
                             ta1_x_xxyz_xxz_0, \
                             ta1_x_xxyz_xyy_0, \
                             ta1_x_xxyz_xyz_0, \
                             ta1_x_xxyz_xzz_0, \
                             ta1_x_xxyz_yyy_0, \
                             ta1_x_xxyz_yyz_0, \
                             ta1_x_xxyz_yzz_0, \
                             ta1_x_xxyz_zzz_0, \
                             ta1_x_xxz_xxx_0,  \
                             ta1_x_xxz_xxx_1,  \
                             ta1_x_xxz_xxz_0,  \
                             ta1_x_xxz_xxz_1,  \
                             ta1_x_xxz_xyz_0,  \
                             ta1_x_xxz_xyz_1,  \
                             ta1_x_xxz_xz_0,   \
                             ta1_x_xxz_xz_1,   \
                             ta1_x_xxz_xzz_0,  \
                             ta1_x_xxz_xzz_1,  \
                             ta1_x_xxz_yyz_0,  \
                             ta1_x_xxz_yyz_1,  \
                             ta1_x_xxz_yz_0,   \
                             ta1_x_xxz_yz_1,   \
                             ta1_x_xxz_yzz_0,  \
                             ta1_x_xxz_yzz_1,  \
                             ta1_x_xxz_zz_0,   \
                             ta1_x_xxz_zz_1,   \
                             ta1_x_xxz_zzz_0,  \
                             ta1_x_xxz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyz_xxx_0[i] = ta1_x_xxz_xxx_0[i] * pa_y[i] - ta1_x_xxz_xxx_1[i] * pc_y[i];

        ta1_x_xxyz_xxy_0[i] = ta1_x_xxy_xxy_0[i] * pa_z[i] - ta1_x_xxy_xxy_1[i] * pc_z[i];

        ta1_x_xxyz_xxz_0[i] = ta1_x_xxz_xxz_0[i] * pa_y[i] - ta1_x_xxz_xxz_1[i] * pc_y[i];

        ta1_x_xxyz_xyy_0[i] = ta1_x_xxy_xyy_0[i] * pa_z[i] - ta1_x_xxy_xyy_1[i] * pc_z[i];

        ta1_x_xxyz_xyz_0[i] = ta1_x_xxz_xz_0[i] * fe_0 - ta1_x_xxz_xz_1[i] * fe_0 + ta1_x_xxz_xyz_0[i] * pa_y[i] - ta1_x_xxz_xyz_1[i] * pc_y[i];

        ta1_x_xxyz_xzz_0[i] = ta1_x_xxz_xzz_0[i] * pa_y[i] - ta1_x_xxz_xzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyy_0[i] = ta1_x_xxy_yyy_0[i] * pa_z[i] - ta1_x_xxy_yyy_1[i] * pc_z[i];

        ta1_x_xxyz_yyz_0[i] =
            2.0 * ta1_x_xxz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxz_yz_1[i] * fe_0 + ta1_x_xxz_yyz_0[i] * pa_y[i] - ta1_x_xxz_yyz_1[i] * pc_y[i];

        ta1_x_xxyz_yzz_0[i] = ta1_x_xxz_zz_0[i] * fe_0 - ta1_x_xxz_zz_1[i] * fe_0 + ta1_x_xxz_yzz_0[i] * pa_y[i] - ta1_x_xxz_yzz_1[i] * pc_y[i];

        ta1_x_xxyz_zzz_0[i] = ta1_x_xxz_zzz_0[i] * pa_y[i] - ta1_x_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto ta1_x_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 50);

    auto ta1_x_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 51);

    auto ta1_x_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 52);

    auto ta1_x_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 53);

    auto ta1_x_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 54);

    auto ta1_x_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 55);

    auto ta1_x_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 56);

    auto ta1_x_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 57);

    auto ta1_x_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 58);

    auto ta1_x_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_yyy_0,   \
                             ta1_x_xx_yyy_1,   \
                             ta1_x_xxz_xx_0,   \
                             ta1_x_xxz_xx_1,   \
                             ta1_x_xxz_xxx_0,  \
                             ta1_x_xxz_xxx_1,  \
                             ta1_x_xxz_xxy_0,  \
                             ta1_x_xxz_xxy_1,  \
                             ta1_x_xxz_xxz_0,  \
                             ta1_x_xxz_xxz_1,  \
                             ta1_x_xxz_xy_0,   \
                             ta1_x_xxz_xy_1,   \
                             ta1_x_xxz_xyy_0,  \
                             ta1_x_xxz_xyy_1,  \
                             ta1_x_xxz_xyz_0,  \
                             ta1_x_xxz_xyz_1,  \
                             ta1_x_xxz_xz_0,   \
                             ta1_x_xxz_xz_1,   \
                             ta1_x_xxz_xzz_0,  \
                             ta1_x_xxz_xzz_1,  \
                             ta1_x_xxz_yyy_0,  \
                             ta1_x_xxz_yyy_1,  \
                             ta1_x_xxzz_xxx_0, \
                             ta1_x_xxzz_xxy_0, \
                             ta1_x_xxzz_xxz_0, \
                             ta1_x_xxzz_xyy_0, \
                             ta1_x_xxzz_xyz_0, \
                             ta1_x_xxzz_xzz_0, \
                             ta1_x_xxzz_yyy_0, \
                             ta1_x_xxzz_yyz_0, \
                             ta1_x_xxzz_yzz_0, \
                             ta1_x_xxzz_zzz_0, \
                             ta1_x_xzz_yyz_0,  \
                             ta1_x_xzz_yyz_1,  \
                             ta1_x_xzz_yzz_0,  \
                             ta1_x_xzz_yzz_1,  \
                             ta1_x_xzz_zzz_0,  \
                             ta1_x_xzz_zzz_1,  \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             ta_xzz_yyz_1,     \
                             ta_xzz_yzz_1,     \
                             ta_xzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_xxx_0[i] = ta1_x_xx_xxx_0[i] * fe_0 - ta1_x_xx_xxx_1[i] * fe_0 + ta1_x_xxz_xxx_0[i] * pa_z[i] - ta1_x_xxz_xxx_1[i] * pc_z[i];

        ta1_x_xxzz_xxy_0[i] = ta1_x_xx_xxy_0[i] * fe_0 - ta1_x_xx_xxy_1[i] * fe_0 + ta1_x_xxz_xxy_0[i] * pa_z[i] - ta1_x_xxz_xxy_1[i] * pc_z[i];

        ta1_x_xxzz_xxz_0[i] = ta1_x_xx_xxz_0[i] * fe_0 - ta1_x_xx_xxz_1[i] * fe_0 + ta1_x_xxz_xx_0[i] * fe_0 - ta1_x_xxz_xx_1[i] * fe_0 +
                              ta1_x_xxz_xxz_0[i] * pa_z[i] - ta1_x_xxz_xxz_1[i] * pc_z[i];

        ta1_x_xxzz_xyy_0[i] = ta1_x_xx_xyy_0[i] * fe_0 - ta1_x_xx_xyy_1[i] * fe_0 + ta1_x_xxz_xyy_0[i] * pa_z[i] - ta1_x_xxz_xyy_1[i] * pc_z[i];

        ta1_x_xxzz_xyz_0[i] = ta1_x_xx_xyz_0[i] * fe_0 - ta1_x_xx_xyz_1[i] * fe_0 + ta1_x_xxz_xy_0[i] * fe_0 - ta1_x_xxz_xy_1[i] * fe_0 +
                              ta1_x_xxz_xyz_0[i] * pa_z[i] - ta1_x_xxz_xyz_1[i] * pc_z[i];

        ta1_x_xxzz_xzz_0[i] = ta1_x_xx_xzz_0[i] * fe_0 - ta1_x_xx_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xz_1[i] * fe_0 +
                              ta1_x_xxz_xzz_0[i] * pa_z[i] - ta1_x_xxz_xzz_1[i] * pc_z[i];

        ta1_x_xxzz_yyy_0[i] = ta1_x_xx_yyy_0[i] * fe_0 - ta1_x_xx_yyy_1[i] * fe_0 + ta1_x_xxz_yyy_0[i] * pa_z[i] - ta1_x_xxz_yyy_1[i] * pc_z[i];

        ta1_x_xxzz_yyz_0[i] =
            ta1_x_zz_yyz_0[i] * fe_0 - ta1_x_zz_yyz_1[i] * fe_0 + ta_xzz_yyz_1[i] + ta1_x_xzz_yyz_0[i] * pa_x[i] - ta1_x_xzz_yyz_1[i] * pc_x[i];

        ta1_x_xxzz_yzz_0[i] =
            ta1_x_zz_yzz_0[i] * fe_0 - ta1_x_zz_yzz_1[i] * fe_0 + ta_xzz_yzz_1[i] + ta1_x_xzz_yzz_0[i] * pa_x[i] - ta1_x_xzz_yzz_1[i] * pc_x[i];

        ta1_x_xxzz_zzz_0[i] =
            ta1_x_zz_zzz_0[i] * fe_0 - ta1_x_zz_zzz_1[i] * fe_0 + ta_xzz_zzz_1[i] + ta1_x_xzz_zzz_0[i] * pa_x[i] - ta1_x_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto ta1_x_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 60);

    auto ta1_x_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 61);

    auto ta1_x_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 62);

    auto ta1_x_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 63);

    auto ta1_x_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 64);

    auto ta1_x_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 65);

    auto ta1_x_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 66);

    auto ta1_x_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 67);

    auto ta1_x_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 68);

    auto ta1_x_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 69);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xy_xxx_0,   \
                             ta1_x_xy_xxx_1,   \
                             ta1_x_xy_xxz_0,   \
                             ta1_x_xy_xxz_1,   \
                             ta1_x_xy_xzz_0,   \
                             ta1_x_xy_xzz_1,   \
                             ta1_x_xyy_xxx_0,  \
                             ta1_x_xyy_xxx_1,  \
                             ta1_x_xyy_xxz_0,  \
                             ta1_x_xyy_xxz_1,  \
                             ta1_x_xyy_xzz_0,  \
                             ta1_x_xyy_xzz_1,  \
                             ta1_x_xyyy_xxx_0, \
                             ta1_x_xyyy_xxy_0, \
                             ta1_x_xyyy_xxz_0, \
                             ta1_x_xyyy_xyy_0, \
                             ta1_x_xyyy_xyz_0, \
                             ta1_x_xyyy_xzz_0, \
                             ta1_x_xyyy_yyy_0, \
                             ta1_x_xyyy_yyz_0, \
                             ta1_x_xyyy_yzz_0, \
                             ta1_x_xyyy_zzz_0, \
                             ta1_x_yyy_xxy_0,  \
                             ta1_x_yyy_xxy_1,  \
                             ta1_x_yyy_xy_0,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_xyy_0,  \
                             ta1_x_yyy_xyy_1,  \
                             ta1_x_yyy_xyz_0,  \
                             ta1_x_yyy_xyz_1,  \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yyy_0,  \
                             ta1_x_yyy_yyy_1,  \
                             ta1_x_yyy_yyz_0,  \
                             ta1_x_yyy_yyz_1,  \
                             ta1_x_yyy_yz_0,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_yzz_0,  \
                             ta1_x_yyy_yzz_1,  \
                             ta1_x_yyy_zzz_0,  \
                             ta1_x_yyy_zzz_1,  \
                             ta_yyy_xxy_1,     \
                             ta_yyy_xyy_1,     \
                             ta_yyy_xyz_1,     \
                             ta_yyy_yyy_1,     \
                             ta_yyy_yyz_1,     \
                             ta_yyy_yzz_1,     \
                             ta_yyy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_xxx_0[i] =
            2.0 * ta1_x_xy_xxx_0[i] * fe_0 - 2.0 * ta1_x_xy_xxx_1[i] * fe_0 + ta1_x_xyy_xxx_0[i] * pa_y[i] - ta1_x_xyy_xxx_1[i] * pc_y[i];

        ta1_x_xyyy_xxy_0[i] = 2.0 * ta1_x_yyy_xy_0[i] * fe_0 - 2.0 * ta1_x_yyy_xy_1[i] * fe_0 + ta_yyy_xxy_1[i] + ta1_x_yyy_xxy_0[i] * pa_x[i] -
                              ta1_x_yyy_xxy_1[i] * pc_x[i];

        ta1_x_xyyy_xxz_0[i] =
            2.0 * ta1_x_xy_xxz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxz_1[i] * fe_0 + ta1_x_xyy_xxz_0[i] * pa_y[i] - ta1_x_xyy_xxz_1[i] * pc_y[i];

        ta1_x_xyyy_xyy_0[i] =
            ta1_x_yyy_yy_0[i] * fe_0 - ta1_x_yyy_yy_1[i] * fe_0 + ta_yyy_xyy_1[i] + ta1_x_yyy_xyy_0[i] * pa_x[i] - ta1_x_yyy_xyy_1[i] * pc_x[i];

        ta1_x_xyyy_xyz_0[i] =
            ta1_x_yyy_yz_0[i] * fe_0 - ta1_x_yyy_yz_1[i] * fe_0 + ta_yyy_xyz_1[i] + ta1_x_yyy_xyz_0[i] * pa_x[i] - ta1_x_yyy_xyz_1[i] * pc_x[i];

        ta1_x_xyyy_xzz_0[i] =
            2.0 * ta1_x_xy_xzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xzz_1[i] * fe_0 + ta1_x_xyy_xzz_0[i] * pa_y[i] - ta1_x_xyy_xzz_1[i] * pc_y[i];

        ta1_x_xyyy_yyy_0[i] = ta_yyy_yyy_1[i] + ta1_x_yyy_yyy_0[i] * pa_x[i] - ta1_x_yyy_yyy_1[i] * pc_x[i];

        ta1_x_xyyy_yyz_0[i] = ta_yyy_yyz_1[i] + ta1_x_yyy_yyz_0[i] * pa_x[i] - ta1_x_yyy_yyz_1[i] * pc_x[i];

        ta1_x_xyyy_yzz_0[i] = ta_yyy_yzz_1[i] + ta1_x_yyy_yzz_0[i] * pa_x[i] - ta1_x_yyy_yzz_1[i] * pc_x[i];

        ta1_x_xyyy_zzz_0[i] = ta_yyy_zzz_1[i] + ta1_x_yyy_zzz_0[i] * pa_x[i] - ta1_x_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto ta1_x_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 70);

    auto ta1_x_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 71);

    auto ta1_x_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 72);

    auto ta1_x_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 73);

    auto ta1_x_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 74);

    auto ta1_x_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 75);

    auto ta1_x_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 76);

    auto ta1_x_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 77);

    auto ta1_x_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 78);

    auto ta1_x_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 79);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xyy_xxx_0,  \
                             ta1_x_xyy_xxx_1,  \
                             ta1_x_xyy_xxy_0,  \
                             ta1_x_xyy_xxy_1,  \
                             ta1_x_xyy_xy_0,   \
                             ta1_x_xyy_xy_1,   \
                             ta1_x_xyy_xyy_0,  \
                             ta1_x_xyy_xyy_1,  \
                             ta1_x_xyy_xyz_0,  \
                             ta1_x_xyy_xyz_1,  \
                             ta1_x_xyy_yyy_0,  \
                             ta1_x_xyy_yyy_1,  \
                             ta1_x_xyyz_xxx_0, \
                             ta1_x_xyyz_xxy_0, \
                             ta1_x_xyyz_xxz_0, \
                             ta1_x_xyyz_xyy_0, \
                             ta1_x_xyyz_xyz_0, \
                             ta1_x_xyyz_xzz_0, \
                             ta1_x_xyyz_yyy_0, \
                             ta1_x_xyyz_yyz_0, \
                             ta1_x_xyyz_yzz_0, \
                             ta1_x_xyyz_zzz_0, \
                             ta1_x_xyz_xxz_0,  \
                             ta1_x_xyz_xxz_1,  \
                             ta1_x_xyz_xzz_0,  \
                             ta1_x_xyz_xzz_1,  \
                             ta1_x_xz_xxz_0,   \
                             ta1_x_xz_xxz_1,   \
                             ta1_x_xz_xzz_0,   \
                             ta1_x_xz_xzz_1,   \
                             ta1_x_yyz_yyz_0,  \
                             ta1_x_yyz_yyz_1,  \
                             ta1_x_yyz_yzz_0,  \
                             ta1_x_yyz_yzz_1,  \
                             ta1_x_yyz_zzz_0,  \
                             ta1_x_yyz_zzz_1,  \
                             ta_yyz_yyz_1,     \
                             ta_yyz_yzz_1,     \
                             ta_yyz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyz_xxx_0[i] = ta1_x_xyy_xxx_0[i] * pa_z[i] - ta1_x_xyy_xxx_1[i] * pc_z[i];

        ta1_x_xyyz_xxy_0[i] = ta1_x_xyy_xxy_0[i] * pa_z[i] - ta1_x_xyy_xxy_1[i] * pc_z[i];

        ta1_x_xyyz_xxz_0[i] = ta1_x_xz_xxz_0[i] * fe_0 - ta1_x_xz_xxz_1[i] * fe_0 + ta1_x_xyz_xxz_0[i] * pa_y[i] - ta1_x_xyz_xxz_1[i] * pc_y[i];

        ta1_x_xyyz_xyy_0[i] = ta1_x_xyy_xyy_0[i] * pa_z[i] - ta1_x_xyy_xyy_1[i] * pc_z[i];

        ta1_x_xyyz_xyz_0[i] = ta1_x_xyy_xy_0[i] * fe_0 - ta1_x_xyy_xy_1[i] * fe_0 + ta1_x_xyy_xyz_0[i] * pa_z[i] - ta1_x_xyy_xyz_1[i] * pc_z[i];

        ta1_x_xyyz_xzz_0[i] = ta1_x_xz_xzz_0[i] * fe_0 - ta1_x_xz_xzz_1[i] * fe_0 + ta1_x_xyz_xzz_0[i] * pa_y[i] - ta1_x_xyz_xzz_1[i] * pc_y[i];

        ta1_x_xyyz_yyy_0[i] = ta1_x_xyy_yyy_0[i] * pa_z[i] - ta1_x_xyy_yyy_1[i] * pc_z[i];

        ta1_x_xyyz_yyz_0[i] = ta_yyz_yyz_1[i] + ta1_x_yyz_yyz_0[i] * pa_x[i] - ta1_x_yyz_yyz_1[i] * pc_x[i];

        ta1_x_xyyz_yzz_0[i] = ta_yyz_yzz_1[i] + ta1_x_yyz_yzz_0[i] * pa_x[i] - ta1_x_yyz_yzz_1[i] * pc_x[i];

        ta1_x_xyyz_zzz_0[i] = ta_yyz_zzz_1[i] + ta1_x_yyz_zzz_0[i] * pa_x[i] - ta1_x_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto ta1_x_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 80);

    auto ta1_x_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 81);

    auto ta1_x_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 82);

    auto ta1_x_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 83);

    auto ta1_x_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 84);

    auto ta1_x_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 85);

    auto ta1_x_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 86);

    auto ta1_x_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 87);

    auto ta1_x_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 88);

    auto ta1_x_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 89);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xyzz_xxx_0, \
                             ta1_x_xyzz_xxy_0, \
                             ta1_x_xyzz_xxz_0, \
                             ta1_x_xyzz_xyy_0, \
                             ta1_x_xyzz_xyz_0, \
                             ta1_x_xyzz_xzz_0, \
                             ta1_x_xyzz_yyy_0, \
                             ta1_x_xyzz_yyz_0, \
                             ta1_x_xyzz_yzz_0, \
                             ta1_x_xyzz_zzz_0, \
                             ta1_x_xzz_xx_0,   \
                             ta1_x_xzz_xx_1,   \
                             ta1_x_xzz_xxx_0,  \
                             ta1_x_xzz_xxx_1,  \
                             ta1_x_xzz_xxy_0,  \
                             ta1_x_xzz_xxy_1,  \
                             ta1_x_xzz_xxz_0,  \
                             ta1_x_xzz_xxz_1,  \
                             ta1_x_xzz_xy_0,   \
                             ta1_x_xzz_xy_1,   \
                             ta1_x_xzz_xyy_0,  \
                             ta1_x_xzz_xyy_1,  \
                             ta1_x_xzz_xyz_0,  \
                             ta1_x_xzz_xyz_1,  \
                             ta1_x_xzz_xz_0,   \
                             ta1_x_xzz_xz_1,   \
                             ta1_x_xzz_xzz_0,  \
                             ta1_x_xzz_xzz_1,  \
                             ta1_x_xzz_zzz_0,  \
                             ta1_x_xzz_zzz_1,  \
                             ta1_x_yzz_yyy_0,  \
                             ta1_x_yzz_yyy_1,  \
                             ta1_x_yzz_yyz_0,  \
                             ta1_x_yzz_yyz_1,  \
                             ta1_x_yzz_yzz_0,  \
                             ta1_x_yzz_yzz_1,  \
                             ta_yzz_yyy_1,     \
                             ta_yzz_yyz_1,     \
                             ta_yzz_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzz_xxx_0[i] = ta1_x_xzz_xxx_0[i] * pa_y[i] - ta1_x_xzz_xxx_1[i] * pc_y[i];

        ta1_x_xyzz_xxy_0[i] = ta1_x_xzz_xx_0[i] * fe_0 - ta1_x_xzz_xx_1[i] * fe_0 + ta1_x_xzz_xxy_0[i] * pa_y[i] - ta1_x_xzz_xxy_1[i] * pc_y[i];

        ta1_x_xyzz_xxz_0[i] = ta1_x_xzz_xxz_0[i] * pa_y[i] - ta1_x_xzz_xxz_1[i] * pc_y[i];

        ta1_x_xyzz_xyy_0[i] =
            2.0 * ta1_x_xzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xzz_xy_1[i] * fe_0 + ta1_x_xzz_xyy_0[i] * pa_y[i] - ta1_x_xzz_xyy_1[i] * pc_y[i];

        ta1_x_xyzz_xyz_0[i] = ta1_x_xzz_xz_0[i] * fe_0 - ta1_x_xzz_xz_1[i] * fe_0 + ta1_x_xzz_xyz_0[i] * pa_y[i] - ta1_x_xzz_xyz_1[i] * pc_y[i];

        ta1_x_xyzz_xzz_0[i] = ta1_x_xzz_xzz_0[i] * pa_y[i] - ta1_x_xzz_xzz_1[i] * pc_y[i];

        ta1_x_xyzz_yyy_0[i] = ta_yzz_yyy_1[i] + ta1_x_yzz_yyy_0[i] * pa_x[i] - ta1_x_yzz_yyy_1[i] * pc_x[i];

        ta1_x_xyzz_yyz_0[i] = ta_yzz_yyz_1[i] + ta1_x_yzz_yyz_0[i] * pa_x[i] - ta1_x_yzz_yyz_1[i] * pc_x[i];

        ta1_x_xyzz_yzz_0[i] = ta_yzz_yzz_1[i] + ta1_x_yzz_yzz_0[i] * pa_x[i] - ta1_x_yzz_yzz_1[i] * pc_x[i];

        ta1_x_xyzz_zzz_0[i] = ta1_x_xzz_zzz_0[i] * pa_y[i] - ta1_x_xzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto ta1_x_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 90);

    auto ta1_x_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 91);

    auto ta1_x_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 92);

    auto ta1_x_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 93);

    auto ta1_x_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 94);

    auto ta1_x_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 95);

    auto ta1_x_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 96);

    auto ta1_x_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 97);

    auto ta1_x_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 98);

    auto ta1_x_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 99);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xz_xxx_0,   \
                             ta1_x_xz_xxx_1,   \
                             ta1_x_xz_xxy_0,   \
                             ta1_x_xz_xxy_1,   \
                             ta1_x_xz_xyy_0,   \
                             ta1_x_xz_xyy_1,   \
                             ta1_x_xzz_xxx_0,  \
                             ta1_x_xzz_xxx_1,  \
                             ta1_x_xzz_xxy_0,  \
                             ta1_x_xzz_xxy_1,  \
                             ta1_x_xzz_xyy_0,  \
                             ta1_x_xzz_xyy_1,  \
                             ta1_x_xzzz_xxx_0, \
                             ta1_x_xzzz_xxy_0, \
                             ta1_x_xzzz_xxz_0, \
                             ta1_x_xzzz_xyy_0, \
                             ta1_x_xzzz_xyz_0, \
                             ta1_x_xzzz_xzz_0, \
                             ta1_x_xzzz_yyy_0, \
                             ta1_x_xzzz_yyz_0, \
                             ta1_x_xzzz_yzz_0, \
                             ta1_x_xzzz_zzz_0, \
                             ta1_x_zzz_xxz_0,  \
                             ta1_x_zzz_xxz_1,  \
                             ta1_x_zzz_xyz_0,  \
                             ta1_x_zzz_xyz_1,  \
                             ta1_x_zzz_xz_0,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_xzz_0,  \
                             ta1_x_zzz_xzz_1,  \
                             ta1_x_zzz_yyy_0,  \
                             ta1_x_zzz_yyy_1,  \
                             ta1_x_zzz_yyz_0,  \
                             ta1_x_zzz_yyz_1,  \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_yzz_0,  \
                             ta1_x_zzz_yzz_1,  \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             ta1_x_zzz_zzz_0,  \
                             ta1_x_zzz_zzz_1,  \
                             ta_zzz_xxz_1,     \
                             ta_zzz_xyz_1,     \
                             ta_zzz_xzz_1,     \
                             ta_zzz_yyy_1,     \
                             ta_zzz_yyz_1,     \
                             ta_zzz_yzz_1,     \
                             ta_zzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_xxx_0[i] =
            2.0 * ta1_x_xz_xxx_0[i] * fe_0 - 2.0 * ta1_x_xz_xxx_1[i] * fe_0 + ta1_x_xzz_xxx_0[i] * pa_z[i] - ta1_x_xzz_xxx_1[i] * pc_z[i];

        ta1_x_xzzz_xxy_0[i] =
            2.0 * ta1_x_xz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxy_1[i] * fe_0 + ta1_x_xzz_xxy_0[i] * pa_z[i] - ta1_x_xzz_xxy_1[i] * pc_z[i];

        ta1_x_xzzz_xxz_0[i] = 2.0 * ta1_x_zzz_xz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xz_1[i] * fe_0 + ta_zzz_xxz_1[i] + ta1_x_zzz_xxz_0[i] * pa_x[i] -
                              ta1_x_zzz_xxz_1[i] * pc_x[i];

        ta1_x_xzzz_xyy_0[i] =
            2.0 * ta1_x_xz_xyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xyy_1[i] * fe_0 + ta1_x_xzz_xyy_0[i] * pa_z[i] - ta1_x_xzz_xyy_1[i] * pc_z[i];

        ta1_x_xzzz_xyz_0[i] =
            ta1_x_zzz_yz_0[i] * fe_0 - ta1_x_zzz_yz_1[i] * fe_0 + ta_zzz_xyz_1[i] + ta1_x_zzz_xyz_0[i] * pa_x[i] - ta1_x_zzz_xyz_1[i] * pc_x[i];

        ta1_x_xzzz_xzz_0[i] =
            ta1_x_zzz_zz_0[i] * fe_0 - ta1_x_zzz_zz_1[i] * fe_0 + ta_zzz_xzz_1[i] + ta1_x_zzz_xzz_0[i] * pa_x[i] - ta1_x_zzz_xzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyy_0[i] = ta_zzz_yyy_1[i] + ta1_x_zzz_yyy_0[i] * pa_x[i] - ta1_x_zzz_yyy_1[i] * pc_x[i];

        ta1_x_xzzz_yyz_0[i] = ta_zzz_yyz_1[i] + ta1_x_zzz_yyz_0[i] * pa_x[i] - ta1_x_zzz_yyz_1[i] * pc_x[i];

        ta1_x_xzzz_yzz_0[i] = ta_zzz_yzz_1[i] + ta1_x_zzz_yzz_0[i] * pa_x[i] - ta1_x_zzz_yzz_1[i] * pc_x[i];

        ta1_x_xzzz_zzz_0[i] = ta_zzz_zzz_1[i] + ta1_x_zzz_zzz_0[i] * pa_x[i] - ta1_x_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto ta1_x_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 100);

    auto ta1_x_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 101);

    auto ta1_x_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 102);

    auto ta1_x_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 103);

    auto ta1_x_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 104);

    auto ta1_x_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 105);

    auto ta1_x_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 106);

    auto ta1_x_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 107);

    auto ta1_x_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 108);

    auto ta1_x_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 109);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yy_xxx_0,   \
                             ta1_x_yy_xxx_1,   \
                             ta1_x_yy_xxy_0,   \
                             ta1_x_yy_xxy_1,   \
                             ta1_x_yy_xxz_0,   \
                             ta1_x_yy_xxz_1,   \
                             ta1_x_yy_xyy_0,   \
                             ta1_x_yy_xyy_1,   \
                             ta1_x_yy_xyz_0,   \
                             ta1_x_yy_xyz_1,   \
                             ta1_x_yy_xzz_0,   \
                             ta1_x_yy_xzz_1,   \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yy_yyz_0,   \
                             ta1_x_yy_yyz_1,   \
                             ta1_x_yy_yzz_0,   \
                             ta1_x_yy_yzz_1,   \
                             ta1_x_yy_zzz_0,   \
                             ta1_x_yy_zzz_1,   \
                             ta1_x_yyy_xx_0,   \
                             ta1_x_yyy_xx_1,   \
                             ta1_x_yyy_xxx_0,  \
                             ta1_x_yyy_xxx_1,  \
                             ta1_x_yyy_xxy_0,  \
                             ta1_x_yyy_xxy_1,  \
                             ta1_x_yyy_xxz_0,  \
                             ta1_x_yyy_xxz_1,  \
                             ta1_x_yyy_xy_0,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_xyy_0,  \
                             ta1_x_yyy_xyy_1,  \
                             ta1_x_yyy_xyz_0,  \
                             ta1_x_yyy_xyz_1,  \
                             ta1_x_yyy_xz_0,   \
                             ta1_x_yyy_xz_1,   \
                             ta1_x_yyy_xzz_0,  \
                             ta1_x_yyy_xzz_1,  \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yyy_0,  \
                             ta1_x_yyy_yyy_1,  \
                             ta1_x_yyy_yyz_0,  \
                             ta1_x_yyy_yyz_1,  \
                             ta1_x_yyy_yz_0,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_yzz_0,  \
                             ta1_x_yyy_yzz_1,  \
                             ta1_x_yyy_zz_0,   \
                             ta1_x_yyy_zz_1,   \
                             ta1_x_yyy_zzz_0,  \
                             ta1_x_yyy_zzz_1,  \
                             ta1_x_yyyy_xxx_0, \
                             ta1_x_yyyy_xxy_0, \
                             ta1_x_yyyy_xxz_0, \
                             ta1_x_yyyy_xyy_0, \
                             ta1_x_yyyy_xyz_0, \
                             ta1_x_yyyy_xzz_0, \
                             ta1_x_yyyy_yyy_0, \
                             ta1_x_yyyy_yyz_0, \
                             ta1_x_yyyy_yzz_0, \
                             ta1_x_yyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_xxx_0[i] =
            3.0 * ta1_x_yy_xxx_0[i] * fe_0 - 3.0 * ta1_x_yy_xxx_1[i] * fe_0 + ta1_x_yyy_xxx_0[i] * pa_y[i] - ta1_x_yyy_xxx_1[i] * pc_y[i];

        ta1_x_yyyy_xxy_0[i] = 3.0 * ta1_x_yy_xxy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxy_1[i] * fe_0 + ta1_x_yyy_xx_0[i] * fe_0 - ta1_x_yyy_xx_1[i] * fe_0 +
                              ta1_x_yyy_xxy_0[i] * pa_y[i] - ta1_x_yyy_xxy_1[i] * pc_y[i];

        ta1_x_yyyy_xxz_0[i] =
            3.0 * ta1_x_yy_xxz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxz_1[i] * fe_0 + ta1_x_yyy_xxz_0[i] * pa_y[i] - ta1_x_yyy_xxz_1[i] * pc_y[i];

        ta1_x_yyyy_xyy_0[i] = 3.0 * ta1_x_yy_xyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xy_0[i] * fe_0 -
                              2.0 * ta1_x_yyy_xy_1[i] * fe_0 + ta1_x_yyy_xyy_0[i] * pa_y[i] - ta1_x_yyy_xyy_1[i] * pc_y[i];

        ta1_x_yyyy_xyz_0[i] = 3.0 * ta1_x_yy_xyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyz_1[i] * fe_0 + ta1_x_yyy_xz_0[i] * fe_0 - ta1_x_yyy_xz_1[i] * fe_0 +
                              ta1_x_yyy_xyz_0[i] * pa_y[i] - ta1_x_yyy_xyz_1[i] * pc_y[i];

        ta1_x_yyyy_xzz_0[i] =
            3.0 * ta1_x_yy_xzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xzz_1[i] * fe_0 + ta1_x_yyy_xzz_0[i] * pa_y[i] - ta1_x_yyy_xzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyy_0[i] = 3.0 * ta1_x_yy_yyy_0[i] * fe_0 - 3.0 * ta1_x_yy_yyy_1[i] * fe_0 + 3.0 * ta1_x_yyy_yy_0[i] * fe_0 -
                              3.0 * ta1_x_yyy_yy_1[i] * fe_0 + ta1_x_yyy_yyy_0[i] * pa_y[i] - ta1_x_yyy_yyy_1[i] * pc_y[i];

        ta1_x_yyyy_yyz_0[i] = 3.0 * ta1_x_yy_yyz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_yz_0[i] * fe_0 -
                              2.0 * ta1_x_yyy_yz_1[i] * fe_0 + ta1_x_yyy_yyz_0[i] * pa_y[i] - ta1_x_yyy_yyz_1[i] * pc_y[i];

        ta1_x_yyyy_yzz_0[i] = 3.0 * ta1_x_yy_yzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yzz_1[i] * fe_0 + ta1_x_yyy_zz_0[i] * fe_0 - ta1_x_yyy_zz_1[i] * fe_0 +
                              ta1_x_yyy_yzz_0[i] * pa_y[i] - ta1_x_yyy_yzz_1[i] * pc_y[i];

        ta1_x_yyyy_zzz_0[i] =
            3.0 * ta1_x_yy_zzz_0[i] * fe_0 - 3.0 * ta1_x_yy_zzz_1[i] * fe_0 + ta1_x_yyy_zzz_0[i] * pa_y[i] - ta1_x_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto ta1_x_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 110);

    auto ta1_x_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 111);

    auto ta1_x_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 112);

    auto ta1_x_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 113);

    auto ta1_x_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 114);

    auto ta1_x_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 115);

    auto ta1_x_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 116);

    auto ta1_x_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 117);

    auto ta1_x_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 118);

    auto ta1_x_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 119);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyy_xxx_0,  \
                             ta1_x_yyy_xxx_1,  \
                             ta1_x_yyy_xxy_0,  \
                             ta1_x_yyy_xxy_1,  \
                             ta1_x_yyy_xy_0,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_xyy_0,  \
                             ta1_x_yyy_xyy_1,  \
                             ta1_x_yyy_xyz_0,  \
                             ta1_x_yyy_xyz_1,  \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yyy_0,  \
                             ta1_x_yyy_yyy_1,  \
                             ta1_x_yyy_yyz_0,  \
                             ta1_x_yyy_yyz_1,  \
                             ta1_x_yyy_yz_0,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_yzz_0,  \
                             ta1_x_yyy_yzz_1,  \
                             ta1_x_yyyz_xxx_0, \
                             ta1_x_yyyz_xxy_0, \
                             ta1_x_yyyz_xxz_0, \
                             ta1_x_yyyz_xyy_0, \
                             ta1_x_yyyz_xyz_0, \
                             ta1_x_yyyz_xzz_0, \
                             ta1_x_yyyz_yyy_0, \
                             ta1_x_yyyz_yyz_0, \
                             ta1_x_yyyz_yzz_0, \
                             ta1_x_yyyz_zzz_0, \
                             ta1_x_yyz_xxz_0,  \
                             ta1_x_yyz_xxz_1,  \
                             ta1_x_yyz_xzz_0,  \
                             ta1_x_yyz_xzz_1,  \
                             ta1_x_yyz_zzz_0,  \
                             ta1_x_yyz_zzz_1,  \
                             ta1_x_yz_xxz_0,   \
                             ta1_x_yz_xxz_1,   \
                             ta1_x_yz_xzz_0,   \
                             ta1_x_yz_xzz_1,   \
                             ta1_x_yz_zzz_0,   \
                             ta1_x_yz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_xxx_0[i] = ta1_x_yyy_xxx_0[i] * pa_z[i] - ta1_x_yyy_xxx_1[i] * pc_z[i];

        ta1_x_yyyz_xxy_0[i] = ta1_x_yyy_xxy_0[i] * pa_z[i] - ta1_x_yyy_xxy_1[i] * pc_z[i];

        ta1_x_yyyz_xxz_0[i] =
            2.0 * ta1_x_yz_xxz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxz_1[i] * fe_0 + ta1_x_yyz_xxz_0[i] * pa_y[i] - ta1_x_yyz_xxz_1[i] * pc_y[i];

        ta1_x_yyyz_xyy_0[i] = ta1_x_yyy_xyy_0[i] * pa_z[i] - ta1_x_yyy_xyy_1[i] * pc_z[i];

        ta1_x_yyyz_xyz_0[i] = ta1_x_yyy_xy_0[i] * fe_0 - ta1_x_yyy_xy_1[i] * fe_0 + ta1_x_yyy_xyz_0[i] * pa_z[i] - ta1_x_yyy_xyz_1[i] * pc_z[i];

        ta1_x_yyyz_xzz_0[i] =
            2.0 * ta1_x_yz_xzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xzz_1[i] * fe_0 + ta1_x_yyz_xzz_0[i] * pa_y[i] - ta1_x_yyz_xzz_1[i] * pc_y[i];

        ta1_x_yyyz_yyy_0[i] = ta1_x_yyy_yyy_0[i] * pa_z[i] - ta1_x_yyy_yyy_1[i] * pc_z[i];

        ta1_x_yyyz_yyz_0[i] = ta1_x_yyy_yy_0[i] * fe_0 - ta1_x_yyy_yy_1[i] * fe_0 + ta1_x_yyy_yyz_0[i] * pa_z[i] - ta1_x_yyy_yyz_1[i] * pc_z[i];

        ta1_x_yyyz_yzz_0[i] =
            2.0 * ta1_x_yyy_yz_0[i] * fe_0 - 2.0 * ta1_x_yyy_yz_1[i] * fe_0 + ta1_x_yyy_yzz_0[i] * pa_z[i] - ta1_x_yyy_yzz_1[i] * pc_z[i];

        ta1_x_yyyz_zzz_0[i] =
            2.0 * ta1_x_yz_zzz_0[i] * fe_0 - 2.0 * ta1_x_yz_zzz_1[i] * fe_0 + ta1_x_yyz_zzz_0[i] * pa_y[i] - ta1_x_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto ta1_x_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 120);

    auto ta1_x_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 121);

    auto ta1_x_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 122);

    auto ta1_x_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 123);

    auto ta1_x_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 124);

    auto ta1_x_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 125);

    auto ta1_x_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 126);

    auto ta1_x_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 127);

    auto ta1_x_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 128);

    auto ta1_x_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 129);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yy_xxy_0,   \
                             ta1_x_yy_xxy_1,   \
                             ta1_x_yy_xyy_0,   \
                             ta1_x_yy_xyy_1,   \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yyz_xxy_0,  \
                             ta1_x_yyz_xxy_1,  \
                             ta1_x_yyz_xyy_0,  \
                             ta1_x_yyz_xyy_1,  \
                             ta1_x_yyz_yyy_0,  \
                             ta1_x_yyz_yyy_1,  \
                             ta1_x_yyzz_xxx_0, \
                             ta1_x_yyzz_xxy_0, \
                             ta1_x_yyzz_xxz_0, \
                             ta1_x_yyzz_xyy_0, \
                             ta1_x_yyzz_xyz_0, \
                             ta1_x_yyzz_xzz_0, \
                             ta1_x_yyzz_yyy_0, \
                             ta1_x_yyzz_yyz_0, \
                             ta1_x_yyzz_yzz_0, \
                             ta1_x_yyzz_zzz_0, \
                             ta1_x_yzz_xxx_0,  \
                             ta1_x_yzz_xxx_1,  \
                             ta1_x_yzz_xxz_0,  \
                             ta1_x_yzz_xxz_1,  \
                             ta1_x_yzz_xyz_0,  \
                             ta1_x_yzz_xyz_1,  \
                             ta1_x_yzz_xz_0,   \
                             ta1_x_yzz_xz_1,   \
                             ta1_x_yzz_xzz_0,  \
                             ta1_x_yzz_xzz_1,  \
                             ta1_x_yzz_yyz_0,  \
                             ta1_x_yzz_yyz_1,  \
                             ta1_x_yzz_yz_0,   \
                             ta1_x_yzz_yz_1,   \
                             ta1_x_yzz_yzz_0,  \
                             ta1_x_yzz_yzz_1,  \
                             ta1_x_yzz_zz_0,   \
                             ta1_x_yzz_zz_1,   \
                             ta1_x_yzz_zzz_0,  \
                             ta1_x_yzz_zzz_1,  \
                             ta1_x_zz_xxx_0,   \
                             ta1_x_zz_xxx_1,   \
                             ta1_x_zz_xxz_0,   \
                             ta1_x_zz_xxz_1,   \
                             ta1_x_zz_xyz_0,   \
                             ta1_x_zz_xyz_1,   \
                             ta1_x_zz_xzz_0,   \
                             ta1_x_zz_xzz_1,   \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_xxx_0[i] = ta1_x_zz_xxx_0[i] * fe_0 - ta1_x_zz_xxx_1[i] * fe_0 + ta1_x_yzz_xxx_0[i] * pa_y[i] - ta1_x_yzz_xxx_1[i] * pc_y[i];

        ta1_x_yyzz_xxy_0[i] = ta1_x_yy_xxy_0[i] * fe_0 - ta1_x_yy_xxy_1[i] * fe_0 + ta1_x_yyz_xxy_0[i] * pa_z[i] - ta1_x_yyz_xxy_1[i] * pc_z[i];

        ta1_x_yyzz_xxz_0[i] = ta1_x_zz_xxz_0[i] * fe_0 - ta1_x_zz_xxz_1[i] * fe_0 + ta1_x_yzz_xxz_0[i] * pa_y[i] - ta1_x_yzz_xxz_1[i] * pc_y[i];

        ta1_x_yyzz_xyy_0[i] = ta1_x_yy_xyy_0[i] * fe_0 - ta1_x_yy_xyy_1[i] * fe_0 + ta1_x_yyz_xyy_0[i] * pa_z[i] - ta1_x_yyz_xyy_1[i] * pc_z[i];

        ta1_x_yyzz_xyz_0[i] = ta1_x_zz_xyz_0[i] * fe_0 - ta1_x_zz_xyz_1[i] * fe_0 + ta1_x_yzz_xz_0[i] * fe_0 - ta1_x_yzz_xz_1[i] * fe_0 +
                              ta1_x_yzz_xyz_0[i] * pa_y[i] - ta1_x_yzz_xyz_1[i] * pc_y[i];

        ta1_x_yyzz_xzz_0[i] = ta1_x_zz_xzz_0[i] * fe_0 - ta1_x_zz_xzz_1[i] * fe_0 + ta1_x_yzz_xzz_0[i] * pa_y[i] - ta1_x_yzz_xzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyy_0[i] = ta1_x_yy_yyy_0[i] * fe_0 - ta1_x_yy_yyy_1[i] * fe_0 + ta1_x_yyz_yyy_0[i] * pa_z[i] - ta1_x_yyz_yyy_1[i] * pc_z[i];

        ta1_x_yyzz_yyz_0[i] = ta1_x_zz_yyz_0[i] * fe_0 - ta1_x_zz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yzz_yz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yz_1[i] * fe_0 +
                              ta1_x_yzz_yyz_0[i] * pa_y[i] - ta1_x_yzz_yyz_1[i] * pc_y[i];

        ta1_x_yyzz_yzz_0[i] = ta1_x_zz_yzz_0[i] * fe_0 - ta1_x_zz_yzz_1[i] * fe_0 + ta1_x_yzz_zz_0[i] * fe_0 - ta1_x_yzz_zz_1[i] * fe_0 +
                              ta1_x_yzz_yzz_0[i] * pa_y[i] - ta1_x_yzz_yzz_1[i] * pc_y[i];

        ta1_x_yyzz_zzz_0[i] = ta1_x_zz_zzz_0[i] * fe_0 - ta1_x_zz_zzz_1[i] * fe_0 + ta1_x_yzz_zzz_0[i] * pa_y[i] - ta1_x_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto ta1_x_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 130);

    auto ta1_x_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 131);

    auto ta1_x_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 132);

    auto ta1_x_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 133);

    auto ta1_x_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 134);

    auto ta1_x_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 135);

    auto ta1_x_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 136);

    auto ta1_x_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 137);

    auto ta1_x_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 138);

    auto ta1_x_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 139);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yzzz_xxx_0, \
                             ta1_x_yzzz_xxy_0, \
                             ta1_x_yzzz_xxz_0, \
                             ta1_x_yzzz_xyy_0, \
                             ta1_x_yzzz_xyz_0, \
                             ta1_x_yzzz_xzz_0, \
                             ta1_x_yzzz_yyy_0, \
                             ta1_x_yzzz_yyz_0, \
                             ta1_x_yzzz_yzz_0, \
                             ta1_x_yzzz_zzz_0, \
                             ta1_x_zzz_xx_0,   \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xxx_0,  \
                             ta1_x_zzz_xxx_1,  \
                             ta1_x_zzz_xxy_0,  \
                             ta1_x_zzz_xxy_1,  \
                             ta1_x_zzz_xxz_0,  \
                             ta1_x_zzz_xxz_1,  \
                             ta1_x_zzz_xy_0,   \
                             ta1_x_zzz_xy_1,   \
                             ta1_x_zzz_xyy_0,  \
                             ta1_x_zzz_xyy_1,  \
                             ta1_x_zzz_xyz_0,  \
                             ta1_x_zzz_xyz_1,  \
                             ta1_x_zzz_xz_0,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_xzz_0,  \
                             ta1_x_zzz_xzz_1,  \
                             ta1_x_zzz_yy_0,   \
                             ta1_x_zzz_yy_1,   \
                             ta1_x_zzz_yyy_0,  \
                             ta1_x_zzz_yyy_1,  \
                             ta1_x_zzz_yyz_0,  \
                             ta1_x_zzz_yyz_1,  \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_yzz_0,  \
                             ta1_x_zzz_yzz_1,  \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             ta1_x_zzz_zzz_0,  \
                             ta1_x_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_xxx_0[i] = ta1_x_zzz_xxx_0[i] * pa_y[i] - ta1_x_zzz_xxx_1[i] * pc_y[i];

        ta1_x_yzzz_xxy_0[i] = ta1_x_zzz_xx_0[i] * fe_0 - ta1_x_zzz_xx_1[i] * fe_0 + ta1_x_zzz_xxy_0[i] * pa_y[i] - ta1_x_zzz_xxy_1[i] * pc_y[i];

        ta1_x_yzzz_xxz_0[i] = ta1_x_zzz_xxz_0[i] * pa_y[i] - ta1_x_zzz_xxz_1[i] * pc_y[i];

        ta1_x_yzzz_xyy_0[i] =
            2.0 * ta1_x_zzz_xy_0[i] * fe_0 - 2.0 * ta1_x_zzz_xy_1[i] * fe_0 + ta1_x_zzz_xyy_0[i] * pa_y[i] - ta1_x_zzz_xyy_1[i] * pc_y[i];

        ta1_x_yzzz_xyz_0[i] = ta1_x_zzz_xz_0[i] * fe_0 - ta1_x_zzz_xz_1[i] * fe_0 + ta1_x_zzz_xyz_0[i] * pa_y[i] - ta1_x_zzz_xyz_1[i] * pc_y[i];

        ta1_x_yzzz_xzz_0[i] = ta1_x_zzz_xzz_0[i] * pa_y[i] - ta1_x_zzz_xzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyy_0[i] =
            3.0 * ta1_x_zzz_yy_0[i] * fe_0 - 3.0 * ta1_x_zzz_yy_1[i] * fe_0 + ta1_x_zzz_yyy_0[i] * pa_y[i] - ta1_x_zzz_yyy_1[i] * pc_y[i];

        ta1_x_yzzz_yyz_0[i] =
            2.0 * ta1_x_zzz_yz_0[i] * fe_0 - 2.0 * ta1_x_zzz_yz_1[i] * fe_0 + ta1_x_zzz_yyz_0[i] * pa_y[i] - ta1_x_zzz_yyz_1[i] * pc_y[i];

        ta1_x_yzzz_yzz_0[i] = ta1_x_zzz_zz_0[i] * fe_0 - ta1_x_zzz_zz_1[i] * fe_0 + ta1_x_zzz_yzz_0[i] * pa_y[i] - ta1_x_zzz_yzz_1[i] * pc_y[i];

        ta1_x_yzzz_zzz_0[i] = ta1_x_zzz_zzz_0[i] * pa_y[i] - ta1_x_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto ta1_x_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 140);

    auto ta1_x_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 141);

    auto ta1_x_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 142);

    auto ta1_x_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 143);

    auto ta1_x_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 144);

    auto ta1_x_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 145);

    auto ta1_x_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 146);

    auto ta1_x_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 147);

    auto ta1_x_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 148);

    auto ta1_x_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 149);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_zz_xxx_0,   \
                             ta1_x_zz_xxx_1,   \
                             ta1_x_zz_xxy_0,   \
                             ta1_x_zz_xxy_1,   \
                             ta1_x_zz_xxz_0,   \
                             ta1_x_zz_xxz_1,   \
                             ta1_x_zz_xyy_0,   \
                             ta1_x_zz_xyy_1,   \
                             ta1_x_zz_xyz_0,   \
                             ta1_x_zz_xyz_1,   \
                             ta1_x_zz_xzz_0,   \
                             ta1_x_zz_xzz_1,   \
                             ta1_x_zz_yyy_0,   \
                             ta1_x_zz_yyy_1,   \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             ta1_x_zzz_xx_0,   \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xxx_0,  \
                             ta1_x_zzz_xxx_1,  \
                             ta1_x_zzz_xxy_0,  \
                             ta1_x_zzz_xxy_1,  \
                             ta1_x_zzz_xxz_0,  \
                             ta1_x_zzz_xxz_1,  \
                             ta1_x_zzz_xy_0,   \
                             ta1_x_zzz_xy_1,   \
                             ta1_x_zzz_xyy_0,  \
                             ta1_x_zzz_xyy_1,  \
                             ta1_x_zzz_xyz_0,  \
                             ta1_x_zzz_xyz_1,  \
                             ta1_x_zzz_xz_0,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_xzz_0,  \
                             ta1_x_zzz_xzz_1,  \
                             ta1_x_zzz_yy_0,   \
                             ta1_x_zzz_yy_1,   \
                             ta1_x_zzz_yyy_0,  \
                             ta1_x_zzz_yyy_1,  \
                             ta1_x_zzz_yyz_0,  \
                             ta1_x_zzz_yyz_1,  \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_yzz_0,  \
                             ta1_x_zzz_yzz_1,  \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             ta1_x_zzz_zzz_0,  \
                             ta1_x_zzz_zzz_1,  \
                             ta1_x_zzzz_xxx_0, \
                             ta1_x_zzzz_xxy_0, \
                             ta1_x_zzzz_xxz_0, \
                             ta1_x_zzzz_xyy_0, \
                             ta1_x_zzzz_xyz_0, \
                             ta1_x_zzzz_xzz_0, \
                             ta1_x_zzzz_yyy_0, \
                             ta1_x_zzzz_yyz_0, \
                             ta1_x_zzzz_yzz_0, \
                             ta1_x_zzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_xxx_0[i] =
            3.0 * ta1_x_zz_xxx_0[i] * fe_0 - 3.0 * ta1_x_zz_xxx_1[i] * fe_0 + ta1_x_zzz_xxx_0[i] * pa_z[i] - ta1_x_zzz_xxx_1[i] * pc_z[i];

        ta1_x_zzzz_xxy_0[i] =
            3.0 * ta1_x_zz_xxy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxy_1[i] * fe_0 + ta1_x_zzz_xxy_0[i] * pa_z[i] - ta1_x_zzz_xxy_1[i] * pc_z[i];

        ta1_x_zzzz_xxz_0[i] = 3.0 * ta1_x_zz_xxz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxz_1[i] * fe_0 + ta1_x_zzz_xx_0[i] * fe_0 - ta1_x_zzz_xx_1[i] * fe_0 +
                              ta1_x_zzz_xxz_0[i] * pa_z[i] - ta1_x_zzz_xxz_1[i] * pc_z[i];

        ta1_x_zzzz_xyy_0[i] =
            3.0 * ta1_x_zz_xyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xyy_1[i] * fe_0 + ta1_x_zzz_xyy_0[i] * pa_z[i] - ta1_x_zzz_xyy_1[i] * pc_z[i];

        ta1_x_zzzz_xyz_0[i] = 3.0 * ta1_x_zz_xyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyz_1[i] * fe_0 + ta1_x_zzz_xy_0[i] * fe_0 - ta1_x_zzz_xy_1[i] * fe_0 +
                              ta1_x_zzz_xyz_0[i] * pa_z[i] - ta1_x_zzz_xyz_1[i] * pc_z[i];

        ta1_x_zzzz_xzz_0[i] = 3.0 * ta1_x_zz_xzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xz_0[i] * fe_0 -
                              2.0 * ta1_x_zzz_xz_1[i] * fe_0 + ta1_x_zzz_xzz_0[i] * pa_z[i] - ta1_x_zzz_xzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyy_0[i] =
            3.0 * ta1_x_zz_yyy_0[i] * fe_0 - 3.0 * ta1_x_zz_yyy_1[i] * fe_0 + ta1_x_zzz_yyy_0[i] * pa_z[i] - ta1_x_zzz_yyy_1[i] * pc_z[i];

        ta1_x_zzzz_yyz_0[i] = 3.0 * ta1_x_zz_yyz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyz_1[i] * fe_0 + ta1_x_zzz_yy_0[i] * fe_0 - ta1_x_zzz_yy_1[i] * fe_0 +
                              ta1_x_zzz_yyz_0[i] * pa_z[i] - ta1_x_zzz_yyz_1[i] * pc_z[i];

        ta1_x_zzzz_yzz_0[i] = 3.0 * ta1_x_zz_yzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_yz_0[i] * fe_0 -
                              2.0 * ta1_x_zzz_yz_1[i] * fe_0 + ta1_x_zzz_yzz_0[i] * pa_z[i] - ta1_x_zzz_yzz_1[i] * pc_z[i];

        ta1_x_zzzz_zzz_0[i] = 3.0 * ta1_x_zz_zzz_0[i] * fe_0 - 3.0 * ta1_x_zz_zzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_zz_0[i] * fe_0 -
                              3.0 * ta1_x_zzz_zz_1[i] * fe_0 + ta1_x_zzz_zzz_0[i] * pa_z[i] - ta1_x_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 150-160 components of targeted buffer : GF

    auto ta1_y_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 150);

    auto ta1_y_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 151);

    auto ta1_y_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 152);

    auto ta1_y_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 153);

    auto ta1_y_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 154);

    auto ta1_y_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 155);

    auto ta1_y_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 156);

    auto ta1_y_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 157);

    auto ta1_y_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 158);

    auto ta1_y_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 159);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxy_0,   \
                             ta1_y_xx_xxy_1,   \
                             ta1_y_xx_xxz_0,   \
                             ta1_y_xx_xxz_1,   \
                             ta1_y_xx_xyy_0,   \
                             ta1_y_xx_xyy_1,   \
                             ta1_y_xx_xyz_0,   \
                             ta1_y_xx_xyz_1,   \
                             ta1_y_xx_xzz_0,   \
                             ta1_y_xx_xzz_1,   \
                             ta1_y_xx_yyy_0,   \
                             ta1_y_xx_yyy_1,   \
                             ta1_y_xx_yyz_0,   \
                             ta1_y_xx_yyz_1,   \
                             ta1_y_xx_yzz_0,   \
                             ta1_y_xx_yzz_1,   \
                             ta1_y_xx_zzz_0,   \
                             ta1_y_xx_zzz_1,   \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xxx_0,  \
                             ta1_y_xxx_xxx_1,  \
                             ta1_y_xxx_xxy_0,  \
                             ta1_y_xxx_xxy_1,  \
                             ta1_y_xxx_xxz_0,  \
                             ta1_y_xxx_xxz_1,  \
                             ta1_y_xxx_xy_0,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xyy_0,  \
                             ta1_y_xxx_xyy_1,  \
                             ta1_y_xxx_xyz_0,  \
                             ta1_y_xxx_xyz_1,  \
                             ta1_y_xxx_xz_0,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_xzz_0,  \
                             ta1_y_xxx_xzz_1,  \
                             ta1_y_xxx_yy_0,   \
                             ta1_y_xxx_yy_1,   \
                             ta1_y_xxx_yyy_0,  \
                             ta1_y_xxx_yyy_1,  \
                             ta1_y_xxx_yyz_0,  \
                             ta1_y_xxx_yyz_1,  \
                             ta1_y_xxx_yz_0,   \
                             ta1_y_xxx_yz_1,   \
                             ta1_y_xxx_yzz_0,  \
                             ta1_y_xxx_yzz_1,  \
                             ta1_y_xxx_zz_0,   \
                             ta1_y_xxx_zz_1,   \
                             ta1_y_xxx_zzz_0,  \
                             ta1_y_xxx_zzz_1,  \
                             ta1_y_xxxx_xxx_0, \
                             ta1_y_xxxx_xxy_0, \
                             ta1_y_xxxx_xxz_0, \
                             ta1_y_xxxx_xyy_0, \
                             ta1_y_xxxx_xyz_0, \
                             ta1_y_xxxx_xzz_0, \
                             ta1_y_xxxx_yyy_0, \
                             ta1_y_xxxx_yyz_0, \
                             ta1_y_xxxx_yzz_0, \
                             ta1_y_xxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_xxx_0[i] = 3.0 * ta1_y_xx_xxx_0[i] * fe_0 - 3.0 * ta1_y_xx_xxx_1[i] * fe_0 + 3.0 * ta1_y_xxx_xx_0[i] * fe_0 -
                              3.0 * ta1_y_xxx_xx_1[i] * fe_0 + ta1_y_xxx_xxx_0[i] * pa_x[i] - ta1_y_xxx_xxx_1[i] * pc_x[i];

        ta1_y_xxxx_xxy_0[i] = 3.0 * ta1_y_xx_xxy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xy_0[i] * fe_0 -
                              2.0 * ta1_y_xxx_xy_1[i] * fe_0 + ta1_y_xxx_xxy_0[i] * pa_x[i] - ta1_y_xxx_xxy_1[i] * pc_x[i];

        ta1_y_xxxx_xxz_0[i] = 3.0 * ta1_y_xx_xxz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xz_0[i] * fe_0 -
                              2.0 * ta1_y_xxx_xz_1[i] * fe_0 + ta1_y_xxx_xxz_0[i] * pa_x[i] - ta1_y_xxx_xxz_1[i] * pc_x[i];

        ta1_y_xxxx_xyy_0[i] = 3.0 * ta1_y_xx_xyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xyy_1[i] * fe_0 + ta1_y_xxx_yy_0[i] * fe_0 - ta1_y_xxx_yy_1[i] * fe_0 +
                              ta1_y_xxx_xyy_0[i] * pa_x[i] - ta1_y_xxx_xyy_1[i] * pc_x[i];

        ta1_y_xxxx_xyz_0[i] = 3.0 * ta1_y_xx_xyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyz_1[i] * fe_0 + ta1_y_xxx_yz_0[i] * fe_0 - ta1_y_xxx_yz_1[i] * fe_0 +
                              ta1_y_xxx_xyz_0[i] * pa_x[i] - ta1_y_xxx_xyz_1[i] * pc_x[i];

        ta1_y_xxxx_xzz_0[i] = 3.0 * ta1_y_xx_xzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xzz_1[i] * fe_0 + ta1_y_xxx_zz_0[i] * fe_0 - ta1_y_xxx_zz_1[i] * fe_0 +
                              ta1_y_xxx_xzz_0[i] * pa_x[i] - ta1_y_xxx_xzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyy_0[i] =
            3.0 * ta1_y_xx_yyy_0[i] * fe_0 - 3.0 * ta1_y_xx_yyy_1[i] * fe_0 + ta1_y_xxx_yyy_0[i] * pa_x[i] - ta1_y_xxx_yyy_1[i] * pc_x[i];

        ta1_y_xxxx_yyz_0[i] =
            3.0 * ta1_y_xx_yyz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyz_1[i] * fe_0 + ta1_y_xxx_yyz_0[i] * pa_x[i] - ta1_y_xxx_yyz_1[i] * pc_x[i];

        ta1_y_xxxx_yzz_0[i] =
            3.0 * ta1_y_xx_yzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yzz_1[i] * fe_0 + ta1_y_xxx_yzz_0[i] * pa_x[i] - ta1_y_xxx_yzz_1[i] * pc_x[i];

        ta1_y_xxxx_zzz_0[i] =
            3.0 * ta1_y_xx_zzz_0[i] * fe_0 - 3.0 * ta1_y_xx_zzz_1[i] * fe_0 + ta1_y_xxx_zzz_0[i] * pa_x[i] - ta1_y_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : GF

    auto ta1_y_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 160);

    auto ta1_y_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 161);

    auto ta1_y_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 162);

    auto ta1_y_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 163);

    auto ta1_y_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 164);

    auto ta1_y_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 165);

    auto ta1_y_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 166);

    auto ta1_y_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 167);

    auto ta1_y_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 168);

    auto ta1_y_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 169);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xxx_0,  \
                             ta1_y_xxx_xxx_1,  \
                             ta1_y_xxx_xxy_0,  \
                             ta1_y_xxx_xxy_1,  \
                             ta1_y_xxx_xxz_0,  \
                             ta1_y_xxx_xxz_1,  \
                             ta1_y_xxx_xy_0,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xyy_0,  \
                             ta1_y_xxx_xyy_1,  \
                             ta1_y_xxx_xyz_0,  \
                             ta1_y_xxx_xyz_1,  \
                             ta1_y_xxx_xz_0,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_xzz_0,  \
                             ta1_y_xxx_xzz_1,  \
                             ta1_y_xxx_zzz_0,  \
                             ta1_y_xxx_zzz_1,  \
                             ta1_y_xxxy_xxx_0, \
                             ta1_y_xxxy_xxy_0, \
                             ta1_y_xxxy_xxz_0, \
                             ta1_y_xxxy_xyy_0, \
                             ta1_y_xxxy_xyz_0, \
                             ta1_y_xxxy_xzz_0, \
                             ta1_y_xxxy_yyy_0, \
                             ta1_y_xxxy_yyz_0, \
                             ta1_y_xxxy_yzz_0, \
                             ta1_y_xxxy_zzz_0, \
                             ta1_y_xxy_yyy_0,  \
                             ta1_y_xxy_yyy_1,  \
                             ta1_y_xxy_yyz_0,  \
                             ta1_y_xxy_yyz_1,  \
                             ta1_y_xxy_yzz_0,  \
                             ta1_y_xxy_yzz_1,  \
                             ta1_y_xy_yyy_0,   \
                             ta1_y_xy_yyy_1,   \
                             ta1_y_xy_yyz_0,   \
                             ta1_y_xy_yyz_1,   \
                             ta1_y_xy_yzz_0,   \
                             ta1_y_xy_yzz_1,   \
                             ta_xxx_xxx_1,     \
                             ta_xxx_xxy_1,     \
                             ta_xxx_xxz_1,     \
                             ta_xxx_xyy_1,     \
                             ta_xxx_xyz_1,     \
                             ta_xxx_xzz_1,     \
                             ta_xxx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_xxx_0[i] = ta_xxx_xxx_1[i] + ta1_y_xxx_xxx_0[i] * pa_y[i] - ta1_y_xxx_xxx_1[i] * pc_y[i];

        ta1_y_xxxy_xxy_0[i] =
            ta1_y_xxx_xx_0[i] * fe_0 - ta1_y_xxx_xx_1[i] * fe_0 + ta_xxx_xxy_1[i] + ta1_y_xxx_xxy_0[i] * pa_y[i] - ta1_y_xxx_xxy_1[i] * pc_y[i];

        ta1_y_xxxy_xxz_0[i] = ta_xxx_xxz_1[i] + ta1_y_xxx_xxz_0[i] * pa_y[i] - ta1_y_xxx_xxz_1[i] * pc_y[i];

        ta1_y_xxxy_xyy_0[i] = 2.0 * ta1_y_xxx_xy_0[i] * fe_0 - 2.0 * ta1_y_xxx_xy_1[i] * fe_0 + ta_xxx_xyy_1[i] + ta1_y_xxx_xyy_0[i] * pa_y[i] -
                              ta1_y_xxx_xyy_1[i] * pc_y[i];

        ta1_y_xxxy_xyz_0[i] =
            ta1_y_xxx_xz_0[i] * fe_0 - ta1_y_xxx_xz_1[i] * fe_0 + ta_xxx_xyz_1[i] + ta1_y_xxx_xyz_0[i] * pa_y[i] - ta1_y_xxx_xyz_1[i] * pc_y[i];

        ta1_y_xxxy_xzz_0[i] = ta_xxx_xzz_1[i] + ta1_y_xxx_xzz_0[i] * pa_y[i] - ta1_y_xxx_xzz_1[i] * pc_y[i];

        ta1_y_xxxy_yyy_0[i] =
            2.0 * ta1_y_xy_yyy_0[i] * fe_0 - 2.0 * ta1_y_xy_yyy_1[i] * fe_0 + ta1_y_xxy_yyy_0[i] * pa_x[i] - ta1_y_xxy_yyy_1[i] * pc_x[i];

        ta1_y_xxxy_yyz_0[i] =
            2.0 * ta1_y_xy_yyz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyz_1[i] * fe_0 + ta1_y_xxy_yyz_0[i] * pa_x[i] - ta1_y_xxy_yyz_1[i] * pc_x[i];

        ta1_y_xxxy_yzz_0[i] =
            2.0 * ta1_y_xy_yzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yzz_1[i] * fe_0 + ta1_y_xxy_yzz_0[i] * pa_x[i] - ta1_y_xxy_yzz_1[i] * pc_x[i];

        ta1_y_xxxy_zzz_0[i] = ta_xxx_zzz_1[i] + ta1_y_xxx_zzz_0[i] * pa_y[i] - ta1_y_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : GF

    auto ta1_y_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 170);

    auto ta1_y_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 171);

    auto ta1_y_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 172);

    auto ta1_y_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 173);

    auto ta1_y_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 174);

    auto ta1_y_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 175);

    auto ta1_y_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 176);

    auto ta1_y_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 177);

    auto ta1_y_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 178);

    auto ta1_y_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 179);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xxx_0,  \
                             ta1_y_xxx_xxx_1,  \
                             ta1_y_xxx_xxy_0,  \
                             ta1_y_xxx_xxy_1,  \
                             ta1_y_xxx_xxz_0,  \
                             ta1_y_xxx_xxz_1,  \
                             ta1_y_xxx_xy_0,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xyy_0,  \
                             ta1_y_xxx_xyy_1,  \
                             ta1_y_xxx_xyz_0,  \
                             ta1_y_xxx_xyz_1,  \
                             ta1_y_xxx_xz_0,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_xzz_0,  \
                             ta1_y_xxx_xzz_1,  \
                             ta1_y_xxx_yyy_0,  \
                             ta1_y_xxx_yyy_1,  \
                             ta1_y_xxxz_xxx_0, \
                             ta1_y_xxxz_xxy_0, \
                             ta1_y_xxxz_xxz_0, \
                             ta1_y_xxxz_xyy_0, \
                             ta1_y_xxxz_xyz_0, \
                             ta1_y_xxxz_xzz_0, \
                             ta1_y_xxxz_yyy_0, \
                             ta1_y_xxxz_yyz_0, \
                             ta1_y_xxxz_yzz_0, \
                             ta1_y_xxxz_zzz_0, \
                             ta1_y_xxz_yyz_0,  \
                             ta1_y_xxz_yyz_1,  \
                             ta1_y_xxz_yzz_0,  \
                             ta1_y_xxz_yzz_1,  \
                             ta1_y_xxz_zzz_0,  \
                             ta1_y_xxz_zzz_1,  \
                             ta1_y_xz_yyz_0,   \
                             ta1_y_xz_yyz_1,   \
                             ta1_y_xz_yzz_0,   \
                             ta1_y_xz_yzz_1,   \
                             ta1_y_xz_zzz_0,   \
                             ta1_y_xz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_xxx_0[i] = ta1_y_xxx_xxx_0[i] * pa_z[i] - ta1_y_xxx_xxx_1[i] * pc_z[i];

        ta1_y_xxxz_xxy_0[i] = ta1_y_xxx_xxy_0[i] * pa_z[i] - ta1_y_xxx_xxy_1[i] * pc_z[i];

        ta1_y_xxxz_xxz_0[i] = ta1_y_xxx_xx_0[i] * fe_0 - ta1_y_xxx_xx_1[i] * fe_0 + ta1_y_xxx_xxz_0[i] * pa_z[i] - ta1_y_xxx_xxz_1[i] * pc_z[i];

        ta1_y_xxxz_xyy_0[i] = ta1_y_xxx_xyy_0[i] * pa_z[i] - ta1_y_xxx_xyy_1[i] * pc_z[i];

        ta1_y_xxxz_xyz_0[i] = ta1_y_xxx_xy_0[i] * fe_0 - ta1_y_xxx_xy_1[i] * fe_0 + ta1_y_xxx_xyz_0[i] * pa_z[i] - ta1_y_xxx_xyz_1[i] * pc_z[i];

        ta1_y_xxxz_xzz_0[i] =
            2.0 * ta1_y_xxx_xz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xz_1[i] * fe_0 + ta1_y_xxx_xzz_0[i] * pa_z[i] - ta1_y_xxx_xzz_1[i] * pc_z[i];

        ta1_y_xxxz_yyy_0[i] = ta1_y_xxx_yyy_0[i] * pa_z[i] - ta1_y_xxx_yyy_1[i] * pc_z[i];

        ta1_y_xxxz_yyz_0[i] =
            2.0 * ta1_y_xz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyz_1[i] * fe_0 + ta1_y_xxz_yyz_0[i] * pa_x[i] - ta1_y_xxz_yyz_1[i] * pc_x[i];

        ta1_y_xxxz_yzz_0[i] =
            2.0 * ta1_y_xz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yzz_1[i] * fe_0 + ta1_y_xxz_yzz_0[i] * pa_x[i] - ta1_y_xxz_yzz_1[i] * pc_x[i];

        ta1_y_xxxz_zzz_0[i] =
            2.0 * ta1_y_xz_zzz_0[i] * fe_0 - 2.0 * ta1_y_xz_zzz_1[i] * fe_0 + ta1_y_xxz_zzz_0[i] * pa_x[i] - ta1_y_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 180-190 components of targeted buffer : GF

    auto ta1_y_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 180);

    auto ta1_y_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 181);

    auto ta1_y_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 182);

    auto ta1_y_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 183);

    auto ta1_y_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 184);

    auto ta1_y_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 185);

    auto ta1_y_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 186);

    auto ta1_y_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 187);

    auto ta1_y_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 188);

    auto ta1_y_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 189);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxz_0,   \
                             ta1_y_xx_xxz_1,   \
                             ta1_y_xx_xzz_0,   \
                             ta1_y_xx_xzz_1,   \
                             ta1_y_xxy_xxx_0,  \
                             ta1_y_xxy_xxx_1,  \
                             ta1_y_xxy_xxz_0,  \
                             ta1_y_xxy_xxz_1,  \
                             ta1_y_xxy_xzz_0,  \
                             ta1_y_xxy_xzz_1,  \
                             ta1_y_xxyy_xxx_0, \
                             ta1_y_xxyy_xxy_0, \
                             ta1_y_xxyy_xxz_0, \
                             ta1_y_xxyy_xyy_0, \
                             ta1_y_xxyy_xyz_0, \
                             ta1_y_xxyy_xzz_0, \
                             ta1_y_xxyy_yyy_0, \
                             ta1_y_xxyy_yyz_0, \
                             ta1_y_xxyy_yzz_0, \
                             ta1_y_xxyy_zzz_0, \
                             ta1_y_xyy_xxy_0,  \
                             ta1_y_xyy_xxy_1,  \
                             ta1_y_xyy_xy_0,   \
                             ta1_y_xyy_xy_1,   \
                             ta1_y_xyy_xyy_0,  \
                             ta1_y_xyy_xyy_1,  \
                             ta1_y_xyy_xyz_0,  \
                             ta1_y_xyy_xyz_1,  \
                             ta1_y_xyy_yy_0,   \
                             ta1_y_xyy_yy_1,   \
                             ta1_y_xyy_yyy_0,  \
                             ta1_y_xyy_yyy_1,  \
                             ta1_y_xyy_yyz_0,  \
                             ta1_y_xyy_yyz_1,  \
                             ta1_y_xyy_yz_0,   \
                             ta1_y_xyy_yz_1,   \
                             ta1_y_xyy_yzz_0,  \
                             ta1_y_xyy_yzz_1,  \
                             ta1_y_xyy_zzz_0,  \
                             ta1_y_xyy_zzz_1,  \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yy_zzz_0,   \
                             ta1_y_yy_zzz_1,   \
                             ta_xxy_xxx_1,     \
                             ta_xxy_xxz_1,     \
                             ta_xxy_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_xxx_0[i] =
            ta1_y_xx_xxx_0[i] * fe_0 - ta1_y_xx_xxx_1[i] * fe_0 + ta_xxy_xxx_1[i] + ta1_y_xxy_xxx_0[i] * pa_y[i] - ta1_y_xxy_xxx_1[i] * pc_y[i];

        ta1_y_xxyy_xxy_0[i] = ta1_y_yy_xxy_0[i] * fe_0 - ta1_y_yy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xyy_xy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xy_1[i] * fe_0 +
                              ta1_y_xyy_xxy_0[i] * pa_x[i] - ta1_y_xyy_xxy_1[i] * pc_x[i];

        ta1_y_xxyy_xxz_0[i] =
            ta1_y_xx_xxz_0[i] * fe_0 - ta1_y_xx_xxz_1[i] * fe_0 + ta_xxy_xxz_1[i] + ta1_y_xxy_xxz_0[i] * pa_y[i] - ta1_y_xxy_xxz_1[i] * pc_y[i];

        ta1_y_xxyy_xyy_0[i] = ta1_y_yy_xyy_0[i] * fe_0 - ta1_y_yy_xyy_1[i] * fe_0 + ta1_y_xyy_yy_0[i] * fe_0 - ta1_y_xyy_yy_1[i] * fe_0 +
                              ta1_y_xyy_xyy_0[i] * pa_x[i] - ta1_y_xyy_xyy_1[i] * pc_x[i];

        ta1_y_xxyy_xyz_0[i] = ta1_y_yy_xyz_0[i] * fe_0 - ta1_y_yy_xyz_1[i] * fe_0 + ta1_y_xyy_yz_0[i] * fe_0 - ta1_y_xyy_yz_1[i] * fe_0 +
                              ta1_y_xyy_xyz_0[i] * pa_x[i] - ta1_y_xyy_xyz_1[i] * pc_x[i];

        ta1_y_xxyy_xzz_0[i] =
            ta1_y_xx_xzz_0[i] * fe_0 - ta1_y_xx_xzz_1[i] * fe_0 + ta_xxy_xzz_1[i] + ta1_y_xxy_xzz_0[i] * pa_y[i] - ta1_y_xxy_xzz_1[i] * pc_y[i];

        ta1_y_xxyy_yyy_0[i] = ta1_y_yy_yyy_0[i] * fe_0 - ta1_y_yy_yyy_1[i] * fe_0 + ta1_y_xyy_yyy_0[i] * pa_x[i] - ta1_y_xyy_yyy_1[i] * pc_x[i];

        ta1_y_xxyy_yyz_0[i] = ta1_y_yy_yyz_0[i] * fe_0 - ta1_y_yy_yyz_1[i] * fe_0 + ta1_y_xyy_yyz_0[i] * pa_x[i] - ta1_y_xyy_yyz_1[i] * pc_x[i];

        ta1_y_xxyy_yzz_0[i] = ta1_y_yy_yzz_0[i] * fe_0 - ta1_y_yy_yzz_1[i] * fe_0 + ta1_y_xyy_yzz_0[i] * pa_x[i] - ta1_y_xyy_yzz_1[i] * pc_x[i];

        ta1_y_xxyy_zzz_0[i] = ta1_y_yy_zzz_0[i] * fe_0 - ta1_y_yy_zzz_1[i] * fe_0 + ta1_y_xyy_zzz_0[i] * pa_x[i] - ta1_y_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 190-200 components of targeted buffer : GF

    auto ta1_y_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 190);

    auto ta1_y_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 191);

    auto ta1_y_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 192);

    auto ta1_y_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 193);

    auto ta1_y_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 194);

    auto ta1_y_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 195);

    auto ta1_y_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 196);

    auto ta1_y_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 197);

    auto ta1_y_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 198);

    auto ta1_y_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 199);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xxy_xxx_0,  \
                             ta1_y_xxy_xxx_1,  \
                             ta1_y_xxy_xxy_0,  \
                             ta1_y_xxy_xxy_1,  \
                             ta1_y_xxy_xy_0,   \
                             ta1_y_xxy_xy_1,   \
                             ta1_y_xxy_xyy_0,  \
                             ta1_y_xxy_xyy_1,  \
                             ta1_y_xxy_xyz_0,  \
                             ta1_y_xxy_xyz_1,  \
                             ta1_y_xxy_yyy_0,  \
                             ta1_y_xxy_yyy_1,  \
                             ta1_y_xxyz_xxx_0, \
                             ta1_y_xxyz_xxy_0, \
                             ta1_y_xxyz_xxz_0, \
                             ta1_y_xxyz_xyy_0, \
                             ta1_y_xxyz_xyz_0, \
                             ta1_y_xxyz_xzz_0, \
                             ta1_y_xxyz_yyy_0, \
                             ta1_y_xxyz_yyz_0, \
                             ta1_y_xxyz_yzz_0, \
                             ta1_y_xxyz_zzz_0, \
                             ta1_y_xxz_xxz_0,  \
                             ta1_y_xxz_xxz_1,  \
                             ta1_y_xxz_xzz_0,  \
                             ta1_y_xxz_xzz_1,  \
                             ta1_y_xxz_zzz_0,  \
                             ta1_y_xxz_zzz_1,  \
                             ta1_y_xyz_yyz_0,  \
                             ta1_y_xyz_yyz_1,  \
                             ta1_y_xyz_yzz_0,  \
                             ta1_y_xyz_yzz_1,  \
                             ta1_y_yz_yyz_0,   \
                             ta1_y_yz_yyz_1,   \
                             ta1_y_yz_yzz_0,   \
                             ta1_y_yz_yzz_1,   \
                             ta_xxz_xxz_1,     \
                             ta_xxz_xzz_1,     \
                             ta_xxz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyz_xxx_0[i] = ta1_y_xxy_xxx_0[i] * pa_z[i] - ta1_y_xxy_xxx_1[i] * pc_z[i];

        ta1_y_xxyz_xxy_0[i] = ta1_y_xxy_xxy_0[i] * pa_z[i] - ta1_y_xxy_xxy_1[i] * pc_z[i];

        ta1_y_xxyz_xxz_0[i] = ta_xxz_xxz_1[i] + ta1_y_xxz_xxz_0[i] * pa_y[i] - ta1_y_xxz_xxz_1[i] * pc_y[i];

        ta1_y_xxyz_xyy_0[i] = ta1_y_xxy_xyy_0[i] * pa_z[i] - ta1_y_xxy_xyy_1[i] * pc_z[i];

        ta1_y_xxyz_xyz_0[i] = ta1_y_xxy_xy_0[i] * fe_0 - ta1_y_xxy_xy_1[i] * fe_0 + ta1_y_xxy_xyz_0[i] * pa_z[i] - ta1_y_xxy_xyz_1[i] * pc_z[i];

        ta1_y_xxyz_xzz_0[i] = ta_xxz_xzz_1[i] + ta1_y_xxz_xzz_0[i] * pa_y[i] - ta1_y_xxz_xzz_1[i] * pc_y[i];

        ta1_y_xxyz_yyy_0[i] = ta1_y_xxy_yyy_0[i] * pa_z[i] - ta1_y_xxy_yyy_1[i] * pc_z[i];

        ta1_y_xxyz_yyz_0[i] = ta1_y_yz_yyz_0[i] * fe_0 - ta1_y_yz_yyz_1[i] * fe_0 + ta1_y_xyz_yyz_0[i] * pa_x[i] - ta1_y_xyz_yyz_1[i] * pc_x[i];

        ta1_y_xxyz_yzz_0[i] = ta1_y_yz_yzz_0[i] * fe_0 - ta1_y_yz_yzz_1[i] * fe_0 + ta1_y_xyz_yzz_0[i] * pa_x[i] - ta1_y_xyz_yzz_1[i] * pc_x[i];

        ta1_y_xxyz_zzz_0[i] = ta_xxz_zzz_1[i] + ta1_y_xxz_zzz_0[i] * pa_y[i] - ta1_y_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 200-210 components of targeted buffer : GF

    auto ta1_y_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 200);

    auto ta1_y_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 201);

    auto ta1_y_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 202);

    auto ta1_y_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 203);

    auto ta1_y_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 204);

    auto ta1_y_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 205);

    auto ta1_y_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 206);

    auto ta1_y_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 207);

    auto ta1_y_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 208);

    auto ta1_y_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 209);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxy_0,   \
                             ta1_y_xx_xxy_1,   \
                             ta1_y_xx_xyy_0,   \
                             ta1_y_xx_xyy_1,   \
                             ta1_y_xxz_xxx_0,  \
                             ta1_y_xxz_xxx_1,  \
                             ta1_y_xxz_xxy_0,  \
                             ta1_y_xxz_xxy_1,  \
                             ta1_y_xxz_xyy_0,  \
                             ta1_y_xxz_xyy_1,  \
                             ta1_y_xxzz_xxx_0, \
                             ta1_y_xxzz_xxy_0, \
                             ta1_y_xxzz_xxz_0, \
                             ta1_y_xxzz_xyy_0, \
                             ta1_y_xxzz_xyz_0, \
                             ta1_y_xxzz_xzz_0, \
                             ta1_y_xxzz_yyy_0, \
                             ta1_y_xxzz_yyz_0, \
                             ta1_y_xxzz_yzz_0, \
                             ta1_y_xxzz_zzz_0, \
                             ta1_y_xzz_xxz_0,  \
                             ta1_y_xzz_xxz_1,  \
                             ta1_y_xzz_xyz_0,  \
                             ta1_y_xzz_xyz_1,  \
                             ta1_y_xzz_xz_0,   \
                             ta1_y_xzz_xz_1,   \
                             ta1_y_xzz_xzz_0,  \
                             ta1_y_xzz_xzz_1,  \
                             ta1_y_xzz_yyy_0,  \
                             ta1_y_xzz_yyy_1,  \
                             ta1_y_xzz_yyz_0,  \
                             ta1_y_xzz_yyz_1,  \
                             ta1_y_xzz_yz_0,   \
                             ta1_y_xzz_yz_1,   \
                             ta1_y_xzz_yzz_0,  \
                             ta1_y_xzz_yzz_1,  \
                             ta1_y_xzz_zz_0,   \
                             ta1_y_xzz_zz_1,   \
                             ta1_y_xzz_zzz_0,  \
                             ta1_y_xzz_zzz_1,  \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xyz_0,   \
                             ta1_y_zz_xyz_1,   \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_yyy_0,   \
                             ta1_y_zz_yyy_1,   \
                             ta1_y_zz_yyz_0,   \
                             ta1_y_zz_yyz_1,   \
                             ta1_y_zz_yzz_0,   \
                             ta1_y_zz_yzz_1,   \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_xxx_0[i] = ta1_y_xx_xxx_0[i] * fe_0 - ta1_y_xx_xxx_1[i] * fe_0 + ta1_y_xxz_xxx_0[i] * pa_z[i] - ta1_y_xxz_xxx_1[i] * pc_z[i];

        ta1_y_xxzz_xxy_0[i] = ta1_y_xx_xxy_0[i] * fe_0 - ta1_y_xx_xxy_1[i] * fe_0 + ta1_y_xxz_xxy_0[i] * pa_z[i] - ta1_y_xxz_xxy_1[i] * pc_z[i];

        ta1_y_xxzz_xxz_0[i] = ta1_y_zz_xxz_0[i] * fe_0 - ta1_y_zz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xz_1[i] * fe_0 +
                              ta1_y_xzz_xxz_0[i] * pa_x[i] - ta1_y_xzz_xxz_1[i] * pc_x[i];

        ta1_y_xxzz_xyy_0[i] = ta1_y_xx_xyy_0[i] * fe_0 - ta1_y_xx_xyy_1[i] * fe_0 + ta1_y_xxz_xyy_0[i] * pa_z[i] - ta1_y_xxz_xyy_1[i] * pc_z[i];

        ta1_y_xxzz_xyz_0[i] = ta1_y_zz_xyz_0[i] * fe_0 - ta1_y_zz_xyz_1[i] * fe_0 + ta1_y_xzz_yz_0[i] * fe_0 - ta1_y_xzz_yz_1[i] * fe_0 +
                              ta1_y_xzz_xyz_0[i] * pa_x[i] - ta1_y_xzz_xyz_1[i] * pc_x[i];

        ta1_y_xxzz_xzz_0[i] = ta1_y_zz_xzz_0[i] * fe_0 - ta1_y_zz_xzz_1[i] * fe_0 + ta1_y_xzz_zz_0[i] * fe_0 - ta1_y_xzz_zz_1[i] * fe_0 +
                              ta1_y_xzz_xzz_0[i] * pa_x[i] - ta1_y_xzz_xzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyy_0[i] = ta1_y_zz_yyy_0[i] * fe_0 - ta1_y_zz_yyy_1[i] * fe_0 + ta1_y_xzz_yyy_0[i] * pa_x[i] - ta1_y_xzz_yyy_1[i] * pc_x[i];

        ta1_y_xxzz_yyz_0[i] = ta1_y_zz_yyz_0[i] * fe_0 - ta1_y_zz_yyz_1[i] * fe_0 + ta1_y_xzz_yyz_0[i] * pa_x[i] - ta1_y_xzz_yyz_1[i] * pc_x[i];

        ta1_y_xxzz_yzz_0[i] = ta1_y_zz_yzz_0[i] * fe_0 - ta1_y_zz_yzz_1[i] * fe_0 + ta1_y_xzz_yzz_0[i] * pa_x[i] - ta1_y_xzz_yzz_1[i] * pc_x[i];

        ta1_y_xxzz_zzz_0[i] = ta1_y_zz_zzz_0[i] * fe_0 - ta1_y_zz_zzz_1[i] * fe_0 + ta1_y_xzz_zzz_0[i] * pa_x[i] - ta1_y_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : GF

    auto ta1_y_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 210);

    auto ta1_y_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 211);

    auto ta1_y_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 212);

    auto ta1_y_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 213);

    auto ta1_y_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 214);

    auto ta1_y_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 215);

    auto ta1_y_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 216);

    auto ta1_y_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 217);

    auto ta1_y_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 218);

    auto ta1_y_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 219);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xyyy_xxx_0, \
                             ta1_y_xyyy_xxy_0, \
                             ta1_y_xyyy_xxz_0, \
                             ta1_y_xyyy_xyy_0, \
                             ta1_y_xyyy_xyz_0, \
                             ta1_y_xyyy_xzz_0, \
                             ta1_y_xyyy_yyy_0, \
                             ta1_y_xyyy_yyz_0, \
                             ta1_y_xyyy_yzz_0, \
                             ta1_y_xyyy_zzz_0, \
                             ta1_y_yyy_xx_0,   \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xxx_0,  \
                             ta1_y_yyy_xxx_1,  \
                             ta1_y_yyy_xxy_0,  \
                             ta1_y_yyy_xxy_1,  \
                             ta1_y_yyy_xxz_0,  \
                             ta1_y_yyy_xxz_1,  \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xyy_0,  \
                             ta1_y_yyy_xyy_1,  \
                             ta1_y_yyy_xyz_0,  \
                             ta1_y_yyy_xyz_1,  \
                             ta1_y_yyy_xz_0,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_xzz_0,  \
                             ta1_y_yyy_xzz_1,  \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yyy_0,  \
                             ta1_y_yyy_yyy_1,  \
                             ta1_y_yyy_yyz_0,  \
                             ta1_y_yyy_yyz_1,  \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_yzz_0,  \
                             ta1_y_yyy_yzz_1,  \
                             ta1_y_yyy_zz_0,   \
                             ta1_y_yyy_zz_1,   \
                             ta1_y_yyy_zzz_0,  \
                             ta1_y_yyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_xxx_0[i] =
            3.0 * ta1_y_yyy_xx_0[i] * fe_0 - 3.0 * ta1_y_yyy_xx_1[i] * fe_0 + ta1_y_yyy_xxx_0[i] * pa_x[i] - ta1_y_yyy_xxx_1[i] * pc_x[i];

        ta1_y_xyyy_xxy_0[i] =
            2.0 * ta1_y_yyy_xy_0[i] * fe_0 - 2.0 * ta1_y_yyy_xy_1[i] * fe_0 + ta1_y_yyy_xxy_0[i] * pa_x[i] - ta1_y_yyy_xxy_1[i] * pc_x[i];

        ta1_y_xyyy_xxz_0[i] =
            2.0 * ta1_y_yyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xz_1[i] * fe_0 + ta1_y_yyy_xxz_0[i] * pa_x[i] - ta1_y_yyy_xxz_1[i] * pc_x[i];

        ta1_y_xyyy_xyy_0[i] = ta1_y_yyy_yy_0[i] * fe_0 - ta1_y_yyy_yy_1[i] * fe_0 + ta1_y_yyy_xyy_0[i] * pa_x[i] - ta1_y_yyy_xyy_1[i] * pc_x[i];

        ta1_y_xyyy_xyz_0[i] = ta1_y_yyy_yz_0[i] * fe_0 - ta1_y_yyy_yz_1[i] * fe_0 + ta1_y_yyy_xyz_0[i] * pa_x[i] - ta1_y_yyy_xyz_1[i] * pc_x[i];

        ta1_y_xyyy_xzz_0[i] = ta1_y_yyy_zz_0[i] * fe_0 - ta1_y_yyy_zz_1[i] * fe_0 + ta1_y_yyy_xzz_0[i] * pa_x[i] - ta1_y_yyy_xzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyy_0[i] = ta1_y_yyy_yyy_0[i] * pa_x[i] - ta1_y_yyy_yyy_1[i] * pc_x[i];

        ta1_y_xyyy_yyz_0[i] = ta1_y_yyy_yyz_0[i] * pa_x[i] - ta1_y_yyy_yyz_1[i] * pc_x[i];

        ta1_y_xyyy_yzz_0[i] = ta1_y_yyy_yzz_0[i] * pa_x[i] - ta1_y_yyy_yzz_1[i] * pc_x[i];

        ta1_y_xyyy_zzz_0[i] = ta1_y_yyy_zzz_0[i] * pa_x[i] - ta1_y_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 220-230 components of targeted buffer : GF

    auto ta1_y_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 220);

    auto ta1_y_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 221);

    auto ta1_y_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 222);

    auto ta1_y_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 223);

    auto ta1_y_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 224);

    auto ta1_y_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 225);

    auto ta1_y_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 226);

    auto ta1_y_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 227);

    auto ta1_y_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 228);

    auto ta1_y_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 229);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xyy_xxx_0,  \
                             ta1_y_xyy_xxx_1,  \
                             ta1_y_xyy_xxy_0,  \
                             ta1_y_xyy_xxy_1,  \
                             ta1_y_xyy_xyy_0,  \
                             ta1_y_xyy_xyy_1,  \
                             ta1_y_xyyz_xxx_0, \
                             ta1_y_xyyz_xxy_0, \
                             ta1_y_xyyz_xxz_0, \
                             ta1_y_xyyz_xyy_0, \
                             ta1_y_xyyz_xyz_0, \
                             ta1_y_xyyz_xzz_0, \
                             ta1_y_xyyz_yyy_0, \
                             ta1_y_xyyz_yyz_0, \
                             ta1_y_xyyz_yzz_0, \
                             ta1_y_xyyz_zzz_0, \
                             ta1_y_yyz_xxz_0,  \
                             ta1_y_yyz_xxz_1,  \
                             ta1_y_yyz_xyz_0,  \
                             ta1_y_yyz_xyz_1,  \
                             ta1_y_yyz_xz_0,   \
                             ta1_y_yyz_xz_1,   \
                             ta1_y_yyz_xzz_0,  \
                             ta1_y_yyz_xzz_1,  \
                             ta1_y_yyz_yyy_0,  \
                             ta1_y_yyz_yyy_1,  \
                             ta1_y_yyz_yyz_0,  \
                             ta1_y_yyz_yyz_1,  \
                             ta1_y_yyz_yz_0,   \
                             ta1_y_yyz_yz_1,   \
                             ta1_y_yyz_yzz_0,  \
                             ta1_y_yyz_yzz_1,  \
                             ta1_y_yyz_zz_0,   \
                             ta1_y_yyz_zz_1,   \
                             ta1_y_yyz_zzz_0,  \
                             ta1_y_yyz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyz_xxx_0[i] = ta1_y_xyy_xxx_0[i] * pa_z[i] - ta1_y_xyy_xxx_1[i] * pc_z[i];

        ta1_y_xyyz_xxy_0[i] = ta1_y_xyy_xxy_0[i] * pa_z[i] - ta1_y_xyy_xxy_1[i] * pc_z[i];

        ta1_y_xyyz_xxz_0[i] =
            2.0 * ta1_y_yyz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xz_1[i] * fe_0 + ta1_y_yyz_xxz_0[i] * pa_x[i] - ta1_y_yyz_xxz_1[i] * pc_x[i];

        ta1_y_xyyz_xyy_0[i] = ta1_y_xyy_xyy_0[i] * pa_z[i] - ta1_y_xyy_xyy_1[i] * pc_z[i];

        ta1_y_xyyz_xyz_0[i] = ta1_y_yyz_yz_0[i] * fe_0 - ta1_y_yyz_yz_1[i] * fe_0 + ta1_y_yyz_xyz_0[i] * pa_x[i] - ta1_y_yyz_xyz_1[i] * pc_x[i];

        ta1_y_xyyz_xzz_0[i] = ta1_y_yyz_zz_0[i] * fe_0 - ta1_y_yyz_zz_1[i] * fe_0 + ta1_y_yyz_xzz_0[i] * pa_x[i] - ta1_y_yyz_xzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyy_0[i] = ta1_y_yyz_yyy_0[i] * pa_x[i] - ta1_y_yyz_yyy_1[i] * pc_x[i];

        ta1_y_xyyz_yyz_0[i] = ta1_y_yyz_yyz_0[i] * pa_x[i] - ta1_y_yyz_yyz_1[i] * pc_x[i];

        ta1_y_xyyz_yzz_0[i] = ta1_y_yyz_yzz_0[i] * pa_x[i] - ta1_y_yyz_yzz_1[i] * pc_x[i];

        ta1_y_xyyz_zzz_0[i] = ta1_y_yyz_zzz_0[i] * pa_x[i] - ta1_y_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 230-240 components of targeted buffer : GF

    auto ta1_y_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 230);

    auto ta1_y_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 231);

    auto ta1_y_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 232);

    auto ta1_y_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 233);

    auto ta1_y_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 234);

    auto ta1_y_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 235);

    auto ta1_y_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 236);

    auto ta1_y_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 237);

    auto ta1_y_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 238);

    auto ta1_y_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 239);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xyzz_xxx_0, \
                             ta1_y_xyzz_xxy_0, \
                             ta1_y_xyzz_xxz_0, \
                             ta1_y_xyzz_xyy_0, \
                             ta1_y_xyzz_xyz_0, \
                             ta1_y_xyzz_xzz_0, \
                             ta1_y_xyzz_yyy_0, \
                             ta1_y_xyzz_yyz_0, \
                             ta1_y_xyzz_yzz_0, \
                             ta1_y_xyzz_zzz_0, \
                             ta1_y_xzz_xxx_0,  \
                             ta1_y_xzz_xxx_1,  \
                             ta1_y_xzz_xxz_0,  \
                             ta1_y_xzz_xxz_1,  \
                             ta1_y_xzz_xzz_0,  \
                             ta1_y_xzz_xzz_1,  \
                             ta1_y_yzz_xxy_0,  \
                             ta1_y_yzz_xxy_1,  \
                             ta1_y_yzz_xy_0,   \
                             ta1_y_yzz_xy_1,   \
                             ta1_y_yzz_xyy_0,  \
                             ta1_y_yzz_xyy_1,  \
                             ta1_y_yzz_xyz_0,  \
                             ta1_y_yzz_xyz_1,  \
                             ta1_y_yzz_yy_0,   \
                             ta1_y_yzz_yy_1,   \
                             ta1_y_yzz_yyy_0,  \
                             ta1_y_yzz_yyy_1,  \
                             ta1_y_yzz_yyz_0,  \
                             ta1_y_yzz_yyz_1,  \
                             ta1_y_yzz_yz_0,   \
                             ta1_y_yzz_yz_1,   \
                             ta1_y_yzz_yzz_0,  \
                             ta1_y_yzz_yzz_1,  \
                             ta1_y_yzz_zzz_0,  \
                             ta1_y_yzz_zzz_1,  \
                             ta_xzz_xxx_1,     \
                             ta_xzz_xxz_1,     \
                             ta_xzz_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzz_xxx_0[i] = ta_xzz_xxx_1[i] + ta1_y_xzz_xxx_0[i] * pa_y[i] - ta1_y_xzz_xxx_1[i] * pc_y[i];

        ta1_y_xyzz_xxy_0[i] =
            2.0 * ta1_y_yzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yzz_xy_1[i] * fe_0 + ta1_y_yzz_xxy_0[i] * pa_x[i] - ta1_y_yzz_xxy_1[i] * pc_x[i];

        ta1_y_xyzz_xxz_0[i] = ta_xzz_xxz_1[i] + ta1_y_xzz_xxz_0[i] * pa_y[i] - ta1_y_xzz_xxz_1[i] * pc_y[i];

        ta1_y_xyzz_xyy_0[i] = ta1_y_yzz_yy_0[i] * fe_0 - ta1_y_yzz_yy_1[i] * fe_0 + ta1_y_yzz_xyy_0[i] * pa_x[i] - ta1_y_yzz_xyy_1[i] * pc_x[i];

        ta1_y_xyzz_xyz_0[i] = ta1_y_yzz_yz_0[i] * fe_0 - ta1_y_yzz_yz_1[i] * fe_0 + ta1_y_yzz_xyz_0[i] * pa_x[i] - ta1_y_yzz_xyz_1[i] * pc_x[i];

        ta1_y_xyzz_xzz_0[i] = ta_xzz_xzz_1[i] + ta1_y_xzz_xzz_0[i] * pa_y[i] - ta1_y_xzz_xzz_1[i] * pc_y[i];

        ta1_y_xyzz_yyy_0[i] = ta1_y_yzz_yyy_0[i] * pa_x[i] - ta1_y_yzz_yyy_1[i] * pc_x[i];

        ta1_y_xyzz_yyz_0[i] = ta1_y_yzz_yyz_0[i] * pa_x[i] - ta1_y_yzz_yyz_1[i] * pc_x[i];

        ta1_y_xyzz_yzz_0[i] = ta1_y_yzz_yzz_0[i] * pa_x[i] - ta1_y_yzz_yzz_1[i] * pc_x[i];

        ta1_y_xyzz_zzz_0[i] = ta1_y_yzz_zzz_0[i] * pa_x[i] - ta1_y_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 240-250 components of targeted buffer : GF

    auto ta1_y_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 240);

    auto ta1_y_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 241);

    auto ta1_y_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 242);

    auto ta1_y_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 243);

    auto ta1_y_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 244);

    auto ta1_y_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 245);

    auto ta1_y_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 246);

    auto ta1_y_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 247);

    auto ta1_y_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 248);

    auto ta1_y_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 249);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xzzz_xxx_0, \
                             ta1_y_xzzz_xxy_0, \
                             ta1_y_xzzz_xxz_0, \
                             ta1_y_xzzz_xyy_0, \
                             ta1_y_xzzz_xyz_0, \
                             ta1_y_xzzz_xzz_0, \
                             ta1_y_xzzz_yyy_0, \
                             ta1_y_xzzz_yyz_0, \
                             ta1_y_xzzz_yzz_0, \
                             ta1_y_xzzz_zzz_0, \
                             ta1_y_zzz_xx_0,   \
                             ta1_y_zzz_xx_1,   \
                             ta1_y_zzz_xxx_0,  \
                             ta1_y_zzz_xxx_1,  \
                             ta1_y_zzz_xxy_0,  \
                             ta1_y_zzz_xxy_1,  \
                             ta1_y_zzz_xxz_0,  \
                             ta1_y_zzz_xxz_1,  \
                             ta1_y_zzz_xy_0,   \
                             ta1_y_zzz_xy_1,   \
                             ta1_y_zzz_xyy_0,  \
                             ta1_y_zzz_xyy_1,  \
                             ta1_y_zzz_xyz_0,  \
                             ta1_y_zzz_xyz_1,  \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_xzz_0,  \
                             ta1_y_zzz_xzz_1,  \
                             ta1_y_zzz_yy_0,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yyy_0,  \
                             ta1_y_zzz_yyy_1,  \
                             ta1_y_zzz_yyz_0,  \
                             ta1_y_zzz_yyz_1,  \
                             ta1_y_zzz_yz_0,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_yzz_0,  \
                             ta1_y_zzz_yzz_1,  \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             ta1_y_zzz_zzz_0,  \
                             ta1_y_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_xxx_0[i] =
            3.0 * ta1_y_zzz_xx_0[i] * fe_0 - 3.0 * ta1_y_zzz_xx_1[i] * fe_0 + ta1_y_zzz_xxx_0[i] * pa_x[i] - ta1_y_zzz_xxx_1[i] * pc_x[i];

        ta1_y_xzzz_xxy_0[i] =
            2.0 * ta1_y_zzz_xy_0[i] * fe_0 - 2.0 * ta1_y_zzz_xy_1[i] * fe_0 + ta1_y_zzz_xxy_0[i] * pa_x[i] - ta1_y_zzz_xxy_1[i] * pc_x[i];

        ta1_y_xzzz_xxz_0[i] =
            2.0 * ta1_y_zzz_xz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xz_1[i] * fe_0 + ta1_y_zzz_xxz_0[i] * pa_x[i] - ta1_y_zzz_xxz_1[i] * pc_x[i];

        ta1_y_xzzz_xyy_0[i] = ta1_y_zzz_yy_0[i] * fe_0 - ta1_y_zzz_yy_1[i] * fe_0 + ta1_y_zzz_xyy_0[i] * pa_x[i] - ta1_y_zzz_xyy_1[i] * pc_x[i];

        ta1_y_xzzz_xyz_0[i] = ta1_y_zzz_yz_0[i] * fe_0 - ta1_y_zzz_yz_1[i] * fe_0 + ta1_y_zzz_xyz_0[i] * pa_x[i] - ta1_y_zzz_xyz_1[i] * pc_x[i];

        ta1_y_xzzz_xzz_0[i] = ta1_y_zzz_zz_0[i] * fe_0 - ta1_y_zzz_zz_1[i] * fe_0 + ta1_y_zzz_xzz_0[i] * pa_x[i] - ta1_y_zzz_xzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyy_0[i] = ta1_y_zzz_yyy_0[i] * pa_x[i] - ta1_y_zzz_yyy_1[i] * pc_x[i];

        ta1_y_xzzz_yyz_0[i] = ta1_y_zzz_yyz_0[i] * pa_x[i] - ta1_y_zzz_yyz_1[i] * pc_x[i];

        ta1_y_xzzz_yzz_0[i] = ta1_y_zzz_yzz_0[i] * pa_x[i] - ta1_y_zzz_yzz_1[i] * pc_x[i];

        ta1_y_xzzz_zzz_0[i] = ta1_y_zzz_zzz_0[i] * pa_x[i] - ta1_y_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 250-260 components of targeted buffer : GF

    auto ta1_y_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 250);

    auto ta1_y_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 251);

    auto ta1_y_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 252);

    auto ta1_y_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 253);

    auto ta1_y_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 254);

    auto ta1_y_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 255);

    auto ta1_y_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 256);

    auto ta1_y_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 257);

    auto ta1_y_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 258);

    auto ta1_y_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 259);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_yy_xxx_0,   \
                             ta1_y_yy_xxx_1,   \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xxz_0,   \
                             ta1_y_yy_xxz_1,   \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_xzz_0,   \
                             ta1_y_yy_xzz_1,   \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yy_zzz_0,   \
                             ta1_y_yy_zzz_1,   \
                             ta1_y_yyy_xx_0,   \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xxx_0,  \
                             ta1_y_yyy_xxx_1,  \
                             ta1_y_yyy_xxy_0,  \
                             ta1_y_yyy_xxy_1,  \
                             ta1_y_yyy_xxz_0,  \
                             ta1_y_yyy_xxz_1,  \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xyy_0,  \
                             ta1_y_yyy_xyy_1,  \
                             ta1_y_yyy_xyz_0,  \
                             ta1_y_yyy_xyz_1,  \
                             ta1_y_yyy_xz_0,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_xzz_0,  \
                             ta1_y_yyy_xzz_1,  \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yyy_0,  \
                             ta1_y_yyy_yyy_1,  \
                             ta1_y_yyy_yyz_0,  \
                             ta1_y_yyy_yyz_1,  \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_yzz_0,  \
                             ta1_y_yyy_yzz_1,  \
                             ta1_y_yyy_zz_0,   \
                             ta1_y_yyy_zz_1,   \
                             ta1_y_yyy_zzz_0,  \
                             ta1_y_yyy_zzz_1,  \
                             ta1_y_yyyy_xxx_0, \
                             ta1_y_yyyy_xxy_0, \
                             ta1_y_yyyy_xxz_0, \
                             ta1_y_yyyy_xyy_0, \
                             ta1_y_yyyy_xyz_0, \
                             ta1_y_yyyy_xzz_0, \
                             ta1_y_yyyy_yyy_0, \
                             ta1_y_yyyy_yyz_0, \
                             ta1_y_yyyy_yzz_0, \
                             ta1_y_yyyy_zzz_0, \
                             ta_yyy_xxx_1,     \
                             ta_yyy_xxy_1,     \
                             ta_yyy_xxz_1,     \
                             ta_yyy_xyy_1,     \
                             ta_yyy_xyz_1,     \
                             ta_yyy_xzz_1,     \
                             ta_yyy_yyy_1,     \
                             ta_yyy_yyz_1,     \
                             ta_yyy_yzz_1,     \
                             ta_yyy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_xxx_0[i] = 3.0 * ta1_y_yy_xxx_0[i] * fe_0 - 3.0 * ta1_y_yy_xxx_1[i] * fe_0 + ta_yyy_xxx_1[i] + ta1_y_yyy_xxx_0[i] * pa_y[i] -
                              ta1_y_yyy_xxx_1[i] * pc_y[i];

        ta1_y_yyyy_xxy_0[i] = 3.0 * ta1_y_yy_xxy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxy_1[i] * fe_0 + ta1_y_yyy_xx_0[i] * fe_0 - ta1_y_yyy_xx_1[i] * fe_0 +
                              ta_yyy_xxy_1[i] + ta1_y_yyy_xxy_0[i] * pa_y[i] - ta1_y_yyy_xxy_1[i] * pc_y[i];

        ta1_y_yyyy_xxz_0[i] = 3.0 * ta1_y_yy_xxz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxz_1[i] * fe_0 + ta_yyy_xxz_1[i] + ta1_y_yyy_xxz_0[i] * pa_y[i] -
                              ta1_y_yyy_xxz_1[i] * pc_y[i];

        ta1_y_yyyy_xyy_0[i] = 3.0 * ta1_y_yy_xyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xy_0[i] * fe_0 -
                              2.0 * ta1_y_yyy_xy_1[i] * fe_0 + ta_yyy_xyy_1[i] + ta1_y_yyy_xyy_0[i] * pa_y[i] - ta1_y_yyy_xyy_1[i] * pc_y[i];

        ta1_y_yyyy_xyz_0[i] = 3.0 * ta1_y_yy_xyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyz_1[i] * fe_0 + ta1_y_yyy_xz_0[i] * fe_0 - ta1_y_yyy_xz_1[i] * fe_0 +
                              ta_yyy_xyz_1[i] + ta1_y_yyy_xyz_0[i] * pa_y[i] - ta1_y_yyy_xyz_1[i] * pc_y[i];

        ta1_y_yyyy_xzz_0[i] = 3.0 * ta1_y_yy_xzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xzz_1[i] * fe_0 + ta_yyy_xzz_1[i] + ta1_y_yyy_xzz_0[i] * pa_y[i] -
                              ta1_y_yyy_xzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyy_0[i] = 3.0 * ta1_y_yy_yyy_0[i] * fe_0 - 3.0 * ta1_y_yy_yyy_1[i] * fe_0 + 3.0 * ta1_y_yyy_yy_0[i] * fe_0 -
                              3.0 * ta1_y_yyy_yy_1[i] * fe_0 + ta_yyy_yyy_1[i] + ta1_y_yyy_yyy_0[i] * pa_y[i] - ta1_y_yyy_yyy_1[i] * pc_y[i];

        ta1_y_yyyy_yyz_0[i] = 3.0 * ta1_y_yy_yyz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yz_0[i] * fe_0 -
                              2.0 * ta1_y_yyy_yz_1[i] * fe_0 + ta_yyy_yyz_1[i] + ta1_y_yyy_yyz_0[i] * pa_y[i] - ta1_y_yyy_yyz_1[i] * pc_y[i];

        ta1_y_yyyy_yzz_0[i] = 3.0 * ta1_y_yy_yzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yzz_1[i] * fe_0 + ta1_y_yyy_zz_0[i] * fe_0 - ta1_y_yyy_zz_1[i] * fe_0 +
                              ta_yyy_yzz_1[i] + ta1_y_yyy_yzz_0[i] * pa_y[i] - ta1_y_yyy_yzz_1[i] * pc_y[i];

        ta1_y_yyyy_zzz_0[i] = 3.0 * ta1_y_yy_zzz_0[i] * fe_0 - 3.0 * ta1_y_yy_zzz_1[i] * fe_0 + ta_yyy_zzz_1[i] + ta1_y_yyy_zzz_0[i] * pa_y[i] -
                              ta1_y_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 260-270 components of targeted buffer : GF

    auto ta1_y_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 260);

    auto ta1_y_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 261);

    auto ta1_y_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 262);

    auto ta1_y_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 263);

    auto ta1_y_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 264);

    auto ta1_y_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 265);

    auto ta1_y_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 266);

    auto ta1_y_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 267);

    auto ta1_y_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 268);

    auto ta1_y_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 269);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_yyy_xx_0,   \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xxx_0,  \
                             ta1_y_yyy_xxx_1,  \
                             ta1_y_yyy_xxy_0,  \
                             ta1_y_yyy_xxy_1,  \
                             ta1_y_yyy_xxz_0,  \
                             ta1_y_yyy_xxz_1,  \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xyy_0,  \
                             ta1_y_yyy_xyy_1,  \
                             ta1_y_yyy_xyz_0,  \
                             ta1_y_yyy_xyz_1,  \
                             ta1_y_yyy_xz_0,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_xzz_0,  \
                             ta1_y_yyy_xzz_1,  \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yyy_0,  \
                             ta1_y_yyy_yyy_1,  \
                             ta1_y_yyy_yyz_0,  \
                             ta1_y_yyy_yyz_1,  \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_yzz_0,  \
                             ta1_y_yyy_yzz_1,  \
                             ta1_y_yyy_zz_0,   \
                             ta1_y_yyy_zz_1,   \
                             ta1_y_yyy_zzz_0,  \
                             ta1_y_yyy_zzz_1,  \
                             ta1_y_yyyz_xxx_0, \
                             ta1_y_yyyz_xxy_0, \
                             ta1_y_yyyz_xxz_0, \
                             ta1_y_yyyz_xyy_0, \
                             ta1_y_yyyz_xyz_0, \
                             ta1_y_yyyz_xzz_0, \
                             ta1_y_yyyz_yyy_0, \
                             ta1_y_yyyz_yyz_0, \
                             ta1_y_yyyz_yzz_0, \
                             ta1_y_yyyz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_xxx_0[i] = ta1_y_yyy_xxx_0[i] * pa_z[i] - ta1_y_yyy_xxx_1[i] * pc_z[i];

        ta1_y_yyyz_xxy_0[i] = ta1_y_yyy_xxy_0[i] * pa_z[i] - ta1_y_yyy_xxy_1[i] * pc_z[i];

        ta1_y_yyyz_xxz_0[i] = ta1_y_yyy_xx_0[i] * fe_0 - ta1_y_yyy_xx_1[i] * fe_0 + ta1_y_yyy_xxz_0[i] * pa_z[i] - ta1_y_yyy_xxz_1[i] * pc_z[i];

        ta1_y_yyyz_xyy_0[i] = ta1_y_yyy_xyy_0[i] * pa_z[i] - ta1_y_yyy_xyy_1[i] * pc_z[i];

        ta1_y_yyyz_xyz_0[i] = ta1_y_yyy_xy_0[i] * fe_0 - ta1_y_yyy_xy_1[i] * fe_0 + ta1_y_yyy_xyz_0[i] * pa_z[i] - ta1_y_yyy_xyz_1[i] * pc_z[i];

        ta1_y_yyyz_xzz_0[i] =
            2.0 * ta1_y_yyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xz_1[i] * fe_0 + ta1_y_yyy_xzz_0[i] * pa_z[i] - ta1_y_yyy_xzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyy_0[i] = ta1_y_yyy_yyy_0[i] * pa_z[i] - ta1_y_yyy_yyy_1[i] * pc_z[i];

        ta1_y_yyyz_yyz_0[i] = ta1_y_yyy_yy_0[i] * fe_0 - ta1_y_yyy_yy_1[i] * fe_0 + ta1_y_yyy_yyz_0[i] * pa_z[i] - ta1_y_yyy_yyz_1[i] * pc_z[i];

        ta1_y_yyyz_yzz_0[i] =
            2.0 * ta1_y_yyy_yz_0[i] * fe_0 - 2.0 * ta1_y_yyy_yz_1[i] * fe_0 + ta1_y_yyy_yzz_0[i] * pa_z[i] - ta1_y_yyy_yzz_1[i] * pc_z[i];

        ta1_y_yyyz_zzz_0[i] =
            3.0 * ta1_y_yyy_zz_0[i] * fe_0 - 3.0 * ta1_y_yyy_zz_1[i] * fe_0 + ta1_y_yyy_zzz_0[i] * pa_z[i] - ta1_y_yyy_zzz_1[i] * pc_z[i];
    }

    // Set up 270-280 components of targeted buffer : GF

    auto ta1_y_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 270);

    auto ta1_y_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 271);

    auto ta1_y_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 272);

    auto ta1_y_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 273);

    auto ta1_y_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 274);

    auto ta1_y_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 275);

    auto ta1_y_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 276);

    auto ta1_y_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 277);

    auto ta1_y_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 278);

    auto ta1_y_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 279);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yy_xxx_0,   \
                             ta1_y_yy_xxx_1,   \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yyz_xxx_0,  \
                             ta1_y_yyz_xxx_1,  \
                             ta1_y_yyz_xxy_0,  \
                             ta1_y_yyz_xxy_1,  \
                             ta1_y_yyz_xy_0,   \
                             ta1_y_yyz_xy_1,   \
                             ta1_y_yyz_xyy_0,  \
                             ta1_y_yyz_xyy_1,  \
                             ta1_y_yyz_xyz_0,  \
                             ta1_y_yyz_xyz_1,  \
                             ta1_y_yyz_yy_0,   \
                             ta1_y_yyz_yy_1,   \
                             ta1_y_yyz_yyy_0,  \
                             ta1_y_yyz_yyy_1,  \
                             ta1_y_yyz_yyz_0,  \
                             ta1_y_yyz_yyz_1,  \
                             ta1_y_yyz_yz_0,   \
                             ta1_y_yyz_yz_1,   \
                             ta1_y_yyz_yzz_0,  \
                             ta1_y_yyz_yzz_1,  \
                             ta1_y_yyzz_xxx_0, \
                             ta1_y_yyzz_xxy_0, \
                             ta1_y_yyzz_xxz_0, \
                             ta1_y_yyzz_xyy_0, \
                             ta1_y_yyzz_xyz_0, \
                             ta1_y_yyzz_xzz_0, \
                             ta1_y_yyzz_yyy_0, \
                             ta1_y_yyzz_yyz_0, \
                             ta1_y_yyzz_yzz_0, \
                             ta1_y_yyzz_zzz_0, \
                             ta1_y_yzz_xxz_0,  \
                             ta1_y_yzz_xxz_1,  \
                             ta1_y_yzz_xzz_0,  \
                             ta1_y_yzz_xzz_1,  \
                             ta1_y_yzz_zzz_0,  \
                             ta1_y_yzz_zzz_1,  \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             ta_yzz_xxz_1,     \
                             ta_yzz_xzz_1,     \
                             ta_yzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_xxx_0[i] = ta1_y_yy_xxx_0[i] * fe_0 - ta1_y_yy_xxx_1[i] * fe_0 + ta1_y_yyz_xxx_0[i] * pa_z[i] - ta1_y_yyz_xxx_1[i] * pc_z[i];

        ta1_y_yyzz_xxy_0[i] = ta1_y_yy_xxy_0[i] * fe_0 - ta1_y_yy_xxy_1[i] * fe_0 + ta1_y_yyz_xxy_0[i] * pa_z[i] - ta1_y_yyz_xxy_1[i] * pc_z[i];

        ta1_y_yyzz_xxz_0[i] =
            ta1_y_zz_xxz_0[i] * fe_0 - ta1_y_zz_xxz_1[i] * fe_0 + ta_yzz_xxz_1[i] + ta1_y_yzz_xxz_0[i] * pa_y[i] - ta1_y_yzz_xxz_1[i] * pc_y[i];

        ta1_y_yyzz_xyy_0[i] = ta1_y_yy_xyy_0[i] * fe_0 - ta1_y_yy_xyy_1[i] * fe_0 + ta1_y_yyz_xyy_0[i] * pa_z[i] - ta1_y_yyz_xyy_1[i] * pc_z[i];

        ta1_y_yyzz_xyz_0[i] = ta1_y_yy_xyz_0[i] * fe_0 - ta1_y_yy_xyz_1[i] * fe_0 + ta1_y_yyz_xy_0[i] * fe_0 - ta1_y_yyz_xy_1[i] * fe_0 +
                              ta1_y_yyz_xyz_0[i] * pa_z[i] - ta1_y_yyz_xyz_1[i] * pc_z[i];

        ta1_y_yyzz_xzz_0[i] =
            ta1_y_zz_xzz_0[i] * fe_0 - ta1_y_zz_xzz_1[i] * fe_0 + ta_yzz_xzz_1[i] + ta1_y_yzz_xzz_0[i] * pa_y[i] - ta1_y_yzz_xzz_1[i] * pc_y[i];

        ta1_y_yyzz_yyy_0[i] = ta1_y_yy_yyy_0[i] * fe_0 - ta1_y_yy_yyy_1[i] * fe_0 + ta1_y_yyz_yyy_0[i] * pa_z[i] - ta1_y_yyz_yyy_1[i] * pc_z[i];

        ta1_y_yyzz_yyz_0[i] = ta1_y_yy_yyz_0[i] * fe_0 - ta1_y_yy_yyz_1[i] * fe_0 + ta1_y_yyz_yy_0[i] * fe_0 - ta1_y_yyz_yy_1[i] * fe_0 +
                              ta1_y_yyz_yyz_0[i] * pa_z[i] - ta1_y_yyz_yyz_1[i] * pc_z[i];

        ta1_y_yyzz_yzz_0[i] = ta1_y_yy_yzz_0[i] * fe_0 - ta1_y_yy_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_yz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yz_1[i] * fe_0 +
                              ta1_y_yyz_yzz_0[i] * pa_z[i] - ta1_y_yyz_yzz_1[i] * pc_z[i];

        ta1_y_yyzz_zzz_0[i] =
            ta1_y_zz_zzz_0[i] * fe_0 - ta1_y_zz_zzz_1[i] * fe_0 + ta_yzz_zzz_1[i] + ta1_y_yzz_zzz_0[i] * pa_y[i] - ta1_y_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 280-290 components of targeted buffer : GF

    auto ta1_y_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 280);

    auto ta1_y_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 281);

    auto ta1_y_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 282);

    auto ta1_y_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 283);

    auto ta1_y_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 284);

    auto ta1_y_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 285);

    auto ta1_y_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 286);

    auto ta1_y_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 287);

    auto ta1_y_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 288);

    auto ta1_y_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 289);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yz_xxy_0,   \
                             ta1_y_yz_xxy_1,   \
                             ta1_y_yz_xyy_0,   \
                             ta1_y_yz_xyy_1,   \
                             ta1_y_yz_yyy_0,   \
                             ta1_y_yz_yyy_1,   \
                             ta1_y_yzz_xxy_0,  \
                             ta1_y_yzz_xxy_1,  \
                             ta1_y_yzz_xyy_0,  \
                             ta1_y_yzz_xyy_1,  \
                             ta1_y_yzz_yyy_0,  \
                             ta1_y_yzz_yyy_1,  \
                             ta1_y_yzzz_xxx_0, \
                             ta1_y_yzzz_xxy_0, \
                             ta1_y_yzzz_xxz_0, \
                             ta1_y_yzzz_xyy_0, \
                             ta1_y_yzzz_xyz_0, \
                             ta1_y_yzzz_xzz_0, \
                             ta1_y_yzzz_yyy_0, \
                             ta1_y_yzzz_yyz_0, \
                             ta1_y_yzzz_yzz_0, \
                             ta1_y_yzzz_zzz_0, \
                             ta1_y_zzz_xxx_0,  \
                             ta1_y_zzz_xxx_1,  \
                             ta1_y_zzz_xxz_0,  \
                             ta1_y_zzz_xxz_1,  \
                             ta1_y_zzz_xyz_0,  \
                             ta1_y_zzz_xyz_1,  \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_xzz_0,  \
                             ta1_y_zzz_xzz_1,  \
                             ta1_y_zzz_yyz_0,  \
                             ta1_y_zzz_yyz_1,  \
                             ta1_y_zzz_yz_0,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_yzz_0,  \
                             ta1_y_zzz_yzz_1,  \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             ta1_y_zzz_zzz_0,  \
                             ta1_y_zzz_zzz_1,  \
                             ta_zzz_xxx_1,     \
                             ta_zzz_xxz_1,     \
                             ta_zzz_xyz_1,     \
                             ta_zzz_xzz_1,     \
                             ta_zzz_yyz_1,     \
                             ta_zzz_yzz_1,     \
                             ta_zzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_xxx_0[i] = ta_zzz_xxx_1[i] + ta1_y_zzz_xxx_0[i] * pa_y[i] - ta1_y_zzz_xxx_1[i] * pc_y[i];

        ta1_y_yzzz_xxy_0[i] =
            2.0 * ta1_y_yz_xxy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxy_1[i] * fe_0 + ta1_y_yzz_xxy_0[i] * pa_z[i] - ta1_y_yzz_xxy_1[i] * pc_z[i];

        ta1_y_yzzz_xxz_0[i] = ta_zzz_xxz_1[i] + ta1_y_zzz_xxz_0[i] * pa_y[i] - ta1_y_zzz_xxz_1[i] * pc_y[i];

        ta1_y_yzzz_xyy_0[i] =
            2.0 * ta1_y_yz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xyy_1[i] * fe_0 + ta1_y_yzz_xyy_0[i] * pa_z[i] - ta1_y_yzz_xyy_1[i] * pc_z[i];

        ta1_y_yzzz_xyz_0[i] =
            ta1_y_zzz_xz_0[i] * fe_0 - ta1_y_zzz_xz_1[i] * fe_0 + ta_zzz_xyz_1[i] + ta1_y_zzz_xyz_0[i] * pa_y[i] - ta1_y_zzz_xyz_1[i] * pc_y[i];

        ta1_y_yzzz_xzz_0[i] = ta_zzz_xzz_1[i] + ta1_y_zzz_xzz_0[i] * pa_y[i] - ta1_y_zzz_xzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyy_0[i] =
            2.0 * ta1_y_yz_yyy_0[i] * fe_0 - 2.0 * ta1_y_yz_yyy_1[i] * fe_0 + ta1_y_yzz_yyy_0[i] * pa_z[i] - ta1_y_yzz_yyy_1[i] * pc_z[i];

        ta1_y_yzzz_yyz_0[i] = 2.0 * ta1_y_zzz_yz_0[i] * fe_0 - 2.0 * ta1_y_zzz_yz_1[i] * fe_0 + ta_zzz_yyz_1[i] + ta1_y_zzz_yyz_0[i] * pa_y[i] -
                              ta1_y_zzz_yyz_1[i] * pc_y[i];

        ta1_y_yzzz_yzz_0[i] =
            ta1_y_zzz_zz_0[i] * fe_0 - ta1_y_zzz_zz_1[i] * fe_0 + ta_zzz_yzz_1[i] + ta1_y_zzz_yzz_0[i] * pa_y[i] - ta1_y_zzz_yzz_1[i] * pc_y[i];

        ta1_y_yzzz_zzz_0[i] = ta_zzz_zzz_1[i] + ta1_y_zzz_zzz_0[i] * pa_y[i] - ta1_y_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 290-300 components of targeted buffer : GF

    auto ta1_y_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 290);

    auto ta1_y_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 291);

    auto ta1_y_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 292);

    auto ta1_y_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 293);

    auto ta1_y_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 294);

    auto ta1_y_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 295);

    auto ta1_y_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 296);

    auto ta1_y_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 297);

    auto ta1_y_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 298);

    auto ta1_y_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 299);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_zz_xxx_0,   \
                             ta1_y_zz_xxx_1,   \
                             ta1_y_zz_xxy_0,   \
                             ta1_y_zz_xxy_1,   \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xyy_0,   \
                             ta1_y_zz_xyy_1,   \
                             ta1_y_zz_xyz_0,   \
                             ta1_y_zz_xyz_1,   \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_yyy_0,   \
                             ta1_y_zz_yyy_1,   \
                             ta1_y_zz_yyz_0,   \
                             ta1_y_zz_yyz_1,   \
                             ta1_y_zz_yzz_0,   \
                             ta1_y_zz_yzz_1,   \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             ta1_y_zzz_xx_0,   \
                             ta1_y_zzz_xx_1,   \
                             ta1_y_zzz_xxx_0,  \
                             ta1_y_zzz_xxx_1,  \
                             ta1_y_zzz_xxy_0,  \
                             ta1_y_zzz_xxy_1,  \
                             ta1_y_zzz_xxz_0,  \
                             ta1_y_zzz_xxz_1,  \
                             ta1_y_zzz_xy_0,   \
                             ta1_y_zzz_xy_1,   \
                             ta1_y_zzz_xyy_0,  \
                             ta1_y_zzz_xyy_1,  \
                             ta1_y_zzz_xyz_0,  \
                             ta1_y_zzz_xyz_1,  \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_xzz_0,  \
                             ta1_y_zzz_xzz_1,  \
                             ta1_y_zzz_yy_0,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yyy_0,  \
                             ta1_y_zzz_yyy_1,  \
                             ta1_y_zzz_yyz_0,  \
                             ta1_y_zzz_yyz_1,  \
                             ta1_y_zzz_yz_0,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_yzz_0,  \
                             ta1_y_zzz_yzz_1,  \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             ta1_y_zzz_zzz_0,  \
                             ta1_y_zzz_zzz_1,  \
                             ta1_y_zzzz_xxx_0, \
                             ta1_y_zzzz_xxy_0, \
                             ta1_y_zzzz_xxz_0, \
                             ta1_y_zzzz_xyy_0, \
                             ta1_y_zzzz_xyz_0, \
                             ta1_y_zzzz_xzz_0, \
                             ta1_y_zzzz_yyy_0, \
                             ta1_y_zzzz_yyz_0, \
                             ta1_y_zzzz_yzz_0, \
                             ta1_y_zzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_xxx_0[i] =
            3.0 * ta1_y_zz_xxx_0[i] * fe_0 - 3.0 * ta1_y_zz_xxx_1[i] * fe_0 + ta1_y_zzz_xxx_0[i] * pa_z[i] - ta1_y_zzz_xxx_1[i] * pc_z[i];

        ta1_y_zzzz_xxy_0[i] =
            3.0 * ta1_y_zz_xxy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxy_1[i] * fe_0 + ta1_y_zzz_xxy_0[i] * pa_z[i] - ta1_y_zzz_xxy_1[i] * pc_z[i];

        ta1_y_zzzz_xxz_0[i] = 3.0 * ta1_y_zz_xxz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxz_1[i] * fe_0 + ta1_y_zzz_xx_0[i] * fe_0 - ta1_y_zzz_xx_1[i] * fe_0 +
                              ta1_y_zzz_xxz_0[i] * pa_z[i] - ta1_y_zzz_xxz_1[i] * pc_z[i];

        ta1_y_zzzz_xyy_0[i] =
            3.0 * ta1_y_zz_xyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xyy_1[i] * fe_0 + ta1_y_zzz_xyy_0[i] * pa_z[i] - ta1_y_zzz_xyy_1[i] * pc_z[i];

        ta1_y_zzzz_xyz_0[i] = 3.0 * ta1_y_zz_xyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyz_1[i] * fe_0 + ta1_y_zzz_xy_0[i] * fe_0 - ta1_y_zzz_xy_1[i] * fe_0 +
                              ta1_y_zzz_xyz_0[i] * pa_z[i] - ta1_y_zzz_xyz_1[i] * pc_z[i];

        ta1_y_zzzz_xzz_0[i] = 3.0 * ta1_y_zz_xzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xz_0[i] * fe_0 -
                              2.0 * ta1_y_zzz_xz_1[i] * fe_0 + ta1_y_zzz_xzz_0[i] * pa_z[i] - ta1_y_zzz_xzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyy_0[i] =
            3.0 * ta1_y_zz_yyy_0[i] * fe_0 - 3.0 * ta1_y_zz_yyy_1[i] * fe_0 + ta1_y_zzz_yyy_0[i] * pa_z[i] - ta1_y_zzz_yyy_1[i] * pc_z[i];

        ta1_y_zzzz_yyz_0[i] = 3.0 * ta1_y_zz_yyz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyz_1[i] * fe_0 + ta1_y_zzz_yy_0[i] * fe_0 - ta1_y_zzz_yy_1[i] * fe_0 +
                              ta1_y_zzz_yyz_0[i] * pa_z[i] - ta1_y_zzz_yyz_1[i] * pc_z[i];

        ta1_y_zzzz_yzz_0[i] = 3.0 * ta1_y_zz_yzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yz_0[i] * fe_0 -
                              2.0 * ta1_y_zzz_yz_1[i] * fe_0 + ta1_y_zzz_yzz_0[i] * pa_z[i] - ta1_y_zzz_yzz_1[i] * pc_z[i];

        ta1_y_zzzz_zzz_0[i] = 3.0 * ta1_y_zz_zzz_0[i] * fe_0 - 3.0 * ta1_y_zz_zzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_zz_0[i] * fe_0 -
                              3.0 * ta1_y_zzz_zz_1[i] * fe_0 + ta1_y_zzz_zzz_0[i] * pa_z[i] - ta1_y_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 300-310 components of targeted buffer : GF

    auto ta1_z_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 300);

    auto ta1_z_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 301);

    auto ta1_z_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 302);

    auto ta1_z_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 303);

    auto ta1_z_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 304);

    auto ta1_z_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 305);

    auto ta1_z_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 306);

    auto ta1_z_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 307);

    auto ta1_z_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 308);

    auto ta1_z_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 309);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxy_0,   \
                             ta1_z_xx_xxy_1,   \
                             ta1_z_xx_xxz_0,   \
                             ta1_z_xx_xxz_1,   \
                             ta1_z_xx_xyy_0,   \
                             ta1_z_xx_xyy_1,   \
                             ta1_z_xx_xyz_0,   \
                             ta1_z_xx_xyz_1,   \
                             ta1_z_xx_xzz_0,   \
                             ta1_z_xx_xzz_1,   \
                             ta1_z_xx_yyy_0,   \
                             ta1_z_xx_yyy_1,   \
                             ta1_z_xx_yyz_0,   \
                             ta1_z_xx_yyz_1,   \
                             ta1_z_xx_yzz_0,   \
                             ta1_z_xx_yzz_1,   \
                             ta1_z_xx_zzz_0,   \
                             ta1_z_xx_zzz_1,   \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xxx_0,  \
                             ta1_z_xxx_xxx_1,  \
                             ta1_z_xxx_xxy_0,  \
                             ta1_z_xxx_xxy_1,  \
                             ta1_z_xxx_xxz_0,  \
                             ta1_z_xxx_xxz_1,  \
                             ta1_z_xxx_xy_0,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xyy_0,  \
                             ta1_z_xxx_xyy_1,  \
                             ta1_z_xxx_xyz_0,  \
                             ta1_z_xxx_xyz_1,  \
                             ta1_z_xxx_xz_0,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_xzz_0,  \
                             ta1_z_xxx_xzz_1,  \
                             ta1_z_xxx_yy_0,   \
                             ta1_z_xxx_yy_1,   \
                             ta1_z_xxx_yyy_0,  \
                             ta1_z_xxx_yyy_1,  \
                             ta1_z_xxx_yyz_0,  \
                             ta1_z_xxx_yyz_1,  \
                             ta1_z_xxx_yz_0,   \
                             ta1_z_xxx_yz_1,   \
                             ta1_z_xxx_yzz_0,  \
                             ta1_z_xxx_yzz_1,  \
                             ta1_z_xxx_zz_0,   \
                             ta1_z_xxx_zz_1,   \
                             ta1_z_xxx_zzz_0,  \
                             ta1_z_xxx_zzz_1,  \
                             ta1_z_xxxx_xxx_0, \
                             ta1_z_xxxx_xxy_0, \
                             ta1_z_xxxx_xxz_0, \
                             ta1_z_xxxx_xyy_0, \
                             ta1_z_xxxx_xyz_0, \
                             ta1_z_xxxx_xzz_0, \
                             ta1_z_xxxx_yyy_0, \
                             ta1_z_xxxx_yyz_0, \
                             ta1_z_xxxx_yzz_0, \
                             ta1_z_xxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_xxx_0[i] = 3.0 * ta1_z_xx_xxx_0[i] * fe_0 - 3.0 * ta1_z_xx_xxx_1[i] * fe_0 + 3.0 * ta1_z_xxx_xx_0[i] * fe_0 -
                              3.0 * ta1_z_xxx_xx_1[i] * fe_0 + ta1_z_xxx_xxx_0[i] * pa_x[i] - ta1_z_xxx_xxx_1[i] * pc_x[i];

        ta1_z_xxxx_xxy_0[i] = 3.0 * ta1_z_xx_xxy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xy_0[i] * fe_0 -
                              2.0 * ta1_z_xxx_xy_1[i] * fe_0 + ta1_z_xxx_xxy_0[i] * pa_x[i] - ta1_z_xxx_xxy_1[i] * pc_x[i];

        ta1_z_xxxx_xxz_0[i] = 3.0 * ta1_z_xx_xxz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xz_0[i] * fe_0 -
                              2.0 * ta1_z_xxx_xz_1[i] * fe_0 + ta1_z_xxx_xxz_0[i] * pa_x[i] - ta1_z_xxx_xxz_1[i] * pc_x[i];

        ta1_z_xxxx_xyy_0[i] = 3.0 * ta1_z_xx_xyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xyy_1[i] * fe_0 + ta1_z_xxx_yy_0[i] * fe_0 - ta1_z_xxx_yy_1[i] * fe_0 +
                              ta1_z_xxx_xyy_0[i] * pa_x[i] - ta1_z_xxx_xyy_1[i] * pc_x[i];

        ta1_z_xxxx_xyz_0[i] = 3.0 * ta1_z_xx_xyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyz_1[i] * fe_0 + ta1_z_xxx_yz_0[i] * fe_0 - ta1_z_xxx_yz_1[i] * fe_0 +
                              ta1_z_xxx_xyz_0[i] * pa_x[i] - ta1_z_xxx_xyz_1[i] * pc_x[i];

        ta1_z_xxxx_xzz_0[i] = 3.0 * ta1_z_xx_xzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xzz_1[i] * fe_0 + ta1_z_xxx_zz_0[i] * fe_0 - ta1_z_xxx_zz_1[i] * fe_0 +
                              ta1_z_xxx_xzz_0[i] * pa_x[i] - ta1_z_xxx_xzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyy_0[i] =
            3.0 * ta1_z_xx_yyy_0[i] * fe_0 - 3.0 * ta1_z_xx_yyy_1[i] * fe_0 + ta1_z_xxx_yyy_0[i] * pa_x[i] - ta1_z_xxx_yyy_1[i] * pc_x[i];

        ta1_z_xxxx_yyz_0[i] =
            3.0 * ta1_z_xx_yyz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyz_1[i] * fe_0 + ta1_z_xxx_yyz_0[i] * pa_x[i] - ta1_z_xxx_yyz_1[i] * pc_x[i];

        ta1_z_xxxx_yzz_0[i] =
            3.0 * ta1_z_xx_yzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yzz_1[i] * fe_0 + ta1_z_xxx_yzz_0[i] * pa_x[i] - ta1_z_xxx_yzz_1[i] * pc_x[i];

        ta1_z_xxxx_zzz_0[i] =
            3.0 * ta1_z_xx_zzz_0[i] * fe_0 - 3.0 * ta1_z_xx_zzz_1[i] * fe_0 + ta1_z_xxx_zzz_0[i] * pa_x[i] - ta1_z_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : GF

    auto ta1_z_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 310);

    auto ta1_z_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 311);

    auto ta1_z_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 312);

    auto ta1_z_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 313);

    auto ta1_z_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 314);

    auto ta1_z_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 315);

    auto ta1_z_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 316);

    auto ta1_z_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 317);

    auto ta1_z_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 318);

    auto ta1_z_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 319);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xxx_0,  \
                             ta1_z_xxx_xxx_1,  \
                             ta1_z_xxx_xxy_0,  \
                             ta1_z_xxx_xxy_1,  \
                             ta1_z_xxx_xxz_0,  \
                             ta1_z_xxx_xxz_1,  \
                             ta1_z_xxx_xy_0,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xyy_0,  \
                             ta1_z_xxx_xyy_1,  \
                             ta1_z_xxx_xyz_0,  \
                             ta1_z_xxx_xyz_1,  \
                             ta1_z_xxx_xz_0,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_xzz_0,  \
                             ta1_z_xxx_xzz_1,  \
                             ta1_z_xxx_zzz_0,  \
                             ta1_z_xxx_zzz_1,  \
                             ta1_z_xxxy_xxx_0, \
                             ta1_z_xxxy_xxy_0, \
                             ta1_z_xxxy_xxz_0, \
                             ta1_z_xxxy_xyy_0, \
                             ta1_z_xxxy_xyz_0, \
                             ta1_z_xxxy_xzz_0, \
                             ta1_z_xxxy_yyy_0, \
                             ta1_z_xxxy_yyz_0, \
                             ta1_z_xxxy_yzz_0, \
                             ta1_z_xxxy_zzz_0, \
                             ta1_z_xxy_yyy_0,  \
                             ta1_z_xxy_yyy_1,  \
                             ta1_z_xxy_yyz_0,  \
                             ta1_z_xxy_yyz_1,  \
                             ta1_z_xxy_yzz_0,  \
                             ta1_z_xxy_yzz_1,  \
                             ta1_z_xy_yyy_0,   \
                             ta1_z_xy_yyy_1,   \
                             ta1_z_xy_yyz_0,   \
                             ta1_z_xy_yyz_1,   \
                             ta1_z_xy_yzz_0,   \
                             ta1_z_xy_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_xxx_0[i] = ta1_z_xxx_xxx_0[i] * pa_y[i] - ta1_z_xxx_xxx_1[i] * pc_y[i];

        ta1_z_xxxy_xxy_0[i] = ta1_z_xxx_xx_0[i] * fe_0 - ta1_z_xxx_xx_1[i] * fe_0 + ta1_z_xxx_xxy_0[i] * pa_y[i] - ta1_z_xxx_xxy_1[i] * pc_y[i];

        ta1_z_xxxy_xxz_0[i] = ta1_z_xxx_xxz_0[i] * pa_y[i] - ta1_z_xxx_xxz_1[i] * pc_y[i];

        ta1_z_xxxy_xyy_0[i] =
            2.0 * ta1_z_xxx_xy_0[i] * fe_0 - 2.0 * ta1_z_xxx_xy_1[i] * fe_0 + ta1_z_xxx_xyy_0[i] * pa_y[i] - ta1_z_xxx_xyy_1[i] * pc_y[i];

        ta1_z_xxxy_xyz_0[i] = ta1_z_xxx_xz_0[i] * fe_0 - ta1_z_xxx_xz_1[i] * fe_0 + ta1_z_xxx_xyz_0[i] * pa_y[i] - ta1_z_xxx_xyz_1[i] * pc_y[i];

        ta1_z_xxxy_xzz_0[i] = ta1_z_xxx_xzz_0[i] * pa_y[i] - ta1_z_xxx_xzz_1[i] * pc_y[i];

        ta1_z_xxxy_yyy_0[i] =
            2.0 * ta1_z_xy_yyy_0[i] * fe_0 - 2.0 * ta1_z_xy_yyy_1[i] * fe_0 + ta1_z_xxy_yyy_0[i] * pa_x[i] - ta1_z_xxy_yyy_1[i] * pc_x[i];

        ta1_z_xxxy_yyz_0[i] =
            2.0 * ta1_z_xy_yyz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyz_1[i] * fe_0 + ta1_z_xxy_yyz_0[i] * pa_x[i] - ta1_z_xxy_yyz_1[i] * pc_x[i];

        ta1_z_xxxy_yzz_0[i] =
            2.0 * ta1_z_xy_yzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yzz_1[i] * fe_0 + ta1_z_xxy_yzz_0[i] * pa_x[i] - ta1_z_xxy_yzz_1[i] * pc_x[i];

        ta1_z_xxxy_zzz_0[i] = ta1_z_xxx_zzz_0[i] * pa_y[i] - ta1_z_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 320-330 components of targeted buffer : GF

    auto ta1_z_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 320);

    auto ta1_z_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 321);

    auto ta1_z_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 322);

    auto ta1_z_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 323);

    auto ta1_z_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 324);

    auto ta1_z_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 325);

    auto ta1_z_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 326);

    auto ta1_z_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 327);

    auto ta1_z_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 328);

    auto ta1_z_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 329);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xxx_0,  \
                             ta1_z_xxx_xxx_1,  \
                             ta1_z_xxx_xxy_0,  \
                             ta1_z_xxx_xxy_1,  \
                             ta1_z_xxx_xxz_0,  \
                             ta1_z_xxx_xxz_1,  \
                             ta1_z_xxx_xy_0,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xyy_0,  \
                             ta1_z_xxx_xyy_1,  \
                             ta1_z_xxx_xyz_0,  \
                             ta1_z_xxx_xyz_1,  \
                             ta1_z_xxx_xz_0,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_xzz_0,  \
                             ta1_z_xxx_xzz_1,  \
                             ta1_z_xxx_yyy_0,  \
                             ta1_z_xxx_yyy_1,  \
                             ta1_z_xxxz_xxx_0, \
                             ta1_z_xxxz_xxy_0, \
                             ta1_z_xxxz_xxz_0, \
                             ta1_z_xxxz_xyy_0, \
                             ta1_z_xxxz_xyz_0, \
                             ta1_z_xxxz_xzz_0, \
                             ta1_z_xxxz_yyy_0, \
                             ta1_z_xxxz_yyz_0, \
                             ta1_z_xxxz_yzz_0, \
                             ta1_z_xxxz_zzz_0, \
                             ta1_z_xxz_yyz_0,  \
                             ta1_z_xxz_yyz_1,  \
                             ta1_z_xxz_yzz_0,  \
                             ta1_z_xxz_yzz_1,  \
                             ta1_z_xxz_zzz_0,  \
                             ta1_z_xxz_zzz_1,  \
                             ta1_z_xz_yyz_0,   \
                             ta1_z_xz_yyz_1,   \
                             ta1_z_xz_yzz_0,   \
                             ta1_z_xz_yzz_1,   \
                             ta1_z_xz_zzz_0,   \
                             ta1_z_xz_zzz_1,   \
                             ta_xxx_xxx_1,     \
                             ta_xxx_xxy_1,     \
                             ta_xxx_xxz_1,     \
                             ta_xxx_xyy_1,     \
                             ta_xxx_xyz_1,     \
                             ta_xxx_xzz_1,     \
                             ta_xxx_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_xxx_0[i] = ta_xxx_xxx_1[i] + ta1_z_xxx_xxx_0[i] * pa_z[i] - ta1_z_xxx_xxx_1[i] * pc_z[i];

        ta1_z_xxxz_xxy_0[i] = ta_xxx_xxy_1[i] + ta1_z_xxx_xxy_0[i] * pa_z[i] - ta1_z_xxx_xxy_1[i] * pc_z[i];

        ta1_z_xxxz_xxz_0[i] =
            ta1_z_xxx_xx_0[i] * fe_0 - ta1_z_xxx_xx_1[i] * fe_0 + ta_xxx_xxz_1[i] + ta1_z_xxx_xxz_0[i] * pa_z[i] - ta1_z_xxx_xxz_1[i] * pc_z[i];

        ta1_z_xxxz_xyy_0[i] = ta_xxx_xyy_1[i] + ta1_z_xxx_xyy_0[i] * pa_z[i] - ta1_z_xxx_xyy_1[i] * pc_z[i];

        ta1_z_xxxz_xyz_0[i] =
            ta1_z_xxx_xy_0[i] * fe_0 - ta1_z_xxx_xy_1[i] * fe_0 + ta_xxx_xyz_1[i] + ta1_z_xxx_xyz_0[i] * pa_z[i] - ta1_z_xxx_xyz_1[i] * pc_z[i];

        ta1_z_xxxz_xzz_0[i] = 2.0 * ta1_z_xxx_xz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xz_1[i] * fe_0 + ta_xxx_xzz_1[i] + ta1_z_xxx_xzz_0[i] * pa_z[i] -
                              ta1_z_xxx_xzz_1[i] * pc_z[i];

        ta1_z_xxxz_yyy_0[i] = ta_xxx_yyy_1[i] + ta1_z_xxx_yyy_0[i] * pa_z[i] - ta1_z_xxx_yyy_1[i] * pc_z[i];

        ta1_z_xxxz_yyz_0[i] =
            2.0 * ta1_z_xz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyz_1[i] * fe_0 + ta1_z_xxz_yyz_0[i] * pa_x[i] - ta1_z_xxz_yyz_1[i] * pc_x[i];

        ta1_z_xxxz_yzz_0[i] =
            2.0 * ta1_z_xz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yzz_1[i] * fe_0 + ta1_z_xxz_yzz_0[i] * pa_x[i] - ta1_z_xxz_yzz_1[i] * pc_x[i];

        ta1_z_xxxz_zzz_0[i] =
            2.0 * ta1_z_xz_zzz_0[i] * fe_0 - 2.0 * ta1_z_xz_zzz_1[i] * fe_0 + ta1_z_xxz_zzz_0[i] * pa_x[i] - ta1_z_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 330-340 components of targeted buffer : GF

    auto ta1_z_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 330);

    auto ta1_z_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 331);

    auto ta1_z_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 332);

    auto ta1_z_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 333);

    auto ta1_z_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 334);

    auto ta1_z_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 335);

    auto ta1_z_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 336);

    auto ta1_z_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 337);

    auto ta1_z_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 338);

    auto ta1_z_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 339);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxz_0,   \
                             ta1_z_xx_xxz_1,   \
                             ta1_z_xx_xzz_0,   \
                             ta1_z_xx_xzz_1,   \
                             ta1_z_xxy_xxx_0,  \
                             ta1_z_xxy_xxx_1,  \
                             ta1_z_xxy_xxz_0,  \
                             ta1_z_xxy_xxz_1,  \
                             ta1_z_xxy_xzz_0,  \
                             ta1_z_xxy_xzz_1,  \
                             ta1_z_xxyy_xxx_0, \
                             ta1_z_xxyy_xxy_0, \
                             ta1_z_xxyy_xxz_0, \
                             ta1_z_xxyy_xyy_0, \
                             ta1_z_xxyy_xyz_0, \
                             ta1_z_xxyy_xzz_0, \
                             ta1_z_xxyy_yyy_0, \
                             ta1_z_xxyy_yyz_0, \
                             ta1_z_xxyy_yzz_0, \
                             ta1_z_xxyy_zzz_0, \
                             ta1_z_xyy_xxy_0,  \
                             ta1_z_xyy_xxy_1,  \
                             ta1_z_xyy_xy_0,   \
                             ta1_z_xyy_xy_1,   \
                             ta1_z_xyy_xyy_0,  \
                             ta1_z_xyy_xyy_1,  \
                             ta1_z_xyy_xyz_0,  \
                             ta1_z_xyy_xyz_1,  \
                             ta1_z_xyy_yy_0,   \
                             ta1_z_xyy_yy_1,   \
                             ta1_z_xyy_yyy_0,  \
                             ta1_z_xyy_yyy_1,  \
                             ta1_z_xyy_yyz_0,  \
                             ta1_z_xyy_yyz_1,  \
                             ta1_z_xyy_yz_0,   \
                             ta1_z_xyy_yz_1,   \
                             ta1_z_xyy_yzz_0,  \
                             ta1_z_xyy_yzz_1,  \
                             ta1_z_xyy_zzz_0,  \
                             ta1_z_xyy_zzz_1,  \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_xyz_0,   \
                             ta1_z_yy_xyz_1,   \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yy_yyz_0,   \
                             ta1_z_yy_yyz_1,   \
                             ta1_z_yy_yzz_0,   \
                             ta1_z_yy_yzz_1,   \
                             ta1_z_yy_zzz_0,   \
                             ta1_z_yy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_xxx_0[i] = ta1_z_xx_xxx_0[i] * fe_0 - ta1_z_xx_xxx_1[i] * fe_0 + ta1_z_xxy_xxx_0[i] * pa_y[i] - ta1_z_xxy_xxx_1[i] * pc_y[i];

        ta1_z_xxyy_xxy_0[i] = ta1_z_yy_xxy_0[i] * fe_0 - ta1_z_yy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xyy_xy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xy_1[i] * fe_0 +
                              ta1_z_xyy_xxy_0[i] * pa_x[i] - ta1_z_xyy_xxy_1[i] * pc_x[i];

        ta1_z_xxyy_xxz_0[i] = ta1_z_xx_xxz_0[i] * fe_0 - ta1_z_xx_xxz_1[i] * fe_0 + ta1_z_xxy_xxz_0[i] * pa_y[i] - ta1_z_xxy_xxz_1[i] * pc_y[i];

        ta1_z_xxyy_xyy_0[i] = ta1_z_yy_xyy_0[i] * fe_0 - ta1_z_yy_xyy_1[i] * fe_0 + ta1_z_xyy_yy_0[i] * fe_0 - ta1_z_xyy_yy_1[i] * fe_0 +
                              ta1_z_xyy_xyy_0[i] * pa_x[i] - ta1_z_xyy_xyy_1[i] * pc_x[i];

        ta1_z_xxyy_xyz_0[i] = ta1_z_yy_xyz_0[i] * fe_0 - ta1_z_yy_xyz_1[i] * fe_0 + ta1_z_xyy_yz_0[i] * fe_0 - ta1_z_xyy_yz_1[i] * fe_0 +
                              ta1_z_xyy_xyz_0[i] * pa_x[i] - ta1_z_xyy_xyz_1[i] * pc_x[i];

        ta1_z_xxyy_xzz_0[i] = ta1_z_xx_xzz_0[i] * fe_0 - ta1_z_xx_xzz_1[i] * fe_0 + ta1_z_xxy_xzz_0[i] * pa_y[i] - ta1_z_xxy_xzz_1[i] * pc_y[i];

        ta1_z_xxyy_yyy_0[i] = ta1_z_yy_yyy_0[i] * fe_0 - ta1_z_yy_yyy_1[i] * fe_0 + ta1_z_xyy_yyy_0[i] * pa_x[i] - ta1_z_xyy_yyy_1[i] * pc_x[i];

        ta1_z_xxyy_yyz_0[i] = ta1_z_yy_yyz_0[i] * fe_0 - ta1_z_yy_yyz_1[i] * fe_0 + ta1_z_xyy_yyz_0[i] * pa_x[i] - ta1_z_xyy_yyz_1[i] * pc_x[i];

        ta1_z_xxyy_yzz_0[i] = ta1_z_yy_yzz_0[i] * fe_0 - ta1_z_yy_yzz_1[i] * fe_0 + ta1_z_xyy_yzz_0[i] * pa_x[i] - ta1_z_xyy_yzz_1[i] * pc_x[i];

        ta1_z_xxyy_zzz_0[i] = ta1_z_yy_zzz_0[i] * fe_0 - ta1_z_yy_zzz_1[i] * fe_0 + ta1_z_xyy_zzz_0[i] * pa_x[i] - ta1_z_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 340-350 components of targeted buffer : GF

    auto ta1_z_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 340);

    auto ta1_z_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 341);

    auto ta1_z_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 342);

    auto ta1_z_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 343);

    auto ta1_z_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 344);

    auto ta1_z_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 345);

    auto ta1_z_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 346);

    auto ta1_z_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 347);

    auto ta1_z_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 348);

    auto ta1_z_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 349);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_xxy_xxy_0,  \
                             ta1_z_xxy_xxy_1,  \
                             ta1_z_xxy_xyy_0,  \
                             ta1_z_xxy_xyy_1,  \
                             ta1_z_xxy_yyy_0,  \
                             ta1_z_xxy_yyy_1,  \
                             ta1_z_xxyz_xxx_0, \
                             ta1_z_xxyz_xxy_0, \
                             ta1_z_xxyz_xxz_0, \
                             ta1_z_xxyz_xyy_0, \
                             ta1_z_xxyz_xyz_0, \
                             ta1_z_xxyz_xzz_0, \
                             ta1_z_xxyz_yyy_0, \
                             ta1_z_xxyz_yyz_0, \
                             ta1_z_xxyz_yzz_0, \
                             ta1_z_xxyz_zzz_0, \
                             ta1_z_xxz_xxx_0,  \
                             ta1_z_xxz_xxx_1,  \
                             ta1_z_xxz_xxz_0,  \
                             ta1_z_xxz_xxz_1,  \
                             ta1_z_xxz_xyz_0,  \
                             ta1_z_xxz_xyz_1,  \
                             ta1_z_xxz_xz_0,   \
                             ta1_z_xxz_xz_1,   \
                             ta1_z_xxz_xzz_0,  \
                             ta1_z_xxz_xzz_1,  \
                             ta1_z_xxz_zzz_0,  \
                             ta1_z_xxz_zzz_1,  \
                             ta1_z_xyz_yyz_0,  \
                             ta1_z_xyz_yyz_1,  \
                             ta1_z_xyz_yzz_0,  \
                             ta1_z_xyz_yzz_1,  \
                             ta1_z_yz_yyz_0,   \
                             ta1_z_yz_yyz_1,   \
                             ta1_z_yz_yzz_0,   \
                             ta1_z_yz_yzz_1,   \
                             ta_xxy_xxy_1,     \
                             ta_xxy_xyy_1,     \
                             ta_xxy_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyz_xxx_0[i] = ta1_z_xxz_xxx_0[i] * pa_y[i] - ta1_z_xxz_xxx_1[i] * pc_y[i];

        ta1_z_xxyz_xxy_0[i] = ta_xxy_xxy_1[i] + ta1_z_xxy_xxy_0[i] * pa_z[i] - ta1_z_xxy_xxy_1[i] * pc_z[i];

        ta1_z_xxyz_xxz_0[i] = ta1_z_xxz_xxz_0[i] * pa_y[i] - ta1_z_xxz_xxz_1[i] * pc_y[i];

        ta1_z_xxyz_xyy_0[i] = ta_xxy_xyy_1[i] + ta1_z_xxy_xyy_0[i] * pa_z[i] - ta1_z_xxy_xyy_1[i] * pc_z[i];

        ta1_z_xxyz_xyz_0[i] = ta1_z_xxz_xz_0[i] * fe_0 - ta1_z_xxz_xz_1[i] * fe_0 + ta1_z_xxz_xyz_0[i] * pa_y[i] - ta1_z_xxz_xyz_1[i] * pc_y[i];

        ta1_z_xxyz_xzz_0[i] = ta1_z_xxz_xzz_0[i] * pa_y[i] - ta1_z_xxz_xzz_1[i] * pc_y[i];

        ta1_z_xxyz_yyy_0[i] = ta_xxy_yyy_1[i] + ta1_z_xxy_yyy_0[i] * pa_z[i] - ta1_z_xxy_yyy_1[i] * pc_z[i];

        ta1_z_xxyz_yyz_0[i] = ta1_z_yz_yyz_0[i] * fe_0 - ta1_z_yz_yyz_1[i] * fe_0 + ta1_z_xyz_yyz_0[i] * pa_x[i] - ta1_z_xyz_yyz_1[i] * pc_x[i];

        ta1_z_xxyz_yzz_0[i] = ta1_z_yz_yzz_0[i] * fe_0 - ta1_z_yz_yzz_1[i] * fe_0 + ta1_z_xyz_yzz_0[i] * pa_x[i] - ta1_z_xyz_yzz_1[i] * pc_x[i];

        ta1_z_xxyz_zzz_0[i] = ta1_z_xxz_zzz_0[i] * pa_y[i] - ta1_z_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 350-360 components of targeted buffer : GF

    auto ta1_z_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 350);

    auto ta1_z_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 351);

    auto ta1_z_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 352);

    auto ta1_z_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 353);

    auto ta1_z_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 354);

    auto ta1_z_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 355);

    auto ta1_z_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 356);

    auto ta1_z_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 357);

    auto ta1_z_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 358);

    auto ta1_z_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 359);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxy_0,   \
                             ta1_z_xx_xxy_1,   \
                             ta1_z_xx_xyy_0,   \
                             ta1_z_xx_xyy_1,   \
                             ta1_z_xxz_xxx_0,  \
                             ta1_z_xxz_xxx_1,  \
                             ta1_z_xxz_xxy_0,  \
                             ta1_z_xxz_xxy_1,  \
                             ta1_z_xxz_xyy_0,  \
                             ta1_z_xxz_xyy_1,  \
                             ta1_z_xxzz_xxx_0, \
                             ta1_z_xxzz_xxy_0, \
                             ta1_z_xxzz_xxz_0, \
                             ta1_z_xxzz_xyy_0, \
                             ta1_z_xxzz_xyz_0, \
                             ta1_z_xxzz_xzz_0, \
                             ta1_z_xxzz_yyy_0, \
                             ta1_z_xxzz_yyz_0, \
                             ta1_z_xxzz_yzz_0, \
                             ta1_z_xxzz_zzz_0, \
                             ta1_z_xzz_xxz_0,  \
                             ta1_z_xzz_xxz_1,  \
                             ta1_z_xzz_xyz_0,  \
                             ta1_z_xzz_xyz_1,  \
                             ta1_z_xzz_xz_0,   \
                             ta1_z_xzz_xz_1,   \
                             ta1_z_xzz_xzz_0,  \
                             ta1_z_xzz_xzz_1,  \
                             ta1_z_xzz_yyy_0,  \
                             ta1_z_xzz_yyy_1,  \
                             ta1_z_xzz_yyz_0,  \
                             ta1_z_xzz_yyz_1,  \
                             ta1_z_xzz_yz_0,   \
                             ta1_z_xzz_yz_1,   \
                             ta1_z_xzz_yzz_0,  \
                             ta1_z_xzz_yzz_1,  \
                             ta1_z_xzz_zz_0,   \
                             ta1_z_xzz_zz_1,   \
                             ta1_z_xzz_zzz_0,  \
                             ta1_z_xzz_zzz_1,  \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_yyy_0,   \
                             ta1_z_zz_yyy_1,   \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta_xxz_xxx_1,     \
                             ta_xxz_xxy_1,     \
                             ta_xxz_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_xxx_0[i] =
            ta1_z_xx_xxx_0[i] * fe_0 - ta1_z_xx_xxx_1[i] * fe_0 + ta_xxz_xxx_1[i] + ta1_z_xxz_xxx_0[i] * pa_z[i] - ta1_z_xxz_xxx_1[i] * pc_z[i];

        ta1_z_xxzz_xxy_0[i] =
            ta1_z_xx_xxy_0[i] * fe_0 - ta1_z_xx_xxy_1[i] * fe_0 + ta_xxz_xxy_1[i] + ta1_z_xxz_xxy_0[i] * pa_z[i] - ta1_z_xxz_xxy_1[i] * pc_z[i];

        ta1_z_xxzz_xxz_0[i] = ta1_z_zz_xxz_0[i] * fe_0 - ta1_z_zz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xz_1[i] * fe_0 +
                              ta1_z_xzz_xxz_0[i] * pa_x[i] - ta1_z_xzz_xxz_1[i] * pc_x[i];

        ta1_z_xxzz_xyy_0[i] =
            ta1_z_xx_xyy_0[i] * fe_0 - ta1_z_xx_xyy_1[i] * fe_0 + ta_xxz_xyy_1[i] + ta1_z_xxz_xyy_0[i] * pa_z[i] - ta1_z_xxz_xyy_1[i] * pc_z[i];

        ta1_z_xxzz_xyz_0[i] = ta1_z_zz_xyz_0[i] * fe_0 - ta1_z_zz_xyz_1[i] * fe_0 + ta1_z_xzz_yz_0[i] * fe_0 - ta1_z_xzz_yz_1[i] * fe_0 +
                              ta1_z_xzz_xyz_0[i] * pa_x[i] - ta1_z_xzz_xyz_1[i] * pc_x[i];

        ta1_z_xxzz_xzz_0[i] = ta1_z_zz_xzz_0[i] * fe_0 - ta1_z_zz_xzz_1[i] * fe_0 + ta1_z_xzz_zz_0[i] * fe_0 - ta1_z_xzz_zz_1[i] * fe_0 +
                              ta1_z_xzz_xzz_0[i] * pa_x[i] - ta1_z_xzz_xzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyy_0[i] = ta1_z_zz_yyy_0[i] * fe_0 - ta1_z_zz_yyy_1[i] * fe_0 + ta1_z_xzz_yyy_0[i] * pa_x[i] - ta1_z_xzz_yyy_1[i] * pc_x[i];

        ta1_z_xxzz_yyz_0[i] = ta1_z_zz_yyz_0[i] * fe_0 - ta1_z_zz_yyz_1[i] * fe_0 + ta1_z_xzz_yyz_0[i] * pa_x[i] - ta1_z_xzz_yyz_1[i] * pc_x[i];

        ta1_z_xxzz_yzz_0[i] = ta1_z_zz_yzz_0[i] * fe_0 - ta1_z_zz_yzz_1[i] * fe_0 + ta1_z_xzz_yzz_0[i] * pa_x[i] - ta1_z_xzz_yzz_1[i] * pc_x[i];

        ta1_z_xxzz_zzz_0[i] = ta1_z_zz_zzz_0[i] * fe_0 - ta1_z_zz_zzz_1[i] * fe_0 + ta1_z_xzz_zzz_0[i] * pa_x[i] - ta1_z_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 360-370 components of targeted buffer : GF

    auto ta1_z_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 360);

    auto ta1_z_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 361);

    auto ta1_z_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 362);

    auto ta1_z_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 363);

    auto ta1_z_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 364);

    auto ta1_z_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 365);

    auto ta1_z_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 366);

    auto ta1_z_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 367);

    auto ta1_z_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 368);

    auto ta1_z_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 369);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xyyy_xxx_0, \
                             ta1_z_xyyy_xxy_0, \
                             ta1_z_xyyy_xxz_0, \
                             ta1_z_xyyy_xyy_0, \
                             ta1_z_xyyy_xyz_0, \
                             ta1_z_xyyy_xzz_0, \
                             ta1_z_xyyy_yyy_0, \
                             ta1_z_xyyy_yyz_0, \
                             ta1_z_xyyy_yzz_0, \
                             ta1_z_xyyy_zzz_0, \
                             ta1_z_yyy_xx_0,   \
                             ta1_z_yyy_xx_1,   \
                             ta1_z_yyy_xxx_0,  \
                             ta1_z_yyy_xxx_1,  \
                             ta1_z_yyy_xxy_0,  \
                             ta1_z_yyy_xxy_1,  \
                             ta1_z_yyy_xxz_0,  \
                             ta1_z_yyy_xxz_1,  \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_xyy_0,  \
                             ta1_z_yyy_xyy_1,  \
                             ta1_z_yyy_xyz_0,  \
                             ta1_z_yyy_xyz_1,  \
                             ta1_z_yyy_xz_0,   \
                             ta1_z_yyy_xz_1,   \
                             ta1_z_yyy_xzz_0,  \
                             ta1_z_yyy_xzz_1,  \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yyy_0,  \
                             ta1_z_yyy_yyy_1,  \
                             ta1_z_yyy_yyz_0,  \
                             ta1_z_yyy_yyz_1,  \
                             ta1_z_yyy_yz_0,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_yzz_0,  \
                             ta1_z_yyy_yzz_1,  \
                             ta1_z_yyy_zz_0,   \
                             ta1_z_yyy_zz_1,   \
                             ta1_z_yyy_zzz_0,  \
                             ta1_z_yyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_xxx_0[i] =
            3.0 * ta1_z_yyy_xx_0[i] * fe_0 - 3.0 * ta1_z_yyy_xx_1[i] * fe_0 + ta1_z_yyy_xxx_0[i] * pa_x[i] - ta1_z_yyy_xxx_1[i] * pc_x[i];

        ta1_z_xyyy_xxy_0[i] =
            2.0 * ta1_z_yyy_xy_0[i] * fe_0 - 2.0 * ta1_z_yyy_xy_1[i] * fe_0 + ta1_z_yyy_xxy_0[i] * pa_x[i] - ta1_z_yyy_xxy_1[i] * pc_x[i];

        ta1_z_xyyy_xxz_0[i] =
            2.0 * ta1_z_yyy_xz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xz_1[i] * fe_0 + ta1_z_yyy_xxz_0[i] * pa_x[i] - ta1_z_yyy_xxz_1[i] * pc_x[i];

        ta1_z_xyyy_xyy_0[i] = ta1_z_yyy_yy_0[i] * fe_0 - ta1_z_yyy_yy_1[i] * fe_0 + ta1_z_yyy_xyy_0[i] * pa_x[i] - ta1_z_yyy_xyy_1[i] * pc_x[i];

        ta1_z_xyyy_xyz_0[i] = ta1_z_yyy_yz_0[i] * fe_0 - ta1_z_yyy_yz_1[i] * fe_0 + ta1_z_yyy_xyz_0[i] * pa_x[i] - ta1_z_yyy_xyz_1[i] * pc_x[i];

        ta1_z_xyyy_xzz_0[i] = ta1_z_yyy_zz_0[i] * fe_0 - ta1_z_yyy_zz_1[i] * fe_0 + ta1_z_yyy_xzz_0[i] * pa_x[i] - ta1_z_yyy_xzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyy_0[i] = ta1_z_yyy_yyy_0[i] * pa_x[i] - ta1_z_yyy_yyy_1[i] * pc_x[i];

        ta1_z_xyyy_yyz_0[i] = ta1_z_yyy_yyz_0[i] * pa_x[i] - ta1_z_yyy_yyz_1[i] * pc_x[i];

        ta1_z_xyyy_yzz_0[i] = ta1_z_yyy_yzz_0[i] * pa_x[i] - ta1_z_yyy_yzz_1[i] * pc_x[i];

        ta1_z_xyyy_zzz_0[i] = ta1_z_yyy_zzz_0[i] * pa_x[i] - ta1_z_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 370-380 components of targeted buffer : GF

    auto ta1_z_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 370);

    auto ta1_z_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 371);

    auto ta1_z_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 372);

    auto ta1_z_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 373);

    auto ta1_z_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 374);

    auto ta1_z_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 375);

    auto ta1_z_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 376);

    auto ta1_z_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 377);

    auto ta1_z_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 378);

    auto ta1_z_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 379);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xyy_xxx_0,  \
                             ta1_z_xyy_xxx_1,  \
                             ta1_z_xyy_xxy_0,  \
                             ta1_z_xyy_xxy_1,  \
                             ta1_z_xyy_xyy_0,  \
                             ta1_z_xyy_xyy_1,  \
                             ta1_z_xyyz_xxx_0, \
                             ta1_z_xyyz_xxy_0, \
                             ta1_z_xyyz_xxz_0, \
                             ta1_z_xyyz_xyy_0, \
                             ta1_z_xyyz_xyz_0, \
                             ta1_z_xyyz_xzz_0, \
                             ta1_z_xyyz_yyy_0, \
                             ta1_z_xyyz_yyz_0, \
                             ta1_z_xyyz_yzz_0, \
                             ta1_z_xyyz_zzz_0, \
                             ta1_z_yyz_xxz_0,  \
                             ta1_z_yyz_xxz_1,  \
                             ta1_z_yyz_xyz_0,  \
                             ta1_z_yyz_xyz_1,  \
                             ta1_z_yyz_xz_0,   \
                             ta1_z_yyz_xz_1,   \
                             ta1_z_yyz_xzz_0,  \
                             ta1_z_yyz_xzz_1,  \
                             ta1_z_yyz_yyy_0,  \
                             ta1_z_yyz_yyy_1,  \
                             ta1_z_yyz_yyz_0,  \
                             ta1_z_yyz_yyz_1,  \
                             ta1_z_yyz_yz_0,   \
                             ta1_z_yyz_yz_1,   \
                             ta1_z_yyz_yzz_0,  \
                             ta1_z_yyz_yzz_1,  \
                             ta1_z_yyz_zz_0,   \
                             ta1_z_yyz_zz_1,   \
                             ta1_z_yyz_zzz_0,  \
                             ta1_z_yyz_zzz_1,  \
                             ta_xyy_xxx_1,     \
                             ta_xyy_xxy_1,     \
                             ta_xyy_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyz_xxx_0[i] = ta_xyy_xxx_1[i] + ta1_z_xyy_xxx_0[i] * pa_z[i] - ta1_z_xyy_xxx_1[i] * pc_z[i];

        ta1_z_xyyz_xxy_0[i] = ta_xyy_xxy_1[i] + ta1_z_xyy_xxy_0[i] * pa_z[i] - ta1_z_xyy_xxy_1[i] * pc_z[i];

        ta1_z_xyyz_xxz_0[i] =
            2.0 * ta1_z_yyz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xz_1[i] * fe_0 + ta1_z_yyz_xxz_0[i] * pa_x[i] - ta1_z_yyz_xxz_1[i] * pc_x[i];

        ta1_z_xyyz_xyy_0[i] = ta_xyy_xyy_1[i] + ta1_z_xyy_xyy_0[i] * pa_z[i] - ta1_z_xyy_xyy_1[i] * pc_z[i];

        ta1_z_xyyz_xyz_0[i] = ta1_z_yyz_yz_0[i] * fe_0 - ta1_z_yyz_yz_1[i] * fe_0 + ta1_z_yyz_xyz_0[i] * pa_x[i] - ta1_z_yyz_xyz_1[i] * pc_x[i];

        ta1_z_xyyz_xzz_0[i] = ta1_z_yyz_zz_0[i] * fe_0 - ta1_z_yyz_zz_1[i] * fe_0 + ta1_z_yyz_xzz_0[i] * pa_x[i] - ta1_z_yyz_xzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyy_0[i] = ta1_z_yyz_yyy_0[i] * pa_x[i] - ta1_z_yyz_yyy_1[i] * pc_x[i];

        ta1_z_xyyz_yyz_0[i] = ta1_z_yyz_yyz_0[i] * pa_x[i] - ta1_z_yyz_yyz_1[i] * pc_x[i];

        ta1_z_xyyz_yzz_0[i] = ta1_z_yyz_yzz_0[i] * pa_x[i] - ta1_z_yyz_yzz_1[i] * pc_x[i];

        ta1_z_xyyz_zzz_0[i] = ta1_z_yyz_zzz_0[i] * pa_x[i] - ta1_z_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 380-390 components of targeted buffer : GF

    auto ta1_z_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 380);

    auto ta1_z_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 381);

    auto ta1_z_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 382);

    auto ta1_z_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 383);

    auto ta1_z_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 384);

    auto ta1_z_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 385);

    auto ta1_z_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 386);

    auto ta1_z_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 387);

    auto ta1_z_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 388);

    auto ta1_z_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 389);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xyzz_xxx_0, \
                             ta1_z_xyzz_xxy_0, \
                             ta1_z_xyzz_xxz_0, \
                             ta1_z_xyzz_xyy_0, \
                             ta1_z_xyzz_xyz_0, \
                             ta1_z_xyzz_xzz_0, \
                             ta1_z_xyzz_yyy_0, \
                             ta1_z_xyzz_yyz_0, \
                             ta1_z_xyzz_yzz_0, \
                             ta1_z_xyzz_zzz_0, \
                             ta1_z_xzz_xxx_0,  \
                             ta1_z_xzz_xxx_1,  \
                             ta1_z_xzz_xxz_0,  \
                             ta1_z_xzz_xxz_1,  \
                             ta1_z_xzz_xzz_0,  \
                             ta1_z_xzz_xzz_1,  \
                             ta1_z_yzz_xxy_0,  \
                             ta1_z_yzz_xxy_1,  \
                             ta1_z_yzz_xy_0,   \
                             ta1_z_yzz_xy_1,   \
                             ta1_z_yzz_xyy_0,  \
                             ta1_z_yzz_xyy_1,  \
                             ta1_z_yzz_xyz_0,  \
                             ta1_z_yzz_xyz_1,  \
                             ta1_z_yzz_yy_0,   \
                             ta1_z_yzz_yy_1,   \
                             ta1_z_yzz_yyy_0,  \
                             ta1_z_yzz_yyy_1,  \
                             ta1_z_yzz_yyz_0,  \
                             ta1_z_yzz_yyz_1,  \
                             ta1_z_yzz_yz_0,   \
                             ta1_z_yzz_yz_1,   \
                             ta1_z_yzz_yzz_0,  \
                             ta1_z_yzz_yzz_1,  \
                             ta1_z_yzz_zzz_0,  \
                             ta1_z_yzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzz_xxx_0[i] = ta1_z_xzz_xxx_0[i] * pa_y[i] - ta1_z_xzz_xxx_1[i] * pc_y[i];

        ta1_z_xyzz_xxy_0[i] =
            2.0 * ta1_z_yzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yzz_xy_1[i] * fe_0 + ta1_z_yzz_xxy_0[i] * pa_x[i] - ta1_z_yzz_xxy_1[i] * pc_x[i];

        ta1_z_xyzz_xxz_0[i] = ta1_z_xzz_xxz_0[i] * pa_y[i] - ta1_z_xzz_xxz_1[i] * pc_y[i];

        ta1_z_xyzz_xyy_0[i] = ta1_z_yzz_yy_0[i] * fe_0 - ta1_z_yzz_yy_1[i] * fe_0 + ta1_z_yzz_xyy_0[i] * pa_x[i] - ta1_z_yzz_xyy_1[i] * pc_x[i];

        ta1_z_xyzz_xyz_0[i] = ta1_z_yzz_yz_0[i] * fe_0 - ta1_z_yzz_yz_1[i] * fe_0 + ta1_z_yzz_xyz_0[i] * pa_x[i] - ta1_z_yzz_xyz_1[i] * pc_x[i];

        ta1_z_xyzz_xzz_0[i] = ta1_z_xzz_xzz_0[i] * pa_y[i] - ta1_z_xzz_xzz_1[i] * pc_y[i];

        ta1_z_xyzz_yyy_0[i] = ta1_z_yzz_yyy_0[i] * pa_x[i] - ta1_z_yzz_yyy_1[i] * pc_x[i];

        ta1_z_xyzz_yyz_0[i] = ta1_z_yzz_yyz_0[i] * pa_x[i] - ta1_z_yzz_yyz_1[i] * pc_x[i];

        ta1_z_xyzz_yzz_0[i] = ta1_z_yzz_yzz_0[i] * pa_x[i] - ta1_z_yzz_yzz_1[i] * pc_x[i];

        ta1_z_xyzz_zzz_0[i] = ta1_z_yzz_zzz_0[i] * pa_x[i] - ta1_z_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 390-400 components of targeted buffer : GF

    auto ta1_z_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 390);

    auto ta1_z_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 391);

    auto ta1_z_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 392);

    auto ta1_z_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 393);

    auto ta1_z_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 394);

    auto ta1_z_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 395);

    auto ta1_z_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 396);

    auto ta1_z_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 397);

    auto ta1_z_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 398);

    auto ta1_z_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 399);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xzzz_xxx_0, \
                             ta1_z_xzzz_xxy_0, \
                             ta1_z_xzzz_xxz_0, \
                             ta1_z_xzzz_xyy_0, \
                             ta1_z_xzzz_xyz_0, \
                             ta1_z_xzzz_xzz_0, \
                             ta1_z_xzzz_yyy_0, \
                             ta1_z_xzzz_yyz_0, \
                             ta1_z_xzzz_yzz_0, \
                             ta1_z_xzzz_zzz_0, \
                             ta1_z_zzz_xx_0,   \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xxx_0,  \
                             ta1_z_zzz_xxx_1,  \
                             ta1_z_zzz_xxy_0,  \
                             ta1_z_zzz_xxy_1,  \
                             ta1_z_zzz_xxz_0,  \
                             ta1_z_zzz_xxz_1,  \
                             ta1_z_zzz_xy_0,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xyy_0,  \
                             ta1_z_zzz_xyy_1,  \
                             ta1_z_zzz_xyz_0,  \
                             ta1_z_zzz_xyz_1,  \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_xzz_0,  \
                             ta1_z_zzz_xzz_1,  \
                             ta1_z_zzz_yy_0,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yyy_0,  \
                             ta1_z_zzz_yyy_1,  \
                             ta1_z_zzz_yyz_0,  \
                             ta1_z_zzz_yyz_1,  \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_yzz_0,  \
                             ta1_z_zzz_yzz_1,  \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta1_z_zzz_zzz_0,  \
                             ta1_z_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_xxx_0[i] =
            3.0 * ta1_z_zzz_xx_0[i] * fe_0 - 3.0 * ta1_z_zzz_xx_1[i] * fe_0 + ta1_z_zzz_xxx_0[i] * pa_x[i] - ta1_z_zzz_xxx_1[i] * pc_x[i];

        ta1_z_xzzz_xxy_0[i] =
            2.0 * ta1_z_zzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xy_1[i] * fe_0 + ta1_z_zzz_xxy_0[i] * pa_x[i] - ta1_z_zzz_xxy_1[i] * pc_x[i];

        ta1_z_xzzz_xxz_0[i] =
            2.0 * ta1_z_zzz_xz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xz_1[i] * fe_0 + ta1_z_zzz_xxz_0[i] * pa_x[i] - ta1_z_zzz_xxz_1[i] * pc_x[i];

        ta1_z_xzzz_xyy_0[i] = ta1_z_zzz_yy_0[i] * fe_0 - ta1_z_zzz_yy_1[i] * fe_0 + ta1_z_zzz_xyy_0[i] * pa_x[i] - ta1_z_zzz_xyy_1[i] * pc_x[i];

        ta1_z_xzzz_xyz_0[i] = ta1_z_zzz_yz_0[i] * fe_0 - ta1_z_zzz_yz_1[i] * fe_0 + ta1_z_zzz_xyz_0[i] * pa_x[i] - ta1_z_zzz_xyz_1[i] * pc_x[i];

        ta1_z_xzzz_xzz_0[i] = ta1_z_zzz_zz_0[i] * fe_0 - ta1_z_zzz_zz_1[i] * fe_0 + ta1_z_zzz_xzz_0[i] * pa_x[i] - ta1_z_zzz_xzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyy_0[i] = ta1_z_zzz_yyy_0[i] * pa_x[i] - ta1_z_zzz_yyy_1[i] * pc_x[i];

        ta1_z_xzzz_yyz_0[i] = ta1_z_zzz_yyz_0[i] * pa_x[i] - ta1_z_zzz_yyz_1[i] * pc_x[i];

        ta1_z_xzzz_yzz_0[i] = ta1_z_zzz_yzz_0[i] * pa_x[i] - ta1_z_zzz_yzz_1[i] * pc_x[i];

        ta1_z_xzzz_zzz_0[i] = ta1_z_zzz_zzz_0[i] * pa_x[i] - ta1_z_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 400-410 components of targeted buffer : GF

    auto ta1_z_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 400);

    auto ta1_z_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 401);

    auto ta1_z_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 402);

    auto ta1_z_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 403);

    auto ta1_z_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 404);

    auto ta1_z_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 405);

    auto ta1_z_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 406);

    auto ta1_z_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 407);

    auto ta1_z_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 408);

    auto ta1_z_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 409);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yy_xxx_0,   \
                             ta1_z_yy_xxx_1,   \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xxz_0,   \
                             ta1_z_yy_xxz_1,   \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_xyz_0,   \
                             ta1_z_yy_xyz_1,   \
                             ta1_z_yy_xzz_0,   \
                             ta1_z_yy_xzz_1,   \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yy_yyz_0,   \
                             ta1_z_yy_yyz_1,   \
                             ta1_z_yy_yzz_0,   \
                             ta1_z_yy_yzz_1,   \
                             ta1_z_yy_zzz_0,   \
                             ta1_z_yy_zzz_1,   \
                             ta1_z_yyy_xx_0,   \
                             ta1_z_yyy_xx_1,   \
                             ta1_z_yyy_xxx_0,  \
                             ta1_z_yyy_xxx_1,  \
                             ta1_z_yyy_xxy_0,  \
                             ta1_z_yyy_xxy_1,  \
                             ta1_z_yyy_xxz_0,  \
                             ta1_z_yyy_xxz_1,  \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_xyy_0,  \
                             ta1_z_yyy_xyy_1,  \
                             ta1_z_yyy_xyz_0,  \
                             ta1_z_yyy_xyz_1,  \
                             ta1_z_yyy_xz_0,   \
                             ta1_z_yyy_xz_1,   \
                             ta1_z_yyy_xzz_0,  \
                             ta1_z_yyy_xzz_1,  \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yyy_0,  \
                             ta1_z_yyy_yyy_1,  \
                             ta1_z_yyy_yyz_0,  \
                             ta1_z_yyy_yyz_1,  \
                             ta1_z_yyy_yz_0,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_yzz_0,  \
                             ta1_z_yyy_yzz_1,  \
                             ta1_z_yyy_zz_0,   \
                             ta1_z_yyy_zz_1,   \
                             ta1_z_yyy_zzz_0,  \
                             ta1_z_yyy_zzz_1,  \
                             ta1_z_yyyy_xxx_0, \
                             ta1_z_yyyy_xxy_0, \
                             ta1_z_yyyy_xxz_0, \
                             ta1_z_yyyy_xyy_0, \
                             ta1_z_yyyy_xyz_0, \
                             ta1_z_yyyy_xzz_0, \
                             ta1_z_yyyy_yyy_0, \
                             ta1_z_yyyy_yyz_0, \
                             ta1_z_yyyy_yzz_0, \
                             ta1_z_yyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_xxx_0[i] =
            3.0 * ta1_z_yy_xxx_0[i] * fe_0 - 3.0 * ta1_z_yy_xxx_1[i] * fe_0 + ta1_z_yyy_xxx_0[i] * pa_y[i] - ta1_z_yyy_xxx_1[i] * pc_y[i];

        ta1_z_yyyy_xxy_0[i] = 3.0 * ta1_z_yy_xxy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxy_1[i] * fe_0 + ta1_z_yyy_xx_0[i] * fe_0 - ta1_z_yyy_xx_1[i] * fe_0 +
                              ta1_z_yyy_xxy_0[i] * pa_y[i] - ta1_z_yyy_xxy_1[i] * pc_y[i];

        ta1_z_yyyy_xxz_0[i] =
            3.0 * ta1_z_yy_xxz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxz_1[i] * fe_0 + ta1_z_yyy_xxz_0[i] * pa_y[i] - ta1_z_yyy_xxz_1[i] * pc_y[i];

        ta1_z_yyyy_xyy_0[i] = 3.0 * ta1_z_yy_xyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xy_0[i] * fe_0 -
                              2.0 * ta1_z_yyy_xy_1[i] * fe_0 + ta1_z_yyy_xyy_0[i] * pa_y[i] - ta1_z_yyy_xyy_1[i] * pc_y[i];

        ta1_z_yyyy_xyz_0[i] = 3.0 * ta1_z_yy_xyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyz_1[i] * fe_0 + ta1_z_yyy_xz_0[i] * fe_0 - ta1_z_yyy_xz_1[i] * fe_0 +
                              ta1_z_yyy_xyz_0[i] * pa_y[i] - ta1_z_yyy_xyz_1[i] * pc_y[i];

        ta1_z_yyyy_xzz_0[i] =
            3.0 * ta1_z_yy_xzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xzz_1[i] * fe_0 + ta1_z_yyy_xzz_0[i] * pa_y[i] - ta1_z_yyy_xzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyy_0[i] = 3.0 * ta1_z_yy_yyy_0[i] * fe_0 - 3.0 * ta1_z_yy_yyy_1[i] * fe_0 + 3.0 * ta1_z_yyy_yy_0[i] * fe_0 -
                              3.0 * ta1_z_yyy_yy_1[i] * fe_0 + ta1_z_yyy_yyy_0[i] * pa_y[i] - ta1_z_yyy_yyy_1[i] * pc_y[i];

        ta1_z_yyyy_yyz_0[i] = 3.0 * ta1_z_yy_yyz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yz_0[i] * fe_0 -
                              2.0 * ta1_z_yyy_yz_1[i] * fe_0 + ta1_z_yyy_yyz_0[i] * pa_y[i] - ta1_z_yyy_yyz_1[i] * pc_y[i];

        ta1_z_yyyy_yzz_0[i] = 3.0 * ta1_z_yy_yzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yzz_1[i] * fe_0 + ta1_z_yyy_zz_0[i] * fe_0 - ta1_z_yyy_zz_1[i] * fe_0 +
                              ta1_z_yyy_yzz_0[i] * pa_y[i] - ta1_z_yyy_yzz_1[i] * pc_y[i];

        ta1_z_yyyy_zzz_0[i] =
            3.0 * ta1_z_yy_zzz_0[i] * fe_0 - 3.0 * ta1_z_yy_zzz_1[i] * fe_0 + ta1_z_yyy_zzz_0[i] * pa_y[i] - ta1_z_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 410-420 components of targeted buffer : GF

    auto ta1_z_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 410);

    auto ta1_z_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 411);

    auto ta1_z_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 412);

    auto ta1_z_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 413);

    auto ta1_z_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 414);

    auto ta1_z_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 415);

    auto ta1_z_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 416);

    auto ta1_z_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 417);

    auto ta1_z_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 418);

    auto ta1_z_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 419);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyy_xxx_0,  \
                             ta1_z_yyy_xxx_1,  \
                             ta1_z_yyy_xxy_0,  \
                             ta1_z_yyy_xxy_1,  \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_xyy_0,  \
                             ta1_z_yyy_xyy_1,  \
                             ta1_z_yyy_xyz_0,  \
                             ta1_z_yyy_xyz_1,  \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yyy_0,  \
                             ta1_z_yyy_yyy_1,  \
                             ta1_z_yyy_yyz_0,  \
                             ta1_z_yyy_yyz_1,  \
                             ta1_z_yyy_yz_0,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_yzz_0,  \
                             ta1_z_yyy_yzz_1,  \
                             ta1_z_yyyz_xxx_0, \
                             ta1_z_yyyz_xxy_0, \
                             ta1_z_yyyz_xxz_0, \
                             ta1_z_yyyz_xyy_0, \
                             ta1_z_yyyz_xyz_0, \
                             ta1_z_yyyz_xzz_0, \
                             ta1_z_yyyz_yyy_0, \
                             ta1_z_yyyz_yyz_0, \
                             ta1_z_yyyz_yzz_0, \
                             ta1_z_yyyz_zzz_0, \
                             ta1_z_yyz_xxz_0,  \
                             ta1_z_yyz_xxz_1,  \
                             ta1_z_yyz_xzz_0,  \
                             ta1_z_yyz_xzz_1,  \
                             ta1_z_yyz_zzz_0,  \
                             ta1_z_yyz_zzz_1,  \
                             ta1_z_yz_xxz_0,   \
                             ta1_z_yz_xxz_1,   \
                             ta1_z_yz_xzz_0,   \
                             ta1_z_yz_xzz_1,   \
                             ta1_z_yz_zzz_0,   \
                             ta1_z_yz_zzz_1,   \
                             ta_yyy_xxx_1,     \
                             ta_yyy_xxy_1,     \
                             ta_yyy_xyy_1,     \
                             ta_yyy_xyz_1,     \
                             ta_yyy_yyy_1,     \
                             ta_yyy_yyz_1,     \
                             ta_yyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_xxx_0[i] = ta_yyy_xxx_1[i] + ta1_z_yyy_xxx_0[i] * pa_z[i] - ta1_z_yyy_xxx_1[i] * pc_z[i];

        ta1_z_yyyz_xxy_0[i] = ta_yyy_xxy_1[i] + ta1_z_yyy_xxy_0[i] * pa_z[i] - ta1_z_yyy_xxy_1[i] * pc_z[i];

        ta1_z_yyyz_xxz_0[i] =
            2.0 * ta1_z_yz_xxz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxz_1[i] * fe_0 + ta1_z_yyz_xxz_0[i] * pa_y[i] - ta1_z_yyz_xxz_1[i] * pc_y[i];

        ta1_z_yyyz_xyy_0[i] = ta_yyy_xyy_1[i] + ta1_z_yyy_xyy_0[i] * pa_z[i] - ta1_z_yyy_xyy_1[i] * pc_z[i];

        ta1_z_yyyz_xyz_0[i] =
            ta1_z_yyy_xy_0[i] * fe_0 - ta1_z_yyy_xy_1[i] * fe_0 + ta_yyy_xyz_1[i] + ta1_z_yyy_xyz_0[i] * pa_z[i] - ta1_z_yyy_xyz_1[i] * pc_z[i];

        ta1_z_yyyz_xzz_0[i] =
            2.0 * ta1_z_yz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xzz_1[i] * fe_0 + ta1_z_yyz_xzz_0[i] * pa_y[i] - ta1_z_yyz_xzz_1[i] * pc_y[i];

        ta1_z_yyyz_yyy_0[i] = ta_yyy_yyy_1[i] + ta1_z_yyy_yyy_0[i] * pa_z[i] - ta1_z_yyy_yyy_1[i] * pc_z[i];

        ta1_z_yyyz_yyz_0[i] =
            ta1_z_yyy_yy_0[i] * fe_0 - ta1_z_yyy_yy_1[i] * fe_0 + ta_yyy_yyz_1[i] + ta1_z_yyy_yyz_0[i] * pa_z[i] - ta1_z_yyy_yyz_1[i] * pc_z[i];

        ta1_z_yyyz_yzz_0[i] = 2.0 * ta1_z_yyy_yz_0[i] * fe_0 - 2.0 * ta1_z_yyy_yz_1[i] * fe_0 + ta_yyy_yzz_1[i] + ta1_z_yyy_yzz_0[i] * pa_z[i] -
                              ta1_z_yyy_yzz_1[i] * pc_z[i];

        ta1_z_yyyz_zzz_0[i] =
            2.0 * ta1_z_yz_zzz_0[i] * fe_0 - 2.0 * ta1_z_yz_zzz_1[i] * fe_0 + ta1_z_yyz_zzz_0[i] * pa_y[i] - ta1_z_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 420-430 components of targeted buffer : GF

    auto ta1_z_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 420);

    auto ta1_z_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 421);

    auto ta1_z_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 422);

    auto ta1_z_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 423);

    auto ta1_z_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 424);

    auto ta1_z_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 425);

    auto ta1_z_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 426);

    auto ta1_z_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 427);

    auto ta1_z_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 428);

    auto ta1_z_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 429);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yyz_xxy_0,  \
                             ta1_z_yyz_xxy_1,  \
                             ta1_z_yyz_xyy_0,  \
                             ta1_z_yyz_xyy_1,  \
                             ta1_z_yyz_yyy_0,  \
                             ta1_z_yyz_yyy_1,  \
                             ta1_z_yyzz_xxx_0, \
                             ta1_z_yyzz_xxy_0, \
                             ta1_z_yyzz_xxz_0, \
                             ta1_z_yyzz_xyy_0, \
                             ta1_z_yyzz_xyz_0, \
                             ta1_z_yyzz_xzz_0, \
                             ta1_z_yyzz_yyy_0, \
                             ta1_z_yyzz_yyz_0, \
                             ta1_z_yyzz_yzz_0, \
                             ta1_z_yyzz_zzz_0, \
                             ta1_z_yzz_xxx_0,  \
                             ta1_z_yzz_xxx_1,  \
                             ta1_z_yzz_xxz_0,  \
                             ta1_z_yzz_xxz_1,  \
                             ta1_z_yzz_xyz_0,  \
                             ta1_z_yzz_xyz_1,  \
                             ta1_z_yzz_xz_0,   \
                             ta1_z_yzz_xz_1,   \
                             ta1_z_yzz_xzz_0,  \
                             ta1_z_yzz_xzz_1,  \
                             ta1_z_yzz_yyz_0,  \
                             ta1_z_yzz_yyz_1,  \
                             ta1_z_yzz_yz_0,   \
                             ta1_z_yzz_yz_1,   \
                             ta1_z_yzz_yzz_0,  \
                             ta1_z_yzz_yzz_1,  \
                             ta1_z_yzz_zz_0,   \
                             ta1_z_yzz_zz_1,   \
                             ta1_z_yzz_zzz_0,  \
                             ta1_z_yzz_zzz_1,  \
                             ta1_z_zz_xxx_0,   \
                             ta1_z_zz_xxx_1,   \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta_yyz_xxy_1,     \
                             ta_yyz_xyy_1,     \
                             ta_yyz_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_xxx_0[i] = ta1_z_zz_xxx_0[i] * fe_0 - ta1_z_zz_xxx_1[i] * fe_0 + ta1_z_yzz_xxx_0[i] * pa_y[i] - ta1_z_yzz_xxx_1[i] * pc_y[i];

        ta1_z_yyzz_xxy_0[i] =
            ta1_z_yy_xxy_0[i] * fe_0 - ta1_z_yy_xxy_1[i] * fe_0 + ta_yyz_xxy_1[i] + ta1_z_yyz_xxy_0[i] * pa_z[i] - ta1_z_yyz_xxy_1[i] * pc_z[i];

        ta1_z_yyzz_xxz_0[i] = ta1_z_zz_xxz_0[i] * fe_0 - ta1_z_zz_xxz_1[i] * fe_0 + ta1_z_yzz_xxz_0[i] * pa_y[i] - ta1_z_yzz_xxz_1[i] * pc_y[i];

        ta1_z_yyzz_xyy_0[i] =
            ta1_z_yy_xyy_0[i] * fe_0 - ta1_z_yy_xyy_1[i] * fe_0 + ta_yyz_xyy_1[i] + ta1_z_yyz_xyy_0[i] * pa_z[i] - ta1_z_yyz_xyy_1[i] * pc_z[i];

        ta1_z_yyzz_xyz_0[i] = ta1_z_zz_xyz_0[i] * fe_0 - ta1_z_zz_xyz_1[i] * fe_0 + ta1_z_yzz_xz_0[i] * fe_0 - ta1_z_yzz_xz_1[i] * fe_0 +
                              ta1_z_yzz_xyz_0[i] * pa_y[i] - ta1_z_yzz_xyz_1[i] * pc_y[i];

        ta1_z_yyzz_xzz_0[i] = ta1_z_zz_xzz_0[i] * fe_0 - ta1_z_zz_xzz_1[i] * fe_0 + ta1_z_yzz_xzz_0[i] * pa_y[i] - ta1_z_yzz_xzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyy_0[i] =
            ta1_z_yy_yyy_0[i] * fe_0 - ta1_z_yy_yyy_1[i] * fe_0 + ta_yyz_yyy_1[i] + ta1_z_yyz_yyy_0[i] * pa_z[i] - ta1_z_yyz_yyy_1[i] * pc_z[i];

        ta1_z_yyzz_yyz_0[i] = ta1_z_zz_yyz_0[i] * fe_0 - ta1_z_zz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yzz_yz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yz_1[i] * fe_0 +
                              ta1_z_yzz_yyz_0[i] * pa_y[i] - ta1_z_yzz_yyz_1[i] * pc_y[i];

        ta1_z_yyzz_yzz_0[i] = ta1_z_zz_yzz_0[i] * fe_0 - ta1_z_zz_yzz_1[i] * fe_0 + ta1_z_yzz_zz_0[i] * fe_0 - ta1_z_yzz_zz_1[i] * fe_0 +
                              ta1_z_yzz_yzz_0[i] * pa_y[i] - ta1_z_yzz_yzz_1[i] * pc_y[i];

        ta1_z_yyzz_zzz_0[i] = ta1_z_zz_zzz_0[i] * fe_0 - ta1_z_zz_zzz_1[i] * fe_0 + ta1_z_yzz_zzz_0[i] * pa_y[i] - ta1_z_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 430-440 components of targeted buffer : GF

    auto ta1_z_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 430);

    auto ta1_z_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 431);

    auto ta1_z_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 432);

    auto ta1_z_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 433);

    auto ta1_z_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 434);

    auto ta1_z_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 435);

    auto ta1_z_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 436);

    auto ta1_z_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 437);

    auto ta1_z_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 438);

    auto ta1_z_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 439);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yzzz_xxx_0, \
                             ta1_z_yzzz_xxy_0, \
                             ta1_z_yzzz_xxz_0, \
                             ta1_z_yzzz_xyy_0, \
                             ta1_z_yzzz_xyz_0, \
                             ta1_z_yzzz_xzz_0, \
                             ta1_z_yzzz_yyy_0, \
                             ta1_z_yzzz_yyz_0, \
                             ta1_z_yzzz_yzz_0, \
                             ta1_z_yzzz_zzz_0, \
                             ta1_z_zzz_xx_0,   \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xxx_0,  \
                             ta1_z_zzz_xxx_1,  \
                             ta1_z_zzz_xxy_0,  \
                             ta1_z_zzz_xxy_1,  \
                             ta1_z_zzz_xxz_0,  \
                             ta1_z_zzz_xxz_1,  \
                             ta1_z_zzz_xy_0,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xyy_0,  \
                             ta1_z_zzz_xyy_1,  \
                             ta1_z_zzz_xyz_0,  \
                             ta1_z_zzz_xyz_1,  \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_xzz_0,  \
                             ta1_z_zzz_xzz_1,  \
                             ta1_z_zzz_yy_0,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yyy_0,  \
                             ta1_z_zzz_yyy_1,  \
                             ta1_z_zzz_yyz_0,  \
                             ta1_z_zzz_yyz_1,  \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_yzz_0,  \
                             ta1_z_zzz_yzz_1,  \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta1_z_zzz_zzz_0,  \
                             ta1_z_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_xxx_0[i] = ta1_z_zzz_xxx_0[i] * pa_y[i] - ta1_z_zzz_xxx_1[i] * pc_y[i];

        ta1_z_yzzz_xxy_0[i] = ta1_z_zzz_xx_0[i] * fe_0 - ta1_z_zzz_xx_1[i] * fe_0 + ta1_z_zzz_xxy_0[i] * pa_y[i] - ta1_z_zzz_xxy_1[i] * pc_y[i];

        ta1_z_yzzz_xxz_0[i] = ta1_z_zzz_xxz_0[i] * pa_y[i] - ta1_z_zzz_xxz_1[i] * pc_y[i];

        ta1_z_yzzz_xyy_0[i] =
            2.0 * ta1_z_zzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xy_1[i] * fe_0 + ta1_z_zzz_xyy_0[i] * pa_y[i] - ta1_z_zzz_xyy_1[i] * pc_y[i];

        ta1_z_yzzz_xyz_0[i] = ta1_z_zzz_xz_0[i] * fe_0 - ta1_z_zzz_xz_1[i] * fe_0 + ta1_z_zzz_xyz_0[i] * pa_y[i] - ta1_z_zzz_xyz_1[i] * pc_y[i];

        ta1_z_yzzz_xzz_0[i] = ta1_z_zzz_xzz_0[i] * pa_y[i] - ta1_z_zzz_xzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyy_0[i] =
            3.0 * ta1_z_zzz_yy_0[i] * fe_0 - 3.0 * ta1_z_zzz_yy_1[i] * fe_0 + ta1_z_zzz_yyy_0[i] * pa_y[i] - ta1_z_zzz_yyy_1[i] * pc_y[i];

        ta1_z_yzzz_yyz_0[i] =
            2.0 * ta1_z_zzz_yz_0[i] * fe_0 - 2.0 * ta1_z_zzz_yz_1[i] * fe_0 + ta1_z_zzz_yyz_0[i] * pa_y[i] - ta1_z_zzz_yyz_1[i] * pc_y[i];

        ta1_z_yzzz_yzz_0[i] = ta1_z_zzz_zz_0[i] * fe_0 - ta1_z_zzz_zz_1[i] * fe_0 + ta1_z_zzz_yzz_0[i] * pa_y[i] - ta1_z_zzz_yzz_1[i] * pc_y[i];

        ta1_z_yzzz_zzz_0[i] = ta1_z_zzz_zzz_0[i] * pa_y[i] - ta1_z_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 440-450 components of targeted buffer : GF

    auto ta1_z_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 440);

    auto ta1_z_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 441);

    auto ta1_z_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 442);

    auto ta1_z_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 443);

    auto ta1_z_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 444);

    auto ta1_z_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 445);

    auto ta1_z_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 446);

    auto ta1_z_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 447);

    auto ta1_z_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 448);

    auto ta1_z_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 449);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_zz_xxx_0,   \
                             ta1_z_zz_xxx_1,   \
                             ta1_z_zz_xxy_0,   \
                             ta1_z_zz_xxy_1,   \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xyy_0,   \
                             ta1_z_zz_xyy_1,   \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_yyy_0,   \
                             ta1_z_zz_yyy_1,   \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta1_z_zzz_xx_0,   \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xxx_0,  \
                             ta1_z_zzz_xxx_1,  \
                             ta1_z_zzz_xxy_0,  \
                             ta1_z_zzz_xxy_1,  \
                             ta1_z_zzz_xxz_0,  \
                             ta1_z_zzz_xxz_1,  \
                             ta1_z_zzz_xy_0,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xyy_0,  \
                             ta1_z_zzz_xyy_1,  \
                             ta1_z_zzz_xyz_0,  \
                             ta1_z_zzz_xyz_1,  \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_xzz_0,  \
                             ta1_z_zzz_xzz_1,  \
                             ta1_z_zzz_yy_0,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yyy_0,  \
                             ta1_z_zzz_yyy_1,  \
                             ta1_z_zzz_yyz_0,  \
                             ta1_z_zzz_yyz_1,  \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_yzz_0,  \
                             ta1_z_zzz_yzz_1,  \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta1_z_zzz_zzz_0,  \
                             ta1_z_zzz_zzz_1,  \
                             ta1_z_zzzz_xxx_0, \
                             ta1_z_zzzz_xxy_0, \
                             ta1_z_zzzz_xxz_0, \
                             ta1_z_zzzz_xyy_0, \
                             ta1_z_zzzz_xyz_0, \
                             ta1_z_zzzz_xzz_0, \
                             ta1_z_zzzz_yyy_0, \
                             ta1_z_zzzz_yyz_0, \
                             ta1_z_zzzz_yzz_0, \
                             ta1_z_zzzz_zzz_0, \
                             ta_zzz_xxx_1,     \
                             ta_zzz_xxy_1,     \
                             ta_zzz_xxz_1,     \
                             ta_zzz_xyy_1,     \
                             ta_zzz_xyz_1,     \
                             ta_zzz_xzz_1,     \
                             ta_zzz_yyy_1,     \
                             ta_zzz_yyz_1,     \
                             ta_zzz_yzz_1,     \
                             ta_zzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_xxx_0[i] = 3.0 * ta1_z_zz_xxx_0[i] * fe_0 - 3.0 * ta1_z_zz_xxx_1[i] * fe_0 + ta_zzz_xxx_1[i] + ta1_z_zzz_xxx_0[i] * pa_z[i] -
                              ta1_z_zzz_xxx_1[i] * pc_z[i];

        ta1_z_zzzz_xxy_0[i] = 3.0 * ta1_z_zz_xxy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxy_1[i] * fe_0 + ta_zzz_xxy_1[i] + ta1_z_zzz_xxy_0[i] * pa_z[i] -
                              ta1_z_zzz_xxy_1[i] * pc_z[i];

        ta1_z_zzzz_xxz_0[i] = 3.0 * ta1_z_zz_xxz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxz_1[i] * fe_0 + ta1_z_zzz_xx_0[i] * fe_0 - ta1_z_zzz_xx_1[i] * fe_0 +
                              ta_zzz_xxz_1[i] + ta1_z_zzz_xxz_0[i] * pa_z[i] - ta1_z_zzz_xxz_1[i] * pc_z[i];

        ta1_z_zzzz_xyy_0[i] = 3.0 * ta1_z_zz_xyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xyy_1[i] * fe_0 + ta_zzz_xyy_1[i] + ta1_z_zzz_xyy_0[i] * pa_z[i] -
                              ta1_z_zzz_xyy_1[i] * pc_z[i];

        ta1_z_zzzz_xyz_0[i] = 3.0 * ta1_z_zz_xyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyz_1[i] * fe_0 + ta1_z_zzz_xy_0[i] * fe_0 - ta1_z_zzz_xy_1[i] * fe_0 +
                              ta_zzz_xyz_1[i] + ta1_z_zzz_xyz_0[i] * pa_z[i] - ta1_z_zzz_xyz_1[i] * pc_z[i];

        ta1_z_zzzz_xzz_0[i] = 3.0 * ta1_z_zz_xzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xz_0[i] * fe_0 -
                              2.0 * ta1_z_zzz_xz_1[i] * fe_0 + ta_zzz_xzz_1[i] + ta1_z_zzz_xzz_0[i] * pa_z[i] - ta1_z_zzz_xzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyy_0[i] = 3.0 * ta1_z_zz_yyy_0[i] * fe_0 - 3.0 * ta1_z_zz_yyy_1[i] * fe_0 + ta_zzz_yyy_1[i] + ta1_z_zzz_yyy_0[i] * pa_z[i] -
                              ta1_z_zzz_yyy_1[i] * pc_z[i];

        ta1_z_zzzz_yyz_0[i] = 3.0 * ta1_z_zz_yyz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyz_1[i] * fe_0 + ta1_z_zzz_yy_0[i] * fe_0 - ta1_z_zzz_yy_1[i] * fe_0 +
                              ta_zzz_yyz_1[i] + ta1_z_zzz_yyz_0[i] * pa_z[i] - ta1_z_zzz_yyz_1[i] * pc_z[i];

        ta1_z_zzzz_yzz_0[i] = 3.0 * ta1_z_zz_yzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yz_0[i] * fe_0 -
                              2.0 * ta1_z_zzz_yz_1[i] * fe_0 + ta_zzz_yzz_1[i] + ta1_z_zzz_yzz_0[i] * pa_z[i] - ta1_z_zzz_yzz_1[i] * pc_z[i];

        ta1_z_zzzz_zzz_0[i] = 3.0 * ta1_z_zz_zzz_0[i] * fe_0 - 3.0 * ta1_z_zz_zzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_zz_0[i] * fe_0 -
                              3.0 * ta1_z_zzz_zz_1[i] * fe_0 + ta_zzz_zzz_1[i] + ta1_z_zzz_zzz_0[i] * pa_z[i] - ta1_z_zzz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
