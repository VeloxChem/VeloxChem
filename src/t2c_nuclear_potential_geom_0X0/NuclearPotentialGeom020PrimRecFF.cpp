#include "NuclearPotentialGeom020PrimRecFF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_ff(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_ff,
                                        const size_t idx_npot_geom_020_0_pf,
                                        const size_t idx_npot_geom_020_1_pf,
                                        const size_t idx_npot_geom_020_0_dd,
                                        const size_t idx_npot_geom_020_1_dd,
                                        const size_t idx_npot_geom_010_1_df,
                                        const size_t idx_npot_geom_020_0_df,
                                        const size_t idx_npot_geom_020_1_df,
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

    // Set up components of auxiliary buffer : PF

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

    // Set up components of auxiliary buffer : PF

    auto ta2_xx_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf);

    auto ta2_xx_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 1);

    auto ta2_xx_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 2);

    auto ta2_xx_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 3);

    auto ta2_xx_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 4);

    auto ta2_xx_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 5);

    auto ta2_xx_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 6);

    auto ta2_xx_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 7);

    auto ta2_xx_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 8);

    auto ta2_xx_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 9);

    auto ta2_xx_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 10);

    auto ta2_xx_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 11);

    auto ta2_xx_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 12);

    auto ta2_xx_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 13);

    auto ta2_xx_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 14);

    auto ta2_xx_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 15);

    auto ta2_xx_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 16);

    auto ta2_xx_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 17);

    auto ta2_xx_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 18);

    auto ta2_xx_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 19);

    auto ta2_xx_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 20);

    auto ta2_xx_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 21);

    auto ta2_xx_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 22);

    auto ta2_xx_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 23);

    auto ta2_xx_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 24);

    auto ta2_xx_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 25);

    auto ta2_xx_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 26);

    auto ta2_xx_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 27);

    auto ta2_xx_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 28);

    auto ta2_xx_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 29);

    auto ta2_xy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 30);

    auto ta2_xy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 31);

    auto ta2_xy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 32);

    auto ta2_xy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 33);

    auto ta2_xy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 34);

    auto ta2_xy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 35);

    auto ta2_xy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 36);

    auto ta2_xy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 37);

    auto ta2_xy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 38);

    auto ta2_xy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 39);

    auto ta2_xy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 40);

    auto ta2_xy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 41);

    auto ta2_xy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 42);

    auto ta2_xy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 43);

    auto ta2_xy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 44);

    auto ta2_xy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 45);

    auto ta2_xy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 46);

    auto ta2_xy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 47);

    auto ta2_xy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 48);

    auto ta2_xy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 49);

    auto ta2_xy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 50);

    auto ta2_xy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 51);

    auto ta2_xy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 52);

    auto ta2_xy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 53);

    auto ta2_xy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 54);

    auto ta2_xy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 55);

    auto ta2_xy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 56);

    auto ta2_xy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 57);

    auto ta2_xy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 58);

    auto ta2_xy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 59);

    auto ta2_xz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 60);

    auto ta2_xz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 61);

    auto ta2_xz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 62);

    auto ta2_xz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 63);

    auto ta2_xz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 64);

    auto ta2_xz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 65);

    auto ta2_xz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 66);

    auto ta2_xz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 67);

    auto ta2_xz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 68);

    auto ta2_xz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 69);

    auto ta2_xz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 70);

    auto ta2_xz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 71);

    auto ta2_xz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 72);

    auto ta2_xz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 73);

    auto ta2_xz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 74);

    auto ta2_xz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 75);

    auto ta2_xz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 76);

    auto ta2_xz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 77);

    auto ta2_xz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 78);

    auto ta2_xz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 79);

    auto ta2_xz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 80);

    auto ta2_xz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 81);

    auto ta2_xz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 82);

    auto ta2_xz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 83);

    auto ta2_xz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 84);

    auto ta2_xz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 85);

    auto ta2_xz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 86);

    auto ta2_xz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 87);

    auto ta2_xz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 88);

    auto ta2_xz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 89);

    auto ta2_yy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 90);

    auto ta2_yy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 91);

    auto ta2_yy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 92);

    auto ta2_yy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 93);

    auto ta2_yy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 94);

    auto ta2_yy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 95);

    auto ta2_yy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 96);

    auto ta2_yy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 97);

    auto ta2_yy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 98);

    auto ta2_yy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 99);

    auto ta2_yy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 100);

    auto ta2_yy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 101);

    auto ta2_yy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 102);

    auto ta2_yy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 103);

    auto ta2_yy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 104);

    auto ta2_yy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 105);

    auto ta2_yy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 106);

    auto ta2_yy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 107);

    auto ta2_yy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 108);

    auto ta2_yy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 109);

    auto ta2_yy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 110);

    auto ta2_yy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 111);

    auto ta2_yy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 112);

    auto ta2_yy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 113);

    auto ta2_yy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 114);

    auto ta2_yy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 115);

    auto ta2_yy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 116);

    auto ta2_yy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 117);

    auto ta2_yy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 118);

    auto ta2_yy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 119);

    auto ta2_yz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 120);

    auto ta2_yz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 121);

    auto ta2_yz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 122);

    auto ta2_yz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 123);

    auto ta2_yz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 124);

    auto ta2_yz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 125);

    auto ta2_yz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 126);

    auto ta2_yz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 127);

    auto ta2_yz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 128);

    auto ta2_yz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 129);

    auto ta2_yz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 130);

    auto ta2_yz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 131);

    auto ta2_yz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 132);

    auto ta2_yz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 133);

    auto ta2_yz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 134);

    auto ta2_yz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 135);

    auto ta2_yz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 136);

    auto ta2_yz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 137);

    auto ta2_yz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 138);

    auto ta2_yz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 139);

    auto ta2_yz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 140);

    auto ta2_yz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 141);

    auto ta2_yz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 142);

    auto ta2_yz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 143);

    auto ta2_yz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 144);

    auto ta2_yz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 145);

    auto ta2_yz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 146);

    auto ta2_yz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 147);

    auto ta2_yz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 148);

    auto ta2_yz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 149);

    auto ta2_zz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 150);

    auto ta2_zz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 151);

    auto ta2_zz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 152);

    auto ta2_zz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 153);

    auto ta2_zz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 154);

    auto ta2_zz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 155);

    auto ta2_zz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 156);

    auto ta2_zz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 157);

    auto ta2_zz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 158);

    auto ta2_zz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 159);

    auto ta2_zz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 160);

    auto ta2_zz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 161);

    auto ta2_zz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 162);

    auto ta2_zz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 163);

    auto ta2_zz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 164);

    auto ta2_zz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 165);

    auto ta2_zz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 166);

    auto ta2_zz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 167);

    auto ta2_zz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 168);

    auto ta2_zz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 169);

    auto ta2_zz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 170);

    auto ta2_zz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 171);

    auto ta2_zz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 172);

    auto ta2_zz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 173);

    auto ta2_zz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 174);

    auto ta2_zz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 175);

    auto ta2_zz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 176);

    auto ta2_zz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 177);

    auto ta2_zz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 178);

    auto ta2_zz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 179);

    // Set up components of auxiliary buffer : DD

    auto ta2_xx_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd);

    auto ta2_xx_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 1);

    auto ta2_xx_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 2);

    auto ta2_xx_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 3);

    auto ta2_xx_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 4);

    auto ta2_xx_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 5);

    auto ta2_xx_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 14);

    auto ta2_xx_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 18);

    auto ta2_xx_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 19);

    auto ta2_xx_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 20);

    auto ta2_xx_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 21);

    auto ta2_xx_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 22);

    auto ta2_xx_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 23);

    auto ta2_xx_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 30);

    auto ta2_xx_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 31);

    auto ta2_xx_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 32);

    auto ta2_xx_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 33);

    auto ta2_xx_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 34);

    auto ta2_xx_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 35);

    auto ta2_xy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 36);

    auto ta2_xy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 37);

    auto ta2_xy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 38);

    auto ta2_xy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 39);

    auto ta2_xy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 40);

    auto ta2_xy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 41);

    auto ta2_xy_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 43);

    auto ta2_xy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 54);

    auto ta2_xy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 55);

    auto ta2_xy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 56);

    auto ta2_xy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 57);

    auto ta2_xy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 58);

    auto ta2_xy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 59);

    auto ta2_xy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 66);

    auto ta2_xy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 67);

    auto ta2_xy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 68);

    auto ta2_xy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 69);

    auto ta2_xy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 70);

    auto ta2_xy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 71);

    auto ta2_xz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 72);

    auto ta2_xz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 73);

    auto ta2_xz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 74);

    auto ta2_xz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 75);

    auto ta2_xz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 76);

    auto ta2_xz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 77);

    auto ta2_xz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 86);

    auto ta2_xz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 90);

    auto ta2_xz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 91);

    auto ta2_xz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 92);

    auto ta2_xz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 93);

    auto ta2_xz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 94);

    auto ta2_xz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 95);

    auto ta2_xz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 102);

    auto ta2_xz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 103);

    auto ta2_xz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 104);

    auto ta2_xz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 105);

    auto ta2_xz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 106);

    auto ta2_xz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 107);

    auto ta2_yy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 108);

    auto ta2_yy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 109);

    auto ta2_yy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 110);

    auto ta2_yy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 111);

    auto ta2_yy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 112);

    auto ta2_yy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 113);

    auto ta2_yy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 126);

    auto ta2_yy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 127);

    auto ta2_yy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 128);

    auto ta2_yy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 129);

    auto ta2_yy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 130);

    auto ta2_yy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 131);

    auto ta2_yy_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 136);

    auto ta2_yy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 138);

    auto ta2_yy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 139);

    auto ta2_yy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 140);

    auto ta2_yy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 141);

    auto ta2_yy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 142);

    auto ta2_yy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 143);

    auto ta2_yz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 144);

    auto ta2_yz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 145);

    auto ta2_yz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 146);

    auto ta2_yz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 147);

    auto ta2_yz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 148);

    auto ta2_yz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 149);

    auto ta2_yz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 162);

    auto ta2_yz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 163);

    auto ta2_yz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 164);

    auto ta2_yz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 165);

    auto ta2_yz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 166);

    auto ta2_yz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 167);

    auto ta2_yz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 172);

    auto ta2_yz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 174);

    auto ta2_yz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 175);

    auto ta2_yz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 176);

    auto ta2_yz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 177);

    auto ta2_yz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 178);

    auto ta2_yz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 179);

    auto ta2_zz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 180);

    auto ta2_zz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 181);

    auto ta2_zz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 182);

    auto ta2_zz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 183);

    auto ta2_zz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 184);

    auto ta2_zz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 185);

    auto ta2_zz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 198);

    auto ta2_zz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 199);

    auto ta2_zz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 200);

    auto ta2_zz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 201);

    auto ta2_zz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 202);

    auto ta2_zz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 203);

    auto ta2_zz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 208);

    auto ta2_zz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 210);

    auto ta2_zz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 211);

    auto ta2_zz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 212);

    auto ta2_zz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 213);

    auto ta2_zz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 214);

    auto ta2_zz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 215);

    // Set up components of auxiliary buffer : DD

    auto ta2_xx_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd);

    auto ta2_xx_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 1);

    auto ta2_xx_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 2);

    auto ta2_xx_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 3);

    auto ta2_xx_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 4);

    auto ta2_xx_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 5);

    auto ta2_xx_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 14);

    auto ta2_xx_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 18);

    auto ta2_xx_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 19);

    auto ta2_xx_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 20);

    auto ta2_xx_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 21);

    auto ta2_xx_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 22);

    auto ta2_xx_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 23);

    auto ta2_xx_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 30);

    auto ta2_xx_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 31);

    auto ta2_xx_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 32);

    auto ta2_xx_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 33);

    auto ta2_xx_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 34);

    auto ta2_xx_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 35);

    auto ta2_xy_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 36);

    auto ta2_xy_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 37);

    auto ta2_xy_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 38);

    auto ta2_xy_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 39);

    auto ta2_xy_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 40);

    auto ta2_xy_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 41);

    auto ta2_xy_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 43);

    auto ta2_xy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 54);

    auto ta2_xy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 55);

    auto ta2_xy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 56);

    auto ta2_xy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 57);

    auto ta2_xy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 58);

    auto ta2_xy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 59);

    auto ta2_xy_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 66);

    auto ta2_xy_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 67);

    auto ta2_xy_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 68);

    auto ta2_xy_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 69);

    auto ta2_xy_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 70);

    auto ta2_xy_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 71);

    auto ta2_xz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 72);

    auto ta2_xz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 73);

    auto ta2_xz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 74);

    auto ta2_xz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 75);

    auto ta2_xz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 76);

    auto ta2_xz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 77);

    auto ta2_xz_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 86);

    auto ta2_xz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 90);

    auto ta2_xz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 91);

    auto ta2_xz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 92);

    auto ta2_xz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 93);

    auto ta2_xz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 94);

    auto ta2_xz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 95);

    auto ta2_xz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 102);

    auto ta2_xz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 103);

    auto ta2_xz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 104);

    auto ta2_xz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 105);

    auto ta2_xz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 106);

    auto ta2_xz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 107);

    auto ta2_yy_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 108);

    auto ta2_yy_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 109);

    auto ta2_yy_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 110);

    auto ta2_yy_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 111);

    auto ta2_yy_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 112);

    auto ta2_yy_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 113);

    auto ta2_yy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 126);

    auto ta2_yy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 127);

    auto ta2_yy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 128);

    auto ta2_yy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 129);

    auto ta2_yy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 130);

    auto ta2_yy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 131);

    auto ta2_yy_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 136);

    auto ta2_yy_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 138);

    auto ta2_yy_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 139);

    auto ta2_yy_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 140);

    auto ta2_yy_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 141);

    auto ta2_yy_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 142);

    auto ta2_yy_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 143);

    auto ta2_yz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 144);

    auto ta2_yz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 145);

    auto ta2_yz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 146);

    auto ta2_yz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 147);

    auto ta2_yz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 148);

    auto ta2_yz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 149);

    auto ta2_yz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 162);

    auto ta2_yz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 163);

    auto ta2_yz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 164);

    auto ta2_yz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 165);

    auto ta2_yz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 166);

    auto ta2_yz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 167);

    auto ta2_yz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 172);

    auto ta2_yz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 174);

    auto ta2_yz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 175);

    auto ta2_yz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 176);

    auto ta2_yz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 177);

    auto ta2_yz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 178);

    auto ta2_yz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 179);

    auto ta2_zz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 180);

    auto ta2_zz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 181);

    auto ta2_zz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 182);

    auto ta2_zz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 183);

    auto ta2_zz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 184);

    auto ta2_zz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 185);

    auto ta2_zz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 198);

    auto ta2_zz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 199);

    auto ta2_zz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 200);

    auto ta2_zz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 201);

    auto ta2_zz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 202);

    auto ta2_zz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 203);

    auto ta2_zz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 208);

    auto ta2_zz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 210);

    auto ta2_zz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 211);

    auto ta2_zz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 212);

    auto ta2_zz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 213);

    auto ta2_zz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 214);

    auto ta2_zz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 215);

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

    auto ta1_x_xy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 11);

    auto ta1_x_xy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 13);

    auto ta1_x_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 22);

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

    auto ta1_x_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 47);

    auto ta1_x_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 48);

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

    auto ta1_y_xy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 71);

    auto ta1_y_xy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 73);

    auto ta1_y_xy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 76);

    auto ta1_y_xy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 77);

    auto ta1_y_xy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 78);

    auto ta1_y_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 82);

    auto ta1_y_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 85);

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

    // Set up components of auxiliary buffer : DF

    auto ta2_xx_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df);

    auto ta2_xx_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 1);

    auto ta2_xx_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 2);

    auto ta2_xx_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 3);

    auto ta2_xx_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 4);

    auto ta2_xx_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 5);

    auto ta2_xx_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 6);

    auto ta2_xx_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 7);

    auto ta2_xx_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 8);

    auto ta2_xx_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 9);

    auto ta2_xx_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 10);

    auto ta2_xx_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 11);

    auto ta2_xx_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 12);

    auto ta2_xx_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 13);

    auto ta2_xx_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 15);

    auto ta2_xx_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 16);

    auto ta2_xx_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 20);

    auto ta2_xx_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 21);

    auto ta2_xx_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 22);

    auto ta2_xx_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 23);

    auto ta2_xx_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 24);

    auto ta2_xx_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 25);

    auto ta2_xx_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 29);

    auto ta2_xx_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 30);

    auto ta2_xx_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 31);

    auto ta2_xx_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 32);

    auto ta2_xx_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 33);

    auto ta2_xx_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 34);

    auto ta2_xx_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 35);

    auto ta2_xx_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 36);

    auto ta2_xx_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 37);

    auto ta2_xx_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 38);

    auto ta2_xx_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 39);

    auto ta2_xx_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 42);

    auto ta2_xx_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 45);

    auto ta2_xx_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 47);

    auto ta2_xx_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 48);

    auto ta2_xx_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 49);

    auto ta2_xx_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 50);

    auto ta2_xx_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 51);

    auto ta2_xx_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 52);

    auto ta2_xx_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 53);

    auto ta2_xx_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 54);

    auto ta2_xx_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 55);

    auto ta2_xx_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 56);

    auto ta2_xx_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 57);

    auto ta2_xx_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 58);

    auto ta2_xx_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 59);

    auto ta2_xy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 60);

    auto ta2_xy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 61);

    auto ta2_xy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 62);

    auto ta2_xy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 63);

    auto ta2_xy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 64);

    auto ta2_xy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 65);

    auto ta2_xy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 66);

    auto ta2_xy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 67);

    auto ta2_xy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 68);

    auto ta2_xy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 69);

    auto ta2_xy_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 70);

    auto ta2_xy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 71);

    auto ta2_xy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 73);

    auto ta2_xy_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 74);

    auto ta2_xy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 76);

    auto ta2_xy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 77);

    auto ta2_xy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 78);

    auto ta2_xy_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 80);

    auto ta2_xy_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 81);

    auto ta2_xy_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 82);

    auto ta2_xy_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 83);

    auto ta2_xy_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 85);

    auto ta2_xy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 90);

    auto ta2_xy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 91);

    auto ta2_xy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 92);

    auto ta2_xy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 93);

    auto ta2_xy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 94);

    auto ta2_xy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 95);

    auto ta2_xy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 96);

    auto ta2_xy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 97);

    auto ta2_xy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 98);

    auto ta2_xy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 99);

    auto ta2_xy_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 101);

    auto ta2_xy_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 103);

    auto ta2_xy_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 106);

    auto ta2_xy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 107);

    auto ta2_xy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 108);

    auto ta2_xy_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 109);

    auto ta2_xy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 110);

    auto ta2_xy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 111);

    auto ta2_xy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 112);

    auto ta2_xy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 113);

    auto ta2_xy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 114);

    auto ta2_xy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 115);

    auto ta2_xy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 116);

    auto ta2_xy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 117);

    auto ta2_xy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 118);

    auto ta2_xy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 119);

    auto ta2_xz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 120);

    auto ta2_xz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 121);

    auto ta2_xz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 122);

    auto ta2_xz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 123);

    auto ta2_xz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 124);

    auto ta2_xz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 125);

    auto ta2_xz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 126);

    auto ta2_xz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 127);

    auto ta2_xz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 128);

    auto ta2_xz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 129);

    auto ta2_xz_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 130);

    auto ta2_xz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 131);

    auto ta2_xz_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 132);

    auto ta2_xz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 133);

    auto ta2_xz_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 135);

    auto ta2_xz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 140);

    auto ta2_xz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 142);

    auto ta2_xz_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 144);

    auto ta2_xz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 145);

    auto ta2_xz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 147);

    auto ta2_xz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 148);

    auto ta2_xz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 149);

    auto ta2_xz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 150);

    auto ta2_xz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 151);

    auto ta2_xz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 152);

    auto ta2_xz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 153);

    auto ta2_xz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 154);

    auto ta2_xz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 155);

    auto ta2_xz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 156);

    auto ta2_xz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 157);

    auto ta2_xz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 158);

    auto ta2_xz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 159);

    auto ta2_xz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 162);

    auto ta2_xz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 165);

    auto ta2_xz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 166);

    auto ta2_xz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 167);

    auto ta2_xz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 168);

    auto ta2_xz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 169);

    auto ta2_xz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 170);

    auto ta2_xz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 171);

    auto ta2_xz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 172);

    auto ta2_xz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 173);

    auto ta2_xz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 174);

    auto ta2_xz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 175);

    auto ta2_xz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 176);

    auto ta2_xz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 177);

    auto ta2_xz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 178);

    auto ta2_xz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 179);

    auto ta2_yy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 180);

    auto ta2_yy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 181);

    auto ta2_yy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 182);

    auto ta2_yy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 183);

    auto ta2_yy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 184);

    auto ta2_yy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 185);

    auto ta2_yy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 186);

    auto ta2_yy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 187);

    auto ta2_yy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 188);

    auto ta2_yy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 189);

    auto ta2_yy_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 190);

    auto ta2_yy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 191);

    auto ta2_yy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 193);

    auto ta2_yy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 196);

    auto ta2_yy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 197);

    auto ta2_yy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 198);

    auto ta2_yy_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 202);

    auto ta2_yy_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 205);

    auto ta2_yy_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 207);

    auto ta2_yy_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 208);

    auto ta2_yy_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 209);

    auto ta2_yy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 210);

    auto ta2_yy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 211);

    auto ta2_yy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 212);

    auto ta2_yy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 213);

    auto ta2_yy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 214);

    auto ta2_yy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 215);

    auto ta2_yy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 216);

    auto ta2_yy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 217);

    auto ta2_yy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 218);

    auto ta2_yy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 219);

    auto ta2_yy_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 221);

    auto ta2_yy_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 223);

    auto ta2_yy_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 224);

    auto ta2_yy_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 226);

    auto ta2_yy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 227);

    auto ta2_yy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 228);

    auto ta2_yy_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 229);

    auto ta2_yy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 230);

    auto ta2_yy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 231);

    auto ta2_yy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 232);

    auto ta2_yy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 233);

    auto ta2_yy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 234);

    auto ta2_yy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 235);

    auto ta2_yy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 236);

    auto ta2_yy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 237);

    auto ta2_yy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 238);

    auto ta2_yy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 239);

    auto ta2_yz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 240);

    auto ta2_yz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 241);

    auto ta2_yz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 242);

    auto ta2_yz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 243);

    auto ta2_yz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 244);

    auto ta2_yz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 245);

    auto ta2_yz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 246);

    auto ta2_yz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 247);

    auto ta2_yz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 248);

    auto ta2_yz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 249);

    auto ta2_yz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 251);

    auto ta2_yz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 253);

    auto ta2_yz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 256);

    auto ta2_yz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 257);

    auto ta2_yz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 258);

    auto ta2_yz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 260);

    auto ta2_yz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 262);

    auto ta2_yz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 265);

    auto ta2_yz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 267);

    auto ta2_yz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 268);

    auto ta2_yz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 269);

    auto ta2_yz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 270);

    auto ta2_yz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 271);

    auto ta2_yz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 272);

    auto ta2_yz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 273);

    auto ta2_yz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 274);

    auto ta2_yz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 275);

    auto ta2_yz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 276);

    auto ta2_yz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 277);

    auto ta2_yz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 278);

    auto ta2_yz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 279);

    auto ta2_yz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 282);

    auto ta2_yz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 284);

    auto ta2_yz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 285);

    auto ta2_yz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 286);

    auto ta2_yz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 287);

    auto ta2_yz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 288);

    auto ta2_yz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 289);

    auto ta2_yz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 290);

    auto ta2_yz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 291);

    auto ta2_yz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 292);

    auto ta2_yz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 293);

    auto ta2_yz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 294);

    auto ta2_yz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 295);

    auto ta2_yz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 296);

    auto ta2_yz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 297);

    auto ta2_yz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 298);

    auto ta2_yz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 299);

    auto ta2_zz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 300);

    auto ta2_zz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 301);

    auto ta2_zz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 302);

    auto ta2_zz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 303);

    auto ta2_zz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 304);

    auto ta2_zz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 305);

    auto ta2_zz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 306);

    auto ta2_zz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 307);

    auto ta2_zz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 308);

    auto ta2_zz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 309);

    auto ta2_zz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 311);

    auto ta2_zz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 313);

    auto ta2_zz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 316);

    auto ta2_zz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 317);

    auto ta2_zz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 318);

    auto ta2_zz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 320);

    auto ta2_zz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 322);

    auto ta2_zz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 325);

    auto ta2_zz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 327);

    auto ta2_zz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 328);

    auto ta2_zz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 329);

    auto ta2_zz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 330);

    auto ta2_zz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 331);

    auto ta2_zz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 332);

    auto ta2_zz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 333);

    auto ta2_zz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 334);

    auto ta2_zz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 335);

    auto ta2_zz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 336);

    auto ta2_zz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 337);

    auto ta2_zz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 338);

    auto ta2_zz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 339);

    auto ta2_zz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 342);

    auto ta2_zz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 344);

    auto ta2_zz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 345);

    auto ta2_zz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 346);

    auto ta2_zz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 347);

    auto ta2_zz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 348);

    auto ta2_zz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 349);

    auto ta2_zz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 350);

    auto ta2_zz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 351);

    auto ta2_zz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 352);

    auto ta2_zz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 353);

    auto ta2_zz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 354);

    auto ta2_zz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 355);

    auto ta2_zz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 356);

    auto ta2_zz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 357);

    auto ta2_zz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 358);

    auto ta2_zz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 359);

    // Set up components of auxiliary buffer : DF

    auto ta2_xx_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df);

    auto ta2_xx_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 1);

    auto ta2_xx_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 2);

    auto ta2_xx_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 3);

    auto ta2_xx_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 4);

    auto ta2_xx_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 5);

    auto ta2_xx_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 6);

    auto ta2_xx_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 7);

    auto ta2_xx_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 8);

    auto ta2_xx_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 9);

    auto ta2_xx_xy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 10);

    auto ta2_xx_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 11);

    auto ta2_xx_xy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 12);

    auto ta2_xx_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 13);

    auto ta2_xx_xy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 15);

    auto ta2_xx_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 16);

    auto ta2_xx_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 20);

    auto ta2_xx_xz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 21);

    auto ta2_xx_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 22);

    auto ta2_xx_xz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 23);

    auto ta2_xx_xz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 24);

    auto ta2_xx_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 25);

    auto ta2_xx_xz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 29);

    auto ta2_xx_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 30);

    auto ta2_xx_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 31);

    auto ta2_xx_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 32);

    auto ta2_xx_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 33);

    auto ta2_xx_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 34);

    auto ta2_xx_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 35);

    auto ta2_xx_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 36);

    auto ta2_xx_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 37);

    auto ta2_xx_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 38);

    auto ta2_xx_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 39);

    auto ta2_xx_yz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 42);

    auto ta2_xx_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 45);

    auto ta2_xx_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 47);

    auto ta2_xx_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 48);

    auto ta2_xx_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 49);

    auto ta2_xx_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 50);

    auto ta2_xx_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 51);

    auto ta2_xx_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 52);

    auto ta2_xx_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 53);

    auto ta2_xx_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 54);

    auto ta2_xx_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 55);

    auto ta2_xx_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 56);

    auto ta2_xx_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 57);

    auto ta2_xx_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 58);

    auto ta2_xx_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 59);

    auto ta2_xy_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 60);

    auto ta2_xy_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 61);

    auto ta2_xy_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 62);

    auto ta2_xy_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 63);

    auto ta2_xy_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 64);

    auto ta2_xy_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 65);

    auto ta2_xy_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 66);

    auto ta2_xy_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 67);

    auto ta2_xy_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 68);

    auto ta2_xy_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 69);

    auto ta2_xy_xy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 70);

    auto ta2_xy_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 71);

    auto ta2_xy_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 73);

    auto ta2_xy_xy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 74);

    auto ta2_xy_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 76);

    auto ta2_xy_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 77);

    auto ta2_xy_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 78);

    auto ta2_xy_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 80);

    auto ta2_xy_xz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 81);

    auto ta2_xy_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 82);

    auto ta2_xy_xz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 83);

    auto ta2_xy_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 85);

    auto ta2_xy_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 90);

    auto ta2_xy_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 91);

    auto ta2_xy_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 92);

    auto ta2_xy_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 93);

    auto ta2_xy_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 94);

    auto ta2_xy_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 95);

    auto ta2_xy_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 96);

    auto ta2_xy_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 97);

    auto ta2_xy_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 98);

    auto ta2_xy_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 99);

    auto ta2_xy_yz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 101);

    auto ta2_xy_yz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 103);

    auto ta2_xy_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 106);

    auto ta2_xy_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 107);

    auto ta2_xy_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 108);

    auto ta2_xy_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 109);

    auto ta2_xy_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 110);

    auto ta2_xy_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 111);

    auto ta2_xy_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 112);

    auto ta2_xy_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 113);

    auto ta2_xy_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 114);

    auto ta2_xy_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 115);

    auto ta2_xy_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 116);

    auto ta2_xy_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 117);

    auto ta2_xy_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 118);

    auto ta2_xy_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 119);

    auto ta2_xz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 120);

    auto ta2_xz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 121);

    auto ta2_xz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 122);

    auto ta2_xz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 123);

    auto ta2_xz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 124);

    auto ta2_xz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 125);

    auto ta2_xz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 126);

    auto ta2_xz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 127);

    auto ta2_xz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 128);

    auto ta2_xz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 129);

    auto ta2_xz_xy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 130);

    auto ta2_xz_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 131);

    auto ta2_xz_xy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 132);

    auto ta2_xz_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 133);

    auto ta2_xz_xy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 135);

    auto ta2_xz_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 140);

    auto ta2_xz_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 142);

    auto ta2_xz_xz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 144);

    auto ta2_xz_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 145);

    auto ta2_xz_xz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 147);

    auto ta2_xz_xz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 148);

    auto ta2_xz_xz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 149);

    auto ta2_xz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 150);

    auto ta2_xz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 151);

    auto ta2_xz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 152);

    auto ta2_xz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 153);

    auto ta2_xz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 154);

    auto ta2_xz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 155);

    auto ta2_xz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 156);

    auto ta2_xz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 157);

    auto ta2_xz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 158);

    auto ta2_xz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 159);

    auto ta2_xz_yz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 162);

    auto ta2_xz_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 165);

    auto ta2_xz_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 166);

    auto ta2_xz_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 167);

    auto ta2_xz_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 168);

    auto ta2_xz_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 169);

    auto ta2_xz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 170);

    auto ta2_xz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 171);

    auto ta2_xz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 172);

    auto ta2_xz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 173);

    auto ta2_xz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 174);

    auto ta2_xz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 175);

    auto ta2_xz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 176);

    auto ta2_xz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 177);

    auto ta2_xz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 178);

    auto ta2_xz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 179);

    auto ta2_yy_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 180);

    auto ta2_yy_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 181);

    auto ta2_yy_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 182);

    auto ta2_yy_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 183);

    auto ta2_yy_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 184);

    auto ta2_yy_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 185);

    auto ta2_yy_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 186);

    auto ta2_yy_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 187);

    auto ta2_yy_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 188);

    auto ta2_yy_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 189);

    auto ta2_yy_xy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 190);

    auto ta2_yy_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 191);

    auto ta2_yy_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 193);

    auto ta2_yy_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 196);

    auto ta2_yy_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 197);

    auto ta2_yy_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 198);

    auto ta2_yy_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 202);

    auto ta2_yy_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 205);

    auto ta2_yy_xz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 207);

    auto ta2_yy_xz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 208);

    auto ta2_yy_xz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 209);

    auto ta2_yy_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 210);

    auto ta2_yy_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 211);

    auto ta2_yy_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 212);

    auto ta2_yy_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 213);

    auto ta2_yy_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 214);

    auto ta2_yy_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 215);

    auto ta2_yy_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 216);

    auto ta2_yy_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 217);

    auto ta2_yy_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 218);

    auto ta2_yy_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 219);

    auto ta2_yy_yz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 221);

    auto ta2_yy_yz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 223);

    auto ta2_yy_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 224);

    auto ta2_yy_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 226);

    auto ta2_yy_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 227);

    auto ta2_yy_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 228);

    auto ta2_yy_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 229);

    auto ta2_yy_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 230);

    auto ta2_yy_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 231);

    auto ta2_yy_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 232);

    auto ta2_yy_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 233);

    auto ta2_yy_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 234);

    auto ta2_yy_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 235);

    auto ta2_yy_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 236);

    auto ta2_yy_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 237);

    auto ta2_yy_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 238);

    auto ta2_yy_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 239);

    auto ta2_yz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 240);

    auto ta2_yz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 241);

    auto ta2_yz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 242);

    auto ta2_yz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 243);

    auto ta2_yz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 244);

    auto ta2_yz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 245);

    auto ta2_yz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 246);

    auto ta2_yz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 247);

    auto ta2_yz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 248);

    auto ta2_yz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 249);

    auto ta2_yz_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 251);

    auto ta2_yz_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 253);

    auto ta2_yz_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 256);

    auto ta2_yz_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 257);

    auto ta2_yz_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 258);

    auto ta2_yz_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 260);

    auto ta2_yz_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 262);

    auto ta2_yz_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 265);

    auto ta2_yz_xz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 267);

    auto ta2_yz_xz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 268);

    auto ta2_yz_xz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 269);

    auto ta2_yz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 270);

    auto ta2_yz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 271);

    auto ta2_yz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 272);

    auto ta2_yz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 273);

    auto ta2_yz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 274);

    auto ta2_yz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 275);

    auto ta2_yz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 276);

    auto ta2_yz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 277);

    auto ta2_yz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 278);

    auto ta2_yz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 279);

    auto ta2_yz_yz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 282);

    auto ta2_yz_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 284);

    auto ta2_yz_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 285);

    auto ta2_yz_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 286);

    auto ta2_yz_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 287);

    auto ta2_yz_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 288);

    auto ta2_yz_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 289);

    auto ta2_yz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 290);

    auto ta2_yz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 291);

    auto ta2_yz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 292);

    auto ta2_yz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 293);

    auto ta2_yz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 294);

    auto ta2_yz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 295);

    auto ta2_yz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 296);

    auto ta2_yz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 297);

    auto ta2_yz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 298);

    auto ta2_yz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 299);

    auto ta2_zz_xx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 300);

    auto ta2_zz_xx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 301);

    auto ta2_zz_xx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 302);

    auto ta2_zz_xx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 303);

    auto ta2_zz_xx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 304);

    auto ta2_zz_xx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 305);

    auto ta2_zz_xx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 306);

    auto ta2_zz_xx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 307);

    auto ta2_zz_xx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 308);

    auto ta2_zz_xx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 309);

    auto ta2_zz_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 311);

    auto ta2_zz_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 313);

    auto ta2_zz_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 316);

    auto ta2_zz_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 317);

    auto ta2_zz_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 318);

    auto ta2_zz_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 320);

    auto ta2_zz_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 322);

    auto ta2_zz_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 325);

    auto ta2_zz_xz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 327);

    auto ta2_zz_xz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 328);

    auto ta2_zz_xz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 329);

    auto ta2_zz_yy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 330);

    auto ta2_zz_yy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 331);

    auto ta2_zz_yy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 332);

    auto ta2_zz_yy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 333);

    auto ta2_zz_yy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 334);

    auto ta2_zz_yy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 335);

    auto ta2_zz_yy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 336);

    auto ta2_zz_yy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 337);

    auto ta2_zz_yy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 338);

    auto ta2_zz_yy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 339);

    auto ta2_zz_yz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 342);

    auto ta2_zz_yz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 344);

    auto ta2_zz_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 345);

    auto ta2_zz_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 346);

    auto ta2_zz_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 347);

    auto ta2_zz_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 348);

    auto ta2_zz_yz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 349);

    auto ta2_zz_zz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 350);

    auto ta2_zz_zz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 351);

    auto ta2_zz_zz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 352);

    auto ta2_zz_zz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 353);

    auto ta2_zz_zz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 354);

    auto ta2_zz_zz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 355);

    auto ta2_zz_zz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 356);

    auto ta2_zz_zz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 357);

    auto ta2_zz_zz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 358);

    auto ta2_zz_zz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 359);

    // Set up 0-10 components of targeted buffer : FF

    auto ta2_xx_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff);

    auto ta2_xx_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 1);

    auto ta2_xx_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 2);

    auto ta2_xx_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 3);

    auto ta2_xx_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 4);

    auto ta2_xx_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 5);

    auto ta2_xx_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 6);

    auto ta2_xx_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 7);

    auto ta2_xx_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 8);

    auto ta2_xx_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 9);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xx_xxx_1, ta1_x_xx_xxy_1, ta1_x_xx_xxz_1, ta1_x_xx_xyy_1, ta1_x_xx_xyz_1, ta1_x_xx_xzz_1, ta1_x_xx_yyy_1, ta1_x_xx_yyz_1, ta1_x_xx_yzz_1, ta1_x_xx_zzz_1, ta2_xx_x_xxx_0, ta2_xx_x_xxx_1, ta2_xx_x_xxy_0, ta2_xx_x_xxy_1, ta2_xx_x_xxz_0, ta2_xx_x_xxz_1, ta2_xx_x_xyy_0, ta2_xx_x_xyy_1, ta2_xx_x_xyz_0, ta2_xx_x_xyz_1, ta2_xx_x_xzz_0, ta2_xx_x_xzz_1, ta2_xx_x_yyy_0, ta2_xx_x_yyy_1, ta2_xx_x_yyz_0, ta2_xx_x_yyz_1, ta2_xx_x_yzz_0, ta2_xx_x_yzz_1, ta2_xx_x_zzz_0, ta2_xx_x_zzz_1, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xxx_xxx_0, ta2_xx_xxx_xxy_0, ta2_xx_xxx_xxz_0, ta2_xx_xxx_xyy_0, ta2_xx_xxx_xyz_0, ta2_xx_xxx_xzz_0, ta2_xx_xxx_yyy_0, ta2_xx_xxx_yyz_0, ta2_xx_xxx_yzz_0, ta2_xx_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxx_xxx_0[i] = 2.0 * ta2_xx_x_xxx_0[i] * fe_0 - 2.0 * ta2_xx_x_xxx_1[i] * fe_0 + 3.0 * ta2_xx_xx_xx_0[i] * fe_0 - 3.0 * ta2_xx_xx_xx_1[i] * fe_0 + 2.0 * ta1_x_xx_xxx_1[i] + ta2_xx_xx_xxx_0[i] * pa_x[i] - ta2_xx_xx_xxx_1[i] * pc_x[i];

        ta2_xx_xxx_xxy_0[i] = 2.0 * ta2_xx_x_xxy_0[i] * fe_0 - 2.0 * ta2_xx_x_xxy_1[i] * fe_0 + 2.0 * ta2_xx_xx_xy_0[i] * fe_0 - 2.0 * ta2_xx_xx_xy_1[i] * fe_0 + 2.0 * ta1_x_xx_xxy_1[i] + ta2_xx_xx_xxy_0[i] * pa_x[i] - ta2_xx_xx_xxy_1[i] * pc_x[i];

        ta2_xx_xxx_xxz_0[i] = 2.0 * ta2_xx_x_xxz_0[i] * fe_0 - 2.0 * ta2_xx_x_xxz_1[i] * fe_0 + 2.0 * ta2_xx_xx_xz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xz_1[i] * fe_0 + 2.0 * ta1_x_xx_xxz_1[i] + ta2_xx_xx_xxz_0[i] * pa_x[i] - ta2_xx_xx_xxz_1[i] * pc_x[i];

        ta2_xx_xxx_xyy_0[i] = 2.0 * ta2_xx_x_xyy_0[i] * fe_0 - 2.0 * ta2_xx_x_xyy_1[i] * fe_0 + ta2_xx_xx_yy_0[i] * fe_0 - ta2_xx_xx_yy_1[i] * fe_0 + 2.0 * ta1_x_xx_xyy_1[i] + ta2_xx_xx_xyy_0[i] * pa_x[i] - ta2_xx_xx_xyy_1[i] * pc_x[i];

        ta2_xx_xxx_xyz_0[i] = 2.0 * ta2_xx_x_xyz_0[i] * fe_0 - 2.0 * ta2_xx_x_xyz_1[i] * fe_0 + ta2_xx_xx_yz_0[i] * fe_0 - ta2_xx_xx_yz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyz_1[i] + ta2_xx_xx_xyz_0[i] * pa_x[i] - ta2_xx_xx_xyz_1[i] * pc_x[i];

        ta2_xx_xxx_xzz_0[i] = 2.0 * ta2_xx_x_xzz_0[i] * fe_0 - 2.0 * ta2_xx_x_xzz_1[i] * fe_0 + ta2_xx_xx_zz_0[i] * fe_0 - ta2_xx_xx_zz_1[i] * fe_0 + 2.0 * ta1_x_xx_xzz_1[i] + ta2_xx_xx_xzz_0[i] * pa_x[i] - ta2_xx_xx_xzz_1[i] * pc_x[i];

        ta2_xx_xxx_yyy_0[i] = 2.0 * ta2_xx_x_yyy_0[i] * fe_0 - 2.0 * ta2_xx_x_yyy_1[i] * fe_0 + 2.0 * ta1_x_xx_yyy_1[i] + ta2_xx_xx_yyy_0[i] * pa_x[i] - ta2_xx_xx_yyy_1[i] * pc_x[i];

        ta2_xx_xxx_yyz_0[i] = 2.0 * ta2_xx_x_yyz_0[i] * fe_0 - 2.0 * ta2_xx_x_yyz_1[i] * fe_0 + 2.0 * ta1_x_xx_yyz_1[i] + ta2_xx_xx_yyz_0[i] * pa_x[i] - ta2_xx_xx_yyz_1[i] * pc_x[i];

        ta2_xx_xxx_yzz_0[i] = 2.0 * ta2_xx_x_yzz_0[i] * fe_0 - 2.0 * ta2_xx_x_yzz_1[i] * fe_0 + 2.0 * ta1_x_xx_yzz_1[i] + ta2_xx_xx_yzz_0[i] * pa_x[i] - ta2_xx_xx_yzz_1[i] * pc_x[i];

        ta2_xx_xxx_zzz_0[i] = 2.0 * ta2_xx_x_zzz_0[i] * fe_0 - 2.0 * ta2_xx_x_zzz_1[i] * fe_0 + 2.0 * ta1_x_xx_zzz_1[i] + ta2_xx_xx_zzz_0[i] * pa_x[i] - ta2_xx_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto ta2_xx_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 10);

    auto ta2_xx_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 11);

    auto ta2_xx_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 12);

    auto ta2_xx_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 13);

    auto ta2_xx_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 14);

    auto ta2_xx_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 15);

    auto ta2_xx_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 16);

    auto ta2_xx_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 17);

    auto ta2_xx_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 18);

    auto ta2_xx_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 19);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xxy_xxx_0, ta2_xx_xxy_xxy_0, ta2_xx_xxy_xxz_0, ta2_xx_xxy_xyy_0, ta2_xx_xxy_xyz_0, ta2_xx_xxy_xzz_0, ta2_xx_xxy_yyy_0, ta2_xx_xxy_yyz_0, ta2_xx_xxy_yzz_0, ta2_xx_xxy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxy_xxx_0[i] = ta2_xx_xx_xxx_0[i] * pa_y[i] - ta2_xx_xx_xxx_1[i] * pc_y[i];

        ta2_xx_xxy_xxy_0[i] = ta2_xx_xx_xx_0[i] * fe_0 - ta2_xx_xx_xx_1[i] * fe_0 + ta2_xx_xx_xxy_0[i] * pa_y[i] - ta2_xx_xx_xxy_1[i] * pc_y[i];

        ta2_xx_xxy_xxz_0[i] = ta2_xx_xx_xxz_0[i] * pa_y[i] - ta2_xx_xx_xxz_1[i] * pc_y[i];

        ta2_xx_xxy_xyy_0[i] = 2.0 * ta2_xx_xx_xy_0[i] * fe_0 - 2.0 * ta2_xx_xx_xy_1[i] * fe_0 + ta2_xx_xx_xyy_0[i] * pa_y[i] - ta2_xx_xx_xyy_1[i] * pc_y[i];

        ta2_xx_xxy_xyz_0[i] = ta2_xx_xx_xz_0[i] * fe_0 - ta2_xx_xx_xz_1[i] * fe_0 + ta2_xx_xx_xyz_0[i] * pa_y[i] - ta2_xx_xx_xyz_1[i] * pc_y[i];

        ta2_xx_xxy_xzz_0[i] = ta2_xx_xx_xzz_0[i] * pa_y[i] - ta2_xx_xx_xzz_1[i] * pc_y[i];

        ta2_xx_xxy_yyy_0[i] = 3.0 * ta2_xx_xx_yy_0[i] * fe_0 - 3.0 * ta2_xx_xx_yy_1[i] * fe_0 + ta2_xx_xx_yyy_0[i] * pa_y[i] - ta2_xx_xx_yyy_1[i] * pc_y[i];

        ta2_xx_xxy_yyz_0[i] = 2.0 * ta2_xx_xx_yz_0[i] * fe_0 - 2.0 * ta2_xx_xx_yz_1[i] * fe_0 + ta2_xx_xx_yyz_0[i] * pa_y[i] - ta2_xx_xx_yyz_1[i] * pc_y[i];

        ta2_xx_xxy_yzz_0[i] = ta2_xx_xx_zz_0[i] * fe_0 - ta2_xx_xx_zz_1[i] * fe_0 + ta2_xx_xx_yzz_0[i] * pa_y[i] - ta2_xx_xx_yzz_1[i] * pc_y[i];

        ta2_xx_xxy_zzz_0[i] = ta2_xx_xx_zzz_0[i] * pa_y[i] - ta2_xx_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto ta2_xx_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 20);

    auto ta2_xx_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 21);

    auto ta2_xx_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 22);

    auto ta2_xx_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 23);

    auto ta2_xx_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 24);

    auto ta2_xx_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 25);

    auto ta2_xx_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 26);

    auto ta2_xx_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 27);

    auto ta2_xx_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 28);

    auto ta2_xx_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 29);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xxz_xxx_0, ta2_xx_xxz_xxy_0, ta2_xx_xxz_xxz_0, ta2_xx_xxz_xyy_0, ta2_xx_xxz_xyz_0, ta2_xx_xxz_xzz_0, ta2_xx_xxz_yyy_0, ta2_xx_xxz_yyz_0, ta2_xx_xxz_yzz_0, ta2_xx_xxz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxz_xxx_0[i] = ta2_xx_xx_xxx_0[i] * pa_z[i] - ta2_xx_xx_xxx_1[i] * pc_z[i];

        ta2_xx_xxz_xxy_0[i] = ta2_xx_xx_xxy_0[i] * pa_z[i] - ta2_xx_xx_xxy_1[i] * pc_z[i];

        ta2_xx_xxz_xxz_0[i] = ta2_xx_xx_xx_0[i] * fe_0 - ta2_xx_xx_xx_1[i] * fe_0 + ta2_xx_xx_xxz_0[i] * pa_z[i] - ta2_xx_xx_xxz_1[i] * pc_z[i];

        ta2_xx_xxz_xyy_0[i] = ta2_xx_xx_xyy_0[i] * pa_z[i] - ta2_xx_xx_xyy_1[i] * pc_z[i];

        ta2_xx_xxz_xyz_0[i] = ta2_xx_xx_xy_0[i] * fe_0 - ta2_xx_xx_xy_1[i] * fe_0 + ta2_xx_xx_xyz_0[i] * pa_z[i] - ta2_xx_xx_xyz_1[i] * pc_z[i];

        ta2_xx_xxz_xzz_0[i] = 2.0 * ta2_xx_xx_xz_0[i] * fe_0 - 2.0 * ta2_xx_xx_xz_1[i] * fe_0 + ta2_xx_xx_xzz_0[i] * pa_z[i] - ta2_xx_xx_xzz_1[i] * pc_z[i];

        ta2_xx_xxz_yyy_0[i] = ta2_xx_xx_yyy_0[i] * pa_z[i] - ta2_xx_xx_yyy_1[i] * pc_z[i];

        ta2_xx_xxz_yyz_0[i] = ta2_xx_xx_yy_0[i] * fe_0 - ta2_xx_xx_yy_1[i] * fe_0 + ta2_xx_xx_yyz_0[i] * pa_z[i] - ta2_xx_xx_yyz_1[i] * pc_z[i];

        ta2_xx_xxz_yzz_0[i] = 2.0 * ta2_xx_xx_yz_0[i] * fe_0 - 2.0 * ta2_xx_xx_yz_1[i] * fe_0 + ta2_xx_xx_yzz_0[i] * pa_z[i] - ta2_xx_xx_yzz_1[i] * pc_z[i];

        ta2_xx_xxz_zzz_0[i] = 3.0 * ta2_xx_xx_zz_0[i] * fe_0 - 3.0 * ta2_xx_xx_zz_1[i] * fe_0 + ta2_xx_xx_zzz_0[i] * pa_z[i] - ta2_xx_xx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto ta2_xx_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 30);

    auto ta2_xx_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 31);

    auto ta2_xx_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 32);

    auto ta2_xx_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 33);

    auto ta2_xx_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 34);

    auto ta2_xx_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 35);

    auto ta2_xx_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 36);

    auto ta2_xx_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 37);

    auto ta2_xx_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 38);

    auto ta2_xx_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 39);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_yy_xxy_1, ta1_x_yy_xyy_1, ta1_x_yy_xyz_1, ta1_x_yy_yyy_1, ta1_x_yy_yyz_1, ta1_x_yy_yzz_1, ta1_x_yy_zzz_1, ta2_xx_x_xxx_0, ta2_xx_x_xxx_1, ta2_xx_x_xxz_0, ta2_xx_x_xxz_1, ta2_xx_x_xzz_0, ta2_xx_x_xzz_1, ta2_xx_xy_xxx_0, ta2_xx_xy_xxx_1, ta2_xx_xy_xxz_0, ta2_xx_xy_xxz_1, ta2_xx_xy_xzz_0, ta2_xx_xy_xzz_1, ta2_xx_xyy_xxx_0, ta2_xx_xyy_xxy_0, ta2_xx_xyy_xxz_0, ta2_xx_xyy_xyy_0, ta2_xx_xyy_xyz_0, ta2_xx_xyy_xzz_0, ta2_xx_xyy_yyy_0, ta2_xx_xyy_yyz_0, ta2_xx_xyy_yzz_0, ta2_xx_xyy_zzz_0, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_zzz_0, ta2_xx_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyy_xxx_0[i] = ta2_xx_x_xxx_0[i] * fe_0 - ta2_xx_x_xxx_1[i] * fe_0 + ta2_xx_xy_xxx_0[i] * pa_y[i] - ta2_xx_xy_xxx_1[i] * pc_y[i];

        ta2_xx_xyy_xxy_0[i] = 2.0 * ta2_xx_yy_xy_0[i] * fe_0 - 2.0 * ta2_xx_yy_xy_1[i] * fe_0 + 2.0 * ta1_x_yy_xxy_1[i] + ta2_xx_yy_xxy_0[i] * pa_x[i] - ta2_xx_yy_xxy_1[i] * pc_x[i];

        ta2_xx_xyy_xxz_0[i] = ta2_xx_x_xxz_0[i] * fe_0 - ta2_xx_x_xxz_1[i] * fe_0 + ta2_xx_xy_xxz_0[i] * pa_y[i] - ta2_xx_xy_xxz_1[i] * pc_y[i];

        ta2_xx_xyy_xyy_0[i] = ta2_xx_yy_yy_0[i] * fe_0 - ta2_xx_yy_yy_1[i] * fe_0 + 2.0 * ta1_x_yy_xyy_1[i] + ta2_xx_yy_xyy_0[i] * pa_x[i] - ta2_xx_yy_xyy_1[i] * pc_x[i];

        ta2_xx_xyy_xyz_0[i] = ta2_xx_yy_yz_0[i] * fe_0 - ta2_xx_yy_yz_1[i] * fe_0 + 2.0 * ta1_x_yy_xyz_1[i] + ta2_xx_yy_xyz_0[i] * pa_x[i] - ta2_xx_yy_xyz_1[i] * pc_x[i];

        ta2_xx_xyy_xzz_0[i] = ta2_xx_x_xzz_0[i] * fe_0 - ta2_xx_x_xzz_1[i] * fe_0 + ta2_xx_xy_xzz_0[i] * pa_y[i] - ta2_xx_xy_xzz_1[i] * pc_y[i];

        ta2_xx_xyy_yyy_0[i] = 2.0 * ta1_x_yy_yyy_1[i] + ta2_xx_yy_yyy_0[i] * pa_x[i] - ta2_xx_yy_yyy_1[i] * pc_x[i];

        ta2_xx_xyy_yyz_0[i] = 2.0 * ta1_x_yy_yyz_1[i] + ta2_xx_yy_yyz_0[i] * pa_x[i] - ta2_xx_yy_yyz_1[i] * pc_x[i];

        ta2_xx_xyy_yzz_0[i] = 2.0 * ta1_x_yy_yzz_1[i] + ta2_xx_yy_yzz_0[i] * pa_x[i] - ta2_xx_yy_yzz_1[i] * pc_x[i];

        ta2_xx_xyy_zzz_0[i] = 2.0 * ta1_x_yy_zzz_1[i] + ta2_xx_yy_zzz_0[i] * pa_x[i] - ta2_xx_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto ta2_xx_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 40);

    auto ta2_xx_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 41);

    auto ta2_xx_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 42);

    auto ta2_xx_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 43);

    auto ta2_xx_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 44);

    auto ta2_xx_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 45);

    auto ta2_xx_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 46);

    auto ta2_xx_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 47);

    auto ta2_xx_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 48);

    auto ta2_xx_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 49);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_yz_yyz_1, ta1_x_yz_yzz_1, ta2_xx_xy_xxy_0, ta2_xx_xy_xxy_1, ta2_xx_xy_xyy_0, ta2_xx_xy_xyy_1, ta2_xx_xy_yyy_0, ta2_xx_xy_yyy_1, ta2_xx_xyz_xxx_0, ta2_xx_xyz_xxy_0, ta2_xx_xyz_xxz_0, ta2_xx_xyz_xyy_0, ta2_xx_xyz_xyz_0, ta2_xx_xyz_xzz_0, ta2_xx_xyz_yyy_0, ta2_xx_xyz_yyz_0, ta2_xx_xyz_yzz_0, ta2_xx_xyz_zzz_0, ta2_xx_xz_xxx_0, ta2_xx_xz_xxx_1, ta2_xx_xz_xxz_0, ta2_xx_xz_xxz_1, ta2_xx_xz_xyz_0, ta2_xx_xz_xyz_1, ta2_xx_xz_xz_0, ta2_xx_xz_xz_1, ta2_xx_xz_xzz_0, ta2_xx_xz_xzz_1, ta2_xx_xz_zzz_0, ta2_xx_xz_zzz_1, ta2_xx_yz_yyz_0, ta2_xx_yz_yyz_1, ta2_xx_yz_yzz_0, ta2_xx_yz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyz_xxx_0[i] = ta2_xx_xz_xxx_0[i] * pa_y[i] - ta2_xx_xz_xxx_1[i] * pc_y[i];

        ta2_xx_xyz_xxy_0[i] = ta2_xx_xy_xxy_0[i] * pa_z[i] - ta2_xx_xy_xxy_1[i] * pc_z[i];

        ta2_xx_xyz_xxz_0[i] = ta2_xx_xz_xxz_0[i] * pa_y[i] - ta2_xx_xz_xxz_1[i] * pc_y[i];

        ta2_xx_xyz_xyy_0[i] = ta2_xx_xy_xyy_0[i] * pa_z[i] - ta2_xx_xy_xyy_1[i] * pc_z[i];

        ta2_xx_xyz_xyz_0[i] = ta2_xx_xz_xz_0[i] * fe_0 - ta2_xx_xz_xz_1[i] * fe_0 + ta2_xx_xz_xyz_0[i] * pa_y[i] - ta2_xx_xz_xyz_1[i] * pc_y[i];

        ta2_xx_xyz_xzz_0[i] = ta2_xx_xz_xzz_0[i] * pa_y[i] - ta2_xx_xz_xzz_1[i] * pc_y[i];

        ta2_xx_xyz_yyy_0[i] = ta2_xx_xy_yyy_0[i] * pa_z[i] - ta2_xx_xy_yyy_1[i] * pc_z[i];

        ta2_xx_xyz_yyz_0[i] = 2.0 * ta1_x_yz_yyz_1[i] + ta2_xx_yz_yyz_0[i] * pa_x[i] - ta2_xx_yz_yyz_1[i] * pc_x[i];

        ta2_xx_xyz_yzz_0[i] = 2.0 * ta1_x_yz_yzz_1[i] + ta2_xx_yz_yzz_0[i] * pa_x[i] - ta2_xx_yz_yzz_1[i] * pc_x[i];

        ta2_xx_xyz_zzz_0[i] = ta2_xx_xz_zzz_0[i] * pa_y[i] - ta2_xx_xz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto ta2_xx_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 50);

    auto ta2_xx_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 51);

    auto ta2_xx_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 52);

    auto ta2_xx_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 53);

    auto ta2_xx_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 54);

    auto ta2_xx_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 55);

    auto ta2_xx_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 56);

    auto ta2_xx_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 57);

    auto ta2_xx_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 58);

    auto ta2_xx_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 59);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_zz_xxz_1, ta1_x_zz_xyz_1, ta1_x_zz_xzz_1, ta1_x_zz_yyy_1, ta1_x_zz_yyz_1, ta1_x_zz_yzz_1, ta1_x_zz_zzz_1, ta2_xx_x_xxx_0, ta2_xx_x_xxx_1, ta2_xx_x_xxy_0, ta2_xx_x_xxy_1, ta2_xx_x_xyy_0, ta2_xx_x_xyy_1, ta2_xx_xz_xxx_0, ta2_xx_xz_xxx_1, ta2_xx_xz_xxy_0, ta2_xx_xz_xxy_1, ta2_xx_xz_xyy_0, ta2_xx_xz_xyy_1, ta2_xx_xzz_xxx_0, ta2_xx_xzz_xxy_0, ta2_xx_xzz_xxz_0, ta2_xx_xzz_xyy_0, ta2_xx_xzz_xyz_0, ta2_xx_xzz_xzz_0, ta2_xx_xzz_yyy_0, ta2_xx_xzz_yyz_0, ta2_xx_xzz_yzz_0, ta2_xx_xzz_zzz_0, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzz_xxx_0[i] = ta2_xx_x_xxx_0[i] * fe_0 - ta2_xx_x_xxx_1[i] * fe_0 + ta2_xx_xz_xxx_0[i] * pa_z[i] - ta2_xx_xz_xxx_1[i] * pc_z[i];

        ta2_xx_xzz_xxy_0[i] = ta2_xx_x_xxy_0[i] * fe_0 - ta2_xx_x_xxy_1[i] * fe_0 + ta2_xx_xz_xxy_0[i] * pa_z[i] - ta2_xx_xz_xxy_1[i] * pc_z[i];

        ta2_xx_xzz_xxz_0[i] = 2.0 * ta2_xx_zz_xz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxz_1[i] + ta2_xx_zz_xxz_0[i] * pa_x[i] - ta2_xx_zz_xxz_1[i] * pc_x[i];

        ta2_xx_xzz_xyy_0[i] = ta2_xx_x_xyy_0[i] * fe_0 - ta2_xx_x_xyy_1[i] * fe_0 + ta2_xx_xz_xyy_0[i] * pa_z[i] - ta2_xx_xz_xyy_1[i] * pc_z[i];

        ta2_xx_xzz_xyz_0[i] = ta2_xx_zz_yz_0[i] * fe_0 - ta2_xx_zz_yz_1[i] * fe_0 + 2.0 * ta1_x_zz_xyz_1[i] + ta2_xx_zz_xyz_0[i] * pa_x[i] - ta2_xx_zz_xyz_1[i] * pc_x[i];

        ta2_xx_xzz_xzz_0[i] = ta2_xx_zz_zz_0[i] * fe_0 - ta2_xx_zz_zz_1[i] * fe_0 + 2.0 * ta1_x_zz_xzz_1[i] + ta2_xx_zz_xzz_0[i] * pa_x[i] - ta2_xx_zz_xzz_1[i] * pc_x[i];

        ta2_xx_xzz_yyy_0[i] = 2.0 * ta1_x_zz_yyy_1[i] + ta2_xx_zz_yyy_0[i] * pa_x[i] - ta2_xx_zz_yyy_1[i] * pc_x[i];

        ta2_xx_xzz_yyz_0[i] = 2.0 * ta1_x_zz_yyz_1[i] + ta2_xx_zz_yyz_0[i] * pa_x[i] - ta2_xx_zz_yyz_1[i] * pc_x[i];

        ta2_xx_xzz_yzz_0[i] = 2.0 * ta1_x_zz_yzz_1[i] + ta2_xx_zz_yzz_0[i] * pa_x[i] - ta2_xx_zz_yzz_1[i] * pc_x[i];

        ta2_xx_xzz_zzz_0[i] = 2.0 * ta1_x_zz_zzz_1[i] + ta2_xx_zz_zzz_0[i] * pa_x[i] - ta2_xx_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto ta2_xx_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 60);

    auto ta2_xx_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 61);

    auto ta2_xx_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 62);

    auto ta2_xx_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 63);

    auto ta2_xx_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 64);

    auto ta2_xx_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 65);

    auto ta2_xx_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 66);

    auto ta2_xx_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 67);

    auto ta2_xx_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 68);

    auto ta2_xx_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 69);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_y_xxx_0, ta2_xx_y_xxx_1, ta2_xx_y_xxy_0, ta2_xx_y_xxy_1, ta2_xx_y_xxz_0, ta2_xx_y_xxz_1, ta2_xx_y_xyy_0, ta2_xx_y_xyy_1, ta2_xx_y_xyz_0, ta2_xx_y_xyz_1, ta2_xx_y_xzz_0, ta2_xx_y_xzz_1, ta2_xx_y_yyy_0, ta2_xx_y_yyy_1, ta2_xx_y_yyz_0, ta2_xx_y_yyz_1, ta2_xx_y_yzz_0, ta2_xx_y_yzz_1, ta2_xx_y_zzz_0, ta2_xx_y_zzz_1, ta2_xx_yy_xx_0, ta2_xx_yy_xx_1, ta2_xx_yy_xxx_0, ta2_xx_yy_xxx_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xxz_0, ta2_xx_yy_xxz_1, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_xz_0, ta2_xx_yy_xz_1, ta2_xx_yy_xzz_0, ta2_xx_yy_xzz_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_zz_0, ta2_xx_yy_zz_1, ta2_xx_yy_zzz_0, ta2_xx_yy_zzz_1, ta2_xx_yyy_xxx_0, ta2_xx_yyy_xxy_0, ta2_xx_yyy_xxz_0, ta2_xx_yyy_xyy_0, ta2_xx_yyy_xyz_0, ta2_xx_yyy_xzz_0, ta2_xx_yyy_yyy_0, ta2_xx_yyy_yyz_0, ta2_xx_yyy_yzz_0, ta2_xx_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyy_xxx_0[i] = 2.0 * ta2_xx_y_xxx_0[i] * fe_0 - 2.0 * ta2_xx_y_xxx_1[i] * fe_0 + ta2_xx_yy_xxx_0[i] * pa_y[i] - ta2_xx_yy_xxx_1[i] * pc_y[i];

        ta2_xx_yyy_xxy_0[i] = 2.0 * ta2_xx_y_xxy_0[i] * fe_0 - 2.0 * ta2_xx_y_xxy_1[i] * fe_0 + ta2_xx_yy_xx_0[i] * fe_0 - ta2_xx_yy_xx_1[i] * fe_0 + ta2_xx_yy_xxy_0[i] * pa_y[i] - ta2_xx_yy_xxy_1[i] * pc_y[i];

        ta2_xx_yyy_xxz_0[i] = 2.0 * ta2_xx_y_xxz_0[i] * fe_0 - 2.0 * ta2_xx_y_xxz_1[i] * fe_0 + ta2_xx_yy_xxz_0[i] * pa_y[i] - ta2_xx_yy_xxz_1[i] * pc_y[i];

        ta2_xx_yyy_xyy_0[i] = 2.0 * ta2_xx_y_xyy_0[i] * fe_0 - 2.0 * ta2_xx_y_xyy_1[i] * fe_0 + 2.0 * ta2_xx_yy_xy_0[i] * fe_0 - 2.0 * ta2_xx_yy_xy_1[i] * fe_0 + ta2_xx_yy_xyy_0[i] * pa_y[i] - ta2_xx_yy_xyy_1[i] * pc_y[i];

        ta2_xx_yyy_xyz_0[i] = 2.0 * ta2_xx_y_xyz_0[i] * fe_0 - 2.0 * ta2_xx_y_xyz_1[i] * fe_0 + ta2_xx_yy_xz_0[i] * fe_0 - ta2_xx_yy_xz_1[i] * fe_0 + ta2_xx_yy_xyz_0[i] * pa_y[i] - ta2_xx_yy_xyz_1[i] * pc_y[i];

        ta2_xx_yyy_xzz_0[i] = 2.0 * ta2_xx_y_xzz_0[i] * fe_0 - 2.0 * ta2_xx_y_xzz_1[i] * fe_0 + ta2_xx_yy_xzz_0[i] * pa_y[i] - ta2_xx_yy_xzz_1[i] * pc_y[i];

        ta2_xx_yyy_yyy_0[i] = 2.0 * ta2_xx_y_yyy_0[i] * fe_0 - 2.0 * ta2_xx_y_yyy_1[i] * fe_0 + 3.0 * ta2_xx_yy_yy_0[i] * fe_0 - 3.0 * ta2_xx_yy_yy_1[i] * fe_0 + ta2_xx_yy_yyy_0[i] * pa_y[i] - ta2_xx_yy_yyy_1[i] * pc_y[i];

        ta2_xx_yyy_yyz_0[i] = 2.0 * ta2_xx_y_yyz_0[i] * fe_0 - 2.0 * ta2_xx_y_yyz_1[i] * fe_0 + 2.0 * ta2_xx_yy_yz_0[i] * fe_0 - 2.0 * ta2_xx_yy_yz_1[i] * fe_0 + ta2_xx_yy_yyz_0[i] * pa_y[i] - ta2_xx_yy_yyz_1[i] * pc_y[i];

        ta2_xx_yyy_yzz_0[i] = 2.0 * ta2_xx_y_yzz_0[i] * fe_0 - 2.0 * ta2_xx_y_yzz_1[i] * fe_0 + ta2_xx_yy_zz_0[i] * fe_0 - ta2_xx_yy_zz_1[i] * fe_0 + ta2_xx_yy_yzz_0[i] * pa_y[i] - ta2_xx_yy_yzz_1[i] * pc_y[i];

        ta2_xx_yyy_zzz_0[i] = 2.0 * ta2_xx_y_zzz_0[i] * fe_0 - 2.0 * ta2_xx_y_zzz_1[i] * fe_0 + ta2_xx_yy_zzz_0[i] * pa_y[i] - ta2_xx_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto ta2_xx_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 70);

    auto ta2_xx_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 71);

    auto ta2_xx_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 72);

    auto ta2_xx_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 73);

    auto ta2_xx_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 74);

    auto ta2_xx_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 75);

    auto ta2_xx_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 76);

    auto ta2_xx_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 77);

    auto ta2_xx_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 78);

    auto ta2_xx_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 79);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_yy_xxx_0, ta2_xx_yy_xxx_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yyz_xxx_0, ta2_xx_yyz_xxy_0, ta2_xx_yyz_xxz_0, ta2_xx_yyz_xyy_0, ta2_xx_yyz_xyz_0, ta2_xx_yyz_xzz_0, ta2_xx_yyz_yyy_0, ta2_xx_yyz_yyz_0, ta2_xx_yyz_yzz_0, ta2_xx_yyz_zzz_0, ta2_xx_yz_xxz_0, ta2_xx_yz_xxz_1, ta2_xx_yz_xzz_0, ta2_xx_yz_xzz_1, ta2_xx_yz_zzz_0, ta2_xx_yz_zzz_1, ta2_xx_z_xxz_0, ta2_xx_z_xxz_1, ta2_xx_z_xzz_0, ta2_xx_z_xzz_1, ta2_xx_z_zzz_0, ta2_xx_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyz_xxx_0[i] = ta2_xx_yy_xxx_0[i] * pa_z[i] - ta2_xx_yy_xxx_1[i] * pc_z[i];

        ta2_xx_yyz_xxy_0[i] = ta2_xx_yy_xxy_0[i] * pa_z[i] - ta2_xx_yy_xxy_1[i] * pc_z[i];

        ta2_xx_yyz_xxz_0[i] = ta2_xx_z_xxz_0[i] * fe_0 - ta2_xx_z_xxz_1[i] * fe_0 + ta2_xx_yz_xxz_0[i] * pa_y[i] - ta2_xx_yz_xxz_1[i] * pc_y[i];

        ta2_xx_yyz_xyy_0[i] = ta2_xx_yy_xyy_0[i] * pa_z[i] - ta2_xx_yy_xyy_1[i] * pc_z[i];

        ta2_xx_yyz_xyz_0[i] = ta2_xx_yy_xy_0[i] * fe_0 - ta2_xx_yy_xy_1[i] * fe_0 + ta2_xx_yy_xyz_0[i] * pa_z[i] - ta2_xx_yy_xyz_1[i] * pc_z[i];

        ta2_xx_yyz_xzz_0[i] = ta2_xx_z_xzz_0[i] * fe_0 - ta2_xx_z_xzz_1[i] * fe_0 + ta2_xx_yz_xzz_0[i] * pa_y[i] - ta2_xx_yz_xzz_1[i] * pc_y[i];

        ta2_xx_yyz_yyy_0[i] = ta2_xx_yy_yyy_0[i] * pa_z[i] - ta2_xx_yy_yyy_1[i] * pc_z[i];

        ta2_xx_yyz_yyz_0[i] = ta2_xx_yy_yy_0[i] * fe_0 - ta2_xx_yy_yy_1[i] * fe_0 + ta2_xx_yy_yyz_0[i] * pa_z[i] - ta2_xx_yy_yyz_1[i] * pc_z[i];

        ta2_xx_yyz_yzz_0[i] = 2.0 * ta2_xx_yy_yz_0[i] * fe_0 - 2.0 * ta2_xx_yy_yz_1[i] * fe_0 + ta2_xx_yy_yzz_0[i] * pa_z[i] - ta2_xx_yy_yzz_1[i] * pc_z[i];

        ta2_xx_yyz_zzz_0[i] = ta2_xx_z_zzz_0[i] * fe_0 - ta2_xx_z_zzz_1[i] * fe_0 + ta2_xx_yz_zzz_0[i] * pa_y[i] - ta2_xx_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto ta2_xx_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 80);

    auto ta2_xx_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 81);

    auto ta2_xx_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 82);

    auto ta2_xx_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 83);

    auto ta2_xx_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 84);

    auto ta2_xx_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 85);

    auto ta2_xx_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 86);

    auto ta2_xx_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 87);

    auto ta2_xx_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 88);

    auto ta2_xx_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 89);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_yzz_xxx_0, ta2_xx_yzz_xxy_0, ta2_xx_yzz_xxz_0, ta2_xx_yzz_xyy_0, ta2_xx_yzz_xyz_0, ta2_xx_yzz_xzz_0, ta2_xx_yzz_yyy_0, ta2_xx_yzz_yyz_0, ta2_xx_yzz_yzz_0, ta2_xx_yzz_zzz_0, ta2_xx_zz_xx_0, ta2_xx_zz_xx_1, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxy_0, ta2_xx_zz_xxy_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xy_0, ta2_xx_zz_xy_1, ta2_xx_zz_xyy_0, ta2_xx_zz_xyy_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_yy_0, ta2_xx_zz_yy_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzz_xxx_0[i] = ta2_xx_zz_xxx_0[i] * pa_y[i] - ta2_xx_zz_xxx_1[i] * pc_y[i];

        ta2_xx_yzz_xxy_0[i] = ta2_xx_zz_xx_0[i] * fe_0 - ta2_xx_zz_xx_1[i] * fe_0 + ta2_xx_zz_xxy_0[i] * pa_y[i] - ta2_xx_zz_xxy_1[i] * pc_y[i];

        ta2_xx_yzz_xxz_0[i] = ta2_xx_zz_xxz_0[i] * pa_y[i] - ta2_xx_zz_xxz_1[i] * pc_y[i];

        ta2_xx_yzz_xyy_0[i] = 2.0 * ta2_xx_zz_xy_0[i] * fe_0 - 2.0 * ta2_xx_zz_xy_1[i] * fe_0 + ta2_xx_zz_xyy_0[i] * pa_y[i] - ta2_xx_zz_xyy_1[i] * pc_y[i];

        ta2_xx_yzz_xyz_0[i] = ta2_xx_zz_xz_0[i] * fe_0 - ta2_xx_zz_xz_1[i] * fe_0 + ta2_xx_zz_xyz_0[i] * pa_y[i] - ta2_xx_zz_xyz_1[i] * pc_y[i];

        ta2_xx_yzz_xzz_0[i] = ta2_xx_zz_xzz_0[i] * pa_y[i] - ta2_xx_zz_xzz_1[i] * pc_y[i];

        ta2_xx_yzz_yyy_0[i] = 3.0 * ta2_xx_zz_yy_0[i] * fe_0 - 3.0 * ta2_xx_zz_yy_1[i] * fe_0 + ta2_xx_zz_yyy_0[i] * pa_y[i] - ta2_xx_zz_yyy_1[i] * pc_y[i];

        ta2_xx_yzz_yyz_0[i] = 2.0 * ta2_xx_zz_yz_0[i] * fe_0 - 2.0 * ta2_xx_zz_yz_1[i] * fe_0 + ta2_xx_zz_yyz_0[i] * pa_y[i] - ta2_xx_zz_yyz_1[i] * pc_y[i];

        ta2_xx_yzz_yzz_0[i] = ta2_xx_zz_zz_0[i] * fe_0 - ta2_xx_zz_zz_1[i] * fe_0 + ta2_xx_zz_yzz_0[i] * pa_y[i] - ta2_xx_zz_yzz_1[i] * pc_y[i];

        ta2_xx_yzz_zzz_0[i] = ta2_xx_zz_zzz_0[i] * pa_y[i] - ta2_xx_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto ta2_xx_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 90);

    auto ta2_xx_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 91);

    auto ta2_xx_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 92);

    auto ta2_xx_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 93);

    auto ta2_xx_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 94);

    auto ta2_xx_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 95);

    auto ta2_xx_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 96);

    auto ta2_xx_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 97);

    auto ta2_xx_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 98);

    auto ta2_xx_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 99);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_z_xxx_0, ta2_xx_z_xxx_1, ta2_xx_z_xxy_0, ta2_xx_z_xxy_1, ta2_xx_z_xxz_0, ta2_xx_z_xxz_1, ta2_xx_z_xyy_0, ta2_xx_z_xyy_1, ta2_xx_z_xyz_0, ta2_xx_z_xyz_1, ta2_xx_z_xzz_0, ta2_xx_z_xzz_1, ta2_xx_z_yyy_0, ta2_xx_z_yyy_1, ta2_xx_z_yyz_0, ta2_xx_z_yyz_1, ta2_xx_z_yzz_0, ta2_xx_z_yzz_1, ta2_xx_z_zzz_0, ta2_xx_z_zzz_1, ta2_xx_zz_xx_0, ta2_xx_zz_xx_1, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxy_0, ta2_xx_zz_xxy_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xy_0, ta2_xx_zz_xy_1, ta2_xx_zz_xyy_0, ta2_xx_zz_xyy_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_yy_0, ta2_xx_zz_yy_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, ta2_xx_zzz_xxx_0, ta2_xx_zzz_xxy_0, ta2_xx_zzz_xxz_0, ta2_xx_zzz_xyy_0, ta2_xx_zzz_xyz_0, ta2_xx_zzz_xzz_0, ta2_xx_zzz_yyy_0, ta2_xx_zzz_yyz_0, ta2_xx_zzz_yzz_0, ta2_xx_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzz_xxx_0[i] = 2.0 * ta2_xx_z_xxx_0[i] * fe_0 - 2.0 * ta2_xx_z_xxx_1[i] * fe_0 + ta2_xx_zz_xxx_0[i] * pa_z[i] - ta2_xx_zz_xxx_1[i] * pc_z[i];

        ta2_xx_zzz_xxy_0[i] = 2.0 * ta2_xx_z_xxy_0[i] * fe_0 - 2.0 * ta2_xx_z_xxy_1[i] * fe_0 + ta2_xx_zz_xxy_0[i] * pa_z[i] - ta2_xx_zz_xxy_1[i] * pc_z[i];

        ta2_xx_zzz_xxz_0[i] = 2.0 * ta2_xx_z_xxz_0[i] * fe_0 - 2.0 * ta2_xx_z_xxz_1[i] * fe_0 + ta2_xx_zz_xx_0[i] * fe_0 - ta2_xx_zz_xx_1[i] * fe_0 + ta2_xx_zz_xxz_0[i] * pa_z[i] - ta2_xx_zz_xxz_1[i] * pc_z[i];

        ta2_xx_zzz_xyy_0[i] = 2.0 * ta2_xx_z_xyy_0[i] * fe_0 - 2.0 * ta2_xx_z_xyy_1[i] * fe_0 + ta2_xx_zz_xyy_0[i] * pa_z[i] - ta2_xx_zz_xyy_1[i] * pc_z[i];

        ta2_xx_zzz_xyz_0[i] = 2.0 * ta2_xx_z_xyz_0[i] * fe_0 - 2.0 * ta2_xx_z_xyz_1[i] * fe_0 + ta2_xx_zz_xy_0[i] * fe_0 - ta2_xx_zz_xy_1[i] * fe_0 + ta2_xx_zz_xyz_0[i] * pa_z[i] - ta2_xx_zz_xyz_1[i] * pc_z[i];

        ta2_xx_zzz_xzz_0[i] = 2.0 * ta2_xx_z_xzz_0[i] * fe_0 - 2.0 * ta2_xx_z_xzz_1[i] * fe_0 + 2.0 * ta2_xx_zz_xz_0[i] * fe_0 - 2.0 * ta2_xx_zz_xz_1[i] * fe_0 + ta2_xx_zz_xzz_0[i] * pa_z[i] - ta2_xx_zz_xzz_1[i] * pc_z[i];

        ta2_xx_zzz_yyy_0[i] = 2.0 * ta2_xx_z_yyy_0[i] * fe_0 - 2.0 * ta2_xx_z_yyy_1[i] * fe_0 + ta2_xx_zz_yyy_0[i] * pa_z[i] - ta2_xx_zz_yyy_1[i] * pc_z[i];

        ta2_xx_zzz_yyz_0[i] = 2.0 * ta2_xx_z_yyz_0[i] * fe_0 - 2.0 * ta2_xx_z_yyz_1[i] * fe_0 + ta2_xx_zz_yy_0[i] * fe_0 - ta2_xx_zz_yy_1[i] * fe_0 + ta2_xx_zz_yyz_0[i] * pa_z[i] - ta2_xx_zz_yyz_1[i] * pc_z[i];

        ta2_xx_zzz_yzz_0[i] = 2.0 * ta2_xx_z_yzz_0[i] * fe_0 - 2.0 * ta2_xx_z_yzz_1[i] * fe_0 + 2.0 * ta2_xx_zz_yz_0[i] * fe_0 - 2.0 * ta2_xx_zz_yz_1[i] * fe_0 + ta2_xx_zz_yzz_0[i] * pa_z[i] - ta2_xx_zz_yzz_1[i] * pc_z[i];

        ta2_xx_zzz_zzz_0[i] = 2.0 * ta2_xx_z_zzz_0[i] * fe_0 - 2.0 * ta2_xx_z_zzz_1[i] * fe_0 + 3.0 * ta2_xx_zz_zz_0[i] * fe_0 - 3.0 * ta2_xx_zz_zz_1[i] * fe_0 + ta2_xx_zz_zzz_0[i] * pa_z[i] - ta2_xx_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 100-110 components of targeted buffer : FF

    auto ta2_xy_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 100);

    auto ta2_xy_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 101);

    auto ta2_xy_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 102);

    auto ta2_xy_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 103);

    auto ta2_xy_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 104);

    auto ta2_xy_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 105);

    auto ta2_xy_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 106);

    auto ta2_xy_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 107);

    auto ta2_xy_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 108);

    auto ta2_xy_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 109);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xx_xxx_1, ta1_y_xx_xxy_1, ta1_y_xx_xxz_1, ta1_y_xx_xyy_1, ta1_y_xx_xyz_1, ta1_y_xx_xzz_1, ta1_y_xx_yyy_1, ta1_y_xx_yyz_1, ta1_y_xx_yzz_1, ta1_y_xx_zzz_1, ta2_xy_x_xxx_0, ta2_xy_x_xxx_1, ta2_xy_x_xxy_0, ta2_xy_x_xxy_1, ta2_xy_x_xxz_0, ta2_xy_x_xxz_1, ta2_xy_x_xyy_0, ta2_xy_x_xyy_1, ta2_xy_x_xyz_0, ta2_xy_x_xyz_1, ta2_xy_x_xzz_0, ta2_xy_x_xzz_1, ta2_xy_x_yyy_0, ta2_xy_x_yyy_1, ta2_xy_x_yyz_0, ta2_xy_x_yyz_1, ta2_xy_x_yzz_0, ta2_xy_x_yzz_1, ta2_xy_x_zzz_0, ta2_xy_x_zzz_1, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_yy_0, ta2_xy_xx_yy_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xx_yyz_0, ta2_xy_xx_yyz_1, ta2_xy_xx_yz_0, ta2_xy_xx_yz_1, ta2_xy_xx_yzz_0, ta2_xy_xx_yzz_1, ta2_xy_xx_zz_0, ta2_xy_xx_zz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xxx_xxx_0, ta2_xy_xxx_xxy_0, ta2_xy_xxx_xxz_0, ta2_xy_xxx_xyy_0, ta2_xy_xxx_xyz_0, ta2_xy_xxx_xzz_0, ta2_xy_xxx_yyy_0, ta2_xy_xxx_yyz_0, ta2_xy_xxx_yzz_0, ta2_xy_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxx_xxx_0[i] = 2.0 * ta2_xy_x_xxx_0[i] * fe_0 - 2.0 * ta2_xy_x_xxx_1[i] * fe_0 + 3.0 * ta2_xy_xx_xx_0[i] * fe_0 - 3.0 * ta2_xy_xx_xx_1[i] * fe_0 + ta1_y_xx_xxx_1[i] + ta2_xy_xx_xxx_0[i] * pa_x[i] - ta2_xy_xx_xxx_1[i] * pc_x[i];

        ta2_xy_xxx_xxy_0[i] = 2.0 * ta2_xy_x_xxy_0[i] * fe_0 - 2.0 * ta2_xy_x_xxy_1[i] * fe_0 + 2.0 * ta2_xy_xx_xy_0[i] * fe_0 - 2.0 * ta2_xy_xx_xy_1[i] * fe_0 + ta1_y_xx_xxy_1[i] + ta2_xy_xx_xxy_0[i] * pa_x[i] - ta2_xy_xx_xxy_1[i] * pc_x[i];

        ta2_xy_xxx_xxz_0[i] = 2.0 * ta2_xy_x_xxz_0[i] * fe_0 - 2.0 * ta2_xy_x_xxz_1[i] * fe_0 + 2.0 * ta2_xy_xx_xz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xz_1[i] * fe_0 + ta1_y_xx_xxz_1[i] + ta2_xy_xx_xxz_0[i] * pa_x[i] - ta2_xy_xx_xxz_1[i] * pc_x[i];

        ta2_xy_xxx_xyy_0[i] = 2.0 * ta2_xy_x_xyy_0[i] * fe_0 - 2.0 * ta2_xy_x_xyy_1[i] * fe_0 + ta2_xy_xx_yy_0[i] * fe_0 - ta2_xy_xx_yy_1[i] * fe_0 + ta1_y_xx_xyy_1[i] + ta2_xy_xx_xyy_0[i] * pa_x[i] - ta2_xy_xx_xyy_1[i] * pc_x[i];

        ta2_xy_xxx_xyz_0[i] = 2.0 * ta2_xy_x_xyz_0[i] * fe_0 - 2.0 * ta2_xy_x_xyz_1[i] * fe_0 + ta2_xy_xx_yz_0[i] * fe_0 - ta2_xy_xx_yz_1[i] * fe_0 + ta1_y_xx_xyz_1[i] + ta2_xy_xx_xyz_0[i] * pa_x[i] - ta2_xy_xx_xyz_1[i] * pc_x[i];

        ta2_xy_xxx_xzz_0[i] = 2.0 * ta2_xy_x_xzz_0[i] * fe_0 - 2.0 * ta2_xy_x_xzz_1[i] * fe_0 + ta2_xy_xx_zz_0[i] * fe_0 - ta2_xy_xx_zz_1[i] * fe_0 + ta1_y_xx_xzz_1[i] + ta2_xy_xx_xzz_0[i] * pa_x[i] - ta2_xy_xx_xzz_1[i] * pc_x[i];

        ta2_xy_xxx_yyy_0[i] = 2.0 * ta2_xy_x_yyy_0[i] * fe_0 - 2.0 * ta2_xy_x_yyy_1[i] * fe_0 + ta1_y_xx_yyy_1[i] + ta2_xy_xx_yyy_0[i] * pa_x[i] - ta2_xy_xx_yyy_1[i] * pc_x[i];

        ta2_xy_xxx_yyz_0[i] = 2.0 * ta2_xy_x_yyz_0[i] * fe_0 - 2.0 * ta2_xy_x_yyz_1[i] * fe_0 + ta1_y_xx_yyz_1[i] + ta2_xy_xx_yyz_0[i] * pa_x[i] - ta2_xy_xx_yyz_1[i] * pc_x[i];

        ta2_xy_xxx_yzz_0[i] = 2.0 * ta2_xy_x_yzz_0[i] * fe_0 - 2.0 * ta2_xy_x_yzz_1[i] * fe_0 + ta1_y_xx_yzz_1[i] + ta2_xy_xx_yzz_0[i] * pa_x[i] - ta2_xy_xx_yzz_1[i] * pc_x[i];

        ta2_xy_xxx_zzz_0[i] = 2.0 * ta2_xy_x_zzz_0[i] * fe_0 - 2.0 * ta2_xy_x_zzz_1[i] * fe_0 + ta1_y_xx_zzz_1[i] + ta2_xy_xx_zzz_0[i] * pa_x[i] - ta2_xy_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 110-120 components of targeted buffer : FF

    auto ta2_xy_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 110);

    auto ta2_xy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 111);

    auto ta2_xy_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 112);

    auto ta2_xy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 113);

    auto ta2_xy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 114);

    auto ta2_xy_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 115);

    auto ta2_xy_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 116);

    auto ta2_xy_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 117);

    auto ta2_xy_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 118);

    auto ta2_xy_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 119);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xx_xxx_1, ta1_x_xx_xxy_1, ta1_x_xx_xxz_1, ta1_x_xx_xyy_1, ta1_x_xx_xyz_1, ta1_x_xx_xzz_1, ta1_x_xx_zzz_1, ta1_y_xy_yyy_1, ta1_y_xy_yyz_1, ta1_y_xy_yzz_1, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xxy_xxx_0, ta2_xy_xxy_xxy_0, ta2_xy_xxy_xxz_0, ta2_xy_xxy_xyy_0, ta2_xy_xxy_xyz_0, ta2_xy_xxy_xzz_0, ta2_xy_xxy_yyy_0, ta2_xy_xxy_yyz_0, ta2_xy_xxy_yzz_0, ta2_xy_xxy_zzz_0, ta2_xy_xy_yyy_0, ta2_xy_xy_yyy_1, ta2_xy_xy_yyz_0, ta2_xy_xy_yyz_1, ta2_xy_xy_yzz_0, ta2_xy_xy_yzz_1, ta2_xy_y_yyy_0, ta2_xy_y_yyy_1, ta2_xy_y_yyz_0, ta2_xy_y_yyz_1, ta2_xy_y_yzz_0, ta2_xy_y_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxy_xxx_0[i] = ta1_x_xx_xxx_1[i] + ta2_xy_xx_xxx_0[i] * pa_y[i] - ta2_xy_xx_xxx_1[i] * pc_y[i];

        ta2_xy_xxy_xxy_0[i] = ta2_xy_xx_xx_0[i] * fe_0 - ta2_xy_xx_xx_1[i] * fe_0 + ta1_x_xx_xxy_1[i] + ta2_xy_xx_xxy_0[i] * pa_y[i] - ta2_xy_xx_xxy_1[i] * pc_y[i];

        ta2_xy_xxy_xxz_0[i] = ta1_x_xx_xxz_1[i] + ta2_xy_xx_xxz_0[i] * pa_y[i] - ta2_xy_xx_xxz_1[i] * pc_y[i];

        ta2_xy_xxy_xyy_0[i] = 2.0 * ta2_xy_xx_xy_0[i] * fe_0 - 2.0 * ta2_xy_xx_xy_1[i] * fe_0 + ta1_x_xx_xyy_1[i] + ta2_xy_xx_xyy_0[i] * pa_y[i] - ta2_xy_xx_xyy_1[i] * pc_y[i];

        ta2_xy_xxy_xyz_0[i] = ta2_xy_xx_xz_0[i] * fe_0 - ta2_xy_xx_xz_1[i] * fe_0 + ta1_x_xx_xyz_1[i] + ta2_xy_xx_xyz_0[i] * pa_y[i] - ta2_xy_xx_xyz_1[i] * pc_y[i];

        ta2_xy_xxy_xzz_0[i] = ta1_x_xx_xzz_1[i] + ta2_xy_xx_xzz_0[i] * pa_y[i] - ta2_xy_xx_xzz_1[i] * pc_y[i];

        ta2_xy_xxy_yyy_0[i] = ta2_xy_y_yyy_0[i] * fe_0 - ta2_xy_y_yyy_1[i] * fe_0 + ta1_y_xy_yyy_1[i] + ta2_xy_xy_yyy_0[i] * pa_x[i] - ta2_xy_xy_yyy_1[i] * pc_x[i];

        ta2_xy_xxy_yyz_0[i] = ta2_xy_y_yyz_0[i] * fe_0 - ta2_xy_y_yyz_1[i] * fe_0 + ta1_y_xy_yyz_1[i] + ta2_xy_xy_yyz_0[i] * pa_x[i] - ta2_xy_xy_yyz_1[i] * pc_x[i];

        ta2_xy_xxy_yzz_0[i] = ta2_xy_y_yzz_0[i] * fe_0 - ta2_xy_y_yzz_1[i] * fe_0 + ta1_y_xy_yzz_1[i] + ta2_xy_xy_yzz_0[i] * pa_x[i] - ta2_xy_xy_yzz_1[i] * pc_x[i];

        ta2_xy_xxy_zzz_0[i] = ta1_x_xx_zzz_1[i] + ta2_xy_xx_zzz_0[i] * pa_y[i] - ta2_xy_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : FF

    auto ta2_xy_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 120);

    auto ta2_xy_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 121);

    auto ta2_xy_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 122);

    auto ta2_xy_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 123);

    auto ta2_xy_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 124);

    auto ta2_xy_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 125);

    auto ta2_xy_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 126);

    auto ta2_xy_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 127);

    auto ta2_xy_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 128);

    auto ta2_xy_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 129);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_yy_0, ta2_xy_xx_yy_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xx_yyz_0, ta2_xy_xx_yyz_1, ta2_xy_xx_yz_0, ta2_xy_xx_yz_1, ta2_xy_xx_yzz_0, ta2_xy_xx_yzz_1, ta2_xy_xx_zz_0, ta2_xy_xx_zz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xxz_xxx_0, ta2_xy_xxz_xxy_0, ta2_xy_xxz_xxz_0, ta2_xy_xxz_xyy_0, ta2_xy_xxz_xyz_0, ta2_xy_xxz_xzz_0, ta2_xy_xxz_yyy_0, ta2_xy_xxz_yyz_0, ta2_xy_xxz_yzz_0, ta2_xy_xxz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxz_xxx_0[i] = ta2_xy_xx_xxx_0[i] * pa_z[i] - ta2_xy_xx_xxx_1[i] * pc_z[i];

        ta2_xy_xxz_xxy_0[i] = ta2_xy_xx_xxy_0[i] * pa_z[i] - ta2_xy_xx_xxy_1[i] * pc_z[i];

        ta2_xy_xxz_xxz_0[i] = ta2_xy_xx_xx_0[i] * fe_0 - ta2_xy_xx_xx_1[i] * fe_0 + ta2_xy_xx_xxz_0[i] * pa_z[i] - ta2_xy_xx_xxz_1[i] * pc_z[i];

        ta2_xy_xxz_xyy_0[i] = ta2_xy_xx_xyy_0[i] * pa_z[i] - ta2_xy_xx_xyy_1[i] * pc_z[i];

        ta2_xy_xxz_xyz_0[i] = ta2_xy_xx_xy_0[i] * fe_0 - ta2_xy_xx_xy_1[i] * fe_0 + ta2_xy_xx_xyz_0[i] * pa_z[i] - ta2_xy_xx_xyz_1[i] * pc_z[i];

        ta2_xy_xxz_xzz_0[i] = 2.0 * ta2_xy_xx_xz_0[i] * fe_0 - 2.0 * ta2_xy_xx_xz_1[i] * fe_0 + ta2_xy_xx_xzz_0[i] * pa_z[i] - ta2_xy_xx_xzz_1[i] * pc_z[i];

        ta2_xy_xxz_yyy_0[i] = ta2_xy_xx_yyy_0[i] * pa_z[i] - ta2_xy_xx_yyy_1[i] * pc_z[i];

        ta2_xy_xxz_yyz_0[i] = ta2_xy_xx_yy_0[i] * fe_0 - ta2_xy_xx_yy_1[i] * fe_0 + ta2_xy_xx_yyz_0[i] * pa_z[i] - ta2_xy_xx_yyz_1[i] * pc_z[i];

        ta2_xy_xxz_yzz_0[i] = 2.0 * ta2_xy_xx_yz_0[i] * fe_0 - 2.0 * ta2_xy_xx_yz_1[i] * fe_0 + ta2_xy_xx_yzz_0[i] * pa_z[i] - ta2_xy_xx_yzz_1[i] * pc_z[i];

        ta2_xy_xxz_zzz_0[i] = 3.0 * ta2_xy_xx_zz_0[i] * fe_0 - 3.0 * ta2_xy_xx_zz_1[i] * fe_0 + ta2_xy_xx_zzz_0[i] * pa_z[i] - ta2_xy_xx_zzz_1[i] * pc_z[i];
    }

    // Set up 130-140 components of targeted buffer : FF

    auto ta2_xy_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 130);

    auto ta2_xy_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 131);

    auto ta2_xy_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 132);

    auto ta2_xy_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 133);

    auto ta2_xy_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 134);

    auto ta2_xy_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 135);

    auto ta2_xy_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 136);

    auto ta2_xy_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 137);

    auto ta2_xy_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 138);

    auto ta2_xy_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 139);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_yy_xxx_1, ta1_y_yy_xxy_1, ta1_y_yy_xxz_1, ta1_y_yy_xyy_1, ta1_y_yy_xyz_1, ta1_y_yy_xzz_1, ta1_y_yy_yyy_1, ta1_y_yy_yyz_1, ta1_y_yy_yzz_1, ta1_y_yy_zzz_1, ta2_xy_xyy_xxx_0, ta2_xy_xyy_xxy_0, ta2_xy_xyy_xxz_0, ta2_xy_xyy_xyy_0, ta2_xy_xyy_xyz_0, ta2_xy_xyy_xzz_0, ta2_xy_xyy_yyy_0, ta2_xy_xyy_yyz_0, ta2_xy_xyy_yzz_0, ta2_xy_xyy_zzz_0, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyy_xxx_0[i] = 3.0 * ta2_xy_yy_xx_0[i] * fe_0 - 3.0 * ta2_xy_yy_xx_1[i] * fe_0 + ta1_y_yy_xxx_1[i] + ta2_xy_yy_xxx_0[i] * pa_x[i] - ta2_xy_yy_xxx_1[i] * pc_x[i];

        ta2_xy_xyy_xxy_0[i] = 2.0 * ta2_xy_yy_xy_0[i] * fe_0 - 2.0 * ta2_xy_yy_xy_1[i] * fe_0 + ta1_y_yy_xxy_1[i] + ta2_xy_yy_xxy_0[i] * pa_x[i] - ta2_xy_yy_xxy_1[i] * pc_x[i];

        ta2_xy_xyy_xxz_0[i] = 2.0 * ta2_xy_yy_xz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xz_1[i] * fe_0 + ta1_y_yy_xxz_1[i] + ta2_xy_yy_xxz_0[i] * pa_x[i] - ta2_xy_yy_xxz_1[i] * pc_x[i];

        ta2_xy_xyy_xyy_0[i] = ta2_xy_yy_yy_0[i] * fe_0 - ta2_xy_yy_yy_1[i] * fe_0 + ta1_y_yy_xyy_1[i] + ta2_xy_yy_xyy_0[i] * pa_x[i] - ta2_xy_yy_xyy_1[i] * pc_x[i];

        ta2_xy_xyy_xyz_0[i] = ta2_xy_yy_yz_0[i] * fe_0 - ta2_xy_yy_yz_1[i] * fe_0 + ta1_y_yy_xyz_1[i] + ta2_xy_yy_xyz_0[i] * pa_x[i] - ta2_xy_yy_xyz_1[i] * pc_x[i];

        ta2_xy_xyy_xzz_0[i] = ta2_xy_yy_zz_0[i] * fe_0 - ta2_xy_yy_zz_1[i] * fe_0 + ta1_y_yy_xzz_1[i] + ta2_xy_yy_xzz_0[i] * pa_x[i] - ta2_xy_yy_xzz_1[i] * pc_x[i];

        ta2_xy_xyy_yyy_0[i] = ta1_y_yy_yyy_1[i] + ta2_xy_yy_yyy_0[i] * pa_x[i] - ta2_xy_yy_yyy_1[i] * pc_x[i];

        ta2_xy_xyy_yyz_0[i] = ta1_y_yy_yyz_1[i] + ta2_xy_yy_yyz_0[i] * pa_x[i] - ta2_xy_yy_yyz_1[i] * pc_x[i];

        ta2_xy_xyy_yzz_0[i] = ta1_y_yy_yzz_1[i] + ta2_xy_yy_yzz_0[i] * pa_x[i] - ta2_xy_yy_yzz_1[i] * pc_x[i];

        ta2_xy_xyy_zzz_0[i] = ta1_y_yy_zzz_1[i] + ta2_xy_yy_zzz_0[i] * pa_x[i] - ta2_xy_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 140-150 components of targeted buffer : FF

    auto ta2_xy_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 140);

    auto ta2_xy_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 141);

    auto ta2_xy_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 142);

    auto ta2_xy_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 143);

    auto ta2_xy_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 144);

    auto ta2_xy_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 145);

    auto ta2_xy_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 146);

    auto ta2_xy_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 147);

    auto ta2_xy_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 148);

    auto ta2_xy_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 149);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xz_xxz_1, ta1_x_xz_xzz_1, ta1_y_yz_yyz_1, ta1_y_yz_yzz_1, ta1_y_yz_zzz_1, ta2_xy_xy_xxx_0, ta2_xy_xy_xxx_1, ta2_xy_xy_xxy_0, ta2_xy_xy_xxy_1, ta2_xy_xy_xy_0, ta2_xy_xy_xy_1, ta2_xy_xy_xyy_0, ta2_xy_xy_xyy_1, ta2_xy_xy_xyz_0, ta2_xy_xy_xyz_1, ta2_xy_xy_yyy_0, ta2_xy_xy_yyy_1, ta2_xy_xyz_xxx_0, ta2_xy_xyz_xxy_0, ta2_xy_xyz_xxz_0, ta2_xy_xyz_xyy_0, ta2_xy_xyz_xyz_0, ta2_xy_xyz_xzz_0, ta2_xy_xyz_yyy_0, ta2_xy_xyz_yyz_0, ta2_xy_xyz_yzz_0, ta2_xy_xyz_zzz_0, ta2_xy_xz_xxz_0, ta2_xy_xz_xxz_1, ta2_xy_xz_xzz_0, ta2_xy_xz_xzz_1, ta2_xy_yz_yyz_0, ta2_xy_yz_yyz_1, ta2_xy_yz_yzz_0, ta2_xy_yz_yzz_1, ta2_xy_yz_zzz_0, ta2_xy_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyz_xxx_0[i] = ta2_xy_xy_xxx_0[i] * pa_z[i] - ta2_xy_xy_xxx_1[i] * pc_z[i];

        ta2_xy_xyz_xxy_0[i] = ta2_xy_xy_xxy_0[i] * pa_z[i] - ta2_xy_xy_xxy_1[i] * pc_z[i];

        ta2_xy_xyz_xxz_0[i] = ta1_x_xz_xxz_1[i] + ta2_xy_xz_xxz_0[i] * pa_y[i] - ta2_xy_xz_xxz_1[i] * pc_y[i];

        ta2_xy_xyz_xyy_0[i] = ta2_xy_xy_xyy_0[i] * pa_z[i] - ta2_xy_xy_xyy_1[i] * pc_z[i];

        ta2_xy_xyz_xyz_0[i] = ta2_xy_xy_xy_0[i] * fe_0 - ta2_xy_xy_xy_1[i] * fe_0 + ta2_xy_xy_xyz_0[i] * pa_z[i] - ta2_xy_xy_xyz_1[i] * pc_z[i];

        ta2_xy_xyz_xzz_0[i] = ta1_x_xz_xzz_1[i] + ta2_xy_xz_xzz_0[i] * pa_y[i] - ta2_xy_xz_xzz_1[i] * pc_y[i];

        ta2_xy_xyz_yyy_0[i] = ta2_xy_xy_yyy_0[i] * pa_z[i] - ta2_xy_xy_yyy_1[i] * pc_z[i];

        ta2_xy_xyz_yyz_0[i] = ta1_y_yz_yyz_1[i] + ta2_xy_yz_yyz_0[i] * pa_x[i] - ta2_xy_yz_yyz_1[i] * pc_x[i];

        ta2_xy_xyz_yzz_0[i] = ta1_y_yz_yzz_1[i] + ta2_xy_yz_yzz_0[i] * pa_x[i] - ta2_xy_yz_yzz_1[i] * pc_x[i];

        ta2_xy_xyz_zzz_0[i] = ta1_y_yz_zzz_1[i] + ta2_xy_yz_zzz_0[i] * pa_x[i] - ta2_xy_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : FF

    auto ta2_xy_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 150);

    auto ta2_xy_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 151);

    auto ta2_xy_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 152);

    auto ta2_xy_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 153);

    auto ta2_xy_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 154);

    auto ta2_xy_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 155);

    auto ta2_xy_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 156);

    auto ta2_xy_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 157);

    auto ta2_xy_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 158);

    auto ta2_xy_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 159);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_zz_xxz_1, ta1_y_zz_xyz_1, ta1_y_zz_xzz_1, ta1_y_zz_yyy_1, ta1_y_zz_yyz_1, ta1_y_zz_yzz_1, ta1_y_zz_zzz_1, ta2_xy_x_xxx_0, ta2_xy_x_xxx_1, ta2_xy_x_xxy_0, ta2_xy_x_xxy_1, ta2_xy_x_xyy_0, ta2_xy_x_xyy_1, ta2_xy_xz_xxx_0, ta2_xy_xz_xxx_1, ta2_xy_xz_xxy_0, ta2_xy_xz_xxy_1, ta2_xy_xz_xyy_0, ta2_xy_xz_xyy_1, ta2_xy_xzz_xxx_0, ta2_xy_xzz_xxy_0, ta2_xy_xzz_xxz_0, ta2_xy_xzz_xyy_0, ta2_xy_xzz_xyz_0, ta2_xy_xzz_xzz_0, ta2_xy_xzz_yyy_0, ta2_xy_xzz_yyz_0, ta2_xy_xzz_yzz_0, ta2_xy_xzz_zzz_0, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_yyy_0, ta2_xy_zz_yyy_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzz_xxx_0[i] = ta2_xy_x_xxx_0[i] * fe_0 - ta2_xy_x_xxx_1[i] * fe_0 + ta2_xy_xz_xxx_0[i] * pa_z[i] - ta2_xy_xz_xxx_1[i] * pc_z[i];

        ta2_xy_xzz_xxy_0[i] = ta2_xy_x_xxy_0[i] * fe_0 - ta2_xy_x_xxy_1[i] * fe_0 + ta2_xy_xz_xxy_0[i] * pa_z[i] - ta2_xy_xz_xxy_1[i] * pc_z[i];

        ta2_xy_xzz_xxz_0[i] = 2.0 * ta2_xy_zz_xz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xz_1[i] * fe_0 + ta1_y_zz_xxz_1[i] + ta2_xy_zz_xxz_0[i] * pa_x[i] - ta2_xy_zz_xxz_1[i] * pc_x[i];

        ta2_xy_xzz_xyy_0[i] = ta2_xy_x_xyy_0[i] * fe_0 - ta2_xy_x_xyy_1[i] * fe_0 + ta2_xy_xz_xyy_0[i] * pa_z[i] - ta2_xy_xz_xyy_1[i] * pc_z[i];

        ta2_xy_xzz_xyz_0[i] = ta2_xy_zz_yz_0[i] * fe_0 - ta2_xy_zz_yz_1[i] * fe_0 + ta1_y_zz_xyz_1[i] + ta2_xy_zz_xyz_0[i] * pa_x[i] - ta2_xy_zz_xyz_1[i] * pc_x[i];

        ta2_xy_xzz_xzz_0[i] = ta2_xy_zz_zz_0[i] * fe_0 - ta2_xy_zz_zz_1[i] * fe_0 + ta1_y_zz_xzz_1[i] + ta2_xy_zz_xzz_0[i] * pa_x[i] - ta2_xy_zz_xzz_1[i] * pc_x[i];

        ta2_xy_xzz_yyy_0[i] = ta1_y_zz_yyy_1[i] + ta2_xy_zz_yyy_0[i] * pa_x[i] - ta2_xy_zz_yyy_1[i] * pc_x[i];

        ta2_xy_xzz_yyz_0[i] = ta1_y_zz_yyz_1[i] + ta2_xy_zz_yyz_0[i] * pa_x[i] - ta2_xy_zz_yyz_1[i] * pc_x[i];

        ta2_xy_xzz_yzz_0[i] = ta1_y_zz_yzz_1[i] + ta2_xy_zz_yzz_0[i] * pa_x[i] - ta2_xy_zz_yzz_1[i] * pc_x[i];

        ta2_xy_xzz_zzz_0[i] = ta1_y_zz_zzz_1[i] + ta2_xy_zz_zzz_0[i] * pa_x[i] - ta2_xy_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : FF

    auto ta2_xy_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 160);

    auto ta2_xy_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 161);

    auto ta2_xy_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 162);

    auto ta2_xy_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 163);

    auto ta2_xy_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 164);

    auto ta2_xy_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 165);

    auto ta2_xy_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 166);

    auto ta2_xy_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 167);

    auto ta2_xy_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 168);

    auto ta2_xy_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 169);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yy_xxx_1, ta1_x_yy_xxy_1, ta1_x_yy_xxz_1, ta1_x_yy_xyy_1, ta1_x_yy_xyz_1, ta1_x_yy_xzz_1, ta1_x_yy_yyy_1, ta1_x_yy_yyz_1, ta1_x_yy_yzz_1, ta1_x_yy_zzz_1, ta2_xy_y_xxx_0, ta2_xy_y_xxx_1, ta2_xy_y_xxy_0, ta2_xy_y_xxy_1, ta2_xy_y_xxz_0, ta2_xy_y_xxz_1, ta2_xy_y_xyy_0, ta2_xy_y_xyy_1, ta2_xy_y_xyz_0, ta2_xy_y_xyz_1, ta2_xy_y_xzz_0, ta2_xy_y_xzz_1, ta2_xy_y_yyy_0, ta2_xy_y_yyy_1, ta2_xy_y_yyz_0, ta2_xy_y_yyz_1, ta2_xy_y_yzz_0, ta2_xy_y_yzz_1, ta2_xy_y_zzz_0, ta2_xy_y_zzz_1, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yyy_xxx_0, ta2_xy_yyy_xxy_0, ta2_xy_yyy_xxz_0, ta2_xy_yyy_xyy_0, ta2_xy_yyy_xyz_0, ta2_xy_yyy_xzz_0, ta2_xy_yyy_yyy_0, ta2_xy_yyy_yyz_0, ta2_xy_yyy_yzz_0, ta2_xy_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyy_xxx_0[i] = 2.0 * ta2_xy_y_xxx_0[i] * fe_0 - 2.0 * ta2_xy_y_xxx_1[i] * fe_0 + ta1_x_yy_xxx_1[i] + ta2_xy_yy_xxx_0[i] * pa_y[i] - ta2_xy_yy_xxx_1[i] * pc_y[i];

        ta2_xy_yyy_xxy_0[i] = 2.0 * ta2_xy_y_xxy_0[i] * fe_0 - 2.0 * ta2_xy_y_xxy_1[i] * fe_0 + ta2_xy_yy_xx_0[i] * fe_0 - ta2_xy_yy_xx_1[i] * fe_0 + ta1_x_yy_xxy_1[i] + ta2_xy_yy_xxy_0[i] * pa_y[i] - ta2_xy_yy_xxy_1[i] * pc_y[i];

        ta2_xy_yyy_xxz_0[i] = 2.0 * ta2_xy_y_xxz_0[i] * fe_0 - 2.0 * ta2_xy_y_xxz_1[i] * fe_0 + ta1_x_yy_xxz_1[i] + ta2_xy_yy_xxz_0[i] * pa_y[i] - ta2_xy_yy_xxz_1[i] * pc_y[i];

        ta2_xy_yyy_xyy_0[i] = 2.0 * ta2_xy_y_xyy_0[i] * fe_0 - 2.0 * ta2_xy_y_xyy_1[i] * fe_0 + 2.0 * ta2_xy_yy_xy_0[i] * fe_0 - 2.0 * ta2_xy_yy_xy_1[i] * fe_0 + ta1_x_yy_xyy_1[i] + ta2_xy_yy_xyy_0[i] * pa_y[i] - ta2_xy_yy_xyy_1[i] * pc_y[i];

        ta2_xy_yyy_xyz_0[i] = 2.0 * ta2_xy_y_xyz_0[i] * fe_0 - 2.0 * ta2_xy_y_xyz_1[i] * fe_0 + ta2_xy_yy_xz_0[i] * fe_0 - ta2_xy_yy_xz_1[i] * fe_0 + ta1_x_yy_xyz_1[i] + ta2_xy_yy_xyz_0[i] * pa_y[i] - ta2_xy_yy_xyz_1[i] * pc_y[i];

        ta2_xy_yyy_xzz_0[i] = 2.0 * ta2_xy_y_xzz_0[i] * fe_0 - 2.0 * ta2_xy_y_xzz_1[i] * fe_0 + ta1_x_yy_xzz_1[i] + ta2_xy_yy_xzz_0[i] * pa_y[i] - ta2_xy_yy_xzz_1[i] * pc_y[i];

        ta2_xy_yyy_yyy_0[i] = 2.0 * ta2_xy_y_yyy_0[i] * fe_0 - 2.0 * ta2_xy_y_yyy_1[i] * fe_0 + 3.0 * ta2_xy_yy_yy_0[i] * fe_0 - 3.0 * ta2_xy_yy_yy_1[i] * fe_0 + ta1_x_yy_yyy_1[i] + ta2_xy_yy_yyy_0[i] * pa_y[i] - ta2_xy_yy_yyy_1[i] * pc_y[i];

        ta2_xy_yyy_yyz_0[i] = 2.0 * ta2_xy_y_yyz_0[i] * fe_0 - 2.0 * ta2_xy_y_yyz_1[i] * fe_0 + 2.0 * ta2_xy_yy_yz_0[i] * fe_0 - 2.0 * ta2_xy_yy_yz_1[i] * fe_0 + ta1_x_yy_yyz_1[i] + ta2_xy_yy_yyz_0[i] * pa_y[i] - ta2_xy_yy_yyz_1[i] * pc_y[i];

        ta2_xy_yyy_yzz_0[i] = 2.0 * ta2_xy_y_yzz_0[i] * fe_0 - 2.0 * ta2_xy_y_yzz_1[i] * fe_0 + ta2_xy_yy_zz_0[i] * fe_0 - ta2_xy_yy_zz_1[i] * fe_0 + ta1_x_yy_yzz_1[i] + ta2_xy_yy_yzz_0[i] * pa_y[i] - ta2_xy_yy_yzz_1[i] * pc_y[i];

        ta2_xy_yyy_zzz_0[i] = 2.0 * ta2_xy_y_zzz_0[i] * fe_0 - 2.0 * ta2_xy_y_zzz_1[i] * fe_0 + ta1_x_yy_zzz_1[i] + ta2_xy_yy_zzz_0[i] * pa_y[i] - ta2_xy_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : FF

    auto ta2_xy_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 170);

    auto ta2_xy_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 171);

    auto ta2_xy_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 172);

    auto ta2_xy_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 173);

    auto ta2_xy_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 174);

    auto ta2_xy_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 175);

    auto ta2_xy_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 176);

    auto ta2_xy_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 177);

    auto ta2_xy_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 178);

    auto ta2_xy_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 179);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yyz_xxx_0, ta2_xy_yyz_xxy_0, ta2_xy_yyz_xxz_0, ta2_xy_yyz_xyy_0, ta2_xy_yyz_xyz_0, ta2_xy_yyz_xzz_0, ta2_xy_yyz_yyy_0, ta2_xy_yyz_yyz_0, ta2_xy_yyz_yzz_0, ta2_xy_yyz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyz_xxx_0[i] = ta2_xy_yy_xxx_0[i] * pa_z[i] - ta2_xy_yy_xxx_1[i] * pc_z[i];

        ta2_xy_yyz_xxy_0[i] = ta2_xy_yy_xxy_0[i] * pa_z[i] - ta2_xy_yy_xxy_1[i] * pc_z[i];

        ta2_xy_yyz_xxz_0[i] = ta2_xy_yy_xx_0[i] * fe_0 - ta2_xy_yy_xx_1[i] * fe_0 + ta2_xy_yy_xxz_0[i] * pa_z[i] - ta2_xy_yy_xxz_1[i] * pc_z[i];

        ta2_xy_yyz_xyy_0[i] = ta2_xy_yy_xyy_0[i] * pa_z[i] - ta2_xy_yy_xyy_1[i] * pc_z[i];

        ta2_xy_yyz_xyz_0[i] = ta2_xy_yy_xy_0[i] * fe_0 - ta2_xy_yy_xy_1[i] * fe_0 + ta2_xy_yy_xyz_0[i] * pa_z[i] - ta2_xy_yy_xyz_1[i] * pc_z[i];

        ta2_xy_yyz_xzz_0[i] = 2.0 * ta2_xy_yy_xz_0[i] * fe_0 - 2.0 * ta2_xy_yy_xz_1[i] * fe_0 + ta2_xy_yy_xzz_0[i] * pa_z[i] - ta2_xy_yy_xzz_1[i] * pc_z[i];

        ta2_xy_yyz_yyy_0[i] = ta2_xy_yy_yyy_0[i] * pa_z[i] - ta2_xy_yy_yyy_1[i] * pc_z[i];

        ta2_xy_yyz_yyz_0[i] = ta2_xy_yy_yy_0[i] * fe_0 - ta2_xy_yy_yy_1[i] * fe_0 + ta2_xy_yy_yyz_0[i] * pa_z[i] - ta2_xy_yy_yyz_1[i] * pc_z[i];

        ta2_xy_yyz_yzz_0[i] = 2.0 * ta2_xy_yy_yz_0[i] * fe_0 - 2.0 * ta2_xy_yy_yz_1[i] * fe_0 + ta2_xy_yy_yzz_0[i] * pa_z[i] - ta2_xy_yy_yzz_1[i] * pc_z[i];

        ta2_xy_yyz_zzz_0[i] = 3.0 * ta2_xy_yy_zz_0[i] * fe_0 - 3.0 * ta2_xy_yy_zz_1[i] * fe_0 + ta2_xy_yy_zzz_0[i] * pa_z[i] - ta2_xy_yy_zzz_1[i] * pc_z[i];
    }

    // Set up 180-190 components of targeted buffer : FF

    auto ta2_xy_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 180);

    auto ta2_xy_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 181);

    auto ta2_xy_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 182);

    auto ta2_xy_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 183);

    auto ta2_xy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 184);

    auto ta2_xy_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 185);

    auto ta2_xy_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 186);

    auto ta2_xy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 187);

    auto ta2_xy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 188);

    auto ta2_xy_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 189);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_zz_xxx_1, ta1_x_zz_xxz_1, ta1_x_zz_xyz_1, ta1_x_zz_xzz_1, ta1_x_zz_yyz_1, ta1_x_zz_yzz_1, ta1_x_zz_zzz_1, ta2_xy_y_xxy_0, ta2_xy_y_xxy_1, ta2_xy_y_xyy_0, ta2_xy_y_xyy_1, ta2_xy_y_yyy_0, ta2_xy_y_yyy_1, ta2_xy_yz_xxy_0, ta2_xy_yz_xxy_1, ta2_xy_yz_xyy_0, ta2_xy_yz_xyy_1, ta2_xy_yz_yyy_0, ta2_xy_yz_yyy_1, ta2_xy_yzz_xxx_0, ta2_xy_yzz_xxy_0, ta2_xy_yzz_xxz_0, ta2_xy_yzz_xyy_0, ta2_xy_yzz_xyz_0, ta2_xy_yzz_xzz_0, ta2_xy_yzz_yyy_0, ta2_xy_yzz_yyz_0, ta2_xy_yzz_yzz_0, ta2_xy_yzz_zzz_0, ta2_xy_zz_xxx_0, ta2_xy_zz_xxx_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzz_xxx_0[i] = ta1_x_zz_xxx_1[i] + ta2_xy_zz_xxx_0[i] * pa_y[i] - ta2_xy_zz_xxx_1[i] * pc_y[i];

        ta2_xy_yzz_xxy_0[i] = ta2_xy_y_xxy_0[i] * fe_0 - ta2_xy_y_xxy_1[i] * fe_0 + ta2_xy_yz_xxy_0[i] * pa_z[i] - ta2_xy_yz_xxy_1[i] * pc_z[i];

        ta2_xy_yzz_xxz_0[i] = ta1_x_zz_xxz_1[i] + ta2_xy_zz_xxz_0[i] * pa_y[i] - ta2_xy_zz_xxz_1[i] * pc_y[i];

        ta2_xy_yzz_xyy_0[i] = ta2_xy_y_xyy_0[i] * fe_0 - ta2_xy_y_xyy_1[i] * fe_0 + ta2_xy_yz_xyy_0[i] * pa_z[i] - ta2_xy_yz_xyy_1[i] * pc_z[i];

        ta2_xy_yzz_xyz_0[i] = ta2_xy_zz_xz_0[i] * fe_0 - ta2_xy_zz_xz_1[i] * fe_0 + ta1_x_zz_xyz_1[i] + ta2_xy_zz_xyz_0[i] * pa_y[i] - ta2_xy_zz_xyz_1[i] * pc_y[i];

        ta2_xy_yzz_xzz_0[i] = ta1_x_zz_xzz_1[i] + ta2_xy_zz_xzz_0[i] * pa_y[i] - ta2_xy_zz_xzz_1[i] * pc_y[i];

        ta2_xy_yzz_yyy_0[i] = ta2_xy_y_yyy_0[i] * fe_0 - ta2_xy_y_yyy_1[i] * fe_0 + ta2_xy_yz_yyy_0[i] * pa_z[i] - ta2_xy_yz_yyy_1[i] * pc_z[i];

        ta2_xy_yzz_yyz_0[i] = 2.0 * ta2_xy_zz_yz_0[i] * fe_0 - 2.0 * ta2_xy_zz_yz_1[i] * fe_0 + ta1_x_zz_yyz_1[i] + ta2_xy_zz_yyz_0[i] * pa_y[i] - ta2_xy_zz_yyz_1[i] * pc_y[i];

        ta2_xy_yzz_yzz_0[i] = ta2_xy_zz_zz_0[i] * fe_0 - ta2_xy_zz_zz_1[i] * fe_0 + ta1_x_zz_yzz_1[i] + ta2_xy_zz_yzz_0[i] * pa_y[i] - ta2_xy_zz_yzz_1[i] * pc_y[i];

        ta2_xy_yzz_zzz_0[i] = ta1_x_zz_zzz_1[i] + ta2_xy_zz_zzz_0[i] * pa_y[i] - ta2_xy_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 190-200 components of targeted buffer : FF

    auto ta2_xy_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 190);

    auto ta2_xy_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 191);

    auto ta2_xy_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 192);

    auto ta2_xy_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 193);

    auto ta2_xy_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 194);

    auto ta2_xy_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 195);

    auto ta2_xy_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 196);

    auto ta2_xy_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 197);

    auto ta2_xy_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 198);

    auto ta2_xy_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 199);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_z_xxx_0, ta2_xy_z_xxx_1, ta2_xy_z_xxy_0, ta2_xy_z_xxy_1, ta2_xy_z_xxz_0, ta2_xy_z_xxz_1, ta2_xy_z_xyy_0, ta2_xy_z_xyy_1, ta2_xy_z_xyz_0, ta2_xy_z_xyz_1, ta2_xy_z_xzz_0, ta2_xy_z_xzz_1, ta2_xy_z_yyy_0, ta2_xy_z_yyy_1, ta2_xy_z_yyz_0, ta2_xy_z_yyz_1, ta2_xy_z_yzz_0, ta2_xy_z_yzz_1, ta2_xy_z_zzz_0, ta2_xy_z_zzz_1, ta2_xy_zz_xx_0, ta2_xy_zz_xx_1, ta2_xy_zz_xxx_0, ta2_xy_zz_xxx_1, ta2_xy_zz_xxy_0, ta2_xy_zz_xxy_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xy_0, ta2_xy_zz_xy_1, ta2_xy_zz_xyy_0, ta2_xy_zz_xyy_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_yy_0, ta2_xy_zz_yy_1, ta2_xy_zz_yyy_0, ta2_xy_zz_yyy_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, ta2_xy_zzz_xxx_0, ta2_xy_zzz_xxy_0, ta2_xy_zzz_xxz_0, ta2_xy_zzz_xyy_0, ta2_xy_zzz_xyz_0, ta2_xy_zzz_xzz_0, ta2_xy_zzz_yyy_0, ta2_xy_zzz_yyz_0, ta2_xy_zzz_yzz_0, ta2_xy_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzz_xxx_0[i] = 2.0 * ta2_xy_z_xxx_0[i] * fe_0 - 2.0 * ta2_xy_z_xxx_1[i] * fe_0 + ta2_xy_zz_xxx_0[i] * pa_z[i] - ta2_xy_zz_xxx_1[i] * pc_z[i];

        ta2_xy_zzz_xxy_0[i] = 2.0 * ta2_xy_z_xxy_0[i] * fe_0 - 2.0 * ta2_xy_z_xxy_1[i] * fe_0 + ta2_xy_zz_xxy_0[i] * pa_z[i] - ta2_xy_zz_xxy_1[i] * pc_z[i];

        ta2_xy_zzz_xxz_0[i] = 2.0 * ta2_xy_z_xxz_0[i] * fe_0 - 2.0 * ta2_xy_z_xxz_1[i] * fe_0 + ta2_xy_zz_xx_0[i] * fe_0 - ta2_xy_zz_xx_1[i] * fe_0 + ta2_xy_zz_xxz_0[i] * pa_z[i] - ta2_xy_zz_xxz_1[i] * pc_z[i];

        ta2_xy_zzz_xyy_0[i] = 2.0 * ta2_xy_z_xyy_0[i] * fe_0 - 2.0 * ta2_xy_z_xyy_1[i] * fe_0 + ta2_xy_zz_xyy_0[i] * pa_z[i] - ta2_xy_zz_xyy_1[i] * pc_z[i];

        ta2_xy_zzz_xyz_0[i] = 2.0 * ta2_xy_z_xyz_0[i] * fe_0 - 2.0 * ta2_xy_z_xyz_1[i] * fe_0 + ta2_xy_zz_xy_0[i] * fe_0 - ta2_xy_zz_xy_1[i] * fe_0 + ta2_xy_zz_xyz_0[i] * pa_z[i] - ta2_xy_zz_xyz_1[i] * pc_z[i];

        ta2_xy_zzz_xzz_0[i] = 2.0 * ta2_xy_z_xzz_0[i] * fe_0 - 2.0 * ta2_xy_z_xzz_1[i] * fe_0 + 2.0 * ta2_xy_zz_xz_0[i] * fe_0 - 2.0 * ta2_xy_zz_xz_1[i] * fe_0 + ta2_xy_zz_xzz_0[i] * pa_z[i] - ta2_xy_zz_xzz_1[i] * pc_z[i];

        ta2_xy_zzz_yyy_0[i] = 2.0 * ta2_xy_z_yyy_0[i] * fe_0 - 2.0 * ta2_xy_z_yyy_1[i] * fe_0 + ta2_xy_zz_yyy_0[i] * pa_z[i] - ta2_xy_zz_yyy_1[i] * pc_z[i];

        ta2_xy_zzz_yyz_0[i] = 2.0 * ta2_xy_z_yyz_0[i] * fe_0 - 2.0 * ta2_xy_z_yyz_1[i] * fe_0 + ta2_xy_zz_yy_0[i] * fe_0 - ta2_xy_zz_yy_1[i] * fe_0 + ta2_xy_zz_yyz_0[i] * pa_z[i] - ta2_xy_zz_yyz_1[i] * pc_z[i];

        ta2_xy_zzz_yzz_0[i] = 2.0 * ta2_xy_z_yzz_0[i] * fe_0 - 2.0 * ta2_xy_z_yzz_1[i] * fe_0 + 2.0 * ta2_xy_zz_yz_0[i] * fe_0 - 2.0 * ta2_xy_zz_yz_1[i] * fe_0 + ta2_xy_zz_yzz_0[i] * pa_z[i] - ta2_xy_zz_yzz_1[i] * pc_z[i];

        ta2_xy_zzz_zzz_0[i] = 2.0 * ta2_xy_z_zzz_0[i] * fe_0 - 2.0 * ta2_xy_z_zzz_1[i] * fe_0 + 3.0 * ta2_xy_zz_zz_0[i] * fe_0 - 3.0 * ta2_xy_zz_zz_1[i] * fe_0 + ta2_xy_zz_zzz_0[i] * pa_z[i] - ta2_xy_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 200-210 components of targeted buffer : FF

    auto ta2_xz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 200);

    auto ta2_xz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 201);

    auto ta2_xz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 202);

    auto ta2_xz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 203);

    auto ta2_xz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 204);

    auto ta2_xz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 205);

    auto ta2_xz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 206);

    auto ta2_xz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 207);

    auto ta2_xz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 208);

    auto ta2_xz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 209);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xx_xxx_1, ta1_z_xx_xxy_1, ta1_z_xx_xxz_1, ta1_z_xx_xyy_1, ta1_z_xx_xyz_1, ta1_z_xx_xzz_1, ta1_z_xx_yyy_1, ta1_z_xx_yyz_1, ta1_z_xx_yzz_1, ta1_z_xx_zzz_1, ta2_xz_x_xxx_0, ta2_xz_x_xxx_1, ta2_xz_x_xxy_0, ta2_xz_x_xxy_1, ta2_xz_x_xxz_0, ta2_xz_x_xxz_1, ta2_xz_x_xyy_0, ta2_xz_x_xyy_1, ta2_xz_x_xyz_0, ta2_xz_x_xyz_1, ta2_xz_x_xzz_0, ta2_xz_x_xzz_1, ta2_xz_x_yyy_0, ta2_xz_x_yyy_1, ta2_xz_x_yyz_0, ta2_xz_x_yyz_1, ta2_xz_x_yzz_0, ta2_xz_x_yzz_1, ta2_xz_x_zzz_0, ta2_xz_x_zzz_1, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_yy_0, ta2_xz_xx_yy_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xx_yyz_0, ta2_xz_xx_yyz_1, ta2_xz_xx_yz_0, ta2_xz_xx_yz_1, ta2_xz_xx_yzz_0, ta2_xz_xx_yzz_1, ta2_xz_xx_zz_0, ta2_xz_xx_zz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xxx_xxx_0, ta2_xz_xxx_xxy_0, ta2_xz_xxx_xxz_0, ta2_xz_xxx_xyy_0, ta2_xz_xxx_xyz_0, ta2_xz_xxx_xzz_0, ta2_xz_xxx_yyy_0, ta2_xz_xxx_yyz_0, ta2_xz_xxx_yzz_0, ta2_xz_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxx_xxx_0[i] = 2.0 * ta2_xz_x_xxx_0[i] * fe_0 - 2.0 * ta2_xz_x_xxx_1[i] * fe_0 + 3.0 * ta2_xz_xx_xx_0[i] * fe_0 - 3.0 * ta2_xz_xx_xx_1[i] * fe_0 + ta1_z_xx_xxx_1[i] + ta2_xz_xx_xxx_0[i] * pa_x[i] - ta2_xz_xx_xxx_1[i] * pc_x[i];

        ta2_xz_xxx_xxy_0[i] = 2.0 * ta2_xz_x_xxy_0[i] * fe_0 - 2.0 * ta2_xz_x_xxy_1[i] * fe_0 + 2.0 * ta2_xz_xx_xy_0[i] * fe_0 - 2.0 * ta2_xz_xx_xy_1[i] * fe_0 + ta1_z_xx_xxy_1[i] + ta2_xz_xx_xxy_0[i] * pa_x[i] - ta2_xz_xx_xxy_1[i] * pc_x[i];

        ta2_xz_xxx_xxz_0[i] = 2.0 * ta2_xz_x_xxz_0[i] * fe_0 - 2.0 * ta2_xz_x_xxz_1[i] * fe_0 + 2.0 * ta2_xz_xx_xz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xz_1[i] * fe_0 + ta1_z_xx_xxz_1[i] + ta2_xz_xx_xxz_0[i] * pa_x[i] - ta2_xz_xx_xxz_1[i] * pc_x[i];

        ta2_xz_xxx_xyy_0[i] = 2.0 * ta2_xz_x_xyy_0[i] * fe_0 - 2.0 * ta2_xz_x_xyy_1[i] * fe_0 + ta2_xz_xx_yy_0[i] * fe_0 - ta2_xz_xx_yy_1[i] * fe_0 + ta1_z_xx_xyy_1[i] + ta2_xz_xx_xyy_0[i] * pa_x[i] - ta2_xz_xx_xyy_1[i] * pc_x[i];

        ta2_xz_xxx_xyz_0[i] = 2.0 * ta2_xz_x_xyz_0[i] * fe_0 - 2.0 * ta2_xz_x_xyz_1[i] * fe_0 + ta2_xz_xx_yz_0[i] * fe_0 - ta2_xz_xx_yz_1[i] * fe_0 + ta1_z_xx_xyz_1[i] + ta2_xz_xx_xyz_0[i] * pa_x[i] - ta2_xz_xx_xyz_1[i] * pc_x[i];

        ta2_xz_xxx_xzz_0[i] = 2.0 * ta2_xz_x_xzz_0[i] * fe_0 - 2.0 * ta2_xz_x_xzz_1[i] * fe_0 + ta2_xz_xx_zz_0[i] * fe_0 - ta2_xz_xx_zz_1[i] * fe_0 + ta1_z_xx_xzz_1[i] + ta2_xz_xx_xzz_0[i] * pa_x[i] - ta2_xz_xx_xzz_1[i] * pc_x[i];

        ta2_xz_xxx_yyy_0[i] = 2.0 * ta2_xz_x_yyy_0[i] * fe_0 - 2.0 * ta2_xz_x_yyy_1[i] * fe_0 + ta1_z_xx_yyy_1[i] + ta2_xz_xx_yyy_0[i] * pa_x[i] - ta2_xz_xx_yyy_1[i] * pc_x[i];

        ta2_xz_xxx_yyz_0[i] = 2.0 * ta2_xz_x_yyz_0[i] * fe_0 - 2.0 * ta2_xz_x_yyz_1[i] * fe_0 + ta1_z_xx_yyz_1[i] + ta2_xz_xx_yyz_0[i] * pa_x[i] - ta2_xz_xx_yyz_1[i] * pc_x[i];

        ta2_xz_xxx_yzz_0[i] = 2.0 * ta2_xz_x_yzz_0[i] * fe_0 - 2.0 * ta2_xz_x_yzz_1[i] * fe_0 + ta1_z_xx_yzz_1[i] + ta2_xz_xx_yzz_0[i] * pa_x[i] - ta2_xz_xx_yzz_1[i] * pc_x[i];

        ta2_xz_xxx_zzz_0[i] = 2.0 * ta2_xz_x_zzz_0[i] * fe_0 - 2.0 * ta2_xz_x_zzz_1[i] * fe_0 + ta1_z_xx_zzz_1[i] + ta2_xz_xx_zzz_0[i] * pa_x[i] - ta2_xz_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : FF

    auto ta2_xz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 210);

    auto ta2_xz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 211);

    auto ta2_xz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 212);

    auto ta2_xz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 213);

    auto ta2_xz_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 214);

    auto ta2_xz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 215);

    auto ta2_xz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 216);

    auto ta2_xz_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 217);

    auto ta2_xz_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 218);

    auto ta2_xz_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 219);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_yy_0, ta2_xz_xx_yy_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xx_yyz_0, ta2_xz_xx_yyz_1, ta2_xz_xx_yz_0, ta2_xz_xx_yz_1, ta2_xz_xx_yzz_0, ta2_xz_xx_yzz_1, ta2_xz_xx_zz_0, ta2_xz_xx_zz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xxy_xxx_0, ta2_xz_xxy_xxy_0, ta2_xz_xxy_xxz_0, ta2_xz_xxy_xyy_0, ta2_xz_xxy_xyz_0, ta2_xz_xxy_xzz_0, ta2_xz_xxy_yyy_0, ta2_xz_xxy_yyz_0, ta2_xz_xxy_yzz_0, ta2_xz_xxy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxy_xxx_0[i] = ta2_xz_xx_xxx_0[i] * pa_y[i] - ta2_xz_xx_xxx_1[i] * pc_y[i];

        ta2_xz_xxy_xxy_0[i] = ta2_xz_xx_xx_0[i] * fe_0 - ta2_xz_xx_xx_1[i] * fe_0 + ta2_xz_xx_xxy_0[i] * pa_y[i] - ta2_xz_xx_xxy_1[i] * pc_y[i];

        ta2_xz_xxy_xxz_0[i] = ta2_xz_xx_xxz_0[i] * pa_y[i] - ta2_xz_xx_xxz_1[i] * pc_y[i];

        ta2_xz_xxy_xyy_0[i] = 2.0 * ta2_xz_xx_xy_0[i] * fe_0 - 2.0 * ta2_xz_xx_xy_1[i] * fe_0 + ta2_xz_xx_xyy_0[i] * pa_y[i] - ta2_xz_xx_xyy_1[i] * pc_y[i];

        ta2_xz_xxy_xyz_0[i] = ta2_xz_xx_xz_0[i] * fe_0 - ta2_xz_xx_xz_1[i] * fe_0 + ta2_xz_xx_xyz_0[i] * pa_y[i] - ta2_xz_xx_xyz_1[i] * pc_y[i];

        ta2_xz_xxy_xzz_0[i] = ta2_xz_xx_xzz_0[i] * pa_y[i] - ta2_xz_xx_xzz_1[i] * pc_y[i];

        ta2_xz_xxy_yyy_0[i] = 3.0 * ta2_xz_xx_yy_0[i] * fe_0 - 3.0 * ta2_xz_xx_yy_1[i] * fe_0 + ta2_xz_xx_yyy_0[i] * pa_y[i] - ta2_xz_xx_yyy_1[i] * pc_y[i];

        ta2_xz_xxy_yyz_0[i] = 2.0 * ta2_xz_xx_yz_0[i] * fe_0 - 2.0 * ta2_xz_xx_yz_1[i] * fe_0 + ta2_xz_xx_yyz_0[i] * pa_y[i] - ta2_xz_xx_yyz_1[i] * pc_y[i];

        ta2_xz_xxy_yzz_0[i] = ta2_xz_xx_zz_0[i] * fe_0 - ta2_xz_xx_zz_1[i] * fe_0 + ta2_xz_xx_yzz_0[i] * pa_y[i] - ta2_xz_xx_yzz_1[i] * pc_y[i];

        ta2_xz_xxy_zzz_0[i] = ta2_xz_xx_zzz_0[i] * pa_y[i] - ta2_xz_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 220-230 components of targeted buffer : FF

    auto ta2_xz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 220);

    auto ta2_xz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 221);

    auto ta2_xz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 222);

    auto ta2_xz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 223);

    auto ta2_xz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 224);

    auto ta2_xz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 225);

    auto ta2_xz_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 226);

    auto ta2_xz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 227);

    auto ta2_xz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 228);

    auto ta2_xz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 229);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xx_xxx_1, ta1_x_xx_xxy_1, ta1_x_xx_xxz_1, ta1_x_xx_xyy_1, ta1_x_xx_xyz_1, ta1_x_xx_xzz_1, ta1_x_xx_yyy_1, ta1_z_xz_yyz_1, ta1_z_xz_yzz_1, ta1_z_xz_zzz_1, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xxz_xxx_0, ta2_xz_xxz_xxy_0, ta2_xz_xxz_xxz_0, ta2_xz_xxz_xyy_0, ta2_xz_xxz_xyz_0, ta2_xz_xxz_xzz_0, ta2_xz_xxz_yyy_0, ta2_xz_xxz_yyz_0, ta2_xz_xxz_yzz_0, ta2_xz_xxz_zzz_0, ta2_xz_xz_yyz_0, ta2_xz_xz_yyz_1, ta2_xz_xz_yzz_0, ta2_xz_xz_yzz_1, ta2_xz_xz_zzz_0, ta2_xz_xz_zzz_1, ta2_xz_z_yyz_0, ta2_xz_z_yyz_1, ta2_xz_z_yzz_0, ta2_xz_z_yzz_1, ta2_xz_z_zzz_0, ta2_xz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxz_xxx_0[i] = ta1_x_xx_xxx_1[i] + ta2_xz_xx_xxx_0[i] * pa_z[i] - ta2_xz_xx_xxx_1[i] * pc_z[i];

        ta2_xz_xxz_xxy_0[i] = ta1_x_xx_xxy_1[i] + ta2_xz_xx_xxy_0[i] * pa_z[i] - ta2_xz_xx_xxy_1[i] * pc_z[i];

        ta2_xz_xxz_xxz_0[i] = ta2_xz_xx_xx_0[i] * fe_0 - ta2_xz_xx_xx_1[i] * fe_0 + ta1_x_xx_xxz_1[i] + ta2_xz_xx_xxz_0[i] * pa_z[i] - ta2_xz_xx_xxz_1[i] * pc_z[i];

        ta2_xz_xxz_xyy_0[i] = ta1_x_xx_xyy_1[i] + ta2_xz_xx_xyy_0[i] * pa_z[i] - ta2_xz_xx_xyy_1[i] * pc_z[i];

        ta2_xz_xxz_xyz_0[i] = ta2_xz_xx_xy_0[i] * fe_0 - ta2_xz_xx_xy_1[i] * fe_0 + ta1_x_xx_xyz_1[i] + ta2_xz_xx_xyz_0[i] * pa_z[i] - ta2_xz_xx_xyz_1[i] * pc_z[i];

        ta2_xz_xxz_xzz_0[i] = 2.0 * ta2_xz_xx_xz_0[i] * fe_0 - 2.0 * ta2_xz_xx_xz_1[i] * fe_0 + ta1_x_xx_xzz_1[i] + ta2_xz_xx_xzz_0[i] * pa_z[i] - ta2_xz_xx_xzz_1[i] * pc_z[i];

        ta2_xz_xxz_yyy_0[i] = ta1_x_xx_yyy_1[i] + ta2_xz_xx_yyy_0[i] * pa_z[i] - ta2_xz_xx_yyy_1[i] * pc_z[i];

        ta2_xz_xxz_yyz_0[i] = ta2_xz_z_yyz_0[i] * fe_0 - ta2_xz_z_yyz_1[i] * fe_0 + ta1_z_xz_yyz_1[i] + ta2_xz_xz_yyz_0[i] * pa_x[i] - ta2_xz_xz_yyz_1[i] * pc_x[i];

        ta2_xz_xxz_yzz_0[i] = ta2_xz_z_yzz_0[i] * fe_0 - ta2_xz_z_yzz_1[i] * fe_0 + ta1_z_xz_yzz_1[i] + ta2_xz_xz_yzz_0[i] * pa_x[i] - ta2_xz_xz_yzz_1[i] * pc_x[i];

        ta2_xz_xxz_zzz_0[i] = ta2_xz_z_zzz_0[i] * fe_0 - ta2_xz_z_zzz_1[i] * fe_0 + ta1_z_xz_zzz_1[i] + ta2_xz_xz_zzz_0[i] * pa_x[i] - ta2_xz_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 230-240 components of targeted buffer : FF

    auto ta2_xz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 230);

    auto ta2_xz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 231);

    auto ta2_xz_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 232);

    auto ta2_xz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 233);

    auto ta2_xz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 234);

    auto ta2_xz_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 235);

    auto ta2_xz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 236);

    auto ta2_xz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 237);

    auto ta2_xz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 238);

    auto ta2_xz_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 239);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_yy_xxy_1, ta1_z_yy_xyy_1, ta1_z_yy_xyz_1, ta1_z_yy_yyy_1, ta1_z_yy_yyz_1, ta1_z_yy_yzz_1, ta1_z_yy_zzz_1, ta2_xz_x_xxx_0, ta2_xz_x_xxx_1, ta2_xz_x_xxz_0, ta2_xz_x_xxz_1, ta2_xz_x_xzz_0, ta2_xz_x_xzz_1, ta2_xz_xy_xxx_0, ta2_xz_xy_xxx_1, ta2_xz_xy_xxz_0, ta2_xz_xy_xxz_1, ta2_xz_xy_xzz_0, ta2_xz_xy_xzz_1, ta2_xz_xyy_xxx_0, ta2_xz_xyy_xxy_0, ta2_xz_xyy_xxz_0, ta2_xz_xyy_xyy_0, ta2_xz_xyy_xyz_0, ta2_xz_xyy_xzz_0, ta2_xz_xyy_yyy_0, ta2_xz_xyy_yyz_0, ta2_xz_xyy_yzz_0, ta2_xz_xyy_zzz_0, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_zzz_0, ta2_xz_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyy_xxx_0[i] = ta2_xz_x_xxx_0[i] * fe_0 - ta2_xz_x_xxx_1[i] * fe_0 + ta2_xz_xy_xxx_0[i] * pa_y[i] - ta2_xz_xy_xxx_1[i] * pc_y[i];

        ta2_xz_xyy_xxy_0[i] = 2.0 * ta2_xz_yy_xy_0[i] * fe_0 - 2.0 * ta2_xz_yy_xy_1[i] * fe_0 + ta1_z_yy_xxy_1[i] + ta2_xz_yy_xxy_0[i] * pa_x[i] - ta2_xz_yy_xxy_1[i] * pc_x[i];

        ta2_xz_xyy_xxz_0[i] = ta2_xz_x_xxz_0[i] * fe_0 - ta2_xz_x_xxz_1[i] * fe_0 + ta2_xz_xy_xxz_0[i] * pa_y[i] - ta2_xz_xy_xxz_1[i] * pc_y[i];

        ta2_xz_xyy_xyy_0[i] = ta2_xz_yy_yy_0[i] * fe_0 - ta2_xz_yy_yy_1[i] * fe_0 + ta1_z_yy_xyy_1[i] + ta2_xz_yy_xyy_0[i] * pa_x[i] - ta2_xz_yy_xyy_1[i] * pc_x[i];

        ta2_xz_xyy_xyz_0[i] = ta2_xz_yy_yz_0[i] * fe_0 - ta2_xz_yy_yz_1[i] * fe_0 + ta1_z_yy_xyz_1[i] + ta2_xz_yy_xyz_0[i] * pa_x[i] - ta2_xz_yy_xyz_1[i] * pc_x[i];

        ta2_xz_xyy_xzz_0[i] = ta2_xz_x_xzz_0[i] * fe_0 - ta2_xz_x_xzz_1[i] * fe_0 + ta2_xz_xy_xzz_0[i] * pa_y[i] - ta2_xz_xy_xzz_1[i] * pc_y[i];

        ta2_xz_xyy_yyy_0[i] = ta1_z_yy_yyy_1[i] + ta2_xz_yy_yyy_0[i] * pa_x[i] - ta2_xz_yy_yyy_1[i] * pc_x[i];

        ta2_xz_xyy_yyz_0[i] = ta1_z_yy_yyz_1[i] + ta2_xz_yy_yyz_0[i] * pa_x[i] - ta2_xz_yy_yyz_1[i] * pc_x[i];

        ta2_xz_xyy_yzz_0[i] = ta1_z_yy_yzz_1[i] + ta2_xz_yy_yzz_0[i] * pa_x[i] - ta2_xz_yy_yzz_1[i] * pc_x[i];

        ta2_xz_xyy_zzz_0[i] = ta1_z_yy_zzz_1[i] + ta2_xz_yy_zzz_0[i] * pa_x[i] - ta2_xz_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 240-250 components of targeted buffer : FF

    auto ta2_xz_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 240);

    auto ta2_xz_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 241);

    auto ta2_xz_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 242);

    auto ta2_xz_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 243);

    auto ta2_xz_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 244);

    auto ta2_xz_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 245);

    auto ta2_xz_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 246);

    auto ta2_xz_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 247);

    auto ta2_xz_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 248);

    auto ta2_xz_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 249);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xy_xxy_1, ta1_x_xy_xyy_1, ta1_z_yz_yyy_1, ta1_z_yz_yyz_1, ta1_z_yz_yzz_1, ta2_xz_xy_xxy_0, ta2_xz_xy_xxy_1, ta2_xz_xy_xyy_0, ta2_xz_xy_xyy_1, ta2_xz_xyz_xxx_0, ta2_xz_xyz_xxy_0, ta2_xz_xyz_xxz_0, ta2_xz_xyz_xyy_0, ta2_xz_xyz_xyz_0, ta2_xz_xyz_xzz_0, ta2_xz_xyz_yyy_0, ta2_xz_xyz_yyz_0, ta2_xz_xyz_yzz_0, ta2_xz_xyz_zzz_0, ta2_xz_xz_xxx_0, ta2_xz_xz_xxx_1, ta2_xz_xz_xxz_0, ta2_xz_xz_xxz_1, ta2_xz_xz_xyz_0, ta2_xz_xz_xyz_1, ta2_xz_xz_xz_0, ta2_xz_xz_xz_1, ta2_xz_xz_xzz_0, ta2_xz_xz_xzz_1, ta2_xz_xz_zzz_0, ta2_xz_xz_zzz_1, ta2_xz_yz_yyy_0, ta2_xz_yz_yyy_1, ta2_xz_yz_yyz_0, ta2_xz_yz_yyz_1, ta2_xz_yz_yzz_0, ta2_xz_yz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyz_xxx_0[i] = ta2_xz_xz_xxx_0[i] * pa_y[i] - ta2_xz_xz_xxx_1[i] * pc_y[i];

        ta2_xz_xyz_xxy_0[i] = ta1_x_xy_xxy_1[i] + ta2_xz_xy_xxy_0[i] * pa_z[i] - ta2_xz_xy_xxy_1[i] * pc_z[i];

        ta2_xz_xyz_xxz_0[i] = ta2_xz_xz_xxz_0[i] * pa_y[i] - ta2_xz_xz_xxz_1[i] * pc_y[i];

        ta2_xz_xyz_xyy_0[i] = ta1_x_xy_xyy_1[i] + ta2_xz_xy_xyy_0[i] * pa_z[i] - ta2_xz_xy_xyy_1[i] * pc_z[i];

        ta2_xz_xyz_xyz_0[i] = ta2_xz_xz_xz_0[i] * fe_0 - ta2_xz_xz_xz_1[i] * fe_0 + ta2_xz_xz_xyz_0[i] * pa_y[i] - ta2_xz_xz_xyz_1[i] * pc_y[i];

        ta2_xz_xyz_xzz_0[i] = ta2_xz_xz_xzz_0[i] * pa_y[i] - ta2_xz_xz_xzz_1[i] * pc_y[i];

        ta2_xz_xyz_yyy_0[i] = ta1_z_yz_yyy_1[i] + ta2_xz_yz_yyy_0[i] * pa_x[i] - ta2_xz_yz_yyy_1[i] * pc_x[i];

        ta2_xz_xyz_yyz_0[i] = ta1_z_yz_yyz_1[i] + ta2_xz_yz_yyz_0[i] * pa_x[i] - ta2_xz_yz_yyz_1[i] * pc_x[i];

        ta2_xz_xyz_yzz_0[i] = ta1_z_yz_yzz_1[i] + ta2_xz_yz_yzz_0[i] * pa_x[i] - ta2_xz_yz_yzz_1[i] * pc_x[i];

        ta2_xz_xyz_zzz_0[i] = ta2_xz_xz_zzz_0[i] * pa_y[i] - ta2_xz_xz_zzz_1[i] * pc_y[i];
    }

    // Set up 250-260 components of targeted buffer : FF

    auto ta2_xz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 250);

    auto ta2_xz_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 251);

    auto ta2_xz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 252);

    auto ta2_xz_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 253);

    auto ta2_xz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 254);

    auto ta2_xz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 255);

    auto ta2_xz_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 256);

    auto ta2_xz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 257);

    auto ta2_xz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 258);

    auto ta2_xz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 259);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_zz_xxx_1, ta1_z_zz_xxy_1, ta1_z_zz_xxz_1, ta1_z_zz_xyy_1, ta1_z_zz_xyz_1, ta1_z_zz_xzz_1, ta1_z_zz_yyy_1, ta1_z_zz_yyz_1, ta1_z_zz_yzz_1, ta1_z_zz_zzz_1, ta2_xz_xzz_xxx_0, ta2_xz_xzz_xxy_0, ta2_xz_xzz_xxz_0, ta2_xz_xzz_xyy_0, ta2_xz_xzz_xyz_0, ta2_xz_xzz_xzz_0, ta2_xz_xzz_yyy_0, ta2_xz_xzz_yyz_0, ta2_xz_xzz_yzz_0, ta2_xz_xzz_zzz_0, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzz_xxx_0[i] = 3.0 * ta2_xz_zz_xx_0[i] * fe_0 - 3.0 * ta2_xz_zz_xx_1[i] * fe_0 + ta1_z_zz_xxx_1[i] + ta2_xz_zz_xxx_0[i] * pa_x[i] - ta2_xz_zz_xxx_1[i] * pc_x[i];

        ta2_xz_xzz_xxy_0[i] = 2.0 * ta2_xz_zz_xy_0[i] * fe_0 - 2.0 * ta2_xz_zz_xy_1[i] * fe_0 + ta1_z_zz_xxy_1[i] + ta2_xz_zz_xxy_0[i] * pa_x[i] - ta2_xz_zz_xxy_1[i] * pc_x[i];

        ta2_xz_xzz_xxz_0[i] = 2.0 * ta2_xz_zz_xz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xz_1[i] * fe_0 + ta1_z_zz_xxz_1[i] + ta2_xz_zz_xxz_0[i] * pa_x[i] - ta2_xz_zz_xxz_1[i] * pc_x[i];

        ta2_xz_xzz_xyy_0[i] = ta2_xz_zz_yy_0[i] * fe_0 - ta2_xz_zz_yy_1[i] * fe_0 + ta1_z_zz_xyy_1[i] + ta2_xz_zz_xyy_0[i] * pa_x[i] - ta2_xz_zz_xyy_1[i] * pc_x[i];

        ta2_xz_xzz_xyz_0[i] = ta2_xz_zz_yz_0[i] * fe_0 - ta2_xz_zz_yz_1[i] * fe_0 + ta1_z_zz_xyz_1[i] + ta2_xz_zz_xyz_0[i] * pa_x[i] - ta2_xz_zz_xyz_1[i] * pc_x[i];

        ta2_xz_xzz_xzz_0[i] = ta2_xz_zz_zz_0[i] * fe_0 - ta2_xz_zz_zz_1[i] * fe_0 + ta1_z_zz_xzz_1[i] + ta2_xz_zz_xzz_0[i] * pa_x[i] - ta2_xz_zz_xzz_1[i] * pc_x[i];

        ta2_xz_xzz_yyy_0[i] = ta1_z_zz_yyy_1[i] + ta2_xz_zz_yyy_0[i] * pa_x[i] - ta2_xz_zz_yyy_1[i] * pc_x[i];

        ta2_xz_xzz_yyz_0[i] = ta1_z_zz_yyz_1[i] + ta2_xz_zz_yyz_0[i] * pa_x[i] - ta2_xz_zz_yyz_1[i] * pc_x[i];

        ta2_xz_xzz_yzz_0[i] = ta1_z_zz_yzz_1[i] + ta2_xz_zz_yzz_0[i] * pa_x[i] - ta2_xz_zz_yzz_1[i] * pc_x[i];

        ta2_xz_xzz_zzz_0[i] = ta1_z_zz_zzz_1[i] + ta2_xz_zz_zzz_0[i] * pa_x[i] - ta2_xz_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 260-270 components of targeted buffer : FF

    auto ta2_xz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 260);

    auto ta2_xz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 261);

    auto ta2_xz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 262);

    auto ta2_xz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 263);

    auto ta2_xz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 264);

    auto ta2_xz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 265);

    auto ta2_xz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 266);

    auto ta2_xz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 267);

    auto ta2_xz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 268);

    auto ta2_xz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 269);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_y_xxx_0, ta2_xz_y_xxx_1, ta2_xz_y_xxy_0, ta2_xz_y_xxy_1, ta2_xz_y_xxz_0, ta2_xz_y_xxz_1, ta2_xz_y_xyy_0, ta2_xz_y_xyy_1, ta2_xz_y_xyz_0, ta2_xz_y_xyz_1, ta2_xz_y_xzz_0, ta2_xz_y_xzz_1, ta2_xz_y_yyy_0, ta2_xz_y_yyy_1, ta2_xz_y_yyz_0, ta2_xz_y_yyz_1, ta2_xz_y_yzz_0, ta2_xz_y_yzz_1, ta2_xz_y_zzz_0, ta2_xz_y_zzz_1, ta2_xz_yy_xx_0, ta2_xz_yy_xx_1, ta2_xz_yy_xxx_0, ta2_xz_yy_xxx_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xxz_0, ta2_xz_yy_xxz_1, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_xz_0, ta2_xz_yy_xz_1, ta2_xz_yy_xzz_0, ta2_xz_yy_xzz_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_zz_0, ta2_xz_yy_zz_1, ta2_xz_yy_zzz_0, ta2_xz_yy_zzz_1, ta2_xz_yyy_xxx_0, ta2_xz_yyy_xxy_0, ta2_xz_yyy_xxz_0, ta2_xz_yyy_xyy_0, ta2_xz_yyy_xyz_0, ta2_xz_yyy_xzz_0, ta2_xz_yyy_yyy_0, ta2_xz_yyy_yyz_0, ta2_xz_yyy_yzz_0, ta2_xz_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyy_xxx_0[i] = 2.0 * ta2_xz_y_xxx_0[i] * fe_0 - 2.0 * ta2_xz_y_xxx_1[i] * fe_0 + ta2_xz_yy_xxx_0[i] * pa_y[i] - ta2_xz_yy_xxx_1[i] * pc_y[i];

        ta2_xz_yyy_xxy_0[i] = 2.0 * ta2_xz_y_xxy_0[i] * fe_0 - 2.0 * ta2_xz_y_xxy_1[i] * fe_0 + ta2_xz_yy_xx_0[i] * fe_0 - ta2_xz_yy_xx_1[i] * fe_0 + ta2_xz_yy_xxy_0[i] * pa_y[i] - ta2_xz_yy_xxy_1[i] * pc_y[i];

        ta2_xz_yyy_xxz_0[i] = 2.0 * ta2_xz_y_xxz_0[i] * fe_0 - 2.0 * ta2_xz_y_xxz_1[i] * fe_0 + ta2_xz_yy_xxz_0[i] * pa_y[i] - ta2_xz_yy_xxz_1[i] * pc_y[i];

        ta2_xz_yyy_xyy_0[i] = 2.0 * ta2_xz_y_xyy_0[i] * fe_0 - 2.0 * ta2_xz_y_xyy_1[i] * fe_0 + 2.0 * ta2_xz_yy_xy_0[i] * fe_0 - 2.0 * ta2_xz_yy_xy_1[i] * fe_0 + ta2_xz_yy_xyy_0[i] * pa_y[i] - ta2_xz_yy_xyy_1[i] * pc_y[i];

        ta2_xz_yyy_xyz_0[i] = 2.0 * ta2_xz_y_xyz_0[i] * fe_0 - 2.0 * ta2_xz_y_xyz_1[i] * fe_0 + ta2_xz_yy_xz_0[i] * fe_0 - ta2_xz_yy_xz_1[i] * fe_0 + ta2_xz_yy_xyz_0[i] * pa_y[i] - ta2_xz_yy_xyz_1[i] * pc_y[i];

        ta2_xz_yyy_xzz_0[i] = 2.0 * ta2_xz_y_xzz_0[i] * fe_0 - 2.0 * ta2_xz_y_xzz_1[i] * fe_0 + ta2_xz_yy_xzz_0[i] * pa_y[i] - ta2_xz_yy_xzz_1[i] * pc_y[i];

        ta2_xz_yyy_yyy_0[i] = 2.0 * ta2_xz_y_yyy_0[i] * fe_0 - 2.0 * ta2_xz_y_yyy_1[i] * fe_0 + 3.0 * ta2_xz_yy_yy_0[i] * fe_0 - 3.0 * ta2_xz_yy_yy_1[i] * fe_0 + ta2_xz_yy_yyy_0[i] * pa_y[i] - ta2_xz_yy_yyy_1[i] * pc_y[i];

        ta2_xz_yyy_yyz_0[i] = 2.0 * ta2_xz_y_yyz_0[i] * fe_0 - 2.0 * ta2_xz_y_yyz_1[i] * fe_0 + 2.0 * ta2_xz_yy_yz_0[i] * fe_0 - 2.0 * ta2_xz_yy_yz_1[i] * fe_0 + ta2_xz_yy_yyz_0[i] * pa_y[i] - ta2_xz_yy_yyz_1[i] * pc_y[i];

        ta2_xz_yyy_yzz_0[i] = 2.0 * ta2_xz_y_yzz_0[i] * fe_0 - 2.0 * ta2_xz_y_yzz_1[i] * fe_0 + ta2_xz_yy_zz_0[i] * fe_0 - ta2_xz_yy_zz_1[i] * fe_0 + ta2_xz_yy_yzz_0[i] * pa_y[i] - ta2_xz_yy_yzz_1[i] * pc_y[i];

        ta2_xz_yyy_zzz_0[i] = 2.0 * ta2_xz_y_zzz_0[i] * fe_0 - 2.0 * ta2_xz_y_zzz_1[i] * fe_0 + ta2_xz_yy_zzz_0[i] * pa_y[i] - ta2_xz_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 270-280 components of targeted buffer : FF

    auto ta2_xz_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 270);

    auto ta2_xz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 271);

    auto ta2_xz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 272);

    auto ta2_xz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 273);

    auto ta2_xz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 274);

    auto ta2_xz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 275);

    auto ta2_xz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 276);

    auto ta2_xz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 277);

    auto ta2_xz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 278);

    auto ta2_xz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 279);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yy_xxx_1, ta1_x_yy_xxy_1, ta1_x_yy_xyy_1, ta1_x_yy_xyz_1, ta1_x_yy_yyy_1, ta1_x_yy_yyz_1, ta1_x_yy_yzz_1, ta2_xz_yy_xxx_0, ta2_xz_yy_xxx_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yyz_xxx_0, ta2_xz_yyz_xxy_0, ta2_xz_yyz_xxz_0, ta2_xz_yyz_xyy_0, ta2_xz_yyz_xyz_0, ta2_xz_yyz_xzz_0, ta2_xz_yyz_yyy_0, ta2_xz_yyz_yyz_0, ta2_xz_yyz_yzz_0, ta2_xz_yyz_zzz_0, ta2_xz_yz_xxz_0, ta2_xz_yz_xxz_1, ta2_xz_yz_xzz_0, ta2_xz_yz_xzz_1, ta2_xz_yz_zzz_0, ta2_xz_yz_zzz_1, ta2_xz_z_xxz_0, ta2_xz_z_xxz_1, ta2_xz_z_xzz_0, ta2_xz_z_xzz_1, ta2_xz_z_zzz_0, ta2_xz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyz_xxx_0[i] = ta1_x_yy_xxx_1[i] + ta2_xz_yy_xxx_0[i] * pa_z[i] - ta2_xz_yy_xxx_1[i] * pc_z[i];

        ta2_xz_yyz_xxy_0[i] = ta1_x_yy_xxy_1[i] + ta2_xz_yy_xxy_0[i] * pa_z[i] - ta2_xz_yy_xxy_1[i] * pc_z[i];

        ta2_xz_yyz_xxz_0[i] = ta2_xz_z_xxz_0[i] * fe_0 - ta2_xz_z_xxz_1[i] * fe_0 + ta2_xz_yz_xxz_0[i] * pa_y[i] - ta2_xz_yz_xxz_1[i] * pc_y[i];

        ta2_xz_yyz_xyy_0[i] = ta1_x_yy_xyy_1[i] + ta2_xz_yy_xyy_0[i] * pa_z[i] - ta2_xz_yy_xyy_1[i] * pc_z[i];

        ta2_xz_yyz_xyz_0[i] = ta2_xz_yy_xy_0[i] * fe_0 - ta2_xz_yy_xy_1[i] * fe_0 + ta1_x_yy_xyz_1[i] + ta2_xz_yy_xyz_0[i] * pa_z[i] - ta2_xz_yy_xyz_1[i] * pc_z[i];

        ta2_xz_yyz_xzz_0[i] = ta2_xz_z_xzz_0[i] * fe_0 - ta2_xz_z_xzz_1[i] * fe_0 + ta2_xz_yz_xzz_0[i] * pa_y[i] - ta2_xz_yz_xzz_1[i] * pc_y[i];

        ta2_xz_yyz_yyy_0[i] = ta1_x_yy_yyy_1[i] + ta2_xz_yy_yyy_0[i] * pa_z[i] - ta2_xz_yy_yyy_1[i] * pc_z[i];

        ta2_xz_yyz_yyz_0[i] = ta2_xz_yy_yy_0[i] * fe_0 - ta2_xz_yy_yy_1[i] * fe_0 + ta1_x_yy_yyz_1[i] + ta2_xz_yy_yyz_0[i] * pa_z[i] - ta2_xz_yy_yyz_1[i] * pc_z[i];

        ta2_xz_yyz_yzz_0[i] = 2.0 * ta2_xz_yy_yz_0[i] * fe_0 - 2.0 * ta2_xz_yy_yz_1[i] * fe_0 + ta1_x_yy_yzz_1[i] + ta2_xz_yy_yzz_0[i] * pa_z[i] - ta2_xz_yy_yzz_1[i] * pc_z[i];

        ta2_xz_yyz_zzz_0[i] = ta2_xz_z_zzz_0[i] * fe_0 - ta2_xz_z_zzz_1[i] * fe_0 + ta2_xz_yz_zzz_0[i] * pa_y[i] - ta2_xz_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 280-290 components of targeted buffer : FF

    auto ta2_xz_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 280);

    auto ta2_xz_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 281);

    auto ta2_xz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 282);

    auto ta2_xz_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 283);

    auto ta2_xz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 284);

    auto ta2_xz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 285);

    auto ta2_xz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 286);

    auto ta2_xz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 287);

    auto ta2_xz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 288);

    auto ta2_xz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 289);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_yzz_xxx_0, ta2_xz_yzz_xxy_0, ta2_xz_yzz_xxz_0, ta2_xz_yzz_xyy_0, ta2_xz_yzz_xyz_0, ta2_xz_yzz_xzz_0, ta2_xz_yzz_yyy_0, ta2_xz_yzz_yyz_0, ta2_xz_yzz_yzz_0, ta2_xz_yzz_zzz_0, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzz_xxx_0[i] = ta2_xz_zz_xxx_0[i] * pa_y[i] - ta2_xz_zz_xxx_1[i] * pc_y[i];

        ta2_xz_yzz_xxy_0[i] = ta2_xz_zz_xx_0[i] * fe_0 - ta2_xz_zz_xx_1[i] * fe_0 + ta2_xz_zz_xxy_0[i] * pa_y[i] - ta2_xz_zz_xxy_1[i] * pc_y[i];

        ta2_xz_yzz_xxz_0[i] = ta2_xz_zz_xxz_0[i] * pa_y[i] - ta2_xz_zz_xxz_1[i] * pc_y[i];

        ta2_xz_yzz_xyy_0[i] = 2.0 * ta2_xz_zz_xy_0[i] * fe_0 - 2.0 * ta2_xz_zz_xy_1[i] * fe_0 + ta2_xz_zz_xyy_0[i] * pa_y[i] - ta2_xz_zz_xyy_1[i] * pc_y[i];

        ta2_xz_yzz_xyz_0[i] = ta2_xz_zz_xz_0[i] * fe_0 - ta2_xz_zz_xz_1[i] * fe_0 + ta2_xz_zz_xyz_0[i] * pa_y[i] - ta2_xz_zz_xyz_1[i] * pc_y[i];

        ta2_xz_yzz_xzz_0[i] = ta2_xz_zz_xzz_0[i] * pa_y[i] - ta2_xz_zz_xzz_1[i] * pc_y[i];

        ta2_xz_yzz_yyy_0[i] = 3.0 * ta2_xz_zz_yy_0[i] * fe_0 - 3.0 * ta2_xz_zz_yy_1[i] * fe_0 + ta2_xz_zz_yyy_0[i] * pa_y[i] - ta2_xz_zz_yyy_1[i] * pc_y[i];

        ta2_xz_yzz_yyz_0[i] = 2.0 * ta2_xz_zz_yz_0[i] * fe_0 - 2.0 * ta2_xz_zz_yz_1[i] * fe_0 + ta2_xz_zz_yyz_0[i] * pa_y[i] - ta2_xz_zz_yyz_1[i] * pc_y[i];

        ta2_xz_yzz_yzz_0[i] = ta2_xz_zz_zz_0[i] * fe_0 - ta2_xz_zz_zz_1[i] * fe_0 + ta2_xz_zz_yzz_0[i] * pa_y[i] - ta2_xz_zz_yzz_1[i] * pc_y[i];

        ta2_xz_yzz_zzz_0[i] = ta2_xz_zz_zzz_0[i] * pa_y[i] - ta2_xz_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 290-300 components of targeted buffer : FF

    auto ta2_xz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 290);

    auto ta2_xz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 291);

    auto ta2_xz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 292);

    auto ta2_xz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 293);

    auto ta2_xz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 294);

    auto ta2_xz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 295);

    auto ta2_xz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 296);

    auto ta2_xz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 297);

    auto ta2_xz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 298);

    auto ta2_xz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 299);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zz_xxx_1, ta1_x_zz_xxy_1, ta1_x_zz_xxz_1, ta1_x_zz_xyy_1, ta1_x_zz_xyz_1, ta1_x_zz_xzz_1, ta1_x_zz_yyy_1, ta1_x_zz_yyz_1, ta1_x_zz_yzz_1, ta1_x_zz_zzz_1, ta2_xz_z_xxx_0, ta2_xz_z_xxx_1, ta2_xz_z_xxy_0, ta2_xz_z_xxy_1, ta2_xz_z_xxz_0, ta2_xz_z_xxz_1, ta2_xz_z_xyy_0, ta2_xz_z_xyy_1, ta2_xz_z_xyz_0, ta2_xz_z_xyz_1, ta2_xz_z_xzz_0, ta2_xz_z_xzz_1, ta2_xz_z_yyy_0, ta2_xz_z_yyy_1, ta2_xz_z_yyz_0, ta2_xz_z_yyz_1, ta2_xz_z_yzz_0, ta2_xz_z_yzz_1, ta2_xz_z_zzz_0, ta2_xz_z_zzz_1, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, ta2_xz_zzz_xxx_0, ta2_xz_zzz_xxy_0, ta2_xz_zzz_xxz_0, ta2_xz_zzz_xyy_0, ta2_xz_zzz_xyz_0, ta2_xz_zzz_xzz_0, ta2_xz_zzz_yyy_0, ta2_xz_zzz_yyz_0, ta2_xz_zzz_yzz_0, ta2_xz_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzz_xxx_0[i] = 2.0 * ta2_xz_z_xxx_0[i] * fe_0 - 2.0 * ta2_xz_z_xxx_1[i] * fe_0 + ta1_x_zz_xxx_1[i] + ta2_xz_zz_xxx_0[i] * pa_z[i] - ta2_xz_zz_xxx_1[i] * pc_z[i];

        ta2_xz_zzz_xxy_0[i] = 2.0 * ta2_xz_z_xxy_0[i] * fe_0 - 2.0 * ta2_xz_z_xxy_1[i] * fe_0 + ta1_x_zz_xxy_1[i] + ta2_xz_zz_xxy_0[i] * pa_z[i] - ta2_xz_zz_xxy_1[i] * pc_z[i];

        ta2_xz_zzz_xxz_0[i] = 2.0 * ta2_xz_z_xxz_0[i] * fe_0 - 2.0 * ta2_xz_z_xxz_1[i] * fe_0 + ta2_xz_zz_xx_0[i] * fe_0 - ta2_xz_zz_xx_1[i] * fe_0 + ta1_x_zz_xxz_1[i] + ta2_xz_zz_xxz_0[i] * pa_z[i] - ta2_xz_zz_xxz_1[i] * pc_z[i];

        ta2_xz_zzz_xyy_0[i] = 2.0 * ta2_xz_z_xyy_0[i] * fe_0 - 2.0 * ta2_xz_z_xyy_1[i] * fe_0 + ta1_x_zz_xyy_1[i] + ta2_xz_zz_xyy_0[i] * pa_z[i] - ta2_xz_zz_xyy_1[i] * pc_z[i];

        ta2_xz_zzz_xyz_0[i] = 2.0 * ta2_xz_z_xyz_0[i] * fe_0 - 2.0 * ta2_xz_z_xyz_1[i] * fe_0 + ta2_xz_zz_xy_0[i] * fe_0 - ta2_xz_zz_xy_1[i] * fe_0 + ta1_x_zz_xyz_1[i] + ta2_xz_zz_xyz_0[i] * pa_z[i] - ta2_xz_zz_xyz_1[i] * pc_z[i];

        ta2_xz_zzz_xzz_0[i] = 2.0 * ta2_xz_z_xzz_0[i] * fe_0 - 2.0 * ta2_xz_z_xzz_1[i] * fe_0 + 2.0 * ta2_xz_zz_xz_0[i] * fe_0 - 2.0 * ta2_xz_zz_xz_1[i] * fe_0 + ta1_x_zz_xzz_1[i] + ta2_xz_zz_xzz_0[i] * pa_z[i] - ta2_xz_zz_xzz_1[i] * pc_z[i];

        ta2_xz_zzz_yyy_0[i] = 2.0 * ta2_xz_z_yyy_0[i] * fe_0 - 2.0 * ta2_xz_z_yyy_1[i] * fe_0 + ta1_x_zz_yyy_1[i] + ta2_xz_zz_yyy_0[i] * pa_z[i] - ta2_xz_zz_yyy_1[i] * pc_z[i];

        ta2_xz_zzz_yyz_0[i] = 2.0 * ta2_xz_z_yyz_0[i] * fe_0 - 2.0 * ta2_xz_z_yyz_1[i] * fe_0 + ta2_xz_zz_yy_0[i] * fe_0 - ta2_xz_zz_yy_1[i] * fe_0 + ta1_x_zz_yyz_1[i] + ta2_xz_zz_yyz_0[i] * pa_z[i] - ta2_xz_zz_yyz_1[i] * pc_z[i];

        ta2_xz_zzz_yzz_0[i] = 2.0 * ta2_xz_z_yzz_0[i] * fe_0 - 2.0 * ta2_xz_z_yzz_1[i] * fe_0 + 2.0 * ta2_xz_zz_yz_0[i] * fe_0 - 2.0 * ta2_xz_zz_yz_1[i] * fe_0 + ta1_x_zz_yzz_1[i] + ta2_xz_zz_yzz_0[i] * pa_z[i] - ta2_xz_zz_yzz_1[i] * pc_z[i];

        ta2_xz_zzz_zzz_0[i] = 2.0 * ta2_xz_z_zzz_0[i] * fe_0 - 2.0 * ta2_xz_z_zzz_1[i] * fe_0 + 3.0 * ta2_xz_zz_zz_0[i] * fe_0 - 3.0 * ta2_xz_zz_zz_1[i] * fe_0 + ta1_x_zz_zzz_1[i] + ta2_xz_zz_zzz_0[i] * pa_z[i] - ta2_xz_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 300-310 components of targeted buffer : FF

    auto ta2_yy_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 300);

    auto ta2_yy_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 301);

    auto ta2_yy_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 302);

    auto ta2_yy_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 303);

    auto ta2_yy_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 304);

    auto ta2_yy_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 305);

    auto ta2_yy_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 306);

    auto ta2_yy_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 307);

    auto ta2_yy_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 308);

    auto ta2_yy_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 309);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_x_xxx_0, ta2_yy_x_xxx_1, ta2_yy_x_xxy_0, ta2_yy_x_xxy_1, ta2_yy_x_xxz_0, ta2_yy_x_xxz_1, ta2_yy_x_xyy_0, ta2_yy_x_xyy_1, ta2_yy_x_xyz_0, ta2_yy_x_xyz_1, ta2_yy_x_xzz_0, ta2_yy_x_xzz_1, ta2_yy_x_yyy_0, ta2_yy_x_yyy_1, ta2_yy_x_yyz_0, ta2_yy_x_yyz_1, ta2_yy_x_yzz_0, ta2_yy_x_yzz_1, ta2_yy_x_zzz_0, ta2_yy_x_zzz_1, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_yy_0, ta2_yy_xx_yy_1, ta2_yy_xx_yyy_0, ta2_yy_xx_yyy_1, ta2_yy_xx_yyz_0, ta2_yy_xx_yyz_1, ta2_yy_xx_yz_0, ta2_yy_xx_yz_1, ta2_yy_xx_yzz_0, ta2_yy_xx_yzz_1, ta2_yy_xx_zz_0, ta2_yy_xx_zz_1, ta2_yy_xx_zzz_0, ta2_yy_xx_zzz_1, ta2_yy_xxx_xxx_0, ta2_yy_xxx_xxy_0, ta2_yy_xxx_xxz_0, ta2_yy_xxx_xyy_0, ta2_yy_xxx_xyz_0, ta2_yy_xxx_xzz_0, ta2_yy_xxx_yyy_0, ta2_yy_xxx_yyz_0, ta2_yy_xxx_yzz_0, ta2_yy_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxx_xxx_0[i] = 2.0 * ta2_yy_x_xxx_0[i] * fe_0 - 2.0 * ta2_yy_x_xxx_1[i] * fe_0 + 3.0 * ta2_yy_xx_xx_0[i] * fe_0 - 3.0 * ta2_yy_xx_xx_1[i] * fe_0 + ta2_yy_xx_xxx_0[i] * pa_x[i] - ta2_yy_xx_xxx_1[i] * pc_x[i];

        ta2_yy_xxx_xxy_0[i] = 2.0 * ta2_yy_x_xxy_0[i] * fe_0 - 2.0 * ta2_yy_x_xxy_1[i] * fe_0 + 2.0 * ta2_yy_xx_xy_0[i] * fe_0 - 2.0 * ta2_yy_xx_xy_1[i] * fe_0 + ta2_yy_xx_xxy_0[i] * pa_x[i] - ta2_yy_xx_xxy_1[i] * pc_x[i];

        ta2_yy_xxx_xxz_0[i] = 2.0 * ta2_yy_x_xxz_0[i] * fe_0 - 2.0 * ta2_yy_x_xxz_1[i] * fe_0 + 2.0 * ta2_yy_xx_xz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xz_1[i] * fe_0 + ta2_yy_xx_xxz_0[i] * pa_x[i] - ta2_yy_xx_xxz_1[i] * pc_x[i];

        ta2_yy_xxx_xyy_0[i] = 2.0 * ta2_yy_x_xyy_0[i] * fe_0 - 2.0 * ta2_yy_x_xyy_1[i] * fe_0 + ta2_yy_xx_yy_0[i] * fe_0 - ta2_yy_xx_yy_1[i] * fe_0 + ta2_yy_xx_xyy_0[i] * pa_x[i] - ta2_yy_xx_xyy_1[i] * pc_x[i];

        ta2_yy_xxx_xyz_0[i] = 2.0 * ta2_yy_x_xyz_0[i] * fe_0 - 2.0 * ta2_yy_x_xyz_1[i] * fe_0 + ta2_yy_xx_yz_0[i] * fe_0 - ta2_yy_xx_yz_1[i] * fe_0 + ta2_yy_xx_xyz_0[i] * pa_x[i] - ta2_yy_xx_xyz_1[i] * pc_x[i];

        ta2_yy_xxx_xzz_0[i] = 2.0 * ta2_yy_x_xzz_0[i] * fe_0 - 2.0 * ta2_yy_x_xzz_1[i] * fe_0 + ta2_yy_xx_zz_0[i] * fe_0 - ta2_yy_xx_zz_1[i] * fe_0 + ta2_yy_xx_xzz_0[i] * pa_x[i] - ta2_yy_xx_xzz_1[i] * pc_x[i];

        ta2_yy_xxx_yyy_0[i] = 2.0 * ta2_yy_x_yyy_0[i] * fe_0 - 2.0 * ta2_yy_x_yyy_1[i] * fe_0 + ta2_yy_xx_yyy_0[i] * pa_x[i] - ta2_yy_xx_yyy_1[i] * pc_x[i];

        ta2_yy_xxx_yyz_0[i] = 2.0 * ta2_yy_x_yyz_0[i] * fe_0 - 2.0 * ta2_yy_x_yyz_1[i] * fe_0 + ta2_yy_xx_yyz_0[i] * pa_x[i] - ta2_yy_xx_yyz_1[i] * pc_x[i];

        ta2_yy_xxx_yzz_0[i] = 2.0 * ta2_yy_x_yzz_0[i] * fe_0 - 2.0 * ta2_yy_x_yzz_1[i] * fe_0 + ta2_yy_xx_yzz_0[i] * pa_x[i] - ta2_yy_xx_yzz_1[i] * pc_x[i];

        ta2_yy_xxx_zzz_0[i] = 2.0 * ta2_yy_x_zzz_0[i] * fe_0 - 2.0 * ta2_yy_x_zzz_1[i] * fe_0 + ta2_yy_xx_zzz_0[i] * pa_x[i] - ta2_yy_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : FF

    auto ta2_yy_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 310);

    auto ta2_yy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 311);

    auto ta2_yy_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 312);

    auto ta2_yy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 313);

    auto ta2_yy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 314);

    auto ta2_yy_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 315);

    auto ta2_yy_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 316);

    auto ta2_yy_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 317);

    auto ta2_yy_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 318);

    auto ta2_yy_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 319);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xx_xxx_1, ta1_y_xx_xxy_1, ta1_y_xx_xxz_1, ta1_y_xx_xyy_1, ta1_y_xx_xyz_1, ta1_y_xx_xzz_1, ta1_y_xx_zzz_1, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_zzz_0, ta2_yy_xx_zzz_1, ta2_yy_xxy_xxx_0, ta2_yy_xxy_xxy_0, ta2_yy_xxy_xxz_0, ta2_yy_xxy_xyy_0, ta2_yy_xxy_xyz_0, ta2_yy_xxy_xzz_0, ta2_yy_xxy_yyy_0, ta2_yy_xxy_yyz_0, ta2_yy_xxy_yzz_0, ta2_yy_xxy_zzz_0, ta2_yy_xy_yyy_0, ta2_yy_xy_yyy_1, ta2_yy_xy_yyz_0, ta2_yy_xy_yyz_1, ta2_yy_xy_yzz_0, ta2_yy_xy_yzz_1, ta2_yy_y_yyy_0, ta2_yy_y_yyy_1, ta2_yy_y_yyz_0, ta2_yy_y_yyz_1, ta2_yy_y_yzz_0, ta2_yy_y_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxy_xxx_0[i] = 2.0 * ta1_y_xx_xxx_1[i] + ta2_yy_xx_xxx_0[i] * pa_y[i] - ta2_yy_xx_xxx_1[i] * pc_y[i];

        ta2_yy_xxy_xxy_0[i] = ta2_yy_xx_xx_0[i] * fe_0 - ta2_yy_xx_xx_1[i] * fe_0 + 2.0 * ta1_y_xx_xxy_1[i] + ta2_yy_xx_xxy_0[i] * pa_y[i] - ta2_yy_xx_xxy_1[i] * pc_y[i];

        ta2_yy_xxy_xxz_0[i] = 2.0 * ta1_y_xx_xxz_1[i] + ta2_yy_xx_xxz_0[i] * pa_y[i] - ta2_yy_xx_xxz_1[i] * pc_y[i];

        ta2_yy_xxy_xyy_0[i] = 2.0 * ta2_yy_xx_xy_0[i] * fe_0 - 2.0 * ta2_yy_xx_xy_1[i] * fe_0 + 2.0 * ta1_y_xx_xyy_1[i] + ta2_yy_xx_xyy_0[i] * pa_y[i] - ta2_yy_xx_xyy_1[i] * pc_y[i];

        ta2_yy_xxy_xyz_0[i] = ta2_yy_xx_xz_0[i] * fe_0 - ta2_yy_xx_xz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyz_1[i] + ta2_yy_xx_xyz_0[i] * pa_y[i] - ta2_yy_xx_xyz_1[i] * pc_y[i];

        ta2_yy_xxy_xzz_0[i] = 2.0 * ta1_y_xx_xzz_1[i] + ta2_yy_xx_xzz_0[i] * pa_y[i] - ta2_yy_xx_xzz_1[i] * pc_y[i];

        ta2_yy_xxy_yyy_0[i] = ta2_yy_y_yyy_0[i] * fe_0 - ta2_yy_y_yyy_1[i] * fe_0 + ta2_yy_xy_yyy_0[i] * pa_x[i] - ta2_yy_xy_yyy_1[i] * pc_x[i];

        ta2_yy_xxy_yyz_0[i] = ta2_yy_y_yyz_0[i] * fe_0 - ta2_yy_y_yyz_1[i] * fe_0 + ta2_yy_xy_yyz_0[i] * pa_x[i] - ta2_yy_xy_yyz_1[i] * pc_x[i];

        ta2_yy_xxy_yzz_0[i] = ta2_yy_y_yzz_0[i] * fe_0 - ta2_yy_y_yzz_1[i] * fe_0 + ta2_yy_xy_yzz_0[i] * pa_x[i] - ta2_yy_xy_yzz_1[i] * pc_x[i];

        ta2_yy_xxy_zzz_0[i] = 2.0 * ta1_y_xx_zzz_1[i] + ta2_yy_xx_zzz_0[i] * pa_y[i] - ta2_yy_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 320-330 components of targeted buffer : FF

    auto ta2_yy_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 320);

    auto ta2_yy_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 321);

    auto ta2_yy_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 322);

    auto ta2_yy_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 323);

    auto ta2_yy_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 324);

    auto ta2_yy_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 325);

    auto ta2_yy_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 326);

    auto ta2_yy_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 327);

    auto ta2_yy_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 328);

    auto ta2_yy_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 329);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_yyy_0, ta2_yy_xx_yyy_1, ta2_yy_xxz_xxx_0, ta2_yy_xxz_xxy_0, ta2_yy_xxz_xxz_0, ta2_yy_xxz_xyy_0, ta2_yy_xxz_xyz_0, ta2_yy_xxz_xzz_0, ta2_yy_xxz_yyy_0, ta2_yy_xxz_yyz_0, ta2_yy_xxz_yzz_0, ta2_yy_xxz_zzz_0, ta2_yy_xz_yyz_0, ta2_yy_xz_yyz_1, ta2_yy_xz_yzz_0, ta2_yy_xz_yzz_1, ta2_yy_xz_zzz_0, ta2_yy_xz_zzz_1, ta2_yy_z_yyz_0, ta2_yy_z_yyz_1, ta2_yy_z_yzz_0, ta2_yy_z_yzz_1, ta2_yy_z_zzz_0, ta2_yy_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxz_xxx_0[i] = ta2_yy_xx_xxx_0[i] * pa_z[i] - ta2_yy_xx_xxx_1[i] * pc_z[i];

        ta2_yy_xxz_xxy_0[i] = ta2_yy_xx_xxy_0[i] * pa_z[i] - ta2_yy_xx_xxy_1[i] * pc_z[i];

        ta2_yy_xxz_xxz_0[i] = ta2_yy_xx_xx_0[i] * fe_0 - ta2_yy_xx_xx_1[i] * fe_0 + ta2_yy_xx_xxz_0[i] * pa_z[i] - ta2_yy_xx_xxz_1[i] * pc_z[i];

        ta2_yy_xxz_xyy_0[i] = ta2_yy_xx_xyy_0[i] * pa_z[i] - ta2_yy_xx_xyy_1[i] * pc_z[i];

        ta2_yy_xxz_xyz_0[i] = ta2_yy_xx_xy_0[i] * fe_0 - ta2_yy_xx_xy_1[i] * fe_0 + ta2_yy_xx_xyz_0[i] * pa_z[i] - ta2_yy_xx_xyz_1[i] * pc_z[i];

        ta2_yy_xxz_xzz_0[i] = 2.0 * ta2_yy_xx_xz_0[i] * fe_0 - 2.0 * ta2_yy_xx_xz_1[i] * fe_0 + ta2_yy_xx_xzz_0[i] * pa_z[i] - ta2_yy_xx_xzz_1[i] * pc_z[i];

        ta2_yy_xxz_yyy_0[i] = ta2_yy_xx_yyy_0[i] * pa_z[i] - ta2_yy_xx_yyy_1[i] * pc_z[i];

        ta2_yy_xxz_yyz_0[i] = ta2_yy_z_yyz_0[i] * fe_0 - ta2_yy_z_yyz_1[i] * fe_0 + ta2_yy_xz_yyz_0[i] * pa_x[i] - ta2_yy_xz_yyz_1[i] * pc_x[i];

        ta2_yy_xxz_yzz_0[i] = ta2_yy_z_yzz_0[i] * fe_0 - ta2_yy_z_yzz_1[i] * fe_0 + ta2_yy_xz_yzz_0[i] * pa_x[i] - ta2_yy_xz_yzz_1[i] * pc_x[i];

        ta2_yy_xxz_zzz_0[i] = ta2_yy_z_zzz_0[i] * fe_0 - ta2_yy_z_zzz_1[i] * fe_0 + ta2_yy_xz_zzz_0[i] * pa_x[i] - ta2_yy_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 330-340 components of targeted buffer : FF

    auto ta2_yy_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 330);

    auto ta2_yy_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 331);

    auto ta2_yy_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 332);

    auto ta2_yy_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 333);

    auto ta2_yy_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 334);

    auto ta2_yy_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 335);

    auto ta2_yy_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 336);

    auto ta2_yy_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 337);

    auto ta2_yy_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 338);

    auto ta2_yy_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 339);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xyy_xxx_0, ta2_yy_xyy_xxy_0, ta2_yy_xyy_xxz_0, ta2_yy_xyy_xyy_0, ta2_yy_xyy_xyz_0, ta2_yy_xyy_xzz_0, ta2_yy_xyy_yyy_0, ta2_yy_xyy_yyz_0, ta2_yy_xyy_yzz_0, ta2_yy_xyy_zzz_0, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyy_xxx_0[i] = 3.0 * ta2_yy_yy_xx_0[i] * fe_0 - 3.0 * ta2_yy_yy_xx_1[i] * fe_0 + ta2_yy_yy_xxx_0[i] * pa_x[i] - ta2_yy_yy_xxx_1[i] * pc_x[i];

        ta2_yy_xyy_xxy_0[i] = 2.0 * ta2_yy_yy_xy_0[i] * fe_0 - 2.0 * ta2_yy_yy_xy_1[i] * fe_0 + ta2_yy_yy_xxy_0[i] * pa_x[i] - ta2_yy_yy_xxy_1[i] * pc_x[i];

        ta2_yy_xyy_xxz_0[i] = 2.0 * ta2_yy_yy_xz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xz_1[i] * fe_0 + ta2_yy_yy_xxz_0[i] * pa_x[i] - ta2_yy_yy_xxz_1[i] * pc_x[i];

        ta2_yy_xyy_xyy_0[i] = ta2_yy_yy_yy_0[i] * fe_0 - ta2_yy_yy_yy_1[i] * fe_0 + ta2_yy_yy_xyy_0[i] * pa_x[i] - ta2_yy_yy_xyy_1[i] * pc_x[i];

        ta2_yy_xyy_xyz_0[i] = ta2_yy_yy_yz_0[i] * fe_0 - ta2_yy_yy_yz_1[i] * fe_0 + ta2_yy_yy_xyz_0[i] * pa_x[i] - ta2_yy_yy_xyz_1[i] * pc_x[i];

        ta2_yy_xyy_xzz_0[i] = ta2_yy_yy_zz_0[i] * fe_0 - ta2_yy_yy_zz_1[i] * fe_0 + ta2_yy_yy_xzz_0[i] * pa_x[i] - ta2_yy_yy_xzz_1[i] * pc_x[i];

        ta2_yy_xyy_yyy_0[i] = ta2_yy_yy_yyy_0[i] * pa_x[i] - ta2_yy_yy_yyy_1[i] * pc_x[i];

        ta2_yy_xyy_yyz_0[i] = ta2_yy_yy_yyz_0[i] * pa_x[i] - ta2_yy_yy_yyz_1[i] * pc_x[i];

        ta2_yy_xyy_yzz_0[i] = ta2_yy_yy_yzz_0[i] * pa_x[i] - ta2_yy_yy_yzz_1[i] * pc_x[i];

        ta2_yy_xyy_zzz_0[i] = ta2_yy_yy_zzz_0[i] * pa_x[i] - ta2_yy_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 340-350 components of targeted buffer : FF

    auto ta2_yy_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 340);

    auto ta2_yy_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 341);

    auto ta2_yy_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 342);

    auto ta2_yy_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 343);

    auto ta2_yy_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 344);

    auto ta2_yy_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 345);

    auto ta2_yy_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 346);

    auto ta2_yy_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 347);

    auto ta2_yy_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 348);

    auto ta2_yy_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 349);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xz_xxz_1, ta1_y_xz_xzz_1, ta2_yy_xy_xxx_0, ta2_yy_xy_xxx_1, ta2_yy_xy_xxy_0, ta2_yy_xy_xxy_1, ta2_yy_xy_xyy_0, ta2_yy_xy_xyy_1, ta2_yy_xyz_xxx_0, ta2_yy_xyz_xxy_0, ta2_yy_xyz_xxz_0, ta2_yy_xyz_xyy_0, ta2_yy_xyz_xyz_0, ta2_yy_xyz_xzz_0, ta2_yy_xyz_yyy_0, ta2_yy_xyz_yyz_0, ta2_yy_xyz_yzz_0, ta2_yy_xyz_zzz_0, ta2_yy_xz_xxz_0, ta2_yy_xz_xxz_1, ta2_yy_xz_xzz_0, ta2_yy_xz_xzz_1, ta2_yy_yz_xyz_0, ta2_yy_yz_xyz_1, ta2_yy_yz_yyy_0, ta2_yy_yz_yyy_1, ta2_yy_yz_yyz_0, ta2_yy_yz_yyz_1, ta2_yy_yz_yz_0, ta2_yy_yz_yz_1, ta2_yy_yz_yzz_0, ta2_yy_yz_yzz_1, ta2_yy_yz_zzz_0, ta2_yy_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyz_xxx_0[i] = ta2_yy_xy_xxx_0[i] * pa_z[i] - ta2_yy_xy_xxx_1[i] * pc_z[i];

        ta2_yy_xyz_xxy_0[i] = ta2_yy_xy_xxy_0[i] * pa_z[i] - ta2_yy_xy_xxy_1[i] * pc_z[i];

        ta2_yy_xyz_xxz_0[i] = 2.0 * ta1_y_xz_xxz_1[i] + ta2_yy_xz_xxz_0[i] * pa_y[i] - ta2_yy_xz_xxz_1[i] * pc_y[i];

        ta2_yy_xyz_xyy_0[i] = ta2_yy_xy_xyy_0[i] * pa_z[i] - ta2_yy_xy_xyy_1[i] * pc_z[i];

        ta2_yy_xyz_xyz_0[i] = ta2_yy_yz_yz_0[i] * fe_0 - ta2_yy_yz_yz_1[i] * fe_0 + ta2_yy_yz_xyz_0[i] * pa_x[i] - ta2_yy_yz_xyz_1[i] * pc_x[i];

        ta2_yy_xyz_xzz_0[i] = 2.0 * ta1_y_xz_xzz_1[i] + ta2_yy_xz_xzz_0[i] * pa_y[i] - ta2_yy_xz_xzz_1[i] * pc_y[i];

        ta2_yy_xyz_yyy_0[i] = ta2_yy_yz_yyy_0[i] * pa_x[i] - ta2_yy_yz_yyy_1[i] * pc_x[i];

        ta2_yy_xyz_yyz_0[i] = ta2_yy_yz_yyz_0[i] * pa_x[i] - ta2_yy_yz_yyz_1[i] * pc_x[i];

        ta2_yy_xyz_yzz_0[i] = ta2_yy_yz_yzz_0[i] * pa_x[i] - ta2_yy_yz_yzz_1[i] * pc_x[i];

        ta2_yy_xyz_zzz_0[i] = ta2_yy_yz_zzz_0[i] * pa_x[i] - ta2_yy_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 350-360 components of targeted buffer : FF

    auto ta2_yy_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 350);

    auto ta2_yy_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 351);

    auto ta2_yy_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 352);

    auto ta2_yy_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 353);

    auto ta2_yy_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 354);

    auto ta2_yy_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 355);

    auto ta2_yy_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 356);

    auto ta2_yy_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 357);

    auto ta2_yy_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 358);

    auto ta2_yy_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 359);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xzz_xxx_0, ta2_yy_xzz_xxy_0, ta2_yy_xzz_xxz_0, ta2_yy_xzz_xyy_0, ta2_yy_xzz_xyz_0, ta2_yy_xzz_xzz_0, ta2_yy_xzz_yyy_0, ta2_yy_xzz_yyz_0, ta2_yy_xzz_yzz_0, ta2_yy_xzz_zzz_0, ta2_yy_zz_xx_0, ta2_yy_zz_xx_1, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxy_0, ta2_yy_zz_xxy_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xy_0, ta2_yy_zz_xy_1, ta2_yy_zz_xyy_0, ta2_yy_zz_xyy_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_yy_0, ta2_yy_zz_yy_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzz_xxx_0[i] = 3.0 * ta2_yy_zz_xx_0[i] * fe_0 - 3.0 * ta2_yy_zz_xx_1[i] * fe_0 + ta2_yy_zz_xxx_0[i] * pa_x[i] - ta2_yy_zz_xxx_1[i] * pc_x[i];

        ta2_yy_xzz_xxy_0[i] = 2.0 * ta2_yy_zz_xy_0[i] * fe_0 - 2.0 * ta2_yy_zz_xy_1[i] * fe_0 + ta2_yy_zz_xxy_0[i] * pa_x[i] - ta2_yy_zz_xxy_1[i] * pc_x[i];

        ta2_yy_xzz_xxz_0[i] = 2.0 * ta2_yy_zz_xz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xz_1[i] * fe_0 + ta2_yy_zz_xxz_0[i] * pa_x[i] - ta2_yy_zz_xxz_1[i] * pc_x[i];

        ta2_yy_xzz_xyy_0[i] = ta2_yy_zz_yy_0[i] * fe_0 - ta2_yy_zz_yy_1[i] * fe_0 + ta2_yy_zz_xyy_0[i] * pa_x[i] - ta2_yy_zz_xyy_1[i] * pc_x[i];

        ta2_yy_xzz_xyz_0[i] = ta2_yy_zz_yz_0[i] * fe_0 - ta2_yy_zz_yz_1[i] * fe_0 + ta2_yy_zz_xyz_0[i] * pa_x[i] - ta2_yy_zz_xyz_1[i] * pc_x[i];

        ta2_yy_xzz_xzz_0[i] = ta2_yy_zz_zz_0[i] * fe_0 - ta2_yy_zz_zz_1[i] * fe_0 + ta2_yy_zz_xzz_0[i] * pa_x[i] - ta2_yy_zz_xzz_1[i] * pc_x[i];

        ta2_yy_xzz_yyy_0[i] = ta2_yy_zz_yyy_0[i] * pa_x[i] - ta2_yy_zz_yyy_1[i] * pc_x[i];

        ta2_yy_xzz_yyz_0[i] = ta2_yy_zz_yyz_0[i] * pa_x[i] - ta2_yy_zz_yyz_1[i] * pc_x[i];

        ta2_yy_xzz_yzz_0[i] = ta2_yy_zz_yzz_0[i] * pa_x[i] - ta2_yy_zz_yzz_1[i] * pc_x[i];

        ta2_yy_xzz_zzz_0[i] = ta2_yy_zz_zzz_0[i] * pa_x[i] - ta2_yy_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 360-370 components of targeted buffer : FF

    auto ta2_yy_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 360);

    auto ta2_yy_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 361);

    auto ta2_yy_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 362);

    auto ta2_yy_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 363);

    auto ta2_yy_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 364);

    auto ta2_yy_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 365);

    auto ta2_yy_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 366);

    auto ta2_yy_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 367);

    auto ta2_yy_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 368);

    auto ta2_yy_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 369);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yy_xxx_1, ta1_y_yy_xxy_1, ta1_y_yy_xxz_1, ta1_y_yy_xyy_1, ta1_y_yy_xyz_1, ta1_y_yy_xzz_1, ta1_y_yy_yyy_1, ta1_y_yy_yyz_1, ta1_y_yy_yzz_1, ta1_y_yy_zzz_1, ta2_yy_y_xxx_0, ta2_yy_y_xxx_1, ta2_yy_y_xxy_0, ta2_yy_y_xxy_1, ta2_yy_y_xxz_0, ta2_yy_y_xxz_1, ta2_yy_y_xyy_0, ta2_yy_y_xyy_1, ta2_yy_y_xyz_0, ta2_yy_y_xyz_1, ta2_yy_y_xzz_0, ta2_yy_y_xzz_1, ta2_yy_y_yyy_0, ta2_yy_y_yyy_1, ta2_yy_y_yyz_0, ta2_yy_y_yyz_1, ta2_yy_y_yzz_0, ta2_yy_y_yzz_1, ta2_yy_y_zzz_0, ta2_yy_y_zzz_1, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yyy_xxx_0, ta2_yy_yyy_xxy_0, ta2_yy_yyy_xxz_0, ta2_yy_yyy_xyy_0, ta2_yy_yyy_xyz_0, ta2_yy_yyy_xzz_0, ta2_yy_yyy_yyy_0, ta2_yy_yyy_yyz_0, ta2_yy_yyy_yzz_0, ta2_yy_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyy_xxx_0[i] = 2.0 * ta2_yy_y_xxx_0[i] * fe_0 - 2.0 * ta2_yy_y_xxx_1[i] * fe_0 + 2.0 * ta1_y_yy_xxx_1[i] + ta2_yy_yy_xxx_0[i] * pa_y[i] - ta2_yy_yy_xxx_1[i] * pc_y[i];

        ta2_yy_yyy_xxy_0[i] = 2.0 * ta2_yy_y_xxy_0[i] * fe_0 - 2.0 * ta2_yy_y_xxy_1[i] * fe_0 + ta2_yy_yy_xx_0[i] * fe_0 - ta2_yy_yy_xx_1[i] * fe_0 + 2.0 * ta1_y_yy_xxy_1[i] + ta2_yy_yy_xxy_0[i] * pa_y[i] - ta2_yy_yy_xxy_1[i] * pc_y[i];

        ta2_yy_yyy_xxz_0[i] = 2.0 * ta2_yy_y_xxz_0[i] * fe_0 - 2.0 * ta2_yy_y_xxz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxz_1[i] + ta2_yy_yy_xxz_0[i] * pa_y[i] - ta2_yy_yy_xxz_1[i] * pc_y[i];

        ta2_yy_yyy_xyy_0[i] = 2.0 * ta2_yy_y_xyy_0[i] * fe_0 - 2.0 * ta2_yy_y_xyy_1[i] * fe_0 + 2.0 * ta2_yy_yy_xy_0[i] * fe_0 - 2.0 * ta2_yy_yy_xy_1[i] * fe_0 + 2.0 * ta1_y_yy_xyy_1[i] + ta2_yy_yy_xyy_0[i] * pa_y[i] - ta2_yy_yy_xyy_1[i] * pc_y[i];

        ta2_yy_yyy_xyz_0[i] = 2.0 * ta2_yy_y_xyz_0[i] * fe_0 - 2.0 * ta2_yy_y_xyz_1[i] * fe_0 + ta2_yy_yy_xz_0[i] * fe_0 - ta2_yy_yy_xz_1[i] * fe_0 + 2.0 * ta1_y_yy_xyz_1[i] + ta2_yy_yy_xyz_0[i] * pa_y[i] - ta2_yy_yy_xyz_1[i] * pc_y[i];

        ta2_yy_yyy_xzz_0[i] = 2.0 * ta2_yy_y_xzz_0[i] * fe_0 - 2.0 * ta2_yy_y_xzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xzz_1[i] + ta2_yy_yy_xzz_0[i] * pa_y[i] - ta2_yy_yy_xzz_1[i] * pc_y[i];

        ta2_yy_yyy_yyy_0[i] = 2.0 * ta2_yy_y_yyy_0[i] * fe_0 - 2.0 * ta2_yy_y_yyy_1[i] * fe_0 + 3.0 * ta2_yy_yy_yy_0[i] * fe_0 - 3.0 * ta2_yy_yy_yy_1[i] * fe_0 + 2.0 * ta1_y_yy_yyy_1[i] + ta2_yy_yy_yyy_0[i] * pa_y[i] - ta2_yy_yy_yyy_1[i] * pc_y[i];

        ta2_yy_yyy_yyz_0[i] = 2.0 * ta2_yy_y_yyz_0[i] * fe_0 - 2.0 * ta2_yy_y_yyz_1[i] * fe_0 + 2.0 * ta2_yy_yy_yz_0[i] * fe_0 - 2.0 * ta2_yy_yy_yz_1[i] * fe_0 + 2.0 * ta1_y_yy_yyz_1[i] + ta2_yy_yy_yyz_0[i] * pa_y[i] - ta2_yy_yy_yyz_1[i] * pc_y[i];

        ta2_yy_yyy_yzz_0[i] = 2.0 * ta2_yy_y_yzz_0[i] * fe_0 - 2.0 * ta2_yy_y_yzz_1[i] * fe_0 + ta2_yy_yy_zz_0[i] * fe_0 - ta2_yy_yy_zz_1[i] * fe_0 + 2.0 * ta1_y_yy_yzz_1[i] + ta2_yy_yy_yzz_0[i] * pa_y[i] - ta2_yy_yy_yzz_1[i] * pc_y[i];

        ta2_yy_yyy_zzz_0[i] = 2.0 * ta2_yy_y_zzz_0[i] * fe_0 - 2.0 * ta2_yy_y_zzz_1[i] * fe_0 + 2.0 * ta1_y_yy_zzz_1[i] + ta2_yy_yy_zzz_0[i] * pa_y[i] - ta2_yy_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 370-380 components of targeted buffer : FF

    auto ta2_yy_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 370);

    auto ta2_yy_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 371);

    auto ta2_yy_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 372);

    auto ta2_yy_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 373);

    auto ta2_yy_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 374);

    auto ta2_yy_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 375);

    auto ta2_yy_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 376);

    auto ta2_yy_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 377);

    auto ta2_yy_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 378);

    auto ta2_yy_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 379);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yyz_xxx_0, ta2_yy_yyz_xxy_0, ta2_yy_yyz_xxz_0, ta2_yy_yyz_xyy_0, ta2_yy_yyz_xyz_0, ta2_yy_yyz_xzz_0, ta2_yy_yyz_yyy_0, ta2_yy_yyz_yyz_0, ta2_yy_yyz_yzz_0, ta2_yy_yyz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyz_xxx_0[i] = ta2_yy_yy_xxx_0[i] * pa_z[i] - ta2_yy_yy_xxx_1[i] * pc_z[i];

        ta2_yy_yyz_xxy_0[i] = ta2_yy_yy_xxy_0[i] * pa_z[i] - ta2_yy_yy_xxy_1[i] * pc_z[i];

        ta2_yy_yyz_xxz_0[i] = ta2_yy_yy_xx_0[i] * fe_0 - ta2_yy_yy_xx_1[i] * fe_0 + ta2_yy_yy_xxz_0[i] * pa_z[i] - ta2_yy_yy_xxz_1[i] * pc_z[i];

        ta2_yy_yyz_xyy_0[i] = ta2_yy_yy_xyy_0[i] * pa_z[i] - ta2_yy_yy_xyy_1[i] * pc_z[i];

        ta2_yy_yyz_xyz_0[i] = ta2_yy_yy_xy_0[i] * fe_0 - ta2_yy_yy_xy_1[i] * fe_0 + ta2_yy_yy_xyz_0[i] * pa_z[i] - ta2_yy_yy_xyz_1[i] * pc_z[i];

        ta2_yy_yyz_xzz_0[i] = 2.0 * ta2_yy_yy_xz_0[i] * fe_0 - 2.0 * ta2_yy_yy_xz_1[i] * fe_0 + ta2_yy_yy_xzz_0[i] * pa_z[i] - ta2_yy_yy_xzz_1[i] * pc_z[i];

        ta2_yy_yyz_yyy_0[i] = ta2_yy_yy_yyy_0[i] * pa_z[i] - ta2_yy_yy_yyy_1[i] * pc_z[i];

        ta2_yy_yyz_yyz_0[i] = ta2_yy_yy_yy_0[i] * fe_0 - ta2_yy_yy_yy_1[i] * fe_0 + ta2_yy_yy_yyz_0[i] * pa_z[i] - ta2_yy_yy_yyz_1[i] * pc_z[i];

        ta2_yy_yyz_yzz_0[i] = 2.0 * ta2_yy_yy_yz_0[i] * fe_0 - 2.0 * ta2_yy_yy_yz_1[i] * fe_0 + ta2_yy_yy_yzz_0[i] * pa_z[i] - ta2_yy_yy_yzz_1[i] * pc_z[i];

        ta2_yy_yyz_zzz_0[i] = 3.0 * ta2_yy_yy_zz_0[i] * fe_0 - 3.0 * ta2_yy_yy_zz_1[i] * fe_0 + ta2_yy_yy_zzz_0[i] * pa_z[i] - ta2_yy_yy_zzz_1[i] * pc_z[i];
    }

    // Set up 380-390 components of targeted buffer : FF

    auto ta2_yy_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 380);

    auto ta2_yy_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 381);

    auto ta2_yy_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 382);

    auto ta2_yy_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 383);

    auto ta2_yy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 384);

    auto ta2_yy_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 385);

    auto ta2_yy_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 386);

    auto ta2_yy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 387);

    auto ta2_yy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 388);

    auto ta2_yy_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 389);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_zz_xxx_1, ta1_y_zz_xxz_1, ta1_y_zz_xyz_1, ta1_y_zz_xzz_1, ta1_y_zz_yyz_1, ta1_y_zz_yzz_1, ta1_y_zz_zzz_1, ta2_yy_y_xxy_0, ta2_yy_y_xxy_1, ta2_yy_y_xyy_0, ta2_yy_y_xyy_1, ta2_yy_y_yyy_0, ta2_yy_y_yyy_1, ta2_yy_yz_xxy_0, ta2_yy_yz_xxy_1, ta2_yy_yz_xyy_0, ta2_yy_yz_xyy_1, ta2_yy_yz_yyy_0, ta2_yy_yz_yyy_1, ta2_yy_yzz_xxx_0, ta2_yy_yzz_xxy_0, ta2_yy_yzz_xxz_0, ta2_yy_yzz_xyy_0, ta2_yy_yzz_xyz_0, ta2_yy_yzz_xzz_0, ta2_yy_yzz_yyy_0, ta2_yy_yzz_yyz_0, ta2_yy_yzz_yzz_0, ta2_yy_yzz_zzz_0, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzz_xxx_0[i] = 2.0 * ta1_y_zz_xxx_1[i] + ta2_yy_zz_xxx_0[i] * pa_y[i] - ta2_yy_zz_xxx_1[i] * pc_y[i];

        ta2_yy_yzz_xxy_0[i] = ta2_yy_y_xxy_0[i] * fe_0 - ta2_yy_y_xxy_1[i] * fe_0 + ta2_yy_yz_xxy_0[i] * pa_z[i] - ta2_yy_yz_xxy_1[i] * pc_z[i];

        ta2_yy_yzz_xxz_0[i] = 2.0 * ta1_y_zz_xxz_1[i] + ta2_yy_zz_xxz_0[i] * pa_y[i] - ta2_yy_zz_xxz_1[i] * pc_y[i];

        ta2_yy_yzz_xyy_0[i] = ta2_yy_y_xyy_0[i] * fe_0 - ta2_yy_y_xyy_1[i] * fe_0 + ta2_yy_yz_xyy_0[i] * pa_z[i] - ta2_yy_yz_xyy_1[i] * pc_z[i];

        ta2_yy_yzz_xyz_0[i] = ta2_yy_zz_xz_0[i] * fe_0 - ta2_yy_zz_xz_1[i] * fe_0 + 2.0 * ta1_y_zz_xyz_1[i] + ta2_yy_zz_xyz_0[i] * pa_y[i] - ta2_yy_zz_xyz_1[i] * pc_y[i];

        ta2_yy_yzz_xzz_0[i] = 2.0 * ta1_y_zz_xzz_1[i] + ta2_yy_zz_xzz_0[i] * pa_y[i] - ta2_yy_zz_xzz_1[i] * pc_y[i];

        ta2_yy_yzz_yyy_0[i] = ta2_yy_y_yyy_0[i] * fe_0 - ta2_yy_y_yyy_1[i] * fe_0 + ta2_yy_yz_yyy_0[i] * pa_z[i] - ta2_yy_yz_yyy_1[i] * pc_z[i];

        ta2_yy_yzz_yyz_0[i] = 2.0 * ta2_yy_zz_yz_0[i] * fe_0 - 2.0 * ta2_yy_zz_yz_1[i] * fe_0 + 2.0 * ta1_y_zz_yyz_1[i] + ta2_yy_zz_yyz_0[i] * pa_y[i] - ta2_yy_zz_yyz_1[i] * pc_y[i];

        ta2_yy_yzz_yzz_0[i] = ta2_yy_zz_zz_0[i] * fe_0 - ta2_yy_zz_zz_1[i] * fe_0 + 2.0 * ta1_y_zz_yzz_1[i] + ta2_yy_zz_yzz_0[i] * pa_y[i] - ta2_yy_zz_yzz_1[i] * pc_y[i];

        ta2_yy_yzz_zzz_0[i] = 2.0 * ta1_y_zz_zzz_1[i] + ta2_yy_zz_zzz_0[i] * pa_y[i] - ta2_yy_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 390-400 components of targeted buffer : FF

    auto ta2_yy_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 390);

    auto ta2_yy_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 391);

    auto ta2_yy_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 392);

    auto ta2_yy_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 393);

    auto ta2_yy_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 394);

    auto ta2_yy_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 395);

    auto ta2_yy_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 396);

    auto ta2_yy_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 397);

    auto ta2_yy_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 398);

    auto ta2_yy_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 399);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_z_xxx_0, ta2_yy_z_xxx_1, ta2_yy_z_xxy_0, ta2_yy_z_xxy_1, ta2_yy_z_xxz_0, ta2_yy_z_xxz_1, ta2_yy_z_xyy_0, ta2_yy_z_xyy_1, ta2_yy_z_xyz_0, ta2_yy_z_xyz_1, ta2_yy_z_xzz_0, ta2_yy_z_xzz_1, ta2_yy_z_yyy_0, ta2_yy_z_yyy_1, ta2_yy_z_yyz_0, ta2_yy_z_yyz_1, ta2_yy_z_yzz_0, ta2_yy_z_yzz_1, ta2_yy_z_zzz_0, ta2_yy_z_zzz_1, ta2_yy_zz_xx_0, ta2_yy_zz_xx_1, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxy_0, ta2_yy_zz_xxy_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xy_0, ta2_yy_zz_xy_1, ta2_yy_zz_xyy_0, ta2_yy_zz_xyy_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_yy_0, ta2_yy_zz_yy_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, ta2_yy_zzz_xxx_0, ta2_yy_zzz_xxy_0, ta2_yy_zzz_xxz_0, ta2_yy_zzz_xyy_0, ta2_yy_zzz_xyz_0, ta2_yy_zzz_xzz_0, ta2_yy_zzz_yyy_0, ta2_yy_zzz_yyz_0, ta2_yy_zzz_yzz_0, ta2_yy_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzz_xxx_0[i] = 2.0 * ta2_yy_z_xxx_0[i] * fe_0 - 2.0 * ta2_yy_z_xxx_1[i] * fe_0 + ta2_yy_zz_xxx_0[i] * pa_z[i] - ta2_yy_zz_xxx_1[i] * pc_z[i];

        ta2_yy_zzz_xxy_0[i] = 2.0 * ta2_yy_z_xxy_0[i] * fe_0 - 2.0 * ta2_yy_z_xxy_1[i] * fe_0 + ta2_yy_zz_xxy_0[i] * pa_z[i] - ta2_yy_zz_xxy_1[i] * pc_z[i];

        ta2_yy_zzz_xxz_0[i] = 2.0 * ta2_yy_z_xxz_0[i] * fe_0 - 2.0 * ta2_yy_z_xxz_1[i] * fe_0 + ta2_yy_zz_xx_0[i] * fe_0 - ta2_yy_zz_xx_1[i] * fe_0 + ta2_yy_zz_xxz_0[i] * pa_z[i] - ta2_yy_zz_xxz_1[i] * pc_z[i];

        ta2_yy_zzz_xyy_0[i] = 2.0 * ta2_yy_z_xyy_0[i] * fe_0 - 2.0 * ta2_yy_z_xyy_1[i] * fe_0 + ta2_yy_zz_xyy_0[i] * pa_z[i] - ta2_yy_zz_xyy_1[i] * pc_z[i];

        ta2_yy_zzz_xyz_0[i] = 2.0 * ta2_yy_z_xyz_0[i] * fe_0 - 2.0 * ta2_yy_z_xyz_1[i] * fe_0 + ta2_yy_zz_xy_0[i] * fe_0 - ta2_yy_zz_xy_1[i] * fe_0 + ta2_yy_zz_xyz_0[i] * pa_z[i] - ta2_yy_zz_xyz_1[i] * pc_z[i];

        ta2_yy_zzz_xzz_0[i] = 2.0 * ta2_yy_z_xzz_0[i] * fe_0 - 2.0 * ta2_yy_z_xzz_1[i] * fe_0 + 2.0 * ta2_yy_zz_xz_0[i] * fe_0 - 2.0 * ta2_yy_zz_xz_1[i] * fe_0 + ta2_yy_zz_xzz_0[i] * pa_z[i] - ta2_yy_zz_xzz_1[i] * pc_z[i];

        ta2_yy_zzz_yyy_0[i] = 2.0 * ta2_yy_z_yyy_0[i] * fe_0 - 2.0 * ta2_yy_z_yyy_1[i] * fe_0 + ta2_yy_zz_yyy_0[i] * pa_z[i] - ta2_yy_zz_yyy_1[i] * pc_z[i];

        ta2_yy_zzz_yyz_0[i] = 2.0 * ta2_yy_z_yyz_0[i] * fe_0 - 2.0 * ta2_yy_z_yyz_1[i] * fe_0 + ta2_yy_zz_yy_0[i] * fe_0 - ta2_yy_zz_yy_1[i] * fe_0 + ta2_yy_zz_yyz_0[i] * pa_z[i] - ta2_yy_zz_yyz_1[i] * pc_z[i];

        ta2_yy_zzz_yzz_0[i] = 2.0 * ta2_yy_z_yzz_0[i] * fe_0 - 2.0 * ta2_yy_z_yzz_1[i] * fe_0 + 2.0 * ta2_yy_zz_yz_0[i] * fe_0 - 2.0 * ta2_yy_zz_yz_1[i] * fe_0 + ta2_yy_zz_yzz_0[i] * pa_z[i] - ta2_yy_zz_yzz_1[i] * pc_z[i];

        ta2_yy_zzz_zzz_0[i] = 2.0 * ta2_yy_z_zzz_0[i] * fe_0 - 2.0 * ta2_yy_z_zzz_1[i] * fe_0 + 3.0 * ta2_yy_zz_zz_0[i] * fe_0 - 3.0 * ta2_yy_zz_zz_1[i] * fe_0 + ta2_yy_zz_zzz_0[i] * pa_z[i] - ta2_yy_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 400-410 components of targeted buffer : FF

    auto ta2_yz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 400);

    auto ta2_yz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 401);

    auto ta2_yz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 402);

    auto ta2_yz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 403);

    auto ta2_yz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 404);

    auto ta2_yz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 405);

    auto ta2_yz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 406);

    auto ta2_yz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 407);

    auto ta2_yz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 408);

    auto ta2_yz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 409);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_x_xxx_0, ta2_yz_x_xxx_1, ta2_yz_x_xxy_0, ta2_yz_x_xxy_1, ta2_yz_x_xxz_0, ta2_yz_x_xxz_1, ta2_yz_x_xyy_0, ta2_yz_x_xyy_1, ta2_yz_x_xyz_0, ta2_yz_x_xyz_1, ta2_yz_x_xzz_0, ta2_yz_x_xzz_1, ta2_yz_x_yyy_0, ta2_yz_x_yyy_1, ta2_yz_x_yyz_0, ta2_yz_x_yyz_1, ta2_yz_x_yzz_0, ta2_yz_x_yzz_1, ta2_yz_x_zzz_0, ta2_yz_x_zzz_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_yy_0, ta2_yz_xx_yy_1, ta2_yz_xx_yyy_0, ta2_yz_xx_yyy_1, ta2_yz_xx_yyz_0, ta2_yz_xx_yyz_1, ta2_yz_xx_yz_0, ta2_yz_xx_yz_1, ta2_yz_xx_yzz_0, ta2_yz_xx_yzz_1, ta2_yz_xx_zz_0, ta2_yz_xx_zz_1, ta2_yz_xx_zzz_0, ta2_yz_xx_zzz_1, ta2_yz_xxx_xxx_0, ta2_yz_xxx_xxy_0, ta2_yz_xxx_xxz_0, ta2_yz_xxx_xyy_0, ta2_yz_xxx_xyz_0, ta2_yz_xxx_xzz_0, ta2_yz_xxx_yyy_0, ta2_yz_xxx_yyz_0, ta2_yz_xxx_yzz_0, ta2_yz_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxx_xxx_0[i] = 2.0 * ta2_yz_x_xxx_0[i] * fe_0 - 2.0 * ta2_yz_x_xxx_1[i] * fe_0 + 3.0 * ta2_yz_xx_xx_0[i] * fe_0 - 3.0 * ta2_yz_xx_xx_1[i] * fe_0 + ta2_yz_xx_xxx_0[i] * pa_x[i] - ta2_yz_xx_xxx_1[i] * pc_x[i];

        ta2_yz_xxx_xxy_0[i] = 2.0 * ta2_yz_x_xxy_0[i] * fe_0 - 2.0 * ta2_yz_x_xxy_1[i] * fe_0 + 2.0 * ta2_yz_xx_xy_0[i] * fe_0 - 2.0 * ta2_yz_xx_xy_1[i] * fe_0 + ta2_yz_xx_xxy_0[i] * pa_x[i] - ta2_yz_xx_xxy_1[i] * pc_x[i];

        ta2_yz_xxx_xxz_0[i] = 2.0 * ta2_yz_x_xxz_0[i] * fe_0 - 2.0 * ta2_yz_x_xxz_1[i] * fe_0 + 2.0 * ta2_yz_xx_xz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xz_1[i] * fe_0 + ta2_yz_xx_xxz_0[i] * pa_x[i] - ta2_yz_xx_xxz_1[i] * pc_x[i];

        ta2_yz_xxx_xyy_0[i] = 2.0 * ta2_yz_x_xyy_0[i] * fe_0 - 2.0 * ta2_yz_x_xyy_1[i] * fe_0 + ta2_yz_xx_yy_0[i] * fe_0 - ta2_yz_xx_yy_1[i] * fe_0 + ta2_yz_xx_xyy_0[i] * pa_x[i] - ta2_yz_xx_xyy_1[i] * pc_x[i];

        ta2_yz_xxx_xyz_0[i] = 2.0 * ta2_yz_x_xyz_0[i] * fe_0 - 2.0 * ta2_yz_x_xyz_1[i] * fe_0 + ta2_yz_xx_yz_0[i] * fe_0 - ta2_yz_xx_yz_1[i] * fe_0 + ta2_yz_xx_xyz_0[i] * pa_x[i] - ta2_yz_xx_xyz_1[i] * pc_x[i];

        ta2_yz_xxx_xzz_0[i] = 2.0 * ta2_yz_x_xzz_0[i] * fe_0 - 2.0 * ta2_yz_x_xzz_1[i] * fe_0 + ta2_yz_xx_zz_0[i] * fe_0 - ta2_yz_xx_zz_1[i] * fe_0 + ta2_yz_xx_xzz_0[i] * pa_x[i] - ta2_yz_xx_xzz_1[i] * pc_x[i];

        ta2_yz_xxx_yyy_0[i] = 2.0 * ta2_yz_x_yyy_0[i] * fe_0 - 2.0 * ta2_yz_x_yyy_1[i] * fe_0 + ta2_yz_xx_yyy_0[i] * pa_x[i] - ta2_yz_xx_yyy_1[i] * pc_x[i];

        ta2_yz_xxx_yyz_0[i] = 2.0 * ta2_yz_x_yyz_0[i] * fe_0 - 2.0 * ta2_yz_x_yyz_1[i] * fe_0 + ta2_yz_xx_yyz_0[i] * pa_x[i] - ta2_yz_xx_yyz_1[i] * pc_x[i];

        ta2_yz_xxx_yzz_0[i] = 2.0 * ta2_yz_x_yzz_0[i] * fe_0 - 2.0 * ta2_yz_x_yzz_1[i] * fe_0 + ta2_yz_xx_yzz_0[i] * pa_x[i] - ta2_yz_xx_yzz_1[i] * pc_x[i];

        ta2_yz_xxx_zzz_0[i] = 2.0 * ta2_yz_x_zzz_0[i] * fe_0 - 2.0 * ta2_yz_x_zzz_1[i] * fe_0 + ta2_yz_xx_zzz_0[i] * pa_x[i] - ta2_yz_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 410-420 components of targeted buffer : FF

    auto ta2_yz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 410);

    auto ta2_yz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 411);

    auto ta2_yz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 412);

    auto ta2_yz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 413);

    auto ta2_yz_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 414);

    auto ta2_yz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 415);

    auto ta2_yz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 416);

    auto ta2_yz_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 417);

    auto ta2_yz_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 418);

    auto ta2_yz_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 419);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xx_xxx_1, ta1_z_xx_xxy_1, ta1_z_xx_xxz_1, ta1_z_xx_xyy_1, ta1_z_xx_xyz_1, ta1_z_xx_xzz_1, ta1_z_xx_zzz_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_zzz_0, ta2_yz_xx_zzz_1, ta2_yz_xxy_xxx_0, ta2_yz_xxy_xxy_0, ta2_yz_xxy_xxz_0, ta2_yz_xxy_xyy_0, ta2_yz_xxy_xyz_0, ta2_yz_xxy_xzz_0, ta2_yz_xxy_yyy_0, ta2_yz_xxy_yyz_0, ta2_yz_xxy_yzz_0, ta2_yz_xxy_zzz_0, ta2_yz_xy_yyy_0, ta2_yz_xy_yyy_1, ta2_yz_xy_yyz_0, ta2_yz_xy_yyz_1, ta2_yz_xy_yzz_0, ta2_yz_xy_yzz_1, ta2_yz_y_yyy_0, ta2_yz_y_yyy_1, ta2_yz_y_yyz_0, ta2_yz_y_yyz_1, ta2_yz_y_yzz_0, ta2_yz_y_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxy_xxx_0[i] = ta1_z_xx_xxx_1[i] + ta2_yz_xx_xxx_0[i] * pa_y[i] - ta2_yz_xx_xxx_1[i] * pc_y[i];

        ta2_yz_xxy_xxy_0[i] = ta2_yz_xx_xx_0[i] * fe_0 - ta2_yz_xx_xx_1[i] * fe_0 + ta1_z_xx_xxy_1[i] + ta2_yz_xx_xxy_0[i] * pa_y[i] - ta2_yz_xx_xxy_1[i] * pc_y[i];

        ta2_yz_xxy_xxz_0[i] = ta1_z_xx_xxz_1[i] + ta2_yz_xx_xxz_0[i] * pa_y[i] - ta2_yz_xx_xxz_1[i] * pc_y[i];

        ta2_yz_xxy_xyy_0[i] = 2.0 * ta2_yz_xx_xy_0[i] * fe_0 - 2.0 * ta2_yz_xx_xy_1[i] * fe_0 + ta1_z_xx_xyy_1[i] + ta2_yz_xx_xyy_0[i] * pa_y[i] - ta2_yz_xx_xyy_1[i] * pc_y[i];

        ta2_yz_xxy_xyz_0[i] = ta2_yz_xx_xz_0[i] * fe_0 - ta2_yz_xx_xz_1[i] * fe_0 + ta1_z_xx_xyz_1[i] + ta2_yz_xx_xyz_0[i] * pa_y[i] - ta2_yz_xx_xyz_1[i] * pc_y[i];

        ta2_yz_xxy_xzz_0[i] = ta1_z_xx_xzz_1[i] + ta2_yz_xx_xzz_0[i] * pa_y[i] - ta2_yz_xx_xzz_1[i] * pc_y[i];

        ta2_yz_xxy_yyy_0[i] = ta2_yz_y_yyy_0[i] * fe_0 - ta2_yz_y_yyy_1[i] * fe_0 + ta2_yz_xy_yyy_0[i] * pa_x[i] - ta2_yz_xy_yyy_1[i] * pc_x[i];

        ta2_yz_xxy_yyz_0[i] = ta2_yz_y_yyz_0[i] * fe_0 - ta2_yz_y_yyz_1[i] * fe_0 + ta2_yz_xy_yyz_0[i] * pa_x[i] - ta2_yz_xy_yyz_1[i] * pc_x[i];

        ta2_yz_xxy_yzz_0[i] = ta2_yz_y_yzz_0[i] * fe_0 - ta2_yz_y_yzz_1[i] * fe_0 + ta2_yz_xy_yzz_0[i] * pa_x[i] - ta2_yz_xy_yzz_1[i] * pc_x[i];

        ta2_yz_xxy_zzz_0[i] = ta1_z_xx_zzz_1[i] + ta2_yz_xx_zzz_0[i] * pa_y[i] - ta2_yz_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 420-430 components of targeted buffer : FF

    auto ta2_yz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 420);

    auto ta2_yz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 421);

    auto ta2_yz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 422);

    auto ta2_yz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 423);

    auto ta2_yz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 424);

    auto ta2_yz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 425);

    auto ta2_yz_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 426);

    auto ta2_yz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 427);

    auto ta2_yz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 428);

    auto ta2_yz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 429);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xx_xxx_1, ta1_y_xx_xxy_1, ta1_y_xx_xxz_1, ta1_y_xx_xyy_1, ta1_y_xx_xyz_1, ta1_y_xx_xzz_1, ta1_y_xx_yyy_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_yyy_0, ta2_yz_xx_yyy_1, ta2_yz_xxz_xxx_0, ta2_yz_xxz_xxy_0, ta2_yz_xxz_xxz_0, ta2_yz_xxz_xyy_0, ta2_yz_xxz_xyz_0, ta2_yz_xxz_xzz_0, ta2_yz_xxz_yyy_0, ta2_yz_xxz_yyz_0, ta2_yz_xxz_yzz_0, ta2_yz_xxz_zzz_0, ta2_yz_xz_yyz_0, ta2_yz_xz_yyz_1, ta2_yz_xz_yzz_0, ta2_yz_xz_yzz_1, ta2_yz_xz_zzz_0, ta2_yz_xz_zzz_1, ta2_yz_z_yyz_0, ta2_yz_z_yyz_1, ta2_yz_z_yzz_0, ta2_yz_z_yzz_1, ta2_yz_z_zzz_0, ta2_yz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxz_xxx_0[i] = ta1_y_xx_xxx_1[i] + ta2_yz_xx_xxx_0[i] * pa_z[i] - ta2_yz_xx_xxx_1[i] * pc_z[i];

        ta2_yz_xxz_xxy_0[i] = ta1_y_xx_xxy_1[i] + ta2_yz_xx_xxy_0[i] * pa_z[i] - ta2_yz_xx_xxy_1[i] * pc_z[i];

        ta2_yz_xxz_xxz_0[i] = ta2_yz_xx_xx_0[i] * fe_0 - ta2_yz_xx_xx_1[i] * fe_0 + ta1_y_xx_xxz_1[i] + ta2_yz_xx_xxz_0[i] * pa_z[i] - ta2_yz_xx_xxz_1[i] * pc_z[i];

        ta2_yz_xxz_xyy_0[i] = ta1_y_xx_xyy_1[i] + ta2_yz_xx_xyy_0[i] * pa_z[i] - ta2_yz_xx_xyy_1[i] * pc_z[i];

        ta2_yz_xxz_xyz_0[i] = ta2_yz_xx_xy_0[i] * fe_0 - ta2_yz_xx_xy_1[i] * fe_0 + ta1_y_xx_xyz_1[i] + ta2_yz_xx_xyz_0[i] * pa_z[i] - ta2_yz_xx_xyz_1[i] * pc_z[i];

        ta2_yz_xxz_xzz_0[i] = 2.0 * ta2_yz_xx_xz_0[i] * fe_0 - 2.0 * ta2_yz_xx_xz_1[i] * fe_0 + ta1_y_xx_xzz_1[i] + ta2_yz_xx_xzz_0[i] * pa_z[i] - ta2_yz_xx_xzz_1[i] * pc_z[i];

        ta2_yz_xxz_yyy_0[i] = ta1_y_xx_yyy_1[i] + ta2_yz_xx_yyy_0[i] * pa_z[i] - ta2_yz_xx_yyy_1[i] * pc_z[i];

        ta2_yz_xxz_yyz_0[i] = ta2_yz_z_yyz_0[i] * fe_0 - ta2_yz_z_yyz_1[i] * fe_0 + ta2_yz_xz_yyz_0[i] * pa_x[i] - ta2_yz_xz_yyz_1[i] * pc_x[i];

        ta2_yz_xxz_yzz_0[i] = ta2_yz_z_yzz_0[i] * fe_0 - ta2_yz_z_yzz_1[i] * fe_0 + ta2_yz_xz_yzz_0[i] * pa_x[i] - ta2_yz_xz_yzz_1[i] * pc_x[i];

        ta2_yz_xxz_zzz_0[i] = ta2_yz_z_zzz_0[i] * fe_0 - ta2_yz_z_zzz_1[i] * fe_0 + ta2_yz_xz_zzz_0[i] * pa_x[i] - ta2_yz_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 430-440 components of targeted buffer : FF

    auto ta2_yz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 430);

    auto ta2_yz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 431);

    auto ta2_yz_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 432);

    auto ta2_yz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 433);

    auto ta2_yz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 434);

    auto ta2_yz_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 435);

    auto ta2_yz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 436);

    auto ta2_yz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 437);

    auto ta2_yz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 438);

    auto ta2_yz_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 439);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xyy_xxx_0, ta2_yz_xyy_xxy_0, ta2_yz_xyy_xxz_0, ta2_yz_xyy_xyy_0, ta2_yz_xyy_xyz_0, ta2_yz_xyy_xzz_0, ta2_yz_xyy_yyy_0, ta2_yz_xyy_yyz_0, ta2_yz_xyy_yzz_0, ta2_yz_xyy_zzz_0, ta2_yz_yy_xx_0, ta2_yz_yy_xx_1, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxz_0, ta2_yz_yy_xxz_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xz_0, ta2_yz_yy_xz_1, ta2_yz_yy_xzz_0, ta2_yz_yy_xzz_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_zz_0, ta2_yz_yy_zz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyy_xxx_0[i] = 3.0 * ta2_yz_yy_xx_0[i] * fe_0 - 3.0 * ta2_yz_yy_xx_1[i] * fe_0 + ta2_yz_yy_xxx_0[i] * pa_x[i] - ta2_yz_yy_xxx_1[i] * pc_x[i];

        ta2_yz_xyy_xxy_0[i] = 2.0 * ta2_yz_yy_xy_0[i] * fe_0 - 2.0 * ta2_yz_yy_xy_1[i] * fe_0 + ta2_yz_yy_xxy_0[i] * pa_x[i] - ta2_yz_yy_xxy_1[i] * pc_x[i];

        ta2_yz_xyy_xxz_0[i] = 2.0 * ta2_yz_yy_xz_0[i] * fe_0 - 2.0 * ta2_yz_yy_xz_1[i] * fe_0 + ta2_yz_yy_xxz_0[i] * pa_x[i] - ta2_yz_yy_xxz_1[i] * pc_x[i];

        ta2_yz_xyy_xyy_0[i] = ta2_yz_yy_yy_0[i] * fe_0 - ta2_yz_yy_yy_1[i] * fe_0 + ta2_yz_yy_xyy_0[i] * pa_x[i] - ta2_yz_yy_xyy_1[i] * pc_x[i];

        ta2_yz_xyy_xyz_0[i] = ta2_yz_yy_yz_0[i] * fe_0 - ta2_yz_yy_yz_1[i] * fe_0 + ta2_yz_yy_xyz_0[i] * pa_x[i] - ta2_yz_yy_xyz_1[i] * pc_x[i];

        ta2_yz_xyy_xzz_0[i] = ta2_yz_yy_zz_0[i] * fe_0 - ta2_yz_yy_zz_1[i] * fe_0 + ta2_yz_yy_xzz_0[i] * pa_x[i] - ta2_yz_yy_xzz_1[i] * pc_x[i];

        ta2_yz_xyy_yyy_0[i] = ta2_yz_yy_yyy_0[i] * pa_x[i] - ta2_yz_yy_yyy_1[i] * pc_x[i];

        ta2_yz_xyy_yyz_0[i] = ta2_yz_yy_yyz_0[i] * pa_x[i] - ta2_yz_yy_yyz_1[i] * pc_x[i];

        ta2_yz_xyy_yzz_0[i] = ta2_yz_yy_yzz_0[i] * pa_x[i] - ta2_yz_yy_yzz_1[i] * pc_x[i];

        ta2_yz_xyy_zzz_0[i] = ta2_yz_yy_zzz_0[i] * pa_x[i] - ta2_yz_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 440-450 components of targeted buffer : FF

    auto ta2_yz_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 440);

    auto ta2_yz_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 441);

    auto ta2_yz_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 442);

    auto ta2_yz_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 443);

    auto ta2_yz_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 444);

    auto ta2_yz_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 445);

    auto ta2_yz_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 446);

    auto ta2_yz_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 447);

    auto ta2_yz_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 448);

    auto ta2_yz_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 449);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xy_xxy_1, ta1_y_xy_xyy_1, ta1_z_xz_xxx_1, ta1_z_xz_xxz_1, ta1_z_xz_xzz_1, ta2_yz_xy_xxy_0, ta2_yz_xy_xxy_1, ta2_yz_xy_xyy_0, ta2_yz_xy_xyy_1, ta2_yz_xyz_xxx_0, ta2_yz_xyz_xxy_0, ta2_yz_xyz_xxz_0, ta2_yz_xyz_xyy_0, ta2_yz_xyz_xyz_0, ta2_yz_xyz_xzz_0, ta2_yz_xyz_yyy_0, ta2_yz_xyz_yyz_0, ta2_yz_xyz_yzz_0, ta2_yz_xyz_zzz_0, ta2_yz_xz_xxx_0, ta2_yz_xz_xxx_1, ta2_yz_xz_xxz_0, ta2_yz_xz_xxz_1, ta2_yz_xz_xzz_0, ta2_yz_xz_xzz_1, ta2_yz_yz_xyz_0, ta2_yz_yz_xyz_1, ta2_yz_yz_yyy_0, ta2_yz_yz_yyy_1, ta2_yz_yz_yyz_0, ta2_yz_yz_yyz_1, ta2_yz_yz_yz_0, ta2_yz_yz_yz_1, ta2_yz_yz_yzz_0, ta2_yz_yz_yzz_1, ta2_yz_yz_zzz_0, ta2_yz_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyz_xxx_0[i] = ta1_z_xz_xxx_1[i] + ta2_yz_xz_xxx_0[i] * pa_y[i] - ta2_yz_xz_xxx_1[i] * pc_y[i];

        ta2_yz_xyz_xxy_0[i] = ta1_y_xy_xxy_1[i] + ta2_yz_xy_xxy_0[i] * pa_z[i] - ta2_yz_xy_xxy_1[i] * pc_z[i];

        ta2_yz_xyz_xxz_0[i] = ta1_z_xz_xxz_1[i] + ta2_yz_xz_xxz_0[i] * pa_y[i] - ta2_yz_xz_xxz_1[i] * pc_y[i];

        ta2_yz_xyz_xyy_0[i] = ta1_y_xy_xyy_1[i] + ta2_yz_xy_xyy_0[i] * pa_z[i] - ta2_yz_xy_xyy_1[i] * pc_z[i];

        ta2_yz_xyz_xyz_0[i] = ta2_yz_yz_yz_0[i] * fe_0 - ta2_yz_yz_yz_1[i] * fe_0 + ta2_yz_yz_xyz_0[i] * pa_x[i] - ta2_yz_yz_xyz_1[i] * pc_x[i];

        ta2_yz_xyz_xzz_0[i] = ta1_z_xz_xzz_1[i] + ta2_yz_xz_xzz_0[i] * pa_y[i] - ta2_yz_xz_xzz_1[i] * pc_y[i];

        ta2_yz_xyz_yyy_0[i] = ta2_yz_yz_yyy_0[i] * pa_x[i] - ta2_yz_yz_yyy_1[i] * pc_x[i];

        ta2_yz_xyz_yyz_0[i] = ta2_yz_yz_yyz_0[i] * pa_x[i] - ta2_yz_yz_yyz_1[i] * pc_x[i];

        ta2_yz_xyz_yzz_0[i] = ta2_yz_yz_yzz_0[i] * pa_x[i] - ta2_yz_yz_yzz_1[i] * pc_x[i];

        ta2_yz_xyz_zzz_0[i] = ta2_yz_yz_zzz_0[i] * pa_x[i] - ta2_yz_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 450-460 components of targeted buffer : FF

    auto ta2_yz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 450);

    auto ta2_yz_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 451);

    auto ta2_yz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 452);

    auto ta2_yz_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 453);

    auto ta2_yz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 454);

    auto ta2_yz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 455);

    auto ta2_yz_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 456);

    auto ta2_yz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 457);

    auto ta2_yz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 458);

    auto ta2_yz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 459);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xzz_xxx_0, ta2_yz_xzz_xxy_0, ta2_yz_xzz_xxz_0, ta2_yz_xzz_xyy_0, ta2_yz_xzz_xyz_0, ta2_yz_xzz_xzz_0, ta2_yz_xzz_yyy_0, ta2_yz_xzz_yyz_0, ta2_yz_xzz_yzz_0, ta2_yz_xzz_zzz_0, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzz_xxx_0[i] = 3.0 * ta2_yz_zz_xx_0[i] * fe_0 - 3.0 * ta2_yz_zz_xx_1[i] * fe_0 + ta2_yz_zz_xxx_0[i] * pa_x[i] - ta2_yz_zz_xxx_1[i] * pc_x[i];

        ta2_yz_xzz_xxy_0[i] = 2.0 * ta2_yz_zz_xy_0[i] * fe_0 - 2.0 * ta2_yz_zz_xy_1[i] * fe_0 + ta2_yz_zz_xxy_0[i] * pa_x[i] - ta2_yz_zz_xxy_1[i] * pc_x[i];

        ta2_yz_xzz_xxz_0[i] = 2.0 * ta2_yz_zz_xz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xz_1[i] * fe_0 + ta2_yz_zz_xxz_0[i] * pa_x[i] - ta2_yz_zz_xxz_1[i] * pc_x[i];

        ta2_yz_xzz_xyy_0[i] = ta2_yz_zz_yy_0[i] * fe_0 - ta2_yz_zz_yy_1[i] * fe_0 + ta2_yz_zz_xyy_0[i] * pa_x[i] - ta2_yz_zz_xyy_1[i] * pc_x[i];

        ta2_yz_xzz_xyz_0[i] = ta2_yz_zz_yz_0[i] * fe_0 - ta2_yz_zz_yz_1[i] * fe_0 + ta2_yz_zz_xyz_0[i] * pa_x[i] - ta2_yz_zz_xyz_1[i] * pc_x[i];

        ta2_yz_xzz_xzz_0[i] = ta2_yz_zz_zz_0[i] * fe_0 - ta2_yz_zz_zz_1[i] * fe_0 + ta2_yz_zz_xzz_0[i] * pa_x[i] - ta2_yz_zz_xzz_1[i] * pc_x[i];

        ta2_yz_xzz_yyy_0[i] = ta2_yz_zz_yyy_0[i] * pa_x[i] - ta2_yz_zz_yyy_1[i] * pc_x[i];

        ta2_yz_xzz_yyz_0[i] = ta2_yz_zz_yyz_0[i] * pa_x[i] - ta2_yz_zz_yyz_1[i] * pc_x[i];

        ta2_yz_xzz_yzz_0[i] = ta2_yz_zz_yzz_0[i] * pa_x[i] - ta2_yz_zz_yzz_1[i] * pc_x[i];

        ta2_yz_xzz_zzz_0[i] = ta2_yz_zz_zzz_0[i] * pa_x[i] - ta2_yz_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 460-470 components of targeted buffer : FF

    auto ta2_yz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 460);

    auto ta2_yz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 461);

    auto ta2_yz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 462);

    auto ta2_yz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 463);

    auto ta2_yz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 464);

    auto ta2_yz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 465);

    auto ta2_yz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 466);

    auto ta2_yz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 467);

    auto ta2_yz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 468);

    auto ta2_yz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 469);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yy_xxx_1, ta1_z_yy_xxy_1, ta1_z_yy_xxz_1, ta1_z_yy_xyy_1, ta1_z_yy_xyz_1, ta1_z_yy_xzz_1, ta1_z_yy_yyy_1, ta1_z_yy_yyz_1, ta1_z_yy_yzz_1, ta1_z_yy_zzz_1, ta2_yz_y_xxx_0, ta2_yz_y_xxx_1, ta2_yz_y_xxy_0, ta2_yz_y_xxy_1, ta2_yz_y_xxz_0, ta2_yz_y_xxz_1, ta2_yz_y_xyy_0, ta2_yz_y_xyy_1, ta2_yz_y_xyz_0, ta2_yz_y_xyz_1, ta2_yz_y_xzz_0, ta2_yz_y_xzz_1, ta2_yz_y_yyy_0, ta2_yz_y_yyy_1, ta2_yz_y_yyz_0, ta2_yz_y_yyz_1, ta2_yz_y_yzz_0, ta2_yz_y_yzz_1, ta2_yz_y_zzz_0, ta2_yz_y_zzz_1, ta2_yz_yy_xx_0, ta2_yz_yy_xx_1, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxz_0, ta2_yz_yy_xxz_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xz_0, ta2_yz_yy_xz_1, ta2_yz_yy_xzz_0, ta2_yz_yy_xzz_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_zz_0, ta2_yz_yy_zz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, ta2_yz_yyy_xxx_0, ta2_yz_yyy_xxy_0, ta2_yz_yyy_xxz_0, ta2_yz_yyy_xyy_0, ta2_yz_yyy_xyz_0, ta2_yz_yyy_xzz_0, ta2_yz_yyy_yyy_0, ta2_yz_yyy_yyz_0, ta2_yz_yyy_yzz_0, ta2_yz_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyy_xxx_0[i] = 2.0 * ta2_yz_y_xxx_0[i] * fe_0 - 2.0 * ta2_yz_y_xxx_1[i] * fe_0 + ta1_z_yy_xxx_1[i] + ta2_yz_yy_xxx_0[i] * pa_y[i] - ta2_yz_yy_xxx_1[i] * pc_y[i];

        ta2_yz_yyy_xxy_0[i] = 2.0 * ta2_yz_y_xxy_0[i] * fe_0 - 2.0 * ta2_yz_y_xxy_1[i] * fe_0 + ta2_yz_yy_xx_0[i] * fe_0 - ta2_yz_yy_xx_1[i] * fe_0 + ta1_z_yy_xxy_1[i] + ta2_yz_yy_xxy_0[i] * pa_y[i] - ta2_yz_yy_xxy_1[i] * pc_y[i];

        ta2_yz_yyy_xxz_0[i] = 2.0 * ta2_yz_y_xxz_0[i] * fe_0 - 2.0 * ta2_yz_y_xxz_1[i] * fe_0 + ta1_z_yy_xxz_1[i] + ta2_yz_yy_xxz_0[i] * pa_y[i] - ta2_yz_yy_xxz_1[i] * pc_y[i];

        ta2_yz_yyy_xyy_0[i] = 2.0 * ta2_yz_y_xyy_0[i] * fe_0 - 2.0 * ta2_yz_y_xyy_1[i] * fe_0 + 2.0 * ta2_yz_yy_xy_0[i] * fe_0 - 2.0 * ta2_yz_yy_xy_1[i] * fe_0 + ta1_z_yy_xyy_1[i] + ta2_yz_yy_xyy_0[i] * pa_y[i] - ta2_yz_yy_xyy_1[i] * pc_y[i];

        ta2_yz_yyy_xyz_0[i] = 2.0 * ta2_yz_y_xyz_0[i] * fe_0 - 2.0 * ta2_yz_y_xyz_1[i] * fe_0 + ta2_yz_yy_xz_0[i] * fe_0 - ta2_yz_yy_xz_1[i] * fe_0 + ta1_z_yy_xyz_1[i] + ta2_yz_yy_xyz_0[i] * pa_y[i] - ta2_yz_yy_xyz_1[i] * pc_y[i];

        ta2_yz_yyy_xzz_0[i] = 2.0 * ta2_yz_y_xzz_0[i] * fe_0 - 2.0 * ta2_yz_y_xzz_1[i] * fe_0 + ta1_z_yy_xzz_1[i] + ta2_yz_yy_xzz_0[i] * pa_y[i] - ta2_yz_yy_xzz_1[i] * pc_y[i];

        ta2_yz_yyy_yyy_0[i] = 2.0 * ta2_yz_y_yyy_0[i] * fe_0 - 2.0 * ta2_yz_y_yyy_1[i] * fe_0 + 3.0 * ta2_yz_yy_yy_0[i] * fe_0 - 3.0 * ta2_yz_yy_yy_1[i] * fe_0 + ta1_z_yy_yyy_1[i] + ta2_yz_yy_yyy_0[i] * pa_y[i] - ta2_yz_yy_yyy_1[i] * pc_y[i];

        ta2_yz_yyy_yyz_0[i] = 2.0 * ta2_yz_y_yyz_0[i] * fe_0 - 2.0 * ta2_yz_y_yyz_1[i] * fe_0 + 2.0 * ta2_yz_yy_yz_0[i] * fe_0 - 2.0 * ta2_yz_yy_yz_1[i] * fe_0 + ta1_z_yy_yyz_1[i] + ta2_yz_yy_yyz_0[i] * pa_y[i] - ta2_yz_yy_yyz_1[i] * pc_y[i];

        ta2_yz_yyy_yzz_0[i] = 2.0 * ta2_yz_y_yzz_0[i] * fe_0 - 2.0 * ta2_yz_y_yzz_1[i] * fe_0 + ta2_yz_yy_zz_0[i] * fe_0 - ta2_yz_yy_zz_1[i] * fe_0 + ta1_z_yy_yzz_1[i] + ta2_yz_yy_yzz_0[i] * pa_y[i] - ta2_yz_yy_yzz_1[i] * pc_y[i];

        ta2_yz_yyy_zzz_0[i] = 2.0 * ta2_yz_y_zzz_0[i] * fe_0 - 2.0 * ta2_yz_y_zzz_1[i] * fe_0 + ta1_z_yy_zzz_1[i] + ta2_yz_yy_zzz_0[i] * pa_y[i] - ta2_yz_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 470-480 components of targeted buffer : FF

    auto ta2_yz_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 470);

    auto ta2_yz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 471);

    auto ta2_yz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 472);

    auto ta2_yz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 473);

    auto ta2_yz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 474);

    auto ta2_yz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 475);

    auto ta2_yz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 476);

    auto ta2_yz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 477);

    auto ta2_yz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 478);

    auto ta2_yz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 479);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yy_xxx_1, ta1_y_yy_xxy_1, ta1_y_yy_xyy_1, ta1_y_yy_xyz_1, ta1_y_yy_yyy_1, ta1_y_yy_yyz_1, ta1_y_yy_yzz_1, ta1_z_yz_xxz_1, ta1_z_yz_xzz_1, ta1_z_yz_zzz_1, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yyz_xxx_0, ta2_yz_yyz_xxy_0, ta2_yz_yyz_xxz_0, ta2_yz_yyz_xyy_0, ta2_yz_yyz_xyz_0, ta2_yz_yyz_xzz_0, ta2_yz_yyz_yyy_0, ta2_yz_yyz_yyz_0, ta2_yz_yyz_yzz_0, ta2_yz_yyz_zzz_0, ta2_yz_yz_xxz_0, ta2_yz_yz_xxz_1, ta2_yz_yz_xzz_0, ta2_yz_yz_xzz_1, ta2_yz_yz_zzz_0, ta2_yz_yz_zzz_1, ta2_yz_z_xxz_0, ta2_yz_z_xxz_1, ta2_yz_z_xzz_0, ta2_yz_z_xzz_1, ta2_yz_z_zzz_0, ta2_yz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyz_xxx_0[i] = ta1_y_yy_xxx_1[i] + ta2_yz_yy_xxx_0[i] * pa_z[i] - ta2_yz_yy_xxx_1[i] * pc_z[i];

        ta2_yz_yyz_xxy_0[i] = ta1_y_yy_xxy_1[i] + ta2_yz_yy_xxy_0[i] * pa_z[i] - ta2_yz_yy_xxy_1[i] * pc_z[i];

        ta2_yz_yyz_xxz_0[i] = ta2_yz_z_xxz_0[i] * fe_0 - ta2_yz_z_xxz_1[i] * fe_0 + ta1_z_yz_xxz_1[i] + ta2_yz_yz_xxz_0[i] * pa_y[i] - ta2_yz_yz_xxz_1[i] * pc_y[i];

        ta2_yz_yyz_xyy_0[i] = ta1_y_yy_xyy_1[i] + ta2_yz_yy_xyy_0[i] * pa_z[i] - ta2_yz_yy_xyy_1[i] * pc_z[i];

        ta2_yz_yyz_xyz_0[i] = ta2_yz_yy_xy_0[i] * fe_0 - ta2_yz_yy_xy_1[i] * fe_0 + ta1_y_yy_xyz_1[i] + ta2_yz_yy_xyz_0[i] * pa_z[i] - ta2_yz_yy_xyz_1[i] * pc_z[i];

        ta2_yz_yyz_xzz_0[i] = ta2_yz_z_xzz_0[i] * fe_0 - ta2_yz_z_xzz_1[i] * fe_0 + ta1_z_yz_xzz_1[i] + ta2_yz_yz_xzz_0[i] * pa_y[i] - ta2_yz_yz_xzz_1[i] * pc_y[i];

        ta2_yz_yyz_yyy_0[i] = ta1_y_yy_yyy_1[i] + ta2_yz_yy_yyy_0[i] * pa_z[i] - ta2_yz_yy_yyy_1[i] * pc_z[i];

        ta2_yz_yyz_yyz_0[i] = ta2_yz_yy_yy_0[i] * fe_0 - ta2_yz_yy_yy_1[i] * fe_0 + ta1_y_yy_yyz_1[i] + ta2_yz_yy_yyz_0[i] * pa_z[i] - ta2_yz_yy_yyz_1[i] * pc_z[i];

        ta2_yz_yyz_yzz_0[i] = 2.0 * ta2_yz_yy_yz_0[i] * fe_0 - 2.0 * ta2_yz_yy_yz_1[i] * fe_0 + ta1_y_yy_yzz_1[i] + ta2_yz_yy_yzz_0[i] * pa_z[i] - ta2_yz_yy_yzz_1[i] * pc_z[i];

        ta2_yz_yyz_zzz_0[i] = ta2_yz_z_zzz_0[i] * fe_0 - ta2_yz_z_zzz_1[i] * fe_0 + ta1_z_yz_zzz_1[i] + ta2_yz_yz_zzz_0[i] * pa_y[i] - ta2_yz_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 480-490 components of targeted buffer : FF

    auto ta2_yz_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 480);

    auto ta2_yz_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 481);

    auto ta2_yz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 482);

    auto ta2_yz_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 483);

    auto ta2_yz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 484);

    auto ta2_yz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 485);

    auto ta2_yz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 486);

    auto ta2_yz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 487);

    auto ta2_yz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 488);

    auto ta2_yz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 489);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_zz_xxx_1, ta1_z_zz_xxy_1, ta1_z_zz_xxz_1, ta1_z_zz_xyy_1, ta1_z_zz_xyz_1, ta1_z_zz_xzz_1, ta1_z_zz_yyy_1, ta1_z_zz_yyz_1, ta1_z_zz_yzz_1, ta1_z_zz_zzz_1, ta2_yz_yzz_xxx_0, ta2_yz_yzz_xxy_0, ta2_yz_yzz_xxz_0, ta2_yz_yzz_xyy_0, ta2_yz_yzz_xyz_0, ta2_yz_yzz_xzz_0, ta2_yz_yzz_yyy_0, ta2_yz_yzz_yyz_0, ta2_yz_yzz_yzz_0, ta2_yz_yzz_zzz_0, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzz_xxx_0[i] = ta1_z_zz_xxx_1[i] + ta2_yz_zz_xxx_0[i] * pa_y[i] - ta2_yz_zz_xxx_1[i] * pc_y[i];

        ta2_yz_yzz_xxy_0[i] = ta2_yz_zz_xx_0[i] * fe_0 - ta2_yz_zz_xx_1[i] * fe_0 + ta1_z_zz_xxy_1[i] + ta2_yz_zz_xxy_0[i] * pa_y[i] - ta2_yz_zz_xxy_1[i] * pc_y[i];

        ta2_yz_yzz_xxz_0[i] = ta1_z_zz_xxz_1[i] + ta2_yz_zz_xxz_0[i] * pa_y[i] - ta2_yz_zz_xxz_1[i] * pc_y[i];

        ta2_yz_yzz_xyy_0[i] = 2.0 * ta2_yz_zz_xy_0[i] * fe_0 - 2.0 * ta2_yz_zz_xy_1[i] * fe_0 + ta1_z_zz_xyy_1[i] + ta2_yz_zz_xyy_0[i] * pa_y[i] - ta2_yz_zz_xyy_1[i] * pc_y[i];

        ta2_yz_yzz_xyz_0[i] = ta2_yz_zz_xz_0[i] * fe_0 - ta2_yz_zz_xz_1[i] * fe_0 + ta1_z_zz_xyz_1[i] + ta2_yz_zz_xyz_0[i] * pa_y[i] - ta2_yz_zz_xyz_1[i] * pc_y[i];

        ta2_yz_yzz_xzz_0[i] = ta1_z_zz_xzz_1[i] + ta2_yz_zz_xzz_0[i] * pa_y[i] - ta2_yz_zz_xzz_1[i] * pc_y[i];

        ta2_yz_yzz_yyy_0[i] = 3.0 * ta2_yz_zz_yy_0[i] * fe_0 - 3.0 * ta2_yz_zz_yy_1[i] * fe_0 + ta1_z_zz_yyy_1[i] + ta2_yz_zz_yyy_0[i] * pa_y[i] - ta2_yz_zz_yyy_1[i] * pc_y[i];

        ta2_yz_yzz_yyz_0[i] = 2.0 * ta2_yz_zz_yz_0[i] * fe_0 - 2.0 * ta2_yz_zz_yz_1[i] * fe_0 + ta1_z_zz_yyz_1[i] + ta2_yz_zz_yyz_0[i] * pa_y[i] - ta2_yz_zz_yyz_1[i] * pc_y[i];

        ta2_yz_yzz_yzz_0[i] = ta2_yz_zz_zz_0[i] * fe_0 - ta2_yz_zz_zz_1[i] * fe_0 + ta1_z_zz_yzz_1[i] + ta2_yz_zz_yzz_0[i] * pa_y[i] - ta2_yz_zz_yzz_1[i] * pc_y[i];

        ta2_yz_yzz_zzz_0[i] = ta1_z_zz_zzz_1[i] + ta2_yz_zz_zzz_0[i] * pa_y[i] - ta2_yz_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 490-500 components of targeted buffer : FF

    auto ta2_yz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 490);

    auto ta2_yz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 491);

    auto ta2_yz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 492);

    auto ta2_yz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 493);

    auto ta2_yz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 494);

    auto ta2_yz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 495);

    auto ta2_yz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 496);

    auto ta2_yz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 497);

    auto ta2_yz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 498);

    auto ta2_yz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 499);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zz_xxx_1, ta1_y_zz_xxy_1, ta1_y_zz_xxz_1, ta1_y_zz_xyy_1, ta1_y_zz_xyz_1, ta1_y_zz_xzz_1, ta1_y_zz_yyy_1, ta1_y_zz_yyz_1, ta1_y_zz_yzz_1, ta1_y_zz_zzz_1, ta2_yz_z_xxx_0, ta2_yz_z_xxx_1, ta2_yz_z_xxy_0, ta2_yz_z_xxy_1, ta2_yz_z_xxz_0, ta2_yz_z_xxz_1, ta2_yz_z_xyy_0, ta2_yz_z_xyy_1, ta2_yz_z_xyz_0, ta2_yz_z_xyz_1, ta2_yz_z_xzz_0, ta2_yz_z_xzz_1, ta2_yz_z_yyy_0, ta2_yz_z_yyy_1, ta2_yz_z_yyz_0, ta2_yz_z_yyz_1, ta2_yz_z_yzz_0, ta2_yz_z_yzz_1, ta2_yz_z_zzz_0, ta2_yz_z_zzz_1, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, ta2_yz_zzz_xxx_0, ta2_yz_zzz_xxy_0, ta2_yz_zzz_xxz_0, ta2_yz_zzz_xyy_0, ta2_yz_zzz_xyz_0, ta2_yz_zzz_xzz_0, ta2_yz_zzz_yyy_0, ta2_yz_zzz_yyz_0, ta2_yz_zzz_yzz_0, ta2_yz_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzz_xxx_0[i] = 2.0 * ta2_yz_z_xxx_0[i] * fe_0 - 2.0 * ta2_yz_z_xxx_1[i] * fe_0 + ta1_y_zz_xxx_1[i] + ta2_yz_zz_xxx_0[i] * pa_z[i] - ta2_yz_zz_xxx_1[i] * pc_z[i];

        ta2_yz_zzz_xxy_0[i] = 2.0 * ta2_yz_z_xxy_0[i] * fe_0 - 2.0 * ta2_yz_z_xxy_1[i] * fe_0 + ta1_y_zz_xxy_1[i] + ta2_yz_zz_xxy_0[i] * pa_z[i] - ta2_yz_zz_xxy_1[i] * pc_z[i];

        ta2_yz_zzz_xxz_0[i] = 2.0 * ta2_yz_z_xxz_0[i] * fe_0 - 2.0 * ta2_yz_z_xxz_1[i] * fe_0 + ta2_yz_zz_xx_0[i] * fe_0 - ta2_yz_zz_xx_1[i] * fe_0 + ta1_y_zz_xxz_1[i] + ta2_yz_zz_xxz_0[i] * pa_z[i] - ta2_yz_zz_xxz_1[i] * pc_z[i];

        ta2_yz_zzz_xyy_0[i] = 2.0 * ta2_yz_z_xyy_0[i] * fe_0 - 2.0 * ta2_yz_z_xyy_1[i] * fe_0 + ta1_y_zz_xyy_1[i] + ta2_yz_zz_xyy_0[i] * pa_z[i] - ta2_yz_zz_xyy_1[i] * pc_z[i];

        ta2_yz_zzz_xyz_0[i] = 2.0 * ta2_yz_z_xyz_0[i] * fe_0 - 2.0 * ta2_yz_z_xyz_1[i] * fe_0 + ta2_yz_zz_xy_0[i] * fe_0 - ta2_yz_zz_xy_1[i] * fe_0 + ta1_y_zz_xyz_1[i] + ta2_yz_zz_xyz_0[i] * pa_z[i] - ta2_yz_zz_xyz_1[i] * pc_z[i];

        ta2_yz_zzz_xzz_0[i] = 2.0 * ta2_yz_z_xzz_0[i] * fe_0 - 2.0 * ta2_yz_z_xzz_1[i] * fe_0 + 2.0 * ta2_yz_zz_xz_0[i] * fe_0 - 2.0 * ta2_yz_zz_xz_1[i] * fe_0 + ta1_y_zz_xzz_1[i] + ta2_yz_zz_xzz_0[i] * pa_z[i] - ta2_yz_zz_xzz_1[i] * pc_z[i];

        ta2_yz_zzz_yyy_0[i] = 2.0 * ta2_yz_z_yyy_0[i] * fe_0 - 2.0 * ta2_yz_z_yyy_1[i] * fe_0 + ta1_y_zz_yyy_1[i] + ta2_yz_zz_yyy_0[i] * pa_z[i] - ta2_yz_zz_yyy_1[i] * pc_z[i];

        ta2_yz_zzz_yyz_0[i] = 2.0 * ta2_yz_z_yyz_0[i] * fe_0 - 2.0 * ta2_yz_z_yyz_1[i] * fe_0 + ta2_yz_zz_yy_0[i] * fe_0 - ta2_yz_zz_yy_1[i] * fe_0 + ta1_y_zz_yyz_1[i] + ta2_yz_zz_yyz_0[i] * pa_z[i] - ta2_yz_zz_yyz_1[i] * pc_z[i];

        ta2_yz_zzz_yzz_0[i] = 2.0 * ta2_yz_z_yzz_0[i] * fe_0 - 2.0 * ta2_yz_z_yzz_1[i] * fe_0 + 2.0 * ta2_yz_zz_yz_0[i] * fe_0 - 2.0 * ta2_yz_zz_yz_1[i] * fe_0 + ta1_y_zz_yzz_1[i] + ta2_yz_zz_yzz_0[i] * pa_z[i] - ta2_yz_zz_yzz_1[i] * pc_z[i];

        ta2_yz_zzz_zzz_0[i] = 2.0 * ta2_yz_z_zzz_0[i] * fe_0 - 2.0 * ta2_yz_z_zzz_1[i] * fe_0 + 3.0 * ta2_yz_zz_zz_0[i] * fe_0 - 3.0 * ta2_yz_zz_zz_1[i] * fe_0 + ta1_y_zz_zzz_1[i] + ta2_yz_zz_zzz_0[i] * pa_z[i] - ta2_yz_zz_zzz_1[i] * pc_z[i];
    }

    // Set up 500-510 components of targeted buffer : FF

    auto ta2_zz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 500);

    auto ta2_zz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 501);

    auto ta2_zz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 502);

    auto ta2_zz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 503);

    auto ta2_zz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 504);

    auto ta2_zz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 505);

    auto ta2_zz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 506);

    auto ta2_zz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 507);

    auto ta2_zz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 508);

    auto ta2_zz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 509);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_x_xxx_0, ta2_zz_x_xxx_1, ta2_zz_x_xxy_0, ta2_zz_x_xxy_1, ta2_zz_x_xxz_0, ta2_zz_x_xxz_1, ta2_zz_x_xyy_0, ta2_zz_x_xyy_1, ta2_zz_x_xyz_0, ta2_zz_x_xyz_1, ta2_zz_x_xzz_0, ta2_zz_x_xzz_1, ta2_zz_x_yyy_0, ta2_zz_x_yyy_1, ta2_zz_x_yyz_0, ta2_zz_x_yyz_1, ta2_zz_x_yzz_0, ta2_zz_x_yzz_1, ta2_zz_x_zzz_0, ta2_zz_x_zzz_1, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_yy_0, ta2_zz_xx_yy_1, ta2_zz_xx_yyy_0, ta2_zz_xx_yyy_1, ta2_zz_xx_yyz_0, ta2_zz_xx_yyz_1, ta2_zz_xx_yz_0, ta2_zz_xx_yz_1, ta2_zz_xx_yzz_0, ta2_zz_xx_yzz_1, ta2_zz_xx_zz_0, ta2_zz_xx_zz_1, ta2_zz_xx_zzz_0, ta2_zz_xx_zzz_1, ta2_zz_xxx_xxx_0, ta2_zz_xxx_xxy_0, ta2_zz_xxx_xxz_0, ta2_zz_xxx_xyy_0, ta2_zz_xxx_xyz_0, ta2_zz_xxx_xzz_0, ta2_zz_xxx_yyy_0, ta2_zz_xxx_yyz_0, ta2_zz_xxx_yzz_0, ta2_zz_xxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxx_xxx_0[i] = 2.0 * ta2_zz_x_xxx_0[i] * fe_0 - 2.0 * ta2_zz_x_xxx_1[i] * fe_0 + 3.0 * ta2_zz_xx_xx_0[i] * fe_0 - 3.0 * ta2_zz_xx_xx_1[i] * fe_0 + ta2_zz_xx_xxx_0[i] * pa_x[i] - ta2_zz_xx_xxx_1[i] * pc_x[i];

        ta2_zz_xxx_xxy_0[i] = 2.0 * ta2_zz_x_xxy_0[i] * fe_0 - 2.0 * ta2_zz_x_xxy_1[i] * fe_0 + 2.0 * ta2_zz_xx_xy_0[i] * fe_0 - 2.0 * ta2_zz_xx_xy_1[i] * fe_0 + ta2_zz_xx_xxy_0[i] * pa_x[i] - ta2_zz_xx_xxy_1[i] * pc_x[i];

        ta2_zz_xxx_xxz_0[i] = 2.0 * ta2_zz_x_xxz_0[i] * fe_0 - 2.0 * ta2_zz_x_xxz_1[i] * fe_0 + 2.0 * ta2_zz_xx_xz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xz_1[i] * fe_0 + ta2_zz_xx_xxz_0[i] * pa_x[i] - ta2_zz_xx_xxz_1[i] * pc_x[i];

        ta2_zz_xxx_xyy_0[i] = 2.0 * ta2_zz_x_xyy_0[i] * fe_0 - 2.0 * ta2_zz_x_xyy_1[i] * fe_0 + ta2_zz_xx_yy_0[i] * fe_0 - ta2_zz_xx_yy_1[i] * fe_0 + ta2_zz_xx_xyy_0[i] * pa_x[i] - ta2_zz_xx_xyy_1[i] * pc_x[i];

        ta2_zz_xxx_xyz_0[i] = 2.0 * ta2_zz_x_xyz_0[i] * fe_0 - 2.0 * ta2_zz_x_xyz_1[i] * fe_0 + ta2_zz_xx_yz_0[i] * fe_0 - ta2_zz_xx_yz_1[i] * fe_0 + ta2_zz_xx_xyz_0[i] * pa_x[i] - ta2_zz_xx_xyz_1[i] * pc_x[i];

        ta2_zz_xxx_xzz_0[i] = 2.0 * ta2_zz_x_xzz_0[i] * fe_0 - 2.0 * ta2_zz_x_xzz_1[i] * fe_0 + ta2_zz_xx_zz_0[i] * fe_0 - ta2_zz_xx_zz_1[i] * fe_0 + ta2_zz_xx_xzz_0[i] * pa_x[i] - ta2_zz_xx_xzz_1[i] * pc_x[i];

        ta2_zz_xxx_yyy_0[i] = 2.0 * ta2_zz_x_yyy_0[i] * fe_0 - 2.0 * ta2_zz_x_yyy_1[i] * fe_0 + ta2_zz_xx_yyy_0[i] * pa_x[i] - ta2_zz_xx_yyy_1[i] * pc_x[i];

        ta2_zz_xxx_yyz_0[i] = 2.0 * ta2_zz_x_yyz_0[i] * fe_0 - 2.0 * ta2_zz_x_yyz_1[i] * fe_0 + ta2_zz_xx_yyz_0[i] * pa_x[i] - ta2_zz_xx_yyz_1[i] * pc_x[i];

        ta2_zz_xxx_yzz_0[i] = 2.0 * ta2_zz_x_yzz_0[i] * fe_0 - 2.0 * ta2_zz_x_yzz_1[i] * fe_0 + ta2_zz_xx_yzz_0[i] * pa_x[i] - ta2_zz_xx_yzz_1[i] * pc_x[i];

        ta2_zz_xxx_zzz_0[i] = 2.0 * ta2_zz_x_zzz_0[i] * fe_0 - 2.0 * ta2_zz_x_zzz_1[i] * fe_0 + ta2_zz_xx_zzz_0[i] * pa_x[i] - ta2_zz_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 510-520 components of targeted buffer : FF

    auto ta2_zz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 510);

    auto ta2_zz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 511);

    auto ta2_zz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 512);

    auto ta2_zz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 513);

    auto ta2_zz_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 514);

    auto ta2_zz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 515);

    auto ta2_zz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 516);

    auto ta2_zz_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 517);

    auto ta2_zz_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 518);

    auto ta2_zz_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 519);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_zzz_0, ta2_zz_xx_zzz_1, ta2_zz_xxy_xxx_0, ta2_zz_xxy_xxy_0, ta2_zz_xxy_xxz_0, ta2_zz_xxy_xyy_0, ta2_zz_xxy_xyz_0, ta2_zz_xxy_xzz_0, ta2_zz_xxy_yyy_0, ta2_zz_xxy_yyz_0, ta2_zz_xxy_yzz_0, ta2_zz_xxy_zzz_0, ta2_zz_xy_yyy_0, ta2_zz_xy_yyy_1, ta2_zz_xy_yyz_0, ta2_zz_xy_yyz_1, ta2_zz_xy_yzz_0, ta2_zz_xy_yzz_1, ta2_zz_y_yyy_0, ta2_zz_y_yyy_1, ta2_zz_y_yyz_0, ta2_zz_y_yyz_1, ta2_zz_y_yzz_0, ta2_zz_y_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxy_xxx_0[i] = ta2_zz_xx_xxx_0[i] * pa_y[i] - ta2_zz_xx_xxx_1[i] * pc_y[i];

        ta2_zz_xxy_xxy_0[i] = ta2_zz_xx_xx_0[i] * fe_0 - ta2_zz_xx_xx_1[i] * fe_0 + ta2_zz_xx_xxy_0[i] * pa_y[i] - ta2_zz_xx_xxy_1[i] * pc_y[i];

        ta2_zz_xxy_xxz_0[i] = ta2_zz_xx_xxz_0[i] * pa_y[i] - ta2_zz_xx_xxz_1[i] * pc_y[i];

        ta2_zz_xxy_xyy_0[i] = 2.0 * ta2_zz_xx_xy_0[i] * fe_0 - 2.0 * ta2_zz_xx_xy_1[i] * fe_0 + ta2_zz_xx_xyy_0[i] * pa_y[i] - ta2_zz_xx_xyy_1[i] * pc_y[i];

        ta2_zz_xxy_xyz_0[i] = ta2_zz_xx_xz_0[i] * fe_0 - ta2_zz_xx_xz_1[i] * fe_0 + ta2_zz_xx_xyz_0[i] * pa_y[i] - ta2_zz_xx_xyz_1[i] * pc_y[i];

        ta2_zz_xxy_xzz_0[i] = ta2_zz_xx_xzz_0[i] * pa_y[i] - ta2_zz_xx_xzz_1[i] * pc_y[i];

        ta2_zz_xxy_yyy_0[i] = ta2_zz_y_yyy_0[i] * fe_0 - ta2_zz_y_yyy_1[i] * fe_0 + ta2_zz_xy_yyy_0[i] * pa_x[i] - ta2_zz_xy_yyy_1[i] * pc_x[i];

        ta2_zz_xxy_yyz_0[i] = ta2_zz_y_yyz_0[i] * fe_0 - ta2_zz_y_yyz_1[i] * fe_0 + ta2_zz_xy_yyz_0[i] * pa_x[i] - ta2_zz_xy_yyz_1[i] * pc_x[i];

        ta2_zz_xxy_yzz_0[i] = ta2_zz_y_yzz_0[i] * fe_0 - ta2_zz_y_yzz_1[i] * fe_0 + ta2_zz_xy_yzz_0[i] * pa_x[i] - ta2_zz_xy_yzz_1[i] * pc_x[i];

        ta2_zz_xxy_zzz_0[i] = ta2_zz_xx_zzz_0[i] * pa_y[i] - ta2_zz_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 520-530 components of targeted buffer : FF

    auto ta2_zz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 520);

    auto ta2_zz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 521);

    auto ta2_zz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 522);

    auto ta2_zz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 523);

    auto ta2_zz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 524);

    auto ta2_zz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 525);

    auto ta2_zz_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 526);

    auto ta2_zz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 527);

    auto ta2_zz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 528);

    auto ta2_zz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 529);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xx_xxx_1, ta1_z_xx_xxy_1, ta1_z_xx_xxz_1, ta1_z_xx_xyy_1, ta1_z_xx_xyz_1, ta1_z_xx_xzz_1, ta1_z_xx_yyy_1, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_yyy_0, ta2_zz_xx_yyy_1, ta2_zz_xxz_xxx_0, ta2_zz_xxz_xxy_0, ta2_zz_xxz_xxz_0, ta2_zz_xxz_xyy_0, ta2_zz_xxz_xyz_0, ta2_zz_xxz_xzz_0, ta2_zz_xxz_yyy_0, ta2_zz_xxz_yyz_0, ta2_zz_xxz_yzz_0, ta2_zz_xxz_zzz_0, ta2_zz_xz_yyz_0, ta2_zz_xz_yyz_1, ta2_zz_xz_yzz_0, ta2_zz_xz_yzz_1, ta2_zz_xz_zzz_0, ta2_zz_xz_zzz_1, ta2_zz_z_yyz_0, ta2_zz_z_yyz_1, ta2_zz_z_yzz_0, ta2_zz_z_yzz_1, ta2_zz_z_zzz_0, ta2_zz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxz_xxx_0[i] = 2.0 * ta1_z_xx_xxx_1[i] + ta2_zz_xx_xxx_0[i] * pa_z[i] - ta2_zz_xx_xxx_1[i] * pc_z[i];

        ta2_zz_xxz_xxy_0[i] = 2.0 * ta1_z_xx_xxy_1[i] + ta2_zz_xx_xxy_0[i] * pa_z[i] - ta2_zz_xx_xxy_1[i] * pc_z[i];

        ta2_zz_xxz_xxz_0[i] = ta2_zz_xx_xx_0[i] * fe_0 - ta2_zz_xx_xx_1[i] * fe_0 + 2.0 * ta1_z_xx_xxz_1[i] + ta2_zz_xx_xxz_0[i] * pa_z[i] - ta2_zz_xx_xxz_1[i] * pc_z[i];

        ta2_zz_xxz_xyy_0[i] = 2.0 * ta1_z_xx_xyy_1[i] + ta2_zz_xx_xyy_0[i] * pa_z[i] - ta2_zz_xx_xyy_1[i] * pc_z[i];

        ta2_zz_xxz_xyz_0[i] = ta2_zz_xx_xy_0[i] * fe_0 - ta2_zz_xx_xy_1[i] * fe_0 + 2.0 * ta1_z_xx_xyz_1[i] + ta2_zz_xx_xyz_0[i] * pa_z[i] - ta2_zz_xx_xyz_1[i] * pc_z[i];

        ta2_zz_xxz_xzz_0[i] = 2.0 * ta2_zz_xx_xz_0[i] * fe_0 - 2.0 * ta2_zz_xx_xz_1[i] * fe_0 + 2.0 * ta1_z_xx_xzz_1[i] + ta2_zz_xx_xzz_0[i] * pa_z[i] - ta2_zz_xx_xzz_1[i] * pc_z[i];

        ta2_zz_xxz_yyy_0[i] = 2.0 * ta1_z_xx_yyy_1[i] + ta2_zz_xx_yyy_0[i] * pa_z[i] - ta2_zz_xx_yyy_1[i] * pc_z[i];

        ta2_zz_xxz_yyz_0[i] = ta2_zz_z_yyz_0[i] * fe_0 - ta2_zz_z_yyz_1[i] * fe_0 + ta2_zz_xz_yyz_0[i] * pa_x[i] - ta2_zz_xz_yyz_1[i] * pc_x[i];

        ta2_zz_xxz_yzz_0[i] = ta2_zz_z_yzz_0[i] * fe_0 - ta2_zz_z_yzz_1[i] * fe_0 + ta2_zz_xz_yzz_0[i] * pa_x[i] - ta2_zz_xz_yzz_1[i] * pc_x[i];

        ta2_zz_xxz_zzz_0[i] = ta2_zz_z_zzz_0[i] * fe_0 - ta2_zz_z_zzz_1[i] * fe_0 + ta2_zz_xz_zzz_0[i] * pa_x[i] - ta2_zz_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 530-540 components of targeted buffer : FF

    auto ta2_zz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 530);

    auto ta2_zz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 531);

    auto ta2_zz_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 532);

    auto ta2_zz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 533);

    auto ta2_zz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 534);

    auto ta2_zz_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 535);

    auto ta2_zz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 536);

    auto ta2_zz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 537);

    auto ta2_zz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 538);

    auto ta2_zz_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 539);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xyy_xxx_0, ta2_zz_xyy_xxy_0, ta2_zz_xyy_xxz_0, ta2_zz_xyy_xyy_0, ta2_zz_xyy_xyz_0, ta2_zz_xyy_xzz_0, ta2_zz_xyy_yyy_0, ta2_zz_xyy_yyz_0, ta2_zz_xyy_yzz_0, ta2_zz_xyy_zzz_0, ta2_zz_yy_xx_0, ta2_zz_yy_xx_1, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxz_0, ta2_zz_yy_xxz_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xz_0, ta2_zz_yy_xz_1, ta2_zz_yy_xzz_0, ta2_zz_yy_xzz_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_zz_0, ta2_zz_yy_zz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyy_xxx_0[i] = 3.0 * ta2_zz_yy_xx_0[i] * fe_0 - 3.0 * ta2_zz_yy_xx_1[i] * fe_0 + ta2_zz_yy_xxx_0[i] * pa_x[i] - ta2_zz_yy_xxx_1[i] * pc_x[i];

        ta2_zz_xyy_xxy_0[i] = 2.0 * ta2_zz_yy_xy_0[i] * fe_0 - 2.0 * ta2_zz_yy_xy_1[i] * fe_0 + ta2_zz_yy_xxy_0[i] * pa_x[i] - ta2_zz_yy_xxy_1[i] * pc_x[i];

        ta2_zz_xyy_xxz_0[i] = 2.0 * ta2_zz_yy_xz_0[i] * fe_0 - 2.0 * ta2_zz_yy_xz_1[i] * fe_0 + ta2_zz_yy_xxz_0[i] * pa_x[i] - ta2_zz_yy_xxz_1[i] * pc_x[i];

        ta2_zz_xyy_xyy_0[i] = ta2_zz_yy_yy_0[i] * fe_0 - ta2_zz_yy_yy_1[i] * fe_0 + ta2_zz_yy_xyy_0[i] * pa_x[i] - ta2_zz_yy_xyy_1[i] * pc_x[i];

        ta2_zz_xyy_xyz_0[i] = ta2_zz_yy_yz_0[i] * fe_0 - ta2_zz_yy_yz_1[i] * fe_0 + ta2_zz_yy_xyz_0[i] * pa_x[i] - ta2_zz_yy_xyz_1[i] * pc_x[i];

        ta2_zz_xyy_xzz_0[i] = ta2_zz_yy_zz_0[i] * fe_0 - ta2_zz_yy_zz_1[i] * fe_0 + ta2_zz_yy_xzz_0[i] * pa_x[i] - ta2_zz_yy_xzz_1[i] * pc_x[i];

        ta2_zz_xyy_yyy_0[i] = ta2_zz_yy_yyy_0[i] * pa_x[i] - ta2_zz_yy_yyy_1[i] * pc_x[i];

        ta2_zz_xyy_yyz_0[i] = ta2_zz_yy_yyz_0[i] * pa_x[i] - ta2_zz_yy_yyz_1[i] * pc_x[i];

        ta2_zz_xyy_yzz_0[i] = ta2_zz_yy_yzz_0[i] * pa_x[i] - ta2_zz_yy_yzz_1[i] * pc_x[i];

        ta2_zz_xyy_zzz_0[i] = ta2_zz_yy_zzz_0[i] * pa_x[i] - ta2_zz_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 540-550 components of targeted buffer : FF

    auto ta2_zz_xyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 540);

    auto ta2_zz_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 541);

    auto ta2_zz_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 542);

    auto ta2_zz_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 543);

    auto ta2_zz_xyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 544);

    auto ta2_zz_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 545);

    auto ta2_zz_xyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 546);

    auto ta2_zz_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 547);

    auto ta2_zz_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 548);

    auto ta2_zz_xyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 549);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xy_xxy_1, ta1_z_xy_xyy_1, ta2_zz_xy_xxy_0, ta2_zz_xy_xxy_1, ta2_zz_xy_xyy_0, ta2_zz_xy_xyy_1, ta2_zz_xyz_xxx_0, ta2_zz_xyz_xxy_0, ta2_zz_xyz_xxz_0, ta2_zz_xyz_xyy_0, ta2_zz_xyz_xyz_0, ta2_zz_xyz_xzz_0, ta2_zz_xyz_yyy_0, ta2_zz_xyz_yyz_0, ta2_zz_xyz_yzz_0, ta2_zz_xyz_zzz_0, ta2_zz_xz_xxx_0, ta2_zz_xz_xxx_1, ta2_zz_xz_xxz_0, ta2_zz_xz_xxz_1, ta2_zz_xz_xzz_0, ta2_zz_xz_xzz_1, ta2_zz_yz_xyz_0, ta2_zz_yz_xyz_1, ta2_zz_yz_yyy_0, ta2_zz_yz_yyy_1, ta2_zz_yz_yyz_0, ta2_zz_yz_yyz_1, ta2_zz_yz_yz_0, ta2_zz_yz_yz_1, ta2_zz_yz_yzz_0, ta2_zz_yz_yzz_1, ta2_zz_yz_zzz_0, ta2_zz_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyz_xxx_0[i] = ta2_zz_xz_xxx_0[i] * pa_y[i] - ta2_zz_xz_xxx_1[i] * pc_y[i];

        ta2_zz_xyz_xxy_0[i] = 2.0 * ta1_z_xy_xxy_1[i] + ta2_zz_xy_xxy_0[i] * pa_z[i] - ta2_zz_xy_xxy_1[i] * pc_z[i];

        ta2_zz_xyz_xxz_0[i] = ta2_zz_xz_xxz_0[i] * pa_y[i] - ta2_zz_xz_xxz_1[i] * pc_y[i];

        ta2_zz_xyz_xyy_0[i] = 2.0 * ta1_z_xy_xyy_1[i] + ta2_zz_xy_xyy_0[i] * pa_z[i] - ta2_zz_xy_xyy_1[i] * pc_z[i];

        ta2_zz_xyz_xyz_0[i] = ta2_zz_yz_yz_0[i] * fe_0 - ta2_zz_yz_yz_1[i] * fe_0 + ta2_zz_yz_xyz_0[i] * pa_x[i] - ta2_zz_yz_xyz_1[i] * pc_x[i];

        ta2_zz_xyz_xzz_0[i] = ta2_zz_xz_xzz_0[i] * pa_y[i] - ta2_zz_xz_xzz_1[i] * pc_y[i];

        ta2_zz_xyz_yyy_0[i] = ta2_zz_yz_yyy_0[i] * pa_x[i] - ta2_zz_yz_yyy_1[i] * pc_x[i];

        ta2_zz_xyz_yyz_0[i] = ta2_zz_yz_yyz_0[i] * pa_x[i] - ta2_zz_yz_yyz_1[i] * pc_x[i];

        ta2_zz_xyz_yzz_0[i] = ta2_zz_yz_yzz_0[i] * pa_x[i] - ta2_zz_yz_yzz_1[i] * pc_x[i];

        ta2_zz_xyz_zzz_0[i] = ta2_zz_yz_zzz_0[i] * pa_x[i] - ta2_zz_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 550-560 components of targeted buffer : FF

    auto ta2_zz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 550);

    auto ta2_zz_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 551);

    auto ta2_zz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 552);

    auto ta2_zz_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 553);

    auto ta2_zz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 554);

    auto ta2_zz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 555);

    auto ta2_zz_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 556);

    auto ta2_zz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 557);

    auto ta2_zz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 558);

    auto ta2_zz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 559);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xzz_xxx_0, ta2_zz_xzz_xxy_0, ta2_zz_xzz_xxz_0, ta2_zz_xzz_xyy_0, ta2_zz_xzz_xyz_0, ta2_zz_xzz_xzz_0, ta2_zz_xzz_yyy_0, ta2_zz_xzz_yyz_0, ta2_zz_xzz_yzz_0, ta2_zz_xzz_zzz_0, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzz_xxx_0[i] = 3.0 * ta2_zz_zz_xx_0[i] * fe_0 - 3.0 * ta2_zz_zz_xx_1[i] * fe_0 + ta2_zz_zz_xxx_0[i] * pa_x[i] - ta2_zz_zz_xxx_1[i] * pc_x[i];

        ta2_zz_xzz_xxy_0[i] = 2.0 * ta2_zz_zz_xy_0[i] * fe_0 - 2.0 * ta2_zz_zz_xy_1[i] * fe_0 + ta2_zz_zz_xxy_0[i] * pa_x[i] - ta2_zz_zz_xxy_1[i] * pc_x[i];

        ta2_zz_xzz_xxz_0[i] = 2.0 * ta2_zz_zz_xz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xz_1[i] * fe_0 + ta2_zz_zz_xxz_0[i] * pa_x[i] - ta2_zz_zz_xxz_1[i] * pc_x[i];

        ta2_zz_xzz_xyy_0[i] = ta2_zz_zz_yy_0[i] * fe_0 - ta2_zz_zz_yy_1[i] * fe_0 + ta2_zz_zz_xyy_0[i] * pa_x[i] - ta2_zz_zz_xyy_1[i] * pc_x[i];

        ta2_zz_xzz_xyz_0[i] = ta2_zz_zz_yz_0[i] * fe_0 - ta2_zz_zz_yz_1[i] * fe_0 + ta2_zz_zz_xyz_0[i] * pa_x[i] - ta2_zz_zz_xyz_1[i] * pc_x[i];

        ta2_zz_xzz_xzz_0[i] = ta2_zz_zz_zz_0[i] * fe_0 - ta2_zz_zz_zz_1[i] * fe_0 + ta2_zz_zz_xzz_0[i] * pa_x[i] - ta2_zz_zz_xzz_1[i] * pc_x[i];

        ta2_zz_xzz_yyy_0[i] = ta2_zz_zz_yyy_0[i] * pa_x[i] - ta2_zz_zz_yyy_1[i] * pc_x[i];

        ta2_zz_xzz_yyz_0[i] = ta2_zz_zz_yyz_0[i] * pa_x[i] - ta2_zz_zz_yyz_1[i] * pc_x[i];

        ta2_zz_xzz_yzz_0[i] = ta2_zz_zz_yzz_0[i] * pa_x[i] - ta2_zz_zz_yzz_1[i] * pc_x[i];

        ta2_zz_xzz_zzz_0[i] = ta2_zz_zz_zzz_0[i] * pa_x[i] - ta2_zz_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 560-570 components of targeted buffer : FF

    auto ta2_zz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 560);

    auto ta2_zz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 561);

    auto ta2_zz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 562);

    auto ta2_zz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 563);

    auto ta2_zz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 564);

    auto ta2_zz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 565);

    auto ta2_zz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 566);

    auto ta2_zz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 567);

    auto ta2_zz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 568);

    auto ta2_zz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 569);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_y_xxx_0, ta2_zz_y_xxx_1, ta2_zz_y_xxy_0, ta2_zz_y_xxy_1, ta2_zz_y_xxz_0, ta2_zz_y_xxz_1, ta2_zz_y_xyy_0, ta2_zz_y_xyy_1, ta2_zz_y_xyz_0, ta2_zz_y_xyz_1, ta2_zz_y_xzz_0, ta2_zz_y_xzz_1, ta2_zz_y_yyy_0, ta2_zz_y_yyy_1, ta2_zz_y_yyz_0, ta2_zz_y_yyz_1, ta2_zz_y_yzz_0, ta2_zz_y_yzz_1, ta2_zz_y_zzz_0, ta2_zz_y_zzz_1, ta2_zz_yy_xx_0, ta2_zz_yy_xx_1, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxz_0, ta2_zz_yy_xxz_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xz_0, ta2_zz_yy_xz_1, ta2_zz_yy_xzz_0, ta2_zz_yy_xzz_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_zz_0, ta2_zz_yy_zz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, ta2_zz_yyy_xxx_0, ta2_zz_yyy_xxy_0, ta2_zz_yyy_xxz_0, ta2_zz_yyy_xyy_0, ta2_zz_yyy_xyz_0, ta2_zz_yyy_xzz_0, ta2_zz_yyy_yyy_0, ta2_zz_yyy_yyz_0, ta2_zz_yyy_yzz_0, ta2_zz_yyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyy_xxx_0[i] = 2.0 * ta2_zz_y_xxx_0[i] * fe_0 - 2.0 * ta2_zz_y_xxx_1[i] * fe_0 + ta2_zz_yy_xxx_0[i] * pa_y[i] - ta2_zz_yy_xxx_1[i] * pc_y[i];

        ta2_zz_yyy_xxy_0[i] = 2.0 * ta2_zz_y_xxy_0[i] * fe_0 - 2.0 * ta2_zz_y_xxy_1[i] * fe_0 + ta2_zz_yy_xx_0[i] * fe_0 - ta2_zz_yy_xx_1[i] * fe_0 + ta2_zz_yy_xxy_0[i] * pa_y[i] - ta2_zz_yy_xxy_1[i] * pc_y[i];

        ta2_zz_yyy_xxz_0[i] = 2.0 * ta2_zz_y_xxz_0[i] * fe_0 - 2.0 * ta2_zz_y_xxz_1[i] * fe_0 + ta2_zz_yy_xxz_0[i] * pa_y[i] - ta2_zz_yy_xxz_1[i] * pc_y[i];

        ta2_zz_yyy_xyy_0[i] = 2.0 * ta2_zz_y_xyy_0[i] * fe_0 - 2.0 * ta2_zz_y_xyy_1[i] * fe_0 + 2.0 * ta2_zz_yy_xy_0[i] * fe_0 - 2.0 * ta2_zz_yy_xy_1[i] * fe_0 + ta2_zz_yy_xyy_0[i] * pa_y[i] - ta2_zz_yy_xyy_1[i] * pc_y[i];

        ta2_zz_yyy_xyz_0[i] = 2.0 * ta2_zz_y_xyz_0[i] * fe_0 - 2.0 * ta2_zz_y_xyz_1[i] * fe_0 + ta2_zz_yy_xz_0[i] * fe_0 - ta2_zz_yy_xz_1[i] * fe_0 + ta2_zz_yy_xyz_0[i] * pa_y[i] - ta2_zz_yy_xyz_1[i] * pc_y[i];

        ta2_zz_yyy_xzz_0[i] = 2.0 * ta2_zz_y_xzz_0[i] * fe_0 - 2.0 * ta2_zz_y_xzz_1[i] * fe_0 + ta2_zz_yy_xzz_0[i] * pa_y[i] - ta2_zz_yy_xzz_1[i] * pc_y[i];

        ta2_zz_yyy_yyy_0[i] = 2.0 * ta2_zz_y_yyy_0[i] * fe_0 - 2.0 * ta2_zz_y_yyy_1[i] * fe_0 + 3.0 * ta2_zz_yy_yy_0[i] * fe_0 - 3.0 * ta2_zz_yy_yy_1[i] * fe_0 + ta2_zz_yy_yyy_0[i] * pa_y[i] - ta2_zz_yy_yyy_1[i] * pc_y[i];

        ta2_zz_yyy_yyz_0[i] = 2.0 * ta2_zz_y_yyz_0[i] * fe_0 - 2.0 * ta2_zz_y_yyz_1[i] * fe_0 + 2.0 * ta2_zz_yy_yz_0[i] * fe_0 - 2.0 * ta2_zz_yy_yz_1[i] * fe_0 + ta2_zz_yy_yyz_0[i] * pa_y[i] - ta2_zz_yy_yyz_1[i] * pc_y[i];

        ta2_zz_yyy_yzz_0[i] = 2.0 * ta2_zz_y_yzz_0[i] * fe_0 - 2.0 * ta2_zz_y_yzz_1[i] * fe_0 + ta2_zz_yy_zz_0[i] * fe_0 - ta2_zz_yy_zz_1[i] * fe_0 + ta2_zz_yy_yzz_0[i] * pa_y[i] - ta2_zz_yy_yzz_1[i] * pc_y[i];

        ta2_zz_yyy_zzz_0[i] = 2.0 * ta2_zz_y_zzz_0[i] * fe_0 - 2.0 * ta2_zz_y_zzz_1[i] * fe_0 + ta2_zz_yy_zzz_0[i] * pa_y[i] - ta2_zz_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 570-580 components of targeted buffer : FF

    auto ta2_zz_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 570);

    auto ta2_zz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 571);

    auto ta2_zz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 572);

    auto ta2_zz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 573);

    auto ta2_zz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 574);

    auto ta2_zz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 575);

    auto ta2_zz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 576);

    auto ta2_zz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 577);

    auto ta2_zz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 578);

    auto ta2_zz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 579);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yy_xxx_1, ta1_z_yy_xxy_1, ta1_z_yy_xyy_1, ta1_z_yy_xyz_1, ta1_z_yy_yyy_1, ta1_z_yy_yyz_1, ta1_z_yy_yzz_1, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yyz_xxx_0, ta2_zz_yyz_xxy_0, ta2_zz_yyz_xxz_0, ta2_zz_yyz_xyy_0, ta2_zz_yyz_xyz_0, ta2_zz_yyz_xzz_0, ta2_zz_yyz_yyy_0, ta2_zz_yyz_yyz_0, ta2_zz_yyz_yzz_0, ta2_zz_yyz_zzz_0, ta2_zz_yz_xxz_0, ta2_zz_yz_xxz_1, ta2_zz_yz_xzz_0, ta2_zz_yz_xzz_1, ta2_zz_yz_zzz_0, ta2_zz_yz_zzz_1, ta2_zz_z_xxz_0, ta2_zz_z_xxz_1, ta2_zz_z_xzz_0, ta2_zz_z_xzz_1, ta2_zz_z_zzz_0, ta2_zz_z_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyz_xxx_0[i] = 2.0 * ta1_z_yy_xxx_1[i] + ta2_zz_yy_xxx_0[i] * pa_z[i] - ta2_zz_yy_xxx_1[i] * pc_z[i];

        ta2_zz_yyz_xxy_0[i] = 2.0 * ta1_z_yy_xxy_1[i] + ta2_zz_yy_xxy_0[i] * pa_z[i] - ta2_zz_yy_xxy_1[i] * pc_z[i];

        ta2_zz_yyz_xxz_0[i] = ta2_zz_z_xxz_0[i] * fe_0 - ta2_zz_z_xxz_1[i] * fe_0 + ta2_zz_yz_xxz_0[i] * pa_y[i] - ta2_zz_yz_xxz_1[i] * pc_y[i];

        ta2_zz_yyz_xyy_0[i] = 2.0 * ta1_z_yy_xyy_1[i] + ta2_zz_yy_xyy_0[i] * pa_z[i] - ta2_zz_yy_xyy_1[i] * pc_z[i];

        ta2_zz_yyz_xyz_0[i] = ta2_zz_yy_xy_0[i] * fe_0 - ta2_zz_yy_xy_1[i] * fe_0 + 2.0 * ta1_z_yy_xyz_1[i] + ta2_zz_yy_xyz_0[i] * pa_z[i] - ta2_zz_yy_xyz_1[i] * pc_z[i];

        ta2_zz_yyz_xzz_0[i] = ta2_zz_z_xzz_0[i] * fe_0 - ta2_zz_z_xzz_1[i] * fe_0 + ta2_zz_yz_xzz_0[i] * pa_y[i] - ta2_zz_yz_xzz_1[i] * pc_y[i];

        ta2_zz_yyz_yyy_0[i] = 2.0 * ta1_z_yy_yyy_1[i] + ta2_zz_yy_yyy_0[i] * pa_z[i] - ta2_zz_yy_yyy_1[i] * pc_z[i];

        ta2_zz_yyz_yyz_0[i] = ta2_zz_yy_yy_0[i] * fe_0 - ta2_zz_yy_yy_1[i] * fe_0 + 2.0 * ta1_z_yy_yyz_1[i] + ta2_zz_yy_yyz_0[i] * pa_z[i] - ta2_zz_yy_yyz_1[i] * pc_z[i];

        ta2_zz_yyz_yzz_0[i] = 2.0 * ta2_zz_yy_yz_0[i] * fe_0 - 2.0 * ta2_zz_yy_yz_1[i] * fe_0 + 2.0 * ta1_z_yy_yzz_1[i] + ta2_zz_yy_yzz_0[i] * pa_z[i] - ta2_zz_yy_yzz_1[i] * pc_z[i];

        ta2_zz_yyz_zzz_0[i] = ta2_zz_z_zzz_0[i] * fe_0 - ta2_zz_z_zzz_1[i] * fe_0 + ta2_zz_yz_zzz_0[i] * pa_y[i] - ta2_zz_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 580-590 components of targeted buffer : FF

    auto ta2_zz_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 580);

    auto ta2_zz_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 581);

    auto ta2_zz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 582);

    auto ta2_zz_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 583);

    auto ta2_zz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 584);

    auto ta2_zz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 585);

    auto ta2_zz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 586);

    auto ta2_zz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 587);

    auto ta2_zz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 588);

    auto ta2_zz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 589);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_yzz_xxx_0, ta2_zz_yzz_xxy_0, ta2_zz_yzz_xxz_0, ta2_zz_yzz_xyy_0, ta2_zz_yzz_xyz_0, ta2_zz_yzz_xzz_0, ta2_zz_yzz_yyy_0, ta2_zz_yzz_yyz_0, ta2_zz_yzz_yzz_0, ta2_zz_yzz_zzz_0, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzz_xxx_0[i] = ta2_zz_zz_xxx_0[i] * pa_y[i] - ta2_zz_zz_xxx_1[i] * pc_y[i];

        ta2_zz_yzz_xxy_0[i] = ta2_zz_zz_xx_0[i] * fe_0 - ta2_zz_zz_xx_1[i] * fe_0 + ta2_zz_zz_xxy_0[i] * pa_y[i] - ta2_zz_zz_xxy_1[i] * pc_y[i];

        ta2_zz_yzz_xxz_0[i] = ta2_zz_zz_xxz_0[i] * pa_y[i] - ta2_zz_zz_xxz_1[i] * pc_y[i];

        ta2_zz_yzz_xyy_0[i] = 2.0 * ta2_zz_zz_xy_0[i] * fe_0 - 2.0 * ta2_zz_zz_xy_1[i] * fe_0 + ta2_zz_zz_xyy_0[i] * pa_y[i] - ta2_zz_zz_xyy_1[i] * pc_y[i];

        ta2_zz_yzz_xyz_0[i] = ta2_zz_zz_xz_0[i] * fe_0 - ta2_zz_zz_xz_1[i] * fe_0 + ta2_zz_zz_xyz_0[i] * pa_y[i] - ta2_zz_zz_xyz_1[i] * pc_y[i];

        ta2_zz_yzz_xzz_0[i] = ta2_zz_zz_xzz_0[i] * pa_y[i] - ta2_zz_zz_xzz_1[i] * pc_y[i];

        ta2_zz_yzz_yyy_0[i] = 3.0 * ta2_zz_zz_yy_0[i] * fe_0 - 3.0 * ta2_zz_zz_yy_1[i] * fe_0 + ta2_zz_zz_yyy_0[i] * pa_y[i] - ta2_zz_zz_yyy_1[i] * pc_y[i];

        ta2_zz_yzz_yyz_0[i] = 2.0 * ta2_zz_zz_yz_0[i] * fe_0 - 2.0 * ta2_zz_zz_yz_1[i] * fe_0 + ta2_zz_zz_yyz_0[i] * pa_y[i] - ta2_zz_zz_yyz_1[i] * pc_y[i];

        ta2_zz_yzz_yzz_0[i] = ta2_zz_zz_zz_0[i] * fe_0 - ta2_zz_zz_zz_1[i] * fe_0 + ta2_zz_zz_yzz_0[i] * pa_y[i] - ta2_zz_zz_yzz_1[i] * pc_y[i];

        ta2_zz_yzz_zzz_0[i] = ta2_zz_zz_zzz_0[i] * pa_y[i] - ta2_zz_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 590-600 components of targeted buffer : FF

    auto ta2_zz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 590);

    auto ta2_zz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 591);

    auto ta2_zz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 592);

    auto ta2_zz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 593);

    auto ta2_zz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 594);

    auto ta2_zz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 595);

    auto ta2_zz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 596);

    auto ta2_zz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 597);

    auto ta2_zz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 598);

    auto ta2_zz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 599);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zz_xxx_1, ta1_z_zz_xxy_1, ta1_z_zz_xxz_1, ta1_z_zz_xyy_1, ta1_z_zz_xyz_1, ta1_z_zz_xzz_1, ta1_z_zz_yyy_1, ta1_z_zz_yyz_1, ta1_z_zz_yzz_1, ta1_z_zz_zzz_1, ta2_zz_z_xxx_0, ta2_zz_z_xxx_1, ta2_zz_z_xxy_0, ta2_zz_z_xxy_1, ta2_zz_z_xxz_0, ta2_zz_z_xxz_1, ta2_zz_z_xyy_0, ta2_zz_z_xyy_1, ta2_zz_z_xyz_0, ta2_zz_z_xyz_1, ta2_zz_z_xzz_0, ta2_zz_z_xzz_1, ta2_zz_z_yyy_0, ta2_zz_z_yyy_1, ta2_zz_z_yyz_0, ta2_zz_z_yyz_1, ta2_zz_z_yzz_0, ta2_zz_z_yzz_1, ta2_zz_z_zzz_0, ta2_zz_z_zzz_1, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, ta2_zz_zzz_xxx_0, ta2_zz_zzz_xxy_0, ta2_zz_zzz_xxz_0, ta2_zz_zzz_xyy_0, ta2_zz_zzz_xyz_0, ta2_zz_zzz_xzz_0, ta2_zz_zzz_yyy_0, ta2_zz_zzz_yyz_0, ta2_zz_zzz_yzz_0, ta2_zz_zzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzz_xxx_0[i] = 2.0 * ta2_zz_z_xxx_0[i] * fe_0 - 2.0 * ta2_zz_z_xxx_1[i] * fe_0 + 2.0 * ta1_z_zz_xxx_1[i] + ta2_zz_zz_xxx_0[i] * pa_z[i] - ta2_zz_zz_xxx_1[i] * pc_z[i];

        ta2_zz_zzz_xxy_0[i] = 2.0 * ta2_zz_z_xxy_0[i] * fe_0 - 2.0 * ta2_zz_z_xxy_1[i] * fe_0 + 2.0 * ta1_z_zz_xxy_1[i] + ta2_zz_zz_xxy_0[i] * pa_z[i] - ta2_zz_zz_xxy_1[i] * pc_z[i];

        ta2_zz_zzz_xxz_0[i] = 2.0 * ta2_zz_z_xxz_0[i] * fe_0 - 2.0 * ta2_zz_z_xxz_1[i] * fe_0 + ta2_zz_zz_xx_0[i] * fe_0 - ta2_zz_zz_xx_1[i] * fe_0 + 2.0 * ta1_z_zz_xxz_1[i] + ta2_zz_zz_xxz_0[i] * pa_z[i] - ta2_zz_zz_xxz_1[i] * pc_z[i];

        ta2_zz_zzz_xyy_0[i] = 2.0 * ta2_zz_z_xyy_0[i] * fe_0 - 2.0 * ta2_zz_z_xyy_1[i] * fe_0 + 2.0 * ta1_z_zz_xyy_1[i] + ta2_zz_zz_xyy_0[i] * pa_z[i] - ta2_zz_zz_xyy_1[i] * pc_z[i];

        ta2_zz_zzz_xyz_0[i] = 2.0 * ta2_zz_z_xyz_0[i] * fe_0 - 2.0 * ta2_zz_z_xyz_1[i] * fe_0 + ta2_zz_zz_xy_0[i] * fe_0 - ta2_zz_zz_xy_1[i] * fe_0 + 2.0 * ta1_z_zz_xyz_1[i] + ta2_zz_zz_xyz_0[i] * pa_z[i] - ta2_zz_zz_xyz_1[i] * pc_z[i];

        ta2_zz_zzz_xzz_0[i] = 2.0 * ta2_zz_z_xzz_0[i] * fe_0 - 2.0 * ta2_zz_z_xzz_1[i] * fe_0 + 2.0 * ta2_zz_zz_xz_0[i] * fe_0 - 2.0 * ta2_zz_zz_xz_1[i] * fe_0 + 2.0 * ta1_z_zz_xzz_1[i] + ta2_zz_zz_xzz_0[i] * pa_z[i] - ta2_zz_zz_xzz_1[i] * pc_z[i];

        ta2_zz_zzz_yyy_0[i] = 2.0 * ta2_zz_z_yyy_0[i] * fe_0 - 2.0 * ta2_zz_z_yyy_1[i] * fe_0 + 2.0 * ta1_z_zz_yyy_1[i] + ta2_zz_zz_yyy_0[i] * pa_z[i] - ta2_zz_zz_yyy_1[i] * pc_z[i];

        ta2_zz_zzz_yyz_0[i] = 2.0 * ta2_zz_z_yyz_0[i] * fe_0 - 2.0 * ta2_zz_z_yyz_1[i] * fe_0 + ta2_zz_zz_yy_0[i] * fe_0 - ta2_zz_zz_yy_1[i] * fe_0 + 2.0 * ta1_z_zz_yyz_1[i] + ta2_zz_zz_yyz_0[i] * pa_z[i] - ta2_zz_zz_yyz_1[i] * pc_z[i];

        ta2_zz_zzz_yzz_0[i] = 2.0 * ta2_zz_z_yzz_0[i] * fe_0 - 2.0 * ta2_zz_z_yzz_1[i] * fe_0 + 2.0 * ta2_zz_zz_yz_0[i] * fe_0 - 2.0 * ta2_zz_zz_yz_1[i] * fe_0 + 2.0 * ta1_z_zz_yzz_1[i] + ta2_zz_zz_yzz_0[i] * pa_z[i] - ta2_zz_zz_yzz_1[i] * pc_z[i];

        ta2_zz_zzz_zzz_0[i] = 2.0 * ta2_zz_z_zzz_0[i] * fe_0 - 2.0 * ta2_zz_z_zzz_1[i] * fe_0 + 3.0 * ta2_zz_zz_zz_0[i] * fe_0 - 3.0 * ta2_zz_zz_zz_1[i] * fe_0 + 2.0 * ta1_z_zz_zzz_1[i] + ta2_zz_zz_zzz_0[i] * pa_z[i] - ta2_zz_zz_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

