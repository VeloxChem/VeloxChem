#include "NuclearPotentialGeom020PrimRecGF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_gf(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_gf,
                                        const size_t idx_npot_geom_020_0_df,
                                        const size_t idx_npot_geom_020_1_df,
                                        const size_t idx_npot_geom_020_0_fd,
                                        const size_t idx_npot_geom_020_1_fd,
                                        const size_t idx_npot_geom_010_1_ff,
                                        const size_t idx_npot_geom_020_0_ff,
                                        const size_t idx_npot_geom_020_1_ff,
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

    auto ta2_xx_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 12);

    auto ta2_xx_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 15);

    auto ta2_xx_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 20);

    auto ta2_xx_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 21);

    auto ta2_xx_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 22);

    auto ta2_xx_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 23);

    auto ta2_xx_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 25);

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

    auto ta2_xy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 71);

    auto ta2_xy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 73);

    auto ta2_xy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 76);

    auto ta2_xy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 77);

    auto ta2_xy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 78);

    auto ta2_xy_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 80);

    auto ta2_xy_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 81);

    auto ta2_xy_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 83);

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

    auto ta2_xz_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 132);

    auto ta2_xz_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 135);

    auto ta2_xz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 142);

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

    auto ta2_yy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 196);

    auto ta2_yy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 197);

    auto ta2_yy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 198);

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

    auto ta2_yy_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 226);

    auto ta2_yy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 227);

    auto ta2_yy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 228);

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

    auto ta2_yz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 256);

    auto ta2_yz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 257);

    auto ta2_yz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 258);

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

    auto ta2_yz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 285);

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

    auto ta2_zz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 316);

    auto ta2_zz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 317);

    auto ta2_zz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 318);

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

    auto ta2_zz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 345);

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

    auto ta2_xx_xy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 12);

    auto ta2_xx_xy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 15);

    auto ta2_xx_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 20);

    auto ta2_xx_xz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 21);

    auto ta2_xx_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 22);

    auto ta2_xx_xz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 23);

    auto ta2_xx_xz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 25);

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

    auto ta2_xy_xy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 71);

    auto ta2_xy_xy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 73);

    auto ta2_xy_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 76);

    auto ta2_xy_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 77);

    auto ta2_xy_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 78);

    auto ta2_xy_xz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_df + 80);

    auto ta2_xy_xz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_df + 81);

    auto ta2_xy_xz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 83);

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

    auto ta2_xz_xy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 132);

    auto ta2_xz_xy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 135);

    auto ta2_xz_xz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_df + 142);

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

    auto ta2_yy_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 196);

    auto ta2_yy_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 197);

    auto ta2_yy_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 198);

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

    auto ta2_yy_yz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 226);

    auto ta2_yy_yz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 227);

    auto ta2_yy_yz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 228);

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

    auto ta2_yz_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 256);

    auto ta2_yz_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 257);

    auto ta2_yz_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 258);

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

    auto ta2_yz_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 285);

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

    auto ta2_zz_xy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_df + 316);

    auto ta2_zz_xy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_df + 317);

    auto ta2_zz_xy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 318);

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

    auto ta2_zz_yz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_df + 345);

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

    // Set up components of auxiliary buffer : FD

    auto ta2_xx_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd);

    auto ta2_xx_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 1);

    auto ta2_xx_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 2);

    auto ta2_xx_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 3);

    auto ta2_xx_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 4);

    auto ta2_xx_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 5);

    auto ta2_xx_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 6);

    auto ta2_xx_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 7);

    auto ta2_xx_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 8);

    auto ta2_xx_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 12);

    auto ta2_xx_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 13);

    auto ta2_xx_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 14);

    auto ta2_xx_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 16);

    auto ta2_xx_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 17);

    auto ta2_xx_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 19);

    auto ta2_xx_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 30);

    auto ta2_xx_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 31);

    auto ta2_xx_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 32);

    auto ta2_xx_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 36);

    auto ta2_xx_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 37);

    auto ta2_xx_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 38);

    auto ta2_xx_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 39);

    auto ta2_xx_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 40);

    auto ta2_xx_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 41);

    auto ta2_xx_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 50);

    auto ta2_xx_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 52);

    auto ta2_xx_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 53);

    auto ta2_xx_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 54);

    auto ta2_xx_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 55);

    auto ta2_xx_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 56);

    auto ta2_xx_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 57);

    auto ta2_xx_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 58);

    auto ta2_xx_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 59);

    auto ta2_xy_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 60);

    auto ta2_xy_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 61);

    auto ta2_xy_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 62);

    auto ta2_xy_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 63);

    auto ta2_xy_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 64);

    auto ta2_xy_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 65);

    auto ta2_xy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 67);

    auto ta2_xy_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 69);

    auto ta2_xy_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 70);

    auto ta2_xy_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 72);

    auto ta2_xy_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 73);

    auto ta2_xy_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 74);

    auto ta2_xy_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 78);

    auto ta2_xy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 79);

    auto ta2_xy_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 80);

    auto ta2_xy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 81);

    auto ta2_xy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 82);

    auto ta2_xy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 96);

    auto ta2_xy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 97);

    auto ta2_xy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 98);

    auto ta2_xy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 99);

    auto ta2_xy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 100);

    auto ta2_xy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 101);

    auto ta2_xy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 103);

    auto ta2_xy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 105);

    auto ta2_xy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 106);

    auto ta2_xy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 112);

    auto ta2_xy_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 114);

    auto ta2_xy_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 115);

    auto ta2_xy_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 116);

    auto ta2_xy_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 117);

    auto ta2_xy_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 118);

    auto ta2_xy_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 119);

    auto ta2_xz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 120);

    auto ta2_xz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 121);

    auto ta2_xz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 122);

    auto ta2_xz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 123);

    auto ta2_xz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 124);

    auto ta2_xz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 125);

    auto ta2_xz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 126);

    auto ta2_xz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 127);

    auto ta2_xz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 128);

    auto ta2_xz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 134);

    auto ta2_xz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 136);

    auto ta2_xz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 137);

    auto ta2_xz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 150);

    auto ta2_xz_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 151);

    auto ta2_xz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 152);

    auto ta2_xz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 154);

    auto ta2_xz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 155);

    auto ta2_xz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 156);

    auto ta2_xz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 157);

    auto ta2_xz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 158);

    auto ta2_xz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 159);

    auto ta2_xz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 160);

    auto ta2_xz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 161);

    auto ta2_xz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 166);

    auto ta2_xz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 170);

    auto ta2_xz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 172);

    auto ta2_xz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 173);

    auto ta2_xz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 174);

    auto ta2_xz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 175);

    auto ta2_xz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 176);

    auto ta2_xz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 177);

    auto ta2_xz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 178);

    auto ta2_xz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 179);

    auto ta2_yy_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 180);

    auto ta2_yy_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 181);

    auto ta2_yy_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 182);

    auto ta2_yy_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 183);

    auto ta2_yy_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 184);

    auto ta2_yy_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 185);

    auto ta2_yy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 187);

    auto ta2_yy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 199);

    auto ta2_yy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 201);

    auto ta2_yy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 202);

    auto ta2_yy_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 212);

    auto ta2_yy_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 214);

    auto ta2_yy_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 215);

    auto ta2_yy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 216);

    auto ta2_yy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 217);

    auto ta2_yy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 218);

    auto ta2_yy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 219);

    auto ta2_yy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 220);

    auto ta2_yy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 221);

    auto ta2_yy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 223);

    auto ta2_yy_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 224);

    auto ta2_yy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 225);

    auto ta2_yy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 226);

    auto ta2_yy_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 227);

    auto ta2_yy_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 229);

    auto ta2_yy_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 231);

    auto ta2_yy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 232);

    auto ta2_yy_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 234);

    auto ta2_yy_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 235);

    auto ta2_yy_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 236);

    auto ta2_yy_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 237);

    auto ta2_yy_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 238);

    auto ta2_yy_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 239);

    auto ta2_yz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 240);

    auto ta2_yz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 241);

    auto ta2_yz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 242);

    auto ta2_yz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 243);

    auto ta2_yz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 244);

    auto ta2_yz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 245);

    auto ta2_yz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 254);

    auto ta2_yz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 259);

    auto ta2_yz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 261);

    auto ta2_yz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 262);

    auto ta2_yz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 272);

    auto ta2_yz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 274);

    auto ta2_yz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 275);

    auto ta2_yz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 276);

    auto ta2_yz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 277);

    auto ta2_yz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 278);

    auto ta2_yz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 279);

    auto ta2_yz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 280);

    auto ta2_yz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 281);

    auto ta2_yz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 284);

    auto ta2_yz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 286);

    auto ta2_yz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 287);

    auto ta2_yz_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 289);

    auto ta2_yz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 290);

    auto ta2_yz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 291);

    auto ta2_yz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 292);

    auto ta2_yz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 293);

    auto ta2_yz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 294);

    auto ta2_yz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 295);

    auto ta2_yz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 296);

    auto ta2_yz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 297);

    auto ta2_yz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 298);

    auto ta2_yz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 299);

    auto ta2_zz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 300);

    auto ta2_zz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 301);

    auto ta2_zz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 302);

    auto ta2_zz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 303);

    auto ta2_zz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 304);

    auto ta2_zz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 305);

    auto ta2_zz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 314);

    auto ta2_zz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 319);

    auto ta2_zz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 321);

    auto ta2_zz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 322);

    auto ta2_zz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 332);

    auto ta2_zz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 334);

    auto ta2_zz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 335);

    auto ta2_zz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 336);

    auto ta2_zz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 337);

    auto ta2_zz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 338);

    auto ta2_zz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 339);

    auto ta2_zz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 340);

    auto ta2_zz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 341);

    auto ta2_zz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 344);

    auto ta2_zz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 346);

    auto ta2_zz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 347);

    auto ta2_zz_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 349);

    auto ta2_zz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 350);

    auto ta2_zz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 351);

    auto ta2_zz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 352);

    auto ta2_zz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 353);

    auto ta2_zz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 354);

    auto ta2_zz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 355);

    auto ta2_zz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 356);

    auto ta2_zz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 357);

    auto ta2_zz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 358);

    auto ta2_zz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 359);

    // Set up components of auxiliary buffer : FD

    auto ta2_xx_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd);

    auto ta2_xx_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 1);

    auto ta2_xx_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 2);

    auto ta2_xx_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 3);

    auto ta2_xx_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 4);

    auto ta2_xx_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 5);

    auto ta2_xx_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 6);

    auto ta2_xx_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 7);

    auto ta2_xx_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 8);

    auto ta2_xx_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 12);

    auto ta2_xx_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 13);

    auto ta2_xx_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 14);

    auto ta2_xx_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 16);

    auto ta2_xx_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 17);

    auto ta2_xx_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 19);

    auto ta2_xx_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 30);

    auto ta2_xx_xzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 31);

    auto ta2_xx_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 32);

    auto ta2_xx_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 36);

    auto ta2_xx_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 37);

    auto ta2_xx_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 38);

    auto ta2_xx_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 39);

    auto ta2_xx_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 40);

    auto ta2_xx_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 41);

    auto ta2_xx_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 50);

    auto ta2_xx_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 52);

    auto ta2_xx_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 53);

    auto ta2_xx_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 54);

    auto ta2_xx_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 55);

    auto ta2_xx_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 56);

    auto ta2_xx_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 57);

    auto ta2_xx_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 58);

    auto ta2_xx_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 59);

    auto ta2_xy_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 60);

    auto ta2_xy_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 61);

    auto ta2_xy_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 62);

    auto ta2_xy_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 63);

    auto ta2_xy_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 64);

    auto ta2_xy_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 65);

    auto ta2_xy_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 67);

    auto ta2_xy_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 69);

    auto ta2_xy_xxy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 70);

    auto ta2_xy_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 72);

    auto ta2_xy_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 73);

    auto ta2_xy_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 74);

    auto ta2_xy_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 78);

    auto ta2_xy_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 79);

    auto ta2_xy_xyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 80);

    auto ta2_xy_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 81);

    auto ta2_xy_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 82);

    auto ta2_xy_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 96);

    auto ta2_xy_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 97);

    auto ta2_xy_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 98);

    auto ta2_xy_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 99);

    auto ta2_xy_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 100);

    auto ta2_xy_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 101);

    auto ta2_xy_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 103);

    auto ta2_xy_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 105);

    auto ta2_xy_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 106);

    auto ta2_xy_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 112);

    auto ta2_xy_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 114);

    auto ta2_xy_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 115);

    auto ta2_xy_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 116);

    auto ta2_xy_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 117);

    auto ta2_xy_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 118);

    auto ta2_xy_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 119);

    auto ta2_xz_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 120);

    auto ta2_xz_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 121);

    auto ta2_xz_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 122);

    auto ta2_xz_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 123);

    auto ta2_xz_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 124);

    auto ta2_xz_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 125);

    auto ta2_xz_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 126);

    auto ta2_xz_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 127);

    auto ta2_xz_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 128);

    auto ta2_xz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 134);

    auto ta2_xz_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 136);

    auto ta2_xz_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 137);

    auto ta2_xz_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 150);

    auto ta2_xz_xzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 151);

    auto ta2_xz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 152);

    auto ta2_xz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 154);

    auto ta2_xz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 155);

    auto ta2_xz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 156);

    auto ta2_xz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 157);

    auto ta2_xz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 158);

    auto ta2_xz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 159);

    auto ta2_xz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 160);

    auto ta2_xz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 161);

    auto ta2_xz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 166);

    auto ta2_xz_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 170);

    auto ta2_xz_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 172);

    auto ta2_xz_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 173);

    auto ta2_xz_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 174);

    auto ta2_xz_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 175);

    auto ta2_xz_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 176);

    auto ta2_xz_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 177);

    auto ta2_xz_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 178);

    auto ta2_xz_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 179);

    auto ta2_yy_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 180);

    auto ta2_yy_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 181);

    auto ta2_yy_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 182);

    auto ta2_yy_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 183);

    auto ta2_yy_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 184);

    auto ta2_yy_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 185);

    auto ta2_yy_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 187);

    auto ta2_yy_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 199);

    auto ta2_yy_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 201);

    auto ta2_yy_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 202);

    auto ta2_yy_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 212);

    auto ta2_yy_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 214);

    auto ta2_yy_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 215);

    auto ta2_yy_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 216);

    auto ta2_yy_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 217);

    auto ta2_yy_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 218);

    auto ta2_yy_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 219);

    auto ta2_yy_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 220);

    auto ta2_yy_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 221);

    auto ta2_yy_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 223);

    auto ta2_yy_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 224);

    auto ta2_yy_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 225);

    auto ta2_yy_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 226);

    auto ta2_yy_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 227);

    auto ta2_yy_yzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 229);

    auto ta2_yy_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 231);

    auto ta2_yy_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 232);

    auto ta2_yy_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 234);

    auto ta2_yy_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 235);

    auto ta2_yy_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 236);

    auto ta2_yy_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 237);

    auto ta2_yy_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 238);

    auto ta2_yy_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 239);

    auto ta2_yz_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 240);

    auto ta2_yz_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 241);

    auto ta2_yz_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 242);

    auto ta2_yz_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 243);

    auto ta2_yz_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 244);

    auto ta2_yz_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 245);

    auto ta2_yz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 254);

    auto ta2_yz_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 259);

    auto ta2_yz_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 261);

    auto ta2_yz_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 262);

    auto ta2_yz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 272);

    auto ta2_yz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 274);

    auto ta2_yz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 275);

    auto ta2_yz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 276);

    auto ta2_yz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 277);

    auto ta2_yz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 278);

    auto ta2_yz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 279);

    auto ta2_yz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 280);

    auto ta2_yz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 281);

    auto ta2_yz_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 284);

    auto ta2_yz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 286);

    auto ta2_yz_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 287);

    auto ta2_yz_yzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 289);

    auto ta2_yz_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 290);

    auto ta2_yz_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 291);

    auto ta2_yz_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 292);

    auto ta2_yz_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 293);

    auto ta2_yz_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 294);

    auto ta2_yz_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 295);

    auto ta2_yz_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 296);

    auto ta2_yz_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 297);

    auto ta2_yz_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 298);

    auto ta2_yz_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 299);

    auto ta2_zz_xxx_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 300);

    auto ta2_zz_xxx_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 301);

    auto ta2_zz_xxx_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 302);

    auto ta2_zz_xxx_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 303);

    auto ta2_zz_xxx_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 304);

    auto ta2_zz_xxx_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 305);

    auto ta2_zz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 314);

    auto ta2_zz_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 319);

    auto ta2_zz_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 321);

    auto ta2_zz_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 322);

    auto ta2_zz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 332);

    auto ta2_zz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 334);

    auto ta2_zz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 335);

    auto ta2_zz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 336);

    auto ta2_zz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 337);

    auto ta2_zz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 338);

    auto ta2_zz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 339);

    auto ta2_zz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 340);

    auto ta2_zz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 341);

    auto ta2_zz_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 344);

    auto ta2_zz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 346);

    auto ta2_zz_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 347);

    auto ta2_zz_yzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 349);

    auto ta2_zz_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 350);

    auto ta2_zz_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 351);

    auto ta2_zz_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 352);

    auto ta2_zz_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 353);

    auto ta2_zz_zzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 354);

    auto ta2_zz_zzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 355);

    auto ta2_zz_zzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 356);

    auto ta2_zz_zzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 357);

    auto ta2_zz_zzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 358);

    auto ta2_zz_zzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 359);

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

    auto ta1_x_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 15);

    auto ta1_x_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 16);

    auto ta1_x_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 20);

    auto ta1_x_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 21);

    auto ta1_x_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 22);

    auto ta1_x_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 23);

    auto ta1_x_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 25);

    auto ta1_x_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 29);

    auto ta1_x_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 30);

    auto ta1_x_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 31);

    auto ta1_x_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 33);

    auto ta1_x_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 36);

    auto ta1_x_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 37);

    auto ta1_x_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 38);

    auto ta1_x_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 50);

    auto ta1_x_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 52);

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

    auto ta1_x_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 73);

    auto ta1_x_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 76);

    auto ta1_x_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 77);

    auto ta1_x_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 78);

    auto ta1_x_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 79);

    auto ta1_x_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 82);

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

    auto ta1_y_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 115);

    auto ta1_y_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 116);

    auto ta1_y_xxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 117);

    auto ta1_y_xxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 118);

    auto ta1_y_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 120);

    auto ta1_y_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 121);

    auto ta1_y_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 122);

    auto ta1_y_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 123);

    auto ta1_y_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 125);

    auto ta1_y_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 129);

    auto ta1_y_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 130);

    auto ta1_y_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 131);

    auto ta1_y_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 133);

    auto ta1_y_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 134);

    auto ta1_y_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 136);

    auto ta1_y_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 137);

    auto ta1_y_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 138);

    auto ta1_y_xyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 139);

    auto ta1_y_xzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 150);

    auto ta1_y_xzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 152);

    auto ta1_y_xzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 155);

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

    auto ta1_y_yyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 171);

    auto ta1_y_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 173);

    auto ta1_y_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 176);

    auto ta1_y_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 177);

    auto ta1_y_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 178);

    auto ta1_y_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 179);

    auto ta1_y_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 182);

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

    auto ta1_z_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 236);

    auto ta1_z_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 237);

    auto ta1_z_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 238);

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

    auto ta1_z_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 282);

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

    // Set up components of auxiliary buffer : FF

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

    auto ta2_xx_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 10);

    auto ta2_xx_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 11);

    auto ta2_xx_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 12);

    auto ta2_xx_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 13);

    auto ta2_xx_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 14);

    auto ta2_xx_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 15);

    auto ta2_xx_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 16);

    auto ta2_xx_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 19);

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

    auto ta2_xx_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 30);

    auto ta2_xx_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 31);

    auto ta2_xx_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 32);

    auto ta2_xx_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 33);

    auto ta2_xx_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 34);

    auto ta2_xx_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 35);

    auto ta2_xx_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 36);

    auto ta2_xx_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 37);

    auto ta2_xx_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 38);

    auto ta2_xx_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 42);

    auto ta2_xx_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 45);

    auto ta2_xx_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 50);

    auto ta2_xx_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 51);

    auto ta2_xx_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 52);

    auto ta2_xx_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 53);

    auto ta2_xx_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 54);

    auto ta2_xx_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 55);

    auto ta2_xx_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 57);

    auto ta2_xx_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 58);

    auto ta2_xx_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 59);

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

    auto ta2_xx_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 71);

    auto ta2_xx_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 72);

    auto ta2_xx_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 73);

    auto ta2_xx_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 75);

    auto ta2_xx_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 76);

    auto ta2_xx_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 77);

    auto ta2_xx_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 78);

    auto ta2_xx_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 79);

    auto ta2_xx_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 80);

    auto ta2_xx_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 82);

    auto ta2_xx_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 84);

    auto ta2_xx_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 85);

    auto ta2_xx_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 86);

    auto ta2_xx_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 87);

    auto ta2_xx_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 88);

    auto ta2_xx_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 89);

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

    auto ta2_xy_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 110);

    auto ta2_xy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 111);

    auto ta2_xy_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 112);

    auto ta2_xy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 113);

    auto ta2_xy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 114);

    auto ta2_xy_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 115);

    auto ta2_xy_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 116);

    auto ta2_xy_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 117);

    auto ta2_xy_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 118);

    auto ta2_xy_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 120);

    auto ta2_xy_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 121);

    auto ta2_xy_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 122);

    auto ta2_xy_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 123);

    auto ta2_xy_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 124);

    auto ta2_xy_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 125);

    auto ta2_xy_xxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 126);

    auto ta2_xy_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 129);

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

    auto ta2_xy_xyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 141);

    auto ta2_xy_xyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 143);

    auto ta2_xy_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 150);

    auto ta2_xy_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 151);

    auto ta2_xy_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 152);

    auto ta2_xy_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 153);

    auto ta2_xy_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 155);

    auto ta2_xy_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 157);

    auto ta2_xy_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 158);

    auto ta2_xy_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 159);

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

    auto ta2_xy_yyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 170);

    auto ta2_xy_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 171);

    auto ta2_xy_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 173);

    auto ta2_xy_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 174);

    auto ta2_xy_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 176);

    auto ta2_xy_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 177);

    auto ta2_xy_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 178);

    auto ta2_xy_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 179);

    auto ta2_xy_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 181);

    auto ta2_xy_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 182);

    auto ta2_xy_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 183);

    auto ta2_xy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 184);

    auto ta2_xy_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 185);

    auto ta2_xy_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 186);

    auto ta2_xy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 187);

    auto ta2_xy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 188);

    auto ta2_xy_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 189);

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

    auto ta2_xz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 210);

    auto ta2_xz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 211);

    auto ta2_xz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 212);

    auto ta2_xz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 213);

    auto ta2_xz_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 214);

    auto ta2_xz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 215);

    auto ta2_xz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 216);

    auto ta2_xz_xxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 219);

    auto ta2_xz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 220);

    auto ta2_xz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 221);

    auto ta2_xz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 222);

    auto ta2_xz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 223);

    auto ta2_xz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 224);

    auto ta2_xz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 225);

    auto ta2_xz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 227);

    auto ta2_xz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 228);

    auto ta2_xz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 229);

    auto ta2_xz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 230);

    auto ta2_xz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 231);

    auto ta2_xz_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 232);

    auto ta2_xz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 233);

    auto ta2_xz_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 235);

    auto ta2_xz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 236);

    auto ta2_xz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 237);

    auto ta2_xz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 238);

    auto ta2_xz_xyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 242);

    auto ta2_xz_xyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 245);

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

    auto ta2_xz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 271);

    auto ta2_xz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 272);

    auto ta2_xz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 273);

    auto ta2_xz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 274);

    auto ta2_xz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 275);

    auto ta2_xz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 276);

    auto ta2_xz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 277);

    auto ta2_xz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 278);

    auto ta2_xz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 279);

    auto ta2_xz_yzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 280);

    auto ta2_xz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 282);

    auto ta2_xz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 284);

    auto ta2_xz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 285);

    auto ta2_xz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 286);

    auto ta2_xz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 287);

    auto ta2_xz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 288);

    auto ta2_xz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 289);

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

    auto ta2_yy_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 310);

    auto ta2_yy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 311);

    auto ta2_yy_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 312);

    auto ta2_yy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 313);

    auto ta2_yy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 314);

    auto ta2_yy_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 315);

    auto ta2_yy_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 316);

    auto ta2_yy_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 317);

    auto ta2_yy_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 318);

    auto ta2_yy_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 320);

    auto ta2_yy_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 321);

    auto ta2_yy_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 322);

    auto ta2_yy_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 323);

    auto ta2_yy_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 325);

    auto ta2_yy_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 327);

    auto ta2_yy_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 328);

    auto ta2_yy_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 329);

    auto ta2_yy_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 330);

    auto ta2_yy_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 331);

    auto ta2_yy_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 333);

    auto ta2_yy_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 334);

    auto ta2_yy_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 336);

    auto ta2_yy_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 337);

    auto ta2_yy_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 338);

    auto ta2_yy_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 339);

    auto ta2_yy_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 347);

    auto ta2_yy_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 348);

    auto ta2_yy_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 350);

    auto ta2_yy_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 352);

    auto ta2_yy_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 354);

    auto ta2_yy_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 355);

    auto ta2_yy_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 356);

    auto ta2_yy_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 357);

    auto ta2_yy_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 358);

    auto ta2_yy_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 359);

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

    auto ta2_yy_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 381);

    auto ta2_yy_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 382);

    auto ta2_yy_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 383);

    auto ta2_yy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 384);

    auto ta2_yy_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 385);

    auto ta2_yy_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 386);

    auto ta2_yy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 387);

    auto ta2_yy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 388);

    auto ta2_yy_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 389);

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

    auto ta2_yz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 410);

    auto ta2_yz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 411);

    auto ta2_yz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 412);

    auto ta2_yz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 413);

    auto ta2_yz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 415);

    auto ta2_yz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 416);

    auto ta2_yz_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 417);

    auto ta2_yz_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 418);

    auto ta2_yz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 420);

    auto ta2_yz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 421);

    auto ta2_yz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 422);

    auto ta2_yz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 423);

    auto ta2_yz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 424);

    auto ta2_yz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 425);

    auto ta2_yz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 427);

    auto ta2_yz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 428);

    auto ta2_yz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 429);

    auto ta2_yz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 430);

    auto ta2_yz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 431);

    auto ta2_yz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 433);

    auto ta2_yz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 434);

    auto ta2_yz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 436);

    auto ta2_yz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 437);

    auto ta2_yz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 438);

    auto ta2_yz_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 439);

    auto ta2_yz_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 447);

    auto ta2_yz_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 448);

    auto ta2_yz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 450);

    auto ta2_yz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 452);

    auto ta2_yz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 454);

    auto ta2_yz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 455);

    auto ta2_yz_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 456);

    auto ta2_yz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 457);

    auto ta2_yz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 458);

    auto ta2_yz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 459);

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

    auto ta2_yz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 471);

    auto ta2_yz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 472);

    auto ta2_yz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 473);

    auto ta2_yz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 474);

    auto ta2_yz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 475);

    auto ta2_yz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 476);

    auto ta2_yz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 477);

    auto ta2_yz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 478);

    auto ta2_yz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 479);

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

    auto ta2_zz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 510);

    auto ta2_zz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 511);

    auto ta2_zz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 512);

    auto ta2_zz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 513);

    auto ta2_zz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 515);

    auto ta2_zz_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 516);

    auto ta2_zz_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 517);

    auto ta2_zz_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 518);

    auto ta2_zz_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 520);

    auto ta2_zz_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 521);

    auto ta2_zz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 522);

    auto ta2_zz_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 523);

    auto ta2_zz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 524);

    auto ta2_zz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 525);

    auto ta2_zz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 527);

    auto ta2_zz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 528);

    auto ta2_zz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 529);

    auto ta2_zz_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 530);

    auto ta2_zz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 531);

    auto ta2_zz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 533);

    auto ta2_zz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 534);

    auto ta2_zz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 536);

    auto ta2_zz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 537);

    auto ta2_zz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 538);

    auto ta2_zz_xyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 539);

    auto ta2_zz_xyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 547);

    auto ta2_zz_xyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 548);

    auto ta2_zz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 550);

    auto ta2_zz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 552);

    auto ta2_zz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 554);

    auto ta2_zz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 555);

    auto ta2_zz_xzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 556);

    auto ta2_zz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 557);

    auto ta2_zz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 558);

    auto ta2_zz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 559);

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

    auto ta2_zz_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 571);

    auto ta2_zz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 572);

    auto ta2_zz_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 573);

    auto ta2_zz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 574);

    auto ta2_zz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 575);

    auto ta2_zz_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 576);

    auto ta2_zz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 577);

    auto ta2_zz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 578);

    auto ta2_zz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 579);

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

    // Set up components of auxiliary buffer : FF

    auto ta2_xx_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff);

    auto ta2_xx_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 1);

    auto ta2_xx_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 2);

    auto ta2_xx_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 3);

    auto ta2_xx_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 4);

    auto ta2_xx_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 5);

    auto ta2_xx_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 6);

    auto ta2_xx_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 7);

    auto ta2_xx_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 8);

    auto ta2_xx_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 9);

    auto ta2_xx_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 10);

    auto ta2_xx_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 11);

    auto ta2_xx_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 12);

    auto ta2_xx_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 13);

    auto ta2_xx_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 14);

    auto ta2_xx_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 15);

    auto ta2_xx_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 16);

    auto ta2_xx_xxy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 19);

    auto ta2_xx_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 20);

    auto ta2_xx_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 21);

    auto ta2_xx_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 22);

    auto ta2_xx_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 23);

    auto ta2_xx_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 24);

    auto ta2_xx_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 25);

    auto ta2_xx_xxz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 26);

    auto ta2_xx_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 27);

    auto ta2_xx_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 28);

    auto ta2_xx_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 29);

    auto ta2_xx_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 30);

    auto ta2_xx_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 31);

    auto ta2_xx_xyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 32);

    auto ta2_xx_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 33);

    auto ta2_xx_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 34);

    auto ta2_xx_xyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 35);

    auto ta2_xx_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 36);

    auto ta2_xx_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 37);

    auto ta2_xx_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 38);

    auto ta2_xx_xyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 42);

    auto ta2_xx_xyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 45);

    auto ta2_xx_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 50);

    auto ta2_xx_xzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 51);

    auto ta2_xx_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 52);

    auto ta2_xx_xzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 53);

    auto ta2_xx_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 54);

    auto ta2_xx_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 55);

    auto ta2_xx_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 57);

    auto ta2_xx_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 58);

    auto ta2_xx_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 59);

    auto ta2_xx_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 60);

    auto ta2_xx_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 61);

    auto ta2_xx_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 62);

    auto ta2_xx_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 63);

    auto ta2_xx_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 64);

    auto ta2_xx_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 65);

    auto ta2_xx_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 66);

    auto ta2_xx_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 67);

    auto ta2_xx_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 68);

    auto ta2_xx_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 69);

    auto ta2_xx_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 71);

    auto ta2_xx_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 72);

    auto ta2_xx_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 73);

    auto ta2_xx_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 75);

    auto ta2_xx_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 76);

    auto ta2_xx_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 77);

    auto ta2_xx_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 78);

    auto ta2_xx_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 79);

    auto ta2_xx_yzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 80);

    auto ta2_xx_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 82);

    auto ta2_xx_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 84);

    auto ta2_xx_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 85);

    auto ta2_xx_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 86);

    auto ta2_xx_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 87);

    auto ta2_xx_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 88);

    auto ta2_xx_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 89);

    auto ta2_xx_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 90);

    auto ta2_xx_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 91);

    auto ta2_xx_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 92);

    auto ta2_xx_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 93);

    auto ta2_xx_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 94);

    auto ta2_xx_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 95);

    auto ta2_xx_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 96);

    auto ta2_xx_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 97);

    auto ta2_xx_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 98);

    auto ta2_xx_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 99);

    auto ta2_xy_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 100);

    auto ta2_xy_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 101);

    auto ta2_xy_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 102);

    auto ta2_xy_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 103);

    auto ta2_xy_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 104);

    auto ta2_xy_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 105);

    auto ta2_xy_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 106);

    auto ta2_xy_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 107);

    auto ta2_xy_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 108);

    auto ta2_xy_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 109);

    auto ta2_xy_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 110);

    auto ta2_xy_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 111);

    auto ta2_xy_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 112);

    auto ta2_xy_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 113);

    auto ta2_xy_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 114);

    auto ta2_xy_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 115);

    auto ta2_xy_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 116);

    auto ta2_xy_xxy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 117);

    auto ta2_xy_xxy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 118);

    auto ta2_xy_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 120);

    auto ta2_xy_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 121);

    auto ta2_xy_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 122);

    auto ta2_xy_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 123);

    auto ta2_xy_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 124);

    auto ta2_xy_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 125);

    auto ta2_xy_xxz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 126);

    auto ta2_xy_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 129);

    auto ta2_xy_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 130);

    auto ta2_xy_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 131);

    auto ta2_xy_xyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 132);

    auto ta2_xy_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 133);

    auto ta2_xy_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 134);

    auto ta2_xy_xyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 135);

    auto ta2_xy_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 136);

    auto ta2_xy_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 137);

    auto ta2_xy_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 138);

    auto ta2_xy_xyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 139);

    auto ta2_xy_xyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 141);

    auto ta2_xy_xyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 143);

    auto ta2_xy_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 150);

    auto ta2_xy_xzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 151);

    auto ta2_xy_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 152);

    auto ta2_xy_xzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 153);

    auto ta2_xy_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 155);

    auto ta2_xy_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 157);

    auto ta2_xy_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 158);

    auto ta2_xy_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 159);

    auto ta2_xy_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 160);

    auto ta2_xy_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 161);

    auto ta2_xy_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 162);

    auto ta2_xy_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 163);

    auto ta2_xy_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 164);

    auto ta2_xy_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 165);

    auto ta2_xy_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 166);

    auto ta2_xy_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 167);

    auto ta2_xy_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 168);

    auto ta2_xy_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 169);

    auto ta2_xy_yyz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 170);

    auto ta2_xy_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 171);

    auto ta2_xy_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 173);

    auto ta2_xy_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 174);

    auto ta2_xy_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 176);

    auto ta2_xy_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 177);

    auto ta2_xy_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 178);

    auto ta2_xy_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 179);

    auto ta2_xy_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 181);

    auto ta2_xy_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 182);

    auto ta2_xy_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 183);

    auto ta2_xy_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 184);

    auto ta2_xy_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 185);

    auto ta2_xy_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 186);

    auto ta2_xy_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 187);

    auto ta2_xy_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 188);

    auto ta2_xy_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 189);

    auto ta2_xy_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 190);

    auto ta2_xy_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 191);

    auto ta2_xy_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 192);

    auto ta2_xy_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 193);

    auto ta2_xy_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 194);

    auto ta2_xy_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 195);

    auto ta2_xy_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 196);

    auto ta2_xy_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 197);

    auto ta2_xy_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 198);

    auto ta2_xy_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 199);

    auto ta2_xz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 200);

    auto ta2_xz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 201);

    auto ta2_xz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 202);

    auto ta2_xz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 203);

    auto ta2_xz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 204);

    auto ta2_xz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 205);

    auto ta2_xz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 206);

    auto ta2_xz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 207);

    auto ta2_xz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 208);

    auto ta2_xz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 209);

    auto ta2_xz_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 210);

    auto ta2_xz_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 211);

    auto ta2_xz_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 212);

    auto ta2_xz_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 213);

    auto ta2_xz_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 214);

    auto ta2_xz_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 215);

    auto ta2_xz_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 216);

    auto ta2_xz_xxy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 219);

    auto ta2_xz_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 220);

    auto ta2_xz_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 221);

    auto ta2_xz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 222);

    auto ta2_xz_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 223);

    auto ta2_xz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 224);

    auto ta2_xz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 225);

    auto ta2_xz_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 227);

    auto ta2_xz_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 228);

    auto ta2_xz_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 229);

    auto ta2_xz_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 230);

    auto ta2_xz_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 231);

    auto ta2_xz_xyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 232);

    auto ta2_xz_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 233);

    auto ta2_xz_xyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 235);

    auto ta2_xz_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 236);

    auto ta2_xz_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 237);

    auto ta2_xz_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 238);

    auto ta2_xz_xyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 242);

    auto ta2_xz_xyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 245);

    auto ta2_xz_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 250);

    auto ta2_xz_xzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 251);

    auto ta2_xz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 252);

    auto ta2_xz_xzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 253);

    auto ta2_xz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 254);

    auto ta2_xz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 255);

    auto ta2_xz_xzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 256);

    auto ta2_xz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 257);

    auto ta2_xz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 258);

    auto ta2_xz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 259);

    auto ta2_xz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 260);

    auto ta2_xz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 261);

    auto ta2_xz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 262);

    auto ta2_xz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 263);

    auto ta2_xz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 264);

    auto ta2_xz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 265);

    auto ta2_xz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 266);

    auto ta2_xz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 267);

    auto ta2_xz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 268);

    auto ta2_xz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 269);

    auto ta2_xz_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 271);

    auto ta2_xz_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 272);

    auto ta2_xz_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 273);

    auto ta2_xz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 274);

    auto ta2_xz_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 275);

    auto ta2_xz_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 276);

    auto ta2_xz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 277);

    auto ta2_xz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 278);

    auto ta2_xz_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 279);

    auto ta2_xz_yzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 280);

    auto ta2_xz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 282);

    auto ta2_xz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 284);

    auto ta2_xz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 285);

    auto ta2_xz_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 286);

    auto ta2_xz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 287);

    auto ta2_xz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 288);

    auto ta2_xz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 289);

    auto ta2_xz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 290);

    auto ta2_xz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 291);

    auto ta2_xz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 292);

    auto ta2_xz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 293);

    auto ta2_xz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 294);

    auto ta2_xz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 295);

    auto ta2_xz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 296);

    auto ta2_xz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 297);

    auto ta2_xz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 298);

    auto ta2_xz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 299);

    auto ta2_yy_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 300);

    auto ta2_yy_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 301);

    auto ta2_yy_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 302);

    auto ta2_yy_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 303);

    auto ta2_yy_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 304);

    auto ta2_yy_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 305);

    auto ta2_yy_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 306);

    auto ta2_yy_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 307);

    auto ta2_yy_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 308);

    auto ta2_yy_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 309);

    auto ta2_yy_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 310);

    auto ta2_yy_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 311);

    auto ta2_yy_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 312);

    auto ta2_yy_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 313);

    auto ta2_yy_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 314);

    auto ta2_yy_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 315);

    auto ta2_yy_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 316);

    auto ta2_yy_xxy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 317);

    auto ta2_yy_xxy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 318);

    auto ta2_yy_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 320);

    auto ta2_yy_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 321);

    auto ta2_yy_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 322);

    auto ta2_yy_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 323);

    auto ta2_yy_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 325);

    auto ta2_yy_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 327);

    auto ta2_yy_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 328);

    auto ta2_yy_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 329);

    auto ta2_yy_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 330);

    auto ta2_yy_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 331);

    auto ta2_yy_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 333);

    auto ta2_yy_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 334);

    auto ta2_yy_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 336);

    auto ta2_yy_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 337);

    auto ta2_yy_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 338);

    auto ta2_yy_xyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 339);

    auto ta2_yy_xyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 347);

    auto ta2_yy_xyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 348);

    auto ta2_yy_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 350);

    auto ta2_yy_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 352);

    auto ta2_yy_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 354);

    auto ta2_yy_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 355);

    auto ta2_yy_xzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 356);

    auto ta2_yy_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 357);

    auto ta2_yy_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 358);

    auto ta2_yy_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 359);

    auto ta2_yy_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 360);

    auto ta2_yy_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 361);

    auto ta2_yy_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 362);

    auto ta2_yy_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 363);

    auto ta2_yy_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 364);

    auto ta2_yy_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 365);

    auto ta2_yy_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 366);

    auto ta2_yy_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 367);

    auto ta2_yy_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 368);

    auto ta2_yy_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 369);

    auto ta2_yy_yyz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 370);

    auto ta2_yy_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 371);

    auto ta2_yy_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 372);

    auto ta2_yy_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 373);

    auto ta2_yy_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 374);

    auto ta2_yy_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 375);

    auto ta2_yy_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 376);

    auto ta2_yy_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 377);

    auto ta2_yy_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 378);

    auto ta2_yy_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 379);

    auto ta2_yy_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 381);

    auto ta2_yy_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 382);

    auto ta2_yy_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 383);

    auto ta2_yy_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 384);

    auto ta2_yy_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 385);

    auto ta2_yy_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 386);

    auto ta2_yy_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 387);

    auto ta2_yy_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 388);

    auto ta2_yy_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 389);

    auto ta2_yy_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 390);

    auto ta2_yy_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 391);

    auto ta2_yy_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 392);

    auto ta2_yy_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 393);

    auto ta2_yy_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 394);

    auto ta2_yy_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 395);

    auto ta2_yy_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 396);

    auto ta2_yy_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 397);

    auto ta2_yy_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 398);

    auto ta2_yy_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 399);

    auto ta2_yz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 400);

    auto ta2_yz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 401);

    auto ta2_yz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 402);

    auto ta2_yz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 403);

    auto ta2_yz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 404);

    auto ta2_yz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 405);

    auto ta2_yz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 406);

    auto ta2_yz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 407);

    auto ta2_yz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 408);

    auto ta2_yz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 409);

    auto ta2_yz_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 410);

    auto ta2_yz_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 411);

    auto ta2_yz_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 412);

    auto ta2_yz_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 413);

    auto ta2_yz_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 415);

    auto ta2_yz_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 416);

    auto ta2_yz_xxy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 417);

    auto ta2_yz_xxy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 418);

    auto ta2_yz_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 420);

    auto ta2_yz_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 421);

    auto ta2_yz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 422);

    auto ta2_yz_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 423);

    auto ta2_yz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 424);

    auto ta2_yz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 425);

    auto ta2_yz_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 427);

    auto ta2_yz_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 428);

    auto ta2_yz_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 429);

    auto ta2_yz_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 430);

    auto ta2_yz_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 431);

    auto ta2_yz_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 433);

    auto ta2_yz_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 434);

    auto ta2_yz_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 436);

    auto ta2_yz_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 437);

    auto ta2_yz_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 438);

    auto ta2_yz_xyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 439);

    auto ta2_yz_xyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 447);

    auto ta2_yz_xyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 448);

    auto ta2_yz_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 450);

    auto ta2_yz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 452);

    auto ta2_yz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 454);

    auto ta2_yz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 455);

    auto ta2_yz_xzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 456);

    auto ta2_yz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 457);

    auto ta2_yz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 458);

    auto ta2_yz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 459);

    auto ta2_yz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 460);

    auto ta2_yz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 461);

    auto ta2_yz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 462);

    auto ta2_yz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 463);

    auto ta2_yz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 464);

    auto ta2_yz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 465);

    auto ta2_yz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 466);

    auto ta2_yz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 467);

    auto ta2_yz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 468);

    auto ta2_yz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 469);

    auto ta2_yz_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 471);

    auto ta2_yz_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 472);

    auto ta2_yz_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 473);

    auto ta2_yz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 474);

    auto ta2_yz_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 475);

    auto ta2_yz_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 476);

    auto ta2_yz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 477);

    auto ta2_yz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 478);

    auto ta2_yz_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 479);

    auto ta2_yz_yzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 480);

    auto ta2_yz_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 481);

    auto ta2_yz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 482);

    auto ta2_yz_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 483);

    auto ta2_yz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 484);

    auto ta2_yz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 485);

    auto ta2_yz_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 486);

    auto ta2_yz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 487);

    auto ta2_yz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 488);

    auto ta2_yz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 489);

    auto ta2_yz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 490);

    auto ta2_yz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 491);

    auto ta2_yz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 492);

    auto ta2_yz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 493);

    auto ta2_yz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 494);

    auto ta2_yz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 495);

    auto ta2_yz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 496);

    auto ta2_yz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 497);

    auto ta2_yz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 498);

    auto ta2_yz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 499);

    auto ta2_zz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 500);

    auto ta2_zz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 501);

    auto ta2_zz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 502);

    auto ta2_zz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 503);

    auto ta2_zz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 504);

    auto ta2_zz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 505);

    auto ta2_zz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 506);

    auto ta2_zz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 507);

    auto ta2_zz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 508);

    auto ta2_zz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 509);

    auto ta2_zz_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 510);

    auto ta2_zz_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 511);

    auto ta2_zz_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 512);

    auto ta2_zz_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 513);

    auto ta2_zz_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 515);

    auto ta2_zz_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 516);

    auto ta2_zz_xxy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 517);

    auto ta2_zz_xxy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 518);

    auto ta2_zz_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 520);

    auto ta2_zz_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 521);

    auto ta2_zz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 522);

    auto ta2_zz_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 523);

    auto ta2_zz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 524);

    auto ta2_zz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 525);

    auto ta2_zz_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 527);

    auto ta2_zz_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 528);

    auto ta2_zz_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 529);

    auto ta2_zz_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 530);

    auto ta2_zz_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 531);

    auto ta2_zz_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 533);

    auto ta2_zz_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 534);

    auto ta2_zz_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 536);

    auto ta2_zz_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 537);

    auto ta2_zz_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 538);

    auto ta2_zz_xyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 539);

    auto ta2_zz_xyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 547);

    auto ta2_zz_xyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 548);

    auto ta2_zz_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 550);

    auto ta2_zz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 552);

    auto ta2_zz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 554);

    auto ta2_zz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 555);

    auto ta2_zz_xzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 556);

    auto ta2_zz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 557);

    auto ta2_zz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 558);

    auto ta2_zz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 559);

    auto ta2_zz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 560);

    auto ta2_zz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 561);

    auto ta2_zz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 562);

    auto ta2_zz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 563);

    auto ta2_zz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 564);

    auto ta2_zz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 565);

    auto ta2_zz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 566);

    auto ta2_zz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 567);

    auto ta2_zz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 568);

    auto ta2_zz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 569);

    auto ta2_zz_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 571);

    auto ta2_zz_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 572);

    auto ta2_zz_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 573);

    auto ta2_zz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 574);

    auto ta2_zz_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 575);

    auto ta2_zz_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 576);

    auto ta2_zz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 577);

    auto ta2_zz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 578);

    auto ta2_zz_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 579);

    auto ta2_zz_yzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 580);

    auto ta2_zz_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 581);

    auto ta2_zz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 582);

    auto ta2_zz_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 583);

    auto ta2_zz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 584);

    auto ta2_zz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 585);

    auto ta2_zz_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 586);

    auto ta2_zz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 587);

    auto ta2_zz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 588);

    auto ta2_zz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 589);

    auto ta2_zz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 590);

    auto ta2_zz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 591);

    auto ta2_zz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 592);

    auto ta2_zz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 593);

    auto ta2_zz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 594);

    auto ta2_zz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 595);

    auto ta2_zz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 596);

    auto ta2_zz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 597);

    auto ta2_zz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 598);

    auto ta2_zz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 599);

    // Set up 0-10 components of targeted buffer : GF

    auto ta2_xx_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf);

    auto ta2_xx_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 1);

    auto ta2_xx_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 2);

    auto ta2_xx_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 3);

    auto ta2_xx_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 4);

    auto ta2_xx_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 5);

    auto ta2_xx_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 6);

    auto ta2_xx_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 7);

    auto ta2_xx_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 8);

    auto ta2_xx_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 9);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_1, ta1_x_xxx_yyy_1, ta1_x_xxx_yyz_1, ta1_x_xxx_yzz_1, ta1_x_xxx_zzz_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xx_yyz_0, ta2_xx_xx_yyz_1, ta2_xx_xx_yzz_0, ta2_xx_xx_yzz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xxx_xx_0, ta2_xx_xxx_xx_1, ta2_xx_xxx_xxx_0, ta2_xx_xxx_xxx_1, ta2_xx_xxx_xxy_0, ta2_xx_xxx_xxy_1, ta2_xx_xxx_xxz_0, ta2_xx_xxx_xxz_1, ta2_xx_xxx_xy_0, ta2_xx_xxx_xy_1, ta2_xx_xxx_xyy_0, ta2_xx_xxx_xyy_1, ta2_xx_xxx_xyz_0, ta2_xx_xxx_xyz_1, ta2_xx_xxx_xz_0, ta2_xx_xxx_xz_1, ta2_xx_xxx_xzz_0, ta2_xx_xxx_xzz_1, ta2_xx_xxx_yy_0, ta2_xx_xxx_yy_1, ta2_xx_xxx_yyy_0, ta2_xx_xxx_yyy_1, ta2_xx_xxx_yyz_0, ta2_xx_xxx_yyz_1, ta2_xx_xxx_yz_0, ta2_xx_xxx_yz_1, ta2_xx_xxx_yzz_0, ta2_xx_xxx_yzz_1, ta2_xx_xxx_zz_0, ta2_xx_xxx_zz_1, ta2_xx_xxx_zzz_0, ta2_xx_xxx_zzz_1, ta2_xx_xxxx_xxx_0, ta2_xx_xxxx_xxy_0, ta2_xx_xxxx_xxz_0, ta2_xx_xxxx_xyy_0, ta2_xx_xxxx_xyz_0, ta2_xx_xxxx_xzz_0, ta2_xx_xxxx_yyy_0, ta2_xx_xxxx_yyz_0, ta2_xx_xxxx_yzz_0, ta2_xx_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxx_xxx_0[i] = 3.0 * ta2_xx_xx_xxx_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxx_1[i] * fe_0 + 3.0 * ta2_xx_xxx_xx_0[i] * fe_0 - 3.0 * ta2_xx_xxx_xx_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxx_1[i] + ta2_xx_xxx_xxx_0[i] * pa_x[i] - ta2_xx_xxx_xxx_1[i] * pc_x[i];

        ta2_xx_xxxx_xxy_0[i] = 3.0 * ta2_xx_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxy_1[i] * fe_0 + 2.0 * ta2_xx_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxy_1[i] + ta2_xx_xxx_xxy_0[i] * pa_x[i] - ta2_xx_xxx_xxy_1[i] * pc_x[i];

        ta2_xx_xxxx_xxz_0[i] = 3.0 * ta2_xx_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxz_1[i] * fe_0 + 2.0 * ta2_xx_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxz_1[i] + ta2_xx_xxx_xxz_0[i] * pa_x[i] - ta2_xx_xxx_xxz_1[i] * pc_x[i];

        ta2_xx_xxxx_xyy_0[i] = 3.0 * ta2_xx_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyy_1[i] * fe_0 + ta2_xx_xxx_yy_0[i] * fe_0 - ta2_xx_xxx_yy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyy_1[i] + ta2_xx_xxx_xyy_0[i] * pa_x[i] - ta2_xx_xxx_xyy_1[i] * pc_x[i];

        ta2_xx_xxxx_xyz_0[i] = 3.0 * ta2_xx_xx_xyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyz_1[i] * fe_0 + ta2_xx_xxx_yz_0[i] * fe_0 - ta2_xx_xxx_yz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyz_1[i] + ta2_xx_xxx_xyz_0[i] * pa_x[i] - ta2_xx_xxx_xyz_1[i] * pc_x[i];

        ta2_xx_xxxx_xzz_0[i] = 3.0 * ta2_xx_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xzz_1[i] * fe_0 + ta2_xx_xxx_zz_0[i] * fe_0 - ta2_xx_xxx_zz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xzz_1[i] + ta2_xx_xxx_xzz_0[i] * pa_x[i] - ta2_xx_xxx_xzz_1[i] * pc_x[i];

        ta2_xx_xxxx_yyy_0[i] = 3.0 * ta2_xx_xx_yyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_yyy_1[i] + ta2_xx_xxx_yyy_0[i] * pa_x[i] - ta2_xx_xxx_yyy_1[i] * pc_x[i];

        ta2_xx_xxxx_yyz_0[i] = 3.0 * ta2_xx_xx_yyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yyz_1[i] + ta2_xx_xxx_yyz_0[i] * pa_x[i] - ta2_xx_xxx_yyz_1[i] * pc_x[i];

        ta2_xx_xxxx_yzz_0[i] = 3.0 * ta2_xx_xx_yzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yzz_1[i] + ta2_xx_xxx_yzz_0[i] * pa_x[i] - ta2_xx_xxx_yzz_1[i] * pc_x[i];

        ta2_xx_xxxx_zzz_0[i] = 3.0 * ta2_xx_xx_zzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_zzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_zzz_1[i] + ta2_xx_xxx_zzz_0[i] * pa_x[i] - ta2_xx_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto ta2_xx_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 10);

    auto ta2_xx_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 11);

    auto ta2_xx_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 12);

    auto ta2_xx_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 13);

    auto ta2_xx_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 14);

    auto ta2_xx_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 15);

    auto ta2_xx_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 16);

    auto ta2_xx_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 17);

    auto ta2_xx_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 18);

    auto ta2_xx_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 19);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_xxx_xx_0, ta2_xx_xxx_xx_1, ta2_xx_xxx_xxx_0, ta2_xx_xxx_xxx_1, ta2_xx_xxx_xxy_0, ta2_xx_xxx_xxy_1, ta2_xx_xxx_xxz_0, ta2_xx_xxx_xxz_1, ta2_xx_xxx_xy_0, ta2_xx_xxx_xy_1, ta2_xx_xxx_xyy_0, ta2_xx_xxx_xyy_1, ta2_xx_xxx_xyz_0, ta2_xx_xxx_xyz_1, ta2_xx_xxx_xz_0, ta2_xx_xxx_xz_1, ta2_xx_xxx_xzz_0, ta2_xx_xxx_xzz_1, ta2_xx_xxx_yy_0, ta2_xx_xxx_yy_1, ta2_xx_xxx_yyy_0, ta2_xx_xxx_yyy_1, ta2_xx_xxx_yyz_0, ta2_xx_xxx_yyz_1, ta2_xx_xxx_yz_0, ta2_xx_xxx_yz_1, ta2_xx_xxx_yzz_0, ta2_xx_xxx_yzz_1, ta2_xx_xxx_zz_0, ta2_xx_xxx_zz_1, ta2_xx_xxx_zzz_0, ta2_xx_xxx_zzz_1, ta2_xx_xxxy_xxx_0, ta2_xx_xxxy_xxy_0, ta2_xx_xxxy_xxz_0, ta2_xx_xxxy_xyy_0, ta2_xx_xxxy_xyz_0, ta2_xx_xxxy_xzz_0, ta2_xx_xxxy_yyy_0, ta2_xx_xxxy_yyz_0, ta2_xx_xxxy_yzz_0, ta2_xx_xxxy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxy_xxx_0[i] = ta2_xx_xxx_xxx_0[i] * pa_y[i] - ta2_xx_xxx_xxx_1[i] * pc_y[i];

        ta2_xx_xxxy_xxy_0[i] = ta2_xx_xxx_xx_0[i] * fe_0 - ta2_xx_xxx_xx_1[i] * fe_0 + ta2_xx_xxx_xxy_0[i] * pa_y[i] - ta2_xx_xxx_xxy_1[i] * pc_y[i];

        ta2_xx_xxxy_xxz_0[i] = ta2_xx_xxx_xxz_0[i] * pa_y[i] - ta2_xx_xxx_xxz_1[i] * pc_y[i];

        ta2_xx_xxxy_xyy_0[i] = 2.0 * ta2_xx_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xy_1[i] * fe_0 + ta2_xx_xxx_xyy_0[i] * pa_y[i] - ta2_xx_xxx_xyy_1[i] * pc_y[i];

        ta2_xx_xxxy_xyz_0[i] = ta2_xx_xxx_xz_0[i] * fe_0 - ta2_xx_xxx_xz_1[i] * fe_0 + ta2_xx_xxx_xyz_0[i] * pa_y[i] - ta2_xx_xxx_xyz_1[i] * pc_y[i];

        ta2_xx_xxxy_xzz_0[i] = ta2_xx_xxx_xzz_0[i] * pa_y[i] - ta2_xx_xxx_xzz_1[i] * pc_y[i];

        ta2_xx_xxxy_yyy_0[i] = 3.0 * ta2_xx_xxx_yy_0[i] * fe_0 - 3.0 * ta2_xx_xxx_yy_1[i] * fe_0 + ta2_xx_xxx_yyy_0[i] * pa_y[i] - ta2_xx_xxx_yyy_1[i] * pc_y[i];

        ta2_xx_xxxy_yyz_0[i] = 2.0 * ta2_xx_xxx_yz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_yz_1[i] * fe_0 + ta2_xx_xxx_yyz_0[i] * pa_y[i] - ta2_xx_xxx_yyz_1[i] * pc_y[i];

        ta2_xx_xxxy_yzz_0[i] = ta2_xx_xxx_zz_0[i] * fe_0 - ta2_xx_xxx_zz_1[i] * fe_0 + ta2_xx_xxx_yzz_0[i] * pa_y[i] - ta2_xx_xxx_yzz_1[i] * pc_y[i];

        ta2_xx_xxxy_zzz_0[i] = ta2_xx_xxx_zzz_0[i] * pa_y[i] - ta2_xx_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto ta2_xx_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 20);

    auto ta2_xx_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 21);

    auto ta2_xx_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 22);

    auto ta2_xx_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 23);

    auto ta2_xx_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 24);

    auto ta2_xx_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 25);

    auto ta2_xx_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 26);

    auto ta2_xx_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 27);

    auto ta2_xx_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 28);

    auto ta2_xx_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 29);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_xxx_xx_0, ta2_xx_xxx_xx_1, ta2_xx_xxx_xxx_0, ta2_xx_xxx_xxx_1, ta2_xx_xxx_xxy_0, ta2_xx_xxx_xxy_1, ta2_xx_xxx_xxz_0, ta2_xx_xxx_xxz_1, ta2_xx_xxx_xy_0, ta2_xx_xxx_xy_1, ta2_xx_xxx_xyy_0, ta2_xx_xxx_xyy_1, ta2_xx_xxx_xyz_0, ta2_xx_xxx_xyz_1, ta2_xx_xxx_xz_0, ta2_xx_xxx_xz_1, ta2_xx_xxx_xzz_0, ta2_xx_xxx_xzz_1, ta2_xx_xxx_yy_0, ta2_xx_xxx_yy_1, ta2_xx_xxx_yyy_0, ta2_xx_xxx_yyy_1, ta2_xx_xxx_yyz_0, ta2_xx_xxx_yyz_1, ta2_xx_xxx_yz_0, ta2_xx_xxx_yz_1, ta2_xx_xxx_yzz_0, ta2_xx_xxx_yzz_1, ta2_xx_xxx_zz_0, ta2_xx_xxx_zz_1, ta2_xx_xxx_zzz_0, ta2_xx_xxx_zzz_1, ta2_xx_xxxz_xxx_0, ta2_xx_xxxz_xxy_0, ta2_xx_xxxz_xxz_0, ta2_xx_xxxz_xyy_0, ta2_xx_xxxz_xyz_0, ta2_xx_xxxz_xzz_0, ta2_xx_xxxz_yyy_0, ta2_xx_xxxz_yyz_0, ta2_xx_xxxz_yzz_0, ta2_xx_xxxz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxz_xxx_0[i] = ta2_xx_xxx_xxx_0[i] * pa_z[i] - ta2_xx_xxx_xxx_1[i] * pc_z[i];

        ta2_xx_xxxz_xxy_0[i] = ta2_xx_xxx_xxy_0[i] * pa_z[i] - ta2_xx_xxx_xxy_1[i] * pc_z[i];

        ta2_xx_xxxz_xxz_0[i] = ta2_xx_xxx_xx_0[i] * fe_0 - ta2_xx_xxx_xx_1[i] * fe_0 + ta2_xx_xxx_xxz_0[i] * pa_z[i] - ta2_xx_xxx_xxz_1[i] * pc_z[i];

        ta2_xx_xxxz_xyy_0[i] = ta2_xx_xxx_xyy_0[i] * pa_z[i] - ta2_xx_xxx_xyy_1[i] * pc_z[i];

        ta2_xx_xxxz_xyz_0[i] = ta2_xx_xxx_xy_0[i] * fe_0 - ta2_xx_xxx_xy_1[i] * fe_0 + ta2_xx_xxx_xyz_0[i] * pa_z[i] - ta2_xx_xxx_xyz_1[i] * pc_z[i];

        ta2_xx_xxxz_xzz_0[i] = 2.0 * ta2_xx_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xz_1[i] * fe_0 + ta2_xx_xxx_xzz_0[i] * pa_z[i] - ta2_xx_xxx_xzz_1[i] * pc_z[i];

        ta2_xx_xxxz_yyy_0[i] = ta2_xx_xxx_yyy_0[i] * pa_z[i] - ta2_xx_xxx_yyy_1[i] * pc_z[i];

        ta2_xx_xxxz_yyz_0[i] = ta2_xx_xxx_yy_0[i] * fe_0 - ta2_xx_xxx_yy_1[i] * fe_0 + ta2_xx_xxx_yyz_0[i] * pa_z[i] - ta2_xx_xxx_yyz_1[i] * pc_z[i];

        ta2_xx_xxxz_yzz_0[i] = 2.0 * ta2_xx_xxx_yz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_yz_1[i] * fe_0 + ta2_xx_xxx_yzz_0[i] * pa_z[i] - ta2_xx_xxx_yzz_1[i] * pc_z[i];

        ta2_xx_xxxz_zzz_0[i] = 3.0 * ta2_xx_xxx_zz_0[i] * fe_0 - 3.0 * ta2_xx_xxx_zz_1[i] * fe_0 + ta2_xx_xxx_zzz_0[i] * pa_z[i] - ta2_xx_xxx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto ta2_xx_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 30);

    auto ta2_xx_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 31);

    auto ta2_xx_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 32);

    auto ta2_xx_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 33);

    auto ta2_xx_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 34);

    auto ta2_xx_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 35);

    auto ta2_xx_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 36);

    auto ta2_xx_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 37);

    auto ta2_xx_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 38);

    auto ta2_xx_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 39);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyy_yyy_1, ta1_x_xyy_yyz_1, ta1_x_xyy_yzz_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_zzz_0, ta2_xx_xx_zzz_1, ta2_xx_xxy_xx_0, ta2_xx_xxy_xx_1, ta2_xx_xxy_xxx_0, ta2_xx_xxy_xxx_1, ta2_xx_xxy_xxy_0, ta2_xx_xxy_xxy_1, ta2_xx_xxy_xxz_0, ta2_xx_xxy_xxz_1, ta2_xx_xxy_xy_0, ta2_xx_xxy_xy_1, ta2_xx_xxy_xyy_0, ta2_xx_xxy_xyy_1, ta2_xx_xxy_xyz_0, ta2_xx_xxy_xyz_1, ta2_xx_xxy_xz_0, ta2_xx_xxy_xz_1, ta2_xx_xxy_xzz_0, ta2_xx_xxy_xzz_1, ta2_xx_xxy_zzz_0, ta2_xx_xxy_zzz_1, ta2_xx_xxyy_xxx_0, ta2_xx_xxyy_xxy_0, ta2_xx_xxyy_xxz_0, ta2_xx_xxyy_xyy_0, ta2_xx_xxyy_xyz_0, ta2_xx_xxyy_xzz_0, ta2_xx_xxyy_yyy_0, ta2_xx_xxyy_yyz_0, ta2_xx_xxyy_yzz_0, ta2_xx_xxyy_zzz_0, ta2_xx_xyy_yyy_0, ta2_xx_xyy_yyy_1, ta2_xx_xyy_yyz_0, ta2_xx_xyy_yyz_1, ta2_xx_xyy_yzz_0, ta2_xx_xyy_yzz_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyy_xxx_0[i] = ta2_xx_xx_xxx_0[i] * fe_0 - ta2_xx_xx_xxx_1[i] * fe_0 + ta2_xx_xxy_xxx_0[i] * pa_y[i] - ta2_xx_xxy_xxx_1[i] * pc_y[i];

        ta2_xx_xxyy_xxy_0[i] = ta2_xx_xx_xxy_0[i] * fe_0 - ta2_xx_xx_xxy_1[i] * fe_0 + ta2_xx_xxy_xx_0[i] * fe_0 - ta2_xx_xxy_xx_1[i] * fe_0 + ta2_xx_xxy_xxy_0[i] * pa_y[i] - ta2_xx_xxy_xxy_1[i] * pc_y[i];

        ta2_xx_xxyy_xxz_0[i] = ta2_xx_xx_xxz_0[i] * fe_0 - ta2_xx_xx_xxz_1[i] * fe_0 + ta2_xx_xxy_xxz_0[i] * pa_y[i] - ta2_xx_xxy_xxz_1[i] * pc_y[i];

        ta2_xx_xxyy_xyy_0[i] = ta2_xx_xx_xyy_0[i] * fe_0 - ta2_xx_xx_xyy_1[i] * fe_0 + 2.0 * ta2_xx_xxy_xy_0[i] * fe_0 - 2.0 * ta2_xx_xxy_xy_1[i] * fe_0 + ta2_xx_xxy_xyy_0[i] * pa_y[i] - ta2_xx_xxy_xyy_1[i] * pc_y[i];

        ta2_xx_xxyy_xyz_0[i] = ta2_xx_xx_xyz_0[i] * fe_0 - ta2_xx_xx_xyz_1[i] * fe_0 + ta2_xx_xxy_xz_0[i] * fe_0 - ta2_xx_xxy_xz_1[i] * fe_0 + ta2_xx_xxy_xyz_0[i] * pa_y[i] - ta2_xx_xxy_xyz_1[i] * pc_y[i];

        ta2_xx_xxyy_xzz_0[i] = ta2_xx_xx_xzz_0[i] * fe_0 - ta2_xx_xx_xzz_1[i] * fe_0 + ta2_xx_xxy_xzz_0[i] * pa_y[i] - ta2_xx_xxy_xzz_1[i] * pc_y[i];

        ta2_xx_xxyy_yyy_0[i] = ta2_xx_yy_yyy_0[i] * fe_0 - ta2_xx_yy_yyy_1[i] * fe_0 + 2.0 * ta1_x_xyy_yyy_1[i] + ta2_xx_xyy_yyy_0[i] * pa_x[i] - ta2_xx_xyy_yyy_1[i] * pc_x[i];

        ta2_xx_xxyy_yyz_0[i] = ta2_xx_yy_yyz_0[i] * fe_0 - ta2_xx_yy_yyz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yyz_1[i] + ta2_xx_xyy_yyz_0[i] * pa_x[i] - ta2_xx_xyy_yyz_1[i] * pc_x[i];

        ta2_xx_xxyy_yzz_0[i] = ta2_xx_yy_yzz_0[i] * fe_0 - ta2_xx_yy_yzz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yzz_1[i] + ta2_xx_xyy_yzz_0[i] * pa_x[i] - ta2_xx_xyy_yzz_1[i] * pc_x[i];

        ta2_xx_xxyy_zzz_0[i] = ta2_xx_xx_zzz_0[i] * fe_0 - ta2_xx_xx_zzz_1[i] * fe_0 + ta2_xx_xxy_zzz_0[i] * pa_y[i] - ta2_xx_xxy_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto ta2_xx_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 40);

    auto ta2_xx_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 41);

    auto ta2_xx_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 42);

    auto ta2_xx_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 43);

    auto ta2_xx_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 44);

    auto ta2_xx_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 45);

    auto ta2_xx_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 46);

    auto ta2_xx_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 47);

    auto ta2_xx_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 48);

    auto ta2_xx_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 49);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_xxy_xxy_0, ta2_xx_xxy_xxy_1, ta2_xx_xxy_xyy_0, ta2_xx_xxy_xyy_1, ta2_xx_xxy_yyy_0, ta2_xx_xxy_yyy_1, ta2_xx_xxyz_xxx_0, ta2_xx_xxyz_xxy_0, ta2_xx_xxyz_xxz_0, ta2_xx_xxyz_xyy_0, ta2_xx_xxyz_xyz_0, ta2_xx_xxyz_xzz_0, ta2_xx_xxyz_yyy_0, ta2_xx_xxyz_yyz_0, ta2_xx_xxyz_yzz_0, ta2_xx_xxyz_zzz_0, ta2_xx_xxz_xxx_0, ta2_xx_xxz_xxx_1, ta2_xx_xxz_xxz_0, ta2_xx_xxz_xxz_1, ta2_xx_xxz_xyz_0, ta2_xx_xxz_xyz_1, ta2_xx_xxz_xz_0, ta2_xx_xxz_xz_1, ta2_xx_xxz_xzz_0, ta2_xx_xxz_xzz_1, ta2_xx_xxz_yyz_0, ta2_xx_xxz_yyz_1, ta2_xx_xxz_yz_0, ta2_xx_xxz_yz_1, ta2_xx_xxz_yzz_0, ta2_xx_xxz_yzz_1, ta2_xx_xxz_zz_0, ta2_xx_xxz_zz_1, ta2_xx_xxz_zzz_0, ta2_xx_xxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyz_xxx_0[i] = ta2_xx_xxz_xxx_0[i] * pa_y[i] - ta2_xx_xxz_xxx_1[i] * pc_y[i];

        ta2_xx_xxyz_xxy_0[i] = ta2_xx_xxy_xxy_0[i] * pa_z[i] - ta2_xx_xxy_xxy_1[i] * pc_z[i];

        ta2_xx_xxyz_xxz_0[i] = ta2_xx_xxz_xxz_0[i] * pa_y[i] - ta2_xx_xxz_xxz_1[i] * pc_y[i];

        ta2_xx_xxyz_xyy_0[i] = ta2_xx_xxy_xyy_0[i] * pa_z[i] - ta2_xx_xxy_xyy_1[i] * pc_z[i];

        ta2_xx_xxyz_xyz_0[i] = ta2_xx_xxz_xz_0[i] * fe_0 - ta2_xx_xxz_xz_1[i] * fe_0 + ta2_xx_xxz_xyz_0[i] * pa_y[i] - ta2_xx_xxz_xyz_1[i] * pc_y[i];

        ta2_xx_xxyz_xzz_0[i] = ta2_xx_xxz_xzz_0[i] * pa_y[i] - ta2_xx_xxz_xzz_1[i] * pc_y[i];

        ta2_xx_xxyz_yyy_0[i] = ta2_xx_xxy_yyy_0[i] * pa_z[i] - ta2_xx_xxy_yyy_1[i] * pc_z[i];

        ta2_xx_xxyz_yyz_0[i] = 2.0 * ta2_xx_xxz_yz_0[i] * fe_0 - 2.0 * ta2_xx_xxz_yz_1[i] * fe_0 + ta2_xx_xxz_yyz_0[i] * pa_y[i] - ta2_xx_xxz_yyz_1[i] * pc_y[i];

        ta2_xx_xxyz_yzz_0[i] = ta2_xx_xxz_zz_0[i] * fe_0 - ta2_xx_xxz_zz_1[i] * fe_0 + ta2_xx_xxz_yzz_0[i] * pa_y[i] - ta2_xx_xxz_yzz_1[i] * pc_y[i];

        ta2_xx_xxyz_zzz_0[i] = ta2_xx_xxz_zzz_0[i] * pa_y[i] - ta2_xx_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto ta2_xx_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 50);

    auto ta2_xx_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 51);

    auto ta2_xx_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 52);

    auto ta2_xx_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 53);

    auto ta2_xx_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 54);

    auto ta2_xx_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 55);

    auto ta2_xx_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 56);

    auto ta2_xx_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 57);

    auto ta2_xx_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 58);

    auto ta2_xx_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 59);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xzz_yyz_1, ta1_x_xzz_yzz_1, ta1_x_xzz_zzz_1, ta2_xx_xx_xxx_0, ta2_xx_xx_xxx_1, ta2_xx_xx_xxy_0, ta2_xx_xx_xxy_1, ta2_xx_xx_xxz_0, ta2_xx_xx_xxz_1, ta2_xx_xx_xyy_0, ta2_xx_xx_xyy_1, ta2_xx_xx_xyz_0, ta2_xx_xx_xyz_1, ta2_xx_xx_xzz_0, ta2_xx_xx_xzz_1, ta2_xx_xx_yyy_0, ta2_xx_xx_yyy_1, ta2_xx_xxz_xx_0, ta2_xx_xxz_xx_1, ta2_xx_xxz_xxx_0, ta2_xx_xxz_xxx_1, ta2_xx_xxz_xxy_0, ta2_xx_xxz_xxy_1, ta2_xx_xxz_xxz_0, ta2_xx_xxz_xxz_1, ta2_xx_xxz_xy_0, ta2_xx_xxz_xy_1, ta2_xx_xxz_xyy_0, ta2_xx_xxz_xyy_1, ta2_xx_xxz_xyz_0, ta2_xx_xxz_xyz_1, ta2_xx_xxz_xz_0, ta2_xx_xxz_xz_1, ta2_xx_xxz_xzz_0, ta2_xx_xxz_xzz_1, ta2_xx_xxz_yyy_0, ta2_xx_xxz_yyy_1, ta2_xx_xxzz_xxx_0, ta2_xx_xxzz_xxy_0, ta2_xx_xxzz_xxz_0, ta2_xx_xxzz_xyy_0, ta2_xx_xxzz_xyz_0, ta2_xx_xxzz_xzz_0, ta2_xx_xxzz_yyy_0, ta2_xx_xxzz_yyz_0, ta2_xx_xxzz_yzz_0, ta2_xx_xxzz_zzz_0, ta2_xx_xzz_yyz_0, ta2_xx_xzz_yyz_1, ta2_xx_xzz_yzz_0, ta2_xx_xzz_yzz_1, ta2_xx_xzz_zzz_0, ta2_xx_xzz_zzz_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxzz_xxx_0[i] = ta2_xx_xx_xxx_0[i] * fe_0 - ta2_xx_xx_xxx_1[i] * fe_0 + ta2_xx_xxz_xxx_0[i] * pa_z[i] - ta2_xx_xxz_xxx_1[i] * pc_z[i];

        ta2_xx_xxzz_xxy_0[i] = ta2_xx_xx_xxy_0[i] * fe_0 - ta2_xx_xx_xxy_1[i] * fe_0 + ta2_xx_xxz_xxy_0[i] * pa_z[i] - ta2_xx_xxz_xxy_1[i] * pc_z[i];

        ta2_xx_xxzz_xxz_0[i] = ta2_xx_xx_xxz_0[i] * fe_0 - ta2_xx_xx_xxz_1[i] * fe_0 + ta2_xx_xxz_xx_0[i] * fe_0 - ta2_xx_xxz_xx_1[i] * fe_0 + ta2_xx_xxz_xxz_0[i] * pa_z[i] - ta2_xx_xxz_xxz_1[i] * pc_z[i];

        ta2_xx_xxzz_xyy_0[i] = ta2_xx_xx_xyy_0[i] * fe_0 - ta2_xx_xx_xyy_1[i] * fe_0 + ta2_xx_xxz_xyy_0[i] * pa_z[i] - ta2_xx_xxz_xyy_1[i] * pc_z[i];

        ta2_xx_xxzz_xyz_0[i] = ta2_xx_xx_xyz_0[i] * fe_0 - ta2_xx_xx_xyz_1[i] * fe_0 + ta2_xx_xxz_xy_0[i] * fe_0 - ta2_xx_xxz_xy_1[i] * fe_0 + ta2_xx_xxz_xyz_0[i] * pa_z[i] - ta2_xx_xxz_xyz_1[i] * pc_z[i];

        ta2_xx_xxzz_xzz_0[i] = ta2_xx_xx_xzz_0[i] * fe_0 - ta2_xx_xx_xzz_1[i] * fe_0 + 2.0 * ta2_xx_xxz_xz_0[i] * fe_0 - 2.0 * ta2_xx_xxz_xz_1[i] * fe_0 + ta2_xx_xxz_xzz_0[i] * pa_z[i] - ta2_xx_xxz_xzz_1[i] * pc_z[i];

        ta2_xx_xxzz_yyy_0[i] = ta2_xx_xx_yyy_0[i] * fe_0 - ta2_xx_xx_yyy_1[i] * fe_0 + ta2_xx_xxz_yyy_0[i] * pa_z[i] - ta2_xx_xxz_yyy_1[i] * pc_z[i];

        ta2_xx_xxzz_yyz_0[i] = ta2_xx_zz_yyz_0[i] * fe_0 - ta2_xx_zz_yyz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yyz_1[i] + ta2_xx_xzz_yyz_0[i] * pa_x[i] - ta2_xx_xzz_yyz_1[i] * pc_x[i];

        ta2_xx_xxzz_yzz_0[i] = ta2_xx_zz_yzz_0[i] * fe_0 - ta2_xx_zz_yzz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yzz_1[i] + ta2_xx_xzz_yzz_0[i] * pa_x[i] - ta2_xx_xzz_yzz_1[i] * pc_x[i];

        ta2_xx_xxzz_zzz_0[i] = ta2_xx_zz_zzz_0[i] * fe_0 - ta2_xx_zz_zzz_1[i] * fe_0 + 2.0 * ta1_x_xzz_zzz_1[i] + ta2_xx_xzz_zzz_0[i] * pa_x[i] - ta2_xx_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto ta2_xx_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 60);

    auto ta2_xx_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 61);

    auto ta2_xx_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 62);

    auto ta2_xx_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 63);

    auto ta2_xx_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 64);

    auto ta2_xx_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 65);

    auto ta2_xx_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 66);

    auto ta2_xx_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 67);

    auto ta2_xx_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 68);

    auto ta2_xx_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 69);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_yyy_xxy_1, ta1_x_yyy_xyy_1, ta1_x_yyy_xyz_1, ta1_x_yyy_yyy_1, ta1_x_yyy_yyz_1, ta1_x_yyy_yzz_1, ta1_x_yyy_zzz_1, ta2_xx_xy_xxx_0, ta2_xx_xy_xxx_1, ta2_xx_xy_xxz_0, ta2_xx_xy_xxz_1, ta2_xx_xy_xzz_0, ta2_xx_xy_xzz_1, ta2_xx_xyy_xxx_0, ta2_xx_xyy_xxx_1, ta2_xx_xyy_xxz_0, ta2_xx_xyy_xxz_1, ta2_xx_xyy_xzz_0, ta2_xx_xyy_xzz_1, ta2_xx_xyyy_xxx_0, ta2_xx_xyyy_xxy_0, ta2_xx_xyyy_xxz_0, ta2_xx_xyyy_xyy_0, ta2_xx_xyyy_xyz_0, ta2_xx_xyyy_xzz_0, ta2_xx_xyyy_yyy_0, ta2_xx_xyyy_yyz_0, ta2_xx_xyyy_yzz_0, ta2_xx_xyyy_zzz_0, ta2_xx_yyy_xxy_0, ta2_xx_yyy_xxy_1, ta2_xx_yyy_xy_0, ta2_xx_yyy_xy_1, ta2_xx_yyy_xyy_0, ta2_xx_yyy_xyy_1, ta2_xx_yyy_xyz_0, ta2_xx_yyy_xyz_1, ta2_xx_yyy_yy_0, ta2_xx_yyy_yy_1, ta2_xx_yyy_yyy_0, ta2_xx_yyy_yyy_1, ta2_xx_yyy_yyz_0, ta2_xx_yyy_yyz_1, ta2_xx_yyy_yz_0, ta2_xx_yyy_yz_1, ta2_xx_yyy_yzz_0, ta2_xx_yyy_yzz_1, ta2_xx_yyy_zzz_0, ta2_xx_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyy_xxx_0[i] = 2.0 * ta2_xx_xy_xxx_0[i] * fe_0 - 2.0 * ta2_xx_xy_xxx_1[i] * fe_0 + ta2_xx_xyy_xxx_0[i] * pa_y[i] - ta2_xx_xyy_xxx_1[i] * pc_y[i];

        ta2_xx_xyyy_xxy_0[i] = 2.0 * ta2_xx_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xx_yyy_xy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxy_1[i] + ta2_xx_yyy_xxy_0[i] * pa_x[i] - ta2_xx_yyy_xxy_1[i] * pc_x[i];

        ta2_xx_xyyy_xxz_0[i] = 2.0 * ta2_xx_xy_xxz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xxz_1[i] * fe_0 + ta2_xx_xyy_xxz_0[i] * pa_y[i] - ta2_xx_xyy_xxz_1[i] * pc_y[i];

        ta2_xx_xyyy_xyy_0[i] = ta2_xx_yyy_yy_0[i] * fe_0 - ta2_xx_yyy_yy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyy_1[i] + ta2_xx_yyy_xyy_0[i] * pa_x[i] - ta2_xx_yyy_xyy_1[i] * pc_x[i];

        ta2_xx_xyyy_xyz_0[i] = ta2_xx_yyy_yz_0[i] * fe_0 - ta2_xx_yyy_yz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyz_1[i] + ta2_xx_yyy_xyz_0[i] * pa_x[i] - ta2_xx_yyy_xyz_1[i] * pc_x[i];

        ta2_xx_xyyy_xzz_0[i] = 2.0 * ta2_xx_xy_xzz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xzz_1[i] * fe_0 + ta2_xx_xyy_xzz_0[i] * pa_y[i] - ta2_xx_xyy_xzz_1[i] * pc_y[i];

        ta2_xx_xyyy_yyy_0[i] = 2.0 * ta1_x_yyy_yyy_1[i] + ta2_xx_yyy_yyy_0[i] * pa_x[i] - ta2_xx_yyy_yyy_1[i] * pc_x[i];

        ta2_xx_xyyy_yyz_0[i] = 2.0 * ta1_x_yyy_yyz_1[i] + ta2_xx_yyy_yyz_0[i] * pa_x[i] - ta2_xx_yyy_yyz_1[i] * pc_x[i];

        ta2_xx_xyyy_yzz_0[i] = 2.0 * ta1_x_yyy_yzz_1[i] + ta2_xx_yyy_yzz_0[i] * pa_x[i] - ta2_xx_yyy_yzz_1[i] * pc_x[i];

        ta2_xx_xyyy_zzz_0[i] = 2.0 * ta1_x_yyy_zzz_1[i] + ta2_xx_yyy_zzz_0[i] * pa_x[i] - ta2_xx_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto ta2_xx_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 70);

    auto ta2_xx_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 71);

    auto ta2_xx_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 72);

    auto ta2_xx_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 73);

    auto ta2_xx_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 74);

    auto ta2_xx_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 75);

    auto ta2_xx_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 76);

    auto ta2_xx_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 77);

    auto ta2_xx_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 78);

    auto ta2_xx_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 79);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_yyz_yyz_1, ta1_x_yyz_yzz_1, ta1_x_yyz_zzz_1, ta2_xx_xyy_xxx_0, ta2_xx_xyy_xxx_1, ta2_xx_xyy_xxy_0, ta2_xx_xyy_xxy_1, ta2_xx_xyy_xy_0, ta2_xx_xyy_xy_1, ta2_xx_xyy_xyy_0, ta2_xx_xyy_xyy_1, ta2_xx_xyy_xyz_0, ta2_xx_xyy_xyz_1, ta2_xx_xyy_yyy_0, ta2_xx_xyy_yyy_1, ta2_xx_xyyz_xxx_0, ta2_xx_xyyz_xxy_0, ta2_xx_xyyz_xxz_0, ta2_xx_xyyz_xyy_0, ta2_xx_xyyz_xyz_0, ta2_xx_xyyz_xzz_0, ta2_xx_xyyz_yyy_0, ta2_xx_xyyz_yyz_0, ta2_xx_xyyz_yzz_0, ta2_xx_xyyz_zzz_0, ta2_xx_xyz_xxz_0, ta2_xx_xyz_xxz_1, ta2_xx_xyz_xzz_0, ta2_xx_xyz_xzz_1, ta2_xx_xz_xxz_0, ta2_xx_xz_xxz_1, ta2_xx_xz_xzz_0, ta2_xx_xz_xzz_1, ta2_xx_yyz_yyz_0, ta2_xx_yyz_yyz_1, ta2_xx_yyz_yzz_0, ta2_xx_yyz_yzz_1, ta2_xx_yyz_zzz_0, ta2_xx_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyz_xxx_0[i] = ta2_xx_xyy_xxx_0[i] * pa_z[i] - ta2_xx_xyy_xxx_1[i] * pc_z[i];

        ta2_xx_xyyz_xxy_0[i] = ta2_xx_xyy_xxy_0[i] * pa_z[i] - ta2_xx_xyy_xxy_1[i] * pc_z[i];

        ta2_xx_xyyz_xxz_0[i] = ta2_xx_xz_xxz_0[i] * fe_0 - ta2_xx_xz_xxz_1[i] * fe_0 + ta2_xx_xyz_xxz_0[i] * pa_y[i] - ta2_xx_xyz_xxz_1[i] * pc_y[i];

        ta2_xx_xyyz_xyy_0[i] = ta2_xx_xyy_xyy_0[i] * pa_z[i] - ta2_xx_xyy_xyy_1[i] * pc_z[i];

        ta2_xx_xyyz_xyz_0[i] = ta2_xx_xyy_xy_0[i] * fe_0 - ta2_xx_xyy_xy_1[i] * fe_0 + ta2_xx_xyy_xyz_0[i] * pa_z[i] - ta2_xx_xyy_xyz_1[i] * pc_z[i];

        ta2_xx_xyyz_xzz_0[i] = ta2_xx_xz_xzz_0[i] * fe_0 - ta2_xx_xz_xzz_1[i] * fe_0 + ta2_xx_xyz_xzz_0[i] * pa_y[i] - ta2_xx_xyz_xzz_1[i] * pc_y[i];

        ta2_xx_xyyz_yyy_0[i] = ta2_xx_xyy_yyy_0[i] * pa_z[i] - ta2_xx_xyy_yyy_1[i] * pc_z[i];

        ta2_xx_xyyz_yyz_0[i] = 2.0 * ta1_x_yyz_yyz_1[i] + ta2_xx_yyz_yyz_0[i] * pa_x[i] - ta2_xx_yyz_yyz_1[i] * pc_x[i];

        ta2_xx_xyyz_yzz_0[i] = 2.0 * ta1_x_yyz_yzz_1[i] + ta2_xx_yyz_yzz_0[i] * pa_x[i] - ta2_xx_yyz_yzz_1[i] * pc_x[i];

        ta2_xx_xyyz_zzz_0[i] = 2.0 * ta1_x_yyz_zzz_1[i] + ta2_xx_yyz_zzz_0[i] * pa_x[i] - ta2_xx_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto ta2_xx_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 80);

    auto ta2_xx_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 81);

    auto ta2_xx_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 82);

    auto ta2_xx_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 83);

    auto ta2_xx_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 84);

    auto ta2_xx_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 85);

    auto ta2_xx_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 86);

    auto ta2_xx_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 87);

    auto ta2_xx_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 88);

    auto ta2_xx_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 89);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_yzz_yyy_1, ta1_x_yzz_yyz_1, ta1_x_yzz_yzz_1, ta2_xx_xyzz_xxx_0, ta2_xx_xyzz_xxy_0, ta2_xx_xyzz_xxz_0, ta2_xx_xyzz_xyy_0, ta2_xx_xyzz_xyz_0, ta2_xx_xyzz_xzz_0, ta2_xx_xyzz_yyy_0, ta2_xx_xyzz_yyz_0, ta2_xx_xyzz_yzz_0, ta2_xx_xyzz_zzz_0, ta2_xx_xzz_xx_0, ta2_xx_xzz_xx_1, ta2_xx_xzz_xxx_0, ta2_xx_xzz_xxx_1, ta2_xx_xzz_xxy_0, ta2_xx_xzz_xxy_1, ta2_xx_xzz_xxz_0, ta2_xx_xzz_xxz_1, ta2_xx_xzz_xy_0, ta2_xx_xzz_xy_1, ta2_xx_xzz_xyy_0, ta2_xx_xzz_xyy_1, ta2_xx_xzz_xyz_0, ta2_xx_xzz_xyz_1, ta2_xx_xzz_xz_0, ta2_xx_xzz_xz_1, ta2_xx_xzz_xzz_0, ta2_xx_xzz_xzz_1, ta2_xx_xzz_zzz_0, ta2_xx_xzz_zzz_1, ta2_xx_yzz_yyy_0, ta2_xx_yzz_yyy_1, ta2_xx_yzz_yyz_0, ta2_xx_yzz_yyz_1, ta2_xx_yzz_yzz_0, ta2_xx_yzz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyzz_xxx_0[i] = ta2_xx_xzz_xxx_0[i] * pa_y[i] - ta2_xx_xzz_xxx_1[i] * pc_y[i];

        ta2_xx_xyzz_xxy_0[i] = ta2_xx_xzz_xx_0[i] * fe_0 - ta2_xx_xzz_xx_1[i] * fe_0 + ta2_xx_xzz_xxy_0[i] * pa_y[i] - ta2_xx_xzz_xxy_1[i] * pc_y[i];

        ta2_xx_xyzz_xxz_0[i] = ta2_xx_xzz_xxz_0[i] * pa_y[i] - ta2_xx_xzz_xxz_1[i] * pc_y[i];

        ta2_xx_xyzz_xyy_0[i] = 2.0 * ta2_xx_xzz_xy_0[i] * fe_0 - 2.0 * ta2_xx_xzz_xy_1[i] * fe_0 + ta2_xx_xzz_xyy_0[i] * pa_y[i] - ta2_xx_xzz_xyy_1[i] * pc_y[i];

        ta2_xx_xyzz_xyz_0[i] = ta2_xx_xzz_xz_0[i] * fe_0 - ta2_xx_xzz_xz_1[i] * fe_0 + ta2_xx_xzz_xyz_0[i] * pa_y[i] - ta2_xx_xzz_xyz_1[i] * pc_y[i];

        ta2_xx_xyzz_xzz_0[i] = ta2_xx_xzz_xzz_0[i] * pa_y[i] - ta2_xx_xzz_xzz_1[i] * pc_y[i];

        ta2_xx_xyzz_yyy_0[i] = 2.0 * ta1_x_yzz_yyy_1[i] + ta2_xx_yzz_yyy_0[i] * pa_x[i] - ta2_xx_yzz_yyy_1[i] * pc_x[i];

        ta2_xx_xyzz_yyz_0[i] = 2.0 * ta1_x_yzz_yyz_1[i] + ta2_xx_yzz_yyz_0[i] * pa_x[i] - ta2_xx_yzz_yyz_1[i] * pc_x[i];

        ta2_xx_xyzz_yzz_0[i] = 2.0 * ta1_x_yzz_yzz_1[i] + ta2_xx_yzz_yzz_0[i] * pa_x[i] - ta2_xx_yzz_yzz_1[i] * pc_x[i];

        ta2_xx_xyzz_zzz_0[i] = ta2_xx_xzz_zzz_0[i] * pa_y[i] - ta2_xx_xzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto ta2_xx_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 90);

    auto ta2_xx_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 91);

    auto ta2_xx_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 92);

    auto ta2_xx_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 93);

    auto ta2_xx_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 94);

    auto ta2_xx_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 95);

    auto ta2_xx_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 96);

    auto ta2_xx_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 97);

    auto ta2_xx_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 98);

    auto ta2_xx_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 99);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_zzz_xxz_1, ta1_x_zzz_xyz_1, ta1_x_zzz_xzz_1, ta1_x_zzz_yyy_1, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_1, ta2_xx_xz_xxx_0, ta2_xx_xz_xxx_1, ta2_xx_xz_xxy_0, ta2_xx_xz_xxy_1, ta2_xx_xz_xyy_0, ta2_xx_xz_xyy_1, ta2_xx_xzz_xxx_0, ta2_xx_xzz_xxx_1, ta2_xx_xzz_xxy_0, ta2_xx_xzz_xxy_1, ta2_xx_xzz_xyy_0, ta2_xx_xzz_xyy_1, ta2_xx_xzzz_xxx_0, ta2_xx_xzzz_xxy_0, ta2_xx_xzzz_xxz_0, ta2_xx_xzzz_xyy_0, ta2_xx_xzzz_xyz_0, ta2_xx_xzzz_xzz_0, ta2_xx_xzzz_yyy_0, ta2_xx_xzzz_yyz_0, ta2_xx_xzzz_yzz_0, ta2_xx_xzzz_zzz_0, ta2_xx_zzz_xxz_0, ta2_xx_zzz_xxz_1, ta2_xx_zzz_xyz_0, ta2_xx_zzz_xyz_1, ta2_xx_zzz_xz_0, ta2_xx_zzz_xz_1, ta2_xx_zzz_xzz_0, ta2_xx_zzz_xzz_1, ta2_xx_zzz_yyy_0, ta2_xx_zzz_yyy_1, ta2_xx_zzz_yyz_0, ta2_xx_zzz_yyz_1, ta2_xx_zzz_yz_0, ta2_xx_zzz_yz_1, ta2_xx_zzz_yzz_0, ta2_xx_zzz_yzz_1, ta2_xx_zzz_zz_0, ta2_xx_zzz_zz_1, ta2_xx_zzz_zzz_0, ta2_xx_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzzz_xxx_0[i] = 2.0 * ta2_xx_xz_xxx_0[i] * fe_0 - 2.0 * ta2_xx_xz_xxx_1[i] * fe_0 + ta2_xx_xzz_xxx_0[i] * pa_z[i] - ta2_xx_xzz_xxx_1[i] * pc_z[i];

        ta2_xx_xzzz_xxy_0[i] = 2.0 * ta2_xx_xz_xxy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xxy_1[i] * fe_0 + ta2_xx_xzz_xxy_0[i] * pa_z[i] - ta2_xx_xzz_xxy_1[i] * pc_z[i];

        ta2_xx_xzzz_xxz_0[i] = 2.0 * ta2_xx_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxz_1[i] + ta2_xx_zzz_xxz_0[i] * pa_x[i] - ta2_xx_zzz_xxz_1[i] * pc_x[i];

        ta2_xx_xzzz_xyy_0[i] = 2.0 * ta2_xx_xz_xyy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xyy_1[i] * fe_0 + ta2_xx_xzz_xyy_0[i] * pa_z[i] - ta2_xx_xzz_xyy_1[i] * pc_z[i];

        ta2_xx_xzzz_xyz_0[i] = ta2_xx_zzz_yz_0[i] * fe_0 - ta2_xx_zzz_yz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyz_1[i] + ta2_xx_zzz_xyz_0[i] * pa_x[i] - ta2_xx_zzz_xyz_1[i] * pc_x[i];

        ta2_xx_xzzz_xzz_0[i] = ta2_xx_zzz_zz_0[i] * fe_0 - ta2_xx_zzz_zz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xzz_1[i] + ta2_xx_zzz_xzz_0[i] * pa_x[i] - ta2_xx_zzz_xzz_1[i] * pc_x[i];

        ta2_xx_xzzz_yyy_0[i] = 2.0 * ta1_x_zzz_yyy_1[i] + ta2_xx_zzz_yyy_0[i] * pa_x[i] - ta2_xx_zzz_yyy_1[i] * pc_x[i];

        ta2_xx_xzzz_yyz_0[i] = 2.0 * ta1_x_zzz_yyz_1[i] + ta2_xx_zzz_yyz_0[i] * pa_x[i] - ta2_xx_zzz_yyz_1[i] * pc_x[i];

        ta2_xx_xzzz_yzz_0[i] = 2.0 * ta1_x_zzz_yzz_1[i] + ta2_xx_zzz_yzz_0[i] * pa_x[i] - ta2_xx_zzz_yzz_1[i] * pc_x[i];

        ta2_xx_xzzz_zzz_0[i] = 2.0 * ta1_x_zzz_zzz_1[i] + ta2_xx_zzz_zzz_0[i] * pa_x[i] - ta2_xx_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto ta2_xx_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 100);

    auto ta2_xx_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 101);

    auto ta2_xx_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 102);

    auto ta2_xx_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 103);

    auto ta2_xx_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 104);

    auto ta2_xx_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 105);

    auto ta2_xx_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 106);

    auto ta2_xx_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 107);

    auto ta2_xx_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 108);

    auto ta2_xx_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 109);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_yy_xxx_0, ta2_xx_yy_xxx_1, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xxz_0, ta2_xx_yy_xxz_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_xyz_0, ta2_xx_yy_xyz_1, ta2_xx_yy_xzz_0, ta2_xx_yy_xzz_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yy_yyz_0, ta2_xx_yy_yyz_1, ta2_xx_yy_yzz_0, ta2_xx_yy_yzz_1, ta2_xx_yy_zzz_0, ta2_xx_yy_zzz_1, ta2_xx_yyy_xx_0, ta2_xx_yyy_xx_1, ta2_xx_yyy_xxx_0, ta2_xx_yyy_xxx_1, ta2_xx_yyy_xxy_0, ta2_xx_yyy_xxy_1, ta2_xx_yyy_xxz_0, ta2_xx_yyy_xxz_1, ta2_xx_yyy_xy_0, ta2_xx_yyy_xy_1, ta2_xx_yyy_xyy_0, ta2_xx_yyy_xyy_1, ta2_xx_yyy_xyz_0, ta2_xx_yyy_xyz_1, ta2_xx_yyy_xz_0, ta2_xx_yyy_xz_1, ta2_xx_yyy_xzz_0, ta2_xx_yyy_xzz_1, ta2_xx_yyy_yy_0, ta2_xx_yyy_yy_1, ta2_xx_yyy_yyy_0, ta2_xx_yyy_yyy_1, ta2_xx_yyy_yyz_0, ta2_xx_yyy_yyz_1, ta2_xx_yyy_yz_0, ta2_xx_yyy_yz_1, ta2_xx_yyy_yzz_0, ta2_xx_yyy_yzz_1, ta2_xx_yyy_zz_0, ta2_xx_yyy_zz_1, ta2_xx_yyy_zzz_0, ta2_xx_yyy_zzz_1, ta2_xx_yyyy_xxx_0, ta2_xx_yyyy_xxy_0, ta2_xx_yyyy_xxz_0, ta2_xx_yyyy_xyy_0, ta2_xx_yyyy_xyz_0, ta2_xx_yyyy_xzz_0, ta2_xx_yyyy_yyy_0, ta2_xx_yyyy_yyz_0, ta2_xx_yyyy_yzz_0, ta2_xx_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyy_xxx_0[i] = 3.0 * ta2_xx_yy_xxx_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxx_1[i] * fe_0 + ta2_xx_yyy_xxx_0[i] * pa_y[i] - ta2_xx_yyy_xxx_1[i] * pc_y[i];

        ta2_xx_yyyy_xxy_0[i] = 3.0 * ta2_xx_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxy_1[i] * fe_0 + ta2_xx_yyy_xx_0[i] * fe_0 - ta2_xx_yyy_xx_1[i] * fe_0 + ta2_xx_yyy_xxy_0[i] * pa_y[i] - ta2_xx_yyy_xxy_1[i] * pc_y[i];

        ta2_xx_yyyy_xxz_0[i] = 3.0 * ta2_xx_yy_xxz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxz_1[i] * fe_0 + ta2_xx_yyy_xxz_0[i] * pa_y[i] - ta2_xx_yyy_xxz_1[i] * pc_y[i];

        ta2_xx_yyyy_xyy_0[i] = 3.0 * ta2_xx_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyy_1[i] * fe_0 + 2.0 * ta2_xx_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xx_yyy_xy_1[i] * fe_0 + ta2_xx_yyy_xyy_0[i] * pa_y[i] - ta2_xx_yyy_xyy_1[i] * pc_y[i];

        ta2_xx_yyyy_xyz_0[i] = 3.0 * ta2_xx_yy_xyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyz_1[i] * fe_0 + ta2_xx_yyy_xz_0[i] * fe_0 - ta2_xx_yyy_xz_1[i] * fe_0 + ta2_xx_yyy_xyz_0[i] * pa_y[i] - ta2_xx_yyy_xyz_1[i] * pc_y[i];

        ta2_xx_yyyy_xzz_0[i] = 3.0 * ta2_xx_yy_xzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xzz_1[i] * fe_0 + ta2_xx_yyy_xzz_0[i] * pa_y[i] - ta2_xx_yyy_xzz_1[i] * pc_y[i];

        ta2_xx_yyyy_yyy_0[i] = 3.0 * ta2_xx_yy_yyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyy_1[i] * fe_0 + 3.0 * ta2_xx_yyy_yy_0[i] * fe_0 - 3.0 * ta2_xx_yyy_yy_1[i] * fe_0 + ta2_xx_yyy_yyy_0[i] * pa_y[i] - ta2_xx_yyy_yyy_1[i] * pc_y[i];

        ta2_xx_yyyy_yyz_0[i] = 3.0 * ta2_xx_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyz_1[i] * fe_0 + 2.0 * ta2_xx_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xx_yyy_yz_1[i] * fe_0 + ta2_xx_yyy_yyz_0[i] * pa_y[i] - ta2_xx_yyy_yyz_1[i] * pc_y[i];

        ta2_xx_yyyy_yzz_0[i] = 3.0 * ta2_xx_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yzz_1[i] * fe_0 + ta2_xx_yyy_zz_0[i] * fe_0 - ta2_xx_yyy_zz_1[i] * fe_0 + ta2_xx_yyy_yzz_0[i] * pa_y[i] - ta2_xx_yyy_yzz_1[i] * pc_y[i];

        ta2_xx_yyyy_zzz_0[i] = 3.0 * ta2_xx_yy_zzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_zzz_1[i] * fe_0 + ta2_xx_yyy_zzz_0[i] * pa_y[i] - ta2_xx_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto ta2_xx_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 110);

    auto ta2_xx_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 111);

    auto ta2_xx_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 112);

    auto ta2_xx_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 113);

    auto ta2_xx_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 114);

    auto ta2_xx_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 115);

    auto ta2_xx_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 116);

    auto ta2_xx_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 117);

    auto ta2_xx_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 118);

    auto ta2_xx_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 119);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_yyy_xxx_0, ta2_xx_yyy_xxx_1, ta2_xx_yyy_xxy_0, ta2_xx_yyy_xxy_1, ta2_xx_yyy_xy_0, ta2_xx_yyy_xy_1, ta2_xx_yyy_xyy_0, ta2_xx_yyy_xyy_1, ta2_xx_yyy_xyz_0, ta2_xx_yyy_xyz_1, ta2_xx_yyy_yy_0, ta2_xx_yyy_yy_1, ta2_xx_yyy_yyy_0, ta2_xx_yyy_yyy_1, ta2_xx_yyy_yyz_0, ta2_xx_yyy_yyz_1, ta2_xx_yyy_yz_0, ta2_xx_yyy_yz_1, ta2_xx_yyy_yzz_0, ta2_xx_yyy_yzz_1, ta2_xx_yyyz_xxx_0, ta2_xx_yyyz_xxy_0, ta2_xx_yyyz_xxz_0, ta2_xx_yyyz_xyy_0, ta2_xx_yyyz_xyz_0, ta2_xx_yyyz_xzz_0, ta2_xx_yyyz_yyy_0, ta2_xx_yyyz_yyz_0, ta2_xx_yyyz_yzz_0, ta2_xx_yyyz_zzz_0, ta2_xx_yyz_xxz_0, ta2_xx_yyz_xxz_1, ta2_xx_yyz_xzz_0, ta2_xx_yyz_xzz_1, ta2_xx_yyz_zzz_0, ta2_xx_yyz_zzz_1, ta2_xx_yz_xxz_0, ta2_xx_yz_xxz_1, ta2_xx_yz_xzz_0, ta2_xx_yz_xzz_1, ta2_xx_yz_zzz_0, ta2_xx_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyz_xxx_0[i] = ta2_xx_yyy_xxx_0[i] * pa_z[i] - ta2_xx_yyy_xxx_1[i] * pc_z[i];

        ta2_xx_yyyz_xxy_0[i] = ta2_xx_yyy_xxy_0[i] * pa_z[i] - ta2_xx_yyy_xxy_1[i] * pc_z[i];

        ta2_xx_yyyz_xxz_0[i] = 2.0 * ta2_xx_yz_xxz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xxz_1[i] * fe_0 + ta2_xx_yyz_xxz_0[i] * pa_y[i] - ta2_xx_yyz_xxz_1[i] * pc_y[i];

        ta2_xx_yyyz_xyy_0[i] = ta2_xx_yyy_xyy_0[i] * pa_z[i] - ta2_xx_yyy_xyy_1[i] * pc_z[i];

        ta2_xx_yyyz_xyz_0[i] = ta2_xx_yyy_xy_0[i] * fe_0 - ta2_xx_yyy_xy_1[i] * fe_0 + ta2_xx_yyy_xyz_0[i] * pa_z[i] - ta2_xx_yyy_xyz_1[i] * pc_z[i];

        ta2_xx_yyyz_xzz_0[i] = 2.0 * ta2_xx_yz_xzz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xzz_1[i] * fe_0 + ta2_xx_yyz_xzz_0[i] * pa_y[i] - ta2_xx_yyz_xzz_1[i] * pc_y[i];

        ta2_xx_yyyz_yyy_0[i] = ta2_xx_yyy_yyy_0[i] * pa_z[i] - ta2_xx_yyy_yyy_1[i] * pc_z[i];

        ta2_xx_yyyz_yyz_0[i] = ta2_xx_yyy_yy_0[i] * fe_0 - ta2_xx_yyy_yy_1[i] * fe_0 + ta2_xx_yyy_yyz_0[i] * pa_z[i] - ta2_xx_yyy_yyz_1[i] * pc_z[i];

        ta2_xx_yyyz_yzz_0[i] = 2.0 * ta2_xx_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xx_yyy_yz_1[i] * fe_0 + ta2_xx_yyy_yzz_0[i] * pa_z[i] - ta2_xx_yyy_yzz_1[i] * pc_z[i];

        ta2_xx_yyyz_zzz_0[i] = 2.0 * ta2_xx_yz_zzz_0[i] * fe_0 - 2.0 * ta2_xx_yz_zzz_1[i] * fe_0 + ta2_xx_yyz_zzz_0[i] * pa_y[i] - ta2_xx_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto ta2_xx_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 120);

    auto ta2_xx_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 121);

    auto ta2_xx_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 122);

    auto ta2_xx_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 123);

    auto ta2_xx_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 124);

    auto ta2_xx_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 125);

    auto ta2_xx_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 126);

    auto ta2_xx_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 127);

    auto ta2_xx_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 128);

    auto ta2_xx_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 129);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_yy_xxy_0, ta2_xx_yy_xxy_1, ta2_xx_yy_xyy_0, ta2_xx_yy_xyy_1, ta2_xx_yy_yyy_0, ta2_xx_yy_yyy_1, ta2_xx_yyz_xxy_0, ta2_xx_yyz_xxy_1, ta2_xx_yyz_xyy_0, ta2_xx_yyz_xyy_1, ta2_xx_yyz_yyy_0, ta2_xx_yyz_yyy_1, ta2_xx_yyzz_xxx_0, ta2_xx_yyzz_xxy_0, ta2_xx_yyzz_xxz_0, ta2_xx_yyzz_xyy_0, ta2_xx_yyzz_xyz_0, ta2_xx_yyzz_xzz_0, ta2_xx_yyzz_yyy_0, ta2_xx_yyzz_yyz_0, ta2_xx_yyzz_yzz_0, ta2_xx_yyzz_zzz_0, ta2_xx_yzz_xxx_0, ta2_xx_yzz_xxx_1, ta2_xx_yzz_xxz_0, ta2_xx_yzz_xxz_1, ta2_xx_yzz_xyz_0, ta2_xx_yzz_xyz_1, ta2_xx_yzz_xz_0, ta2_xx_yzz_xz_1, ta2_xx_yzz_xzz_0, ta2_xx_yzz_xzz_1, ta2_xx_yzz_yyz_0, ta2_xx_yzz_yyz_1, ta2_xx_yzz_yz_0, ta2_xx_yzz_yz_1, ta2_xx_yzz_yzz_0, ta2_xx_yzz_yzz_1, ta2_xx_yzz_zz_0, ta2_xx_yzz_zz_1, ta2_xx_yzz_zzz_0, ta2_xx_yzz_zzz_1, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyzz_xxx_0[i] = ta2_xx_zz_xxx_0[i] * fe_0 - ta2_xx_zz_xxx_1[i] * fe_0 + ta2_xx_yzz_xxx_0[i] * pa_y[i] - ta2_xx_yzz_xxx_1[i] * pc_y[i];

        ta2_xx_yyzz_xxy_0[i] = ta2_xx_yy_xxy_0[i] * fe_0 - ta2_xx_yy_xxy_1[i] * fe_0 + ta2_xx_yyz_xxy_0[i] * pa_z[i] - ta2_xx_yyz_xxy_1[i] * pc_z[i];

        ta2_xx_yyzz_xxz_0[i] = ta2_xx_zz_xxz_0[i] * fe_0 - ta2_xx_zz_xxz_1[i] * fe_0 + ta2_xx_yzz_xxz_0[i] * pa_y[i] - ta2_xx_yzz_xxz_1[i] * pc_y[i];

        ta2_xx_yyzz_xyy_0[i] = ta2_xx_yy_xyy_0[i] * fe_0 - ta2_xx_yy_xyy_1[i] * fe_0 + ta2_xx_yyz_xyy_0[i] * pa_z[i] - ta2_xx_yyz_xyy_1[i] * pc_z[i];

        ta2_xx_yyzz_xyz_0[i] = ta2_xx_zz_xyz_0[i] * fe_0 - ta2_xx_zz_xyz_1[i] * fe_0 + ta2_xx_yzz_xz_0[i] * fe_0 - ta2_xx_yzz_xz_1[i] * fe_0 + ta2_xx_yzz_xyz_0[i] * pa_y[i] - ta2_xx_yzz_xyz_1[i] * pc_y[i];

        ta2_xx_yyzz_xzz_0[i] = ta2_xx_zz_xzz_0[i] * fe_0 - ta2_xx_zz_xzz_1[i] * fe_0 + ta2_xx_yzz_xzz_0[i] * pa_y[i] - ta2_xx_yzz_xzz_1[i] * pc_y[i];

        ta2_xx_yyzz_yyy_0[i] = ta2_xx_yy_yyy_0[i] * fe_0 - ta2_xx_yy_yyy_1[i] * fe_0 + ta2_xx_yyz_yyy_0[i] * pa_z[i] - ta2_xx_yyz_yyy_1[i] * pc_z[i];

        ta2_xx_yyzz_yyz_0[i] = ta2_xx_zz_yyz_0[i] * fe_0 - ta2_xx_zz_yyz_1[i] * fe_0 + 2.0 * ta2_xx_yzz_yz_0[i] * fe_0 - 2.0 * ta2_xx_yzz_yz_1[i] * fe_0 + ta2_xx_yzz_yyz_0[i] * pa_y[i] - ta2_xx_yzz_yyz_1[i] * pc_y[i];

        ta2_xx_yyzz_yzz_0[i] = ta2_xx_zz_yzz_0[i] * fe_0 - ta2_xx_zz_yzz_1[i] * fe_0 + ta2_xx_yzz_zz_0[i] * fe_0 - ta2_xx_yzz_zz_1[i] * fe_0 + ta2_xx_yzz_yzz_0[i] * pa_y[i] - ta2_xx_yzz_yzz_1[i] * pc_y[i];

        ta2_xx_yyzz_zzz_0[i] = ta2_xx_zz_zzz_0[i] * fe_0 - ta2_xx_zz_zzz_1[i] * fe_0 + ta2_xx_yzz_zzz_0[i] * pa_y[i] - ta2_xx_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto ta2_xx_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 130);

    auto ta2_xx_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 131);

    auto ta2_xx_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 132);

    auto ta2_xx_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 133);

    auto ta2_xx_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 134);

    auto ta2_xx_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 135);

    auto ta2_xx_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 136);

    auto ta2_xx_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 137);

    auto ta2_xx_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 138);

    auto ta2_xx_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 139);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_yzzz_xxx_0, ta2_xx_yzzz_xxy_0, ta2_xx_yzzz_xxz_0, ta2_xx_yzzz_xyy_0, ta2_xx_yzzz_xyz_0, ta2_xx_yzzz_xzz_0, ta2_xx_yzzz_yyy_0, ta2_xx_yzzz_yyz_0, ta2_xx_yzzz_yzz_0, ta2_xx_yzzz_zzz_0, ta2_xx_zzz_xx_0, ta2_xx_zzz_xx_1, ta2_xx_zzz_xxx_0, ta2_xx_zzz_xxx_1, ta2_xx_zzz_xxy_0, ta2_xx_zzz_xxy_1, ta2_xx_zzz_xxz_0, ta2_xx_zzz_xxz_1, ta2_xx_zzz_xy_0, ta2_xx_zzz_xy_1, ta2_xx_zzz_xyy_0, ta2_xx_zzz_xyy_1, ta2_xx_zzz_xyz_0, ta2_xx_zzz_xyz_1, ta2_xx_zzz_xz_0, ta2_xx_zzz_xz_1, ta2_xx_zzz_xzz_0, ta2_xx_zzz_xzz_1, ta2_xx_zzz_yy_0, ta2_xx_zzz_yy_1, ta2_xx_zzz_yyy_0, ta2_xx_zzz_yyy_1, ta2_xx_zzz_yyz_0, ta2_xx_zzz_yyz_1, ta2_xx_zzz_yz_0, ta2_xx_zzz_yz_1, ta2_xx_zzz_yzz_0, ta2_xx_zzz_yzz_1, ta2_xx_zzz_zz_0, ta2_xx_zzz_zz_1, ta2_xx_zzz_zzz_0, ta2_xx_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzzz_xxx_0[i] = ta2_xx_zzz_xxx_0[i] * pa_y[i] - ta2_xx_zzz_xxx_1[i] * pc_y[i];

        ta2_xx_yzzz_xxy_0[i] = ta2_xx_zzz_xx_0[i] * fe_0 - ta2_xx_zzz_xx_1[i] * fe_0 + ta2_xx_zzz_xxy_0[i] * pa_y[i] - ta2_xx_zzz_xxy_1[i] * pc_y[i];

        ta2_xx_yzzz_xxz_0[i] = ta2_xx_zzz_xxz_0[i] * pa_y[i] - ta2_xx_zzz_xxz_1[i] * pc_y[i];

        ta2_xx_yzzz_xyy_0[i] = 2.0 * ta2_xx_zzz_xy_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xy_1[i] * fe_0 + ta2_xx_zzz_xyy_0[i] * pa_y[i] - ta2_xx_zzz_xyy_1[i] * pc_y[i];

        ta2_xx_yzzz_xyz_0[i] = ta2_xx_zzz_xz_0[i] * fe_0 - ta2_xx_zzz_xz_1[i] * fe_0 + ta2_xx_zzz_xyz_0[i] * pa_y[i] - ta2_xx_zzz_xyz_1[i] * pc_y[i];

        ta2_xx_yzzz_xzz_0[i] = ta2_xx_zzz_xzz_0[i] * pa_y[i] - ta2_xx_zzz_xzz_1[i] * pc_y[i];

        ta2_xx_yzzz_yyy_0[i] = 3.0 * ta2_xx_zzz_yy_0[i] * fe_0 - 3.0 * ta2_xx_zzz_yy_1[i] * fe_0 + ta2_xx_zzz_yyy_0[i] * pa_y[i] - ta2_xx_zzz_yyy_1[i] * pc_y[i];

        ta2_xx_yzzz_yyz_0[i] = 2.0 * ta2_xx_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_yz_1[i] * fe_0 + ta2_xx_zzz_yyz_0[i] * pa_y[i] - ta2_xx_zzz_yyz_1[i] * pc_y[i];

        ta2_xx_yzzz_yzz_0[i] = ta2_xx_zzz_zz_0[i] * fe_0 - ta2_xx_zzz_zz_1[i] * fe_0 + ta2_xx_zzz_yzz_0[i] * pa_y[i] - ta2_xx_zzz_yzz_1[i] * pc_y[i];

        ta2_xx_yzzz_zzz_0[i] = ta2_xx_zzz_zzz_0[i] * pa_y[i] - ta2_xx_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto ta2_xx_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 140);

    auto ta2_xx_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 141);

    auto ta2_xx_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 142);

    auto ta2_xx_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 143);

    auto ta2_xx_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 144);

    auto ta2_xx_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 145);

    auto ta2_xx_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 146);

    auto ta2_xx_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 147);

    auto ta2_xx_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 148);

    auto ta2_xx_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 149);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_zz_xxx_0, ta2_xx_zz_xxx_1, ta2_xx_zz_xxy_0, ta2_xx_zz_xxy_1, ta2_xx_zz_xxz_0, ta2_xx_zz_xxz_1, ta2_xx_zz_xyy_0, ta2_xx_zz_xyy_1, ta2_xx_zz_xyz_0, ta2_xx_zz_xyz_1, ta2_xx_zz_xzz_0, ta2_xx_zz_xzz_1, ta2_xx_zz_yyy_0, ta2_xx_zz_yyy_1, ta2_xx_zz_yyz_0, ta2_xx_zz_yyz_1, ta2_xx_zz_yzz_0, ta2_xx_zz_yzz_1, ta2_xx_zz_zzz_0, ta2_xx_zz_zzz_1, ta2_xx_zzz_xx_0, ta2_xx_zzz_xx_1, ta2_xx_zzz_xxx_0, ta2_xx_zzz_xxx_1, ta2_xx_zzz_xxy_0, ta2_xx_zzz_xxy_1, ta2_xx_zzz_xxz_0, ta2_xx_zzz_xxz_1, ta2_xx_zzz_xy_0, ta2_xx_zzz_xy_1, ta2_xx_zzz_xyy_0, ta2_xx_zzz_xyy_1, ta2_xx_zzz_xyz_0, ta2_xx_zzz_xyz_1, ta2_xx_zzz_xz_0, ta2_xx_zzz_xz_1, ta2_xx_zzz_xzz_0, ta2_xx_zzz_xzz_1, ta2_xx_zzz_yy_0, ta2_xx_zzz_yy_1, ta2_xx_zzz_yyy_0, ta2_xx_zzz_yyy_1, ta2_xx_zzz_yyz_0, ta2_xx_zzz_yyz_1, ta2_xx_zzz_yz_0, ta2_xx_zzz_yz_1, ta2_xx_zzz_yzz_0, ta2_xx_zzz_yzz_1, ta2_xx_zzz_zz_0, ta2_xx_zzz_zz_1, ta2_xx_zzz_zzz_0, ta2_xx_zzz_zzz_1, ta2_xx_zzzz_xxx_0, ta2_xx_zzzz_xxy_0, ta2_xx_zzzz_xxz_0, ta2_xx_zzzz_xyy_0, ta2_xx_zzzz_xyz_0, ta2_xx_zzzz_xzz_0, ta2_xx_zzzz_yyy_0, ta2_xx_zzzz_yyz_0, ta2_xx_zzzz_yzz_0, ta2_xx_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzzz_xxx_0[i] = 3.0 * ta2_xx_zz_xxx_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxx_1[i] * fe_0 + ta2_xx_zzz_xxx_0[i] * pa_z[i] - ta2_xx_zzz_xxx_1[i] * pc_z[i];

        ta2_xx_zzzz_xxy_0[i] = 3.0 * ta2_xx_zz_xxy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxy_1[i] * fe_0 + ta2_xx_zzz_xxy_0[i] * pa_z[i] - ta2_xx_zzz_xxy_1[i] * pc_z[i];

        ta2_xx_zzzz_xxz_0[i] = 3.0 * ta2_xx_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxz_1[i] * fe_0 + ta2_xx_zzz_xx_0[i] * fe_0 - ta2_xx_zzz_xx_1[i] * fe_0 + ta2_xx_zzz_xxz_0[i] * pa_z[i] - ta2_xx_zzz_xxz_1[i] * pc_z[i];

        ta2_xx_zzzz_xyy_0[i] = 3.0 * ta2_xx_zz_xyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyy_1[i] * fe_0 + ta2_xx_zzz_xyy_0[i] * pa_z[i] - ta2_xx_zzz_xyy_1[i] * pc_z[i];

        ta2_xx_zzzz_xyz_0[i] = 3.0 * ta2_xx_zz_xyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyz_1[i] * fe_0 + ta2_xx_zzz_xy_0[i] * fe_0 - ta2_xx_zzz_xy_1[i] * fe_0 + ta2_xx_zzz_xyz_0[i] * pa_z[i] - ta2_xx_zzz_xyz_1[i] * pc_z[i];

        ta2_xx_zzzz_xzz_0[i] = 3.0 * ta2_xx_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xzz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xz_1[i] * fe_0 + ta2_xx_zzz_xzz_0[i] * pa_z[i] - ta2_xx_zzz_xzz_1[i] * pc_z[i];

        ta2_xx_zzzz_yyy_0[i] = 3.0 * ta2_xx_zz_yyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyy_1[i] * fe_0 + ta2_xx_zzz_yyy_0[i] * pa_z[i] - ta2_xx_zzz_yyy_1[i] * pc_z[i];

        ta2_xx_zzzz_yyz_0[i] = 3.0 * ta2_xx_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyz_1[i] * fe_0 + ta2_xx_zzz_yy_0[i] * fe_0 - ta2_xx_zzz_yy_1[i] * fe_0 + ta2_xx_zzz_yyz_0[i] * pa_z[i] - ta2_xx_zzz_yyz_1[i] * pc_z[i];

        ta2_xx_zzzz_yzz_0[i] = 3.0 * ta2_xx_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yzz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_yz_1[i] * fe_0 + ta2_xx_zzz_yzz_0[i] * pa_z[i] - ta2_xx_zzz_yzz_1[i] * pc_z[i];

        ta2_xx_zzzz_zzz_0[i] = 3.0 * ta2_xx_zz_zzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_zzz_1[i] * fe_0 + 3.0 * ta2_xx_zzz_zz_0[i] * fe_0 - 3.0 * ta2_xx_zzz_zz_1[i] * fe_0 + ta2_xx_zzz_zzz_0[i] * pa_z[i] - ta2_xx_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 150-160 components of targeted buffer : GF

    auto ta2_xy_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 150);

    auto ta2_xy_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 151);

    auto ta2_xy_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 152);

    auto ta2_xy_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 153);

    auto ta2_xy_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 154);

    auto ta2_xy_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 155);

    auto ta2_xy_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 156);

    auto ta2_xy_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 157);

    auto ta2_xy_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 158);

    auto ta2_xy_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 159);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xxx_xxx_1, ta1_y_xxx_xxy_1, ta1_y_xxx_xxz_1, ta1_y_xxx_xyy_1, ta1_y_xxx_xyz_1, ta1_y_xxx_xzz_1, ta1_y_xxx_yyy_1, ta1_y_xxx_yyz_1, ta1_y_xxx_yzz_1, ta1_y_xxx_zzz_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xx_yyz_0, ta2_xy_xx_yyz_1, ta2_xy_xx_yzz_0, ta2_xy_xx_yzz_1, ta2_xy_xx_zzz_0, ta2_xy_xx_zzz_1, ta2_xy_xxx_xx_0, ta2_xy_xxx_xx_1, ta2_xy_xxx_xxx_0, ta2_xy_xxx_xxx_1, ta2_xy_xxx_xxy_0, ta2_xy_xxx_xxy_1, ta2_xy_xxx_xxz_0, ta2_xy_xxx_xxz_1, ta2_xy_xxx_xy_0, ta2_xy_xxx_xy_1, ta2_xy_xxx_xyy_0, ta2_xy_xxx_xyy_1, ta2_xy_xxx_xyz_0, ta2_xy_xxx_xyz_1, ta2_xy_xxx_xz_0, ta2_xy_xxx_xz_1, ta2_xy_xxx_xzz_0, ta2_xy_xxx_xzz_1, ta2_xy_xxx_yy_0, ta2_xy_xxx_yy_1, ta2_xy_xxx_yyy_0, ta2_xy_xxx_yyy_1, ta2_xy_xxx_yyz_0, ta2_xy_xxx_yyz_1, ta2_xy_xxx_yz_0, ta2_xy_xxx_yz_1, ta2_xy_xxx_yzz_0, ta2_xy_xxx_yzz_1, ta2_xy_xxx_zz_0, ta2_xy_xxx_zz_1, ta2_xy_xxx_zzz_0, ta2_xy_xxx_zzz_1, ta2_xy_xxxx_xxx_0, ta2_xy_xxxx_xxy_0, ta2_xy_xxxx_xxz_0, ta2_xy_xxxx_xyy_0, ta2_xy_xxxx_xyz_0, ta2_xy_xxxx_xzz_0, ta2_xy_xxxx_yyy_0, ta2_xy_xxxx_yyz_0, ta2_xy_xxxx_yzz_0, ta2_xy_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxx_xxx_0[i] = 3.0 * ta2_xy_xx_xxx_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxx_1[i] * fe_0 + 3.0 * ta2_xy_xxx_xx_0[i] * fe_0 - 3.0 * ta2_xy_xxx_xx_1[i] * fe_0 + ta1_y_xxx_xxx_1[i] + ta2_xy_xxx_xxx_0[i] * pa_x[i] - ta2_xy_xxx_xxx_1[i] * pc_x[i];

        ta2_xy_xxxx_xxy_0[i] = 3.0 * ta2_xy_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxy_1[i] * fe_0 + 2.0 * ta2_xy_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xy_1[i] * fe_0 + ta1_y_xxx_xxy_1[i] + ta2_xy_xxx_xxy_0[i] * pa_x[i] - ta2_xy_xxx_xxy_1[i] * pc_x[i];

        ta2_xy_xxxx_xxz_0[i] = 3.0 * ta2_xy_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxz_1[i] * fe_0 + 2.0 * ta2_xy_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xz_1[i] * fe_0 + ta1_y_xxx_xxz_1[i] + ta2_xy_xxx_xxz_0[i] * pa_x[i] - ta2_xy_xxx_xxz_1[i] * pc_x[i];

        ta2_xy_xxxx_xyy_0[i] = 3.0 * ta2_xy_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyy_1[i] * fe_0 + ta2_xy_xxx_yy_0[i] * fe_0 - ta2_xy_xxx_yy_1[i] * fe_0 + ta1_y_xxx_xyy_1[i] + ta2_xy_xxx_xyy_0[i] * pa_x[i] - ta2_xy_xxx_xyy_1[i] * pc_x[i];

        ta2_xy_xxxx_xyz_0[i] = 3.0 * ta2_xy_xx_xyz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyz_1[i] * fe_0 + ta2_xy_xxx_yz_0[i] * fe_0 - ta2_xy_xxx_yz_1[i] * fe_0 + ta1_y_xxx_xyz_1[i] + ta2_xy_xxx_xyz_0[i] * pa_x[i] - ta2_xy_xxx_xyz_1[i] * pc_x[i];

        ta2_xy_xxxx_xzz_0[i] = 3.0 * ta2_xy_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xzz_1[i] * fe_0 + ta2_xy_xxx_zz_0[i] * fe_0 - ta2_xy_xxx_zz_1[i] * fe_0 + ta1_y_xxx_xzz_1[i] + ta2_xy_xxx_xzz_0[i] * pa_x[i] - ta2_xy_xxx_xzz_1[i] * pc_x[i];

        ta2_xy_xxxx_yyy_0[i] = 3.0 * ta2_xy_xx_yyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_yyy_1[i] * fe_0 + ta1_y_xxx_yyy_1[i] + ta2_xy_xxx_yyy_0[i] * pa_x[i] - ta2_xy_xxx_yyy_1[i] * pc_x[i];

        ta2_xy_xxxx_yyz_0[i] = 3.0 * ta2_xy_xx_yyz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yyz_1[i] * fe_0 + ta1_y_xxx_yyz_1[i] + ta2_xy_xxx_yyz_0[i] * pa_x[i] - ta2_xy_xxx_yyz_1[i] * pc_x[i];

        ta2_xy_xxxx_yzz_0[i] = 3.0 * ta2_xy_xx_yzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yzz_1[i] * fe_0 + ta1_y_xxx_yzz_1[i] + ta2_xy_xxx_yzz_0[i] * pa_x[i] - ta2_xy_xxx_yzz_1[i] * pc_x[i];

        ta2_xy_xxxx_zzz_0[i] = 3.0 * ta2_xy_xx_zzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_zzz_1[i] * fe_0 + ta1_y_xxx_zzz_1[i] + ta2_xy_xxx_zzz_0[i] * pa_x[i] - ta2_xy_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : GF

    auto ta2_xy_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 160);

    auto ta2_xy_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 161);

    auto ta2_xy_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 162);

    auto ta2_xy_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 163);

    auto ta2_xy_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 164);

    auto ta2_xy_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 165);

    auto ta2_xy_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 166);

    auto ta2_xy_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 167);

    auto ta2_xy_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 168);

    auto ta2_xy_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 169);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_1, ta1_x_xxx_zzz_1, ta1_y_xxy_yyy_1, ta1_y_xxy_yyz_1, ta1_y_xxy_yzz_1, ta2_xy_xxx_xx_0, ta2_xy_xxx_xx_1, ta2_xy_xxx_xxx_0, ta2_xy_xxx_xxx_1, ta2_xy_xxx_xxy_0, ta2_xy_xxx_xxy_1, ta2_xy_xxx_xxz_0, ta2_xy_xxx_xxz_1, ta2_xy_xxx_xy_0, ta2_xy_xxx_xy_1, ta2_xy_xxx_xyy_0, ta2_xy_xxx_xyy_1, ta2_xy_xxx_xyz_0, ta2_xy_xxx_xyz_1, ta2_xy_xxx_xz_0, ta2_xy_xxx_xz_1, ta2_xy_xxx_xzz_0, ta2_xy_xxx_xzz_1, ta2_xy_xxx_zzz_0, ta2_xy_xxx_zzz_1, ta2_xy_xxxy_xxx_0, ta2_xy_xxxy_xxy_0, ta2_xy_xxxy_xxz_0, ta2_xy_xxxy_xyy_0, ta2_xy_xxxy_xyz_0, ta2_xy_xxxy_xzz_0, ta2_xy_xxxy_yyy_0, ta2_xy_xxxy_yyz_0, ta2_xy_xxxy_yzz_0, ta2_xy_xxxy_zzz_0, ta2_xy_xxy_yyy_0, ta2_xy_xxy_yyy_1, ta2_xy_xxy_yyz_0, ta2_xy_xxy_yyz_1, ta2_xy_xxy_yzz_0, ta2_xy_xxy_yzz_1, ta2_xy_xy_yyy_0, ta2_xy_xy_yyy_1, ta2_xy_xy_yyz_0, ta2_xy_xy_yyz_1, ta2_xy_xy_yzz_0, ta2_xy_xy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxy_xxx_0[i] = ta1_x_xxx_xxx_1[i] + ta2_xy_xxx_xxx_0[i] * pa_y[i] - ta2_xy_xxx_xxx_1[i] * pc_y[i];

        ta2_xy_xxxy_xxy_0[i] = ta2_xy_xxx_xx_0[i] * fe_0 - ta2_xy_xxx_xx_1[i] * fe_0 + ta1_x_xxx_xxy_1[i] + ta2_xy_xxx_xxy_0[i] * pa_y[i] - ta2_xy_xxx_xxy_1[i] * pc_y[i];

        ta2_xy_xxxy_xxz_0[i] = ta1_x_xxx_xxz_1[i] + ta2_xy_xxx_xxz_0[i] * pa_y[i] - ta2_xy_xxx_xxz_1[i] * pc_y[i];

        ta2_xy_xxxy_xyy_0[i] = 2.0 * ta2_xy_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xy_1[i] * fe_0 + ta1_x_xxx_xyy_1[i] + ta2_xy_xxx_xyy_0[i] * pa_y[i] - ta2_xy_xxx_xyy_1[i] * pc_y[i];

        ta2_xy_xxxy_xyz_0[i] = ta2_xy_xxx_xz_0[i] * fe_0 - ta2_xy_xxx_xz_1[i] * fe_0 + ta1_x_xxx_xyz_1[i] + ta2_xy_xxx_xyz_0[i] * pa_y[i] - ta2_xy_xxx_xyz_1[i] * pc_y[i];

        ta2_xy_xxxy_xzz_0[i] = ta1_x_xxx_xzz_1[i] + ta2_xy_xxx_xzz_0[i] * pa_y[i] - ta2_xy_xxx_xzz_1[i] * pc_y[i];

        ta2_xy_xxxy_yyy_0[i] = 2.0 * ta2_xy_xy_yyy_0[i] * fe_0 - 2.0 * ta2_xy_xy_yyy_1[i] * fe_0 + ta1_y_xxy_yyy_1[i] + ta2_xy_xxy_yyy_0[i] * pa_x[i] - ta2_xy_xxy_yyy_1[i] * pc_x[i];

        ta2_xy_xxxy_yyz_0[i] = 2.0 * ta2_xy_xy_yyz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yyz_1[i] * fe_0 + ta1_y_xxy_yyz_1[i] + ta2_xy_xxy_yyz_0[i] * pa_x[i] - ta2_xy_xxy_yyz_1[i] * pc_x[i];

        ta2_xy_xxxy_yzz_0[i] = 2.0 * ta2_xy_xy_yzz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yzz_1[i] * fe_0 + ta1_y_xxy_yzz_1[i] + ta2_xy_xxy_yzz_0[i] * pa_x[i] - ta2_xy_xxy_yzz_1[i] * pc_x[i];

        ta2_xy_xxxy_zzz_0[i] = ta1_x_xxx_zzz_1[i] + ta2_xy_xxx_zzz_0[i] * pa_y[i] - ta2_xy_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : GF

    auto ta2_xy_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 170);

    auto ta2_xy_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 171);

    auto ta2_xy_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 172);

    auto ta2_xy_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 173);

    auto ta2_xy_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 174);

    auto ta2_xy_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 175);

    auto ta2_xy_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 176);

    auto ta2_xy_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 177);

    auto ta2_xy_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 178);

    auto ta2_xy_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 179);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_xxx_xx_0, ta2_xy_xxx_xx_1, ta2_xy_xxx_xxx_0, ta2_xy_xxx_xxx_1, ta2_xy_xxx_xxy_0, ta2_xy_xxx_xxy_1, ta2_xy_xxx_xxz_0, ta2_xy_xxx_xxz_1, ta2_xy_xxx_xy_0, ta2_xy_xxx_xy_1, ta2_xy_xxx_xyy_0, ta2_xy_xxx_xyy_1, ta2_xy_xxx_xyz_0, ta2_xy_xxx_xyz_1, ta2_xy_xxx_xz_0, ta2_xy_xxx_xz_1, ta2_xy_xxx_xzz_0, ta2_xy_xxx_xzz_1, ta2_xy_xxx_yy_0, ta2_xy_xxx_yy_1, ta2_xy_xxx_yyy_0, ta2_xy_xxx_yyy_1, ta2_xy_xxx_yyz_0, ta2_xy_xxx_yyz_1, ta2_xy_xxx_yz_0, ta2_xy_xxx_yz_1, ta2_xy_xxx_yzz_0, ta2_xy_xxx_yzz_1, ta2_xy_xxx_zz_0, ta2_xy_xxx_zz_1, ta2_xy_xxx_zzz_0, ta2_xy_xxx_zzz_1, ta2_xy_xxxz_xxx_0, ta2_xy_xxxz_xxy_0, ta2_xy_xxxz_xxz_0, ta2_xy_xxxz_xyy_0, ta2_xy_xxxz_xyz_0, ta2_xy_xxxz_xzz_0, ta2_xy_xxxz_yyy_0, ta2_xy_xxxz_yyz_0, ta2_xy_xxxz_yzz_0, ta2_xy_xxxz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxz_xxx_0[i] = ta2_xy_xxx_xxx_0[i] * pa_z[i] - ta2_xy_xxx_xxx_1[i] * pc_z[i];

        ta2_xy_xxxz_xxy_0[i] = ta2_xy_xxx_xxy_0[i] * pa_z[i] - ta2_xy_xxx_xxy_1[i] * pc_z[i];

        ta2_xy_xxxz_xxz_0[i] = ta2_xy_xxx_xx_0[i] * fe_0 - ta2_xy_xxx_xx_1[i] * fe_0 + ta2_xy_xxx_xxz_0[i] * pa_z[i] - ta2_xy_xxx_xxz_1[i] * pc_z[i];

        ta2_xy_xxxz_xyy_0[i] = ta2_xy_xxx_xyy_0[i] * pa_z[i] - ta2_xy_xxx_xyy_1[i] * pc_z[i];

        ta2_xy_xxxz_xyz_0[i] = ta2_xy_xxx_xy_0[i] * fe_0 - ta2_xy_xxx_xy_1[i] * fe_0 + ta2_xy_xxx_xyz_0[i] * pa_z[i] - ta2_xy_xxx_xyz_1[i] * pc_z[i];

        ta2_xy_xxxz_xzz_0[i] = 2.0 * ta2_xy_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xz_1[i] * fe_0 + ta2_xy_xxx_xzz_0[i] * pa_z[i] - ta2_xy_xxx_xzz_1[i] * pc_z[i];

        ta2_xy_xxxz_yyy_0[i] = ta2_xy_xxx_yyy_0[i] * pa_z[i] - ta2_xy_xxx_yyy_1[i] * pc_z[i];

        ta2_xy_xxxz_yyz_0[i] = ta2_xy_xxx_yy_0[i] * fe_0 - ta2_xy_xxx_yy_1[i] * fe_0 + ta2_xy_xxx_yyz_0[i] * pa_z[i] - ta2_xy_xxx_yyz_1[i] * pc_z[i];

        ta2_xy_xxxz_yzz_0[i] = 2.0 * ta2_xy_xxx_yz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_yz_1[i] * fe_0 + ta2_xy_xxx_yzz_0[i] * pa_z[i] - ta2_xy_xxx_yzz_1[i] * pc_z[i];

        ta2_xy_xxxz_zzz_0[i] = 3.0 * ta2_xy_xxx_zz_0[i] * fe_0 - 3.0 * ta2_xy_xxx_zz_1[i] * fe_0 + ta2_xy_xxx_zzz_0[i] * pa_z[i] - ta2_xy_xxx_zzz_1[i] * pc_z[i];
    }

    // Set up 180-190 components of targeted buffer : GF

    auto ta2_xy_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 180);

    auto ta2_xy_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 181);

    auto ta2_xy_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 182);

    auto ta2_xy_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 183);

    auto ta2_xy_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 184);

    auto ta2_xy_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 185);

    auto ta2_xy_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 186);

    auto ta2_xy_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 187);

    auto ta2_xy_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 188);

    auto ta2_xy_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 189);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxy_xxx_1, ta1_x_xxy_xxz_1, ta1_x_xxy_xzz_1, ta1_y_xyy_xxy_1, ta1_y_xyy_xyy_1, ta1_y_xyy_xyz_1, ta1_y_xyy_yyy_1, ta1_y_xyy_yyz_1, ta1_y_xyy_yzz_1, ta1_y_xyy_zzz_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xxy_xxx_0, ta2_xy_xxy_xxx_1, ta2_xy_xxy_xxz_0, ta2_xy_xxy_xxz_1, ta2_xy_xxy_xzz_0, ta2_xy_xxy_xzz_1, ta2_xy_xxyy_xxx_0, ta2_xy_xxyy_xxy_0, ta2_xy_xxyy_xxz_0, ta2_xy_xxyy_xyy_0, ta2_xy_xxyy_xyz_0, ta2_xy_xxyy_xzz_0, ta2_xy_xxyy_yyy_0, ta2_xy_xxyy_yyz_0, ta2_xy_xxyy_yzz_0, ta2_xy_xxyy_zzz_0, ta2_xy_xyy_xxy_0, ta2_xy_xyy_xxy_1, ta2_xy_xyy_xy_0, ta2_xy_xyy_xy_1, ta2_xy_xyy_xyy_0, ta2_xy_xyy_xyy_1, ta2_xy_xyy_xyz_0, ta2_xy_xyy_xyz_1, ta2_xy_xyy_yy_0, ta2_xy_xyy_yy_1, ta2_xy_xyy_yyy_0, ta2_xy_xyy_yyy_1, ta2_xy_xyy_yyz_0, ta2_xy_xyy_yyz_1, ta2_xy_xyy_yz_0, ta2_xy_xyy_yz_1, ta2_xy_xyy_yzz_0, ta2_xy_xyy_yzz_1, ta2_xy_xyy_zzz_0, ta2_xy_xyy_zzz_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyy_xxx_0[i] = ta2_xy_xx_xxx_0[i] * fe_0 - ta2_xy_xx_xxx_1[i] * fe_0 + ta1_x_xxy_xxx_1[i] + ta2_xy_xxy_xxx_0[i] * pa_y[i] - ta2_xy_xxy_xxx_1[i] * pc_y[i];

        ta2_xy_xxyy_xxy_0[i] = ta2_xy_yy_xxy_0[i] * fe_0 - ta2_xy_yy_xxy_1[i] * fe_0 + 2.0 * ta2_xy_xyy_xy_0[i] * fe_0 - 2.0 * ta2_xy_xyy_xy_1[i] * fe_0 + ta1_y_xyy_xxy_1[i] + ta2_xy_xyy_xxy_0[i] * pa_x[i] - ta2_xy_xyy_xxy_1[i] * pc_x[i];

        ta2_xy_xxyy_xxz_0[i] = ta2_xy_xx_xxz_0[i] * fe_0 - ta2_xy_xx_xxz_1[i] * fe_0 + ta1_x_xxy_xxz_1[i] + ta2_xy_xxy_xxz_0[i] * pa_y[i] - ta2_xy_xxy_xxz_1[i] * pc_y[i];

        ta2_xy_xxyy_xyy_0[i] = ta2_xy_yy_xyy_0[i] * fe_0 - ta2_xy_yy_xyy_1[i] * fe_0 + ta2_xy_xyy_yy_0[i] * fe_0 - ta2_xy_xyy_yy_1[i] * fe_0 + ta1_y_xyy_xyy_1[i] + ta2_xy_xyy_xyy_0[i] * pa_x[i] - ta2_xy_xyy_xyy_1[i] * pc_x[i];

        ta2_xy_xxyy_xyz_0[i] = ta2_xy_yy_xyz_0[i] * fe_0 - ta2_xy_yy_xyz_1[i] * fe_0 + ta2_xy_xyy_yz_0[i] * fe_0 - ta2_xy_xyy_yz_1[i] * fe_0 + ta1_y_xyy_xyz_1[i] + ta2_xy_xyy_xyz_0[i] * pa_x[i] - ta2_xy_xyy_xyz_1[i] * pc_x[i];

        ta2_xy_xxyy_xzz_0[i] = ta2_xy_xx_xzz_0[i] * fe_0 - ta2_xy_xx_xzz_1[i] * fe_0 + ta1_x_xxy_xzz_1[i] + ta2_xy_xxy_xzz_0[i] * pa_y[i] - ta2_xy_xxy_xzz_1[i] * pc_y[i];

        ta2_xy_xxyy_yyy_0[i] = ta2_xy_yy_yyy_0[i] * fe_0 - ta2_xy_yy_yyy_1[i] * fe_0 + ta1_y_xyy_yyy_1[i] + ta2_xy_xyy_yyy_0[i] * pa_x[i] - ta2_xy_xyy_yyy_1[i] * pc_x[i];

        ta2_xy_xxyy_yyz_0[i] = ta2_xy_yy_yyz_0[i] * fe_0 - ta2_xy_yy_yyz_1[i] * fe_0 + ta1_y_xyy_yyz_1[i] + ta2_xy_xyy_yyz_0[i] * pa_x[i] - ta2_xy_xyy_yyz_1[i] * pc_x[i];

        ta2_xy_xxyy_yzz_0[i] = ta2_xy_yy_yzz_0[i] * fe_0 - ta2_xy_yy_yzz_1[i] * fe_0 + ta1_y_xyy_yzz_1[i] + ta2_xy_xyy_yzz_0[i] * pa_x[i] - ta2_xy_xyy_yzz_1[i] * pc_x[i];

        ta2_xy_xxyy_zzz_0[i] = ta2_xy_yy_zzz_0[i] * fe_0 - ta2_xy_yy_zzz_1[i] * fe_0 + ta1_y_xyy_zzz_1[i] + ta2_xy_xyy_zzz_0[i] * pa_x[i] - ta2_xy_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 190-200 components of targeted buffer : GF

    auto ta2_xy_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 190);

    auto ta2_xy_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 191);

    auto ta2_xy_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 192);

    auto ta2_xy_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 193);

    auto ta2_xy_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 194);

    auto ta2_xy_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 195);

    auto ta2_xy_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 196);

    auto ta2_xy_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 197);

    auto ta2_xy_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 198);

    auto ta2_xy_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 199);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxz_xxz_1, ta1_x_xxz_xzz_1, ta1_x_xxz_zzz_1, ta2_xy_xxy_xxx_0, ta2_xy_xxy_xxx_1, ta2_xy_xxy_xxy_0, ta2_xy_xxy_xxy_1, ta2_xy_xxy_xy_0, ta2_xy_xxy_xy_1, ta2_xy_xxy_xyy_0, ta2_xy_xxy_xyy_1, ta2_xy_xxy_xyz_0, ta2_xy_xxy_xyz_1, ta2_xy_xxy_yy_0, ta2_xy_xxy_yy_1, ta2_xy_xxy_yyy_0, ta2_xy_xxy_yyy_1, ta2_xy_xxy_yyz_0, ta2_xy_xxy_yyz_1, ta2_xy_xxy_yz_0, ta2_xy_xxy_yz_1, ta2_xy_xxy_yzz_0, ta2_xy_xxy_yzz_1, ta2_xy_xxyz_xxx_0, ta2_xy_xxyz_xxy_0, ta2_xy_xxyz_xxz_0, ta2_xy_xxyz_xyy_0, ta2_xy_xxyz_xyz_0, ta2_xy_xxyz_xzz_0, ta2_xy_xxyz_yyy_0, ta2_xy_xxyz_yyz_0, ta2_xy_xxyz_yzz_0, ta2_xy_xxyz_zzz_0, ta2_xy_xxz_xxz_0, ta2_xy_xxz_xxz_1, ta2_xy_xxz_xzz_0, ta2_xy_xxz_xzz_1, ta2_xy_xxz_zzz_0, ta2_xy_xxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyz_xxx_0[i] = ta2_xy_xxy_xxx_0[i] * pa_z[i] - ta2_xy_xxy_xxx_1[i] * pc_z[i];

        ta2_xy_xxyz_xxy_0[i] = ta2_xy_xxy_xxy_0[i] * pa_z[i] - ta2_xy_xxy_xxy_1[i] * pc_z[i];

        ta2_xy_xxyz_xxz_0[i] = ta1_x_xxz_xxz_1[i] + ta2_xy_xxz_xxz_0[i] * pa_y[i] - ta2_xy_xxz_xxz_1[i] * pc_y[i];

        ta2_xy_xxyz_xyy_0[i] = ta2_xy_xxy_xyy_0[i] * pa_z[i] - ta2_xy_xxy_xyy_1[i] * pc_z[i];

        ta2_xy_xxyz_xyz_0[i] = ta2_xy_xxy_xy_0[i] * fe_0 - ta2_xy_xxy_xy_1[i] * fe_0 + ta2_xy_xxy_xyz_0[i] * pa_z[i] - ta2_xy_xxy_xyz_1[i] * pc_z[i];

        ta2_xy_xxyz_xzz_0[i] = ta1_x_xxz_xzz_1[i] + ta2_xy_xxz_xzz_0[i] * pa_y[i] - ta2_xy_xxz_xzz_1[i] * pc_y[i];

        ta2_xy_xxyz_yyy_0[i] = ta2_xy_xxy_yyy_0[i] * pa_z[i] - ta2_xy_xxy_yyy_1[i] * pc_z[i];

        ta2_xy_xxyz_yyz_0[i] = ta2_xy_xxy_yy_0[i] * fe_0 - ta2_xy_xxy_yy_1[i] * fe_0 + ta2_xy_xxy_yyz_0[i] * pa_z[i] - ta2_xy_xxy_yyz_1[i] * pc_z[i];

        ta2_xy_xxyz_yzz_0[i] = 2.0 * ta2_xy_xxy_yz_0[i] * fe_0 - 2.0 * ta2_xy_xxy_yz_1[i] * fe_0 + ta2_xy_xxy_yzz_0[i] * pa_z[i] - ta2_xy_xxy_yzz_1[i] * pc_z[i];

        ta2_xy_xxyz_zzz_0[i] = ta1_x_xxz_zzz_1[i] + ta2_xy_xxz_zzz_0[i] * pa_y[i] - ta2_xy_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 200-210 components of targeted buffer : GF

    auto ta2_xy_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 200);

    auto ta2_xy_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 201);

    auto ta2_xy_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 202);

    auto ta2_xy_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 203);

    auto ta2_xy_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 204);

    auto ta2_xy_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 205);

    auto ta2_xy_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 206);

    auto ta2_xy_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 207);

    auto ta2_xy_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 208);

    auto ta2_xy_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 209);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xzz_yyz_1, ta1_y_xzz_yzz_1, ta1_y_xzz_zzz_1, ta2_xy_xx_xxx_0, ta2_xy_xx_xxx_1, ta2_xy_xx_xxy_0, ta2_xy_xx_xxy_1, ta2_xy_xx_xxz_0, ta2_xy_xx_xxz_1, ta2_xy_xx_xyy_0, ta2_xy_xx_xyy_1, ta2_xy_xx_xyz_0, ta2_xy_xx_xyz_1, ta2_xy_xx_xzz_0, ta2_xy_xx_xzz_1, ta2_xy_xx_yyy_0, ta2_xy_xx_yyy_1, ta2_xy_xxz_xx_0, ta2_xy_xxz_xx_1, ta2_xy_xxz_xxx_0, ta2_xy_xxz_xxx_1, ta2_xy_xxz_xxy_0, ta2_xy_xxz_xxy_1, ta2_xy_xxz_xxz_0, ta2_xy_xxz_xxz_1, ta2_xy_xxz_xy_0, ta2_xy_xxz_xy_1, ta2_xy_xxz_xyy_0, ta2_xy_xxz_xyy_1, ta2_xy_xxz_xyz_0, ta2_xy_xxz_xyz_1, ta2_xy_xxz_xz_0, ta2_xy_xxz_xz_1, ta2_xy_xxz_xzz_0, ta2_xy_xxz_xzz_1, ta2_xy_xxz_yyy_0, ta2_xy_xxz_yyy_1, ta2_xy_xxzz_xxx_0, ta2_xy_xxzz_xxy_0, ta2_xy_xxzz_xxz_0, ta2_xy_xxzz_xyy_0, ta2_xy_xxzz_xyz_0, ta2_xy_xxzz_xzz_0, ta2_xy_xxzz_yyy_0, ta2_xy_xxzz_yyz_0, ta2_xy_xxzz_yzz_0, ta2_xy_xxzz_zzz_0, ta2_xy_xzz_yyz_0, ta2_xy_xzz_yyz_1, ta2_xy_xzz_yzz_0, ta2_xy_xzz_yzz_1, ta2_xy_xzz_zzz_0, ta2_xy_xzz_zzz_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxzz_xxx_0[i] = ta2_xy_xx_xxx_0[i] * fe_0 - ta2_xy_xx_xxx_1[i] * fe_0 + ta2_xy_xxz_xxx_0[i] * pa_z[i] - ta2_xy_xxz_xxx_1[i] * pc_z[i];

        ta2_xy_xxzz_xxy_0[i] = ta2_xy_xx_xxy_0[i] * fe_0 - ta2_xy_xx_xxy_1[i] * fe_0 + ta2_xy_xxz_xxy_0[i] * pa_z[i] - ta2_xy_xxz_xxy_1[i] * pc_z[i];

        ta2_xy_xxzz_xxz_0[i] = ta2_xy_xx_xxz_0[i] * fe_0 - ta2_xy_xx_xxz_1[i] * fe_0 + ta2_xy_xxz_xx_0[i] * fe_0 - ta2_xy_xxz_xx_1[i] * fe_0 + ta2_xy_xxz_xxz_0[i] * pa_z[i] - ta2_xy_xxz_xxz_1[i] * pc_z[i];

        ta2_xy_xxzz_xyy_0[i] = ta2_xy_xx_xyy_0[i] * fe_0 - ta2_xy_xx_xyy_1[i] * fe_0 + ta2_xy_xxz_xyy_0[i] * pa_z[i] - ta2_xy_xxz_xyy_1[i] * pc_z[i];

        ta2_xy_xxzz_xyz_0[i] = ta2_xy_xx_xyz_0[i] * fe_0 - ta2_xy_xx_xyz_1[i] * fe_0 + ta2_xy_xxz_xy_0[i] * fe_0 - ta2_xy_xxz_xy_1[i] * fe_0 + ta2_xy_xxz_xyz_0[i] * pa_z[i] - ta2_xy_xxz_xyz_1[i] * pc_z[i];

        ta2_xy_xxzz_xzz_0[i] = ta2_xy_xx_xzz_0[i] * fe_0 - ta2_xy_xx_xzz_1[i] * fe_0 + 2.0 * ta2_xy_xxz_xz_0[i] * fe_0 - 2.0 * ta2_xy_xxz_xz_1[i] * fe_0 + ta2_xy_xxz_xzz_0[i] * pa_z[i] - ta2_xy_xxz_xzz_1[i] * pc_z[i];

        ta2_xy_xxzz_yyy_0[i] = ta2_xy_xx_yyy_0[i] * fe_0 - ta2_xy_xx_yyy_1[i] * fe_0 + ta2_xy_xxz_yyy_0[i] * pa_z[i] - ta2_xy_xxz_yyy_1[i] * pc_z[i];

        ta2_xy_xxzz_yyz_0[i] = ta2_xy_zz_yyz_0[i] * fe_0 - ta2_xy_zz_yyz_1[i] * fe_0 + ta1_y_xzz_yyz_1[i] + ta2_xy_xzz_yyz_0[i] * pa_x[i] - ta2_xy_xzz_yyz_1[i] * pc_x[i];

        ta2_xy_xxzz_yzz_0[i] = ta2_xy_zz_yzz_0[i] * fe_0 - ta2_xy_zz_yzz_1[i] * fe_0 + ta1_y_xzz_yzz_1[i] + ta2_xy_xzz_yzz_0[i] * pa_x[i] - ta2_xy_xzz_yzz_1[i] * pc_x[i];

        ta2_xy_xxzz_zzz_0[i] = ta2_xy_zz_zzz_0[i] * fe_0 - ta2_xy_zz_zzz_1[i] * fe_0 + ta1_y_xzz_zzz_1[i] + ta2_xy_xzz_zzz_0[i] * pa_x[i] - ta2_xy_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : GF

    auto ta2_xy_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 210);

    auto ta2_xy_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 211);

    auto ta2_xy_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 212);

    auto ta2_xy_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 213);

    auto ta2_xy_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 214);

    auto ta2_xy_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 215);

    auto ta2_xy_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 216);

    auto ta2_xy_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 217);

    auto ta2_xy_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 218);

    auto ta2_xy_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 219);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_yyy_xxx_1, ta1_y_yyy_xxy_1, ta1_y_yyy_xxz_1, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_1, ta1_y_yyy_xzz_1, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_1, ta1_y_yyy_zzz_1, ta2_xy_xyyy_xxx_0, ta2_xy_xyyy_xxy_0, ta2_xy_xyyy_xxz_0, ta2_xy_xyyy_xyy_0, ta2_xy_xyyy_xyz_0, ta2_xy_xyyy_xzz_0, ta2_xy_xyyy_yyy_0, ta2_xy_xyyy_yyz_0, ta2_xy_xyyy_yzz_0, ta2_xy_xyyy_zzz_0, ta2_xy_yyy_xx_0, ta2_xy_yyy_xx_1, ta2_xy_yyy_xxx_0, ta2_xy_yyy_xxx_1, ta2_xy_yyy_xxy_0, ta2_xy_yyy_xxy_1, ta2_xy_yyy_xxz_0, ta2_xy_yyy_xxz_1, ta2_xy_yyy_xy_0, ta2_xy_yyy_xy_1, ta2_xy_yyy_xyy_0, ta2_xy_yyy_xyy_1, ta2_xy_yyy_xyz_0, ta2_xy_yyy_xyz_1, ta2_xy_yyy_xz_0, ta2_xy_yyy_xz_1, ta2_xy_yyy_xzz_0, ta2_xy_yyy_xzz_1, ta2_xy_yyy_yy_0, ta2_xy_yyy_yy_1, ta2_xy_yyy_yyy_0, ta2_xy_yyy_yyy_1, ta2_xy_yyy_yyz_0, ta2_xy_yyy_yyz_1, ta2_xy_yyy_yz_0, ta2_xy_yyy_yz_1, ta2_xy_yyy_yzz_0, ta2_xy_yyy_yzz_1, ta2_xy_yyy_zz_0, ta2_xy_yyy_zz_1, ta2_xy_yyy_zzz_0, ta2_xy_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyy_xxx_0[i] = 3.0 * ta2_xy_yyy_xx_0[i] * fe_0 - 3.0 * ta2_xy_yyy_xx_1[i] * fe_0 + ta1_y_yyy_xxx_1[i] + ta2_xy_yyy_xxx_0[i] * pa_x[i] - ta2_xy_yyy_xxx_1[i] * pc_x[i];

        ta2_xy_xyyy_xxy_0[i] = 2.0 * ta2_xy_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xy_1[i] * fe_0 + ta1_y_yyy_xxy_1[i] + ta2_xy_yyy_xxy_0[i] * pa_x[i] - ta2_xy_yyy_xxy_1[i] * pc_x[i];

        ta2_xy_xyyy_xxz_0[i] = 2.0 * ta2_xy_yyy_xz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xz_1[i] * fe_0 + ta1_y_yyy_xxz_1[i] + ta2_xy_yyy_xxz_0[i] * pa_x[i] - ta2_xy_yyy_xxz_1[i] * pc_x[i];

        ta2_xy_xyyy_xyy_0[i] = ta2_xy_yyy_yy_0[i] * fe_0 - ta2_xy_yyy_yy_1[i] * fe_0 + ta1_y_yyy_xyy_1[i] + ta2_xy_yyy_xyy_0[i] * pa_x[i] - ta2_xy_yyy_xyy_1[i] * pc_x[i];

        ta2_xy_xyyy_xyz_0[i] = ta2_xy_yyy_yz_0[i] * fe_0 - ta2_xy_yyy_yz_1[i] * fe_0 + ta1_y_yyy_xyz_1[i] + ta2_xy_yyy_xyz_0[i] * pa_x[i] - ta2_xy_yyy_xyz_1[i] * pc_x[i];

        ta2_xy_xyyy_xzz_0[i] = ta2_xy_yyy_zz_0[i] * fe_0 - ta2_xy_yyy_zz_1[i] * fe_0 + ta1_y_yyy_xzz_1[i] + ta2_xy_yyy_xzz_0[i] * pa_x[i] - ta2_xy_yyy_xzz_1[i] * pc_x[i];

        ta2_xy_xyyy_yyy_0[i] = ta1_y_yyy_yyy_1[i] + ta2_xy_yyy_yyy_0[i] * pa_x[i] - ta2_xy_yyy_yyy_1[i] * pc_x[i];

        ta2_xy_xyyy_yyz_0[i] = ta1_y_yyy_yyz_1[i] + ta2_xy_yyy_yyz_0[i] * pa_x[i] - ta2_xy_yyy_yyz_1[i] * pc_x[i];

        ta2_xy_xyyy_yzz_0[i] = ta1_y_yyy_yzz_1[i] + ta2_xy_yyy_yzz_0[i] * pa_x[i] - ta2_xy_yyy_yzz_1[i] * pc_x[i];

        ta2_xy_xyyy_zzz_0[i] = ta1_y_yyy_zzz_1[i] + ta2_xy_yyy_zzz_0[i] * pa_x[i] - ta2_xy_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 220-230 components of targeted buffer : GF

    auto ta2_xy_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 220);

    auto ta2_xy_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 221);

    auto ta2_xy_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 222);

    auto ta2_xy_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 223);

    auto ta2_xy_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 224);

    auto ta2_xy_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 225);

    auto ta2_xy_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 226);

    auto ta2_xy_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 227);

    auto ta2_xy_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 228);

    auto ta2_xy_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 229);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_yyz_yyz_1, ta1_y_yyz_yzz_1, ta1_y_yyz_zzz_1, ta2_xy_xyy_xx_0, ta2_xy_xyy_xx_1, ta2_xy_xyy_xxx_0, ta2_xy_xyy_xxx_1, ta2_xy_xyy_xxy_0, ta2_xy_xyy_xxy_1, ta2_xy_xyy_xxz_0, ta2_xy_xyy_xxz_1, ta2_xy_xyy_xy_0, ta2_xy_xyy_xy_1, ta2_xy_xyy_xyy_0, ta2_xy_xyy_xyy_1, ta2_xy_xyy_xyz_0, ta2_xy_xyy_xyz_1, ta2_xy_xyy_xz_0, ta2_xy_xyy_xz_1, ta2_xy_xyy_xzz_0, ta2_xy_xyy_xzz_1, ta2_xy_xyy_yyy_0, ta2_xy_xyy_yyy_1, ta2_xy_xyyz_xxx_0, ta2_xy_xyyz_xxy_0, ta2_xy_xyyz_xxz_0, ta2_xy_xyyz_xyy_0, ta2_xy_xyyz_xyz_0, ta2_xy_xyyz_xzz_0, ta2_xy_xyyz_yyy_0, ta2_xy_xyyz_yyz_0, ta2_xy_xyyz_yzz_0, ta2_xy_xyyz_zzz_0, ta2_xy_yyz_yyz_0, ta2_xy_yyz_yyz_1, ta2_xy_yyz_yzz_0, ta2_xy_yyz_yzz_1, ta2_xy_yyz_zzz_0, ta2_xy_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyz_xxx_0[i] = ta2_xy_xyy_xxx_0[i] * pa_z[i] - ta2_xy_xyy_xxx_1[i] * pc_z[i];

        ta2_xy_xyyz_xxy_0[i] = ta2_xy_xyy_xxy_0[i] * pa_z[i] - ta2_xy_xyy_xxy_1[i] * pc_z[i];

        ta2_xy_xyyz_xxz_0[i] = ta2_xy_xyy_xx_0[i] * fe_0 - ta2_xy_xyy_xx_1[i] * fe_0 + ta2_xy_xyy_xxz_0[i] * pa_z[i] - ta2_xy_xyy_xxz_1[i] * pc_z[i];

        ta2_xy_xyyz_xyy_0[i] = ta2_xy_xyy_xyy_0[i] * pa_z[i] - ta2_xy_xyy_xyy_1[i] * pc_z[i];

        ta2_xy_xyyz_xyz_0[i] = ta2_xy_xyy_xy_0[i] * fe_0 - ta2_xy_xyy_xy_1[i] * fe_0 + ta2_xy_xyy_xyz_0[i] * pa_z[i] - ta2_xy_xyy_xyz_1[i] * pc_z[i];

        ta2_xy_xyyz_xzz_0[i] = 2.0 * ta2_xy_xyy_xz_0[i] * fe_0 - 2.0 * ta2_xy_xyy_xz_1[i] * fe_0 + ta2_xy_xyy_xzz_0[i] * pa_z[i] - ta2_xy_xyy_xzz_1[i] * pc_z[i];

        ta2_xy_xyyz_yyy_0[i] = ta2_xy_xyy_yyy_0[i] * pa_z[i] - ta2_xy_xyy_yyy_1[i] * pc_z[i];

        ta2_xy_xyyz_yyz_0[i] = ta1_y_yyz_yyz_1[i] + ta2_xy_yyz_yyz_0[i] * pa_x[i] - ta2_xy_yyz_yyz_1[i] * pc_x[i];

        ta2_xy_xyyz_yzz_0[i] = ta1_y_yyz_yzz_1[i] + ta2_xy_yyz_yzz_0[i] * pa_x[i] - ta2_xy_yyz_yzz_1[i] * pc_x[i];

        ta2_xy_xyyz_zzz_0[i] = ta1_y_yyz_zzz_1[i] + ta2_xy_yyz_zzz_0[i] * pa_x[i] - ta2_xy_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 230-240 components of targeted buffer : GF

    auto ta2_xy_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 230);

    auto ta2_xy_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 231);

    auto ta2_xy_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 232);

    auto ta2_xy_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 233);

    auto ta2_xy_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 234);

    auto ta2_xy_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 235);

    auto ta2_xy_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 236);

    auto ta2_xy_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 237);

    auto ta2_xy_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 238);

    auto ta2_xy_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 239);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xzz_xxx_1, ta1_x_xzz_xxz_1, ta1_x_xzz_xzz_1, ta1_y_yzz_xyz_1, ta1_y_yzz_yyy_1, ta1_y_yzz_yyz_1, ta1_y_yzz_yzz_1, ta1_y_yzz_zzz_1, ta2_xy_xy_xxy_0, ta2_xy_xy_xxy_1, ta2_xy_xy_xyy_0, ta2_xy_xy_xyy_1, ta2_xy_xyz_xxy_0, ta2_xy_xyz_xxy_1, ta2_xy_xyz_xyy_0, ta2_xy_xyz_xyy_1, ta2_xy_xyzz_xxx_0, ta2_xy_xyzz_xxy_0, ta2_xy_xyzz_xxz_0, ta2_xy_xyzz_xyy_0, ta2_xy_xyzz_xyz_0, ta2_xy_xyzz_xzz_0, ta2_xy_xyzz_yyy_0, ta2_xy_xyzz_yyz_0, ta2_xy_xyzz_yzz_0, ta2_xy_xyzz_zzz_0, ta2_xy_xzz_xxx_0, ta2_xy_xzz_xxx_1, ta2_xy_xzz_xxz_0, ta2_xy_xzz_xxz_1, ta2_xy_xzz_xzz_0, ta2_xy_xzz_xzz_1, ta2_xy_yzz_xyz_0, ta2_xy_yzz_xyz_1, ta2_xy_yzz_yyy_0, ta2_xy_yzz_yyy_1, ta2_xy_yzz_yyz_0, ta2_xy_yzz_yyz_1, ta2_xy_yzz_yz_0, ta2_xy_yzz_yz_1, ta2_xy_yzz_yzz_0, ta2_xy_yzz_yzz_1, ta2_xy_yzz_zzz_0, ta2_xy_yzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyzz_xxx_0[i] = ta1_x_xzz_xxx_1[i] + ta2_xy_xzz_xxx_0[i] * pa_y[i] - ta2_xy_xzz_xxx_1[i] * pc_y[i];

        ta2_xy_xyzz_xxy_0[i] = ta2_xy_xy_xxy_0[i] * fe_0 - ta2_xy_xy_xxy_1[i] * fe_0 + ta2_xy_xyz_xxy_0[i] * pa_z[i] - ta2_xy_xyz_xxy_1[i] * pc_z[i];

        ta2_xy_xyzz_xxz_0[i] = ta1_x_xzz_xxz_1[i] + ta2_xy_xzz_xxz_0[i] * pa_y[i] - ta2_xy_xzz_xxz_1[i] * pc_y[i];

        ta2_xy_xyzz_xyy_0[i] = ta2_xy_xy_xyy_0[i] * fe_0 - ta2_xy_xy_xyy_1[i] * fe_0 + ta2_xy_xyz_xyy_0[i] * pa_z[i] - ta2_xy_xyz_xyy_1[i] * pc_z[i];

        ta2_xy_xyzz_xyz_0[i] = ta2_xy_yzz_yz_0[i] * fe_0 - ta2_xy_yzz_yz_1[i] * fe_0 + ta1_y_yzz_xyz_1[i] + ta2_xy_yzz_xyz_0[i] * pa_x[i] - ta2_xy_yzz_xyz_1[i] * pc_x[i];

        ta2_xy_xyzz_xzz_0[i] = ta1_x_xzz_xzz_1[i] + ta2_xy_xzz_xzz_0[i] * pa_y[i] - ta2_xy_xzz_xzz_1[i] * pc_y[i];

        ta2_xy_xyzz_yyy_0[i] = ta1_y_yzz_yyy_1[i] + ta2_xy_yzz_yyy_0[i] * pa_x[i] - ta2_xy_yzz_yyy_1[i] * pc_x[i];

        ta2_xy_xyzz_yyz_0[i] = ta1_y_yzz_yyz_1[i] + ta2_xy_yzz_yyz_0[i] * pa_x[i] - ta2_xy_yzz_yyz_1[i] * pc_x[i];

        ta2_xy_xyzz_yzz_0[i] = ta1_y_yzz_yzz_1[i] + ta2_xy_yzz_yzz_0[i] * pa_x[i] - ta2_xy_yzz_yzz_1[i] * pc_x[i];

        ta2_xy_xyzz_zzz_0[i] = ta1_y_yzz_zzz_1[i] + ta2_xy_yzz_zzz_0[i] * pa_x[i] - ta2_xy_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 240-250 components of targeted buffer : GF

    auto ta2_xy_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 240);

    auto ta2_xy_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 241);

    auto ta2_xy_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 242);

    auto ta2_xy_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 243);

    auto ta2_xy_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 244);

    auto ta2_xy_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 245);

    auto ta2_xy_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 246);

    auto ta2_xy_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 247);

    auto ta2_xy_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 248);

    auto ta2_xy_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 249);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_zzz_xxz_1, ta1_y_zzz_xyz_1, ta1_y_zzz_xzz_1, ta1_y_zzz_yyy_1, ta1_y_zzz_yyz_1, ta1_y_zzz_yzz_1, ta1_y_zzz_zzz_1, ta2_xy_xz_xxx_0, ta2_xy_xz_xxx_1, ta2_xy_xz_xxy_0, ta2_xy_xz_xxy_1, ta2_xy_xz_xyy_0, ta2_xy_xz_xyy_1, ta2_xy_xzz_xxx_0, ta2_xy_xzz_xxx_1, ta2_xy_xzz_xxy_0, ta2_xy_xzz_xxy_1, ta2_xy_xzz_xyy_0, ta2_xy_xzz_xyy_1, ta2_xy_xzzz_xxx_0, ta2_xy_xzzz_xxy_0, ta2_xy_xzzz_xxz_0, ta2_xy_xzzz_xyy_0, ta2_xy_xzzz_xyz_0, ta2_xy_xzzz_xzz_0, ta2_xy_xzzz_yyy_0, ta2_xy_xzzz_yyz_0, ta2_xy_xzzz_yzz_0, ta2_xy_xzzz_zzz_0, ta2_xy_zzz_xxz_0, ta2_xy_zzz_xxz_1, ta2_xy_zzz_xyz_0, ta2_xy_zzz_xyz_1, ta2_xy_zzz_xz_0, ta2_xy_zzz_xz_1, ta2_xy_zzz_xzz_0, ta2_xy_zzz_xzz_1, ta2_xy_zzz_yyy_0, ta2_xy_zzz_yyy_1, ta2_xy_zzz_yyz_0, ta2_xy_zzz_yyz_1, ta2_xy_zzz_yz_0, ta2_xy_zzz_yz_1, ta2_xy_zzz_yzz_0, ta2_xy_zzz_yzz_1, ta2_xy_zzz_zz_0, ta2_xy_zzz_zz_1, ta2_xy_zzz_zzz_0, ta2_xy_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzzz_xxx_0[i] = 2.0 * ta2_xy_xz_xxx_0[i] * fe_0 - 2.0 * ta2_xy_xz_xxx_1[i] * fe_0 + ta2_xy_xzz_xxx_0[i] * pa_z[i] - ta2_xy_xzz_xxx_1[i] * pc_z[i];

        ta2_xy_xzzz_xxy_0[i] = 2.0 * ta2_xy_xz_xxy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xxy_1[i] * fe_0 + ta2_xy_xzz_xxy_0[i] * pa_z[i] - ta2_xy_xzz_xxy_1[i] * pc_z[i];

        ta2_xy_xzzz_xxz_0[i] = 2.0 * ta2_xy_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_xz_1[i] * fe_0 + ta1_y_zzz_xxz_1[i] + ta2_xy_zzz_xxz_0[i] * pa_x[i] - ta2_xy_zzz_xxz_1[i] * pc_x[i];

        ta2_xy_xzzz_xyy_0[i] = 2.0 * ta2_xy_xz_xyy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xyy_1[i] * fe_0 + ta2_xy_xzz_xyy_0[i] * pa_z[i] - ta2_xy_xzz_xyy_1[i] * pc_z[i];

        ta2_xy_xzzz_xyz_0[i] = ta2_xy_zzz_yz_0[i] * fe_0 - ta2_xy_zzz_yz_1[i] * fe_0 + ta1_y_zzz_xyz_1[i] + ta2_xy_zzz_xyz_0[i] * pa_x[i] - ta2_xy_zzz_xyz_1[i] * pc_x[i];

        ta2_xy_xzzz_xzz_0[i] = ta2_xy_zzz_zz_0[i] * fe_0 - ta2_xy_zzz_zz_1[i] * fe_0 + ta1_y_zzz_xzz_1[i] + ta2_xy_zzz_xzz_0[i] * pa_x[i] - ta2_xy_zzz_xzz_1[i] * pc_x[i];

        ta2_xy_xzzz_yyy_0[i] = ta1_y_zzz_yyy_1[i] + ta2_xy_zzz_yyy_0[i] * pa_x[i] - ta2_xy_zzz_yyy_1[i] * pc_x[i];

        ta2_xy_xzzz_yyz_0[i] = ta1_y_zzz_yyz_1[i] + ta2_xy_zzz_yyz_0[i] * pa_x[i] - ta2_xy_zzz_yyz_1[i] * pc_x[i];

        ta2_xy_xzzz_yzz_0[i] = ta1_y_zzz_yzz_1[i] + ta2_xy_zzz_yzz_0[i] * pa_x[i] - ta2_xy_zzz_yzz_1[i] * pc_x[i];

        ta2_xy_xzzz_zzz_0[i] = ta1_y_zzz_zzz_1[i] + ta2_xy_zzz_zzz_0[i] * pa_x[i] - ta2_xy_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 250-260 components of targeted buffer : GF

    auto ta2_xy_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 250);

    auto ta2_xy_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 251);

    auto ta2_xy_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 252);

    auto ta2_xy_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 253);

    auto ta2_xy_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 254);

    auto ta2_xy_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 255);

    auto ta2_xy_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 256);

    auto ta2_xy_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 257);

    auto ta2_xy_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 258);

    auto ta2_xy_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 259);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yyy_xxx_1, ta1_x_yyy_xxy_1, ta1_x_yyy_xxz_1, ta1_x_yyy_xyy_1, ta1_x_yyy_xyz_1, ta1_x_yyy_xzz_1, ta1_x_yyy_yyy_1, ta1_x_yyy_yyz_1, ta1_x_yyy_yzz_1, ta1_x_yyy_zzz_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xxz_0, ta2_xy_yy_xxz_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_xzz_0, ta2_xy_yy_xzz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yy_zzz_0, ta2_xy_yy_zzz_1, ta2_xy_yyy_xx_0, ta2_xy_yyy_xx_1, ta2_xy_yyy_xxx_0, ta2_xy_yyy_xxx_1, ta2_xy_yyy_xxy_0, ta2_xy_yyy_xxy_1, ta2_xy_yyy_xxz_0, ta2_xy_yyy_xxz_1, ta2_xy_yyy_xy_0, ta2_xy_yyy_xy_1, ta2_xy_yyy_xyy_0, ta2_xy_yyy_xyy_1, ta2_xy_yyy_xyz_0, ta2_xy_yyy_xyz_1, ta2_xy_yyy_xz_0, ta2_xy_yyy_xz_1, ta2_xy_yyy_xzz_0, ta2_xy_yyy_xzz_1, ta2_xy_yyy_yy_0, ta2_xy_yyy_yy_1, ta2_xy_yyy_yyy_0, ta2_xy_yyy_yyy_1, ta2_xy_yyy_yyz_0, ta2_xy_yyy_yyz_1, ta2_xy_yyy_yz_0, ta2_xy_yyy_yz_1, ta2_xy_yyy_yzz_0, ta2_xy_yyy_yzz_1, ta2_xy_yyy_zz_0, ta2_xy_yyy_zz_1, ta2_xy_yyy_zzz_0, ta2_xy_yyy_zzz_1, ta2_xy_yyyy_xxx_0, ta2_xy_yyyy_xxy_0, ta2_xy_yyyy_xxz_0, ta2_xy_yyyy_xyy_0, ta2_xy_yyyy_xyz_0, ta2_xy_yyyy_xzz_0, ta2_xy_yyyy_yyy_0, ta2_xy_yyyy_yyz_0, ta2_xy_yyyy_yzz_0, ta2_xy_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyy_xxx_0[i] = 3.0 * ta2_xy_yy_xxx_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxx_1[i] * fe_0 + ta1_x_yyy_xxx_1[i] + ta2_xy_yyy_xxx_0[i] * pa_y[i] - ta2_xy_yyy_xxx_1[i] * pc_y[i];

        ta2_xy_yyyy_xxy_0[i] = 3.0 * ta2_xy_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxy_1[i] * fe_0 + ta2_xy_yyy_xx_0[i] * fe_0 - ta2_xy_yyy_xx_1[i] * fe_0 + ta1_x_yyy_xxy_1[i] + ta2_xy_yyy_xxy_0[i] * pa_y[i] - ta2_xy_yyy_xxy_1[i] * pc_y[i];

        ta2_xy_yyyy_xxz_0[i] = 3.0 * ta2_xy_yy_xxz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxz_1[i] * fe_0 + ta1_x_yyy_xxz_1[i] + ta2_xy_yyy_xxz_0[i] * pa_y[i] - ta2_xy_yyy_xxz_1[i] * pc_y[i];

        ta2_xy_yyyy_xyy_0[i] = 3.0 * ta2_xy_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyy_1[i] * fe_0 + 2.0 * ta2_xy_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xy_1[i] * fe_0 + ta1_x_yyy_xyy_1[i] + ta2_xy_yyy_xyy_0[i] * pa_y[i] - ta2_xy_yyy_xyy_1[i] * pc_y[i];

        ta2_xy_yyyy_xyz_0[i] = 3.0 * ta2_xy_yy_xyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyz_1[i] * fe_0 + ta2_xy_yyy_xz_0[i] * fe_0 - ta2_xy_yyy_xz_1[i] * fe_0 + ta1_x_yyy_xyz_1[i] + ta2_xy_yyy_xyz_0[i] * pa_y[i] - ta2_xy_yyy_xyz_1[i] * pc_y[i];

        ta2_xy_yyyy_xzz_0[i] = 3.0 * ta2_xy_yy_xzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xzz_1[i] * fe_0 + ta1_x_yyy_xzz_1[i] + ta2_xy_yyy_xzz_0[i] * pa_y[i] - ta2_xy_yyy_xzz_1[i] * pc_y[i];

        ta2_xy_yyyy_yyy_0[i] = 3.0 * ta2_xy_yy_yyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyy_1[i] * fe_0 + 3.0 * ta2_xy_yyy_yy_0[i] * fe_0 - 3.0 * ta2_xy_yyy_yy_1[i] * fe_0 + ta1_x_yyy_yyy_1[i] + ta2_xy_yyy_yyy_0[i] * pa_y[i] - ta2_xy_yyy_yyy_1[i] * pc_y[i];

        ta2_xy_yyyy_yyz_0[i] = 3.0 * ta2_xy_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyz_1[i] * fe_0 + 2.0 * ta2_xy_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_yz_1[i] * fe_0 + ta1_x_yyy_yyz_1[i] + ta2_xy_yyy_yyz_0[i] * pa_y[i] - ta2_xy_yyy_yyz_1[i] * pc_y[i];

        ta2_xy_yyyy_yzz_0[i] = 3.0 * ta2_xy_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yzz_1[i] * fe_0 + ta2_xy_yyy_zz_0[i] * fe_0 - ta2_xy_yyy_zz_1[i] * fe_0 + ta1_x_yyy_yzz_1[i] + ta2_xy_yyy_yzz_0[i] * pa_y[i] - ta2_xy_yyy_yzz_1[i] * pc_y[i];

        ta2_xy_yyyy_zzz_0[i] = 3.0 * ta2_xy_yy_zzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_zzz_1[i] * fe_0 + ta1_x_yyy_zzz_1[i] + ta2_xy_yyy_zzz_0[i] * pa_y[i] - ta2_xy_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 260-270 components of targeted buffer : GF

    auto ta2_xy_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 260);

    auto ta2_xy_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 261);

    auto ta2_xy_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 262);

    auto ta2_xy_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 263);

    auto ta2_xy_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 264);

    auto ta2_xy_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 265);

    auto ta2_xy_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 266);

    auto ta2_xy_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 267);

    auto ta2_xy_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 268);

    auto ta2_xy_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 269);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_yyy_xx_0, ta2_xy_yyy_xx_1, ta2_xy_yyy_xxx_0, ta2_xy_yyy_xxx_1, ta2_xy_yyy_xxy_0, ta2_xy_yyy_xxy_1, ta2_xy_yyy_xxz_0, ta2_xy_yyy_xxz_1, ta2_xy_yyy_xy_0, ta2_xy_yyy_xy_1, ta2_xy_yyy_xyy_0, ta2_xy_yyy_xyy_1, ta2_xy_yyy_xyz_0, ta2_xy_yyy_xyz_1, ta2_xy_yyy_xz_0, ta2_xy_yyy_xz_1, ta2_xy_yyy_xzz_0, ta2_xy_yyy_xzz_1, ta2_xy_yyy_yy_0, ta2_xy_yyy_yy_1, ta2_xy_yyy_yyy_0, ta2_xy_yyy_yyy_1, ta2_xy_yyy_yyz_0, ta2_xy_yyy_yyz_1, ta2_xy_yyy_yz_0, ta2_xy_yyy_yz_1, ta2_xy_yyy_yzz_0, ta2_xy_yyy_yzz_1, ta2_xy_yyy_zz_0, ta2_xy_yyy_zz_1, ta2_xy_yyy_zzz_0, ta2_xy_yyy_zzz_1, ta2_xy_yyyz_xxx_0, ta2_xy_yyyz_xxy_0, ta2_xy_yyyz_xxz_0, ta2_xy_yyyz_xyy_0, ta2_xy_yyyz_xyz_0, ta2_xy_yyyz_xzz_0, ta2_xy_yyyz_yyy_0, ta2_xy_yyyz_yyz_0, ta2_xy_yyyz_yzz_0, ta2_xy_yyyz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyz_xxx_0[i] = ta2_xy_yyy_xxx_0[i] * pa_z[i] - ta2_xy_yyy_xxx_1[i] * pc_z[i];

        ta2_xy_yyyz_xxy_0[i] = ta2_xy_yyy_xxy_0[i] * pa_z[i] - ta2_xy_yyy_xxy_1[i] * pc_z[i];

        ta2_xy_yyyz_xxz_0[i] = ta2_xy_yyy_xx_0[i] * fe_0 - ta2_xy_yyy_xx_1[i] * fe_0 + ta2_xy_yyy_xxz_0[i] * pa_z[i] - ta2_xy_yyy_xxz_1[i] * pc_z[i];

        ta2_xy_yyyz_xyy_0[i] = ta2_xy_yyy_xyy_0[i] * pa_z[i] - ta2_xy_yyy_xyy_1[i] * pc_z[i];

        ta2_xy_yyyz_xyz_0[i] = ta2_xy_yyy_xy_0[i] * fe_0 - ta2_xy_yyy_xy_1[i] * fe_0 + ta2_xy_yyy_xyz_0[i] * pa_z[i] - ta2_xy_yyy_xyz_1[i] * pc_z[i];

        ta2_xy_yyyz_xzz_0[i] = 2.0 * ta2_xy_yyy_xz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xz_1[i] * fe_0 + ta2_xy_yyy_xzz_0[i] * pa_z[i] - ta2_xy_yyy_xzz_1[i] * pc_z[i];

        ta2_xy_yyyz_yyy_0[i] = ta2_xy_yyy_yyy_0[i] * pa_z[i] - ta2_xy_yyy_yyy_1[i] * pc_z[i];

        ta2_xy_yyyz_yyz_0[i] = ta2_xy_yyy_yy_0[i] * fe_0 - ta2_xy_yyy_yy_1[i] * fe_0 + ta2_xy_yyy_yyz_0[i] * pa_z[i] - ta2_xy_yyy_yyz_1[i] * pc_z[i];

        ta2_xy_yyyz_yzz_0[i] = 2.0 * ta2_xy_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_yz_1[i] * fe_0 + ta2_xy_yyy_yzz_0[i] * pa_z[i] - ta2_xy_yyy_yzz_1[i] * pc_z[i];

        ta2_xy_yyyz_zzz_0[i] = 3.0 * ta2_xy_yyy_zz_0[i] * fe_0 - 3.0 * ta2_xy_yyy_zz_1[i] * fe_0 + ta2_xy_yyy_zzz_0[i] * pa_z[i] - ta2_xy_yyy_zzz_1[i] * pc_z[i];
    }

    // Set up 270-280 components of targeted buffer : GF

    auto ta2_xy_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 270);

    auto ta2_xy_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 271);

    auto ta2_xy_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 272);

    auto ta2_xy_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 273);

    auto ta2_xy_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 274);

    auto ta2_xy_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 275);

    auto ta2_xy_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 276);

    auto ta2_xy_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 277);

    auto ta2_xy_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 278);

    auto ta2_xy_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 279);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yzz_xxz_1, ta1_x_yzz_xzz_1, ta1_x_yzz_zzz_1, ta2_xy_yy_xxx_0, ta2_xy_yy_xxx_1, ta2_xy_yy_xxy_0, ta2_xy_yy_xxy_1, ta2_xy_yy_xyy_0, ta2_xy_yy_xyy_1, ta2_xy_yy_xyz_0, ta2_xy_yy_xyz_1, ta2_xy_yy_yyy_0, ta2_xy_yy_yyy_1, ta2_xy_yy_yyz_0, ta2_xy_yy_yyz_1, ta2_xy_yy_yzz_0, ta2_xy_yy_yzz_1, ta2_xy_yyz_xxx_0, ta2_xy_yyz_xxx_1, ta2_xy_yyz_xxy_0, ta2_xy_yyz_xxy_1, ta2_xy_yyz_xy_0, ta2_xy_yyz_xy_1, ta2_xy_yyz_xyy_0, ta2_xy_yyz_xyy_1, ta2_xy_yyz_xyz_0, ta2_xy_yyz_xyz_1, ta2_xy_yyz_yy_0, ta2_xy_yyz_yy_1, ta2_xy_yyz_yyy_0, ta2_xy_yyz_yyy_1, ta2_xy_yyz_yyz_0, ta2_xy_yyz_yyz_1, ta2_xy_yyz_yz_0, ta2_xy_yyz_yz_1, ta2_xy_yyz_yzz_0, ta2_xy_yyz_yzz_1, ta2_xy_yyzz_xxx_0, ta2_xy_yyzz_xxy_0, ta2_xy_yyzz_xxz_0, ta2_xy_yyzz_xyy_0, ta2_xy_yyzz_xyz_0, ta2_xy_yyzz_xzz_0, ta2_xy_yyzz_yyy_0, ta2_xy_yyzz_yyz_0, ta2_xy_yyzz_yzz_0, ta2_xy_yyzz_zzz_0, ta2_xy_yzz_xxz_0, ta2_xy_yzz_xxz_1, ta2_xy_yzz_xzz_0, ta2_xy_yzz_xzz_1, ta2_xy_yzz_zzz_0, ta2_xy_yzz_zzz_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyzz_xxx_0[i] = ta2_xy_yy_xxx_0[i] * fe_0 - ta2_xy_yy_xxx_1[i] * fe_0 + ta2_xy_yyz_xxx_0[i] * pa_z[i] - ta2_xy_yyz_xxx_1[i] * pc_z[i];

        ta2_xy_yyzz_xxy_0[i] = ta2_xy_yy_xxy_0[i] * fe_0 - ta2_xy_yy_xxy_1[i] * fe_0 + ta2_xy_yyz_xxy_0[i] * pa_z[i] - ta2_xy_yyz_xxy_1[i] * pc_z[i];

        ta2_xy_yyzz_xxz_0[i] = ta2_xy_zz_xxz_0[i] * fe_0 - ta2_xy_zz_xxz_1[i] * fe_0 + ta1_x_yzz_xxz_1[i] + ta2_xy_yzz_xxz_0[i] * pa_y[i] - ta2_xy_yzz_xxz_1[i] * pc_y[i];

        ta2_xy_yyzz_xyy_0[i] = ta2_xy_yy_xyy_0[i] * fe_0 - ta2_xy_yy_xyy_1[i] * fe_0 + ta2_xy_yyz_xyy_0[i] * pa_z[i] - ta2_xy_yyz_xyy_1[i] * pc_z[i];

        ta2_xy_yyzz_xyz_0[i] = ta2_xy_yy_xyz_0[i] * fe_0 - ta2_xy_yy_xyz_1[i] * fe_0 + ta2_xy_yyz_xy_0[i] * fe_0 - ta2_xy_yyz_xy_1[i] * fe_0 + ta2_xy_yyz_xyz_0[i] * pa_z[i] - ta2_xy_yyz_xyz_1[i] * pc_z[i];

        ta2_xy_yyzz_xzz_0[i] = ta2_xy_zz_xzz_0[i] * fe_0 - ta2_xy_zz_xzz_1[i] * fe_0 + ta1_x_yzz_xzz_1[i] + ta2_xy_yzz_xzz_0[i] * pa_y[i] - ta2_xy_yzz_xzz_1[i] * pc_y[i];

        ta2_xy_yyzz_yyy_0[i] = ta2_xy_yy_yyy_0[i] * fe_0 - ta2_xy_yy_yyy_1[i] * fe_0 + ta2_xy_yyz_yyy_0[i] * pa_z[i] - ta2_xy_yyz_yyy_1[i] * pc_z[i];

        ta2_xy_yyzz_yyz_0[i] = ta2_xy_yy_yyz_0[i] * fe_0 - ta2_xy_yy_yyz_1[i] * fe_0 + ta2_xy_yyz_yy_0[i] * fe_0 - ta2_xy_yyz_yy_1[i] * fe_0 + ta2_xy_yyz_yyz_0[i] * pa_z[i] - ta2_xy_yyz_yyz_1[i] * pc_z[i];

        ta2_xy_yyzz_yzz_0[i] = ta2_xy_yy_yzz_0[i] * fe_0 - ta2_xy_yy_yzz_1[i] * fe_0 + 2.0 * ta2_xy_yyz_yz_0[i] * fe_0 - 2.0 * ta2_xy_yyz_yz_1[i] * fe_0 + ta2_xy_yyz_yzz_0[i] * pa_z[i] - ta2_xy_yyz_yzz_1[i] * pc_z[i];

        ta2_xy_yyzz_zzz_0[i] = ta2_xy_zz_zzz_0[i] * fe_0 - ta2_xy_zz_zzz_1[i] * fe_0 + ta1_x_yzz_zzz_1[i] + ta2_xy_yzz_zzz_0[i] * pa_y[i] - ta2_xy_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 280-290 components of targeted buffer : GF

    auto ta2_xy_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 280);

    auto ta2_xy_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 281);

    auto ta2_xy_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 282);

    auto ta2_xy_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 283);

    auto ta2_xy_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 284);

    auto ta2_xy_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 285);

    auto ta2_xy_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 286);

    auto ta2_xy_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 287);

    auto ta2_xy_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 288);

    auto ta2_xy_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 289);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_zzz_xxx_1, ta1_x_zzz_xxz_1, ta1_x_zzz_xyz_1, ta1_x_zzz_xzz_1, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_1, ta2_xy_yz_xxy_0, ta2_xy_yz_xxy_1, ta2_xy_yz_xyy_0, ta2_xy_yz_xyy_1, ta2_xy_yz_yyy_0, ta2_xy_yz_yyy_1, ta2_xy_yzz_xxy_0, ta2_xy_yzz_xxy_1, ta2_xy_yzz_xyy_0, ta2_xy_yzz_xyy_1, ta2_xy_yzz_yyy_0, ta2_xy_yzz_yyy_1, ta2_xy_yzzz_xxx_0, ta2_xy_yzzz_xxy_0, ta2_xy_yzzz_xxz_0, ta2_xy_yzzz_xyy_0, ta2_xy_yzzz_xyz_0, ta2_xy_yzzz_xzz_0, ta2_xy_yzzz_yyy_0, ta2_xy_yzzz_yyz_0, ta2_xy_yzzz_yzz_0, ta2_xy_yzzz_zzz_0, ta2_xy_zzz_xxx_0, ta2_xy_zzz_xxx_1, ta2_xy_zzz_xxz_0, ta2_xy_zzz_xxz_1, ta2_xy_zzz_xyz_0, ta2_xy_zzz_xyz_1, ta2_xy_zzz_xz_0, ta2_xy_zzz_xz_1, ta2_xy_zzz_xzz_0, ta2_xy_zzz_xzz_1, ta2_xy_zzz_yyz_0, ta2_xy_zzz_yyz_1, ta2_xy_zzz_yz_0, ta2_xy_zzz_yz_1, ta2_xy_zzz_yzz_0, ta2_xy_zzz_yzz_1, ta2_xy_zzz_zz_0, ta2_xy_zzz_zz_1, ta2_xy_zzz_zzz_0, ta2_xy_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzzz_xxx_0[i] = ta1_x_zzz_xxx_1[i] + ta2_xy_zzz_xxx_0[i] * pa_y[i] - ta2_xy_zzz_xxx_1[i] * pc_y[i];

        ta2_xy_yzzz_xxy_0[i] = 2.0 * ta2_xy_yz_xxy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xxy_1[i] * fe_0 + ta2_xy_yzz_xxy_0[i] * pa_z[i] - ta2_xy_yzz_xxy_1[i] * pc_z[i];

        ta2_xy_yzzz_xxz_0[i] = ta1_x_zzz_xxz_1[i] + ta2_xy_zzz_xxz_0[i] * pa_y[i] - ta2_xy_zzz_xxz_1[i] * pc_y[i];

        ta2_xy_yzzz_xyy_0[i] = 2.0 * ta2_xy_yz_xyy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xyy_1[i] * fe_0 + ta2_xy_yzz_xyy_0[i] * pa_z[i] - ta2_xy_yzz_xyy_1[i] * pc_z[i];

        ta2_xy_yzzz_xyz_0[i] = ta2_xy_zzz_xz_0[i] * fe_0 - ta2_xy_zzz_xz_1[i] * fe_0 + ta1_x_zzz_xyz_1[i] + ta2_xy_zzz_xyz_0[i] * pa_y[i] - ta2_xy_zzz_xyz_1[i] * pc_y[i];

        ta2_xy_yzzz_xzz_0[i] = ta1_x_zzz_xzz_1[i] + ta2_xy_zzz_xzz_0[i] * pa_y[i] - ta2_xy_zzz_xzz_1[i] * pc_y[i];

        ta2_xy_yzzz_yyy_0[i] = 2.0 * ta2_xy_yz_yyy_0[i] * fe_0 - 2.0 * ta2_xy_yz_yyy_1[i] * fe_0 + ta2_xy_yzz_yyy_0[i] * pa_z[i] - ta2_xy_yzz_yyy_1[i] * pc_z[i];

        ta2_xy_yzzz_yyz_0[i] = 2.0 * ta2_xy_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_yz_1[i] * fe_0 + ta1_x_zzz_yyz_1[i] + ta2_xy_zzz_yyz_0[i] * pa_y[i] - ta2_xy_zzz_yyz_1[i] * pc_y[i];

        ta2_xy_yzzz_yzz_0[i] = ta2_xy_zzz_zz_0[i] * fe_0 - ta2_xy_zzz_zz_1[i] * fe_0 + ta1_x_zzz_yzz_1[i] + ta2_xy_zzz_yzz_0[i] * pa_y[i] - ta2_xy_zzz_yzz_1[i] * pc_y[i];

        ta2_xy_yzzz_zzz_0[i] = ta1_x_zzz_zzz_1[i] + ta2_xy_zzz_zzz_0[i] * pa_y[i] - ta2_xy_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 290-300 components of targeted buffer : GF

    auto ta2_xy_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 290);

    auto ta2_xy_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 291);

    auto ta2_xy_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 292);

    auto ta2_xy_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 293);

    auto ta2_xy_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 294);

    auto ta2_xy_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 295);

    auto ta2_xy_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 296);

    auto ta2_xy_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 297);

    auto ta2_xy_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 298);

    auto ta2_xy_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 299);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_zz_xxx_0, ta2_xy_zz_xxx_1, ta2_xy_zz_xxy_0, ta2_xy_zz_xxy_1, ta2_xy_zz_xxz_0, ta2_xy_zz_xxz_1, ta2_xy_zz_xyy_0, ta2_xy_zz_xyy_1, ta2_xy_zz_xyz_0, ta2_xy_zz_xyz_1, ta2_xy_zz_xzz_0, ta2_xy_zz_xzz_1, ta2_xy_zz_yyy_0, ta2_xy_zz_yyy_1, ta2_xy_zz_yyz_0, ta2_xy_zz_yyz_1, ta2_xy_zz_yzz_0, ta2_xy_zz_yzz_1, ta2_xy_zz_zzz_0, ta2_xy_zz_zzz_1, ta2_xy_zzz_xx_0, ta2_xy_zzz_xx_1, ta2_xy_zzz_xxx_0, ta2_xy_zzz_xxx_1, ta2_xy_zzz_xxy_0, ta2_xy_zzz_xxy_1, ta2_xy_zzz_xxz_0, ta2_xy_zzz_xxz_1, ta2_xy_zzz_xy_0, ta2_xy_zzz_xy_1, ta2_xy_zzz_xyy_0, ta2_xy_zzz_xyy_1, ta2_xy_zzz_xyz_0, ta2_xy_zzz_xyz_1, ta2_xy_zzz_xz_0, ta2_xy_zzz_xz_1, ta2_xy_zzz_xzz_0, ta2_xy_zzz_xzz_1, ta2_xy_zzz_yy_0, ta2_xy_zzz_yy_1, ta2_xy_zzz_yyy_0, ta2_xy_zzz_yyy_1, ta2_xy_zzz_yyz_0, ta2_xy_zzz_yyz_1, ta2_xy_zzz_yz_0, ta2_xy_zzz_yz_1, ta2_xy_zzz_yzz_0, ta2_xy_zzz_yzz_1, ta2_xy_zzz_zz_0, ta2_xy_zzz_zz_1, ta2_xy_zzz_zzz_0, ta2_xy_zzz_zzz_1, ta2_xy_zzzz_xxx_0, ta2_xy_zzzz_xxy_0, ta2_xy_zzzz_xxz_0, ta2_xy_zzzz_xyy_0, ta2_xy_zzzz_xyz_0, ta2_xy_zzzz_xzz_0, ta2_xy_zzzz_yyy_0, ta2_xy_zzzz_yyz_0, ta2_xy_zzzz_yzz_0, ta2_xy_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzzz_xxx_0[i] = 3.0 * ta2_xy_zz_xxx_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxx_1[i] * fe_0 + ta2_xy_zzz_xxx_0[i] * pa_z[i] - ta2_xy_zzz_xxx_1[i] * pc_z[i];

        ta2_xy_zzzz_xxy_0[i] = 3.0 * ta2_xy_zz_xxy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxy_1[i] * fe_0 + ta2_xy_zzz_xxy_0[i] * pa_z[i] - ta2_xy_zzz_xxy_1[i] * pc_z[i];

        ta2_xy_zzzz_xxz_0[i] = 3.0 * ta2_xy_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxz_1[i] * fe_0 + ta2_xy_zzz_xx_0[i] * fe_0 - ta2_xy_zzz_xx_1[i] * fe_0 + ta2_xy_zzz_xxz_0[i] * pa_z[i] - ta2_xy_zzz_xxz_1[i] * pc_z[i];

        ta2_xy_zzzz_xyy_0[i] = 3.0 * ta2_xy_zz_xyy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xyy_1[i] * fe_0 + ta2_xy_zzz_xyy_0[i] * pa_z[i] - ta2_xy_zzz_xyy_1[i] * pc_z[i];

        ta2_xy_zzzz_xyz_0[i] = 3.0 * ta2_xy_zz_xyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xyz_1[i] * fe_0 + ta2_xy_zzz_xy_0[i] * fe_0 - ta2_xy_zzz_xy_1[i] * fe_0 + ta2_xy_zzz_xyz_0[i] * pa_z[i] - ta2_xy_zzz_xyz_1[i] * pc_z[i];

        ta2_xy_zzzz_xzz_0[i] = 3.0 * ta2_xy_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xzz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_xz_1[i] * fe_0 + ta2_xy_zzz_xzz_0[i] * pa_z[i] - ta2_xy_zzz_xzz_1[i] * pc_z[i];

        ta2_xy_zzzz_yyy_0[i] = 3.0 * ta2_xy_zz_yyy_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyy_1[i] * fe_0 + ta2_xy_zzz_yyy_0[i] * pa_z[i] - ta2_xy_zzz_yyy_1[i] * pc_z[i];

        ta2_xy_zzzz_yyz_0[i] = 3.0 * ta2_xy_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyz_1[i] * fe_0 + ta2_xy_zzz_yy_0[i] * fe_0 - ta2_xy_zzz_yy_1[i] * fe_0 + ta2_xy_zzz_yyz_0[i] * pa_z[i] - ta2_xy_zzz_yyz_1[i] * pc_z[i];

        ta2_xy_zzzz_yzz_0[i] = 3.0 * ta2_xy_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yzz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_yz_1[i] * fe_0 + ta2_xy_zzz_yzz_0[i] * pa_z[i] - ta2_xy_zzz_yzz_1[i] * pc_z[i];

        ta2_xy_zzzz_zzz_0[i] = 3.0 * ta2_xy_zz_zzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_zzz_1[i] * fe_0 + 3.0 * ta2_xy_zzz_zz_0[i] * fe_0 - 3.0 * ta2_xy_zzz_zz_1[i] * fe_0 + ta2_xy_zzz_zzz_0[i] * pa_z[i] - ta2_xy_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 300-310 components of targeted buffer : GF

    auto ta2_xz_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 300);

    auto ta2_xz_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 301);

    auto ta2_xz_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 302);

    auto ta2_xz_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 303);

    auto ta2_xz_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 304);

    auto ta2_xz_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 305);

    auto ta2_xz_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 306);

    auto ta2_xz_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 307);

    auto ta2_xz_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 308);

    auto ta2_xz_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 309);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xxx_xxx_1, ta1_z_xxx_xxy_1, ta1_z_xxx_xxz_1, ta1_z_xxx_xyy_1, ta1_z_xxx_xyz_1, ta1_z_xxx_xzz_1, ta1_z_xxx_yyy_1, ta1_z_xxx_yyz_1, ta1_z_xxx_yzz_1, ta1_z_xxx_zzz_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_yyy_0, ta2_xz_xx_yyy_1, ta2_xz_xx_yyz_0, ta2_xz_xx_yyz_1, ta2_xz_xx_yzz_0, ta2_xz_xx_yzz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xxx_xx_0, ta2_xz_xxx_xx_1, ta2_xz_xxx_xxx_0, ta2_xz_xxx_xxx_1, ta2_xz_xxx_xxy_0, ta2_xz_xxx_xxy_1, ta2_xz_xxx_xxz_0, ta2_xz_xxx_xxz_1, ta2_xz_xxx_xy_0, ta2_xz_xxx_xy_1, ta2_xz_xxx_xyy_0, ta2_xz_xxx_xyy_1, ta2_xz_xxx_xyz_0, ta2_xz_xxx_xyz_1, ta2_xz_xxx_xz_0, ta2_xz_xxx_xz_1, ta2_xz_xxx_xzz_0, ta2_xz_xxx_xzz_1, ta2_xz_xxx_yy_0, ta2_xz_xxx_yy_1, ta2_xz_xxx_yyy_0, ta2_xz_xxx_yyy_1, ta2_xz_xxx_yyz_0, ta2_xz_xxx_yyz_1, ta2_xz_xxx_yz_0, ta2_xz_xxx_yz_1, ta2_xz_xxx_yzz_0, ta2_xz_xxx_yzz_1, ta2_xz_xxx_zz_0, ta2_xz_xxx_zz_1, ta2_xz_xxx_zzz_0, ta2_xz_xxx_zzz_1, ta2_xz_xxxx_xxx_0, ta2_xz_xxxx_xxy_0, ta2_xz_xxxx_xxz_0, ta2_xz_xxxx_xyy_0, ta2_xz_xxxx_xyz_0, ta2_xz_xxxx_xzz_0, ta2_xz_xxxx_yyy_0, ta2_xz_xxxx_yyz_0, ta2_xz_xxxx_yzz_0, ta2_xz_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxx_xxx_0[i] = 3.0 * ta2_xz_xx_xxx_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxx_1[i] * fe_0 + 3.0 * ta2_xz_xxx_xx_0[i] * fe_0 - 3.0 * ta2_xz_xxx_xx_1[i] * fe_0 + ta1_z_xxx_xxx_1[i] + ta2_xz_xxx_xxx_0[i] * pa_x[i] - ta2_xz_xxx_xxx_1[i] * pc_x[i];

        ta2_xz_xxxx_xxy_0[i] = 3.0 * ta2_xz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxy_1[i] * fe_0 + 2.0 * ta2_xz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xy_1[i] * fe_0 + ta1_z_xxx_xxy_1[i] + ta2_xz_xxx_xxy_0[i] * pa_x[i] - ta2_xz_xxx_xxy_1[i] * pc_x[i];

        ta2_xz_xxxx_xxz_0[i] = 3.0 * ta2_xz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxz_1[i] * fe_0 + 2.0 * ta2_xz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xz_1[i] * fe_0 + ta1_z_xxx_xxz_1[i] + ta2_xz_xxx_xxz_0[i] * pa_x[i] - ta2_xz_xxx_xxz_1[i] * pc_x[i];

        ta2_xz_xxxx_xyy_0[i] = 3.0 * ta2_xz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyy_1[i] * fe_0 + ta2_xz_xxx_yy_0[i] * fe_0 - ta2_xz_xxx_yy_1[i] * fe_0 + ta1_z_xxx_xyy_1[i] + ta2_xz_xxx_xyy_0[i] * pa_x[i] - ta2_xz_xxx_xyy_1[i] * pc_x[i];

        ta2_xz_xxxx_xyz_0[i] = 3.0 * ta2_xz_xx_xyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyz_1[i] * fe_0 + ta2_xz_xxx_yz_0[i] * fe_0 - ta2_xz_xxx_yz_1[i] * fe_0 + ta1_z_xxx_xyz_1[i] + ta2_xz_xxx_xyz_0[i] * pa_x[i] - ta2_xz_xxx_xyz_1[i] * pc_x[i];

        ta2_xz_xxxx_xzz_0[i] = 3.0 * ta2_xz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xzz_1[i] * fe_0 + ta2_xz_xxx_zz_0[i] * fe_0 - ta2_xz_xxx_zz_1[i] * fe_0 + ta1_z_xxx_xzz_1[i] + ta2_xz_xxx_xzz_0[i] * pa_x[i] - ta2_xz_xxx_xzz_1[i] * pc_x[i];

        ta2_xz_xxxx_yyy_0[i] = 3.0 * ta2_xz_xx_yyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyy_1[i] * fe_0 + ta1_z_xxx_yyy_1[i] + ta2_xz_xxx_yyy_0[i] * pa_x[i] - ta2_xz_xxx_yyy_1[i] * pc_x[i];

        ta2_xz_xxxx_yyz_0[i] = 3.0 * ta2_xz_xx_yyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyz_1[i] * fe_0 + ta1_z_xxx_yyz_1[i] + ta2_xz_xxx_yyz_0[i] * pa_x[i] - ta2_xz_xxx_yyz_1[i] * pc_x[i];

        ta2_xz_xxxx_yzz_0[i] = 3.0 * ta2_xz_xx_yzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yzz_1[i] * fe_0 + ta1_z_xxx_yzz_1[i] + ta2_xz_xxx_yzz_0[i] * pa_x[i] - ta2_xz_xxx_yzz_1[i] * pc_x[i];

        ta2_xz_xxxx_zzz_0[i] = 3.0 * ta2_xz_xx_zzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_zzz_1[i] * fe_0 + ta1_z_xxx_zzz_1[i] + ta2_xz_xxx_zzz_0[i] * pa_x[i] - ta2_xz_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : GF

    auto ta2_xz_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 310);

    auto ta2_xz_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 311);

    auto ta2_xz_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 312);

    auto ta2_xz_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 313);

    auto ta2_xz_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 314);

    auto ta2_xz_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 315);

    auto ta2_xz_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 316);

    auto ta2_xz_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 317);

    auto ta2_xz_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 318);

    auto ta2_xz_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 319);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_xxx_xx_0, ta2_xz_xxx_xx_1, ta2_xz_xxx_xxx_0, ta2_xz_xxx_xxx_1, ta2_xz_xxx_xxy_0, ta2_xz_xxx_xxy_1, ta2_xz_xxx_xxz_0, ta2_xz_xxx_xxz_1, ta2_xz_xxx_xy_0, ta2_xz_xxx_xy_1, ta2_xz_xxx_xyy_0, ta2_xz_xxx_xyy_1, ta2_xz_xxx_xyz_0, ta2_xz_xxx_xyz_1, ta2_xz_xxx_xz_0, ta2_xz_xxx_xz_1, ta2_xz_xxx_xzz_0, ta2_xz_xxx_xzz_1, ta2_xz_xxx_yy_0, ta2_xz_xxx_yy_1, ta2_xz_xxx_yyy_0, ta2_xz_xxx_yyy_1, ta2_xz_xxx_yyz_0, ta2_xz_xxx_yyz_1, ta2_xz_xxx_yz_0, ta2_xz_xxx_yz_1, ta2_xz_xxx_yzz_0, ta2_xz_xxx_yzz_1, ta2_xz_xxx_zz_0, ta2_xz_xxx_zz_1, ta2_xz_xxx_zzz_0, ta2_xz_xxx_zzz_1, ta2_xz_xxxy_xxx_0, ta2_xz_xxxy_xxy_0, ta2_xz_xxxy_xxz_0, ta2_xz_xxxy_xyy_0, ta2_xz_xxxy_xyz_0, ta2_xz_xxxy_xzz_0, ta2_xz_xxxy_yyy_0, ta2_xz_xxxy_yyz_0, ta2_xz_xxxy_yzz_0, ta2_xz_xxxy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxy_xxx_0[i] = ta2_xz_xxx_xxx_0[i] * pa_y[i] - ta2_xz_xxx_xxx_1[i] * pc_y[i];

        ta2_xz_xxxy_xxy_0[i] = ta2_xz_xxx_xx_0[i] * fe_0 - ta2_xz_xxx_xx_1[i] * fe_0 + ta2_xz_xxx_xxy_0[i] * pa_y[i] - ta2_xz_xxx_xxy_1[i] * pc_y[i];

        ta2_xz_xxxy_xxz_0[i] = ta2_xz_xxx_xxz_0[i] * pa_y[i] - ta2_xz_xxx_xxz_1[i] * pc_y[i];

        ta2_xz_xxxy_xyy_0[i] = 2.0 * ta2_xz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xy_1[i] * fe_0 + ta2_xz_xxx_xyy_0[i] * pa_y[i] - ta2_xz_xxx_xyy_1[i] * pc_y[i];

        ta2_xz_xxxy_xyz_0[i] = ta2_xz_xxx_xz_0[i] * fe_0 - ta2_xz_xxx_xz_1[i] * fe_0 + ta2_xz_xxx_xyz_0[i] * pa_y[i] - ta2_xz_xxx_xyz_1[i] * pc_y[i];

        ta2_xz_xxxy_xzz_0[i] = ta2_xz_xxx_xzz_0[i] * pa_y[i] - ta2_xz_xxx_xzz_1[i] * pc_y[i];

        ta2_xz_xxxy_yyy_0[i] = 3.0 * ta2_xz_xxx_yy_0[i] * fe_0 - 3.0 * ta2_xz_xxx_yy_1[i] * fe_0 + ta2_xz_xxx_yyy_0[i] * pa_y[i] - ta2_xz_xxx_yyy_1[i] * pc_y[i];

        ta2_xz_xxxy_yyz_0[i] = 2.0 * ta2_xz_xxx_yz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_yz_1[i] * fe_0 + ta2_xz_xxx_yyz_0[i] * pa_y[i] - ta2_xz_xxx_yyz_1[i] * pc_y[i];

        ta2_xz_xxxy_yzz_0[i] = ta2_xz_xxx_zz_0[i] * fe_0 - ta2_xz_xxx_zz_1[i] * fe_0 + ta2_xz_xxx_yzz_0[i] * pa_y[i] - ta2_xz_xxx_yzz_1[i] * pc_y[i];

        ta2_xz_xxxy_zzz_0[i] = ta2_xz_xxx_zzz_0[i] * pa_y[i] - ta2_xz_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 320-330 components of targeted buffer : GF

    auto ta2_xz_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 320);

    auto ta2_xz_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 321);

    auto ta2_xz_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 322);

    auto ta2_xz_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 323);

    auto ta2_xz_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 324);

    auto ta2_xz_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 325);

    auto ta2_xz_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 326);

    auto ta2_xz_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 327);

    auto ta2_xz_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 328);

    auto ta2_xz_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 329);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_1, ta1_x_xxx_yyy_1, ta1_z_xxz_yyz_1, ta1_z_xxz_yzz_1, ta1_z_xxz_zzz_1, ta2_xz_xxx_xx_0, ta2_xz_xxx_xx_1, ta2_xz_xxx_xxx_0, ta2_xz_xxx_xxx_1, ta2_xz_xxx_xxy_0, ta2_xz_xxx_xxy_1, ta2_xz_xxx_xxz_0, ta2_xz_xxx_xxz_1, ta2_xz_xxx_xy_0, ta2_xz_xxx_xy_1, ta2_xz_xxx_xyy_0, ta2_xz_xxx_xyy_1, ta2_xz_xxx_xyz_0, ta2_xz_xxx_xyz_1, ta2_xz_xxx_xz_0, ta2_xz_xxx_xz_1, ta2_xz_xxx_xzz_0, ta2_xz_xxx_xzz_1, ta2_xz_xxx_yyy_0, ta2_xz_xxx_yyy_1, ta2_xz_xxxz_xxx_0, ta2_xz_xxxz_xxy_0, ta2_xz_xxxz_xxz_0, ta2_xz_xxxz_xyy_0, ta2_xz_xxxz_xyz_0, ta2_xz_xxxz_xzz_0, ta2_xz_xxxz_yyy_0, ta2_xz_xxxz_yyz_0, ta2_xz_xxxz_yzz_0, ta2_xz_xxxz_zzz_0, ta2_xz_xxz_yyz_0, ta2_xz_xxz_yyz_1, ta2_xz_xxz_yzz_0, ta2_xz_xxz_yzz_1, ta2_xz_xxz_zzz_0, ta2_xz_xxz_zzz_1, ta2_xz_xz_yyz_0, ta2_xz_xz_yyz_1, ta2_xz_xz_yzz_0, ta2_xz_xz_yzz_1, ta2_xz_xz_zzz_0, ta2_xz_xz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxz_xxx_0[i] = ta1_x_xxx_xxx_1[i] + ta2_xz_xxx_xxx_0[i] * pa_z[i] - ta2_xz_xxx_xxx_1[i] * pc_z[i];

        ta2_xz_xxxz_xxy_0[i] = ta1_x_xxx_xxy_1[i] + ta2_xz_xxx_xxy_0[i] * pa_z[i] - ta2_xz_xxx_xxy_1[i] * pc_z[i];

        ta2_xz_xxxz_xxz_0[i] = ta2_xz_xxx_xx_0[i] * fe_0 - ta2_xz_xxx_xx_1[i] * fe_0 + ta1_x_xxx_xxz_1[i] + ta2_xz_xxx_xxz_0[i] * pa_z[i] - ta2_xz_xxx_xxz_1[i] * pc_z[i];

        ta2_xz_xxxz_xyy_0[i] = ta1_x_xxx_xyy_1[i] + ta2_xz_xxx_xyy_0[i] * pa_z[i] - ta2_xz_xxx_xyy_1[i] * pc_z[i];

        ta2_xz_xxxz_xyz_0[i] = ta2_xz_xxx_xy_0[i] * fe_0 - ta2_xz_xxx_xy_1[i] * fe_0 + ta1_x_xxx_xyz_1[i] + ta2_xz_xxx_xyz_0[i] * pa_z[i] - ta2_xz_xxx_xyz_1[i] * pc_z[i];

        ta2_xz_xxxz_xzz_0[i] = 2.0 * ta2_xz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xz_1[i] * fe_0 + ta1_x_xxx_xzz_1[i] + ta2_xz_xxx_xzz_0[i] * pa_z[i] - ta2_xz_xxx_xzz_1[i] * pc_z[i];

        ta2_xz_xxxz_yyy_0[i] = ta1_x_xxx_yyy_1[i] + ta2_xz_xxx_yyy_0[i] * pa_z[i] - ta2_xz_xxx_yyy_1[i] * pc_z[i];

        ta2_xz_xxxz_yyz_0[i] = 2.0 * ta2_xz_xz_yyz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yyz_1[i] * fe_0 + ta1_z_xxz_yyz_1[i] + ta2_xz_xxz_yyz_0[i] * pa_x[i] - ta2_xz_xxz_yyz_1[i] * pc_x[i];

        ta2_xz_xxxz_yzz_0[i] = 2.0 * ta2_xz_xz_yzz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yzz_1[i] * fe_0 + ta1_z_xxz_yzz_1[i] + ta2_xz_xxz_yzz_0[i] * pa_x[i] - ta2_xz_xxz_yzz_1[i] * pc_x[i];

        ta2_xz_xxxz_zzz_0[i] = 2.0 * ta2_xz_xz_zzz_0[i] * fe_0 - 2.0 * ta2_xz_xz_zzz_1[i] * fe_0 + ta1_z_xxz_zzz_1[i] + ta2_xz_xxz_zzz_0[i] * pa_x[i] - ta2_xz_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 330-340 components of targeted buffer : GF

    auto ta2_xz_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 330);

    auto ta2_xz_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 331);

    auto ta2_xz_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 332);

    auto ta2_xz_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 333);

    auto ta2_xz_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 334);

    auto ta2_xz_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 335);

    auto ta2_xz_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 336);

    auto ta2_xz_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 337);

    auto ta2_xz_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 338);

    auto ta2_xz_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 339);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xyy_yyy_1, ta1_z_xyy_yyz_1, ta1_z_xyy_yzz_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xxz_0, ta2_xz_xx_xxz_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xx_xyz_0, ta2_xz_xx_xyz_1, ta2_xz_xx_xzz_0, ta2_xz_xx_xzz_1, ta2_xz_xx_zzz_0, ta2_xz_xx_zzz_1, ta2_xz_xxy_xx_0, ta2_xz_xxy_xx_1, ta2_xz_xxy_xxx_0, ta2_xz_xxy_xxx_1, ta2_xz_xxy_xxy_0, ta2_xz_xxy_xxy_1, ta2_xz_xxy_xxz_0, ta2_xz_xxy_xxz_1, ta2_xz_xxy_xy_0, ta2_xz_xxy_xy_1, ta2_xz_xxy_xyy_0, ta2_xz_xxy_xyy_1, ta2_xz_xxy_xyz_0, ta2_xz_xxy_xyz_1, ta2_xz_xxy_xz_0, ta2_xz_xxy_xz_1, ta2_xz_xxy_xzz_0, ta2_xz_xxy_xzz_1, ta2_xz_xxy_zzz_0, ta2_xz_xxy_zzz_1, ta2_xz_xxyy_xxx_0, ta2_xz_xxyy_xxy_0, ta2_xz_xxyy_xxz_0, ta2_xz_xxyy_xyy_0, ta2_xz_xxyy_xyz_0, ta2_xz_xxyy_xzz_0, ta2_xz_xxyy_yyy_0, ta2_xz_xxyy_yyz_0, ta2_xz_xxyy_yzz_0, ta2_xz_xxyy_zzz_0, ta2_xz_xyy_yyy_0, ta2_xz_xyy_yyy_1, ta2_xz_xyy_yyz_0, ta2_xz_xyy_yyz_1, ta2_xz_xyy_yzz_0, ta2_xz_xyy_yzz_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyy_xxx_0[i] = ta2_xz_xx_xxx_0[i] * fe_0 - ta2_xz_xx_xxx_1[i] * fe_0 + ta2_xz_xxy_xxx_0[i] * pa_y[i] - ta2_xz_xxy_xxx_1[i] * pc_y[i];

        ta2_xz_xxyy_xxy_0[i] = ta2_xz_xx_xxy_0[i] * fe_0 - ta2_xz_xx_xxy_1[i] * fe_0 + ta2_xz_xxy_xx_0[i] * fe_0 - ta2_xz_xxy_xx_1[i] * fe_0 + ta2_xz_xxy_xxy_0[i] * pa_y[i] - ta2_xz_xxy_xxy_1[i] * pc_y[i];

        ta2_xz_xxyy_xxz_0[i] = ta2_xz_xx_xxz_0[i] * fe_0 - ta2_xz_xx_xxz_1[i] * fe_0 + ta2_xz_xxy_xxz_0[i] * pa_y[i] - ta2_xz_xxy_xxz_1[i] * pc_y[i];

        ta2_xz_xxyy_xyy_0[i] = ta2_xz_xx_xyy_0[i] * fe_0 - ta2_xz_xx_xyy_1[i] * fe_0 + 2.0 * ta2_xz_xxy_xy_0[i] * fe_0 - 2.0 * ta2_xz_xxy_xy_1[i] * fe_0 + ta2_xz_xxy_xyy_0[i] * pa_y[i] - ta2_xz_xxy_xyy_1[i] * pc_y[i];

        ta2_xz_xxyy_xyz_0[i] = ta2_xz_xx_xyz_0[i] * fe_0 - ta2_xz_xx_xyz_1[i] * fe_0 + ta2_xz_xxy_xz_0[i] * fe_0 - ta2_xz_xxy_xz_1[i] * fe_0 + ta2_xz_xxy_xyz_0[i] * pa_y[i] - ta2_xz_xxy_xyz_1[i] * pc_y[i];

        ta2_xz_xxyy_xzz_0[i] = ta2_xz_xx_xzz_0[i] * fe_0 - ta2_xz_xx_xzz_1[i] * fe_0 + ta2_xz_xxy_xzz_0[i] * pa_y[i] - ta2_xz_xxy_xzz_1[i] * pc_y[i];

        ta2_xz_xxyy_yyy_0[i] = ta2_xz_yy_yyy_0[i] * fe_0 - ta2_xz_yy_yyy_1[i] * fe_0 + ta1_z_xyy_yyy_1[i] + ta2_xz_xyy_yyy_0[i] * pa_x[i] - ta2_xz_xyy_yyy_1[i] * pc_x[i];

        ta2_xz_xxyy_yyz_0[i] = ta2_xz_yy_yyz_0[i] * fe_0 - ta2_xz_yy_yyz_1[i] * fe_0 + ta1_z_xyy_yyz_1[i] + ta2_xz_xyy_yyz_0[i] * pa_x[i] - ta2_xz_xyy_yyz_1[i] * pc_x[i];

        ta2_xz_xxyy_yzz_0[i] = ta2_xz_yy_yzz_0[i] * fe_0 - ta2_xz_yy_yzz_1[i] * fe_0 + ta1_z_xyy_yzz_1[i] + ta2_xz_xyy_yzz_0[i] * pa_x[i] - ta2_xz_xyy_yzz_1[i] * pc_x[i];

        ta2_xz_xxyy_zzz_0[i] = ta2_xz_xx_zzz_0[i] * fe_0 - ta2_xz_xx_zzz_1[i] * fe_0 + ta2_xz_xxy_zzz_0[i] * pa_y[i] - ta2_xz_xxy_zzz_1[i] * pc_y[i];
    }

    // Set up 340-350 components of targeted buffer : GF

    auto ta2_xz_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 340);

    auto ta2_xz_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 341);

    auto ta2_xz_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 342);

    auto ta2_xz_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 343);

    auto ta2_xz_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 344);

    auto ta2_xz_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 345);

    auto ta2_xz_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 346);

    auto ta2_xz_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 347);

    auto ta2_xz_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 348);

    auto ta2_xz_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 349);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxy_xxy_1, ta1_x_xxy_xyy_1, ta1_x_xxy_yyy_1, ta2_xz_xxy_xxy_0, ta2_xz_xxy_xxy_1, ta2_xz_xxy_xyy_0, ta2_xz_xxy_xyy_1, ta2_xz_xxy_yyy_0, ta2_xz_xxy_yyy_1, ta2_xz_xxyz_xxx_0, ta2_xz_xxyz_xxy_0, ta2_xz_xxyz_xxz_0, ta2_xz_xxyz_xyy_0, ta2_xz_xxyz_xyz_0, ta2_xz_xxyz_xzz_0, ta2_xz_xxyz_yyy_0, ta2_xz_xxyz_yyz_0, ta2_xz_xxyz_yzz_0, ta2_xz_xxyz_zzz_0, ta2_xz_xxz_xxx_0, ta2_xz_xxz_xxx_1, ta2_xz_xxz_xxz_0, ta2_xz_xxz_xxz_1, ta2_xz_xxz_xyz_0, ta2_xz_xxz_xyz_1, ta2_xz_xxz_xz_0, ta2_xz_xxz_xz_1, ta2_xz_xxz_xzz_0, ta2_xz_xxz_xzz_1, ta2_xz_xxz_yyz_0, ta2_xz_xxz_yyz_1, ta2_xz_xxz_yz_0, ta2_xz_xxz_yz_1, ta2_xz_xxz_yzz_0, ta2_xz_xxz_yzz_1, ta2_xz_xxz_zz_0, ta2_xz_xxz_zz_1, ta2_xz_xxz_zzz_0, ta2_xz_xxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyz_xxx_0[i] = ta2_xz_xxz_xxx_0[i] * pa_y[i] - ta2_xz_xxz_xxx_1[i] * pc_y[i];

        ta2_xz_xxyz_xxy_0[i] = ta1_x_xxy_xxy_1[i] + ta2_xz_xxy_xxy_0[i] * pa_z[i] - ta2_xz_xxy_xxy_1[i] * pc_z[i];

        ta2_xz_xxyz_xxz_0[i] = ta2_xz_xxz_xxz_0[i] * pa_y[i] - ta2_xz_xxz_xxz_1[i] * pc_y[i];

        ta2_xz_xxyz_xyy_0[i] = ta1_x_xxy_xyy_1[i] + ta2_xz_xxy_xyy_0[i] * pa_z[i] - ta2_xz_xxy_xyy_1[i] * pc_z[i];

        ta2_xz_xxyz_xyz_0[i] = ta2_xz_xxz_xz_0[i] * fe_0 - ta2_xz_xxz_xz_1[i] * fe_0 + ta2_xz_xxz_xyz_0[i] * pa_y[i] - ta2_xz_xxz_xyz_1[i] * pc_y[i];

        ta2_xz_xxyz_xzz_0[i] = ta2_xz_xxz_xzz_0[i] * pa_y[i] - ta2_xz_xxz_xzz_1[i] * pc_y[i];

        ta2_xz_xxyz_yyy_0[i] = ta1_x_xxy_yyy_1[i] + ta2_xz_xxy_yyy_0[i] * pa_z[i] - ta2_xz_xxy_yyy_1[i] * pc_z[i];

        ta2_xz_xxyz_yyz_0[i] = 2.0 * ta2_xz_xxz_yz_0[i] * fe_0 - 2.0 * ta2_xz_xxz_yz_1[i] * fe_0 + ta2_xz_xxz_yyz_0[i] * pa_y[i] - ta2_xz_xxz_yyz_1[i] * pc_y[i];

        ta2_xz_xxyz_yzz_0[i] = ta2_xz_xxz_zz_0[i] * fe_0 - ta2_xz_xxz_zz_1[i] * fe_0 + ta2_xz_xxz_yzz_0[i] * pa_y[i] - ta2_xz_xxz_yzz_1[i] * pc_y[i];

        ta2_xz_xxyz_zzz_0[i] = ta2_xz_xxz_zzz_0[i] * pa_y[i] - ta2_xz_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 350-360 components of targeted buffer : GF

    auto ta2_xz_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 350);

    auto ta2_xz_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 351);

    auto ta2_xz_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 352);

    auto ta2_xz_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 353);

    auto ta2_xz_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 354);

    auto ta2_xz_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 355);

    auto ta2_xz_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 356);

    auto ta2_xz_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 357);

    auto ta2_xz_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 358);

    auto ta2_xz_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 359);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxz_xxx_1, ta1_x_xxz_xxy_1, ta1_x_xxz_xyy_1, ta1_z_xzz_xxz_1, ta1_z_xzz_xyz_1, ta1_z_xzz_xzz_1, ta1_z_xzz_yyy_1, ta1_z_xzz_yyz_1, ta1_z_xzz_yzz_1, ta1_z_xzz_zzz_1, ta2_xz_xx_xxx_0, ta2_xz_xx_xxx_1, ta2_xz_xx_xxy_0, ta2_xz_xx_xxy_1, ta2_xz_xx_xyy_0, ta2_xz_xx_xyy_1, ta2_xz_xxz_xxx_0, ta2_xz_xxz_xxx_1, ta2_xz_xxz_xxy_0, ta2_xz_xxz_xxy_1, ta2_xz_xxz_xyy_0, ta2_xz_xxz_xyy_1, ta2_xz_xxzz_xxx_0, ta2_xz_xxzz_xxy_0, ta2_xz_xxzz_xxz_0, ta2_xz_xxzz_xyy_0, ta2_xz_xxzz_xyz_0, ta2_xz_xxzz_xzz_0, ta2_xz_xxzz_yyy_0, ta2_xz_xxzz_yyz_0, ta2_xz_xxzz_yzz_0, ta2_xz_xxzz_zzz_0, ta2_xz_xzz_xxz_0, ta2_xz_xzz_xxz_1, ta2_xz_xzz_xyz_0, ta2_xz_xzz_xyz_1, ta2_xz_xzz_xz_0, ta2_xz_xzz_xz_1, ta2_xz_xzz_xzz_0, ta2_xz_xzz_xzz_1, ta2_xz_xzz_yyy_0, ta2_xz_xzz_yyy_1, ta2_xz_xzz_yyz_0, ta2_xz_xzz_yyz_1, ta2_xz_xzz_yz_0, ta2_xz_xzz_yz_1, ta2_xz_xzz_yzz_0, ta2_xz_xzz_yzz_1, ta2_xz_xzz_zz_0, ta2_xz_xzz_zz_1, ta2_xz_xzz_zzz_0, ta2_xz_xzz_zzz_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxzz_xxx_0[i] = ta2_xz_xx_xxx_0[i] * fe_0 - ta2_xz_xx_xxx_1[i] * fe_0 + ta1_x_xxz_xxx_1[i] + ta2_xz_xxz_xxx_0[i] * pa_z[i] - ta2_xz_xxz_xxx_1[i] * pc_z[i];

        ta2_xz_xxzz_xxy_0[i] = ta2_xz_xx_xxy_0[i] * fe_0 - ta2_xz_xx_xxy_1[i] * fe_0 + ta1_x_xxz_xxy_1[i] + ta2_xz_xxz_xxy_0[i] * pa_z[i] - ta2_xz_xxz_xxy_1[i] * pc_z[i];

        ta2_xz_xxzz_xxz_0[i] = ta2_xz_zz_xxz_0[i] * fe_0 - ta2_xz_zz_xxz_1[i] * fe_0 + 2.0 * ta2_xz_xzz_xz_0[i] * fe_0 - 2.0 * ta2_xz_xzz_xz_1[i] * fe_0 + ta1_z_xzz_xxz_1[i] + ta2_xz_xzz_xxz_0[i] * pa_x[i] - ta2_xz_xzz_xxz_1[i] * pc_x[i];

        ta2_xz_xxzz_xyy_0[i] = ta2_xz_xx_xyy_0[i] * fe_0 - ta2_xz_xx_xyy_1[i] * fe_0 + ta1_x_xxz_xyy_1[i] + ta2_xz_xxz_xyy_0[i] * pa_z[i] - ta2_xz_xxz_xyy_1[i] * pc_z[i];

        ta2_xz_xxzz_xyz_0[i] = ta2_xz_zz_xyz_0[i] * fe_0 - ta2_xz_zz_xyz_1[i] * fe_0 + ta2_xz_xzz_yz_0[i] * fe_0 - ta2_xz_xzz_yz_1[i] * fe_0 + ta1_z_xzz_xyz_1[i] + ta2_xz_xzz_xyz_0[i] * pa_x[i] - ta2_xz_xzz_xyz_1[i] * pc_x[i];

        ta2_xz_xxzz_xzz_0[i] = ta2_xz_zz_xzz_0[i] * fe_0 - ta2_xz_zz_xzz_1[i] * fe_0 + ta2_xz_xzz_zz_0[i] * fe_0 - ta2_xz_xzz_zz_1[i] * fe_0 + ta1_z_xzz_xzz_1[i] + ta2_xz_xzz_xzz_0[i] * pa_x[i] - ta2_xz_xzz_xzz_1[i] * pc_x[i];

        ta2_xz_xxzz_yyy_0[i] = ta2_xz_zz_yyy_0[i] * fe_0 - ta2_xz_zz_yyy_1[i] * fe_0 + ta1_z_xzz_yyy_1[i] + ta2_xz_xzz_yyy_0[i] * pa_x[i] - ta2_xz_xzz_yyy_1[i] * pc_x[i];

        ta2_xz_xxzz_yyz_0[i] = ta2_xz_zz_yyz_0[i] * fe_0 - ta2_xz_zz_yyz_1[i] * fe_0 + ta1_z_xzz_yyz_1[i] + ta2_xz_xzz_yyz_0[i] * pa_x[i] - ta2_xz_xzz_yyz_1[i] * pc_x[i];

        ta2_xz_xxzz_yzz_0[i] = ta2_xz_zz_yzz_0[i] * fe_0 - ta2_xz_zz_yzz_1[i] * fe_0 + ta1_z_xzz_yzz_1[i] + ta2_xz_xzz_yzz_0[i] * pa_x[i] - ta2_xz_xzz_yzz_1[i] * pc_x[i];

        ta2_xz_xxzz_zzz_0[i] = ta2_xz_zz_zzz_0[i] * fe_0 - ta2_xz_zz_zzz_1[i] * fe_0 + ta1_z_xzz_zzz_1[i] + ta2_xz_xzz_zzz_0[i] * pa_x[i] - ta2_xz_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 360-370 components of targeted buffer : GF

    auto ta2_xz_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 360);

    auto ta2_xz_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 361);

    auto ta2_xz_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 362);

    auto ta2_xz_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 363);

    auto ta2_xz_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 364);

    auto ta2_xz_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 365);

    auto ta2_xz_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 366);

    auto ta2_xz_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 367);

    auto ta2_xz_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 368);

    auto ta2_xz_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 369);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_yyy_xxy_1, ta1_z_yyy_xyy_1, ta1_z_yyy_xyz_1, ta1_z_yyy_yyy_1, ta1_z_yyy_yyz_1, ta1_z_yyy_yzz_1, ta1_z_yyy_zzz_1, ta2_xz_xy_xxx_0, ta2_xz_xy_xxx_1, ta2_xz_xy_xxz_0, ta2_xz_xy_xxz_1, ta2_xz_xy_xzz_0, ta2_xz_xy_xzz_1, ta2_xz_xyy_xxx_0, ta2_xz_xyy_xxx_1, ta2_xz_xyy_xxz_0, ta2_xz_xyy_xxz_1, ta2_xz_xyy_xzz_0, ta2_xz_xyy_xzz_1, ta2_xz_xyyy_xxx_0, ta2_xz_xyyy_xxy_0, ta2_xz_xyyy_xxz_0, ta2_xz_xyyy_xyy_0, ta2_xz_xyyy_xyz_0, ta2_xz_xyyy_xzz_0, ta2_xz_xyyy_yyy_0, ta2_xz_xyyy_yyz_0, ta2_xz_xyyy_yzz_0, ta2_xz_xyyy_zzz_0, ta2_xz_yyy_xxy_0, ta2_xz_yyy_xxy_1, ta2_xz_yyy_xy_0, ta2_xz_yyy_xy_1, ta2_xz_yyy_xyy_0, ta2_xz_yyy_xyy_1, ta2_xz_yyy_xyz_0, ta2_xz_yyy_xyz_1, ta2_xz_yyy_yy_0, ta2_xz_yyy_yy_1, ta2_xz_yyy_yyy_0, ta2_xz_yyy_yyy_1, ta2_xz_yyy_yyz_0, ta2_xz_yyy_yyz_1, ta2_xz_yyy_yz_0, ta2_xz_yyy_yz_1, ta2_xz_yyy_yzz_0, ta2_xz_yyy_yzz_1, ta2_xz_yyy_zzz_0, ta2_xz_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyy_xxx_0[i] = 2.0 * ta2_xz_xy_xxx_0[i] * fe_0 - 2.0 * ta2_xz_xy_xxx_1[i] * fe_0 + ta2_xz_xyy_xxx_0[i] * pa_y[i] - ta2_xz_xyy_xxx_1[i] * pc_y[i];

        ta2_xz_xyyy_xxy_0[i] = 2.0 * ta2_xz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xz_yyy_xy_1[i] * fe_0 + ta1_z_yyy_xxy_1[i] + ta2_xz_yyy_xxy_0[i] * pa_x[i] - ta2_xz_yyy_xxy_1[i] * pc_x[i];

        ta2_xz_xyyy_xxz_0[i] = 2.0 * ta2_xz_xy_xxz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xxz_1[i] * fe_0 + ta2_xz_xyy_xxz_0[i] * pa_y[i] - ta2_xz_xyy_xxz_1[i] * pc_y[i];

        ta2_xz_xyyy_xyy_0[i] = ta2_xz_yyy_yy_0[i] * fe_0 - ta2_xz_yyy_yy_1[i] * fe_0 + ta1_z_yyy_xyy_1[i] + ta2_xz_yyy_xyy_0[i] * pa_x[i] - ta2_xz_yyy_xyy_1[i] * pc_x[i];

        ta2_xz_xyyy_xyz_0[i] = ta2_xz_yyy_yz_0[i] * fe_0 - ta2_xz_yyy_yz_1[i] * fe_0 + ta1_z_yyy_xyz_1[i] + ta2_xz_yyy_xyz_0[i] * pa_x[i] - ta2_xz_yyy_xyz_1[i] * pc_x[i];

        ta2_xz_xyyy_xzz_0[i] = 2.0 * ta2_xz_xy_xzz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xzz_1[i] * fe_0 + ta2_xz_xyy_xzz_0[i] * pa_y[i] - ta2_xz_xyy_xzz_1[i] * pc_y[i];

        ta2_xz_xyyy_yyy_0[i] = ta1_z_yyy_yyy_1[i] + ta2_xz_yyy_yyy_0[i] * pa_x[i] - ta2_xz_yyy_yyy_1[i] * pc_x[i];

        ta2_xz_xyyy_yyz_0[i] = ta1_z_yyy_yyz_1[i] + ta2_xz_yyy_yyz_0[i] * pa_x[i] - ta2_xz_yyy_yyz_1[i] * pc_x[i];

        ta2_xz_xyyy_yzz_0[i] = ta1_z_yyy_yzz_1[i] + ta2_xz_yyy_yzz_0[i] * pa_x[i] - ta2_xz_yyy_yzz_1[i] * pc_x[i];

        ta2_xz_xyyy_zzz_0[i] = ta1_z_yyy_zzz_1[i] + ta2_xz_yyy_zzz_0[i] * pa_x[i] - ta2_xz_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 370-380 components of targeted buffer : GF

    auto ta2_xz_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 370);

    auto ta2_xz_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 371);

    auto ta2_xz_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 372);

    auto ta2_xz_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 373);

    auto ta2_xz_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 374);

    auto ta2_xz_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 375);

    auto ta2_xz_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 376);

    auto ta2_xz_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 377);

    auto ta2_xz_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 378);

    auto ta2_xz_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 379);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyy_xxx_1, ta1_x_xyy_xxy_1, ta1_x_xyy_xyy_1, ta1_z_yyz_xyz_1, ta1_z_yyz_yyy_1, ta1_z_yyz_yyz_1, ta1_z_yyz_yzz_1, ta1_z_yyz_zzz_1, ta2_xz_xyy_xxx_0, ta2_xz_xyy_xxx_1, ta2_xz_xyy_xxy_0, ta2_xz_xyy_xxy_1, ta2_xz_xyy_xyy_0, ta2_xz_xyy_xyy_1, ta2_xz_xyyz_xxx_0, ta2_xz_xyyz_xxy_0, ta2_xz_xyyz_xxz_0, ta2_xz_xyyz_xyy_0, ta2_xz_xyyz_xyz_0, ta2_xz_xyyz_xzz_0, ta2_xz_xyyz_yyy_0, ta2_xz_xyyz_yyz_0, ta2_xz_xyyz_yzz_0, ta2_xz_xyyz_zzz_0, ta2_xz_xyz_xxz_0, ta2_xz_xyz_xxz_1, ta2_xz_xyz_xzz_0, ta2_xz_xyz_xzz_1, ta2_xz_xz_xxz_0, ta2_xz_xz_xxz_1, ta2_xz_xz_xzz_0, ta2_xz_xz_xzz_1, ta2_xz_yyz_xyz_0, ta2_xz_yyz_xyz_1, ta2_xz_yyz_yyy_0, ta2_xz_yyz_yyy_1, ta2_xz_yyz_yyz_0, ta2_xz_yyz_yyz_1, ta2_xz_yyz_yz_0, ta2_xz_yyz_yz_1, ta2_xz_yyz_yzz_0, ta2_xz_yyz_yzz_1, ta2_xz_yyz_zzz_0, ta2_xz_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyz_xxx_0[i] = ta1_x_xyy_xxx_1[i] + ta2_xz_xyy_xxx_0[i] * pa_z[i] - ta2_xz_xyy_xxx_1[i] * pc_z[i];

        ta2_xz_xyyz_xxy_0[i] = ta1_x_xyy_xxy_1[i] + ta2_xz_xyy_xxy_0[i] * pa_z[i] - ta2_xz_xyy_xxy_1[i] * pc_z[i];

        ta2_xz_xyyz_xxz_0[i] = ta2_xz_xz_xxz_0[i] * fe_0 - ta2_xz_xz_xxz_1[i] * fe_0 + ta2_xz_xyz_xxz_0[i] * pa_y[i] - ta2_xz_xyz_xxz_1[i] * pc_y[i];

        ta2_xz_xyyz_xyy_0[i] = ta1_x_xyy_xyy_1[i] + ta2_xz_xyy_xyy_0[i] * pa_z[i] - ta2_xz_xyy_xyy_1[i] * pc_z[i];

        ta2_xz_xyyz_xyz_0[i] = ta2_xz_yyz_yz_0[i] * fe_0 - ta2_xz_yyz_yz_1[i] * fe_0 + ta1_z_yyz_xyz_1[i] + ta2_xz_yyz_xyz_0[i] * pa_x[i] - ta2_xz_yyz_xyz_1[i] * pc_x[i];

        ta2_xz_xyyz_xzz_0[i] = ta2_xz_xz_xzz_0[i] * fe_0 - ta2_xz_xz_xzz_1[i] * fe_0 + ta2_xz_xyz_xzz_0[i] * pa_y[i] - ta2_xz_xyz_xzz_1[i] * pc_y[i];

        ta2_xz_xyyz_yyy_0[i] = ta1_z_yyz_yyy_1[i] + ta2_xz_yyz_yyy_0[i] * pa_x[i] - ta2_xz_yyz_yyy_1[i] * pc_x[i];

        ta2_xz_xyyz_yyz_0[i] = ta1_z_yyz_yyz_1[i] + ta2_xz_yyz_yyz_0[i] * pa_x[i] - ta2_xz_yyz_yyz_1[i] * pc_x[i];

        ta2_xz_xyyz_yzz_0[i] = ta1_z_yyz_yzz_1[i] + ta2_xz_yyz_yzz_0[i] * pa_x[i] - ta2_xz_yyz_yzz_1[i] * pc_x[i];

        ta2_xz_xyyz_zzz_0[i] = ta1_z_yyz_zzz_1[i] + ta2_xz_yyz_zzz_0[i] * pa_x[i] - ta2_xz_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 380-390 components of targeted buffer : GF

    auto ta2_xz_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 380);

    auto ta2_xz_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 381);

    auto ta2_xz_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 382);

    auto ta2_xz_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 383);

    auto ta2_xz_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 384);

    auto ta2_xz_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 385);

    auto ta2_xz_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 386);

    auto ta2_xz_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 387);

    auto ta2_xz_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 388);

    auto ta2_xz_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 389);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_yzz_yyy_1, ta1_z_yzz_yyz_1, ta1_z_yzz_yzz_1, ta2_xz_xyzz_xxx_0, ta2_xz_xyzz_xxy_0, ta2_xz_xyzz_xxz_0, ta2_xz_xyzz_xyy_0, ta2_xz_xyzz_xyz_0, ta2_xz_xyzz_xzz_0, ta2_xz_xyzz_yyy_0, ta2_xz_xyzz_yyz_0, ta2_xz_xyzz_yzz_0, ta2_xz_xyzz_zzz_0, ta2_xz_xzz_xx_0, ta2_xz_xzz_xx_1, ta2_xz_xzz_xxx_0, ta2_xz_xzz_xxx_1, ta2_xz_xzz_xxy_0, ta2_xz_xzz_xxy_1, ta2_xz_xzz_xxz_0, ta2_xz_xzz_xxz_1, ta2_xz_xzz_xy_0, ta2_xz_xzz_xy_1, ta2_xz_xzz_xyy_0, ta2_xz_xzz_xyy_1, ta2_xz_xzz_xyz_0, ta2_xz_xzz_xyz_1, ta2_xz_xzz_xz_0, ta2_xz_xzz_xz_1, ta2_xz_xzz_xzz_0, ta2_xz_xzz_xzz_1, ta2_xz_xzz_zzz_0, ta2_xz_xzz_zzz_1, ta2_xz_yzz_yyy_0, ta2_xz_yzz_yyy_1, ta2_xz_yzz_yyz_0, ta2_xz_yzz_yyz_1, ta2_xz_yzz_yzz_0, ta2_xz_yzz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyzz_xxx_0[i] = ta2_xz_xzz_xxx_0[i] * pa_y[i] - ta2_xz_xzz_xxx_1[i] * pc_y[i];

        ta2_xz_xyzz_xxy_0[i] = ta2_xz_xzz_xx_0[i] * fe_0 - ta2_xz_xzz_xx_1[i] * fe_0 + ta2_xz_xzz_xxy_0[i] * pa_y[i] - ta2_xz_xzz_xxy_1[i] * pc_y[i];

        ta2_xz_xyzz_xxz_0[i] = ta2_xz_xzz_xxz_0[i] * pa_y[i] - ta2_xz_xzz_xxz_1[i] * pc_y[i];

        ta2_xz_xyzz_xyy_0[i] = 2.0 * ta2_xz_xzz_xy_0[i] * fe_0 - 2.0 * ta2_xz_xzz_xy_1[i] * fe_0 + ta2_xz_xzz_xyy_0[i] * pa_y[i] - ta2_xz_xzz_xyy_1[i] * pc_y[i];

        ta2_xz_xyzz_xyz_0[i] = ta2_xz_xzz_xz_0[i] * fe_0 - ta2_xz_xzz_xz_1[i] * fe_0 + ta2_xz_xzz_xyz_0[i] * pa_y[i] - ta2_xz_xzz_xyz_1[i] * pc_y[i];

        ta2_xz_xyzz_xzz_0[i] = ta2_xz_xzz_xzz_0[i] * pa_y[i] - ta2_xz_xzz_xzz_1[i] * pc_y[i];

        ta2_xz_xyzz_yyy_0[i] = ta1_z_yzz_yyy_1[i] + ta2_xz_yzz_yyy_0[i] * pa_x[i] - ta2_xz_yzz_yyy_1[i] * pc_x[i];

        ta2_xz_xyzz_yyz_0[i] = ta1_z_yzz_yyz_1[i] + ta2_xz_yzz_yyz_0[i] * pa_x[i] - ta2_xz_yzz_yyz_1[i] * pc_x[i];

        ta2_xz_xyzz_yzz_0[i] = ta1_z_yzz_yzz_1[i] + ta2_xz_yzz_yzz_0[i] * pa_x[i] - ta2_xz_yzz_yzz_1[i] * pc_x[i];

        ta2_xz_xyzz_zzz_0[i] = ta2_xz_xzz_zzz_0[i] * pa_y[i] - ta2_xz_xzz_zzz_1[i] * pc_y[i];
    }

    // Set up 390-400 components of targeted buffer : GF

    auto ta2_xz_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 390);

    auto ta2_xz_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 391);

    auto ta2_xz_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 392);

    auto ta2_xz_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 393);

    auto ta2_xz_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 394);

    auto ta2_xz_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 395);

    auto ta2_xz_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 396);

    auto ta2_xz_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 397);

    auto ta2_xz_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 398);

    auto ta2_xz_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 399);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_zzz_xxx_1, ta1_z_zzz_xxy_1, ta1_z_zzz_xxz_1, ta1_z_zzz_xyy_1, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_1, ta1_z_zzz_yyy_1, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_1, ta2_xz_xzzz_xxx_0, ta2_xz_xzzz_xxy_0, ta2_xz_xzzz_xxz_0, ta2_xz_xzzz_xyy_0, ta2_xz_xzzz_xyz_0, ta2_xz_xzzz_xzz_0, ta2_xz_xzzz_yyy_0, ta2_xz_xzzz_yyz_0, ta2_xz_xzzz_yzz_0, ta2_xz_xzzz_zzz_0, ta2_xz_zzz_xx_0, ta2_xz_zzz_xx_1, ta2_xz_zzz_xxx_0, ta2_xz_zzz_xxx_1, ta2_xz_zzz_xxy_0, ta2_xz_zzz_xxy_1, ta2_xz_zzz_xxz_0, ta2_xz_zzz_xxz_1, ta2_xz_zzz_xy_0, ta2_xz_zzz_xy_1, ta2_xz_zzz_xyy_0, ta2_xz_zzz_xyy_1, ta2_xz_zzz_xyz_0, ta2_xz_zzz_xyz_1, ta2_xz_zzz_xz_0, ta2_xz_zzz_xz_1, ta2_xz_zzz_xzz_0, ta2_xz_zzz_xzz_1, ta2_xz_zzz_yy_0, ta2_xz_zzz_yy_1, ta2_xz_zzz_yyy_0, ta2_xz_zzz_yyy_1, ta2_xz_zzz_yyz_0, ta2_xz_zzz_yyz_1, ta2_xz_zzz_yz_0, ta2_xz_zzz_yz_1, ta2_xz_zzz_yzz_0, ta2_xz_zzz_yzz_1, ta2_xz_zzz_zz_0, ta2_xz_zzz_zz_1, ta2_xz_zzz_zzz_0, ta2_xz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzzz_xxx_0[i] = 3.0 * ta2_xz_zzz_xx_0[i] * fe_0 - 3.0 * ta2_xz_zzz_xx_1[i] * fe_0 + ta1_z_zzz_xxx_1[i] + ta2_xz_zzz_xxx_0[i] * pa_x[i] - ta2_xz_zzz_xxx_1[i] * pc_x[i];

        ta2_xz_xzzz_xxy_0[i] = 2.0 * ta2_xz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xy_1[i] * fe_0 + ta1_z_zzz_xxy_1[i] + ta2_xz_zzz_xxy_0[i] * pa_x[i] - ta2_xz_zzz_xxy_1[i] * pc_x[i];

        ta2_xz_xzzz_xxz_0[i] = 2.0 * ta2_xz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xz_1[i] * fe_0 + ta1_z_zzz_xxz_1[i] + ta2_xz_zzz_xxz_0[i] * pa_x[i] - ta2_xz_zzz_xxz_1[i] * pc_x[i];

        ta2_xz_xzzz_xyy_0[i] = ta2_xz_zzz_yy_0[i] * fe_0 - ta2_xz_zzz_yy_1[i] * fe_0 + ta1_z_zzz_xyy_1[i] + ta2_xz_zzz_xyy_0[i] * pa_x[i] - ta2_xz_zzz_xyy_1[i] * pc_x[i];

        ta2_xz_xzzz_xyz_0[i] = ta2_xz_zzz_yz_0[i] * fe_0 - ta2_xz_zzz_yz_1[i] * fe_0 + ta1_z_zzz_xyz_1[i] + ta2_xz_zzz_xyz_0[i] * pa_x[i] - ta2_xz_zzz_xyz_1[i] * pc_x[i];

        ta2_xz_xzzz_xzz_0[i] = ta2_xz_zzz_zz_0[i] * fe_0 - ta2_xz_zzz_zz_1[i] * fe_0 + ta1_z_zzz_xzz_1[i] + ta2_xz_zzz_xzz_0[i] * pa_x[i] - ta2_xz_zzz_xzz_1[i] * pc_x[i];

        ta2_xz_xzzz_yyy_0[i] = ta1_z_zzz_yyy_1[i] + ta2_xz_zzz_yyy_0[i] * pa_x[i] - ta2_xz_zzz_yyy_1[i] * pc_x[i];

        ta2_xz_xzzz_yyz_0[i] = ta1_z_zzz_yyz_1[i] + ta2_xz_zzz_yyz_0[i] * pa_x[i] - ta2_xz_zzz_yyz_1[i] * pc_x[i];

        ta2_xz_xzzz_yzz_0[i] = ta1_z_zzz_yzz_1[i] + ta2_xz_zzz_yzz_0[i] * pa_x[i] - ta2_xz_zzz_yzz_1[i] * pc_x[i];

        ta2_xz_xzzz_zzz_0[i] = ta1_z_zzz_zzz_1[i] + ta2_xz_zzz_zzz_0[i] * pa_x[i] - ta2_xz_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 400-410 components of targeted buffer : GF

    auto ta2_xz_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 400);

    auto ta2_xz_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 401);

    auto ta2_xz_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 402);

    auto ta2_xz_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 403);

    auto ta2_xz_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 404);

    auto ta2_xz_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 405);

    auto ta2_xz_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 406);

    auto ta2_xz_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 407);

    auto ta2_xz_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 408);

    auto ta2_xz_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 409);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_yy_xxx_0, ta2_xz_yy_xxx_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xxz_0, ta2_xz_yy_xxz_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_xyz_0, ta2_xz_yy_xyz_1, ta2_xz_yy_xzz_0, ta2_xz_yy_xzz_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yy_yyz_0, ta2_xz_yy_yyz_1, ta2_xz_yy_yzz_0, ta2_xz_yy_yzz_1, ta2_xz_yy_zzz_0, ta2_xz_yy_zzz_1, ta2_xz_yyy_xx_0, ta2_xz_yyy_xx_1, ta2_xz_yyy_xxx_0, ta2_xz_yyy_xxx_1, ta2_xz_yyy_xxy_0, ta2_xz_yyy_xxy_1, ta2_xz_yyy_xxz_0, ta2_xz_yyy_xxz_1, ta2_xz_yyy_xy_0, ta2_xz_yyy_xy_1, ta2_xz_yyy_xyy_0, ta2_xz_yyy_xyy_1, ta2_xz_yyy_xyz_0, ta2_xz_yyy_xyz_1, ta2_xz_yyy_xz_0, ta2_xz_yyy_xz_1, ta2_xz_yyy_xzz_0, ta2_xz_yyy_xzz_1, ta2_xz_yyy_yy_0, ta2_xz_yyy_yy_1, ta2_xz_yyy_yyy_0, ta2_xz_yyy_yyy_1, ta2_xz_yyy_yyz_0, ta2_xz_yyy_yyz_1, ta2_xz_yyy_yz_0, ta2_xz_yyy_yz_1, ta2_xz_yyy_yzz_0, ta2_xz_yyy_yzz_1, ta2_xz_yyy_zz_0, ta2_xz_yyy_zz_1, ta2_xz_yyy_zzz_0, ta2_xz_yyy_zzz_1, ta2_xz_yyyy_xxx_0, ta2_xz_yyyy_xxy_0, ta2_xz_yyyy_xxz_0, ta2_xz_yyyy_xyy_0, ta2_xz_yyyy_xyz_0, ta2_xz_yyyy_xzz_0, ta2_xz_yyyy_yyy_0, ta2_xz_yyyy_yyz_0, ta2_xz_yyyy_yzz_0, ta2_xz_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyy_xxx_0[i] = 3.0 * ta2_xz_yy_xxx_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxx_1[i] * fe_0 + ta2_xz_yyy_xxx_0[i] * pa_y[i] - ta2_xz_yyy_xxx_1[i] * pc_y[i];

        ta2_xz_yyyy_xxy_0[i] = 3.0 * ta2_xz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxy_1[i] * fe_0 + ta2_xz_yyy_xx_0[i] * fe_0 - ta2_xz_yyy_xx_1[i] * fe_0 + ta2_xz_yyy_xxy_0[i] * pa_y[i] - ta2_xz_yyy_xxy_1[i] * pc_y[i];

        ta2_xz_yyyy_xxz_0[i] = 3.0 * ta2_xz_yy_xxz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxz_1[i] * fe_0 + ta2_xz_yyy_xxz_0[i] * pa_y[i] - ta2_xz_yyy_xxz_1[i] * pc_y[i];

        ta2_xz_yyyy_xyy_0[i] = 3.0 * ta2_xz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyy_1[i] * fe_0 + 2.0 * ta2_xz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_xz_yyy_xy_1[i] * fe_0 + ta2_xz_yyy_xyy_0[i] * pa_y[i] - ta2_xz_yyy_xyy_1[i] * pc_y[i];

        ta2_xz_yyyy_xyz_0[i] = 3.0 * ta2_xz_yy_xyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyz_1[i] * fe_0 + ta2_xz_yyy_xz_0[i] * fe_0 - ta2_xz_yyy_xz_1[i] * fe_0 + ta2_xz_yyy_xyz_0[i] * pa_y[i] - ta2_xz_yyy_xyz_1[i] * pc_y[i];

        ta2_xz_yyyy_xzz_0[i] = 3.0 * ta2_xz_yy_xzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xzz_1[i] * fe_0 + ta2_xz_yyy_xzz_0[i] * pa_y[i] - ta2_xz_yyy_xzz_1[i] * pc_y[i];

        ta2_xz_yyyy_yyy_0[i] = 3.0 * ta2_xz_yy_yyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyy_1[i] * fe_0 + 3.0 * ta2_xz_yyy_yy_0[i] * fe_0 - 3.0 * ta2_xz_yyy_yy_1[i] * fe_0 + ta2_xz_yyy_yyy_0[i] * pa_y[i] - ta2_xz_yyy_yyy_1[i] * pc_y[i];

        ta2_xz_yyyy_yyz_0[i] = 3.0 * ta2_xz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyz_1[i] * fe_0 + 2.0 * ta2_xz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xz_yyy_yz_1[i] * fe_0 + ta2_xz_yyy_yyz_0[i] * pa_y[i] - ta2_xz_yyy_yyz_1[i] * pc_y[i];

        ta2_xz_yyyy_yzz_0[i] = 3.0 * ta2_xz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yzz_1[i] * fe_0 + ta2_xz_yyy_zz_0[i] * fe_0 - ta2_xz_yyy_zz_1[i] * fe_0 + ta2_xz_yyy_yzz_0[i] * pa_y[i] - ta2_xz_yyy_yzz_1[i] * pc_y[i];

        ta2_xz_yyyy_zzz_0[i] = 3.0 * ta2_xz_yy_zzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_zzz_1[i] * fe_0 + ta2_xz_yyy_zzz_0[i] * pa_y[i] - ta2_xz_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 410-420 components of targeted buffer : GF

    auto ta2_xz_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 410);

    auto ta2_xz_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 411);

    auto ta2_xz_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 412);

    auto ta2_xz_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 413);

    auto ta2_xz_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 414);

    auto ta2_xz_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 415);

    auto ta2_xz_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 416);

    auto ta2_xz_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 417);

    auto ta2_xz_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 418);

    auto ta2_xz_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 419);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyy_xxx_1, ta1_x_yyy_xxy_1, ta1_x_yyy_xyy_1, ta1_x_yyy_xyz_1, ta1_x_yyy_yyy_1, ta1_x_yyy_yyz_1, ta1_x_yyy_yzz_1, ta2_xz_yyy_xxx_0, ta2_xz_yyy_xxx_1, ta2_xz_yyy_xxy_0, ta2_xz_yyy_xxy_1, ta2_xz_yyy_xy_0, ta2_xz_yyy_xy_1, ta2_xz_yyy_xyy_0, ta2_xz_yyy_xyy_1, ta2_xz_yyy_xyz_0, ta2_xz_yyy_xyz_1, ta2_xz_yyy_yy_0, ta2_xz_yyy_yy_1, ta2_xz_yyy_yyy_0, ta2_xz_yyy_yyy_1, ta2_xz_yyy_yyz_0, ta2_xz_yyy_yyz_1, ta2_xz_yyy_yz_0, ta2_xz_yyy_yz_1, ta2_xz_yyy_yzz_0, ta2_xz_yyy_yzz_1, ta2_xz_yyyz_xxx_0, ta2_xz_yyyz_xxy_0, ta2_xz_yyyz_xxz_0, ta2_xz_yyyz_xyy_0, ta2_xz_yyyz_xyz_0, ta2_xz_yyyz_xzz_0, ta2_xz_yyyz_yyy_0, ta2_xz_yyyz_yyz_0, ta2_xz_yyyz_yzz_0, ta2_xz_yyyz_zzz_0, ta2_xz_yyz_xxz_0, ta2_xz_yyz_xxz_1, ta2_xz_yyz_xzz_0, ta2_xz_yyz_xzz_1, ta2_xz_yyz_zzz_0, ta2_xz_yyz_zzz_1, ta2_xz_yz_xxz_0, ta2_xz_yz_xxz_1, ta2_xz_yz_xzz_0, ta2_xz_yz_xzz_1, ta2_xz_yz_zzz_0, ta2_xz_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyz_xxx_0[i] = ta1_x_yyy_xxx_1[i] + ta2_xz_yyy_xxx_0[i] * pa_z[i] - ta2_xz_yyy_xxx_1[i] * pc_z[i];

        ta2_xz_yyyz_xxy_0[i] = ta1_x_yyy_xxy_1[i] + ta2_xz_yyy_xxy_0[i] * pa_z[i] - ta2_xz_yyy_xxy_1[i] * pc_z[i];

        ta2_xz_yyyz_xxz_0[i] = 2.0 * ta2_xz_yz_xxz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xxz_1[i] * fe_0 + ta2_xz_yyz_xxz_0[i] * pa_y[i] - ta2_xz_yyz_xxz_1[i] * pc_y[i];

        ta2_xz_yyyz_xyy_0[i] = ta1_x_yyy_xyy_1[i] + ta2_xz_yyy_xyy_0[i] * pa_z[i] - ta2_xz_yyy_xyy_1[i] * pc_z[i];

        ta2_xz_yyyz_xyz_0[i] = ta2_xz_yyy_xy_0[i] * fe_0 - ta2_xz_yyy_xy_1[i] * fe_0 + ta1_x_yyy_xyz_1[i] + ta2_xz_yyy_xyz_0[i] * pa_z[i] - ta2_xz_yyy_xyz_1[i] * pc_z[i];

        ta2_xz_yyyz_xzz_0[i] = 2.0 * ta2_xz_yz_xzz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xzz_1[i] * fe_0 + ta2_xz_yyz_xzz_0[i] * pa_y[i] - ta2_xz_yyz_xzz_1[i] * pc_y[i];

        ta2_xz_yyyz_yyy_0[i] = ta1_x_yyy_yyy_1[i] + ta2_xz_yyy_yyy_0[i] * pa_z[i] - ta2_xz_yyy_yyy_1[i] * pc_z[i];

        ta2_xz_yyyz_yyz_0[i] = ta2_xz_yyy_yy_0[i] * fe_0 - ta2_xz_yyy_yy_1[i] * fe_0 + ta1_x_yyy_yyz_1[i] + ta2_xz_yyy_yyz_0[i] * pa_z[i] - ta2_xz_yyy_yyz_1[i] * pc_z[i];

        ta2_xz_yyyz_yzz_0[i] = 2.0 * ta2_xz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_xz_yyy_yz_1[i] * fe_0 + ta1_x_yyy_yzz_1[i] + ta2_xz_yyy_yzz_0[i] * pa_z[i] - ta2_xz_yyy_yzz_1[i] * pc_z[i];

        ta2_xz_yyyz_zzz_0[i] = 2.0 * ta2_xz_yz_zzz_0[i] * fe_0 - 2.0 * ta2_xz_yz_zzz_1[i] * fe_0 + ta2_xz_yyz_zzz_0[i] * pa_y[i] - ta2_xz_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 420-430 components of targeted buffer : GF

    auto ta2_xz_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 420);

    auto ta2_xz_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 421);

    auto ta2_xz_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 422);

    auto ta2_xz_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 423);

    auto ta2_xz_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 424);

    auto ta2_xz_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 425);

    auto ta2_xz_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 426);

    auto ta2_xz_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 427);

    auto ta2_xz_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 428);

    auto ta2_xz_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 429);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyz_xxy_1, ta1_x_yyz_xyy_1, ta1_x_yyz_yyy_1, ta2_xz_yy_xxy_0, ta2_xz_yy_xxy_1, ta2_xz_yy_xyy_0, ta2_xz_yy_xyy_1, ta2_xz_yy_yyy_0, ta2_xz_yy_yyy_1, ta2_xz_yyz_xxy_0, ta2_xz_yyz_xxy_1, ta2_xz_yyz_xyy_0, ta2_xz_yyz_xyy_1, ta2_xz_yyz_yyy_0, ta2_xz_yyz_yyy_1, ta2_xz_yyzz_xxx_0, ta2_xz_yyzz_xxy_0, ta2_xz_yyzz_xxz_0, ta2_xz_yyzz_xyy_0, ta2_xz_yyzz_xyz_0, ta2_xz_yyzz_xzz_0, ta2_xz_yyzz_yyy_0, ta2_xz_yyzz_yyz_0, ta2_xz_yyzz_yzz_0, ta2_xz_yyzz_zzz_0, ta2_xz_yzz_xxx_0, ta2_xz_yzz_xxx_1, ta2_xz_yzz_xxz_0, ta2_xz_yzz_xxz_1, ta2_xz_yzz_xyz_0, ta2_xz_yzz_xyz_1, ta2_xz_yzz_xz_0, ta2_xz_yzz_xz_1, ta2_xz_yzz_xzz_0, ta2_xz_yzz_xzz_1, ta2_xz_yzz_yyz_0, ta2_xz_yzz_yyz_1, ta2_xz_yzz_yz_0, ta2_xz_yzz_yz_1, ta2_xz_yzz_yzz_0, ta2_xz_yzz_yzz_1, ta2_xz_yzz_zz_0, ta2_xz_yzz_zz_1, ta2_xz_yzz_zzz_0, ta2_xz_yzz_zzz_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyzz_xxx_0[i] = ta2_xz_zz_xxx_0[i] * fe_0 - ta2_xz_zz_xxx_1[i] * fe_0 + ta2_xz_yzz_xxx_0[i] * pa_y[i] - ta2_xz_yzz_xxx_1[i] * pc_y[i];

        ta2_xz_yyzz_xxy_0[i] = ta2_xz_yy_xxy_0[i] * fe_0 - ta2_xz_yy_xxy_1[i] * fe_0 + ta1_x_yyz_xxy_1[i] + ta2_xz_yyz_xxy_0[i] * pa_z[i] - ta2_xz_yyz_xxy_1[i] * pc_z[i];

        ta2_xz_yyzz_xxz_0[i] = ta2_xz_zz_xxz_0[i] * fe_0 - ta2_xz_zz_xxz_1[i] * fe_0 + ta2_xz_yzz_xxz_0[i] * pa_y[i] - ta2_xz_yzz_xxz_1[i] * pc_y[i];

        ta2_xz_yyzz_xyy_0[i] = ta2_xz_yy_xyy_0[i] * fe_0 - ta2_xz_yy_xyy_1[i] * fe_0 + ta1_x_yyz_xyy_1[i] + ta2_xz_yyz_xyy_0[i] * pa_z[i] - ta2_xz_yyz_xyy_1[i] * pc_z[i];

        ta2_xz_yyzz_xyz_0[i] = ta2_xz_zz_xyz_0[i] * fe_0 - ta2_xz_zz_xyz_1[i] * fe_0 + ta2_xz_yzz_xz_0[i] * fe_0 - ta2_xz_yzz_xz_1[i] * fe_0 + ta2_xz_yzz_xyz_0[i] * pa_y[i] - ta2_xz_yzz_xyz_1[i] * pc_y[i];

        ta2_xz_yyzz_xzz_0[i] = ta2_xz_zz_xzz_0[i] * fe_0 - ta2_xz_zz_xzz_1[i] * fe_0 + ta2_xz_yzz_xzz_0[i] * pa_y[i] - ta2_xz_yzz_xzz_1[i] * pc_y[i];

        ta2_xz_yyzz_yyy_0[i] = ta2_xz_yy_yyy_0[i] * fe_0 - ta2_xz_yy_yyy_1[i] * fe_0 + ta1_x_yyz_yyy_1[i] + ta2_xz_yyz_yyy_0[i] * pa_z[i] - ta2_xz_yyz_yyy_1[i] * pc_z[i];

        ta2_xz_yyzz_yyz_0[i] = ta2_xz_zz_yyz_0[i] * fe_0 - ta2_xz_zz_yyz_1[i] * fe_0 + 2.0 * ta2_xz_yzz_yz_0[i] * fe_0 - 2.0 * ta2_xz_yzz_yz_1[i] * fe_0 + ta2_xz_yzz_yyz_0[i] * pa_y[i] - ta2_xz_yzz_yyz_1[i] * pc_y[i];

        ta2_xz_yyzz_yzz_0[i] = ta2_xz_zz_yzz_0[i] * fe_0 - ta2_xz_zz_yzz_1[i] * fe_0 + ta2_xz_yzz_zz_0[i] * fe_0 - ta2_xz_yzz_zz_1[i] * fe_0 + ta2_xz_yzz_yzz_0[i] * pa_y[i] - ta2_xz_yzz_yzz_1[i] * pc_y[i];

        ta2_xz_yyzz_zzz_0[i] = ta2_xz_zz_zzz_0[i] * fe_0 - ta2_xz_zz_zzz_1[i] * fe_0 + ta2_xz_yzz_zzz_0[i] * pa_y[i] - ta2_xz_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 430-440 components of targeted buffer : GF

    auto ta2_xz_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 430);

    auto ta2_xz_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 431);

    auto ta2_xz_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 432);

    auto ta2_xz_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 433);

    auto ta2_xz_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 434);

    auto ta2_xz_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 435);

    auto ta2_xz_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 436);

    auto ta2_xz_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 437);

    auto ta2_xz_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 438);

    auto ta2_xz_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 439);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_yzzz_xxx_0, ta2_xz_yzzz_xxy_0, ta2_xz_yzzz_xxz_0, ta2_xz_yzzz_xyy_0, ta2_xz_yzzz_xyz_0, ta2_xz_yzzz_xzz_0, ta2_xz_yzzz_yyy_0, ta2_xz_yzzz_yyz_0, ta2_xz_yzzz_yzz_0, ta2_xz_yzzz_zzz_0, ta2_xz_zzz_xx_0, ta2_xz_zzz_xx_1, ta2_xz_zzz_xxx_0, ta2_xz_zzz_xxx_1, ta2_xz_zzz_xxy_0, ta2_xz_zzz_xxy_1, ta2_xz_zzz_xxz_0, ta2_xz_zzz_xxz_1, ta2_xz_zzz_xy_0, ta2_xz_zzz_xy_1, ta2_xz_zzz_xyy_0, ta2_xz_zzz_xyy_1, ta2_xz_zzz_xyz_0, ta2_xz_zzz_xyz_1, ta2_xz_zzz_xz_0, ta2_xz_zzz_xz_1, ta2_xz_zzz_xzz_0, ta2_xz_zzz_xzz_1, ta2_xz_zzz_yy_0, ta2_xz_zzz_yy_1, ta2_xz_zzz_yyy_0, ta2_xz_zzz_yyy_1, ta2_xz_zzz_yyz_0, ta2_xz_zzz_yyz_1, ta2_xz_zzz_yz_0, ta2_xz_zzz_yz_1, ta2_xz_zzz_yzz_0, ta2_xz_zzz_yzz_1, ta2_xz_zzz_zz_0, ta2_xz_zzz_zz_1, ta2_xz_zzz_zzz_0, ta2_xz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzzz_xxx_0[i] = ta2_xz_zzz_xxx_0[i] * pa_y[i] - ta2_xz_zzz_xxx_1[i] * pc_y[i];

        ta2_xz_yzzz_xxy_0[i] = ta2_xz_zzz_xx_0[i] * fe_0 - ta2_xz_zzz_xx_1[i] * fe_0 + ta2_xz_zzz_xxy_0[i] * pa_y[i] - ta2_xz_zzz_xxy_1[i] * pc_y[i];

        ta2_xz_yzzz_xxz_0[i] = ta2_xz_zzz_xxz_0[i] * pa_y[i] - ta2_xz_zzz_xxz_1[i] * pc_y[i];

        ta2_xz_yzzz_xyy_0[i] = 2.0 * ta2_xz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xy_1[i] * fe_0 + ta2_xz_zzz_xyy_0[i] * pa_y[i] - ta2_xz_zzz_xyy_1[i] * pc_y[i];

        ta2_xz_yzzz_xyz_0[i] = ta2_xz_zzz_xz_0[i] * fe_0 - ta2_xz_zzz_xz_1[i] * fe_0 + ta2_xz_zzz_xyz_0[i] * pa_y[i] - ta2_xz_zzz_xyz_1[i] * pc_y[i];

        ta2_xz_yzzz_xzz_0[i] = ta2_xz_zzz_xzz_0[i] * pa_y[i] - ta2_xz_zzz_xzz_1[i] * pc_y[i];

        ta2_xz_yzzz_yyy_0[i] = 3.0 * ta2_xz_zzz_yy_0[i] * fe_0 - 3.0 * ta2_xz_zzz_yy_1[i] * fe_0 + ta2_xz_zzz_yyy_0[i] * pa_y[i] - ta2_xz_zzz_yyy_1[i] * pc_y[i];

        ta2_xz_yzzz_yyz_0[i] = 2.0 * ta2_xz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_yz_1[i] * fe_0 + ta2_xz_zzz_yyz_0[i] * pa_y[i] - ta2_xz_zzz_yyz_1[i] * pc_y[i];

        ta2_xz_yzzz_yzz_0[i] = ta2_xz_zzz_zz_0[i] * fe_0 - ta2_xz_zzz_zz_1[i] * fe_0 + ta2_xz_zzz_yzz_0[i] * pa_y[i] - ta2_xz_zzz_yzz_1[i] * pc_y[i];

        ta2_xz_yzzz_zzz_0[i] = ta2_xz_zzz_zzz_0[i] * pa_y[i] - ta2_xz_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 440-450 components of targeted buffer : GF

    auto ta2_xz_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 440);

    auto ta2_xz_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 441);

    auto ta2_xz_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 442);

    auto ta2_xz_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 443);

    auto ta2_xz_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 444);

    auto ta2_xz_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 445);

    auto ta2_xz_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 446);

    auto ta2_xz_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 447);

    auto ta2_xz_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 448);

    auto ta2_xz_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 449);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zzz_xxx_1, ta1_x_zzz_xxy_1, ta1_x_zzz_xxz_1, ta1_x_zzz_xyy_1, ta1_x_zzz_xyz_1, ta1_x_zzz_xzz_1, ta1_x_zzz_yyy_1, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_1, ta2_xz_zz_xxx_0, ta2_xz_zz_xxx_1, ta2_xz_zz_xxy_0, ta2_xz_zz_xxy_1, ta2_xz_zz_xxz_0, ta2_xz_zz_xxz_1, ta2_xz_zz_xyy_0, ta2_xz_zz_xyy_1, ta2_xz_zz_xyz_0, ta2_xz_zz_xyz_1, ta2_xz_zz_xzz_0, ta2_xz_zz_xzz_1, ta2_xz_zz_yyy_0, ta2_xz_zz_yyy_1, ta2_xz_zz_yyz_0, ta2_xz_zz_yyz_1, ta2_xz_zz_yzz_0, ta2_xz_zz_yzz_1, ta2_xz_zz_zzz_0, ta2_xz_zz_zzz_1, ta2_xz_zzz_xx_0, ta2_xz_zzz_xx_1, ta2_xz_zzz_xxx_0, ta2_xz_zzz_xxx_1, ta2_xz_zzz_xxy_0, ta2_xz_zzz_xxy_1, ta2_xz_zzz_xxz_0, ta2_xz_zzz_xxz_1, ta2_xz_zzz_xy_0, ta2_xz_zzz_xy_1, ta2_xz_zzz_xyy_0, ta2_xz_zzz_xyy_1, ta2_xz_zzz_xyz_0, ta2_xz_zzz_xyz_1, ta2_xz_zzz_xz_0, ta2_xz_zzz_xz_1, ta2_xz_zzz_xzz_0, ta2_xz_zzz_xzz_1, ta2_xz_zzz_yy_0, ta2_xz_zzz_yy_1, ta2_xz_zzz_yyy_0, ta2_xz_zzz_yyy_1, ta2_xz_zzz_yyz_0, ta2_xz_zzz_yyz_1, ta2_xz_zzz_yz_0, ta2_xz_zzz_yz_1, ta2_xz_zzz_yzz_0, ta2_xz_zzz_yzz_1, ta2_xz_zzz_zz_0, ta2_xz_zzz_zz_1, ta2_xz_zzz_zzz_0, ta2_xz_zzz_zzz_1, ta2_xz_zzzz_xxx_0, ta2_xz_zzzz_xxy_0, ta2_xz_zzzz_xxz_0, ta2_xz_zzzz_xyy_0, ta2_xz_zzzz_xyz_0, ta2_xz_zzzz_xzz_0, ta2_xz_zzzz_yyy_0, ta2_xz_zzzz_yyz_0, ta2_xz_zzzz_yzz_0, ta2_xz_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzzz_xxx_0[i] = 3.0 * ta2_xz_zz_xxx_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxx_1[i] * fe_0 + ta1_x_zzz_xxx_1[i] + ta2_xz_zzz_xxx_0[i] * pa_z[i] - ta2_xz_zzz_xxx_1[i] * pc_z[i];

        ta2_xz_zzzz_xxy_0[i] = 3.0 * ta2_xz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxy_1[i] * fe_0 + ta1_x_zzz_xxy_1[i] + ta2_xz_zzz_xxy_0[i] * pa_z[i] - ta2_xz_zzz_xxy_1[i] * pc_z[i];

        ta2_xz_zzzz_xxz_0[i] = 3.0 * ta2_xz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxz_1[i] * fe_0 + ta2_xz_zzz_xx_0[i] * fe_0 - ta2_xz_zzz_xx_1[i] * fe_0 + ta1_x_zzz_xxz_1[i] + ta2_xz_zzz_xxz_0[i] * pa_z[i] - ta2_xz_zzz_xxz_1[i] * pc_z[i];

        ta2_xz_zzzz_xyy_0[i] = 3.0 * ta2_xz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyy_1[i] * fe_0 + ta1_x_zzz_xyy_1[i] + ta2_xz_zzz_xyy_0[i] * pa_z[i] - ta2_xz_zzz_xyy_1[i] * pc_z[i];

        ta2_xz_zzzz_xyz_0[i] = 3.0 * ta2_xz_zz_xyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyz_1[i] * fe_0 + ta2_xz_zzz_xy_0[i] * fe_0 - ta2_xz_zzz_xy_1[i] * fe_0 + ta1_x_zzz_xyz_1[i] + ta2_xz_zzz_xyz_0[i] * pa_z[i] - ta2_xz_zzz_xyz_1[i] * pc_z[i];

        ta2_xz_zzzz_xzz_0[i] = 3.0 * ta2_xz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xzz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xz_1[i] * fe_0 + ta1_x_zzz_xzz_1[i] + ta2_xz_zzz_xzz_0[i] * pa_z[i] - ta2_xz_zzz_xzz_1[i] * pc_z[i];

        ta2_xz_zzzz_yyy_0[i] = 3.0 * ta2_xz_zz_yyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyy_1[i] * fe_0 + ta1_x_zzz_yyy_1[i] + ta2_xz_zzz_yyy_0[i] * pa_z[i] - ta2_xz_zzz_yyy_1[i] * pc_z[i];

        ta2_xz_zzzz_yyz_0[i] = 3.0 * ta2_xz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyz_1[i] * fe_0 + ta2_xz_zzz_yy_0[i] * fe_0 - ta2_xz_zzz_yy_1[i] * fe_0 + ta1_x_zzz_yyz_1[i] + ta2_xz_zzz_yyz_0[i] * pa_z[i] - ta2_xz_zzz_yyz_1[i] * pc_z[i];

        ta2_xz_zzzz_yzz_0[i] = 3.0 * ta2_xz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yzz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_yz_1[i] * fe_0 + ta1_x_zzz_yzz_1[i] + ta2_xz_zzz_yzz_0[i] * pa_z[i] - ta2_xz_zzz_yzz_1[i] * pc_z[i];

        ta2_xz_zzzz_zzz_0[i] = 3.0 * ta2_xz_zz_zzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_zzz_1[i] * fe_0 + 3.0 * ta2_xz_zzz_zz_0[i] * fe_0 - 3.0 * ta2_xz_zzz_zz_1[i] * fe_0 + ta1_x_zzz_zzz_1[i] + ta2_xz_zzz_zzz_0[i] * pa_z[i] - ta2_xz_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 450-460 components of targeted buffer : GF

    auto ta2_yy_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 450);

    auto ta2_yy_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 451);

    auto ta2_yy_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 452);

    auto ta2_yy_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 453);

    auto ta2_yy_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 454);

    auto ta2_yy_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 455);

    auto ta2_yy_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 456);

    auto ta2_yy_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 457);

    auto ta2_yy_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 458);

    auto ta2_yy_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 459);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xx_xyz_0, ta2_yy_xx_xyz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xx_yyy_0, ta2_yy_xx_yyy_1, ta2_yy_xx_yyz_0, ta2_yy_xx_yyz_1, ta2_yy_xx_yzz_0, ta2_yy_xx_yzz_1, ta2_yy_xx_zzz_0, ta2_yy_xx_zzz_1, ta2_yy_xxx_xx_0, ta2_yy_xxx_xx_1, ta2_yy_xxx_xxx_0, ta2_yy_xxx_xxx_1, ta2_yy_xxx_xxy_0, ta2_yy_xxx_xxy_1, ta2_yy_xxx_xxz_0, ta2_yy_xxx_xxz_1, ta2_yy_xxx_xy_0, ta2_yy_xxx_xy_1, ta2_yy_xxx_xyy_0, ta2_yy_xxx_xyy_1, ta2_yy_xxx_xyz_0, ta2_yy_xxx_xyz_1, ta2_yy_xxx_xz_0, ta2_yy_xxx_xz_1, ta2_yy_xxx_xzz_0, ta2_yy_xxx_xzz_1, ta2_yy_xxx_yy_0, ta2_yy_xxx_yy_1, ta2_yy_xxx_yyy_0, ta2_yy_xxx_yyy_1, ta2_yy_xxx_yyz_0, ta2_yy_xxx_yyz_1, ta2_yy_xxx_yz_0, ta2_yy_xxx_yz_1, ta2_yy_xxx_yzz_0, ta2_yy_xxx_yzz_1, ta2_yy_xxx_zz_0, ta2_yy_xxx_zz_1, ta2_yy_xxx_zzz_0, ta2_yy_xxx_zzz_1, ta2_yy_xxxx_xxx_0, ta2_yy_xxxx_xxy_0, ta2_yy_xxxx_xxz_0, ta2_yy_xxxx_xyy_0, ta2_yy_xxxx_xyz_0, ta2_yy_xxxx_xzz_0, ta2_yy_xxxx_yyy_0, ta2_yy_xxxx_yyz_0, ta2_yy_xxxx_yzz_0, ta2_yy_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxx_xxx_0[i] = 3.0 * ta2_yy_xx_xxx_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxx_1[i] * fe_0 + 3.0 * ta2_yy_xxx_xx_0[i] * fe_0 - 3.0 * ta2_yy_xxx_xx_1[i] * fe_0 + ta2_yy_xxx_xxx_0[i] * pa_x[i] - ta2_yy_xxx_xxx_1[i] * pc_x[i];

        ta2_yy_xxxx_xxy_0[i] = 3.0 * ta2_yy_xx_xxy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxy_1[i] * fe_0 + 2.0 * ta2_yy_xxx_xy_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xy_1[i] * fe_0 + ta2_yy_xxx_xxy_0[i] * pa_x[i] - ta2_yy_xxx_xxy_1[i] * pc_x[i];

        ta2_yy_xxxx_xxz_0[i] = 3.0 * ta2_yy_xx_xxz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxz_1[i] * fe_0 + 2.0 * ta2_yy_xxx_xz_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xz_1[i] * fe_0 + ta2_yy_xxx_xxz_0[i] * pa_x[i] - ta2_yy_xxx_xxz_1[i] * pc_x[i];

        ta2_yy_xxxx_xyy_0[i] = 3.0 * ta2_yy_xx_xyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyy_1[i] * fe_0 + ta2_yy_xxx_yy_0[i] * fe_0 - ta2_yy_xxx_yy_1[i] * fe_0 + ta2_yy_xxx_xyy_0[i] * pa_x[i] - ta2_yy_xxx_xyy_1[i] * pc_x[i];

        ta2_yy_xxxx_xyz_0[i] = 3.0 * ta2_yy_xx_xyz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyz_1[i] * fe_0 + ta2_yy_xxx_yz_0[i] * fe_0 - ta2_yy_xxx_yz_1[i] * fe_0 + ta2_yy_xxx_xyz_0[i] * pa_x[i] - ta2_yy_xxx_xyz_1[i] * pc_x[i];

        ta2_yy_xxxx_xzz_0[i] = 3.0 * ta2_yy_xx_xzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xzz_1[i] * fe_0 + ta2_yy_xxx_zz_0[i] * fe_0 - ta2_yy_xxx_zz_1[i] * fe_0 + ta2_yy_xxx_xzz_0[i] * pa_x[i] - ta2_yy_xxx_xzz_1[i] * pc_x[i];

        ta2_yy_xxxx_yyy_0[i] = 3.0 * ta2_yy_xx_yyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_yyy_1[i] * fe_0 + ta2_yy_xxx_yyy_0[i] * pa_x[i] - ta2_yy_xxx_yyy_1[i] * pc_x[i];

        ta2_yy_xxxx_yyz_0[i] = 3.0 * ta2_yy_xx_yyz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yyz_1[i] * fe_0 + ta2_yy_xxx_yyz_0[i] * pa_x[i] - ta2_yy_xxx_yyz_1[i] * pc_x[i];

        ta2_yy_xxxx_yzz_0[i] = 3.0 * ta2_yy_xx_yzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yzz_1[i] * fe_0 + ta2_yy_xxx_yzz_0[i] * pa_x[i] - ta2_yy_xxx_yzz_1[i] * pc_x[i];

        ta2_yy_xxxx_zzz_0[i] = 3.0 * ta2_yy_xx_zzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_zzz_1[i] * fe_0 + ta2_yy_xxx_zzz_0[i] * pa_x[i] - ta2_yy_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 460-470 components of targeted buffer : GF

    auto ta2_yy_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 460);

    auto ta2_yy_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 461);

    auto ta2_yy_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 462);

    auto ta2_yy_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 463);

    auto ta2_yy_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 464);

    auto ta2_yy_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 465);

    auto ta2_yy_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 466);

    auto ta2_yy_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 467);

    auto ta2_yy_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 468);

    auto ta2_yy_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 469);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxx_xxx_1, ta1_y_xxx_xxy_1, ta1_y_xxx_xxz_1, ta1_y_xxx_xyy_1, ta1_y_xxx_xyz_1, ta1_y_xxx_xzz_1, ta1_y_xxx_zzz_1, ta2_yy_xxx_xx_0, ta2_yy_xxx_xx_1, ta2_yy_xxx_xxx_0, ta2_yy_xxx_xxx_1, ta2_yy_xxx_xxy_0, ta2_yy_xxx_xxy_1, ta2_yy_xxx_xxz_0, ta2_yy_xxx_xxz_1, ta2_yy_xxx_xy_0, ta2_yy_xxx_xy_1, ta2_yy_xxx_xyy_0, ta2_yy_xxx_xyy_1, ta2_yy_xxx_xyz_0, ta2_yy_xxx_xyz_1, ta2_yy_xxx_xz_0, ta2_yy_xxx_xz_1, ta2_yy_xxx_xzz_0, ta2_yy_xxx_xzz_1, ta2_yy_xxx_zzz_0, ta2_yy_xxx_zzz_1, ta2_yy_xxxy_xxx_0, ta2_yy_xxxy_xxy_0, ta2_yy_xxxy_xxz_0, ta2_yy_xxxy_xyy_0, ta2_yy_xxxy_xyz_0, ta2_yy_xxxy_xzz_0, ta2_yy_xxxy_yyy_0, ta2_yy_xxxy_yyz_0, ta2_yy_xxxy_yzz_0, ta2_yy_xxxy_zzz_0, ta2_yy_xxy_yyy_0, ta2_yy_xxy_yyy_1, ta2_yy_xxy_yyz_0, ta2_yy_xxy_yyz_1, ta2_yy_xxy_yzz_0, ta2_yy_xxy_yzz_1, ta2_yy_xy_yyy_0, ta2_yy_xy_yyy_1, ta2_yy_xy_yyz_0, ta2_yy_xy_yyz_1, ta2_yy_xy_yzz_0, ta2_yy_xy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxy_xxx_0[i] = 2.0 * ta1_y_xxx_xxx_1[i] + ta2_yy_xxx_xxx_0[i] * pa_y[i] - ta2_yy_xxx_xxx_1[i] * pc_y[i];

        ta2_yy_xxxy_xxy_0[i] = ta2_yy_xxx_xx_0[i] * fe_0 - ta2_yy_xxx_xx_1[i] * fe_0 + 2.0 * ta1_y_xxx_xxy_1[i] + ta2_yy_xxx_xxy_0[i] * pa_y[i] - ta2_yy_xxx_xxy_1[i] * pc_y[i];

        ta2_yy_xxxy_xxz_0[i] = 2.0 * ta1_y_xxx_xxz_1[i] + ta2_yy_xxx_xxz_0[i] * pa_y[i] - ta2_yy_xxx_xxz_1[i] * pc_y[i];

        ta2_yy_xxxy_xyy_0[i] = 2.0 * ta2_yy_xxx_xy_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyy_1[i] + ta2_yy_xxx_xyy_0[i] * pa_y[i] - ta2_yy_xxx_xyy_1[i] * pc_y[i];

        ta2_yy_xxxy_xyz_0[i] = ta2_yy_xxx_xz_0[i] * fe_0 - ta2_yy_xxx_xz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyz_1[i] + ta2_yy_xxx_xyz_0[i] * pa_y[i] - ta2_yy_xxx_xyz_1[i] * pc_y[i];

        ta2_yy_xxxy_xzz_0[i] = 2.0 * ta1_y_xxx_xzz_1[i] + ta2_yy_xxx_xzz_0[i] * pa_y[i] - ta2_yy_xxx_xzz_1[i] * pc_y[i];

        ta2_yy_xxxy_yyy_0[i] = 2.0 * ta2_yy_xy_yyy_0[i] * fe_0 - 2.0 * ta2_yy_xy_yyy_1[i] * fe_0 + ta2_yy_xxy_yyy_0[i] * pa_x[i] - ta2_yy_xxy_yyy_1[i] * pc_x[i];

        ta2_yy_xxxy_yyz_0[i] = 2.0 * ta2_yy_xy_yyz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yyz_1[i] * fe_0 + ta2_yy_xxy_yyz_0[i] * pa_x[i] - ta2_yy_xxy_yyz_1[i] * pc_x[i];

        ta2_yy_xxxy_yzz_0[i] = 2.0 * ta2_yy_xy_yzz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yzz_1[i] * fe_0 + ta2_yy_xxy_yzz_0[i] * pa_x[i] - ta2_yy_xxy_yzz_1[i] * pc_x[i];

        ta2_yy_xxxy_zzz_0[i] = 2.0 * ta1_y_xxx_zzz_1[i] + ta2_yy_xxx_zzz_0[i] * pa_y[i] - ta2_yy_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 470-480 components of targeted buffer : GF

    auto ta2_yy_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 470);

    auto ta2_yy_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 471);

    auto ta2_yy_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 472);

    auto ta2_yy_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 473);

    auto ta2_yy_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 474);

    auto ta2_yy_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 475);

    auto ta2_yy_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 476);

    auto ta2_yy_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 477);

    auto ta2_yy_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 478);

    auto ta2_yy_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 479);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xxx_xx_0, ta2_yy_xxx_xx_1, ta2_yy_xxx_xxx_0, ta2_yy_xxx_xxx_1, ta2_yy_xxx_xxy_0, ta2_yy_xxx_xxy_1, ta2_yy_xxx_xxz_0, ta2_yy_xxx_xxz_1, ta2_yy_xxx_xy_0, ta2_yy_xxx_xy_1, ta2_yy_xxx_xyy_0, ta2_yy_xxx_xyy_1, ta2_yy_xxx_xyz_0, ta2_yy_xxx_xyz_1, ta2_yy_xxx_xz_0, ta2_yy_xxx_xz_1, ta2_yy_xxx_xzz_0, ta2_yy_xxx_xzz_1, ta2_yy_xxx_yyy_0, ta2_yy_xxx_yyy_1, ta2_yy_xxxz_xxx_0, ta2_yy_xxxz_xxy_0, ta2_yy_xxxz_xxz_0, ta2_yy_xxxz_xyy_0, ta2_yy_xxxz_xyz_0, ta2_yy_xxxz_xzz_0, ta2_yy_xxxz_yyy_0, ta2_yy_xxxz_yyz_0, ta2_yy_xxxz_yzz_0, ta2_yy_xxxz_zzz_0, ta2_yy_xxz_yyz_0, ta2_yy_xxz_yyz_1, ta2_yy_xxz_yzz_0, ta2_yy_xxz_yzz_1, ta2_yy_xxz_zzz_0, ta2_yy_xxz_zzz_1, ta2_yy_xz_yyz_0, ta2_yy_xz_yyz_1, ta2_yy_xz_yzz_0, ta2_yy_xz_yzz_1, ta2_yy_xz_zzz_0, ta2_yy_xz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxz_xxx_0[i] = ta2_yy_xxx_xxx_0[i] * pa_z[i] - ta2_yy_xxx_xxx_1[i] * pc_z[i];

        ta2_yy_xxxz_xxy_0[i] = ta2_yy_xxx_xxy_0[i] * pa_z[i] - ta2_yy_xxx_xxy_1[i] * pc_z[i];

        ta2_yy_xxxz_xxz_0[i] = ta2_yy_xxx_xx_0[i] * fe_0 - ta2_yy_xxx_xx_1[i] * fe_0 + ta2_yy_xxx_xxz_0[i] * pa_z[i] - ta2_yy_xxx_xxz_1[i] * pc_z[i];

        ta2_yy_xxxz_xyy_0[i] = ta2_yy_xxx_xyy_0[i] * pa_z[i] - ta2_yy_xxx_xyy_1[i] * pc_z[i];

        ta2_yy_xxxz_xyz_0[i] = ta2_yy_xxx_xy_0[i] * fe_0 - ta2_yy_xxx_xy_1[i] * fe_0 + ta2_yy_xxx_xyz_0[i] * pa_z[i] - ta2_yy_xxx_xyz_1[i] * pc_z[i];

        ta2_yy_xxxz_xzz_0[i] = 2.0 * ta2_yy_xxx_xz_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xz_1[i] * fe_0 + ta2_yy_xxx_xzz_0[i] * pa_z[i] - ta2_yy_xxx_xzz_1[i] * pc_z[i];

        ta2_yy_xxxz_yyy_0[i] = ta2_yy_xxx_yyy_0[i] * pa_z[i] - ta2_yy_xxx_yyy_1[i] * pc_z[i];

        ta2_yy_xxxz_yyz_0[i] = 2.0 * ta2_yy_xz_yyz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yyz_1[i] * fe_0 + ta2_yy_xxz_yyz_0[i] * pa_x[i] - ta2_yy_xxz_yyz_1[i] * pc_x[i];

        ta2_yy_xxxz_yzz_0[i] = 2.0 * ta2_yy_xz_yzz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yzz_1[i] * fe_0 + ta2_yy_xxz_yzz_0[i] * pa_x[i] - ta2_yy_xxz_yzz_1[i] * pc_x[i];

        ta2_yy_xxxz_zzz_0[i] = 2.0 * ta2_yy_xz_zzz_0[i] * fe_0 - 2.0 * ta2_yy_xz_zzz_1[i] * fe_0 + ta2_yy_xxz_zzz_0[i] * pa_x[i] - ta2_yy_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 480-490 components of targeted buffer : GF

    auto ta2_yy_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 480);

    auto ta2_yy_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 481);

    auto ta2_yy_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 482);

    auto ta2_yy_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 483);

    auto ta2_yy_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 484);

    auto ta2_yy_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 485);

    auto ta2_yy_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 486);

    auto ta2_yy_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 487);

    auto ta2_yy_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 488);

    auto ta2_yy_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 489);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxy_xxx_1, ta1_y_xxy_xxz_1, ta1_y_xxy_xzz_1, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxz_0, ta2_yy_xx_xxz_1, ta2_yy_xx_xzz_0, ta2_yy_xx_xzz_1, ta2_yy_xxy_xxx_0, ta2_yy_xxy_xxx_1, ta2_yy_xxy_xxz_0, ta2_yy_xxy_xxz_1, ta2_yy_xxy_xzz_0, ta2_yy_xxy_xzz_1, ta2_yy_xxyy_xxx_0, ta2_yy_xxyy_xxy_0, ta2_yy_xxyy_xxz_0, ta2_yy_xxyy_xyy_0, ta2_yy_xxyy_xyz_0, ta2_yy_xxyy_xzz_0, ta2_yy_xxyy_yyy_0, ta2_yy_xxyy_yyz_0, ta2_yy_xxyy_yzz_0, ta2_yy_xxyy_zzz_0, ta2_yy_xyy_xxy_0, ta2_yy_xyy_xxy_1, ta2_yy_xyy_xy_0, ta2_yy_xyy_xy_1, ta2_yy_xyy_xyy_0, ta2_yy_xyy_xyy_1, ta2_yy_xyy_xyz_0, ta2_yy_xyy_xyz_1, ta2_yy_xyy_yy_0, ta2_yy_xyy_yy_1, ta2_yy_xyy_yyy_0, ta2_yy_xyy_yyy_1, ta2_yy_xyy_yyz_0, ta2_yy_xyy_yyz_1, ta2_yy_xyy_yz_0, ta2_yy_xyy_yz_1, ta2_yy_xyy_yzz_0, ta2_yy_xyy_yzz_1, ta2_yy_xyy_zzz_0, ta2_yy_xyy_zzz_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyy_xxx_0[i] = ta2_yy_xx_xxx_0[i] * fe_0 - ta2_yy_xx_xxx_1[i] * fe_0 + 2.0 * ta1_y_xxy_xxx_1[i] + ta2_yy_xxy_xxx_0[i] * pa_y[i] - ta2_yy_xxy_xxx_1[i] * pc_y[i];

        ta2_yy_xxyy_xxy_0[i] = ta2_yy_yy_xxy_0[i] * fe_0 - ta2_yy_yy_xxy_1[i] * fe_0 + 2.0 * ta2_yy_xyy_xy_0[i] * fe_0 - 2.0 * ta2_yy_xyy_xy_1[i] * fe_0 + ta2_yy_xyy_xxy_0[i] * pa_x[i] - ta2_yy_xyy_xxy_1[i] * pc_x[i];

        ta2_yy_xxyy_xxz_0[i] = ta2_yy_xx_xxz_0[i] * fe_0 - ta2_yy_xx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xxz_1[i] + ta2_yy_xxy_xxz_0[i] * pa_y[i] - ta2_yy_xxy_xxz_1[i] * pc_y[i];

        ta2_yy_xxyy_xyy_0[i] = ta2_yy_yy_xyy_0[i] * fe_0 - ta2_yy_yy_xyy_1[i] * fe_0 + ta2_yy_xyy_yy_0[i] * fe_0 - ta2_yy_xyy_yy_1[i] * fe_0 + ta2_yy_xyy_xyy_0[i] * pa_x[i] - ta2_yy_xyy_xyy_1[i] * pc_x[i];

        ta2_yy_xxyy_xyz_0[i] = ta2_yy_yy_xyz_0[i] * fe_0 - ta2_yy_yy_xyz_1[i] * fe_0 + ta2_yy_xyy_yz_0[i] * fe_0 - ta2_yy_xyy_yz_1[i] * fe_0 + ta2_yy_xyy_xyz_0[i] * pa_x[i] - ta2_yy_xyy_xyz_1[i] * pc_x[i];

        ta2_yy_xxyy_xzz_0[i] = ta2_yy_xx_xzz_0[i] * fe_0 - ta2_yy_xx_xzz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xzz_1[i] + ta2_yy_xxy_xzz_0[i] * pa_y[i] - ta2_yy_xxy_xzz_1[i] * pc_y[i];

        ta2_yy_xxyy_yyy_0[i] = ta2_yy_yy_yyy_0[i] * fe_0 - ta2_yy_yy_yyy_1[i] * fe_0 + ta2_yy_xyy_yyy_0[i] * pa_x[i] - ta2_yy_xyy_yyy_1[i] * pc_x[i];

        ta2_yy_xxyy_yyz_0[i] = ta2_yy_yy_yyz_0[i] * fe_0 - ta2_yy_yy_yyz_1[i] * fe_0 + ta2_yy_xyy_yyz_0[i] * pa_x[i] - ta2_yy_xyy_yyz_1[i] * pc_x[i];

        ta2_yy_xxyy_yzz_0[i] = ta2_yy_yy_yzz_0[i] * fe_0 - ta2_yy_yy_yzz_1[i] * fe_0 + ta2_yy_xyy_yzz_0[i] * pa_x[i] - ta2_yy_xyy_yzz_1[i] * pc_x[i];

        ta2_yy_xxyy_zzz_0[i] = ta2_yy_yy_zzz_0[i] * fe_0 - ta2_yy_yy_zzz_1[i] * fe_0 + ta2_yy_xyy_zzz_0[i] * pa_x[i] - ta2_yy_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 490-500 components of targeted buffer : GF

    auto ta2_yy_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 490);

    auto ta2_yy_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 491);

    auto ta2_yy_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 492);

    auto ta2_yy_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 493);

    auto ta2_yy_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 494);

    auto ta2_yy_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 495);

    auto ta2_yy_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 496);

    auto ta2_yy_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 497);

    auto ta2_yy_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 498);

    auto ta2_yy_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 499);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxz_xxz_1, ta1_y_xxz_xzz_1, ta1_y_xxz_zzz_1, ta2_yy_xxy_xxx_0, ta2_yy_xxy_xxx_1, ta2_yy_xxy_xxy_0, ta2_yy_xxy_xxy_1, ta2_yy_xxy_xy_0, ta2_yy_xxy_xy_1, ta2_yy_xxy_xyy_0, ta2_yy_xxy_xyy_1, ta2_yy_xxy_xyz_0, ta2_yy_xxy_xyz_1, ta2_yy_xxy_yyy_0, ta2_yy_xxy_yyy_1, ta2_yy_xxyz_xxx_0, ta2_yy_xxyz_xxy_0, ta2_yy_xxyz_xxz_0, ta2_yy_xxyz_xyy_0, ta2_yy_xxyz_xyz_0, ta2_yy_xxyz_xzz_0, ta2_yy_xxyz_yyy_0, ta2_yy_xxyz_yyz_0, ta2_yy_xxyz_yzz_0, ta2_yy_xxyz_zzz_0, ta2_yy_xxz_xxz_0, ta2_yy_xxz_xxz_1, ta2_yy_xxz_xzz_0, ta2_yy_xxz_xzz_1, ta2_yy_xxz_zzz_0, ta2_yy_xxz_zzz_1, ta2_yy_xyz_yyz_0, ta2_yy_xyz_yyz_1, ta2_yy_xyz_yzz_0, ta2_yy_xyz_yzz_1, ta2_yy_yz_yyz_0, ta2_yy_yz_yyz_1, ta2_yy_yz_yzz_0, ta2_yy_yz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyz_xxx_0[i] = ta2_yy_xxy_xxx_0[i] * pa_z[i] - ta2_yy_xxy_xxx_1[i] * pc_z[i];

        ta2_yy_xxyz_xxy_0[i] = ta2_yy_xxy_xxy_0[i] * pa_z[i] - ta2_yy_xxy_xxy_1[i] * pc_z[i];

        ta2_yy_xxyz_xxz_0[i] = 2.0 * ta1_y_xxz_xxz_1[i] + ta2_yy_xxz_xxz_0[i] * pa_y[i] - ta2_yy_xxz_xxz_1[i] * pc_y[i];

        ta2_yy_xxyz_xyy_0[i] = ta2_yy_xxy_xyy_0[i] * pa_z[i] - ta2_yy_xxy_xyy_1[i] * pc_z[i];

        ta2_yy_xxyz_xyz_0[i] = ta2_yy_xxy_xy_0[i] * fe_0 - ta2_yy_xxy_xy_1[i] * fe_0 + ta2_yy_xxy_xyz_0[i] * pa_z[i] - ta2_yy_xxy_xyz_1[i] * pc_z[i];

        ta2_yy_xxyz_xzz_0[i] = 2.0 * ta1_y_xxz_xzz_1[i] + ta2_yy_xxz_xzz_0[i] * pa_y[i] - ta2_yy_xxz_xzz_1[i] * pc_y[i];

        ta2_yy_xxyz_yyy_0[i] = ta2_yy_xxy_yyy_0[i] * pa_z[i] - ta2_yy_xxy_yyy_1[i] * pc_z[i];

        ta2_yy_xxyz_yyz_0[i] = ta2_yy_yz_yyz_0[i] * fe_0 - ta2_yy_yz_yyz_1[i] * fe_0 + ta2_yy_xyz_yyz_0[i] * pa_x[i] - ta2_yy_xyz_yyz_1[i] * pc_x[i];

        ta2_yy_xxyz_yzz_0[i] = ta2_yy_yz_yzz_0[i] * fe_0 - ta2_yy_yz_yzz_1[i] * fe_0 + ta2_yy_xyz_yzz_0[i] * pa_x[i] - ta2_yy_xyz_yzz_1[i] * pc_x[i];

        ta2_yy_xxyz_zzz_0[i] = 2.0 * ta1_y_xxz_zzz_1[i] + ta2_yy_xxz_zzz_0[i] * pa_y[i] - ta2_yy_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 500-510 components of targeted buffer : GF

    auto ta2_yy_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 500);

    auto ta2_yy_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 501);

    auto ta2_yy_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 502);

    auto ta2_yy_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 503);

    auto ta2_yy_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 504);

    auto ta2_yy_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 505);

    auto ta2_yy_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 506);

    auto ta2_yy_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 507);

    auto ta2_yy_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 508);

    auto ta2_yy_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 509);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xx_xxx_0, ta2_yy_xx_xxx_1, ta2_yy_xx_xxy_0, ta2_yy_xx_xxy_1, ta2_yy_xx_xyy_0, ta2_yy_xx_xyy_1, ta2_yy_xxz_xxx_0, ta2_yy_xxz_xxx_1, ta2_yy_xxz_xxy_0, ta2_yy_xxz_xxy_1, ta2_yy_xxz_xyy_0, ta2_yy_xxz_xyy_1, ta2_yy_xxzz_xxx_0, ta2_yy_xxzz_xxy_0, ta2_yy_xxzz_xxz_0, ta2_yy_xxzz_xyy_0, ta2_yy_xxzz_xyz_0, ta2_yy_xxzz_xzz_0, ta2_yy_xxzz_yyy_0, ta2_yy_xxzz_yyz_0, ta2_yy_xxzz_yzz_0, ta2_yy_xxzz_zzz_0, ta2_yy_xzz_xxz_0, ta2_yy_xzz_xxz_1, ta2_yy_xzz_xyz_0, ta2_yy_xzz_xyz_1, ta2_yy_xzz_xz_0, ta2_yy_xzz_xz_1, ta2_yy_xzz_xzz_0, ta2_yy_xzz_xzz_1, ta2_yy_xzz_yyy_0, ta2_yy_xzz_yyy_1, ta2_yy_xzz_yyz_0, ta2_yy_xzz_yyz_1, ta2_yy_xzz_yz_0, ta2_yy_xzz_yz_1, ta2_yy_xzz_yzz_0, ta2_yy_xzz_yzz_1, ta2_yy_xzz_zz_0, ta2_yy_xzz_zz_1, ta2_yy_xzz_zzz_0, ta2_yy_xzz_zzz_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxzz_xxx_0[i] = ta2_yy_xx_xxx_0[i] * fe_0 - ta2_yy_xx_xxx_1[i] * fe_0 + ta2_yy_xxz_xxx_0[i] * pa_z[i] - ta2_yy_xxz_xxx_1[i] * pc_z[i];

        ta2_yy_xxzz_xxy_0[i] = ta2_yy_xx_xxy_0[i] * fe_0 - ta2_yy_xx_xxy_1[i] * fe_0 + ta2_yy_xxz_xxy_0[i] * pa_z[i] - ta2_yy_xxz_xxy_1[i] * pc_z[i];

        ta2_yy_xxzz_xxz_0[i] = ta2_yy_zz_xxz_0[i] * fe_0 - ta2_yy_zz_xxz_1[i] * fe_0 + 2.0 * ta2_yy_xzz_xz_0[i] * fe_0 - 2.0 * ta2_yy_xzz_xz_1[i] * fe_0 + ta2_yy_xzz_xxz_0[i] * pa_x[i] - ta2_yy_xzz_xxz_1[i] * pc_x[i];

        ta2_yy_xxzz_xyy_0[i] = ta2_yy_xx_xyy_0[i] * fe_0 - ta2_yy_xx_xyy_1[i] * fe_0 + ta2_yy_xxz_xyy_0[i] * pa_z[i] - ta2_yy_xxz_xyy_1[i] * pc_z[i];

        ta2_yy_xxzz_xyz_0[i] = ta2_yy_zz_xyz_0[i] * fe_0 - ta2_yy_zz_xyz_1[i] * fe_0 + ta2_yy_xzz_yz_0[i] * fe_0 - ta2_yy_xzz_yz_1[i] * fe_0 + ta2_yy_xzz_xyz_0[i] * pa_x[i] - ta2_yy_xzz_xyz_1[i] * pc_x[i];

        ta2_yy_xxzz_xzz_0[i] = ta2_yy_zz_xzz_0[i] * fe_0 - ta2_yy_zz_xzz_1[i] * fe_0 + ta2_yy_xzz_zz_0[i] * fe_0 - ta2_yy_xzz_zz_1[i] * fe_0 + ta2_yy_xzz_xzz_0[i] * pa_x[i] - ta2_yy_xzz_xzz_1[i] * pc_x[i];

        ta2_yy_xxzz_yyy_0[i] = ta2_yy_zz_yyy_0[i] * fe_0 - ta2_yy_zz_yyy_1[i] * fe_0 + ta2_yy_xzz_yyy_0[i] * pa_x[i] - ta2_yy_xzz_yyy_1[i] * pc_x[i];

        ta2_yy_xxzz_yyz_0[i] = ta2_yy_zz_yyz_0[i] * fe_0 - ta2_yy_zz_yyz_1[i] * fe_0 + ta2_yy_xzz_yyz_0[i] * pa_x[i] - ta2_yy_xzz_yyz_1[i] * pc_x[i];

        ta2_yy_xxzz_yzz_0[i] = ta2_yy_zz_yzz_0[i] * fe_0 - ta2_yy_zz_yzz_1[i] * fe_0 + ta2_yy_xzz_yzz_0[i] * pa_x[i] - ta2_yy_xzz_yzz_1[i] * pc_x[i];

        ta2_yy_xxzz_zzz_0[i] = ta2_yy_zz_zzz_0[i] * fe_0 - ta2_yy_zz_zzz_1[i] * fe_0 + ta2_yy_xzz_zzz_0[i] * pa_x[i] - ta2_yy_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 510-520 components of targeted buffer : GF

    auto ta2_yy_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 510);

    auto ta2_yy_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 511);

    auto ta2_yy_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 512);

    auto ta2_yy_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 513);

    auto ta2_yy_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 514);

    auto ta2_yy_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 515);

    auto ta2_yy_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 516);

    auto ta2_yy_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 517);

    auto ta2_yy_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 518);

    auto ta2_yy_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 519);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xyyy_xxx_0, ta2_yy_xyyy_xxy_0, ta2_yy_xyyy_xxz_0, ta2_yy_xyyy_xyy_0, ta2_yy_xyyy_xyz_0, ta2_yy_xyyy_xzz_0, ta2_yy_xyyy_yyy_0, ta2_yy_xyyy_yyz_0, ta2_yy_xyyy_yzz_0, ta2_yy_xyyy_zzz_0, ta2_yy_yyy_xx_0, ta2_yy_yyy_xx_1, ta2_yy_yyy_xxx_0, ta2_yy_yyy_xxx_1, ta2_yy_yyy_xxy_0, ta2_yy_yyy_xxy_1, ta2_yy_yyy_xxz_0, ta2_yy_yyy_xxz_1, ta2_yy_yyy_xy_0, ta2_yy_yyy_xy_1, ta2_yy_yyy_xyy_0, ta2_yy_yyy_xyy_1, ta2_yy_yyy_xyz_0, ta2_yy_yyy_xyz_1, ta2_yy_yyy_xz_0, ta2_yy_yyy_xz_1, ta2_yy_yyy_xzz_0, ta2_yy_yyy_xzz_1, ta2_yy_yyy_yy_0, ta2_yy_yyy_yy_1, ta2_yy_yyy_yyy_0, ta2_yy_yyy_yyy_1, ta2_yy_yyy_yyz_0, ta2_yy_yyy_yyz_1, ta2_yy_yyy_yz_0, ta2_yy_yyy_yz_1, ta2_yy_yyy_yzz_0, ta2_yy_yyy_yzz_1, ta2_yy_yyy_zz_0, ta2_yy_yyy_zz_1, ta2_yy_yyy_zzz_0, ta2_yy_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyy_xxx_0[i] = 3.0 * ta2_yy_yyy_xx_0[i] * fe_0 - 3.0 * ta2_yy_yyy_xx_1[i] * fe_0 + ta2_yy_yyy_xxx_0[i] * pa_x[i] - ta2_yy_yyy_xxx_1[i] * pc_x[i];

        ta2_yy_xyyy_xxy_0[i] = 2.0 * ta2_yy_yyy_xy_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xy_1[i] * fe_0 + ta2_yy_yyy_xxy_0[i] * pa_x[i] - ta2_yy_yyy_xxy_1[i] * pc_x[i];

        ta2_yy_xyyy_xxz_0[i] = 2.0 * ta2_yy_yyy_xz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xz_1[i] * fe_0 + ta2_yy_yyy_xxz_0[i] * pa_x[i] - ta2_yy_yyy_xxz_1[i] * pc_x[i];

        ta2_yy_xyyy_xyy_0[i] = ta2_yy_yyy_yy_0[i] * fe_0 - ta2_yy_yyy_yy_1[i] * fe_0 + ta2_yy_yyy_xyy_0[i] * pa_x[i] - ta2_yy_yyy_xyy_1[i] * pc_x[i];

        ta2_yy_xyyy_xyz_0[i] = ta2_yy_yyy_yz_0[i] * fe_0 - ta2_yy_yyy_yz_1[i] * fe_0 + ta2_yy_yyy_xyz_0[i] * pa_x[i] - ta2_yy_yyy_xyz_1[i] * pc_x[i];

        ta2_yy_xyyy_xzz_0[i] = ta2_yy_yyy_zz_0[i] * fe_0 - ta2_yy_yyy_zz_1[i] * fe_0 + ta2_yy_yyy_xzz_0[i] * pa_x[i] - ta2_yy_yyy_xzz_1[i] * pc_x[i];

        ta2_yy_xyyy_yyy_0[i] = ta2_yy_yyy_yyy_0[i] * pa_x[i] - ta2_yy_yyy_yyy_1[i] * pc_x[i];

        ta2_yy_xyyy_yyz_0[i] = ta2_yy_yyy_yyz_0[i] * pa_x[i] - ta2_yy_yyy_yyz_1[i] * pc_x[i];

        ta2_yy_xyyy_yzz_0[i] = ta2_yy_yyy_yzz_0[i] * pa_x[i] - ta2_yy_yyy_yzz_1[i] * pc_x[i];

        ta2_yy_xyyy_zzz_0[i] = ta2_yy_yyy_zzz_0[i] * pa_x[i] - ta2_yy_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 520-530 components of targeted buffer : GF

    auto ta2_yy_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 520);

    auto ta2_yy_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 521);

    auto ta2_yy_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 522);

    auto ta2_yy_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 523);

    auto ta2_yy_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 524);

    auto ta2_yy_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 525);

    auto ta2_yy_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 526);

    auto ta2_yy_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 527);

    auto ta2_yy_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 528);

    auto ta2_yy_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 529);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xyy_xxx_0, ta2_yy_xyy_xxx_1, ta2_yy_xyy_xxy_0, ta2_yy_xyy_xxy_1, ta2_yy_xyy_xyy_0, ta2_yy_xyy_xyy_1, ta2_yy_xyyz_xxx_0, ta2_yy_xyyz_xxy_0, ta2_yy_xyyz_xxz_0, ta2_yy_xyyz_xyy_0, ta2_yy_xyyz_xyz_0, ta2_yy_xyyz_xzz_0, ta2_yy_xyyz_yyy_0, ta2_yy_xyyz_yyz_0, ta2_yy_xyyz_yzz_0, ta2_yy_xyyz_zzz_0, ta2_yy_yyz_xxz_0, ta2_yy_yyz_xxz_1, ta2_yy_yyz_xyz_0, ta2_yy_yyz_xyz_1, ta2_yy_yyz_xz_0, ta2_yy_yyz_xz_1, ta2_yy_yyz_xzz_0, ta2_yy_yyz_xzz_1, ta2_yy_yyz_yyy_0, ta2_yy_yyz_yyy_1, ta2_yy_yyz_yyz_0, ta2_yy_yyz_yyz_1, ta2_yy_yyz_yz_0, ta2_yy_yyz_yz_1, ta2_yy_yyz_yzz_0, ta2_yy_yyz_yzz_1, ta2_yy_yyz_zz_0, ta2_yy_yyz_zz_1, ta2_yy_yyz_zzz_0, ta2_yy_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyz_xxx_0[i] = ta2_yy_xyy_xxx_0[i] * pa_z[i] - ta2_yy_xyy_xxx_1[i] * pc_z[i];

        ta2_yy_xyyz_xxy_0[i] = ta2_yy_xyy_xxy_0[i] * pa_z[i] - ta2_yy_xyy_xxy_1[i] * pc_z[i];

        ta2_yy_xyyz_xxz_0[i] = 2.0 * ta2_yy_yyz_xz_0[i] * fe_0 - 2.0 * ta2_yy_yyz_xz_1[i] * fe_0 + ta2_yy_yyz_xxz_0[i] * pa_x[i] - ta2_yy_yyz_xxz_1[i] * pc_x[i];

        ta2_yy_xyyz_xyy_0[i] = ta2_yy_xyy_xyy_0[i] * pa_z[i] - ta2_yy_xyy_xyy_1[i] * pc_z[i];

        ta2_yy_xyyz_xyz_0[i] = ta2_yy_yyz_yz_0[i] * fe_0 - ta2_yy_yyz_yz_1[i] * fe_0 + ta2_yy_yyz_xyz_0[i] * pa_x[i] - ta2_yy_yyz_xyz_1[i] * pc_x[i];

        ta2_yy_xyyz_xzz_0[i] = ta2_yy_yyz_zz_0[i] * fe_0 - ta2_yy_yyz_zz_1[i] * fe_0 + ta2_yy_yyz_xzz_0[i] * pa_x[i] - ta2_yy_yyz_xzz_1[i] * pc_x[i];

        ta2_yy_xyyz_yyy_0[i] = ta2_yy_yyz_yyy_0[i] * pa_x[i] - ta2_yy_yyz_yyy_1[i] * pc_x[i];

        ta2_yy_xyyz_yyz_0[i] = ta2_yy_yyz_yyz_0[i] * pa_x[i] - ta2_yy_yyz_yyz_1[i] * pc_x[i];

        ta2_yy_xyyz_yzz_0[i] = ta2_yy_yyz_yzz_0[i] * pa_x[i] - ta2_yy_yyz_yzz_1[i] * pc_x[i];

        ta2_yy_xyyz_zzz_0[i] = ta2_yy_yyz_zzz_0[i] * pa_x[i] - ta2_yy_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 530-540 components of targeted buffer : GF

    auto ta2_yy_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 530);

    auto ta2_yy_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 531);

    auto ta2_yy_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 532);

    auto ta2_yy_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 533);

    auto ta2_yy_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 534);

    auto ta2_yy_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 535);

    auto ta2_yy_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 536);

    auto ta2_yy_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 537);

    auto ta2_yy_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 538);

    auto ta2_yy_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 539);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xzz_xxx_1, ta1_y_xzz_xxz_1, ta1_y_xzz_xzz_1, ta2_yy_xyzz_xxx_0, ta2_yy_xyzz_xxy_0, ta2_yy_xyzz_xxz_0, ta2_yy_xyzz_xyy_0, ta2_yy_xyzz_xyz_0, ta2_yy_xyzz_xzz_0, ta2_yy_xyzz_yyy_0, ta2_yy_xyzz_yyz_0, ta2_yy_xyzz_yzz_0, ta2_yy_xyzz_zzz_0, ta2_yy_xzz_xxx_0, ta2_yy_xzz_xxx_1, ta2_yy_xzz_xxz_0, ta2_yy_xzz_xxz_1, ta2_yy_xzz_xzz_0, ta2_yy_xzz_xzz_1, ta2_yy_yzz_xxy_0, ta2_yy_yzz_xxy_1, ta2_yy_yzz_xy_0, ta2_yy_yzz_xy_1, ta2_yy_yzz_xyy_0, ta2_yy_yzz_xyy_1, ta2_yy_yzz_xyz_0, ta2_yy_yzz_xyz_1, ta2_yy_yzz_yy_0, ta2_yy_yzz_yy_1, ta2_yy_yzz_yyy_0, ta2_yy_yzz_yyy_1, ta2_yy_yzz_yyz_0, ta2_yy_yzz_yyz_1, ta2_yy_yzz_yz_0, ta2_yy_yzz_yz_1, ta2_yy_yzz_yzz_0, ta2_yy_yzz_yzz_1, ta2_yy_yzz_zzz_0, ta2_yy_yzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyzz_xxx_0[i] = 2.0 * ta1_y_xzz_xxx_1[i] + ta2_yy_xzz_xxx_0[i] * pa_y[i] - ta2_yy_xzz_xxx_1[i] * pc_y[i];

        ta2_yy_xyzz_xxy_0[i] = 2.0 * ta2_yy_yzz_xy_0[i] * fe_0 - 2.0 * ta2_yy_yzz_xy_1[i] * fe_0 + ta2_yy_yzz_xxy_0[i] * pa_x[i] - ta2_yy_yzz_xxy_1[i] * pc_x[i];

        ta2_yy_xyzz_xxz_0[i] = 2.0 * ta1_y_xzz_xxz_1[i] + ta2_yy_xzz_xxz_0[i] * pa_y[i] - ta2_yy_xzz_xxz_1[i] * pc_y[i];

        ta2_yy_xyzz_xyy_0[i] = ta2_yy_yzz_yy_0[i] * fe_0 - ta2_yy_yzz_yy_1[i] * fe_0 + ta2_yy_yzz_xyy_0[i] * pa_x[i] - ta2_yy_yzz_xyy_1[i] * pc_x[i];

        ta2_yy_xyzz_xyz_0[i] = ta2_yy_yzz_yz_0[i] * fe_0 - ta2_yy_yzz_yz_1[i] * fe_0 + ta2_yy_yzz_xyz_0[i] * pa_x[i] - ta2_yy_yzz_xyz_1[i] * pc_x[i];

        ta2_yy_xyzz_xzz_0[i] = 2.0 * ta1_y_xzz_xzz_1[i] + ta2_yy_xzz_xzz_0[i] * pa_y[i] - ta2_yy_xzz_xzz_1[i] * pc_y[i];

        ta2_yy_xyzz_yyy_0[i] = ta2_yy_yzz_yyy_0[i] * pa_x[i] - ta2_yy_yzz_yyy_1[i] * pc_x[i];

        ta2_yy_xyzz_yyz_0[i] = ta2_yy_yzz_yyz_0[i] * pa_x[i] - ta2_yy_yzz_yyz_1[i] * pc_x[i];

        ta2_yy_xyzz_yzz_0[i] = ta2_yy_yzz_yzz_0[i] * pa_x[i] - ta2_yy_yzz_yzz_1[i] * pc_x[i];

        ta2_yy_xyzz_zzz_0[i] = ta2_yy_yzz_zzz_0[i] * pa_x[i] - ta2_yy_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 540-550 components of targeted buffer : GF

    auto ta2_yy_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 540);

    auto ta2_yy_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 541);

    auto ta2_yy_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 542);

    auto ta2_yy_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 543);

    auto ta2_yy_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 544);

    auto ta2_yy_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 545);

    auto ta2_yy_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 546);

    auto ta2_yy_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 547);

    auto ta2_yy_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 548);

    auto ta2_yy_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 549);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xzzz_xxx_0, ta2_yy_xzzz_xxy_0, ta2_yy_xzzz_xxz_0, ta2_yy_xzzz_xyy_0, ta2_yy_xzzz_xyz_0, ta2_yy_xzzz_xzz_0, ta2_yy_xzzz_yyy_0, ta2_yy_xzzz_yyz_0, ta2_yy_xzzz_yzz_0, ta2_yy_xzzz_zzz_0, ta2_yy_zzz_xx_0, ta2_yy_zzz_xx_1, ta2_yy_zzz_xxx_0, ta2_yy_zzz_xxx_1, ta2_yy_zzz_xxy_0, ta2_yy_zzz_xxy_1, ta2_yy_zzz_xxz_0, ta2_yy_zzz_xxz_1, ta2_yy_zzz_xy_0, ta2_yy_zzz_xy_1, ta2_yy_zzz_xyy_0, ta2_yy_zzz_xyy_1, ta2_yy_zzz_xyz_0, ta2_yy_zzz_xyz_1, ta2_yy_zzz_xz_0, ta2_yy_zzz_xz_1, ta2_yy_zzz_xzz_0, ta2_yy_zzz_xzz_1, ta2_yy_zzz_yy_0, ta2_yy_zzz_yy_1, ta2_yy_zzz_yyy_0, ta2_yy_zzz_yyy_1, ta2_yy_zzz_yyz_0, ta2_yy_zzz_yyz_1, ta2_yy_zzz_yz_0, ta2_yy_zzz_yz_1, ta2_yy_zzz_yzz_0, ta2_yy_zzz_yzz_1, ta2_yy_zzz_zz_0, ta2_yy_zzz_zz_1, ta2_yy_zzz_zzz_0, ta2_yy_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzzz_xxx_0[i] = 3.0 * ta2_yy_zzz_xx_0[i] * fe_0 - 3.0 * ta2_yy_zzz_xx_1[i] * fe_0 + ta2_yy_zzz_xxx_0[i] * pa_x[i] - ta2_yy_zzz_xxx_1[i] * pc_x[i];

        ta2_yy_xzzz_xxy_0[i] = 2.0 * ta2_yy_zzz_xy_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xy_1[i] * fe_0 + ta2_yy_zzz_xxy_0[i] * pa_x[i] - ta2_yy_zzz_xxy_1[i] * pc_x[i];

        ta2_yy_xzzz_xxz_0[i] = 2.0 * ta2_yy_zzz_xz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xz_1[i] * fe_0 + ta2_yy_zzz_xxz_0[i] * pa_x[i] - ta2_yy_zzz_xxz_1[i] * pc_x[i];

        ta2_yy_xzzz_xyy_0[i] = ta2_yy_zzz_yy_0[i] * fe_0 - ta2_yy_zzz_yy_1[i] * fe_0 + ta2_yy_zzz_xyy_0[i] * pa_x[i] - ta2_yy_zzz_xyy_1[i] * pc_x[i];

        ta2_yy_xzzz_xyz_0[i] = ta2_yy_zzz_yz_0[i] * fe_0 - ta2_yy_zzz_yz_1[i] * fe_0 + ta2_yy_zzz_xyz_0[i] * pa_x[i] - ta2_yy_zzz_xyz_1[i] * pc_x[i];

        ta2_yy_xzzz_xzz_0[i] = ta2_yy_zzz_zz_0[i] * fe_0 - ta2_yy_zzz_zz_1[i] * fe_0 + ta2_yy_zzz_xzz_0[i] * pa_x[i] - ta2_yy_zzz_xzz_1[i] * pc_x[i];

        ta2_yy_xzzz_yyy_0[i] = ta2_yy_zzz_yyy_0[i] * pa_x[i] - ta2_yy_zzz_yyy_1[i] * pc_x[i];

        ta2_yy_xzzz_yyz_0[i] = ta2_yy_zzz_yyz_0[i] * pa_x[i] - ta2_yy_zzz_yyz_1[i] * pc_x[i];

        ta2_yy_xzzz_yzz_0[i] = ta2_yy_zzz_yzz_0[i] * pa_x[i] - ta2_yy_zzz_yzz_1[i] * pc_x[i];

        ta2_yy_xzzz_zzz_0[i] = ta2_yy_zzz_zzz_0[i] * pa_x[i] - ta2_yy_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 550-560 components of targeted buffer : GF

    auto ta2_yy_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 550);

    auto ta2_yy_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 551);

    auto ta2_yy_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 552);

    auto ta2_yy_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 553);

    auto ta2_yy_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 554);

    auto ta2_yy_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 555);

    auto ta2_yy_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 556);

    auto ta2_yy_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 557);

    auto ta2_yy_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 558);

    auto ta2_yy_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 559);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yyy_xxx_1, ta1_y_yyy_xxy_1, ta1_y_yyy_xxz_1, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_1, ta1_y_yyy_xzz_1, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_1, ta1_y_yyy_zzz_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xxz_0, ta2_yy_yy_xxz_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_xzz_0, ta2_yy_yy_xzz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yy_zzz_0, ta2_yy_yy_zzz_1, ta2_yy_yyy_xx_0, ta2_yy_yyy_xx_1, ta2_yy_yyy_xxx_0, ta2_yy_yyy_xxx_1, ta2_yy_yyy_xxy_0, ta2_yy_yyy_xxy_1, ta2_yy_yyy_xxz_0, ta2_yy_yyy_xxz_1, ta2_yy_yyy_xy_0, ta2_yy_yyy_xy_1, ta2_yy_yyy_xyy_0, ta2_yy_yyy_xyy_1, ta2_yy_yyy_xyz_0, ta2_yy_yyy_xyz_1, ta2_yy_yyy_xz_0, ta2_yy_yyy_xz_1, ta2_yy_yyy_xzz_0, ta2_yy_yyy_xzz_1, ta2_yy_yyy_yy_0, ta2_yy_yyy_yy_1, ta2_yy_yyy_yyy_0, ta2_yy_yyy_yyy_1, ta2_yy_yyy_yyz_0, ta2_yy_yyy_yyz_1, ta2_yy_yyy_yz_0, ta2_yy_yyy_yz_1, ta2_yy_yyy_yzz_0, ta2_yy_yyy_yzz_1, ta2_yy_yyy_zz_0, ta2_yy_yyy_zz_1, ta2_yy_yyy_zzz_0, ta2_yy_yyy_zzz_1, ta2_yy_yyyy_xxx_0, ta2_yy_yyyy_xxy_0, ta2_yy_yyyy_xxz_0, ta2_yy_yyyy_xyy_0, ta2_yy_yyyy_xyz_0, ta2_yy_yyyy_xzz_0, ta2_yy_yyyy_yyy_0, ta2_yy_yyyy_yyz_0, ta2_yy_yyyy_yzz_0, ta2_yy_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyy_xxx_0[i] = 3.0 * ta2_yy_yy_xxx_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxx_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxx_1[i] + ta2_yy_yyy_xxx_0[i] * pa_y[i] - ta2_yy_yyy_xxx_1[i] * pc_y[i];

        ta2_yy_yyyy_xxy_0[i] = 3.0 * ta2_yy_yy_xxy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxy_1[i] * fe_0 + ta2_yy_yyy_xx_0[i] * fe_0 - ta2_yy_yyy_xx_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxy_1[i] + ta2_yy_yyy_xxy_0[i] * pa_y[i] - ta2_yy_yyy_xxy_1[i] * pc_y[i];

        ta2_yy_yyyy_xxz_0[i] = 3.0 * ta2_yy_yy_xxz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxz_1[i] + ta2_yy_yyy_xxz_0[i] * pa_y[i] - ta2_yy_yyy_xxz_1[i] * pc_y[i];

        ta2_yy_yyyy_xyy_0[i] = 3.0 * ta2_yy_yy_xyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyy_1[i] * fe_0 + 2.0 * ta2_yy_yyy_xy_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyy_1[i] + ta2_yy_yyy_xyy_0[i] * pa_y[i] - ta2_yy_yyy_xyy_1[i] * pc_y[i];

        ta2_yy_yyyy_xyz_0[i] = 3.0 * ta2_yy_yy_xyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyz_1[i] * fe_0 + ta2_yy_yyy_xz_0[i] * fe_0 - ta2_yy_yyy_xz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyz_1[i] + ta2_yy_yyy_xyz_0[i] * pa_y[i] - ta2_yy_yyy_xyz_1[i] * pc_y[i];

        ta2_yy_yyyy_xzz_0[i] = 3.0 * ta2_yy_yy_xzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xzz_1[i] + ta2_yy_yyy_xzz_0[i] * pa_y[i] - ta2_yy_yyy_xzz_1[i] * pc_y[i];

        ta2_yy_yyyy_yyy_0[i] = 3.0 * ta2_yy_yy_yyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyy_1[i] * fe_0 + 3.0 * ta2_yy_yyy_yy_0[i] * fe_0 - 3.0 * ta2_yy_yyy_yy_1[i] * fe_0 + 2.0 * ta1_y_yyy_yyy_1[i] + ta2_yy_yyy_yyy_0[i] * pa_y[i] - ta2_yy_yyy_yyy_1[i] * pc_y[i];

        ta2_yy_yyyy_yyz_0[i] = 3.0 * ta2_yy_yy_yyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyz_1[i] * fe_0 + 2.0 * ta2_yy_yyy_yz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_yz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yyz_1[i] + ta2_yy_yyy_yyz_0[i] * pa_y[i] - ta2_yy_yyy_yyz_1[i] * pc_y[i];

        ta2_yy_yyyy_yzz_0[i] = 3.0 * ta2_yy_yy_yzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yzz_1[i] * fe_0 + ta2_yy_yyy_zz_0[i] * fe_0 - ta2_yy_yyy_zz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yzz_1[i] + ta2_yy_yyy_yzz_0[i] * pa_y[i] - ta2_yy_yyy_yzz_1[i] * pc_y[i];

        ta2_yy_yyyy_zzz_0[i] = 3.0 * ta2_yy_yy_zzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_zzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_zzz_1[i] + ta2_yy_yyy_zzz_0[i] * pa_y[i] - ta2_yy_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 560-570 components of targeted buffer : GF

    auto ta2_yy_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 560);

    auto ta2_yy_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 561);

    auto ta2_yy_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 562);

    auto ta2_yy_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 563);

    auto ta2_yy_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 564);

    auto ta2_yy_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 565);

    auto ta2_yy_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 566);

    auto ta2_yy_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 567);

    auto ta2_yy_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 568);

    auto ta2_yy_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 569);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_yyy_xx_0, ta2_yy_yyy_xx_1, ta2_yy_yyy_xxx_0, ta2_yy_yyy_xxx_1, ta2_yy_yyy_xxy_0, ta2_yy_yyy_xxy_1, ta2_yy_yyy_xxz_0, ta2_yy_yyy_xxz_1, ta2_yy_yyy_xy_0, ta2_yy_yyy_xy_1, ta2_yy_yyy_xyy_0, ta2_yy_yyy_xyy_1, ta2_yy_yyy_xyz_0, ta2_yy_yyy_xyz_1, ta2_yy_yyy_xz_0, ta2_yy_yyy_xz_1, ta2_yy_yyy_xzz_0, ta2_yy_yyy_xzz_1, ta2_yy_yyy_yy_0, ta2_yy_yyy_yy_1, ta2_yy_yyy_yyy_0, ta2_yy_yyy_yyy_1, ta2_yy_yyy_yyz_0, ta2_yy_yyy_yyz_1, ta2_yy_yyy_yz_0, ta2_yy_yyy_yz_1, ta2_yy_yyy_yzz_0, ta2_yy_yyy_yzz_1, ta2_yy_yyy_zz_0, ta2_yy_yyy_zz_1, ta2_yy_yyy_zzz_0, ta2_yy_yyy_zzz_1, ta2_yy_yyyz_xxx_0, ta2_yy_yyyz_xxy_0, ta2_yy_yyyz_xxz_0, ta2_yy_yyyz_xyy_0, ta2_yy_yyyz_xyz_0, ta2_yy_yyyz_xzz_0, ta2_yy_yyyz_yyy_0, ta2_yy_yyyz_yyz_0, ta2_yy_yyyz_yzz_0, ta2_yy_yyyz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyz_xxx_0[i] = ta2_yy_yyy_xxx_0[i] * pa_z[i] - ta2_yy_yyy_xxx_1[i] * pc_z[i];

        ta2_yy_yyyz_xxy_0[i] = ta2_yy_yyy_xxy_0[i] * pa_z[i] - ta2_yy_yyy_xxy_1[i] * pc_z[i];

        ta2_yy_yyyz_xxz_0[i] = ta2_yy_yyy_xx_0[i] * fe_0 - ta2_yy_yyy_xx_1[i] * fe_0 + ta2_yy_yyy_xxz_0[i] * pa_z[i] - ta2_yy_yyy_xxz_1[i] * pc_z[i];

        ta2_yy_yyyz_xyy_0[i] = ta2_yy_yyy_xyy_0[i] * pa_z[i] - ta2_yy_yyy_xyy_1[i] * pc_z[i];

        ta2_yy_yyyz_xyz_0[i] = ta2_yy_yyy_xy_0[i] * fe_0 - ta2_yy_yyy_xy_1[i] * fe_0 + ta2_yy_yyy_xyz_0[i] * pa_z[i] - ta2_yy_yyy_xyz_1[i] * pc_z[i];

        ta2_yy_yyyz_xzz_0[i] = 2.0 * ta2_yy_yyy_xz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xz_1[i] * fe_0 + ta2_yy_yyy_xzz_0[i] * pa_z[i] - ta2_yy_yyy_xzz_1[i] * pc_z[i];

        ta2_yy_yyyz_yyy_0[i] = ta2_yy_yyy_yyy_0[i] * pa_z[i] - ta2_yy_yyy_yyy_1[i] * pc_z[i];

        ta2_yy_yyyz_yyz_0[i] = ta2_yy_yyy_yy_0[i] * fe_0 - ta2_yy_yyy_yy_1[i] * fe_0 + ta2_yy_yyy_yyz_0[i] * pa_z[i] - ta2_yy_yyy_yyz_1[i] * pc_z[i];

        ta2_yy_yyyz_yzz_0[i] = 2.0 * ta2_yy_yyy_yz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_yz_1[i] * fe_0 + ta2_yy_yyy_yzz_0[i] * pa_z[i] - ta2_yy_yyy_yzz_1[i] * pc_z[i];

        ta2_yy_yyyz_zzz_0[i] = 3.0 * ta2_yy_yyy_zz_0[i] * fe_0 - 3.0 * ta2_yy_yyy_zz_1[i] * fe_0 + ta2_yy_yyy_zzz_0[i] * pa_z[i] - ta2_yy_yyy_zzz_1[i] * pc_z[i];
    }

    // Set up 570-580 components of targeted buffer : GF

    auto ta2_yy_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 570);

    auto ta2_yy_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 571);

    auto ta2_yy_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 572);

    auto ta2_yy_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 573);

    auto ta2_yy_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 574);

    auto ta2_yy_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 575);

    auto ta2_yy_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 576);

    auto ta2_yy_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 577);

    auto ta2_yy_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 578);

    auto ta2_yy_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 579);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yzz_xxz_1, ta1_y_yzz_xzz_1, ta1_y_yzz_zzz_1, ta2_yy_yy_xxx_0, ta2_yy_yy_xxx_1, ta2_yy_yy_xxy_0, ta2_yy_yy_xxy_1, ta2_yy_yy_xyy_0, ta2_yy_yy_xyy_1, ta2_yy_yy_xyz_0, ta2_yy_yy_xyz_1, ta2_yy_yy_yyy_0, ta2_yy_yy_yyy_1, ta2_yy_yy_yyz_0, ta2_yy_yy_yyz_1, ta2_yy_yy_yzz_0, ta2_yy_yy_yzz_1, ta2_yy_yyz_xxx_0, ta2_yy_yyz_xxx_1, ta2_yy_yyz_xxy_0, ta2_yy_yyz_xxy_1, ta2_yy_yyz_xy_0, ta2_yy_yyz_xy_1, ta2_yy_yyz_xyy_0, ta2_yy_yyz_xyy_1, ta2_yy_yyz_xyz_0, ta2_yy_yyz_xyz_1, ta2_yy_yyz_yy_0, ta2_yy_yyz_yy_1, ta2_yy_yyz_yyy_0, ta2_yy_yyz_yyy_1, ta2_yy_yyz_yyz_0, ta2_yy_yyz_yyz_1, ta2_yy_yyz_yz_0, ta2_yy_yyz_yz_1, ta2_yy_yyz_yzz_0, ta2_yy_yyz_yzz_1, ta2_yy_yyzz_xxx_0, ta2_yy_yyzz_xxy_0, ta2_yy_yyzz_xxz_0, ta2_yy_yyzz_xyy_0, ta2_yy_yyzz_xyz_0, ta2_yy_yyzz_xzz_0, ta2_yy_yyzz_yyy_0, ta2_yy_yyzz_yyz_0, ta2_yy_yyzz_yzz_0, ta2_yy_yyzz_zzz_0, ta2_yy_yzz_xxz_0, ta2_yy_yzz_xxz_1, ta2_yy_yzz_xzz_0, ta2_yy_yzz_xzz_1, ta2_yy_yzz_zzz_0, ta2_yy_yzz_zzz_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyzz_xxx_0[i] = ta2_yy_yy_xxx_0[i] * fe_0 - ta2_yy_yy_xxx_1[i] * fe_0 + ta2_yy_yyz_xxx_0[i] * pa_z[i] - ta2_yy_yyz_xxx_1[i] * pc_z[i];

        ta2_yy_yyzz_xxy_0[i] = ta2_yy_yy_xxy_0[i] * fe_0 - ta2_yy_yy_xxy_1[i] * fe_0 + ta2_yy_yyz_xxy_0[i] * pa_z[i] - ta2_yy_yyz_xxy_1[i] * pc_z[i];

        ta2_yy_yyzz_xxz_0[i] = ta2_yy_zz_xxz_0[i] * fe_0 - ta2_yy_zz_xxz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xxz_1[i] + ta2_yy_yzz_xxz_0[i] * pa_y[i] - ta2_yy_yzz_xxz_1[i] * pc_y[i];

        ta2_yy_yyzz_xyy_0[i] = ta2_yy_yy_xyy_0[i] * fe_0 - ta2_yy_yy_xyy_1[i] * fe_0 + ta2_yy_yyz_xyy_0[i] * pa_z[i] - ta2_yy_yyz_xyy_1[i] * pc_z[i];

        ta2_yy_yyzz_xyz_0[i] = ta2_yy_yy_xyz_0[i] * fe_0 - ta2_yy_yy_xyz_1[i] * fe_0 + ta2_yy_yyz_xy_0[i] * fe_0 - ta2_yy_yyz_xy_1[i] * fe_0 + ta2_yy_yyz_xyz_0[i] * pa_z[i] - ta2_yy_yyz_xyz_1[i] * pc_z[i];

        ta2_yy_yyzz_xzz_0[i] = ta2_yy_zz_xzz_0[i] * fe_0 - ta2_yy_zz_xzz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xzz_1[i] + ta2_yy_yzz_xzz_0[i] * pa_y[i] - ta2_yy_yzz_xzz_1[i] * pc_y[i];

        ta2_yy_yyzz_yyy_0[i] = ta2_yy_yy_yyy_0[i] * fe_0 - ta2_yy_yy_yyy_1[i] * fe_0 + ta2_yy_yyz_yyy_0[i] * pa_z[i] - ta2_yy_yyz_yyy_1[i] * pc_z[i];

        ta2_yy_yyzz_yyz_0[i] = ta2_yy_yy_yyz_0[i] * fe_0 - ta2_yy_yy_yyz_1[i] * fe_0 + ta2_yy_yyz_yy_0[i] * fe_0 - ta2_yy_yyz_yy_1[i] * fe_0 + ta2_yy_yyz_yyz_0[i] * pa_z[i] - ta2_yy_yyz_yyz_1[i] * pc_z[i];

        ta2_yy_yyzz_yzz_0[i] = ta2_yy_yy_yzz_0[i] * fe_0 - ta2_yy_yy_yzz_1[i] * fe_0 + 2.0 * ta2_yy_yyz_yz_0[i] * fe_0 - 2.0 * ta2_yy_yyz_yz_1[i] * fe_0 + ta2_yy_yyz_yzz_0[i] * pa_z[i] - ta2_yy_yyz_yzz_1[i] * pc_z[i];

        ta2_yy_yyzz_zzz_0[i] = ta2_yy_zz_zzz_0[i] * fe_0 - ta2_yy_zz_zzz_1[i] * fe_0 + 2.0 * ta1_y_yzz_zzz_1[i] + ta2_yy_yzz_zzz_0[i] * pa_y[i] - ta2_yy_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 580-590 components of targeted buffer : GF

    auto ta2_yy_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 580);

    auto ta2_yy_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 581);

    auto ta2_yy_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 582);

    auto ta2_yy_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 583);

    auto ta2_yy_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 584);

    auto ta2_yy_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 585);

    auto ta2_yy_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 586);

    auto ta2_yy_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 587);

    auto ta2_yy_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 588);

    auto ta2_yy_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 589);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_zzz_xxx_1, ta1_y_zzz_xxz_1, ta1_y_zzz_xyz_1, ta1_y_zzz_xzz_1, ta1_y_zzz_yyz_1, ta1_y_zzz_yzz_1, ta1_y_zzz_zzz_1, ta2_yy_yz_xxy_0, ta2_yy_yz_xxy_1, ta2_yy_yz_xyy_0, ta2_yy_yz_xyy_1, ta2_yy_yz_yyy_0, ta2_yy_yz_yyy_1, ta2_yy_yzz_xxy_0, ta2_yy_yzz_xxy_1, ta2_yy_yzz_xyy_0, ta2_yy_yzz_xyy_1, ta2_yy_yzz_yyy_0, ta2_yy_yzz_yyy_1, ta2_yy_yzzz_xxx_0, ta2_yy_yzzz_xxy_0, ta2_yy_yzzz_xxz_0, ta2_yy_yzzz_xyy_0, ta2_yy_yzzz_xyz_0, ta2_yy_yzzz_xzz_0, ta2_yy_yzzz_yyy_0, ta2_yy_yzzz_yyz_0, ta2_yy_yzzz_yzz_0, ta2_yy_yzzz_zzz_0, ta2_yy_zzz_xxx_0, ta2_yy_zzz_xxx_1, ta2_yy_zzz_xxz_0, ta2_yy_zzz_xxz_1, ta2_yy_zzz_xyz_0, ta2_yy_zzz_xyz_1, ta2_yy_zzz_xz_0, ta2_yy_zzz_xz_1, ta2_yy_zzz_xzz_0, ta2_yy_zzz_xzz_1, ta2_yy_zzz_yyz_0, ta2_yy_zzz_yyz_1, ta2_yy_zzz_yz_0, ta2_yy_zzz_yz_1, ta2_yy_zzz_yzz_0, ta2_yy_zzz_yzz_1, ta2_yy_zzz_zz_0, ta2_yy_zzz_zz_1, ta2_yy_zzz_zzz_0, ta2_yy_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzzz_xxx_0[i] = 2.0 * ta1_y_zzz_xxx_1[i] + ta2_yy_zzz_xxx_0[i] * pa_y[i] - ta2_yy_zzz_xxx_1[i] * pc_y[i];

        ta2_yy_yzzz_xxy_0[i] = 2.0 * ta2_yy_yz_xxy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xxy_1[i] * fe_0 + ta2_yy_yzz_xxy_0[i] * pa_z[i] - ta2_yy_yzz_xxy_1[i] * pc_z[i];

        ta2_yy_yzzz_xxz_0[i] = 2.0 * ta1_y_zzz_xxz_1[i] + ta2_yy_zzz_xxz_0[i] * pa_y[i] - ta2_yy_zzz_xxz_1[i] * pc_y[i];

        ta2_yy_yzzz_xyy_0[i] = 2.0 * ta2_yy_yz_xyy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xyy_1[i] * fe_0 + ta2_yy_yzz_xyy_0[i] * pa_z[i] - ta2_yy_yzz_xyy_1[i] * pc_z[i];

        ta2_yy_yzzz_xyz_0[i] = ta2_yy_zzz_xz_0[i] * fe_0 - ta2_yy_zzz_xz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyz_1[i] + ta2_yy_zzz_xyz_0[i] * pa_y[i] - ta2_yy_zzz_xyz_1[i] * pc_y[i];

        ta2_yy_yzzz_xzz_0[i] = 2.0 * ta1_y_zzz_xzz_1[i] + ta2_yy_zzz_xzz_0[i] * pa_y[i] - ta2_yy_zzz_xzz_1[i] * pc_y[i];

        ta2_yy_yzzz_yyy_0[i] = 2.0 * ta2_yy_yz_yyy_0[i] * fe_0 - 2.0 * ta2_yy_yz_yyy_1[i] * fe_0 + ta2_yy_yzz_yyy_0[i] * pa_z[i] - ta2_yy_yzz_yyy_1[i] * pc_z[i];

        ta2_yy_yzzz_yyz_0[i] = 2.0 * ta2_yy_zzz_yz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_yz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyz_1[i] + ta2_yy_zzz_yyz_0[i] * pa_y[i] - ta2_yy_zzz_yyz_1[i] * pc_y[i];

        ta2_yy_yzzz_yzz_0[i] = ta2_yy_zzz_zz_0[i] * fe_0 - ta2_yy_zzz_zz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yzz_1[i] + ta2_yy_zzz_yzz_0[i] * pa_y[i] - ta2_yy_zzz_yzz_1[i] * pc_y[i];

        ta2_yy_yzzz_zzz_0[i] = 2.0 * ta1_y_zzz_zzz_1[i] + ta2_yy_zzz_zzz_0[i] * pa_y[i] - ta2_yy_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 590-600 components of targeted buffer : GF

    auto ta2_yy_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 590);

    auto ta2_yy_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 591);

    auto ta2_yy_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 592);

    auto ta2_yy_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 593);

    auto ta2_yy_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 594);

    auto ta2_yy_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 595);

    auto ta2_yy_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 596);

    auto ta2_yy_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 597);

    auto ta2_yy_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 598);

    auto ta2_yy_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 599);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_zz_xxx_0, ta2_yy_zz_xxx_1, ta2_yy_zz_xxy_0, ta2_yy_zz_xxy_1, ta2_yy_zz_xxz_0, ta2_yy_zz_xxz_1, ta2_yy_zz_xyy_0, ta2_yy_zz_xyy_1, ta2_yy_zz_xyz_0, ta2_yy_zz_xyz_1, ta2_yy_zz_xzz_0, ta2_yy_zz_xzz_1, ta2_yy_zz_yyy_0, ta2_yy_zz_yyy_1, ta2_yy_zz_yyz_0, ta2_yy_zz_yyz_1, ta2_yy_zz_yzz_0, ta2_yy_zz_yzz_1, ta2_yy_zz_zzz_0, ta2_yy_zz_zzz_1, ta2_yy_zzz_xx_0, ta2_yy_zzz_xx_1, ta2_yy_zzz_xxx_0, ta2_yy_zzz_xxx_1, ta2_yy_zzz_xxy_0, ta2_yy_zzz_xxy_1, ta2_yy_zzz_xxz_0, ta2_yy_zzz_xxz_1, ta2_yy_zzz_xy_0, ta2_yy_zzz_xy_1, ta2_yy_zzz_xyy_0, ta2_yy_zzz_xyy_1, ta2_yy_zzz_xyz_0, ta2_yy_zzz_xyz_1, ta2_yy_zzz_xz_0, ta2_yy_zzz_xz_1, ta2_yy_zzz_xzz_0, ta2_yy_zzz_xzz_1, ta2_yy_zzz_yy_0, ta2_yy_zzz_yy_1, ta2_yy_zzz_yyy_0, ta2_yy_zzz_yyy_1, ta2_yy_zzz_yyz_0, ta2_yy_zzz_yyz_1, ta2_yy_zzz_yz_0, ta2_yy_zzz_yz_1, ta2_yy_zzz_yzz_0, ta2_yy_zzz_yzz_1, ta2_yy_zzz_zz_0, ta2_yy_zzz_zz_1, ta2_yy_zzz_zzz_0, ta2_yy_zzz_zzz_1, ta2_yy_zzzz_xxx_0, ta2_yy_zzzz_xxy_0, ta2_yy_zzzz_xxz_0, ta2_yy_zzzz_xyy_0, ta2_yy_zzzz_xyz_0, ta2_yy_zzzz_xzz_0, ta2_yy_zzzz_yyy_0, ta2_yy_zzzz_yyz_0, ta2_yy_zzzz_yzz_0, ta2_yy_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzzz_xxx_0[i] = 3.0 * ta2_yy_zz_xxx_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxx_1[i] * fe_0 + ta2_yy_zzz_xxx_0[i] * pa_z[i] - ta2_yy_zzz_xxx_1[i] * pc_z[i];

        ta2_yy_zzzz_xxy_0[i] = 3.0 * ta2_yy_zz_xxy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxy_1[i] * fe_0 + ta2_yy_zzz_xxy_0[i] * pa_z[i] - ta2_yy_zzz_xxy_1[i] * pc_z[i];

        ta2_yy_zzzz_xxz_0[i] = 3.0 * ta2_yy_zz_xxz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxz_1[i] * fe_0 + ta2_yy_zzz_xx_0[i] * fe_0 - ta2_yy_zzz_xx_1[i] * fe_0 + ta2_yy_zzz_xxz_0[i] * pa_z[i] - ta2_yy_zzz_xxz_1[i] * pc_z[i];

        ta2_yy_zzzz_xyy_0[i] = 3.0 * ta2_yy_zz_xyy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xyy_1[i] * fe_0 + ta2_yy_zzz_xyy_0[i] * pa_z[i] - ta2_yy_zzz_xyy_1[i] * pc_z[i];

        ta2_yy_zzzz_xyz_0[i] = 3.0 * ta2_yy_zz_xyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xyz_1[i] * fe_0 + ta2_yy_zzz_xy_0[i] * fe_0 - ta2_yy_zzz_xy_1[i] * fe_0 + ta2_yy_zzz_xyz_0[i] * pa_z[i] - ta2_yy_zzz_xyz_1[i] * pc_z[i];

        ta2_yy_zzzz_xzz_0[i] = 3.0 * ta2_yy_zz_xzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xzz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_xz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xz_1[i] * fe_0 + ta2_yy_zzz_xzz_0[i] * pa_z[i] - ta2_yy_zzz_xzz_1[i] * pc_z[i];

        ta2_yy_zzzz_yyy_0[i] = 3.0 * ta2_yy_zz_yyy_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyy_1[i] * fe_0 + ta2_yy_zzz_yyy_0[i] * pa_z[i] - ta2_yy_zzz_yyy_1[i] * pc_z[i];

        ta2_yy_zzzz_yyz_0[i] = 3.0 * ta2_yy_zz_yyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyz_1[i] * fe_0 + ta2_yy_zzz_yy_0[i] * fe_0 - ta2_yy_zzz_yy_1[i] * fe_0 + ta2_yy_zzz_yyz_0[i] * pa_z[i] - ta2_yy_zzz_yyz_1[i] * pc_z[i];

        ta2_yy_zzzz_yzz_0[i] = 3.0 * ta2_yy_zz_yzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yzz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_yz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_yz_1[i] * fe_0 + ta2_yy_zzz_yzz_0[i] * pa_z[i] - ta2_yy_zzz_yzz_1[i] * pc_z[i];

        ta2_yy_zzzz_zzz_0[i] = 3.0 * ta2_yy_zz_zzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_zzz_1[i] * fe_0 + 3.0 * ta2_yy_zzz_zz_0[i] * fe_0 - 3.0 * ta2_yy_zzz_zz_1[i] * fe_0 + ta2_yy_zzz_zzz_0[i] * pa_z[i] - ta2_yy_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 600-610 components of targeted buffer : GF

    auto ta2_yz_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 600);

    auto ta2_yz_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 601);

    auto ta2_yz_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 602);

    auto ta2_yz_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 603);

    auto ta2_yz_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 604);

    auto ta2_yz_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 605);

    auto ta2_yz_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 606);

    auto ta2_yz_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 607);

    auto ta2_yz_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 608);

    auto ta2_yz_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 609);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xx_xyz_0, ta2_yz_xx_xyz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xx_yyy_0, ta2_yz_xx_yyy_1, ta2_yz_xx_yyz_0, ta2_yz_xx_yyz_1, ta2_yz_xx_yzz_0, ta2_yz_xx_yzz_1, ta2_yz_xx_zzz_0, ta2_yz_xx_zzz_1, ta2_yz_xxx_xx_0, ta2_yz_xxx_xx_1, ta2_yz_xxx_xxx_0, ta2_yz_xxx_xxx_1, ta2_yz_xxx_xxy_0, ta2_yz_xxx_xxy_1, ta2_yz_xxx_xxz_0, ta2_yz_xxx_xxz_1, ta2_yz_xxx_xy_0, ta2_yz_xxx_xy_1, ta2_yz_xxx_xyy_0, ta2_yz_xxx_xyy_1, ta2_yz_xxx_xyz_0, ta2_yz_xxx_xyz_1, ta2_yz_xxx_xz_0, ta2_yz_xxx_xz_1, ta2_yz_xxx_xzz_0, ta2_yz_xxx_xzz_1, ta2_yz_xxx_yy_0, ta2_yz_xxx_yy_1, ta2_yz_xxx_yyy_0, ta2_yz_xxx_yyy_1, ta2_yz_xxx_yyz_0, ta2_yz_xxx_yyz_1, ta2_yz_xxx_yz_0, ta2_yz_xxx_yz_1, ta2_yz_xxx_yzz_0, ta2_yz_xxx_yzz_1, ta2_yz_xxx_zz_0, ta2_yz_xxx_zz_1, ta2_yz_xxx_zzz_0, ta2_yz_xxx_zzz_1, ta2_yz_xxxx_xxx_0, ta2_yz_xxxx_xxy_0, ta2_yz_xxxx_xxz_0, ta2_yz_xxxx_xyy_0, ta2_yz_xxxx_xyz_0, ta2_yz_xxxx_xzz_0, ta2_yz_xxxx_yyy_0, ta2_yz_xxxx_yyz_0, ta2_yz_xxxx_yzz_0, ta2_yz_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxx_xxx_0[i] = 3.0 * ta2_yz_xx_xxx_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxx_1[i] * fe_0 + 3.0 * ta2_yz_xxx_xx_0[i] * fe_0 - 3.0 * ta2_yz_xxx_xx_1[i] * fe_0 + ta2_yz_xxx_xxx_0[i] * pa_x[i] - ta2_yz_xxx_xxx_1[i] * pc_x[i];

        ta2_yz_xxxx_xxy_0[i] = 3.0 * ta2_yz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxy_1[i] * fe_0 + 2.0 * ta2_yz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xy_1[i] * fe_0 + ta2_yz_xxx_xxy_0[i] * pa_x[i] - ta2_yz_xxx_xxy_1[i] * pc_x[i];

        ta2_yz_xxxx_xxz_0[i] = 3.0 * ta2_yz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxz_1[i] * fe_0 + 2.0 * ta2_yz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xz_1[i] * fe_0 + ta2_yz_xxx_xxz_0[i] * pa_x[i] - ta2_yz_xxx_xxz_1[i] * pc_x[i];

        ta2_yz_xxxx_xyy_0[i] = 3.0 * ta2_yz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyy_1[i] * fe_0 + ta2_yz_xxx_yy_0[i] * fe_0 - ta2_yz_xxx_yy_1[i] * fe_0 + ta2_yz_xxx_xyy_0[i] * pa_x[i] - ta2_yz_xxx_xyy_1[i] * pc_x[i];

        ta2_yz_xxxx_xyz_0[i] = 3.0 * ta2_yz_xx_xyz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyz_1[i] * fe_0 + ta2_yz_xxx_yz_0[i] * fe_0 - ta2_yz_xxx_yz_1[i] * fe_0 + ta2_yz_xxx_xyz_0[i] * pa_x[i] - ta2_yz_xxx_xyz_1[i] * pc_x[i];

        ta2_yz_xxxx_xzz_0[i] = 3.0 * ta2_yz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xzz_1[i] * fe_0 + ta2_yz_xxx_zz_0[i] * fe_0 - ta2_yz_xxx_zz_1[i] * fe_0 + ta2_yz_xxx_xzz_0[i] * pa_x[i] - ta2_yz_xxx_xzz_1[i] * pc_x[i];

        ta2_yz_xxxx_yyy_0[i] = 3.0 * ta2_yz_xx_yyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_yyy_1[i] * fe_0 + ta2_yz_xxx_yyy_0[i] * pa_x[i] - ta2_yz_xxx_yyy_1[i] * pc_x[i];

        ta2_yz_xxxx_yyz_0[i] = 3.0 * ta2_yz_xx_yyz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yyz_1[i] * fe_0 + ta2_yz_xxx_yyz_0[i] * pa_x[i] - ta2_yz_xxx_yyz_1[i] * pc_x[i];

        ta2_yz_xxxx_yzz_0[i] = 3.0 * ta2_yz_xx_yzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yzz_1[i] * fe_0 + ta2_yz_xxx_yzz_0[i] * pa_x[i] - ta2_yz_xxx_yzz_1[i] * pc_x[i];

        ta2_yz_xxxx_zzz_0[i] = 3.0 * ta2_yz_xx_zzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_zzz_1[i] * fe_0 + ta2_yz_xxx_zzz_0[i] * pa_x[i] - ta2_yz_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 610-620 components of targeted buffer : GF

    auto ta2_yz_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 610);

    auto ta2_yz_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 611);

    auto ta2_yz_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 612);

    auto ta2_yz_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 613);

    auto ta2_yz_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 614);

    auto ta2_yz_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 615);

    auto ta2_yz_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 616);

    auto ta2_yz_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 617);

    auto ta2_yz_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 618);

    auto ta2_yz_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 619);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxx_xxx_1, ta1_z_xxx_xxy_1, ta1_z_xxx_xxz_1, ta1_z_xxx_xyy_1, ta1_z_xxx_xyz_1, ta1_z_xxx_xzz_1, ta1_z_xxx_zzz_1, ta2_yz_xxx_xx_0, ta2_yz_xxx_xx_1, ta2_yz_xxx_xxx_0, ta2_yz_xxx_xxx_1, ta2_yz_xxx_xxy_0, ta2_yz_xxx_xxy_1, ta2_yz_xxx_xxz_0, ta2_yz_xxx_xxz_1, ta2_yz_xxx_xy_0, ta2_yz_xxx_xy_1, ta2_yz_xxx_xyy_0, ta2_yz_xxx_xyy_1, ta2_yz_xxx_xyz_0, ta2_yz_xxx_xyz_1, ta2_yz_xxx_xz_0, ta2_yz_xxx_xz_1, ta2_yz_xxx_xzz_0, ta2_yz_xxx_xzz_1, ta2_yz_xxx_zzz_0, ta2_yz_xxx_zzz_1, ta2_yz_xxxy_xxx_0, ta2_yz_xxxy_xxy_0, ta2_yz_xxxy_xxz_0, ta2_yz_xxxy_xyy_0, ta2_yz_xxxy_xyz_0, ta2_yz_xxxy_xzz_0, ta2_yz_xxxy_yyy_0, ta2_yz_xxxy_yyz_0, ta2_yz_xxxy_yzz_0, ta2_yz_xxxy_zzz_0, ta2_yz_xxy_yyy_0, ta2_yz_xxy_yyy_1, ta2_yz_xxy_yyz_0, ta2_yz_xxy_yyz_1, ta2_yz_xxy_yzz_0, ta2_yz_xxy_yzz_1, ta2_yz_xy_yyy_0, ta2_yz_xy_yyy_1, ta2_yz_xy_yyz_0, ta2_yz_xy_yyz_1, ta2_yz_xy_yzz_0, ta2_yz_xy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxy_xxx_0[i] = ta1_z_xxx_xxx_1[i] + ta2_yz_xxx_xxx_0[i] * pa_y[i] - ta2_yz_xxx_xxx_1[i] * pc_y[i];

        ta2_yz_xxxy_xxy_0[i] = ta2_yz_xxx_xx_0[i] * fe_0 - ta2_yz_xxx_xx_1[i] * fe_0 + ta1_z_xxx_xxy_1[i] + ta2_yz_xxx_xxy_0[i] * pa_y[i] - ta2_yz_xxx_xxy_1[i] * pc_y[i];

        ta2_yz_xxxy_xxz_0[i] = ta1_z_xxx_xxz_1[i] + ta2_yz_xxx_xxz_0[i] * pa_y[i] - ta2_yz_xxx_xxz_1[i] * pc_y[i];

        ta2_yz_xxxy_xyy_0[i] = 2.0 * ta2_yz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xy_1[i] * fe_0 + ta1_z_xxx_xyy_1[i] + ta2_yz_xxx_xyy_0[i] * pa_y[i] - ta2_yz_xxx_xyy_1[i] * pc_y[i];

        ta2_yz_xxxy_xyz_0[i] = ta2_yz_xxx_xz_0[i] * fe_0 - ta2_yz_xxx_xz_1[i] * fe_0 + ta1_z_xxx_xyz_1[i] + ta2_yz_xxx_xyz_0[i] * pa_y[i] - ta2_yz_xxx_xyz_1[i] * pc_y[i];

        ta2_yz_xxxy_xzz_0[i] = ta1_z_xxx_xzz_1[i] + ta2_yz_xxx_xzz_0[i] * pa_y[i] - ta2_yz_xxx_xzz_1[i] * pc_y[i];

        ta2_yz_xxxy_yyy_0[i] = 2.0 * ta2_yz_xy_yyy_0[i] * fe_0 - 2.0 * ta2_yz_xy_yyy_1[i] * fe_0 + ta2_yz_xxy_yyy_0[i] * pa_x[i] - ta2_yz_xxy_yyy_1[i] * pc_x[i];

        ta2_yz_xxxy_yyz_0[i] = 2.0 * ta2_yz_xy_yyz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yyz_1[i] * fe_0 + ta2_yz_xxy_yyz_0[i] * pa_x[i] - ta2_yz_xxy_yyz_1[i] * pc_x[i];

        ta2_yz_xxxy_yzz_0[i] = 2.0 * ta2_yz_xy_yzz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yzz_1[i] * fe_0 + ta2_yz_xxy_yzz_0[i] * pa_x[i] - ta2_yz_xxy_yzz_1[i] * pc_x[i];

        ta2_yz_xxxy_zzz_0[i] = ta1_z_xxx_zzz_1[i] + ta2_yz_xxx_zzz_0[i] * pa_y[i] - ta2_yz_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 620-630 components of targeted buffer : GF

    auto ta2_yz_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 620);

    auto ta2_yz_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 621);

    auto ta2_yz_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 622);

    auto ta2_yz_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 623);

    auto ta2_yz_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 624);

    auto ta2_yz_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 625);

    auto ta2_yz_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 626);

    auto ta2_yz_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 627);

    auto ta2_yz_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 628);

    auto ta2_yz_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 629);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxx_xxx_1, ta1_y_xxx_xxy_1, ta1_y_xxx_xxz_1, ta1_y_xxx_xyy_1, ta1_y_xxx_xyz_1, ta1_y_xxx_xzz_1, ta1_y_xxx_yyy_1, ta2_yz_xxx_xx_0, ta2_yz_xxx_xx_1, ta2_yz_xxx_xxx_0, ta2_yz_xxx_xxx_1, ta2_yz_xxx_xxy_0, ta2_yz_xxx_xxy_1, ta2_yz_xxx_xxz_0, ta2_yz_xxx_xxz_1, ta2_yz_xxx_xy_0, ta2_yz_xxx_xy_1, ta2_yz_xxx_xyy_0, ta2_yz_xxx_xyy_1, ta2_yz_xxx_xyz_0, ta2_yz_xxx_xyz_1, ta2_yz_xxx_xz_0, ta2_yz_xxx_xz_1, ta2_yz_xxx_xzz_0, ta2_yz_xxx_xzz_1, ta2_yz_xxx_yyy_0, ta2_yz_xxx_yyy_1, ta2_yz_xxxz_xxx_0, ta2_yz_xxxz_xxy_0, ta2_yz_xxxz_xxz_0, ta2_yz_xxxz_xyy_0, ta2_yz_xxxz_xyz_0, ta2_yz_xxxz_xzz_0, ta2_yz_xxxz_yyy_0, ta2_yz_xxxz_yyz_0, ta2_yz_xxxz_yzz_0, ta2_yz_xxxz_zzz_0, ta2_yz_xxz_yyz_0, ta2_yz_xxz_yyz_1, ta2_yz_xxz_yzz_0, ta2_yz_xxz_yzz_1, ta2_yz_xxz_zzz_0, ta2_yz_xxz_zzz_1, ta2_yz_xz_yyz_0, ta2_yz_xz_yyz_1, ta2_yz_xz_yzz_0, ta2_yz_xz_yzz_1, ta2_yz_xz_zzz_0, ta2_yz_xz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxz_xxx_0[i] = ta1_y_xxx_xxx_1[i] + ta2_yz_xxx_xxx_0[i] * pa_z[i] - ta2_yz_xxx_xxx_1[i] * pc_z[i];

        ta2_yz_xxxz_xxy_0[i] = ta1_y_xxx_xxy_1[i] + ta2_yz_xxx_xxy_0[i] * pa_z[i] - ta2_yz_xxx_xxy_1[i] * pc_z[i];

        ta2_yz_xxxz_xxz_0[i] = ta2_yz_xxx_xx_0[i] * fe_0 - ta2_yz_xxx_xx_1[i] * fe_0 + ta1_y_xxx_xxz_1[i] + ta2_yz_xxx_xxz_0[i] * pa_z[i] - ta2_yz_xxx_xxz_1[i] * pc_z[i];

        ta2_yz_xxxz_xyy_0[i] = ta1_y_xxx_xyy_1[i] + ta2_yz_xxx_xyy_0[i] * pa_z[i] - ta2_yz_xxx_xyy_1[i] * pc_z[i];

        ta2_yz_xxxz_xyz_0[i] = ta2_yz_xxx_xy_0[i] * fe_0 - ta2_yz_xxx_xy_1[i] * fe_0 + ta1_y_xxx_xyz_1[i] + ta2_yz_xxx_xyz_0[i] * pa_z[i] - ta2_yz_xxx_xyz_1[i] * pc_z[i];

        ta2_yz_xxxz_xzz_0[i] = 2.0 * ta2_yz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xz_1[i] * fe_0 + ta1_y_xxx_xzz_1[i] + ta2_yz_xxx_xzz_0[i] * pa_z[i] - ta2_yz_xxx_xzz_1[i] * pc_z[i];

        ta2_yz_xxxz_yyy_0[i] = ta1_y_xxx_yyy_1[i] + ta2_yz_xxx_yyy_0[i] * pa_z[i] - ta2_yz_xxx_yyy_1[i] * pc_z[i];

        ta2_yz_xxxz_yyz_0[i] = 2.0 * ta2_yz_xz_yyz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yyz_1[i] * fe_0 + ta2_yz_xxz_yyz_0[i] * pa_x[i] - ta2_yz_xxz_yyz_1[i] * pc_x[i];

        ta2_yz_xxxz_yzz_0[i] = 2.0 * ta2_yz_xz_yzz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yzz_1[i] * fe_0 + ta2_yz_xxz_yzz_0[i] * pa_x[i] - ta2_yz_xxz_yzz_1[i] * pc_x[i];

        ta2_yz_xxxz_zzz_0[i] = 2.0 * ta2_yz_xz_zzz_0[i] * fe_0 - 2.0 * ta2_yz_xz_zzz_1[i] * fe_0 + ta2_yz_xxz_zzz_0[i] * pa_x[i] - ta2_yz_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 630-640 components of targeted buffer : GF

    auto ta2_yz_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 630);

    auto ta2_yz_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 631);

    auto ta2_yz_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 632);

    auto ta2_yz_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 633);

    auto ta2_yz_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 634);

    auto ta2_yz_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 635);

    auto ta2_yz_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 636);

    auto ta2_yz_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 637);

    auto ta2_yz_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 638);

    auto ta2_yz_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 639);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxy_xxx_1, ta1_z_xxy_xxz_1, ta1_z_xxy_xzz_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxz_0, ta2_yz_xx_xxz_1, ta2_yz_xx_xzz_0, ta2_yz_xx_xzz_1, ta2_yz_xxy_xxx_0, ta2_yz_xxy_xxx_1, ta2_yz_xxy_xxz_0, ta2_yz_xxy_xxz_1, ta2_yz_xxy_xzz_0, ta2_yz_xxy_xzz_1, ta2_yz_xxyy_xxx_0, ta2_yz_xxyy_xxy_0, ta2_yz_xxyy_xxz_0, ta2_yz_xxyy_xyy_0, ta2_yz_xxyy_xyz_0, ta2_yz_xxyy_xzz_0, ta2_yz_xxyy_yyy_0, ta2_yz_xxyy_yyz_0, ta2_yz_xxyy_yzz_0, ta2_yz_xxyy_zzz_0, ta2_yz_xyy_xxy_0, ta2_yz_xyy_xxy_1, ta2_yz_xyy_xy_0, ta2_yz_xyy_xy_1, ta2_yz_xyy_xyy_0, ta2_yz_xyy_xyy_1, ta2_yz_xyy_xyz_0, ta2_yz_xyy_xyz_1, ta2_yz_xyy_yy_0, ta2_yz_xyy_yy_1, ta2_yz_xyy_yyy_0, ta2_yz_xyy_yyy_1, ta2_yz_xyy_yyz_0, ta2_yz_xyy_yyz_1, ta2_yz_xyy_yz_0, ta2_yz_xyy_yz_1, ta2_yz_xyy_yzz_0, ta2_yz_xyy_yzz_1, ta2_yz_xyy_zzz_0, ta2_yz_xyy_zzz_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyy_xxx_0[i] = ta2_yz_xx_xxx_0[i] * fe_0 - ta2_yz_xx_xxx_1[i] * fe_0 + ta1_z_xxy_xxx_1[i] + ta2_yz_xxy_xxx_0[i] * pa_y[i] - ta2_yz_xxy_xxx_1[i] * pc_y[i];

        ta2_yz_xxyy_xxy_0[i] = ta2_yz_yy_xxy_0[i] * fe_0 - ta2_yz_yy_xxy_1[i] * fe_0 + 2.0 * ta2_yz_xyy_xy_0[i] * fe_0 - 2.0 * ta2_yz_xyy_xy_1[i] * fe_0 + ta2_yz_xyy_xxy_0[i] * pa_x[i] - ta2_yz_xyy_xxy_1[i] * pc_x[i];

        ta2_yz_xxyy_xxz_0[i] = ta2_yz_xx_xxz_0[i] * fe_0 - ta2_yz_xx_xxz_1[i] * fe_0 + ta1_z_xxy_xxz_1[i] + ta2_yz_xxy_xxz_0[i] * pa_y[i] - ta2_yz_xxy_xxz_1[i] * pc_y[i];

        ta2_yz_xxyy_xyy_0[i] = ta2_yz_yy_xyy_0[i] * fe_0 - ta2_yz_yy_xyy_1[i] * fe_0 + ta2_yz_xyy_yy_0[i] * fe_0 - ta2_yz_xyy_yy_1[i] * fe_0 + ta2_yz_xyy_xyy_0[i] * pa_x[i] - ta2_yz_xyy_xyy_1[i] * pc_x[i];

        ta2_yz_xxyy_xyz_0[i] = ta2_yz_yy_xyz_0[i] * fe_0 - ta2_yz_yy_xyz_1[i] * fe_0 + ta2_yz_xyy_yz_0[i] * fe_0 - ta2_yz_xyy_yz_1[i] * fe_0 + ta2_yz_xyy_xyz_0[i] * pa_x[i] - ta2_yz_xyy_xyz_1[i] * pc_x[i];

        ta2_yz_xxyy_xzz_0[i] = ta2_yz_xx_xzz_0[i] * fe_0 - ta2_yz_xx_xzz_1[i] * fe_0 + ta1_z_xxy_xzz_1[i] + ta2_yz_xxy_xzz_0[i] * pa_y[i] - ta2_yz_xxy_xzz_1[i] * pc_y[i];

        ta2_yz_xxyy_yyy_0[i] = ta2_yz_yy_yyy_0[i] * fe_0 - ta2_yz_yy_yyy_1[i] * fe_0 + ta2_yz_xyy_yyy_0[i] * pa_x[i] - ta2_yz_xyy_yyy_1[i] * pc_x[i];

        ta2_yz_xxyy_yyz_0[i] = ta2_yz_yy_yyz_0[i] * fe_0 - ta2_yz_yy_yyz_1[i] * fe_0 + ta2_yz_xyy_yyz_0[i] * pa_x[i] - ta2_yz_xyy_yyz_1[i] * pc_x[i];

        ta2_yz_xxyy_yzz_0[i] = ta2_yz_yy_yzz_0[i] * fe_0 - ta2_yz_yy_yzz_1[i] * fe_0 + ta2_yz_xyy_yzz_0[i] * pa_x[i] - ta2_yz_xyy_yzz_1[i] * pc_x[i];

        ta2_yz_xxyy_zzz_0[i] = ta2_yz_yy_zzz_0[i] * fe_0 - ta2_yz_yy_zzz_1[i] * fe_0 + ta2_yz_xyy_zzz_0[i] * pa_x[i] - ta2_yz_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 640-650 components of targeted buffer : GF

    auto ta2_yz_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 640);

    auto ta2_yz_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 641);

    auto ta2_yz_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 642);

    auto ta2_yz_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 643);

    auto ta2_yz_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 644);

    auto ta2_yz_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 645);

    auto ta2_yz_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 646);

    auto ta2_yz_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 647);

    auto ta2_yz_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 648);

    auto ta2_yz_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 649);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxy_xxy_1, ta1_y_xxy_xyy_1, ta1_y_xxy_yyy_1, ta1_z_xxz_xxx_1, ta1_z_xxz_xxz_1, ta1_z_xxz_xyz_1, ta1_z_xxz_xzz_1, ta1_z_xxz_zzz_1, ta2_yz_xxy_xxy_0, ta2_yz_xxy_xxy_1, ta2_yz_xxy_xyy_0, ta2_yz_xxy_xyy_1, ta2_yz_xxy_yyy_0, ta2_yz_xxy_yyy_1, ta2_yz_xxyz_xxx_0, ta2_yz_xxyz_xxy_0, ta2_yz_xxyz_xxz_0, ta2_yz_xxyz_xyy_0, ta2_yz_xxyz_xyz_0, ta2_yz_xxyz_xzz_0, ta2_yz_xxyz_yyy_0, ta2_yz_xxyz_yyz_0, ta2_yz_xxyz_yzz_0, ta2_yz_xxyz_zzz_0, ta2_yz_xxz_xxx_0, ta2_yz_xxz_xxx_1, ta2_yz_xxz_xxz_0, ta2_yz_xxz_xxz_1, ta2_yz_xxz_xyz_0, ta2_yz_xxz_xyz_1, ta2_yz_xxz_xz_0, ta2_yz_xxz_xz_1, ta2_yz_xxz_xzz_0, ta2_yz_xxz_xzz_1, ta2_yz_xxz_zzz_0, ta2_yz_xxz_zzz_1, ta2_yz_xyz_yyz_0, ta2_yz_xyz_yyz_1, ta2_yz_xyz_yzz_0, ta2_yz_xyz_yzz_1, ta2_yz_yz_yyz_0, ta2_yz_yz_yyz_1, ta2_yz_yz_yzz_0, ta2_yz_yz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyz_xxx_0[i] = ta1_z_xxz_xxx_1[i] + ta2_yz_xxz_xxx_0[i] * pa_y[i] - ta2_yz_xxz_xxx_1[i] * pc_y[i];

        ta2_yz_xxyz_xxy_0[i] = ta1_y_xxy_xxy_1[i] + ta2_yz_xxy_xxy_0[i] * pa_z[i] - ta2_yz_xxy_xxy_1[i] * pc_z[i];

        ta2_yz_xxyz_xxz_0[i] = ta1_z_xxz_xxz_1[i] + ta2_yz_xxz_xxz_0[i] * pa_y[i] - ta2_yz_xxz_xxz_1[i] * pc_y[i];

        ta2_yz_xxyz_xyy_0[i] = ta1_y_xxy_xyy_1[i] + ta2_yz_xxy_xyy_0[i] * pa_z[i] - ta2_yz_xxy_xyy_1[i] * pc_z[i];

        ta2_yz_xxyz_xyz_0[i] = ta2_yz_xxz_xz_0[i] * fe_0 - ta2_yz_xxz_xz_1[i] * fe_0 + ta1_z_xxz_xyz_1[i] + ta2_yz_xxz_xyz_0[i] * pa_y[i] - ta2_yz_xxz_xyz_1[i] * pc_y[i];

        ta2_yz_xxyz_xzz_0[i] = ta1_z_xxz_xzz_1[i] + ta2_yz_xxz_xzz_0[i] * pa_y[i] - ta2_yz_xxz_xzz_1[i] * pc_y[i];

        ta2_yz_xxyz_yyy_0[i] = ta1_y_xxy_yyy_1[i] + ta2_yz_xxy_yyy_0[i] * pa_z[i] - ta2_yz_xxy_yyy_1[i] * pc_z[i];

        ta2_yz_xxyz_yyz_0[i] = ta2_yz_yz_yyz_0[i] * fe_0 - ta2_yz_yz_yyz_1[i] * fe_0 + ta2_yz_xyz_yyz_0[i] * pa_x[i] - ta2_yz_xyz_yyz_1[i] * pc_x[i];

        ta2_yz_xxyz_yzz_0[i] = ta2_yz_yz_yzz_0[i] * fe_0 - ta2_yz_yz_yzz_1[i] * fe_0 + ta2_yz_xyz_yzz_0[i] * pa_x[i] - ta2_yz_xyz_yzz_1[i] * pc_x[i];

        ta2_yz_xxyz_zzz_0[i] = ta1_z_xxz_zzz_1[i] + ta2_yz_xxz_zzz_0[i] * pa_y[i] - ta2_yz_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 650-660 components of targeted buffer : GF

    auto ta2_yz_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 650);

    auto ta2_yz_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 651);

    auto ta2_yz_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 652);

    auto ta2_yz_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 653);

    auto ta2_yz_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 654);

    auto ta2_yz_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 655);

    auto ta2_yz_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 656);

    auto ta2_yz_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 657);

    auto ta2_yz_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 658);

    auto ta2_yz_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 659);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxz_xxx_1, ta1_y_xxz_xxy_1, ta1_y_xxz_xyy_1, ta2_yz_xx_xxx_0, ta2_yz_xx_xxx_1, ta2_yz_xx_xxy_0, ta2_yz_xx_xxy_1, ta2_yz_xx_xyy_0, ta2_yz_xx_xyy_1, ta2_yz_xxz_xxx_0, ta2_yz_xxz_xxx_1, ta2_yz_xxz_xxy_0, ta2_yz_xxz_xxy_1, ta2_yz_xxz_xyy_0, ta2_yz_xxz_xyy_1, ta2_yz_xxzz_xxx_0, ta2_yz_xxzz_xxy_0, ta2_yz_xxzz_xxz_0, ta2_yz_xxzz_xyy_0, ta2_yz_xxzz_xyz_0, ta2_yz_xxzz_xzz_0, ta2_yz_xxzz_yyy_0, ta2_yz_xxzz_yyz_0, ta2_yz_xxzz_yzz_0, ta2_yz_xxzz_zzz_0, ta2_yz_xzz_xxz_0, ta2_yz_xzz_xxz_1, ta2_yz_xzz_xyz_0, ta2_yz_xzz_xyz_1, ta2_yz_xzz_xz_0, ta2_yz_xzz_xz_1, ta2_yz_xzz_xzz_0, ta2_yz_xzz_xzz_1, ta2_yz_xzz_yyy_0, ta2_yz_xzz_yyy_1, ta2_yz_xzz_yyz_0, ta2_yz_xzz_yyz_1, ta2_yz_xzz_yz_0, ta2_yz_xzz_yz_1, ta2_yz_xzz_yzz_0, ta2_yz_xzz_yzz_1, ta2_yz_xzz_zz_0, ta2_yz_xzz_zz_1, ta2_yz_xzz_zzz_0, ta2_yz_xzz_zzz_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxzz_xxx_0[i] = ta2_yz_xx_xxx_0[i] * fe_0 - ta2_yz_xx_xxx_1[i] * fe_0 + ta1_y_xxz_xxx_1[i] + ta2_yz_xxz_xxx_0[i] * pa_z[i] - ta2_yz_xxz_xxx_1[i] * pc_z[i];

        ta2_yz_xxzz_xxy_0[i] = ta2_yz_xx_xxy_0[i] * fe_0 - ta2_yz_xx_xxy_1[i] * fe_0 + ta1_y_xxz_xxy_1[i] + ta2_yz_xxz_xxy_0[i] * pa_z[i] - ta2_yz_xxz_xxy_1[i] * pc_z[i];

        ta2_yz_xxzz_xxz_0[i] = ta2_yz_zz_xxz_0[i] * fe_0 - ta2_yz_zz_xxz_1[i] * fe_0 + 2.0 * ta2_yz_xzz_xz_0[i] * fe_0 - 2.0 * ta2_yz_xzz_xz_1[i] * fe_0 + ta2_yz_xzz_xxz_0[i] * pa_x[i] - ta2_yz_xzz_xxz_1[i] * pc_x[i];

        ta2_yz_xxzz_xyy_0[i] = ta2_yz_xx_xyy_0[i] * fe_0 - ta2_yz_xx_xyy_1[i] * fe_0 + ta1_y_xxz_xyy_1[i] + ta2_yz_xxz_xyy_0[i] * pa_z[i] - ta2_yz_xxz_xyy_1[i] * pc_z[i];

        ta2_yz_xxzz_xyz_0[i] = ta2_yz_zz_xyz_0[i] * fe_0 - ta2_yz_zz_xyz_1[i] * fe_0 + ta2_yz_xzz_yz_0[i] * fe_0 - ta2_yz_xzz_yz_1[i] * fe_0 + ta2_yz_xzz_xyz_0[i] * pa_x[i] - ta2_yz_xzz_xyz_1[i] * pc_x[i];

        ta2_yz_xxzz_xzz_0[i] = ta2_yz_zz_xzz_0[i] * fe_0 - ta2_yz_zz_xzz_1[i] * fe_0 + ta2_yz_xzz_zz_0[i] * fe_0 - ta2_yz_xzz_zz_1[i] * fe_0 + ta2_yz_xzz_xzz_0[i] * pa_x[i] - ta2_yz_xzz_xzz_1[i] * pc_x[i];

        ta2_yz_xxzz_yyy_0[i] = ta2_yz_zz_yyy_0[i] * fe_0 - ta2_yz_zz_yyy_1[i] * fe_0 + ta2_yz_xzz_yyy_0[i] * pa_x[i] - ta2_yz_xzz_yyy_1[i] * pc_x[i];

        ta2_yz_xxzz_yyz_0[i] = ta2_yz_zz_yyz_0[i] * fe_0 - ta2_yz_zz_yyz_1[i] * fe_0 + ta2_yz_xzz_yyz_0[i] * pa_x[i] - ta2_yz_xzz_yyz_1[i] * pc_x[i];

        ta2_yz_xxzz_yzz_0[i] = ta2_yz_zz_yzz_0[i] * fe_0 - ta2_yz_zz_yzz_1[i] * fe_0 + ta2_yz_xzz_yzz_0[i] * pa_x[i] - ta2_yz_xzz_yzz_1[i] * pc_x[i];

        ta2_yz_xxzz_zzz_0[i] = ta2_yz_zz_zzz_0[i] * fe_0 - ta2_yz_zz_zzz_1[i] * fe_0 + ta2_yz_xzz_zzz_0[i] * pa_x[i] - ta2_yz_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 660-670 components of targeted buffer : GF

    auto ta2_yz_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 660);

    auto ta2_yz_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 661);

    auto ta2_yz_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 662);

    auto ta2_yz_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 663);

    auto ta2_yz_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 664);

    auto ta2_yz_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 665);

    auto ta2_yz_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 666);

    auto ta2_yz_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 667);

    auto ta2_yz_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 668);

    auto ta2_yz_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 669);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xyyy_xxx_0, ta2_yz_xyyy_xxy_0, ta2_yz_xyyy_xxz_0, ta2_yz_xyyy_xyy_0, ta2_yz_xyyy_xyz_0, ta2_yz_xyyy_xzz_0, ta2_yz_xyyy_yyy_0, ta2_yz_xyyy_yyz_0, ta2_yz_xyyy_yzz_0, ta2_yz_xyyy_zzz_0, ta2_yz_yyy_xx_0, ta2_yz_yyy_xx_1, ta2_yz_yyy_xxx_0, ta2_yz_yyy_xxx_1, ta2_yz_yyy_xxy_0, ta2_yz_yyy_xxy_1, ta2_yz_yyy_xxz_0, ta2_yz_yyy_xxz_1, ta2_yz_yyy_xy_0, ta2_yz_yyy_xy_1, ta2_yz_yyy_xyy_0, ta2_yz_yyy_xyy_1, ta2_yz_yyy_xyz_0, ta2_yz_yyy_xyz_1, ta2_yz_yyy_xz_0, ta2_yz_yyy_xz_1, ta2_yz_yyy_xzz_0, ta2_yz_yyy_xzz_1, ta2_yz_yyy_yy_0, ta2_yz_yyy_yy_1, ta2_yz_yyy_yyy_0, ta2_yz_yyy_yyy_1, ta2_yz_yyy_yyz_0, ta2_yz_yyy_yyz_1, ta2_yz_yyy_yz_0, ta2_yz_yyy_yz_1, ta2_yz_yyy_yzz_0, ta2_yz_yyy_yzz_1, ta2_yz_yyy_zz_0, ta2_yz_yyy_zz_1, ta2_yz_yyy_zzz_0, ta2_yz_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyy_xxx_0[i] = 3.0 * ta2_yz_yyy_xx_0[i] * fe_0 - 3.0 * ta2_yz_yyy_xx_1[i] * fe_0 + ta2_yz_yyy_xxx_0[i] * pa_x[i] - ta2_yz_yyy_xxx_1[i] * pc_x[i];

        ta2_yz_xyyy_xxy_0[i] = 2.0 * ta2_yz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xy_1[i] * fe_0 + ta2_yz_yyy_xxy_0[i] * pa_x[i] - ta2_yz_yyy_xxy_1[i] * pc_x[i];

        ta2_yz_xyyy_xxz_0[i] = 2.0 * ta2_yz_yyy_xz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xz_1[i] * fe_0 + ta2_yz_yyy_xxz_0[i] * pa_x[i] - ta2_yz_yyy_xxz_1[i] * pc_x[i];

        ta2_yz_xyyy_xyy_0[i] = ta2_yz_yyy_yy_0[i] * fe_0 - ta2_yz_yyy_yy_1[i] * fe_0 + ta2_yz_yyy_xyy_0[i] * pa_x[i] - ta2_yz_yyy_xyy_1[i] * pc_x[i];

        ta2_yz_xyyy_xyz_0[i] = ta2_yz_yyy_yz_0[i] * fe_0 - ta2_yz_yyy_yz_1[i] * fe_0 + ta2_yz_yyy_xyz_0[i] * pa_x[i] - ta2_yz_yyy_xyz_1[i] * pc_x[i];

        ta2_yz_xyyy_xzz_0[i] = ta2_yz_yyy_zz_0[i] * fe_0 - ta2_yz_yyy_zz_1[i] * fe_0 + ta2_yz_yyy_xzz_0[i] * pa_x[i] - ta2_yz_yyy_xzz_1[i] * pc_x[i];

        ta2_yz_xyyy_yyy_0[i] = ta2_yz_yyy_yyy_0[i] * pa_x[i] - ta2_yz_yyy_yyy_1[i] * pc_x[i];

        ta2_yz_xyyy_yyz_0[i] = ta2_yz_yyy_yyz_0[i] * pa_x[i] - ta2_yz_yyy_yyz_1[i] * pc_x[i];

        ta2_yz_xyyy_yzz_0[i] = ta2_yz_yyy_yzz_0[i] * pa_x[i] - ta2_yz_yyy_yzz_1[i] * pc_x[i];

        ta2_yz_xyyy_zzz_0[i] = ta2_yz_yyy_zzz_0[i] * pa_x[i] - ta2_yz_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 670-680 components of targeted buffer : GF

    auto ta2_yz_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 670);

    auto ta2_yz_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 671);

    auto ta2_yz_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 672);

    auto ta2_yz_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 673);

    auto ta2_yz_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 674);

    auto ta2_yz_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 675);

    auto ta2_yz_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 676);

    auto ta2_yz_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 677);

    auto ta2_yz_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 678);

    auto ta2_yz_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 679);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xyy_xxx_1, ta1_y_xyy_xxy_1, ta1_y_xyy_xyy_1, ta2_yz_xyy_xxx_0, ta2_yz_xyy_xxx_1, ta2_yz_xyy_xxy_0, ta2_yz_xyy_xxy_1, ta2_yz_xyy_xyy_0, ta2_yz_xyy_xyy_1, ta2_yz_xyyz_xxx_0, ta2_yz_xyyz_xxy_0, ta2_yz_xyyz_xxz_0, ta2_yz_xyyz_xyy_0, ta2_yz_xyyz_xyz_0, ta2_yz_xyyz_xzz_0, ta2_yz_xyyz_yyy_0, ta2_yz_xyyz_yyz_0, ta2_yz_xyyz_yzz_0, ta2_yz_xyyz_zzz_0, ta2_yz_yyz_xxz_0, ta2_yz_yyz_xxz_1, ta2_yz_yyz_xyz_0, ta2_yz_yyz_xyz_1, ta2_yz_yyz_xz_0, ta2_yz_yyz_xz_1, ta2_yz_yyz_xzz_0, ta2_yz_yyz_xzz_1, ta2_yz_yyz_yyy_0, ta2_yz_yyz_yyy_1, ta2_yz_yyz_yyz_0, ta2_yz_yyz_yyz_1, ta2_yz_yyz_yz_0, ta2_yz_yyz_yz_1, ta2_yz_yyz_yzz_0, ta2_yz_yyz_yzz_1, ta2_yz_yyz_zz_0, ta2_yz_yyz_zz_1, ta2_yz_yyz_zzz_0, ta2_yz_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyz_xxx_0[i] = ta1_y_xyy_xxx_1[i] + ta2_yz_xyy_xxx_0[i] * pa_z[i] - ta2_yz_xyy_xxx_1[i] * pc_z[i];

        ta2_yz_xyyz_xxy_0[i] = ta1_y_xyy_xxy_1[i] + ta2_yz_xyy_xxy_0[i] * pa_z[i] - ta2_yz_xyy_xxy_1[i] * pc_z[i];

        ta2_yz_xyyz_xxz_0[i] = 2.0 * ta2_yz_yyz_xz_0[i] * fe_0 - 2.0 * ta2_yz_yyz_xz_1[i] * fe_0 + ta2_yz_yyz_xxz_0[i] * pa_x[i] - ta2_yz_yyz_xxz_1[i] * pc_x[i];

        ta2_yz_xyyz_xyy_0[i] = ta1_y_xyy_xyy_1[i] + ta2_yz_xyy_xyy_0[i] * pa_z[i] - ta2_yz_xyy_xyy_1[i] * pc_z[i];

        ta2_yz_xyyz_xyz_0[i] = ta2_yz_yyz_yz_0[i] * fe_0 - ta2_yz_yyz_yz_1[i] * fe_0 + ta2_yz_yyz_xyz_0[i] * pa_x[i] - ta2_yz_yyz_xyz_1[i] * pc_x[i];

        ta2_yz_xyyz_xzz_0[i] = ta2_yz_yyz_zz_0[i] * fe_0 - ta2_yz_yyz_zz_1[i] * fe_0 + ta2_yz_yyz_xzz_0[i] * pa_x[i] - ta2_yz_yyz_xzz_1[i] * pc_x[i];

        ta2_yz_xyyz_yyy_0[i] = ta2_yz_yyz_yyy_0[i] * pa_x[i] - ta2_yz_yyz_yyy_1[i] * pc_x[i];

        ta2_yz_xyyz_yyz_0[i] = ta2_yz_yyz_yyz_0[i] * pa_x[i] - ta2_yz_yyz_yyz_1[i] * pc_x[i];

        ta2_yz_xyyz_yzz_0[i] = ta2_yz_yyz_yzz_0[i] * pa_x[i] - ta2_yz_yyz_yzz_1[i] * pc_x[i];

        ta2_yz_xyyz_zzz_0[i] = ta2_yz_yyz_zzz_0[i] * pa_x[i] - ta2_yz_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 680-690 components of targeted buffer : GF

    auto ta2_yz_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 680);

    auto ta2_yz_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 681);

    auto ta2_yz_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 682);

    auto ta2_yz_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 683);

    auto ta2_yz_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 684);

    auto ta2_yz_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 685);

    auto ta2_yz_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 686);

    auto ta2_yz_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 687);

    auto ta2_yz_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 688);

    auto ta2_yz_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 689);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xzz_xxx_1, ta1_z_xzz_xxz_1, ta1_z_xzz_xzz_1, ta2_yz_xyzz_xxx_0, ta2_yz_xyzz_xxy_0, ta2_yz_xyzz_xxz_0, ta2_yz_xyzz_xyy_0, ta2_yz_xyzz_xyz_0, ta2_yz_xyzz_xzz_0, ta2_yz_xyzz_yyy_0, ta2_yz_xyzz_yyz_0, ta2_yz_xyzz_yzz_0, ta2_yz_xyzz_zzz_0, ta2_yz_xzz_xxx_0, ta2_yz_xzz_xxx_1, ta2_yz_xzz_xxz_0, ta2_yz_xzz_xxz_1, ta2_yz_xzz_xzz_0, ta2_yz_xzz_xzz_1, ta2_yz_yzz_xxy_0, ta2_yz_yzz_xxy_1, ta2_yz_yzz_xy_0, ta2_yz_yzz_xy_1, ta2_yz_yzz_xyy_0, ta2_yz_yzz_xyy_1, ta2_yz_yzz_xyz_0, ta2_yz_yzz_xyz_1, ta2_yz_yzz_yy_0, ta2_yz_yzz_yy_1, ta2_yz_yzz_yyy_0, ta2_yz_yzz_yyy_1, ta2_yz_yzz_yyz_0, ta2_yz_yzz_yyz_1, ta2_yz_yzz_yz_0, ta2_yz_yzz_yz_1, ta2_yz_yzz_yzz_0, ta2_yz_yzz_yzz_1, ta2_yz_yzz_zzz_0, ta2_yz_yzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyzz_xxx_0[i] = ta1_z_xzz_xxx_1[i] + ta2_yz_xzz_xxx_0[i] * pa_y[i] - ta2_yz_xzz_xxx_1[i] * pc_y[i];

        ta2_yz_xyzz_xxy_0[i] = 2.0 * ta2_yz_yzz_xy_0[i] * fe_0 - 2.0 * ta2_yz_yzz_xy_1[i] * fe_0 + ta2_yz_yzz_xxy_0[i] * pa_x[i] - ta2_yz_yzz_xxy_1[i] * pc_x[i];

        ta2_yz_xyzz_xxz_0[i] = ta1_z_xzz_xxz_1[i] + ta2_yz_xzz_xxz_0[i] * pa_y[i] - ta2_yz_xzz_xxz_1[i] * pc_y[i];

        ta2_yz_xyzz_xyy_0[i] = ta2_yz_yzz_yy_0[i] * fe_0 - ta2_yz_yzz_yy_1[i] * fe_0 + ta2_yz_yzz_xyy_0[i] * pa_x[i] - ta2_yz_yzz_xyy_1[i] * pc_x[i];

        ta2_yz_xyzz_xyz_0[i] = ta2_yz_yzz_yz_0[i] * fe_0 - ta2_yz_yzz_yz_1[i] * fe_0 + ta2_yz_yzz_xyz_0[i] * pa_x[i] - ta2_yz_yzz_xyz_1[i] * pc_x[i];

        ta2_yz_xyzz_xzz_0[i] = ta1_z_xzz_xzz_1[i] + ta2_yz_xzz_xzz_0[i] * pa_y[i] - ta2_yz_xzz_xzz_1[i] * pc_y[i];

        ta2_yz_xyzz_yyy_0[i] = ta2_yz_yzz_yyy_0[i] * pa_x[i] - ta2_yz_yzz_yyy_1[i] * pc_x[i];

        ta2_yz_xyzz_yyz_0[i] = ta2_yz_yzz_yyz_0[i] * pa_x[i] - ta2_yz_yzz_yyz_1[i] * pc_x[i];

        ta2_yz_xyzz_yzz_0[i] = ta2_yz_yzz_yzz_0[i] * pa_x[i] - ta2_yz_yzz_yzz_1[i] * pc_x[i];

        ta2_yz_xyzz_zzz_0[i] = ta2_yz_yzz_zzz_0[i] * pa_x[i] - ta2_yz_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 690-700 components of targeted buffer : GF

    auto ta2_yz_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 690);

    auto ta2_yz_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 691);

    auto ta2_yz_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 692);

    auto ta2_yz_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 693);

    auto ta2_yz_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 694);

    auto ta2_yz_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 695);

    auto ta2_yz_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 696);

    auto ta2_yz_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 697);

    auto ta2_yz_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 698);

    auto ta2_yz_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 699);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xzzz_xxx_0, ta2_yz_xzzz_xxy_0, ta2_yz_xzzz_xxz_0, ta2_yz_xzzz_xyy_0, ta2_yz_xzzz_xyz_0, ta2_yz_xzzz_xzz_0, ta2_yz_xzzz_yyy_0, ta2_yz_xzzz_yyz_0, ta2_yz_xzzz_yzz_0, ta2_yz_xzzz_zzz_0, ta2_yz_zzz_xx_0, ta2_yz_zzz_xx_1, ta2_yz_zzz_xxx_0, ta2_yz_zzz_xxx_1, ta2_yz_zzz_xxy_0, ta2_yz_zzz_xxy_1, ta2_yz_zzz_xxz_0, ta2_yz_zzz_xxz_1, ta2_yz_zzz_xy_0, ta2_yz_zzz_xy_1, ta2_yz_zzz_xyy_0, ta2_yz_zzz_xyy_1, ta2_yz_zzz_xyz_0, ta2_yz_zzz_xyz_1, ta2_yz_zzz_xz_0, ta2_yz_zzz_xz_1, ta2_yz_zzz_xzz_0, ta2_yz_zzz_xzz_1, ta2_yz_zzz_yy_0, ta2_yz_zzz_yy_1, ta2_yz_zzz_yyy_0, ta2_yz_zzz_yyy_1, ta2_yz_zzz_yyz_0, ta2_yz_zzz_yyz_1, ta2_yz_zzz_yz_0, ta2_yz_zzz_yz_1, ta2_yz_zzz_yzz_0, ta2_yz_zzz_yzz_1, ta2_yz_zzz_zz_0, ta2_yz_zzz_zz_1, ta2_yz_zzz_zzz_0, ta2_yz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzzz_xxx_0[i] = 3.0 * ta2_yz_zzz_xx_0[i] * fe_0 - 3.0 * ta2_yz_zzz_xx_1[i] * fe_0 + ta2_yz_zzz_xxx_0[i] * pa_x[i] - ta2_yz_zzz_xxx_1[i] * pc_x[i];

        ta2_yz_xzzz_xxy_0[i] = 2.0 * ta2_yz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xy_1[i] * fe_0 + ta2_yz_zzz_xxy_0[i] * pa_x[i] - ta2_yz_zzz_xxy_1[i] * pc_x[i];

        ta2_yz_xzzz_xxz_0[i] = 2.0 * ta2_yz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xz_1[i] * fe_0 + ta2_yz_zzz_xxz_0[i] * pa_x[i] - ta2_yz_zzz_xxz_1[i] * pc_x[i];

        ta2_yz_xzzz_xyy_0[i] = ta2_yz_zzz_yy_0[i] * fe_0 - ta2_yz_zzz_yy_1[i] * fe_0 + ta2_yz_zzz_xyy_0[i] * pa_x[i] - ta2_yz_zzz_xyy_1[i] * pc_x[i];

        ta2_yz_xzzz_xyz_0[i] = ta2_yz_zzz_yz_0[i] * fe_0 - ta2_yz_zzz_yz_1[i] * fe_0 + ta2_yz_zzz_xyz_0[i] * pa_x[i] - ta2_yz_zzz_xyz_1[i] * pc_x[i];

        ta2_yz_xzzz_xzz_0[i] = ta2_yz_zzz_zz_0[i] * fe_0 - ta2_yz_zzz_zz_1[i] * fe_0 + ta2_yz_zzz_xzz_0[i] * pa_x[i] - ta2_yz_zzz_xzz_1[i] * pc_x[i];

        ta2_yz_xzzz_yyy_0[i] = ta2_yz_zzz_yyy_0[i] * pa_x[i] - ta2_yz_zzz_yyy_1[i] * pc_x[i];

        ta2_yz_xzzz_yyz_0[i] = ta2_yz_zzz_yyz_0[i] * pa_x[i] - ta2_yz_zzz_yyz_1[i] * pc_x[i];

        ta2_yz_xzzz_yzz_0[i] = ta2_yz_zzz_yzz_0[i] * pa_x[i] - ta2_yz_zzz_yzz_1[i] * pc_x[i];

        ta2_yz_xzzz_zzz_0[i] = ta2_yz_zzz_zzz_0[i] * pa_x[i] - ta2_yz_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 700-710 components of targeted buffer : GF

    auto ta2_yz_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 700);

    auto ta2_yz_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 701);

    auto ta2_yz_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 702);

    auto ta2_yz_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 703);

    auto ta2_yz_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 704);

    auto ta2_yz_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 705);

    auto ta2_yz_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 706);

    auto ta2_yz_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 707);

    auto ta2_yz_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 708);

    auto ta2_yz_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 709);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yyy_xxx_1, ta1_z_yyy_xxy_1, ta1_z_yyy_xxz_1, ta1_z_yyy_xyy_1, ta1_z_yyy_xyz_1, ta1_z_yyy_xzz_1, ta1_z_yyy_yyy_1, ta1_z_yyy_yyz_1, ta1_z_yyy_yzz_1, ta1_z_yyy_zzz_1, ta2_yz_yy_xxx_0, ta2_yz_yy_xxx_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xxz_0, ta2_yz_yy_xxz_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_xyz_0, ta2_yz_yy_xyz_1, ta2_yz_yy_xzz_0, ta2_yz_yy_xzz_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yy_yyz_0, ta2_yz_yy_yyz_1, ta2_yz_yy_yzz_0, ta2_yz_yy_yzz_1, ta2_yz_yy_zzz_0, ta2_yz_yy_zzz_1, ta2_yz_yyy_xx_0, ta2_yz_yyy_xx_1, ta2_yz_yyy_xxx_0, ta2_yz_yyy_xxx_1, ta2_yz_yyy_xxy_0, ta2_yz_yyy_xxy_1, ta2_yz_yyy_xxz_0, ta2_yz_yyy_xxz_1, ta2_yz_yyy_xy_0, ta2_yz_yyy_xy_1, ta2_yz_yyy_xyy_0, ta2_yz_yyy_xyy_1, ta2_yz_yyy_xyz_0, ta2_yz_yyy_xyz_1, ta2_yz_yyy_xz_0, ta2_yz_yyy_xz_1, ta2_yz_yyy_xzz_0, ta2_yz_yyy_xzz_1, ta2_yz_yyy_yy_0, ta2_yz_yyy_yy_1, ta2_yz_yyy_yyy_0, ta2_yz_yyy_yyy_1, ta2_yz_yyy_yyz_0, ta2_yz_yyy_yyz_1, ta2_yz_yyy_yz_0, ta2_yz_yyy_yz_1, ta2_yz_yyy_yzz_0, ta2_yz_yyy_yzz_1, ta2_yz_yyy_zz_0, ta2_yz_yyy_zz_1, ta2_yz_yyy_zzz_0, ta2_yz_yyy_zzz_1, ta2_yz_yyyy_xxx_0, ta2_yz_yyyy_xxy_0, ta2_yz_yyyy_xxz_0, ta2_yz_yyyy_xyy_0, ta2_yz_yyyy_xyz_0, ta2_yz_yyyy_xzz_0, ta2_yz_yyyy_yyy_0, ta2_yz_yyyy_yyz_0, ta2_yz_yyyy_yzz_0, ta2_yz_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyy_xxx_0[i] = 3.0 * ta2_yz_yy_xxx_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxx_1[i] * fe_0 + ta1_z_yyy_xxx_1[i] + ta2_yz_yyy_xxx_0[i] * pa_y[i] - ta2_yz_yyy_xxx_1[i] * pc_y[i];

        ta2_yz_yyyy_xxy_0[i] = 3.0 * ta2_yz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxy_1[i] * fe_0 + ta2_yz_yyy_xx_0[i] * fe_0 - ta2_yz_yyy_xx_1[i] * fe_0 + ta1_z_yyy_xxy_1[i] + ta2_yz_yyy_xxy_0[i] * pa_y[i] - ta2_yz_yyy_xxy_1[i] * pc_y[i];

        ta2_yz_yyyy_xxz_0[i] = 3.0 * ta2_yz_yy_xxz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxz_1[i] * fe_0 + ta1_z_yyy_xxz_1[i] + ta2_yz_yyy_xxz_0[i] * pa_y[i] - ta2_yz_yyy_xxz_1[i] * pc_y[i];

        ta2_yz_yyyy_xyy_0[i] = 3.0 * ta2_yz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyy_1[i] * fe_0 + 2.0 * ta2_yz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xy_1[i] * fe_0 + ta1_z_yyy_xyy_1[i] + ta2_yz_yyy_xyy_0[i] * pa_y[i] - ta2_yz_yyy_xyy_1[i] * pc_y[i];

        ta2_yz_yyyy_xyz_0[i] = 3.0 * ta2_yz_yy_xyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyz_1[i] * fe_0 + ta2_yz_yyy_xz_0[i] * fe_0 - ta2_yz_yyy_xz_1[i] * fe_0 + ta1_z_yyy_xyz_1[i] + ta2_yz_yyy_xyz_0[i] * pa_y[i] - ta2_yz_yyy_xyz_1[i] * pc_y[i];

        ta2_yz_yyyy_xzz_0[i] = 3.0 * ta2_yz_yy_xzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xzz_1[i] * fe_0 + ta1_z_yyy_xzz_1[i] + ta2_yz_yyy_xzz_0[i] * pa_y[i] - ta2_yz_yyy_xzz_1[i] * pc_y[i];

        ta2_yz_yyyy_yyy_0[i] = 3.0 * ta2_yz_yy_yyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyy_1[i] * fe_0 + 3.0 * ta2_yz_yyy_yy_0[i] * fe_0 - 3.0 * ta2_yz_yyy_yy_1[i] * fe_0 + ta1_z_yyy_yyy_1[i] + ta2_yz_yyy_yyy_0[i] * pa_y[i] - ta2_yz_yyy_yyy_1[i] * pc_y[i];

        ta2_yz_yyyy_yyz_0[i] = 3.0 * ta2_yz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyz_1[i] * fe_0 + 2.0 * ta2_yz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_yz_1[i] * fe_0 + ta1_z_yyy_yyz_1[i] + ta2_yz_yyy_yyz_0[i] * pa_y[i] - ta2_yz_yyy_yyz_1[i] * pc_y[i];

        ta2_yz_yyyy_yzz_0[i] = 3.0 * ta2_yz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yzz_1[i] * fe_0 + ta2_yz_yyy_zz_0[i] * fe_0 - ta2_yz_yyy_zz_1[i] * fe_0 + ta1_z_yyy_yzz_1[i] + ta2_yz_yyy_yzz_0[i] * pa_y[i] - ta2_yz_yyy_yzz_1[i] * pc_y[i];

        ta2_yz_yyyy_zzz_0[i] = 3.0 * ta2_yz_yy_zzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_zzz_1[i] * fe_0 + ta1_z_yyy_zzz_1[i] + ta2_yz_yyy_zzz_0[i] * pa_y[i] - ta2_yz_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 710-720 components of targeted buffer : GF

    auto ta2_yz_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 710);

    auto ta2_yz_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 711);

    auto ta2_yz_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 712);

    auto ta2_yz_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 713);

    auto ta2_yz_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 714);

    auto ta2_yz_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 715);

    auto ta2_yz_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 716);

    auto ta2_yz_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 717);

    auto ta2_yz_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 718);

    auto ta2_yz_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 719);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyy_xxx_1, ta1_y_yyy_xxy_1, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_1, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_1, ta1_z_yyz_xxz_1, ta1_z_yyz_xzz_1, ta1_z_yyz_zzz_1, ta2_yz_yyy_xxx_0, ta2_yz_yyy_xxx_1, ta2_yz_yyy_xxy_0, ta2_yz_yyy_xxy_1, ta2_yz_yyy_xy_0, ta2_yz_yyy_xy_1, ta2_yz_yyy_xyy_0, ta2_yz_yyy_xyy_1, ta2_yz_yyy_xyz_0, ta2_yz_yyy_xyz_1, ta2_yz_yyy_yy_0, ta2_yz_yyy_yy_1, ta2_yz_yyy_yyy_0, ta2_yz_yyy_yyy_1, ta2_yz_yyy_yyz_0, ta2_yz_yyy_yyz_1, ta2_yz_yyy_yz_0, ta2_yz_yyy_yz_1, ta2_yz_yyy_yzz_0, ta2_yz_yyy_yzz_1, ta2_yz_yyyz_xxx_0, ta2_yz_yyyz_xxy_0, ta2_yz_yyyz_xxz_0, ta2_yz_yyyz_xyy_0, ta2_yz_yyyz_xyz_0, ta2_yz_yyyz_xzz_0, ta2_yz_yyyz_yyy_0, ta2_yz_yyyz_yyz_0, ta2_yz_yyyz_yzz_0, ta2_yz_yyyz_zzz_0, ta2_yz_yyz_xxz_0, ta2_yz_yyz_xxz_1, ta2_yz_yyz_xzz_0, ta2_yz_yyz_xzz_1, ta2_yz_yyz_zzz_0, ta2_yz_yyz_zzz_1, ta2_yz_yz_xxz_0, ta2_yz_yz_xxz_1, ta2_yz_yz_xzz_0, ta2_yz_yz_xzz_1, ta2_yz_yz_zzz_0, ta2_yz_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyz_xxx_0[i] = ta1_y_yyy_xxx_1[i] + ta2_yz_yyy_xxx_0[i] * pa_z[i] - ta2_yz_yyy_xxx_1[i] * pc_z[i];

        ta2_yz_yyyz_xxy_0[i] = ta1_y_yyy_xxy_1[i] + ta2_yz_yyy_xxy_0[i] * pa_z[i] - ta2_yz_yyy_xxy_1[i] * pc_z[i];

        ta2_yz_yyyz_xxz_0[i] = 2.0 * ta2_yz_yz_xxz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xxz_1[i] * fe_0 + ta1_z_yyz_xxz_1[i] + ta2_yz_yyz_xxz_0[i] * pa_y[i] - ta2_yz_yyz_xxz_1[i] * pc_y[i];

        ta2_yz_yyyz_xyy_0[i] = ta1_y_yyy_xyy_1[i] + ta2_yz_yyy_xyy_0[i] * pa_z[i] - ta2_yz_yyy_xyy_1[i] * pc_z[i];

        ta2_yz_yyyz_xyz_0[i] = ta2_yz_yyy_xy_0[i] * fe_0 - ta2_yz_yyy_xy_1[i] * fe_0 + ta1_y_yyy_xyz_1[i] + ta2_yz_yyy_xyz_0[i] * pa_z[i] - ta2_yz_yyy_xyz_1[i] * pc_z[i];

        ta2_yz_yyyz_xzz_0[i] = 2.0 * ta2_yz_yz_xzz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xzz_1[i] * fe_0 + ta1_z_yyz_xzz_1[i] + ta2_yz_yyz_xzz_0[i] * pa_y[i] - ta2_yz_yyz_xzz_1[i] * pc_y[i];

        ta2_yz_yyyz_yyy_0[i] = ta1_y_yyy_yyy_1[i] + ta2_yz_yyy_yyy_0[i] * pa_z[i] - ta2_yz_yyy_yyy_1[i] * pc_z[i];

        ta2_yz_yyyz_yyz_0[i] = ta2_yz_yyy_yy_0[i] * fe_0 - ta2_yz_yyy_yy_1[i] * fe_0 + ta1_y_yyy_yyz_1[i] + ta2_yz_yyy_yyz_0[i] * pa_z[i] - ta2_yz_yyy_yyz_1[i] * pc_z[i];

        ta2_yz_yyyz_yzz_0[i] = 2.0 * ta2_yz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_yz_1[i] * fe_0 + ta1_y_yyy_yzz_1[i] + ta2_yz_yyy_yzz_0[i] * pa_z[i] - ta2_yz_yyy_yzz_1[i] * pc_z[i];

        ta2_yz_yyyz_zzz_0[i] = 2.0 * ta2_yz_yz_zzz_0[i] * fe_0 - 2.0 * ta2_yz_yz_zzz_1[i] * fe_0 + ta1_z_yyz_zzz_1[i] + ta2_yz_yyz_zzz_0[i] * pa_y[i] - ta2_yz_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 720-730 components of targeted buffer : GF

    auto ta2_yz_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 720);

    auto ta2_yz_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 721);

    auto ta2_yz_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 722);

    auto ta2_yz_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 723);

    auto ta2_yz_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 724);

    auto ta2_yz_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 725);

    auto ta2_yz_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 726);

    auto ta2_yz_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 727);

    auto ta2_yz_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 728);

    auto ta2_yz_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 729);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyz_xxy_1, ta1_y_yyz_xyy_1, ta1_y_yyz_yyy_1, ta1_z_yzz_xxx_1, ta1_z_yzz_xxz_1, ta1_z_yzz_xyz_1, ta1_z_yzz_xzz_1, ta1_z_yzz_yyz_1, ta1_z_yzz_yzz_1, ta1_z_yzz_zzz_1, ta2_yz_yy_xxy_0, ta2_yz_yy_xxy_1, ta2_yz_yy_xyy_0, ta2_yz_yy_xyy_1, ta2_yz_yy_yyy_0, ta2_yz_yy_yyy_1, ta2_yz_yyz_xxy_0, ta2_yz_yyz_xxy_1, ta2_yz_yyz_xyy_0, ta2_yz_yyz_xyy_1, ta2_yz_yyz_yyy_0, ta2_yz_yyz_yyy_1, ta2_yz_yyzz_xxx_0, ta2_yz_yyzz_xxy_0, ta2_yz_yyzz_xxz_0, ta2_yz_yyzz_xyy_0, ta2_yz_yyzz_xyz_0, ta2_yz_yyzz_xzz_0, ta2_yz_yyzz_yyy_0, ta2_yz_yyzz_yyz_0, ta2_yz_yyzz_yzz_0, ta2_yz_yyzz_zzz_0, ta2_yz_yzz_xxx_0, ta2_yz_yzz_xxx_1, ta2_yz_yzz_xxz_0, ta2_yz_yzz_xxz_1, ta2_yz_yzz_xyz_0, ta2_yz_yzz_xyz_1, ta2_yz_yzz_xz_0, ta2_yz_yzz_xz_1, ta2_yz_yzz_xzz_0, ta2_yz_yzz_xzz_1, ta2_yz_yzz_yyz_0, ta2_yz_yzz_yyz_1, ta2_yz_yzz_yz_0, ta2_yz_yzz_yz_1, ta2_yz_yzz_yzz_0, ta2_yz_yzz_yzz_1, ta2_yz_yzz_zz_0, ta2_yz_yzz_zz_1, ta2_yz_yzz_zzz_0, ta2_yz_yzz_zzz_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyzz_xxx_0[i] = ta2_yz_zz_xxx_0[i] * fe_0 - ta2_yz_zz_xxx_1[i] * fe_0 + ta1_z_yzz_xxx_1[i] + ta2_yz_yzz_xxx_0[i] * pa_y[i] - ta2_yz_yzz_xxx_1[i] * pc_y[i];

        ta2_yz_yyzz_xxy_0[i] = ta2_yz_yy_xxy_0[i] * fe_0 - ta2_yz_yy_xxy_1[i] * fe_0 + ta1_y_yyz_xxy_1[i] + ta2_yz_yyz_xxy_0[i] * pa_z[i] - ta2_yz_yyz_xxy_1[i] * pc_z[i];

        ta2_yz_yyzz_xxz_0[i] = ta2_yz_zz_xxz_0[i] * fe_0 - ta2_yz_zz_xxz_1[i] * fe_0 + ta1_z_yzz_xxz_1[i] + ta2_yz_yzz_xxz_0[i] * pa_y[i] - ta2_yz_yzz_xxz_1[i] * pc_y[i];

        ta2_yz_yyzz_xyy_0[i] = ta2_yz_yy_xyy_0[i] * fe_0 - ta2_yz_yy_xyy_1[i] * fe_0 + ta1_y_yyz_xyy_1[i] + ta2_yz_yyz_xyy_0[i] * pa_z[i] - ta2_yz_yyz_xyy_1[i] * pc_z[i];

        ta2_yz_yyzz_xyz_0[i] = ta2_yz_zz_xyz_0[i] * fe_0 - ta2_yz_zz_xyz_1[i] * fe_0 + ta2_yz_yzz_xz_0[i] * fe_0 - ta2_yz_yzz_xz_1[i] * fe_0 + ta1_z_yzz_xyz_1[i] + ta2_yz_yzz_xyz_0[i] * pa_y[i] - ta2_yz_yzz_xyz_1[i] * pc_y[i];

        ta2_yz_yyzz_xzz_0[i] = ta2_yz_zz_xzz_0[i] * fe_0 - ta2_yz_zz_xzz_1[i] * fe_0 + ta1_z_yzz_xzz_1[i] + ta2_yz_yzz_xzz_0[i] * pa_y[i] - ta2_yz_yzz_xzz_1[i] * pc_y[i];

        ta2_yz_yyzz_yyy_0[i] = ta2_yz_yy_yyy_0[i] * fe_0 - ta2_yz_yy_yyy_1[i] * fe_0 + ta1_y_yyz_yyy_1[i] + ta2_yz_yyz_yyy_0[i] * pa_z[i] - ta2_yz_yyz_yyy_1[i] * pc_z[i];

        ta2_yz_yyzz_yyz_0[i] = ta2_yz_zz_yyz_0[i] * fe_0 - ta2_yz_zz_yyz_1[i] * fe_0 + 2.0 * ta2_yz_yzz_yz_0[i] * fe_0 - 2.0 * ta2_yz_yzz_yz_1[i] * fe_0 + ta1_z_yzz_yyz_1[i] + ta2_yz_yzz_yyz_0[i] * pa_y[i] - ta2_yz_yzz_yyz_1[i] * pc_y[i];

        ta2_yz_yyzz_yzz_0[i] = ta2_yz_zz_yzz_0[i] * fe_0 - ta2_yz_zz_yzz_1[i] * fe_0 + ta2_yz_yzz_zz_0[i] * fe_0 - ta2_yz_yzz_zz_1[i] * fe_0 + ta1_z_yzz_yzz_1[i] + ta2_yz_yzz_yzz_0[i] * pa_y[i] - ta2_yz_yzz_yzz_1[i] * pc_y[i];

        ta2_yz_yyzz_zzz_0[i] = ta2_yz_zz_zzz_0[i] * fe_0 - ta2_yz_zz_zzz_1[i] * fe_0 + ta1_z_yzz_zzz_1[i] + ta2_yz_yzz_zzz_0[i] * pa_y[i] - ta2_yz_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 730-740 components of targeted buffer : GF

    auto ta2_yz_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 730);

    auto ta2_yz_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 731);

    auto ta2_yz_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 732);

    auto ta2_yz_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 733);

    auto ta2_yz_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 734);

    auto ta2_yz_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 735);

    auto ta2_yz_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 736);

    auto ta2_yz_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 737);

    auto ta2_yz_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 738);

    auto ta2_yz_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 739);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_zzz_xxx_1, ta1_z_zzz_xxy_1, ta1_z_zzz_xxz_1, ta1_z_zzz_xyy_1, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_1, ta1_z_zzz_yyy_1, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_1, ta2_yz_yzzz_xxx_0, ta2_yz_yzzz_xxy_0, ta2_yz_yzzz_xxz_0, ta2_yz_yzzz_xyy_0, ta2_yz_yzzz_xyz_0, ta2_yz_yzzz_xzz_0, ta2_yz_yzzz_yyy_0, ta2_yz_yzzz_yyz_0, ta2_yz_yzzz_yzz_0, ta2_yz_yzzz_zzz_0, ta2_yz_zzz_xx_0, ta2_yz_zzz_xx_1, ta2_yz_zzz_xxx_0, ta2_yz_zzz_xxx_1, ta2_yz_zzz_xxy_0, ta2_yz_zzz_xxy_1, ta2_yz_zzz_xxz_0, ta2_yz_zzz_xxz_1, ta2_yz_zzz_xy_0, ta2_yz_zzz_xy_1, ta2_yz_zzz_xyy_0, ta2_yz_zzz_xyy_1, ta2_yz_zzz_xyz_0, ta2_yz_zzz_xyz_1, ta2_yz_zzz_xz_0, ta2_yz_zzz_xz_1, ta2_yz_zzz_xzz_0, ta2_yz_zzz_xzz_1, ta2_yz_zzz_yy_0, ta2_yz_zzz_yy_1, ta2_yz_zzz_yyy_0, ta2_yz_zzz_yyy_1, ta2_yz_zzz_yyz_0, ta2_yz_zzz_yyz_1, ta2_yz_zzz_yz_0, ta2_yz_zzz_yz_1, ta2_yz_zzz_yzz_0, ta2_yz_zzz_yzz_1, ta2_yz_zzz_zz_0, ta2_yz_zzz_zz_1, ta2_yz_zzz_zzz_0, ta2_yz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzzz_xxx_0[i] = ta1_z_zzz_xxx_1[i] + ta2_yz_zzz_xxx_0[i] * pa_y[i] - ta2_yz_zzz_xxx_1[i] * pc_y[i];

        ta2_yz_yzzz_xxy_0[i] = ta2_yz_zzz_xx_0[i] * fe_0 - ta2_yz_zzz_xx_1[i] * fe_0 + ta1_z_zzz_xxy_1[i] + ta2_yz_zzz_xxy_0[i] * pa_y[i] - ta2_yz_zzz_xxy_1[i] * pc_y[i];

        ta2_yz_yzzz_xxz_0[i] = ta1_z_zzz_xxz_1[i] + ta2_yz_zzz_xxz_0[i] * pa_y[i] - ta2_yz_zzz_xxz_1[i] * pc_y[i];

        ta2_yz_yzzz_xyy_0[i] = 2.0 * ta2_yz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xy_1[i] * fe_0 + ta1_z_zzz_xyy_1[i] + ta2_yz_zzz_xyy_0[i] * pa_y[i] - ta2_yz_zzz_xyy_1[i] * pc_y[i];

        ta2_yz_yzzz_xyz_0[i] = ta2_yz_zzz_xz_0[i] * fe_0 - ta2_yz_zzz_xz_1[i] * fe_0 + ta1_z_zzz_xyz_1[i] + ta2_yz_zzz_xyz_0[i] * pa_y[i] - ta2_yz_zzz_xyz_1[i] * pc_y[i];

        ta2_yz_yzzz_xzz_0[i] = ta1_z_zzz_xzz_1[i] + ta2_yz_zzz_xzz_0[i] * pa_y[i] - ta2_yz_zzz_xzz_1[i] * pc_y[i];

        ta2_yz_yzzz_yyy_0[i] = 3.0 * ta2_yz_zzz_yy_0[i] * fe_0 - 3.0 * ta2_yz_zzz_yy_1[i] * fe_0 + ta1_z_zzz_yyy_1[i] + ta2_yz_zzz_yyy_0[i] * pa_y[i] - ta2_yz_zzz_yyy_1[i] * pc_y[i];

        ta2_yz_yzzz_yyz_0[i] = 2.0 * ta2_yz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_yz_1[i] * fe_0 + ta1_z_zzz_yyz_1[i] + ta2_yz_zzz_yyz_0[i] * pa_y[i] - ta2_yz_zzz_yyz_1[i] * pc_y[i];

        ta2_yz_yzzz_yzz_0[i] = ta2_yz_zzz_zz_0[i] * fe_0 - ta2_yz_zzz_zz_1[i] * fe_0 + ta1_z_zzz_yzz_1[i] + ta2_yz_zzz_yzz_0[i] * pa_y[i] - ta2_yz_zzz_yzz_1[i] * pc_y[i];

        ta2_yz_yzzz_zzz_0[i] = ta1_z_zzz_zzz_1[i] + ta2_yz_zzz_zzz_0[i] * pa_y[i] - ta2_yz_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 740-750 components of targeted buffer : GF

    auto ta2_yz_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 740);

    auto ta2_yz_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 741);

    auto ta2_yz_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 742);

    auto ta2_yz_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 743);

    auto ta2_yz_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 744);

    auto ta2_yz_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 745);

    auto ta2_yz_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 746);

    auto ta2_yz_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 747);

    auto ta2_yz_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 748);

    auto ta2_yz_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 749);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zzz_xxx_1, ta1_y_zzz_xxy_1, ta1_y_zzz_xxz_1, ta1_y_zzz_xyy_1, ta1_y_zzz_xyz_1, ta1_y_zzz_xzz_1, ta1_y_zzz_yyy_1, ta1_y_zzz_yyz_1, ta1_y_zzz_yzz_1, ta1_y_zzz_zzz_1, ta2_yz_zz_xxx_0, ta2_yz_zz_xxx_1, ta2_yz_zz_xxy_0, ta2_yz_zz_xxy_1, ta2_yz_zz_xxz_0, ta2_yz_zz_xxz_1, ta2_yz_zz_xyy_0, ta2_yz_zz_xyy_1, ta2_yz_zz_xyz_0, ta2_yz_zz_xyz_1, ta2_yz_zz_xzz_0, ta2_yz_zz_xzz_1, ta2_yz_zz_yyy_0, ta2_yz_zz_yyy_1, ta2_yz_zz_yyz_0, ta2_yz_zz_yyz_1, ta2_yz_zz_yzz_0, ta2_yz_zz_yzz_1, ta2_yz_zz_zzz_0, ta2_yz_zz_zzz_1, ta2_yz_zzz_xx_0, ta2_yz_zzz_xx_1, ta2_yz_zzz_xxx_0, ta2_yz_zzz_xxx_1, ta2_yz_zzz_xxy_0, ta2_yz_zzz_xxy_1, ta2_yz_zzz_xxz_0, ta2_yz_zzz_xxz_1, ta2_yz_zzz_xy_0, ta2_yz_zzz_xy_1, ta2_yz_zzz_xyy_0, ta2_yz_zzz_xyy_1, ta2_yz_zzz_xyz_0, ta2_yz_zzz_xyz_1, ta2_yz_zzz_xz_0, ta2_yz_zzz_xz_1, ta2_yz_zzz_xzz_0, ta2_yz_zzz_xzz_1, ta2_yz_zzz_yy_0, ta2_yz_zzz_yy_1, ta2_yz_zzz_yyy_0, ta2_yz_zzz_yyy_1, ta2_yz_zzz_yyz_0, ta2_yz_zzz_yyz_1, ta2_yz_zzz_yz_0, ta2_yz_zzz_yz_1, ta2_yz_zzz_yzz_0, ta2_yz_zzz_yzz_1, ta2_yz_zzz_zz_0, ta2_yz_zzz_zz_1, ta2_yz_zzz_zzz_0, ta2_yz_zzz_zzz_1, ta2_yz_zzzz_xxx_0, ta2_yz_zzzz_xxy_0, ta2_yz_zzzz_xxz_0, ta2_yz_zzzz_xyy_0, ta2_yz_zzzz_xyz_0, ta2_yz_zzzz_xzz_0, ta2_yz_zzzz_yyy_0, ta2_yz_zzzz_yyz_0, ta2_yz_zzzz_yzz_0, ta2_yz_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzzz_xxx_0[i] = 3.0 * ta2_yz_zz_xxx_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxx_1[i] * fe_0 + ta1_y_zzz_xxx_1[i] + ta2_yz_zzz_xxx_0[i] * pa_z[i] - ta2_yz_zzz_xxx_1[i] * pc_z[i];

        ta2_yz_zzzz_xxy_0[i] = 3.0 * ta2_yz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxy_1[i] * fe_0 + ta1_y_zzz_xxy_1[i] + ta2_yz_zzz_xxy_0[i] * pa_z[i] - ta2_yz_zzz_xxy_1[i] * pc_z[i];

        ta2_yz_zzzz_xxz_0[i] = 3.0 * ta2_yz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxz_1[i] * fe_0 + ta2_yz_zzz_xx_0[i] * fe_0 - ta2_yz_zzz_xx_1[i] * fe_0 + ta1_y_zzz_xxz_1[i] + ta2_yz_zzz_xxz_0[i] * pa_z[i] - ta2_yz_zzz_xxz_1[i] * pc_z[i];

        ta2_yz_zzzz_xyy_0[i] = 3.0 * ta2_yz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyy_1[i] * fe_0 + ta1_y_zzz_xyy_1[i] + ta2_yz_zzz_xyy_0[i] * pa_z[i] - ta2_yz_zzz_xyy_1[i] * pc_z[i];

        ta2_yz_zzzz_xyz_0[i] = 3.0 * ta2_yz_zz_xyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyz_1[i] * fe_0 + ta2_yz_zzz_xy_0[i] * fe_0 - ta2_yz_zzz_xy_1[i] * fe_0 + ta1_y_zzz_xyz_1[i] + ta2_yz_zzz_xyz_0[i] * pa_z[i] - ta2_yz_zzz_xyz_1[i] * pc_z[i];

        ta2_yz_zzzz_xzz_0[i] = 3.0 * ta2_yz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xzz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xz_1[i] * fe_0 + ta1_y_zzz_xzz_1[i] + ta2_yz_zzz_xzz_0[i] * pa_z[i] - ta2_yz_zzz_xzz_1[i] * pc_z[i];

        ta2_yz_zzzz_yyy_0[i] = 3.0 * ta2_yz_zz_yyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyy_1[i] * fe_0 + ta1_y_zzz_yyy_1[i] + ta2_yz_zzz_yyy_0[i] * pa_z[i] - ta2_yz_zzz_yyy_1[i] * pc_z[i];

        ta2_yz_zzzz_yyz_0[i] = 3.0 * ta2_yz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyz_1[i] * fe_0 + ta2_yz_zzz_yy_0[i] * fe_0 - ta2_yz_zzz_yy_1[i] * fe_0 + ta1_y_zzz_yyz_1[i] + ta2_yz_zzz_yyz_0[i] * pa_z[i] - ta2_yz_zzz_yyz_1[i] * pc_z[i];

        ta2_yz_zzzz_yzz_0[i] = 3.0 * ta2_yz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yzz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_yz_1[i] * fe_0 + ta1_y_zzz_yzz_1[i] + ta2_yz_zzz_yzz_0[i] * pa_z[i] - ta2_yz_zzz_yzz_1[i] * pc_z[i];

        ta2_yz_zzzz_zzz_0[i] = 3.0 * ta2_yz_zz_zzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_zzz_1[i] * fe_0 + 3.0 * ta2_yz_zzz_zz_0[i] * fe_0 - 3.0 * ta2_yz_zzz_zz_1[i] * fe_0 + ta1_y_zzz_zzz_1[i] + ta2_yz_zzz_zzz_0[i] * pa_z[i] - ta2_yz_zzz_zzz_1[i] * pc_z[i];
    }

    // Set up 750-760 components of targeted buffer : GF

    auto ta2_zz_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 750);

    auto ta2_zz_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 751);

    auto ta2_zz_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 752);

    auto ta2_zz_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 753);

    auto ta2_zz_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 754);

    auto ta2_zz_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 755);

    auto ta2_zz_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 756);

    auto ta2_zz_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 757);

    auto ta2_zz_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 758);

    auto ta2_zz_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 759);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xx_xyz_0, ta2_zz_xx_xyz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xx_yyy_0, ta2_zz_xx_yyy_1, ta2_zz_xx_yyz_0, ta2_zz_xx_yyz_1, ta2_zz_xx_yzz_0, ta2_zz_xx_yzz_1, ta2_zz_xx_zzz_0, ta2_zz_xx_zzz_1, ta2_zz_xxx_xx_0, ta2_zz_xxx_xx_1, ta2_zz_xxx_xxx_0, ta2_zz_xxx_xxx_1, ta2_zz_xxx_xxy_0, ta2_zz_xxx_xxy_1, ta2_zz_xxx_xxz_0, ta2_zz_xxx_xxz_1, ta2_zz_xxx_xy_0, ta2_zz_xxx_xy_1, ta2_zz_xxx_xyy_0, ta2_zz_xxx_xyy_1, ta2_zz_xxx_xyz_0, ta2_zz_xxx_xyz_1, ta2_zz_xxx_xz_0, ta2_zz_xxx_xz_1, ta2_zz_xxx_xzz_0, ta2_zz_xxx_xzz_1, ta2_zz_xxx_yy_0, ta2_zz_xxx_yy_1, ta2_zz_xxx_yyy_0, ta2_zz_xxx_yyy_1, ta2_zz_xxx_yyz_0, ta2_zz_xxx_yyz_1, ta2_zz_xxx_yz_0, ta2_zz_xxx_yz_1, ta2_zz_xxx_yzz_0, ta2_zz_xxx_yzz_1, ta2_zz_xxx_zz_0, ta2_zz_xxx_zz_1, ta2_zz_xxx_zzz_0, ta2_zz_xxx_zzz_1, ta2_zz_xxxx_xxx_0, ta2_zz_xxxx_xxy_0, ta2_zz_xxxx_xxz_0, ta2_zz_xxxx_xyy_0, ta2_zz_xxxx_xyz_0, ta2_zz_xxxx_xzz_0, ta2_zz_xxxx_yyy_0, ta2_zz_xxxx_yyz_0, ta2_zz_xxxx_yzz_0, ta2_zz_xxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxx_xxx_0[i] = 3.0 * ta2_zz_xx_xxx_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxx_1[i] * fe_0 + 3.0 * ta2_zz_xxx_xx_0[i] * fe_0 - 3.0 * ta2_zz_xxx_xx_1[i] * fe_0 + ta2_zz_xxx_xxx_0[i] * pa_x[i] - ta2_zz_xxx_xxx_1[i] * pc_x[i];

        ta2_zz_xxxx_xxy_0[i] = 3.0 * ta2_zz_xx_xxy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxy_1[i] * fe_0 + 2.0 * ta2_zz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xy_1[i] * fe_0 + ta2_zz_xxx_xxy_0[i] * pa_x[i] - ta2_zz_xxx_xxy_1[i] * pc_x[i];

        ta2_zz_xxxx_xxz_0[i] = 3.0 * ta2_zz_xx_xxz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxz_1[i] * fe_0 + 2.0 * ta2_zz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xz_1[i] * fe_0 + ta2_zz_xxx_xxz_0[i] * pa_x[i] - ta2_zz_xxx_xxz_1[i] * pc_x[i];

        ta2_zz_xxxx_xyy_0[i] = 3.0 * ta2_zz_xx_xyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyy_1[i] * fe_0 + ta2_zz_xxx_yy_0[i] * fe_0 - ta2_zz_xxx_yy_1[i] * fe_0 + ta2_zz_xxx_xyy_0[i] * pa_x[i] - ta2_zz_xxx_xyy_1[i] * pc_x[i];

        ta2_zz_xxxx_xyz_0[i] = 3.0 * ta2_zz_xx_xyz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyz_1[i] * fe_0 + ta2_zz_xxx_yz_0[i] * fe_0 - ta2_zz_xxx_yz_1[i] * fe_0 + ta2_zz_xxx_xyz_0[i] * pa_x[i] - ta2_zz_xxx_xyz_1[i] * pc_x[i];

        ta2_zz_xxxx_xzz_0[i] = 3.0 * ta2_zz_xx_xzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xzz_1[i] * fe_0 + ta2_zz_xxx_zz_0[i] * fe_0 - ta2_zz_xxx_zz_1[i] * fe_0 + ta2_zz_xxx_xzz_0[i] * pa_x[i] - ta2_zz_xxx_xzz_1[i] * pc_x[i];

        ta2_zz_xxxx_yyy_0[i] = 3.0 * ta2_zz_xx_yyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_yyy_1[i] * fe_0 + ta2_zz_xxx_yyy_0[i] * pa_x[i] - ta2_zz_xxx_yyy_1[i] * pc_x[i];

        ta2_zz_xxxx_yyz_0[i] = 3.0 * ta2_zz_xx_yyz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yyz_1[i] * fe_0 + ta2_zz_xxx_yyz_0[i] * pa_x[i] - ta2_zz_xxx_yyz_1[i] * pc_x[i];

        ta2_zz_xxxx_yzz_0[i] = 3.0 * ta2_zz_xx_yzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yzz_1[i] * fe_0 + ta2_zz_xxx_yzz_0[i] * pa_x[i] - ta2_zz_xxx_yzz_1[i] * pc_x[i];

        ta2_zz_xxxx_zzz_0[i] = 3.0 * ta2_zz_xx_zzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_zzz_1[i] * fe_0 + ta2_zz_xxx_zzz_0[i] * pa_x[i] - ta2_zz_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 760-770 components of targeted buffer : GF

    auto ta2_zz_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 760);

    auto ta2_zz_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 761);

    auto ta2_zz_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 762);

    auto ta2_zz_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 763);

    auto ta2_zz_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 764);

    auto ta2_zz_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 765);

    auto ta2_zz_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 766);

    auto ta2_zz_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 767);

    auto ta2_zz_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 768);

    auto ta2_zz_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 769);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xxx_xx_0, ta2_zz_xxx_xx_1, ta2_zz_xxx_xxx_0, ta2_zz_xxx_xxx_1, ta2_zz_xxx_xxy_0, ta2_zz_xxx_xxy_1, ta2_zz_xxx_xxz_0, ta2_zz_xxx_xxz_1, ta2_zz_xxx_xy_0, ta2_zz_xxx_xy_1, ta2_zz_xxx_xyy_0, ta2_zz_xxx_xyy_1, ta2_zz_xxx_xyz_0, ta2_zz_xxx_xyz_1, ta2_zz_xxx_xz_0, ta2_zz_xxx_xz_1, ta2_zz_xxx_xzz_0, ta2_zz_xxx_xzz_1, ta2_zz_xxx_zzz_0, ta2_zz_xxx_zzz_1, ta2_zz_xxxy_xxx_0, ta2_zz_xxxy_xxy_0, ta2_zz_xxxy_xxz_0, ta2_zz_xxxy_xyy_0, ta2_zz_xxxy_xyz_0, ta2_zz_xxxy_xzz_0, ta2_zz_xxxy_yyy_0, ta2_zz_xxxy_yyz_0, ta2_zz_xxxy_yzz_0, ta2_zz_xxxy_zzz_0, ta2_zz_xxy_yyy_0, ta2_zz_xxy_yyy_1, ta2_zz_xxy_yyz_0, ta2_zz_xxy_yyz_1, ta2_zz_xxy_yzz_0, ta2_zz_xxy_yzz_1, ta2_zz_xy_yyy_0, ta2_zz_xy_yyy_1, ta2_zz_xy_yyz_0, ta2_zz_xy_yyz_1, ta2_zz_xy_yzz_0, ta2_zz_xy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxy_xxx_0[i] = ta2_zz_xxx_xxx_0[i] * pa_y[i] - ta2_zz_xxx_xxx_1[i] * pc_y[i];

        ta2_zz_xxxy_xxy_0[i] = ta2_zz_xxx_xx_0[i] * fe_0 - ta2_zz_xxx_xx_1[i] * fe_0 + ta2_zz_xxx_xxy_0[i] * pa_y[i] - ta2_zz_xxx_xxy_1[i] * pc_y[i];

        ta2_zz_xxxy_xxz_0[i] = ta2_zz_xxx_xxz_0[i] * pa_y[i] - ta2_zz_xxx_xxz_1[i] * pc_y[i];

        ta2_zz_xxxy_xyy_0[i] = 2.0 * ta2_zz_xxx_xy_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xy_1[i] * fe_0 + ta2_zz_xxx_xyy_0[i] * pa_y[i] - ta2_zz_xxx_xyy_1[i] * pc_y[i];

        ta2_zz_xxxy_xyz_0[i] = ta2_zz_xxx_xz_0[i] * fe_0 - ta2_zz_xxx_xz_1[i] * fe_0 + ta2_zz_xxx_xyz_0[i] * pa_y[i] - ta2_zz_xxx_xyz_1[i] * pc_y[i];

        ta2_zz_xxxy_xzz_0[i] = ta2_zz_xxx_xzz_0[i] * pa_y[i] - ta2_zz_xxx_xzz_1[i] * pc_y[i];

        ta2_zz_xxxy_yyy_0[i] = 2.0 * ta2_zz_xy_yyy_0[i] * fe_0 - 2.0 * ta2_zz_xy_yyy_1[i] * fe_0 + ta2_zz_xxy_yyy_0[i] * pa_x[i] - ta2_zz_xxy_yyy_1[i] * pc_x[i];

        ta2_zz_xxxy_yyz_0[i] = 2.0 * ta2_zz_xy_yyz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yyz_1[i] * fe_0 + ta2_zz_xxy_yyz_0[i] * pa_x[i] - ta2_zz_xxy_yyz_1[i] * pc_x[i];

        ta2_zz_xxxy_yzz_0[i] = 2.0 * ta2_zz_xy_yzz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yzz_1[i] * fe_0 + ta2_zz_xxy_yzz_0[i] * pa_x[i] - ta2_zz_xxy_yzz_1[i] * pc_x[i];

        ta2_zz_xxxy_zzz_0[i] = ta2_zz_xxx_zzz_0[i] * pa_y[i] - ta2_zz_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 770-780 components of targeted buffer : GF

    auto ta2_zz_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 770);

    auto ta2_zz_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 771);

    auto ta2_zz_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 772);

    auto ta2_zz_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 773);

    auto ta2_zz_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 774);

    auto ta2_zz_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 775);

    auto ta2_zz_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 776);

    auto ta2_zz_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 777);

    auto ta2_zz_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 778);

    auto ta2_zz_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 779);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxx_xxx_1, ta1_z_xxx_xxy_1, ta1_z_xxx_xxz_1, ta1_z_xxx_xyy_1, ta1_z_xxx_xyz_1, ta1_z_xxx_xzz_1, ta1_z_xxx_yyy_1, ta2_zz_xxx_xx_0, ta2_zz_xxx_xx_1, ta2_zz_xxx_xxx_0, ta2_zz_xxx_xxx_1, ta2_zz_xxx_xxy_0, ta2_zz_xxx_xxy_1, ta2_zz_xxx_xxz_0, ta2_zz_xxx_xxz_1, ta2_zz_xxx_xy_0, ta2_zz_xxx_xy_1, ta2_zz_xxx_xyy_0, ta2_zz_xxx_xyy_1, ta2_zz_xxx_xyz_0, ta2_zz_xxx_xyz_1, ta2_zz_xxx_xz_0, ta2_zz_xxx_xz_1, ta2_zz_xxx_xzz_0, ta2_zz_xxx_xzz_1, ta2_zz_xxx_yyy_0, ta2_zz_xxx_yyy_1, ta2_zz_xxxz_xxx_0, ta2_zz_xxxz_xxy_0, ta2_zz_xxxz_xxz_0, ta2_zz_xxxz_xyy_0, ta2_zz_xxxz_xyz_0, ta2_zz_xxxz_xzz_0, ta2_zz_xxxz_yyy_0, ta2_zz_xxxz_yyz_0, ta2_zz_xxxz_yzz_0, ta2_zz_xxxz_zzz_0, ta2_zz_xxz_yyz_0, ta2_zz_xxz_yyz_1, ta2_zz_xxz_yzz_0, ta2_zz_xxz_yzz_1, ta2_zz_xxz_zzz_0, ta2_zz_xxz_zzz_1, ta2_zz_xz_yyz_0, ta2_zz_xz_yyz_1, ta2_zz_xz_yzz_0, ta2_zz_xz_yzz_1, ta2_zz_xz_zzz_0, ta2_zz_xz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxz_xxx_0[i] = 2.0 * ta1_z_xxx_xxx_1[i] + ta2_zz_xxx_xxx_0[i] * pa_z[i] - ta2_zz_xxx_xxx_1[i] * pc_z[i];

        ta2_zz_xxxz_xxy_0[i] = 2.0 * ta1_z_xxx_xxy_1[i] + ta2_zz_xxx_xxy_0[i] * pa_z[i] - ta2_zz_xxx_xxy_1[i] * pc_z[i];

        ta2_zz_xxxz_xxz_0[i] = ta2_zz_xxx_xx_0[i] * fe_0 - ta2_zz_xxx_xx_1[i] * fe_0 + 2.0 * ta1_z_xxx_xxz_1[i] + ta2_zz_xxx_xxz_0[i] * pa_z[i] - ta2_zz_xxx_xxz_1[i] * pc_z[i];

        ta2_zz_xxxz_xyy_0[i] = 2.0 * ta1_z_xxx_xyy_1[i] + ta2_zz_xxx_xyy_0[i] * pa_z[i] - ta2_zz_xxx_xyy_1[i] * pc_z[i];

        ta2_zz_xxxz_xyz_0[i] = ta2_zz_xxx_xy_0[i] * fe_0 - ta2_zz_xxx_xy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyz_1[i] + ta2_zz_xxx_xyz_0[i] * pa_z[i] - ta2_zz_xxx_xyz_1[i] * pc_z[i];

        ta2_zz_xxxz_xzz_0[i] = 2.0 * ta2_zz_xxx_xz_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xzz_1[i] + ta2_zz_xxx_xzz_0[i] * pa_z[i] - ta2_zz_xxx_xzz_1[i] * pc_z[i];

        ta2_zz_xxxz_yyy_0[i] = 2.0 * ta1_z_xxx_yyy_1[i] + ta2_zz_xxx_yyy_0[i] * pa_z[i] - ta2_zz_xxx_yyy_1[i] * pc_z[i];

        ta2_zz_xxxz_yyz_0[i] = 2.0 * ta2_zz_xz_yyz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yyz_1[i] * fe_0 + ta2_zz_xxz_yyz_0[i] * pa_x[i] - ta2_zz_xxz_yyz_1[i] * pc_x[i];

        ta2_zz_xxxz_yzz_0[i] = 2.0 * ta2_zz_xz_yzz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yzz_1[i] * fe_0 + ta2_zz_xxz_yzz_0[i] * pa_x[i] - ta2_zz_xxz_yzz_1[i] * pc_x[i];

        ta2_zz_xxxz_zzz_0[i] = 2.0 * ta2_zz_xz_zzz_0[i] * fe_0 - 2.0 * ta2_zz_xz_zzz_1[i] * fe_0 + ta2_zz_xxz_zzz_0[i] * pa_x[i] - ta2_zz_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 780-790 components of targeted buffer : GF

    auto ta2_zz_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 780);

    auto ta2_zz_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 781);

    auto ta2_zz_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 782);

    auto ta2_zz_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 783);

    auto ta2_zz_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 784);

    auto ta2_zz_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 785);

    auto ta2_zz_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 786);

    auto ta2_zz_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 787);

    auto ta2_zz_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 788);

    auto ta2_zz_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 789);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxz_0, ta2_zz_xx_xxz_1, ta2_zz_xx_xzz_0, ta2_zz_xx_xzz_1, ta2_zz_xxy_xxx_0, ta2_zz_xxy_xxx_1, ta2_zz_xxy_xxz_0, ta2_zz_xxy_xxz_1, ta2_zz_xxy_xzz_0, ta2_zz_xxy_xzz_1, ta2_zz_xxyy_xxx_0, ta2_zz_xxyy_xxy_0, ta2_zz_xxyy_xxz_0, ta2_zz_xxyy_xyy_0, ta2_zz_xxyy_xyz_0, ta2_zz_xxyy_xzz_0, ta2_zz_xxyy_yyy_0, ta2_zz_xxyy_yyz_0, ta2_zz_xxyy_yzz_0, ta2_zz_xxyy_zzz_0, ta2_zz_xyy_xxy_0, ta2_zz_xyy_xxy_1, ta2_zz_xyy_xy_0, ta2_zz_xyy_xy_1, ta2_zz_xyy_xyy_0, ta2_zz_xyy_xyy_1, ta2_zz_xyy_xyz_0, ta2_zz_xyy_xyz_1, ta2_zz_xyy_yy_0, ta2_zz_xyy_yy_1, ta2_zz_xyy_yyy_0, ta2_zz_xyy_yyy_1, ta2_zz_xyy_yyz_0, ta2_zz_xyy_yyz_1, ta2_zz_xyy_yz_0, ta2_zz_xyy_yz_1, ta2_zz_xyy_yzz_0, ta2_zz_xyy_yzz_1, ta2_zz_xyy_zzz_0, ta2_zz_xyy_zzz_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyy_xxx_0[i] = ta2_zz_xx_xxx_0[i] * fe_0 - ta2_zz_xx_xxx_1[i] * fe_0 + ta2_zz_xxy_xxx_0[i] * pa_y[i] - ta2_zz_xxy_xxx_1[i] * pc_y[i];

        ta2_zz_xxyy_xxy_0[i] = ta2_zz_yy_xxy_0[i] * fe_0 - ta2_zz_yy_xxy_1[i] * fe_0 + 2.0 * ta2_zz_xyy_xy_0[i] * fe_0 - 2.0 * ta2_zz_xyy_xy_1[i] * fe_0 + ta2_zz_xyy_xxy_0[i] * pa_x[i] - ta2_zz_xyy_xxy_1[i] * pc_x[i];

        ta2_zz_xxyy_xxz_0[i] = ta2_zz_xx_xxz_0[i] * fe_0 - ta2_zz_xx_xxz_1[i] * fe_0 + ta2_zz_xxy_xxz_0[i] * pa_y[i] - ta2_zz_xxy_xxz_1[i] * pc_y[i];

        ta2_zz_xxyy_xyy_0[i] = ta2_zz_yy_xyy_0[i] * fe_0 - ta2_zz_yy_xyy_1[i] * fe_0 + ta2_zz_xyy_yy_0[i] * fe_0 - ta2_zz_xyy_yy_1[i] * fe_0 + ta2_zz_xyy_xyy_0[i] * pa_x[i] - ta2_zz_xyy_xyy_1[i] * pc_x[i];

        ta2_zz_xxyy_xyz_0[i] = ta2_zz_yy_xyz_0[i] * fe_0 - ta2_zz_yy_xyz_1[i] * fe_0 + ta2_zz_xyy_yz_0[i] * fe_0 - ta2_zz_xyy_yz_1[i] * fe_0 + ta2_zz_xyy_xyz_0[i] * pa_x[i] - ta2_zz_xyy_xyz_1[i] * pc_x[i];

        ta2_zz_xxyy_xzz_0[i] = ta2_zz_xx_xzz_0[i] * fe_0 - ta2_zz_xx_xzz_1[i] * fe_0 + ta2_zz_xxy_xzz_0[i] * pa_y[i] - ta2_zz_xxy_xzz_1[i] * pc_y[i];

        ta2_zz_xxyy_yyy_0[i] = ta2_zz_yy_yyy_0[i] * fe_0 - ta2_zz_yy_yyy_1[i] * fe_0 + ta2_zz_xyy_yyy_0[i] * pa_x[i] - ta2_zz_xyy_yyy_1[i] * pc_x[i];

        ta2_zz_xxyy_yyz_0[i] = ta2_zz_yy_yyz_0[i] * fe_0 - ta2_zz_yy_yyz_1[i] * fe_0 + ta2_zz_xyy_yyz_0[i] * pa_x[i] - ta2_zz_xyy_yyz_1[i] * pc_x[i];

        ta2_zz_xxyy_yzz_0[i] = ta2_zz_yy_yzz_0[i] * fe_0 - ta2_zz_yy_yzz_1[i] * fe_0 + ta2_zz_xyy_yzz_0[i] * pa_x[i] - ta2_zz_xyy_yzz_1[i] * pc_x[i];

        ta2_zz_xxyy_zzz_0[i] = ta2_zz_yy_zzz_0[i] * fe_0 - ta2_zz_yy_zzz_1[i] * fe_0 + ta2_zz_xyy_zzz_0[i] * pa_x[i] - ta2_zz_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 790-800 components of targeted buffer : GF

    auto ta2_zz_xxyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 790);

    auto ta2_zz_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 791);

    auto ta2_zz_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 792);

    auto ta2_zz_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 793);

    auto ta2_zz_xxyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 794);

    auto ta2_zz_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 795);

    auto ta2_zz_xxyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 796);

    auto ta2_zz_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 797);

    auto ta2_zz_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 798);

    auto ta2_zz_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 799);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxy_xxy_1, ta1_z_xxy_xyy_1, ta1_z_xxy_yyy_1, ta2_zz_xxy_xxy_0, ta2_zz_xxy_xxy_1, ta2_zz_xxy_xyy_0, ta2_zz_xxy_xyy_1, ta2_zz_xxy_yyy_0, ta2_zz_xxy_yyy_1, ta2_zz_xxyz_xxx_0, ta2_zz_xxyz_xxy_0, ta2_zz_xxyz_xxz_0, ta2_zz_xxyz_xyy_0, ta2_zz_xxyz_xyz_0, ta2_zz_xxyz_xzz_0, ta2_zz_xxyz_yyy_0, ta2_zz_xxyz_yyz_0, ta2_zz_xxyz_yzz_0, ta2_zz_xxyz_zzz_0, ta2_zz_xxz_xxx_0, ta2_zz_xxz_xxx_1, ta2_zz_xxz_xxz_0, ta2_zz_xxz_xxz_1, ta2_zz_xxz_xyz_0, ta2_zz_xxz_xyz_1, ta2_zz_xxz_xz_0, ta2_zz_xxz_xz_1, ta2_zz_xxz_xzz_0, ta2_zz_xxz_xzz_1, ta2_zz_xxz_zzz_0, ta2_zz_xxz_zzz_1, ta2_zz_xyz_yyz_0, ta2_zz_xyz_yyz_1, ta2_zz_xyz_yzz_0, ta2_zz_xyz_yzz_1, ta2_zz_yz_yyz_0, ta2_zz_yz_yyz_1, ta2_zz_yz_yzz_0, ta2_zz_yz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyz_xxx_0[i] = ta2_zz_xxz_xxx_0[i] * pa_y[i] - ta2_zz_xxz_xxx_1[i] * pc_y[i];

        ta2_zz_xxyz_xxy_0[i] = 2.0 * ta1_z_xxy_xxy_1[i] + ta2_zz_xxy_xxy_0[i] * pa_z[i] - ta2_zz_xxy_xxy_1[i] * pc_z[i];

        ta2_zz_xxyz_xxz_0[i] = ta2_zz_xxz_xxz_0[i] * pa_y[i] - ta2_zz_xxz_xxz_1[i] * pc_y[i];

        ta2_zz_xxyz_xyy_0[i] = 2.0 * ta1_z_xxy_xyy_1[i] + ta2_zz_xxy_xyy_0[i] * pa_z[i] - ta2_zz_xxy_xyy_1[i] * pc_z[i];

        ta2_zz_xxyz_xyz_0[i] = ta2_zz_xxz_xz_0[i] * fe_0 - ta2_zz_xxz_xz_1[i] * fe_0 + ta2_zz_xxz_xyz_0[i] * pa_y[i] - ta2_zz_xxz_xyz_1[i] * pc_y[i];

        ta2_zz_xxyz_xzz_0[i] = ta2_zz_xxz_xzz_0[i] * pa_y[i] - ta2_zz_xxz_xzz_1[i] * pc_y[i];

        ta2_zz_xxyz_yyy_0[i] = 2.0 * ta1_z_xxy_yyy_1[i] + ta2_zz_xxy_yyy_0[i] * pa_z[i] - ta2_zz_xxy_yyy_1[i] * pc_z[i];

        ta2_zz_xxyz_yyz_0[i] = ta2_zz_yz_yyz_0[i] * fe_0 - ta2_zz_yz_yyz_1[i] * fe_0 + ta2_zz_xyz_yyz_0[i] * pa_x[i] - ta2_zz_xyz_yyz_1[i] * pc_x[i];

        ta2_zz_xxyz_yzz_0[i] = ta2_zz_yz_yzz_0[i] * fe_0 - ta2_zz_yz_yzz_1[i] * fe_0 + ta2_zz_xyz_yzz_0[i] * pa_x[i] - ta2_zz_xyz_yzz_1[i] * pc_x[i];

        ta2_zz_xxyz_zzz_0[i] = ta2_zz_xxz_zzz_0[i] * pa_y[i] - ta2_zz_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 800-810 components of targeted buffer : GF

    auto ta2_zz_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 800);

    auto ta2_zz_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 801);

    auto ta2_zz_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 802);

    auto ta2_zz_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 803);

    auto ta2_zz_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 804);

    auto ta2_zz_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 805);

    auto ta2_zz_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 806);

    auto ta2_zz_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 807);

    auto ta2_zz_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 808);

    auto ta2_zz_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 809);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxz_xxx_1, ta1_z_xxz_xxy_1, ta1_z_xxz_xyy_1, ta2_zz_xx_xxx_0, ta2_zz_xx_xxx_1, ta2_zz_xx_xxy_0, ta2_zz_xx_xxy_1, ta2_zz_xx_xyy_0, ta2_zz_xx_xyy_1, ta2_zz_xxz_xxx_0, ta2_zz_xxz_xxx_1, ta2_zz_xxz_xxy_0, ta2_zz_xxz_xxy_1, ta2_zz_xxz_xyy_0, ta2_zz_xxz_xyy_1, ta2_zz_xxzz_xxx_0, ta2_zz_xxzz_xxy_0, ta2_zz_xxzz_xxz_0, ta2_zz_xxzz_xyy_0, ta2_zz_xxzz_xyz_0, ta2_zz_xxzz_xzz_0, ta2_zz_xxzz_yyy_0, ta2_zz_xxzz_yyz_0, ta2_zz_xxzz_yzz_0, ta2_zz_xxzz_zzz_0, ta2_zz_xzz_xxz_0, ta2_zz_xzz_xxz_1, ta2_zz_xzz_xyz_0, ta2_zz_xzz_xyz_1, ta2_zz_xzz_xz_0, ta2_zz_xzz_xz_1, ta2_zz_xzz_xzz_0, ta2_zz_xzz_xzz_1, ta2_zz_xzz_yyy_0, ta2_zz_xzz_yyy_1, ta2_zz_xzz_yyz_0, ta2_zz_xzz_yyz_1, ta2_zz_xzz_yz_0, ta2_zz_xzz_yz_1, ta2_zz_xzz_yzz_0, ta2_zz_xzz_yzz_1, ta2_zz_xzz_zz_0, ta2_zz_xzz_zz_1, ta2_zz_xzz_zzz_0, ta2_zz_xzz_zzz_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxzz_xxx_0[i] = ta2_zz_xx_xxx_0[i] * fe_0 - ta2_zz_xx_xxx_1[i] * fe_0 + 2.0 * ta1_z_xxz_xxx_1[i] + ta2_zz_xxz_xxx_0[i] * pa_z[i] - ta2_zz_xxz_xxx_1[i] * pc_z[i];

        ta2_zz_xxzz_xxy_0[i] = ta2_zz_xx_xxy_0[i] * fe_0 - ta2_zz_xx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xxy_1[i] + ta2_zz_xxz_xxy_0[i] * pa_z[i] - ta2_zz_xxz_xxy_1[i] * pc_z[i];

        ta2_zz_xxzz_xxz_0[i] = ta2_zz_zz_xxz_0[i] * fe_0 - ta2_zz_zz_xxz_1[i] * fe_0 + 2.0 * ta2_zz_xzz_xz_0[i] * fe_0 - 2.0 * ta2_zz_xzz_xz_1[i] * fe_0 + ta2_zz_xzz_xxz_0[i] * pa_x[i] - ta2_zz_xzz_xxz_1[i] * pc_x[i];

        ta2_zz_xxzz_xyy_0[i] = ta2_zz_xx_xyy_0[i] * fe_0 - ta2_zz_xx_xyy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xyy_1[i] + ta2_zz_xxz_xyy_0[i] * pa_z[i] - ta2_zz_xxz_xyy_1[i] * pc_z[i];

        ta2_zz_xxzz_xyz_0[i] = ta2_zz_zz_xyz_0[i] * fe_0 - ta2_zz_zz_xyz_1[i] * fe_0 + ta2_zz_xzz_yz_0[i] * fe_0 - ta2_zz_xzz_yz_1[i] * fe_0 + ta2_zz_xzz_xyz_0[i] * pa_x[i] - ta2_zz_xzz_xyz_1[i] * pc_x[i];

        ta2_zz_xxzz_xzz_0[i] = ta2_zz_zz_xzz_0[i] * fe_0 - ta2_zz_zz_xzz_1[i] * fe_0 + ta2_zz_xzz_zz_0[i] * fe_0 - ta2_zz_xzz_zz_1[i] * fe_0 + ta2_zz_xzz_xzz_0[i] * pa_x[i] - ta2_zz_xzz_xzz_1[i] * pc_x[i];

        ta2_zz_xxzz_yyy_0[i] = ta2_zz_zz_yyy_0[i] * fe_0 - ta2_zz_zz_yyy_1[i] * fe_0 + ta2_zz_xzz_yyy_0[i] * pa_x[i] - ta2_zz_xzz_yyy_1[i] * pc_x[i];

        ta2_zz_xxzz_yyz_0[i] = ta2_zz_zz_yyz_0[i] * fe_0 - ta2_zz_zz_yyz_1[i] * fe_0 + ta2_zz_xzz_yyz_0[i] * pa_x[i] - ta2_zz_xzz_yyz_1[i] * pc_x[i];

        ta2_zz_xxzz_yzz_0[i] = ta2_zz_zz_yzz_0[i] * fe_0 - ta2_zz_zz_yzz_1[i] * fe_0 + ta2_zz_xzz_yzz_0[i] * pa_x[i] - ta2_zz_xzz_yzz_1[i] * pc_x[i];

        ta2_zz_xxzz_zzz_0[i] = ta2_zz_zz_zzz_0[i] * fe_0 - ta2_zz_zz_zzz_1[i] * fe_0 + ta2_zz_xzz_zzz_0[i] * pa_x[i] - ta2_zz_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 810-820 components of targeted buffer : GF

    auto ta2_zz_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 810);

    auto ta2_zz_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 811);

    auto ta2_zz_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 812);

    auto ta2_zz_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 813);

    auto ta2_zz_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 814);

    auto ta2_zz_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 815);

    auto ta2_zz_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 816);

    auto ta2_zz_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 817);

    auto ta2_zz_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 818);

    auto ta2_zz_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 819);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xyyy_xxx_0, ta2_zz_xyyy_xxy_0, ta2_zz_xyyy_xxz_0, ta2_zz_xyyy_xyy_0, ta2_zz_xyyy_xyz_0, ta2_zz_xyyy_xzz_0, ta2_zz_xyyy_yyy_0, ta2_zz_xyyy_yyz_0, ta2_zz_xyyy_yzz_0, ta2_zz_xyyy_zzz_0, ta2_zz_yyy_xx_0, ta2_zz_yyy_xx_1, ta2_zz_yyy_xxx_0, ta2_zz_yyy_xxx_1, ta2_zz_yyy_xxy_0, ta2_zz_yyy_xxy_1, ta2_zz_yyy_xxz_0, ta2_zz_yyy_xxz_1, ta2_zz_yyy_xy_0, ta2_zz_yyy_xy_1, ta2_zz_yyy_xyy_0, ta2_zz_yyy_xyy_1, ta2_zz_yyy_xyz_0, ta2_zz_yyy_xyz_1, ta2_zz_yyy_xz_0, ta2_zz_yyy_xz_1, ta2_zz_yyy_xzz_0, ta2_zz_yyy_xzz_1, ta2_zz_yyy_yy_0, ta2_zz_yyy_yy_1, ta2_zz_yyy_yyy_0, ta2_zz_yyy_yyy_1, ta2_zz_yyy_yyz_0, ta2_zz_yyy_yyz_1, ta2_zz_yyy_yz_0, ta2_zz_yyy_yz_1, ta2_zz_yyy_yzz_0, ta2_zz_yyy_yzz_1, ta2_zz_yyy_zz_0, ta2_zz_yyy_zz_1, ta2_zz_yyy_zzz_0, ta2_zz_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyy_xxx_0[i] = 3.0 * ta2_zz_yyy_xx_0[i] * fe_0 - 3.0 * ta2_zz_yyy_xx_1[i] * fe_0 + ta2_zz_yyy_xxx_0[i] * pa_x[i] - ta2_zz_yyy_xxx_1[i] * pc_x[i];

        ta2_zz_xyyy_xxy_0[i] = 2.0 * ta2_zz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xy_1[i] * fe_0 + ta2_zz_yyy_xxy_0[i] * pa_x[i] - ta2_zz_yyy_xxy_1[i] * pc_x[i];

        ta2_zz_xyyy_xxz_0[i] = 2.0 * ta2_zz_yyy_xz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xz_1[i] * fe_0 + ta2_zz_yyy_xxz_0[i] * pa_x[i] - ta2_zz_yyy_xxz_1[i] * pc_x[i];

        ta2_zz_xyyy_xyy_0[i] = ta2_zz_yyy_yy_0[i] * fe_0 - ta2_zz_yyy_yy_1[i] * fe_0 + ta2_zz_yyy_xyy_0[i] * pa_x[i] - ta2_zz_yyy_xyy_1[i] * pc_x[i];

        ta2_zz_xyyy_xyz_0[i] = ta2_zz_yyy_yz_0[i] * fe_0 - ta2_zz_yyy_yz_1[i] * fe_0 + ta2_zz_yyy_xyz_0[i] * pa_x[i] - ta2_zz_yyy_xyz_1[i] * pc_x[i];

        ta2_zz_xyyy_xzz_0[i] = ta2_zz_yyy_zz_0[i] * fe_0 - ta2_zz_yyy_zz_1[i] * fe_0 + ta2_zz_yyy_xzz_0[i] * pa_x[i] - ta2_zz_yyy_xzz_1[i] * pc_x[i];

        ta2_zz_xyyy_yyy_0[i] = ta2_zz_yyy_yyy_0[i] * pa_x[i] - ta2_zz_yyy_yyy_1[i] * pc_x[i];

        ta2_zz_xyyy_yyz_0[i] = ta2_zz_yyy_yyz_0[i] * pa_x[i] - ta2_zz_yyy_yyz_1[i] * pc_x[i];

        ta2_zz_xyyy_yzz_0[i] = ta2_zz_yyy_yzz_0[i] * pa_x[i] - ta2_zz_yyy_yzz_1[i] * pc_x[i];

        ta2_zz_xyyy_zzz_0[i] = ta2_zz_yyy_zzz_0[i] * pa_x[i] - ta2_zz_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 820-830 components of targeted buffer : GF

    auto ta2_zz_xyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 820);

    auto ta2_zz_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 821);

    auto ta2_zz_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 822);

    auto ta2_zz_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 823);

    auto ta2_zz_xyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 824);

    auto ta2_zz_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 825);

    auto ta2_zz_xyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 826);

    auto ta2_zz_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 827);

    auto ta2_zz_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 828);

    auto ta2_zz_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 829);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xyy_xxx_1, ta1_z_xyy_xxy_1, ta1_z_xyy_xyy_1, ta2_zz_xyy_xxx_0, ta2_zz_xyy_xxx_1, ta2_zz_xyy_xxy_0, ta2_zz_xyy_xxy_1, ta2_zz_xyy_xyy_0, ta2_zz_xyy_xyy_1, ta2_zz_xyyz_xxx_0, ta2_zz_xyyz_xxy_0, ta2_zz_xyyz_xxz_0, ta2_zz_xyyz_xyy_0, ta2_zz_xyyz_xyz_0, ta2_zz_xyyz_xzz_0, ta2_zz_xyyz_yyy_0, ta2_zz_xyyz_yyz_0, ta2_zz_xyyz_yzz_0, ta2_zz_xyyz_zzz_0, ta2_zz_yyz_xxz_0, ta2_zz_yyz_xxz_1, ta2_zz_yyz_xyz_0, ta2_zz_yyz_xyz_1, ta2_zz_yyz_xz_0, ta2_zz_yyz_xz_1, ta2_zz_yyz_xzz_0, ta2_zz_yyz_xzz_1, ta2_zz_yyz_yyy_0, ta2_zz_yyz_yyy_1, ta2_zz_yyz_yyz_0, ta2_zz_yyz_yyz_1, ta2_zz_yyz_yz_0, ta2_zz_yyz_yz_1, ta2_zz_yyz_yzz_0, ta2_zz_yyz_yzz_1, ta2_zz_yyz_zz_0, ta2_zz_yyz_zz_1, ta2_zz_yyz_zzz_0, ta2_zz_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyz_xxx_0[i] = 2.0 * ta1_z_xyy_xxx_1[i] + ta2_zz_xyy_xxx_0[i] * pa_z[i] - ta2_zz_xyy_xxx_1[i] * pc_z[i];

        ta2_zz_xyyz_xxy_0[i] = 2.0 * ta1_z_xyy_xxy_1[i] + ta2_zz_xyy_xxy_0[i] * pa_z[i] - ta2_zz_xyy_xxy_1[i] * pc_z[i];

        ta2_zz_xyyz_xxz_0[i] = 2.0 * ta2_zz_yyz_xz_0[i] * fe_0 - 2.0 * ta2_zz_yyz_xz_1[i] * fe_0 + ta2_zz_yyz_xxz_0[i] * pa_x[i] - ta2_zz_yyz_xxz_1[i] * pc_x[i];

        ta2_zz_xyyz_xyy_0[i] = 2.0 * ta1_z_xyy_xyy_1[i] + ta2_zz_xyy_xyy_0[i] * pa_z[i] - ta2_zz_xyy_xyy_1[i] * pc_z[i];

        ta2_zz_xyyz_xyz_0[i] = ta2_zz_yyz_yz_0[i] * fe_0 - ta2_zz_yyz_yz_1[i] * fe_0 + ta2_zz_yyz_xyz_0[i] * pa_x[i] - ta2_zz_yyz_xyz_1[i] * pc_x[i];

        ta2_zz_xyyz_xzz_0[i] = ta2_zz_yyz_zz_0[i] * fe_0 - ta2_zz_yyz_zz_1[i] * fe_0 + ta2_zz_yyz_xzz_0[i] * pa_x[i] - ta2_zz_yyz_xzz_1[i] * pc_x[i];

        ta2_zz_xyyz_yyy_0[i] = ta2_zz_yyz_yyy_0[i] * pa_x[i] - ta2_zz_yyz_yyy_1[i] * pc_x[i];

        ta2_zz_xyyz_yyz_0[i] = ta2_zz_yyz_yyz_0[i] * pa_x[i] - ta2_zz_yyz_yyz_1[i] * pc_x[i];

        ta2_zz_xyyz_yzz_0[i] = ta2_zz_yyz_yzz_0[i] * pa_x[i] - ta2_zz_yyz_yzz_1[i] * pc_x[i];

        ta2_zz_xyyz_zzz_0[i] = ta2_zz_yyz_zzz_0[i] * pa_x[i] - ta2_zz_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 830-840 components of targeted buffer : GF

    auto ta2_zz_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 830);

    auto ta2_zz_xyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 831);

    auto ta2_zz_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 832);

    auto ta2_zz_xyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 833);

    auto ta2_zz_xyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 834);

    auto ta2_zz_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 835);

    auto ta2_zz_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 836);

    auto ta2_zz_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 837);

    auto ta2_zz_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 838);

    auto ta2_zz_xyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 839);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xyzz_xxx_0, ta2_zz_xyzz_xxy_0, ta2_zz_xyzz_xxz_0, ta2_zz_xyzz_xyy_0, ta2_zz_xyzz_xyz_0, ta2_zz_xyzz_xzz_0, ta2_zz_xyzz_yyy_0, ta2_zz_xyzz_yyz_0, ta2_zz_xyzz_yzz_0, ta2_zz_xyzz_zzz_0, ta2_zz_xzz_xxx_0, ta2_zz_xzz_xxx_1, ta2_zz_xzz_xxz_0, ta2_zz_xzz_xxz_1, ta2_zz_xzz_xzz_0, ta2_zz_xzz_xzz_1, ta2_zz_yzz_xxy_0, ta2_zz_yzz_xxy_1, ta2_zz_yzz_xy_0, ta2_zz_yzz_xy_1, ta2_zz_yzz_xyy_0, ta2_zz_yzz_xyy_1, ta2_zz_yzz_xyz_0, ta2_zz_yzz_xyz_1, ta2_zz_yzz_yy_0, ta2_zz_yzz_yy_1, ta2_zz_yzz_yyy_0, ta2_zz_yzz_yyy_1, ta2_zz_yzz_yyz_0, ta2_zz_yzz_yyz_1, ta2_zz_yzz_yz_0, ta2_zz_yzz_yz_1, ta2_zz_yzz_yzz_0, ta2_zz_yzz_yzz_1, ta2_zz_yzz_zzz_0, ta2_zz_yzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyzz_xxx_0[i] = ta2_zz_xzz_xxx_0[i] * pa_y[i] - ta2_zz_xzz_xxx_1[i] * pc_y[i];

        ta2_zz_xyzz_xxy_0[i] = 2.0 * ta2_zz_yzz_xy_0[i] * fe_0 - 2.0 * ta2_zz_yzz_xy_1[i] * fe_0 + ta2_zz_yzz_xxy_0[i] * pa_x[i] - ta2_zz_yzz_xxy_1[i] * pc_x[i];

        ta2_zz_xyzz_xxz_0[i] = ta2_zz_xzz_xxz_0[i] * pa_y[i] - ta2_zz_xzz_xxz_1[i] * pc_y[i];

        ta2_zz_xyzz_xyy_0[i] = ta2_zz_yzz_yy_0[i] * fe_0 - ta2_zz_yzz_yy_1[i] * fe_0 + ta2_zz_yzz_xyy_0[i] * pa_x[i] - ta2_zz_yzz_xyy_1[i] * pc_x[i];

        ta2_zz_xyzz_xyz_0[i] = ta2_zz_yzz_yz_0[i] * fe_0 - ta2_zz_yzz_yz_1[i] * fe_0 + ta2_zz_yzz_xyz_0[i] * pa_x[i] - ta2_zz_yzz_xyz_1[i] * pc_x[i];

        ta2_zz_xyzz_xzz_0[i] = ta2_zz_xzz_xzz_0[i] * pa_y[i] - ta2_zz_xzz_xzz_1[i] * pc_y[i];

        ta2_zz_xyzz_yyy_0[i] = ta2_zz_yzz_yyy_0[i] * pa_x[i] - ta2_zz_yzz_yyy_1[i] * pc_x[i];

        ta2_zz_xyzz_yyz_0[i] = ta2_zz_yzz_yyz_0[i] * pa_x[i] - ta2_zz_yzz_yyz_1[i] * pc_x[i];

        ta2_zz_xyzz_yzz_0[i] = ta2_zz_yzz_yzz_0[i] * pa_x[i] - ta2_zz_yzz_yzz_1[i] * pc_x[i];

        ta2_zz_xyzz_zzz_0[i] = ta2_zz_yzz_zzz_0[i] * pa_x[i] - ta2_zz_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 840-850 components of targeted buffer : GF

    auto ta2_zz_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 840);

    auto ta2_zz_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 841);

    auto ta2_zz_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 842);

    auto ta2_zz_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 843);

    auto ta2_zz_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 844);

    auto ta2_zz_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 845);

    auto ta2_zz_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 846);

    auto ta2_zz_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 847);

    auto ta2_zz_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 848);

    auto ta2_zz_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 849);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xzzz_xxx_0, ta2_zz_xzzz_xxy_0, ta2_zz_xzzz_xxz_0, ta2_zz_xzzz_xyy_0, ta2_zz_xzzz_xyz_0, ta2_zz_xzzz_xzz_0, ta2_zz_xzzz_yyy_0, ta2_zz_xzzz_yyz_0, ta2_zz_xzzz_yzz_0, ta2_zz_xzzz_zzz_0, ta2_zz_zzz_xx_0, ta2_zz_zzz_xx_1, ta2_zz_zzz_xxx_0, ta2_zz_zzz_xxx_1, ta2_zz_zzz_xxy_0, ta2_zz_zzz_xxy_1, ta2_zz_zzz_xxz_0, ta2_zz_zzz_xxz_1, ta2_zz_zzz_xy_0, ta2_zz_zzz_xy_1, ta2_zz_zzz_xyy_0, ta2_zz_zzz_xyy_1, ta2_zz_zzz_xyz_0, ta2_zz_zzz_xyz_1, ta2_zz_zzz_xz_0, ta2_zz_zzz_xz_1, ta2_zz_zzz_xzz_0, ta2_zz_zzz_xzz_1, ta2_zz_zzz_yy_0, ta2_zz_zzz_yy_1, ta2_zz_zzz_yyy_0, ta2_zz_zzz_yyy_1, ta2_zz_zzz_yyz_0, ta2_zz_zzz_yyz_1, ta2_zz_zzz_yz_0, ta2_zz_zzz_yz_1, ta2_zz_zzz_yzz_0, ta2_zz_zzz_yzz_1, ta2_zz_zzz_zz_0, ta2_zz_zzz_zz_1, ta2_zz_zzz_zzz_0, ta2_zz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzzz_xxx_0[i] = 3.0 * ta2_zz_zzz_xx_0[i] * fe_0 - 3.0 * ta2_zz_zzz_xx_1[i] * fe_0 + ta2_zz_zzz_xxx_0[i] * pa_x[i] - ta2_zz_zzz_xxx_1[i] * pc_x[i];

        ta2_zz_xzzz_xxy_0[i] = 2.0 * ta2_zz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xy_1[i] * fe_0 + ta2_zz_zzz_xxy_0[i] * pa_x[i] - ta2_zz_zzz_xxy_1[i] * pc_x[i];

        ta2_zz_xzzz_xxz_0[i] = 2.0 * ta2_zz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xz_1[i] * fe_0 + ta2_zz_zzz_xxz_0[i] * pa_x[i] - ta2_zz_zzz_xxz_1[i] * pc_x[i];

        ta2_zz_xzzz_xyy_0[i] = ta2_zz_zzz_yy_0[i] * fe_0 - ta2_zz_zzz_yy_1[i] * fe_0 + ta2_zz_zzz_xyy_0[i] * pa_x[i] - ta2_zz_zzz_xyy_1[i] * pc_x[i];

        ta2_zz_xzzz_xyz_0[i] = ta2_zz_zzz_yz_0[i] * fe_0 - ta2_zz_zzz_yz_1[i] * fe_0 + ta2_zz_zzz_xyz_0[i] * pa_x[i] - ta2_zz_zzz_xyz_1[i] * pc_x[i];

        ta2_zz_xzzz_xzz_0[i] = ta2_zz_zzz_zz_0[i] * fe_0 - ta2_zz_zzz_zz_1[i] * fe_0 + ta2_zz_zzz_xzz_0[i] * pa_x[i] - ta2_zz_zzz_xzz_1[i] * pc_x[i];

        ta2_zz_xzzz_yyy_0[i] = ta2_zz_zzz_yyy_0[i] * pa_x[i] - ta2_zz_zzz_yyy_1[i] * pc_x[i];

        ta2_zz_xzzz_yyz_0[i] = ta2_zz_zzz_yyz_0[i] * pa_x[i] - ta2_zz_zzz_yyz_1[i] * pc_x[i];

        ta2_zz_xzzz_yzz_0[i] = ta2_zz_zzz_yzz_0[i] * pa_x[i] - ta2_zz_zzz_yzz_1[i] * pc_x[i];

        ta2_zz_xzzz_zzz_0[i] = ta2_zz_zzz_zzz_0[i] * pa_x[i] - ta2_zz_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 850-860 components of targeted buffer : GF

    auto ta2_zz_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 850);

    auto ta2_zz_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 851);

    auto ta2_zz_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 852);

    auto ta2_zz_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 853);

    auto ta2_zz_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 854);

    auto ta2_zz_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 855);

    auto ta2_zz_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 856);

    auto ta2_zz_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 857);

    auto ta2_zz_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 858);

    auto ta2_zz_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 859);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_yy_xxx_0, ta2_zz_yy_xxx_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xxz_0, ta2_zz_yy_xxz_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_xyz_0, ta2_zz_yy_xyz_1, ta2_zz_yy_xzz_0, ta2_zz_yy_xzz_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yy_yyz_0, ta2_zz_yy_yyz_1, ta2_zz_yy_yzz_0, ta2_zz_yy_yzz_1, ta2_zz_yy_zzz_0, ta2_zz_yy_zzz_1, ta2_zz_yyy_xx_0, ta2_zz_yyy_xx_1, ta2_zz_yyy_xxx_0, ta2_zz_yyy_xxx_1, ta2_zz_yyy_xxy_0, ta2_zz_yyy_xxy_1, ta2_zz_yyy_xxz_0, ta2_zz_yyy_xxz_1, ta2_zz_yyy_xy_0, ta2_zz_yyy_xy_1, ta2_zz_yyy_xyy_0, ta2_zz_yyy_xyy_1, ta2_zz_yyy_xyz_0, ta2_zz_yyy_xyz_1, ta2_zz_yyy_xz_0, ta2_zz_yyy_xz_1, ta2_zz_yyy_xzz_0, ta2_zz_yyy_xzz_1, ta2_zz_yyy_yy_0, ta2_zz_yyy_yy_1, ta2_zz_yyy_yyy_0, ta2_zz_yyy_yyy_1, ta2_zz_yyy_yyz_0, ta2_zz_yyy_yyz_1, ta2_zz_yyy_yz_0, ta2_zz_yyy_yz_1, ta2_zz_yyy_yzz_0, ta2_zz_yyy_yzz_1, ta2_zz_yyy_zz_0, ta2_zz_yyy_zz_1, ta2_zz_yyy_zzz_0, ta2_zz_yyy_zzz_1, ta2_zz_yyyy_xxx_0, ta2_zz_yyyy_xxy_0, ta2_zz_yyyy_xxz_0, ta2_zz_yyyy_xyy_0, ta2_zz_yyyy_xyz_0, ta2_zz_yyyy_xzz_0, ta2_zz_yyyy_yyy_0, ta2_zz_yyyy_yyz_0, ta2_zz_yyyy_yzz_0, ta2_zz_yyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyy_xxx_0[i] = 3.0 * ta2_zz_yy_xxx_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxx_1[i] * fe_0 + ta2_zz_yyy_xxx_0[i] * pa_y[i] - ta2_zz_yyy_xxx_1[i] * pc_y[i];

        ta2_zz_yyyy_xxy_0[i] = 3.0 * ta2_zz_yy_xxy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxy_1[i] * fe_0 + ta2_zz_yyy_xx_0[i] * fe_0 - ta2_zz_yyy_xx_1[i] * fe_0 + ta2_zz_yyy_xxy_0[i] * pa_y[i] - ta2_zz_yyy_xxy_1[i] * pc_y[i];

        ta2_zz_yyyy_xxz_0[i] = 3.0 * ta2_zz_yy_xxz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxz_1[i] * fe_0 + ta2_zz_yyy_xxz_0[i] * pa_y[i] - ta2_zz_yyy_xxz_1[i] * pc_y[i];

        ta2_zz_yyyy_xyy_0[i] = 3.0 * ta2_zz_yy_xyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyy_1[i] * fe_0 + 2.0 * ta2_zz_yyy_xy_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xy_1[i] * fe_0 + ta2_zz_yyy_xyy_0[i] * pa_y[i] - ta2_zz_yyy_xyy_1[i] * pc_y[i];

        ta2_zz_yyyy_xyz_0[i] = 3.0 * ta2_zz_yy_xyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyz_1[i] * fe_0 + ta2_zz_yyy_xz_0[i] * fe_0 - ta2_zz_yyy_xz_1[i] * fe_0 + ta2_zz_yyy_xyz_0[i] * pa_y[i] - ta2_zz_yyy_xyz_1[i] * pc_y[i];

        ta2_zz_yyyy_xzz_0[i] = 3.0 * ta2_zz_yy_xzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xzz_1[i] * fe_0 + ta2_zz_yyy_xzz_0[i] * pa_y[i] - ta2_zz_yyy_xzz_1[i] * pc_y[i];

        ta2_zz_yyyy_yyy_0[i] = 3.0 * ta2_zz_yy_yyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyy_1[i] * fe_0 + 3.0 * ta2_zz_yyy_yy_0[i] * fe_0 - 3.0 * ta2_zz_yyy_yy_1[i] * fe_0 + ta2_zz_yyy_yyy_0[i] * pa_y[i] - ta2_zz_yyy_yyy_1[i] * pc_y[i];

        ta2_zz_yyyy_yyz_0[i] = 3.0 * ta2_zz_yy_yyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyz_1[i] * fe_0 + 2.0 * ta2_zz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_yz_1[i] * fe_0 + ta2_zz_yyy_yyz_0[i] * pa_y[i] - ta2_zz_yyy_yyz_1[i] * pc_y[i];

        ta2_zz_yyyy_yzz_0[i] = 3.0 * ta2_zz_yy_yzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yzz_1[i] * fe_0 + ta2_zz_yyy_zz_0[i] * fe_0 - ta2_zz_yyy_zz_1[i] * fe_0 + ta2_zz_yyy_yzz_0[i] * pa_y[i] - ta2_zz_yyy_yzz_1[i] * pc_y[i];

        ta2_zz_yyyy_zzz_0[i] = 3.0 * ta2_zz_yy_zzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_zzz_1[i] * fe_0 + ta2_zz_yyy_zzz_0[i] * pa_y[i] - ta2_zz_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 860-870 components of targeted buffer : GF

    auto ta2_zz_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 860);

    auto ta2_zz_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 861);

    auto ta2_zz_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 862);

    auto ta2_zz_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 863);

    auto ta2_zz_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 864);

    auto ta2_zz_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 865);

    auto ta2_zz_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 866);

    auto ta2_zz_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 867);

    auto ta2_zz_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 868);

    auto ta2_zz_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 869);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyy_xxx_1, ta1_z_yyy_xxy_1, ta1_z_yyy_xyy_1, ta1_z_yyy_xyz_1, ta1_z_yyy_yyy_1, ta1_z_yyy_yyz_1, ta1_z_yyy_yzz_1, ta2_zz_yyy_xxx_0, ta2_zz_yyy_xxx_1, ta2_zz_yyy_xxy_0, ta2_zz_yyy_xxy_1, ta2_zz_yyy_xy_0, ta2_zz_yyy_xy_1, ta2_zz_yyy_xyy_0, ta2_zz_yyy_xyy_1, ta2_zz_yyy_xyz_0, ta2_zz_yyy_xyz_1, ta2_zz_yyy_yy_0, ta2_zz_yyy_yy_1, ta2_zz_yyy_yyy_0, ta2_zz_yyy_yyy_1, ta2_zz_yyy_yyz_0, ta2_zz_yyy_yyz_1, ta2_zz_yyy_yz_0, ta2_zz_yyy_yz_1, ta2_zz_yyy_yzz_0, ta2_zz_yyy_yzz_1, ta2_zz_yyyz_xxx_0, ta2_zz_yyyz_xxy_0, ta2_zz_yyyz_xxz_0, ta2_zz_yyyz_xyy_0, ta2_zz_yyyz_xyz_0, ta2_zz_yyyz_xzz_0, ta2_zz_yyyz_yyy_0, ta2_zz_yyyz_yyz_0, ta2_zz_yyyz_yzz_0, ta2_zz_yyyz_zzz_0, ta2_zz_yyz_xxz_0, ta2_zz_yyz_xxz_1, ta2_zz_yyz_xzz_0, ta2_zz_yyz_xzz_1, ta2_zz_yyz_zzz_0, ta2_zz_yyz_zzz_1, ta2_zz_yz_xxz_0, ta2_zz_yz_xxz_1, ta2_zz_yz_xzz_0, ta2_zz_yz_xzz_1, ta2_zz_yz_zzz_0, ta2_zz_yz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyz_xxx_0[i] = 2.0 * ta1_z_yyy_xxx_1[i] + ta2_zz_yyy_xxx_0[i] * pa_z[i] - ta2_zz_yyy_xxx_1[i] * pc_z[i];

        ta2_zz_yyyz_xxy_0[i] = 2.0 * ta1_z_yyy_xxy_1[i] + ta2_zz_yyy_xxy_0[i] * pa_z[i] - ta2_zz_yyy_xxy_1[i] * pc_z[i];

        ta2_zz_yyyz_xxz_0[i] = 2.0 * ta2_zz_yz_xxz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xxz_1[i] * fe_0 + ta2_zz_yyz_xxz_0[i] * pa_y[i] - ta2_zz_yyz_xxz_1[i] * pc_y[i];

        ta2_zz_yyyz_xyy_0[i] = 2.0 * ta1_z_yyy_xyy_1[i] + ta2_zz_yyy_xyy_0[i] * pa_z[i] - ta2_zz_yyy_xyy_1[i] * pc_z[i];

        ta2_zz_yyyz_xyz_0[i] = ta2_zz_yyy_xy_0[i] * fe_0 - ta2_zz_yyy_xy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyz_1[i] + ta2_zz_yyy_xyz_0[i] * pa_z[i] - ta2_zz_yyy_xyz_1[i] * pc_z[i];

        ta2_zz_yyyz_xzz_0[i] = 2.0 * ta2_zz_yz_xzz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xzz_1[i] * fe_0 + ta2_zz_yyz_xzz_0[i] * pa_y[i] - ta2_zz_yyz_xzz_1[i] * pc_y[i];

        ta2_zz_yyyz_yyy_0[i] = 2.0 * ta1_z_yyy_yyy_1[i] + ta2_zz_yyy_yyy_0[i] * pa_z[i] - ta2_zz_yyy_yyy_1[i] * pc_z[i];

        ta2_zz_yyyz_yyz_0[i] = ta2_zz_yyy_yy_0[i] * fe_0 - ta2_zz_yyy_yy_1[i] * fe_0 + 2.0 * ta1_z_yyy_yyz_1[i] + ta2_zz_yyy_yyz_0[i] * pa_z[i] - ta2_zz_yyy_yyz_1[i] * pc_z[i];

        ta2_zz_yyyz_yzz_0[i] = 2.0 * ta2_zz_yyy_yz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_yz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yzz_1[i] + ta2_zz_yyy_yzz_0[i] * pa_z[i] - ta2_zz_yyy_yzz_1[i] * pc_z[i];

        ta2_zz_yyyz_zzz_0[i] = 2.0 * ta2_zz_yz_zzz_0[i] * fe_0 - 2.0 * ta2_zz_yz_zzz_1[i] * fe_0 + ta2_zz_yyz_zzz_0[i] * pa_y[i] - ta2_zz_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 870-880 components of targeted buffer : GF

    auto ta2_zz_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 870);

    auto ta2_zz_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 871);

    auto ta2_zz_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 872);

    auto ta2_zz_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 873);

    auto ta2_zz_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 874);

    auto ta2_zz_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 875);

    auto ta2_zz_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 876);

    auto ta2_zz_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 877);

    auto ta2_zz_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 878);

    auto ta2_zz_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 879);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyz_xxy_1, ta1_z_yyz_xyy_1, ta1_z_yyz_yyy_1, ta2_zz_yy_xxy_0, ta2_zz_yy_xxy_1, ta2_zz_yy_xyy_0, ta2_zz_yy_xyy_1, ta2_zz_yy_yyy_0, ta2_zz_yy_yyy_1, ta2_zz_yyz_xxy_0, ta2_zz_yyz_xxy_1, ta2_zz_yyz_xyy_0, ta2_zz_yyz_xyy_1, ta2_zz_yyz_yyy_0, ta2_zz_yyz_yyy_1, ta2_zz_yyzz_xxx_0, ta2_zz_yyzz_xxy_0, ta2_zz_yyzz_xxz_0, ta2_zz_yyzz_xyy_0, ta2_zz_yyzz_xyz_0, ta2_zz_yyzz_xzz_0, ta2_zz_yyzz_yyy_0, ta2_zz_yyzz_yyz_0, ta2_zz_yyzz_yzz_0, ta2_zz_yyzz_zzz_0, ta2_zz_yzz_xxx_0, ta2_zz_yzz_xxx_1, ta2_zz_yzz_xxz_0, ta2_zz_yzz_xxz_1, ta2_zz_yzz_xyz_0, ta2_zz_yzz_xyz_1, ta2_zz_yzz_xz_0, ta2_zz_yzz_xz_1, ta2_zz_yzz_xzz_0, ta2_zz_yzz_xzz_1, ta2_zz_yzz_yyz_0, ta2_zz_yzz_yyz_1, ta2_zz_yzz_yz_0, ta2_zz_yzz_yz_1, ta2_zz_yzz_yzz_0, ta2_zz_yzz_yzz_1, ta2_zz_yzz_zz_0, ta2_zz_yzz_zz_1, ta2_zz_yzz_zzz_0, ta2_zz_yzz_zzz_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyzz_xxx_0[i] = ta2_zz_zz_xxx_0[i] * fe_0 - ta2_zz_zz_xxx_1[i] * fe_0 + ta2_zz_yzz_xxx_0[i] * pa_y[i] - ta2_zz_yzz_xxx_1[i] * pc_y[i];

        ta2_zz_yyzz_xxy_0[i] = ta2_zz_yy_xxy_0[i] * fe_0 - ta2_zz_yy_xxy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xxy_1[i] + ta2_zz_yyz_xxy_0[i] * pa_z[i] - ta2_zz_yyz_xxy_1[i] * pc_z[i];

        ta2_zz_yyzz_xxz_0[i] = ta2_zz_zz_xxz_0[i] * fe_0 - ta2_zz_zz_xxz_1[i] * fe_0 + ta2_zz_yzz_xxz_0[i] * pa_y[i] - ta2_zz_yzz_xxz_1[i] * pc_y[i];

        ta2_zz_yyzz_xyy_0[i] = ta2_zz_yy_xyy_0[i] * fe_0 - ta2_zz_yy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xyy_1[i] + ta2_zz_yyz_xyy_0[i] * pa_z[i] - ta2_zz_yyz_xyy_1[i] * pc_z[i];

        ta2_zz_yyzz_xyz_0[i] = ta2_zz_zz_xyz_0[i] * fe_0 - ta2_zz_zz_xyz_1[i] * fe_0 + ta2_zz_yzz_xz_0[i] * fe_0 - ta2_zz_yzz_xz_1[i] * fe_0 + ta2_zz_yzz_xyz_0[i] * pa_y[i] - ta2_zz_yzz_xyz_1[i] * pc_y[i];

        ta2_zz_yyzz_xzz_0[i] = ta2_zz_zz_xzz_0[i] * fe_0 - ta2_zz_zz_xzz_1[i] * fe_0 + ta2_zz_yzz_xzz_0[i] * pa_y[i] - ta2_zz_yzz_xzz_1[i] * pc_y[i];

        ta2_zz_yyzz_yyy_0[i] = ta2_zz_yy_yyy_0[i] * fe_0 - ta2_zz_yy_yyy_1[i] * fe_0 + 2.0 * ta1_z_yyz_yyy_1[i] + ta2_zz_yyz_yyy_0[i] * pa_z[i] - ta2_zz_yyz_yyy_1[i] * pc_z[i];

        ta2_zz_yyzz_yyz_0[i] = ta2_zz_zz_yyz_0[i] * fe_0 - ta2_zz_zz_yyz_1[i] * fe_0 + 2.0 * ta2_zz_yzz_yz_0[i] * fe_0 - 2.0 * ta2_zz_yzz_yz_1[i] * fe_0 + ta2_zz_yzz_yyz_0[i] * pa_y[i] - ta2_zz_yzz_yyz_1[i] * pc_y[i];

        ta2_zz_yyzz_yzz_0[i] = ta2_zz_zz_yzz_0[i] * fe_0 - ta2_zz_zz_yzz_1[i] * fe_0 + ta2_zz_yzz_zz_0[i] * fe_0 - ta2_zz_yzz_zz_1[i] * fe_0 + ta2_zz_yzz_yzz_0[i] * pa_y[i] - ta2_zz_yzz_yzz_1[i] * pc_y[i];

        ta2_zz_yyzz_zzz_0[i] = ta2_zz_zz_zzz_0[i] * fe_0 - ta2_zz_zz_zzz_1[i] * fe_0 + ta2_zz_yzz_zzz_0[i] * pa_y[i] - ta2_zz_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 880-890 components of targeted buffer : GF

    auto ta2_zz_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 880);

    auto ta2_zz_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 881);

    auto ta2_zz_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 882);

    auto ta2_zz_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 883);

    auto ta2_zz_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 884);

    auto ta2_zz_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 885);

    auto ta2_zz_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 886);

    auto ta2_zz_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 887);

    auto ta2_zz_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 888);

    auto ta2_zz_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 889);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_yzzz_xxx_0, ta2_zz_yzzz_xxy_0, ta2_zz_yzzz_xxz_0, ta2_zz_yzzz_xyy_0, ta2_zz_yzzz_xyz_0, ta2_zz_yzzz_xzz_0, ta2_zz_yzzz_yyy_0, ta2_zz_yzzz_yyz_0, ta2_zz_yzzz_yzz_0, ta2_zz_yzzz_zzz_0, ta2_zz_zzz_xx_0, ta2_zz_zzz_xx_1, ta2_zz_zzz_xxx_0, ta2_zz_zzz_xxx_1, ta2_zz_zzz_xxy_0, ta2_zz_zzz_xxy_1, ta2_zz_zzz_xxz_0, ta2_zz_zzz_xxz_1, ta2_zz_zzz_xy_0, ta2_zz_zzz_xy_1, ta2_zz_zzz_xyy_0, ta2_zz_zzz_xyy_1, ta2_zz_zzz_xyz_0, ta2_zz_zzz_xyz_1, ta2_zz_zzz_xz_0, ta2_zz_zzz_xz_1, ta2_zz_zzz_xzz_0, ta2_zz_zzz_xzz_1, ta2_zz_zzz_yy_0, ta2_zz_zzz_yy_1, ta2_zz_zzz_yyy_0, ta2_zz_zzz_yyy_1, ta2_zz_zzz_yyz_0, ta2_zz_zzz_yyz_1, ta2_zz_zzz_yz_0, ta2_zz_zzz_yz_1, ta2_zz_zzz_yzz_0, ta2_zz_zzz_yzz_1, ta2_zz_zzz_zz_0, ta2_zz_zzz_zz_1, ta2_zz_zzz_zzz_0, ta2_zz_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzzz_xxx_0[i] = ta2_zz_zzz_xxx_0[i] * pa_y[i] - ta2_zz_zzz_xxx_1[i] * pc_y[i];

        ta2_zz_yzzz_xxy_0[i] = ta2_zz_zzz_xx_0[i] * fe_0 - ta2_zz_zzz_xx_1[i] * fe_0 + ta2_zz_zzz_xxy_0[i] * pa_y[i] - ta2_zz_zzz_xxy_1[i] * pc_y[i];

        ta2_zz_yzzz_xxz_0[i] = ta2_zz_zzz_xxz_0[i] * pa_y[i] - ta2_zz_zzz_xxz_1[i] * pc_y[i];

        ta2_zz_yzzz_xyy_0[i] = 2.0 * ta2_zz_zzz_xy_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xy_1[i] * fe_0 + ta2_zz_zzz_xyy_0[i] * pa_y[i] - ta2_zz_zzz_xyy_1[i] * pc_y[i];

        ta2_zz_yzzz_xyz_0[i] = ta2_zz_zzz_xz_0[i] * fe_0 - ta2_zz_zzz_xz_1[i] * fe_0 + ta2_zz_zzz_xyz_0[i] * pa_y[i] - ta2_zz_zzz_xyz_1[i] * pc_y[i];

        ta2_zz_yzzz_xzz_0[i] = ta2_zz_zzz_xzz_0[i] * pa_y[i] - ta2_zz_zzz_xzz_1[i] * pc_y[i];

        ta2_zz_yzzz_yyy_0[i] = 3.0 * ta2_zz_zzz_yy_0[i] * fe_0 - 3.0 * ta2_zz_zzz_yy_1[i] * fe_0 + ta2_zz_zzz_yyy_0[i] * pa_y[i] - ta2_zz_zzz_yyy_1[i] * pc_y[i];

        ta2_zz_yzzz_yyz_0[i] = 2.0 * ta2_zz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_yz_1[i] * fe_0 + ta2_zz_zzz_yyz_0[i] * pa_y[i] - ta2_zz_zzz_yyz_1[i] * pc_y[i];

        ta2_zz_yzzz_yzz_0[i] = ta2_zz_zzz_zz_0[i] * fe_0 - ta2_zz_zzz_zz_1[i] * fe_0 + ta2_zz_zzz_yzz_0[i] * pa_y[i] - ta2_zz_zzz_yzz_1[i] * pc_y[i];

        ta2_zz_yzzz_zzz_0[i] = ta2_zz_zzz_zzz_0[i] * pa_y[i] - ta2_zz_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 890-900 components of targeted buffer : GF

    auto ta2_zz_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_gf + 890);

    auto ta2_zz_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 891);

    auto ta2_zz_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 892);

    auto ta2_zz_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 893);

    auto ta2_zz_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 894);

    auto ta2_zz_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 895);

    auto ta2_zz_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_gf + 896);

    auto ta2_zz_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 897);

    auto ta2_zz_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 898);

    auto ta2_zz_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_gf + 899);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zzz_xxx_1, ta1_z_zzz_xxy_1, ta1_z_zzz_xxz_1, ta1_z_zzz_xyy_1, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_1, ta1_z_zzz_yyy_1, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_1, ta2_zz_zz_xxx_0, ta2_zz_zz_xxx_1, ta2_zz_zz_xxy_0, ta2_zz_zz_xxy_1, ta2_zz_zz_xxz_0, ta2_zz_zz_xxz_1, ta2_zz_zz_xyy_0, ta2_zz_zz_xyy_1, ta2_zz_zz_xyz_0, ta2_zz_zz_xyz_1, ta2_zz_zz_xzz_0, ta2_zz_zz_xzz_1, ta2_zz_zz_yyy_0, ta2_zz_zz_yyy_1, ta2_zz_zz_yyz_0, ta2_zz_zz_yyz_1, ta2_zz_zz_yzz_0, ta2_zz_zz_yzz_1, ta2_zz_zz_zzz_0, ta2_zz_zz_zzz_1, ta2_zz_zzz_xx_0, ta2_zz_zzz_xx_1, ta2_zz_zzz_xxx_0, ta2_zz_zzz_xxx_1, ta2_zz_zzz_xxy_0, ta2_zz_zzz_xxy_1, ta2_zz_zzz_xxz_0, ta2_zz_zzz_xxz_1, ta2_zz_zzz_xy_0, ta2_zz_zzz_xy_1, ta2_zz_zzz_xyy_0, ta2_zz_zzz_xyy_1, ta2_zz_zzz_xyz_0, ta2_zz_zzz_xyz_1, ta2_zz_zzz_xz_0, ta2_zz_zzz_xz_1, ta2_zz_zzz_xzz_0, ta2_zz_zzz_xzz_1, ta2_zz_zzz_yy_0, ta2_zz_zzz_yy_1, ta2_zz_zzz_yyy_0, ta2_zz_zzz_yyy_1, ta2_zz_zzz_yyz_0, ta2_zz_zzz_yyz_1, ta2_zz_zzz_yz_0, ta2_zz_zzz_yz_1, ta2_zz_zzz_yzz_0, ta2_zz_zzz_yzz_1, ta2_zz_zzz_zz_0, ta2_zz_zzz_zz_1, ta2_zz_zzz_zzz_0, ta2_zz_zzz_zzz_1, ta2_zz_zzzz_xxx_0, ta2_zz_zzzz_xxy_0, ta2_zz_zzzz_xxz_0, ta2_zz_zzzz_xyy_0, ta2_zz_zzzz_xyz_0, ta2_zz_zzzz_xzz_0, ta2_zz_zzzz_yyy_0, ta2_zz_zzzz_yyz_0, ta2_zz_zzzz_yzz_0, ta2_zz_zzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzzz_xxx_0[i] = 3.0 * ta2_zz_zz_xxx_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxx_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxx_1[i] + ta2_zz_zzz_xxx_0[i] * pa_z[i] - ta2_zz_zzz_xxx_1[i] * pc_z[i];

        ta2_zz_zzzz_xxy_0[i] = 3.0 * ta2_zz_zz_xxy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxy_1[i] + ta2_zz_zzz_xxy_0[i] * pa_z[i] - ta2_zz_zzz_xxy_1[i] * pc_z[i];

        ta2_zz_zzzz_xxz_0[i] = 3.0 * ta2_zz_zz_xxz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxz_1[i] * fe_0 + ta2_zz_zzz_xx_0[i] * fe_0 - ta2_zz_zzz_xx_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxz_1[i] + ta2_zz_zzz_xxz_0[i] * pa_z[i] - ta2_zz_zzz_xxz_1[i] * pc_z[i];

        ta2_zz_zzzz_xyy_0[i] = 3.0 * ta2_zz_zz_xyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyy_1[i] + ta2_zz_zzz_xyy_0[i] * pa_z[i] - ta2_zz_zzz_xyy_1[i] * pc_z[i];

        ta2_zz_zzzz_xyz_0[i] = 3.0 * ta2_zz_zz_xyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyz_1[i] * fe_0 + ta2_zz_zzz_xy_0[i] * fe_0 - ta2_zz_zzz_xy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyz_1[i] + ta2_zz_zzz_xyz_0[i] * pa_z[i] - ta2_zz_zzz_xyz_1[i] * pc_z[i];

        ta2_zz_zzzz_xzz_0[i] = 3.0 * ta2_zz_zz_xzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xzz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_xz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xzz_1[i] + ta2_zz_zzz_xzz_0[i] * pa_z[i] - ta2_zz_zzz_xzz_1[i] * pc_z[i];

        ta2_zz_zzzz_yyy_0[i] = 3.0 * ta2_zz_zz_yyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyy_1[i] + ta2_zz_zzz_yyy_0[i] * pa_z[i] - ta2_zz_zzz_yyy_1[i] * pc_z[i];

        ta2_zz_zzzz_yyz_0[i] = 3.0 * ta2_zz_zz_yyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyz_1[i] * fe_0 + ta2_zz_zzz_yy_0[i] * fe_0 - ta2_zz_zzz_yy_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyz_1[i] + ta2_zz_zzz_yyz_0[i] * pa_z[i] - ta2_zz_zzz_yyz_1[i] * pc_z[i];

        ta2_zz_zzzz_yzz_0[i] = 3.0 * ta2_zz_zz_yzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yzz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_yz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_yz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yzz_1[i] + ta2_zz_zzz_yzz_0[i] * pa_z[i] - ta2_zz_zzz_yzz_1[i] * pc_z[i];

        ta2_zz_zzzz_zzz_0[i] = 3.0 * ta2_zz_zz_zzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_zzz_1[i] * fe_0 + 3.0 * ta2_zz_zzz_zz_0[i] * fe_0 - 3.0 * ta2_zz_zzz_zz_1[i] * fe_0 + 2.0 * ta1_z_zzz_zzz_1[i] + ta2_zz_zzz_zzz_0[i] * pa_z[i] - ta2_zz_zzz_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

