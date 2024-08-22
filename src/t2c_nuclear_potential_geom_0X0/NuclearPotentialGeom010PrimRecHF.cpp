#include "NuclearPotentialGeom010PrimRecHF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_hf(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_hf,
                                        const size_t idx_npot_geom_010_0_ff,
                                        const size_t idx_npot_geom_010_1_ff,
                                        const size_t idx_npot_geom_010_0_gd,
                                        const size_t idx_npot_geom_010_1_gd,
                                        const size_t idx_npot_1_gf,
                                        const size_t idx_npot_geom_010_0_gf,
                                        const size_t idx_npot_geom_010_1_gf,
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

    auto ta1_x_xxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 19);

    auto ta1_x_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 20);

    auto ta1_x_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 21);

    auto ta1_x_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 22);

    auto ta1_x_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 23);

    auto ta1_x_xxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 24);

    auto ta1_x_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 25);

    auto ta1_x_xxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 26);

    auto ta1_x_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 29);

    auto ta1_x_xyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 30);

    auto ta1_x_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 31);

    auto ta1_x_xyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 32);

    auto ta1_x_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 33);

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

    auto ta1_x_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 79);

    auto ta1_x_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 80);

    auto ta1_x_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 82);

    auto ta1_x_yzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 84);

    auto ta1_x_yzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 85);

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

    auto ta1_y_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 115);

    auto ta1_y_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 116);

    auto ta1_y_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 117);

    auto ta1_y_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 118);

    auto ta1_y_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 120);

    auto ta1_y_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 121);

    auto ta1_y_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 123);

    auto ta1_y_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 127);

    auto ta1_y_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 128);

    auto ta1_y_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 129);

    auto ta1_y_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 131);

    auto ta1_y_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 133);

    auto ta1_y_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 134);

    auto ta1_y_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 136);

    auto ta1_y_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 137);

    auto ta1_y_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 138);

    auto ta1_y_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 139);

    auto ta1_y_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 147);

    auto ta1_y_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 148);

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

    auto ta1_y_yyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 173);

    auto ta1_y_yyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 174);

    auto ta1_y_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 176);

    auto ta1_y_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 177);

    auto ta1_y_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 178);

    auto ta1_y_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 179);

    auto ta1_y_yzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 181);

    auto ta1_y_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 182);

    auto ta1_y_yzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 183);

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

    auto ta1_z_xxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 212);

    auto ta1_z_xxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 215);

    auto ta1_z_xxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 216);

    auto ta1_z_xxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 217);

    auto ta1_z_xxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 218);

    auto ta1_z_xxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 220);

    auto ta1_z_xxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 221);

    auto ta1_z_xxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 222);

    auto ta1_z_xxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 223);

    auto ta1_z_xxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 225);

    auto ta1_z_xxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 227);

    auto ta1_z_xxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 228);

    auto ta1_z_xxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 229);

    auto ta1_z_xyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 231);

    auto ta1_z_xyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 233);

    auto ta1_z_xyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 234);

    auto ta1_z_xyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 236);

    auto ta1_z_xyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 237);

    auto ta1_z_xyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 238);

    auto ta1_z_xyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 239);

    auto ta1_z_xyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 247);

    auto ta1_z_xyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 248);

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

    auto ta1_z_yyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 275);

    auto ta1_z_yyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_ff + 276);

    auto ta1_z_yyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 277);

    auto ta1_z_yyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 278);

    auto ta1_z_yyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 279);

    auto ta1_z_yzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_ff + 280);

    auto ta1_z_yzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_ff + 282);

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

    auto ta1_x_xxy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 19);

    auto ta1_x_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 20);

    auto ta1_x_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 21);

    auto ta1_x_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 22);

    auto ta1_x_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 23);

    auto ta1_x_xxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 24);

    auto ta1_x_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 25);

    auto ta1_x_xxz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 26);

    auto ta1_x_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 29);

    auto ta1_x_xyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 30);

    auto ta1_x_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 31);

    auto ta1_x_xyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 32);

    auto ta1_x_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 33);

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

    auto ta1_x_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 79);

    auto ta1_x_yzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 80);

    auto ta1_x_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 82);

    auto ta1_x_yzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 84);

    auto ta1_x_yzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 85);

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

    auto ta1_y_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 123);

    auto ta1_y_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 127);

    auto ta1_y_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 128);

    auto ta1_y_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 129);

    auto ta1_y_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 131);

    auto ta1_y_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 133);

    auto ta1_y_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 134);

    auto ta1_y_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 136);

    auto ta1_y_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 137);

    auto ta1_y_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 138);

    auto ta1_y_xyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 139);

    auto ta1_y_xyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 147);

    auto ta1_y_xyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 148);

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

    auto ta1_y_yyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 173);

    auto ta1_y_yyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 174);

    auto ta1_y_yyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 176);

    auto ta1_y_yyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 177);

    auto ta1_y_yyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 178);

    auto ta1_y_yyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 179);

    auto ta1_y_yzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 181);

    auto ta1_y_yzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 182);

    auto ta1_y_yzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 183);

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

    auto ta1_z_xxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 212);

    auto ta1_z_xxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 215);

    auto ta1_z_xxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 216);

    auto ta1_z_xxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 217);

    auto ta1_z_xxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 218);

    auto ta1_z_xxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_ff + 220);

    auto ta1_z_xxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 221);

    auto ta1_z_xxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 222);

    auto ta1_z_xxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 223);

    auto ta1_z_xxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 225);

    auto ta1_z_xxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 227);

    auto ta1_z_xxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 228);

    auto ta1_z_xxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 229);

    auto ta1_z_xyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 231);

    auto ta1_z_xyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 233);

    auto ta1_z_xyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 234);

    auto ta1_z_xyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_ff + 236);

    auto ta1_z_xyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 237);

    auto ta1_z_xyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 238);

    auto ta1_z_xyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 239);

    auto ta1_z_xyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 247);

    auto ta1_z_xyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_ff + 248);

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

    // Set up components of auxiliary buffer : GD

    auto ta1_x_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd);

    auto ta1_x_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 1);

    auto ta1_x_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 2);

    auto ta1_x_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 3);

    auto ta1_x_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 4);

    auto ta1_x_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 5);

    auto ta1_x_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 6);

    auto ta1_x_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 7);

    auto ta1_x_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 8);

    auto ta1_x_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 12);

    auto ta1_x_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 13);

    auto ta1_x_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 14);

    auto ta1_x_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 16);

    auto ta1_x_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 17);

    auto ta1_x_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 18);

    auto ta1_x_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 19);

    auto ta1_x_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 20);

    auto ta1_x_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 21);

    auto ta1_x_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 22);

    auto ta1_x_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 30);

    auto ta1_x_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 31);

    auto ta1_x_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 32);

    auto ta1_x_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 33);

    auto ta1_x_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 34);

    auto ta1_x_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 35);

    auto ta1_x_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 37);

    auto ta1_x_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 54);

    auto ta1_x_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 55);

    auto ta1_x_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 56);

    auto ta1_x_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 60);

    auto ta1_x_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 61);

    auto ta1_x_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 62);

    auto ta1_x_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 63);

    auto ta1_x_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 64);

    auto ta1_x_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 65);

    auto ta1_x_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 74);

    auto ta1_x_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 76);

    auto ta1_x_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 77);

    auto ta1_x_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 80);

    auto ta1_x_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 82);

    auto ta1_x_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 83);

    auto ta1_x_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 84);

    auto ta1_x_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 85);

    auto ta1_x_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 86);

    auto ta1_x_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 87);

    auto ta1_x_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 88);

    auto ta1_x_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 89);

    auto ta1_y_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 90);

    auto ta1_y_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 91);

    auto ta1_y_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 92);

    auto ta1_y_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 93);

    auto ta1_y_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 94);

    auto ta1_y_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 95);

    auto ta1_y_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 97);

    auto ta1_y_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 108);

    auto ta1_y_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 109);

    auto ta1_y_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 110);

    auto ta1_y_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 111);

    auto ta1_y_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 112);

    auto ta1_y_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 122);

    auto ta1_y_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 124);

    auto ta1_y_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 125);

    auto ta1_y_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 127);

    auto ta1_y_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 129);

    auto ta1_y_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 130);

    auto ta1_y_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 146);

    auto ta1_y_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 148);

    auto ta1_y_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 149);

    auto ta1_y_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 150);

    auto ta1_y_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 151);

    auto ta1_y_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 152);

    auto ta1_y_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 153);

    auto ta1_y_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 154);

    auto ta1_y_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 155);

    auto ta1_y_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 157);

    auto ta1_y_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 158);

    auto ta1_y_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 159);

    auto ta1_y_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 160);

    auto ta1_y_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 161);

    auto ta1_y_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 162);

    auto ta1_y_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 163);

    auto ta1_y_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 164);

    auto ta1_y_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 165);

    auto ta1_y_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 166);

    auto ta1_y_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 167);

    auto ta1_y_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 169);

    auto ta1_y_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 171);

    auto ta1_y_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 172);

    auto ta1_y_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 174);

    auto ta1_y_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 175);

    auto ta1_y_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 176);

    auto ta1_y_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 177);

    auto ta1_y_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 178);

    auto ta1_y_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 179);

    auto ta1_z_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 180);

    auto ta1_z_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 181);

    auto ta1_z_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 182);

    auto ta1_z_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 183);

    auto ta1_z_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 184);

    auto ta1_z_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 185);

    auto ta1_z_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 194);

    auto ta1_z_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 199);

    auto ta1_z_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 201);

    auto ta1_z_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 202);

    auto ta1_z_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 210);

    auto ta1_z_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 211);

    auto ta1_z_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 212);

    auto ta1_z_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 214);

    auto ta1_z_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 215);

    auto ta1_z_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 217);

    auto ta1_z_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 219);

    auto ta1_z_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 220);

    auto ta1_z_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 236);

    auto ta1_z_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 238);

    auto ta1_z_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 239);

    auto ta1_z_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 240);

    auto ta1_z_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 241);

    auto ta1_z_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 242);

    auto ta1_z_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 243);

    auto ta1_z_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 244);

    auto ta1_z_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 245);

    auto ta1_z_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 248);

    auto ta1_z_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 250);

    auto ta1_z_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 251);

    auto ta1_z_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 252);

    auto ta1_z_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 253);

    auto ta1_z_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 254);

    auto ta1_z_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 255);

    auto ta1_z_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 256);

    auto ta1_z_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 257);

    auto ta1_z_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 259);

    auto ta1_z_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 260);

    auto ta1_z_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 261);

    auto ta1_z_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 262);

    auto ta1_z_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 263);

    auto ta1_z_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 264);

    auto ta1_z_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 265);

    auto ta1_z_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 266);

    auto ta1_z_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 267);

    auto ta1_z_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 268);

    auto ta1_z_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 269);

    // Set up components of auxiliary buffer : GD

    auto ta1_x_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd);

    auto ta1_x_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 1);

    auto ta1_x_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 2);

    auto ta1_x_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 3);

    auto ta1_x_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 4);

    auto ta1_x_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 5);

    auto ta1_x_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 6);

    auto ta1_x_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 7);

    auto ta1_x_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 8);

    auto ta1_x_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 12);

    auto ta1_x_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 13);

    auto ta1_x_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 14);

    auto ta1_x_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 16);

    auto ta1_x_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 17);

    auto ta1_x_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 18);

    auto ta1_x_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 19);

    auto ta1_x_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 20);

    auto ta1_x_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 21);

    auto ta1_x_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 22);

    auto ta1_x_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 30);

    auto ta1_x_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 31);

    auto ta1_x_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 32);

    auto ta1_x_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 33);

    auto ta1_x_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 34);

    auto ta1_x_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 35);

    auto ta1_x_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 37);

    auto ta1_x_xzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 54);

    auto ta1_x_xzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 55);

    auto ta1_x_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 56);

    auto ta1_x_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 60);

    auto ta1_x_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 61);

    auto ta1_x_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 62);

    auto ta1_x_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 63);

    auto ta1_x_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 64);

    auto ta1_x_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 65);

    auto ta1_x_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 74);

    auto ta1_x_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 76);

    auto ta1_x_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 77);

    auto ta1_x_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 80);

    auto ta1_x_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 82);

    auto ta1_x_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 83);

    auto ta1_x_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 84);

    auto ta1_x_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 85);

    auto ta1_x_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 86);

    auto ta1_x_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 87);

    auto ta1_x_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 88);

    auto ta1_x_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 89);

    auto ta1_y_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 90);

    auto ta1_y_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 91);

    auto ta1_y_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 92);

    auto ta1_y_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 93);

    auto ta1_y_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 94);

    auto ta1_y_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 95);

    auto ta1_y_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 97);

    auto ta1_y_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 108);

    auto ta1_y_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 109);

    auto ta1_y_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 110);

    auto ta1_y_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 111);

    auto ta1_y_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 112);

    auto ta1_y_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 122);

    auto ta1_y_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 124);

    auto ta1_y_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 125);

    auto ta1_y_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 127);

    auto ta1_y_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 129);

    auto ta1_y_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 130);

    auto ta1_y_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 146);

    auto ta1_y_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 148);

    auto ta1_y_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 149);

    auto ta1_y_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 150);

    auto ta1_y_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 151);

    auto ta1_y_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 152);

    auto ta1_y_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 153);

    auto ta1_y_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 154);

    auto ta1_y_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 155);

    auto ta1_y_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 157);

    auto ta1_y_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 158);

    auto ta1_y_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 159);

    auto ta1_y_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 160);

    auto ta1_y_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 161);

    auto ta1_y_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 162);

    auto ta1_y_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 163);

    auto ta1_y_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 164);

    auto ta1_y_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 165);

    auto ta1_y_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 166);

    auto ta1_y_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 167);

    auto ta1_y_yzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 169);

    auto ta1_y_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 171);

    auto ta1_y_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 172);

    auto ta1_y_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 174);

    auto ta1_y_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 175);

    auto ta1_y_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 176);

    auto ta1_y_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 177);

    auto ta1_y_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 178);

    auto ta1_y_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 179);

    auto ta1_z_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 180);

    auto ta1_z_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 181);

    auto ta1_z_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 182);

    auto ta1_z_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 183);

    auto ta1_z_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 184);

    auto ta1_z_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 185);

    auto ta1_z_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 194);

    auto ta1_z_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 199);

    auto ta1_z_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 201);

    auto ta1_z_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 202);

    auto ta1_z_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 210);

    auto ta1_z_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 211);

    auto ta1_z_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 212);

    auto ta1_z_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 214);

    auto ta1_z_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 215);

    auto ta1_z_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 217);

    auto ta1_z_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 219);

    auto ta1_z_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 220);

    auto ta1_z_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 236);

    auto ta1_z_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 238);

    auto ta1_z_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 239);

    auto ta1_z_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 240);

    auto ta1_z_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 241);

    auto ta1_z_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 242);

    auto ta1_z_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 243);

    auto ta1_z_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 244);

    auto ta1_z_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 245);

    auto ta1_z_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 248);

    auto ta1_z_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 250);

    auto ta1_z_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 251);

    auto ta1_z_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 252);

    auto ta1_z_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 253);

    auto ta1_z_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 254);

    auto ta1_z_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 255);

    auto ta1_z_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 256);

    auto ta1_z_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 257);

    auto ta1_z_yzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 259);

    auto ta1_z_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 260);

    auto ta1_z_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 261);

    auto ta1_z_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 262);

    auto ta1_z_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 263);

    auto ta1_z_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 264);

    auto ta1_z_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 265);

    auto ta1_z_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 266);

    auto ta1_z_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 267);

    auto ta1_z_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 268);

    auto ta1_z_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 269);

    // Set up components of auxiliary buffer : GF

    auto ta_xxxx_xxx_1 = pbuffer.data(idx_npot_1_gf);

    auto ta_xxxx_xxy_1 = pbuffer.data(idx_npot_1_gf + 1);

    auto ta_xxxx_xxz_1 = pbuffer.data(idx_npot_1_gf + 2);

    auto ta_xxxx_xyy_1 = pbuffer.data(idx_npot_1_gf + 3);

    auto ta_xxxx_xyz_1 = pbuffer.data(idx_npot_1_gf + 4);

    auto ta_xxxx_xzz_1 = pbuffer.data(idx_npot_1_gf + 5);

    auto ta_xxxx_yyy_1 = pbuffer.data(idx_npot_1_gf + 6);

    auto ta_xxxx_yyz_1 = pbuffer.data(idx_npot_1_gf + 7);

    auto ta_xxxx_yzz_1 = pbuffer.data(idx_npot_1_gf + 8);

    auto ta_xxxx_zzz_1 = pbuffer.data(idx_npot_1_gf + 9);

    auto ta_xxxy_xxx_1 = pbuffer.data(idx_npot_1_gf + 10);

    auto ta_xxxy_xxy_1 = pbuffer.data(idx_npot_1_gf + 11);

    auto ta_xxxy_xxz_1 = pbuffer.data(idx_npot_1_gf + 12);

    auto ta_xxxy_xyy_1 = pbuffer.data(idx_npot_1_gf + 13);

    auto ta_xxxy_xzz_1 = pbuffer.data(idx_npot_1_gf + 15);

    auto ta_xxxy_yyy_1 = pbuffer.data(idx_npot_1_gf + 16);

    auto ta_xxxz_xxx_1 = pbuffer.data(idx_npot_1_gf + 20);

    auto ta_xxxz_xxy_1 = pbuffer.data(idx_npot_1_gf + 21);

    auto ta_xxxz_xxz_1 = pbuffer.data(idx_npot_1_gf + 22);

    auto ta_xxxz_xyy_1 = pbuffer.data(idx_npot_1_gf + 23);

    auto ta_xxxz_xzz_1 = pbuffer.data(idx_npot_1_gf + 25);

    auto ta_xxxz_zzz_1 = pbuffer.data(idx_npot_1_gf + 29);

    auto ta_xxyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 30);

    auto ta_xxyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 31);

    auto ta_xxyy_xxz_1 = pbuffer.data(idx_npot_1_gf + 32);

    auto ta_xxyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 33);

    auto ta_xxyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 34);

    auto ta_xxyy_xzz_1 = pbuffer.data(idx_npot_1_gf + 35);

    auto ta_xxyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 36);

    auto ta_xxyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 37);

    auto ta_xxyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 38);

    auto ta_xxzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 50);

    auto ta_xxzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 51);

    auto ta_xxzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 52);

    auto ta_xxzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 53);

    auto ta_xxzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 54);

    auto ta_xxzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 55);

    auto ta_xxzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 57);

    auto ta_xxzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 58);

    auto ta_xxzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 59);

    auto ta_xyyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 60);

    auto ta_xyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 61);

    auto ta_xyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 63);

    auto ta_xyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 66);

    auto ta_xyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 67);

    auto ta_xyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 68);

    auto ta_xzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 90);

    auto ta_xzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 92);

    auto ta_xzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 95);

    auto ta_xzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 97);

    auto ta_xzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 98);

    auto ta_xzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 99);

    auto ta_yyyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 100);

    auto ta_yyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 101);

    auto ta_yyyy_xxz_1 = pbuffer.data(idx_npot_1_gf + 102);

    auto ta_yyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 103);

    auto ta_yyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 104);

    auto ta_yyyy_xzz_1 = pbuffer.data(idx_npot_1_gf + 105);

    auto ta_yyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 106);

    auto ta_yyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 107);

    auto ta_yyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 108);

    auto ta_yyyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 109);

    auto ta_yyyz_xxy_1 = pbuffer.data(idx_npot_1_gf + 111);

    auto ta_yyyz_xyy_1 = pbuffer.data(idx_npot_1_gf + 113);

    auto ta_yyyz_yyy_1 = pbuffer.data(idx_npot_1_gf + 116);

    auto ta_yyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 117);

    auto ta_yyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 118);

    auto ta_yyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 119);

    auto ta_yyzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 121);

    auto ta_yyzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 122);

    auto ta_yyzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 123);

    auto ta_yyzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 124);

    auto ta_yyzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 125);

    auto ta_yyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 126);

    auto ta_yyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 127);

    auto ta_yyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 128);

    auto ta_yyzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 129);

    auto ta_yzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 132);

    auto ta_yzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 135);

    auto ta_yzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 136);

    auto ta_yzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 137);

    auto ta_yzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 138);

    auto ta_yzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 139);

    auto ta_zzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 140);

    auto ta_zzzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 141);

    auto ta_zzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 142);

    auto ta_zzzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 143);

    auto ta_zzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 144);

    auto ta_zzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 145);

    auto ta_zzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 146);

    auto ta_zzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 147);

    auto ta_zzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 148);

    auto ta_zzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 149);

    // Set up components of auxiliary buffer : GF

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

    auto ta1_x_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 10);

    auto ta1_x_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 11);

    auto ta1_x_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 12);

    auto ta1_x_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 13);

    auto ta1_x_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 14);

    auto ta1_x_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 15);

    auto ta1_x_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 16);

    auto ta1_x_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 19);

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

    auto ta1_x_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 42);

    auto ta1_x_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 45);

    auto ta1_x_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 49);

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

    auto ta1_x_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 60);

    auto ta1_x_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 61);

    auto ta1_x_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 62);

    auto ta1_x_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 63);

    auto ta1_x_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 64);

    auto ta1_x_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 65);

    auto ta1_x_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 66);

    auto ta1_x_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 67);

    auto ta1_x_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 68);

    auto ta1_x_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 71);

    auto ta1_x_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 72);

    auto ta1_x_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 73);

    auto ta1_x_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 75);

    auto ta1_x_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 80);

    auto ta1_x_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 82);

    auto ta1_x_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 85);

    auto ta1_x_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 90);

    auto ta1_x_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 91);

    auto ta1_x_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 92);

    auto ta1_x_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 93);

    auto ta1_x_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 94);

    auto ta1_x_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 95);

    auto ta1_x_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 97);

    auto ta1_x_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 98);

    auto ta1_x_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 99);

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

    auto ta1_x_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 111);

    auto ta1_x_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 112);

    auto ta1_x_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 113);

    auto ta1_x_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 115);

    auto ta1_x_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 116);

    auto ta1_x_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 117);

    auto ta1_x_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 118);

    auto ta1_x_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 119);

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

    auto ta1_x_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 130);

    auto ta1_x_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 132);

    auto ta1_x_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 134);

    auto ta1_x_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 135);

    auto ta1_x_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 136);

    auto ta1_x_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 137);

    auto ta1_x_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 138);

    auto ta1_x_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 139);

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

    auto ta1_y_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 160);

    auto ta1_y_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 161);

    auto ta1_y_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 162);

    auto ta1_y_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 163);

    auto ta1_y_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 164);

    auto ta1_y_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 165);

    auto ta1_y_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 166);

    auto ta1_y_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 167);

    auto ta1_y_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 168);

    auto ta1_y_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 170);

    auto ta1_y_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 171);

    auto ta1_y_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 172);

    auto ta1_y_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 173);

    auto ta1_y_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 175);

    auto ta1_y_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 177);

    auto ta1_y_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 178);

    auto ta1_y_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 179);

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

    auto ta1_y_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 191);

    auto ta1_y_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 193);

    auto ta1_y_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 197);

    auto ta1_y_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 198);

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

    auto ta1_y_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 210);

    auto ta1_y_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 211);

    auto ta1_y_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 213);

    auto ta1_y_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 214);

    auto ta1_y_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 216);

    auto ta1_y_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 217);

    auto ta1_y_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 218);

    auto ta1_y_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 219);

    auto ta1_y_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 227);

    auto ta1_y_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 228);

    auto ta1_y_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 229);

    auto ta1_y_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 236);

    auto ta1_y_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 237);

    auto ta1_y_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 238);

    auto ta1_y_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 240);

    auto ta1_y_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 242);

    auto ta1_y_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 244);

    auto ta1_y_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 245);

    auto ta1_y_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 246);

    auto ta1_y_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 247);

    auto ta1_y_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 248);

    auto ta1_y_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 249);

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

    auto ta1_y_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 281);

    auto ta1_y_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 282);

    auto ta1_y_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 283);

    auto ta1_y_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 284);

    auto ta1_y_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 285);

    auto ta1_y_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 286);

    auto ta1_y_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 287);

    auto ta1_y_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 288);

    auto ta1_y_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 289);

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

    auto ta1_z_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 310);

    auto ta1_z_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 311);

    auto ta1_z_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 312);

    auto ta1_z_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 313);

    auto ta1_z_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 315);

    auto ta1_z_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 316);

    auto ta1_z_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 317);

    auto ta1_z_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 318);

    auto ta1_z_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 320);

    auto ta1_z_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 321);

    auto ta1_z_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 322);

    auto ta1_z_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 323);

    auto ta1_z_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 324);

    auto ta1_z_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 325);

    auto ta1_z_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 327);

    auto ta1_z_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 328);

    auto ta1_z_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 329);

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

    auto ta1_z_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 342);

    auto ta1_z_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 345);

    auto ta1_z_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 347);

    auto ta1_z_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 348);

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

    auto ta1_z_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 360);

    auto ta1_z_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 361);

    auto ta1_z_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 363);

    auto ta1_z_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 364);

    auto ta1_z_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 366);

    auto ta1_z_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 367);

    auto ta1_z_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 368);

    auto ta1_z_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 369);

    auto ta1_z_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 377);

    auto ta1_z_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 378);

    auto ta1_z_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 379);

    auto ta1_z_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 386);

    auto ta1_z_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 387);

    auto ta1_z_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 388);

    auto ta1_z_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 390);

    auto ta1_z_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 392);

    auto ta1_z_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 394);

    auto ta1_z_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 395);

    auto ta1_z_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 396);

    auto ta1_z_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 397);

    auto ta1_z_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 398);

    auto ta1_z_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 399);

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

    auto ta1_z_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 411);

    auto ta1_z_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 412);

    auto ta1_z_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 413);

    auto ta1_z_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 414);

    auto ta1_z_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 415);

    auto ta1_z_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 416);

    auto ta1_z_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 417);

    auto ta1_z_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 418);

    auto ta1_z_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 419);

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

    // Set up components of auxiliary buffer : GF

    auto ta1_x_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf);

    auto ta1_x_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 1);

    auto ta1_x_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 2);

    auto ta1_x_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 3);

    auto ta1_x_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 4);

    auto ta1_x_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 5);

    auto ta1_x_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 6);

    auto ta1_x_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 7);

    auto ta1_x_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 8);

    auto ta1_x_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 9);

    auto ta1_x_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 10);

    auto ta1_x_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 11);

    auto ta1_x_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 12);

    auto ta1_x_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 13);

    auto ta1_x_xxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 14);

    auto ta1_x_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 15);

    auto ta1_x_xxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 16);

    auto ta1_x_xxxy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 19);

    auto ta1_x_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 20);

    auto ta1_x_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 21);

    auto ta1_x_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 22);

    auto ta1_x_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 23);

    auto ta1_x_xxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 24);

    auto ta1_x_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 25);

    auto ta1_x_xxxz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 26);

    auto ta1_x_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 27);

    auto ta1_x_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 28);

    auto ta1_x_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 29);

    auto ta1_x_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 30);

    auto ta1_x_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 31);

    auto ta1_x_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 32);

    auto ta1_x_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 33);

    auto ta1_x_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 34);

    auto ta1_x_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 35);

    auto ta1_x_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 36);

    auto ta1_x_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 37);

    auto ta1_x_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 38);

    auto ta1_x_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 39);

    auto ta1_x_xxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 42);

    auto ta1_x_xxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 45);

    auto ta1_x_xxyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 49);

    auto ta1_x_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 50);

    auto ta1_x_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 51);

    auto ta1_x_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 52);

    auto ta1_x_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 53);

    auto ta1_x_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 54);

    auto ta1_x_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 55);

    auto ta1_x_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 56);

    auto ta1_x_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 57);

    auto ta1_x_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 58);

    auto ta1_x_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 59);

    auto ta1_x_xyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 60);

    auto ta1_x_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 61);

    auto ta1_x_xyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 62);

    auto ta1_x_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 63);

    auto ta1_x_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 64);

    auto ta1_x_xyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 65);

    auto ta1_x_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 66);

    auto ta1_x_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 67);

    auto ta1_x_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 68);

    auto ta1_x_xyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 71);

    auto ta1_x_xyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 72);

    auto ta1_x_xyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 73);

    auto ta1_x_xyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 75);

    auto ta1_x_xyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 80);

    auto ta1_x_xyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 82);

    auto ta1_x_xyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 85);

    auto ta1_x_xzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 90);

    auto ta1_x_xzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 91);

    auto ta1_x_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 92);

    auto ta1_x_xzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 93);

    auto ta1_x_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 94);

    auto ta1_x_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 95);

    auto ta1_x_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 97);

    auto ta1_x_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 98);

    auto ta1_x_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 99);

    auto ta1_x_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 100);

    auto ta1_x_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 101);

    auto ta1_x_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 102);

    auto ta1_x_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 103);

    auto ta1_x_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 104);

    auto ta1_x_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 105);

    auto ta1_x_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 106);

    auto ta1_x_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 107);

    auto ta1_x_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 108);

    auto ta1_x_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 109);

    auto ta1_x_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 111);

    auto ta1_x_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 112);

    auto ta1_x_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 113);

    auto ta1_x_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 115);

    auto ta1_x_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 116);

    auto ta1_x_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 117);

    auto ta1_x_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 118);

    auto ta1_x_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 119);

    auto ta1_x_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 120);

    auto ta1_x_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 121);

    auto ta1_x_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 122);

    auto ta1_x_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 123);

    auto ta1_x_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 124);

    auto ta1_x_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 125);

    auto ta1_x_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 126);

    auto ta1_x_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 127);

    auto ta1_x_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 128);

    auto ta1_x_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 129);

    auto ta1_x_yzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 130);

    auto ta1_x_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 132);

    auto ta1_x_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 134);

    auto ta1_x_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 135);

    auto ta1_x_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 136);

    auto ta1_x_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 137);

    auto ta1_x_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 138);

    auto ta1_x_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 139);

    auto ta1_x_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 140);

    auto ta1_x_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 141);

    auto ta1_x_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 142);

    auto ta1_x_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 143);

    auto ta1_x_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 144);

    auto ta1_x_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 145);

    auto ta1_x_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 146);

    auto ta1_x_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 147);

    auto ta1_x_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 148);

    auto ta1_x_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 149);

    auto ta1_y_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 150);

    auto ta1_y_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 151);

    auto ta1_y_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 152);

    auto ta1_y_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 153);

    auto ta1_y_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 154);

    auto ta1_y_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 155);

    auto ta1_y_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 156);

    auto ta1_y_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 157);

    auto ta1_y_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 158);

    auto ta1_y_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 159);

    auto ta1_y_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 160);

    auto ta1_y_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 161);

    auto ta1_y_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 162);

    auto ta1_y_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 163);

    auto ta1_y_xxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 164);

    auto ta1_y_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 165);

    auto ta1_y_xxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 166);

    auto ta1_y_xxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 167);

    auto ta1_y_xxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 168);

    auto ta1_y_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 170);

    auto ta1_y_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 171);

    auto ta1_y_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 172);

    auto ta1_y_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 173);

    auto ta1_y_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 175);

    auto ta1_y_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 177);

    auto ta1_y_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 178);

    auto ta1_y_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 179);

    auto ta1_y_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 180);

    auto ta1_y_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 181);

    auto ta1_y_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 182);

    auto ta1_y_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 183);

    auto ta1_y_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 184);

    auto ta1_y_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 185);

    auto ta1_y_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 186);

    auto ta1_y_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 187);

    auto ta1_y_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 188);

    auto ta1_y_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 189);

    auto ta1_y_xxyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 191);

    auto ta1_y_xxyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 193);

    auto ta1_y_xxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 197);

    auto ta1_y_xxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 198);

    auto ta1_y_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 200);

    auto ta1_y_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 201);

    auto ta1_y_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 202);

    auto ta1_y_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 203);

    auto ta1_y_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 204);

    auto ta1_y_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 205);

    auto ta1_y_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 206);

    auto ta1_y_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 207);

    auto ta1_y_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 208);

    auto ta1_y_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 209);

    auto ta1_y_xyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 210);

    auto ta1_y_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 211);

    auto ta1_y_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 213);

    auto ta1_y_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 214);

    auto ta1_y_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 216);

    auto ta1_y_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 217);

    auto ta1_y_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 218);

    auto ta1_y_xyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 219);

    auto ta1_y_xyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 227);

    auto ta1_y_xyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 228);

    auto ta1_y_xyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 229);

    auto ta1_y_xyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 236);

    auto ta1_y_xyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 237);

    auto ta1_y_xyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 238);

    auto ta1_y_xzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 240);

    auto ta1_y_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 242);

    auto ta1_y_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 244);

    auto ta1_y_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 245);

    auto ta1_y_xzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 246);

    auto ta1_y_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 247);

    auto ta1_y_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 248);

    auto ta1_y_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 249);

    auto ta1_y_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 250);

    auto ta1_y_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 251);

    auto ta1_y_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 252);

    auto ta1_y_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 253);

    auto ta1_y_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 254);

    auto ta1_y_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 255);

    auto ta1_y_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 256);

    auto ta1_y_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 257);

    auto ta1_y_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 258);

    auto ta1_y_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 259);

    auto ta1_y_yyyz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 260);

    auto ta1_y_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 261);

    auto ta1_y_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 262);

    auto ta1_y_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 263);

    auto ta1_y_yyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 264);

    auto ta1_y_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 265);

    auto ta1_y_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 266);

    auto ta1_y_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 267);

    auto ta1_y_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 268);

    auto ta1_y_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 269);

    auto ta1_y_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 270);

    auto ta1_y_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 271);

    auto ta1_y_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 272);

    auto ta1_y_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 273);

    auto ta1_y_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 274);

    auto ta1_y_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 275);

    auto ta1_y_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 276);

    auto ta1_y_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 277);

    auto ta1_y_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 278);

    auto ta1_y_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 279);

    auto ta1_y_yzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 281);

    auto ta1_y_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 282);

    auto ta1_y_yzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 283);

    auto ta1_y_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 284);

    auto ta1_y_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 285);

    auto ta1_y_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 286);

    auto ta1_y_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 287);

    auto ta1_y_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 288);

    auto ta1_y_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 289);

    auto ta1_y_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 290);

    auto ta1_y_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 291);

    auto ta1_y_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 292);

    auto ta1_y_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 293);

    auto ta1_y_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 294);

    auto ta1_y_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 295);

    auto ta1_y_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 296);

    auto ta1_y_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 297);

    auto ta1_y_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 298);

    auto ta1_y_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 299);

    auto ta1_z_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 300);

    auto ta1_z_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 301);

    auto ta1_z_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 302);

    auto ta1_z_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 303);

    auto ta1_z_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 304);

    auto ta1_z_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 305);

    auto ta1_z_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 306);

    auto ta1_z_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 307);

    auto ta1_z_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 308);

    auto ta1_z_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 309);

    auto ta1_z_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 310);

    auto ta1_z_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 311);

    auto ta1_z_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 312);

    auto ta1_z_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 313);

    auto ta1_z_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 315);

    auto ta1_z_xxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 316);

    auto ta1_z_xxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 317);

    auto ta1_z_xxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 318);

    auto ta1_z_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 320);

    auto ta1_z_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 321);

    auto ta1_z_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 322);

    auto ta1_z_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 323);

    auto ta1_z_xxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 324);

    auto ta1_z_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 325);

    auto ta1_z_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 327);

    auto ta1_z_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 328);

    auto ta1_z_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 329);

    auto ta1_z_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 330);

    auto ta1_z_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 331);

    auto ta1_z_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 332);

    auto ta1_z_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 333);

    auto ta1_z_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 334);

    auto ta1_z_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 335);

    auto ta1_z_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 336);

    auto ta1_z_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 337);

    auto ta1_z_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 338);

    auto ta1_z_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 339);

    auto ta1_z_xxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 342);

    auto ta1_z_xxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 345);

    auto ta1_z_xxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 347);

    auto ta1_z_xxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 348);

    auto ta1_z_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 350);

    auto ta1_z_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 351);

    auto ta1_z_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 352);

    auto ta1_z_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 353);

    auto ta1_z_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 354);

    auto ta1_z_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 355);

    auto ta1_z_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 356);

    auto ta1_z_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 357);

    auto ta1_z_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 358);

    auto ta1_z_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 359);

    auto ta1_z_xyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 360);

    auto ta1_z_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 361);

    auto ta1_z_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 363);

    auto ta1_z_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 364);

    auto ta1_z_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 366);

    auto ta1_z_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 367);

    auto ta1_z_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 368);

    auto ta1_z_xyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 369);

    auto ta1_z_xyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 377);

    auto ta1_z_xyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 378);

    auto ta1_z_xyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 379);

    auto ta1_z_xyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 386);

    auto ta1_z_xyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 387);

    auto ta1_z_xyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 388);

    auto ta1_z_xzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 390);

    auto ta1_z_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 392);

    auto ta1_z_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 394);

    auto ta1_z_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 395);

    auto ta1_z_xzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 396);

    auto ta1_z_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 397);

    auto ta1_z_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 398);

    auto ta1_z_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 399);

    auto ta1_z_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 400);

    auto ta1_z_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 401);

    auto ta1_z_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 402);

    auto ta1_z_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 403);

    auto ta1_z_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 404);

    auto ta1_z_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 405);

    auto ta1_z_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 406);

    auto ta1_z_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 407);

    auto ta1_z_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 408);

    auto ta1_z_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 409);

    auto ta1_z_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 411);

    auto ta1_z_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 412);

    auto ta1_z_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 413);

    auto ta1_z_yyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 414);

    auto ta1_z_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 415);

    auto ta1_z_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 416);

    auto ta1_z_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 417);

    auto ta1_z_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 418);

    auto ta1_z_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 419);

    auto ta1_z_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 420);

    auto ta1_z_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 421);

    auto ta1_z_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 422);

    auto ta1_z_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 423);

    auto ta1_z_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 424);

    auto ta1_z_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 425);

    auto ta1_z_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 426);

    auto ta1_z_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 427);

    auto ta1_z_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 428);

    auto ta1_z_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 429);

    auto ta1_z_yzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 430);

    auto ta1_z_yzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 431);

    auto ta1_z_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 432);

    auto ta1_z_yzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 433);

    auto ta1_z_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 434);

    auto ta1_z_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 435);

    auto ta1_z_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 436);

    auto ta1_z_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 437);

    auto ta1_z_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 438);

    auto ta1_z_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 439);

    auto ta1_z_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 440);

    auto ta1_z_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 441);

    auto ta1_z_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 442);

    auto ta1_z_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 443);

    auto ta1_z_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 444);

    auto ta1_z_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 445);

    auto ta1_z_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 446);

    auto ta1_z_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 447);

    auto ta1_z_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 448);

    auto ta1_z_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 449);

    // Set up 0-10 components of targeted buffer : HF

    auto ta1_x_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf);

    auto ta1_x_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 1);

    auto ta1_x_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 2);

    auto ta1_x_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 3);

    auto ta1_x_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 4);

    auto ta1_x_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 5);

    auto ta1_x_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 6);

    auto ta1_x_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 7);

    auto ta1_x_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 8);

    auto ta1_x_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 9);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xxx_xxx_0, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_0, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_0, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_0, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_0, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_0, ta1_x_xxx_xzz_1, ta1_x_xxx_yyy_0, ta1_x_xxx_yyy_1, ta1_x_xxx_yyz_0, ta1_x_xxx_yyz_1, ta1_x_xxx_yzz_0, ta1_x_xxx_yzz_1, ta1_x_xxx_zzz_0, ta1_x_xxx_zzz_1, ta1_x_xxxx_xx_0, ta1_x_xxxx_xx_1, ta1_x_xxxx_xxx_0, ta1_x_xxxx_xxx_1, ta1_x_xxxx_xxy_0, ta1_x_xxxx_xxy_1, ta1_x_xxxx_xxz_0, ta1_x_xxxx_xxz_1, ta1_x_xxxx_xy_0, ta1_x_xxxx_xy_1, ta1_x_xxxx_xyy_0, ta1_x_xxxx_xyy_1, ta1_x_xxxx_xyz_0, ta1_x_xxxx_xyz_1, ta1_x_xxxx_xz_0, ta1_x_xxxx_xz_1, ta1_x_xxxx_xzz_0, ta1_x_xxxx_xzz_1, ta1_x_xxxx_yy_0, ta1_x_xxxx_yy_1, ta1_x_xxxx_yyy_0, ta1_x_xxxx_yyy_1, ta1_x_xxxx_yyz_0, ta1_x_xxxx_yyz_1, ta1_x_xxxx_yz_0, ta1_x_xxxx_yz_1, ta1_x_xxxx_yzz_0, ta1_x_xxxx_yzz_1, ta1_x_xxxx_zz_0, ta1_x_xxxx_zz_1, ta1_x_xxxx_zzz_0, ta1_x_xxxx_zzz_1, ta1_x_xxxxx_xxx_0, ta1_x_xxxxx_xxy_0, ta1_x_xxxxx_xxz_0, ta1_x_xxxxx_xyy_0, ta1_x_xxxxx_xyz_0, ta1_x_xxxxx_xzz_0, ta1_x_xxxxx_yyy_0, ta1_x_xxxxx_yyz_0, ta1_x_xxxxx_yzz_0, ta1_x_xxxxx_zzz_0, ta_xxxx_xxx_1, ta_xxxx_xxy_1, ta_xxxx_xxz_1, ta_xxxx_xyy_1, ta_xxxx_xyz_1, ta_xxxx_xzz_1, ta_xxxx_yyy_1, ta_xxxx_yyz_1, ta_xxxx_yzz_1, ta_xxxx_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxx_xxx_0[i] = 4.0 * ta1_x_xxx_xxx_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxx_1[i] * fe_0 + 3.0 * ta1_x_xxxx_xx_0[i] * fe_0 - 3.0 * ta1_x_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxx_1[i] + ta1_x_xxxx_xxx_0[i] * pa_x[i] - ta1_x_xxxx_xxx_1[i] * pc_x[i];

        ta1_x_xxxxx_xxy_0[i] = 4.0 * ta1_x_xxx_xxy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_x_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xy_1[i] * fe_0 + ta_xxxx_xxy_1[i] + ta1_x_xxxx_xxy_0[i] * pa_x[i] - ta1_x_xxxx_xxy_1[i] * pc_x[i];

        ta1_x_xxxxx_xxz_0[i] = 4.0 * ta1_x_xxx_xxz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_x_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xz_1[i] * fe_0 + ta_xxxx_xxz_1[i] + ta1_x_xxxx_xxz_0[i] * pa_x[i] - ta1_x_xxxx_xxz_1[i] * pc_x[i];

        ta1_x_xxxxx_xyy_0[i] = 4.0 * ta1_x_xxx_xyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyy_1[i] * fe_0 + ta1_x_xxxx_yy_0[i] * fe_0 - ta1_x_xxxx_yy_1[i] * fe_0 + ta_xxxx_xyy_1[i] + ta1_x_xxxx_xyy_0[i] * pa_x[i] - ta1_x_xxxx_xyy_1[i] * pc_x[i];

        ta1_x_xxxxx_xyz_0[i] = 4.0 * ta1_x_xxx_xyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyz_1[i] * fe_0 + ta1_x_xxxx_yz_0[i] * fe_0 - ta1_x_xxxx_yz_1[i] * fe_0 + ta_xxxx_xyz_1[i] + ta1_x_xxxx_xyz_0[i] * pa_x[i] - ta1_x_xxxx_xyz_1[i] * pc_x[i];

        ta1_x_xxxxx_xzz_0[i] = 4.0 * ta1_x_xxx_xzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xzz_1[i] * fe_0 + ta1_x_xxxx_zz_0[i] * fe_0 - ta1_x_xxxx_zz_1[i] * fe_0 + ta_xxxx_xzz_1[i] + ta1_x_xxxx_xzz_0[i] * pa_x[i] - ta1_x_xxxx_xzz_1[i] * pc_x[i];

        ta1_x_xxxxx_yyy_0[i] = 4.0 * ta1_x_xxx_yyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyy_1[i] * fe_0 + ta_xxxx_yyy_1[i] + ta1_x_xxxx_yyy_0[i] * pa_x[i] - ta1_x_xxxx_yyy_1[i] * pc_x[i];

        ta1_x_xxxxx_yyz_0[i] = 4.0 * ta1_x_xxx_yyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyz_1[i] * fe_0 + ta_xxxx_yyz_1[i] + ta1_x_xxxx_yyz_0[i] * pa_x[i] - ta1_x_xxxx_yyz_1[i] * pc_x[i];

        ta1_x_xxxxx_yzz_0[i] = 4.0 * ta1_x_xxx_yzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yzz_1[i] * fe_0 + ta_xxxx_yzz_1[i] + ta1_x_xxxx_yzz_0[i] * pa_x[i] - ta1_x_xxxx_yzz_1[i] * pc_x[i];

        ta1_x_xxxxx_zzz_0[i] = 4.0 * ta1_x_xxx_zzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_zzz_1[i] * fe_0 + ta_xxxx_zzz_1[i] + ta1_x_xxxx_zzz_0[i] * pa_x[i] - ta1_x_xxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : HF

    auto ta1_x_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 10);

    auto ta1_x_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 11);

    auto ta1_x_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 12);

    auto ta1_x_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 13);

    auto ta1_x_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 14);

    auto ta1_x_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 15);

    auto ta1_x_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 16);

    auto ta1_x_xxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 17);

    auto ta1_x_xxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 18);

    auto ta1_x_xxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 19);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxxx_xx_0, ta1_x_xxxx_xx_1, ta1_x_xxxx_xxx_0, ta1_x_xxxx_xxx_1, ta1_x_xxxx_xxy_0, ta1_x_xxxx_xxy_1, ta1_x_xxxx_xxz_0, ta1_x_xxxx_xxz_1, ta1_x_xxxx_xy_0, ta1_x_xxxx_xy_1, ta1_x_xxxx_xyy_0, ta1_x_xxxx_xyy_1, ta1_x_xxxx_xyz_0, ta1_x_xxxx_xyz_1, ta1_x_xxxx_xz_0, ta1_x_xxxx_xz_1, ta1_x_xxxx_xzz_0, ta1_x_xxxx_xzz_1, ta1_x_xxxx_yy_0, ta1_x_xxxx_yy_1, ta1_x_xxxx_yyy_0, ta1_x_xxxx_yyy_1, ta1_x_xxxx_yyz_0, ta1_x_xxxx_yyz_1, ta1_x_xxxx_yz_0, ta1_x_xxxx_yz_1, ta1_x_xxxx_yzz_0, ta1_x_xxxx_yzz_1, ta1_x_xxxx_zz_0, ta1_x_xxxx_zz_1, ta1_x_xxxx_zzz_0, ta1_x_xxxx_zzz_1, ta1_x_xxxxy_xxx_0, ta1_x_xxxxy_xxy_0, ta1_x_xxxxy_xxz_0, ta1_x_xxxxy_xyy_0, ta1_x_xxxxy_xyz_0, ta1_x_xxxxy_xzz_0, ta1_x_xxxxy_yyy_0, ta1_x_xxxxy_yyz_0, ta1_x_xxxxy_yzz_0, ta1_x_xxxxy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxy_xxx_0[i] = ta1_x_xxxx_xxx_0[i] * pa_y[i] - ta1_x_xxxx_xxx_1[i] * pc_y[i];

        ta1_x_xxxxy_xxy_0[i] = ta1_x_xxxx_xx_0[i] * fe_0 - ta1_x_xxxx_xx_1[i] * fe_0 + ta1_x_xxxx_xxy_0[i] * pa_y[i] - ta1_x_xxxx_xxy_1[i] * pc_y[i];

        ta1_x_xxxxy_xxz_0[i] = ta1_x_xxxx_xxz_0[i] * pa_y[i] - ta1_x_xxxx_xxz_1[i] * pc_y[i];

        ta1_x_xxxxy_xyy_0[i] = 2.0 * ta1_x_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xy_1[i] * fe_0 + ta1_x_xxxx_xyy_0[i] * pa_y[i] - ta1_x_xxxx_xyy_1[i] * pc_y[i];

        ta1_x_xxxxy_xyz_0[i] = ta1_x_xxxx_xz_0[i] * fe_0 - ta1_x_xxxx_xz_1[i] * fe_0 + ta1_x_xxxx_xyz_0[i] * pa_y[i] - ta1_x_xxxx_xyz_1[i] * pc_y[i];

        ta1_x_xxxxy_xzz_0[i] = ta1_x_xxxx_xzz_0[i] * pa_y[i] - ta1_x_xxxx_xzz_1[i] * pc_y[i];

        ta1_x_xxxxy_yyy_0[i] = 3.0 * ta1_x_xxxx_yy_0[i] * fe_0 - 3.0 * ta1_x_xxxx_yy_1[i] * fe_0 + ta1_x_xxxx_yyy_0[i] * pa_y[i] - ta1_x_xxxx_yyy_1[i] * pc_y[i];

        ta1_x_xxxxy_yyz_0[i] = 2.0 * ta1_x_xxxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_yz_1[i] * fe_0 + ta1_x_xxxx_yyz_0[i] * pa_y[i] - ta1_x_xxxx_yyz_1[i] * pc_y[i];

        ta1_x_xxxxy_yzz_0[i] = ta1_x_xxxx_zz_0[i] * fe_0 - ta1_x_xxxx_zz_1[i] * fe_0 + ta1_x_xxxx_yzz_0[i] * pa_y[i] - ta1_x_xxxx_yzz_1[i] * pc_y[i];

        ta1_x_xxxxy_zzz_0[i] = ta1_x_xxxx_zzz_0[i] * pa_y[i] - ta1_x_xxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : HF

    auto ta1_x_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 20);

    auto ta1_x_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 21);

    auto ta1_x_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 22);

    auto ta1_x_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 23);

    auto ta1_x_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 24);

    auto ta1_x_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 25);

    auto ta1_x_xxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 26);

    auto ta1_x_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 27);

    auto ta1_x_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 28);

    auto ta1_x_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 29);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_xxxx_xx_0, ta1_x_xxxx_xx_1, ta1_x_xxxx_xxx_0, ta1_x_xxxx_xxx_1, ta1_x_xxxx_xxy_0, ta1_x_xxxx_xxy_1, ta1_x_xxxx_xxz_0, ta1_x_xxxx_xxz_1, ta1_x_xxxx_xy_0, ta1_x_xxxx_xy_1, ta1_x_xxxx_xyy_0, ta1_x_xxxx_xyy_1, ta1_x_xxxx_xyz_0, ta1_x_xxxx_xyz_1, ta1_x_xxxx_xz_0, ta1_x_xxxx_xz_1, ta1_x_xxxx_xzz_0, ta1_x_xxxx_xzz_1, ta1_x_xxxx_yy_0, ta1_x_xxxx_yy_1, ta1_x_xxxx_yyy_0, ta1_x_xxxx_yyy_1, ta1_x_xxxx_yyz_0, ta1_x_xxxx_yyz_1, ta1_x_xxxx_yz_0, ta1_x_xxxx_yz_1, ta1_x_xxxx_yzz_0, ta1_x_xxxx_yzz_1, ta1_x_xxxx_zz_0, ta1_x_xxxx_zz_1, ta1_x_xxxx_zzz_0, ta1_x_xxxx_zzz_1, ta1_x_xxxxz_xxx_0, ta1_x_xxxxz_xxy_0, ta1_x_xxxxz_xxz_0, ta1_x_xxxxz_xyy_0, ta1_x_xxxxz_xyz_0, ta1_x_xxxxz_xzz_0, ta1_x_xxxxz_yyy_0, ta1_x_xxxxz_yyz_0, ta1_x_xxxxz_yzz_0, ta1_x_xxxxz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxz_xxx_0[i] = ta1_x_xxxx_xxx_0[i] * pa_z[i] - ta1_x_xxxx_xxx_1[i] * pc_z[i];

        ta1_x_xxxxz_xxy_0[i] = ta1_x_xxxx_xxy_0[i] * pa_z[i] - ta1_x_xxxx_xxy_1[i] * pc_z[i];

        ta1_x_xxxxz_xxz_0[i] = ta1_x_xxxx_xx_0[i] * fe_0 - ta1_x_xxxx_xx_1[i] * fe_0 + ta1_x_xxxx_xxz_0[i] * pa_z[i] - ta1_x_xxxx_xxz_1[i] * pc_z[i];

        ta1_x_xxxxz_xyy_0[i] = ta1_x_xxxx_xyy_0[i] * pa_z[i] - ta1_x_xxxx_xyy_1[i] * pc_z[i];

        ta1_x_xxxxz_xyz_0[i] = ta1_x_xxxx_xy_0[i] * fe_0 - ta1_x_xxxx_xy_1[i] * fe_0 + ta1_x_xxxx_xyz_0[i] * pa_z[i] - ta1_x_xxxx_xyz_1[i] * pc_z[i];

        ta1_x_xxxxz_xzz_0[i] = 2.0 * ta1_x_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_xz_1[i] * fe_0 + ta1_x_xxxx_xzz_0[i] * pa_z[i] - ta1_x_xxxx_xzz_1[i] * pc_z[i];

        ta1_x_xxxxz_yyy_0[i] = ta1_x_xxxx_yyy_0[i] * pa_z[i] - ta1_x_xxxx_yyy_1[i] * pc_z[i];

        ta1_x_xxxxz_yyz_0[i] = ta1_x_xxxx_yy_0[i] * fe_0 - ta1_x_xxxx_yy_1[i] * fe_0 + ta1_x_xxxx_yyz_0[i] * pa_z[i] - ta1_x_xxxx_yyz_1[i] * pc_z[i];

        ta1_x_xxxxz_yzz_0[i] = 2.0 * ta1_x_xxxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxx_yz_1[i] * fe_0 + ta1_x_xxxx_yzz_0[i] * pa_z[i] - ta1_x_xxxx_yzz_1[i] * pc_z[i];

        ta1_x_xxxxz_zzz_0[i] = 3.0 * ta1_x_xxxx_zz_0[i] * fe_0 - 3.0 * ta1_x_xxxx_zz_1[i] * fe_0 + ta1_x_xxxx_zzz_0[i] * pa_z[i] - ta1_x_xxxx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : HF

    auto ta1_x_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 30);

    auto ta1_x_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 31);

    auto ta1_x_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 32);

    auto ta1_x_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 33);

    auto ta1_x_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 34);

    auto ta1_x_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 35);

    auto ta1_x_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 36);

    auto ta1_x_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 37);

    auto ta1_x_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 38);

    auto ta1_x_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 39);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxx_xxx_0, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_0, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_0, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_0, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_0, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_0, ta1_x_xxx_xzz_1, ta1_x_xxx_zzz_0, ta1_x_xxx_zzz_1, ta1_x_xxxy_xx_0, ta1_x_xxxy_xx_1, ta1_x_xxxy_xxx_0, ta1_x_xxxy_xxx_1, ta1_x_xxxy_xxy_0, ta1_x_xxxy_xxy_1, ta1_x_xxxy_xxz_0, ta1_x_xxxy_xxz_1, ta1_x_xxxy_xy_0, ta1_x_xxxy_xy_1, ta1_x_xxxy_xyy_0, ta1_x_xxxy_xyy_1, ta1_x_xxxy_xyz_0, ta1_x_xxxy_xyz_1, ta1_x_xxxy_xz_0, ta1_x_xxxy_xz_1, ta1_x_xxxy_xzz_0, ta1_x_xxxy_xzz_1, ta1_x_xxxy_zzz_0, ta1_x_xxxy_zzz_1, ta1_x_xxxyy_xxx_0, ta1_x_xxxyy_xxy_0, ta1_x_xxxyy_xxz_0, ta1_x_xxxyy_xyy_0, ta1_x_xxxyy_xyz_0, ta1_x_xxxyy_xzz_0, ta1_x_xxxyy_yyy_0, ta1_x_xxxyy_yyz_0, ta1_x_xxxyy_yzz_0, ta1_x_xxxyy_zzz_0, ta1_x_xxyy_yyy_0, ta1_x_xxyy_yyy_1, ta1_x_xxyy_yyz_0, ta1_x_xxyy_yyz_1, ta1_x_xxyy_yzz_0, ta1_x_xxyy_yzz_1, ta1_x_xyy_yyy_0, ta1_x_xyy_yyy_1, ta1_x_xyy_yyz_0, ta1_x_xyy_yyz_1, ta1_x_xyy_yzz_0, ta1_x_xyy_yzz_1, ta_xxyy_yyy_1, ta_xxyy_yyz_1, ta_xxyy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyy_xxx_0[i] = ta1_x_xxx_xxx_0[i] * fe_0 - ta1_x_xxx_xxx_1[i] * fe_0 + ta1_x_xxxy_xxx_0[i] * pa_y[i] - ta1_x_xxxy_xxx_1[i] * pc_y[i];

        ta1_x_xxxyy_xxy_0[i] = ta1_x_xxx_xxy_0[i] * fe_0 - ta1_x_xxx_xxy_1[i] * fe_0 + ta1_x_xxxy_xx_0[i] * fe_0 - ta1_x_xxxy_xx_1[i] * fe_0 + ta1_x_xxxy_xxy_0[i] * pa_y[i] - ta1_x_xxxy_xxy_1[i] * pc_y[i];

        ta1_x_xxxyy_xxz_0[i] = ta1_x_xxx_xxz_0[i] * fe_0 - ta1_x_xxx_xxz_1[i] * fe_0 + ta1_x_xxxy_xxz_0[i] * pa_y[i] - ta1_x_xxxy_xxz_1[i] * pc_y[i];

        ta1_x_xxxyy_xyy_0[i] = ta1_x_xxx_xyy_0[i] * fe_0 - ta1_x_xxx_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxxy_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xy_1[i] * fe_0 + ta1_x_xxxy_xyy_0[i] * pa_y[i] - ta1_x_xxxy_xyy_1[i] * pc_y[i];

        ta1_x_xxxyy_xyz_0[i] = ta1_x_xxx_xyz_0[i] * fe_0 - ta1_x_xxx_xyz_1[i] * fe_0 + ta1_x_xxxy_xz_0[i] * fe_0 - ta1_x_xxxy_xz_1[i] * fe_0 + ta1_x_xxxy_xyz_0[i] * pa_y[i] - ta1_x_xxxy_xyz_1[i] * pc_y[i];

        ta1_x_xxxyy_xzz_0[i] = ta1_x_xxx_xzz_0[i] * fe_0 - ta1_x_xxx_xzz_1[i] * fe_0 + ta1_x_xxxy_xzz_0[i] * pa_y[i] - ta1_x_xxxy_xzz_1[i] * pc_y[i];

        ta1_x_xxxyy_yyy_0[i] = 2.0 * ta1_x_xyy_yyy_0[i] * fe_0 - 2.0 * ta1_x_xyy_yyy_1[i] * fe_0 + ta_xxyy_yyy_1[i] + ta1_x_xxyy_yyy_0[i] * pa_x[i] - ta1_x_xxyy_yyy_1[i] * pc_x[i];

        ta1_x_xxxyy_yyz_0[i] = 2.0 * ta1_x_xyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yyz_1[i] * fe_0 + ta_xxyy_yyz_1[i] + ta1_x_xxyy_yyz_0[i] * pa_x[i] - ta1_x_xxyy_yyz_1[i] * pc_x[i];

        ta1_x_xxxyy_yzz_0[i] = 2.0 * ta1_x_xyy_yzz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yzz_1[i] * fe_0 + ta_xxyy_yzz_1[i] + ta1_x_xxyy_yzz_0[i] * pa_x[i] - ta1_x_xxyy_yzz_1[i] * pc_x[i];

        ta1_x_xxxyy_zzz_0[i] = ta1_x_xxx_zzz_0[i] * fe_0 - ta1_x_xxx_zzz_1[i] * fe_0 + ta1_x_xxxy_zzz_0[i] * pa_y[i] - ta1_x_xxxy_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : HF

    auto ta1_x_xxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 40);

    auto ta1_x_xxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 41);

    auto ta1_x_xxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 42);

    auto ta1_x_xxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 43);

    auto ta1_x_xxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 44);

    auto ta1_x_xxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 45);

    auto ta1_x_xxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 46);

    auto ta1_x_xxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 47);

    auto ta1_x_xxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 48);

    auto ta1_x_xxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 49);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxxy_xxy_0, ta1_x_xxxy_xxy_1, ta1_x_xxxy_xyy_0, ta1_x_xxxy_xyy_1, ta1_x_xxxy_yyy_0, ta1_x_xxxy_yyy_1, ta1_x_xxxyz_xxx_0, ta1_x_xxxyz_xxy_0, ta1_x_xxxyz_xxz_0, ta1_x_xxxyz_xyy_0, ta1_x_xxxyz_xyz_0, ta1_x_xxxyz_xzz_0, ta1_x_xxxyz_yyy_0, ta1_x_xxxyz_yyz_0, ta1_x_xxxyz_yzz_0, ta1_x_xxxyz_zzz_0, ta1_x_xxxz_xxx_0, ta1_x_xxxz_xxx_1, ta1_x_xxxz_xxz_0, ta1_x_xxxz_xxz_1, ta1_x_xxxz_xyz_0, ta1_x_xxxz_xyz_1, ta1_x_xxxz_xz_0, ta1_x_xxxz_xz_1, ta1_x_xxxz_xzz_0, ta1_x_xxxz_xzz_1, ta1_x_xxxz_yyz_0, ta1_x_xxxz_yyz_1, ta1_x_xxxz_yz_0, ta1_x_xxxz_yz_1, ta1_x_xxxz_yzz_0, ta1_x_xxxz_yzz_1, ta1_x_xxxz_zz_0, ta1_x_xxxz_zz_1, ta1_x_xxxz_zzz_0, ta1_x_xxxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyz_xxx_0[i] = ta1_x_xxxz_xxx_0[i] * pa_y[i] - ta1_x_xxxz_xxx_1[i] * pc_y[i];

        ta1_x_xxxyz_xxy_0[i] = ta1_x_xxxy_xxy_0[i] * pa_z[i] - ta1_x_xxxy_xxy_1[i] * pc_z[i];

        ta1_x_xxxyz_xxz_0[i] = ta1_x_xxxz_xxz_0[i] * pa_y[i] - ta1_x_xxxz_xxz_1[i] * pc_y[i];

        ta1_x_xxxyz_xyy_0[i] = ta1_x_xxxy_xyy_0[i] * pa_z[i] - ta1_x_xxxy_xyy_1[i] * pc_z[i];

        ta1_x_xxxyz_xyz_0[i] = ta1_x_xxxz_xz_0[i] * fe_0 - ta1_x_xxxz_xz_1[i] * fe_0 + ta1_x_xxxz_xyz_0[i] * pa_y[i] - ta1_x_xxxz_xyz_1[i] * pc_y[i];

        ta1_x_xxxyz_xzz_0[i] = ta1_x_xxxz_xzz_0[i] * pa_y[i] - ta1_x_xxxz_xzz_1[i] * pc_y[i];

        ta1_x_xxxyz_yyy_0[i] = ta1_x_xxxy_yyy_0[i] * pa_z[i] - ta1_x_xxxy_yyy_1[i] * pc_z[i];

        ta1_x_xxxyz_yyz_0[i] = 2.0 * ta1_x_xxxz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_yz_1[i] * fe_0 + ta1_x_xxxz_yyz_0[i] * pa_y[i] - ta1_x_xxxz_yyz_1[i] * pc_y[i];

        ta1_x_xxxyz_yzz_0[i] = ta1_x_xxxz_zz_0[i] * fe_0 - ta1_x_xxxz_zz_1[i] * fe_0 + ta1_x_xxxz_yzz_0[i] * pa_y[i] - ta1_x_xxxz_yzz_1[i] * pc_y[i];

        ta1_x_xxxyz_zzz_0[i] = ta1_x_xxxz_zzz_0[i] * pa_y[i] - ta1_x_xxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : HF

    auto ta1_x_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 50);

    auto ta1_x_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 51);

    auto ta1_x_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 52);

    auto ta1_x_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 53);

    auto ta1_x_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 54);

    auto ta1_x_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 55);

    auto ta1_x_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 56);

    auto ta1_x_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 57);

    auto ta1_x_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 58);

    auto ta1_x_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 59);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxx_xxx_0, ta1_x_xxx_xxx_1, ta1_x_xxx_xxy_0, ta1_x_xxx_xxy_1, ta1_x_xxx_xxz_0, ta1_x_xxx_xxz_1, ta1_x_xxx_xyy_0, ta1_x_xxx_xyy_1, ta1_x_xxx_xyz_0, ta1_x_xxx_xyz_1, ta1_x_xxx_xzz_0, ta1_x_xxx_xzz_1, ta1_x_xxx_yyy_0, ta1_x_xxx_yyy_1, ta1_x_xxxz_xx_0, ta1_x_xxxz_xx_1, ta1_x_xxxz_xxx_0, ta1_x_xxxz_xxx_1, ta1_x_xxxz_xxy_0, ta1_x_xxxz_xxy_1, ta1_x_xxxz_xxz_0, ta1_x_xxxz_xxz_1, ta1_x_xxxz_xy_0, ta1_x_xxxz_xy_1, ta1_x_xxxz_xyy_0, ta1_x_xxxz_xyy_1, ta1_x_xxxz_xyz_0, ta1_x_xxxz_xyz_1, ta1_x_xxxz_xz_0, ta1_x_xxxz_xz_1, ta1_x_xxxz_xzz_0, ta1_x_xxxz_xzz_1, ta1_x_xxxz_yyy_0, ta1_x_xxxz_yyy_1, ta1_x_xxxzz_xxx_0, ta1_x_xxxzz_xxy_0, ta1_x_xxxzz_xxz_0, ta1_x_xxxzz_xyy_0, ta1_x_xxxzz_xyz_0, ta1_x_xxxzz_xzz_0, ta1_x_xxxzz_yyy_0, ta1_x_xxxzz_yyz_0, ta1_x_xxxzz_yzz_0, ta1_x_xxxzz_zzz_0, ta1_x_xxzz_yyz_0, ta1_x_xxzz_yyz_1, ta1_x_xxzz_yzz_0, ta1_x_xxzz_yzz_1, ta1_x_xxzz_zzz_0, ta1_x_xxzz_zzz_1, ta1_x_xzz_yyz_0, ta1_x_xzz_yyz_1, ta1_x_xzz_yzz_0, ta1_x_xzz_yzz_1, ta1_x_xzz_zzz_0, ta1_x_xzz_zzz_1, ta_xxzz_yyz_1, ta_xxzz_yzz_1, ta_xxzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzz_xxx_0[i] = ta1_x_xxx_xxx_0[i] * fe_0 - ta1_x_xxx_xxx_1[i] * fe_0 + ta1_x_xxxz_xxx_0[i] * pa_z[i] - ta1_x_xxxz_xxx_1[i] * pc_z[i];

        ta1_x_xxxzz_xxy_0[i] = ta1_x_xxx_xxy_0[i] * fe_0 - ta1_x_xxx_xxy_1[i] * fe_0 + ta1_x_xxxz_xxy_0[i] * pa_z[i] - ta1_x_xxxz_xxy_1[i] * pc_z[i];

        ta1_x_xxxzz_xxz_0[i] = ta1_x_xxx_xxz_0[i] * fe_0 - ta1_x_xxx_xxz_1[i] * fe_0 + ta1_x_xxxz_xx_0[i] * fe_0 - ta1_x_xxxz_xx_1[i] * fe_0 + ta1_x_xxxz_xxz_0[i] * pa_z[i] - ta1_x_xxxz_xxz_1[i] * pc_z[i];

        ta1_x_xxxzz_xyy_0[i] = ta1_x_xxx_xyy_0[i] * fe_0 - ta1_x_xxx_xyy_1[i] * fe_0 + ta1_x_xxxz_xyy_0[i] * pa_z[i] - ta1_x_xxxz_xyy_1[i] * pc_z[i];

        ta1_x_xxxzz_xyz_0[i] = ta1_x_xxx_xyz_0[i] * fe_0 - ta1_x_xxx_xyz_1[i] * fe_0 + ta1_x_xxxz_xy_0[i] * fe_0 - ta1_x_xxxz_xy_1[i] * fe_0 + ta1_x_xxxz_xyz_0[i] * pa_z[i] - ta1_x_xxxz_xyz_1[i] * pc_z[i];

        ta1_x_xxxzz_xzz_0[i] = ta1_x_xxx_xzz_0[i] * fe_0 - ta1_x_xxx_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxxz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xz_1[i] * fe_0 + ta1_x_xxxz_xzz_0[i] * pa_z[i] - ta1_x_xxxz_xzz_1[i] * pc_z[i];

        ta1_x_xxxzz_yyy_0[i] = ta1_x_xxx_yyy_0[i] * fe_0 - ta1_x_xxx_yyy_1[i] * fe_0 + ta1_x_xxxz_yyy_0[i] * pa_z[i] - ta1_x_xxxz_yyy_1[i] * pc_z[i];

        ta1_x_xxxzz_yyz_0[i] = 2.0 * ta1_x_xzz_yyz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yyz_1[i] * fe_0 + ta_xxzz_yyz_1[i] + ta1_x_xxzz_yyz_0[i] * pa_x[i] - ta1_x_xxzz_yyz_1[i] * pc_x[i];

        ta1_x_xxxzz_yzz_0[i] = 2.0 * ta1_x_xzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yzz_1[i] * fe_0 + ta_xxzz_yzz_1[i] + ta1_x_xxzz_yzz_0[i] * pa_x[i] - ta1_x_xxzz_yzz_1[i] * pc_x[i];

        ta1_x_xxxzz_zzz_0[i] = 2.0 * ta1_x_xzz_zzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_zzz_1[i] * fe_0 + ta_xxzz_zzz_1[i] + ta1_x_xxzz_zzz_0[i] * pa_x[i] - ta1_x_xxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : HF

    auto ta1_x_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 60);

    auto ta1_x_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 61);

    auto ta1_x_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 62);

    auto ta1_x_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 63);

    auto ta1_x_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 64);

    auto ta1_x_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 65);

    auto ta1_x_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 66);

    auto ta1_x_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 67);

    auto ta1_x_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 68);

    auto ta1_x_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 69);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxy_xxx_0, ta1_x_xxy_xxx_1, ta1_x_xxy_xxy_0, ta1_x_xxy_xxy_1, ta1_x_xxy_xxz_0, ta1_x_xxy_xxz_1, ta1_x_xxy_xyy_0, ta1_x_xxy_xyy_1, ta1_x_xxy_xyz_0, ta1_x_xxy_xyz_1, ta1_x_xxy_xzz_0, ta1_x_xxy_xzz_1, ta1_x_xxy_zzz_0, ta1_x_xxy_zzz_1, ta1_x_xxyy_xx_0, ta1_x_xxyy_xx_1, ta1_x_xxyy_xxx_0, ta1_x_xxyy_xxx_1, ta1_x_xxyy_xxy_0, ta1_x_xxyy_xxy_1, ta1_x_xxyy_xxz_0, ta1_x_xxyy_xxz_1, ta1_x_xxyy_xy_0, ta1_x_xxyy_xy_1, ta1_x_xxyy_xyy_0, ta1_x_xxyy_xyy_1, ta1_x_xxyy_xyz_0, ta1_x_xxyy_xyz_1, ta1_x_xxyy_xz_0, ta1_x_xxyy_xz_1, ta1_x_xxyy_xzz_0, ta1_x_xxyy_xzz_1, ta1_x_xxyy_zzz_0, ta1_x_xxyy_zzz_1, ta1_x_xxyyy_xxx_0, ta1_x_xxyyy_xxy_0, ta1_x_xxyyy_xxz_0, ta1_x_xxyyy_xyy_0, ta1_x_xxyyy_xyz_0, ta1_x_xxyyy_xzz_0, ta1_x_xxyyy_yyy_0, ta1_x_xxyyy_yyz_0, ta1_x_xxyyy_yzz_0, ta1_x_xxyyy_zzz_0, ta1_x_xyyy_yyy_0, ta1_x_xyyy_yyy_1, ta1_x_xyyy_yyz_0, ta1_x_xyyy_yyz_1, ta1_x_xyyy_yzz_0, ta1_x_xyyy_yzz_1, ta1_x_yyy_yyy_0, ta1_x_yyy_yyy_1, ta1_x_yyy_yyz_0, ta1_x_yyy_yyz_1, ta1_x_yyy_yzz_0, ta1_x_yyy_yzz_1, ta_xyyy_yyy_1, ta_xyyy_yyz_1, ta_xyyy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyy_xxx_0[i] = 2.0 * ta1_x_xxy_xxx_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxx_1[i] * fe_0 + ta1_x_xxyy_xxx_0[i] * pa_y[i] - ta1_x_xxyy_xxx_1[i] * pc_y[i];

        ta1_x_xxyyy_xxy_0[i] = 2.0 * ta1_x_xxy_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxy_1[i] * fe_0 + ta1_x_xxyy_xx_0[i] * fe_0 - ta1_x_xxyy_xx_1[i] * fe_0 + ta1_x_xxyy_xxy_0[i] * pa_y[i] - ta1_x_xxyy_xxy_1[i] * pc_y[i];

        ta1_x_xxyyy_xxz_0[i] = 2.0 * ta1_x_xxy_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xxz_1[i] * fe_0 + ta1_x_xxyy_xxz_0[i] * pa_y[i] - ta1_x_xxyy_xxz_1[i] * pc_y[i];

        ta1_x_xxyyy_xyy_0[i] = 2.0 * ta1_x_xxy_xyy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxyy_xy_0[i] * fe_0 - 2.0 * ta1_x_xxyy_xy_1[i] * fe_0 + ta1_x_xxyy_xyy_0[i] * pa_y[i] - ta1_x_xxyy_xyy_1[i] * pc_y[i];

        ta1_x_xxyyy_xyz_0[i] = 2.0 * ta1_x_xxy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xyz_1[i] * fe_0 + ta1_x_xxyy_xz_0[i] * fe_0 - ta1_x_xxyy_xz_1[i] * fe_0 + ta1_x_xxyy_xyz_0[i] * pa_y[i] - ta1_x_xxyy_xyz_1[i] * pc_y[i];

        ta1_x_xxyyy_xzz_0[i] = 2.0 * ta1_x_xxy_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xzz_1[i] * fe_0 + ta1_x_xxyy_xzz_0[i] * pa_y[i] - ta1_x_xxyy_xzz_1[i] * pc_y[i];

        ta1_x_xxyyy_yyy_0[i] = ta1_x_yyy_yyy_0[i] * fe_0 - ta1_x_yyy_yyy_1[i] * fe_0 + ta_xyyy_yyy_1[i] + ta1_x_xyyy_yyy_0[i] * pa_x[i] - ta1_x_xyyy_yyy_1[i] * pc_x[i];

        ta1_x_xxyyy_yyz_0[i] = ta1_x_yyy_yyz_0[i] * fe_0 - ta1_x_yyy_yyz_1[i] * fe_0 + ta_xyyy_yyz_1[i] + ta1_x_xyyy_yyz_0[i] * pa_x[i] - ta1_x_xyyy_yyz_1[i] * pc_x[i];

        ta1_x_xxyyy_yzz_0[i] = ta1_x_yyy_yzz_0[i] * fe_0 - ta1_x_yyy_yzz_1[i] * fe_0 + ta_xyyy_yzz_1[i] + ta1_x_xyyy_yzz_0[i] * pa_x[i] - ta1_x_xyyy_yzz_1[i] * pc_x[i];

        ta1_x_xxyyy_zzz_0[i] = 2.0 * ta1_x_xxy_zzz_0[i] * fe_0 - 2.0 * ta1_x_xxy_zzz_1[i] * fe_0 + ta1_x_xxyy_zzz_0[i] * pa_y[i] - ta1_x_xxyy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : HF

    auto ta1_x_xxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 70);

    auto ta1_x_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 71);

    auto ta1_x_xxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 72);

    auto ta1_x_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 73);

    auto ta1_x_xxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 74);

    auto ta1_x_xxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 75);

    auto ta1_x_xxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 76);

    auto ta1_x_xxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 77);

    auto ta1_x_xxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 78);

    auto ta1_x_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 79);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxyy_xxx_0, ta1_x_xxyy_xxx_1, ta1_x_xxyy_xxy_0, ta1_x_xxyy_xxy_1, ta1_x_xxyy_xy_0, ta1_x_xxyy_xy_1, ta1_x_xxyy_xyy_0, ta1_x_xxyy_xyy_1, ta1_x_xxyy_xyz_0, ta1_x_xxyy_xyz_1, ta1_x_xxyy_yy_0, ta1_x_xxyy_yy_1, ta1_x_xxyy_yyy_0, ta1_x_xxyy_yyy_1, ta1_x_xxyy_yyz_0, ta1_x_xxyy_yyz_1, ta1_x_xxyy_yz_0, ta1_x_xxyy_yz_1, ta1_x_xxyy_yzz_0, ta1_x_xxyy_yzz_1, ta1_x_xxyyz_xxx_0, ta1_x_xxyyz_xxy_0, ta1_x_xxyyz_xxz_0, ta1_x_xxyyz_xyy_0, ta1_x_xxyyz_xyz_0, ta1_x_xxyyz_xzz_0, ta1_x_xxyyz_yyy_0, ta1_x_xxyyz_yyz_0, ta1_x_xxyyz_yzz_0, ta1_x_xxyyz_zzz_0, ta1_x_xxyz_xxz_0, ta1_x_xxyz_xxz_1, ta1_x_xxyz_xzz_0, ta1_x_xxyz_xzz_1, ta1_x_xxyz_zzz_0, ta1_x_xxyz_zzz_1, ta1_x_xxz_xxz_0, ta1_x_xxz_xxz_1, ta1_x_xxz_xzz_0, ta1_x_xxz_xzz_1, ta1_x_xxz_zzz_0, ta1_x_xxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyz_xxx_0[i] = ta1_x_xxyy_xxx_0[i] * pa_z[i] - ta1_x_xxyy_xxx_1[i] * pc_z[i];

        ta1_x_xxyyz_xxy_0[i] = ta1_x_xxyy_xxy_0[i] * pa_z[i] - ta1_x_xxyy_xxy_1[i] * pc_z[i];

        ta1_x_xxyyz_xxz_0[i] = ta1_x_xxz_xxz_0[i] * fe_0 - ta1_x_xxz_xxz_1[i] * fe_0 + ta1_x_xxyz_xxz_0[i] * pa_y[i] - ta1_x_xxyz_xxz_1[i] * pc_y[i];

        ta1_x_xxyyz_xyy_0[i] = ta1_x_xxyy_xyy_0[i] * pa_z[i] - ta1_x_xxyy_xyy_1[i] * pc_z[i];

        ta1_x_xxyyz_xyz_0[i] = ta1_x_xxyy_xy_0[i] * fe_0 - ta1_x_xxyy_xy_1[i] * fe_0 + ta1_x_xxyy_xyz_0[i] * pa_z[i] - ta1_x_xxyy_xyz_1[i] * pc_z[i];

        ta1_x_xxyyz_xzz_0[i] = ta1_x_xxz_xzz_0[i] * fe_0 - ta1_x_xxz_xzz_1[i] * fe_0 + ta1_x_xxyz_xzz_0[i] * pa_y[i] - ta1_x_xxyz_xzz_1[i] * pc_y[i];

        ta1_x_xxyyz_yyy_0[i] = ta1_x_xxyy_yyy_0[i] * pa_z[i] - ta1_x_xxyy_yyy_1[i] * pc_z[i];

        ta1_x_xxyyz_yyz_0[i] = ta1_x_xxyy_yy_0[i] * fe_0 - ta1_x_xxyy_yy_1[i] * fe_0 + ta1_x_xxyy_yyz_0[i] * pa_z[i] - ta1_x_xxyy_yyz_1[i] * pc_z[i];

        ta1_x_xxyyz_yzz_0[i] = 2.0 * ta1_x_xxyy_yz_0[i] * fe_0 - 2.0 * ta1_x_xxyy_yz_1[i] * fe_0 + ta1_x_xxyy_yzz_0[i] * pa_z[i] - ta1_x_xxyy_yzz_1[i] * pc_z[i];

        ta1_x_xxyyz_zzz_0[i] = ta1_x_xxz_zzz_0[i] * fe_0 - ta1_x_xxz_zzz_1[i] * fe_0 + ta1_x_xxyz_zzz_0[i] * pa_y[i] - ta1_x_xxyz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : HF

    auto ta1_x_xxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 80);

    auto ta1_x_xxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 81);

    auto ta1_x_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 82);

    auto ta1_x_xxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 83);

    auto ta1_x_xxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 84);

    auto ta1_x_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 85);

    auto ta1_x_xxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 86);

    auto ta1_x_xxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 87);

    auto ta1_x_xxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 88);

    auto ta1_x_xxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 89);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxyzz_xxx_0, ta1_x_xxyzz_xxy_0, ta1_x_xxyzz_xxz_0, ta1_x_xxyzz_xyy_0, ta1_x_xxyzz_xyz_0, ta1_x_xxyzz_xzz_0, ta1_x_xxyzz_yyy_0, ta1_x_xxyzz_yyz_0, ta1_x_xxyzz_yzz_0, ta1_x_xxyzz_zzz_0, ta1_x_xxzz_xx_0, ta1_x_xxzz_xx_1, ta1_x_xxzz_xxx_0, ta1_x_xxzz_xxx_1, ta1_x_xxzz_xxy_0, ta1_x_xxzz_xxy_1, ta1_x_xxzz_xxz_0, ta1_x_xxzz_xxz_1, ta1_x_xxzz_xy_0, ta1_x_xxzz_xy_1, ta1_x_xxzz_xyy_0, ta1_x_xxzz_xyy_1, ta1_x_xxzz_xyz_0, ta1_x_xxzz_xyz_1, ta1_x_xxzz_xz_0, ta1_x_xxzz_xz_1, ta1_x_xxzz_xzz_0, ta1_x_xxzz_xzz_1, ta1_x_xxzz_yy_0, ta1_x_xxzz_yy_1, ta1_x_xxzz_yyy_0, ta1_x_xxzz_yyy_1, ta1_x_xxzz_yyz_0, ta1_x_xxzz_yyz_1, ta1_x_xxzz_yz_0, ta1_x_xxzz_yz_1, ta1_x_xxzz_yzz_0, ta1_x_xxzz_yzz_1, ta1_x_xxzz_zz_0, ta1_x_xxzz_zz_1, ta1_x_xxzz_zzz_0, ta1_x_xxzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzz_xxx_0[i] = ta1_x_xxzz_xxx_0[i] * pa_y[i] - ta1_x_xxzz_xxx_1[i] * pc_y[i];

        ta1_x_xxyzz_xxy_0[i] = ta1_x_xxzz_xx_0[i] * fe_0 - ta1_x_xxzz_xx_1[i] * fe_0 + ta1_x_xxzz_xxy_0[i] * pa_y[i] - ta1_x_xxzz_xxy_1[i] * pc_y[i];

        ta1_x_xxyzz_xxz_0[i] = ta1_x_xxzz_xxz_0[i] * pa_y[i] - ta1_x_xxzz_xxz_1[i] * pc_y[i];

        ta1_x_xxyzz_xyy_0[i] = 2.0 * ta1_x_xxzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xxzz_xy_1[i] * fe_0 + ta1_x_xxzz_xyy_0[i] * pa_y[i] - ta1_x_xxzz_xyy_1[i] * pc_y[i];

        ta1_x_xxyzz_xyz_0[i] = ta1_x_xxzz_xz_0[i] * fe_0 - ta1_x_xxzz_xz_1[i] * fe_0 + ta1_x_xxzz_xyz_0[i] * pa_y[i] - ta1_x_xxzz_xyz_1[i] * pc_y[i];

        ta1_x_xxyzz_xzz_0[i] = ta1_x_xxzz_xzz_0[i] * pa_y[i] - ta1_x_xxzz_xzz_1[i] * pc_y[i];

        ta1_x_xxyzz_yyy_0[i] = 3.0 * ta1_x_xxzz_yy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yy_1[i] * fe_0 + ta1_x_xxzz_yyy_0[i] * pa_y[i] - ta1_x_xxzz_yyy_1[i] * pc_y[i];

        ta1_x_xxyzz_yyz_0[i] = 2.0 * ta1_x_xxzz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxzz_yz_1[i] * fe_0 + ta1_x_xxzz_yyz_0[i] * pa_y[i] - ta1_x_xxzz_yyz_1[i] * pc_y[i];

        ta1_x_xxyzz_yzz_0[i] = ta1_x_xxzz_zz_0[i] * fe_0 - ta1_x_xxzz_zz_1[i] * fe_0 + ta1_x_xxzz_yzz_0[i] * pa_y[i] - ta1_x_xxzz_yzz_1[i] * pc_y[i];

        ta1_x_xxyzz_zzz_0[i] = ta1_x_xxzz_zzz_0[i] * pa_y[i] - ta1_x_xxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : HF

    auto ta1_x_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 90);

    auto ta1_x_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 91);

    auto ta1_x_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 92);

    auto ta1_x_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 93);

    auto ta1_x_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 94);

    auto ta1_x_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 95);

    auto ta1_x_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 96);

    auto ta1_x_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 97);

    auto ta1_x_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 98);

    auto ta1_x_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 99);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxz_xxx_0, ta1_x_xxz_xxx_1, ta1_x_xxz_xxy_0, ta1_x_xxz_xxy_1, ta1_x_xxz_xxz_0, ta1_x_xxz_xxz_1, ta1_x_xxz_xyy_0, ta1_x_xxz_xyy_1, ta1_x_xxz_xyz_0, ta1_x_xxz_xyz_1, ta1_x_xxz_xzz_0, ta1_x_xxz_xzz_1, ta1_x_xxz_yyy_0, ta1_x_xxz_yyy_1, ta1_x_xxzz_xx_0, ta1_x_xxzz_xx_1, ta1_x_xxzz_xxx_0, ta1_x_xxzz_xxx_1, ta1_x_xxzz_xxy_0, ta1_x_xxzz_xxy_1, ta1_x_xxzz_xxz_0, ta1_x_xxzz_xxz_1, ta1_x_xxzz_xy_0, ta1_x_xxzz_xy_1, ta1_x_xxzz_xyy_0, ta1_x_xxzz_xyy_1, ta1_x_xxzz_xyz_0, ta1_x_xxzz_xyz_1, ta1_x_xxzz_xz_0, ta1_x_xxzz_xz_1, ta1_x_xxzz_xzz_0, ta1_x_xxzz_xzz_1, ta1_x_xxzz_yyy_0, ta1_x_xxzz_yyy_1, ta1_x_xxzzz_xxx_0, ta1_x_xxzzz_xxy_0, ta1_x_xxzzz_xxz_0, ta1_x_xxzzz_xyy_0, ta1_x_xxzzz_xyz_0, ta1_x_xxzzz_xzz_0, ta1_x_xxzzz_yyy_0, ta1_x_xxzzz_yyz_0, ta1_x_xxzzz_yzz_0, ta1_x_xxzzz_zzz_0, ta1_x_xzzz_yyz_0, ta1_x_xzzz_yyz_1, ta1_x_xzzz_yzz_0, ta1_x_xzzz_yzz_1, ta1_x_xzzz_zzz_0, ta1_x_xzzz_zzz_1, ta1_x_zzz_yyz_0, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_0, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_0, ta1_x_zzz_zzz_1, ta_xzzz_yyz_1, ta_xzzz_yzz_1, ta_xzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzz_xxx_0[i] = 2.0 * ta1_x_xxz_xxx_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxx_1[i] * fe_0 + ta1_x_xxzz_xxx_0[i] * pa_z[i] - ta1_x_xxzz_xxx_1[i] * pc_z[i];

        ta1_x_xxzzz_xxy_0[i] = 2.0 * ta1_x_xxz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxy_1[i] * fe_0 + ta1_x_xxzz_xxy_0[i] * pa_z[i] - ta1_x_xxzz_xxy_1[i] * pc_z[i];

        ta1_x_xxzzz_xxz_0[i] = 2.0 * ta1_x_xxz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxz_1[i] * fe_0 + ta1_x_xxzz_xx_0[i] * fe_0 - ta1_x_xxzz_xx_1[i] * fe_0 + ta1_x_xxzz_xxz_0[i] * pa_z[i] - ta1_x_xxzz_xxz_1[i] * pc_z[i];

        ta1_x_xxzzz_xyy_0[i] = 2.0 * ta1_x_xxz_xyy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyy_1[i] * fe_0 + ta1_x_xxzz_xyy_0[i] * pa_z[i] - ta1_x_xxzz_xyy_1[i] * pc_z[i];

        ta1_x_xxzzz_xyz_0[i] = 2.0 * ta1_x_xxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyz_1[i] * fe_0 + ta1_x_xxzz_xy_0[i] * fe_0 - ta1_x_xxzz_xy_1[i] * fe_0 + ta1_x_xxzz_xyz_0[i] * pa_z[i] - ta1_x_xxzz_xyz_1[i] * pc_z[i];

        ta1_x_xxzzz_xzz_0[i] = 2.0 * ta1_x_xxz_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxzz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxzz_xz_1[i] * fe_0 + ta1_x_xxzz_xzz_0[i] * pa_z[i] - ta1_x_xxzz_xzz_1[i] * pc_z[i];

        ta1_x_xxzzz_yyy_0[i] = 2.0 * ta1_x_xxz_yyy_0[i] * fe_0 - 2.0 * ta1_x_xxz_yyy_1[i] * fe_0 + ta1_x_xxzz_yyy_0[i] * pa_z[i] - ta1_x_xxzz_yyy_1[i] * pc_z[i];

        ta1_x_xxzzz_yyz_0[i] = ta1_x_zzz_yyz_0[i] * fe_0 - ta1_x_zzz_yyz_1[i] * fe_0 + ta_xzzz_yyz_1[i] + ta1_x_xzzz_yyz_0[i] * pa_x[i] - ta1_x_xzzz_yyz_1[i] * pc_x[i];

        ta1_x_xxzzz_yzz_0[i] = ta1_x_zzz_yzz_0[i] * fe_0 - ta1_x_zzz_yzz_1[i] * fe_0 + ta_xzzz_yzz_1[i] + ta1_x_xzzz_yzz_0[i] * pa_x[i] - ta1_x_xzzz_yzz_1[i] * pc_x[i];

        ta1_x_xxzzz_zzz_0[i] = ta1_x_zzz_zzz_0[i] * fe_0 - ta1_x_zzz_zzz_1[i] * fe_0 + ta_xzzz_zzz_1[i] + ta1_x_xzzz_zzz_0[i] * pa_x[i] - ta1_x_xzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : HF

    auto ta1_x_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 100);

    auto ta1_x_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 101);

    auto ta1_x_xyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 102);

    auto ta1_x_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 103);

    auto ta1_x_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 104);

    auto ta1_x_xyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 105);

    auto ta1_x_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 106);

    auto ta1_x_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 107);

    auto ta1_x_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 108);

    auto ta1_x_xyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 109);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyy_xxx_0, ta1_x_xyy_xxx_1, ta1_x_xyy_xxz_0, ta1_x_xyy_xxz_1, ta1_x_xyy_xzz_0, ta1_x_xyy_xzz_1, ta1_x_xyyy_xxx_0, ta1_x_xyyy_xxx_1, ta1_x_xyyy_xxz_0, ta1_x_xyyy_xxz_1, ta1_x_xyyy_xzz_0, ta1_x_xyyy_xzz_1, ta1_x_xyyyy_xxx_0, ta1_x_xyyyy_xxy_0, ta1_x_xyyyy_xxz_0, ta1_x_xyyyy_xyy_0, ta1_x_xyyyy_xyz_0, ta1_x_xyyyy_xzz_0, ta1_x_xyyyy_yyy_0, ta1_x_xyyyy_yyz_0, ta1_x_xyyyy_yzz_0, ta1_x_xyyyy_zzz_0, ta1_x_yyyy_xxy_0, ta1_x_yyyy_xxy_1, ta1_x_yyyy_xy_0, ta1_x_yyyy_xy_1, ta1_x_yyyy_xyy_0, ta1_x_yyyy_xyy_1, ta1_x_yyyy_xyz_0, ta1_x_yyyy_xyz_1, ta1_x_yyyy_yy_0, ta1_x_yyyy_yy_1, ta1_x_yyyy_yyy_0, ta1_x_yyyy_yyy_1, ta1_x_yyyy_yyz_0, ta1_x_yyyy_yyz_1, ta1_x_yyyy_yz_0, ta1_x_yyyy_yz_1, ta1_x_yyyy_yzz_0, ta1_x_yyyy_yzz_1, ta1_x_yyyy_zzz_0, ta1_x_yyyy_zzz_1, ta_yyyy_xxy_1, ta_yyyy_xyy_1, ta_yyyy_xyz_1, ta_yyyy_yyy_1, ta_yyyy_yyz_1, ta_yyyy_yzz_1, ta_yyyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyy_xxx_0[i] = 3.0 * ta1_x_xyy_xxx_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxx_1[i] * fe_0 + ta1_x_xyyy_xxx_0[i] * pa_y[i] - ta1_x_xyyy_xxx_1[i] * pc_y[i];

        ta1_x_xyyyy_xxy_0[i] = 2.0 * ta1_x_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_x_yyyy_xy_1[i] * fe_0 + ta_yyyy_xxy_1[i] + ta1_x_yyyy_xxy_0[i] * pa_x[i] - ta1_x_yyyy_xxy_1[i] * pc_x[i];

        ta1_x_xyyyy_xxz_0[i] = 3.0 * ta1_x_xyy_xxz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxz_1[i] * fe_0 + ta1_x_xyyy_xxz_0[i] * pa_y[i] - ta1_x_xyyy_xxz_1[i] * pc_y[i];

        ta1_x_xyyyy_xyy_0[i] = ta1_x_yyyy_yy_0[i] * fe_0 - ta1_x_yyyy_yy_1[i] * fe_0 + ta_yyyy_xyy_1[i] + ta1_x_yyyy_xyy_0[i] * pa_x[i] - ta1_x_yyyy_xyy_1[i] * pc_x[i];

        ta1_x_xyyyy_xyz_0[i] = ta1_x_yyyy_yz_0[i] * fe_0 - ta1_x_yyyy_yz_1[i] * fe_0 + ta_yyyy_xyz_1[i] + ta1_x_yyyy_xyz_0[i] * pa_x[i] - ta1_x_yyyy_xyz_1[i] * pc_x[i];

        ta1_x_xyyyy_xzz_0[i] = 3.0 * ta1_x_xyy_xzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xzz_1[i] * fe_0 + ta1_x_xyyy_xzz_0[i] * pa_y[i] - ta1_x_xyyy_xzz_1[i] * pc_y[i];

        ta1_x_xyyyy_yyy_0[i] = ta_yyyy_yyy_1[i] + ta1_x_yyyy_yyy_0[i] * pa_x[i] - ta1_x_yyyy_yyy_1[i] * pc_x[i];

        ta1_x_xyyyy_yyz_0[i] = ta_yyyy_yyz_1[i] + ta1_x_yyyy_yyz_0[i] * pa_x[i] - ta1_x_yyyy_yyz_1[i] * pc_x[i];

        ta1_x_xyyyy_yzz_0[i] = ta_yyyy_yzz_1[i] + ta1_x_yyyy_yzz_0[i] * pa_x[i] - ta1_x_yyyy_yzz_1[i] * pc_x[i];

        ta1_x_xyyyy_zzz_0[i] = ta_yyyy_zzz_1[i] + ta1_x_yyyy_zzz_0[i] * pa_x[i] - ta1_x_yyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 110-120 components of targeted buffer : HF

    auto ta1_x_xyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 110);

    auto ta1_x_xyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 111);

    auto ta1_x_xyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 112);

    auto ta1_x_xyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 113);

    auto ta1_x_xyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 114);

    auto ta1_x_xyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 115);

    auto ta1_x_xyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 116);

    auto ta1_x_xyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 117);

    auto ta1_x_xyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 118);

    auto ta1_x_xyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 119);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyyy_xxx_0, ta1_x_xyyy_xxx_1, ta1_x_xyyy_xxy_0, ta1_x_xyyy_xxy_1, ta1_x_xyyy_xy_0, ta1_x_xyyy_xy_1, ta1_x_xyyy_xyy_0, ta1_x_xyyy_xyy_1, ta1_x_xyyy_xyz_0, ta1_x_xyyy_xyz_1, ta1_x_xyyy_yyy_0, ta1_x_xyyy_yyy_1, ta1_x_xyyyz_xxx_0, ta1_x_xyyyz_xxy_0, ta1_x_xyyyz_xxz_0, ta1_x_xyyyz_xyy_0, ta1_x_xyyyz_xyz_0, ta1_x_xyyyz_xzz_0, ta1_x_xyyyz_yyy_0, ta1_x_xyyyz_yyz_0, ta1_x_xyyyz_yzz_0, ta1_x_xyyyz_zzz_0, ta1_x_xyyz_xxz_0, ta1_x_xyyz_xxz_1, ta1_x_xyyz_xzz_0, ta1_x_xyyz_xzz_1, ta1_x_xyz_xxz_0, ta1_x_xyz_xxz_1, ta1_x_xyz_xzz_0, ta1_x_xyz_xzz_1, ta1_x_yyyz_yyz_0, ta1_x_yyyz_yyz_1, ta1_x_yyyz_yzz_0, ta1_x_yyyz_yzz_1, ta1_x_yyyz_zzz_0, ta1_x_yyyz_zzz_1, ta_yyyz_yyz_1, ta_yyyz_yzz_1, ta_yyyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyz_xxx_0[i] = ta1_x_xyyy_xxx_0[i] * pa_z[i] - ta1_x_xyyy_xxx_1[i] * pc_z[i];

        ta1_x_xyyyz_xxy_0[i] = ta1_x_xyyy_xxy_0[i] * pa_z[i] - ta1_x_xyyy_xxy_1[i] * pc_z[i];

        ta1_x_xyyyz_xxz_0[i] = 2.0 * ta1_x_xyz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xxz_1[i] * fe_0 + ta1_x_xyyz_xxz_0[i] * pa_y[i] - ta1_x_xyyz_xxz_1[i] * pc_y[i];

        ta1_x_xyyyz_xyy_0[i] = ta1_x_xyyy_xyy_0[i] * pa_z[i] - ta1_x_xyyy_xyy_1[i] * pc_z[i];

        ta1_x_xyyyz_xyz_0[i] = ta1_x_xyyy_xy_0[i] * fe_0 - ta1_x_xyyy_xy_1[i] * fe_0 + ta1_x_xyyy_xyz_0[i] * pa_z[i] - ta1_x_xyyy_xyz_1[i] * pc_z[i];

        ta1_x_xyyyz_xzz_0[i] = 2.0 * ta1_x_xyz_xzz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xzz_1[i] * fe_0 + ta1_x_xyyz_xzz_0[i] * pa_y[i] - ta1_x_xyyz_xzz_1[i] * pc_y[i];

        ta1_x_xyyyz_yyy_0[i] = ta1_x_xyyy_yyy_0[i] * pa_z[i] - ta1_x_xyyy_yyy_1[i] * pc_z[i];

        ta1_x_xyyyz_yyz_0[i] = ta_yyyz_yyz_1[i] + ta1_x_yyyz_yyz_0[i] * pa_x[i] - ta1_x_yyyz_yyz_1[i] * pc_x[i];

        ta1_x_xyyyz_yzz_0[i] = ta_yyyz_yzz_1[i] + ta1_x_yyyz_yzz_0[i] * pa_x[i] - ta1_x_yyyz_yzz_1[i] * pc_x[i];

        ta1_x_xyyyz_zzz_0[i] = ta_yyyz_zzz_1[i] + ta1_x_yyyz_zzz_0[i] * pa_x[i] - ta1_x_yyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 120-130 components of targeted buffer : HF

    auto ta1_x_xyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 120);

    auto ta1_x_xyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 121);

    auto ta1_x_xyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 122);

    auto ta1_x_xyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 123);

    auto ta1_x_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 124);

    auto ta1_x_xyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 125);

    auto ta1_x_xyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 126);

    auto ta1_x_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 127);

    auto ta1_x_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 128);

    auto ta1_x_xyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 129);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xyy_xxy_0, ta1_x_xyy_xxy_1, ta1_x_xyy_xyy_0, ta1_x_xyy_xyy_1, ta1_x_xyyz_xxy_0, ta1_x_xyyz_xxy_1, ta1_x_xyyz_xyy_0, ta1_x_xyyz_xyy_1, ta1_x_xyyzz_xxx_0, ta1_x_xyyzz_xxy_0, ta1_x_xyyzz_xxz_0, ta1_x_xyyzz_xyy_0, ta1_x_xyyzz_xyz_0, ta1_x_xyyzz_xzz_0, ta1_x_xyyzz_yyy_0, ta1_x_xyyzz_yyz_0, ta1_x_xyyzz_yzz_0, ta1_x_xyyzz_zzz_0, ta1_x_xyzz_xxx_0, ta1_x_xyzz_xxx_1, ta1_x_xyzz_xxz_0, ta1_x_xyzz_xxz_1, ta1_x_xyzz_xzz_0, ta1_x_xyzz_xzz_1, ta1_x_xzz_xxx_0, ta1_x_xzz_xxx_1, ta1_x_xzz_xxz_0, ta1_x_xzz_xxz_1, ta1_x_xzz_xzz_0, ta1_x_xzz_xzz_1, ta1_x_yyzz_xyz_0, ta1_x_yyzz_xyz_1, ta1_x_yyzz_yyy_0, ta1_x_yyzz_yyy_1, ta1_x_yyzz_yyz_0, ta1_x_yyzz_yyz_1, ta1_x_yyzz_yz_0, ta1_x_yyzz_yz_1, ta1_x_yyzz_yzz_0, ta1_x_yyzz_yzz_1, ta1_x_yyzz_zzz_0, ta1_x_yyzz_zzz_1, ta_yyzz_xyz_1, ta_yyzz_yyy_1, ta_yyzz_yyz_1, ta_yyzz_yzz_1, ta_yyzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzz_xxx_0[i] = ta1_x_xzz_xxx_0[i] * fe_0 - ta1_x_xzz_xxx_1[i] * fe_0 + ta1_x_xyzz_xxx_0[i] * pa_y[i] - ta1_x_xyzz_xxx_1[i] * pc_y[i];

        ta1_x_xyyzz_xxy_0[i] = ta1_x_xyy_xxy_0[i] * fe_0 - ta1_x_xyy_xxy_1[i] * fe_0 + ta1_x_xyyz_xxy_0[i] * pa_z[i] - ta1_x_xyyz_xxy_1[i] * pc_z[i];

        ta1_x_xyyzz_xxz_0[i] = ta1_x_xzz_xxz_0[i] * fe_0 - ta1_x_xzz_xxz_1[i] * fe_0 + ta1_x_xyzz_xxz_0[i] * pa_y[i] - ta1_x_xyzz_xxz_1[i] * pc_y[i];

        ta1_x_xyyzz_xyy_0[i] = ta1_x_xyy_xyy_0[i] * fe_0 - ta1_x_xyy_xyy_1[i] * fe_0 + ta1_x_xyyz_xyy_0[i] * pa_z[i] - ta1_x_xyyz_xyy_1[i] * pc_z[i];

        ta1_x_xyyzz_xyz_0[i] = ta1_x_yyzz_yz_0[i] * fe_0 - ta1_x_yyzz_yz_1[i] * fe_0 + ta_yyzz_xyz_1[i] + ta1_x_yyzz_xyz_0[i] * pa_x[i] - ta1_x_yyzz_xyz_1[i] * pc_x[i];

        ta1_x_xyyzz_xzz_0[i] = ta1_x_xzz_xzz_0[i] * fe_0 - ta1_x_xzz_xzz_1[i] * fe_0 + ta1_x_xyzz_xzz_0[i] * pa_y[i] - ta1_x_xyzz_xzz_1[i] * pc_y[i];

        ta1_x_xyyzz_yyy_0[i] = ta_yyzz_yyy_1[i] + ta1_x_yyzz_yyy_0[i] * pa_x[i] - ta1_x_yyzz_yyy_1[i] * pc_x[i];

        ta1_x_xyyzz_yyz_0[i] = ta_yyzz_yyz_1[i] + ta1_x_yyzz_yyz_0[i] * pa_x[i] - ta1_x_yyzz_yyz_1[i] * pc_x[i];

        ta1_x_xyyzz_yzz_0[i] = ta_yyzz_yzz_1[i] + ta1_x_yyzz_yzz_0[i] * pa_x[i] - ta1_x_yyzz_yzz_1[i] * pc_x[i];

        ta1_x_xyyzz_zzz_0[i] = ta_yyzz_zzz_1[i] + ta1_x_yyzz_zzz_0[i] * pa_x[i] - ta1_x_yyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : HF

    auto ta1_x_xyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 130);

    auto ta1_x_xyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 131);

    auto ta1_x_xyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 132);

    auto ta1_x_xyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 133);

    auto ta1_x_xyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 134);

    auto ta1_x_xyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 135);

    auto ta1_x_xyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 136);

    auto ta1_x_xyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 137);

    auto ta1_x_xyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 138);

    auto ta1_x_xyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 139);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyzzz_xxx_0, ta1_x_xyzzz_xxy_0, ta1_x_xyzzz_xxz_0, ta1_x_xyzzz_xyy_0, ta1_x_xyzzz_xyz_0, ta1_x_xyzzz_xzz_0, ta1_x_xyzzz_yyy_0, ta1_x_xyzzz_yyz_0, ta1_x_xyzzz_yzz_0, ta1_x_xyzzz_zzz_0, ta1_x_xzzz_xx_0, ta1_x_xzzz_xx_1, ta1_x_xzzz_xxx_0, ta1_x_xzzz_xxx_1, ta1_x_xzzz_xxy_0, ta1_x_xzzz_xxy_1, ta1_x_xzzz_xxz_0, ta1_x_xzzz_xxz_1, ta1_x_xzzz_xy_0, ta1_x_xzzz_xy_1, ta1_x_xzzz_xyy_0, ta1_x_xzzz_xyy_1, ta1_x_xzzz_xyz_0, ta1_x_xzzz_xyz_1, ta1_x_xzzz_xz_0, ta1_x_xzzz_xz_1, ta1_x_xzzz_xzz_0, ta1_x_xzzz_xzz_1, ta1_x_xzzz_zzz_0, ta1_x_xzzz_zzz_1, ta1_x_yzzz_yyy_0, ta1_x_yzzz_yyy_1, ta1_x_yzzz_yyz_0, ta1_x_yzzz_yyz_1, ta1_x_yzzz_yzz_0, ta1_x_yzzz_yzz_1, ta_yzzz_yyy_1, ta_yzzz_yyz_1, ta_yzzz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzz_xxx_0[i] = ta1_x_xzzz_xxx_0[i] * pa_y[i] - ta1_x_xzzz_xxx_1[i] * pc_y[i];

        ta1_x_xyzzz_xxy_0[i] = ta1_x_xzzz_xx_0[i] * fe_0 - ta1_x_xzzz_xx_1[i] * fe_0 + ta1_x_xzzz_xxy_0[i] * pa_y[i] - ta1_x_xzzz_xxy_1[i] * pc_y[i];

        ta1_x_xyzzz_xxz_0[i] = ta1_x_xzzz_xxz_0[i] * pa_y[i] - ta1_x_xzzz_xxz_1[i] * pc_y[i];

        ta1_x_xyzzz_xyy_0[i] = 2.0 * ta1_x_xzzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xzzz_xy_1[i] * fe_0 + ta1_x_xzzz_xyy_0[i] * pa_y[i] - ta1_x_xzzz_xyy_1[i] * pc_y[i];

        ta1_x_xyzzz_xyz_0[i] = ta1_x_xzzz_xz_0[i] * fe_0 - ta1_x_xzzz_xz_1[i] * fe_0 + ta1_x_xzzz_xyz_0[i] * pa_y[i] - ta1_x_xzzz_xyz_1[i] * pc_y[i];

        ta1_x_xyzzz_xzz_0[i] = ta1_x_xzzz_xzz_0[i] * pa_y[i] - ta1_x_xzzz_xzz_1[i] * pc_y[i];

        ta1_x_xyzzz_yyy_0[i] = ta_yzzz_yyy_1[i] + ta1_x_yzzz_yyy_0[i] * pa_x[i] - ta1_x_yzzz_yyy_1[i] * pc_x[i];

        ta1_x_xyzzz_yyz_0[i] = ta_yzzz_yyz_1[i] + ta1_x_yzzz_yyz_0[i] * pa_x[i] - ta1_x_yzzz_yyz_1[i] * pc_x[i];

        ta1_x_xyzzz_yzz_0[i] = ta_yzzz_yzz_1[i] + ta1_x_yzzz_yzz_0[i] * pa_x[i] - ta1_x_yzzz_yzz_1[i] * pc_x[i];

        ta1_x_xyzzz_zzz_0[i] = ta1_x_xzzz_zzz_0[i] * pa_y[i] - ta1_x_xzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : HF

    auto ta1_x_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 140);

    auto ta1_x_xzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 141);

    auto ta1_x_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 142);

    auto ta1_x_xzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 143);

    auto ta1_x_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 144);

    auto ta1_x_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 145);

    auto ta1_x_xzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 146);

    auto ta1_x_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 147);

    auto ta1_x_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 148);

    auto ta1_x_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 149);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xzz_xxx_0, ta1_x_xzz_xxx_1, ta1_x_xzz_xxy_0, ta1_x_xzz_xxy_1, ta1_x_xzz_xyy_0, ta1_x_xzz_xyy_1, ta1_x_xzzz_xxx_0, ta1_x_xzzz_xxx_1, ta1_x_xzzz_xxy_0, ta1_x_xzzz_xxy_1, ta1_x_xzzz_xyy_0, ta1_x_xzzz_xyy_1, ta1_x_xzzzz_xxx_0, ta1_x_xzzzz_xxy_0, ta1_x_xzzzz_xxz_0, ta1_x_xzzzz_xyy_0, ta1_x_xzzzz_xyz_0, ta1_x_xzzzz_xzz_0, ta1_x_xzzzz_yyy_0, ta1_x_xzzzz_yyz_0, ta1_x_xzzzz_yzz_0, ta1_x_xzzzz_zzz_0, ta1_x_zzzz_xxz_0, ta1_x_zzzz_xxz_1, ta1_x_zzzz_xyz_0, ta1_x_zzzz_xyz_1, ta1_x_zzzz_xz_0, ta1_x_zzzz_xz_1, ta1_x_zzzz_xzz_0, ta1_x_zzzz_xzz_1, ta1_x_zzzz_yyy_0, ta1_x_zzzz_yyy_1, ta1_x_zzzz_yyz_0, ta1_x_zzzz_yyz_1, ta1_x_zzzz_yz_0, ta1_x_zzzz_yz_1, ta1_x_zzzz_yzz_0, ta1_x_zzzz_yzz_1, ta1_x_zzzz_zz_0, ta1_x_zzzz_zz_1, ta1_x_zzzz_zzz_0, ta1_x_zzzz_zzz_1, ta_zzzz_xxz_1, ta_zzzz_xyz_1, ta_zzzz_xzz_1, ta_zzzz_yyy_1, ta_zzzz_yyz_1, ta_zzzz_yzz_1, ta_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzz_xxx_0[i] = 3.0 * ta1_x_xzz_xxx_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxx_1[i] * fe_0 + ta1_x_xzzz_xxx_0[i] * pa_z[i] - ta1_x_xzzz_xxx_1[i] * pc_z[i];

        ta1_x_xzzzz_xxy_0[i] = 3.0 * ta1_x_xzz_xxy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxy_1[i] * fe_0 + ta1_x_xzzz_xxy_0[i] * pa_z[i] - ta1_x_xzzz_xxy_1[i] * pc_z[i];

        ta1_x_xzzzz_xxz_0[i] = 2.0 * ta1_x_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xz_1[i] * fe_0 + ta_zzzz_xxz_1[i] + ta1_x_zzzz_xxz_0[i] * pa_x[i] - ta1_x_zzzz_xxz_1[i] * pc_x[i];

        ta1_x_xzzzz_xyy_0[i] = 3.0 * ta1_x_xzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xyy_1[i] * fe_0 + ta1_x_xzzz_xyy_0[i] * pa_z[i] - ta1_x_xzzz_xyy_1[i] * pc_z[i];

        ta1_x_xzzzz_xyz_0[i] = ta1_x_zzzz_yz_0[i] * fe_0 - ta1_x_zzzz_yz_1[i] * fe_0 + ta_zzzz_xyz_1[i] + ta1_x_zzzz_xyz_0[i] * pa_x[i] - ta1_x_zzzz_xyz_1[i] * pc_x[i];

        ta1_x_xzzzz_xzz_0[i] = ta1_x_zzzz_zz_0[i] * fe_0 - ta1_x_zzzz_zz_1[i] * fe_0 + ta_zzzz_xzz_1[i] + ta1_x_zzzz_xzz_0[i] * pa_x[i] - ta1_x_zzzz_xzz_1[i] * pc_x[i];

        ta1_x_xzzzz_yyy_0[i] = ta_zzzz_yyy_1[i] + ta1_x_zzzz_yyy_0[i] * pa_x[i] - ta1_x_zzzz_yyy_1[i] * pc_x[i];

        ta1_x_xzzzz_yyz_0[i] = ta_zzzz_yyz_1[i] + ta1_x_zzzz_yyz_0[i] * pa_x[i] - ta1_x_zzzz_yyz_1[i] * pc_x[i];

        ta1_x_xzzzz_yzz_0[i] = ta_zzzz_yzz_1[i] + ta1_x_zzzz_yzz_0[i] * pa_x[i] - ta1_x_zzzz_yzz_1[i] * pc_x[i];

        ta1_x_xzzzz_zzz_0[i] = ta_zzzz_zzz_1[i] + ta1_x_zzzz_zzz_0[i] * pa_x[i] - ta1_x_zzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : HF

    auto ta1_x_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 150);

    auto ta1_x_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 151);

    auto ta1_x_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 152);

    auto ta1_x_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 153);

    auto ta1_x_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 154);

    auto ta1_x_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 155);

    auto ta1_x_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 156);

    auto ta1_x_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 157);

    auto ta1_x_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 158);

    auto ta1_x_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 159);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yyy_xxx_0, ta1_x_yyy_xxx_1, ta1_x_yyy_xxy_0, ta1_x_yyy_xxy_1, ta1_x_yyy_xxz_0, ta1_x_yyy_xxz_1, ta1_x_yyy_xyy_0, ta1_x_yyy_xyy_1, ta1_x_yyy_xyz_0, ta1_x_yyy_xyz_1, ta1_x_yyy_xzz_0, ta1_x_yyy_xzz_1, ta1_x_yyy_yyy_0, ta1_x_yyy_yyy_1, ta1_x_yyy_yyz_0, ta1_x_yyy_yyz_1, ta1_x_yyy_yzz_0, ta1_x_yyy_yzz_1, ta1_x_yyy_zzz_0, ta1_x_yyy_zzz_1, ta1_x_yyyy_xx_0, ta1_x_yyyy_xx_1, ta1_x_yyyy_xxx_0, ta1_x_yyyy_xxx_1, ta1_x_yyyy_xxy_0, ta1_x_yyyy_xxy_1, ta1_x_yyyy_xxz_0, ta1_x_yyyy_xxz_1, ta1_x_yyyy_xy_0, ta1_x_yyyy_xy_1, ta1_x_yyyy_xyy_0, ta1_x_yyyy_xyy_1, ta1_x_yyyy_xyz_0, ta1_x_yyyy_xyz_1, ta1_x_yyyy_xz_0, ta1_x_yyyy_xz_1, ta1_x_yyyy_xzz_0, ta1_x_yyyy_xzz_1, ta1_x_yyyy_yy_0, ta1_x_yyyy_yy_1, ta1_x_yyyy_yyy_0, ta1_x_yyyy_yyy_1, ta1_x_yyyy_yyz_0, ta1_x_yyyy_yyz_1, ta1_x_yyyy_yz_0, ta1_x_yyyy_yz_1, ta1_x_yyyy_yzz_0, ta1_x_yyyy_yzz_1, ta1_x_yyyy_zz_0, ta1_x_yyyy_zz_1, ta1_x_yyyy_zzz_0, ta1_x_yyyy_zzz_1, ta1_x_yyyyy_xxx_0, ta1_x_yyyyy_xxy_0, ta1_x_yyyyy_xxz_0, ta1_x_yyyyy_xyy_0, ta1_x_yyyyy_xyz_0, ta1_x_yyyyy_xzz_0, ta1_x_yyyyy_yyy_0, ta1_x_yyyyy_yyz_0, ta1_x_yyyyy_yzz_0, ta1_x_yyyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyy_xxx_0[i] = 4.0 * ta1_x_yyy_xxx_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxx_1[i] * fe_0 + ta1_x_yyyy_xxx_0[i] * pa_y[i] - ta1_x_yyyy_xxx_1[i] * pc_y[i];

        ta1_x_yyyyy_xxy_0[i] = 4.0 * ta1_x_yyy_xxy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxy_1[i] * fe_0 + ta1_x_yyyy_xx_0[i] * fe_0 - ta1_x_yyyy_xx_1[i] * fe_0 + ta1_x_yyyy_xxy_0[i] * pa_y[i] - ta1_x_yyyy_xxy_1[i] * pc_y[i];

        ta1_x_yyyyy_xxz_0[i] = 4.0 * ta1_x_yyy_xxz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxz_1[i] * fe_0 + ta1_x_yyyy_xxz_0[i] * pa_y[i] - ta1_x_yyyy_xxz_1[i] * pc_y[i];

        ta1_x_yyyyy_xyy_0[i] = 4.0 * ta1_x_yyy_xyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_x_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_x_yyyy_xy_1[i] * fe_0 + ta1_x_yyyy_xyy_0[i] * pa_y[i] - ta1_x_yyyy_xyy_1[i] * pc_y[i];

        ta1_x_yyyyy_xyz_0[i] = 4.0 * ta1_x_yyy_xyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyz_1[i] * fe_0 + ta1_x_yyyy_xz_0[i] * fe_0 - ta1_x_yyyy_xz_1[i] * fe_0 + ta1_x_yyyy_xyz_0[i] * pa_y[i] - ta1_x_yyyy_xyz_1[i] * pc_y[i];

        ta1_x_yyyyy_xzz_0[i] = 4.0 * ta1_x_yyy_xzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xzz_1[i] * fe_0 + ta1_x_yyyy_xzz_0[i] * pa_y[i] - ta1_x_yyyy_xzz_1[i] * pc_y[i];

        ta1_x_yyyyy_yyy_0[i] = 4.0 * ta1_x_yyy_yyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyy_1[i] * fe_0 + 3.0 * ta1_x_yyyy_yy_0[i] * fe_0 - 3.0 * ta1_x_yyyy_yy_1[i] * fe_0 + ta1_x_yyyy_yyy_0[i] * pa_y[i] - ta1_x_yyyy_yyy_1[i] * pc_y[i];

        ta1_x_yyyyy_yyz_0[i] = 4.0 * ta1_x_yyy_yyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_x_yyyy_yz_1[i] * fe_0 + ta1_x_yyyy_yyz_0[i] * pa_y[i] - ta1_x_yyyy_yyz_1[i] * pc_y[i];

        ta1_x_yyyyy_yzz_0[i] = 4.0 * ta1_x_yyy_yzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yzz_1[i] * fe_0 + ta1_x_yyyy_zz_0[i] * fe_0 - ta1_x_yyyy_zz_1[i] * fe_0 + ta1_x_yyyy_yzz_0[i] * pa_y[i] - ta1_x_yyyy_yzz_1[i] * pc_y[i];

        ta1_x_yyyyy_zzz_0[i] = 4.0 * ta1_x_yyy_zzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_zzz_1[i] * fe_0 + ta1_x_yyyy_zzz_0[i] * pa_y[i] - ta1_x_yyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 160-170 components of targeted buffer : HF

    auto ta1_x_yyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 160);

    auto ta1_x_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 161);

    auto ta1_x_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 162);

    auto ta1_x_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 163);

    auto ta1_x_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 164);

    auto ta1_x_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 165);

    auto ta1_x_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 166);

    auto ta1_x_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 167);

    auto ta1_x_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 168);

    auto ta1_x_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 169);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyy_xxx_0, ta1_x_yyyy_xxx_1, ta1_x_yyyy_xxy_0, ta1_x_yyyy_xxy_1, ta1_x_yyyy_xy_0, ta1_x_yyyy_xy_1, ta1_x_yyyy_xyy_0, ta1_x_yyyy_xyy_1, ta1_x_yyyy_xyz_0, ta1_x_yyyy_xyz_1, ta1_x_yyyy_yy_0, ta1_x_yyyy_yy_1, ta1_x_yyyy_yyy_0, ta1_x_yyyy_yyy_1, ta1_x_yyyy_yyz_0, ta1_x_yyyy_yyz_1, ta1_x_yyyy_yz_0, ta1_x_yyyy_yz_1, ta1_x_yyyy_yzz_0, ta1_x_yyyy_yzz_1, ta1_x_yyyyz_xxx_0, ta1_x_yyyyz_xxy_0, ta1_x_yyyyz_xxz_0, ta1_x_yyyyz_xyy_0, ta1_x_yyyyz_xyz_0, ta1_x_yyyyz_xzz_0, ta1_x_yyyyz_yyy_0, ta1_x_yyyyz_yyz_0, ta1_x_yyyyz_yzz_0, ta1_x_yyyyz_zzz_0, ta1_x_yyyz_xxz_0, ta1_x_yyyz_xxz_1, ta1_x_yyyz_xzz_0, ta1_x_yyyz_xzz_1, ta1_x_yyyz_zzz_0, ta1_x_yyyz_zzz_1, ta1_x_yyz_xxz_0, ta1_x_yyz_xxz_1, ta1_x_yyz_xzz_0, ta1_x_yyz_xzz_1, ta1_x_yyz_zzz_0, ta1_x_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyz_xxx_0[i] = ta1_x_yyyy_xxx_0[i] * pa_z[i] - ta1_x_yyyy_xxx_1[i] * pc_z[i];

        ta1_x_yyyyz_xxy_0[i] = ta1_x_yyyy_xxy_0[i] * pa_z[i] - ta1_x_yyyy_xxy_1[i] * pc_z[i];

        ta1_x_yyyyz_xxz_0[i] = 3.0 * ta1_x_yyz_xxz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xxz_1[i] * fe_0 + ta1_x_yyyz_xxz_0[i] * pa_y[i] - ta1_x_yyyz_xxz_1[i] * pc_y[i];

        ta1_x_yyyyz_xyy_0[i] = ta1_x_yyyy_xyy_0[i] * pa_z[i] - ta1_x_yyyy_xyy_1[i] * pc_z[i];

        ta1_x_yyyyz_xyz_0[i] = ta1_x_yyyy_xy_0[i] * fe_0 - ta1_x_yyyy_xy_1[i] * fe_0 + ta1_x_yyyy_xyz_0[i] * pa_z[i] - ta1_x_yyyy_xyz_1[i] * pc_z[i];

        ta1_x_yyyyz_xzz_0[i] = 3.0 * ta1_x_yyz_xzz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xzz_1[i] * fe_0 + ta1_x_yyyz_xzz_0[i] * pa_y[i] - ta1_x_yyyz_xzz_1[i] * pc_y[i];

        ta1_x_yyyyz_yyy_0[i] = ta1_x_yyyy_yyy_0[i] * pa_z[i] - ta1_x_yyyy_yyy_1[i] * pc_z[i];

        ta1_x_yyyyz_yyz_0[i] = ta1_x_yyyy_yy_0[i] * fe_0 - ta1_x_yyyy_yy_1[i] * fe_0 + ta1_x_yyyy_yyz_0[i] * pa_z[i] - ta1_x_yyyy_yyz_1[i] * pc_z[i];

        ta1_x_yyyyz_yzz_0[i] = 2.0 * ta1_x_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_x_yyyy_yz_1[i] * fe_0 + ta1_x_yyyy_yzz_0[i] * pa_z[i] - ta1_x_yyyy_yzz_1[i] * pc_z[i];

        ta1_x_yyyyz_zzz_0[i] = 3.0 * ta1_x_yyz_zzz_0[i] * fe_0 - 3.0 * ta1_x_yyz_zzz_1[i] * fe_0 + ta1_x_yyyz_zzz_0[i] * pa_y[i] - ta1_x_yyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : HF

    auto ta1_x_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 170);

    auto ta1_x_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 171);

    auto ta1_x_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 172);

    auto ta1_x_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 173);

    auto ta1_x_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 174);

    auto ta1_x_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 175);

    auto ta1_x_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 176);

    auto ta1_x_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 177);

    auto ta1_x_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 178);

    auto ta1_x_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 179);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyy_xxy_0, ta1_x_yyy_xxy_1, ta1_x_yyy_xyy_0, ta1_x_yyy_xyy_1, ta1_x_yyy_yyy_0, ta1_x_yyy_yyy_1, ta1_x_yyyz_xxy_0, ta1_x_yyyz_xxy_1, ta1_x_yyyz_xyy_0, ta1_x_yyyz_xyy_1, ta1_x_yyyz_yyy_0, ta1_x_yyyz_yyy_1, ta1_x_yyyzz_xxx_0, ta1_x_yyyzz_xxy_0, ta1_x_yyyzz_xxz_0, ta1_x_yyyzz_xyy_0, ta1_x_yyyzz_xyz_0, ta1_x_yyyzz_xzz_0, ta1_x_yyyzz_yyy_0, ta1_x_yyyzz_yyz_0, ta1_x_yyyzz_yzz_0, ta1_x_yyyzz_zzz_0, ta1_x_yyzz_xxx_0, ta1_x_yyzz_xxx_1, ta1_x_yyzz_xxz_0, ta1_x_yyzz_xxz_1, ta1_x_yyzz_xyz_0, ta1_x_yyzz_xyz_1, ta1_x_yyzz_xz_0, ta1_x_yyzz_xz_1, ta1_x_yyzz_xzz_0, ta1_x_yyzz_xzz_1, ta1_x_yyzz_yyz_0, ta1_x_yyzz_yyz_1, ta1_x_yyzz_yz_0, ta1_x_yyzz_yz_1, ta1_x_yyzz_yzz_0, ta1_x_yyzz_yzz_1, ta1_x_yyzz_zz_0, ta1_x_yyzz_zz_1, ta1_x_yyzz_zzz_0, ta1_x_yyzz_zzz_1, ta1_x_yzz_xxx_0, ta1_x_yzz_xxx_1, ta1_x_yzz_xxz_0, ta1_x_yzz_xxz_1, ta1_x_yzz_xyz_0, ta1_x_yzz_xyz_1, ta1_x_yzz_xzz_0, ta1_x_yzz_xzz_1, ta1_x_yzz_yyz_0, ta1_x_yzz_yyz_1, ta1_x_yzz_yzz_0, ta1_x_yzz_yzz_1, ta1_x_yzz_zzz_0, ta1_x_yzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzz_xxx_0[i] = 2.0 * ta1_x_yzz_xxx_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxx_1[i] * fe_0 + ta1_x_yyzz_xxx_0[i] * pa_y[i] - ta1_x_yyzz_xxx_1[i] * pc_y[i];

        ta1_x_yyyzz_xxy_0[i] = ta1_x_yyy_xxy_0[i] * fe_0 - ta1_x_yyy_xxy_1[i] * fe_0 + ta1_x_yyyz_xxy_0[i] * pa_z[i] - ta1_x_yyyz_xxy_1[i] * pc_z[i];

        ta1_x_yyyzz_xxz_0[i] = 2.0 * ta1_x_yzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xxz_1[i] * fe_0 + ta1_x_yyzz_xxz_0[i] * pa_y[i] - ta1_x_yyzz_xxz_1[i] * pc_y[i];

        ta1_x_yyyzz_xyy_0[i] = ta1_x_yyy_xyy_0[i] * fe_0 - ta1_x_yyy_xyy_1[i] * fe_0 + ta1_x_yyyz_xyy_0[i] * pa_z[i] - ta1_x_yyyz_xyy_1[i] * pc_z[i];

        ta1_x_yyyzz_xyz_0[i] = 2.0 * ta1_x_yzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xyz_1[i] * fe_0 + ta1_x_yyzz_xz_0[i] * fe_0 - ta1_x_yyzz_xz_1[i] * fe_0 + ta1_x_yyzz_xyz_0[i] * pa_y[i] - ta1_x_yyzz_xyz_1[i] * pc_y[i];

        ta1_x_yyyzz_xzz_0[i] = 2.0 * ta1_x_yzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xzz_1[i] * fe_0 + ta1_x_yyzz_xzz_0[i] * pa_y[i] - ta1_x_yyzz_xzz_1[i] * pc_y[i];

        ta1_x_yyyzz_yyy_0[i] = ta1_x_yyy_yyy_0[i] * fe_0 - ta1_x_yyy_yyy_1[i] * fe_0 + ta1_x_yyyz_yyy_0[i] * pa_z[i] - ta1_x_yyyz_yyy_1[i] * pc_z[i];

        ta1_x_yyyzz_yyz_0[i] = 2.0 * ta1_x_yzz_yyz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyzz_yz_0[i] * fe_0 - 2.0 * ta1_x_yyzz_yz_1[i] * fe_0 + ta1_x_yyzz_yyz_0[i] * pa_y[i] - ta1_x_yyzz_yyz_1[i] * pc_y[i];

        ta1_x_yyyzz_yzz_0[i] = 2.0 * ta1_x_yzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yzz_1[i] * fe_0 + ta1_x_yyzz_zz_0[i] * fe_0 - ta1_x_yyzz_zz_1[i] * fe_0 + ta1_x_yyzz_yzz_0[i] * pa_y[i] - ta1_x_yyzz_yzz_1[i] * pc_y[i];

        ta1_x_yyyzz_zzz_0[i] = 2.0 * ta1_x_yzz_zzz_0[i] * fe_0 - 2.0 * ta1_x_yzz_zzz_1[i] * fe_0 + ta1_x_yyzz_zzz_0[i] * pa_y[i] - ta1_x_yyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 180-190 components of targeted buffer : HF

    auto ta1_x_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 180);

    auto ta1_x_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 181);

    auto ta1_x_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 182);

    auto ta1_x_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 183);

    auto ta1_x_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 184);

    auto ta1_x_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 185);

    auto ta1_x_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 186);

    auto ta1_x_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 187);

    auto ta1_x_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 188);

    auto ta1_x_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 189);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyz_xxy_0, ta1_x_yyz_xxy_1, ta1_x_yyz_xyy_0, ta1_x_yyz_xyy_1, ta1_x_yyz_yyy_0, ta1_x_yyz_yyy_1, ta1_x_yyzz_xxy_0, ta1_x_yyzz_xxy_1, ta1_x_yyzz_xyy_0, ta1_x_yyzz_xyy_1, ta1_x_yyzz_yyy_0, ta1_x_yyzz_yyy_1, ta1_x_yyzzz_xxx_0, ta1_x_yyzzz_xxy_0, ta1_x_yyzzz_xxz_0, ta1_x_yyzzz_xyy_0, ta1_x_yyzzz_xyz_0, ta1_x_yyzzz_xzz_0, ta1_x_yyzzz_yyy_0, ta1_x_yyzzz_yyz_0, ta1_x_yyzzz_yzz_0, ta1_x_yyzzz_zzz_0, ta1_x_yzzz_xxx_0, ta1_x_yzzz_xxx_1, ta1_x_yzzz_xxz_0, ta1_x_yzzz_xxz_1, ta1_x_yzzz_xyz_0, ta1_x_yzzz_xyz_1, ta1_x_yzzz_xz_0, ta1_x_yzzz_xz_1, ta1_x_yzzz_xzz_0, ta1_x_yzzz_xzz_1, ta1_x_yzzz_yyz_0, ta1_x_yzzz_yyz_1, ta1_x_yzzz_yz_0, ta1_x_yzzz_yz_1, ta1_x_yzzz_yzz_0, ta1_x_yzzz_yzz_1, ta1_x_yzzz_zz_0, ta1_x_yzzz_zz_1, ta1_x_yzzz_zzz_0, ta1_x_yzzz_zzz_1, ta1_x_zzz_xxx_0, ta1_x_zzz_xxx_1, ta1_x_zzz_xxz_0, ta1_x_zzz_xxz_1, ta1_x_zzz_xyz_0, ta1_x_zzz_xyz_1, ta1_x_zzz_xzz_0, ta1_x_zzz_xzz_1, ta1_x_zzz_yyz_0, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_0, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_0, ta1_x_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzz_xxx_0[i] = ta1_x_zzz_xxx_0[i] * fe_0 - ta1_x_zzz_xxx_1[i] * fe_0 + ta1_x_yzzz_xxx_0[i] * pa_y[i] - ta1_x_yzzz_xxx_1[i] * pc_y[i];

        ta1_x_yyzzz_xxy_0[i] = 2.0 * ta1_x_yyz_xxy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xxy_1[i] * fe_0 + ta1_x_yyzz_xxy_0[i] * pa_z[i] - ta1_x_yyzz_xxy_1[i] * pc_z[i];

        ta1_x_yyzzz_xxz_0[i] = ta1_x_zzz_xxz_0[i] * fe_0 - ta1_x_zzz_xxz_1[i] * fe_0 + ta1_x_yzzz_xxz_0[i] * pa_y[i] - ta1_x_yzzz_xxz_1[i] * pc_y[i];

        ta1_x_yyzzz_xyy_0[i] = 2.0 * ta1_x_yyz_xyy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xyy_1[i] * fe_0 + ta1_x_yyzz_xyy_0[i] * pa_z[i] - ta1_x_yyzz_xyy_1[i] * pc_z[i];

        ta1_x_yyzzz_xyz_0[i] = ta1_x_zzz_xyz_0[i] * fe_0 - ta1_x_zzz_xyz_1[i] * fe_0 + ta1_x_yzzz_xz_0[i] * fe_0 - ta1_x_yzzz_xz_1[i] * fe_0 + ta1_x_yzzz_xyz_0[i] * pa_y[i] - ta1_x_yzzz_xyz_1[i] * pc_y[i];

        ta1_x_yyzzz_xzz_0[i] = ta1_x_zzz_xzz_0[i] * fe_0 - ta1_x_zzz_xzz_1[i] * fe_0 + ta1_x_yzzz_xzz_0[i] * pa_y[i] - ta1_x_yzzz_xzz_1[i] * pc_y[i];

        ta1_x_yyzzz_yyy_0[i] = 2.0 * ta1_x_yyz_yyy_0[i] * fe_0 - 2.0 * ta1_x_yyz_yyy_1[i] * fe_0 + ta1_x_yyzz_yyy_0[i] * pa_z[i] - ta1_x_yyzz_yyy_1[i] * pc_z[i];

        ta1_x_yyzzz_yyz_0[i] = ta1_x_zzz_yyz_0[i] * fe_0 - ta1_x_zzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yz_1[i] * fe_0 + ta1_x_yzzz_yyz_0[i] * pa_y[i] - ta1_x_yzzz_yyz_1[i] * pc_y[i];

        ta1_x_yyzzz_yzz_0[i] = ta1_x_zzz_yzz_0[i] * fe_0 - ta1_x_zzz_yzz_1[i] * fe_0 + ta1_x_yzzz_zz_0[i] * fe_0 - ta1_x_yzzz_zz_1[i] * fe_0 + ta1_x_yzzz_yzz_0[i] * pa_y[i] - ta1_x_yzzz_yzz_1[i] * pc_y[i];

        ta1_x_yyzzz_zzz_0[i] = ta1_x_zzz_zzz_0[i] * fe_0 - ta1_x_zzz_zzz_1[i] * fe_0 + ta1_x_yzzz_zzz_0[i] * pa_y[i] - ta1_x_yzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 190-200 components of targeted buffer : HF

    auto ta1_x_yzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 190);

    auto ta1_x_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 191);

    auto ta1_x_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 192);

    auto ta1_x_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 193);

    auto ta1_x_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 194);

    auto ta1_x_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 195);

    auto ta1_x_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 196);

    auto ta1_x_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 197);

    auto ta1_x_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 198);

    auto ta1_x_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 199);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yzzzz_xxx_0, ta1_x_yzzzz_xxy_0, ta1_x_yzzzz_xxz_0, ta1_x_yzzzz_xyy_0, ta1_x_yzzzz_xyz_0, ta1_x_yzzzz_xzz_0, ta1_x_yzzzz_yyy_0, ta1_x_yzzzz_yyz_0, ta1_x_yzzzz_yzz_0, ta1_x_yzzzz_zzz_0, ta1_x_zzzz_xx_0, ta1_x_zzzz_xx_1, ta1_x_zzzz_xxx_0, ta1_x_zzzz_xxx_1, ta1_x_zzzz_xxy_0, ta1_x_zzzz_xxy_1, ta1_x_zzzz_xxz_0, ta1_x_zzzz_xxz_1, ta1_x_zzzz_xy_0, ta1_x_zzzz_xy_1, ta1_x_zzzz_xyy_0, ta1_x_zzzz_xyy_1, ta1_x_zzzz_xyz_0, ta1_x_zzzz_xyz_1, ta1_x_zzzz_xz_0, ta1_x_zzzz_xz_1, ta1_x_zzzz_xzz_0, ta1_x_zzzz_xzz_1, ta1_x_zzzz_yy_0, ta1_x_zzzz_yy_1, ta1_x_zzzz_yyy_0, ta1_x_zzzz_yyy_1, ta1_x_zzzz_yyz_0, ta1_x_zzzz_yyz_1, ta1_x_zzzz_yz_0, ta1_x_zzzz_yz_1, ta1_x_zzzz_yzz_0, ta1_x_zzzz_yzz_1, ta1_x_zzzz_zz_0, ta1_x_zzzz_zz_1, ta1_x_zzzz_zzz_0, ta1_x_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzz_xxx_0[i] = ta1_x_zzzz_xxx_0[i] * pa_y[i] - ta1_x_zzzz_xxx_1[i] * pc_y[i];

        ta1_x_yzzzz_xxy_0[i] = ta1_x_zzzz_xx_0[i] * fe_0 - ta1_x_zzzz_xx_1[i] * fe_0 + ta1_x_zzzz_xxy_0[i] * pa_y[i] - ta1_x_zzzz_xxy_1[i] * pc_y[i];

        ta1_x_yzzzz_xxz_0[i] = ta1_x_zzzz_xxz_0[i] * pa_y[i] - ta1_x_zzzz_xxz_1[i] * pc_y[i];

        ta1_x_yzzzz_xyy_0[i] = 2.0 * ta1_x_zzzz_xy_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xy_1[i] * fe_0 + ta1_x_zzzz_xyy_0[i] * pa_y[i] - ta1_x_zzzz_xyy_1[i] * pc_y[i];

        ta1_x_yzzzz_xyz_0[i] = ta1_x_zzzz_xz_0[i] * fe_0 - ta1_x_zzzz_xz_1[i] * fe_0 + ta1_x_zzzz_xyz_0[i] * pa_y[i] - ta1_x_zzzz_xyz_1[i] * pc_y[i];

        ta1_x_yzzzz_xzz_0[i] = ta1_x_zzzz_xzz_0[i] * pa_y[i] - ta1_x_zzzz_xzz_1[i] * pc_y[i];

        ta1_x_yzzzz_yyy_0[i] = 3.0 * ta1_x_zzzz_yy_0[i] * fe_0 - 3.0 * ta1_x_zzzz_yy_1[i] * fe_0 + ta1_x_zzzz_yyy_0[i] * pa_y[i] - ta1_x_zzzz_yyy_1[i] * pc_y[i];

        ta1_x_yzzzz_yyz_0[i] = 2.0 * ta1_x_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_yz_1[i] * fe_0 + ta1_x_zzzz_yyz_0[i] * pa_y[i] - ta1_x_zzzz_yyz_1[i] * pc_y[i];

        ta1_x_yzzzz_yzz_0[i] = ta1_x_zzzz_zz_0[i] * fe_0 - ta1_x_zzzz_zz_1[i] * fe_0 + ta1_x_zzzz_yzz_0[i] * pa_y[i] - ta1_x_zzzz_yzz_1[i] * pc_y[i];

        ta1_x_yzzzz_zzz_0[i] = ta1_x_zzzz_zzz_0[i] * pa_y[i] - ta1_x_zzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 200-210 components of targeted buffer : HF

    auto ta1_x_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 200);

    auto ta1_x_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 201);

    auto ta1_x_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 202);

    auto ta1_x_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 203);

    auto ta1_x_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 204);

    auto ta1_x_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 205);

    auto ta1_x_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 206);

    auto ta1_x_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 207);

    auto ta1_x_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 208);

    auto ta1_x_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 209);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zzz_xxx_0, ta1_x_zzz_xxx_1, ta1_x_zzz_xxy_0, ta1_x_zzz_xxy_1, ta1_x_zzz_xxz_0, ta1_x_zzz_xxz_1, ta1_x_zzz_xyy_0, ta1_x_zzz_xyy_1, ta1_x_zzz_xyz_0, ta1_x_zzz_xyz_1, ta1_x_zzz_xzz_0, ta1_x_zzz_xzz_1, ta1_x_zzz_yyy_0, ta1_x_zzz_yyy_1, ta1_x_zzz_yyz_0, ta1_x_zzz_yyz_1, ta1_x_zzz_yzz_0, ta1_x_zzz_yzz_1, ta1_x_zzz_zzz_0, ta1_x_zzz_zzz_1, ta1_x_zzzz_xx_0, ta1_x_zzzz_xx_1, ta1_x_zzzz_xxx_0, ta1_x_zzzz_xxx_1, ta1_x_zzzz_xxy_0, ta1_x_zzzz_xxy_1, ta1_x_zzzz_xxz_0, ta1_x_zzzz_xxz_1, ta1_x_zzzz_xy_0, ta1_x_zzzz_xy_1, ta1_x_zzzz_xyy_0, ta1_x_zzzz_xyy_1, ta1_x_zzzz_xyz_0, ta1_x_zzzz_xyz_1, ta1_x_zzzz_xz_0, ta1_x_zzzz_xz_1, ta1_x_zzzz_xzz_0, ta1_x_zzzz_xzz_1, ta1_x_zzzz_yy_0, ta1_x_zzzz_yy_1, ta1_x_zzzz_yyy_0, ta1_x_zzzz_yyy_1, ta1_x_zzzz_yyz_0, ta1_x_zzzz_yyz_1, ta1_x_zzzz_yz_0, ta1_x_zzzz_yz_1, ta1_x_zzzz_yzz_0, ta1_x_zzzz_yzz_1, ta1_x_zzzz_zz_0, ta1_x_zzzz_zz_1, ta1_x_zzzz_zzz_0, ta1_x_zzzz_zzz_1, ta1_x_zzzzz_xxx_0, ta1_x_zzzzz_xxy_0, ta1_x_zzzzz_xxz_0, ta1_x_zzzzz_xyy_0, ta1_x_zzzzz_xyz_0, ta1_x_zzzzz_xzz_0, ta1_x_zzzzz_yyy_0, ta1_x_zzzzz_yyz_0, ta1_x_zzzzz_yzz_0, ta1_x_zzzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzz_xxx_0[i] = 4.0 * ta1_x_zzz_xxx_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxx_1[i] * fe_0 + ta1_x_zzzz_xxx_0[i] * pa_z[i] - ta1_x_zzzz_xxx_1[i] * pc_z[i];

        ta1_x_zzzzz_xxy_0[i] = 4.0 * ta1_x_zzz_xxy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxy_1[i] * fe_0 + ta1_x_zzzz_xxy_0[i] * pa_z[i] - ta1_x_zzzz_xxy_1[i] * pc_z[i];

        ta1_x_zzzzz_xxz_0[i] = 4.0 * ta1_x_zzz_xxz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxz_1[i] * fe_0 + ta1_x_zzzz_xx_0[i] * fe_0 - ta1_x_zzzz_xx_1[i] * fe_0 + ta1_x_zzzz_xxz_0[i] * pa_z[i] - ta1_x_zzzz_xxz_1[i] * pc_z[i];

        ta1_x_zzzzz_xyy_0[i] = 4.0 * ta1_x_zzz_xyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyy_1[i] * fe_0 + ta1_x_zzzz_xyy_0[i] * pa_z[i] - ta1_x_zzzz_xyy_1[i] * pc_z[i];

        ta1_x_zzzzz_xyz_0[i] = 4.0 * ta1_x_zzz_xyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyz_1[i] * fe_0 + ta1_x_zzzz_xy_0[i] * fe_0 - ta1_x_zzzz_xy_1[i] * fe_0 + ta1_x_zzzz_xyz_0[i] * pa_z[i] - ta1_x_zzzz_xyz_1[i] * pc_z[i];

        ta1_x_zzzzz_xzz_0[i] = 4.0 * ta1_x_zzz_xzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_xz_1[i] * fe_0 + ta1_x_zzzz_xzz_0[i] * pa_z[i] - ta1_x_zzzz_xzz_1[i] * pc_z[i];

        ta1_x_zzzzz_yyy_0[i] = 4.0 * ta1_x_zzz_yyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyy_1[i] * fe_0 + ta1_x_zzzz_yyy_0[i] * pa_z[i] - ta1_x_zzzz_yyy_1[i] * pc_z[i];

        ta1_x_zzzzz_yyz_0[i] = 4.0 * ta1_x_zzz_yyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyz_1[i] * fe_0 + ta1_x_zzzz_yy_0[i] * fe_0 - ta1_x_zzzz_yy_1[i] * fe_0 + ta1_x_zzzz_yyz_0[i] * pa_z[i] - ta1_x_zzzz_yyz_1[i] * pc_z[i];

        ta1_x_zzzzz_yzz_0[i] = 4.0 * ta1_x_zzz_yzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_zzzz_yz_1[i] * fe_0 + ta1_x_zzzz_yzz_0[i] * pa_z[i] - ta1_x_zzzz_yzz_1[i] * pc_z[i];

        ta1_x_zzzzz_zzz_0[i] = 4.0 * ta1_x_zzz_zzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_zzz_1[i] * fe_0 + 3.0 * ta1_x_zzzz_zz_0[i] * fe_0 - 3.0 * ta1_x_zzzz_zz_1[i] * fe_0 + ta1_x_zzzz_zzz_0[i] * pa_z[i] - ta1_x_zzzz_zzz_1[i] * pc_z[i];
    }

    // Set up 210-220 components of targeted buffer : HF

    auto ta1_y_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 210);

    auto ta1_y_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 211);

    auto ta1_y_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 212);

    auto ta1_y_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 213);

    auto ta1_y_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 214);

    auto ta1_y_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 215);

    auto ta1_y_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 216);

    auto ta1_y_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 217);

    auto ta1_y_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 218);

    auto ta1_y_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 219);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xxx_xxx_0, ta1_y_xxx_xxx_1, ta1_y_xxx_xxy_0, ta1_y_xxx_xxy_1, ta1_y_xxx_xxz_0, ta1_y_xxx_xxz_1, ta1_y_xxx_xyy_0, ta1_y_xxx_xyy_1, ta1_y_xxx_xyz_0, ta1_y_xxx_xyz_1, ta1_y_xxx_xzz_0, ta1_y_xxx_xzz_1, ta1_y_xxx_yyy_0, ta1_y_xxx_yyy_1, ta1_y_xxx_yyz_0, ta1_y_xxx_yyz_1, ta1_y_xxx_yzz_0, ta1_y_xxx_yzz_1, ta1_y_xxx_zzz_0, ta1_y_xxx_zzz_1, ta1_y_xxxx_xx_0, ta1_y_xxxx_xx_1, ta1_y_xxxx_xxx_0, ta1_y_xxxx_xxx_1, ta1_y_xxxx_xxy_0, ta1_y_xxxx_xxy_1, ta1_y_xxxx_xxz_0, ta1_y_xxxx_xxz_1, ta1_y_xxxx_xy_0, ta1_y_xxxx_xy_1, ta1_y_xxxx_xyy_0, ta1_y_xxxx_xyy_1, ta1_y_xxxx_xyz_0, ta1_y_xxxx_xyz_1, ta1_y_xxxx_xz_0, ta1_y_xxxx_xz_1, ta1_y_xxxx_xzz_0, ta1_y_xxxx_xzz_1, ta1_y_xxxx_yy_0, ta1_y_xxxx_yy_1, ta1_y_xxxx_yyy_0, ta1_y_xxxx_yyy_1, ta1_y_xxxx_yyz_0, ta1_y_xxxx_yyz_1, ta1_y_xxxx_yz_0, ta1_y_xxxx_yz_1, ta1_y_xxxx_yzz_0, ta1_y_xxxx_yzz_1, ta1_y_xxxx_zz_0, ta1_y_xxxx_zz_1, ta1_y_xxxx_zzz_0, ta1_y_xxxx_zzz_1, ta1_y_xxxxx_xxx_0, ta1_y_xxxxx_xxy_0, ta1_y_xxxxx_xxz_0, ta1_y_xxxxx_xyy_0, ta1_y_xxxxx_xyz_0, ta1_y_xxxxx_xzz_0, ta1_y_xxxxx_yyy_0, ta1_y_xxxxx_yyz_0, ta1_y_xxxxx_yzz_0, ta1_y_xxxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxx_xxx_0[i] = 4.0 * ta1_y_xxx_xxx_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxx_1[i] * fe_0 + 3.0 * ta1_y_xxxx_xx_0[i] * fe_0 - 3.0 * ta1_y_xxxx_xx_1[i] * fe_0 + ta1_y_xxxx_xxx_0[i] * pa_x[i] - ta1_y_xxxx_xxx_1[i] * pc_x[i];

        ta1_y_xxxxx_xxy_0[i] = 4.0 * ta1_y_xxx_xxy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xy_1[i] * fe_0 + ta1_y_xxxx_xxy_0[i] * pa_x[i] - ta1_y_xxxx_xxy_1[i] * pc_x[i];

        ta1_y_xxxxx_xxz_0[i] = 4.0 * ta1_y_xxx_xxz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xz_1[i] * fe_0 + ta1_y_xxxx_xxz_0[i] * pa_x[i] - ta1_y_xxxx_xxz_1[i] * pc_x[i];

        ta1_y_xxxxx_xyy_0[i] = 4.0 * ta1_y_xxx_xyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyy_1[i] * fe_0 + ta1_y_xxxx_yy_0[i] * fe_0 - ta1_y_xxxx_yy_1[i] * fe_0 + ta1_y_xxxx_xyy_0[i] * pa_x[i] - ta1_y_xxxx_xyy_1[i] * pc_x[i];

        ta1_y_xxxxx_xyz_0[i] = 4.0 * ta1_y_xxx_xyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyz_1[i] * fe_0 + ta1_y_xxxx_yz_0[i] * fe_0 - ta1_y_xxxx_yz_1[i] * fe_0 + ta1_y_xxxx_xyz_0[i] * pa_x[i] - ta1_y_xxxx_xyz_1[i] * pc_x[i];

        ta1_y_xxxxx_xzz_0[i] = 4.0 * ta1_y_xxx_xzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xzz_1[i] * fe_0 + ta1_y_xxxx_zz_0[i] * fe_0 - ta1_y_xxxx_zz_1[i] * fe_0 + ta1_y_xxxx_xzz_0[i] * pa_x[i] - ta1_y_xxxx_xzz_1[i] * pc_x[i];

        ta1_y_xxxxx_yyy_0[i] = 4.0 * ta1_y_xxx_yyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_yyy_1[i] * fe_0 + ta1_y_xxxx_yyy_0[i] * pa_x[i] - ta1_y_xxxx_yyy_1[i] * pc_x[i];

        ta1_y_xxxxx_yyz_0[i] = 4.0 * ta1_y_xxx_yyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yyz_1[i] * fe_0 + ta1_y_xxxx_yyz_0[i] * pa_x[i] - ta1_y_xxxx_yyz_1[i] * pc_x[i];

        ta1_y_xxxxx_yzz_0[i] = 4.0 * ta1_y_xxx_yzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yzz_1[i] * fe_0 + ta1_y_xxxx_yzz_0[i] * pa_x[i] - ta1_y_xxxx_yzz_1[i] * pc_x[i];

        ta1_y_xxxxx_zzz_0[i] = 4.0 * ta1_y_xxx_zzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_zzz_1[i] * fe_0 + ta1_y_xxxx_zzz_0[i] * pa_x[i] - ta1_y_xxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 220-230 components of targeted buffer : HF

    auto ta1_y_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 220);

    auto ta1_y_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 221);

    auto ta1_y_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 222);

    auto ta1_y_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 223);

    auto ta1_y_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 224);

    auto ta1_y_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 225);

    auto ta1_y_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 226);

    auto ta1_y_xxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 227);

    auto ta1_y_xxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 228);

    auto ta1_y_xxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 229);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxx_xx_0, ta1_y_xxxx_xx_1, ta1_y_xxxx_xxx_0, ta1_y_xxxx_xxx_1, ta1_y_xxxx_xxy_0, ta1_y_xxxx_xxy_1, ta1_y_xxxx_xxz_0, ta1_y_xxxx_xxz_1, ta1_y_xxxx_xy_0, ta1_y_xxxx_xy_1, ta1_y_xxxx_xyy_0, ta1_y_xxxx_xyy_1, ta1_y_xxxx_xyz_0, ta1_y_xxxx_xyz_1, ta1_y_xxxx_xz_0, ta1_y_xxxx_xz_1, ta1_y_xxxx_xzz_0, ta1_y_xxxx_xzz_1, ta1_y_xxxx_zzz_0, ta1_y_xxxx_zzz_1, ta1_y_xxxxy_xxx_0, ta1_y_xxxxy_xxy_0, ta1_y_xxxxy_xxz_0, ta1_y_xxxxy_xyy_0, ta1_y_xxxxy_xyz_0, ta1_y_xxxxy_xzz_0, ta1_y_xxxxy_yyy_0, ta1_y_xxxxy_yyz_0, ta1_y_xxxxy_yzz_0, ta1_y_xxxxy_zzz_0, ta1_y_xxxy_yyy_0, ta1_y_xxxy_yyy_1, ta1_y_xxxy_yyz_0, ta1_y_xxxy_yyz_1, ta1_y_xxxy_yzz_0, ta1_y_xxxy_yzz_1, ta1_y_xxy_yyy_0, ta1_y_xxy_yyy_1, ta1_y_xxy_yyz_0, ta1_y_xxy_yyz_1, ta1_y_xxy_yzz_0, ta1_y_xxy_yzz_1, ta_xxxx_xxx_1, ta_xxxx_xxy_1, ta_xxxx_xxz_1, ta_xxxx_xyy_1, ta_xxxx_xyz_1, ta_xxxx_xzz_1, ta_xxxx_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxy_xxx_0[i] = ta_xxxx_xxx_1[i] + ta1_y_xxxx_xxx_0[i] * pa_y[i] - ta1_y_xxxx_xxx_1[i] * pc_y[i];

        ta1_y_xxxxy_xxy_0[i] = ta1_y_xxxx_xx_0[i] * fe_0 - ta1_y_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxy_1[i] + ta1_y_xxxx_xxy_0[i] * pa_y[i] - ta1_y_xxxx_xxy_1[i] * pc_y[i];

        ta1_y_xxxxy_xxz_0[i] = ta_xxxx_xxz_1[i] + ta1_y_xxxx_xxz_0[i] * pa_y[i] - ta1_y_xxxx_xxz_1[i] * pc_y[i];

        ta1_y_xxxxy_xyy_0[i] = 2.0 * ta1_y_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xy_1[i] * fe_0 + ta_xxxx_xyy_1[i] + ta1_y_xxxx_xyy_0[i] * pa_y[i] - ta1_y_xxxx_xyy_1[i] * pc_y[i];

        ta1_y_xxxxy_xyz_0[i] = ta1_y_xxxx_xz_0[i] * fe_0 - ta1_y_xxxx_xz_1[i] * fe_0 + ta_xxxx_xyz_1[i] + ta1_y_xxxx_xyz_0[i] * pa_y[i] - ta1_y_xxxx_xyz_1[i] * pc_y[i];

        ta1_y_xxxxy_xzz_0[i] = ta_xxxx_xzz_1[i] + ta1_y_xxxx_xzz_0[i] * pa_y[i] - ta1_y_xxxx_xzz_1[i] * pc_y[i];

        ta1_y_xxxxy_yyy_0[i] = 3.0 * ta1_y_xxy_yyy_0[i] * fe_0 - 3.0 * ta1_y_xxy_yyy_1[i] * fe_0 + ta1_y_xxxy_yyy_0[i] * pa_x[i] - ta1_y_xxxy_yyy_1[i] * pc_x[i];

        ta1_y_xxxxy_yyz_0[i] = 3.0 * ta1_y_xxy_yyz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yyz_1[i] * fe_0 + ta1_y_xxxy_yyz_0[i] * pa_x[i] - ta1_y_xxxy_yyz_1[i] * pc_x[i];

        ta1_y_xxxxy_yzz_0[i] = 3.0 * ta1_y_xxy_yzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yzz_1[i] * fe_0 + ta1_y_xxxy_yzz_0[i] * pa_x[i] - ta1_y_xxxy_yzz_1[i] * pc_x[i];

        ta1_y_xxxxy_zzz_0[i] = ta_xxxx_zzz_1[i] + ta1_y_xxxx_zzz_0[i] * pa_y[i] - ta1_y_xxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 230-240 components of targeted buffer : HF

    auto ta1_y_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 230);

    auto ta1_y_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 231);

    auto ta1_y_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 232);

    auto ta1_y_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 233);

    auto ta1_y_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 234);

    auto ta1_y_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 235);

    auto ta1_y_xxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 236);

    auto ta1_y_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 237);

    auto ta1_y_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 238);

    auto ta1_y_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 239);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxx_xx_0, ta1_y_xxxx_xx_1, ta1_y_xxxx_xxx_0, ta1_y_xxxx_xxx_1, ta1_y_xxxx_xxy_0, ta1_y_xxxx_xxy_1, ta1_y_xxxx_xxz_0, ta1_y_xxxx_xxz_1, ta1_y_xxxx_xy_0, ta1_y_xxxx_xy_1, ta1_y_xxxx_xyy_0, ta1_y_xxxx_xyy_1, ta1_y_xxxx_xyz_0, ta1_y_xxxx_xyz_1, ta1_y_xxxx_xz_0, ta1_y_xxxx_xz_1, ta1_y_xxxx_xzz_0, ta1_y_xxxx_xzz_1, ta1_y_xxxx_yyy_0, ta1_y_xxxx_yyy_1, ta1_y_xxxxz_xxx_0, ta1_y_xxxxz_xxy_0, ta1_y_xxxxz_xxz_0, ta1_y_xxxxz_xyy_0, ta1_y_xxxxz_xyz_0, ta1_y_xxxxz_xzz_0, ta1_y_xxxxz_yyy_0, ta1_y_xxxxz_yyz_0, ta1_y_xxxxz_yzz_0, ta1_y_xxxxz_zzz_0, ta1_y_xxxz_yyz_0, ta1_y_xxxz_yyz_1, ta1_y_xxxz_yzz_0, ta1_y_xxxz_yzz_1, ta1_y_xxxz_zzz_0, ta1_y_xxxz_zzz_1, ta1_y_xxz_yyz_0, ta1_y_xxz_yyz_1, ta1_y_xxz_yzz_0, ta1_y_xxz_yzz_1, ta1_y_xxz_zzz_0, ta1_y_xxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxz_xxx_0[i] = ta1_y_xxxx_xxx_0[i] * pa_z[i] - ta1_y_xxxx_xxx_1[i] * pc_z[i];

        ta1_y_xxxxz_xxy_0[i] = ta1_y_xxxx_xxy_0[i] * pa_z[i] - ta1_y_xxxx_xxy_1[i] * pc_z[i];

        ta1_y_xxxxz_xxz_0[i] = ta1_y_xxxx_xx_0[i] * fe_0 - ta1_y_xxxx_xx_1[i] * fe_0 + ta1_y_xxxx_xxz_0[i] * pa_z[i] - ta1_y_xxxx_xxz_1[i] * pc_z[i];

        ta1_y_xxxxz_xyy_0[i] = ta1_y_xxxx_xyy_0[i] * pa_z[i] - ta1_y_xxxx_xyy_1[i] * pc_z[i];

        ta1_y_xxxxz_xyz_0[i] = ta1_y_xxxx_xy_0[i] * fe_0 - ta1_y_xxxx_xy_1[i] * fe_0 + ta1_y_xxxx_xyz_0[i] * pa_z[i] - ta1_y_xxxx_xyz_1[i] * pc_z[i];

        ta1_y_xxxxz_xzz_0[i] = 2.0 * ta1_y_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_y_xxxx_xz_1[i] * fe_0 + ta1_y_xxxx_xzz_0[i] * pa_z[i] - ta1_y_xxxx_xzz_1[i] * pc_z[i];

        ta1_y_xxxxz_yyy_0[i] = ta1_y_xxxx_yyy_0[i] * pa_z[i] - ta1_y_xxxx_yyy_1[i] * pc_z[i];

        ta1_y_xxxxz_yyz_0[i] = 3.0 * ta1_y_xxz_yyz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yyz_1[i] * fe_0 + ta1_y_xxxz_yyz_0[i] * pa_x[i] - ta1_y_xxxz_yyz_1[i] * pc_x[i];

        ta1_y_xxxxz_yzz_0[i] = 3.0 * ta1_y_xxz_yzz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yzz_1[i] * fe_0 + ta1_y_xxxz_yzz_0[i] * pa_x[i] - ta1_y_xxxz_yzz_1[i] * pc_x[i];

        ta1_y_xxxxz_zzz_0[i] = 3.0 * ta1_y_xxz_zzz_0[i] * fe_0 - 3.0 * ta1_y_xxz_zzz_1[i] * fe_0 + ta1_y_xxxz_zzz_0[i] * pa_x[i] - ta1_y_xxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 240-250 components of targeted buffer : HF

    auto ta1_y_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 240);

    auto ta1_y_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 241);

    auto ta1_y_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 242);

    auto ta1_y_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 243);

    auto ta1_y_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 244);

    auto ta1_y_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 245);

    auto ta1_y_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 246);

    auto ta1_y_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 247);

    auto ta1_y_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 248);

    auto ta1_y_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 249);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxx_xxx_0, ta1_y_xxx_xxx_1, ta1_y_xxx_xxz_0, ta1_y_xxx_xxz_1, ta1_y_xxx_xzz_0, ta1_y_xxx_xzz_1, ta1_y_xxxy_xxx_0, ta1_y_xxxy_xxx_1, ta1_y_xxxy_xxz_0, ta1_y_xxxy_xxz_1, ta1_y_xxxy_xzz_0, ta1_y_xxxy_xzz_1, ta1_y_xxxyy_xxx_0, ta1_y_xxxyy_xxy_0, ta1_y_xxxyy_xxz_0, ta1_y_xxxyy_xyy_0, ta1_y_xxxyy_xyz_0, ta1_y_xxxyy_xzz_0, ta1_y_xxxyy_yyy_0, ta1_y_xxxyy_yyz_0, ta1_y_xxxyy_yzz_0, ta1_y_xxxyy_zzz_0, ta1_y_xxyy_xxy_0, ta1_y_xxyy_xxy_1, ta1_y_xxyy_xy_0, ta1_y_xxyy_xy_1, ta1_y_xxyy_xyy_0, ta1_y_xxyy_xyy_1, ta1_y_xxyy_xyz_0, ta1_y_xxyy_xyz_1, ta1_y_xxyy_yy_0, ta1_y_xxyy_yy_1, ta1_y_xxyy_yyy_0, ta1_y_xxyy_yyy_1, ta1_y_xxyy_yyz_0, ta1_y_xxyy_yyz_1, ta1_y_xxyy_yz_0, ta1_y_xxyy_yz_1, ta1_y_xxyy_yzz_0, ta1_y_xxyy_yzz_1, ta1_y_xxyy_zzz_0, ta1_y_xxyy_zzz_1, ta1_y_xyy_xxy_0, ta1_y_xyy_xxy_1, ta1_y_xyy_xyy_0, ta1_y_xyy_xyy_1, ta1_y_xyy_xyz_0, ta1_y_xyy_xyz_1, ta1_y_xyy_yyy_0, ta1_y_xyy_yyy_1, ta1_y_xyy_yyz_0, ta1_y_xyy_yyz_1, ta1_y_xyy_yzz_0, ta1_y_xyy_yzz_1, ta1_y_xyy_zzz_0, ta1_y_xyy_zzz_1, ta_xxxy_xxx_1, ta_xxxy_xxz_1, ta_xxxy_xzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyy_xxx_0[i] = ta1_y_xxx_xxx_0[i] * fe_0 - ta1_y_xxx_xxx_1[i] * fe_0 + ta_xxxy_xxx_1[i] + ta1_y_xxxy_xxx_0[i] * pa_y[i] - ta1_y_xxxy_xxx_1[i] * pc_y[i];

        ta1_y_xxxyy_xxy_0[i] = 2.0 * ta1_y_xyy_xxy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxyy_xy_0[i] * fe_0 - 2.0 * ta1_y_xxyy_xy_1[i] * fe_0 + ta1_y_xxyy_xxy_0[i] * pa_x[i] - ta1_y_xxyy_xxy_1[i] * pc_x[i];

        ta1_y_xxxyy_xxz_0[i] = ta1_y_xxx_xxz_0[i] * fe_0 - ta1_y_xxx_xxz_1[i] * fe_0 + ta_xxxy_xxz_1[i] + ta1_y_xxxy_xxz_0[i] * pa_y[i] - ta1_y_xxxy_xxz_1[i] * pc_y[i];

        ta1_y_xxxyy_xyy_0[i] = 2.0 * ta1_y_xyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xyy_1[i] * fe_0 + ta1_y_xxyy_yy_0[i] * fe_0 - ta1_y_xxyy_yy_1[i] * fe_0 + ta1_y_xxyy_xyy_0[i] * pa_x[i] - ta1_y_xxyy_xyy_1[i] * pc_x[i];

        ta1_y_xxxyy_xyz_0[i] = 2.0 * ta1_y_xyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xyy_xyz_1[i] * fe_0 + ta1_y_xxyy_yz_0[i] * fe_0 - ta1_y_xxyy_yz_1[i] * fe_0 + ta1_y_xxyy_xyz_0[i] * pa_x[i] - ta1_y_xxyy_xyz_1[i] * pc_x[i];

        ta1_y_xxxyy_xzz_0[i] = ta1_y_xxx_xzz_0[i] * fe_0 - ta1_y_xxx_xzz_1[i] * fe_0 + ta_xxxy_xzz_1[i] + ta1_y_xxxy_xzz_0[i] * pa_y[i] - ta1_y_xxxy_xzz_1[i] * pc_y[i];

        ta1_y_xxxyy_yyy_0[i] = 2.0 * ta1_y_xyy_yyy_0[i] * fe_0 - 2.0 * ta1_y_xyy_yyy_1[i] * fe_0 + ta1_y_xxyy_yyy_0[i] * pa_x[i] - ta1_y_xxyy_yyy_1[i] * pc_x[i];

        ta1_y_xxxyy_yyz_0[i] = 2.0 * ta1_y_xyy_yyz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yyz_1[i] * fe_0 + ta1_y_xxyy_yyz_0[i] * pa_x[i] - ta1_y_xxyy_yyz_1[i] * pc_x[i];

        ta1_y_xxxyy_yzz_0[i] = 2.0 * ta1_y_xyy_yzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yzz_1[i] * fe_0 + ta1_y_xxyy_yzz_0[i] * pa_x[i] - ta1_y_xxyy_yzz_1[i] * pc_x[i];

        ta1_y_xxxyy_zzz_0[i] = 2.0 * ta1_y_xyy_zzz_0[i] * fe_0 - 2.0 * ta1_y_xyy_zzz_1[i] * fe_0 + ta1_y_xxyy_zzz_0[i] * pa_x[i] - ta1_y_xxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 250-260 components of targeted buffer : HF

    auto ta1_y_xxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 250);

    auto ta1_y_xxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 251);

    auto ta1_y_xxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 252);

    auto ta1_y_xxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 253);

    auto ta1_y_xxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 254);

    auto ta1_y_xxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 255);

    auto ta1_y_xxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 256);

    auto ta1_y_xxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 257);

    auto ta1_y_xxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 258);

    auto ta1_y_xxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 259);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxxy_xxx_0, ta1_y_xxxy_xxx_1, ta1_y_xxxy_xxy_0, ta1_y_xxxy_xxy_1, ta1_y_xxxy_xy_0, ta1_y_xxxy_xy_1, ta1_y_xxxy_xyy_0, ta1_y_xxxy_xyy_1, ta1_y_xxxy_xyz_0, ta1_y_xxxy_xyz_1, ta1_y_xxxy_yyy_0, ta1_y_xxxy_yyy_1, ta1_y_xxxyz_xxx_0, ta1_y_xxxyz_xxy_0, ta1_y_xxxyz_xxz_0, ta1_y_xxxyz_xyy_0, ta1_y_xxxyz_xyz_0, ta1_y_xxxyz_xzz_0, ta1_y_xxxyz_yyy_0, ta1_y_xxxyz_yyz_0, ta1_y_xxxyz_yzz_0, ta1_y_xxxyz_zzz_0, ta1_y_xxxz_xxz_0, ta1_y_xxxz_xxz_1, ta1_y_xxxz_xzz_0, ta1_y_xxxz_xzz_1, ta1_y_xxxz_zzz_0, ta1_y_xxxz_zzz_1, ta1_y_xxyz_yyz_0, ta1_y_xxyz_yyz_1, ta1_y_xxyz_yzz_0, ta1_y_xxyz_yzz_1, ta1_y_xyz_yyz_0, ta1_y_xyz_yyz_1, ta1_y_xyz_yzz_0, ta1_y_xyz_yzz_1, ta_xxxz_xxz_1, ta_xxxz_xzz_1, ta_xxxz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyz_xxx_0[i] = ta1_y_xxxy_xxx_0[i] * pa_z[i] - ta1_y_xxxy_xxx_1[i] * pc_z[i];

        ta1_y_xxxyz_xxy_0[i] = ta1_y_xxxy_xxy_0[i] * pa_z[i] - ta1_y_xxxy_xxy_1[i] * pc_z[i];

        ta1_y_xxxyz_xxz_0[i] = ta_xxxz_xxz_1[i] + ta1_y_xxxz_xxz_0[i] * pa_y[i] - ta1_y_xxxz_xxz_1[i] * pc_y[i];

        ta1_y_xxxyz_xyy_0[i] = ta1_y_xxxy_xyy_0[i] * pa_z[i] - ta1_y_xxxy_xyy_1[i] * pc_z[i];

        ta1_y_xxxyz_xyz_0[i] = ta1_y_xxxy_xy_0[i] * fe_0 - ta1_y_xxxy_xy_1[i] * fe_0 + ta1_y_xxxy_xyz_0[i] * pa_z[i] - ta1_y_xxxy_xyz_1[i] * pc_z[i];

        ta1_y_xxxyz_xzz_0[i] = ta_xxxz_xzz_1[i] + ta1_y_xxxz_xzz_0[i] * pa_y[i] - ta1_y_xxxz_xzz_1[i] * pc_y[i];

        ta1_y_xxxyz_yyy_0[i] = ta1_y_xxxy_yyy_0[i] * pa_z[i] - ta1_y_xxxy_yyy_1[i] * pc_z[i];

        ta1_y_xxxyz_yyz_0[i] = 2.0 * ta1_y_xyz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yyz_1[i] * fe_0 + ta1_y_xxyz_yyz_0[i] * pa_x[i] - ta1_y_xxyz_yyz_1[i] * pc_x[i];

        ta1_y_xxxyz_yzz_0[i] = 2.0 * ta1_y_xyz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yzz_1[i] * fe_0 + ta1_y_xxyz_yzz_0[i] * pa_x[i] - ta1_y_xxyz_yzz_1[i] * pc_x[i];

        ta1_y_xxxyz_zzz_0[i] = ta_xxxz_zzz_1[i] + ta1_y_xxxz_zzz_0[i] * pa_y[i] - ta1_y_xxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 260-270 components of targeted buffer : HF

    auto ta1_y_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 260);

    auto ta1_y_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 261);

    auto ta1_y_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 262);

    auto ta1_y_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 263);

    auto ta1_y_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 264);

    auto ta1_y_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 265);

    auto ta1_y_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 266);

    auto ta1_y_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 267);

    auto ta1_y_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 268);

    auto ta1_y_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 269);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxx_xxx_0, ta1_y_xxx_xxx_1, ta1_y_xxx_xxy_0, ta1_y_xxx_xxy_1, ta1_y_xxx_xyy_0, ta1_y_xxx_xyy_1, ta1_y_xxxz_xxx_0, ta1_y_xxxz_xxx_1, ta1_y_xxxz_xxy_0, ta1_y_xxxz_xxy_1, ta1_y_xxxz_xyy_0, ta1_y_xxxz_xyy_1, ta1_y_xxxzz_xxx_0, ta1_y_xxxzz_xxy_0, ta1_y_xxxzz_xxz_0, ta1_y_xxxzz_xyy_0, ta1_y_xxxzz_xyz_0, ta1_y_xxxzz_xzz_0, ta1_y_xxxzz_yyy_0, ta1_y_xxxzz_yyz_0, ta1_y_xxxzz_yzz_0, ta1_y_xxxzz_zzz_0, ta1_y_xxzz_xxz_0, ta1_y_xxzz_xxz_1, ta1_y_xxzz_xyz_0, ta1_y_xxzz_xyz_1, ta1_y_xxzz_xz_0, ta1_y_xxzz_xz_1, ta1_y_xxzz_xzz_0, ta1_y_xxzz_xzz_1, ta1_y_xxzz_yyy_0, ta1_y_xxzz_yyy_1, ta1_y_xxzz_yyz_0, ta1_y_xxzz_yyz_1, ta1_y_xxzz_yz_0, ta1_y_xxzz_yz_1, ta1_y_xxzz_yzz_0, ta1_y_xxzz_yzz_1, ta1_y_xxzz_zz_0, ta1_y_xxzz_zz_1, ta1_y_xxzz_zzz_0, ta1_y_xxzz_zzz_1, ta1_y_xzz_xxz_0, ta1_y_xzz_xxz_1, ta1_y_xzz_xyz_0, ta1_y_xzz_xyz_1, ta1_y_xzz_xzz_0, ta1_y_xzz_xzz_1, ta1_y_xzz_yyy_0, ta1_y_xzz_yyy_1, ta1_y_xzz_yyz_0, ta1_y_xzz_yyz_1, ta1_y_xzz_yzz_0, ta1_y_xzz_yzz_1, ta1_y_xzz_zzz_0, ta1_y_xzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzz_xxx_0[i] = ta1_y_xxx_xxx_0[i] * fe_0 - ta1_y_xxx_xxx_1[i] * fe_0 + ta1_y_xxxz_xxx_0[i] * pa_z[i] - ta1_y_xxxz_xxx_1[i] * pc_z[i];

        ta1_y_xxxzz_xxy_0[i] = ta1_y_xxx_xxy_0[i] * fe_0 - ta1_y_xxx_xxy_1[i] * fe_0 + ta1_y_xxxz_xxy_0[i] * pa_z[i] - ta1_y_xxxz_xxy_1[i] * pc_z[i];

        ta1_y_xxxzz_xxz_0[i] = 2.0 * ta1_y_xzz_xxz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxzz_xz_0[i] * fe_0 - 2.0 * ta1_y_xxzz_xz_1[i] * fe_0 + ta1_y_xxzz_xxz_0[i] * pa_x[i] - ta1_y_xxzz_xxz_1[i] * pc_x[i];

        ta1_y_xxxzz_xyy_0[i] = ta1_y_xxx_xyy_0[i] * fe_0 - ta1_y_xxx_xyy_1[i] * fe_0 + ta1_y_xxxz_xyy_0[i] * pa_z[i] - ta1_y_xxxz_xyy_1[i] * pc_z[i];

        ta1_y_xxxzz_xyz_0[i] = 2.0 * ta1_y_xzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xyz_1[i] * fe_0 + ta1_y_xxzz_yz_0[i] * fe_0 - ta1_y_xxzz_yz_1[i] * fe_0 + ta1_y_xxzz_xyz_0[i] * pa_x[i] - ta1_y_xxzz_xyz_1[i] * pc_x[i];

        ta1_y_xxxzz_xzz_0[i] = 2.0 * ta1_y_xzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xzz_1[i] * fe_0 + ta1_y_xxzz_zz_0[i] * fe_0 - ta1_y_xxzz_zz_1[i] * fe_0 + ta1_y_xxzz_xzz_0[i] * pa_x[i] - ta1_y_xxzz_xzz_1[i] * pc_x[i];

        ta1_y_xxxzz_yyy_0[i] = 2.0 * ta1_y_xzz_yyy_0[i] * fe_0 - 2.0 * ta1_y_xzz_yyy_1[i] * fe_0 + ta1_y_xxzz_yyy_0[i] * pa_x[i] - ta1_y_xxzz_yyy_1[i] * pc_x[i];

        ta1_y_xxxzz_yyz_0[i] = 2.0 * ta1_y_xzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yyz_1[i] * fe_0 + ta1_y_xxzz_yyz_0[i] * pa_x[i] - ta1_y_xxzz_yyz_1[i] * pc_x[i];

        ta1_y_xxxzz_yzz_0[i] = 2.0 * ta1_y_xzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yzz_1[i] * fe_0 + ta1_y_xxzz_yzz_0[i] * pa_x[i] - ta1_y_xxzz_yzz_1[i] * pc_x[i];

        ta1_y_xxxzz_zzz_0[i] = 2.0 * ta1_y_xzz_zzz_0[i] * fe_0 - 2.0 * ta1_y_xzz_zzz_1[i] * fe_0 + ta1_y_xxzz_zzz_0[i] * pa_x[i] - ta1_y_xxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 270-280 components of targeted buffer : HF

    auto ta1_y_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 270);

    auto ta1_y_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 271);

    auto ta1_y_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 272);

    auto ta1_y_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 273);

    auto ta1_y_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 274);

    auto ta1_y_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 275);

    auto ta1_y_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 276);

    auto ta1_y_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 277);

    auto ta1_y_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 278);

    auto ta1_y_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 279);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxy_xxx_0, ta1_y_xxy_xxx_1, ta1_y_xxy_xxz_0, ta1_y_xxy_xxz_1, ta1_y_xxy_xzz_0, ta1_y_xxy_xzz_1, ta1_y_xxyy_xxx_0, ta1_y_xxyy_xxx_1, ta1_y_xxyy_xxz_0, ta1_y_xxyy_xxz_1, ta1_y_xxyy_xzz_0, ta1_y_xxyy_xzz_1, ta1_y_xxyyy_xxx_0, ta1_y_xxyyy_xxy_0, ta1_y_xxyyy_xxz_0, ta1_y_xxyyy_xyy_0, ta1_y_xxyyy_xyz_0, ta1_y_xxyyy_xzz_0, ta1_y_xxyyy_yyy_0, ta1_y_xxyyy_yyz_0, ta1_y_xxyyy_yzz_0, ta1_y_xxyyy_zzz_0, ta1_y_xyyy_xxy_0, ta1_y_xyyy_xxy_1, ta1_y_xyyy_xy_0, ta1_y_xyyy_xy_1, ta1_y_xyyy_xyy_0, ta1_y_xyyy_xyy_1, ta1_y_xyyy_xyz_0, ta1_y_xyyy_xyz_1, ta1_y_xyyy_yy_0, ta1_y_xyyy_yy_1, ta1_y_xyyy_yyy_0, ta1_y_xyyy_yyy_1, ta1_y_xyyy_yyz_0, ta1_y_xyyy_yyz_1, ta1_y_xyyy_yz_0, ta1_y_xyyy_yz_1, ta1_y_xyyy_yzz_0, ta1_y_xyyy_yzz_1, ta1_y_xyyy_zzz_0, ta1_y_xyyy_zzz_1, ta1_y_yyy_xxy_0, ta1_y_yyy_xxy_1, ta1_y_yyy_xyy_0, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_0, ta1_y_yyy_xyz_1, ta1_y_yyy_yyy_0, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_0, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_0, ta1_y_yyy_yzz_1, ta1_y_yyy_zzz_0, ta1_y_yyy_zzz_1, ta_xxyy_xxx_1, ta_xxyy_xxz_1, ta_xxyy_xzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyy_xxx_0[i] = 2.0 * ta1_y_xxy_xxx_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxx_1[i] * fe_0 + ta_xxyy_xxx_1[i] + ta1_y_xxyy_xxx_0[i] * pa_y[i] - ta1_y_xxyy_xxx_1[i] * pc_y[i];

        ta1_y_xxyyy_xxy_0[i] = ta1_y_yyy_xxy_0[i] * fe_0 - ta1_y_yyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xyyy_xy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xy_1[i] * fe_0 + ta1_y_xyyy_xxy_0[i] * pa_x[i] - ta1_y_xyyy_xxy_1[i] * pc_x[i];

        ta1_y_xxyyy_xxz_0[i] = 2.0 * ta1_y_xxy_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxz_1[i] * fe_0 + ta_xxyy_xxz_1[i] + ta1_y_xxyy_xxz_0[i] * pa_y[i] - ta1_y_xxyy_xxz_1[i] * pc_y[i];

        ta1_y_xxyyy_xyy_0[i] = ta1_y_yyy_xyy_0[i] * fe_0 - ta1_y_yyy_xyy_1[i] * fe_0 + ta1_y_xyyy_yy_0[i] * fe_0 - ta1_y_xyyy_yy_1[i] * fe_0 + ta1_y_xyyy_xyy_0[i] * pa_x[i] - ta1_y_xyyy_xyy_1[i] * pc_x[i];

        ta1_y_xxyyy_xyz_0[i] = ta1_y_yyy_xyz_0[i] * fe_0 - ta1_y_yyy_xyz_1[i] * fe_0 + ta1_y_xyyy_yz_0[i] * fe_0 - ta1_y_xyyy_yz_1[i] * fe_0 + ta1_y_xyyy_xyz_0[i] * pa_x[i] - ta1_y_xyyy_xyz_1[i] * pc_x[i];

        ta1_y_xxyyy_xzz_0[i] = 2.0 * ta1_y_xxy_xzz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xzz_1[i] * fe_0 + ta_xxyy_xzz_1[i] + ta1_y_xxyy_xzz_0[i] * pa_y[i] - ta1_y_xxyy_xzz_1[i] * pc_y[i];

        ta1_y_xxyyy_yyy_0[i] = ta1_y_yyy_yyy_0[i] * fe_0 - ta1_y_yyy_yyy_1[i] * fe_0 + ta1_y_xyyy_yyy_0[i] * pa_x[i] - ta1_y_xyyy_yyy_1[i] * pc_x[i];

        ta1_y_xxyyy_yyz_0[i] = ta1_y_yyy_yyz_0[i] * fe_0 - ta1_y_yyy_yyz_1[i] * fe_0 + ta1_y_xyyy_yyz_0[i] * pa_x[i] - ta1_y_xyyy_yyz_1[i] * pc_x[i];

        ta1_y_xxyyy_yzz_0[i] = ta1_y_yyy_yzz_0[i] * fe_0 - ta1_y_yyy_yzz_1[i] * fe_0 + ta1_y_xyyy_yzz_0[i] * pa_x[i] - ta1_y_xyyy_yzz_1[i] * pc_x[i];

        ta1_y_xxyyy_zzz_0[i] = ta1_y_yyy_zzz_0[i] * fe_0 - ta1_y_yyy_zzz_1[i] * fe_0 + ta1_y_xyyy_zzz_0[i] * pa_x[i] - ta1_y_xyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 280-290 components of targeted buffer : HF

    auto ta1_y_xxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 280);

    auto ta1_y_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 281);

    auto ta1_y_xxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 282);

    auto ta1_y_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 283);

    auto ta1_y_xxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 284);

    auto ta1_y_xxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 285);

    auto ta1_y_xxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 286);

    auto ta1_y_xxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 287);

    auto ta1_y_xxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 288);

    auto ta1_y_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 289);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxyy_xx_0, ta1_y_xxyy_xx_1, ta1_y_xxyy_xxx_0, ta1_y_xxyy_xxx_1, ta1_y_xxyy_xxy_0, ta1_y_xxyy_xxy_1, ta1_y_xxyy_xxz_0, ta1_y_xxyy_xxz_1, ta1_y_xxyy_xy_0, ta1_y_xxyy_xy_1, ta1_y_xxyy_xyy_0, ta1_y_xxyy_xyy_1, ta1_y_xxyy_xyz_0, ta1_y_xxyy_xyz_1, ta1_y_xxyy_xz_0, ta1_y_xxyy_xz_1, ta1_y_xxyy_xzz_0, ta1_y_xxyy_xzz_1, ta1_y_xxyy_yyy_0, ta1_y_xxyy_yyy_1, ta1_y_xxyyz_xxx_0, ta1_y_xxyyz_xxy_0, ta1_y_xxyyz_xxz_0, ta1_y_xxyyz_xyy_0, ta1_y_xxyyz_xyz_0, ta1_y_xxyyz_xzz_0, ta1_y_xxyyz_yyy_0, ta1_y_xxyyz_yyz_0, ta1_y_xxyyz_yzz_0, ta1_y_xxyyz_zzz_0, ta1_y_xyyz_yyz_0, ta1_y_xyyz_yyz_1, ta1_y_xyyz_yzz_0, ta1_y_xyyz_yzz_1, ta1_y_xyyz_zzz_0, ta1_y_xyyz_zzz_1, ta1_y_yyz_yyz_0, ta1_y_yyz_yyz_1, ta1_y_yyz_yzz_0, ta1_y_yyz_yzz_1, ta1_y_yyz_zzz_0, ta1_y_yyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyz_xxx_0[i] = ta1_y_xxyy_xxx_0[i] * pa_z[i] - ta1_y_xxyy_xxx_1[i] * pc_z[i];

        ta1_y_xxyyz_xxy_0[i] = ta1_y_xxyy_xxy_0[i] * pa_z[i] - ta1_y_xxyy_xxy_1[i] * pc_z[i];

        ta1_y_xxyyz_xxz_0[i] = ta1_y_xxyy_xx_0[i] * fe_0 - ta1_y_xxyy_xx_1[i] * fe_0 + ta1_y_xxyy_xxz_0[i] * pa_z[i] - ta1_y_xxyy_xxz_1[i] * pc_z[i];

        ta1_y_xxyyz_xyy_0[i] = ta1_y_xxyy_xyy_0[i] * pa_z[i] - ta1_y_xxyy_xyy_1[i] * pc_z[i];

        ta1_y_xxyyz_xyz_0[i] = ta1_y_xxyy_xy_0[i] * fe_0 - ta1_y_xxyy_xy_1[i] * fe_0 + ta1_y_xxyy_xyz_0[i] * pa_z[i] - ta1_y_xxyy_xyz_1[i] * pc_z[i];

        ta1_y_xxyyz_xzz_0[i] = 2.0 * ta1_y_xxyy_xz_0[i] * fe_0 - 2.0 * ta1_y_xxyy_xz_1[i] * fe_0 + ta1_y_xxyy_xzz_0[i] * pa_z[i] - ta1_y_xxyy_xzz_1[i] * pc_z[i];

        ta1_y_xxyyz_yyy_0[i] = ta1_y_xxyy_yyy_0[i] * pa_z[i] - ta1_y_xxyy_yyy_1[i] * pc_z[i];

        ta1_y_xxyyz_yyz_0[i] = ta1_y_yyz_yyz_0[i] * fe_0 - ta1_y_yyz_yyz_1[i] * fe_0 + ta1_y_xyyz_yyz_0[i] * pa_x[i] - ta1_y_xyyz_yyz_1[i] * pc_x[i];

        ta1_y_xxyyz_yzz_0[i] = ta1_y_yyz_yzz_0[i] * fe_0 - ta1_y_yyz_yzz_1[i] * fe_0 + ta1_y_xyyz_yzz_0[i] * pa_x[i] - ta1_y_xyyz_yzz_1[i] * pc_x[i];

        ta1_y_xxyyz_zzz_0[i] = ta1_y_yyz_zzz_0[i] * fe_0 - ta1_y_yyz_zzz_1[i] * fe_0 + ta1_y_xyyz_zzz_0[i] * pa_x[i] - ta1_y_xyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 290-300 components of targeted buffer : HF

    auto ta1_y_xxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 290);

    auto ta1_y_xxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 291);

    auto ta1_y_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 292);

    auto ta1_y_xxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 293);

    auto ta1_y_xxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 294);

    auto ta1_y_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 295);

    auto ta1_y_xxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 296);

    auto ta1_y_xxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 297);

    auto ta1_y_xxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 298);

    auto ta1_y_xxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 299);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xxy_xxy_0, ta1_y_xxy_xxy_1, ta1_y_xxy_xyy_0, ta1_y_xxy_xyy_1, ta1_y_xxyz_xxy_0, ta1_y_xxyz_xxy_1, ta1_y_xxyz_xyy_0, ta1_y_xxyz_xyy_1, ta1_y_xxyzz_xxx_0, ta1_y_xxyzz_xxy_0, ta1_y_xxyzz_xxz_0, ta1_y_xxyzz_xyy_0, ta1_y_xxyzz_xyz_0, ta1_y_xxyzz_xzz_0, ta1_y_xxyzz_yyy_0, ta1_y_xxyzz_yyz_0, ta1_y_xxyzz_yzz_0, ta1_y_xxyzz_zzz_0, ta1_y_xxzz_xxx_0, ta1_y_xxzz_xxx_1, ta1_y_xxzz_xxz_0, ta1_y_xxzz_xxz_1, ta1_y_xxzz_xyz_0, ta1_y_xxzz_xyz_1, ta1_y_xxzz_xz_0, ta1_y_xxzz_xz_1, ta1_y_xxzz_xzz_0, ta1_y_xxzz_xzz_1, ta1_y_xxzz_zzz_0, ta1_y_xxzz_zzz_1, ta1_y_xyzz_yyy_0, ta1_y_xyzz_yyy_1, ta1_y_xyzz_yyz_0, ta1_y_xyzz_yyz_1, ta1_y_xyzz_yzz_0, ta1_y_xyzz_yzz_1, ta1_y_yzz_yyy_0, ta1_y_yzz_yyy_1, ta1_y_yzz_yyz_0, ta1_y_yzz_yyz_1, ta1_y_yzz_yzz_0, ta1_y_yzz_yzz_1, ta_xxzz_xxx_1, ta_xxzz_xxz_1, ta_xxzz_xyz_1, ta_xxzz_xzz_1, ta_xxzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzz_xxx_0[i] = ta_xxzz_xxx_1[i] + ta1_y_xxzz_xxx_0[i] * pa_y[i] - ta1_y_xxzz_xxx_1[i] * pc_y[i];

        ta1_y_xxyzz_xxy_0[i] = ta1_y_xxy_xxy_0[i] * fe_0 - ta1_y_xxy_xxy_1[i] * fe_0 + ta1_y_xxyz_xxy_0[i] * pa_z[i] - ta1_y_xxyz_xxy_1[i] * pc_z[i];

        ta1_y_xxyzz_xxz_0[i] = ta_xxzz_xxz_1[i] + ta1_y_xxzz_xxz_0[i] * pa_y[i] - ta1_y_xxzz_xxz_1[i] * pc_y[i];

        ta1_y_xxyzz_xyy_0[i] = ta1_y_xxy_xyy_0[i] * fe_0 - ta1_y_xxy_xyy_1[i] * fe_0 + ta1_y_xxyz_xyy_0[i] * pa_z[i] - ta1_y_xxyz_xyy_1[i] * pc_z[i];

        ta1_y_xxyzz_xyz_0[i] = ta1_y_xxzz_xz_0[i] * fe_0 - ta1_y_xxzz_xz_1[i] * fe_0 + ta_xxzz_xyz_1[i] + ta1_y_xxzz_xyz_0[i] * pa_y[i] - ta1_y_xxzz_xyz_1[i] * pc_y[i];

        ta1_y_xxyzz_xzz_0[i] = ta_xxzz_xzz_1[i] + ta1_y_xxzz_xzz_0[i] * pa_y[i] - ta1_y_xxzz_xzz_1[i] * pc_y[i];

        ta1_y_xxyzz_yyy_0[i] = ta1_y_yzz_yyy_0[i] * fe_0 - ta1_y_yzz_yyy_1[i] * fe_0 + ta1_y_xyzz_yyy_0[i] * pa_x[i] - ta1_y_xyzz_yyy_1[i] * pc_x[i];

        ta1_y_xxyzz_yyz_0[i] = ta1_y_yzz_yyz_0[i] * fe_0 - ta1_y_yzz_yyz_1[i] * fe_0 + ta1_y_xyzz_yyz_0[i] * pa_x[i] - ta1_y_xyzz_yyz_1[i] * pc_x[i];

        ta1_y_xxyzz_yzz_0[i] = ta1_y_yzz_yzz_0[i] * fe_0 - ta1_y_yzz_yzz_1[i] * fe_0 + ta1_y_xyzz_yzz_0[i] * pa_x[i] - ta1_y_xyzz_yzz_1[i] * pc_x[i];

        ta1_y_xxyzz_zzz_0[i] = ta_xxzz_zzz_1[i] + ta1_y_xxzz_zzz_0[i] * pa_y[i] - ta1_y_xxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 300-310 components of targeted buffer : HF

    auto ta1_y_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 300);

    auto ta1_y_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 301);

    auto ta1_y_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 302);

    auto ta1_y_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 303);

    auto ta1_y_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 304);

    auto ta1_y_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 305);

    auto ta1_y_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 306);

    auto ta1_y_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 307);

    auto ta1_y_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 308);

    auto ta1_y_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 309);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxz_xxx_0, ta1_y_xxz_xxx_1, ta1_y_xxz_xxy_0, ta1_y_xxz_xxy_1, ta1_y_xxz_xyy_0, ta1_y_xxz_xyy_1, ta1_y_xxzz_xxx_0, ta1_y_xxzz_xxx_1, ta1_y_xxzz_xxy_0, ta1_y_xxzz_xxy_1, ta1_y_xxzz_xyy_0, ta1_y_xxzz_xyy_1, ta1_y_xxzzz_xxx_0, ta1_y_xxzzz_xxy_0, ta1_y_xxzzz_xxz_0, ta1_y_xxzzz_xyy_0, ta1_y_xxzzz_xyz_0, ta1_y_xxzzz_xzz_0, ta1_y_xxzzz_yyy_0, ta1_y_xxzzz_yyz_0, ta1_y_xxzzz_yzz_0, ta1_y_xxzzz_zzz_0, ta1_y_xzzz_xxz_0, ta1_y_xzzz_xxz_1, ta1_y_xzzz_xyz_0, ta1_y_xzzz_xyz_1, ta1_y_xzzz_xz_0, ta1_y_xzzz_xz_1, ta1_y_xzzz_xzz_0, ta1_y_xzzz_xzz_1, ta1_y_xzzz_yyy_0, ta1_y_xzzz_yyy_1, ta1_y_xzzz_yyz_0, ta1_y_xzzz_yyz_1, ta1_y_xzzz_yz_0, ta1_y_xzzz_yz_1, ta1_y_xzzz_yzz_0, ta1_y_xzzz_yzz_1, ta1_y_xzzz_zz_0, ta1_y_xzzz_zz_1, ta1_y_xzzz_zzz_0, ta1_y_xzzz_zzz_1, ta1_y_zzz_xxz_0, ta1_y_zzz_xxz_1, ta1_y_zzz_xyz_0, ta1_y_zzz_xyz_1, ta1_y_zzz_xzz_0, ta1_y_zzz_xzz_1, ta1_y_zzz_yyy_0, ta1_y_zzz_yyy_1, ta1_y_zzz_yyz_0, ta1_y_zzz_yyz_1, ta1_y_zzz_yzz_0, ta1_y_zzz_yzz_1, ta1_y_zzz_zzz_0, ta1_y_zzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzz_xxx_0[i] = 2.0 * ta1_y_xxz_xxx_0[i] * fe_0 - 2.0 * ta1_y_xxz_xxx_1[i] * fe_0 + ta1_y_xxzz_xxx_0[i] * pa_z[i] - ta1_y_xxzz_xxx_1[i] * pc_z[i];

        ta1_y_xxzzz_xxy_0[i] = 2.0 * ta1_y_xxz_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xxy_1[i] * fe_0 + ta1_y_xxzz_xxy_0[i] * pa_z[i] - ta1_y_xxzz_xxy_1[i] * pc_z[i];

        ta1_y_xxzzz_xxz_0[i] = ta1_y_zzz_xxz_0[i] * fe_0 - ta1_y_zzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xz_1[i] * fe_0 + ta1_y_xzzz_xxz_0[i] * pa_x[i] - ta1_y_xzzz_xxz_1[i] * pc_x[i];

        ta1_y_xxzzz_xyy_0[i] = 2.0 * ta1_y_xxz_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xyy_1[i] * fe_0 + ta1_y_xxzz_xyy_0[i] * pa_z[i] - ta1_y_xxzz_xyy_1[i] * pc_z[i];

        ta1_y_xxzzz_xyz_0[i] = ta1_y_zzz_xyz_0[i] * fe_0 - ta1_y_zzz_xyz_1[i] * fe_0 + ta1_y_xzzz_yz_0[i] * fe_0 - ta1_y_xzzz_yz_1[i] * fe_0 + ta1_y_xzzz_xyz_0[i] * pa_x[i] - ta1_y_xzzz_xyz_1[i] * pc_x[i];

        ta1_y_xxzzz_xzz_0[i] = ta1_y_zzz_xzz_0[i] * fe_0 - ta1_y_zzz_xzz_1[i] * fe_0 + ta1_y_xzzz_zz_0[i] * fe_0 - ta1_y_xzzz_zz_1[i] * fe_0 + ta1_y_xzzz_xzz_0[i] * pa_x[i] - ta1_y_xzzz_xzz_1[i] * pc_x[i];

        ta1_y_xxzzz_yyy_0[i] = ta1_y_zzz_yyy_0[i] * fe_0 - ta1_y_zzz_yyy_1[i] * fe_0 + ta1_y_xzzz_yyy_0[i] * pa_x[i] - ta1_y_xzzz_yyy_1[i] * pc_x[i];

        ta1_y_xxzzz_yyz_0[i] = ta1_y_zzz_yyz_0[i] * fe_0 - ta1_y_zzz_yyz_1[i] * fe_0 + ta1_y_xzzz_yyz_0[i] * pa_x[i] - ta1_y_xzzz_yyz_1[i] * pc_x[i];

        ta1_y_xxzzz_yzz_0[i] = ta1_y_zzz_yzz_0[i] * fe_0 - ta1_y_zzz_yzz_1[i] * fe_0 + ta1_y_xzzz_yzz_0[i] * pa_x[i] - ta1_y_xzzz_yzz_1[i] * pc_x[i];

        ta1_y_xxzzz_zzz_0[i] = ta1_y_zzz_zzz_0[i] * fe_0 - ta1_y_zzz_zzz_1[i] * fe_0 + ta1_y_xzzz_zzz_0[i] * pa_x[i] - ta1_y_xzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : HF

    auto ta1_y_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 310);

    auto ta1_y_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 311);

    auto ta1_y_xyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 312);

    auto ta1_y_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 313);

    auto ta1_y_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 314);

    auto ta1_y_xyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 315);

    auto ta1_y_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 316);

    auto ta1_y_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 317);

    auto ta1_y_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 318);

    auto ta1_y_xyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 319);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyyy_xxx_0, ta1_y_xyyyy_xxy_0, ta1_y_xyyyy_xxz_0, ta1_y_xyyyy_xyy_0, ta1_y_xyyyy_xyz_0, ta1_y_xyyyy_xzz_0, ta1_y_xyyyy_yyy_0, ta1_y_xyyyy_yyz_0, ta1_y_xyyyy_yzz_0, ta1_y_xyyyy_zzz_0, ta1_y_yyyy_xx_0, ta1_y_yyyy_xx_1, ta1_y_yyyy_xxx_0, ta1_y_yyyy_xxx_1, ta1_y_yyyy_xxy_0, ta1_y_yyyy_xxy_1, ta1_y_yyyy_xxz_0, ta1_y_yyyy_xxz_1, ta1_y_yyyy_xy_0, ta1_y_yyyy_xy_1, ta1_y_yyyy_xyy_0, ta1_y_yyyy_xyy_1, ta1_y_yyyy_xyz_0, ta1_y_yyyy_xyz_1, ta1_y_yyyy_xz_0, ta1_y_yyyy_xz_1, ta1_y_yyyy_xzz_0, ta1_y_yyyy_xzz_1, ta1_y_yyyy_yy_0, ta1_y_yyyy_yy_1, ta1_y_yyyy_yyy_0, ta1_y_yyyy_yyy_1, ta1_y_yyyy_yyz_0, ta1_y_yyyy_yyz_1, ta1_y_yyyy_yz_0, ta1_y_yyyy_yz_1, ta1_y_yyyy_yzz_0, ta1_y_yyyy_yzz_1, ta1_y_yyyy_zz_0, ta1_y_yyyy_zz_1, ta1_y_yyyy_zzz_0, ta1_y_yyyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyy_xxx_0[i] = 3.0 * ta1_y_yyyy_xx_0[i] * fe_0 - 3.0 * ta1_y_yyyy_xx_1[i] * fe_0 + ta1_y_yyyy_xxx_0[i] * pa_x[i] - ta1_y_yyyy_xxx_1[i] * pc_x[i];

        ta1_y_xyyyy_xxy_0[i] = 2.0 * ta1_y_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xy_1[i] * fe_0 + ta1_y_yyyy_xxy_0[i] * pa_x[i] - ta1_y_yyyy_xxy_1[i] * pc_x[i];

        ta1_y_xyyyy_xxz_0[i] = 2.0 * ta1_y_yyyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xz_1[i] * fe_0 + ta1_y_yyyy_xxz_0[i] * pa_x[i] - ta1_y_yyyy_xxz_1[i] * pc_x[i];

        ta1_y_xyyyy_xyy_0[i] = ta1_y_yyyy_yy_0[i] * fe_0 - ta1_y_yyyy_yy_1[i] * fe_0 + ta1_y_yyyy_xyy_0[i] * pa_x[i] - ta1_y_yyyy_xyy_1[i] * pc_x[i];

        ta1_y_xyyyy_xyz_0[i] = ta1_y_yyyy_yz_0[i] * fe_0 - ta1_y_yyyy_yz_1[i] * fe_0 + ta1_y_yyyy_xyz_0[i] * pa_x[i] - ta1_y_yyyy_xyz_1[i] * pc_x[i];

        ta1_y_xyyyy_xzz_0[i] = ta1_y_yyyy_zz_0[i] * fe_0 - ta1_y_yyyy_zz_1[i] * fe_0 + ta1_y_yyyy_xzz_0[i] * pa_x[i] - ta1_y_yyyy_xzz_1[i] * pc_x[i];

        ta1_y_xyyyy_yyy_0[i] = ta1_y_yyyy_yyy_0[i] * pa_x[i] - ta1_y_yyyy_yyy_1[i] * pc_x[i];

        ta1_y_xyyyy_yyz_0[i] = ta1_y_yyyy_yyz_0[i] * pa_x[i] - ta1_y_yyyy_yyz_1[i] * pc_x[i];

        ta1_y_xyyyy_yzz_0[i] = ta1_y_yyyy_yzz_0[i] * pa_x[i] - ta1_y_yyyy_yzz_1[i] * pc_x[i];

        ta1_y_xyyyy_zzz_0[i] = ta1_y_yyyy_zzz_0[i] * pa_x[i] - ta1_y_yyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 320-330 components of targeted buffer : HF

    auto ta1_y_xyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 320);

    auto ta1_y_xyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 321);

    auto ta1_y_xyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 322);

    auto ta1_y_xyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 323);

    auto ta1_y_xyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 324);

    auto ta1_y_xyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 325);

    auto ta1_y_xyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 326);

    auto ta1_y_xyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 327);

    auto ta1_y_xyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 328);

    auto ta1_y_xyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 329);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xyyy_xxx_0, ta1_y_xyyy_xxx_1, ta1_y_xyyy_xxy_0, ta1_y_xyyy_xxy_1, ta1_y_xyyy_xyy_0, ta1_y_xyyy_xyy_1, ta1_y_xyyyz_xxx_0, ta1_y_xyyyz_xxy_0, ta1_y_xyyyz_xxz_0, ta1_y_xyyyz_xyy_0, ta1_y_xyyyz_xyz_0, ta1_y_xyyyz_xzz_0, ta1_y_xyyyz_yyy_0, ta1_y_xyyyz_yyz_0, ta1_y_xyyyz_yzz_0, ta1_y_xyyyz_zzz_0, ta1_y_yyyz_xxz_0, ta1_y_yyyz_xxz_1, ta1_y_yyyz_xyz_0, ta1_y_yyyz_xyz_1, ta1_y_yyyz_xz_0, ta1_y_yyyz_xz_1, ta1_y_yyyz_xzz_0, ta1_y_yyyz_xzz_1, ta1_y_yyyz_yyy_0, ta1_y_yyyz_yyy_1, ta1_y_yyyz_yyz_0, ta1_y_yyyz_yyz_1, ta1_y_yyyz_yz_0, ta1_y_yyyz_yz_1, ta1_y_yyyz_yzz_0, ta1_y_yyyz_yzz_1, ta1_y_yyyz_zz_0, ta1_y_yyyz_zz_1, ta1_y_yyyz_zzz_0, ta1_y_yyyz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyz_xxx_0[i] = ta1_y_xyyy_xxx_0[i] * pa_z[i] - ta1_y_xyyy_xxx_1[i] * pc_z[i];

        ta1_y_xyyyz_xxy_0[i] = ta1_y_xyyy_xxy_0[i] * pa_z[i] - ta1_y_xyyy_xxy_1[i] * pc_z[i];

        ta1_y_xyyyz_xxz_0[i] = 2.0 * ta1_y_yyyz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xz_1[i] * fe_0 + ta1_y_yyyz_xxz_0[i] * pa_x[i] - ta1_y_yyyz_xxz_1[i] * pc_x[i];

        ta1_y_xyyyz_xyy_0[i] = ta1_y_xyyy_xyy_0[i] * pa_z[i] - ta1_y_xyyy_xyy_1[i] * pc_z[i];

        ta1_y_xyyyz_xyz_0[i] = ta1_y_yyyz_yz_0[i] * fe_0 - ta1_y_yyyz_yz_1[i] * fe_0 + ta1_y_yyyz_xyz_0[i] * pa_x[i] - ta1_y_yyyz_xyz_1[i] * pc_x[i];

        ta1_y_xyyyz_xzz_0[i] = ta1_y_yyyz_zz_0[i] * fe_0 - ta1_y_yyyz_zz_1[i] * fe_0 + ta1_y_yyyz_xzz_0[i] * pa_x[i] - ta1_y_yyyz_xzz_1[i] * pc_x[i];

        ta1_y_xyyyz_yyy_0[i] = ta1_y_yyyz_yyy_0[i] * pa_x[i] - ta1_y_yyyz_yyy_1[i] * pc_x[i];

        ta1_y_xyyyz_yyz_0[i] = ta1_y_yyyz_yyz_0[i] * pa_x[i] - ta1_y_yyyz_yyz_1[i] * pc_x[i];

        ta1_y_xyyyz_yzz_0[i] = ta1_y_yyyz_yzz_0[i] * pa_x[i] - ta1_y_yyyz_yzz_1[i] * pc_x[i];

        ta1_y_xyyyz_zzz_0[i] = ta1_y_yyyz_zzz_0[i] * pa_x[i] - ta1_y_yyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 330-340 components of targeted buffer : HF

    auto ta1_y_xyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 330);

    auto ta1_y_xyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 331);

    auto ta1_y_xyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 332);

    auto ta1_y_xyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 333);

    auto ta1_y_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 334);

    auto ta1_y_xyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 335);

    auto ta1_y_xyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 336);

    auto ta1_y_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 337);

    auto ta1_y_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 338);

    auto ta1_y_xyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 339);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyzz_xxx_0, ta1_y_xyyzz_xxy_0, ta1_y_xyyzz_xxz_0, ta1_y_xyyzz_xyy_0, ta1_y_xyyzz_xyz_0, ta1_y_xyyzz_xzz_0, ta1_y_xyyzz_yyy_0, ta1_y_xyyzz_yyz_0, ta1_y_xyyzz_yzz_0, ta1_y_xyyzz_zzz_0, ta1_y_yyzz_xx_0, ta1_y_yyzz_xx_1, ta1_y_yyzz_xxx_0, ta1_y_yyzz_xxx_1, ta1_y_yyzz_xxy_0, ta1_y_yyzz_xxy_1, ta1_y_yyzz_xxz_0, ta1_y_yyzz_xxz_1, ta1_y_yyzz_xy_0, ta1_y_yyzz_xy_1, ta1_y_yyzz_xyy_0, ta1_y_yyzz_xyy_1, ta1_y_yyzz_xyz_0, ta1_y_yyzz_xyz_1, ta1_y_yyzz_xz_0, ta1_y_yyzz_xz_1, ta1_y_yyzz_xzz_0, ta1_y_yyzz_xzz_1, ta1_y_yyzz_yy_0, ta1_y_yyzz_yy_1, ta1_y_yyzz_yyy_0, ta1_y_yyzz_yyy_1, ta1_y_yyzz_yyz_0, ta1_y_yyzz_yyz_1, ta1_y_yyzz_yz_0, ta1_y_yyzz_yz_1, ta1_y_yyzz_yzz_0, ta1_y_yyzz_yzz_1, ta1_y_yyzz_zz_0, ta1_y_yyzz_zz_1, ta1_y_yyzz_zzz_0, ta1_y_yyzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzz_xxx_0[i] = 3.0 * ta1_y_yyzz_xx_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xx_1[i] * fe_0 + ta1_y_yyzz_xxx_0[i] * pa_x[i] - ta1_y_yyzz_xxx_1[i] * pc_x[i];

        ta1_y_xyyzz_xxy_0[i] = 2.0 * ta1_y_yyzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yyzz_xy_1[i] * fe_0 + ta1_y_yyzz_xxy_0[i] * pa_x[i] - ta1_y_yyzz_xxy_1[i] * pc_x[i];

        ta1_y_xyyzz_xxz_0[i] = 2.0 * ta1_y_yyzz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyzz_xz_1[i] * fe_0 + ta1_y_yyzz_xxz_0[i] * pa_x[i] - ta1_y_yyzz_xxz_1[i] * pc_x[i];

        ta1_y_xyyzz_xyy_0[i] = ta1_y_yyzz_yy_0[i] * fe_0 - ta1_y_yyzz_yy_1[i] * fe_0 + ta1_y_yyzz_xyy_0[i] * pa_x[i] - ta1_y_yyzz_xyy_1[i] * pc_x[i];

        ta1_y_xyyzz_xyz_0[i] = ta1_y_yyzz_yz_0[i] * fe_0 - ta1_y_yyzz_yz_1[i] * fe_0 + ta1_y_yyzz_xyz_0[i] * pa_x[i] - ta1_y_yyzz_xyz_1[i] * pc_x[i];

        ta1_y_xyyzz_xzz_0[i] = ta1_y_yyzz_zz_0[i] * fe_0 - ta1_y_yyzz_zz_1[i] * fe_0 + ta1_y_yyzz_xzz_0[i] * pa_x[i] - ta1_y_yyzz_xzz_1[i] * pc_x[i];

        ta1_y_xyyzz_yyy_0[i] = ta1_y_yyzz_yyy_0[i] * pa_x[i] - ta1_y_yyzz_yyy_1[i] * pc_x[i];

        ta1_y_xyyzz_yyz_0[i] = ta1_y_yyzz_yyz_0[i] * pa_x[i] - ta1_y_yyzz_yyz_1[i] * pc_x[i];

        ta1_y_xyyzz_yzz_0[i] = ta1_y_yyzz_yzz_0[i] * pa_x[i] - ta1_y_yyzz_yzz_1[i] * pc_x[i];

        ta1_y_xyyzz_zzz_0[i] = ta1_y_yyzz_zzz_0[i] * pa_x[i] - ta1_y_yyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 340-350 components of targeted buffer : HF

    auto ta1_y_xyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 340);

    auto ta1_y_xyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 341);

    auto ta1_y_xyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 342);

    auto ta1_y_xyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 343);

    auto ta1_y_xyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 344);

    auto ta1_y_xyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 345);

    auto ta1_y_xyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 346);

    auto ta1_y_xyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 347);

    auto ta1_y_xyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 348);

    auto ta1_y_xyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 349);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xyzzz_xxx_0, ta1_y_xyzzz_xxy_0, ta1_y_xyzzz_xxz_0, ta1_y_xyzzz_xyy_0, ta1_y_xyzzz_xyz_0, ta1_y_xyzzz_xzz_0, ta1_y_xyzzz_yyy_0, ta1_y_xyzzz_yyz_0, ta1_y_xyzzz_yzz_0, ta1_y_xyzzz_zzz_0, ta1_y_xzzz_xxx_0, ta1_y_xzzz_xxx_1, ta1_y_xzzz_xxz_0, ta1_y_xzzz_xxz_1, ta1_y_xzzz_xzz_0, ta1_y_xzzz_xzz_1, ta1_y_yzzz_xxy_0, ta1_y_yzzz_xxy_1, ta1_y_yzzz_xy_0, ta1_y_yzzz_xy_1, ta1_y_yzzz_xyy_0, ta1_y_yzzz_xyy_1, ta1_y_yzzz_xyz_0, ta1_y_yzzz_xyz_1, ta1_y_yzzz_yy_0, ta1_y_yzzz_yy_1, ta1_y_yzzz_yyy_0, ta1_y_yzzz_yyy_1, ta1_y_yzzz_yyz_0, ta1_y_yzzz_yyz_1, ta1_y_yzzz_yz_0, ta1_y_yzzz_yz_1, ta1_y_yzzz_yzz_0, ta1_y_yzzz_yzz_1, ta1_y_yzzz_zzz_0, ta1_y_yzzz_zzz_1, ta_xzzz_xxx_1, ta_xzzz_xxz_1, ta_xzzz_xzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzz_xxx_0[i] = ta_xzzz_xxx_1[i] + ta1_y_xzzz_xxx_0[i] * pa_y[i] - ta1_y_xzzz_xxx_1[i] * pc_y[i];

        ta1_y_xyzzz_xxy_0[i] = 2.0 * ta1_y_yzzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xy_1[i] * fe_0 + ta1_y_yzzz_xxy_0[i] * pa_x[i] - ta1_y_yzzz_xxy_1[i] * pc_x[i];

        ta1_y_xyzzz_xxz_0[i] = ta_xzzz_xxz_1[i] + ta1_y_xzzz_xxz_0[i] * pa_y[i] - ta1_y_xzzz_xxz_1[i] * pc_y[i];

        ta1_y_xyzzz_xyy_0[i] = ta1_y_yzzz_yy_0[i] * fe_0 - ta1_y_yzzz_yy_1[i] * fe_0 + ta1_y_yzzz_xyy_0[i] * pa_x[i] - ta1_y_yzzz_xyy_1[i] * pc_x[i];

        ta1_y_xyzzz_xyz_0[i] = ta1_y_yzzz_yz_0[i] * fe_0 - ta1_y_yzzz_yz_1[i] * fe_0 + ta1_y_yzzz_xyz_0[i] * pa_x[i] - ta1_y_yzzz_xyz_1[i] * pc_x[i];

        ta1_y_xyzzz_xzz_0[i] = ta_xzzz_xzz_1[i] + ta1_y_xzzz_xzz_0[i] * pa_y[i] - ta1_y_xzzz_xzz_1[i] * pc_y[i];

        ta1_y_xyzzz_yyy_0[i] = ta1_y_yzzz_yyy_0[i] * pa_x[i] - ta1_y_yzzz_yyy_1[i] * pc_x[i];

        ta1_y_xyzzz_yyz_0[i] = ta1_y_yzzz_yyz_0[i] * pa_x[i] - ta1_y_yzzz_yyz_1[i] * pc_x[i];

        ta1_y_xyzzz_yzz_0[i] = ta1_y_yzzz_yzz_0[i] * pa_x[i] - ta1_y_yzzz_yzz_1[i] * pc_x[i];

        ta1_y_xyzzz_zzz_0[i] = ta1_y_yzzz_zzz_0[i] * pa_x[i] - ta1_y_yzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 350-360 components of targeted buffer : HF

    auto ta1_y_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 350);

    auto ta1_y_xzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 351);

    auto ta1_y_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 352);

    auto ta1_y_xzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 353);

    auto ta1_y_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 354);

    auto ta1_y_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 355);

    auto ta1_y_xzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 356);

    auto ta1_y_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 357);

    auto ta1_y_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 358);

    auto ta1_y_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 359);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xzzzz_xxx_0, ta1_y_xzzzz_xxy_0, ta1_y_xzzzz_xxz_0, ta1_y_xzzzz_xyy_0, ta1_y_xzzzz_xyz_0, ta1_y_xzzzz_xzz_0, ta1_y_xzzzz_yyy_0, ta1_y_xzzzz_yyz_0, ta1_y_xzzzz_yzz_0, ta1_y_xzzzz_zzz_0, ta1_y_zzzz_xx_0, ta1_y_zzzz_xx_1, ta1_y_zzzz_xxx_0, ta1_y_zzzz_xxx_1, ta1_y_zzzz_xxy_0, ta1_y_zzzz_xxy_1, ta1_y_zzzz_xxz_0, ta1_y_zzzz_xxz_1, ta1_y_zzzz_xy_0, ta1_y_zzzz_xy_1, ta1_y_zzzz_xyy_0, ta1_y_zzzz_xyy_1, ta1_y_zzzz_xyz_0, ta1_y_zzzz_xyz_1, ta1_y_zzzz_xz_0, ta1_y_zzzz_xz_1, ta1_y_zzzz_xzz_0, ta1_y_zzzz_xzz_1, ta1_y_zzzz_yy_0, ta1_y_zzzz_yy_1, ta1_y_zzzz_yyy_0, ta1_y_zzzz_yyy_1, ta1_y_zzzz_yyz_0, ta1_y_zzzz_yyz_1, ta1_y_zzzz_yz_0, ta1_y_zzzz_yz_1, ta1_y_zzzz_yzz_0, ta1_y_zzzz_yzz_1, ta1_y_zzzz_zz_0, ta1_y_zzzz_zz_1, ta1_y_zzzz_zzz_0, ta1_y_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzz_xxx_0[i] = 3.0 * ta1_y_zzzz_xx_0[i] * fe_0 - 3.0 * ta1_y_zzzz_xx_1[i] * fe_0 + ta1_y_zzzz_xxx_0[i] * pa_x[i] - ta1_y_zzzz_xxx_1[i] * pc_x[i];

        ta1_y_xzzzz_xxy_0[i] = 2.0 * ta1_y_zzzz_xy_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xy_1[i] * fe_0 + ta1_y_zzzz_xxy_0[i] * pa_x[i] - ta1_y_zzzz_xxy_1[i] * pc_x[i];

        ta1_y_xzzzz_xxz_0[i] = 2.0 * ta1_y_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xz_1[i] * fe_0 + ta1_y_zzzz_xxz_0[i] * pa_x[i] - ta1_y_zzzz_xxz_1[i] * pc_x[i];

        ta1_y_xzzzz_xyy_0[i] = ta1_y_zzzz_yy_0[i] * fe_0 - ta1_y_zzzz_yy_1[i] * fe_0 + ta1_y_zzzz_xyy_0[i] * pa_x[i] - ta1_y_zzzz_xyy_1[i] * pc_x[i];

        ta1_y_xzzzz_xyz_0[i] = ta1_y_zzzz_yz_0[i] * fe_0 - ta1_y_zzzz_yz_1[i] * fe_0 + ta1_y_zzzz_xyz_0[i] * pa_x[i] - ta1_y_zzzz_xyz_1[i] * pc_x[i];

        ta1_y_xzzzz_xzz_0[i] = ta1_y_zzzz_zz_0[i] * fe_0 - ta1_y_zzzz_zz_1[i] * fe_0 + ta1_y_zzzz_xzz_0[i] * pa_x[i] - ta1_y_zzzz_xzz_1[i] * pc_x[i];

        ta1_y_xzzzz_yyy_0[i] = ta1_y_zzzz_yyy_0[i] * pa_x[i] - ta1_y_zzzz_yyy_1[i] * pc_x[i];

        ta1_y_xzzzz_yyz_0[i] = ta1_y_zzzz_yyz_0[i] * pa_x[i] - ta1_y_zzzz_yyz_1[i] * pc_x[i];

        ta1_y_xzzzz_yzz_0[i] = ta1_y_zzzz_yzz_0[i] * pa_x[i] - ta1_y_zzzz_yzz_1[i] * pc_x[i];

        ta1_y_xzzzz_zzz_0[i] = ta1_y_zzzz_zzz_0[i] * pa_x[i] - ta1_y_zzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 360-370 components of targeted buffer : HF

    auto ta1_y_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 360);

    auto ta1_y_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 361);

    auto ta1_y_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 362);

    auto ta1_y_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 363);

    auto ta1_y_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 364);

    auto ta1_y_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 365);

    auto ta1_y_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 366);

    auto ta1_y_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 367);

    auto ta1_y_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 368);

    auto ta1_y_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 369);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yyy_xxx_0, ta1_y_yyy_xxx_1, ta1_y_yyy_xxy_0, ta1_y_yyy_xxy_1, ta1_y_yyy_xxz_0, ta1_y_yyy_xxz_1, ta1_y_yyy_xyy_0, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_0, ta1_y_yyy_xyz_1, ta1_y_yyy_xzz_0, ta1_y_yyy_xzz_1, ta1_y_yyy_yyy_0, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_0, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_0, ta1_y_yyy_yzz_1, ta1_y_yyy_zzz_0, ta1_y_yyy_zzz_1, ta1_y_yyyy_xx_0, ta1_y_yyyy_xx_1, ta1_y_yyyy_xxx_0, ta1_y_yyyy_xxx_1, ta1_y_yyyy_xxy_0, ta1_y_yyyy_xxy_1, ta1_y_yyyy_xxz_0, ta1_y_yyyy_xxz_1, ta1_y_yyyy_xy_0, ta1_y_yyyy_xy_1, ta1_y_yyyy_xyy_0, ta1_y_yyyy_xyy_1, ta1_y_yyyy_xyz_0, ta1_y_yyyy_xyz_1, ta1_y_yyyy_xz_0, ta1_y_yyyy_xz_1, ta1_y_yyyy_xzz_0, ta1_y_yyyy_xzz_1, ta1_y_yyyy_yy_0, ta1_y_yyyy_yy_1, ta1_y_yyyy_yyy_0, ta1_y_yyyy_yyy_1, ta1_y_yyyy_yyz_0, ta1_y_yyyy_yyz_1, ta1_y_yyyy_yz_0, ta1_y_yyyy_yz_1, ta1_y_yyyy_yzz_0, ta1_y_yyyy_yzz_1, ta1_y_yyyy_zz_0, ta1_y_yyyy_zz_1, ta1_y_yyyy_zzz_0, ta1_y_yyyy_zzz_1, ta1_y_yyyyy_xxx_0, ta1_y_yyyyy_xxy_0, ta1_y_yyyyy_xxz_0, ta1_y_yyyyy_xyy_0, ta1_y_yyyyy_xyz_0, ta1_y_yyyyy_xzz_0, ta1_y_yyyyy_yyy_0, ta1_y_yyyyy_yyz_0, ta1_y_yyyyy_yzz_0, ta1_y_yyyyy_zzz_0, ta_yyyy_xxx_1, ta_yyyy_xxy_1, ta_yyyy_xxz_1, ta_yyyy_xyy_1, ta_yyyy_xyz_1, ta_yyyy_xzz_1, ta_yyyy_yyy_1, ta_yyyy_yyz_1, ta_yyyy_yzz_1, ta_yyyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyy_xxx_0[i] = 4.0 * ta1_y_yyy_xxx_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxx_1[i] * fe_0 + ta_yyyy_xxx_1[i] + ta1_y_yyyy_xxx_0[i] * pa_y[i] - ta1_y_yyyy_xxx_1[i] * pc_y[i];

        ta1_y_yyyyy_xxy_0[i] = 4.0 * ta1_y_yyy_xxy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxy_1[i] * fe_0 + ta1_y_yyyy_xx_0[i] * fe_0 - ta1_y_yyyy_xx_1[i] * fe_0 + ta_yyyy_xxy_1[i] + ta1_y_yyyy_xxy_0[i] * pa_y[i] - ta1_y_yyyy_xxy_1[i] * pc_y[i];

        ta1_y_yyyyy_xxz_0[i] = 4.0 * ta1_y_yyy_xxz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxz_1[i] * fe_0 + ta_yyyy_xxz_1[i] + ta1_y_yyyy_xxz_0[i] * pa_y[i] - ta1_y_yyyy_xxz_1[i] * pc_y[i];

        ta1_y_yyyyy_xyy_0[i] = 4.0 * ta1_y_yyy_xyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_y_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xy_1[i] * fe_0 + ta_yyyy_xyy_1[i] + ta1_y_yyyy_xyy_0[i] * pa_y[i] - ta1_y_yyyy_xyy_1[i] * pc_y[i];

        ta1_y_yyyyy_xyz_0[i] = 4.0 * ta1_y_yyy_xyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyz_1[i] * fe_0 + ta1_y_yyyy_xz_0[i] * fe_0 - ta1_y_yyyy_xz_1[i] * fe_0 + ta_yyyy_xyz_1[i] + ta1_y_yyyy_xyz_0[i] * pa_y[i] - ta1_y_yyyy_xyz_1[i] * pc_y[i];

        ta1_y_yyyyy_xzz_0[i] = 4.0 * ta1_y_yyy_xzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xzz_1[i] * fe_0 + ta_yyyy_xzz_1[i] + ta1_y_yyyy_xzz_0[i] * pa_y[i] - ta1_y_yyyy_xzz_1[i] * pc_y[i];

        ta1_y_yyyyy_yyy_0[i] = 4.0 * ta1_y_yyy_yyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyy_1[i] * fe_0 + 3.0 * ta1_y_yyyy_yy_0[i] * fe_0 - 3.0 * ta1_y_yyyy_yy_1[i] * fe_0 + ta_yyyy_yyy_1[i] + ta1_y_yyyy_yyy_0[i] * pa_y[i] - ta1_y_yyyy_yyy_1[i] * pc_y[i];

        ta1_y_yyyyy_yyz_0[i] = 4.0 * ta1_y_yyy_yyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_y_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_yz_1[i] * fe_0 + ta_yyyy_yyz_1[i] + ta1_y_yyyy_yyz_0[i] * pa_y[i] - ta1_y_yyyy_yyz_1[i] * pc_y[i];

        ta1_y_yyyyy_yzz_0[i] = 4.0 * ta1_y_yyy_yzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yzz_1[i] * fe_0 + ta1_y_yyyy_zz_0[i] * fe_0 - ta1_y_yyyy_zz_1[i] * fe_0 + ta_yyyy_yzz_1[i] + ta1_y_yyyy_yzz_0[i] * pa_y[i] - ta1_y_yyyy_yzz_1[i] * pc_y[i];

        ta1_y_yyyyy_zzz_0[i] = 4.0 * ta1_y_yyy_zzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_zzz_1[i] * fe_0 + ta_yyyy_zzz_1[i] + ta1_y_yyyy_zzz_0[i] * pa_y[i] - ta1_y_yyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 370-380 components of targeted buffer : HF

    auto ta1_y_yyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 370);

    auto ta1_y_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 371);

    auto ta1_y_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 372);

    auto ta1_y_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 373);

    auto ta1_y_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 374);

    auto ta1_y_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 375);

    auto ta1_y_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 376);

    auto ta1_y_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 377);

    auto ta1_y_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 378);

    auto ta1_y_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 379);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_yyyy_xx_0, ta1_y_yyyy_xx_1, ta1_y_yyyy_xxx_0, ta1_y_yyyy_xxx_1, ta1_y_yyyy_xxy_0, ta1_y_yyyy_xxy_1, ta1_y_yyyy_xxz_0, ta1_y_yyyy_xxz_1, ta1_y_yyyy_xy_0, ta1_y_yyyy_xy_1, ta1_y_yyyy_xyy_0, ta1_y_yyyy_xyy_1, ta1_y_yyyy_xyz_0, ta1_y_yyyy_xyz_1, ta1_y_yyyy_xz_0, ta1_y_yyyy_xz_1, ta1_y_yyyy_xzz_0, ta1_y_yyyy_xzz_1, ta1_y_yyyy_yy_0, ta1_y_yyyy_yy_1, ta1_y_yyyy_yyy_0, ta1_y_yyyy_yyy_1, ta1_y_yyyy_yyz_0, ta1_y_yyyy_yyz_1, ta1_y_yyyy_yz_0, ta1_y_yyyy_yz_1, ta1_y_yyyy_yzz_0, ta1_y_yyyy_yzz_1, ta1_y_yyyy_zz_0, ta1_y_yyyy_zz_1, ta1_y_yyyy_zzz_0, ta1_y_yyyy_zzz_1, ta1_y_yyyyz_xxx_0, ta1_y_yyyyz_xxy_0, ta1_y_yyyyz_xxz_0, ta1_y_yyyyz_xyy_0, ta1_y_yyyyz_xyz_0, ta1_y_yyyyz_xzz_0, ta1_y_yyyyz_yyy_0, ta1_y_yyyyz_yyz_0, ta1_y_yyyyz_yzz_0, ta1_y_yyyyz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyz_xxx_0[i] = ta1_y_yyyy_xxx_0[i] * pa_z[i] - ta1_y_yyyy_xxx_1[i] * pc_z[i];

        ta1_y_yyyyz_xxy_0[i] = ta1_y_yyyy_xxy_0[i] * pa_z[i] - ta1_y_yyyy_xxy_1[i] * pc_z[i];

        ta1_y_yyyyz_xxz_0[i] = ta1_y_yyyy_xx_0[i] * fe_0 - ta1_y_yyyy_xx_1[i] * fe_0 + ta1_y_yyyy_xxz_0[i] * pa_z[i] - ta1_y_yyyy_xxz_1[i] * pc_z[i];

        ta1_y_yyyyz_xyy_0[i] = ta1_y_yyyy_xyy_0[i] * pa_z[i] - ta1_y_yyyy_xyy_1[i] * pc_z[i];

        ta1_y_yyyyz_xyz_0[i] = ta1_y_yyyy_xy_0[i] * fe_0 - ta1_y_yyyy_xy_1[i] * fe_0 + ta1_y_yyyy_xyz_0[i] * pa_z[i] - ta1_y_yyyy_xyz_1[i] * pc_z[i];

        ta1_y_yyyyz_xzz_0[i] = 2.0 * ta1_y_yyyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_xz_1[i] * fe_0 + ta1_y_yyyy_xzz_0[i] * pa_z[i] - ta1_y_yyyy_xzz_1[i] * pc_z[i];

        ta1_y_yyyyz_yyy_0[i] = ta1_y_yyyy_yyy_0[i] * pa_z[i] - ta1_y_yyyy_yyy_1[i] * pc_z[i];

        ta1_y_yyyyz_yyz_0[i] = ta1_y_yyyy_yy_0[i] * fe_0 - ta1_y_yyyy_yy_1[i] * fe_0 + ta1_y_yyyy_yyz_0[i] * pa_z[i] - ta1_y_yyyy_yyz_1[i] * pc_z[i];

        ta1_y_yyyyz_yzz_0[i] = 2.0 * ta1_y_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_y_yyyy_yz_1[i] * fe_0 + ta1_y_yyyy_yzz_0[i] * pa_z[i] - ta1_y_yyyy_yzz_1[i] * pc_z[i];

        ta1_y_yyyyz_zzz_0[i] = 3.0 * ta1_y_yyyy_zz_0[i] * fe_0 - 3.0 * ta1_y_yyyy_zz_1[i] * fe_0 + ta1_y_yyyy_zzz_0[i] * pa_z[i] - ta1_y_yyyy_zzz_1[i] * pc_z[i];
    }

    // Set up 380-390 components of targeted buffer : HF

    auto ta1_y_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 380);

    auto ta1_y_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 381);

    auto ta1_y_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 382);

    auto ta1_y_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 383);

    auto ta1_y_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 384);

    auto ta1_y_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 385);

    auto ta1_y_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 386);

    auto ta1_y_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 387);

    auto ta1_y_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 388);

    auto ta1_y_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 389);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyy_xxx_0, ta1_y_yyy_xxx_1, ta1_y_yyy_xxy_0, ta1_y_yyy_xxy_1, ta1_y_yyy_xyy_0, ta1_y_yyy_xyy_1, ta1_y_yyy_xyz_0, ta1_y_yyy_xyz_1, ta1_y_yyy_yyy_0, ta1_y_yyy_yyy_1, ta1_y_yyy_yyz_0, ta1_y_yyy_yyz_1, ta1_y_yyy_yzz_0, ta1_y_yyy_yzz_1, ta1_y_yyyz_xxx_0, ta1_y_yyyz_xxx_1, ta1_y_yyyz_xxy_0, ta1_y_yyyz_xxy_1, ta1_y_yyyz_xy_0, ta1_y_yyyz_xy_1, ta1_y_yyyz_xyy_0, ta1_y_yyyz_xyy_1, ta1_y_yyyz_xyz_0, ta1_y_yyyz_xyz_1, ta1_y_yyyz_yy_0, ta1_y_yyyz_yy_1, ta1_y_yyyz_yyy_0, ta1_y_yyyz_yyy_1, ta1_y_yyyz_yyz_0, ta1_y_yyyz_yyz_1, ta1_y_yyyz_yz_0, ta1_y_yyyz_yz_1, ta1_y_yyyz_yzz_0, ta1_y_yyyz_yzz_1, ta1_y_yyyzz_xxx_0, ta1_y_yyyzz_xxy_0, ta1_y_yyyzz_xxz_0, ta1_y_yyyzz_xyy_0, ta1_y_yyyzz_xyz_0, ta1_y_yyyzz_xzz_0, ta1_y_yyyzz_yyy_0, ta1_y_yyyzz_yyz_0, ta1_y_yyyzz_yzz_0, ta1_y_yyyzz_zzz_0, ta1_y_yyzz_xxz_0, ta1_y_yyzz_xxz_1, ta1_y_yyzz_xzz_0, ta1_y_yyzz_xzz_1, ta1_y_yyzz_zzz_0, ta1_y_yyzz_zzz_1, ta1_y_yzz_xxz_0, ta1_y_yzz_xxz_1, ta1_y_yzz_xzz_0, ta1_y_yzz_xzz_1, ta1_y_yzz_zzz_0, ta1_y_yzz_zzz_1, ta_yyzz_xxz_1, ta_yyzz_xzz_1, ta_yyzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzz_xxx_0[i] = ta1_y_yyy_xxx_0[i] * fe_0 - ta1_y_yyy_xxx_1[i] * fe_0 + ta1_y_yyyz_xxx_0[i] * pa_z[i] - ta1_y_yyyz_xxx_1[i] * pc_z[i];

        ta1_y_yyyzz_xxy_0[i] = ta1_y_yyy_xxy_0[i] * fe_0 - ta1_y_yyy_xxy_1[i] * fe_0 + ta1_y_yyyz_xxy_0[i] * pa_z[i] - ta1_y_yyyz_xxy_1[i] * pc_z[i];

        ta1_y_yyyzz_xxz_0[i] = 2.0 * ta1_y_yzz_xxz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xxz_1[i] * fe_0 + ta_yyzz_xxz_1[i] + ta1_y_yyzz_xxz_0[i] * pa_y[i] - ta1_y_yyzz_xxz_1[i] * pc_y[i];

        ta1_y_yyyzz_xyy_0[i] = ta1_y_yyy_xyy_0[i] * fe_0 - ta1_y_yyy_xyy_1[i] * fe_0 + ta1_y_yyyz_xyy_0[i] * pa_z[i] - ta1_y_yyyz_xyy_1[i] * pc_z[i];

        ta1_y_yyyzz_xyz_0[i] = ta1_y_yyy_xyz_0[i] * fe_0 - ta1_y_yyy_xyz_1[i] * fe_0 + ta1_y_yyyz_xy_0[i] * fe_0 - ta1_y_yyyz_xy_1[i] * fe_0 + ta1_y_yyyz_xyz_0[i] * pa_z[i] - ta1_y_yyyz_xyz_1[i] * pc_z[i];

        ta1_y_yyyzz_xzz_0[i] = 2.0 * ta1_y_yzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xzz_1[i] * fe_0 + ta_yyzz_xzz_1[i] + ta1_y_yyzz_xzz_0[i] * pa_y[i] - ta1_y_yyzz_xzz_1[i] * pc_y[i];

        ta1_y_yyyzz_yyy_0[i] = ta1_y_yyy_yyy_0[i] * fe_0 - ta1_y_yyy_yyy_1[i] * fe_0 + ta1_y_yyyz_yyy_0[i] * pa_z[i] - ta1_y_yyyz_yyy_1[i] * pc_z[i];

        ta1_y_yyyzz_yyz_0[i] = ta1_y_yyy_yyz_0[i] * fe_0 - ta1_y_yyy_yyz_1[i] * fe_0 + ta1_y_yyyz_yy_0[i] * fe_0 - ta1_y_yyyz_yy_1[i] * fe_0 + ta1_y_yyyz_yyz_0[i] * pa_z[i] - ta1_y_yyyz_yyz_1[i] * pc_z[i];

        ta1_y_yyyzz_yzz_0[i] = ta1_y_yyy_yzz_0[i] * fe_0 - ta1_y_yyy_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyyz_yz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yz_1[i] * fe_0 + ta1_y_yyyz_yzz_0[i] * pa_z[i] - ta1_y_yyyz_yzz_1[i] * pc_z[i];

        ta1_y_yyyzz_zzz_0[i] = 2.0 * ta1_y_yzz_zzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_zzz_1[i] * fe_0 + ta_yyzz_zzz_1[i] + ta1_y_yyzz_zzz_0[i] * pa_y[i] - ta1_y_yyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 390-400 components of targeted buffer : HF

    auto ta1_y_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 390);

    auto ta1_y_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 391);

    auto ta1_y_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 392);

    auto ta1_y_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 393);

    auto ta1_y_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 394);

    auto ta1_y_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 395);

    auto ta1_y_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 396);

    auto ta1_y_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 397);

    auto ta1_y_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 398);

    auto ta1_y_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 399);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyz_xxx_0, ta1_y_yyz_xxx_1, ta1_y_yyz_xxy_0, ta1_y_yyz_xxy_1, ta1_y_yyz_xyy_0, ta1_y_yyz_xyy_1, ta1_y_yyz_xyz_0, ta1_y_yyz_xyz_1, ta1_y_yyz_yyy_0, ta1_y_yyz_yyy_1, ta1_y_yyz_yyz_0, ta1_y_yyz_yyz_1, ta1_y_yyz_yzz_0, ta1_y_yyz_yzz_1, ta1_y_yyzz_xxx_0, ta1_y_yyzz_xxx_1, ta1_y_yyzz_xxy_0, ta1_y_yyzz_xxy_1, ta1_y_yyzz_xy_0, ta1_y_yyzz_xy_1, ta1_y_yyzz_xyy_0, ta1_y_yyzz_xyy_1, ta1_y_yyzz_xyz_0, ta1_y_yyzz_xyz_1, ta1_y_yyzz_yy_0, ta1_y_yyzz_yy_1, ta1_y_yyzz_yyy_0, ta1_y_yyzz_yyy_1, ta1_y_yyzz_yyz_0, ta1_y_yyzz_yyz_1, ta1_y_yyzz_yz_0, ta1_y_yyzz_yz_1, ta1_y_yyzz_yzz_0, ta1_y_yyzz_yzz_1, ta1_y_yyzzz_xxx_0, ta1_y_yyzzz_xxy_0, ta1_y_yyzzz_xxz_0, ta1_y_yyzzz_xyy_0, ta1_y_yyzzz_xyz_0, ta1_y_yyzzz_xzz_0, ta1_y_yyzzz_yyy_0, ta1_y_yyzzz_yyz_0, ta1_y_yyzzz_yzz_0, ta1_y_yyzzz_zzz_0, ta1_y_yzzz_xxz_0, ta1_y_yzzz_xxz_1, ta1_y_yzzz_xzz_0, ta1_y_yzzz_xzz_1, ta1_y_yzzz_zzz_0, ta1_y_yzzz_zzz_1, ta1_y_zzz_xxz_0, ta1_y_zzz_xxz_1, ta1_y_zzz_xzz_0, ta1_y_zzz_xzz_1, ta1_y_zzz_zzz_0, ta1_y_zzz_zzz_1, ta_yzzz_xxz_1, ta_yzzz_xzz_1, ta_yzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzz_xxx_0[i] = 2.0 * ta1_y_yyz_xxx_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxx_1[i] * fe_0 + ta1_y_yyzz_xxx_0[i] * pa_z[i] - ta1_y_yyzz_xxx_1[i] * pc_z[i];

        ta1_y_yyzzz_xxy_0[i] = 2.0 * ta1_y_yyz_xxy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xxy_1[i] * fe_0 + ta1_y_yyzz_xxy_0[i] * pa_z[i] - ta1_y_yyzz_xxy_1[i] * pc_z[i];

        ta1_y_yyzzz_xxz_0[i] = ta1_y_zzz_xxz_0[i] * fe_0 - ta1_y_zzz_xxz_1[i] * fe_0 + ta_yzzz_xxz_1[i] + ta1_y_yzzz_xxz_0[i] * pa_y[i] - ta1_y_yzzz_xxz_1[i] * pc_y[i];

        ta1_y_yyzzz_xyy_0[i] = 2.0 * ta1_y_yyz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyy_1[i] * fe_0 + ta1_y_yyzz_xyy_0[i] * pa_z[i] - ta1_y_yyzz_xyy_1[i] * pc_z[i];

        ta1_y_yyzzz_xyz_0[i] = 2.0 * ta1_y_yyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyz_1[i] * fe_0 + ta1_y_yyzz_xy_0[i] * fe_0 - ta1_y_yyzz_xy_1[i] * fe_0 + ta1_y_yyzz_xyz_0[i] * pa_z[i] - ta1_y_yyzz_xyz_1[i] * pc_z[i];

        ta1_y_yyzzz_xzz_0[i] = ta1_y_zzz_xzz_0[i] * fe_0 - ta1_y_zzz_xzz_1[i] * fe_0 + ta_yzzz_xzz_1[i] + ta1_y_yzzz_xzz_0[i] * pa_y[i] - ta1_y_yzzz_xzz_1[i] * pc_y[i];

        ta1_y_yyzzz_yyy_0[i] = 2.0 * ta1_y_yyz_yyy_0[i] * fe_0 - 2.0 * ta1_y_yyz_yyy_1[i] * fe_0 + ta1_y_yyzz_yyy_0[i] * pa_z[i] - ta1_y_yyzz_yyy_1[i] * pc_z[i];

        ta1_y_yyzzz_yyz_0[i] = 2.0 * ta1_y_yyz_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yyz_1[i] * fe_0 + ta1_y_yyzz_yy_0[i] * fe_0 - ta1_y_yyzz_yy_1[i] * fe_0 + ta1_y_yyzz_yyz_0[i] * pa_z[i] - ta1_y_yyzz_yyz_1[i] * pc_z[i];

        ta1_y_yyzzz_yzz_0[i] = 2.0 * ta1_y_yyz_yzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyzz_yz_0[i] * fe_0 - 2.0 * ta1_y_yyzz_yz_1[i] * fe_0 + ta1_y_yyzz_yzz_0[i] * pa_z[i] - ta1_y_yyzz_yzz_1[i] * pc_z[i];

        ta1_y_yyzzz_zzz_0[i] = ta1_y_zzz_zzz_0[i] * fe_0 - ta1_y_zzz_zzz_1[i] * fe_0 + ta_yzzz_zzz_1[i] + ta1_y_yzzz_zzz_0[i] * pa_y[i] - ta1_y_yzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 400-410 components of targeted buffer : HF

    auto ta1_y_yzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 400);

    auto ta1_y_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 401);

    auto ta1_y_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 402);

    auto ta1_y_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 403);

    auto ta1_y_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 404);

    auto ta1_y_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 405);

    auto ta1_y_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 406);

    auto ta1_y_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 407);

    auto ta1_y_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 408);

    auto ta1_y_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 409);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yzz_xxy_0, ta1_y_yzz_xxy_1, ta1_y_yzz_xyy_0, ta1_y_yzz_xyy_1, ta1_y_yzz_yyy_0, ta1_y_yzz_yyy_1, ta1_y_yzzz_xxy_0, ta1_y_yzzz_xxy_1, ta1_y_yzzz_xyy_0, ta1_y_yzzz_xyy_1, ta1_y_yzzz_yyy_0, ta1_y_yzzz_yyy_1, ta1_y_yzzzz_xxx_0, ta1_y_yzzzz_xxy_0, ta1_y_yzzzz_xxz_0, ta1_y_yzzzz_xyy_0, ta1_y_yzzzz_xyz_0, ta1_y_yzzzz_xzz_0, ta1_y_yzzzz_yyy_0, ta1_y_yzzzz_yyz_0, ta1_y_yzzzz_yzz_0, ta1_y_yzzzz_zzz_0, ta1_y_zzzz_xxx_0, ta1_y_zzzz_xxx_1, ta1_y_zzzz_xxz_0, ta1_y_zzzz_xxz_1, ta1_y_zzzz_xyz_0, ta1_y_zzzz_xyz_1, ta1_y_zzzz_xz_0, ta1_y_zzzz_xz_1, ta1_y_zzzz_xzz_0, ta1_y_zzzz_xzz_1, ta1_y_zzzz_yyz_0, ta1_y_zzzz_yyz_1, ta1_y_zzzz_yz_0, ta1_y_zzzz_yz_1, ta1_y_zzzz_yzz_0, ta1_y_zzzz_yzz_1, ta1_y_zzzz_zz_0, ta1_y_zzzz_zz_1, ta1_y_zzzz_zzz_0, ta1_y_zzzz_zzz_1, ta_zzzz_xxx_1, ta_zzzz_xxz_1, ta_zzzz_xyz_1, ta_zzzz_xzz_1, ta_zzzz_yyz_1, ta_zzzz_yzz_1, ta_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzz_xxx_0[i] = ta_zzzz_xxx_1[i] + ta1_y_zzzz_xxx_0[i] * pa_y[i] - ta1_y_zzzz_xxx_1[i] * pc_y[i];

        ta1_y_yzzzz_xxy_0[i] = 3.0 * ta1_y_yzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxy_1[i] * fe_0 + ta1_y_yzzz_xxy_0[i] * pa_z[i] - ta1_y_yzzz_xxy_1[i] * pc_z[i];

        ta1_y_yzzzz_xxz_0[i] = ta_zzzz_xxz_1[i] + ta1_y_zzzz_xxz_0[i] * pa_y[i] - ta1_y_zzzz_xxz_1[i] * pc_y[i];

        ta1_y_yzzzz_xyy_0[i] = 3.0 * ta1_y_yzz_xyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xyy_1[i] * fe_0 + ta1_y_yzzz_xyy_0[i] * pa_z[i] - ta1_y_yzzz_xyy_1[i] * pc_z[i];

        ta1_y_yzzzz_xyz_0[i] = ta1_y_zzzz_xz_0[i] * fe_0 - ta1_y_zzzz_xz_1[i] * fe_0 + ta_zzzz_xyz_1[i] + ta1_y_zzzz_xyz_0[i] * pa_y[i] - ta1_y_zzzz_xyz_1[i] * pc_y[i];

        ta1_y_yzzzz_xzz_0[i] = ta_zzzz_xzz_1[i] + ta1_y_zzzz_xzz_0[i] * pa_y[i] - ta1_y_zzzz_xzz_1[i] * pc_y[i];

        ta1_y_yzzzz_yyy_0[i] = 3.0 * ta1_y_yzz_yyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_yyy_1[i] * fe_0 + ta1_y_yzzz_yyy_0[i] * pa_z[i] - ta1_y_yzzz_yyy_1[i] * pc_z[i];

        ta1_y_yzzzz_yyz_0[i] = 2.0 * ta1_y_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_yz_1[i] * fe_0 + ta_zzzz_yyz_1[i] + ta1_y_zzzz_yyz_0[i] * pa_y[i] - ta1_y_zzzz_yyz_1[i] * pc_y[i];

        ta1_y_yzzzz_yzz_0[i] = ta1_y_zzzz_zz_0[i] * fe_0 - ta1_y_zzzz_zz_1[i] * fe_0 + ta_zzzz_yzz_1[i] + ta1_y_zzzz_yzz_0[i] * pa_y[i] - ta1_y_zzzz_yzz_1[i] * pc_y[i];

        ta1_y_yzzzz_zzz_0[i] = ta_zzzz_zzz_1[i] + ta1_y_zzzz_zzz_0[i] * pa_y[i] - ta1_y_zzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 410-420 components of targeted buffer : HF

    auto ta1_y_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 410);

    auto ta1_y_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 411);

    auto ta1_y_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 412);

    auto ta1_y_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 413);

    auto ta1_y_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 414);

    auto ta1_y_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 415);

    auto ta1_y_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 416);

    auto ta1_y_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 417);

    auto ta1_y_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 418);

    auto ta1_y_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 419);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zzz_xxx_0, ta1_y_zzz_xxx_1, ta1_y_zzz_xxy_0, ta1_y_zzz_xxy_1, ta1_y_zzz_xxz_0, ta1_y_zzz_xxz_1, ta1_y_zzz_xyy_0, ta1_y_zzz_xyy_1, ta1_y_zzz_xyz_0, ta1_y_zzz_xyz_1, ta1_y_zzz_xzz_0, ta1_y_zzz_xzz_1, ta1_y_zzz_yyy_0, ta1_y_zzz_yyy_1, ta1_y_zzz_yyz_0, ta1_y_zzz_yyz_1, ta1_y_zzz_yzz_0, ta1_y_zzz_yzz_1, ta1_y_zzz_zzz_0, ta1_y_zzz_zzz_1, ta1_y_zzzz_xx_0, ta1_y_zzzz_xx_1, ta1_y_zzzz_xxx_0, ta1_y_zzzz_xxx_1, ta1_y_zzzz_xxy_0, ta1_y_zzzz_xxy_1, ta1_y_zzzz_xxz_0, ta1_y_zzzz_xxz_1, ta1_y_zzzz_xy_0, ta1_y_zzzz_xy_1, ta1_y_zzzz_xyy_0, ta1_y_zzzz_xyy_1, ta1_y_zzzz_xyz_0, ta1_y_zzzz_xyz_1, ta1_y_zzzz_xz_0, ta1_y_zzzz_xz_1, ta1_y_zzzz_xzz_0, ta1_y_zzzz_xzz_1, ta1_y_zzzz_yy_0, ta1_y_zzzz_yy_1, ta1_y_zzzz_yyy_0, ta1_y_zzzz_yyy_1, ta1_y_zzzz_yyz_0, ta1_y_zzzz_yyz_1, ta1_y_zzzz_yz_0, ta1_y_zzzz_yz_1, ta1_y_zzzz_yzz_0, ta1_y_zzzz_yzz_1, ta1_y_zzzz_zz_0, ta1_y_zzzz_zz_1, ta1_y_zzzz_zzz_0, ta1_y_zzzz_zzz_1, ta1_y_zzzzz_xxx_0, ta1_y_zzzzz_xxy_0, ta1_y_zzzzz_xxz_0, ta1_y_zzzzz_xyy_0, ta1_y_zzzzz_xyz_0, ta1_y_zzzzz_xzz_0, ta1_y_zzzzz_yyy_0, ta1_y_zzzzz_yyz_0, ta1_y_zzzzz_yzz_0, ta1_y_zzzzz_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzz_xxx_0[i] = 4.0 * ta1_y_zzz_xxx_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxx_1[i] * fe_0 + ta1_y_zzzz_xxx_0[i] * pa_z[i] - ta1_y_zzzz_xxx_1[i] * pc_z[i];

        ta1_y_zzzzz_xxy_0[i] = 4.0 * ta1_y_zzz_xxy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxy_1[i] * fe_0 + ta1_y_zzzz_xxy_0[i] * pa_z[i] - ta1_y_zzzz_xxy_1[i] * pc_z[i];

        ta1_y_zzzzz_xxz_0[i] = 4.0 * ta1_y_zzz_xxz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxz_1[i] * fe_0 + ta1_y_zzzz_xx_0[i] * fe_0 - ta1_y_zzzz_xx_1[i] * fe_0 + ta1_y_zzzz_xxz_0[i] * pa_z[i] - ta1_y_zzzz_xxz_1[i] * pc_z[i];

        ta1_y_zzzzz_xyy_0[i] = 4.0 * ta1_y_zzz_xyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyy_1[i] * fe_0 + ta1_y_zzzz_xyy_0[i] * pa_z[i] - ta1_y_zzzz_xyy_1[i] * pc_z[i];

        ta1_y_zzzzz_xyz_0[i] = 4.0 * ta1_y_zzz_xyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyz_1[i] * fe_0 + ta1_y_zzzz_xy_0[i] * fe_0 - ta1_y_zzzz_xy_1[i] * fe_0 + ta1_y_zzzz_xyz_0[i] * pa_z[i] - ta1_y_zzzz_xyz_1[i] * pc_z[i];

        ta1_y_zzzzz_xzz_0[i] = 4.0 * ta1_y_zzz_xzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_xz_1[i] * fe_0 + ta1_y_zzzz_xzz_0[i] * pa_z[i] - ta1_y_zzzz_xzz_1[i] * pc_z[i];

        ta1_y_zzzzz_yyy_0[i] = 4.0 * ta1_y_zzz_yyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyy_1[i] * fe_0 + ta1_y_zzzz_yyy_0[i] * pa_z[i] - ta1_y_zzzz_yyy_1[i] * pc_z[i];

        ta1_y_zzzzz_yyz_0[i] = 4.0 * ta1_y_zzz_yyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyz_1[i] * fe_0 + ta1_y_zzzz_yy_0[i] * fe_0 - ta1_y_zzzz_yy_1[i] * fe_0 + ta1_y_zzzz_yyz_0[i] * pa_z[i] - ta1_y_zzzz_yyz_1[i] * pc_z[i];

        ta1_y_zzzzz_yzz_0[i] = 4.0 * ta1_y_zzz_yzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_y_zzzz_yz_1[i] * fe_0 + ta1_y_zzzz_yzz_0[i] * pa_z[i] - ta1_y_zzzz_yzz_1[i] * pc_z[i];

        ta1_y_zzzzz_zzz_0[i] = 4.0 * ta1_y_zzz_zzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_zzz_1[i] * fe_0 + 3.0 * ta1_y_zzzz_zz_0[i] * fe_0 - 3.0 * ta1_y_zzzz_zz_1[i] * fe_0 + ta1_y_zzzz_zzz_0[i] * pa_z[i] - ta1_y_zzzz_zzz_1[i] * pc_z[i];
    }

    // Set up 420-430 components of targeted buffer : HF

    auto ta1_z_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 420);

    auto ta1_z_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 421);

    auto ta1_z_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 422);

    auto ta1_z_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 423);

    auto ta1_z_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 424);

    auto ta1_z_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 425);

    auto ta1_z_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 426);

    auto ta1_z_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 427);

    auto ta1_z_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 428);

    auto ta1_z_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 429);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xxx_xxx_0, ta1_z_xxx_xxx_1, ta1_z_xxx_xxy_0, ta1_z_xxx_xxy_1, ta1_z_xxx_xxz_0, ta1_z_xxx_xxz_1, ta1_z_xxx_xyy_0, ta1_z_xxx_xyy_1, ta1_z_xxx_xyz_0, ta1_z_xxx_xyz_1, ta1_z_xxx_xzz_0, ta1_z_xxx_xzz_1, ta1_z_xxx_yyy_0, ta1_z_xxx_yyy_1, ta1_z_xxx_yyz_0, ta1_z_xxx_yyz_1, ta1_z_xxx_yzz_0, ta1_z_xxx_yzz_1, ta1_z_xxx_zzz_0, ta1_z_xxx_zzz_1, ta1_z_xxxx_xx_0, ta1_z_xxxx_xx_1, ta1_z_xxxx_xxx_0, ta1_z_xxxx_xxx_1, ta1_z_xxxx_xxy_0, ta1_z_xxxx_xxy_1, ta1_z_xxxx_xxz_0, ta1_z_xxxx_xxz_1, ta1_z_xxxx_xy_0, ta1_z_xxxx_xy_1, ta1_z_xxxx_xyy_0, ta1_z_xxxx_xyy_1, ta1_z_xxxx_xyz_0, ta1_z_xxxx_xyz_1, ta1_z_xxxx_xz_0, ta1_z_xxxx_xz_1, ta1_z_xxxx_xzz_0, ta1_z_xxxx_xzz_1, ta1_z_xxxx_yy_0, ta1_z_xxxx_yy_1, ta1_z_xxxx_yyy_0, ta1_z_xxxx_yyy_1, ta1_z_xxxx_yyz_0, ta1_z_xxxx_yyz_1, ta1_z_xxxx_yz_0, ta1_z_xxxx_yz_1, ta1_z_xxxx_yzz_0, ta1_z_xxxx_yzz_1, ta1_z_xxxx_zz_0, ta1_z_xxxx_zz_1, ta1_z_xxxx_zzz_0, ta1_z_xxxx_zzz_1, ta1_z_xxxxx_xxx_0, ta1_z_xxxxx_xxy_0, ta1_z_xxxxx_xxz_0, ta1_z_xxxxx_xyy_0, ta1_z_xxxxx_xyz_0, ta1_z_xxxxx_xzz_0, ta1_z_xxxxx_yyy_0, ta1_z_xxxxx_yyz_0, ta1_z_xxxxx_yzz_0, ta1_z_xxxxx_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxx_xxx_0[i] = 4.0 * ta1_z_xxx_xxx_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxx_1[i] * fe_0 + 3.0 * ta1_z_xxxx_xx_0[i] * fe_0 - 3.0 * ta1_z_xxxx_xx_1[i] * fe_0 + ta1_z_xxxx_xxx_0[i] * pa_x[i] - ta1_z_xxxx_xxx_1[i] * pc_x[i];

        ta1_z_xxxxx_xxy_0[i] = 4.0 * ta1_z_xxx_xxy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xy_1[i] * fe_0 + ta1_z_xxxx_xxy_0[i] * pa_x[i] - ta1_z_xxxx_xxy_1[i] * pc_x[i];

        ta1_z_xxxxx_xxz_0[i] = 4.0 * ta1_z_xxx_xxz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xz_1[i] * fe_0 + ta1_z_xxxx_xxz_0[i] * pa_x[i] - ta1_z_xxxx_xxz_1[i] * pc_x[i];

        ta1_z_xxxxx_xyy_0[i] = 4.0 * ta1_z_xxx_xyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyy_1[i] * fe_0 + ta1_z_xxxx_yy_0[i] * fe_0 - ta1_z_xxxx_yy_1[i] * fe_0 + ta1_z_xxxx_xyy_0[i] * pa_x[i] - ta1_z_xxxx_xyy_1[i] * pc_x[i];

        ta1_z_xxxxx_xyz_0[i] = 4.0 * ta1_z_xxx_xyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyz_1[i] * fe_0 + ta1_z_xxxx_yz_0[i] * fe_0 - ta1_z_xxxx_yz_1[i] * fe_0 + ta1_z_xxxx_xyz_0[i] * pa_x[i] - ta1_z_xxxx_xyz_1[i] * pc_x[i];

        ta1_z_xxxxx_xzz_0[i] = 4.0 * ta1_z_xxx_xzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xzz_1[i] * fe_0 + ta1_z_xxxx_zz_0[i] * fe_0 - ta1_z_xxxx_zz_1[i] * fe_0 + ta1_z_xxxx_xzz_0[i] * pa_x[i] - ta1_z_xxxx_xzz_1[i] * pc_x[i];

        ta1_z_xxxxx_yyy_0[i] = 4.0 * ta1_z_xxx_yyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_yyy_1[i] * fe_0 + ta1_z_xxxx_yyy_0[i] * pa_x[i] - ta1_z_xxxx_yyy_1[i] * pc_x[i];

        ta1_z_xxxxx_yyz_0[i] = 4.0 * ta1_z_xxx_yyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yyz_1[i] * fe_0 + ta1_z_xxxx_yyz_0[i] * pa_x[i] - ta1_z_xxxx_yyz_1[i] * pc_x[i];

        ta1_z_xxxxx_yzz_0[i] = 4.0 * ta1_z_xxx_yzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yzz_1[i] * fe_0 + ta1_z_xxxx_yzz_0[i] * pa_x[i] - ta1_z_xxxx_yzz_1[i] * pc_x[i];

        ta1_z_xxxxx_zzz_0[i] = 4.0 * ta1_z_xxx_zzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_zzz_1[i] * fe_0 + ta1_z_xxxx_zzz_0[i] * pa_x[i] - ta1_z_xxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 430-440 components of targeted buffer : HF

    auto ta1_z_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 430);

    auto ta1_z_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 431);

    auto ta1_z_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 432);

    auto ta1_z_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 433);

    auto ta1_z_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 434);

    auto ta1_z_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 435);

    auto ta1_z_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 436);

    auto ta1_z_xxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 437);

    auto ta1_z_xxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 438);

    auto ta1_z_xxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 439);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxx_xx_0, ta1_z_xxxx_xx_1, ta1_z_xxxx_xxx_0, ta1_z_xxxx_xxx_1, ta1_z_xxxx_xxy_0, ta1_z_xxxx_xxy_1, ta1_z_xxxx_xxz_0, ta1_z_xxxx_xxz_1, ta1_z_xxxx_xy_0, ta1_z_xxxx_xy_1, ta1_z_xxxx_xyy_0, ta1_z_xxxx_xyy_1, ta1_z_xxxx_xyz_0, ta1_z_xxxx_xyz_1, ta1_z_xxxx_xz_0, ta1_z_xxxx_xz_1, ta1_z_xxxx_xzz_0, ta1_z_xxxx_xzz_1, ta1_z_xxxx_zzz_0, ta1_z_xxxx_zzz_1, ta1_z_xxxxy_xxx_0, ta1_z_xxxxy_xxy_0, ta1_z_xxxxy_xxz_0, ta1_z_xxxxy_xyy_0, ta1_z_xxxxy_xyz_0, ta1_z_xxxxy_xzz_0, ta1_z_xxxxy_yyy_0, ta1_z_xxxxy_yyz_0, ta1_z_xxxxy_yzz_0, ta1_z_xxxxy_zzz_0, ta1_z_xxxy_yyy_0, ta1_z_xxxy_yyy_1, ta1_z_xxxy_yyz_0, ta1_z_xxxy_yyz_1, ta1_z_xxxy_yzz_0, ta1_z_xxxy_yzz_1, ta1_z_xxy_yyy_0, ta1_z_xxy_yyy_1, ta1_z_xxy_yyz_0, ta1_z_xxy_yyz_1, ta1_z_xxy_yzz_0, ta1_z_xxy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxy_xxx_0[i] = ta1_z_xxxx_xxx_0[i] * pa_y[i] - ta1_z_xxxx_xxx_1[i] * pc_y[i];

        ta1_z_xxxxy_xxy_0[i] = ta1_z_xxxx_xx_0[i] * fe_0 - ta1_z_xxxx_xx_1[i] * fe_0 + ta1_z_xxxx_xxy_0[i] * pa_y[i] - ta1_z_xxxx_xxy_1[i] * pc_y[i];

        ta1_z_xxxxy_xxz_0[i] = ta1_z_xxxx_xxz_0[i] * pa_y[i] - ta1_z_xxxx_xxz_1[i] * pc_y[i];

        ta1_z_xxxxy_xyy_0[i] = 2.0 * ta1_z_xxxx_xy_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xy_1[i] * fe_0 + ta1_z_xxxx_xyy_0[i] * pa_y[i] - ta1_z_xxxx_xyy_1[i] * pc_y[i];

        ta1_z_xxxxy_xyz_0[i] = ta1_z_xxxx_xz_0[i] * fe_0 - ta1_z_xxxx_xz_1[i] * fe_0 + ta1_z_xxxx_xyz_0[i] * pa_y[i] - ta1_z_xxxx_xyz_1[i] * pc_y[i];

        ta1_z_xxxxy_xzz_0[i] = ta1_z_xxxx_xzz_0[i] * pa_y[i] - ta1_z_xxxx_xzz_1[i] * pc_y[i];

        ta1_z_xxxxy_yyy_0[i] = 3.0 * ta1_z_xxy_yyy_0[i] * fe_0 - 3.0 * ta1_z_xxy_yyy_1[i] * fe_0 + ta1_z_xxxy_yyy_0[i] * pa_x[i] - ta1_z_xxxy_yyy_1[i] * pc_x[i];

        ta1_z_xxxxy_yyz_0[i] = 3.0 * ta1_z_xxy_yyz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yyz_1[i] * fe_0 + ta1_z_xxxy_yyz_0[i] * pa_x[i] - ta1_z_xxxy_yyz_1[i] * pc_x[i];

        ta1_z_xxxxy_yzz_0[i] = 3.0 * ta1_z_xxy_yzz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yzz_1[i] * fe_0 + ta1_z_xxxy_yzz_0[i] * pa_x[i] - ta1_z_xxxy_yzz_1[i] * pc_x[i];

        ta1_z_xxxxy_zzz_0[i] = ta1_z_xxxx_zzz_0[i] * pa_y[i] - ta1_z_xxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 440-450 components of targeted buffer : HF

    auto ta1_z_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 440);

    auto ta1_z_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 441);

    auto ta1_z_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 442);

    auto ta1_z_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 443);

    auto ta1_z_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 444);

    auto ta1_z_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 445);

    auto ta1_z_xxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 446);

    auto ta1_z_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 447);

    auto ta1_z_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 448);

    auto ta1_z_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 449);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxx_xx_0, ta1_z_xxxx_xx_1, ta1_z_xxxx_xxx_0, ta1_z_xxxx_xxx_1, ta1_z_xxxx_xxy_0, ta1_z_xxxx_xxy_1, ta1_z_xxxx_xxz_0, ta1_z_xxxx_xxz_1, ta1_z_xxxx_xy_0, ta1_z_xxxx_xy_1, ta1_z_xxxx_xyy_0, ta1_z_xxxx_xyy_1, ta1_z_xxxx_xyz_0, ta1_z_xxxx_xyz_1, ta1_z_xxxx_xz_0, ta1_z_xxxx_xz_1, ta1_z_xxxx_xzz_0, ta1_z_xxxx_xzz_1, ta1_z_xxxx_yyy_0, ta1_z_xxxx_yyy_1, ta1_z_xxxxz_xxx_0, ta1_z_xxxxz_xxy_0, ta1_z_xxxxz_xxz_0, ta1_z_xxxxz_xyy_0, ta1_z_xxxxz_xyz_0, ta1_z_xxxxz_xzz_0, ta1_z_xxxxz_yyy_0, ta1_z_xxxxz_yyz_0, ta1_z_xxxxz_yzz_0, ta1_z_xxxxz_zzz_0, ta1_z_xxxz_yyz_0, ta1_z_xxxz_yyz_1, ta1_z_xxxz_yzz_0, ta1_z_xxxz_yzz_1, ta1_z_xxxz_zzz_0, ta1_z_xxxz_zzz_1, ta1_z_xxz_yyz_0, ta1_z_xxz_yyz_1, ta1_z_xxz_yzz_0, ta1_z_xxz_yzz_1, ta1_z_xxz_zzz_0, ta1_z_xxz_zzz_1, ta_xxxx_xxx_1, ta_xxxx_xxy_1, ta_xxxx_xxz_1, ta_xxxx_xyy_1, ta_xxxx_xyz_1, ta_xxxx_xzz_1, ta_xxxx_yyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxz_xxx_0[i] = ta_xxxx_xxx_1[i] + ta1_z_xxxx_xxx_0[i] * pa_z[i] - ta1_z_xxxx_xxx_1[i] * pc_z[i];

        ta1_z_xxxxz_xxy_0[i] = ta_xxxx_xxy_1[i] + ta1_z_xxxx_xxy_0[i] * pa_z[i] - ta1_z_xxxx_xxy_1[i] * pc_z[i];

        ta1_z_xxxxz_xxz_0[i] = ta1_z_xxxx_xx_0[i] * fe_0 - ta1_z_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxz_1[i] + ta1_z_xxxx_xxz_0[i] * pa_z[i] - ta1_z_xxxx_xxz_1[i] * pc_z[i];

        ta1_z_xxxxz_xyy_0[i] = ta_xxxx_xyy_1[i] + ta1_z_xxxx_xyy_0[i] * pa_z[i] - ta1_z_xxxx_xyy_1[i] * pc_z[i];

        ta1_z_xxxxz_xyz_0[i] = ta1_z_xxxx_xy_0[i] * fe_0 - ta1_z_xxxx_xy_1[i] * fe_0 + ta_xxxx_xyz_1[i] + ta1_z_xxxx_xyz_0[i] * pa_z[i] - ta1_z_xxxx_xyz_1[i] * pc_z[i];

        ta1_z_xxxxz_xzz_0[i] = 2.0 * ta1_z_xxxx_xz_0[i] * fe_0 - 2.0 * ta1_z_xxxx_xz_1[i] * fe_0 + ta_xxxx_xzz_1[i] + ta1_z_xxxx_xzz_0[i] * pa_z[i] - ta1_z_xxxx_xzz_1[i] * pc_z[i];

        ta1_z_xxxxz_yyy_0[i] = ta_xxxx_yyy_1[i] + ta1_z_xxxx_yyy_0[i] * pa_z[i] - ta1_z_xxxx_yyy_1[i] * pc_z[i];

        ta1_z_xxxxz_yyz_0[i] = 3.0 * ta1_z_xxz_yyz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yyz_1[i] * fe_0 + ta1_z_xxxz_yyz_0[i] * pa_x[i] - ta1_z_xxxz_yyz_1[i] * pc_x[i];

        ta1_z_xxxxz_yzz_0[i] = 3.0 * ta1_z_xxz_yzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yzz_1[i] * fe_0 + ta1_z_xxxz_yzz_0[i] * pa_x[i] - ta1_z_xxxz_yzz_1[i] * pc_x[i];

        ta1_z_xxxxz_zzz_0[i] = 3.0 * ta1_z_xxz_zzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_zzz_1[i] * fe_0 + ta1_z_xxxz_zzz_0[i] * pa_x[i] - ta1_z_xxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 450-460 components of targeted buffer : HF

    auto ta1_z_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 450);

    auto ta1_z_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 451);

    auto ta1_z_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 452);

    auto ta1_z_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 453);

    auto ta1_z_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 454);

    auto ta1_z_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 455);

    auto ta1_z_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 456);

    auto ta1_z_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 457);

    auto ta1_z_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 458);

    auto ta1_z_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 459);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxx_xxx_0, ta1_z_xxx_xxx_1, ta1_z_xxx_xxz_0, ta1_z_xxx_xxz_1, ta1_z_xxx_xzz_0, ta1_z_xxx_xzz_1, ta1_z_xxxy_xxx_0, ta1_z_xxxy_xxx_1, ta1_z_xxxy_xxz_0, ta1_z_xxxy_xxz_1, ta1_z_xxxy_xzz_0, ta1_z_xxxy_xzz_1, ta1_z_xxxyy_xxx_0, ta1_z_xxxyy_xxy_0, ta1_z_xxxyy_xxz_0, ta1_z_xxxyy_xyy_0, ta1_z_xxxyy_xyz_0, ta1_z_xxxyy_xzz_0, ta1_z_xxxyy_yyy_0, ta1_z_xxxyy_yyz_0, ta1_z_xxxyy_yzz_0, ta1_z_xxxyy_zzz_0, ta1_z_xxyy_xxy_0, ta1_z_xxyy_xxy_1, ta1_z_xxyy_xy_0, ta1_z_xxyy_xy_1, ta1_z_xxyy_xyy_0, ta1_z_xxyy_xyy_1, ta1_z_xxyy_xyz_0, ta1_z_xxyy_xyz_1, ta1_z_xxyy_yy_0, ta1_z_xxyy_yy_1, ta1_z_xxyy_yyy_0, ta1_z_xxyy_yyy_1, ta1_z_xxyy_yyz_0, ta1_z_xxyy_yyz_1, ta1_z_xxyy_yz_0, ta1_z_xxyy_yz_1, ta1_z_xxyy_yzz_0, ta1_z_xxyy_yzz_1, ta1_z_xxyy_zzz_0, ta1_z_xxyy_zzz_1, ta1_z_xyy_xxy_0, ta1_z_xyy_xxy_1, ta1_z_xyy_xyy_0, ta1_z_xyy_xyy_1, ta1_z_xyy_xyz_0, ta1_z_xyy_xyz_1, ta1_z_xyy_yyy_0, ta1_z_xyy_yyy_1, ta1_z_xyy_yyz_0, ta1_z_xyy_yyz_1, ta1_z_xyy_yzz_0, ta1_z_xyy_yzz_1, ta1_z_xyy_zzz_0, ta1_z_xyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyy_xxx_0[i] = ta1_z_xxx_xxx_0[i] * fe_0 - ta1_z_xxx_xxx_1[i] * fe_0 + ta1_z_xxxy_xxx_0[i] * pa_y[i] - ta1_z_xxxy_xxx_1[i] * pc_y[i];

        ta1_z_xxxyy_xxy_0[i] = 2.0 * ta1_z_xyy_xxy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxyy_xy_0[i] * fe_0 - 2.0 * ta1_z_xxyy_xy_1[i] * fe_0 + ta1_z_xxyy_xxy_0[i] * pa_x[i] - ta1_z_xxyy_xxy_1[i] * pc_x[i];

        ta1_z_xxxyy_xxz_0[i] = ta1_z_xxx_xxz_0[i] * fe_0 - ta1_z_xxx_xxz_1[i] * fe_0 + ta1_z_xxxy_xxz_0[i] * pa_y[i] - ta1_z_xxxy_xxz_1[i] * pc_y[i];

        ta1_z_xxxyy_xyy_0[i] = 2.0 * ta1_z_xyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xyy_1[i] * fe_0 + ta1_z_xxyy_yy_0[i] * fe_0 - ta1_z_xxyy_yy_1[i] * fe_0 + ta1_z_xxyy_xyy_0[i] * pa_x[i] - ta1_z_xxyy_xyy_1[i] * pc_x[i];

        ta1_z_xxxyy_xyz_0[i] = 2.0 * ta1_z_xyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xyy_xyz_1[i] * fe_0 + ta1_z_xxyy_yz_0[i] * fe_0 - ta1_z_xxyy_yz_1[i] * fe_0 + ta1_z_xxyy_xyz_0[i] * pa_x[i] - ta1_z_xxyy_xyz_1[i] * pc_x[i];

        ta1_z_xxxyy_xzz_0[i] = ta1_z_xxx_xzz_0[i] * fe_0 - ta1_z_xxx_xzz_1[i] * fe_0 + ta1_z_xxxy_xzz_0[i] * pa_y[i] - ta1_z_xxxy_xzz_1[i] * pc_y[i];

        ta1_z_xxxyy_yyy_0[i] = 2.0 * ta1_z_xyy_yyy_0[i] * fe_0 - 2.0 * ta1_z_xyy_yyy_1[i] * fe_0 + ta1_z_xxyy_yyy_0[i] * pa_x[i] - ta1_z_xxyy_yyy_1[i] * pc_x[i];

        ta1_z_xxxyy_yyz_0[i] = 2.0 * ta1_z_xyy_yyz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yyz_1[i] * fe_0 + ta1_z_xxyy_yyz_0[i] * pa_x[i] - ta1_z_xxyy_yyz_1[i] * pc_x[i];

        ta1_z_xxxyy_yzz_0[i] = 2.0 * ta1_z_xyy_yzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yzz_1[i] * fe_0 + ta1_z_xxyy_yzz_0[i] * pa_x[i] - ta1_z_xxyy_yzz_1[i] * pc_x[i];

        ta1_z_xxxyy_zzz_0[i] = 2.0 * ta1_z_xyy_zzz_0[i] * fe_0 - 2.0 * ta1_z_xyy_zzz_1[i] * fe_0 + ta1_z_xxyy_zzz_0[i] * pa_x[i] - ta1_z_xxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 460-470 components of targeted buffer : HF

    auto ta1_z_xxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 460);

    auto ta1_z_xxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 461);

    auto ta1_z_xxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 462);

    auto ta1_z_xxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 463);

    auto ta1_z_xxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 464);

    auto ta1_z_xxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 465);

    auto ta1_z_xxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 466);

    auto ta1_z_xxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 467);

    auto ta1_z_xxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 468);

    auto ta1_z_xxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 469);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxxy_xxy_0, ta1_z_xxxy_xxy_1, ta1_z_xxxy_xyy_0, ta1_z_xxxy_xyy_1, ta1_z_xxxy_yyy_0, ta1_z_xxxy_yyy_1, ta1_z_xxxyz_xxx_0, ta1_z_xxxyz_xxy_0, ta1_z_xxxyz_xxz_0, ta1_z_xxxyz_xyy_0, ta1_z_xxxyz_xyz_0, ta1_z_xxxyz_xzz_0, ta1_z_xxxyz_yyy_0, ta1_z_xxxyz_yyz_0, ta1_z_xxxyz_yzz_0, ta1_z_xxxyz_zzz_0, ta1_z_xxxz_xxx_0, ta1_z_xxxz_xxx_1, ta1_z_xxxz_xxz_0, ta1_z_xxxz_xxz_1, ta1_z_xxxz_xyz_0, ta1_z_xxxz_xyz_1, ta1_z_xxxz_xz_0, ta1_z_xxxz_xz_1, ta1_z_xxxz_xzz_0, ta1_z_xxxz_xzz_1, ta1_z_xxxz_zzz_0, ta1_z_xxxz_zzz_1, ta1_z_xxyz_yyz_0, ta1_z_xxyz_yyz_1, ta1_z_xxyz_yzz_0, ta1_z_xxyz_yzz_1, ta1_z_xyz_yyz_0, ta1_z_xyz_yyz_1, ta1_z_xyz_yzz_0, ta1_z_xyz_yzz_1, ta_xxxy_xxy_1, ta_xxxy_xyy_1, ta_xxxy_yyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyz_xxx_0[i] = ta1_z_xxxz_xxx_0[i] * pa_y[i] - ta1_z_xxxz_xxx_1[i] * pc_y[i];

        ta1_z_xxxyz_xxy_0[i] = ta_xxxy_xxy_1[i] + ta1_z_xxxy_xxy_0[i] * pa_z[i] - ta1_z_xxxy_xxy_1[i] * pc_z[i];

        ta1_z_xxxyz_xxz_0[i] = ta1_z_xxxz_xxz_0[i] * pa_y[i] - ta1_z_xxxz_xxz_1[i] * pc_y[i];

        ta1_z_xxxyz_xyy_0[i] = ta_xxxy_xyy_1[i] + ta1_z_xxxy_xyy_0[i] * pa_z[i] - ta1_z_xxxy_xyy_1[i] * pc_z[i];

        ta1_z_xxxyz_xyz_0[i] = ta1_z_xxxz_xz_0[i] * fe_0 - ta1_z_xxxz_xz_1[i] * fe_0 + ta1_z_xxxz_xyz_0[i] * pa_y[i] - ta1_z_xxxz_xyz_1[i] * pc_y[i];

        ta1_z_xxxyz_xzz_0[i] = ta1_z_xxxz_xzz_0[i] * pa_y[i] - ta1_z_xxxz_xzz_1[i] * pc_y[i];

        ta1_z_xxxyz_yyy_0[i] = ta_xxxy_yyy_1[i] + ta1_z_xxxy_yyy_0[i] * pa_z[i] - ta1_z_xxxy_yyy_1[i] * pc_z[i];

        ta1_z_xxxyz_yyz_0[i] = 2.0 * ta1_z_xyz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yyz_1[i] * fe_0 + ta1_z_xxyz_yyz_0[i] * pa_x[i] - ta1_z_xxyz_yyz_1[i] * pc_x[i];

        ta1_z_xxxyz_yzz_0[i] = 2.0 * ta1_z_xyz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yzz_1[i] * fe_0 + ta1_z_xxyz_yzz_0[i] * pa_x[i] - ta1_z_xxyz_yzz_1[i] * pc_x[i];

        ta1_z_xxxyz_zzz_0[i] = ta1_z_xxxz_zzz_0[i] * pa_y[i] - ta1_z_xxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 470-480 components of targeted buffer : HF

    auto ta1_z_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 470);

    auto ta1_z_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 471);

    auto ta1_z_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 472);

    auto ta1_z_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 473);

    auto ta1_z_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 474);

    auto ta1_z_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 475);

    auto ta1_z_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 476);

    auto ta1_z_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 477);

    auto ta1_z_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 478);

    auto ta1_z_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 479);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxx_xxx_0, ta1_z_xxx_xxx_1, ta1_z_xxx_xxy_0, ta1_z_xxx_xxy_1, ta1_z_xxx_xyy_0, ta1_z_xxx_xyy_1, ta1_z_xxxz_xxx_0, ta1_z_xxxz_xxx_1, ta1_z_xxxz_xxy_0, ta1_z_xxxz_xxy_1, ta1_z_xxxz_xyy_0, ta1_z_xxxz_xyy_1, ta1_z_xxxzz_xxx_0, ta1_z_xxxzz_xxy_0, ta1_z_xxxzz_xxz_0, ta1_z_xxxzz_xyy_0, ta1_z_xxxzz_xyz_0, ta1_z_xxxzz_xzz_0, ta1_z_xxxzz_yyy_0, ta1_z_xxxzz_yyz_0, ta1_z_xxxzz_yzz_0, ta1_z_xxxzz_zzz_0, ta1_z_xxzz_xxz_0, ta1_z_xxzz_xxz_1, ta1_z_xxzz_xyz_0, ta1_z_xxzz_xyz_1, ta1_z_xxzz_xz_0, ta1_z_xxzz_xz_1, ta1_z_xxzz_xzz_0, ta1_z_xxzz_xzz_1, ta1_z_xxzz_yyy_0, ta1_z_xxzz_yyy_1, ta1_z_xxzz_yyz_0, ta1_z_xxzz_yyz_1, ta1_z_xxzz_yz_0, ta1_z_xxzz_yz_1, ta1_z_xxzz_yzz_0, ta1_z_xxzz_yzz_1, ta1_z_xxzz_zz_0, ta1_z_xxzz_zz_1, ta1_z_xxzz_zzz_0, ta1_z_xxzz_zzz_1, ta1_z_xzz_xxz_0, ta1_z_xzz_xxz_1, ta1_z_xzz_xyz_0, ta1_z_xzz_xyz_1, ta1_z_xzz_xzz_0, ta1_z_xzz_xzz_1, ta1_z_xzz_yyy_0, ta1_z_xzz_yyy_1, ta1_z_xzz_yyz_0, ta1_z_xzz_yyz_1, ta1_z_xzz_yzz_0, ta1_z_xzz_yzz_1, ta1_z_xzz_zzz_0, ta1_z_xzz_zzz_1, ta_xxxz_xxx_1, ta_xxxz_xxy_1, ta_xxxz_xyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzz_xxx_0[i] = ta1_z_xxx_xxx_0[i] * fe_0 - ta1_z_xxx_xxx_1[i] * fe_0 + ta_xxxz_xxx_1[i] + ta1_z_xxxz_xxx_0[i] * pa_z[i] - ta1_z_xxxz_xxx_1[i] * pc_z[i];

        ta1_z_xxxzz_xxy_0[i] = ta1_z_xxx_xxy_0[i] * fe_0 - ta1_z_xxx_xxy_1[i] * fe_0 + ta_xxxz_xxy_1[i] + ta1_z_xxxz_xxy_0[i] * pa_z[i] - ta1_z_xxxz_xxy_1[i] * pc_z[i];

        ta1_z_xxxzz_xxz_0[i] = 2.0 * ta1_z_xzz_xxz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxzz_xz_0[i] * fe_0 - 2.0 * ta1_z_xxzz_xz_1[i] * fe_0 + ta1_z_xxzz_xxz_0[i] * pa_x[i] - ta1_z_xxzz_xxz_1[i] * pc_x[i];

        ta1_z_xxxzz_xyy_0[i] = ta1_z_xxx_xyy_0[i] * fe_0 - ta1_z_xxx_xyy_1[i] * fe_0 + ta_xxxz_xyy_1[i] + ta1_z_xxxz_xyy_0[i] * pa_z[i] - ta1_z_xxxz_xyy_1[i] * pc_z[i];

        ta1_z_xxxzz_xyz_0[i] = 2.0 * ta1_z_xzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xyz_1[i] * fe_0 + ta1_z_xxzz_yz_0[i] * fe_0 - ta1_z_xxzz_yz_1[i] * fe_0 + ta1_z_xxzz_xyz_0[i] * pa_x[i] - ta1_z_xxzz_xyz_1[i] * pc_x[i];

        ta1_z_xxxzz_xzz_0[i] = 2.0 * ta1_z_xzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xzz_1[i] * fe_0 + ta1_z_xxzz_zz_0[i] * fe_0 - ta1_z_xxzz_zz_1[i] * fe_0 + ta1_z_xxzz_xzz_0[i] * pa_x[i] - ta1_z_xxzz_xzz_1[i] * pc_x[i];

        ta1_z_xxxzz_yyy_0[i] = 2.0 * ta1_z_xzz_yyy_0[i] * fe_0 - 2.0 * ta1_z_xzz_yyy_1[i] * fe_0 + ta1_z_xxzz_yyy_0[i] * pa_x[i] - ta1_z_xxzz_yyy_1[i] * pc_x[i];

        ta1_z_xxxzz_yyz_0[i] = 2.0 * ta1_z_xzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yyz_1[i] * fe_0 + ta1_z_xxzz_yyz_0[i] * pa_x[i] - ta1_z_xxzz_yyz_1[i] * pc_x[i];

        ta1_z_xxxzz_yzz_0[i] = 2.0 * ta1_z_xzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yzz_1[i] * fe_0 + ta1_z_xxzz_yzz_0[i] * pa_x[i] - ta1_z_xxzz_yzz_1[i] * pc_x[i];

        ta1_z_xxxzz_zzz_0[i] = 2.0 * ta1_z_xzz_zzz_0[i] * fe_0 - 2.0 * ta1_z_xzz_zzz_1[i] * fe_0 + ta1_z_xxzz_zzz_0[i] * pa_x[i] - ta1_z_xxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 480-490 components of targeted buffer : HF

    auto ta1_z_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 480);

    auto ta1_z_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 481);

    auto ta1_z_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 482);

    auto ta1_z_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 483);

    auto ta1_z_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 484);

    auto ta1_z_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 485);

    auto ta1_z_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 486);

    auto ta1_z_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 487);

    auto ta1_z_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 488);

    auto ta1_z_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 489);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxy_xxx_0, ta1_z_xxy_xxx_1, ta1_z_xxy_xxz_0, ta1_z_xxy_xxz_1, ta1_z_xxy_xzz_0, ta1_z_xxy_xzz_1, ta1_z_xxyy_xxx_0, ta1_z_xxyy_xxx_1, ta1_z_xxyy_xxz_0, ta1_z_xxyy_xxz_1, ta1_z_xxyy_xzz_0, ta1_z_xxyy_xzz_1, ta1_z_xxyyy_xxx_0, ta1_z_xxyyy_xxy_0, ta1_z_xxyyy_xxz_0, ta1_z_xxyyy_xyy_0, ta1_z_xxyyy_xyz_0, ta1_z_xxyyy_xzz_0, ta1_z_xxyyy_yyy_0, ta1_z_xxyyy_yyz_0, ta1_z_xxyyy_yzz_0, ta1_z_xxyyy_zzz_0, ta1_z_xyyy_xxy_0, ta1_z_xyyy_xxy_1, ta1_z_xyyy_xy_0, ta1_z_xyyy_xy_1, ta1_z_xyyy_xyy_0, ta1_z_xyyy_xyy_1, ta1_z_xyyy_xyz_0, ta1_z_xyyy_xyz_1, ta1_z_xyyy_yy_0, ta1_z_xyyy_yy_1, ta1_z_xyyy_yyy_0, ta1_z_xyyy_yyy_1, ta1_z_xyyy_yyz_0, ta1_z_xyyy_yyz_1, ta1_z_xyyy_yz_0, ta1_z_xyyy_yz_1, ta1_z_xyyy_yzz_0, ta1_z_xyyy_yzz_1, ta1_z_xyyy_zzz_0, ta1_z_xyyy_zzz_1, ta1_z_yyy_xxy_0, ta1_z_yyy_xxy_1, ta1_z_yyy_xyy_0, ta1_z_yyy_xyy_1, ta1_z_yyy_xyz_0, ta1_z_yyy_xyz_1, ta1_z_yyy_yyy_0, ta1_z_yyy_yyy_1, ta1_z_yyy_yyz_0, ta1_z_yyy_yyz_1, ta1_z_yyy_yzz_0, ta1_z_yyy_yzz_1, ta1_z_yyy_zzz_0, ta1_z_yyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyy_xxx_0[i] = 2.0 * ta1_z_xxy_xxx_0[i] * fe_0 - 2.0 * ta1_z_xxy_xxx_1[i] * fe_0 + ta1_z_xxyy_xxx_0[i] * pa_y[i] - ta1_z_xxyy_xxx_1[i] * pc_y[i];

        ta1_z_xxyyy_xxy_0[i] = ta1_z_yyy_xxy_0[i] * fe_0 - ta1_z_yyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xyyy_xy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xy_1[i] * fe_0 + ta1_z_xyyy_xxy_0[i] * pa_x[i] - ta1_z_xyyy_xxy_1[i] * pc_x[i];

        ta1_z_xxyyy_xxz_0[i] = 2.0 * ta1_z_xxy_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xxz_1[i] * fe_0 + ta1_z_xxyy_xxz_0[i] * pa_y[i] - ta1_z_xxyy_xxz_1[i] * pc_y[i];

        ta1_z_xxyyy_xyy_0[i] = ta1_z_yyy_xyy_0[i] * fe_0 - ta1_z_yyy_xyy_1[i] * fe_0 + ta1_z_xyyy_yy_0[i] * fe_0 - ta1_z_xyyy_yy_1[i] * fe_0 + ta1_z_xyyy_xyy_0[i] * pa_x[i] - ta1_z_xyyy_xyy_1[i] * pc_x[i];

        ta1_z_xxyyy_xyz_0[i] = ta1_z_yyy_xyz_0[i] * fe_0 - ta1_z_yyy_xyz_1[i] * fe_0 + ta1_z_xyyy_yz_0[i] * fe_0 - ta1_z_xyyy_yz_1[i] * fe_0 + ta1_z_xyyy_xyz_0[i] * pa_x[i] - ta1_z_xyyy_xyz_1[i] * pc_x[i];

        ta1_z_xxyyy_xzz_0[i] = 2.0 * ta1_z_xxy_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xzz_1[i] * fe_0 + ta1_z_xxyy_xzz_0[i] * pa_y[i] - ta1_z_xxyy_xzz_1[i] * pc_y[i];

        ta1_z_xxyyy_yyy_0[i] = ta1_z_yyy_yyy_0[i] * fe_0 - ta1_z_yyy_yyy_1[i] * fe_0 + ta1_z_xyyy_yyy_0[i] * pa_x[i] - ta1_z_xyyy_yyy_1[i] * pc_x[i];

        ta1_z_xxyyy_yyz_0[i] = ta1_z_yyy_yyz_0[i] * fe_0 - ta1_z_yyy_yyz_1[i] * fe_0 + ta1_z_xyyy_yyz_0[i] * pa_x[i] - ta1_z_xyyy_yyz_1[i] * pc_x[i];

        ta1_z_xxyyy_yzz_0[i] = ta1_z_yyy_yzz_0[i] * fe_0 - ta1_z_yyy_yzz_1[i] * fe_0 + ta1_z_xyyy_yzz_0[i] * pa_x[i] - ta1_z_xyyy_yzz_1[i] * pc_x[i];

        ta1_z_xxyyy_zzz_0[i] = ta1_z_yyy_zzz_0[i] * fe_0 - ta1_z_yyy_zzz_1[i] * fe_0 + ta1_z_xyyy_zzz_0[i] * pa_x[i] - ta1_z_xyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 490-500 components of targeted buffer : HF

    auto ta1_z_xxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 490);

    auto ta1_z_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 491);

    auto ta1_z_xxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 492);

    auto ta1_z_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 493);

    auto ta1_z_xxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 494);

    auto ta1_z_xxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 495);

    auto ta1_z_xxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 496);

    auto ta1_z_xxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 497);

    auto ta1_z_xxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 498);

    auto ta1_z_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 499);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xxyy_xxx_0, ta1_z_xxyy_xxx_1, ta1_z_xxyy_xxy_0, ta1_z_xxyy_xxy_1, ta1_z_xxyy_xy_0, ta1_z_xxyy_xy_1, ta1_z_xxyy_xyy_0, ta1_z_xxyy_xyy_1, ta1_z_xxyy_xyz_0, ta1_z_xxyy_xyz_1, ta1_z_xxyy_yyy_0, ta1_z_xxyy_yyy_1, ta1_z_xxyyz_xxx_0, ta1_z_xxyyz_xxy_0, ta1_z_xxyyz_xxz_0, ta1_z_xxyyz_xyy_0, ta1_z_xxyyz_xyz_0, ta1_z_xxyyz_xzz_0, ta1_z_xxyyz_yyy_0, ta1_z_xxyyz_yyz_0, ta1_z_xxyyz_yzz_0, ta1_z_xxyyz_zzz_0, ta1_z_xxyz_xxz_0, ta1_z_xxyz_xxz_1, ta1_z_xxyz_xzz_0, ta1_z_xxyz_xzz_1, ta1_z_xxz_xxz_0, ta1_z_xxz_xxz_1, ta1_z_xxz_xzz_0, ta1_z_xxz_xzz_1, ta1_z_xyyz_yyz_0, ta1_z_xyyz_yyz_1, ta1_z_xyyz_yzz_0, ta1_z_xyyz_yzz_1, ta1_z_xyyz_zzz_0, ta1_z_xyyz_zzz_1, ta1_z_yyz_yyz_0, ta1_z_yyz_yyz_1, ta1_z_yyz_yzz_0, ta1_z_yyz_yzz_1, ta1_z_yyz_zzz_0, ta1_z_yyz_zzz_1, ta_xxyy_xxx_1, ta_xxyy_xxy_1, ta_xxyy_xyy_1, ta_xxyy_xyz_1, ta_xxyy_yyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyz_xxx_0[i] = ta_xxyy_xxx_1[i] + ta1_z_xxyy_xxx_0[i] * pa_z[i] - ta1_z_xxyy_xxx_1[i] * pc_z[i];

        ta1_z_xxyyz_xxy_0[i] = ta_xxyy_xxy_1[i] + ta1_z_xxyy_xxy_0[i] * pa_z[i] - ta1_z_xxyy_xxy_1[i] * pc_z[i];

        ta1_z_xxyyz_xxz_0[i] = ta1_z_xxz_xxz_0[i] * fe_0 - ta1_z_xxz_xxz_1[i] * fe_0 + ta1_z_xxyz_xxz_0[i] * pa_y[i] - ta1_z_xxyz_xxz_1[i] * pc_y[i];

        ta1_z_xxyyz_xyy_0[i] = ta_xxyy_xyy_1[i] + ta1_z_xxyy_xyy_0[i] * pa_z[i] - ta1_z_xxyy_xyy_1[i] * pc_z[i];

        ta1_z_xxyyz_xyz_0[i] = ta1_z_xxyy_xy_0[i] * fe_0 - ta1_z_xxyy_xy_1[i] * fe_0 + ta_xxyy_xyz_1[i] + ta1_z_xxyy_xyz_0[i] * pa_z[i] - ta1_z_xxyy_xyz_1[i] * pc_z[i];

        ta1_z_xxyyz_xzz_0[i] = ta1_z_xxz_xzz_0[i] * fe_0 - ta1_z_xxz_xzz_1[i] * fe_0 + ta1_z_xxyz_xzz_0[i] * pa_y[i] - ta1_z_xxyz_xzz_1[i] * pc_y[i];

        ta1_z_xxyyz_yyy_0[i] = ta_xxyy_yyy_1[i] + ta1_z_xxyy_yyy_0[i] * pa_z[i] - ta1_z_xxyy_yyy_1[i] * pc_z[i];

        ta1_z_xxyyz_yyz_0[i] = ta1_z_yyz_yyz_0[i] * fe_0 - ta1_z_yyz_yyz_1[i] * fe_0 + ta1_z_xyyz_yyz_0[i] * pa_x[i] - ta1_z_xyyz_yyz_1[i] * pc_x[i];

        ta1_z_xxyyz_yzz_0[i] = ta1_z_yyz_yzz_0[i] * fe_0 - ta1_z_yyz_yzz_1[i] * fe_0 + ta1_z_xyyz_yzz_0[i] * pa_x[i] - ta1_z_xyyz_yzz_1[i] * pc_x[i];

        ta1_z_xxyyz_zzz_0[i] = ta1_z_yyz_zzz_0[i] * fe_0 - ta1_z_yyz_zzz_1[i] * fe_0 + ta1_z_xyyz_zzz_0[i] * pa_x[i] - ta1_z_xyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 500-510 components of targeted buffer : HF

    auto ta1_z_xxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 500);

    auto ta1_z_xxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 501);

    auto ta1_z_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 502);

    auto ta1_z_xxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 503);

    auto ta1_z_xxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 504);

    auto ta1_z_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 505);

    auto ta1_z_xxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 506);

    auto ta1_z_xxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 507);

    auto ta1_z_xxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 508);

    auto ta1_z_xxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 509);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyzz_xxx_0, ta1_z_xxyzz_xxy_0, ta1_z_xxyzz_xxz_0, ta1_z_xxyzz_xyy_0, ta1_z_xxyzz_xyz_0, ta1_z_xxyzz_xzz_0, ta1_z_xxyzz_yyy_0, ta1_z_xxyzz_yyz_0, ta1_z_xxyzz_yzz_0, ta1_z_xxyzz_zzz_0, ta1_z_xxzz_xx_0, ta1_z_xxzz_xx_1, ta1_z_xxzz_xxx_0, ta1_z_xxzz_xxx_1, ta1_z_xxzz_xxy_0, ta1_z_xxzz_xxy_1, ta1_z_xxzz_xxz_0, ta1_z_xxzz_xxz_1, ta1_z_xxzz_xy_0, ta1_z_xxzz_xy_1, ta1_z_xxzz_xyy_0, ta1_z_xxzz_xyy_1, ta1_z_xxzz_xyz_0, ta1_z_xxzz_xyz_1, ta1_z_xxzz_xz_0, ta1_z_xxzz_xz_1, ta1_z_xxzz_xzz_0, ta1_z_xxzz_xzz_1, ta1_z_xxzz_zzz_0, ta1_z_xxzz_zzz_1, ta1_z_xyzz_yyy_0, ta1_z_xyzz_yyy_1, ta1_z_xyzz_yyz_0, ta1_z_xyzz_yyz_1, ta1_z_xyzz_yzz_0, ta1_z_xyzz_yzz_1, ta1_z_yzz_yyy_0, ta1_z_yzz_yyy_1, ta1_z_yzz_yyz_0, ta1_z_yzz_yyz_1, ta1_z_yzz_yzz_0, ta1_z_yzz_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzz_xxx_0[i] = ta1_z_xxzz_xxx_0[i] * pa_y[i] - ta1_z_xxzz_xxx_1[i] * pc_y[i];

        ta1_z_xxyzz_xxy_0[i] = ta1_z_xxzz_xx_0[i] * fe_0 - ta1_z_xxzz_xx_1[i] * fe_0 + ta1_z_xxzz_xxy_0[i] * pa_y[i] - ta1_z_xxzz_xxy_1[i] * pc_y[i];

        ta1_z_xxyzz_xxz_0[i] = ta1_z_xxzz_xxz_0[i] * pa_y[i] - ta1_z_xxzz_xxz_1[i] * pc_y[i];

        ta1_z_xxyzz_xyy_0[i] = 2.0 * ta1_z_xxzz_xy_0[i] * fe_0 - 2.0 * ta1_z_xxzz_xy_1[i] * fe_0 + ta1_z_xxzz_xyy_0[i] * pa_y[i] - ta1_z_xxzz_xyy_1[i] * pc_y[i];

        ta1_z_xxyzz_xyz_0[i] = ta1_z_xxzz_xz_0[i] * fe_0 - ta1_z_xxzz_xz_1[i] * fe_0 + ta1_z_xxzz_xyz_0[i] * pa_y[i] - ta1_z_xxzz_xyz_1[i] * pc_y[i];

        ta1_z_xxyzz_xzz_0[i] = ta1_z_xxzz_xzz_0[i] * pa_y[i] - ta1_z_xxzz_xzz_1[i] * pc_y[i];

        ta1_z_xxyzz_yyy_0[i] = ta1_z_yzz_yyy_0[i] * fe_0 - ta1_z_yzz_yyy_1[i] * fe_0 + ta1_z_xyzz_yyy_0[i] * pa_x[i] - ta1_z_xyzz_yyy_1[i] * pc_x[i];

        ta1_z_xxyzz_yyz_0[i] = ta1_z_yzz_yyz_0[i] * fe_0 - ta1_z_yzz_yyz_1[i] * fe_0 + ta1_z_xyzz_yyz_0[i] * pa_x[i] - ta1_z_xyzz_yyz_1[i] * pc_x[i];

        ta1_z_xxyzz_yzz_0[i] = ta1_z_yzz_yzz_0[i] * fe_0 - ta1_z_yzz_yzz_1[i] * fe_0 + ta1_z_xyzz_yzz_0[i] * pa_x[i] - ta1_z_xyzz_yzz_1[i] * pc_x[i];

        ta1_z_xxyzz_zzz_0[i] = ta1_z_xxzz_zzz_0[i] * pa_y[i] - ta1_z_xxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 510-520 components of targeted buffer : HF

    auto ta1_z_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 510);

    auto ta1_z_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 511);

    auto ta1_z_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 512);

    auto ta1_z_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 513);

    auto ta1_z_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 514);

    auto ta1_z_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 515);

    auto ta1_z_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 516);

    auto ta1_z_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 517);

    auto ta1_z_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 518);

    auto ta1_z_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 519);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxz_xxx_0, ta1_z_xxz_xxx_1, ta1_z_xxz_xxy_0, ta1_z_xxz_xxy_1, ta1_z_xxz_xyy_0, ta1_z_xxz_xyy_1, ta1_z_xxzz_xxx_0, ta1_z_xxzz_xxx_1, ta1_z_xxzz_xxy_0, ta1_z_xxzz_xxy_1, ta1_z_xxzz_xyy_0, ta1_z_xxzz_xyy_1, ta1_z_xxzzz_xxx_0, ta1_z_xxzzz_xxy_0, ta1_z_xxzzz_xxz_0, ta1_z_xxzzz_xyy_0, ta1_z_xxzzz_xyz_0, ta1_z_xxzzz_xzz_0, ta1_z_xxzzz_yyy_0, ta1_z_xxzzz_yyz_0, ta1_z_xxzzz_yzz_0, ta1_z_xxzzz_zzz_0, ta1_z_xzzz_xxz_0, ta1_z_xzzz_xxz_1, ta1_z_xzzz_xyz_0, ta1_z_xzzz_xyz_1, ta1_z_xzzz_xz_0, ta1_z_xzzz_xz_1, ta1_z_xzzz_xzz_0, ta1_z_xzzz_xzz_1, ta1_z_xzzz_yyy_0, ta1_z_xzzz_yyy_1, ta1_z_xzzz_yyz_0, ta1_z_xzzz_yyz_1, ta1_z_xzzz_yz_0, ta1_z_xzzz_yz_1, ta1_z_xzzz_yzz_0, ta1_z_xzzz_yzz_1, ta1_z_xzzz_zz_0, ta1_z_xzzz_zz_1, ta1_z_xzzz_zzz_0, ta1_z_xzzz_zzz_1, ta1_z_zzz_xxz_0, ta1_z_zzz_xxz_1, ta1_z_zzz_xyz_0, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_0, ta1_z_zzz_xzz_1, ta1_z_zzz_yyy_0, ta1_z_zzz_yyy_1, ta1_z_zzz_yyz_0, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_0, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_0, ta1_z_zzz_zzz_1, ta_xxzz_xxx_1, ta_xxzz_xxy_1, ta_xxzz_xyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzz_xxx_0[i] = 2.0 * ta1_z_xxz_xxx_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxx_1[i] * fe_0 + ta_xxzz_xxx_1[i] + ta1_z_xxzz_xxx_0[i] * pa_z[i] - ta1_z_xxzz_xxx_1[i] * pc_z[i];

        ta1_z_xxzzz_xxy_0[i] = 2.0 * ta1_z_xxz_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxy_1[i] * fe_0 + ta_xxzz_xxy_1[i] + ta1_z_xxzz_xxy_0[i] * pa_z[i] - ta1_z_xxzz_xxy_1[i] * pc_z[i];

        ta1_z_xxzzz_xxz_0[i] = ta1_z_zzz_xxz_0[i] * fe_0 - ta1_z_zzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xz_1[i] * fe_0 + ta1_z_xzzz_xxz_0[i] * pa_x[i] - ta1_z_xzzz_xxz_1[i] * pc_x[i];

        ta1_z_xxzzz_xyy_0[i] = 2.0 * ta1_z_xxz_xyy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xyy_1[i] * fe_0 + ta_xxzz_xyy_1[i] + ta1_z_xxzz_xyy_0[i] * pa_z[i] - ta1_z_xxzz_xyy_1[i] * pc_z[i];

        ta1_z_xxzzz_xyz_0[i] = ta1_z_zzz_xyz_0[i] * fe_0 - ta1_z_zzz_xyz_1[i] * fe_0 + ta1_z_xzzz_yz_0[i] * fe_0 - ta1_z_xzzz_yz_1[i] * fe_0 + ta1_z_xzzz_xyz_0[i] * pa_x[i] - ta1_z_xzzz_xyz_1[i] * pc_x[i];

        ta1_z_xxzzz_xzz_0[i] = ta1_z_zzz_xzz_0[i] * fe_0 - ta1_z_zzz_xzz_1[i] * fe_0 + ta1_z_xzzz_zz_0[i] * fe_0 - ta1_z_xzzz_zz_1[i] * fe_0 + ta1_z_xzzz_xzz_0[i] * pa_x[i] - ta1_z_xzzz_xzz_1[i] * pc_x[i];

        ta1_z_xxzzz_yyy_0[i] = ta1_z_zzz_yyy_0[i] * fe_0 - ta1_z_zzz_yyy_1[i] * fe_0 + ta1_z_xzzz_yyy_0[i] * pa_x[i] - ta1_z_xzzz_yyy_1[i] * pc_x[i];

        ta1_z_xxzzz_yyz_0[i] = ta1_z_zzz_yyz_0[i] * fe_0 - ta1_z_zzz_yyz_1[i] * fe_0 + ta1_z_xzzz_yyz_0[i] * pa_x[i] - ta1_z_xzzz_yyz_1[i] * pc_x[i];

        ta1_z_xxzzz_yzz_0[i] = ta1_z_zzz_yzz_0[i] * fe_0 - ta1_z_zzz_yzz_1[i] * fe_0 + ta1_z_xzzz_yzz_0[i] * pa_x[i] - ta1_z_xzzz_yzz_1[i] * pc_x[i];

        ta1_z_xxzzz_zzz_0[i] = ta1_z_zzz_zzz_0[i] * fe_0 - ta1_z_zzz_zzz_1[i] * fe_0 + ta1_z_xzzz_zzz_0[i] * pa_x[i] - ta1_z_xzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 520-530 components of targeted buffer : HF

    auto ta1_z_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 520);

    auto ta1_z_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 521);

    auto ta1_z_xyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 522);

    auto ta1_z_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 523);

    auto ta1_z_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 524);

    auto ta1_z_xyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 525);

    auto ta1_z_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 526);

    auto ta1_z_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 527);

    auto ta1_z_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 528);

    auto ta1_z_xyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 529);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyyy_xxx_0, ta1_z_xyyyy_xxy_0, ta1_z_xyyyy_xxz_0, ta1_z_xyyyy_xyy_0, ta1_z_xyyyy_xyz_0, ta1_z_xyyyy_xzz_0, ta1_z_xyyyy_yyy_0, ta1_z_xyyyy_yyz_0, ta1_z_xyyyy_yzz_0, ta1_z_xyyyy_zzz_0, ta1_z_yyyy_xx_0, ta1_z_yyyy_xx_1, ta1_z_yyyy_xxx_0, ta1_z_yyyy_xxx_1, ta1_z_yyyy_xxy_0, ta1_z_yyyy_xxy_1, ta1_z_yyyy_xxz_0, ta1_z_yyyy_xxz_1, ta1_z_yyyy_xy_0, ta1_z_yyyy_xy_1, ta1_z_yyyy_xyy_0, ta1_z_yyyy_xyy_1, ta1_z_yyyy_xyz_0, ta1_z_yyyy_xyz_1, ta1_z_yyyy_xz_0, ta1_z_yyyy_xz_1, ta1_z_yyyy_xzz_0, ta1_z_yyyy_xzz_1, ta1_z_yyyy_yy_0, ta1_z_yyyy_yy_1, ta1_z_yyyy_yyy_0, ta1_z_yyyy_yyy_1, ta1_z_yyyy_yyz_0, ta1_z_yyyy_yyz_1, ta1_z_yyyy_yz_0, ta1_z_yyyy_yz_1, ta1_z_yyyy_yzz_0, ta1_z_yyyy_yzz_1, ta1_z_yyyy_zz_0, ta1_z_yyyy_zz_1, ta1_z_yyyy_zzz_0, ta1_z_yyyy_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyy_xxx_0[i] = 3.0 * ta1_z_yyyy_xx_0[i] * fe_0 - 3.0 * ta1_z_yyyy_xx_1[i] * fe_0 + ta1_z_yyyy_xxx_0[i] * pa_x[i] - ta1_z_yyyy_xxx_1[i] * pc_x[i];

        ta1_z_xyyyy_xxy_0[i] = 2.0 * ta1_z_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xy_1[i] * fe_0 + ta1_z_yyyy_xxy_0[i] * pa_x[i] - ta1_z_yyyy_xxy_1[i] * pc_x[i];

        ta1_z_xyyyy_xxz_0[i] = 2.0 * ta1_z_yyyy_xz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xz_1[i] * fe_0 + ta1_z_yyyy_xxz_0[i] * pa_x[i] - ta1_z_yyyy_xxz_1[i] * pc_x[i];

        ta1_z_xyyyy_xyy_0[i] = ta1_z_yyyy_yy_0[i] * fe_0 - ta1_z_yyyy_yy_1[i] * fe_0 + ta1_z_yyyy_xyy_0[i] * pa_x[i] - ta1_z_yyyy_xyy_1[i] * pc_x[i];

        ta1_z_xyyyy_xyz_0[i] = ta1_z_yyyy_yz_0[i] * fe_0 - ta1_z_yyyy_yz_1[i] * fe_0 + ta1_z_yyyy_xyz_0[i] * pa_x[i] - ta1_z_yyyy_xyz_1[i] * pc_x[i];

        ta1_z_xyyyy_xzz_0[i] = ta1_z_yyyy_zz_0[i] * fe_0 - ta1_z_yyyy_zz_1[i] * fe_0 + ta1_z_yyyy_xzz_0[i] * pa_x[i] - ta1_z_yyyy_xzz_1[i] * pc_x[i];

        ta1_z_xyyyy_yyy_0[i] = ta1_z_yyyy_yyy_0[i] * pa_x[i] - ta1_z_yyyy_yyy_1[i] * pc_x[i];

        ta1_z_xyyyy_yyz_0[i] = ta1_z_yyyy_yyz_0[i] * pa_x[i] - ta1_z_yyyy_yyz_1[i] * pc_x[i];

        ta1_z_xyyyy_yzz_0[i] = ta1_z_yyyy_yzz_0[i] * pa_x[i] - ta1_z_yyyy_yzz_1[i] * pc_x[i];

        ta1_z_xyyyy_zzz_0[i] = ta1_z_yyyy_zzz_0[i] * pa_x[i] - ta1_z_yyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 530-540 components of targeted buffer : HF

    auto ta1_z_xyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 530);

    auto ta1_z_xyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 531);

    auto ta1_z_xyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 532);

    auto ta1_z_xyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 533);

    auto ta1_z_xyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 534);

    auto ta1_z_xyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 535);

    auto ta1_z_xyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 536);

    auto ta1_z_xyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 537);

    auto ta1_z_xyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 538);

    auto ta1_z_xyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 539);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xyyy_xxx_0, ta1_z_xyyy_xxx_1, ta1_z_xyyy_xxy_0, ta1_z_xyyy_xxy_1, ta1_z_xyyy_xyy_0, ta1_z_xyyy_xyy_1, ta1_z_xyyyz_xxx_0, ta1_z_xyyyz_xxy_0, ta1_z_xyyyz_xxz_0, ta1_z_xyyyz_xyy_0, ta1_z_xyyyz_xyz_0, ta1_z_xyyyz_xzz_0, ta1_z_xyyyz_yyy_0, ta1_z_xyyyz_yyz_0, ta1_z_xyyyz_yzz_0, ta1_z_xyyyz_zzz_0, ta1_z_yyyz_xxz_0, ta1_z_yyyz_xxz_1, ta1_z_yyyz_xyz_0, ta1_z_yyyz_xyz_1, ta1_z_yyyz_xz_0, ta1_z_yyyz_xz_1, ta1_z_yyyz_xzz_0, ta1_z_yyyz_xzz_1, ta1_z_yyyz_yyy_0, ta1_z_yyyz_yyy_1, ta1_z_yyyz_yyz_0, ta1_z_yyyz_yyz_1, ta1_z_yyyz_yz_0, ta1_z_yyyz_yz_1, ta1_z_yyyz_yzz_0, ta1_z_yyyz_yzz_1, ta1_z_yyyz_zz_0, ta1_z_yyyz_zz_1, ta1_z_yyyz_zzz_0, ta1_z_yyyz_zzz_1, ta_xyyy_xxx_1, ta_xyyy_xxy_1, ta_xyyy_xyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyz_xxx_0[i] = ta_xyyy_xxx_1[i] + ta1_z_xyyy_xxx_0[i] * pa_z[i] - ta1_z_xyyy_xxx_1[i] * pc_z[i];

        ta1_z_xyyyz_xxy_0[i] = ta_xyyy_xxy_1[i] + ta1_z_xyyy_xxy_0[i] * pa_z[i] - ta1_z_xyyy_xxy_1[i] * pc_z[i];

        ta1_z_xyyyz_xxz_0[i] = 2.0 * ta1_z_yyyz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xz_1[i] * fe_0 + ta1_z_yyyz_xxz_0[i] * pa_x[i] - ta1_z_yyyz_xxz_1[i] * pc_x[i];

        ta1_z_xyyyz_xyy_0[i] = ta_xyyy_xyy_1[i] + ta1_z_xyyy_xyy_0[i] * pa_z[i] - ta1_z_xyyy_xyy_1[i] * pc_z[i];

        ta1_z_xyyyz_xyz_0[i] = ta1_z_yyyz_yz_0[i] * fe_0 - ta1_z_yyyz_yz_1[i] * fe_0 + ta1_z_yyyz_xyz_0[i] * pa_x[i] - ta1_z_yyyz_xyz_1[i] * pc_x[i];

        ta1_z_xyyyz_xzz_0[i] = ta1_z_yyyz_zz_0[i] * fe_0 - ta1_z_yyyz_zz_1[i] * fe_0 + ta1_z_yyyz_xzz_0[i] * pa_x[i] - ta1_z_yyyz_xzz_1[i] * pc_x[i];

        ta1_z_xyyyz_yyy_0[i] = ta1_z_yyyz_yyy_0[i] * pa_x[i] - ta1_z_yyyz_yyy_1[i] * pc_x[i];

        ta1_z_xyyyz_yyz_0[i] = ta1_z_yyyz_yyz_0[i] * pa_x[i] - ta1_z_yyyz_yyz_1[i] * pc_x[i];

        ta1_z_xyyyz_yzz_0[i] = ta1_z_yyyz_yzz_0[i] * pa_x[i] - ta1_z_yyyz_yzz_1[i] * pc_x[i];

        ta1_z_xyyyz_zzz_0[i] = ta1_z_yyyz_zzz_0[i] * pa_x[i] - ta1_z_yyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 540-550 components of targeted buffer : HF

    auto ta1_z_xyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 540);

    auto ta1_z_xyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 541);

    auto ta1_z_xyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 542);

    auto ta1_z_xyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 543);

    auto ta1_z_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 544);

    auto ta1_z_xyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 545);

    auto ta1_z_xyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 546);

    auto ta1_z_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 547);

    auto ta1_z_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 548);

    auto ta1_z_xyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 549);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyzz_xxx_0, ta1_z_xyyzz_xxy_0, ta1_z_xyyzz_xxz_0, ta1_z_xyyzz_xyy_0, ta1_z_xyyzz_xyz_0, ta1_z_xyyzz_xzz_0, ta1_z_xyyzz_yyy_0, ta1_z_xyyzz_yyz_0, ta1_z_xyyzz_yzz_0, ta1_z_xyyzz_zzz_0, ta1_z_yyzz_xx_0, ta1_z_yyzz_xx_1, ta1_z_yyzz_xxx_0, ta1_z_yyzz_xxx_1, ta1_z_yyzz_xxy_0, ta1_z_yyzz_xxy_1, ta1_z_yyzz_xxz_0, ta1_z_yyzz_xxz_1, ta1_z_yyzz_xy_0, ta1_z_yyzz_xy_1, ta1_z_yyzz_xyy_0, ta1_z_yyzz_xyy_1, ta1_z_yyzz_xyz_0, ta1_z_yyzz_xyz_1, ta1_z_yyzz_xz_0, ta1_z_yyzz_xz_1, ta1_z_yyzz_xzz_0, ta1_z_yyzz_xzz_1, ta1_z_yyzz_yy_0, ta1_z_yyzz_yy_1, ta1_z_yyzz_yyy_0, ta1_z_yyzz_yyy_1, ta1_z_yyzz_yyz_0, ta1_z_yyzz_yyz_1, ta1_z_yyzz_yz_0, ta1_z_yyzz_yz_1, ta1_z_yyzz_yzz_0, ta1_z_yyzz_yzz_1, ta1_z_yyzz_zz_0, ta1_z_yyzz_zz_1, ta1_z_yyzz_zzz_0, ta1_z_yyzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzz_xxx_0[i] = 3.0 * ta1_z_yyzz_xx_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xx_1[i] * fe_0 + ta1_z_yyzz_xxx_0[i] * pa_x[i] - ta1_z_yyzz_xxx_1[i] * pc_x[i];

        ta1_z_xyyzz_xxy_0[i] = 2.0 * ta1_z_yyzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yyzz_xy_1[i] * fe_0 + ta1_z_yyzz_xxy_0[i] * pa_x[i] - ta1_z_yyzz_xxy_1[i] * pc_x[i];

        ta1_z_xyyzz_xxz_0[i] = 2.0 * ta1_z_yyzz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyzz_xz_1[i] * fe_0 + ta1_z_yyzz_xxz_0[i] * pa_x[i] - ta1_z_yyzz_xxz_1[i] * pc_x[i];

        ta1_z_xyyzz_xyy_0[i] = ta1_z_yyzz_yy_0[i] * fe_0 - ta1_z_yyzz_yy_1[i] * fe_0 + ta1_z_yyzz_xyy_0[i] * pa_x[i] - ta1_z_yyzz_xyy_1[i] * pc_x[i];

        ta1_z_xyyzz_xyz_0[i] = ta1_z_yyzz_yz_0[i] * fe_0 - ta1_z_yyzz_yz_1[i] * fe_0 + ta1_z_yyzz_xyz_0[i] * pa_x[i] - ta1_z_yyzz_xyz_1[i] * pc_x[i];

        ta1_z_xyyzz_xzz_0[i] = ta1_z_yyzz_zz_0[i] * fe_0 - ta1_z_yyzz_zz_1[i] * fe_0 + ta1_z_yyzz_xzz_0[i] * pa_x[i] - ta1_z_yyzz_xzz_1[i] * pc_x[i];

        ta1_z_xyyzz_yyy_0[i] = ta1_z_yyzz_yyy_0[i] * pa_x[i] - ta1_z_yyzz_yyy_1[i] * pc_x[i];

        ta1_z_xyyzz_yyz_0[i] = ta1_z_yyzz_yyz_0[i] * pa_x[i] - ta1_z_yyzz_yyz_1[i] * pc_x[i];

        ta1_z_xyyzz_yzz_0[i] = ta1_z_yyzz_yzz_0[i] * pa_x[i] - ta1_z_yyzz_yzz_1[i] * pc_x[i];

        ta1_z_xyyzz_zzz_0[i] = ta1_z_yyzz_zzz_0[i] * pa_x[i] - ta1_z_yyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 550-560 components of targeted buffer : HF

    auto ta1_z_xyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 550);

    auto ta1_z_xyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 551);

    auto ta1_z_xyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 552);

    auto ta1_z_xyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 553);

    auto ta1_z_xyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 554);

    auto ta1_z_xyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 555);

    auto ta1_z_xyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 556);

    auto ta1_z_xyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 557);

    auto ta1_z_xyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 558);

    auto ta1_z_xyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 559);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xyzzz_xxx_0, ta1_z_xyzzz_xxy_0, ta1_z_xyzzz_xxz_0, ta1_z_xyzzz_xyy_0, ta1_z_xyzzz_xyz_0, ta1_z_xyzzz_xzz_0, ta1_z_xyzzz_yyy_0, ta1_z_xyzzz_yyz_0, ta1_z_xyzzz_yzz_0, ta1_z_xyzzz_zzz_0, ta1_z_xzzz_xxx_0, ta1_z_xzzz_xxx_1, ta1_z_xzzz_xxz_0, ta1_z_xzzz_xxz_1, ta1_z_xzzz_xzz_0, ta1_z_xzzz_xzz_1, ta1_z_yzzz_xxy_0, ta1_z_yzzz_xxy_1, ta1_z_yzzz_xy_0, ta1_z_yzzz_xy_1, ta1_z_yzzz_xyy_0, ta1_z_yzzz_xyy_1, ta1_z_yzzz_xyz_0, ta1_z_yzzz_xyz_1, ta1_z_yzzz_yy_0, ta1_z_yzzz_yy_1, ta1_z_yzzz_yyy_0, ta1_z_yzzz_yyy_1, ta1_z_yzzz_yyz_0, ta1_z_yzzz_yyz_1, ta1_z_yzzz_yz_0, ta1_z_yzzz_yz_1, ta1_z_yzzz_yzz_0, ta1_z_yzzz_yzz_1, ta1_z_yzzz_zzz_0, ta1_z_yzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzz_xxx_0[i] = ta1_z_xzzz_xxx_0[i] * pa_y[i] - ta1_z_xzzz_xxx_1[i] * pc_y[i];

        ta1_z_xyzzz_xxy_0[i] = 2.0 * ta1_z_yzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xy_1[i] * fe_0 + ta1_z_yzzz_xxy_0[i] * pa_x[i] - ta1_z_yzzz_xxy_1[i] * pc_x[i];

        ta1_z_xyzzz_xxz_0[i] = ta1_z_xzzz_xxz_0[i] * pa_y[i] - ta1_z_xzzz_xxz_1[i] * pc_y[i];

        ta1_z_xyzzz_xyy_0[i] = ta1_z_yzzz_yy_0[i] * fe_0 - ta1_z_yzzz_yy_1[i] * fe_0 + ta1_z_yzzz_xyy_0[i] * pa_x[i] - ta1_z_yzzz_xyy_1[i] * pc_x[i];

        ta1_z_xyzzz_xyz_0[i] = ta1_z_yzzz_yz_0[i] * fe_0 - ta1_z_yzzz_yz_1[i] * fe_0 + ta1_z_yzzz_xyz_0[i] * pa_x[i] - ta1_z_yzzz_xyz_1[i] * pc_x[i];

        ta1_z_xyzzz_xzz_0[i] = ta1_z_xzzz_xzz_0[i] * pa_y[i] - ta1_z_xzzz_xzz_1[i] * pc_y[i];

        ta1_z_xyzzz_yyy_0[i] = ta1_z_yzzz_yyy_0[i] * pa_x[i] - ta1_z_yzzz_yyy_1[i] * pc_x[i];

        ta1_z_xyzzz_yyz_0[i] = ta1_z_yzzz_yyz_0[i] * pa_x[i] - ta1_z_yzzz_yyz_1[i] * pc_x[i];

        ta1_z_xyzzz_yzz_0[i] = ta1_z_yzzz_yzz_0[i] * pa_x[i] - ta1_z_yzzz_yzz_1[i] * pc_x[i];

        ta1_z_xyzzz_zzz_0[i] = ta1_z_yzzz_zzz_0[i] * pa_x[i] - ta1_z_yzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 560-570 components of targeted buffer : HF

    auto ta1_z_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 560);

    auto ta1_z_xzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 561);

    auto ta1_z_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 562);

    auto ta1_z_xzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 563);

    auto ta1_z_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 564);

    auto ta1_z_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 565);

    auto ta1_z_xzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 566);

    auto ta1_z_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 567);

    auto ta1_z_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 568);

    auto ta1_z_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 569);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xzzzz_xxx_0, ta1_z_xzzzz_xxy_0, ta1_z_xzzzz_xxz_0, ta1_z_xzzzz_xyy_0, ta1_z_xzzzz_xyz_0, ta1_z_xzzzz_xzz_0, ta1_z_xzzzz_yyy_0, ta1_z_xzzzz_yyz_0, ta1_z_xzzzz_yzz_0, ta1_z_xzzzz_zzz_0, ta1_z_zzzz_xx_0, ta1_z_zzzz_xx_1, ta1_z_zzzz_xxx_0, ta1_z_zzzz_xxx_1, ta1_z_zzzz_xxy_0, ta1_z_zzzz_xxy_1, ta1_z_zzzz_xxz_0, ta1_z_zzzz_xxz_1, ta1_z_zzzz_xy_0, ta1_z_zzzz_xy_1, ta1_z_zzzz_xyy_0, ta1_z_zzzz_xyy_1, ta1_z_zzzz_xyz_0, ta1_z_zzzz_xyz_1, ta1_z_zzzz_xz_0, ta1_z_zzzz_xz_1, ta1_z_zzzz_xzz_0, ta1_z_zzzz_xzz_1, ta1_z_zzzz_yy_0, ta1_z_zzzz_yy_1, ta1_z_zzzz_yyy_0, ta1_z_zzzz_yyy_1, ta1_z_zzzz_yyz_0, ta1_z_zzzz_yyz_1, ta1_z_zzzz_yz_0, ta1_z_zzzz_yz_1, ta1_z_zzzz_yzz_0, ta1_z_zzzz_yzz_1, ta1_z_zzzz_zz_0, ta1_z_zzzz_zz_1, ta1_z_zzzz_zzz_0, ta1_z_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzz_xxx_0[i] = 3.0 * ta1_z_zzzz_xx_0[i] * fe_0 - 3.0 * ta1_z_zzzz_xx_1[i] * fe_0 + ta1_z_zzzz_xxx_0[i] * pa_x[i] - ta1_z_zzzz_xxx_1[i] * pc_x[i];

        ta1_z_xzzzz_xxy_0[i] = 2.0 * ta1_z_zzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xy_1[i] * fe_0 + ta1_z_zzzz_xxy_0[i] * pa_x[i] - ta1_z_zzzz_xxy_1[i] * pc_x[i];

        ta1_z_xzzzz_xxz_0[i] = 2.0 * ta1_z_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xz_1[i] * fe_0 + ta1_z_zzzz_xxz_0[i] * pa_x[i] - ta1_z_zzzz_xxz_1[i] * pc_x[i];

        ta1_z_xzzzz_xyy_0[i] = ta1_z_zzzz_yy_0[i] * fe_0 - ta1_z_zzzz_yy_1[i] * fe_0 + ta1_z_zzzz_xyy_0[i] * pa_x[i] - ta1_z_zzzz_xyy_1[i] * pc_x[i];

        ta1_z_xzzzz_xyz_0[i] = ta1_z_zzzz_yz_0[i] * fe_0 - ta1_z_zzzz_yz_1[i] * fe_0 + ta1_z_zzzz_xyz_0[i] * pa_x[i] - ta1_z_zzzz_xyz_1[i] * pc_x[i];

        ta1_z_xzzzz_xzz_0[i] = ta1_z_zzzz_zz_0[i] * fe_0 - ta1_z_zzzz_zz_1[i] * fe_0 + ta1_z_zzzz_xzz_0[i] * pa_x[i] - ta1_z_zzzz_xzz_1[i] * pc_x[i];

        ta1_z_xzzzz_yyy_0[i] = ta1_z_zzzz_yyy_0[i] * pa_x[i] - ta1_z_zzzz_yyy_1[i] * pc_x[i];

        ta1_z_xzzzz_yyz_0[i] = ta1_z_zzzz_yyz_0[i] * pa_x[i] - ta1_z_zzzz_yyz_1[i] * pc_x[i];

        ta1_z_xzzzz_yzz_0[i] = ta1_z_zzzz_yzz_0[i] * pa_x[i] - ta1_z_zzzz_yzz_1[i] * pc_x[i];

        ta1_z_xzzzz_zzz_0[i] = ta1_z_zzzz_zzz_0[i] * pa_x[i] - ta1_z_zzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 570-580 components of targeted buffer : HF

    auto ta1_z_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 570);

    auto ta1_z_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 571);

    auto ta1_z_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 572);

    auto ta1_z_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 573);

    auto ta1_z_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 574);

    auto ta1_z_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 575);

    auto ta1_z_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 576);

    auto ta1_z_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 577);

    auto ta1_z_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 578);

    auto ta1_z_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 579);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yyy_xxx_0, ta1_z_yyy_xxx_1, ta1_z_yyy_xxy_0, ta1_z_yyy_xxy_1, ta1_z_yyy_xxz_0, ta1_z_yyy_xxz_1, ta1_z_yyy_xyy_0, ta1_z_yyy_xyy_1, ta1_z_yyy_xyz_0, ta1_z_yyy_xyz_1, ta1_z_yyy_xzz_0, ta1_z_yyy_xzz_1, ta1_z_yyy_yyy_0, ta1_z_yyy_yyy_1, ta1_z_yyy_yyz_0, ta1_z_yyy_yyz_1, ta1_z_yyy_yzz_0, ta1_z_yyy_yzz_1, ta1_z_yyy_zzz_0, ta1_z_yyy_zzz_1, ta1_z_yyyy_xx_0, ta1_z_yyyy_xx_1, ta1_z_yyyy_xxx_0, ta1_z_yyyy_xxx_1, ta1_z_yyyy_xxy_0, ta1_z_yyyy_xxy_1, ta1_z_yyyy_xxz_0, ta1_z_yyyy_xxz_1, ta1_z_yyyy_xy_0, ta1_z_yyyy_xy_1, ta1_z_yyyy_xyy_0, ta1_z_yyyy_xyy_1, ta1_z_yyyy_xyz_0, ta1_z_yyyy_xyz_1, ta1_z_yyyy_xz_0, ta1_z_yyyy_xz_1, ta1_z_yyyy_xzz_0, ta1_z_yyyy_xzz_1, ta1_z_yyyy_yy_0, ta1_z_yyyy_yy_1, ta1_z_yyyy_yyy_0, ta1_z_yyyy_yyy_1, ta1_z_yyyy_yyz_0, ta1_z_yyyy_yyz_1, ta1_z_yyyy_yz_0, ta1_z_yyyy_yz_1, ta1_z_yyyy_yzz_0, ta1_z_yyyy_yzz_1, ta1_z_yyyy_zz_0, ta1_z_yyyy_zz_1, ta1_z_yyyy_zzz_0, ta1_z_yyyy_zzz_1, ta1_z_yyyyy_xxx_0, ta1_z_yyyyy_xxy_0, ta1_z_yyyyy_xxz_0, ta1_z_yyyyy_xyy_0, ta1_z_yyyyy_xyz_0, ta1_z_yyyyy_xzz_0, ta1_z_yyyyy_yyy_0, ta1_z_yyyyy_yyz_0, ta1_z_yyyyy_yzz_0, ta1_z_yyyyy_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyy_xxx_0[i] = 4.0 * ta1_z_yyy_xxx_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxx_1[i] * fe_0 + ta1_z_yyyy_xxx_0[i] * pa_y[i] - ta1_z_yyyy_xxx_1[i] * pc_y[i];

        ta1_z_yyyyy_xxy_0[i] = 4.0 * ta1_z_yyy_xxy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxy_1[i] * fe_0 + ta1_z_yyyy_xx_0[i] * fe_0 - ta1_z_yyyy_xx_1[i] * fe_0 + ta1_z_yyyy_xxy_0[i] * pa_y[i] - ta1_z_yyyy_xxy_1[i] * pc_y[i];

        ta1_z_yyyyy_xxz_0[i] = 4.0 * ta1_z_yyy_xxz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxz_1[i] * fe_0 + ta1_z_yyyy_xxz_0[i] * pa_y[i] - ta1_z_yyyy_xxz_1[i] * pc_y[i];

        ta1_z_yyyyy_xyy_0[i] = 4.0 * ta1_z_yyy_xyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yyyy_xy_0[i] * fe_0 - 2.0 * ta1_z_yyyy_xy_1[i] * fe_0 + ta1_z_yyyy_xyy_0[i] * pa_y[i] - ta1_z_yyyy_xyy_1[i] * pc_y[i];

        ta1_z_yyyyy_xyz_0[i] = 4.0 * ta1_z_yyy_xyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyz_1[i] * fe_0 + ta1_z_yyyy_xz_0[i] * fe_0 - ta1_z_yyyy_xz_1[i] * fe_0 + ta1_z_yyyy_xyz_0[i] * pa_y[i] - ta1_z_yyyy_xyz_1[i] * pc_y[i];

        ta1_z_yyyyy_xzz_0[i] = 4.0 * ta1_z_yyy_xzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xzz_1[i] * fe_0 + ta1_z_yyyy_xzz_0[i] * pa_y[i] - ta1_z_yyyy_xzz_1[i] * pc_y[i];

        ta1_z_yyyyy_yyy_0[i] = 4.0 * ta1_z_yyy_yyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyy_1[i] * fe_0 + 3.0 * ta1_z_yyyy_yy_0[i] * fe_0 - 3.0 * ta1_z_yyyy_yy_1[i] * fe_0 + ta1_z_yyyy_yyy_0[i] * pa_y[i] - ta1_z_yyyy_yyy_1[i] * pc_y[i];

        ta1_z_yyyyy_yyz_0[i] = 4.0 * ta1_z_yyy_yyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_yz_1[i] * fe_0 + ta1_z_yyyy_yyz_0[i] * pa_y[i] - ta1_z_yyyy_yyz_1[i] * pc_y[i];

        ta1_z_yyyyy_yzz_0[i] = 4.0 * ta1_z_yyy_yzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yzz_1[i] * fe_0 + ta1_z_yyyy_zz_0[i] * fe_0 - ta1_z_yyyy_zz_1[i] * fe_0 + ta1_z_yyyy_yzz_0[i] * pa_y[i] - ta1_z_yyyy_yzz_1[i] * pc_y[i];

        ta1_z_yyyyy_zzz_0[i] = 4.0 * ta1_z_yyy_zzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_zzz_1[i] * fe_0 + ta1_z_yyyy_zzz_0[i] * pa_y[i] - ta1_z_yyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 580-590 components of targeted buffer : HF

    auto ta1_z_yyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 580);

    auto ta1_z_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 581);

    auto ta1_z_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 582);

    auto ta1_z_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 583);

    auto ta1_z_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 584);

    auto ta1_z_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 585);

    auto ta1_z_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 586);

    auto ta1_z_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 587);

    auto ta1_z_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 588);

    auto ta1_z_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 589);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyy_xxx_0, ta1_z_yyyy_xxx_1, ta1_z_yyyy_xxy_0, ta1_z_yyyy_xxy_1, ta1_z_yyyy_xy_0, ta1_z_yyyy_xy_1, ta1_z_yyyy_xyy_0, ta1_z_yyyy_xyy_1, ta1_z_yyyy_xyz_0, ta1_z_yyyy_xyz_1, ta1_z_yyyy_yy_0, ta1_z_yyyy_yy_1, ta1_z_yyyy_yyy_0, ta1_z_yyyy_yyy_1, ta1_z_yyyy_yyz_0, ta1_z_yyyy_yyz_1, ta1_z_yyyy_yz_0, ta1_z_yyyy_yz_1, ta1_z_yyyy_yzz_0, ta1_z_yyyy_yzz_1, ta1_z_yyyyz_xxx_0, ta1_z_yyyyz_xxy_0, ta1_z_yyyyz_xxz_0, ta1_z_yyyyz_xyy_0, ta1_z_yyyyz_xyz_0, ta1_z_yyyyz_xzz_0, ta1_z_yyyyz_yyy_0, ta1_z_yyyyz_yyz_0, ta1_z_yyyyz_yzz_0, ta1_z_yyyyz_zzz_0, ta1_z_yyyz_xxz_0, ta1_z_yyyz_xxz_1, ta1_z_yyyz_xzz_0, ta1_z_yyyz_xzz_1, ta1_z_yyyz_zzz_0, ta1_z_yyyz_zzz_1, ta1_z_yyz_xxz_0, ta1_z_yyz_xxz_1, ta1_z_yyz_xzz_0, ta1_z_yyz_xzz_1, ta1_z_yyz_zzz_0, ta1_z_yyz_zzz_1, ta_yyyy_xxx_1, ta_yyyy_xxy_1, ta_yyyy_xyy_1, ta_yyyy_xyz_1, ta_yyyy_yyy_1, ta_yyyy_yyz_1, ta_yyyy_yzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyz_xxx_0[i] = ta_yyyy_xxx_1[i] + ta1_z_yyyy_xxx_0[i] * pa_z[i] - ta1_z_yyyy_xxx_1[i] * pc_z[i];

        ta1_z_yyyyz_xxy_0[i] = ta_yyyy_xxy_1[i] + ta1_z_yyyy_xxy_0[i] * pa_z[i] - ta1_z_yyyy_xxy_1[i] * pc_z[i];

        ta1_z_yyyyz_xxz_0[i] = 3.0 * ta1_z_yyz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxz_1[i] * fe_0 + ta1_z_yyyz_xxz_0[i] * pa_y[i] - ta1_z_yyyz_xxz_1[i] * pc_y[i];

        ta1_z_yyyyz_xyy_0[i] = ta_yyyy_xyy_1[i] + ta1_z_yyyy_xyy_0[i] * pa_z[i] - ta1_z_yyyy_xyy_1[i] * pc_z[i];

        ta1_z_yyyyz_xyz_0[i] = ta1_z_yyyy_xy_0[i] * fe_0 - ta1_z_yyyy_xy_1[i] * fe_0 + ta_yyyy_xyz_1[i] + ta1_z_yyyy_xyz_0[i] * pa_z[i] - ta1_z_yyyy_xyz_1[i] * pc_z[i];

        ta1_z_yyyyz_xzz_0[i] = 3.0 * ta1_z_yyz_xzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xzz_1[i] * fe_0 + ta1_z_yyyz_xzz_0[i] * pa_y[i] - ta1_z_yyyz_xzz_1[i] * pc_y[i];

        ta1_z_yyyyz_yyy_0[i] = ta_yyyy_yyy_1[i] + ta1_z_yyyy_yyy_0[i] * pa_z[i] - ta1_z_yyyy_yyy_1[i] * pc_z[i];

        ta1_z_yyyyz_yyz_0[i] = ta1_z_yyyy_yy_0[i] * fe_0 - ta1_z_yyyy_yy_1[i] * fe_0 + ta_yyyy_yyz_1[i] + ta1_z_yyyy_yyz_0[i] * pa_z[i] - ta1_z_yyyy_yyz_1[i] * pc_z[i];

        ta1_z_yyyyz_yzz_0[i] = 2.0 * ta1_z_yyyy_yz_0[i] * fe_0 - 2.0 * ta1_z_yyyy_yz_1[i] * fe_0 + ta_yyyy_yzz_1[i] + ta1_z_yyyy_yzz_0[i] * pa_z[i] - ta1_z_yyyy_yzz_1[i] * pc_z[i];

        ta1_z_yyyyz_zzz_0[i] = 3.0 * ta1_z_yyz_zzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_zzz_1[i] * fe_0 + ta1_z_yyyz_zzz_0[i] * pa_y[i] - ta1_z_yyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 590-600 components of targeted buffer : HF

    auto ta1_z_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 590);

    auto ta1_z_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 591);

    auto ta1_z_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 592);

    auto ta1_z_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 593);

    auto ta1_z_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 594);

    auto ta1_z_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 595);

    auto ta1_z_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 596);

    auto ta1_z_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 597);

    auto ta1_z_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 598);

    auto ta1_z_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 599);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyy_xxy_0, ta1_z_yyy_xxy_1, ta1_z_yyy_xyy_0, ta1_z_yyy_xyy_1, ta1_z_yyy_yyy_0, ta1_z_yyy_yyy_1, ta1_z_yyyz_xxy_0, ta1_z_yyyz_xxy_1, ta1_z_yyyz_xyy_0, ta1_z_yyyz_xyy_1, ta1_z_yyyz_yyy_0, ta1_z_yyyz_yyy_1, ta1_z_yyyzz_xxx_0, ta1_z_yyyzz_xxy_0, ta1_z_yyyzz_xxz_0, ta1_z_yyyzz_xyy_0, ta1_z_yyyzz_xyz_0, ta1_z_yyyzz_xzz_0, ta1_z_yyyzz_yyy_0, ta1_z_yyyzz_yyz_0, ta1_z_yyyzz_yzz_0, ta1_z_yyyzz_zzz_0, ta1_z_yyzz_xxx_0, ta1_z_yyzz_xxx_1, ta1_z_yyzz_xxz_0, ta1_z_yyzz_xxz_1, ta1_z_yyzz_xyz_0, ta1_z_yyzz_xyz_1, ta1_z_yyzz_xz_0, ta1_z_yyzz_xz_1, ta1_z_yyzz_xzz_0, ta1_z_yyzz_xzz_1, ta1_z_yyzz_yyz_0, ta1_z_yyzz_yyz_1, ta1_z_yyzz_yz_0, ta1_z_yyzz_yz_1, ta1_z_yyzz_yzz_0, ta1_z_yyzz_yzz_1, ta1_z_yyzz_zz_0, ta1_z_yyzz_zz_1, ta1_z_yyzz_zzz_0, ta1_z_yyzz_zzz_1, ta1_z_yzz_xxx_0, ta1_z_yzz_xxx_1, ta1_z_yzz_xxz_0, ta1_z_yzz_xxz_1, ta1_z_yzz_xyz_0, ta1_z_yzz_xyz_1, ta1_z_yzz_xzz_0, ta1_z_yzz_xzz_1, ta1_z_yzz_yyz_0, ta1_z_yzz_yyz_1, ta1_z_yzz_yzz_0, ta1_z_yzz_yzz_1, ta1_z_yzz_zzz_0, ta1_z_yzz_zzz_1, ta_yyyz_xxy_1, ta_yyyz_xyy_1, ta_yyyz_yyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzz_xxx_0[i] = 2.0 * ta1_z_yzz_xxx_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxx_1[i] * fe_0 + ta1_z_yyzz_xxx_0[i] * pa_y[i] - ta1_z_yyzz_xxx_1[i] * pc_y[i];

        ta1_z_yyyzz_xxy_0[i] = ta1_z_yyy_xxy_0[i] * fe_0 - ta1_z_yyy_xxy_1[i] * fe_0 + ta_yyyz_xxy_1[i] + ta1_z_yyyz_xxy_0[i] * pa_z[i] - ta1_z_yyyz_xxy_1[i] * pc_z[i];

        ta1_z_yyyzz_xxz_0[i] = 2.0 * ta1_z_yzz_xxz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xxz_1[i] * fe_0 + ta1_z_yyzz_xxz_0[i] * pa_y[i] - ta1_z_yyzz_xxz_1[i] * pc_y[i];

        ta1_z_yyyzz_xyy_0[i] = ta1_z_yyy_xyy_0[i] * fe_0 - ta1_z_yyy_xyy_1[i] * fe_0 + ta_yyyz_xyy_1[i] + ta1_z_yyyz_xyy_0[i] * pa_z[i] - ta1_z_yyyz_xyy_1[i] * pc_z[i];

        ta1_z_yyyzz_xyz_0[i] = 2.0 * ta1_z_yzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyz_1[i] * fe_0 + ta1_z_yyzz_xz_0[i] * fe_0 - ta1_z_yyzz_xz_1[i] * fe_0 + ta1_z_yyzz_xyz_0[i] * pa_y[i] - ta1_z_yyzz_xyz_1[i] * pc_y[i];

        ta1_z_yyyzz_xzz_0[i] = 2.0 * ta1_z_yzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xzz_1[i] * fe_0 + ta1_z_yyzz_xzz_0[i] * pa_y[i] - ta1_z_yyzz_xzz_1[i] * pc_y[i];

        ta1_z_yyyzz_yyy_0[i] = ta1_z_yyy_yyy_0[i] * fe_0 - ta1_z_yyy_yyy_1[i] * fe_0 + ta_yyyz_yyy_1[i] + ta1_z_yyyz_yyy_0[i] * pa_z[i] - ta1_z_yyyz_yyy_1[i] * pc_z[i];

        ta1_z_yyyzz_yyz_0[i] = 2.0 * ta1_z_yzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyzz_yz_0[i] * fe_0 - 2.0 * ta1_z_yyzz_yz_1[i] * fe_0 + ta1_z_yyzz_yyz_0[i] * pa_y[i] - ta1_z_yyzz_yyz_1[i] * pc_y[i];

        ta1_z_yyyzz_yzz_0[i] = 2.0 * ta1_z_yzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yzz_1[i] * fe_0 + ta1_z_yyzz_zz_0[i] * fe_0 - ta1_z_yyzz_zz_1[i] * fe_0 + ta1_z_yyzz_yzz_0[i] * pa_y[i] - ta1_z_yyzz_yzz_1[i] * pc_y[i];

        ta1_z_yyyzz_zzz_0[i] = 2.0 * ta1_z_yzz_zzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_zzz_1[i] * fe_0 + ta1_z_yyzz_zzz_0[i] * pa_y[i] - ta1_z_yyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 600-610 components of targeted buffer : HF

    auto ta1_z_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 600);

    auto ta1_z_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 601);

    auto ta1_z_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 602);

    auto ta1_z_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 603);

    auto ta1_z_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 604);

    auto ta1_z_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 605);

    auto ta1_z_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 606);

    auto ta1_z_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 607);

    auto ta1_z_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 608);

    auto ta1_z_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 609);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyz_xxy_0, ta1_z_yyz_xxy_1, ta1_z_yyz_xyy_0, ta1_z_yyz_xyy_1, ta1_z_yyz_yyy_0, ta1_z_yyz_yyy_1, ta1_z_yyzz_xxy_0, ta1_z_yyzz_xxy_1, ta1_z_yyzz_xyy_0, ta1_z_yyzz_xyy_1, ta1_z_yyzz_yyy_0, ta1_z_yyzz_yyy_1, ta1_z_yyzzz_xxx_0, ta1_z_yyzzz_xxy_0, ta1_z_yyzzz_xxz_0, ta1_z_yyzzz_xyy_0, ta1_z_yyzzz_xyz_0, ta1_z_yyzzz_xzz_0, ta1_z_yyzzz_yyy_0, ta1_z_yyzzz_yyz_0, ta1_z_yyzzz_yzz_0, ta1_z_yyzzz_zzz_0, ta1_z_yzzz_xxx_0, ta1_z_yzzz_xxx_1, ta1_z_yzzz_xxz_0, ta1_z_yzzz_xxz_1, ta1_z_yzzz_xyz_0, ta1_z_yzzz_xyz_1, ta1_z_yzzz_xz_0, ta1_z_yzzz_xz_1, ta1_z_yzzz_xzz_0, ta1_z_yzzz_xzz_1, ta1_z_yzzz_yyz_0, ta1_z_yzzz_yyz_1, ta1_z_yzzz_yz_0, ta1_z_yzzz_yz_1, ta1_z_yzzz_yzz_0, ta1_z_yzzz_yzz_1, ta1_z_yzzz_zz_0, ta1_z_yzzz_zz_1, ta1_z_yzzz_zzz_0, ta1_z_yzzz_zzz_1, ta1_z_zzz_xxx_0, ta1_z_zzz_xxx_1, ta1_z_zzz_xxz_0, ta1_z_zzz_xxz_1, ta1_z_zzz_xyz_0, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_0, ta1_z_zzz_xzz_1, ta1_z_zzz_yyz_0, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_0, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_0, ta1_z_zzz_zzz_1, ta_yyzz_xxy_1, ta_yyzz_xyy_1, ta_yyzz_yyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzz_xxx_0[i] = ta1_z_zzz_xxx_0[i] * fe_0 - ta1_z_zzz_xxx_1[i] * fe_0 + ta1_z_yzzz_xxx_0[i] * pa_y[i] - ta1_z_yzzz_xxx_1[i] * pc_y[i];

        ta1_z_yyzzz_xxy_0[i] = 2.0 * ta1_z_yyz_xxy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xxy_1[i] * fe_0 + ta_yyzz_xxy_1[i] + ta1_z_yyzz_xxy_0[i] * pa_z[i] - ta1_z_yyzz_xxy_1[i] * pc_z[i];

        ta1_z_yyzzz_xxz_0[i] = ta1_z_zzz_xxz_0[i] * fe_0 - ta1_z_zzz_xxz_1[i] * fe_0 + ta1_z_yzzz_xxz_0[i] * pa_y[i] - ta1_z_yzzz_xxz_1[i] * pc_y[i];

        ta1_z_yyzzz_xyy_0[i] = 2.0 * ta1_z_yyz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyy_1[i] * fe_0 + ta_yyzz_xyy_1[i] + ta1_z_yyzz_xyy_0[i] * pa_z[i] - ta1_z_yyzz_xyy_1[i] * pc_z[i];

        ta1_z_yyzzz_xyz_0[i] = ta1_z_zzz_xyz_0[i] * fe_0 - ta1_z_zzz_xyz_1[i] * fe_0 + ta1_z_yzzz_xz_0[i] * fe_0 - ta1_z_yzzz_xz_1[i] * fe_0 + ta1_z_yzzz_xyz_0[i] * pa_y[i] - ta1_z_yzzz_xyz_1[i] * pc_y[i];

        ta1_z_yyzzz_xzz_0[i] = ta1_z_zzz_xzz_0[i] * fe_0 - ta1_z_zzz_xzz_1[i] * fe_0 + ta1_z_yzzz_xzz_0[i] * pa_y[i] - ta1_z_yzzz_xzz_1[i] * pc_y[i];

        ta1_z_yyzzz_yyy_0[i] = 2.0 * ta1_z_yyz_yyy_0[i] * fe_0 - 2.0 * ta1_z_yyz_yyy_1[i] * fe_0 + ta_yyzz_yyy_1[i] + ta1_z_yyzz_yyy_0[i] * pa_z[i] - ta1_z_yyzz_yyy_1[i] * pc_z[i];

        ta1_z_yyzzz_yyz_0[i] = ta1_z_zzz_yyz_0[i] * fe_0 - ta1_z_zzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yz_1[i] * fe_0 + ta1_z_yzzz_yyz_0[i] * pa_y[i] - ta1_z_yzzz_yyz_1[i] * pc_y[i];

        ta1_z_yyzzz_yzz_0[i] = ta1_z_zzz_yzz_0[i] * fe_0 - ta1_z_zzz_yzz_1[i] * fe_0 + ta1_z_yzzz_zz_0[i] * fe_0 - ta1_z_yzzz_zz_1[i] * fe_0 + ta1_z_yzzz_yzz_0[i] * pa_y[i] - ta1_z_yzzz_yzz_1[i] * pc_y[i];

        ta1_z_yyzzz_zzz_0[i] = ta1_z_zzz_zzz_0[i] * fe_0 - ta1_z_zzz_zzz_1[i] * fe_0 + ta1_z_yzzz_zzz_0[i] * pa_y[i] - ta1_z_yzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 610-620 components of targeted buffer : HF

    auto ta1_z_yzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 610);

    auto ta1_z_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 611);

    auto ta1_z_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 612);

    auto ta1_z_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 613);

    auto ta1_z_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 614);

    auto ta1_z_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 615);

    auto ta1_z_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 616);

    auto ta1_z_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 617);

    auto ta1_z_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 618);

    auto ta1_z_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 619);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yzzzz_xxx_0, ta1_z_yzzzz_xxy_0, ta1_z_yzzzz_xxz_0, ta1_z_yzzzz_xyy_0, ta1_z_yzzzz_xyz_0, ta1_z_yzzzz_xzz_0, ta1_z_yzzzz_yyy_0, ta1_z_yzzzz_yyz_0, ta1_z_yzzzz_yzz_0, ta1_z_yzzzz_zzz_0, ta1_z_zzzz_xx_0, ta1_z_zzzz_xx_1, ta1_z_zzzz_xxx_0, ta1_z_zzzz_xxx_1, ta1_z_zzzz_xxy_0, ta1_z_zzzz_xxy_1, ta1_z_zzzz_xxz_0, ta1_z_zzzz_xxz_1, ta1_z_zzzz_xy_0, ta1_z_zzzz_xy_1, ta1_z_zzzz_xyy_0, ta1_z_zzzz_xyy_1, ta1_z_zzzz_xyz_0, ta1_z_zzzz_xyz_1, ta1_z_zzzz_xz_0, ta1_z_zzzz_xz_1, ta1_z_zzzz_xzz_0, ta1_z_zzzz_xzz_1, ta1_z_zzzz_yy_0, ta1_z_zzzz_yy_1, ta1_z_zzzz_yyy_0, ta1_z_zzzz_yyy_1, ta1_z_zzzz_yyz_0, ta1_z_zzzz_yyz_1, ta1_z_zzzz_yz_0, ta1_z_zzzz_yz_1, ta1_z_zzzz_yzz_0, ta1_z_zzzz_yzz_1, ta1_z_zzzz_zz_0, ta1_z_zzzz_zz_1, ta1_z_zzzz_zzz_0, ta1_z_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzz_xxx_0[i] = ta1_z_zzzz_xxx_0[i] * pa_y[i] - ta1_z_zzzz_xxx_1[i] * pc_y[i];

        ta1_z_yzzzz_xxy_0[i] = ta1_z_zzzz_xx_0[i] * fe_0 - ta1_z_zzzz_xx_1[i] * fe_0 + ta1_z_zzzz_xxy_0[i] * pa_y[i] - ta1_z_zzzz_xxy_1[i] * pc_y[i];

        ta1_z_yzzzz_xxz_0[i] = ta1_z_zzzz_xxz_0[i] * pa_y[i] - ta1_z_zzzz_xxz_1[i] * pc_y[i];

        ta1_z_yzzzz_xyy_0[i] = 2.0 * ta1_z_zzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xy_1[i] * fe_0 + ta1_z_zzzz_xyy_0[i] * pa_y[i] - ta1_z_zzzz_xyy_1[i] * pc_y[i];

        ta1_z_yzzzz_xyz_0[i] = ta1_z_zzzz_xz_0[i] * fe_0 - ta1_z_zzzz_xz_1[i] * fe_0 + ta1_z_zzzz_xyz_0[i] * pa_y[i] - ta1_z_zzzz_xyz_1[i] * pc_y[i];

        ta1_z_yzzzz_xzz_0[i] = ta1_z_zzzz_xzz_0[i] * pa_y[i] - ta1_z_zzzz_xzz_1[i] * pc_y[i];

        ta1_z_yzzzz_yyy_0[i] = 3.0 * ta1_z_zzzz_yy_0[i] * fe_0 - 3.0 * ta1_z_zzzz_yy_1[i] * fe_0 + ta1_z_zzzz_yyy_0[i] * pa_y[i] - ta1_z_zzzz_yyy_1[i] * pc_y[i];

        ta1_z_yzzzz_yyz_0[i] = 2.0 * ta1_z_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_yz_1[i] * fe_0 + ta1_z_zzzz_yyz_0[i] * pa_y[i] - ta1_z_zzzz_yyz_1[i] * pc_y[i];

        ta1_z_yzzzz_yzz_0[i] = ta1_z_zzzz_zz_0[i] * fe_0 - ta1_z_zzzz_zz_1[i] * fe_0 + ta1_z_zzzz_yzz_0[i] * pa_y[i] - ta1_z_zzzz_yzz_1[i] * pc_y[i];

        ta1_z_yzzzz_zzz_0[i] = ta1_z_zzzz_zzz_0[i] * pa_y[i] - ta1_z_zzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 620-630 components of targeted buffer : HF

    auto ta1_z_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 620);

    auto ta1_z_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 621);

    auto ta1_z_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 622);

    auto ta1_z_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 623);

    auto ta1_z_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 624);

    auto ta1_z_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 625);

    auto ta1_z_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 626);

    auto ta1_z_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 627);

    auto ta1_z_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 628);

    auto ta1_z_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 629);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zzz_xxx_0, ta1_z_zzz_xxx_1, ta1_z_zzz_xxy_0, ta1_z_zzz_xxy_1, ta1_z_zzz_xxz_0, ta1_z_zzz_xxz_1, ta1_z_zzz_xyy_0, ta1_z_zzz_xyy_1, ta1_z_zzz_xyz_0, ta1_z_zzz_xyz_1, ta1_z_zzz_xzz_0, ta1_z_zzz_xzz_1, ta1_z_zzz_yyy_0, ta1_z_zzz_yyy_1, ta1_z_zzz_yyz_0, ta1_z_zzz_yyz_1, ta1_z_zzz_yzz_0, ta1_z_zzz_yzz_1, ta1_z_zzz_zzz_0, ta1_z_zzz_zzz_1, ta1_z_zzzz_xx_0, ta1_z_zzzz_xx_1, ta1_z_zzzz_xxx_0, ta1_z_zzzz_xxx_1, ta1_z_zzzz_xxy_0, ta1_z_zzzz_xxy_1, ta1_z_zzzz_xxz_0, ta1_z_zzzz_xxz_1, ta1_z_zzzz_xy_0, ta1_z_zzzz_xy_1, ta1_z_zzzz_xyy_0, ta1_z_zzzz_xyy_1, ta1_z_zzzz_xyz_0, ta1_z_zzzz_xyz_1, ta1_z_zzzz_xz_0, ta1_z_zzzz_xz_1, ta1_z_zzzz_xzz_0, ta1_z_zzzz_xzz_1, ta1_z_zzzz_yy_0, ta1_z_zzzz_yy_1, ta1_z_zzzz_yyy_0, ta1_z_zzzz_yyy_1, ta1_z_zzzz_yyz_0, ta1_z_zzzz_yyz_1, ta1_z_zzzz_yz_0, ta1_z_zzzz_yz_1, ta1_z_zzzz_yzz_0, ta1_z_zzzz_yzz_1, ta1_z_zzzz_zz_0, ta1_z_zzzz_zz_1, ta1_z_zzzz_zzz_0, ta1_z_zzzz_zzz_1, ta1_z_zzzzz_xxx_0, ta1_z_zzzzz_xxy_0, ta1_z_zzzzz_xxz_0, ta1_z_zzzzz_xyy_0, ta1_z_zzzzz_xyz_0, ta1_z_zzzzz_xzz_0, ta1_z_zzzzz_yyy_0, ta1_z_zzzzz_yyz_0, ta1_z_zzzzz_yzz_0, ta1_z_zzzzz_zzz_0, ta_zzzz_xxx_1, ta_zzzz_xxy_1, ta_zzzz_xxz_1, ta_zzzz_xyy_1, ta_zzzz_xyz_1, ta_zzzz_xzz_1, ta_zzzz_yyy_1, ta_zzzz_yyz_1, ta_zzzz_yzz_1, ta_zzzz_zzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzz_xxx_0[i] = 4.0 * ta1_z_zzz_xxx_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxx_1[i] * fe_0 + ta_zzzz_xxx_1[i] + ta1_z_zzzz_xxx_0[i] * pa_z[i] - ta1_z_zzzz_xxx_1[i] * pc_z[i];

        ta1_z_zzzzz_xxy_0[i] = 4.0 * ta1_z_zzz_xxy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxy_1[i] * fe_0 + ta_zzzz_xxy_1[i] + ta1_z_zzzz_xxy_0[i] * pa_z[i] - ta1_z_zzzz_xxy_1[i] * pc_z[i];

        ta1_z_zzzzz_xxz_0[i] = 4.0 * ta1_z_zzz_xxz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxz_1[i] * fe_0 + ta1_z_zzzz_xx_0[i] * fe_0 - ta1_z_zzzz_xx_1[i] * fe_0 + ta_zzzz_xxz_1[i] + ta1_z_zzzz_xxz_0[i] * pa_z[i] - ta1_z_zzzz_xxz_1[i] * pc_z[i];

        ta1_z_zzzzz_xyy_0[i] = 4.0 * ta1_z_zzz_xyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyy_1[i] * fe_0 + ta_zzzz_xyy_1[i] + ta1_z_zzzz_xyy_0[i] * pa_z[i] - ta1_z_zzzz_xyy_1[i] * pc_z[i];

        ta1_z_zzzzz_xyz_0[i] = 4.0 * ta1_z_zzz_xyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyz_1[i] * fe_0 + ta1_z_zzzz_xy_0[i] * fe_0 - ta1_z_zzzz_xy_1[i] * fe_0 + ta_zzzz_xyz_1[i] + ta1_z_zzzz_xyz_0[i] * pa_z[i] - ta1_z_zzzz_xyz_1[i] * pc_z[i];

        ta1_z_zzzzz_xzz_0[i] = 4.0 * ta1_z_zzz_xzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_xz_1[i] * fe_0 + ta_zzzz_xzz_1[i] + ta1_z_zzzz_xzz_0[i] * pa_z[i] - ta1_z_zzzz_xzz_1[i] * pc_z[i];

        ta1_z_zzzzz_yyy_0[i] = 4.0 * ta1_z_zzz_yyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyy_1[i] * fe_0 + ta_zzzz_yyy_1[i] + ta1_z_zzzz_yyy_0[i] * pa_z[i] - ta1_z_zzzz_yyy_1[i] * pc_z[i];

        ta1_z_zzzzz_yyz_0[i] = 4.0 * ta1_z_zzz_yyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyz_1[i] * fe_0 + ta1_z_zzzz_yy_0[i] * fe_0 - ta1_z_zzzz_yy_1[i] * fe_0 + ta_zzzz_yyz_1[i] + ta1_z_zzzz_yyz_0[i] * pa_z[i] - ta1_z_zzzz_yyz_1[i] * pc_z[i];

        ta1_z_zzzzz_yzz_0[i] = 4.0 * ta1_z_zzz_yzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_zzzz_yz_1[i] * fe_0 + ta_zzzz_yzz_1[i] + ta1_z_zzzz_yzz_0[i] * pa_z[i] - ta1_z_zzzz_yzz_1[i] * pc_z[i];

        ta1_z_zzzzz_zzz_0[i] = 4.0 * ta1_z_zzz_zzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_zzz_1[i] * fe_0 + 3.0 * ta1_z_zzzz_zz_0[i] * fe_0 - 3.0 * ta1_z_zzzz_zz_1[i] * fe_0 + ta_zzzz_zzz_1[i] + ta1_z_zzzz_zzz_0[i] * pa_z[i] - ta1_z_zzzz_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

