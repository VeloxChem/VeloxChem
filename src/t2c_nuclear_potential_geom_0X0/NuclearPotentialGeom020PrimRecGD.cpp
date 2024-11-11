#include "NuclearPotentialGeom020PrimRecGD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_gd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_gd,
                                        const size_t              idx_npot_geom_020_0_dd,
                                        const size_t              idx_npot_geom_020_1_dd,
                                        const size_t              idx_npot_geom_020_0_fp,
                                        const size_t              idx_npot_geom_020_1_fp,
                                        const size_t              idx_npot_geom_010_1_fd,
                                        const size_t              idx_npot_geom_020_0_fd,
                                        const size_t              idx_npot_geom_020_1_fd,
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

    auto ta2_xx_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd);

    auto ta2_xx_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 1);

    auto ta2_xx_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 2);

    auto ta2_xx_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 3);

    auto ta2_xx_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 4);

    auto ta2_xx_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 5);

    auto ta2_xx_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 6);

    auto ta2_xx_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 8);

    auto ta2_xx_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 12);

    auto ta2_xx_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 13);

    auto ta2_xx_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 14);

    auto ta2_xx_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 18);

    auto ta2_xx_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 19);

    auto ta2_xx_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 20);

    auto ta2_xx_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 21);

    auto ta2_xx_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 22);

    auto ta2_xx_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 23);

    auto ta2_xx_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 26);

    auto ta2_xx_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 29);

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

    auto ta2_xy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 45);

    auto ta2_xy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 46);

    auto ta2_xy_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 48);

    auto ta2_xy_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 49);

    auto ta2_xy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 54);

    auto ta2_xy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 55);

    auto ta2_xy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 56);

    auto ta2_xy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 57);

    auto ta2_xy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 58);

    auto ta2_xy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 59);

    auto ta2_xy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 61);

    auto ta2_xy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 63);

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

    auto ta2_xz_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 78);

    auto ta2_xz_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 80);

    auto ta2_xz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 86);

    auto ta2_xz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 88);

    auto ta2_xz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 89);

    auto ta2_xz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 90);

    auto ta2_xz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 91);

    auto ta2_xz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 92);

    auto ta2_xz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 93);

    auto ta2_xz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 94);

    auto ta2_xz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 95);

    auto ta2_xz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 98);

    auto ta2_xz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 101);

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

    auto ta2_yy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 117);

    auto ta2_yy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 118);

    auto ta2_yy_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 124);

    auto ta2_yy_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 125);

    auto ta2_yy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 126);

    auto ta2_yy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 127);

    auto ta2_yy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 128);

    auto ta2_yy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 129);

    auto ta2_yy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 130);

    auto ta2_yy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 131);

    auto ta2_yy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 133);

    auto ta2_yy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 135);

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

    auto ta2_yz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 153);

    auto ta2_yz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 154);

    auto ta2_yz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 160);

    auto ta2_yz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 161);

    auto ta2_yz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 162);

    auto ta2_yz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 163);

    auto ta2_yz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 164);

    auto ta2_yz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 165);

    auto ta2_yz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 166);

    auto ta2_yz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 167);

    auto ta2_yz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 170);

    auto ta2_yz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 172);

    auto ta2_yz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 173);

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

    auto ta2_zz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 189);

    auto ta2_zz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 190);

    auto ta2_zz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 196);

    auto ta2_zz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 197);

    auto ta2_zz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 198);

    auto ta2_zz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 199);

    auto ta2_zz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 200);

    auto ta2_zz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 201);

    auto ta2_zz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 202);

    auto ta2_zz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 203);

    auto ta2_zz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 206);

    auto ta2_zz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 208);

    auto ta2_zz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 209);

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

    auto ta2_xx_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 6);

    auto ta2_xx_xy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 8);

    auto ta2_xx_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 12);

    auto ta2_xx_xz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 13);

    auto ta2_xx_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 14);

    auto ta2_xx_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 18);

    auto ta2_xx_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 19);

    auto ta2_xx_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 20);

    auto ta2_xx_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 21);

    auto ta2_xx_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 22);

    auto ta2_xx_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 23);

    auto ta2_xx_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 26);

    auto ta2_xx_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 29);

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

    auto ta2_xy_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 45);

    auto ta2_xy_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 46);

    auto ta2_xy_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 48);

    auto ta2_xy_xz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 49);

    auto ta2_xy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 54);

    auto ta2_xy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 55);

    auto ta2_xy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 56);

    auto ta2_xy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 57);

    auto ta2_xy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 58);

    auto ta2_xy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 59);

    auto ta2_xy_yz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 61);

    auto ta2_xy_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 63);

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

    auto ta2_xz_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 78);

    auto ta2_xz_xy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 80);

    auto ta2_xz_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 86);

    auto ta2_xz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 88);

    auto ta2_xz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 89);

    auto ta2_xz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 90);

    auto ta2_xz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 91);

    auto ta2_xz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 92);

    auto ta2_xz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 93);

    auto ta2_xz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 94);

    auto ta2_xz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 95);

    auto ta2_xz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 98);

    auto ta2_xz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 101);

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

    auto ta2_yy_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 117);

    auto ta2_yy_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 118);

    auto ta2_yy_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 124);

    auto ta2_yy_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 125);

    auto ta2_yy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 126);

    auto ta2_yy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 127);

    auto ta2_yy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 128);

    auto ta2_yy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 129);

    auto ta2_yy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 130);

    auto ta2_yy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 131);

    auto ta2_yy_yz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 133);

    auto ta2_yy_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 135);

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

    auto ta2_yz_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 153);

    auto ta2_yz_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 154);

    auto ta2_yz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 160);

    auto ta2_yz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 161);

    auto ta2_yz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 162);

    auto ta2_yz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 163);

    auto ta2_yz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 164);

    auto ta2_yz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 165);

    auto ta2_yz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 166);

    auto ta2_yz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 167);

    auto ta2_yz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 170);

    auto ta2_yz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 172);

    auto ta2_yz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 173);

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

    auto ta2_zz_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 189);

    auto ta2_zz_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 190);

    auto ta2_zz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 196);

    auto ta2_zz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 197);

    auto ta2_zz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 198);

    auto ta2_zz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 199);

    auto ta2_zz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 200);

    auto ta2_zz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 201);

    auto ta2_zz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 202);

    auto ta2_zz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 203);

    auto ta2_zz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 206);

    auto ta2_zz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 208);

    auto ta2_zz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 209);

    auto ta2_zz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 210);

    auto ta2_zz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 211);

    auto ta2_zz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 212);

    auto ta2_zz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 213);

    auto ta2_zz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 214);

    auto ta2_zz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 215);

    // Set up components of auxiliary buffer : FP

    auto ta2_xx_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp);

    auto ta2_xx_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 1);

    auto ta2_xx_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 2);

    auto ta2_xx_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 3);

    auto ta2_xx_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 6);

    auto ta2_xx_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 8);

    auto ta2_xx_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 15);

    auto ta2_xx_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 18);

    auto ta2_xx_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 19);

    auto ta2_xx_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 20);

    auto ta2_xx_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 26);

    auto ta2_xx_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 27);

    auto ta2_xx_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 28);

    auto ta2_xx_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 29);

    auto ta2_xy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 30);

    auto ta2_xy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 31);

    auto ta2_xy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 32);

    auto ta2_xy_xxy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 34);

    auto ta2_xy_xxz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 36);

    auto ta2_xy_xyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 39);

    auto ta2_xy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 40);

    auto ta2_xy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 48);

    auto ta2_xy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 49);

    auto ta2_xy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 50);

    auto ta2_xy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 52);

    auto ta2_xy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 57);

    auto ta2_xy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 58);

    auto ta2_xy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 59);

    auto ta2_xz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 60);

    auto ta2_xz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 61);

    auto ta2_xz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 62);

    auto ta2_xz_xxy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 63);

    auto ta2_xz_xxz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 68);

    auto ta2_xz_xzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 75);

    auto ta2_xz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 77);

    auto ta2_xz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 78);

    auto ta2_xz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 79);

    auto ta2_xz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 80);

    auto ta2_xz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 86);

    auto ta2_xz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 87);

    auto ta2_xz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 88);

    auto ta2_xz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 89);

    auto ta2_yy_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 90);

    auto ta2_yy_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 91);

    auto ta2_yy_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 92);

    auto ta2_yy_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 100);

    auto ta2_yy_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 107);

    auto ta2_yy_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 108);

    auto ta2_yy_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 109);

    auto ta2_yy_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 110);

    auto ta2_yy_yyz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 112);

    auto ta2_yy_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 113);

    auto ta2_yy_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 115);

    auto ta2_yy_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 117);

    auto ta2_yy_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 118);

    auto ta2_yy_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 119);

    auto ta2_yz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 120);

    auto ta2_yz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 121);

    auto ta2_yz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 122);

    auto ta2_yz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 130);

    auto ta2_yz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 137);

    auto ta2_yz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 138);

    auto ta2_yz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 139);

    auto ta2_yz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 140);

    auto ta2_yz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 143);

    auto ta2_yz_yzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 145);

    auto ta2_yz_yzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 146);

    auto ta2_yz_zzz_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 147);

    auto ta2_yz_zzz_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 148);

    auto ta2_yz_zzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 149);

    auto ta2_zz_xxx_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 150);

    auto ta2_zz_xxx_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 151);

    auto ta2_zz_xxx_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 152);

    auto ta2_zz_xyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 160);

    auto ta2_zz_xzz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 167);

    auto ta2_zz_yyy_x_0 = pbuffer.data(idx_npot_geom_020_0_fp + 168);

    auto ta2_zz_yyy_y_0 = pbuffer.data(idx_npot_geom_020_0_fp + 169);

    auto ta2_zz_yyy_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 170);

    auto ta2_zz_yyz_z_0 = pbuffer.data(idx_npot_geom_020_0_fp + 173);

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

    auto ta2_xx_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 6);

    auto ta2_xx_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 8);

    auto ta2_xx_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 15);

    auto ta2_xx_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 18);

    auto ta2_xx_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 19);

    auto ta2_xx_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 20);

    auto ta2_xx_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 26);

    auto ta2_xx_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 27);

    auto ta2_xx_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 28);

    auto ta2_xx_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 29);

    auto ta2_xy_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 30);

    auto ta2_xy_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 31);

    auto ta2_xy_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 32);

    auto ta2_xy_xxy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 34);

    auto ta2_xy_xxz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 36);

    auto ta2_xy_xyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 39);

    auto ta2_xy_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 40);

    auto ta2_xy_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 48);

    auto ta2_xy_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 49);

    auto ta2_xy_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 50);

    auto ta2_xy_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 52);

    auto ta2_xy_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 57);

    auto ta2_xy_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 58);

    auto ta2_xy_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 59);

    auto ta2_xz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 60);

    auto ta2_xz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 61);

    auto ta2_xz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 62);

    auto ta2_xz_xxy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 63);

    auto ta2_xz_xxz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 68);

    auto ta2_xz_xzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 75);

    auto ta2_xz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 77);

    auto ta2_xz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 78);

    auto ta2_xz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 79);

    auto ta2_xz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 80);

    auto ta2_xz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 86);

    auto ta2_xz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 87);

    auto ta2_xz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 88);

    auto ta2_xz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 89);

    auto ta2_yy_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 90);

    auto ta2_yy_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 91);

    auto ta2_yy_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 92);

    auto ta2_yy_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 100);

    auto ta2_yy_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 107);

    auto ta2_yy_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 108);

    auto ta2_yy_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 109);

    auto ta2_yy_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 110);

    auto ta2_yy_yyz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 112);

    auto ta2_yy_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 113);

    auto ta2_yy_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 115);

    auto ta2_yy_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 117);

    auto ta2_yy_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 118);

    auto ta2_yy_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 119);

    auto ta2_yz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 120);

    auto ta2_yz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 121);

    auto ta2_yz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 122);

    auto ta2_yz_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 130);

    auto ta2_yz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 137);

    auto ta2_yz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 138);

    auto ta2_yz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 139);

    auto ta2_yz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 140);

    auto ta2_yz_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 143);

    auto ta2_yz_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 145);

    auto ta2_yz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 146);

    auto ta2_yz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 147);

    auto ta2_yz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 148);

    auto ta2_yz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 149);

    auto ta2_zz_xxx_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 150);

    auto ta2_zz_xxx_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 151);

    auto ta2_zz_xxx_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 152);

    auto ta2_zz_xyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 160);

    auto ta2_zz_xzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 167);

    auto ta2_zz_yyy_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 168);

    auto ta2_zz_yyy_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 169);

    auto ta2_zz_yyy_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 170);

    auto ta2_zz_yyz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 173);

    auto ta2_zz_yzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 175);

    auto ta2_zz_yzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 176);

    auto ta2_zz_zzz_x_1 = pbuffer.data(idx_npot_geom_020_1_fp + 177);

    auto ta2_zz_zzz_y_1 = pbuffer.data(idx_npot_geom_020_1_fp + 178);

    auto ta2_zz_zzz_z_1 = pbuffer.data(idx_npot_geom_020_1_fp + 179);

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

    auto ta1_x_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 12);

    auto ta1_x_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 13);

    auto ta1_x_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 14);

    auto ta1_x_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 17);

    auto ta1_x_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 18);

    auto ta1_x_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 19);

    auto ta1_x_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 21);

    auto ta1_x_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 22);

    auto ta1_x_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 30);

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

    auto ta1_x_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 45);

    auto ta1_x_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 46);

    auto ta1_x_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 47);

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

    auto ta1_y_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 77);

    auto ta1_y_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 78);

    auto ta1_y_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 79);

    auto ta1_y_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 81);

    auto ta1_y_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 82);

    auto ta1_y_xyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 83);

    auto ta1_y_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 90);

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

    auto ta1_y_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 105);

    auto ta1_y_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 106);

    auto ta1_y_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 107);

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

    auto ta1_z_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 132);

    auto ta1_z_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 133);

    auto ta1_z_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 134);

    auto ta1_z_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 136);

    auto ta1_z_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 137);

    auto ta1_z_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 138);

    auto ta1_z_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 139);

    auto ta1_z_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 141);

    auto ta1_z_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 142);

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

    auto ta2_xx_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 9);

    auto ta2_xx_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 11);

    auto ta2_xx_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 12);

    auto ta2_xx_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 13);

    auto ta2_xx_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 14);

    auto ta2_xx_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 15);

    auto ta2_xx_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 16);

    auto ta2_xx_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 17);

    auto ta2_xx_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 18);

    auto ta2_xx_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 19);

    auto ta2_xx_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 20);

    auto ta2_xx_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 21);

    auto ta2_xx_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 22);

    auto ta2_xx_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 26);

    auto ta2_xx_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 30);

    auto ta2_xx_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 31);

    auto ta2_xx_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 32);

    auto ta2_xx_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 34);

    auto ta2_xx_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 35);

    auto ta2_xx_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 36);

    auto ta2_xx_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 37);

    auto ta2_xx_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 38);

    auto ta2_xx_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 39);

    auto ta2_xx_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 40);

    auto ta2_xx_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 41);

    auto ta2_xx_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 43);

    auto ta2_xx_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 44);

    auto ta2_xx_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 45);

    auto ta2_xx_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 46);

    auto ta2_xx_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 47);

    auto ta2_xx_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 48);

    auto ta2_xx_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 50);

    auto ta2_xx_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 51);

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

    auto ta2_xy_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 66);

    auto ta2_xy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 67);

    auto ta2_xy_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 68);

    auto ta2_xy_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 69);

    auto ta2_xy_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 70);

    auto ta2_xy_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 72);

    auto ta2_xy_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 73);

    auto ta2_xy_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 74);

    auto ta2_xy_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 75);

    auto ta2_xy_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 77);

    auto ta2_xy_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 78);

    auto ta2_xy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 79);

    auto ta2_xy_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 80);

    auto ta2_xy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 81);

    auto ta2_xy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 82);

    auto ta2_xy_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 83);

    auto ta2_xy_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 85);

    auto ta2_xy_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 90);

    auto ta2_xy_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 91);

    auto ta2_xy_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 92);

    auto ta2_xy_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 94);

    auto ta2_xy_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 95);

    auto ta2_xy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 96);

    auto ta2_xy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 97);

    auto ta2_xy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 98);

    auto ta2_xy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 99);

    auto ta2_xy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 100);

    auto ta2_xy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 101);

    auto ta2_xy_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 102);

    auto ta2_xy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 103);

    auto ta2_xy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 105);

    auto ta2_xy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 106);

    auto ta2_xy_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 107);

    auto ta2_xy_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 109);

    auto ta2_xy_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 110);

    auto ta2_xy_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 111);

    auto ta2_xy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 112);

    auto ta2_xy_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 113);

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

    auto ta2_xz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 129);

    auto ta2_xz_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 131);

    auto ta2_xz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 132);

    auto ta2_xz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 133);

    auto ta2_xz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 134);

    auto ta2_xz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 136);

    auto ta2_xz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 137);

    auto ta2_xz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 138);

    auto ta2_xz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 139);

    auto ta2_xz_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 140);

    auto ta2_xz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 141);

    auto ta2_xz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 142);

    auto ta2_xz_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 146);

    auto ta2_xz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 150);

    auto ta2_xz_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 151);

    auto ta2_xz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 152);

    auto ta2_xz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 153);

    auto ta2_xz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 154);

    auto ta2_xz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 155);

    auto ta2_xz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 156);

    auto ta2_xz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 157);

    auto ta2_xz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 158);

    auto ta2_xz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 159);

    auto ta2_xz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 160);

    auto ta2_xz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 161);

    auto ta2_xz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 163);

    auto ta2_xz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 164);

    auto ta2_xz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 165);

    auto ta2_xz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 166);

    auto ta2_xz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 167);

    auto ta2_xz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 168);

    auto ta2_xz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 170);

    auto ta2_xz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 171);

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

    auto ta2_yy_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 186);

    auto ta2_yy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 187);

    auto ta2_yy_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 188);

    auto ta2_yy_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 189);

    auto ta2_yy_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 190);

    auto ta2_yy_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 192);

    auto ta2_yy_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 193);

    auto ta2_yy_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 194);

    auto ta2_yy_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 196);

    auto ta2_yy_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 197);

    auto ta2_yy_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 198);

    auto ta2_yy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 199);

    auto ta2_yy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 201);

    auto ta2_yy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 202);

    auto ta2_yy_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 203);

    auto ta2_yy_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 208);

    auto ta2_yy_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 210);

    auto ta2_yy_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 212);

    auto ta2_yy_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 213);

    auto ta2_yy_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 214);

    auto ta2_yy_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 215);

    auto ta2_yy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 216);

    auto ta2_yy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 217);

    auto ta2_yy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 218);

    auto ta2_yy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 219);

    auto ta2_yy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 220);

    auto ta2_yy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 221);

    auto ta2_yy_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 222);

    auto ta2_yy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 223);

    auto ta2_yy_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 224);

    auto ta2_yy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 225);

    auto ta2_yy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 226);

    auto ta2_yy_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 227);

    auto ta2_yy_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 229);

    auto ta2_yy_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 230);

    auto ta2_yy_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 231);

    auto ta2_yy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 232);

    auto ta2_yy_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 233);

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

    auto ta2_yz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 246);

    auto ta2_yz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 247);

    auto ta2_yz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 248);

    auto ta2_yz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 249);

    auto ta2_yz_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 250);

    auto ta2_yz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 252);

    auto ta2_yz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 253);

    auto ta2_yz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 254);

    auto ta2_yz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 256);

    auto ta2_yz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 257);

    auto ta2_yz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 258);

    auto ta2_yz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 259);

    auto ta2_yz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 261);

    auto ta2_yz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 262);

    auto ta2_yz_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 263);

    auto ta2_yz_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 268);

    auto ta2_yz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 270);

    auto ta2_yz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 272);

    auto ta2_yz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 273);

    auto ta2_yz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 274);

    auto ta2_yz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 275);

    auto ta2_yz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 276);

    auto ta2_yz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 277);

    auto ta2_yz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 278);

    auto ta2_yz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 279);

    auto ta2_yz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 280);

    auto ta2_yz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 281);

    auto ta2_yz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 283);

    auto ta2_yz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 284);

    auto ta2_yz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 285);

    auto ta2_yz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 286);

    auto ta2_yz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 287);

    auto ta2_yz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 288);

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

    auto ta2_zz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 306);

    auto ta2_zz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 307);

    auto ta2_zz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 308);

    auto ta2_zz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 309);

    auto ta2_zz_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 310);

    auto ta2_zz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 312);

    auto ta2_zz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 313);

    auto ta2_zz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 314);

    auto ta2_zz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 316);

    auto ta2_zz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 317);

    auto ta2_zz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 318);

    auto ta2_zz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 319);

    auto ta2_zz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 321);

    auto ta2_zz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 322);

    auto ta2_zz_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 323);

    auto ta2_zz_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 328);

    auto ta2_zz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 330);

    auto ta2_zz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 332);

    auto ta2_zz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 333);

    auto ta2_zz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 334);

    auto ta2_zz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 335);

    auto ta2_zz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 336);

    auto ta2_zz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 337);

    auto ta2_zz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 338);

    auto ta2_zz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 339);

    auto ta2_zz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 340);

    auto ta2_zz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 341);

    auto ta2_zz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 343);

    auto ta2_zz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 344);

    auto ta2_zz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 345);

    auto ta2_zz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 346);

    auto ta2_zz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 347);

    auto ta2_zz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 348);

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

    auto ta2_xx_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 9);

    auto ta2_xx_xxy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 11);

    auto ta2_xx_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 12);

    auto ta2_xx_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 13);

    auto ta2_xx_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 14);

    auto ta2_xx_xxz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 15);

    auto ta2_xx_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 16);

    auto ta2_xx_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 17);

    auto ta2_xx_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 18);

    auto ta2_xx_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 19);

    auto ta2_xx_xyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 20);

    auto ta2_xx_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 21);

    auto ta2_xx_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 22);

    auto ta2_xx_xyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 26);

    auto ta2_xx_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 30);

    auto ta2_xx_xzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 31);

    auto ta2_xx_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 32);

    auto ta2_xx_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 34);

    auto ta2_xx_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 35);

    auto ta2_xx_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 36);

    auto ta2_xx_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 37);

    auto ta2_xx_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 38);

    auto ta2_xx_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 39);

    auto ta2_xx_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 40);

    auto ta2_xx_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 41);

    auto ta2_xx_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 43);

    auto ta2_xx_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 44);

    auto ta2_xx_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 45);

    auto ta2_xx_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 46);

    auto ta2_xx_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 47);

    auto ta2_xx_yzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 48);

    auto ta2_xx_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 50);

    auto ta2_xx_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 51);

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

    auto ta2_xy_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 66);

    auto ta2_xy_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 67);

    auto ta2_xy_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 68);

    auto ta2_xy_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 69);

    auto ta2_xy_xxy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 70);

    auto ta2_xy_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 72);

    auto ta2_xy_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 73);

    auto ta2_xy_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 74);

    auto ta2_xy_xxz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 75);

    auto ta2_xy_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 77);

    auto ta2_xy_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 78);

    auto ta2_xy_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 79);

    auto ta2_xy_xyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 80);

    auto ta2_xy_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 81);

    auto ta2_xy_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 82);

    auto ta2_xy_xyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 83);

    auto ta2_xy_xyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 85);

    auto ta2_xy_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 90);

    auto ta2_xy_xzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 91);

    auto ta2_xy_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 92);

    auto ta2_xy_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 94);

    auto ta2_xy_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 95);

    auto ta2_xy_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 96);

    auto ta2_xy_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 97);

    auto ta2_xy_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 98);

    auto ta2_xy_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 99);

    auto ta2_xy_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 100);

    auto ta2_xy_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 101);

    auto ta2_xy_yyz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 102);

    auto ta2_xy_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 103);

    auto ta2_xy_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 105);

    auto ta2_xy_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 106);

    auto ta2_xy_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 107);

    auto ta2_xy_yzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 109);

    auto ta2_xy_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 110);

    auto ta2_xy_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 111);

    auto ta2_xy_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 112);

    auto ta2_xy_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 113);

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

    auto ta2_xz_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 129);

    auto ta2_xz_xxy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 131);

    auto ta2_xz_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 132);

    auto ta2_xz_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 133);

    auto ta2_xz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 134);

    auto ta2_xz_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 136);

    auto ta2_xz_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 137);

    auto ta2_xz_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 138);

    auto ta2_xz_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 139);

    auto ta2_xz_xyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 140);

    auto ta2_xz_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 141);

    auto ta2_xz_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 142);

    auto ta2_xz_xyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 146);

    auto ta2_xz_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 150);

    auto ta2_xz_xzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 151);

    auto ta2_xz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 152);

    auto ta2_xz_xzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 153);

    auto ta2_xz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 154);

    auto ta2_xz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 155);

    auto ta2_xz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 156);

    auto ta2_xz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 157);

    auto ta2_xz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 158);

    auto ta2_xz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 159);

    auto ta2_xz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 160);

    auto ta2_xz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 161);

    auto ta2_xz_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 163);

    auto ta2_xz_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 164);

    auto ta2_xz_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 165);

    auto ta2_xz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 166);

    auto ta2_xz_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 167);

    auto ta2_xz_yzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 168);

    auto ta2_xz_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 170);

    auto ta2_xz_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 171);

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

    auto ta2_yy_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 186);

    auto ta2_yy_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 187);

    auto ta2_yy_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 188);

    auto ta2_yy_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 189);

    auto ta2_yy_xxy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 190);

    auto ta2_yy_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 192);

    auto ta2_yy_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 193);

    auto ta2_yy_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 194);

    auto ta2_yy_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 196);

    auto ta2_yy_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 197);

    auto ta2_yy_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 198);

    auto ta2_yy_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 199);

    auto ta2_yy_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 201);

    auto ta2_yy_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 202);

    auto ta2_yy_xyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 203);

    auto ta2_yy_xyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 208);

    auto ta2_yy_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 210);

    auto ta2_yy_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 212);

    auto ta2_yy_xzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 213);

    auto ta2_yy_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 214);

    auto ta2_yy_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 215);

    auto ta2_yy_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 216);

    auto ta2_yy_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 217);

    auto ta2_yy_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 218);

    auto ta2_yy_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 219);

    auto ta2_yy_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 220);

    auto ta2_yy_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 221);

    auto ta2_yy_yyz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 222);

    auto ta2_yy_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 223);

    auto ta2_yy_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 224);

    auto ta2_yy_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 225);

    auto ta2_yy_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 226);

    auto ta2_yy_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 227);

    auto ta2_yy_yzz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 229);

    auto ta2_yy_yzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 230);

    auto ta2_yy_yzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 231);

    auto ta2_yy_yzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 232);

    auto ta2_yy_yzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 233);

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

    auto ta2_yz_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 246);

    auto ta2_yz_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 247);

    auto ta2_yz_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 248);

    auto ta2_yz_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 249);

    auto ta2_yz_xxy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 250);

    auto ta2_yz_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 252);

    auto ta2_yz_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 253);

    auto ta2_yz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 254);

    auto ta2_yz_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 256);

    auto ta2_yz_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 257);

    auto ta2_yz_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 258);

    auto ta2_yz_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 259);

    auto ta2_yz_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 261);

    auto ta2_yz_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 262);

    auto ta2_yz_xyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 263);

    auto ta2_yz_xyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 268);

    auto ta2_yz_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 270);

    auto ta2_yz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 272);

    auto ta2_yz_xzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 273);

    auto ta2_yz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 274);

    auto ta2_yz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 275);

    auto ta2_yz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 276);

    auto ta2_yz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 277);

    auto ta2_yz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 278);

    auto ta2_yz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 279);

    auto ta2_yz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 280);

    auto ta2_yz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 281);

    auto ta2_yz_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 283);

    auto ta2_yz_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 284);

    auto ta2_yz_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 285);

    auto ta2_yz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 286);

    auto ta2_yz_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 287);

    auto ta2_yz_yzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 288);

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

    auto ta2_zz_xxy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 306);

    auto ta2_zz_xxy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 307);

    auto ta2_zz_xxy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 308);

    auto ta2_zz_xxy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 309);

    auto ta2_zz_xxy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 310);

    auto ta2_zz_xxz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 312);

    auto ta2_zz_xxz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 313);

    auto ta2_zz_xxz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 314);

    auto ta2_zz_xxz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 316);

    auto ta2_zz_xxz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 317);

    auto ta2_zz_xyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 318);

    auto ta2_zz_xyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 319);

    auto ta2_zz_xyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 321);

    auto ta2_zz_xyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 322);

    auto ta2_zz_xyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 323);

    auto ta2_zz_xyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 328);

    auto ta2_zz_xzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 330);

    auto ta2_zz_xzz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 332);

    auto ta2_zz_xzz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 333);

    auto ta2_zz_xzz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 334);

    auto ta2_zz_xzz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 335);

    auto ta2_zz_yyy_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 336);

    auto ta2_zz_yyy_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 337);

    auto ta2_zz_yyy_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 338);

    auto ta2_zz_yyy_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 339);

    auto ta2_zz_yyy_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 340);

    auto ta2_zz_yyy_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 341);

    auto ta2_zz_yyz_xy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 343);

    auto ta2_zz_yyz_xz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 344);

    auto ta2_zz_yyz_yy_1 = pbuffer.data(idx_npot_geom_020_1_fd + 345);

    auto ta2_zz_yyz_yz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 346);

    auto ta2_zz_yyz_zz_1 = pbuffer.data(idx_npot_geom_020_1_fd + 347);

    auto ta2_zz_yzz_xx_1 = pbuffer.data(idx_npot_geom_020_1_fd + 348);

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

    // Set up 0-6 components of targeted buffer : GD

    auto ta2_xx_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd);

    auto ta2_xx_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 1);

    auto ta2_xx_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 2);

    auto ta2_xx_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 3);

    auto ta2_xx_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 4);

    auto ta2_xx_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 5);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxx_yz_1,   \
                             ta1_x_xxx_zz_1,   \
                             ta2_xx_xx_xx_0,   \
                             ta2_xx_xx_xx_1,   \
                             ta2_xx_xx_xy_0,   \
                             ta2_xx_xx_xy_1,   \
                             ta2_xx_xx_xz_0,   \
                             ta2_xx_xx_xz_1,   \
                             ta2_xx_xx_yy_0,   \
                             ta2_xx_xx_yy_1,   \
                             ta2_xx_xx_yz_0,   \
                             ta2_xx_xx_yz_1,   \
                             ta2_xx_xx_zz_0,   \
                             ta2_xx_xx_zz_1,   \
                             ta2_xx_xxx_x_0,   \
                             ta2_xx_xxx_x_1,   \
                             ta2_xx_xxx_xx_0,  \
                             ta2_xx_xxx_xx_1,  \
                             ta2_xx_xxx_xy_0,  \
                             ta2_xx_xxx_xy_1,  \
                             ta2_xx_xxx_xz_0,  \
                             ta2_xx_xxx_xz_1,  \
                             ta2_xx_xxx_y_0,   \
                             ta2_xx_xxx_y_1,   \
                             ta2_xx_xxx_yy_0,  \
                             ta2_xx_xxx_yy_1,  \
                             ta2_xx_xxx_yz_0,  \
                             ta2_xx_xxx_yz_1,  \
                             ta2_xx_xxx_z_0,   \
                             ta2_xx_xxx_z_1,   \
                             ta2_xx_xxx_zz_0,  \
                             ta2_xx_xxx_zz_1,  \
                             ta2_xx_xxxx_xx_0, \
                             ta2_xx_xxxx_xy_0, \
                             ta2_xx_xxxx_xz_0, \
                             ta2_xx_xxxx_yy_0, \
                             ta2_xx_xxxx_yz_0, \
                             ta2_xx_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxx_xx_0[i] = 3.0 * ta2_xx_xx_xx_0[i] * fe_0 - 3.0 * ta2_xx_xx_xx_1[i] * fe_0 + 2.0 * ta2_xx_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_xx_xxx_x_1[i] * fe_0 + 2.0 * ta1_x_xxx_xx_1[i] + ta2_xx_xxx_xx_0[i] * pa_x[i] - ta2_xx_xxx_xx_1[i] * pc_x[i];

        ta2_xx_xxxx_xy_0[i] = 3.0 * ta2_xx_xx_xy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xy_1[i] * fe_0 + ta2_xx_xxx_y_0[i] * fe_0 - ta2_xx_xxx_y_1[i] * fe_0 +
                              2.0 * ta1_x_xxx_xy_1[i] + ta2_xx_xxx_xy_0[i] * pa_x[i] - ta2_xx_xxx_xy_1[i] * pc_x[i];

        ta2_xx_xxxx_xz_0[i] = 3.0 * ta2_xx_xx_xz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xz_1[i] * fe_0 + ta2_xx_xxx_z_0[i] * fe_0 - ta2_xx_xxx_z_1[i] * fe_0 +
                              2.0 * ta1_x_xxx_xz_1[i] + ta2_xx_xxx_xz_0[i] * pa_x[i] - ta2_xx_xxx_xz_1[i] * pc_x[i];

        ta2_xx_xxxx_yy_0[i] = 3.0 * ta2_xx_xx_yy_0[i] * fe_0 - 3.0 * ta2_xx_xx_yy_1[i] * fe_0 + 2.0 * ta1_x_xxx_yy_1[i] +
                              ta2_xx_xxx_yy_0[i] * pa_x[i] - ta2_xx_xxx_yy_1[i] * pc_x[i];

        ta2_xx_xxxx_yz_0[i] = 3.0 * ta2_xx_xx_yz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yz_1[i] +
                              ta2_xx_xxx_yz_0[i] * pa_x[i] - ta2_xx_xxx_yz_1[i] * pc_x[i];

        ta2_xx_xxxx_zz_0[i] = 3.0 * ta2_xx_xx_zz_0[i] * fe_0 - 3.0 * ta2_xx_xx_zz_1[i] * fe_0 + 2.0 * ta1_x_xxx_zz_1[i] +
                              ta2_xx_xxx_zz_0[i] * pa_x[i] - ta2_xx_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto ta2_xx_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 6);

    auto ta2_xx_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 7);

    auto ta2_xx_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 8);

    auto ta2_xx_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 9);

    auto ta2_xx_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 10);

    auto ta2_xx_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 11);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xx_xxx_x_0,   \
                             ta2_xx_xxx_x_1,   \
                             ta2_xx_xxx_xx_0,  \
                             ta2_xx_xxx_xx_1,  \
                             ta2_xx_xxx_xy_0,  \
                             ta2_xx_xxx_xy_1,  \
                             ta2_xx_xxx_xz_0,  \
                             ta2_xx_xxx_xz_1,  \
                             ta2_xx_xxx_y_0,   \
                             ta2_xx_xxx_y_1,   \
                             ta2_xx_xxx_yy_0,  \
                             ta2_xx_xxx_yy_1,  \
                             ta2_xx_xxx_yz_0,  \
                             ta2_xx_xxx_yz_1,  \
                             ta2_xx_xxx_z_0,   \
                             ta2_xx_xxx_z_1,   \
                             ta2_xx_xxx_zz_0,  \
                             ta2_xx_xxx_zz_1,  \
                             ta2_xx_xxxy_xx_0, \
                             ta2_xx_xxxy_xy_0, \
                             ta2_xx_xxxy_xz_0, \
                             ta2_xx_xxxy_yy_0, \
                             ta2_xx_xxxy_yz_0, \
                             ta2_xx_xxxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxy_xx_0[i] = ta2_xx_xxx_xx_0[i] * pa_y[i] - ta2_xx_xxx_xx_1[i] * pc_y[i];

        ta2_xx_xxxy_xy_0[i] = ta2_xx_xxx_x_0[i] * fe_0 - ta2_xx_xxx_x_1[i] * fe_0 + ta2_xx_xxx_xy_0[i] * pa_y[i] - ta2_xx_xxx_xy_1[i] * pc_y[i];

        ta2_xx_xxxy_xz_0[i] = ta2_xx_xxx_xz_0[i] * pa_y[i] - ta2_xx_xxx_xz_1[i] * pc_y[i];

        ta2_xx_xxxy_yy_0[i] =
            2.0 * ta2_xx_xxx_y_0[i] * fe_0 - 2.0 * ta2_xx_xxx_y_1[i] * fe_0 + ta2_xx_xxx_yy_0[i] * pa_y[i] - ta2_xx_xxx_yy_1[i] * pc_y[i];

        ta2_xx_xxxy_yz_0[i] = ta2_xx_xxx_z_0[i] * fe_0 - ta2_xx_xxx_z_1[i] * fe_0 + ta2_xx_xxx_yz_0[i] * pa_y[i] - ta2_xx_xxx_yz_1[i] * pc_y[i];

        ta2_xx_xxxy_zz_0[i] = ta2_xx_xxx_zz_0[i] * pa_y[i] - ta2_xx_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto ta2_xx_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 12);

    auto ta2_xx_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 13);

    auto ta2_xx_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 14);

    auto ta2_xx_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 15);

    auto ta2_xx_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 16);

    auto ta2_xx_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 17);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xx_xxx_x_0,   \
                             ta2_xx_xxx_x_1,   \
                             ta2_xx_xxx_xx_0,  \
                             ta2_xx_xxx_xx_1,  \
                             ta2_xx_xxx_xy_0,  \
                             ta2_xx_xxx_xy_1,  \
                             ta2_xx_xxx_xz_0,  \
                             ta2_xx_xxx_xz_1,  \
                             ta2_xx_xxx_y_0,   \
                             ta2_xx_xxx_y_1,   \
                             ta2_xx_xxx_yy_0,  \
                             ta2_xx_xxx_yy_1,  \
                             ta2_xx_xxx_yz_0,  \
                             ta2_xx_xxx_yz_1,  \
                             ta2_xx_xxx_z_0,   \
                             ta2_xx_xxx_z_1,   \
                             ta2_xx_xxx_zz_0,  \
                             ta2_xx_xxx_zz_1,  \
                             ta2_xx_xxxz_xx_0, \
                             ta2_xx_xxxz_xy_0, \
                             ta2_xx_xxxz_xz_0, \
                             ta2_xx_xxxz_yy_0, \
                             ta2_xx_xxxz_yz_0, \
                             ta2_xx_xxxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxz_xx_0[i] = ta2_xx_xxx_xx_0[i] * pa_z[i] - ta2_xx_xxx_xx_1[i] * pc_z[i];

        ta2_xx_xxxz_xy_0[i] = ta2_xx_xxx_xy_0[i] * pa_z[i] - ta2_xx_xxx_xy_1[i] * pc_z[i];

        ta2_xx_xxxz_xz_0[i] = ta2_xx_xxx_x_0[i] * fe_0 - ta2_xx_xxx_x_1[i] * fe_0 + ta2_xx_xxx_xz_0[i] * pa_z[i] - ta2_xx_xxx_xz_1[i] * pc_z[i];

        ta2_xx_xxxz_yy_0[i] = ta2_xx_xxx_yy_0[i] * pa_z[i] - ta2_xx_xxx_yy_1[i] * pc_z[i];

        ta2_xx_xxxz_yz_0[i] = ta2_xx_xxx_y_0[i] * fe_0 - ta2_xx_xxx_y_1[i] * fe_0 + ta2_xx_xxx_yz_0[i] * pa_z[i] - ta2_xx_xxx_yz_1[i] * pc_z[i];

        ta2_xx_xxxz_zz_0[i] =
            2.0 * ta2_xx_xxx_z_0[i] * fe_0 - 2.0 * ta2_xx_xxx_z_1[i] * fe_0 + ta2_xx_xxx_zz_0[i] * pa_z[i] - ta2_xx_xxx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto ta2_xx_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 18);

    auto ta2_xx_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 19);

    auto ta2_xx_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 20);

    auto ta2_xx_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 21);

    auto ta2_xx_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 22);

    auto ta2_xx_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 23);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xyy_yy_1,   \
                             ta1_x_xyy_yz_1,   \
                             ta2_xx_xx_xx_0,   \
                             ta2_xx_xx_xx_1,   \
                             ta2_xx_xx_xy_0,   \
                             ta2_xx_xx_xy_1,   \
                             ta2_xx_xx_xz_0,   \
                             ta2_xx_xx_xz_1,   \
                             ta2_xx_xx_zz_0,   \
                             ta2_xx_xx_zz_1,   \
                             ta2_xx_xxy_x_0,   \
                             ta2_xx_xxy_x_1,   \
                             ta2_xx_xxy_xx_0,  \
                             ta2_xx_xxy_xx_1,  \
                             ta2_xx_xxy_xy_0,  \
                             ta2_xx_xxy_xy_1,  \
                             ta2_xx_xxy_xz_0,  \
                             ta2_xx_xxy_xz_1,  \
                             ta2_xx_xxy_zz_0,  \
                             ta2_xx_xxy_zz_1,  \
                             ta2_xx_xxyy_xx_0, \
                             ta2_xx_xxyy_xy_0, \
                             ta2_xx_xxyy_xz_0, \
                             ta2_xx_xxyy_yy_0, \
                             ta2_xx_xxyy_yz_0, \
                             ta2_xx_xxyy_zz_0, \
                             ta2_xx_xyy_yy_0,  \
                             ta2_xx_xyy_yy_1,  \
                             ta2_xx_xyy_yz_0,  \
                             ta2_xx_xyy_yz_1,  \
                             ta2_xx_yy_yy_0,   \
                             ta2_xx_yy_yy_1,   \
                             ta2_xx_yy_yz_0,   \
                             ta2_xx_yy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyy_xx_0[i] = ta2_xx_xx_xx_0[i] * fe_0 - ta2_xx_xx_xx_1[i] * fe_0 + ta2_xx_xxy_xx_0[i] * pa_y[i] - ta2_xx_xxy_xx_1[i] * pc_y[i];

        ta2_xx_xxyy_xy_0[i] = ta2_xx_xx_xy_0[i] * fe_0 - ta2_xx_xx_xy_1[i] * fe_0 + ta2_xx_xxy_x_0[i] * fe_0 - ta2_xx_xxy_x_1[i] * fe_0 +
                              ta2_xx_xxy_xy_0[i] * pa_y[i] - ta2_xx_xxy_xy_1[i] * pc_y[i];

        ta2_xx_xxyy_xz_0[i] = ta2_xx_xx_xz_0[i] * fe_0 - ta2_xx_xx_xz_1[i] * fe_0 + ta2_xx_xxy_xz_0[i] * pa_y[i] - ta2_xx_xxy_xz_1[i] * pc_y[i];

        ta2_xx_xxyy_yy_0[i] = ta2_xx_yy_yy_0[i] * fe_0 - ta2_xx_yy_yy_1[i] * fe_0 + 2.0 * ta1_x_xyy_yy_1[i] + ta2_xx_xyy_yy_0[i] * pa_x[i] -
                              ta2_xx_xyy_yy_1[i] * pc_x[i];

        ta2_xx_xxyy_yz_0[i] = ta2_xx_yy_yz_0[i] * fe_0 - ta2_xx_yy_yz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yz_1[i] + ta2_xx_xyy_yz_0[i] * pa_x[i] -
                              ta2_xx_xyy_yz_1[i] * pc_x[i];

        ta2_xx_xxyy_zz_0[i] = ta2_xx_xx_zz_0[i] * fe_0 - ta2_xx_xx_zz_1[i] * fe_0 + ta2_xx_xxy_zz_0[i] * pa_y[i] - ta2_xx_xxy_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto ta2_xx_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 24);

    auto ta2_xx_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 25);

    auto ta2_xx_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 26);

    auto ta2_xx_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 27);

    auto ta2_xx_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 28);

    auto ta2_xx_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 29);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta2_xx_xxy_xy_0,  \
                             ta2_xx_xxy_xy_1,  \
                             ta2_xx_xxy_yy_0,  \
                             ta2_xx_xxy_yy_1,  \
                             ta2_xx_xxyz_xx_0, \
                             ta2_xx_xxyz_xy_0, \
                             ta2_xx_xxyz_xz_0, \
                             ta2_xx_xxyz_yy_0, \
                             ta2_xx_xxyz_yz_0, \
                             ta2_xx_xxyz_zz_0, \
                             ta2_xx_xxz_xx_0,  \
                             ta2_xx_xxz_xx_1,  \
                             ta2_xx_xxz_xz_0,  \
                             ta2_xx_xxz_xz_1,  \
                             ta2_xx_xxz_yz_0,  \
                             ta2_xx_xxz_yz_1,  \
                             ta2_xx_xxz_z_0,   \
                             ta2_xx_xxz_z_1,   \
                             ta2_xx_xxz_zz_0,  \
                             ta2_xx_xxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyz_xx_0[i] = ta2_xx_xxz_xx_0[i] * pa_y[i] - ta2_xx_xxz_xx_1[i] * pc_y[i];

        ta2_xx_xxyz_xy_0[i] = ta2_xx_xxy_xy_0[i] * pa_z[i] - ta2_xx_xxy_xy_1[i] * pc_z[i];

        ta2_xx_xxyz_xz_0[i] = ta2_xx_xxz_xz_0[i] * pa_y[i] - ta2_xx_xxz_xz_1[i] * pc_y[i];

        ta2_xx_xxyz_yy_0[i] = ta2_xx_xxy_yy_0[i] * pa_z[i] - ta2_xx_xxy_yy_1[i] * pc_z[i];

        ta2_xx_xxyz_yz_0[i] = ta2_xx_xxz_z_0[i] * fe_0 - ta2_xx_xxz_z_1[i] * fe_0 + ta2_xx_xxz_yz_0[i] * pa_y[i] - ta2_xx_xxz_yz_1[i] * pc_y[i];

        ta2_xx_xxyz_zz_0[i] = ta2_xx_xxz_zz_0[i] * pa_y[i] - ta2_xx_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto ta2_xx_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 30);

    auto ta2_xx_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 31);

    auto ta2_xx_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 32);

    auto ta2_xx_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 33);

    auto ta2_xx_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 34);

    auto ta2_xx_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 35);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xzz_yz_1,   \
                             ta1_x_xzz_zz_1,   \
                             ta2_xx_xx_xx_0,   \
                             ta2_xx_xx_xx_1,   \
                             ta2_xx_xx_xy_0,   \
                             ta2_xx_xx_xy_1,   \
                             ta2_xx_xx_xz_0,   \
                             ta2_xx_xx_xz_1,   \
                             ta2_xx_xx_yy_0,   \
                             ta2_xx_xx_yy_1,   \
                             ta2_xx_xxz_x_0,   \
                             ta2_xx_xxz_x_1,   \
                             ta2_xx_xxz_xx_0,  \
                             ta2_xx_xxz_xx_1,  \
                             ta2_xx_xxz_xy_0,  \
                             ta2_xx_xxz_xy_1,  \
                             ta2_xx_xxz_xz_0,  \
                             ta2_xx_xxz_xz_1,  \
                             ta2_xx_xxz_yy_0,  \
                             ta2_xx_xxz_yy_1,  \
                             ta2_xx_xxzz_xx_0, \
                             ta2_xx_xxzz_xy_0, \
                             ta2_xx_xxzz_xz_0, \
                             ta2_xx_xxzz_yy_0, \
                             ta2_xx_xxzz_yz_0, \
                             ta2_xx_xxzz_zz_0, \
                             ta2_xx_xzz_yz_0,  \
                             ta2_xx_xzz_yz_1,  \
                             ta2_xx_xzz_zz_0,  \
                             ta2_xx_xzz_zz_1,  \
                             ta2_xx_zz_yz_0,   \
                             ta2_xx_zz_yz_1,   \
                             ta2_xx_zz_zz_0,   \
                             ta2_xx_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxzz_xx_0[i] = ta2_xx_xx_xx_0[i] * fe_0 - ta2_xx_xx_xx_1[i] * fe_0 + ta2_xx_xxz_xx_0[i] * pa_z[i] - ta2_xx_xxz_xx_1[i] * pc_z[i];

        ta2_xx_xxzz_xy_0[i] = ta2_xx_xx_xy_0[i] * fe_0 - ta2_xx_xx_xy_1[i] * fe_0 + ta2_xx_xxz_xy_0[i] * pa_z[i] - ta2_xx_xxz_xy_1[i] * pc_z[i];

        ta2_xx_xxzz_xz_0[i] = ta2_xx_xx_xz_0[i] * fe_0 - ta2_xx_xx_xz_1[i] * fe_0 + ta2_xx_xxz_x_0[i] * fe_0 - ta2_xx_xxz_x_1[i] * fe_0 +
                              ta2_xx_xxz_xz_0[i] * pa_z[i] - ta2_xx_xxz_xz_1[i] * pc_z[i];

        ta2_xx_xxzz_yy_0[i] = ta2_xx_xx_yy_0[i] * fe_0 - ta2_xx_xx_yy_1[i] * fe_0 + ta2_xx_xxz_yy_0[i] * pa_z[i] - ta2_xx_xxz_yy_1[i] * pc_z[i];

        ta2_xx_xxzz_yz_0[i] = ta2_xx_zz_yz_0[i] * fe_0 - ta2_xx_zz_yz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yz_1[i] + ta2_xx_xzz_yz_0[i] * pa_x[i] -
                              ta2_xx_xzz_yz_1[i] * pc_x[i];

        ta2_xx_xxzz_zz_0[i] = ta2_xx_zz_zz_0[i] * fe_0 - ta2_xx_zz_zz_1[i] * fe_0 + 2.0 * ta1_x_xzz_zz_1[i] + ta2_xx_xzz_zz_0[i] * pa_x[i] -
                              ta2_xx_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto ta2_xx_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 36);

    auto ta2_xx_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 37);

    auto ta2_xx_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 38);

    auto ta2_xx_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 39);

    auto ta2_xx_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 40);

    auto ta2_xx_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 41);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_zz_1,   \
                             ta2_xx_xy_xx_0,   \
                             ta2_xx_xy_xx_1,   \
                             ta2_xx_xy_xz_0,   \
                             ta2_xx_xy_xz_1,   \
                             ta2_xx_xyy_xx_0,  \
                             ta2_xx_xyy_xx_1,  \
                             ta2_xx_xyy_xz_0,  \
                             ta2_xx_xyy_xz_1,  \
                             ta2_xx_xyyy_xx_0, \
                             ta2_xx_xyyy_xy_0, \
                             ta2_xx_xyyy_xz_0, \
                             ta2_xx_xyyy_yy_0, \
                             ta2_xx_xyyy_yz_0, \
                             ta2_xx_xyyy_zz_0, \
                             ta2_xx_yyy_xy_0,  \
                             ta2_xx_yyy_xy_1,  \
                             ta2_xx_yyy_y_0,   \
                             ta2_xx_yyy_y_1,   \
                             ta2_xx_yyy_yy_0,  \
                             ta2_xx_yyy_yy_1,  \
                             ta2_xx_yyy_yz_0,  \
                             ta2_xx_yyy_yz_1,  \
                             ta2_xx_yyy_zz_0,  \
                             ta2_xx_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyy_xx_0[i] =
            2.0 * ta2_xx_xy_xx_0[i] * fe_0 - 2.0 * ta2_xx_xy_xx_1[i] * fe_0 + ta2_xx_xyy_xx_0[i] * pa_y[i] - ta2_xx_xyy_xx_1[i] * pc_y[i];

        ta2_xx_xyyy_xy_0[i] = ta2_xx_yyy_y_0[i] * fe_0 - ta2_xx_yyy_y_1[i] * fe_0 + 2.0 * ta1_x_yyy_xy_1[i] + ta2_xx_yyy_xy_0[i] * pa_x[i] -
                              ta2_xx_yyy_xy_1[i] * pc_x[i];

        ta2_xx_xyyy_xz_0[i] =
            2.0 * ta2_xx_xy_xz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xz_1[i] * fe_0 + ta2_xx_xyy_xz_0[i] * pa_y[i] - ta2_xx_xyy_xz_1[i] * pc_y[i];

        ta2_xx_xyyy_yy_0[i] = 2.0 * ta1_x_yyy_yy_1[i] + ta2_xx_yyy_yy_0[i] * pa_x[i] - ta2_xx_yyy_yy_1[i] * pc_x[i];

        ta2_xx_xyyy_yz_0[i] = 2.0 * ta1_x_yyy_yz_1[i] + ta2_xx_yyy_yz_0[i] * pa_x[i] - ta2_xx_yyy_yz_1[i] * pc_x[i];

        ta2_xx_xyyy_zz_0[i] = 2.0 * ta1_x_yyy_zz_1[i] + ta2_xx_yyy_zz_0[i] * pa_x[i] - ta2_xx_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto ta2_xx_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 42);

    auto ta2_xx_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 43);

    auto ta2_xx_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 44);

    auto ta2_xx_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 45);

    auto ta2_xx_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 46);

    auto ta2_xx_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 47);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyz_yz_1,   \
                             ta1_x_yyz_zz_1,   \
                             ta2_xx_xyy_xx_0,  \
                             ta2_xx_xyy_xx_1,  \
                             ta2_xx_xyy_xy_0,  \
                             ta2_xx_xyy_xy_1,  \
                             ta2_xx_xyy_yy_0,  \
                             ta2_xx_xyy_yy_1,  \
                             ta2_xx_xyyz_xx_0, \
                             ta2_xx_xyyz_xy_0, \
                             ta2_xx_xyyz_xz_0, \
                             ta2_xx_xyyz_yy_0, \
                             ta2_xx_xyyz_yz_0, \
                             ta2_xx_xyyz_zz_0, \
                             ta2_xx_xyz_xz_0,  \
                             ta2_xx_xyz_xz_1,  \
                             ta2_xx_xz_xz_0,   \
                             ta2_xx_xz_xz_1,   \
                             ta2_xx_yyz_yz_0,  \
                             ta2_xx_yyz_yz_1,  \
                             ta2_xx_yyz_zz_0,  \
                             ta2_xx_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyz_xx_0[i] = ta2_xx_xyy_xx_0[i] * pa_z[i] - ta2_xx_xyy_xx_1[i] * pc_z[i];

        ta2_xx_xyyz_xy_0[i] = ta2_xx_xyy_xy_0[i] * pa_z[i] - ta2_xx_xyy_xy_1[i] * pc_z[i];

        ta2_xx_xyyz_xz_0[i] = ta2_xx_xz_xz_0[i] * fe_0 - ta2_xx_xz_xz_1[i] * fe_0 + ta2_xx_xyz_xz_0[i] * pa_y[i] - ta2_xx_xyz_xz_1[i] * pc_y[i];

        ta2_xx_xyyz_yy_0[i] = ta2_xx_xyy_yy_0[i] * pa_z[i] - ta2_xx_xyy_yy_1[i] * pc_z[i];

        ta2_xx_xyyz_yz_0[i] = 2.0 * ta1_x_yyz_yz_1[i] + ta2_xx_yyz_yz_0[i] * pa_x[i] - ta2_xx_yyz_yz_1[i] * pc_x[i];

        ta2_xx_xyyz_zz_0[i] = 2.0 * ta1_x_yyz_zz_1[i] + ta2_xx_yyz_zz_0[i] * pa_x[i] - ta2_xx_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto ta2_xx_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 48);

    auto ta2_xx_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 49);

    auto ta2_xx_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 50);

    auto ta2_xx_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 51);

    auto ta2_xx_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 52);

    auto ta2_xx_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 53);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_yzz_yy_1,   \
                             ta1_x_yzz_yz_1,   \
                             ta2_xx_xyzz_xx_0, \
                             ta2_xx_xyzz_xy_0, \
                             ta2_xx_xyzz_xz_0, \
                             ta2_xx_xyzz_yy_0, \
                             ta2_xx_xyzz_yz_0, \
                             ta2_xx_xyzz_zz_0, \
                             ta2_xx_xzz_x_0,   \
                             ta2_xx_xzz_x_1,   \
                             ta2_xx_xzz_xx_0,  \
                             ta2_xx_xzz_xx_1,  \
                             ta2_xx_xzz_xy_0,  \
                             ta2_xx_xzz_xy_1,  \
                             ta2_xx_xzz_xz_0,  \
                             ta2_xx_xzz_xz_1,  \
                             ta2_xx_xzz_zz_0,  \
                             ta2_xx_xzz_zz_1,  \
                             ta2_xx_yzz_yy_0,  \
                             ta2_xx_yzz_yy_1,  \
                             ta2_xx_yzz_yz_0,  \
                             ta2_xx_yzz_yz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyzz_xx_0[i] = ta2_xx_xzz_xx_0[i] * pa_y[i] - ta2_xx_xzz_xx_1[i] * pc_y[i];

        ta2_xx_xyzz_xy_0[i] = ta2_xx_xzz_x_0[i] * fe_0 - ta2_xx_xzz_x_1[i] * fe_0 + ta2_xx_xzz_xy_0[i] * pa_y[i] - ta2_xx_xzz_xy_1[i] * pc_y[i];

        ta2_xx_xyzz_xz_0[i] = ta2_xx_xzz_xz_0[i] * pa_y[i] - ta2_xx_xzz_xz_1[i] * pc_y[i];

        ta2_xx_xyzz_yy_0[i] = 2.0 * ta1_x_yzz_yy_1[i] + ta2_xx_yzz_yy_0[i] * pa_x[i] - ta2_xx_yzz_yy_1[i] * pc_x[i];

        ta2_xx_xyzz_yz_0[i] = 2.0 * ta1_x_yzz_yz_1[i] + ta2_xx_yzz_yz_0[i] * pa_x[i] - ta2_xx_yzz_yz_1[i] * pc_x[i];

        ta2_xx_xyzz_zz_0[i] = ta2_xx_xzz_zz_0[i] * pa_y[i] - ta2_xx_xzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto ta2_xx_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 54);

    auto ta2_xx_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 55);

    auto ta2_xx_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 56);

    auto ta2_xx_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 57);

    auto ta2_xx_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 58);

    auto ta2_xx_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_yy_1,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_1,   \
                             ta2_xx_xz_xx_0,   \
                             ta2_xx_xz_xx_1,   \
                             ta2_xx_xz_xy_0,   \
                             ta2_xx_xz_xy_1,   \
                             ta2_xx_xzz_xx_0,  \
                             ta2_xx_xzz_xx_1,  \
                             ta2_xx_xzz_xy_0,  \
                             ta2_xx_xzz_xy_1,  \
                             ta2_xx_xzzz_xx_0, \
                             ta2_xx_xzzz_xy_0, \
                             ta2_xx_xzzz_xz_0, \
                             ta2_xx_xzzz_yy_0, \
                             ta2_xx_xzzz_yz_0, \
                             ta2_xx_xzzz_zz_0, \
                             ta2_xx_zzz_xz_0,  \
                             ta2_xx_zzz_xz_1,  \
                             ta2_xx_zzz_yy_0,  \
                             ta2_xx_zzz_yy_1,  \
                             ta2_xx_zzz_yz_0,  \
                             ta2_xx_zzz_yz_1,  \
                             ta2_xx_zzz_z_0,   \
                             ta2_xx_zzz_z_1,   \
                             ta2_xx_zzz_zz_0,  \
                             ta2_xx_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzzz_xx_0[i] =
            2.0 * ta2_xx_xz_xx_0[i] * fe_0 - 2.0 * ta2_xx_xz_xx_1[i] * fe_0 + ta2_xx_xzz_xx_0[i] * pa_z[i] - ta2_xx_xzz_xx_1[i] * pc_z[i];

        ta2_xx_xzzz_xy_0[i] =
            2.0 * ta2_xx_xz_xy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xy_1[i] * fe_0 + ta2_xx_xzz_xy_0[i] * pa_z[i] - ta2_xx_xzz_xy_1[i] * pc_z[i];

        ta2_xx_xzzz_xz_0[i] = ta2_xx_zzz_z_0[i] * fe_0 - ta2_xx_zzz_z_1[i] * fe_0 + 2.0 * ta1_x_zzz_xz_1[i] + ta2_xx_zzz_xz_0[i] * pa_x[i] -
                              ta2_xx_zzz_xz_1[i] * pc_x[i];

        ta2_xx_xzzz_yy_0[i] = 2.0 * ta1_x_zzz_yy_1[i] + ta2_xx_zzz_yy_0[i] * pa_x[i] - ta2_xx_zzz_yy_1[i] * pc_x[i];

        ta2_xx_xzzz_yz_0[i] = 2.0 * ta1_x_zzz_yz_1[i] + ta2_xx_zzz_yz_0[i] * pa_x[i] - ta2_xx_zzz_yz_1[i] * pc_x[i];

        ta2_xx_xzzz_zz_0[i] = 2.0 * ta1_x_zzz_zz_1[i] + ta2_xx_zzz_zz_0[i] * pa_x[i] - ta2_xx_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto ta2_xx_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 60);

    auto ta2_xx_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 61);

    auto ta2_xx_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 62);

    auto ta2_xx_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 63);

    auto ta2_xx_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 64);

    auto ta2_xx_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 65);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xx_yy_xx_0,   \
                             ta2_xx_yy_xx_1,   \
                             ta2_xx_yy_xy_0,   \
                             ta2_xx_yy_xy_1,   \
                             ta2_xx_yy_xz_0,   \
                             ta2_xx_yy_xz_1,   \
                             ta2_xx_yy_yy_0,   \
                             ta2_xx_yy_yy_1,   \
                             ta2_xx_yy_yz_0,   \
                             ta2_xx_yy_yz_1,   \
                             ta2_xx_yy_zz_0,   \
                             ta2_xx_yy_zz_1,   \
                             ta2_xx_yyy_x_0,   \
                             ta2_xx_yyy_x_1,   \
                             ta2_xx_yyy_xx_0,  \
                             ta2_xx_yyy_xx_1,  \
                             ta2_xx_yyy_xy_0,  \
                             ta2_xx_yyy_xy_1,  \
                             ta2_xx_yyy_xz_0,  \
                             ta2_xx_yyy_xz_1,  \
                             ta2_xx_yyy_y_0,   \
                             ta2_xx_yyy_y_1,   \
                             ta2_xx_yyy_yy_0,  \
                             ta2_xx_yyy_yy_1,  \
                             ta2_xx_yyy_yz_0,  \
                             ta2_xx_yyy_yz_1,  \
                             ta2_xx_yyy_z_0,   \
                             ta2_xx_yyy_z_1,   \
                             ta2_xx_yyy_zz_0,  \
                             ta2_xx_yyy_zz_1,  \
                             ta2_xx_yyyy_xx_0, \
                             ta2_xx_yyyy_xy_0, \
                             ta2_xx_yyyy_xz_0, \
                             ta2_xx_yyyy_yy_0, \
                             ta2_xx_yyyy_yz_0, \
                             ta2_xx_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyy_xx_0[i] =
            3.0 * ta2_xx_yy_xx_0[i] * fe_0 - 3.0 * ta2_xx_yy_xx_1[i] * fe_0 + ta2_xx_yyy_xx_0[i] * pa_y[i] - ta2_xx_yyy_xx_1[i] * pc_y[i];

        ta2_xx_yyyy_xy_0[i] = 3.0 * ta2_xx_yy_xy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xy_1[i] * fe_0 + ta2_xx_yyy_x_0[i] * fe_0 - ta2_xx_yyy_x_1[i] * fe_0 +
                              ta2_xx_yyy_xy_0[i] * pa_y[i] - ta2_xx_yyy_xy_1[i] * pc_y[i];

        ta2_xx_yyyy_xz_0[i] =
            3.0 * ta2_xx_yy_xz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xz_1[i] * fe_0 + ta2_xx_yyy_xz_0[i] * pa_y[i] - ta2_xx_yyy_xz_1[i] * pc_y[i];

        ta2_xx_yyyy_yy_0[i] = 3.0 * ta2_xx_yy_yy_0[i] * fe_0 - 3.0 * ta2_xx_yy_yy_1[i] * fe_0 + 2.0 * ta2_xx_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_xx_yyy_y_1[i] * fe_0 + ta2_xx_yyy_yy_0[i] * pa_y[i] - ta2_xx_yyy_yy_1[i] * pc_y[i];

        ta2_xx_yyyy_yz_0[i] = 3.0 * ta2_xx_yy_yz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yz_1[i] * fe_0 + ta2_xx_yyy_z_0[i] * fe_0 - ta2_xx_yyy_z_1[i] * fe_0 +
                              ta2_xx_yyy_yz_0[i] * pa_y[i] - ta2_xx_yyy_yz_1[i] * pc_y[i];

        ta2_xx_yyyy_zz_0[i] =
            3.0 * ta2_xx_yy_zz_0[i] * fe_0 - 3.0 * ta2_xx_yy_zz_1[i] * fe_0 + ta2_xx_yyy_zz_0[i] * pa_y[i] - ta2_xx_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto ta2_xx_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 66);

    auto ta2_xx_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 67);

    auto ta2_xx_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 68);

    auto ta2_xx_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 69);

    auto ta2_xx_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 70);

    auto ta2_xx_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 71);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta2_xx_yyy_xx_0,  \
                             ta2_xx_yyy_xx_1,  \
                             ta2_xx_yyy_xy_0,  \
                             ta2_xx_yyy_xy_1,  \
                             ta2_xx_yyy_y_0,   \
                             ta2_xx_yyy_y_1,   \
                             ta2_xx_yyy_yy_0,  \
                             ta2_xx_yyy_yy_1,  \
                             ta2_xx_yyy_yz_0,  \
                             ta2_xx_yyy_yz_1,  \
                             ta2_xx_yyyz_xx_0, \
                             ta2_xx_yyyz_xy_0, \
                             ta2_xx_yyyz_xz_0, \
                             ta2_xx_yyyz_yy_0, \
                             ta2_xx_yyyz_yz_0, \
                             ta2_xx_yyyz_zz_0, \
                             ta2_xx_yyz_xz_0,  \
                             ta2_xx_yyz_xz_1,  \
                             ta2_xx_yyz_zz_0,  \
                             ta2_xx_yyz_zz_1,  \
                             ta2_xx_yz_xz_0,   \
                             ta2_xx_yz_xz_1,   \
                             ta2_xx_yz_zz_0,   \
                             ta2_xx_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyz_xx_0[i] = ta2_xx_yyy_xx_0[i] * pa_z[i] - ta2_xx_yyy_xx_1[i] * pc_z[i];

        ta2_xx_yyyz_xy_0[i] = ta2_xx_yyy_xy_0[i] * pa_z[i] - ta2_xx_yyy_xy_1[i] * pc_z[i];

        ta2_xx_yyyz_xz_0[i] =
            2.0 * ta2_xx_yz_xz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xz_1[i] * fe_0 + ta2_xx_yyz_xz_0[i] * pa_y[i] - ta2_xx_yyz_xz_1[i] * pc_y[i];

        ta2_xx_yyyz_yy_0[i] = ta2_xx_yyy_yy_0[i] * pa_z[i] - ta2_xx_yyy_yy_1[i] * pc_z[i];

        ta2_xx_yyyz_yz_0[i] = ta2_xx_yyy_y_0[i] * fe_0 - ta2_xx_yyy_y_1[i] * fe_0 + ta2_xx_yyy_yz_0[i] * pa_z[i] - ta2_xx_yyy_yz_1[i] * pc_z[i];

        ta2_xx_yyyz_zz_0[i] =
            2.0 * ta2_xx_yz_zz_0[i] * fe_0 - 2.0 * ta2_xx_yz_zz_1[i] * fe_0 + ta2_xx_yyz_zz_0[i] * pa_y[i] - ta2_xx_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto ta2_xx_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 72);

    auto ta2_xx_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 73);

    auto ta2_xx_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 74);

    auto ta2_xx_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 75);

    auto ta2_xx_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 76);

    auto ta2_xx_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 77);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta2_xx_yy_xy_0,   \
                             ta2_xx_yy_xy_1,   \
                             ta2_xx_yy_yy_0,   \
                             ta2_xx_yy_yy_1,   \
                             ta2_xx_yyz_xy_0,  \
                             ta2_xx_yyz_xy_1,  \
                             ta2_xx_yyz_yy_0,  \
                             ta2_xx_yyz_yy_1,  \
                             ta2_xx_yyzz_xx_0, \
                             ta2_xx_yyzz_xy_0, \
                             ta2_xx_yyzz_xz_0, \
                             ta2_xx_yyzz_yy_0, \
                             ta2_xx_yyzz_yz_0, \
                             ta2_xx_yyzz_zz_0, \
                             ta2_xx_yzz_xx_0,  \
                             ta2_xx_yzz_xx_1,  \
                             ta2_xx_yzz_xz_0,  \
                             ta2_xx_yzz_xz_1,  \
                             ta2_xx_yzz_yz_0,  \
                             ta2_xx_yzz_yz_1,  \
                             ta2_xx_yzz_z_0,   \
                             ta2_xx_yzz_z_1,   \
                             ta2_xx_yzz_zz_0,  \
                             ta2_xx_yzz_zz_1,  \
                             ta2_xx_zz_xx_0,   \
                             ta2_xx_zz_xx_1,   \
                             ta2_xx_zz_xz_0,   \
                             ta2_xx_zz_xz_1,   \
                             ta2_xx_zz_yz_0,   \
                             ta2_xx_zz_yz_1,   \
                             ta2_xx_zz_zz_0,   \
                             ta2_xx_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyzz_xx_0[i] = ta2_xx_zz_xx_0[i] * fe_0 - ta2_xx_zz_xx_1[i] * fe_0 + ta2_xx_yzz_xx_0[i] * pa_y[i] - ta2_xx_yzz_xx_1[i] * pc_y[i];

        ta2_xx_yyzz_xy_0[i] = ta2_xx_yy_xy_0[i] * fe_0 - ta2_xx_yy_xy_1[i] * fe_0 + ta2_xx_yyz_xy_0[i] * pa_z[i] - ta2_xx_yyz_xy_1[i] * pc_z[i];

        ta2_xx_yyzz_xz_0[i] = ta2_xx_zz_xz_0[i] * fe_0 - ta2_xx_zz_xz_1[i] * fe_0 + ta2_xx_yzz_xz_0[i] * pa_y[i] - ta2_xx_yzz_xz_1[i] * pc_y[i];

        ta2_xx_yyzz_yy_0[i] = ta2_xx_yy_yy_0[i] * fe_0 - ta2_xx_yy_yy_1[i] * fe_0 + ta2_xx_yyz_yy_0[i] * pa_z[i] - ta2_xx_yyz_yy_1[i] * pc_z[i];

        ta2_xx_yyzz_yz_0[i] = ta2_xx_zz_yz_0[i] * fe_0 - ta2_xx_zz_yz_1[i] * fe_0 + ta2_xx_yzz_z_0[i] * fe_0 - ta2_xx_yzz_z_1[i] * fe_0 +
                              ta2_xx_yzz_yz_0[i] * pa_y[i] - ta2_xx_yzz_yz_1[i] * pc_y[i];

        ta2_xx_yyzz_zz_0[i] = ta2_xx_zz_zz_0[i] * fe_0 - ta2_xx_zz_zz_1[i] * fe_0 + ta2_xx_yzz_zz_0[i] * pa_y[i] - ta2_xx_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto ta2_xx_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 78);

    auto ta2_xx_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 79);

    auto ta2_xx_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 80);

    auto ta2_xx_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 81);

    auto ta2_xx_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 82);

    auto ta2_xx_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 83);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xx_yzzz_xx_0, \
                             ta2_xx_yzzz_xy_0, \
                             ta2_xx_yzzz_xz_0, \
                             ta2_xx_yzzz_yy_0, \
                             ta2_xx_yzzz_yz_0, \
                             ta2_xx_yzzz_zz_0, \
                             ta2_xx_zzz_x_0,   \
                             ta2_xx_zzz_x_1,   \
                             ta2_xx_zzz_xx_0,  \
                             ta2_xx_zzz_xx_1,  \
                             ta2_xx_zzz_xy_0,  \
                             ta2_xx_zzz_xy_1,  \
                             ta2_xx_zzz_xz_0,  \
                             ta2_xx_zzz_xz_1,  \
                             ta2_xx_zzz_y_0,   \
                             ta2_xx_zzz_y_1,   \
                             ta2_xx_zzz_yy_0,  \
                             ta2_xx_zzz_yy_1,  \
                             ta2_xx_zzz_yz_0,  \
                             ta2_xx_zzz_yz_1,  \
                             ta2_xx_zzz_z_0,   \
                             ta2_xx_zzz_z_1,   \
                             ta2_xx_zzz_zz_0,  \
                             ta2_xx_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzzz_xx_0[i] = ta2_xx_zzz_xx_0[i] * pa_y[i] - ta2_xx_zzz_xx_1[i] * pc_y[i];

        ta2_xx_yzzz_xy_0[i] = ta2_xx_zzz_x_0[i] * fe_0 - ta2_xx_zzz_x_1[i] * fe_0 + ta2_xx_zzz_xy_0[i] * pa_y[i] - ta2_xx_zzz_xy_1[i] * pc_y[i];

        ta2_xx_yzzz_xz_0[i] = ta2_xx_zzz_xz_0[i] * pa_y[i] - ta2_xx_zzz_xz_1[i] * pc_y[i];

        ta2_xx_yzzz_yy_0[i] =
            2.0 * ta2_xx_zzz_y_0[i] * fe_0 - 2.0 * ta2_xx_zzz_y_1[i] * fe_0 + ta2_xx_zzz_yy_0[i] * pa_y[i] - ta2_xx_zzz_yy_1[i] * pc_y[i];

        ta2_xx_yzzz_yz_0[i] = ta2_xx_zzz_z_0[i] * fe_0 - ta2_xx_zzz_z_1[i] * fe_0 + ta2_xx_zzz_yz_0[i] * pa_y[i] - ta2_xx_zzz_yz_1[i] * pc_y[i];

        ta2_xx_yzzz_zz_0[i] = ta2_xx_zzz_zz_0[i] * pa_y[i] - ta2_xx_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto ta2_xx_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 84);

    auto ta2_xx_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 85);

    auto ta2_xx_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 86);

    auto ta2_xx_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 87);

    auto ta2_xx_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 88);

    auto ta2_xx_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 89);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xx_zz_xx_0,   \
                             ta2_xx_zz_xx_1,   \
                             ta2_xx_zz_xy_0,   \
                             ta2_xx_zz_xy_1,   \
                             ta2_xx_zz_xz_0,   \
                             ta2_xx_zz_xz_1,   \
                             ta2_xx_zz_yy_0,   \
                             ta2_xx_zz_yy_1,   \
                             ta2_xx_zz_yz_0,   \
                             ta2_xx_zz_yz_1,   \
                             ta2_xx_zz_zz_0,   \
                             ta2_xx_zz_zz_1,   \
                             ta2_xx_zzz_x_0,   \
                             ta2_xx_zzz_x_1,   \
                             ta2_xx_zzz_xx_0,  \
                             ta2_xx_zzz_xx_1,  \
                             ta2_xx_zzz_xy_0,  \
                             ta2_xx_zzz_xy_1,  \
                             ta2_xx_zzz_xz_0,  \
                             ta2_xx_zzz_xz_1,  \
                             ta2_xx_zzz_y_0,   \
                             ta2_xx_zzz_y_1,   \
                             ta2_xx_zzz_yy_0,  \
                             ta2_xx_zzz_yy_1,  \
                             ta2_xx_zzz_yz_0,  \
                             ta2_xx_zzz_yz_1,  \
                             ta2_xx_zzz_z_0,   \
                             ta2_xx_zzz_z_1,   \
                             ta2_xx_zzz_zz_0,  \
                             ta2_xx_zzz_zz_1,  \
                             ta2_xx_zzzz_xx_0, \
                             ta2_xx_zzzz_xy_0, \
                             ta2_xx_zzzz_xz_0, \
                             ta2_xx_zzzz_yy_0, \
                             ta2_xx_zzzz_yz_0, \
                             ta2_xx_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzzz_xx_0[i] =
            3.0 * ta2_xx_zz_xx_0[i] * fe_0 - 3.0 * ta2_xx_zz_xx_1[i] * fe_0 + ta2_xx_zzz_xx_0[i] * pa_z[i] - ta2_xx_zzz_xx_1[i] * pc_z[i];

        ta2_xx_zzzz_xy_0[i] =
            3.0 * ta2_xx_zz_xy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xy_1[i] * fe_0 + ta2_xx_zzz_xy_0[i] * pa_z[i] - ta2_xx_zzz_xy_1[i] * pc_z[i];

        ta2_xx_zzzz_xz_0[i] = 3.0 * ta2_xx_zz_xz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xz_1[i] * fe_0 + ta2_xx_zzz_x_0[i] * fe_0 - ta2_xx_zzz_x_1[i] * fe_0 +
                              ta2_xx_zzz_xz_0[i] * pa_z[i] - ta2_xx_zzz_xz_1[i] * pc_z[i];

        ta2_xx_zzzz_yy_0[i] =
            3.0 * ta2_xx_zz_yy_0[i] * fe_0 - 3.0 * ta2_xx_zz_yy_1[i] * fe_0 + ta2_xx_zzz_yy_0[i] * pa_z[i] - ta2_xx_zzz_yy_1[i] * pc_z[i];

        ta2_xx_zzzz_yz_0[i] = 3.0 * ta2_xx_zz_yz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yz_1[i] * fe_0 + ta2_xx_zzz_y_0[i] * fe_0 - ta2_xx_zzz_y_1[i] * fe_0 +
                              ta2_xx_zzz_yz_0[i] * pa_z[i] - ta2_xx_zzz_yz_1[i] * pc_z[i];

        ta2_xx_zzzz_zz_0[i] = 3.0 * ta2_xx_zz_zz_0[i] * fe_0 - 3.0 * ta2_xx_zz_zz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_xx_zzz_z_1[i] * fe_0 + ta2_xx_zzz_zz_0[i] * pa_z[i] - ta2_xx_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 90-96 components of targeted buffer : GD

    auto ta2_xy_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 90);

    auto ta2_xy_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 91);

    auto ta2_xy_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 92);

    auto ta2_xy_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 93);

    auto ta2_xy_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 94);

    auto ta2_xy_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 95);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_yy_1,   \
                             ta1_y_xxx_yz_1,   \
                             ta1_y_xxx_zz_1,   \
                             ta2_xy_xx_xx_0,   \
                             ta2_xy_xx_xx_1,   \
                             ta2_xy_xx_xy_0,   \
                             ta2_xy_xx_xy_1,   \
                             ta2_xy_xx_xz_0,   \
                             ta2_xy_xx_xz_1,   \
                             ta2_xy_xx_yy_0,   \
                             ta2_xy_xx_yy_1,   \
                             ta2_xy_xx_yz_0,   \
                             ta2_xy_xx_yz_1,   \
                             ta2_xy_xx_zz_0,   \
                             ta2_xy_xx_zz_1,   \
                             ta2_xy_xxx_x_0,   \
                             ta2_xy_xxx_x_1,   \
                             ta2_xy_xxx_xx_0,  \
                             ta2_xy_xxx_xx_1,  \
                             ta2_xy_xxx_xy_0,  \
                             ta2_xy_xxx_xy_1,  \
                             ta2_xy_xxx_xz_0,  \
                             ta2_xy_xxx_xz_1,  \
                             ta2_xy_xxx_y_0,   \
                             ta2_xy_xxx_y_1,   \
                             ta2_xy_xxx_yy_0,  \
                             ta2_xy_xxx_yy_1,  \
                             ta2_xy_xxx_yz_0,  \
                             ta2_xy_xxx_yz_1,  \
                             ta2_xy_xxx_z_0,   \
                             ta2_xy_xxx_z_1,   \
                             ta2_xy_xxx_zz_0,  \
                             ta2_xy_xxx_zz_1,  \
                             ta2_xy_xxxx_xx_0, \
                             ta2_xy_xxxx_xy_0, \
                             ta2_xy_xxxx_xz_0, \
                             ta2_xy_xxxx_yy_0, \
                             ta2_xy_xxxx_yz_0, \
                             ta2_xy_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxx_xx_0[i] = 3.0 * ta2_xy_xx_xx_0[i] * fe_0 - 3.0 * ta2_xy_xx_xx_1[i] * fe_0 + 2.0 * ta2_xy_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_xy_xxx_x_1[i] * fe_0 + ta1_y_xxx_xx_1[i] + ta2_xy_xxx_xx_0[i] * pa_x[i] - ta2_xy_xxx_xx_1[i] * pc_x[i];

        ta2_xy_xxxx_xy_0[i] = 3.0 * ta2_xy_xx_xy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xy_1[i] * fe_0 + ta2_xy_xxx_y_0[i] * fe_0 - ta2_xy_xxx_y_1[i] * fe_0 +
                              ta1_y_xxx_xy_1[i] + ta2_xy_xxx_xy_0[i] * pa_x[i] - ta2_xy_xxx_xy_1[i] * pc_x[i];

        ta2_xy_xxxx_xz_0[i] = 3.0 * ta2_xy_xx_xz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xz_1[i] * fe_0 + ta2_xy_xxx_z_0[i] * fe_0 - ta2_xy_xxx_z_1[i] * fe_0 +
                              ta1_y_xxx_xz_1[i] + ta2_xy_xxx_xz_0[i] * pa_x[i] - ta2_xy_xxx_xz_1[i] * pc_x[i];

        ta2_xy_xxxx_yy_0[i] = 3.0 * ta2_xy_xx_yy_0[i] * fe_0 - 3.0 * ta2_xy_xx_yy_1[i] * fe_0 + ta1_y_xxx_yy_1[i] + ta2_xy_xxx_yy_0[i] * pa_x[i] -
                              ta2_xy_xxx_yy_1[i] * pc_x[i];

        ta2_xy_xxxx_yz_0[i] = 3.0 * ta2_xy_xx_yz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yz_1[i] * fe_0 + ta1_y_xxx_yz_1[i] + ta2_xy_xxx_yz_0[i] * pa_x[i] -
                              ta2_xy_xxx_yz_1[i] * pc_x[i];

        ta2_xy_xxxx_zz_0[i] = 3.0 * ta2_xy_xx_zz_0[i] * fe_0 - 3.0 * ta2_xy_xx_zz_1[i] * fe_0 + ta1_y_xxx_zz_1[i] + ta2_xy_xxx_zz_0[i] * pa_x[i] -
                              ta2_xy_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : GD

    auto ta2_xy_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 96);

    auto ta2_xy_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 97);

    auto ta2_xy_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 98);

    auto ta2_xy_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 99);

    auto ta2_xy_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 100);

    auto ta2_xy_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 101);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_y_xxy_yy_1,   \
                             ta1_y_xxy_yz_1,   \
                             ta2_xy_xxx_x_0,   \
                             ta2_xy_xxx_x_1,   \
                             ta2_xy_xxx_xx_0,  \
                             ta2_xy_xxx_xx_1,  \
                             ta2_xy_xxx_xy_0,  \
                             ta2_xy_xxx_xy_1,  \
                             ta2_xy_xxx_xz_0,  \
                             ta2_xy_xxx_xz_1,  \
                             ta2_xy_xxx_zz_0,  \
                             ta2_xy_xxx_zz_1,  \
                             ta2_xy_xxxy_xx_0, \
                             ta2_xy_xxxy_xy_0, \
                             ta2_xy_xxxy_xz_0, \
                             ta2_xy_xxxy_yy_0, \
                             ta2_xy_xxxy_yz_0, \
                             ta2_xy_xxxy_zz_0, \
                             ta2_xy_xxy_yy_0,  \
                             ta2_xy_xxy_yy_1,  \
                             ta2_xy_xxy_yz_0,  \
                             ta2_xy_xxy_yz_1,  \
                             ta2_xy_xy_yy_0,   \
                             ta2_xy_xy_yy_1,   \
                             ta2_xy_xy_yz_0,   \
                             ta2_xy_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxy_xx_0[i] = ta1_x_xxx_xx_1[i] + ta2_xy_xxx_xx_0[i] * pa_y[i] - ta2_xy_xxx_xx_1[i] * pc_y[i];

        ta2_xy_xxxy_xy_0[i] =
            ta2_xy_xxx_x_0[i] * fe_0 - ta2_xy_xxx_x_1[i] * fe_0 + ta1_x_xxx_xy_1[i] + ta2_xy_xxx_xy_0[i] * pa_y[i] - ta2_xy_xxx_xy_1[i] * pc_y[i];

        ta2_xy_xxxy_xz_0[i] = ta1_x_xxx_xz_1[i] + ta2_xy_xxx_xz_0[i] * pa_y[i] - ta2_xy_xxx_xz_1[i] * pc_y[i];

        ta2_xy_xxxy_yy_0[i] = 2.0 * ta2_xy_xy_yy_0[i] * fe_0 - 2.0 * ta2_xy_xy_yy_1[i] * fe_0 + ta1_y_xxy_yy_1[i] + ta2_xy_xxy_yy_0[i] * pa_x[i] -
                              ta2_xy_xxy_yy_1[i] * pc_x[i];

        ta2_xy_xxxy_yz_0[i] = 2.0 * ta2_xy_xy_yz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yz_1[i] * fe_0 + ta1_y_xxy_yz_1[i] + ta2_xy_xxy_yz_0[i] * pa_x[i] -
                              ta2_xy_xxy_yz_1[i] * pc_x[i];

        ta2_xy_xxxy_zz_0[i] = ta1_x_xxx_zz_1[i] + ta2_xy_xxx_zz_0[i] * pa_y[i] - ta2_xy_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : GD

    auto ta2_xy_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 102);

    auto ta2_xy_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 103);

    auto ta2_xy_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 104);

    auto ta2_xy_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 105);

    auto ta2_xy_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 106);

    auto ta2_xy_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 107);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xy_xxx_x_0,   \
                             ta2_xy_xxx_x_1,   \
                             ta2_xy_xxx_xx_0,  \
                             ta2_xy_xxx_xx_1,  \
                             ta2_xy_xxx_xy_0,  \
                             ta2_xy_xxx_xy_1,  \
                             ta2_xy_xxx_xz_0,  \
                             ta2_xy_xxx_xz_1,  \
                             ta2_xy_xxx_y_0,   \
                             ta2_xy_xxx_y_1,   \
                             ta2_xy_xxx_yy_0,  \
                             ta2_xy_xxx_yy_1,  \
                             ta2_xy_xxx_yz_0,  \
                             ta2_xy_xxx_yz_1,  \
                             ta2_xy_xxx_z_0,   \
                             ta2_xy_xxx_z_1,   \
                             ta2_xy_xxx_zz_0,  \
                             ta2_xy_xxx_zz_1,  \
                             ta2_xy_xxxz_xx_0, \
                             ta2_xy_xxxz_xy_0, \
                             ta2_xy_xxxz_xz_0, \
                             ta2_xy_xxxz_yy_0, \
                             ta2_xy_xxxz_yz_0, \
                             ta2_xy_xxxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxz_xx_0[i] = ta2_xy_xxx_xx_0[i] * pa_z[i] - ta2_xy_xxx_xx_1[i] * pc_z[i];

        ta2_xy_xxxz_xy_0[i] = ta2_xy_xxx_xy_0[i] * pa_z[i] - ta2_xy_xxx_xy_1[i] * pc_z[i];

        ta2_xy_xxxz_xz_0[i] = ta2_xy_xxx_x_0[i] * fe_0 - ta2_xy_xxx_x_1[i] * fe_0 + ta2_xy_xxx_xz_0[i] * pa_z[i] - ta2_xy_xxx_xz_1[i] * pc_z[i];

        ta2_xy_xxxz_yy_0[i] = ta2_xy_xxx_yy_0[i] * pa_z[i] - ta2_xy_xxx_yy_1[i] * pc_z[i];

        ta2_xy_xxxz_yz_0[i] = ta2_xy_xxx_y_0[i] * fe_0 - ta2_xy_xxx_y_1[i] * fe_0 + ta2_xy_xxx_yz_0[i] * pa_z[i] - ta2_xy_xxx_yz_1[i] * pc_z[i];

        ta2_xy_xxxz_zz_0[i] =
            2.0 * ta2_xy_xxx_z_0[i] * fe_0 - 2.0 * ta2_xy_xxx_z_1[i] * fe_0 + ta2_xy_xxx_zz_0[i] * pa_z[i] - ta2_xy_xxx_zz_1[i] * pc_z[i];
    }

    // Set up 108-114 components of targeted buffer : GD

    auto ta2_xy_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 108);

    auto ta2_xy_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 109);

    auto ta2_xy_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 110);

    auto ta2_xy_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 111);

    auto ta2_xy_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 112);

    auto ta2_xy_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 113);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xxy_xx_1,   \
                             ta1_x_xxy_xz_1,   \
                             ta1_y_xyy_xy_1,   \
                             ta1_y_xyy_yy_1,   \
                             ta1_y_xyy_yz_1,   \
                             ta1_y_xyy_zz_1,   \
                             ta2_xy_xx_xx_0,   \
                             ta2_xy_xx_xx_1,   \
                             ta2_xy_xx_xz_0,   \
                             ta2_xy_xx_xz_1,   \
                             ta2_xy_xxy_xx_0,  \
                             ta2_xy_xxy_xx_1,  \
                             ta2_xy_xxy_xz_0,  \
                             ta2_xy_xxy_xz_1,  \
                             ta2_xy_xxyy_xx_0, \
                             ta2_xy_xxyy_xy_0, \
                             ta2_xy_xxyy_xz_0, \
                             ta2_xy_xxyy_yy_0, \
                             ta2_xy_xxyy_yz_0, \
                             ta2_xy_xxyy_zz_0, \
                             ta2_xy_xyy_xy_0,  \
                             ta2_xy_xyy_xy_1,  \
                             ta2_xy_xyy_y_0,   \
                             ta2_xy_xyy_y_1,   \
                             ta2_xy_xyy_yy_0,  \
                             ta2_xy_xyy_yy_1,  \
                             ta2_xy_xyy_yz_0,  \
                             ta2_xy_xyy_yz_1,  \
                             ta2_xy_xyy_zz_0,  \
                             ta2_xy_xyy_zz_1,  \
                             ta2_xy_yy_xy_0,   \
                             ta2_xy_yy_xy_1,   \
                             ta2_xy_yy_yy_0,   \
                             ta2_xy_yy_yy_1,   \
                             ta2_xy_yy_yz_0,   \
                             ta2_xy_yy_yz_1,   \
                             ta2_xy_yy_zz_0,   \
                             ta2_xy_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyy_xx_0[i] =
            ta2_xy_xx_xx_0[i] * fe_0 - ta2_xy_xx_xx_1[i] * fe_0 + ta1_x_xxy_xx_1[i] + ta2_xy_xxy_xx_0[i] * pa_y[i] - ta2_xy_xxy_xx_1[i] * pc_y[i];

        ta2_xy_xxyy_xy_0[i] = ta2_xy_yy_xy_0[i] * fe_0 - ta2_xy_yy_xy_1[i] * fe_0 + ta2_xy_xyy_y_0[i] * fe_0 - ta2_xy_xyy_y_1[i] * fe_0 +
                              ta1_y_xyy_xy_1[i] + ta2_xy_xyy_xy_0[i] * pa_x[i] - ta2_xy_xyy_xy_1[i] * pc_x[i];

        ta2_xy_xxyy_xz_0[i] =
            ta2_xy_xx_xz_0[i] * fe_0 - ta2_xy_xx_xz_1[i] * fe_0 + ta1_x_xxy_xz_1[i] + ta2_xy_xxy_xz_0[i] * pa_y[i] - ta2_xy_xxy_xz_1[i] * pc_y[i];

        ta2_xy_xxyy_yy_0[i] =
            ta2_xy_yy_yy_0[i] * fe_0 - ta2_xy_yy_yy_1[i] * fe_0 + ta1_y_xyy_yy_1[i] + ta2_xy_xyy_yy_0[i] * pa_x[i] - ta2_xy_xyy_yy_1[i] * pc_x[i];

        ta2_xy_xxyy_yz_0[i] =
            ta2_xy_yy_yz_0[i] * fe_0 - ta2_xy_yy_yz_1[i] * fe_0 + ta1_y_xyy_yz_1[i] + ta2_xy_xyy_yz_0[i] * pa_x[i] - ta2_xy_xyy_yz_1[i] * pc_x[i];

        ta2_xy_xxyy_zz_0[i] =
            ta2_xy_yy_zz_0[i] * fe_0 - ta2_xy_yy_zz_1[i] * fe_0 + ta1_y_xyy_zz_1[i] + ta2_xy_xyy_zz_0[i] * pa_x[i] - ta2_xy_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 114-120 components of targeted buffer : GD

    auto ta2_xy_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 114);

    auto ta2_xy_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 115);

    auto ta2_xy_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 116);

    auto ta2_xy_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 117);

    auto ta2_xy_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 118);

    auto ta2_xy_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 119);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxz_xz_1,   \
                             ta1_x_xxz_zz_1,   \
                             ta2_xy_xxy_xx_0,  \
                             ta2_xy_xxy_xx_1,  \
                             ta2_xy_xxy_xy_0,  \
                             ta2_xy_xxy_xy_1,  \
                             ta2_xy_xxy_y_0,   \
                             ta2_xy_xxy_y_1,   \
                             ta2_xy_xxy_yy_0,  \
                             ta2_xy_xxy_yy_1,  \
                             ta2_xy_xxy_yz_0,  \
                             ta2_xy_xxy_yz_1,  \
                             ta2_xy_xxyz_xx_0, \
                             ta2_xy_xxyz_xy_0, \
                             ta2_xy_xxyz_xz_0, \
                             ta2_xy_xxyz_yy_0, \
                             ta2_xy_xxyz_yz_0, \
                             ta2_xy_xxyz_zz_0, \
                             ta2_xy_xxz_xz_0,  \
                             ta2_xy_xxz_xz_1,  \
                             ta2_xy_xxz_zz_0,  \
                             ta2_xy_xxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyz_xx_0[i] = ta2_xy_xxy_xx_0[i] * pa_z[i] - ta2_xy_xxy_xx_1[i] * pc_z[i];

        ta2_xy_xxyz_xy_0[i] = ta2_xy_xxy_xy_0[i] * pa_z[i] - ta2_xy_xxy_xy_1[i] * pc_z[i];

        ta2_xy_xxyz_xz_0[i] = ta1_x_xxz_xz_1[i] + ta2_xy_xxz_xz_0[i] * pa_y[i] - ta2_xy_xxz_xz_1[i] * pc_y[i];

        ta2_xy_xxyz_yy_0[i] = ta2_xy_xxy_yy_0[i] * pa_z[i] - ta2_xy_xxy_yy_1[i] * pc_z[i];

        ta2_xy_xxyz_yz_0[i] = ta2_xy_xxy_y_0[i] * fe_0 - ta2_xy_xxy_y_1[i] * fe_0 + ta2_xy_xxy_yz_0[i] * pa_z[i] - ta2_xy_xxy_yz_1[i] * pc_z[i];

        ta2_xy_xxyz_zz_0[i] = ta1_x_xxz_zz_1[i] + ta2_xy_xxz_zz_0[i] * pa_y[i] - ta2_xy_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 120-126 components of targeted buffer : GD

    auto ta2_xy_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 120);

    auto ta2_xy_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 121);

    auto ta2_xy_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 122);

    auto ta2_xy_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 123);

    auto ta2_xy_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 124);

    auto ta2_xy_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 125);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xzz_yz_1,   \
                             ta1_y_xzz_zz_1,   \
                             ta2_xy_xx_xx_0,   \
                             ta2_xy_xx_xx_1,   \
                             ta2_xy_xx_xy_0,   \
                             ta2_xy_xx_xy_1,   \
                             ta2_xy_xx_xz_0,   \
                             ta2_xy_xx_xz_1,   \
                             ta2_xy_xx_yy_0,   \
                             ta2_xy_xx_yy_1,   \
                             ta2_xy_xxz_x_0,   \
                             ta2_xy_xxz_x_1,   \
                             ta2_xy_xxz_xx_0,  \
                             ta2_xy_xxz_xx_1,  \
                             ta2_xy_xxz_xy_0,  \
                             ta2_xy_xxz_xy_1,  \
                             ta2_xy_xxz_xz_0,  \
                             ta2_xy_xxz_xz_1,  \
                             ta2_xy_xxz_yy_0,  \
                             ta2_xy_xxz_yy_1,  \
                             ta2_xy_xxzz_xx_0, \
                             ta2_xy_xxzz_xy_0, \
                             ta2_xy_xxzz_xz_0, \
                             ta2_xy_xxzz_yy_0, \
                             ta2_xy_xxzz_yz_0, \
                             ta2_xy_xxzz_zz_0, \
                             ta2_xy_xzz_yz_0,  \
                             ta2_xy_xzz_yz_1,  \
                             ta2_xy_xzz_zz_0,  \
                             ta2_xy_xzz_zz_1,  \
                             ta2_xy_zz_yz_0,   \
                             ta2_xy_zz_yz_1,   \
                             ta2_xy_zz_zz_0,   \
                             ta2_xy_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxzz_xx_0[i] = ta2_xy_xx_xx_0[i] * fe_0 - ta2_xy_xx_xx_1[i] * fe_0 + ta2_xy_xxz_xx_0[i] * pa_z[i] - ta2_xy_xxz_xx_1[i] * pc_z[i];

        ta2_xy_xxzz_xy_0[i] = ta2_xy_xx_xy_0[i] * fe_0 - ta2_xy_xx_xy_1[i] * fe_0 + ta2_xy_xxz_xy_0[i] * pa_z[i] - ta2_xy_xxz_xy_1[i] * pc_z[i];

        ta2_xy_xxzz_xz_0[i] = ta2_xy_xx_xz_0[i] * fe_0 - ta2_xy_xx_xz_1[i] * fe_0 + ta2_xy_xxz_x_0[i] * fe_0 - ta2_xy_xxz_x_1[i] * fe_0 +
                              ta2_xy_xxz_xz_0[i] * pa_z[i] - ta2_xy_xxz_xz_1[i] * pc_z[i];

        ta2_xy_xxzz_yy_0[i] = ta2_xy_xx_yy_0[i] * fe_0 - ta2_xy_xx_yy_1[i] * fe_0 + ta2_xy_xxz_yy_0[i] * pa_z[i] - ta2_xy_xxz_yy_1[i] * pc_z[i];

        ta2_xy_xxzz_yz_0[i] =
            ta2_xy_zz_yz_0[i] * fe_0 - ta2_xy_zz_yz_1[i] * fe_0 + ta1_y_xzz_yz_1[i] + ta2_xy_xzz_yz_0[i] * pa_x[i] - ta2_xy_xzz_yz_1[i] * pc_x[i];

        ta2_xy_xxzz_zz_0[i] =
            ta2_xy_zz_zz_0[i] * fe_0 - ta2_xy_zz_zz_1[i] * fe_0 + ta1_y_xzz_zz_1[i] + ta2_xy_xzz_zz_0[i] * pa_x[i] - ta2_xy_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : GD

    auto ta2_xy_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 126);

    auto ta2_xy_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 127);

    auto ta2_xy_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 128);

    auto ta2_xy_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 129);

    auto ta2_xy_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 130);

    auto ta2_xy_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 131);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_zz_1,   \
                             ta2_xy_xyyy_xx_0, \
                             ta2_xy_xyyy_xy_0, \
                             ta2_xy_xyyy_xz_0, \
                             ta2_xy_xyyy_yy_0, \
                             ta2_xy_xyyy_yz_0, \
                             ta2_xy_xyyy_zz_0, \
                             ta2_xy_yyy_x_0,   \
                             ta2_xy_yyy_x_1,   \
                             ta2_xy_yyy_xx_0,  \
                             ta2_xy_yyy_xx_1,  \
                             ta2_xy_yyy_xy_0,  \
                             ta2_xy_yyy_xy_1,  \
                             ta2_xy_yyy_xz_0,  \
                             ta2_xy_yyy_xz_1,  \
                             ta2_xy_yyy_y_0,   \
                             ta2_xy_yyy_y_1,   \
                             ta2_xy_yyy_yy_0,  \
                             ta2_xy_yyy_yy_1,  \
                             ta2_xy_yyy_yz_0,  \
                             ta2_xy_yyy_yz_1,  \
                             ta2_xy_yyy_z_0,   \
                             ta2_xy_yyy_z_1,   \
                             ta2_xy_yyy_zz_0,  \
                             ta2_xy_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyy_xx_0[i] = 2.0 * ta2_xy_yyy_x_0[i] * fe_0 - 2.0 * ta2_xy_yyy_x_1[i] * fe_0 + ta1_y_yyy_xx_1[i] + ta2_xy_yyy_xx_0[i] * pa_x[i] -
                              ta2_xy_yyy_xx_1[i] * pc_x[i];

        ta2_xy_xyyy_xy_0[i] =
            ta2_xy_yyy_y_0[i] * fe_0 - ta2_xy_yyy_y_1[i] * fe_0 + ta1_y_yyy_xy_1[i] + ta2_xy_yyy_xy_0[i] * pa_x[i] - ta2_xy_yyy_xy_1[i] * pc_x[i];

        ta2_xy_xyyy_xz_0[i] =
            ta2_xy_yyy_z_0[i] * fe_0 - ta2_xy_yyy_z_1[i] * fe_0 + ta1_y_yyy_xz_1[i] + ta2_xy_yyy_xz_0[i] * pa_x[i] - ta2_xy_yyy_xz_1[i] * pc_x[i];

        ta2_xy_xyyy_yy_0[i] = ta1_y_yyy_yy_1[i] + ta2_xy_yyy_yy_0[i] * pa_x[i] - ta2_xy_yyy_yy_1[i] * pc_x[i];

        ta2_xy_xyyy_yz_0[i] = ta1_y_yyy_yz_1[i] + ta2_xy_yyy_yz_0[i] * pa_x[i] - ta2_xy_yyy_yz_1[i] * pc_x[i];

        ta2_xy_xyyy_zz_0[i] = ta1_y_yyy_zz_1[i] + ta2_xy_yyy_zz_0[i] * pa_x[i] - ta2_xy_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 132-138 components of targeted buffer : GD

    auto ta2_xy_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 132);

    auto ta2_xy_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 133);

    auto ta2_xy_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 134);

    auto ta2_xy_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 135);

    auto ta2_xy_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 136);

    auto ta2_xy_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 137);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_yyz_yz_1,   \
                             ta1_y_yyz_zz_1,   \
                             ta2_xy_xyy_x_0,   \
                             ta2_xy_xyy_x_1,   \
                             ta2_xy_xyy_xx_0,  \
                             ta2_xy_xyy_xx_1,  \
                             ta2_xy_xyy_xy_0,  \
                             ta2_xy_xyy_xy_1,  \
                             ta2_xy_xyy_xz_0,  \
                             ta2_xy_xyy_xz_1,  \
                             ta2_xy_xyy_yy_0,  \
                             ta2_xy_xyy_yy_1,  \
                             ta2_xy_xyyz_xx_0, \
                             ta2_xy_xyyz_xy_0, \
                             ta2_xy_xyyz_xz_0, \
                             ta2_xy_xyyz_yy_0, \
                             ta2_xy_xyyz_yz_0, \
                             ta2_xy_xyyz_zz_0, \
                             ta2_xy_yyz_yz_0,  \
                             ta2_xy_yyz_yz_1,  \
                             ta2_xy_yyz_zz_0,  \
                             ta2_xy_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyz_xx_0[i] = ta2_xy_xyy_xx_0[i] * pa_z[i] - ta2_xy_xyy_xx_1[i] * pc_z[i];

        ta2_xy_xyyz_xy_0[i] = ta2_xy_xyy_xy_0[i] * pa_z[i] - ta2_xy_xyy_xy_1[i] * pc_z[i];

        ta2_xy_xyyz_xz_0[i] = ta2_xy_xyy_x_0[i] * fe_0 - ta2_xy_xyy_x_1[i] * fe_0 + ta2_xy_xyy_xz_0[i] * pa_z[i] - ta2_xy_xyy_xz_1[i] * pc_z[i];

        ta2_xy_xyyz_yy_0[i] = ta2_xy_xyy_yy_0[i] * pa_z[i] - ta2_xy_xyy_yy_1[i] * pc_z[i];

        ta2_xy_xyyz_yz_0[i] = ta1_y_yyz_yz_1[i] + ta2_xy_yyz_yz_0[i] * pa_x[i] - ta2_xy_yyz_yz_1[i] * pc_x[i];

        ta2_xy_xyyz_zz_0[i] = ta1_y_yyz_zz_1[i] + ta2_xy_yyz_zz_0[i] * pa_x[i] - ta2_xy_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 138-144 components of targeted buffer : GD

    auto ta2_xy_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 138);

    auto ta2_xy_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 139);

    auto ta2_xy_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 140);

    auto ta2_xy_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 141);

    auto ta2_xy_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 142);

    auto ta2_xy_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 143);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xzz_xx_1,   \
                             ta1_x_xzz_xz_1,   \
                             ta1_y_yzz_yy_1,   \
                             ta1_y_yzz_yz_1,   \
                             ta1_y_yzz_zz_1,   \
                             ta2_xy_xy_xy_0,   \
                             ta2_xy_xy_xy_1,   \
                             ta2_xy_xyz_xy_0,  \
                             ta2_xy_xyz_xy_1,  \
                             ta2_xy_xyzz_xx_0, \
                             ta2_xy_xyzz_xy_0, \
                             ta2_xy_xyzz_xz_0, \
                             ta2_xy_xyzz_yy_0, \
                             ta2_xy_xyzz_yz_0, \
                             ta2_xy_xyzz_zz_0, \
                             ta2_xy_xzz_xx_0,  \
                             ta2_xy_xzz_xx_1,  \
                             ta2_xy_xzz_xz_0,  \
                             ta2_xy_xzz_xz_1,  \
                             ta2_xy_yzz_yy_0,  \
                             ta2_xy_yzz_yy_1,  \
                             ta2_xy_yzz_yz_0,  \
                             ta2_xy_yzz_yz_1,  \
                             ta2_xy_yzz_zz_0,  \
                             ta2_xy_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyzz_xx_0[i] = ta1_x_xzz_xx_1[i] + ta2_xy_xzz_xx_0[i] * pa_y[i] - ta2_xy_xzz_xx_1[i] * pc_y[i];

        ta2_xy_xyzz_xy_0[i] = ta2_xy_xy_xy_0[i] * fe_0 - ta2_xy_xy_xy_1[i] * fe_0 + ta2_xy_xyz_xy_0[i] * pa_z[i] - ta2_xy_xyz_xy_1[i] * pc_z[i];

        ta2_xy_xyzz_xz_0[i] = ta1_x_xzz_xz_1[i] + ta2_xy_xzz_xz_0[i] * pa_y[i] - ta2_xy_xzz_xz_1[i] * pc_y[i];

        ta2_xy_xyzz_yy_0[i] = ta1_y_yzz_yy_1[i] + ta2_xy_yzz_yy_0[i] * pa_x[i] - ta2_xy_yzz_yy_1[i] * pc_x[i];

        ta2_xy_xyzz_yz_0[i] = ta1_y_yzz_yz_1[i] + ta2_xy_yzz_yz_0[i] * pa_x[i] - ta2_xy_yzz_yz_1[i] * pc_x[i];

        ta2_xy_xyzz_zz_0[i] = ta1_y_yzz_zz_1[i] + ta2_xy_yzz_zz_0[i] * pa_x[i] - ta2_xy_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 144-150 components of targeted buffer : GD

    auto ta2_xy_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 144);

    auto ta2_xy_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 145);

    auto ta2_xy_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 146);

    auto ta2_xy_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 147);

    auto ta2_xy_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 148);

    auto ta2_xy_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 149);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_zz_1,   \
                             ta2_xy_xz_xx_0,   \
                             ta2_xy_xz_xx_1,   \
                             ta2_xy_xz_xy_0,   \
                             ta2_xy_xz_xy_1,   \
                             ta2_xy_xzz_xx_0,  \
                             ta2_xy_xzz_xx_1,  \
                             ta2_xy_xzz_xy_0,  \
                             ta2_xy_xzz_xy_1,  \
                             ta2_xy_xzzz_xx_0, \
                             ta2_xy_xzzz_xy_0, \
                             ta2_xy_xzzz_xz_0, \
                             ta2_xy_xzzz_yy_0, \
                             ta2_xy_xzzz_yz_0, \
                             ta2_xy_xzzz_zz_0, \
                             ta2_xy_zzz_xz_0,  \
                             ta2_xy_zzz_xz_1,  \
                             ta2_xy_zzz_yy_0,  \
                             ta2_xy_zzz_yy_1,  \
                             ta2_xy_zzz_yz_0,  \
                             ta2_xy_zzz_yz_1,  \
                             ta2_xy_zzz_z_0,   \
                             ta2_xy_zzz_z_1,   \
                             ta2_xy_zzz_zz_0,  \
                             ta2_xy_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzzz_xx_0[i] =
            2.0 * ta2_xy_xz_xx_0[i] * fe_0 - 2.0 * ta2_xy_xz_xx_1[i] * fe_0 + ta2_xy_xzz_xx_0[i] * pa_z[i] - ta2_xy_xzz_xx_1[i] * pc_z[i];

        ta2_xy_xzzz_xy_0[i] =
            2.0 * ta2_xy_xz_xy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xy_1[i] * fe_0 + ta2_xy_xzz_xy_0[i] * pa_z[i] - ta2_xy_xzz_xy_1[i] * pc_z[i];

        ta2_xy_xzzz_xz_0[i] =
            ta2_xy_zzz_z_0[i] * fe_0 - ta2_xy_zzz_z_1[i] * fe_0 + ta1_y_zzz_xz_1[i] + ta2_xy_zzz_xz_0[i] * pa_x[i] - ta2_xy_zzz_xz_1[i] * pc_x[i];

        ta2_xy_xzzz_yy_0[i] = ta1_y_zzz_yy_1[i] + ta2_xy_zzz_yy_0[i] * pa_x[i] - ta2_xy_zzz_yy_1[i] * pc_x[i];

        ta2_xy_xzzz_yz_0[i] = ta1_y_zzz_yz_1[i] + ta2_xy_zzz_yz_0[i] * pa_x[i] - ta2_xy_zzz_yz_1[i] * pc_x[i];

        ta2_xy_xzzz_zz_0[i] = ta1_y_zzz_zz_1[i] + ta2_xy_zzz_zz_0[i] * pa_x[i] - ta2_xy_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 150-156 components of targeted buffer : GD

    auto ta2_xy_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 150);

    auto ta2_xy_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 151);

    auto ta2_xy_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 152);

    auto ta2_xy_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 153);

    auto ta2_xy_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 154);

    auto ta2_xy_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 155);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yyy_xx_1,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_xz_1,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_zz_1,   \
                             ta2_xy_yy_xx_0,   \
                             ta2_xy_yy_xx_1,   \
                             ta2_xy_yy_xy_0,   \
                             ta2_xy_yy_xy_1,   \
                             ta2_xy_yy_xz_0,   \
                             ta2_xy_yy_xz_1,   \
                             ta2_xy_yy_yy_0,   \
                             ta2_xy_yy_yy_1,   \
                             ta2_xy_yy_yz_0,   \
                             ta2_xy_yy_yz_1,   \
                             ta2_xy_yy_zz_0,   \
                             ta2_xy_yy_zz_1,   \
                             ta2_xy_yyy_x_0,   \
                             ta2_xy_yyy_x_1,   \
                             ta2_xy_yyy_xx_0,  \
                             ta2_xy_yyy_xx_1,  \
                             ta2_xy_yyy_xy_0,  \
                             ta2_xy_yyy_xy_1,  \
                             ta2_xy_yyy_xz_0,  \
                             ta2_xy_yyy_xz_1,  \
                             ta2_xy_yyy_y_0,   \
                             ta2_xy_yyy_y_1,   \
                             ta2_xy_yyy_yy_0,  \
                             ta2_xy_yyy_yy_1,  \
                             ta2_xy_yyy_yz_0,  \
                             ta2_xy_yyy_yz_1,  \
                             ta2_xy_yyy_z_0,   \
                             ta2_xy_yyy_z_1,   \
                             ta2_xy_yyy_zz_0,  \
                             ta2_xy_yyy_zz_1,  \
                             ta2_xy_yyyy_xx_0, \
                             ta2_xy_yyyy_xy_0, \
                             ta2_xy_yyyy_xz_0, \
                             ta2_xy_yyyy_yy_0, \
                             ta2_xy_yyyy_yz_0, \
                             ta2_xy_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyy_xx_0[i] = 3.0 * ta2_xy_yy_xx_0[i] * fe_0 - 3.0 * ta2_xy_yy_xx_1[i] * fe_0 + ta1_x_yyy_xx_1[i] + ta2_xy_yyy_xx_0[i] * pa_y[i] -
                              ta2_xy_yyy_xx_1[i] * pc_y[i];

        ta2_xy_yyyy_xy_0[i] = 3.0 * ta2_xy_yy_xy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xy_1[i] * fe_0 + ta2_xy_yyy_x_0[i] * fe_0 - ta2_xy_yyy_x_1[i] * fe_0 +
                              ta1_x_yyy_xy_1[i] + ta2_xy_yyy_xy_0[i] * pa_y[i] - ta2_xy_yyy_xy_1[i] * pc_y[i];

        ta2_xy_yyyy_xz_0[i] = 3.0 * ta2_xy_yy_xz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xz_1[i] * fe_0 + ta1_x_yyy_xz_1[i] + ta2_xy_yyy_xz_0[i] * pa_y[i] -
                              ta2_xy_yyy_xz_1[i] * pc_y[i];

        ta2_xy_yyyy_yy_0[i] = 3.0 * ta2_xy_yy_yy_0[i] * fe_0 - 3.0 * ta2_xy_yy_yy_1[i] * fe_0 + 2.0 * ta2_xy_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_xy_yyy_y_1[i] * fe_0 + ta1_x_yyy_yy_1[i] + ta2_xy_yyy_yy_0[i] * pa_y[i] - ta2_xy_yyy_yy_1[i] * pc_y[i];

        ta2_xy_yyyy_yz_0[i] = 3.0 * ta2_xy_yy_yz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yz_1[i] * fe_0 + ta2_xy_yyy_z_0[i] * fe_0 - ta2_xy_yyy_z_1[i] * fe_0 +
                              ta1_x_yyy_yz_1[i] + ta2_xy_yyy_yz_0[i] * pa_y[i] - ta2_xy_yyy_yz_1[i] * pc_y[i];

        ta2_xy_yyyy_zz_0[i] = 3.0 * ta2_xy_yy_zz_0[i] * fe_0 - 3.0 * ta2_xy_yy_zz_1[i] * fe_0 + ta1_x_yyy_zz_1[i] + ta2_xy_yyy_zz_0[i] * pa_y[i] -
                              ta2_xy_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 156-162 components of targeted buffer : GD

    auto ta2_xy_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 156);

    auto ta2_xy_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 157);

    auto ta2_xy_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 158);

    auto ta2_xy_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 159);

    auto ta2_xy_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 160);

    auto ta2_xy_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 161);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xy_yyy_x_0,   \
                             ta2_xy_yyy_x_1,   \
                             ta2_xy_yyy_xx_0,  \
                             ta2_xy_yyy_xx_1,  \
                             ta2_xy_yyy_xy_0,  \
                             ta2_xy_yyy_xy_1,  \
                             ta2_xy_yyy_xz_0,  \
                             ta2_xy_yyy_xz_1,  \
                             ta2_xy_yyy_y_0,   \
                             ta2_xy_yyy_y_1,   \
                             ta2_xy_yyy_yy_0,  \
                             ta2_xy_yyy_yy_1,  \
                             ta2_xy_yyy_yz_0,  \
                             ta2_xy_yyy_yz_1,  \
                             ta2_xy_yyy_z_0,   \
                             ta2_xy_yyy_z_1,   \
                             ta2_xy_yyy_zz_0,  \
                             ta2_xy_yyy_zz_1,  \
                             ta2_xy_yyyz_xx_0, \
                             ta2_xy_yyyz_xy_0, \
                             ta2_xy_yyyz_xz_0, \
                             ta2_xy_yyyz_yy_0, \
                             ta2_xy_yyyz_yz_0, \
                             ta2_xy_yyyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyz_xx_0[i] = ta2_xy_yyy_xx_0[i] * pa_z[i] - ta2_xy_yyy_xx_1[i] * pc_z[i];

        ta2_xy_yyyz_xy_0[i] = ta2_xy_yyy_xy_0[i] * pa_z[i] - ta2_xy_yyy_xy_1[i] * pc_z[i];

        ta2_xy_yyyz_xz_0[i] = ta2_xy_yyy_x_0[i] * fe_0 - ta2_xy_yyy_x_1[i] * fe_0 + ta2_xy_yyy_xz_0[i] * pa_z[i] - ta2_xy_yyy_xz_1[i] * pc_z[i];

        ta2_xy_yyyz_yy_0[i] = ta2_xy_yyy_yy_0[i] * pa_z[i] - ta2_xy_yyy_yy_1[i] * pc_z[i];

        ta2_xy_yyyz_yz_0[i] = ta2_xy_yyy_y_0[i] * fe_0 - ta2_xy_yyy_y_1[i] * fe_0 + ta2_xy_yyy_yz_0[i] * pa_z[i] - ta2_xy_yyy_yz_1[i] * pc_z[i];

        ta2_xy_yyyz_zz_0[i] =
            2.0 * ta2_xy_yyy_z_0[i] * fe_0 - 2.0 * ta2_xy_yyy_z_1[i] * fe_0 + ta2_xy_yyy_zz_0[i] * pa_z[i] - ta2_xy_yyy_zz_1[i] * pc_z[i];
    }

    // Set up 162-168 components of targeted buffer : GD

    auto ta2_xy_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 162);

    auto ta2_xy_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 163);

    auto ta2_xy_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 164);

    auto ta2_xy_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 165);

    auto ta2_xy_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 166);

    auto ta2_xy_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 167);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yzz_xz_1,   \
                             ta1_x_yzz_zz_1,   \
                             ta2_xy_yy_xx_0,   \
                             ta2_xy_yy_xx_1,   \
                             ta2_xy_yy_xy_0,   \
                             ta2_xy_yy_xy_1,   \
                             ta2_xy_yy_yy_0,   \
                             ta2_xy_yy_yy_1,   \
                             ta2_xy_yy_yz_0,   \
                             ta2_xy_yy_yz_1,   \
                             ta2_xy_yyz_xx_0,  \
                             ta2_xy_yyz_xx_1,  \
                             ta2_xy_yyz_xy_0,  \
                             ta2_xy_yyz_xy_1,  \
                             ta2_xy_yyz_y_0,   \
                             ta2_xy_yyz_y_1,   \
                             ta2_xy_yyz_yy_0,  \
                             ta2_xy_yyz_yy_1,  \
                             ta2_xy_yyz_yz_0,  \
                             ta2_xy_yyz_yz_1,  \
                             ta2_xy_yyzz_xx_0, \
                             ta2_xy_yyzz_xy_0, \
                             ta2_xy_yyzz_xz_0, \
                             ta2_xy_yyzz_yy_0, \
                             ta2_xy_yyzz_yz_0, \
                             ta2_xy_yyzz_zz_0, \
                             ta2_xy_yzz_xz_0,  \
                             ta2_xy_yzz_xz_1,  \
                             ta2_xy_yzz_zz_0,  \
                             ta2_xy_yzz_zz_1,  \
                             ta2_xy_zz_xz_0,   \
                             ta2_xy_zz_xz_1,   \
                             ta2_xy_zz_zz_0,   \
                             ta2_xy_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyzz_xx_0[i] = ta2_xy_yy_xx_0[i] * fe_0 - ta2_xy_yy_xx_1[i] * fe_0 + ta2_xy_yyz_xx_0[i] * pa_z[i] - ta2_xy_yyz_xx_1[i] * pc_z[i];

        ta2_xy_yyzz_xy_0[i] = ta2_xy_yy_xy_0[i] * fe_0 - ta2_xy_yy_xy_1[i] * fe_0 + ta2_xy_yyz_xy_0[i] * pa_z[i] - ta2_xy_yyz_xy_1[i] * pc_z[i];

        ta2_xy_yyzz_xz_0[i] =
            ta2_xy_zz_xz_0[i] * fe_0 - ta2_xy_zz_xz_1[i] * fe_0 + ta1_x_yzz_xz_1[i] + ta2_xy_yzz_xz_0[i] * pa_y[i] - ta2_xy_yzz_xz_1[i] * pc_y[i];

        ta2_xy_yyzz_yy_0[i] = ta2_xy_yy_yy_0[i] * fe_0 - ta2_xy_yy_yy_1[i] * fe_0 + ta2_xy_yyz_yy_0[i] * pa_z[i] - ta2_xy_yyz_yy_1[i] * pc_z[i];

        ta2_xy_yyzz_yz_0[i] = ta2_xy_yy_yz_0[i] * fe_0 - ta2_xy_yy_yz_1[i] * fe_0 + ta2_xy_yyz_y_0[i] * fe_0 - ta2_xy_yyz_y_1[i] * fe_0 +
                              ta2_xy_yyz_yz_0[i] * pa_z[i] - ta2_xy_yyz_yz_1[i] * pc_z[i];

        ta2_xy_yyzz_zz_0[i] =
            ta2_xy_zz_zz_0[i] * fe_0 - ta2_xy_zz_zz_1[i] * fe_0 + ta1_x_yzz_zz_1[i] + ta2_xy_yzz_zz_0[i] * pa_y[i] - ta2_xy_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 168-174 components of targeted buffer : GD

    auto ta2_xy_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 168);

    auto ta2_xy_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 169);

    auto ta2_xy_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 170);

    auto ta2_xy_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 171);

    auto ta2_xy_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 172);

    auto ta2_xy_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 173);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_1,   \
                             ta2_xy_yz_xy_0,   \
                             ta2_xy_yz_xy_1,   \
                             ta2_xy_yz_yy_0,   \
                             ta2_xy_yz_yy_1,   \
                             ta2_xy_yzz_xy_0,  \
                             ta2_xy_yzz_xy_1,  \
                             ta2_xy_yzz_yy_0,  \
                             ta2_xy_yzz_yy_1,  \
                             ta2_xy_yzzz_xx_0, \
                             ta2_xy_yzzz_xy_0, \
                             ta2_xy_yzzz_xz_0, \
                             ta2_xy_yzzz_yy_0, \
                             ta2_xy_yzzz_yz_0, \
                             ta2_xy_yzzz_zz_0, \
                             ta2_xy_zzz_xx_0,  \
                             ta2_xy_zzz_xx_1,  \
                             ta2_xy_zzz_xz_0,  \
                             ta2_xy_zzz_xz_1,  \
                             ta2_xy_zzz_yz_0,  \
                             ta2_xy_zzz_yz_1,  \
                             ta2_xy_zzz_z_0,   \
                             ta2_xy_zzz_z_1,   \
                             ta2_xy_zzz_zz_0,  \
                             ta2_xy_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzzz_xx_0[i] = ta1_x_zzz_xx_1[i] + ta2_xy_zzz_xx_0[i] * pa_y[i] - ta2_xy_zzz_xx_1[i] * pc_y[i];

        ta2_xy_yzzz_xy_0[i] =
            2.0 * ta2_xy_yz_xy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xy_1[i] * fe_0 + ta2_xy_yzz_xy_0[i] * pa_z[i] - ta2_xy_yzz_xy_1[i] * pc_z[i];

        ta2_xy_yzzz_xz_0[i] = ta1_x_zzz_xz_1[i] + ta2_xy_zzz_xz_0[i] * pa_y[i] - ta2_xy_zzz_xz_1[i] * pc_y[i];

        ta2_xy_yzzz_yy_0[i] =
            2.0 * ta2_xy_yz_yy_0[i] * fe_0 - 2.0 * ta2_xy_yz_yy_1[i] * fe_0 + ta2_xy_yzz_yy_0[i] * pa_z[i] - ta2_xy_yzz_yy_1[i] * pc_z[i];

        ta2_xy_yzzz_yz_0[i] =
            ta2_xy_zzz_z_0[i] * fe_0 - ta2_xy_zzz_z_1[i] * fe_0 + ta1_x_zzz_yz_1[i] + ta2_xy_zzz_yz_0[i] * pa_y[i] - ta2_xy_zzz_yz_1[i] * pc_y[i];

        ta2_xy_yzzz_zz_0[i] = ta1_x_zzz_zz_1[i] + ta2_xy_zzz_zz_0[i] * pa_y[i] - ta2_xy_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 174-180 components of targeted buffer : GD

    auto ta2_xy_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 174);

    auto ta2_xy_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 175);

    auto ta2_xy_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 176);

    auto ta2_xy_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 177);

    auto ta2_xy_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 178);

    auto ta2_xy_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 179);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xy_zz_xx_0,   \
                             ta2_xy_zz_xx_1,   \
                             ta2_xy_zz_xy_0,   \
                             ta2_xy_zz_xy_1,   \
                             ta2_xy_zz_xz_0,   \
                             ta2_xy_zz_xz_1,   \
                             ta2_xy_zz_yy_0,   \
                             ta2_xy_zz_yy_1,   \
                             ta2_xy_zz_yz_0,   \
                             ta2_xy_zz_yz_1,   \
                             ta2_xy_zz_zz_0,   \
                             ta2_xy_zz_zz_1,   \
                             ta2_xy_zzz_x_0,   \
                             ta2_xy_zzz_x_1,   \
                             ta2_xy_zzz_xx_0,  \
                             ta2_xy_zzz_xx_1,  \
                             ta2_xy_zzz_xy_0,  \
                             ta2_xy_zzz_xy_1,  \
                             ta2_xy_zzz_xz_0,  \
                             ta2_xy_zzz_xz_1,  \
                             ta2_xy_zzz_y_0,   \
                             ta2_xy_zzz_y_1,   \
                             ta2_xy_zzz_yy_0,  \
                             ta2_xy_zzz_yy_1,  \
                             ta2_xy_zzz_yz_0,  \
                             ta2_xy_zzz_yz_1,  \
                             ta2_xy_zzz_z_0,   \
                             ta2_xy_zzz_z_1,   \
                             ta2_xy_zzz_zz_0,  \
                             ta2_xy_zzz_zz_1,  \
                             ta2_xy_zzzz_xx_0, \
                             ta2_xy_zzzz_xy_0, \
                             ta2_xy_zzzz_xz_0, \
                             ta2_xy_zzzz_yy_0, \
                             ta2_xy_zzzz_yz_0, \
                             ta2_xy_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzzz_xx_0[i] =
            3.0 * ta2_xy_zz_xx_0[i] * fe_0 - 3.0 * ta2_xy_zz_xx_1[i] * fe_0 + ta2_xy_zzz_xx_0[i] * pa_z[i] - ta2_xy_zzz_xx_1[i] * pc_z[i];

        ta2_xy_zzzz_xy_0[i] =
            3.0 * ta2_xy_zz_xy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xy_1[i] * fe_0 + ta2_xy_zzz_xy_0[i] * pa_z[i] - ta2_xy_zzz_xy_1[i] * pc_z[i];

        ta2_xy_zzzz_xz_0[i] = 3.0 * ta2_xy_zz_xz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xz_1[i] * fe_0 + ta2_xy_zzz_x_0[i] * fe_0 - ta2_xy_zzz_x_1[i] * fe_0 +
                              ta2_xy_zzz_xz_0[i] * pa_z[i] - ta2_xy_zzz_xz_1[i] * pc_z[i];

        ta2_xy_zzzz_yy_0[i] =
            3.0 * ta2_xy_zz_yy_0[i] * fe_0 - 3.0 * ta2_xy_zz_yy_1[i] * fe_0 + ta2_xy_zzz_yy_0[i] * pa_z[i] - ta2_xy_zzz_yy_1[i] * pc_z[i];

        ta2_xy_zzzz_yz_0[i] = 3.0 * ta2_xy_zz_yz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yz_1[i] * fe_0 + ta2_xy_zzz_y_0[i] * fe_0 - ta2_xy_zzz_y_1[i] * fe_0 +
                              ta2_xy_zzz_yz_0[i] * pa_z[i] - ta2_xy_zzz_yz_1[i] * pc_z[i];

        ta2_xy_zzzz_zz_0[i] = 3.0 * ta2_xy_zz_zz_0[i] * fe_0 - 3.0 * ta2_xy_zz_zz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_xy_zzz_z_1[i] * fe_0 + ta2_xy_zzz_zz_0[i] * pa_z[i] - ta2_xy_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 180-186 components of targeted buffer : GD

    auto ta2_xz_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 180);

    auto ta2_xz_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 181);

    auto ta2_xz_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 182);

    auto ta2_xz_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 183);

    auto ta2_xz_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 184);

    auto ta2_xz_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 185);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_yy_1,   \
                             ta1_z_xxx_yz_1,   \
                             ta1_z_xxx_zz_1,   \
                             ta2_xz_xx_xx_0,   \
                             ta2_xz_xx_xx_1,   \
                             ta2_xz_xx_xy_0,   \
                             ta2_xz_xx_xy_1,   \
                             ta2_xz_xx_xz_0,   \
                             ta2_xz_xx_xz_1,   \
                             ta2_xz_xx_yy_0,   \
                             ta2_xz_xx_yy_1,   \
                             ta2_xz_xx_yz_0,   \
                             ta2_xz_xx_yz_1,   \
                             ta2_xz_xx_zz_0,   \
                             ta2_xz_xx_zz_1,   \
                             ta2_xz_xxx_x_0,   \
                             ta2_xz_xxx_x_1,   \
                             ta2_xz_xxx_xx_0,  \
                             ta2_xz_xxx_xx_1,  \
                             ta2_xz_xxx_xy_0,  \
                             ta2_xz_xxx_xy_1,  \
                             ta2_xz_xxx_xz_0,  \
                             ta2_xz_xxx_xz_1,  \
                             ta2_xz_xxx_y_0,   \
                             ta2_xz_xxx_y_1,   \
                             ta2_xz_xxx_yy_0,  \
                             ta2_xz_xxx_yy_1,  \
                             ta2_xz_xxx_yz_0,  \
                             ta2_xz_xxx_yz_1,  \
                             ta2_xz_xxx_z_0,   \
                             ta2_xz_xxx_z_1,   \
                             ta2_xz_xxx_zz_0,  \
                             ta2_xz_xxx_zz_1,  \
                             ta2_xz_xxxx_xx_0, \
                             ta2_xz_xxxx_xy_0, \
                             ta2_xz_xxxx_xz_0, \
                             ta2_xz_xxxx_yy_0, \
                             ta2_xz_xxxx_yz_0, \
                             ta2_xz_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxx_xx_0[i] = 3.0 * ta2_xz_xx_xx_0[i] * fe_0 - 3.0 * ta2_xz_xx_xx_1[i] * fe_0 + 2.0 * ta2_xz_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_xz_xxx_x_1[i] * fe_0 + ta1_z_xxx_xx_1[i] + ta2_xz_xxx_xx_0[i] * pa_x[i] - ta2_xz_xxx_xx_1[i] * pc_x[i];

        ta2_xz_xxxx_xy_0[i] = 3.0 * ta2_xz_xx_xy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xy_1[i] * fe_0 + ta2_xz_xxx_y_0[i] * fe_0 - ta2_xz_xxx_y_1[i] * fe_0 +
                              ta1_z_xxx_xy_1[i] + ta2_xz_xxx_xy_0[i] * pa_x[i] - ta2_xz_xxx_xy_1[i] * pc_x[i];

        ta2_xz_xxxx_xz_0[i] = 3.0 * ta2_xz_xx_xz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xz_1[i] * fe_0 + ta2_xz_xxx_z_0[i] * fe_0 - ta2_xz_xxx_z_1[i] * fe_0 +
                              ta1_z_xxx_xz_1[i] + ta2_xz_xxx_xz_0[i] * pa_x[i] - ta2_xz_xxx_xz_1[i] * pc_x[i];

        ta2_xz_xxxx_yy_0[i] = 3.0 * ta2_xz_xx_yy_0[i] * fe_0 - 3.0 * ta2_xz_xx_yy_1[i] * fe_0 + ta1_z_xxx_yy_1[i] + ta2_xz_xxx_yy_0[i] * pa_x[i] -
                              ta2_xz_xxx_yy_1[i] * pc_x[i];

        ta2_xz_xxxx_yz_0[i] = 3.0 * ta2_xz_xx_yz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yz_1[i] * fe_0 + ta1_z_xxx_yz_1[i] + ta2_xz_xxx_yz_0[i] * pa_x[i] -
                              ta2_xz_xxx_yz_1[i] * pc_x[i];

        ta2_xz_xxxx_zz_0[i] = 3.0 * ta2_xz_xx_zz_0[i] * fe_0 - 3.0 * ta2_xz_xx_zz_1[i] * fe_0 + ta1_z_xxx_zz_1[i] + ta2_xz_xxx_zz_0[i] * pa_x[i] -
                              ta2_xz_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : GD

    auto ta2_xz_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 186);

    auto ta2_xz_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 187);

    auto ta2_xz_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 188);

    auto ta2_xz_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 189);

    auto ta2_xz_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 190);

    auto ta2_xz_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 191);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xz_xxx_x_0,   \
                             ta2_xz_xxx_x_1,   \
                             ta2_xz_xxx_xx_0,  \
                             ta2_xz_xxx_xx_1,  \
                             ta2_xz_xxx_xy_0,  \
                             ta2_xz_xxx_xy_1,  \
                             ta2_xz_xxx_xz_0,  \
                             ta2_xz_xxx_xz_1,  \
                             ta2_xz_xxx_y_0,   \
                             ta2_xz_xxx_y_1,   \
                             ta2_xz_xxx_yy_0,  \
                             ta2_xz_xxx_yy_1,  \
                             ta2_xz_xxx_yz_0,  \
                             ta2_xz_xxx_yz_1,  \
                             ta2_xz_xxx_z_0,   \
                             ta2_xz_xxx_z_1,   \
                             ta2_xz_xxx_zz_0,  \
                             ta2_xz_xxx_zz_1,  \
                             ta2_xz_xxxy_xx_0, \
                             ta2_xz_xxxy_xy_0, \
                             ta2_xz_xxxy_xz_0, \
                             ta2_xz_xxxy_yy_0, \
                             ta2_xz_xxxy_yz_0, \
                             ta2_xz_xxxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxy_xx_0[i] = ta2_xz_xxx_xx_0[i] * pa_y[i] - ta2_xz_xxx_xx_1[i] * pc_y[i];

        ta2_xz_xxxy_xy_0[i] = ta2_xz_xxx_x_0[i] * fe_0 - ta2_xz_xxx_x_1[i] * fe_0 + ta2_xz_xxx_xy_0[i] * pa_y[i] - ta2_xz_xxx_xy_1[i] * pc_y[i];

        ta2_xz_xxxy_xz_0[i] = ta2_xz_xxx_xz_0[i] * pa_y[i] - ta2_xz_xxx_xz_1[i] * pc_y[i];

        ta2_xz_xxxy_yy_0[i] =
            2.0 * ta2_xz_xxx_y_0[i] * fe_0 - 2.0 * ta2_xz_xxx_y_1[i] * fe_0 + ta2_xz_xxx_yy_0[i] * pa_y[i] - ta2_xz_xxx_yy_1[i] * pc_y[i];

        ta2_xz_xxxy_yz_0[i] = ta2_xz_xxx_z_0[i] * fe_0 - ta2_xz_xxx_z_1[i] * fe_0 + ta2_xz_xxx_yz_0[i] * pa_y[i] - ta2_xz_xxx_yz_1[i] * pc_y[i];

        ta2_xz_xxxy_zz_0[i] = ta2_xz_xxx_zz_0[i] * pa_y[i] - ta2_xz_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 192-198 components of targeted buffer : GD

    auto ta2_xz_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 192);

    auto ta2_xz_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 193);

    auto ta2_xz_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 194);

    auto ta2_xz_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 195);

    auto ta2_xz_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 196);

    auto ta2_xz_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 197);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_z_xxz_yz_1,   \
                             ta1_z_xxz_zz_1,   \
                             ta2_xz_xxx_x_0,   \
                             ta2_xz_xxx_x_1,   \
                             ta2_xz_xxx_xx_0,  \
                             ta2_xz_xxx_xx_1,  \
                             ta2_xz_xxx_xy_0,  \
                             ta2_xz_xxx_xy_1,  \
                             ta2_xz_xxx_xz_0,  \
                             ta2_xz_xxx_xz_1,  \
                             ta2_xz_xxx_yy_0,  \
                             ta2_xz_xxx_yy_1,  \
                             ta2_xz_xxxz_xx_0, \
                             ta2_xz_xxxz_xy_0, \
                             ta2_xz_xxxz_xz_0, \
                             ta2_xz_xxxz_yy_0, \
                             ta2_xz_xxxz_yz_0, \
                             ta2_xz_xxxz_zz_0, \
                             ta2_xz_xxz_yz_0,  \
                             ta2_xz_xxz_yz_1,  \
                             ta2_xz_xxz_zz_0,  \
                             ta2_xz_xxz_zz_1,  \
                             ta2_xz_xz_yz_0,   \
                             ta2_xz_xz_yz_1,   \
                             ta2_xz_xz_zz_0,   \
                             ta2_xz_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxz_xx_0[i] = ta1_x_xxx_xx_1[i] + ta2_xz_xxx_xx_0[i] * pa_z[i] - ta2_xz_xxx_xx_1[i] * pc_z[i];

        ta2_xz_xxxz_xy_0[i] = ta1_x_xxx_xy_1[i] + ta2_xz_xxx_xy_0[i] * pa_z[i] - ta2_xz_xxx_xy_1[i] * pc_z[i];

        ta2_xz_xxxz_xz_0[i] =
            ta2_xz_xxx_x_0[i] * fe_0 - ta2_xz_xxx_x_1[i] * fe_0 + ta1_x_xxx_xz_1[i] + ta2_xz_xxx_xz_0[i] * pa_z[i] - ta2_xz_xxx_xz_1[i] * pc_z[i];

        ta2_xz_xxxz_yy_0[i] = ta1_x_xxx_yy_1[i] + ta2_xz_xxx_yy_0[i] * pa_z[i] - ta2_xz_xxx_yy_1[i] * pc_z[i];

        ta2_xz_xxxz_yz_0[i] = 2.0 * ta2_xz_xz_yz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yz_1[i] * fe_0 + ta1_z_xxz_yz_1[i] + ta2_xz_xxz_yz_0[i] * pa_x[i] -
                              ta2_xz_xxz_yz_1[i] * pc_x[i];

        ta2_xz_xxxz_zz_0[i] = 2.0 * ta2_xz_xz_zz_0[i] * fe_0 - 2.0 * ta2_xz_xz_zz_1[i] * fe_0 + ta1_z_xxz_zz_1[i] + ta2_xz_xxz_zz_0[i] * pa_x[i] -
                              ta2_xz_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 198-204 components of targeted buffer : GD

    auto ta2_xz_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 198);

    auto ta2_xz_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 199);

    auto ta2_xz_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 200);

    auto ta2_xz_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 201);

    auto ta2_xz_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 202);

    auto ta2_xz_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 203);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xyy_yy_1,   \
                             ta1_z_xyy_yz_1,   \
                             ta2_xz_xx_xx_0,   \
                             ta2_xz_xx_xx_1,   \
                             ta2_xz_xx_xy_0,   \
                             ta2_xz_xx_xy_1,   \
                             ta2_xz_xx_xz_0,   \
                             ta2_xz_xx_xz_1,   \
                             ta2_xz_xx_zz_0,   \
                             ta2_xz_xx_zz_1,   \
                             ta2_xz_xxy_x_0,   \
                             ta2_xz_xxy_x_1,   \
                             ta2_xz_xxy_xx_0,  \
                             ta2_xz_xxy_xx_1,  \
                             ta2_xz_xxy_xy_0,  \
                             ta2_xz_xxy_xy_1,  \
                             ta2_xz_xxy_xz_0,  \
                             ta2_xz_xxy_xz_1,  \
                             ta2_xz_xxy_zz_0,  \
                             ta2_xz_xxy_zz_1,  \
                             ta2_xz_xxyy_xx_0, \
                             ta2_xz_xxyy_xy_0, \
                             ta2_xz_xxyy_xz_0, \
                             ta2_xz_xxyy_yy_0, \
                             ta2_xz_xxyy_yz_0, \
                             ta2_xz_xxyy_zz_0, \
                             ta2_xz_xyy_yy_0,  \
                             ta2_xz_xyy_yy_1,  \
                             ta2_xz_xyy_yz_0,  \
                             ta2_xz_xyy_yz_1,  \
                             ta2_xz_yy_yy_0,   \
                             ta2_xz_yy_yy_1,   \
                             ta2_xz_yy_yz_0,   \
                             ta2_xz_yy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyy_xx_0[i] = ta2_xz_xx_xx_0[i] * fe_0 - ta2_xz_xx_xx_1[i] * fe_0 + ta2_xz_xxy_xx_0[i] * pa_y[i] - ta2_xz_xxy_xx_1[i] * pc_y[i];

        ta2_xz_xxyy_xy_0[i] = ta2_xz_xx_xy_0[i] * fe_0 - ta2_xz_xx_xy_1[i] * fe_0 + ta2_xz_xxy_x_0[i] * fe_0 - ta2_xz_xxy_x_1[i] * fe_0 +
                              ta2_xz_xxy_xy_0[i] * pa_y[i] - ta2_xz_xxy_xy_1[i] * pc_y[i];

        ta2_xz_xxyy_xz_0[i] = ta2_xz_xx_xz_0[i] * fe_0 - ta2_xz_xx_xz_1[i] * fe_0 + ta2_xz_xxy_xz_0[i] * pa_y[i] - ta2_xz_xxy_xz_1[i] * pc_y[i];

        ta2_xz_xxyy_yy_0[i] =
            ta2_xz_yy_yy_0[i] * fe_0 - ta2_xz_yy_yy_1[i] * fe_0 + ta1_z_xyy_yy_1[i] + ta2_xz_xyy_yy_0[i] * pa_x[i] - ta2_xz_xyy_yy_1[i] * pc_x[i];

        ta2_xz_xxyy_yz_0[i] =
            ta2_xz_yy_yz_0[i] * fe_0 - ta2_xz_yy_yz_1[i] * fe_0 + ta1_z_xyy_yz_1[i] + ta2_xz_xyy_yz_0[i] * pa_x[i] - ta2_xz_xyy_yz_1[i] * pc_x[i];

        ta2_xz_xxyy_zz_0[i] = ta2_xz_xx_zz_0[i] * fe_0 - ta2_xz_xx_zz_1[i] * fe_0 + ta2_xz_xxy_zz_0[i] * pa_y[i] - ta2_xz_xxy_zz_1[i] * pc_y[i];
    }

    // Set up 204-210 components of targeted buffer : GD

    auto ta2_xz_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 204);

    auto ta2_xz_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 205);

    auto ta2_xz_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 206);

    auto ta2_xz_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 207);

    auto ta2_xz_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 208);

    auto ta2_xz_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 209);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxy_xy_1,   \
                             ta1_x_xxy_yy_1,   \
                             ta2_xz_xxy_xy_0,  \
                             ta2_xz_xxy_xy_1,  \
                             ta2_xz_xxy_yy_0,  \
                             ta2_xz_xxy_yy_1,  \
                             ta2_xz_xxyz_xx_0, \
                             ta2_xz_xxyz_xy_0, \
                             ta2_xz_xxyz_xz_0, \
                             ta2_xz_xxyz_yy_0, \
                             ta2_xz_xxyz_yz_0, \
                             ta2_xz_xxyz_zz_0, \
                             ta2_xz_xxz_xx_0,  \
                             ta2_xz_xxz_xx_1,  \
                             ta2_xz_xxz_xz_0,  \
                             ta2_xz_xxz_xz_1,  \
                             ta2_xz_xxz_yz_0,  \
                             ta2_xz_xxz_yz_1,  \
                             ta2_xz_xxz_z_0,   \
                             ta2_xz_xxz_z_1,   \
                             ta2_xz_xxz_zz_0,  \
                             ta2_xz_xxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyz_xx_0[i] = ta2_xz_xxz_xx_0[i] * pa_y[i] - ta2_xz_xxz_xx_1[i] * pc_y[i];

        ta2_xz_xxyz_xy_0[i] = ta1_x_xxy_xy_1[i] + ta2_xz_xxy_xy_0[i] * pa_z[i] - ta2_xz_xxy_xy_1[i] * pc_z[i];

        ta2_xz_xxyz_xz_0[i] = ta2_xz_xxz_xz_0[i] * pa_y[i] - ta2_xz_xxz_xz_1[i] * pc_y[i];

        ta2_xz_xxyz_yy_0[i] = ta1_x_xxy_yy_1[i] + ta2_xz_xxy_yy_0[i] * pa_z[i] - ta2_xz_xxy_yy_1[i] * pc_z[i];

        ta2_xz_xxyz_yz_0[i] = ta2_xz_xxz_z_0[i] * fe_0 - ta2_xz_xxz_z_1[i] * fe_0 + ta2_xz_xxz_yz_0[i] * pa_y[i] - ta2_xz_xxz_yz_1[i] * pc_y[i];

        ta2_xz_xxyz_zz_0[i] = ta2_xz_xxz_zz_0[i] * pa_y[i] - ta2_xz_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 210-216 components of targeted buffer : GD

    auto ta2_xz_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 210);

    auto ta2_xz_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 211);

    auto ta2_xz_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 212);

    auto ta2_xz_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 213);

    auto ta2_xz_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 214);

    auto ta2_xz_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 215);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xxz_xx_1,   \
                             ta1_x_xxz_xy_1,   \
                             ta1_z_xzz_xz_1,   \
                             ta1_z_xzz_yy_1,   \
                             ta1_z_xzz_yz_1,   \
                             ta1_z_xzz_zz_1,   \
                             ta2_xz_xx_xx_0,   \
                             ta2_xz_xx_xx_1,   \
                             ta2_xz_xx_xy_0,   \
                             ta2_xz_xx_xy_1,   \
                             ta2_xz_xxz_xx_0,  \
                             ta2_xz_xxz_xx_1,  \
                             ta2_xz_xxz_xy_0,  \
                             ta2_xz_xxz_xy_1,  \
                             ta2_xz_xxzz_xx_0, \
                             ta2_xz_xxzz_xy_0, \
                             ta2_xz_xxzz_xz_0, \
                             ta2_xz_xxzz_yy_0, \
                             ta2_xz_xxzz_yz_0, \
                             ta2_xz_xxzz_zz_0, \
                             ta2_xz_xzz_xz_0,  \
                             ta2_xz_xzz_xz_1,  \
                             ta2_xz_xzz_yy_0,  \
                             ta2_xz_xzz_yy_1,  \
                             ta2_xz_xzz_yz_0,  \
                             ta2_xz_xzz_yz_1,  \
                             ta2_xz_xzz_z_0,   \
                             ta2_xz_xzz_z_1,   \
                             ta2_xz_xzz_zz_0,  \
                             ta2_xz_xzz_zz_1,  \
                             ta2_xz_zz_xz_0,   \
                             ta2_xz_zz_xz_1,   \
                             ta2_xz_zz_yy_0,   \
                             ta2_xz_zz_yy_1,   \
                             ta2_xz_zz_yz_0,   \
                             ta2_xz_zz_yz_1,   \
                             ta2_xz_zz_zz_0,   \
                             ta2_xz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxzz_xx_0[i] =
            ta2_xz_xx_xx_0[i] * fe_0 - ta2_xz_xx_xx_1[i] * fe_0 + ta1_x_xxz_xx_1[i] + ta2_xz_xxz_xx_0[i] * pa_z[i] - ta2_xz_xxz_xx_1[i] * pc_z[i];

        ta2_xz_xxzz_xy_0[i] =
            ta2_xz_xx_xy_0[i] * fe_0 - ta2_xz_xx_xy_1[i] * fe_0 + ta1_x_xxz_xy_1[i] + ta2_xz_xxz_xy_0[i] * pa_z[i] - ta2_xz_xxz_xy_1[i] * pc_z[i];

        ta2_xz_xxzz_xz_0[i] = ta2_xz_zz_xz_0[i] * fe_0 - ta2_xz_zz_xz_1[i] * fe_0 + ta2_xz_xzz_z_0[i] * fe_0 - ta2_xz_xzz_z_1[i] * fe_0 +
                              ta1_z_xzz_xz_1[i] + ta2_xz_xzz_xz_0[i] * pa_x[i] - ta2_xz_xzz_xz_1[i] * pc_x[i];

        ta2_xz_xxzz_yy_0[i] =
            ta2_xz_zz_yy_0[i] * fe_0 - ta2_xz_zz_yy_1[i] * fe_0 + ta1_z_xzz_yy_1[i] + ta2_xz_xzz_yy_0[i] * pa_x[i] - ta2_xz_xzz_yy_1[i] * pc_x[i];

        ta2_xz_xxzz_yz_0[i] =
            ta2_xz_zz_yz_0[i] * fe_0 - ta2_xz_zz_yz_1[i] * fe_0 + ta1_z_xzz_yz_1[i] + ta2_xz_xzz_yz_0[i] * pa_x[i] - ta2_xz_xzz_yz_1[i] * pc_x[i];

        ta2_xz_xxzz_zz_0[i] =
            ta2_xz_zz_zz_0[i] * fe_0 - ta2_xz_zz_zz_1[i] * fe_0 + ta1_z_xzz_zz_1[i] + ta2_xz_xzz_zz_0[i] * pa_x[i] - ta2_xz_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 216-222 components of targeted buffer : GD

    auto ta2_xz_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 216);

    auto ta2_xz_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 217);

    auto ta2_xz_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 218);

    auto ta2_xz_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 219);

    auto ta2_xz_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 220);

    auto ta2_xz_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 221);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_zz_1,   \
                             ta2_xz_xy_xx_0,   \
                             ta2_xz_xy_xx_1,   \
                             ta2_xz_xy_xz_0,   \
                             ta2_xz_xy_xz_1,   \
                             ta2_xz_xyy_xx_0,  \
                             ta2_xz_xyy_xx_1,  \
                             ta2_xz_xyy_xz_0,  \
                             ta2_xz_xyy_xz_1,  \
                             ta2_xz_xyyy_xx_0, \
                             ta2_xz_xyyy_xy_0, \
                             ta2_xz_xyyy_xz_0, \
                             ta2_xz_xyyy_yy_0, \
                             ta2_xz_xyyy_yz_0, \
                             ta2_xz_xyyy_zz_0, \
                             ta2_xz_yyy_xy_0,  \
                             ta2_xz_yyy_xy_1,  \
                             ta2_xz_yyy_y_0,   \
                             ta2_xz_yyy_y_1,   \
                             ta2_xz_yyy_yy_0,  \
                             ta2_xz_yyy_yy_1,  \
                             ta2_xz_yyy_yz_0,  \
                             ta2_xz_yyy_yz_1,  \
                             ta2_xz_yyy_zz_0,  \
                             ta2_xz_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyy_xx_0[i] =
            2.0 * ta2_xz_xy_xx_0[i] * fe_0 - 2.0 * ta2_xz_xy_xx_1[i] * fe_0 + ta2_xz_xyy_xx_0[i] * pa_y[i] - ta2_xz_xyy_xx_1[i] * pc_y[i];

        ta2_xz_xyyy_xy_0[i] =
            ta2_xz_yyy_y_0[i] * fe_0 - ta2_xz_yyy_y_1[i] * fe_0 + ta1_z_yyy_xy_1[i] + ta2_xz_yyy_xy_0[i] * pa_x[i] - ta2_xz_yyy_xy_1[i] * pc_x[i];

        ta2_xz_xyyy_xz_0[i] =
            2.0 * ta2_xz_xy_xz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xz_1[i] * fe_0 + ta2_xz_xyy_xz_0[i] * pa_y[i] - ta2_xz_xyy_xz_1[i] * pc_y[i];

        ta2_xz_xyyy_yy_0[i] = ta1_z_yyy_yy_1[i] + ta2_xz_yyy_yy_0[i] * pa_x[i] - ta2_xz_yyy_yy_1[i] * pc_x[i];

        ta2_xz_xyyy_yz_0[i] = ta1_z_yyy_yz_1[i] + ta2_xz_yyy_yz_0[i] * pa_x[i] - ta2_xz_yyy_yz_1[i] * pc_x[i];

        ta2_xz_xyyy_zz_0[i] = ta1_z_yyy_zz_1[i] + ta2_xz_yyy_zz_0[i] * pa_x[i] - ta2_xz_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 222-228 components of targeted buffer : GD

    auto ta2_xz_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 222);

    auto ta2_xz_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 223);

    auto ta2_xz_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 224);

    auto ta2_xz_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 225);

    auto ta2_xz_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 226);

    auto ta2_xz_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 227);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xyy_xx_1,   \
                             ta1_x_xyy_xy_1,   \
                             ta1_z_yyz_yy_1,   \
                             ta1_z_yyz_yz_1,   \
                             ta1_z_yyz_zz_1,   \
                             ta2_xz_xyy_xx_0,  \
                             ta2_xz_xyy_xx_1,  \
                             ta2_xz_xyy_xy_0,  \
                             ta2_xz_xyy_xy_1,  \
                             ta2_xz_xyyz_xx_0, \
                             ta2_xz_xyyz_xy_0, \
                             ta2_xz_xyyz_xz_0, \
                             ta2_xz_xyyz_yy_0, \
                             ta2_xz_xyyz_yz_0, \
                             ta2_xz_xyyz_zz_0, \
                             ta2_xz_xyz_xz_0,  \
                             ta2_xz_xyz_xz_1,  \
                             ta2_xz_xz_xz_0,   \
                             ta2_xz_xz_xz_1,   \
                             ta2_xz_yyz_yy_0,  \
                             ta2_xz_yyz_yy_1,  \
                             ta2_xz_yyz_yz_0,  \
                             ta2_xz_yyz_yz_1,  \
                             ta2_xz_yyz_zz_0,  \
                             ta2_xz_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyz_xx_0[i] = ta1_x_xyy_xx_1[i] + ta2_xz_xyy_xx_0[i] * pa_z[i] - ta2_xz_xyy_xx_1[i] * pc_z[i];

        ta2_xz_xyyz_xy_0[i] = ta1_x_xyy_xy_1[i] + ta2_xz_xyy_xy_0[i] * pa_z[i] - ta2_xz_xyy_xy_1[i] * pc_z[i];

        ta2_xz_xyyz_xz_0[i] = ta2_xz_xz_xz_0[i] * fe_0 - ta2_xz_xz_xz_1[i] * fe_0 + ta2_xz_xyz_xz_0[i] * pa_y[i] - ta2_xz_xyz_xz_1[i] * pc_y[i];

        ta2_xz_xyyz_yy_0[i] = ta1_z_yyz_yy_1[i] + ta2_xz_yyz_yy_0[i] * pa_x[i] - ta2_xz_yyz_yy_1[i] * pc_x[i];

        ta2_xz_xyyz_yz_0[i] = ta1_z_yyz_yz_1[i] + ta2_xz_yyz_yz_0[i] * pa_x[i] - ta2_xz_yyz_yz_1[i] * pc_x[i];

        ta2_xz_xyyz_zz_0[i] = ta1_z_yyz_zz_1[i] + ta2_xz_yyz_zz_0[i] * pa_x[i] - ta2_xz_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 228-234 components of targeted buffer : GD

    auto ta2_xz_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 228);

    auto ta2_xz_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 229);

    auto ta2_xz_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 230);

    auto ta2_xz_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 231);

    auto ta2_xz_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 232);

    auto ta2_xz_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 233);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_yzz_yy_1,   \
                             ta1_z_yzz_yz_1,   \
                             ta2_xz_xyzz_xx_0, \
                             ta2_xz_xyzz_xy_0, \
                             ta2_xz_xyzz_xz_0, \
                             ta2_xz_xyzz_yy_0, \
                             ta2_xz_xyzz_yz_0, \
                             ta2_xz_xyzz_zz_0, \
                             ta2_xz_xzz_x_0,   \
                             ta2_xz_xzz_x_1,   \
                             ta2_xz_xzz_xx_0,  \
                             ta2_xz_xzz_xx_1,  \
                             ta2_xz_xzz_xy_0,  \
                             ta2_xz_xzz_xy_1,  \
                             ta2_xz_xzz_xz_0,  \
                             ta2_xz_xzz_xz_1,  \
                             ta2_xz_xzz_zz_0,  \
                             ta2_xz_xzz_zz_1,  \
                             ta2_xz_yzz_yy_0,  \
                             ta2_xz_yzz_yy_1,  \
                             ta2_xz_yzz_yz_0,  \
                             ta2_xz_yzz_yz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyzz_xx_0[i] = ta2_xz_xzz_xx_0[i] * pa_y[i] - ta2_xz_xzz_xx_1[i] * pc_y[i];

        ta2_xz_xyzz_xy_0[i] = ta2_xz_xzz_x_0[i] * fe_0 - ta2_xz_xzz_x_1[i] * fe_0 + ta2_xz_xzz_xy_0[i] * pa_y[i] - ta2_xz_xzz_xy_1[i] * pc_y[i];

        ta2_xz_xyzz_xz_0[i] = ta2_xz_xzz_xz_0[i] * pa_y[i] - ta2_xz_xzz_xz_1[i] * pc_y[i];

        ta2_xz_xyzz_yy_0[i] = ta1_z_yzz_yy_1[i] + ta2_xz_yzz_yy_0[i] * pa_x[i] - ta2_xz_yzz_yy_1[i] * pc_x[i];

        ta2_xz_xyzz_yz_0[i] = ta1_z_yzz_yz_1[i] + ta2_xz_yzz_yz_0[i] * pa_x[i] - ta2_xz_yzz_yz_1[i] * pc_x[i];

        ta2_xz_xyzz_zz_0[i] = ta2_xz_xzz_zz_0[i] * pa_y[i] - ta2_xz_xzz_zz_1[i] * pc_y[i];
    }

    // Set up 234-240 components of targeted buffer : GD

    auto ta2_xz_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 234);

    auto ta2_xz_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 235);

    auto ta2_xz_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 236);

    auto ta2_xz_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 237);

    auto ta2_xz_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 238);

    auto ta2_xz_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 239);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_1,   \
                             ta2_xz_xzzz_xx_0, \
                             ta2_xz_xzzz_xy_0, \
                             ta2_xz_xzzz_xz_0, \
                             ta2_xz_xzzz_yy_0, \
                             ta2_xz_xzzz_yz_0, \
                             ta2_xz_xzzz_zz_0, \
                             ta2_xz_zzz_x_0,   \
                             ta2_xz_zzz_x_1,   \
                             ta2_xz_zzz_xx_0,  \
                             ta2_xz_zzz_xx_1,  \
                             ta2_xz_zzz_xy_0,  \
                             ta2_xz_zzz_xy_1,  \
                             ta2_xz_zzz_xz_0,  \
                             ta2_xz_zzz_xz_1,  \
                             ta2_xz_zzz_y_0,   \
                             ta2_xz_zzz_y_1,   \
                             ta2_xz_zzz_yy_0,  \
                             ta2_xz_zzz_yy_1,  \
                             ta2_xz_zzz_yz_0,  \
                             ta2_xz_zzz_yz_1,  \
                             ta2_xz_zzz_z_0,   \
                             ta2_xz_zzz_z_1,   \
                             ta2_xz_zzz_zz_0,  \
                             ta2_xz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzzz_xx_0[i] = 2.0 * ta2_xz_zzz_x_0[i] * fe_0 - 2.0 * ta2_xz_zzz_x_1[i] * fe_0 + ta1_z_zzz_xx_1[i] + ta2_xz_zzz_xx_0[i] * pa_x[i] -
                              ta2_xz_zzz_xx_1[i] * pc_x[i];

        ta2_xz_xzzz_xy_0[i] =
            ta2_xz_zzz_y_0[i] * fe_0 - ta2_xz_zzz_y_1[i] * fe_0 + ta1_z_zzz_xy_1[i] + ta2_xz_zzz_xy_0[i] * pa_x[i] - ta2_xz_zzz_xy_1[i] * pc_x[i];

        ta2_xz_xzzz_xz_0[i] =
            ta2_xz_zzz_z_0[i] * fe_0 - ta2_xz_zzz_z_1[i] * fe_0 + ta1_z_zzz_xz_1[i] + ta2_xz_zzz_xz_0[i] * pa_x[i] - ta2_xz_zzz_xz_1[i] * pc_x[i];

        ta2_xz_xzzz_yy_0[i] = ta1_z_zzz_yy_1[i] + ta2_xz_zzz_yy_0[i] * pa_x[i] - ta2_xz_zzz_yy_1[i] * pc_x[i];

        ta2_xz_xzzz_yz_0[i] = ta1_z_zzz_yz_1[i] + ta2_xz_zzz_yz_0[i] * pa_x[i] - ta2_xz_zzz_yz_1[i] * pc_x[i];

        ta2_xz_xzzz_zz_0[i] = ta1_z_zzz_zz_1[i] + ta2_xz_zzz_zz_0[i] * pa_x[i] - ta2_xz_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 240-246 components of targeted buffer : GD

    auto ta2_xz_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 240);

    auto ta2_xz_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 241);

    auto ta2_xz_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 242);

    auto ta2_xz_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 243);

    auto ta2_xz_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 244);

    auto ta2_xz_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 245);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xz_yy_xx_0,   \
                             ta2_xz_yy_xx_1,   \
                             ta2_xz_yy_xy_0,   \
                             ta2_xz_yy_xy_1,   \
                             ta2_xz_yy_xz_0,   \
                             ta2_xz_yy_xz_1,   \
                             ta2_xz_yy_yy_0,   \
                             ta2_xz_yy_yy_1,   \
                             ta2_xz_yy_yz_0,   \
                             ta2_xz_yy_yz_1,   \
                             ta2_xz_yy_zz_0,   \
                             ta2_xz_yy_zz_1,   \
                             ta2_xz_yyy_x_0,   \
                             ta2_xz_yyy_x_1,   \
                             ta2_xz_yyy_xx_0,  \
                             ta2_xz_yyy_xx_1,  \
                             ta2_xz_yyy_xy_0,  \
                             ta2_xz_yyy_xy_1,  \
                             ta2_xz_yyy_xz_0,  \
                             ta2_xz_yyy_xz_1,  \
                             ta2_xz_yyy_y_0,   \
                             ta2_xz_yyy_y_1,   \
                             ta2_xz_yyy_yy_0,  \
                             ta2_xz_yyy_yy_1,  \
                             ta2_xz_yyy_yz_0,  \
                             ta2_xz_yyy_yz_1,  \
                             ta2_xz_yyy_z_0,   \
                             ta2_xz_yyy_z_1,   \
                             ta2_xz_yyy_zz_0,  \
                             ta2_xz_yyy_zz_1,  \
                             ta2_xz_yyyy_xx_0, \
                             ta2_xz_yyyy_xy_0, \
                             ta2_xz_yyyy_xz_0, \
                             ta2_xz_yyyy_yy_0, \
                             ta2_xz_yyyy_yz_0, \
                             ta2_xz_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyy_xx_0[i] =
            3.0 * ta2_xz_yy_xx_0[i] * fe_0 - 3.0 * ta2_xz_yy_xx_1[i] * fe_0 + ta2_xz_yyy_xx_0[i] * pa_y[i] - ta2_xz_yyy_xx_1[i] * pc_y[i];

        ta2_xz_yyyy_xy_0[i] = 3.0 * ta2_xz_yy_xy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xy_1[i] * fe_0 + ta2_xz_yyy_x_0[i] * fe_0 - ta2_xz_yyy_x_1[i] * fe_0 +
                              ta2_xz_yyy_xy_0[i] * pa_y[i] - ta2_xz_yyy_xy_1[i] * pc_y[i];

        ta2_xz_yyyy_xz_0[i] =
            3.0 * ta2_xz_yy_xz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xz_1[i] * fe_0 + ta2_xz_yyy_xz_0[i] * pa_y[i] - ta2_xz_yyy_xz_1[i] * pc_y[i];

        ta2_xz_yyyy_yy_0[i] = 3.0 * ta2_xz_yy_yy_0[i] * fe_0 - 3.0 * ta2_xz_yy_yy_1[i] * fe_0 + 2.0 * ta2_xz_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_xz_yyy_y_1[i] * fe_0 + ta2_xz_yyy_yy_0[i] * pa_y[i] - ta2_xz_yyy_yy_1[i] * pc_y[i];

        ta2_xz_yyyy_yz_0[i] = 3.0 * ta2_xz_yy_yz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yz_1[i] * fe_0 + ta2_xz_yyy_z_0[i] * fe_0 - ta2_xz_yyy_z_1[i] * fe_0 +
                              ta2_xz_yyy_yz_0[i] * pa_y[i] - ta2_xz_yyy_yz_1[i] * pc_y[i];

        ta2_xz_yyyy_zz_0[i] =
            3.0 * ta2_xz_yy_zz_0[i] * fe_0 - 3.0 * ta2_xz_yy_zz_1[i] * fe_0 + ta2_xz_yyy_zz_0[i] * pa_y[i] - ta2_xz_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 246-252 components of targeted buffer : GD

    auto ta2_xz_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 246);

    auto ta2_xz_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 247);

    auto ta2_xz_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 248);

    auto ta2_xz_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 249);

    auto ta2_xz_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 250);

    auto ta2_xz_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 251);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyy_xx_1,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yz_1,   \
                             ta2_xz_yyy_xx_0,  \
                             ta2_xz_yyy_xx_1,  \
                             ta2_xz_yyy_xy_0,  \
                             ta2_xz_yyy_xy_1,  \
                             ta2_xz_yyy_y_0,   \
                             ta2_xz_yyy_y_1,   \
                             ta2_xz_yyy_yy_0,  \
                             ta2_xz_yyy_yy_1,  \
                             ta2_xz_yyy_yz_0,  \
                             ta2_xz_yyy_yz_1,  \
                             ta2_xz_yyyz_xx_0, \
                             ta2_xz_yyyz_xy_0, \
                             ta2_xz_yyyz_xz_0, \
                             ta2_xz_yyyz_yy_0, \
                             ta2_xz_yyyz_yz_0, \
                             ta2_xz_yyyz_zz_0, \
                             ta2_xz_yyz_xz_0,  \
                             ta2_xz_yyz_xz_1,  \
                             ta2_xz_yyz_zz_0,  \
                             ta2_xz_yyz_zz_1,  \
                             ta2_xz_yz_xz_0,   \
                             ta2_xz_yz_xz_1,   \
                             ta2_xz_yz_zz_0,   \
                             ta2_xz_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyz_xx_0[i] = ta1_x_yyy_xx_1[i] + ta2_xz_yyy_xx_0[i] * pa_z[i] - ta2_xz_yyy_xx_1[i] * pc_z[i];

        ta2_xz_yyyz_xy_0[i] = ta1_x_yyy_xy_1[i] + ta2_xz_yyy_xy_0[i] * pa_z[i] - ta2_xz_yyy_xy_1[i] * pc_z[i];

        ta2_xz_yyyz_xz_0[i] =
            2.0 * ta2_xz_yz_xz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xz_1[i] * fe_0 + ta2_xz_yyz_xz_0[i] * pa_y[i] - ta2_xz_yyz_xz_1[i] * pc_y[i];

        ta2_xz_yyyz_yy_0[i] = ta1_x_yyy_yy_1[i] + ta2_xz_yyy_yy_0[i] * pa_z[i] - ta2_xz_yyy_yy_1[i] * pc_z[i];

        ta2_xz_yyyz_yz_0[i] =
            ta2_xz_yyy_y_0[i] * fe_0 - ta2_xz_yyy_y_1[i] * fe_0 + ta1_x_yyy_yz_1[i] + ta2_xz_yyy_yz_0[i] * pa_z[i] - ta2_xz_yyy_yz_1[i] * pc_z[i];

        ta2_xz_yyyz_zz_0[i] =
            2.0 * ta2_xz_yz_zz_0[i] * fe_0 - 2.0 * ta2_xz_yz_zz_1[i] * fe_0 + ta2_xz_yyz_zz_0[i] * pa_y[i] - ta2_xz_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 252-258 components of targeted buffer : GD

    auto ta2_xz_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 252);

    auto ta2_xz_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 253);

    auto ta2_xz_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 254);

    auto ta2_xz_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 255);

    auto ta2_xz_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 256);

    auto ta2_xz_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 257);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyz_xy_1,   \
                             ta1_x_yyz_yy_1,   \
                             ta2_xz_yy_xy_0,   \
                             ta2_xz_yy_xy_1,   \
                             ta2_xz_yy_yy_0,   \
                             ta2_xz_yy_yy_1,   \
                             ta2_xz_yyz_xy_0,  \
                             ta2_xz_yyz_xy_1,  \
                             ta2_xz_yyz_yy_0,  \
                             ta2_xz_yyz_yy_1,  \
                             ta2_xz_yyzz_xx_0, \
                             ta2_xz_yyzz_xy_0, \
                             ta2_xz_yyzz_xz_0, \
                             ta2_xz_yyzz_yy_0, \
                             ta2_xz_yyzz_yz_0, \
                             ta2_xz_yyzz_zz_0, \
                             ta2_xz_yzz_xx_0,  \
                             ta2_xz_yzz_xx_1,  \
                             ta2_xz_yzz_xz_0,  \
                             ta2_xz_yzz_xz_1,  \
                             ta2_xz_yzz_yz_0,  \
                             ta2_xz_yzz_yz_1,  \
                             ta2_xz_yzz_z_0,   \
                             ta2_xz_yzz_z_1,   \
                             ta2_xz_yzz_zz_0,  \
                             ta2_xz_yzz_zz_1,  \
                             ta2_xz_zz_xx_0,   \
                             ta2_xz_zz_xx_1,   \
                             ta2_xz_zz_xz_0,   \
                             ta2_xz_zz_xz_1,   \
                             ta2_xz_zz_yz_0,   \
                             ta2_xz_zz_yz_1,   \
                             ta2_xz_zz_zz_0,   \
                             ta2_xz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyzz_xx_0[i] = ta2_xz_zz_xx_0[i] * fe_0 - ta2_xz_zz_xx_1[i] * fe_0 + ta2_xz_yzz_xx_0[i] * pa_y[i] - ta2_xz_yzz_xx_1[i] * pc_y[i];

        ta2_xz_yyzz_xy_0[i] =
            ta2_xz_yy_xy_0[i] * fe_0 - ta2_xz_yy_xy_1[i] * fe_0 + ta1_x_yyz_xy_1[i] + ta2_xz_yyz_xy_0[i] * pa_z[i] - ta2_xz_yyz_xy_1[i] * pc_z[i];

        ta2_xz_yyzz_xz_0[i] = ta2_xz_zz_xz_0[i] * fe_0 - ta2_xz_zz_xz_1[i] * fe_0 + ta2_xz_yzz_xz_0[i] * pa_y[i] - ta2_xz_yzz_xz_1[i] * pc_y[i];

        ta2_xz_yyzz_yy_0[i] =
            ta2_xz_yy_yy_0[i] * fe_0 - ta2_xz_yy_yy_1[i] * fe_0 + ta1_x_yyz_yy_1[i] + ta2_xz_yyz_yy_0[i] * pa_z[i] - ta2_xz_yyz_yy_1[i] * pc_z[i];

        ta2_xz_yyzz_yz_0[i] = ta2_xz_zz_yz_0[i] * fe_0 - ta2_xz_zz_yz_1[i] * fe_0 + ta2_xz_yzz_z_0[i] * fe_0 - ta2_xz_yzz_z_1[i] * fe_0 +
                              ta2_xz_yzz_yz_0[i] * pa_y[i] - ta2_xz_yzz_yz_1[i] * pc_y[i];

        ta2_xz_yyzz_zz_0[i] = ta2_xz_zz_zz_0[i] * fe_0 - ta2_xz_zz_zz_1[i] * fe_0 + ta2_xz_yzz_zz_0[i] * pa_y[i] - ta2_xz_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 258-264 components of targeted buffer : GD

    auto ta2_xz_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 258);

    auto ta2_xz_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 259);

    auto ta2_xz_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 260);

    auto ta2_xz_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 261);

    auto ta2_xz_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 262);

    auto ta2_xz_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 263);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xz_yzzz_xx_0, \
                             ta2_xz_yzzz_xy_0, \
                             ta2_xz_yzzz_xz_0, \
                             ta2_xz_yzzz_yy_0, \
                             ta2_xz_yzzz_yz_0, \
                             ta2_xz_yzzz_zz_0, \
                             ta2_xz_zzz_x_0,   \
                             ta2_xz_zzz_x_1,   \
                             ta2_xz_zzz_xx_0,  \
                             ta2_xz_zzz_xx_1,  \
                             ta2_xz_zzz_xy_0,  \
                             ta2_xz_zzz_xy_1,  \
                             ta2_xz_zzz_xz_0,  \
                             ta2_xz_zzz_xz_1,  \
                             ta2_xz_zzz_y_0,   \
                             ta2_xz_zzz_y_1,   \
                             ta2_xz_zzz_yy_0,  \
                             ta2_xz_zzz_yy_1,  \
                             ta2_xz_zzz_yz_0,  \
                             ta2_xz_zzz_yz_1,  \
                             ta2_xz_zzz_z_0,   \
                             ta2_xz_zzz_z_1,   \
                             ta2_xz_zzz_zz_0,  \
                             ta2_xz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzzz_xx_0[i] = ta2_xz_zzz_xx_0[i] * pa_y[i] - ta2_xz_zzz_xx_1[i] * pc_y[i];

        ta2_xz_yzzz_xy_0[i] = ta2_xz_zzz_x_0[i] * fe_0 - ta2_xz_zzz_x_1[i] * fe_0 + ta2_xz_zzz_xy_0[i] * pa_y[i] - ta2_xz_zzz_xy_1[i] * pc_y[i];

        ta2_xz_yzzz_xz_0[i] = ta2_xz_zzz_xz_0[i] * pa_y[i] - ta2_xz_zzz_xz_1[i] * pc_y[i];

        ta2_xz_yzzz_yy_0[i] =
            2.0 * ta2_xz_zzz_y_0[i] * fe_0 - 2.0 * ta2_xz_zzz_y_1[i] * fe_0 + ta2_xz_zzz_yy_0[i] * pa_y[i] - ta2_xz_zzz_yy_1[i] * pc_y[i];

        ta2_xz_yzzz_yz_0[i] = ta2_xz_zzz_z_0[i] * fe_0 - ta2_xz_zzz_z_1[i] * fe_0 + ta2_xz_zzz_yz_0[i] * pa_y[i] - ta2_xz_zzz_yz_1[i] * pc_y[i];

        ta2_xz_yzzz_zz_0[i] = ta2_xz_zzz_zz_0[i] * pa_y[i] - ta2_xz_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 264-270 components of targeted buffer : GD

    auto ta2_xz_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 264);

    auto ta2_xz_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 265);

    auto ta2_xz_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 266);

    auto ta2_xz_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 267);

    auto ta2_xz_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 268);

    auto ta2_xz_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 269);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xy_1,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_yy_1,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_1,   \
                             ta2_xz_zz_xx_0,   \
                             ta2_xz_zz_xx_1,   \
                             ta2_xz_zz_xy_0,   \
                             ta2_xz_zz_xy_1,   \
                             ta2_xz_zz_xz_0,   \
                             ta2_xz_zz_xz_1,   \
                             ta2_xz_zz_yy_0,   \
                             ta2_xz_zz_yy_1,   \
                             ta2_xz_zz_yz_0,   \
                             ta2_xz_zz_yz_1,   \
                             ta2_xz_zz_zz_0,   \
                             ta2_xz_zz_zz_1,   \
                             ta2_xz_zzz_x_0,   \
                             ta2_xz_zzz_x_1,   \
                             ta2_xz_zzz_xx_0,  \
                             ta2_xz_zzz_xx_1,  \
                             ta2_xz_zzz_xy_0,  \
                             ta2_xz_zzz_xy_1,  \
                             ta2_xz_zzz_xz_0,  \
                             ta2_xz_zzz_xz_1,  \
                             ta2_xz_zzz_y_0,   \
                             ta2_xz_zzz_y_1,   \
                             ta2_xz_zzz_yy_0,  \
                             ta2_xz_zzz_yy_1,  \
                             ta2_xz_zzz_yz_0,  \
                             ta2_xz_zzz_yz_1,  \
                             ta2_xz_zzz_z_0,   \
                             ta2_xz_zzz_z_1,   \
                             ta2_xz_zzz_zz_0,  \
                             ta2_xz_zzz_zz_1,  \
                             ta2_xz_zzzz_xx_0, \
                             ta2_xz_zzzz_xy_0, \
                             ta2_xz_zzzz_xz_0, \
                             ta2_xz_zzzz_yy_0, \
                             ta2_xz_zzzz_yz_0, \
                             ta2_xz_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzzz_xx_0[i] = 3.0 * ta2_xz_zz_xx_0[i] * fe_0 - 3.0 * ta2_xz_zz_xx_1[i] * fe_0 + ta1_x_zzz_xx_1[i] + ta2_xz_zzz_xx_0[i] * pa_z[i] -
                              ta2_xz_zzz_xx_1[i] * pc_z[i];

        ta2_xz_zzzz_xy_0[i] = 3.0 * ta2_xz_zz_xy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xy_1[i] * fe_0 + ta1_x_zzz_xy_1[i] + ta2_xz_zzz_xy_0[i] * pa_z[i] -
                              ta2_xz_zzz_xy_1[i] * pc_z[i];

        ta2_xz_zzzz_xz_0[i] = 3.0 * ta2_xz_zz_xz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xz_1[i] * fe_0 + ta2_xz_zzz_x_0[i] * fe_0 - ta2_xz_zzz_x_1[i] * fe_0 +
                              ta1_x_zzz_xz_1[i] + ta2_xz_zzz_xz_0[i] * pa_z[i] - ta2_xz_zzz_xz_1[i] * pc_z[i];

        ta2_xz_zzzz_yy_0[i] = 3.0 * ta2_xz_zz_yy_0[i] * fe_0 - 3.0 * ta2_xz_zz_yy_1[i] * fe_0 + ta1_x_zzz_yy_1[i] + ta2_xz_zzz_yy_0[i] * pa_z[i] -
                              ta2_xz_zzz_yy_1[i] * pc_z[i];

        ta2_xz_zzzz_yz_0[i] = 3.0 * ta2_xz_zz_yz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yz_1[i] * fe_0 + ta2_xz_zzz_y_0[i] * fe_0 - ta2_xz_zzz_y_1[i] * fe_0 +
                              ta1_x_zzz_yz_1[i] + ta2_xz_zzz_yz_0[i] * pa_z[i] - ta2_xz_zzz_yz_1[i] * pc_z[i];

        ta2_xz_zzzz_zz_0[i] = 3.0 * ta2_xz_zz_zz_0[i] * fe_0 - 3.0 * ta2_xz_zz_zz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_xz_zzz_z_1[i] * fe_0 + ta1_x_zzz_zz_1[i] + ta2_xz_zzz_zz_0[i] * pa_z[i] - ta2_xz_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 270-276 components of targeted buffer : GD

    auto ta2_yy_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 270);

    auto ta2_yy_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 271);

    auto ta2_yy_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 272);

    auto ta2_yy_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 273);

    auto ta2_yy_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 274);

    auto ta2_yy_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 275);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yy_xx_xx_0,   \
                             ta2_yy_xx_xx_1,   \
                             ta2_yy_xx_xy_0,   \
                             ta2_yy_xx_xy_1,   \
                             ta2_yy_xx_xz_0,   \
                             ta2_yy_xx_xz_1,   \
                             ta2_yy_xx_yy_0,   \
                             ta2_yy_xx_yy_1,   \
                             ta2_yy_xx_yz_0,   \
                             ta2_yy_xx_yz_1,   \
                             ta2_yy_xx_zz_0,   \
                             ta2_yy_xx_zz_1,   \
                             ta2_yy_xxx_x_0,   \
                             ta2_yy_xxx_x_1,   \
                             ta2_yy_xxx_xx_0,  \
                             ta2_yy_xxx_xx_1,  \
                             ta2_yy_xxx_xy_0,  \
                             ta2_yy_xxx_xy_1,  \
                             ta2_yy_xxx_xz_0,  \
                             ta2_yy_xxx_xz_1,  \
                             ta2_yy_xxx_y_0,   \
                             ta2_yy_xxx_y_1,   \
                             ta2_yy_xxx_yy_0,  \
                             ta2_yy_xxx_yy_1,  \
                             ta2_yy_xxx_yz_0,  \
                             ta2_yy_xxx_yz_1,  \
                             ta2_yy_xxx_z_0,   \
                             ta2_yy_xxx_z_1,   \
                             ta2_yy_xxx_zz_0,  \
                             ta2_yy_xxx_zz_1,  \
                             ta2_yy_xxxx_xx_0, \
                             ta2_yy_xxxx_xy_0, \
                             ta2_yy_xxxx_xz_0, \
                             ta2_yy_xxxx_yy_0, \
                             ta2_yy_xxxx_yz_0, \
                             ta2_yy_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxx_xx_0[i] = 3.0 * ta2_yy_xx_xx_0[i] * fe_0 - 3.0 * ta2_yy_xx_xx_1[i] * fe_0 + 2.0 * ta2_yy_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_yy_xxx_x_1[i] * fe_0 + ta2_yy_xxx_xx_0[i] * pa_x[i] - ta2_yy_xxx_xx_1[i] * pc_x[i];

        ta2_yy_xxxx_xy_0[i] = 3.0 * ta2_yy_xx_xy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xy_1[i] * fe_0 + ta2_yy_xxx_y_0[i] * fe_0 - ta2_yy_xxx_y_1[i] * fe_0 +
                              ta2_yy_xxx_xy_0[i] * pa_x[i] - ta2_yy_xxx_xy_1[i] * pc_x[i];

        ta2_yy_xxxx_xz_0[i] = 3.0 * ta2_yy_xx_xz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xz_1[i] * fe_0 + ta2_yy_xxx_z_0[i] * fe_0 - ta2_yy_xxx_z_1[i] * fe_0 +
                              ta2_yy_xxx_xz_0[i] * pa_x[i] - ta2_yy_xxx_xz_1[i] * pc_x[i];

        ta2_yy_xxxx_yy_0[i] =
            3.0 * ta2_yy_xx_yy_0[i] * fe_0 - 3.0 * ta2_yy_xx_yy_1[i] * fe_0 + ta2_yy_xxx_yy_0[i] * pa_x[i] - ta2_yy_xxx_yy_1[i] * pc_x[i];

        ta2_yy_xxxx_yz_0[i] =
            3.0 * ta2_yy_xx_yz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yz_1[i] * fe_0 + ta2_yy_xxx_yz_0[i] * pa_x[i] - ta2_yy_xxx_yz_1[i] * pc_x[i];

        ta2_yy_xxxx_zz_0[i] =
            3.0 * ta2_yy_xx_zz_0[i] * fe_0 - 3.0 * ta2_yy_xx_zz_1[i] * fe_0 + ta2_yy_xxx_zz_0[i] * pa_x[i] - ta2_yy_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 276-282 components of targeted buffer : GD

    auto ta2_yy_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 276);

    auto ta2_yy_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 277);

    auto ta2_yy_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 278);

    auto ta2_yy_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 279);

    auto ta2_yy_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 280);

    auto ta2_yy_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 281);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_zz_1,   \
                             ta2_yy_xxx_x_0,   \
                             ta2_yy_xxx_x_1,   \
                             ta2_yy_xxx_xx_0,  \
                             ta2_yy_xxx_xx_1,  \
                             ta2_yy_xxx_xy_0,  \
                             ta2_yy_xxx_xy_1,  \
                             ta2_yy_xxx_xz_0,  \
                             ta2_yy_xxx_xz_1,  \
                             ta2_yy_xxx_zz_0,  \
                             ta2_yy_xxx_zz_1,  \
                             ta2_yy_xxxy_xx_0, \
                             ta2_yy_xxxy_xy_0, \
                             ta2_yy_xxxy_xz_0, \
                             ta2_yy_xxxy_yy_0, \
                             ta2_yy_xxxy_yz_0, \
                             ta2_yy_xxxy_zz_0, \
                             ta2_yy_xxy_yy_0,  \
                             ta2_yy_xxy_yy_1,  \
                             ta2_yy_xxy_yz_0,  \
                             ta2_yy_xxy_yz_1,  \
                             ta2_yy_xy_yy_0,   \
                             ta2_yy_xy_yy_1,   \
                             ta2_yy_xy_yz_0,   \
                             ta2_yy_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxy_xx_0[i] = 2.0 * ta1_y_xxx_xx_1[i] + ta2_yy_xxx_xx_0[i] * pa_y[i] - ta2_yy_xxx_xx_1[i] * pc_y[i];

        ta2_yy_xxxy_xy_0[i] = ta2_yy_xxx_x_0[i] * fe_0 - ta2_yy_xxx_x_1[i] * fe_0 + 2.0 * ta1_y_xxx_xy_1[i] + ta2_yy_xxx_xy_0[i] * pa_y[i] -
                              ta2_yy_xxx_xy_1[i] * pc_y[i];

        ta2_yy_xxxy_xz_0[i] = 2.0 * ta1_y_xxx_xz_1[i] + ta2_yy_xxx_xz_0[i] * pa_y[i] - ta2_yy_xxx_xz_1[i] * pc_y[i];

        ta2_yy_xxxy_yy_0[i] =
            2.0 * ta2_yy_xy_yy_0[i] * fe_0 - 2.0 * ta2_yy_xy_yy_1[i] * fe_0 + ta2_yy_xxy_yy_0[i] * pa_x[i] - ta2_yy_xxy_yy_1[i] * pc_x[i];

        ta2_yy_xxxy_yz_0[i] =
            2.0 * ta2_yy_xy_yz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yz_1[i] * fe_0 + ta2_yy_xxy_yz_0[i] * pa_x[i] - ta2_yy_xxy_yz_1[i] * pc_x[i];

        ta2_yy_xxxy_zz_0[i] = 2.0 * ta1_y_xxx_zz_1[i] + ta2_yy_xxx_zz_0[i] * pa_y[i] - ta2_yy_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 282-288 components of targeted buffer : GD

    auto ta2_yy_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 282);

    auto ta2_yy_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 283);

    auto ta2_yy_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 284);

    auto ta2_yy_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 285);

    auto ta2_yy_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 286);

    auto ta2_yy_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 287);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta2_yy_xxx_x_0,   \
                             ta2_yy_xxx_x_1,   \
                             ta2_yy_xxx_xx_0,  \
                             ta2_yy_xxx_xx_1,  \
                             ta2_yy_xxx_xy_0,  \
                             ta2_yy_xxx_xy_1,  \
                             ta2_yy_xxx_xz_0,  \
                             ta2_yy_xxx_xz_1,  \
                             ta2_yy_xxx_yy_0,  \
                             ta2_yy_xxx_yy_1,  \
                             ta2_yy_xxxz_xx_0, \
                             ta2_yy_xxxz_xy_0, \
                             ta2_yy_xxxz_xz_0, \
                             ta2_yy_xxxz_yy_0, \
                             ta2_yy_xxxz_yz_0, \
                             ta2_yy_xxxz_zz_0, \
                             ta2_yy_xxz_yz_0,  \
                             ta2_yy_xxz_yz_1,  \
                             ta2_yy_xxz_zz_0,  \
                             ta2_yy_xxz_zz_1,  \
                             ta2_yy_xz_yz_0,   \
                             ta2_yy_xz_yz_1,   \
                             ta2_yy_xz_zz_0,   \
                             ta2_yy_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxz_xx_0[i] = ta2_yy_xxx_xx_0[i] * pa_z[i] - ta2_yy_xxx_xx_1[i] * pc_z[i];

        ta2_yy_xxxz_xy_0[i] = ta2_yy_xxx_xy_0[i] * pa_z[i] - ta2_yy_xxx_xy_1[i] * pc_z[i];

        ta2_yy_xxxz_xz_0[i] = ta2_yy_xxx_x_0[i] * fe_0 - ta2_yy_xxx_x_1[i] * fe_0 + ta2_yy_xxx_xz_0[i] * pa_z[i] - ta2_yy_xxx_xz_1[i] * pc_z[i];

        ta2_yy_xxxz_yy_0[i] = ta2_yy_xxx_yy_0[i] * pa_z[i] - ta2_yy_xxx_yy_1[i] * pc_z[i];

        ta2_yy_xxxz_yz_0[i] =
            2.0 * ta2_yy_xz_yz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yz_1[i] * fe_0 + ta2_yy_xxz_yz_0[i] * pa_x[i] - ta2_yy_xxz_yz_1[i] * pc_x[i];

        ta2_yy_xxxz_zz_0[i] =
            2.0 * ta2_yy_xz_zz_0[i] * fe_0 - 2.0 * ta2_yy_xz_zz_1[i] * fe_0 + ta2_yy_xxz_zz_0[i] * pa_x[i] - ta2_yy_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 288-294 components of targeted buffer : GD

    auto ta2_yy_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 288);

    auto ta2_yy_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 289);

    auto ta2_yy_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 290);

    auto ta2_yy_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 291);

    auto ta2_yy_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 292);

    auto ta2_yy_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 293);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxy_xx_1,   \
                             ta1_y_xxy_xz_1,   \
                             ta2_yy_xx_xx_0,   \
                             ta2_yy_xx_xx_1,   \
                             ta2_yy_xx_xz_0,   \
                             ta2_yy_xx_xz_1,   \
                             ta2_yy_xxy_xx_0,  \
                             ta2_yy_xxy_xx_1,  \
                             ta2_yy_xxy_xz_0,  \
                             ta2_yy_xxy_xz_1,  \
                             ta2_yy_xxyy_xx_0, \
                             ta2_yy_xxyy_xy_0, \
                             ta2_yy_xxyy_xz_0, \
                             ta2_yy_xxyy_yy_0, \
                             ta2_yy_xxyy_yz_0, \
                             ta2_yy_xxyy_zz_0, \
                             ta2_yy_xyy_xy_0,  \
                             ta2_yy_xyy_xy_1,  \
                             ta2_yy_xyy_y_0,   \
                             ta2_yy_xyy_y_1,   \
                             ta2_yy_xyy_yy_0,  \
                             ta2_yy_xyy_yy_1,  \
                             ta2_yy_xyy_yz_0,  \
                             ta2_yy_xyy_yz_1,  \
                             ta2_yy_xyy_zz_0,  \
                             ta2_yy_xyy_zz_1,  \
                             ta2_yy_yy_xy_0,   \
                             ta2_yy_yy_xy_1,   \
                             ta2_yy_yy_yy_0,   \
                             ta2_yy_yy_yy_1,   \
                             ta2_yy_yy_yz_0,   \
                             ta2_yy_yy_yz_1,   \
                             ta2_yy_yy_zz_0,   \
                             ta2_yy_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyy_xx_0[i] = ta2_yy_xx_xx_0[i] * fe_0 - ta2_yy_xx_xx_1[i] * fe_0 + 2.0 * ta1_y_xxy_xx_1[i] + ta2_yy_xxy_xx_0[i] * pa_y[i] -
                              ta2_yy_xxy_xx_1[i] * pc_y[i];

        ta2_yy_xxyy_xy_0[i] = ta2_yy_yy_xy_0[i] * fe_0 - ta2_yy_yy_xy_1[i] * fe_0 + ta2_yy_xyy_y_0[i] * fe_0 - ta2_yy_xyy_y_1[i] * fe_0 +
                              ta2_yy_xyy_xy_0[i] * pa_x[i] - ta2_yy_xyy_xy_1[i] * pc_x[i];

        ta2_yy_xxyy_xz_0[i] = ta2_yy_xx_xz_0[i] * fe_0 - ta2_yy_xx_xz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xz_1[i] + ta2_yy_xxy_xz_0[i] * pa_y[i] -
                              ta2_yy_xxy_xz_1[i] * pc_y[i];

        ta2_yy_xxyy_yy_0[i] = ta2_yy_yy_yy_0[i] * fe_0 - ta2_yy_yy_yy_1[i] * fe_0 + ta2_yy_xyy_yy_0[i] * pa_x[i] - ta2_yy_xyy_yy_1[i] * pc_x[i];

        ta2_yy_xxyy_yz_0[i] = ta2_yy_yy_yz_0[i] * fe_0 - ta2_yy_yy_yz_1[i] * fe_0 + ta2_yy_xyy_yz_0[i] * pa_x[i] - ta2_yy_xyy_yz_1[i] * pc_x[i];

        ta2_yy_xxyy_zz_0[i] = ta2_yy_yy_zz_0[i] * fe_0 - ta2_yy_yy_zz_1[i] * fe_0 + ta2_yy_xyy_zz_0[i] * pa_x[i] - ta2_yy_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 294-300 components of targeted buffer : GD

    auto ta2_yy_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 294);

    auto ta2_yy_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 295);

    auto ta2_yy_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 296);

    auto ta2_yy_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 297);

    auto ta2_yy_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 298);

    auto ta2_yy_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 299);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xxz_xz_1,   \
                             ta1_y_xxz_zz_1,   \
                             ta2_yy_xxy_xx_0,  \
                             ta2_yy_xxy_xx_1,  \
                             ta2_yy_xxy_xy_0,  \
                             ta2_yy_xxy_xy_1,  \
                             ta2_yy_xxy_yy_0,  \
                             ta2_yy_xxy_yy_1,  \
                             ta2_yy_xxyz_xx_0, \
                             ta2_yy_xxyz_xy_0, \
                             ta2_yy_xxyz_xz_0, \
                             ta2_yy_xxyz_yy_0, \
                             ta2_yy_xxyz_yz_0, \
                             ta2_yy_xxyz_zz_0, \
                             ta2_yy_xxz_xz_0,  \
                             ta2_yy_xxz_xz_1,  \
                             ta2_yy_xxz_zz_0,  \
                             ta2_yy_xxz_zz_1,  \
                             ta2_yy_xyz_yz_0,  \
                             ta2_yy_xyz_yz_1,  \
                             ta2_yy_yz_yz_0,   \
                             ta2_yy_yz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyz_xx_0[i] = ta2_yy_xxy_xx_0[i] * pa_z[i] - ta2_yy_xxy_xx_1[i] * pc_z[i];

        ta2_yy_xxyz_xy_0[i] = ta2_yy_xxy_xy_0[i] * pa_z[i] - ta2_yy_xxy_xy_1[i] * pc_z[i];

        ta2_yy_xxyz_xz_0[i] = 2.0 * ta1_y_xxz_xz_1[i] + ta2_yy_xxz_xz_0[i] * pa_y[i] - ta2_yy_xxz_xz_1[i] * pc_y[i];

        ta2_yy_xxyz_yy_0[i] = ta2_yy_xxy_yy_0[i] * pa_z[i] - ta2_yy_xxy_yy_1[i] * pc_z[i];

        ta2_yy_xxyz_yz_0[i] = ta2_yy_yz_yz_0[i] * fe_0 - ta2_yy_yz_yz_1[i] * fe_0 + ta2_yy_xyz_yz_0[i] * pa_x[i] - ta2_yy_xyz_yz_1[i] * pc_x[i];

        ta2_yy_xxyz_zz_0[i] = 2.0 * ta1_y_xxz_zz_1[i] + ta2_yy_xxz_zz_0[i] * pa_y[i] - ta2_yy_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 300-306 components of targeted buffer : GD

    auto ta2_yy_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 300);

    auto ta2_yy_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 301);

    auto ta2_yy_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 302);

    auto ta2_yy_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 303);

    auto ta2_yy_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 304);

    auto ta2_yy_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 305);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta2_yy_xx_xx_0,   \
                             ta2_yy_xx_xx_1,   \
                             ta2_yy_xx_xy_0,   \
                             ta2_yy_xx_xy_1,   \
                             ta2_yy_xxz_xx_0,  \
                             ta2_yy_xxz_xx_1,  \
                             ta2_yy_xxz_xy_0,  \
                             ta2_yy_xxz_xy_1,  \
                             ta2_yy_xxzz_xx_0, \
                             ta2_yy_xxzz_xy_0, \
                             ta2_yy_xxzz_xz_0, \
                             ta2_yy_xxzz_yy_0, \
                             ta2_yy_xxzz_yz_0, \
                             ta2_yy_xxzz_zz_0, \
                             ta2_yy_xzz_xz_0,  \
                             ta2_yy_xzz_xz_1,  \
                             ta2_yy_xzz_yy_0,  \
                             ta2_yy_xzz_yy_1,  \
                             ta2_yy_xzz_yz_0,  \
                             ta2_yy_xzz_yz_1,  \
                             ta2_yy_xzz_z_0,   \
                             ta2_yy_xzz_z_1,   \
                             ta2_yy_xzz_zz_0,  \
                             ta2_yy_xzz_zz_1,  \
                             ta2_yy_zz_xz_0,   \
                             ta2_yy_zz_xz_1,   \
                             ta2_yy_zz_yy_0,   \
                             ta2_yy_zz_yy_1,   \
                             ta2_yy_zz_yz_0,   \
                             ta2_yy_zz_yz_1,   \
                             ta2_yy_zz_zz_0,   \
                             ta2_yy_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxzz_xx_0[i] = ta2_yy_xx_xx_0[i] * fe_0 - ta2_yy_xx_xx_1[i] * fe_0 + ta2_yy_xxz_xx_0[i] * pa_z[i] - ta2_yy_xxz_xx_1[i] * pc_z[i];

        ta2_yy_xxzz_xy_0[i] = ta2_yy_xx_xy_0[i] * fe_0 - ta2_yy_xx_xy_1[i] * fe_0 + ta2_yy_xxz_xy_0[i] * pa_z[i] - ta2_yy_xxz_xy_1[i] * pc_z[i];

        ta2_yy_xxzz_xz_0[i] = ta2_yy_zz_xz_0[i] * fe_0 - ta2_yy_zz_xz_1[i] * fe_0 + ta2_yy_xzz_z_0[i] * fe_0 - ta2_yy_xzz_z_1[i] * fe_0 +
                              ta2_yy_xzz_xz_0[i] * pa_x[i] - ta2_yy_xzz_xz_1[i] * pc_x[i];

        ta2_yy_xxzz_yy_0[i] = ta2_yy_zz_yy_0[i] * fe_0 - ta2_yy_zz_yy_1[i] * fe_0 + ta2_yy_xzz_yy_0[i] * pa_x[i] - ta2_yy_xzz_yy_1[i] * pc_x[i];

        ta2_yy_xxzz_yz_0[i] = ta2_yy_zz_yz_0[i] * fe_0 - ta2_yy_zz_yz_1[i] * fe_0 + ta2_yy_xzz_yz_0[i] * pa_x[i] - ta2_yy_xzz_yz_1[i] * pc_x[i];

        ta2_yy_xxzz_zz_0[i] = ta2_yy_zz_zz_0[i] * fe_0 - ta2_yy_zz_zz_1[i] * fe_0 + ta2_yy_xzz_zz_0[i] * pa_x[i] - ta2_yy_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 306-312 components of targeted buffer : GD

    auto ta2_yy_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 306);

    auto ta2_yy_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 307);

    auto ta2_yy_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 308);

    auto ta2_yy_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 309);

    auto ta2_yy_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 310);

    auto ta2_yy_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 311);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yy_xyyy_xx_0, \
                             ta2_yy_xyyy_xy_0, \
                             ta2_yy_xyyy_xz_0, \
                             ta2_yy_xyyy_yy_0, \
                             ta2_yy_xyyy_yz_0, \
                             ta2_yy_xyyy_zz_0, \
                             ta2_yy_yyy_x_0,   \
                             ta2_yy_yyy_x_1,   \
                             ta2_yy_yyy_xx_0,  \
                             ta2_yy_yyy_xx_1,  \
                             ta2_yy_yyy_xy_0,  \
                             ta2_yy_yyy_xy_1,  \
                             ta2_yy_yyy_xz_0,  \
                             ta2_yy_yyy_xz_1,  \
                             ta2_yy_yyy_y_0,   \
                             ta2_yy_yyy_y_1,   \
                             ta2_yy_yyy_yy_0,  \
                             ta2_yy_yyy_yy_1,  \
                             ta2_yy_yyy_yz_0,  \
                             ta2_yy_yyy_yz_1,  \
                             ta2_yy_yyy_z_0,   \
                             ta2_yy_yyy_z_1,   \
                             ta2_yy_yyy_zz_0,  \
                             ta2_yy_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyy_xx_0[i] =
            2.0 * ta2_yy_yyy_x_0[i] * fe_0 - 2.0 * ta2_yy_yyy_x_1[i] * fe_0 + ta2_yy_yyy_xx_0[i] * pa_x[i] - ta2_yy_yyy_xx_1[i] * pc_x[i];

        ta2_yy_xyyy_xy_0[i] = ta2_yy_yyy_y_0[i] * fe_0 - ta2_yy_yyy_y_1[i] * fe_0 + ta2_yy_yyy_xy_0[i] * pa_x[i] - ta2_yy_yyy_xy_1[i] * pc_x[i];

        ta2_yy_xyyy_xz_0[i] = ta2_yy_yyy_z_0[i] * fe_0 - ta2_yy_yyy_z_1[i] * fe_0 + ta2_yy_yyy_xz_0[i] * pa_x[i] - ta2_yy_yyy_xz_1[i] * pc_x[i];

        ta2_yy_xyyy_yy_0[i] = ta2_yy_yyy_yy_0[i] * pa_x[i] - ta2_yy_yyy_yy_1[i] * pc_x[i];

        ta2_yy_xyyy_yz_0[i] = ta2_yy_yyy_yz_0[i] * pa_x[i] - ta2_yy_yyy_yz_1[i] * pc_x[i];

        ta2_yy_xyyy_zz_0[i] = ta2_yy_yyy_zz_0[i] * pa_x[i] - ta2_yy_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 312-318 components of targeted buffer : GD

    auto ta2_yy_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 312);

    auto ta2_yy_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 313);

    auto ta2_yy_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 314);

    auto ta2_yy_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 315);

    auto ta2_yy_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 316);

    auto ta2_yy_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 317);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta2_yy_xyy_xx_0,  \
                             ta2_yy_xyy_xx_1,  \
                             ta2_yy_xyy_xy_0,  \
                             ta2_yy_xyy_xy_1,  \
                             ta2_yy_xyyz_xx_0, \
                             ta2_yy_xyyz_xy_0, \
                             ta2_yy_xyyz_xz_0, \
                             ta2_yy_xyyz_yy_0, \
                             ta2_yy_xyyz_yz_0, \
                             ta2_yy_xyyz_zz_0, \
                             ta2_yy_yyz_xz_0,  \
                             ta2_yy_yyz_xz_1,  \
                             ta2_yy_yyz_yy_0,  \
                             ta2_yy_yyz_yy_1,  \
                             ta2_yy_yyz_yz_0,  \
                             ta2_yy_yyz_yz_1,  \
                             ta2_yy_yyz_z_0,   \
                             ta2_yy_yyz_z_1,   \
                             ta2_yy_yyz_zz_0,  \
                             ta2_yy_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyz_xx_0[i] = ta2_yy_xyy_xx_0[i] * pa_z[i] - ta2_yy_xyy_xx_1[i] * pc_z[i];

        ta2_yy_xyyz_xy_0[i] = ta2_yy_xyy_xy_0[i] * pa_z[i] - ta2_yy_xyy_xy_1[i] * pc_z[i];

        ta2_yy_xyyz_xz_0[i] = ta2_yy_yyz_z_0[i] * fe_0 - ta2_yy_yyz_z_1[i] * fe_0 + ta2_yy_yyz_xz_0[i] * pa_x[i] - ta2_yy_yyz_xz_1[i] * pc_x[i];

        ta2_yy_xyyz_yy_0[i] = ta2_yy_yyz_yy_0[i] * pa_x[i] - ta2_yy_yyz_yy_1[i] * pc_x[i];

        ta2_yy_xyyz_yz_0[i] = ta2_yy_yyz_yz_0[i] * pa_x[i] - ta2_yy_yyz_yz_1[i] * pc_x[i];

        ta2_yy_xyyz_zz_0[i] = ta2_yy_yyz_zz_0[i] * pa_x[i] - ta2_yy_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 318-324 components of targeted buffer : GD

    auto ta2_yy_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 318);

    auto ta2_yy_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 319);

    auto ta2_yy_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 320);

    auto ta2_yy_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 321);

    auto ta2_yy_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 322);

    auto ta2_yy_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 323);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xzz_xx_1,   \
                             ta1_y_xzz_xz_1,   \
                             ta2_yy_xyzz_xx_0, \
                             ta2_yy_xyzz_xy_0, \
                             ta2_yy_xyzz_xz_0, \
                             ta2_yy_xyzz_yy_0, \
                             ta2_yy_xyzz_yz_0, \
                             ta2_yy_xyzz_zz_0, \
                             ta2_yy_xzz_xx_0,  \
                             ta2_yy_xzz_xx_1,  \
                             ta2_yy_xzz_xz_0,  \
                             ta2_yy_xzz_xz_1,  \
                             ta2_yy_yzz_xy_0,  \
                             ta2_yy_yzz_xy_1,  \
                             ta2_yy_yzz_y_0,   \
                             ta2_yy_yzz_y_1,   \
                             ta2_yy_yzz_yy_0,  \
                             ta2_yy_yzz_yy_1,  \
                             ta2_yy_yzz_yz_0,  \
                             ta2_yy_yzz_yz_1,  \
                             ta2_yy_yzz_zz_0,  \
                             ta2_yy_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyzz_xx_0[i] = 2.0 * ta1_y_xzz_xx_1[i] + ta2_yy_xzz_xx_0[i] * pa_y[i] - ta2_yy_xzz_xx_1[i] * pc_y[i];

        ta2_yy_xyzz_xy_0[i] = ta2_yy_yzz_y_0[i] * fe_0 - ta2_yy_yzz_y_1[i] * fe_0 + ta2_yy_yzz_xy_0[i] * pa_x[i] - ta2_yy_yzz_xy_1[i] * pc_x[i];

        ta2_yy_xyzz_xz_0[i] = 2.0 * ta1_y_xzz_xz_1[i] + ta2_yy_xzz_xz_0[i] * pa_y[i] - ta2_yy_xzz_xz_1[i] * pc_y[i];

        ta2_yy_xyzz_yy_0[i] = ta2_yy_yzz_yy_0[i] * pa_x[i] - ta2_yy_yzz_yy_1[i] * pc_x[i];

        ta2_yy_xyzz_yz_0[i] = ta2_yy_yzz_yz_0[i] * pa_x[i] - ta2_yy_yzz_yz_1[i] * pc_x[i];

        ta2_yy_xyzz_zz_0[i] = ta2_yy_yzz_zz_0[i] * pa_x[i] - ta2_yy_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 324-330 components of targeted buffer : GD

    auto ta2_yy_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 324);

    auto ta2_yy_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 325);

    auto ta2_yy_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 326);

    auto ta2_yy_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 327);

    auto ta2_yy_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 328);

    auto ta2_yy_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 329);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yy_xzzz_xx_0, \
                             ta2_yy_xzzz_xy_0, \
                             ta2_yy_xzzz_xz_0, \
                             ta2_yy_xzzz_yy_0, \
                             ta2_yy_xzzz_yz_0, \
                             ta2_yy_xzzz_zz_0, \
                             ta2_yy_zzz_x_0,   \
                             ta2_yy_zzz_x_1,   \
                             ta2_yy_zzz_xx_0,  \
                             ta2_yy_zzz_xx_1,  \
                             ta2_yy_zzz_xy_0,  \
                             ta2_yy_zzz_xy_1,  \
                             ta2_yy_zzz_xz_0,  \
                             ta2_yy_zzz_xz_1,  \
                             ta2_yy_zzz_y_0,   \
                             ta2_yy_zzz_y_1,   \
                             ta2_yy_zzz_yy_0,  \
                             ta2_yy_zzz_yy_1,  \
                             ta2_yy_zzz_yz_0,  \
                             ta2_yy_zzz_yz_1,  \
                             ta2_yy_zzz_z_0,   \
                             ta2_yy_zzz_z_1,   \
                             ta2_yy_zzz_zz_0,  \
                             ta2_yy_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzzz_xx_0[i] =
            2.0 * ta2_yy_zzz_x_0[i] * fe_0 - 2.0 * ta2_yy_zzz_x_1[i] * fe_0 + ta2_yy_zzz_xx_0[i] * pa_x[i] - ta2_yy_zzz_xx_1[i] * pc_x[i];

        ta2_yy_xzzz_xy_0[i] = ta2_yy_zzz_y_0[i] * fe_0 - ta2_yy_zzz_y_1[i] * fe_0 + ta2_yy_zzz_xy_0[i] * pa_x[i] - ta2_yy_zzz_xy_1[i] * pc_x[i];

        ta2_yy_xzzz_xz_0[i] = ta2_yy_zzz_z_0[i] * fe_0 - ta2_yy_zzz_z_1[i] * fe_0 + ta2_yy_zzz_xz_0[i] * pa_x[i] - ta2_yy_zzz_xz_1[i] * pc_x[i];

        ta2_yy_xzzz_yy_0[i] = ta2_yy_zzz_yy_0[i] * pa_x[i] - ta2_yy_zzz_yy_1[i] * pc_x[i];

        ta2_yy_xzzz_yz_0[i] = ta2_yy_zzz_yz_0[i] * pa_x[i] - ta2_yy_zzz_yz_1[i] * pc_x[i];

        ta2_yy_xzzz_zz_0[i] = ta2_yy_zzz_zz_0[i] * pa_x[i] - ta2_yy_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 330-336 components of targeted buffer : GD

    auto ta2_yy_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 330);

    auto ta2_yy_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 331);

    auto ta2_yy_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 332);

    auto ta2_yy_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 333);

    auto ta2_yy_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 334);

    auto ta2_yy_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 335);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_zz_1,   \
                             ta2_yy_yy_xx_0,   \
                             ta2_yy_yy_xx_1,   \
                             ta2_yy_yy_xy_0,   \
                             ta2_yy_yy_xy_1,   \
                             ta2_yy_yy_xz_0,   \
                             ta2_yy_yy_xz_1,   \
                             ta2_yy_yy_yy_0,   \
                             ta2_yy_yy_yy_1,   \
                             ta2_yy_yy_yz_0,   \
                             ta2_yy_yy_yz_1,   \
                             ta2_yy_yy_zz_0,   \
                             ta2_yy_yy_zz_1,   \
                             ta2_yy_yyy_x_0,   \
                             ta2_yy_yyy_x_1,   \
                             ta2_yy_yyy_xx_0,  \
                             ta2_yy_yyy_xx_1,  \
                             ta2_yy_yyy_xy_0,  \
                             ta2_yy_yyy_xy_1,  \
                             ta2_yy_yyy_xz_0,  \
                             ta2_yy_yyy_xz_1,  \
                             ta2_yy_yyy_y_0,   \
                             ta2_yy_yyy_y_1,   \
                             ta2_yy_yyy_yy_0,  \
                             ta2_yy_yyy_yy_1,  \
                             ta2_yy_yyy_yz_0,  \
                             ta2_yy_yyy_yz_1,  \
                             ta2_yy_yyy_z_0,   \
                             ta2_yy_yyy_z_1,   \
                             ta2_yy_yyy_zz_0,  \
                             ta2_yy_yyy_zz_1,  \
                             ta2_yy_yyyy_xx_0, \
                             ta2_yy_yyyy_xy_0, \
                             ta2_yy_yyyy_xz_0, \
                             ta2_yy_yyyy_yy_0, \
                             ta2_yy_yyyy_yz_0, \
                             ta2_yy_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyy_xx_0[i] = 3.0 * ta2_yy_yy_xx_0[i] * fe_0 - 3.0 * ta2_yy_yy_xx_1[i] * fe_0 + 2.0 * ta1_y_yyy_xx_1[i] +
                              ta2_yy_yyy_xx_0[i] * pa_y[i] - ta2_yy_yyy_xx_1[i] * pc_y[i];

        ta2_yy_yyyy_xy_0[i] = 3.0 * ta2_yy_yy_xy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xy_1[i] * fe_0 + ta2_yy_yyy_x_0[i] * fe_0 - ta2_yy_yyy_x_1[i] * fe_0 +
                              2.0 * ta1_y_yyy_xy_1[i] + ta2_yy_yyy_xy_0[i] * pa_y[i] - ta2_yy_yyy_xy_1[i] * pc_y[i];

        ta2_yy_yyyy_xz_0[i] = 3.0 * ta2_yy_yy_xz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xz_1[i] +
                              ta2_yy_yyy_xz_0[i] * pa_y[i] - ta2_yy_yyy_xz_1[i] * pc_y[i];

        ta2_yy_yyyy_yy_0[i] = 3.0 * ta2_yy_yy_yy_0[i] * fe_0 - 3.0 * ta2_yy_yy_yy_1[i] * fe_0 + 2.0 * ta2_yy_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_yy_yyy_y_1[i] * fe_0 + 2.0 * ta1_y_yyy_yy_1[i] + ta2_yy_yyy_yy_0[i] * pa_y[i] - ta2_yy_yyy_yy_1[i] * pc_y[i];

        ta2_yy_yyyy_yz_0[i] = 3.0 * ta2_yy_yy_yz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yz_1[i] * fe_0 + ta2_yy_yyy_z_0[i] * fe_0 - ta2_yy_yyy_z_1[i] * fe_0 +
                              2.0 * ta1_y_yyy_yz_1[i] + ta2_yy_yyy_yz_0[i] * pa_y[i] - ta2_yy_yyy_yz_1[i] * pc_y[i];

        ta2_yy_yyyy_zz_0[i] = 3.0 * ta2_yy_yy_zz_0[i] * fe_0 - 3.0 * ta2_yy_yy_zz_1[i] * fe_0 + 2.0 * ta1_y_yyy_zz_1[i] +
                              ta2_yy_yyy_zz_0[i] * pa_y[i] - ta2_yy_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 336-342 components of targeted buffer : GD

    auto ta2_yy_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 336);

    auto ta2_yy_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 337);

    auto ta2_yy_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 338);

    auto ta2_yy_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 339);

    auto ta2_yy_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 340);

    auto ta2_yy_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 341);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_yy_yyy_x_0,   \
                             ta2_yy_yyy_x_1,   \
                             ta2_yy_yyy_xx_0,  \
                             ta2_yy_yyy_xx_1,  \
                             ta2_yy_yyy_xy_0,  \
                             ta2_yy_yyy_xy_1,  \
                             ta2_yy_yyy_xz_0,  \
                             ta2_yy_yyy_xz_1,  \
                             ta2_yy_yyy_y_0,   \
                             ta2_yy_yyy_y_1,   \
                             ta2_yy_yyy_yy_0,  \
                             ta2_yy_yyy_yy_1,  \
                             ta2_yy_yyy_yz_0,  \
                             ta2_yy_yyy_yz_1,  \
                             ta2_yy_yyy_z_0,   \
                             ta2_yy_yyy_z_1,   \
                             ta2_yy_yyy_zz_0,  \
                             ta2_yy_yyy_zz_1,  \
                             ta2_yy_yyyz_xx_0, \
                             ta2_yy_yyyz_xy_0, \
                             ta2_yy_yyyz_xz_0, \
                             ta2_yy_yyyz_yy_0, \
                             ta2_yy_yyyz_yz_0, \
                             ta2_yy_yyyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyz_xx_0[i] = ta2_yy_yyy_xx_0[i] * pa_z[i] - ta2_yy_yyy_xx_1[i] * pc_z[i];

        ta2_yy_yyyz_xy_0[i] = ta2_yy_yyy_xy_0[i] * pa_z[i] - ta2_yy_yyy_xy_1[i] * pc_z[i];

        ta2_yy_yyyz_xz_0[i] = ta2_yy_yyy_x_0[i] * fe_0 - ta2_yy_yyy_x_1[i] * fe_0 + ta2_yy_yyy_xz_0[i] * pa_z[i] - ta2_yy_yyy_xz_1[i] * pc_z[i];

        ta2_yy_yyyz_yy_0[i] = ta2_yy_yyy_yy_0[i] * pa_z[i] - ta2_yy_yyy_yy_1[i] * pc_z[i];

        ta2_yy_yyyz_yz_0[i] = ta2_yy_yyy_y_0[i] * fe_0 - ta2_yy_yyy_y_1[i] * fe_0 + ta2_yy_yyy_yz_0[i] * pa_z[i] - ta2_yy_yyy_yz_1[i] * pc_z[i];

        ta2_yy_yyyz_zz_0[i] =
            2.0 * ta2_yy_yyy_z_0[i] * fe_0 - 2.0 * ta2_yy_yyy_z_1[i] * fe_0 + ta2_yy_yyy_zz_0[i] * pa_z[i] - ta2_yy_yyy_zz_1[i] * pc_z[i];
    }

    // Set up 342-348 components of targeted buffer : GD

    auto ta2_yy_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 342);

    auto ta2_yy_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 343);

    auto ta2_yy_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 344);

    auto ta2_yy_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 345);

    auto ta2_yy_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 346);

    auto ta2_yy_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 347);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yzz_xz_1,   \
                             ta1_y_yzz_zz_1,   \
                             ta2_yy_yy_xx_0,   \
                             ta2_yy_yy_xx_1,   \
                             ta2_yy_yy_xy_0,   \
                             ta2_yy_yy_xy_1,   \
                             ta2_yy_yy_yy_0,   \
                             ta2_yy_yy_yy_1,   \
                             ta2_yy_yy_yz_0,   \
                             ta2_yy_yy_yz_1,   \
                             ta2_yy_yyz_xx_0,  \
                             ta2_yy_yyz_xx_1,  \
                             ta2_yy_yyz_xy_0,  \
                             ta2_yy_yyz_xy_1,  \
                             ta2_yy_yyz_y_0,   \
                             ta2_yy_yyz_y_1,   \
                             ta2_yy_yyz_yy_0,  \
                             ta2_yy_yyz_yy_1,  \
                             ta2_yy_yyz_yz_0,  \
                             ta2_yy_yyz_yz_1,  \
                             ta2_yy_yyzz_xx_0, \
                             ta2_yy_yyzz_xy_0, \
                             ta2_yy_yyzz_xz_0, \
                             ta2_yy_yyzz_yy_0, \
                             ta2_yy_yyzz_yz_0, \
                             ta2_yy_yyzz_zz_0, \
                             ta2_yy_yzz_xz_0,  \
                             ta2_yy_yzz_xz_1,  \
                             ta2_yy_yzz_zz_0,  \
                             ta2_yy_yzz_zz_1,  \
                             ta2_yy_zz_xz_0,   \
                             ta2_yy_zz_xz_1,   \
                             ta2_yy_zz_zz_0,   \
                             ta2_yy_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyzz_xx_0[i] = ta2_yy_yy_xx_0[i] * fe_0 - ta2_yy_yy_xx_1[i] * fe_0 + ta2_yy_yyz_xx_0[i] * pa_z[i] - ta2_yy_yyz_xx_1[i] * pc_z[i];

        ta2_yy_yyzz_xy_0[i] = ta2_yy_yy_xy_0[i] * fe_0 - ta2_yy_yy_xy_1[i] * fe_0 + ta2_yy_yyz_xy_0[i] * pa_z[i] - ta2_yy_yyz_xy_1[i] * pc_z[i];

        ta2_yy_yyzz_xz_0[i] = ta2_yy_zz_xz_0[i] * fe_0 - ta2_yy_zz_xz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xz_1[i] + ta2_yy_yzz_xz_0[i] * pa_y[i] -
                              ta2_yy_yzz_xz_1[i] * pc_y[i];

        ta2_yy_yyzz_yy_0[i] = ta2_yy_yy_yy_0[i] * fe_0 - ta2_yy_yy_yy_1[i] * fe_0 + ta2_yy_yyz_yy_0[i] * pa_z[i] - ta2_yy_yyz_yy_1[i] * pc_z[i];

        ta2_yy_yyzz_yz_0[i] = ta2_yy_yy_yz_0[i] * fe_0 - ta2_yy_yy_yz_1[i] * fe_0 + ta2_yy_yyz_y_0[i] * fe_0 - ta2_yy_yyz_y_1[i] * fe_0 +
                              ta2_yy_yyz_yz_0[i] * pa_z[i] - ta2_yy_yyz_yz_1[i] * pc_z[i];

        ta2_yy_yyzz_zz_0[i] = ta2_yy_zz_zz_0[i] * fe_0 - ta2_yy_zz_zz_1[i] * fe_0 + 2.0 * ta1_y_yzz_zz_1[i] + ta2_yy_yzz_zz_0[i] * pa_y[i] -
                              ta2_yy_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 348-354 components of targeted buffer : GD

    auto ta2_yy_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 348);

    auto ta2_yy_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 349);

    auto ta2_yy_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 350);

    auto ta2_yy_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 351);

    auto ta2_yy_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 352);

    auto ta2_yy_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 353);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_zzz_xx_1,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_zz_1,   \
                             ta2_yy_yz_xy_0,   \
                             ta2_yy_yz_xy_1,   \
                             ta2_yy_yz_yy_0,   \
                             ta2_yy_yz_yy_1,   \
                             ta2_yy_yzz_xy_0,  \
                             ta2_yy_yzz_xy_1,  \
                             ta2_yy_yzz_yy_0,  \
                             ta2_yy_yzz_yy_1,  \
                             ta2_yy_yzzz_xx_0, \
                             ta2_yy_yzzz_xy_0, \
                             ta2_yy_yzzz_xz_0, \
                             ta2_yy_yzzz_yy_0, \
                             ta2_yy_yzzz_yz_0, \
                             ta2_yy_yzzz_zz_0, \
                             ta2_yy_zzz_xx_0,  \
                             ta2_yy_zzz_xx_1,  \
                             ta2_yy_zzz_xz_0,  \
                             ta2_yy_zzz_xz_1,  \
                             ta2_yy_zzz_yz_0,  \
                             ta2_yy_zzz_yz_1,  \
                             ta2_yy_zzz_z_0,   \
                             ta2_yy_zzz_z_1,   \
                             ta2_yy_zzz_zz_0,  \
                             ta2_yy_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzzz_xx_0[i] = 2.0 * ta1_y_zzz_xx_1[i] + ta2_yy_zzz_xx_0[i] * pa_y[i] - ta2_yy_zzz_xx_1[i] * pc_y[i];

        ta2_yy_yzzz_xy_0[i] =
            2.0 * ta2_yy_yz_xy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xy_1[i] * fe_0 + ta2_yy_yzz_xy_0[i] * pa_z[i] - ta2_yy_yzz_xy_1[i] * pc_z[i];

        ta2_yy_yzzz_xz_0[i] = 2.0 * ta1_y_zzz_xz_1[i] + ta2_yy_zzz_xz_0[i] * pa_y[i] - ta2_yy_zzz_xz_1[i] * pc_y[i];

        ta2_yy_yzzz_yy_0[i] =
            2.0 * ta2_yy_yz_yy_0[i] * fe_0 - 2.0 * ta2_yy_yz_yy_1[i] * fe_0 + ta2_yy_yzz_yy_0[i] * pa_z[i] - ta2_yy_yzz_yy_1[i] * pc_z[i];

        ta2_yy_yzzz_yz_0[i] = ta2_yy_zzz_z_0[i] * fe_0 - ta2_yy_zzz_z_1[i] * fe_0 + 2.0 * ta1_y_zzz_yz_1[i] + ta2_yy_zzz_yz_0[i] * pa_y[i] -
                              ta2_yy_zzz_yz_1[i] * pc_y[i];

        ta2_yy_yzzz_zz_0[i] = 2.0 * ta1_y_zzz_zz_1[i] + ta2_yy_zzz_zz_0[i] * pa_y[i] - ta2_yy_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 354-360 components of targeted buffer : GD

    auto ta2_yy_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 354);

    auto ta2_yy_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 355);

    auto ta2_yy_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 356);

    auto ta2_yy_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 357);

    auto ta2_yy_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 358);

    auto ta2_yy_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 359);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_yy_zz_xx_0,   \
                             ta2_yy_zz_xx_1,   \
                             ta2_yy_zz_xy_0,   \
                             ta2_yy_zz_xy_1,   \
                             ta2_yy_zz_xz_0,   \
                             ta2_yy_zz_xz_1,   \
                             ta2_yy_zz_yy_0,   \
                             ta2_yy_zz_yy_1,   \
                             ta2_yy_zz_yz_0,   \
                             ta2_yy_zz_yz_1,   \
                             ta2_yy_zz_zz_0,   \
                             ta2_yy_zz_zz_1,   \
                             ta2_yy_zzz_x_0,   \
                             ta2_yy_zzz_x_1,   \
                             ta2_yy_zzz_xx_0,  \
                             ta2_yy_zzz_xx_1,  \
                             ta2_yy_zzz_xy_0,  \
                             ta2_yy_zzz_xy_1,  \
                             ta2_yy_zzz_xz_0,  \
                             ta2_yy_zzz_xz_1,  \
                             ta2_yy_zzz_y_0,   \
                             ta2_yy_zzz_y_1,   \
                             ta2_yy_zzz_yy_0,  \
                             ta2_yy_zzz_yy_1,  \
                             ta2_yy_zzz_yz_0,  \
                             ta2_yy_zzz_yz_1,  \
                             ta2_yy_zzz_z_0,   \
                             ta2_yy_zzz_z_1,   \
                             ta2_yy_zzz_zz_0,  \
                             ta2_yy_zzz_zz_1,  \
                             ta2_yy_zzzz_xx_0, \
                             ta2_yy_zzzz_xy_0, \
                             ta2_yy_zzzz_xz_0, \
                             ta2_yy_zzzz_yy_0, \
                             ta2_yy_zzzz_yz_0, \
                             ta2_yy_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzzz_xx_0[i] =
            3.0 * ta2_yy_zz_xx_0[i] * fe_0 - 3.0 * ta2_yy_zz_xx_1[i] * fe_0 + ta2_yy_zzz_xx_0[i] * pa_z[i] - ta2_yy_zzz_xx_1[i] * pc_z[i];

        ta2_yy_zzzz_xy_0[i] =
            3.0 * ta2_yy_zz_xy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xy_1[i] * fe_0 + ta2_yy_zzz_xy_0[i] * pa_z[i] - ta2_yy_zzz_xy_1[i] * pc_z[i];

        ta2_yy_zzzz_xz_0[i] = 3.0 * ta2_yy_zz_xz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xz_1[i] * fe_0 + ta2_yy_zzz_x_0[i] * fe_0 - ta2_yy_zzz_x_1[i] * fe_0 +
                              ta2_yy_zzz_xz_0[i] * pa_z[i] - ta2_yy_zzz_xz_1[i] * pc_z[i];

        ta2_yy_zzzz_yy_0[i] =
            3.0 * ta2_yy_zz_yy_0[i] * fe_0 - 3.0 * ta2_yy_zz_yy_1[i] * fe_0 + ta2_yy_zzz_yy_0[i] * pa_z[i] - ta2_yy_zzz_yy_1[i] * pc_z[i];

        ta2_yy_zzzz_yz_0[i] = 3.0 * ta2_yy_zz_yz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yz_1[i] * fe_0 + ta2_yy_zzz_y_0[i] * fe_0 - ta2_yy_zzz_y_1[i] * fe_0 +
                              ta2_yy_zzz_yz_0[i] * pa_z[i] - ta2_yy_zzz_yz_1[i] * pc_z[i];

        ta2_yy_zzzz_zz_0[i] = 3.0 * ta2_yy_zz_zz_0[i] * fe_0 - 3.0 * ta2_yy_zz_zz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_yy_zzz_z_1[i] * fe_0 + ta2_yy_zzz_zz_0[i] * pa_z[i] - ta2_yy_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 360-366 components of targeted buffer : GD

    auto ta2_yz_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 360);

    auto ta2_yz_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 361);

    auto ta2_yz_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 362);

    auto ta2_yz_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 363);

    auto ta2_yz_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 364);

    auto ta2_yz_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 365);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yz_xx_xx_0,   \
                             ta2_yz_xx_xx_1,   \
                             ta2_yz_xx_xy_0,   \
                             ta2_yz_xx_xy_1,   \
                             ta2_yz_xx_xz_0,   \
                             ta2_yz_xx_xz_1,   \
                             ta2_yz_xx_yy_0,   \
                             ta2_yz_xx_yy_1,   \
                             ta2_yz_xx_yz_0,   \
                             ta2_yz_xx_yz_1,   \
                             ta2_yz_xx_zz_0,   \
                             ta2_yz_xx_zz_1,   \
                             ta2_yz_xxx_x_0,   \
                             ta2_yz_xxx_x_1,   \
                             ta2_yz_xxx_xx_0,  \
                             ta2_yz_xxx_xx_1,  \
                             ta2_yz_xxx_xy_0,  \
                             ta2_yz_xxx_xy_1,  \
                             ta2_yz_xxx_xz_0,  \
                             ta2_yz_xxx_xz_1,  \
                             ta2_yz_xxx_y_0,   \
                             ta2_yz_xxx_y_1,   \
                             ta2_yz_xxx_yy_0,  \
                             ta2_yz_xxx_yy_1,  \
                             ta2_yz_xxx_yz_0,  \
                             ta2_yz_xxx_yz_1,  \
                             ta2_yz_xxx_z_0,   \
                             ta2_yz_xxx_z_1,   \
                             ta2_yz_xxx_zz_0,  \
                             ta2_yz_xxx_zz_1,  \
                             ta2_yz_xxxx_xx_0, \
                             ta2_yz_xxxx_xy_0, \
                             ta2_yz_xxxx_xz_0, \
                             ta2_yz_xxxx_yy_0, \
                             ta2_yz_xxxx_yz_0, \
                             ta2_yz_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxx_xx_0[i] = 3.0 * ta2_yz_xx_xx_0[i] * fe_0 - 3.0 * ta2_yz_xx_xx_1[i] * fe_0 + 2.0 * ta2_yz_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_yz_xxx_x_1[i] * fe_0 + ta2_yz_xxx_xx_0[i] * pa_x[i] - ta2_yz_xxx_xx_1[i] * pc_x[i];

        ta2_yz_xxxx_xy_0[i] = 3.0 * ta2_yz_xx_xy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xy_1[i] * fe_0 + ta2_yz_xxx_y_0[i] * fe_0 - ta2_yz_xxx_y_1[i] * fe_0 +
                              ta2_yz_xxx_xy_0[i] * pa_x[i] - ta2_yz_xxx_xy_1[i] * pc_x[i];

        ta2_yz_xxxx_xz_0[i] = 3.0 * ta2_yz_xx_xz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xz_1[i] * fe_0 + ta2_yz_xxx_z_0[i] * fe_0 - ta2_yz_xxx_z_1[i] * fe_0 +
                              ta2_yz_xxx_xz_0[i] * pa_x[i] - ta2_yz_xxx_xz_1[i] * pc_x[i];

        ta2_yz_xxxx_yy_0[i] =
            3.0 * ta2_yz_xx_yy_0[i] * fe_0 - 3.0 * ta2_yz_xx_yy_1[i] * fe_0 + ta2_yz_xxx_yy_0[i] * pa_x[i] - ta2_yz_xxx_yy_1[i] * pc_x[i];

        ta2_yz_xxxx_yz_0[i] =
            3.0 * ta2_yz_xx_yz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yz_1[i] * fe_0 + ta2_yz_xxx_yz_0[i] * pa_x[i] - ta2_yz_xxx_yz_1[i] * pc_x[i];

        ta2_yz_xxxx_zz_0[i] =
            3.0 * ta2_yz_xx_zz_0[i] * fe_0 - 3.0 * ta2_yz_xx_zz_1[i] * fe_0 + ta2_yz_xxx_zz_0[i] * pa_x[i] - ta2_yz_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 366-372 components of targeted buffer : GD

    auto ta2_yz_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 366);

    auto ta2_yz_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 367);

    auto ta2_yz_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 368);

    auto ta2_yz_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 369);

    auto ta2_yz_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 370);

    auto ta2_yz_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 371);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_zz_1,   \
                             ta2_yz_xxx_x_0,   \
                             ta2_yz_xxx_x_1,   \
                             ta2_yz_xxx_xx_0,  \
                             ta2_yz_xxx_xx_1,  \
                             ta2_yz_xxx_xy_0,  \
                             ta2_yz_xxx_xy_1,  \
                             ta2_yz_xxx_xz_0,  \
                             ta2_yz_xxx_xz_1,  \
                             ta2_yz_xxx_zz_0,  \
                             ta2_yz_xxx_zz_1,  \
                             ta2_yz_xxxy_xx_0, \
                             ta2_yz_xxxy_xy_0, \
                             ta2_yz_xxxy_xz_0, \
                             ta2_yz_xxxy_yy_0, \
                             ta2_yz_xxxy_yz_0, \
                             ta2_yz_xxxy_zz_0, \
                             ta2_yz_xxy_yy_0,  \
                             ta2_yz_xxy_yy_1,  \
                             ta2_yz_xxy_yz_0,  \
                             ta2_yz_xxy_yz_1,  \
                             ta2_yz_xy_yy_0,   \
                             ta2_yz_xy_yy_1,   \
                             ta2_yz_xy_yz_0,   \
                             ta2_yz_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxy_xx_0[i] = ta1_z_xxx_xx_1[i] + ta2_yz_xxx_xx_0[i] * pa_y[i] - ta2_yz_xxx_xx_1[i] * pc_y[i];

        ta2_yz_xxxy_xy_0[i] =
            ta2_yz_xxx_x_0[i] * fe_0 - ta2_yz_xxx_x_1[i] * fe_0 + ta1_z_xxx_xy_1[i] + ta2_yz_xxx_xy_0[i] * pa_y[i] - ta2_yz_xxx_xy_1[i] * pc_y[i];

        ta2_yz_xxxy_xz_0[i] = ta1_z_xxx_xz_1[i] + ta2_yz_xxx_xz_0[i] * pa_y[i] - ta2_yz_xxx_xz_1[i] * pc_y[i];

        ta2_yz_xxxy_yy_0[i] =
            2.0 * ta2_yz_xy_yy_0[i] * fe_0 - 2.0 * ta2_yz_xy_yy_1[i] * fe_0 + ta2_yz_xxy_yy_0[i] * pa_x[i] - ta2_yz_xxy_yy_1[i] * pc_x[i];

        ta2_yz_xxxy_yz_0[i] =
            2.0 * ta2_yz_xy_yz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yz_1[i] * fe_0 + ta2_yz_xxy_yz_0[i] * pa_x[i] - ta2_yz_xxy_yz_1[i] * pc_x[i];

        ta2_yz_xxxy_zz_0[i] = ta1_z_xxx_zz_1[i] + ta2_yz_xxx_zz_0[i] * pa_y[i] - ta2_yz_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 372-378 components of targeted buffer : GD

    auto ta2_yz_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 372);

    auto ta2_yz_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 373);

    auto ta2_yz_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 374);

    auto ta2_yz_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 375);

    auto ta2_yz_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 376);

    auto ta2_yz_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 377);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_yy_1,   \
                             ta2_yz_xxx_x_0,   \
                             ta2_yz_xxx_x_1,   \
                             ta2_yz_xxx_xx_0,  \
                             ta2_yz_xxx_xx_1,  \
                             ta2_yz_xxx_xy_0,  \
                             ta2_yz_xxx_xy_1,  \
                             ta2_yz_xxx_xz_0,  \
                             ta2_yz_xxx_xz_1,  \
                             ta2_yz_xxx_yy_0,  \
                             ta2_yz_xxx_yy_1,  \
                             ta2_yz_xxxz_xx_0, \
                             ta2_yz_xxxz_xy_0, \
                             ta2_yz_xxxz_xz_0, \
                             ta2_yz_xxxz_yy_0, \
                             ta2_yz_xxxz_yz_0, \
                             ta2_yz_xxxz_zz_0, \
                             ta2_yz_xxz_yz_0,  \
                             ta2_yz_xxz_yz_1,  \
                             ta2_yz_xxz_zz_0,  \
                             ta2_yz_xxz_zz_1,  \
                             ta2_yz_xz_yz_0,   \
                             ta2_yz_xz_yz_1,   \
                             ta2_yz_xz_zz_0,   \
                             ta2_yz_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxz_xx_0[i] = ta1_y_xxx_xx_1[i] + ta2_yz_xxx_xx_0[i] * pa_z[i] - ta2_yz_xxx_xx_1[i] * pc_z[i];

        ta2_yz_xxxz_xy_0[i] = ta1_y_xxx_xy_1[i] + ta2_yz_xxx_xy_0[i] * pa_z[i] - ta2_yz_xxx_xy_1[i] * pc_z[i];

        ta2_yz_xxxz_xz_0[i] =
            ta2_yz_xxx_x_0[i] * fe_0 - ta2_yz_xxx_x_1[i] * fe_0 + ta1_y_xxx_xz_1[i] + ta2_yz_xxx_xz_0[i] * pa_z[i] - ta2_yz_xxx_xz_1[i] * pc_z[i];

        ta2_yz_xxxz_yy_0[i] = ta1_y_xxx_yy_1[i] + ta2_yz_xxx_yy_0[i] * pa_z[i] - ta2_yz_xxx_yy_1[i] * pc_z[i];

        ta2_yz_xxxz_yz_0[i] =
            2.0 * ta2_yz_xz_yz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yz_1[i] * fe_0 + ta2_yz_xxz_yz_0[i] * pa_x[i] - ta2_yz_xxz_yz_1[i] * pc_x[i];

        ta2_yz_xxxz_zz_0[i] =
            2.0 * ta2_yz_xz_zz_0[i] * fe_0 - 2.0 * ta2_yz_xz_zz_1[i] * fe_0 + ta2_yz_xxz_zz_0[i] * pa_x[i] - ta2_yz_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 378-384 components of targeted buffer : GD

    auto ta2_yz_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 378);

    auto ta2_yz_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 379);

    auto ta2_yz_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 380);

    auto ta2_yz_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 381);

    auto ta2_yz_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 382);

    auto ta2_yz_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 383);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxy_xx_1,   \
                             ta1_z_xxy_xz_1,   \
                             ta2_yz_xx_xx_0,   \
                             ta2_yz_xx_xx_1,   \
                             ta2_yz_xx_xz_0,   \
                             ta2_yz_xx_xz_1,   \
                             ta2_yz_xxy_xx_0,  \
                             ta2_yz_xxy_xx_1,  \
                             ta2_yz_xxy_xz_0,  \
                             ta2_yz_xxy_xz_1,  \
                             ta2_yz_xxyy_xx_0, \
                             ta2_yz_xxyy_xy_0, \
                             ta2_yz_xxyy_xz_0, \
                             ta2_yz_xxyy_yy_0, \
                             ta2_yz_xxyy_yz_0, \
                             ta2_yz_xxyy_zz_0, \
                             ta2_yz_xyy_xy_0,  \
                             ta2_yz_xyy_xy_1,  \
                             ta2_yz_xyy_y_0,   \
                             ta2_yz_xyy_y_1,   \
                             ta2_yz_xyy_yy_0,  \
                             ta2_yz_xyy_yy_1,  \
                             ta2_yz_xyy_yz_0,  \
                             ta2_yz_xyy_yz_1,  \
                             ta2_yz_xyy_zz_0,  \
                             ta2_yz_xyy_zz_1,  \
                             ta2_yz_yy_xy_0,   \
                             ta2_yz_yy_xy_1,   \
                             ta2_yz_yy_yy_0,   \
                             ta2_yz_yy_yy_1,   \
                             ta2_yz_yy_yz_0,   \
                             ta2_yz_yy_yz_1,   \
                             ta2_yz_yy_zz_0,   \
                             ta2_yz_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyy_xx_0[i] =
            ta2_yz_xx_xx_0[i] * fe_0 - ta2_yz_xx_xx_1[i] * fe_0 + ta1_z_xxy_xx_1[i] + ta2_yz_xxy_xx_0[i] * pa_y[i] - ta2_yz_xxy_xx_1[i] * pc_y[i];

        ta2_yz_xxyy_xy_0[i] = ta2_yz_yy_xy_0[i] * fe_0 - ta2_yz_yy_xy_1[i] * fe_0 + ta2_yz_xyy_y_0[i] * fe_0 - ta2_yz_xyy_y_1[i] * fe_0 +
                              ta2_yz_xyy_xy_0[i] * pa_x[i] - ta2_yz_xyy_xy_1[i] * pc_x[i];

        ta2_yz_xxyy_xz_0[i] =
            ta2_yz_xx_xz_0[i] * fe_0 - ta2_yz_xx_xz_1[i] * fe_0 + ta1_z_xxy_xz_1[i] + ta2_yz_xxy_xz_0[i] * pa_y[i] - ta2_yz_xxy_xz_1[i] * pc_y[i];

        ta2_yz_xxyy_yy_0[i] = ta2_yz_yy_yy_0[i] * fe_0 - ta2_yz_yy_yy_1[i] * fe_0 + ta2_yz_xyy_yy_0[i] * pa_x[i] - ta2_yz_xyy_yy_1[i] * pc_x[i];

        ta2_yz_xxyy_yz_0[i] = ta2_yz_yy_yz_0[i] * fe_0 - ta2_yz_yy_yz_1[i] * fe_0 + ta2_yz_xyy_yz_0[i] * pa_x[i] - ta2_yz_xyy_yz_1[i] * pc_x[i];

        ta2_yz_xxyy_zz_0[i] = ta2_yz_yy_zz_0[i] * fe_0 - ta2_yz_yy_zz_1[i] * fe_0 + ta2_yz_xyy_zz_0[i] * pa_x[i] - ta2_yz_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 384-390 components of targeted buffer : GD

    auto ta2_yz_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 384);

    auto ta2_yz_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 385);

    auto ta2_yz_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 386);

    auto ta2_yz_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 387);

    auto ta2_yz_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 388);

    auto ta2_yz_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 389);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xxy_xy_1,   \
                             ta1_y_xxy_yy_1,   \
                             ta1_z_xxz_xx_1,   \
                             ta1_z_xxz_xz_1,   \
                             ta1_z_xxz_zz_1,   \
                             ta2_yz_xxy_xy_0,  \
                             ta2_yz_xxy_xy_1,  \
                             ta2_yz_xxy_yy_0,  \
                             ta2_yz_xxy_yy_1,  \
                             ta2_yz_xxyz_xx_0, \
                             ta2_yz_xxyz_xy_0, \
                             ta2_yz_xxyz_xz_0, \
                             ta2_yz_xxyz_yy_0, \
                             ta2_yz_xxyz_yz_0, \
                             ta2_yz_xxyz_zz_0, \
                             ta2_yz_xxz_xx_0,  \
                             ta2_yz_xxz_xx_1,  \
                             ta2_yz_xxz_xz_0,  \
                             ta2_yz_xxz_xz_1,  \
                             ta2_yz_xxz_zz_0,  \
                             ta2_yz_xxz_zz_1,  \
                             ta2_yz_xyz_yz_0,  \
                             ta2_yz_xyz_yz_1,  \
                             ta2_yz_yz_yz_0,   \
                             ta2_yz_yz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyz_xx_0[i] = ta1_z_xxz_xx_1[i] + ta2_yz_xxz_xx_0[i] * pa_y[i] - ta2_yz_xxz_xx_1[i] * pc_y[i];

        ta2_yz_xxyz_xy_0[i] = ta1_y_xxy_xy_1[i] + ta2_yz_xxy_xy_0[i] * pa_z[i] - ta2_yz_xxy_xy_1[i] * pc_z[i];

        ta2_yz_xxyz_xz_0[i] = ta1_z_xxz_xz_1[i] + ta2_yz_xxz_xz_0[i] * pa_y[i] - ta2_yz_xxz_xz_1[i] * pc_y[i];

        ta2_yz_xxyz_yy_0[i] = ta1_y_xxy_yy_1[i] + ta2_yz_xxy_yy_0[i] * pa_z[i] - ta2_yz_xxy_yy_1[i] * pc_z[i];

        ta2_yz_xxyz_yz_0[i] = ta2_yz_yz_yz_0[i] * fe_0 - ta2_yz_yz_yz_1[i] * fe_0 + ta2_yz_xyz_yz_0[i] * pa_x[i] - ta2_yz_xyz_yz_1[i] * pc_x[i];

        ta2_yz_xxyz_zz_0[i] = ta1_z_xxz_zz_1[i] + ta2_yz_xxz_zz_0[i] * pa_y[i] - ta2_yz_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 390-396 components of targeted buffer : GD

    auto ta2_yz_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 390);

    auto ta2_yz_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 391);

    auto ta2_yz_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 392);

    auto ta2_yz_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 393);

    auto ta2_yz_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 394);

    auto ta2_yz_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 395);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxz_xx_1,   \
                             ta1_y_xxz_xy_1,   \
                             ta2_yz_xx_xx_0,   \
                             ta2_yz_xx_xx_1,   \
                             ta2_yz_xx_xy_0,   \
                             ta2_yz_xx_xy_1,   \
                             ta2_yz_xxz_xx_0,  \
                             ta2_yz_xxz_xx_1,  \
                             ta2_yz_xxz_xy_0,  \
                             ta2_yz_xxz_xy_1,  \
                             ta2_yz_xxzz_xx_0, \
                             ta2_yz_xxzz_xy_0, \
                             ta2_yz_xxzz_xz_0, \
                             ta2_yz_xxzz_yy_0, \
                             ta2_yz_xxzz_yz_0, \
                             ta2_yz_xxzz_zz_0, \
                             ta2_yz_xzz_xz_0,  \
                             ta2_yz_xzz_xz_1,  \
                             ta2_yz_xzz_yy_0,  \
                             ta2_yz_xzz_yy_1,  \
                             ta2_yz_xzz_yz_0,  \
                             ta2_yz_xzz_yz_1,  \
                             ta2_yz_xzz_z_0,   \
                             ta2_yz_xzz_z_1,   \
                             ta2_yz_xzz_zz_0,  \
                             ta2_yz_xzz_zz_1,  \
                             ta2_yz_zz_xz_0,   \
                             ta2_yz_zz_xz_1,   \
                             ta2_yz_zz_yy_0,   \
                             ta2_yz_zz_yy_1,   \
                             ta2_yz_zz_yz_0,   \
                             ta2_yz_zz_yz_1,   \
                             ta2_yz_zz_zz_0,   \
                             ta2_yz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxzz_xx_0[i] =
            ta2_yz_xx_xx_0[i] * fe_0 - ta2_yz_xx_xx_1[i] * fe_0 + ta1_y_xxz_xx_1[i] + ta2_yz_xxz_xx_0[i] * pa_z[i] - ta2_yz_xxz_xx_1[i] * pc_z[i];

        ta2_yz_xxzz_xy_0[i] =
            ta2_yz_xx_xy_0[i] * fe_0 - ta2_yz_xx_xy_1[i] * fe_0 + ta1_y_xxz_xy_1[i] + ta2_yz_xxz_xy_0[i] * pa_z[i] - ta2_yz_xxz_xy_1[i] * pc_z[i];

        ta2_yz_xxzz_xz_0[i] = ta2_yz_zz_xz_0[i] * fe_0 - ta2_yz_zz_xz_1[i] * fe_0 + ta2_yz_xzz_z_0[i] * fe_0 - ta2_yz_xzz_z_1[i] * fe_0 +
                              ta2_yz_xzz_xz_0[i] * pa_x[i] - ta2_yz_xzz_xz_1[i] * pc_x[i];

        ta2_yz_xxzz_yy_0[i] = ta2_yz_zz_yy_0[i] * fe_0 - ta2_yz_zz_yy_1[i] * fe_0 + ta2_yz_xzz_yy_0[i] * pa_x[i] - ta2_yz_xzz_yy_1[i] * pc_x[i];

        ta2_yz_xxzz_yz_0[i] = ta2_yz_zz_yz_0[i] * fe_0 - ta2_yz_zz_yz_1[i] * fe_0 + ta2_yz_xzz_yz_0[i] * pa_x[i] - ta2_yz_xzz_yz_1[i] * pc_x[i];

        ta2_yz_xxzz_zz_0[i] = ta2_yz_zz_zz_0[i] * fe_0 - ta2_yz_zz_zz_1[i] * fe_0 + ta2_yz_xzz_zz_0[i] * pa_x[i] - ta2_yz_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 396-402 components of targeted buffer : GD

    auto ta2_yz_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 396);

    auto ta2_yz_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 397);

    auto ta2_yz_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 398);

    auto ta2_yz_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 399);

    auto ta2_yz_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 400);

    auto ta2_yz_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 401);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yz_xyyy_xx_0, \
                             ta2_yz_xyyy_xy_0, \
                             ta2_yz_xyyy_xz_0, \
                             ta2_yz_xyyy_yy_0, \
                             ta2_yz_xyyy_yz_0, \
                             ta2_yz_xyyy_zz_0, \
                             ta2_yz_yyy_x_0,   \
                             ta2_yz_yyy_x_1,   \
                             ta2_yz_yyy_xx_0,  \
                             ta2_yz_yyy_xx_1,  \
                             ta2_yz_yyy_xy_0,  \
                             ta2_yz_yyy_xy_1,  \
                             ta2_yz_yyy_xz_0,  \
                             ta2_yz_yyy_xz_1,  \
                             ta2_yz_yyy_y_0,   \
                             ta2_yz_yyy_y_1,   \
                             ta2_yz_yyy_yy_0,  \
                             ta2_yz_yyy_yy_1,  \
                             ta2_yz_yyy_yz_0,  \
                             ta2_yz_yyy_yz_1,  \
                             ta2_yz_yyy_z_0,   \
                             ta2_yz_yyy_z_1,   \
                             ta2_yz_yyy_zz_0,  \
                             ta2_yz_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyy_xx_0[i] =
            2.0 * ta2_yz_yyy_x_0[i] * fe_0 - 2.0 * ta2_yz_yyy_x_1[i] * fe_0 + ta2_yz_yyy_xx_0[i] * pa_x[i] - ta2_yz_yyy_xx_1[i] * pc_x[i];

        ta2_yz_xyyy_xy_0[i] = ta2_yz_yyy_y_0[i] * fe_0 - ta2_yz_yyy_y_1[i] * fe_0 + ta2_yz_yyy_xy_0[i] * pa_x[i] - ta2_yz_yyy_xy_1[i] * pc_x[i];

        ta2_yz_xyyy_xz_0[i] = ta2_yz_yyy_z_0[i] * fe_0 - ta2_yz_yyy_z_1[i] * fe_0 + ta2_yz_yyy_xz_0[i] * pa_x[i] - ta2_yz_yyy_xz_1[i] * pc_x[i];

        ta2_yz_xyyy_yy_0[i] = ta2_yz_yyy_yy_0[i] * pa_x[i] - ta2_yz_yyy_yy_1[i] * pc_x[i];

        ta2_yz_xyyy_yz_0[i] = ta2_yz_yyy_yz_0[i] * pa_x[i] - ta2_yz_yyy_yz_1[i] * pc_x[i];

        ta2_yz_xyyy_zz_0[i] = ta2_yz_yyy_zz_0[i] * pa_x[i] - ta2_yz_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 402-408 components of targeted buffer : GD

    auto ta2_yz_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 402);

    auto ta2_yz_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 403);

    auto ta2_yz_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 404);

    auto ta2_yz_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 405);

    auto ta2_yz_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 406);

    auto ta2_yz_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 407);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xyy_xx_1,   \
                             ta1_y_xyy_xy_1,   \
                             ta2_yz_xyy_xx_0,  \
                             ta2_yz_xyy_xx_1,  \
                             ta2_yz_xyy_xy_0,  \
                             ta2_yz_xyy_xy_1,  \
                             ta2_yz_xyyz_xx_0, \
                             ta2_yz_xyyz_xy_0, \
                             ta2_yz_xyyz_xz_0, \
                             ta2_yz_xyyz_yy_0, \
                             ta2_yz_xyyz_yz_0, \
                             ta2_yz_xyyz_zz_0, \
                             ta2_yz_yyz_xz_0,  \
                             ta2_yz_yyz_xz_1,  \
                             ta2_yz_yyz_yy_0,  \
                             ta2_yz_yyz_yy_1,  \
                             ta2_yz_yyz_yz_0,  \
                             ta2_yz_yyz_yz_1,  \
                             ta2_yz_yyz_z_0,   \
                             ta2_yz_yyz_z_1,   \
                             ta2_yz_yyz_zz_0,  \
                             ta2_yz_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyz_xx_0[i] = ta1_y_xyy_xx_1[i] + ta2_yz_xyy_xx_0[i] * pa_z[i] - ta2_yz_xyy_xx_1[i] * pc_z[i];

        ta2_yz_xyyz_xy_0[i] = ta1_y_xyy_xy_1[i] + ta2_yz_xyy_xy_0[i] * pa_z[i] - ta2_yz_xyy_xy_1[i] * pc_z[i];

        ta2_yz_xyyz_xz_0[i] = ta2_yz_yyz_z_0[i] * fe_0 - ta2_yz_yyz_z_1[i] * fe_0 + ta2_yz_yyz_xz_0[i] * pa_x[i] - ta2_yz_yyz_xz_1[i] * pc_x[i];

        ta2_yz_xyyz_yy_0[i] = ta2_yz_yyz_yy_0[i] * pa_x[i] - ta2_yz_yyz_yy_1[i] * pc_x[i];

        ta2_yz_xyyz_yz_0[i] = ta2_yz_yyz_yz_0[i] * pa_x[i] - ta2_yz_yyz_yz_1[i] * pc_x[i];

        ta2_yz_xyyz_zz_0[i] = ta2_yz_yyz_zz_0[i] * pa_x[i] - ta2_yz_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 408-414 components of targeted buffer : GD

    auto ta2_yz_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 408);

    auto ta2_yz_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 409);

    auto ta2_yz_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 410);

    auto ta2_yz_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 411);

    auto ta2_yz_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 412);

    auto ta2_yz_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 413);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xzz_xx_1,   \
                             ta1_z_xzz_xz_1,   \
                             ta2_yz_xyzz_xx_0, \
                             ta2_yz_xyzz_xy_0, \
                             ta2_yz_xyzz_xz_0, \
                             ta2_yz_xyzz_yy_0, \
                             ta2_yz_xyzz_yz_0, \
                             ta2_yz_xyzz_zz_0, \
                             ta2_yz_xzz_xx_0,  \
                             ta2_yz_xzz_xx_1,  \
                             ta2_yz_xzz_xz_0,  \
                             ta2_yz_xzz_xz_1,  \
                             ta2_yz_yzz_xy_0,  \
                             ta2_yz_yzz_xy_1,  \
                             ta2_yz_yzz_y_0,   \
                             ta2_yz_yzz_y_1,   \
                             ta2_yz_yzz_yy_0,  \
                             ta2_yz_yzz_yy_1,  \
                             ta2_yz_yzz_yz_0,  \
                             ta2_yz_yzz_yz_1,  \
                             ta2_yz_yzz_zz_0,  \
                             ta2_yz_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyzz_xx_0[i] = ta1_z_xzz_xx_1[i] + ta2_yz_xzz_xx_0[i] * pa_y[i] - ta2_yz_xzz_xx_1[i] * pc_y[i];

        ta2_yz_xyzz_xy_0[i] = ta2_yz_yzz_y_0[i] * fe_0 - ta2_yz_yzz_y_1[i] * fe_0 + ta2_yz_yzz_xy_0[i] * pa_x[i] - ta2_yz_yzz_xy_1[i] * pc_x[i];

        ta2_yz_xyzz_xz_0[i] = ta1_z_xzz_xz_1[i] + ta2_yz_xzz_xz_0[i] * pa_y[i] - ta2_yz_xzz_xz_1[i] * pc_y[i];

        ta2_yz_xyzz_yy_0[i] = ta2_yz_yzz_yy_0[i] * pa_x[i] - ta2_yz_yzz_yy_1[i] * pc_x[i];

        ta2_yz_xyzz_yz_0[i] = ta2_yz_yzz_yz_0[i] * pa_x[i] - ta2_yz_yzz_yz_1[i] * pc_x[i];

        ta2_yz_xyzz_zz_0[i] = ta2_yz_yzz_zz_0[i] * pa_x[i] - ta2_yz_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 414-420 components of targeted buffer : GD

    auto ta2_yz_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 414);

    auto ta2_yz_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 415);

    auto ta2_yz_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 416);

    auto ta2_yz_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 417);

    auto ta2_yz_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 418);

    auto ta2_yz_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 419);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yz_xzzz_xx_0, \
                             ta2_yz_xzzz_xy_0, \
                             ta2_yz_xzzz_xz_0, \
                             ta2_yz_xzzz_yy_0, \
                             ta2_yz_xzzz_yz_0, \
                             ta2_yz_xzzz_zz_0, \
                             ta2_yz_zzz_x_0,   \
                             ta2_yz_zzz_x_1,   \
                             ta2_yz_zzz_xx_0,  \
                             ta2_yz_zzz_xx_1,  \
                             ta2_yz_zzz_xy_0,  \
                             ta2_yz_zzz_xy_1,  \
                             ta2_yz_zzz_xz_0,  \
                             ta2_yz_zzz_xz_1,  \
                             ta2_yz_zzz_y_0,   \
                             ta2_yz_zzz_y_1,   \
                             ta2_yz_zzz_yy_0,  \
                             ta2_yz_zzz_yy_1,  \
                             ta2_yz_zzz_yz_0,  \
                             ta2_yz_zzz_yz_1,  \
                             ta2_yz_zzz_z_0,   \
                             ta2_yz_zzz_z_1,   \
                             ta2_yz_zzz_zz_0,  \
                             ta2_yz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzzz_xx_0[i] =
            2.0 * ta2_yz_zzz_x_0[i] * fe_0 - 2.0 * ta2_yz_zzz_x_1[i] * fe_0 + ta2_yz_zzz_xx_0[i] * pa_x[i] - ta2_yz_zzz_xx_1[i] * pc_x[i];

        ta2_yz_xzzz_xy_0[i] = ta2_yz_zzz_y_0[i] * fe_0 - ta2_yz_zzz_y_1[i] * fe_0 + ta2_yz_zzz_xy_0[i] * pa_x[i] - ta2_yz_zzz_xy_1[i] * pc_x[i];

        ta2_yz_xzzz_xz_0[i] = ta2_yz_zzz_z_0[i] * fe_0 - ta2_yz_zzz_z_1[i] * fe_0 + ta2_yz_zzz_xz_0[i] * pa_x[i] - ta2_yz_zzz_xz_1[i] * pc_x[i];

        ta2_yz_xzzz_yy_0[i] = ta2_yz_zzz_yy_0[i] * pa_x[i] - ta2_yz_zzz_yy_1[i] * pc_x[i];

        ta2_yz_xzzz_yz_0[i] = ta2_yz_zzz_yz_0[i] * pa_x[i] - ta2_yz_zzz_yz_1[i] * pc_x[i];

        ta2_yz_xzzz_zz_0[i] = ta2_yz_zzz_zz_0[i] * pa_x[i] - ta2_yz_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 420-426 components of targeted buffer : GD

    auto ta2_yz_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 420);

    auto ta2_yz_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 421);

    auto ta2_yz_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 422);

    auto ta2_yz_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 423);

    auto ta2_yz_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 424);

    auto ta2_yz_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 425);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yyy_xx_1,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_xz_1,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_zz_1,   \
                             ta2_yz_yy_xx_0,   \
                             ta2_yz_yy_xx_1,   \
                             ta2_yz_yy_xy_0,   \
                             ta2_yz_yy_xy_1,   \
                             ta2_yz_yy_xz_0,   \
                             ta2_yz_yy_xz_1,   \
                             ta2_yz_yy_yy_0,   \
                             ta2_yz_yy_yy_1,   \
                             ta2_yz_yy_yz_0,   \
                             ta2_yz_yy_yz_1,   \
                             ta2_yz_yy_zz_0,   \
                             ta2_yz_yy_zz_1,   \
                             ta2_yz_yyy_x_0,   \
                             ta2_yz_yyy_x_1,   \
                             ta2_yz_yyy_xx_0,  \
                             ta2_yz_yyy_xx_1,  \
                             ta2_yz_yyy_xy_0,  \
                             ta2_yz_yyy_xy_1,  \
                             ta2_yz_yyy_xz_0,  \
                             ta2_yz_yyy_xz_1,  \
                             ta2_yz_yyy_y_0,   \
                             ta2_yz_yyy_y_1,   \
                             ta2_yz_yyy_yy_0,  \
                             ta2_yz_yyy_yy_1,  \
                             ta2_yz_yyy_yz_0,  \
                             ta2_yz_yyy_yz_1,  \
                             ta2_yz_yyy_z_0,   \
                             ta2_yz_yyy_z_1,   \
                             ta2_yz_yyy_zz_0,  \
                             ta2_yz_yyy_zz_1,  \
                             ta2_yz_yyyy_xx_0, \
                             ta2_yz_yyyy_xy_0, \
                             ta2_yz_yyyy_xz_0, \
                             ta2_yz_yyyy_yy_0, \
                             ta2_yz_yyyy_yz_0, \
                             ta2_yz_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyy_xx_0[i] = 3.0 * ta2_yz_yy_xx_0[i] * fe_0 - 3.0 * ta2_yz_yy_xx_1[i] * fe_0 + ta1_z_yyy_xx_1[i] + ta2_yz_yyy_xx_0[i] * pa_y[i] -
                              ta2_yz_yyy_xx_1[i] * pc_y[i];

        ta2_yz_yyyy_xy_0[i] = 3.0 * ta2_yz_yy_xy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xy_1[i] * fe_0 + ta2_yz_yyy_x_0[i] * fe_0 - ta2_yz_yyy_x_1[i] * fe_0 +
                              ta1_z_yyy_xy_1[i] + ta2_yz_yyy_xy_0[i] * pa_y[i] - ta2_yz_yyy_xy_1[i] * pc_y[i];

        ta2_yz_yyyy_xz_0[i] = 3.0 * ta2_yz_yy_xz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xz_1[i] * fe_0 + ta1_z_yyy_xz_1[i] + ta2_yz_yyy_xz_0[i] * pa_y[i] -
                              ta2_yz_yyy_xz_1[i] * pc_y[i];

        ta2_yz_yyyy_yy_0[i] = 3.0 * ta2_yz_yy_yy_0[i] * fe_0 - 3.0 * ta2_yz_yy_yy_1[i] * fe_0 + 2.0 * ta2_yz_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_yz_yyy_y_1[i] * fe_0 + ta1_z_yyy_yy_1[i] + ta2_yz_yyy_yy_0[i] * pa_y[i] - ta2_yz_yyy_yy_1[i] * pc_y[i];

        ta2_yz_yyyy_yz_0[i] = 3.0 * ta2_yz_yy_yz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yz_1[i] * fe_0 + ta2_yz_yyy_z_0[i] * fe_0 - ta2_yz_yyy_z_1[i] * fe_0 +
                              ta1_z_yyy_yz_1[i] + ta2_yz_yyy_yz_0[i] * pa_y[i] - ta2_yz_yyy_yz_1[i] * pc_y[i];

        ta2_yz_yyyy_zz_0[i] = 3.0 * ta2_yz_yy_zz_0[i] * fe_0 - 3.0 * ta2_yz_yy_zz_1[i] * fe_0 + ta1_z_yyy_zz_1[i] + ta2_yz_yyy_zz_0[i] * pa_y[i] -
                              ta2_yz_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 426-432 components of targeted buffer : GD

    auto ta2_yz_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 426);

    auto ta2_yz_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 427);

    auto ta2_yz_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 428);

    auto ta2_yz_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 429);

    auto ta2_yz_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 430);

    auto ta2_yz_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 431);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_z_yyz_xz_1,   \
                             ta1_z_yyz_zz_1,   \
                             ta2_yz_yyy_xx_0,  \
                             ta2_yz_yyy_xx_1,  \
                             ta2_yz_yyy_xy_0,  \
                             ta2_yz_yyy_xy_1,  \
                             ta2_yz_yyy_y_0,   \
                             ta2_yz_yyy_y_1,   \
                             ta2_yz_yyy_yy_0,  \
                             ta2_yz_yyy_yy_1,  \
                             ta2_yz_yyy_yz_0,  \
                             ta2_yz_yyy_yz_1,  \
                             ta2_yz_yyyz_xx_0, \
                             ta2_yz_yyyz_xy_0, \
                             ta2_yz_yyyz_xz_0, \
                             ta2_yz_yyyz_yy_0, \
                             ta2_yz_yyyz_yz_0, \
                             ta2_yz_yyyz_zz_0, \
                             ta2_yz_yyz_xz_0,  \
                             ta2_yz_yyz_xz_1,  \
                             ta2_yz_yyz_zz_0,  \
                             ta2_yz_yyz_zz_1,  \
                             ta2_yz_yz_xz_0,   \
                             ta2_yz_yz_xz_1,   \
                             ta2_yz_yz_zz_0,   \
                             ta2_yz_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyz_xx_0[i] = ta1_y_yyy_xx_1[i] + ta2_yz_yyy_xx_0[i] * pa_z[i] - ta2_yz_yyy_xx_1[i] * pc_z[i];

        ta2_yz_yyyz_xy_0[i] = ta1_y_yyy_xy_1[i] + ta2_yz_yyy_xy_0[i] * pa_z[i] - ta2_yz_yyy_xy_1[i] * pc_z[i];

        ta2_yz_yyyz_xz_0[i] = 2.0 * ta2_yz_yz_xz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xz_1[i] * fe_0 + ta1_z_yyz_xz_1[i] + ta2_yz_yyz_xz_0[i] * pa_y[i] -
                              ta2_yz_yyz_xz_1[i] * pc_y[i];

        ta2_yz_yyyz_yy_0[i] = ta1_y_yyy_yy_1[i] + ta2_yz_yyy_yy_0[i] * pa_z[i] - ta2_yz_yyy_yy_1[i] * pc_z[i];

        ta2_yz_yyyz_yz_0[i] =
            ta2_yz_yyy_y_0[i] * fe_0 - ta2_yz_yyy_y_1[i] * fe_0 + ta1_y_yyy_yz_1[i] + ta2_yz_yyy_yz_0[i] * pa_z[i] - ta2_yz_yyy_yz_1[i] * pc_z[i];

        ta2_yz_yyyz_zz_0[i] = 2.0 * ta2_yz_yz_zz_0[i] * fe_0 - 2.0 * ta2_yz_yz_zz_1[i] * fe_0 + ta1_z_yyz_zz_1[i] + ta2_yz_yyz_zz_0[i] * pa_y[i] -
                              ta2_yz_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 432-438 components of targeted buffer : GD

    auto ta2_yz_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 432);

    auto ta2_yz_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 433);

    auto ta2_yz_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 434);

    auto ta2_yz_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 435);

    auto ta2_yz_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 436);

    auto ta2_yz_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 437);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yyz_xy_1,   \
                             ta1_y_yyz_yy_1,   \
                             ta1_z_yzz_xx_1,   \
                             ta1_z_yzz_xz_1,   \
                             ta1_z_yzz_yz_1,   \
                             ta1_z_yzz_zz_1,   \
                             ta2_yz_yy_xy_0,   \
                             ta2_yz_yy_xy_1,   \
                             ta2_yz_yy_yy_0,   \
                             ta2_yz_yy_yy_1,   \
                             ta2_yz_yyz_xy_0,  \
                             ta2_yz_yyz_xy_1,  \
                             ta2_yz_yyz_yy_0,  \
                             ta2_yz_yyz_yy_1,  \
                             ta2_yz_yyzz_xx_0, \
                             ta2_yz_yyzz_xy_0, \
                             ta2_yz_yyzz_xz_0, \
                             ta2_yz_yyzz_yy_0, \
                             ta2_yz_yyzz_yz_0, \
                             ta2_yz_yyzz_zz_0, \
                             ta2_yz_yzz_xx_0,  \
                             ta2_yz_yzz_xx_1,  \
                             ta2_yz_yzz_xz_0,  \
                             ta2_yz_yzz_xz_1,  \
                             ta2_yz_yzz_yz_0,  \
                             ta2_yz_yzz_yz_1,  \
                             ta2_yz_yzz_z_0,   \
                             ta2_yz_yzz_z_1,   \
                             ta2_yz_yzz_zz_0,  \
                             ta2_yz_yzz_zz_1,  \
                             ta2_yz_zz_xx_0,   \
                             ta2_yz_zz_xx_1,   \
                             ta2_yz_zz_xz_0,   \
                             ta2_yz_zz_xz_1,   \
                             ta2_yz_zz_yz_0,   \
                             ta2_yz_zz_yz_1,   \
                             ta2_yz_zz_zz_0,   \
                             ta2_yz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyzz_xx_0[i] =
            ta2_yz_zz_xx_0[i] * fe_0 - ta2_yz_zz_xx_1[i] * fe_0 + ta1_z_yzz_xx_1[i] + ta2_yz_yzz_xx_0[i] * pa_y[i] - ta2_yz_yzz_xx_1[i] * pc_y[i];

        ta2_yz_yyzz_xy_0[i] =
            ta2_yz_yy_xy_0[i] * fe_0 - ta2_yz_yy_xy_1[i] * fe_0 + ta1_y_yyz_xy_1[i] + ta2_yz_yyz_xy_0[i] * pa_z[i] - ta2_yz_yyz_xy_1[i] * pc_z[i];

        ta2_yz_yyzz_xz_0[i] =
            ta2_yz_zz_xz_0[i] * fe_0 - ta2_yz_zz_xz_1[i] * fe_0 + ta1_z_yzz_xz_1[i] + ta2_yz_yzz_xz_0[i] * pa_y[i] - ta2_yz_yzz_xz_1[i] * pc_y[i];

        ta2_yz_yyzz_yy_0[i] =
            ta2_yz_yy_yy_0[i] * fe_0 - ta2_yz_yy_yy_1[i] * fe_0 + ta1_y_yyz_yy_1[i] + ta2_yz_yyz_yy_0[i] * pa_z[i] - ta2_yz_yyz_yy_1[i] * pc_z[i];

        ta2_yz_yyzz_yz_0[i] = ta2_yz_zz_yz_0[i] * fe_0 - ta2_yz_zz_yz_1[i] * fe_0 + ta2_yz_yzz_z_0[i] * fe_0 - ta2_yz_yzz_z_1[i] * fe_0 +
                              ta1_z_yzz_yz_1[i] + ta2_yz_yzz_yz_0[i] * pa_y[i] - ta2_yz_yzz_yz_1[i] * pc_y[i];

        ta2_yz_yyzz_zz_0[i] =
            ta2_yz_zz_zz_0[i] * fe_0 - ta2_yz_zz_zz_1[i] * fe_0 + ta1_z_yzz_zz_1[i] + ta2_yz_yzz_zz_0[i] * pa_y[i] - ta2_yz_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 438-444 components of targeted buffer : GD

    auto ta2_yz_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 438);

    auto ta2_yz_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 439);

    auto ta2_yz_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 440);

    auto ta2_yz_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 441);

    auto ta2_yz_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 442);

    auto ta2_yz_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 443);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_1,   \
                             ta2_yz_yzzz_xx_0, \
                             ta2_yz_yzzz_xy_0, \
                             ta2_yz_yzzz_xz_0, \
                             ta2_yz_yzzz_yy_0, \
                             ta2_yz_yzzz_yz_0, \
                             ta2_yz_yzzz_zz_0, \
                             ta2_yz_zzz_x_0,   \
                             ta2_yz_zzz_x_1,   \
                             ta2_yz_zzz_xx_0,  \
                             ta2_yz_zzz_xx_1,  \
                             ta2_yz_zzz_xy_0,  \
                             ta2_yz_zzz_xy_1,  \
                             ta2_yz_zzz_xz_0,  \
                             ta2_yz_zzz_xz_1,  \
                             ta2_yz_zzz_y_0,   \
                             ta2_yz_zzz_y_1,   \
                             ta2_yz_zzz_yy_0,  \
                             ta2_yz_zzz_yy_1,  \
                             ta2_yz_zzz_yz_0,  \
                             ta2_yz_zzz_yz_1,  \
                             ta2_yz_zzz_z_0,   \
                             ta2_yz_zzz_z_1,   \
                             ta2_yz_zzz_zz_0,  \
                             ta2_yz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzzz_xx_0[i] = ta1_z_zzz_xx_1[i] + ta2_yz_zzz_xx_0[i] * pa_y[i] - ta2_yz_zzz_xx_1[i] * pc_y[i];

        ta2_yz_yzzz_xy_0[i] =
            ta2_yz_zzz_x_0[i] * fe_0 - ta2_yz_zzz_x_1[i] * fe_0 + ta1_z_zzz_xy_1[i] + ta2_yz_zzz_xy_0[i] * pa_y[i] - ta2_yz_zzz_xy_1[i] * pc_y[i];

        ta2_yz_yzzz_xz_0[i] = ta1_z_zzz_xz_1[i] + ta2_yz_zzz_xz_0[i] * pa_y[i] - ta2_yz_zzz_xz_1[i] * pc_y[i];

        ta2_yz_yzzz_yy_0[i] = 2.0 * ta2_yz_zzz_y_0[i] * fe_0 - 2.0 * ta2_yz_zzz_y_1[i] * fe_0 + ta1_z_zzz_yy_1[i] + ta2_yz_zzz_yy_0[i] * pa_y[i] -
                              ta2_yz_zzz_yy_1[i] * pc_y[i];

        ta2_yz_yzzz_yz_0[i] =
            ta2_yz_zzz_z_0[i] * fe_0 - ta2_yz_zzz_z_1[i] * fe_0 + ta1_z_zzz_yz_1[i] + ta2_yz_zzz_yz_0[i] * pa_y[i] - ta2_yz_zzz_yz_1[i] * pc_y[i];

        ta2_yz_yzzz_zz_0[i] = ta1_z_zzz_zz_1[i] + ta2_yz_zzz_zz_0[i] * pa_y[i] - ta2_yz_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 444-450 components of targeted buffer : GD

    auto ta2_yz_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 444);

    auto ta2_yz_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 445);

    auto ta2_yz_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 446);

    auto ta2_yz_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 447);

    auto ta2_yz_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 448);

    auto ta2_yz_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 449);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_zzz_xx_1,   \
                             ta1_y_zzz_xy_1,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_zz_1,   \
                             ta2_yz_zz_xx_0,   \
                             ta2_yz_zz_xx_1,   \
                             ta2_yz_zz_xy_0,   \
                             ta2_yz_zz_xy_1,   \
                             ta2_yz_zz_xz_0,   \
                             ta2_yz_zz_xz_1,   \
                             ta2_yz_zz_yy_0,   \
                             ta2_yz_zz_yy_1,   \
                             ta2_yz_zz_yz_0,   \
                             ta2_yz_zz_yz_1,   \
                             ta2_yz_zz_zz_0,   \
                             ta2_yz_zz_zz_1,   \
                             ta2_yz_zzz_x_0,   \
                             ta2_yz_zzz_x_1,   \
                             ta2_yz_zzz_xx_0,  \
                             ta2_yz_zzz_xx_1,  \
                             ta2_yz_zzz_xy_0,  \
                             ta2_yz_zzz_xy_1,  \
                             ta2_yz_zzz_xz_0,  \
                             ta2_yz_zzz_xz_1,  \
                             ta2_yz_zzz_y_0,   \
                             ta2_yz_zzz_y_1,   \
                             ta2_yz_zzz_yy_0,  \
                             ta2_yz_zzz_yy_1,  \
                             ta2_yz_zzz_yz_0,  \
                             ta2_yz_zzz_yz_1,  \
                             ta2_yz_zzz_z_0,   \
                             ta2_yz_zzz_z_1,   \
                             ta2_yz_zzz_zz_0,  \
                             ta2_yz_zzz_zz_1,  \
                             ta2_yz_zzzz_xx_0, \
                             ta2_yz_zzzz_xy_0, \
                             ta2_yz_zzzz_xz_0, \
                             ta2_yz_zzzz_yy_0, \
                             ta2_yz_zzzz_yz_0, \
                             ta2_yz_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzzz_xx_0[i] = 3.0 * ta2_yz_zz_xx_0[i] * fe_0 - 3.0 * ta2_yz_zz_xx_1[i] * fe_0 + ta1_y_zzz_xx_1[i] + ta2_yz_zzz_xx_0[i] * pa_z[i] -
                              ta2_yz_zzz_xx_1[i] * pc_z[i];

        ta2_yz_zzzz_xy_0[i] = 3.0 * ta2_yz_zz_xy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xy_1[i] * fe_0 + ta1_y_zzz_xy_1[i] + ta2_yz_zzz_xy_0[i] * pa_z[i] -
                              ta2_yz_zzz_xy_1[i] * pc_z[i];

        ta2_yz_zzzz_xz_0[i] = 3.0 * ta2_yz_zz_xz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xz_1[i] * fe_0 + ta2_yz_zzz_x_0[i] * fe_0 - ta2_yz_zzz_x_1[i] * fe_0 +
                              ta1_y_zzz_xz_1[i] + ta2_yz_zzz_xz_0[i] * pa_z[i] - ta2_yz_zzz_xz_1[i] * pc_z[i];

        ta2_yz_zzzz_yy_0[i] = 3.0 * ta2_yz_zz_yy_0[i] * fe_0 - 3.0 * ta2_yz_zz_yy_1[i] * fe_0 + ta1_y_zzz_yy_1[i] + ta2_yz_zzz_yy_0[i] * pa_z[i] -
                              ta2_yz_zzz_yy_1[i] * pc_z[i];

        ta2_yz_zzzz_yz_0[i] = 3.0 * ta2_yz_zz_yz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yz_1[i] * fe_0 + ta2_yz_zzz_y_0[i] * fe_0 - ta2_yz_zzz_y_1[i] * fe_0 +
                              ta1_y_zzz_yz_1[i] + ta2_yz_zzz_yz_0[i] * pa_z[i] - ta2_yz_zzz_yz_1[i] * pc_z[i];

        ta2_yz_zzzz_zz_0[i] = 3.0 * ta2_yz_zz_zz_0[i] * fe_0 - 3.0 * ta2_yz_zz_zz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_yz_zzz_z_1[i] * fe_0 + ta1_y_zzz_zz_1[i] + ta2_yz_zzz_zz_0[i] * pa_z[i] - ta2_yz_zzz_zz_1[i] * pc_z[i];
    }

    // Set up 450-456 components of targeted buffer : GD

    auto ta2_zz_xxxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 450);

    auto ta2_zz_xxxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 451);

    auto ta2_zz_xxxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 452);

    auto ta2_zz_xxxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 453);

    auto ta2_zz_xxxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 454);

    auto ta2_zz_xxxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 455);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_zz_xx_xx_0,   \
                             ta2_zz_xx_xx_1,   \
                             ta2_zz_xx_xy_0,   \
                             ta2_zz_xx_xy_1,   \
                             ta2_zz_xx_xz_0,   \
                             ta2_zz_xx_xz_1,   \
                             ta2_zz_xx_yy_0,   \
                             ta2_zz_xx_yy_1,   \
                             ta2_zz_xx_yz_0,   \
                             ta2_zz_xx_yz_1,   \
                             ta2_zz_xx_zz_0,   \
                             ta2_zz_xx_zz_1,   \
                             ta2_zz_xxx_x_0,   \
                             ta2_zz_xxx_x_1,   \
                             ta2_zz_xxx_xx_0,  \
                             ta2_zz_xxx_xx_1,  \
                             ta2_zz_xxx_xy_0,  \
                             ta2_zz_xxx_xy_1,  \
                             ta2_zz_xxx_xz_0,  \
                             ta2_zz_xxx_xz_1,  \
                             ta2_zz_xxx_y_0,   \
                             ta2_zz_xxx_y_1,   \
                             ta2_zz_xxx_yy_0,  \
                             ta2_zz_xxx_yy_1,  \
                             ta2_zz_xxx_yz_0,  \
                             ta2_zz_xxx_yz_1,  \
                             ta2_zz_xxx_z_0,   \
                             ta2_zz_xxx_z_1,   \
                             ta2_zz_xxx_zz_0,  \
                             ta2_zz_xxx_zz_1,  \
                             ta2_zz_xxxx_xx_0, \
                             ta2_zz_xxxx_xy_0, \
                             ta2_zz_xxxx_xz_0, \
                             ta2_zz_xxxx_yy_0, \
                             ta2_zz_xxxx_yz_0, \
                             ta2_zz_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxx_xx_0[i] = 3.0 * ta2_zz_xx_xx_0[i] * fe_0 - 3.0 * ta2_zz_xx_xx_1[i] * fe_0 + 2.0 * ta2_zz_xxx_x_0[i] * fe_0 -
                              2.0 * ta2_zz_xxx_x_1[i] * fe_0 + ta2_zz_xxx_xx_0[i] * pa_x[i] - ta2_zz_xxx_xx_1[i] * pc_x[i];

        ta2_zz_xxxx_xy_0[i] = 3.0 * ta2_zz_xx_xy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xy_1[i] * fe_0 + ta2_zz_xxx_y_0[i] * fe_0 - ta2_zz_xxx_y_1[i] * fe_0 +
                              ta2_zz_xxx_xy_0[i] * pa_x[i] - ta2_zz_xxx_xy_1[i] * pc_x[i];

        ta2_zz_xxxx_xz_0[i] = 3.0 * ta2_zz_xx_xz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xz_1[i] * fe_0 + ta2_zz_xxx_z_0[i] * fe_0 - ta2_zz_xxx_z_1[i] * fe_0 +
                              ta2_zz_xxx_xz_0[i] * pa_x[i] - ta2_zz_xxx_xz_1[i] * pc_x[i];

        ta2_zz_xxxx_yy_0[i] =
            3.0 * ta2_zz_xx_yy_0[i] * fe_0 - 3.0 * ta2_zz_xx_yy_1[i] * fe_0 + ta2_zz_xxx_yy_0[i] * pa_x[i] - ta2_zz_xxx_yy_1[i] * pc_x[i];

        ta2_zz_xxxx_yz_0[i] =
            3.0 * ta2_zz_xx_yz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yz_1[i] * fe_0 + ta2_zz_xxx_yz_0[i] * pa_x[i] - ta2_zz_xxx_yz_1[i] * pc_x[i];

        ta2_zz_xxxx_zz_0[i] =
            3.0 * ta2_zz_xx_zz_0[i] * fe_0 - 3.0 * ta2_zz_xx_zz_1[i] * fe_0 + ta2_zz_xxx_zz_0[i] * pa_x[i] - ta2_zz_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 456-462 components of targeted buffer : GD

    auto ta2_zz_xxxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 456);

    auto ta2_zz_xxxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 457);

    auto ta2_zz_xxxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 458);

    auto ta2_zz_xxxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 459);

    auto ta2_zz_xxxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 460);

    auto ta2_zz_xxxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 461);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta2_zz_xxx_x_0,   \
                             ta2_zz_xxx_x_1,   \
                             ta2_zz_xxx_xx_0,  \
                             ta2_zz_xxx_xx_1,  \
                             ta2_zz_xxx_xy_0,  \
                             ta2_zz_xxx_xy_1,  \
                             ta2_zz_xxx_xz_0,  \
                             ta2_zz_xxx_xz_1,  \
                             ta2_zz_xxx_zz_0,  \
                             ta2_zz_xxx_zz_1,  \
                             ta2_zz_xxxy_xx_0, \
                             ta2_zz_xxxy_xy_0, \
                             ta2_zz_xxxy_xz_0, \
                             ta2_zz_xxxy_yy_0, \
                             ta2_zz_xxxy_yz_0, \
                             ta2_zz_xxxy_zz_0, \
                             ta2_zz_xxy_yy_0,  \
                             ta2_zz_xxy_yy_1,  \
                             ta2_zz_xxy_yz_0,  \
                             ta2_zz_xxy_yz_1,  \
                             ta2_zz_xy_yy_0,   \
                             ta2_zz_xy_yy_1,   \
                             ta2_zz_xy_yz_0,   \
                             ta2_zz_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxy_xx_0[i] = ta2_zz_xxx_xx_0[i] * pa_y[i] - ta2_zz_xxx_xx_1[i] * pc_y[i];

        ta2_zz_xxxy_xy_0[i] = ta2_zz_xxx_x_0[i] * fe_0 - ta2_zz_xxx_x_1[i] * fe_0 + ta2_zz_xxx_xy_0[i] * pa_y[i] - ta2_zz_xxx_xy_1[i] * pc_y[i];

        ta2_zz_xxxy_xz_0[i] = ta2_zz_xxx_xz_0[i] * pa_y[i] - ta2_zz_xxx_xz_1[i] * pc_y[i];

        ta2_zz_xxxy_yy_0[i] =
            2.0 * ta2_zz_xy_yy_0[i] * fe_0 - 2.0 * ta2_zz_xy_yy_1[i] * fe_0 + ta2_zz_xxy_yy_0[i] * pa_x[i] - ta2_zz_xxy_yy_1[i] * pc_x[i];

        ta2_zz_xxxy_yz_0[i] =
            2.0 * ta2_zz_xy_yz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yz_1[i] * fe_0 + ta2_zz_xxy_yz_0[i] * pa_x[i] - ta2_zz_xxy_yz_1[i] * pc_x[i];

        ta2_zz_xxxy_zz_0[i] = ta2_zz_xxx_zz_0[i] * pa_y[i] - ta2_zz_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 462-468 components of targeted buffer : GD

    auto ta2_zz_xxxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 462);

    auto ta2_zz_xxxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 463);

    auto ta2_zz_xxxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 464);

    auto ta2_zz_xxxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 465);

    auto ta2_zz_xxxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 466);

    auto ta2_zz_xxxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 467);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_yy_1,   \
                             ta2_zz_xxx_x_0,   \
                             ta2_zz_xxx_x_1,   \
                             ta2_zz_xxx_xx_0,  \
                             ta2_zz_xxx_xx_1,  \
                             ta2_zz_xxx_xy_0,  \
                             ta2_zz_xxx_xy_1,  \
                             ta2_zz_xxx_xz_0,  \
                             ta2_zz_xxx_xz_1,  \
                             ta2_zz_xxx_yy_0,  \
                             ta2_zz_xxx_yy_1,  \
                             ta2_zz_xxxz_xx_0, \
                             ta2_zz_xxxz_xy_0, \
                             ta2_zz_xxxz_xz_0, \
                             ta2_zz_xxxz_yy_0, \
                             ta2_zz_xxxz_yz_0, \
                             ta2_zz_xxxz_zz_0, \
                             ta2_zz_xxz_yz_0,  \
                             ta2_zz_xxz_yz_1,  \
                             ta2_zz_xxz_zz_0,  \
                             ta2_zz_xxz_zz_1,  \
                             ta2_zz_xz_yz_0,   \
                             ta2_zz_xz_yz_1,   \
                             ta2_zz_xz_zz_0,   \
                             ta2_zz_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxz_xx_0[i] = 2.0 * ta1_z_xxx_xx_1[i] + ta2_zz_xxx_xx_0[i] * pa_z[i] - ta2_zz_xxx_xx_1[i] * pc_z[i];

        ta2_zz_xxxz_xy_0[i] = 2.0 * ta1_z_xxx_xy_1[i] + ta2_zz_xxx_xy_0[i] * pa_z[i] - ta2_zz_xxx_xy_1[i] * pc_z[i];

        ta2_zz_xxxz_xz_0[i] = ta2_zz_xxx_x_0[i] * fe_0 - ta2_zz_xxx_x_1[i] * fe_0 + 2.0 * ta1_z_xxx_xz_1[i] + ta2_zz_xxx_xz_0[i] * pa_z[i] -
                              ta2_zz_xxx_xz_1[i] * pc_z[i];

        ta2_zz_xxxz_yy_0[i] = 2.0 * ta1_z_xxx_yy_1[i] + ta2_zz_xxx_yy_0[i] * pa_z[i] - ta2_zz_xxx_yy_1[i] * pc_z[i];

        ta2_zz_xxxz_yz_0[i] =
            2.0 * ta2_zz_xz_yz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yz_1[i] * fe_0 + ta2_zz_xxz_yz_0[i] * pa_x[i] - ta2_zz_xxz_yz_1[i] * pc_x[i];

        ta2_zz_xxxz_zz_0[i] =
            2.0 * ta2_zz_xz_zz_0[i] * fe_0 - 2.0 * ta2_zz_xz_zz_1[i] * fe_0 + ta2_zz_xxz_zz_0[i] * pa_x[i] - ta2_zz_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 468-474 components of targeted buffer : GD

    auto ta2_zz_xxyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 468);

    auto ta2_zz_xxyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 469);

    auto ta2_zz_xxyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 470);

    auto ta2_zz_xxyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 471);

    auto ta2_zz_xxyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 472);

    auto ta2_zz_xxyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 473);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta2_zz_xx_xx_0,   \
                             ta2_zz_xx_xx_1,   \
                             ta2_zz_xx_xz_0,   \
                             ta2_zz_xx_xz_1,   \
                             ta2_zz_xxy_xx_0,  \
                             ta2_zz_xxy_xx_1,  \
                             ta2_zz_xxy_xz_0,  \
                             ta2_zz_xxy_xz_1,  \
                             ta2_zz_xxyy_xx_0, \
                             ta2_zz_xxyy_xy_0, \
                             ta2_zz_xxyy_xz_0, \
                             ta2_zz_xxyy_yy_0, \
                             ta2_zz_xxyy_yz_0, \
                             ta2_zz_xxyy_zz_0, \
                             ta2_zz_xyy_xy_0,  \
                             ta2_zz_xyy_xy_1,  \
                             ta2_zz_xyy_y_0,   \
                             ta2_zz_xyy_y_1,   \
                             ta2_zz_xyy_yy_0,  \
                             ta2_zz_xyy_yy_1,  \
                             ta2_zz_xyy_yz_0,  \
                             ta2_zz_xyy_yz_1,  \
                             ta2_zz_xyy_zz_0,  \
                             ta2_zz_xyy_zz_1,  \
                             ta2_zz_yy_xy_0,   \
                             ta2_zz_yy_xy_1,   \
                             ta2_zz_yy_yy_0,   \
                             ta2_zz_yy_yy_1,   \
                             ta2_zz_yy_yz_0,   \
                             ta2_zz_yy_yz_1,   \
                             ta2_zz_yy_zz_0,   \
                             ta2_zz_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyy_xx_0[i] = ta2_zz_xx_xx_0[i] * fe_0 - ta2_zz_xx_xx_1[i] * fe_0 + ta2_zz_xxy_xx_0[i] * pa_y[i] - ta2_zz_xxy_xx_1[i] * pc_y[i];

        ta2_zz_xxyy_xy_0[i] = ta2_zz_yy_xy_0[i] * fe_0 - ta2_zz_yy_xy_1[i] * fe_0 + ta2_zz_xyy_y_0[i] * fe_0 - ta2_zz_xyy_y_1[i] * fe_0 +
                              ta2_zz_xyy_xy_0[i] * pa_x[i] - ta2_zz_xyy_xy_1[i] * pc_x[i];

        ta2_zz_xxyy_xz_0[i] = ta2_zz_xx_xz_0[i] * fe_0 - ta2_zz_xx_xz_1[i] * fe_0 + ta2_zz_xxy_xz_0[i] * pa_y[i] - ta2_zz_xxy_xz_1[i] * pc_y[i];

        ta2_zz_xxyy_yy_0[i] = ta2_zz_yy_yy_0[i] * fe_0 - ta2_zz_yy_yy_1[i] * fe_0 + ta2_zz_xyy_yy_0[i] * pa_x[i] - ta2_zz_xyy_yy_1[i] * pc_x[i];

        ta2_zz_xxyy_yz_0[i] = ta2_zz_yy_yz_0[i] * fe_0 - ta2_zz_yy_yz_1[i] * fe_0 + ta2_zz_xyy_yz_0[i] * pa_x[i] - ta2_zz_xyy_yz_1[i] * pc_x[i];

        ta2_zz_xxyy_zz_0[i] = ta2_zz_yy_zz_0[i] * fe_0 - ta2_zz_yy_zz_1[i] * fe_0 + ta2_zz_xyy_zz_0[i] * pa_x[i] - ta2_zz_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 474-480 components of targeted buffer : GD

    auto ta2_zz_xxyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 474);

    auto ta2_zz_xxyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 475);

    auto ta2_zz_xxyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 476);

    auto ta2_zz_xxyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 477);

    auto ta2_zz_xxyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 478);

    auto ta2_zz_xxyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 479);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_xxy_xy_1,   \
                             ta1_z_xxy_yy_1,   \
                             ta2_zz_xxy_xy_0,  \
                             ta2_zz_xxy_xy_1,  \
                             ta2_zz_xxy_yy_0,  \
                             ta2_zz_xxy_yy_1,  \
                             ta2_zz_xxyz_xx_0, \
                             ta2_zz_xxyz_xy_0, \
                             ta2_zz_xxyz_xz_0, \
                             ta2_zz_xxyz_yy_0, \
                             ta2_zz_xxyz_yz_0, \
                             ta2_zz_xxyz_zz_0, \
                             ta2_zz_xxz_xx_0,  \
                             ta2_zz_xxz_xx_1,  \
                             ta2_zz_xxz_xz_0,  \
                             ta2_zz_xxz_xz_1,  \
                             ta2_zz_xxz_zz_0,  \
                             ta2_zz_xxz_zz_1,  \
                             ta2_zz_xyz_yz_0,  \
                             ta2_zz_xyz_yz_1,  \
                             ta2_zz_yz_yz_0,   \
                             ta2_zz_yz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyz_xx_0[i] = ta2_zz_xxz_xx_0[i] * pa_y[i] - ta2_zz_xxz_xx_1[i] * pc_y[i];

        ta2_zz_xxyz_xy_0[i] = 2.0 * ta1_z_xxy_xy_1[i] + ta2_zz_xxy_xy_0[i] * pa_z[i] - ta2_zz_xxy_xy_1[i] * pc_z[i];

        ta2_zz_xxyz_xz_0[i] = ta2_zz_xxz_xz_0[i] * pa_y[i] - ta2_zz_xxz_xz_1[i] * pc_y[i];

        ta2_zz_xxyz_yy_0[i] = 2.0 * ta1_z_xxy_yy_1[i] + ta2_zz_xxy_yy_0[i] * pa_z[i] - ta2_zz_xxy_yy_1[i] * pc_z[i];

        ta2_zz_xxyz_yz_0[i] = ta2_zz_yz_yz_0[i] * fe_0 - ta2_zz_yz_yz_1[i] * fe_0 + ta2_zz_xyz_yz_0[i] * pa_x[i] - ta2_zz_xyz_yz_1[i] * pc_x[i];

        ta2_zz_xxyz_zz_0[i] = ta2_zz_xxz_zz_0[i] * pa_y[i] - ta2_zz_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 480-486 components of targeted buffer : GD

    auto ta2_zz_xxzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 480);

    auto ta2_zz_xxzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 481);

    auto ta2_zz_xxzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 482);

    auto ta2_zz_xxzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 483);

    auto ta2_zz_xxzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 484);

    auto ta2_zz_xxzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 485);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxz_xx_1,   \
                             ta1_z_xxz_xy_1,   \
                             ta2_zz_xx_xx_0,   \
                             ta2_zz_xx_xx_1,   \
                             ta2_zz_xx_xy_0,   \
                             ta2_zz_xx_xy_1,   \
                             ta2_zz_xxz_xx_0,  \
                             ta2_zz_xxz_xx_1,  \
                             ta2_zz_xxz_xy_0,  \
                             ta2_zz_xxz_xy_1,  \
                             ta2_zz_xxzz_xx_0, \
                             ta2_zz_xxzz_xy_0, \
                             ta2_zz_xxzz_xz_0, \
                             ta2_zz_xxzz_yy_0, \
                             ta2_zz_xxzz_yz_0, \
                             ta2_zz_xxzz_zz_0, \
                             ta2_zz_xzz_xz_0,  \
                             ta2_zz_xzz_xz_1,  \
                             ta2_zz_xzz_yy_0,  \
                             ta2_zz_xzz_yy_1,  \
                             ta2_zz_xzz_yz_0,  \
                             ta2_zz_xzz_yz_1,  \
                             ta2_zz_xzz_z_0,   \
                             ta2_zz_xzz_z_1,   \
                             ta2_zz_xzz_zz_0,  \
                             ta2_zz_xzz_zz_1,  \
                             ta2_zz_zz_xz_0,   \
                             ta2_zz_zz_xz_1,   \
                             ta2_zz_zz_yy_0,   \
                             ta2_zz_zz_yy_1,   \
                             ta2_zz_zz_yz_0,   \
                             ta2_zz_zz_yz_1,   \
                             ta2_zz_zz_zz_0,   \
                             ta2_zz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxzz_xx_0[i] = ta2_zz_xx_xx_0[i] * fe_0 - ta2_zz_xx_xx_1[i] * fe_0 + 2.0 * ta1_z_xxz_xx_1[i] + ta2_zz_xxz_xx_0[i] * pa_z[i] -
                              ta2_zz_xxz_xx_1[i] * pc_z[i];

        ta2_zz_xxzz_xy_0[i] = ta2_zz_xx_xy_0[i] * fe_0 - ta2_zz_xx_xy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xy_1[i] + ta2_zz_xxz_xy_0[i] * pa_z[i] -
                              ta2_zz_xxz_xy_1[i] * pc_z[i];

        ta2_zz_xxzz_xz_0[i] = ta2_zz_zz_xz_0[i] * fe_0 - ta2_zz_zz_xz_1[i] * fe_0 + ta2_zz_xzz_z_0[i] * fe_0 - ta2_zz_xzz_z_1[i] * fe_0 +
                              ta2_zz_xzz_xz_0[i] * pa_x[i] - ta2_zz_xzz_xz_1[i] * pc_x[i];

        ta2_zz_xxzz_yy_0[i] = ta2_zz_zz_yy_0[i] * fe_0 - ta2_zz_zz_yy_1[i] * fe_0 + ta2_zz_xzz_yy_0[i] * pa_x[i] - ta2_zz_xzz_yy_1[i] * pc_x[i];

        ta2_zz_xxzz_yz_0[i] = ta2_zz_zz_yz_0[i] * fe_0 - ta2_zz_zz_yz_1[i] * fe_0 + ta2_zz_xzz_yz_0[i] * pa_x[i] - ta2_zz_xzz_yz_1[i] * pc_x[i];

        ta2_zz_xxzz_zz_0[i] = ta2_zz_zz_zz_0[i] * fe_0 - ta2_zz_zz_zz_1[i] * fe_0 + ta2_zz_xzz_zz_0[i] * pa_x[i] - ta2_zz_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 486-492 components of targeted buffer : GD

    auto ta2_zz_xyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 486);

    auto ta2_zz_xyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 487);

    auto ta2_zz_xyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 488);

    auto ta2_zz_xyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 489);

    auto ta2_zz_xyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 490);

    auto ta2_zz_xyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 491);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_zz_xyyy_xx_0, \
                             ta2_zz_xyyy_xy_0, \
                             ta2_zz_xyyy_xz_0, \
                             ta2_zz_xyyy_yy_0, \
                             ta2_zz_xyyy_yz_0, \
                             ta2_zz_xyyy_zz_0, \
                             ta2_zz_yyy_x_0,   \
                             ta2_zz_yyy_x_1,   \
                             ta2_zz_yyy_xx_0,  \
                             ta2_zz_yyy_xx_1,  \
                             ta2_zz_yyy_xy_0,  \
                             ta2_zz_yyy_xy_1,  \
                             ta2_zz_yyy_xz_0,  \
                             ta2_zz_yyy_xz_1,  \
                             ta2_zz_yyy_y_0,   \
                             ta2_zz_yyy_y_1,   \
                             ta2_zz_yyy_yy_0,  \
                             ta2_zz_yyy_yy_1,  \
                             ta2_zz_yyy_yz_0,  \
                             ta2_zz_yyy_yz_1,  \
                             ta2_zz_yyy_z_0,   \
                             ta2_zz_yyy_z_1,   \
                             ta2_zz_yyy_zz_0,  \
                             ta2_zz_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyy_xx_0[i] =
            2.0 * ta2_zz_yyy_x_0[i] * fe_0 - 2.0 * ta2_zz_yyy_x_1[i] * fe_0 + ta2_zz_yyy_xx_0[i] * pa_x[i] - ta2_zz_yyy_xx_1[i] * pc_x[i];

        ta2_zz_xyyy_xy_0[i] = ta2_zz_yyy_y_0[i] * fe_0 - ta2_zz_yyy_y_1[i] * fe_0 + ta2_zz_yyy_xy_0[i] * pa_x[i] - ta2_zz_yyy_xy_1[i] * pc_x[i];

        ta2_zz_xyyy_xz_0[i] = ta2_zz_yyy_z_0[i] * fe_0 - ta2_zz_yyy_z_1[i] * fe_0 + ta2_zz_yyy_xz_0[i] * pa_x[i] - ta2_zz_yyy_xz_1[i] * pc_x[i];

        ta2_zz_xyyy_yy_0[i] = ta2_zz_yyy_yy_0[i] * pa_x[i] - ta2_zz_yyy_yy_1[i] * pc_x[i];

        ta2_zz_xyyy_yz_0[i] = ta2_zz_yyy_yz_0[i] * pa_x[i] - ta2_zz_yyy_yz_1[i] * pc_x[i];

        ta2_zz_xyyy_zz_0[i] = ta2_zz_yyy_zz_0[i] * pa_x[i] - ta2_zz_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 492-498 components of targeted buffer : GD

    auto ta2_zz_xyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 492);

    auto ta2_zz_xyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 493);

    auto ta2_zz_xyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 494);

    auto ta2_zz_xyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 495);

    auto ta2_zz_xyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 496);

    auto ta2_zz_xyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 497);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xyy_xx_1,   \
                             ta1_z_xyy_xy_1,   \
                             ta2_zz_xyy_xx_0,  \
                             ta2_zz_xyy_xx_1,  \
                             ta2_zz_xyy_xy_0,  \
                             ta2_zz_xyy_xy_1,  \
                             ta2_zz_xyyz_xx_0, \
                             ta2_zz_xyyz_xy_0, \
                             ta2_zz_xyyz_xz_0, \
                             ta2_zz_xyyz_yy_0, \
                             ta2_zz_xyyz_yz_0, \
                             ta2_zz_xyyz_zz_0, \
                             ta2_zz_yyz_xz_0,  \
                             ta2_zz_yyz_xz_1,  \
                             ta2_zz_yyz_yy_0,  \
                             ta2_zz_yyz_yy_1,  \
                             ta2_zz_yyz_yz_0,  \
                             ta2_zz_yyz_yz_1,  \
                             ta2_zz_yyz_z_0,   \
                             ta2_zz_yyz_z_1,   \
                             ta2_zz_yyz_zz_0,  \
                             ta2_zz_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyz_xx_0[i] = 2.0 * ta1_z_xyy_xx_1[i] + ta2_zz_xyy_xx_0[i] * pa_z[i] - ta2_zz_xyy_xx_1[i] * pc_z[i];

        ta2_zz_xyyz_xy_0[i] = 2.0 * ta1_z_xyy_xy_1[i] + ta2_zz_xyy_xy_0[i] * pa_z[i] - ta2_zz_xyy_xy_1[i] * pc_z[i];

        ta2_zz_xyyz_xz_0[i] = ta2_zz_yyz_z_0[i] * fe_0 - ta2_zz_yyz_z_1[i] * fe_0 + ta2_zz_yyz_xz_0[i] * pa_x[i] - ta2_zz_yyz_xz_1[i] * pc_x[i];

        ta2_zz_xyyz_yy_0[i] = ta2_zz_yyz_yy_0[i] * pa_x[i] - ta2_zz_yyz_yy_1[i] * pc_x[i];

        ta2_zz_xyyz_yz_0[i] = ta2_zz_yyz_yz_0[i] * pa_x[i] - ta2_zz_yyz_yz_1[i] * pc_x[i];

        ta2_zz_xyyz_zz_0[i] = ta2_zz_yyz_zz_0[i] * pa_x[i] - ta2_zz_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 498-504 components of targeted buffer : GD

    auto ta2_zz_xyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 498);

    auto ta2_zz_xyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 499);

    auto ta2_zz_xyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 500);

    auto ta2_zz_xyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 501);

    auto ta2_zz_xyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 502);

    auto ta2_zz_xyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 503);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta2_zz_xyzz_xx_0, \
                             ta2_zz_xyzz_xy_0, \
                             ta2_zz_xyzz_xz_0, \
                             ta2_zz_xyzz_yy_0, \
                             ta2_zz_xyzz_yz_0, \
                             ta2_zz_xyzz_zz_0, \
                             ta2_zz_xzz_xx_0,  \
                             ta2_zz_xzz_xx_1,  \
                             ta2_zz_xzz_xz_0,  \
                             ta2_zz_xzz_xz_1,  \
                             ta2_zz_yzz_xy_0,  \
                             ta2_zz_yzz_xy_1,  \
                             ta2_zz_yzz_y_0,   \
                             ta2_zz_yzz_y_1,   \
                             ta2_zz_yzz_yy_0,  \
                             ta2_zz_yzz_yy_1,  \
                             ta2_zz_yzz_yz_0,  \
                             ta2_zz_yzz_yz_1,  \
                             ta2_zz_yzz_zz_0,  \
                             ta2_zz_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyzz_xx_0[i] = ta2_zz_xzz_xx_0[i] * pa_y[i] - ta2_zz_xzz_xx_1[i] * pc_y[i];

        ta2_zz_xyzz_xy_0[i] = ta2_zz_yzz_y_0[i] * fe_0 - ta2_zz_yzz_y_1[i] * fe_0 + ta2_zz_yzz_xy_0[i] * pa_x[i] - ta2_zz_yzz_xy_1[i] * pc_x[i];

        ta2_zz_xyzz_xz_0[i] = ta2_zz_xzz_xz_0[i] * pa_y[i] - ta2_zz_xzz_xz_1[i] * pc_y[i];

        ta2_zz_xyzz_yy_0[i] = ta2_zz_yzz_yy_0[i] * pa_x[i] - ta2_zz_yzz_yy_1[i] * pc_x[i];

        ta2_zz_xyzz_yz_0[i] = ta2_zz_yzz_yz_0[i] * pa_x[i] - ta2_zz_yzz_yz_1[i] * pc_x[i];

        ta2_zz_xyzz_zz_0[i] = ta2_zz_yzz_zz_0[i] * pa_x[i] - ta2_zz_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 504-510 components of targeted buffer : GD

    auto ta2_zz_xzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 504);

    auto ta2_zz_xzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 505);

    auto ta2_zz_xzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 506);

    auto ta2_zz_xzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 507);

    auto ta2_zz_xzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 508);

    auto ta2_zz_xzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 509);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_zz_xzzz_xx_0, \
                             ta2_zz_xzzz_xy_0, \
                             ta2_zz_xzzz_xz_0, \
                             ta2_zz_xzzz_yy_0, \
                             ta2_zz_xzzz_yz_0, \
                             ta2_zz_xzzz_zz_0, \
                             ta2_zz_zzz_x_0,   \
                             ta2_zz_zzz_x_1,   \
                             ta2_zz_zzz_xx_0,  \
                             ta2_zz_zzz_xx_1,  \
                             ta2_zz_zzz_xy_0,  \
                             ta2_zz_zzz_xy_1,  \
                             ta2_zz_zzz_xz_0,  \
                             ta2_zz_zzz_xz_1,  \
                             ta2_zz_zzz_y_0,   \
                             ta2_zz_zzz_y_1,   \
                             ta2_zz_zzz_yy_0,  \
                             ta2_zz_zzz_yy_1,  \
                             ta2_zz_zzz_yz_0,  \
                             ta2_zz_zzz_yz_1,  \
                             ta2_zz_zzz_z_0,   \
                             ta2_zz_zzz_z_1,   \
                             ta2_zz_zzz_zz_0,  \
                             ta2_zz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzzz_xx_0[i] =
            2.0 * ta2_zz_zzz_x_0[i] * fe_0 - 2.0 * ta2_zz_zzz_x_1[i] * fe_0 + ta2_zz_zzz_xx_0[i] * pa_x[i] - ta2_zz_zzz_xx_1[i] * pc_x[i];

        ta2_zz_xzzz_xy_0[i] = ta2_zz_zzz_y_0[i] * fe_0 - ta2_zz_zzz_y_1[i] * fe_0 + ta2_zz_zzz_xy_0[i] * pa_x[i] - ta2_zz_zzz_xy_1[i] * pc_x[i];

        ta2_zz_xzzz_xz_0[i] = ta2_zz_zzz_z_0[i] * fe_0 - ta2_zz_zzz_z_1[i] * fe_0 + ta2_zz_zzz_xz_0[i] * pa_x[i] - ta2_zz_zzz_xz_1[i] * pc_x[i];

        ta2_zz_xzzz_yy_0[i] = ta2_zz_zzz_yy_0[i] * pa_x[i] - ta2_zz_zzz_yy_1[i] * pc_x[i];

        ta2_zz_xzzz_yz_0[i] = ta2_zz_zzz_yz_0[i] * pa_x[i] - ta2_zz_zzz_yz_1[i] * pc_x[i];

        ta2_zz_xzzz_zz_0[i] = ta2_zz_zzz_zz_0[i] * pa_x[i] - ta2_zz_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 510-516 components of targeted buffer : GD

    auto ta2_zz_yyyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 510);

    auto ta2_zz_yyyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 511);

    auto ta2_zz_yyyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 512);

    auto ta2_zz_yyyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 513);

    auto ta2_zz_yyyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 514);

    auto ta2_zz_yyyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 515);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_zz_yy_xx_0,   \
                             ta2_zz_yy_xx_1,   \
                             ta2_zz_yy_xy_0,   \
                             ta2_zz_yy_xy_1,   \
                             ta2_zz_yy_xz_0,   \
                             ta2_zz_yy_xz_1,   \
                             ta2_zz_yy_yy_0,   \
                             ta2_zz_yy_yy_1,   \
                             ta2_zz_yy_yz_0,   \
                             ta2_zz_yy_yz_1,   \
                             ta2_zz_yy_zz_0,   \
                             ta2_zz_yy_zz_1,   \
                             ta2_zz_yyy_x_0,   \
                             ta2_zz_yyy_x_1,   \
                             ta2_zz_yyy_xx_0,  \
                             ta2_zz_yyy_xx_1,  \
                             ta2_zz_yyy_xy_0,  \
                             ta2_zz_yyy_xy_1,  \
                             ta2_zz_yyy_xz_0,  \
                             ta2_zz_yyy_xz_1,  \
                             ta2_zz_yyy_y_0,   \
                             ta2_zz_yyy_y_1,   \
                             ta2_zz_yyy_yy_0,  \
                             ta2_zz_yyy_yy_1,  \
                             ta2_zz_yyy_yz_0,  \
                             ta2_zz_yyy_yz_1,  \
                             ta2_zz_yyy_z_0,   \
                             ta2_zz_yyy_z_1,   \
                             ta2_zz_yyy_zz_0,  \
                             ta2_zz_yyy_zz_1,  \
                             ta2_zz_yyyy_xx_0, \
                             ta2_zz_yyyy_xy_0, \
                             ta2_zz_yyyy_xz_0, \
                             ta2_zz_yyyy_yy_0, \
                             ta2_zz_yyyy_yz_0, \
                             ta2_zz_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyy_xx_0[i] =
            3.0 * ta2_zz_yy_xx_0[i] * fe_0 - 3.0 * ta2_zz_yy_xx_1[i] * fe_0 + ta2_zz_yyy_xx_0[i] * pa_y[i] - ta2_zz_yyy_xx_1[i] * pc_y[i];

        ta2_zz_yyyy_xy_0[i] = 3.0 * ta2_zz_yy_xy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xy_1[i] * fe_0 + ta2_zz_yyy_x_0[i] * fe_0 - ta2_zz_yyy_x_1[i] * fe_0 +
                              ta2_zz_yyy_xy_0[i] * pa_y[i] - ta2_zz_yyy_xy_1[i] * pc_y[i];

        ta2_zz_yyyy_xz_0[i] =
            3.0 * ta2_zz_yy_xz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xz_1[i] * fe_0 + ta2_zz_yyy_xz_0[i] * pa_y[i] - ta2_zz_yyy_xz_1[i] * pc_y[i];

        ta2_zz_yyyy_yy_0[i] = 3.0 * ta2_zz_yy_yy_0[i] * fe_0 - 3.0 * ta2_zz_yy_yy_1[i] * fe_0 + 2.0 * ta2_zz_yyy_y_0[i] * fe_0 -
                              2.0 * ta2_zz_yyy_y_1[i] * fe_0 + ta2_zz_yyy_yy_0[i] * pa_y[i] - ta2_zz_yyy_yy_1[i] * pc_y[i];

        ta2_zz_yyyy_yz_0[i] = 3.0 * ta2_zz_yy_yz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yz_1[i] * fe_0 + ta2_zz_yyy_z_0[i] * fe_0 - ta2_zz_yyy_z_1[i] * fe_0 +
                              ta2_zz_yyy_yz_0[i] * pa_y[i] - ta2_zz_yyy_yz_1[i] * pc_y[i];

        ta2_zz_yyyy_zz_0[i] =
            3.0 * ta2_zz_yy_zz_0[i] * fe_0 - 3.0 * ta2_zz_yy_zz_1[i] * fe_0 + ta2_zz_yyy_zz_0[i] * pa_y[i] - ta2_zz_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 516-522 components of targeted buffer : GD

    auto ta2_zz_yyyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 516);

    auto ta2_zz_yyyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 517);

    auto ta2_zz_yyyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 518);

    auto ta2_zz_yyyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 519);

    auto ta2_zz_yyyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 520);

    auto ta2_zz_yyyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 521);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyy_xx_1,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yz_1,   \
                             ta2_zz_yyy_xx_0,  \
                             ta2_zz_yyy_xx_1,  \
                             ta2_zz_yyy_xy_0,  \
                             ta2_zz_yyy_xy_1,  \
                             ta2_zz_yyy_y_0,   \
                             ta2_zz_yyy_y_1,   \
                             ta2_zz_yyy_yy_0,  \
                             ta2_zz_yyy_yy_1,  \
                             ta2_zz_yyy_yz_0,  \
                             ta2_zz_yyy_yz_1,  \
                             ta2_zz_yyyz_xx_0, \
                             ta2_zz_yyyz_xy_0, \
                             ta2_zz_yyyz_xz_0, \
                             ta2_zz_yyyz_yy_0, \
                             ta2_zz_yyyz_yz_0, \
                             ta2_zz_yyyz_zz_0, \
                             ta2_zz_yyz_xz_0,  \
                             ta2_zz_yyz_xz_1,  \
                             ta2_zz_yyz_zz_0,  \
                             ta2_zz_yyz_zz_1,  \
                             ta2_zz_yz_xz_0,   \
                             ta2_zz_yz_xz_1,   \
                             ta2_zz_yz_zz_0,   \
                             ta2_zz_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyz_xx_0[i] = 2.0 * ta1_z_yyy_xx_1[i] + ta2_zz_yyy_xx_0[i] * pa_z[i] - ta2_zz_yyy_xx_1[i] * pc_z[i];

        ta2_zz_yyyz_xy_0[i] = 2.0 * ta1_z_yyy_xy_1[i] + ta2_zz_yyy_xy_0[i] * pa_z[i] - ta2_zz_yyy_xy_1[i] * pc_z[i];

        ta2_zz_yyyz_xz_0[i] =
            2.0 * ta2_zz_yz_xz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xz_1[i] * fe_0 + ta2_zz_yyz_xz_0[i] * pa_y[i] - ta2_zz_yyz_xz_1[i] * pc_y[i];

        ta2_zz_yyyz_yy_0[i] = 2.0 * ta1_z_yyy_yy_1[i] + ta2_zz_yyy_yy_0[i] * pa_z[i] - ta2_zz_yyy_yy_1[i] * pc_z[i];

        ta2_zz_yyyz_yz_0[i] = ta2_zz_yyy_y_0[i] * fe_0 - ta2_zz_yyy_y_1[i] * fe_0 + 2.0 * ta1_z_yyy_yz_1[i] + ta2_zz_yyy_yz_0[i] * pa_z[i] -
                              ta2_zz_yyy_yz_1[i] * pc_z[i];

        ta2_zz_yyyz_zz_0[i] =
            2.0 * ta2_zz_yz_zz_0[i] * fe_0 - 2.0 * ta2_zz_yz_zz_1[i] * fe_0 + ta2_zz_yyz_zz_0[i] * pa_y[i] - ta2_zz_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 522-528 components of targeted buffer : GD

    auto ta2_zz_yyzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 522);

    auto ta2_zz_yyzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 523);

    auto ta2_zz_yyzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 524);

    auto ta2_zz_yyzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 525);

    auto ta2_zz_yyzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 526);

    auto ta2_zz_yyzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 527);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyz_xy_1,   \
                             ta1_z_yyz_yy_1,   \
                             ta2_zz_yy_xy_0,   \
                             ta2_zz_yy_xy_1,   \
                             ta2_zz_yy_yy_0,   \
                             ta2_zz_yy_yy_1,   \
                             ta2_zz_yyz_xy_0,  \
                             ta2_zz_yyz_xy_1,  \
                             ta2_zz_yyz_yy_0,  \
                             ta2_zz_yyz_yy_1,  \
                             ta2_zz_yyzz_xx_0, \
                             ta2_zz_yyzz_xy_0, \
                             ta2_zz_yyzz_xz_0, \
                             ta2_zz_yyzz_yy_0, \
                             ta2_zz_yyzz_yz_0, \
                             ta2_zz_yyzz_zz_0, \
                             ta2_zz_yzz_xx_0,  \
                             ta2_zz_yzz_xx_1,  \
                             ta2_zz_yzz_xz_0,  \
                             ta2_zz_yzz_xz_1,  \
                             ta2_zz_yzz_yz_0,  \
                             ta2_zz_yzz_yz_1,  \
                             ta2_zz_yzz_z_0,   \
                             ta2_zz_yzz_z_1,   \
                             ta2_zz_yzz_zz_0,  \
                             ta2_zz_yzz_zz_1,  \
                             ta2_zz_zz_xx_0,   \
                             ta2_zz_zz_xx_1,   \
                             ta2_zz_zz_xz_0,   \
                             ta2_zz_zz_xz_1,   \
                             ta2_zz_zz_yz_0,   \
                             ta2_zz_zz_yz_1,   \
                             ta2_zz_zz_zz_0,   \
                             ta2_zz_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyzz_xx_0[i] = ta2_zz_zz_xx_0[i] * fe_0 - ta2_zz_zz_xx_1[i] * fe_0 + ta2_zz_yzz_xx_0[i] * pa_y[i] - ta2_zz_yzz_xx_1[i] * pc_y[i];

        ta2_zz_yyzz_xy_0[i] = ta2_zz_yy_xy_0[i] * fe_0 - ta2_zz_yy_xy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xy_1[i] + ta2_zz_yyz_xy_0[i] * pa_z[i] -
                              ta2_zz_yyz_xy_1[i] * pc_z[i];

        ta2_zz_yyzz_xz_0[i] = ta2_zz_zz_xz_0[i] * fe_0 - ta2_zz_zz_xz_1[i] * fe_0 + ta2_zz_yzz_xz_0[i] * pa_y[i] - ta2_zz_yzz_xz_1[i] * pc_y[i];

        ta2_zz_yyzz_yy_0[i] = ta2_zz_yy_yy_0[i] * fe_0 - ta2_zz_yy_yy_1[i] * fe_0 + 2.0 * ta1_z_yyz_yy_1[i] + ta2_zz_yyz_yy_0[i] * pa_z[i] -
                              ta2_zz_yyz_yy_1[i] * pc_z[i];

        ta2_zz_yyzz_yz_0[i] = ta2_zz_zz_yz_0[i] * fe_0 - ta2_zz_zz_yz_1[i] * fe_0 + ta2_zz_yzz_z_0[i] * fe_0 - ta2_zz_yzz_z_1[i] * fe_0 +
                              ta2_zz_yzz_yz_0[i] * pa_y[i] - ta2_zz_yzz_yz_1[i] * pc_y[i];

        ta2_zz_yyzz_zz_0[i] = ta2_zz_zz_zz_0[i] * fe_0 - ta2_zz_zz_zz_1[i] * fe_0 + ta2_zz_yzz_zz_0[i] * pa_y[i] - ta2_zz_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 528-534 components of targeted buffer : GD

    auto ta2_zz_yzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 528);

    auto ta2_zz_yzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 529);

    auto ta2_zz_yzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 530);

    auto ta2_zz_yzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 531);

    auto ta2_zz_yzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 532);

    auto ta2_zz_yzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 533);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_zz_yzzz_xx_0, \
                             ta2_zz_yzzz_xy_0, \
                             ta2_zz_yzzz_xz_0, \
                             ta2_zz_yzzz_yy_0, \
                             ta2_zz_yzzz_yz_0, \
                             ta2_zz_yzzz_zz_0, \
                             ta2_zz_zzz_x_0,   \
                             ta2_zz_zzz_x_1,   \
                             ta2_zz_zzz_xx_0,  \
                             ta2_zz_zzz_xx_1,  \
                             ta2_zz_zzz_xy_0,  \
                             ta2_zz_zzz_xy_1,  \
                             ta2_zz_zzz_xz_0,  \
                             ta2_zz_zzz_xz_1,  \
                             ta2_zz_zzz_y_0,   \
                             ta2_zz_zzz_y_1,   \
                             ta2_zz_zzz_yy_0,  \
                             ta2_zz_zzz_yy_1,  \
                             ta2_zz_zzz_yz_0,  \
                             ta2_zz_zzz_yz_1,  \
                             ta2_zz_zzz_z_0,   \
                             ta2_zz_zzz_z_1,   \
                             ta2_zz_zzz_zz_0,  \
                             ta2_zz_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzzz_xx_0[i] = ta2_zz_zzz_xx_0[i] * pa_y[i] - ta2_zz_zzz_xx_1[i] * pc_y[i];

        ta2_zz_yzzz_xy_0[i] = ta2_zz_zzz_x_0[i] * fe_0 - ta2_zz_zzz_x_1[i] * fe_0 + ta2_zz_zzz_xy_0[i] * pa_y[i] - ta2_zz_zzz_xy_1[i] * pc_y[i];

        ta2_zz_yzzz_xz_0[i] = ta2_zz_zzz_xz_0[i] * pa_y[i] - ta2_zz_zzz_xz_1[i] * pc_y[i];

        ta2_zz_yzzz_yy_0[i] =
            2.0 * ta2_zz_zzz_y_0[i] * fe_0 - 2.0 * ta2_zz_zzz_y_1[i] * fe_0 + ta2_zz_zzz_yy_0[i] * pa_y[i] - ta2_zz_zzz_yy_1[i] * pc_y[i];

        ta2_zz_yzzz_yz_0[i] = ta2_zz_zzz_z_0[i] * fe_0 - ta2_zz_zzz_z_1[i] * fe_0 + ta2_zz_zzz_yz_0[i] * pa_y[i] - ta2_zz_zzz_yz_1[i] * pc_y[i];

        ta2_zz_yzzz_zz_0[i] = ta2_zz_zzz_zz_0[i] * pa_y[i] - ta2_zz_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 534-540 components of targeted buffer : GD

    auto ta2_zz_zzzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_gd + 534);

    auto ta2_zz_zzzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 535);

    auto ta2_zz_zzzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 536);

    auto ta2_zz_zzzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_gd + 537);

    auto ta2_zz_zzzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 538);

    auto ta2_zz_zzzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_gd + 539);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_1,   \
                             ta2_zz_zz_xx_0,   \
                             ta2_zz_zz_xx_1,   \
                             ta2_zz_zz_xy_0,   \
                             ta2_zz_zz_xy_1,   \
                             ta2_zz_zz_xz_0,   \
                             ta2_zz_zz_xz_1,   \
                             ta2_zz_zz_yy_0,   \
                             ta2_zz_zz_yy_1,   \
                             ta2_zz_zz_yz_0,   \
                             ta2_zz_zz_yz_1,   \
                             ta2_zz_zz_zz_0,   \
                             ta2_zz_zz_zz_1,   \
                             ta2_zz_zzz_x_0,   \
                             ta2_zz_zzz_x_1,   \
                             ta2_zz_zzz_xx_0,  \
                             ta2_zz_zzz_xx_1,  \
                             ta2_zz_zzz_xy_0,  \
                             ta2_zz_zzz_xy_1,  \
                             ta2_zz_zzz_xz_0,  \
                             ta2_zz_zzz_xz_1,  \
                             ta2_zz_zzz_y_0,   \
                             ta2_zz_zzz_y_1,   \
                             ta2_zz_zzz_yy_0,  \
                             ta2_zz_zzz_yy_1,  \
                             ta2_zz_zzz_yz_0,  \
                             ta2_zz_zzz_yz_1,  \
                             ta2_zz_zzz_z_0,   \
                             ta2_zz_zzz_z_1,   \
                             ta2_zz_zzz_zz_0,  \
                             ta2_zz_zzz_zz_1,  \
                             ta2_zz_zzzz_xx_0, \
                             ta2_zz_zzzz_xy_0, \
                             ta2_zz_zzzz_xz_0, \
                             ta2_zz_zzzz_yy_0, \
                             ta2_zz_zzzz_yz_0, \
                             ta2_zz_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzzz_xx_0[i] = 3.0 * ta2_zz_zz_xx_0[i] * fe_0 - 3.0 * ta2_zz_zz_xx_1[i] * fe_0 + 2.0 * ta1_z_zzz_xx_1[i] +
                              ta2_zz_zzz_xx_0[i] * pa_z[i] - ta2_zz_zzz_xx_1[i] * pc_z[i];

        ta2_zz_zzzz_xy_0[i] = 3.0 * ta2_zz_zz_xy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xy_1[i] +
                              ta2_zz_zzz_xy_0[i] * pa_z[i] - ta2_zz_zzz_xy_1[i] * pc_z[i];

        ta2_zz_zzzz_xz_0[i] = 3.0 * ta2_zz_zz_xz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xz_1[i] * fe_0 + ta2_zz_zzz_x_0[i] * fe_0 - ta2_zz_zzz_x_1[i] * fe_0 +
                              2.0 * ta1_z_zzz_xz_1[i] + ta2_zz_zzz_xz_0[i] * pa_z[i] - ta2_zz_zzz_xz_1[i] * pc_z[i];

        ta2_zz_zzzz_yy_0[i] = 3.0 * ta2_zz_zz_yy_0[i] * fe_0 - 3.0 * ta2_zz_zz_yy_1[i] * fe_0 + 2.0 * ta1_z_zzz_yy_1[i] +
                              ta2_zz_zzz_yy_0[i] * pa_z[i] - ta2_zz_zzz_yy_1[i] * pc_z[i];

        ta2_zz_zzzz_yz_0[i] = 3.0 * ta2_zz_zz_yz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yz_1[i] * fe_0 + ta2_zz_zzz_y_0[i] * fe_0 - ta2_zz_zzz_y_1[i] * fe_0 +
                              2.0 * ta1_z_zzz_yz_1[i] + ta2_zz_zzz_yz_0[i] * pa_z[i] - ta2_zz_zzz_yz_1[i] * pc_z[i];

        ta2_zz_zzzz_zz_0[i] = 3.0 * ta2_zz_zz_zz_0[i] * fe_0 - 3.0 * ta2_zz_zz_zz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_z_0[i] * fe_0 -
                              2.0 * ta2_zz_zzz_z_1[i] * fe_0 + 2.0 * ta1_z_zzz_zz_1[i] + ta2_zz_zzz_zz_0[i] * pa_z[i] - ta2_zz_zzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
