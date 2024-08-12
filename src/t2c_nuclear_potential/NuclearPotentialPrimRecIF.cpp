#include "NuclearPotentialPrimRecIF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_if(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_if,
                               const size_t              idx_npot_0_gf,
                               const size_t              idx_npot_1_gf,
                               const size_t              idx_npot_0_hd,
                               const size_t              idx_npot_1_hd,
                               const size_t              idx_npot_0_hf,
                               const size_t              idx_npot_1_hf,
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

    // Set up components of auxiliary buffer : GF

    auto ta_xxxx_xxx_0 = pbuffer.data(idx_npot_0_gf);

    auto ta_xxxx_xxy_0 = pbuffer.data(idx_npot_0_gf + 1);

    auto ta_xxxx_xxz_0 = pbuffer.data(idx_npot_0_gf + 2);

    auto ta_xxxx_xyy_0 = pbuffer.data(idx_npot_0_gf + 3);

    auto ta_xxxx_xyz_0 = pbuffer.data(idx_npot_0_gf + 4);

    auto ta_xxxx_xzz_0 = pbuffer.data(idx_npot_0_gf + 5);

    auto ta_xxxx_yyy_0 = pbuffer.data(idx_npot_0_gf + 6);

    auto ta_xxxx_yyz_0 = pbuffer.data(idx_npot_0_gf + 7);

    auto ta_xxxx_yzz_0 = pbuffer.data(idx_npot_0_gf + 8);

    auto ta_xxxx_zzz_0 = pbuffer.data(idx_npot_0_gf + 9);

    auto ta_xxxy_xxx_0 = pbuffer.data(idx_npot_0_gf + 10);

    auto ta_xxxy_xxz_0 = pbuffer.data(idx_npot_0_gf + 12);

    auto ta_xxxy_xzz_0 = pbuffer.data(idx_npot_0_gf + 15);

    auto ta_xxxy_yyy_0 = pbuffer.data(idx_npot_0_gf + 16);

    auto ta_xxxy_yyz_0 = pbuffer.data(idx_npot_0_gf + 17);

    auto ta_xxxy_yzz_0 = pbuffer.data(idx_npot_0_gf + 18);

    auto ta_xxxz_xxx_0 = pbuffer.data(idx_npot_0_gf + 20);

    auto ta_xxxz_xxy_0 = pbuffer.data(idx_npot_0_gf + 21);

    auto ta_xxxz_xxz_0 = pbuffer.data(idx_npot_0_gf + 22);

    auto ta_xxxz_xyy_0 = pbuffer.data(idx_npot_0_gf + 23);

    auto ta_xxxz_xzz_0 = pbuffer.data(idx_npot_0_gf + 25);

    auto ta_xxxz_yyz_0 = pbuffer.data(idx_npot_0_gf + 27);

    auto ta_xxxz_yzz_0 = pbuffer.data(idx_npot_0_gf + 28);

    auto ta_xxxz_zzz_0 = pbuffer.data(idx_npot_0_gf + 29);

    auto ta_xxyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 30);

    auto ta_xxyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 31);

    auto ta_xxyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 32);

    auto ta_xxyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 33);

    auto ta_xxyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 34);

    auto ta_xxyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 35);

    auto ta_xxyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 36);

    auto ta_xxyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 37);

    auto ta_xxyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 38);

    auto ta_xxyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 39);

    auto ta_xxyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 42);

    auto ta_xxyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 45);

    auto ta_xxyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 47);

    auto ta_xxyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 48);

    auto ta_xxzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 50);

    auto ta_xxzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 51);

    auto ta_xxzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 52);

    auto ta_xxzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 53);

    auto ta_xxzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 54);

    auto ta_xxzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 55);

    auto ta_xxzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 56);

    auto ta_xxzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 57);

    auto ta_xxzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 58);

    auto ta_xxzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 59);

    auto ta_xyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 61);

    auto ta_xyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 63);

    auto ta_xyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 64);

    auto ta_xyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 66);

    auto ta_xyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 67);

    auto ta_xyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 68);

    auto ta_xyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 69);

    auto ta_xyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 77);

    auto ta_xyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 78);

    auto ta_xyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 79);

    auto ta_xyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 86);

    auto ta_xyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 87);

    auto ta_xyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 88);

    auto ta_xzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 92);

    auto ta_xzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 94);

    auto ta_xzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 95);

    auto ta_xzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 96);

    auto ta_xzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 97);

    auto ta_xzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 98);

    auto ta_xzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 99);

    auto ta_yyyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 100);

    auto ta_yyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 101);

    auto ta_yyyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 102);

    auto ta_yyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 103);

    auto ta_yyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 104);

    auto ta_yyyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 105);

    auto ta_yyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 106);

    auto ta_yyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 107);

    auto ta_yyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 108);

    auto ta_yyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 109);

    auto ta_yyyz_xxy_0 = pbuffer.data(idx_npot_0_gf + 111);

    auto ta_yyyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 112);

    auto ta_yyyz_xyy_0 = pbuffer.data(idx_npot_0_gf + 113);

    auto ta_yyyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 115);

    auto ta_yyyz_yyy_0 = pbuffer.data(idx_npot_0_gf + 116);

    auto ta_yyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 117);

    auto ta_yyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 118);

    auto ta_yyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 119);

    auto ta_yyzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 120);

    auto ta_yyzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 121);

    auto ta_yyzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 122);

    auto ta_yyzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 123);

    auto ta_yyzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 124);

    auto ta_yyzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 125);

    auto ta_yyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 126);

    auto ta_yyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 127);

    auto ta_yyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 128);

    auto ta_yyzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 129);

    auto ta_yzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 130);

    auto ta_yzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 132);

    auto ta_yzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 134);

    auto ta_yzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 135);

    auto ta_yzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 136);

    auto ta_yzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 137);

    auto ta_yzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 138);

    auto ta_yzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 139);

    auto ta_zzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 140);

    auto ta_zzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 141);

    auto ta_zzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 142);

    auto ta_zzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 143);

    auto ta_zzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 144);

    auto ta_zzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 145);

    auto ta_zzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 146);

    auto ta_zzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 147);

    auto ta_zzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 148);

    auto ta_zzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 149);

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

    auto ta_xxxy_xxz_1 = pbuffer.data(idx_npot_1_gf + 12);

    auto ta_xxxy_xzz_1 = pbuffer.data(idx_npot_1_gf + 15);

    auto ta_xxxy_yyy_1 = pbuffer.data(idx_npot_1_gf + 16);

    auto ta_xxxy_yyz_1 = pbuffer.data(idx_npot_1_gf + 17);

    auto ta_xxxy_yzz_1 = pbuffer.data(idx_npot_1_gf + 18);

    auto ta_xxxz_xxx_1 = pbuffer.data(idx_npot_1_gf + 20);

    auto ta_xxxz_xxy_1 = pbuffer.data(idx_npot_1_gf + 21);

    auto ta_xxxz_xxz_1 = pbuffer.data(idx_npot_1_gf + 22);

    auto ta_xxxz_xyy_1 = pbuffer.data(idx_npot_1_gf + 23);

    auto ta_xxxz_xzz_1 = pbuffer.data(idx_npot_1_gf + 25);

    auto ta_xxxz_yyz_1 = pbuffer.data(idx_npot_1_gf + 27);

    auto ta_xxxz_yzz_1 = pbuffer.data(idx_npot_1_gf + 28);

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

    auto ta_xxyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 39);

    auto ta_xxyz_xxz_1 = pbuffer.data(idx_npot_1_gf + 42);

    auto ta_xxyz_xzz_1 = pbuffer.data(idx_npot_1_gf + 45);

    auto ta_xxyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 47);

    auto ta_xxyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 48);

    auto ta_xxzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 50);

    auto ta_xxzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 51);

    auto ta_xxzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 52);

    auto ta_xxzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 53);

    auto ta_xxzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 54);

    auto ta_xxzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 55);

    auto ta_xxzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 56);

    auto ta_xxzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 57);

    auto ta_xxzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 58);

    auto ta_xxzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 59);

    auto ta_xyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 61);

    auto ta_xyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 63);

    auto ta_xyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 64);

    auto ta_xyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 66);

    auto ta_xyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 67);

    auto ta_xyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 68);

    auto ta_xyyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 69);

    auto ta_xyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 77);

    auto ta_xyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 78);

    auto ta_xyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 79);

    auto ta_xyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 86);

    auto ta_xyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 87);

    auto ta_xyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 88);

    auto ta_xzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 92);

    auto ta_xzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 94);

    auto ta_xzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 95);

    auto ta_xzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 96);

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

    auto ta_yyyz_xxz_1 = pbuffer.data(idx_npot_1_gf + 112);

    auto ta_yyyz_xyy_1 = pbuffer.data(idx_npot_1_gf + 113);

    auto ta_yyyz_xzz_1 = pbuffer.data(idx_npot_1_gf + 115);

    auto ta_yyyz_yyy_1 = pbuffer.data(idx_npot_1_gf + 116);

    auto ta_yyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 117);

    auto ta_yyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 118);

    auto ta_yyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 119);

    auto ta_yyzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 120);

    auto ta_yyzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 121);

    auto ta_yyzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 122);

    auto ta_yyzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 123);

    auto ta_yyzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 124);

    auto ta_yyzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 125);

    auto ta_yyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 126);

    auto ta_yyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 127);

    auto ta_yyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 128);

    auto ta_yyzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 129);

    auto ta_yzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 130);

    auto ta_yzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 132);

    auto ta_yzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 134);

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

    // Set up components of auxiliary buffer : HD

    auto ta_xxxxx_xx_0 = pbuffer.data(idx_npot_0_hd);

    auto ta_xxxxx_xy_0 = pbuffer.data(idx_npot_0_hd + 1);

    auto ta_xxxxx_xz_0 = pbuffer.data(idx_npot_0_hd + 2);

    auto ta_xxxxx_yy_0 = pbuffer.data(idx_npot_0_hd + 3);

    auto ta_xxxxx_yz_0 = pbuffer.data(idx_npot_0_hd + 4);

    auto ta_xxxxx_zz_0 = pbuffer.data(idx_npot_0_hd + 5);

    auto ta_xxxxz_xz_0 = pbuffer.data(idx_npot_0_hd + 14);

    auto ta_xxxyy_xy_0 = pbuffer.data(idx_npot_0_hd + 19);

    auto ta_xxxyy_yy_0 = pbuffer.data(idx_npot_0_hd + 21);

    auto ta_xxxyy_yz_0 = pbuffer.data(idx_npot_0_hd + 22);

    auto ta_xxxzz_xx_0 = pbuffer.data(idx_npot_0_hd + 30);

    auto ta_xxxzz_xy_0 = pbuffer.data(idx_npot_0_hd + 31);

    auto ta_xxxzz_xz_0 = pbuffer.data(idx_npot_0_hd + 32);

    auto ta_xxxzz_yz_0 = pbuffer.data(idx_npot_0_hd + 34);

    auto ta_xxxzz_zz_0 = pbuffer.data(idx_npot_0_hd + 35);

    auto ta_xxyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 37);

    auto ta_xxyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 39);

    auto ta_xxyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 40);

    auto ta_xxzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 54);

    auto ta_xxzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 55);

    auto ta_xxzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 56);

    auto ta_xxzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 58);

    auto ta_xxzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 59);

    auto ta_xyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 61);

    auto ta_xyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 63);

    auto ta_xyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 64);

    auto ta_xyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 76);

    auto ta_xzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 86);

    auto ta_xzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 88);

    auto ta_xzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 89);

    auto ta_yyyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 90);

    auto ta_yyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 91);

    auto ta_yyyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 92);

    auto ta_yyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 93);

    auto ta_yyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 94);

    auto ta_yyyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 95);

    auto ta_yyyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 98);

    auto ta_yyyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 100);

    auto ta_yyyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 101);

    auto ta_yyyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 102);

    auto ta_yyyzz_xy_0 = pbuffer.data(idx_npot_0_hd + 103);

    auto ta_yyyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 104);

    auto ta_yyyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 105);

    auto ta_yyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 106);

    auto ta_yyyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 107);

    auto ta_yyzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 108);

    auto ta_yyzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 109);

    auto ta_yyzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 110);

    auto ta_yyzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 111);

    auto ta_yyzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 112);

    auto ta_yyzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 113);

    auto ta_yzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 115);

    auto ta_yzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 116);

    auto ta_yzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 117);

    auto ta_yzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 118);

    auto ta_yzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 119);

    auto ta_zzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 120);

    auto ta_zzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 121);

    auto ta_zzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 122);

    auto ta_zzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 123);

    auto ta_zzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 124);

    auto ta_zzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 125);

    // Set up components of auxiliary buffer : HD

    auto ta_xxxxx_xx_1 = pbuffer.data(idx_npot_1_hd);

    auto ta_xxxxx_xy_1 = pbuffer.data(idx_npot_1_hd + 1);

    auto ta_xxxxx_xz_1 = pbuffer.data(idx_npot_1_hd + 2);

    auto ta_xxxxx_yy_1 = pbuffer.data(idx_npot_1_hd + 3);

    auto ta_xxxxx_yz_1 = pbuffer.data(idx_npot_1_hd + 4);

    auto ta_xxxxx_zz_1 = pbuffer.data(idx_npot_1_hd + 5);

    auto ta_xxxxz_xz_1 = pbuffer.data(idx_npot_1_hd + 14);

    auto ta_xxxyy_xy_1 = pbuffer.data(idx_npot_1_hd + 19);

    auto ta_xxxyy_yy_1 = pbuffer.data(idx_npot_1_hd + 21);

    auto ta_xxxyy_yz_1 = pbuffer.data(idx_npot_1_hd + 22);

    auto ta_xxxzz_xx_1 = pbuffer.data(idx_npot_1_hd + 30);

    auto ta_xxxzz_xy_1 = pbuffer.data(idx_npot_1_hd + 31);

    auto ta_xxxzz_xz_1 = pbuffer.data(idx_npot_1_hd + 32);

    auto ta_xxxzz_yz_1 = pbuffer.data(idx_npot_1_hd + 34);

    auto ta_xxxzz_zz_1 = pbuffer.data(idx_npot_1_hd + 35);

    auto ta_xxyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 37);

    auto ta_xxyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 39);

    auto ta_xxyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 40);

    auto ta_xxzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 54);

    auto ta_xxzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 55);

    auto ta_xxzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 56);

    auto ta_xxzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 58);

    auto ta_xxzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 59);

    auto ta_xyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 61);

    auto ta_xyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 63);

    auto ta_xyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 64);

    auto ta_xyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 76);

    auto ta_xzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 86);

    auto ta_xzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 88);

    auto ta_xzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 89);

    auto ta_yyyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 90);

    auto ta_yyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 91);

    auto ta_yyyyy_xz_1 = pbuffer.data(idx_npot_1_hd + 92);

    auto ta_yyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 93);

    auto ta_yyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 94);

    auto ta_yyyyy_zz_1 = pbuffer.data(idx_npot_1_hd + 95);

    auto ta_yyyyz_xz_1 = pbuffer.data(idx_npot_1_hd + 98);

    auto ta_yyyyz_yz_1 = pbuffer.data(idx_npot_1_hd + 100);

    auto ta_yyyyz_zz_1 = pbuffer.data(idx_npot_1_hd + 101);

    auto ta_yyyzz_xx_1 = pbuffer.data(idx_npot_1_hd + 102);

    auto ta_yyyzz_xy_1 = pbuffer.data(idx_npot_1_hd + 103);

    auto ta_yyyzz_xz_1 = pbuffer.data(idx_npot_1_hd + 104);

    auto ta_yyyzz_yy_1 = pbuffer.data(idx_npot_1_hd + 105);

    auto ta_yyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 106);

    auto ta_yyyzz_zz_1 = pbuffer.data(idx_npot_1_hd + 107);

    auto ta_yyzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 108);

    auto ta_yyzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 109);

    auto ta_yyzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 110);

    auto ta_yyzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 111);

    auto ta_yyzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 112);

    auto ta_yyzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 113);

    auto ta_yzzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 115);

    auto ta_yzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 116);

    auto ta_yzzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 117);

    auto ta_yzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 118);

    auto ta_yzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 119);

    auto ta_zzzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 120);

    auto ta_zzzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 121);

    auto ta_zzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 122);

    auto ta_zzzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 123);

    auto ta_zzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 124);

    auto ta_zzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 125);

    // Set up components of auxiliary buffer : HF

    auto ta_xxxxx_xxx_0 = pbuffer.data(idx_npot_0_hf);

    auto ta_xxxxx_xxy_0 = pbuffer.data(idx_npot_0_hf + 1);

    auto ta_xxxxx_xxz_0 = pbuffer.data(idx_npot_0_hf + 2);

    auto ta_xxxxx_xyy_0 = pbuffer.data(idx_npot_0_hf + 3);

    auto ta_xxxxx_xyz_0 = pbuffer.data(idx_npot_0_hf + 4);

    auto ta_xxxxx_xzz_0 = pbuffer.data(idx_npot_0_hf + 5);

    auto ta_xxxxx_yyy_0 = pbuffer.data(idx_npot_0_hf + 6);

    auto ta_xxxxx_yyz_0 = pbuffer.data(idx_npot_0_hf + 7);

    auto ta_xxxxx_yzz_0 = pbuffer.data(idx_npot_0_hf + 8);

    auto ta_xxxxx_zzz_0 = pbuffer.data(idx_npot_0_hf + 9);

    auto ta_xxxxy_xxx_0 = pbuffer.data(idx_npot_0_hf + 10);

    auto ta_xxxxy_xxy_0 = pbuffer.data(idx_npot_0_hf + 11);

    auto ta_xxxxy_xxz_0 = pbuffer.data(idx_npot_0_hf + 12);

    auto ta_xxxxy_xyy_0 = pbuffer.data(idx_npot_0_hf + 13);

    auto ta_xxxxy_xzz_0 = pbuffer.data(idx_npot_0_hf + 15);

    auto ta_xxxxy_yyy_0 = pbuffer.data(idx_npot_0_hf + 16);

    auto ta_xxxxy_yyz_0 = pbuffer.data(idx_npot_0_hf + 17);

    auto ta_xxxxy_yzz_0 = pbuffer.data(idx_npot_0_hf + 18);

    auto ta_xxxxz_xxx_0 = pbuffer.data(idx_npot_0_hf + 20);

    auto ta_xxxxz_xxy_0 = pbuffer.data(idx_npot_0_hf + 21);

    auto ta_xxxxz_xxz_0 = pbuffer.data(idx_npot_0_hf + 22);

    auto ta_xxxxz_xyy_0 = pbuffer.data(idx_npot_0_hf + 23);

    auto ta_xxxxz_xyz_0 = pbuffer.data(idx_npot_0_hf + 24);

    auto ta_xxxxz_xzz_0 = pbuffer.data(idx_npot_0_hf + 25);

    auto ta_xxxxz_yyz_0 = pbuffer.data(idx_npot_0_hf + 27);

    auto ta_xxxxz_yzz_0 = pbuffer.data(idx_npot_0_hf + 28);

    auto ta_xxxxz_zzz_0 = pbuffer.data(idx_npot_0_hf + 29);

    auto ta_xxxyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 30);

    auto ta_xxxyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 31);

    auto ta_xxxyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 32);

    auto ta_xxxyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 33);

    auto ta_xxxyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 34);

    auto ta_xxxyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 35);

    auto ta_xxxyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 36);

    auto ta_xxxyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 37);

    auto ta_xxxyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 38);

    auto ta_xxxyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 39);

    auto ta_xxxyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 42);

    auto ta_xxxyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 45);

    auto ta_xxxyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 47);

    auto ta_xxxyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 48);

    auto ta_xxxzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 50);

    auto ta_xxxzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 51);

    auto ta_xxxzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 52);

    auto ta_xxxzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 53);

    auto ta_xxxzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 54);

    auto ta_xxxzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 55);

    auto ta_xxxzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 56);

    auto ta_xxxzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 57);

    auto ta_xxxzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 58);

    auto ta_xxxzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 59);

    auto ta_xxyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 60);

    auto ta_xxyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 61);

    auto ta_xxyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 62);

    auto ta_xxyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 63);

    auto ta_xxyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 64);

    auto ta_xxyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 65);

    auto ta_xxyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 66);

    auto ta_xxyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 67);

    auto ta_xxyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 68);

    auto ta_xxyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 69);

    auto ta_xxyyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 71);

    auto ta_xxyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 72);

    auto ta_xxyyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 73);

    auto ta_xxyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 75);

    auto ta_xxyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 77);

    auto ta_xxyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 78);

    auto ta_xxyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 79);

    auto ta_xxyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 80);

    auto ta_xxyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 82);

    auto ta_xxyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 85);

    auto ta_xxyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 86);

    auto ta_xxyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 87);

    auto ta_xxyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 88);

    auto ta_xxzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 90);

    auto ta_xxzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 91);

    auto ta_xxzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 92);

    auto ta_xxzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 93);

    auto ta_xxzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 94);

    auto ta_xxzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 95);

    auto ta_xxzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 96);

    auto ta_xxzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 97);

    auto ta_xxzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 98);

    auto ta_xxzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 99);

    auto ta_xyyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 100);

    auto ta_xyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 101);

    auto ta_xyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 103);

    auto ta_xyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 104);

    auto ta_xyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 106);

    auto ta_xyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 107);

    auto ta_xyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 108);

    auto ta_xyyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 109);

    auto ta_xyyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 117);

    auto ta_xyyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 118);

    auto ta_xyyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 119);

    auto ta_xyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 124);

    auto ta_xyyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 126);

    auto ta_xyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 127);

    auto ta_xyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 128);

    auto ta_xyyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 129);

    auto ta_xyzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 136);

    auto ta_xyzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 137);

    auto ta_xyzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 138);

    auto ta_xzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 140);

    auto ta_xzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 142);

    auto ta_xzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 144);

    auto ta_xzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 145);

    auto ta_xzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 146);

    auto ta_xzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 147);

    auto ta_xzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 148);

    auto ta_xzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 149);

    auto ta_yyyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 150);

    auto ta_yyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 151);

    auto ta_yyyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 152);

    auto ta_yyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 153);

    auto ta_yyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 154);

    auto ta_yyyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 155);

    auto ta_yyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 156);

    auto ta_yyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 157);

    auto ta_yyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 158);

    auto ta_yyyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 159);

    auto ta_yyyyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 161);

    auto ta_yyyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 162);

    auto ta_yyyyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 163);

    auto ta_yyyyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 164);

    auto ta_yyyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 165);

    auto ta_yyyyz_yyy_0 = pbuffer.data(idx_npot_0_hf + 166);

    auto ta_yyyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 167);

    auto ta_yyyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 168);

    auto ta_yyyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 169);

    auto ta_yyyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 170);

    auto ta_yyyzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 171);

    auto ta_yyyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 172);

    auto ta_yyyzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 173);

    auto ta_yyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 174);

    auto ta_yyyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 175);

    auto ta_yyyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 176);

    auto ta_yyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 177);

    auto ta_yyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 178);

    auto ta_yyyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 179);

    auto ta_yyzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 180);

    auto ta_yyzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 181);

    auto ta_yyzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 182);

    auto ta_yyzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 183);

    auto ta_yyzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 184);

    auto ta_yyzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 185);

    auto ta_yyzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 186);

    auto ta_yyzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 187);

    auto ta_yyzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 188);

    auto ta_yyzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 189);

    auto ta_yzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 190);

    auto ta_yzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 191);

    auto ta_yzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 192);

    auto ta_yzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 193);

    auto ta_yzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 194);

    auto ta_yzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 195);

    auto ta_yzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 196);

    auto ta_yzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 197);

    auto ta_yzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 198);

    auto ta_yzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 199);

    auto ta_zzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 200);

    auto ta_zzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 201);

    auto ta_zzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 202);

    auto ta_zzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 203);

    auto ta_zzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 204);

    auto ta_zzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 205);

    auto ta_zzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 206);

    auto ta_zzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 207);

    auto ta_zzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 208);

    auto ta_zzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 209);

    // Set up components of auxiliary buffer : HF

    auto ta_xxxxx_xxx_1 = pbuffer.data(idx_npot_1_hf);

    auto ta_xxxxx_xxy_1 = pbuffer.data(idx_npot_1_hf + 1);

    auto ta_xxxxx_xxz_1 = pbuffer.data(idx_npot_1_hf + 2);

    auto ta_xxxxx_xyy_1 = pbuffer.data(idx_npot_1_hf + 3);

    auto ta_xxxxx_xyz_1 = pbuffer.data(idx_npot_1_hf + 4);

    auto ta_xxxxx_xzz_1 = pbuffer.data(idx_npot_1_hf + 5);

    auto ta_xxxxx_yyy_1 = pbuffer.data(idx_npot_1_hf + 6);

    auto ta_xxxxx_yyz_1 = pbuffer.data(idx_npot_1_hf + 7);

    auto ta_xxxxx_yzz_1 = pbuffer.data(idx_npot_1_hf + 8);

    auto ta_xxxxx_zzz_1 = pbuffer.data(idx_npot_1_hf + 9);

    auto ta_xxxxy_xxx_1 = pbuffer.data(idx_npot_1_hf + 10);

    auto ta_xxxxy_xxy_1 = pbuffer.data(idx_npot_1_hf + 11);

    auto ta_xxxxy_xxz_1 = pbuffer.data(idx_npot_1_hf + 12);

    auto ta_xxxxy_xyy_1 = pbuffer.data(idx_npot_1_hf + 13);

    auto ta_xxxxy_xzz_1 = pbuffer.data(idx_npot_1_hf + 15);

    auto ta_xxxxy_yyy_1 = pbuffer.data(idx_npot_1_hf + 16);

    auto ta_xxxxy_yyz_1 = pbuffer.data(idx_npot_1_hf + 17);

    auto ta_xxxxy_yzz_1 = pbuffer.data(idx_npot_1_hf + 18);

    auto ta_xxxxz_xxx_1 = pbuffer.data(idx_npot_1_hf + 20);

    auto ta_xxxxz_xxy_1 = pbuffer.data(idx_npot_1_hf + 21);

    auto ta_xxxxz_xxz_1 = pbuffer.data(idx_npot_1_hf + 22);

    auto ta_xxxxz_xyy_1 = pbuffer.data(idx_npot_1_hf + 23);

    auto ta_xxxxz_xyz_1 = pbuffer.data(idx_npot_1_hf + 24);

    auto ta_xxxxz_xzz_1 = pbuffer.data(idx_npot_1_hf + 25);

    auto ta_xxxxz_yyz_1 = pbuffer.data(idx_npot_1_hf + 27);

    auto ta_xxxxz_yzz_1 = pbuffer.data(idx_npot_1_hf + 28);

    auto ta_xxxxz_zzz_1 = pbuffer.data(idx_npot_1_hf + 29);

    auto ta_xxxyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 30);

    auto ta_xxxyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 31);

    auto ta_xxxyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 32);

    auto ta_xxxyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 33);

    auto ta_xxxyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 34);

    auto ta_xxxyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 35);

    auto ta_xxxyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 36);

    auto ta_xxxyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 37);

    auto ta_xxxyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 38);

    auto ta_xxxyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 39);

    auto ta_xxxyz_xxz_1 = pbuffer.data(idx_npot_1_hf + 42);

    auto ta_xxxyz_xzz_1 = pbuffer.data(idx_npot_1_hf + 45);

    auto ta_xxxyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 47);

    auto ta_xxxyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 48);

    auto ta_xxxzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 50);

    auto ta_xxxzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 51);

    auto ta_xxxzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 52);

    auto ta_xxxzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 53);

    auto ta_xxxzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 54);

    auto ta_xxxzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 55);

    auto ta_xxxzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 56);

    auto ta_xxxzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 57);

    auto ta_xxxzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 58);

    auto ta_xxxzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 59);

    auto ta_xxyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 60);

    auto ta_xxyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 61);

    auto ta_xxyyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 62);

    auto ta_xxyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 63);

    auto ta_xxyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 64);

    auto ta_xxyyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 65);

    auto ta_xxyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 66);

    auto ta_xxyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 67);

    auto ta_xxyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 68);

    auto ta_xxyyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 69);

    auto ta_xxyyz_xxy_1 = pbuffer.data(idx_npot_1_hf + 71);

    auto ta_xxyyz_xxz_1 = pbuffer.data(idx_npot_1_hf + 72);

    auto ta_xxyyz_xyy_1 = pbuffer.data(idx_npot_1_hf + 73);

    auto ta_xxyyz_xzz_1 = pbuffer.data(idx_npot_1_hf + 75);

    auto ta_xxyyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 77);

    auto ta_xxyyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 78);

    auto ta_xxyyz_zzz_1 = pbuffer.data(idx_npot_1_hf + 79);

    auto ta_xxyzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 80);

    auto ta_xxyzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 82);

    auto ta_xxyzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 85);

    auto ta_xxyzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 86);

    auto ta_xxyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 87);

    auto ta_xxyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 88);

    auto ta_xxzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 90);

    auto ta_xxzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 91);

    auto ta_xxzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 92);

    auto ta_xxzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 93);

    auto ta_xxzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 94);

    auto ta_xxzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 95);

    auto ta_xxzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 96);

    auto ta_xxzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 97);

    auto ta_xxzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 98);

    auto ta_xxzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 99);

    auto ta_xyyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 100);

    auto ta_xyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 101);

    auto ta_xyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 103);

    auto ta_xyyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 104);

    auto ta_xyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 106);

    auto ta_xyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 107);

    auto ta_xyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 108);

    auto ta_xyyyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 109);

    auto ta_xyyyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 117);

    auto ta_xyyyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 118);

    auto ta_xyyyz_zzz_1 = pbuffer.data(idx_npot_1_hf + 119);

    auto ta_xyyzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 124);

    auto ta_xyyzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 126);

    auto ta_xyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 127);

    auto ta_xyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 128);

    auto ta_xyyzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 129);

    auto ta_xyzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 136);

    auto ta_xyzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 137);

    auto ta_xyzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 138);

    auto ta_xzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 140);

    auto ta_xzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 142);

    auto ta_xzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 144);

    auto ta_xzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 145);

    auto ta_xzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 146);

    auto ta_xzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 147);

    auto ta_xzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 148);

    auto ta_xzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 149);

    auto ta_yyyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 150);

    auto ta_yyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 151);

    auto ta_yyyyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 152);

    auto ta_yyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 153);

    auto ta_yyyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 154);

    auto ta_yyyyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 155);

    auto ta_yyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 156);

    auto ta_yyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 157);

    auto ta_yyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 158);

    auto ta_yyyyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 159);

    auto ta_yyyyz_xxy_1 = pbuffer.data(idx_npot_1_hf + 161);

    auto ta_yyyyz_xxz_1 = pbuffer.data(idx_npot_1_hf + 162);

    auto ta_yyyyz_xyy_1 = pbuffer.data(idx_npot_1_hf + 163);

    auto ta_yyyyz_xyz_1 = pbuffer.data(idx_npot_1_hf + 164);

    auto ta_yyyyz_xzz_1 = pbuffer.data(idx_npot_1_hf + 165);

    auto ta_yyyyz_yyy_1 = pbuffer.data(idx_npot_1_hf + 166);

    auto ta_yyyyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 167);

    auto ta_yyyyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 168);

    auto ta_yyyyz_zzz_1 = pbuffer.data(idx_npot_1_hf + 169);

    auto ta_yyyzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 170);

    auto ta_yyyzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 171);

    auto ta_yyyzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 172);

    auto ta_yyyzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 173);

    auto ta_yyyzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 174);

    auto ta_yyyzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 175);

    auto ta_yyyzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 176);

    auto ta_yyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 177);

    auto ta_yyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 178);

    auto ta_yyyzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 179);

    auto ta_yyzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 180);

    auto ta_yyzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 181);

    auto ta_yyzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 182);

    auto ta_yyzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 183);

    auto ta_yyzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 184);

    auto ta_yyzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 185);

    auto ta_yyzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 186);

    auto ta_yyzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 187);

    auto ta_yyzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 188);

    auto ta_yyzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 189);

    auto ta_yzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 190);

    auto ta_yzzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 191);

    auto ta_yzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 192);

    auto ta_yzzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 193);

    auto ta_yzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 194);

    auto ta_yzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 195);

    auto ta_yzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 196);

    auto ta_yzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 197);

    auto ta_yzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 198);

    auto ta_yzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 199);

    auto ta_zzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 200);

    auto ta_zzzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 201);

    auto ta_zzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 202);

    auto ta_zzzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 203);

    auto ta_zzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 204);

    auto ta_zzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 205);

    auto ta_zzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 206);

    auto ta_zzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 207);

    auto ta_zzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 208);

    auto ta_zzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 209);

    // Set up 0-10 components of targeted buffer : IF

    auto ta_xxxxxx_xxx_0 = pbuffer.data(idx_npot_0_if);

    auto ta_xxxxxx_xxy_0 = pbuffer.data(idx_npot_0_if + 1);

    auto ta_xxxxxx_xxz_0 = pbuffer.data(idx_npot_0_if + 2);

    auto ta_xxxxxx_xyy_0 = pbuffer.data(idx_npot_0_if + 3);

    auto ta_xxxxxx_xyz_0 = pbuffer.data(idx_npot_0_if + 4);

    auto ta_xxxxxx_xzz_0 = pbuffer.data(idx_npot_0_if + 5);

    auto ta_xxxxxx_yyy_0 = pbuffer.data(idx_npot_0_if + 6);

    auto ta_xxxxxx_yyz_0 = pbuffer.data(idx_npot_0_if + 7);

    auto ta_xxxxxx_yzz_0 = pbuffer.data(idx_npot_0_if + 8);

    auto ta_xxxxxx_zzz_0 = pbuffer.data(idx_npot_0_if + 9);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxy_0,   \
                             ta_xxxx_xxy_1,   \
                             ta_xxxx_xxz_0,   \
                             ta_xxxx_xxz_1,   \
                             ta_xxxx_xyy_0,   \
                             ta_xxxx_xyy_1,   \
                             ta_xxxx_xyz_0,   \
                             ta_xxxx_xyz_1,   \
                             ta_xxxx_xzz_0,   \
                             ta_xxxx_xzz_1,   \
                             ta_xxxx_yyy_0,   \
                             ta_xxxx_yyy_1,   \
                             ta_xxxx_yyz_0,   \
                             ta_xxxx_yyz_1,   \
                             ta_xxxx_yzz_0,   \
                             ta_xxxx_yzz_1,   \
                             ta_xxxx_zzz_0,   \
                             ta_xxxx_zzz_1,   \
                             ta_xxxxx_xx_0,   \
                             ta_xxxxx_xx_1,   \
                             ta_xxxxx_xxx_0,  \
                             ta_xxxxx_xxx_1,  \
                             ta_xxxxx_xxy_0,  \
                             ta_xxxxx_xxy_1,  \
                             ta_xxxxx_xxz_0,  \
                             ta_xxxxx_xxz_1,  \
                             ta_xxxxx_xy_0,   \
                             ta_xxxxx_xy_1,   \
                             ta_xxxxx_xyy_0,  \
                             ta_xxxxx_xyy_1,  \
                             ta_xxxxx_xyz_0,  \
                             ta_xxxxx_xyz_1,  \
                             ta_xxxxx_xz_0,   \
                             ta_xxxxx_xz_1,   \
                             ta_xxxxx_xzz_0,  \
                             ta_xxxxx_xzz_1,  \
                             ta_xxxxx_yy_0,   \
                             ta_xxxxx_yy_1,   \
                             ta_xxxxx_yyy_0,  \
                             ta_xxxxx_yyy_1,  \
                             ta_xxxxx_yyz_0,  \
                             ta_xxxxx_yyz_1,  \
                             ta_xxxxx_yz_0,   \
                             ta_xxxxx_yz_1,   \
                             ta_xxxxx_yzz_0,  \
                             ta_xxxxx_yzz_1,  \
                             ta_xxxxx_zz_0,   \
                             ta_xxxxx_zz_1,   \
                             ta_xxxxx_zzz_0,  \
                             ta_xxxxx_zzz_1,  \
                             ta_xxxxxx_xxx_0, \
                             ta_xxxxxx_xxy_0, \
                             ta_xxxxxx_xxz_0, \
                             ta_xxxxxx_xyy_0, \
                             ta_xxxxxx_xyz_0, \
                             ta_xxxxxx_xzz_0, \
                             ta_xxxxxx_yyy_0, \
                             ta_xxxxxx_yyz_0, \
                             ta_xxxxxx_yzz_0, \
                             ta_xxxxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_xxx_0[i] = 5.0 * ta_xxxx_xxx_0[i] * fe_0 - 5.0 * ta_xxxx_xxx_1[i] * fe_0 + 3.0 * ta_xxxxx_xx_0[i] * fe_0 -
                             3.0 * ta_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxx_0[i] * pa_x[i] - ta_xxxxx_xxx_1[i] * pc_x[i];

        ta_xxxxxx_xxy_0[i] = 5.0 * ta_xxxx_xxy_0[i] * fe_0 - 5.0 * ta_xxxx_xxy_1[i] * fe_0 + 2.0 * ta_xxxxx_xy_0[i] * fe_0 -
                             2.0 * ta_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xxy_0[i] * pa_x[i] - ta_xxxxx_xxy_1[i] * pc_x[i];

        ta_xxxxxx_xxz_0[i] = 5.0 * ta_xxxx_xxz_0[i] * fe_0 - 5.0 * ta_xxxx_xxz_1[i] * fe_0 + 2.0 * ta_xxxxx_xz_0[i] * fe_0 -
                             2.0 * ta_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xxz_0[i] * pa_x[i] - ta_xxxxx_xxz_1[i] * pc_x[i];

        ta_xxxxxx_xyy_0[i] = 5.0 * ta_xxxx_xyy_0[i] * fe_0 - 5.0 * ta_xxxx_xyy_1[i] * fe_0 + ta_xxxxx_yy_0[i] * fe_0 - ta_xxxxx_yy_1[i] * fe_0 +
                             ta_xxxxx_xyy_0[i] * pa_x[i] - ta_xxxxx_xyy_1[i] * pc_x[i];

        ta_xxxxxx_xyz_0[i] = 5.0 * ta_xxxx_xyz_0[i] * fe_0 - 5.0 * ta_xxxx_xyz_1[i] * fe_0 + ta_xxxxx_yz_0[i] * fe_0 - ta_xxxxx_yz_1[i] * fe_0 +
                             ta_xxxxx_xyz_0[i] * pa_x[i] - ta_xxxxx_xyz_1[i] * pc_x[i];

        ta_xxxxxx_xzz_0[i] = 5.0 * ta_xxxx_xzz_0[i] * fe_0 - 5.0 * ta_xxxx_xzz_1[i] * fe_0 + ta_xxxxx_zz_0[i] * fe_0 - ta_xxxxx_zz_1[i] * fe_0 +
                             ta_xxxxx_xzz_0[i] * pa_x[i] - ta_xxxxx_xzz_1[i] * pc_x[i];

        ta_xxxxxx_yyy_0[i] =
            5.0 * ta_xxxx_yyy_0[i] * fe_0 - 5.0 * ta_xxxx_yyy_1[i] * fe_0 + ta_xxxxx_yyy_0[i] * pa_x[i] - ta_xxxxx_yyy_1[i] * pc_x[i];

        ta_xxxxxx_yyz_0[i] =
            5.0 * ta_xxxx_yyz_0[i] * fe_0 - 5.0 * ta_xxxx_yyz_1[i] * fe_0 + ta_xxxxx_yyz_0[i] * pa_x[i] - ta_xxxxx_yyz_1[i] * pc_x[i];

        ta_xxxxxx_yzz_0[i] =
            5.0 * ta_xxxx_yzz_0[i] * fe_0 - 5.0 * ta_xxxx_yzz_1[i] * fe_0 + ta_xxxxx_yzz_0[i] * pa_x[i] - ta_xxxxx_yzz_1[i] * pc_x[i];

        ta_xxxxxx_zzz_0[i] =
            5.0 * ta_xxxx_zzz_0[i] * fe_0 - 5.0 * ta_xxxx_zzz_1[i] * fe_0 + ta_xxxxx_zzz_0[i] * pa_x[i] - ta_xxxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : IF

    auto ta_xxxxxy_xxx_0 = pbuffer.data(idx_npot_0_if + 10);

    auto ta_xxxxxy_xxy_0 = pbuffer.data(idx_npot_0_if + 11);

    auto ta_xxxxxy_xxz_0 = pbuffer.data(idx_npot_0_if + 12);

    auto ta_xxxxxy_xyy_0 = pbuffer.data(idx_npot_0_if + 13);

    auto ta_xxxxxy_xyz_0 = pbuffer.data(idx_npot_0_if + 14);

    auto ta_xxxxxy_xzz_0 = pbuffer.data(idx_npot_0_if + 15);

    auto ta_xxxxxy_yyy_0 = pbuffer.data(idx_npot_0_if + 16);

    auto ta_xxxxxy_yyz_0 = pbuffer.data(idx_npot_0_if + 17);

    auto ta_xxxxxy_yzz_0 = pbuffer.data(idx_npot_0_if + 18);

    auto ta_xxxxxy_zzz_0 = pbuffer.data(idx_npot_0_if + 19);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxxxx_xx_0,   \
                             ta_xxxxx_xx_1,   \
                             ta_xxxxx_xxx_0,  \
                             ta_xxxxx_xxx_1,  \
                             ta_xxxxx_xxy_0,  \
                             ta_xxxxx_xxy_1,  \
                             ta_xxxxx_xxz_0,  \
                             ta_xxxxx_xxz_1,  \
                             ta_xxxxx_xy_0,   \
                             ta_xxxxx_xy_1,   \
                             ta_xxxxx_xyy_0,  \
                             ta_xxxxx_xyy_1,  \
                             ta_xxxxx_xyz_0,  \
                             ta_xxxxx_xyz_1,  \
                             ta_xxxxx_xz_0,   \
                             ta_xxxxx_xz_1,   \
                             ta_xxxxx_xzz_0,  \
                             ta_xxxxx_xzz_1,  \
                             ta_xxxxx_zzz_0,  \
                             ta_xxxxx_zzz_1,  \
                             ta_xxxxxy_xxx_0, \
                             ta_xxxxxy_xxy_0, \
                             ta_xxxxxy_xxz_0, \
                             ta_xxxxxy_xyy_0, \
                             ta_xxxxxy_xyz_0, \
                             ta_xxxxxy_xzz_0, \
                             ta_xxxxxy_yyy_0, \
                             ta_xxxxxy_yyz_0, \
                             ta_xxxxxy_yzz_0, \
                             ta_xxxxxy_zzz_0, \
                             ta_xxxxy_yyy_0,  \
                             ta_xxxxy_yyy_1,  \
                             ta_xxxxy_yyz_0,  \
                             ta_xxxxy_yyz_1,  \
                             ta_xxxxy_yzz_0,  \
                             ta_xxxxy_yzz_1,  \
                             ta_xxxy_yyy_0,   \
                             ta_xxxy_yyy_1,   \
                             ta_xxxy_yyz_0,   \
                             ta_xxxy_yyz_1,   \
                             ta_xxxy_yzz_0,   \
                             ta_xxxy_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_xxx_0[i] = ta_xxxxx_xxx_0[i] * pa_y[i] - ta_xxxxx_xxx_1[i] * pc_y[i];

        ta_xxxxxy_xxy_0[i] = ta_xxxxx_xx_0[i] * fe_0 - ta_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxy_0[i] * pa_y[i] - ta_xxxxx_xxy_1[i] * pc_y[i];

        ta_xxxxxy_xxz_0[i] = ta_xxxxx_xxz_0[i] * pa_y[i] - ta_xxxxx_xxz_1[i] * pc_y[i];

        ta_xxxxxy_xyy_0[i] =
            2.0 * ta_xxxxx_xy_0[i] * fe_0 - 2.0 * ta_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xyy_0[i] * pa_y[i] - ta_xxxxx_xyy_1[i] * pc_y[i];

        ta_xxxxxy_xyz_0[i] = ta_xxxxx_xz_0[i] * fe_0 - ta_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xyz_0[i] * pa_y[i] - ta_xxxxx_xyz_1[i] * pc_y[i];

        ta_xxxxxy_xzz_0[i] = ta_xxxxx_xzz_0[i] * pa_y[i] - ta_xxxxx_xzz_1[i] * pc_y[i];

        ta_xxxxxy_yyy_0[i] =
            4.0 * ta_xxxy_yyy_0[i] * fe_0 - 4.0 * ta_xxxy_yyy_1[i] * fe_0 + ta_xxxxy_yyy_0[i] * pa_x[i] - ta_xxxxy_yyy_1[i] * pc_x[i];

        ta_xxxxxy_yyz_0[i] =
            4.0 * ta_xxxy_yyz_0[i] * fe_0 - 4.0 * ta_xxxy_yyz_1[i] * fe_0 + ta_xxxxy_yyz_0[i] * pa_x[i] - ta_xxxxy_yyz_1[i] * pc_x[i];

        ta_xxxxxy_yzz_0[i] =
            4.0 * ta_xxxy_yzz_0[i] * fe_0 - 4.0 * ta_xxxy_yzz_1[i] * fe_0 + ta_xxxxy_yzz_0[i] * pa_x[i] - ta_xxxxy_yzz_1[i] * pc_x[i];

        ta_xxxxxy_zzz_0[i] = ta_xxxxx_zzz_0[i] * pa_y[i] - ta_xxxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : IF

    auto ta_xxxxxz_xxx_0 = pbuffer.data(idx_npot_0_if + 20);

    auto ta_xxxxxz_xxy_0 = pbuffer.data(idx_npot_0_if + 21);

    auto ta_xxxxxz_xxz_0 = pbuffer.data(idx_npot_0_if + 22);

    auto ta_xxxxxz_xyy_0 = pbuffer.data(idx_npot_0_if + 23);

    auto ta_xxxxxz_xyz_0 = pbuffer.data(idx_npot_0_if + 24);

    auto ta_xxxxxz_xzz_0 = pbuffer.data(idx_npot_0_if + 25);

    auto ta_xxxxxz_yyy_0 = pbuffer.data(idx_npot_0_if + 26);

    auto ta_xxxxxz_yyz_0 = pbuffer.data(idx_npot_0_if + 27);

    auto ta_xxxxxz_yzz_0 = pbuffer.data(idx_npot_0_if + 28);

    auto ta_xxxxxz_zzz_0 = pbuffer.data(idx_npot_0_if + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxxxx_xx_0,   \
                             ta_xxxxx_xx_1,   \
                             ta_xxxxx_xxx_0,  \
                             ta_xxxxx_xxx_1,  \
                             ta_xxxxx_xxy_0,  \
                             ta_xxxxx_xxy_1,  \
                             ta_xxxxx_xxz_0,  \
                             ta_xxxxx_xxz_1,  \
                             ta_xxxxx_xy_0,   \
                             ta_xxxxx_xy_1,   \
                             ta_xxxxx_xyy_0,  \
                             ta_xxxxx_xyy_1,  \
                             ta_xxxxx_xyz_0,  \
                             ta_xxxxx_xyz_1,  \
                             ta_xxxxx_xz_0,   \
                             ta_xxxxx_xz_1,   \
                             ta_xxxxx_xzz_0,  \
                             ta_xxxxx_xzz_1,  \
                             ta_xxxxx_yyy_0,  \
                             ta_xxxxx_yyy_1,  \
                             ta_xxxxxz_xxx_0, \
                             ta_xxxxxz_xxy_0, \
                             ta_xxxxxz_xxz_0, \
                             ta_xxxxxz_xyy_0, \
                             ta_xxxxxz_xyz_0, \
                             ta_xxxxxz_xzz_0, \
                             ta_xxxxxz_yyy_0, \
                             ta_xxxxxz_yyz_0, \
                             ta_xxxxxz_yzz_0, \
                             ta_xxxxxz_zzz_0, \
                             ta_xxxxz_yyz_0,  \
                             ta_xxxxz_yyz_1,  \
                             ta_xxxxz_yzz_0,  \
                             ta_xxxxz_yzz_1,  \
                             ta_xxxxz_zzz_0,  \
                             ta_xxxxz_zzz_1,  \
                             ta_xxxz_yyz_0,   \
                             ta_xxxz_yyz_1,   \
                             ta_xxxz_yzz_0,   \
                             ta_xxxz_yzz_1,   \
                             ta_xxxz_zzz_0,   \
                             ta_xxxz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_xxx_0[i] = ta_xxxxx_xxx_0[i] * pa_z[i] - ta_xxxxx_xxx_1[i] * pc_z[i];

        ta_xxxxxz_xxy_0[i] = ta_xxxxx_xxy_0[i] * pa_z[i] - ta_xxxxx_xxy_1[i] * pc_z[i];

        ta_xxxxxz_xxz_0[i] = ta_xxxxx_xx_0[i] * fe_0 - ta_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxz_0[i] * pa_z[i] - ta_xxxxx_xxz_1[i] * pc_z[i];

        ta_xxxxxz_xyy_0[i] = ta_xxxxx_xyy_0[i] * pa_z[i] - ta_xxxxx_xyy_1[i] * pc_z[i];

        ta_xxxxxz_xyz_0[i] = ta_xxxxx_xy_0[i] * fe_0 - ta_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xyz_0[i] * pa_z[i] - ta_xxxxx_xyz_1[i] * pc_z[i];

        ta_xxxxxz_xzz_0[i] =
            2.0 * ta_xxxxx_xz_0[i] * fe_0 - 2.0 * ta_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xzz_0[i] * pa_z[i] - ta_xxxxx_xzz_1[i] * pc_z[i];

        ta_xxxxxz_yyy_0[i] = ta_xxxxx_yyy_0[i] * pa_z[i] - ta_xxxxx_yyy_1[i] * pc_z[i];

        ta_xxxxxz_yyz_0[i] =
            4.0 * ta_xxxz_yyz_0[i] * fe_0 - 4.0 * ta_xxxz_yyz_1[i] * fe_0 + ta_xxxxz_yyz_0[i] * pa_x[i] - ta_xxxxz_yyz_1[i] * pc_x[i];

        ta_xxxxxz_yzz_0[i] =
            4.0 * ta_xxxz_yzz_0[i] * fe_0 - 4.0 * ta_xxxz_yzz_1[i] * fe_0 + ta_xxxxz_yzz_0[i] * pa_x[i] - ta_xxxxz_yzz_1[i] * pc_x[i];

        ta_xxxxxz_zzz_0[i] =
            4.0 * ta_xxxz_zzz_0[i] * fe_0 - 4.0 * ta_xxxz_zzz_1[i] * fe_0 + ta_xxxxz_zzz_0[i] * pa_x[i] - ta_xxxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : IF

    auto ta_xxxxyy_xxx_0 = pbuffer.data(idx_npot_0_if + 30);

    auto ta_xxxxyy_xxy_0 = pbuffer.data(idx_npot_0_if + 31);

    auto ta_xxxxyy_xxz_0 = pbuffer.data(idx_npot_0_if + 32);

    auto ta_xxxxyy_xyy_0 = pbuffer.data(idx_npot_0_if + 33);

    auto ta_xxxxyy_xyz_0 = pbuffer.data(idx_npot_0_if + 34);

    auto ta_xxxxyy_xzz_0 = pbuffer.data(idx_npot_0_if + 35);

    auto ta_xxxxyy_yyy_0 = pbuffer.data(idx_npot_0_if + 36);

    auto ta_xxxxyy_yyz_0 = pbuffer.data(idx_npot_0_if + 37);

    auto ta_xxxxyy_yzz_0 = pbuffer.data(idx_npot_0_if + 38);

    auto ta_xxxxyy_zzz_0 = pbuffer.data(idx_npot_0_if + 39);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxz_0,   \
                             ta_xxxx_xxz_1,   \
                             ta_xxxx_xzz_0,   \
                             ta_xxxx_xzz_1,   \
                             ta_xxxxy_xxx_0,  \
                             ta_xxxxy_xxx_1,  \
                             ta_xxxxy_xxz_0,  \
                             ta_xxxxy_xxz_1,  \
                             ta_xxxxy_xzz_0,  \
                             ta_xxxxy_xzz_1,  \
                             ta_xxxxyy_xxx_0, \
                             ta_xxxxyy_xxy_0, \
                             ta_xxxxyy_xxz_0, \
                             ta_xxxxyy_xyy_0, \
                             ta_xxxxyy_xyz_0, \
                             ta_xxxxyy_xzz_0, \
                             ta_xxxxyy_yyy_0, \
                             ta_xxxxyy_yyz_0, \
                             ta_xxxxyy_yzz_0, \
                             ta_xxxxyy_zzz_0, \
                             ta_xxxyy_xxy_0,  \
                             ta_xxxyy_xxy_1,  \
                             ta_xxxyy_xy_0,   \
                             ta_xxxyy_xy_1,   \
                             ta_xxxyy_xyy_0,  \
                             ta_xxxyy_xyy_1,  \
                             ta_xxxyy_xyz_0,  \
                             ta_xxxyy_xyz_1,  \
                             ta_xxxyy_yy_0,   \
                             ta_xxxyy_yy_1,   \
                             ta_xxxyy_yyy_0,  \
                             ta_xxxyy_yyy_1,  \
                             ta_xxxyy_yyz_0,  \
                             ta_xxxyy_yyz_1,  \
                             ta_xxxyy_yz_0,   \
                             ta_xxxyy_yz_1,   \
                             ta_xxxyy_yzz_0,  \
                             ta_xxxyy_yzz_1,  \
                             ta_xxxyy_zzz_0,  \
                             ta_xxxyy_zzz_1,  \
                             ta_xxyy_xxy_0,   \
                             ta_xxyy_xxy_1,   \
                             ta_xxyy_xyy_0,   \
                             ta_xxyy_xyy_1,   \
                             ta_xxyy_xyz_0,   \
                             ta_xxyy_xyz_1,   \
                             ta_xxyy_yyy_0,   \
                             ta_xxyy_yyy_1,   \
                             ta_xxyy_yyz_0,   \
                             ta_xxyy_yyz_1,   \
                             ta_xxyy_yzz_0,   \
                             ta_xxyy_yzz_1,   \
                             ta_xxyy_zzz_0,   \
                             ta_xxyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_xxx_0[i] = ta_xxxx_xxx_0[i] * fe_0 - ta_xxxx_xxx_1[i] * fe_0 + ta_xxxxy_xxx_0[i] * pa_y[i] - ta_xxxxy_xxx_1[i] * pc_y[i];

        ta_xxxxyy_xxy_0[i] = 3.0 * ta_xxyy_xxy_0[i] * fe_0 - 3.0 * ta_xxyy_xxy_1[i] * fe_0 + 2.0 * ta_xxxyy_xy_0[i] * fe_0 -
                             2.0 * ta_xxxyy_xy_1[i] * fe_0 + ta_xxxyy_xxy_0[i] * pa_x[i] - ta_xxxyy_xxy_1[i] * pc_x[i];

        ta_xxxxyy_xxz_0[i] = ta_xxxx_xxz_0[i] * fe_0 - ta_xxxx_xxz_1[i] * fe_0 + ta_xxxxy_xxz_0[i] * pa_y[i] - ta_xxxxy_xxz_1[i] * pc_y[i];

        ta_xxxxyy_xyy_0[i] = 3.0 * ta_xxyy_xyy_0[i] * fe_0 - 3.0 * ta_xxyy_xyy_1[i] * fe_0 + ta_xxxyy_yy_0[i] * fe_0 - ta_xxxyy_yy_1[i] * fe_0 +
                             ta_xxxyy_xyy_0[i] * pa_x[i] - ta_xxxyy_xyy_1[i] * pc_x[i];

        ta_xxxxyy_xyz_0[i] = 3.0 * ta_xxyy_xyz_0[i] * fe_0 - 3.0 * ta_xxyy_xyz_1[i] * fe_0 + ta_xxxyy_yz_0[i] * fe_0 - ta_xxxyy_yz_1[i] * fe_0 +
                             ta_xxxyy_xyz_0[i] * pa_x[i] - ta_xxxyy_xyz_1[i] * pc_x[i];

        ta_xxxxyy_xzz_0[i] = ta_xxxx_xzz_0[i] * fe_0 - ta_xxxx_xzz_1[i] * fe_0 + ta_xxxxy_xzz_0[i] * pa_y[i] - ta_xxxxy_xzz_1[i] * pc_y[i];

        ta_xxxxyy_yyy_0[i] =
            3.0 * ta_xxyy_yyy_0[i] * fe_0 - 3.0 * ta_xxyy_yyy_1[i] * fe_0 + ta_xxxyy_yyy_0[i] * pa_x[i] - ta_xxxyy_yyy_1[i] * pc_x[i];

        ta_xxxxyy_yyz_0[i] =
            3.0 * ta_xxyy_yyz_0[i] * fe_0 - 3.0 * ta_xxyy_yyz_1[i] * fe_0 + ta_xxxyy_yyz_0[i] * pa_x[i] - ta_xxxyy_yyz_1[i] * pc_x[i];

        ta_xxxxyy_yzz_0[i] =
            3.0 * ta_xxyy_yzz_0[i] * fe_0 - 3.0 * ta_xxyy_yzz_1[i] * fe_0 + ta_xxxyy_yzz_0[i] * pa_x[i] - ta_xxxyy_yzz_1[i] * pc_x[i];

        ta_xxxxyy_zzz_0[i] =
            3.0 * ta_xxyy_zzz_0[i] * fe_0 - 3.0 * ta_xxyy_zzz_1[i] * fe_0 + ta_xxxyy_zzz_0[i] * pa_x[i] - ta_xxxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : IF

    auto ta_xxxxyz_xxx_0 = pbuffer.data(idx_npot_0_if + 40);

    auto ta_xxxxyz_xxy_0 = pbuffer.data(idx_npot_0_if + 41);

    auto ta_xxxxyz_xxz_0 = pbuffer.data(idx_npot_0_if + 42);

    auto ta_xxxxyz_xyy_0 = pbuffer.data(idx_npot_0_if + 43);

    auto ta_xxxxyz_xyz_0 = pbuffer.data(idx_npot_0_if + 44);

    auto ta_xxxxyz_xzz_0 = pbuffer.data(idx_npot_0_if + 45);

    auto ta_xxxxyz_yyy_0 = pbuffer.data(idx_npot_0_if + 46);

    auto ta_xxxxyz_yyz_0 = pbuffer.data(idx_npot_0_if + 47);

    auto ta_xxxxyz_yzz_0 = pbuffer.data(idx_npot_0_if + 48);

    auto ta_xxxxyz_zzz_0 = pbuffer.data(idx_npot_0_if + 49);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxxxy_xxy_0,  \
                             ta_xxxxy_xxy_1,  \
                             ta_xxxxy_xyy_0,  \
                             ta_xxxxy_xyy_1,  \
                             ta_xxxxy_yyy_0,  \
                             ta_xxxxy_yyy_1,  \
                             ta_xxxxyz_xxx_0, \
                             ta_xxxxyz_xxy_0, \
                             ta_xxxxyz_xxz_0, \
                             ta_xxxxyz_xyy_0, \
                             ta_xxxxyz_xyz_0, \
                             ta_xxxxyz_xzz_0, \
                             ta_xxxxyz_yyy_0, \
                             ta_xxxxyz_yyz_0, \
                             ta_xxxxyz_yzz_0, \
                             ta_xxxxyz_zzz_0, \
                             ta_xxxxz_xxx_0,  \
                             ta_xxxxz_xxx_1,  \
                             ta_xxxxz_xxz_0,  \
                             ta_xxxxz_xxz_1,  \
                             ta_xxxxz_xyz_0,  \
                             ta_xxxxz_xyz_1,  \
                             ta_xxxxz_xz_0,   \
                             ta_xxxxz_xz_1,   \
                             ta_xxxxz_xzz_0,  \
                             ta_xxxxz_xzz_1,  \
                             ta_xxxxz_zzz_0,  \
                             ta_xxxxz_zzz_1,  \
                             ta_xxxyz_yyz_0,  \
                             ta_xxxyz_yyz_1,  \
                             ta_xxxyz_yzz_0,  \
                             ta_xxxyz_yzz_1,  \
                             ta_xxyz_yyz_0,   \
                             ta_xxyz_yyz_1,   \
                             ta_xxyz_yzz_0,   \
                             ta_xxyz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyz_xxx_0[i] = ta_xxxxz_xxx_0[i] * pa_y[i] - ta_xxxxz_xxx_1[i] * pc_y[i];

        ta_xxxxyz_xxy_0[i] = ta_xxxxy_xxy_0[i] * pa_z[i] - ta_xxxxy_xxy_1[i] * pc_z[i];

        ta_xxxxyz_xxz_0[i] = ta_xxxxz_xxz_0[i] * pa_y[i] - ta_xxxxz_xxz_1[i] * pc_y[i];

        ta_xxxxyz_xyy_0[i] = ta_xxxxy_xyy_0[i] * pa_z[i] - ta_xxxxy_xyy_1[i] * pc_z[i];

        ta_xxxxyz_xyz_0[i] = ta_xxxxz_xz_0[i] * fe_0 - ta_xxxxz_xz_1[i] * fe_0 + ta_xxxxz_xyz_0[i] * pa_y[i] - ta_xxxxz_xyz_1[i] * pc_y[i];

        ta_xxxxyz_xzz_0[i] = ta_xxxxz_xzz_0[i] * pa_y[i] - ta_xxxxz_xzz_1[i] * pc_y[i];

        ta_xxxxyz_yyy_0[i] = ta_xxxxy_yyy_0[i] * pa_z[i] - ta_xxxxy_yyy_1[i] * pc_z[i];

        ta_xxxxyz_yyz_0[i] =
            3.0 * ta_xxyz_yyz_0[i] * fe_0 - 3.0 * ta_xxyz_yyz_1[i] * fe_0 + ta_xxxyz_yyz_0[i] * pa_x[i] - ta_xxxyz_yyz_1[i] * pc_x[i];

        ta_xxxxyz_yzz_0[i] =
            3.0 * ta_xxyz_yzz_0[i] * fe_0 - 3.0 * ta_xxyz_yzz_1[i] * fe_0 + ta_xxxyz_yzz_0[i] * pa_x[i] - ta_xxxyz_yzz_1[i] * pc_x[i];

        ta_xxxxyz_zzz_0[i] = ta_xxxxz_zzz_0[i] * pa_y[i] - ta_xxxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : IF

    auto ta_xxxxzz_xxx_0 = pbuffer.data(idx_npot_0_if + 50);

    auto ta_xxxxzz_xxy_0 = pbuffer.data(idx_npot_0_if + 51);

    auto ta_xxxxzz_xxz_0 = pbuffer.data(idx_npot_0_if + 52);

    auto ta_xxxxzz_xyy_0 = pbuffer.data(idx_npot_0_if + 53);

    auto ta_xxxxzz_xyz_0 = pbuffer.data(idx_npot_0_if + 54);

    auto ta_xxxxzz_xzz_0 = pbuffer.data(idx_npot_0_if + 55);

    auto ta_xxxxzz_yyy_0 = pbuffer.data(idx_npot_0_if + 56);

    auto ta_xxxxzz_yyz_0 = pbuffer.data(idx_npot_0_if + 57);

    auto ta_xxxxzz_yzz_0 = pbuffer.data(idx_npot_0_if + 58);

    auto ta_xxxxzz_zzz_0 = pbuffer.data(idx_npot_0_if + 59);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxy_0,   \
                             ta_xxxx_xxy_1,   \
                             ta_xxxx_xyy_0,   \
                             ta_xxxx_xyy_1,   \
                             ta_xxxxz_xxx_0,  \
                             ta_xxxxz_xxx_1,  \
                             ta_xxxxz_xxy_0,  \
                             ta_xxxxz_xxy_1,  \
                             ta_xxxxz_xyy_0,  \
                             ta_xxxxz_xyy_1,  \
                             ta_xxxxzz_xxx_0, \
                             ta_xxxxzz_xxy_0, \
                             ta_xxxxzz_xxz_0, \
                             ta_xxxxzz_xyy_0, \
                             ta_xxxxzz_xyz_0, \
                             ta_xxxxzz_xzz_0, \
                             ta_xxxxzz_yyy_0, \
                             ta_xxxxzz_yyz_0, \
                             ta_xxxxzz_yzz_0, \
                             ta_xxxxzz_zzz_0, \
                             ta_xxxzz_xxz_0,  \
                             ta_xxxzz_xxz_1,  \
                             ta_xxxzz_xyz_0,  \
                             ta_xxxzz_xyz_1,  \
                             ta_xxxzz_xz_0,   \
                             ta_xxxzz_xz_1,   \
                             ta_xxxzz_xzz_0,  \
                             ta_xxxzz_xzz_1,  \
                             ta_xxxzz_yyy_0,  \
                             ta_xxxzz_yyy_1,  \
                             ta_xxxzz_yyz_0,  \
                             ta_xxxzz_yyz_1,  \
                             ta_xxxzz_yz_0,   \
                             ta_xxxzz_yz_1,   \
                             ta_xxxzz_yzz_0,  \
                             ta_xxxzz_yzz_1,  \
                             ta_xxxzz_zz_0,   \
                             ta_xxxzz_zz_1,   \
                             ta_xxxzz_zzz_0,  \
                             ta_xxxzz_zzz_1,  \
                             ta_xxzz_xxz_0,   \
                             ta_xxzz_xxz_1,   \
                             ta_xxzz_xyz_0,   \
                             ta_xxzz_xyz_1,   \
                             ta_xxzz_xzz_0,   \
                             ta_xxzz_xzz_1,   \
                             ta_xxzz_yyy_0,   \
                             ta_xxzz_yyy_1,   \
                             ta_xxzz_yyz_0,   \
                             ta_xxzz_yyz_1,   \
                             ta_xxzz_yzz_0,   \
                             ta_xxzz_yzz_1,   \
                             ta_xxzz_zzz_0,   \
                             ta_xxzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_xxx_0[i] = ta_xxxx_xxx_0[i] * fe_0 - ta_xxxx_xxx_1[i] * fe_0 + ta_xxxxz_xxx_0[i] * pa_z[i] - ta_xxxxz_xxx_1[i] * pc_z[i];

        ta_xxxxzz_xxy_0[i] = ta_xxxx_xxy_0[i] * fe_0 - ta_xxxx_xxy_1[i] * fe_0 + ta_xxxxz_xxy_0[i] * pa_z[i] - ta_xxxxz_xxy_1[i] * pc_z[i];

        ta_xxxxzz_xxz_0[i] = 3.0 * ta_xxzz_xxz_0[i] * fe_0 - 3.0 * ta_xxzz_xxz_1[i] * fe_0 + 2.0 * ta_xxxzz_xz_0[i] * fe_0 -
                             2.0 * ta_xxxzz_xz_1[i] * fe_0 + ta_xxxzz_xxz_0[i] * pa_x[i] - ta_xxxzz_xxz_1[i] * pc_x[i];

        ta_xxxxzz_xyy_0[i] = ta_xxxx_xyy_0[i] * fe_0 - ta_xxxx_xyy_1[i] * fe_0 + ta_xxxxz_xyy_0[i] * pa_z[i] - ta_xxxxz_xyy_1[i] * pc_z[i];

        ta_xxxxzz_xyz_0[i] = 3.0 * ta_xxzz_xyz_0[i] * fe_0 - 3.0 * ta_xxzz_xyz_1[i] * fe_0 + ta_xxxzz_yz_0[i] * fe_0 - ta_xxxzz_yz_1[i] * fe_0 +
                             ta_xxxzz_xyz_0[i] * pa_x[i] - ta_xxxzz_xyz_1[i] * pc_x[i];

        ta_xxxxzz_xzz_0[i] = 3.0 * ta_xxzz_xzz_0[i] * fe_0 - 3.0 * ta_xxzz_xzz_1[i] * fe_0 + ta_xxxzz_zz_0[i] * fe_0 - ta_xxxzz_zz_1[i] * fe_0 +
                             ta_xxxzz_xzz_0[i] * pa_x[i] - ta_xxxzz_xzz_1[i] * pc_x[i];

        ta_xxxxzz_yyy_0[i] =
            3.0 * ta_xxzz_yyy_0[i] * fe_0 - 3.0 * ta_xxzz_yyy_1[i] * fe_0 + ta_xxxzz_yyy_0[i] * pa_x[i] - ta_xxxzz_yyy_1[i] * pc_x[i];

        ta_xxxxzz_yyz_0[i] =
            3.0 * ta_xxzz_yyz_0[i] * fe_0 - 3.0 * ta_xxzz_yyz_1[i] * fe_0 + ta_xxxzz_yyz_0[i] * pa_x[i] - ta_xxxzz_yyz_1[i] * pc_x[i];

        ta_xxxxzz_yzz_0[i] =
            3.0 * ta_xxzz_yzz_0[i] * fe_0 - 3.0 * ta_xxzz_yzz_1[i] * fe_0 + ta_xxxzz_yzz_0[i] * pa_x[i] - ta_xxxzz_yzz_1[i] * pc_x[i];

        ta_xxxxzz_zzz_0[i] =
            3.0 * ta_xxzz_zzz_0[i] * fe_0 - 3.0 * ta_xxzz_zzz_1[i] * fe_0 + ta_xxxzz_zzz_0[i] * pa_x[i] - ta_xxxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : IF

    auto ta_xxxyyy_xxx_0 = pbuffer.data(idx_npot_0_if + 60);

    auto ta_xxxyyy_xxy_0 = pbuffer.data(idx_npot_0_if + 61);

    auto ta_xxxyyy_xxz_0 = pbuffer.data(idx_npot_0_if + 62);

    auto ta_xxxyyy_xyy_0 = pbuffer.data(idx_npot_0_if + 63);

    auto ta_xxxyyy_xyz_0 = pbuffer.data(idx_npot_0_if + 64);

    auto ta_xxxyyy_xzz_0 = pbuffer.data(idx_npot_0_if + 65);

    auto ta_xxxyyy_yyy_0 = pbuffer.data(idx_npot_0_if + 66);

    auto ta_xxxyyy_yyz_0 = pbuffer.data(idx_npot_0_if + 67);

    auto ta_xxxyyy_yzz_0 = pbuffer.data(idx_npot_0_if + 68);

    auto ta_xxxyyy_zzz_0 = pbuffer.data(idx_npot_0_if + 69);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxxy_xxx_0,   \
                             ta_xxxy_xxx_1,   \
                             ta_xxxy_xxz_0,   \
                             ta_xxxy_xxz_1,   \
                             ta_xxxy_xzz_0,   \
                             ta_xxxy_xzz_1,   \
                             ta_xxxyy_xxx_0,  \
                             ta_xxxyy_xxx_1,  \
                             ta_xxxyy_xxz_0,  \
                             ta_xxxyy_xxz_1,  \
                             ta_xxxyy_xzz_0,  \
                             ta_xxxyy_xzz_1,  \
                             ta_xxxyyy_xxx_0, \
                             ta_xxxyyy_xxy_0, \
                             ta_xxxyyy_xxz_0, \
                             ta_xxxyyy_xyy_0, \
                             ta_xxxyyy_xyz_0, \
                             ta_xxxyyy_xzz_0, \
                             ta_xxxyyy_yyy_0, \
                             ta_xxxyyy_yyz_0, \
                             ta_xxxyyy_yzz_0, \
                             ta_xxxyyy_zzz_0, \
                             ta_xxyyy_xxy_0,  \
                             ta_xxyyy_xxy_1,  \
                             ta_xxyyy_xy_0,   \
                             ta_xxyyy_xy_1,   \
                             ta_xxyyy_xyy_0,  \
                             ta_xxyyy_xyy_1,  \
                             ta_xxyyy_xyz_0,  \
                             ta_xxyyy_xyz_1,  \
                             ta_xxyyy_yy_0,   \
                             ta_xxyyy_yy_1,   \
                             ta_xxyyy_yyy_0,  \
                             ta_xxyyy_yyy_1,  \
                             ta_xxyyy_yyz_0,  \
                             ta_xxyyy_yyz_1,  \
                             ta_xxyyy_yz_0,   \
                             ta_xxyyy_yz_1,   \
                             ta_xxyyy_yzz_0,  \
                             ta_xxyyy_yzz_1,  \
                             ta_xxyyy_zzz_0,  \
                             ta_xxyyy_zzz_1,  \
                             ta_xyyy_xxy_0,   \
                             ta_xyyy_xxy_1,   \
                             ta_xyyy_xyy_0,   \
                             ta_xyyy_xyy_1,   \
                             ta_xyyy_xyz_0,   \
                             ta_xyyy_xyz_1,   \
                             ta_xyyy_yyy_0,   \
                             ta_xyyy_yyy_1,   \
                             ta_xyyy_yyz_0,   \
                             ta_xyyy_yyz_1,   \
                             ta_xyyy_yzz_0,   \
                             ta_xyyy_yzz_1,   \
                             ta_xyyy_zzz_0,   \
                             ta_xyyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_xxx_0[i] =
            2.0 * ta_xxxy_xxx_0[i] * fe_0 - 2.0 * ta_xxxy_xxx_1[i] * fe_0 + ta_xxxyy_xxx_0[i] * pa_y[i] - ta_xxxyy_xxx_1[i] * pc_y[i];

        ta_xxxyyy_xxy_0[i] = 2.0 * ta_xyyy_xxy_0[i] * fe_0 - 2.0 * ta_xyyy_xxy_1[i] * fe_0 + 2.0 * ta_xxyyy_xy_0[i] * fe_0 -
                             2.0 * ta_xxyyy_xy_1[i] * fe_0 + ta_xxyyy_xxy_0[i] * pa_x[i] - ta_xxyyy_xxy_1[i] * pc_x[i];

        ta_xxxyyy_xxz_0[i] =
            2.0 * ta_xxxy_xxz_0[i] * fe_0 - 2.0 * ta_xxxy_xxz_1[i] * fe_0 + ta_xxxyy_xxz_0[i] * pa_y[i] - ta_xxxyy_xxz_1[i] * pc_y[i];

        ta_xxxyyy_xyy_0[i] = 2.0 * ta_xyyy_xyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyy_1[i] * fe_0 + ta_xxyyy_yy_0[i] * fe_0 - ta_xxyyy_yy_1[i] * fe_0 +
                             ta_xxyyy_xyy_0[i] * pa_x[i] - ta_xxyyy_xyy_1[i] * pc_x[i];

        ta_xxxyyy_xyz_0[i] = 2.0 * ta_xyyy_xyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyz_1[i] * fe_0 + ta_xxyyy_yz_0[i] * fe_0 - ta_xxyyy_yz_1[i] * fe_0 +
                             ta_xxyyy_xyz_0[i] * pa_x[i] - ta_xxyyy_xyz_1[i] * pc_x[i];

        ta_xxxyyy_xzz_0[i] =
            2.0 * ta_xxxy_xzz_0[i] * fe_0 - 2.0 * ta_xxxy_xzz_1[i] * fe_0 + ta_xxxyy_xzz_0[i] * pa_y[i] - ta_xxxyy_xzz_1[i] * pc_y[i];

        ta_xxxyyy_yyy_0[i] =
            2.0 * ta_xyyy_yyy_0[i] * fe_0 - 2.0 * ta_xyyy_yyy_1[i] * fe_0 + ta_xxyyy_yyy_0[i] * pa_x[i] - ta_xxyyy_yyy_1[i] * pc_x[i];

        ta_xxxyyy_yyz_0[i] =
            2.0 * ta_xyyy_yyz_0[i] * fe_0 - 2.0 * ta_xyyy_yyz_1[i] * fe_0 + ta_xxyyy_yyz_0[i] * pa_x[i] - ta_xxyyy_yyz_1[i] * pc_x[i];

        ta_xxxyyy_yzz_0[i] =
            2.0 * ta_xyyy_yzz_0[i] * fe_0 - 2.0 * ta_xyyy_yzz_1[i] * fe_0 + ta_xxyyy_yzz_0[i] * pa_x[i] - ta_xxyyy_yzz_1[i] * pc_x[i];

        ta_xxxyyy_zzz_0[i] =
            2.0 * ta_xyyy_zzz_0[i] * fe_0 - 2.0 * ta_xyyy_zzz_1[i] * fe_0 + ta_xxyyy_zzz_0[i] * pa_x[i] - ta_xxyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : IF

    auto ta_xxxyyz_xxx_0 = pbuffer.data(idx_npot_0_if + 70);

    auto ta_xxxyyz_xxy_0 = pbuffer.data(idx_npot_0_if + 71);

    auto ta_xxxyyz_xxz_0 = pbuffer.data(idx_npot_0_if + 72);

    auto ta_xxxyyz_xyy_0 = pbuffer.data(idx_npot_0_if + 73);

    auto ta_xxxyyz_xyz_0 = pbuffer.data(idx_npot_0_if + 74);

    auto ta_xxxyyz_xzz_0 = pbuffer.data(idx_npot_0_if + 75);

    auto ta_xxxyyz_yyy_0 = pbuffer.data(idx_npot_0_if + 76);

    auto ta_xxxyyz_yyz_0 = pbuffer.data(idx_npot_0_if + 77);

    auto ta_xxxyyz_yzz_0 = pbuffer.data(idx_npot_0_if + 78);

    auto ta_xxxyyz_zzz_0 = pbuffer.data(idx_npot_0_if + 79);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxxyy_xxx_0,  \
                             ta_xxxyy_xxx_1,  \
                             ta_xxxyy_xxy_0,  \
                             ta_xxxyy_xxy_1,  \
                             ta_xxxyy_xy_0,   \
                             ta_xxxyy_xy_1,   \
                             ta_xxxyy_xyy_0,  \
                             ta_xxxyy_xyy_1,  \
                             ta_xxxyy_xyz_0,  \
                             ta_xxxyy_xyz_1,  \
                             ta_xxxyy_yyy_0,  \
                             ta_xxxyy_yyy_1,  \
                             ta_xxxyyz_xxx_0, \
                             ta_xxxyyz_xxy_0, \
                             ta_xxxyyz_xxz_0, \
                             ta_xxxyyz_xyy_0, \
                             ta_xxxyyz_xyz_0, \
                             ta_xxxyyz_xzz_0, \
                             ta_xxxyyz_yyy_0, \
                             ta_xxxyyz_yyz_0, \
                             ta_xxxyyz_yzz_0, \
                             ta_xxxyyz_zzz_0, \
                             ta_xxxyz_xxz_0,  \
                             ta_xxxyz_xxz_1,  \
                             ta_xxxyz_xzz_0,  \
                             ta_xxxyz_xzz_1,  \
                             ta_xxxz_xxz_0,   \
                             ta_xxxz_xxz_1,   \
                             ta_xxxz_xzz_0,   \
                             ta_xxxz_xzz_1,   \
                             ta_xxyyz_yyz_0,  \
                             ta_xxyyz_yyz_1,  \
                             ta_xxyyz_yzz_0,  \
                             ta_xxyyz_yzz_1,  \
                             ta_xxyyz_zzz_0,  \
                             ta_xxyyz_zzz_1,  \
                             ta_xyyz_yyz_0,   \
                             ta_xyyz_yyz_1,   \
                             ta_xyyz_yzz_0,   \
                             ta_xyyz_yzz_1,   \
                             ta_xyyz_zzz_0,   \
                             ta_xyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_xxx_0[i] = ta_xxxyy_xxx_0[i] * pa_z[i] - ta_xxxyy_xxx_1[i] * pc_z[i];

        ta_xxxyyz_xxy_0[i] = ta_xxxyy_xxy_0[i] * pa_z[i] - ta_xxxyy_xxy_1[i] * pc_z[i];

        ta_xxxyyz_xxz_0[i] = ta_xxxz_xxz_0[i] * fe_0 - ta_xxxz_xxz_1[i] * fe_0 + ta_xxxyz_xxz_0[i] * pa_y[i] - ta_xxxyz_xxz_1[i] * pc_y[i];

        ta_xxxyyz_xyy_0[i] = ta_xxxyy_xyy_0[i] * pa_z[i] - ta_xxxyy_xyy_1[i] * pc_z[i];

        ta_xxxyyz_xyz_0[i] = ta_xxxyy_xy_0[i] * fe_0 - ta_xxxyy_xy_1[i] * fe_0 + ta_xxxyy_xyz_0[i] * pa_z[i] - ta_xxxyy_xyz_1[i] * pc_z[i];

        ta_xxxyyz_xzz_0[i] = ta_xxxz_xzz_0[i] * fe_0 - ta_xxxz_xzz_1[i] * fe_0 + ta_xxxyz_xzz_0[i] * pa_y[i] - ta_xxxyz_xzz_1[i] * pc_y[i];

        ta_xxxyyz_yyy_0[i] = ta_xxxyy_yyy_0[i] * pa_z[i] - ta_xxxyy_yyy_1[i] * pc_z[i];

        ta_xxxyyz_yyz_0[i] =
            2.0 * ta_xyyz_yyz_0[i] * fe_0 - 2.0 * ta_xyyz_yyz_1[i] * fe_0 + ta_xxyyz_yyz_0[i] * pa_x[i] - ta_xxyyz_yyz_1[i] * pc_x[i];

        ta_xxxyyz_yzz_0[i] =
            2.0 * ta_xyyz_yzz_0[i] * fe_0 - 2.0 * ta_xyyz_yzz_1[i] * fe_0 + ta_xxyyz_yzz_0[i] * pa_x[i] - ta_xxyyz_yzz_1[i] * pc_x[i];

        ta_xxxyyz_zzz_0[i] =
            2.0 * ta_xyyz_zzz_0[i] * fe_0 - 2.0 * ta_xyyz_zzz_1[i] * fe_0 + ta_xxyyz_zzz_0[i] * pa_x[i] - ta_xxyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : IF

    auto ta_xxxyzz_xxx_0 = pbuffer.data(idx_npot_0_if + 80);

    auto ta_xxxyzz_xxy_0 = pbuffer.data(idx_npot_0_if + 81);

    auto ta_xxxyzz_xxz_0 = pbuffer.data(idx_npot_0_if + 82);

    auto ta_xxxyzz_xyy_0 = pbuffer.data(idx_npot_0_if + 83);

    auto ta_xxxyzz_xyz_0 = pbuffer.data(idx_npot_0_if + 84);

    auto ta_xxxyzz_xzz_0 = pbuffer.data(idx_npot_0_if + 85);

    auto ta_xxxyzz_yyy_0 = pbuffer.data(idx_npot_0_if + 86);

    auto ta_xxxyzz_yyz_0 = pbuffer.data(idx_npot_0_if + 87);

    auto ta_xxxyzz_yzz_0 = pbuffer.data(idx_npot_0_if + 88);

    auto ta_xxxyzz_zzz_0 = pbuffer.data(idx_npot_0_if + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxxyzz_xxx_0, \
                             ta_xxxyzz_xxy_0, \
                             ta_xxxyzz_xxz_0, \
                             ta_xxxyzz_xyy_0, \
                             ta_xxxyzz_xyz_0, \
                             ta_xxxyzz_xzz_0, \
                             ta_xxxyzz_yyy_0, \
                             ta_xxxyzz_yyz_0, \
                             ta_xxxyzz_yzz_0, \
                             ta_xxxyzz_zzz_0, \
                             ta_xxxzz_xx_0,   \
                             ta_xxxzz_xx_1,   \
                             ta_xxxzz_xxx_0,  \
                             ta_xxxzz_xxx_1,  \
                             ta_xxxzz_xxy_0,  \
                             ta_xxxzz_xxy_1,  \
                             ta_xxxzz_xxz_0,  \
                             ta_xxxzz_xxz_1,  \
                             ta_xxxzz_xy_0,   \
                             ta_xxxzz_xy_1,   \
                             ta_xxxzz_xyy_0,  \
                             ta_xxxzz_xyy_1,  \
                             ta_xxxzz_xyz_0,  \
                             ta_xxxzz_xyz_1,  \
                             ta_xxxzz_xz_0,   \
                             ta_xxxzz_xz_1,   \
                             ta_xxxzz_xzz_0,  \
                             ta_xxxzz_xzz_1,  \
                             ta_xxxzz_zzz_0,  \
                             ta_xxxzz_zzz_1,  \
                             ta_xxyzz_yyy_0,  \
                             ta_xxyzz_yyy_1,  \
                             ta_xxyzz_yyz_0,  \
                             ta_xxyzz_yyz_1,  \
                             ta_xxyzz_yzz_0,  \
                             ta_xxyzz_yzz_1,  \
                             ta_xyzz_yyy_0,   \
                             ta_xyzz_yyy_1,   \
                             ta_xyzz_yyz_0,   \
                             ta_xyzz_yyz_1,   \
                             ta_xyzz_yzz_0,   \
                             ta_xyzz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_xxx_0[i] = ta_xxxzz_xxx_0[i] * pa_y[i] - ta_xxxzz_xxx_1[i] * pc_y[i];

        ta_xxxyzz_xxy_0[i] = ta_xxxzz_xx_0[i] * fe_0 - ta_xxxzz_xx_1[i] * fe_0 + ta_xxxzz_xxy_0[i] * pa_y[i] - ta_xxxzz_xxy_1[i] * pc_y[i];

        ta_xxxyzz_xxz_0[i] = ta_xxxzz_xxz_0[i] * pa_y[i] - ta_xxxzz_xxz_1[i] * pc_y[i];

        ta_xxxyzz_xyy_0[i] =
            2.0 * ta_xxxzz_xy_0[i] * fe_0 - 2.0 * ta_xxxzz_xy_1[i] * fe_0 + ta_xxxzz_xyy_0[i] * pa_y[i] - ta_xxxzz_xyy_1[i] * pc_y[i];

        ta_xxxyzz_xyz_0[i] = ta_xxxzz_xz_0[i] * fe_0 - ta_xxxzz_xz_1[i] * fe_0 + ta_xxxzz_xyz_0[i] * pa_y[i] - ta_xxxzz_xyz_1[i] * pc_y[i];

        ta_xxxyzz_xzz_0[i] = ta_xxxzz_xzz_0[i] * pa_y[i] - ta_xxxzz_xzz_1[i] * pc_y[i];

        ta_xxxyzz_yyy_0[i] =
            2.0 * ta_xyzz_yyy_0[i] * fe_0 - 2.0 * ta_xyzz_yyy_1[i] * fe_0 + ta_xxyzz_yyy_0[i] * pa_x[i] - ta_xxyzz_yyy_1[i] * pc_x[i];

        ta_xxxyzz_yyz_0[i] =
            2.0 * ta_xyzz_yyz_0[i] * fe_0 - 2.0 * ta_xyzz_yyz_1[i] * fe_0 + ta_xxyzz_yyz_0[i] * pa_x[i] - ta_xxyzz_yyz_1[i] * pc_x[i];

        ta_xxxyzz_yzz_0[i] =
            2.0 * ta_xyzz_yzz_0[i] * fe_0 - 2.0 * ta_xyzz_yzz_1[i] * fe_0 + ta_xxyzz_yzz_0[i] * pa_x[i] - ta_xxyzz_yzz_1[i] * pc_x[i];

        ta_xxxyzz_zzz_0[i] = ta_xxxzz_zzz_0[i] * pa_y[i] - ta_xxxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : IF

    auto ta_xxxzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 90);

    auto ta_xxxzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 91);

    auto ta_xxxzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 92);

    auto ta_xxxzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 93);

    auto ta_xxxzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 94);

    auto ta_xxxzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 95);

    auto ta_xxxzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 96);

    auto ta_xxxzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 97);

    auto ta_xxxzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 98);

    auto ta_xxxzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 99);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxxz_xxx_0,   \
                             ta_xxxz_xxx_1,   \
                             ta_xxxz_xxy_0,   \
                             ta_xxxz_xxy_1,   \
                             ta_xxxz_xyy_0,   \
                             ta_xxxz_xyy_1,   \
                             ta_xxxzz_xxx_0,  \
                             ta_xxxzz_xxx_1,  \
                             ta_xxxzz_xxy_0,  \
                             ta_xxxzz_xxy_1,  \
                             ta_xxxzz_xyy_0,  \
                             ta_xxxzz_xyy_1,  \
                             ta_xxxzzz_xxx_0, \
                             ta_xxxzzz_xxy_0, \
                             ta_xxxzzz_xxz_0, \
                             ta_xxxzzz_xyy_0, \
                             ta_xxxzzz_xyz_0, \
                             ta_xxxzzz_xzz_0, \
                             ta_xxxzzz_yyy_0, \
                             ta_xxxzzz_yyz_0, \
                             ta_xxxzzz_yzz_0, \
                             ta_xxxzzz_zzz_0, \
                             ta_xxzzz_xxz_0,  \
                             ta_xxzzz_xxz_1,  \
                             ta_xxzzz_xyz_0,  \
                             ta_xxzzz_xyz_1,  \
                             ta_xxzzz_xz_0,   \
                             ta_xxzzz_xz_1,   \
                             ta_xxzzz_xzz_0,  \
                             ta_xxzzz_xzz_1,  \
                             ta_xxzzz_yyy_0,  \
                             ta_xxzzz_yyy_1,  \
                             ta_xxzzz_yyz_0,  \
                             ta_xxzzz_yyz_1,  \
                             ta_xxzzz_yz_0,   \
                             ta_xxzzz_yz_1,   \
                             ta_xxzzz_yzz_0,  \
                             ta_xxzzz_yzz_1,  \
                             ta_xxzzz_zz_0,   \
                             ta_xxzzz_zz_1,   \
                             ta_xxzzz_zzz_0,  \
                             ta_xxzzz_zzz_1,  \
                             ta_xzzz_xxz_0,   \
                             ta_xzzz_xxz_1,   \
                             ta_xzzz_xyz_0,   \
                             ta_xzzz_xyz_1,   \
                             ta_xzzz_xzz_0,   \
                             ta_xzzz_xzz_1,   \
                             ta_xzzz_yyy_0,   \
                             ta_xzzz_yyy_1,   \
                             ta_xzzz_yyz_0,   \
                             ta_xzzz_yyz_1,   \
                             ta_xzzz_yzz_0,   \
                             ta_xzzz_yzz_1,   \
                             ta_xzzz_zzz_0,   \
                             ta_xzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_xxx_0[i] =
            2.0 * ta_xxxz_xxx_0[i] * fe_0 - 2.0 * ta_xxxz_xxx_1[i] * fe_0 + ta_xxxzz_xxx_0[i] * pa_z[i] - ta_xxxzz_xxx_1[i] * pc_z[i];

        ta_xxxzzz_xxy_0[i] =
            2.0 * ta_xxxz_xxy_0[i] * fe_0 - 2.0 * ta_xxxz_xxy_1[i] * fe_0 + ta_xxxzz_xxy_0[i] * pa_z[i] - ta_xxxzz_xxy_1[i] * pc_z[i];

        ta_xxxzzz_xxz_0[i] = 2.0 * ta_xzzz_xxz_0[i] * fe_0 - 2.0 * ta_xzzz_xxz_1[i] * fe_0 + 2.0 * ta_xxzzz_xz_0[i] * fe_0 -
                             2.0 * ta_xxzzz_xz_1[i] * fe_0 + ta_xxzzz_xxz_0[i] * pa_x[i] - ta_xxzzz_xxz_1[i] * pc_x[i];

        ta_xxxzzz_xyy_0[i] =
            2.0 * ta_xxxz_xyy_0[i] * fe_0 - 2.0 * ta_xxxz_xyy_1[i] * fe_0 + ta_xxxzz_xyy_0[i] * pa_z[i] - ta_xxxzz_xyy_1[i] * pc_z[i];

        ta_xxxzzz_xyz_0[i] = 2.0 * ta_xzzz_xyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyz_1[i] * fe_0 + ta_xxzzz_yz_0[i] * fe_0 - ta_xxzzz_yz_1[i] * fe_0 +
                             ta_xxzzz_xyz_0[i] * pa_x[i] - ta_xxzzz_xyz_1[i] * pc_x[i];

        ta_xxxzzz_xzz_0[i] = 2.0 * ta_xzzz_xzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzz_1[i] * fe_0 + ta_xxzzz_zz_0[i] * fe_0 - ta_xxzzz_zz_1[i] * fe_0 +
                             ta_xxzzz_xzz_0[i] * pa_x[i] - ta_xxzzz_xzz_1[i] * pc_x[i];

        ta_xxxzzz_yyy_0[i] =
            2.0 * ta_xzzz_yyy_0[i] * fe_0 - 2.0 * ta_xzzz_yyy_1[i] * fe_0 + ta_xxzzz_yyy_0[i] * pa_x[i] - ta_xxzzz_yyy_1[i] * pc_x[i];

        ta_xxxzzz_yyz_0[i] =
            2.0 * ta_xzzz_yyz_0[i] * fe_0 - 2.0 * ta_xzzz_yyz_1[i] * fe_0 + ta_xxzzz_yyz_0[i] * pa_x[i] - ta_xxzzz_yyz_1[i] * pc_x[i];

        ta_xxxzzz_yzz_0[i] =
            2.0 * ta_xzzz_yzz_0[i] * fe_0 - 2.0 * ta_xzzz_yzz_1[i] * fe_0 + ta_xxzzz_yzz_0[i] * pa_x[i] - ta_xxzzz_yzz_1[i] * pc_x[i];

        ta_xxxzzz_zzz_0[i] =
            2.0 * ta_xzzz_zzz_0[i] * fe_0 - 2.0 * ta_xzzz_zzz_1[i] * fe_0 + ta_xxzzz_zzz_0[i] * pa_x[i] - ta_xxzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : IF

    auto ta_xxyyyy_xxx_0 = pbuffer.data(idx_npot_0_if + 100);

    auto ta_xxyyyy_xxy_0 = pbuffer.data(idx_npot_0_if + 101);

    auto ta_xxyyyy_xxz_0 = pbuffer.data(idx_npot_0_if + 102);

    auto ta_xxyyyy_xyy_0 = pbuffer.data(idx_npot_0_if + 103);

    auto ta_xxyyyy_xyz_0 = pbuffer.data(idx_npot_0_if + 104);

    auto ta_xxyyyy_xzz_0 = pbuffer.data(idx_npot_0_if + 105);

    auto ta_xxyyyy_yyy_0 = pbuffer.data(idx_npot_0_if + 106);

    auto ta_xxyyyy_yyz_0 = pbuffer.data(idx_npot_0_if + 107);

    auto ta_xxyyyy_yzz_0 = pbuffer.data(idx_npot_0_if + 108);

    auto ta_xxyyyy_zzz_0 = pbuffer.data(idx_npot_0_if + 109);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxyy_xxx_0,   \
                             ta_xxyy_xxx_1,   \
                             ta_xxyy_xxz_0,   \
                             ta_xxyy_xxz_1,   \
                             ta_xxyy_xzz_0,   \
                             ta_xxyy_xzz_1,   \
                             ta_xxyyy_xxx_0,  \
                             ta_xxyyy_xxx_1,  \
                             ta_xxyyy_xxz_0,  \
                             ta_xxyyy_xxz_1,  \
                             ta_xxyyy_xzz_0,  \
                             ta_xxyyy_xzz_1,  \
                             ta_xxyyyy_xxx_0, \
                             ta_xxyyyy_xxy_0, \
                             ta_xxyyyy_xxz_0, \
                             ta_xxyyyy_xyy_0, \
                             ta_xxyyyy_xyz_0, \
                             ta_xxyyyy_xzz_0, \
                             ta_xxyyyy_yyy_0, \
                             ta_xxyyyy_yyz_0, \
                             ta_xxyyyy_yzz_0, \
                             ta_xxyyyy_zzz_0, \
                             ta_xyyyy_xxy_0,  \
                             ta_xyyyy_xxy_1,  \
                             ta_xyyyy_xy_0,   \
                             ta_xyyyy_xy_1,   \
                             ta_xyyyy_xyy_0,  \
                             ta_xyyyy_xyy_1,  \
                             ta_xyyyy_xyz_0,  \
                             ta_xyyyy_xyz_1,  \
                             ta_xyyyy_yy_0,   \
                             ta_xyyyy_yy_1,   \
                             ta_xyyyy_yyy_0,  \
                             ta_xyyyy_yyy_1,  \
                             ta_xyyyy_yyz_0,  \
                             ta_xyyyy_yyz_1,  \
                             ta_xyyyy_yz_0,   \
                             ta_xyyyy_yz_1,   \
                             ta_xyyyy_yzz_0,  \
                             ta_xyyyy_yzz_1,  \
                             ta_xyyyy_zzz_0,  \
                             ta_xyyyy_zzz_1,  \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_xyz_0,   \
                             ta_yyyy_xyz_1,   \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyy_yyz_0,   \
                             ta_yyyy_yyz_1,   \
                             ta_yyyy_yzz_0,   \
                             ta_yyyy_yzz_1,   \
                             ta_yyyy_zzz_0,   \
                             ta_yyyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_xxx_0[i] =
            3.0 * ta_xxyy_xxx_0[i] * fe_0 - 3.0 * ta_xxyy_xxx_1[i] * fe_0 + ta_xxyyy_xxx_0[i] * pa_y[i] - ta_xxyyy_xxx_1[i] * pc_y[i];

        ta_xxyyyy_xxy_0[i] = ta_yyyy_xxy_0[i] * fe_0 - ta_yyyy_xxy_1[i] * fe_0 + 2.0 * ta_xyyyy_xy_0[i] * fe_0 - 2.0 * ta_xyyyy_xy_1[i] * fe_0 +
                             ta_xyyyy_xxy_0[i] * pa_x[i] - ta_xyyyy_xxy_1[i] * pc_x[i];

        ta_xxyyyy_xxz_0[i] =
            3.0 * ta_xxyy_xxz_0[i] * fe_0 - 3.0 * ta_xxyy_xxz_1[i] * fe_0 + ta_xxyyy_xxz_0[i] * pa_y[i] - ta_xxyyy_xxz_1[i] * pc_y[i];

        ta_xxyyyy_xyy_0[i] = ta_yyyy_xyy_0[i] * fe_0 - ta_yyyy_xyy_1[i] * fe_0 + ta_xyyyy_yy_0[i] * fe_0 - ta_xyyyy_yy_1[i] * fe_0 +
                             ta_xyyyy_xyy_0[i] * pa_x[i] - ta_xyyyy_xyy_1[i] * pc_x[i];

        ta_xxyyyy_xyz_0[i] = ta_yyyy_xyz_0[i] * fe_0 - ta_yyyy_xyz_1[i] * fe_0 + ta_xyyyy_yz_0[i] * fe_0 - ta_xyyyy_yz_1[i] * fe_0 +
                             ta_xyyyy_xyz_0[i] * pa_x[i] - ta_xyyyy_xyz_1[i] * pc_x[i];

        ta_xxyyyy_xzz_0[i] =
            3.0 * ta_xxyy_xzz_0[i] * fe_0 - 3.0 * ta_xxyy_xzz_1[i] * fe_0 + ta_xxyyy_xzz_0[i] * pa_y[i] - ta_xxyyy_xzz_1[i] * pc_y[i];

        ta_xxyyyy_yyy_0[i] = ta_yyyy_yyy_0[i] * fe_0 - ta_yyyy_yyy_1[i] * fe_0 + ta_xyyyy_yyy_0[i] * pa_x[i] - ta_xyyyy_yyy_1[i] * pc_x[i];

        ta_xxyyyy_yyz_0[i] = ta_yyyy_yyz_0[i] * fe_0 - ta_yyyy_yyz_1[i] * fe_0 + ta_xyyyy_yyz_0[i] * pa_x[i] - ta_xyyyy_yyz_1[i] * pc_x[i];

        ta_xxyyyy_yzz_0[i] = ta_yyyy_yzz_0[i] * fe_0 - ta_yyyy_yzz_1[i] * fe_0 + ta_xyyyy_yzz_0[i] * pa_x[i] - ta_xyyyy_yzz_1[i] * pc_x[i];

        ta_xxyyyy_zzz_0[i] = ta_yyyy_zzz_0[i] * fe_0 - ta_yyyy_zzz_1[i] * fe_0 + ta_xyyyy_zzz_0[i] * pa_x[i] - ta_xyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 110-120 components of targeted buffer : IF

    auto ta_xxyyyz_xxx_0 = pbuffer.data(idx_npot_0_if + 110);

    auto ta_xxyyyz_xxy_0 = pbuffer.data(idx_npot_0_if + 111);

    auto ta_xxyyyz_xxz_0 = pbuffer.data(idx_npot_0_if + 112);

    auto ta_xxyyyz_xyy_0 = pbuffer.data(idx_npot_0_if + 113);

    auto ta_xxyyyz_xyz_0 = pbuffer.data(idx_npot_0_if + 114);

    auto ta_xxyyyz_xzz_0 = pbuffer.data(idx_npot_0_if + 115);

    auto ta_xxyyyz_yyy_0 = pbuffer.data(idx_npot_0_if + 116);

    auto ta_xxyyyz_yyz_0 = pbuffer.data(idx_npot_0_if + 117);

    auto ta_xxyyyz_yzz_0 = pbuffer.data(idx_npot_0_if + 118);

    auto ta_xxyyyz_zzz_0 = pbuffer.data(idx_npot_0_if + 119);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxyyy_xxx_0,  \
                             ta_xxyyy_xxx_1,  \
                             ta_xxyyy_xxy_0,  \
                             ta_xxyyy_xxy_1,  \
                             ta_xxyyy_xy_0,   \
                             ta_xxyyy_xy_1,   \
                             ta_xxyyy_xyy_0,  \
                             ta_xxyyy_xyy_1,  \
                             ta_xxyyy_xyz_0,  \
                             ta_xxyyy_xyz_1,  \
                             ta_xxyyy_yyy_0,  \
                             ta_xxyyy_yyy_1,  \
                             ta_xxyyyz_xxx_0, \
                             ta_xxyyyz_xxy_0, \
                             ta_xxyyyz_xxz_0, \
                             ta_xxyyyz_xyy_0, \
                             ta_xxyyyz_xyz_0, \
                             ta_xxyyyz_xzz_0, \
                             ta_xxyyyz_yyy_0, \
                             ta_xxyyyz_yyz_0, \
                             ta_xxyyyz_yzz_0, \
                             ta_xxyyyz_zzz_0, \
                             ta_xxyyz_xxz_0,  \
                             ta_xxyyz_xxz_1,  \
                             ta_xxyyz_xzz_0,  \
                             ta_xxyyz_xzz_1,  \
                             ta_xxyz_xxz_0,   \
                             ta_xxyz_xxz_1,   \
                             ta_xxyz_xzz_0,   \
                             ta_xxyz_xzz_1,   \
                             ta_xyyyz_yyz_0,  \
                             ta_xyyyz_yyz_1,  \
                             ta_xyyyz_yzz_0,  \
                             ta_xyyyz_yzz_1,  \
                             ta_xyyyz_zzz_0,  \
                             ta_xyyyz_zzz_1,  \
                             ta_yyyz_yyz_0,   \
                             ta_yyyz_yyz_1,   \
                             ta_yyyz_yzz_0,   \
                             ta_yyyz_yzz_1,   \
                             ta_yyyz_zzz_0,   \
                             ta_yyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_xxx_0[i] = ta_xxyyy_xxx_0[i] * pa_z[i] - ta_xxyyy_xxx_1[i] * pc_z[i];

        ta_xxyyyz_xxy_0[i] = ta_xxyyy_xxy_0[i] * pa_z[i] - ta_xxyyy_xxy_1[i] * pc_z[i];

        ta_xxyyyz_xxz_0[i] =
            2.0 * ta_xxyz_xxz_0[i] * fe_0 - 2.0 * ta_xxyz_xxz_1[i] * fe_0 + ta_xxyyz_xxz_0[i] * pa_y[i] - ta_xxyyz_xxz_1[i] * pc_y[i];

        ta_xxyyyz_xyy_0[i] = ta_xxyyy_xyy_0[i] * pa_z[i] - ta_xxyyy_xyy_1[i] * pc_z[i];

        ta_xxyyyz_xyz_0[i] = ta_xxyyy_xy_0[i] * fe_0 - ta_xxyyy_xy_1[i] * fe_0 + ta_xxyyy_xyz_0[i] * pa_z[i] - ta_xxyyy_xyz_1[i] * pc_z[i];

        ta_xxyyyz_xzz_0[i] =
            2.0 * ta_xxyz_xzz_0[i] * fe_0 - 2.0 * ta_xxyz_xzz_1[i] * fe_0 + ta_xxyyz_xzz_0[i] * pa_y[i] - ta_xxyyz_xzz_1[i] * pc_y[i];

        ta_xxyyyz_yyy_0[i] = ta_xxyyy_yyy_0[i] * pa_z[i] - ta_xxyyy_yyy_1[i] * pc_z[i];

        ta_xxyyyz_yyz_0[i] = ta_yyyz_yyz_0[i] * fe_0 - ta_yyyz_yyz_1[i] * fe_0 + ta_xyyyz_yyz_0[i] * pa_x[i] - ta_xyyyz_yyz_1[i] * pc_x[i];

        ta_xxyyyz_yzz_0[i] = ta_yyyz_yzz_0[i] * fe_0 - ta_yyyz_yzz_1[i] * fe_0 + ta_xyyyz_yzz_0[i] * pa_x[i] - ta_xyyyz_yzz_1[i] * pc_x[i];

        ta_xxyyyz_zzz_0[i] = ta_yyyz_zzz_0[i] * fe_0 - ta_yyyz_zzz_1[i] * fe_0 + ta_xyyyz_zzz_0[i] * pa_x[i] - ta_xyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 120-130 components of targeted buffer : IF

    auto ta_xxyyzz_xxx_0 = pbuffer.data(idx_npot_0_if + 120);

    auto ta_xxyyzz_xxy_0 = pbuffer.data(idx_npot_0_if + 121);

    auto ta_xxyyzz_xxz_0 = pbuffer.data(idx_npot_0_if + 122);

    auto ta_xxyyzz_xyy_0 = pbuffer.data(idx_npot_0_if + 123);

    auto ta_xxyyzz_xyz_0 = pbuffer.data(idx_npot_0_if + 124);

    auto ta_xxyyzz_xzz_0 = pbuffer.data(idx_npot_0_if + 125);

    auto ta_xxyyzz_yyy_0 = pbuffer.data(idx_npot_0_if + 126);

    auto ta_xxyyzz_yyz_0 = pbuffer.data(idx_npot_0_if + 127);

    auto ta_xxyyzz_yzz_0 = pbuffer.data(idx_npot_0_if + 128);

    auto ta_xxyyzz_zzz_0 = pbuffer.data(idx_npot_0_if + 129);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxyy_xxy_0,   \
                             ta_xxyy_xxy_1,   \
                             ta_xxyy_xyy_0,   \
                             ta_xxyy_xyy_1,   \
                             ta_xxyyz_xxy_0,  \
                             ta_xxyyz_xxy_1,  \
                             ta_xxyyz_xyy_0,  \
                             ta_xxyyz_xyy_1,  \
                             ta_xxyyzz_xxx_0, \
                             ta_xxyyzz_xxy_0, \
                             ta_xxyyzz_xxz_0, \
                             ta_xxyyzz_xyy_0, \
                             ta_xxyyzz_xyz_0, \
                             ta_xxyyzz_xzz_0, \
                             ta_xxyyzz_yyy_0, \
                             ta_xxyyzz_yyz_0, \
                             ta_xxyyzz_yzz_0, \
                             ta_xxyyzz_zzz_0, \
                             ta_xxyzz_xxx_0,  \
                             ta_xxyzz_xxx_1,  \
                             ta_xxyzz_xxz_0,  \
                             ta_xxyzz_xxz_1,  \
                             ta_xxyzz_xzz_0,  \
                             ta_xxyzz_xzz_1,  \
                             ta_xxzz_xxx_0,   \
                             ta_xxzz_xxx_1,   \
                             ta_xxzz_xxz_0,   \
                             ta_xxzz_xxz_1,   \
                             ta_xxzz_xzz_0,   \
                             ta_xxzz_xzz_1,   \
                             ta_xyyzz_xyz_0,  \
                             ta_xyyzz_xyz_1,  \
                             ta_xyyzz_yyy_0,  \
                             ta_xyyzz_yyy_1,  \
                             ta_xyyzz_yyz_0,  \
                             ta_xyyzz_yyz_1,  \
                             ta_xyyzz_yz_0,   \
                             ta_xyyzz_yz_1,   \
                             ta_xyyzz_yzz_0,  \
                             ta_xyyzz_yzz_1,  \
                             ta_xyyzz_zzz_0,  \
                             ta_xyyzz_zzz_1,  \
                             ta_yyzz_xyz_0,   \
                             ta_yyzz_xyz_1,   \
                             ta_yyzz_yyy_0,   \
                             ta_yyzz_yyy_1,   \
                             ta_yyzz_yyz_0,   \
                             ta_yyzz_yyz_1,   \
                             ta_yyzz_yzz_0,   \
                             ta_yyzz_yzz_1,   \
                             ta_yyzz_zzz_0,   \
                             ta_yyzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_xxx_0[i] = ta_xxzz_xxx_0[i] * fe_0 - ta_xxzz_xxx_1[i] * fe_0 + ta_xxyzz_xxx_0[i] * pa_y[i] - ta_xxyzz_xxx_1[i] * pc_y[i];

        ta_xxyyzz_xxy_0[i] = ta_xxyy_xxy_0[i] * fe_0 - ta_xxyy_xxy_1[i] * fe_0 + ta_xxyyz_xxy_0[i] * pa_z[i] - ta_xxyyz_xxy_1[i] * pc_z[i];

        ta_xxyyzz_xxz_0[i] = ta_xxzz_xxz_0[i] * fe_0 - ta_xxzz_xxz_1[i] * fe_0 + ta_xxyzz_xxz_0[i] * pa_y[i] - ta_xxyzz_xxz_1[i] * pc_y[i];

        ta_xxyyzz_xyy_0[i] = ta_xxyy_xyy_0[i] * fe_0 - ta_xxyy_xyy_1[i] * fe_0 + ta_xxyyz_xyy_0[i] * pa_z[i] - ta_xxyyz_xyy_1[i] * pc_z[i];

        ta_xxyyzz_xyz_0[i] = ta_yyzz_xyz_0[i] * fe_0 - ta_yyzz_xyz_1[i] * fe_0 + ta_xyyzz_yz_0[i] * fe_0 - ta_xyyzz_yz_1[i] * fe_0 +
                             ta_xyyzz_xyz_0[i] * pa_x[i] - ta_xyyzz_xyz_1[i] * pc_x[i];

        ta_xxyyzz_xzz_0[i] = ta_xxzz_xzz_0[i] * fe_0 - ta_xxzz_xzz_1[i] * fe_0 + ta_xxyzz_xzz_0[i] * pa_y[i] - ta_xxyzz_xzz_1[i] * pc_y[i];

        ta_xxyyzz_yyy_0[i] = ta_yyzz_yyy_0[i] * fe_0 - ta_yyzz_yyy_1[i] * fe_0 + ta_xyyzz_yyy_0[i] * pa_x[i] - ta_xyyzz_yyy_1[i] * pc_x[i];

        ta_xxyyzz_yyz_0[i] = ta_yyzz_yyz_0[i] * fe_0 - ta_yyzz_yyz_1[i] * fe_0 + ta_xyyzz_yyz_0[i] * pa_x[i] - ta_xyyzz_yyz_1[i] * pc_x[i];

        ta_xxyyzz_yzz_0[i] = ta_yyzz_yzz_0[i] * fe_0 - ta_yyzz_yzz_1[i] * fe_0 + ta_xyyzz_yzz_0[i] * pa_x[i] - ta_xyyzz_yzz_1[i] * pc_x[i];

        ta_xxyyzz_zzz_0[i] = ta_yyzz_zzz_0[i] * fe_0 - ta_yyzz_zzz_1[i] * fe_0 + ta_xyyzz_zzz_0[i] * pa_x[i] - ta_xyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : IF

    auto ta_xxyzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 130);

    auto ta_xxyzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 131);

    auto ta_xxyzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 132);

    auto ta_xxyzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 133);

    auto ta_xxyzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 134);

    auto ta_xxyzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 135);

    auto ta_xxyzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 136);

    auto ta_xxyzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 137);

    auto ta_xxyzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 138);

    auto ta_xxyzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 139);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxyzzz_xxx_0, \
                             ta_xxyzzz_xxy_0, \
                             ta_xxyzzz_xxz_0, \
                             ta_xxyzzz_xyy_0, \
                             ta_xxyzzz_xyz_0, \
                             ta_xxyzzz_xzz_0, \
                             ta_xxyzzz_yyy_0, \
                             ta_xxyzzz_yyz_0, \
                             ta_xxyzzz_yzz_0, \
                             ta_xxyzzz_zzz_0, \
                             ta_xxzzz_xx_0,   \
                             ta_xxzzz_xx_1,   \
                             ta_xxzzz_xxx_0,  \
                             ta_xxzzz_xxx_1,  \
                             ta_xxzzz_xxy_0,  \
                             ta_xxzzz_xxy_1,  \
                             ta_xxzzz_xxz_0,  \
                             ta_xxzzz_xxz_1,  \
                             ta_xxzzz_xy_0,   \
                             ta_xxzzz_xy_1,   \
                             ta_xxzzz_xyy_0,  \
                             ta_xxzzz_xyy_1,  \
                             ta_xxzzz_xyz_0,  \
                             ta_xxzzz_xyz_1,  \
                             ta_xxzzz_xz_0,   \
                             ta_xxzzz_xz_1,   \
                             ta_xxzzz_xzz_0,  \
                             ta_xxzzz_xzz_1,  \
                             ta_xxzzz_zzz_0,  \
                             ta_xxzzz_zzz_1,  \
                             ta_xyzzz_yyy_0,  \
                             ta_xyzzz_yyy_1,  \
                             ta_xyzzz_yyz_0,  \
                             ta_xyzzz_yyz_1,  \
                             ta_xyzzz_yzz_0,  \
                             ta_xyzzz_yzz_1,  \
                             ta_yzzz_yyy_0,   \
                             ta_yzzz_yyy_1,   \
                             ta_yzzz_yyz_0,   \
                             ta_yzzz_yyz_1,   \
                             ta_yzzz_yzz_0,   \
                             ta_yzzz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_xxx_0[i] = ta_xxzzz_xxx_0[i] * pa_y[i] - ta_xxzzz_xxx_1[i] * pc_y[i];

        ta_xxyzzz_xxy_0[i] = ta_xxzzz_xx_0[i] * fe_0 - ta_xxzzz_xx_1[i] * fe_0 + ta_xxzzz_xxy_0[i] * pa_y[i] - ta_xxzzz_xxy_1[i] * pc_y[i];

        ta_xxyzzz_xxz_0[i] = ta_xxzzz_xxz_0[i] * pa_y[i] - ta_xxzzz_xxz_1[i] * pc_y[i];

        ta_xxyzzz_xyy_0[i] =
            2.0 * ta_xxzzz_xy_0[i] * fe_0 - 2.0 * ta_xxzzz_xy_1[i] * fe_0 + ta_xxzzz_xyy_0[i] * pa_y[i] - ta_xxzzz_xyy_1[i] * pc_y[i];

        ta_xxyzzz_xyz_0[i] = ta_xxzzz_xz_0[i] * fe_0 - ta_xxzzz_xz_1[i] * fe_0 + ta_xxzzz_xyz_0[i] * pa_y[i] - ta_xxzzz_xyz_1[i] * pc_y[i];

        ta_xxyzzz_xzz_0[i] = ta_xxzzz_xzz_0[i] * pa_y[i] - ta_xxzzz_xzz_1[i] * pc_y[i];

        ta_xxyzzz_yyy_0[i] = ta_yzzz_yyy_0[i] * fe_0 - ta_yzzz_yyy_1[i] * fe_0 + ta_xyzzz_yyy_0[i] * pa_x[i] - ta_xyzzz_yyy_1[i] * pc_x[i];

        ta_xxyzzz_yyz_0[i] = ta_yzzz_yyz_0[i] * fe_0 - ta_yzzz_yyz_1[i] * fe_0 + ta_xyzzz_yyz_0[i] * pa_x[i] - ta_xyzzz_yyz_1[i] * pc_x[i];

        ta_xxyzzz_yzz_0[i] = ta_yzzz_yzz_0[i] * fe_0 - ta_yzzz_yzz_1[i] * fe_0 + ta_xyzzz_yzz_0[i] * pa_x[i] - ta_xyzzz_yzz_1[i] * pc_x[i];

        ta_xxyzzz_zzz_0[i] = ta_xxzzz_zzz_0[i] * pa_y[i] - ta_xxzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : IF

    auto ta_xxzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 140);

    auto ta_xxzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 141);

    auto ta_xxzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 142);

    auto ta_xxzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 143);

    auto ta_xxzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 144);

    auto ta_xxzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 145);

    auto ta_xxzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 146);

    auto ta_xxzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 147);

    auto ta_xxzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 148);

    auto ta_xxzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxzz_xxx_0,   \
                             ta_xxzz_xxx_1,   \
                             ta_xxzz_xxy_0,   \
                             ta_xxzz_xxy_1,   \
                             ta_xxzz_xyy_0,   \
                             ta_xxzz_xyy_1,   \
                             ta_xxzzz_xxx_0,  \
                             ta_xxzzz_xxx_1,  \
                             ta_xxzzz_xxy_0,  \
                             ta_xxzzz_xxy_1,  \
                             ta_xxzzz_xyy_0,  \
                             ta_xxzzz_xyy_1,  \
                             ta_xxzzzz_xxx_0, \
                             ta_xxzzzz_xxy_0, \
                             ta_xxzzzz_xxz_0, \
                             ta_xxzzzz_xyy_0, \
                             ta_xxzzzz_xyz_0, \
                             ta_xxzzzz_xzz_0, \
                             ta_xxzzzz_yyy_0, \
                             ta_xxzzzz_yyz_0, \
                             ta_xxzzzz_yzz_0, \
                             ta_xxzzzz_zzz_0, \
                             ta_xzzzz_xxz_0,  \
                             ta_xzzzz_xxz_1,  \
                             ta_xzzzz_xyz_0,  \
                             ta_xzzzz_xyz_1,  \
                             ta_xzzzz_xz_0,   \
                             ta_xzzzz_xz_1,   \
                             ta_xzzzz_xzz_0,  \
                             ta_xzzzz_xzz_1,  \
                             ta_xzzzz_yyy_0,  \
                             ta_xzzzz_yyy_1,  \
                             ta_xzzzz_yyz_0,  \
                             ta_xzzzz_yyz_1,  \
                             ta_xzzzz_yz_0,   \
                             ta_xzzzz_yz_1,   \
                             ta_xzzzz_yzz_0,  \
                             ta_xzzzz_yzz_1,  \
                             ta_xzzzz_zz_0,   \
                             ta_xzzzz_zz_1,   \
                             ta_xzzzz_zzz_0,  \
                             ta_xzzzz_zzz_1,  \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_yyy_0,   \
                             ta_zzzz_yyy_1,   \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_xxx_0[i] =
            3.0 * ta_xxzz_xxx_0[i] * fe_0 - 3.0 * ta_xxzz_xxx_1[i] * fe_0 + ta_xxzzz_xxx_0[i] * pa_z[i] - ta_xxzzz_xxx_1[i] * pc_z[i];

        ta_xxzzzz_xxy_0[i] =
            3.0 * ta_xxzz_xxy_0[i] * fe_0 - 3.0 * ta_xxzz_xxy_1[i] * fe_0 + ta_xxzzz_xxy_0[i] * pa_z[i] - ta_xxzzz_xxy_1[i] * pc_z[i];

        ta_xxzzzz_xxz_0[i] = ta_zzzz_xxz_0[i] * fe_0 - ta_zzzz_xxz_1[i] * fe_0 + 2.0 * ta_xzzzz_xz_0[i] * fe_0 - 2.0 * ta_xzzzz_xz_1[i] * fe_0 +
                             ta_xzzzz_xxz_0[i] * pa_x[i] - ta_xzzzz_xxz_1[i] * pc_x[i];

        ta_xxzzzz_xyy_0[i] =
            3.0 * ta_xxzz_xyy_0[i] * fe_0 - 3.0 * ta_xxzz_xyy_1[i] * fe_0 + ta_xxzzz_xyy_0[i] * pa_z[i] - ta_xxzzz_xyy_1[i] * pc_z[i];

        ta_xxzzzz_xyz_0[i] = ta_zzzz_xyz_0[i] * fe_0 - ta_zzzz_xyz_1[i] * fe_0 + ta_xzzzz_yz_0[i] * fe_0 - ta_xzzzz_yz_1[i] * fe_0 +
                             ta_xzzzz_xyz_0[i] * pa_x[i] - ta_xzzzz_xyz_1[i] * pc_x[i];

        ta_xxzzzz_xzz_0[i] = ta_zzzz_xzz_0[i] * fe_0 - ta_zzzz_xzz_1[i] * fe_0 + ta_xzzzz_zz_0[i] * fe_0 - ta_xzzzz_zz_1[i] * fe_0 +
                             ta_xzzzz_xzz_0[i] * pa_x[i] - ta_xzzzz_xzz_1[i] * pc_x[i];

        ta_xxzzzz_yyy_0[i] = ta_zzzz_yyy_0[i] * fe_0 - ta_zzzz_yyy_1[i] * fe_0 + ta_xzzzz_yyy_0[i] * pa_x[i] - ta_xzzzz_yyy_1[i] * pc_x[i];

        ta_xxzzzz_yyz_0[i] = ta_zzzz_yyz_0[i] * fe_0 - ta_zzzz_yyz_1[i] * fe_0 + ta_xzzzz_yyz_0[i] * pa_x[i] - ta_xzzzz_yyz_1[i] * pc_x[i];

        ta_xxzzzz_yzz_0[i] = ta_zzzz_yzz_0[i] * fe_0 - ta_zzzz_yzz_1[i] * fe_0 + ta_xzzzz_yzz_0[i] * pa_x[i] - ta_xzzzz_yzz_1[i] * pc_x[i];

        ta_xxzzzz_zzz_0[i] = ta_zzzz_zzz_0[i] * fe_0 - ta_zzzz_zzz_1[i] * fe_0 + ta_xzzzz_zzz_0[i] * pa_x[i] - ta_xzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : IF

    auto ta_xyyyyy_xxx_0 = pbuffer.data(idx_npot_0_if + 150);

    auto ta_xyyyyy_xxy_0 = pbuffer.data(idx_npot_0_if + 151);

    auto ta_xyyyyy_xxz_0 = pbuffer.data(idx_npot_0_if + 152);

    auto ta_xyyyyy_xyy_0 = pbuffer.data(idx_npot_0_if + 153);

    auto ta_xyyyyy_xyz_0 = pbuffer.data(idx_npot_0_if + 154);

    auto ta_xyyyyy_xzz_0 = pbuffer.data(idx_npot_0_if + 155);

    auto ta_xyyyyy_yyy_0 = pbuffer.data(idx_npot_0_if + 156);

    auto ta_xyyyyy_yyz_0 = pbuffer.data(idx_npot_0_if + 157);

    auto ta_xyyyyy_yzz_0 = pbuffer.data(idx_npot_0_if + 158);

    auto ta_xyyyyy_zzz_0 = pbuffer.data(idx_npot_0_if + 159);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xyyyyy_xxx_0, \
                             ta_xyyyyy_xxy_0, \
                             ta_xyyyyy_xxz_0, \
                             ta_xyyyyy_xyy_0, \
                             ta_xyyyyy_xyz_0, \
                             ta_xyyyyy_xzz_0, \
                             ta_xyyyyy_yyy_0, \
                             ta_xyyyyy_yyz_0, \
                             ta_xyyyyy_yzz_0, \
                             ta_xyyyyy_zzz_0, \
                             ta_yyyyy_xx_0,   \
                             ta_yyyyy_xx_1,   \
                             ta_yyyyy_xxx_0,  \
                             ta_yyyyy_xxx_1,  \
                             ta_yyyyy_xxy_0,  \
                             ta_yyyyy_xxy_1,  \
                             ta_yyyyy_xxz_0,  \
                             ta_yyyyy_xxz_1,  \
                             ta_yyyyy_xy_0,   \
                             ta_yyyyy_xy_1,   \
                             ta_yyyyy_xyy_0,  \
                             ta_yyyyy_xyy_1,  \
                             ta_yyyyy_xyz_0,  \
                             ta_yyyyy_xyz_1,  \
                             ta_yyyyy_xz_0,   \
                             ta_yyyyy_xz_1,   \
                             ta_yyyyy_xzz_0,  \
                             ta_yyyyy_xzz_1,  \
                             ta_yyyyy_yy_0,   \
                             ta_yyyyy_yy_1,   \
                             ta_yyyyy_yyy_0,  \
                             ta_yyyyy_yyy_1,  \
                             ta_yyyyy_yyz_0,  \
                             ta_yyyyy_yyz_1,  \
                             ta_yyyyy_yz_0,   \
                             ta_yyyyy_yz_1,   \
                             ta_yyyyy_yzz_0,  \
                             ta_yyyyy_yzz_1,  \
                             ta_yyyyy_zz_0,   \
                             ta_yyyyy_zz_1,   \
                             ta_yyyyy_zzz_0,  \
                             ta_yyyyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_xxx_0[i] =
            3.0 * ta_yyyyy_xx_0[i] * fe_0 - 3.0 * ta_yyyyy_xx_1[i] * fe_0 + ta_yyyyy_xxx_0[i] * pa_x[i] - ta_yyyyy_xxx_1[i] * pc_x[i];

        ta_xyyyyy_xxy_0[i] =
            2.0 * ta_yyyyy_xy_0[i] * fe_0 - 2.0 * ta_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xxy_0[i] * pa_x[i] - ta_yyyyy_xxy_1[i] * pc_x[i];

        ta_xyyyyy_xxz_0[i] =
            2.0 * ta_yyyyy_xz_0[i] * fe_0 - 2.0 * ta_yyyyy_xz_1[i] * fe_0 + ta_yyyyy_xxz_0[i] * pa_x[i] - ta_yyyyy_xxz_1[i] * pc_x[i];

        ta_xyyyyy_xyy_0[i] = ta_yyyyy_yy_0[i] * fe_0 - ta_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_xyy_0[i] * pa_x[i] - ta_yyyyy_xyy_1[i] * pc_x[i];

        ta_xyyyyy_xyz_0[i] = ta_yyyyy_yz_0[i] * fe_0 - ta_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_xyz_0[i] * pa_x[i] - ta_yyyyy_xyz_1[i] * pc_x[i];

        ta_xyyyyy_xzz_0[i] = ta_yyyyy_zz_0[i] * fe_0 - ta_yyyyy_zz_1[i] * fe_0 + ta_yyyyy_xzz_0[i] * pa_x[i] - ta_yyyyy_xzz_1[i] * pc_x[i];

        ta_xyyyyy_yyy_0[i] = ta_yyyyy_yyy_0[i] * pa_x[i] - ta_yyyyy_yyy_1[i] * pc_x[i];

        ta_xyyyyy_yyz_0[i] = ta_yyyyy_yyz_0[i] * pa_x[i] - ta_yyyyy_yyz_1[i] * pc_x[i];

        ta_xyyyyy_yzz_0[i] = ta_yyyyy_yzz_0[i] * pa_x[i] - ta_yyyyy_yzz_1[i] * pc_x[i];

        ta_xyyyyy_zzz_0[i] = ta_yyyyy_zzz_0[i] * pa_x[i] - ta_yyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : IF

    auto ta_xyyyyz_xxx_0 = pbuffer.data(idx_npot_0_if + 160);

    auto ta_xyyyyz_xxy_0 = pbuffer.data(idx_npot_0_if + 161);

    auto ta_xyyyyz_xxz_0 = pbuffer.data(idx_npot_0_if + 162);

    auto ta_xyyyyz_xyy_0 = pbuffer.data(idx_npot_0_if + 163);

    auto ta_xyyyyz_xyz_0 = pbuffer.data(idx_npot_0_if + 164);

    auto ta_xyyyyz_xzz_0 = pbuffer.data(idx_npot_0_if + 165);

    auto ta_xyyyyz_yyy_0 = pbuffer.data(idx_npot_0_if + 166);

    auto ta_xyyyyz_yyz_0 = pbuffer.data(idx_npot_0_if + 167);

    auto ta_xyyyyz_yzz_0 = pbuffer.data(idx_npot_0_if + 168);

    auto ta_xyyyyz_zzz_0 = pbuffer.data(idx_npot_0_if + 169);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xyyyy_xxx_0,  \
                             ta_xyyyy_xxx_1,  \
                             ta_xyyyy_xxy_0,  \
                             ta_xyyyy_xxy_1,  \
                             ta_xyyyy_xyy_0,  \
                             ta_xyyyy_xyy_1,  \
                             ta_xyyyyz_xxx_0, \
                             ta_xyyyyz_xxy_0, \
                             ta_xyyyyz_xxz_0, \
                             ta_xyyyyz_xyy_0, \
                             ta_xyyyyz_xyz_0, \
                             ta_xyyyyz_xzz_0, \
                             ta_xyyyyz_yyy_0, \
                             ta_xyyyyz_yyz_0, \
                             ta_xyyyyz_yzz_0, \
                             ta_xyyyyz_zzz_0, \
                             ta_yyyyz_xxz_0,  \
                             ta_yyyyz_xxz_1,  \
                             ta_yyyyz_xyz_0,  \
                             ta_yyyyz_xyz_1,  \
                             ta_yyyyz_xz_0,   \
                             ta_yyyyz_xz_1,   \
                             ta_yyyyz_xzz_0,  \
                             ta_yyyyz_xzz_1,  \
                             ta_yyyyz_yyy_0,  \
                             ta_yyyyz_yyy_1,  \
                             ta_yyyyz_yyz_0,  \
                             ta_yyyyz_yyz_1,  \
                             ta_yyyyz_yz_0,   \
                             ta_yyyyz_yz_1,   \
                             ta_yyyyz_yzz_0,  \
                             ta_yyyyz_yzz_1,  \
                             ta_yyyyz_zz_0,   \
                             ta_yyyyz_zz_1,   \
                             ta_yyyyz_zzz_0,  \
                             ta_yyyyz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyz_xxx_0[i] = ta_xyyyy_xxx_0[i] * pa_z[i] - ta_xyyyy_xxx_1[i] * pc_z[i];

        ta_xyyyyz_xxy_0[i] = ta_xyyyy_xxy_0[i] * pa_z[i] - ta_xyyyy_xxy_1[i] * pc_z[i];

        ta_xyyyyz_xxz_0[i] =
            2.0 * ta_yyyyz_xz_0[i] * fe_0 - 2.0 * ta_yyyyz_xz_1[i] * fe_0 + ta_yyyyz_xxz_0[i] * pa_x[i] - ta_yyyyz_xxz_1[i] * pc_x[i];

        ta_xyyyyz_xyy_0[i] = ta_xyyyy_xyy_0[i] * pa_z[i] - ta_xyyyy_xyy_1[i] * pc_z[i];

        ta_xyyyyz_xyz_0[i] = ta_yyyyz_yz_0[i] * fe_0 - ta_yyyyz_yz_1[i] * fe_0 + ta_yyyyz_xyz_0[i] * pa_x[i] - ta_yyyyz_xyz_1[i] * pc_x[i];

        ta_xyyyyz_xzz_0[i] = ta_yyyyz_zz_0[i] * fe_0 - ta_yyyyz_zz_1[i] * fe_0 + ta_yyyyz_xzz_0[i] * pa_x[i] - ta_yyyyz_xzz_1[i] * pc_x[i];

        ta_xyyyyz_yyy_0[i] = ta_yyyyz_yyy_0[i] * pa_x[i] - ta_yyyyz_yyy_1[i] * pc_x[i];

        ta_xyyyyz_yyz_0[i] = ta_yyyyz_yyz_0[i] * pa_x[i] - ta_yyyyz_yyz_1[i] * pc_x[i];

        ta_xyyyyz_yzz_0[i] = ta_yyyyz_yzz_0[i] * pa_x[i] - ta_yyyyz_yzz_1[i] * pc_x[i];

        ta_xyyyyz_zzz_0[i] = ta_yyyyz_zzz_0[i] * pa_x[i] - ta_yyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 170-180 components of targeted buffer : IF

    auto ta_xyyyzz_xxx_0 = pbuffer.data(idx_npot_0_if + 170);

    auto ta_xyyyzz_xxy_0 = pbuffer.data(idx_npot_0_if + 171);

    auto ta_xyyyzz_xxz_0 = pbuffer.data(idx_npot_0_if + 172);

    auto ta_xyyyzz_xyy_0 = pbuffer.data(idx_npot_0_if + 173);

    auto ta_xyyyzz_xyz_0 = pbuffer.data(idx_npot_0_if + 174);

    auto ta_xyyyzz_xzz_0 = pbuffer.data(idx_npot_0_if + 175);

    auto ta_xyyyzz_yyy_0 = pbuffer.data(idx_npot_0_if + 176);

    auto ta_xyyyzz_yyz_0 = pbuffer.data(idx_npot_0_if + 177);

    auto ta_xyyyzz_yzz_0 = pbuffer.data(idx_npot_0_if + 178);

    auto ta_xyyyzz_zzz_0 = pbuffer.data(idx_npot_0_if + 179);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xyyyzz_xxx_0, \
                             ta_xyyyzz_xxy_0, \
                             ta_xyyyzz_xxz_0, \
                             ta_xyyyzz_xyy_0, \
                             ta_xyyyzz_xyz_0, \
                             ta_xyyyzz_xzz_0, \
                             ta_xyyyzz_yyy_0, \
                             ta_xyyyzz_yyz_0, \
                             ta_xyyyzz_yzz_0, \
                             ta_xyyyzz_zzz_0, \
                             ta_yyyzz_xx_0,   \
                             ta_yyyzz_xx_1,   \
                             ta_yyyzz_xxx_0,  \
                             ta_yyyzz_xxx_1,  \
                             ta_yyyzz_xxy_0,  \
                             ta_yyyzz_xxy_1,  \
                             ta_yyyzz_xxz_0,  \
                             ta_yyyzz_xxz_1,  \
                             ta_yyyzz_xy_0,   \
                             ta_yyyzz_xy_1,   \
                             ta_yyyzz_xyy_0,  \
                             ta_yyyzz_xyy_1,  \
                             ta_yyyzz_xyz_0,  \
                             ta_yyyzz_xyz_1,  \
                             ta_yyyzz_xz_0,   \
                             ta_yyyzz_xz_1,   \
                             ta_yyyzz_xzz_0,  \
                             ta_yyyzz_xzz_1,  \
                             ta_yyyzz_yy_0,   \
                             ta_yyyzz_yy_1,   \
                             ta_yyyzz_yyy_0,  \
                             ta_yyyzz_yyy_1,  \
                             ta_yyyzz_yyz_0,  \
                             ta_yyyzz_yyz_1,  \
                             ta_yyyzz_yz_0,   \
                             ta_yyyzz_yz_1,   \
                             ta_yyyzz_yzz_0,  \
                             ta_yyyzz_yzz_1,  \
                             ta_yyyzz_zz_0,   \
                             ta_yyyzz_zz_1,   \
                             ta_yyyzz_zzz_0,  \
                             ta_yyyzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_xxx_0[i] =
            3.0 * ta_yyyzz_xx_0[i] * fe_0 - 3.0 * ta_yyyzz_xx_1[i] * fe_0 + ta_yyyzz_xxx_0[i] * pa_x[i] - ta_yyyzz_xxx_1[i] * pc_x[i];

        ta_xyyyzz_xxy_0[i] =
            2.0 * ta_yyyzz_xy_0[i] * fe_0 - 2.0 * ta_yyyzz_xy_1[i] * fe_0 + ta_yyyzz_xxy_0[i] * pa_x[i] - ta_yyyzz_xxy_1[i] * pc_x[i];

        ta_xyyyzz_xxz_0[i] =
            2.0 * ta_yyyzz_xz_0[i] * fe_0 - 2.0 * ta_yyyzz_xz_1[i] * fe_0 + ta_yyyzz_xxz_0[i] * pa_x[i] - ta_yyyzz_xxz_1[i] * pc_x[i];

        ta_xyyyzz_xyy_0[i] = ta_yyyzz_yy_0[i] * fe_0 - ta_yyyzz_yy_1[i] * fe_0 + ta_yyyzz_xyy_0[i] * pa_x[i] - ta_yyyzz_xyy_1[i] * pc_x[i];

        ta_xyyyzz_xyz_0[i] = ta_yyyzz_yz_0[i] * fe_0 - ta_yyyzz_yz_1[i] * fe_0 + ta_yyyzz_xyz_0[i] * pa_x[i] - ta_yyyzz_xyz_1[i] * pc_x[i];

        ta_xyyyzz_xzz_0[i] = ta_yyyzz_zz_0[i] * fe_0 - ta_yyyzz_zz_1[i] * fe_0 + ta_yyyzz_xzz_0[i] * pa_x[i] - ta_yyyzz_xzz_1[i] * pc_x[i];

        ta_xyyyzz_yyy_0[i] = ta_yyyzz_yyy_0[i] * pa_x[i] - ta_yyyzz_yyy_1[i] * pc_x[i];

        ta_xyyyzz_yyz_0[i] = ta_yyyzz_yyz_0[i] * pa_x[i] - ta_yyyzz_yyz_1[i] * pc_x[i];

        ta_xyyyzz_yzz_0[i] = ta_yyyzz_yzz_0[i] * pa_x[i] - ta_yyyzz_yzz_1[i] * pc_x[i];

        ta_xyyyzz_zzz_0[i] = ta_yyyzz_zzz_0[i] * pa_x[i] - ta_yyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 180-190 components of targeted buffer : IF

    auto ta_xyyzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 180);

    auto ta_xyyzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 181);

    auto ta_xyyzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 182);

    auto ta_xyyzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 183);

    auto ta_xyyzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 184);

    auto ta_xyyzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 185);

    auto ta_xyyzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 186);

    auto ta_xyyzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 187);

    auto ta_xyyzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 188);

    auto ta_xyyzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 189);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xyyzzz_xxx_0, \
                             ta_xyyzzz_xxy_0, \
                             ta_xyyzzz_xxz_0, \
                             ta_xyyzzz_xyy_0, \
                             ta_xyyzzz_xyz_0, \
                             ta_xyyzzz_xzz_0, \
                             ta_xyyzzz_yyy_0, \
                             ta_xyyzzz_yyz_0, \
                             ta_xyyzzz_yzz_0, \
                             ta_xyyzzz_zzz_0, \
                             ta_yyzzz_xx_0,   \
                             ta_yyzzz_xx_1,   \
                             ta_yyzzz_xxx_0,  \
                             ta_yyzzz_xxx_1,  \
                             ta_yyzzz_xxy_0,  \
                             ta_yyzzz_xxy_1,  \
                             ta_yyzzz_xxz_0,  \
                             ta_yyzzz_xxz_1,  \
                             ta_yyzzz_xy_0,   \
                             ta_yyzzz_xy_1,   \
                             ta_yyzzz_xyy_0,  \
                             ta_yyzzz_xyy_1,  \
                             ta_yyzzz_xyz_0,  \
                             ta_yyzzz_xyz_1,  \
                             ta_yyzzz_xz_0,   \
                             ta_yyzzz_xz_1,   \
                             ta_yyzzz_xzz_0,  \
                             ta_yyzzz_xzz_1,  \
                             ta_yyzzz_yy_0,   \
                             ta_yyzzz_yy_1,   \
                             ta_yyzzz_yyy_0,  \
                             ta_yyzzz_yyy_1,  \
                             ta_yyzzz_yyz_0,  \
                             ta_yyzzz_yyz_1,  \
                             ta_yyzzz_yz_0,   \
                             ta_yyzzz_yz_1,   \
                             ta_yyzzz_yzz_0,  \
                             ta_yyzzz_yzz_1,  \
                             ta_yyzzz_zz_0,   \
                             ta_yyzzz_zz_1,   \
                             ta_yyzzz_zzz_0,  \
                             ta_yyzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_xxx_0[i] =
            3.0 * ta_yyzzz_xx_0[i] * fe_0 - 3.0 * ta_yyzzz_xx_1[i] * fe_0 + ta_yyzzz_xxx_0[i] * pa_x[i] - ta_yyzzz_xxx_1[i] * pc_x[i];

        ta_xyyzzz_xxy_0[i] =
            2.0 * ta_yyzzz_xy_0[i] * fe_0 - 2.0 * ta_yyzzz_xy_1[i] * fe_0 + ta_yyzzz_xxy_0[i] * pa_x[i] - ta_yyzzz_xxy_1[i] * pc_x[i];

        ta_xyyzzz_xxz_0[i] =
            2.0 * ta_yyzzz_xz_0[i] * fe_0 - 2.0 * ta_yyzzz_xz_1[i] * fe_0 + ta_yyzzz_xxz_0[i] * pa_x[i] - ta_yyzzz_xxz_1[i] * pc_x[i];

        ta_xyyzzz_xyy_0[i] = ta_yyzzz_yy_0[i] * fe_0 - ta_yyzzz_yy_1[i] * fe_0 + ta_yyzzz_xyy_0[i] * pa_x[i] - ta_yyzzz_xyy_1[i] * pc_x[i];

        ta_xyyzzz_xyz_0[i] = ta_yyzzz_yz_0[i] * fe_0 - ta_yyzzz_yz_1[i] * fe_0 + ta_yyzzz_xyz_0[i] * pa_x[i] - ta_yyzzz_xyz_1[i] * pc_x[i];

        ta_xyyzzz_xzz_0[i] = ta_yyzzz_zz_0[i] * fe_0 - ta_yyzzz_zz_1[i] * fe_0 + ta_yyzzz_xzz_0[i] * pa_x[i] - ta_yyzzz_xzz_1[i] * pc_x[i];

        ta_xyyzzz_yyy_0[i] = ta_yyzzz_yyy_0[i] * pa_x[i] - ta_yyzzz_yyy_1[i] * pc_x[i];

        ta_xyyzzz_yyz_0[i] = ta_yyzzz_yyz_0[i] * pa_x[i] - ta_yyzzz_yyz_1[i] * pc_x[i];

        ta_xyyzzz_yzz_0[i] = ta_yyzzz_yzz_0[i] * pa_x[i] - ta_yyzzz_yzz_1[i] * pc_x[i];

        ta_xyyzzz_zzz_0[i] = ta_yyzzz_zzz_0[i] * pa_x[i] - ta_yyzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 190-200 components of targeted buffer : IF

    auto ta_xyzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 190);

    auto ta_xyzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 191);

    auto ta_xyzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 192);

    auto ta_xyzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 193);

    auto ta_xyzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 194);

    auto ta_xyzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 195);

    auto ta_xyzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 196);

    auto ta_xyzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 197);

    auto ta_xyzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 198);

    auto ta_xyzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 199);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xyzzzz_xxx_0, \
                             ta_xyzzzz_xxy_0, \
                             ta_xyzzzz_xxz_0, \
                             ta_xyzzzz_xyy_0, \
                             ta_xyzzzz_xyz_0, \
                             ta_xyzzzz_xzz_0, \
                             ta_xyzzzz_yyy_0, \
                             ta_xyzzzz_yyz_0, \
                             ta_xyzzzz_yzz_0, \
                             ta_xyzzzz_zzz_0, \
                             ta_xzzzz_xxx_0,  \
                             ta_xzzzz_xxx_1,  \
                             ta_xzzzz_xxz_0,  \
                             ta_xzzzz_xxz_1,  \
                             ta_xzzzz_xzz_0,  \
                             ta_xzzzz_xzz_1,  \
                             ta_yzzzz_xxy_0,  \
                             ta_yzzzz_xxy_1,  \
                             ta_yzzzz_xy_0,   \
                             ta_yzzzz_xy_1,   \
                             ta_yzzzz_xyy_0,  \
                             ta_yzzzz_xyy_1,  \
                             ta_yzzzz_xyz_0,  \
                             ta_yzzzz_xyz_1,  \
                             ta_yzzzz_yy_0,   \
                             ta_yzzzz_yy_1,   \
                             ta_yzzzz_yyy_0,  \
                             ta_yzzzz_yyy_1,  \
                             ta_yzzzz_yyz_0,  \
                             ta_yzzzz_yyz_1,  \
                             ta_yzzzz_yz_0,   \
                             ta_yzzzz_yz_1,   \
                             ta_yzzzz_yzz_0,  \
                             ta_yzzzz_yzz_1,  \
                             ta_yzzzz_zzz_0,  \
                             ta_yzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzzz_xxx_0[i] = ta_xzzzz_xxx_0[i] * pa_y[i] - ta_xzzzz_xxx_1[i] * pc_y[i];

        ta_xyzzzz_xxy_0[i] =
            2.0 * ta_yzzzz_xy_0[i] * fe_0 - 2.0 * ta_yzzzz_xy_1[i] * fe_0 + ta_yzzzz_xxy_0[i] * pa_x[i] - ta_yzzzz_xxy_1[i] * pc_x[i];

        ta_xyzzzz_xxz_0[i] = ta_xzzzz_xxz_0[i] * pa_y[i] - ta_xzzzz_xxz_1[i] * pc_y[i];

        ta_xyzzzz_xyy_0[i] = ta_yzzzz_yy_0[i] * fe_0 - ta_yzzzz_yy_1[i] * fe_0 + ta_yzzzz_xyy_0[i] * pa_x[i] - ta_yzzzz_xyy_1[i] * pc_x[i];

        ta_xyzzzz_xyz_0[i] = ta_yzzzz_yz_0[i] * fe_0 - ta_yzzzz_yz_1[i] * fe_0 + ta_yzzzz_xyz_0[i] * pa_x[i] - ta_yzzzz_xyz_1[i] * pc_x[i];

        ta_xyzzzz_xzz_0[i] = ta_xzzzz_xzz_0[i] * pa_y[i] - ta_xzzzz_xzz_1[i] * pc_y[i];

        ta_xyzzzz_yyy_0[i] = ta_yzzzz_yyy_0[i] * pa_x[i] - ta_yzzzz_yyy_1[i] * pc_x[i];

        ta_xyzzzz_yyz_0[i] = ta_yzzzz_yyz_0[i] * pa_x[i] - ta_yzzzz_yyz_1[i] * pc_x[i];

        ta_xyzzzz_yzz_0[i] = ta_yzzzz_yzz_0[i] * pa_x[i] - ta_yzzzz_yzz_1[i] * pc_x[i];

        ta_xyzzzz_zzz_0[i] = ta_yzzzz_zzz_0[i] * pa_x[i] - ta_yzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 200-210 components of targeted buffer : IF

    auto ta_xzzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 200);

    auto ta_xzzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 201);

    auto ta_xzzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 202);

    auto ta_xzzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 203);

    auto ta_xzzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 204);

    auto ta_xzzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 205);

    auto ta_xzzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 206);

    auto ta_xzzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 207);

    auto ta_xzzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 208);

    auto ta_xzzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 209);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xzzzzz_xxx_0, \
                             ta_xzzzzz_xxy_0, \
                             ta_xzzzzz_xxz_0, \
                             ta_xzzzzz_xyy_0, \
                             ta_xzzzzz_xyz_0, \
                             ta_xzzzzz_xzz_0, \
                             ta_xzzzzz_yyy_0, \
                             ta_xzzzzz_yyz_0, \
                             ta_xzzzzz_yzz_0, \
                             ta_xzzzzz_zzz_0, \
                             ta_zzzzz_xx_0,   \
                             ta_zzzzz_xx_1,   \
                             ta_zzzzz_xxx_0,  \
                             ta_zzzzz_xxx_1,  \
                             ta_zzzzz_xxy_0,  \
                             ta_zzzzz_xxy_1,  \
                             ta_zzzzz_xxz_0,  \
                             ta_zzzzz_xxz_1,  \
                             ta_zzzzz_xy_0,   \
                             ta_zzzzz_xy_1,   \
                             ta_zzzzz_xyy_0,  \
                             ta_zzzzz_xyy_1,  \
                             ta_zzzzz_xyz_0,  \
                             ta_zzzzz_xyz_1,  \
                             ta_zzzzz_xz_0,   \
                             ta_zzzzz_xz_1,   \
                             ta_zzzzz_xzz_0,  \
                             ta_zzzzz_xzz_1,  \
                             ta_zzzzz_yy_0,   \
                             ta_zzzzz_yy_1,   \
                             ta_zzzzz_yyy_0,  \
                             ta_zzzzz_yyy_1,  \
                             ta_zzzzz_yyz_0,  \
                             ta_zzzzz_yyz_1,  \
                             ta_zzzzz_yz_0,   \
                             ta_zzzzz_yz_1,   \
                             ta_zzzzz_yzz_0,  \
                             ta_zzzzz_yzz_1,  \
                             ta_zzzzz_zz_0,   \
                             ta_zzzzz_zz_1,   \
                             ta_zzzzz_zzz_0,  \
                             ta_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_xxx_0[i] =
            3.0 * ta_zzzzz_xx_0[i] * fe_0 - 3.0 * ta_zzzzz_xx_1[i] * fe_0 + ta_zzzzz_xxx_0[i] * pa_x[i] - ta_zzzzz_xxx_1[i] * pc_x[i];

        ta_xzzzzz_xxy_0[i] =
            2.0 * ta_zzzzz_xy_0[i] * fe_0 - 2.0 * ta_zzzzz_xy_1[i] * fe_0 + ta_zzzzz_xxy_0[i] * pa_x[i] - ta_zzzzz_xxy_1[i] * pc_x[i];

        ta_xzzzzz_xxz_0[i] =
            2.0 * ta_zzzzz_xz_0[i] * fe_0 - 2.0 * ta_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xxz_0[i] * pa_x[i] - ta_zzzzz_xxz_1[i] * pc_x[i];

        ta_xzzzzz_xyy_0[i] = ta_zzzzz_yy_0[i] * fe_0 - ta_zzzzz_yy_1[i] * fe_0 + ta_zzzzz_xyy_0[i] * pa_x[i] - ta_zzzzz_xyy_1[i] * pc_x[i];

        ta_xzzzzz_xyz_0[i] = ta_zzzzz_yz_0[i] * fe_0 - ta_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_xyz_0[i] * pa_x[i] - ta_zzzzz_xyz_1[i] * pc_x[i];

        ta_xzzzzz_xzz_0[i] = ta_zzzzz_zz_0[i] * fe_0 - ta_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_xzz_0[i] * pa_x[i] - ta_zzzzz_xzz_1[i] * pc_x[i];

        ta_xzzzzz_yyy_0[i] = ta_zzzzz_yyy_0[i] * pa_x[i] - ta_zzzzz_yyy_1[i] * pc_x[i];

        ta_xzzzzz_yyz_0[i] = ta_zzzzz_yyz_0[i] * pa_x[i] - ta_zzzzz_yyz_1[i] * pc_x[i];

        ta_xzzzzz_yzz_0[i] = ta_zzzzz_yzz_0[i] * pa_x[i] - ta_zzzzz_yzz_1[i] * pc_x[i];

        ta_xzzzzz_zzz_0[i] = ta_zzzzz_zzz_0[i] * pa_x[i] - ta_zzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : IF

    auto ta_yyyyyy_xxx_0 = pbuffer.data(idx_npot_0_if + 210);

    auto ta_yyyyyy_xxy_0 = pbuffer.data(idx_npot_0_if + 211);

    auto ta_yyyyyy_xxz_0 = pbuffer.data(idx_npot_0_if + 212);

    auto ta_yyyyyy_xyy_0 = pbuffer.data(idx_npot_0_if + 213);

    auto ta_yyyyyy_xyz_0 = pbuffer.data(idx_npot_0_if + 214);

    auto ta_yyyyyy_xzz_0 = pbuffer.data(idx_npot_0_if + 215);

    auto ta_yyyyyy_yyy_0 = pbuffer.data(idx_npot_0_if + 216);

    auto ta_yyyyyy_yyz_0 = pbuffer.data(idx_npot_0_if + 217);

    auto ta_yyyyyy_yzz_0 = pbuffer.data(idx_npot_0_if + 218);

    auto ta_yyyyyy_zzz_0 = pbuffer.data(idx_npot_0_if + 219);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta_yyyy_xxx_0,   \
                             ta_yyyy_xxx_1,   \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xxz_0,   \
                             ta_yyyy_xxz_1,   \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_xyz_0,   \
                             ta_yyyy_xyz_1,   \
                             ta_yyyy_xzz_0,   \
                             ta_yyyy_xzz_1,   \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyy_yyz_0,   \
                             ta_yyyy_yyz_1,   \
                             ta_yyyy_yzz_0,   \
                             ta_yyyy_yzz_1,   \
                             ta_yyyy_zzz_0,   \
                             ta_yyyy_zzz_1,   \
                             ta_yyyyy_xx_0,   \
                             ta_yyyyy_xx_1,   \
                             ta_yyyyy_xxx_0,  \
                             ta_yyyyy_xxx_1,  \
                             ta_yyyyy_xxy_0,  \
                             ta_yyyyy_xxy_1,  \
                             ta_yyyyy_xxz_0,  \
                             ta_yyyyy_xxz_1,  \
                             ta_yyyyy_xy_0,   \
                             ta_yyyyy_xy_1,   \
                             ta_yyyyy_xyy_0,  \
                             ta_yyyyy_xyy_1,  \
                             ta_yyyyy_xyz_0,  \
                             ta_yyyyy_xyz_1,  \
                             ta_yyyyy_xz_0,   \
                             ta_yyyyy_xz_1,   \
                             ta_yyyyy_xzz_0,  \
                             ta_yyyyy_xzz_1,  \
                             ta_yyyyy_yy_0,   \
                             ta_yyyyy_yy_1,   \
                             ta_yyyyy_yyy_0,  \
                             ta_yyyyy_yyy_1,  \
                             ta_yyyyy_yyz_0,  \
                             ta_yyyyy_yyz_1,  \
                             ta_yyyyy_yz_0,   \
                             ta_yyyyy_yz_1,   \
                             ta_yyyyy_yzz_0,  \
                             ta_yyyyy_yzz_1,  \
                             ta_yyyyy_zz_0,   \
                             ta_yyyyy_zz_1,   \
                             ta_yyyyy_zzz_0,  \
                             ta_yyyyy_zzz_1,  \
                             ta_yyyyyy_xxx_0, \
                             ta_yyyyyy_xxy_0, \
                             ta_yyyyyy_xxz_0, \
                             ta_yyyyyy_xyy_0, \
                             ta_yyyyyy_xyz_0, \
                             ta_yyyyyy_xzz_0, \
                             ta_yyyyyy_yyy_0, \
                             ta_yyyyyy_yyz_0, \
                             ta_yyyyyy_yzz_0, \
                             ta_yyyyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_xxx_0[i] =
            5.0 * ta_yyyy_xxx_0[i] * fe_0 - 5.0 * ta_yyyy_xxx_1[i] * fe_0 + ta_yyyyy_xxx_0[i] * pa_y[i] - ta_yyyyy_xxx_1[i] * pc_y[i];

        ta_yyyyyy_xxy_0[i] = 5.0 * ta_yyyy_xxy_0[i] * fe_0 - 5.0 * ta_yyyy_xxy_1[i] * fe_0 + ta_yyyyy_xx_0[i] * fe_0 - ta_yyyyy_xx_1[i] * fe_0 +
                             ta_yyyyy_xxy_0[i] * pa_y[i] - ta_yyyyy_xxy_1[i] * pc_y[i];

        ta_yyyyyy_xxz_0[i] =
            5.0 * ta_yyyy_xxz_0[i] * fe_0 - 5.0 * ta_yyyy_xxz_1[i] * fe_0 + ta_yyyyy_xxz_0[i] * pa_y[i] - ta_yyyyy_xxz_1[i] * pc_y[i];

        ta_yyyyyy_xyy_0[i] = 5.0 * ta_yyyy_xyy_0[i] * fe_0 - 5.0 * ta_yyyy_xyy_1[i] * fe_0 + 2.0 * ta_yyyyy_xy_0[i] * fe_0 -
                             2.0 * ta_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xyy_0[i] * pa_y[i] - ta_yyyyy_xyy_1[i] * pc_y[i];

        ta_yyyyyy_xyz_0[i] = 5.0 * ta_yyyy_xyz_0[i] * fe_0 - 5.0 * ta_yyyy_xyz_1[i] * fe_0 + ta_yyyyy_xz_0[i] * fe_0 - ta_yyyyy_xz_1[i] * fe_0 +
                             ta_yyyyy_xyz_0[i] * pa_y[i] - ta_yyyyy_xyz_1[i] * pc_y[i];

        ta_yyyyyy_xzz_0[i] =
            5.0 * ta_yyyy_xzz_0[i] * fe_0 - 5.0 * ta_yyyy_xzz_1[i] * fe_0 + ta_yyyyy_xzz_0[i] * pa_y[i] - ta_yyyyy_xzz_1[i] * pc_y[i];

        ta_yyyyyy_yyy_0[i] = 5.0 * ta_yyyy_yyy_0[i] * fe_0 - 5.0 * ta_yyyy_yyy_1[i] * fe_0 + 3.0 * ta_yyyyy_yy_0[i] * fe_0 -
                             3.0 * ta_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_yyy_0[i] * pa_y[i] - ta_yyyyy_yyy_1[i] * pc_y[i];

        ta_yyyyyy_yyz_0[i] = 5.0 * ta_yyyy_yyz_0[i] * fe_0 - 5.0 * ta_yyyy_yyz_1[i] * fe_0 + 2.0 * ta_yyyyy_yz_0[i] * fe_0 -
                             2.0 * ta_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_yyz_0[i] * pa_y[i] - ta_yyyyy_yyz_1[i] * pc_y[i];

        ta_yyyyyy_yzz_0[i] = 5.0 * ta_yyyy_yzz_0[i] * fe_0 - 5.0 * ta_yyyy_yzz_1[i] * fe_0 + ta_yyyyy_zz_0[i] * fe_0 - ta_yyyyy_zz_1[i] * fe_0 +
                             ta_yyyyy_yzz_0[i] * pa_y[i] - ta_yyyyy_yzz_1[i] * pc_y[i];

        ta_yyyyyy_zzz_0[i] =
            5.0 * ta_yyyy_zzz_0[i] * fe_0 - 5.0 * ta_yyyy_zzz_1[i] * fe_0 + ta_yyyyy_zzz_0[i] * pa_y[i] - ta_yyyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 220-230 components of targeted buffer : IF

    auto ta_yyyyyz_xxx_0 = pbuffer.data(idx_npot_0_if + 220);

    auto ta_yyyyyz_xxy_0 = pbuffer.data(idx_npot_0_if + 221);

    auto ta_yyyyyz_xxz_0 = pbuffer.data(idx_npot_0_if + 222);

    auto ta_yyyyyz_xyy_0 = pbuffer.data(idx_npot_0_if + 223);

    auto ta_yyyyyz_xyz_0 = pbuffer.data(idx_npot_0_if + 224);

    auto ta_yyyyyz_xzz_0 = pbuffer.data(idx_npot_0_if + 225);

    auto ta_yyyyyz_yyy_0 = pbuffer.data(idx_npot_0_if + 226);

    auto ta_yyyyyz_yyz_0 = pbuffer.data(idx_npot_0_if + 227);

    auto ta_yyyyyz_yzz_0 = pbuffer.data(idx_npot_0_if + 228);

    auto ta_yyyyyz_zzz_0 = pbuffer.data(idx_npot_0_if + 229);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyyyy_xxx_0,  \
                             ta_yyyyy_xxx_1,  \
                             ta_yyyyy_xxy_0,  \
                             ta_yyyyy_xxy_1,  \
                             ta_yyyyy_xy_0,   \
                             ta_yyyyy_xy_1,   \
                             ta_yyyyy_xyy_0,  \
                             ta_yyyyy_xyy_1,  \
                             ta_yyyyy_xyz_0,  \
                             ta_yyyyy_xyz_1,  \
                             ta_yyyyy_yy_0,   \
                             ta_yyyyy_yy_1,   \
                             ta_yyyyy_yyy_0,  \
                             ta_yyyyy_yyy_1,  \
                             ta_yyyyy_yyz_0,  \
                             ta_yyyyy_yyz_1,  \
                             ta_yyyyy_yz_0,   \
                             ta_yyyyy_yz_1,   \
                             ta_yyyyy_yzz_0,  \
                             ta_yyyyy_yzz_1,  \
                             ta_yyyyyz_xxx_0, \
                             ta_yyyyyz_xxy_0, \
                             ta_yyyyyz_xxz_0, \
                             ta_yyyyyz_xyy_0, \
                             ta_yyyyyz_xyz_0, \
                             ta_yyyyyz_xzz_0, \
                             ta_yyyyyz_yyy_0, \
                             ta_yyyyyz_yyz_0, \
                             ta_yyyyyz_yzz_0, \
                             ta_yyyyyz_zzz_0, \
                             ta_yyyyz_xxz_0,  \
                             ta_yyyyz_xxz_1,  \
                             ta_yyyyz_xzz_0,  \
                             ta_yyyyz_xzz_1,  \
                             ta_yyyyz_zzz_0,  \
                             ta_yyyyz_zzz_1,  \
                             ta_yyyz_xxz_0,   \
                             ta_yyyz_xxz_1,   \
                             ta_yyyz_xzz_0,   \
                             ta_yyyz_xzz_1,   \
                             ta_yyyz_zzz_0,   \
                             ta_yyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_xxx_0[i] = ta_yyyyy_xxx_0[i] * pa_z[i] - ta_yyyyy_xxx_1[i] * pc_z[i];

        ta_yyyyyz_xxy_0[i] = ta_yyyyy_xxy_0[i] * pa_z[i] - ta_yyyyy_xxy_1[i] * pc_z[i];

        ta_yyyyyz_xxz_0[i] =
            4.0 * ta_yyyz_xxz_0[i] * fe_0 - 4.0 * ta_yyyz_xxz_1[i] * fe_0 + ta_yyyyz_xxz_0[i] * pa_y[i] - ta_yyyyz_xxz_1[i] * pc_y[i];

        ta_yyyyyz_xyy_0[i] = ta_yyyyy_xyy_0[i] * pa_z[i] - ta_yyyyy_xyy_1[i] * pc_z[i];

        ta_yyyyyz_xyz_0[i] = ta_yyyyy_xy_0[i] * fe_0 - ta_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xyz_0[i] * pa_z[i] - ta_yyyyy_xyz_1[i] * pc_z[i];

        ta_yyyyyz_xzz_0[i] =
            4.0 * ta_yyyz_xzz_0[i] * fe_0 - 4.0 * ta_yyyz_xzz_1[i] * fe_0 + ta_yyyyz_xzz_0[i] * pa_y[i] - ta_yyyyz_xzz_1[i] * pc_y[i];

        ta_yyyyyz_yyy_0[i] = ta_yyyyy_yyy_0[i] * pa_z[i] - ta_yyyyy_yyy_1[i] * pc_z[i];

        ta_yyyyyz_yyz_0[i] = ta_yyyyy_yy_0[i] * fe_0 - ta_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_yyz_0[i] * pa_z[i] - ta_yyyyy_yyz_1[i] * pc_z[i];

        ta_yyyyyz_yzz_0[i] =
            2.0 * ta_yyyyy_yz_0[i] * fe_0 - 2.0 * ta_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_yzz_0[i] * pa_z[i] - ta_yyyyy_yzz_1[i] * pc_z[i];

        ta_yyyyyz_zzz_0[i] =
            4.0 * ta_yyyz_zzz_0[i] * fe_0 - 4.0 * ta_yyyz_zzz_1[i] * fe_0 + ta_yyyyz_zzz_0[i] * pa_y[i] - ta_yyyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 230-240 components of targeted buffer : IF

    auto ta_yyyyzz_xxx_0 = pbuffer.data(idx_npot_0_if + 230);

    auto ta_yyyyzz_xxy_0 = pbuffer.data(idx_npot_0_if + 231);

    auto ta_yyyyzz_xxz_0 = pbuffer.data(idx_npot_0_if + 232);

    auto ta_yyyyzz_xyy_0 = pbuffer.data(idx_npot_0_if + 233);

    auto ta_yyyyzz_xyz_0 = pbuffer.data(idx_npot_0_if + 234);

    auto ta_yyyyzz_xzz_0 = pbuffer.data(idx_npot_0_if + 235);

    auto ta_yyyyzz_yyy_0 = pbuffer.data(idx_npot_0_if + 236);

    auto ta_yyyyzz_yyz_0 = pbuffer.data(idx_npot_0_if + 237);

    auto ta_yyyyzz_yzz_0 = pbuffer.data(idx_npot_0_if + 238);

    auto ta_yyyyzz_zzz_0 = pbuffer.data(idx_npot_0_if + 239);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyyz_xxy_0,  \
                             ta_yyyyz_xxy_1,  \
                             ta_yyyyz_xyy_0,  \
                             ta_yyyyz_xyy_1,  \
                             ta_yyyyz_yyy_0,  \
                             ta_yyyyz_yyy_1,  \
                             ta_yyyyzz_xxx_0, \
                             ta_yyyyzz_xxy_0, \
                             ta_yyyyzz_xxz_0, \
                             ta_yyyyzz_xyy_0, \
                             ta_yyyyzz_xyz_0, \
                             ta_yyyyzz_xzz_0, \
                             ta_yyyyzz_yyy_0, \
                             ta_yyyyzz_yyz_0, \
                             ta_yyyyzz_yzz_0, \
                             ta_yyyyzz_zzz_0, \
                             ta_yyyzz_xxx_0,  \
                             ta_yyyzz_xxx_1,  \
                             ta_yyyzz_xxz_0,  \
                             ta_yyyzz_xxz_1,  \
                             ta_yyyzz_xyz_0,  \
                             ta_yyyzz_xyz_1,  \
                             ta_yyyzz_xz_0,   \
                             ta_yyyzz_xz_1,   \
                             ta_yyyzz_xzz_0,  \
                             ta_yyyzz_xzz_1,  \
                             ta_yyyzz_yyz_0,  \
                             ta_yyyzz_yyz_1,  \
                             ta_yyyzz_yz_0,   \
                             ta_yyyzz_yz_1,   \
                             ta_yyyzz_yzz_0,  \
                             ta_yyyzz_yzz_1,  \
                             ta_yyyzz_zz_0,   \
                             ta_yyyzz_zz_1,   \
                             ta_yyyzz_zzz_0,  \
                             ta_yyyzz_zzz_1,  \
                             ta_yyzz_xxx_0,   \
                             ta_yyzz_xxx_1,   \
                             ta_yyzz_xxz_0,   \
                             ta_yyzz_xxz_1,   \
                             ta_yyzz_xyz_0,   \
                             ta_yyzz_xyz_1,   \
                             ta_yyzz_xzz_0,   \
                             ta_yyzz_xzz_1,   \
                             ta_yyzz_yyz_0,   \
                             ta_yyzz_yyz_1,   \
                             ta_yyzz_yzz_0,   \
                             ta_yyzz_yzz_1,   \
                             ta_yyzz_zzz_0,   \
                             ta_yyzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_xxx_0[i] =
            3.0 * ta_yyzz_xxx_0[i] * fe_0 - 3.0 * ta_yyzz_xxx_1[i] * fe_0 + ta_yyyzz_xxx_0[i] * pa_y[i] - ta_yyyzz_xxx_1[i] * pc_y[i];

        ta_yyyyzz_xxy_0[i] = ta_yyyy_xxy_0[i] * fe_0 - ta_yyyy_xxy_1[i] * fe_0 + ta_yyyyz_xxy_0[i] * pa_z[i] - ta_yyyyz_xxy_1[i] * pc_z[i];

        ta_yyyyzz_xxz_0[i] =
            3.0 * ta_yyzz_xxz_0[i] * fe_0 - 3.0 * ta_yyzz_xxz_1[i] * fe_0 + ta_yyyzz_xxz_0[i] * pa_y[i] - ta_yyyzz_xxz_1[i] * pc_y[i];

        ta_yyyyzz_xyy_0[i] = ta_yyyy_xyy_0[i] * fe_0 - ta_yyyy_xyy_1[i] * fe_0 + ta_yyyyz_xyy_0[i] * pa_z[i] - ta_yyyyz_xyy_1[i] * pc_z[i];

        ta_yyyyzz_xyz_0[i] = 3.0 * ta_yyzz_xyz_0[i] * fe_0 - 3.0 * ta_yyzz_xyz_1[i] * fe_0 + ta_yyyzz_xz_0[i] * fe_0 - ta_yyyzz_xz_1[i] * fe_0 +
                             ta_yyyzz_xyz_0[i] * pa_y[i] - ta_yyyzz_xyz_1[i] * pc_y[i];

        ta_yyyyzz_xzz_0[i] =
            3.0 * ta_yyzz_xzz_0[i] * fe_0 - 3.0 * ta_yyzz_xzz_1[i] * fe_0 + ta_yyyzz_xzz_0[i] * pa_y[i] - ta_yyyzz_xzz_1[i] * pc_y[i];

        ta_yyyyzz_yyy_0[i] = ta_yyyy_yyy_0[i] * fe_0 - ta_yyyy_yyy_1[i] * fe_0 + ta_yyyyz_yyy_0[i] * pa_z[i] - ta_yyyyz_yyy_1[i] * pc_z[i];

        ta_yyyyzz_yyz_0[i] = 3.0 * ta_yyzz_yyz_0[i] * fe_0 - 3.0 * ta_yyzz_yyz_1[i] * fe_0 + 2.0 * ta_yyyzz_yz_0[i] * fe_0 -
                             2.0 * ta_yyyzz_yz_1[i] * fe_0 + ta_yyyzz_yyz_0[i] * pa_y[i] - ta_yyyzz_yyz_1[i] * pc_y[i];

        ta_yyyyzz_yzz_0[i] = 3.0 * ta_yyzz_yzz_0[i] * fe_0 - 3.0 * ta_yyzz_yzz_1[i] * fe_0 + ta_yyyzz_zz_0[i] * fe_0 - ta_yyyzz_zz_1[i] * fe_0 +
                             ta_yyyzz_yzz_0[i] * pa_y[i] - ta_yyyzz_yzz_1[i] * pc_y[i];

        ta_yyyyzz_zzz_0[i] =
            3.0 * ta_yyzz_zzz_0[i] * fe_0 - 3.0 * ta_yyzz_zzz_1[i] * fe_0 + ta_yyyzz_zzz_0[i] * pa_y[i] - ta_yyyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 240-250 components of targeted buffer : IF

    auto ta_yyyzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 240);

    auto ta_yyyzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 241);

    auto ta_yyyzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 242);

    auto ta_yyyzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 243);

    auto ta_yyyzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 244);

    auto ta_yyyzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 245);

    auto ta_yyyzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 246);

    auto ta_yyyzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 247);

    auto ta_yyyzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 248);

    auto ta_yyyzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 249);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyyz_xxy_0,   \
                             ta_yyyz_xxy_1,   \
                             ta_yyyz_xyy_0,   \
                             ta_yyyz_xyy_1,   \
                             ta_yyyz_yyy_0,   \
                             ta_yyyz_yyy_1,   \
                             ta_yyyzz_xxy_0,  \
                             ta_yyyzz_xxy_1,  \
                             ta_yyyzz_xyy_0,  \
                             ta_yyyzz_xyy_1,  \
                             ta_yyyzz_yyy_0,  \
                             ta_yyyzz_yyy_1,  \
                             ta_yyyzzz_xxx_0, \
                             ta_yyyzzz_xxy_0, \
                             ta_yyyzzz_xxz_0, \
                             ta_yyyzzz_xyy_0, \
                             ta_yyyzzz_xyz_0, \
                             ta_yyyzzz_xzz_0, \
                             ta_yyyzzz_yyy_0, \
                             ta_yyyzzz_yyz_0, \
                             ta_yyyzzz_yzz_0, \
                             ta_yyyzzz_zzz_0, \
                             ta_yyzzz_xxx_0,  \
                             ta_yyzzz_xxx_1,  \
                             ta_yyzzz_xxz_0,  \
                             ta_yyzzz_xxz_1,  \
                             ta_yyzzz_xyz_0,  \
                             ta_yyzzz_xyz_1,  \
                             ta_yyzzz_xz_0,   \
                             ta_yyzzz_xz_1,   \
                             ta_yyzzz_xzz_0,  \
                             ta_yyzzz_xzz_1,  \
                             ta_yyzzz_yyz_0,  \
                             ta_yyzzz_yyz_1,  \
                             ta_yyzzz_yz_0,   \
                             ta_yyzzz_yz_1,   \
                             ta_yyzzz_yzz_0,  \
                             ta_yyzzz_yzz_1,  \
                             ta_yyzzz_zz_0,   \
                             ta_yyzzz_zz_1,   \
                             ta_yyzzz_zzz_0,  \
                             ta_yyzzz_zzz_1,  \
                             ta_yzzz_xxx_0,   \
                             ta_yzzz_xxx_1,   \
                             ta_yzzz_xxz_0,   \
                             ta_yzzz_xxz_1,   \
                             ta_yzzz_xyz_0,   \
                             ta_yzzz_xyz_1,   \
                             ta_yzzz_xzz_0,   \
                             ta_yzzz_xzz_1,   \
                             ta_yzzz_yyz_0,   \
                             ta_yzzz_yyz_1,   \
                             ta_yzzz_yzz_0,   \
                             ta_yzzz_yzz_1,   \
                             ta_yzzz_zzz_0,   \
                             ta_yzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_xxx_0[i] =
            2.0 * ta_yzzz_xxx_0[i] * fe_0 - 2.0 * ta_yzzz_xxx_1[i] * fe_0 + ta_yyzzz_xxx_0[i] * pa_y[i] - ta_yyzzz_xxx_1[i] * pc_y[i];

        ta_yyyzzz_xxy_0[i] =
            2.0 * ta_yyyz_xxy_0[i] * fe_0 - 2.0 * ta_yyyz_xxy_1[i] * fe_0 + ta_yyyzz_xxy_0[i] * pa_z[i] - ta_yyyzz_xxy_1[i] * pc_z[i];

        ta_yyyzzz_xxz_0[i] =
            2.0 * ta_yzzz_xxz_0[i] * fe_0 - 2.0 * ta_yzzz_xxz_1[i] * fe_0 + ta_yyzzz_xxz_0[i] * pa_y[i] - ta_yyzzz_xxz_1[i] * pc_y[i];

        ta_yyyzzz_xyy_0[i] =
            2.0 * ta_yyyz_xyy_0[i] * fe_0 - 2.0 * ta_yyyz_xyy_1[i] * fe_0 + ta_yyyzz_xyy_0[i] * pa_z[i] - ta_yyyzz_xyy_1[i] * pc_z[i];

        ta_yyyzzz_xyz_0[i] = 2.0 * ta_yzzz_xyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyz_1[i] * fe_0 + ta_yyzzz_xz_0[i] * fe_0 - ta_yyzzz_xz_1[i] * fe_0 +
                             ta_yyzzz_xyz_0[i] * pa_y[i] - ta_yyzzz_xyz_1[i] * pc_y[i];

        ta_yyyzzz_xzz_0[i] =
            2.0 * ta_yzzz_xzz_0[i] * fe_0 - 2.0 * ta_yzzz_xzz_1[i] * fe_0 + ta_yyzzz_xzz_0[i] * pa_y[i] - ta_yyzzz_xzz_1[i] * pc_y[i];

        ta_yyyzzz_yyy_0[i] =
            2.0 * ta_yyyz_yyy_0[i] * fe_0 - 2.0 * ta_yyyz_yyy_1[i] * fe_0 + ta_yyyzz_yyy_0[i] * pa_z[i] - ta_yyyzz_yyy_1[i] * pc_z[i];

        ta_yyyzzz_yyz_0[i] = 2.0 * ta_yzzz_yyz_0[i] * fe_0 - 2.0 * ta_yzzz_yyz_1[i] * fe_0 + 2.0 * ta_yyzzz_yz_0[i] * fe_0 -
                             2.0 * ta_yyzzz_yz_1[i] * fe_0 + ta_yyzzz_yyz_0[i] * pa_y[i] - ta_yyzzz_yyz_1[i] * pc_y[i];

        ta_yyyzzz_yzz_0[i] = 2.0 * ta_yzzz_yzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzz_1[i] * fe_0 + ta_yyzzz_zz_0[i] * fe_0 - ta_yyzzz_zz_1[i] * fe_0 +
                             ta_yyzzz_yzz_0[i] * pa_y[i] - ta_yyzzz_yzz_1[i] * pc_y[i];

        ta_yyyzzz_zzz_0[i] =
            2.0 * ta_yzzz_zzz_0[i] * fe_0 - 2.0 * ta_yzzz_zzz_1[i] * fe_0 + ta_yyzzz_zzz_0[i] * pa_y[i] - ta_yyzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 250-260 components of targeted buffer : IF

    auto ta_yyzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 250);

    auto ta_yyzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 251);

    auto ta_yyzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 252);

    auto ta_yyzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 253);

    auto ta_yyzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 254);

    auto ta_yyzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 255);

    auto ta_yyzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 256);

    auto ta_yyzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 257);

    auto ta_yyzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 258);

    auto ta_yyzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 259);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyzz_xxy_0,   \
                             ta_yyzz_xxy_1,   \
                             ta_yyzz_xyy_0,   \
                             ta_yyzz_xyy_1,   \
                             ta_yyzz_yyy_0,   \
                             ta_yyzz_yyy_1,   \
                             ta_yyzzz_xxy_0,  \
                             ta_yyzzz_xxy_1,  \
                             ta_yyzzz_xyy_0,  \
                             ta_yyzzz_xyy_1,  \
                             ta_yyzzz_yyy_0,  \
                             ta_yyzzz_yyy_1,  \
                             ta_yyzzzz_xxx_0, \
                             ta_yyzzzz_xxy_0, \
                             ta_yyzzzz_xxz_0, \
                             ta_yyzzzz_xyy_0, \
                             ta_yyzzzz_xyz_0, \
                             ta_yyzzzz_xzz_0, \
                             ta_yyzzzz_yyy_0, \
                             ta_yyzzzz_yyz_0, \
                             ta_yyzzzz_yzz_0, \
                             ta_yyzzzz_zzz_0, \
                             ta_yzzzz_xxx_0,  \
                             ta_yzzzz_xxx_1,  \
                             ta_yzzzz_xxz_0,  \
                             ta_yzzzz_xxz_1,  \
                             ta_yzzzz_xyz_0,  \
                             ta_yzzzz_xyz_1,  \
                             ta_yzzzz_xz_0,   \
                             ta_yzzzz_xz_1,   \
                             ta_yzzzz_xzz_0,  \
                             ta_yzzzz_xzz_1,  \
                             ta_yzzzz_yyz_0,  \
                             ta_yzzzz_yyz_1,  \
                             ta_yzzzz_yz_0,   \
                             ta_yzzzz_yz_1,   \
                             ta_yzzzz_yzz_0,  \
                             ta_yzzzz_yzz_1,  \
                             ta_yzzzz_zz_0,   \
                             ta_yzzzz_zz_1,   \
                             ta_yzzzz_zzz_0,  \
                             ta_yzzzz_zzz_1,  \
                             ta_zzzz_xxx_0,   \
                             ta_zzzz_xxx_1,   \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_xxx_0[i] = ta_zzzz_xxx_0[i] * fe_0 - ta_zzzz_xxx_1[i] * fe_0 + ta_yzzzz_xxx_0[i] * pa_y[i] - ta_yzzzz_xxx_1[i] * pc_y[i];

        ta_yyzzzz_xxy_0[i] =
            3.0 * ta_yyzz_xxy_0[i] * fe_0 - 3.0 * ta_yyzz_xxy_1[i] * fe_0 + ta_yyzzz_xxy_0[i] * pa_z[i] - ta_yyzzz_xxy_1[i] * pc_z[i];

        ta_yyzzzz_xxz_0[i] = ta_zzzz_xxz_0[i] * fe_0 - ta_zzzz_xxz_1[i] * fe_0 + ta_yzzzz_xxz_0[i] * pa_y[i] - ta_yzzzz_xxz_1[i] * pc_y[i];

        ta_yyzzzz_xyy_0[i] =
            3.0 * ta_yyzz_xyy_0[i] * fe_0 - 3.0 * ta_yyzz_xyy_1[i] * fe_0 + ta_yyzzz_xyy_0[i] * pa_z[i] - ta_yyzzz_xyy_1[i] * pc_z[i];

        ta_yyzzzz_xyz_0[i] = ta_zzzz_xyz_0[i] * fe_0 - ta_zzzz_xyz_1[i] * fe_0 + ta_yzzzz_xz_0[i] * fe_0 - ta_yzzzz_xz_1[i] * fe_0 +
                             ta_yzzzz_xyz_0[i] * pa_y[i] - ta_yzzzz_xyz_1[i] * pc_y[i];

        ta_yyzzzz_xzz_0[i] = ta_zzzz_xzz_0[i] * fe_0 - ta_zzzz_xzz_1[i] * fe_0 + ta_yzzzz_xzz_0[i] * pa_y[i] - ta_yzzzz_xzz_1[i] * pc_y[i];

        ta_yyzzzz_yyy_0[i] =
            3.0 * ta_yyzz_yyy_0[i] * fe_0 - 3.0 * ta_yyzz_yyy_1[i] * fe_0 + ta_yyzzz_yyy_0[i] * pa_z[i] - ta_yyzzz_yyy_1[i] * pc_z[i];

        ta_yyzzzz_yyz_0[i] = ta_zzzz_yyz_0[i] * fe_0 - ta_zzzz_yyz_1[i] * fe_0 + 2.0 * ta_yzzzz_yz_0[i] * fe_0 - 2.0 * ta_yzzzz_yz_1[i] * fe_0 +
                             ta_yzzzz_yyz_0[i] * pa_y[i] - ta_yzzzz_yyz_1[i] * pc_y[i];

        ta_yyzzzz_yzz_0[i] = ta_zzzz_yzz_0[i] * fe_0 - ta_zzzz_yzz_1[i] * fe_0 + ta_yzzzz_zz_0[i] * fe_0 - ta_yzzzz_zz_1[i] * fe_0 +
                             ta_yzzzz_yzz_0[i] * pa_y[i] - ta_yzzzz_yzz_1[i] * pc_y[i];

        ta_yyzzzz_zzz_0[i] = ta_zzzz_zzz_0[i] * fe_0 - ta_zzzz_zzz_1[i] * fe_0 + ta_yzzzz_zzz_0[i] * pa_y[i] - ta_yzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 260-270 components of targeted buffer : IF

    auto ta_yzzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 260);

    auto ta_yzzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 261);

    auto ta_yzzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 262);

    auto ta_yzzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 263);

    auto ta_yzzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 264);

    auto ta_yzzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 265);

    auto ta_yzzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 266);

    auto ta_yzzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 267);

    auto ta_yzzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 268);

    auto ta_yzzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 269);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta_yzzzzz_xxx_0, \
                             ta_yzzzzz_xxy_0, \
                             ta_yzzzzz_xxz_0, \
                             ta_yzzzzz_xyy_0, \
                             ta_yzzzzz_xyz_0, \
                             ta_yzzzzz_xzz_0, \
                             ta_yzzzzz_yyy_0, \
                             ta_yzzzzz_yyz_0, \
                             ta_yzzzzz_yzz_0, \
                             ta_yzzzzz_zzz_0, \
                             ta_zzzzz_xx_0,   \
                             ta_zzzzz_xx_1,   \
                             ta_zzzzz_xxx_0,  \
                             ta_zzzzz_xxx_1,  \
                             ta_zzzzz_xxy_0,  \
                             ta_zzzzz_xxy_1,  \
                             ta_zzzzz_xxz_0,  \
                             ta_zzzzz_xxz_1,  \
                             ta_zzzzz_xy_0,   \
                             ta_zzzzz_xy_1,   \
                             ta_zzzzz_xyy_0,  \
                             ta_zzzzz_xyy_1,  \
                             ta_zzzzz_xyz_0,  \
                             ta_zzzzz_xyz_1,  \
                             ta_zzzzz_xz_0,   \
                             ta_zzzzz_xz_1,   \
                             ta_zzzzz_xzz_0,  \
                             ta_zzzzz_xzz_1,  \
                             ta_zzzzz_yy_0,   \
                             ta_zzzzz_yy_1,   \
                             ta_zzzzz_yyy_0,  \
                             ta_zzzzz_yyy_1,  \
                             ta_zzzzz_yyz_0,  \
                             ta_zzzzz_yyz_1,  \
                             ta_zzzzz_yz_0,   \
                             ta_zzzzz_yz_1,   \
                             ta_zzzzz_yzz_0,  \
                             ta_zzzzz_yzz_1,  \
                             ta_zzzzz_zz_0,   \
                             ta_zzzzz_zz_1,   \
                             ta_zzzzz_zzz_0,  \
                             ta_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_xxx_0[i] = ta_zzzzz_xxx_0[i] * pa_y[i] - ta_zzzzz_xxx_1[i] * pc_y[i];

        ta_yzzzzz_xxy_0[i] = ta_zzzzz_xx_0[i] * fe_0 - ta_zzzzz_xx_1[i] * fe_0 + ta_zzzzz_xxy_0[i] * pa_y[i] - ta_zzzzz_xxy_1[i] * pc_y[i];

        ta_yzzzzz_xxz_0[i] = ta_zzzzz_xxz_0[i] * pa_y[i] - ta_zzzzz_xxz_1[i] * pc_y[i];

        ta_yzzzzz_xyy_0[i] =
            2.0 * ta_zzzzz_xy_0[i] * fe_0 - 2.0 * ta_zzzzz_xy_1[i] * fe_0 + ta_zzzzz_xyy_0[i] * pa_y[i] - ta_zzzzz_xyy_1[i] * pc_y[i];

        ta_yzzzzz_xyz_0[i] = ta_zzzzz_xz_0[i] * fe_0 - ta_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xyz_0[i] * pa_y[i] - ta_zzzzz_xyz_1[i] * pc_y[i];

        ta_yzzzzz_xzz_0[i] = ta_zzzzz_xzz_0[i] * pa_y[i] - ta_zzzzz_xzz_1[i] * pc_y[i];

        ta_yzzzzz_yyy_0[i] =
            3.0 * ta_zzzzz_yy_0[i] * fe_0 - 3.0 * ta_zzzzz_yy_1[i] * fe_0 + ta_zzzzz_yyy_0[i] * pa_y[i] - ta_zzzzz_yyy_1[i] * pc_y[i];

        ta_yzzzzz_yyz_0[i] =
            2.0 * ta_zzzzz_yz_0[i] * fe_0 - 2.0 * ta_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_yyz_0[i] * pa_y[i] - ta_zzzzz_yyz_1[i] * pc_y[i];

        ta_yzzzzz_yzz_0[i] = ta_zzzzz_zz_0[i] * fe_0 - ta_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_yzz_0[i] * pa_y[i] - ta_zzzzz_yzz_1[i] * pc_y[i];

        ta_yzzzzz_zzz_0[i] = ta_zzzzz_zzz_0[i] * pa_y[i] - ta_zzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 270-280 components of targeted buffer : IF

    auto ta_zzzzzz_xxx_0 = pbuffer.data(idx_npot_0_if + 270);

    auto ta_zzzzzz_xxy_0 = pbuffer.data(idx_npot_0_if + 271);

    auto ta_zzzzzz_xxz_0 = pbuffer.data(idx_npot_0_if + 272);

    auto ta_zzzzzz_xyy_0 = pbuffer.data(idx_npot_0_if + 273);

    auto ta_zzzzzz_xyz_0 = pbuffer.data(idx_npot_0_if + 274);

    auto ta_zzzzzz_xzz_0 = pbuffer.data(idx_npot_0_if + 275);

    auto ta_zzzzzz_yyy_0 = pbuffer.data(idx_npot_0_if + 276);

    auto ta_zzzzzz_yyz_0 = pbuffer.data(idx_npot_0_if + 277);

    auto ta_zzzzzz_yzz_0 = pbuffer.data(idx_npot_0_if + 278);

    auto ta_zzzzzz_zzz_0 = pbuffer.data(idx_npot_0_if + 279);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta_zzzz_xxx_0,   \
                             ta_zzzz_xxx_1,   \
                             ta_zzzz_xxy_0,   \
                             ta_zzzz_xxy_1,   \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xyy_0,   \
                             ta_zzzz_xyy_1,   \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_yyy_0,   \
                             ta_zzzz_yyy_1,   \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             ta_zzzzz_xx_0,   \
                             ta_zzzzz_xx_1,   \
                             ta_zzzzz_xxx_0,  \
                             ta_zzzzz_xxx_1,  \
                             ta_zzzzz_xxy_0,  \
                             ta_zzzzz_xxy_1,  \
                             ta_zzzzz_xxz_0,  \
                             ta_zzzzz_xxz_1,  \
                             ta_zzzzz_xy_0,   \
                             ta_zzzzz_xy_1,   \
                             ta_zzzzz_xyy_0,  \
                             ta_zzzzz_xyy_1,  \
                             ta_zzzzz_xyz_0,  \
                             ta_zzzzz_xyz_1,  \
                             ta_zzzzz_xz_0,   \
                             ta_zzzzz_xz_1,   \
                             ta_zzzzz_xzz_0,  \
                             ta_zzzzz_xzz_1,  \
                             ta_zzzzz_yy_0,   \
                             ta_zzzzz_yy_1,   \
                             ta_zzzzz_yyy_0,  \
                             ta_zzzzz_yyy_1,  \
                             ta_zzzzz_yyz_0,  \
                             ta_zzzzz_yyz_1,  \
                             ta_zzzzz_yz_0,   \
                             ta_zzzzz_yz_1,   \
                             ta_zzzzz_yzz_0,  \
                             ta_zzzzz_yzz_1,  \
                             ta_zzzzz_zz_0,   \
                             ta_zzzzz_zz_1,   \
                             ta_zzzzz_zzz_0,  \
                             ta_zzzzz_zzz_1,  \
                             ta_zzzzzz_xxx_0, \
                             ta_zzzzzz_xxy_0, \
                             ta_zzzzzz_xxz_0, \
                             ta_zzzzzz_xyy_0, \
                             ta_zzzzzz_xyz_0, \
                             ta_zzzzzz_xzz_0, \
                             ta_zzzzzz_yyy_0, \
                             ta_zzzzzz_yyz_0, \
                             ta_zzzzzz_yzz_0, \
                             ta_zzzzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_xxx_0[i] =
            5.0 * ta_zzzz_xxx_0[i] * fe_0 - 5.0 * ta_zzzz_xxx_1[i] * fe_0 + ta_zzzzz_xxx_0[i] * pa_z[i] - ta_zzzzz_xxx_1[i] * pc_z[i];

        ta_zzzzzz_xxy_0[i] =
            5.0 * ta_zzzz_xxy_0[i] * fe_0 - 5.0 * ta_zzzz_xxy_1[i] * fe_0 + ta_zzzzz_xxy_0[i] * pa_z[i] - ta_zzzzz_xxy_1[i] * pc_z[i];

        ta_zzzzzz_xxz_0[i] = 5.0 * ta_zzzz_xxz_0[i] * fe_0 - 5.0 * ta_zzzz_xxz_1[i] * fe_0 + ta_zzzzz_xx_0[i] * fe_0 - ta_zzzzz_xx_1[i] * fe_0 +
                             ta_zzzzz_xxz_0[i] * pa_z[i] - ta_zzzzz_xxz_1[i] * pc_z[i];

        ta_zzzzzz_xyy_0[i] =
            5.0 * ta_zzzz_xyy_0[i] * fe_0 - 5.0 * ta_zzzz_xyy_1[i] * fe_0 + ta_zzzzz_xyy_0[i] * pa_z[i] - ta_zzzzz_xyy_1[i] * pc_z[i];

        ta_zzzzzz_xyz_0[i] = 5.0 * ta_zzzz_xyz_0[i] * fe_0 - 5.0 * ta_zzzz_xyz_1[i] * fe_0 + ta_zzzzz_xy_0[i] * fe_0 - ta_zzzzz_xy_1[i] * fe_0 +
                             ta_zzzzz_xyz_0[i] * pa_z[i] - ta_zzzzz_xyz_1[i] * pc_z[i];

        ta_zzzzzz_xzz_0[i] = 5.0 * ta_zzzz_xzz_0[i] * fe_0 - 5.0 * ta_zzzz_xzz_1[i] * fe_0 + 2.0 * ta_zzzzz_xz_0[i] * fe_0 -
                             2.0 * ta_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xzz_0[i] * pa_z[i] - ta_zzzzz_xzz_1[i] * pc_z[i];

        ta_zzzzzz_yyy_0[i] =
            5.0 * ta_zzzz_yyy_0[i] * fe_0 - 5.0 * ta_zzzz_yyy_1[i] * fe_0 + ta_zzzzz_yyy_0[i] * pa_z[i] - ta_zzzzz_yyy_1[i] * pc_z[i];

        ta_zzzzzz_yyz_0[i] = 5.0 * ta_zzzz_yyz_0[i] * fe_0 - 5.0 * ta_zzzz_yyz_1[i] * fe_0 + ta_zzzzz_yy_0[i] * fe_0 - ta_zzzzz_yy_1[i] * fe_0 +
                             ta_zzzzz_yyz_0[i] * pa_z[i] - ta_zzzzz_yyz_1[i] * pc_z[i];

        ta_zzzzzz_yzz_0[i] = 5.0 * ta_zzzz_yzz_0[i] * fe_0 - 5.0 * ta_zzzz_yzz_1[i] * fe_0 + 2.0 * ta_zzzzz_yz_0[i] * fe_0 -
                             2.0 * ta_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_yzz_0[i] * pa_z[i] - ta_zzzzz_yzz_1[i] * pc_z[i];

        ta_zzzzzz_zzz_0[i] = 5.0 * ta_zzzz_zzz_0[i] * fe_0 - 5.0 * ta_zzzz_zzz_1[i] * fe_0 + 3.0 * ta_zzzzz_zz_0[i] * fe_0 -
                             3.0 * ta_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_zzz_0[i] * pa_z[i] - ta_zzzzz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
