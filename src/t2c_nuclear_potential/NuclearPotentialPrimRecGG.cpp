#include "NuclearPotentialPrimRecGG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_gg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_gg,
                               const size_t              idx_npot_0_dg,
                               const size_t              idx_npot_1_dg,
                               const size_t              idx_npot_0_ff,
                               const size_t              idx_npot_1_ff,
                               const size_t              idx_npot_0_fg,
                               const size_t              idx_npot_1_fg,
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

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_0 = pbuffer.data(idx_npot_0_dg);

    auto ta_xx_xxxy_0 = pbuffer.data(idx_npot_0_dg + 1);

    auto ta_xx_xxxz_0 = pbuffer.data(idx_npot_0_dg + 2);

    auto ta_xx_xxyy_0 = pbuffer.data(idx_npot_0_dg + 3);

    auto ta_xx_xxyz_0 = pbuffer.data(idx_npot_0_dg + 4);

    auto ta_xx_xxzz_0 = pbuffer.data(idx_npot_0_dg + 5);

    auto ta_xx_xyyy_0 = pbuffer.data(idx_npot_0_dg + 6);

    auto ta_xx_xyyz_0 = pbuffer.data(idx_npot_0_dg + 7);

    auto ta_xx_xyzz_0 = pbuffer.data(idx_npot_0_dg + 8);

    auto ta_xx_xzzz_0 = pbuffer.data(idx_npot_0_dg + 9);

    auto ta_xx_yyyy_0 = pbuffer.data(idx_npot_0_dg + 10);

    auto ta_xx_yyyz_0 = pbuffer.data(idx_npot_0_dg + 11);

    auto ta_xx_yyzz_0 = pbuffer.data(idx_npot_0_dg + 12);

    auto ta_xx_yzzz_0 = pbuffer.data(idx_npot_0_dg + 13);

    auto ta_xx_zzzz_0 = pbuffer.data(idx_npot_0_dg + 14);

    auto ta_xy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 25);

    auto ta_xy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 26);

    auto ta_xy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 27);

    auto ta_xy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 28);

    auto ta_xz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 41);

    auto ta_xz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 42);

    auto ta_xz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 43);

    auto ta_xz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 44);

    auto ta_yy_xxxx_0 = pbuffer.data(idx_npot_0_dg + 45);

    auto ta_yy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 46);

    auto ta_yy_xxxz_0 = pbuffer.data(idx_npot_0_dg + 47);

    auto ta_yy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 48);

    auto ta_yy_xxyz_0 = pbuffer.data(idx_npot_0_dg + 49);

    auto ta_yy_xxzz_0 = pbuffer.data(idx_npot_0_dg + 50);

    auto ta_yy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 51);

    auto ta_yy_xyyz_0 = pbuffer.data(idx_npot_0_dg + 52);

    auto ta_yy_xyzz_0 = pbuffer.data(idx_npot_0_dg + 53);

    auto ta_yy_xzzz_0 = pbuffer.data(idx_npot_0_dg + 54);

    auto ta_yy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 55);

    auto ta_yy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 56);

    auto ta_yy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 57);

    auto ta_yy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 58);

    auto ta_yy_zzzz_0 = pbuffer.data(idx_npot_0_dg + 59);

    auto ta_yz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 62);

    auto ta_yz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 65);

    auto ta_yz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 69);

    auto ta_yz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 71);

    auto ta_yz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 72);

    auto ta_yz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 73);

    auto ta_yz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 74);

    auto ta_zz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 75);

    auto ta_zz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 76);

    auto ta_zz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 77);

    auto ta_zz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 78);

    auto ta_zz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 79);

    auto ta_zz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 80);

    auto ta_zz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 81);

    auto ta_zz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 82);

    auto ta_zz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 83);

    auto ta_zz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 84);

    auto ta_zz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 85);

    auto ta_zz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 86);

    auto ta_zz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 87);

    auto ta_zz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 88);

    auto ta_zz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_1 = pbuffer.data(idx_npot_1_dg);

    auto ta_xx_xxxy_1 = pbuffer.data(idx_npot_1_dg + 1);

    auto ta_xx_xxxz_1 = pbuffer.data(idx_npot_1_dg + 2);

    auto ta_xx_xxyy_1 = pbuffer.data(idx_npot_1_dg + 3);

    auto ta_xx_xxyz_1 = pbuffer.data(idx_npot_1_dg + 4);

    auto ta_xx_xxzz_1 = pbuffer.data(idx_npot_1_dg + 5);

    auto ta_xx_xyyy_1 = pbuffer.data(idx_npot_1_dg + 6);

    auto ta_xx_xyyz_1 = pbuffer.data(idx_npot_1_dg + 7);

    auto ta_xx_xyzz_1 = pbuffer.data(idx_npot_1_dg + 8);

    auto ta_xx_xzzz_1 = pbuffer.data(idx_npot_1_dg + 9);

    auto ta_xx_yyyy_1 = pbuffer.data(idx_npot_1_dg + 10);

    auto ta_xx_yyyz_1 = pbuffer.data(idx_npot_1_dg + 11);

    auto ta_xx_yyzz_1 = pbuffer.data(idx_npot_1_dg + 12);

    auto ta_xx_yzzz_1 = pbuffer.data(idx_npot_1_dg + 13);

    auto ta_xx_zzzz_1 = pbuffer.data(idx_npot_1_dg + 14);

    auto ta_xy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 25);

    auto ta_xy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 26);

    auto ta_xy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 27);

    auto ta_xy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 28);

    auto ta_xz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 41);

    auto ta_xz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 42);

    auto ta_xz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 43);

    auto ta_xz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 44);

    auto ta_yy_xxxx_1 = pbuffer.data(idx_npot_1_dg + 45);

    auto ta_yy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 46);

    auto ta_yy_xxxz_1 = pbuffer.data(idx_npot_1_dg + 47);

    auto ta_yy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 48);

    auto ta_yy_xxyz_1 = pbuffer.data(idx_npot_1_dg + 49);

    auto ta_yy_xxzz_1 = pbuffer.data(idx_npot_1_dg + 50);

    auto ta_yy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 51);

    auto ta_yy_xyyz_1 = pbuffer.data(idx_npot_1_dg + 52);

    auto ta_yy_xyzz_1 = pbuffer.data(idx_npot_1_dg + 53);

    auto ta_yy_xzzz_1 = pbuffer.data(idx_npot_1_dg + 54);

    auto ta_yy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 55);

    auto ta_yy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 56);

    auto ta_yy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 57);

    auto ta_yy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 58);

    auto ta_yy_zzzz_1 = pbuffer.data(idx_npot_1_dg + 59);

    auto ta_yz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 62);

    auto ta_yz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 65);

    auto ta_yz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 69);

    auto ta_yz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 71);

    auto ta_yz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 72);

    auto ta_yz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 73);

    auto ta_yz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 74);

    auto ta_zz_xxxx_1 = pbuffer.data(idx_npot_1_dg + 75);

    auto ta_zz_xxxy_1 = pbuffer.data(idx_npot_1_dg + 76);

    auto ta_zz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 77);

    auto ta_zz_xxyy_1 = pbuffer.data(idx_npot_1_dg + 78);

    auto ta_zz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 79);

    auto ta_zz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 80);

    auto ta_zz_xyyy_1 = pbuffer.data(idx_npot_1_dg + 81);

    auto ta_zz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 82);

    auto ta_zz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 83);

    auto ta_zz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 84);

    auto ta_zz_yyyy_1 = pbuffer.data(idx_npot_1_dg + 85);

    auto ta_zz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 86);

    auto ta_zz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 87);

    auto ta_zz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 88);

    auto ta_zz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 89);

    // Set up components of auxiliary buffer : FF

    auto ta_xxx_xxx_0 = pbuffer.data(idx_npot_0_ff);

    auto ta_xxx_xxy_0 = pbuffer.data(idx_npot_0_ff + 1);

    auto ta_xxx_xxz_0 = pbuffer.data(idx_npot_0_ff + 2);

    auto ta_xxx_xyy_0 = pbuffer.data(idx_npot_0_ff + 3);

    auto ta_xxx_xyz_0 = pbuffer.data(idx_npot_0_ff + 4);

    auto ta_xxx_xzz_0 = pbuffer.data(idx_npot_0_ff + 5);

    auto ta_xxx_yyy_0 = pbuffer.data(idx_npot_0_ff + 6);

    auto ta_xxx_yyz_0 = pbuffer.data(idx_npot_0_ff + 7);

    auto ta_xxx_yzz_0 = pbuffer.data(idx_npot_0_ff + 8);

    auto ta_xxx_zzz_0 = pbuffer.data(idx_npot_0_ff + 9);

    auto ta_xxz_xxz_0 = pbuffer.data(idx_npot_0_ff + 22);

    auto ta_xxz_xyz_0 = pbuffer.data(idx_npot_0_ff + 24);

    auto ta_xxz_xzz_0 = pbuffer.data(idx_npot_0_ff + 25);

    auto ta_xyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 31);

    auto ta_xyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 33);

    auto ta_xyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 34);

    auto ta_xyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 36);

    auto ta_xyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 37);

    auto ta_xyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 38);

    auto ta_xzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 52);

    auto ta_xzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 54);

    auto ta_xzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 55);

    auto ta_xzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 57);

    auto ta_xzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 58);

    auto ta_xzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 59);

    auto ta_yyy_xxx_0 = pbuffer.data(idx_npot_0_ff + 60);

    auto ta_yyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 61);

    auto ta_yyy_xxz_0 = pbuffer.data(idx_npot_0_ff + 62);

    auto ta_yyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 63);

    auto ta_yyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 64);

    auto ta_yyy_xzz_0 = pbuffer.data(idx_npot_0_ff + 65);

    auto ta_yyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 66);

    auto ta_yyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 67);

    auto ta_yyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 68);

    auto ta_yyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 69);

    auto ta_yyz_xxz_0 = pbuffer.data(idx_npot_0_ff + 72);

    auto ta_yyz_xyz_0 = pbuffer.data(idx_npot_0_ff + 74);

    auto ta_yyz_xzz_0 = pbuffer.data(idx_npot_0_ff + 75);

    auto ta_yyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 77);

    auto ta_yyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 78);

    auto ta_yyz_zzz_0 = pbuffer.data(idx_npot_0_ff + 79);

    auto ta_yzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 81);

    auto ta_yzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 82);

    auto ta_yzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 83);

    auto ta_yzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 84);

    auto ta_yzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 85);

    auto ta_yzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 86);

    auto ta_yzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 87);

    auto ta_yzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 88);

    auto ta_yzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 89);

    auto ta_zzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 90);

    auto ta_zzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 91);

    auto ta_zzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 92);

    auto ta_zzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 93);

    auto ta_zzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 94);

    auto ta_zzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 95);

    auto ta_zzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 96);

    auto ta_zzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 97);

    auto ta_zzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 98);

    auto ta_zzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 99);

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

    auto ta_xxz_xxz_1 = pbuffer.data(idx_npot_1_ff + 22);

    auto ta_xxz_xyz_1 = pbuffer.data(idx_npot_1_ff + 24);

    auto ta_xxz_xzz_1 = pbuffer.data(idx_npot_1_ff + 25);

    auto ta_xyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 31);

    auto ta_xyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 33);

    auto ta_xyy_xyz_1 = pbuffer.data(idx_npot_1_ff + 34);

    auto ta_xyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 36);

    auto ta_xyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 37);

    auto ta_xyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 38);

    auto ta_xzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 52);

    auto ta_xzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 54);

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

    auto ta_yyz_xxz_1 = pbuffer.data(idx_npot_1_ff + 72);

    auto ta_yyz_xyz_1 = pbuffer.data(idx_npot_1_ff + 74);

    auto ta_yyz_xzz_1 = pbuffer.data(idx_npot_1_ff + 75);

    auto ta_yyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 77);

    auto ta_yyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 78);

    auto ta_yyz_zzz_1 = pbuffer.data(idx_npot_1_ff + 79);

    auto ta_yzz_xxy_1 = pbuffer.data(idx_npot_1_ff + 81);

    auto ta_yzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 82);

    auto ta_yzz_xyy_1 = pbuffer.data(idx_npot_1_ff + 83);

    auto ta_yzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 84);

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

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_0 = pbuffer.data(idx_npot_0_fg);

    auto ta_xxx_xxxy_0 = pbuffer.data(idx_npot_0_fg + 1);

    auto ta_xxx_xxxz_0 = pbuffer.data(idx_npot_0_fg + 2);

    auto ta_xxx_xxyy_0 = pbuffer.data(idx_npot_0_fg + 3);

    auto ta_xxx_xxyz_0 = pbuffer.data(idx_npot_0_fg + 4);

    auto ta_xxx_xxzz_0 = pbuffer.data(idx_npot_0_fg + 5);

    auto ta_xxx_xyyy_0 = pbuffer.data(idx_npot_0_fg + 6);

    auto ta_xxx_xyyz_0 = pbuffer.data(idx_npot_0_fg + 7);

    auto ta_xxx_xyzz_0 = pbuffer.data(idx_npot_0_fg + 8);

    auto ta_xxx_xzzz_0 = pbuffer.data(idx_npot_0_fg + 9);

    auto ta_xxx_yyyy_0 = pbuffer.data(idx_npot_0_fg + 10);

    auto ta_xxx_yyyz_0 = pbuffer.data(idx_npot_0_fg + 11);

    auto ta_xxx_yyzz_0 = pbuffer.data(idx_npot_0_fg + 12);

    auto ta_xxx_yzzz_0 = pbuffer.data(idx_npot_0_fg + 13);

    auto ta_xxx_zzzz_0 = pbuffer.data(idx_npot_0_fg + 14);

    auto ta_xxy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 15);

    auto ta_xxy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 16);

    auto ta_xxy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 17);

    auto ta_xxy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 18);

    auto ta_xxy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 20);

    auto ta_xxy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 21);

    auto ta_xxy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 24);

    auto ta_xxy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 25);

    auto ta_xxy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 26);

    auto ta_xxy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 27);

    auto ta_xxy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 28);

    auto ta_xxz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 30);

    auto ta_xxz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 31);

    auto ta_xxz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 32);

    auto ta_xxz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 33);

    auto ta_xxz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 34);

    auto ta_xxz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 35);

    auto ta_xxz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 36);

    auto ta_xxz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 37);

    auto ta_xxz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 38);

    auto ta_xxz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 39);

    auto ta_xxz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 41);

    auto ta_xxz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 42);

    auto ta_xxz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 43);

    auto ta_xxz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 44);

    auto ta_xyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 45);

    auto ta_xyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 46);

    auto ta_xyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 48);

    auto ta_xyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 49);

    auto ta_xyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 51);

    auto ta_xyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 52);

    auto ta_xyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 53);

    auto ta_xyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 55);

    auto ta_xyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 56);

    auto ta_xyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 57);

    auto ta_xyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 58);

    auto ta_xyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 59);

    auto ta_xyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 71);

    auto ta_xyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 72);

    auto ta_xyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 73);

    auto ta_xzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 75);

    auto ta_xzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 77);

    auto ta_xzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 79);

    auto ta_xzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 80);

    auto ta_xzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 82);

    auto ta_xzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 83);

    auto ta_xzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 84);

    auto ta_xzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 85);

    auto ta_xzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 86);

    auto ta_xzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 87);

    auto ta_xzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 88);

    auto ta_xzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 89);

    auto ta_yyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 90);

    auto ta_yyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 91);

    auto ta_yyy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 92);

    auto ta_yyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 93);

    auto ta_yyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 94);

    auto ta_yyy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 95);

    auto ta_yyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 96);

    auto ta_yyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 97);

    auto ta_yyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 98);

    auto ta_yyy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 99);

    auto ta_yyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 100);

    auto ta_yyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 101);

    auto ta_yyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 102);

    auto ta_yyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 103);

    auto ta_yyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 104);

    auto ta_yyz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 106);

    auto ta_yyz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 107);

    auto ta_yyz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 108);

    auto ta_yyz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 109);

    auto ta_yyz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 110);

    auto ta_yyz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 111);

    auto ta_yyz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 112);

    auto ta_yyz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 113);

    auto ta_yyz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 114);

    auto ta_yyz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 115);

    auto ta_yyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 116);

    auto ta_yyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 117);

    auto ta_yyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 118);

    auto ta_yyz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 119);

    auto ta_yzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 120);

    auto ta_yzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 121);

    auto ta_yzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 122);

    auto ta_yzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 123);

    auto ta_yzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 124);

    auto ta_yzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 125);

    auto ta_yzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 126);

    auto ta_yzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 127);

    auto ta_yzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 128);

    auto ta_yzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 129);

    auto ta_yzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 130);

    auto ta_yzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 131);

    auto ta_yzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 132);

    auto ta_yzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 133);

    auto ta_yzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 134);

    auto ta_zzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 135);

    auto ta_zzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 136);

    auto ta_zzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 137);

    auto ta_zzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 138);

    auto ta_zzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 139);

    auto ta_zzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 140);

    auto ta_zzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 141);

    auto ta_zzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 142);

    auto ta_zzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 143);

    auto ta_zzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 144);

    auto ta_zzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 145);

    auto ta_zzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 146);

    auto ta_zzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 147);

    auto ta_zzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 148);

    auto ta_zzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 149);

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_1 = pbuffer.data(idx_npot_1_fg);

    auto ta_xxx_xxxy_1 = pbuffer.data(idx_npot_1_fg + 1);

    auto ta_xxx_xxxz_1 = pbuffer.data(idx_npot_1_fg + 2);

    auto ta_xxx_xxyy_1 = pbuffer.data(idx_npot_1_fg + 3);

    auto ta_xxx_xxyz_1 = pbuffer.data(idx_npot_1_fg + 4);

    auto ta_xxx_xxzz_1 = pbuffer.data(idx_npot_1_fg + 5);

    auto ta_xxx_xyyy_1 = pbuffer.data(idx_npot_1_fg + 6);

    auto ta_xxx_xyyz_1 = pbuffer.data(idx_npot_1_fg + 7);

    auto ta_xxx_xyzz_1 = pbuffer.data(idx_npot_1_fg + 8);

    auto ta_xxx_xzzz_1 = pbuffer.data(idx_npot_1_fg + 9);

    auto ta_xxx_yyyy_1 = pbuffer.data(idx_npot_1_fg + 10);

    auto ta_xxx_yyyz_1 = pbuffer.data(idx_npot_1_fg + 11);

    auto ta_xxx_yyzz_1 = pbuffer.data(idx_npot_1_fg + 12);

    auto ta_xxx_yzzz_1 = pbuffer.data(idx_npot_1_fg + 13);

    auto ta_xxx_zzzz_1 = pbuffer.data(idx_npot_1_fg + 14);

    auto ta_xxy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 15);

    auto ta_xxy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 16);

    auto ta_xxy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 17);

    auto ta_xxy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 18);

    auto ta_xxy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 20);

    auto ta_xxy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 21);

    auto ta_xxy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 24);

    auto ta_xxy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 25);

    auto ta_xxy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 26);

    auto ta_xxy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 27);

    auto ta_xxy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 28);

    auto ta_xxz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 30);

    auto ta_xxz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 31);

    auto ta_xxz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 32);

    auto ta_xxz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 33);

    auto ta_xxz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 34);

    auto ta_xxz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 35);

    auto ta_xxz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 36);

    auto ta_xxz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 37);

    auto ta_xxz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 38);

    auto ta_xxz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 39);

    auto ta_xxz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 41);

    auto ta_xxz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 42);

    auto ta_xxz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 43);

    auto ta_xxz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 44);

    auto ta_xyy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 45);

    auto ta_xyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 46);

    auto ta_xyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 48);

    auto ta_xyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 49);

    auto ta_xyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 51);

    auto ta_xyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 52);

    auto ta_xyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 53);

    auto ta_xyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 55);

    auto ta_xyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 56);

    auto ta_xyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 57);

    auto ta_xyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 58);

    auto ta_xyy_zzzz_1 = pbuffer.data(idx_npot_1_fg + 59);

    auto ta_xyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 71);

    auto ta_xyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 72);

    auto ta_xyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 73);

    auto ta_xzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 75);

    auto ta_xzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 77);

    auto ta_xzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 79);

    auto ta_xzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 80);

    auto ta_xzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 82);

    auto ta_xzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 83);

    auto ta_xzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 84);

    auto ta_xzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 85);

    auto ta_xzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 86);

    auto ta_xzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 87);

    auto ta_xzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 88);

    auto ta_xzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 89);

    auto ta_yyy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 90);

    auto ta_yyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 91);

    auto ta_yyy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 92);

    auto ta_yyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 93);

    auto ta_yyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 94);

    auto ta_yyy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 95);

    auto ta_yyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 96);

    auto ta_yyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 97);

    auto ta_yyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 98);

    auto ta_yyy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 99);

    auto ta_yyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 100);

    auto ta_yyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 101);

    auto ta_yyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 102);

    auto ta_yyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 103);

    auto ta_yyy_zzzz_1 = pbuffer.data(idx_npot_1_fg + 104);

    auto ta_yyz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 106);

    auto ta_yyz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 107);

    auto ta_yyz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 108);

    auto ta_yyz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 109);

    auto ta_yyz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 110);

    auto ta_yyz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 111);

    auto ta_yyz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 112);

    auto ta_yyz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 113);

    auto ta_yyz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 114);

    auto ta_yyz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 115);

    auto ta_yyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 116);

    auto ta_yyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 117);

    auto ta_yyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 118);

    auto ta_yyz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 119);

    auto ta_yzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 120);

    auto ta_yzz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 121);

    auto ta_yzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 122);

    auto ta_yzz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 123);

    auto ta_yzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 124);

    auto ta_yzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 125);

    auto ta_yzz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 126);

    auto ta_yzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 127);

    auto ta_yzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 128);

    auto ta_yzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 129);

    auto ta_yzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 130);

    auto ta_yzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 131);

    auto ta_yzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 132);

    auto ta_yzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 133);

    auto ta_yzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 134);

    auto ta_zzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 135);

    auto ta_zzz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 136);

    auto ta_zzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 137);

    auto ta_zzz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 138);

    auto ta_zzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 139);

    auto ta_zzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 140);

    auto ta_zzz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 141);

    auto ta_zzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 142);

    auto ta_zzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 143);

    auto ta_zzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 144);

    auto ta_zzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 145);

    auto ta_zzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 146);

    auto ta_zzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 147);

    auto ta_zzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 148);

    auto ta_zzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 149);

    // Set up 0-15 components of targeted buffer : GG

    auto ta_xxxx_xxxx_0 = pbuffer.data(idx_npot_0_gg);

    auto ta_xxxx_xxxy_0 = pbuffer.data(idx_npot_0_gg + 1);

    auto ta_xxxx_xxxz_0 = pbuffer.data(idx_npot_0_gg + 2);

    auto ta_xxxx_xxyy_0 = pbuffer.data(idx_npot_0_gg + 3);

    auto ta_xxxx_xxyz_0 = pbuffer.data(idx_npot_0_gg + 4);

    auto ta_xxxx_xxzz_0 = pbuffer.data(idx_npot_0_gg + 5);

    auto ta_xxxx_xyyy_0 = pbuffer.data(idx_npot_0_gg + 6);

    auto ta_xxxx_xyyz_0 = pbuffer.data(idx_npot_0_gg + 7);

    auto ta_xxxx_xyzz_0 = pbuffer.data(idx_npot_0_gg + 8);

    auto ta_xxxx_xzzz_0 = pbuffer.data(idx_npot_0_gg + 9);

    auto ta_xxxx_yyyy_0 = pbuffer.data(idx_npot_0_gg + 10);

    auto ta_xxxx_yyyz_0 = pbuffer.data(idx_npot_0_gg + 11);

    auto ta_xxxx_yyzz_0 = pbuffer.data(idx_npot_0_gg + 12);

    auto ta_xxxx_yzzz_0 = pbuffer.data(idx_npot_0_gg + 13);

    auto ta_xxxx_zzzz_0 = pbuffer.data(idx_npot_0_gg + 14);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xx_xxxx_0,   \
                             ta_xx_xxxx_1,   \
                             ta_xx_xxxy_0,   \
                             ta_xx_xxxy_1,   \
                             ta_xx_xxxz_0,   \
                             ta_xx_xxxz_1,   \
                             ta_xx_xxyy_0,   \
                             ta_xx_xxyy_1,   \
                             ta_xx_xxyz_0,   \
                             ta_xx_xxyz_1,   \
                             ta_xx_xxzz_0,   \
                             ta_xx_xxzz_1,   \
                             ta_xx_xyyy_0,   \
                             ta_xx_xyyy_1,   \
                             ta_xx_xyyz_0,   \
                             ta_xx_xyyz_1,   \
                             ta_xx_xyzz_0,   \
                             ta_xx_xyzz_1,   \
                             ta_xx_xzzz_0,   \
                             ta_xx_xzzz_1,   \
                             ta_xx_yyyy_0,   \
                             ta_xx_yyyy_1,   \
                             ta_xx_yyyz_0,   \
                             ta_xx_yyyz_1,   \
                             ta_xx_yyzz_0,   \
                             ta_xx_yyzz_1,   \
                             ta_xx_yzzz_0,   \
                             ta_xx_yzzz_1,   \
                             ta_xx_zzzz_0,   \
                             ta_xx_zzzz_1,   \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxxx_0,  \
                             ta_xxx_xxxx_1,  \
                             ta_xxx_xxxy_0,  \
                             ta_xxx_xxxy_1,  \
                             ta_xxx_xxxz_0,  \
                             ta_xxx_xxxz_1,  \
                             ta_xxx_xxy_0,   \
                             ta_xxx_xxy_1,   \
                             ta_xxx_xxyy_0,  \
                             ta_xxx_xxyy_1,  \
                             ta_xxx_xxyz_0,  \
                             ta_xxx_xxyz_1,  \
                             ta_xxx_xxz_0,   \
                             ta_xxx_xxz_1,   \
                             ta_xxx_xxzz_0,  \
                             ta_xxx_xxzz_1,  \
                             ta_xxx_xyy_0,   \
                             ta_xxx_xyy_1,   \
                             ta_xxx_xyyy_0,  \
                             ta_xxx_xyyy_1,  \
                             ta_xxx_xyyz_0,  \
                             ta_xxx_xyyz_1,  \
                             ta_xxx_xyz_0,   \
                             ta_xxx_xyz_1,   \
                             ta_xxx_xyzz_0,  \
                             ta_xxx_xyzz_1,  \
                             ta_xxx_xzz_0,   \
                             ta_xxx_xzz_1,   \
                             ta_xxx_xzzz_0,  \
                             ta_xxx_xzzz_1,  \
                             ta_xxx_yyy_0,   \
                             ta_xxx_yyy_1,   \
                             ta_xxx_yyyy_0,  \
                             ta_xxx_yyyy_1,  \
                             ta_xxx_yyyz_0,  \
                             ta_xxx_yyyz_1,  \
                             ta_xxx_yyz_0,   \
                             ta_xxx_yyz_1,   \
                             ta_xxx_yyzz_0,  \
                             ta_xxx_yyzz_1,  \
                             ta_xxx_yzz_0,   \
                             ta_xxx_yzz_1,   \
                             ta_xxx_yzzz_0,  \
                             ta_xxx_yzzz_1,  \
                             ta_xxx_zzz_0,   \
                             ta_xxx_zzz_1,   \
                             ta_xxx_zzzz_0,  \
                             ta_xxx_zzzz_1,  \
                             ta_xxxx_xxxx_0, \
                             ta_xxxx_xxxy_0, \
                             ta_xxxx_xxxz_0, \
                             ta_xxxx_xxyy_0, \
                             ta_xxxx_xxyz_0, \
                             ta_xxxx_xxzz_0, \
                             ta_xxxx_xyyy_0, \
                             ta_xxxx_xyyz_0, \
                             ta_xxxx_xyzz_0, \
                             ta_xxxx_xzzz_0, \
                             ta_xxxx_yyyy_0, \
                             ta_xxxx_yyyz_0, \
                             ta_xxxx_yyzz_0, \
                             ta_xxxx_yzzz_0, \
                             ta_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_xxxx_0[i] = 3.0 * ta_xx_xxxx_0[i] * fe_0 - 3.0 * ta_xx_xxxx_1[i] * fe_0 + 4.0 * ta_xxx_xxx_0[i] * fe_0 -
                            4.0 * ta_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxx_0[i] * pa_x[i] - ta_xxx_xxxx_1[i] * pc_x[i];

        ta_xxxx_xxxy_0[i] = 3.0 * ta_xx_xxxy_0[i] * fe_0 - 3.0 * ta_xx_xxxy_1[i] * fe_0 + 3.0 * ta_xxx_xxy_0[i] * fe_0 -
                            3.0 * ta_xxx_xxy_1[i] * fe_0 + ta_xxx_xxxy_0[i] * pa_x[i] - ta_xxx_xxxy_1[i] * pc_x[i];

        ta_xxxx_xxxz_0[i] = 3.0 * ta_xx_xxxz_0[i] * fe_0 - 3.0 * ta_xx_xxxz_1[i] * fe_0 + 3.0 * ta_xxx_xxz_0[i] * fe_0 -
                            3.0 * ta_xxx_xxz_1[i] * fe_0 + ta_xxx_xxxz_0[i] * pa_x[i] - ta_xxx_xxxz_1[i] * pc_x[i];

        ta_xxxx_xxyy_0[i] = 3.0 * ta_xx_xxyy_0[i] * fe_0 - 3.0 * ta_xx_xxyy_1[i] * fe_0 + 2.0 * ta_xxx_xyy_0[i] * fe_0 -
                            2.0 * ta_xxx_xyy_1[i] * fe_0 + ta_xxx_xxyy_0[i] * pa_x[i] - ta_xxx_xxyy_1[i] * pc_x[i];

        ta_xxxx_xxyz_0[i] = 3.0 * ta_xx_xxyz_0[i] * fe_0 - 3.0 * ta_xx_xxyz_1[i] * fe_0 + 2.0 * ta_xxx_xyz_0[i] * fe_0 -
                            2.0 * ta_xxx_xyz_1[i] * fe_0 + ta_xxx_xxyz_0[i] * pa_x[i] - ta_xxx_xxyz_1[i] * pc_x[i];

        ta_xxxx_xxzz_0[i] = 3.0 * ta_xx_xxzz_0[i] * fe_0 - 3.0 * ta_xx_xxzz_1[i] * fe_0 + 2.0 * ta_xxx_xzz_0[i] * fe_0 -
                            2.0 * ta_xxx_xzz_1[i] * fe_0 + ta_xxx_xxzz_0[i] * pa_x[i] - ta_xxx_xxzz_1[i] * pc_x[i];

        ta_xxxx_xyyy_0[i] = 3.0 * ta_xx_xyyy_0[i] * fe_0 - 3.0 * ta_xx_xyyy_1[i] * fe_0 + ta_xxx_yyy_0[i] * fe_0 - ta_xxx_yyy_1[i] * fe_0 +
                            ta_xxx_xyyy_0[i] * pa_x[i] - ta_xxx_xyyy_1[i] * pc_x[i];

        ta_xxxx_xyyz_0[i] = 3.0 * ta_xx_xyyz_0[i] * fe_0 - 3.0 * ta_xx_xyyz_1[i] * fe_0 + ta_xxx_yyz_0[i] * fe_0 - ta_xxx_yyz_1[i] * fe_0 +
                            ta_xxx_xyyz_0[i] * pa_x[i] - ta_xxx_xyyz_1[i] * pc_x[i];

        ta_xxxx_xyzz_0[i] = 3.0 * ta_xx_xyzz_0[i] * fe_0 - 3.0 * ta_xx_xyzz_1[i] * fe_0 + ta_xxx_yzz_0[i] * fe_0 - ta_xxx_yzz_1[i] * fe_0 +
                            ta_xxx_xyzz_0[i] * pa_x[i] - ta_xxx_xyzz_1[i] * pc_x[i];

        ta_xxxx_xzzz_0[i] = 3.0 * ta_xx_xzzz_0[i] * fe_0 - 3.0 * ta_xx_xzzz_1[i] * fe_0 + ta_xxx_zzz_0[i] * fe_0 - ta_xxx_zzz_1[i] * fe_0 +
                            ta_xxx_xzzz_0[i] * pa_x[i] - ta_xxx_xzzz_1[i] * pc_x[i];

        ta_xxxx_yyyy_0[i] = 3.0 * ta_xx_yyyy_0[i] * fe_0 - 3.0 * ta_xx_yyyy_1[i] * fe_0 + ta_xxx_yyyy_0[i] * pa_x[i] - ta_xxx_yyyy_1[i] * pc_x[i];

        ta_xxxx_yyyz_0[i] = 3.0 * ta_xx_yyyz_0[i] * fe_0 - 3.0 * ta_xx_yyyz_1[i] * fe_0 + ta_xxx_yyyz_0[i] * pa_x[i] - ta_xxx_yyyz_1[i] * pc_x[i];

        ta_xxxx_yyzz_0[i] = 3.0 * ta_xx_yyzz_0[i] * fe_0 - 3.0 * ta_xx_yyzz_1[i] * fe_0 + ta_xxx_yyzz_0[i] * pa_x[i] - ta_xxx_yyzz_1[i] * pc_x[i];

        ta_xxxx_yzzz_0[i] = 3.0 * ta_xx_yzzz_0[i] * fe_0 - 3.0 * ta_xx_yzzz_1[i] * fe_0 + ta_xxx_yzzz_0[i] * pa_x[i] - ta_xxx_yzzz_1[i] * pc_x[i];

        ta_xxxx_zzzz_0[i] = 3.0 * ta_xx_zzzz_0[i] * fe_0 - 3.0 * ta_xx_zzzz_1[i] * fe_0 + ta_xxx_zzzz_0[i] * pa_x[i] - ta_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

    auto ta_xxxy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 15);

    auto ta_xxxy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 16);

    auto ta_xxxy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 17);

    auto ta_xxxy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 18);

    auto ta_xxxy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 19);

    auto ta_xxxy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 20);

    auto ta_xxxy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 21);

    auto ta_xxxy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 22);

    auto ta_xxxy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 23);

    auto ta_xxxy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 24);

    auto ta_xxxy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 25);

    auto ta_xxxy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 26);

    auto ta_xxxy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 27);

    auto ta_xxxy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 28);

    auto ta_xxxy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxxx_0,  \
                             ta_xxx_xxxx_1,  \
                             ta_xxx_xxxy_0,  \
                             ta_xxx_xxxy_1,  \
                             ta_xxx_xxxz_0,  \
                             ta_xxx_xxxz_1,  \
                             ta_xxx_xxy_0,   \
                             ta_xxx_xxy_1,   \
                             ta_xxx_xxyy_0,  \
                             ta_xxx_xxyy_1,  \
                             ta_xxx_xxyz_0,  \
                             ta_xxx_xxyz_1,  \
                             ta_xxx_xxz_0,   \
                             ta_xxx_xxz_1,   \
                             ta_xxx_xxzz_0,  \
                             ta_xxx_xxzz_1,  \
                             ta_xxx_xyy_0,   \
                             ta_xxx_xyy_1,   \
                             ta_xxx_xyyy_0,  \
                             ta_xxx_xyyy_1,  \
                             ta_xxx_xyyz_0,  \
                             ta_xxx_xyyz_1,  \
                             ta_xxx_xyz_0,   \
                             ta_xxx_xyz_1,   \
                             ta_xxx_xyzz_0,  \
                             ta_xxx_xyzz_1,  \
                             ta_xxx_xzz_0,   \
                             ta_xxx_xzz_1,   \
                             ta_xxx_xzzz_0,  \
                             ta_xxx_xzzz_1,  \
                             ta_xxx_zzzz_0,  \
                             ta_xxx_zzzz_1,  \
                             ta_xxxy_xxxx_0, \
                             ta_xxxy_xxxy_0, \
                             ta_xxxy_xxxz_0, \
                             ta_xxxy_xxyy_0, \
                             ta_xxxy_xxyz_0, \
                             ta_xxxy_xxzz_0, \
                             ta_xxxy_xyyy_0, \
                             ta_xxxy_xyyz_0, \
                             ta_xxxy_xyzz_0, \
                             ta_xxxy_xzzz_0, \
                             ta_xxxy_yyyy_0, \
                             ta_xxxy_yyyz_0, \
                             ta_xxxy_yyzz_0, \
                             ta_xxxy_yzzz_0, \
                             ta_xxxy_zzzz_0, \
                             ta_xxy_yyyy_0,  \
                             ta_xxy_yyyy_1,  \
                             ta_xxy_yyyz_0,  \
                             ta_xxy_yyyz_1,  \
                             ta_xxy_yyzz_0,  \
                             ta_xxy_yyzz_1,  \
                             ta_xxy_yzzz_0,  \
                             ta_xxy_yzzz_1,  \
                             ta_xy_yyyy_0,   \
                             ta_xy_yyyy_1,   \
                             ta_xy_yyyz_0,   \
                             ta_xy_yyyz_1,   \
                             ta_xy_yyzz_0,   \
                             ta_xy_yyzz_1,   \
                             ta_xy_yzzz_0,   \
                             ta_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_xxxx_0[i] = ta_xxx_xxxx_0[i] * pa_y[i] - ta_xxx_xxxx_1[i] * pc_y[i];

        ta_xxxy_xxxy_0[i] = ta_xxx_xxx_0[i] * fe_0 - ta_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxy_0[i] * pa_y[i] - ta_xxx_xxxy_1[i] * pc_y[i];

        ta_xxxy_xxxz_0[i] = ta_xxx_xxxz_0[i] * pa_y[i] - ta_xxx_xxxz_1[i] * pc_y[i];

        ta_xxxy_xxyy_0[i] = 2.0 * ta_xxx_xxy_0[i] * fe_0 - 2.0 * ta_xxx_xxy_1[i] * fe_0 + ta_xxx_xxyy_0[i] * pa_y[i] - ta_xxx_xxyy_1[i] * pc_y[i];

        ta_xxxy_xxyz_0[i] = ta_xxx_xxz_0[i] * fe_0 - ta_xxx_xxz_1[i] * fe_0 + ta_xxx_xxyz_0[i] * pa_y[i] - ta_xxx_xxyz_1[i] * pc_y[i];

        ta_xxxy_xxzz_0[i] = ta_xxx_xxzz_0[i] * pa_y[i] - ta_xxx_xxzz_1[i] * pc_y[i];

        ta_xxxy_xyyy_0[i] = 3.0 * ta_xxx_xyy_0[i] * fe_0 - 3.0 * ta_xxx_xyy_1[i] * fe_0 + ta_xxx_xyyy_0[i] * pa_y[i] - ta_xxx_xyyy_1[i] * pc_y[i];

        ta_xxxy_xyyz_0[i] = 2.0 * ta_xxx_xyz_0[i] * fe_0 - 2.0 * ta_xxx_xyz_1[i] * fe_0 + ta_xxx_xyyz_0[i] * pa_y[i] - ta_xxx_xyyz_1[i] * pc_y[i];

        ta_xxxy_xyzz_0[i] = ta_xxx_xzz_0[i] * fe_0 - ta_xxx_xzz_1[i] * fe_0 + ta_xxx_xyzz_0[i] * pa_y[i] - ta_xxx_xyzz_1[i] * pc_y[i];

        ta_xxxy_xzzz_0[i] = ta_xxx_xzzz_0[i] * pa_y[i] - ta_xxx_xzzz_1[i] * pc_y[i];

        ta_xxxy_yyyy_0[i] = 2.0 * ta_xy_yyyy_0[i] * fe_0 - 2.0 * ta_xy_yyyy_1[i] * fe_0 + ta_xxy_yyyy_0[i] * pa_x[i] - ta_xxy_yyyy_1[i] * pc_x[i];

        ta_xxxy_yyyz_0[i] = 2.0 * ta_xy_yyyz_0[i] * fe_0 - 2.0 * ta_xy_yyyz_1[i] * fe_0 + ta_xxy_yyyz_0[i] * pa_x[i] - ta_xxy_yyyz_1[i] * pc_x[i];

        ta_xxxy_yyzz_0[i] = 2.0 * ta_xy_yyzz_0[i] * fe_0 - 2.0 * ta_xy_yyzz_1[i] * fe_0 + ta_xxy_yyzz_0[i] * pa_x[i] - ta_xxy_yyzz_1[i] * pc_x[i];

        ta_xxxy_yzzz_0[i] = 2.0 * ta_xy_yzzz_0[i] * fe_0 - 2.0 * ta_xy_yzzz_1[i] * fe_0 + ta_xxy_yzzz_0[i] * pa_x[i] - ta_xxy_yzzz_1[i] * pc_x[i];

        ta_xxxy_zzzz_0[i] = ta_xxx_zzzz_0[i] * pa_y[i] - ta_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

    auto ta_xxxz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 30);

    auto ta_xxxz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 31);

    auto ta_xxxz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 32);

    auto ta_xxxz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 33);

    auto ta_xxxz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 34);

    auto ta_xxxz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 35);

    auto ta_xxxz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 36);

    auto ta_xxxz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 37);

    auto ta_xxxz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 38);

    auto ta_xxxz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 39);

    auto ta_xxxz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 40);

    auto ta_xxxz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 41);

    auto ta_xxxz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 42);

    auto ta_xxxz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 43);

    auto ta_xxxz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 44);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxxx_0,  \
                             ta_xxx_xxxx_1,  \
                             ta_xxx_xxxy_0,  \
                             ta_xxx_xxxy_1,  \
                             ta_xxx_xxxz_0,  \
                             ta_xxx_xxxz_1,  \
                             ta_xxx_xxy_0,   \
                             ta_xxx_xxy_1,   \
                             ta_xxx_xxyy_0,  \
                             ta_xxx_xxyy_1,  \
                             ta_xxx_xxyz_0,  \
                             ta_xxx_xxyz_1,  \
                             ta_xxx_xxz_0,   \
                             ta_xxx_xxz_1,   \
                             ta_xxx_xxzz_0,  \
                             ta_xxx_xxzz_1,  \
                             ta_xxx_xyy_0,   \
                             ta_xxx_xyy_1,   \
                             ta_xxx_xyyy_0,  \
                             ta_xxx_xyyy_1,  \
                             ta_xxx_xyyz_0,  \
                             ta_xxx_xyyz_1,  \
                             ta_xxx_xyz_0,   \
                             ta_xxx_xyz_1,   \
                             ta_xxx_xyzz_0,  \
                             ta_xxx_xyzz_1,  \
                             ta_xxx_xzz_0,   \
                             ta_xxx_xzz_1,   \
                             ta_xxx_xzzz_0,  \
                             ta_xxx_xzzz_1,  \
                             ta_xxx_yyyy_0,  \
                             ta_xxx_yyyy_1,  \
                             ta_xxxz_xxxx_0, \
                             ta_xxxz_xxxy_0, \
                             ta_xxxz_xxxz_0, \
                             ta_xxxz_xxyy_0, \
                             ta_xxxz_xxyz_0, \
                             ta_xxxz_xxzz_0, \
                             ta_xxxz_xyyy_0, \
                             ta_xxxz_xyyz_0, \
                             ta_xxxz_xyzz_0, \
                             ta_xxxz_xzzz_0, \
                             ta_xxxz_yyyy_0, \
                             ta_xxxz_yyyz_0, \
                             ta_xxxz_yyzz_0, \
                             ta_xxxz_yzzz_0, \
                             ta_xxxz_zzzz_0, \
                             ta_xxz_yyyz_0,  \
                             ta_xxz_yyyz_1,  \
                             ta_xxz_yyzz_0,  \
                             ta_xxz_yyzz_1,  \
                             ta_xxz_yzzz_0,  \
                             ta_xxz_yzzz_1,  \
                             ta_xxz_zzzz_0,  \
                             ta_xxz_zzzz_1,  \
                             ta_xz_yyyz_0,   \
                             ta_xz_yyyz_1,   \
                             ta_xz_yyzz_0,   \
                             ta_xz_yyzz_1,   \
                             ta_xz_yzzz_0,   \
                             ta_xz_yzzz_1,   \
                             ta_xz_zzzz_0,   \
                             ta_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_xxxx_0[i] = ta_xxx_xxxx_0[i] * pa_z[i] - ta_xxx_xxxx_1[i] * pc_z[i];

        ta_xxxz_xxxy_0[i] = ta_xxx_xxxy_0[i] * pa_z[i] - ta_xxx_xxxy_1[i] * pc_z[i];

        ta_xxxz_xxxz_0[i] = ta_xxx_xxx_0[i] * fe_0 - ta_xxx_xxx_1[i] * fe_0 + ta_xxx_xxxz_0[i] * pa_z[i] - ta_xxx_xxxz_1[i] * pc_z[i];

        ta_xxxz_xxyy_0[i] = ta_xxx_xxyy_0[i] * pa_z[i] - ta_xxx_xxyy_1[i] * pc_z[i];

        ta_xxxz_xxyz_0[i] = ta_xxx_xxy_0[i] * fe_0 - ta_xxx_xxy_1[i] * fe_0 + ta_xxx_xxyz_0[i] * pa_z[i] - ta_xxx_xxyz_1[i] * pc_z[i];

        ta_xxxz_xxzz_0[i] = 2.0 * ta_xxx_xxz_0[i] * fe_0 - 2.0 * ta_xxx_xxz_1[i] * fe_0 + ta_xxx_xxzz_0[i] * pa_z[i] - ta_xxx_xxzz_1[i] * pc_z[i];

        ta_xxxz_xyyy_0[i] = ta_xxx_xyyy_0[i] * pa_z[i] - ta_xxx_xyyy_1[i] * pc_z[i];

        ta_xxxz_xyyz_0[i] = ta_xxx_xyy_0[i] * fe_0 - ta_xxx_xyy_1[i] * fe_0 + ta_xxx_xyyz_0[i] * pa_z[i] - ta_xxx_xyyz_1[i] * pc_z[i];

        ta_xxxz_xyzz_0[i] = 2.0 * ta_xxx_xyz_0[i] * fe_0 - 2.0 * ta_xxx_xyz_1[i] * fe_0 + ta_xxx_xyzz_0[i] * pa_z[i] - ta_xxx_xyzz_1[i] * pc_z[i];

        ta_xxxz_xzzz_0[i] = 3.0 * ta_xxx_xzz_0[i] * fe_0 - 3.0 * ta_xxx_xzz_1[i] * fe_0 + ta_xxx_xzzz_0[i] * pa_z[i] - ta_xxx_xzzz_1[i] * pc_z[i];

        ta_xxxz_yyyy_0[i] = ta_xxx_yyyy_0[i] * pa_z[i] - ta_xxx_yyyy_1[i] * pc_z[i];

        ta_xxxz_yyyz_0[i] = 2.0 * ta_xz_yyyz_0[i] * fe_0 - 2.0 * ta_xz_yyyz_1[i] * fe_0 + ta_xxz_yyyz_0[i] * pa_x[i] - ta_xxz_yyyz_1[i] * pc_x[i];

        ta_xxxz_yyzz_0[i] = 2.0 * ta_xz_yyzz_0[i] * fe_0 - 2.0 * ta_xz_yyzz_1[i] * fe_0 + ta_xxz_yyzz_0[i] * pa_x[i] - ta_xxz_yyzz_1[i] * pc_x[i];

        ta_xxxz_yzzz_0[i] = 2.0 * ta_xz_yzzz_0[i] * fe_0 - 2.0 * ta_xz_yzzz_1[i] * fe_0 + ta_xxz_yzzz_0[i] * pa_x[i] - ta_xxz_yzzz_1[i] * pc_x[i];

        ta_xxxz_zzzz_0[i] = 2.0 * ta_xz_zzzz_0[i] * fe_0 - 2.0 * ta_xz_zzzz_1[i] * fe_0 + ta_xxz_zzzz_0[i] * pa_x[i] - ta_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : GG

    auto ta_xxyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 45);

    auto ta_xxyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 46);

    auto ta_xxyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 47);

    auto ta_xxyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 48);

    auto ta_xxyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 49);

    auto ta_xxyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 50);

    auto ta_xxyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 51);

    auto ta_xxyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 52);

    auto ta_xxyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 53);

    auto ta_xxyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 54);

    auto ta_xxyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 55);

    auto ta_xxyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 56);

    auto ta_xxyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 57);

    auto ta_xxyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 58);

    auto ta_xxyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xx_xxxx_0,   \
                             ta_xx_xxxx_1,   \
                             ta_xx_xxxz_0,   \
                             ta_xx_xxxz_1,   \
                             ta_xx_xxzz_0,   \
                             ta_xx_xxzz_1,   \
                             ta_xx_xzzz_0,   \
                             ta_xx_xzzz_1,   \
                             ta_xxy_xxxx_0,  \
                             ta_xxy_xxxx_1,  \
                             ta_xxy_xxxz_0,  \
                             ta_xxy_xxxz_1,  \
                             ta_xxy_xxzz_0,  \
                             ta_xxy_xxzz_1,  \
                             ta_xxy_xzzz_0,  \
                             ta_xxy_xzzz_1,  \
                             ta_xxyy_xxxx_0, \
                             ta_xxyy_xxxy_0, \
                             ta_xxyy_xxxz_0, \
                             ta_xxyy_xxyy_0, \
                             ta_xxyy_xxyz_0, \
                             ta_xxyy_xxzz_0, \
                             ta_xxyy_xyyy_0, \
                             ta_xxyy_xyyz_0, \
                             ta_xxyy_xyzz_0, \
                             ta_xxyy_xzzz_0, \
                             ta_xxyy_yyyy_0, \
                             ta_xxyy_yyyz_0, \
                             ta_xxyy_yyzz_0, \
                             ta_xxyy_yzzz_0, \
                             ta_xxyy_zzzz_0, \
                             ta_xyy_xxxy_0,  \
                             ta_xyy_xxxy_1,  \
                             ta_xyy_xxy_0,   \
                             ta_xyy_xxy_1,   \
                             ta_xyy_xxyy_0,  \
                             ta_xyy_xxyy_1,  \
                             ta_xyy_xxyz_0,  \
                             ta_xyy_xxyz_1,  \
                             ta_xyy_xyy_0,   \
                             ta_xyy_xyy_1,   \
                             ta_xyy_xyyy_0,  \
                             ta_xyy_xyyy_1,  \
                             ta_xyy_xyyz_0,  \
                             ta_xyy_xyyz_1,  \
                             ta_xyy_xyz_0,   \
                             ta_xyy_xyz_1,   \
                             ta_xyy_xyzz_0,  \
                             ta_xyy_xyzz_1,  \
                             ta_xyy_yyy_0,   \
                             ta_xyy_yyy_1,   \
                             ta_xyy_yyyy_0,  \
                             ta_xyy_yyyy_1,  \
                             ta_xyy_yyyz_0,  \
                             ta_xyy_yyyz_1,  \
                             ta_xyy_yyz_0,   \
                             ta_xyy_yyz_1,   \
                             ta_xyy_yyzz_0,  \
                             ta_xyy_yyzz_1,  \
                             ta_xyy_yzz_0,   \
                             ta_xyy_yzz_1,   \
                             ta_xyy_yzzz_0,  \
                             ta_xyy_yzzz_1,  \
                             ta_xyy_zzzz_0,  \
                             ta_xyy_zzzz_1,  \
                             ta_yy_xxxy_0,   \
                             ta_yy_xxxy_1,   \
                             ta_yy_xxyy_0,   \
                             ta_yy_xxyy_1,   \
                             ta_yy_xxyz_0,   \
                             ta_yy_xxyz_1,   \
                             ta_yy_xyyy_0,   \
                             ta_yy_xyyy_1,   \
                             ta_yy_xyyz_0,   \
                             ta_yy_xyyz_1,   \
                             ta_yy_xyzz_0,   \
                             ta_yy_xyzz_1,   \
                             ta_yy_yyyy_0,   \
                             ta_yy_yyyy_1,   \
                             ta_yy_yyyz_0,   \
                             ta_yy_yyyz_1,   \
                             ta_yy_yyzz_0,   \
                             ta_yy_yyzz_1,   \
                             ta_yy_yzzz_0,   \
                             ta_yy_yzzz_1,   \
                             ta_yy_zzzz_0,   \
                             ta_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_xxxx_0[i] = ta_xx_xxxx_0[i] * fe_0 - ta_xx_xxxx_1[i] * fe_0 + ta_xxy_xxxx_0[i] * pa_y[i] - ta_xxy_xxxx_1[i] * pc_y[i];

        ta_xxyy_xxxy_0[i] = ta_yy_xxxy_0[i] * fe_0 - ta_yy_xxxy_1[i] * fe_0 + 3.0 * ta_xyy_xxy_0[i] * fe_0 - 3.0 * ta_xyy_xxy_1[i] * fe_0 +
                            ta_xyy_xxxy_0[i] * pa_x[i] - ta_xyy_xxxy_1[i] * pc_x[i];

        ta_xxyy_xxxz_0[i] = ta_xx_xxxz_0[i] * fe_0 - ta_xx_xxxz_1[i] * fe_0 + ta_xxy_xxxz_0[i] * pa_y[i] - ta_xxy_xxxz_1[i] * pc_y[i];

        ta_xxyy_xxyy_0[i] = ta_yy_xxyy_0[i] * fe_0 - ta_yy_xxyy_1[i] * fe_0 + 2.0 * ta_xyy_xyy_0[i] * fe_0 - 2.0 * ta_xyy_xyy_1[i] * fe_0 +
                            ta_xyy_xxyy_0[i] * pa_x[i] - ta_xyy_xxyy_1[i] * pc_x[i];

        ta_xxyy_xxyz_0[i] = ta_yy_xxyz_0[i] * fe_0 - ta_yy_xxyz_1[i] * fe_0 + 2.0 * ta_xyy_xyz_0[i] * fe_0 - 2.0 * ta_xyy_xyz_1[i] * fe_0 +
                            ta_xyy_xxyz_0[i] * pa_x[i] - ta_xyy_xxyz_1[i] * pc_x[i];

        ta_xxyy_xxzz_0[i] = ta_xx_xxzz_0[i] * fe_0 - ta_xx_xxzz_1[i] * fe_0 + ta_xxy_xxzz_0[i] * pa_y[i] - ta_xxy_xxzz_1[i] * pc_y[i];

        ta_xxyy_xyyy_0[i] = ta_yy_xyyy_0[i] * fe_0 - ta_yy_xyyy_1[i] * fe_0 + ta_xyy_yyy_0[i] * fe_0 - ta_xyy_yyy_1[i] * fe_0 +
                            ta_xyy_xyyy_0[i] * pa_x[i] - ta_xyy_xyyy_1[i] * pc_x[i];

        ta_xxyy_xyyz_0[i] = ta_yy_xyyz_0[i] * fe_0 - ta_yy_xyyz_1[i] * fe_0 + ta_xyy_yyz_0[i] * fe_0 - ta_xyy_yyz_1[i] * fe_0 +
                            ta_xyy_xyyz_0[i] * pa_x[i] - ta_xyy_xyyz_1[i] * pc_x[i];

        ta_xxyy_xyzz_0[i] = ta_yy_xyzz_0[i] * fe_0 - ta_yy_xyzz_1[i] * fe_0 + ta_xyy_yzz_0[i] * fe_0 - ta_xyy_yzz_1[i] * fe_0 +
                            ta_xyy_xyzz_0[i] * pa_x[i] - ta_xyy_xyzz_1[i] * pc_x[i];

        ta_xxyy_xzzz_0[i] = ta_xx_xzzz_0[i] * fe_0 - ta_xx_xzzz_1[i] * fe_0 + ta_xxy_xzzz_0[i] * pa_y[i] - ta_xxy_xzzz_1[i] * pc_y[i];

        ta_xxyy_yyyy_0[i] = ta_yy_yyyy_0[i] * fe_0 - ta_yy_yyyy_1[i] * fe_0 + ta_xyy_yyyy_0[i] * pa_x[i] - ta_xyy_yyyy_1[i] * pc_x[i];

        ta_xxyy_yyyz_0[i] = ta_yy_yyyz_0[i] * fe_0 - ta_yy_yyyz_1[i] * fe_0 + ta_xyy_yyyz_0[i] * pa_x[i] - ta_xyy_yyyz_1[i] * pc_x[i];

        ta_xxyy_yyzz_0[i] = ta_yy_yyzz_0[i] * fe_0 - ta_yy_yyzz_1[i] * fe_0 + ta_xyy_yyzz_0[i] * pa_x[i] - ta_xyy_yyzz_1[i] * pc_x[i];

        ta_xxyy_yzzz_0[i] = ta_yy_yzzz_0[i] * fe_0 - ta_yy_yzzz_1[i] * fe_0 + ta_xyy_yzzz_0[i] * pa_x[i] - ta_xyy_yzzz_1[i] * pc_x[i];

        ta_xxyy_zzzz_0[i] = ta_yy_zzzz_0[i] * fe_0 - ta_yy_zzzz_1[i] * fe_0 + ta_xyy_zzzz_0[i] * pa_x[i] - ta_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto ta_xxyz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 60);

    auto ta_xxyz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 61);

    auto ta_xxyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 62);

    auto ta_xxyz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 63);

    auto ta_xxyz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 64);

    auto ta_xxyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 65);

    auto ta_xxyz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 66);

    auto ta_xxyz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 67);

    auto ta_xxyz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 68);

    auto ta_xxyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 69);

    auto ta_xxyz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 70);

    auto ta_xxyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 71);

    auto ta_xxyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 72);

    auto ta_xxyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 73);

    auto ta_xxyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 74);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxy_xxxy_0,  \
                             ta_xxy_xxxy_1,  \
                             ta_xxy_xxyy_0,  \
                             ta_xxy_xxyy_1,  \
                             ta_xxy_xyyy_0,  \
                             ta_xxy_xyyy_1,  \
                             ta_xxy_yyyy_0,  \
                             ta_xxy_yyyy_1,  \
                             ta_xxyz_xxxx_0, \
                             ta_xxyz_xxxy_0, \
                             ta_xxyz_xxxz_0, \
                             ta_xxyz_xxyy_0, \
                             ta_xxyz_xxyz_0, \
                             ta_xxyz_xxzz_0, \
                             ta_xxyz_xyyy_0, \
                             ta_xxyz_xyyz_0, \
                             ta_xxyz_xyzz_0, \
                             ta_xxyz_xzzz_0, \
                             ta_xxyz_yyyy_0, \
                             ta_xxyz_yyyz_0, \
                             ta_xxyz_yyzz_0, \
                             ta_xxyz_yzzz_0, \
                             ta_xxyz_zzzz_0, \
                             ta_xxz_xxxx_0,  \
                             ta_xxz_xxxx_1,  \
                             ta_xxz_xxxz_0,  \
                             ta_xxz_xxxz_1,  \
                             ta_xxz_xxyz_0,  \
                             ta_xxz_xxyz_1,  \
                             ta_xxz_xxz_0,   \
                             ta_xxz_xxz_1,   \
                             ta_xxz_xxzz_0,  \
                             ta_xxz_xxzz_1,  \
                             ta_xxz_xyyz_0,  \
                             ta_xxz_xyyz_1,  \
                             ta_xxz_xyz_0,   \
                             ta_xxz_xyz_1,   \
                             ta_xxz_xyzz_0,  \
                             ta_xxz_xyzz_1,  \
                             ta_xxz_xzz_0,   \
                             ta_xxz_xzz_1,   \
                             ta_xxz_xzzz_0,  \
                             ta_xxz_xzzz_1,  \
                             ta_xxz_zzzz_0,  \
                             ta_xxz_zzzz_1,  \
                             ta_xyz_yyyz_0,  \
                             ta_xyz_yyyz_1,  \
                             ta_xyz_yyzz_0,  \
                             ta_xyz_yyzz_1,  \
                             ta_xyz_yzzz_0,  \
                             ta_xyz_yzzz_1,  \
                             ta_yz_yyyz_0,   \
                             ta_yz_yyyz_1,   \
                             ta_yz_yyzz_0,   \
                             ta_yz_yyzz_1,   \
                             ta_yz_yzzz_0,   \
                             ta_yz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyz_xxxx_0[i] = ta_xxz_xxxx_0[i] * pa_y[i] - ta_xxz_xxxx_1[i] * pc_y[i];

        ta_xxyz_xxxy_0[i] = ta_xxy_xxxy_0[i] * pa_z[i] - ta_xxy_xxxy_1[i] * pc_z[i];

        ta_xxyz_xxxz_0[i] = ta_xxz_xxxz_0[i] * pa_y[i] - ta_xxz_xxxz_1[i] * pc_y[i];

        ta_xxyz_xxyy_0[i] = ta_xxy_xxyy_0[i] * pa_z[i] - ta_xxy_xxyy_1[i] * pc_z[i];

        ta_xxyz_xxyz_0[i] = ta_xxz_xxz_0[i] * fe_0 - ta_xxz_xxz_1[i] * fe_0 + ta_xxz_xxyz_0[i] * pa_y[i] - ta_xxz_xxyz_1[i] * pc_y[i];

        ta_xxyz_xxzz_0[i] = ta_xxz_xxzz_0[i] * pa_y[i] - ta_xxz_xxzz_1[i] * pc_y[i];

        ta_xxyz_xyyy_0[i] = ta_xxy_xyyy_0[i] * pa_z[i] - ta_xxy_xyyy_1[i] * pc_z[i];

        ta_xxyz_xyyz_0[i] = 2.0 * ta_xxz_xyz_0[i] * fe_0 - 2.0 * ta_xxz_xyz_1[i] * fe_0 + ta_xxz_xyyz_0[i] * pa_y[i] - ta_xxz_xyyz_1[i] * pc_y[i];

        ta_xxyz_xyzz_0[i] = ta_xxz_xzz_0[i] * fe_0 - ta_xxz_xzz_1[i] * fe_0 + ta_xxz_xyzz_0[i] * pa_y[i] - ta_xxz_xyzz_1[i] * pc_y[i];

        ta_xxyz_xzzz_0[i] = ta_xxz_xzzz_0[i] * pa_y[i] - ta_xxz_xzzz_1[i] * pc_y[i];

        ta_xxyz_yyyy_0[i] = ta_xxy_yyyy_0[i] * pa_z[i] - ta_xxy_yyyy_1[i] * pc_z[i];

        ta_xxyz_yyyz_0[i] = ta_yz_yyyz_0[i] * fe_0 - ta_yz_yyyz_1[i] * fe_0 + ta_xyz_yyyz_0[i] * pa_x[i] - ta_xyz_yyyz_1[i] * pc_x[i];

        ta_xxyz_yyzz_0[i] = ta_yz_yyzz_0[i] * fe_0 - ta_yz_yyzz_1[i] * fe_0 + ta_xyz_yyzz_0[i] * pa_x[i] - ta_xyz_yyzz_1[i] * pc_x[i];

        ta_xxyz_yzzz_0[i] = ta_yz_yzzz_0[i] * fe_0 - ta_yz_yzzz_1[i] * fe_0 + ta_xyz_yzzz_0[i] * pa_x[i] - ta_xyz_yzzz_1[i] * pc_x[i];

        ta_xxyz_zzzz_0[i] = ta_xxz_zzzz_0[i] * pa_y[i] - ta_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

    auto ta_xxzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 75);

    auto ta_xxzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 76);

    auto ta_xxzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 77);

    auto ta_xxzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 78);

    auto ta_xxzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 79);

    auto ta_xxzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 80);

    auto ta_xxzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 81);

    auto ta_xxzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 82);

    auto ta_xxzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 83);

    auto ta_xxzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 84);

    auto ta_xxzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 85);

    auto ta_xxzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 86);

    auto ta_xxzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 87);

    auto ta_xxzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 88);

    auto ta_xxzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xx_xxxx_0,   \
                             ta_xx_xxxx_1,   \
                             ta_xx_xxxy_0,   \
                             ta_xx_xxxy_1,   \
                             ta_xx_xxyy_0,   \
                             ta_xx_xxyy_1,   \
                             ta_xx_xyyy_0,   \
                             ta_xx_xyyy_1,   \
                             ta_xxz_xxxx_0,  \
                             ta_xxz_xxxx_1,  \
                             ta_xxz_xxxy_0,  \
                             ta_xxz_xxxy_1,  \
                             ta_xxz_xxyy_0,  \
                             ta_xxz_xxyy_1,  \
                             ta_xxz_xyyy_0,  \
                             ta_xxz_xyyy_1,  \
                             ta_xxzz_xxxx_0, \
                             ta_xxzz_xxxy_0, \
                             ta_xxzz_xxxz_0, \
                             ta_xxzz_xxyy_0, \
                             ta_xxzz_xxyz_0, \
                             ta_xxzz_xxzz_0, \
                             ta_xxzz_xyyy_0, \
                             ta_xxzz_xyyz_0, \
                             ta_xxzz_xyzz_0, \
                             ta_xxzz_xzzz_0, \
                             ta_xxzz_yyyy_0, \
                             ta_xxzz_yyyz_0, \
                             ta_xxzz_yyzz_0, \
                             ta_xxzz_yzzz_0, \
                             ta_xxzz_zzzz_0, \
                             ta_xzz_xxxz_0,  \
                             ta_xzz_xxxz_1,  \
                             ta_xzz_xxyz_0,  \
                             ta_xzz_xxyz_1,  \
                             ta_xzz_xxz_0,   \
                             ta_xzz_xxz_1,   \
                             ta_xzz_xxzz_0,  \
                             ta_xzz_xxzz_1,  \
                             ta_xzz_xyyz_0,  \
                             ta_xzz_xyyz_1,  \
                             ta_xzz_xyz_0,   \
                             ta_xzz_xyz_1,   \
                             ta_xzz_xyzz_0,  \
                             ta_xzz_xyzz_1,  \
                             ta_xzz_xzz_0,   \
                             ta_xzz_xzz_1,   \
                             ta_xzz_xzzz_0,  \
                             ta_xzz_xzzz_1,  \
                             ta_xzz_yyyy_0,  \
                             ta_xzz_yyyy_1,  \
                             ta_xzz_yyyz_0,  \
                             ta_xzz_yyyz_1,  \
                             ta_xzz_yyz_0,   \
                             ta_xzz_yyz_1,   \
                             ta_xzz_yyzz_0,  \
                             ta_xzz_yyzz_1,  \
                             ta_xzz_yzz_0,   \
                             ta_xzz_yzz_1,   \
                             ta_xzz_yzzz_0,  \
                             ta_xzz_yzzz_1,  \
                             ta_xzz_zzz_0,   \
                             ta_xzz_zzz_1,   \
                             ta_xzz_zzzz_0,  \
                             ta_xzz_zzzz_1,  \
                             ta_zz_xxxz_0,   \
                             ta_zz_xxxz_1,   \
                             ta_zz_xxyz_0,   \
                             ta_zz_xxyz_1,   \
                             ta_zz_xxzz_0,   \
                             ta_zz_xxzz_1,   \
                             ta_zz_xyyz_0,   \
                             ta_zz_xyyz_1,   \
                             ta_zz_xyzz_0,   \
                             ta_zz_xyzz_1,   \
                             ta_zz_xzzz_0,   \
                             ta_zz_xzzz_1,   \
                             ta_zz_yyyy_0,   \
                             ta_zz_yyyy_1,   \
                             ta_zz_yyyz_0,   \
                             ta_zz_yyyz_1,   \
                             ta_zz_yyzz_0,   \
                             ta_zz_yyzz_1,   \
                             ta_zz_yzzz_0,   \
                             ta_zz_yzzz_1,   \
                             ta_zz_zzzz_0,   \
                             ta_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_xxxx_0[i] = ta_xx_xxxx_0[i] * fe_0 - ta_xx_xxxx_1[i] * fe_0 + ta_xxz_xxxx_0[i] * pa_z[i] - ta_xxz_xxxx_1[i] * pc_z[i];

        ta_xxzz_xxxy_0[i] = ta_xx_xxxy_0[i] * fe_0 - ta_xx_xxxy_1[i] * fe_0 + ta_xxz_xxxy_0[i] * pa_z[i] - ta_xxz_xxxy_1[i] * pc_z[i];

        ta_xxzz_xxxz_0[i] = ta_zz_xxxz_0[i] * fe_0 - ta_zz_xxxz_1[i] * fe_0 + 3.0 * ta_xzz_xxz_0[i] * fe_0 - 3.0 * ta_xzz_xxz_1[i] * fe_0 +
                            ta_xzz_xxxz_0[i] * pa_x[i] - ta_xzz_xxxz_1[i] * pc_x[i];

        ta_xxzz_xxyy_0[i] = ta_xx_xxyy_0[i] * fe_0 - ta_xx_xxyy_1[i] * fe_0 + ta_xxz_xxyy_0[i] * pa_z[i] - ta_xxz_xxyy_1[i] * pc_z[i];

        ta_xxzz_xxyz_0[i] = ta_zz_xxyz_0[i] * fe_0 - ta_zz_xxyz_1[i] * fe_0 + 2.0 * ta_xzz_xyz_0[i] * fe_0 - 2.0 * ta_xzz_xyz_1[i] * fe_0 +
                            ta_xzz_xxyz_0[i] * pa_x[i] - ta_xzz_xxyz_1[i] * pc_x[i];

        ta_xxzz_xxzz_0[i] = ta_zz_xxzz_0[i] * fe_0 - ta_zz_xxzz_1[i] * fe_0 + 2.0 * ta_xzz_xzz_0[i] * fe_0 - 2.0 * ta_xzz_xzz_1[i] * fe_0 +
                            ta_xzz_xxzz_0[i] * pa_x[i] - ta_xzz_xxzz_1[i] * pc_x[i];

        ta_xxzz_xyyy_0[i] = ta_xx_xyyy_0[i] * fe_0 - ta_xx_xyyy_1[i] * fe_0 + ta_xxz_xyyy_0[i] * pa_z[i] - ta_xxz_xyyy_1[i] * pc_z[i];

        ta_xxzz_xyyz_0[i] = ta_zz_xyyz_0[i] * fe_0 - ta_zz_xyyz_1[i] * fe_0 + ta_xzz_yyz_0[i] * fe_0 - ta_xzz_yyz_1[i] * fe_0 +
                            ta_xzz_xyyz_0[i] * pa_x[i] - ta_xzz_xyyz_1[i] * pc_x[i];

        ta_xxzz_xyzz_0[i] = ta_zz_xyzz_0[i] * fe_0 - ta_zz_xyzz_1[i] * fe_0 + ta_xzz_yzz_0[i] * fe_0 - ta_xzz_yzz_1[i] * fe_0 +
                            ta_xzz_xyzz_0[i] * pa_x[i] - ta_xzz_xyzz_1[i] * pc_x[i];

        ta_xxzz_xzzz_0[i] = ta_zz_xzzz_0[i] * fe_0 - ta_zz_xzzz_1[i] * fe_0 + ta_xzz_zzz_0[i] * fe_0 - ta_xzz_zzz_1[i] * fe_0 +
                            ta_xzz_xzzz_0[i] * pa_x[i] - ta_xzz_xzzz_1[i] * pc_x[i];

        ta_xxzz_yyyy_0[i] = ta_zz_yyyy_0[i] * fe_0 - ta_zz_yyyy_1[i] * fe_0 + ta_xzz_yyyy_0[i] * pa_x[i] - ta_xzz_yyyy_1[i] * pc_x[i];

        ta_xxzz_yyyz_0[i] = ta_zz_yyyz_0[i] * fe_0 - ta_zz_yyyz_1[i] * fe_0 + ta_xzz_yyyz_0[i] * pa_x[i] - ta_xzz_yyyz_1[i] * pc_x[i];

        ta_xxzz_yyzz_0[i] = ta_zz_yyzz_0[i] * fe_0 - ta_zz_yyzz_1[i] * fe_0 + ta_xzz_yyzz_0[i] * pa_x[i] - ta_xzz_yyzz_1[i] * pc_x[i];

        ta_xxzz_yzzz_0[i] = ta_zz_yzzz_0[i] * fe_0 - ta_zz_yzzz_1[i] * fe_0 + ta_xzz_yzzz_0[i] * pa_x[i] - ta_xzz_yzzz_1[i] * pc_x[i];

        ta_xxzz_zzzz_0[i] = ta_zz_zzzz_0[i] * fe_0 - ta_zz_zzzz_1[i] * fe_0 + ta_xzz_zzzz_0[i] * pa_x[i] - ta_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

    auto ta_xyyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 90);

    auto ta_xyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 91);

    auto ta_xyyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 92);

    auto ta_xyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 93);

    auto ta_xyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 94);

    auto ta_xyyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 95);

    auto ta_xyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 96);

    auto ta_xyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 97);

    auto ta_xyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 98);

    auto ta_xyyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 99);

    auto ta_xyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 100);

    auto ta_xyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 101);

    auto ta_xyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 102);

    auto ta_xyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 103);

    auto ta_xyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 104);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyy_xxxx_0, \
                             ta_xyyy_xxxy_0, \
                             ta_xyyy_xxxz_0, \
                             ta_xyyy_xxyy_0, \
                             ta_xyyy_xxyz_0, \
                             ta_xyyy_xxzz_0, \
                             ta_xyyy_xyyy_0, \
                             ta_xyyy_xyyz_0, \
                             ta_xyyy_xyzz_0, \
                             ta_xyyy_xzzz_0, \
                             ta_xyyy_yyyy_0, \
                             ta_xyyy_yyyz_0, \
                             ta_xyyy_yyzz_0, \
                             ta_xyyy_yzzz_0, \
                             ta_xyyy_zzzz_0, \
                             ta_yyy_xxx_0,   \
                             ta_yyy_xxx_1,   \
                             ta_yyy_xxxx_0,  \
                             ta_yyy_xxxx_1,  \
                             ta_yyy_xxxy_0,  \
                             ta_yyy_xxxy_1,  \
                             ta_yyy_xxxz_0,  \
                             ta_yyy_xxxz_1,  \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xxyy_0,  \
                             ta_yyy_xxyy_1,  \
                             ta_yyy_xxyz_0,  \
                             ta_yyy_xxyz_1,  \
                             ta_yyy_xxz_0,   \
                             ta_yyy_xxz_1,   \
                             ta_yyy_xxzz_0,  \
                             ta_yyy_xxzz_1,  \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_xyyy_0,  \
                             ta_yyy_xyyy_1,  \
                             ta_yyy_xyyz_0,  \
                             ta_yyy_xyyz_1,  \
                             ta_yyy_xyz_0,   \
                             ta_yyy_xyz_1,   \
                             ta_yyy_xyzz_0,  \
                             ta_yyy_xyzz_1,  \
                             ta_yyy_xzz_0,   \
                             ta_yyy_xzz_1,   \
                             ta_yyy_xzzz_0,  \
                             ta_yyy_xzzz_1,  \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyy_yyyy_0,  \
                             ta_yyy_yyyy_1,  \
                             ta_yyy_yyyz_0,  \
                             ta_yyy_yyyz_1,  \
                             ta_yyy_yyz_0,   \
                             ta_yyy_yyz_1,   \
                             ta_yyy_yyzz_0,  \
                             ta_yyy_yyzz_1,  \
                             ta_yyy_yzz_0,   \
                             ta_yyy_yzz_1,   \
                             ta_yyy_yzzz_0,  \
                             ta_yyy_yzzz_1,  \
                             ta_yyy_zzz_0,   \
                             ta_yyy_zzz_1,   \
                             ta_yyy_zzzz_0,  \
                             ta_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_xxxx_0[i] = 4.0 * ta_yyy_xxx_0[i] * fe_0 - 4.0 * ta_yyy_xxx_1[i] * fe_0 + ta_yyy_xxxx_0[i] * pa_x[i] - ta_yyy_xxxx_1[i] * pc_x[i];

        ta_xyyy_xxxy_0[i] = 3.0 * ta_yyy_xxy_0[i] * fe_0 - 3.0 * ta_yyy_xxy_1[i] * fe_0 + ta_yyy_xxxy_0[i] * pa_x[i] - ta_yyy_xxxy_1[i] * pc_x[i];

        ta_xyyy_xxxz_0[i] = 3.0 * ta_yyy_xxz_0[i] * fe_0 - 3.0 * ta_yyy_xxz_1[i] * fe_0 + ta_yyy_xxxz_0[i] * pa_x[i] - ta_yyy_xxxz_1[i] * pc_x[i];

        ta_xyyy_xxyy_0[i] = 2.0 * ta_yyy_xyy_0[i] * fe_0 - 2.0 * ta_yyy_xyy_1[i] * fe_0 + ta_yyy_xxyy_0[i] * pa_x[i] - ta_yyy_xxyy_1[i] * pc_x[i];

        ta_xyyy_xxyz_0[i] = 2.0 * ta_yyy_xyz_0[i] * fe_0 - 2.0 * ta_yyy_xyz_1[i] * fe_0 + ta_yyy_xxyz_0[i] * pa_x[i] - ta_yyy_xxyz_1[i] * pc_x[i];

        ta_xyyy_xxzz_0[i] = 2.0 * ta_yyy_xzz_0[i] * fe_0 - 2.0 * ta_yyy_xzz_1[i] * fe_0 + ta_yyy_xxzz_0[i] * pa_x[i] - ta_yyy_xxzz_1[i] * pc_x[i];

        ta_xyyy_xyyy_0[i] = ta_yyy_yyy_0[i] * fe_0 - ta_yyy_yyy_1[i] * fe_0 + ta_yyy_xyyy_0[i] * pa_x[i] - ta_yyy_xyyy_1[i] * pc_x[i];

        ta_xyyy_xyyz_0[i] = ta_yyy_yyz_0[i] * fe_0 - ta_yyy_yyz_1[i] * fe_0 + ta_yyy_xyyz_0[i] * pa_x[i] - ta_yyy_xyyz_1[i] * pc_x[i];

        ta_xyyy_xyzz_0[i] = ta_yyy_yzz_0[i] * fe_0 - ta_yyy_yzz_1[i] * fe_0 + ta_yyy_xyzz_0[i] * pa_x[i] - ta_yyy_xyzz_1[i] * pc_x[i];

        ta_xyyy_xzzz_0[i] = ta_yyy_zzz_0[i] * fe_0 - ta_yyy_zzz_1[i] * fe_0 + ta_yyy_xzzz_0[i] * pa_x[i] - ta_yyy_xzzz_1[i] * pc_x[i];

        ta_xyyy_yyyy_0[i] = ta_yyy_yyyy_0[i] * pa_x[i] - ta_yyy_yyyy_1[i] * pc_x[i];

        ta_xyyy_yyyz_0[i] = ta_yyy_yyyz_0[i] * pa_x[i] - ta_yyy_yyyz_1[i] * pc_x[i];

        ta_xyyy_yyzz_0[i] = ta_yyy_yyzz_0[i] * pa_x[i] - ta_yyy_yyzz_1[i] * pc_x[i];

        ta_xyyy_yzzz_0[i] = ta_yyy_yzzz_0[i] * pa_x[i] - ta_yyy_yzzz_1[i] * pc_x[i];

        ta_xyyy_zzzz_0[i] = ta_yyy_zzzz_0[i] * pa_x[i] - ta_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto ta_xyyz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 105);

    auto ta_xyyz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 106);

    auto ta_xyyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 107);

    auto ta_xyyz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 108);

    auto ta_xyyz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 109);

    auto ta_xyyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 110);

    auto ta_xyyz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 111);

    auto ta_xyyz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 112);

    auto ta_xyyz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 113);

    auto ta_xyyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 114);

    auto ta_xyyz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 115);

    auto ta_xyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 116);

    auto ta_xyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 117);

    auto ta_xyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 118);

    auto ta_xyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xyy_xxxx_0,  \
                             ta_xyy_xxxx_1,  \
                             ta_xyy_xxxy_0,  \
                             ta_xyy_xxxy_1,  \
                             ta_xyy_xxyy_0,  \
                             ta_xyy_xxyy_1,  \
                             ta_xyy_xyyy_0,  \
                             ta_xyy_xyyy_1,  \
                             ta_xyyz_xxxx_0, \
                             ta_xyyz_xxxy_0, \
                             ta_xyyz_xxxz_0, \
                             ta_xyyz_xxyy_0, \
                             ta_xyyz_xxyz_0, \
                             ta_xyyz_xxzz_0, \
                             ta_xyyz_xyyy_0, \
                             ta_xyyz_xyyz_0, \
                             ta_xyyz_xyzz_0, \
                             ta_xyyz_xzzz_0, \
                             ta_xyyz_yyyy_0, \
                             ta_xyyz_yyyz_0, \
                             ta_xyyz_yyzz_0, \
                             ta_xyyz_yzzz_0, \
                             ta_xyyz_zzzz_0, \
                             ta_yyz_xxxz_0,  \
                             ta_yyz_xxxz_1,  \
                             ta_yyz_xxyz_0,  \
                             ta_yyz_xxyz_1,  \
                             ta_yyz_xxz_0,   \
                             ta_yyz_xxz_1,   \
                             ta_yyz_xxzz_0,  \
                             ta_yyz_xxzz_1,  \
                             ta_yyz_xyyz_0,  \
                             ta_yyz_xyyz_1,  \
                             ta_yyz_xyz_0,   \
                             ta_yyz_xyz_1,   \
                             ta_yyz_xyzz_0,  \
                             ta_yyz_xyzz_1,  \
                             ta_yyz_xzz_0,   \
                             ta_yyz_xzz_1,   \
                             ta_yyz_xzzz_0,  \
                             ta_yyz_xzzz_1,  \
                             ta_yyz_yyyy_0,  \
                             ta_yyz_yyyy_1,  \
                             ta_yyz_yyyz_0,  \
                             ta_yyz_yyyz_1,  \
                             ta_yyz_yyz_0,   \
                             ta_yyz_yyz_1,   \
                             ta_yyz_yyzz_0,  \
                             ta_yyz_yyzz_1,  \
                             ta_yyz_yzz_0,   \
                             ta_yyz_yzz_1,   \
                             ta_yyz_yzzz_0,  \
                             ta_yyz_yzzz_1,  \
                             ta_yyz_zzz_0,   \
                             ta_yyz_zzz_1,   \
                             ta_yyz_zzzz_0,  \
                             ta_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyz_xxxx_0[i] = ta_xyy_xxxx_0[i] * pa_z[i] - ta_xyy_xxxx_1[i] * pc_z[i];

        ta_xyyz_xxxy_0[i] = ta_xyy_xxxy_0[i] * pa_z[i] - ta_xyy_xxxy_1[i] * pc_z[i];

        ta_xyyz_xxxz_0[i] = 3.0 * ta_yyz_xxz_0[i] * fe_0 - 3.0 * ta_yyz_xxz_1[i] * fe_0 + ta_yyz_xxxz_0[i] * pa_x[i] - ta_yyz_xxxz_1[i] * pc_x[i];

        ta_xyyz_xxyy_0[i] = ta_xyy_xxyy_0[i] * pa_z[i] - ta_xyy_xxyy_1[i] * pc_z[i];

        ta_xyyz_xxyz_0[i] = 2.0 * ta_yyz_xyz_0[i] * fe_0 - 2.0 * ta_yyz_xyz_1[i] * fe_0 + ta_yyz_xxyz_0[i] * pa_x[i] - ta_yyz_xxyz_1[i] * pc_x[i];

        ta_xyyz_xxzz_0[i] = 2.0 * ta_yyz_xzz_0[i] * fe_0 - 2.0 * ta_yyz_xzz_1[i] * fe_0 + ta_yyz_xxzz_0[i] * pa_x[i] - ta_yyz_xxzz_1[i] * pc_x[i];

        ta_xyyz_xyyy_0[i] = ta_xyy_xyyy_0[i] * pa_z[i] - ta_xyy_xyyy_1[i] * pc_z[i];

        ta_xyyz_xyyz_0[i] = ta_yyz_yyz_0[i] * fe_0 - ta_yyz_yyz_1[i] * fe_0 + ta_yyz_xyyz_0[i] * pa_x[i] - ta_yyz_xyyz_1[i] * pc_x[i];

        ta_xyyz_xyzz_0[i] = ta_yyz_yzz_0[i] * fe_0 - ta_yyz_yzz_1[i] * fe_0 + ta_yyz_xyzz_0[i] * pa_x[i] - ta_yyz_xyzz_1[i] * pc_x[i];

        ta_xyyz_xzzz_0[i] = ta_yyz_zzz_0[i] * fe_0 - ta_yyz_zzz_1[i] * fe_0 + ta_yyz_xzzz_0[i] * pa_x[i] - ta_yyz_xzzz_1[i] * pc_x[i];

        ta_xyyz_yyyy_0[i] = ta_yyz_yyyy_0[i] * pa_x[i] - ta_yyz_yyyy_1[i] * pc_x[i];

        ta_xyyz_yyyz_0[i] = ta_yyz_yyyz_0[i] * pa_x[i] - ta_yyz_yyyz_1[i] * pc_x[i];

        ta_xyyz_yyzz_0[i] = ta_yyz_yyzz_0[i] * pa_x[i] - ta_yyz_yyzz_1[i] * pc_x[i];

        ta_xyyz_yzzz_0[i] = ta_yyz_yzzz_0[i] * pa_x[i] - ta_yyz_yzzz_1[i] * pc_x[i];

        ta_xyyz_zzzz_0[i] = ta_yyz_zzzz_0[i] * pa_x[i] - ta_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto ta_xyzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 120);

    auto ta_xyzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 121);

    auto ta_xyzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 122);

    auto ta_xyzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 123);

    auto ta_xyzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 124);

    auto ta_xyzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 125);

    auto ta_xyzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 126);

    auto ta_xyzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 127);

    auto ta_xyzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 128);

    auto ta_xyzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 129);

    auto ta_xyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 130);

    auto ta_xyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 131);

    auto ta_xyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 132);

    auto ta_xyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 133);

    auto ta_xyzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 134);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xyzz_xxxx_0, \
                             ta_xyzz_xxxy_0, \
                             ta_xyzz_xxxz_0, \
                             ta_xyzz_xxyy_0, \
                             ta_xyzz_xxyz_0, \
                             ta_xyzz_xxzz_0, \
                             ta_xyzz_xyyy_0, \
                             ta_xyzz_xyyz_0, \
                             ta_xyzz_xyzz_0, \
                             ta_xyzz_xzzz_0, \
                             ta_xyzz_yyyy_0, \
                             ta_xyzz_yyyz_0, \
                             ta_xyzz_yyzz_0, \
                             ta_xyzz_yzzz_0, \
                             ta_xyzz_zzzz_0, \
                             ta_xzz_xxxx_0,  \
                             ta_xzz_xxxx_1,  \
                             ta_xzz_xxxz_0,  \
                             ta_xzz_xxxz_1,  \
                             ta_xzz_xxzz_0,  \
                             ta_xzz_xxzz_1,  \
                             ta_xzz_xzzz_0,  \
                             ta_xzz_xzzz_1,  \
                             ta_yzz_xxxy_0,  \
                             ta_yzz_xxxy_1,  \
                             ta_yzz_xxy_0,   \
                             ta_yzz_xxy_1,   \
                             ta_yzz_xxyy_0,  \
                             ta_yzz_xxyy_1,  \
                             ta_yzz_xxyz_0,  \
                             ta_yzz_xxyz_1,  \
                             ta_yzz_xyy_0,   \
                             ta_yzz_xyy_1,   \
                             ta_yzz_xyyy_0,  \
                             ta_yzz_xyyy_1,  \
                             ta_yzz_xyyz_0,  \
                             ta_yzz_xyyz_1,  \
                             ta_yzz_xyz_0,   \
                             ta_yzz_xyz_1,   \
                             ta_yzz_xyzz_0,  \
                             ta_yzz_xyzz_1,  \
                             ta_yzz_yyy_0,   \
                             ta_yzz_yyy_1,   \
                             ta_yzz_yyyy_0,  \
                             ta_yzz_yyyy_1,  \
                             ta_yzz_yyyz_0,  \
                             ta_yzz_yyyz_1,  \
                             ta_yzz_yyz_0,   \
                             ta_yzz_yyz_1,   \
                             ta_yzz_yyzz_0,  \
                             ta_yzz_yyzz_1,  \
                             ta_yzz_yzz_0,   \
                             ta_yzz_yzz_1,   \
                             ta_yzz_yzzz_0,  \
                             ta_yzz_yzzz_1,  \
                             ta_yzz_zzzz_0,  \
                             ta_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzz_xxxx_0[i] = ta_xzz_xxxx_0[i] * pa_y[i] - ta_xzz_xxxx_1[i] * pc_y[i];

        ta_xyzz_xxxy_0[i] = 3.0 * ta_yzz_xxy_0[i] * fe_0 - 3.0 * ta_yzz_xxy_1[i] * fe_0 + ta_yzz_xxxy_0[i] * pa_x[i] - ta_yzz_xxxy_1[i] * pc_x[i];

        ta_xyzz_xxxz_0[i] = ta_xzz_xxxz_0[i] * pa_y[i] - ta_xzz_xxxz_1[i] * pc_y[i];

        ta_xyzz_xxyy_0[i] = 2.0 * ta_yzz_xyy_0[i] * fe_0 - 2.0 * ta_yzz_xyy_1[i] * fe_0 + ta_yzz_xxyy_0[i] * pa_x[i] - ta_yzz_xxyy_1[i] * pc_x[i];

        ta_xyzz_xxyz_0[i] = 2.0 * ta_yzz_xyz_0[i] * fe_0 - 2.0 * ta_yzz_xyz_1[i] * fe_0 + ta_yzz_xxyz_0[i] * pa_x[i] - ta_yzz_xxyz_1[i] * pc_x[i];

        ta_xyzz_xxzz_0[i] = ta_xzz_xxzz_0[i] * pa_y[i] - ta_xzz_xxzz_1[i] * pc_y[i];

        ta_xyzz_xyyy_0[i] = ta_yzz_yyy_0[i] * fe_0 - ta_yzz_yyy_1[i] * fe_0 + ta_yzz_xyyy_0[i] * pa_x[i] - ta_yzz_xyyy_1[i] * pc_x[i];

        ta_xyzz_xyyz_0[i] = ta_yzz_yyz_0[i] * fe_0 - ta_yzz_yyz_1[i] * fe_0 + ta_yzz_xyyz_0[i] * pa_x[i] - ta_yzz_xyyz_1[i] * pc_x[i];

        ta_xyzz_xyzz_0[i] = ta_yzz_yzz_0[i] * fe_0 - ta_yzz_yzz_1[i] * fe_0 + ta_yzz_xyzz_0[i] * pa_x[i] - ta_yzz_xyzz_1[i] * pc_x[i];

        ta_xyzz_xzzz_0[i] = ta_xzz_xzzz_0[i] * pa_y[i] - ta_xzz_xzzz_1[i] * pc_y[i];

        ta_xyzz_yyyy_0[i] = ta_yzz_yyyy_0[i] * pa_x[i] - ta_yzz_yyyy_1[i] * pc_x[i];

        ta_xyzz_yyyz_0[i] = ta_yzz_yyyz_0[i] * pa_x[i] - ta_yzz_yyyz_1[i] * pc_x[i];

        ta_xyzz_yyzz_0[i] = ta_yzz_yyzz_0[i] * pa_x[i] - ta_yzz_yyzz_1[i] * pc_x[i];

        ta_xyzz_yzzz_0[i] = ta_yzz_yzzz_0[i] * pa_x[i] - ta_yzz_yzzz_1[i] * pc_x[i];

        ta_xyzz_zzzz_0[i] = ta_yzz_zzzz_0[i] * pa_x[i] - ta_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 135-150 components of targeted buffer : GG

    auto ta_xzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 135);

    auto ta_xzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 136);

    auto ta_xzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 137);

    auto ta_xzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 138);

    auto ta_xzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 139);

    auto ta_xzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 140);

    auto ta_xzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 141);

    auto ta_xzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 142);

    auto ta_xzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 143);

    auto ta_xzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 144);

    auto ta_xzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 145);

    auto ta_xzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 146);

    auto ta_xzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 147);

    auto ta_xzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 148);

    auto ta_xzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 149);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xzzz_xxxx_0, \
                             ta_xzzz_xxxy_0, \
                             ta_xzzz_xxxz_0, \
                             ta_xzzz_xxyy_0, \
                             ta_xzzz_xxyz_0, \
                             ta_xzzz_xxzz_0, \
                             ta_xzzz_xyyy_0, \
                             ta_xzzz_xyyz_0, \
                             ta_xzzz_xyzz_0, \
                             ta_xzzz_xzzz_0, \
                             ta_xzzz_yyyy_0, \
                             ta_xzzz_yyyz_0, \
                             ta_xzzz_yyzz_0, \
                             ta_xzzz_yzzz_0, \
                             ta_xzzz_zzzz_0, \
                             ta_zzz_xxx_0,   \
                             ta_zzz_xxx_1,   \
                             ta_zzz_xxxx_0,  \
                             ta_zzz_xxxx_1,  \
                             ta_zzz_xxxy_0,  \
                             ta_zzz_xxxy_1,  \
                             ta_zzz_xxxz_0,  \
                             ta_zzz_xxxz_1,  \
                             ta_zzz_xxy_0,   \
                             ta_zzz_xxy_1,   \
                             ta_zzz_xxyy_0,  \
                             ta_zzz_xxyy_1,  \
                             ta_zzz_xxyz_0,  \
                             ta_zzz_xxyz_1,  \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xxzz_0,  \
                             ta_zzz_xxzz_1,  \
                             ta_zzz_xyy_0,   \
                             ta_zzz_xyy_1,   \
                             ta_zzz_xyyy_0,  \
                             ta_zzz_xyyy_1,  \
                             ta_zzz_xyyz_0,  \
                             ta_zzz_xyyz_1,  \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xyzz_0,  \
                             ta_zzz_xyzz_1,  \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_xzzz_0,  \
                             ta_zzz_xzzz_1,  \
                             ta_zzz_yyy_0,   \
                             ta_zzz_yyy_1,   \
                             ta_zzz_yyyy_0,  \
                             ta_zzz_yyyy_1,  \
                             ta_zzz_yyyz_0,  \
                             ta_zzz_yyyz_1,  \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yyzz_0,  \
                             ta_zzz_yyzz_1,  \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_yzzz_0,  \
                             ta_zzz_yzzz_1,  \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             ta_zzz_zzzz_0,  \
                             ta_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_xxxx_0[i] = 4.0 * ta_zzz_xxx_0[i] * fe_0 - 4.0 * ta_zzz_xxx_1[i] * fe_0 + ta_zzz_xxxx_0[i] * pa_x[i] - ta_zzz_xxxx_1[i] * pc_x[i];

        ta_xzzz_xxxy_0[i] = 3.0 * ta_zzz_xxy_0[i] * fe_0 - 3.0 * ta_zzz_xxy_1[i] * fe_0 + ta_zzz_xxxy_0[i] * pa_x[i] - ta_zzz_xxxy_1[i] * pc_x[i];

        ta_xzzz_xxxz_0[i] = 3.0 * ta_zzz_xxz_0[i] * fe_0 - 3.0 * ta_zzz_xxz_1[i] * fe_0 + ta_zzz_xxxz_0[i] * pa_x[i] - ta_zzz_xxxz_1[i] * pc_x[i];

        ta_xzzz_xxyy_0[i] = 2.0 * ta_zzz_xyy_0[i] * fe_0 - 2.0 * ta_zzz_xyy_1[i] * fe_0 + ta_zzz_xxyy_0[i] * pa_x[i] - ta_zzz_xxyy_1[i] * pc_x[i];

        ta_xzzz_xxyz_0[i] = 2.0 * ta_zzz_xyz_0[i] * fe_0 - 2.0 * ta_zzz_xyz_1[i] * fe_0 + ta_zzz_xxyz_0[i] * pa_x[i] - ta_zzz_xxyz_1[i] * pc_x[i];

        ta_xzzz_xxzz_0[i] = 2.0 * ta_zzz_xzz_0[i] * fe_0 - 2.0 * ta_zzz_xzz_1[i] * fe_0 + ta_zzz_xxzz_0[i] * pa_x[i] - ta_zzz_xxzz_1[i] * pc_x[i];

        ta_xzzz_xyyy_0[i] = ta_zzz_yyy_0[i] * fe_0 - ta_zzz_yyy_1[i] * fe_0 + ta_zzz_xyyy_0[i] * pa_x[i] - ta_zzz_xyyy_1[i] * pc_x[i];

        ta_xzzz_xyyz_0[i] = ta_zzz_yyz_0[i] * fe_0 - ta_zzz_yyz_1[i] * fe_0 + ta_zzz_xyyz_0[i] * pa_x[i] - ta_zzz_xyyz_1[i] * pc_x[i];

        ta_xzzz_xyzz_0[i] = ta_zzz_yzz_0[i] * fe_0 - ta_zzz_yzz_1[i] * fe_0 + ta_zzz_xyzz_0[i] * pa_x[i] - ta_zzz_xyzz_1[i] * pc_x[i];

        ta_xzzz_xzzz_0[i] = ta_zzz_zzz_0[i] * fe_0 - ta_zzz_zzz_1[i] * fe_0 + ta_zzz_xzzz_0[i] * pa_x[i] - ta_zzz_xzzz_1[i] * pc_x[i];

        ta_xzzz_yyyy_0[i] = ta_zzz_yyyy_0[i] * pa_x[i] - ta_zzz_yyyy_1[i] * pc_x[i];

        ta_xzzz_yyyz_0[i] = ta_zzz_yyyz_0[i] * pa_x[i] - ta_zzz_yyyz_1[i] * pc_x[i];

        ta_xzzz_yyzz_0[i] = ta_zzz_yyzz_0[i] * pa_x[i] - ta_zzz_yyzz_1[i] * pc_x[i];

        ta_xzzz_yzzz_0[i] = ta_zzz_yzzz_0[i] * pa_x[i] - ta_zzz_yzzz_1[i] * pc_x[i];

        ta_xzzz_zzzz_0[i] = ta_zzz_zzzz_0[i] * pa_x[i] - ta_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

    auto ta_yyyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 150);

    auto ta_yyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 151);

    auto ta_yyyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 152);

    auto ta_yyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 153);

    auto ta_yyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 154);

    auto ta_yyyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 155);

    auto ta_yyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 156);

    auto ta_yyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 157);

    auto ta_yyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 158);

    auto ta_yyyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 159);

    auto ta_yyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 160);

    auto ta_yyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 161);

    auto ta_yyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 162);

    auto ta_yyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 163);

    auto ta_yyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 164);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yy_xxxx_0,   \
                             ta_yy_xxxx_1,   \
                             ta_yy_xxxy_0,   \
                             ta_yy_xxxy_1,   \
                             ta_yy_xxxz_0,   \
                             ta_yy_xxxz_1,   \
                             ta_yy_xxyy_0,   \
                             ta_yy_xxyy_1,   \
                             ta_yy_xxyz_0,   \
                             ta_yy_xxyz_1,   \
                             ta_yy_xxzz_0,   \
                             ta_yy_xxzz_1,   \
                             ta_yy_xyyy_0,   \
                             ta_yy_xyyy_1,   \
                             ta_yy_xyyz_0,   \
                             ta_yy_xyyz_1,   \
                             ta_yy_xyzz_0,   \
                             ta_yy_xyzz_1,   \
                             ta_yy_xzzz_0,   \
                             ta_yy_xzzz_1,   \
                             ta_yy_yyyy_0,   \
                             ta_yy_yyyy_1,   \
                             ta_yy_yyyz_0,   \
                             ta_yy_yyyz_1,   \
                             ta_yy_yyzz_0,   \
                             ta_yy_yyzz_1,   \
                             ta_yy_yzzz_0,   \
                             ta_yy_yzzz_1,   \
                             ta_yy_zzzz_0,   \
                             ta_yy_zzzz_1,   \
                             ta_yyy_xxx_0,   \
                             ta_yyy_xxx_1,   \
                             ta_yyy_xxxx_0,  \
                             ta_yyy_xxxx_1,  \
                             ta_yyy_xxxy_0,  \
                             ta_yyy_xxxy_1,  \
                             ta_yyy_xxxz_0,  \
                             ta_yyy_xxxz_1,  \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xxyy_0,  \
                             ta_yyy_xxyy_1,  \
                             ta_yyy_xxyz_0,  \
                             ta_yyy_xxyz_1,  \
                             ta_yyy_xxz_0,   \
                             ta_yyy_xxz_1,   \
                             ta_yyy_xxzz_0,  \
                             ta_yyy_xxzz_1,  \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_xyyy_0,  \
                             ta_yyy_xyyy_1,  \
                             ta_yyy_xyyz_0,  \
                             ta_yyy_xyyz_1,  \
                             ta_yyy_xyz_0,   \
                             ta_yyy_xyz_1,   \
                             ta_yyy_xyzz_0,  \
                             ta_yyy_xyzz_1,  \
                             ta_yyy_xzz_0,   \
                             ta_yyy_xzz_1,   \
                             ta_yyy_xzzz_0,  \
                             ta_yyy_xzzz_1,  \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyy_yyyy_0,  \
                             ta_yyy_yyyy_1,  \
                             ta_yyy_yyyz_0,  \
                             ta_yyy_yyyz_1,  \
                             ta_yyy_yyz_0,   \
                             ta_yyy_yyz_1,   \
                             ta_yyy_yyzz_0,  \
                             ta_yyy_yyzz_1,  \
                             ta_yyy_yzz_0,   \
                             ta_yyy_yzz_1,   \
                             ta_yyy_yzzz_0,  \
                             ta_yyy_yzzz_1,  \
                             ta_yyy_zzz_0,   \
                             ta_yyy_zzz_1,   \
                             ta_yyy_zzzz_0,  \
                             ta_yyy_zzzz_1,  \
                             ta_yyyy_xxxx_0, \
                             ta_yyyy_xxxy_0, \
                             ta_yyyy_xxxz_0, \
                             ta_yyyy_xxyy_0, \
                             ta_yyyy_xxyz_0, \
                             ta_yyyy_xxzz_0, \
                             ta_yyyy_xyyy_0, \
                             ta_yyyy_xyyz_0, \
                             ta_yyyy_xyzz_0, \
                             ta_yyyy_xzzz_0, \
                             ta_yyyy_yyyy_0, \
                             ta_yyyy_yyyz_0, \
                             ta_yyyy_yyzz_0, \
                             ta_yyyy_yzzz_0, \
                             ta_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_xxxx_0[i] = 3.0 * ta_yy_xxxx_0[i] * fe_0 - 3.0 * ta_yy_xxxx_1[i] * fe_0 + ta_yyy_xxxx_0[i] * pa_y[i] - ta_yyy_xxxx_1[i] * pc_y[i];

        ta_yyyy_xxxy_0[i] = 3.0 * ta_yy_xxxy_0[i] * fe_0 - 3.0 * ta_yy_xxxy_1[i] * fe_0 + ta_yyy_xxx_0[i] * fe_0 - ta_yyy_xxx_1[i] * fe_0 +
                            ta_yyy_xxxy_0[i] * pa_y[i] - ta_yyy_xxxy_1[i] * pc_y[i];

        ta_yyyy_xxxz_0[i] = 3.0 * ta_yy_xxxz_0[i] * fe_0 - 3.0 * ta_yy_xxxz_1[i] * fe_0 + ta_yyy_xxxz_0[i] * pa_y[i] - ta_yyy_xxxz_1[i] * pc_y[i];

        ta_yyyy_xxyy_0[i] = 3.0 * ta_yy_xxyy_0[i] * fe_0 - 3.0 * ta_yy_xxyy_1[i] * fe_0 + 2.0 * ta_yyy_xxy_0[i] * fe_0 -
                            2.0 * ta_yyy_xxy_1[i] * fe_0 + ta_yyy_xxyy_0[i] * pa_y[i] - ta_yyy_xxyy_1[i] * pc_y[i];

        ta_yyyy_xxyz_0[i] = 3.0 * ta_yy_xxyz_0[i] * fe_0 - 3.0 * ta_yy_xxyz_1[i] * fe_0 + ta_yyy_xxz_0[i] * fe_0 - ta_yyy_xxz_1[i] * fe_0 +
                            ta_yyy_xxyz_0[i] * pa_y[i] - ta_yyy_xxyz_1[i] * pc_y[i];

        ta_yyyy_xxzz_0[i] = 3.0 * ta_yy_xxzz_0[i] * fe_0 - 3.0 * ta_yy_xxzz_1[i] * fe_0 + ta_yyy_xxzz_0[i] * pa_y[i] - ta_yyy_xxzz_1[i] * pc_y[i];

        ta_yyyy_xyyy_0[i] = 3.0 * ta_yy_xyyy_0[i] * fe_0 - 3.0 * ta_yy_xyyy_1[i] * fe_0 + 3.0 * ta_yyy_xyy_0[i] * fe_0 -
                            3.0 * ta_yyy_xyy_1[i] * fe_0 + ta_yyy_xyyy_0[i] * pa_y[i] - ta_yyy_xyyy_1[i] * pc_y[i];

        ta_yyyy_xyyz_0[i] = 3.0 * ta_yy_xyyz_0[i] * fe_0 - 3.0 * ta_yy_xyyz_1[i] * fe_0 + 2.0 * ta_yyy_xyz_0[i] * fe_0 -
                            2.0 * ta_yyy_xyz_1[i] * fe_0 + ta_yyy_xyyz_0[i] * pa_y[i] - ta_yyy_xyyz_1[i] * pc_y[i];

        ta_yyyy_xyzz_0[i] = 3.0 * ta_yy_xyzz_0[i] * fe_0 - 3.0 * ta_yy_xyzz_1[i] * fe_0 + ta_yyy_xzz_0[i] * fe_0 - ta_yyy_xzz_1[i] * fe_0 +
                            ta_yyy_xyzz_0[i] * pa_y[i] - ta_yyy_xyzz_1[i] * pc_y[i];

        ta_yyyy_xzzz_0[i] = 3.0 * ta_yy_xzzz_0[i] * fe_0 - 3.0 * ta_yy_xzzz_1[i] * fe_0 + ta_yyy_xzzz_0[i] * pa_y[i] - ta_yyy_xzzz_1[i] * pc_y[i];

        ta_yyyy_yyyy_0[i] = 3.0 * ta_yy_yyyy_0[i] * fe_0 - 3.0 * ta_yy_yyyy_1[i] * fe_0 + 4.0 * ta_yyy_yyy_0[i] * fe_0 -
                            4.0 * ta_yyy_yyy_1[i] * fe_0 + ta_yyy_yyyy_0[i] * pa_y[i] - ta_yyy_yyyy_1[i] * pc_y[i];

        ta_yyyy_yyyz_0[i] = 3.0 * ta_yy_yyyz_0[i] * fe_0 - 3.0 * ta_yy_yyyz_1[i] * fe_0 + 3.0 * ta_yyy_yyz_0[i] * fe_0 -
                            3.0 * ta_yyy_yyz_1[i] * fe_0 + ta_yyy_yyyz_0[i] * pa_y[i] - ta_yyy_yyyz_1[i] * pc_y[i];

        ta_yyyy_yyzz_0[i] = 3.0 * ta_yy_yyzz_0[i] * fe_0 - 3.0 * ta_yy_yyzz_1[i] * fe_0 + 2.0 * ta_yyy_yzz_0[i] * fe_0 -
                            2.0 * ta_yyy_yzz_1[i] * fe_0 + ta_yyy_yyzz_0[i] * pa_y[i] - ta_yyy_yyzz_1[i] * pc_y[i];

        ta_yyyy_yzzz_0[i] = 3.0 * ta_yy_yzzz_0[i] * fe_0 - 3.0 * ta_yy_yzzz_1[i] * fe_0 + ta_yyy_zzz_0[i] * fe_0 - ta_yyy_zzz_1[i] * fe_0 +
                            ta_yyy_yzzz_0[i] * pa_y[i] - ta_yyy_yzzz_1[i] * pc_y[i];

        ta_yyyy_zzzz_0[i] = 3.0 * ta_yy_zzzz_0[i] * fe_0 - 3.0 * ta_yy_zzzz_1[i] * fe_0 + ta_yyy_zzzz_0[i] * pa_y[i] - ta_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto ta_yyyz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 165);

    auto ta_yyyz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 166);

    auto ta_yyyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 167);

    auto ta_yyyz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 168);

    auto ta_yyyz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 169);

    auto ta_yyyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 170);

    auto ta_yyyz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 171);

    auto ta_yyyz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 172);

    auto ta_yyyz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 173);

    auto ta_yyyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 174);

    auto ta_yyyz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 175);

    auto ta_yyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 176);

    auto ta_yyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 177);

    auto ta_yyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 178);

    auto ta_yyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 179);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyy_xxxx_0,  \
                             ta_yyy_xxxx_1,  \
                             ta_yyy_xxxy_0,  \
                             ta_yyy_xxxy_1,  \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xxyy_0,  \
                             ta_yyy_xxyy_1,  \
                             ta_yyy_xxyz_0,  \
                             ta_yyy_xxyz_1,  \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_xyyy_0,  \
                             ta_yyy_xyyy_1,  \
                             ta_yyy_xyyz_0,  \
                             ta_yyy_xyyz_1,  \
                             ta_yyy_xyz_0,   \
                             ta_yyy_xyz_1,   \
                             ta_yyy_xyzz_0,  \
                             ta_yyy_xyzz_1,  \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyy_yyyy_0,  \
                             ta_yyy_yyyy_1,  \
                             ta_yyy_yyyz_0,  \
                             ta_yyy_yyyz_1,  \
                             ta_yyy_yyz_0,   \
                             ta_yyy_yyz_1,   \
                             ta_yyy_yyzz_0,  \
                             ta_yyy_yyzz_1,  \
                             ta_yyy_yzz_0,   \
                             ta_yyy_yzz_1,   \
                             ta_yyy_yzzz_0,  \
                             ta_yyy_yzzz_1,  \
                             ta_yyyz_xxxx_0, \
                             ta_yyyz_xxxy_0, \
                             ta_yyyz_xxxz_0, \
                             ta_yyyz_xxyy_0, \
                             ta_yyyz_xxyz_0, \
                             ta_yyyz_xxzz_0, \
                             ta_yyyz_xyyy_0, \
                             ta_yyyz_xyyz_0, \
                             ta_yyyz_xyzz_0, \
                             ta_yyyz_xzzz_0, \
                             ta_yyyz_yyyy_0, \
                             ta_yyyz_yyyz_0, \
                             ta_yyyz_yyzz_0, \
                             ta_yyyz_yzzz_0, \
                             ta_yyyz_zzzz_0, \
                             ta_yyz_xxxz_0,  \
                             ta_yyz_xxxz_1,  \
                             ta_yyz_xxzz_0,  \
                             ta_yyz_xxzz_1,  \
                             ta_yyz_xzzz_0,  \
                             ta_yyz_xzzz_1,  \
                             ta_yyz_zzzz_0,  \
                             ta_yyz_zzzz_1,  \
                             ta_yz_xxxz_0,   \
                             ta_yz_xxxz_1,   \
                             ta_yz_xxzz_0,   \
                             ta_yz_xxzz_1,   \
                             ta_yz_xzzz_0,   \
                             ta_yz_xzzz_1,   \
                             ta_yz_zzzz_0,   \
                             ta_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_xxxx_0[i] = ta_yyy_xxxx_0[i] * pa_z[i] - ta_yyy_xxxx_1[i] * pc_z[i];

        ta_yyyz_xxxy_0[i] = ta_yyy_xxxy_0[i] * pa_z[i] - ta_yyy_xxxy_1[i] * pc_z[i];

        ta_yyyz_xxxz_0[i] = 2.0 * ta_yz_xxxz_0[i] * fe_0 - 2.0 * ta_yz_xxxz_1[i] * fe_0 + ta_yyz_xxxz_0[i] * pa_y[i] - ta_yyz_xxxz_1[i] * pc_y[i];

        ta_yyyz_xxyy_0[i] = ta_yyy_xxyy_0[i] * pa_z[i] - ta_yyy_xxyy_1[i] * pc_z[i];

        ta_yyyz_xxyz_0[i] = ta_yyy_xxy_0[i] * fe_0 - ta_yyy_xxy_1[i] * fe_0 + ta_yyy_xxyz_0[i] * pa_z[i] - ta_yyy_xxyz_1[i] * pc_z[i];

        ta_yyyz_xxzz_0[i] = 2.0 * ta_yz_xxzz_0[i] * fe_0 - 2.0 * ta_yz_xxzz_1[i] * fe_0 + ta_yyz_xxzz_0[i] * pa_y[i] - ta_yyz_xxzz_1[i] * pc_y[i];

        ta_yyyz_xyyy_0[i] = ta_yyy_xyyy_0[i] * pa_z[i] - ta_yyy_xyyy_1[i] * pc_z[i];

        ta_yyyz_xyyz_0[i] = ta_yyy_xyy_0[i] * fe_0 - ta_yyy_xyy_1[i] * fe_0 + ta_yyy_xyyz_0[i] * pa_z[i] - ta_yyy_xyyz_1[i] * pc_z[i];

        ta_yyyz_xyzz_0[i] = 2.0 * ta_yyy_xyz_0[i] * fe_0 - 2.0 * ta_yyy_xyz_1[i] * fe_0 + ta_yyy_xyzz_0[i] * pa_z[i] - ta_yyy_xyzz_1[i] * pc_z[i];

        ta_yyyz_xzzz_0[i] = 2.0 * ta_yz_xzzz_0[i] * fe_0 - 2.0 * ta_yz_xzzz_1[i] * fe_0 + ta_yyz_xzzz_0[i] * pa_y[i] - ta_yyz_xzzz_1[i] * pc_y[i];

        ta_yyyz_yyyy_0[i] = ta_yyy_yyyy_0[i] * pa_z[i] - ta_yyy_yyyy_1[i] * pc_z[i];

        ta_yyyz_yyyz_0[i] = ta_yyy_yyy_0[i] * fe_0 - ta_yyy_yyy_1[i] * fe_0 + ta_yyy_yyyz_0[i] * pa_z[i] - ta_yyy_yyyz_1[i] * pc_z[i];

        ta_yyyz_yyzz_0[i] = 2.0 * ta_yyy_yyz_0[i] * fe_0 - 2.0 * ta_yyy_yyz_1[i] * fe_0 + ta_yyy_yyzz_0[i] * pa_z[i] - ta_yyy_yyzz_1[i] * pc_z[i];

        ta_yyyz_yzzz_0[i] = 3.0 * ta_yyy_yzz_0[i] * fe_0 - 3.0 * ta_yyy_yzz_1[i] * fe_0 + ta_yyy_yzzz_0[i] * pa_z[i] - ta_yyy_yzzz_1[i] * pc_z[i];

        ta_yyyz_zzzz_0[i] = 2.0 * ta_yz_zzzz_0[i] * fe_0 - 2.0 * ta_yz_zzzz_1[i] * fe_0 + ta_yyz_zzzz_0[i] * pa_y[i] - ta_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : GG

    auto ta_yyzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 180);

    auto ta_yyzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 181);

    auto ta_yyzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 182);

    auto ta_yyzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 183);

    auto ta_yyzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 184);

    auto ta_yyzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 185);

    auto ta_yyzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 186);

    auto ta_yyzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 187);

    auto ta_yyzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 188);

    auto ta_yyzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 189);

    auto ta_yyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 190);

    auto ta_yyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 191);

    auto ta_yyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 192);

    auto ta_yyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 193);

    auto ta_yyzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 194);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yy_xxxy_0,   \
                             ta_yy_xxxy_1,   \
                             ta_yy_xxyy_0,   \
                             ta_yy_xxyy_1,   \
                             ta_yy_xyyy_0,   \
                             ta_yy_xyyy_1,   \
                             ta_yy_yyyy_0,   \
                             ta_yy_yyyy_1,   \
                             ta_yyz_xxxy_0,  \
                             ta_yyz_xxxy_1,  \
                             ta_yyz_xxyy_0,  \
                             ta_yyz_xxyy_1,  \
                             ta_yyz_xyyy_0,  \
                             ta_yyz_xyyy_1,  \
                             ta_yyz_yyyy_0,  \
                             ta_yyz_yyyy_1,  \
                             ta_yyzz_xxxx_0, \
                             ta_yyzz_xxxy_0, \
                             ta_yyzz_xxxz_0, \
                             ta_yyzz_xxyy_0, \
                             ta_yyzz_xxyz_0, \
                             ta_yyzz_xxzz_0, \
                             ta_yyzz_xyyy_0, \
                             ta_yyzz_xyyz_0, \
                             ta_yyzz_xyzz_0, \
                             ta_yyzz_xzzz_0, \
                             ta_yyzz_yyyy_0, \
                             ta_yyzz_yyyz_0, \
                             ta_yyzz_yyzz_0, \
                             ta_yyzz_yzzz_0, \
                             ta_yyzz_zzzz_0, \
                             ta_yzz_xxxx_0,  \
                             ta_yzz_xxxx_1,  \
                             ta_yzz_xxxz_0,  \
                             ta_yzz_xxxz_1,  \
                             ta_yzz_xxyz_0,  \
                             ta_yzz_xxyz_1,  \
                             ta_yzz_xxz_0,   \
                             ta_yzz_xxz_1,   \
                             ta_yzz_xxzz_0,  \
                             ta_yzz_xxzz_1,  \
                             ta_yzz_xyyz_0,  \
                             ta_yzz_xyyz_1,  \
                             ta_yzz_xyz_0,   \
                             ta_yzz_xyz_1,   \
                             ta_yzz_xyzz_0,  \
                             ta_yzz_xyzz_1,  \
                             ta_yzz_xzz_0,   \
                             ta_yzz_xzz_1,   \
                             ta_yzz_xzzz_0,  \
                             ta_yzz_xzzz_1,  \
                             ta_yzz_yyyz_0,  \
                             ta_yzz_yyyz_1,  \
                             ta_yzz_yyz_0,   \
                             ta_yzz_yyz_1,   \
                             ta_yzz_yyzz_0,  \
                             ta_yzz_yyzz_1,  \
                             ta_yzz_yzz_0,   \
                             ta_yzz_yzz_1,   \
                             ta_yzz_yzzz_0,  \
                             ta_yzz_yzzz_1,  \
                             ta_yzz_zzz_0,   \
                             ta_yzz_zzz_1,   \
                             ta_yzz_zzzz_0,  \
                             ta_yzz_zzzz_1,  \
                             ta_zz_xxxx_0,   \
                             ta_zz_xxxx_1,   \
                             ta_zz_xxxz_0,   \
                             ta_zz_xxxz_1,   \
                             ta_zz_xxyz_0,   \
                             ta_zz_xxyz_1,   \
                             ta_zz_xxzz_0,   \
                             ta_zz_xxzz_1,   \
                             ta_zz_xyyz_0,   \
                             ta_zz_xyyz_1,   \
                             ta_zz_xyzz_0,   \
                             ta_zz_xyzz_1,   \
                             ta_zz_xzzz_0,   \
                             ta_zz_xzzz_1,   \
                             ta_zz_yyyz_0,   \
                             ta_zz_yyyz_1,   \
                             ta_zz_yyzz_0,   \
                             ta_zz_yyzz_1,   \
                             ta_zz_yzzz_0,   \
                             ta_zz_yzzz_1,   \
                             ta_zz_zzzz_0,   \
                             ta_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_xxxx_0[i] = ta_zz_xxxx_0[i] * fe_0 - ta_zz_xxxx_1[i] * fe_0 + ta_yzz_xxxx_0[i] * pa_y[i] - ta_yzz_xxxx_1[i] * pc_y[i];

        ta_yyzz_xxxy_0[i] = ta_yy_xxxy_0[i] * fe_0 - ta_yy_xxxy_1[i] * fe_0 + ta_yyz_xxxy_0[i] * pa_z[i] - ta_yyz_xxxy_1[i] * pc_z[i];

        ta_yyzz_xxxz_0[i] = ta_zz_xxxz_0[i] * fe_0 - ta_zz_xxxz_1[i] * fe_0 + ta_yzz_xxxz_0[i] * pa_y[i] - ta_yzz_xxxz_1[i] * pc_y[i];

        ta_yyzz_xxyy_0[i] = ta_yy_xxyy_0[i] * fe_0 - ta_yy_xxyy_1[i] * fe_0 + ta_yyz_xxyy_0[i] * pa_z[i] - ta_yyz_xxyy_1[i] * pc_z[i];

        ta_yyzz_xxyz_0[i] = ta_zz_xxyz_0[i] * fe_0 - ta_zz_xxyz_1[i] * fe_0 + ta_yzz_xxz_0[i] * fe_0 - ta_yzz_xxz_1[i] * fe_0 +
                            ta_yzz_xxyz_0[i] * pa_y[i] - ta_yzz_xxyz_1[i] * pc_y[i];

        ta_yyzz_xxzz_0[i] = ta_zz_xxzz_0[i] * fe_0 - ta_zz_xxzz_1[i] * fe_0 + ta_yzz_xxzz_0[i] * pa_y[i] - ta_yzz_xxzz_1[i] * pc_y[i];

        ta_yyzz_xyyy_0[i] = ta_yy_xyyy_0[i] * fe_0 - ta_yy_xyyy_1[i] * fe_0 + ta_yyz_xyyy_0[i] * pa_z[i] - ta_yyz_xyyy_1[i] * pc_z[i];

        ta_yyzz_xyyz_0[i] = ta_zz_xyyz_0[i] * fe_0 - ta_zz_xyyz_1[i] * fe_0 + 2.0 * ta_yzz_xyz_0[i] * fe_0 - 2.0 * ta_yzz_xyz_1[i] * fe_0 +
                            ta_yzz_xyyz_0[i] * pa_y[i] - ta_yzz_xyyz_1[i] * pc_y[i];

        ta_yyzz_xyzz_0[i] = ta_zz_xyzz_0[i] * fe_0 - ta_zz_xyzz_1[i] * fe_0 + ta_yzz_xzz_0[i] * fe_0 - ta_yzz_xzz_1[i] * fe_0 +
                            ta_yzz_xyzz_0[i] * pa_y[i] - ta_yzz_xyzz_1[i] * pc_y[i];

        ta_yyzz_xzzz_0[i] = ta_zz_xzzz_0[i] * fe_0 - ta_zz_xzzz_1[i] * fe_0 + ta_yzz_xzzz_0[i] * pa_y[i] - ta_yzz_xzzz_1[i] * pc_y[i];

        ta_yyzz_yyyy_0[i] = ta_yy_yyyy_0[i] * fe_0 - ta_yy_yyyy_1[i] * fe_0 + ta_yyz_yyyy_0[i] * pa_z[i] - ta_yyz_yyyy_1[i] * pc_z[i];

        ta_yyzz_yyyz_0[i] = ta_zz_yyyz_0[i] * fe_0 - ta_zz_yyyz_1[i] * fe_0 + 3.0 * ta_yzz_yyz_0[i] * fe_0 - 3.0 * ta_yzz_yyz_1[i] * fe_0 +
                            ta_yzz_yyyz_0[i] * pa_y[i] - ta_yzz_yyyz_1[i] * pc_y[i];

        ta_yyzz_yyzz_0[i] = ta_zz_yyzz_0[i] * fe_0 - ta_zz_yyzz_1[i] * fe_0 + 2.0 * ta_yzz_yzz_0[i] * fe_0 - 2.0 * ta_yzz_yzz_1[i] * fe_0 +
                            ta_yzz_yyzz_0[i] * pa_y[i] - ta_yzz_yyzz_1[i] * pc_y[i];

        ta_yyzz_yzzz_0[i] = ta_zz_yzzz_0[i] * fe_0 - ta_zz_yzzz_1[i] * fe_0 + ta_yzz_zzz_0[i] * fe_0 - ta_yzz_zzz_1[i] * fe_0 +
                            ta_yzz_yzzz_0[i] * pa_y[i] - ta_yzz_yzzz_1[i] * pc_y[i];

        ta_yyzz_zzzz_0[i] = ta_zz_zzzz_0[i] * fe_0 - ta_zz_zzzz_1[i] * fe_0 + ta_yzz_zzzz_0[i] * pa_y[i] - ta_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto ta_yzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 195);

    auto ta_yzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 196);

    auto ta_yzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 197);

    auto ta_yzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 198);

    auto ta_yzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 199);

    auto ta_yzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 200);

    auto ta_yzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 201);

    auto ta_yzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 202);

    auto ta_yzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 203);

    auto ta_yzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 204);

    auto ta_yzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 205);

    auto ta_yzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 206);

    auto ta_yzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 207);

    auto ta_yzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 208);

    auto ta_yzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 209);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yzzz_xxxx_0, \
                             ta_yzzz_xxxy_0, \
                             ta_yzzz_xxxz_0, \
                             ta_yzzz_xxyy_0, \
                             ta_yzzz_xxyz_0, \
                             ta_yzzz_xxzz_0, \
                             ta_yzzz_xyyy_0, \
                             ta_yzzz_xyyz_0, \
                             ta_yzzz_xyzz_0, \
                             ta_yzzz_xzzz_0, \
                             ta_yzzz_yyyy_0, \
                             ta_yzzz_yyyz_0, \
                             ta_yzzz_yyzz_0, \
                             ta_yzzz_yzzz_0, \
                             ta_yzzz_zzzz_0, \
                             ta_zzz_xxx_0,   \
                             ta_zzz_xxx_1,   \
                             ta_zzz_xxxx_0,  \
                             ta_zzz_xxxx_1,  \
                             ta_zzz_xxxy_0,  \
                             ta_zzz_xxxy_1,  \
                             ta_zzz_xxxz_0,  \
                             ta_zzz_xxxz_1,  \
                             ta_zzz_xxy_0,   \
                             ta_zzz_xxy_1,   \
                             ta_zzz_xxyy_0,  \
                             ta_zzz_xxyy_1,  \
                             ta_zzz_xxyz_0,  \
                             ta_zzz_xxyz_1,  \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xxzz_0,  \
                             ta_zzz_xxzz_1,  \
                             ta_zzz_xyy_0,   \
                             ta_zzz_xyy_1,   \
                             ta_zzz_xyyy_0,  \
                             ta_zzz_xyyy_1,  \
                             ta_zzz_xyyz_0,  \
                             ta_zzz_xyyz_1,  \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xyzz_0,  \
                             ta_zzz_xyzz_1,  \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_xzzz_0,  \
                             ta_zzz_xzzz_1,  \
                             ta_zzz_yyy_0,   \
                             ta_zzz_yyy_1,   \
                             ta_zzz_yyyy_0,  \
                             ta_zzz_yyyy_1,  \
                             ta_zzz_yyyz_0,  \
                             ta_zzz_yyyz_1,  \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yyzz_0,  \
                             ta_zzz_yyzz_1,  \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_yzzz_0,  \
                             ta_zzz_yzzz_1,  \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             ta_zzz_zzzz_0,  \
                             ta_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_xxxx_0[i] = ta_zzz_xxxx_0[i] * pa_y[i] - ta_zzz_xxxx_1[i] * pc_y[i];

        ta_yzzz_xxxy_0[i] = ta_zzz_xxx_0[i] * fe_0 - ta_zzz_xxx_1[i] * fe_0 + ta_zzz_xxxy_0[i] * pa_y[i] - ta_zzz_xxxy_1[i] * pc_y[i];

        ta_yzzz_xxxz_0[i] = ta_zzz_xxxz_0[i] * pa_y[i] - ta_zzz_xxxz_1[i] * pc_y[i];

        ta_yzzz_xxyy_0[i] = 2.0 * ta_zzz_xxy_0[i] * fe_0 - 2.0 * ta_zzz_xxy_1[i] * fe_0 + ta_zzz_xxyy_0[i] * pa_y[i] - ta_zzz_xxyy_1[i] * pc_y[i];

        ta_yzzz_xxyz_0[i] = ta_zzz_xxz_0[i] * fe_0 - ta_zzz_xxz_1[i] * fe_0 + ta_zzz_xxyz_0[i] * pa_y[i] - ta_zzz_xxyz_1[i] * pc_y[i];

        ta_yzzz_xxzz_0[i] = ta_zzz_xxzz_0[i] * pa_y[i] - ta_zzz_xxzz_1[i] * pc_y[i];

        ta_yzzz_xyyy_0[i] = 3.0 * ta_zzz_xyy_0[i] * fe_0 - 3.0 * ta_zzz_xyy_1[i] * fe_0 + ta_zzz_xyyy_0[i] * pa_y[i] - ta_zzz_xyyy_1[i] * pc_y[i];

        ta_yzzz_xyyz_0[i] = 2.0 * ta_zzz_xyz_0[i] * fe_0 - 2.0 * ta_zzz_xyz_1[i] * fe_0 + ta_zzz_xyyz_0[i] * pa_y[i] - ta_zzz_xyyz_1[i] * pc_y[i];

        ta_yzzz_xyzz_0[i] = ta_zzz_xzz_0[i] * fe_0 - ta_zzz_xzz_1[i] * fe_0 + ta_zzz_xyzz_0[i] * pa_y[i] - ta_zzz_xyzz_1[i] * pc_y[i];

        ta_yzzz_xzzz_0[i] = ta_zzz_xzzz_0[i] * pa_y[i] - ta_zzz_xzzz_1[i] * pc_y[i];

        ta_yzzz_yyyy_0[i] = 4.0 * ta_zzz_yyy_0[i] * fe_0 - 4.0 * ta_zzz_yyy_1[i] * fe_0 + ta_zzz_yyyy_0[i] * pa_y[i] - ta_zzz_yyyy_1[i] * pc_y[i];

        ta_yzzz_yyyz_0[i] = 3.0 * ta_zzz_yyz_0[i] * fe_0 - 3.0 * ta_zzz_yyz_1[i] * fe_0 + ta_zzz_yyyz_0[i] * pa_y[i] - ta_zzz_yyyz_1[i] * pc_y[i];

        ta_yzzz_yyzz_0[i] = 2.0 * ta_zzz_yzz_0[i] * fe_0 - 2.0 * ta_zzz_yzz_1[i] * fe_0 + ta_zzz_yyzz_0[i] * pa_y[i] - ta_zzz_yyzz_1[i] * pc_y[i];

        ta_yzzz_yzzz_0[i] = ta_zzz_zzz_0[i] * fe_0 - ta_zzz_zzz_1[i] * fe_0 + ta_zzz_yzzz_0[i] * pa_y[i] - ta_zzz_yzzz_1[i] * pc_y[i];

        ta_yzzz_zzzz_0[i] = ta_zzz_zzzz_0[i] * pa_y[i] - ta_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

    auto ta_zzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 210);

    auto ta_zzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 211);

    auto ta_zzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 212);

    auto ta_zzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 213);

    auto ta_zzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 214);

    auto ta_zzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 215);

    auto ta_zzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 216);

    auto ta_zzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 217);

    auto ta_zzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 218);

    auto ta_zzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 219);

    auto ta_zzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 220);

    auto ta_zzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 221);

    auto ta_zzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 222);

    auto ta_zzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 223);

    auto ta_zzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 224);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta_zz_xxxx_0,   \
                             ta_zz_xxxx_1,   \
                             ta_zz_xxxy_0,   \
                             ta_zz_xxxy_1,   \
                             ta_zz_xxxz_0,   \
                             ta_zz_xxxz_1,   \
                             ta_zz_xxyy_0,   \
                             ta_zz_xxyy_1,   \
                             ta_zz_xxyz_0,   \
                             ta_zz_xxyz_1,   \
                             ta_zz_xxzz_0,   \
                             ta_zz_xxzz_1,   \
                             ta_zz_xyyy_0,   \
                             ta_zz_xyyy_1,   \
                             ta_zz_xyyz_0,   \
                             ta_zz_xyyz_1,   \
                             ta_zz_xyzz_0,   \
                             ta_zz_xyzz_1,   \
                             ta_zz_xzzz_0,   \
                             ta_zz_xzzz_1,   \
                             ta_zz_yyyy_0,   \
                             ta_zz_yyyy_1,   \
                             ta_zz_yyyz_0,   \
                             ta_zz_yyyz_1,   \
                             ta_zz_yyzz_0,   \
                             ta_zz_yyzz_1,   \
                             ta_zz_yzzz_0,   \
                             ta_zz_yzzz_1,   \
                             ta_zz_zzzz_0,   \
                             ta_zz_zzzz_1,   \
                             ta_zzz_xxx_0,   \
                             ta_zzz_xxx_1,   \
                             ta_zzz_xxxx_0,  \
                             ta_zzz_xxxx_1,  \
                             ta_zzz_xxxy_0,  \
                             ta_zzz_xxxy_1,  \
                             ta_zzz_xxxz_0,  \
                             ta_zzz_xxxz_1,  \
                             ta_zzz_xxy_0,   \
                             ta_zzz_xxy_1,   \
                             ta_zzz_xxyy_0,  \
                             ta_zzz_xxyy_1,  \
                             ta_zzz_xxyz_0,  \
                             ta_zzz_xxyz_1,  \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xxzz_0,  \
                             ta_zzz_xxzz_1,  \
                             ta_zzz_xyy_0,   \
                             ta_zzz_xyy_1,   \
                             ta_zzz_xyyy_0,  \
                             ta_zzz_xyyy_1,  \
                             ta_zzz_xyyz_0,  \
                             ta_zzz_xyyz_1,  \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xyzz_0,  \
                             ta_zzz_xyzz_1,  \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_xzzz_0,  \
                             ta_zzz_xzzz_1,  \
                             ta_zzz_yyy_0,   \
                             ta_zzz_yyy_1,   \
                             ta_zzz_yyyy_0,  \
                             ta_zzz_yyyy_1,  \
                             ta_zzz_yyyz_0,  \
                             ta_zzz_yyyz_1,  \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yyzz_0,  \
                             ta_zzz_yyzz_1,  \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_yzzz_0,  \
                             ta_zzz_yzzz_1,  \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             ta_zzz_zzzz_0,  \
                             ta_zzz_zzzz_1,  \
                             ta_zzzz_xxxx_0, \
                             ta_zzzz_xxxy_0, \
                             ta_zzzz_xxxz_0, \
                             ta_zzzz_xxyy_0, \
                             ta_zzzz_xxyz_0, \
                             ta_zzzz_xxzz_0, \
                             ta_zzzz_xyyy_0, \
                             ta_zzzz_xyyz_0, \
                             ta_zzzz_xyzz_0, \
                             ta_zzzz_xzzz_0, \
                             ta_zzzz_yyyy_0, \
                             ta_zzzz_yyyz_0, \
                             ta_zzzz_yyzz_0, \
                             ta_zzzz_yzzz_0, \
                             ta_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_xxxx_0[i] = 3.0 * ta_zz_xxxx_0[i] * fe_0 - 3.0 * ta_zz_xxxx_1[i] * fe_0 + ta_zzz_xxxx_0[i] * pa_z[i] - ta_zzz_xxxx_1[i] * pc_z[i];

        ta_zzzz_xxxy_0[i] = 3.0 * ta_zz_xxxy_0[i] * fe_0 - 3.0 * ta_zz_xxxy_1[i] * fe_0 + ta_zzz_xxxy_0[i] * pa_z[i] - ta_zzz_xxxy_1[i] * pc_z[i];

        ta_zzzz_xxxz_0[i] = 3.0 * ta_zz_xxxz_0[i] * fe_0 - 3.0 * ta_zz_xxxz_1[i] * fe_0 + ta_zzz_xxx_0[i] * fe_0 - ta_zzz_xxx_1[i] * fe_0 +
                            ta_zzz_xxxz_0[i] * pa_z[i] - ta_zzz_xxxz_1[i] * pc_z[i];

        ta_zzzz_xxyy_0[i] = 3.0 * ta_zz_xxyy_0[i] * fe_0 - 3.0 * ta_zz_xxyy_1[i] * fe_0 + ta_zzz_xxyy_0[i] * pa_z[i] - ta_zzz_xxyy_1[i] * pc_z[i];

        ta_zzzz_xxyz_0[i] = 3.0 * ta_zz_xxyz_0[i] * fe_0 - 3.0 * ta_zz_xxyz_1[i] * fe_0 + ta_zzz_xxy_0[i] * fe_0 - ta_zzz_xxy_1[i] * fe_0 +
                            ta_zzz_xxyz_0[i] * pa_z[i] - ta_zzz_xxyz_1[i] * pc_z[i];

        ta_zzzz_xxzz_0[i] = 3.0 * ta_zz_xxzz_0[i] * fe_0 - 3.0 * ta_zz_xxzz_1[i] * fe_0 + 2.0 * ta_zzz_xxz_0[i] * fe_0 -
                            2.0 * ta_zzz_xxz_1[i] * fe_0 + ta_zzz_xxzz_0[i] * pa_z[i] - ta_zzz_xxzz_1[i] * pc_z[i];

        ta_zzzz_xyyy_0[i] = 3.0 * ta_zz_xyyy_0[i] * fe_0 - 3.0 * ta_zz_xyyy_1[i] * fe_0 + ta_zzz_xyyy_0[i] * pa_z[i] - ta_zzz_xyyy_1[i] * pc_z[i];

        ta_zzzz_xyyz_0[i] = 3.0 * ta_zz_xyyz_0[i] * fe_0 - 3.0 * ta_zz_xyyz_1[i] * fe_0 + ta_zzz_xyy_0[i] * fe_0 - ta_zzz_xyy_1[i] * fe_0 +
                            ta_zzz_xyyz_0[i] * pa_z[i] - ta_zzz_xyyz_1[i] * pc_z[i];

        ta_zzzz_xyzz_0[i] = 3.0 * ta_zz_xyzz_0[i] * fe_0 - 3.0 * ta_zz_xyzz_1[i] * fe_0 + 2.0 * ta_zzz_xyz_0[i] * fe_0 -
                            2.0 * ta_zzz_xyz_1[i] * fe_0 + ta_zzz_xyzz_0[i] * pa_z[i] - ta_zzz_xyzz_1[i] * pc_z[i];

        ta_zzzz_xzzz_0[i] = 3.0 * ta_zz_xzzz_0[i] * fe_0 - 3.0 * ta_zz_xzzz_1[i] * fe_0 + 3.0 * ta_zzz_xzz_0[i] * fe_0 -
                            3.0 * ta_zzz_xzz_1[i] * fe_0 + ta_zzz_xzzz_0[i] * pa_z[i] - ta_zzz_xzzz_1[i] * pc_z[i];

        ta_zzzz_yyyy_0[i] = 3.0 * ta_zz_yyyy_0[i] * fe_0 - 3.0 * ta_zz_yyyy_1[i] * fe_0 + ta_zzz_yyyy_0[i] * pa_z[i] - ta_zzz_yyyy_1[i] * pc_z[i];

        ta_zzzz_yyyz_0[i] = 3.0 * ta_zz_yyyz_0[i] * fe_0 - 3.0 * ta_zz_yyyz_1[i] * fe_0 + ta_zzz_yyy_0[i] * fe_0 - ta_zzz_yyy_1[i] * fe_0 +
                            ta_zzz_yyyz_0[i] * pa_z[i] - ta_zzz_yyyz_1[i] * pc_z[i];

        ta_zzzz_yyzz_0[i] = 3.0 * ta_zz_yyzz_0[i] * fe_0 - 3.0 * ta_zz_yyzz_1[i] * fe_0 + 2.0 * ta_zzz_yyz_0[i] * fe_0 -
                            2.0 * ta_zzz_yyz_1[i] * fe_0 + ta_zzz_yyzz_0[i] * pa_z[i] - ta_zzz_yyzz_1[i] * pc_z[i];

        ta_zzzz_yzzz_0[i] = 3.0 * ta_zz_yzzz_0[i] * fe_0 - 3.0 * ta_zz_yzzz_1[i] * fe_0 + 3.0 * ta_zzz_yzz_0[i] * fe_0 -
                            3.0 * ta_zzz_yzz_1[i] * fe_0 + ta_zzz_yzzz_0[i] * pa_z[i] - ta_zzz_yzzz_1[i] * pc_z[i];

        ta_zzzz_zzzz_0[i] = 3.0 * ta_zz_zzzz_0[i] * fe_0 - 3.0 * ta_zz_zzzz_1[i] * fe_0 + 4.0 * ta_zzz_zzz_0[i] * fe_0 -
                            4.0 * ta_zzz_zzz_1[i] * fe_0 + ta_zzz_zzzz_0[i] * pa_z[i] - ta_zzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
