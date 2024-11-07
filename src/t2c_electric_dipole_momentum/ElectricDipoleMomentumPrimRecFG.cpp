#include "ElectricDipoleMomentumPrimRecFG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_fg(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_fg,
                                      const size_t              idx_dip_pg,
                                      const size_t              idx_dip_df,
                                      const size_t              idx_ovl_dg,
                                      const size_t              idx_dip_dg,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : PG

    auto tr_x_x_xxxx = pbuffer.data(idx_dip_pg);

    auto tr_x_x_xxxy = pbuffer.data(idx_dip_pg + 1);

    auto tr_x_x_xxxz = pbuffer.data(idx_dip_pg + 2);

    auto tr_x_x_xxyy = pbuffer.data(idx_dip_pg + 3);

    auto tr_x_x_xxyz = pbuffer.data(idx_dip_pg + 4);

    auto tr_x_x_xxzz = pbuffer.data(idx_dip_pg + 5);

    auto tr_x_x_xyyy = pbuffer.data(idx_dip_pg + 6);

    auto tr_x_x_xyyz = pbuffer.data(idx_dip_pg + 7);

    auto tr_x_x_xyzz = pbuffer.data(idx_dip_pg + 8);

    auto tr_x_x_xzzz = pbuffer.data(idx_dip_pg + 9);

    auto tr_x_x_yyyy = pbuffer.data(idx_dip_pg + 10);

    auto tr_x_x_yyyz = pbuffer.data(idx_dip_pg + 11);

    auto tr_x_x_yyzz = pbuffer.data(idx_dip_pg + 12);

    auto tr_x_x_yzzz = pbuffer.data(idx_dip_pg + 13);

    auto tr_x_x_zzzz = pbuffer.data(idx_dip_pg + 14);

    auto tr_x_y_xxxx = pbuffer.data(idx_dip_pg + 15);

    auto tr_x_y_xxxy = pbuffer.data(idx_dip_pg + 16);

    auto tr_x_y_xxxz = pbuffer.data(idx_dip_pg + 17);

    auto tr_x_y_xxyy = pbuffer.data(idx_dip_pg + 18);

    auto tr_x_y_xxyz = pbuffer.data(idx_dip_pg + 19);

    auto tr_x_y_xxzz = pbuffer.data(idx_dip_pg + 20);

    auto tr_x_y_xyyy = pbuffer.data(idx_dip_pg + 21);

    auto tr_x_y_xyyz = pbuffer.data(idx_dip_pg + 22);

    auto tr_x_y_xyzz = pbuffer.data(idx_dip_pg + 23);

    auto tr_x_y_xzzz = pbuffer.data(idx_dip_pg + 24);

    auto tr_x_y_yyyy = pbuffer.data(idx_dip_pg + 25);

    auto tr_x_y_yyyz = pbuffer.data(idx_dip_pg + 26);

    auto tr_x_y_yyzz = pbuffer.data(idx_dip_pg + 27);

    auto tr_x_y_yzzz = pbuffer.data(idx_dip_pg + 28);

    auto tr_x_y_zzzz = pbuffer.data(idx_dip_pg + 29);

    auto tr_x_z_xxxx = pbuffer.data(idx_dip_pg + 30);

    auto tr_x_z_xxxy = pbuffer.data(idx_dip_pg + 31);

    auto tr_x_z_xxxz = pbuffer.data(idx_dip_pg + 32);

    auto tr_x_z_xxyy = pbuffer.data(idx_dip_pg + 33);

    auto tr_x_z_xxyz = pbuffer.data(idx_dip_pg + 34);

    auto tr_x_z_xxzz = pbuffer.data(idx_dip_pg + 35);

    auto tr_x_z_xyyy = pbuffer.data(idx_dip_pg + 36);

    auto tr_x_z_xyyz = pbuffer.data(idx_dip_pg + 37);

    auto tr_x_z_xyzz = pbuffer.data(idx_dip_pg + 38);

    auto tr_x_z_xzzz = pbuffer.data(idx_dip_pg + 39);

    auto tr_x_z_yyyy = pbuffer.data(idx_dip_pg + 40);

    auto tr_x_z_yyyz = pbuffer.data(idx_dip_pg + 41);

    auto tr_x_z_yyzz = pbuffer.data(idx_dip_pg + 42);

    auto tr_x_z_yzzz = pbuffer.data(idx_dip_pg + 43);

    auto tr_x_z_zzzz = pbuffer.data(idx_dip_pg + 44);

    auto tr_y_x_xxxx = pbuffer.data(idx_dip_pg + 45);

    auto tr_y_x_xxxy = pbuffer.data(idx_dip_pg + 46);

    auto tr_y_x_xxxz = pbuffer.data(idx_dip_pg + 47);

    auto tr_y_x_xxyy = pbuffer.data(idx_dip_pg + 48);

    auto tr_y_x_xxyz = pbuffer.data(idx_dip_pg + 49);

    auto tr_y_x_xxzz = pbuffer.data(idx_dip_pg + 50);

    auto tr_y_x_xyyy = pbuffer.data(idx_dip_pg + 51);

    auto tr_y_x_xyyz = pbuffer.data(idx_dip_pg + 52);

    auto tr_y_x_xyzz = pbuffer.data(idx_dip_pg + 53);

    auto tr_y_x_xzzz = pbuffer.data(idx_dip_pg + 54);

    auto tr_y_x_yyyy = pbuffer.data(idx_dip_pg + 55);

    auto tr_y_x_yyyz = pbuffer.data(idx_dip_pg + 56);

    auto tr_y_x_yyzz = pbuffer.data(idx_dip_pg + 57);

    auto tr_y_x_yzzz = pbuffer.data(idx_dip_pg + 58);

    auto tr_y_x_zzzz = pbuffer.data(idx_dip_pg + 59);

    auto tr_y_y_xxxx = pbuffer.data(idx_dip_pg + 60);

    auto tr_y_y_xxxy = pbuffer.data(idx_dip_pg + 61);

    auto tr_y_y_xxxz = pbuffer.data(idx_dip_pg + 62);

    auto tr_y_y_xxyy = pbuffer.data(idx_dip_pg + 63);

    auto tr_y_y_xxyz = pbuffer.data(idx_dip_pg + 64);

    auto tr_y_y_xxzz = pbuffer.data(idx_dip_pg + 65);

    auto tr_y_y_xyyy = pbuffer.data(idx_dip_pg + 66);

    auto tr_y_y_xyyz = pbuffer.data(idx_dip_pg + 67);

    auto tr_y_y_xyzz = pbuffer.data(idx_dip_pg + 68);

    auto tr_y_y_xzzz = pbuffer.data(idx_dip_pg + 69);

    auto tr_y_y_yyyy = pbuffer.data(idx_dip_pg + 70);

    auto tr_y_y_yyyz = pbuffer.data(idx_dip_pg + 71);

    auto tr_y_y_yyzz = pbuffer.data(idx_dip_pg + 72);

    auto tr_y_y_yzzz = pbuffer.data(idx_dip_pg + 73);

    auto tr_y_y_zzzz = pbuffer.data(idx_dip_pg + 74);

    auto tr_y_z_xxxx = pbuffer.data(idx_dip_pg + 75);

    auto tr_y_z_xxxy = pbuffer.data(idx_dip_pg + 76);

    auto tr_y_z_xxxz = pbuffer.data(idx_dip_pg + 77);

    auto tr_y_z_xxyy = pbuffer.data(idx_dip_pg + 78);

    auto tr_y_z_xxyz = pbuffer.data(idx_dip_pg + 79);

    auto tr_y_z_xxzz = pbuffer.data(idx_dip_pg + 80);

    auto tr_y_z_xyyy = pbuffer.data(idx_dip_pg + 81);

    auto tr_y_z_xyyz = pbuffer.data(idx_dip_pg + 82);

    auto tr_y_z_xyzz = pbuffer.data(idx_dip_pg + 83);

    auto tr_y_z_xzzz = pbuffer.data(idx_dip_pg + 84);

    auto tr_y_z_yyyy = pbuffer.data(idx_dip_pg + 85);

    auto tr_y_z_yyyz = pbuffer.data(idx_dip_pg + 86);

    auto tr_y_z_yyzz = pbuffer.data(idx_dip_pg + 87);

    auto tr_y_z_yzzz = pbuffer.data(idx_dip_pg + 88);

    auto tr_y_z_zzzz = pbuffer.data(idx_dip_pg + 89);

    auto tr_z_x_xxxx = pbuffer.data(idx_dip_pg + 90);

    auto tr_z_x_xxxy = pbuffer.data(idx_dip_pg + 91);

    auto tr_z_x_xxxz = pbuffer.data(idx_dip_pg + 92);

    auto tr_z_x_xxyy = pbuffer.data(idx_dip_pg + 93);

    auto tr_z_x_xxyz = pbuffer.data(idx_dip_pg + 94);

    auto tr_z_x_xxzz = pbuffer.data(idx_dip_pg + 95);

    auto tr_z_x_xyyy = pbuffer.data(idx_dip_pg + 96);

    auto tr_z_x_xyyz = pbuffer.data(idx_dip_pg + 97);

    auto tr_z_x_xyzz = pbuffer.data(idx_dip_pg + 98);

    auto tr_z_x_xzzz = pbuffer.data(idx_dip_pg + 99);

    auto tr_z_x_yyyy = pbuffer.data(idx_dip_pg + 100);

    auto tr_z_x_yyyz = pbuffer.data(idx_dip_pg + 101);

    auto tr_z_x_yyzz = pbuffer.data(idx_dip_pg + 102);

    auto tr_z_x_yzzz = pbuffer.data(idx_dip_pg + 103);

    auto tr_z_x_zzzz = pbuffer.data(idx_dip_pg + 104);

    auto tr_z_y_xxxx = pbuffer.data(idx_dip_pg + 105);

    auto tr_z_y_xxxy = pbuffer.data(idx_dip_pg + 106);

    auto tr_z_y_xxxz = pbuffer.data(idx_dip_pg + 107);

    auto tr_z_y_xxyy = pbuffer.data(idx_dip_pg + 108);

    auto tr_z_y_xxyz = pbuffer.data(idx_dip_pg + 109);

    auto tr_z_y_xxzz = pbuffer.data(idx_dip_pg + 110);

    auto tr_z_y_xyyy = pbuffer.data(idx_dip_pg + 111);

    auto tr_z_y_xyyz = pbuffer.data(idx_dip_pg + 112);

    auto tr_z_y_xyzz = pbuffer.data(idx_dip_pg + 113);

    auto tr_z_y_xzzz = pbuffer.data(idx_dip_pg + 114);

    auto tr_z_y_yyyy = pbuffer.data(idx_dip_pg + 115);

    auto tr_z_y_yyyz = pbuffer.data(idx_dip_pg + 116);

    auto tr_z_y_yyzz = pbuffer.data(idx_dip_pg + 117);

    auto tr_z_y_yzzz = pbuffer.data(idx_dip_pg + 118);

    auto tr_z_y_zzzz = pbuffer.data(idx_dip_pg + 119);

    auto tr_z_z_xxxx = pbuffer.data(idx_dip_pg + 120);

    auto tr_z_z_xxxy = pbuffer.data(idx_dip_pg + 121);

    auto tr_z_z_xxxz = pbuffer.data(idx_dip_pg + 122);

    auto tr_z_z_xxyy = pbuffer.data(idx_dip_pg + 123);

    auto tr_z_z_xxyz = pbuffer.data(idx_dip_pg + 124);

    auto tr_z_z_xxzz = pbuffer.data(idx_dip_pg + 125);

    auto tr_z_z_xyyy = pbuffer.data(idx_dip_pg + 126);

    auto tr_z_z_xyyz = pbuffer.data(idx_dip_pg + 127);

    auto tr_z_z_xyzz = pbuffer.data(idx_dip_pg + 128);

    auto tr_z_z_xzzz = pbuffer.data(idx_dip_pg + 129);

    auto tr_z_z_yyyy = pbuffer.data(idx_dip_pg + 130);

    auto tr_z_z_yyyz = pbuffer.data(idx_dip_pg + 131);

    auto tr_z_z_yyzz = pbuffer.data(idx_dip_pg + 132);

    auto tr_z_z_yzzz = pbuffer.data(idx_dip_pg + 133);

    auto tr_z_z_zzzz = pbuffer.data(idx_dip_pg + 134);

    // Set up components of auxiliary buffer : DF

    auto tr_x_xx_xxx = pbuffer.data(idx_dip_df);

    auto tr_x_xx_xxy = pbuffer.data(idx_dip_df + 1);

    auto tr_x_xx_xxz = pbuffer.data(idx_dip_df + 2);

    auto tr_x_xx_xyy = pbuffer.data(idx_dip_df + 3);

    auto tr_x_xx_xyz = pbuffer.data(idx_dip_df + 4);

    auto tr_x_xx_xzz = pbuffer.data(idx_dip_df + 5);

    auto tr_x_xx_yyy = pbuffer.data(idx_dip_df + 6);

    auto tr_x_xx_yyz = pbuffer.data(idx_dip_df + 7);

    auto tr_x_xx_yzz = pbuffer.data(idx_dip_df + 8);

    auto tr_x_xx_zzz = pbuffer.data(idx_dip_df + 9);

    auto tr_x_xz_xxz = pbuffer.data(idx_dip_df + 22);

    auto tr_x_xz_xyz = pbuffer.data(idx_dip_df + 24);

    auto tr_x_xz_xzz = pbuffer.data(idx_dip_df + 25);

    auto tr_x_yy_xxx = pbuffer.data(idx_dip_df + 30);

    auto tr_x_yy_xxy = pbuffer.data(idx_dip_df + 31);

    auto tr_x_yy_xxz = pbuffer.data(idx_dip_df + 32);

    auto tr_x_yy_xyy = pbuffer.data(idx_dip_df + 33);

    auto tr_x_yy_xyz = pbuffer.data(idx_dip_df + 34);

    auto tr_x_yy_xzz = pbuffer.data(idx_dip_df + 35);

    auto tr_x_yy_yyy = pbuffer.data(idx_dip_df + 36);

    auto tr_x_yy_yyz = pbuffer.data(idx_dip_df + 37);

    auto tr_x_yy_yzz = pbuffer.data(idx_dip_df + 38);

    auto tr_x_yy_zzz = pbuffer.data(idx_dip_df + 39);

    auto tr_x_zz_xxx = pbuffer.data(idx_dip_df + 50);

    auto tr_x_zz_xxy = pbuffer.data(idx_dip_df + 51);

    auto tr_x_zz_xxz = pbuffer.data(idx_dip_df + 52);

    auto tr_x_zz_xyy = pbuffer.data(idx_dip_df + 53);

    auto tr_x_zz_xyz = pbuffer.data(idx_dip_df + 54);

    auto tr_x_zz_xzz = pbuffer.data(idx_dip_df + 55);

    auto tr_x_zz_yyy = pbuffer.data(idx_dip_df + 56);

    auto tr_x_zz_yyz = pbuffer.data(idx_dip_df + 57);

    auto tr_x_zz_yzz = pbuffer.data(idx_dip_df + 58);

    auto tr_x_zz_zzz = pbuffer.data(idx_dip_df + 59);

    auto tr_y_xx_xxx = pbuffer.data(idx_dip_df + 60);

    auto tr_y_xx_xxy = pbuffer.data(idx_dip_df + 61);

    auto tr_y_xx_xxz = pbuffer.data(idx_dip_df + 62);

    auto tr_y_xx_xyy = pbuffer.data(idx_dip_df + 63);

    auto tr_y_xx_xyz = pbuffer.data(idx_dip_df + 64);

    auto tr_y_xx_xzz = pbuffer.data(idx_dip_df + 65);

    auto tr_y_xx_yyy = pbuffer.data(idx_dip_df + 66);

    auto tr_y_xx_yyz = pbuffer.data(idx_dip_df + 67);

    auto tr_y_xx_yzz = pbuffer.data(idx_dip_df + 68);

    auto tr_y_xx_zzz = pbuffer.data(idx_dip_df + 69);

    auto tr_y_xy_xxy = pbuffer.data(idx_dip_df + 71);

    auto tr_y_xy_xyy = pbuffer.data(idx_dip_df + 73);

    auto tr_y_xy_xyz = pbuffer.data(idx_dip_df + 74);

    auto tr_y_xy_yyy = pbuffer.data(idx_dip_df + 76);

    auto tr_y_xy_yyz = pbuffer.data(idx_dip_df + 77);

    auto tr_y_xy_yzz = pbuffer.data(idx_dip_df + 78);

    auto tr_y_yy_xxx = pbuffer.data(idx_dip_df + 90);

    auto tr_y_yy_xxy = pbuffer.data(idx_dip_df + 91);

    auto tr_y_yy_xxz = pbuffer.data(idx_dip_df + 92);

    auto tr_y_yy_xyy = pbuffer.data(idx_dip_df + 93);

    auto tr_y_yy_xyz = pbuffer.data(idx_dip_df + 94);

    auto tr_y_yy_xzz = pbuffer.data(idx_dip_df + 95);

    auto tr_y_yy_yyy = pbuffer.data(idx_dip_df + 96);

    auto tr_y_yy_yyz = pbuffer.data(idx_dip_df + 97);

    auto tr_y_yy_yzz = pbuffer.data(idx_dip_df + 98);

    auto tr_y_yy_zzz = pbuffer.data(idx_dip_df + 99);

    auto tr_y_yz_xxz = pbuffer.data(idx_dip_df + 102);

    auto tr_y_yz_xyz = pbuffer.data(idx_dip_df + 104);

    auto tr_y_yz_xzz = pbuffer.data(idx_dip_df + 105);

    auto tr_y_yz_yyz = pbuffer.data(idx_dip_df + 107);

    auto tr_y_yz_yzz = pbuffer.data(idx_dip_df + 108);

    auto tr_y_yz_zzz = pbuffer.data(idx_dip_df + 109);

    auto tr_y_zz_xxx = pbuffer.data(idx_dip_df + 110);

    auto tr_y_zz_xxy = pbuffer.data(idx_dip_df + 111);

    auto tr_y_zz_xxz = pbuffer.data(idx_dip_df + 112);

    auto tr_y_zz_xyy = pbuffer.data(idx_dip_df + 113);

    auto tr_y_zz_xyz = pbuffer.data(idx_dip_df + 114);

    auto tr_y_zz_xzz = pbuffer.data(idx_dip_df + 115);

    auto tr_y_zz_yyy = pbuffer.data(idx_dip_df + 116);

    auto tr_y_zz_yyz = pbuffer.data(idx_dip_df + 117);

    auto tr_y_zz_yzz = pbuffer.data(idx_dip_df + 118);

    auto tr_y_zz_zzz = pbuffer.data(idx_dip_df + 119);

    auto tr_z_xx_xxx = pbuffer.data(idx_dip_df + 120);

    auto tr_z_xx_xxy = pbuffer.data(idx_dip_df + 121);

    auto tr_z_xx_xxz = pbuffer.data(idx_dip_df + 122);

    auto tr_z_xx_xyy = pbuffer.data(idx_dip_df + 123);

    auto tr_z_xx_xyz = pbuffer.data(idx_dip_df + 124);

    auto tr_z_xx_xzz = pbuffer.data(idx_dip_df + 125);

    auto tr_z_xx_yyy = pbuffer.data(idx_dip_df + 126);

    auto tr_z_xx_yyz = pbuffer.data(idx_dip_df + 127);

    auto tr_z_xx_yzz = pbuffer.data(idx_dip_df + 128);

    auto tr_z_xx_zzz = pbuffer.data(idx_dip_df + 129);

    auto tr_z_xz_xxz = pbuffer.data(idx_dip_df + 142);

    auto tr_z_xz_xyz = pbuffer.data(idx_dip_df + 144);

    auto tr_z_xz_xzz = pbuffer.data(idx_dip_df + 145);

    auto tr_z_xz_yyz = pbuffer.data(idx_dip_df + 147);

    auto tr_z_xz_yzz = pbuffer.data(idx_dip_df + 148);

    auto tr_z_xz_zzz = pbuffer.data(idx_dip_df + 149);

    auto tr_z_yy_xxx = pbuffer.data(idx_dip_df + 150);

    auto tr_z_yy_xxy = pbuffer.data(idx_dip_df + 151);

    auto tr_z_yy_xxz = pbuffer.data(idx_dip_df + 152);

    auto tr_z_yy_xyy = pbuffer.data(idx_dip_df + 153);

    auto tr_z_yy_xyz = pbuffer.data(idx_dip_df + 154);

    auto tr_z_yy_xzz = pbuffer.data(idx_dip_df + 155);

    auto tr_z_yy_yyy = pbuffer.data(idx_dip_df + 156);

    auto tr_z_yy_yyz = pbuffer.data(idx_dip_df + 157);

    auto tr_z_yy_yzz = pbuffer.data(idx_dip_df + 158);

    auto tr_z_yy_zzz = pbuffer.data(idx_dip_df + 159);

    auto tr_z_yz_xxy = pbuffer.data(idx_dip_df + 161);

    auto tr_z_yz_xxz = pbuffer.data(idx_dip_df + 162);

    auto tr_z_yz_xyy = pbuffer.data(idx_dip_df + 163);

    auto tr_z_yz_xyz = pbuffer.data(idx_dip_df + 164);

    auto tr_z_yz_xzz = pbuffer.data(idx_dip_df + 165);

    auto tr_z_yz_yyy = pbuffer.data(idx_dip_df + 166);

    auto tr_z_yz_yyz = pbuffer.data(idx_dip_df + 167);

    auto tr_z_yz_yzz = pbuffer.data(idx_dip_df + 168);

    auto tr_z_yz_zzz = pbuffer.data(idx_dip_df + 169);

    auto tr_z_zz_xxx = pbuffer.data(idx_dip_df + 170);

    auto tr_z_zz_xxy = pbuffer.data(idx_dip_df + 171);

    auto tr_z_zz_xxz = pbuffer.data(idx_dip_df + 172);

    auto tr_z_zz_xyy = pbuffer.data(idx_dip_df + 173);

    auto tr_z_zz_xyz = pbuffer.data(idx_dip_df + 174);

    auto tr_z_zz_xzz = pbuffer.data(idx_dip_df + 175);

    auto tr_z_zz_yyy = pbuffer.data(idx_dip_df + 176);

    auto tr_z_zz_yyz = pbuffer.data(idx_dip_df + 177);

    auto tr_z_zz_yzz = pbuffer.data(idx_dip_df + 178);

    auto tr_z_zz_zzz = pbuffer.data(idx_dip_df + 179);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_ovl_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_ovl_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_ovl_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_ovl_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_ovl_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_ovl_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_ovl_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_ovl_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_ovl_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_ovl_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_ovl_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_ovl_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_ovl_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_ovl_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_ovl_dg + 14);

    auto ts_yy_xxxx = pbuffer.data(idx_ovl_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_ovl_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_ovl_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_ovl_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_ovl_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_ovl_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_ovl_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_ovl_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_ovl_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_ovl_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_ovl_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_ovl_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_ovl_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_ovl_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_ovl_dg + 59);

    auto ts_yz_yyyz = pbuffer.data(idx_ovl_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_ovl_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_ovl_dg + 73);

    auto ts_zz_xxxx = pbuffer.data(idx_ovl_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_ovl_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_ovl_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_ovl_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_ovl_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_ovl_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_ovl_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_ovl_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_ovl_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_ovl_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_ovl_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_ovl_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_ovl_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_ovl_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_ovl_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto tr_x_xx_xxxx = pbuffer.data(idx_dip_dg);

    auto tr_x_xx_xxxy = pbuffer.data(idx_dip_dg + 1);

    auto tr_x_xx_xxxz = pbuffer.data(idx_dip_dg + 2);

    auto tr_x_xx_xxyy = pbuffer.data(idx_dip_dg + 3);

    auto tr_x_xx_xxyz = pbuffer.data(idx_dip_dg + 4);

    auto tr_x_xx_xxzz = pbuffer.data(idx_dip_dg + 5);

    auto tr_x_xx_xyyy = pbuffer.data(idx_dip_dg + 6);

    auto tr_x_xx_xyyz = pbuffer.data(idx_dip_dg + 7);

    auto tr_x_xx_xyzz = pbuffer.data(idx_dip_dg + 8);

    auto tr_x_xx_xzzz = pbuffer.data(idx_dip_dg + 9);

    auto tr_x_xx_yyyy = pbuffer.data(idx_dip_dg + 10);

    auto tr_x_xx_yyyz = pbuffer.data(idx_dip_dg + 11);

    auto tr_x_xx_yyzz = pbuffer.data(idx_dip_dg + 12);

    auto tr_x_xx_yzzz = pbuffer.data(idx_dip_dg + 13);

    auto tr_x_xx_zzzz = pbuffer.data(idx_dip_dg + 14);

    auto tr_x_xy_xxxx = pbuffer.data(idx_dip_dg + 15);

    auto tr_x_xy_xxxy = pbuffer.data(idx_dip_dg + 16);

    auto tr_x_xy_xxxz = pbuffer.data(idx_dip_dg + 17);

    auto tr_x_xy_xxyy = pbuffer.data(idx_dip_dg + 18);

    auto tr_x_xy_xxzz = pbuffer.data(idx_dip_dg + 20);

    auto tr_x_xy_xyyy = pbuffer.data(idx_dip_dg + 21);

    auto tr_x_xy_xzzz = pbuffer.data(idx_dip_dg + 24);

    auto tr_x_xy_yyyy = pbuffer.data(idx_dip_dg + 25);

    auto tr_x_xz_xxxx = pbuffer.data(idx_dip_dg + 30);

    auto tr_x_xz_xxxy = pbuffer.data(idx_dip_dg + 31);

    auto tr_x_xz_xxxz = pbuffer.data(idx_dip_dg + 32);

    auto tr_x_xz_xxyy = pbuffer.data(idx_dip_dg + 33);

    auto tr_x_xz_xxyz = pbuffer.data(idx_dip_dg + 34);

    auto tr_x_xz_xxzz = pbuffer.data(idx_dip_dg + 35);

    auto tr_x_xz_xyyy = pbuffer.data(idx_dip_dg + 36);

    auto tr_x_xz_xyyz = pbuffer.data(idx_dip_dg + 37);

    auto tr_x_xz_xyzz = pbuffer.data(idx_dip_dg + 38);

    auto tr_x_xz_xzzz = pbuffer.data(idx_dip_dg + 39);

    auto tr_x_xz_zzzz = pbuffer.data(idx_dip_dg + 44);

    auto tr_x_yy_xxxx = pbuffer.data(idx_dip_dg + 45);

    auto tr_x_yy_xxxy = pbuffer.data(idx_dip_dg + 46);

    auto tr_x_yy_xxxz = pbuffer.data(idx_dip_dg + 47);

    auto tr_x_yy_xxyy = pbuffer.data(idx_dip_dg + 48);

    auto tr_x_yy_xxyz = pbuffer.data(idx_dip_dg + 49);

    auto tr_x_yy_xxzz = pbuffer.data(idx_dip_dg + 50);

    auto tr_x_yy_xyyy = pbuffer.data(idx_dip_dg + 51);

    auto tr_x_yy_xyyz = pbuffer.data(idx_dip_dg + 52);

    auto tr_x_yy_xyzz = pbuffer.data(idx_dip_dg + 53);

    auto tr_x_yy_xzzz = pbuffer.data(idx_dip_dg + 54);

    auto tr_x_yy_yyyy = pbuffer.data(idx_dip_dg + 55);

    auto tr_x_yy_yyyz = pbuffer.data(idx_dip_dg + 56);

    auto tr_x_yy_yyzz = pbuffer.data(idx_dip_dg + 57);

    auto tr_x_yy_yzzz = pbuffer.data(idx_dip_dg + 58);

    auto tr_x_yy_zzzz = pbuffer.data(idx_dip_dg + 59);

    auto tr_x_yz_xxxz = pbuffer.data(idx_dip_dg + 62);

    auto tr_x_yz_xxzz = pbuffer.data(idx_dip_dg + 65);

    auto tr_x_yz_xzzz = pbuffer.data(idx_dip_dg + 69);

    auto tr_x_yz_yyyz = pbuffer.data(idx_dip_dg + 71);

    auto tr_x_yz_yyzz = pbuffer.data(idx_dip_dg + 72);

    auto tr_x_yz_yzzz = pbuffer.data(idx_dip_dg + 73);

    auto tr_x_yz_zzzz = pbuffer.data(idx_dip_dg + 74);

    auto tr_x_zz_xxxx = pbuffer.data(idx_dip_dg + 75);

    auto tr_x_zz_xxxy = pbuffer.data(idx_dip_dg + 76);

    auto tr_x_zz_xxxz = pbuffer.data(idx_dip_dg + 77);

    auto tr_x_zz_xxyy = pbuffer.data(idx_dip_dg + 78);

    auto tr_x_zz_xxyz = pbuffer.data(idx_dip_dg + 79);

    auto tr_x_zz_xxzz = pbuffer.data(idx_dip_dg + 80);

    auto tr_x_zz_xyyy = pbuffer.data(idx_dip_dg + 81);

    auto tr_x_zz_xyyz = pbuffer.data(idx_dip_dg + 82);

    auto tr_x_zz_xyzz = pbuffer.data(idx_dip_dg + 83);

    auto tr_x_zz_xzzz = pbuffer.data(idx_dip_dg + 84);

    auto tr_x_zz_yyyy = pbuffer.data(idx_dip_dg + 85);

    auto tr_x_zz_yyyz = pbuffer.data(idx_dip_dg + 86);

    auto tr_x_zz_yyzz = pbuffer.data(idx_dip_dg + 87);

    auto tr_x_zz_yzzz = pbuffer.data(idx_dip_dg + 88);

    auto tr_x_zz_zzzz = pbuffer.data(idx_dip_dg + 89);

    auto tr_y_xx_xxxx = pbuffer.data(idx_dip_dg + 90);

    auto tr_y_xx_xxxy = pbuffer.data(idx_dip_dg + 91);

    auto tr_y_xx_xxxz = pbuffer.data(idx_dip_dg + 92);

    auto tr_y_xx_xxyy = pbuffer.data(idx_dip_dg + 93);

    auto tr_y_xx_xxyz = pbuffer.data(idx_dip_dg + 94);

    auto tr_y_xx_xxzz = pbuffer.data(idx_dip_dg + 95);

    auto tr_y_xx_xyyy = pbuffer.data(idx_dip_dg + 96);

    auto tr_y_xx_xyyz = pbuffer.data(idx_dip_dg + 97);

    auto tr_y_xx_xyzz = pbuffer.data(idx_dip_dg + 98);

    auto tr_y_xx_xzzz = pbuffer.data(idx_dip_dg + 99);

    auto tr_y_xx_yyyy = pbuffer.data(idx_dip_dg + 100);

    auto tr_y_xx_yyyz = pbuffer.data(idx_dip_dg + 101);

    auto tr_y_xx_yyzz = pbuffer.data(idx_dip_dg + 102);

    auto tr_y_xx_yzzz = pbuffer.data(idx_dip_dg + 103);

    auto tr_y_xx_zzzz = pbuffer.data(idx_dip_dg + 104);

    auto tr_y_xy_xxxx = pbuffer.data(idx_dip_dg + 105);

    auto tr_y_xy_xxxy = pbuffer.data(idx_dip_dg + 106);

    auto tr_y_xy_xxyy = pbuffer.data(idx_dip_dg + 108);

    auto tr_y_xy_xxyz = pbuffer.data(idx_dip_dg + 109);

    auto tr_y_xy_xyyy = pbuffer.data(idx_dip_dg + 111);

    auto tr_y_xy_xyyz = pbuffer.data(idx_dip_dg + 112);

    auto tr_y_xy_xyzz = pbuffer.data(idx_dip_dg + 113);

    auto tr_y_xy_yyyy = pbuffer.data(idx_dip_dg + 115);

    auto tr_y_xy_yyyz = pbuffer.data(idx_dip_dg + 116);

    auto tr_y_xy_yyzz = pbuffer.data(idx_dip_dg + 117);

    auto tr_y_xy_yzzz = pbuffer.data(idx_dip_dg + 118);

    auto tr_y_xy_zzzz = pbuffer.data(idx_dip_dg + 119);

    auto tr_y_xz_yyyz = pbuffer.data(idx_dip_dg + 131);

    auto tr_y_xz_yyzz = pbuffer.data(idx_dip_dg + 132);

    auto tr_y_xz_yzzz = pbuffer.data(idx_dip_dg + 133);

    auto tr_y_xz_zzzz = pbuffer.data(idx_dip_dg + 134);

    auto tr_y_yy_xxxx = pbuffer.data(idx_dip_dg + 135);

    auto tr_y_yy_xxxy = pbuffer.data(idx_dip_dg + 136);

    auto tr_y_yy_xxxz = pbuffer.data(idx_dip_dg + 137);

    auto tr_y_yy_xxyy = pbuffer.data(idx_dip_dg + 138);

    auto tr_y_yy_xxyz = pbuffer.data(idx_dip_dg + 139);

    auto tr_y_yy_xxzz = pbuffer.data(idx_dip_dg + 140);

    auto tr_y_yy_xyyy = pbuffer.data(idx_dip_dg + 141);

    auto tr_y_yy_xyyz = pbuffer.data(idx_dip_dg + 142);

    auto tr_y_yy_xyzz = pbuffer.data(idx_dip_dg + 143);

    auto tr_y_yy_xzzz = pbuffer.data(idx_dip_dg + 144);

    auto tr_y_yy_yyyy = pbuffer.data(idx_dip_dg + 145);

    auto tr_y_yy_yyyz = pbuffer.data(idx_dip_dg + 146);

    auto tr_y_yy_yyzz = pbuffer.data(idx_dip_dg + 147);

    auto tr_y_yy_yzzz = pbuffer.data(idx_dip_dg + 148);

    auto tr_y_yy_zzzz = pbuffer.data(idx_dip_dg + 149);

    auto tr_y_yz_xxxy = pbuffer.data(idx_dip_dg + 151);

    auto tr_y_yz_xxxz = pbuffer.data(idx_dip_dg + 152);

    auto tr_y_yz_xxyy = pbuffer.data(idx_dip_dg + 153);

    auto tr_y_yz_xxyz = pbuffer.data(idx_dip_dg + 154);

    auto tr_y_yz_xxzz = pbuffer.data(idx_dip_dg + 155);

    auto tr_y_yz_xyyy = pbuffer.data(idx_dip_dg + 156);

    auto tr_y_yz_xyyz = pbuffer.data(idx_dip_dg + 157);

    auto tr_y_yz_xyzz = pbuffer.data(idx_dip_dg + 158);

    auto tr_y_yz_xzzz = pbuffer.data(idx_dip_dg + 159);

    auto tr_y_yz_yyyy = pbuffer.data(idx_dip_dg + 160);

    auto tr_y_yz_yyyz = pbuffer.data(idx_dip_dg + 161);

    auto tr_y_yz_yyzz = pbuffer.data(idx_dip_dg + 162);

    auto tr_y_yz_yzzz = pbuffer.data(idx_dip_dg + 163);

    auto tr_y_yz_zzzz = pbuffer.data(idx_dip_dg + 164);

    auto tr_y_zz_xxxx = pbuffer.data(idx_dip_dg + 165);

    auto tr_y_zz_xxxy = pbuffer.data(idx_dip_dg + 166);

    auto tr_y_zz_xxxz = pbuffer.data(idx_dip_dg + 167);

    auto tr_y_zz_xxyy = pbuffer.data(idx_dip_dg + 168);

    auto tr_y_zz_xxyz = pbuffer.data(idx_dip_dg + 169);

    auto tr_y_zz_xxzz = pbuffer.data(idx_dip_dg + 170);

    auto tr_y_zz_xyyy = pbuffer.data(idx_dip_dg + 171);

    auto tr_y_zz_xyyz = pbuffer.data(idx_dip_dg + 172);

    auto tr_y_zz_xyzz = pbuffer.data(idx_dip_dg + 173);

    auto tr_y_zz_xzzz = pbuffer.data(idx_dip_dg + 174);

    auto tr_y_zz_yyyy = pbuffer.data(idx_dip_dg + 175);

    auto tr_y_zz_yyyz = pbuffer.data(idx_dip_dg + 176);

    auto tr_y_zz_yyzz = pbuffer.data(idx_dip_dg + 177);

    auto tr_y_zz_yzzz = pbuffer.data(idx_dip_dg + 178);

    auto tr_y_zz_zzzz = pbuffer.data(idx_dip_dg + 179);

    auto tr_z_xx_xxxx = pbuffer.data(idx_dip_dg + 180);

    auto tr_z_xx_xxxy = pbuffer.data(idx_dip_dg + 181);

    auto tr_z_xx_xxxz = pbuffer.data(idx_dip_dg + 182);

    auto tr_z_xx_xxyy = pbuffer.data(idx_dip_dg + 183);

    auto tr_z_xx_xxyz = pbuffer.data(idx_dip_dg + 184);

    auto tr_z_xx_xxzz = pbuffer.data(idx_dip_dg + 185);

    auto tr_z_xx_xyyy = pbuffer.data(idx_dip_dg + 186);

    auto tr_z_xx_xyyz = pbuffer.data(idx_dip_dg + 187);

    auto tr_z_xx_xyzz = pbuffer.data(idx_dip_dg + 188);

    auto tr_z_xx_xzzz = pbuffer.data(idx_dip_dg + 189);

    auto tr_z_xx_yyyy = pbuffer.data(idx_dip_dg + 190);

    auto tr_z_xx_yyyz = pbuffer.data(idx_dip_dg + 191);

    auto tr_z_xx_yyzz = pbuffer.data(idx_dip_dg + 192);

    auto tr_z_xx_yzzz = pbuffer.data(idx_dip_dg + 193);

    auto tr_z_xx_zzzz = pbuffer.data(idx_dip_dg + 194);

    auto tr_z_xy_yyyy = pbuffer.data(idx_dip_dg + 205);

    auto tr_z_xy_yyyz = pbuffer.data(idx_dip_dg + 206);

    auto tr_z_xy_yyzz = pbuffer.data(idx_dip_dg + 207);

    auto tr_z_xy_yzzz = pbuffer.data(idx_dip_dg + 208);

    auto tr_z_xz_xxxx = pbuffer.data(idx_dip_dg + 210);

    auto tr_z_xz_xxxz = pbuffer.data(idx_dip_dg + 212);

    auto tr_z_xz_xxyz = pbuffer.data(idx_dip_dg + 214);

    auto tr_z_xz_xxzz = pbuffer.data(idx_dip_dg + 215);

    auto tr_z_xz_xyyz = pbuffer.data(idx_dip_dg + 217);

    auto tr_z_xz_xyzz = pbuffer.data(idx_dip_dg + 218);

    auto tr_z_xz_xzzz = pbuffer.data(idx_dip_dg + 219);

    auto tr_z_xz_yyyy = pbuffer.data(idx_dip_dg + 220);

    auto tr_z_xz_yyyz = pbuffer.data(idx_dip_dg + 221);

    auto tr_z_xz_yyzz = pbuffer.data(idx_dip_dg + 222);

    auto tr_z_xz_yzzz = pbuffer.data(idx_dip_dg + 223);

    auto tr_z_xz_zzzz = pbuffer.data(idx_dip_dg + 224);

    auto tr_z_yy_xxxx = pbuffer.data(idx_dip_dg + 225);

    auto tr_z_yy_xxxy = pbuffer.data(idx_dip_dg + 226);

    auto tr_z_yy_xxxz = pbuffer.data(idx_dip_dg + 227);

    auto tr_z_yy_xxyy = pbuffer.data(idx_dip_dg + 228);

    auto tr_z_yy_xxyz = pbuffer.data(idx_dip_dg + 229);

    auto tr_z_yy_xxzz = pbuffer.data(idx_dip_dg + 230);

    auto tr_z_yy_xyyy = pbuffer.data(idx_dip_dg + 231);

    auto tr_z_yy_xyyz = pbuffer.data(idx_dip_dg + 232);

    auto tr_z_yy_xyzz = pbuffer.data(idx_dip_dg + 233);

    auto tr_z_yy_xzzz = pbuffer.data(idx_dip_dg + 234);

    auto tr_z_yy_yyyy = pbuffer.data(idx_dip_dg + 235);

    auto tr_z_yy_yyyz = pbuffer.data(idx_dip_dg + 236);

    auto tr_z_yy_yyzz = pbuffer.data(idx_dip_dg + 237);

    auto tr_z_yy_yzzz = pbuffer.data(idx_dip_dg + 238);

    auto tr_z_yy_zzzz = pbuffer.data(idx_dip_dg + 239);

    auto tr_z_yz_xxxx = pbuffer.data(idx_dip_dg + 240);

    auto tr_z_yz_xxxy = pbuffer.data(idx_dip_dg + 241);

    auto tr_z_yz_xxxz = pbuffer.data(idx_dip_dg + 242);

    auto tr_z_yz_xxyy = pbuffer.data(idx_dip_dg + 243);

    auto tr_z_yz_xxyz = pbuffer.data(idx_dip_dg + 244);

    auto tr_z_yz_xxzz = pbuffer.data(idx_dip_dg + 245);

    auto tr_z_yz_xyyy = pbuffer.data(idx_dip_dg + 246);

    auto tr_z_yz_xyyz = pbuffer.data(idx_dip_dg + 247);

    auto tr_z_yz_xyzz = pbuffer.data(idx_dip_dg + 248);

    auto tr_z_yz_xzzz = pbuffer.data(idx_dip_dg + 249);

    auto tr_z_yz_yyyy = pbuffer.data(idx_dip_dg + 250);

    auto tr_z_yz_yyyz = pbuffer.data(idx_dip_dg + 251);

    auto tr_z_yz_yyzz = pbuffer.data(idx_dip_dg + 252);

    auto tr_z_yz_yzzz = pbuffer.data(idx_dip_dg + 253);

    auto tr_z_yz_zzzz = pbuffer.data(idx_dip_dg + 254);

    auto tr_z_zz_xxxx = pbuffer.data(idx_dip_dg + 255);

    auto tr_z_zz_xxxy = pbuffer.data(idx_dip_dg + 256);

    auto tr_z_zz_xxxz = pbuffer.data(idx_dip_dg + 257);

    auto tr_z_zz_xxyy = pbuffer.data(idx_dip_dg + 258);

    auto tr_z_zz_xxyz = pbuffer.data(idx_dip_dg + 259);

    auto tr_z_zz_xxzz = pbuffer.data(idx_dip_dg + 260);

    auto tr_z_zz_xyyy = pbuffer.data(idx_dip_dg + 261);

    auto tr_z_zz_xyyz = pbuffer.data(idx_dip_dg + 262);

    auto tr_z_zz_xyzz = pbuffer.data(idx_dip_dg + 263);

    auto tr_z_zz_xzzz = pbuffer.data(idx_dip_dg + 264);

    auto tr_z_zz_yyyy = pbuffer.data(idx_dip_dg + 265);

    auto tr_z_zz_yyyz = pbuffer.data(idx_dip_dg + 266);

    auto tr_z_zz_yyzz = pbuffer.data(idx_dip_dg + 267);

    auto tr_z_zz_yzzz = pbuffer.data(idx_dip_dg + 268);

    auto tr_z_zz_zzzz = pbuffer.data(idx_dip_dg + 269);

    // Set up 0-15 components of targeted buffer : FG

    auto tr_x_xxx_xxxx = pbuffer.data(idx_dip_fg);

    auto tr_x_xxx_xxxy = pbuffer.data(idx_dip_fg + 1);

    auto tr_x_xxx_xxxz = pbuffer.data(idx_dip_fg + 2);

    auto tr_x_xxx_xxyy = pbuffer.data(idx_dip_fg + 3);

    auto tr_x_xxx_xxyz = pbuffer.data(idx_dip_fg + 4);

    auto tr_x_xxx_xxzz = pbuffer.data(idx_dip_fg + 5);

    auto tr_x_xxx_xyyy = pbuffer.data(idx_dip_fg + 6);

    auto tr_x_xxx_xyyz = pbuffer.data(idx_dip_fg + 7);

    auto tr_x_xxx_xyzz = pbuffer.data(idx_dip_fg + 8);

    auto tr_x_xxx_xzzz = pbuffer.data(idx_dip_fg + 9);

    auto tr_x_xxx_yyyy = pbuffer.data(idx_dip_fg + 10);

    auto tr_x_xxx_yyyz = pbuffer.data(idx_dip_fg + 11);

    auto tr_x_xxx_yyzz = pbuffer.data(idx_dip_fg + 12);

    auto tr_x_xxx_yzzz = pbuffer.data(idx_dip_fg + 13);

    auto tr_x_xxx_zzzz = pbuffer.data(idx_dip_fg + 14);

#pragma omp simd aligned(pa_x,              \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxy,   \
                             tr_x_x_xxxz,   \
                             tr_x_x_xxyy,   \
                             tr_x_x_xxyz,   \
                             tr_x_x_xxzz,   \
                             tr_x_x_xyyy,   \
                             tr_x_x_xyyz,   \
                             tr_x_x_xyzz,   \
                             tr_x_x_xzzz,   \
                             tr_x_x_yyyy,   \
                             tr_x_x_yyyz,   \
                             tr_x_x_yyzz,   \
                             tr_x_x_yzzz,   \
                             tr_x_x_zzzz,   \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxxx,  \
                             tr_x_xx_xxxy,  \
                             tr_x_xx_xxxz,  \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxyy,  \
                             tr_x_xx_xxyz,  \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xxzz,  \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyyy,  \
                             tr_x_xx_xyyz,  \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xyzz,  \
                             tr_x_xx_xzz,   \
                             tr_x_xx_xzzz,  \
                             tr_x_xx_yyy,   \
                             tr_x_xx_yyyy,  \
                             tr_x_xx_yyyz,  \
                             tr_x_xx_yyz,   \
                             tr_x_xx_yyzz,  \
                             tr_x_xx_yzz,   \
                             tr_x_xx_yzzz,  \
                             tr_x_xx_zzz,   \
                             tr_x_xx_zzzz,  \
                             tr_x_xxx_xxxx, \
                             tr_x_xxx_xxxy, \
                             tr_x_xxx_xxxz, \
                             tr_x_xxx_xxyy, \
                             tr_x_xxx_xxyz, \
                             tr_x_xxx_xxzz, \
                             tr_x_xxx_xyyy, \
                             tr_x_xxx_xyyz, \
                             tr_x_xxx_xyzz, \
                             tr_x_xxx_xzzz, \
                             tr_x_xxx_yyyy, \
                             tr_x_xxx_yyyz, \
                             tr_x_xxx_yyzz, \
                             tr_x_xxx_yzzz, \
                             tr_x_xxx_zzzz, \
                             ts_xx_xxxx,    \
                             ts_xx_xxxy,    \
                             ts_xx_xxxz,    \
                             ts_xx_xxyy,    \
                             ts_xx_xxyz,    \
                             ts_xx_xxzz,    \
                             ts_xx_xyyy,    \
                             ts_xx_xyyz,    \
                             ts_xx_xyzz,    \
                             ts_xx_xzzz,    \
                             ts_xx_yyyy,    \
                             ts_xx_yyyz,    \
                             ts_xx_yyzz,    \
                             ts_xx_yzzz,    \
                             ts_xx_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_xxxx[i] = 2.0 * tr_x_x_xxxx[i] * fe_0 + 4.0 * tr_x_xx_xxx[i] * fe_0 + ts_xx_xxxx[i] * fe_0 + tr_x_xx_xxxx[i] * pa_x[i];

        tr_x_xxx_xxxy[i] = 2.0 * tr_x_x_xxxy[i] * fe_0 + 3.0 * tr_x_xx_xxy[i] * fe_0 + ts_xx_xxxy[i] * fe_0 + tr_x_xx_xxxy[i] * pa_x[i];

        tr_x_xxx_xxxz[i] = 2.0 * tr_x_x_xxxz[i] * fe_0 + 3.0 * tr_x_xx_xxz[i] * fe_0 + ts_xx_xxxz[i] * fe_0 + tr_x_xx_xxxz[i] * pa_x[i];

        tr_x_xxx_xxyy[i] = 2.0 * tr_x_x_xxyy[i] * fe_0 + 2.0 * tr_x_xx_xyy[i] * fe_0 + ts_xx_xxyy[i] * fe_0 + tr_x_xx_xxyy[i] * pa_x[i];

        tr_x_xxx_xxyz[i] = 2.0 * tr_x_x_xxyz[i] * fe_0 + 2.0 * tr_x_xx_xyz[i] * fe_0 + ts_xx_xxyz[i] * fe_0 + tr_x_xx_xxyz[i] * pa_x[i];

        tr_x_xxx_xxzz[i] = 2.0 * tr_x_x_xxzz[i] * fe_0 + 2.0 * tr_x_xx_xzz[i] * fe_0 + ts_xx_xxzz[i] * fe_0 + tr_x_xx_xxzz[i] * pa_x[i];

        tr_x_xxx_xyyy[i] = 2.0 * tr_x_x_xyyy[i] * fe_0 + tr_x_xx_yyy[i] * fe_0 + ts_xx_xyyy[i] * fe_0 + tr_x_xx_xyyy[i] * pa_x[i];

        tr_x_xxx_xyyz[i] = 2.0 * tr_x_x_xyyz[i] * fe_0 + tr_x_xx_yyz[i] * fe_0 + ts_xx_xyyz[i] * fe_0 + tr_x_xx_xyyz[i] * pa_x[i];

        tr_x_xxx_xyzz[i] = 2.0 * tr_x_x_xyzz[i] * fe_0 + tr_x_xx_yzz[i] * fe_0 + ts_xx_xyzz[i] * fe_0 + tr_x_xx_xyzz[i] * pa_x[i];

        tr_x_xxx_xzzz[i] = 2.0 * tr_x_x_xzzz[i] * fe_0 + tr_x_xx_zzz[i] * fe_0 + ts_xx_xzzz[i] * fe_0 + tr_x_xx_xzzz[i] * pa_x[i];

        tr_x_xxx_yyyy[i] = 2.0 * tr_x_x_yyyy[i] * fe_0 + ts_xx_yyyy[i] * fe_0 + tr_x_xx_yyyy[i] * pa_x[i];

        tr_x_xxx_yyyz[i] = 2.0 * tr_x_x_yyyz[i] * fe_0 + ts_xx_yyyz[i] * fe_0 + tr_x_xx_yyyz[i] * pa_x[i];

        tr_x_xxx_yyzz[i] = 2.0 * tr_x_x_yyzz[i] * fe_0 + ts_xx_yyzz[i] * fe_0 + tr_x_xx_yyzz[i] * pa_x[i];

        tr_x_xxx_yzzz[i] = 2.0 * tr_x_x_yzzz[i] * fe_0 + ts_xx_yzzz[i] * fe_0 + tr_x_xx_yzzz[i] * pa_x[i];

        tr_x_xxx_zzzz[i] = 2.0 * tr_x_x_zzzz[i] * fe_0 + ts_xx_zzzz[i] * fe_0 + tr_x_xx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto tr_x_xxy_xxxx = pbuffer.data(idx_dip_fg + 15);

    auto tr_x_xxy_xxxy = pbuffer.data(idx_dip_fg + 16);

    auto tr_x_xxy_xxxz = pbuffer.data(idx_dip_fg + 17);

    auto tr_x_xxy_xxyy = pbuffer.data(idx_dip_fg + 18);

    auto tr_x_xxy_xxyz = pbuffer.data(idx_dip_fg + 19);

    auto tr_x_xxy_xxzz = pbuffer.data(idx_dip_fg + 20);

    auto tr_x_xxy_xyyy = pbuffer.data(idx_dip_fg + 21);

    auto tr_x_xxy_xyyz = pbuffer.data(idx_dip_fg + 22);

    auto tr_x_xxy_xyzz = pbuffer.data(idx_dip_fg + 23);

    auto tr_x_xxy_xzzz = pbuffer.data(idx_dip_fg + 24);

    auto tr_x_xxy_yyyy = pbuffer.data(idx_dip_fg + 25);

    auto tr_x_xxy_yyyz = pbuffer.data(idx_dip_fg + 26);

    auto tr_x_xxy_yyzz = pbuffer.data(idx_dip_fg + 27);

    auto tr_x_xxy_yzzz = pbuffer.data(idx_dip_fg + 28);

    auto tr_x_xxy_zzzz = pbuffer.data(idx_dip_fg + 29);

#pragma omp simd aligned(pa_y,              \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxxx,  \
                             tr_x_xx_xxxy,  \
                             tr_x_xx_xxxz,  \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxyy,  \
                             tr_x_xx_xxyz,  \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xxzz,  \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyyy,  \
                             tr_x_xx_xyyz,  \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xyzz,  \
                             tr_x_xx_xzz,   \
                             tr_x_xx_xzzz,  \
                             tr_x_xx_yyy,   \
                             tr_x_xx_yyyy,  \
                             tr_x_xx_yyyz,  \
                             tr_x_xx_yyz,   \
                             tr_x_xx_yyzz,  \
                             tr_x_xx_yzz,   \
                             tr_x_xx_yzzz,  \
                             tr_x_xx_zzz,   \
                             tr_x_xx_zzzz,  \
                             tr_x_xxy_xxxx, \
                             tr_x_xxy_xxxy, \
                             tr_x_xxy_xxxz, \
                             tr_x_xxy_xxyy, \
                             tr_x_xxy_xxyz, \
                             tr_x_xxy_xxzz, \
                             tr_x_xxy_xyyy, \
                             tr_x_xxy_xyyz, \
                             tr_x_xxy_xyzz, \
                             tr_x_xxy_xzzz, \
                             tr_x_xxy_yyyy, \
                             tr_x_xxy_yyyz, \
                             tr_x_xxy_yyzz, \
                             tr_x_xxy_yzzz, \
                             tr_x_xxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_xxxx[i] = tr_x_xx_xxxx[i] * pa_y[i];

        tr_x_xxy_xxxy[i] = tr_x_xx_xxx[i] * fe_0 + tr_x_xx_xxxy[i] * pa_y[i];

        tr_x_xxy_xxxz[i] = tr_x_xx_xxxz[i] * pa_y[i];

        tr_x_xxy_xxyy[i] = 2.0 * tr_x_xx_xxy[i] * fe_0 + tr_x_xx_xxyy[i] * pa_y[i];

        tr_x_xxy_xxyz[i] = tr_x_xx_xxz[i] * fe_0 + tr_x_xx_xxyz[i] * pa_y[i];

        tr_x_xxy_xxzz[i] = tr_x_xx_xxzz[i] * pa_y[i];

        tr_x_xxy_xyyy[i] = 3.0 * tr_x_xx_xyy[i] * fe_0 + tr_x_xx_xyyy[i] * pa_y[i];

        tr_x_xxy_xyyz[i] = 2.0 * tr_x_xx_xyz[i] * fe_0 + tr_x_xx_xyyz[i] * pa_y[i];

        tr_x_xxy_xyzz[i] = tr_x_xx_xzz[i] * fe_0 + tr_x_xx_xyzz[i] * pa_y[i];

        tr_x_xxy_xzzz[i] = tr_x_xx_xzzz[i] * pa_y[i];

        tr_x_xxy_yyyy[i] = 4.0 * tr_x_xx_yyy[i] * fe_0 + tr_x_xx_yyyy[i] * pa_y[i];

        tr_x_xxy_yyyz[i] = 3.0 * tr_x_xx_yyz[i] * fe_0 + tr_x_xx_yyyz[i] * pa_y[i];

        tr_x_xxy_yyzz[i] = 2.0 * tr_x_xx_yzz[i] * fe_0 + tr_x_xx_yyzz[i] * pa_y[i];

        tr_x_xxy_yzzz[i] = tr_x_xx_zzz[i] * fe_0 + tr_x_xx_yzzz[i] * pa_y[i];

        tr_x_xxy_zzzz[i] = tr_x_xx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto tr_x_xxz_xxxx = pbuffer.data(idx_dip_fg + 30);

    auto tr_x_xxz_xxxy = pbuffer.data(idx_dip_fg + 31);

    auto tr_x_xxz_xxxz = pbuffer.data(idx_dip_fg + 32);

    auto tr_x_xxz_xxyy = pbuffer.data(idx_dip_fg + 33);

    auto tr_x_xxz_xxyz = pbuffer.data(idx_dip_fg + 34);

    auto tr_x_xxz_xxzz = pbuffer.data(idx_dip_fg + 35);

    auto tr_x_xxz_xyyy = pbuffer.data(idx_dip_fg + 36);

    auto tr_x_xxz_xyyz = pbuffer.data(idx_dip_fg + 37);

    auto tr_x_xxz_xyzz = pbuffer.data(idx_dip_fg + 38);

    auto tr_x_xxz_xzzz = pbuffer.data(idx_dip_fg + 39);

    auto tr_x_xxz_yyyy = pbuffer.data(idx_dip_fg + 40);

    auto tr_x_xxz_yyyz = pbuffer.data(idx_dip_fg + 41);

    auto tr_x_xxz_yyzz = pbuffer.data(idx_dip_fg + 42);

    auto tr_x_xxz_yzzz = pbuffer.data(idx_dip_fg + 43);

    auto tr_x_xxz_zzzz = pbuffer.data(idx_dip_fg + 44);

#pragma omp simd aligned(pa_z,              \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxxx,  \
                             tr_x_xx_xxxy,  \
                             tr_x_xx_xxxz,  \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxyy,  \
                             tr_x_xx_xxyz,  \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xxzz,  \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyyy,  \
                             tr_x_xx_xyyz,  \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xyzz,  \
                             tr_x_xx_xzz,   \
                             tr_x_xx_xzzz,  \
                             tr_x_xx_yyy,   \
                             tr_x_xx_yyyy,  \
                             tr_x_xx_yyyz,  \
                             tr_x_xx_yyz,   \
                             tr_x_xx_yyzz,  \
                             tr_x_xx_yzz,   \
                             tr_x_xx_yzzz,  \
                             tr_x_xx_zzz,   \
                             tr_x_xx_zzzz,  \
                             tr_x_xxz_xxxx, \
                             tr_x_xxz_xxxy, \
                             tr_x_xxz_xxxz, \
                             tr_x_xxz_xxyy, \
                             tr_x_xxz_xxyz, \
                             tr_x_xxz_xxzz, \
                             tr_x_xxz_xyyy, \
                             tr_x_xxz_xyyz, \
                             tr_x_xxz_xyzz, \
                             tr_x_xxz_xzzz, \
                             tr_x_xxz_yyyy, \
                             tr_x_xxz_yyyz, \
                             tr_x_xxz_yyzz, \
                             tr_x_xxz_yzzz, \
                             tr_x_xxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_xxxx[i] = tr_x_xx_xxxx[i] * pa_z[i];

        tr_x_xxz_xxxy[i] = tr_x_xx_xxxy[i] * pa_z[i];

        tr_x_xxz_xxxz[i] = tr_x_xx_xxx[i] * fe_0 + tr_x_xx_xxxz[i] * pa_z[i];

        tr_x_xxz_xxyy[i] = tr_x_xx_xxyy[i] * pa_z[i];

        tr_x_xxz_xxyz[i] = tr_x_xx_xxy[i] * fe_0 + tr_x_xx_xxyz[i] * pa_z[i];

        tr_x_xxz_xxzz[i] = 2.0 * tr_x_xx_xxz[i] * fe_0 + tr_x_xx_xxzz[i] * pa_z[i];

        tr_x_xxz_xyyy[i] = tr_x_xx_xyyy[i] * pa_z[i];

        tr_x_xxz_xyyz[i] = tr_x_xx_xyy[i] * fe_0 + tr_x_xx_xyyz[i] * pa_z[i];

        tr_x_xxz_xyzz[i] = 2.0 * tr_x_xx_xyz[i] * fe_0 + tr_x_xx_xyzz[i] * pa_z[i];

        tr_x_xxz_xzzz[i] = 3.0 * tr_x_xx_xzz[i] * fe_0 + tr_x_xx_xzzz[i] * pa_z[i];

        tr_x_xxz_yyyy[i] = tr_x_xx_yyyy[i] * pa_z[i];

        tr_x_xxz_yyyz[i] = tr_x_xx_yyy[i] * fe_0 + tr_x_xx_yyyz[i] * pa_z[i];

        tr_x_xxz_yyzz[i] = 2.0 * tr_x_xx_yyz[i] * fe_0 + tr_x_xx_yyzz[i] * pa_z[i];

        tr_x_xxz_yzzz[i] = 3.0 * tr_x_xx_yzz[i] * fe_0 + tr_x_xx_yzzz[i] * pa_z[i];

        tr_x_xxz_zzzz[i] = 4.0 * tr_x_xx_zzz[i] * fe_0 + tr_x_xx_zzzz[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto tr_x_xyy_xxxx = pbuffer.data(idx_dip_fg + 45);

    auto tr_x_xyy_xxxy = pbuffer.data(idx_dip_fg + 46);

    auto tr_x_xyy_xxxz = pbuffer.data(idx_dip_fg + 47);

    auto tr_x_xyy_xxyy = pbuffer.data(idx_dip_fg + 48);

    auto tr_x_xyy_xxyz = pbuffer.data(idx_dip_fg + 49);

    auto tr_x_xyy_xxzz = pbuffer.data(idx_dip_fg + 50);

    auto tr_x_xyy_xyyy = pbuffer.data(idx_dip_fg + 51);

    auto tr_x_xyy_xyyz = pbuffer.data(idx_dip_fg + 52);

    auto tr_x_xyy_xyzz = pbuffer.data(idx_dip_fg + 53);

    auto tr_x_xyy_xzzz = pbuffer.data(idx_dip_fg + 54);

    auto tr_x_xyy_yyyy = pbuffer.data(idx_dip_fg + 55);

    auto tr_x_xyy_yyyz = pbuffer.data(idx_dip_fg + 56);

    auto tr_x_xyy_yyzz = pbuffer.data(idx_dip_fg + 57);

    auto tr_x_xyy_yzzz = pbuffer.data(idx_dip_fg + 58);

    auto tr_x_xyy_zzzz = pbuffer.data(idx_dip_fg + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxz,   \
                             tr_x_x_xxzz,   \
                             tr_x_x_xzzz,   \
                             tr_x_xy_xxxx,  \
                             tr_x_xy_xxxz,  \
                             tr_x_xy_xxzz,  \
                             tr_x_xy_xzzz,  \
                             tr_x_xyy_xxxx, \
                             tr_x_xyy_xxxy, \
                             tr_x_xyy_xxxz, \
                             tr_x_xyy_xxyy, \
                             tr_x_xyy_xxyz, \
                             tr_x_xyy_xxzz, \
                             tr_x_xyy_xyyy, \
                             tr_x_xyy_xyyz, \
                             tr_x_xyy_xyzz, \
                             tr_x_xyy_xzzz, \
                             tr_x_xyy_yyyy, \
                             tr_x_xyy_yyyz, \
                             tr_x_xyy_yyzz, \
                             tr_x_xyy_yzzz, \
                             tr_x_xyy_zzzz, \
                             tr_x_yy_xxxy,  \
                             tr_x_yy_xxy,   \
                             tr_x_yy_xxyy,  \
                             tr_x_yy_xxyz,  \
                             tr_x_yy_xyy,   \
                             tr_x_yy_xyyy,  \
                             tr_x_yy_xyyz,  \
                             tr_x_yy_xyz,   \
                             tr_x_yy_xyzz,  \
                             tr_x_yy_yyy,   \
                             tr_x_yy_yyyy,  \
                             tr_x_yy_yyyz,  \
                             tr_x_yy_yyz,   \
                             tr_x_yy_yyzz,  \
                             tr_x_yy_yzz,   \
                             tr_x_yy_yzzz,  \
                             tr_x_yy_zzzz,  \
                             ts_yy_xxxy,    \
                             ts_yy_xxyy,    \
                             ts_yy_xxyz,    \
                             ts_yy_xyyy,    \
                             ts_yy_xyyz,    \
                             ts_yy_xyzz,    \
                             ts_yy_yyyy,    \
                             ts_yy_yyyz,    \
                             ts_yy_yyzz,    \
                             ts_yy_yzzz,    \
                             ts_yy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_xxxx[i] = tr_x_x_xxxx[i] * fe_0 + tr_x_xy_xxxx[i] * pa_y[i];

        tr_x_xyy_xxxy[i] = 3.0 * tr_x_yy_xxy[i] * fe_0 + ts_yy_xxxy[i] * fe_0 + tr_x_yy_xxxy[i] * pa_x[i];

        tr_x_xyy_xxxz[i] = tr_x_x_xxxz[i] * fe_0 + tr_x_xy_xxxz[i] * pa_y[i];

        tr_x_xyy_xxyy[i] = 2.0 * tr_x_yy_xyy[i] * fe_0 + ts_yy_xxyy[i] * fe_0 + tr_x_yy_xxyy[i] * pa_x[i];

        tr_x_xyy_xxyz[i] = 2.0 * tr_x_yy_xyz[i] * fe_0 + ts_yy_xxyz[i] * fe_0 + tr_x_yy_xxyz[i] * pa_x[i];

        tr_x_xyy_xxzz[i] = tr_x_x_xxzz[i] * fe_0 + tr_x_xy_xxzz[i] * pa_y[i];

        tr_x_xyy_xyyy[i] = tr_x_yy_yyy[i] * fe_0 + ts_yy_xyyy[i] * fe_0 + tr_x_yy_xyyy[i] * pa_x[i];

        tr_x_xyy_xyyz[i] = tr_x_yy_yyz[i] * fe_0 + ts_yy_xyyz[i] * fe_0 + tr_x_yy_xyyz[i] * pa_x[i];

        tr_x_xyy_xyzz[i] = tr_x_yy_yzz[i] * fe_0 + ts_yy_xyzz[i] * fe_0 + tr_x_yy_xyzz[i] * pa_x[i];

        tr_x_xyy_xzzz[i] = tr_x_x_xzzz[i] * fe_0 + tr_x_xy_xzzz[i] * pa_y[i];

        tr_x_xyy_yyyy[i] = ts_yy_yyyy[i] * fe_0 + tr_x_yy_yyyy[i] * pa_x[i];

        tr_x_xyy_yyyz[i] = ts_yy_yyyz[i] * fe_0 + tr_x_yy_yyyz[i] * pa_x[i];

        tr_x_xyy_yyzz[i] = ts_yy_yyzz[i] * fe_0 + tr_x_yy_yyzz[i] * pa_x[i];

        tr_x_xyy_yzzz[i] = ts_yy_yzzz[i] * fe_0 + tr_x_yy_yzzz[i] * pa_x[i];

        tr_x_xyy_zzzz[i] = ts_yy_zzzz[i] * fe_0 + tr_x_yy_zzzz[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto tr_x_xyz_xxxx = pbuffer.data(idx_dip_fg + 60);

    auto tr_x_xyz_xxxy = pbuffer.data(idx_dip_fg + 61);

    auto tr_x_xyz_xxxz = pbuffer.data(idx_dip_fg + 62);

    auto tr_x_xyz_xxyy = pbuffer.data(idx_dip_fg + 63);

    auto tr_x_xyz_xxyz = pbuffer.data(idx_dip_fg + 64);

    auto tr_x_xyz_xxzz = pbuffer.data(idx_dip_fg + 65);

    auto tr_x_xyz_xyyy = pbuffer.data(idx_dip_fg + 66);

    auto tr_x_xyz_xyyz = pbuffer.data(idx_dip_fg + 67);

    auto tr_x_xyz_xyzz = pbuffer.data(idx_dip_fg + 68);

    auto tr_x_xyz_xzzz = pbuffer.data(idx_dip_fg + 69);

    auto tr_x_xyz_yyyy = pbuffer.data(idx_dip_fg + 70);

    auto tr_x_xyz_yyyz = pbuffer.data(idx_dip_fg + 71);

    auto tr_x_xyz_yyzz = pbuffer.data(idx_dip_fg + 72);

    auto tr_x_xyz_yzzz = pbuffer.data(idx_dip_fg + 73);

    auto tr_x_xyz_zzzz = pbuffer.data(idx_dip_fg + 74);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             tr_x_xy_xxxy,  \
                             tr_x_xy_xxyy,  \
                             tr_x_xy_xyyy,  \
                             tr_x_xy_yyyy,  \
                             tr_x_xyz_xxxx, \
                             tr_x_xyz_xxxy, \
                             tr_x_xyz_xxxz, \
                             tr_x_xyz_xxyy, \
                             tr_x_xyz_xxyz, \
                             tr_x_xyz_xxzz, \
                             tr_x_xyz_xyyy, \
                             tr_x_xyz_xyyz, \
                             tr_x_xyz_xyzz, \
                             tr_x_xyz_xzzz, \
                             tr_x_xyz_yyyy, \
                             tr_x_xyz_yyyz, \
                             tr_x_xyz_yyzz, \
                             tr_x_xyz_yzzz, \
                             tr_x_xyz_zzzz, \
                             tr_x_xz_xxxx,  \
                             tr_x_xz_xxxz,  \
                             tr_x_xz_xxyz,  \
                             tr_x_xz_xxz,   \
                             tr_x_xz_xxzz,  \
                             tr_x_xz_xyyz,  \
                             tr_x_xz_xyz,   \
                             tr_x_xz_xyzz,  \
                             tr_x_xz_xzz,   \
                             tr_x_xz_xzzz,  \
                             tr_x_xz_zzzz,  \
                             tr_x_yz_yyyz,  \
                             tr_x_yz_yyzz,  \
                             tr_x_yz_yzzz,  \
                             ts_yz_yyyz,    \
                             ts_yz_yyzz,    \
                             ts_yz_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyz_xxxx[i] = tr_x_xz_xxxx[i] * pa_y[i];

        tr_x_xyz_xxxy[i] = tr_x_xy_xxxy[i] * pa_z[i];

        tr_x_xyz_xxxz[i] = tr_x_xz_xxxz[i] * pa_y[i];

        tr_x_xyz_xxyy[i] = tr_x_xy_xxyy[i] * pa_z[i];

        tr_x_xyz_xxyz[i] = tr_x_xz_xxz[i] * fe_0 + tr_x_xz_xxyz[i] * pa_y[i];

        tr_x_xyz_xxzz[i] = tr_x_xz_xxzz[i] * pa_y[i];

        tr_x_xyz_xyyy[i] = tr_x_xy_xyyy[i] * pa_z[i];

        tr_x_xyz_xyyz[i] = 2.0 * tr_x_xz_xyz[i] * fe_0 + tr_x_xz_xyyz[i] * pa_y[i];

        tr_x_xyz_xyzz[i] = tr_x_xz_xzz[i] * fe_0 + tr_x_xz_xyzz[i] * pa_y[i];

        tr_x_xyz_xzzz[i] = tr_x_xz_xzzz[i] * pa_y[i];

        tr_x_xyz_yyyy[i] = tr_x_xy_yyyy[i] * pa_z[i];

        tr_x_xyz_yyyz[i] = ts_yz_yyyz[i] * fe_0 + tr_x_yz_yyyz[i] * pa_x[i];

        tr_x_xyz_yyzz[i] = ts_yz_yyzz[i] * fe_0 + tr_x_yz_yyzz[i] * pa_x[i];

        tr_x_xyz_yzzz[i] = ts_yz_yzzz[i] * fe_0 + tr_x_yz_yzzz[i] * pa_x[i];

        tr_x_xyz_zzzz[i] = tr_x_xz_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto tr_x_xzz_xxxx = pbuffer.data(idx_dip_fg + 75);

    auto tr_x_xzz_xxxy = pbuffer.data(idx_dip_fg + 76);

    auto tr_x_xzz_xxxz = pbuffer.data(idx_dip_fg + 77);

    auto tr_x_xzz_xxyy = pbuffer.data(idx_dip_fg + 78);

    auto tr_x_xzz_xxyz = pbuffer.data(idx_dip_fg + 79);

    auto tr_x_xzz_xxzz = pbuffer.data(idx_dip_fg + 80);

    auto tr_x_xzz_xyyy = pbuffer.data(idx_dip_fg + 81);

    auto tr_x_xzz_xyyz = pbuffer.data(idx_dip_fg + 82);

    auto tr_x_xzz_xyzz = pbuffer.data(idx_dip_fg + 83);

    auto tr_x_xzz_xzzz = pbuffer.data(idx_dip_fg + 84);

    auto tr_x_xzz_yyyy = pbuffer.data(idx_dip_fg + 85);

    auto tr_x_xzz_yyyz = pbuffer.data(idx_dip_fg + 86);

    auto tr_x_xzz_yyzz = pbuffer.data(idx_dip_fg + 87);

    auto tr_x_xzz_yzzz = pbuffer.data(idx_dip_fg + 88);

    auto tr_x_xzz_zzzz = pbuffer.data(idx_dip_fg + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxy,   \
                             tr_x_x_xxyy,   \
                             tr_x_x_xyyy,   \
                             tr_x_xz_xxxx,  \
                             tr_x_xz_xxxy,  \
                             tr_x_xz_xxyy,  \
                             tr_x_xz_xyyy,  \
                             tr_x_xzz_xxxx, \
                             tr_x_xzz_xxxy, \
                             tr_x_xzz_xxxz, \
                             tr_x_xzz_xxyy, \
                             tr_x_xzz_xxyz, \
                             tr_x_xzz_xxzz, \
                             tr_x_xzz_xyyy, \
                             tr_x_xzz_xyyz, \
                             tr_x_xzz_xyzz, \
                             tr_x_xzz_xzzz, \
                             tr_x_xzz_yyyy, \
                             tr_x_xzz_yyyz, \
                             tr_x_xzz_yyzz, \
                             tr_x_xzz_yzzz, \
                             tr_x_xzz_zzzz, \
                             tr_x_zz_xxxz,  \
                             tr_x_zz_xxyz,  \
                             tr_x_zz_xxz,   \
                             tr_x_zz_xxzz,  \
                             tr_x_zz_xyyz,  \
                             tr_x_zz_xyz,   \
                             tr_x_zz_xyzz,  \
                             tr_x_zz_xzz,   \
                             tr_x_zz_xzzz,  \
                             tr_x_zz_yyyy,  \
                             tr_x_zz_yyyz,  \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yyzz,  \
                             tr_x_zz_yzz,   \
                             tr_x_zz_yzzz,  \
                             tr_x_zz_zzz,   \
                             tr_x_zz_zzzz,  \
                             ts_zz_xxxz,    \
                             ts_zz_xxyz,    \
                             ts_zz_xxzz,    \
                             ts_zz_xyyz,    \
                             ts_zz_xyzz,    \
                             ts_zz_xzzz,    \
                             ts_zz_yyyy,    \
                             ts_zz_yyyz,    \
                             ts_zz_yyzz,    \
                             ts_zz_yzzz,    \
                             ts_zz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_xxxx[i] = tr_x_x_xxxx[i] * fe_0 + tr_x_xz_xxxx[i] * pa_z[i];

        tr_x_xzz_xxxy[i] = tr_x_x_xxxy[i] * fe_0 + tr_x_xz_xxxy[i] * pa_z[i];

        tr_x_xzz_xxxz[i] = 3.0 * tr_x_zz_xxz[i] * fe_0 + ts_zz_xxxz[i] * fe_0 + tr_x_zz_xxxz[i] * pa_x[i];

        tr_x_xzz_xxyy[i] = tr_x_x_xxyy[i] * fe_0 + tr_x_xz_xxyy[i] * pa_z[i];

        tr_x_xzz_xxyz[i] = 2.0 * tr_x_zz_xyz[i] * fe_0 + ts_zz_xxyz[i] * fe_0 + tr_x_zz_xxyz[i] * pa_x[i];

        tr_x_xzz_xxzz[i] = 2.0 * tr_x_zz_xzz[i] * fe_0 + ts_zz_xxzz[i] * fe_0 + tr_x_zz_xxzz[i] * pa_x[i];

        tr_x_xzz_xyyy[i] = tr_x_x_xyyy[i] * fe_0 + tr_x_xz_xyyy[i] * pa_z[i];

        tr_x_xzz_xyyz[i] = tr_x_zz_yyz[i] * fe_0 + ts_zz_xyyz[i] * fe_0 + tr_x_zz_xyyz[i] * pa_x[i];

        tr_x_xzz_xyzz[i] = tr_x_zz_yzz[i] * fe_0 + ts_zz_xyzz[i] * fe_0 + tr_x_zz_xyzz[i] * pa_x[i];

        tr_x_xzz_xzzz[i] = tr_x_zz_zzz[i] * fe_0 + ts_zz_xzzz[i] * fe_0 + tr_x_zz_xzzz[i] * pa_x[i];

        tr_x_xzz_yyyy[i] = ts_zz_yyyy[i] * fe_0 + tr_x_zz_yyyy[i] * pa_x[i];

        tr_x_xzz_yyyz[i] = ts_zz_yyyz[i] * fe_0 + tr_x_zz_yyyz[i] * pa_x[i];

        tr_x_xzz_yyzz[i] = ts_zz_yyzz[i] * fe_0 + tr_x_zz_yyzz[i] * pa_x[i];

        tr_x_xzz_yzzz[i] = ts_zz_yzzz[i] * fe_0 + tr_x_zz_yzzz[i] * pa_x[i];

        tr_x_xzz_zzzz[i] = ts_zz_zzzz[i] * fe_0 + tr_x_zz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto tr_x_yyy_xxxx = pbuffer.data(idx_dip_fg + 90);

    auto tr_x_yyy_xxxy = pbuffer.data(idx_dip_fg + 91);

    auto tr_x_yyy_xxxz = pbuffer.data(idx_dip_fg + 92);

    auto tr_x_yyy_xxyy = pbuffer.data(idx_dip_fg + 93);

    auto tr_x_yyy_xxyz = pbuffer.data(idx_dip_fg + 94);

    auto tr_x_yyy_xxzz = pbuffer.data(idx_dip_fg + 95);

    auto tr_x_yyy_xyyy = pbuffer.data(idx_dip_fg + 96);

    auto tr_x_yyy_xyyz = pbuffer.data(idx_dip_fg + 97);

    auto tr_x_yyy_xyzz = pbuffer.data(idx_dip_fg + 98);

    auto tr_x_yyy_xzzz = pbuffer.data(idx_dip_fg + 99);

    auto tr_x_yyy_yyyy = pbuffer.data(idx_dip_fg + 100);

    auto tr_x_yyy_yyyz = pbuffer.data(idx_dip_fg + 101);

    auto tr_x_yyy_yyzz = pbuffer.data(idx_dip_fg + 102);

    auto tr_x_yyy_yzzz = pbuffer.data(idx_dip_fg + 103);

    auto tr_x_yyy_zzzz = pbuffer.data(idx_dip_fg + 104);

#pragma omp simd aligned(pa_y,              \
                             tr_x_y_xxxx,   \
                             tr_x_y_xxxy,   \
                             tr_x_y_xxxz,   \
                             tr_x_y_xxyy,   \
                             tr_x_y_xxyz,   \
                             tr_x_y_xxzz,   \
                             tr_x_y_xyyy,   \
                             tr_x_y_xyyz,   \
                             tr_x_y_xyzz,   \
                             tr_x_y_xzzz,   \
                             tr_x_y_yyyy,   \
                             tr_x_y_yyyz,   \
                             tr_x_y_yyzz,   \
                             tr_x_y_yzzz,   \
                             tr_x_y_zzzz,   \
                             tr_x_yy_xxx,   \
                             tr_x_yy_xxxx,  \
                             tr_x_yy_xxxy,  \
                             tr_x_yy_xxxz,  \
                             tr_x_yy_xxy,   \
                             tr_x_yy_xxyy,  \
                             tr_x_yy_xxyz,  \
                             tr_x_yy_xxz,   \
                             tr_x_yy_xxzz,  \
                             tr_x_yy_xyy,   \
                             tr_x_yy_xyyy,  \
                             tr_x_yy_xyyz,  \
                             tr_x_yy_xyz,   \
                             tr_x_yy_xyzz,  \
                             tr_x_yy_xzz,   \
                             tr_x_yy_xzzz,  \
                             tr_x_yy_yyy,   \
                             tr_x_yy_yyyy,  \
                             tr_x_yy_yyyz,  \
                             tr_x_yy_yyz,   \
                             tr_x_yy_yyzz,  \
                             tr_x_yy_yzz,   \
                             tr_x_yy_yzzz,  \
                             tr_x_yy_zzz,   \
                             tr_x_yy_zzzz,  \
                             tr_x_yyy_xxxx, \
                             tr_x_yyy_xxxy, \
                             tr_x_yyy_xxxz, \
                             tr_x_yyy_xxyy, \
                             tr_x_yyy_xxyz, \
                             tr_x_yyy_xxzz, \
                             tr_x_yyy_xyyy, \
                             tr_x_yyy_xyyz, \
                             tr_x_yyy_xyzz, \
                             tr_x_yyy_xzzz, \
                             tr_x_yyy_yyyy, \
                             tr_x_yyy_yyyz, \
                             tr_x_yyy_yyzz, \
                             tr_x_yyy_yzzz, \
                             tr_x_yyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_xxxx[i] = 2.0 * tr_x_y_xxxx[i] * fe_0 + tr_x_yy_xxxx[i] * pa_y[i];

        tr_x_yyy_xxxy[i] = 2.0 * tr_x_y_xxxy[i] * fe_0 + tr_x_yy_xxx[i] * fe_0 + tr_x_yy_xxxy[i] * pa_y[i];

        tr_x_yyy_xxxz[i] = 2.0 * tr_x_y_xxxz[i] * fe_0 + tr_x_yy_xxxz[i] * pa_y[i];

        tr_x_yyy_xxyy[i] = 2.0 * tr_x_y_xxyy[i] * fe_0 + 2.0 * tr_x_yy_xxy[i] * fe_0 + tr_x_yy_xxyy[i] * pa_y[i];

        tr_x_yyy_xxyz[i] = 2.0 * tr_x_y_xxyz[i] * fe_0 + tr_x_yy_xxz[i] * fe_0 + tr_x_yy_xxyz[i] * pa_y[i];

        tr_x_yyy_xxzz[i] = 2.0 * tr_x_y_xxzz[i] * fe_0 + tr_x_yy_xxzz[i] * pa_y[i];

        tr_x_yyy_xyyy[i] = 2.0 * tr_x_y_xyyy[i] * fe_0 + 3.0 * tr_x_yy_xyy[i] * fe_0 + tr_x_yy_xyyy[i] * pa_y[i];

        tr_x_yyy_xyyz[i] = 2.0 * tr_x_y_xyyz[i] * fe_0 + 2.0 * tr_x_yy_xyz[i] * fe_0 + tr_x_yy_xyyz[i] * pa_y[i];

        tr_x_yyy_xyzz[i] = 2.0 * tr_x_y_xyzz[i] * fe_0 + tr_x_yy_xzz[i] * fe_0 + tr_x_yy_xyzz[i] * pa_y[i];

        tr_x_yyy_xzzz[i] = 2.0 * tr_x_y_xzzz[i] * fe_0 + tr_x_yy_xzzz[i] * pa_y[i];

        tr_x_yyy_yyyy[i] = 2.0 * tr_x_y_yyyy[i] * fe_0 + 4.0 * tr_x_yy_yyy[i] * fe_0 + tr_x_yy_yyyy[i] * pa_y[i];

        tr_x_yyy_yyyz[i] = 2.0 * tr_x_y_yyyz[i] * fe_0 + 3.0 * tr_x_yy_yyz[i] * fe_0 + tr_x_yy_yyyz[i] * pa_y[i];

        tr_x_yyy_yyzz[i] = 2.0 * tr_x_y_yyzz[i] * fe_0 + 2.0 * tr_x_yy_yzz[i] * fe_0 + tr_x_yy_yyzz[i] * pa_y[i];

        tr_x_yyy_yzzz[i] = 2.0 * tr_x_y_yzzz[i] * fe_0 + tr_x_yy_zzz[i] * fe_0 + tr_x_yy_yzzz[i] * pa_y[i];

        tr_x_yyy_zzzz[i] = 2.0 * tr_x_y_zzzz[i] * fe_0 + tr_x_yy_zzzz[i] * pa_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto tr_x_yyz_xxxx = pbuffer.data(idx_dip_fg + 105);

    auto tr_x_yyz_xxxy = pbuffer.data(idx_dip_fg + 106);

    auto tr_x_yyz_xxxz = pbuffer.data(idx_dip_fg + 107);

    auto tr_x_yyz_xxyy = pbuffer.data(idx_dip_fg + 108);

    auto tr_x_yyz_xxyz = pbuffer.data(idx_dip_fg + 109);

    auto tr_x_yyz_xxzz = pbuffer.data(idx_dip_fg + 110);

    auto tr_x_yyz_xyyy = pbuffer.data(idx_dip_fg + 111);

    auto tr_x_yyz_xyyz = pbuffer.data(idx_dip_fg + 112);

    auto tr_x_yyz_xyzz = pbuffer.data(idx_dip_fg + 113);

    auto tr_x_yyz_xzzz = pbuffer.data(idx_dip_fg + 114);

    auto tr_x_yyz_yyyy = pbuffer.data(idx_dip_fg + 115);

    auto tr_x_yyz_yyyz = pbuffer.data(idx_dip_fg + 116);

    auto tr_x_yyz_yyzz = pbuffer.data(idx_dip_fg + 117);

    auto tr_x_yyz_yzzz = pbuffer.data(idx_dip_fg + 118);

    auto tr_x_yyz_zzzz = pbuffer.data(idx_dip_fg + 119);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yy_xxxx,  \
                             tr_x_yy_xxxy,  \
                             tr_x_yy_xxy,   \
                             tr_x_yy_xxyy,  \
                             tr_x_yy_xxyz,  \
                             tr_x_yy_xyy,   \
                             tr_x_yy_xyyy,  \
                             tr_x_yy_xyyz,  \
                             tr_x_yy_xyz,   \
                             tr_x_yy_xyzz,  \
                             tr_x_yy_yyy,   \
                             tr_x_yy_yyyy,  \
                             tr_x_yy_yyyz,  \
                             tr_x_yy_yyz,   \
                             tr_x_yy_yyzz,  \
                             tr_x_yy_yzz,   \
                             tr_x_yy_yzzz,  \
                             tr_x_yyz_xxxx, \
                             tr_x_yyz_xxxy, \
                             tr_x_yyz_xxxz, \
                             tr_x_yyz_xxyy, \
                             tr_x_yyz_xxyz, \
                             tr_x_yyz_xxzz, \
                             tr_x_yyz_xyyy, \
                             tr_x_yyz_xyyz, \
                             tr_x_yyz_xyzz, \
                             tr_x_yyz_xzzz, \
                             tr_x_yyz_yyyy, \
                             tr_x_yyz_yyyz, \
                             tr_x_yyz_yyzz, \
                             tr_x_yyz_yzzz, \
                             tr_x_yyz_zzzz, \
                             tr_x_yz_xxxz,  \
                             tr_x_yz_xxzz,  \
                             tr_x_yz_xzzz,  \
                             tr_x_yz_zzzz,  \
                             tr_x_z_xxxz,   \
                             tr_x_z_xxzz,   \
                             tr_x_z_xzzz,   \
                             tr_x_z_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_xxxx[i] = tr_x_yy_xxxx[i] * pa_z[i];

        tr_x_yyz_xxxy[i] = tr_x_yy_xxxy[i] * pa_z[i];

        tr_x_yyz_xxxz[i] = tr_x_z_xxxz[i] * fe_0 + tr_x_yz_xxxz[i] * pa_y[i];

        tr_x_yyz_xxyy[i] = tr_x_yy_xxyy[i] * pa_z[i];

        tr_x_yyz_xxyz[i] = tr_x_yy_xxy[i] * fe_0 + tr_x_yy_xxyz[i] * pa_z[i];

        tr_x_yyz_xxzz[i] = tr_x_z_xxzz[i] * fe_0 + tr_x_yz_xxzz[i] * pa_y[i];

        tr_x_yyz_xyyy[i] = tr_x_yy_xyyy[i] * pa_z[i];

        tr_x_yyz_xyyz[i] = tr_x_yy_xyy[i] * fe_0 + tr_x_yy_xyyz[i] * pa_z[i];

        tr_x_yyz_xyzz[i] = 2.0 * tr_x_yy_xyz[i] * fe_0 + tr_x_yy_xyzz[i] * pa_z[i];

        tr_x_yyz_xzzz[i] = tr_x_z_xzzz[i] * fe_0 + tr_x_yz_xzzz[i] * pa_y[i];

        tr_x_yyz_yyyy[i] = tr_x_yy_yyyy[i] * pa_z[i];

        tr_x_yyz_yyyz[i] = tr_x_yy_yyy[i] * fe_0 + tr_x_yy_yyyz[i] * pa_z[i];

        tr_x_yyz_yyzz[i] = 2.0 * tr_x_yy_yyz[i] * fe_0 + tr_x_yy_yyzz[i] * pa_z[i];

        tr_x_yyz_yzzz[i] = 3.0 * tr_x_yy_yzz[i] * fe_0 + tr_x_yy_yzzz[i] * pa_z[i];

        tr_x_yyz_zzzz[i] = tr_x_z_zzzz[i] * fe_0 + tr_x_yz_zzzz[i] * pa_y[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto tr_x_yzz_xxxx = pbuffer.data(idx_dip_fg + 120);

    auto tr_x_yzz_xxxy = pbuffer.data(idx_dip_fg + 121);

    auto tr_x_yzz_xxxz = pbuffer.data(idx_dip_fg + 122);

    auto tr_x_yzz_xxyy = pbuffer.data(idx_dip_fg + 123);

    auto tr_x_yzz_xxyz = pbuffer.data(idx_dip_fg + 124);

    auto tr_x_yzz_xxzz = pbuffer.data(idx_dip_fg + 125);

    auto tr_x_yzz_xyyy = pbuffer.data(idx_dip_fg + 126);

    auto tr_x_yzz_xyyz = pbuffer.data(idx_dip_fg + 127);

    auto tr_x_yzz_xyzz = pbuffer.data(idx_dip_fg + 128);

    auto tr_x_yzz_xzzz = pbuffer.data(idx_dip_fg + 129);

    auto tr_x_yzz_yyyy = pbuffer.data(idx_dip_fg + 130);

    auto tr_x_yzz_yyyz = pbuffer.data(idx_dip_fg + 131);

    auto tr_x_yzz_yyzz = pbuffer.data(idx_dip_fg + 132);

    auto tr_x_yzz_yzzz = pbuffer.data(idx_dip_fg + 133);

    auto tr_x_yzz_zzzz = pbuffer.data(idx_dip_fg + 134);

#pragma omp simd aligned(pa_y,              \
                             tr_x_yzz_xxxx, \
                             tr_x_yzz_xxxy, \
                             tr_x_yzz_xxxz, \
                             tr_x_yzz_xxyy, \
                             tr_x_yzz_xxyz, \
                             tr_x_yzz_xxzz, \
                             tr_x_yzz_xyyy, \
                             tr_x_yzz_xyyz, \
                             tr_x_yzz_xyzz, \
                             tr_x_yzz_xzzz, \
                             tr_x_yzz_yyyy, \
                             tr_x_yzz_yyyz, \
                             tr_x_yzz_yyzz, \
                             tr_x_yzz_yzzz, \
                             tr_x_yzz_zzzz, \
                             tr_x_zz_xxx,   \
                             tr_x_zz_xxxx,  \
                             tr_x_zz_xxxy,  \
                             tr_x_zz_xxxz,  \
                             tr_x_zz_xxy,   \
                             tr_x_zz_xxyy,  \
                             tr_x_zz_xxyz,  \
                             tr_x_zz_xxz,   \
                             tr_x_zz_xxzz,  \
                             tr_x_zz_xyy,   \
                             tr_x_zz_xyyy,  \
                             tr_x_zz_xyyz,  \
                             tr_x_zz_xyz,   \
                             tr_x_zz_xyzz,  \
                             tr_x_zz_xzz,   \
                             tr_x_zz_xzzz,  \
                             tr_x_zz_yyy,   \
                             tr_x_zz_yyyy,  \
                             tr_x_zz_yyyz,  \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yyzz,  \
                             tr_x_zz_yzz,   \
                             tr_x_zz_yzzz,  \
                             tr_x_zz_zzz,   \
                             tr_x_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_xxxx[i] = tr_x_zz_xxxx[i] * pa_y[i];

        tr_x_yzz_xxxy[i] = tr_x_zz_xxx[i] * fe_0 + tr_x_zz_xxxy[i] * pa_y[i];

        tr_x_yzz_xxxz[i] = tr_x_zz_xxxz[i] * pa_y[i];

        tr_x_yzz_xxyy[i] = 2.0 * tr_x_zz_xxy[i] * fe_0 + tr_x_zz_xxyy[i] * pa_y[i];

        tr_x_yzz_xxyz[i] = tr_x_zz_xxz[i] * fe_0 + tr_x_zz_xxyz[i] * pa_y[i];

        tr_x_yzz_xxzz[i] = tr_x_zz_xxzz[i] * pa_y[i];

        tr_x_yzz_xyyy[i] = 3.0 * tr_x_zz_xyy[i] * fe_0 + tr_x_zz_xyyy[i] * pa_y[i];

        tr_x_yzz_xyyz[i] = 2.0 * tr_x_zz_xyz[i] * fe_0 + tr_x_zz_xyyz[i] * pa_y[i];

        tr_x_yzz_xyzz[i] = tr_x_zz_xzz[i] * fe_0 + tr_x_zz_xyzz[i] * pa_y[i];

        tr_x_yzz_xzzz[i] = tr_x_zz_xzzz[i] * pa_y[i];

        tr_x_yzz_yyyy[i] = 4.0 * tr_x_zz_yyy[i] * fe_0 + tr_x_zz_yyyy[i] * pa_y[i];

        tr_x_yzz_yyyz[i] = 3.0 * tr_x_zz_yyz[i] * fe_0 + tr_x_zz_yyyz[i] * pa_y[i];

        tr_x_yzz_yyzz[i] = 2.0 * tr_x_zz_yzz[i] * fe_0 + tr_x_zz_yyzz[i] * pa_y[i];

        tr_x_yzz_yzzz[i] = tr_x_zz_zzz[i] * fe_0 + tr_x_zz_yzzz[i] * pa_y[i];

        tr_x_yzz_zzzz[i] = tr_x_zz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto tr_x_zzz_xxxx = pbuffer.data(idx_dip_fg + 135);

    auto tr_x_zzz_xxxy = pbuffer.data(idx_dip_fg + 136);

    auto tr_x_zzz_xxxz = pbuffer.data(idx_dip_fg + 137);

    auto tr_x_zzz_xxyy = pbuffer.data(idx_dip_fg + 138);

    auto tr_x_zzz_xxyz = pbuffer.data(idx_dip_fg + 139);

    auto tr_x_zzz_xxzz = pbuffer.data(idx_dip_fg + 140);

    auto tr_x_zzz_xyyy = pbuffer.data(idx_dip_fg + 141);

    auto tr_x_zzz_xyyz = pbuffer.data(idx_dip_fg + 142);

    auto tr_x_zzz_xyzz = pbuffer.data(idx_dip_fg + 143);

    auto tr_x_zzz_xzzz = pbuffer.data(idx_dip_fg + 144);

    auto tr_x_zzz_yyyy = pbuffer.data(idx_dip_fg + 145);

    auto tr_x_zzz_yyyz = pbuffer.data(idx_dip_fg + 146);

    auto tr_x_zzz_yyzz = pbuffer.data(idx_dip_fg + 147);

    auto tr_x_zzz_yzzz = pbuffer.data(idx_dip_fg + 148);

    auto tr_x_zzz_zzzz = pbuffer.data(idx_dip_fg + 149);

#pragma omp simd aligned(pa_z,              \
                             tr_x_z_xxxx,   \
                             tr_x_z_xxxy,   \
                             tr_x_z_xxxz,   \
                             tr_x_z_xxyy,   \
                             tr_x_z_xxyz,   \
                             tr_x_z_xxzz,   \
                             tr_x_z_xyyy,   \
                             tr_x_z_xyyz,   \
                             tr_x_z_xyzz,   \
                             tr_x_z_xzzz,   \
                             tr_x_z_yyyy,   \
                             tr_x_z_yyyz,   \
                             tr_x_z_yyzz,   \
                             tr_x_z_yzzz,   \
                             tr_x_z_zzzz,   \
                             tr_x_zz_xxx,   \
                             tr_x_zz_xxxx,  \
                             tr_x_zz_xxxy,  \
                             tr_x_zz_xxxz,  \
                             tr_x_zz_xxy,   \
                             tr_x_zz_xxyy,  \
                             tr_x_zz_xxyz,  \
                             tr_x_zz_xxz,   \
                             tr_x_zz_xxzz,  \
                             tr_x_zz_xyy,   \
                             tr_x_zz_xyyy,  \
                             tr_x_zz_xyyz,  \
                             tr_x_zz_xyz,   \
                             tr_x_zz_xyzz,  \
                             tr_x_zz_xzz,   \
                             tr_x_zz_xzzz,  \
                             tr_x_zz_yyy,   \
                             tr_x_zz_yyyy,  \
                             tr_x_zz_yyyz,  \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yyzz,  \
                             tr_x_zz_yzz,   \
                             tr_x_zz_yzzz,  \
                             tr_x_zz_zzz,   \
                             tr_x_zz_zzzz,  \
                             tr_x_zzz_xxxx, \
                             tr_x_zzz_xxxy, \
                             tr_x_zzz_xxxz, \
                             tr_x_zzz_xxyy, \
                             tr_x_zzz_xxyz, \
                             tr_x_zzz_xxzz, \
                             tr_x_zzz_xyyy, \
                             tr_x_zzz_xyyz, \
                             tr_x_zzz_xyzz, \
                             tr_x_zzz_xzzz, \
                             tr_x_zzz_yyyy, \
                             tr_x_zzz_yyyz, \
                             tr_x_zzz_yyzz, \
                             tr_x_zzz_yzzz, \
                             tr_x_zzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_xxxx[i] = 2.0 * tr_x_z_xxxx[i] * fe_0 + tr_x_zz_xxxx[i] * pa_z[i];

        tr_x_zzz_xxxy[i] = 2.0 * tr_x_z_xxxy[i] * fe_0 + tr_x_zz_xxxy[i] * pa_z[i];

        tr_x_zzz_xxxz[i] = 2.0 * tr_x_z_xxxz[i] * fe_0 + tr_x_zz_xxx[i] * fe_0 + tr_x_zz_xxxz[i] * pa_z[i];

        tr_x_zzz_xxyy[i] = 2.0 * tr_x_z_xxyy[i] * fe_0 + tr_x_zz_xxyy[i] * pa_z[i];

        tr_x_zzz_xxyz[i] = 2.0 * tr_x_z_xxyz[i] * fe_0 + tr_x_zz_xxy[i] * fe_0 + tr_x_zz_xxyz[i] * pa_z[i];

        tr_x_zzz_xxzz[i] = 2.0 * tr_x_z_xxzz[i] * fe_0 + 2.0 * tr_x_zz_xxz[i] * fe_0 + tr_x_zz_xxzz[i] * pa_z[i];

        tr_x_zzz_xyyy[i] = 2.0 * tr_x_z_xyyy[i] * fe_0 + tr_x_zz_xyyy[i] * pa_z[i];

        tr_x_zzz_xyyz[i] = 2.0 * tr_x_z_xyyz[i] * fe_0 + tr_x_zz_xyy[i] * fe_0 + tr_x_zz_xyyz[i] * pa_z[i];

        tr_x_zzz_xyzz[i] = 2.0 * tr_x_z_xyzz[i] * fe_0 + 2.0 * tr_x_zz_xyz[i] * fe_0 + tr_x_zz_xyzz[i] * pa_z[i];

        tr_x_zzz_xzzz[i] = 2.0 * tr_x_z_xzzz[i] * fe_0 + 3.0 * tr_x_zz_xzz[i] * fe_0 + tr_x_zz_xzzz[i] * pa_z[i];

        tr_x_zzz_yyyy[i] = 2.0 * tr_x_z_yyyy[i] * fe_0 + tr_x_zz_yyyy[i] * pa_z[i];

        tr_x_zzz_yyyz[i] = 2.0 * tr_x_z_yyyz[i] * fe_0 + tr_x_zz_yyy[i] * fe_0 + tr_x_zz_yyyz[i] * pa_z[i];

        tr_x_zzz_yyzz[i] = 2.0 * tr_x_z_yyzz[i] * fe_0 + 2.0 * tr_x_zz_yyz[i] * fe_0 + tr_x_zz_yyzz[i] * pa_z[i];

        tr_x_zzz_yzzz[i] = 2.0 * tr_x_z_yzzz[i] * fe_0 + 3.0 * tr_x_zz_yzz[i] * fe_0 + tr_x_zz_yzzz[i] * pa_z[i];

        tr_x_zzz_zzzz[i] = 2.0 * tr_x_z_zzzz[i] * fe_0 + 4.0 * tr_x_zz_zzz[i] * fe_0 + tr_x_zz_zzzz[i] * pa_z[i];
    }

    // Set up 150-165 components of targeted buffer : FG

    auto tr_y_xxx_xxxx = pbuffer.data(idx_dip_fg + 150);

    auto tr_y_xxx_xxxy = pbuffer.data(idx_dip_fg + 151);

    auto tr_y_xxx_xxxz = pbuffer.data(idx_dip_fg + 152);

    auto tr_y_xxx_xxyy = pbuffer.data(idx_dip_fg + 153);

    auto tr_y_xxx_xxyz = pbuffer.data(idx_dip_fg + 154);

    auto tr_y_xxx_xxzz = pbuffer.data(idx_dip_fg + 155);

    auto tr_y_xxx_xyyy = pbuffer.data(idx_dip_fg + 156);

    auto tr_y_xxx_xyyz = pbuffer.data(idx_dip_fg + 157);

    auto tr_y_xxx_xyzz = pbuffer.data(idx_dip_fg + 158);

    auto tr_y_xxx_xzzz = pbuffer.data(idx_dip_fg + 159);

    auto tr_y_xxx_yyyy = pbuffer.data(idx_dip_fg + 160);

    auto tr_y_xxx_yyyz = pbuffer.data(idx_dip_fg + 161);

    auto tr_y_xxx_yyzz = pbuffer.data(idx_dip_fg + 162);

    auto tr_y_xxx_yzzz = pbuffer.data(idx_dip_fg + 163);

    auto tr_y_xxx_zzzz = pbuffer.data(idx_dip_fg + 164);

#pragma omp simd aligned(pa_x,              \
                             tr_y_x_xxxx,   \
                             tr_y_x_xxxy,   \
                             tr_y_x_xxxz,   \
                             tr_y_x_xxyy,   \
                             tr_y_x_xxyz,   \
                             tr_y_x_xxzz,   \
                             tr_y_x_xyyy,   \
                             tr_y_x_xyyz,   \
                             tr_y_x_xyzz,   \
                             tr_y_x_xzzz,   \
                             tr_y_x_yyyy,   \
                             tr_y_x_yyyz,   \
                             tr_y_x_yyzz,   \
                             tr_y_x_yzzz,   \
                             tr_y_x_zzzz,   \
                             tr_y_xx_xxx,   \
                             tr_y_xx_xxxx,  \
                             tr_y_xx_xxxy,  \
                             tr_y_xx_xxxz,  \
                             tr_y_xx_xxy,   \
                             tr_y_xx_xxyy,  \
                             tr_y_xx_xxyz,  \
                             tr_y_xx_xxz,   \
                             tr_y_xx_xxzz,  \
                             tr_y_xx_xyy,   \
                             tr_y_xx_xyyy,  \
                             tr_y_xx_xyyz,  \
                             tr_y_xx_xyz,   \
                             tr_y_xx_xyzz,  \
                             tr_y_xx_xzz,   \
                             tr_y_xx_xzzz,  \
                             tr_y_xx_yyy,   \
                             tr_y_xx_yyyy,  \
                             tr_y_xx_yyyz,  \
                             tr_y_xx_yyz,   \
                             tr_y_xx_yyzz,  \
                             tr_y_xx_yzz,   \
                             tr_y_xx_yzzz,  \
                             tr_y_xx_zzz,   \
                             tr_y_xx_zzzz,  \
                             tr_y_xxx_xxxx, \
                             tr_y_xxx_xxxy, \
                             tr_y_xxx_xxxz, \
                             tr_y_xxx_xxyy, \
                             tr_y_xxx_xxyz, \
                             tr_y_xxx_xxzz, \
                             tr_y_xxx_xyyy, \
                             tr_y_xxx_xyyz, \
                             tr_y_xxx_xyzz, \
                             tr_y_xxx_xzzz, \
                             tr_y_xxx_yyyy, \
                             tr_y_xxx_yyyz, \
                             tr_y_xxx_yyzz, \
                             tr_y_xxx_yzzz, \
                             tr_y_xxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_xxxx[i] = 2.0 * tr_y_x_xxxx[i] * fe_0 + 4.0 * tr_y_xx_xxx[i] * fe_0 + tr_y_xx_xxxx[i] * pa_x[i];

        tr_y_xxx_xxxy[i] = 2.0 * tr_y_x_xxxy[i] * fe_0 + 3.0 * tr_y_xx_xxy[i] * fe_0 + tr_y_xx_xxxy[i] * pa_x[i];

        tr_y_xxx_xxxz[i] = 2.0 * tr_y_x_xxxz[i] * fe_0 + 3.0 * tr_y_xx_xxz[i] * fe_0 + tr_y_xx_xxxz[i] * pa_x[i];

        tr_y_xxx_xxyy[i] = 2.0 * tr_y_x_xxyy[i] * fe_0 + 2.0 * tr_y_xx_xyy[i] * fe_0 + tr_y_xx_xxyy[i] * pa_x[i];

        tr_y_xxx_xxyz[i] = 2.0 * tr_y_x_xxyz[i] * fe_0 + 2.0 * tr_y_xx_xyz[i] * fe_0 + tr_y_xx_xxyz[i] * pa_x[i];

        tr_y_xxx_xxzz[i] = 2.0 * tr_y_x_xxzz[i] * fe_0 + 2.0 * tr_y_xx_xzz[i] * fe_0 + tr_y_xx_xxzz[i] * pa_x[i];

        tr_y_xxx_xyyy[i] = 2.0 * tr_y_x_xyyy[i] * fe_0 + tr_y_xx_yyy[i] * fe_0 + tr_y_xx_xyyy[i] * pa_x[i];

        tr_y_xxx_xyyz[i] = 2.0 * tr_y_x_xyyz[i] * fe_0 + tr_y_xx_yyz[i] * fe_0 + tr_y_xx_xyyz[i] * pa_x[i];

        tr_y_xxx_xyzz[i] = 2.0 * tr_y_x_xyzz[i] * fe_0 + tr_y_xx_yzz[i] * fe_0 + tr_y_xx_xyzz[i] * pa_x[i];

        tr_y_xxx_xzzz[i] = 2.0 * tr_y_x_xzzz[i] * fe_0 + tr_y_xx_zzz[i] * fe_0 + tr_y_xx_xzzz[i] * pa_x[i];

        tr_y_xxx_yyyy[i] = 2.0 * tr_y_x_yyyy[i] * fe_0 + tr_y_xx_yyyy[i] * pa_x[i];

        tr_y_xxx_yyyz[i] = 2.0 * tr_y_x_yyyz[i] * fe_0 + tr_y_xx_yyyz[i] * pa_x[i];

        tr_y_xxx_yyzz[i] = 2.0 * tr_y_x_yyzz[i] * fe_0 + tr_y_xx_yyzz[i] * pa_x[i];

        tr_y_xxx_yzzz[i] = 2.0 * tr_y_x_yzzz[i] * fe_0 + tr_y_xx_yzzz[i] * pa_x[i];

        tr_y_xxx_zzzz[i] = 2.0 * tr_y_x_zzzz[i] * fe_0 + tr_y_xx_zzzz[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : FG

    auto tr_y_xxy_xxxx = pbuffer.data(idx_dip_fg + 165);

    auto tr_y_xxy_xxxy = pbuffer.data(idx_dip_fg + 166);

    auto tr_y_xxy_xxxz = pbuffer.data(idx_dip_fg + 167);

    auto tr_y_xxy_xxyy = pbuffer.data(idx_dip_fg + 168);

    auto tr_y_xxy_xxyz = pbuffer.data(idx_dip_fg + 169);

    auto tr_y_xxy_xxzz = pbuffer.data(idx_dip_fg + 170);

    auto tr_y_xxy_xyyy = pbuffer.data(idx_dip_fg + 171);

    auto tr_y_xxy_xyyz = pbuffer.data(idx_dip_fg + 172);

    auto tr_y_xxy_xyzz = pbuffer.data(idx_dip_fg + 173);

    auto tr_y_xxy_xzzz = pbuffer.data(idx_dip_fg + 174);

    auto tr_y_xxy_yyyy = pbuffer.data(idx_dip_fg + 175);

    auto tr_y_xxy_yyyz = pbuffer.data(idx_dip_fg + 176);

    auto tr_y_xxy_yyzz = pbuffer.data(idx_dip_fg + 177);

    auto tr_y_xxy_yzzz = pbuffer.data(idx_dip_fg + 178);

    auto tr_y_xxy_zzzz = pbuffer.data(idx_dip_fg + 179);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_y_xx_xxxx,  \
                             tr_y_xx_xxxz,  \
                             tr_y_xx_xxzz,  \
                             tr_y_xx_xzzz,  \
                             tr_y_xxy_xxxx, \
                             tr_y_xxy_xxxy, \
                             tr_y_xxy_xxxz, \
                             tr_y_xxy_xxyy, \
                             tr_y_xxy_xxyz, \
                             tr_y_xxy_xxzz, \
                             tr_y_xxy_xyyy, \
                             tr_y_xxy_xyyz, \
                             tr_y_xxy_xyzz, \
                             tr_y_xxy_xzzz, \
                             tr_y_xxy_yyyy, \
                             tr_y_xxy_yyyz, \
                             tr_y_xxy_yyzz, \
                             tr_y_xxy_yzzz, \
                             tr_y_xxy_zzzz, \
                             tr_y_xy_xxxy,  \
                             tr_y_xy_xxy,   \
                             tr_y_xy_xxyy,  \
                             tr_y_xy_xxyz,  \
                             tr_y_xy_xyy,   \
                             tr_y_xy_xyyy,  \
                             tr_y_xy_xyyz,  \
                             tr_y_xy_xyz,   \
                             tr_y_xy_xyzz,  \
                             tr_y_xy_yyy,   \
                             tr_y_xy_yyyy,  \
                             tr_y_xy_yyyz,  \
                             tr_y_xy_yyz,   \
                             tr_y_xy_yyzz,  \
                             tr_y_xy_yzz,   \
                             tr_y_xy_yzzz,  \
                             tr_y_xy_zzzz,  \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxyy,   \
                             tr_y_y_xxyz,   \
                             tr_y_y_xyyy,   \
                             tr_y_y_xyyz,   \
                             tr_y_y_xyzz,   \
                             tr_y_y_yyyy,   \
                             tr_y_y_yyyz,   \
                             tr_y_y_yyzz,   \
                             tr_y_y_yzzz,   \
                             tr_y_y_zzzz,   \
                             ts_xx_xxxx,    \
                             ts_xx_xxxz,    \
                             ts_xx_xxzz,    \
                             ts_xx_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_xxxx[i] = ts_xx_xxxx[i] * fe_0 + tr_y_xx_xxxx[i] * pa_y[i];

        tr_y_xxy_xxxy[i] = tr_y_y_xxxy[i] * fe_0 + 3.0 * tr_y_xy_xxy[i] * fe_0 + tr_y_xy_xxxy[i] * pa_x[i];

        tr_y_xxy_xxxz[i] = ts_xx_xxxz[i] * fe_0 + tr_y_xx_xxxz[i] * pa_y[i];

        tr_y_xxy_xxyy[i] = tr_y_y_xxyy[i] * fe_0 + 2.0 * tr_y_xy_xyy[i] * fe_0 + tr_y_xy_xxyy[i] * pa_x[i];

        tr_y_xxy_xxyz[i] = tr_y_y_xxyz[i] * fe_0 + 2.0 * tr_y_xy_xyz[i] * fe_0 + tr_y_xy_xxyz[i] * pa_x[i];

        tr_y_xxy_xxzz[i] = ts_xx_xxzz[i] * fe_0 + tr_y_xx_xxzz[i] * pa_y[i];

        tr_y_xxy_xyyy[i] = tr_y_y_xyyy[i] * fe_0 + tr_y_xy_yyy[i] * fe_0 + tr_y_xy_xyyy[i] * pa_x[i];

        tr_y_xxy_xyyz[i] = tr_y_y_xyyz[i] * fe_0 + tr_y_xy_yyz[i] * fe_0 + tr_y_xy_xyyz[i] * pa_x[i];

        tr_y_xxy_xyzz[i] = tr_y_y_xyzz[i] * fe_0 + tr_y_xy_yzz[i] * fe_0 + tr_y_xy_xyzz[i] * pa_x[i];

        tr_y_xxy_xzzz[i] = ts_xx_xzzz[i] * fe_0 + tr_y_xx_xzzz[i] * pa_y[i];

        tr_y_xxy_yyyy[i] = tr_y_y_yyyy[i] * fe_0 + tr_y_xy_yyyy[i] * pa_x[i];

        tr_y_xxy_yyyz[i] = tr_y_y_yyyz[i] * fe_0 + tr_y_xy_yyyz[i] * pa_x[i];

        tr_y_xxy_yyzz[i] = tr_y_y_yyzz[i] * fe_0 + tr_y_xy_yyzz[i] * pa_x[i];

        tr_y_xxy_yzzz[i] = tr_y_y_yzzz[i] * fe_0 + tr_y_xy_yzzz[i] * pa_x[i];

        tr_y_xxy_zzzz[i] = tr_y_y_zzzz[i] * fe_0 + tr_y_xy_zzzz[i] * pa_x[i];
    }

    // Set up 180-195 components of targeted buffer : FG

    auto tr_y_xxz_xxxx = pbuffer.data(idx_dip_fg + 180);

    auto tr_y_xxz_xxxy = pbuffer.data(idx_dip_fg + 181);

    auto tr_y_xxz_xxxz = pbuffer.data(idx_dip_fg + 182);

    auto tr_y_xxz_xxyy = pbuffer.data(idx_dip_fg + 183);

    auto tr_y_xxz_xxyz = pbuffer.data(idx_dip_fg + 184);

    auto tr_y_xxz_xxzz = pbuffer.data(idx_dip_fg + 185);

    auto tr_y_xxz_xyyy = pbuffer.data(idx_dip_fg + 186);

    auto tr_y_xxz_xyyz = pbuffer.data(idx_dip_fg + 187);

    auto tr_y_xxz_xyzz = pbuffer.data(idx_dip_fg + 188);

    auto tr_y_xxz_xzzz = pbuffer.data(idx_dip_fg + 189);

    auto tr_y_xxz_yyyy = pbuffer.data(idx_dip_fg + 190);

    auto tr_y_xxz_yyyz = pbuffer.data(idx_dip_fg + 191);

    auto tr_y_xxz_yyzz = pbuffer.data(idx_dip_fg + 192);

    auto tr_y_xxz_yzzz = pbuffer.data(idx_dip_fg + 193);

    auto tr_y_xxz_zzzz = pbuffer.data(idx_dip_fg + 194);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xx_xxx,   \
                             tr_y_xx_xxxx,  \
                             tr_y_xx_xxxy,  \
                             tr_y_xx_xxxz,  \
                             tr_y_xx_xxy,   \
                             tr_y_xx_xxyy,  \
                             tr_y_xx_xxyz,  \
                             tr_y_xx_xxz,   \
                             tr_y_xx_xxzz,  \
                             tr_y_xx_xyy,   \
                             tr_y_xx_xyyy,  \
                             tr_y_xx_xyyz,  \
                             tr_y_xx_xyz,   \
                             tr_y_xx_xyzz,  \
                             tr_y_xx_xzz,   \
                             tr_y_xx_xzzz,  \
                             tr_y_xx_yyyy,  \
                             tr_y_xxz_xxxx, \
                             tr_y_xxz_xxxy, \
                             tr_y_xxz_xxxz, \
                             tr_y_xxz_xxyy, \
                             tr_y_xxz_xxyz, \
                             tr_y_xxz_xxzz, \
                             tr_y_xxz_xyyy, \
                             tr_y_xxz_xyyz, \
                             tr_y_xxz_xyzz, \
                             tr_y_xxz_xzzz, \
                             tr_y_xxz_yyyy, \
                             tr_y_xxz_yyyz, \
                             tr_y_xxz_yyzz, \
                             tr_y_xxz_yzzz, \
                             tr_y_xxz_zzzz, \
                             tr_y_xz_yyyz,  \
                             tr_y_xz_yyzz,  \
                             tr_y_xz_yzzz,  \
                             tr_y_xz_zzzz,  \
                             tr_y_z_yyyz,   \
                             tr_y_z_yyzz,   \
                             tr_y_z_yzzz,   \
                             tr_y_z_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_xxxx[i] = tr_y_xx_xxxx[i] * pa_z[i];

        tr_y_xxz_xxxy[i] = tr_y_xx_xxxy[i] * pa_z[i];

        tr_y_xxz_xxxz[i] = tr_y_xx_xxx[i] * fe_0 + tr_y_xx_xxxz[i] * pa_z[i];

        tr_y_xxz_xxyy[i] = tr_y_xx_xxyy[i] * pa_z[i];

        tr_y_xxz_xxyz[i] = tr_y_xx_xxy[i] * fe_0 + tr_y_xx_xxyz[i] * pa_z[i];

        tr_y_xxz_xxzz[i] = 2.0 * tr_y_xx_xxz[i] * fe_0 + tr_y_xx_xxzz[i] * pa_z[i];

        tr_y_xxz_xyyy[i] = tr_y_xx_xyyy[i] * pa_z[i];

        tr_y_xxz_xyyz[i] = tr_y_xx_xyy[i] * fe_0 + tr_y_xx_xyyz[i] * pa_z[i];

        tr_y_xxz_xyzz[i] = 2.0 * tr_y_xx_xyz[i] * fe_0 + tr_y_xx_xyzz[i] * pa_z[i];

        tr_y_xxz_xzzz[i] = 3.0 * tr_y_xx_xzz[i] * fe_0 + tr_y_xx_xzzz[i] * pa_z[i];

        tr_y_xxz_yyyy[i] = tr_y_xx_yyyy[i] * pa_z[i];

        tr_y_xxz_yyyz[i] = tr_y_z_yyyz[i] * fe_0 + tr_y_xz_yyyz[i] * pa_x[i];

        tr_y_xxz_yyzz[i] = tr_y_z_yyzz[i] * fe_0 + tr_y_xz_yyzz[i] * pa_x[i];

        tr_y_xxz_yzzz[i] = tr_y_z_yzzz[i] * fe_0 + tr_y_xz_yzzz[i] * pa_x[i];

        tr_y_xxz_zzzz[i] = tr_y_z_zzzz[i] * fe_0 + tr_y_xz_zzzz[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : FG

    auto tr_y_xyy_xxxx = pbuffer.data(idx_dip_fg + 195);

    auto tr_y_xyy_xxxy = pbuffer.data(idx_dip_fg + 196);

    auto tr_y_xyy_xxxz = pbuffer.data(idx_dip_fg + 197);

    auto tr_y_xyy_xxyy = pbuffer.data(idx_dip_fg + 198);

    auto tr_y_xyy_xxyz = pbuffer.data(idx_dip_fg + 199);

    auto tr_y_xyy_xxzz = pbuffer.data(idx_dip_fg + 200);

    auto tr_y_xyy_xyyy = pbuffer.data(idx_dip_fg + 201);

    auto tr_y_xyy_xyyz = pbuffer.data(idx_dip_fg + 202);

    auto tr_y_xyy_xyzz = pbuffer.data(idx_dip_fg + 203);

    auto tr_y_xyy_xzzz = pbuffer.data(idx_dip_fg + 204);

    auto tr_y_xyy_yyyy = pbuffer.data(idx_dip_fg + 205);

    auto tr_y_xyy_yyyz = pbuffer.data(idx_dip_fg + 206);

    auto tr_y_xyy_yyzz = pbuffer.data(idx_dip_fg + 207);

    auto tr_y_xyy_yzzz = pbuffer.data(idx_dip_fg + 208);

    auto tr_y_xyy_zzzz = pbuffer.data(idx_dip_fg + 209);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xyy_xxxx, \
                             tr_y_xyy_xxxy, \
                             tr_y_xyy_xxxz, \
                             tr_y_xyy_xxyy, \
                             tr_y_xyy_xxyz, \
                             tr_y_xyy_xxzz, \
                             tr_y_xyy_xyyy, \
                             tr_y_xyy_xyyz, \
                             tr_y_xyy_xyzz, \
                             tr_y_xyy_xzzz, \
                             tr_y_xyy_yyyy, \
                             tr_y_xyy_yyyz, \
                             tr_y_xyy_yyzz, \
                             tr_y_xyy_yzzz, \
                             tr_y_xyy_zzzz, \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxxx,  \
                             tr_y_yy_xxxy,  \
                             tr_y_yy_xxxz,  \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xxyy,  \
                             tr_y_yy_xxyz,  \
                             tr_y_yy_xxz,   \
                             tr_y_yy_xxzz,  \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyyy,  \
                             tr_y_yy_xyyz,  \
                             tr_y_yy_xyz,   \
                             tr_y_yy_xyzz,  \
                             tr_y_yy_xzz,   \
                             tr_y_yy_xzzz,  \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyyy,  \
                             tr_y_yy_yyyz,  \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yyzz,  \
                             tr_y_yy_yzz,   \
                             tr_y_yy_yzzz,  \
                             tr_y_yy_zzz,   \
                             tr_y_yy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_xxxx[i] = 4.0 * tr_y_yy_xxx[i] * fe_0 + tr_y_yy_xxxx[i] * pa_x[i];

        tr_y_xyy_xxxy[i] = 3.0 * tr_y_yy_xxy[i] * fe_0 + tr_y_yy_xxxy[i] * pa_x[i];

        tr_y_xyy_xxxz[i] = 3.0 * tr_y_yy_xxz[i] * fe_0 + tr_y_yy_xxxz[i] * pa_x[i];

        tr_y_xyy_xxyy[i] = 2.0 * tr_y_yy_xyy[i] * fe_0 + tr_y_yy_xxyy[i] * pa_x[i];

        tr_y_xyy_xxyz[i] = 2.0 * tr_y_yy_xyz[i] * fe_0 + tr_y_yy_xxyz[i] * pa_x[i];

        tr_y_xyy_xxzz[i] = 2.0 * tr_y_yy_xzz[i] * fe_0 + tr_y_yy_xxzz[i] * pa_x[i];

        tr_y_xyy_xyyy[i] = tr_y_yy_yyy[i] * fe_0 + tr_y_yy_xyyy[i] * pa_x[i];

        tr_y_xyy_xyyz[i] = tr_y_yy_yyz[i] * fe_0 + tr_y_yy_xyyz[i] * pa_x[i];

        tr_y_xyy_xyzz[i] = tr_y_yy_yzz[i] * fe_0 + tr_y_yy_xyzz[i] * pa_x[i];

        tr_y_xyy_xzzz[i] = tr_y_yy_zzz[i] * fe_0 + tr_y_yy_xzzz[i] * pa_x[i];

        tr_y_xyy_yyyy[i] = tr_y_yy_yyyy[i] * pa_x[i];

        tr_y_xyy_yyyz[i] = tr_y_yy_yyyz[i] * pa_x[i];

        tr_y_xyy_yyzz[i] = tr_y_yy_yyzz[i] * pa_x[i];

        tr_y_xyy_yzzz[i] = tr_y_yy_yzzz[i] * pa_x[i];

        tr_y_xyy_zzzz[i] = tr_y_yy_zzzz[i] * pa_x[i];
    }

    // Set up 210-225 components of targeted buffer : FG

    auto tr_y_xyz_xxxx = pbuffer.data(idx_dip_fg + 210);

    auto tr_y_xyz_xxxy = pbuffer.data(idx_dip_fg + 211);

    auto tr_y_xyz_xxxz = pbuffer.data(idx_dip_fg + 212);

    auto tr_y_xyz_xxyy = pbuffer.data(idx_dip_fg + 213);

    auto tr_y_xyz_xxyz = pbuffer.data(idx_dip_fg + 214);

    auto tr_y_xyz_xxzz = pbuffer.data(idx_dip_fg + 215);

    auto tr_y_xyz_xyyy = pbuffer.data(idx_dip_fg + 216);

    auto tr_y_xyz_xyyz = pbuffer.data(idx_dip_fg + 217);

    auto tr_y_xyz_xyzz = pbuffer.data(idx_dip_fg + 218);

    auto tr_y_xyz_xzzz = pbuffer.data(idx_dip_fg + 219);

    auto tr_y_xyz_yyyy = pbuffer.data(idx_dip_fg + 220);

    auto tr_y_xyz_yyyz = pbuffer.data(idx_dip_fg + 221);

    auto tr_y_xyz_yyzz = pbuffer.data(idx_dip_fg + 222);

    auto tr_y_xyz_yzzz = pbuffer.data(idx_dip_fg + 223);

    auto tr_y_xyz_zzzz = pbuffer.data(idx_dip_fg + 224);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xy_xxxx,  \
                             tr_y_xy_xxxy,  \
                             tr_y_xy_xxyy,  \
                             tr_y_xy_xyyy,  \
                             tr_y_xyz_xxxx, \
                             tr_y_xyz_xxxy, \
                             tr_y_xyz_xxxz, \
                             tr_y_xyz_xxyy, \
                             tr_y_xyz_xxyz, \
                             tr_y_xyz_xxzz, \
                             tr_y_xyz_xyyy, \
                             tr_y_xyz_xyyz, \
                             tr_y_xyz_xyzz, \
                             tr_y_xyz_xzzz, \
                             tr_y_xyz_yyyy, \
                             tr_y_xyz_yyyz, \
                             tr_y_xyz_yyzz, \
                             tr_y_xyz_yzzz, \
                             tr_y_xyz_zzzz, \
                             tr_y_yz_xxxz,  \
                             tr_y_yz_xxyz,  \
                             tr_y_yz_xxz,   \
                             tr_y_yz_xxzz,  \
                             tr_y_yz_xyyz,  \
                             tr_y_yz_xyz,   \
                             tr_y_yz_xyzz,  \
                             tr_y_yz_xzz,   \
                             tr_y_yz_xzzz,  \
                             tr_y_yz_yyyy,  \
                             tr_y_yz_yyyz,  \
                             tr_y_yz_yyz,   \
                             tr_y_yz_yyzz,  \
                             tr_y_yz_yzz,   \
                             tr_y_yz_yzzz,  \
                             tr_y_yz_zzz,   \
                             tr_y_yz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyz_xxxx[i] = tr_y_xy_xxxx[i] * pa_z[i];

        tr_y_xyz_xxxy[i] = tr_y_xy_xxxy[i] * pa_z[i];

        tr_y_xyz_xxxz[i] = 3.0 * tr_y_yz_xxz[i] * fe_0 + tr_y_yz_xxxz[i] * pa_x[i];

        tr_y_xyz_xxyy[i] = tr_y_xy_xxyy[i] * pa_z[i];

        tr_y_xyz_xxyz[i] = 2.0 * tr_y_yz_xyz[i] * fe_0 + tr_y_yz_xxyz[i] * pa_x[i];

        tr_y_xyz_xxzz[i] = 2.0 * tr_y_yz_xzz[i] * fe_0 + tr_y_yz_xxzz[i] * pa_x[i];

        tr_y_xyz_xyyy[i] = tr_y_xy_xyyy[i] * pa_z[i];

        tr_y_xyz_xyyz[i] = tr_y_yz_yyz[i] * fe_0 + tr_y_yz_xyyz[i] * pa_x[i];

        tr_y_xyz_xyzz[i] = tr_y_yz_yzz[i] * fe_0 + tr_y_yz_xyzz[i] * pa_x[i];

        tr_y_xyz_xzzz[i] = tr_y_yz_zzz[i] * fe_0 + tr_y_yz_xzzz[i] * pa_x[i];

        tr_y_xyz_yyyy[i] = tr_y_yz_yyyy[i] * pa_x[i];

        tr_y_xyz_yyyz[i] = tr_y_yz_yyyz[i] * pa_x[i];

        tr_y_xyz_yyzz[i] = tr_y_yz_yyzz[i] * pa_x[i];

        tr_y_xyz_yzzz[i] = tr_y_yz_yzzz[i] * pa_x[i];

        tr_y_xyz_zzzz[i] = tr_y_yz_zzzz[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : FG

    auto tr_y_xzz_xxxx = pbuffer.data(idx_dip_fg + 225);

    auto tr_y_xzz_xxxy = pbuffer.data(idx_dip_fg + 226);

    auto tr_y_xzz_xxxz = pbuffer.data(idx_dip_fg + 227);

    auto tr_y_xzz_xxyy = pbuffer.data(idx_dip_fg + 228);

    auto tr_y_xzz_xxyz = pbuffer.data(idx_dip_fg + 229);

    auto tr_y_xzz_xxzz = pbuffer.data(idx_dip_fg + 230);

    auto tr_y_xzz_xyyy = pbuffer.data(idx_dip_fg + 231);

    auto tr_y_xzz_xyyz = pbuffer.data(idx_dip_fg + 232);

    auto tr_y_xzz_xyzz = pbuffer.data(idx_dip_fg + 233);

    auto tr_y_xzz_xzzz = pbuffer.data(idx_dip_fg + 234);

    auto tr_y_xzz_yyyy = pbuffer.data(idx_dip_fg + 235);

    auto tr_y_xzz_yyyz = pbuffer.data(idx_dip_fg + 236);

    auto tr_y_xzz_yyzz = pbuffer.data(idx_dip_fg + 237);

    auto tr_y_xzz_yzzz = pbuffer.data(idx_dip_fg + 238);

    auto tr_y_xzz_zzzz = pbuffer.data(idx_dip_fg + 239);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xzz_xxxx, \
                             tr_y_xzz_xxxy, \
                             tr_y_xzz_xxxz, \
                             tr_y_xzz_xxyy, \
                             tr_y_xzz_xxyz, \
                             tr_y_xzz_xxzz, \
                             tr_y_xzz_xyyy, \
                             tr_y_xzz_xyyz, \
                             tr_y_xzz_xyzz, \
                             tr_y_xzz_xzzz, \
                             tr_y_xzz_yyyy, \
                             tr_y_xzz_yyyz, \
                             tr_y_xzz_yyzz, \
                             tr_y_xzz_yzzz, \
                             tr_y_xzz_zzzz, \
                             tr_y_zz_xxx,   \
                             tr_y_zz_xxxx,  \
                             tr_y_zz_xxxy,  \
                             tr_y_zz_xxxz,  \
                             tr_y_zz_xxy,   \
                             tr_y_zz_xxyy,  \
                             tr_y_zz_xxyz,  \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xxzz,  \
                             tr_y_zz_xyy,   \
                             tr_y_zz_xyyy,  \
                             tr_y_zz_xyyz,  \
                             tr_y_zz_xyz,   \
                             tr_y_zz_xyzz,  \
                             tr_y_zz_xzz,   \
                             tr_y_zz_xzzz,  \
                             tr_y_zz_yyy,   \
                             tr_y_zz_yyyy,  \
                             tr_y_zz_yyyz,  \
                             tr_y_zz_yyz,   \
                             tr_y_zz_yyzz,  \
                             tr_y_zz_yzz,   \
                             tr_y_zz_yzzz,  \
                             tr_y_zz_zzz,   \
                             tr_y_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_xxxx[i] = 4.0 * tr_y_zz_xxx[i] * fe_0 + tr_y_zz_xxxx[i] * pa_x[i];

        tr_y_xzz_xxxy[i] = 3.0 * tr_y_zz_xxy[i] * fe_0 + tr_y_zz_xxxy[i] * pa_x[i];

        tr_y_xzz_xxxz[i] = 3.0 * tr_y_zz_xxz[i] * fe_0 + tr_y_zz_xxxz[i] * pa_x[i];

        tr_y_xzz_xxyy[i] = 2.0 * tr_y_zz_xyy[i] * fe_0 + tr_y_zz_xxyy[i] * pa_x[i];

        tr_y_xzz_xxyz[i] = 2.0 * tr_y_zz_xyz[i] * fe_0 + tr_y_zz_xxyz[i] * pa_x[i];

        tr_y_xzz_xxzz[i] = 2.0 * tr_y_zz_xzz[i] * fe_0 + tr_y_zz_xxzz[i] * pa_x[i];

        tr_y_xzz_xyyy[i] = tr_y_zz_yyy[i] * fe_0 + tr_y_zz_xyyy[i] * pa_x[i];

        tr_y_xzz_xyyz[i] = tr_y_zz_yyz[i] * fe_0 + tr_y_zz_xyyz[i] * pa_x[i];

        tr_y_xzz_xyzz[i] = tr_y_zz_yzz[i] * fe_0 + tr_y_zz_xyzz[i] * pa_x[i];

        tr_y_xzz_xzzz[i] = tr_y_zz_zzz[i] * fe_0 + tr_y_zz_xzzz[i] * pa_x[i];

        tr_y_xzz_yyyy[i] = tr_y_zz_yyyy[i] * pa_x[i];

        tr_y_xzz_yyyz[i] = tr_y_zz_yyyz[i] * pa_x[i];

        tr_y_xzz_yyzz[i] = tr_y_zz_yyzz[i] * pa_x[i];

        tr_y_xzz_yzzz[i] = tr_y_zz_yzzz[i] * pa_x[i];

        tr_y_xzz_zzzz[i] = tr_y_zz_zzzz[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : FG

    auto tr_y_yyy_xxxx = pbuffer.data(idx_dip_fg + 240);

    auto tr_y_yyy_xxxy = pbuffer.data(idx_dip_fg + 241);

    auto tr_y_yyy_xxxz = pbuffer.data(idx_dip_fg + 242);

    auto tr_y_yyy_xxyy = pbuffer.data(idx_dip_fg + 243);

    auto tr_y_yyy_xxyz = pbuffer.data(idx_dip_fg + 244);

    auto tr_y_yyy_xxzz = pbuffer.data(idx_dip_fg + 245);

    auto tr_y_yyy_xyyy = pbuffer.data(idx_dip_fg + 246);

    auto tr_y_yyy_xyyz = pbuffer.data(idx_dip_fg + 247);

    auto tr_y_yyy_xyzz = pbuffer.data(idx_dip_fg + 248);

    auto tr_y_yyy_xzzz = pbuffer.data(idx_dip_fg + 249);

    auto tr_y_yyy_yyyy = pbuffer.data(idx_dip_fg + 250);

    auto tr_y_yyy_yyyz = pbuffer.data(idx_dip_fg + 251);

    auto tr_y_yyy_yyzz = pbuffer.data(idx_dip_fg + 252);

    auto tr_y_yyy_yzzz = pbuffer.data(idx_dip_fg + 253);

    auto tr_y_yyy_zzzz = pbuffer.data(idx_dip_fg + 254);

#pragma omp simd aligned(pa_y,              \
                             tr_y_y_xxxx,   \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxxz,   \
                             tr_y_y_xxyy,   \
                             tr_y_y_xxyz,   \
                             tr_y_y_xxzz,   \
                             tr_y_y_xyyy,   \
                             tr_y_y_xyyz,   \
                             tr_y_y_xyzz,   \
                             tr_y_y_xzzz,   \
                             tr_y_y_yyyy,   \
                             tr_y_y_yyyz,   \
                             tr_y_y_yyzz,   \
                             tr_y_y_yzzz,   \
                             tr_y_y_zzzz,   \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxxx,  \
                             tr_y_yy_xxxy,  \
                             tr_y_yy_xxxz,  \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xxyy,  \
                             tr_y_yy_xxyz,  \
                             tr_y_yy_xxz,   \
                             tr_y_yy_xxzz,  \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyyy,  \
                             tr_y_yy_xyyz,  \
                             tr_y_yy_xyz,   \
                             tr_y_yy_xyzz,  \
                             tr_y_yy_xzz,   \
                             tr_y_yy_xzzz,  \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyyy,  \
                             tr_y_yy_yyyz,  \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yyzz,  \
                             tr_y_yy_yzz,   \
                             tr_y_yy_yzzz,  \
                             tr_y_yy_zzz,   \
                             tr_y_yy_zzzz,  \
                             tr_y_yyy_xxxx, \
                             tr_y_yyy_xxxy, \
                             tr_y_yyy_xxxz, \
                             tr_y_yyy_xxyy, \
                             tr_y_yyy_xxyz, \
                             tr_y_yyy_xxzz, \
                             tr_y_yyy_xyyy, \
                             tr_y_yyy_xyyz, \
                             tr_y_yyy_xyzz, \
                             tr_y_yyy_xzzz, \
                             tr_y_yyy_yyyy, \
                             tr_y_yyy_yyyz, \
                             tr_y_yyy_yyzz, \
                             tr_y_yyy_yzzz, \
                             tr_y_yyy_zzzz, \
                             ts_yy_xxxx,    \
                             ts_yy_xxxy,    \
                             ts_yy_xxxz,    \
                             ts_yy_xxyy,    \
                             ts_yy_xxyz,    \
                             ts_yy_xxzz,    \
                             ts_yy_xyyy,    \
                             ts_yy_xyyz,    \
                             ts_yy_xyzz,    \
                             ts_yy_xzzz,    \
                             ts_yy_yyyy,    \
                             ts_yy_yyyz,    \
                             ts_yy_yyzz,    \
                             ts_yy_yzzz,    \
                             ts_yy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_xxxx[i] = 2.0 * tr_y_y_xxxx[i] * fe_0 + ts_yy_xxxx[i] * fe_0 + tr_y_yy_xxxx[i] * pa_y[i];

        tr_y_yyy_xxxy[i] = 2.0 * tr_y_y_xxxy[i] * fe_0 + tr_y_yy_xxx[i] * fe_0 + ts_yy_xxxy[i] * fe_0 + tr_y_yy_xxxy[i] * pa_y[i];

        tr_y_yyy_xxxz[i] = 2.0 * tr_y_y_xxxz[i] * fe_0 + ts_yy_xxxz[i] * fe_0 + tr_y_yy_xxxz[i] * pa_y[i];

        tr_y_yyy_xxyy[i] = 2.0 * tr_y_y_xxyy[i] * fe_0 + 2.0 * tr_y_yy_xxy[i] * fe_0 + ts_yy_xxyy[i] * fe_0 + tr_y_yy_xxyy[i] * pa_y[i];

        tr_y_yyy_xxyz[i] = 2.0 * tr_y_y_xxyz[i] * fe_0 + tr_y_yy_xxz[i] * fe_0 + ts_yy_xxyz[i] * fe_0 + tr_y_yy_xxyz[i] * pa_y[i];

        tr_y_yyy_xxzz[i] = 2.0 * tr_y_y_xxzz[i] * fe_0 + ts_yy_xxzz[i] * fe_0 + tr_y_yy_xxzz[i] * pa_y[i];

        tr_y_yyy_xyyy[i] = 2.0 * tr_y_y_xyyy[i] * fe_0 + 3.0 * tr_y_yy_xyy[i] * fe_0 + ts_yy_xyyy[i] * fe_0 + tr_y_yy_xyyy[i] * pa_y[i];

        tr_y_yyy_xyyz[i] = 2.0 * tr_y_y_xyyz[i] * fe_0 + 2.0 * tr_y_yy_xyz[i] * fe_0 + ts_yy_xyyz[i] * fe_0 + tr_y_yy_xyyz[i] * pa_y[i];

        tr_y_yyy_xyzz[i] = 2.0 * tr_y_y_xyzz[i] * fe_0 + tr_y_yy_xzz[i] * fe_0 + ts_yy_xyzz[i] * fe_0 + tr_y_yy_xyzz[i] * pa_y[i];

        tr_y_yyy_xzzz[i] = 2.0 * tr_y_y_xzzz[i] * fe_0 + ts_yy_xzzz[i] * fe_0 + tr_y_yy_xzzz[i] * pa_y[i];

        tr_y_yyy_yyyy[i] = 2.0 * tr_y_y_yyyy[i] * fe_0 + 4.0 * tr_y_yy_yyy[i] * fe_0 + ts_yy_yyyy[i] * fe_0 + tr_y_yy_yyyy[i] * pa_y[i];

        tr_y_yyy_yyyz[i] = 2.0 * tr_y_y_yyyz[i] * fe_0 + 3.0 * tr_y_yy_yyz[i] * fe_0 + ts_yy_yyyz[i] * fe_0 + tr_y_yy_yyyz[i] * pa_y[i];

        tr_y_yyy_yyzz[i] = 2.0 * tr_y_y_yyzz[i] * fe_0 + 2.0 * tr_y_yy_yzz[i] * fe_0 + ts_yy_yyzz[i] * fe_0 + tr_y_yy_yyzz[i] * pa_y[i];

        tr_y_yyy_yzzz[i] = 2.0 * tr_y_y_yzzz[i] * fe_0 + tr_y_yy_zzz[i] * fe_0 + ts_yy_yzzz[i] * fe_0 + tr_y_yy_yzzz[i] * pa_y[i];

        tr_y_yyy_zzzz[i] = 2.0 * tr_y_y_zzzz[i] * fe_0 + ts_yy_zzzz[i] * fe_0 + tr_y_yy_zzzz[i] * pa_y[i];
    }

    // Set up 255-270 components of targeted buffer : FG

    auto tr_y_yyz_xxxx = pbuffer.data(idx_dip_fg + 255);

    auto tr_y_yyz_xxxy = pbuffer.data(idx_dip_fg + 256);

    auto tr_y_yyz_xxxz = pbuffer.data(idx_dip_fg + 257);

    auto tr_y_yyz_xxyy = pbuffer.data(idx_dip_fg + 258);

    auto tr_y_yyz_xxyz = pbuffer.data(idx_dip_fg + 259);

    auto tr_y_yyz_xxzz = pbuffer.data(idx_dip_fg + 260);

    auto tr_y_yyz_xyyy = pbuffer.data(idx_dip_fg + 261);

    auto tr_y_yyz_xyyz = pbuffer.data(idx_dip_fg + 262);

    auto tr_y_yyz_xyzz = pbuffer.data(idx_dip_fg + 263);

    auto tr_y_yyz_xzzz = pbuffer.data(idx_dip_fg + 264);

    auto tr_y_yyz_yyyy = pbuffer.data(idx_dip_fg + 265);

    auto tr_y_yyz_yyyz = pbuffer.data(idx_dip_fg + 266);

    auto tr_y_yyz_yyzz = pbuffer.data(idx_dip_fg + 267);

    auto tr_y_yyz_yzzz = pbuffer.data(idx_dip_fg + 268);

    auto tr_y_yyz_zzzz = pbuffer.data(idx_dip_fg + 269);

#pragma omp simd aligned(pa_z,              \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxxx,  \
                             tr_y_yy_xxxy,  \
                             tr_y_yy_xxxz,  \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xxyy,  \
                             tr_y_yy_xxyz,  \
                             tr_y_yy_xxz,   \
                             tr_y_yy_xxzz,  \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyyy,  \
                             tr_y_yy_xyyz,  \
                             tr_y_yy_xyz,   \
                             tr_y_yy_xyzz,  \
                             tr_y_yy_xzz,   \
                             tr_y_yy_xzzz,  \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyyy,  \
                             tr_y_yy_yyyz,  \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yyzz,  \
                             tr_y_yy_yzz,   \
                             tr_y_yy_yzzz,  \
                             tr_y_yy_zzz,   \
                             tr_y_yy_zzzz,  \
                             tr_y_yyz_xxxx, \
                             tr_y_yyz_xxxy, \
                             tr_y_yyz_xxxz, \
                             tr_y_yyz_xxyy, \
                             tr_y_yyz_xxyz, \
                             tr_y_yyz_xxzz, \
                             tr_y_yyz_xyyy, \
                             tr_y_yyz_xyyz, \
                             tr_y_yyz_xyzz, \
                             tr_y_yyz_xzzz, \
                             tr_y_yyz_yyyy, \
                             tr_y_yyz_yyyz, \
                             tr_y_yyz_yyzz, \
                             tr_y_yyz_yzzz, \
                             tr_y_yyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_xxxx[i] = tr_y_yy_xxxx[i] * pa_z[i];

        tr_y_yyz_xxxy[i] = tr_y_yy_xxxy[i] * pa_z[i];

        tr_y_yyz_xxxz[i] = tr_y_yy_xxx[i] * fe_0 + tr_y_yy_xxxz[i] * pa_z[i];

        tr_y_yyz_xxyy[i] = tr_y_yy_xxyy[i] * pa_z[i];

        tr_y_yyz_xxyz[i] = tr_y_yy_xxy[i] * fe_0 + tr_y_yy_xxyz[i] * pa_z[i];

        tr_y_yyz_xxzz[i] = 2.0 * tr_y_yy_xxz[i] * fe_0 + tr_y_yy_xxzz[i] * pa_z[i];

        tr_y_yyz_xyyy[i] = tr_y_yy_xyyy[i] * pa_z[i];

        tr_y_yyz_xyyz[i] = tr_y_yy_xyy[i] * fe_0 + tr_y_yy_xyyz[i] * pa_z[i];

        tr_y_yyz_xyzz[i] = 2.0 * tr_y_yy_xyz[i] * fe_0 + tr_y_yy_xyzz[i] * pa_z[i];

        tr_y_yyz_xzzz[i] = 3.0 * tr_y_yy_xzz[i] * fe_0 + tr_y_yy_xzzz[i] * pa_z[i];

        tr_y_yyz_yyyy[i] = tr_y_yy_yyyy[i] * pa_z[i];

        tr_y_yyz_yyyz[i] = tr_y_yy_yyy[i] * fe_0 + tr_y_yy_yyyz[i] * pa_z[i];

        tr_y_yyz_yyzz[i] = 2.0 * tr_y_yy_yyz[i] * fe_0 + tr_y_yy_yyzz[i] * pa_z[i];

        tr_y_yyz_yzzz[i] = 3.0 * tr_y_yy_yzz[i] * fe_0 + tr_y_yy_yzzz[i] * pa_z[i];

        tr_y_yyz_zzzz[i] = 4.0 * tr_y_yy_zzz[i] * fe_0 + tr_y_yy_zzzz[i] * pa_z[i];
    }

    // Set up 270-285 components of targeted buffer : FG

    auto tr_y_yzz_xxxx = pbuffer.data(idx_dip_fg + 270);

    auto tr_y_yzz_xxxy = pbuffer.data(idx_dip_fg + 271);

    auto tr_y_yzz_xxxz = pbuffer.data(idx_dip_fg + 272);

    auto tr_y_yzz_xxyy = pbuffer.data(idx_dip_fg + 273);

    auto tr_y_yzz_xxyz = pbuffer.data(idx_dip_fg + 274);

    auto tr_y_yzz_xxzz = pbuffer.data(idx_dip_fg + 275);

    auto tr_y_yzz_xyyy = pbuffer.data(idx_dip_fg + 276);

    auto tr_y_yzz_xyyz = pbuffer.data(idx_dip_fg + 277);

    auto tr_y_yzz_xyzz = pbuffer.data(idx_dip_fg + 278);

    auto tr_y_yzz_xzzz = pbuffer.data(idx_dip_fg + 279);

    auto tr_y_yzz_yyyy = pbuffer.data(idx_dip_fg + 280);

    auto tr_y_yzz_yyyz = pbuffer.data(idx_dip_fg + 281);

    auto tr_y_yzz_yyzz = pbuffer.data(idx_dip_fg + 282);

    auto tr_y_yzz_yzzz = pbuffer.data(idx_dip_fg + 283);

    auto tr_y_yzz_zzzz = pbuffer.data(idx_dip_fg + 284);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxyy,   \
                             tr_y_y_xyyy,   \
                             tr_y_y_yyyy,   \
                             tr_y_yz_xxxy,  \
                             tr_y_yz_xxyy,  \
                             tr_y_yz_xyyy,  \
                             tr_y_yz_yyyy,  \
                             tr_y_yzz_xxxx, \
                             tr_y_yzz_xxxy, \
                             tr_y_yzz_xxxz, \
                             tr_y_yzz_xxyy, \
                             tr_y_yzz_xxyz, \
                             tr_y_yzz_xxzz, \
                             tr_y_yzz_xyyy, \
                             tr_y_yzz_xyyz, \
                             tr_y_yzz_xyzz, \
                             tr_y_yzz_xzzz, \
                             tr_y_yzz_yyyy, \
                             tr_y_yzz_yyyz, \
                             tr_y_yzz_yyzz, \
                             tr_y_yzz_yzzz, \
                             tr_y_yzz_zzzz, \
                             tr_y_zz_xxxx,  \
                             tr_y_zz_xxxz,  \
                             tr_y_zz_xxyz,  \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xxzz,  \
                             tr_y_zz_xyyz,  \
                             tr_y_zz_xyz,   \
                             tr_y_zz_xyzz,  \
                             tr_y_zz_xzz,   \
                             tr_y_zz_xzzz,  \
                             tr_y_zz_yyyz,  \
                             tr_y_zz_yyz,   \
                             tr_y_zz_yyzz,  \
                             tr_y_zz_yzz,   \
                             tr_y_zz_yzzz,  \
                             tr_y_zz_zzz,   \
                             tr_y_zz_zzzz,  \
                             ts_zz_xxxx,    \
                             ts_zz_xxxz,    \
                             ts_zz_xxyz,    \
                             ts_zz_xxzz,    \
                             ts_zz_xyyz,    \
                             ts_zz_xyzz,    \
                             ts_zz_xzzz,    \
                             ts_zz_yyyz,    \
                             ts_zz_yyzz,    \
                             ts_zz_yzzz,    \
                             ts_zz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_xxxx[i] = ts_zz_xxxx[i] * fe_0 + tr_y_zz_xxxx[i] * pa_y[i];

        tr_y_yzz_xxxy[i] = tr_y_y_xxxy[i] * fe_0 + tr_y_yz_xxxy[i] * pa_z[i];

        tr_y_yzz_xxxz[i] = ts_zz_xxxz[i] * fe_0 + tr_y_zz_xxxz[i] * pa_y[i];

        tr_y_yzz_xxyy[i] = tr_y_y_xxyy[i] * fe_0 + tr_y_yz_xxyy[i] * pa_z[i];

        tr_y_yzz_xxyz[i] = tr_y_zz_xxz[i] * fe_0 + ts_zz_xxyz[i] * fe_0 + tr_y_zz_xxyz[i] * pa_y[i];

        tr_y_yzz_xxzz[i] = ts_zz_xxzz[i] * fe_0 + tr_y_zz_xxzz[i] * pa_y[i];

        tr_y_yzz_xyyy[i] = tr_y_y_xyyy[i] * fe_0 + tr_y_yz_xyyy[i] * pa_z[i];

        tr_y_yzz_xyyz[i] = 2.0 * tr_y_zz_xyz[i] * fe_0 + ts_zz_xyyz[i] * fe_0 + tr_y_zz_xyyz[i] * pa_y[i];

        tr_y_yzz_xyzz[i] = tr_y_zz_xzz[i] * fe_0 + ts_zz_xyzz[i] * fe_0 + tr_y_zz_xyzz[i] * pa_y[i];

        tr_y_yzz_xzzz[i] = ts_zz_xzzz[i] * fe_0 + tr_y_zz_xzzz[i] * pa_y[i];

        tr_y_yzz_yyyy[i] = tr_y_y_yyyy[i] * fe_0 + tr_y_yz_yyyy[i] * pa_z[i];

        tr_y_yzz_yyyz[i] = 3.0 * tr_y_zz_yyz[i] * fe_0 + ts_zz_yyyz[i] * fe_0 + tr_y_zz_yyyz[i] * pa_y[i];

        tr_y_yzz_yyzz[i] = 2.0 * tr_y_zz_yzz[i] * fe_0 + ts_zz_yyzz[i] * fe_0 + tr_y_zz_yyzz[i] * pa_y[i];

        tr_y_yzz_yzzz[i] = tr_y_zz_zzz[i] * fe_0 + ts_zz_yzzz[i] * fe_0 + tr_y_zz_yzzz[i] * pa_y[i];

        tr_y_yzz_zzzz[i] = ts_zz_zzzz[i] * fe_0 + tr_y_zz_zzzz[i] * pa_y[i];
    }

    // Set up 285-300 components of targeted buffer : FG

    auto tr_y_zzz_xxxx = pbuffer.data(idx_dip_fg + 285);

    auto tr_y_zzz_xxxy = pbuffer.data(idx_dip_fg + 286);

    auto tr_y_zzz_xxxz = pbuffer.data(idx_dip_fg + 287);

    auto tr_y_zzz_xxyy = pbuffer.data(idx_dip_fg + 288);

    auto tr_y_zzz_xxyz = pbuffer.data(idx_dip_fg + 289);

    auto tr_y_zzz_xxzz = pbuffer.data(idx_dip_fg + 290);

    auto tr_y_zzz_xyyy = pbuffer.data(idx_dip_fg + 291);

    auto tr_y_zzz_xyyz = pbuffer.data(idx_dip_fg + 292);

    auto tr_y_zzz_xyzz = pbuffer.data(idx_dip_fg + 293);

    auto tr_y_zzz_xzzz = pbuffer.data(idx_dip_fg + 294);

    auto tr_y_zzz_yyyy = pbuffer.data(idx_dip_fg + 295);

    auto tr_y_zzz_yyyz = pbuffer.data(idx_dip_fg + 296);

    auto tr_y_zzz_yyzz = pbuffer.data(idx_dip_fg + 297);

    auto tr_y_zzz_yzzz = pbuffer.data(idx_dip_fg + 298);

    auto tr_y_zzz_zzzz = pbuffer.data(idx_dip_fg + 299);

#pragma omp simd aligned(pa_z,              \
                             tr_y_z_xxxx,   \
                             tr_y_z_xxxy,   \
                             tr_y_z_xxxz,   \
                             tr_y_z_xxyy,   \
                             tr_y_z_xxyz,   \
                             tr_y_z_xxzz,   \
                             tr_y_z_xyyy,   \
                             tr_y_z_xyyz,   \
                             tr_y_z_xyzz,   \
                             tr_y_z_xzzz,   \
                             tr_y_z_yyyy,   \
                             tr_y_z_yyyz,   \
                             tr_y_z_yyzz,   \
                             tr_y_z_yzzz,   \
                             tr_y_z_zzzz,   \
                             tr_y_zz_xxx,   \
                             tr_y_zz_xxxx,  \
                             tr_y_zz_xxxy,  \
                             tr_y_zz_xxxz,  \
                             tr_y_zz_xxy,   \
                             tr_y_zz_xxyy,  \
                             tr_y_zz_xxyz,  \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xxzz,  \
                             tr_y_zz_xyy,   \
                             tr_y_zz_xyyy,  \
                             tr_y_zz_xyyz,  \
                             tr_y_zz_xyz,   \
                             tr_y_zz_xyzz,  \
                             tr_y_zz_xzz,   \
                             tr_y_zz_xzzz,  \
                             tr_y_zz_yyy,   \
                             tr_y_zz_yyyy,  \
                             tr_y_zz_yyyz,  \
                             tr_y_zz_yyz,   \
                             tr_y_zz_yyzz,  \
                             tr_y_zz_yzz,   \
                             tr_y_zz_yzzz,  \
                             tr_y_zz_zzz,   \
                             tr_y_zz_zzzz,  \
                             tr_y_zzz_xxxx, \
                             tr_y_zzz_xxxy, \
                             tr_y_zzz_xxxz, \
                             tr_y_zzz_xxyy, \
                             tr_y_zzz_xxyz, \
                             tr_y_zzz_xxzz, \
                             tr_y_zzz_xyyy, \
                             tr_y_zzz_xyyz, \
                             tr_y_zzz_xyzz, \
                             tr_y_zzz_xzzz, \
                             tr_y_zzz_yyyy, \
                             tr_y_zzz_yyyz, \
                             tr_y_zzz_yyzz, \
                             tr_y_zzz_yzzz, \
                             tr_y_zzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_xxxx[i] = 2.0 * tr_y_z_xxxx[i] * fe_0 + tr_y_zz_xxxx[i] * pa_z[i];

        tr_y_zzz_xxxy[i] = 2.0 * tr_y_z_xxxy[i] * fe_0 + tr_y_zz_xxxy[i] * pa_z[i];

        tr_y_zzz_xxxz[i] = 2.0 * tr_y_z_xxxz[i] * fe_0 + tr_y_zz_xxx[i] * fe_0 + tr_y_zz_xxxz[i] * pa_z[i];

        tr_y_zzz_xxyy[i] = 2.0 * tr_y_z_xxyy[i] * fe_0 + tr_y_zz_xxyy[i] * pa_z[i];

        tr_y_zzz_xxyz[i] = 2.0 * tr_y_z_xxyz[i] * fe_0 + tr_y_zz_xxy[i] * fe_0 + tr_y_zz_xxyz[i] * pa_z[i];

        tr_y_zzz_xxzz[i] = 2.0 * tr_y_z_xxzz[i] * fe_0 + 2.0 * tr_y_zz_xxz[i] * fe_0 + tr_y_zz_xxzz[i] * pa_z[i];

        tr_y_zzz_xyyy[i] = 2.0 * tr_y_z_xyyy[i] * fe_0 + tr_y_zz_xyyy[i] * pa_z[i];

        tr_y_zzz_xyyz[i] = 2.0 * tr_y_z_xyyz[i] * fe_0 + tr_y_zz_xyy[i] * fe_0 + tr_y_zz_xyyz[i] * pa_z[i];

        tr_y_zzz_xyzz[i] = 2.0 * tr_y_z_xyzz[i] * fe_0 + 2.0 * tr_y_zz_xyz[i] * fe_0 + tr_y_zz_xyzz[i] * pa_z[i];

        tr_y_zzz_xzzz[i] = 2.0 * tr_y_z_xzzz[i] * fe_0 + 3.0 * tr_y_zz_xzz[i] * fe_0 + tr_y_zz_xzzz[i] * pa_z[i];

        tr_y_zzz_yyyy[i] = 2.0 * tr_y_z_yyyy[i] * fe_0 + tr_y_zz_yyyy[i] * pa_z[i];

        tr_y_zzz_yyyz[i] = 2.0 * tr_y_z_yyyz[i] * fe_0 + tr_y_zz_yyy[i] * fe_0 + tr_y_zz_yyyz[i] * pa_z[i];

        tr_y_zzz_yyzz[i] = 2.0 * tr_y_z_yyzz[i] * fe_0 + 2.0 * tr_y_zz_yyz[i] * fe_0 + tr_y_zz_yyzz[i] * pa_z[i];

        tr_y_zzz_yzzz[i] = 2.0 * tr_y_z_yzzz[i] * fe_0 + 3.0 * tr_y_zz_yzz[i] * fe_0 + tr_y_zz_yzzz[i] * pa_z[i];

        tr_y_zzz_zzzz[i] = 2.0 * tr_y_z_zzzz[i] * fe_0 + 4.0 * tr_y_zz_zzz[i] * fe_0 + tr_y_zz_zzzz[i] * pa_z[i];
    }

    // Set up 300-315 components of targeted buffer : FG

    auto tr_z_xxx_xxxx = pbuffer.data(idx_dip_fg + 300);

    auto tr_z_xxx_xxxy = pbuffer.data(idx_dip_fg + 301);

    auto tr_z_xxx_xxxz = pbuffer.data(idx_dip_fg + 302);

    auto tr_z_xxx_xxyy = pbuffer.data(idx_dip_fg + 303);

    auto tr_z_xxx_xxyz = pbuffer.data(idx_dip_fg + 304);

    auto tr_z_xxx_xxzz = pbuffer.data(idx_dip_fg + 305);

    auto tr_z_xxx_xyyy = pbuffer.data(idx_dip_fg + 306);

    auto tr_z_xxx_xyyz = pbuffer.data(idx_dip_fg + 307);

    auto tr_z_xxx_xyzz = pbuffer.data(idx_dip_fg + 308);

    auto tr_z_xxx_xzzz = pbuffer.data(idx_dip_fg + 309);

    auto tr_z_xxx_yyyy = pbuffer.data(idx_dip_fg + 310);

    auto tr_z_xxx_yyyz = pbuffer.data(idx_dip_fg + 311);

    auto tr_z_xxx_yyzz = pbuffer.data(idx_dip_fg + 312);

    auto tr_z_xxx_yzzz = pbuffer.data(idx_dip_fg + 313);

    auto tr_z_xxx_zzzz = pbuffer.data(idx_dip_fg + 314);

#pragma omp simd aligned(pa_x,              \
                             tr_z_x_xxxx,   \
                             tr_z_x_xxxy,   \
                             tr_z_x_xxxz,   \
                             tr_z_x_xxyy,   \
                             tr_z_x_xxyz,   \
                             tr_z_x_xxzz,   \
                             tr_z_x_xyyy,   \
                             tr_z_x_xyyz,   \
                             tr_z_x_xyzz,   \
                             tr_z_x_xzzz,   \
                             tr_z_x_yyyy,   \
                             tr_z_x_yyyz,   \
                             tr_z_x_yyzz,   \
                             tr_z_x_yzzz,   \
                             tr_z_x_zzzz,   \
                             tr_z_xx_xxx,   \
                             tr_z_xx_xxxx,  \
                             tr_z_xx_xxxy,  \
                             tr_z_xx_xxxz,  \
                             tr_z_xx_xxy,   \
                             tr_z_xx_xxyy,  \
                             tr_z_xx_xxyz,  \
                             tr_z_xx_xxz,   \
                             tr_z_xx_xxzz,  \
                             tr_z_xx_xyy,   \
                             tr_z_xx_xyyy,  \
                             tr_z_xx_xyyz,  \
                             tr_z_xx_xyz,   \
                             tr_z_xx_xyzz,  \
                             tr_z_xx_xzz,   \
                             tr_z_xx_xzzz,  \
                             tr_z_xx_yyy,   \
                             tr_z_xx_yyyy,  \
                             tr_z_xx_yyyz,  \
                             tr_z_xx_yyz,   \
                             tr_z_xx_yyzz,  \
                             tr_z_xx_yzz,   \
                             tr_z_xx_yzzz,  \
                             tr_z_xx_zzz,   \
                             tr_z_xx_zzzz,  \
                             tr_z_xxx_xxxx, \
                             tr_z_xxx_xxxy, \
                             tr_z_xxx_xxxz, \
                             tr_z_xxx_xxyy, \
                             tr_z_xxx_xxyz, \
                             tr_z_xxx_xxzz, \
                             tr_z_xxx_xyyy, \
                             tr_z_xxx_xyyz, \
                             tr_z_xxx_xyzz, \
                             tr_z_xxx_xzzz, \
                             tr_z_xxx_yyyy, \
                             tr_z_xxx_yyyz, \
                             tr_z_xxx_yyzz, \
                             tr_z_xxx_yzzz, \
                             tr_z_xxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_xxxx[i] = 2.0 * tr_z_x_xxxx[i] * fe_0 + 4.0 * tr_z_xx_xxx[i] * fe_0 + tr_z_xx_xxxx[i] * pa_x[i];

        tr_z_xxx_xxxy[i] = 2.0 * tr_z_x_xxxy[i] * fe_0 + 3.0 * tr_z_xx_xxy[i] * fe_0 + tr_z_xx_xxxy[i] * pa_x[i];

        tr_z_xxx_xxxz[i] = 2.0 * tr_z_x_xxxz[i] * fe_0 + 3.0 * tr_z_xx_xxz[i] * fe_0 + tr_z_xx_xxxz[i] * pa_x[i];

        tr_z_xxx_xxyy[i] = 2.0 * tr_z_x_xxyy[i] * fe_0 + 2.0 * tr_z_xx_xyy[i] * fe_0 + tr_z_xx_xxyy[i] * pa_x[i];

        tr_z_xxx_xxyz[i] = 2.0 * tr_z_x_xxyz[i] * fe_0 + 2.0 * tr_z_xx_xyz[i] * fe_0 + tr_z_xx_xxyz[i] * pa_x[i];

        tr_z_xxx_xxzz[i] = 2.0 * tr_z_x_xxzz[i] * fe_0 + 2.0 * tr_z_xx_xzz[i] * fe_0 + tr_z_xx_xxzz[i] * pa_x[i];

        tr_z_xxx_xyyy[i] = 2.0 * tr_z_x_xyyy[i] * fe_0 + tr_z_xx_yyy[i] * fe_0 + tr_z_xx_xyyy[i] * pa_x[i];

        tr_z_xxx_xyyz[i] = 2.0 * tr_z_x_xyyz[i] * fe_0 + tr_z_xx_yyz[i] * fe_0 + tr_z_xx_xyyz[i] * pa_x[i];

        tr_z_xxx_xyzz[i] = 2.0 * tr_z_x_xyzz[i] * fe_0 + tr_z_xx_yzz[i] * fe_0 + tr_z_xx_xyzz[i] * pa_x[i];

        tr_z_xxx_xzzz[i] = 2.0 * tr_z_x_xzzz[i] * fe_0 + tr_z_xx_zzz[i] * fe_0 + tr_z_xx_xzzz[i] * pa_x[i];

        tr_z_xxx_yyyy[i] = 2.0 * tr_z_x_yyyy[i] * fe_0 + tr_z_xx_yyyy[i] * pa_x[i];

        tr_z_xxx_yyyz[i] = 2.0 * tr_z_x_yyyz[i] * fe_0 + tr_z_xx_yyyz[i] * pa_x[i];

        tr_z_xxx_yyzz[i] = 2.0 * tr_z_x_yyzz[i] * fe_0 + tr_z_xx_yyzz[i] * pa_x[i];

        tr_z_xxx_yzzz[i] = 2.0 * tr_z_x_yzzz[i] * fe_0 + tr_z_xx_yzzz[i] * pa_x[i];

        tr_z_xxx_zzzz[i] = 2.0 * tr_z_x_zzzz[i] * fe_0 + tr_z_xx_zzzz[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : FG

    auto tr_z_xxy_xxxx = pbuffer.data(idx_dip_fg + 315);

    auto tr_z_xxy_xxxy = pbuffer.data(idx_dip_fg + 316);

    auto tr_z_xxy_xxxz = pbuffer.data(idx_dip_fg + 317);

    auto tr_z_xxy_xxyy = pbuffer.data(idx_dip_fg + 318);

    auto tr_z_xxy_xxyz = pbuffer.data(idx_dip_fg + 319);

    auto tr_z_xxy_xxzz = pbuffer.data(idx_dip_fg + 320);

    auto tr_z_xxy_xyyy = pbuffer.data(idx_dip_fg + 321);

    auto tr_z_xxy_xyyz = pbuffer.data(idx_dip_fg + 322);

    auto tr_z_xxy_xyzz = pbuffer.data(idx_dip_fg + 323);

    auto tr_z_xxy_xzzz = pbuffer.data(idx_dip_fg + 324);

    auto tr_z_xxy_yyyy = pbuffer.data(idx_dip_fg + 325);

    auto tr_z_xxy_yyyz = pbuffer.data(idx_dip_fg + 326);

    auto tr_z_xxy_yyzz = pbuffer.data(idx_dip_fg + 327);

    auto tr_z_xxy_yzzz = pbuffer.data(idx_dip_fg + 328);

    auto tr_z_xxy_zzzz = pbuffer.data(idx_dip_fg + 329);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xx_xxx,   \
                             tr_z_xx_xxxx,  \
                             tr_z_xx_xxxy,  \
                             tr_z_xx_xxxz,  \
                             tr_z_xx_xxy,   \
                             tr_z_xx_xxyy,  \
                             tr_z_xx_xxyz,  \
                             tr_z_xx_xxz,   \
                             tr_z_xx_xxzz,  \
                             tr_z_xx_xyy,   \
                             tr_z_xx_xyyy,  \
                             tr_z_xx_xyyz,  \
                             tr_z_xx_xyz,   \
                             tr_z_xx_xyzz,  \
                             tr_z_xx_xzz,   \
                             tr_z_xx_xzzz,  \
                             tr_z_xx_zzzz,  \
                             tr_z_xxy_xxxx, \
                             tr_z_xxy_xxxy, \
                             tr_z_xxy_xxxz, \
                             tr_z_xxy_xxyy, \
                             tr_z_xxy_xxyz, \
                             tr_z_xxy_xxzz, \
                             tr_z_xxy_xyyy, \
                             tr_z_xxy_xyyz, \
                             tr_z_xxy_xyzz, \
                             tr_z_xxy_xzzz, \
                             tr_z_xxy_yyyy, \
                             tr_z_xxy_yyyz, \
                             tr_z_xxy_yyzz, \
                             tr_z_xxy_yzzz, \
                             tr_z_xxy_zzzz, \
                             tr_z_xy_yyyy,  \
                             tr_z_xy_yyyz,  \
                             tr_z_xy_yyzz,  \
                             tr_z_xy_yzzz,  \
                             tr_z_y_yyyy,   \
                             tr_z_y_yyyz,   \
                             tr_z_y_yyzz,   \
                             tr_z_y_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_xxxx[i] = tr_z_xx_xxxx[i] * pa_y[i];

        tr_z_xxy_xxxy[i] = tr_z_xx_xxx[i] * fe_0 + tr_z_xx_xxxy[i] * pa_y[i];

        tr_z_xxy_xxxz[i] = tr_z_xx_xxxz[i] * pa_y[i];

        tr_z_xxy_xxyy[i] = 2.0 * tr_z_xx_xxy[i] * fe_0 + tr_z_xx_xxyy[i] * pa_y[i];

        tr_z_xxy_xxyz[i] = tr_z_xx_xxz[i] * fe_0 + tr_z_xx_xxyz[i] * pa_y[i];

        tr_z_xxy_xxzz[i] = tr_z_xx_xxzz[i] * pa_y[i];

        tr_z_xxy_xyyy[i] = 3.0 * tr_z_xx_xyy[i] * fe_0 + tr_z_xx_xyyy[i] * pa_y[i];

        tr_z_xxy_xyyz[i] = 2.0 * tr_z_xx_xyz[i] * fe_0 + tr_z_xx_xyyz[i] * pa_y[i];

        tr_z_xxy_xyzz[i] = tr_z_xx_xzz[i] * fe_0 + tr_z_xx_xyzz[i] * pa_y[i];

        tr_z_xxy_xzzz[i] = tr_z_xx_xzzz[i] * pa_y[i];

        tr_z_xxy_yyyy[i] = tr_z_y_yyyy[i] * fe_0 + tr_z_xy_yyyy[i] * pa_x[i];

        tr_z_xxy_yyyz[i] = tr_z_y_yyyz[i] * fe_0 + tr_z_xy_yyyz[i] * pa_x[i];

        tr_z_xxy_yyzz[i] = tr_z_y_yyzz[i] * fe_0 + tr_z_xy_yyzz[i] * pa_x[i];

        tr_z_xxy_yzzz[i] = tr_z_y_yzzz[i] * fe_0 + tr_z_xy_yzzz[i] * pa_x[i];

        tr_z_xxy_zzzz[i] = tr_z_xx_zzzz[i] * pa_y[i];
    }

    // Set up 330-345 components of targeted buffer : FG

    auto tr_z_xxz_xxxx = pbuffer.data(idx_dip_fg + 330);

    auto tr_z_xxz_xxxy = pbuffer.data(idx_dip_fg + 331);

    auto tr_z_xxz_xxxz = pbuffer.data(idx_dip_fg + 332);

    auto tr_z_xxz_xxyy = pbuffer.data(idx_dip_fg + 333);

    auto tr_z_xxz_xxyz = pbuffer.data(idx_dip_fg + 334);

    auto tr_z_xxz_xxzz = pbuffer.data(idx_dip_fg + 335);

    auto tr_z_xxz_xyyy = pbuffer.data(idx_dip_fg + 336);

    auto tr_z_xxz_xyyz = pbuffer.data(idx_dip_fg + 337);

    auto tr_z_xxz_xyzz = pbuffer.data(idx_dip_fg + 338);

    auto tr_z_xxz_xzzz = pbuffer.data(idx_dip_fg + 339);

    auto tr_z_xxz_yyyy = pbuffer.data(idx_dip_fg + 340);

    auto tr_z_xxz_yyyz = pbuffer.data(idx_dip_fg + 341);

    auto tr_z_xxz_yyzz = pbuffer.data(idx_dip_fg + 342);

    auto tr_z_xxz_yzzz = pbuffer.data(idx_dip_fg + 343);

    auto tr_z_xxz_zzzz = pbuffer.data(idx_dip_fg + 344);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_z_xx_xxxx,  \
                             tr_z_xx_xxxy,  \
                             tr_z_xx_xxyy,  \
                             tr_z_xx_xyyy,  \
                             tr_z_xxz_xxxx, \
                             tr_z_xxz_xxxy, \
                             tr_z_xxz_xxxz, \
                             tr_z_xxz_xxyy, \
                             tr_z_xxz_xxyz, \
                             tr_z_xxz_xxzz, \
                             tr_z_xxz_xyyy, \
                             tr_z_xxz_xyyz, \
                             tr_z_xxz_xyzz, \
                             tr_z_xxz_xzzz, \
                             tr_z_xxz_yyyy, \
                             tr_z_xxz_yyyz, \
                             tr_z_xxz_yyzz, \
                             tr_z_xxz_yzzz, \
                             tr_z_xxz_zzzz, \
                             tr_z_xz_xxxz,  \
                             tr_z_xz_xxyz,  \
                             tr_z_xz_xxz,   \
                             tr_z_xz_xxzz,  \
                             tr_z_xz_xyyz,  \
                             tr_z_xz_xyz,   \
                             tr_z_xz_xyzz,  \
                             tr_z_xz_xzz,   \
                             tr_z_xz_xzzz,  \
                             tr_z_xz_yyyy,  \
                             tr_z_xz_yyyz,  \
                             tr_z_xz_yyz,   \
                             tr_z_xz_yyzz,  \
                             tr_z_xz_yzz,   \
                             tr_z_xz_yzzz,  \
                             tr_z_xz_zzz,   \
                             tr_z_xz_zzzz,  \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxzz,   \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyzz,   \
                             tr_z_z_xzzz,   \
                             tr_z_z_yyyy,   \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyzz,   \
                             tr_z_z_yzzz,   \
                             tr_z_z_zzzz,   \
                             ts_xx_xxxx,    \
                             ts_xx_xxxy,    \
                             ts_xx_xxyy,    \
                             ts_xx_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_xxxx[i] = ts_xx_xxxx[i] * fe_0 + tr_z_xx_xxxx[i] * pa_z[i];

        tr_z_xxz_xxxy[i] = ts_xx_xxxy[i] * fe_0 + tr_z_xx_xxxy[i] * pa_z[i];

        tr_z_xxz_xxxz[i] = tr_z_z_xxxz[i] * fe_0 + 3.0 * tr_z_xz_xxz[i] * fe_0 + tr_z_xz_xxxz[i] * pa_x[i];

        tr_z_xxz_xxyy[i] = ts_xx_xxyy[i] * fe_0 + tr_z_xx_xxyy[i] * pa_z[i];

        tr_z_xxz_xxyz[i] = tr_z_z_xxyz[i] * fe_0 + 2.0 * tr_z_xz_xyz[i] * fe_0 + tr_z_xz_xxyz[i] * pa_x[i];

        tr_z_xxz_xxzz[i] = tr_z_z_xxzz[i] * fe_0 + 2.0 * tr_z_xz_xzz[i] * fe_0 + tr_z_xz_xxzz[i] * pa_x[i];

        tr_z_xxz_xyyy[i] = ts_xx_xyyy[i] * fe_0 + tr_z_xx_xyyy[i] * pa_z[i];

        tr_z_xxz_xyyz[i] = tr_z_z_xyyz[i] * fe_0 + tr_z_xz_yyz[i] * fe_0 + tr_z_xz_xyyz[i] * pa_x[i];

        tr_z_xxz_xyzz[i] = tr_z_z_xyzz[i] * fe_0 + tr_z_xz_yzz[i] * fe_0 + tr_z_xz_xyzz[i] * pa_x[i];

        tr_z_xxz_xzzz[i] = tr_z_z_xzzz[i] * fe_0 + tr_z_xz_zzz[i] * fe_0 + tr_z_xz_xzzz[i] * pa_x[i];

        tr_z_xxz_yyyy[i] = tr_z_z_yyyy[i] * fe_0 + tr_z_xz_yyyy[i] * pa_x[i];

        tr_z_xxz_yyyz[i] = tr_z_z_yyyz[i] * fe_0 + tr_z_xz_yyyz[i] * pa_x[i];

        tr_z_xxz_yyzz[i] = tr_z_z_yyzz[i] * fe_0 + tr_z_xz_yyzz[i] * pa_x[i];

        tr_z_xxz_yzzz[i] = tr_z_z_yzzz[i] * fe_0 + tr_z_xz_yzzz[i] * pa_x[i];

        tr_z_xxz_zzzz[i] = tr_z_z_zzzz[i] * fe_0 + tr_z_xz_zzzz[i] * pa_x[i];
    }

    // Set up 345-360 components of targeted buffer : FG

    auto tr_z_xyy_xxxx = pbuffer.data(idx_dip_fg + 345);

    auto tr_z_xyy_xxxy = pbuffer.data(idx_dip_fg + 346);

    auto tr_z_xyy_xxxz = pbuffer.data(idx_dip_fg + 347);

    auto tr_z_xyy_xxyy = pbuffer.data(idx_dip_fg + 348);

    auto tr_z_xyy_xxyz = pbuffer.data(idx_dip_fg + 349);

    auto tr_z_xyy_xxzz = pbuffer.data(idx_dip_fg + 350);

    auto tr_z_xyy_xyyy = pbuffer.data(idx_dip_fg + 351);

    auto tr_z_xyy_xyyz = pbuffer.data(idx_dip_fg + 352);

    auto tr_z_xyy_xyzz = pbuffer.data(idx_dip_fg + 353);

    auto tr_z_xyy_xzzz = pbuffer.data(idx_dip_fg + 354);

    auto tr_z_xyy_yyyy = pbuffer.data(idx_dip_fg + 355);

    auto tr_z_xyy_yyyz = pbuffer.data(idx_dip_fg + 356);

    auto tr_z_xyy_yyzz = pbuffer.data(idx_dip_fg + 357);

    auto tr_z_xyy_yzzz = pbuffer.data(idx_dip_fg + 358);

    auto tr_z_xyy_zzzz = pbuffer.data(idx_dip_fg + 359);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xyy_xxxx, \
                             tr_z_xyy_xxxy, \
                             tr_z_xyy_xxxz, \
                             tr_z_xyy_xxyy, \
                             tr_z_xyy_xxyz, \
                             tr_z_xyy_xxzz, \
                             tr_z_xyy_xyyy, \
                             tr_z_xyy_xyyz, \
                             tr_z_xyy_xyzz, \
                             tr_z_xyy_xzzz, \
                             tr_z_xyy_yyyy, \
                             tr_z_xyy_yyyz, \
                             tr_z_xyy_yyzz, \
                             tr_z_xyy_yzzz, \
                             tr_z_xyy_zzzz, \
                             tr_z_yy_xxx,   \
                             tr_z_yy_xxxx,  \
                             tr_z_yy_xxxy,  \
                             tr_z_yy_xxxz,  \
                             tr_z_yy_xxy,   \
                             tr_z_yy_xxyy,  \
                             tr_z_yy_xxyz,  \
                             tr_z_yy_xxz,   \
                             tr_z_yy_xxzz,  \
                             tr_z_yy_xyy,   \
                             tr_z_yy_xyyy,  \
                             tr_z_yy_xyyz,  \
                             tr_z_yy_xyz,   \
                             tr_z_yy_xyzz,  \
                             tr_z_yy_xzz,   \
                             tr_z_yy_xzzz,  \
                             tr_z_yy_yyy,   \
                             tr_z_yy_yyyy,  \
                             tr_z_yy_yyyz,  \
                             tr_z_yy_yyz,   \
                             tr_z_yy_yyzz,  \
                             tr_z_yy_yzz,   \
                             tr_z_yy_yzzz,  \
                             tr_z_yy_zzz,   \
                             tr_z_yy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_xxxx[i] = 4.0 * tr_z_yy_xxx[i] * fe_0 + tr_z_yy_xxxx[i] * pa_x[i];

        tr_z_xyy_xxxy[i] = 3.0 * tr_z_yy_xxy[i] * fe_0 + tr_z_yy_xxxy[i] * pa_x[i];

        tr_z_xyy_xxxz[i] = 3.0 * tr_z_yy_xxz[i] * fe_0 + tr_z_yy_xxxz[i] * pa_x[i];

        tr_z_xyy_xxyy[i] = 2.0 * tr_z_yy_xyy[i] * fe_0 + tr_z_yy_xxyy[i] * pa_x[i];

        tr_z_xyy_xxyz[i] = 2.0 * tr_z_yy_xyz[i] * fe_0 + tr_z_yy_xxyz[i] * pa_x[i];

        tr_z_xyy_xxzz[i] = 2.0 * tr_z_yy_xzz[i] * fe_0 + tr_z_yy_xxzz[i] * pa_x[i];

        tr_z_xyy_xyyy[i] = tr_z_yy_yyy[i] * fe_0 + tr_z_yy_xyyy[i] * pa_x[i];

        tr_z_xyy_xyyz[i] = tr_z_yy_yyz[i] * fe_0 + tr_z_yy_xyyz[i] * pa_x[i];

        tr_z_xyy_xyzz[i] = tr_z_yy_yzz[i] * fe_0 + tr_z_yy_xyzz[i] * pa_x[i];

        tr_z_xyy_xzzz[i] = tr_z_yy_zzz[i] * fe_0 + tr_z_yy_xzzz[i] * pa_x[i];

        tr_z_xyy_yyyy[i] = tr_z_yy_yyyy[i] * pa_x[i];

        tr_z_xyy_yyyz[i] = tr_z_yy_yyyz[i] * pa_x[i];

        tr_z_xyy_yyzz[i] = tr_z_yy_yyzz[i] * pa_x[i];

        tr_z_xyy_yzzz[i] = tr_z_yy_yzzz[i] * pa_x[i];

        tr_z_xyy_zzzz[i] = tr_z_yy_zzzz[i] * pa_x[i];
    }

    // Set up 360-375 components of targeted buffer : FG

    auto tr_z_xyz_xxxx = pbuffer.data(idx_dip_fg + 360);

    auto tr_z_xyz_xxxy = pbuffer.data(idx_dip_fg + 361);

    auto tr_z_xyz_xxxz = pbuffer.data(idx_dip_fg + 362);

    auto tr_z_xyz_xxyy = pbuffer.data(idx_dip_fg + 363);

    auto tr_z_xyz_xxyz = pbuffer.data(idx_dip_fg + 364);

    auto tr_z_xyz_xxzz = pbuffer.data(idx_dip_fg + 365);

    auto tr_z_xyz_xyyy = pbuffer.data(idx_dip_fg + 366);

    auto tr_z_xyz_xyyz = pbuffer.data(idx_dip_fg + 367);

    auto tr_z_xyz_xyzz = pbuffer.data(idx_dip_fg + 368);

    auto tr_z_xyz_xzzz = pbuffer.data(idx_dip_fg + 369);

    auto tr_z_xyz_yyyy = pbuffer.data(idx_dip_fg + 370);

    auto tr_z_xyz_yyyz = pbuffer.data(idx_dip_fg + 371);

    auto tr_z_xyz_yyzz = pbuffer.data(idx_dip_fg + 372);

    auto tr_z_xyz_yzzz = pbuffer.data(idx_dip_fg + 373);

    auto tr_z_xyz_zzzz = pbuffer.data(idx_dip_fg + 374);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xyz_xxxx, \
                             tr_z_xyz_xxxy, \
                             tr_z_xyz_xxxz, \
                             tr_z_xyz_xxyy, \
                             tr_z_xyz_xxyz, \
                             tr_z_xyz_xxzz, \
                             tr_z_xyz_xyyy, \
                             tr_z_xyz_xyyz, \
                             tr_z_xyz_xyzz, \
                             tr_z_xyz_xzzz, \
                             tr_z_xyz_yyyy, \
                             tr_z_xyz_yyyz, \
                             tr_z_xyz_yyzz, \
                             tr_z_xyz_yzzz, \
                             tr_z_xyz_zzzz, \
                             tr_z_xz_xxxx,  \
                             tr_z_xz_xxxz,  \
                             tr_z_xz_xxzz,  \
                             tr_z_xz_xzzz,  \
                             tr_z_yz_xxxy,  \
                             tr_z_yz_xxy,   \
                             tr_z_yz_xxyy,  \
                             tr_z_yz_xxyz,  \
                             tr_z_yz_xyy,   \
                             tr_z_yz_xyyy,  \
                             tr_z_yz_xyyz,  \
                             tr_z_yz_xyz,   \
                             tr_z_yz_xyzz,  \
                             tr_z_yz_yyy,   \
                             tr_z_yz_yyyy,  \
                             tr_z_yz_yyyz,  \
                             tr_z_yz_yyz,   \
                             tr_z_yz_yyzz,  \
                             tr_z_yz_yzz,   \
                             tr_z_yz_yzzz,  \
                             tr_z_yz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyz_xxxx[i] = tr_z_xz_xxxx[i] * pa_y[i];

        tr_z_xyz_xxxy[i] = 3.0 * tr_z_yz_xxy[i] * fe_0 + tr_z_yz_xxxy[i] * pa_x[i];

        tr_z_xyz_xxxz[i] = tr_z_xz_xxxz[i] * pa_y[i];

        tr_z_xyz_xxyy[i] = 2.0 * tr_z_yz_xyy[i] * fe_0 + tr_z_yz_xxyy[i] * pa_x[i];

        tr_z_xyz_xxyz[i] = 2.0 * tr_z_yz_xyz[i] * fe_0 + tr_z_yz_xxyz[i] * pa_x[i];

        tr_z_xyz_xxzz[i] = tr_z_xz_xxzz[i] * pa_y[i];

        tr_z_xyz_xyyy[i] = tr_z_yz_yyy[i] * fe_0 + tr_z_yz_xyyy[i] * pa_x[i];

        tr_z_xyz_xyyz[i] = tr_z_yz_yyz[i] * fe_0 + tr_z_yz_xyyz[i] * pa_x[i];

        tr_z_xyz_xyzz[i] = tr_z_yz_yzz[i] * fe_0 + tr_z_yz_xyzz[i] * pa_x[i];

        tr_z_xyz_xzzz[i] = tr_z_xz_xzzz[i] * pa_y[i];

        tr_z_xyz_yyyy[i] = tr_z_yz_yyyy[i] * pa_x[i];

        tr_z_xyz_yyyz[i] = tr_z_yz_yyyz[i] * pa_x[i];

        tr_z_xyz_yyzz[i] = tr_z_yz_yyzz[i] * pa_x[i];

        tr_z_xyz_yzzz[i] = tr_z_yz_yzzz[i] * pa_x[i];

        tr_z_xyz_zzzz[i] = tr_z_yz_zzzz[i] * pa_x[i];
    }

    // Set up 375-390 components of targeted buffer : FG

    auto tr_z_xzz_xxxx = pbuffer.data(idx_dip_fg + 375);

    auto tr_z_xzz_xxxy = pbuffer.data(idx_dip_fg + 376);

    auto tr_z_xzz_xxxz = pbuffer.data(idx_dip_fg + 377);

    auto tr_z_xzz_xxyy = pbuffer.data(idx_dip_fg + 378);

    auto tr_z_xzz_xxyz = pbuffer.data(idx_dip_fg + 379);

    auto tr_z_xzz_xxzz = pbuffer.data(idx_dip_fg + 380);

    auto tr_z_xzz_xyyy = pbuffer.data(idx_dip_fg + 381);

    auto tr_z_xzz_xyyz = pbuffer.data(idx_dip_fg + 382);

    auto tr_z_xzz_xyzz = pbuffer.data(idx_dip_fg + 383);

    auto tr_z_xzz_xzzz = pbuffer.data(idx_dip_fg + 384);

    auto tr_z_xzz_yyyy = pbuffer.data(idx_dip_fg + 385);

    auto tr_z_xzz_yyyz = pbuffer.data(idx_dip_fg + 386);

    auto tr_z_xzz_yyzz = pbuffer.data(idx_dip_fg + 387);

    auto tr_z_xzz_yzzz = pbuffer.data(idx_dip_fg + 388);

    auto tr_z_xzz_zzzz = pbuffer.data(idx_dip_fg + 389);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xzz_xxxx, \
                             tr_z_xzz_xxxy, \
                             tr_z_xzz_xxxz, \
                             tr_z_xzz_xxyy, \
                             tr_z_xzz_xxyz, \
                             tr_z_xzz_xxzz, \
                             tr_z_xzz_xyyy, \
                             tr_z_xzz_xyyz, \
                             tr_z_xzz_xyzz, \
                             tr_z_xzz_xzzz, \
                             tr_z_xzz_yyyy, \
                             tr_z_xzz_yyyz, \
                             tr_z_xzz_yyzz, \
                             tr_z_xzz_yzzz, \
                             tr_z_xzz_zzzz, \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxxx,  \
                             tr_z_zz_xxxy,  \
                             tr_z_zz_xxxz,  \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxyy,  \
                             tr_z_zz_xxyz,  \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xxzz,  \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyyy,  \
                             tr_z_zz_xyyz,  \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xyzz,  \
                             tr_z_zz_xzz,   \
                             tr_z_zz_xzzz,  \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyyy,  \
                             tr_z_zz_yyyz,  \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yyzz,  \
                             tr_z_zz_yzz,   \
                             tr_z_zz_yzzz,  \
                             tr_z_zz_zzz,   \
                             tr_z_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_xxxx[i] = 4.0 * tr_z_zz_xxx[i] * fe_0 + tr_z_zz_xxxx[i] * pa_x[i];

        tr_z_xzz_xxxy[i] = 3.0 * tr_z_zz_xxy[i] * fe_0 + tr_z_zz_xxxy[i] * pa_x[i];

        tr_z_xzz_xxxz[i] = 3.0 * tr_z_zz_xxz[i] * fe_0 + tr_z_zz_xxxz[i] * pa_x[i];

        tr_z_xzz_xxyy[i] = 2.0 * tr_z_zz_xyy[i] * fe_0 + tr_z_zz_xxyy[i] * pa_x[i];

        tr_z_xzz_xxyz[i] = 2.0 * tr_z_zz_xyz[i] * fe_0 + tr_z_zz_xxyz[i] * pa_x[i];

        tr_z_xzz_xxzz[i] = 2.0 * tr_z_zz_xzz[i] * fe_0 + tr_z_zz_xxzz[i] * pa_x[i];

        tr_z_xzz_xyyy[i] = tr_z_zz_yyy[i] * fe_0 + tr_z_zz_xyyy[i] * pa_x[i];

        tr_z_xzz_xyyz[i] = tr_z_zz_yyz[i] * fe_0 + tr_z_zz_xyyz[i] * pa_x[i];

        tr_z_xzz_xyzz[i] = tr_z_zz_yzz[i] * fe_0 + tr_z_zz_xyzz[i] * pa_x[i];

        tr_z_xzz_xzzz[i] = tr_z_zz_zzz[i] * fe_0 + tr_z_zz_xzzz[i] * pa_x[i];

        tr_z_xzz_yyyy[i] = tr_z_zz_yyyy[i] * pa_x[i];

        tr_z_xzz_yyyz[i] = tr_z_zz_yyyz[i] * pa_x[i];

        tr_z_xzz_yyzz[i] = tr_z_zz_yyzz[i] * pa_x[i];

        tr_z_xzz_yzzz[i] = tr_z_zz_yzzz[i] * pa_x[i];

        tr_z_xzz_zzzz[i] = tr_z_zz_zzzz[i] * pa_x[i];
    }

    // Set up 390-405 components of targeted buffer : FG

    auto tr_z_yyy_xxxx = pbuffer.data(idx_dip_fg + 390);

    auto tr_z_yyy_xxxy = pbuffer.data(idx_dip_fg + 391);

    auto tr_z_yyy_xxxz = pbuffer.data(idx_dip_fg + 392);

    auto tr_z_yyy_xxyy = pbuffer.data(idx_dip_fg + 393);

    auto tr_z_yyy_xxyz = pbuffer.data(idx_dip_fg + 394);

    auto tr_z_yyy_xxzz = pbuffer.data(idx_dip_fg + 395);

    auto tr_z_yyy_xyyy = pbuffer.data(idx_dip_fg + 396);

    auto tr_z_yyy_xyyz = pbuffer.data(idx_dip_fg + 397);

    auto tr_z_yyy_xyzz = pbuffer.data(idx_dip_fg + 398);

    auto tr_z_yyy_xzzz = pbuffer.data(idx_dip_fg + 399);

    auto tr_z_yyy_yyyy = pbuffer.data(idx_dip_fg + 400);

    auto tr_z_yyy_yyyz = pbuffer.data(idx_dip_fg + 401);

    auto tr_z_yyy_yyzz = pbuffer.data(idx_dip_fg + 402);

    auto tr_z_yyy_yzzz = pbuffer.data(idx_dip_fg + 403);

    auto tr_z_yyy_zzzz = pbuffer.data(idx_dip_fg + 404);

#pragma omp simd aligned(pa_y,              \
                             tr_z_y_xxxx,   \
                             tr_z_y_xxxy,   \
                             tr_z_y_xxxz,   \
                             tr_z_y_xxyy,   \
                             tr_z_y_xxyz,   \
                             tr_z_y_xxzz,   \
                             tr_z_y_xyyy,   \
                             tr_z_y_xyyz,   \
                             tr_z_y_xyzz,   \
                             tr_z_y_xzzz,   \
                             tr_z_y_yyyy,   \
                             tr_z_y_yyyz,   \
                             tr_z_y_yyzz,   \
                             tr_z_y_yzzz,   \
                             tr_z_y_zzzz,   \
                             tr_z_yy_xxx,   \
                             tr_z_yy_xxxx,  \
                             tr_z_yy_xxxy,  \
                             tr_z_yy_xxxz,  \
                             tr_z_yy_xxy,   \
                             tr_z_yy_xxyy,  \
                             tr_z_yy_xxyz,  \
                             tr_z_yy_xxz,   \
                             tr_z_yy_xxzz,  \
                             tr_z_yy_xyy,   \
                             tr_z_yy_xyyy,  \
                             tr_z_yy_xyyz,  \
                             tr_z_yy_xyz,   \
                             tr_z_yy_xyzz,  \
                             tr_z_yy_xzz,   \
                             tr_z_yy_xzzz,  \
                             tr_z_yy_yyy,   \
                             tr_z_yy_yyyy,  \
                             tr_z_yy_yyyz,  \
                             tr_z_yy_yyz,   \
                             tr_z_yy_yyzz,  \
                             tr_z_yy_yzz,   \
                             tr_z_yy_yzzz,  \
                             tr_z_yy_zzz,   \
                             tr_z_yy_zzzz,  \
                             tr_z_yyy_xxxx, \
                             tr_z_yyy_xxxy, \
                             tr_z_yyy_xxxz, \
                             tr_z_yyy_xxyy, \
                             tr_z_yyy_xxyz, \
                             tr_z_yyy_xxzz, \
                             tr_z_yyy_xyyy, \
                             tr_z_yyy_xyyz, \
                             tr_z_yyy_xyzz, \
                             tr_z_yyy_xzzz, \
                             tr_z_yyy_yyyy, \
                             tr_z_yyy_yyyz, \
                             tr_z_yyy_yyzz, \
                             tr_z_yyy_yzzz, \
                             tr_z_yyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_xxxx[i] = 2.0 * tr_z_y_xxxx[i] * fe_0 + tr_z_yy_xxxx[i] * pa_y[i];

        tr_z_yyy_xxxy[i] = 2.0 * tr_z_y_xxxy[i] * fe_0 + tr_z_yy_xxx[i] * fe_0 + tr_z_yy_xxxy[i] * pa_y[i];

        tr_z_yyy_xxxz[i] = 2.0 * tr_z_y_xxxz[i] * fe_0 + tr_z_yy_xxxz[i] * pa_y[i];

        tr_z_yyy_xxyy[i] = 2.0 * tr_z_y_xxyy[i] * fe_0 + 2.0 * tr_z_yy_xxy[i] * fe_0 + tr_z_yy_xxyy[i] * pa_y[i];

        tr_z_yyy_xxyz[i] = 2.0 * tr_z_y_xxyz[i] * fe_0 + tr_z_yy_xxz[i] * fe_0 + tr_z_yy_xxyz[i] * pa_y[i];

        tr_z_yyy_xxzz[i] = 2.0 * tr_z_y_xxzz[i] * fe_0 + tr_z_yy_xxzz[i] * pa_y[i];

        tr_z_yyy_xyyy[i] = 2.0 * tr_z_y_xyyy[i] * fe_0 + 3.0 * tr_z_yy_xyy[i] * fe_0 + tr_z_yy_xyyy[i] * pa_y[i];

        tr_z_yyy_xyyz[i] = 2.0 * tr_z_y_xyyz[i] * fe_0 + 2.0 * tr_z_yy_xyz[i] * fe_0 + tr_z_yy_xyyz[i] * pa_y[i];

        tr_z_yyy_xyzz[i] = 2.0 * tr_z_y_xyzz[i] * fe_0 + tr_z_yy_xzz[i] * fe_0 + tr_z_yy_xyzz[i] * pa_y[i];

        tr_z_yyy_xzzz[i] = 2.0 * tr_z_y_xzzz[i] * fe_0 + tr_z_yy_xzzz[i] * pa_y[i];

        tr_z_yyy_yyyy[i] = 2.0 * tr_z_y_yyyy[i] * fe_0 + 4.0 * tr_z_yy_yyy[i] * fe_0 + tr_z_yy_yyyy[i] * pa_y[i];

        tr_z_yyy_yyyz[i] = 2.0 * tr_z_y_yyyz[i] * fe_0 + 3.0 * tr_z_yy_yyz[i] * fe_0 + tr_z_yy_yyyz[i] * pa_y[i];

        tr_z_yyy_yyzz[i] = 2.0 * tr_z_y_yyzz[i] * fe_0 + 2.0 * tr_z_yy_yzz[i] * fe_0 + tr_z_yy_yyzz[i] * pa_y[i];

        tr_z_yyy_yzzz[i] = 2.0 * tr_z_y_yzzz[i] * fe_0 + tr_z_yy_zzz[i] * fe_0 + tr_z_yy_yzzz[i] * pa_y[i];

        tr_z_yyy_zzzz[i] = 2.0 * tr_z_y_zzzz[i] * fe_0 + tr_z_yy_zzzz[i] * pa_y[i];
    }

    // Set up 405-420 components of targeted buffer : FG

    auto tr_z_yyz_xxxx = pbuffer.data(idx_dip_fg + 405);

    auto tr_z_yyz_xxxy = pbuffer.data(idx_dip_fg + 406);

    auto tr_z_yyz_xxxz = pbuffer.data(idx_dip_fg + 407);

    auto tr_z_yyz_xxyy = pbuffer.data(idx_dip_fg + 408);

    auto tr_z_yyz_xxyz = pbuffer.data(idx_dip_fg + 409);

    auto tr_z_yyz_xxzz = pbuffer.data(idx_dip_fg + 410);

    auto tr_z_yyz_xyyy = pbuffer.data(idx_dip_fg + 411);

    auto tr_z_yyz_xyyz = pbuffer.data(idx_dip_fg + 412);

    auto tr_z_yyz_xyzz = pbuffer.data(idx_dip_fg + 413);

    auto tr_z_yyz_xzzz = pbuffer.data(idx_dip_fg + 414);

    auto tr_z_yyz_yyyy = pbuffer.data(idx_dip_fg + 415);

    auto tr_z_yyz_yyyz = pbuffer.data(idx_dip_fg + 416);

    auto tr_z_yyz_yyzz = pbuffer.data(idx_dip_fg + 417);

    auto tr_z_yyz_yzzz = pbuffer.data(idx_dip_fg + 418);

    auto tr_z_yyz_zzzz = pbuffer.data(idx_dip_fg + 419);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_z_yy_xxxy,  \
                             tr_z_yy_xxyy,  \
                             tr_z_yy_xyyy,  \
                             tr_z_yy_yyyy,  \
                             tr_z_yyz_xxxx, \
                             tr_z_yyz_xxxy, \
                             tr_z_yyz_xxxz, \
                             tr_z_yyz_xxyy, \
                             tr_z_yyz_xxyz, \
                             tr_z_yyz_xxzz, \
                             tr_z_yyz_xyyy, \
                             tr_z_yyz_xyyz, \
                             tr_z_yyz_xyzz, \
                             tr_z_yyz_xzzz, \
                             tr_z_yyz_yyyy, \
                             tr_z_yyz_yyyz, \
                             tr_z_yyz_yyzz, \
                             tr_z_yyz_yzzz, \
                             tr_z_yyz_zzzz, \
                             tr_z_yz_xxxx,  \
                             tr_z_yz_xxxz,  \
                             tr_z_yz_xxyz,  \
                             tr_z_yz_xxz,   \
                             tr_z_yz_xxzz,  \
                             tr_z_yz_xyyz,  \
                             tr_z_yz_xyz,   \
                             tr_z_yz_xyzz,  \
                             tr_z_yz_xzz,   \
                             tr_z_yz_xzzz,  \
                             tr_z_yz_yyyz,  \
                             tr_z_yz_yyz,   \
                             tr_z_yz_yyzz,  \
                             tr_z_yz_yzz,   \
                             tr_z_yz_yzzz,  \
                             tr_z_yz_zzz,   \
                             tr_z_yz_zzzz,  \
                             tr_z_z_xxxx,   \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxzz,   \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyzz,   \
                             tr_z_z_xzzz,   \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyzz,   \
                             tr_z_z_yzzz,   \
                             tr_z_z_zzzz,   \
                             ts_yy_xxxy,    \
                             ts_yy_xxyy,    \
                             ts_yy_xyyy,    \
                             ts_yy_yyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_xxxx[i] = tr_z_z_xxxx[i] * fe_0 + tr_z_yz_xxxx[i] * pa_y[i];

        tr_z_yyz_xxxy[i] = ts_yy_xxxy[i] * fe_0 + tr_z_yy_xxxy[i] * pa_z[i];

        tr_z_yyz_xxxz[i] = tr_z_z_xxxz[i] * fe_0 + tr_z_yz_xxxz[i] * pa_y[i];

        tr_z_yyz_xxyy[i] = ts_yy_xxyy[i] * fe_0 + tr_z_yy_xxyy[i] * pa_z[i];

        tr_z_yyz_xxyz[i] = tr_z_z_xxyz[i] * fe_0 + tr_z_yz_xxz[i] * fe_0 + tr_z_yz_xxyz[i] * pa_y[i];

        tr_z_yyz_xxzz[i] = tr_z_z_xxzz[i] * fe_0 + tr_z_yz_xxzz[i] * pa_y[i];

        tr_z_yyz_xyyy[i] = ts_yy_xyyy[i] * fe_0 + tr_z_yy_xyyy[i] * pa_z[i];

        tr_z_yyz_xyyz[i] = tr_z_z_xyyz[i] * fe_0 + 2.0 * tr_z_yz_xyz[i] * fe_0 + tr_z_yz_xyyz[i] * pa_y[i];

        tr_z_yyz_xyzz[i] = tr_z_z_xyzz[i] * fe_0 + tr_z_yz_xzz[i] * fe_0 + tr_z_yz_xyzz[i] * pa_y[i];

        tr_z_yyz_xzzz[i] = tr_z_z_xzzz[i] * fe_0 + tr_z_yz_xzzz[i] * pa_y[i];

        tr_z_yyz_yyyy[i] = ts_yy_yyyy[i] * fe_0 + tr_z_yy_yyyy[i] * pa_z[i];

        tr_z_yyz_yyyz[i] = tr_z_z_yyyz[i] * fe_0 + 3.0 * tr_z_yz_yyz[i] * fe_0 + tr_z_yz_yyyz[i] * pa_y[i];

        tr_z_yyz_yyzz[i] = tr_z_z_yyzz[i] * fe_0 + 2.0 * tr_z_yz_yzz[i] * fe_0 + tr_z_yz_yyzz[i] * pa_y[i];

        tr_z_yyz_yzzz[i] = tr_z_z_yzzz[i] * fe_0 + tr_z_yz_zzz[i] * fe_0 + tr_z_yz_yzzz[i] * pa_y[i];

        tr_z_yyz_zzzz[i] = tr_z_z_zzzz[i] * fe_0 + tr_z_yz_zzzz[i] * pa_y[i];
    }

    // Set up 420-435 components of targeted buffer : FG

    auto tr_z_yzz_xxxx = pbuffer.data(idx_dip_fg + 420);

    auto tr_z_yzz_xxxy = pbuffer.data(idx_dip_fg + 421);

    auto tr_z_yzz_xxxz = pbuffer.data(idx_dip_fg + 422);

    auto tr_z_yzz_xxyy = pbuffer.data(idx_dip_fg + 423);

    auto tr_z_yzz_xxyz = pbuffer.data(idx_dip_fg + 424);

    auto tr_z_yzz_xxzz = pbuffer.data(idx_dip_fg + 425);

    auto tr_z_yzz_xyyy = pbuffer.data(idx_dip_fg + 426);

    auto tr_z_yzz_xyyz = pbuffer.data(idx_dip_fg + 427);

    auto tr_z_yzz_xyzz = pbuffer.data(idx_dip_fg + 428);

    auto tr_z_yzz_xzzz = pbuffer.data(idx_dip_fg + 429);

    auto tr_z_yzz_yyyy = pbuffer.data(idx_dip_fg + 430);

    auto tr_z_yzz_yyyz = pbuffer.data(idx_dip_fg + 431);

    auto tr_z_yzz_yyzz = pbuffer.data(idx_dip_fg + 432);

    auto tr_z_yzz_yzzz = pbuffer.data(idx_dip_fg + 433);

    auto tr_z_yzz_zzzz = pbuffer.data(idx_dip_fg + 434);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yzz_xxxx, \
                             tr_z_yzz_xxxy, \
                             tr_z_yzz_xxxz, \
                             tr_z_yzz_xxyy, \
                             tr_z_yzz_xxyz, \
                             tr_z_yzz_xxzz, \
                             tr_z_yzz_xyyy, \
                             tr_z_yzz_xyyz, \
                             tr_z_yzz_xyzz, \
                             tr_z_yzz_xzzz, \
                             tr_z_yzz_yyyy, \
                             tr_z_yzz_yyyz, \
                             tr_z_yzz_yyzz, \
                             tr_z_yzz_yzzz, \
                             tr_z_yzz_zzzz, \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxxx,  \
                             tr_z_zz_xxxy,  \
                             tr_z_zz_xxxz,  \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxyy,  \
                             tr_z_zz_xxyz,  \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xxzz,  \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyyy,  \
                             tr_z_zz_xyyz,  \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xyzz,  \
                             tr_z_zz_xzz,   \
                             tr_z_zz_xzzz,  \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyyy,  \
                             tr_z_zz_yyyz,  \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yyzz,  \
                             tr_z_zz_yzz,   \
                             tr_z_zz_yzzz,  \
                             tr_z_zz_zzz,   \
                             tr_z_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_xxxx[i] = tr_z_zz_xxxx[i] * pa_y[i];

        tr_z_yzz_xxxy[i] = tr_z_zz_xxx[i] * fe_0 + tr_z_zz_xxxy[i] * pa_y[i];

        tr_z_yzz_xxxz[i] = tr_z_zz_xxxz[i] * pa_y[i];

        tr_z_yzz_xxyy[i] = 2.0 * tr_z_zz_xxy[i] * fe_0 + tr_z_zz_xxyy[i] * pa_y[i];

        tr_z_yzz_xxyz[i] = tr_z_zz_xxz[i] * fe_0 + tr_z_zz_xxyz[i] * pa_y[i];

        tr_z_yzz_xxzz[i] = tr_z_zz_xxzz[i] * pa_y[i];

        tr_z_yzz_xyyy[i] = 3.0 * tr_z_zz_xyy[i] * fe_0 + tr_z_zz_xyyy[i] * pa_y[i];

        tr_z_yzz_xyyz[i] = 2.0 * tr_z_zz_xyz[i] * fe_0 + tr_z_zz_xyyz[i] * pa_y[i];

        tr_z_yzz_xyzz[i] = tr_z_zz_xzz[i] * fe_0 + tr_z_zz_xyzz[i] * pa_y[i];

        tr_z_yzz_xzzz[i] = tr_z_zz_xzzz[i] * pa_y[i];

        tr_z_yzz_yyyy[i] = 4.0 * tr_z_zz_yyy[i] * fe_0 + tr_z_zz_yyyy[i] * pa_y[i];

        tr_z_yzz_yyyz[i] = 3.0 * tr_z_zz_yyz[i] * fe_0 + tr_z_zz_yyyz[i] * pa_y[i];

        tr_z_yzz_yyzz[i] = 2.0 * tr_z_zz_yzz[i] * fe_0 + tr_z_zz_yyzz[i] * pa_y[i];

        tr_z_yzz_yzzz[i] = tr_z_zz_zzz[i] * fe_0 + tr_z_zz_yzzz[i] * pa_y[i];

        tr_z_yzz_zzzz[i] = tr_z_zz_zzzz[i] * pa_y[i];
    }

    // Set up 435-450 components of targeted buffer : FG

    auto tr_z_zzz_xxxx = pbuffer.data(idx_dip_fg + 435);

    auto tr_z_zzz_xxxy = pbuffer.data(idx_dip_fg + 436);

    auto tr_z_zzz_xxxz = pbuffer.data(idx_dip_fg + 437);

    auto tr_z_zzz_xxyy = pbuffer.data(idx_dip_fg + 438);

    auto tr_z_zzz_xxyz = pbuffer.data(idx_dip_fg + 439);

    auto tr_z_zzz_xxzz = pbuffer.data(idx_dip_fg + 440);

    auto tr_z_zzz_xyyy = pbuffer.data(idx_dip_fg + 441);

    auto tr_z_zzz_xyyz = pbuffer.data(idx_dip_fg + 442);

    auto tr_z_zzz_xyzz = pbuffer.data(idx_dip_fg + 443);

    auto tr_z_zzz_xzzz = pbuffer.data(idx_dip_fg + 444);

    auto tr_z_zzz_yyyy = pbuffer.data(idx_dip_fg + 445);

    auto tr_z_zzz_yyyz = pbuffer.data(idx_dip_fg + 446);

    auto tr_z_zzz_yyzz = pbuffer.data(idx_dip_fg + 447);

    auto tr_z_zzz_yzzz = pbuffer.data(idx_dip_fg + 448);

    auto tr_z_zzz_zzzz = pbuffer.data(idx_dip_fg + 449);

#pragma omp simd aligned(pa_z,              \
                             tr_z_z_xxxx,   \
                             tr_z_z_xxxy,   \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxyy,   \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxzz,   \
                             tr_z_z_xyyy,   \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyzz,   \
                             tr_z_z_xzzz,   \
                             tr_z_z_yyyy,   \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyzz,   \
                             tr_z_z_yzzz,   \
                             tr_z_z_zzzz,   \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxxx,  \
                             tr_z_zz_xxxy,  \
                             tr_z_zz_xxxz,  \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxyy,  \
                             tr_z_zz_xxyz,  \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xxzz,  \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyyy,  \
                             tr_z_zz_xyyz,  \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xyzz,  \
                             tr_z_zz_xzz,   \
                             tr_z_zz_xzzz,  \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyyy,  \
                             tr_z_zz_yyyz,  \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yyzz,  \
                             tr_z_zz_yzz,   \
                             tr_z_zz_yzzz,  \
                             tr_z_zz_zzz,   \
                             tr_z_zz_zzzz,  \
                             tr_z_zzz_xxxx, \
                             tr_z_zzz_xxxy, \
                             tr_z_zzz_xxxz, \
                             tr_z_zzz_xxyy, \
                             tr_z_zzz_xxyz, \
                             tr_z_zzz_xxzz, \
                             tr_z_zzz_xyyy, \
                             tr_z_zzz_xyyz, \
                             tr_z_zzz_xyzz, \
                             tr_z_zzz_xzzz, \
                             tr_z_zzz_yyyy, \
                             tr_z_zzz_yyyz, \
                             tr_z_zzz_yyzz, \
                             tr_z_zzz_yzzz, \
                             tr_z_zzz_zzzz, \
                             ts_zz_xxxx,    \
                             ts_zz_xxxy,    \
                             ts_zz_xxxz,    \
                             ts_zz_xxyy,    \
                             ts_zz_xxyz,    \
                             ts_zz_xxzz,    \
                             ts_zz_xyyy,    \
                             ts_zz_xyyz,    \
                             ts_zz_xyzz,    \
                             ts_zz_xzzz,    \
                             ts_zz_yyyy,    \
                             ts_zz_yyyz,    \
                             ts_zz_yyzz,    \
                             ts_zz_yzzz,    \
                             ts_zz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_xxxx[i] = 2.0 * tr_z_z_xxxx[i] * fe_0 + ts_zz_xxxx[i] * fe_0 + tr_z_zz_xxxx[i] * pa_z[i];

        tr_z_zzz_xxxy[i] = 2.0 * tr_z_z_xxxy[i] * fe_0 + ts_zz_xxxy[i] * fe_0 + tr_z_zz_xxxy[i] * pa_z[i];

        tr_z_zzz_xxxz[i] = 2.0 * tr_z_z_xxxz[i] * fe_0 + tr_z_zz_xxx[i] * fe_0 + ts_zz_xxxz[i] * fe_0 + tr_z_zz_xxxz[i] * pa_z[i];

        tr_z_zzz_xxyy[i] = 2.0 * tr_z_z_xxyy[i] * fe_0 + ts_zz_xxyy[i] * fe_0 + tr_z_zz_xxyy[i] * pa_z[i];

        tr_z_zzz_xxyz[i] = 2.0 * tr_z_z_xxyz[i] * fe_0 + tr_z_zz_xxy[i] * fe_0 + ts_zz_xxyz[i] * fe_0 + tr_z_zz_xxyz[i] * pa_z[i];

        tr_z_zzz_xxzz[i] = 2.0 * tr_z_z_xxzz[i] * fe_0 + 2.0 * tr_z_zz_xxz[i] * fe_0 + ts_zz_xxzz[i] * fe_0 + tr_z_zz_xxzz[i] * pa_z[i];

        tr_z_zzz_xyyy[i] = 2.0 * tr_z_z_xyyy[i] * fe_0 + ts_zz_xyyy[i] * fe_0 + tr_z_zz_xyyy[i] * pa_z[i];

        tr_z_zzz_xyyz[i] = 2.0 * tr_z_z_xyyz[i] * fe_0 + tr_z_zz_xyy[i] * fe_0 + ts_zz_xyyz[i] * fe_0 + tr_z_zz_xyyz[i] * pa_z[i];

        tr_z_zzz_xyzz[i] = 2.0 * tr_z_z_xyzz[i] * fe_0 + 2.0 * tr_z_zz_xyz[i] * fe_0 + ts_zz_xyzz[i] * fe_0 + tr_z_zz_xyzz[i] * pa_z[i];

        tr_z_zzz_xzzz[i] = 2.0 * tr_z_z_xzzz[i] * fe_0 + 3.0 * tr_z_zz_xzz[i] * fe_0 + ts_zz_xzzz[i] * fe_0 + tr_z_zz_xzzz[i] * pa_z[i];

        tr_z_zzz_yyyy[i] = 2.0 * tr_z_z_yyyy[i] * fe_0 + ts_zz_yyyy[i] * fe_0 + tr_z_zz_yyyy[i] * pa_z[i];

        tr_z_zzz_yyyz[i] = 2.0 * tr_z_z_yyyz[i] * fe_0 + tr_z_zz_yyy[i] * fe_0 + ts_zz_yyyz[i] * fe_0 + tr_z_zz_yyyz[i] * pa_z[i];

        tr_z_zzz_yyzz[i] = 2.0 * tr_z_z_yyzz[i] * fe_0 + 2.0 * tr_z_zz_yyz[i] * fe_0 + ts_zz_yyzz[i] * fe_0 + tr_z_zz_yyzz[i] * pa_z[i];

        tr_z_zzz_yzzz[i] = 2.0 * tr_z_z_yzzz[i] * fe_0 + 3.0 * tr_z_zz_yzz[i] * fe_0 + ts_zz_yzzz[i] * fe_0 + tr_z_zz_yzzz[i] * pa_z[i];

        tr_z_zzz_zzzz[i] = 2.0 * tr_z_z_zzzz[i] * fe_0 + 4.0 * tr_z_zz_zzz[i] * fe_0 + ts_zz_zzzz[i] * fe_0 + tr_z_zz_zzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
