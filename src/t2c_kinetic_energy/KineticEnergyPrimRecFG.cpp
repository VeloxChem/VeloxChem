#include "KineticEnergyPrimRecFG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_fg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_fg,
                            const size_t              idx_ovl_pg,
                            const size_t              idx_kin_pg,
                            const size_t              idx_kin_df,
                            const size_t              idx_kin_dg,
                            const size_t              idx_ovl_fg,
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

    auto ts_x_xxxx = pbuffer.data(idx_ovl_pg);

    auto ts_x_xxxy = pbuffer.data(idx_ovl_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_ovl_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_ovl_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_ovl_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_ovl_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_ovl_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_ovl_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_ovl_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_ovl_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_ovl_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_ovl_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_ovl_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_ovl_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_ovl_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_ovl_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_ovl_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_ovl_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_ovl_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_ovl_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_ovl_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_ovl_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_ovl_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_ovl_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_ovl_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_ovl_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_ovl_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_ovl_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_ovl_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_ovl_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_ovl_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_ovl_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_ovl_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_ovl_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_ovl_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_ovl_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_ovl_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_ovl_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_ovl_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_ovl_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_ovl_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_ovl_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_ovl_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_ovl_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_ovl_pg + 44);

    // Set up components of auxiliary buffer : PG

    auto tk_x_xxxx = pbuffer.data(idx_kin_pg);

    auto tk_x_xxxy = pbuffer.data(idx_kin_pg + 1);

    auto tk_x_xxxz = pbuffer.data(idx_kin_pg + 2);

    auto tk_x_xxyy = pbuffer.data(idx_kin_pg + 3);

    auto tk_x_xxyz = pbuffer.data(idx_kin_pg + 4);

    auto tk_x_xxzz = pbuffer.data(idx_kin_pg + 5);

    auto tk_x_xyyy = pbuffer.data(idx_kin_pg + 6);

    auto tk_x_xyyz = pbuffer.data(idx_kin_pg + 7);

    auto tk_x_xyzz = pbuffer.data(idx_kin_pg + 8);

    auto tk_x_xzzz = pbuffer.data(idx_kin_pg + 9);

    auto tk_x_yyyy = pbuffer.data(idx_kin_pg + 10);

    auto tk_x_yyyz = pbuffer.data(idx_kin_pg + 11);

    auto tk_x_yyzz = pbuffer.data(idx_kin_pg + 12);

    auto tk_x_yzzz = pbuffer.data(idx_kin_pg + 13);

    auto tk_x_zzzz = pbuffer.data(idx_kin_pg + 14);

    auto tk_y_xxxx = pbuffer.data(idx_kin_pg + 15);

    auto tk_y_xxxy = pbuffer.data(idx_kin_pg + 16);

    auto tk_y_xxxz = pbuffer.data(idx_kin_pg + 17);

    auto tk_y_xxyy = pbuffer.data(idx_kin_pg + 18);

    auto tk_y_xxyz = pbuffer.data(idx_kin_pg + 19);

    auto tk_y_xxzz = pbuffer.data(idx_kin_pg + 20);

    auto tk_y_xyyy = pbuffer.data(idx_kin_pg + 21);

    auto tk_y_xyyz = pbuffer.data(idx_kin_pg + 22);

    auto tk_y_xyzz = pbuffer.data(idx_kin_pg + 23);

    auto tk_y_xzzz = pbuffer.data(idx_kin_pg + 24);

    auto tk_y_yyyy = pbuffer.data(idx_kin_pg + 25);

    auto tk_y_yyyz = pbuffer.data(idx_kin_pg + 26);

    auto tk_y_yyzz = pbuffer.data(idx_kin_pg + 27);

    auto tk_y_yzzz = pbuffer.data(idx_kin_pg + 28);

    auto tk_y_zzzz = pbuffer.data(idx_kin_pg + 29);

    auto tk_z_xxxx = pbuffer.data(idx_kin_pg + 30);

    auto tk_z_xxxy = pbuffer.data(idx_kin_pg + 31);

    auto tk_z_xxxz = pbuffer.data(idx_kin_pg + 32);

    auto tk_z_xxyy = pbuffer.data(idx_kin_pg + 33);

    auto tk_z_xxyz = pbuffer.data(idx_kin_pg + 34);

    auto tk_z_xxzz = pbuffer.data(idx_kin_pg + 35);

    auto tk_z_xyyy = pbuffer.data(idx_kin_pg + 36);

    auto tk_z_xyyz = pbuffer.data(idx_kin_pg + 37);

    auto tk_z_xyzz = pbuffer.data(idx_kin_pg + 38);

    auto tk_z_xzzz = pbuffer.data(idx_kin_pg + 39);

    auto tk_z_yyyy = pbuffer.data(idx_kin_pg + 40);

    auto tk_z_yyyz = pbuffer.data(idx_kin_pg + 41);

    auto tk_z_yyzz = pbuffer.data(idx_kin_pg + 42);

    auto tk_z_yzzz = pbuffer.data(idx_kin_pg + 43);

    auto tk_z_zzzz = pbuffer.data(idx_kin_pg + 44);

    // Set up components of auxiliary buffer : DF

    auto tk_xx_xxx = pbuffer.data(idx_kin_df);

    auto tk_xx_xxy = pbuffer.data(idx_kin_df + 1);

    auto tk_xx_xxz = pbuffer.data(idx_kin_df + 2);

    auto tk_xx_xyy = pbuffer.data(idx_kin_df + 3);

    auto tk_xx_xyz = pbuffer.data(idx_kin_df + 4);

    auto tk_xx_xzz = pbuffer.data(idx_kin_df + 5);

    auto tk_xx_yyy = pbuffer.data(idx_kin_df + 6);

    auto tk_xx_yyz = pbuffer.data(idx_kin_df + 7);

    auto tk_xx_yzz = pbuffer.data(idx_kin_df + 8);

    auto tk_xx_zzz = pbuffer.data(idx_kin_df + 9);

    auto tk_yy_xxx = pbuffer.data(idx_kin_df + 30);

    auto tk_yy_xxy = pbuffer.data(idx_kin_df + 31);

    auto tk_yy_xxz = pbuffer.data(idx_kin_df + 32);

    auto tk_yy_xyy = pbuffer.data(idx_kin_df + 33);

    auto tk_yy_xyz = pbuffer.data(idx_kin_df + 34);

    auto tk_yy_xzz = pbuffer.data(idx_kin_df + 35);

    auto tk_yy_yyy = pbuffer.data(idx_kin_df + 36);

    auto tk_yy_yyz = pbuffer.data(idx_kin_df + 37);

    auto tk_yy_yzz = pbuffer.data(idx_kin_df + 38);

    auto tk_yy_zzz = pbuffer.data(idx_kin_df + 39);

    auto tk_yz_xyz = pbuffer.data(idx_kin_df + 44);

    auto tk_yz_yyz = pbuffer.data(idx_kin_df + 47);

    auto tk_yz_yzz = pbuffer.data(idx_kin_df + 48);

    auto tk_zz_xxx = pbuffer.data(idx_kin_df + 50);

    auto tk_zz_xxy = pbuffer.data(idx_kin_df + 51);

    auto tk_zz_xxz = pbuffer.data(idx_kin_df + 52);

    auto tk_zz_xyy = pbuffer.data(idx_kin_df + 53);

    auto tk_zz_xyz = pbuffer.data(idx_kin_df + 54);

    auto tk_zz_xzz = pbuffer.data(idx_kin_df + 55);

    auto tk_zz_yyy = pbuffer.data(idx_kin_df + 56);

    auto tk_zz_yyz = pbuffer.data(idx_kin_df + 57);

    auto tk_zz_yzz = pbuffer.data(idx_kin_df + 58);

    auto tk_zz_zzz = pbuffer.data(idx_kin_df + 59);

    // Set up components of auxiliary buffer : DG

    auto tk_xx_xxxx = pbuffer.data(idx_kin_dg);

    auto tk_xx_xxxy = pbuffer.data(idx_kin_dg + 1);

    auto tk_xx_xxxz = pbuffer.data(idx_kin_dg + 2);

    auto tk_xx_xxyy = pbuffer.data(idx_kin_dg + 3);

    auto tk_xx_xxyz = pbuffer.data(idx_kin_dg + 4);

    auto tk_xx_xxzz = pbuffer.data(idx_kin_dg + 5);

    auto tk_xx_xyyy = pbuffer.data(idx_kin_dg + 6);

    auto tk_xx_xyyz = pbuffer.data(idx_kin_dg + 7);

    auto tk_xx_xyzz = pbuffer.data(idx_kin_dg + 8);

    auto tk_xx_xzzz = pbuffer.data(idx_kin_dg + 9);

    auto tk_xx_yyyy = pbuffer.data(idx_kin_dg + 10);

    auto tk_xx_yyyz = pbuffer.data(idx_kin_dg + 11);

    auto tk_xx_yyzz = pbuffer.data(idx_kin_dg + 12);

    auto tk_xx_yzzz = pbuffer.data(idx_kin_dg + 13);

    auto tk_xx_zzzz = pbuffer.data(idx_kin_dg + 14);

    auto tk_xy_xxxy = pbuffer.data(idx_kin_dg + 16);

    auto tk_xy_xxyy = pbuffer.data(idx_kin_dg + 18);

    auto tk_xy_xyyy = pbuffer.data(idx_kin_dg + 21);

    auto tk_xz_xxxx = pbuffer.data(idx_kin_dg + 30);

    auto tk_xz_xxxz = pbuffer.data(idx_kin_dg + 32);

    auto tk_xz_xxzz = pbuffer.data(idx_kin_dg + 35);

    auto tk_xz_xzzz = pbuffer.data(idx_kin_dg + 39);

    auto tk_yy_xxxx = pbuffer.data(idx_kin_dg + 45);

    auto tk_yy_xxxy = pbuffer.data(idx_kin_dg + 46);

    auto tk_yy_xxxz = pbuffer.data(idx_kin_dg + 47);

    auto tk_yy_xxyy = pbuffer.data(idx_kin_dg + 48);

    auto tk_yy_xxyz = pbuffer.data(idx_kin_dg + 49);

    auto tk_yy_xxzz = pbuffer.data(idx_kin_dg + 50);

    auto tk_yy_xyyy = pbuffer.data(idx_kin_dg + 51);

    auto tk_yy_xyyz = pbuffer.data(idx_kin_dg + 52);

    auto tk_yy_xyzz = pbuffer.data(idx_kin_dg + 53);

    auto tk_yy_xzzz = pbuffer.data(idx_kin_dg + 54);

    auto tk_yy_yyyy = pbuffer.data(idx_kin_dg + 55);

    auto tk_yy_yyyz = pbuffer.data(idx_kin_dg + 56);

    auto tk_yy_yyzz = pbuffer.data(idx_kin_dg + 57);

    auto tk_yy_yzzz = pbuffer.data(idx_kin_dg + 58);

    auto tk_yy_zzzz = pbuffer.data(idx_kin_dg + 59);

    auto tk_yz_xxyz = pbuffer.data(idx_kin_dg + 64);

    auto tk_yz_xyyz = pbuffer.data(idx_kin_dg + 67);

    auto tk_yz_xyzz = pbuffer.data(idx_kin_dg + 68);

    auto tk_yz_yyyy = pbuffer.data(idx_kin_dg + 70);

    auto tk_yz_yyyz = pbuffer.data(idx_kin_dg + 71);

    auto tk_yz_yyzz = pbuffer.data(idx_kin_dg + 72);

    auto tk_yz_yzzz = pbuffer.data(idx_kin_dg + 73);

    auto tk_yz_zzzz = pbuffer.data(idx_kin_dg + 74);

    auto tk_zz_xxxx = pbuffer.data(idx_kin_dg + 75);

    auto tk_zz_xxxy = pbuffer.data(idx_kin_dg + 76);

    auto tk_zz_xxxz = pbuffer.data(idx_kin_dg + 77);

    auto tk_zz_xxyy = pbuffer.data(idx_kin_dg + 78);

    auto tk_zz_xxyz = pbuffer.data(idx_kin_dg + 79);

    auto tk_zz_xxzz = pbuffer.data(idx_kin_dg + 80);

    auto tk_zz_xyyy = pbuffer.data(idx_kin_dg + 81);

    auto tk_zz_xyyz = pbuffer.data(idx_kin_dg + 82);

    auto tk_zz_xyzz = pbuffer.data(idx_kin_dg + 83);

    auto tk_zz_xzzz = pbuffer.data(idx_kin_dg + 84);

    auto tk_zz_yyyy = pbuffer.data(idx_kin_dg + 85);

    auto tk_zz_yyyz = pbuffer.data(idx_kin_dg + 86);

    auto tk_zz_yyzz = pbuffer.data(idx_kin_dg + 87);

    auto tk_zz_yzzz = pbuffer.data(idx_kin_dg + 88);

    auto tk_zz_zzzz = pbuffer.data(idx_kin_dg + 89);

    // Set up components of auxiliary buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_ovl_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_ovl_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_ovl_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_ovl_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_ovl_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_ovl_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_ovl_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_ovl_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_ovl_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_ovl_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_ovl_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_ovl_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_ovl_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_ovl_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_ovl_fg + 14);

    auto ts_xxy_xxxx = pbuffer.data(idx_ovl_fg + 15);

    auto ts_xxy_xxxy = pbuffer.data(idx_ovl_fg + 16);

    auto ts_xxy_xxxz = pbuffer.data(idx_ovl_fg + 17);

    auto ts_xxy_xxyy = pbuffer.data(idx_ovl_fg + 18);

    auto ts_xxy_xxyz = pbuffer.data(idx_ovl_fg + 19);

    auto ts_xxy_xxzz = pbuffer.data(idx_ovl_fg + 20);

    auto ts_xxy_xyyy = pbuffer.data(idx_ovl_fg + 21);

    auto ts_xxy_xyyz = pbuffer.data(idx_ovl_fg + 22);

    auto ts_xxy_xyzz = pbuffer.data(idx_ovl_fg + 23);

    auto ts_xxy_xzzz = pbuffer.data(idx_ovl_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_ovl_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_ovl_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_ovl_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_ovl_fg + 28);

    auto ts_xxy_zzzz = pbuffer.data(idx_ovl_fg + 29);

    auto ts_xxz_xxxx = pbuffer.data(idx_ovl_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_ovl_fg + 31);

    auto ts_xxz_xxxz = pbuffer.data(idx_ovl_fg + 32);

    auto ts_xxz_xxyy = pbuffer.data(idx_ovl_fg + 33);

    auto ts_xxz_xxyz = pbuffer.data(idx_ovl_fg + 34);

    auto ts_xxz_xxzz = pbuffer.data(idx_ovl_fg + 35);

    auto ts_xxz_xyyy = pbuffer.data(idx_ovl_fg + 36);

    auto ts_xxz_xyyz = pbuffer.data(idx_ovl_fg + 37);

    auto ts_xxz_xyzz = pbuffer.data(idx_ovl_fg + 38);

    auto ts_xxz_xzzz = pbuffer.data(idx_ovl_fg + 39);

    auto ts_xxz_yyyy = pbuffer.data(idx_ovl_fg + 40);

    auto ts_xxz_yyyz = pbuffer.data(idx_ovl_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_ovl_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_ovl_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_ovl_fg + 44);

    auto ts_xyy_xxxx = pbuffer.data(idx_ovl_fg + 45);

    auto ts_xyy_xxxy = pbuffer.data(idx_ovl_fg + 46);

    auto ts_xyy_xxxz = pbuffer.data(idx_ovl_fg + 47);

    auto ts_xyy_xxyy = pbuffer.data(idx_ovl_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_ovl_fg + 49);

    auto ts_xyy_xxzz = pbuffer.data(idx_ovl_fg + 50);

    auto ts_xyy_xyyy = pbuffer.data(idx_ovl_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_ovl_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_ovl_fg + 53);

    auto ts_xyy_xzzz = pbuffer.data(idx_ovl_fg + 54);

    auto ts_xyy_yyyy = pbuffer.data(idx_ovl_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_ovl_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_ovl_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_ovl_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_ovl_fg + 59);

    auto ts_xyz_xxxx = pbuffer.data(idx_ovl_fg + 60);

    auto ts_xyz_xxxy = pbuffer.data(idx_ovl_fg + 61);

    auto ts_xyz_xxxz = pbuffer.data(idx_ovl_fg + 62);

    auto ts_xyz_xxyy = pbuffer.data(idx_ovl_fg + 63);

    auto ts_xyz_xxyz = pbuffer.data(idx_ovl_fg + 64);

    auto ts_xyz_xxzz = pbuffer.data(idx_ovl_fg + 65);

    auto ts_xyz_xyyy = pbuffer.data(idx_ovl_fg + 66);

    auto ts_xyz_xyyz = pbuffer.data(idx_ovl_fg + 67);

    auto ts_xyz_xyzz = pbuffer.data(idx_ovl_fg + 68);

    auto ts_xyz_xzzz = pbuffer.data(idx_ovl_fg + 69);

    auto ts_xyz_yyyy = pbuffer.data(idx_ovl_fg + 70);

    auto ts_xyz_yyyz = pbuffer.data(idx_ovl_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_ovl_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_ovl_fg + 73);

    auto ts_xyz_zzzz = pbuffer.data(idx_ovl_fg + 74);

    auto ts_xzz_xxxx = pbuffer.data(idx_ovl_fg + 75);

    auto ts_xzz_xxxy = pbuffer.data(idx_ovl_fg + 76);

    auto ts_xzz_xxxz = pbuffer.data(idx_ovl_fg + 77);

    auto ts_xzz_xxyy = pbuffer.data(idx_ovl_fg + 78);

    auto ts_xzz_xxyz = pbuffer.data(idx_ovl_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_ovl_fg + 80);

    auto ts_xzz_xyyy = pbuffer.data(idx_ovl_fg + 81);

    auto ts_xzz_xyyz = pbuffer.data(idx_ovl_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_ovl_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_ovl_fg + 84);

    auto ts_xzz_yyyy = pbuffer.data(idx_ovl_fg + 85);

    auto ts_xzz_yyyz = pbuffer.data(idx_ovl_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_ovl_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_ovl_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_ovl_fg + 89);

    auto ts_yyy_xxxx = pbuffer.data(idx_ovl_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_ovl_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_ovl_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_ovl_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_ovl_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_ovl_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_ovl_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_ovl_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_ovl_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_ovl_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_ovl_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_ovl_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_ovl_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_ovl_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_ovl_fg + 104);

    auto ts_yyz_xxxx = pbuffer.data(idx_ovl_fg + 105);

    auto ts_yyz_xxxy = pbuffer.data(idx_ovl_fg + 106);

    auto ts_yyz_xxxz = pbuffer.data(idx_ovl_fg + 107);

    auto ts_yyz_xxyy = pbuffer.data(idx_ovl_fg + 108);

    auto ts_yyz_xxyz = pbuffer.data(idx_ovl_fg + 109);

    auto ts_yyz_xxzz = pbuffer.data(idx_ovl_fg + 110);

    auto ts_yyz_xyyy = pbuffer.data(idx_ovl_fg + 111);

    auto ts_yyz_xyyz = pbuffer.data(idx_ovl_fg + 112);

    auto ts_yyz_xyzz = pbuffer.data(idx_ovl_fg + 113);

    auto ts_yyz_xzzz = pbuffer.data(idx_ovl_fg + 114);

    auto ts_yyz_yyyy = pbuffer.data(idx_ovl_fg + 115);

    auto ts_yyz_yyyz = pbuffer.data(idx_ovl_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_ovl_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_ovl_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_ovl_fg + 119);

    auto ts_yzz_xxxx = pbuffer.data(idx_ovl_fg + 120);

    auto ts_yzz_xxxy = pbuffer.data(idx_ovl_fg + 121);

    auto ts_yzz_xxxz = pbuffer.data(idx_ovl_fg + 122);

    auto ts_yzz_xxyy = pbuffer.data(idx_ovl_fg + 123);

    auto ts_yzz_xxyz = pbuffer.data(idx_ovl_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_ovl_fg + 125);

    auto ts_yzz_xyyy = pbuffer.data(idx_ovl_fg + 126);

    auto ts_yzz_xyyz = pbuffer.data(idx_ovl_fg + 127);

    auto ts_yzz_xyzz = pbuffer.data(idx_ovl_fg + 128);

    auto ts_yzz_xzzz = pbuffer.data(idx_ovl_fg + 129);

    auto ts_yzz_yyyy = pbuffer.data(idx_ovl_fg + 130);

    auto ts_yzz_yyyz = pbuffer.data(idx_ovl_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_ovl_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_ovl_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_ovl_fg + 134);

    auto ts_zzz_xxxx = pbuffer.data(idx_ovl_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_ovl_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_ovl_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_ovl_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_ovl_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_ovl_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_ovl_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_ovl_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_ovl_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_ovl_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_ovl_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_ovl_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_ovl_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_ovl_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_ovl_fg + 149);

    // Set up 0-15 components of targeted buffer : FG

    auto tk_xxx_xxxx = pbuffer.data(idx_kin_fg);

    auto tk_xxx_xxxy = pbuffer.data(idx_kin_fg + 1);

    auto tk_xxx_xxxz = pbuffer.data(idx_kin_fg + 2);

    auto tk_xxx_xxyy = pbuffer.data(idx_kin_fg + 3);

    auto tk_xxx_xxyz = pbuffer.data(idx_kin_fg + 4);

    auto tk_xxx_xxzz = pbuffer.data(idx_kin_fg + 5);

    auto tk_xxx_xyyy = pbuffer.data(idx_kin_fg + 6);

    auto tk_xxx_xyyz = pbuffer.data(idx_kin_fg + 7);

    auto tk_xxx_xyzz = pbuffer.data(idx_kin_fg + 8);

    auto tk_xxx_xzzz = pbuffer.data(idx_kin_fg + 9);

    auto tk_xxx_yyyy = pbuffer.data(idx_kin_fg + 10);

    auto tk_xxx_yyyz = pbuffer.data(idx_kin_fg + 11);

    auto tk_xxx_yyzz = pbuffer.data(idx_kin_fg + 12);

    auto tk_xxx_yzzz = pbuffer.data(idx_kin_fg + 13);

    auto tk_xxx_zzzz = pbuffer.data(idx_kin_fg + 14);

#pragma omp simd aligned(pa_x,            \
                             tk_x_xxxx,   \
                             tk_x_xxxy,   \
                             tk_x_xxxz,   \
                             tk_x_xxyy,   \
                             tk_x_xxyz,   \
                             tk_x_xxzz,   \
                             tk_x_xyyy,   \
                             tk_x_xyyz,   \
                             tk_x_xyzz,   \
                             tk_x_xzzz,   \
                             tk_x_yyyy,   \
                             tk_x_yyyz,   \
                             tk_x_yyzz,   \
                             tk_x_yzzz,   \
                             tk_x_zzzz,   \
                             tk_xx_xxx,   \
                             tk_xx_xxxx,  \
                             tk_xx_xxxy,  \
                             tk_xx_xxxz,  \
                             tk_xx_xxy,   \
                             tk_xx_xxyy,  \
                             tk_xx_xxyz,  \
                             tk_xx_xxz,   \
                             tk_xx_xxzz,  \
                             tk_xx_xyy,   \
                             tk_xx_xyyy,  \
                             tk_xx_xyyz,  \
                             tk_xx_xyz,   \
                             tk_xx_xyzz,  \
                             tk_xx_xzz,   \
                             tk_xx_xzzz,  \
                             tk_xx_yyy,   \
                             tk_xx_yyyy,  \
                             tk_xx_yyyz,  \
                             tk_xx_yyz,   \
                             tk_xx_yyzz,  \
                             tk_xx_yzz,   \
                             tk_xx_yzzz,  \
                             tk_xx_zzz,   \
                             tk_xx_zzzz,  \
                             tk_xxx_xxxx, \
                             tk_xxx_xxxy, \
                             tk_xxx_xxxz, \
                             tk_xxx_xxyy, \
                             tk_xxx_xxyz, \
                             tk_xxx_xxzz, \
                             tk_xxx_xyyy, \
                             tk_xxx_xyyz, \
                             tk_xxx_xyzz, \
                             tk_xxx_xzzz, \
                             tk_xxx_yyyy, \
                             tk_xxx_yyyz, \
                             tk_xxx_yyzz, \
                             tk_xxx_yzzz, \
                             tk_xxx_zzzz, \
                             ts_x_xxxx,   \
                             ts_x_xxxy,   \
                             ts_x_xxxz,   \
                             ts_x_xxyy,   \
                             ts_x_xxyz,   \
                             ts_x_xxzz,   \
                             ts_x_xyyy,   \
                             ts_x_xyyz,   \
                             ts_x_xyzz,   \
                             ts_x_xzzz,   \
                             ts_x_yyyy,   \
                             ts_x_yyyz,   \
                             ts_x_yyzz,   \
                             ts_x_yzzz,   \
                             ts_x_zzzz,   \
                             ts_xxx_xxxx, \
                             ts_xxx_xxxy, \
                             ts_xxx_xxxz, \
                             ts_xxx_xxyy, \
                             ts_xxx_xxyz, \
                             ts_xxx_xxzz, \
                             ts_xxx_xyyy, \
                             ts_xxx_xyyz, \
                             ts_xxx_xyzz, \
                             ts_xxx_xzzz, \
                             ts_xxx_yyyy, \
                             ts_xxx_yyyz, \
                             ts_xxx_yyzz, \
                             ts_xxx_yzzz, \
                             ts_xxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_xxxx[i] = -4.0 * ts_x_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxx[i] * fe_0 + 4.0 * tk_xx_xxx[i] * fe_0 + tk_xx_xxxx[i] * pa_x[i] +
                         2.0 * ts_xxx_xxxx[i] * fz_0;

        tk_xxx_xxxy[i] = -4.0 * ts_x_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxy[i] * fe_0 + 3.0 * tk_xx_xxy[i] * fe_0 + tk_xx_xxxy[i] * pa_x[i] +
                         2.0 * ts_xxx_xxxy[i] * fz_0;

        tk_xxx_xxxz[i] = -4.0 * ts_x_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxz[i] * fe_0 + 3.0 * tk_xx_xxz[i] * fe_0 + tk_xx_xxxz[i] * pa_x[i] +
                         2.0 * ts_xxx_xxxz[i] * fz_0;

        tk_xxx_xxyy[i] = -4.0 * ts_x_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyy[i] * fe_0 + 2.0 * tk_xx_xyy[i] * fe_0 + tk_xx_xxyy[i] * pa_x[i] +
                         2.0 * ts_xxx_xxyy[i] * fz_0;

        tk_xxx_xxyz[i] = -4.0 * ts_x_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyz[i] * fe_0 + 2.0 * tk_xx_xyz[i] * fe_0 + tk_xx_xxyz[i] * pa_x[i] +
                         2.0 * ts_xxx_xxyz[i] * fz_0;

        tk_xxx_xxzz[i] = -4.0 * ts_x_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxzz[i] * fe_0 + 2.0 * tk_xx_xzz[i] * fe_0 + tk_xx_xxzz[i] * pa_x[i] +
                         2.0 * ts_xxx_xxzz[i] * fz_0;

        tk_xxx_xyyy[i] = -4.0 * ts_x_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyy[i] * fe_0 + tk_xx_yyy[i] * fe_0 + tk_xx_xyyy[i] * pa_x[i] +
                         2.0 * ts_xxx_xyyy[i] * fz_0;

        tk_xxx_xyyz[i] = -4.0 * ts_x_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyz[i] * fe_0 + tk_xx_yyz[i] * fe_0 + tk_xx_xyyz[i] * pa_x[i] +
                         2.0 * ts_xxx_xyyz[i] * fz_0;

        tk_xxx_xyzz[i] = -4.0 * ts_x_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyzz[i] * fe_0 + tk_xx_yzz[i] * fe_0 + tk_xx_xyzz[i] * pa_x[i] +
                         2.0 * ts_xxx_xyzz[i] * fz_0;

        tk_xxx_xzzz[i] = -4.0 * ts_x_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xzzz[i] * fe_0 + tk_xx_zzz[i] * fe_0 + tk_xx_xzzz[i] * pa_x[i] +
                         2.0 * ts_xxx_xzzz[i] * fz_0;

        tk_xxx_yyyy[i] = -4.0 * ts_x_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyy[i] * fe_0 + tk_xx_yyyy[i] * pa_x[i] + 2.0 * ts_xxx_yyyy[i] * fz_0;

        tk_xxx_yyyz[i] = -4.0 * ts_x_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyz[i] * fe_0 + tk_xx_yyyz[i] * pa_x[i] + 2.0 * ts_xxx_yyyz[i] * fz_0;

        tk_xxx_yyzz[i] = -4.0 * ts_x_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyzz[i] * fe_0 + tk_xx_yyzz[i] * pa_x[i] + 2.0 * ts_xxx_yyzz[i] * fz_0;

        tk_xxx_yzzz[i] = -4.0 * ts_x_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yzzz[i] * fe_0 + tk_xx_yzzz[i] * pa_x[i] + 2.0 * ts_xxx_yzzz[i] * fz_0;

        tk_xxx_zzzz[i] = -4.0 * ts_x_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_zzzz[i] * fe_0 + tk_xx_zzzz[i] * pa_x[i] + 2.0 * ts_xxx_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : FG

    auto tk_xxy_xxxx = pbuffer.data(idx_kin_fg + 15);

    auto tk_xxy_xxxy = pbuffer.data(idx_kin_fg + 16);

    auto tk_xxy_xxxz = pbuffer.data(idx_kin_fg + 17);

    auto tk_xxy_xxyy = pbuffer.data(idx_kin_fg + 18);

    auto tk_xxy_xxyz = pbuffer.data(idx_kin_fg + 19);

    auto tk_xxy_xxzz = pbuffer.data(idx_kin_fg + 20);

    auto tk_xxy_xyyy = pbuffer.data(idx_kin_fg + 21);

    auto tk_xxy_xyyz = pbuffer.data(idx_kin_fg + 22);

    auto tk_xxy_xyzz = pbuffer.data(idx_kin_fg + 23);

    auto tk_xxy_xzzz = pbuffer.data(idx_kin_fg + 24);

    auto tk_xxy_yyyy = pbuffer.data(idx_kin_fg + 25);

    auto tk_xxy_yyyz = pbuffer.data(idx_kin_fg + 26);

    auto tk_xxy_yyzz = pbuffer.data(idx_kin_fg + 27);

    auto tk_xxy_yzzz = pbuffer.data(idx_kin_fg + 28);

    auto tk_xxy_zzzz = pbuffer.data(idx_kin_fg + 29);

#pragma omp simd aligned(pa_y,            \
                             tk_xx_xxx,   \
                             tk_xx_xxxx,  \
                             tk_xx_xxxy,  \
                             tk_xx_xxxz,  \
                             tk_xx_xxy,   \
                             tk_xx_xxyy,  \
                             tk_xx_xxyz,  \
                             tk_xx_xxz,   \
                             tk_xx_xxzz,  \
                             tk_xx_xyy,   \
                             tk_xx_xyyy,  \
                             tk_xx_xyyz,  \
                             tk_xx_xyz,   \
                             tk_xx_xyzz,  \
                             tk_xx_xzz,   \
                             tk_xx_xzzz,  \
                             tk_xx_yyy,   \
                             tk_xx_yyyy,  \
                             tk_xx_yyyz,  \
                             tk_xx_yyz,   \
                             tk_xx_yyzz,  \
                             tk_xx_yzz,   \
                             tk_xx_yzzz,  \
                             tk_xx_zzz,   \
                             tk_xx_zzzz,  \
                             tk_xxy_xxxx, \
                             tk_xxy_xxxy, \
                             tk_xxy_xxxz, \
                             tk_xxy_xxyy, \
                             tk_xxy_xxyz, \
                             tk_xxy_xxzz, \
                             tk_xxy_xyyy, \
                             tk_xxy_xyyz, \
                             tk_xxy_xyzz, \
                             tk_xxy_xzzz, \
                             tk_xxy_yyyy, \
                             tk_xxy_yyyz, \
                             tk_xxy_yyzz, \
                             tk_xxy_yzzz, \
                             tk_xxy_zzzz, \
                             ts_xxy_xxxx, \
                             ts_xxy_xxxy, \
                             ts_xxy_xxxz, \
                             ts_xxy_xxyy, \
                             ts_xxy_xxyz, \
                             ts_xxy_xxzz, \
                             ts_xxy_xyyy, \
                             ts_xxy_xyyz, \
                             ts_xxy_xyzz, \
                             ts_xxy_xzzz, \
                             ts_xxy_yyyy, \
                             ts_xxy_yyyz, \
                             ts_xxy_yyzz, \
                             ts_xxy_yzzz, \
                             ts_xxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_xxxx[i] = tk_xx_xxxx[i] * pa_y[i] + 2.0 * ts_xxy_xxxx[i] * fz_0;

        tk_xxy_xxxy[i] = tk_xx_xxx[i] * fe_0 + tk_xx_xxxy[i] * pa_y[i] + 2.0 * ts_xxy_xxxy[i] * fz_0;

        tk_xxy_xxxz[i] = tk_xx_xxxz[i] * pa_y[i] + 2.0 * ts_xxy_xxxz[i] * fz_0;

        tk_xxy_xxyy[i] = 2.0 * tk_xx_xxy[i] * fe_0 + tk_xx_xxyy[i] * pa_y[i] + 2.0 * ts_xxy_xxyy[i] * fz_0;

        tk_xxy_xxyz[i] = tk_xx_xxz[i] * fe_0 + tk_xx_xxyz[i] * pa_y[i] + 2.0 * ts_xxy_xxyz[i] * fz_0;

        tk_xxy_xxzz[i] = tk_xx_xxzz[i] * pa_y[i] + 2.0 * ts_xxy_xxzz[i] * fz_0;

        tk_xxy_xyyy[i] = 3.0 * tk_xx_xyy[i] * fe_0 + tk_xx_xyyy[i] * pa_y[i] + 2.0 * ts_xxy_xyyy[i] * fz_0;

        tk_xxy_xyyz[i] = 2.0 * tk_xx_xyz[i] * fe_0 + tk_xx_xyyz[i] * pa_y[i] + 2.0 * ts_xxy_xyyz[i] * fz_0;

        tk_xxy_xyzz[i] = tk_xx_xzz[i] * fe_0 + tk_xx_xyzz[i] * pa_y[i] + 2.0 * ts_xxy_xyzz[i] * fz_0;

        tk_xxy_xzzz[i] = tk_xx_xzzz[i] * pa_y[i] + 2.0 * ts_xxy_xzzz[i] * fz_0;

        tk_xxy_yyyy[i] = 4.0 * tk_xx_yyy[i] * fe_0 + tk_xx_yyyy[i] * pa_y[i] + 2.0 * ts_xxy_yyyy[i] * fz_0;

        tk_xxy_yyyz[i] = 3.0 * tk_xx_yyz[i] * fe_0 + tk_xx_yyyz[i] * pa_y[i] + 2.0 * ts_xxy_yyyz[i] * fz_0;

        tk_xxy_yyzz[i] = 2.0 * tk_xx_yzz[i] * fe_0 + tk_xx_yyzz[i] * pa_y[i] + 2.0 * ts_xxy_yyzz[i] * fz_0;

        tk_xxy_yzzz[i] = tk_xx_zzz[i] * fe_0 + tk_xx_yzzz[i] * pa_y[i] + 2.0 * ts_xxy_yzzz[i] * fz_0;

        tk_xxy_zzzz[i] = tk_xx_zzzz[i] * pa_y[i] + 2.0 * ts_xxy_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : FG

    auto tk_xxz_xxxx = pbuffer.data(idx_kin_fg + 30);

    auto tk_xxz_xxxy = pbuffer.data(idx_kin_fg + 31);

    auto tk_xxz_xxxz = pbuffer.data(idx_kin_fg + 32);

    auto tk_xxz_xxyy = pbuffer.data(idx_kin_fg + 33);

    auto tk_xxz_xxyz = pbuffer.data(idx_kin_fg + 34);

    auto tk_xxz_xxzz = pbuffer.data(idx_kin_fg + 35);

    auto tk_xxz_xyyy = pbuffer.data(idx_kin_fg + 36);

    auto tk_xxz_xyyz = pbuffer.data(idx_kin_fg + 37);

    auto tk_xxz_xyzz = pbuffer.data(idx_kin_fg + 38);

    auto tk_xxz_xzzz = pbuffer.data(idx_kin_fg + 39);

    auto tk_xxz_yyyy = pbuffer.data(idx_kin_fg + 40);

    auto tk_xxz_yyyz = pbuffer.data(idx_kin_fg + 41);

    auto tk_xxz_yyzz = pbuffer.data(idx_kin_fg + 42);

    auto tk_xxz_yzzz = pbuffer.data(idx_kin_fg + 43);

    auto tk_xxz_zzzz = pbuffer.data(idx_kin_fg + 44);

#pragma omp simd aligned(pa_z,            \
                             tk_xx_xxx,   \
                             tk_xx_xxxx,  \
                             tk_xx_xxxy,  \
                             tk_xx_xxxz,  \
                             tk_xx_xxy,   \
                             tk_xx_xxyy,  \
                             tk_xx_xxyz,  \
                             tk_xx_xxz,   \
                             tk_xx_xxzz,  \
                             tk_xx_xyy,   \
                             tk_xx_xyyy,  \
                             tk_xx_xyyz,  \
                             tk_xx_xyz,   \
                             tk_xx_xyzz,  \
                             tk_xx_xzz,   \
                             tk_xx_xzzz,  \
                             tk_xx_yyy,   \
                             tk_xx_yyyy,  \
                             tk_xx_yyyz,  \
                             tk_xx_yyz,   \
                             tk_xx_yyzz,  \
                             tk_xx_yzz,   \
                             tk_xx_yzzz,  \
                             tk_xx_zzz,   \
                             tk_xx_zzzz,  \
                             tk_xxz_xxxx, \
                             tk_xxz_xxxy, \
                             tk_xxz_xxxz, \
                             tk_xxz_xxyy, \
                             tk_xxz_xxyz, \
                             tk_xxz_xxzz, \
                             tk_xxz_xyyy, \
                             tk_xxz_xyyz, \
                             tk_xxz_xyzz, \
                             tk_xxz_xzzz, \
                             tk_xxz_yyyy, \
                             tk_xxz_yyyz, \
                             tk_xxz_yyzz, \
                             tk_xxz_yzzz, \
                             tk_xxz_zzzz, \
                             ts_xxz_xxxx, \
                             ts_xxz_xxxy, \
                             ts_xxz_xxxz, \
                             ts_xxz_xxyy, \
                             ts_xxz_xxyz, \
                             ts_xxz_xxzz, \
                             ts_xxz_xyyy, \
                             ts_xxz_xyyz, \
                             ts_xxz_xyzz, \
                             ts_xxz_xzzz, \
                             ts_xxz_yyyy, \
                             ts_xxz_yyyz, \
                             ts_xxz_yyzz, \
                             ts_xxz_yzzz, \
                             ts_xxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_xxxx[i] = tk_xx_xxxx[i] * pa_z[i] + 2.0 * ts_xxz_xxxx[i] * fz_0;

        tk_xxz_xxxy[i] = tk_xx_xxxy[i] * pa_z[i] + 2.0 * ts_xxz_xxxy[i] * fz_0;

        tk_xxz_xxxz[i] = tk_xx_xxx[i] * fe_0 + tk_xx_xxxz[i] * pa_z[i] + 2.0 * ts_xxz_xxxz[i] * fz_0;

        tk_xxz_xxyy[i] = tk_xx_xxyy[i] * pa_z[i] + 2.0 * ts_xxz_xxyy[i] * fz_0;

        tk_xxz_xxyz[i] = tk_xx_xxy[i] * fe_0 + tk_xx_xxyz[i] * pa_z[i] + 2.0 * ts_xxz_xxyz[i] * fz_0;

        tk_xxz_xxzz[i] = 2.0 * tk_xx_xxz[i] * fe_0 + tk_xx_xxzz[i] * pa_z[i] + 2.0 * ts_xxz_xxzz[i] * fz_0;

        tk_xxz_xyyy[i] = tk_xx_xyyy[i] * pa_z[i] + 2.0 * ts_xxz_xyyy[i] * fz_0;

        tk_xxz_xyyz[i] = tk_xx_xyy[i] * fe_0 + tk_xx_xyyz[i] * pa_z[i] + 2.0 * ts_xxz_xyyz[i] * fz_0;

        tk_xxz_xyzz[i] = 2.0 * tk_xx_xyz[i] * fe_0 + tk_xx_xyzz[i] * pa_z[i] + 2.0 * ts_xxz_xyzz[i] * fz_0;

        tk_xxz_xzzz[i] = 3.0 * tk_xx_xzz[i] * fe_0 + tk_xx_xzzz[i] * pa_z[i] + 2.0 * ts_xxz_xzzz[i] * fz_0;

        tk_xxz_yyyy[i] = tk_xx_yyyy[i] * pa_z[i] + 2.0 * ts_xxz_yyyy[i] * fz_0;

        tk_xxz_yyyz[i] = tk_xx_yyy[i] * fe_0 + tk_xx_yyyz[i] * pa_z[i] + 2.0 * ts_xxz_yyyz[i] * fz_0;

        tk_xxz_yyzz[i] = 2.0 * tk_xx_yyz[i] * fe_0 + tk_xx_yyzz[i] * pa_z[i] + 2.0 * ts_xxz_yyzz[i] * fz_0;

        tk_xxz_yzzz[i] = 3.0 * tk_xx_yzz[i] * fe_0 + tk_xx_yzzz[i] * pa_z[i] + 2.0 * ts_xxz_yzzz[i] * fz_0;

        tk_xxz_zzzz[i] = 4.0 * tk_xx_zzz[i] * fe_0 + tk_xx_zzzz[i] * pa_z[i] + 2.0 * ts_xxz_zzzz[i] * fz_0;
    }

    // Set up 45-60 components of targeted buffer : FG

    auto tk_xyy_xxxx = pbuffer.data(idx_kin_fg + 45);

    auto tk_xyy_xxxy = pbuffer.data(idx_kin_fg + 46);

    auto tk_xyy_xxxz = pbuffer.data(idx_kin_fg + 47);

    auto tk_xyy_xxyy = pbuffer.data(idx_kin_fg + 48);

    auto tk_xyy_xxyz = pbuffer.data(idx_kin_fg + 49);

    auto tk_xyy_xxzz = pbuffer.data(idx_kin_fg + 50);

    auto tk_xyy_xyyy = pbuffer.data(idx_kin_fg + 51);

    auto tk_xyy_xyyz = pbuffer.data(idx_kin_fg + 52);

    auto tk_xyy_xyzz = pbuffer.data(idx_kin_fg + 53);

    auto tk_xyy_xzzz = pbuffer.data(idx_kin_fg + 54);

    auto tk_xyy_yyyy = pbuffer.data(idx_kin_fg + 55);

    auto tk_xyy_yyyz = pbuffer.data(idx_kin_fg + 56);

    auto tk_xyy_yyzz = pbuffer.data(idx_kin_fg + 57);

    auto tk_xyy_yzzz = pbuffer.data(idx_kin_fg + 58);

    auto tk_xyy_zzzz = pbuffer.data(idx_kin_fg + 59);

#pragma omp simd aligned(pa_x,            \
                             tk_xyy_xxxx, \
                             tk_xyy_xxxy, \
                             tk_xyy_xxxz, \
                             tk_xyy_xxyy, \
                             tk_xyy_xxyz, \
                             tk_xyy_xxzz, \
                             tk_xyy_xyyy, \
                             tk_xyy_xyyz, \
                             tk_xyy_xyzz, \
                             tk_xyy_xzzz, \
                             tk_xyy_yyyy, \
                             tk_xyy_yyyz, \
                             tk_xyy_yyzz, \
                             tk_xyy_yzzz, \
                             tk_xyy_zzzz, \
                             tk_yy_xxx,   \
                             tk_yy_xxxx,  \
                             tk_yy_xxxy,  \
                             tk_yy_xxxz,  \
                             tk_yy_xxy,   \
                             tk_yy_xxyy,  \
                             tk_yy_xxyz,  \
                             tk_yy_xxz,   \
                             tk_yy_xxzz,  \
                             tk_yy_xyy,   \
                             tk_yy_xyyy,  \
                             tk_yy_xyyz,  \
                             tk_yy_xyz,   \
                             tk_yy_xyzz,  \
                             tk_yy_xzz,   \
                             tk_yy_xzzz,  \
                             tk_yy_yyy,   \
                             tk_yy_yyyy,  \
                             tk_yy_yyyz,  \
                             tk_yy_yyz,   \
                             tk_yy_yyzz,  \
                             tk_yy_yzz,   \
                             tk_yy_yzzz,  \
                             tk_yy_zzz,   \
                             tk_yy_zzzz,  \
                             ts_xyy_xxxx, \
                             ts_xyy_xxxy, \
                             ts_xyy_xxxz, \
                             ts_xyy_xxyy, \
                             ts_xyy_xxyz, \
                             ts_xyy_xxzz, \
                             ts_xyy_xyyy, \
                             ts_xyy_xyyz, \
                             ts_xyy_xyzz, \
                             ts_xyy_xzzz, \
                             ts_xyy_yyyy, \
                             ts_xyy_yyyz, \
                             ts_xyy_yyzz, \
                             ts_xyy_yzzz, \
                             ts_xyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_xxxx[i] = 4.0 * tk_yy_xxx[i] * fe_0 + tk_yy_xxxx[i] * pa_x[i] + 2.0 * ts_xyy_xxxx[i] * fz_0;

        tk_xyy_xxxy[i] = 3.0 * tk_yy_xxy[i] * fe_0 + tk_yy_xxxy[i] * pa_x[i] + 2.0 * ts_xyy_xxxy[i] * fz_0;

        tk_xyy_xxxz[i] = 3.0 * tk_yy_xxz[i] * fe_0 + tk_yy_xxxz[i] * pa_x[i] + 2.0 * ts_xyy_xxxz[i] * fz_0;

        tk_xyy_xxyy[i] = 2.0 * tk_yy_xyy[i] * fe_0 + tk_yy_xxyy[i] * pa_x[i] + 2.0 * ts_xyy_xxyy[i] * fz_0;

        tk_xyy_xxyz[i] = 2.0 * tk_yy_xyz[i] * fe_0 + tk_yy_xxyz[i] * pa_x[i] + 2.0 * ts_xyy_xxyz[i] * fz_0;

        tk_xyy_xxzz[i] = 2.0 * tk_yy_xzz[i] * fe_0 + tk_yy_xxzz[i] * pa_x[i] + 2.0 * ts_xyy_xxzz[i] * fz_0;

        tk_xyy_xyyy[i] = tk_yy_yyy[i] * fe_0 + tk_yy_xyyy[i] * pa_x[i] + 2.0 * ts_xyy_xyyy[i] * fz_0;

        tk_xyy_xyyz[i] = tk_yy_yyz[i] * fe_0 + tk_yy_xyyz[i] * pa_x[i] + 2.0 * ts_xyy_xyyz[i] * fz_0;

        tk_xyy_xyzz[i] = tk_yy_yzz[i] * fe_0 + tk_yy_xyzz[i] * pa_x[i] + 2.0 * ts_xyy_xyzz[i] * fz_0;

        tk_xyy_xzzz[i] = tk_yy_zzz[i] * fe_0 + tk_yy_xzzz[i] * pa_x[i] + 2.0 * ts_xyy_xzzz[i] * fz_0;

        tk_xyy_yyyy[i] = tk_yy_yyyy[i] * pa_x[i] + 2.0 * ts_xyy_yyyy[i] * fz_0;

        tk_xyy_yyyz[i] = tk_yy_yyyz[i] * pa_x[i] + 2.0 * ts_xyy_yyyz[i] * fz_0;

        tk_xyy_yyzz[i] = tk_yy_yyzz[i] * pa_x[i] + 2.0 * ts_xyy_yyzz[i] * fz_0;

        tk_xyy_yzzz[i] = tk_yy_yzzz[i] * pa_x[i] + 2.0 * ts_xyy_yzzz[i] * fz_0;

        tk_xyy_zzzz[i] = tk_yy_zzzz[i] * pa_x[i] + 2.0 * ts_xyy_zzzz[i] * fz_0;
    }

    // Set up 60-75 components of targeted buffer : FG

    auto tk_xyz_xxxx = pbuffer.data(idx_kin_fg + 60);

    auto tk_xyz_xxxy = pbuffer.data(idx_kin_fg + 61);

    auto tk_xyz_xxxz = pbuffer.data(idx_kin_fg + 62);

    auto tk_xyz_xxyy = pbuffer.data(idx_kin_fg + 63);

    auto tk_xyz_xxyz = pbuffer.data(idx_kin_fg + 64);

    auto tk_xyz_xxzz = pbuffer.data(idx_kin_fg + 65);

    auto tk_xyz_xyyy = pbuffer.data(idx_kin_fg + 66);

    auto tk_xyz_xyyz = pbuffer.data(idx_kin_fg + 67);

    auto tk_xyz_xyzz = pbuffer.data(idx_kin_fg + 68);

    auto tk_xyz_xzzz = pbuffer.data(idx_kin_fg + 69);

    auto tk_xyz_yyyy = pbuffer.data(idx_kin_fg + 70);

    auto tk_xyz_yyyz = pbuffer.data(idx_kin_fg + 71);

    auto tk_xyz_yyzz = pbuffer.data(idx_kin_fg + 72);

    auto tk_xyz_yzzz = pbuffer.data(idx_kin_fg + 73);

    auto tk_xyz_zzzz = pbuffer.data(idx_kin_fg + 74);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             tk_xy_xxxy,  \
                             tk_xy_xxyy,  \
                             tk_xy_xyyy,  \
                             tk_xyz_xxxx, \
                             tk_xyz_xxxy, \
                             tk_xyz_xxxz, \
                             tk_xyz_xxyy, \
                             tk_xyz_xxyz, \
                             tk_xyz_xxzz, \
                             tk_xyz_xyyy, \
                             tk_xyz_xyyz, \
                             tk_xyz_xyzz, \
                             tk_xyz_xzzz, \
                             tk_xyz_yyyy, \
                             tk_xyz_yyyz, \
                             tk_xyz_yyzz, \
                             tk_xyz_yzzz, \
                             tk_xyz_zzzz, \
                             tk_xz_xxxx,  \
                             tk_xz_xxxz,  \
                             tk_xz_xxzz,  \
                             tk_xz_xzzz,  \
                             tk_yz_xxyz,  \
                             tk_yz_xyyz,  \
                             tk_yz_xyz,   \
                             tk_yz_xyzz,  \
                             tk_yz_yyyy,  \
                             tk_yz_yyyz,  \
                             tk_yz_yyz,   \
                             tk_yz_yyzz,  \
                             tk_yz_yzz,   \
                             tk_yz_yzzz,  \
                             tk_yz_zzzz,  \
                             ts_xyz_xxxx, \
                             ts_xyz_xxxy, \
                             ts_xyz_xxxz, \
                             ts_xyz_xxyy, \
                             ts_xyz_xxyz, \
                             ts_xyz_xxzz, \
                             ts_xyz_xyyy, \
                             ts_xyz_xyyz, \
                             ts_xyz_xyzz, \
                             ts_xyz_xzzz, \
                             ts_xyz_yyyy, \
                             ts_xyz_yyyz, \
                             ts_xyz_yyzz, \
                             ts_xyz_yzzz, \
                             ts_xyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyz_xxxx[i] = tk_xz_xxxx[i] * pa_y[i] + 2.0 * ts_xyz_xxxx[i] * fz_0;

        tk_xyz_xxxy[i] = tk_xy_xxxy[i] * pa_z[i] + 2.0 * ts_xyz_xxxy[i] * fz_0;

        tk_xyz_xxxz[i] = tk_xz_xxxz[i] * pa_y[i] + 2.0 * ts_xyz_xxxz[i] * fz_0;

        tk_xyz_xxyy[i] = tk_xy_xxyy[i] * pa_z[i] + 2.0 * ts_xyz_xxyy[i] * fz_0;

        tk_xyz_xxyz[i] = 2.0 * tk_yz_xyz[i] * fe_0 + tk_yz_xxyz[i] * pa_x[i] + 2.0 * ts_xyz_xxyz[i] * fz_0;

        tk_xyz_xxzz[i] = tk_xz_xxzz[i] * pa_y[i] + 2.0 * ts_xyz_xxzz[i] * fz_0;

        tk_xyz_xyyy[i] = tk_xy_xyyy[i] * pa_z[i] + 2.0 * ts_xyz_xyyy[i] * fz_0;

        tk_xyz_xyyz[i] = tk_yz_yyz[i] * fe_0 + tk_yz_xyyz[i] * pa_x[i] + 2.0 * ts_xyz_xyyz[i] * fz_0;

        tk_xyz_xyzz[i] = tk_yz_yzz[i] * fe_0 + tk_yz_xyzz[i] * pa_x[i] + 2.0 * ts_xyz_xyzz[i] * fz_0;

        tk_xyz_xzzz[i] = tk_xz_xzzz[i] * pa_y[i] + 2.0 * ts_xyz_xzzz[i] * fz_0;

        tk_xyz_yyyy[i] = tk_yz_yyyy[i] * pa_x[i] + 2.0 * ts_xyz_yyyy[i] * fz_0;

        tk_xyz_yyyz[i] = tk_yz_yyyz[i] * pa_x[i] + 2.0 * ts_xyz_yyyz[i] * fz_0;

        tk_xyz_yyzz[i] = tk_yz_yyzz[i] * pa_x[i] + 2.0 * ts_xyz_yyzz[i] * fz_0;

        tk_xyz_yzzz[i] = tk_yz_yzzz[i] * pa_x[i] + 2.0 * ts_xyz_yzzz[i] * fz_0;

        tk_xyz_zzzz[i] = tk_yz_zzzz[i] * pa_x[i] + 2.0 * ts_xyz_zzzz[i] * fz_0;
    }

    // Set up 75-90 components of targeted buffer : FG

    auto tk_xzz_xxxx = pbuffer.data(idx_kin_fg + 75);

    auto tk_xzz_xxxy = pbuffer.data(idx_kin_fg + 76);

    auto tk_xzz_xxxz = pbuffer.data(idx_kin_fg + 77);

    auto tk_xzz_xxyy = pbuffer.data(idx_kin_fg + 78);

    auto tk_xzz_xxyz = pbuffer.data(idx_kin_fg + 79);

    auto tk_xzz_xxzz = pbuffer.data(idx_kin_fg + 80);

    auto tk_xzz_xyyy = pbuffer.data(idx_kin_fg + 81);

    auto tk_xzz_xyyz = pbuffer.data(idx_kin_fg + 82);

    auto tk_xzz_xyzz = pbuffer.data(idx_kin_fg + 83);

    auto tk_xzz_xzzz = pbuffer.data(idx_kin_fg + 84);

    auto tk_xzz_yyyy = pbuffer.data(idx_kin_fg + 85);

    auto tk_xzz_yyyz = pbuffer.data(idx_kin_fg + 86);

    auto tk_xzz_yyzz = pbuffer.data(idx_kin_fg + 87);

    auto tk_xzz_yzzz = pbuffer.data(idx_kin_fg + 88);

    auto tk_xzz_zzzz = pbuffer.data(idx_kin_fg + 89);

#pragma omp simd aligned(pa_x,            \
                             tk_xzz_xxxx, \
                             tk_xzz_xxxy, \
                             tk_xzz_xxxz, \
                             tk_xzz_xxyy, \
                             tk_xzz_xxyz, \
                             tk_xzz_xxzz, \
                             tk_xzz_xyyy, \
                             tk_xzz_xyyz, \
                             tk_xzz_xyzz, \
                             tk_xzz_xzzz, \
                             tk_xzz_yyyy, \
                             tk_xzz_yyyz, \
                             tk_xzz_yyzz, \
                             tk_xzz_yzzz, \
                             tk_xzz_zzzz, \
                             tk_zz_xxx,   \
                             tk_zz_xxxx,  \
                             tk_zz_xxxy,  \
                             tk_zz_xxxz,  \
                             tk_zz_xxy,   \
                             tk_zz_xxyy,  \
                             tk_zz_xxyz,  \
                             tk_zz_xxz,   \
                             tk_zz_xxzz,  \
                             tk_zz_xyy,   \
                             tk_zz_xyyy,  \
                             tk_zz_xyyz,  \
                             tk_zz_xyz,   \
                             tk_zz_xyzz,  \
                             tk_zz_xzz,   \
                             tk_zz_xzzz,  \
                             tk_zz_yyy,   \
                             tk_zz_yyyy,  \
                             tk_zz_yyyz,  \
                             tk_zz_yyz,   \
                             tk_zz_yyzz,  \
                             tk_zz_yzz,   \
                             tk_zz_yzzz,  \
                             tk_zz_zzz,   \
                             tk_zz_zzzz,  \
                             ts_xzz_xxxx, \
                             ts_xzz_xxxy, \
                             ts_xzz_xxxz, \
                             ts_xzz_xxyy, \
                             ts_xzz_xxyz, \
                             ts_xzz_xxzz, \
                             ts_xzz_xyyy, \
                             ts_xzz_xyyz, \
                             ts_xzz_xyzz, \
                             ts_xzz_xzzz, \
                             ts_xzz_yyyy, \
                             ts_xzz_yyyz, \
                             ts_xzz_yyzz, \
                             ts_xzz_yzzz, \
                             ts_xzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_xxxx[i] = 4.0 * tk_zz_xxx[i] * fe_0 + tk_zz_xxxx[i] * pa_x[i] + 2.0 * ts_xzz_xxxx[i] * fz_0;

        tk_xzz_xxxy[i] = 3.0 * tk_zz_xxy[i] * fe_0 + tk_zz_xxxy[i] * pa_x[i] + 2.0 * ts_xzz_xxxy[i] * fz_0;

        tk_xzz_xxxz[i] = 3.0 * tk_zz_xxz[i] * fe_0 + tk_zz_xxxz[i] * pa_x[i] + 2.0 * ts_xzz_xxxz[i] * fz_0;

        tk_xzz_xxyy[i] = 2.0 * tk_zz_xyy[i] * fe_0 + tk_zz_xxyy[i] * pa_x[i] + 2.0 * ts_xzz_xxyy[i] * fz_0;

        tk_xzz_xxyz[i] = 2.0 * tk_zz_xyz[i] * fe_0 + tk_zz_xxyz[i] * pa_x[i] + 2.0 * ts_xzz_xxyz[i] * fz_0;

        tk_xzz_xxzz[i] = 2.0 * tk_zz_xzz[i] * fe_0 + tk_zz_xxzz[i] * pa_x[i] + 2.0 * ts_xzz_xxzz[i] * fz_0;

        tk_xzz_xyyy[i] = tk_zz_yyy[i] * fe_0 + tk_zz_xyyy[i] * pa_x[i] + 2.0 * ts_xzz_xyyy[i] * fz_0;

        tk_xzz_xyyz[i] = tk_zz_yyz[i] * fe_0 + tk_zz_xyyz[i] * pa_x[i] + 2.0 * ts_xzz_xyyz[i] * fz_0;

        tk_xzz_xyzz[i] = tk_zz_yzz[i] * fe_0 + tk_zz_xyzz[i] * pa_x[i] + 2.0 * ts_xzz_xyzz[i] * fz_0;

        tk_xzz_xzzz[i] = tk_zz_zzz[i] * fe_0 + tk_zz_xzzz[i] * pa_x[i] + 2.0 * ts_xzz_xzzz[i] * fz_0;

        tk_xzz_yyyy[i] = tk_zz_yyyy[i] * pa_x[i] + 2.0 * ts_xzz_yyyy[i] * fz_0;

        tk_xzz_yyyz[i] = tk_zz_yyyz[i] * pa_x[i] + 2.0 * ts_xzz_yyyz[i] * fz_0;

        tk_xzz_yyzz[i] = tk_zz_yyzz[i] * pa_x[i] + 2.0 * ts_xzz_yyzz[i] * fz_0;

        tk_xzz_yzzz[i] = tk_zz_yzzz[i] * pa_x[i] + 2.0 * ts_xzz_yzzz[i] * fz_0;

        tk_xzz_zzzz[i] = tk_zz_zzzz[i] * pa_x[i] + 2.0 * ts_xzz_zzzz[i] * fz_0;
    }

    // Set up 90-105 components of targeted buffer : FG

    auto tk_yyy_xxxx = pbuffer.data(idx_kin_fg + 90);

    auto tk_yyy_xxxy = pbuffer.data(idx_kin_fg + 91);

    auto tk_yyy_xxxz = pbuffer.data(idx_kin_fg + 92);

    auto tk_yyy_xxyy = pbuffer.data(idx_kin_fg + 93);

    auto tk_yyy_xxyz = pbuffer.data(idx_kin_fg + 94);

    auto tk_yyy_xxzz = pbuffer.data(idx_kin_fg + 95);

    auto tk_yyy_xyyy = pbuffer.data(idx_kin_fg + 96);

    auto tk_yyy_xyyz = pbuffer.data(idx_kin_fg + 97);

    auto tk_yyy_xyzz = pbuffer.data(idx_kin_fg + 98);

    auto tk_yyy_xzzz = pbuffer.data(idx_kin_fg + 99);

    auto tk_yyy_yyyy = pbuffer.data(idx_kin_fg + 100);

    auto tk_yyy_yyyz = pbuffer.data(idx_kin_fg + 101);

    auto tk_yyy_yyzz = pbuffer.data(idx_kin_fg + 102);

    auto tk_yyy_yzzz = pbuffer.data(idx_kin_fg + 103);

    auto tk_yyy_zzzz = pbuffer.data(idx_kin_fg + 104);

#pragma omp simd aligned(pa_y,            \
                             tk_y_xxxx,   \
                             tk_y_xxxy,   \
                             tk_y_xxxz,   \
                             tk_y_xxyy,   \
                             tk_y_xxyz,   \
                             tk_y_xxzz,   \
                             tk_y_xyyy,   \
                             tk_y_xyyz,   \
                             tk_y_xyzz,   \
                             tk_y_xzzz,   \
                             tk_y_yyyy,   \
                             tk_y_yyyz,   \
                             tk_y_yyzz,   \
                             tk_y_yzzz,   \
                             tk_y_zzzz,   \
                             tk_yy_xxx,   \
                             tk_yy_xxxx,  \
                             tk_yy_xxxy,  \
                             tk_yy_xxxz,  \
                             tk_yy_xxy,   \
                             tk_yy_xxyy,  \
                             tk_yy_xxyz,  \
                             tk_yy_xxz,   \
                             tk_yy_xxzz,  \
                             tk_yy_xyy,   \
                             tk_yy_xyyy,  \
                             tk_yy_xyyz,  \
                             tk_yy_xyz,   \
                             tk_yy_xyzz,  \
                             tk_yy_xzz,   \
                             tk_yy_xzzz,  \
                             tk_yy_yyy,   \
                             tk_yy_yyyy,  \
                             tk_yy_yyyz,  \
                             tk_yy_yyz,   \
                             tk_yy_yyzz,  \
                             tk_yy_yzz,   \
                             tk_yy_yzzz,  \
                             tk_yy_zzz,   \
                             tk_yy_zzzz,  \
                             tk_yyy_xxxx, \
                             tk_yyy_xxxy, \
                             tk_yyy_xxxz, \
                             tk_yyy_xxyy, \
                             tk_yyy_xxyz, \
                             tk_yyy_xxzz, \
                             tk_yyy_xyyy, \
                             tk_yyy_xyyz, \
                             tk_yyy_xyzz, \
                             tk_yyy_xzzz, \
                             tk_yyy_yyyy, \
                             tk_yyy_yyyz, \
                             tk_yyy_yyzz, \
                             tk_yyy_yzzz, \
                             tk_yyy_zzzz, \
                             ts_y_xxxx,   \
                             ts_y_xxxy,   \
                             ts_y_xxxz,   \
                             ts_y_xxyy,   \
                             ts_y_xxyz,   \
                             ts_y_xxzz,   \
                             ts_y_xyyy,   \
                             ts_y_xyyz,   \
                             ts_y_xyzz,   \
                             ts_y_xzzz,   \
                             ts_y_yyyy,   \
                             ts_y_yyyz,   \
                             ts_y_yyzz,   \
                             ts_y_yzzz,   \
                             ts_y_zzzz,   \
                             ts_yyy_xxxx, \
                             ts_yyy_xxxy, \
                             ts_yyy_xxxz, \
                             ts_yyy_xxyy, \
                             ts_yyy_xxyz, \
                             ts_yyy_xxzz, \
                             ts_yyy_xyyy, \
                             ts_yyy_xyyz, \
                             ts_yyy_xyzz, \
                             ts_yyy_xzzz, \
                             ts_yyy_yyyy, \
                             ts_yyy_yyyz, \
                             ts_yyy_yyzz, \
                             ts_yyy_yzzz, \
                             ts_yyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_xxxx[i] = -4.0 * ts_y_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxx[i] * fe_0 + tk_yy_xxxx[i] * pa_y[i] + 2.0 * ts_yyy_xxxx[i] * fz_0;

        tk_yyy_xxxy[i] = -4.0 * ts_y_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxy[i] * fe_0 + tk_yy_xxx[i] * fe_0 + tk_yy_xxxy[i] * pa_y[i] +
                         2.0 * ts_yyy_xxxy[i] * fz_0;

        tk_yyy_xxxz[i] = -4.0 * ts_y_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxz[i] * fe_0 + tk_yy_xxxz[i] * pa_y[i] + 2.0 * ts_yyy_xxxz[i] * fz_0;

        tk_yyy_xxyy[i] = -4.0 * ts_y_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyy[i] * fe_0 + 2.0 * tk_yy_xxy[i] * fe_0 + tk_yy_xxyy[i] * pa_y[i] +
                         2.0 * ts_yyy_xxyy[i] * fz_0;

        tk_yyy_xxyz[i] = -4.0 * ts_y_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyz[i] * fe_0 + tk_yy_xxz[i] * fe_0 + tk_yy_xxyz[i] * pa_y[i] +
                         2.0 * ts_yyy_xxyz[i] * fz_0;

        tk_yyy_xxzz[i] = -4.0 * ts_y_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxzz[i] * fe_0 + tk_yy_xxzz[i] * pa_y[i] + 2.0 * ts_yyy_xxzz[i] * fz_0;

        tk_yyy_xyyy[i] = -4.0 * ts_y_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyy[i] * fe_0 + 3.0 * tk_yy_xyy[i] * fe_0 + tk_yy_xyyy[i] * pa_y[i] +
                         2.0 * ts_yyy_xyyy[i] * fz_0;

        tk_yyy_xyyz[i] = -4.0 * ts_y_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyz[i] * fe_0 + 2.0 * tk_yy_xyz[i] * fe_0 + tk_yy_xyyz[i] * pa_y[i] +
                         2.0 * ts_yyy_xyyz[i] * fz_0;

        tk_yyy_xyzz[i] = -4.0 * ts_y_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyzz[i] * fe_0 + tk_yy_xzz[i] * fe_0 + tk_yy_xyzz[i] * pa_y[i] +
                         2.0 * ts_yyy_xyzz[i] * fz_0;

        tk_yyy_xzzz[i] = -4.0 * ts_y_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xzzz[i] * fe_0 + tk_yy_xzzz[i] * pa_y[i] + 2.0 * ts_yyy_xzzz[i] * fz_0;

        tk_yyy_yyyy[i] = -4.0 * ts_y_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyy[i] * fe_0 + 4.0 * tk_yy_yyy[i] * fe_0 + tk_yy_yyyy[i] * pa_y[i] +
                         2.0 * ts_yyy_yyyy[i] * fz_0;

        tk_yyy_yyyz[i] = -4.0 * ts_y_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyz[i] * fe_0 + 3.0 * tk_yy_yyz[i] * fe_0 + tk_yy_yyyz[i] * pa_y[i] +
                         2.0 * ts_yyy_yyyz[i] * fz_0;

        tk_yyy_yyzz[i] = -4.0 * ts_y_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyzz[i] * fe_0 + 2.0 * tk_yy_yzz[i] * fe_0 + tk_yy_yyzz[i] * pa_y[i] +
                         2.0 * ts_yyy_yyzz[i] * fz_0;

        tk_yyy_yzzz[i] = -4.0 * ts_y_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yzzz[i] * fe_0 + tk_yy_zzz[i] * fe_0 + tk_yy_yzzz[i] * pa_y[i] +
                         2.0 * ts_yyy_yzzz[i] * fz_0;

        tk_yyy_zzzz[i] = -4.0 * ts_y_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_zzzz[i] * fe_0 + tk_yy_zzzz[i] * pa_y[i] + 2.0 * ts_yyy_zzzz[i] * fz_0;
    }

    // Set up 105-120 components of targeted buffer : FG

    auto tk_yyz_xxxx = pbuffer.data(idx_kin_fg + 105);

    auto tk_yyz_xxxy = pbuffer.data(idx_kin_fg + 106);

    auto tk_yyz_xxxz = pbuffer.data(idx_kin_fg + 107);

    auto tk_yyz_xxyy = pbuffer.data(idx_kin_fg + 108);

    auto tk_yyz_xxyz = pbuffer.data(idx_kin_fg + 109);

    auto tk_yyz_xxzz = pbuffer.data(idx_kin_fg + 110);

    auto tk_yyz_xyyy = pbuffer.data(idx_kin_fg + 111);

    auto tk_yyz_xyyz = pbuffer.data(idx_kin_fg + 112);

    auto tk_yyz_xyzz = pbuffer.data(idx_kin_fg + 113);

    auto tk_yyz_xzzz = pbuffer.data(idx_kin_fg + 114);

    auto tk_yyz_yyyy = pbuffer.data(idx_kin_fg + 115);

    auto tk_yyz_yyyz = pbuffer.data(idx_kin_fg + 116);

    auto tk_yyz_yyzz = pbuffer.data(idx_kin_fg + 117);

    auto tk_yyz_yzzz = pbuffer.data(idx_kin_fg + 118);

    auto tk_yyz_zzzz = pbuffer.data(idx_kin_fg + 119);

#pragma omp simd aligned(pa_z,            \
                             tk_yy_xxx,   \
                             tk_yy_xxxx,  \
                             tk_yy_xxxy,  \
                             tk_yy_xxxz,  \
                             tk_yy_xxy,   \
                             tk_yy_xxyy,  \
                             tk_yy_xxyz,  \
                             tk_yy_xxz,   \
                             tk_yy_xxzz,  \
                             tk_yy_xyy,   \
                             tk_yy_xyyy,  \
                             tk_yy_xyyz,  \
                             tk_yy_xyz,   \
                             tk_yy_xyzz,  \
                             tk_yy_xzz,   \
                             tk_yy_xzzz,  \
                             tk_yy_yyy,   \
                             tk_yy_yyyy,  \
                             tk_yy_yyyz,  \
                             tk_yy_yyz,   \
                             tk_yy_yyzz,  \
                             tk_yy_yzz,   \
                             tk_yy_yzzz,  \
                             tk_yy_zzz,   \
                             tk_yy_zzzz,  \
                             tk_yyz_xxxx, \
                             tk_yyz_xxxy, \
                             tk_yyz_xxxz, \
                             tk_yyz_xxyy, \
                             tk_yyz_xxyz, \
                             tk_yyz_xxzz, \
                             tk_yyz_xyyy, \
                             tk_yyz_xyyz, \
                             tk_yyz_xyzz, \
                             tk_yyz_xzzz, \
                             tk_yyz_yyyy, \
                             tk_yyz_yyyz, \
                             tk_yyz_yyzz, \
                             tk_yyz_yzzz, \
                             tk_yyz_zzzz, \
                             ts_yyz_xxxx, \
                             ts_yyz_xxxy, \
                             ts_yyz_xxxz, \
                             ts_yyz_xxyy, \
                             ts_yyz_xxyz, \
                             ts_yyz_xxzz, \
                             ts_yyz_xyyy, \
                             ts_yyz_xyyz, \
                             ts_yyz_xyzz, \
                             ts_yyz_xzzz, \
                             ts_yyz_yyyy, \
                             ts_yyz_yyyz, \
                             ts_yyz_yyzz, \
                             ts_yyz_yzzz, \
                             ts_yyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_xxxx[i] = tk_yy_xxxx[i] * pa_z[i] + 2.0 * ts_yyz_xxxx[i] * fz_0;

        tk_yyz_xxxy[i] = tk_yy_xxxy[i] * pa_z[i] + 2.0 * ts_yyz_xxxy[i] * fz_0;

        tk_yyz_xxxz[i] = tk_yy_xxx[i] * fe_0 + tk_yy_xxxz[i] * pa_z[i] + 2.0 * ts_yyz_xxxz[i] * fz_0;

        tk_yyz_xxyy[i] = tk_yy_xxyy[i] * pa_z[i] + 2.0 * ts_yyz_xxyy[i] * fz_0;

        tk_yyz_xxyz[i] = tk_yy_xxy[i] * fe_0 + tk_yy_xxyz[i] * pa_z[i] + 2.0 * ts_yyz_xxyz[i] * fz_0;

        tk_yyz_xxzz[i] = 2.0 * tk_yy_xxz[i] * fe_0 + tk_yy_xxzz[i] * pa_z[i] + 2.0 * ts_yyz_xxzz[i] * fz_0;

        tk_yyz_xyyy[i] = tk_yy_xyyy[i] * pa_z[i] + 2.0 * ts_yyz_xyyy[i] * fz_0;

        tk_yyz_xyyz[i] = tk_yy_xyy[i] * fe_0 + tk_yy_xyyz[i] * pa_z[i] + 2.0 * ts_yyz_xyyz[i] * fz_0;

        tk_yyz_xyzz[i] = 2.0 * tk_yy_xyz[i] * fe_0 + tk_yy_xyzz[i] * pa_z[i] + 2.0 * ts_yyz_xyzz[i] * fz_0;

        tk_yyz_xzzz[i] = 3.0 * tk_yy_xzz[i] * fe_0 + tk_yy_xzzz[i] * pa_z[i] + 2.0 * ts_yyz_xzzz[i] * fz_0;

        tk_yyz_yyyy[i] = tk_yy_yyyy[i] * pa_z[i] + 2.0 * ts_yyz_yyyy[i] * fz_0;

        tk_yyz_yyyz[i] = tk_yy_yyy[i] * fe_0 + tk_yy_yyyz[i] * pa_z[i] + 2.0 * ts_yyz_yyyz[i] * fz_0;

        tk_yyz_yyzz[i] = 2.0 * tk_yy_yyz[i] * fe_0 + tk_yy_yyzz[i] * pa_z[i] + 2.0 * ts_yyz_yyzz[i] * fz_0;

        tk_yyz_yzzz[i] = 3.0 * tk_yy_yzz[i] * fe_0 + tk_yy_yzzz[i] * pa_z[i] + 2.0 * ts_yyz_yzzz[i] * fz_0;

        tk_yyz_zzzz[i] = 4.0 * tk_yy_zzz[i] * fe_0 + tk_yy_zzzz[i] * pa_z[i] + 2.0 * ts_yyz_zzzz[i] * fz_0;
    }

    // Set up 120-135 components of targeted buffer : FG

    auto tk_yzz_xxxx = pbuffer.data(idx_kin_fg + 120);

    auto tk_yzz_xxxy = pbuffer.data(idx_kin_fg + 121);

    auto tk_yzz_xxxz = pbuffer.data(idx_kin_fg + 122);

    auto tk_yzz_xxyy = pbuffer.data(idx_kin_fg + 123);

    auto tk_yzz_xxyz = pbuffer.data(idx_kin_fg + 124);

    auto tk_yzz_xxzz = pbuffer.data(idx_kin_fg + 125);

    auto tk_yzz_xyyy = pbuffer.data(idx_kin_fg + 126);

    auto tk_yzz_xyyz = pbuffer.data(idx_kin_fg + 127);

    auto tk_yzz_xyzz = pbuffer.data(idx_kin_fg + 128);

    auto tk_yzz_xzzz = pbuffer.data(idx_kin_fg + 129);

    auto tk_yzz_yyyy = pbuffer.data(idx_kin_fg + 130);

    auto tk_yzz_yyyz = pbuffer.data(idx_kin_fg + 131);

    auto tk_yzz_yyzz = pbuffer.data(idx_kin_fg + 132);

    auto tk_yzz_yzzz = pbuffer.data(idx_kin_fg + 133);

    auto tk_yzz_zzzz = pbuffer.data(idx_kin_fg + 134);

#pragma omp simd aligned(pa_y,            \
                             tk_yzz_xxxx, \
                             tk_yzz_xxxy, \
                             tk_yzz_xxxz, \
                             tk_yzz_xxyy, \
                             tk_yzz_xxyz, \
                             tk_yzz_xxzz, \
                             tk_yzz_xyyy, \
                             tk_yzz_xyyz, \
                             tk_yzz_xyzz, \
                             tk_yzz_xzzz, \
                             tk_yzz_yyyy, \
                             tk_yzz_yyyz, \
                             tk_yzz_yyzz, \
                             tk_yzz_yzzz, \
                             tk_yzz_zzzz, \
                             tk_zz_xxx,   \
                             tk_zz_xxxx,  \
                             tk_zz_xxxy,  \
                             tk_zz_xxxz,  \
                             tk_zz_xxy,   \
                             tk_zz_xxyy,  \
                             tk_zz_xxyz,  \
                             tk_zz_xxz,   \
                             tk_zz_xxzz,  \
                             tk_zz_xyy,   \
                             tk_zz_xyyy,  \
                             tk_zz_xyyz,  \
                             tk_zz_xyz,   \
                             tk_zz_xyzz,  \
                             tk_zz_xzz,   \
                             tk_zz_xzzz,  \
                             tk_zz_yyy,   \
                             tk_zz_yyyy,  \
                             tk_zz_yyyz,  \
                             tk_zz_yyz,   \
                             tk_zz_yyzz,  \
                             tk_zz_yzz,   \
                             tk_zz_yzzz,  \
                             tk_zz_zzz,   \
                             tk_zz_zzzz,  \
                             ts_yzz_xxxx, \
                             ts_yzz_xxxy, \
                             ts_yzz_xxxz, \
                             ts_yzz_xxyy, \
                             ts_yzz_xxyz, \
                             ts_yzz_xxzz, \
                             ts_yzz_xyyy, \
                             ts_yzz_xyyz, \
                             ts_yzz_xyzz, \
                             ts_yzz_xzzz, \
                             ts_yzz_yyyy, \
                             ts_yzz_yyyz, \
                             ts_yzz_yyzz, \
                             ts_yzz_yzzz, \
                             ts_yzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_xxxx[i] = tk_zz_xxxx[i] * pa_y[i] + 2.0 * ts_yzz_xxxx[i] * fz_0;

        tk_yzz_xxxy[i] = tk_zz_xxx[i] * fe_0 + tk_zz_xxxy[i] * pa_y[i] + 2.0 * ts_yzz_xxxy[i] * fz_0;

        tk_yzz_xxxz[i] = tk_zz_xxxz[i] * pa_y[i] + 2.0 * ts_yzz_xxxz[i] * fz_0;

        tk_yzz_xxyy[i] = 2.0 * tk_zz_xxy[i] * fe_0 + tk_zz_xxyy[i] * pa_y[i] + 2.0 * ts_yzz_xxyy[i] * fz_0;

        tk_yzz_xxyz[i] = tk_zz_xxz[i] * fe_0 + tk_zz_xxyz[i] * pa_y[i] + 2.0 * ts_yzz_xxyz[i] * fz_0;

        tk_yzz_xxzz[i] = tk_zz_xxzz[i] * pa_y[i] + 2.0 * ts_yzz_xxzz[i] * fz_0;

        tk_yzz_xyyy[i] = 3.0 * tk_zz_xyy[i] * fe_0 + tk_zz_xyyy[i] * pa_y[i] + 2.0 * ts_yzz_xyyy[i] * fz_0;

        tk_yzz_xyyz[i] = 2.0 * tk_zz_xyz[i] * fe_0 + tk_zz_xyyz[i] * pa_y[i] + 2.0 * ts_yzz_xyyz[i] * fz_0;

        tk_yzz_xyzz[i] = tk_zz_xzz[i] * fe_0 + tk_zz_xyzz[i] * pa_y[i] + 2.0 * ts_yzz_xyzz[i] * fz_0;

        tk_yzz_xzzz[i] = tk_zz_xzzz[i] * pa_y[i] + 2.0 * ts_yzz_xzzz[i] * fz_0;

        tk_yzz_yyyy[i] = 4.0 * tk_zz_yyy[i] * fe_0 + tk_zz_yyyy[i] * pa_y[i] + 2.0 * ts_yzz_yyyy[i] * fz_0;

        tk_yzz_yyyz[i] = 3.0 * tk_zz_yyz[i] * fe_0 + tk_zz_yyyz[i] * pa_y[i] + 2.0 * ts_yzz_yyyz[i] * fz_0;

        tk_yzz_yyzz[i] = 2.0 * tk_zz_yzz[i] * fe_0 + tk_zz_yyzz[i] * pa_y[i] + 2.0 * ts_yzz_yyzz[i] * fz_0;

        tk_yzz_yzzz[i] = tk_zz_zzz[i] * fe_0 + tk_zz_yzzz[i] * pa_y[i] + 2.0 * ts_yzz_yzzz[i] * fz_0;

        tk_yzz_zzzz[i] = tk_zz_zzzz[i] * pa_y[i] + 2.0 * ts_yzz_zzzz[i] * fz_0;
    }

    // Set up 135-150 components of targeted buffer : FG

    auto tk_zzz_xxxx = pbuffer.data(idx_kin_fg + 135);

    auto tk_zzz_xxxy = pbuffer.data(idx_kin_fg + 136);

    auto tk_zzz_xxxz = pbuffer.data(idx_kin_fg + 137);

    auto tk_zzz_xxyy = pbuffer.data(idx_kin_fg + 138);

    auto tk_zzz_xxyz = pbuffer.data(idx_kin_fg + 139);

    auto tk_zzz_xxzz = pbuffer.data(idx_kin_fg + 140);

    auto tk_zzz_xyyy = pbuffer.data(idx_kin_fg + 141);

    auto tk_zzz_xyyz = pbuffer.data(idx_kin_fg + 142);

    auto tk_zzz_xyzz = pbuffer.data(idx_kin_fg + 143);

    auto tk_zzz_xzzz = pbuffer.data(idx_kin_fg + 144);

    auto tk_zzz_yyyy = pbuffer.data(idx_kin_fg + 145);

    auto tk_zzz_yyyz = pbuffer.data(idx_kin_fg + 146);

    auto tk_zzz_yyzz = pbuffer.data(idx_kin_fg + 147);

    auto tk_zzz_yzzz = pbuffer.data(idx_kin_fg + 148);

    auto tk_zzz_zzzz = pbuffer.data(idx_kin_fg + 149);

#pragma omp simd aligned(pa_z,            \
                             tk_z_xxxx,   \
                             tk_z_xxxy,   \
                             tk_z_xxxz,   \
                             tk_z_xxyy,   \
                             tk_z_xxyz,   \
                             tk_z_xxzz,   \
                             tk_z_xyyy,   \
                             tk_z_xyyz,   \
                             tk_z_xyzz,   \
                             tk_z_xzzz,   \
                             tk_z_yyyy,   \
                             tk_z_yyyz,   \
                             tk_z_yyzz,   \
                             tk_z_yzzz,   \
                             tk_z_zzzz,   \
                             tk_zz_xxx,   \
                             tk_zz_xxxx,  \
                             tk_zz_xxxy,  \
                             tk_zz_xxxz,  \
                             tk_zz_xxy,   \
                             tk_zz_xxyy,  \
                             tk_zz_xxyz,  \
                             tk_zz_xxz,   \
                             tk_zz_xxzz,  \
                             tk_zz_xyy,   \
                             tk_zz_xyyy,  \
                             tk_zz_xyyz,  \
                             tk_zz_xyz,   \
                             tk_zz_xyzz,  \
                             tk_zz_xzz,   \
                             tk_zz_xzzz,  \
                             tk_zz_yyy,   \
                             tk_zz_yyyy,  \
                             tk_zz_yyyz,  \
                             tk_zz_yyz,   \
                             tk_zz_yyzz,  \
                             tk_zz_yzz,   \
                             tk_zz_yzzz,  \
                             tk_zz_zzz,   \
                             tk_zz_zzzz,  \
                             tk_zzz_xxxx, \
                             tk_zzz_xxxy, \
                             tk_zzz_xxxz, \
                             tk_zzz_xxyy, \
                             tk_zzz_xxyz, \
                             tk_zzz_xxzz, \
                             tk_zzz_xyyy, \
                             tk_zzz_xyyz, \
                             tk_zzz_xyzz, \
                             tk_zzz_xzzz, \
                             tk_zzz_yyyy, \
                             tk_zzz_yyyz, \
                             tk_zzz_yyzz, \
                             tk_zzz_yzzz, \
                             tk_zzz_zzzz, \
                             ts_z_xxxx,   \
                             ts_z_xxxy,   \
                             ts_z_xxxz,   \
                             ts_z_xxyy,   \
                             ts_z_xxyz,   \
                             ts_z_xxzz,   \
                             ts_z_xyyy,   \
                             ts_z_xyyz,   \
                             ts_z_xyzz,   \
                             ts_z_xzzz,   \
                             ts_z_yyyy,   \
                             ts_z_yyyz,   \
                             ts_z_yyzz,   \
                             ts_z_yzzz,   \
                             ts_z_zzzz,   \
                             ts_zzz_xxxx, \
                             ts_zzz_xxxy, \
                             ts_zzz_xxxz, \
                             ts_zzz_xxyy, \
                             ts_zzz_xxyz, \
                             ts_zzz_xxzz, \
                             ts_zzz_xyyy, \
                             ts_zzz_xyyz, \
                             ts_zzz_xyzz, \
                             ts_zzz_xzzz, \
                             ts_zzz_yyyy, \
                             ts_zzz_yyyz, \
                             ts_zzz_yyzz, \
                             ts_zzz_yzzz, \
                             ts_zzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_xxxx[i] = -4.0 * ts_z_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxx[i] * fe_0 + tk_zz_xxxx[i] * pa_z[i] + 2.0 * ts_zzz_xxxx[i] * fz_0;

        tk_zzz_xxxy[i] = -4.0 * ts_z_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxy[i] * fe_0 + tk_zz_xxxy[i] * pa_z[i] + 2.0 * ts_zzz_xxxy[i] * fz_0;

        tk_zzz_xxxz[i] = -4.0 * ts_z_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxz[i] * fe_0 + tk_zz_xxx[i] * fe_0 + tk_zz_xxxz[i] * pa_z[i] +
                         2.0 * ts_zzz_xxxz[i] * fz_0;

        tk_zzz_xxyy[i] = -4.0 * ts_z_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyy[i] * fe_0 + tk_zz_xxyy[i] * pa_z[i] + 2.0 * ts_zzz_xxyy[i] * fz_0;

        tk_zzz_xxyz[i] = -4.0 * ts_z_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyz[i] * fe_0 + tk_zz_xxy[i] * fe_0 + tk_zz_xxyz[i] * pa_z[i] +
                         2.0 * ts_zzz_xxyz[i] * fz_0;

        tk_zzz_xxzz[i] = -4.0 * ts_z_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxzz[i] * fe_0 + 2.0 * tk_zz_xxz[i] * fe_0 + tk_zz_xxzz[i] * pa_z[i] +
                         2.0 * ts_zzz_xxzz[i] * fz_0;

        tk_zzz_xyyy[i] = -4.0 * ts_z_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyy[i] * fe_0 + tk_zz_xyyy[i] * pa_z[i] + 2.0 * ts_zzz_xyyy[i] * fz_0;

        tk_zzz_xyyz[i] = -4.0 * ts_z_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyz[i] * fe_0 + tk_zz_xyy[i] * fe_0 + tk_zz_xyyz[i] * pa_z[i] +
                         2.0 * ts_zzz_xyyz[i] * fz_0;

        tk_zzz_xyzz[i] = -4.0 * ts_z_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyzz[i] * fe_0 + 2.0 * tk_zz_xyz[i] * fe_0 + tk_zz_xyzz[i] * pa_z[i] +
                         2.0 * ts_zzz_xyzz[i] * fz_0;

        tk_zzz_xzzz[i] = -4.0 * ts_z_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xzzz[i] * fe_0 + 3.0 * tk_zz_xzz[i] * fe_0 + tk_zz_xzzz[i] * pa_z[i] +
                         2.0 * ts_zzz_xzzz[i] * fz_0;

        tk_zzz_yyyy[i] = -4.0 * ts_z_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyy[i] * fe_0 + tk_zz_yyyy[i] * pa_z[i] + 2.0 * ts_zzz_yyyy[i] * fz_0;

        tk_zzz_yyyz[i] = -4.0 * ts_z_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyz[i] * fe_0 + tk_zz_yyy[i] * fe_0 + tk_zz_yyyz[i] * pa_z[i] +
                         2.0 * ts_zzz_yyyz[i] * fz_0;

        tk_zzz_yyzz[i] = -4.0 * ts_z_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyzz[i] * fe_0 + 2.0 * tk_zz_yyz[i] * fe_0 + tk_zz_yyzz[i] * pa_z[i] +
                         2.0 * ts_zzz_yyzz[i] * fz_0;

        tk_zzz_yzzz[i] = -4.0 * ts_z_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yzzz[i] * fe_0 + 3.0 * tk_zz_yzz[i] * fe_0 + tk_zz_yzzz[i] * pa_z[i] +
                         2.0 * ts_zzz_yzzz[i] * fz_0;

        tk_zzz_zzzz[i] = -4.0 * ts_z_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_zzzz[i] * fe_0 + 4.0 * tk_zz_zzz[i] * fe_0 + tk_zz_zzzz[i] * pa_z[i] +
                         2.0 * ts_zzz_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
