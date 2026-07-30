#include "LocalCorePotentialPrimRecFG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fg,
                                  const size_t idx_pg,
                                  const size_t idx_df,
                                  const size_t idx_dg,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx = pbuffer.data(idx_pg);

    auto tg_x_xxxy = pbuffer.data(idx_pg + 1);

    auto tg_x_xxxz = pbuffer.data(idx_pg + 2);

    auto tg_x_xxyy = pbuffer.data(idx_pg + 3);

    auto tg_x_xxyz = pbuffer.data(idx_pg + 4);

    auto tg_x_xxzz = pbuffer.data(idx_pg + 5);

    auto tg_x_xyyy = pbuffer.data(idx_pg + 6);

    auto tg_x_xyyz = pbuffer.data(idx_pg + 7);

    auto tg_x_xyzz = pbuffer.data(idx_pg + 8);

    auto tg_x_xzzz = pbuffer.data(idx_pg + 9);

    auto tg_x_yyyy = pbuffer.data(idx_pg + 10);

    auto tg_x_yyyz = pbuffer.data(idx_pg + 11);

    auto tg_x_yyzz = pbuffer.data(idx_pg + 12);

    auto tg_x_yzzz = pbuffer.data(idx_pg + 13);

    auto tg_x_zzzz = pbuffer.data(idx_pg + 14);

    auto tg_y_xxxx = pbuffer.data(idx_pg + 15);

    auto tg_y_xxxy = pbuffer.data(idx_pg + 16);

    auto tg_y_xxxz = pbuffer.data(idx_pg + 17);

    auto tg_y_xxyy = pbuffer.data(idx_pg + 18);

    auto tg_y_xxyz = pbuffer.data(idx_pg + 19);

    auto tg_y_xxzz = pbuffer.data(idx_pg + 20);

    auto tg_y_xyyy = pbuffer.data(idx_pg + 21);

    auto tg_y_xyyz = pbuffer.data(idx_pg + 22);

    auto tg_y_xyzz = pbuffer.data(idx_pg + 23);

    auto tg_y_xzzz = pbuffer.data(idx_pg + 24);

    auto tg_y_yyyy = pbuffer.data(idx_pg + 25);

    auto tg_y_yyyz = pbuffer.data(idx_pg + 26);

    auto tg_y_yyzz = pbuffer.data(idx_pg + 27);

    auto tg_y_yzzz = pbuffer.data(idx_pg + 28);

    auto tg_y_zzzz = pbuffer.data(idx_pg + 29);

    auto tg_z_xxxx = pbuffer.data(idx_pg + 30);

    auto tg_z_xxxy = pbuffer.data(idx_pg + 31);

    auto tg_z_xxxz = pbuffer.data(idx_pg + 32);

    auto tg_z_xxyy = pbuffer.data(idx_pg + 33);

    auto tg_z_xxyz = pbuffer.data(idx_pg + 34);

    auto tg_z_xxzz = pbuffer.data(idx_pg + 35);

    auto tg_z_xyyy = pbuffer.data(idx_pg + 36);

    auto tg_z_xyyz = pbuffer.data(idx_pg + 37);

    auto tg_z_xyzz = pbuffer.data(idx_pg + 38);

    auto tg_z_xzzz = pbuffer.data(idx_pg + 39);

    auto tg_z_yyyy = pbuffer.data(idx_pg + 40);

    auto tg_z_yyyz = pbuffer.data(idx_pg + 41);

    auto tg_z_yyzz = pbuffer.data(idx_pg + 42);

    auto tg_z_yzzz = pbuffer.data(idx_pg + 43);

    auto tg_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx = pbuffer.data(idx_df);

    auto tg_xx_xxy = pbuffer.data(idx_df + 1);

    auto tg_xx_xxz = pbuffer.data(idx_df + 2);

    auto tg_xx_xyy = pbuffer.data(idx_df + 3);

    auto tg_xx_xyz = pbuffer.data(idx_df + 4);

    auto tg_xx_xzz = pbuffer.data(idx_df + 5);

    auto tg_xx_yyy = pbuffer.data(idx_df + 6);

    auto tg_xx_yyz = pbuffer.data(idx_df + 7);

    auto tg_xx_yzz = pbuffer.data(idx_df + 8);

    auto tg_xx_zzz = pbuffer.data(idx_df + 9);

    auto tg_yy_xxx = pbuffer.data(idx_df + 30);

    auto tg_yy_xxy = pbuffer.data(idx_df + 31);

    auto tg_yy_xxz = pbuffer.data(idx_df + 32);

    auto tg_yy_xyy = pbuffer.data(idx_df + 33);

    auto tg_yy_xyz = pbuffer.data(idx_df + 34);

    auto tg_yy_xzz = pbuffer.data(idx_df + 35);

    auto tg_yy_yyy = pbuffer.data(idx_df + 36);

    auto tg_yy_yyz = pbuffer.data(idx_df + 37);

    auto tg_yy_yzz = pbuffer.data(idx_df + 38);

    auto tg_yy_zzz = pbuffer.data(idx_df + 39);

    auto tg_yz_xyz = pbuffer.data(idx_df + 44);

    auto tg_yz_yyz = pbuffer.data(idx_df + 47);

    auto tg_yz_yzz = pbuffer.data(idx_df + 48);

    auto tg_zz_xxx = pbuffer.data(idx_df + 50);

    auto tg_zz_xxy = pbuffer.data(idx_df + 51);

    auto tg_zz_xxz = pbuffer.data(idx_df + 52);

    auto tg_zz_xyy = pbuffer.data(idx_df + 53);

    auto tg_zz_xyz = pbuffer.data(idx_df + 54);

    auto tg_zz_xzz = pbuffer.data(idx_df + 55);

    auto tg_zz_yyy = pbuffer.data(idx_df + 56);

    auto tg_zz_yyz = pbuffer.data(idx_df + 57);

    auto tg_zz_yzz = pbuffer.data(idx_df + 58);

    auto tg_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx = pbuffer.data(idx_dg);

    auto tg_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto tg_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto tg_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto tg_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto tg_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto tg_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto tg_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto tg_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto tg_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto tg_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto tg_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto tg_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto tg_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto tg_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto tg_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto tg_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto tg_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto tg_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto tg_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto tg_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto tg_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto tg_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto tg_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto tg_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto tg_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto tg_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto tg_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto tg_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto tg_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto tg_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto tg_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto tg_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto tg_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto tg_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto tg_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto tg_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto tg_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto tg_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto tg_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto tg_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto tg_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto tg_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto tg_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto tg_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto tg_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto tg_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto tg_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto tg_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto tg_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto tg_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto tg_yz_yyyy = pbuffer.data(idx_dg + 70);

    auto tg_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto tg_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto tg_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto tg_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto tg_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto tg_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto tg_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto tg_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto tg_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto tg_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto tg_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto tg_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto tg_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto tg_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto tg_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto tg_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto tg_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto tg_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto tg_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up components of targeted buffer : FG

    auto tg_xxx_xxxx = pbuffer.data(idx_fg);

    auto tg_xxx_xxxy = pbuffer.data(idx_fg + 1);

    auto tg_xxx_xxxz = pbuffer.data(idx_fg + 2);

    auto tg_xxx_xxyy = pbuffer.data(idx_fg + 3);

    auto tg_xxx_xxyz = pbuffer.data(idx_fg + 4);

    auto tg_xxx_xxzz = pbuffer.data(idx_fg + 5);

    auto tg_xxx_xyyy = pbuffer.data(idx_fg + 6);

    auto tg_xxx_xyyz = pbuffer.data(idx_fg + 7);

    auto tg_xxx_xyzz = pbuffer.data(idx_fg + 8);

    auto tg_xxx_xzzz = pbuffer.data(idx_fg + 9);

    auto tg_xxx_yyyy = pbuffer.data(idx_fg + 10);

    auto tg_xxx_yyyz = pbuffer.data(idx_fg + 11);

    auto tg_xxx_yyzz = pbuffer.data(idx_fg + 12);

    auto tg_xxx_yzzz = pbuffer.data(idx_fg + 13);

    auto tg_xxx_zzzz = pbuffer.data(idx_fg + 14);

    auto tg_xxy_xxxx = pbuffer.data(idx_fg + 15);

    auto tg_xxy_xxxy = pbuffer.data(idx_fg + 16);

    auto tg_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto tg_xxy_xxyy = pbuffer.data(idx_fg + 18);

    auto tg_xxy_xxyz = pbuffer.data(idx_fg + 19);

    auto tg_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto tg_xxy_xyyy = pbuffer.data(idx_fg + 21);

    auto tg_xxy_xyyz = pbuffer.data(idx_fg + 22);

    auto tg_xxy_xyzz = pbuffer.data(idx_fg + 23);

    auto tg_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto tg_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto tg_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto tg_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto tg_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto tg_xxy_zzzz = pbuffer.data(idx_fg + 29);

    auto tg_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto tg_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto tg_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto tg_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto tg_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto tg_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto tg_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto tg_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto tg_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto tg_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto tg_xxz_yyyy = pbuffer.data(idx_fg + 40);

    auto tg_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto tg_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto tg_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto tg_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto tg_xyy_xxxx = pbuffer.data(idx_fg + 45);

    auto tg_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto tg_xyy_xxxz = pbuffer.data(idx_fg + 47);

    auto tg_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto tg_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto tg_xyy_xxzz = pbuffer.data(idx_fg + 50);

    auto tg_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto tg_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto tg_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto tg_xyy_xzzz = pbuffer.data(idx_fg + 54);

    auto tg_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto tg_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto tg_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto tg_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto tg_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto tg_xyz_xxxx = pbuffer.data(idx_fg + 60);

    auto tg_xyz_xxxy = pbuffer.data(idx_fg + 61);

    auto tg_xyz_xxxz = pbuffer.data(idx_fg + 62);

    auto tg_xyz_xxyy = pbuffer.data(idx_fg + 63);

    auto tg_xyz_xxyz = pbuffer.data(idx_fg + 64);

    auto tg_xyz_xxzz = pbuffer.data(idx_fg + 65);

    auto tg_xyz_xyyy = pbuffer.data(idx_fg + 66);

    auto tg_xyz_xyyz = pbuffer.data(idx_fg + 67);

    auto tg_xyz_xyzz = pbuffer.data(idx_fg + 68);

    auto tg_xyz_xzzz = pbuffer.data(idx_fg + 69);

    auto tg_xyz_yyyy = pbuffer.data(idx_fg + 70);

    auto tg_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto tg_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto tg_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto tg_xyz_zzzz = pbuffer.data(idx_fg + 74);

    auto tg_xzz_xxxx = pbuffer.data(idx_fg + 75);

    auto tg_xzz_xxxy = pbuffer.data(idx_fg + 76);

    auto tg_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto tg_xzz_xxyy = pbuffer.data(idx_fg + 78);

    auto tg_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto tg_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto tg_xzz_xyyy = pbuffer.data(idx_fg + 81);

    auto tg_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto tg_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto tg_xzz_xzzz = pbuffer.data(idx_fg + 84);

    auto tg_xzz_yyyy = pbuffer.data(idx_fg + 85);

    auto tg_xzz_yyyz = pbuffer.data(idx_fg + 86);

    auto tg_xzz_yyzz = pbuffer.data(idx_fg + 87);

    auto tg_xzz_yzzz = pbuffer.data(idx_fg + 88);

    auto tg_xzz_zzzz = pbuffer.data(idx_fg + 89);

    auto tg_yyy_xxxx = pbuffer.data(idx_fg + 90);

    auto tg_yyy_xxxy = pbuffer.data(idx_fg + 91);

    auto tg_yyy_xxxz = pbuffer.data(idx_fg + 92);

    auto tg_yyy_xxyy = pbuffer.data(idx_fg + 93);

    auto tg_yyy_xxyz = pbuffer.data(idx_fg + 94);

    auto tg_yyy_xxzz = pbuffer.data(idx_fg + 95);

    auto tg_yyy_xyyy = pbuffer.data(idx_fg + 96);

    auto tg_yyy_xyyz = pbuffer.data(idx_fg + 97);

    auto tg_yyy_xyzz = pbuffer.data(idx_fg + 98);

    auto tg_yyy_xzzz = pbuffer.data(idx_fg + 99);

    auto tg_yyy_yyyy = pbuffer.data(idx_fg + 100);

    auto tg_yyy_yyyz = pbuffer.data(idx_fg + 101);

    auto tg_yyy_yyzz = pbuffer.data(idx_fg + 102);

    auto tg_yyy_yzzz = pbuffer.data(idx_fg + 103);

    auto tg_yyy_zzzz = pbuffer.data(idx_fg + 104);

    auto tg_yyz_xxxx = pbuffer.data(idx_fg + 105);

    auto tg_yyz_xxxy = pbuffer.data(idx_fg + 106);

    auto tg_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto tg_yyz_xxyy = pbuffer.data(idx_fg + 108);

    auto tg_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto tg_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto tg_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto tg_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto tg_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto tg_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto tg_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto tg_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto tg_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto tg_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto tg_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto tg_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto tg_yzz_xxxy = pbuffer.data(idx_fg + 121);

    auto tg_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto tg_yzz_xxyy = pbuffer.data(idx_fg + 123);

    auto tg_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto tg_yzz_xxzz = pbuffer.data(idx_fg + 125);

    auto tg_yzz_xyyy = pbuffer.data(idx_fg + 126);

    auto tg_yzz_xyyz = pbuffer.data(idx_fg + 127);

    auto tg_yzz_xyzz = pbuffer.data(idx_fg + 128);

    auto tg_yzz_xzzz = pbuffer.data(idx_fg + 129);

    auto tg_yzz_yyyy = pbuffer.data(idx_fg + 130);

    auto tg_yzz_yyyz = pbuffer.data(idx_fg + 131);

    auto tg_yzz_yyzz = pbuffer.data(idx_fg + 132);

    auto tg_yzz_yzzz = pbuffer.data(idx_fg + 133);

    auto tg_yzz_zzzz = pbuffer.data(idx_fg + 134);

    auto tg_zzz_xxxx = pbuffer.data(idx_fg + 135);

    auto tg_zzz_xxxy = pbuffer.data(idx_fg + 136);

    auto tg_zzz_xxxz = pbuffer.data(idx_fg + 137);

    auto tg_zzz_xxyy = pbuffer.data(idx_fg + 138);

    auto tg_zzz_xxyz = pbuffer.data(idx_fg + 139);

    auto tg_zzz_xxzz = pbuffer.data(idx_fg + 140);

    auto tg_zzz_xyyy = pbuffer.data(idx_fg + 141);

    auto tg_zzz_xyyz = pbuffer.data(idx_fg + 142);

    auto tg_zzz_xyzz = pbuffer.data(idx_fg + 143);

    auto tg_zzz_xzzz = pbuffer.data(idx_fg + 144);

    auto tg_zzz_yyyy = pbuffer.data(idx_fg + 145);

    auto tg_zzz_yyyz = pbuffer.data(idx_fg + 146);

    auto tg_zzz_yyzz = pbuffer.data(idx_fg + 147);

    auto tg_zzz_yzzz = pbuffer.data(idx_fg + 148);

    auto tg_zzz_zzzz = pbuffer.data(idx_fg + 149);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_xxxx, tg_x_xxxy, tg_x_xxxz, tg_x_xxyy, tg_x_xxyz, tg_x_xxzz, tg_x_xyyy, tg_x_xyyz, tg_x_xyzz, tg_x_xzzz, tg_x_yyyy, tg_x_yyyz, tg_x_yyzz, tg_x_yzzz, tg_x_zzzz, tg_xx_xxx, tg_xx_xxxx, tg_xx_xxxy, tg_xx_xxxz, tg_xx_xxy, tg_xx_xxyy, tg_xx_xxyz, tg_xx_xxz, tg_xx_xxzz, tg_xx_xyy, tg_xx_xyyy, tg_xx_xyyz, tg_xx_xyz, tg_xx_xyzz, tg_xx_xzz, tg_xx_xzzz, tg_xx_yyy, tg_xx_yyyy, tg_xx_yyyz, tg_xx_yyz, tg_xx_yyzz, tg_xx_yzz, tg_xx_yzzz, tg_xx_zzz, tg_xx_zzzz, tg_xxx_xxxx, tg_xxx_xxxy, tg_xxx_xxxz, tg_xxx_xxyy, tg_xxx_xxyz, tg_xxx_xxzz, tg_xxx_xyyy, tg_xxx_xyyz, tg_xxx_xyzz, tg_xxx_xzzz, tg_xxx_yyyy, tg_xxx_yyyz, tg_xxx_yyzz, tg_xxx_yzzz, tg_xxx_zzzz, tg_xxy_xxxx, tg_xxy_xxxy, tg_xxy_xxxz, tg_xxy_xxyy, tg_xxy_xxyz, tg_xxy_xxzz, tg_xxy_xyyy, tg_xxy_xyyz, tg_xxy_xyzz, tg_xxy_xzzz, tg_xxy_yyyy, tg_xxy_yyyz, tg_xxy_yyzz, tg_xxy_yzzz, tg_xxy_zzzz, tg_xxz_xxxx, tg_xxz_xxxy, tg_xxz_xxxz, tg_xxz_xxyy, tg_xxz_xxyz, tg_xxz_xxzz, tg_xxz_xyyy, tg_xxz_xyyz, tg_xxz_xyzz, tg_xxz_xzzz, tg_xxz_yyyy, tg_xxz_yyyz, tg_xxz_yyzz, tg_xxz_yzzz, tg_xxz_zzzz, tg_xy_xxxy, tg_xy_xxyy, tg_xy_xyyy, tg_xy_yyyy, tg_xy_yyyz, tg_xy_yyzz, tg_xy_yzzz, tg_xyy_xxxx, tg_xyy_xxxy, tg_xyy_xxxz, tg_xyy_xxyy, tg_xyy_xxyz, tg_xyy_xxzz, tg_xyy_xyyy, tg_xyy_xyyz, tg_xyy_xyzz, tg_xyy_xzzz, tg_xyy_yyyy, tg_xyy_yyyz, tg_xyy_yyzz, tg_xyy_yzzz, tg_xyy_zzzz, tg_xyz_xxxx, tg_xyz_xxxy, tg_xyz_xxxz, tg_xyz_xxyy, tg_xyz_xxyz, tg_xyz_xxzz, tg_xyz_xyyy, tg_xyz_xyyz, tg_xyz_xyzz, tg_xyz_xzzz, tg_xyz_yyyy, tg_xyz_yyyz, tg_xyz_yyzz, tg_xyz_yzzz, tg_xyz_zzzz, tg_xz_xxxx, tg_xz_xxxz, tg_xz_xxzz, tg_xz_xzzz, tg_xz_yyyz, tg_xz_yyzz, tg_xz_yzzz, tg_xz_zzzz, tg_xzz_xxxx, tg_xzz_xxxy, tg_xzz_xxxz, tg_xzz_xxyy, tg_xzz_xxyz, tg_xzz_xxzz, tg_xzz_xyyy, tg_xzz_xyyz, tg_xzz_xyzz, tg_xzz_xzzz, tg_xzz_yyyy, tg_xzz_yyyz, tg_xzz_yyzz, tg_xzz_yzzz, tg_xzz_zzzz, tg_y_xxxx, tg_y_xxxy, tg_y_xxxz, tg_y_xxyy, tg_y_xxyz, tg_y_xxzz, tg_y_xyyy, tg_y_xyyz, tg_y_xyzz, tg_y_xzzz, tg_y_yyyy, tg_y_yyyz, tg_y_yyzz, tg_y_yzzz, tg_y_zzzz, tg_yy_xxx, tg_yy_xxxx, tg_yy_xxxy, tg_yy_xxxz, tg_yy_xxy, tg_yy_xxyy, tg_yy_xxyz, tg_yy_xxz, tg_yy_xxzz, tg_yy_xyy, tg_yy_xyyy, tg_yy_xyyz, tg_yy_xyz, tg_yy_xyzz, tg_yy_xzz, tg_yy_xzzz, tg_yy_yyy, tg_yy_yyyy, tg_yy_yyyz, tg_yy_yyz, tg_yy_yyzz, tg_yy_yzz, tg_yy_yzzz, tg_yy_zzz, tg_yy_zzzz, tg_yyy_xxxx, tg_yyy_xxxy, tg_yyy_xxxz, tg_yyy_xxyy, tg_yyy_xxyz, tg_yyy_xxzz, tg_yyy_xyyy, tg_yyy_xyyz, tg_yyy_xyzz, tg_yyy_xzzz, tg_yyy_yyyy, tg_yyy_yyyz, tg_yyy_yyzz, tg_yyy_yzzz, tg_yyy_zzzz, tg_yyz_xxxx, tg_yyz_xxxy, tg_yyz_xxxz, tg_yyz_xxyy, tg_yyz_xxyz, tg_yyz_xxzz, tg_yyz_xyyy, tg_yyz_xyyz, tg_yyz_xyzz, tg_yyz_xzzz, tg_yyz_yyyy, tg_yyz_yyyz, tg_yyz_yyzz, tg_yyz_yzzz, tg_yyz_zzzz, tg_yz_xxxz, tg_yz_xxyz, tg_yz_xxzz, tg_yz_xyyz, tg_yz_xyz, tg_yz_xyzz, tg_yz_xzzz, tg_yz_yyyy, tg_yz_yyyz, tg_yz_yyz, tg_yz_yyzz, tg_yz_yzz, tg_yz_yzzz, tg_yz_zzzz, tg_yzz_xxxx, tg_yzz_xxxy, tg_yzz_xxxz, tg_yzz_xxyy, tg_yzz_xxyz, tg_yzz_xxzz, tg_yzz_xyyy, tg_yzz_xyyz, tg_yzz_xyzz, tg_yzz_xzzz, tg_yzz_yyyy, tg_yzz_yyyz, tg_yzz_yyzz, tg_yzz_yzzz, tg_yzz_zzzz, tg_z_xxxx, tg_z_xxxy, tg_z_xxxz, tg_z_xxyy, tg_z_xxyz, tg_z_xxzz, tg_z_xyyy, tg_z_xyyz, tg_z_xyzz, tg_z_xzzz, tg_z_yyyy, tg_z_yyyz, tg_z_yyzz, tg_z_yzzz, tg_z_zzzz, tg_zz_xxx, tg_zz_xxxx, tg_zz_xxxy, tg_zz_xxxz, tg_zz_xxy, tg_zz_xxyy, tg_zz_xxyz, tg_zz_xxz, tg_zz_xxzz, tg_zz_xyy, tg_zz_xyyy, tg_zz_xyyz, tg_zz_xyz, tg_zz_xyzz, tg_zz_xzz, tg_zz_xzzz, tg_zz_yyy, tg_zz_yyyy, tg_zz_yyyz, tg_zz_yyz, tg_zz_yyzz, tg_zz_yzz, tg_zz_yzzz, tg_zz_zzz, tg_zz_zzzz, tg_zzz_xxxx, tg_zzz_xxxy, tg_zzz_xxxz, tg_zzz_xxyy, tg_zzz_xxyz, tg_zzz_xxzz, tg_zzz_xyyy, tg_zzz_xyyz, tg_zzz_xyzz, tg_zzz_xzzz, tg_zzz_yyyy, tg_zzz_yyyz, tg_zzz_yyzz, tg_zzz_yzzz, tg_zzz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_xxxx[i] = 2.0 * tg_x_xxxx[i] * fxi[i] + 4.0 * tg_xx_xxx[i] * fxi[i] + tg_xx_xxxx[i] * ra_x[i];

        tg_xxx_xxxy[i] = 2.0 * tg_x_xxxy[i] * fxi[i] + 3.0 * tg_xx_xxy[i] * fxi[i] + tg_xx_xxxy[i] * ra_x[i];

        tg_xxx_xxxz[i] = 2.0 * tg_x_xxxz[i] * fxi[i] + 3.0 * tg_xx_xxz[i] * fxi[i] + tg_xx_xxxz[i] * ra_x[i];

        tg_xxx_xxyy[i] = 2.0 * tg_x_xxyy[i] * fxi[i] + 2.0 * tg_xx_xyy[i] * fxi[i] + tg_xx_xxyy[i] * ra_x[i];

        tg_xxx_xxyz[i] = 2.0 * tg_x_xxyz[i] * fxi[i] + 2.0 * tg_xx_xyz[i] * fxi[i] + tg_xx_xxyz[i] * ra_x[i];

        tg_xxx_xxzz[i] = 2.0 * tg_x_xxzz[i] * fxi[i] + 2.0 * tg_xx_xzz[i] * fxi[i] + tg_xx_xxzz[i] * ra_x[i];

        tg_xxx_xyyy[i] = 2.0 * tg_x_xyyy[i] * fxi[i] + tg_xx_yyy[i] * fxi[i] + tg_xx_xyyy[i] * ra_x[i];

        tg_xxx_xyyz[i] = 2.0 * tg_x_xyyz[i] * fxi[i] + tg_xx_yyz[i] * fxi[i] + tg_xx_xyyz[i] * ra_x[i];

        tg_xxx_xyzz[i] = 2.0 * tg_x_xyzz[i] * fxi[i] + tg_xx_yzz[i] * fxi[i] + tg_xx_xyzz[i] * ra_x[i];

        tg_xxx_xzzz[i] = 2.0 * tg_x_xzzz[i] * fxi[i] + tg_xx_zzz[i] * fxi[i] + tg_xx_xzzz[i] * ra_x[i];

        tg_xxx_yyyy[i] = 2.0 * tg_x_yyyy[i] * fxi[i] + tg_xx_yyyy[i] * ra_x[i];

        tg_xxx_yyyz[i] = 2.0 * tg_x_yyyz[i] * fxi[i] + tg_xx_yyyz[i] * ra_x[i];

        tg_xxx_yyzz[i] = 2.0 * tg_x_yyzz[i] * fxi[i] + tg_xx_yyzz[i] * ra_x[i];

        tg_xxx_yzzz[i] = 2.0 * tg_x_yzzz[i] * fxi[i] + tg_xx_yzzz[i] * ra_x[i];

        tg_xxx_zzzz[i] = 2.0 * tg_x_zzzz[i] * fxi[i] + tg_xx_zzzz[i] * ra_x[i];

        tg_xxy_xxxx[i] = tg_xx_xxxx[i] * ra_y[i];

        tg_xxy_xxxy[i] = tg_xx_xxx[i] * fxi[i] + tg_xx_xxxy[i] * ra_y[i];

        tg_xxy_xxxz[i] = tg_xx_xxxz[i] * ra_y[i];

        tg_xxy_xxyy[i] = 2.0 * tg_xx_xxy[i] * fxi[i] + tg_xx_xxyy[i] * ra_y[i];

        tg_xxy_xxyz[i] = tg_xx_xxz[i] * fxi[i] + tg_xx_xxyz[i] * ra_y[i];

        tg_xxy_xxzz[i] = tg_xx_xxzz[i] * ra_y[i];

        tg_xxy_xyyy[i] = 3.0 * tg_xx_xyy[i] * fxi[i] + tg_xx_xyyy[i] * ra_y[i];

        tg_xxy_xyyz[i] = 2.0 * tg_xx_xyz[i] * fxi[i] + tg_xx_xyyz[i] * ra_y[i];

        tg_xxy_xyzz[i] = tg_xx_xzz[i] * fxi[i] + tg_xx_xyzz[i] * ra_y[i];

        tg_xxy_xzzz[i] = tg_xx_xzzz[i] * ra_y[i];

        tg_xxy_yyyy[i] = tg_y_yyyy[i] * fxi[i] + tg_xy_yyyy[i] * ra_x[i];

        tg_xxy_yyyz[i] = tg_y_yyyz[i] * fxi[i] + tg_xy_yyyz[i] * ra_x[i];

        tg_xxy_yyzz[i] = tg_y_yyzz[i] * fxi[i] + tg_xy_yyzz[i] * ra_x[i];

        tg_xxy_yzzz[i] = tg_y_yzzz[i] * fxi[i] + tg_xy_yzzz[i] * ra_x[i];

        tg_xxy_zzzz[i] = tg_xx_zzzz[i] * ra_y[i];

        tg_xxz_xxxx[i] = tg_xx_xxxx[i] * ra_z[i];

        tg_xxz_xxxy[i] = tg_xx_xxxy[i] * ra_z[i];

        tg_xxz_xxxz[i] = tg_xx_xxx[i] * fxi[i] + tg_xx_xxxz[i] * ra_z[i];

        tg_xxz_xxyy[i] = tg_xx_xxyy[i] * ra_z[i];

        tg_xxz_xxyz[i] = tg_xx_xxy[i] * fxi[i] + tg_xx_xxyz[i] * ra_z[i];

        tg_xxz_xxzz[i] = 2.0 * tg_xx_xxz[i] * fxi[i] + tg_xx_xxzz[i] * ra_z[i];

        tg_xxz_xyyy[i] = tg_xx_xyyy[i] * ra_z[i];

        tg_xxz_xyyz[i] = tg_xx_xyy[i] * fxi[i] + tg_xx_xyyz[i] * ra_z[i];

        tg_xxz_xyzz[i] = 2.0 * tg_xx_xyz[i] * fxi[i] + tg_xx_xyzz[i] * ra_z[i];

        tg_xxz_xzzz[i] = 3.0 * tg_xx_xzz[i] * fxi[i] + tg_xx_xzzz[i] * ra_z[i];

        tg_xxz_yyyy[i] = tg_xx_yyyy[i] * ra_z[i];

        tg_xxz_yyyz[i] = tg_z_yyyz[i] * fxi[i] + tg_xz_yyyz[i] * ra_x[i];

        tg_xxz_yyzz[i] = tg_z_yyzz[i] * fxi[i] + tg_xz_yyzz[i] * ra_x[i];

        tg_xxz_yzzz[i] = tg_z_yzzz[i] * fxi[i] + tg_xz_yzzz[i] * ra_x[i];

        tg_xxz_zzzz[i] = tg_z_zzzz[i] * fxi[i] + tg_xz_zzzz[i] * ra_x[i];

        tg_xyy_xxxx[i] = 4.0 * tg_yy_xxx[i] * fxi[i] + tg_yy_xxxx[i] * ra_x[i];

        tg_xyy_xxxy[i] = 3.0 * tg_yy_xxy[i] * fxi[i] + tg_yy_xxxy[i] * ra_x[i];

        tg_xyy_xxxz[i] = 3.0 * tg_yy_xxz[i] * fxi[i] + tg_yy_xxxz[i] * ra_x[i];

        tg_xyy_xxyy[i] = 2.0 * tg_yy_xyy[i] * fxi[i] + tg_yy_xxyy[i] * ra_x[i];

        tg_xyy_xxyz[i] = 2.0 * tg_yy_xyz[i] * fxi[i] + tg_yy_xxyz[i] * ra_x[i];

        tg_xyy_xxzz[i] = 2.0 * tg_yy_xzz[i] * fxi[i] + tg_yy_xxzz[i] * ra_x[i];

        tg_xyy_xyyy[i] = tg_yy_yyy[i] * fxi[i] + tg_yy_xyyy[i] * ra_x[i];

        tg_xyy_xyyz[i] = tg_yy_yyz[i] * fxi[i] + tg_yy_xyyz[i] * ra_x[i];

        tg_xyy_xyzz[i] = tg_yy_yzz[i] * fxi[i] + tg_yy_xyzz[i] * ra_x[i];

        tg_xyy_xzzz[i] = tg_yy_zzz[i] * fxi[i] + tg_yy_xzzz[i] * ra_x[i];

        tg_xyy_yyyy[i] = tg_yy_yyyy[i] * ra_x[i];

        tg_xyy_yyyz[i] = tg_yy_yyyz[i] * ra_x[i];

        tg_xyy_yyzz[i] = tg_yy_yyzz[i] * ra_x[i];

        tg_xyy_yzzz[i] = tg_yy_yzzz[i] * ra_x[i];

        tg_xyy_zzzz[i] = tg_yy_zzzz[i] * ra_x[i];

        tg_xyz_xxxx[i] = tg_xz_xxxx[i] * ra_y[i];

        tg_xyz_xxxy[i] = tg_xy_xxxy[i] * ra_z[i];

        tg_xyz_xxxz[i] = tg_xz_xxxz[i] * ra_y[i];

        tg_xyz_xxyy[i] = tg_xy_xxyy[i] * ra_z[i];

        tg_xyz_xxyz[i] = 2.0 * tg_yz_xyz[i] * fxi[i] + tg_yz_xxyz[i] * ra_x[i];

        tg_xyz_xxzz[i] = tg_xz_xxzz[i] * ra_y[i];

        tg_xyz_xyyy[i] = tg_xy_xyyy[i] * ra_z[i];

        tg_xyz_xyyz[i] = tg_yz_yyz[i] * fxi[i] + tg_yz_xyyz[i] * ra_x[i];

        tg_xyz_xyzz[i] = tg_yz_yzz[i] * fxi[i] + tg_yz_xyzz[i] * ra_x[i];

        tg_xyz_xzzz[i] = tg_xz_xzzz[i] * ra_y[i];

        tg_xyz_yyyy[i] = tg_yz_yyyy[i] * ra_x[i];

        tg_xyz_yyyz[i] = tg_yz_yyyz[i] * ra_x[i];

        tg_xyz_yyzz[i] = tg_yz_yyzz[i] * ra_x[i];

        tg_xyz_yzzz[i] = tg_yz_yzzz[i] * ra_x[i];

        tg_xyz_zzzz[i] = tg_yz_zzzz[i] * ra_x[i];

        tg_xzz_xxxx[i] = 4.0 * tg_zz_xxx[i] * fxi[i] + tg_zz_xxxx[i] * ra_x[i];

        tg_xzz_xxxy[i] = 3.0 * tg_zz_xxy[i] * fxi[i] + tg_zz_xxxy[i] * ra_x[i];

        tg_xzz_xxxz[i] = 3.0 * tg_zz_xxz[i] * fxi[i] + tg_zz_xxxz[i] * ra_x[i];

        tg_xzz_xxyy[i] = 2.0 * tg_zz_xyy[i] * fxi[i] + tg_zz_xxyy[i] * ra_x[i];

        tg_xzz_xxyz[i] = 2.0 * tg_zz_xyz[i] * fxi[i] + tg_zz_xxyz[i] * ra_x[i];

        tg_xzz_xxzz[i] = 2.0 * tg_zz_xzz[i] * fxi[i] + tg_zz_xxzz[i] * ra_x[i];

        tg_xzz_xyyy[i] = tg_zz_yyy[i] * fxi[i] + tg_zz_xyyy[i] * ra_x[i];

        tg_xzz_xyyz[i] = tg_zz_yyz[i] * fxi[i] + tg_zz_xyyz[i] * ra_x[i];

        tg_xzz_xyzz[i] = tg_zz_yzz[i] * fxi[i] + tg_zz_xyzz[i] * ra_x[i];

        tg_xzz_xzzz[i] = tg_zz_zzz[i] * fxi[i] + tg_zz_xzzz[i] * ra_x[i];

        tg_xzz_yyyy[i] = tg_zz_yyyy[i] * ra_x[i];

        tg_xzz_yyyz[i] = tg_zz_yyyz[i] * ra_x[i];

        tg_xzz_yyzz[i] = tg_zz_yyzz[i] * ra_x[i];

        tg_xzz_yzzz[i] = tg_zz_yzzz[i] * ra_x[i];

        tg_xzz_zzzz[i] = tg_zz_zzzz[i] * ra_x[i];

        tg_yyy_xxxx[i] = 2.0 * tg_y_xxxx[i] * fxi[i] + tg_yy_xxxx[i] * ra_y[i];

        tg_yyy_xxxy[i] = 2.0 * tg_y_xxxy[i] * fxi[i] + tg_yy_xxx[i] * fxi[i] + tg_yy_xxxy[i] * ra_y[i];

        tg_yyy_xxxz[i] = 2.0 * tg_y_xxxz[i] * fxi[i] + tg_yy_xxxz[i] * ra_y[i];

        tg_yyy_xxyy[i] = 2.0 * tg_y_xxyy[i] * fxi[i] + 2.0 * tg_yy_xxy[i] * fxi[i] + tg_yy_xxyy[i] * ra_y[i];

        tg_yyy_xxyz[i] = 2.0 * tg_y_xxyz[i] * fxi[i] + tg_yy_xxz[i] * fxi[i] + tg_yy_xxyz[i] * ra_y[i];

        tg_yyy_xxzz[i] = 2.0 * tg_y_xxzz[i] * fxi[i] + tg_yy_xxzz[i] * ra_y[i];

        tg_yyy_xyyy[i] = 2.0 * tg_y_xyyy[i] * fxi[i] + 3.0 * tg_yy_xyy[i] * fxi[i] + tg_yy_xyyy[i] * ra_y[i];

        tg_yyy_xyyz[i] = 2.0 * tg_y_xyyz[i] * fxi[i] + 2.0 * tg_yy_xyz[i] * fxi[i] + tg_yy_xyyz[i] * ra_y[i];

        tg_yyy_xyzz[i] = 2.0 * tg_y_xyzz[i] * fxi[i] + tg_yy_xzz[i] * fxi[i] + tg_yy_xyzz[i] * ra_y[i];

        tg_yyy_xzzz[i] = 2.0 * tg_y_xzzz[i] * fxi[i] + tg_yy_xzzz[i] * ra_y[i];

        tg_yyy_yyyy[i] = 2.0 * tg_y_yyyy[i] * fxi[i] + 4.0 * tg_yy_yyy[i] * fxi[i] + tg_yy_yyyy[i] * ra_y[i];

        tg_yyy_yyyz[i] = 2.0 * tg_y_yyyz[i] * fxi[i] + 3.0 * tg_yy_yyz[i] * fxi[i] + tg_yy_yyyz[i] * ra_y[i];

        tg_yyy_yyzz[i] = 2.0 * tg_y_yyzz[i] * fxi[i] + 2.0 * tg_yy_yzz[i] * fxi[i] + tg_yy_yyzz[i] * ra_y[i];

        tg_yyy_yzzz[i] = 2.0 * tg_y_yzzz[i] * fxi[i] + tg_yy_zzz[i] * fxi[i] + tg_yy_yzzz[i] * ra_y[i];

        tg_yyy_zzzz[i] = 2.0 * tg_y_zzzz[i] * fxi[i] + tg_yy_zzzz[i] * ra_y[i];

        tg_yyz_xxxx[i] = tg_yy_xxxx[i] * ra_z[i];

        tg_yyz_xxxy[i] = tg_yy_xxxy[i] * ra_z[i];

        tg_yyz_xxxz[i] = tg_z_xxxz[i] * fxi[i] + tg_yz_xxxz[i] * ra_y[i];

        tg_yyz_xxyy[i] = tg_yy_xxyy[i] * ra_z[i];

        tg_yyz_xxyz[i] = tg_yy_xxy[i] * fxi[i] + tg_yy_xxyz[i] * ra_z[i];

        tg_yyz_xxzz[i] = tg_z_xxzz[i] * fxi[i] + tg_yz_xxzz[i] * ra_y[i];

        tg_yyz_xyyy[i] = tg_yy_xyyy[i] * ra_z[i];

        tg_yyz_xyyz[i] = tg_yy_xyy[i] * fxi[i] + tg_yy_xyyz[i] * ra_z[i];

        tg_yyz_xyzz[i] = 2.0 * tg_yy_xyz[i] * fxi[i] + tg_yy_xyzz[i] * ra_z[i];

        tg_yyz_xzzz[i] = tg_z_xzzz[i] * fxi[i] + tg_yz_xzzz[i] * ra_y[i];

        tg_yyz_yyyy[i] = tg_yy_yyyy[i] * ra_z[i];

        tg_yyz_yyyz[i] = tg_yy_yyy[i] * fxi[i] + tg_yy_yyyz[i] * ra_z[i];

        tg_yyz_yyzz[i] = 2.0 * tg_yy_yyz[i] * fxi[i] + tg_yy_yyzz[i] * ra_z[i];

        tg_yyz_yzzz[i] = 3.0 * tg_yy_yzz[i] * fxi[i] + tg_yy_yzzz[i] * ra_z[i];

        tg_yyz_zzzz[i] = tg_z_zzzz[i] * fxi[i] + tg_yz_zzzz[i] * ra_y[i];

        tg_yzz_xxxx[i] = tg_zz_xxxx[i] * ra_y[i];

        tg_yzz_xxxy[i] = tg_zz_xxx[i] * fxi[i] + tg_zz_xxxy[i] * ra_y[i];

        tg_yzz_xxxz[i] = tg_zz_xxxz[i] * ra_y[i];

        tg_yzz_xxyy[i] = 2.0 * tg_zz_xxy[i] * fxi[i] + tg_zz_xxyy[i] * ra_y[i];

        tg_yzz_xxyz[i] = tg_zz_xxz[i] * fxi[i] + tg_zz_xxyz[i] * ra_y[i];

        tg_yzz_xxzz[i] = tg_zz_xxzz[i] * ra_y[i];

        tg_yzz_xyyy[i] = 3.0 * tg_zz_xyy[i] * fxi[i] + tg_zz_xyyy[i] * ra_y[i];

        tg_yzz_xyyz[i] = 2.0 * tg_zz_xyz[i] * fxi[i] + tg_zz_xyyz[i] * ra_y[i];

        tg_yzz_xyzz[i] = tg_zz_xzz[i] * fxi[i] + tg_zz_xyzz[i] * ra_y[i];

        tg_yzz_xzzz[i] = tg_zz_xzzz[i] * ra_y[i];

        tg_yzz_yyyy[i] = 4.0 * tg_zz_yyy[i] * fxi[i] + tg_zz_yyyy[i] * ra_y[i];

        tg_yzz_yyyz[i] = 3.0 * tg_zz_yyz[i] * fxi[i] + tg_zz_yyyz[i] * ra_y[i];

        tg_yzz_yyzz[i] = 2.0 * tg_zz_yzz[i] * fxi[i] + tg_zz_yyzz[i] * ra_y[i];

        tg_yzz_yzzz[i] = tg_zz_zzz[i] * fxi[i] + tg_zz_yzzz[i] * ra_y[i];

        tg_yzz_zzzz[i] = tg_zz_zzzz[i] * ra_y[i];

        tg_zzz_xxxx[i] = 2.0 * tg_z_xxxx[i] * fxi[i] + tg_zz_xxxx[i] * ra_z[i];

        tg_zzz_xxxy[i] = 2.0 * tg_z_xxxy[i] * fxi[i] + tg_zz_xxxy[i] * ra_z[i];

        tg_zzz_xxxz[i] = 2.0 * tg_z_xxxz[i] * fxi[i] + tg_zz_xxx[i] * fxi[i] + tg_zz_xxxz[i] * ra_z[i];

        tg_zzz_xxyy[i] = 2.0 * tg_z_xxyy[i] * fxi[i] + tg_zz_xxyy[i] * ra_z[i];

        tg_zzz_xxyz[i] = 2.0 * tg_z_xxyz[i] * fxi[i] + tg_zz_xxy[i] * fxi[i] + tg_zz_xxyz[i] * ra_z[i];

        tg_zzz_xxzz[i] = 2.0 * tg_z_xxzz[i] * fxi[i] + 2.0 * tg_zz_xxz[i] * fxi[i] + tg_zz_xxzz[i] * ra_z[i];

        tg_zzz_xyyy[i] = 2.0 * tg_z_xyyy[i] * fxi[i] + tg_zz_xyyy[i] * ra_z[i];

        tg_zzz_xyyz[i] = 2.0 * tg_z_xyyz[i] * fxi[i] + tg_zz_xyy[i] * fxi[i] + tg_zz_xyyz[i] * ra_z[i];

        tg_zzz_xyzz[i] = 2.0 * tg_z_xyzz[i] * fxi[i] + 2.0 * tg_zz_xyz[i] * fxi[i] + tg_zz_xyzz[i] * ra_z[i];

        tg_zzz_xzzz[i] = 2.0 * tg_z_xzzz[i] * fxi[i] + 3.0 * tg_zz_xzz[i] * fxi[i] + tg_zz_xzzz[i] * ra_z[i];

        tg_zzz_yyyy[i] = 2.0 * tg_z_yyyy[i] * fxi[i] + tg_zz_yyyy[i] * ra_z[i];

        tg_zzz_yyyz[i] = 2.0 * tg_z_yyyz[i] * fxi[i] + tg_zz_yyy[i] * fxi[i] + tg_zz_yyyz[i] * ra_z[i];

        tg_zzz_yyzz[i] = 2.0 * tg_z_yyzz[i] * fxi[i] + 2.0 * tg_zz_yyz[i] * fxi[i] + tg_zz_yyzz[i] * ra_z[i];

        tg_zzz_yzzz[i] = 2.0 * tg_z_yzzz[i] * fxi[i] + 3.0 * tg_zz_yzz[i] * fxi[i] + tg_zz_yzzz[i] * ra_z[i];

        tg_zzz_zzzz[i] = 2.0 * tg_z_zzzz[i] * fxi[i] + 4.0 * tg_zz_zzz[i] * fxi[i] + tg_zz_zzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

