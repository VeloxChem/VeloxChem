#include "LocalCorePotentialPrimRecGF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gf(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gf,
                                  const size_t idx_df,
                                  const size_t idx_fd,
                                  const size_t idx_ff,
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

    auto tg_xy_yyy = pbuffer.data(idx_df + 16);

    auto tg_xy_yyz = pbuffer.data(idx_df + 17);

    auto tg_xy_yzz = pbuffer.data(idx_df + 18);

    auto tg_xz_yyz = pbuffer.data(idx_df + 27);

    auto tg_xz_yzz = pbuffer.data(idx_df + 28);

    auto tg_xz_zzz = pbuffer.data(idx_df + 29);

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

    auto tg_yz_xxz = pbuffer.data(idx_df + 42);

    auto tg_yz_xzz = pbuffer.data(idx_df + 45);

    auto tg_yz_yyz = pbuffer.data(idx_df + 47);

    auto tg_yz_yzz = pbuffer.data(idx_df + 48);

    auto tg_yz_zzz = pbuffer.data(idx_df + 49);

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

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx = pbuffer.data(idx_fd);

    auto tg_xxx_xy = pbuffer.data(idx_fd + 1);

    auto tg_xxx_xz = pbuffer.data(idx_fd + 2);

    auto tg_xxx_yy = pbuffer.data(idx_fd + 3);

    auto tg_xxx_yz = pbuffer.data(idx_fd + 4);

    auto tg_xxx_zz = pbuffer.data(idx_fd + 5);

    auto tg_xxz_xz = pbuffer.data(idx_fd + 14);

    auto tg_xyy_xy = pbuffer.data(idx_fd + 19);

    auto tg_xyy_yy = pbuffer.data(idx_fd + 21);

    auto tg_xyy_yz = pbuffer.data(idx_fd + 22);

    auto tg_xzz_xz = pbuffer.data(idx_fd + 32);

    auto tg_xzz_yz = pbuffer.data(idx_fd + 34);

    auto tg_xzz_zz = pbuffer.data(idx_fd + 35);

    auto tg_yyy_xx = pbuffer.data(idx_fd + 36);

    auto tg_yyy_xy = pbuffer.data(idx_fd + 37);

    auto tg_yyy_xz = pbuffer.data(idx_fd + 38);

    auto tg_yyy_yy = pbuffer.data(idx_fd + 39);

    auto tg_yyy_yz = pbuffer.data(idx_fd + 40);

    auto tg_yyy_zz = pbuffer.data(idx_fd + 41);

    auto tg_yyz_xz = pbuffer.data(idx_fd + 44);

    auto tg_yyz_yz = pbuffer.data(idx_fd + 46);

    auto tg_yyz_zz = pbuffer.data(idx_fd + 47);

    auto tg_yzz_xy = pbuffer.data(idx_fd + 49);

    auto tg_yzz_xz = pbuffer.data(idx_fd + 50);

    auto tg_yzz_yy = pbuffer.data(idx_fd + 51);

    auto tg_yzz_yz = pbuffer.data(idx_fd + 52);

    auto tg_yzz_zz = pbuffer.data(idx_fd + 53);

    auto tg_zzz_xx = pbuffer.data(idx_fd + 54);

    auto tg_zzz_xy = pbuffer.data(idx_fd + 55);

    auto tg_zzz_xz = pbuffer.data(idx_fd + 56);

    auto tg_zzz_yy = pbuffer.data(idx_fd + 57);

    auto tg_zzz_yz = pbuffer.data(idx_fd + 58);

    auto tg_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx = pbuffer.data(idx_ff);

    auto tg_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto tg_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto tg_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto tg_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto tg_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto tg_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto tg_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto tg_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto tg_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto tg_xxy_xxx = pbuffer.data(idx_ff + 10);

    auto tg_xxy_xxy = pbuffer.data(idx_ff + 11);

    auto tg_xxy_xxz = pbuffer.data(idx_ff + 12);

    auto tg_xxy_xyy = pbuffer.data(idx_ff + 13);

    auto tg_xxy_xzz = pbuffer.data(idx_ff + 15);

    auto tg_xxy_yyy = pbuffer.data(idx_ff + 16);

    auto tg_xxy_yyz = pbuffer.data(idx_ff + 17);

    auto tg_xxy_yzz = pbuffer.data(idx_ff + 18);

    auto tg_xxz_xxx = pbuffer.data(idx_ff + 20);

    auto tg_xxz_xxy = pbuffer.data(idx_ff + 21);

    auto tg_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto tg_xxz_xyy = pbuffer.data(idx_ff + 23);

    auto tg_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto tg_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto tg_xxz_yyz = pbuffer.data(idx_ff + 27);

    auto tg_xxz_yzz = pbuffer.data(idx_ff + 28);

    auto tg_xxz_zzz = pbuffer.data(idx_ff + 29);

    auto tg_xyy_xxx = pbuffer.data(idx_ff + 30);

    auto tg_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto tg_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto tg_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto tg_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto tg_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto tg_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto tg_xyy_zzz = pbuffer.data(idx_ff + 39);

    auto tg_xyz_yyz = pbuffer.data(idx_ff + 47);

    auto tg_xyz_yzz = pbuffer.data(idx_ff + 48);

    auto tg_xzz_xxx = pbuffer.data(idx_ff + 50);

    auto tg_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto tg_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto tg_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto tg_xzz_yyy = pbuffer.data(idx_ff + 56);

    auto tg_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto tg_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto tg_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto tg_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto tg_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto tg_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto tg_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto tg_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto tg_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto tg_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto tg_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto tg_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto tg_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto tg_yyz_xxy = pbuffer.data(idx_ff + 71);

    auto tg_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto tg_yyz_xyy = pbuffer.data(idx_ff + 73);

    auto tg_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto tg_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto tg_yyz_yyy = pbuffer.data(idx_ff + 76);

    auto tg_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto tg_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto tg_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto tg_yzz_xxx = pbuffer.data(idx_ff + 80);

    auto tg_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto tg_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto tg_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto tg_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto tg_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto tg_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto tg_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto tg_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto tg_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto tg_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto tg_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto tg_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto tg_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto tg_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto tg_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto tg_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto tg_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto tg_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto tg_zzz_zzz = pbuffer.data(idx_ff + 99);

    // Set up components of targeted buffer : GF

    auto tg_xxxx_xxx = pbuffer.data(idx_gf);

    auto tg_xxxx_xxy = pbuffer.data(idx_gf + 1);

    auto tg_xxxx_xxz = pbuffer.data(idx_gf + 2);

    auto tg_xxxx_xyy = pbuffer.data(idx_gf + 3);

    auto tg_xxxx_xyz = pbuffer.data(idx_gf + 4);

    auto tg_xxxx_xzz = pbuffer.data(idx_gf + 5);

    auto tg_xxxx_yyy = pbuffer.data(idx_gf + 6);

    auto tg_xxxx_yyz = pbuffer.data(idx_gf + 7);

    auto tg_xxxx_yzz = pbuffer.data(idx_gf + 8);

    auto tg_xxxx_zzz = pbuffer.data(idx_gf + 9);

    auto tg_xxxy_xxx = pbuffer.data(idx_gf + 10);

    auto tg_xxxy_xxy = pbuffer.data(idx_gf + 11);

    auto tg_xxxy_xxz = pbuffer.data(idx_gf + 12);

    auto tg_xxxy_xyy = pbuffer.data(idx_gf + 13);

    auto tg_xxxy_xyz = pbuffer.data(idx_gf + 14);

    auto tg_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto tg_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto tg_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto tg_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto tg_xxxy_zzz = pbuffer.data(idx_gf + 19);

    auto tg_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto tg_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto tg_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto tg_xxxz_xyy = pbuffer.data(idx_gf + 23);

    auto tg_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto tg_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto tg_xxxz_yyy = pbuffer.data(idx_gf + 26);

    auto tg_xxxz_yyz = pbuffer.data(idx_gf + 27);

    auto tg_xxxz_yzz = pbuffer.data(idx_gf + 28);

    auto tg_xxxz_zzz = pbuffer.data(idx_gf + 29);

    auto tg_xxyy_xxx = pbuffer.data(idx_gf + 30);

    auto tg_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto tg_xxyy_xxz = pbuffer.data(idx_gf + 32);

    auto tg_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto tg_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto tg_xxyy_xzz = pbuffer.data(idx_gf + 35);

    auto tg_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto tg_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto tg_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto tg_xxyy_zzz = pbuffer.data(idx_gf + 39);

    auto tg_xxyz_xxx = pbuffer.data(idx_gf + 40);

    auto tg_xxyz_xxy = pbuffer.data(idx_gf + 41);

    auto tg_xxyz_xxz = pbuffer.data(idx_gf + 42);

    auto tg_xxyz_xyy = pbuffer.data(idx_gf + 43);

    auto tg_xxyz_xyz = pbuffer.data(idx_gf + 44);

    auto tg_xxyz_xzz = pbuffer.data(idx_gf + 45);

    auto tg_xxyz_yyy = pbuffer.data(idx_gf + 46);

    auto tg_xxyz_yyz = pbuffer.data(idx_gf + 47);

    auto tg_xxyz_yzz = pbuffer.data(idx_gf + 48);

    auto tg_xxyz_zzz = pbuffer.data(idx_gf + 49);

    auto tg_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto tg_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto tg_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto tg_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto tg_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto tg_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto tg_xxzz_yyy = pbuffer.data(idx_gf + 56);

    auto tg_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto tg_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto tg_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto tg_xyyy_xxx = pbuffer.data(idx_gf + 60);

    auto tg_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto tg_xyyy_xxz = pbuffer.data(idx_gf + 62);

    auto tg_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto tg_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto tg_xyyy_xzz = pbuffer.data(idx_gf + 65);

    auto tg_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto tg_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto tg_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto tg_xyyy_zzz = pbuffer.data(idx_gf + 69);

    auto tg_xyyz_xxx = pbuffer.data(idx_gf + 70);

    auto tg_xyyz_xxy = pbuffer.data(idx_gf + 71);

    auto tg_xyyz_xxz = pbuffer.data(idx_gf + 72);

    auto tg_xyyz_xyy = pbuffer.data(idx_gf + 73);

    auto tg_xyyz_xyz = pbuffer.data(idx_gf + 74);

    auto tg_xyyz_xzz = pbuffer.data(idx_gf + 75);

    auto tg_xyyz_yyy = pbuffer.data(idx_gf + 76);

    auto tg_xyyz_yyz = pbuffer.data(idx_gf + 77);

    auto tg_xyyz_yzz = pbuffer.data(idx_gf + 78);

    auto tg_xyyz_zzz = pbuffer.data(idx_gf + 79);

    auto tg_xyzz_xxx = pbuffer.data(idx_gf + 80);

    auto tg_xyzz_xxy = pbuffer.data(idx_gf + 81);

    auto tg_xyzz_xxz = pbuffer.data(idx_gf + 82);

    auto tg_xyzz_xyy = pbuffer.data(idx_gf + 83);

    auto tg_xyzz_xyz = pbuffer.data(idx_gf + 84);

    auto tg_xyzz_xzz = pbuffer.data(idx_gf + 85);

    auto tg_xyzz_yyy = pbuffer.data(idx_gf + 86);

    auto tg_xyzz_yyz = pbuffer.data(idx_gf + 87);

    auto tg_xyzz_yzz = pbuffer.data(idx_gf + 88);

    auto tg_xyzz_zzz = pbuffer.data(idx_gf + 89);

    auto tg_xzzz_xxx = pbuffer.data(idx_gf + 90);

    auto tg_xzzz_xxy = pbuffer.data(idx_gf + 91);

    auto tg_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto tg_xzzz_xyy = pbuffer.data(idx_gf + 93);

    auto tg_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto tg_xzzz_xzz = pbuffer.data(idx_gf + 95);

    auto tg_xzzz_yyy = pbuffer.data(idx_gf + 96);

    auto tg_xzzz_yyz = pbuffer.data(idx_gf + 97);

    auto tg_xzzz_yzz = pbuffer.data(idx_gf + 98);

    auto tg_xzzz_zzz = pbuffer.data(idx_gf + 99);

    auto tg_yyyy_xxx = pbuffer.data(idx_gf + 100);

    auto tg_yyyy_xxy = pbuffer.data(idx_gf + 101);

    auto tg_yyyy_xxz = pbuffer.data(idx_gf + 102);

    auto tg_yyyy_xyy = pbuffer.data(idx_gf + 103);

    auto tg_yyyy_xyz = pbuffer.data(idx_gf + 104);

    auto tg_yyyy_xzz = pbuffer.data(idx_gf + 105);

    auto tg_yyyy_yyy = pbuffer.data(idx_gf + 106);

    auto tg_yyyy_yyz = pbuffer.data(idx_gf + 107);

    auto tg_yyyy_yzz = pbuffer.data(idx_gf + 108);

    auto tg_yyyy_zzz = pbuffer.data(idx_gf + 109);

    auto tg_yyyz_xxx = pbuffer.data(idx_gf + 110);

    auto tg_yyyz_xxy = pbuffer.data(idx_gf + 111);

    auto tg_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto tg_yyyz_xyy = pbuffer.data(idx_gf + 113);

    auto tg_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto tg_yyyz_xzz = pbuffer.data(idx_gf + 115);

    auto tg_yyyz_yyy = pbuffer.data(idx_gf + 116);

    auto tg_yyyz_yyz = pbuffer.data(idx_gf + 117);

    auto tg_yyyz_yzz = pbuffer.data(idx_gf + 118);

    auto tg_yyyz_zzz = pbuffer.data(idx_gf + 119);

    auto tg_yyzz_xxx = pbuffer.data(idx_gf + 120);

    auto tg_yyzz_xxy = pbuffer.data(idx_gf + 121);

    auto tg_yyzz_xxz = pbuffer.data(idx_gf + 122);

    auto tg_yyzz_xyy = pbuffer.data(idx_gf + 123);

    auto tg_yyzz_xyz = pbuffer.data(idx_gf + 124);

    auto tg_yyzz_xzz = pbuffer.data(idx_gf + 125);

    auto tg_yyzz_yyy = pbuffer.data(idx_gf + 126);

    auto tg_yyzz_yyz = pbuffer.data(idx_gf + 127);

    auto tg_yyzz_yzz = pbuffer.data(idx_gf + 128);

    auto tg_yyzz_zzz = pbuffer.data(idx_gf + 129);

    auto tg_yzzz_xxx = pbuffer.data(idx_gf + 130);

    auto tg_yzzz_xxy = pbuffer.data(idx_gf + 131);

    auto tg_yzzz_xxz = pbuffer.data(idx_gf + 132);

    auto tg_yzzz_xyy = pbuffer.data(idx_gf + 133);

    auto tg_yzzz_xyz = pbuffer.data(idx_gf + 134);

    auto tg_yzzz_xzz = pbuffer.data(idx_gf + 135);

    auto tg_yzzz_yyy = pbuffer.data(idx_gf + 136);

    auto tg_yzzz_yyz = pbuffer.data(idx_gf + 137);

    auto tg_yzzz_yzz = pbuffer.data(idx_gf + 138);

    auto tg_yzzz_zzz = pbuffer.data(idx_gf + 139);

    auto tg_zzzz_xxx = pbuffer.data(idx_gf + 140);

    auto tg_zzzz_xxy = pbuffer.data(idx_gf + 141);

    auto tg_zzzz_xxz = pbuffer.data(idx_gf + 142);

    auto tg_zzzz_xyy = pbuffer.data(idx_gf + 143);

    auto tg_zzzz_xyz = pbuffer.data(idx_gf + 144);

    auto tg_zzzz_xzz = pbuffer.data(idx_gf + 145);

    auto tg_zzzz_yyy = pbuffer.data(idx_gf + 146);

    auto tg_zzzz_yyz = pbuffer.data(idx_gf + 147);

    auto tg_zzzz_yzz = pbuffer.data(idx_gf + 148);

    auto tg_zzzz_zzz = pbuffer.data(idx_gf + 149);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_xxx, tg_xx_xxy, tg_xx_xxz, tg_xx_xyy, tg_xx_xyz, tg_xx_xzz, tg_xx_yyy, tg_xx_yyz, tg_xx_yzz, tg_xx_zzz, tg_xxx_xx, tg_xxx_xxx, tg_xxx_xxy, tg_xxx_xxz, tg_xxx_xy, tg_xxx_xyy, tg_xxx_xyz, tg_xxx_xz, tg_xxx_xzz, tg_xxx_yy, tg_xxx_yyy, tg_xxx_yyz, tg_xxx_yz, tg_xxx_yzz, tg_xxx_zz, tg_xxx_zzz, tg_xxxx_xxx, tg_xxxx_xxy, tg_xxxx_xxz, tg_xxxx_xyy, tg_xxxx_xyz, tg_xxxx_xzz, tg_xxxx_yyy, tg_xxxx_yyz, tg_xxxx_yzz, tg_xxxx_zzz, tg_xxxy_xxx, tg_xxxy_xxy, tg_xxxy_xxz, tg_xxxy_xyy, tg_xxxy_xyz, tg_xxxy_xzz, tg_xxxy_yyy, tg_xxxy_yyz, tg_xxxy_yzz, tg_xxxy_zzz, tg_xxxz_xxx, tg_xxxz_xxy, tg_xxxz_xxz, tg_xxxz_xyy, tg_xxxz_xyz, tg_xxxz_xzz, tg_xxxz_yyy, tg_xxxz_yyz, tg_xxxz_yzz, tg_xxxz_zzz, tg_xxy_xxx, tg_xxy_xxy, tg_xxy_xxz, tg_xxy_xyy, tg_xxy_xzz, tg_xxy_yyy, tg_xxy_yyz, tg_xxy_yzz, tg_xxyy_xxx, tg_xxyy_xxy, tg_xxyy_xxz, tg_xxyy_xyy, tg_xxyy_xyz, tg_xxyy_xzz, tg_xxyy_yyy, tg_xxyy_yyz, tg_xxyy_yzz, tg_xxyy_zzz, tg_xxyz_xxx, tg_xxyz_xxy, tg_xxyz_xxz, tg_xxyz_xyy, tg_xxyz_xyz, tg_xxyz_xzz, tg_xxyz_yyy, tg_xxyz_yyz, tg_xxyz_yzz, tg_xxyz_zzz, tg_xxz_xxx, tg_xxz_xxy, tg_xxz_xxz, tg_xxz_xyy, tg_xxz_xyz, tg_xxz_xz, tg_xxz_xzz, tg_xxz_yyz, tg_xxz_yzz, tg_xxz_zzz, tg_xxzz_xxx, tg_xxzz_xxy, tg_xxzz_xxz, tg_xxzz_xyy, tg_xxzz_xyz, tg_xxzz_xzz, tg_xxzz_yyy, tg_xxzz_yyz, tg_xxzz_yzz, tg_xxzz_zzz, tg_xy_yyy, tg_xy_yyz, tg_xy_yzz, tg_xyy_xxx, tg_xyy_xxy, tg_xyy_xy, tg_xyy_xyy, tg_xyy_xyz, tg_xyy_yy, tg_xyy_yyy, tg_xyy_yyz, tg_xyy_yz, tg_xyy_yzz, tg_xyy_zzz, tg_xyyy_xxx, tg_xyyy_xxy, tg_xyyy_xxz, tg_xyyy_xyy, tg_xyyy_xyz, tg_xyyy_xzz, tg_xyyy_yyy, tg_xyyy_yyz, tg_xyyy_yzz, tg_xyyy_zzz, tg_xyyz_xxx, tg_xyyz_xxy, tg_xyyz_xxz, tg_xyyz_xyy, tg_xyyz_xyz, tg_xyyz_xzz, tg_xyyz_yyy, tg_xyyz_yyz, tg_xyyz_yzz, tg_xyyz_zzz, tg_xyz_yyz, tg_xyz_yzz, tg_xyzz_xxx, tg_xyzz_xxy, tg_xyzz_xxz, tg_xyzz_xyy, tg_xyzz_xyz, tg_xyzz_xzz, tg_xyzz_yyy, tg_xyzz_yyz, tg_xyzz_yzz, tg_xyzz_zzz, tg_xz_yyz, tg_xz_yzz, tg_xz_zzz, tg_xzz_xxx, tg_xzz_xxz, tg_xzz_xyz, tg_xzz_xz, tg_xzz_xzz, tg_xzz_yyy, tg_xzz_yyz, tg_xzz_yz, tg_xzz_yzz, tg_xzz_zz, tg_xzz_zzz, tg_xzzz_xxx, tg_xzzz_xxy, tg_xzzz_xxz, tg_xzzz_xyy, tg_xzzz_xyz, tg_xzzz_xzz, tg_xzzz_yyy, tg_xzzz_yyz, tg_xzzz_yzz, tg_xzzz_zzz, tg_yy_xxx, tg_yy_xxy, tg_yy_xxz, tg_yy_xyy, tg_yy_xyz, tg_yy_xzz, tg_yy_yyy, tg_yy_yyz, tg_yy_yzz, tg_yy_zzz, tg_yyy_xx, tg_yyy_xxx, tg_yyy_xxy, tg_yyy_xxz, tg_yyy_xy, tg_yyy_xyy, tg_yyy_xyz, tg_yyy_xz, tg_yyy_xzz, tg_yyy_yy, tg_yyy_yyy, tg_yyy_yyz, tg_yyy_yz, tg_yyy_yzz, tg_yyy_zz, tg_yyy_zzz, tg_yyyy_xxx, tg_yyyy_xxy, tg_yyyy_xxz, tg_yyyy_xyy, tg_yyyy_xyz, tg_yyyy_xzz, tg_yyyy_yyy, tg_yyyy_yyz, tg_yyyy_yzz, tg_yyyy_zzz, tg_yyyz_xxx, tg_yyyz_xxy, tg_yyyz_xxz, tg_yyyz_xyy, tg_yyyz_xyz, tg_yyyz_xzz, tg_yyyz_yyy, tg_yyyz_yyz, tg_yyyz_yzz, tg_yyyz_zzz, tg_yyz_xxy, tg_yyz_xxz, tg_yyz_xyy, tg_yyz_xyz, tg_yyz_xz, tg_yyz_xzz, tg_yyz_yyy, tg_yyz_yyz, tg_yyz_yz, tg_yyz_yzz, tg_yyz_zz, tg_yyz_zzz, tg_yyzz_xxx, tg_yyzz_xxy, tg_yyzz_xxz, tg_yyzz_xyy, tg_yyzz_xyz, tg_yyzz_xzz, tg_yyzz_yyy, tg_yyzz_yyz, tg_yyzz_yzz, tg_yyzz_zzz, tg_yz_xxz, tg_yz_xzz, tg_yz_yyz, tg_yz_yzz, tg_yz_zzz, tg_yzz_xxx, tg_yzz_xxy, tg_yzz_xxz, tg_yzz_xy, tg_yzz_xyy, tg_yzz_xyz, tg_yzz_xz, tg_yzz_xzz, tg_yzz_yy, tg_yzz_yyy, tg_yzz_yyz, tg_yzz_yz, tg_yzz_yzz, tg_yzz_zz, tg_yzz_zzz, tg_yzzz_xxx, tg_yzzz_xxy, tg_yzzz_xxz, tg_yzzz_xyy, tg_yzzz_xyz, tg_yzzz_xzz, tg_yzzz_yyy, tg_yzzz_yyz, tg_yzzz_yzz, tg_yzzz_zzz, tg_zz_xxx, tg_zz_xxy, tg_zz_xxz, tg_zz_xyy, tg_zz_xyz, tg_zz_xzz, tg_zz_yyy, tg_zz_yyz, tg_zz_yzz, tg_zz_zzz, tg_zzz_xx, tg_zzz_xxx, tg_zzz_xxy, tg_zzz_xxz, tg_zzz_xy, tg_zzz_xyy, tg_zzz_xyz, tg_zzz_xz, tg_zzz_xzz, tg_zzz_yy, tg_zzz_yyy, tg_zzz_yyz, tg_zzz_yz, tg_zzz_yzz, tg_zzz_zz, tg_zzz_zzz, tg_zzzz_xxx, tg_zzzz_xxy, tg_zzzz_xxz, tg_zzzz_xyy, tg_zzzz_xyz, tg_zzzz_xzz, tg_zzzz_yyy, tg_zzzz_yyz, tg_zzzz_yzz, tg_zzzz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_xxx[i] = 3.0 * tg_xx_xxx[i] * fxi[i] + 3.0 * tg_xxx_xx[i] * fxi[i] + tg_xxx_xxx[i] * ra_x[i];

        tg_xxxx_xxy[i] = 3.0 * tg_xx_xxy[i] * fxi[i] + 2.0 * tg_xxx_xy[i] * fxi[i] + tg_xxx_xxy[i] * ra_x[i];

        tg_xxxx_xxz[i] = 3.0 * tg_xx_xxz[i] * fxi[i] + 2.0 * tg_xxx_xz[i] * fxi[i] + tg_xxx_xxz[i] * ra_x[i];

        tg_xxxx_xyy[i] = 3.0 * tg_xx_xyy[i] * fxi[i] + tg_xxx_yy[i] * fxi[i] + tg_xxx_xyy[i] * ra_x[i];

        tg_xxxx_xyz[i] = 3.0 * tg_xx_xyz[i] * fxi[i] + tg_xxx_yz[i] * fxi[i] + tg_xxx_xyz[i] * ra_x[i];

        tg_xxxx_xzz[i] = 3.0 * tg_xx_xzz[i] * fxi[i] + tg_xxx_zz[i] * fxi[i] + tg_xxx_xzz[i] * ra_x[i];

        tg_xxxx_yyy[i] = 3.0 * tg_xx_yyy[i] * fxi[i] + tg_xxx_yyy[i] * ra_x[i];

        tg_xxxx_yyz[i] = 3.0 * tg_xx_yyz[i] * fxi[i] + tg_xxx_yyz[i] * ra_x[i];

        tg_xxxx_yzz[i] = 3.0 * tg_xx_yzz[i] * fxi[i] + tg_xxx_yzz[i] * ra_x[i];

        tg_xxxx_zzz[i] = 3.0 * tg_xx_zzz[i] * fxi[i] + tg_xxx_zzz[i] * ra_x[i];

        tg_xxxy_xxx[i] = tg_xxx_xxx[i] * ra_y[i];

        tg_xxxy_xxy[i] = tg_xxx_xx[i] * fxi[i] + tg_xxx_xxy[i] * ra_y[i];

        tg_xxxy_xxz[i] = tg_xxx_xxz[i] * ra_y[i];

        tg_xxxy_xyy[i] = 2.0 * tg_xxx_xy[i] * fxi[i] + tg_xxx_xyy[i] * ra_y[i];

        tg_xxxy_xyz[i] = tg_xxx_xz[i] * fxi[i] + tg_xxx_xyz[i] * ra_y[i];

        tg_xxxy_xzz[i] = tg_xxx_xzz[i] * ra_y[i];

        tg_xxxy_yyy[i] = 2.0 * tg_xy_yyy[i] * fxi[i] + tg_xxy_yyy[i] * ra_x[i];

        tg_xxxy_yyz[i] = 2.0 * tg_xy_yyz[i] * fxi[i] + tg_xxy_yyz[i] * ra_x[i];

        tg_xxxy_yzz[i] = 2.0 * tg_xy_yzz[i] * fxi[i] + tg_xxy_yzz[i] * ra_x[i];

        tg_xxxy_zzz[i] = tg_xxx_zzz[i] * ra_y[i];

        tg_xxxz_xxx[i] = tg_xxx_xxx[i] * ra_z[i];

        tg_xxxz_xxy[i] = tg_xxx_xxy[i] * ra_z[i];

        tg_xxxz_xxz[i] = tg_xxx_xx[i] * fxi[i] + tg_xxx_xxz[i] * ra_z[i];

        tg_xxxz_xyy[i] = tg_xxx_xyy[i] * ra_z[i];

        tg_xxxz_xyz[i] = tg_xxx_xy[i] * fxi[i] + tg_xxx_xyz[i] * ra_z[i];

        tg_xxxz_xzz[i] = 2.0 * tg_xxx_xz[i] * fxi[i] + tg_xxx_xzz[i] * ra_z[i];

        tg_xxxz_yyy[i] = tg_xxx_yyy[i] * ra_z[i];

        tg_xxxz_yyz[i] = 2.0 * tg_xz_yyz[i] * fxi[i] + tg_xxz_yyz[i] * ra_x[i];

        tg_xxxz_yzz[i] = 2.0 * tg_xz_yzz[i] * fxi[i] + tg_xxz_yzz[i] * ra_x[i];

        tg_xxxz_zzz[i] = 2.0 * tg_xz_zzz[i] * fxi[i] + tg_xxz_zzz[i] * ra_x[i];

        tg_xxyy_xxx[i] = tg_xx_xxx[i] * fxi[i] + tg_xxy_xxx[i] * ra_y[i];

        tg_xxyy_xxy[i] = tg_yy_xxy[i] * fxi[i] + 2.0 * tg_xyy_xy[i] * fxi[i] + tg_xyy_xxy[i] * ra_x[i];

        tg_xxyy_xxz[i] = tg_xx_xxz[i] * fxi[i] + tg_xxy_xxz[i] * ra_y[i];

        tg_xxyy_xyy[i] = tg_yy_xyy[i] * fxi[i] + tg_xyy_yy[i] * fxi[i] + tg_xyy_xyy[i] * ra_x[i];

        tg_xxyy_xyz[i] = tg_yy_xyz[i] * fxi[i] + tg_xyy_yz[i] * fxi[i] + tg_xyy_xyz[i] * ra_x[i];

        tg_xxyy_xzz[i] = tg_xx_xzz[i] * fxi[i] + tg_xxy_xzz[i] * ra_y[i];

        tg_xxyy_yyy[i] = tg_yy_yyy[i] * fxi[i] + tg_xyy_yyy[i] * ra_x[i];

        tg_xxyy_yyz[i] = tg_yy_yyz[i] * fxi[i] + tg_xyy_yyz[i] * ra_x[i];

        tg_xxyy_yzz[i] = tg_yy_yzz[i] * fxi[i] + tg_xyy_yzz[i] * ra_x[i];

        tg_xxyy_zzz[i] = tg_yy_zzz[i] * fxi[i] + tg_xyy_zzz[i] * ra_x[i];

        tg_xxyz_xxx[i] = tg_xxz_xxx[i] * ra_y[i];

        tg_xxyz_xxy[i] = tg_xxy_xxy[i] * ra_z[i];

        tg_xxyz_xxz[i] = tg_xxz_xxz[i] * ra_y[i];

        tg_xxyz_xyy[i] = tg_xxy_xyy[i] * ra_z[i];

        tg_xxyz_xyz[i] = tg_xxz_xz[i] * fxi[i] + tg_xxz_xyz[i] * ra_y[i];

        tg_xxyz_xzz[i] = tg_xxz_xzz[i] * ra_y[i];

        tg_xxyz_yyy[i] = tg_xxy_yyy[i] * ra_z[i];

        tg_xxyz_yyz[i] = tg_yz_yyz[i] * fxi[i] + tg_xyz_yyz[i] * ra_x[i];

        tg_xxyz_yzz[i] = tg_yz_yzz[i] * fxi[i] + tg_xyz_yzz[i] * ra_x[i];

        tg_xxyz_zzz[i] = tg_xxz_zzz[i] * ra_y[i];

        tg_xxzz_xxx[i] = tg_xx_xxx[i] * fxi[i] + tg_xxz_xxx[i] * ra_z[i];

        tg_xxzz_xxy[i] = tg_xx_xxy[i] * fxi[i] + tg_xxz_xxy[i] * ra_z[i];

        tg_xxzz_xxz[i] = tg_zz_xxz[i] * fxi[i] + 2.0 * tg_xzz_xz[i] * fxi[i] + tg_xzz_xxz[i] * ra_x[i];

        tg_xxzz_xyy[i] = tg_xx_xyy[i] * fxi[i] + tg_xxz_xyy[i] * ra_z[i];

        tg_xxzz_xyz[i] = tg_zz_xyz[i] * fxi[i] + tg_xzz_yz[i] * fxi[i] + tg_xzz_xyz[i] * ra_x[i];

        tg_xxzz_xzz[i] = tg_zz_xzz[i] * fxi[i] + tg_xzz_zz[i] * fxi[i] + tg_xzz_xzz[i] * ra_x[i];

        tg_xxzz_yyy[i] = tg_zz_yyy[i] * fxi[i] + tg_xzz_yyy[i] * ra_x[i];

        tg_xxzz_yyz[i] = tg_zz_yyz[i] * fxi[i] + tg_xzz_yyz[i] * ra_x[i];

        tg_xxzz_yzz[i] = tg_zz_yzz[i] * fxi[i] + tg_xzz_yzz[i] * ra_x[i];

        tg_xxzz_zzz[i] = tg_zz_zzz[i] * fxi[i] + tg_xzz_zzz[i] * ra_x[i];

        tg_xyyy_xxx[i] = 3.0 * tg_yyy_xx[i] * fxi[i] + tg_yyy_xxx[i] * ra_x[i];

        tg_xyyy_xxy[i] = 2.0 * tg_yyy_xy[i] * fxi[i] + tg_yyy_xxy[i] * ra_x[i];

        tg_xyyy_xxz[i] = 2.0 * tg_yyy_xz[i] * fxi[i] + tg_yyy_xxz[i] * ra_x[i];

        tg_xyyy_xyy[i] = tg_yyy_yy[i] * fxi[i] + tg_yyy_xyy[i] * ra_x[i];

        tg_xyyy_xyz[i] = tg_yyy_yz[i] * fxi[i] + tg_yyy_xyz[i] * ra_x[i];

        tg_xyyy_xzz[i] = tg_yyy_zz[i] * fxi[i] + tg_yyy_xzz[i] * ra_x[i];

        tg_xyyy_yyy[i] = tg_yyy_yyy[i] * ra_x[i];

        tg_xyyy_yyz[i] = tg_yyy_yyz[i] * ra_x[i];

        tg_xyyy_yzz[i] = tg_yyy_yzz[i] * ra_x[i];

        tg_xyyy_zzz[i] = tg_yyy_zzz[i] * ra_x[i];

        tg_xyyz_xxx[i] = tg_xyy_xxx[i] * ra_z[i];

        tg_xyyz_xxy[i] = tg_xyy_xxy[i] * ra_z[i];

        tg_xyyz_xxz[i] = 2.0 * tg_yyz_xz[i] * fxi[i] + tg_yyz_xxz[i] * ra_x[i];

        tg_xyyz_xyy[i] = tg_xyy_xyy[i] * ra_z[i];

        tg_xyyz_xyz[i] = tg_yyz_yz[i] * fxi[i] + tg_yyz_xyz[i] * ra_x[i];

        tg_xyyz_xzz[i] = tg_yyz_zz[i] * fxi[i] + tg_yyz_xzz[i] * ra_x[i];

        tg_xyyz_yyy[i] = tg_yyz_yyy[i] * ra_x[i];

        tg_xyyz_yyz[i] = tg_yyz_yyz[i] * ra_x[i];

        tg_xyyz_yzz[i] = tg_yyz_yzz[i] * ra_x[i];

        tg_xyyz_zzz[i] = tg_yyz_zzz[i] * ra_x[i];

        tg_xyzz_xxx[i] = tg_xzz_xxx[i] * ra_y[i];

        tg_xyzz_xxy[i] = 2.0 * tg_yzz_xy[i] * fxi[i] + tg_yzz_xxy[i] * ra_x[i];

        tg_xyzz_xxz[i] = tg_xzz_xxz[i] * ra_y[i];

        tg_xyzz_xyy[i] = tg_yzz_yy[i] * fxi[i] + tg_yzz_xyy[i] * ra_x[i];

        tg_xyzz_xyz[i] = tg_yzz_yz[i] * fxi[i] + tg_yzz_xyz[i] * ra_x[i];

        tg_xyzz_xzz[i] = tg_xzz_xzz[i] * ra_y[i];

        tg_xyzz_yyy[i] = tg_yzz_yyy[i] * ra_x[i];

        tg_xyzz_yyz[i] = tg_yzz_yyz[i] * ra_x[i];

        tg_xyzz_yzz[i] = tg_yzz_yzz[i] * ra_x[i];

        tg_xyzz_zzz[i] = tg_yzz_zzz[i] * ra_x[i];

        tg_xzzz_xxx[i] = 3.0 * tg_zzz_xx[i] * fxi[i] + tg_zzz_xxx[i] * ra_x[i];

        tg_xzzz_xxy[i] = 2.0 * tg_zzz_xy[i] * fxi[i] + tg_zzz_xxy[i] * ra_x[i];

        tg_xzzz_xxz[i] = 2.0 * tg_zzz_xz[i] * fxi[i] + tg_zzz_xxz[i] * ra_x[i];

        tg_xzzz_xyy[i] = tg_zzz_yy[i] * fxi[i] + tg_zzz_xyy[i] * ra_x[i];

        tg_xzzz_xyz[i] = tg_zzz_yz[i] * fxi[i] + tg_zzz_xyz[i] * ra_x[i];

        tg_xzzz_xzz[i] = tg_zzz_zz[i] * fxi[i] + tg_zzz_xzz[i] * ra_x[i];

        tg_xzzz_yyy[i] = tg_zzz_yyy[i] * ra_x[i];

        tg_xzzz_yyz[i] = tg_zzz_yyz[i] * ra_x[i];

        tg_xzzz_yzz[i] = tg_zzz_yzz[i] * ra_x[i];

        tg_xzzz_zzz[i] = tg_zzz_zzz[i] * ra_x[i];

        tg_yyyy_xxx[i] = 3.0 * tg_yy_xxx[i] * fxi[i] + tg_yyy_xxx[i] * ra_y[i];

        tg_yyyy_xxy[i] = 3.0 * tg_yy_xxy[i] * fxi[i] + tg_yyy_xx[i] * fxi[i] + tg_yyy_xxy[i] * ra_y[i];

        tg_yyyy_xxz[i] = 3.0 * tg_yy_xxz[i] * fxi[i] + tg_yyy_xxz[i] * ra_y[i];

        tg_yyyy_xyy[i] = 3.0 * tg_yy_xyy[i] * fxi[i] + 2.0 * tg_yyy_xy[i] * fxi[i] + tg_yyy_xyy[i] * ra_y[i];

        tg_yyyy_xyz[i] = 3.0 * tg_yy_xyz[i] * fxi[i] + tg_yyy_xz[i] * fxi[i] + tg_yyy_xyz[i] * ra_y[i];

        tg_yyyy_xzz[i] = 3.0 * tg_yy_xzz[i] * fxi[i] + tg_yyy_xzz[i] * ra_y[i];

        tg_yyyy_yyy[i] = 3.0 * tg_yy_yyy[i] * fxi[i] + 3.0 * tg_yyy_yy[i] * fxi[i] + tg_yyy_yyy[i] * ra_y[i];

        tg_yyyy_yyz[i] = 3.0 * tg_yy_yyz[i] * fxi[i] + 2.0 * tg_yyy_yz[i] * fxi[i] + tg_yyy_yyz[i] * ra_y[i];

        tg_yyyy_yzz[i] = 3.0 * tg_yy_yzz[i] * fxi[i] + tg_yyy_zz[i] * fxi[i] + tg_yyy_yzz[i] * ra_y[i];

        tg_yyyy_zzz[i] = 3.0 * tg_yy_zzz[i] * fxi[i] + tg_yyy_zzz[i] * ra_y[i];

        tg_yyyz_xxx[i] = tg_yyy_xxx[i] * ra_z[i];

        tg_yyyz_xxy[i] = tg_yyy_xxy[i] * ra_z[i];

        tg_yyyz_xxz[i] = 2.0 * tg_yz_xxz[i] * fxi[i] + tg_yyz_xxz[i] * ra_y[i];

        tg_yyyz_xyy[i] = tg_yyy_xyy[i] * ra_z[i];

        tg_yyyz_xyz[i] = tg_yyy_xy[i] * fxi[i] + tg_yyy_xyz[i] * ra_z[i];

        tg_yyyz_xzz[i] = 2.0 * tg_yz_xzz[i] * fxi[i] + tg_yyz_xzz[i] * ra_y[i];

        tg_yyyz_yyy[i] = tg_yyy_yyy[i] * ra_z[i];

        tg_yyyz_yyz[i] = tg_yyy_yy[i] * fxi[i] + tg_yyy_yyz[i] * ra_z[i];

        tg_yyyz_yzz[i] = 2.0 * tg_yyy_yz[i] * fxi[i] + tg_yyy_yzz[i] * ra_z[i];

        tg_yyyz_zzz[i] = 2.0 * tg_yz_zzz[i] * fxi[i] + tg_yyz_zzz[i] * ra_y[i];

        tg_yyzz_xxx[i] = tg_zz_xxx[i] * fxi[i] + tg_yzz_xxx[i] * ra_y[i];

        tg_yyzz_xxy[i] = tg_yy_xxy[i] * fxi[i] + tg_yyz_xxy[i] * ra_z[i];

        tg_yyzz_xxz[i] = tg_zz_xxz[i] * fxi[i] + tg_yzz_xxz[i] * ra_y[i];

        tg_yyzz_xyy[i] = tg_yy_xyy[i] * fxi[i] + tg_yyz_xyy[i] * ra_z[i];

        tg_yyzz_xyz[i] = tg_zz_xyz[i] * fxi[i] + tg_yzz_xz[i] * fxi[i] + tg_yzz_xyz[i] * ra_y[i];

        tg_yyzz_xzz[i] = tg_zz_xzz[i] * fxi[i] + tg_yzz_xzz[i] * ra_y[i];

        tg_yyzz_yyy[i] = tg_yy_yyy[i] * fxi[i] + tg_yyz_yyy[i] * ra_z[i];

        tg_yyzz_yyz[i] = tg_zz_yyz[i] * fxi[i] + 2.0 * tg_yzz_yz[i] * fxi[i] + tg_yzz_yyz[i] * ra_y[i];

        tg_yyzz_yzz[i] = tg_zz_yzz[i] * fxi[i] + tg_yzz_zz[i] * fxi[i] + tg_yzz_yzz[i] * ra_y[i];

        tg_yyzz_zzz[i] = tg_zz_zzz[i] * fxi[i] + tg_yzz_zzz[i] * ra_y[i];

        tg_yzzz_xxx[i] = tg_zzz_xxx[i] * ra_y[i];

        tg_yzzz_xxy[i] = tg_zzz_xx[i] * fxi[i] + tg_zzz_xxy[i] * ra_y[i];

        tg_yzzz_xxz[i] = tg_zzz_xxz[i] * ra_y[i];

        tg_yzzz_xyy[i] = 2.0 * tg_zzz_xy[i] * fxi[i] + tg_zzz_xyy[i] * ra_y[i];

        tg_yzzz_xyz[i] = tg_zzz_xz[i] * fxi[i] + tg_zzz_xyz[i] * ra_y[i];

        tg_yzzz_xzz[i] = tg_zzz_xzz[i] * ra_y[i];

        tg_yzzz_yyy[i] = 3.0 * tg_zzz_yy[i] * fxi[i] + tg_zzz_yyy[i] * ra_y[i];

        tg_yzzz_yyz[i] = 2.0 * tg_zzz_yz[i] * fxi[i] + tg_zzz_yyz[i] * ra_y[i];

        tg_yzzz_yzz[i] = tg_zzz_zz[i] * fxi[i] + tg_zzz_yzz[i] * ra_y[i];

        tg_yzzz_zzz[i] = tg_zzz_zzz[i] * ra_y[i];

        tg_zzzz_xxx[i] = 3.0 * tg_zz_xxx[i] * fxi[i] + tg_zzz_xxx[i] * ra_z[i];

        tg_zzzz_xxy[i] = 3.0 * tg_zz_xxy[i] * fxi[i] + tg_zzz_xxy[i] * ra_z[i];

        tg_zzzz_xxz[i] = 3.0 * tg_zz_xxz[i] * fxi[i] + tg_zzz_xx[i] * fxi[i] + tg_zzz_xxz[i] * ra_z[i];

        tg_zzzz_xyy[i] = 3.0 * tg_zz_xyy[i] * fxi[i] + tg_zzz_xyy[i] * ra_z[i];

        tg_zzzz_xyz[i] = 3.0 * tg_zz_xyz[i] * fxi[i] + tg_zzz_xy[i] * fxi[i] + tg_zzz_xyz[i] * ra_z[i];

        tg_zzzz_xzz[i] = 3.0 * tg_zz_xzz[i] * fxi[i] + 2.0 * tg_zzz_xz[i] * fxi[i] + tg_zzz_xzz[i] * ra_z[i];

        tg_zzzz_yyy[i] = 3.0 * tg_zz_yyy[i] * fxi[i] + tg_zzz_yyy[i] * ra_z[i];

        tg_zzzz_yyz[i] = 3.0 * tg_zz_yyz[i] * fxi[i] + tg_zzz_yy[i] * fxi[i] + tg_zzz_yyz[i] * ra_z[i];

        tg_zzzz_yzz[i] = 3.0 * tg_zz_yzz[i] * fxi[i] + 2.0 * tg_zzz_yz[i] * fxi[i] + tg_zzz_yzz[i] * ra_z[i];

        tg_zzzz_zzz[i] = 3.0 * tg_zz_zzz[i] * fxi[i] + 3.0 * tg_zzz_zz[i] * fxi[i] + tg_zzz_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

