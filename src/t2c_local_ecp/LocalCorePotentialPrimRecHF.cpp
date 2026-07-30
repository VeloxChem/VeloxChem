#include "LocalCorePotentialPrimRecHF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hf(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hf,
                                  const size_t idx_ff,
                                  const size_t idx_gd,
                                  const size_t idx_gf,
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

    auto tg_xxy_xxz = pbuffer.data(idx_ff + 12);

    auto tg_xxy_xzz = pbuffer.data(idx_ff + 15);

    auto tg_xxy_yyy = pbuffer.data(idx_ff + 16);

    auto tg_xxy_yyz = pbuffer.data(idx_ff + 17);

    auto tg_xxy_yzz = pbuffer.data(idx_ff + 18);

    auto tg_xxz_xxx = pbuffer.data(idx_ff + 20);

    auto tg_xxz_xxy = pbuffer.data(idx_ff + 21);

    auto tg_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto tg_xxz_xyy = pbuffer.data(idx_ff + 23);

    auto tg_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto tg_xxz_yyz = pbuffer.data(idx_ff + 27);

    auto tg_xxz_yzz = pbuffer.data(idx_ff + 28);

    auto tg_xxz_zzz = pbuffer.data(idx_ff + 29);

    auto tg_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto tg_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto tg_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto tg_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto tg_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto tg_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto tg_xyy_zzz = pbuffer.data(idx_ff + 39);

    auto tg_xyz_yyz = pbuffer.data(idx_ff + 47);

    auto tg_xyz_yzz = pbuffer.data(idx_ff + 48);

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

    auto tg_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto tg_yyz_yyy = pbuffer.data(idx_ff + 76);

    auto tg_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto tg_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto tg_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto tg_yzz_xxx = pbuffer.data(idx_ff + 80);

    auto tg_yzz_xxz = pbuffer.data(idx_ff + 82);

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

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx = pbuffer.data(idx_gd);

    auto tg_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto tg_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto tg_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto tg_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto tg_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto tg_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto tg_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto tg_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto tg_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto tg_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto tg_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto tg_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto tg_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto tg_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto tg_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto tg_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto tg_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto tg_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto tg_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto tg_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto tg_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto tg_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto tg_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto tg_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto tg_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto tg_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto tg_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto tg_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto tg_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto tg_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto tg_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto tg_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto tg_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto tg_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto tg_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto tg_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto tg_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto tg_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto tg_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto tg_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto tg_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto tg_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto tg_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto tg_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto tg_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto tg_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : GF

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

    auto tg_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto tg_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto tg_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto tg_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto tg_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto tg_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto tg_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto tg_xxxz_xyy = pbuffer.data(idx_gf + 23);

    auto tg_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto tg_xxxz_xzz = pbuffer.data(idx_gf + 25);

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

    auto tg_xxyz_xxz = pbuffer.data(idx_gf + 42);

    auto tg_xxyz_xzz = pbuffer.data(idx_gf + 45);

    auto tg_xxyz_yyz = pbuffer.data(idx_gf + 47);

    auto tg_xxyz_yzz = pbuffer.data(idx_gf + 48);

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

    auto tg_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto tg_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto tg_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto tg_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto tg_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto tg_xyyy_zzz = pbuffer.data(idx_gf + 69);

    auto tg_xyyz_yyz = pbuffer.data(idx_gf + 77);

    auto tg_xyyz_yzz = pbuffer.data(idx_gf + 78);

    auto tg_xyyz_zzz = pbuffer.data(idx_gf + 79);

    auto tg_xyzz_yyy = pbuffer.data(idx_gf + 86);

    auto tg_xyzz_yyz = pbuffer.data(idx_gf + 87);

    auto tg_xyzz_yzz = pbuffer.data(idx_gf + 88);

    auto tg_xzzz_xxx = pbuffer.data(idx_gf + 90);

    auto tg_xzzz_xxz = pbuffer.data(idx_gf + 92);

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

    // Set up components of targeted buffer : HF

    auto tg_xxxxx_xxx = pbuffer.data(idx_hf);

    auto tg_xxxxx_xxy = pbuffer.data(idx_hf + 1);

    auto tg_xxxxx_xxz = pbuffer.data(idx_hf + 2);

    auto tg_xxxxx_xyy = pbuffer.data(idx_hf + 3);

    auto tg_xxxxx_xyz = pbuffer.data(idx_hf + 4);

    auto tg_xxxxx_xzz = pbuffer.data(idx_hf + 5);

    auto tg_xxxxx_yyy = pbuffer.data(idx_hf + 6);

    auto tg_xxxxx_yyz = pbuffer.data(idx_hf + 7);

    auto tg_xxxxx_yzz = pbuffer.data(idx_hf + 8);

    auto tg_xxxxx_zzz = pbuffer.data(idx_hf + 9);

    auto tg_xxxxy_xxx = pbuffer.data(idx_hf + 10);

    auto tg_xxxxy_xxy = pbuffer.data(idx_hf + 11);

    auto tg_xxxxy_xxz = pbuffer.data(idx_hf + 12);

    auto tg_xxxxy_xyy = pbuffer.data(idx_hf + 13);

    auto tg_xxxxy_xyz = pbuffer.data(idx_hf + 14);

    auto tg_xxxxy_xzz = pbuffer.data(idx_hf + 15);

    auto tg_xxxxy_yyy = pbuffer.data(idx_hf + 16);

    auto tg_xxxxy_yyz = pbuffer.data(idx_hf + 17);

    auto tg_xxxxy_yzz = pbuffer.data(idx_hf + 18);

    auto tg_xxxxy_zzz = pbuffer.data(idx_hf + 19);

    auto tg_xxxxz_xxx = pbuffer.data(idx_hf + 20);

    auto tg_xxxxz_xxy = pbuffer.data(idx_hf + 21);

    auto tg_xxxxz_xxz = pbuffer.data(idx_hf + 22);

    auto tg_xxxxz_xyy = pbuffer.data(idx_hf + 23);

    auto tg_xxxxz_xyz = pbuffer.data(idx_hf + 24);

    auto tg_xxxxz_xzz = pbuffer.data(idx_hf + 25);

    auto tg_xxxxz_yyy = pbuffer.data(idx_hf + 26);

    auto tg_xxxxz_yyz = pbuffer.data(idx_hf + 27);

    auto tg_xxxxz_yzz = pbuffer.data(idx_hf + 28);

    auto tg_xxxxz_zzz = pbuffer.data(idx_hf + 29);

    auto tg_xxxyy_xxx = pbuffer.data(idx_hf + 30);

    auto tg_xxxyy_xxy = pbuffer.data(idx_hf + 31);

    auto tg_xxxyy_xxz = pbuffer.data(idx_hf + 32);

    auto tg_xxxyy_xyy = pbuffer.data(idx_hf + 33);

    auto tg_xxxyy_xyz = pbuffer.data(idx_hf + 34);

    auto tg_xxxyy_xzz = pbuffer.data(idx_hf + 35);

    auto tg_xxxyy_yyy = pbuffer.data(idx_hf + 36);

    auto tg_xxxyy_yyz = pbuffer.data(idx_hf + 37);

    auto tg_xxxyy_yzz = pbuffer.data(idx_hf + 38);

    auto tg_xxxyy_zzz = pbuffer.data(idx_hf + 39);

    auto tg_xxxyz_xxx = pbuffer.data(idx_hf + 40);

    auto tg_xxxyz_xxy = pbuffer.data(idx_hf + 41);

    auto tg_xxxyz_xxz = pbuffer.data(idx_hf + 42);

    auto tg_xxxyz_xyy = pbuffer.data(idx_hf + 43);

    auto tg_xxxyz_xyz = pbuffer.data(idx_hf + 44);

    auto tg_xxxyz_xzz = pbuffer.data(idx_hf + 45);

    auto tg_xxxyz_yyy = pbuffer.data(idx_hf + 46);

    auto tg_xxxyz_yyz = pbuffer.data(idx_hf + 47);

    auto tg_xxxyz_yzz = pbuffer.data(idx_hf + 48);

    auto tg_xxxyz_zzz = pbuffer.data(idx_hf + 49);

    auto tg_xxxzz_xxx = pbuffer.data(idx_hf + 50);

    auto tg_xxxzz_xxy = pbuffer.data(idx_hf + 51);

    auto tg_xxxzz_xxz = pbuffer.data(idx_hf + 52);

    auto tg_xxxzz_xyy = pbuffer.data(idx_hf + 53);

    auto tg_xxxzz_xyz = pbuffer.data(idx_hf + 54);

    auto tg_xxxzz_xzz = pbuffer.data(idx_hf + 55);

    auto tg_xxxzz_yyy = pbuffer.data(idx_hf + 56);

    auto tg_xxxzz_yyz = pbuffer.data(idx_hf + 57);

    auto tg_xxxzz_yzz = pbuffer.data(idx_hf + 58);

    auto tg_xxxzz_zzz = pbuffer.data(idx_hf + 59);

    auto tg_xxyyy_xxx = pbuffer.data(idx_hf + 60);

    auto tg_xxyyy_xxy = pbuffer.data(idx_hf + 61);

    auto tg_xxyyy_xxz = pbuffer.data(idx_hf + 62);

    auto tg_xxyyy_xyy = pbuffer.data(idx_hf + 63);

    auto tg_xxyyy_xyz = pbuffer.data(idx_hf + 64);

    auto tg_xxyyy_xzz = pbuffer.data(idx_hf + 65);

    auto tg_xxyyy_yyy = pbuffer.data(idx_hf + 66);

    auto tg_xxyyy_yyz = pbuffer.data(idx_hf + 67);

    auto tg_xxyyy_yzz = pbuffer.data(idx_hf + 68);

    auto tg_xxyyy_zzz = pbuffer.data(idx_hf + 69);

    auto tg_xxyyz_xxx = pbuffer.data(idx_hf + 70);

    auto tg_xxyyz_xxy = pbuffer.data(idx_hf + 71);

    auto tg_xxyyz_xxz = pbuffer.data(idx_hf + 72);

    auto tg_xxyyz_xyy = pbuffer.data(idx_hf + 73);

    auto tg_xxyyz_xyz = pbuffer.data(idx_hf + 74);

    auto tg_xxyyz_xzz = pbuffer.data(idx_hf + 75);

    auto tg_xxyyz_yyy = pbuffer.data(idx_hf + 76);

    auto tg_xxyyz_yyz = pbuffer.data(idx_hf + 77);

    auto tg_xxyyz_yzz = pbuffer.data(idx_hf + 78);

    auto tg_xxyyz_zzz = pbuffer.data(idx_hf + 79);

    auto tg_xxyzz_xxx = pbuffer.data(idx_hf + 80);

    auto tg_xxyzz_xxy = pbuffer.data(idx_hf + 81);

    auto tg_xxyzz_xxz = pbuffer.data(idx_hf + 82);

    auto tg_xxyzz_xyy = pbuffer.data(idx_hf + 83);

    auto tg_xxyzz_xyz = pbuffer.data(idx_hf + 84);

    auto tg_xxyzz_xzz = pbuffer.data(idx_hf + 85);

    auto tg_xxyzz_yyy = pbuffer.data(idx_hf + 86);

    auto tg_xxyzz_yyz = pbuffer.data(idx_hf + 87);

    auto tg_xxyzz_yzz = pbuffer.data(idx_hf + 88);

    auto tg_xxyzz_zzz = pbuffer.data(idx_hf + 89);

    auto tg_xxzzz_xxx = pbuffer.data(idx_hf + 90);

    auto tg_xxzzz_xxy = pbuffer.data(idx_hf + 91);

    auto tg_xxzzz_xxz = pbuffer.data(idx_hf + 92);

    auto tg_xxzzz_xyy = pbuffer.data(idx_hf + 93);

    auto tg_xxzzz_xyz = pbuffer.data(idx_hf + 94);

    auto tg_xxzzz_xzz = pbuffer.data(idx_hf + 95);

    auto tg_xxzzz_yyy = pbuffer.data(idx_hf + 96);

    auto tg_xxzzz_yyz = pbuffer.data(idx_hf + 97);

    auto tg_xxzzz_yzz = pbuffer.data(idx_hf + 98);

    auto tg_xxzzz_zzz = pbuffer.data(idx_hf + 99);

    auto tg_xyyyy_xxx = pbuffer.data(idx_hf + 100);

    auto tg_xyyyy_xxy = pbuffer.data(idx_hf + 101);

    auto tg_xyyyy_xxz = pbuffer.data(idx_hf + 102);

    auto tg_xyyyy_xyy = pbuffer.data(idx_hf + 103);

    auto tg_xyyyy_xyz = pbuffer.data(idx_hf + 104);

    auto tg_xyyyy_xzz = pbuffer.data(idx_hf + 105);

    auto tg_xyyyy_yyy = pbuffer.data(idx_hf + 106);

    auto tg_xyyyy_yyz = pbuffer.data(idx_hf + 107);

    auto tg_xyyyy_yzz = pbuffer.data(idx_hf + 108);

    auto tg_xyyyy_zzz = pbuffer.data(idx_hf + 109);

    auto tg_xyyyz_xxx = pbuffer.data(idx_hf + 110);

    auto tg_xyyyz_xxy = pbuffer.data(idx_hf + 111);

    auto tg_xyyyz_xxz = pbuffer.data(idx_hf + 112);

    auto tg_xyyyz_xyy = pbuffer.data(idx_hf + 113);

    auto tg_xyyyz_xyz = pbuffer.data(idx_hf + 114);

    auto tg_xyyyz_xzz = pbuffer.data(idx_hf + 115);

    auto tg_xyyyz_yyy = pbuffer.data(idx_hf + 116);

    auto tg_xyyyz_yyz = pbuffer.data(idx_hf + 117);

    auto tg_xyyyz_yzz = pbuffer.data(idx_hf + 118);

    auto tg_xyyyz_zzz = pbuffer.data(idx_hf + 119);

    auto tg_xyyzz_xxx = pbuffer.data(idx_hf + 120);

    auto tg_xyyzz_xxy = pbuffer.data(idx_hf + 121);

    auto tg_xyyzz_xxz = pbuffer.data(idx_hf + 122);

    auto tg_xyyzz_xyy = pbuffer.data(idx_hf + 123);

    auto tg_xyyzz_xyz = pbuffer.data(idx_hf + 124);

    auto tg_xyyzz_xzz = pbuffer.data(idx_hf + 125);

    auto tg_xyyzz_yyy = pbuffer.data(idx_hf + 126);

    auto tg_xyyzz_yyz = pbuffer.data(idx_hf + 127);

    auto tg_xyyzz_yzz = pbuffer.data(idx_hf + 128);

    auto tg_xyyzz_zzz = pbuffer.data(idx_hf + 129);

    auto tg_xyzzz_xxx = pbuffer.data(idx_hf + 130);

    auto tg_xyzzz_xxy = pbuffer.data(idx_hf + 131);

    auto tg_xyzzz_xxz = pbuffer.data(idx_hf + 132);

    auto tg_xyzzz_xyy = pbuffer.data(idx_hf + 133);

    auto tg_xyzzz_xyz = pbuffer.data(idx_hf + 134);

    auto tg_xyzzz_xzz = pbuffer.data(idx_hf + 135);

    auto tg_xyzzz_yyy = pbuffer.data(idx_hf + 136);

    auto tg_xyzzz_yyz = pbuffer.data(idx_hf + 137);

    auto tg_xyzzz_yzz = pbuffer.data(idx_hf + 138);

    auto tg_xyzzz_zzz = pbuffer.data(idx_hf + 139);

    auto tg_xzzzz_xxx = pbuffer.data(idx_hf + 140);

    auto tg_xzzzz_xxy = pbuffer.data(idx_hf + 141);

    auto tg_xzzzz_xxz = pbuffer.data(idx_hf + 142);

    auto tg_xzzzz_xyy = pbuffer.data(idx_hf + 143);

    auto tg_xzzzz_xyz = pbuffer.data(idx_hf + 144);

    auto tg_xzzzz_xzz = pbuffer.data(idx_hf + 145);

    auto tg_xzzzz_yyy = pbuffer.data(idx_hf + 146);

    auto tg_xzzzz_yyz = pbuffer.data(idx_hf + 147);

    auto tg_xzzzz_yzz = pbuffer.data(idx_hf + 148);

    auto tg_xzzzz_zzz = pbuffer.data(idx_hf + 149);

    auto tg_yyyyy_xxx = pbuffer.data(idx_hf + 150);

    auto tg_yyyyy_xxy = pbuffer.data(idx_hf + 151);

    auto tg_yyyyy_xxz = pbuffer.data(idx_hf + 152);

    auto tg_yyyyy_xyy = pbuffer.data(idx_hf + 153);

    auto tg_yyyyy_xyz = pbuffer.data(idx_hf + 154);

    auto tg_yyyyy_xzz = pbuffer.data(idx_hf + 155);

    auto tg_yyyyy_yyy = pbuffer.data(idx_hf + 156);

    auto tg_yyyyy_yyz = pbuffer.data(idx_hf + 157);

    auto tg_yyyyy_yzz = pbuffer.data(idx_hf + 158);

    auto tg_yyyyy_zzz = pbuffer.data(idx_hf + 159);

    auto tg_yyyyz_xxx = pbuffer.data(idx_hf + 160);

    auto tg_yyyyz_xxy = pbuffer.data(idx_hf + 161);

    auto tg_yyyyz_xxz = pbuffer.data(idx_hf + 162);

    auto tg_yyyyz_xyy = pbuffer.data(idx_hf + 163);

    auto tg_yyyyz_xyz = pbuffer.data(idx_hf + 164);

    auto tg_yyyyz_xzz = pbuffer.data(idx_hf + 165);

    auto tg_yyyyz_yyy = pbuffer.data(idx_hf + 166);

    auto tg_yyyyz_yyz = pbuffer.data(idx_hf + 167);

    auto tg_yyyyz_yzz = pbuffer.data(idx_hf + 168);

    auto tg_yyyyz_zzz = pbuffer.data(idx_hf + 169);

    auto tg_yyyzz_xxx = pbuffer.data(idx_hf + 170);

    auto tg_yyyzz_xxy = pbuffer.data(idx_hf + 171);

    auto tg_yyyzz_xxz = pbuffer.data(idx_hf + 172);

    auto tg_yyyzz_xyy = pbuffer.data(idx_hf + 173);

    auto tg_yyyzz_xyz = pbuffer.data(idx_hf + 174);

    auto tg_yyyzz_xzz = pbuffer.data(idx_hf + 175);

    auto tg_yyyzz_yyy = pbuffer.data(idx_hf + 176);

    auto tg_yyyzz_yyz = pbuffer.data(idx_hf + 177);

    auto tg_yyyzz_yzz = pbuffer.data(idx_hf + 178);

    auto tg_yyyzz_zzz = pbuffer.data(idx_hf + 179);

    auto tg_yyzzz_xxx = pbuffer.data(idx_hf + 180);

    auto tg_yyzzz_xxy = pbuffer.data(idx_hf + 181);

    auto tg_yyzzz_xxz = pbuffer.data(idx_hf + 182);

    auto tg_yyzzz_xyy = pbuffer.data(idx_hf + 183);

    auto tg_yyzzz_xyz = pbuffer.data(idx_hf + 184);

    auto tg_yyzzz_xzz = pbuffer.data(idx_hf + 185);

    auto tg_yyzzz_yyy = pbuffer.data(idx_hf + 186);

    auto tg_yyzzz_yyz = pbuffer.data(idx_hf + 187);

    auto tg_yyzzz_yzz = pbuffer.data(idx_hf + 188);

    auto tg_yyzzz_zzz = pbuffer.data(idx_hf + 189);

    auto tg_yzzzz_xxx = pbuffer.data(idx_hf + 190);

    auto tg_yzzzz_xxy = pbuffer.data(idx_hf + 191);

    auto tg_yzzzz_xxz = pbuffer.data(idx_hf + 192);

    auto tg_yzzzz_xyy = pbuffer.data(idx_hf + 193);

    auto tg_yzzzz_xyz = pbuffer.data(idx_hf + 194);

    auto tg_yzzzz_xzz = pbuffer.data(idx_hf + 195);

    auto tg_yzzzz_yyy = pbuffer.data(idx_hf + 196);

    auto tg_yzzzz_yyz = pbuffer.data(idx_hf + 197);

    auto tg_yzzzz_yzz = pbuffer.data(idx_hf + 198);

    auto tg_yzzzz_zzz = pbuffer.data(idx_hf + 199);

    auto tg_zzzzz_xxx = pbuffer.data(idx_hf + 200);

    auto tg_zzzzz_xxy = pbuffer.data(idx_hf + 201);

    auto tg_zzzzz_xxz = pbuffer.data(idx_hf + 202);

    auto tg_zzzzz_xyy = pbuffer.data(idx_hf + 203);

    auto tg_zzzzz_xyz = pbuffer.data(idx_hf + 204);

    auto tg_zzzzz_xzz = pbuffer.data(idx_hf + 205);

    auto tg_zzzzz_yyy = pbuffer.data(idx_hf + 206);

    auto tg_zzzzz_yyz = pbuffer.data(idx_hf + 207);

    auto tg_zzzzz_yzz = pbuffer.data(idx_hf + 208);

    auto tg_zzzzz_zzz = pbuffer.data(idx_hf + 209);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_xxx, tg_xxx_xxy, tg_xxx_xxz, tg_xxx_xyy, tg_xxx_xyz, tg_xxx_xzz, tg_xxx_yyy, tg_xxx_yyz, tg_xxx_yzz, tg_xxx_zzz, tg_xxxx_xx, tg_xxxx_xxx, tg_xxxx_xxy, tg_xxxx_xxz, tg_xxxx_xy, tg_xxxx_xyy, tg_xxxx_xyz, tg_xxxx_xz, tg_xxxx_xzz, tg_xxxx_yy, tg_xxxx_yyy, tg_xxxx_yyz, tg_xxxx_yz, tg_xxxx_yzz, tg_xxxx_zz, tg_xxxx_zzz, tg_xxxxx_xxx, tg_xxxxx_xxy, tg_xxxxx_xxz, tg_xxxxx_xyy, tg_xxxxx_xyz, tg_xxxxx_xzz, tg_xxxxx_yyy, tg_xxxxx_yyz, tg_xxxxx_yzz, tg_xxxxx_zzz, tg_xxxxy_xxx, tg_xxxxy_xxy, tg_xxxxy_xxz, tg_xxxxy_xyy, tg_xxxxy_xyz, tg_xxxxy_xzz, tg_xxxxy_yyy, tg_xxxxy_yyz, tg_xxxxy_yzz, tg_xxxxy_zzz, tg_xxxxz_xxx, tg_xxxxz_xxy, tg_xxxxz_xxz, tg_xxxxz_xyy, tg_xxxxz_xyz, tg_xxxxz_xzz, tg_xxxxz_yyy, tg_xxxxz_yyz, tg_xxxxz_yzz, tg_xxxxz_zzz, tg_xxxy_xxx, tg_xxxy_xxy, tg_xxxy_xxz, tg_xxxy_xyy, tg_xxxy_xzz, tg_xxxy_yyy, tg_xxxy_yyz, tg_xxxy_yzz, tg_xxxyy_xxx, tg_xxxyy_xxy, tg_xxxyy_xxz, tg_xxxyy_xyy, tg_xxxyy_xyz, tg_xxxyy_xzz, tg_xxxyy_yyy, tg_xxxyy_yyz, tg_xxxyy_yzz, tg_xxxyy_zzz, tg_xxxyz_xxx, tg_xxxyz_xxy, tg_xxxyz_xxz, tg_xxxyz_xyy, tg_xxxyz_xyz, tg_xxxyz_xzz, tg_xxxyz_yyy, tg_xxxyz_yyz, tg_xxxyz_yzz, tg_xxxyz_zzz, tg_xxxz_xxx, tg_xxxz_xxy, tg_xxxz_xxz, tg_xxxz_xyy, tg_xxxz_xyz, tg_xxxz_xz, tg_xxxz_xzz, tg_xxxz_yyz, tg_xxxz_yzz, tg_xxxz_zzz, tg_xxxzz_xxx, tg_xxxzz_xxy, tg_xxxzz_xxz, tg_xxxzz_xyy, tg_xxxzz_xyz, tg_xxxzz_xzz, tg_xxxzz_yyy, tg_xxxzz_yyz, tg_xxxzz_yzz, tg_xxxzz_zzz, tg_xxy_xxx, tg_xxy_xxz, tg_xxy_xzz, tg_xxy_yyy, tg_xxy_yyz, tg_xxy_yzz, tg_xxyy_xxx, tg_xxyy_xxy, tg_xxyy_xxz, tg_xxyy_xy, tg_xxyy_xyy, tg_xxyy_xyz, tg_xxyy_xzz, tg_xxyy_yy, tg_xxyy_yyy, tg_xxyy_yyz, tg_xxyy_yz, tg_xxyy_yzz, tg_xxyy_zzz, tg_xxyyy_xxx, tg_xxyyy_xxy, tg_xxyyy_xxz, tg_xxyyy_xyy, tg_xxyyy_xyz, tg_xxyyy_xzz, tg_xxyyy_yyy, tg_xxyyy_yyz, tg_xxyyy_yzz, tg_xxyyy_zzz, tg_xxyyz_xxx, tg_xxyyz_xxy, tg_xxyyz_xxz, tg_xxyyz_xyy, tg_xxyyz_xyz, tg_xxyyz_xzz, tg_xxyyz_yyy, tg_xxyyz_yyz, tg_xxyyz_yzz, tg_xxyyz_zzz, tg_xxyz_xxz, tg_xxyz_xzz, tg_xxyz_yyz, tg_xxyz_yzz, tg_xxyzz_xxx, tg_xxyzz_xxy, tg_xxyzz_xxz, tg_xxyzz_xyy, tg_xxyzz_xyz, tg_xxyzz_xzz, tg_xxyzz_yyy, tg_xxyzz_yyz, tg_xxyzz_yzz, tg_xxyzz_zzz, tg_xxz_xxx, tg_xxz_xxy, tg_xxz_xxz, tg_xxz_xyy, tg_xxz_xzz, tg_xxz_yyz, tg_xxz_yzz, tg_xxz_zzz, tg_xxzz_xx, tg_xxzz_xxx, tg_xxzz_xxy, tg_xxzz_xxz, tg_xxzz_xy, tg_xxzz_xyy, tg_xxzz_xyz, tg_xxzz_xz, tg_xxzz_xzz, tg_xxzz_yyy, tg_xxzz_yyz, tg_xxzz_yz, tg_xxzz_yzz, tg_xxzz_zz, tg_xxzz_zzz, tg_xxzzz_xxx, tg_xxzzz_xxy, tg_xxzzz_xxz, tg_xxzzz_xyy, tg_xxzzz_xyz, tg_xxzzz_xzz, tg_xxzzz_yyy, tg_xxzzz_yyz, tg_xxzzz_yzz, tg_xxzzz_zzz, tg_xyy_xxy, tg_xyy_xyy, tg_xyy_xyz, tg_xyy_yyy, tg_xyy_yyz, tg_xyy_yzz, tg_xyy_zzz, tg_xyyy_xxx, tg_xyyy_xxy, tg_xyyy_xy, tg_xyyy_xyy, tg_xyyy_xyz, tg_xyyy_yy, tg_xyyy_yyy, tg_xyyy_yyz, tg_xyyy_yz, tg_xyyy_yzz, tg_xyyy_zzz, tg_xyyyy_xxx, tg_xyyyy_xxy, tg_xyyyy_xxz, tg_xyyyy_xyy, tg_xyyyy_xyz, tg_xyyyy_xzz, tg_xyyyy_yyy, tg_xyyyy_yyz, tg_xyyyy_yzz, tg_xyyyy_zzz, tg_xyyyz_xxx, tg_xyyyz_xxy, tg_xyyyz_xxz, tg_xyyyz_xyy, tg_xyyyz_xyz, tg_xyyyz_xzz, tg_xyyyz_yyy, tg_xyyyz_yyz, tg_xyyyz_yzz, tg_xyyyz_zzz, tg_xyyz_yyz, tg_xyyz_yzz, tg_xyyz_zzz, tg_xyyzz_xxx, tg_xyyzz_xxy, tg_xyyzz_xxz, tg_xyyzz_xyy, tg_xyyzz_xyz, tg_xyyzz_xzz, tg_xyyzz_yyy, tg_xyyzz_yyz, tg_xyyzz_yzz, tg_xyyzz_zzz, tg_xyz_yyz, tg_xyz_yzz, tg_xyzz_yyy, tg_xyzz_yyz, tg_xyzz_yzz, tg_xyzzz_xxx, tg_xyzzz_xxy, tg_xyzzz_xxz, tg_xyzzz_xyy, tg_xyzzz_xyz, tg_xyzzz_xzz, tg_xyzzz_yyy, tg_xyzzz_yyz, tg_xyzzz_yzz, tg_xyzzz_zzz, tg_xzz_xxz, tg_xzz_xyz, tg_xzz_xzz, tg_xzz_yyy, tg_xzz_yyz, tg_xzz_yzz, tg_xzz_zzz, tg_xzzz_xxx, tg_xzzz_xxz, tg_xzzz_xyz, tg_xzzz_xz, tg_xzzz_xzz, tg_xzzz_yyy, tg_xzzz_yyz, tg_xzzz_yz, tg_xzzz_yzz, tg_xzzz_zz, tg_xzzz_zzz, tg_xzzzz_xxx, tg_xzzzz_xxy, tg_xzzzz_xxz, tg_xzzzz_xyy, tg_xzzzz_xyz, tg_xzzzz_xzz, tg_xzzzz_yyy, tg_xzzzz_yyz, tg_xzzzz_yzz, tg_xzzzz_zzz, tg_yyy_xxx, tg_yyy_xxy, tg_yyy_xxz, tg_yyy_xyy, tg_yyy_xyz, tg_yyy_xzz, tg_yyy_yyy, tg_yyy_yyz, tg_yyy_yzz, tg_yyy_zzz, tg_yyyy_xx, tg_yyyy_xxx, tg_yyyy_xxy, tg_yyyy_xxz, tg_yyyy_xy, tg_yyyy_xyy, tg_yyyy_xyz, tg_yyyy_xz, tg_yyyy_xzz, tg_yyyy_yy, tg_yyyy_yyy, tg_yyyy_yyz, tg_yyyy_yz, tg_yyyy_yzz, tg_yyyy_zz, tg_yyyy_zzz, tg_yyyyy_xxx, tg_yyyyy_xxy, tg_yyyyy_xxz, tg_yyyyy_xyy, tg_yyyyy_xyz, tg_yyyyy_xzz, tg_yyyyy_yyy, tg_yyyyy_yyz, tg_yyyyy_yzz, tg_yyyyy_zzz, tg_yyyyz_xxx, tg_yyyyz_xxy, tg_yyyyz_xxz, tg_yyyyz_xyy, tg_yyyyz_xyz, tg_yyyyz_xzz, tg_yyyyz_yyy, tg_yyyyz_yyz, tg_yyyyz_yzz, tg_yyyyz_zzz, tg_yyyz_xxy, tg_yyyz_xxz, tg_yyyz_xyy, tg_yyyz_xyz, tg_yyyz_xz, tg_yyyz_xzz, tg_yyyz_yyy, tg_yyyz_yyz, tg_yyyz_yz, tg_yyyz_yzz, tg_yyyz_zz, tg_yyyz_zzz, tg_yyyzz_xxx, tg_yyyzz_xxy, tg_yyyzz_xxz, tg_yyyzz_xyy, tg_yyyzz_xyz, tg_yyyzz_xzz, tg_yyyzz_yyy, tg_yyyzz_yyz, tg_yyyzz_yzz, tg_yyyzz_zzz, tg_yyz_xxy, tg_yyz_xxz, tg_yyz_xyy, tg_yyz_xzz, tg_yyz_yyy, tg_yyz_yyz, tg_yyz_yzz, tg_yyz_zzz, tg_yyzz_xx, tg_yyzz_xxx, tg_yyzz_xxy, tg_yyzz_xxz, tg_yyzz_xy, tg_yyzz_xyy, tg_yyzz_xyz, tg_yyzz_xz, tg_yyzz_xzz, tg_yyzz_yy, tg_yyzz_yyy, tg_yyzz_yyz, tg_yyzz_yz, tg_yyzz_yzz, tg_yyzz_zz, tg_yyzz_zzz, tg_yyzzz_xxx, tg_yyzzz_xxy, tg_yyzzz_xxz, tg_yyzzz_xyy, tg_yyzzz_xyz, tg_yyzzz_xzz, tg_yyzzz_yyy, tg_yyzzz_yyz, tg_yyzzz_yzz, tg_yyzzz_zzz, tg_yzz_xxx, tg_yzz_xxz, tg_yzz_xyz, tg_yzz_xzz, tg_yzz_yyy, tg_yzz_yyz, tg_yzz_yzz, tg_yzz_zzz, tg_yzzz_xxx, tg_yzzz_xxy, tg_yzzz_xxz, tg_yzzz_xy, tg_yzzz_xyy, tg_yzzz_xyz, tg_yzzz_xz, tg_yzzz_xzz, tg_yzzz_yy, tg_yzzz_yyy, tg_yzzz_yyz, tg_yzzz_yz, tg_yzzz_yzz, tg_yzzz_zz, tg_yzzz_zzz, tg_yzzzz_xxx, tg_yzzzz_xxy, tg_yzzzz_xxz, tg_yzzzz_xyy, tg_yzzzz_xyz, tg_yzzzz_xzz, tg_yzzzz_yyy, tg_yzzzz_yyz, tg_yzzzz_yzz, tg_yzzzz_zzz, tg_zzz_xxx, tg_zzz_xxy, tg_zzz_xxz, tg_zzz_xyy, tg_zzz_xyz, tg_zzz_xzz, tg_zzz_yyy, tg_zzz_yyz, tg_zzz_yzz, tg_zzz_zzz, tg_zzzz_xx, tg_zzzz_xxx, tg_zzzz_xxy, tg_zzzz_xxz, tg_zzzz_xy, tg_zzzz_xyy, tg_zzzz_xyz, tg_zzzz_xz, tg_zzzz_xzz, tg_zzzz_yy, tg_zzzz_yyy, tg_zzzz_yyz, tg_zzzz_yz, tg_zzzz_yzz, tg_zzzz_zz, tg_zzzz_zzz, tg_zzzzz_xxx, tg_zzzzz_xxy, tg_zzzzz_xxz, tg_zzzzz_xyy, tg_zzzzz_xyz, tg_zzzzz_xzz, tg_zzzzz_yyy, tg_zzzzz_yyz, tg_zzzzz_yzz, tg_zzzzz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_xxx[i] = 4.0 * tg_xxx_xxx[i] * fxi[i] + 3.0 * tg_xxxx_xx[i] * fxi[i] + tg_xxxx_xxx[i] * ra_x[i];

        tg_xxxxx_xxy[i] = 4.0 * tg_xxx_xxy[i] * fxi[i] + 2.0 * tg_xxxx_xy[i] * fxi[i] + tg_xxxx_xxy[i] * ra_x[i];

        tg_xxxxx_xxz[i] = 4.0 * tg_xxx_xxz[i] * fxi[i] + 2.0 * tg_xxxx_xz[i] * fxi[i] + tg_xxxx_xxz[i] * ra_x[i];

        tg_xxxxx_xyy[i] = 4.0 * tg_xxx_xyy[i] * fxi[i] + tg_xxxx_yy[i] * fxi[i] + tg_xxxx_xyy[i] * ra_x[i];

        tg_xxxxx_xyz[i] = 4.0 * tg_xxx_xyz[i] * fxi[i] + tg_xxxx_yz[i] * fxi[i] + tg_xxxx_xyz[i] * ra_x[i];

        tg_xxxxx_xzz[i] = 4.0 * tg_xxx_xzz[i] * fxi[i] + tg_xxxx_zz[i] * fxi[i] + tg_xxxx_xzz[i] * ra_x[i];

        tg_xxxxx_yyy[i] = 4.0 * tg_xxx_yyy[i] * fxi[i] + tg_xxxx_yyy[i] * ra_x[i];

        tg_xxxxx_yyz[i] = 4.0 * tg_xxx_yyz[i] * fxi[i] + tg_xxxx_yyz[i] * ra_x[i];

        tg_xxxxx_yzz[i] = 4.0 * tg_xxx_yzz[i] * fxi[i] + tg_xxxx_yzz[i] * ra_x[i];

        tg_xxxxx_zzz[i] = 4.0 * tg_xxx_zzz[i] * fxi[i] + tg_xxxx_zzz[i] * ra_x[i];

        tg_xxxxy_xxx[i] = tg_xxxx_xxx[i] * ra_y[i];

        tg_xxxxy_xxy[i] = tg_xxxx_xx[i] * fxi[i] + tg_xxxx_xxy[i] * ra_y[i];

        tg_xxxxy_xxz[i] = tg_xxxx_xxz[i] * ra_y[i];

        tg_xxxxy_xyy[i] = 2.0 * tg_xxxx_xy[i] * fxi[i] + tg_xxxx_xyy[i] * ra_y[i];

        tg_xxxxy_xyz[i] = tg_xxxx_xz[i] * fxi[i] + tg_xxxx_xyz[i] * ra_y[i];

        tg_xxxxy_xzz[i] = tg_xxxx_xzz[i] * ra_y[i];

        tg_xxxxy_yyy[i] = 3.0 * tg_xxy_yyy[i] * fxi[i] + tg_xxxy_yyy[i] * ra_x[i];

        tg_xxxxy_yyz[i] = 3.0 * tg_xxy_yyz[i] * fxi[i] + tg_xxxy_yyz[i] * ra_x[i];

        tg_xxxxy_yzz[i] = 3.0 * tg_xxy_yzz[i] * fxi[i] + tg_xxxy_yzz[i] * ra_x[i];

        tg_xxxxy_zzz[i] = tg_xxxx_zzz[i] * ra_y[i];

        tg_xxxxz_xxx[i] = tg_xxxx_xxx[i] * ra_z[i];

        tg_xxxxz_xxy[i] = tg_xxxx_xxy[i] * ra_z[i];

        tg_xxxxz_xxz[i] = tg_xxxx_xx[i] * fxi[i] + tg_xxxx_xxz[i] * ra_z[i];

        tg_xxxxz_xyy[i] = tg_xxxx_xyy[i] * ra_z[i];

        tg_xxxxz_xyz[i] = tg_xxxx_xy[i] * fxi[i] + tg_xxxx_xyz[i] * ra_z[i];

        tg_xxxxz_xzz[i] = 2.0 * tg_xxxx_xz[i] * fxi[i] + tg_xxxx_xzz[i] * ra_z[i];

        tg_xxxxz_yyy[i] = tg_xxxx_yyy[i] * ra_z[i];

        tg_xxxxz_yyz[i] = 3.0 * tg_xxz_yyz[i] * fxi[i] + tg_xxxz_yyz[i] * ra_x[i];

        tg_xxxxz_yzz[i] = 3.0 * tg_xxz_yzz[i] * fxi[i] + tg_xxxz_yzz[i] * ra_x[i];

        tg_xxxxz_zzz[i] = 3.0 * tg_xxz_zzz[i] * fxi[i] + tg_xxxz_zzz[i] * ra_x[i];

        tg_xxxyy_xxx[i] = tg_xxx_xxx[i] * fxi[i] + tg_xxxy_xxx[i] * ra_y[i];

        tg_xxxyy_xxy[i] = 2.0 * tg_xyy_xxy[i] * fxi[i] + 2.0 * tg_xxyy_xy[i] * fxi[i] + tg_xxyy_xxy[i] * ra_x[i];

        tg_xxxyy_xxz[i] = tg_xxx_xxz[i] * fxi[i] + tg_xxxy_xxz[i] * ra_y[i];

        tg_xxxyy_xyy[i] = 2.0 * tg_xyy_xyy[i] * fxi[i] + tg_xxyy_yy[i] * fxi[i] + tg_xxyy_xyy[i] * ra_x[i];

        tg_xxxyy_xyz[i] = 2.0 * tg_xyy_xyz[i] * fxi[i] + tg_xxyy_yz[i] * fxi[i] + tg_xxyy_xyz[i] * ra_x[i];

        tg_xxxyy_xzz[i] = tg_xxx_xzz[i] * fxi[i] + tg_xxxy_xzz[i] * ra_y[i];

        tg_xxxyy_yyy[i] = 2.0 * tg_xyy_yyy[i] * fxi[i] + tg_xxyy_yyy[i] * ra_x[i];

        tg_xxxyy_yyz[i] = 2.0 * tg_xyy_yyz[i] * fxi[i] + tg_xxyy_yyz[i] * ra_x[i];

        tg_xxxyy_yzz[i] = 2.0 * tg_xyy_yzz[i] * fxi[i] + tg_xxyy_yzz[i] * ra_x[i];

        tg_xxxyy_zzz[i] = 2.0 * tg_xyy_zzz[i] * fxi[i] + tg_xxyy_zzz[i] * ra_x[i];

        tg_xxxyz_xxx[i] = tg_xxxz_xxx[i] * ra_y[i];

        tg_xxxyz_xxy[i] = tg_xxxy_xxy[i] * ra_z[i];

        tg_xxxyz_xxz[i] = tg_xxxz_xxz[i] * ra_y[i];

        tg_xxxyz_xyy[i] = tg_xxxy_xyy[i] * ra_z[i];

        tg_xxxyz_xyz[i] = tg_xxxz_xz[i] * fxi[i] + tg_xxxz_xyz[i] * ra_y[i];

        tg_xxxyz_xzz[i] = tg_xxxz_xzz[i] * ra_y[i];

        tg_xxxyz_yyy[i] = tg_xxxy_yyy[i] * ra_z[i];

        tg_xxxyz_yyz[i] = 2.0 * tg_xyz_yyz[i] * fxi[i] + tg_xxyz_yyz[i] * ra_x[i];

        tg_xxxyz_yzz[i] = 2.0 * tg_xyz_yzz[i] * fxi[i] + tg_xxyz_yzz[i] * ra_x[i];

        tg_xxxyz_zzz[i] = tg_xxxz_zzz[i] * ra_y[i];

        tg_xxxzz_xxx[i] = tg_xxx_xxx[i] * fxi[i] + tg_xxxz_xxx[i] * ra_z[i];

        tg_xxxzz_xxy[i] = tg_xxx_xxy[i] * fxi[i] + tg_xxxz_xxy[i] * ra_z[i];

        tg_xxxzz_xxz[i] = 2.0 * tg_xzz_xxz[i] * fxi[i] + 2.0 * tg_xxzz_xz[i] * fxi[i] + tg_xxzz_xxz[i] * ra_x[i];

        tg_xxxzz_xyy[i] = tg_xxx_xyy[i] * fxi[i] + tg_xxxz_xyy[i] * ra_z[i];

        tg_xxxzz_xyz[i] = 2.0 * tg_xzz_xyz[i] * fxi[i] + tg_xxzz_yz[i] * fxi[i] + tg_xxzz_xyz[i] * ra_x[i];

        tg_xxxzz_xzz[i] = 2.0 * tg_xzz_xzz[i] * fxi[i] + tg_xxzz_zz[i] * fxi[i] + tg_xxzz_xzz[i] * ra_x[i];

        tg_xxxzz_yyy[i] = 2.0 * tg_xzz_yyy[i] * fxi[i] + tg_xxzz_yyy[i] * ra_x[i];

        tg_xxxzz_yyz[i] = 2.0 * tg_xzz_yyz[i] * fxi[i] + tg_xxzz_yyz[i] * ra_x[i];

        tg_xxxzz_yzz[i] = 2.0 * tg_xzz_yzz[i] * fxi[i] + tg_xxzz_yzz[i] * ra_x[i];

        tg_xxxzz_zzz[i] = 2.0 * tg_xzz_zzz[i] * fxi[i] + tg_xxzz_zzz[i] * ra_x[i];

        tg_xxyyy_xxx[i] = 2.0 * tg_xxy_xxx[i] * fxi[i] + tg_xxyy_xxx[i] * ra_y[i];

        tg_xxyyy_xxy[i] = tg_yyy_xxy[i] * fxi[i] + 2.0 * tg_xyyy_xy[i] * fxi[i] + tg_xyyy_xxy[i] * ra_x[i];

        tg_xxyyy_xxz[i] = 2.0 * tg_xxy_xxz[i] * fxi[i] + tg_xxyy_xxz[i] * ra_y[i];

        tg_xxyyy_xyy[i] = tg_yyy_xyy[i] * fxi[i] + tg_xyyy_yy[i] * fxi[i] + tg_xyyy_xyy[i] * ra_x[i];

        tg_xxyyy_xyz[i] = tg_yyy_xyz[i] * fxi[i] + tg_xyyy_yz[i] * fxi[i] + tg_xyyy_xyz[i] * ra_x[i];

        tg_xxyyy_xzz[i] = 2.0 * tg_xxy_xzz[i] * fxi[i] + tg_xxyy_xzz[i] * ra_y[i];

        tg_xxyyy_yyy[i] = tg_yyy_yyy[i] * fxi[i] + tg_xyyy_yyy[i] * ra_x[i];

        tg_xxyyy_yyz[i] = tg_yyy_yyz[i] * fxi[i] + tg_xyyy_yyz[i] * ra_x[i];

        tg_xxyyy_yzz[i] = tg_yyy_yzz[i] * fxi[i] + tg_xyyy_yzz[i] * ra_x[i];

        tg_xxyyy_zzz[i] = tg_yyy_zzz[i] * fxi[i] + tg_xyyy_zzz[i] * ra_x[i];

        tg_xxyyz_xxx[i] = tg_xxyy_xxx[i] * ra_z[i];

        tg_xxyyz_xxy[i] = tg_xxyy_xxy[i] * ra_z[i];

        tg_xxyyz_xxz[i] = tg_xxz_xxz[i] * fxi[i] + tg_xxyz_xxz[i] * ra_y[i];

        tg_xxyyz_xyy[i] = tg_xxyy_xyy[i] * ra_z[i];

        tg_xxyyz_xyz[i] = tg_xxyy_xy[i] * fxi[i] + tg_xxyy_xyz[i] * ra_z[i];

        tg_xxyyz_xzz[i] = tg_xxz_xzz[i] * fxi[i] + tg_xxyz_xzz[i] * ra_y[i];

        tg_xxyyz_yyy[i] = tg_xxyy_yyy[i] * ra_z[i];

        tg_xxyyz_yyz[i] = tg_yyz_yyz[i] * fxi[i] + tg_xyyz_yyz[i] * ra_x[i];

        tg_xxyyz_yzz[i] = tg_yyz_yzz[i] * fxi[i] + tg_xyyz_yzz[i] * ra_x[i];

        tg_xxyyz_zzz[i] = tg_yyz_zzz[i] * fxi[i] + tg_xyyz_zzz[i] * ra_x[i];

        tg_xxyzz_xxx[i] = tg_xxzz_xxx[i] * ra_y[i];

        tg_xxyzz_xxy[i] = tg_xxzz_xx[i] * fxi[i] + tg_xxzz_xxy[i] * ra_y[i];

        tg_xxyzz_xxz[i] = tg_xxzz_xxz[i] * ra_y[i];

        tg_xxyzz_xyy[i] = 2.0 * tg_xxzz_xy[i] * fxi[i] + tg_xxzz_xyy[i] * ra_y[i];

        tg_xxyzz_xyz[i] = tg_xxzz_xz[i] * fxi[i] + tg_xxzz_xyz[i] * ra_y[i];

        tg_xxyzz_xzz[i] = tg_xxzz_xzz[i] * ra_y[i];

        tg_xxyzz_yyy[i] = tg_yzz_yyy[i] * fxi[i] + tg_xyzz_yyy[i] * ra_x[i];

        tg_xxyzz_yyz[i] = tg_yzz_yyz[i] * fxi[i] + tg_xyzz_yyz[i] * ra_x[i];

        tg_xxyzz_yzz[i] = tg_yzz_yzz[i] * fxi[i] + tg_xyzz_yzz[i] * ra_x[i];

        tg_xxyzz_zzz[i] = tg_xxzz_zzz[i] * ra_y[i];

        tg_xxzzz_xxx[i] = 2.0 * tg_xxz_xxx[i] * fxi[i] + tg_xxzz_xxx[i] * ra_z[i];

        tg_xxzzz_xxy[i] = 2.0 * tg_xxz_xxy[i] * fxi[i] + tg_xxzz_xxy[i] * ra_z[i];

        tg_xxzzz_xxz[i] = tg_zzz_xxz[i] * fxi[i] + 2.0 * tg_xzzz_xz[i] * fxi[i] + tg_xzzz_xxz[i] * ra_x[i];

        tg_xxzzz_xyy[i] = 2.0 * tg_xxz_xyy[i] * fxi[i] + tg_xxzz_xyy[i] * ra_z[i];

        tg_xxzzz_xyz[i] = tg_zzz_xyz[i] * fxi[i] + tg_xzzz_yz[i] * fxi[i] + tg_xzzz_xyz[i] * ra_x[i];

        tg_xxzzz_xzz[i] = tg_zzz_xzz[i] * fxi[i] + tg_xzzz_zz[i] * fxi[i] + tg_xzzz_xzz[i] * ra_x[i];

        tg_xxzzz_yyy[i] = tg_zzz_yyy[i] * fxi[i] + tg_xzzz_yyy[i] * ra_x[i];

        tg_xxzzz_yyz[i] = tg_zzz_yyz[i] * fxi[i] + tg_xzzz_yyz[i] * ra_x[i];

        tg_xxzzz_yzz[i] = tg_zzz_yzz[i] * fxi[i] + tg_xzzz_yzz[i] * ra_x[i];

        tg_xxzzz_zzz[i] = tg_zzz_zzz[i] * fxi[i] + tg_xzzz_zzz[i] * ra_x[i];

        tg_xyyyy_xxx[i] = 3.0 * tg_yyyy_xx[i] * fxi[i] + tg_yyyy_xxx[i] * ra_x[i];

        tg_xyyyy_xxy[i] = 2.0 * tg_yyyy_xy[i] * fxi[i] + tg_yyyy_xxy[i] * ra_x[i];

        tg_xyyyy_xxz[i] = 2.0 * tg_yyyy_xz[i] * fxi[i] + tg_yyyy_xxz[i] * ra_x[i];

        tg_xyyyy_xyy[i] = tg_yyyy_yy[i] * fxi[i] + tg_yyyy_xyy[i] * ra_x[i];

        tg_xyyyy_xyz[i] = tg_yyyy_yz[i] * fxi[i] + tg_yyyy_xyz[i] * ra_x[i];

        tg_xyyyy_xzz[i] = tg_yyyy_zz[i] * fxi[i] + tg_yyyy_xzz[i] * ra_x[i];

        tg_xyyyy_yyy[i] = tg_yyyy_yyy[i] * ra_x[i];

        tg_xyyyy_yyz[i] = tg_yyyy_yyz[i] * ra_x[i];

        tg_xyyyy_yzz[i] = tg_yyyy_yzz[i] * ra_x[i];

        tg_xyyyy_zzz[i] = tg_yyyy_zzz[i] * ra_x[i];

        tg_xyyyz_xxx[i] = tg_xyyy_xxx[i] * ra_z[i];

        tg_xyyyz_xxy[i] = tg_xyyy_xxy[i] * ra_z[i];

        tg_xyyyz_xxz[i] = 2.0 * tg_yyyz_xz[i] * fxi[i] + tg_yyyz_xxz[i] * ra_x[i];

        tg_xyyyz_xyy[i] = tg_xyyy_xyy[i] * ra_z[i];

        tg_xyyyz_xyz[i] = tg_yyyz_yz[i] * fxi[i] + tg_yyyz_xyz[i] * ra_x[i];

        tg_xyyyz_xzz[i] = tg_yyyz_zz[i] * fxi[i] + tg_yyyz_xzz[i] * ra_x[i];

        tg_xyyyz_yyy[i] = tg_yyyz_yyy[i] * ra_x[i];

        tg_xyyyz_yyz[i] = tg_yyyz_yyz[i] * ra_x[i];

        tg_xyyyz_yzz[i] = tg_yyyz_yzz[i] * ra_x[i];

        tg_xyyyz_zzz[i] = tg_yyyz_zzz[i] * ra_x[i];

        tg_xyyzz_xxx[i] = 3.0 * tg_yyzz_xx[i] * fxi[i] + tg_yyzz_xxx[i] * ra_x[i];

        tg_xyyzz_xxy[i] = 2.0 * tg_yyzz_xy[i] * fxi[i] + tg_yyzz_xxy[i] * ra_x[i];

        tg_xyyzz_xxz[i] = 2.0 * tg_yyzz_xz[i] * fxi[i] + tg_yyzz_xxz[i] * ra_x[i];

        tg_xyyzz_xyy[i] = tg_yyzz_yy[i] * fxi[i] + tg_yyzz_xyy[i] * ra_x[i];

        tg_xyyzz_xyz[i] = tg_yyzz_yz[i] * fxi[i] + tg_yyzz_xyz[i] * ra_x[i];

        tg_xyyzz_xzz[i] = tg_yyzz_zz[i] * fxi[i] + tg_yyzz_xzz[i] * ra_x[i];

        tg_xyyzz_yyy[i] = tg_yyzz_yyy[i] * ra_x[i];

        tg_xyyzz_yyz[i] = tg_yyzz_yyz[i] * ra_x[i];

        tg_xyyzz_yzz[i] = tg_yyzz_yzz[i] * ra_x[i];

        tg_xyyzz_zzz[i] = tg_yyzz_zzz[i] * ra_x[i];

        tg_xyzzz_xxx[i] = tg_xzzz_xxx[i] * ra_y[i];

        tg_xyzzz_xxy[i] = 2.0 * tg_yzzz_xy[i] * fxi[i] + tg_yzzz_xxy[i] * ra_x[i];

        tg_xyzzz_xxz[i] = tg_xzzz_xxz[i] * ra_y[i];

        tg_xyzzz_xyy[i] = tg_yzzz_yy[i] * fxi[i] + tg_yzzz_xyy[i] * ra_x[i];

        tg_xyzzz_xyz[i] = tg_yzzz_yz[i] * fxi[i] + tg_yzzz_xyz[i] * ra_x[i];

        tg_xyzzz_xzz[i] = tg_xzzz_xzz[i] * ra_y[i];

        tg_xyzzz_yyy[i] = tg_yzzz_yyy[i] * ra_x[i];

        tg_xyzzz_yyz[i] = tg_yzzz_yyz[i] * ra_x[i];

        tg_xyzzz_yzz[i] = tg_yzzz_yzz[i] * ra_x[i];

        tg_xyzzz_zzz[i] = tg_yzzz_zzz[i] * ra_x[i];

        tg_xzzzz_xxx[i] = 3.0 * tg_zzzz_xx[i] * fxi[i] + tg_zzzz_xxx[i] * ra_x[i];

        tg_xzzzz_xxy[i] = 2.0 * tg_zzzz_xy[i] * fxi[i] + tg_zzzz_xxy[i] * ra_x[i];

        tg_xzzzz_xxz[i] = 2.0 * tg_zzzz_xz[i] * fxi[i] + tg_zzzz_xxz[i] * ra_x[i];

        tg_xzzzz_xyy[i] = tg_zzzz_yy[i] * fxi[i] + tg_zzzz_xyy[i] * ra_x[i];

        tg_xzzzz_xyz[i] = tg_zzzz_yz[i] * fxi[i] + tg_zzzz_xyz[i] * ra_x[i];

        tg_xzzzz_xzz[i] = tg_zzzz_zz[i] * fxi[i] + tg_zzzz_xzz[i] * ra_x[i];

        tg_xzzzz_yyy[i] = tg_zzzz_yyy[i] * ra_x[i];

        tg_xzzzz_yyz[i] = tg_zzzz_yyz[i] * ra_x[i];

        tg_xzzzz_yzz[i] = tg_zzzz_yzz[i] * ra_x[i];

        tg_xzzzz_zzz[i] = tg_zzzz_zzz[i] * ra_x[i];

        tg_yyyyy_xxx[i] = 4.0 * tg_yyy_xxx[i] * fxi[i] + tg_yyyy_xxx[i] * ra_y[i];

        tg_yyyyy_xxy[i] = 4.0 * tg_yyy_xxy[i] * fxi[i] + tg_yyyy_xx[i] * fxi[i] + tg_yyyy_xxy[i] * ra_y[i];

        tg_yyyyy_xxz[i] = 4.0 * tg_yyy_xxz[i] * fxi[i] + tg_yyyy_xxz[i] * ra_y[i];

        tg_yyyyy_xyy[i] = 4.0 * tg_yyy_xyy[i] * fxi[i] + 2.0 * tg_yyyy_xy[i] * fxi[i] + tg_yyyy_xyy[i] * ra_y[i];

        tg_yyyyy_xyz[i] = 4.0 * tg_yyy_xyz[i] * fxi[i] + tg_yyyy_xz[i] * fxi[i] + tg_yyyy_xyz[i] * ra_y[i];

        tg_yyyyy_xzz[i] = 4.0 * tg_yyy_xzz[i] * fxi[i] + tg_yyyy_xzz[i] * ra_y[i];

        tg_yyyyy_yyy[i] = 4.0 * tg_yyy_yyy[i] * fxi[i] + 3.0 * tg_yyyy_yy[i] * fxi[i] + tg_yyyy_yyy[i] * ra_y[i];

        tg_yyyyy_yyz[i] = 4.0 * tg_yyy_yyz[i] * fxi[i] + 2.0 * tg_yyyy_yz[i] * fxi[i] + tg_yyyy_yyz[i] * ra_y[i];

        tg_yyyyy_yzz[i] = 4.0 * tg_yyy_yzz[i] * fxi[i] + tg_yyyy_zz[i] * fxi[i] + tg_yyyy_yzz[i] * ra_y[i];

        tg_yyyyy_zzz[i] = 4.0 * tg_yyy_zzz[i] * fxi[i] + tg_yyyy_zzz[i] * ra_y[i];

        tg_yyyyz_xxx[i] = tg_yyyy_xxx[i] * ra_z[i];

        tg_yyyyz_xxy[i] = tg_yyyy_xxy[i] * ra_z[i];

        tg_yyyyz_xxz[i] = 3.0 * tg_yyz_xxz[i] * fxi[i] + tg_yyyz_xxz[i] * ra_y[i];

        tg_yyyyz_xyy[i] = tg_yyyy_xyy[i] * ra_z[i];

        tg_yyyyz_xyz[i] = tg_yyyy_xy[i] * fxi[i] + tg_yyyy_xyz[i] * ra_z[i];

        tg_yyyyz_xzz[i] = 3.0 * tg_yyz_xzz[i] * fxi[i] + tg_yyyz_xzz[i] * ra_y[i];

        tg_yyyyz_yyy[i] = tg_yyyy_yyy[i] * ra_z[i];

        tg_yyyyz_yyz[i] = tg_yyyy_yy[i] * fxi[i] + tg_yyyy_yyz[i] * ra_z[i];

        tg_yyyyz_yzz[i] = 2.0 * tg_yyyy_yz[i] * fxi[i] + tg_yyyy_yzz[i] * ra_z[i];

        tg_yyyyz_zzz[i] = 3.0 * tg_yyz_zzz[i] * fxi[i] + tg_yyyz_zzz[i] * ra_y[i];

        tg_yyyzz_xxx[i] = 2.0 * tg_yzz_xxx[i] * fxi[i] + tg_yyzz_xxx[i] * ra_y[i];

        tg_yyyzz_xxy[i] = tg_yyy_xxy[i] * fxi[i] + tg_yyyz_xxy[i] * ra_z[i];

        tg_yyyzz_xxz[i] = 2.0 * tg_yzz_xxz[i] * fxi[i] + tg_yyzz_xxz[i] * ra_y[i];

        tg_yyyzz_xyy[i] = tg_yyy_xyy[i] * fxi[i] + tg_yyyz_xyy[i] * ra_z[i];

        tg_yyyzz_xyz[i] = 2.0 * tg_yzz_xyz[i] * fxi[i] + tg_yyzz_xz[i] * fxi[i] + tg_yyzz_xyz[i] * ra_y[i];

        tg_yyyzz_xzz[i] = 2.0 * tg_yzz_xzz[i] * fxi[i] + tg_yyzz_xzz[i] * ra_y[i];

        tg_yyyzz_yyy[i] = tg_yyy_yyy[i] * fxi[i] + tg_yyyz_yyy[i] * ra_z[i];

        tg_yyyzz_yyz[i] = 2.0 * tg_yzz_yyz[i] * fxi[i] + 2.0 * tg_yyzz_yz[i] * fxi[i] + tg_yyzz_yyz[i] * ra_y[i];

        tg_yyyzz_yzz[i] = 2.0 * tg_yzz_yzz[i] * fxi[i] + tg_yyzz_zz[i] * fxi[i] + tg_yyzz_yzz[i] * ra_y[i];

        tg_yyyzz_zzz[i] = 2.0 * tg_yzz_zzz[i] * fxi[i] + tg_yyzz_zzz[i] * ra_y[i];

        tg_yyzzz_xxx[i] = tg_zzz_xxx[i] * fxi[i] + tg_yzzz_xxx[i] * ra_y[i];

        tg_yyzzz_xxy[i] = 2.0 * tg_yyz_xxy[i] * fxi[i] + tg_yyzz_xxy[i] * ra_z[i];

        tg_yyzzz_xxz[i] = tg_zzz_xxz[i] * fxi[i] + tg_yzzz_xxz[i] * ra_y[i];

        tg_yyzzz_xyy[i] = 2.0 * tg_yyz_xyy[i] * fxi[i] + tg_yyzz_xyy[i] * ra_z[i];

        tg_yyzzz_xyz[i] = tg_zzz_xyz[i] * fxi[i] + tg_yzzz_xz[i] * fxi[i] + tg_yzzz_xyz[i] * ra_y[i];

        tg_yyzzz_xzz[i] = tg_zzz_xzz[i] * fxi[i] + tg_yzzz_xzz[i] * ra_y[i];

        tg_yyzzz_yyy[i] = 2.0 * tg_yyz_yyy[i] * fxi[i] + tg_yyzz_yyy[i] * ra_z[i];

        tg_yyzzz_yyz[i] = tg_zzz_yyz[i] * fxi[i] + 2.0 * tg_yzzz_yz[i] * fxi[i] + tg_yzzz_yyz[i] * ra_y[i];

        tg_yyzzz_yzz[i] = tg_zzz_yzz[i] * fxi[i] + tg_yzzz_zz[i] * fxi[i] + tg_yzzz_yzz[i] * ra_y[i];

        tg_yyzzz_zzz[i] = tg_zzz_zzz[i] * fxi[i] + tg_yzzz_zzz[i] * ra_y[i];

        tg_yzzzz_xxx[i] = tg_zzzz_xxx[i] * ra_y[i];

        tg_yzzzz_xxy[i] = tg_zzzz_xx[i] * fxi[i] + tg_zzzz_xxy[i] * ra_y[i];

        tg_yzzzz_xxz[i] = tg_zzzz_xxz[i] * ra_y[i];

        tg_yzzzz_xyy[i] = 2.0 * tg_zzzz_xy[i] * fxi[i] + tg_zzzz_xyy[i] * ra_y[i];

        tg_yzzzz_xyz[i] = tg_zzzz_xz[i] * fxi[i] + tg_zzzz_xyz[i] * ra_y[i];

        tg_yzzzz_xzz[i] = tg_zzzz_xzz[i] * ra_y[i];

        tg_yzzzz_yyy[i] = 3.0 * tg_zzzz_yy[i] * fxi[i] + tg_zzzz_yyy[i] * ra_y[i];

        tg_yzzzz_yyz[i] = 2.0 * tg_zzzz_yz[i] * fxi[i] + tg_zzzz_yyz[i] * ra_y[i];

        tg_yzzzz_yzz[i] = tg_zzzz_zz[i] * fxi[i] + tg_zzzz_yzz[i] * ra_y[i];

        tg_yzzzz_zzz[i] = tg_zzzz_zzz[i] * ra_y[i];

        tg_zzzzz_xxx[i] = 4.0 * tg_zzz_xxx[i] * fxi[i] + tg_zzzz_xxx[i] * ra_z[i];

        tg_zzzzz_xxy[i] = 4.0 * tg_zzz_xxy[i] * fxi[i] + tg_zzzz_xxy[i] * ra_z[i];

        tg_zzzzz_xxz[i] = 4.0 * tg_zzz_xxz[i] * fxi[i] + tg_zzzz_xx[i] * fxi[i] + tg_zzzz_xxz[i] * ra_z[i];

        tg_zzzzz_xyy[i] = 4.0 * tg_zzz_xyy[i] * fxi[i] + tg_zzzz_xyy[i] * ra_z[i];

        tg_zzzzz_xyz[i] = 4.0 * tg_zzz_xyz[i] * fxi[i] + tg_zzzz_xy[i] * fxi[i] + tg_zzzz_xyz[i] * ra_z[i];

        tg_zzzzz_xzz[i] = 4.0 * tg_zzz_xzz[i] * fxi[i] + 2.0 * tg_zzzz_xz[i] * fxi[i] + tg_zzzz_xzz[i] * ra_z[i];

        tg_zzzzz_yyy[i] = 4.0 * tg_zzz_yyy[i] * fxi[i] + tg_zzzz_yyy[i] * ra_z[i];

        tg_zzzzz_yyz[i] = 4.0 * tg_zzz_yyz[i] * fxi[i] + tg_zzzz_yy[i] * fxi[i] + tg_zzzz_yyz[i] * ra_z[i];

        tg_zzzzz_yzz[i] = 4.0 * tg_zzz_yzz[i] * fxi[i] + 2.0 * tg_zzzz_yz[i] * fxi[i] + tg_zzzz_yzz[i] * ra_z[i];

        tg_zzzzz_zzz[i] = 4.0 * tg_zzz_zzz[i] * fxi[i] + 3.0 * tg_zzzz_zz[i] * fxi[i] + tg_zzzz_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

