#include "LocalCorePotentialPrimRecIF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_if(CSimdArray<double>& pbuffer, 
                                  const size_t idx_if,
                                  const size_t idx_gf,
                                  const size_t idx_hd,
                                  const size_t idx_hf,
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

    auto tg_xxxy_xxz = pbuffer.data(idx_gf + 12);

    auto tg_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto tg_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto tg_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto tg_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto tg_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto tg_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto tg_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto tg_xxxz_xyy = pbuffer.data(idx_gf + 23);

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

    auto tg_yzzz_xxz = pbuffer.data(idx_gf + 132);

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

    // Set up components of auxiliary buffer : HD

    auto tg_xxxxx_xx = pbuffer.data(idx_hd);

    auto tg_xxxxx_xy = pbuffer.data(idx_hd + 1);

    auto tg_xxxxx_xz = pbuffer.data(idx_hd + 2);

    auto tg_xxxxx_yy = pbuffer.data(idx_hd + 3);

    auto tg_xxxxx_yz = pbuffer.data(idx_hd + 4);

    auto tg_xxxxx_zz = pbuffer.data(idx_hd + 5);

    auto tg_xxxxz_xz = pbuffer.data(idx_hd + 14);

    auto tg_xxxyy_xy = pbuffer.data(idx_hd + 19);

    auto tg_xxxyy_yy = pbuffer.data(idx_hd + 21);

    auto tg_xxxyy_yz = pbuffer.data(idx_hd + 22);

    auto tg_xxxzz_xx = pbuffer.data(idx_hd + 30);

    auto tg_xxxzz_xy = pbuffer.data(idx_hd + 31);

    auto tg_xxxzz_xz = pbuffer.data(idx_hd + 32);

    auto tg_xxxzz_yz = pbuffer.data(idx_hd + 34);

    auto tg_xxxzz_zz = pbuffer.data(idx_hd + 35);

    auto tg_xxyyy_xy = pbuffer.data(idx_hd + 37);

    auto tg_xxyyy_yy = pbuffer.data(idx_hd + 39);

    auto tg_xxyyy_yz = pbuffer.data(idx_hd + 40);

    auto tg_xxzzz_xx = pbuffer.data(idx_hd + 54);

    auto tg_xxzzz_xy = pbuffer.data(idx_hd + 55);

    auto tg_xxzzz_xz = pbuffer.data(idx_hd + 56);

    auto tg_xxzzz_yz = pbuffer.data(idx_hd + 58);

    auto tg_xxzzz_zz = pbuffer.data(idx_hd + 59);

    auto tg_xyyyy_xy = pbuffer.data(idx_hd + 61);

    auto tg_xyyyy_yy = pbuffer.data(idx_hd + 63);

    auto tg_xyyyy_yz = pbuffer.data(idx_hd + 64);

    auto tg_xyyzz_yz = pbuffer.data(idx_hd + 76);

    auto tg_xzzzz_xz = pbuffer.data(idx_hd + 86);

    auto tg_xzzzz_yz = pbuffer.data(idx_hd + 88);

    auto tg_xzzzz_zz = pbuffer.data(idx_hd + 89);

    auto tg_yyyyy_xx = pbuffer.data(idx_hd + 90);

    auto tg_yyyyy_xy = pbuffer.data(idx_hd + 91);

    auto tg_yyyyy_xz = pbuffer.data(idx_hd + 92);

    auto tg_yyyyy_yy = pbuffer.data(idx_hd + 93);

    auto tg_yyyyy_yz = pbuffer.data(idx_hd + 94);

    auto tg_yyyyy_zz = pbuffer.data(idx_hd + 95);

    auto tg_yyyyz_xz = pbuffer.data(idx_hd + 98);

    auto tg_yyyyz_yz = pbuffer.data(idx_hd + 100);

    auto tg_yyyyz_zz = pbuffer.data(idx_hd + 101);

    auto tg_yyyzz_xx = pbuffer.data(idx_hd + 102);

    auto tg_yyyzz_xy = pbuffer.data(idx_hd + 103);

    auto tg_yyyzz_xz = pbuffer.data(idx_hd + 104);

    auto tg_yyyzz_yy = pbuffer.data(idx_hd + 105);

    auto tg_yyyzz_yz = pbuffer.data(idx_hd + 106);

    auto tg_yyyzz_zz = pbuffer.data(idx_hd + 107);

    auto tg_yyzzz_xx = pbuffer.data(idx_hd + 108);

    auto tg_yyzzz_xy = pbuffer.data(idx_hd + 109);

    auto tg_yyzzz_xz = pbuffer.data(idx_hd + 110);

    auto tg_yyzzz_yy = pbuffer.data(idx_hd + 111);

    auto tg_yyzzz_yz = pbuffer.data(idx_hd + 112);

    auto tg_yyzzz_zz = pbuffer.data(idx_hd + 113);

    auto tg_yzzzz_xy = pbuffer.data(idx_hd + 115);

    auto tg_yzzzz_xz = pbuffer.data(idx_hd + 116);

    auto tg_yzzzz_yy = pbuffer.data(idx_hd + 117);

    auto tg_yzzzz_yz = pbuffer.data(idx_hd + 118);

    auto tg_yzzzz_zz = pbuffer.data(idx_hd + 119);

    auto tg_zzzzz_xx = pbuffer.data(idx_hd + 120);

    auto tg_zzzzz_xy = pbuffer.data(idx_hd + 121);

    auto tg_zzzzz_xz = pbuffer.data(idx_hd + 122);

    auto tg_zzzzz_yy = pbuffer.data(idx_hd + 123);

    auto tg_zzzzz_yz = pbuffer.data(idx_hd + 124);

    auto tg_zzzzz_zz = pbuffer.data(idx_hd + 125);

    // Set up components of auxiliary buffer : HF

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

    auto tg_xxxxy_xzz = pbuffer.data(idx_hf + 15);

    auto tg_xxxxy_yyy = pbuffer.data(idx_hf + 16);

    auto tg_xxxxy_yyz = pbuffer.data(idx_hf + 17);

    auto tg_xxxxy_yzz = pbuffer.data(idx_hf + 18);

    auto tg_xxxxz_xxx = pbuffer.data(idx_hf + 20);

    auto tg_xxxxz_xxy = pbuffer.data(idx_hf + 21);

    auto tg_xxxxz_xxz = pbuffer.data(idx_hf + 22);

    auto tg_xxxxz_xyy = pbuffer.data(idx_hf + 23);

    auto tg_xxxxz_xyz = pbuffer.data(idx_hf + 24);

    auto tg_xxxxz_xzz = pbuffer.data(idx_hf + 25);

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

    auto tg_xxxyz_xxz = pbuffer.data(idx_hf + 42);

    auto tg_xxxyz_xzz = pbuffer.data(idx_hf + 45);

    auto tg_xxxyz_yyz = pbuffer.data(idx_hf + 47);

    auto tg_xxxyz_yzz = pbuffer.data(idx_hf + 48);

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

    auto tg_xxyyz_xxy = pbuffer.data(idx_hf + 71);

    auto tg_xxyyz_xxz = pbuffer.data(idx_hf + 72);

    auto tg_xxyyz_xyy = pbuffer.data(idx_hf + 73);

    auto tg_xxyyz_xzz = pbuffer.data(idx_hf + 75);

    auto tg_xxyyz_yyz = pbuffer.data(idx_hf + 77);

    auto tg_xxyyz_yzz = pbuffer.data(idx_hf + 78);

    auto tg_xxyyz_zzz = pbuffer.data(idx_hf + 79);

    auto tg_xxyzz_xxx = pbuffer.data(idx_hf + 80);

    auto tg_xxyzz_xxz = pbuffer.data(idx_hf + 82);

    auto tg_xxyzz_xzz = pbuffer.data(idx_hf + 85);

    auto tg_xxyzz_yyy = pbuffer.data(idx_hf + 86);

    auto tg_xxyzz_yyz = pbuffer.data(idx_hf + 87);

    auto tg_xxyzz_yzz = pbuffer.data(idx_hf + 88);

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

    auto tg_xyyyy_xyy = pbuffer.data(idx_hf + 103);

    auto tg_xyyyy_xyz = pbuffer.data(idx_hf + 104);

    auto tg_xyyyy_yyy = pbuffer.data(idx_hf + 106);

    auto tg_xyyyy_yyz = pbuffer.data(idx_hf + 107);

    auto tg_xyyyy_yzz = pbuffer.data(idx_hf + 108);

    auto tg_xyyyy_zzz = pbuffer.data(idx_hf + 109);

    auto tg_xyyyz_yyz = pbuffer.data(idx_hf + 117);

    auto tg_xyyyz_yzz = pbuffer.data(idx_hf + 118);

    auto tg_xyyyz_zzz = pbuffer.data(idx_hf + 119);

    auto tg_xyyzz_xyz = pbuffer.data(idx_hf + 124);

    auto tg_xyyzz_yyy = pbuffer.data(idx_hf + 126);

    auto tg_xyyzz_yyz = pbuffer.data(idx_hf + 127);

    auto tg_xyyzz_yzz = pbuffer.data(idx_hf + 128);

    auto tg_xyyzz_zzz = pbuffer.data(idx_hf + 129);

    auto tg_xyzzz_yyy = pbuffer.data(idx_hf + 136);

    auto tg_xyzzz_yyz = pbuffer.data(idx_hf + 137);

    auto tg_xyzzz_yzz = pbuffer.data(idx_hf + 138);

    auto tg_xzzzz_xxx = pbuffer.data(idx_hf + 140);

    auto tg_xzzzz_xxz = pbuffer.data(idx_hf + 142);

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

    // Set up components of targeted buffer : IF

    auto tg_xxxxxx_xxx = pbuffer.data(idx_if);

    auto tg_xxxxxx_xxy = pbuffer.data(idx_if + 1);

    auto tg_xxxxxx_xxz = pbuffer.data(idx_if + 2);

    auto tg_xxxxxx_xyy = pbuffer.data(idx_if + 3);

    auto tg_xxxxxx_xyz = pbuffer.data(idx_if + 4);

    auto tg_xxxxxx_xzz = pbuffer.data(idx_if + 5);

    auto tg_xxxxxx_yyy = pbuffer.data(idx_if + 6);

    auto tg_xxxxxx_yyz = pbuffer.data(idx_if + 7);

    auto tg_xxxxxx_yzz = pbuffer.data(idx_if + 8);

    auto tg_xxxxxx_zzz = pbuffer.data(idx_if + 9);

    auto tg_xxxxxy_xxx = pbuffer.data(idx_if + 10);

    auto tg_xxxxxy_xxy = pbuffer.data(idx_if + 11);

    auto tg_xxxxxy_xxz = pbuffer.data(idx_if + 12);

    auto tg_xxxxxy_xyy = pbuffer.data(idx_if + 13);

    auto tg_xxxxxy_xyz = pbuffer.data(idx_if + 14);

    auto tg_xxxxxy_xzz = pbuffer.data(idx_if + 15);

    auto tg_xxxxxy_yyy = pbuffer.data(idx_if + 16);

    auto tg_xxxxxy_yyz = pbuffer.data(idx_if + 17);

    auto tg_xxxxxy_yzz = pbuffer.data(idx_if + 18);

    auto tg_xxxxxy_zzz = pbuffer.data(idx_if + 19);

    auto tg_xxxxxz_xxx = pbuffer.data(idx_if + 20);

    auto tg_xxxxxz_xxy = pbuffer.data(idx_if + 21);

    auto tg_xxxxxz_xxz = pbuffer.data(idx_if + 22);

    auto tg_xxxxxz_xyy = pbuffer.data(idx_if + 23);

    auto tg_xxxxxz_xyz = pbuffer.data(idx_if + 24);

    auto tg_xxxxxz_xzz = pbuffer.data(idx_if + 25);

    auto tg_xxxxxz_yyy = pbuffer.data(idx_if + 26);

    auto tg_xxxxxz_yyz = pbuffer.data(idx_if + 27);

    auto tg_xxxxxz_yzz = pbuffer.data(idx_if + 28);

    auto tg_xxxxxz_zzz = pbuffer.data(idx_if + 29);

    auto tg_xxxxyy_xxx = pbuffer.data(idx_if + 30);

    auto tg_xxxxyy_xxy = pbuffer.data(idx_if + 31);

    auto tg_xxxxyy_xxz = pbuffer.data(idx_if + 32);

    auto tg_xxxxyy_xyy = pbuffer.data(idx_if + 33);

    auto tg_xxxxyy_xyz = pbuffer.data(idx_if + 34);

    auto tg_xxxxyy_xzz = pbuffer.data(idx_if + 35);

    auto tg_xxxxyy_yyy = pbuffer.data(idx_if + 36);

    auto tg_xxxxyy_yyz = pbuffer.data(idx_if + 37);

    auto tg_xxxxyy_yzz = pbuffer.data(idx_if + 38);

    auto tg_xxxxyy_zzz = pbuffer.data(idx_if + 39);

    auto tg_xxxxyz_xxx = pbuffer.data(idx_if + 40);

    auto tg_xxxxyz_xxy = pbuffer.data(idx_if + 41);

    auto tg_xxxxyz_xxz = pbuffer.data(idx_if + 42);

    auto tg_xxxxyz_xyy = pbuffer.data(idx_if + 43);

    auto tg_xxxxyz_xyz = pbuffer.data(idx_if + 44);

    auto tg_xxxxyz_xzz = pbuffer.data(idx_if + 45);

    auto tg_xxxxyz_yyy = pbuffer.data(idx_if + 46);

    auto tg_xxxxyz_yyz = pbuffer.data(idx_if + 47);

    auto tg_xxxxyz_yzz = pbuffer.data(idx_if + 48);

    auto tg_xxxxyz_zzz = pbuffer.data(idx_if + 49);

    auto tg_xxxxzz_xxx = pbuffer.data(idx_if + 50);

    auto tg_xxxxzz_xxy = pbuffer.data(idx_if + 51);

    auto tg_xxxxzz_xxz = pbuffer.data(idx_if + 52);

    auto tg_xxxxzz_xyy = pbuffer.data(idx_if + 53);

    auto tg_xxxxzz_xyz = pbuffer.data(idx_if + 54);

    auto tg_xxxxzz_xzz = pbuffer.data(idx_if + 55);

    auto tg_xxxxzz_yyy = pbuffer.data(idx_if + 56);

    auto tg_xxxxzz_yyz = pbuffer.data(idx_if + 57);

    auto tg_xxxxzz_yzz = pbuffer.data(idx_if + 58);

    auto tg_xxxxzz_zzz = pbuffer.data(idx_if + 59);

    auto tg_xxxyyy_xxx = pbuffer.data(idx_if + 60);

    auto tg_xxxyyy_xxy = pbuffer.data(idx_if + 61);

    auto tg_xxxyyy_xxz = pbuffer.data(idx_if + 62);

    auto tg_xxxyyy_xyy = pbuffer.data(idx_if + 63);

    auto tg_xxxyyy_xyz = pbuffer.data(idx_if + 64);

    auto tg_xxxyyy_xzz = pbuffer.data(idx_if + 65);

    auto tg_xxxyyy_yyy = pbuffer.data(idx_if + 66);

    auto tg_xxxyyy_yyz = pbuffer.data(idx_if + 67);

    auto tg_xxxyyy_yzz = pbuffer.data(idx_if + 68);

    auto tg_xxxyyy_zzz = pbuffer.data(idx_if + 69);

    auto tg_xxxyyz_xxx = pbuffer.data(idx_if + 70);

    auto tg_xxxyyz_xxy = pbuffer.data(idx_if + 71);

    auto tg_xxxyyz_xxz = pbuffer.data(idx_if + 72);

    auto tg_xxxyyz_xyy = pbuffer.data(idx_if + 73);

    auto tg_xxxyyz_xyz = pbuffer.data(idx_if + 74);

    auto tg_xxxyyz_xzz = pbuffer.data(idx_if + 75);

    auto tg_xxxyyz_yyy = pbuffer.data(idx_if + 76);

    auto tg_xxxyyz_yyz = pbuffer.data(idx_if + 77);

    auto tg_xxxyyz_yzz = pbuffer.data(idx_if + 78);

    auto tg_xxxyyz_zzz = pbuffer.data(idx_if + 79);

    auto tg_xxxyzz_xxx = pbuffer.data(idx_if + 80);

    auto tg_xxxyzz_xxy = pbuffer.data(idx_if + 81);

    auto tg_xxxyzz_xxz = pbuffer.data(idx_if + 82);

    auto tg_xxxyzz_xyy = pbuffer.data(idx_if + 83);

    auto tg_xxxyzz_xyz = pbuffer.data(idx_if + 84);

    auto tg_xxxyzz_xzz = pbuffer.data(idx_if + 85);

    auto tg_xxxyzz_yyy = pbuffer.data(idx_if + 86);

    auto tg_xxxyzz_yyz = pbuffer.data(idx_if + 87);

    auto tg_xxxyzz_yzz = pbuffer.data(idx_if + 88);

    auto tg_xxxyzz_zzz = pbuffer.data(idx_if + 89);

    auto tg_xxxzzz_xxx = pbuffer.data(idx_if + 90);

    auto tg_xxxzzz_xxy = pbuffer.data(idx_if + 91);

    auto tg_xxxzzz_xxz = pbuffer.data(idx_if + 92);

    auto tg_xxxzzz_xyy = pbuffer.data(idx_if + 93);

    auto tg_xxxzzz_xyz = pbuffer.data(idx_if + 94);

    auto tg_xxxzzz_xzz = pbuffer.data(idx_if + 95);

    auto tg_xxxzzz_yyy = pbuffer.data(idx_if + 96);

    auto tg_xxxzzz_yyz = pbuffer.data(idx_if + 97);

    auto tg_xxxzzz_yzz = pbuffer.data(idx_if + 98);

    auto tg_xxxzzz_zzz = pbuffer.data(idx_if + 99);

    auto tg_xxyyyy_xxx = pbuffer.data(idx_if + 100);

    auto tg_xxyyyy_xxy = pbuffer.data(idx_if + 101);

    auto tg_xxyyyy_xxz = pbuffer.data(idx_if + 102);

    auto tg_xxyyyy_xyy = pbuffer.data(idx_if + 103);

    auto tg_xxyyyy_xyz = pbuffer.data(idx_if + 104);

    auto tg_xxyyyy_xzz = pbuffer.data(idx_if + 105);

    auto tg_xxyyyy_yyy = pbuffer.data(idx_if + 106);

    auto tg_xxyyyy_yyz = pbuffer.data(idx_if + 107);

    auto tg_xxyyyy_yzz = pbuffer.data(idx_if + 108);

    auto tg_xxyyyy_zzz = pbuffer.data(idx_if + 109);

    auto tg_xxyyyz_xxx = pbuffer.data(idx_if + 110);

    auto tg_xxyyyz_xxy = pbuffer.data(idx_if + 111);

    auto tg_xxyyyz_xxz = pbuffer.data(idx_if + 112);

    auto tg_xxyyyz_xyy = pbuffer.data(idx_if + 113);

    auto tg_xxyyyz_xyz = pbuffer.data(idx_if + 114);

    auto tg_xxyyyz_xzz = pbuffer.data(idx_if + 115);

    auto tg_xxyyyz_yyy = pbuffer.data(idx_if + 116);

    auto tg_xxyyyz_yyz = pbuffer.data(idx_if + 117);

    auto tg_xxyyyz_yzz = pbuffer.data(idx_if + 118);

    auto tg_xxyyyz_zzz = pbuffer.data(idx_if + 119);

    auto tg_xxyyzz_xxx = pbuffer.data(idx_if + 120);

    auto tg_xxyyzz_xxy = pbuffer.data(idx_if + 121);

    auto tg_xxyyzz_xxz = pbuffer.data(idx_if + 122);

    auto tg_xxyyzz_xyy = pbuffer.data(idx_if + 123);

    auto tg_xxyyzz_xyz = pbuffer.data(idx_if + 124);

    auto tg_xxyyzz_xzz = pbuffer.data(idx_if + 125);

    auto tg_xxyyzz_yyy = pbuffer.data(idx_if + 126);

    auto tg_xxyyzz_yyz = pbuffer.data(idx_if + 127);

    auto tg_xxyyzz_yzz = pbuffer.data(idx_if + 128);

    auto tg_xxyyzz_zzz = pbuffer.data(idx_if + 129);

    auto tg_xxyzzz_xxx = pbuffer.data(idx_if + 130);

    auto tg_xxyzzz_xxy = pbuffer.data(idx_if + 131);

    auto tg_xxyzzz_xxz = pbuffer.data(idx_if + 132);

    auto tg_xxyzzz_xyy = pbuffer.data(idx_if + 133);

    auto tg_xxyzzz_xyz = pbuffer.data(idx_if + 134);

    auto tg_xxyzzz_xzz = pbuffer.data(idx_if + 135);

    auto tg_xxyzzz_yyy = pbuffer.data(idx_if + 136);

    auto tg_xxyzzz_yyz = pbuffer.data(idx_if + 137);

    auto tg_xxyzzz_yzz = pbuffer.data(idx_if + 138);

    auto tg_xxyzzz_zzz = pbuffer.data(idx_if + 139);

    auto tg_xxzzzz_xxx = pbuffer.data(idx_if + 140);

    auto tg_xxzzzz_xxy = pbuffer.data(idx_if + 141);

    auto tg_xxzzzz_xxz = pbuffer.data(idx_if + 142);

    auto tg_xxzzzz_xyy = pbuffer.data(idx_if + 143);

    auto tg_xxzzzz_xyz = pbuffer.data(idx_if + 144);

    auto tg_xxzzzz_xzz = pbuffer.data(idx_if + 145);

    auto tg_xxzzzz_yyy = pbuffer.data(idx_if + 146);

    auto tg_xxzzzz_yyz = pbuffer.data(idx_if + 147);

    auto tg_xxzzzz_yzz = pbuffer.data(idx_if + 148);

    auto tg_xxzzzz_zzz = pbuffer.data(idx_if + 149);

    auto tg_xyyyyy_xxx = pbuffer.data(idx_if + 150);

    auto tg_xyyyyy_xxy = pbuffer.data(idx_if + 151);

    auto tg_xyyyyy_xxz = pbuffer.data(idx_if + 152);

    auto tg_xyyyyy_xyy = pbuffer.data(idx_if + 153);

    auto tg_xyyyyy_xyz = pbuffer.data(idx_if + 154);

    auto tg_xyyyyy_xzz = pbuffer.data(idx_if + 155);

    auto tg_xyyyyy_yyy = pbuffer.data(idx_if + 156);

    auto tg_xyyyyy_yyz = pbuffer.data(idx_if + 157);

    auto tg_xyyyyy_yzz = pbuffer.data(idx_if + 158);

    auto tg_xyyyyy_zzz = pbuffer.data(idx_if + 159);

    auto tg_xyyyyz_xxx = pbuffer.data(idx_if + 160);

    auto tg_xyyyyz_xxy = pbuffer.data(idx_if + 161);

    auto tg_xyyyyz_xxz = pbuffer.data(idx_if + 162);

    auto tg_xyyyyz_xyy = pbuffer.data(idx_if + 163);

    auto tg_xyyyyz_xyz = pbuffer.data(idx_if + 164);

    auto tg_xyyyyz_xzz = pbuffer.data(idx_if + 165);

    auto tg_xyyyyz_yyy = pbuffer.data(idx_if + 166);

    auto tg_xyyyyz_yyz = pbuffer.data(idx_if + 167);

    auto tg_xyyyyz_yzz = pbuffer.data(idx_if + 168);

    auto tg_xyyyyz_zzz = pbuffer.data(idx_if + 169);

    auto tg_xyyyzz_xxx = pbuffer.data(idx_if + 170);

    auto tg_xyyyzz_xxy = pbuffer.data(idx_if + 171);

    auto tg_xyyyzz_xxz = pbuffer.data(idx_if + 172);

    auto tg_xyyyzz_xyy = pbuffer.data(idx_if + 173);

    auto tg_xyyyzz_xyz = pbuffer.data(idx_if + 174);

    auto tg_xyyyzz_xzz = pbuffer.data(idx_if + 175);

    auto tg_xyyyzz_yyy = pbuffer.data(idx_if + 176);

    auto tg_xyyyzz_yyz = pbuffer.data(idx_if + 177);

    auto tg_xyyyzz_yzz = pbuffer.data(idx_if + 178);

    auto tg_xyyyzz_zzz = pbuffer.data(idx_if + 179);

    auto tg_xyyzzz_xxx = pbuffer.data(idx_if + 180);

    auto tg_xyyzzz_xxy = pbuffer.data(idx_if + 181);

    auto tg_xyyzzz_xxz = pbuffer.data(idx_if + 182);

    auto tg_xyyzzz_xyy = pbuffer.data(idx_if + 183);

    auto tg_xyyzzz_xyz = pbuffer.data(idx_if + 184);

    auto tg_xyyzzz_xzz = pbuffer.data(idx_if + 185);

    auto tg_xyyzzz_yyy = pbuffer.data(idx_if + 186);

    auto tg_xyyzzz_yyz = pbuffer.data(idx_if + 187);

    auto tg_xyyzzz_yzz = pbuffer.data(idx_if + 188);

    auto tg_xyyzzz_zzz = pbuffer.data(idx_if + 189);

    auto tg_xyzzzz_xxx = pbuffer.data(idx_if + 190);

    auto tg_xyzzzz_xxy = pbuffer.data(idx_if + 191);

    auto tg_xyzzzz_xxz = pbuffer.data(idx_if + 192);

    auto tg_xyzzzz_xyy = pbuffer.data(idx_if + 193);

    auto tg_xyzzzz_xyz = pbuffer.data(idx_if + 194);

    auto tg_xyzzzz_xzz = pbuffer.data(idx_if + 195);

    auto tg_xyzzzz_yyy = pbuffer.data(idx_if + 196);

    auto tg_xyzzzz_yyz = pbuffer.data(idx_if + 197);

    auto tg_xyzzzz_yzz = pbuffer.data(idx_if + 198);

    auto tg_xyzzzz_zzz = pbuffer.data(idx_if + 199);

    auto tg_xzzzzz_xxx = pbuffer.data(idx_if + 200);

    auto tg_xzzzzz_xxy = pbuffer.data(idx_if + 201);

    auto tg_xzzzzz_xxz = pbuffer.data(idx_if + 202);

    auto tg_xzzzzz_xyy = pbuffer.data(idx_if + 203);

    auto tg_xzzzzz_xyz = pbuffer.data(idx_if + 204);

    auto tg_xzzzzz_xzz = pbuffer.data(idx_if + 205);

    auto tg_xzzzzz_yyy = pbuffer.data(idx_if + 206);

    auto tg_xzzzzz_yyz = pbuffer.data(idx_if + 207);

    auto tg_xzzzzz_yzz = pbuffer.data(idx_if + 208);

    auto tg_xzzzzz_zzz = pbuffer.data(idx_if + 209);

    auto tg_yyyyyy_xxx = pbuffer.data(idx_if + 210);

    auto tg_yyyyyy_xxy = pbuffer.data(idx_if + 211);

    auto tg_yyyyyy_xxz = pbuffer.data(idx_if + 212);

    auto tg_yyyyyy_xyy = pbuffer.data(idx_if + 213);

    auto tg_yyyyyy_xyz = pbuffer.data(idx_if + 214);

    auto tg_yyyyyy_xzz = pbuffer.data(idx_if + 215);

    auto tg_yyyyyy_yyy = pbuffer.data(idx_if + 216);

    auto tg_yyyyyy_yyz = pbuffer.data(idx_if + 217);

    auto tg_yyyyyy_yzz = pbuffer.data(idx_if + 218);

    auto tg_yyyyyy_zzz = pbuffer.data(idx_if + 219);

    auto tg_yyyyyz_xxx = pbuffer.data(idx_if + 220);

    auto tg_yyyyyz_xxy = pbuffer.data(idx_if + 221);

    auto tg_yyyyyz_xxz = pbuffer.data(idx_if + 222);

    auto tg_yyyyyz_xyy = pbuffer.data(idx_if + 223);

    auto tg_yyyyyz_xyz = pbuffer.data(idx_if + 224);

    auto tg_yyyyyz_xzz = pbuffer.data(idx_if + 225);

    auto tg_yyyyyz_yyy = pbuffer.data(idx_if + 226);

    auto tg_yyyyyz_yyz = pbuffer.data(idx_if + 227);

    auto tg_yyyyyz_yzz = pbuffer.data(idx_if + 228);

    auto tg_yyyyyz_zzz = pbuffer.data(idx_if + 229);

    auto tg_yyyyzz_xxx = pbuffer.data(idx_if + 230);

    auto tg_yyyyzz_xxy = pbuffer.data(idx_if + 231);

    auto tg_yyyyzz_xxz = pbuffer.data(idx_if + 232);

    auto tg_yyyyzz_xyy = pbuffer.data(idx_if + 233);

    auto tg_yyyyzz_xyz = pbuffer.data(idx_if + 234);

    auto tg_yyyyzz_xzz = pbuffer.data(idx_if + 235);

    auto tg_yyyyzz_yyy = pbuffer.data(idx_if + 236);

    auto tg_yyyyzz_yyz = pbuffer.data(idx_if + 237);

    auto tg_yyyyzz_yzz = pbuffer.data(idx_if + 238);

    auto tg_yyyyzz_zzz = pbuffer.data(idx_if + 239);

    auto tg_yyyzzz_xxx = pbuffer.data(idx_if + 240);

    auto tg_yyyzzz_xxy = pbuffer.data(idx_if + 241);

    auto tg_yyyzzz_xxz = pbuffer.data(idx_if + 242);

    auto tg_yyyzzz_xyy = pbuffer.data(idx_if + 243);

    auto tg_yyyzzz_xyz = pbuffer.data(idx_if + 244);

    auto tg_yyyzzz_xzz = pbuffer.data(idx_if + 245);

    auto tg_yyyzzz_yyy = pbuffer.data(idx_if + 246);

    auto tg_yyyzzz_yyz = pbuffer.data(idx_if + 247);

    auto tg_yyyzzz_yzz = pbuffer.data(idx_if + 248);

    auto tg_yyyzzz_zzz = pbuffer.data(idx_if + 249);

    auto tg_yyzzzz_xxx = pbuffer.data(idx_if + 250);

    auto tg_yyzzzz_xxy = pbuffer.data(idx_if + 251);

    auto tg_yyzzzz_xxz = pbuffer.data(idx_if + 252);

    auto tg_yyzzzz_xyy = pbuffer.data(idx_if + 253);

    auto tg_yyzzzz_xyz = pbuffer.data(idx_if + 254);

    auto tg_yyzzzz_xzz = pbuffer.data(idx_if + 255);

    auto tg_yyzzzz_yyy = pbuffer.data(idx_if + 256);

    auto tg_yyzzzz_yyz = pbuffer.data(idx_if + 257);

    auto tg_yyzzzz_yzz = pbuffer.data(idx_if + 258);

    auto tg_yyzzzz_zzz = pbuffer.data(idx_if + 259);

    auto tg_yzzzzz_xxx = pbuffer.data(idx_if + 260);

    auto tg_yzzzzz_xxy = pbuffer.data(idx_if + 261);

    auto tg_yzzzzz_xxz = pbuffer.data(idx_if + 262);

    auto tg_yzzzzz_xyy = pbuffer.data(idx_if + 263);

    auto tg_yzzzzz_xyz = pbuffer.data(idx_if + 264);

    auto tg_yzzzzz_xzz = pbuffer.data(idx_if + 265);

    auto tg_yzzzzz_yyy = pbuffer.data(idx_if + 266);

    auto tg_yzzzzz_yyz = pbuffer.data(idx_if + 267);

    auto tg_yzzzzz_yzz = pbuffer.data(idx_if + 268);

    auto tg_yzzzzz_zzz = pbuffer.data(idx_if + 269);

    auto tg_zzzzzz_xxx = pbuffer.data(idx_if + 270);

    auto tg_zzzzzz_xxy = pbuffer.data(idx_if + 271);

    auto tg_zzzzzz_xxz = pbuffer.data(idx_if + 272);

    auto tg_zzzzzz_xyy = pbuffer.data(idx_if + 273);

    auto tg_zzzzzz_xyz = pbuffer.data(idx_if + 274);

    auto tg_zzzzzz_xzz = pbuffer.data(idx_if + 275);

    auto tg_zzzzzz_yyy = pbuffer.data(idx_if + 276);

    auto tg_zzzzzz_yyz = pbuffer.data(idx_if + 277);

    auto tg_zzzzzz_yzz = pbuffer.data(idx_if + 278);

    auto tg_zzzzzz_zzz = pbuffer.data(idx_if + 279);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxx_xxx, tg_xxxx_xxy, tg_xxxx_xxz, tg_xxxx_xyy, tg_xxxx_xyz, tg_xxxx_xzz, tg_xxxx_yyy, tg_xxxx_yyz, tg_xxxx_yzz, tg_xxxx_zzz, tg_xxxxx_xx, tg_xxxxx_xxx, tg_xxxxx_xxy, tg_xxxxx_xxz, tg_xxxxx_xy, tg_xxxxx_xyy, tg_xxxxx_xyz, tg_xxxxx_xz, tg_xxxxx_xzz, tg_xxxxx_yy, tg_xxxxx_yyy, tg_xxxxx_yyz, tg_xxxxx_yz, tg_xxxxx_yzz, tg_xxxxx_zz, tg_xxxxx_zzz, tg_xxxxxx_xxx, tg_xxxxxx_xxy, tg_xxxxxx_xxz, tg_xxxxxx_xyy, tg_xxxxxx_xyz, tg_xxxxxx_xzz, tg_xxxxxx_yyy, tg_xxxxxx_yyz, tg_xxxxxx_yzz, tg_xxxxxx_zzz, tg_xxxxxy_xxx, tg_xxxxxy_xxy, tg_xxxxxy_xxz, tg_xxxxxy_xyy, tg_xxxxxy_xyz, tg_xxxxxy_xzz, tg_xxxxxy_yyy, tg_xxxxxy_yyz, tg_xxxxxy_yzz, tg_xxxxxy_zzz, tg_xxxxxz_xxx, tg_xxxxxz_xxy, tg_xxxxxz_xxz, tg_xxxxxz_xyy, tg_xxxxxz_xyz, tg_xxxxxz_xzz, tg_xxxxxz_yyy, tg_xxxxxz_yyz, tg_xxxxxz_yzz, tg_xxxxxz_zzz, tg_xxxxy_xxx, tg_xxxxy_xxy, tg_xxxxy_xxz, tg_xxxxy_xyy, tg_xxxxy_xzz, tg_xxxxy_yyy, tg_xxxxy_yyz, tg_xxxxy_yzz, tg_xxxxyy_xxx, tg_xxxxyy_xxy, tg_xxxxyy_xxz, tg_xxxxyy_xyy, tg_xxxxyy_xyz, tg_xxxxyy_xzz, tg_xxxxyy_yyy, tg_xxxxyy_yyz, tg_xxxxyy_yzz, tg_xxxxyy_zzz, tg_xxxxyz_xxx, tg_xxxxyz_xxy, tg_xxxxyz_xxz, tg_xxxxyz_xyy, tg_xxxxyz_xyz, tg_xxxxyz_xzz, tg_xxxxyz_yyy, tg_xxxxyz_yyz, tg_xxxxyz_yzz, tg_xxxxyz_zzz, tg_xxxxz_xxx, tg_xxxxz_xxy, tg_xxxxz_xxz, tg_xxxxz_xyy, tg_xxxxz_xyz, tg_xxxxz_xz, tg_xxxxz_xzz, tg_xxxxz_yyz, tg_xxxxz_yzz, tg_xxxxz_zzz, tg_xxxxzz_xxx, tg_xxxxzz_xxy, tg_xxxxzz_xxz, tg_xxxxzz_xyy, tg_xxxxzz_xyz, tg_xxxxzz_xzz, tg_xxxxzz_yyy, tg_xxxxzz_yyz, tg_xxxxzz_yzz, tg_xxxxzz_zzz, tg_xxxy_xxx, tg_xxxy_xxz, tg_xxxy_xzz, tg_xxxy_yyy, tg_xxxy_yyz, tg_xxxy_yzz, tg_xxxyy_xxx, tg_xxxyy_xxy, tg_xxxyy_xxz, tg_xxxyy_xy, tg_xxxyy_xyy, tg_xxxyy_xyz, tg_xxxyy_xzz, tg_xxxyy_yy, tg_xxxyy_yyy, tg_xxxyy_yyz, tg_xxxyy_yz, tg_xxxyy_yzz, tg_xxxyy_zzz, tg_xxxyyy_xxx, tg_xxxyyy_xxy, tg_xxxyyy_xxz, tg_xxxyyy_xyy, tg_xxxyyy_xyz, tg_xxxyyy_xzz, tg_xxxyyy_yyy, tg_xxxyyy_yyz, tg_xxxyyy_yzz, tg_xxxyyy_zzz, tg_xxxyyz_xxx, tg_xxxyyz_xxy, tg_xxxyyz_xxz, tg_xxxyyz_xyy, tg_xxxyyz_xyz, tg_xxxyyz_xzz, tg_xxxyyz_yyy, tg_xxxyyz_yyz, tg_xxxyyz_yzz, tg_xxxyyz_zzz, tg_xxxyz_xxz, tg_xxxyz_xzz, tg_xxxyz_yyz, tg_xxxyz_yzz, tg_xxxyzz_xxx, tg_xxxyzz_xxy, tg_xxxyzz_xxz, tg_xxxyzz_xyy, tg_xxxyzz_xyz, tg_xxxyzz_xzz, tg_xxxyzz_yyy, tg_xxxyzz_yyz, tg_xxxyzz_yzz, tg_xxxyzz_zzz, tg_xxxz_xxx, tg_xxxz_xxy, tg_xxxz_xxz, tg_xxxz_xyy, tg_xxxz_xzz, tg_xxxz_yyz, tg_xxxz_yzz, tg_xxxz_zzz, tg_xxxzz_xx, tg_xxxzz_xxx, tg_xxxzz_xxy, tg_xxxzz_xxz, tg_xxxzz_xy, tg_xxxzz_xyy, tg_xxxzz_xyz, tg_xxxzz_xz, tg_xxxzz_xzz, tg_xxxzz_yyy, tg_xxxzz_yyz, tg_xxxzz_yz, tg_xxxzz_yzz, tg_xxxzz_zz, tg_xxxzz_zzz, tg_xxxzzz_xxx, tg_xxxzzz_xxy, tg_xxxzzz_xxz, tg_xxxzzz_xyy, tg_xxxzzz_xyz, tg_xxxzzz_xzz, tg_xxxzzz_yyy, tg_xxxzzz_yyz, tg_xxxzzz_yzz, tg_xxxzzz_zzz, tg_xxyy_xxx, tg_xxyy_xxy, tg_xxyy_xxz, tg_xxyy_xyy, tg_xxyy_xyz, tg_xxyy_xzz, tg_xxyy_yyy, tg_xxyy_yyz, tg_xxyy_yzz, tg_xxyy_zzz, tg_xxyyy_xxx, tg_xxyyy_xxy, tg_xxyyy_xxz, tg_xxyyy_xy, tg_xxyyy_xyy, tg_xxyyy_xyz, tg_xxyyy_xzz, tg_xxyyy_yy, tg_xxyyy_yyy, tg_xxyyy_yyz, tg_xxyyy_yz, tg_xxyyy_yzz, tg_xxyyy_zzz, tg_xxyyyy_xxx, tg_xxyyyy_xxy, tg_xxyyyy_xxz, tg_xxyyyy_xyy, tg_xxyyyy_xyz, tg_xxyyyy_xzz, tg_xxyyyy_yyy, tg_xxyyyy_yyz, tg_xxyyyy_yzz, tg_xxyyyy_zzz, tg_xxyyyz_xxx, tg_xxyyyz_xxy, tg_xxyyyz_xxz, tg_xxyyyz_xyy, tg_xxyyyz_xyz, tg_xxyyyz_xzz, tg_xxyyyz_yyy, tg_xxyyyz_yyz, tg_xxyyyz_yzz, tg_xxyyyz_zzz, tg_xxyyz_xxy, tg_xxyyz_xxz, tg_xxyyz_xyy, tg_xxyyz_xzz, tg_xxyyz_yyz, tg_xxyyz_yzz, tg_xxyyz_zzz, tg_xxyyzz_xxx, tg_xxyyzz_xxy, tg_xxyyzz_xxz, tg_xxyyzz_xyy, tg_xxyyzz_xyz, tg_xxyyzz_xzz, tg_xxyyzz_yyy, tg_xxyyzz_yyz, tg_xxyyzz_yzz, tg_xxyyzz_zzz, tg_xxyz_xxz, tg_xxyz_xzz, tg_xxyz_yyz, tg_xxyz_yzz, tg_xxyzz_xxx, tg_xxyzz_xxz, tg_xxyzz_xzz, tg_xxyzz_yyy, tg_xxyzz_yyz, tg_xxyzz_yzz, tg_xxyzzz_xxx, tg_xxyzzz_xxy, tg_xxyzzz_xxz, tg_xxyzzz_xyy, tg_xxyzzz_xyz, tg_xxyzzz_xzz, tg_xxyzzz_yyy, tg_xxyzzz_yyz, tg_xxyzzz_yzz, tg_xxyzzz_zzz, tg_xxzz_xxx, tg_xxzz_xxy, tg_xxzz_xxz, tg_xxzz_xyy, tg_xxzz_xyz, tg_xxzz_xzz, tg_xxzz_yyy, tg_xxzz_yyz, tg_xxzz_yzz, tg_xxzz_zzz, tg_xxzzz_xx, tg_xxzzz_xxx, tg_xxzzz_xxy, tg_xxzzz_xxz, tg_xxzzz_xy, tg_xxzzz_xyy, tg_xxzzz_xyz, tg_xxzzz_xz, tg_xxzzz_xzz, tg_xxzzz_yyy, tg_xxzzz_yyz, tg_xxzzz_yz, tg_xxzzz_yzz, tg_xxzzz_zz, tg_xxzzz_zzz, tg_xxzzzz_xxx, tg_xxzzzz_xxy, tg_xxzzzz_xxz, tg_xxzzzz_xyy, tg_xxzzzz_xyz, tg_xxzzzz_xzz, tg_xxzzzz_yyy, tg_xxzzzz_yyz, tg_xxzzzz_yzz, tg_xxzzzz_zzz, tg_xyyy_xxy, tg_xyyy_xyy, tg_xyyy_xyz, tg_xyyy_yyy, tg_xyyy_yyz, tg_xyyy_yzz, tg_xyyy_zzz, tg_xyyyy_xxx, tg_xyyyy_xxy, tg_xyyyy_xy, tg_xyyyy_xyy, tg_xyyyy_xyz, tg_xyyyy_yy, tg_xyyyy_yyy, tg_xyyyy_yyz, tg_xyyyy_yz, tg_xyyyy_yzz, tg_xyyyy_zzz, tg_xyyyyy_xxx, tg_xyyyyy_xxy, tg_xyyyyy_xxz, tg_xyyyyy_xyy, tg_xyyyyy_xyz, tg_xyyyyy_xzz, tg_xyyyyy_yyy, tg_xyyyyy_yyz, tg_xyyyyy_yzz, tg_xyyyyy_zzz, tg_xyyyyz_xxx, tg_xyyyyz_xxy, tg_xyyyyz_xxz, tg_xyyyyz_xyy, tg_xyyyyz_xyz, tg_xyyyyz_xzz, tg_xyyyyz_yyy, tg_xyyyyz_yyz, tg_xyyyyz_yzz, tg_xyyyyz_zzz, tg_xyyyz_yyz, tg_xyyyz_yzz, tg_xyyyz_zzz, tg_xyyyzz_xxx, tg_xyyyzz_xxy, tg_xyyyzz_xxz, tg_xyyyzz_xyy, tg_xyyyzz_xyz, tg_xyyyzz_xzz, tg_xyyyzz_yyy, tg_xyyyzz_yyz, tg_xyyyzz_yzz, tg_xyyyzz_zzz, tg_xyyz_yyz, tg_xyyz_yzz, tg_xyyz_zzz, tg_xyyzz_xyz, tg_xyyzz_yyy, tg_xyyzz_yyz, tg_xyyzz_yz, tg_xyyzz_yzz, tg_xyyzz_zzz, tg_xyyzzz_xxx, tg_xyyzzz_xxy, tg_xyyzzz_xxz, tg_xyyzzz_xyy, tg_xyyzzz_xyz, tg_xyyzzz_xzz, tg_xyyzzz_yyy, tg_xyyzzz_yyz, tg_xyyzzz_yzz, tg_xyyzzz_zzz, tg_xyzz_yyy, tg_xyzz_yyz, tg_xyzz_yzz, tg_xyzzz_yyy, tg_xyzzz_yyz, tg_xyzzz_yzz, tg_xyzzzz_xxx, tg_xyzzzz_xxy, tg_xyzzzz_xxz, tg_xyzzzz_xyy, tg_xyzzzz_xyz, tg_xyzzzz_xzz, tg_xyzzzz_yyy, tg_xyzzzz_yyz, tg_xyzzzz_yzz, tg_xyzzzz_zzz, tg_xzzz_xxz, tg_xzzz_xyz, tg_xzzz_xzz, tg_xzzz_yyy, tg_xzzz_yyz, tg_xzzz_yzz, tg_xzzz_zzz, tg_xzzzz_xxx, tg_xzzzz_xxz, tg_xzzzz_xyz, tg_xzzzz_xz, tg_xzzzz_xzz, tg_xzzzz_yyy, tg_xzzzz_yyz, tg_xzzzz_yz, tg_xzzzz_yzz, tg_xzzzz_zz, tg_xzzzz_zzz, tg_xzzzzz_xxx, tg_xzzzzz_xxy, tg_xzzzzz_xxz, tg_xzzzzz_xyy, tg_xzzzzz_xyz, tg_xzzzzz_xzz, tg_xzzzzz_yyy, tg_xzzzzz_yyz, tg_xzzzzz_yzz, tg_xzzzzz_zzz, tg_yyyy_xxx, tg_yyyy_xxy, tg_yyyy_xxz, tg_yyyy_xyy, tg_yyyy_xyz, tg_yyyy_xzz, tg_yyyy_yyy, tg_yyyy_yyz, tg_yyyy_yzz, tg_yyyy_zzz, tg_yyyyy_xx, tg_yyyyy_xxx, tg_yyyyy_xxy, tg_yyyyy_xxz, tg_yyyyy_xy, tg_yyyyy_xyy, tg_yyyyy_xyz, tg_yyyyy_xz, tg_yyyyy_xzz, tg_yyyyy_yy, tg_yyyyy_yyy, tg_yyyyy_yyz, tg_yyyyy_yz, tg_yyyyy_yzz, tg_yyyyy_zz, tg_yyyyy_zzz, tg_yyyyyy_xxx, tg_yyyyyy_xxy, tg_yyyyyy_xxz, tg_yyyyyy_xyy, tg_yyyyyy_xyz, tg_yyyyyy_xzz, tg_yyyyyy_yyy, tg_yyyyyy_yyz, tg_yyyyyy_yzz, tg_yyyyyy_zzz, tg_yyyyyz_xxx, tg_yyyyyz_xxy, tg_yyyyyz_xxz, tg_yyyyyz_xyy, tg_yyyyyz_xyz, tg_yyyyyz_xzz, tg_yyyyyz_yyy, tg_yyyyyz_yyz, tg_yyyyyz_yzz, tg_yyyyyz_zzz, tg_yyyyz_xxy, tg_yyyyz_xxz, tg_yyyyz_xyy, tg_yyyyz_xyz, tg_yyyyz_xz, tg_yyyyz_xzz, tg_yyyyz_yyy, tg_yyyyz_yyz, tg_yyyyz_yz, tg_yyyyz_yzz, tg_yyyyz_zz, tg_yyyyz_zzz, tg_yyyyzz_xxx, tg_yyyyzz_xxy, tg_yyyyzz_xxz, tg_yyyyzz_xyy, tg_yyyyzz_xyz, tg_yyyyzz_xzz, tg_yyyyzz_yyy, tg_yyyyzz_yyz, tg_yyyyzz_yzz, tg_yyyyzz_zzz, tg_yyyz_xxy, tg_yyyz_xxz, tg_yyyz_xyy, tg_yyyz_xzz, tg_yyyz_yyy, tg_yyyz_yyz, tg_yyyz_yzz, tg_yyyz_zzz, tg_yyyzz_xx, tg_yyyzz_xxx, tg_yyyzz_xxy, tg_yyyzz_xxz, tg_yyyzz_xy, tg_yyyzz_xyy, tg_yyyzz_xyz, tg_yyyzz_xz, tg_yyyzz_xzz, tg_yyyzz_yy, tg_yyyzz_yyy, tg_yyyzz_yyz, tg_yyyzz_yz, tg_yyyzz_yzz, tg_yyyzz_zz, tg_yyyzz_zzz, tg_yyyzzz_xxx, tg_yyyzzz_xxy, tg_yyyzzz_xxz, tg_yyyzzz_xyy, tg_yyyzzz_xyz, tg_yyyzzz_xzz, tg_yyyzzz_yyy, tg_yyyzzz_yyz, tg_yyyzzz_yzz, tg_yyyzzz_zzz, tg_yyzz_xxx, tg_yyzz_xxy, tg_yyzz_xxz, tg_yyzz_xyy, tg_yyzz_xyz, tg_yyzz_xzz, tg_yyzz_yyy, tg_yyzz_yyz, tg_yyzz_yzz, tg_yyzz_zzz, tg_yyzzz_xx, tg_yyzzz_xxx, tg_yyzzz_xxy, tg_yyzzz_xxz, tg_yyzzz_xy, tg_yyzzz_xyy, tg_yyzzz_xyz, tg_yyzzz_xz, tg_yyzzz_xzz, tg_yyzzz_yy, tg_yyzzz_yyy, tg_yyzzz_yyz, tg_yyzzz_yz, tg_yyzzz_yzz, tg_yyzzz_zz, tg_yyzzz_zzz, tg_yyzzzz_xxx, tg_yyzzzz_xxy, tg_yyzzzz_xxz, tg_yyzzzz_xyy, tg_yyzzzz_xyz, tg_yyzzzz_xzz, tg_yyzzzz_yyy, tg_yyzzzz_yyz, tg_yyzzzz_yzz, tg_yyzzzz_zzz, tg_yzzz_xxx, tg_yzzz_xxz, tg_yzzz_xyz, tg_yzzz_xzz, tg_yzzz_yyy, tg_yzzz_yyz, tg_yzzz_yzz, tg_yzzz_zzz, tg_yzzzz_xxx, tg_yzzzz_xxy, tg_yzzzz_xxz, tg_yzzzz_xy, tg_yzzzz_xyy, tg_yzzzz_xyz, tg_yzzzz_xz, tg_yzzzz_xzz, tg_yzzzz_yy, tg_yzzzz_yyy, tg_yzzzz_yyz, tg_yzzzz_yz, tg_yzzzz_yzz, tg_yzzzz_zz, tg_yzzzz_zzz, tg_yzzzzz_xxx, tg_yzzzzz_xxy, tg_yzzzzz_xxz, tg_yzzzzz_xyy, tg_yzzzzz_xyz, tg_yzzzzz_xzz, tg_yzzzzz_yyy, tg_yzzzzz_yyz, tg_yzzzzz_yzz, tg_yzzzzz_zzz, tg_zzzz_xxx, tg_zzzz_xxy, tg_zzzz_xxz, tg_zzzz_xyy, tg_zzzz_xyz, tg_zzzz_xzz, tg_zzzz_yyy, tg_zzzz_yyz, tg_zzzz_yzz, tg_zzzz_zzz, tg_zzzzz_xx, tg_zzzzz_xxx, tg_zzzzz_xxy, tg_zzzzz_xxz, tg_zzzzz_xy, tg_zzzzz_xyy, tg_zzzzz_xyz, tg_zzzzz_xz, tg_zzzzz_xzz, tg_zzzzz_yy, tg_zzzzz_yyy, tg_zzzzz_yyz, tg_zzzzz_yz, tg_zzzzz_yzz, tg_zzzzz_zz, tg_zzzzz_zzz, tg_zzzzzz_xxx, tg_zzzzzz_xxy, tg_zzzzzz_xxz, tg_zzzzzz_xyy, tg_zzzzzz_xyz, tg_zzzzzz_xzz, tg_zzzzzz_yyy, tg_zzzzzz_yyz, tg_zzzzzz_yzz, tg_zzzzzz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxx_xxx[i] = 5.0 * tg_xxxx_xxx[i] * fxi[i] + 3.0 * tg_xxxxx_xx[i] * fxi[i] + tg_xxxxx_xxx[i] * ra_x[i];

        tg_xxxxxx_xxy[i] = 5.0 * tg_xxxx_xxy[i] * fxi[i] + 2.0 * tg_xxxxx_xy[i] * fxi[i] + tg_xxxxx_xxy[i] * ra_x[i];

        tg_xxxxxx_xxz[i] = 5.0 * tg_xxxx_xxz[i] * fxi[i] + 2.0 * tg_xxxxx_xz[i] * fxi[i] + tg_xxxxx_xxz[i] * ra_x[i];

        tg_xxxxxx_xyy[i] = 5.0 * tg_xxxx_xyy[i] * fxi[i] + tg_xxxxx_yy[i] * fxi[i] + tg_xxxxx_xyy[i] * ra_x[i];

        tg_xxxxxx_xyz[i] = 5.0 * tg_xxxx_xyz[i] * fxi[i] + tg_xxxxx_yz[i] * fxi[i] + tg_xxxxx_xyz[i] * ra_x[i];

        tg_xxxxxx_xzz[i] = 5.0 * tg_xxxx_xzz[i] * fxi[i] + tg_xxxxx_zz[i] * fxi[i] + tg_xxxxx_xzz[i] * ra_x[i];

        tg_xxxxxx_yyy[i] = 5.0 * tg_xxxx_yyy[i] * fxi[i] + tg_xxxxx_yyy[i] * ra_x[i];

        tg_xxxxxx_yyz[i] = 5.0 * tg_xxxx_yyz[i] * fxi[i] + tg_xxxxx_yyz[i] * ra_x[i];

        tg_xxxxxx_yzz[i] = 5.0 * tg_xxxx_yzz[i] * fxi[i] + tg_xxxxx_yzz[i] * ra_x[i];

        tg_xxxxxx_zzz[i] = 5.0 * tg_xxxx_zzz[i] * fxi[i] + tg_xxxxx_zzz[i] * ra_x[i];

        tg_xxxxxy_xxx[i] = tg_xxxxx_xxx[i] * ra_y[i];

        tg_xxxxxy_xxy[i] = tg_xxxxx_xx[i] * fxi[i] + tg_xxxxx_xxy[i] * ra_y[i];

        tg_xxxxxy_xxz[i] = tg_xxxxx_xxz[i] * ra_y[i];

        tg_xxxxxy_xyy[i] = 2.0 * tg_xxxxx_xy[i] * fxi[i] + tg_xxxxx_xyy[i] * ra_y[i];

        tg_xxxxxy_xyz[i] = tg_xxxxx_xz[i] * fxi[i] + tg_xxxxx_xyz[i] * ra_y[i];

        tg_xxxxxy_xzz[i] = tg_xxxxx_xzz[i] * ra_y[i];

        tg_xxxxxy_yyy[i] = 4.0 * tg_xxxy_yyy[i] * fxi[i] + tg_xxxxy_yyy[i] * ra_x[i];

        tg_xxxxxy_yyz[i] = 4.0 * tg_xxxy_yyz[i] * fxi[i] + tg_xxxxy_yyz[i] * ra_x[i];

        tg_xxxxxy_yzz[i] = 4.0 * tg_xxxy_yzz[i] * fxi[i] + tg_xxxxy_yzz[i] * ra_x[i];

        tg_xxxxxy_zzz[i] = tg_xxxxx_zzz[i] * ra_y[i];

        tg_xxxxxz_xxx[i] = tg_xxxxx_xxx[i] * ra_z[i];

        tg_xxxxxz_xxy[i] = tg_xxxxx_xxy[i] * ra_z[i];

        tg_xxxxxz_xxz[i] = tg_xxxxx_xx[i] * fxi[i] + tg_xxxxx_xxz[i] * ra_z[i];

        tg_xxxxxz_xyy[i] = tg_xxxxx_xyy[i] * ra_z[i];

        tg_xxxxxz_xyz[i] = tg_xxxxx_xy[i] * fxi[i] + tg_xxxxx_xyz[i] * ra_z[i];

        tg_xxxxxz_xzz[i] = 2.0 * tg_xxxxx_xz[i] * fxi[i] + tg_xxxxx_xzz[i] * ra_z[i];

        tg_xxxxxz_yyy[i] = tg_xxxxx_yyy[i] * ra_z[i];

        tg_xxxxxz_yyz[i] = 4.0 * tg_xxxz_yyz[i] * fxi[i] + tg_xxxxz_yyz[i] * ra_x[i];

        tg_xxxxxz_yzz[i] = 4.0 * tg_xxxz_yzz[i] * fxi[i] + tg_xxxxz_yzz[i] * ra_x[i];

        tg_xxxxxz_zzz[i] = 4.0 * tg_xxxz_zzz[i] * fxi[i] + tg_xxxxz_zzz[i] * ra_x[i];

        tg_xxxxyy_xxx[i] = tg_xxxx_xxx[i] * fxi[i] + tg_xxxxy_xxx[i] * ra_y[i];

        tg_xxxxyy_xxy[i] = 3.0 * tg_xxyy_xxy[i] * fxi[i] + 2.0 * tg_xxxyy_xy[i] * fxi[i] + tg_xxxyy_xxy[i] * ra_x[i];

        tg_xxxxyy_xxz[i] = tg_xxxx_xxz[i] * fxi[i] + tg_xxxxy_xxz[i] * ra_y[i];

        tg_xxxxyy_xyy[i] = 3.0 * tg_xxyy_xyy[i] * fxi[i] + tg_xxxyy_yy[i] * fxi[i] + tg_xxxyy_xyy[i] * ra_x[i];

        tg_xxxxyy_xyz[i] = 3.0 * tg_xxyy_xyz[i] * fxi[i] + tg_xxxyy_yz[i] * fxi[i] + tg_xxxyy_xyz[i] * ra_x[i];

        tg_xxxxyy_xzz[i] = tg_xxxx_xzz[i] * fxi[i] + tg_xxxxy_xzz[i] * ra_y[i];

        tg_xxxxyy_yyy[i] = 3.0 * tg_xxyy_yyy[i] * fxi[i] + tg_xxxyy_yyy[i] * ra_x[i];

        tg_xxxxyy_yyz[i] = 3.0 * tg_xxyy_yyz[i] * fxi[i] + tg_xxxyy_yyz[i] * ra_x[i];

        tg_xxxxyy_yzz[i] = 3.0 * tg_xxyy_yzz[i] * fxi[i] + tg_xxxyy_yzz[i] * ra_x[i];

        tg_xxxxyy_zzz[i] = 3.0 * tg_xxyy_zzz[i] * fxi[i] + tg_xxxyy_zzz[i] * ra_x[i];

        tg_xxxxyz_xxx[i] = tg_xxxxz_xxx[i] * ra_y[i];

        tg_xxxxyz_xxy[i] = tg_xxxxy_xxy[i] * ra_z[i];

        tg_xxxxyz_xxz[i] = tg_xxxxz_xxz[i] * ra_y[i];

        tg_xxxxyz_xyy[i] = tg_xxxxy_xyy[i] * ra_z[i];

        tg_xxxxyz_xyz[i] = tg_xxxxz_xz[i] * fxi[i] + tg_xxxxz_xyz[i] * ra_y[i];

        tg_xxxxyz_xzz[i] = tg_xxxxz_xzz[i] * ra_y[i];

        tg_xxxxyz_yyy[i] = tg_xxxxy_yyy[i] * ra_z[i];

        tg_xxxxyz_yyz[i] = 3.0 * tg_xxyz_yyz[i] * fxi[i] + tg_xxxyz_yyz[i] * ra_x[i];

        tg_xxxxyz_yzz[i] = 3.0 * tg_xxyz_yzz[i] * fxi[i] + tg_xxxyz_yzz[i] * ra_x[i];

        tg_xxxxyz_zzz[i] = tg_xxxxz_zzz[i] * ra_y[i];

        tg_xxxxzz_xxx[i] = tg_xxxx_xxx[i] * fxi[i] + tg_xxxxz_xxx[i] * ra_z[i];

        tg_xxxxzz_xxy[i] = tg_xxxx_xxy[i] * fxi[i] + tg_xxxxz_xxy[i] * ra_z[i];

        tg_xxxxzz_xxz[i] = 3.0 * tg_xxzz_xxz[i] * fxi[i] + 2.0 * tg_xxxzz_xz[i] * fxi[i] + tg_xxxzz_xxz[i] * ra_x[i];

        tg_xxxxzz_xyy[i] = tg_xxxx_xyy[i] * fxi[i] + tg_xxxxz_xyy[i] * ra_z[i];

        tg_xxxxzz_xyz[i] = 3.0 * tg_xxzz_xyz[i] * fxi[i] + tg_xxxzz_yz[i] * fxi[i] + tg_xxxzz_xyz[i] * ra_x[i];

        tg_xxxxzz_xzz[i] = 3.0 * tg_xxzz_xzz[i] * fxi[i] + tg_xxxzz_zz[i] * fxi[i] + tg_xxxzz_xzz[i] * ra_x[i];

        tg_xxxxzz_yyy[i] = 3.0 * tg_xxzz_yyy[i] * fxi[i] + tg_xxxzz_yyy[i] * ra_x[i];

        tg_xxxxzz_yyz[i] = 3.0 * tg_xxzz_yyz[i] * fxi[i] + tg_xxxzz_yyz[i] * ra_x[i];

        tg_xxxxzz_yzz[i] = 3.0 * tg_xxzz_yzz[i] * fxi[i] + tg_xxxzz_yzz[i] * ra_x[i];

        tg_xxxxzz_zzz[i] = 3.0 * tg_xxzz_zzz[i] * fxi[i] + tg_xxxzz_zzz[i] * ra_x[i];

        tg_xxxyyy_xxx[i] = 2.0 * tg_xxxy_xxx[i] * fxi[i] + tg_xxxyy_xxx[i] * ra_y[i];

        tg_xxxyyy_xxy[i] = 2.0 * tg_xyyy_xxy[i] * fxi[i] + 2.0 * tg_xxyyy_xy[i] * fxi[i] + tg_xxyyy_xxy[i] * ra_x[i];

        tg_xxxyyy_xxz[i] = 2.0 * tg_xxxy_xxz[i] * fxi[i] + tg_xxxyy_xxz[i] * ra_y[i];

        tg_xxxyyy_xyy[i] = 2.0 * tg_xyyy_xyy[i] * fxi[i] + tg_xxyyy_yy[i] * fxi[i] + tg_xxyyy_xyy[i] * ra_x[i];

        tg_xxxyyy_xyz[i] = 2.0 * tg_xyyy_xyz[i] * fxi[i] + tg_xxyyy_yz[i] * fxi[i] + tg_xxyyy_xyz[i] * ra_x[i];

        tg_xxxyyy_xzz[i] = 2.0 * tg_xxxy_xzz[i] * fxi[i] + tg_xxxyy_xzz[i] * ra_y[i];

        tg_xxxyyy_yyy[i] = 2.0 * tg_xyyy_yyy[i] * fxi[i] + tg_xxyyy_yyy[i] * ra_x[i];

        tg_xxxyyy_yyz[i] = 2.0 * tg_xyyy_yyz[i] * fxi[i] + tg_xxyyy_yyz[i] * ra_x[i];

        tg_xxxyyy_yzz[i] = 2.0 * tg_xyyy_yzz[i] * fxi[i] + tg_xxyyy_yzz[i] * ra_x[i];

        tg_xxxyyy_zzz[i] = 2.0 * tg_xyyy_zzz[i] * fxi[i] + tg_xxyyy_zzz[i] * ra_x[i];

        tg_xxxyyz_xxx[i] = tg_xxxyy_xxx[i] * ra_z[i];

        tg_xxxyyz_xxy[i] = tg_xxxyy_xxy[i] * ra_z[i];

        tg_xxxyyz_xxz[i] = tg_xxxz_xxz[i] * fxi[i] + tg_xxxyz_xxz[i] * ra_y[i];

        tg_xxxyyz_xyy[i] = tg_xxxyy_xyy[i] * ra_z[i];

        tg_xxxyyz_xyz[i] = tg_xxxyy_xy[i] * fxi[i] + tg_xxxyy_xyz[i] * ra_z[i];

        tg_xxxyyz_xzz[i] = tg_xxxz_xzz[i] * fxi[i] + tg_xxxyz_xzz[i] * ra_y[i];

        tg_xxxyyz_yyy[i] = tg_xxxyy_yyy[i] * ra_z[i];

        tg_xxxyyz_yyz[i] = 2.0 * tg_xyyz_yyz[i] * fxi[i] + tg_xxyyz_yyz[i] * ra_x[i];

        tg_xxxyyz_yzz[i] = 2.0 * tg_xyyz_yzz[i] * fxi[i] + tg_xxyyz_yzz[i] * ra_x[i];

        tg_xxxyyz_zzz[i] = 2.0 * tg_xyyz_zzz[i] * fxi[i] + tg_xxyyz_zzz[i] * ra_x[i];

        tg_xxxyzz_xxx[i] = tg_xxxzz_xxx[i] * ra_y[i];

        tg_xxxyzz_xxy[i] = tg_xxxzz_xx[i] * fxi[i] + tg_xxxzz_xxy[i] * ra_y[i];

        tg_xxxyzz_xxz[i] = tg_xxxzz_xxz[i] * ra_y[i];

        tg_xxxyzz_xyy[i] = 2.0 * tg_xxxzz_xy[i] * fxi[i] + tg_xxxzz_xyy[i] * ra_y[i];

        tg_xxxyzz_xyz[i] = tg_xxxzz_xz[i] * fxi[i] + tg_xxxzz_xyz[i] * ra_y[i];

        tg_xxxyzz_xzz[i] = tg_xxxzz_xzz[i] * ra_y[i];

        tg_xxxyzz_yyy[i] = 2.0 * tg_xyzz_yyy[i] * fxi[i] + tg_xxyzz_yyy[i] * ra_x[i];

        tg_xxxyzz_yyz[i] = 2.0 * tg_xyzz_yyz[i] * fxi[i] + tg_xxyzz_yyz[i] * ra_x[i];

        tg_xxxyzz_yzz[i] = 2.0 * tg_xyzz_yzz[i] * fxi[i] + tg_xxyzz_yzz[i] * ra_x[i];

        tg_xxxyzz_zzz[i] = tg_xxxzz_zzz[i] * ra_y[i];

        tg_xxxzzz_xxx[i] = 2.0 * tg_xxxz_xxx[i] * fxi[i] + tg_xxxzz_xxx[i] * ra_z[i];

        tg_xxxzzz_xxy[i] = 2.0 * tg_xxxz_xxy[i] * fxi[i] + tg_xxxzz_xxy[i] * ra_z[i];

        tg_xxxzzz_xxz[i] = 2.0 * tg_xzzz_xxz[i] * fxi[i] + 2.0 * tg_xxzzz_xz[i] * fxi[i] + tg_xxzzz_xxz[i] * ra_x[i];

        tg_xxxzzz_xyy[i] = 2.0 * tg_xxxz_xyy[i] * fxi[i] + tg_xxxzz_xyy[i] * ra_z[i];

        tg_xxxzzz_xyz[i] = 2.0 * tg_xzzz_xyz[i] * fxi[i] + tg_xxzzz_yz[i] * fxi[i] + tg_xxzzz_xyz[i] * ra_x[i];

        tg_xxxzzz_xzz[i] = 2.0 * tg_xzzz_xzz[i] * fxi[i] + tg_xxzzz_zz[i] * fxi[i] + tg_xxzzz_xzz[i] * ra_x[i];

        tg_xxxzzz_yyy[i] = 2.0 * tg_xzzz_yyy[i] * fxi[i] + tg_xxzzz_yyy[i] * ra_x[i];

        tg_xxxzzz_yyz[i] = 2.0 * tg_xzzz_yyz[i] * fxi[i] + tg_xxzzz_yyz[i] * ra_x[i];

        tg_xxxzzz_yzz[i] = 2.0 * tg_xzzz_yzz[i] * fxi[i] + tg_xxzzz_yzz[i] * ra_x[i];

        tg_xxxzzz_zzz[i] = 2.0 * tg_xzzz_zzz[i] * fxi[i] + tg_xxzzz_zzz[i] * ra_x[i];

        tg_xxyyyy_xxx[i] = 3.0 * tg_xxyy_xxx[i] * fxi[i] + tg_xxyyy_xxx[i] * ra_y[i];

        tg_xxyyyy_xxy[i] = tg_yyyy_xxy[i] * fxi[i] + 2.0 * tg_xyyyy_xy[i] * fxi[i] + tg_xyyyy_xxy[i] * ra_x[i];

        tg_xxyyyy_xxz[i] = 3.0 * tg_xxyy_xxz[i] * fxi[i] + tg_xxyyy_xxz[i] * ra_y[i];

        tg_xxyyyy_xyy[i] = tg_yyyy_xyy[i] * fxi[i] + tg_xyyyy_yy[i] * fxi[i] + tg_xyyyy_xyy[i] * ra_x[i];

        tg_xxyyyy_xyz[i] = tg_yyyy_xyz[i] * fxi[i] + tg_xyyyy_yz[i] * fxi[i] + tg_xyyyy_xyz[i] * ra_x[i];

        tg_xxyyyy_xzz[i] = 3.0 * tg_xxyy_xzz[i] * fxi[i] + tg_xxyyy_xzz[i] * ra_y[i];

        tg_xxyyyy_yyy[i] = tg_yyyy_yyy[i] * fxi[i] + tg_xyyyy_yyy[i] * ra_x[i];

        tg_xxyyyy_yyz[i] = tg_yyyy_yyz[i] * fxi[i] + tg_xyyyy_yyz[i] * ra_x[i];

        tg_xxyyyy_yzz[i] = tg_yyyy_yzz[i] * fxi[i] + tg_xyyyy_yzz[i] * ra_x[i];

        tg_xxyyyy_zzz[i] = tg_yyyy_zzz[i] * fxi[i] + tg_xyyyy_zzz[i] * ra_x[i];

        tg_xxyyyz_xxx[i] = tg_xxyyy_xxx[i] * ra_z[i];

        tg_xxyyyz_xxy[i] = tg_xxyyy_xxy[i] * ra_z[i];

        tg_xxyyyz_xxz[i] = 2.0 * tg_xxyz_xxz[i] * fxi[i] + tg_xxyyz_xxz[i] * ra_y[i];

        tg_xxyyyz_xyy[i] = tg_xxyyy_xyy[i] * ra_z[i];

        tg_xxyyyz_xyz[i] = tg_xxyyy_xy[i] * fxi[i] + tg_xxyyy_xyz[i] * ra_z[i];

        tg_xxyyyz_xzz[i] = 2.0 * tg_xxyz_xzz[i] * fxi[i] + tg_xxyyz_xzz[i] * ra_y[i];

        tg_xxyyyz_yyy[i] = tg_xxyyy_yyy[i] * ra_z[i];

        tg_xxyyyz_yyz[i] = tg_yyyz_yyz[i] * fxi[i] + tg_xyyyz_yyz[i] * ra_x[i];

        tg_xxyyyz_yzz[i] = tg_yyyz_yzz[i] * fxi[i] + tg_xyyyz_yzz[i] * ra_x[i];

        tg_xxyyyz_zzz[i] = tg_yyyz_zzz[i] * fxi[i] + tg_xyyyz_zzz[i] * ra_x[i];

        tg_xxyyzz_xxx[i] = tg_xxzz_xxx[i] * fxi[i] + tg_xxyzz_xxx[i] * ra_y[i];

        tg_xxyyzz_xxy[i] = tg_xxyy_xxy[i] * fxi[i] + tg_xxyyz_xxy[i] * ra_z[i];

        tg_xxyyzz_xxz[i] = tg_xxzz_xxz[i] * fxi[i] + tg_xxyzz_xxz[i] * ra_y[i];

        tg_xxyyzz_xyy[i] = tg_xxyy_xyy[i] * fxi[i] + tg_xxyyz_xyy[i] * ra_z[i];

        tg_xxyyzz_xyz[i] = tg_yyzz_xyz[i] * fxi[i] + tg_xyyzz_yz[i] * fxi[i] + tg_xyyzz_xyz[i] * ra_x[i];

        tg_xxyyzz_xzz[i] = tg_xxzz_xzz[i] * fxi[i] + tg_xxyzz_xzz[i] * ra_y[i];

        tg_xxyyzz_yyy[i] = tg_yyzz_yyy[i] * fxi[i] + tg_xyyzz_yyy[i] * ra_x[i];

        tg_xxyyzz_yyz[i] = tg_yyzz_yyz[i] * fxi[i] + tg_xyyzz_yyz[i] * ra_x[i];

        tg_xxyyzz_yzz[i] = tg_yyzz_yzz[i] * fxi[i] + tg_xyyzz_yzz[i] * ra_x[i];

        tg_xxyyzz_zzz[i] = tg_yyzz_zzz[i] * fxi[i] + tg_xyyzz_zzz[i] * ra_x[i];

        tg_xxyzzz_xxx[i] = tg_xxzzz_xxx[i] * ra_y[i];

        tg_xxyzzz_xxy[i] = tg_xxzzz_xx[i] * fxi[i] + tg_xxzzz_xxy[i] * ra_y[i];

        tg_xxyzzz_xxz[i] = tg_xxzzz_xxz[i] * ra_y[i];

        tg_xxyzzz_xyy[i] = 2.0 * tg_xxzzz_xy[i] * fxi[i] + tg_xxzzz_xyy[i] * ra_y[i];

        tg_xxyzzz_xyz[i] = tg_xxzzz_xz[i] * fxi[i] + tg_xxzzz_xyz[i] * ra_y[i];

        tg_xxyzzz_xzz[i] = tg_xxzzz_xzz[i] * ra_y[i];

        tg_xxyzzz_yyy[i] = tg_yzzz_yyy[i] * fxi[i] + tg_xyzzz_yyy[i] * ra_x[i];

        tg_xxyzzz_yyz[i] = tg_yzzz_yyz[i] * fxi[i] + tg_xyzzz_yyz[i] * ra_x[i];

        tg_xxyzzz_yzz[i] = tg_yzzz_yzz[i] * fxi[i] + tg_xyzzz_yzz[i] * ra_x[i];

        tg_xxyzzz_zzz[i] = tg_xxzzz_zzz[i] * ra_y[i];

        tg_xxzzzz_xxx[i] = 3.0 * tg_xxzz_xxx[i] * fxi[i] + tg_xxzzz_xxx[i] * ra_z[i];

        tg_xxzzzz_xxy[i] = 3.0 * tg_xxzz_xxy[i] * fxi[i] + tg_xxzzz_xxy[i] * ra_z[i];

        tg_xxzzzz_xxz[i] = tg_zzzz_xxz[i] * fxi[i] + 2.0 * tg_xzzzz_xz[i] * fxi[i] + tg_xzzzz_xxz[i] * ra_x[i];

        tg_xxzzzz_xyy[i] = 3.0 * tg_xxzz_xyy[i] * fxi[i] + tg_xxzzz_xyy[i] * ra_z[i];

        tg_xxzzzz_xyz[i] = tg_zzzz_xyz[i] * fxi[i] + tg_xzzzz_yz[i] * fxi[i] + tg_xzzzz_xyz[i] * ra_x[i];

        tg_xxzzzz_xzz[i] = tg_zzzz_xzz[i] * fxi[i] + tg_xzzzz_zz[i] * fxi[i] + tg_xzzzz_xzz[i] * ra_x[i];

        tg_xxzzzz_yyy[i] = tg_zzzz_yyy[i] * fxi[i] + tg_xzzzz_yyy[i] * ra_x[i];

        tg_xxzzzz_yyz[i] = tg_zzzz_yyz[i] * fxi[i] + tg_xzzzz_yyz[i] * ra_x[i];

        tg_xxzzzz_yzz[i] = tg_zzzz_yzz[i] * fxi[i] + tg_xzzzz_yzz[i] * ra_x[i];

        tg_xxzzzz_zzz[i] = tg_zzzz_zzz[i] * fxi[i] + tg_xzzzz_zzz[i] * ra_x[i];

        tg_xyyyyy_xxx[i] = 3.0 * tg_yyyyy_xx[i] * fxi[i] + tg_yyyyy_xxx[i] * ra_x[i];

        tg_xyyyyy_xxy[i] = 2.0 * tg_yyyyy_xy[i] * fxi[i] + tg_yyyyy_xxy[i] * ra_x[i];

        tg_xyyyyy_xxz[i] = 2.0 * tg_yyyyy_xz[i] * fxi[i] + tg_yyyyy_xxz[i] * ra_x[i];

        tg_xyyyyy_xyy[i] = tg_yyyyy_yy[i] * fxi[i] + tg_yyyyy_xyy[i] * ra_x[i];

        tg_xyyyyy_xyz[i] = tg_yyyyy_yz[i] * fxi[i] + tg_yyyyy_xyz[i] * ra_x[i];

        tg_xyyyyy_xzz[i] = tg_yyyyy_zz[i] * fxi[i] + tg_yyyyy_xzz[i] * ra_x[i];

        tg_xyyyyy_yyy[i] = tg_yyyyy_yyy[i] * ra_x[i];

        tg_xyyyyy_yyz[i] = tg_yyyyy_yyz[i] * ra_x[i];

        tg_xyyyyy_yzz[i] = tg_yyyyy_yzz[i] * ra_x[i];

        tg_xyyyyy_zzz[i] = tg_yyyyy_zzz[i] * ra_x[i];

        tg_xyyyyz_xxx[i] = tg_xyyyy_xxx[i] * ra_z[i];

        tg_xyyyyz_xxy[i] = tg_xyyyy_xxy[i] * ra_z[i];

        tg_xyyyyz_xxz[i] = 2.0 * tg_yyyyz_xz[i] * fxi[i] + tg_yyyyz_xxz[i] * ra_x[i];

        tg_xyyyyz_xyy[i] = tg_xyyyy_xyy[i] * ra_z[i];

        tg_xyyyyz_xyz[i] = tg_yyyyz_yz[i] * fxi[i] + tg_yyyyz_xyz[i] * ra_x[i];

        tg_xyyyyz_xzz[i] = tg_yyyyz_zz[i] * fxi[i] + tg_yyyyz_xzz[i] * ra_x[i];

        tg_xyyyyz_yyy[i] = tg_yyyyz_yyy[i] * ra_x[i];

        tg_xyyyyz_yyz[i] = tg_yyyyz_yyz[i] * ra_x[i];

        tg_xyyyyz_yzz[i] = tg_yyyyz_yzz[i] * ra_x[i];

        tg_xyyyyz_zzz[i] = tg_yyyyz_zzz[i] * ra_x[i];

        tg_xyyyzz_xxx[i] = 3.0 * tg_yyyzz_xx[i] * fxi[i] + tg_yyyzz_xxx[i] * ra_x[i];

        tg_xyyyzz_xxy[i] = 2.0 * tg_yyyzz_xy[i] * fxi[i] + tg_yyyzz_xxy[i] * ra_x[i];

        tg_xyyyzz_xxz[i] = 2.0 * tg_yyyzz_xz[i] * fxi[i] + tg_yyyzz_xxz[i] * ra_x[i];

        tg_xyyyzz_xyy[i] = tg_yyyzz_yy[i] * fxi[i] + tg_yyyzz_xyy[i] * ra_x[i];

        tg_xyyyzz_xyz[i] = tg_yyyzz_yz[i] * fxi[i] + tg_yyyzz_xyz[i] * ra_x[i];

        tg_xyyyzz_xzz[i] = tg_yyyzz_zz[i] * fxi[i] + tg_yyyzz_xzz[i] * ra_x[i];

        tg_xyyyzz_yyy[i] = tg_yyyzz_yyy[i] * ra_x[i];

        tg_xyyyzz_yyz[i] = tg_yyyzz_yyz[i] * ra_x[i];

        tg_xyyyzz_yzz[i] = tg_yyyzz_yzz[i] * ra_x[i];

        tg_xyyyzz_zzz[i] = tg_yyyzz_zzz[i] * ra_x[i];

        tg_xyyzzz_xxx[i] = 3.0 * tg_yyzzz_xx[i] * fxi[i] + tg_yyzzz_xxx[i] * ra_x[i];

        tg_xyyzzz_xxy[i] = 2.0 * tg_yyzzz_xy[i] * fxi[i] + tg_yyzzz_xxy[i] * ra_x[i];

        tg_xyyzzz_xxz[i] = 2.0 * tg_yyzzz_xz[i] * fxi[i] + tg_yyzzz_xxz[i] * ra_x[i];

        tg_xyyzzz_xyy[i] = tg_yyzzz_yy[i] * fxi[i] + tg_yyzzz_xyy[i] * ra_x[i];

        tg_xyyzzz_xyz[i] = tg_yyzzz_yz[i] * fxi[i] + tg_yyzzz_xyz[i] * ra_x[i];

        tg_xyyzzz_xzz[i] = tg_yyzzz_zz[i] * fxi[i] + tg_yyzzz_xzz[i] * ra_x[i];

        tg_xyyzzz_yyy[i] = tg_yyzzz_yyy[i] * ra_x[i];

        tg_xyyzzz_yyz[i] = tg_yyzzz_yyz[i] * ra_x[i];

        tg_xyyzzz_yzz[i] = tg_yyzzz_yzz[i] * ra_x[i];

        tg_xyyzzz_zzz[i] = tg_yyzzz_zzz[i] * ra_x[i];

        tg_xyzzzz_xxx[i] = tg_xzzzz_xxx[i] * ra_y[i];

        tg_xyzzzz_xxy[i] = 2.0 * tg_yzzzz_xy[i] * fxi[i] + tg_yzzzz_xxy[i] * ra_x[i];

        tg_xyzzzz_xxz[i] = tg_xzzzz_xxz[i] * ra_y[i];

        tg_xyzzzz_xyy[i] = tg_yzzzz_yy[i] * fxi[i] + tg_yzzzz_xyy[i] * ra_x[i];

        tg_xyzzzz_xyz[i] = tg_yzzzz_yz[i] * fxi[i] + tg_yzzzz_xyz[i] * ra_x[i];

        tg_xyzzzz_xzz[i] = tg_xzzzz_xzz[i] * ra_y[i];

        tg_xyzzzz_yyy[i] = tg_yzzzz_yyy[i] * ra_x[i];

        tg_xyzzzz_yyz[i] = tg_yzzzz_yyz[i] * ra_x[i];

        tg_xyzzzz_yzz[i] = tg_yzzzz_yzz[i] * ra_x[i];

        tg_xyzzzz_zzz[i] = tg_yzzzz_zzz[i] * ra_x[i];

        tg_xzzzzz_xxx[i] = 3.0 * tg_zzzzz_xx[i] * fxi[i] + tg_zzzzz_xxx[i] * ra_x[i];

        tg_xzzzzz_xxy[i] = 2.0 * tg_zzzzz_xy[i] * fxi[i] + tg_zzzzz_xxy[i] * ra_x[i];

        tg_xzzzzz_xxz[i] = 2.0 * tg_zzzzz_xz[i] * fxi[i] + tg_zzzzz_xxz[i] * ra_x[i];

        tg_xzzzzz_xyy[i] = tg_zzzzz_yy[i] * fxi[i] + tg_zzzzz_xyy[i] * ra_x[i];

        tg_xzzzzz_xyz[i] = tg_zzzzz_yz[i] * fxi[i] + tg_zzzzz_xyz[i] * ra_x[i];

        tg_xzzzzz_xzz[i] = tg_zzzzz_zz[i] * fxi[i] + tg_zzzzz_xzz[i] * ra_x[i];

        tg_xzzzzz_yyy[i] = tg_zzzzz_yyy[i] * ra_x[i];

        tg_xzzzzz_yyz[i] = tg_zzzzz_yyz[i] * ra_x[i];

        tg_xzzzzz_yzz[i] = tg_zzzzz_yzz[i] * ra_x[i];

        tg_xzzzzz_zzz[i] = tg_zzzzz_zzz[i] * ra_x[i];

        tg_yyyyyy_xxx[i] = 5.0 * tg_yyyy_xxx[i] * fxi[i] + tg_yyyyy_xxx[i] * ra_y[i];

        tg_yyyyyy_xxy[i] = 5.0 * tg_yyyy_xxy[i] * fxi[i] + tg_yyyyy_xx[i] * fxi[i] + tg_yyyyy_xxy[i] * ra_y[i];

        tg_yyyyyy_xxz[i] = 5.0 * tg_yyyy_xxz[i] * fxi[i] + tg_yyyyy_xxz[i] * ra_y[i];

        tg_yyyyyy_xyy[i] = 5.0 * tg_yyyy_xyy[i] * fxi[i] + 2.0 * tg_yyyyy_xy[i] * fxi[i] + tg_yyyyy_xyy[i] * ra_y[i];

        tg_yyyyyy_xyz[i] = 5.0 * tg_yyyy_xyz[i] * fxi[i] + tg_yyyyy_xz[i] * fxi[i] + tg_yyyyy_xyz[i] * ra_y[i];

        tg_yyyyyy_xzz[i] = 5.0 * tg_yyyy_xzz[i] * fxi[i] + tg_yyyyy_xzz[i] * ra_y[i];

        tg_yyyyyy_yyy[i] = 5.0 * tg_yyyy_yyy[i] * fxi[i] + 3.0 * tg_yyyyy_yy[i] * fxi[i] + tg_yyyyy_yyy[i] * ra_y[i];

        tg_yyyyyy_yyz[i] = 5.0 * tg_yyyy_yyz[i] * fxi[i] + 2.0 * tg_yyyyy_yz[i] * fxi[i] + tg_yyyyy_yyz[i] * ra_y[i];

        tg_yyyyyy_yzz[i] = 5.0 * tg_yyyy_yzz[i] * fxi[i] + tg_yyyyy_zz[i] * fxi[i] + tg_yyyyy_yzz[i] * ra_y[i];

        tg_yyyyyy_zzz[i] = 5.0 * tg_yyyy_zzz[i] * fxi[i] + tg_yyyyy_zzz[i] * ra_y[i];

        tg_yyyyyz_xxx[i] = tg_yyyyy_xxx[i] * ra_z[i];

        tg_yyyyyz_xxy[i] = tg_yyyyy_xxy[i] * ra_z[i];

        tg_yyyyyz_xxz[i] = 4.0 * tg_yyyz_xxz[i] * fxi[i] + tg_yyyyz_xxz[i] * ra_y[i];

        tg_yyyyyz_xyy[i] = tg_yyyyy_xyy[i] * ra_z[i];

        tg_yyyyyz_xyz[i] = tg_yyyyy_xy[i] * fxi[i] + tg_yyyyy_xyz[i] * ra_z[i];

        tg_yyyyyz_xzz[i] = 4.0 * tg_yyyz_xzz[i] * fxi[i] + tg_yyyyz_xzz[i] * ra_y[i];

        tg_yyyyyz_yyy[i] = tg_yyyyy_yyy[i] * ra_z[i];

        tg_yyyyyz_yyz[i] = tg_yyyyy_yy[i] * fxi[i] + tg_yyyyy_yyz[i] * ra_z[i];

        tg_yyyyyz_yzz[i] = 2.0 * tg_yyyyy_yz[i] * fxi[i] + tg_yyyyy_yzz[i] * ra_z[i];

        tg_yyyyyz_zzz[i] = 4.0 * tg_yyyz_zzz[i] * fxi[i] + tg_yyyyz_zzz[i] * ra_y[i];

        tg_yyyyzz_xxx[i] = 3.0 * tg_yyzz_xxx[i] * fxi[i] + tg_yyyzz_xxx[i] * ra_y[i];

        tg_yyyyzz_xxy[i] = tg_yyyy_xxy[i] * fxi[i] + tg_yyyyz_xxy[i] * ra_z[i];

        tg_yyyyzz_xxz[i] = 3.0 * tg_yyzz_xxz[i] * fxi[i] + tg_yyyzz_xxz[i] * ra_y[i];

        tg_yyyyzz_xyy[i] = tg_yyyy_xyy[i] * fxi[i] + tg_yyyyz_xyy[i] * ra_z[i];

        tg_yyyyzz_xyz[i] = 3.0 * tg_yyzz_xyz[i] * fxi[i] + tg_yyyzz_xz[i] * fxi[i] + tg_yyyzz_xyz[i] * ra_y[i];

        tg_yyyyzz_xzz[i] = 3.0 * tg_yyzz_xzz[i] * fxi[i] + tg_yyyzz_xzz[i] * ra_y[i];

        tg_yyyyzz_yyy[i] = tg_yyyy_yyy[i] * fxi[i] + tg_yyyyz_yyy[i] * ra_z[i];

        tg_yyyyzz_yyz[i] = 3.0 * tg_yyzz_yyz[i] * fxi[i] + 2.0 * tg_yyyzz_yz[i] * fxi[i] + tg_yyyzz_yyz[i] * ra_y[i];

        tg_yyyyzz_yzz[i] = 3.0 * tg_yyzz_yzz[i] * fxi[i] + tg_yyyzz_zz[i] * fxi[i] + tg_yyyzz_yzz[i] * ra_y[i];

        tg_yyyyzz_zzz[i] = 3.0 * tg_yyzz_zzz[i] * fxi[i] + tg_yyyzz_zzz[i] * ra_y[i];

        tg_yyyzzz_xxx[i] = 2.0 * tg_yzzz_xxx[i] * fxi[i] + tg_yyzzz_xxx[i] * ra_y[i];

        tg_yyyzzz_xxy[i] = 2.0 * tg_yyyz_xxy[i] * fxi[i] + tg_yyyzz_xxy[i] * ra_z[i];

        tg_yyyzzz_xxz[i] = 2.0 * tg_yzzz_xxz[i] * fxi[i] + tg_yyzzz_xxz[i] * ra_y[i];

        tg_yyyzzz_xyy[i] = 2.0 * tg_yyyz_xyy[i] * fxi[i] + tg_yyyzz_xyy[i] * ra_z[i];

        tg_yyyzzz_xyz[i] = 2.0 * tg_yzzz_xyz[i] * fxi[i] + tg_yyzzz_xz[i] * fxi[i] + tg_yyzzz_xyz[i] * ra_y[i];

        tg_yyyzzz_xzz[i] = 2.0 * tg_yzzz_xzz[i] * fxi[i] + tg_yyzzz_xzz[i] * ra_y[i];

        tg_yyyzzz_yyy[i] = 2.0 * tg_yyyz_yyy[i] * fxi[i] + tg_yyyzz_yyy[i] * ra_z[i];

        tg_yyyzzz_yyz[i] = 2.0 * tg_yzzz_yyz[i] * fxi[i] + 2.0 * tg_yyzzz_yz[i] * fxi[i] + tg_yyzzz_yyz[i] * ra_y[i];

        tg_yyyzzz_yzz[i] = 2.0 * tg_yzzz_yzz[i] * fxi[i] + tg_yyzzz_zz[i] * fxi[i] + tg_yyzzz_yzz[i] * ra_y[i];

        tg_yyyzzz_zzz[i] = 2.0 * tg_yzzz_zzz[i] * fxi[i] + tg_yyzzz_zzz[i] * ra_y[i];

        tg_yyzzzz_xxx[i] = tg_zzzz_xxx[i] * fxi[i] + tg_yzzzz_xxx[i] * ra_y[i];

        tg_yyzzzz_xxy[i] = 3.0 * tg_yyzz_xxy[i] * fxi[i] + tg_yyzzz_xxy[i] * ra_z[i];

        tg_yyzzzz_xxz[i] = tg_zzzz_xxz[i] * fxi[i] + tg_yzzzz_xxz[i] * ra_y[i];

        tg_yyzzzz_xyy[i] = 3.0 * tg_yyzz_xyy[i] * fxi[i] + tg_yyzzz_xyy[i] * ra_z[i];

        tg_yyzzzz_xyz[i] = tg_zzzz_xyz[i] * fxi[i] + tg_yzzzz_xz[i] * fxi[i] + tg_yzzzz_xyz[i] * ra_y[i];

        tg_yyzzzz_xzz[i] = tg_zzzz_xzz[i] * fxi[i] + tg_yzzzz_xzz[i] * ra_y[i];

        tg_yyzzzz_yyy[i] = 3.0 * tg_yyzz_yyy[i] * fxi[i] + tg_yyzzz_yyy[i] * ra_z[i];

        tg_yyzzzz_yyz[i] = tg_zzzz_yyz[i] * fxi[i] + 2.0 * tg_yzzzz_yz[i] * fxi[i] + tg_yzzzz_yyz[i] * ra_y[i];

        tg_yyzzzz_yzz[i] = tg_zzzz_yzz[i] * fxi[i] + tg_yzzzz_zz[i] * fxi[i] + tg_yzzzz_yzz[i] * ra_y[i];

        tg_yyzzzz_zzz[i] = tg_zzzz_zzz[i] * fxi[i] + tg_yzzzz_zzz[i] * ra_y[i];

        tg_yzzzzz_xxx[i] = tg_zzzzz_xxx[i] * ra_y[i];

        tg_yzzzzz_xxy[i] = tg_zzzzz_xx[i] * fxi[i] + tg_zzzzz_xxy[i] * ra_y[i];

        tg_yzzzzz_xxz[i] = tg_zzzzz_xxz[i] * ra_y[i];

        tg_yzzzzz_xyy[i] = 2.0 * tg_zzzzz_xy[i] * fxi[i] + tg_zzzzz_xyy[i] * ra_y[i];

        tg_yzzzzz_xyz[i] = tg_zzzzz_xz[i] * fxi[i] + tg_zzzzz_xyz[i] * ra_y[i];

        tg_yzzzzz_xzz[i] = tg_zzzzz_xzz[i] * ra_y[i];

        tg_yzzzzz_yyy[i] = 3.0 * tg_zzzzz_yy[i] * fxi[i] + tg_zzzzz_yyy[i] * ra_y[i];

        tg_yzzzzz_yyz[i] = 2.0 * tg_zzzzz_yz[i] * fxi[i] + tg_zzzzz_yyz[i] * ra_y[i];

        tg_yzzzzz_yzz[i] = tg_zzzzz_zz[i] * fxi[i] + tg_zzzzz_yzz[i] * ra_y[i];

        tg_yzzzzz_zzz[i] = tg_zzzzz_zzz[i] * ra_y[i];

        tg_zzzzzz_xxx[i] = 5.0 * tg_zzzz_xxx[i] * fxi[i] + tg_zzzzz_xxx[i] * ra_z[i];

        tg_zzzzzz_xxy[i] = 5.0 * tg_zzzz_xxy[i] * fxi[i] + tg_zzzzz_xxy[i] * ra_z[i];

        tg_zzzzzz_xxz[i] = 5.0 * tg_zzzz_xxz[i] * fxi[i] + tg_zzzzz_xx[i] * fxi[i] + tg_zzzzz_xxz[i] * ra_z[i];

        tg_zzzzzz_xyy[i] = 5.0 * tg_zzzz_xyy[i] * fxi[i] + tg_zzzzz_xyy[i] * ra_z[i];

        tg_zzzzzz_xyz[i] = 5.0 * tg_zzzz_xyz[i] * fxi[i] + tg_zzzzz_xy[i] * fxi[i] + tg_zzzzz_xyz[i] * ra_z[i];

        tg_zzzzzz_xzz[i] = 5.0 * tg_zzzz_xzz[i] * fxi[i] + 2.0 * tg_zzzzz_xz[i] * fxi[i] + tg_zzzzz_xzz[i] * ra_z[i];

        tg_zzzzzz_yyy[i] = 5.0 * tg_zzzz_yyy[i] * fxi[i] + tg_zzzzz_yyy[i] * ra_z[i];

        tg_zzzzzz_yyz[i] = 5.0 * tg_zzzz_yyz[i] * fxi[i] + tg_zzzzz_yy[i] * fxi[i] + tg_zzzzz_yyz[i] * ra_z[i];

        tg_zzzzzz_yzz[i] = 5.0 * tg_zzzz_yzz[i] * fxi[i] + 2.0 * tg_zzzzz_yz[i] * fxi[i] + tg_zzzzz_yzz[i] * ra_z[i];

        tg_zzzzzz_zzz[i] = 5.0 * tg_zzzz_zzz[i] * fxi[i] + 3.0 * tg_zzzzz_zz[i] * fxi[i] + tg_zzzzz_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

