#include "ThreeCenterOverlapGradientPrimRecHF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hf(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hf,
                              const size_t idx_gf,
                              const size_t idx_hd,
                              const size_t idx_hf,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_gf + 13);

    auto ts_xxxy_xyz = pbuffer.data(idx_gf + 14);

    auto ts_xxxy_xzz = pbuffer.data(idx_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_gf + 18);

    auto ts_xxxy_zzz = pbuffer.data(idx_gf + 19);

    auto ts_xxxz_xxx = pbuffer.data(idx_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_gf + 25);

    auto ts_xxxz_yyy = pbuffer.data(idx_gf + 26);

    auto ts_xxxz_yyz = pbuffer.data(idx_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_gf + 39);

    auto ts_xxyz_xxx = pbuffer.data(idx_gf + 40);

    auto ts_xxyz_xxy = pbuffer.data(idx_gf + 41);

    auto ts_xxyz_xxz = pbuffer.data(idx_gf + 42);

    auto ts_xxyz_xyy = pbuffer.data(idx_gf + 43);

    auto ts_xxyz_xyz = pbuffer.data(idx_gf + 44);

    auto ts_xxyz_xzz = pbuffer.data(idx_gf + 45);

    auto ts_xxyz_yyy = pbuffer.data(idx_gf + 46);

    auto ts_xxyz_yyz = pbuffer.data(idx_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_gf + 48);

    auto ts_xxyz_zzz = pbuffer.data(idx_gf + 49);

    auto ts_xxzz_xxx = pbuffer.data(idx_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_gf + 59);

    auto ts_xyyy_xxx = pbuffer.data(idx_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_gf + 61);

    auto ts_xyyy_xxz = pbuffer.data(idx_gf + 62);

    auto ts_xyyy_xyy = pbuffer.data(idx_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_gf + 64);

    auto ts_xyyy_xzz = pbuffer.data(idx_gf + 65);

    auto ts_xyyy_yyy = pbuffer.data(idx_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_gf + 69);

    auto ts_xyyz_xxx = pbuffer.data(idx_gf + 70);

    auto ts_xyyz_xxy = pbuffer.data(idx_gf + 71);

    auto ts_xyyz_xxz = pbuffer.data(idx_gf + 72);

    auto ts_xyyz_xyy = pbuffer.data(idx_gf + 73);

    auto ts_xyyz_xyz = pbuffer.data(idx_gf + 74);

    auto ts_xyyz_xzz = pbuffer.data(idx_gf + 75);

    auto ts_xyyz_yyy = pbuffer.data(idx_gf + 76);

    auto ts_xyyz_yyz = pbuffer.data(idx_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_gf + 79);

    auto ts_xyzz_xxx = pbuffer.data(idx_gf + 80);

    auto ts_xyzz_xxy = pbuffer.data(idx_gf + 81);

    auto ts_xyzz_xxz = pbuffer.data(idx_gf + 82);

    auto ts_xyzz_xyy = pbuffer.data(idx_gf + 83);

    auto ts_xyzz_xyz = pbuffer.data(idx_gf + 84);

    auto ts_xyzz_xzz = pbuffer.data(idx_gf + 85);

    auto ts_xyzz_yyy = pbuffer.data(idx_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_gf + 88);

    auto ts_xyzz_zzz = pbuffer.data(idx_gf + 89);

    auto ts_xzzz_xxx = pbuffer.data(idx_gf + 90);

    auto ts_xzzz_xxy = pbuffer.data(idx_gf + 91);

    auto ts_xzzz_xxz = pbuffer.data(idx_gf + 92);

    auto ts_xzzz_xyy = pbuffer.data(idx_gf + 93);

    auto ts_xzzz_xyz = pbuffer.data(idx_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_gf + 109);

    auto ts_yyyz_xxx = pbuffer.data(idx_gf + 110);

    auto ts_yyyz_xxy = pbuffer.data(idx_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_gf + 149);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_hd + 5);

    auto ts_xxxxy_xx = pbuffer.data(idx_hd + 6);

    auto ts_xxxxy_xy = pbuffer.data(idx_hd + 7);

    auto ts_xxxxy_xz = pbuffer.data(idx_hd + 8);

    auto ts_xxxxy_yy = pbuffer.data(idx_hd + 9);

    auto ts_xxxxy_yz = pbuffer.data(idx_hd + 10);

    auto ts_xxxxy_zz = pbuffer.data(idx_hd + 11);

    auto ts_xxxxz_xx = pbuffer.data(idx_hd + 12);

    auto ts_xxxxz_xy = pbuffer.data(idx_hd + 13);

    auto ts_xxxxz_xz = pbuffer.data(idx_hd + 14);

    auto ts_xxxxz_yy = pbuffer.data(idx_hd + 15);

    auto ts_xxxxz_yz = pbuffer.data(idx_hd + 16);

    auto ts_xxxxz_zz = pbuffer.data(idx_hd + 17);

    auto ts_xxxyy_xx = pbuffer.data(idx_hd + 18);

    auto ts_xxxyy_xy = pbuffer.data(idx_hd + 19);

    auto ts_xxxyy_xz = pbuffer.data(idx_hd + 20);

    auto ts_xxxyy_yy = pbuffer.data(idx_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_hd + 22);

    auto ts_xxxyy_zz = pbuffer.data(idx_hd + 23);

    auto ts_xxxyz_xx = pbuffer.data(idx_hd + 24);

    auto ts_xxxyz_xy = pbuffer.data(idx_hd + 25);

    auto ts_xxxyz_xz = pbuffer.data(idx_hd + 26);

    auto ts_xxxyz_yy = pbuffer.data(idx_hd + 27);

    auto ts_xxxyz_yz = pbuffer.data(idx_hd + 28);

    auto ts_xxxyz_zz = pbuffer.data(idx_hd + 29);

    auto ts_xxxzz_xx = pbuffer.data(idx_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_hd + 32);

    auto ts_xxxzz_yy = pbuffer.data(idx_hd + 33);

    auto ts_xxxzz_yz = pbuffer.data(idx_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_hd + 35);

    auto ts_xxyyy_xx = pbuffer.data(idx_hd + 36);

    auto ts_xxyyy_xy = pbuffer.data(idx_hd + 37);

    auto ts_xxyyy_xz = pbuffer.data(idx_hd + 38);

    auto ts_xxyyy_yy = pbuffer.data(idx_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_hd + 40);

    auto ts_xxyyy_zz = pbuffer.data(idx_hd + 41);

    auto ts_xxyyz_xx = pbuffer.data(idx_hd + 42);

    auto ts_xxyyz_xy = pbuffer.data(idx_hd + 43);

    auto ts_xxyyz_xz = pbuffer.data(idx_hd + 44);

    auto ts_xxyyz_yy = pbuffer.data(idx_hd + 45);

    auto ts_xxyyz_yz = pbuffer.data(idx_hd + 46);

    auto ts_xxyyz_zz = pbuffer.data(idx_hd + 47);

    auto ts_xxyzz_xx = pbuffer.data(idx_hd + 48);

    auto ts_xxyzz_xy = pbuffer.data(idx_hd + 49);

    auto ts_xxyzz_xz = pbuffer.data(idx_hd + 50);

    auto ts_xxyzz_yy = pbuffer.data(idx_hd + 51);

    auto ts_xxyzz_yz = pbuffer.data(idx_hd + 52);

    auto ts_xxyzz_zz = pbuffer.data(idx_hd + 53);

    auto ts_xxzzz_xx = pbuffer.data(idx_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_hd + 56);

    auto ts_xxzzz_yy = pbuffer.data(idx_hd + 57);

    auto ts_xxzzz_yz = pbuffer.data(idx_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_hd + 59);

    auto ts_xyyyy_xx = pbuffer.data(idx_hd + 60);

    auto ts_xyyyy_xy = pbuffer.data(idx_hd + 61);

    auto ts_xyyyy_xz = pbuffer.data(idx_hd + 62);

    auto ts_xyyyy_yy = pbuffer.data(idx_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_hd + 64);

    auto ts_xyyyy_zz = pbuffer.data(idx_hd + 65);

    auto ts_xyyyz_xx = pbuffer.data(idx_hd + 66);

    auto ts_xyyyz_xy = pbuffer.data(idx_hd + 67);

    auto ts_xyyyz_xz = pbuffer.data(idx_hd + 68);

    auto ts_xyyyz_yy = pbuffer.data(idx_hd + 69);

    auto ts_xyyyz_yz = pbuffer.data(idx_hd + 70);

    auto ts_xyyyz_zz = pbuffer.data(idx_hd + 71);

    auto ts_xyyzz_xx = pbuffer.data(idx_hd + 72);

    auto ts_xyyzz_xy = pbuffer.data(idx_hd + 73);

    auto ts_xyyzz_xz = pbuffer.data(idx_hd + 74);

    auto ts_xyyzz_yy = pbuffer.data(idx_hd + 75);

    auto ts_xyyzz_yz = pbuffer.data(idx_hd + 76);

    auto ts_xyyzz_zz = pbuffer.data(idx_hd + 77);

    auto ts_xyzzz_xx = pbuffer.data(idx_hd + 78);

    auto ts_xyzzz_xy = pbuffer.data(idx_hd + 79);

    auto ts_xyzzz_xz = pbuffer.data(idx_hd + 80);

    auto ts_xyzzz_yy = pbuffer.data(idx_hd + 81);

    auto ts_xyzzz_yz = pbuffer.data(idx_hd + 82);

    auto ts_xyzzz_zz = pbuffer.data(idx_hd + 83);

    auto ts_xzzzz_xx = pbuffer.data(idx_hd + 84);

    auto ts_xzzzz_xy = pbuffer.data(idx_hd + 85);

    auto ts_xzzzz_xz = pbuffer.data(idx_hd + 86);

    auto ts_xzzzz_yy = pbuffer.data(idx_hd + 87);

    auto ts_xzzzz_yz = pbuffer.data(idx_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_hd + 95);

    auto ts_yyyyz_xx = pbuffer.data(idx_hd + 96);

    auto ts_yyyyz_xy = pbuffer.data(idx_hd + 97);

    auto ts_yyyyz_xz = pbuffer.data(idx_hd + 98);

    auto ts_yyyyz_yy = pbuffer.data(idx_hd + 99);

    auto ts_yyyyz_yz = pbuffer.data(idx_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_hd + 101);

    auto ts_yyyzz_xx = pbuffer.data(idx_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_hd + 107);

    auto ts_yyzzz_xx = pbuffer.data(idx_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_hd + 113);

    auto ts_yzzzz_xx = pbuffer.data(idx_hd + 114);

    auto ts_yzzzz_xy = pbuffer.data(idx_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_hd + 125);

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_hf + 9);

    auto ts_xxxxy_xxx = pbuffer.data(idx_hf + 10);

    auto ts_xxxxy_xxy = pbuffer.data(idx_hf + 11);

    auto ts_xxxxy_xxz = pbuffer.data(idx_hf + 12);

    auto ts_xxxxy_xyy = pbuffer.data(idx_hf + 13);

    auto ts_xxxxy_xyz = pbuffer.data(idx_hf + 14);

    auto ts_xxxxy_xzz = pbuffer.data(idx_hf + 15);

    auto ts_xxxxy_yyy = pbuffer.data(idx_hf + 16);

    auto ts_xxxxy_yyz = pbuffer.data(idx_hf + 17);

    auto ts_xxxxy_yzz = pbuffer.data(idx_hf + 18);

    auto ts_xxxxy_zzz = pbuffer.data(idx_hf + 19);

    auto ts_xxxxz_xxx = pbuffer.data(idx_hf + 20);

    auto ts_xxxxz_xxy = pbuffer.data(idx_hf + 21);

    auto ts_xxxxz_xxz = pbuffer.data(idx_hf + 22);

    auto ts_xxxxz_xyy = pbuffer.data(idx_hf + 23);

    auto ts_xxxxz_xyz = pbuffer.data(idx_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_hf + 25);

    auto ts_xxxxz_yyy = pbuffer.data(idx_hf + 26);

    auto ts_xxxxz_yyz = pbuffer.data(idx_hf + 27);

    auto ts_xxxxz_yzz = pbuffer.data(idx_hf + 28);

    auto ts_xxxxz_zzz = pbuffer.data(idx_hf + 29);

    auto ts_xxxyy_xxx = pbuffer.data(idx_hf + 30);

    auto ts_xxxyy_xxy = pbuffer.data(idx_hf + 31);

    auto ts_xxxyy_xxz = pbuffer.data(idx_hf + 32);

    auto ts_xxxyy_xyy = pbuffer.data(idx_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_hf + 34);

    auto ts_xxxyy_xzz = pbuffer.data(idx_hf + 35);

    auto ts_xxxyy_yyy = pbuffer.data(idx_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_hf + 38);

    auto ts_xxxyy_zzz = pbuffer.data(idx_hf + 39);

    auto ts_xxxyz_xxx = pbuffer.data(idx_hf + 40);

    auto ts_xxxyz_xxy = pbuffer.data(idx_hf + 41);

    auto ts_xxxyz_xxz = pbuffer.data(idx_hf + 42);

    auto ts_xxxyz_xyy = pbuffer.data(idx_hf + 43);

    auto ts_xxxyz_xyz = pbuffer.data(idx_hf + 44);

    auto ts_xxxyz_xzz = pbuffer.data(idx_hf + 45);

    auto ts_xxxyz_yyy = pbuffer.data(idx_hf + 46);

    auto ts_xxxyz_yyz = pbuffer.data(idx_hf + 47);

    auto ts_xxxyz_yzz = pbuffer.data(idx_hf + 48);

    auto ts_xxxyz_zzz = pbuffer.data(idx_hf + 49);

    auto ts_xxxzz_xxx = pbuffer.data(idx_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_hf + 55);

    auto ts_xxxzz_yyy = pbuffer.data(idx_hf + 56);

    auto ts_xxxzz_yyz = pbuffer.data(idx_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_hf + 59);

    auto ts_xxyyy_xxx = pbuffer.data(idx_hf + 60);

    auto ts_xxyyy_xxy = pbuffer.data(idx_hf + 61);

    auto ts_xxyyy_xxz = pbuffer.data(idx_hf + 62);

    auto ts_xxyyy_xyy = pbuffer.data(idx_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_hf + 64);

    auto ts_xxyyy_xzz = pbuffer.data(idx_hf + 65);

    auto ts_xxyyy_yyy = pbuffer.data(idx_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_hf + 68);

    auto ts_xxyyy_zzz = pbuffer.data(idx_hf + 69);

    auto ts_xxyyz_xxx = pbuffer.data(idx_hf + 70);

    auto ts_xxyyz_xxy = pbuffer.data(idx_hf + 71);

    auto ts_xxyyz_xxz = pbuffer.data(idx_hf + 72);

    auto ts_xxyyz_xyy = pbuffer.data(idx_hf + 73);

    auto ts_xxyyz_xyz = pbuffer.data(idx_hf + 74);

    auto ts_xxyyz_xzz = pbuffer.data(idx_hf + 75);

    auto ts_xxyyz_yyy = pbuffer.data(idx_hf + 76);

    auto ts_xxyyz_yyz = pbuffer.data(idx_hf + 77);

    auto ts_xxyyz_yzz = pbuffer.data(idx_hf + 78);

    auto ts_xxyyz_zzz = pbuffer.data(idx_hf + 79);

    auto ts_xxyzz_xxx = pbuffer.data(idx_hf + 80);

    auto ts_xxyzz_xxy = pbuffer.data(idx_hf + 81);

    auto ts_xxyzz_xxz = pbuffer.data(idx_hf + 82);

    auto ts_xxyzz_xyy = pbuffer.data(idx_hf + 83);

    auto ts_xxyzz_xyz = pbuffer.data(idx_hf + 84);

    auto ts_xxyzz_xzz = pbuffer.data(idx_hf + 85);

    auto ts_xxyzz_yyy = pbuffer.data(idx_hf + 86);

    auto ts_xxyzz_yyz = pbuffer.data(idx_hf + 87);

    auto ts_xxyzz_yzz = pbuffer.data(idx_hf + 88);

    auto ts_xxyzz_zzz = pbuffer.data(idx_hf + 89);

    auto ts_xxzzz_xxx = pbuffer.data(idx_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_hf + 95);

    auto ts_xxzzz_yyy = pbuffer.data(idx_hf + 96);

    auto ts_xxzzz_yyz = pbuffer.data(idx_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_hf + 99);

    auto ts_xyyyy_xxx = pbuffer.data(idx_hf + 100);

    auto ts_xyyyy_xxy = pbuffer.data(idx_hf + 101);

    auto ts_xyyyy_xxz = pbuffer.data(idx_hf + 102);

    auto ts_xyyyy_xyy = pbuffer.data(idx_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_hf + 104);

    auto ts_xyyyy_xzz = pbuffer.data(idx_hf + 105);

    auto ts_xyyyy_yyy = pbuffer.data(idx_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_hf + 108);

    auto ts_xyyyy_zzz = pbuffer.data(idx_hf + 109);

    auto ts_xyyyz_xxx = pbuffer.data(idx_hf + 110);

    auto ts_xyyyz_xxy = pbuffer.data(idx_hf + 111);

    auto ts_xyyyz_xxz = pbuffer.data(idx_hf + 112);

    auto ts_xyyyz_xyy = pbuffer.data(idx_hf + 113);

    auto ts_xyyyz_xyz = pbuffer.data(idx_hf + 114);

    auto ts_xyyyz_xzz = pbuffer.data(idx_hf + 115);

    auto ts_xyyyz_yyy = pbuffer.data(idx_hf + 116);

    auto ts_xyyyz_yyz = pbuffer.data(idx_hf + 117);

    auto ts_xyyyz_yzz = pbuffer.data(idx_hf + 118);

    auto ts_xyyyz_zzz = pbuffer.data(idx_hf + 119);

    auto ts_xyyzz_xxx = pbuffer.data(idx_hf + 120);

    auto ts_xyyzz_xxy = pbuffer.data(idx_hf + 121);

    auto ts_xyyzz_xxz = pbuffer.data(idx_hf + 122);

    auto ts_xyyzz_xyy = pbuffer.data(idx_hf + 123);

    auto ts_xyyzz_xyz = pbuffer.data(idx_hf + 124);

    auto ts_xyyzz_xzz = pbuffer.data(idx_hf + 125);

    auto ts_xyyzz_yyy = pbuffer.data(idx_hf + 126);

    auto ts_xyyzz_yyz = pbuffer.data(idx_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_hf + 128);

    auto ts_xyyzz_zzz = pbuffer.data(idx_hf + 129);

    auto ts_xyzzz_xxx = pbuffer.data(idx_hf + 130);

    auto ts_xyzzz_xxy = pbuffer.data(idx_hf + 131);

    auto ts_xyzzz_xxz = pbuffer.data(idx_hf + 132);

    auto ts_xyzzz_xyy = pbuffer.data(idx_hf + 133);

    auto ts_xyzzz_xyz = pbuffer.data(idx_hf + 134);

    auto ts_xyzzz_xzz = pbuffer.data(idx_hf + 135);

    auto ts_xyzzz_yyy = pbuffer.data(idx_hf + 136);

    auto ts_xyzzz_yyz = pbuffer.data(idx_hf + 137);

    auto ts_xyzzz_yzz = pbuffer.data(idx_hf + 138);

    auto ts_xyzzz_zzz = pbuffer.data(idx_hf + 139);

    auto ts_xzzzz_xxx = pbuffer.data(idx_hf + 140);

    auto ts_xzzzz_xxy = pbuffer.data(idx_hf + 141);

    auto ts_xzzzz_xxz = pbuffer.data(idx_hf + 142);

    auto ts_xzzzz_xyy = pbuffer.data(idx_hf + 143);

    auto ts_xzzzz_xyz = pbuffer.data(idx_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_hf + 145);

    auto ts_xzzzz_yyy = pbuffer.data(idx_hf + 146);

    auto ts_xzzzz_yyz = pbuffer.data(idx_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_hf + 159);

    auto ts_yyyyz_xxx = pbuffer.data(idx_hf + 160);

    auto ts_yyyyz_xxy = pbuffer.data(idx_hf + 161);

    auto ts_yyyyz_xxz = pbuffer.data(idx_hf + 162);

    auto ts_yyyyz_xyy = pbuffer.data(idx_hf + 163);

    auto ts_yyyyz_xyz = pbuffer.data(idx_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_hf + 165);

    auto ts_yyyyz_yyy = pbuffer.data(idx_hf + 166);

    auto ts_yyyyz_yyz = pbuffer.data(idx_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_hf + 169);

    auto ts_yyyzz_xxx = pbuffer.data(idx_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_hf + 179);

    auto ts_yyzzz_xxx = pbuffer.data(idx_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_hf + 189);

    auto ts_yzzzz_xxx = pbuffer.data(idx_hf + 190);

    auto ts_yzzzz_xxy = pbuffer.data(idx_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_hf + 209);

    // Set up 0-10 components of targeted buffer : HF

    auto gs_x_xxxxx_xxx = pbuffer.data(idx_g_hf);

    auto gs_x_xxxxx_xxy = pbuffer.data(idx_g_hf + 1);

    auto gs_x_xxxxx_xxz = pbuffer.data(idx_g_hf + 2);

    auto gs_x_xxxxx_xyy = pbuffer.data(idx_g_hf + 3);

    auto gs_x_xxxxx_xyz = pbuffer.data(idx_g_hf + 4);

    auto gs_x_xxxxx_xzz = pbuffer.data(idx_g_hf + 5);

    auto gs_x_xxxxx_yyy = pbuffer.data(idx_g_hf + 6);

    auto gs_x_xxxxx_yyz = pbuffer.data(idx_g_hf + 7);

    auto gs_x_xxxxx_yzz = pbuffer.data(idx_g_hf + 8);

    auto gs_x_xxxxx_zzz = pbuffer.data(idx_g_hf + 9);

    #pragma omp simd aligned(gc_x, gs_x_xxxxx_xxx, gs_x_xxxxx_xxy, gs_x_xxxxx_xxz, gs_x_xxxxx_xyy, gs_x_xxxxx_xyz, gs_x_xxxxx_xzz, gs_x_xxxxx_yyy, gs_x_xxxxx_yyz, gs_x_xxxxx_yzz, gs_x_xxxxx_zzz, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xzz, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yzz, ts_xxxx_zzz, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_yy, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yz, ts_xxxxx_yzz, ts_xxxxx_zz, ts_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_xxx[i] = 10.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxy[i] = 10.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxz[i] = 10.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyy[i] = 10.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyz[i] = 10.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xzz[i] = 10.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyy[i] = 10.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyz[i] = 10.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yzz[i] = 10.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_zzz[i] = 10.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 10-20 components of targeted buffer : HF

    auto gs_x_xxxxy_xxx = pbuffer.data(idx_g_hf + 10);

    auto gs_x_xxxxy_xxy = pbuffer.data(idx_g_hf + 11);

    auto gs_x_xxxxy_xxz = pbuffer.data(idx_g_hf + 12);

    auto gs_x_xxxxy_xyy = pbuffer.data(idx_g_hf + 13);

    auto gs_x_xxxxy_xyz = pbuffer.data(idx_g_hf + 14);

    auto gs_x_xxxxy_xzz = pbuffer.data(idx_g_hf + 15);

    auto gs_x_xxxxy_yyy = pbuffer.data(idx_g_hf + 16);

    auto gs_x_xxxxy_yyz = pbuffer.data(idx_g_hf + 17);

    auto gs_x_xxxxy_yzz = pbuffer.data(idx_g_hf + 18);

    auto gs_x_xxxxy_zzz = pbuffer.data(idx_g_hf + 19);

    #pragma omp simd aligned(gc_x, gs_x_xxxxy_xxx, gs_x_xxxxy_xxy, gs_x_xxxxy_xxz, gs_x_xxxxy_xyy, gs_x_xxxxy_xyz, gs_x_xxxxy_xzz, gs_x_xxxxy_yyy, gs_x_xxxxy_yyz, gs_x_xxxxy_yzz, gs_x_xxxxy_zzz, ts_xxxxy_xx, ts_xxxxy_xxx, ts_xxxxy_xxy, ts_xxxxy_xxz, ts_xxxxy_xy, ts_xxxxy_xyy, ts_xxxxy_xyz, ts_xxxxy_xz, ts_xxxxy_xzz, ts_xxxxy_yy, ts_xxxxy_yyy, ts_xxxxy_yyz, ts_xxxxy_yz, ts_xxxxy_yzz, ts_xxxxy_zz, ts_xxxxy_zzz, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xzz, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yzz, ts_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxy_xxx[i] = 8.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxy[i] = 8.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxz[i] = 8.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyy[i] = 8.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyz[i] = 8.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xzz[i] = 8.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyy[i] = 8.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyz[i] = 8.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yzz[i] = 8.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_zzz[i] = 8.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 20-30 components of targeted buffer : HF

    auto gs_x_xxxxz_xxx = pbuffer.data(idx_g_hf + 20);

    auto gs_x_xxxxz_xxy = pbuffer.data(idx_g_hf + 21);

    auto gs_x_xxxxz_xxz = pbuffer.data(idx_g_hf + 22);

    auto gs_x_xxxxz_xyy = pbuffer.data(idx_g_hf + 23);

    auto gs_x_xxxxz_xyz = pbuffer.data(idx_g_hf + 24);

    auto gs_x_xxxxz_xzz = pbuffer.data(idx_g_hf + 25);

    auto gs_x_xxxxz_yyy = pbuffer.data(idx_g_hf + 26);

    auto gs_x_xxxxz_yyz = pbuffer.data(idx_g_hf + 27);

    auto gs_x_xxxxz_yzz = pbuffer.data(idx_g_hf + 28);

    auto gs_x_xxxxz_zzz = pbuffer.data(idx_g_hf + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxxxz_xxx, gs_x_xxxxz_xxy, gs_x_xxxxz_xxz, gs_x_xxxxz_xyy, gs_x_xxxxz_xyz, gs_x_xxxxz_xzz, gs_x_xxxxz_yyy, gs_x_xxxxz_yyz, gs_x_xxxxz_yzz, gs_x_xxxxz_zzz, ts_xxxxz_xx, ts_xxxxz_xxx, ts_xxxxz_xxy, ts_xxxxz_xxz, ts_xxxxz_xy, ts_xxxxz_xyy, ts_xxxxz_xyz, ts_xxxxz_xz, ts_xxxxz_xzz, ts_xxxxz_yy, ts_xxxxz_yyy, ts_xxxxz_yyz, ts_xxxxz_yz, ts_xxxxz_yzz, ts_xxxxz_zz, ts_xxxxz_zzz, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xzz, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yzz, ts_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxz_xxx[i] = 8.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxy[i] = 8.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxz[i] = 8.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyy[i] = 8.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyz[i] = 8.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xzz[i] = 8.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyy[i] = 8.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyz[i] = 8.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yzz[i] = 8.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_zzz[i] = 8.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-40 components of targeted buffer : HF

    auto gs_x_xxxyy_xxx = pbuffer.data(idx_g_hf + 30);

    auto gs_x_xxxyy_xxy = pbuffer.data(idx_g_hf + 31);

    auto gs_x_xxxyy_xxz = pbuffer.data(idx_g_hf + 32);

    auto gs_x_xxxyy_xyy = pbuffer.data(idx_g_hf + 33);

    auto gs_x_xxxyy_xyz = pbuffer.data(idx_g_hf + 34);

    auto gs_x_xxxyy_xzz = pbuffer.data(idx_g_hf + 35);

    auto gs_x_xxxyy_yyy = pbuffer.data(idx_g_hf + 36);

    auto gs_x_xxxyy_yyz = pbuffer.data(idx_g_hf + 37);

    auto gs_x_xxxyy_yzz = pbuffer.data(idx_g_hf + 38);

    auto gs_x_xxxyy_zzz = pbuffer.data(idx_g_hf + 39);

    #pragma omp simd aligned(gc_x, gs_x_xxxyy_xxx, gs_x_xxxyy_xxy, gs_x_xxxyy_xxz, gs_x_xxxyy_xyy, gs_x_xxxyy_xyz, gs_x_xxxyy_xzz, gs_x_xxxyy_yyy, gs_x_xxxyy_yyz, gs_x_xxxyy_yzz, gs_x_xxxyy_zzz, ts_xxxyy_xx, ts_xxxyy_xxx, ts_xxxyy_xxy, ts_xxxyy_xxz, ts_xxxyy_xy, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_xz, ts_xxxyy_xzz, ts_xxxyy_yy, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yz, ts_xxxyy_yzz, ts_xxxyy_zz, ts_xxxyy_zzz, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xzz, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yzz, ts_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyy_xxx[i] = 6.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxy[i] = 6.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxz[i] = 6.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyy[i] = 6.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyz[i] = 6.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xzz[i] = 6.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyy[i] = 6.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyz[i] = 6.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yzz[i] = 6.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_zzz[i] = 6.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 40-50 components of targeted buffer : HF

    auto gs_x_xxxyz_xxx = pbuffer.data(idx_g_hf + 40);

    auto gs_x_xxxyz_xxy = pbuffer.data(idx_g_hf + 41);

    auto gs_x_xxxyz_xxz = pbuffer.data(idx_g_hf + 42);

    auto gs_x_xxxyz_xyy = pbuffer.data(idx_g_hf + 43);

    auto gs_x_xxxyz_xyz = pbuffer.data(idx_g_hf + 44);

    auto gs_x_xxxyz_xzz = pbuffer.data(idx_g_hf + 45);

    auto gs_x_xxxyz_yyy = pbuffer.data(idx_g_hf + 46);

    auto gs_x_xxxyz_yyz = pbuffer.data(idx_g_hf + 47);

    auto gs_x_xxxyz_yzz = pbuffer.data(idx_g_hf + 48);

    auto gs_x_xxxyz_zzz = pbuffer.data(idx_g_hf + 49);

    #pragma omp simd aligned(gc_x, gs_x_xxxyz_xxx, gs_x_xxxyz_xxy, gs_x_xxxyz_xxz, gs_x_xxxyz_xyy, gs_x_xxxyz_xyz, gs_x_xxxyz_xzz, gs_x_xxxyz_yyy, gs_x_xxxyz_yyz, gs_x_xxxyz_yzz, gs_x_xxxyz_zzz, ts_xxxyz_xx, ts_xxxyz_xxx, ts_xxxyz_xxy, ts_xxxyz_xxz, ts_xxxyz_xy, ts_xxxyz_xyy, ts_xxxyz_xyz, ts_xxxyz_xz, ts_xxxyz_xzz, ts_xxxyz_yy, ts_xxxyz_yyy, ts_xxxyz_yyz, ts_xxxyz_yz, ts_xxxyz_yzz, ts_xxxyz_zz, ts_xxxyz_zzz, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xzz, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yzz, ts_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyz_xxx[i] = 6.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxy[i] = 6.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxz[i] = 6.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyy[i] = 6.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyz[i] = 6.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xzz[i] = 6.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyy[i] = 6.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyz[i] = 6.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yzz[i] = 6.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_zzz[i] = 6.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 50-60 components of targeted buffer : HF

    auto gs_x_xxxzz_xxx = pbuffer.data(idx_g_hf + 50);

    auto gs_x_xxxzz_xxy = pbuffer.data(idx_g_hf + 51);

    auto gs_x_xxxzz_xxz = pbuffer.data(idx_g_hf + 52);

    auto gs_x_xxxzz_xyy = pbuffer.data(idx_g_hf + 53);

    auto gs_x_xxxzz_xyz = pbuffer.data(idx_g_hf + 54);

    auto gs_x_xxxzz_xzz = pbuffer.data(idx_g_hf + 55);

    auto gs_x_xxxzz_yyy = pbuffer.data(idx_g_hf + 56);

    auto gs_x_xxxzz_yyz = pbuffer.data(idx_g_hf + 57);

    auto gs_x_xxxzz_yzz = pbuffer.data(idx_g_hf + 58);

    auto gs_x_xxxzz_zzz = pbuffer.data(idx_g_hf + 59);

    #pragma omp simd aligned(gc_x, gs_x_xxxzz_xxx, gs_x_xxxzz_xxy, gs_x_xxxzz_xxz, gs_x_xxxzz_xyy, gs_x_xxxzz_xyz, gs_x_xxxzz_xzz, gs_x_xxxzz_yyy, gs_x_xxxzz_yyz, gs_x_xxxzz_yzz, gs_x_xxxzz_zzz, ts_xxxzz_xx, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xxz, ts_xxxzz_xy, ts_xxxzz_xyy, ts_xxxzz_xyz, ts_xxxzz_xz, ts_xxxzz_xzz, ts_xxxzz_yy, ts_xxxzz_yyy, ts_xxxzz_yyz, ts_xxxzz_yz, ts_xxxzz_yzz, ts_xxxzz_zz, ts_xxxzz_zzz, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xzz, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yzz, ts_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxzz_xxx[i] = 6.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxy[i] = 6.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxz[i] = 6.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyy[i] = 6.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyz[i] = 6.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xzz[i] = 6.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyy[i] = 6.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyz[i] = 6.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yzz[i] = 6.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_zzz[i] = 6.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-70 components of targeted buffer : HF

    auto gs_x_xxyyy_xxx = pbuffer.data(idx_g_hf + 60);

    auto gs_x_xxyyy_xxy = pbuffer.data(idx_g_hf + 61);

    auto gs_x_xxyyy_xxz = pbuffer.data(idx_g_hf + 62);

    auto gs_x_xxyyy_xyy = pbuffer.data(idx_g_hf + 63);

    auto gs_x_xxyyy_xyz = pbuffer.data(idx_g_hf + 64);

    auto gs_x_xxyyy_xzz = pbuffer.data(idx_g_hf + 65);

    auto gs_x_xxyyy_yyy = pbuffer.data(idx_g_hf + 66);

    auto gs_x_xxyyy_yyz = pbuffer.data(idx_g_hf + 67);

    auto gs_x_xxyyy_yzz = pbuffer.data(idx_g_hf + 68);

    auto gs_x_xxyyy_zzz = pbuffer.data(idx_g_hf + 69);

    #pragma omp simd aligned(gc_x, gs_x_xxyyy_xxx, gs_x_xxyyy_xxy, gs_x_xxyyy_xxz, gs_x_xxyyy_xyy, gs_x_xxyyy_xyz, gs_x_xxyyy_xzz, gs_x_xxyyy_yyy, gs_x_xxyyy_yyz, gs_x_xxyyy_yzz, gs_x_xxyyy_zzz, ts_xxyyy_xx, ts_xxyyy_xxx, ts_xxyyy_xxy, ts_xxyyy_xxz, ts_xxyyy_xy, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_xz, ts_xxyyy_xzz, ts_xxyyy_yy, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yz, ts_xxyyy_yzz, ts_xxyyy_zz, ts_xxyyy_zzz, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xzz, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yzz, ts_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyy_xxx[i] = 4.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxy[i] = 4.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxz[i] = 4.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyy[i] = 4.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyz[i] = 4.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xzz[i] = 4.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyy[i] = 4.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyz[i] = 4.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yzz[i] = 4.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_zzz[i] = 4.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 70-80 components of targeted buffer : HF

    auto gs_x_xxyyz_xxx = pbuffer.data(idx_g_hf + 70);

    auto gs_x_xxyyz_xxy = pbuffer.data(idx_g_hf + 71);

    auto gs_x_xxyyz_xxz = pbuffer.data(idx_g_hf + 72);

    auto gs_x_xxyyz_xyy = pbuffer.data(idx_g_hf + 73);

    auto gs_x_xxyyz_xyz = pbuffer.data(idx_g_hf + 74);

    auto gs_x_xxyyz_xzz = pbuffer.data(idx_g_hf + 75);

    auto gs_x_xxyyz_yyy = pbuffer.data(idx_g_hf + 76);

    auto gs_x_xxyyz_yyz = pbuffer.data(idx_g_hf + 77);

    auto gs_x_xxyyz_yzz = pbuffer.data(idx_g_hf + 78);

    auto gs_x_xxyyz_zzz = pbuffer.data(idx_g_hf + 79);

    #pragma omp simd aligned(gc_x, gs_x_xxyyz_xxx, gs_x_xxyyz_xxy, gs_x_xxyyz_xxz, gs_x_xxyyz_xyy, gs_x_xxyyz_xyz, gs_x_xxyyz_xzz, gs_x_xxyyz_yyy, gs_x_xxyyz_yyz, gs_x_xxyyz_yzz, gs_x_xxyyz_zzz, ts_xxyyz_xx, ts_xxyyz_xxx, ts_xxyyz_xxy, ts_xxyyz_xxz, ts_xxyyz_xy, ts_xxyyz_xyy, ts_xxyyz_xyz, ts_xxyyz_xz, ts_xxyyz_xzz, ts_xxyyz_yy, ts_xxyyz_yyy, ts_xxyyz_yyz, ts_xxyyz_yz, ts_xxyyz_yzz, ts_xxyyz_zz, ts_xxyyz_zzz, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xzz, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yzz, ts_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyz_xxx[i] = 4.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxy[i] = 4.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxz[i] = 4.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyy[i] = 4.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyz[i] = 4.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xzz[i] = 4.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyy[i] = 4.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyz[i] = 4.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yzz[i] = 4.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_zzz[i] = 4.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 80-90 components of targeted buffer : HF

    auto gs_x_xxyzz_xxx = pbuffer.data(idx_g_hf + 80);

    auto gs_x_xxyzz_xxy = pbuffer.data(idx_g_hf + 81);

    auto gs_x_xxyzz_xxz = pbuffer.data(idx_g_hf + 82);

    auto gs_x_xxyzz_xyy = pbuffer.data(idx_g_hf + 83);

    auto gs_x_xxyzz_xyz = pbuffer.data(idx_g_hf + 84);

    auto gs_x_xxyzz_xzz = pbuffer.data(idx_g_hf + 85);

    auto gs_x_xxyzz_yyy = pbuffer.data(idx_g_hf + 86);

    auto gs_x_xxyzz_yyz = pbuffer.data(idx_g_hf + 87);

    auto gs_x_xxyzz_yzz = pbuffer.data(idx_g_hf + 88);

    auto gs_x_xxyzz_zzz = pbuffer.data(idx_g_hf + 89);

    #pragma omp simd aligned(gc_x, gs_x_xxyzz_xxx, gs_x_xxyzz_xxy, gs_x_xxyzz_xxz, gs_x_xxyzz_xyy, gs_x_xxyzz_xyz, gs_x_xxyzz_xzz, gs_x_xxyzz_yyy, gs_x_xxyzz_yyz, gs_x_xxyzz_yzz, gs_x_xxyzz_zzz, ts_xxyzz_xx, ts_xxyzz_xxx, ts_xxyzz_xxy, ts_xxyzz_xxz, ts_xxyzz_xy, ts_xxyzz_xyy, ts_xxyzz_xyz, ts_xxyzz_xz, ts_xxyzz_xzz, ts_xxyzz_yy, ts_xxyzz_yyy, ts_xxyzz_yyz, ts_xxyzz_yz, ts_xxyzz_yzz, ts_xxyzz_zz, ts_xxyzz_zzz, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xzz, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yzz, ts_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyzz_xxx[i] = 4.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxy[i] = 4.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxz[i] = 4.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyy[i] = 4.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyz[i] = 4.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xzz[i] = 4.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyy[i] = 4.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyz[i] = 4.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yzz[i] = 4.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_zzz[i] = 4.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-100 components of targeted buffer : HF

    auto gs_x_xxzzz_xxx = pbuffer.data(idx_g_hf + 90);

    auto gs_x_xxzzz_xxy = pbuffer.data(idx_g_hf + 91);

    auto gs_x_xxzzz_xxz = pbuffer.data(idx_g_hf + 92);

    auto gs_x_xxzzz_xyy = pbuffer.data(idx_g_hf + 93);

    auto gs_x_xxzzz_xyz = pbuffer.data(idx_g_hf + 94);

    auto gs_x_xxzzz_xzz = pbuffer.data(idx_g_hf + 95);

    auto gs_x_xxzzz_yyy = pbuffer.data(idx_g_hf + 96);

    auto gs_x_xxzzz_yyz = pbuffer.data(idx_g_hf + 97);

    auto gs_x_xxzzz_yzz = pbuffer.data(idx_g_hf + 98);

    auto gs_x_xxzzz_zzz = pbuffer.data(idx_g_hf + 99);

    #pragma omp simd aligned(gc_x, gs_x_xxzzz_xxx, gs_x_xxzzz_xxy, gs_x_xxzzz_xxz, gs_x_xxzzz_xyy, gs_x_xxzzz_xyz, gs_x_xxzzz_xzz, gs_x_xxzzz_yyy, gs_x_xxzzz_yyz, gs_x_xxzzz_yzz, gs_x_xxzzz_zzz, ts_xxzzz_xx, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xxz, ts_xxzzz_xy, ts_xxzzz_xyy, ts_xxzzz_xyz, ts_xxzzz_xz, ts_xxzzz_xzz, ts_xxzzz_yy, ts_xxzzz_yyy, ts_xxzzz_yyz, ts_xxzzz_yz, ts_xxzzz_yzz, ts_xxzzz_zz, ts_xxzzz_zzz, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xzz, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yzz, ts_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzzz_xxx[i] = 4.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxy[i] = 4.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxz[i] = 4.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyy[i] = 4.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyz[i] = 4.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xzz[i] = 4.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyy[i] = 4.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyz[i] = 4.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yzz[i] = 4.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_zzz[i] = 4.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 100-110 components of targeted buffer : HF

    auto gs_x_xyyyy_xxx = pbuffer.data(idx_g_hf + 100);

    auto gs_x_xyyyy_xxy = pbuffer.data(idx_g_hf + 101);

    auto gs_x_xyyyy_xxz = pbuffer.data(idx_g_hf + 102);

    auto gs_x_xyyyy_xyy = pbuffer.data(idx_g_hf + 103);

    auto gs_x_xyyyy_xyz = pbuffer.data(idx_g_hf + 104);

    auto gs_x_xyyyy_xzz = pbuffer.data(idx_g_hf + 105);

    auto gs_x_xyyyy_yyy = pbuffer.data(idx_g_hf + 106);

    auto gs_x_xyyyy_yyz = pbuffer.data(idx_g_hf + 107);

    auto gs_x_xyyyy_yzz = pbuffer.data(idx_g_hf + 108);

    auto gs_x_xyyyy_zzz = pbuffer.data(idx_g_hf + 109);

    #pragma omp simd aligned(gc_x, gs_x_xyyyy_xxx, gs_x_xyyyy_xxy, gs_x_xyyyy_xxz, gs_x_xyyyy_xyy, gs_x_xyyyy_xyz, gs_x_xyyyy_xzz, gs_x_xyyyy_yyy, gs_x_xyyyy_yyz, gs_x_xyyyy_yzz, gs_x_xyyyy_zzz, ts_xyyyy_xx, ts_xyyyy_xxx, ts_xyyyy_xxy, ts_xyyyy_xxz, ts_xyyyy_xy, ts_xyyyy_xyy, ts_xyyyy_xyz, ts_xyyyy_xz, ts_xyyyy_xzz, ts_xyyyy_yy, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yz, ts_xyyyy_yzz, ts_xyyyy_zz, ts_xyyyy_zzz, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xzz, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yzz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyy_xxx[i] = 2.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxy[i] = 2.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxz[i] = 2.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyy[i] = 2.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyz[i] = 2.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xzz[i] = 2.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyy[i] = 2.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyz[i] = 2.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yzz[i] = 2.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_zzz[i] = 2.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 110-120 components of targeted buffer : HF

    auto gs_x_xyyyz_xxx = pbuffer.data(idx_g_hf + 110);

    auto gs_x_xyyyz_xxy = pbuffer.data(idx_g_hf + 111);

    auto gs_x_xyyyz_xxz = pbuffer.data(idx_g_hf + 112);

    auto gs_x_xyyyz_xyy = pbuffer.data(idx_g_hf + 113);

    auto gs_x_xyyyz_xyz = pbuffer.data(idx_g_hf + 114);

    auto gs_x_xyyyz_xzz = pbuffer.data(idx_g_hf + 115);

    auto gs_x_xyyyz_yyy = pbuffer.data(idx_g_hf + 116);

    auto gs_x_xyyyz_yyz = pbuffer.data(idx_g_hf + 117);

    auto gs_x_xyyyz_yzz = pbuffer.data(idx_g_hf + 118);

    auto gs_x_xyyyz_zzz = pbuffer.data(idx_g_hf + 119);

    #pragma omp simd aligned(gc_x, gs_x_xyyyz_xxx, gs_x_xyyyz_xxy, gs_x_xyyyz_xxz, gs_x_xyyyz_xyy, gs_x_xyyyz_xyz, gs_x_xyyyz_xzz, gs_x_xyyyz_yyy, gs_x_xyyyz_yyz, gs_x_xyyyz_yzz, gs_x_xyyyz_zzz, ts_xyyyz_xx, ts_xyyyz_xxx, ts_xyyyz_xxy, ts_xyyyz_xxz, ts_xyyyz_xy, ts_xyyyz_xyy, ts_xyyyz_xyz, ts_xyyyz_xz, ts_xyyyz_xzz, ts_xyyyz_yy, ts_xyyyz_yyy, ts_xyyyz_yyz, ts_xyyyz_yz, ts_xyyyz_yzz, ts_xyyyz_zz, ts_xyyyz_zzz, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xzz, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yzz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyz_xxx[i] = 2.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxy[i] = 2.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxz[i] = 2.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyy[i] = 2.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyz[i] = 2.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xzz[i] = 2.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyy[i] = 2.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyz[i] = 2.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yzz[i] = 2.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_zzz[i] = 2.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-130 components of targeted buffer : HF

    auto gs_x_xyyzz_xxx = pbuffer.data(idx_g_hf + 120);

    auto gs_x_xyyzz_xxy = pbuffer.data(idx_g_hf + 121);

    auto gs_x_xyyzz_xxz = pbuffer.data(idx_g_hf + 122);

    auto gs_x_xyyzz_xyy = pbuffer.data(idx_g_hf + 123);

    auto gs_x_xyyzz_xyz = pbuffer.data(idx_g_hf + 124);

    auto gs_x_xyyzz_xzz = pbuffer.data(idx_g_hf + 125);

    auto gs_x_xyyzz_yyy = pbuffer.data(idx_g_hf + 126);

    auto gs_x_xyyzz_yyz = pbuffer.data(idx_g_hf + 127);

    auto gs_x_xyyzz_yzz = pbuffer.data(idx_g_hf + 128);

    auto gs_x_xyyzz_zzz = pbuffer.data(idx_g_hf + 129);

    #pragma omp simd aligned(gc_x, gs_x_xyyzz_xxx, gs_x_xyyzz_xxy, gs_x_xyyzz_xxz, gs_x_xyyzz_xyy, gs_x_xyyzz_xyz, gs_x_xyyzz_xzz, gs_x_xyyzz_yyy, gs_x_xyyzz_yyz, gs_x_xyyzz_yzz, gs_x_xyyzz_zzz, ts_xyyzz_xx, ts_xyyzz_xxx, ts_xyyzz_xxy, ts_xyyzz_xxz, ts_xyyzz_xy, ts_xyyzz_xyy, ts_xyyzz_xyz, ts_xyyzz_xz, ts_xyyzz_xzz, ts_xyyzz_yy, ts_xyyzz_yyy, ts_xyyzz_yyz, ts_xyyzz_yz, ts_xyyzz_yzz, ts_xyyzz_zz, ts_xyyzz_zzz, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xzz, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yzz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyzz_xxx[i] = 2.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxy[i] = 2.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxz[i] = 2.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyy[i] = 2.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyz[i] = 2.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xzz[i] = 2.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyy[i] = 2.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyz[i] = 2.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yzz[i] = 2.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_zzz[i] = 2.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 130-140 components of targeted buffer : HF

    auto gs_x_xyzzz_xxx = pbuffer.data(idx_g_hf + 130);

    auto gs_x_xyzzz_xxy = pbuffer.data(idx_g_hf + 131);

    auto gs_x_xyzzz_xxz = pbuffer.data(idx_g_hf + 132);

    auto gs_x_xyzzz_xyy = pbuffer.data(idx_g_hf + 133);

    auto gs_x_xyzzz_xyz = pbuffer.data(idx_g_hf + 134);

    auto gs_x_xyzzz_xzz = pbuffer.data(idx_g_hf + 135);

    auto gs_x_xyzzz_yyy = pbuffer.data(idx_g_hf + 136);

    auto gs_x_xyzzz_yyz = pbuffer.data(idx_g_hf + 137);

    auto gs_x_xyzzz_yzz = pbuffer.data(idx_g_hf + 138);

    auto gs_x_xyzzz_zzz = pbuffer.data(idx_g_hf + 139);

    #pragma omp simd aligned(gc_x, gs_x_xyzzz_xxx, gs_x_xyzzz_xxy, gs_x_xyzzz_xxz, gs_x_xyzzz_xyy, gs_x_xyzzz_xyz, gs_x_xyzzz_xzz, gs_x_xyzzz_yyy, gs_x_xyzzz_yyz, gs_x_xyzzz_yzz, gs_x_xyzzz_zzz, ts_xyzzz_xx, ts_xyzzz_xxx, ts_xyzzz_xxy, ts_xyzzz_xxz, ts_xyzzz_xy, ts_xyzzz_xyy, ts_xyzzz_xyz, ts_xyzzz_xz, ts_xyzzz_xzz, ts_xyzzz_yy, ts_xyzzz_yyy, ts_xyzzz_yyz, ts_xyzzz_yz, ts_xyzzz_yzz, ts_xyzzz_zz, ts_xyzzz_zzz, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xzz, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yzz, ts_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzzz_xxx[i] = 2.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxy[i] = 2.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxz[i] = 2.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyy[i] = 2.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyz[i] = 2.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xzz[i] = 2.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyy[i] = 2.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyz[i] = 2.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yzz[i] = 2.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_zzz[i] = 2.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 140-150 components of targeted buffer : HF

    auto gs_x_xzzzz_xxx = pbuffer.data(idx_g_hf + 140);

    auto gs_x_xzzzz_xxy = pbuffer.data(idx_g_hf + 141);

    auto gs_x_xzzzz_xxz = pbuffer.data(idx_g_hf + 142);

    auto gs_x_xzzzz_xyy = pbuffer.data(idx_g_hf + 143);

    auto gs_x_xzzzz_xyz = pbuffer.data(idx_g_hf + 144);

    auto gs_x_xzzzz_xzz = pbuffer.data(idx_g_hf + 145);

    auto gs_x_xzzzz_yyy = pbuffer.data(idx_g_hf + 146);

    auto gs_x_xzzzz_yyz = pbuffer.data(idx_g_hf + 147);

    auto gs_x_xzzzz_yzz = pbuffer.data(idx_g_hf + 148);

    auto gs_x_xzzzz_zzz = pbuffer.data(idx_g_hf + 149);

    #pragma omp simd aligned(gc_x, gs_x_xzzzz_xxx, gs_x_xzzzz_xxy, gs_x_xzzzz_xxz, gs_x_xzzzz_xyy, gs_x_xzzzz_xyz, gs_x_xzzzz_xzz, gs_x_xzzzz_yyy, gs_x_xzzzz_yyz, gs_x_xzzzz_yzz, gs_x_xzzzz_zzz, ts_xzzzz_xx, ts_xzzzz_xxx, ts_xzzzz_xxy, ts_xzzzz_xxz, ts_xzzzz_xy, ts_xzzzz_xyy, ts_xzzzz_xyz, ts_xzzzz_xz, ts_xzzzz_xzz, ts_xzzzz_yy, ts_xzzzz_yyy, ts_xzzzz_yyz, ts_xzzzz_yz, ts_xzzzz_yzz, ts_xzzzz_zz, ts_xzzzz_zzz, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzzz_xxx[i] = 2.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxy[i] = 2.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxz[i] = 2.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyy[i] = 2.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyz[i] = 2.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xzz[i] = 2.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyy[i] = 2.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyz[i] = 2.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yzz[i] = 2.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_zzz[i] = 2.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 150-160 components of targeted buffer : HF

    auto gs_x_yyyyy_xxx = pbuffer.data(idx_g_hf + 150);

    auto gs_x_yyyyy_xxy = pbuffer.data(idx_g_hf + 151);

    auto gs_x_yyyyy_xxz = pbuffer.data(idx_g_hf + 152);

    auto gs_x_yyyyy_xyy = pbuffer.data(idx_g_hf + 153);

    auto gs_x_yyyyy_xyz = pbuffer.data(idx_g_hf + 154);

    auto gs_x_yyyyy_xzz = pbuffer.data(idx_g_hf + 155);

    auto gs_x_yyyyy_yyy = pbuffer.data(idx_g_hf + 156);

    auto gs_x_yyyyy_yyz = pbuffer.data(idx_g_hf + 157);

    auto gs_x_yyyyy_yzz = pbuffer.data(idx_g_hf + 158);

    auto gs_x_yyyyy_zzz = pbuffer.data(idx_g_hf + 159);

    #pragma omp simd aligned(gc_x, gs_x_yyyyy_xxx, gs_x_yyyyy_xxy, gs_x_yyyyy_xxz, gs_x_yyyyy_xyy, gs_x_yyyyy_xyz, gs_x_yyyyy_xzz, gs_x_yyyyy_yyy, gs_x_yyyyy_yyz, gs_x_yyyyy_yzz, gs_x_yyyyy_zzz, ts_yyyyy_xx, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xz, ts_yyyyy_xzz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyy_zz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyy_xxx[i] = 6.0 * ts_yyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxy[i] = 4.0 * ts_yyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxz[i] = 4.0 * ts_yyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyy[i] = 2.0 * ts_yyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyz[i] = 2.0 * ts_yyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xzz[i] = 2.0 * ts_yyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyy[i] = 2.0 * ts_yyyyy_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyz[i] = 2.0 * ts_yyyyy_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yzz[i] = 2.0 * ts_yyyyy_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_zzz[i] = 2.0 * ts_yyyyy_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 160-170 components of targeted buffer : HF

    auto gs_x_yyyyz_xxx = pbuffer.data(idx_g_hf + 160);

    auto gs_x_yyyyz_xxy = pbuffer.data(idx_g_hf + 161);

    auto gs_x_yyyyz_xxz = pbuffer.data(idx_g_hf + 162);

    auto gs_x_yyyyz_xyy = pbuffer.data(idx_g_hf + 163);

    auto gs_x_yyyyz_xyz = pbuffer.data(idx_g_hf + 164);

    auto gs_x_yyyyz_xzz = pbuffer.data(idx_g_hf + 165);

    auto gs_x_yyyyz_yyy = pbuffer.data(idx_g_hf + 166);

    auto gs_x_yyyyz_yyz = pbuffer.data(idx_g_hf + 167);

    auto gs_x_yyyyz_yzz = pbuffer.data(idx_g_hf + 168);

    auto gs_x_yyyyz_zzz = pbuffer.data(idx_g_hf + 169);

    #pragma omp simd aligned(gc_x, gs_x_yyyyz_xxx, gs_x_yyyyz_xxy, gs_x_yyyyz_xxz, gs_x_yyyyz_xyy, gs_x_yyyyz_xyz, gs_x_yyyyz_xzz, gs_x_yyyyz_yyy, gs_x_yyyyz_yyz, gs_x_yyyyz_yzz, gs_x_yyyyz_zzz, ts_yyyyz_xx, ts_yyyyz_xxx, ts_yyyyz_xxy, ts_yyyyz_xxz, ts_yyyyz_xy, ts_yyyyz_xyy, ts_yyyyz_xyz, ts_yyyyz_xz, ts_yyyyz_xzz, ts_yyyyz_yy, ts_yyyyz_yyy, ts_yyyyz_yyz, ts_yyyyz_yz, ts_yyyyz_yzz, ts_yyyyz_zz, ts_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyz_xxx[i] = 6.0 * ts_yyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxy[i] = 4.0 * ts_yyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxz[i] = 4.0 * ts_yyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyy[i] = 2.0 * ts_yyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyz[i] = 2.0 * ts_yyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xzz[i] = 2.0 * ts_yyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyy[i] = 2.0 * ts_yyyyz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyz[i] = 2.0 * ts_yyyyz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yzz[i] = 2.0 * ts_yyyyz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_zzz[i] = 2.0 * ts_yyyyz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 170-180 components of targeted buffer : HF

    auto gs_x_yyyzz_xxx = pbuffer.data(idx_g_hf + 170);

    auto gs_x_yyyzz_xxy = pbuffer.data(idx_g_hf + 171);

    auto gs_x_yyyzz_xxz = pbuffer.data(idx_g_hf + 172);

    auto gs_x_yyyzz_xyy = pbuffer.data(idx_g_hf + 173);

    auto gs_x_yyyzz_xyz = pbuffer.data(idx_g_hf + 174);

    auto gs_x_yyyzz_xzz = pbuffer.data(idx_g_hf + 175);

    auto gs_x_yyyzz_yyy = pbuffer.data(idx_g_hf + 176);

    auto gs_x_yyyzz_yyz = pbuffer.data(idx_g_hf + 177);

    auto gs_x_yyyzz_yzz = pbuffer.data(idx_g_hf + 178);

    auto gs_x_yyyzz_zzz = pbuffer.data(idx_g_hf + 179);

    #pragma omp simd aligned(gc_x, gs_x_yyyzz_xxx, gs_x_yyyzz_xxy, gs_x_yyyzz_xxz, gs_x_yyyzz_xyy, gs_x_yyyzz_xyz, gs_x_yyyzz_xzz, gs_x_yyyzz_yyy, gs_x_yyyzz_yyz, gs_x_yyyzz_yzz, gs_x_yyyzz_zzz, ts_yyyzz_xx, ts_yyyzz_xxx, ts_yyyzz_xxy, ts_yyyzz_xxz, ts_yyyzz_xy, ts_yyyzz_xyy, ts_yyyzz_xyz, ts_yyyzz_xz, ts_yyyzz_xzz, ts_yyyzz_yy, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yz, ts_yyyzz_yzz, ts_yyyzz_zz, ts_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyzz_xxx[i] = 6.0 * ts_yyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxy[i] = 4.0 * ts_yyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxz[i] = 4.0 * ts_yyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyy[i] = 2.0 * ts_yyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyz[i] = 2.0 * ts_yyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xzz[i] = 2.0 * ts_yyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyy[i] = 2.0 * ts_yyyzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyz[i] = 2.0 * ts_yyyzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yzz[i] = 2.0 * ts_yyyzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_zzz[i] = 2.0 * ts_yyyzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 180-190 components of targeted buffer : HF

    auto gs_x_yyzzz_xxx = pbuffer.data(idx_g_hf + 180);

    auto gs_x_yyzzz_xxy = pbuffer.data(idx_g_hf + 181);

    auto gs_x_yyzzz_xxz = pbuffer.data(idx_g_hf + 182);

    auto gs_x_yyzzz_xyy = pbuffer.data(idx_g_hf + 183);

    auto gs_x_yyzzz_xyz = pbuffer.data(idx_g_hf + 184);

    auto gs_x_yyzzz_xzz = pbuffer.data(idx_g_hf + 185);

    auto gs_x_yyzzz_yyy = pbuffer.data(idx_g_hf + 186);

    auto gs_x_yyzzz_yyz = pbuffer.data(idx_g_hf + 187);

    auto gs_x_yyzzz_yzz = pbuffer.data(idx_g_hf + 188);

    auto gs_x_yyzzz_zzz = pbuffer.data(idx_g_hf + 189);

    #pragma omp simd aligned(gc_x, gs_x_yyzzz_xxx, gs_x_yyzzz_xxy, gs_x_yyzzz_xxz, gs_x_yyzzz_xyy, gs_x_yyzzz_xyz, gs_x_yyzzz_xzz, gs_x_yyzzz_yyy, gs_x_yyzzz_yyz, gs_x_yyzzz_yzz, gs_x_yyzzz_zzz, ts_yyzzz_xx, ts_yyzzz_xxx, ts_yyzzz_xxy, ts_yyzzz_xxz, ts_yyzzz_xy, ts_yyzzz_xyy, ts_yyzzz_xyz, ts_yyzzz_xz, ts_yyzzz_xzz, ts_yyzzz_yy, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yz, ts_yyzzz_yzz, ts_yyzzz_zz, ts_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzzz_xxx[i] = 6.0 * ts_yyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxy[i] = 4.0 * ts_yyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxz[i] = 4.0 * ts_yyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyy[i] = 2.0 * ts_yyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyz[i] = 2.0 * ts_yyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xzz[i] = 2.0 * ts_yyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyy[i] = 2.0 * ts_yyzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyz[i] = 2.0 * ts_yyzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yzz[i] = 2.0 * ts_yyzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_zzz[i] = 2.0 * ts_yyzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 190-200 components of targeted buffer : HF

    auto gs_x_yzzzz_xxx = pbuffer.data(idx_g_hf + 190);

    auto gs_x_yzzzz_xxy = pbuffer.data(idx_g_hf + 191);

    auto gs_x_yzzzz_xxz = pbuffer.data(idx_g_hf + 192);

    auto gs_x_yzzzz_xyy = pbuffer.data(idx_g_hf + 193);

    auto gs_x_yzzzz_xyz = pbuffer.data(idx_g_hf + 194);

    auto gs_x_yzzzz_xzz = pbuffer.data(idx_g_hf + 195);

    auto gs_x_yzzzz_yyy = pbuffer.data(idx_g_hf + 196);

    auto gs_x_yzzzz_yyz = pbuffer.data(idx_g_hf + 197);

    auto gs_x_yzzzz_yzz = pbuffer.data(idx_g_hf + 198);

    auto gs_x_yzzzz_zzz = pbuffer.data(idx_g_hf + 199);

    #pragma omp simd aligned(gc_x, gs_x_yzzzz_xxx, gs_x_yzzzz_xxy, gs_x_yzzzz_xxz, gs_x_yzzzz_xyy, gs_x_yzzzz_xyz, gs_x_yzzzz_xzz, gs_x_yzzzz_yyy, gs_x_yzzzz_yyz, gs_x_yzzzz_yzz, gs_x_yzzzz_zzz, ts_yzzzz_xx, ts_yzzzz_xxx, ts_yzzzz_xxy, ts_yzzzz_xxz, ts_yzzzz_xy, ts_yzzzz_xyy, ts_yzzzz_xyz, ts_yzzzz_xz, ts_yzzzz_xzz, ts_yzzzz_yy, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yz, ts_yzzzz_yzz, ts_yzzzz_zz, ts_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzzz_xxx[i] = 6.0 * ts_yzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxy[i] = 4.0 * ts_yzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxz[i] = 4.0 * ts_yzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyy[i] = 2.0 * ts_yzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyz[i] = 2.0 * ts_yzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xzz[i] = 2.0 * ts_yzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyy[i] = 2.0 * ts_yzzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyz[i] = 2.0 * ts_yzzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yzz[i] = 2.0 * ts_yzzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_zzz[i] = 2.0 * ts_yzzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 200-210 components of targeted buffer : HF

    auto gs_x_zzzzz_xxx = pbuffer.data(idx_g_hf + 200);

    auto gs_x_zzzzz_xxy = pbuffer.data(idx_g_hf + 201);

    auto gs_x_zzzzz_xxz = pbuffer.data(idx_g_hf + 202);

    auto gs_x_zzzzz_xyy = pbuffer.data(idx_g_hf + 203);

    auto gs_x_zzzzz_xyz = pbuffer.data(idx_g_hf + 204);

    auto gs_x_zzzzz_xzz = pbuffer.data(idx_g_hf + 205);

    auto gs_x_zzzzz_yyy = pbuffer.data(idx_g_hf + 206);

    auto gs_x_zzzzz_yyz = pbuffer.data(idx_g_hf + 207);

    auto gs_x_zzzzz_yzz = pbuffer.data(idx_g_hf + 208);

    auto gs_x_zzzzz_zzz = pbuffer.data(idx_g_hf + 209);

    #pragma omp simd aligned(gc_x, gs_x_zzzzz_xxx, gs_x_zzzzz_xxy, gs_x_zzzzz_xxz, gs_x_zzzzz_xyy, gs_x_zzzzz_xyz, gs_x_zzzzz_xzz, gs_x_zzzzz_yyy, gs_x_zzzzz_yyz, gs_x_zzzzz_yzz, gs_x_zzzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzzz_xxx[i] = 6.0 * ts_zzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxx[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxy[i] = 4.0 * ts_zzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxz[i] = 4.0 * ts_zzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyy[i] = 2.0 * ts_zzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyz[i] = 2.0 * ts_zzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xzz[i] = 2.0 * ts_zzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyy[i] = 2.0 * ts_zzzzz_yyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyz[i] = 2.0 * ts_zzzzz_yyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yzz[i] = 2.0 * ts_zzzzz_yzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_zzz[i] = 2.0 * ts_zzzzz_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-220 components of targeted buffer : HF

    auto gs_y_xxxxx_xxx = pbuffer.data(idx_g_hf + 210);

    auto gs_y_xxxxx_xxy = pbuffer.data(idx_g_hf + 211);

    auto gs_y_xxxxx_xxz = pbuffer.data(idx_g_hf + 212);

    auto gs_y_xxxxx_xyy = pbuffer.data(idx_g_hf + 213);

    auto gs_y_xxxxx_xyz = pbuffer.data(idx_g_hf + 214);

    auto gs_y_xxxxx_xzz = pbuffer.data(idx_g_hf + 215);

    auto gs_y_xxxxx_yyy = pbuffer.data(idx_g_hf + 216);

    auto gs_y_xxxxx_yyz = pbuffer.data(idx_g_hf + 217);

    auto gs_y_xxxxx_yzz = pbuffer.data(idx_g_hf + 218);

    auto gs_y_xxxxx_zzz = pbuffer.data(idx_g_hf + 219);

    #pragma omp simd aligned(gc_y, gs_y_xxxxx_xxx, gs_y_xxxxx_xxy, gs_y_xxxxx_xxz, gs_y_xxxxx_xyy, gs_y_xxxxx_xyz, gs_y_xxxxx_xzz, gs_y_xxxxx_yyy, gs_y_xxxxx_yyz, gs_y_xxxxx_yzz, gs_y_xxxxx_zzz, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_yy, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yz, ts_xxxxx_yzz, ts_xxxxx_zz, ts_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxx_xxx[i] = 2.0 * ts_xxxxx_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxy[i] = 2.0 * ts_xxxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxz[i] = 2.0 * ts_xxxxx_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyy[i] = 4.0 * ts_xxxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyz[i] = 2.0 * ts_xxxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xzz[i] = 2.0 * ts_xxxxx_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyy[i] = 6.0 * ts_xxxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyz[i] = 4.0 * ts_xxxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yzz[i] = 2.0 * ts_xxxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_zzz[i] = 2.0 * ts_xxxxx_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 220-230 components of targeted buffer : HF

    auto gs_y_xxxxy_xxx = pbuffer.data(idx_g_hf + 220);

    auto gs_y_xxxxy_xxy = pbuffer.data(idx_g_hf + 221);

    auto gs_y_xxxxy_xxz = pbuffer.data(idx_g_hf + 222);

    auto gs_y_xxxxy_xyy = pbuffer.data(idx_g_hf + 223);

    auto gs_y_xxxxy_xyz = pbuffer.data(idx_g_hf + 224);

    auto gs_y_xxxxy_xzz = pbuffer.data(idx_g_hf + 225);

    auto gs_y_xxxxy_yyy = pbuffer.data(idx_g_hf + 226);

    auto gs_y_xxxxy_yyz = pbuffer.data(idx_g_hf + 227);

    auto gs_y_xxxxy_yzz = pbuffer.data(idx_g_hf + 228);

    auto gs_y_xxxxy_zzz = pbuffer.data(idx_g_hf + 229);

    #pragma omp simd aligned(gc_y, gs_y_xxxxy_xxx, gs_y_xxxxy_xxy, gs_y_xxxxy_xxz, gs_y_xxxxy_xyy, gs_y_xxxxy_xyz, gs_y_xxxxy_xzz, gs_y_xxxxy_yyy, gs_y_xxxxy_yyz, gs_y_xxxxy_yzz, gs_y_xxxxy_zzz, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xzz, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yzz, ts_xxxx_zzz, ts_xxxxy_xx, ts_xxxxy_xxx, ts_xxxxy_xxy, ts_xxxxy_xxz, ts_xxxxy_xy, ts_xxxxy_xyy, ts_xxxxy_xyz, ts_xxxxy_xz, ts_xxxxy_xzz, ts_xxxxy_yy, ts_xxxxy_yyy, ts_xxxxy_yyz, ts_xxxxy_yz, ts_xxxxy_yzz, ts_xxxxy_zz, ts_xxxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxy_xxx[i] = 2.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxy[i] = 2.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxz[i] = 2.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyy[i] = 2.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyz[i] = 2.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xzz[i] = 2.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyy[i] = 2.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyz[i] = 2.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yzz[i] = 2.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_zzz[i] = 2.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 230-240 components of targeted buffer : HF

    auto gs_y_xxxxz_xxx = pbuffer.data(idx_g_hf + 230);

    auto gs_y_xxxxz_xxy = pbuffer.data(idx_g_hf + 231);

    auto gs_y_xxxxz_xxz = pbuffer.data(idx_g_hf + 232);

    auto gs_y_xxxxz_xyy = pbuffer.data(idx_g_hf + 233);

    auto gs_y_xxxxz_xyz = pbuffer.data(idx_g_hf + 234);

    auto gs_y_xxxxz_xzz = pbuffer.data(idx_g_hf + 235);

    auto gs_y_xxxxz_yyy = pbuffer.data(idx_g_hf + 236);

    auto gs_y_xxxxz_yyz = pbuffer.data(idx_g_hf + 237);

    auto gs_y_xxxxz_yzz = pbuffer.data(idx_g_hf + 238);

    auto gs_y_xxxxz_zzz = pbuffer.data(idx_g_hf + 239);

    #pragma omp simd aligned(gc_y, gs_y_xxxxz_xxx, gs_y_xxxxz_xxy, gs_y_xxxxz_xxz, gs_y_xxxxz_xyy, gs_y_xxxxz_xyz, gs_y_xxxxz_xzz, gs_y_xxxxz_yyy, gs_y_xxxxz_yyz, gs_y_xxxxz_yzz, gs_y_xxxxz_zzz, ts_xxxxz_xx, ts_xxxxz_xxx, ts_xxxxz_xxy, ts_xxxxz_xxz, ts_xxxxz_xy, ts_xxxxz_xyy, ts_xxxxz_xyz, ts_xxxxz_xz, ts_xxxxz_xzz, ts_xxxxz_yy, ts_xxxxz_yyy, ts_xxxxz_yyz, ts_xxxxz_yz, ts_xxxxz_yzz, ts_xxxxz_zz, ts_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxz_xxx[i] = 2.0 * ts_xxxxz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxy[i] = 2.0 * ts_xxxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxz[i] = 2.0 * ts_xxxxz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyy[i] = 4.0 * ts_xxxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyz[i] = 2.0 * ts_xxxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xzz[i] = 2.0 * ts_xxxxz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyy[i] = 6.0 * ts_xxxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyz[i] = 4.0 * ts_xxxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yzz[i] = 2.0 * ts_xxxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_zzz[i] = 2.0 * ts_xxxxz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 240-250 components of targeted buffer : HF

    auto gs_y_xxxyy_xxx = pbuffer.data(idx_g_hf + 240);

    auto gs_y_xxxyy_xxy = pbuffer.data(idx_g_hf + 241);

    auto gs_y_xxxyy_xxz = pbuffer.data(idx_g_hf + 242);

    auto gs_y_xxxyy_xyy = pbuffer.data(idx_g_hf + 243);

    auto gs_y_xxxyy_xyz = pbuffer.data(idx_g_hf + 244);

    auto gs_y_xxxyy_xzz = pbuffer.data(idx_g_hf + 245);

    auto gs_y_xxxyy_yyy = pbuffer.data(idx_g_hf + 246);

    auto gs_y_xxxyy_yyz = pbuffer.data(idx_g_hf + 247);

    auto gs_y_xxxyy_yzz = pbuffer.data(idx_g_hf + 248);

    auto gs_y_xxxyy_zzz = pbuffer.data(idx_g_hf + 249);

    #pragma omp simd aligned(gc_y, gs_y_xxxyy_xxx, gs_y_xxxyy_xxy, gs_y_xxxyy_xxz, gs_y_xxxyy_xyy, gs_y_xxxyy_xyz, gs_y_xxxyy_xzz, gs_y_xxxyy_yyy, gs_y_xxxyy_yyz, gs_y_xxxyy_yzz, gs_y_xxxyy_zzz, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xzz, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yzz, ts_xxxy_zzz, ts_xxxyy_xx, ts_xxxyy_xxx, ts_xxxyy_xxy, ts_xxxyy_xxz, ts_xxxyy_xy, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_xz, ts_xxxyy_xzz, ts_xxxyy_yy, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yz, ts_xxxyy_yzz, ts_xxxyy_zz, ts_xxxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyy_xxx[i] = 4.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxy[i] = 4.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxz[i] = 4.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyy[i] = 4.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyz[i] = 4.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xzz[i] = 4.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyy[i] = 4.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyz[i] = 4.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yzz[i] = 4.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_zzz[i] = 4.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 250-260 components of targeted buffer : HF

    auto gs_y_xxxyz_xxx = pbuffer.data(idx_g_hf + 250);

    auto gs_y_xxxyz_xxy = pbuffer.data(idx_g_hf + 251);

    auto gs_y_xxxyz_xxz = pbuffer.data(idx_g_hf + 252);

    auto gs_y_xxxyz_xyy = pbuffer.data(idx_g_hf + 253);

    auto gs_y_xxxyz_xyz = pbuffer.data(idx_g_hf + 254);

    auto gs_y_xxxyz_xzz = pbuffer.data(idx_g_hf + 255);

    auto gs_y_xxxyz_yyy = pbuffer.data(idx_g_hf + 256);

    auto gs_y_xxxyz_yyz = pbuffer.data(idx_g_hf + 257);

    auto gs_y_xxxyz_yzz = pbuffer.data(idx_g_hf + 258);

    auto gs_y_xxxyz_zzz = pbuffer.data(idx_g_hf + 259);

    #pragma omp simd aligned(gc_y, gs_y_xxxyz_xxx, gs_y_xxxyz_xxy, gs_y_xxxyz_xxz, gs_y_xxxyz_xyy, gs_y_xxxyz_xyz, gs_y_xxxyz_xzz, gs_y_xxxyz_yyy, gs_y_xxxyz_yyz, gs_y_xxxyz_yzz, gs_y_xxxyz_zzz, ts_xxxyz_xx, ts_xxxyz_xxx, ts_xxxyz_xxy, ts_xxxyz_xxz, ts_xxxyz_xy, ts_xxxyz_xyy, ts_xxxyz_xyz, ts_xxxyz_xz, ts_xxxyz_xzz, ts_xxxyz_yy, ts_xxxyz_yyy, ts_xxxyz_yyz, ts_xxxyz_yz, ts_xxxyz_yzz, ts_xxxyz_zz, ts_xxxyz_zzz, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xzz, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yzz, ts_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyz_xxx[i] = 2.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxy[i] = 2.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxz[i] = 2.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyy[i] = 2.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyz[i] = 2.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xzz[i] = 2.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyy[i] = 2.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyz[i] = 2.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yzz[i] = 2.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_zzz[i] = 2.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 260-270 components of targeted buffer : HF

    auto gs_y_xxxzz_xxx = pbuffer.data(idx_g_hf + 260);

    auto gs_y_xxxzz_xxy = pbuffer.data(idx_g_hf + 261);

    auto gs_y_xxxzz_xxz = pbuffer.data(idx_g_hf + 262);

    auto gs_y_xxxzz_xyy = pbuffer.data(idx_g_hf + 263);

    auto gs_y_xxxzz_xyz = pbuffer.data(idx_g_hf + 264);

    auto gs_y_xxxzz_xzz = pbuffer.data(idx_g_hf + 265);

    auto gs_y_xxxzz_yyy = pbuffer.data(idx_g_hf + 266);

    auto gs_y_xxxzz_yyz = pbuffer.data(idx_g_hf + 267);

    auto gs_y_xxxzz_yzz = pbuffer.data(idx_g_hf + 268);

    auto gs_y_xxxzz_zzz = pbuffer.data(idx_g_hf + 269);

    #pragma omp simd aligned(gc_y, gs_y_xxxzz_xxx, gs_y_xxxzz_xxy, gs_y_xxxzz_xxz, gs_y_xxxzz_xyy, gs_y_xxxzz_xyz, gs_y_xxxzz_xzz, gs_y_xxxzz_yyy, gs_y_xxxzz_yyz, gs_y_xxxzz_yzz, gs_y_xxxzz_zzz, ts_xxxzz_xx, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xxz, ts_xxxzz_xy, ts_xxxzz_xyy, ts_xxxzz_xyz, ts_xxxzz_xz, ts_xxxzz_xzz, ts_xxxzz_yy, ts_xxxzz_yyy, ts_xxxzz_yyz, ts_xxxzz_yz, ts_xxxzz_yzz, ts_xxxzz_zz, ts_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxzz_xxx[i] = 2.0 * ts_xxxzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxy[i] = 2.0 * ts_xxxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxz[i] = 2.0 * ts_xxxzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyy[i] = 4.0 * ts_xxxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyz[i] = 2.0 * ts_xxxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xzz[i] = 2.0 * ts_xxxzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyy[i] = 6.0 * ts_xxxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyz[i] = 4.0 * ts_xxxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yzz[i] = 2.0 * ts_xxxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_zzz[i] = 2.0 * ts_xxxzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 270-280 components of targeted buffer : HF

    auto gs_y_xxyyy_xxx = pbuffer.data(idx_g_hf + 270);

    auto gs_y_xxyyy_xxy = pbuffer.data(idx_g_hf + 271);

    auto gs_y_xxyyy_xxz = pbuffer.data(idx_g_hf + 272);

    auto gs_y_xxyyy_xyy = pbuffer.data(idx_g_hf + 273);

    auto gs_y_xxyyy_xyz = pbuffer.data(idx_g_hf + 274);

    auto gs_y_xxyyy_xzz = pbuffer.data(idx_g_hf + 275);

    auto gs_y_xxyyy_yyy = pbuffer.data(idx_g_hf + 276);

    auto gs_y_xxyyy_yyz = pbuffer.data(idx_g_hf + 277);

    auto gs_y_xxyyy_yzz = pbuffer.data(idx_g_hf + 278);

    auto gs_y_xxyyy_zzz = pbuffer.data(idx_g_hf + 279);

    #pragma omp simd aligned(gc_y, gs_y_xxyyy_xxx, gs_y_xxyyy_xxy, gs_y_xxyyy_xxz, gs_y_xxyyy_xyy, gs_y_xxyyy_xyz, gs_y_xxyyy_xzz, gs_y_xxyyy_yyy, gs_y_xxyyy_yyz, gs_y_xxyyy_yzz, gs_y_xxyyy_zzz, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xzz, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yzz, ts_xxyy_zzz, ts_xxyyy_xx, ts_xxyyy_xxx, ts_xxyyy_xxy, ts_xxyyy_xxz, ts_xxyyy_xy, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_xz, ts_xxyyy_xzz, ts_xxyyy_yy, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yz, ts_xxyyy_yzz, ts_xxyyy_zz, ts_xxyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyy_xxx[i] = 6.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxy[i] = 6.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxz[i] = 6.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyy[i] = 6.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyz[i] = 6.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xzz[i] = 6.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyy[i] = 6.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyz[i] = 6.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yzz[i] = 6.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_zzz[i] = 6.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 280-290 components of targeted buffer : HF

    auto gs_y_xxyyz_xxx = pbuffer.data(idx_g_hf + 280);

    auto gs_y_xxyyz_xxy = pbuffer.data(idx_g_hf + 281);

    auto gs_y_xxyyz_xxz = pbuffer.data(idx_g_hf + 282);

    auto gs_y_xxyyz_xyy = pbuffer.data(idx_g_hf + 283);

    auto gs_y_xxyyz_xyz = pbuffer.data(idx_g_hf + 284);

    auto gs_y_xxyyz_xzz = pbuffer.data(idx_g_hf + 285);

    auto gs_y_xxyyz_yyy = pbuffer.data(idx_g_hf + 286);

    auto gs_y_xxyyz_yyz = pbuffer.data(idx_g_hf + 287);

    auto gs_y_xxyyz_yzz = pbuffer.data(idx_g_hf + 288);

    auto gs_y_xxyyz_zzz = pbuffer.data(idx_g_hf + 289);

    #pragma omp simd aligned(gc_y, gs_y_xxyyz_xxx, gs_y_xxyyz_xxy, gs_y_xxyyz_xxz, gs_y_xxyyz_xyy, gs_y_xxyyz_xyz, gs_y_xxyyz_xzz, gs_y_xxyyz_yyy, gs_y_xxyyz_yyz, gs_y_xxyyz_yzz, gs_y_xxyyz_zzz, ts_xxyyz_xx, ts_xxyyz_xxx, ts_xxyyz_xxy, ts_xxyyz_xxz, ts_xxyyz_xy, ts_xxyyz_xyy, ts_xxyyz_xyz, ts_xxyyz_xz, ts_xxyyz_xzz, ts_xxyyz_yy, ts_xxyyz_yyy, ts_xxyyz_yyz, ts_xxyyz_yz, ts_xxyyz_yzz, ts_xxyyz_zz, ts_xxyyz_zzz, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xzz, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yzz, ts_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyz_xxx[i] = 4.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxy[i] = 4.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxz[i] = 4.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyy[i] = 4.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyz[i] = 4.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xzz[i] = 4.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyy[i] = 4.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyz[i] = 4.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yzz[i] = 4.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_zzz[i] = 4.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 290-300 components of targeted buffer : HF

    auto gs_y_xxyzz_xxx = pbuffer.data(idx_g_hf + 290);

    auto gs_y_xxyzz_xxy = pbuffer.data(idx_g_hf + 291);

    auto gs_y_xxyzz_xxz = pbuffer.data(idx_g_hf + 292);

    auto gs_y_xxyzz_xyy = pbuffer.data(idx_g_hf + 293);

    auto gs_y_xxyzz_xyz = pbuffer.data(idx_g_hf + 294);

    auto gs_y_xxyzz_xzz = pbuffer.data(idx_g_hf + 295);

    auto gs_y_xxyzz_yyy = pbuffer.data(idx_g_hf + 296);

    auto gs_y_xxyzz_yyz = pbuffer.data(idx_g_hf + 297);

    auto gs_y_xxyzz_yzz = pbuffer.data(idx_g_hf + 298);

    auto gs_y_xxyzz_zzz = pbuffer.data(idx_g_hf + 299);

    #pragma omp simd aligned(gc_y, gs_y_xxyzz_xxx, gs_y_xxyzz_xxy, gs_y_xxyzz_xxz, gs_y_xxyzz_xyy, gs_y_xxyzz_xyz, gs_y_xxyzz_xzz, gs_y_xxyzz_yyy, gs_y_xxyzz_yyz, gs_y_xxyzz_yzz, gs_y_xxyzz_zzz, ts_xxyzz_xx, ts_xxyzz_xxx, ts_xxyzz_xxy, ts_xxyzz_xxz, ts_xxyzz_xy, ts_xxyzz_xyy, ts_xxyzz_xyz, ts_xxyzz_xz, ts_xxyzz_xzz, ts_xxyzz_yy, ts_xxyzz_yyy, ts_xxyzz_yyz, ts_xxyzz_yz, ts_xxyzz_yzz, ts_xxyzz_zz, ts_xxyzz_zzz, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xzz, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yzz, ts_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyzz_xxx[i] = 2.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxy[i] = 2.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxz[i] = 2.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyy[i] = 2.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyz[i] = 2.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xzz[i] = 2.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyy[i] = 2.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyz[i] = 2.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yzz[i] = 2.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_zzz[i] = 2.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 300-310 components of targeted buffer : HF

    auto gs_y_xxzzz_xxx = pbuffer.data(idx_g_hf + 300);

    auto gs_y_xxzzz_xxy = pbuffer.data(idx_g_hf + 301);

    auto gs_y_xxzzz_xxz = pbuffer.data(idx_g_hf + 302);

    auto gs_y_xxzzz_xyy = pbuffer.data(idx_g_hf + 303);

    auto gs_y_xxzzz_xyz = pbuffer.data(idx_g_hf + 304);

    auto gs_y_xxzzz_xzz = pbuffer.data(idx_g_hf + 305);

    auto gs_y_xxzzz_yyy = pbuffer.data(idx_g_hf + 306);

    auto gs_y_xxzzz_yyz = pbuffer.data(idx_g_hf + 307);

    auto gs_y_xxzzz_yzz = pbuffer.data(idx_g_hf + 308);

    auto gs_y_xxzzz_zzz = pbuffer.data(idx_g_hf + 309);

    #pragma omp simd aligned(gc_y, gs_y_xxzzz_xxx, gs_y_xxzzz_xxy, gs_y_xxzzz_xxz, gs_y_xxzzz_xyy, gs_y_xxzzz_xyz, gs_y_xxzzz_xzz, gs_y_xxzzz_yyy, gs_y_xxzzz_yyz, gs_y_xxzzz_yzz, gs_y_xxzzz_zzz, ts_xxzzz_xx, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xxz, ts_xxzzz_xy, ts_xxzzz_xyy, ts_xxzzz_xyz, ts_xxzzz_xz, ts_xxzzz_xzz, ts_xxzzz_yy, ts_xxzzz_yyy, ts_xxzzz_yyz, ts_xxzzz_yz, ts_xxzzz_yzz, ts_xxzzz_zz, ts_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzzz_xxx[i] = 2.0 * ts_xxzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxy[i] = 2.0 * ts_xxzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxz[i] = 2.0 * ts_xxzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyy[i] = 4.0 * ts_xxzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyz[i] = 2.0 * ts_xxzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xzz[i] = 2.0 * ts_xxzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyy[i] = 6.0 * ts_xxzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyz[i] = 4.0 * ts_xxzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yzz[i] = 2.0 * ts_xxzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_zzz[i] = 2.0 * ts_xxzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 310-320 components of targeted buffer : HF

    auto gs_y_xyyyy_xxx = pbuffer.data(idx_g_hf + 310);

    auto gs_y_xyyyy_xxy = pbuffer.data(idx_g_hf + 311);

    auto gs_y_xyyyy_xxz = pbuffer.data(idx_g_hf + 312);

    auto gs_y_xyyyy_xyy = pbuffer.data(idx_g_hf + 313);

    auto gs_y_xyyyy_xyz = pbuffer.data(idx_g_hf + 314);

    auto gs_y_xyyyy_xzz = pbuffer.data(idx_g_hf + 315);

    auto gs_y_xyyyy_yyy = pbuffer.data(idx_g_hf + 316);

    auto gs_y_xyyyy_yyz = pbuffer.data(idx_g_hf + 317);

    auto gs_y_xyyyy_yzz = pbuffer.data(idx_g_hf + 318);

    auto gs_y_xyyyy_zzz = pbuffer.data(idx_g_hf + 319);

    #pragma omp simd aligned(gc_y, gs_y_xyyyy_xxx, gs_y_xyyyy_xxy, gs_y_xyyyy_xxz, gs_y_xyyyy_xyy, gs_y_xyyyy_xyz, gs_y_xyyyy_xzz, gs_y_xyyyy_yyy, gs_y_xyyyy_yyz, gs_y_xyyyy_yzz, gs_y_xyyyy_zzz, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xzz, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yzz, ts_xyyy_zzz, ts_xyyyy_xx, ts_xyyyy_xxx, ts_xyyyy_xxy, ts_xyyyy_xxz, ts_xyyyy_xy, ts_xyyyy_xyy, ts_xyyyy_xyz, ts_xyyyy_xz, ts_xyyyy_xzz, ts_xyyyy_yy, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yz, ts_xyyyy_yzz, ts_xyyyy_zz, ts_xyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyy_xxx[i] = 8.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxy[i] = 8.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxz[i] = 8.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyy[i] = 8.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyz[i] = 8.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xzz[i] = 8.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyy[i] = 8.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyz[i] = 8.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yzz[i] = 8.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_zzz[i] = 8.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 320-330 components of targeted buffer : HF

    auto gs_y_xyyyz_xxx = pbuffer.data(idx_g_hf + 320);

    auto gs_y_xyyyz_xxy = pbuffer.data(idx_g_hf + 321);

    auto gs_y_xyyyz_xxz = pbuffer.data(idx_g_hf + 322);

    auto gs_y_xyyyz_xyy = pbuffer.data(idx_g_hf + 323);

    auto gs_y_xyyyz_xyz = pbuffer.data(idx_g_hf + 324);

    auto gs_y_xyyyz_xzz = pbuffer.data(idx_g_hf + 325);

    auto gs_y_xyyyz_yyy = pbuffer.data(idx_g_hf + 326);

    auto gs_y_xyyyz_yyz = pbuffer.data(idx_g_hf + 327);

    auto gs_y_xyyyz_yzz = pbuffer.data(idx_g_hf + 328);

    auto gs_y_xyyyz_zzz = pbuffer.data(idx_g_hf + 329);

    #pragma omp simd aligned(gc_y, gs_y_xyyyz_xxx, gs_y_xyyyz_xxy, gs_y_xyyyz_xxz, gs_y_xyyyz_xyy, gs_y_xyyyz_xyz, gs_y_xyyyz_xzz, gs_y_xyyyz_yyy, gs_y_xyyyz_yyz, gs_y_xyyyz_yzz, gs_y_xyyyz_zzz, ts_xyyyz_xx, ts_xyyyz_xxx, ts_xyyyz_xxy, ts_xyyyz_xxz, ts_xyyyz_xy, ts_xyyyz_xyy, ts_xyyyz_xyz, ts_xyyyz_xz, ts_xyyyz_xzz, ts_xyyyz_yy, ts_xyyyz_yyy, ts_xyyyz_yyz, ts_xyyyz_yz, ts_xyyyz_yzz, ts_xyyyz_zz, ts_xyyyz_zzz, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xzz, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yzz, ts_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyz_xxx[i] = 6.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxy[i] = 6.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxz[i] = 6.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyy[i] = 6.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyz[i] = 6.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xzz[i] = 6.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyy[i] = 6.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyz[i] = 6.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yzz[i] = 6.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_zzz[i] = 6.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 330-340 components of targeted buffer : HF

    auto gs_y_xyyzz_xxx = pbuffer.data(idx_g_hf + 330);

    auto gs_y_xyyzz_xxy = pbuffer.data(idx_g_hf + 331);

    auto gs_y_xyyzz_xxz = pbuffer.data(idx_g_hf + 332);

    auto gs_y_xyyzz_xyy = pbuffer.data(idx_g_hf + 333);

    auto gs_y_xyyzz_xyz = pbuffer.data(idx_g_hf + 334);

    auto gs_y_xyyzz_xzz = pbuffer.data(idx_g_hf + 335);

    auto gs_y_xyyzz_yyy = pbuffer.data(idx_g_hf + 336);

    auto gs_y_xyyzz_yyz = pbuffer.data(idx_g_hf + 337);

    auto gs_y_xyyzz_yzz = pbuffer.data(idx_g_hf + 338);

    auto gs_y_xyyzz_zzz = pbuffer.data(idx_g_hf + 339);

    #pragma omp simd aligned(gc_y, gs_y_xyyzz_xxx, gs_y_xyyzz_xxy, gs_y_xyyzz_xxz, gs_y_xyyzz_xyy, gs_y_xyyzz_xyz, gs_y_xyyzz_xzz, gs_y_xyyzz_yyy, gs_y_xyyzz_yyz, gs_y_xyyzz_yzz, gs_y_xyyzz_zzz, ts_xyyzz_xx, ts_xyyzz_xxx, ts_xyyzz_xxy, ts_xyyzz_xxz, ts_xyyzz_xy, ts_xyyzz_xyy, ts_xyyzz_xyz, ts_xyyzz_xz, ts_xyyzz_xzz, ts_xyyzz_yy, ts_xyyzz_yyy, ts_xyyzz_yyz, ts_xyyzz_yz, ts_xyyzz_yzz, ts_xyyzz_zz, ts_xyyzz_zzz, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xzz, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yzz, ts_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyzz_xxx[i] = 4.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxy[i] = 4.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxz[i] = 4.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyy[i] = 4.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyz[i] = 4.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xzz[i] = 4.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyy[i] = 4.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyz[i] = 4.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yzz[i] = 4.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_zzz[i] = 4.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 340-350 components of targeted buffer : HF

    auto gs_y_xyzzz_xxx = pbuffer.data(idx_g_hf + 340);

    auto gs_y_xyzzz_xxy = pbuffer.data(idx_g_hf + 341);

    auto gs_y_xyzzz_xxz = pbuffer.data(idx_g_hf + 342);

    auto gs_y_xyzzz_xyy = pbuffer.data(idx_g_hf + 343);

    auto gs_y_xyzzz_xyz = pbuffer.data(idx_g_hf + 344);

    auto gs_y_xyzzz_xzz = pbuffer.data(idx_g_hf + 345);

    auto gs_y_xyzzz_yyy = pbuffer.data(idx_g_hf + 346);

    auto gs_y_xyzzz_yyz = pbuffer.data(idx_g_hf + 347);

    auto gs_y_xyzzz_yzz = pbuffer.data(idx_g_hf + 348);

    auto gs_y_xyzzz_zzz = pbuffer.data(idx_g_hf + 349);

    #pragma omp simd aligned(gc_y, gs_y_xyzzz_xxx, gs_y_xyzzz_xxy, gs_y_xyzzz_xxz, gs_y_xyzzz_xyy, gs_y_xyzzz_xyz, gs_y_xyzzz_xzz, gs_y_xyzzz_yyy, gs_y_xyzzz_yyz, gs_y_xyzzz_yzz, gs_y_xyzzz_zzz, ts_xyzzz_xx, ts_xyzzz_xxx, ts_xyzzz_xxy, ts_xyzzz_xxz, ts_xyzzz_xy, ts_xyzzz_xyy, ts_xyzzz_xyz, ts_xyzzz_xz, ts_xyzzz_xzz, ts_xyzzz_yy, ts_xyzzz_yyy, ts_xyzzz_yyz, ts_xyzzz_yz, ts_xyzzz_yzz, ts_xyzzz_zz, ts_xyzzz_zzz, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xzz, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yzz, ts_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzzz_xxx[i] = 2.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxy[i] = 2.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxz[i] = 2.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyy[i] = 2.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyz[i] = 2.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xzz[i] = 2.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyy[i] = 2.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyz[i] = 2.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yzz[i] = 2.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_zzz[i] = 2.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 350-360 components of targeted buffer : HF

    auto gs_y_xzzzz_xxx = pbuffer.data(idx_g_hf + 350);

    auto gs_y_xzzzz_xxy = pbuffer.data(idx_g_hf + 351);

    auto gs_y_xzzzz_xxz = pbuffer.data(idx_g_hf + 352);

    auto gs_y_xzzzz_xyy = pbuffer.data(idx_g_hf + 353);

    auto gs_y_xzzzz_xyz = pbuffer.data(idx_g_hf + 354);

    auto gs_y_xzzzz_xzz = pbuffer.data(idx_g_hf + 355);

    auto gs_y_xzzzz_yyy = pbuffer.data(idx_g_hf + 356);

    auto gs_y_xzzzz_yyz = pbuffer.data(idx_g_hf + 357);

    auto gs_y_xzzzz_yzz = pbuffer.data(idx_g_hf + 358);

    auto gs_y_xzzzz_zzz = pbuffer.data(idx_g_hf + 359);

    #pragma omp simd aligned(gc_y, gs_y_xzzzz_xxx, gs_y_xzzzz_xxy, gs_y_xzzzz_xxz, gs_y_xzzzz_xyy, gs_y_xzzzz_xyz, gs_y_xzzzz_xzz, gs_y_xzzzz_yyy, gs_y_xzzzz_yyz, gs_y_xzzzz_yzz, gs_y_xzzzz_zzz, ts_xzzzz_xx, ts_xzzzz_xxx, ts_xzzzz_xxy, ts_xzzzz_xxz, ts_xzzzz_xy, ts_xzzzz_xyy, ts_xzzzz_xyz, ts_xzzzz_xz, ts_xzzzz_xzz, ts_xzzzz_yy, ts_xzzzz_yyy, ts_xzzzz_yyz, ts_xzzzz_yz, ts_xzzzz_yzz, ts_xzzzz_zz, ts_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzzz_xxx[i] = 2.0 * ts_xzzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxy[i] = 2.0 * ts_xzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxz[i] = 2.0 * ts_xzzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyy[i] = 4.0 * ts_xzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyz[i] = 2.0 * ts_xzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xzz[i] = 2.0 * ts_xzzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyy[i] = 6.0 * ts_xzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyz[i] = 4.0 * ts_xzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yzz[i] = 2.0 * ts_xzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_zzz[i] = 2.0 * ts_xzzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 360-370 components of targeted buffer : HF

    auto gs_y_yyyyy_xxx = pbuffer.data(idx_g_hf + 360);

    auto gs_y_yyyyy_xxy = pbuffer.data(idx_g_hf + 361);

    auto gs_y_yyyyy_xxz = pbuffer.data(idx_g_hf + 362);

    auto gs_y_yyyyy_xyy = pbuffer.data(idx_g_hf + 363);

    auto gs_y_yyyyy_xyz = pbuffer.data(idx_g_hf + 364);

    auto gs_y_yyyyy_xzz = pbuffer.data(idx_g_hf + 365);

    auto gs_y_yyyyy_yyy = pbuffer.data(idx_g_hf + 366);

    auto gs_y_yyyyy_yyz = pbuffer.data(idx_g_hf + 367);

    auto gs_y_yyyyy_yzz = pbuffer.data(idx_g_hf + 368);

    auto gs_y_yyyyy_zzz = pbuffer.data(idx_g_hf + 369);

    #pragma omp simd aligned(gc_y, gs_y_yyyyy_xxx, gs_y_yyyyy_xxy, gs_y_yyyyy_xxz, gs_y_yyyyy_xyy, gs_y_yyyyy_xyz, gs_y_yyyyy_xzz, gs_y_yyyyy_yyy, gs_y_yyyyy_yyz, gs_y_yyyyy_yzz, gs_y_yyyyy_zzz, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xzz, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yzz, ts_yyyy_zzz, ts_yyyyy_xx, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xz, ts_yyyyy_xzz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyy_zz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyy_xxx[i] = 10.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxy[i] = 10.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxz[i] = 10.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyy[i] = 10.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyz[i] = 10.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xzz[i] = 10.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyy[i] = 10.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyz[i] = 10.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yzz[i] = 10.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_zzz[i] = 10.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 370-380 components of targeted buffer : HF

    auto gs_y_yyyyz_xxx = pbuffer.data(idx_g_hf + 370);

    auto gs_y_yyyyz_xxy = pbuffer.data(idx_g_hf + 371);

    auto gs_y_yyyyz_xxz = pbuffer.data(idx_g_hf + 372);

    auto gs_y_yyyyz_xyy = pbuffer.data(idx_g_hf + 373);

    auto gs_y_yyyyz_xyz = pbuffer.data(idx_g_hf + 374);

    auto gs_y_yyyyz_xzz = pbuffer.data(idx_g_hf + 375);

    auto gs_y_yyyyz_yyy = pbuffer.data(idx_g_hf + 376);

    auto gs_y_yyyyz_yyz = pbuffer.data(idx_g_hf + 377);

    auto gs_y_yyyyz_yzz = pbuffer.data(idx_g_hf + 378);

    auto gs_y_yyyyz_zzz = pbuffer.data(idx_g_hf + 379);

    #pragma omp simd aligned(gc_y, gs_y_yyyyz_xxx, gs_y_yyyyz_xxy, gs_y_yyyyz_xxz, gs_y_yyyyz_xyy, gs_y_yyyyz_xyz, gs_y_yyyyz_xzz, gs_y_yyyyz_yyy, gs_y_yyyyz_yyz, gs_y_yyyyz_yzz, gs_y_yyyyz_zzz, ts_yyyyz_xx, ts_yyyyz_xxx, ts_yyyyz_xxy, ts_yyyyz_xxz, ts_yyyyz_xy, ts_yyyyz_xyy, ts_yyyyz_xyz, ts_yyyyz_xz, ts_yyyyz_xzz, ts_yyyyz_yy, ts_yyyyz_yyy, ts_yyyyz_yyz, ts_yyyyz_yz, ts_yyyyz_yzz, ts_yyyyz_zz, ts_yyyyz_zzz, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xzz, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yzz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyz_xxx[i] = 8.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxy[i] = 8.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxz[i] = 8.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyy[i] = 8.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyz[i] = 8.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xzz[i] = 8.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyy[i] = 8.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyz[i] = 8.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yzz[i] = 8.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_zzz[i] = 8.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 380-390 components of targeted buffer : HF

    auto gs_y_yyyzz_xxx = pbuffer.data(idx_g_hf + 380);

    auto gs_y_yyyzz_xxy = pbuffer.data(idx_g_hf + 381);

    auto gs_y_yyyzz_xxz = pbuffer.data(idx_g_hf + 382);

    auto gs_y_yyyzz_xyy = pbuffer.data(idx_g_hf + 383);

    auto gs_y_yyyzz_xyz = pbuffer.data(idx_g_hf + 384);

    auto gs_y_yyyzz_xzz = pbuffer.data(idx_g_hf + 385);

    auto gs_y_yyyzz_yyy = pbuffer.data(idx_g_hf + 386);

    auto gs_y_yyyzz_yyz = pbuffer.data(idx_g_hf + 387);

    auto gs_y_yyyzz_yzz = pbuffer.data(idx_g_hf + 388);

    auto gs_y_yyyzz_zzz = pbuffer.data(idx_g_hf + 389);

    #pragma omp simd aligned(gc_y, gs_y_yyyzz_xxx, gs_y_yyyzz_xxy, gs_y_yyyzz_xxz, gs_y_yyyzz_xyy, gs_y_yyyzz_xyz, gs_y_yyyzz_xzz, gs_y_yyyzz_yyy, gs_y_yyyzz_yyz, gs_y_yyyzz_yzz, gs_y_yyyzz_zzz, ts_yyyzz_xx, ts_yyyzz_xxx, ts_yyyzz_xxy, ts_yyyzz_xxz, ts_yyyzz_xy, ts_yyyzz_xyy, ts_yyyzz_xyz, ts_yyyzz_xz, ts_yyyzz_xzz, ts_yyyzz_yy, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yz, ts_yyyzz_yzz, ts_yyyzz_zz, ts_yyyzz_zzz, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xzz, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yzz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyzz_xxx[i] = 6.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxy[i] = 6.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxz[i] = 6.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyy[i] = 6.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyz[i] = 6.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xzz[i] = 6.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyy[i] = 6.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyz[i] = 6.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yzz[i] = 6.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_zzz[i] = 6.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 390-400 components of targeted buffer : HF

    auto gs_y_yyzzz_xxx = pbuffer.data(idx_g_hf + 390);

    auto gs_y_yyzzz_xxy = pbuffer.data(idx_g_hf + 391);

    auto gs_y_yyzzz_xxz = pbuffer.data(idx_g_hf + 392);

    auto gs_y_yyzzz_xyy = pbuffer.data(idx_g_hf + 393);

    auto gs_y_yyzzz_xyz = pbuffer.data(idx_g_hf + 394);

    auto gs_y_yyzzz_xzz = pbuffer.data(idx_g_hf + 395);

    auto gs_y_yyzzz_yyy = pbuffer.data(idx_g_hf + 396);

    auto gs_y_yyzzz_yyz = pbuffer.data(idx_g_hf + 397);

    auto gs_y_yyzzz_yzz = pbuffer.data(idx_g_hf + 398);

    auto gs_y_yyzzz_zzz = pbuffer.data(idx_g_hf + 399);

    #pragma omp simd aligned(gc_y, gs_y_yyzzz_xxx, gs_y_yyzzz_xxy, gs_y_yyzzz_xxz, gs_y_yyzzz_xyy, gs_y_yyzzz_xyz, gs_y_yyzzz_xzz, gs_y_yyzzz_yyy, gs_y_yyzzz_yyz, gs_y_yyzzz_yzz, gs_y_yyzzz_zzz, ts_yyzzz_xx, ts_yyzzz_xxx, ts_yyzzz_xxy, ts_yyzzz_xxz, ts_yyzzz_xy, ts_yyzzz_xyy, ts_yyzzz_xyz, ts_yyzzz_xz, ts_yyzzz_xzz, ts_yyzzz_yy, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yz, ts_yyzzz_yzz, ts_yyzzz_zz, ts_yyzzz_zzz, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xzz, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yzz, ts_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzzz_xxx[i] = 4.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxy[i] = 4.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxz[i] = 4.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyy[i] = 4.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyz[i] = 4.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xzz[i] = 4.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyy[i] = 4.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyz[i] = 4.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yzz[i] = 4.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_zzz[i] = 4.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 400-410 components of targeted buffer : HF

    auto gs_y_yzzzz_xxx = pbuffer.data(idx_g_hf + 400);

    auto gs_y_yzzzz_xxy = pbuffer.data(idx_g_hf + 401);

    auto gs_y_yzzzz_xxz = pbuffer.data(idx_g_hf + 402);

    auto gs_y_yzzzz_xyy = pbuffer.data(idx_g_hf + 403);

    auto gs_y_yzzzz_xyz = pbuffer.data(idx_g_hf + 404);

    auto gs_y_yzzzz_xzz = pbuffer.data(idx_g_hf + 405);

    auto gs_y_yzzzz_yyy = pbuffer.data(idx_g_hf + 406);

    auto gs_y_yzzzz_yyz = pbuffer.data(idx_g_hf + 407);

    auto gs_y_yzzzz_yzz = pbuffer.data(idx_g_hf + 408);

    auto gs_y_yzzzz_zzz = pbuffer.data(idx_g_hf + 409);

    #pragma omp simd aligned(gc_y, gs_y_yzzzz_xxx, gs_y_yzzzz_xxy, gs_y_yzzzz_xxz, gs_y_yzzzz_xyy, gs_y_yzzzz_xyz, gs_y_yzzzz_xzz, gs_y_yzzzz_yyy, gs_y_yzzzz_yyz, gs_y_yzzzz_yzz, gs_y_yzzzz_zzz, ts_yzzzz_xx, ts_yzzzz_xxx, ts_yzzzz_xxy, ts_yzzzz_xxz, ts_yzzzz_xy, ts_yzzzz_xyy, ts_yzzzz_xyz, ts_yzzzz_xz, ts_yzzzz_xzz, ts_yzzzz_yy, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yz, ts_yzzzz_yzz, ts_yzzzz_zz, ts_yzzzz_zzz, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzzz_xxx[i] = 2.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxy[i] = 2.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxz[i] = 2.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyy[i] = 2.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyz[i] = 2.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xzz[i] = 2.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyy[i] = 2.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyz[i] = 2.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yzz[i] = 2.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_zzz[i] = 2.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 410-420 components of targeted buffer : HF

    auto gs_y_zzzzz_xxx = pbuffer.data(idx_g_hf + 410);

    auto gs_y_zzzzz_xxy = pbuffer.data(idx_g_hf + 411);

    auto gs_y_zzzzz_xxz = pbuffer.data(idx_g_hf + 412);

    auto gs_y_zzzzz_xyy = pbuffer.data(idx_g_hf + 413);

    auto gs_y_zzzzz_xyz = pbuffer.data(idx_g_hf + 414);

    auto gs_y_zzzzz_xzz = pbuffer.data(idx_g_hf + 415);

    auto gs_y_zzzzz_yyy = pbuffer.data(idx_g_hf + 416);

    auto gs_y_zzzzz_yyz = pbuffer.data(idx_g_hf + 417);

    auto gs_y_zzzzz_yzz = pbuffer.data(idx_g_hf + 418);

    auto gs_y_zzzzz_zzz = pbuffer.data(idx_g_hf + 419);

    #pragma omp simd aligned(gc_y, gs_y_zzzzz_xxx, gs_y_zzzzz_xxy, gs_y_zzzzz_xxz, gs_y_zzzzz_xyy, gs_y_zzzzz_xyz, gs_y_zzzzz_xzz, gs_y_zzzzz_yyy, gs_y_zzzzz_yyz, gs_y_zzzzz_yzz, gs_y_zzzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzzz_xxx[i] = 2.0 * ts_zzzzz_xxx[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxy[i] = 2.0 * ts_zzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxz[i] = 2.0 * ts_zzzzz_xxz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyy[i] = 4.0 * ts_zzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyz[i] = 2.0 * ts_zzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xzz[i] = 2.0 * ts_zzzzz_xzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyy[i] = 6.0 * ts_zzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyz[i] = 4.0 * ts_zzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yzz[i] = 2.0 * ts_zzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_zzz[i] = 2.0 * ts_zzzzz_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 420-430 components of targeted buffer : HF

    auto gs_z_xxxxx_xxx = pbuffer.data(idx_g_hf + 420);

    auto gs_z_xxxxx_xxy = pbuffer.data(idx_g_hf + 421);

    auto gs_z_xxxxx_xxz = pbuffer.data(idx_g_hf + 422);

    auto gs_z_xxxxx_xyy = pbuffer.data(idx_g_hf + 423);

    auto gs_z_xxxxx_xyz = pbuffer.data(idx_g_hf + 424);

    auto gs_z_xxxxx_xzz = pbuffer.data(idx_g_hf + 425);

    auto gs_z_xxxxx_yyy = pbuffer.data(idx_g_hf + 426);

    auto gs_z_xxxxx_yyz = pbuffer.data(idx_g_hf + 427);

    auto gs_z_xxxxx_yzz = pbuffer.data(idx_g_hf + 428);

    auto gs_z_xxxxx_zzz = pbuffer.data(idx_g_hf + 429);

    #pragma omp simd aligned(gc_z, gs_z_xxxxx_xxx, gs_z_xxxxx_xxy, gs_z_xxxxx_xxz, gs_z_xxxxx_xyy, gs_z_xxxxx_xyz, gs_z_xxxxx_xzz, gs_z_xxxxx_yyy, gs_z_xxxxx_yyz, gs_z_xxxxx_yzz, gs_z_xxxxx_zzz, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_yy, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yz, ts_xxxxx_yzz, ts_xxxxx_zz, ts_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxx_xxx[i] = 2.0 * ts_xxxxx_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxy[i] = 2.0 * ts_xxxxx_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxz[i] = 2.0 * ts_xxxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyy[i] = 2.0 * ts_xxxxx_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyz[i] = 2.0 * ts_xxxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xzz[i] = 4.0 * ts_xxxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyy[i] = 2.0 * ts_xxxxx_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyz[i] = 2.0 * ts_xxxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yzz[i] = 4.0 * ts_xxxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_zzz[i] = 6.0 * ts_xxxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 430-440 components of targeted buffer : HF

    auto gs_z_xxxxy_xxx = pbuffer.data(idx_g_hf + 430);

    auto gs_z_xxxxy_xxy = pbuffer.data(idx_g_hf + 431);

    auto gs_z_xxxxy_xxz = pbuffer.data(idx_g_hf + 432);

    auto gs_z_xxxxy_xyy = pbuffer.data(idx_g_hf + 433);

    auto gs_z_xxxxy_xyz = pbuffer.data(idx_g_hf + 434);

    auto gs_z_xxxxy_xzz = pbuffer.data(idx_g_hf + 435);

    auto gs_z_xxxxy_yyy = pbuffer.data(idx_g_hf + 436);

    auto gs_z_xxxxy_yyz = pbuffer.data(idx_g_hf + 437);

    auto gs_z_xxxxy_yzz = pbuffer.data(idx_g_hf + 438);

    auto gs_z_xxxxy_zzz = pbuffer.data(idx_g_hf + 439);

    #pragma omp simd aligned(gc_z, gs_z_xxxxy_xxx, gs_z_xxxxy_xxy, gs_z_xxxxy_xxz, gs_z_xxxxy_xyy, gs_z_xxxxy_xyz, gs_z_xxxxy_xzz, gs_z_xxxxy_yyy, gs_z_xxxxy_yyz, gs_z_xxxxy_yzz, gs_z_xxxxy_zzz, ts_xxxxy_xx, ts_xxxxy_xxx, ts_xxxxy_xxy, ts_xxxxy_xxz, ts_xxxxy_xy, ts_xxxxy_xyy, ts_xxxxy_xyz, ts_xxxxy_xz, ts_xxxxy_xzz, ts_xxxxy_yy, ts_xxxxy_yyy, ts_xxxxy_yyz, ts_xxxxy_yz, ts_xxxxy_yzz, ts_xxxxy_zz, ts_xxxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxy_xxx[i] = 2.0 * ts_xxxxy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxy[i] = 2.0 * ts_xxxxy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxz[i] = 2.0 * ts_xxxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyy[i] = 2.0 * ts_xxxxy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyz[i] = 2.0 * ts_xxxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xzz[i] = 4.0 * ts_xxxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyy[i] = 2.0 * ts_xxxxy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyz[i] = 2.0 * ts_xxxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yzz[i] = 4.0 * ts_xxxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_zzz[i] = 6.0 * ts_xxxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 440-450 components of targeted buffer : HF

    auto gs_z_xxxxz_xxx = pbuffer.data(idx_g_hf + 440);

    auto gs_z_xxxxz_xxy = pbuffer.data(idx_g_hf + 441);

    auto gs_z_xxxxz_xxz = pbuffer.data(idx_g_hf + 442);

    auto gs_z_xxxxz_xyy = pbuffer.data(idx_g_hf + 443);

    auto gs_z_xxxxz_xyz = pbuffer.data(idx_g_hf + 444);

    auto gs_z_xxxxz_xzz = pbuffer.data(idx_g_hf + 445);

    auto gs_z_xxxxz_yyy = pbuffer.data(idx_g_hf + 446);

    auto gs_z_xxxxz_yyz = pbuffer.data(idx_g_hf + 447);

    auto gs_z_xxxxz_yzz = pbuffer.data(idx_g_hf + 448);

    auto gs_z_xxxxz_zzz = pbuffer.data(idx_g_hf + 449);

    #pragma omp simd aligned(gc_z, gs_z_xxxxz_xxx, gs_z_xxxxz_xxy, gs_z_xxxxz_xxz, gs_z_xxxxz_xyy, gs_z_xxxxz_xyz, gs_z_xxxxz_xzz, gs_z_xxxxz_yyy, gs_z_xxxxz_yyz, gs_z_xxxxz_yzz, gs_z_xxxxz_zzz, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xzz, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yzz, ts_xxxx_zzz, ts_xxxxz_xx, ts_xxxxz_xxx, ts_xxxxz_xxy, ts_xxxxz_xxz, ts_xxxxz_xy, ts_xxxxz_xyy, ts_xxxxz_xyz, ts_xxxxz_xz, ts_xxxxz_xzz, ts_xxxxz_yy, ts_xxxxz_yyy, ts_xxxxz_yyz, ts_xxxxz_yz, ts_xxxxz_yzz, ts_xxxxz_zz, ts_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxz_xxx[i] = 2.0 * ts_xxxx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxy[i] = 2.0 * ts_xxxx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxz[i] = 2.0 * ts_xxxx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyy[i] = 2.0 * ts_xxxx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyz[i] = 2.0 * ts_xxxx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xzz[i] = 2.0 * ts_xxxx_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyy[i] = 2.0 * ts_xxxx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyz[i] = 2.0 * ts_xxxx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yzz[i] = 2.0 * ts_xxxx_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_zzz[i] = 2.0 * ts_xxxx_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 450-460 components of targeted buffer : HF

    auto gs_z_xxxyy_xxx = pbuffer.data(idx_g_hf + 450);

    auto gs_z_xxxyy_xxy = pbuffer.data(idx_g_hf + 451);

    auto gs_z_xxxyy_xxz = pbuffer.data(idx_g_hf + 452);

    auto gs_z_xxxyy_xyy = pbuffer.data(idx_g_hf + 453);

    auto gs_z_xxxyy_xyz = pbuffer.data(idx_g_hf + 454);

    auto gs_z_xxxyy_xzz = pbuffer.data(idx_g_hf + 455);

    auto gs_z_xxxyy_yyy = pbuffer.data(idx_g_hf + 456);

    auto gs_z_xxxyy_yyz = pbuffer.data(idx_g_hf + 457);

    auto gs_z_xxxyy_yzz = pbuffer.data(idx_g_hf + 458);

    auto gs_z_xxxyy_zzz = pbuffer.data(idx_g_hf + 459);

    #pragma omp simd aligned(gc_z, gs_z_xxxyy_xxx, gs_z_xxxyy_xxy, gs_z_xxxyy_xxz, gs_z_xxxyy_xyy, gs_z_xxxyy_xyz, gs_z_xxxyy_xzz, gs_z_xxxyy_yyy, gs_z_xxxyy_yyz, gs_z_xxxyy_yzz, gs_z_xxxyy_zzz, ts_xxxyy_xx, ts_xxxyy_xxx, ts_xxxyy_xxy, ts_xxxyy_xxz, ts_xxxyy_xy, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_xz, ts_xxxyy_xzz, ts_xxxyy_yy, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yz, ts_xxxyy_yzz, ts_xxxyy_zz, ts_xxxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyy_xxx[i] = 2.0 * ts_xxxyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxy[i] = 2.0 * ts_xxxyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxz[i] = 2.0 * ts_xxxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyy[i] = 2.0 * ts_xxxyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyz[i] = 2.0 * ts_xxxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xzz[i] = 4.0 * ts_xxxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyy[i] = 2.0 * ts_xxxyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyz[i] = 2.0 * ts_xxxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yzz[i] = 4.0 * ts_xxxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_zzz[i] = 6.0 * ts_xxxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 460-470 components of targeted buffer : HF

    auto gs_z_xxxyz_xxx = pbuffer.data(idx_g_hf + 460);

    auto gs_z_xxxyz_xxy = pbuffer.data(idx_g_hf + 461);

    auto gs_z_xxxyz_xxz = pbuffer.data(idx_g_hf + 462);

    auto gs_z_xxxyz_xyy = pbuffer.data(idx_g_hf + 463);

    auto gs_z_xxxyz_xyz = pbuffer.data(idx_g_hf + 464);

    auto gs_z_xxxyz_xzz = pbuffer.data(idx_g_hf + 465);

    auto gs_z_xxxyz_yyy = pbuffer.data(idx_g_hf + 466);

    auto gs_z_xxxyz_yyz = pbuffer.data(idx_g_hf + 467);

    auto gs_z_xxxyz_yzz = pbuffer.data(idx_g_hf + 468);

    auto gs_z_xxxyz_zzz = pbuffer.data(idx_g_hf + 469);

    #pragma omp simd aligned(gc_z, gs_z_xxxyz_xxx, gs_z_xxxyz_xxy, gs_z_xxxyz_xxz, gs_z_xxxyz_xyy, gs_z_xxxyz_xyz, gs_z_xxxyz_xzz, gs_z_xxxyz_yyy, gs_z_xxxyz_yyz, gs_z_xxxyz_yzz, gs_z_xxxyz_zzz, ts_xxxy_xxx, ts_xxxy_xxy, ts_xxxy_xxz, ts_xxxy_xyy, ts_xxxy_xyz, ts_xxxy_xzz, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yzz, ts_xxxy_zzz, ts_xxxyz_xx, ts_xxxyz_xxx, ts_xxxyz_xxy, ts_xxxyz_xxz, ts_xxxyz_xy, ts_xxxyz_xyy, ts_xxxyz_xyz, ts_xxxyz_xz, ts_xxxyz_xzz, ts_xxxyz_yy, ts_xxxyz_yyy, ts_xxxyz_yyz, ts_xxxyz_yz, ts_xxxyz_yzz, ts_xxxyz_zz, ts_xxxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyz_xxx[i] = 2.0 * ts_xxxy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxy[i] = 2.0 * ts_xxxy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxz[i] = 2.0 * ts_xxxy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyy[i] = 2.0 * ts_xxxy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyz[i] = 2.0 * ts_xxxy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xzz[i] = 2.0 * ts_xxxy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyy[i] = 2.0 * ts_xxxy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyz[i] = 2.0 * ts_xxxy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yzz[i] = 2.0 * ts_xxxy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_zzz[i] = 2.0 * ts_xxxy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 470-480 components of targeted buffer : HF

    auto gs_z_xxxzz_xxx = pbuffer.data(idx_g_hf + 470);

    auto gs_z_xxxzz_xxy = pbuffer.data(idx_g_hf + 471);

    auto gs_z_xxxzz_xxz = pbuffer.data(idx_g_hf + 472);

    auto gs_z_xxxzz_xyy = pbuffer.data(idx_g_hf + 473);

    auto gs_z_xxxzz_xyz = pbuffer.data(idx_g_hf + 474);

    auto gs_z_xxxzz_xzz = pbuffer.data(idx_g_hf + 475);

    auto gs_z_xxxzz_yyy = pbuffer.data(idx_g_hf + 476);

    auto gs_z_xxxzz_yyz = pbuffer.data(idx_g_hf + 477);

    auto gs_z_xxxzz_yzz = pbuffer.data(idx_g_hf + 478);

    auto gs_z_xxxzz_zzz = pbuffer.data(idx_g_hf + 479);

    #pragma omp simd aligned(gc_z, gs_z_xxxzz_xxx, gs_z_xxxzz_xxy, gs_z_xxxzz_xxz, gs_z_xxxzz_xyy, gs_z_xxxzz_xyz, gs_z_xxxzz_xzz, gs_z_xxxzz_yyy, gs_z_xxxzz_yyz, gs_z_xxxzz_yzz, gs_z_xxxzz_zzz, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xxz, ts_xxxz_xyy, ts_xxxz_xyz, ts_xxxz_xzz, ts_xxxz_yyy, ts_xxxz_yyz, ts_xxxz_yzz, ts_xxxz_zzz, ts_xxxzz_xx, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xxz, ts_xxxzz_xy, ts_xxxzz_xyy, ts_xxxzz_xyz, ts_xxxzz_xz, ts_xxxzz_xzz, ts_xxxzz_yy, ts_xxxzz_yyy, ts_xxxzz_yyz, ts_xxxzz_yz, ts_xxxzz_yzz, ts_xxxzz_zz, ts_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxzz_xxx[i] = 4.0 * ts_xxxz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxy[i] = 4.0 * ts_xxxz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxz[i] = 4.0 * ts_xxxz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyy[i] = 4.0 * ts_xxxz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyz[i] = 4.0 * ts_xxxz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xzz[i] = 4.0 * ts_xxxz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyy[i] = 4.0 * ts_xxxz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyz[i] = 4.0 * ts_xxxz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yzz[i] = 4.0 * ts_xxxz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_zzz[i] = 4.0 * ts_xxxz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 480-490 components of targeted buffer : HF

    auto gs_z_xxyyy_xxx = pbuffer.data(idx_g_hf + 480);

    auto gs_z_xxyyy_xxy = pbuffer.data(idx_g_hf + 481);

    auto gs_z_xxyyy_xxz = pbuffer.data(idx_g_hf + 482);

    auto gs_z_xxyyy_xyy = pbuffer.data(idx_g_hf + 483);

    auto gs_z_xxyyy_xyz = pbuffer.data(idx_g_hf + 484);

    auto gs_z_xxyyy_xzz = pbuffer.data(idx_g_hf + 485);

    auto gs_z_xxyyy_yyy = pbuffer.data(idx_g_hf + 486);

    auto gs_z_xxyyy_yyz = pbuffer.data(idx_g_hf + 487);

    auto gs_z_xxyyy_yzz = pbuffer.data(idx_g_hf + 488);

    auto gs_z_xxyyy_zzz = pbuffer.data(idx_g_hf + 489);

    #pragma omp simd aligned(gc_z, gs_z_xxyyy_xxx, gs_z_xxyyy_xxy, gs_z_xxyyy_xxz, gs_z_xxyyy_xyy, gs_z_xxyyy_xyz, gs_z_xxyyy_xzz, gs_z_xxyyy_yyy, gs_z_xxyyy_yyz, gs_z_xxyyy_yzz, gs_z_xxyyy_zzz, ts_xxyyy_xx, ts_xxyyy_xxx, ts_xxyyy_xxy, ts_xxyyy_xxz, ts_xxyyy_xy, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_xz, ts_xxyyy_xzz, ts_xxyyy_yy, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yz, ts_xxyyy_yzz, ts_xxyyy_zz, ts_xxyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyy_xxx[i] = 2.0 * ts_xxyyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxy[i] = 2.0 * ts_xxyyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxz[i] = 2.0 * ts_xxyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyy[i] = 2.0 * ts_xxyyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyz[i] = 2.0 * ts_xxyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xzz[i] = 4.0 * ts_xxyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyy[i] = 2.0 * ts_xxyyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyz[i] = 2.0 * ts_xxyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yzz[i] = 4.0 * ts_xxyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_zzz[i] = 6.0 * ts_xxyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 490-500 components of targeted buffer : HF

    auto gs_z_xxyyz_xxx = pbuffer.data(idx_g_hf + 490);

    auto gs_z_xxyyz_xxy = pbuffer.data(idx_g_hf + 491);

    auto gs_z_xxyyz_xxz = pbuffer.data(idx_g_hf + 492);

    auto gs_z_xxyyz_xyy = pbuffer.data(idx_g_hf + 493);

    auto gs_z_xxyyz_xyz = pbuffer.data(idx_g_hf + 494);

    auto gs_z_xxyyz_xzz = pbuffer.data(idx_g_hf + 495);

    auto gs_z_xxyyz_yyy = pbuffer.data(idx_g_hf + 496);

    auto gs_z_xxyyz_yyz = pbuffer.data(idx_g_hf + 497);

    auto gs_z_xxyyz_yzz = pbuffer.data(idx_g_hf + 498);

    auto gs_z_xxyyz_zzz = pbuffer.data(idx_g_hf + 499);

    #pragma omp simd aligned(gc_z, gs_z_xxyyz_xxx, gs_z_xxyyz_xxy, gs_z_xxyyz_xxz, gs_z_xxyyz_xyy, gs_z_xxyyz_xyz, gs_z_xxyyz_xzz, gs_z_xxyyz_yyy, gs_z_xxyyz_yyz, gs_z_xxyyz_yzz, gs_z_xxyyz_zzz, ts_xxyy_xxx, ts_xxyy_xxy, ts_xxyy_xxz, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_xzz, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yzz, ts_xxyy_zzz, ts_xxyyz_xx, ts_xxyyz_xxx, ts_xxyyz_xxy, ts_xxyyz_xxz, ts_xxyyz_xy, ts_xxyyz_xyy, ts_xxyyz_xyz, ts_xxyyz_xz, ts_xxyyz_xzz, ts_xxyyz_yy, ts_xxyyz_yyy, ts_xxyyz_yyz, ts_xxyyz_yz, ts_xxyyz_yzz, ts_xxyyz_zz, ts_xxyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyz_xxx[i] = 2.0 * ts_xxyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxy[i] = 2.0 * ts_xxyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxz[i] = 2.0 * ts_xxyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyy[i] = 2.0 * ts_xxyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyz[i] = 2.0 * ts_xxyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xzz[i] = 2.0 * ts_xxyy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyy[i] = 2.0 * ts_xxyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyz[i] = 2.0 * ts_xxyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yzz[i] = 2.0 * ts_xxyy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_zzz[i] = 2.0 * ts_xxyy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 500-510 components of targeted buffer : HF

    auto gs_z_xxyzz_xxx = pbuffer.data(idx_g_hf + 500);

    auto gs_z_xxyzz_xxy = pbuffer.data(idx_g_hf + 501);

    auto gs_z_xxyzz_xxz = pbuffer.data(idx_g_hf + 502);

    auto gs_z_xxyzz_xyy = pbuffer.data(idx_g_hf + 503);

    auto gs_z_xxyzz_xyz = pbuffer.data(idx_g_hf + 504);

    auto gs_z_xxyzz_xzz = pbuffer.data(idx_g_hf + 505);

    auto gs_z_xxyzz_yyy = pbuffer.data(idx_g_hf + 506);

    auto gs_z_xxyzz_yyz = pbuffer.data(idx_g_hf + 507);

    auto gs_z_xxyzz_yzz = pbuffer.data(idx_g_hf + 508);

    auto gs_z_xxyzz_zzz = pbuffer.data(idx_g_hf + 509);

    #pragma omp simd aligned(gc_z, gs_z_xxyzz_xxx, gs_z_xxyzz_xxy, gs_z_xxyzz_xxz, gs_z_xxyzz_xyy, gs_z_xxyzz_xyz, gs_z_xxyzz_xzz, gs_z_xxyzz_yyy, gs_z_xxyzz_yyz, gs_z_xxyzz_yzz, gs_z_xxyzz_zzz, ts_xxyz_xxx, ts_xxyz_xxy, ts_xxyz_xxz, ts_xxyz_xyy, ts_xxyz_xyz, ts_xxyz_xzz, ts_xxyz_yyy, ts_xxyz_yyz, ts_xxyz_yzz, ts_xxyz_zzz, ts_xxyzz_xx, ts_xxyzz_xxx, ts_xxyzz_xxy, ts_xxyzz_xxz, ts_xxyzz_xy, ts_xxyzz_xyy, ts_xxyzz_xyz, ts_xxyzz_xz, ts_xxyzz_xzz, ts_xxyzz_yy, ts_xxyzz_yyy, ts_xxyzz_yyz, ts_xxyzz_yz, ts_xxyzz_yzz, ts_xxyzz_zz, ts_xxyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyzz_xxx[i] = 4.0 * ts_xxyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxy[i] = 4.0 * ts_xxyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxz[i] = 4.0 * ts_xxyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyy[i] = 4.0 * ts_xxyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyz[i] = 4.0 * ts_xxyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xzz[i] = 4.0 * ts_xxyz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyy[i] = 4.0 * ts_xxyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyz[i] = 4.0 * ts_xxyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yzz[i] = 4.0 * ts_xxyz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_zzz[i] = 4.0 * ts_xxyz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 510-520 components of targeted buffer : HF

    auto gs_z_xxzzz_xxx = pbuffer.data(idx_g_hf + 510);

    auto gs_z_xxzzz_xxy = pbuffer.data(idx_g_hf + 511);

    auto gs_z_xxzzz_xxz = pbuffer.data(idx_g_hf + 512);

    auto gs_z_xxzzz_xyy = pbuffer.data(idx_g_hf + 513);

    auto gs_z_xxzzz_xyz = pbuffer.data(idx_g_hf + 514);

    auto gs_z_xxzzz_xzz = pbuffer.data(idx_g_hf + 515);

    auto gs_z_xxzzz_yyy = pbuffer.data(idx_g_hf + 516);

    auto gs_z_xxzzz_yyz = pbuffer.data(idx_g_hf + 517);

    auto gs_z_xxzzz_yzz = pbuffer.data(idx_g_hf + 518);

    auto gs_z_xxzzz_zzz = pbuffer.data(idx_g_hf + 519);

    #pragma omp simd aligned(gc_z, gs_z_xxzzz_xxx, gs_z_xxzzz_xxy, gs_z_xxzzz_xxz, gs_z_xxzzz_xyy, gs_z_xxzzz_xyz, gs_z_xxzzz_xzz, gs_z_xxzzz_yyy, gs_z_xxzzz_yyz, gs_z_xxzzz_yzz, gs_z_xxzzz_zzz, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xxz, ts_xxzz_xyy, ts_xxzz_xyz, ts_xxzz_xzz, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yzz, ts_xxzz_zzz, ts_xxzzz_xx, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xxz, ts_xxzzz_xy, ts_xxzzz_xyy, ts_xxzzz_xyz, ts_xxzzz_xz, ts_xxzzz_xzz, ts_xxzzz_yy, ts_xxzzz_yyy, ts_xxzzz_yyz, ts_xxzzz_yz, ts_xxzzz_yzz, ts_xxzzz_zz, ts_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzzz_xxx[i] = 6.0 * ts_xxzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxy[i] = 6.0 * ts_xxzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxz[i] = 6.0 * ts_xxzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyy[i] = 6.0 * ts_xxzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyz[i] = 6.0 * ts_xxzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xzz[i] = 6.0 * ts_xxzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyy[i] = 6.0 * ts_xxzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyz[i] = 6.0 * ts_xxzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yzz[i] = 6.0 * ts_xxzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_zzz[i] = 6.0 * ts_xxzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 520-530 components of targeted buffer : HF

    auto gs_z_xyyyy_xxx = pbuffer.data(idx_g_hf + 520);

    auto gs_z_xyyyy_xxy = pbuffer.data(idx_g_hf + 521);

    auto gs_z_xyyyy_xxz = pbuffer.data(idx_g_hf + 522);

    auto gs_z_xyyyy_xyy = pbuffer.data(idx_g_hf + 523);

    auto gs_z_xyyyy_xyz = pbuffer.data(idx_g_hf + 524);

    auto gs_z_xyyyy_xzz = pbuffer.data(idx_g_hf + 525);

    auto gs_z_xyyyy_yyy = pbuffer.data(idx_g_hf + 526);

    auto gs_z_xyyyy_yyz = pbuffer.data(idx_g_hf + 527);

    auto gs_z_xyyyy_yzz = pbuffer.data(idx_g_hf + 528);

    auto gs_z_xyyyy_zzz = pbuffer.data(idx_g_hf + 529);

    #pragma omp simd aligned(gc_z, gs_z_xyyyy_xxx, gs_z_xyyyy_xxy, gs_z_xyyyy_xxz, gs_z_xyyyy_xyy, gs_z_xyyyy_xyz, gs_z_xyyyy_xzz, gs_z_xyyyy_yyy, gs_z_xyyyy_yyz, gs_z_xyyyy_yzz, gs_z_xyyyy_zzz, ts_xyyyy_xx, ts_xyyyy_xxx, ts_xyyyy_xxy, ts_xyyyy_xxz, ts_xyyyy_xy, ts_xyyyy_xyy, ts_xyyyy_xyz, ts_xyyyy_xz, ts_xyyyy_xzz, ts_xyyyy_yy, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yz, ts_xyyyy_yzz, ts_xyyyy_zz, ts_xyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyy_xxx[i] = 2.0 * ts_xyyyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxy[i] = 2.0 * ts_xyyyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxz[i] = 2.0 * ts_xyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyy[i] = 2.0 * ts_xyyyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyz[i] = 2.0 * ts_xyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xzz[i] = 4.0 * ts_xyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyy[i] = 2.0 * ts_xyyyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyz[i] = 2.0 * ts_xyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yzz[i] = 4.0 * ts_xyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_zzz[i] = 6.0 * ts_xyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 530-540 components of targeted buffer : HF

    auto gs_z_xyyyz_xxx = pbuffer.data(idx_g_hf + 530);

    auto gs_z_xyyyz_xxy = pbuffer.data(idx_g_hf + 531);

    auto gs_z_xyyyz_xxz = pbuffer.data(idx_g_hf + 532);

    auto gs_z_xyyyz_xyy = pbuffer.data(idx_g_hf + 533);

    auto gs_z_xyyyz_xyz = pbuffer.data(idx_g_hf + 534);

    auto gs_z_xyyyz_xzz = pbuffer.data(idx_g_hf + 535);

    auto gs_z_xyyyz_yyy = pbuffer.data(idx_g_hf + 536);

    auto gs_z_xyyyz_yyz = pbuffer.data(idx_g_hf + 537);

    auto gs_z_xyyyz_yzz = pbuffer.data(idx_g_hf + 538);

    auto gs_z_xyyyz_zzz = pbuffer.data(idx_g_hf + 539);

    #pragma omp simd aligned(gc_z, gs_z_xyyyz_xxx, gs_z_xyyyz_xxy, gs_z_xyyyz_xxz, gs_z_xyyyz_xyy, gs_z_xyyyz_xyz, gs_z_xyyyz_xzz, gs_z_xyyyz_yyy, gs_z_xyyyz_yyz, gs_z_xyyyz_yzz, gs_z_xyyyz_zzz, ts_xyyy_xxx, ts_xyyy_xxy, ts_xyyy_xxz, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_xzz, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yzz, ts_xyyy_zzz, ts_xyyyz_xx, ts_xyyyz_xxx, ts_xyyyz_xxy, ts_xyyyz_xxz, ts_xyyyz_xy, ts_xyyyz_xyy, ts_xyyyz_xyz, ts_xyyyz_xz, ts_xyyyz_xzz, ts_xyyyz_yy, ts_xyyyz_yyy, ts_xyyyz_yyz, ts_xyyyz_yz, ts_xyyyz_yzz, ts_xyyyz_zz, ts_xyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyz_xxx[i] = 2.0 * ts_xyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxy[i] = 2.0 * ts_xyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxz[i] = 2.0 * ts_xyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyy[i] = 2.0 * ts_xyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyz[i] = 2.0 * ts_xyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xzz[i] = 2.0 * ts_xyyy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyy[i] = 2.0 * ts_xyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyz[i] = 2.0 * ts_xyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yzz[i] = 2.0 * ts_xyyy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_zzz[i] = 2.0 * ts_xyyy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 540-550 components of targeted buffer : HF

    auto gs_z_xyyzz_xxx = pbuffer.data(idx_g_hf + 540);

    auto gs_z_xyyzz_xxy = pbuffer.data(idx_g_hf + 541);

    auto gs_z_xyyzz_xxz = pbuffer.data(idx_g_hf + 542);

    auto gs_z_xyyzz_xyy = pbuffer.data(idx_g_hf + 543);

    auto gs_z_xyyzz_xyz = pbuffer.data(idx_g_hf + 544);

    auto gs_z_xyyzz_xzz = pbuffer.data(idx_g_hf + 545);

    auto gs_z_xyyzz_yyy = pbuffer.data(idx_g_hf + 546);

    auto gs_z_xyyzz_yyz = pbuffer.data(idx_g_hf + 547);

    auto gs_z_xyyzz_yzz = pbuffer.data(idx_g_hf + 548);

    auto gs_z_xyyzz_zzz = pbuffer.data(idx_g_hf + 549);

    #pragma omp simd aligned(gc_z, gs_z_xyyzz_xxx, gs_z_xyyzz_xxy, gs_z_xyyzz_xxz, gs_z_xyyzz_xyy, gs_z_xyyzz_xyz, gs_z_xyyzz_xzz, gs_z_xyyzz_yyy, gs_z_xyyzz_yyz, gs_z_xyyzz_yzz, gs_z_xyyzz_zzz, ts_xyyz_xxx, ts_xyyz_xxy, ts_xyyz_xxz, ts_xyyz_xyy, ts_xyyz_xyz, ts_xyyz_xzz, ts_xyyz_yyy, ts_xyyz_yyz, ts_xyyz_yzz, ts_xyyz_zzz, ts_xyyzz_xx, ts_xyyzz_xxx, ts_xyyzz_xxy, ts_xyyzz_xxz, ts_xyyzz_xy, ts_xyyzz_xyy, ts_xyyzz_xyz, ts_xyyzz_xz, ts_xyyzz_xzz, ts_xyyzz_yy, ts_xyyzz_yyy, ts_xyyzz_yyz, ts_xyyzz_yz, ts_xyyzz_yzz, ts_xyyzz_zz, ts_xyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyzz_xxx[i] = 4.0 * ts_xyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxy[i] = 4.0 * ts_xyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxz[i] = 4.0 * ts_xyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyy[i] = 4.0 * ts_xyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyz[i] = 4.0 * ts_xyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xzz[i] = 4.0 * ts_xyyz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyy[i] = 4.0 * ts_xyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyz[i] = 4.0 * ts_xyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yzz[i] = 4.0 * ts_xyyz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_zzz[i] = 4.0 * ts_xyyz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 550-560 components of targeted buffer : HF

    auto gs_z_xyzzz_xxx = pbuffer.data(idx_g_hf + 550);

    auto gs_z_xyzzz_xxy = pbuffer.data(idx_g_hf + 551);

    auto gs_z_xyzzz_xxz = pbuffer.data(idx_g_hf + 552);

    auto gs_z_xyzzz_xyy = pbuffer.data(idx_g_hf + 553);

    auto gs_z_xyzzz_xyz = pbuffer.data(idx_g_hf + 554);

    auto gs_z_xyzzz_xzz = pbuffer.data(idx_g_hf + 555);

    auto gs_z_xyzzz_yyy = pbuffer.data(idx_g_hf + 556);

    auto gs_z_xyzzz_yyz = pbuffer.data(idx_g_hf + 557);

    auto gs_z_xyzzz_yzz = pbuffer.data(idx_g_hf + 558);

    auto gs_z_xyzzz_zzz = pbuffer.data(idx_g_hf + 559);

    #pragma omp simd aligned(gc_z, gs_z_xyzzz_xxx, gs_z_xyzzz_xxy, gs_z_xyzzz_xxz, gs_z_xyzzz_xyy, gs_z_xyzzz_xyz, gs_z_xyzzz_xzz, gs_z_xyzzz_yyy, gs_z_xyzzz_yyz, gs_z_xyzzz_yzz, gs_z_xyzzz_zzz, ts_xyzz_xxx, ts_xyzz_xxy, ts_xyzz_xxz, ts_xyzz_xyy, ts_xyzz_xyz, ts_xyzz_xzz, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yzz, ts_xyzz_zzz, ts_xyzzz_xx, ts_xyzzz_xxx, ts_xyzzz_xxy, ts_xyzzz_xxz, ts_xyzzz_xy, ts_xyzzz_xyy, ts_xyzzz_xyz, ts_xyzzz_xz, ts_xyzzz_xzz, ts_xyzzz_yy, ts_xyzzz_yyy, ts_xyzzz_yyz, ts_xyzzz_yz, ts_xyzzz_yzz, ts_xyzzz_zz, ts_xyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzzz_xxx[i] = 6.0 * ts_xyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxy[i] = 6.0 * ts_xyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxz[i] = 6.0 * ts_xyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyy[i] = 6.0 * ts_xyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyz[i] = 6.0 * ts_xyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xzz[i] = 6.0 * ts_xyzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyy[i] = 6.0 * ts_xyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyz[i] = 6.0 * ts_xyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yzz[i] = 6.0 * ts_xyzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_zzz[i] = 6.0 * ts_xyzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 560-570 components of targeted buffer : HF

    auto gs_z_xzzzz_xxx = pbuffer.data(idx_g_hf + 560);

    auto gs_z_xzzzz_xxy = pbuffer.data(idx_g_hf + 561);

    auto gs_z_xzzzz_xxz = pbuffer.data(idx_g_hf + 562);

    auto gs_z_xzzzz_xyy = pbuffer.data(idx_g_hf + 563);

    auto gs_z_xzzzz_xyz = pbuffer.data(idx_g_hf + 564);

    auto gs_z_xzzzz_xzz = pbuffer.data(idx_g_hf + 565);

    auto gs_z_xzzzz_yyy = pbuffer.data(idx_g_hf + 566);

    auto gs_z_xzzzz_yyz = pbuffer.data(idx_g_hf + 567);

    auto gs_z_xzzzz_yzz = pbuffer.data(idx_g_hf + 568);

    auto gs_z_xzzzz_zzz = pbuffer.data(idx_g_hf + 569);

    #pragma omp simd aligned(gc_z, gs_z_xzzzz_xxx, gs_z_xzzzz_xxy, gs_z_xzzzz_xxz, gs_z_xzzzz_xyy, gs_z_xzzzz_xyz, gs_z_xzzzz_xzz, gs_z_xzzzz_yyy, gs_z_xzzzz_yyz, gs_z_xzzzz_yzz, gs_z_xzzzz_zzz, ts_xzzz_xxx, ts_xzzz_xxy, ts_xzzz_xxz, ts_xzzz_xyy, ts_xzzz_xyz, ts_xzzz_xzz, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yzz, ts_xzzz_zzz, ts_xzzzz_xx, ts_xzzzz_xxx, ts_xzzzz_xxy, ts_xzzzz_xxz, ts_xzzzz_xy, ts_xzzzz_xyy, ts_xzzzz_xyz, ts_xzzzz_xz, ts_xzzzz_xzz, ts_xzzzz_yy, ts_xzzzz_yyy, ts_xzzzz_yyz, ts_xzzzz_yz, ts_xzzzz_yzz, ts_xzzzz_zz, ts_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzzz_xxx[i] = 8.0 * ts_xzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxy[i] = 8.0 * ts_xzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxz[i] = 8.0 * ts_xzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyy[i] = 8.0 * ts_xzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyz[i] = 8.0 * ts_xzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xzz[i] = 8.0 * ts_xzzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyy[i] = 8.0 * ts_xzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyz[i] = 8.0 * ts_xzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yzz[i] = 8.0 * ts_xzzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_zzz[i] = 8.0 * ts_xzzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 570-580 components of targeted buffer : HF

    auto gs_z_yyyyy_xxx = pbuffer.data(idx_g_hf + 570);

    auto gs_z_yyyyy_xxy = pbuffer.data(idx_g_hf + 571);

    auto gs_z_yyyyy_xxz = pbuffer.data(idx_g_hf + 572);

    auto gs_z_yyyyy_xyy = pbuffer.data(idx_g_hf + 573);

    auto gs_z_yyyyy_xyz = pbuffer.data(idx_g_hf + 574);

    auto gs_z_yyyyy_xzz = pbuffer.data(idx_g_hf + 575);

    auto gs_z_yyyyy_yyy = pbuffer.data(idx_g_hf + 576);

    auto gs_z_yyyyy_yyz = pbuffer.data(idx_g_hf + 577);

    auto gs_z_yyyyy_yzz = pbuffer.data(idx_g_hf + 578);

    auto gs_z_yyyyy_zzz = pbuffer.data(idx_g_hf + 579);

    #pragma omp simd aligned(gc_z, gs_z_yyyyy_xxx, gs_z_yyyyy_xxy, gs_z_yyyyy_xxz, gs_z_yyyyy_xyy, gs_z_yyyyy_xyz, gs_z_yyyyy_xzz, gs_z_yyyyy_yyy, gs_z_yyyyy_yyz, gs_z_yyyyy_yzz, gs_z_yyyyy_zzz, ts_yyyyy_xx, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xz, ts_yyyyy_xzz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyy_zz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyy_xxx[i] = 2.0 * ts_yyyyy_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxy[i] = 2.0 * ts_yyyyy_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxz[i] = 2.0 * ts_yyyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyy[i] = 2.0 * ts_yyyyy_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyz[i] = 2.0 * ts_yyyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xzz[i] = 4.0 * ts_yyyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyy[i] = 2.0 * ts_yyyyy_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyz[i] = 2.0 * ts_yyyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yzz[i] = 4.0 * ts_yyyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_zzz[i] = 6.0 * ts_yyyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 580-590 components of targeted buffer : HF

    auto gs_z_yyyyz_xxx = pbuffer.data(idx_g_hf + 580);

    auto gs_z_yyyyz_xxy = pbuffer.data(idx_g_hf + 581);

    auto gs_z_yyyyz_xxz = pbuffer.data(idx_g_hf + 582);

    auto gs_z_yyyyz_xyy = pbuffer.data(idx_g_hf + 583);

    auto gs_z_yyyyz_xyz = pbuffer.data(idx_g_hf + 584);

    auto gs_z_yyyyz_xzz = pbuffer.data(idx_g_hf + 585);

    auto gs_z_yyyyz_yyy = pbuffer.data(idx_g_hf + 586);

    auto gs_z_yyyyz_yyz = pbuffer.data(idx_g_hf + 587);

    auto gs_z_yyyyz_yzz = pbuffer.data(idx_g_hf + 588);

    auto gs_z_yyyyz_zzz = pbuffer.data(idx_g_hf + 589);

    #pragma omp simd aligned(gc_z, gs_z_yyyyz_xxx, gs_z_yyyyz_xxy, gs_z_yyyyz_xxz, gs_z_yyyyz_xyy, gs_z_yyyyz_xyz, gs_z_yyyyz_xzz, gs_z_yyyyz_yyy, gs_z_yyyyz_yyz, gs_z_yyyyz_yzz, gs_z_yyyyz_zzz, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xzz, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yzz, ts_yyyy_zzz, ts_yyyyz_xx, ts_yyyyz_xxx, ts_yyyyz_xxy, ts_yyyyz_xxz, ts_yyyyz_xy, ts_yyyyz_xyy, ts_yyyyz_xyz, ts_yyyyz_xz, ts_yyyyz_xzz, ts_yyyyz_yy, ts_yyyyz_yyy, ts_yyyyz_yyz, ts_yyyyz_yz, ts_yyyyz_yzz, ts_yyyyz_zz, ts_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyz_xxx[i] = 2.0 * ts_yyyy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxy[i] = 2.0 * ts_yyyy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxz[i] = 2.0 * ts_yyyy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyy[i] = 2.0 * ts_yyyy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyz[i] = 2.0 * ts_yyyy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xzz[i] = 2.0 * ts_yyyy_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyy[i] = 2.0 * ts_yyyy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyz[i] = 2.0 * ts_yyyy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yzz[i] = 2.0 * ts_yyyy_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_zzz[i] = 2.0 * ts_yyyy_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 590-600 components of targeted buffer : HF

    auto gs_z_yyyzz_xxx = pbuffer.data(idx_g_hf + 590);

    auto gs_z_yyyzz_xxy = pbuffer.data(idx_g_hf + 591);

    auto gs_z_yyyzz_xxz = pbuffer.data(idx_g_hf + 592);

    auto gs_z_yyyzz_xyy = pbuffer.data(idx_g_hf + 593);

    auto gs_z_yyyzz_xyz = pbuffer.data(idx_g_hf + 594);

    auto gs_z_yyyzz_xzz = pbuffer.data(idx_g_hf + 595);

    auto gs_z_yyyzz_yyy = pbuffer.data(idx_g_hf + 596);

    auto gs_z_yyyzz_yyz = pbuffer.data(idx_g_hf + 597);

    auto gs_z_yyyzz_yzz = pbuffer.data(idx_g_hf + 598);

    auto gs_z_yyyzz_zzz = pbuffer.data(idx_g_hf + 599);

    #pragma omp simd aligned(gc_z, gs_z_yyyzz_xxx, gs_z_yyyzz_xxy, gs_z_yyyzz_xxz, gs_z_yyyzz_xyy, gs_z_yyyzz_xyz, gs_z_yyyzz_xzz, gs_z_yyyzz_yyy, gs_z_yyyzz_yyz, gs_z_yyyzz_yzz, gs_z_yyyzz_zzz, ts_yyyz_xxx, ts_yyyz_xxy, ts_yyyz_xxz, ts_yyyz_xyy, ts_yyyz_xyz, ts_yyyz_xzz, ts_yyyz_yyy, ts_yyyz_yyz, ts_yyyz_yzz, ts_yyyz_zzz, ts_yyyzz_xx, ts_yyyzz_xxx, ts_yyyzz_xxy, ts_yyyzz_xxz, ts_yyyzz_xy, ts_yyyzz_xyy, ts_yyyzz_xyz, ts_yyyzz_xz, ts_yyyzz_xzz, ts_yyyzz_yy, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yz, ts_yyyzz_yzz, ts_yyyzz_zz, ts_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyzz_xxx[i] = 4.0 * ts_yyyz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxy[i] = 4.0 * ts_yyyz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxz[i] = 4.0 * ts_yyyz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyy[i] = 4.0 * ts_yyyz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyz[i] = 4.0 * ts_yyyz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xzz[i] = 4.0 * ts_yyyz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyy[i] = 4.0 * ts_yyyz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyz[i] = 4.0 * ts_yyyz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yzz[i] = 4.0 * ts_yyyz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_zzz[i] = 4.0 * ts_yyyz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 600-610 components of targeted buffer : HF

    auto gs_z_yyzzz_xxx = pbuffer.data(idx_g_hf + 600);

    auto gs_z_yyzzz_xxy = pbuffer.data(idx_g_hf + 601);

    auto gs_z_yyzzz_xxz = pbuffer.data(idx_g_hf + 602);

    auto gs_z_yyzzz_xyy = pbuffer.data(idx_g_hf + 603);

    auto gs_z_yyzzz_xyz = pbuffer.data(idx_g_hf + 604);

    auto gs_z_yyzzz_xzz = pbuffer.data(idx_g_hf + 605);

    auto gs_z_yyzzz_yyy = pbuffer.data(idx_g_hf + 606);

    auto gs_z_yyzzz_yyz = pbuffer.data(idx_g_hf + 607);

    auto gs_z_yyzzz_yzz = pbuffer.data(idx_g_hf + 608);

    auto gs_z_yyzzz_zzz = pbuffer.data(idx_g_hf + 609);

    #pragma omp simd aligned(gc_z, gs_z_yyzzz_xxx, gs_z_yyzzz_xxy, gs_z_yyzzz_xxz, gs_z_yyzzz_xyy, gs_z_yyzzz_xyz, gs_z_yyzzz_xzz, gs_z_yyzzz_yyy, gs_z_yyzzz_yyz, gs_z_yyzzz_yzz, gs_z_yyzzz_zzz, ts_yyzz_xxx, ts_yyzz_xxy, ts_yyzz_xxz, ts_yyzz_xyy, ts_yyzz_xyz, ts_yyzz_xzz, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yzz, ts_yyzz_zzz, ts_yyzzz_xx, ts_yyzzz_xxx, ts_yyzzz_xxy, ts_yyzzz_xxz, ts_yyzzz_xy, ts_yyzzz_xyy, ts_yyzzz_xyz, ts_yyzzz_xz, ts_yyzzz_xzz, ts_yyzzz_yy, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yz, ts_yyzzz_yzz, ts_yyzzz_zz, ts_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzzz_xxx[i] = 6.0 * ts_yyzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxy[i] = 6.0 * ts_yyzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxz[i] = 6.0 * ts_yyzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyy[i] = 6.0 * ts_yyzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyz[i] = 6.0 * ts_yyzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xzz[i] = 6.0 * ts_yyzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyy[i] = 6.0 * ts_yyzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyz[i] = 6.0 * ts_yyzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yzz[i] = 6.0 * ts_yyzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_zzz[i] = 6.0 * ts_yyzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 610-620 components of targeted buffer : HF

    auto gs_z_yzzzz_xxx = pbuffer.data(idx_g_hf + 610);

    auto gs_z_yzzzz_xxy = pbuffer.data(idx_g_hf + 611);

    auto gs_z_yzzzz_xxz = pbuffer.data(idx_g_hf + 612);

    auto gs_z_yzzzz_xyy = pbuffer.data(idx_g_hf + 613);

    auto gs_z_yzzzz_xyz = pbuffer.data(idx_g_hf + 614);

    auto gs_z_yzzzz_xzz = pbuffer.data(idx_g_hf + 615);

    auto gs_z_yzzzz_yyy = pbuffer.data(idx_g_hf + 616);

    auto gs_z_yzzzz_yyz = pbuffer.data(idx_g_hf + 617);

    auto gs_z_yzzzz_yzz = pbuffer.data(idx_g_hf + 618);

    auto gs_z_yzzzz_zzz = pbuffer.data(idx_g_hf + 619);

    #pragma omp simd aligned(gc_z, gs_z_yzzzz_xxx, gs_z_yzzzz_xxy, gs_z_yzzzz_xxz, gs_z_yzzzz_xyy, gs_z_yzzzz_xyz, gs_z_yzzzz_xzz, gs_z_yzzzz_yyy, gs_z_yzzzz_yyz, gs_z_yzzzz_yzz, gs_z_yzzzz_zzz, ts_yzzz_xxx, ts_yzzz_xxy, ts_yzzz_xxz, ts_yzzz_xyy, ts_yzzz_xyz, ts_yzzz_xzz, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yzz, ts_yzzz_zzz, ts_yzzzz_xx, ts_yzzzz_xxx, ts_yzzzz_xxy, ts_yzzzz_xxz, ts_yzzzz_xy, ts_yzzzz_xyy, ts_yzzzz_xyz, ts_yzzzz_xz, ts_yzzzz_xzz, ts_yzzzz_yy, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yz, ts_yzzzz_yzz, ts_yzzzz_zz, ts_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzzz_xxx[i] = 8.0 * ts_yzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxy[i] = 8.0 * ts_yzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxz[i] = 8.0 * ts_yzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyy[i] = 8.0 * ts_yzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyz[i] = 8.0 * ts_yzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xzz[i] = 8.0 * ts_yzzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyy[i] = 8.0 * ts_yzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyz[i] = 8.0 * ts_yzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yzz[i] = 8.0 * ts_yzzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_zzz[i] = 8.0 * ts_yzzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 620-630 components of targeted buffer : HF

    auto gs_z_zzzzz_xxx = pbuffer.data(idx_g_hf + 620);

    auto gs_z_zzzzz_xxy = pbuffer.data(idx_g_hf + 621);

    auto gs_z_zzzzz_xxz = pbuffer.data(idx_g_hf + 622);

    auto gs_z_zzzzz_xyy = pbuffer.data(idx_g_hf + 623);

    auto gs_z_zzzzz_xyz = pbuffer.data(idx_g_hf + 624);

    auto gs_z_zzzzz_xzz = pbuffer.data(idx_g_hf + 625);

    auto gs_z_zzzzz_yyy = pbuffer.data(idx_g_hf + 626);

    auto gs_z_zzzzz_yyz = pbuffer.data(idx_g_hf + 627);

    auto gs_z_zzzzz_yzz = pbuffer.data(idx_g_hf + 628);

    auto gs_z_zzzzz_zzz = pbuffer.data(idx_g_hf + 629);

    #pragma omp simd aligned(gc_z, gs_z_zzzzz_xxx, gs_z_zzzzz_xxy, gs_z_zzzzz_xxz, gs_z_zzzzz_xyy, gs_z_zzzzz_xyz, gs_z_zzzzz_xzz, gs_z_zzzzz_yyy, gs_z_zzzzz_yyz, gs_z_zzzzz_yzz, gs_z_zzzzz_zzz, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzzz_xxx[i] = 10.0 * ts_zzzz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxx[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxy[i] = 10.0 * ts_zzzz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxz[i] = 10.0 * ts_zzzz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyy[i] = 10.0 * ts_zzzz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyz[i] = 10.0 * ts_zzzz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xzz[i] = 10.0 * ts_zzzz_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyy[i] = 10.0 * ts_zzzz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyz[i] = 10.0 * ts_zzzz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yzz[i] = 10.0 * ts_zzzz_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_zzz[i] = 10.0 * ts_zzzz_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_zzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

