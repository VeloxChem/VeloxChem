#include "KineticEnergyPrimRecIF.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_if(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_if,
                            const size_t              idx_ovl_gf,
                            const size_t              idx_kin_gf,
                            const size_t              idx_kin_hd,
                            const size_t              idx_kin_hf,
                            const size_t              idx_ovl_if,
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

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_ovl_gf + 10);

    auto ts_xxxy_xxz = pbuffer.data(idx_ovl_gf + 12);

    auto ts_xxxy_xzz = pbuffer.data(idx_ovl_gf + 15);

    auto ts_xxxz_xxx = pbuffer.data(idx_ovl_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_ovl_gf + 21);

    auto ts_xxxz_xyy = pbuffer.data(idx_ovl_gf + 23);

    auto ts_xxyy_xxx = pbuffer.data(idx_ovl_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_ovl_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_ovl_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_ovl_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_ovl_gf + 39);

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_ovl_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_ovl_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_ovl_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_ovl_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

    auto ts_xyyy_xxy = pbuffer.data(idx_ovl_gf + 61);

    auto ts_xyyy_xyy = pbuffer.data(idx_ovl_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_ovl_gf + 64);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_ovl_gf + 69);

    auto ts_xzzz_xxz = pbuffer.data(idx_ovl_gf + 92);

    auto ts_xzzz_xyz = pbuffer.data(idx_ovl_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_ovl_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_ovl_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

    auto ts_yyyz_xxy = pbuffer.data(idx_ovl_gf + 111);

    auto ts_yyyz_xyy = pbuffer.data(idx_ovl_gf + 113);

    auto ts_yyyz_yyy = pbuffer.data(idx_ovl_gf + 116);

    auto ts_yyzz_xxx = pbuffer.data(idx_ovl_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_ovl_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_ovl_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_ovl_gf + 130);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xyz = pbuffer.data(idx_ovl_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

    // Set up components of auxiliary buffer : GF

    auto tk_xxxx_xxx = pbuffer.data(idx_kin_gf);

    auto tk_xxxx_xxy = pbuffer.data(idx_kin_gf + 1);

    auto tk_xxxx_xxz = pbuffer.data(idx_kin_gf + 2);

    auto tk_xxxx_xyy = pbuffer.data(idx_kin_gf + 3);

    auto tk_xxxx_xyz = pbuffer.data(idx_kin_gf + 4);

    auto tk_xxxx_xzz = pbuffer.data(idx_kin_gf + 5);

    auto tk_xxxx_yyy = pbuffer.data(idx_kin_gf + 6);

    auto tk_xxxx_yyz = pbuffer.data(idx_kin_gf + 7);

    auto tk_xxxx_yzz = pbuffer.data(idx_kin_gf + 8);

    auto tk_xxxx_zzz = pbuffer.data(idx_kin_gf + 9);

    auto tk_xxxy_xxx = pbuffer.data(idx_kin_gf + 10);

    auto tk_xxxy_xxz = pbuffer.data(idx_kin_gf + 12);

    auto tk_xxxy_xzz = pbuffer.data(idx_kin_gf + 15);

    auto tk_xxxz_xxx = pbuffer.data(idx_kin_gf + 20);

    auto tk_xxxz_xxy = pbuffer.data(idx_kin_gf + 21);

    auto tk_xxxz_xyy = pbuffer.data(idx_kin_gf + 23);

    auto tk_xxyy_xxx = pbuffer.data(idx_kin_gf + 30);

    auto tk_xxyy_xxy = pbuffer.data(idx_kin_gf + 31);

    auto tk_xxyy_xxz = pbuffer.data(idx_kin_gf + 32);

    auto tk_xxyy_xyy = pbuffer.data(idx_kin_gf + 33);

    auto tk_xxyy_xyz = pbuffer.data(idx_kin_gf + 34);

    auto tk_xxyy_xzz = pbuffer.data(idx_kin_gf + 35);

    auto tk_xxyy_yyy = pbuffer.data(idx_kin_gf + 36);

    auto tk_xxyy_yyz = pbuffer.data(idx_kin_gf + 37);

    auto tk_xxyy_yzz = pbuffer.data(idx_kin_gf + 38);

    auto tk_xxyy_zzz = pbuffer.data(idx_kin_gf + 39);

    auto tk_xxzz_xxx = pbuffer.data(idx_kin_gf + 50);

    auto tk_xxzz_xxy = pbuffer.data(idx_kin_gf + 51);

    auto tk_xxzz_xxz = pbuffer.data(idx_kin_gf + 52);

    auto tk_xxzz_xyy = pbuffer.data(idx_kin_gf + 53);

    auto tk_xxzz_xyz = pbuffer.data(idx_kin_gf + 54);

    auto tk_xxzz_xzz = pbuffer.data(idx_kin_gf + 55);

    auto tk_xxzz_yyy = pbuffer.data(idx_kin_gf + 56);

    auto tk_xxzz_yyz = pbuffer.data(idx_kin_gf + 57);

    auto tk_xxzz_yzz = pbuffer.data(idx_kin_gf + 58);

    auto tk_xxzz_zzz = pbuffer.data(idx_kin_gf + 59);

    auto tk_xyyy_xxy = pbuffer.data(idx_kin_gf + 61);

    auto tk_xyyy_xyy = pbuffer.data(idx_kin_gf + 63);

    auto tk_xyyy_xyz = pbuffer.data(idx_kin_gf + 64);

    auto tk_xyyy_yyy = pbuffer.data(idx_kin_gf + 66);

    auto tk_xyyy_yyz = pbuffer.data(idx_kin_gf + 67);

    auto tk_xyyy_yzz = pbuffer.data(idx_kin_gf + 68);

    auto tk_xyyy_zzz = pbuffer.data(idx_kin_gf + 69);

    auto tk_xzzz_xxz = pbuffer.data(idx_kin_gf + 92);

    auto tk_xzzz_xyz = pbuffer.data(idx_kin_gf + 94);

    auto tk_xzzz_xzz = pbuffer.data(idx_kin_gf + 95);

    auto tk_xzzz_yyy = pbuffer.data(idx_kin_gf + 96);

    auto tk_xzzz_yyz = pbuffer.data(idx_kin_gf + 97);

    auto tk_xzzz_yzz = pbuffer.data(idx_kin_gf + 98);

    auto tk_xzzz_zzz = pbuffer.data(idx_kin_gf + 99);

    auto tk_yyyy_xxx = pbuffer.data(idx_kin_gf + 100);

    auto tk_yyyy_xxy = pbuffer.data(idx_kin_gf + 101);

    auto tk_yyyy_xxz = pbuffer.data(idx_kin_gf + 102);

    auto tk_yyyy_xyy = pbuffer.data(idx_kin_gf + 103);

    auto tk_yyyy_xyz = pbuffer.data(idx_kin_gf + 104);

    auto tk_yyyy_xzz = pbuffer.data(idx_kin_gf + 105);

    auto tk_yyyy_yyy = pbuffer.data(idx_kin_gf + 106);

    auto tk_yyyy_yyz = pbuffer.data(idx_kin_gf + 107);

    auto tk_yyyy_yzz = pbuffer.data(idx_kin_gf + 108);

    auto tk_yyyy_zzz = pbuffer.data(idx_kin_gf + 109);

    auto tk_yyyz_xxy = pbuffer.data(idx_kin_gf + 111);

    auto tk_yyyz_xyy = pbuffer.data(idx_kin_gf + 113);

    auto tk_yyyz_yyy = pbuffer.data(idx_kin_gf + 116);

    auto tk_yyzz_xxx = pbuffer.data(idx_kin_gf + 120);

    auto tk_yyzz_xxy = pbuffer.data(idx_kin_gf + 121);

    auto tk_yyzz_xxz = pbuffer.data(idx_kin_gf + 122);

    auto tk_yyzz_xyy = pbuffer.data(idx_kin_gf + 123);

    auto tk_yyzz_xyz = pbuffer.data(idx_kin_gf + 124);

    auto tk_yyzz_xzz = pbuffer.data(idx_kin_gf + 125);

    auto tk_yyzz_yyy = pbuffer.data(idx_kin_gf + 126);

    auto tk_yyzz_yyz = pbuffer.data(idx_kin_gf + 127);

    auto tk_yyzz_yzz = pbuffer.data(idx_kin_gf + 128);

    auto tk_yyzz_zzz = pbuffer.data(idx_kin_gf + 129);

    auto tk_yzzz_xxx = pbuffer.data(idx_kin_gf + 130);

    auto tk_yzzz_xxz = pbuffer.data(idx_kin_gf + 132);

    auto tk_yzzz_xyz = pbuffer.data(idx_kin_gf + 134);

    auto tk_yzzz_xzz = pbuffer.data(idx_kin_gf + 135);

    auto tk_yzzz_yyz = pbuffer.data(idx_kin_gf + 137);

    auto tk_yzzz_yzz = pbuffer.data(idx_kin_gf + 138);

    auto tk_yzzz_zzz = pbuffer.data(idx_kin_gf + 139);

    auto tk_zzzz_xxx = pbuffer.data(idx_kin_gf + 140);

    auto tk_zzzz_xxy = pbuffer.data(idx_kin_gf + 141);

    auto tk_zzzz_xxz = pbuffer.data(idx_kin_gf + 142);

    auto tk_zzzz_xyy = pbuffer.data(idx_kin_gf + 143);

    auto tk_zzzz_xyz = pbuffer.data(idx_kin_gf + 144);

    auto tk_zzzz_xzz = pbuffer.data(idx_kin_gf + 145);

    auto tk_zzzz_yyy = pbuffer.data(idx_kin_gf + 146);

    auto tk_zzzz_yyz = pbuffer.data(idx_kin_gf + 147);

    auto tk_zzzz_yzz = pbuffer.data(idx_kin_gf + 148);

    auto tk_zzzz_zzz = pbuffer.data(idx_kin_gf + 149);

    // Set up components of auxiliary buffer : HD

    auto tk_xxxxx_xx = pbuffer.data(idx_kin_hd);

    auto tk_xxxxx_xy = pbuffer.data(idx_kin_hd + 1);

    auto tk_xxxxx_xz = pbuffer.data(idx_kin_hd + 2);

    auto tk_xxxxx_yy = pbuffer.data(idx_kin_hd + 3);

    auto tk_xxxxx_yz = pbuffer.data(idx_kin_hd + 4);

    auto tk_xxxxx_zz = pbuffer.data(idx_kin_hd + 5);

    auto tk_xxxxz_xz = pbuffer.data(idx_kin_hd + 14);

    auto tk_xxxxz_yz = pbuffer.data(idx_kin_hd + 16);

    auto tk_xxxxz_zz = pbuffer.data(idx_kin_hd + 17);

    auto tk_xxxyy_xx = pbuffer.data(idx_kin_hd + 18);

    auto tk_xxxyy_xy = pbuffer.data(idx_kin_hd + 19);

    auto tk_xxxyy_xz = pbuffer.data(idx_kin_hd + 20);

    auto tk_xxxyy_yy = pbuffer.data(idx_kin_hd + 21);

    auto tk_xxxyy_yz = pbuffer.data(idx_kin_hd + 22);

    auto tk_xxxyy_zz = pbuffer.data(idx_kin_hd + 23);

    auto tk_xxxzz_xx = pbuffer.data(idx_kin_hd + 30);

    auto tk_xxxzz_xy = pbuffer.data(idx_kin_hd + 31);

    auto tk_xxxzz_xz = pbuffer.data(idx_kin_hd + 32);

    auto tk_xxxzz_yy = pbuffer.data(idx_kin_hd + 33);

    auto tk_xxxzz_yz = pbuffer.data(idx_kin_hd + 34);

    auto tk_xxxzz_zz = pbuffer.data(idx_kin_hd + 35);

    auto tk_xxyyy_xx = pbuffer.data(idx_kin_hd + 36);

    auto tk_xxyyy_xy = pbuffer.data(idx_kin_hd + 37);

    auto tk_xxyyy_xz = pbuffer.data(idx_kin_hd + 38);

    auto tk_xxyyy_yy = pbuffer.data(idx_kin_hd + 39);

    auto tk_xxyyy_yz = pbuffer.data(idx_kin_hd + 40);

    auto tk_xxyyy_zz = pbuffer.data(idx_kin_hd + 41);

    auto tk_xxzzz_xx = pbuffer.data(idx_kin_hd + 54);

    auto tk_xxzzz_xy = pbuffer.data(idx_kin_hd + 55);

    auto tk_xxzzz_xz = pbuffer.data(idx_kin_hd + 56);

    auto tk_xxzzz_yy = pbuffer.data(idx_kin_hd + 57);

    auto tk_xxzzz_yz = pbuffer.data(idx_kin_hd + 58);

    auto tk_xxzzz_zz = pbuffer.data(idx_kin_hd + 59);

    auto tk_xyyyy_xy = pbuffer.data(idx_kin_hd + 61);

    auto tk_xyyyy_yy = pbuffer.data(idx_kin_hd + 63);

    auto tk_xyyyy_yz = pbuffer.data(idx_kin_hd + 64);

    auto tk_xyyzz_yz = pbuffer.data(idx_kin_hd + 76);

    auto tk_xzzzz_xz = pbuffer.data(idx_kin_hd + 86);

    auto tk_xzzzz_yz = pbuffer.data(idx_kin_hd + 88);

    auto tk_xzzzz_zz = pbuffer.data(idx_kin_hd + 89);

    auto tk_yyyyy_xx = pbuffer.data(idx_kin_hd + 90);

    auto tk_yyyyy_xy = pbuffer.data(idx_kin_hd + 91);

    auto tk_yyyyy_xz = pbuffer.data(idx_kin_hd + 92);

    auto tk_yyyyy_yy = pbuffer.data(idx_kin_hd + 93);

    auto tk_yyyyy_yz = pbuffer.data(idx_kin_hd + 94);

    auto tk_yyyyy_zz = pbuffer.data(idx_kin_hd + 95);

    auto tk_yyyyz_xz = pbuffer.data(idx_kin_hd + 98);

    auto tk_yyyyz_yz = pbuffer.data(idx_kin_hd + 100);

    auto tk_yyyyz_zz = pbuffer.data(idx_kin_hd + 101);

    auto tk_yyyzz_xx = pbuffer.data(idx_kin_hd + 102);

    auto tk_yyyzz_xy = pbuffer.data(idx_kin_hd + 103);

    auto tk_yyyzz_xz = pbuffer.data(idx_kin_hd + 104);

    auto tk_yyyzz_yy = pbuffer.data(idx_kin_hd + 105);

    auto tk_yyyzz_yz = pbuffer.data(idx_kin_hd + 106);

    auto tk_yyyzz_zz = pbuffer.data(idx_kin_hd + 107);

    auto tk_yyzzz_xx = pbuffer.data(idx_kin_hd + 108);

    auto tk_yyzzz_xy = pbuffer.data(idx_kin_hd + 109);

    auto tk_yyzzz_xz = pbuffer.data(idx_kin_hd + 110);

    auto tk_yyzzz_yy = pbuffer.data(idx_kin_hd + 111);

    auto tk_yyzzz_yz = pbuffer.data(idx_kin_hd + 112);

    auto tk_yyzzz_zz = pbuffer.data(idx_kin_hd + 113);

    auto tk_yzzzz_xy = pbuffer.data(idx_kin_hd + 115);

    auto tk_yzzzz_xz = pbuffer.data(idx_kin_hd + 116);

    auto tk_yzzzz_yy = pbuffer.data(idx_kin_hd + 117);

    auto tk_yzzzz_yz = pbuffer.data(idx_kin_hd + 118);

    auto tk_yzzzz_zz = pbuffer.data(idx_kin_hd + 119);

    auto tk_zzzzz_xx = pbuffer.data(idx_kin_hd + 120);

    auto tk_zzzzz_xy = pbuffer.data(idx_kin_hd + 121);

    auto tk_zzzzz_xz = pbuffer.data(idx_kin_hd + 122);

    auto tk_zzzzz_yy = pbuffer.data(idx_kin_hd + 123);

    auto tk_zzzzz_yz = pbuffer.data(idx_kin_hd + 124);

    auto tk_zzzzz_zz = pbuffer.data(idx_kin_hd + 125);

    // Set up components of auxiliary buffer : HF

    auto tk_xxxxx_xxx = pbuffer.data(idx_kin_hf);

    auto tk_xxxxx_xxy = pbuffer.data(idx_kin_hf + 1);

    auto tk_xxxxx_xxz = pbuffer.data(idx_kin_hf + 2);

    auto tk_xxxxx_xyy = pbuffer.data(idx_kin_hf + 3);

    auto tk_xxxxx_xyz = pbuffer.data(idx_kin_hf + 4);

    auto tk_xxxxx_xzz = pbuffer.data(idx_kin_hf + 5);

    auto tk_xxxxx_yyy = pbuffer.data(idx_kin_hf + 6);

    auto tk_xxxxx_yyz = pbuffer.data(idx_kin_hf + 7);

    auto tk_xxxxx_yzz = pbuffer.data(idx_kin_hf + 8);

    auto tk_xxxxx_zzz = pbuffer.data(idx_kin_hf + 9);

    auto tk_xxxxy_xxx = pbuffer.data(idx_kin_hf + 10);

    auto tk_xxxxy_xxy = pbuffer.data(idx_kin_hf + 11);

    auto tk_xxxxy_xxz = pbuffer.data(idx_kin_hf + 12);

    auto tk_xxxxy_xyy = pbuffer.data(idx_kin_hf + 13);

    auto tk_xxxxy_xzz = pbuffer.data(idx_kin_hf + 15);

    auto tk_xxxxy_yyy = pbuffer.data(idx_kin_hf + 16);

    auto tk_xxxxz_xxx = pbuffer.data(idx_kin_hf + 20);

    auto tk_xxxxz_xxy = pbuffer.data(idx_kin_hf + 21);

    auto tk_xxxxz_xxz = pbuffer.data(idx_kin_hf + 22);

    auto tk_xxxxz_xyy = pbuffer.data(idx_kin_hf + 23);

    auto tk_xxxxz_xyz = pbuffer.data(idx_kin_hf + 24);

    auto tk_xxxxz_xzz = pbuffer.data(idx_kin_hf + 25);

    auto tk_xxxxz_yyz = pbuffer.data(idx_kin_hf + 27);

    auto tk_xxxxz_yzz = pbuffer.data(idx_kin_hf + 28);

    auto tk_xxxxz_zzz = pbuffer.data(idx_kin_hf + 29);

    auto tk_xxxyy_xxx = pbuffer.data(idx_kin_hf + 30);

    auto tk_xxxyy_xxy = pbuffer.data(idx_kin_hf + 31);

    auto tk_xxxyy_xxz = pbuffer.data(idx_kin_hf + 32);

    auto tk_xxxyy_xyy = pbuffer.data(idx_kin_hf + 33);

    auto tk_xxxyy_xyz = pbuffer.data(idx_kin_hf + 34);

    auto tk_xxxyy_xzz = pbuffer.data(idx_kin_hf + 35);

    auto tk_xxxyy_yyy = pbuffer.data(idx_kin_hf + 36);

    auto tk_xxxyy_yyz = pbuffer.data(idx_kin_hf + 37);

    auto tk_xxxyy_yzz = pbuffer.data(idx_kin_hf + 38);

    auto tk_xxxyy_zzz = pbuffer.data(idx_kin_hf + 39);

    auto tk_xxxzz_xxx = pbuffer.data(idx_kin_hf + 50);

    auto tk_xxxzz_xxy = pbuffer.data(idx_kin_hf + 51);

    auto tk_xxxzz_xxz = pbuffer.data(idx_kin_hf + 52);

    auto tk_xxxzz_xyy = pbuffer.data(idx_kin_hf + 53);

    auto tk_xxxzz_xyz = pbuffer.data(idx_kin_hf + 54);

    auto tk_xxxzz_xzz = pbuffer.data(idx_kin_hf + 55);

    auto tk_xxxzz_yyy = pbuffer.data(idx_kin_hf + 56);

    auto tk_xxxzz_yyz = pbuffer.data(idx_kin_hf + 57);

    auto tk_xxxzz_yzz = pbuffer.data(idx_kin_hf + 58);

    auto tk_xxxzz_zzz = pbuffer.data(idx_kin_hf + 59);

    auto tk_xxyyy_xxx = pbuffer.data(idx_kin_hf + 60);

    auto tk_xxyyy_xxy = pbuffer.data(idx_kin_hf + 61);

    auto tk_xxyyy_xxz = pbuffer.data(idx_kin_hf + 62);

    auto tk_xxyyy_xyy = pbuffer.data(idx_kin_hf + 63);

    auto tk_xxyyy_xyz = pbuffer.data(idx_kin_hf + 64);

    auto tk_xxyyy_xzz = pbuffer.data(idx_kin_hf + 65);

    auto tk_xxyyy_yyy = pbuffer.data(idx_kin_hf + 66);

    auto tk_xxyyy_yyz = pbuffer.data(idx_kin_hf + 67);

    auto tk_xxyyy_yzz = pbuffer.data(idx_kin_hf + 68);

    auto tk_xxyyy_zzz = pbuffer.data(idx_kin_hf + 69);

    auto tk_xxyyz_xxy = pbuffer.data(idx_kin_hf + 71);

    auto tk_xxyyz_xyy = pbuffer.data(idx_kin_hf + 73);

    auto tk_xxyzz_xxx = pbuffer.data(idx_kin_hf + 80);

    auto tk_xxyzz_xxz = pbuffer.data(idx_kin_hf + 82);

    auto tk_xxyzz_xzz = pbuffer.data(idx_kin_hf + 85);

    auto tk_xxzzz_xxx = pbuffer.data(idx_kin_hf + 90);

    auto tk_xxzzz_xxy = pbuffer.data(idx_kin_hf + 91);

    auto tk_xxzzz_xxz = pbuffer.data(idx_kin_hf + 92);

    auto tk_xxzzz_xyy = pbuffer.data(idx_kin_hf + 93);

    auto tk_xxzzz_xyz = pbuffer.data(idx_kin_hf + 94);

    auto tk_xxzzz_xzz = pbuffer.data(idx_kin_hf + 95);

    auto tk_xxzzz_yyy = pbuffer.data(idx_kin_hf + 96);

    auto tk_xxzzz_yyz = pbuffer.data(idx_kin_hf + 97);

    auto tk_xxzzz_yzz = pbuffer.data(idx_kin_hf + 98);

    auto tk_xxzzz_zzz = pbuffer.data(idx_kin_hf + 99);

    auto tk_xyyyy_xxx = pbuffer.data(idx_kin_hf + 100);

    auto tk_xyyyy_xxy = pbuffer.data(idx_kin_hf + 101);

    auto tk_xyyyy_xyy = pbuffer.data(idx_kin_hf + 103);

    auto tk_xyyyy_xyz = pbuffer.data(idx_kin_hf + 104);

    auto tk_xyyyy_yyy = pbuffer.data(idx_kin_hf + 106);

    auto tk_xyyyy_yyz = pbuffer.data(idx_kin_hf + 107);

    auto tk_xyyyy_yzz = pbuffer.data(idx_kin_hf + 108);

    auto tk_xyyyy_zzz = pbuffer.data(idx_kin_hf + 109);

    auto tk_xyyzz_xyz = pbuffer.data(idx_kin_hf + 124);

    auto tk_xyyzz_yyy = pbuffer.data(idx_kin_hf + 126);

    auto tk_xyyzz_yyz = pbuffer.data(idx_kin_hf + 127);

    auto tk_xyyzz_yzz = pbuffer.data(idx_kin_hf + 128);

    auto tk_xyyzz_zzz = pbuffer.data(idx_kin_hf + 129);

    auto tk_xzzzz_xxx = pbuffer.data(idx_kin_hf + 140);

    auto tk_xzzzz_xxz = pbuffer.data(idx_kin_hf + 142);

    auto tk_xzzzz_xyz = pbuffer.data(idx_kin_hf + 144);

    auto tk_xzzzz_xzz = pbuffer.data(idx_kin_hf + 145);

    auto tk_xzzzz_yyy = pbuffer.data(idx_kin_hf + 146);

    auto tk_xzzzz_yyz = pbuffer.data(idx_kin_hf + 147);

    auto tk_xzzzz_yzz = pbuffer.data(idx_kin_hf + 148);

    auto tk_xzzzz_zzz = pbuffer.data(idx_kin_hf + 149);

    auto tk_yyyyy_xxx = pbuffer.data(idx_kin_hf + 150);

    auto tk_yyyyy_xxy = pbuffer.data(idx_kin_hf + 151);

    auto tk_yyyyy_xxz = pbuffer.data(idx_kin_hf + 152);

    auto tk_yyyyy_xyy = pbuffer.data(idx_kin_hf + 153);

    auto tk_yyyyy_xyz = pbuffer.data(idx_kin_hf + 154);

    auto tk_yyyyy_xzz = pbuffer.data(idx_kin_hf + 155);

    auto tk_yyyyy_yyy = pbuffer.data(idx_kin_hf + 156);

    auto tk_yyyyy_yyz = pbuffer.data(idx_kin_hf + 157);

    auto tk_yyyyy_yzz = pbuffer.data(idx_kin_hf + 158);

    auto tk_yyyyy_zzz = pbuffer.data(idx_kin_hf + 159);

    auto tk_yyyyz_xxy = pbuffer.data(idx_kin_hf + 161);

    auto tk_yyyyz_xxz = pbuffer.data(idx_kin_hf + 162);

    auto tk_yyyyz_xyy = pbuffer.data(idx_kin_hf + 163);

    auto tk_yyyyz_xyz = pbuffer.data(idx_kin_hf + 164);

    auto tk_yyyyz_xzz = pbuffer.data(idx_kin_hf + 165);

    auto tk_yyyyz_yyy = pbuffer.data(idx_kin_hf + 166);

    auto tk_yyyyz_yyz = pbuffer.data(idx_kin_hf + 167);

    auto tk_yyyyz_yzz = pbuffer.data(idx_kin_hf + 168);

    auto tk_yyyyz_zzz = pbuffer.data(idx_kin_hf + 169);

    auto tk_yyyzz_xxx = pbuffer.data(idx_kin_hf + 170);

    auto tk_yyyzz_xxy = pbuffer.data(idx_kin_hf + 171);

    auto tk_yyyzz_xxz = pbuffer.data(idx_kin_hf + 172);

    auto tk_yyyzz_xyy = pbuffer.data(idx_kin_hf + 173);

    auto tk_yyyzz_xyz = pbuffer.data(idx_kin_hf + 174);

    auto tk_yyyzz_xzz = pbuffer.data(idx_kin_hf + 175);

    auto tk_yyyzz_yyy = pbuffer.data(idx_kin_hf + 176);

    auto tk_yyyzz_yyz = pbuffer.data(idx_kin_hf + 177);

    auto tk_yyyzz_yzz = pbuffer.data(idx_kin_hf + 178);

    auto tk_yyyzz_zzz = pbuffer.data(idx_kin_hf + 179);

    auto tk_yyzzz_xxx = pbuffer.data(idx_kin_hf + 180);

    auto tk_yyzzz_xxy = pbuffer.data(idx_kin_hf + 181);

    auto tk_yyzzz_xxz = pbuffer.data(idx_kin_hf + 182);

    auto tk_yyzzz_xyy = pbuffer.data(idx_kin_hf + 183);

    auto tk_yyzzz_xyz = pbuffer.data(idx_kin_hf + 184);

    auto tk_yyzzz_xzz = pbuffer.data(idx_kin_hf + 185);

    auto tk_yyzzz_yyy = pbuffer.data(idx_kin_hf + 186);

    auto tk_yyzzz_yyz = pbuffer.data(idx_kin_hf + 187);

    auto tk_yyzzz_yzz = pbuffer.data(idx_kin_hf + 188);

    auto tk_yyzzz_zzz = pbuffer.data(idx_kin_hf + 189);

    auto tk_yzzzz_xxx = pbuffer.data(idx_kin_hf + 190);

    auto tk_yzzzz_xxy = pbuffer.data(idx_kin_hf + 191);

    auto tk_yzzzz_xxz = pbuffer.data(idx_kin_hf + 192);

    auto tk_yzzzz_xyy = pbuffer.data(idx_kin_hf + 193);

    auto tk_yzzzz_xyz = pbuffer.data(idx_kin_hf + 194);

    auto tk_yzzzz_xzz = pbuffer.data(idx_kin_hf + 195);

    auto tk_yzzzz_yyy = pbuffer.data(idx_kin_hf + 196);

    auto tk_yzzzz_yyz = pbuffer.data(idx_kin_hf + 197);

    auto tk_yzzzz_yzz = pbuffer.data(idx_kin_hf + 198);

    auto tk_yzzzz_zzz = pbuffer.data(idx_kin_hf + 199);

    auto tk_zzzzz_xxx = pbuffer.data(idx_kin_hf + 200);

    auto tk_zzzzz_xxy = pbuffer.data(idx_kin_hf + 201);

    auto tk_zzzzz_xxz = pbuffer.data(idx_kin_hf + 202);

    auto tk_zzzzz_xyy = pbuffer.data(idx_kin_hf + 203);

    auto tk_zzzzz_xyz = pbuffer.data(idx_kin_hf + 204);

    auto tk_zzzzz_xzz = pbuffer.data(idx_kin_hf + 205);

    auto tk_zzzzz_yyy = pbuffer.data(idx_kin_hf + 206);

    auto tk_zzzzz_yyz = pbuffer.data(idx_kin_hf + 207);

    auto tk_zzzzz_yzz = pbuffer.data(idx_kin_hf + 208);

    auto tk_zzzzz_zzz = pbuffer.data(idx_kin_hf + 209);

    // Set up components of auxiliary buffer : IF

    auto ts_xxxxxx_xxx = pbuffer.data(idx_ovl_if);

    auto ts_xxxxxx_xxy = pbuffer.data(idx_ovl_if + 1);

    auto ts_xxxxxx_xxz = pbuffer.data(idx_ovl_if + 2);

    auto ts_xxxxxx_xyy = pbuffer.data(idx_ovl_if + 3);

    auto ts_xxxxxx_xyz = pbuffer.data(idx_ovl_if + 4);

    auto ts_xxxxxx_xzz = pbuffer.data(idx_ovl_if + 5);

    auto ts_xxxxxx_yyy = pbuffer.data(idx_ovl_if + 6);

    auto ts_xxxxxx_yyz = pbuffer.data(idx_ovl_if + 7);

    auto ts_xxxxxx_yzz = pbuffer.data(idx_ovl_if + 8);

    auto ts_xxxxxx_zzz = pbuffer.data(idx_ovl_if + 9);

    auto ts_xxxxxy_xxx = pbuffer.data(idx_ovl_if + 10);

    auto ts_xxxxxy_xxy = pbuffer.data(idx_ovl_if + 11);

    auto ts_xxxxxy_xxz = pbuffer.data(idx_ovl_if + 12);

    auto ts_xxxxxy_xyy = pbuffer.data(idx_ovl_if + 13);

    auto ts_xxxxxy_xyz = pbuffer.data(idx_ovl_if + 14);

    auto ts_xxxxxy_xzz = pbuffer.data(idx_ovl_if + 15);

    auto ts_xxxxxy_yyy = pbuffer.data(idx_ovl_if + 16);

    auto ts_xxxxxy_yyz = pbuffer.data(idx_ovl_if + 17);

    auto ts_xxxxxy_yzz = pbuffer.data(idx_ovl_if + 18);

    auto ts_xxxxxy_zzz = pbuffer.data(idx_ovl_if + 19);

    auto ts_xxxxxz_xxx = pbuffer.data(idx_ovl_if + 20);

    auto ts_xxxxxz_xxy = pbuffer.data(idx_ovl_if + 21);

    auto ts_xxxxxz_xxz = pbuffer.data(idx_ovl_if + 22);

    auto ts_xxxxxz_xyy = pbuffer.data(idx_ovl_if + 23);

    auto ts_xxxxxz_xyz = pbuffer.data(idx_ovl_if + 24);

    auto ts_xxxxxz_xzz = pbuffer.data(idx_ovl_if + 25);

    auto ts_xxxxxz_yyy = pbuffer.data(idx_ovl_if + 26);

    auto ts_xxxxxz_yyz = pbuffer.data(idx_ovl_if + 27);

    auto ts_xxxxxz_yzz = pbuffer.data(idx_ovl_if + 28);

    auto ts_xxxxxz_zzz = pbuffer.data(idx_ovl_if + 29);

    auto ts_xxxxyy_xxx = pbuffer.data(idx_ovl_if + 30);

    auto ts_xxxxyy_xxy = pbuffer.data(idx_ovl_if + 31);

    auto ts_xxxxyy_xxz = pbuffer.data(idx_ovl_if + 32);

    auto ts_xxxxyy_xyy = pbuffer.data(idx_ovl_if + 33);

    auto ts_xxxxyy_xyz = pbuffer.data(idx_ovl_if + 34);

    auto ts_xxxxyy_xzz = pbuffer.data(idx_ovl_if + 35);

    auto ts_xxxxyy_yyy = pbuffer.data(idx_ovl_if + 36);

    auto ts_xxxxyy_yyz = pbuffer.data(idx_ovl_if + 37);

    auto ts_xxxxyy_yzz = pbuffer.data(idx_ovl_if + 38);

    auto ts_xxxxyy_zzz = pbuffer.data(idx_ovl_if + 39);

    auto ts_xxxxyz_xxx = pbuffer.data(idx_ovl_if + 40);

    auto ts_xxxxyz_xxy = pbuffer.data(idx_ovl_if + 41);

    auto ts_xxxxyz_xxz = pbuffer.data(idx_ovl_if + 42);

    auto ts_xxxxyz_xyy = pbuffer.data(idx_ovl_if + 43);

    auto ts_xxxxyz_xyz = pbuffer.data(idx_ovl_if + 44);

    auto ts_xxxxyz_xzz = pbuffer.data(idx_ovl_if + 45);

    auto ts_xxxxyz_yyy = pbuffer.data(idx_ovl_if + 46);

    auto ts_xxxxyz_yyz = pbuffer.data(idx_ovl_if + 47);

    auto ts_xxxxyz_yzz = pbuffer.data(idx_ovl_if + 48);

    auto ts_xxxxyz_zzz = pbuffer.data(idx_ovl_if + 49);

    auto ts_xxxxzz_xxx = pbuffer.data(idx_ovl_if + 50);

    auto ts_xxxxzz_xxy = pbuffer.data(idx_ovl_if + 51);

    auto ts_xxxxzz_xxz = pbuffer.data(idx_ovl_if + 52);

    auto ts_xxxxzz_xyy = pbuffer.data(idx_ovl_if + 53);

    auto ts_xxxxzz_xyz = pbuffer.data(idx_ovl_if + 54);

    auto ts_xxxxzz_xzz = pbuffer.data(idx_ovl_if + 55);

    auto ts_xxxxzz_yyy = pbuffer.data(idx_ovl_if + 56);

    auto ts_xxxxzz_yyz = pbuffer.data(idx_ovl_if + 57);

    auto ts_xxxxzz_yzz = pbuffer.data(idx_ovl_if + 58);

    auto ts_xxxxzz_zzz = pbuffer.data(idx_ovl_if + 59);

    auto ts_xxxyyy_xxx = pbuffer.data(idx_ovl_if + 60);

    auto ts_xxxyyy_xxy = pbuffer.data(idx_ovl_if + 61);

    auto ts_xxxyyy_xxz = pbuffer.data(idx_ovl_if + 62);

    auto ts_xxxyyy_xyy = pbuffer.data(idx_ovl_if + 63);

    auto ts_xxxyyy_xyz = pbuffer.data(idx_ovl_if + 64);

    auto ts_xxxyyy_xzz = pbuffer.data(idx_ovl_if + 65);

    auto ts_xxxyyy_yyy = pbuffer.data(idx_ovl_if + 66);

    auto ts_xxxyyy_yyz = pbuffer.data(idx_ovl_if + 67);

    auto ts_xxxyyy_yzz = pbuffer.data(idx_ovl_if + 68);

    auto ts_xxxyyy_zzz = pbuffer.data(idx_ovl_if + 69);

    auto ts_xxxyyz_xxx = pbuffer.data(idx_ovl_if + 70);

    auto ts_xxxyyz_xxy = pbuffer.data(idx_ovl_if + 71);

    auto ts_xxxyyz_xxz = pbuffer.data(idx_ovl_if + 72);

    auto ts_xxxyyz_xyy = pbuffer.data(idx_ovl_if + 73);

    auto ts_xxxyyz_xyz = pbuffer.data(idx_ovl_if + 74);

    auto ts_xxxyyz_xzz = pbuffer.data(idx_ovl_if + 75);

    auto ts_xxxyyz_yyy = pbuffer.data(idx_ovl_if + 76);

    auto ts_xxxyyz_yyz = pbuffer.data(idx_ovl_if + 77);

    auto ts_xxxyyz_yzz = pbuffer.data(idx_ovl_if + 78);

    auto ts_xxxyyz_zzz = pbuffer.data(idx_ovl_if + 79);

    auto ts_xxxyzz_xxx = pbuffer.data(idx_ovl_if + 80);

    auto ts_xxxyzz_xxy = pbuffer.data(idx_ovl_if + 81);

    auto ts_xxxyzz_xxz = pbuffer.data(idx_ovl_if + 82);

    auto ts_xxxyzz_xyy = pbuffer.data(idx_ovl_if + 83);

    auto ts_xxxyzz_xyz = pbuffer.data(idx_ovl_if + 84);

    auto ts_xxxyzz_xzz = pbuffer.data(idx_ovl_if + 85);

    auto ts_xxxyzz_yyy = pbuffer.data(idx_ovl_if + 86);

    auto ts_xxxyzz_yyz = pbuffer.data(idx_ovl_if + 87);

    auto ts_xxxyzz_yzz = pbuffer.data(idx_ovl_if + 88);

    auto ts_xxxyzz_zzz = pbuffer.data(idx_ovl_if + 89);

    auto ts_xxxzzz_xxx = pbuffer.data(idx_ovl_if + 90);

    auto ts_xxxzzz_xxy = pbuffer.data(idx_ovl_if + 91);

    auto ts_xxxzzz_xxz = pbuffer.data(idx_ovl_if + 92);

    auto ts_xxxzzz_xyy = pbuffer.data(idx_ovl_if + 93);

    auto ts_xxxzzz_xyz = pbuffer.data(idx_ovl_if + 94);

    auto ts_xxxzzz_xzz = pbuffer.data(idx_ovl_if + 95);

    auto ts_xxxzzz_yyy = pbuffer.data(idx_ovl_if + 96);

    auto ts_xxxzzz_yyz = pbuffer.data(idx_ovl_if + 97);

    auto ts_xxxzzz_yzz = pbuffer.data(idx_ovl_if + 98);

    auto ts_xxxzzz_zzz = pbuffer.data(idx_ovl_if + 99);

    auto ts_xxyyyy_xxx = pbuffer.data(idx_ovl_if + 100);

    auto ts_xxyyyy_xxy = pbuffer.data(idx_ovl_if + 101);

    auto ts_xxyyyy_xxz = pbuffer.data(idx_ovl_if + 102);

    auto ts_xxyyyy_xyy = pbuffer.data(idx_ovl_if + 103);

    auto ts_xxyyyy_xyz = pbuffer.data(idx_ovl_if + 104);

    auto ts_xxyyyy_xzz = pbuffer.data(idx_ovl_if + 105);

    auto ts_xxyyyy_yyy = pbuffer.data(idx_ovl_if + 106);

    auto ts_xxyyyy_yyz = pbuffer.data(idx_ovl_if + 107);

    auto ts_xxyyyy_yzz = pbuffer.data(idx_ovl_if + 108);

    auto ts_xxyyyy_zzz = pbuffer.data(idx_ovl_if + 109);

    auto ts_xxyyyz_xxx = pbuffer.data(idx_ovl_if + 110);

    auto ts_xxyyyz_xxy = pbuffer.data(idx_ovl_if + 111);

    auto ts_xxyyyz_xxz = pbuffer.data(idx_ovl_if + 112);

    auto ts_xxyyyz_xyy = pbuffer.data(idx_ovl_if + 113);

    auto ts_xxyyyz_xyz = pbuffer.data(idx_ovl_if + 114);

    auto ts_xxyyyz_xzz = pbuffer.data(idx_ovl_if + 115);

    auto ts_xxyyyz_yyy = pbuffer.data(idx_ovl_if + 116);

    auto ts_xxyyyz_yyz = pbuffer.data(idx_ovl_if + 117);

    auto ts_xxyyyz_yzz = pbuffer.data(idx_ovl_if + 118);

    auto ts_xxyyyz_zzz = pbuffer.data(idx_ovl_if + 119);

    auto ts_xxyyzz_xxx = pbuffer.data(idx_ovl_if + 120);

    auto ts_xxyyzz_xxy = pbuffer.data(idx_ovl_if + 121);

    auto ts_xxyyzz_xxz = pbuffer.data(idx_ovl_if + 122);

    auto ts_xxyyzz_xyy = pbuffer.data(idx_ovl_if + 123);

    auto ts_xxyyzz_xyz = pbuffer.data(idx_ovl_if + 124);

    auto ts_xxyyzz_xzz = pbuffer.data(idx_ovl_if + 125);

    auto ts_xxyyzz_yyy = pbuffer.data(idx_ovl_if + 126);

    auto ts_xxyyzz_yyz = pbuffer.data(idx_ovl_if + 127);

    auto ts_xxyyzz_yzz = pbuffer.data(idx_ovl_if + 128);

    auto ts_xxyyzz_zzz = pbuffer.data(idx_ovl_if + 129);

    auto ts_xxyzzz_xxx = pbuffer.data(idx_ovl_if + 130);

    auto ts_xxyzzz_xxy = pbuffer.data(idx_ovl_if + 131);

    auto ts_xxyzzz_xxz = pbuffer.data(idx_ovl_if + 132);

    auto ts_xxyzzz_xyy = pbuffer.data(idx_ovl_if + 133);

    auto ts_xxyzzz_xyz = pbuffer.data(idx_ovl_if + 134);

    auto ts_xxyzzz_xzz = pbuffer.data(idx_ovl_if + 135);

    auto ts_xxyzzz_yyy = pbuffer.data(idx_ovl_if + 136);

    auto ts_xxyzzz_yyz = pbuffer.data(idx_ovl_if + 137);

    auto ts_xxyzzz_yzz = pbuffer.data(idx_ovl_if + 138);

    auto ts_xxyzzz_zzz = pbuffer.data(idx_ovl_if + 139);

    auto ts_xxzzzz_xxx = pbuffer.data(idx_ovl_if + 140);

    auto ts_xxzzzz_xxy = pbuffer.data(idx_ovl_if + 141);

    auto ts_xxzzzz_xxz = pbuffer.data(idx_ovl_if + 142);

    auto ts_xxzzzz_xyy = pbuffer.data(idx_ovl_if + 143);

    auto ts_xxzzzz_xyz = pbuffer.data(idx_ovl_if + 144);

    auto ts_xxzzzz_xzz = pbuffer.data(idx_ovl_if + 145);

    auto ts_xxzzzz_yyy = pbuffer.data(idx_ovl_if + 146);

    auto ts_xxzzzz_yyz = pbuffer.data(idx_ovl_if + 147);

    auto ts_xxzzzz_yzz = pbuffer.data(idx_ovl_if + 148);

    auto ts_xxzzzz_zzz = pbuffer.data(idx_ovl_if + 149);

    auto ts_xyyyyy_xxx = pbuffer.data(idx_ovl_if + 150);

    auto ts_xyyyyy_xxy = pbuffer.data(idx_ovl_if + 151);

    auto ts_xyyyyy_xxz = pbuffer.data(idx_ovl_if + 152);

    auto ts_xyyyyy_xyy = pbuffer.data(idx_ovl_if + 153);

    auto ts_xyyyyy_xyz = pbuffer.data(idx_ovl_if + 154);

    auto ts_xyyyyy_xzz = pbuffer.data(idx_ovl_if + 155);

    auto ts_xyyyyy_yyy = pbuffer.data(idx_ovl_if + 156);

    auto ts_xyyyyy_yyz = pbuffer.data(idx_ovl_if + 157);

    auto ts_xyyyyy_yzz = pbuffer.data(idx_ovl_if + 158);

    auto ts_xyyyyy_zzz = pbuffer.data(idx_ovl_if + 159);

    auto ts_xyyyyz_xxx = pbuffer.data(idx_ovl_if + 160);

    auto ts_xyyyyz_xxy = pbuffer.data(idx_ovl_if + 161);

    auto ts_xyyyyz_xxz = pbuffer.data(idx_ovl_if + 162);

    auto ts_xyyyyz_xyy = pbuffer.data(idx_ovl_if + 163);

    auto ts_xyyyyz_xyz = pbuffer.data(idx_ovl_if + 164);

    auto ts_xyyyyz_xzz = pbuffer.data(idx_ovl_if + 165);

    auto ts_xyyyyz_yyy = pbuffer.data(idx_ovl_if + 166);

    auto ts_xyyyyz_yyz = pbuffer.data(idx_ovl_if + 167);

    auto ts_xyyyyz_yzz = pbuffer.data(idx_ovl_if + 168);

    auto ts_xyyyyz_zzz = pbuffer.data(idx_ovl_if + 169);

    auto ts_xyyyzz_xxx = pbuffer.data(idx_ovl_if + 170);

    auto ts_xyyyzz_xxy = pbuffer.data(idx_ovl_if + 171);

    auto ts_xyyyzz_xxz = pbuffer.data(idx_ovl_if + 172);

    auto ts_xyyyzz_xyy = pbuffer.data(idx_ovl_if + 173);

    auto ts_xyyyzz_xyz = pbuffer.data(idx_ovl_if + 174);

    auto ts_xyyyzz_xzz = pbuffer.data(idx_ovl_if + 175);

    auto ts_xyyyzz_yyy = pbuffer.data(idx_ovl_if + 176);

    auto ts_xyyyzz_yyz = pbuffer.data(idx_ovl_if + 177);

    auto ts_xyyyzz_yzz = pbuffer.data(idx_ovl_if + 178);

    auto ts_xyyyzz_zzz = pbuffer.data(idx_ovl_if + 179);

    auto ts_xyyzzz_xxx = pbuffer.data(idx_ovl_if + 180);

    auto ts_xyyzzz_xxy = pbuffer.data(idx_ovl_if + 181);

    auto ts_xyyzzz_xxz = pbuffer.data(idx_ovl_if + 182);

    auto ts_xyyzzz_xyy = pbuffer.data(idx_ovl_if + 183);

    auto ts_xyyzzz_xyz = pbuffer.data(idx_ovl_if + 184);

    auto ts_xyyzzz_xzz = pbuffer.data(idx_ovl_if + 185);

    auto ts_xyyzzz_yyy = pbuffer.data(idx_ovl_if + 186);

    auto ts_xyyzzz_yyz = pbuffer.data(idx_ovl_if + 187);

    auto ts_xyyzzz_yzz = pbuffer.data(idx_ovl_if + 188);

    auto ts_xyyzzz_zzz = pbuffer.data(idx_ovl_if + 189);

    auto ts_xyzzzz_xxx = pbuffer.data(idx_ovl_if + 190);

    auto ts_xyzzzz_xxy = pbuffer.data(idx_ovl_if + 191);

    auto ts_xyzzzz_xxz = pbuffer.data(idx_ovl_if + 192);

    auto ts_xyzzzz_xyy = pbuffer.data(idx_ovl_if + 193);

    auto ts_xyzzzz_xyz = pbuffer.data(idx_ovl_if + 194);

    auto ts_xyzzzz_xzz = pbuffer.data(idx_ovl_if + 195);

    auto ts_xyzzzz_yyy = pbuffer.data(idx_ovl_if + 196);

    auto ts_xyzzzz_yyz = pbuffer.data(idx_ovl_if + 197);

    auto ts_xyzzzz_yzz = pbuffer.data(idx_ovl_if + 198);

    auto ts_xyzzzz_zzz = pbuffer.data(idx_ovl_if + 199);

    auto ts_xzzzzz_xxx = pbuffer.data(idx_ovl_if + 200);

    auto ts_xzzzzz_xxy = pbuffer.data(idx_ovl_if + 201);

    auto ts_xzzzzz_xxz = pbuffer.data(idx_ovl_if + 202);

    auto ts_xzzzzz_xyy = pbuffer.data(idx_ovl_if + 203);

    auto ts_xzzzzz_xyz = pbuffer.data(idx_ovl_if + 204);

    auto ts_xzzzzz_xzz = pbuffer.data(idx_ovl_if + 205);

    auto ts_xzzzzz_yyy = pbuffer.data(idx_ovl_if + 206);

    auto ts_xzzzzz_yyz = pbuffer.data(idx_ovl_if + 207);

    auto ts_xzzzzz_yzz = pbuffer.data(idx_ovl_if + 208);

    auto ts_xzzzzz_zzz = pbuffer.data(idx_ovl_if + 209);

    auto ts_yyyyyy_xxx = pbuffer.data(idx_ovl_if + 210);

    auto ts_yyyyyy_xxy = pbuffer.data(idx_ovl_if + 211);

    auto ts_yyyyyy_xxz = pbuffer.data(idx_ovl_if + 212);

    auto ts_yyyyyy_xyy = pbuffer.data(idx_ovl_if + 213);

    auto ts_yyyyyy_xyz = pbuffer.data(idx_ovl_if + 214);

    auto ts_yyyyyy_xzz = pbuffer.data(idx_ovl_if + 215);

    auto ts_yyyyyy_yyy = pbuffer.data(idx_ovl_if + 216);

    auto ts_yyyyyy_yyz = pbuffer.data(idx_ovl_if + 217);

    auto ts_yyyyyy_yzz = pbuffer.data(idx_ovl_if + 218);

    auto ts_yyyyyy_zzz = pbuffer.data(idx_ovl_if + 219);

    auto ts_yyyyyz_xxx = pbuffer.data(idx_ovl_if + 220);

    auto ts_yyyyyz_xxy = pbuffer.data(idx_ovl_if + 221);

    auto ts_yyyyyz_xxz = pbuffer.data(idx_ovl_if + 222);

    auto ts_yyyyyz_xyy = pbuffer.data(idx_ovl_if + 223);

    auto ts_yyyyyz_xyz = pbuffer.data(idx_ovl_if + 224);

    auto ts_yyyyyz_xzz = pbuffer.data(idx_ovl_if + 225);

    auto ts_yyyyyz_yyy = pbuffer.data(idx_ovl_if + 226);

    auto ts_yyyyyz_yyz = pbuffer.data(idx_ovl_if + 227);

    auto ts_yyyyyz_yzz = pbuffer.data(idx_ovl_if + 228);

    auto ts_yyyyyz_zzz = pbuffer.data(idx_ovl_if + 229);

    auto ts_yyyyzz_xxx = pbuffer.data(idx_ovl_if + 230);

    auto ts_yyyyzz_xxy = pbuffer.data(idx_ovl_if + 231);

    auto ts_yyyyzz_xxz = pbuffer.data(idx_ovl_if + 232);

    auto ts_yyyyzz_xyy = pbuffer.data(idx_ovl_if + 233);

    auto ts_yyyyzz_xyz = pbuffer.data(idx_ovl_if + 234);

    auto ts_yyyyzz_xzz = pbuffer.data(idx_ovl_if + 235);

    auto ts_yyyyzz_yyy = pbuffer.data(idx_ovl_if + 236);

    auto ts_yyyyzz_yyz = pbuffer.data(idx_ovl_if + 237);

    auto ts_yyyyzz_yzz = pbuffer.data(idx_ovl_if + 238);

    auto ts_yyyyzz_zzz = pbuffer.data(idx_ovl_if + 239);

    auto ts_yyyzzz_xxx = pbuffer.data(idx_ovl_if + 240);

    auto ts_yyyzzz_xxy = pbuffer.data(idx_ovl_if + 241);

    auto ts_yyyzzz_xxz = pbuffer.data(idx_ovl_if + 242);

    auto ts_yyyzzz_xyy = pbuffer.data(idx_ovl_if + 243);

    auto ts_yyyzzz_xyz = pbuffer.data(idx_ovl_if + 244);

    auto ts_yyyzzz_xzz = pbuffer.data(idx_ovl_if + 245);

    auto ts_yyyzzz_yyy = pbuffer.data(idx_ovl_if + 246);

    auto ts_yyyzzz_yyz = pbuffer.data(idx_ovl_if + 247);

    auto ts_yyyzzz_yzz = pbuffer.data(idx_ovl_if + 248);

    auto ts_yyyzzz_zzz = pbuffer.data(idx_ovl_if + 249);

    auto ts_yyzzzz_xxx = pbuffer.data(idx_ovl_if + 250);

    auto ts_yyzzzz_xxy = pbuffer.data(idx_ovl_if + 251);

    auto ts_yyzzzz_xxz = pbuffer.data(idx_ovl_if + 252);

    auto ts_yyzzzz_xyy = pbuffer.data(idx_ovl_if + 253);

    auto ts_yyzzzz_xyz = pbuffer.data(idx_ovl_if + 254);

    auto ts_yyzzzz_xzz = pbuffer.data(idx_ovl_if + 255);

    auto ts_yyzzzz_yyy = pbuffer.data(idx_ovl_if + 256);

    auto ts_yyzzzz_yyz = pbuffer.data(idx_ovl_if + 257);

    auto ts_yyzzzz_yzz = pbuffer.data(idx_ovl_if + 258);

    auto ts_yyzzzz_zzz = pbuffer.data(idx_ovl_if + 259);

    auto ts_yzzzzz_xxx = pbuffer.data(idx_ovl_if + 260);

    auto ts_yzzzzz_xxy = pbuffer.data(idx_ovl_if + 261);

    auto ts_yzzzzz_xxz = pbuffer.data(idx_ovl_if + 262);

    auto ts_yzzzzz_xyy = pbuffer.data(idx_ovl_if + 263);

    auto ts_yzzzzz_xyz = pbuffer.data(idx_ovl_if + 264);

    auto ts_yzzzzz_xzz = pbuffer.data(idx_ovl_if + 265);

    auto ts_yzzzzz_yyy = pbuffer.data(idx_ovl_if + 266);

    auto ts_yzzzzz_yyz = pbuffer.data(idx_ovl_if + 267);

    auto ts_yzzzzz_yzz = pbuffer.data(idx_ovl_if + 268);

    auto ts_yzzzzz_zzz = pbuffer.data(idx_ovl_if + 269);

    auto ts_zzzzzz_xxx = pbuffer.data(idx_ovl_if + 270);

    auto ts_zzzzzz_xxy = pbuffer.data(idx_ovl_if + 271);

    auto ts_zzzzzz_xxz = pbuffer.data(idx_ovl_if + 272);

    auto ts_zzzzzz_xyy = pbuffer.data(idx_ovl_if + 273);

    auto ts_zzzzzz_xyz = pbuffer.data(idx_ovl_if + 274);

    auto ts_zzzzzz_xzz = pbuffer.data(idx_ovl_if + 275);

    auto ts_zzzzzz_yyy = pbuffer.data(idx_ovl_if + 276);

    auto ts_zzzzzz_yyz = pbuffer.data(idx_ovl_if + 277);

    auto ts_zzzzzz_yzz = pbuffer.data(idx_ovl_if + 278);

    auto ts_zzzzzz_zzz = pbuffer.data(idx_ovl_if + 279);

    // Set up 0-10 components of targeted buffer : IF

    auto tk_xxxxxx_xxx = pbuffer.data(idx_kin_if);

    auto tk_xxxxxx_xxy = pbuffer.data(idx_kin_if + 1);

    auto tk_xxxxxx_xxz = pbuffer.data(idx_kin_if + 2);

    auto tk_xxxxxx_xyy = pbuffer.data(idx_kin_if + 3);

    auto tk_xxxxxx_xyz = pbuffer.data(idx_kin_if + 4);

    auto tk_xxxxxx_xzz = pbuffer.data(idx_kin_if + 5);

    auto tk_xxxxxx_yyy = pbuffer.data(idx_kin_if + 6);

    auto tk_xxxxxx_yyz = pbuffer.data(idx_kin_if + 7);

    auto tk_xxxxxx_yzz = pbuffer.data(idx_kin_if + 8);

    auto tk_xxxxxx_zzz = pbuffer.data(idx_kin_if + 9);

#pragma omp simd aligned(pa_x,              \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxy,   \
                             tk_xxxx_xxz,   \
                             tk_xxxx_xyy,   \
                             tk_xxxx_xyz,   \
                             tk_xxxx_xzz,   \
                             tk_xxxx_yyy,   \
                             tk_xxxx_yyz,   \
                             tk_xxxx_yzz,   \
                             tk_xxxx_zzz,   \
                             tk_xxxxx_xx,   \
                             tk_xxxxx_xxx,  \
                             tk_xxxxx_xxy,  \
                             tk_xxxxx_xxz,  \
                             tk_xxxxx_xy,   \
                             tk_xxxxx_xyy,  \
                             tk_xxxxx_xyz,  \
                             tk_xxxxx_xz,   \
                             tk_xxxxx_xzz,  \
                             tk_xxxxx_yy,   \
                             tk_xxxxx_yyy,  \
                             tk_xxxxx_yyz,  \
                             tk_xxxxx_yz,   \
                             tk_xxxxx_yzz,  \
                             tk_xxxxx_zz,   \
                             tk_xxxxx_zzz,  \
                             tk_xxxxxx_xxx, \
                             tk_xxxxxx_xxy, \
                             tk_xxxxxx_xxz, \
                             tk_xxxxxx_xyy, \
                             tk_xxxxxx_xyz, \
                             tk_xxxxxx_xzz, \
                             tk_xxxxxx_yyy, \
                             tk_xxxxxx_yyz, \
                             tk_xxxxxx_yzz, \
                             tk_xxxxxx_zzz, \
                             ts_xxxx_xxx,   \
                             ts_xxxx_xxy,   \
                             ts_xxxx_xxz,   \
                             ts_xxxx_xyy,   \
                             ts_xxxx_xyz,   \
                             ts_xxxx_xzz,   \
                             ts_xxxx_yyy,   \
                             ts_xxxx_yyz,   \
                             ts_xxxx_yzz,   \
                             ts_xxxx_zzz,   \
                             ts_xxxxxx_xxx, \
                             ts_xxxxxx_xxy, \
                             ts_xxxxxx_xxz, \
                             ts_xxxxxx_xyy, \
                             ts_xxxxxx_xyz, \
                             ts_xxxxxx_xzz, \
                             ts_xxxxxx_yyy, \
                             ts_xxxxxx_yyz, \
                             ts_xxxxxx_yzz, \
                             ts_xxxxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_xxx[i] = -10.0 * ts_xxxx_xxx[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxx[i] * fe_0 + 3.0 * tk_xxxxx_xx[i] * fe_0 +
                           tk_xxxxx_xxx[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxx[i] * fz_0;

        tk_xxxxxx_xxy[i] = -10.0 * ts_xxxx_xxy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxy[i] * fe_0 + 2.0 * tk_xxxxx_xy[i] * fe_0 +
                           tk_xxxxx_xxy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxy[i] * fz_0;

        tk_xxxxxx_xxz[i] = -10.0 * ts_xxxx_xxz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxz[i] * fe_0 + 2.0 * tk_xxxxx_xz[i] * fe_0 +
                           tk_xxxxx_xxz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxz[i] * fz_0;

        tk_xxxxxx_xyy[i] = -10.0 * ts_xxxx_xyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyy[i] * fe_0 + tk_xxxxx_yy[i] * fe_0 + tk_xxxxx_xyy[i] * pa_x[i] +
                           2.0 * ts_xxxxxx_xyy[i] * fz_0;

        tk_xxxxxx_xyz[i] = -10.0 * ts_xxxx_xyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyz[i] * fe_0 + tk_xxxxx_yz[i] * fe_0 + tk_xxxxx_xyz[i] * pa_x[i] +
                           2.0 * ts_xxxxxx_xyz[i] * fz_0;

        tk_xxxxxx_xzz[i] = -10.0 * ts_xxxx_xzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xzz[i] * fe_0 + tk_xxxxx_zz[i] * fe_0 + tk_xxxxx_xzz[i] * pa_x[i] +
                           2.0 * ts_xxxxxx_xzz[i] * fz_0;

        tk_xxxxxx_yyy[i] =
            -10.0 * ts_xxxx_yyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyy[i] * fe_0 + tk_xxxxx_yyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyy[i] * fz_0;

        tk_xxxxxx_yyz[i] =
            -10.0 * ts_xxxx_yyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyz[i] * fe_0 + tk_xxxxx_yyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyz[i] * fz_0;

        tk_xxxxxx_yzz[i] =
            -10.0 * ts_xxxx_yzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yzz[i] * fe_0 + tk_xxxxx_yzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yzz[i] * fz_0;

        tk_xxxxxx_zzz[i] =
            -10.0 * ts_xxxx_zzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_zzz[i] * fe_0 + tk_xxxxx_zzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : IF

    auto tk_xxxxxy_xxx = pbuffer.data(idx_kin_if + 10);

    auto tk_xxxxxy_xxy = pbuffer.data(idx_kin_if + 11);

    auto tk_xxxxxy_xxz = pbuffer.data(idx_kin_if + 12);

    auto tk_xxxxxy_xyy = pbuffer.data(idx_kin_if + 13);

    auto tk_xxxxxy_xyz = pbuffer.data(idx_kin_if + 14);

    auto tk_xxxxxy_xzz = pbuffer.data(idx_kin_if + 15);

    auto tk_xxxxxy_yyy = pbuffer.data(idx_kin_if + 16);

    auto tk_xxxxxy_yyz = pbuffer.data(idx_kin_if + 17);

    auto tk_xxxxxy_yzz = pbuffer.data(idx_kin_if + 18);

    auto tk_xxxxxy_zzz = pbuffer.data(idx_kin_if + 19);

#pragma omp simd aligned(pa_y,              \
                             tk_xxxxx_xx,   \
                             tk_xxxxx_xxx,  \
                             tk_xxxxx_xxy,  \
                             tk_xxxxx_xxz,  \
                             tk_xxxxx_xy,   \
                             tk_xxxxx_xyy,  \
                             tk_xxxxx_xyz,  \
                             tk_xxxxx_xz,   \
                             tk_xxxxx_xzz,  \
                             tk_xxxxx_yy,   \
                             tk_xxxxx_yyy,  \
                             tk_xxxxx_yyz,  \
                             tk_xxxxx_yz,   \
                             tk_xxxxx_yzz,  \
                             tk_xxxxx_zz,   \
                             tk_xxxxx_zzz,  \
                             tk_xxxxxy_xxx, \
                             tk_xxxxxy_xxy, \
                             tk_xxxxxy_xxz, \
                             tk_xxxxxy_xyy, \
                             tk_xxxxxy_xyz, \
                             tk_xxxxxy_xzz, \
                             tk_xxxxxy_yyy, \
                             tk_xxxxxy_yyz, \
                             tk_xxxxxy_yzz, \
                             tk_xxxxxy_zzz, \
                             ts_xxxxxy_xxx, \
                             ts_xxxxxy_xxy, \
                             ts_xxxxxy_xxz, \
                             ts_xxxxxy_xyy, \
                             ts_xxxxxy_xyz, \
                             ts_xxxxxy_xzz, \
                             ts_xxxxxy_yyy, \
                             ts_xxxxxy_yyz, \
                             ts_xxxxxy_yzz, \
                             ts_xxxxxy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_xxx[i] = tk_xxxxx_xxx[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxx[i] * fz_0;

        tk_xxxxxy_xxy[i] = tk_xxxxx_xx[i] * fe_0 + tk_xxxxx_xxy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxy[i] * fz_0;

        tk_xxxxxy_xxz[i] = tk_xxxxx_xxz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxz[i] * fz_0;

        tk_xxxxxy_xyy[i] = 2.0 * tk_xxxxx_xy[i] * fe_0 + tk_xxxxx_xyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyy[i] * fz_0;

        tk_xxxxxy_xyz[i] = tk_xxxxx_xz[i] * fe_0 + tk_xxxxx_xyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyz[i] * fz_0;

        tk_xxxxxy_xzz[i] = tk_xxxxx_xzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xzz[i] * fz_0;

        tk_xxxxxy_yyy[i] = 3.0 * tk_xxxxx_yy[i] * fe_0 + tk_xxxxx_yyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyy[i] * fz_0;

        tk_xxxxxy_yyz[i] = 2.0 * tk_xxxxx_yz[i] * fe_0 + tk_xxxxx_yyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyz[i] * fz_0;

        tk_xxxxxy_yzz[i] = tk_xxxxx_zz[i] * fe_0 + tk_xxxxx_yzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yzz[i] * fz_0;

        tk_xxxxxy_zzz[i] = tk_xxxxx_zzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : IF

    auto tk_xxxxxz_xxx = pbuffer.data(idx_kin_if + 20);

    auto tk_xxxxxz_xxy = pbuffer.data(idx_kin_if + 21);

    auto tk_xxxxxz_xxz = pbuffer.data(idx_kin_if + 22);

    auto tk_xxxxxz_xyy = pbuffer.data(idx_kin_if + 23);

    auto tk_xxxxxz_xyz = pbuffer.data(idx_kin_if + 24);

    auto tk_xxxxxz_xzz = pbuffer.data(idx_kin_if + 25);

    auto tk_xxxxxz_yyy = pbuffer.data(idx_kin_if + 26);

    auto tk_xxxxxz_yyz = pbuffer.data(idx_kin_if + 27);

    auto tk_xxxxxz_yzz = pbuffer.data(idx_kin_if + 28);

    auto tk_xxxxxz_zzz = pbuffer.data(idx_kin_if + 29);

#pragma omp simd aligned(pa_z,              \
                             tk_xxxxx_xx,   \
                             tk_xxxxx_xxx,  \
                             tk_xxxxx_xxy,  \
                             tk_xxxxx_xxz,  \
                             tk_xxxxx_xy,   \
                             tk_xxxxx_xyy,  \
                             tk_xxxxx_xyz,  \
                             tk_xxxxx_xz,   \
                             tk_xxxxx_xzz,  \
                             tk_xxxxx_yy,   \
                             tk_xxxxx_yyy,  \
                             tk_xxxxx_yyz,  \
                             tk_xxxxx_yz,   \
                             tk_xxxxx_yzz,  \
                             tk_xxxxx_zz,   \
                             tk_xxxxx_zzz,  \
                             tk_xxxxxz_xxx, \
                             tk_xxxxxz_xxy, \
                             tk_xxxxxz_xxz, \
                             tk_xxxxxz_xyy, \
                             tk_xxxxxz_xyz, \
                             tk_xxxxxz_xzz, \
                             tk_xxxxxz_yyy, \
                             tk_xxxxxz_yyz, \
                             tk_xxxxxz_yzz, \
                             tk_xxxxxz_zzz, \
                             ts_xxxxxz_xxx, \
                             ts_xxxxxz_xxy, \
                             ts_xxxxxz_xxz, \
                             ts_xxxxxz_xyy, \
                             ts_xxxxxz_xyz, \
                             ts_xxxxxz_xzz, \
                             ts_xxxxxz_yyy, \
                             ts_xxxxxz_yyz, \
                             ts_xxxxxz_yzz, \
                             ts_xxxxxz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_xxx[i] = tk_xxxxx_xxx[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxx[i] * fz_0;

        tk_xxxxxz_xxy[i] = tk_xxxxx_xxy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxy[i] * fz_0;

        tk_xxxxxz_xxz[i] = tk_xxxxx_xx[i] * fe_0 + tk_xxxxx_xxz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxz[i] * fz_0;

        tk_xxxxxz_xyy[i] = tk_xxxxx_xyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyy[i] * fz_0;

        tk_xxxxxz_xyz[i] = tk_xxxxx_xy[i] * fe_0 + tk_xxxxx_xyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyz[i] * fz_0;

        tk_xxxxxz_xzz[i] = 2.0 * tk_xxxxx_xz[i] * fe_0 + tk_xxxxx_xzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xzz[i] * fz_0;

        tk_xxxxxz_yyy[i] = tk_xxxxx_yyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyy[i] * fz_0;

        tk_xxxxxz_yyz[i] = tk_xxxxx_yy[i] * fe_0 + tk_xxxxx_yyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyz[i] * fz_0;

        tk_xxxxxz_yzz[i] = 2.0 * tk_xxxxx_yz[i] * fe_0 + tk_xxxxx_yzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yzz[i] * fz_0;

        tk_xxxxxz_zzz[i] = 3.0 * tk_xxxxx_zz[i] * fe_0 + tk_xxxxx_zzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_zzz[i] * fz_0;
    }

    // Set up 30-40 components of targeted buffer : IF

    auto tk_xxxxyy_xxx = pbuffer.data(idx_kin_if + 30);

    auto tk_xxxxyy_xxy = pbuffer.data(idx_kin_if + 31);

    auto tk_xxxxyy_xxz = pbuffer.data(idx_kin_if + 32);

    auto tk_xxxxyy_xyy = pbuffer.data(idx_kin_if + 33);

    auto tk_xxxxyy_xyz = pbuffer.data(idx_kin_if + 34);

    auto tk_xxxxyy_xzz = pbuffer.data(idx_kin_if + 35);

    auto tk_xxxxyy_yyy = pbuffer.data(idx_kin_if + 36);

    auto tk_xxxxyy_yyz = pbuffer.data(idx_kin_if + 37);

    auto tk_xxxxyy_yzz = pbuffer.data(idx_kin_if + 38);

    auto tk_xxxxyy_zzz = pbuffer.data(idx_kin_if + 39);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxz,   \
                             tk_xxxx_xzz,   \
                             tk_xxxxy_xxx,  \
                             tk_xxxxy_xxz,  \
                             tk_xxxxy_xzz,  \
                             tk_xxxxyy_xxx, \
                             tk_xxxxyy_xxy, \
                             tk_xxxxyy_xxz, \
                             tk_xxxxyy_xyy, \
                             tk_xxxxyy_xyz, \
                             tk_xxxxyy_xzz, \
                             tk_xxxxyy_yyy, \
                             tk_xxxxyy_yyz, \
                             tk_xxxxyy_yzz, \
                             tk_xxxxyy_zzz, \
                             tk_xxxyy_xxy,  \
                             tk_xxxyy_xy,   \
                             tk_xxxyy_xyy,  \
                             tk_xxxyy_xyz,  \
                             tk_xxxyy_yy,   \
                             tk_xxxyy_yyy,  \
                             tk_xxxyy_yyz,  \
                             tk_xxxyy_yz,   \
                             tk_xxxyy_yzz,  \
                             tk_xxxyy_zzz,  \
                             tk_xxyy_xxy,   \
                             tk_xxyy_xyy,   \
                             tk_xxyy_xyz,   \
                             tk_xxyy_yyy,   \
                             tk_xxyy_yyz,   \
                             tk_xxyy_yzz,   \
                             tk_xxyy_zzz,   \
                             ts_xxxx_xxx,   \
                             ts_xxxx_xxz,   \
                             ts_xxxx_xzz,   \
                             ts_xxxxyy_xxx, \
                             ts_xxxxyy_xxy, \
                             ts_xxxxyy_xxz, \
                             ts_xxxxyy_xyy, \
                             ts_xxxxyy_xyz, \
                             ts_xxxxyy_xzz, \
                             ts_xxxxyy_yyy, \
                             ts_xxxxyy_yyz, \
                             ts_xxxxyy_yzz, \
                             ts_xxxxyy_zzz, \
                             ts_xxyy_xxy,   \
                             ts_xxyy_xyy,   \
                             ts_xxyy_xyz,   \
                             ts_xxyy_yyy,   \
                             ts_xxyy_yyz,   \
                             ts_xxyy_yzz,   \
                             ts_xxyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_xxx[i] = -2.0 * ts_xxxx_xxx[i] * fbe_0 * fz_0 + tk_xxxx_xxx[i] * fe_0 + tk_xxxxy_xxx[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxx[i] * fz_0;

        tk_xxxxyy_xxy[i] = -6.0 * ts_xxyy_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxy[i] * fe_0 + 2.0 * tk_xxxyy_xy[i] * fe_0 +
                           tk_xxxyy_xxy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxy[i] * fz_0;

        tk_xxxxyy_xxz[i] = -2.0 * ts_xxxx_xxz[i] * fbe_0 * fz_0 + tk_xxxx_xxz[i] * fe_0 + tk_xxxxy_xxz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxz[i] * fz_0;

        tk_xxxxyy_xyy[i] = -6.0 * ts_xxyy_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyy[i] * fe_0 + tk_xxxyy_yy[i] * fe_0 + tk_xxxyy_xyy[i] * pa_x[i] +
                           2.0 * ts_xxxxyy_xyy[i] * fz_0;

        tk_xxxxyy_xyz[i] = -6.0 * ts_xxyy_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyz[i] * fe_0 + tk_xxxyy_yz[i] * fe_0 + tk_xxxyy_xyz[i] * pa_x[i] +
                           2.0 * ts_xxxxyy_xyz[i] * fz_0;

        tk_xxxxyy_xzz[i] = -2.0 * ts_xxxx_xzz[i] * fbe_0 * fz_0 + tk_xxxx_xzz[i] * fe_0 + tk_xxxxy_xzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xzz[i] * fz_0;

        tk_xxxxyy_yyy[i] =
            -6.0 * ts_xxyy_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyy[i] * fe_0 + tk_xxxyy_yyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyy[i] * fz_0;

        tk_xxxxyy_yyz[i] =
            -6.0 * ts_xxyy_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyz[i] * fe_0 + tk_xxxyy_yyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyz[i] * fz_0;

        tk_xxxxyy_yzz[i] =
            -6.0 * ts_xxyy_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yzz[i] * fe_0 + tk_xxxyy_yzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yzz[i] * fz_0;

        tk_xxxxyy_zzz[i] =
            -6.0 * ts_xxyy_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_zzz[i] * fe_0 + tk_xxxyy_zzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_zzz[i] * fz_0;
    }

    // Set up 40-50 components of targeted buffer : IF

    auto tk_xxxxyz_xxx = pbuffer.data(idx_kin_if + 40);

    auto tk_xxxxyz_xxy = pbuffer.data(idx_kin_if + 41);

    auto tk_xxxxyz_xxz = pbuffer.data(idx_kin_if + 42);

    auto tk_xxxxyz_xyy = pbuffer.data(idx_kin_if + 43);

    auto tk_xxxxyz_xyz = pbuffer.data(idx_kin_if + 44);

    auto tk_xxxxyz_xzz = pbuffer.data(idx_kin_if + 45);

    auto tk_xxxxyz_yyy = pbuffer.data(idx_kin_if + 46);

    auto tk_xxxxyz_yyz = pbuffer.data(idx_kin_if + 47);

    auto tk_xxxxyz_yzz = pbuffer.data(idx_kin_if + 48);

    auto tk_xxxxyz_zzz = pbuffer.data(idx_kin_if + 49);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_xxxxy_xxy,  \
                             tk_xxxxy_xyy,  \
                             tk_xxxxy_yyy,  \
                             tk_xxxxyz_xxx, \
                             tk_xxxxyz_xxy, \
                             tk_xxxxyz_xxz, \
                             tk_xxxxyz_xyy, \
                             tk_xxxxyz_xyz, \
                             tk_xxxxyz_xzz, \
                             tk_xxxxyz_yyy, \
                             tk_xxxxyz_yyz, \
                             tk_xxxxyz_yzz, \
                             tk_xxxxyz_zzz, \
                             tk_xxxxz_xxx,  \
                             tk_xxxxz_xxz,  \
                             tk_xxxxz_xyz,  \
                             tk_xxxxz_xz,   \
                             tk_xxxxz_xzz,  \
                             tk_xxxxz_yyz,  \
                             tk_xxxxz_yz,   \
                             tk_xxxxz_yzz,  \
                             tk_xxxxz_zz,   \
                             tk_xxxxz_zzz,  \
                             ts_xxxxyz_xxx, \
                             ts_xxxxyz_xxy, \
                             ts_xxxxyz_xxz, \
                             ts_xxxxyz_xyy, \
                             ts_xxxxyz_xyz, \
                             ts_xxxxyz_xzz, \
                             ts_xxxxyz_yyy, \
                             ts_xxxxyz_yyz, \
                             ts_xxxxyz_yzz, \
                             ts_xxxxyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxyz_xxx[i] = tk_xxxxz_xxx[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxx[i] * fz_0;

        tk_xxxxyz_xxy[i] = tk_xxxxy_xxy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxy[i] * fz_0;

        tk_xxxxyz_xxz[i] = tk_xxxxz_xxz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxz[i] * fz_0;

        tk_xxxxyz_xyy[i] = tk_xxxxy_xyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xyy[i] * fz_0;

        tk_xxxxyz_xyz[i] = tk_xxxxz_xz[i] * fe_0 + tk_xxxxz_xyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyz[i] * fz_0;

        tk_xxxxyz_xzz[i] = tk_xxxxz_xzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xzz[i] * fz_0;

        tk_xxxxyz_yyy[i] = tk_xxxxy_yyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_yyy[i] * fz_0;

        tk_xxxxyz_yyz[i] = 2.0 * tk_xxxxz_yz[i] * fe_0 + tk_xxxxz_yyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyz[i] * fz_0;

        tk_xxxxyz_yzz[i] = tk_xxxxz_zz[i] * fe_0 + tk_xxxxz_yzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yzz[i] * fz_0;

        tk_xxxxyz_zzz[i] = tk_xxxxz_zzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_zzz[i] * fz_0;
    }

    // Set up 50-60 components of targeted buffer : IF

    auto tk_xxxxzz_xxx = pbuffer.data(idx_kin_if + 50);

    auto tk_xxxxzz_xxy = pbuffer.data(idx_kin_if + 51);

    auto tk_xxxxzz_xxz = pbuffer.data(idx_kin_if + 52);

    auto tk_xxxxzz_xyy = pbuffer.data(idx_kin_if + 53);

    auto tk_xxxxzz_xyz = pbuffer.data(idx_kin_if + 54);

    auto tk_xxxxzz_xzz = pbuffer.data(idx_kin_if + 55);

    auto tk_xxxxzz_yyy = pbuffer.data(idx_kin_if + 56);

    auto tk_xxxxzz_yyz = pbuffer.data(idx_kin_if + 57);

    auto tk_xxxxzz_yzz = pbuffer.data(idx_kin_if + 58);

    auto tk_xxxxzz_zzz = pbuffer.data(idx_kin_if + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxy,   \
                             tk_xxxx_xyy,   \
                             tk_xxxxz_xxx,  \
                             tk_xxxxz_xxy,  \
                             tk_xxxxz_xyy,  \
                             tk_xxxxzz_xxx, \
                             tk_xxxxzz_xxy, \
                             tk_xxxxzz_xxz, \
                             tk_xxxxzz_xyy, \
                             tk_xxxxzz_xyz, \
                             tk_xxxxzz_xzz, \
                             tk_xxxxzz_yyy, \
                             tk_xxxxzz_yyz, \
                             tk_xxxxzz_yzz, \
                             tk_xxxxzz_zzz, \
                             tk_xxxzz_xxz,  \
                             tk_xxxzz_xyz,  \
                             tk_xxxzz_xz,   \
                             tk_xxxzz_xzz,  \
                             tk_xxxzz_yyy,  \
                             tk_xxxzz_yyz,  \
                             tk_xxxzz_yz,   \
                             tk_xxxzz_yzz,  \
                             tk_xxxzz_zz,   \
                             tk_xxxzz_zzz,  \
                             tk_xxzz_xxz,   \
                             tk_xxzz_xyz,   \
                             tk_xxzz_xzz,   \
                             tk_xxzz_yyy,   \
                             tk_xxzz_yyz,   \
                             tk_xxzz_yzz,   \
                             tk_xxzz_zzz,   \
                             ts_xxxx_xxx,   \
                             ts_xxxx_xxy,   \
                             ts_xxxx_xyy,   \
                             ts_xxxxzz_xxx, \
                             ts_xxxxzz_xxy, \
                             ts_xxxxzz_xxz, \
                             ts_xxxxzz_xyy, \
                             ts_xxxxzz_xyz, \
                             ts_xxxxzz_xzz, \
                             ts_xxxxzz_yyy, \
                             ts_xxxxzz_yyz, \
                             ts_xxxxzz_yzz, \
                             ts_xxxxzz_zzz, \
                             ts_xxzz_xxz,   \
                             ts_xxzz_xyz,   \
                             ts_xxzz_xzz,   \
                             ts_xxzz_yyy,   \
                             ts_xxzz_yyz,   \
                             ts_xxzz_yzz,   \
                             ts_xxzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_xxx[i] = -2.0 * ts_xxxx_xxx[i] * fbe_0 * fz_0 + tk_xxxx_xxx[i] * fe_0 + tk_xxxxz_xxx[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxx[i] * fz_0;

        tk_xxxxzz_xxy[i] = -2.0 * ts_xxxx_xxy[i] * fbe_0 * fz_0 + tk_xxxx_xxy[i] * fe_0 + tk_xxxxz_xxy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxy[i] * fz_0;

        tk_xxxxzz_xxz[i] = -6.0 * ts_xxzz_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxz[i] * fe_0 + 2.0 * tk_xxxzz_xz[i] * fe_0 +
                           tk_xxxzz_xxz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxz[i] * fz_0;

        tk_xxxxzz_xyy[i] = -2.0 * ts_xxxx_xyy[i] * fbe_0 * fz_0 + tk_xxxx_xyy[i] * fe_0 + tk_xxxxz_xyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xyy[i] * fz_0;

        tk_xxxxzz_xyz[i] = -6.0 * ts_xxzz_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyz[i] * fe_0 + tk_xxxzz_yz[i] * fe_0 + tk_xxxzz_xyz[i] * pa_x[i] +
                           2.0 * ts_xxxxzz_xyz[i] * fz_0;

        tk_xxxxzz_xzz[i] = -6.0 * ts_xxzz_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xzz[i] * fe_0 + tk_xxxzz_zz[i] * fe_0 + tk_xxxzz_xzz[i] * pa_x[i] +
                           2.0 * ts_xxxxzz_xzz[i] * fz_0;

        tk_xxxxzz_yyy[i] =
            -6.0 * ts_xxzz_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyy[i] * fe_0 + tk_xxxzz_yyy[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyy[i] * fz_0;

        tk_xxxxzz_yyz[i] =
            -6.0 * ts_xxzz_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyz[i] * fe_0 + tk_xxxzz_yyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyz[i] * fz_0;

        tk_xxxxzz_yzz[i] =
            -6.0 * ts_xxzz_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yzz[i] * fe_0 + tk_xxxzz_yzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yzz[i] * fz_0;

        tk_xxxxzz_zzz[i] =
            -6.0 * ts_xxzz_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_zzz[i] * fe_0 + tk_xxxzz_zzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_zzz[i] * fz_0;
    }

    // Set up 60-70 components of targeted buffer : IF

    auto tk_xxxyyy_xxx = pbuffer.data(idx_kin_if + 60);

    auto tk_xxxyyy_xxy = pbuffer.data(idx_kin_if + 61);

    auto tk_xxxyyy_xxz = pbuffer.data(idx_kin_if + 62);

    auto tk_xxxyyy_xyy = pbuffer.data(idx_kin_if + 63);

    auto tk_xxxyyy_xyz = pbuffer.data(idx_kin_if + 64);

    auto tk_xxxyyy_xzz = pbuffer.data(idx_kin_if + 65);

    auto tk_xxxyyy_yyy = pbuffer.data(idx_kin_if + 66);

    auto tk_xxxyyy_yyz = pbuffer.data(idx_kin_if + 67);

    auto tk_xxxyyy_yzz = pbuffer.data(idx_kin_if + 68);

    auto tk_xxxyyy_zzz = pbuffer.data(idx_kin_if + 69);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xxxy_xxx,   \
                             tk_xxxy_xxz,   \
                             tk_xxxy_xzz,   \
                             tk_xxxyy_xxx,  \
                             tk_xxxyy_xxz,  \
                             tk_xxxyy_xzz,  \
                             tk_xxxyyy_xxx, \
                             tk_xxxyyy_xxy, \
                             tk_xxxyyy_xxz, \
                             tk_xxxyyy_xyy, \
                             tk_xxxyyy_xyz, \
                             tk_xxxyyy_xzz, \
                             tk_xxxyyy_yyy, \
                             tk_xxxyyy_yyz, \
                             tk_xxxyyy_yzz, \
                             tk_xxxyyy_zzz, \
                             tk_xxyyy_xxy,  \
                             tk_xxyyy_xy,   \
                             tk_xxyyy_xyy,  \
                             tk_xxyyy_xyz,  \
                             tk_xxyyy_yy,   \
                             tk_xxyyy_yyy,  \
                             tk_xxyyy_yyz,  \
                             tk_xxyyy_yz,   \
                             tk_xxyyy_yzz,  \
                             tk_xxyyy_zzz,  \
                             tk_xyyy_xxy,   \
                             tk_xyyy_xyy,   \
                             tk_xyyy_xyz,   \
                             tk_xyyy_yyy,   \
                             tk_xyyy_yyz,   \
                             tk_xyyy_yzz,   \
                             tk_xyyy_zzz,   \
                             ts_xxxy_xxx,   \
                             ts_xxxy_xxz,   \
                             ts_xxxy_xzz,   \
                             ts_xxxyyy_xxx, \
                             ts_xxxyyy_xxy, \
                             ts_xxxyyy_xxz, \
                             ts_xxxyyy_xyy, \
                             ts_xxxyyy_xyz, \
                             ts_xxxyyy_xzz, \
                             ts_xxxyyy_yyy, \
                             ts_xxxyyy_yyz, \
                             ts_xxxyyy_yzz, \
                             ts_xxxyyy_zzz, \
                             ts_xyyy_xxy,   \
                             ts_xyyy_xyy,   \
                             ts_xyyy_xyz,   \
                             ts_xyyy_yyy,   \
                             ts_xyyy_yyz,   \
                             ts_xyyy_yzz,   \
                             ts_xyyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_xxx[i] =
            -4.0 * ts_xxxy_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxx[i] * fe_0 + tk_xxxyy_xxx[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxx[i] * fz_0;

        tk_xxxyyy_xxy[i] = -4.0 * ts_xyyy_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxy[i] * fe_0 + 2.0 * tk_xxyyy_xy[i] * fe_0 +
                           tk_xxyyy_xxy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxy[i] * fz_0;

        tk_xxxyyy_xxz[i] =
            -4.0 * ts_xxxy_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxz[i] * fe_0 + tk_xxxyy_xxz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxz[i] * fz_0;

        tk_xxxyyy_xyy[i] = -4.0 * ts_xyyy_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyy[i] * fe_0 + tk_xxyyy_yy[i] * fe_0 + tk_xxyyy_xyy[i] * pa_x[i] +
                           2.0 * ts_xxxyyy_xyy[i] * fz_0;

        tk_xxxyyy_xyz[i] = -4.0 * ts_xyyy_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyz[i] * fe_0 + tk_xxyyy_yz[i] * fe_0 + tk_xxyyy_xyz[i] * pa_x[i] +
                           2.0 * ts_xxxyyy_xyz[i] * fz_0;

        tk_xxxyyy_xzz[i] =
            -4.0 * ts_xxxy_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xzz[i] * fe_0 + tk_xxxyy_xzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xzz[i] * fz_0;

        tk_xxxyyy_yyy[i] =
            -4.0 * ts_xyyy_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyy[i] * fe_0 + tk_xxyyy_yyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyy[i] * fz_0;

        tk_xxxyyy_yyz[i] =
            -4.0 * ts_xyyy_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyz[i] * fe_0 + tk_xxyyy_yyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyz[i] * fz_0;

        tk_xxxyyy_yzz[i] =
            -4.0 * ts_xyyy_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yzz[i] * fe_0 + tk_xxyyy_yzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yzz[i] * fz_0;

        tk_xxxyyy_zzz[i] =
            -4.0 * ts_xyyy_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_zzz[i] * fe_0 + tk_xxyyy_zzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_zzz[i] * fz_0;
    }

    // Set up 70-80 components of targeted buffer : IF

    auto tk_xxxyyz_xxx = pbuffer.data(idx_kin_if + 70);

    auto tk_xxxyyz_xxy = pbuffer.data(idx_kin_if + 71);

    auto tk_xxxyyz_xxz = pbuffer.data(idx_kin_if + 72);

    auto tk_xxxyyz_xyy = pbuffer.data(idx_kin_if + 73);

    auto tk_xxxyyz_xyz = pbuffer.data(idx_kin_if + 74);

    auto tk_xxxyyz_xzz = pbuffer.data(idx_kin_if + 75);

    auto tk_xxxyyz_yyy = pbuffer.data(idx_kin_if + 76);

    auto tk_xxxyyz_yyz = pbuffer.data(idx_kin_if + 77);

    auto tk_xxxyyz_yzz = pbuffer.data(idx_kin_if + 78);

    auto tk_xxxyyz_zzz = pbuffer.data(idx_kin_if + 79);

#pragma omp simd aligned(pa_z,              \
                             tk_xxxyy_xx,   \
                             tk_xxxyy_xxx,  \
                             tk_xxxyy_xxy,  \
                             tk_xxxyy_xxz,  \
                             tk_xxxyy_xy,   \
                             tk_xxxyy_xyy,  \
                             tk_xxxyy_xyz,  \
                             tk_xxxyy_xz,   \
                             tk_xxxyy_xzz,  \
                             tk_xxxyy_yy,   \
                             tk_xxxyy_yyy,  \
                             tk_xxxyy_yyz,  \
                             tk_xxxyy_yz,   \
                             tk_xxxyy_yzz,  \
                             tk_xxxyy_zz,   \
                             tk_xxxyy_zzz,  \
                             tk_xxxyyz_xxx, \
                             tk_xxxyyz_xxy, \
                             tk_xxxyyz_xxz, \
                             tk_xxxyyz_xyy, \
                             tk_xxxyyz_xyz, \
                             tk_xxxyyz_xzz, \
                             tk_xxxyyz_yyy, \
                             tk_xxxyyz_yyz, \
                             tk_xxxyyz_yzz, \
                             tk_xxxyyz_zzz, \
                             ts_xxxyyz_xxx, \
                             ts_xxxyyz_xxy, \
                             ts_xxxyyz_xxz, \
                             ts_xxxyyz_xyy, \
                             ts_xxxyyz_xyz, \
                             ts_xxxyyz_xzz, \
                             ts_xxxyyz_yyy, \
                             ts_xxxyyz_yyz, \
                             ts_xxxyyz_yzz, \
                             ts_xxxyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_xxx[i] = tk_xxxyy_xxx[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxx[i] * fz_0;

        tk_xxxyyz_xxy[i] = tk_xxxyy_xxy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxy[i] * fz_0;

        tk_xxxyyz_xxz[i] = tk_xxxyy_xx[i] * fe_0 + tk_xxxyy_xxz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxz[i] * fz_0;

        tk_xxxyyz_xyy[i] = tk_xxxyy_xyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyy[i] * fz_0;

        tk_xxxyyz_xyz[i] = tk_xxxyy_xy[i] * fe_0 + tk_xxxyy_xyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyz[i] * fz_0;

        tk_xxxyyz_xzz[i] = 2.0 * tk_xxxyy_xz[i] * fe_0 + tk_xxxyy_xzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xzz[i] * fz_0;

        tk_xxxyyz_yyy[i] = tk_xxxyy_yyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyy[i] * fz_0;

        tk_xxxyyz_yyz[i] = tk_xxxyy_yy[i] * fe_0 + tk_xxxyy_yyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyz[i] * fz_0;

        tk_xxxyyz_yzz[i] = 2.0 * tk_xxxyy_yz[i] * fe_0 + tk_xxxyy_yzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yzz[i] * fz_0;

        tk_xxxyyz_zzz[i] = 3.0 * tk_xxxyy_zz[i] * fe_0 + tk_xxxyy_zzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_zzz[i] * fz_0;
    }

    // Set up 80-90 components of targeted buffer : IF

    auto tk_xxxyzz_xxx = pbuffer.data(idx_kin_if + 80);

    auto tk_xxxyzz_xxy = pbuffer.data(idx_kin_if + 81);

    auto tk_xxxyzz_xxz = pbuffer.data(idx_kin_if + 82);

    auto tk_xxxyzz_xyy = pbuffer.data(idx_kin_if + 83);

    auto tk_xxxyzz_xyz = pbuffer.data(idx_kin_if + 84);

    auto tk_xxxyzz_xzz = pbuffer.data(idx_kin_if + 85);

    auto tk_xxxyzz_yyy = pbuffer.data(idx_kin_if + 86);

    auto tk_xxxyzz_yyz = pbuffer.data(idx_kin_if + 87);

    auto tk_xxxyzz_yzz = pbuffer.data(idx_kin_if + 88);

    auto tk_xxxyzz_zzz = pbuffer.data(idx_kin_if + 89);

#pragma omp simd aligned(pa_y,              \
                             tk_xxxyzz_xxx, \
                             tk_xxxyzz_xxy, \
                             tk_xxxyzz_xxz, \
                             tk_xxxyzz_xyy, \
                             tk_xxxyzz_xyz, \
                             tk_xxxyzz_xzz, \
                             tk_xxxyzz_yyy, \
                             tk_xxxyzz_yyz, \
                             tk_xxxyzz_yzz, \
                             tk_xxxyzz_zzz, \
                             tk_xxxzz_xx,   \
                             tk_xxxzz_xxx,  \
                             tk_xxxzz_xxy,  \
                             tk_xxxzz_xxz,  \
                             tk_xxxzz_xy,   \
                             tk_xxxzz_xyy,  \
                             tk_xxxzz_xyz,  \
                             tk_xxxzz_xz,   \
                             tk_xxxzz_xzz,  \
                             tk_xxxzz_yy,   \
                             tk_xxxzz_yyy,  \
                             tk_xxxzz_yyz,  \
                             tk_xxxzz_yz,   \
                             tk_xxxzz_yzz,  \
                             tk_xxxzz_zz,   \
                             tk_xxxzz_zzz,  \
                             ts_xxxyzz_xxx, \
                             ts_xxxyzz_xxy, \
                             ts_xxxyzz_xxz, \
                             ts_xxxyzz_xyy, \
                             ts_xxxyzz_xyz, \
                             ts_xxxyzz_xzz, \
                             ts_xxxyzz_yyy, \
                             ts_xxxyzz_yyz, \
                             ts_xxxyzz_yzz, \
                             ts_xxxyzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_xxx[i] = tk_xxxzz_xxx[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxx[i] * fz_0;

        tk_xxxyzz_xxy[i] = tk_xxxzz_xx[i] * fe_0 + tk_xxxzz_xxy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxy[i] * fz_0;

        tk_xxxyzz_xxz[i] = tk_xxxzz_xxz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxz[i] * fz_0;

        tk_xxxyzz_xyy[i] = 2.0 * tk_xxxzz_xy[i] * fe_0 + tk_xxxzz_xyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyy[i] * fz_0;

        tk_xxxyzz_xyz[i] = tk_xxxzz_xz[i] * fe_0 + tk_xxxzz_xyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyz[i] * fz_0;

        tk_xxxyzz_xzz[i] = tk_xxxzz_xzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xzz[i] * fz_0;

        tk_xxxyzz_yyy[i] = 3.0 * tk_xxxzz_yy[i] * fe_0 + tk_xxxzz_yyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyy[i] * fz_0;

        tk_xxxyzz_yyz[i] = 2.0 * tk_xxxzz_yz[i] * fe_0 + tk_xxxzz_yyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyz[i] * fz_0;

        tk_xxxyzz_yzz[i] = tk_xxxzz_zz[i] * fe_0 + tk_xxxzz_yzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yzz[i] * fz_0;

        tk_xxxyzz_zzz[i] = tk_xxxzz_zzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_zzz[i] * fz_0;
    }

    // Set up 90-100 components of targeted buffer : IF

    auto tk_xxxzzz_xxx = pbuffer.data(idx_kin_if + 90);

    auto tk_xxxzzz_xxy = pbuffer.data(idx_kin_if + 91);

    auto tk_xxxzzz_xxz = pbuffer.data(idx_kin_if + 92);

    auto tk_xxxzzz_xyy = pbuffer.data(idx_kin_if + 93);

    auto tk_xxxzzz_xyz = pbuffer.data(idx_kin_if + 94);

    auto tk_xxxzzz_xzz = pbuffer.data(idx_kin_if + 95);

    auto tk_xxxzzz_yyy = pbuffer.data(idx_kin_if + 96);

    auto tk_xxxzzz_yyz = pbuffer.data(idx_kin_if + 97);

    auto tk_xxxzzz_yzz = pbuffer.data(idx_kin_if + 98);

    auto tk_xxxzzz_zzz = pbuffer.data(idx_kin_if + 99);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xxxz_xxx,   \
                             tk_xxxz_xxy,   \
                             tk_xxxz_xyy,   \
                             tk_xxxzz_xxx,  \
                             tk_xxxzz_xxy,  \
                             tk_xxxzz_xyy,  \
                             tk_xxxzzz_xxx, \
                             tk_xxxzzz_xxy, \
                             tk_xxxzzz_xxz, \
                             tk_xxxzzz_xyy, \
                             tk_xxxzzz_xyz, \
                             tk_xxxzzz_xzz, \
                             tk_xxxzzz_yyy, \
                             tk_xxxzzz_yyz, \
                             tk_xxxzzz_yzz, \
                             tk_xxxzzz_zzz, \
                             tk_xxzzz_xxz,  \
                             tk_xxzzz_xyz,  \
                             tk_xxzzz_xz,   \
                             tk_xxzzz_xzz,  \
                             tk_xxzzz_yyy,  \
                             tk_xxzzz_yyz,  \
                             tk_xxzzz_yz,   \
                             tk_xxzzz_yzz,  \
                             tk_xxzzz_zz,   \
                             tk_xxzzz_zzz,  \
                             tk_xzzz_xxz,   \
                             tk_xzzz_xyz,   \
                             tk_xzzz_xzz,   \
                             tk_xzzz_yyy,   \
                             tk_xzzz_yyz,   \
                             tk_xzzz_yzz,   \
                             tk_xzzz_zzz,   \
                             ts_xxxz_xxx,   \
                             ts_xxxz_xxy,   \
                             ts_xxxz_xyy,   \
                             ts_xxxzzz_xxx, \
                             ts_xxxzzz_xxy, \
                             ts_xxxzzz_xxz, \
                             ts_xxxzzz_xyy, \
                             ts_xxxzzz_xyz, \
                             ts_xxxzzz_xzz, \
                             ts_xxxzzz_yyy, \
                             ts_xxxzzz_yyz, \
                             ts_xxxzzz_yzz, \
                             ts_xxxzzz_zzz, \
                             ts_xzzz_xxz,   \
                             ts_xzzz_xyz,   \
                             ts_xzzz_xzz,   \
                             ts_xzzz_yyy,   \
                             ts_xzzz_yyz,   \
                             ts_xzzz_yzz,   \
                             ts_xzzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_xxx[i] =
            -4.0 * ts_xxxz_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxx[i] * fe_0 + tk_xxxzz_xxx[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxx[i] * fz_0;

        tk_xxxzzz_xxy[i] =
            -4.0 * ts_xxxz_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxy[i] * fe_0 + tk_xxxzz_xxy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxy[i] * fz_0;

        tk_xxxzzz_xxz[i] = -4.0 * ts_xzzz_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxz[i] * fe_0 + 2.0 * tk_xxzzz_xz[i] * fe_0 +
                           tk_xxzzz_xxz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxz[i] * fz_0;

        tk_xxxzzz_xyy[i] =
            -4.0 * ts_xxxz_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xyy[i] * fe_0 + tk_xxxzz_xyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xyy[i] * fz_0;

        tk_xxxzzz_xyz[i] = -4.0 * ts_xzzz_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyz[i] * fe_0 + tk_xxzzz_yz[i] * fe_0 + tk_xxzzz_xyz[i] * pa_x[i] +
                           2.0 * ts_xxxzzz_xyz[i] * fz_0;

        tk_xxxzzz_xzz[i] = -4.0 * ts_xzzz_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xzz[i] * fe_0 + tk_xxzzz_zz[i] * fe_0 + tk_xxzzz_xzz[i] * pa_x[i] +
                           2.0 * ts_xxxzzz_xzz[i] * fz_0;

        tk_xxxzzz_yyy[i] =
            -4.0 * ts_xzzz_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyy[i] * fe_0 + tk_xxzzz_yyy[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyy[i] * fz_0;

        tk_xxxzzz_yyz[i] =
            -4.0 * ts_xzzz_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyz[i] * fe_0 + tk_xxzzz_yyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyz[i] * fz_0;

        tk_xxxzzz_yzz[i] =
            -4.0 * ts_xzzz_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yzz[i] * fe_0 + tk_xxzzz_yzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yzz[i] * fz_0;

        tk_xxxzzz_zzz[i] =
            -4.0 * ts_xzzz_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_zzz[i] * fe_0 + tk_xxzzz_zzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_zzz[i] * fz_0;
    }

    // Set up 100-110 components of targeted buffer : IF

    auto tk_xxyyyy_xxx = pbuffer.data(idx_kin_if + 100);

    auto tk_xxyyyy_xxy = pbuffer.data(idx_kin_if + 101);

    auto tk_xxyyyy_xxz = pbuffer.data(idx_kin_if + 102);

    auto tk_xxyyyy_xyy = pbuffer.data(idx_kin_if + 103);

    auto tk_xxyyyy_xyz = pbuffer.data(idx_kin_if + 104);

    auto tk_xxyyyy_xzz = pbuffer.data(idx_kin_if + 105);

    auto tk_xxyyyy_yyy = pbuffer.data(idx_kin_if + 106);

    auto tk_xxyyyy_yyz = pbuffer.data(idx_kin_if + 107);

    auto tk_xxyyyy_yzz = pbuffer.data(idx_kin_if + 108);

    auto tk_xxyyyy_zzz = pbuffer.data(idx_kin_if + 109);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xxyy_xxx,   \
                             tk_xxyy_xxz,   \
                             tk_xxyy_xzz,   \
                             tk_xxyyy_xxx,  \
                             tk_xxyyy_xxz,  \
                             tk_xxyyy_xzz,  \
                             tk_xxyyyy_xxx, \
                             tk_xxyyyy_xxy, \
                             tk_xxyyyy_xxz, \
                             tk_xxyyyy_xyy, \
                             tk_xxyyyy_xyz, \
                             tk_xxyyyy_xzz, \
                             tk_xxyyyy_yyy, \
                             tk_xxyyyy_yyz, \
                             tk_xxyyyy_yzz, \
                             tk_xxyyyy_zzz, \
                             tk_xyyyy_xxy,  \
                             tk_xyyyy_xy,   \
                             tk_xyyyy_xyy,  \
                             tk_xyyyy_xyz,  \
                             tk_xyyyy_yy,   \
                             tk_xyyyy_yyy,  \
                             tk_xyyyy_yyz,  \
                             tk_xyyyy_yz,   \
                             tk_xyyyy_yzz,  \
                             tk_xyyyy_zzz,  \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xyy,   \
                             tk_yyyy_xyz,   \
                             tk_yyyy_yyy,   \
                             tk_yyyy_yyz,   \
                             tk_yyyy_yzz,   \
                             tk_yyyy_zzz,   \
                             ts_xxyy_xxx,   \
                             ts_xxyy_xxz,   \
                             ts_xxyy_xzz,   \
                             ts_xxyyyy_xxx, \
                             ts_xxyyyy_xxy, \
                             ts_xxyyyy_xxz, \
                             ts_xxyyyy_xyy, \
                             ts_xxyyyy_xyz, \
                             ts_xxyyyy_xzz, \
                             ts_xxyyyy_yyy, \
                             ts_xxyyyy_yyz, \
                             ts_xxyyyy_yzz, \
                             ts_xxyyyy_zzz, \
                             ts_yyyy_xxy,   \
                             ts_yyyy_xyy,   \
                             ts_yyyy_xyz,   \
                             ts_yyyy_yyy,   \
                             ts_yyyy_yyz,   \
                             ts_yyyy_yzz,   \
                             ts_yyyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_xxx[i] =
            -6.0 * ts_xxyy_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxx[i] * fe_0 + tk_xxyyy_xxx[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxx[i] * fz_0;

        tk_xxyyyy_xxy[i] = -2.0 * ts_yyyy_xxy[i] * fbe_0 * fz_0 + tk_yyyy_xxy[i] * fe_0 + 2.0 * tk_xyyyy_xy[i] * fe_0 + tk_xyyyy_xxy[i] * pa_x[i] +
                           2.0 * ts_xxyyyy_xxy[i] * fz_0;

        tk_xxyyyy_xxz[i] =
            -6.0 * ts_xxyy_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxz[i] * fe_0 + tk_xxyyy_xxz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxz[i] * fz_0;

        tk_xxyyyy_xyy[i] = -2.0 * ts_yyyy_xyy[i] * fbe_0 * fz_0 + tk_yyyy_xyy[i] * fe_0 + tk_xyyyy_yy[i] * fe_0 + tk_xyyyy_xyy[i] * pa_x[i] +
                           2.0 * ts_xxyyyy_xyy[i] * fz_0;

        tk_xxyyyy_xyz[i] = -2.0 * ts_yyyy_xyz[i] * fbe_0 * fz_0 + tk_yyyy_xyz[i] * fe_0 + tk_xyyyy_yz[i] * fe_0 + tk_xyyyy_xyz[i] * pa_x[i] +
                           2.0 * ts_xxyyyy_xyz[i] * fz_0;

        tk_xxyyyy_xzz[i] =
            -6.0 * ts_xxyy_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xzz[i] * fe_0 + tk_xxyyy_xzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xzz[i] * fz_0;

        tk_xxyyyy_yyy[i] = -2.0 * ts_yyyy_yyy[i] * fbe_0 * fz_0 + tk_yyyy_yyy[i] * fe_0 + tk_xyyyy_yyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyy[i] * fz_0;

        tk_xxyyyy_yyz[i] = -2.0 * ts_yyyy_yyz[i] * fbe_0 * fz_0 + tk_yyyy_yyz[i] * fe_0 + tk_xyyyy_yyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyz[i] * fz_0;

        tk_xxyyyy_yzz[i] = -2.0 * ts_yyyy_yzz[i] * fbe_0 * fz_0 + tk_yyyy_yzz[i] * fe_0 + tk_xyyyy_yzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yzz[i] * fz_0;

        tk_xxyyyy_zzz[i] = -2.0 * ts_yyyy_zzz[i] * fbe_0 * fz_0 + tk_yyyy_zzz[i] * fe_0 + tk_xyyyy_zzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_zzz[i] * fz_0;
    }

    // Set up 110-120 components of targeted buffer : IF

    auto tk_xxyyyz_xxx = pbuffer.data(idx_kin_if + 110);

    auto tk_xxyyyz_xxy = pbuffer.data(idx_kin_if + 111);

    auto tk_xxyyyz_xxz = pbuffer.data(idx_kin_if + 112);

    auto tk_xxyyyz_xyy = pbuffer.data(idx_kin_if + 113);

    auto tk_xxyyyz_xyz = pbuffer.data(idx_kin_if + 114);

    auto tk_xxyyyz_xzz = pbuffer.data(idx_kin_if + 115);

    auto tk_xxyyyz_yyy = pbuffer.data(idx_kin_if + 116);

    auto tk_xxyyyz_yyz = pbuffer.data(idx_kin_if + 117);

    auto tk_xxyyyz_yzz = pbuffer.data(idx_kin_if + 118);

    auto tk_xxyyyz_zzz = pbuffer.data(idx_kin_if + 119);

#pragma omp simd aligned(pa_z,              \
                             tk_xxyyy_xx,   \
                             tk_xxyyy_xxx,  \
                             tk_xxyyy_xxy,  \
                             tk_xxyyy_xxz,  \
                             tk_xxyyy_xy,   \
                             tk_xxyyy_xyy,  \
                             tk_xxyyy_xyz,  \
                             tk_xxyyy_xz,   \
                             tk_xxyyy_xzz,  \
                             tk_xxyyy_yy,   \
                             tk_xxyyy_yyy,  \
                             tk_xxyyy_yyz,  \
                             tk_xxyyy_yz,   \
                             tk_xxyyy_yzz,  \
                             tk_xxyyy_zz,   \
                             tk_xxyyy_zzz,  \
                             tk_xxyyyz_xxx, \
                             tk_xxyyyz_xxy, \
                             tk_xxyyyz_xxz, \
                             tk_xxyyyz_xyy, \
                             tk_xxyyyz_xyz, \
                             tk_xxyyyz_xzz, \
                             tk_xxyyyz_yyy, \
                             tk_xxyyyz_yyz, \
                             tk_xxyyyz_yzz, \
                             tk_xxyyyz_zzz, \
                             ts_xxyyyz_xxx, \
                             ts_xxyyyz_xxy, \
                             ts_xxyyyz_xxz, \
                             ts_xxyyyz_xyy, \
                             ts_xxyyyz_xyz, \
                             ts_xxyyyz_xzz, \
                             ts_xxyyyz_yyy, \
                             ts_xxyyyz_yyz, \
                             ts_xxyyyz_yzz, \
                             ts_xxyyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_xxx[i] = tk_xxyyy_xxx[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxx[i] * fz_0;

        tk_xxyyyz_xxy[i] = tk_xxyyy_xxy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxy[i] * fz_0;

        tk_xxyyyz_xxz[i] = tk_xxyyy_xx[i] * fe_0 + tk_xxyyy_xxz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxz[i] * fz_0;

        tk_xxyyyz_xyy[i] = tk_xxyyy_xyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyy[i] * fz_0;

        tk_xxyyyz_xyz[i] = tk_xxyyy_xy[i] * fe_0 + tk_xxyyy_xyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyz[i] * fz_0;

        tk_xxyyyz_xzz[i] = 2.0 * tk_xxyyy_xz[i] * fe_0 + tk_xxyyy_xzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xzz[i] * fz_0;

        tk_xxyyyz_yyy[i] = tk_xxyyy_yyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyy[i] * fz_0;

        tk_xxyyyz_yyz[i] = tk_xxyyy_yy[i] * fe_0 + tk_xxyyy_yyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyz[i] * fz_0;

        tk_xxyyyz_yzz[i] = 2.0 * tk_xxyyy_yz[i] * fe_0 + tk_xxyyy_yzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yzz[i] * fz_0;

        tk_xxyyyz_zzz[i] = 3.0 * tk_xxyyy_zz[i] * fe_0 + tk_xxyyy_zzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_zzz[i] * fz_0;
    }

    // Set up 120-130 components of targeted buffer : IF

    auto tk_xxyyzz_xxx = pbuffer.data(idx_kin_if + 120);

    auto tk_xxyyzz_xxy = pbuffer.data(idx_kin_if + 121);

    auto tk_xxyyzz_xxz = pbuffer.data(idx_kin_if + 122);

    auto tk_xxyyzz_xyy = pbuffer.data(idx_kin_if + 123);

    auto tk_xxyyzz_xyz = pbuffer.data(idx_kin_if + 124);

    auto tk_xxyyzz_xzz = pbuffer.data(idx_kin_if + 125);

    auto tk_xxyyzz_yyy = pbuffer.data(idx_kin_if + 126);

    auto tk_xxyyzz_yyz = pbuffer.data(idx_kin_if + 127);

    auto tk_xxyyzz_yzz = pbuffer.data(idx_kin_if + 128);

    auto tk_xxyyzz_zzz = pbuffer.data(idx_kin_if + 129);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             tk_xxyy_xxy,   \
                             tk_xxyy_xyy,   \
                             tk_xxyyz_xxy,  \
                             tk_xxyyz_xyy,  \
                             tk_xxyyzz_xxx, \
                             tk_xxyyzz_xxy, \
                             tk_xxyyzz_xxz, \
                             tk_xxyyzz_xyy, \
                             tk_xxyyzz_xyz, \
                             tk_xxyyzz_xzz, \
                             tk_xxyyzz_yyy, \
                             tk_xxyyzz_yyz, \
                             tk_xxyyzz_yzz, \
                             tk_xxyyzz_zzz, \
                             tk_xxyzz_xxx,  \
                             tk_xxyzz_xxz,  \
                             tk_xxyzz_xzz,  \
                             tk_xxzz_xxx,   \
                             tk_xxzz_xxz,   \
                             tk_xxzz_xzz,   \
                             tk_xyyzz_xyz,  \
                             tk_xyyzz_yyy,  \
                             tk_xyyzz_yyz,  \
                             tk_xyyzz_yz,   \
                             tk_xyyzz_yzz,  \
                             tk_xyyzz_zzz,  \
                             tk_yyzz_xyz,   \
                             tk_yyzz_yyy,   \
                             tk_yyzz_yyz,   \
                             tk_yyzz_yzz,   \
                             tk_yyzz_zzz,   \
                             ts_xxyy_xxy,   \
                             ts_xxyy_xyy,   \
                             ts_xxyyzz_xxx, \
                             ts_xxyyzz_xxy, \
                             ts_xxyyzz_xxz, \
                             ts_xxyyzz_xyy, \
                             ts_xxyyzz_xyz, \
                             ts_xxyyzz_xzz, \
                             ts_xxyyzz_yyy, \
                             ts_xxyyzz_yyz, \
                             ts_xxyyzz_yzz, \
                             ts_xxyyzz_zzz, \
                             ts_xxzz_xxx,   \
                             ts_xxzz_xxz,   \
                             ts_xxzz_xzz,   \
                             ts_yyzz_xyz,   \
                             ts_yyzz_yyy,   \
                             ts_yyzz_yyz,   \
                             ts_yyzz_yzz,   \
                             ts_yyzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_xxx[i] = -2.0 * ts_xxzz_xxx[i] * fbe_0 * fz_0 + tk_xxzz_xxx[i] * fe_0 + tk_xxyzz_xxx[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxx[i] * fz_0;

        tk_xxyyzz_xxy[i] = -2.0 * ts_xxyy_xxy[i] * fbe_0 * fz_0 + tk_xxyy_xxy[i] * fe_0 + tk_xxyyz_xxy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxy[i] * fz_0;

        tk_xxyyzz_xxz[i] = -2.0 * ts_xxzz_xxz[i] * fbe_0 * fz_0 + tk_xxzz_xxz[i] * fe_0 + tk_xxyzz_xxz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxz[i] * fz_0;

        tk_xxyyzz_xyy[i] = -2.0 * ts_xxyy_xyy[i] * fbe_0 * fz_0 + tk_xxyy_xyy[i] * fe_0 + tk_xxyyz_xyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xyy[i] * fz_0;

        tk_xxyyzz_xyz[i] = -2.0 * ts_yyzz_xyz[i] * fbe_0 * fz_0 + tk_yyzz_xyz[i] * fe_0 + tk_xyyzz_yz[i] * fe_0 + tk_xyyzz_xyz[i] * pa_x[i] +
                           2.0 * ts_xxyyzz_xyz[i] * fz_0;

        tk_xxyyzz_xzz[i] = -2.0 * ts_xxzz_xzz[i] * fbe_0 * fz_0 + tk_xxzz_xzz[i] * fe_0 + tk_xxyzz_xzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xzz[i] * fz_0;

        tk_xxyyzz_yyy[i] = -2.0 * ts_yyzz_yyy[i] * fbe_0 * fz_0 + tk_yyzz_yyy[i] * fe_0 + tk_xyyzz_yyy[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyy[i] * fz_0;

        tk_xxyyzz_yyz[i] = -2.0 * ts_yyzz_yyz[i] * fbe_0 * fz_0 + tk_yyzz_yyz[i] * fe_0 + tk_xyyzz_yyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyz[i] * fz_0;

        tk_xxyyzz_yzz[i] = -2.0 * ts_yyzz_yzz[i] * fbe_0 * fz_0 + tk_yyzz_yzz[i] * fe_0 + tk_xyyzz_yzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yzz[i] * fz_0;

        tk_xxyyzz_zzz[i] = -2.0 * ts_yyzz_zzz[i] * fbe_0 * fz_0 + tk_yyzz_zzz[i] * fe_0 + tk_xyyzz_zzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_zzz[i] * fz_0;
    }

    // Set up 130-140 components of targeted buffer : IF

    auto tk_xxyzzz_xxx = pbuffer.data(idx_kin_if + 130);

    auto tk_xxyzzz_xxy = pbuffer.data(idx_kin_if + 131);

    auto tk_xxyzzz_xxz = pbuffer.data(idx_kin_if + 132);

    auto tk_xxyzzz_xyy = pbuffer.data(idx_kin_if + 133);

    auto tk_xxyzzz_xyz = pbuffer.data(idx_kin_if + 134);

    auto tk_xxyzzz_xzz = pbuffer.data(idx_kin_if + 135);

    auto tk_xxyzzz_yyy = pbuffer.data(idx_kin_if + 136);

    auto tk_xxyzzz_yyz = pbuffer.data(idx_kin_if + 137);

    auto tk_xxyzzz_yzz = pbuffer.data(idx_kin_if + 138);

    auto tk_xxyzzz_zzz = pbuffer.data(idx_kin_if + 139);

#pragma omp simd aligned(pa_y,              \
                             tk_xxyzzz_xxx, \
                             tk_xxyzzz_xxy, \
                             tk_xxyzzz_xxz, \
                             tk_xxyzzz_xyy, \
                             tk_xxyzzz_xyz, \
                             tk_xxyzzz_xzz, \
                             tk_xxyzzz_yyy, \
                             tk_xxyzzz_yyz, \
                             tk_xxyzzz_yzz, \
                             tk_xxyzzz_zzz, \
                             tk_xxzzz_xx,   \
                             tk_xxzzz_xxx,  \
                             tk_xxzzz_xxy,  \
                             tk_xxzzz_xxz,  \
                             tk_xxzzz_xy,   \
                             tk_xxzzz_xyy,  \
                             tk_xxzzz_xyz,  \
                             tk_xxzzz_xz,   \
                             tk_xxzzz_xzz,  \
                             tk_xxzzz_yy,   \
                             tk_xxzzz_yyy,  \
                             tk_xxzzz_yyz,  \
                             tk_xxzzz_yz,   \
                             tk_xxzzz_yzz,  \
                             tk_xxzzz_zz,   \
                             tk_xxzzz_zzz,  \
                             ts_xxyzzz_xxx, \
                             ts_xxyzzz_xxy, \
                             ts_xxyzzz_xxz, \
                             ts_xxyzzz_xyy, \
                             ts_xxyzzz_xyz, \
                             ts_xxyzzz_xzz, \
                             ts_xxyzzz_yyy, \
                             ts_xxyzzz_yyz, \
                             ts_xxyzzz_yzz, \
                             ts_xxyzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_xxx[i] = tk_xxzzz_xxx[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxx[i] * fz_0;

        tk_xxyzzz_xxy[i] = tk_xxzzz_xx[i] * fe_0 + tk_xxzzz_xxy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxy[i] * fz_0;

        tk_xxyzzz_xxz[i] = tk_xxzzz_xxz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxz[i] * fz_0;

        tk_xxyzzz_xyy[i] = 2.0 * tk_xxzzz_xy[i] * fe_0 + tk_xxzzz_xyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyy[i] * fz_0;

        tk_xxyzzz_xyz[i] = tk_xxzzz_xz[i] * fe_0 + tk_xxzzz_xyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyz[i] * fz_0;

        tk_xxyzzz_xzz[i] = tk_xxzzz_xzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xzz[i] * fz_0;

        tk_xxyzzz_yyy[i] = 3.0 * tk_xxzzz_yy[i] * fe_0 + tk_xxzzz_yyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyy[i] * fz_0;

        tk_xxyzzz_yyz[i] = 2.0 * tk_xxzzz_yz[i] * fe_0 + tk_xxzzz_yyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyz[i] * fz_0;

        tk_xxyzzz_yzz[i] = tk_xxzzz_zz[i] * fe_0 + tk_xxzzz_yzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yzz[i] * fz_0;

        tk_xxyzzz_zzz[i] = tk_xxzzz_zzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_zzz[i] * fz_0;
    }

    // Set up 140-150 components of targeted buffer : IF

    auto tk_xxzzzz_xxx = pbuffer.data(idx_kin_if + 140);

    auto tk_xxzzzz_xxy = pbuffer.data(idx_kin_if + 141);

    auto tk_xxzzzz_xxz = pbuffer.data(idx_kin_if + 142);

    auto tk_xxzzzz_xyy = pbuffer.data(idx_kin_if + 143);

    auto tk_xxzzzz_xyz = pbuffer.data(idx_kin_if + 144);

    auto tk_xxzzzz_xzz = pbuffer.data(idx_kin_if + 145);

    auto tk_xxzzzz_yyy = pbuffer.data(idx_kin_if + 146);

    auto tk_xxzzzz_yyz = pbuffer.data(idx_kin_if + 147);

    auto tk_xxzzzz_yzz = pbuffer.data(idx_kin_if + 148);

    auto tk_xxzzzz_zzz = pbuffer.data(idx_kin_if + 149);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xxzz_xxx,   \
                             tk_xxzz_xxy,   \
                             tk_xxzz_xyy,   \
                             tk_xxzzz_xxx,  \
                             tk_xxzzz_xxy,  \
                             tk_xxzzz_xyy,  \
                             tk_xxzzzz_xxx, \
                             tk_xxzzzz_xxy, \
                             tk_xxzzzz_xxz, \
                             tk_xxzzzz_xyy, \
                             tk_xxzzzz_xyz, \
                             tk_xxzzzz_xzz, \
                             tk_xxzzzz_yyy, \
                             tk_xxzzzz_yyz, \
                             tk_xxzzzz_yzz, \
                             tk_xxzzzz_zzz, \
                             tk_xzzzz_xxz,  \
                             tk_xzzzz_xyz,  \
                             tk_xzzzz_xz,   \
                             tk_xzzzz_xzz,  \
                             tk_xzzzz_yyy,  \
                             tk_xzzzz_yyz,  \
                             tk_xzzzz_yz,   \
                             tk_xzzzz_yzz,  \
                             tk_xzzzz_zz,   \
                             tk_xzzzz_zzz,  \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xzz,   \
                             tk_zzzz_yyy,   \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yzz,   \
                             tk_zzzz_zzz,   \
                             ts_xxzz_xxx,   \
                             ts_xxzz_xxy,   \
                             ts_xxzz_xyy,   \
                             ts_xxzzzz_xxx, \
                             ts_xxzzzz_xxy, \
                             ts_xxzzzz_xxz, \
                             ts_xxzzzz_xyy, \
                             ts_xxzzzz_xyz, \
                             ts_xxzzzz_xzz, \
                             ts_xxzzzz_yyy, \
                             ts_xxzzzz_yyz, \
                             ts_xxzzzz_yzz, \
                             ts_xxzzzz_zzz, \
                             ts_zzzz_xxz,   \
                             ts_zzzz_xyz,   \
                             ts_zzzz_xzz,   \
                             ts_zzzz_yyy,   \
                             ts_zzzz_yyz,   \
                             ts_zzzz_yzz,   \
                             ts_zzzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_xxx[i] =
            -6.0 * ts_xxzz_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxx[i] * fe_0 + tk_xxzzz_xxx[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxx[i] * fz_0;

        tk_xxzzzz_xxy[i] =
            -6.0 * ts_xxzz_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxy[i] * fe_0 + tk_xxzzz_xxy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxy[i] * fz_0;

        tk_xxzzzz_xxz[i] = -2.0 * ts_zzzz_xxz[i] * fbe_0 * fz_0 + tk_zzzz_xxz[i] * fe_0 + 2.0 * tk_xzzzz_xz[i] * fe_0 + tk_xzzzz_xxz[i] * pa_x[i] +
                           2.0 * ts_xxzzzz_xxz[i] * fz_0;

        tk_xxzzzz_xyy[i] =
            -6.0 * ts_xxzz_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyy[i] * fe_0 + tk_xxzzz_xyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xyy[i] * fz_0;

        tk_xxzzzz_xyz[i] = -2.0 * ts_zzzz_xyz[i] * fbe_0 * fz_0 + tk_zzzz_xyz[i] * fe_0 + tk_xzzzz_yz[i] * fe_0 + tk_xzzzz_xyz[i] * pa_x[i] +
                           2.0 * ts_xxzzzz_xyz[i] * fz_0;

        tk_xxzzzz_xzz[i] = -2.0 * ts_zzzz_xzz[i] * fbe_0 * fz_0 + tk_zzzz_xzz[i] * fe_0 + tk_xzzzz_zz[i] * fe_0 + tk_xzzzz_xzz[i] * pa_x[i] +
                           2.0 * ts_xxzzzz_xzz[i] * fz_0;

        tk_xxzzzz_yyy[i] = -2.0 * ts_zzzz_yyy[i] * fbe_0 * fz_0 + tk_zzzz_yyy[i] * fe_0 + tk_xzzzz_yyy[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyy[i] * fz_0;

        tk_xxzzzz_yyz[i] = -2.0 * ts_zzzz_yyz[i] * fbe_0 * fz_0 + tk_zzzz_yyz[i] * fe_0 + tk_xzzzz_yyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyz[i] * fz_0;

        tk_xxzzzz_yzz[i] = -2.0 * ts_zzzz_yzz[i] * fbe_0 * fz_0 + tk_zzzz_yzz[i] * fe_0 + tk_xzzzz_yzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yzz[i] * fz_0;

        tk_xxzzzz_zzz[i] = -2.0 * ts_zzzz_zzz[i] * fbe_0 * fz_0 + tk_zzzz_zzz[i] * fe_0 + tk_xzzzz_zzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_zzz[i] * fz_0;
    }

    // Set up 150-160 components of targeted buffer : IF

    auto tk_xyyyyy_xxx = pbuffer.data(idx_kin_if + 150);

    auto tk_xyyyyy_xxy = pbuffer.data(idx_kin_if + 151);

    auto tk_xyyyyy_xxz = pbuffer.data(idx_kin_if + 152);

    auto tk_xyyyyy_xyy = pbuffer.data(idx_kin_if + 153);

    auto tk_xyyyyy_xyz = pbuffer.data(idx_kin_if + 154);

    auto tk_xyyyyy_xzz = pbuffer.data(idx_kin_if + 155);

    auto tk_xyyyyy_yyy = pbuffer.data(idx_kin_if + 156);

    auto tk_xyyyyy_yyz = pbuffer.data(idx_kin_if + 157);

    auto tk_xyyyyy_yzz = pbuffer.data(idx_kin_if + 158);

    auto tk_xyyyyy_zzz = pbuffer.data(idx_kin_if + 159);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyyyy_xxx, \
                             tk_xyyyyy_xxy, \
                             tk_xyyyyy_xxz, \
                             tk_xyyyyy_xyy, \
                             tk_xyyyyy_xyz, \
                             tk_xyyyyy_xzz, \
                             tk_xyyyyy_yyy, \
                             tk_xyyyyy_yyz, \
                             tk_xyyyyy_yzz, \
                             tk_xyyyyy_zzz, \
                             tk_yyyyy_xx,   \
                             tk_yyyyy_xxx,  \
                             tk_yyyyy_xxy,  \
                             tk_yyyyy_xxz,  \
                             tk_yyyyy_xy,   \
                             tk_yyyyy_xyy,  \
                             tk_yyyyy_xyz,  \
                             tk_yyyyy_xz,   \
                             tk_yyyyy_xzz,  \
                             tk_yyyyy_yy,   \
                             tk_yyyyy_yyy,  \
                             tk_yyyyy_yyz,  \
                             tk_yyyyy_yz,   \
                             tk_yyyyy_yzz,  \
                             tk_yyyyy_zz,   \
                             tk_yyyyy_zzz,  \
                             ts_xyyyyy_xxx, \
                             ts_xyyyyy_xxy, \
                             ts_xyyyyy_xxz, \
                             ts_xyyyyy_xyy, \
                             ts_xyyyyy_xyz, \
                             ts_xyyyyy_xzz, \
                             ts_xyyyyy_yyy, \
                             ts_xyyyyy_yyz, \
                             ts_xyyyyy_yzz, \
                             ts_xyyyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_xxx[i] = 3.0 * tk_yyyyy_xx[i] * fe_0 + tk_yyyyy_xxx[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxx[i] * fz_0;

        tk_xyyyyy_xxy[i] = 2.0 * tk_yyyyy_xy[i] * fe_0 + tk_yyyyy_xxy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxy[i] * fz_0;

        tk_xyyyyy_xxz[i] = 2.0 * tk_yyyyy_xz[i] * fe_0 + tk_yyyyy_xxz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxz[i] * fz_0;

        tk_xyyyyy_xyy[i] = tk_yyyyy_yy[i] * fe_0 + tk_yyyyy_xyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyy[i] * fz_0;

        tk_xyyyyy_xyz[i] = tk_yyyyy_yz[i] * fe_0 + tk_yyyyy_xyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyz[i] * fz_0;

        tk_xyyyyy_xzz[i] = tk_yyyyy_zz[i] * fe_0 + tk_yyyyy_xzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xzz[i] * fz_0;

        tk_xyyyyy_yyy[i] = tk_yyyyy_yyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyy[i] * fz_0;

        tk_xyyyyy_yyz[i] = tk_yyyyy_yyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyz[i] * fz_0;

        tk_xyyyyy_yzz[i] = tk_yyyyy_yzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yzz[i] * fz_0;

        tk_xyyyyy_zzz[i] = tk_yyyyy_zzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_zzz[i] * fz_0;
    }

    // Set up 160-170 components of targeted buffer : IF

    auto tk_xyyyyz_xxx = pbuffer.data(idx_kin_if + 160);

    auto tk_xyyyyz_xxy = pbuffer.data(idx_kin_if + 161);

    auto tk_xyyyyz_xxz = pbuffer.data(idx_kin_if + 162);

    auto tk_xyyyyz_xyy = pbuffer.data(idx_kin_if + 163);

    auto tk_xyyyyz_xyz = pbuffer.data(idx_kin_if + 164);

    auto tk_xyyyyz_xzz = pbuffer.data(idx_kin_if + 165);

    auto tk_xyyyyz_yyy = pbuffer.data(idx_kin_if + 166);

    auto tk_xyyyyz_yyz = pbuffer.data(idx_kin_if + 167);

    auto tk_xyyyyz_yzz = pbuffer.data(idx_kin_if + 168);

    auto tk_xyyyyz_zzz = pbuffer.data(idx_kin_if + 169);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xyyyy_xxx,  \
                             tk_xyyyy_xxy,  \
                             tk_xyyyy_xyy,  \
                             tk_xyyyyz_xxx, \
                             tk_xyyyyz_xxy, \
                             tk_xyyyyz_xxz, \
                             tk_xyyyyz_xyy, \
                             tk_xyyyyz_xyz, \
                             tk_xyyyyz_xzz, \
                             tk_xyyyyz_yyy, \
                             tk_xyyyyz_yyz, \
                             tk_xyyyyz_yzz, \
                             tk_xyyyyz_zzz, \
                             tk_yyyyz_xxz,  \
                             tk_yyyyz_xyz,  \
                             tk_yyyyz_xz,   \
                             tk_yyyyz_xzz,  \
                             tk_yyyyz_yyy,  \
                             tk_yyyyz_yyz,  \
                             tk_yyyyz_yz,   \
                             tk_yyyyz_yzz,  \
                             tk_yyyyz_zz,   \
                             tk_yyyyz_zzz,  \
                             ts_xyyyyz_xxx, \
                             ts_xyyyyz_xxy, \
                             ts_xyyyyz_xxz, \
                             ts_xyyyyz_xyy, \
                             ts_xyyyyz_xyz, \
                             ts_xyyyyz_xzz, \
                             ts_xyyyyz_yyy, \
                             ts_xyyyyz_yyz, \
                             ts_xyyyyz_yzz, \
                             ts_xyyyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyz_xxx[i] = tk_xyyyy_xxx[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxx[i] * fz_0;

        tk_xyyyyz_xxy[i] = tk_xyyyy_xxy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxy[i] * fz_0;

        tk_xyyyyz_xxz[i] = 2.0 * tk_yyyyz_xz[i] * fe_0 + tk_yyyyz_xxz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxz[i] * fz_0;

        tk_xyyyyz_xyy[i] = tk_xyyyy_xyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xyy[i] * fz_0;

        tk_xyyyyz_xyz[i] = tk_yyyyz_yz[i] * fe_0 + tk_yyyyz_xyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyz[i] * fz_0;

        tk_xyyyyz_xzz[i] = tk_yyyyz_zz[i] * fe_0 + tk_yyyyz_xzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xzz[i] * fz_0;

        tk_xyyyyz_yyy[i] = tk_yyyyz_yyy[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyy[i] * fz_0;

        tk_xyyyyz_yyz[i] = tk_yyyyz_yyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyz[i] * fz_0;

        tk_xyyyyz_yzz[i] = tk_yyyyz_yzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yzz[i] * fz_0;

        tk_xyyyyz_zzz[i] = tk_yyyyz_zzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_zzz[i] * fz_0;
    }

    // Set up 170-180 components of targeted buffer : IF

    auto tk_xyyyzz_xxx = pbuffer.data(idx_kin_if + 170);

    auto tk_xyyyzz_xxy = pbuffer.data(idx_kin_if + 171);

    auto tk_xyyyzz_xxz = pbuffer.data(idx_kin_if + 172);

    auto tk_xyyyzz_xyy = pbuffer.data(idx_kin_if + 173);

    auto tk_xyyyzz_xyz = pbuffer.data(idx_kin_if + 174);

    auto tk_xyyyzz_xzz = pbuffer.data(idx_kin_if + 175);

    auto tk_xyyyzz_yyy = pbuffer.data(idx_kin_if + 176);

    auto tk_xyyyzz_yyz = pbuffer.data(idx_kin_if + 177);

    auto tk_xyyyzz_yzz = pbuffer.data(idx_kin_if + 178);

    auto tk_xyyyzz_zzz = pbuffer.data(idx_kin_if + 179);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyyzz_xxx, \
                             tk_xyyyzz_xxy, \
                             tk_xyyyzz_xxz, \
                             tk_xyyyzz_xyy, \
                             tk_xyyyzz_xyz, \
                             tk_xyyyzz_xzz, \
                             tk_xyyyzz_yyy, \
                             tk_xyyyzz_yyz, \
                             tk_xyyyzz_yzz, \
                             tk_xyyyzz_zzz, \
                             tk_yyyzz_xx,   \
                             tk_yyyzz_xxx,  \
                             tk_yyyzz_xxy,  \
                             tk_yyyzz_xxz,  \
                             tk_yyyzz_xy,   \
                             tk_yyyzz_xyy,  \
                             tk_yyyzz_xyz,  \
                             tk_yyyzz_xz,   \
                             tk_yyyzz_xzz,  \
                             tk_yyyzz_yy,   \
                             tk_yyyzz_yyy,  \
                             tk_yyyzz_yyz,  \
                             tk_yyyzz_yz,   \
                             tk_yyyzz_yzz,  \
                             tk_yyyzz_zz,   \
                             tk_yyyzz_zzz,  \
                             ts_xyyyzz_xxx, \
                             ts_xyyyzz_xxy, \
                             ts_xyyyzz_xxz, \
                             ts_xyyyzz_xyy, \
                             ts_xyyyzz_xyz, \
                             ts_xyyyzz_xzz, \
                             ts_xyyyzz_yyy, \
                             ts_xyyyzz_yyz, \
                             ts_xyyyzz_yzz, \
                             ts_xyyyzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_xxx[i] = 3.0 * tk_yyyzz_xx[i] * fe_0 + tk_yyyzz_xxx[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxx[i] * fz_0;

        tk_xyyyzz_xxy[i] = 2.0 * tk_yyyzz_xy[i] * fe_0 + tk_yyyzz_xxy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxy[i] * fz_0;

        tk_xyyyzz_xxz[i] = 2.0 * tk_yyyzz_xz[i] * fe_0 + tk_yyyzz_xxz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxz[i] * fz_0;

        tk_xyyyzz_xyy[i] = tk_yyyzz_yy[i] * fe_0 + tk_yyyzz_xyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyy[i] * fz_0;

        tk_xyyyzz_xyz[i] = tk_yyyzz_yz[i] * fe_0 + tk_yyyzz_xyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyz[i] * fz_0;

        tk_xyyyzz_xzz[i] = tk_yyyzz_zz[i] * fe_0 + tk_yyyzz_xzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xzz[i] * fz_0;

        tk_xyyyzz_yyy[i] = tk_yyyzz_yyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyy[i] * fz_0;

        tk_xyyyzz_yyz[i] = tk_yyyzz_yyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyz[i] * fz_0;

        tk_xyyyzz_yzz[i] = tk_yyyzz_yzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yzz[i] * fz_0;

        tk_xyyyzz_zzz[i] = tk_yyyzz_zzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_zzz[i] * fz_0;
    }

    // Set up 180-190 components of targeted buffer : IF

    auto tk_xyyzzz_xxx = pbuffer.data(idx_kin_if + 180);

    auto tk_xyyzzz_xxy = pbuffer.data(idx_kin_if + 181);

    auto tk_xyyzzz_xxz = pbuffer.data(idx_kin_if + 182);

    auto tk_xyyzzz_xyy = pbuffer.data(idx_kin_if + 183);

    auto tk_xyyzzz_xyz = pbuffer.data(idx_kin_if + 184);

    auto tk_xyyzzz_xzz = pbuffer.data(idx_kin_if + 185);

    auto tk_xyyzzz_yyy = pbuffer.data(idx_kin_if + 186);

    auto tk_xyyzzz_yyz = pbuffer.data(idx_kin_if + 187);

    auto tk_xyyzzz_yzz = pbuffer.data(idx_kin_if + 188);

    auto tk_xyyzzz_zzz = pbuffer.data(idx_kin_if + 189);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyzzz_xxx, \
                             tk_xyyzzz_xxy, \
                             tk_xyyzzz_xxz, \
                             tk_xyyzzz_xyy, \
                             tk_xyyzzz_xyz, \
                             tk_xyyzzz_xzz, \
                             tk_xyyzzz_yyy, \
                             tk_xyyzzz_yyz, \
                             tk_xyyzzz_yzz, \
                             tk_xyyzzz_zzz, \
                             tk_yyzzz_xx,   \
                             tk_yyzzz_xxx,  \
                             tk_yyzzz_xxy,  \
                             tk_yyzzz_xxz,  \
                             tk_yyzzz_xy,   \
                             tk_yyzzz_xyy,  \
                             tk_yyzzz_xyz,  \
                             tk_yyzzz_xz,   \
                             tk_yyzzz_xzz,  \
                             tk_yyzzz_yy,   \
                             tk_yyzzz_yyy,  \
                             tk_yyzzz_yyz,  \
                             tk_yyzzz_yz,   \
                             tk_yyzzz_yzz,  \
                             tk_yyzzz_zz,   \
                             tk_yyzzz_zzz,  \
                             ts_xyyzzz_xxx, \
                             ts_xyyzzz_xxy, \
                             ts_xyyzzz_xxz, \
                             ts_xyyzzz_xyy, \
                             ts_xyyzzz_xyz, \
                             ts_xyyzzz_xzz, \
                             ts_xyyzzz_yyy, \
                             ts_xyyzzz_yyz, \
                             ts_xyyzzz_yzz, \
                             ts_xyyzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_xxx[i] = 3.0 * tk_yyzzz_xx[i] * fe_0 + tk_yyzzz_xxx[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxx[i] * fz_0;

        tk_xyyzzz_xxy[i] = 2.0 * tk_yyzzz_xy[i] * fe_0 + tk_yyzzz_xxy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxy[i] * fz_0;

        tk_xyyzzz_xxz[i] = 2.0 * tk_yyzzz_xz[i] * fe_0 + tk_yyzzz_xxz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxz[i] * fz_0;

        tk_xyyzzz_xyy[i] = tk_yyzzz_yy[i] * fe_0 + tk_yyzzz_xyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyy[i] * fz_0;

        tk_xyyzzz_xyz[i] = tk_yyzzz_yz[i] * fe_0 + tk_yyzzz_xyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyz[i] * fz_0;

        tk_xyyzzz_xzz[i] = tk_yyzzz_zz[i] * fe_0 + tk_yyzzz_xzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xzz[i] * fz_0;

        tk_xyyzzz_yyy[i] = tk_yyzzz_yyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyy[i] * fz_0;

        tk_xyyzzz_yyz[i] = tk_yyzzz_yyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyz[i] * fz_0;

        tk_xyyzzz_yzz[i] = tk_yyzzz_yzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yzz[i] * fz_0;

        tk_xyyzzz_zzz[i] = tk_yyzzz_zzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_zzz[i] * fz_0;
    }

    // Set up 190-200 components of targeted buffer : IF

    auto tk_xyzzzz_xxx = pbuffer.data(idx_kin_if + 190);

    auto tk_xyzzzz_xxy = pbuffer.data(idx_kin_if + 191);

    auto tk_xyzzzz_xxz = pbuffer.data(idx_kin_if + 192);

    auto tk_xyzzzz_xyy = pbuffer.data(idx_kin_if + 193);

    auto tk_xyzzzz_xyz = pbuffer.data(idx_kin_if + 194);

    auto tk_xyzzzz_xzz = pbuffer.data(idx_kin_if + 195);

    auto tk_xyzzzz_yyy = pbuffer.data(idx_kin_if + 196);

    auto tk_xyzzzz_yyz = pbuffer.data(idx_kin_if + 197);

    auto tk_xyzzzz_yzz = pbuffer.data(idx_kin_if + 198);

    auto tk_xyzzzz_zzz = pbuffer.data(idx_kin_if + 199);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xyzzzz_xxx, \
                             tk_xyzzzz_xxy, \
                             tk_xyzzzz_xxz, \
                             tk_xyzzzz_xyy, \
                             tk_xyzzzz_xyz, \
                             tk_xyzzzz_xzz, \
                             tk_xyzzzz_yyy, \
                             tk_xyzzzz_yyz, \
                             tk_xyzzzz_yzz, \
                             tk_xyzzzz_zzz, \
                             tk_xzzzz_xxx,  \
                             tk_xzzzz_xxz,  \
                             tk_xzzzz_xzz,  \
                             tk_yzzzz_xxy,  \
                             tk_yzzzz_xy,   \
                             tk_yzzzz_xyy,  \
                             tk_yzzzz_xyz,  \
                             tk_yzzzz_yy,   \
                             tk_yzzzz_yyy,  \
                             tk_yzzzz_yyz,  \
                             tk_yzzzz_yz,   \
                             tk_yzzzz_yzz,  \
                             tk_yzzzz_zzz,  \
                             ts_xyzzzz_xxx, \
                             ts_xyzzzz_xxy, \
                             ts_xyzzzz_xxz, \
                             ts_xyzzzz_xyy, \
                             ts_xyzzzz_xyz, \
                             ts_xyzzzz_xzz, \
                             ts_xyzzzz_yyy, \
                             ts_xyzzzz_yyz, \
                             ts_xyzzzz_yzz, \
                             ts_xyzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzzz_xxx[i] = tk_xzzzz_xxx[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxx[i] * fz_0;

        tk_xyzzzz_xxy[i] = 2.0 * tk_yzzzz_xy[i] * fe_0 + tk_yzzzz_xxy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxy[i] * fz_0;

        tk_xyzzzz_xxz[i] = tk_xzzzz_xxz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxz[i] * fz_0;

        tk_xyzzzz_xyy[i] = tk_yzzzz_yy[i] * fe_0 + tk_yzzzz_xyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyy[i] * fz_0;

        tk_xyzzzz_xyz[i] = tk_yzzzz_yz[i] * fe_0 + tk_yzzzz_xyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyz[i] * fz_0;

        tk_xyzzzz_xzz[i] = tk_xzzzz_xzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xzz[i] * fz_0;

        tk_xyzzzz_yyy[i] = tk_yzzzz_yyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyy[i] * fz_0;

        tk_xyzzzz_yyz[i] = tk_yzzzz_yyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyz[i] * fz_0;

        tk_xyzzzz_yzz[i] = tk_yzzzz_yzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yzz[i] * fz_0;

        tk_xyzzzz_zzz[i] = tk_yzzzz_zzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_zzz[i] * fz_0;
    }

    // Set up 200-210 components of targeted buffer : IF

    auto tk_xzzzzz_xxx = pbuffer.data(idx_kin_if + 200);

    auto tk_xzzzzz_xxy = pbuffer.data(idx_kin_if + 201);

    auto tk_xzzzzz_xxz = pbuffer.data(idx_kin_if + 202);

    auto tk_xzzzzz_xyy = pbuffer.data(idx_kin_if + 203);

    auto tk_xzzzzz_xyz = pbuffer.data(idx_kin_if + 204);

    auto tk_xzzzzz_xzz = pbuffer.data(idx_kin_if + 205);

    auto tk_xzzzzz_yyy = pbuffer.data(idx_kin_if + 206);

    auto tk_xzzzzz_yyz = pbuffer.data(idx_kin_if + 207);

    auto tk_xzzzzz_yzz = pbuffer.data(idx_kin_if + 208);

    auto tk_xzzzzz_zzz = pbuffer.data(idx_kin_if + 209);

#pragma omp simd aligned(pa_x,              \
                             tk_xzzzzz_xxx, \
                             tk_xzzzzz_xxy, \
                             tk_xzzzzz_xxz, \
                             tk_xzzzzz_xyy, \
                             tk_xzzzzz_xyz, \
                             tk_xzzzzz_xzz, \
                             tk_xzzzzz_yyy, \
                             tk_xzzzzz_yyz, \
                             tk_xzzzzz_yzz, \
                             tk_xzzzzz_zzz, \
                             tk_zzzzz_xx,   \
                             tk_zzzzz_xxx,  \
                             tk_zzzzz_xxy,  \
                             tk_zzzzz_xxz,  \
                             tk_zzzzz_xy,   \
                             tk_zzzzz_xyy,  \
                             tk_zzzzz_xyz,  \
                             tk_zzzzz_xz,   \
                             tk_zzzzz_xzz,  \
                             tk_zzzzz_yy,   \
                             tk_zzzzz_yyy,  \
                             tk_zzzzz_yyz,  \
                             tk_zzzzz_yz,   \
                             tk_zzzzz_yzz,  \
                             tk_zzzzz_zz,   \
                             tk_zzzzz_zzz,  \
                             ts_xzzzzz_xxx, \
                             ts_xzzzzz_xxy, \
                             ts_xzzzzz_xxz, \
                             ts_xzzzzz_xyy, \
                             ts_xzzzzz_xyz, \
                             ts_xzzzzz_xzz, \
                             ts_xzzzzz_yyy, \
                             ts_xzzzzz_yyz, \
                             ts_xzzzzz_yzz, \
                             ts_xzzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_xxx[i] = 3.0 * tk_zzzzz_xx[i] * fe_0 + tk_zzzzz_xxx[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxx[i] * fz_0;

        tk_xzzzzz_xxy[i] = 2.0 * tk_zzzzz_xy[i] * fe_0 + tk_zzzzz_xxy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxy[i] * fz_0;

        tk_xzzzzz_xxz[i] = 2.0 * tk_zzzzz_xz[i] * fe_0 + tk_zzzzz_xxz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxz[i] * fz_0;

        tk_xzzzzz_xyy[i] = tk_zzzzz_yy[i] * fe_0 + tk_zzzzz_xyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyy[i] * fz_0;

        tk_xzzzzz_xyz[i] = tk_zzzzz_yz[i] * fe_0 + tk_zzzzz_xyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyz[i] * fz_0;

        tk_xzzzzz_xzz[i] = tk_zzzzz_zz[i] * fe_0 + tk_zzzzz_xzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xzz[i] * fz_0;

        tk_xzzzzz_yyy[i] = tk_zzzzz_yyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyy[i] * fz_0;

        tk_xzzzzz_yyz[i] = tk_zzzzz_yyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyz[i] * fz_0;

        tk_xzzzzz_yzz[i] = tk_zzzzz_yzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yzz[i] * fz_0;

        tk_xzzzzz_zzz[i] = tk_zzzzz_zzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_zzz[i] * fz_0;
    }

    // Set up 210-220 components of targeted buffer : IF

    auto tk_yyyyyy_xxx = pbuffer.data(idx_kin_if + 210);

    auto tk_yyyyyy_xxy = pbuffer.data(idx_kin_if + 211);

    auto tk_yyyyyy_xxz = pbuffer.data(idx_kin_if + 212);

    auto tk_yyyyyy_xyy = pbuffer.data(idx_kin_if + 213);

    auto tk_yyyyyy_xyz = pbuffer.data(idx_kin_if + 214);

    auto tk_yyyyyy_xzz = pbuffer.data(idx_kin_if + 215);

    auto tk_yyyyyy_yyy = pbuffer.data(idx_kin_if + 216);

    auto tk_yyyyyy_yyz = pbuffer.data(idx_kin_if + 217);

    auto tk_yyyyyy_yzz = pbuffer.data(idx_kin_if + 218);

    auto tk_yyyyyy_zzz = pbuffer.data(idx_kin_if + 219);

#pragma omp simd aligned(pa_y,              \
                             tk_yyyy_xxx,   \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xxz,   \
                             tk_yyyy_xyy,   \
                             tk_yyyy_xyz,   \
                             tk_yyyy_xzz,   \
                             tk_yyyy_yyy,   \
                             tk_yyyy_yyz,   \
                             tk_yyyy_yzz,   \
                             tk_yyyy_zzz,   \
                             tk_yyyyy_xx,   \
                             tk_yyyyy_xxx,  \
                             tk_yyyyy_xxy,  \
                             tk_yyyyy_xxz,  \
                             tk_yyyyy_xy,   \
                             tk_yyyyy_xyy,  \
                             tk_yyyyy_xyz,  \
                             tk_yyyyy_xz,   \
                             tk_yyyyy_xzz,  \
                             tk_yyyyy_yy,   \
                             tk_yyyyy_yyy,  \
                             tk_yyyyy_yyz,  \
                             tk_yyyyy_yz,   \
                             tk_yyyyy_yzz,  \
                             tk_yyyyy_zz,   \
                             tk_yyyyy_zzz,  \
                             tk_yyyyyy_xxx, \
                             tk_yyyyyy_xxy, \
                             tk_yyyyyy_xxz, \
                             tk_yyyyyy_xyy, \
                             tk_yyyyyy_xyz, \
                             tk_yyyyyy_xzz, \
                             tk_yyyyyy_yyy, \
                             tk_yyyyyy_yyz, \
                             tk_yyyyyy_yzz, \
                             tk_yyyyyy_zzz, \
                             ts_yyyy_xxx,   \
                             ts_yyyy_xxy,   \
                             ts_yyyy_xxz,   \
                             ts_yyyy_xyy,   \
                             ts_yyyy_xyz,   \
                             ts_yyyy_xzz,   \
                             ts_yyyy_yyy,   \
                             ts_yyyy_yyz,   \
                             ts_yyyy_yzz,   \
                             ts_yyyy_zzz,   \
                             ts_yyyyyy_xxx, \
                             ts_yyyyyy_xxy, \
                             ts_yyyyyy_xxz, \
                             ts_yyyyyy_xyy, \
                             ts_yyyyyy_xyz, \
                             ts_yyyyyy_xzz, \
                             ts_yyyyyy_yyy, \
                             ts_yyyyyy_yyz, \
                             ts_yyyyyy_yzz, \
                             ts_yyyyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_xxx[i] =
            -10.0 * ts_yyyy_xxx[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxx[i] * fe_0 + tk_yyyyy_xxx[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxx[i] * fz_0;

        tk_yyyyyy_xxy[i] = -10.0 * ts_yyyy_xxy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxy[i] * fe_0 + tk_yyyyy_xx[i] * fe_0 + tk_yyyyy_xxy[i] * pa_y[i] +
                           2.0 * ts_yyyyyy_xxy[i] * fz_0;

        tk_yyyyyy_xxz[i] =
            -10.0 * ts_yyyy_xxz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxz[i] * fe_0 + tk_yyyyy_xxz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxz[i] * fz_0;

        tk_yyyyyy_xyy[i] = -10.0 * ts_yyyy_xyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyy[i] * fe_0 + 2.0 * tk_yyyyy_xy[i] * fe_0 +
                           tk_yyyyy_xyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyy[i] * fz_0;

        tk_yyyyyy_xyz[i] = -10.0 * ts_yyyy_xyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyz[i] * fe_0 + tk_yyyyy_xz[i] * fe_0 + tk_yyyyy_xyz[i] * pa_y[i] +
                           2.0 * ts_yyyyyy_xyz[i] * fz_0;

        tk_yyyyyy_xzz[i] =
            -10.0 * ts_yyyy_xzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xzz[i] * fe_0 + tk_yyyyy_xzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xzz[i] * fz_0;

        tk_yyyyyy_yyy[i] = -10.0 * ts_yyyy_yyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyy[i] * fe_0 + 3.0 * tk_yyyyy_yy[i] * fe_0 +
                           tk_yyyyy_yyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyy[i] * fz_0;

        tk_yyyyyy_yyz[i] = -10.0 * ts_yyyy_yyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyz[i] * fe_0 + 2.0 * tk_yyyyy_yz[i] * fe_0 +
                           tk_yyyyy_yyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyz[i] * fz_0;

        tk_yyyyyy_yzz[i] = -10.0 * ts_yyyy_yzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yzz[i] * fe_0 + tk_yyyyy_zz[i] * fe_0 + tk_yyyyy_yzz[i] * pa_y[i] +
                           2.0 * ts_yyyyyy_yzz[i] * fz_0;

        tk_yyyyyy_zzz[i] =
            -10.0 * ts_yyyy_zzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_zzz[i] * fe_0 + tk_yyyyy_zzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_zzz[i] * fz_0;
    }

    // Set up 220-230 components of targeted buffer : IF

    auto tk_yyyyyz_xxx = pbuffer.data(idx_kin_if + 220);

    auto tk_yyyyyz_xxy = pbuffer.data(idx_kin_if + 221);

    auto tk_yyyyyz_xxz = pbuffer.data(idx_kin_if + 222);

    auto tk_yyyyyz_xyy = pbuffer.data(idx_kin_if + 223);

    auto tk_yyyyyz_xyz = pbuffer.data(idx_kin_if + 224);

    auto tk_yyyyyz_xzz = pbuffer.data(idx_kin_if + 225);

    auto tk_yyyyyz_yyy = pbuffer.data(idx_kin_if + 226);

    auto tk_yyyyyz_yyz = pbuffer.data(idx_kin_if + 227);

    auto tk_yyyyyz_yzz = pbuffer.data(idx_kin_if + 228);

    auto tk_yyyyyz_zzz = pbuffer.data(idx_kin_if + 229);

#pragma omp simd aligned(pa_z,              \
                             tk_yyyyy_xx,   \
                             tk_yyyyy_xxx,  \
                             tk_yyyyy_xxy,  \
                             tk_yyyyy_xxz,  \
                             tk_yyyyy_xy,   \
                             tk_yyyyy_xyy,  \
                             tk_yyyyy_xyz,  \
                             tk_yyyyy_xz,   \
                             tk_yyyyy_xzz,  \
                             tk_yyyyy_yy,   \
                             tk_yyyyy_yyy,  \
                             tk_yyyyy_yyz,  \
                             tk_yyyyy_yz,   \
                             tk_yyyyy_yzz,  \
                             tk_yyyyy_zz,   \
                             tk_yyyyy_zzz,  \
                             tk_yyyyyz_xxx, \
                             tk_yyyyyz_xxy, \
                             tk_yyyyyz_xxz, \
                             tk_yyyyyz_xyy, \
                             tk_yyyyyz_xyz, \
                             tk_yyyyyz_xzz, \
                             tk_yyyyyz_yyy, \
                             tk_yyyyyz_yyz, \
                             tk_yyyyyz_yzz, \
                             tk_yyyyyz_zzz, \
                             ts_yyyyyz_xxx, \
                             ts_yyyyyz_xxy, \
                             ts_yyyyyz_xxz, \
                             ts_yyyyyz_xyy, \
                             ts_yyyyyz_xyz, \
                             ts_yyyyyz_xzz, \
                             ts_yyyyyz_yyy, \
                             ts_yyyyyz_yyz, \
                             ts_yyyyyz_yzz, \
                             ts_yyyyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_xxx[i] = tk_yyyyy_xxx[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxx[i] * fz_0;

        tk_yyyyyz_xxy[i] = tk_yyyyy_xxy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxy[i] * fz_0;

        tk_yyyyyz_xxz[i] = tk_yyyyy_xx[i] * fe_0 + tk_yyyyy_xxz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxz[i] * fz_0;

        tk_yyyyyz_xyy[i] = tk_yyyyy_xyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyy[i] * fz_0;

        tk_yyyyyz_xyz[i] = tk_yyyyy_xy[i] * fe_0 + tk_yyyyy_xyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyz[i] * fz_0;

        tk_yyyyyz_xzz[i] = 2.0 * tk_yyyyy_xz[i] * fe_0 + tk_yyyyy_xzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xzz[i] * fz_0;

        tk_yyyyyz_yyy[i] = tk_yyyyy_yyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyy[i] * fz_0;

        tk_yyyyyz_yyz[i] = tk_yyyyy_yy[i] * fe_0 + tk_yyyyy_yyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyz[i] * fz_0;

        tk_yyyyyz_yzz[i] = 2.0 * tk_yyyyy_yz[i] * fe_0 + tk_yyyyy_yzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yzz[i] * fz_0;

        tk_yyyyyz_zzz[i] = 3.0 * tk_yyyyy_zz[i] * fe_0 + tk_yyyyy_zzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_zzz[i] * fz_0;
    }

    // Set up 230-240 components of targeted buffer : IF

    auto tk_yyyyzz_xxx = pbuffer.data(idx_kin_if + 230);

    auto tk_yyyyzz_xxy = pbuffer.data(idx_kin_if + 231);

    auto tk_yyyyzz_xxz = pbuffer.data(idx_kin_if + 232);

    auto tk_yyyyzz_xyy = pbuffer.data(idx_kin_if + 233);

    auto tk_yyyyzz_xyz = pbuffer.data(idx_kin_if + 234);

    auto tk_yyyyzz_xzz = pbuffer.data(idx_kin_if + 235);

    auto tk_yyyyzz_yyy = pbuffer.data(idx_kin_if + 236);

    auto tk_yyyyzz_yyz = pbuffer.data(idx_kin_if + 237);

    auto tk_yyyyzz_yzz = pbuffer.data(idx_kin_if + 238);

    auto tk_yyyyzz_zzz = pbuffer.data(idx_kin_if + 239);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xyy,   \
                             tk_yyyy_yyy,   \
                             tk_yyyyz_xxy,  \
                             tk_yyyyz_xyy,  \
                             tk_yyyyz_yyy,  \
                             tk_yyyyzz_xxx, \
                             tk_yyyyzz_xxy, \
                             tk_yyyyzz_xxz, \
                             tk_yyyyzz_xyy, \
                             tk_yyyyzz_xyz, \
                             tk_yyyyzz_xzz, \
                             tk_yyyyzz_yyy, \
                             tk_yyyyzz_yyz, \
                             tk_yyyyzz_yzz, \
                             tk_yyyyzz_zzz, \
                             tk_yyyzz_xxx,  \
                             tk_yyyzz_xxz,  \
                             tk_yyyzz_xyz,  \
                             tk_yyyzz_xz,   \
                             tk_yyyzz_xzz,  \
                             tk_yyyzz_yyz,  \
                             tk_yyyzz_yz,   \
                             tk_yyyzz_yzz,  \
                             tk_yyyzz_zz,   \
                             tk_yyyzz_zzz,  \
                             tk_yyzz_xxx,   \
                             tk_yyzz_xxz,   \
                             tk_yyzz_xyz,   \
                             tk_yyzz_xzz,   \
                             tk_yyzz_yyz,   \
                             tk_yyzz_yzz,   \
                             tk_yyzz_zzz,   \
                             ts_yyyy_xxy,   \
                             ts_yyyy_xyy,   \
                             ts_yyyy_yyy,   \
                             ts_yyyyzz_xxx, \
                             ts_yyyyzz_xxy, \
                             ts_yyyyzz_xxz, \
                             ts_yyyyzz_xyy, \
                             ts_yyyyzz_xyz, \
                             ts_yyyyzz_xzz, \
                             ts_yyyyzz_yyy, \
                             ts_yyyyzz_yyz, \
                             ts_yyyyzz_yzz, \
                             ts_yyyyzz_zzz, \
                             ts_yyzz_xxx,   \
                             ts_yyzz_xxz,   \
                             ts_yyzz_xyz,   \
                             ts_yyzz_xzz,   \
                             ts_yyzz_yyz,   \
                             ts_yyzz_yzz,   \
                             ts_yyzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_xxx[i] =
            -6.0 * ts_yyzz_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxx[i] * fe_0 + tk_yyyzz_xxx[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxx[i] * fz_0;

        tk_yyyyzz_xxy[i] = -2.0 * ts_yyyy_xxy[i] * fbe_0 * fz_0 + tk_yyyy_xxy[i] * fe_0 + tk_yyyyz_xxy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxy[i] * fz_0;

        tk_yyyyzz_xxz[i] =
            -6.0 * ts_yyzz_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxz[i] * fe_0 + tk_yyyzz_xxz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxz[i] * fz_0;

        tk_yyyyzz_xyy[i] = -2.0 * ts_yyyy_xyy[i] * fbe_0 * fz_0 + tk_yyyy_xyy[i] * fe_0 + tk_yyyyz_xyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xyy[i] * fz_0;

        tk_yyyyzz_xyz[i] = -6.0 * ts_yyzz_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyz[i] * fe_0 + tk_yyyzz_xz[i] * fe_0 + tk_yyyzz_xyz[i] * pa_y[i] +
                           2.0 * ts_yyyyzz_xyz[i] * fz_0;

        tk_yyyyzz_xzz[i] =
            -6.0 * ts_yyzz_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xzz[i] * fe_0 + tk_yyyzz_xzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xzz[i] * fz_0;

        tk_yyyyzz_yyy[i] = -2.0 * ts_yyyy_yyy[i] * fbe_0 * fz_0 + tk_yyyy_yyy[i] * fe_0 + tk_yyyyz_yyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_yyy[i] * fz_0;

        tk_yyyyzz_yyz[i] = -6.0 * ts_yyzz_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyz[i] * fe_0 + 2.0 * tk_yyyzz_yz[i] * fe_0 +
                           tk_yyyzz_yyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyz[i] * fz_0;

        tk_yyyyzz_yzz[i] = -6.0 * ts_yyzz_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yzz[i] * fe_0 + tk_yyyzz_zz[i] * fe_0 + tk_yyyzz_yzz[i] * pa_y[i] +
                           2.0 * ts_yyyyzz_yzz[i] * fz_0;

        tk_yyyyzz_zzz[i] =
            -6.0 * ts_yyzz_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_zzz[i] * fe_0 + tk_yyyzz_zzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_zzz[i] * fz_0;
    }

    // Set up 240-250 components of targeted buffer : IF

    auto tk_yyyzzz_xxx = pbuffer.data(idx_kin_if + 240);

    auto tk_yyyzzz_xxy = pbuffer.data(idx_kin_if + 241);

    auto tk_yyyzzz_xxz = pbuffer.data(idx_kin_if + 242);

    auto tk_yyyzzz_xyy = pbuffer.data(idx_kin_if + 243);

    auto tk_yyyzzz_xyz = pbuffer.data(idx_kin_if + 244);

    auto tk_yyyzzz_xzz = pbuffer.data(idx_kin_if + 245);

    auto tk_yyyzzz_yyy = pbuffer.data(idx_kin_if + 246);

    auto tk_yyyzzz_yyz = pbuffer.data(idx_kin_if + 247);

    auto tk_yyyzzz_yzz = pbuffer.data(idx_kin_if + 248);

    auto tk_yyyzzz_zzz = pbuffer.data(idx_kin_if + 249);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yyyz_xxy,   \
                             tk_yyyz_xyy,   \
                             tk_yyyz_yyy,   \
                             tk_yyyzz_xxy,  \
                             tk_yyyzz_xyy,  \
                             tk_yyyzz_yyy,  \
                             tk_yyyzzz_xxx, \
                             tk_yyyzzz_xxy, \
                             tk_yyyzzz_xxz, \
                             tk_yyyzzz_xyy, \
                             tk_yyyzzz_xyz, \
                             tk_yyyzzz_xzz, \
                             tk_yyyzzz_yyy, \
                             tk_yyyzzz_yyz, \
                             tk_yyyzzz_yzz, \
                             tk_yyyzzz_zzz, \
                             tk_yyzzz_xxx,  \
                             tk_yyzzz_xxz,  \
                             tk_yyzzz_xyz,  \
                             tk_yyzzz_xz,   \
                             tk_yyzzz_xzz,  \
                             tk_yyzzz_yyz,  \
                             tk_yyzzz_yz,   \
                             tk_yyzzz_yzz,  \
                             tk_yyzzz_zz,   \
                             tk_yyzzz_zzz,  \
                             tk_yzzz_xxx,   \
                             tk_yzzz_xxz,   \
                             tk_yzzz_xyz,   \
                             tk_yzzz_xzz,   \
                             tk_yzzz_yyz,   \
                             tk_yzzz_yzz,   \
                             tk_yzzz_zzz,   \
                             ts_yyyz_xxy,   \
                             ts_yyyz_xyy,   \
                             ts_yyyz_yyy,   \
                             ts_yyyzzz_xxx, \
                             ts_yyyzzz_xxy, \
                             ts_yyyzzz_xxz, \
                             ts_yyyzzz_xyy, \
                             ts_yyyzzz_xyz, \
                             ts_yyyzzz_xzz, \
                             ts_yyyzzz_yyy, \
                             ts_yyyzzz_yyz, \
                             ts_yyyzzz_yzz, \
                             ts_yyyzzz_zzz, \
                             ts_yzzz_xxx,   \
                             ts_yzzz_xxz,   \
                             ts_yzzz_xyz,   \
                             ts_yzzz_xzz,   \
                             ts_yzzz_yyz,   \
                             ts_yzzz_yzz,   \
                             ts_yzzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_xxx[i] =
            -4.0 * ts_yzzz_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxx[i] * fe_0 + tk_yyzzz_xxx[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxx[i] * fz_0;

        tk_yyyzzz_xxy[i] =
            -4.0 * ts_yyyz_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxy[i] * fe_0 + tk_yyyzz_xxy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxy[i] * fz_0;

        tk_yyyzzz_xxz[i] =
            -4.0 * ts_yzzz_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxz[i] * fe_0 + tk_yyzzz_xxz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxz[i] * fz_0;

        tk_yyyzzz_xyy[i] =
            -4.0 * ts_yyyz_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xyy[i] * fe_0 + tk_yyyzz_xyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xyy[i] * fz_0;

        tk_yyyzzz_xyz[i] = -4.0 * ts_yzzz_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyz[i] * fe_0 + tk_yyzzz_xz[i] * fe_0 + tk_yyzzz_xyz[i] * pa_y[i] +
                           2.0 * ts_yyyzzz_xyz[i] * fz_0;

        tk_yyyzzz_xzz[i] =
            -4.0 * ts_yzzz_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xzz[i] * fe_0 + tk_yyzzz_xzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xzz[i] * fz_0;

        tk_yyyzzz_yyy[i] =
            -4.0 * ts_yyyz_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_yyy[i] * fe_0 + tk_yyyzz_yyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_yyy[i] * fz_0;

        tk_yyyzzz_yyz[i] = -4.0 * ts_yzzz_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyz[i] * fe_0 + 2.0 * tk_yyzzz_yz[i] * fe_0 +
                           tk_yyzzz_yyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyz[i] * fz_0;

        tk_yyyzzz_yzz[i] = -4.0 * ts_yzzz_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yzz[i] * fe_0 + tk_yyzzz_zz[i] * fe_0 + tk_yyzzz_yzz[i] * pa_y[i] +
                           2.0 * ts_yyyzzz_yzz[i] * fz_0;

        tk_yyyzzz_zzz[i] =
            -4.0 * ts_yzzz_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_zzz[i] * fe_0 + tk_yyzzz_zzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_zzz[i] * fz_0;
    }

    // Set up 250-260 components of targeted buffer : IF

    auto tk_yyzzzz_xxx = pbuffer.data(idx_kin_if + 250);

    auto tk_yyzzzz_xxy = pbuffer.data(idx_kin_if + 251);

    auto tk_yyzzzz_xxz = pbuffer.data(idx_kin_if + 252);

    auto tk_yyzzzz_xyy = pbuffer.data(idx_kin_if + 253);

    auto tk_yyzzzz_xyz = pbuffer.data(idx_kin_if + 254);

    auto tk_yyzzzz_xzz = pbuffer.data(idx_kin_if + 255);

    auto tk_yyzzzz_yyy = pbuffer.data(idx_kin_if + 256);

    auto tk_yyzzzz_yyz = pbuffer.data(idx_kin_if + 257);

    auto tk_yyzzzz_yzz = pbuffer.data(idx_kin_if + 258);

    auto tk_yyzzzz_zzz = pbuffer.data(idx_kin_if + 259);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yyzz_xxy,   \
                             tk_yyzz_xyy,   \
                             tk_yyzz_yyy,   \
                             tk_yyzzz_xxy,  \
                             tk_yyzzz_xyy,  \
                             tk_yyzzz_yyy,  \
                             tk_yyzzzz_xxx, \
                             tk_yyzzzz_xxy, \
                             tk_yyzzzz_xxz, \
                             tk_yyzzzz_xyy, \
                             tk_yyzzzz_xyz, \
                             tk_yyzzzz_xzz, \
                             tk_yyzzzz_yyy, \
                             tk_yyzzzz_yyz, \
                             tk_yyzzzz_yzz, \
                             tk_yyzzzz_zzz, \
                             tk_yzzzz_xxx,  \
                             tk_yzzzz_xxz,  \
                             tk_yzzzz_xyz,  \
                             tk_yzzzz_xz,   \
                             tk_yzzzz_xzz,  \
                             tk_yzzzz_yyz,  \
                             tk_yzzzz_yz,   \
                             tk_yzzzz_yzz,  \
                             tk_yzzzz_zz,   \
                             tk_yzzzz_zzz,  \
                             tk_zzzz_xxx,   \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xzz,   \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yzz,   \
                             tk_zzzz_zzz,   \
                             ts_yyzz_xxy,   \
                             ts_yyzz_xyy,   \
                             ts_yyzz_yyy,   \
                             ts_yyzzzz_xxx, \
                             ts_yyzzzz_xxy, \
                             ts_yyzzzz_xxz, \
                             ts_yyzzzz_xyy, \
                             ts_yyzzzz_xyz, \
                             ts_yyzzzz_xzz, \
                             ts_yyzzzz_yyy, \
                             ts_yyzzzz_yyz, \
                             ts_yyzzzz_yzz, \
                             ts_yyzzzz_zzz, \
                             ts_zzzz_xxx,   \
                             ts_zzzz_xxz,   \
                             ts_zzzz_xyz,   \
                             ts_zzzz_xzz,   \
                             ts_zzzz_yyz,   \
                             ts_zzzz_yzz,   \
                             ts_zzzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_xxx[i] = -2.0 * ts_zzzz_xxx[i] * fbe_0 * fz_0 + tk_zzzz_xxx[i] * fe_0 + tk_yzzzz_xxx[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxx[i] * fz_0;

        tk_yyzzzz_xxy[i] =
            -6.0 * ts_yyzz_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxy[i] * fe_0 + tk_yyzzz_xxy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxy[i] * fz_0;

        tk_yyzzzz_xxz[i] = -2.0 * ts_zzzz_xxz[i] * fbe_0 * fz_0 + tk_zzzz_xxz[i] * fe_0 + tk_yzzzz_xxz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxz[i] * fz_0;

        tk_yyzzzz_xyy[i] =
            -6.0 * ts_yyzz_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyy[i] * fe_0 + tk_yyzzz_xyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xyy[i] * fz_0;

        tk_yyzzzz_xyz[i] = -2.0 * ts_zzzz_xyz[i] * fbe_0 * fz_0 + tk_zzzz_xyz[i] * fe_0 + tk_yzzzz_xz[i] * fe_0 + tk_yzzzz_xyz[i] * pa_y[i] +
                           2.0 * ts_yyzzzz_xyz[i] * fz_0;

        tk_yyzzzz_xzz[i] = -2.0 * ts_zzzz_xzz[i] * fbe_0 * fz_0 + tk_zzzz_xzz[i] * fe_0 + tk_yzzzz_xzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xzz[i] * fz_0;

        tk_yyzzzz_yyy[i] =
            -6.0 * ts_yyzz_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyy[i] * fe_0 + tk_yyzzz_yyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_yyy[i] * fz_0;

        tk_yyzzzz_yyz[i] = -2.0 * ts_zzzz_yyz[i] * fbe_0 * fz_0 + tk_zzzz_yyz[i] * fe_0 + 2.0 * tk_yzzzz_yz[i] * fe_0 + tk_yzzzz_yyz[i] * pa_y[i] +
                           2.0 * ts_yyzzzz_yyz[i] * fz_0;

        tk_yyzzzz_yzz[i] = -2.0 * ts_zzzz_yzz[i] * fbe_0 * fz_0 + tk_zzzz_yzz[i] * fe_0 + tk_yzzzz_zz[i] * fe_0 + tk_yzzzz_yzz[i] * pa_y[i] +
                           2.0 * ts_yyzzzz_yzz[i] * fz_0;

        tk_yyzzzz_zzz[i] = -2.0 * ts_zzzz_zzz[i] * fbe_0 * fz_0 + tk_zzzz_zzz[i] * fe_0 + tk_yzzzz_zzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_zzz[i] * fz_0;
    }

    // Set up 260-270 components of targeted buffer : IF

    auto tk_yzzzzz_xxx = pbuffer.data(idx_kin_if + 260);

    auto tk_yzzzzz_xxy = pbuffer.data(idx_kin_if + 261);

    auto tk_yzzzzz_xxz = pbuffer.data(idx_kin_if + 262);

    auto tk_yzzzzz_xyy = pbuffer.data(idx_kin_if + 263);

    auto tk_yzzzzz_xyz = pbuffer.data(idx_kin_if + 264);

    auto tk_yzzzzz_xzz = pbuffer.data(idx_kin_if + 265);

    auto tk_yzzzzz_yyy = pbuffer.data(idx_kin_if + 266);

    auto tk_yzzzzz_yyz = pbuffer.data(idx_kin_if + 267);

    auto tk_yzzzzz_yzz = pbuffer.data(idx_kin_if + 268);

    auto tk_yzzzzz_zzz = pbuffer.data(idx_kin_if + 269);

#pragma omp simd aligned(pa_y,              \
                             tk_yzzzzz_xxx, \
                             tk_yzzzzz_xxy, \
                             tk_yzzzzz_xxz, \
                             tk_yzzzzz_xyy, \
                             tk_yzzzzz_xyz, \
                             tk_yzzzzz_xzz, \
                             tk_yzzzzz_yyy, \
                             tk_yzzzzz_yyz, \
                             tk_yzzzzz_yzz, \
                             tk_yzzzzz_zzz, \
                             tk_zzzzz_xx,   \
                             tk_zzzzz_xxx,  \
                             tk_zzzzz_xxy,  \
                             tk_zzzzz_xxz,  \
                             tk_zzzzz_xy,   \
                             tk_zzzzz_xyy,  \
                             tk_zzzzz_xyz,  \
                             tk_zzzzz_xz,   \
                             tk_zzzzz_xzz,  \
                             tk_zzzzz_yy,   \
                             tk_zzzzz_yyy,  \
                             tk_zzzzz_yyz,  \
                             tk_zzzzz_yz,   \
                             tk_zzzzz_yzz,  \
                             tk_zzzzz_zz,   \
                             tk_zzzzz_zzz,  \
                             ts_yzzzzz_xxx, \
                             ts_yzzzzz_xxy, \
                             ts_yzzzzz_xxz, \
                             ts_yzzzzz_xyy, \
                             ts_yzzzzz_xyz, \
                             ts_yzzzzz_xzz, \
                             ts_yzzzzz_yyy, \
                             ts_yzzzzz_yyz, \
                             ts_yzzzzz_yzz, \
                             ts_yzzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_xxx[i] = tk_zzzzz_xxx[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxx[i] * fz_0;

        tk_yzzzzz_xxy[i] = tk_zzzzz_xx[i] * fe_0 + tk_zzzzz_xxy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxy[i] * fz_0;

        tk_yzzzzz_xxz[i] = tk_zzzzz_xxz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxz[i] * fz_0;

        tk_yzzzzz_xyy[i] = 2.0 * tk_zzzzz_xy[i] * fe_0 + tk_zzzzz_xyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyy[i] * fz_0;

        tk_yzzzzz_xyz[i] = tk_zzzzz_xz[i] * fe_0 + tk_zzzzz_xyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyz[i] * fz_0;

        tk_yzzzzz_xzz[i] = tk_zzzzz_xzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xzz[i] * fz_0;

        tk_yzzzzz_yyy[i] = 3.0 * tk_zzzzz_yy[i] * fe_0 + tk_zzzzz_yyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyy[i] * fz_0;

        tk_yzzzzz_yyz[i] = 2.0 * tk_zzzzz_yz[i] * fe_0 + tk_zzzzz_yyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyz[i] * fz_0;

        tk_yzzzzz_yzz[i] = tk_zzzzz_zz[i] * fe_0 + tk_zzzzz_yzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yzz[i] * fz_0;

        tk_yzzzzz_zzz[i] = tk_zzzzz_zzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_zzz[i] * fz_0;
    }

    // Set up 270-280 components of targeted buffer : IF

    auto tk_zzzzzz_xxx = pbuffer.data(idx_kin_if + 270);

    auto tk_zzzzzz_xxy = pbuffer.data(idx_kin_if + 271);

    auto tk_zzzzzz_xxz = pbuffer.data(idx_kin_if + 272);

    auto tk_zzzzzz_xyy = pbuffer.data(idx_kin_if + 273);

    auto tk_zzzzzz_xyz = pbuffer.data(idx_kin_if + 274);

    auto tk_zzzzzz_xzz = pbuffer.data(idx_kin_if + 275);

    auto tk_zzzzzz_yyy = pbuffer.data(idx_kin_if + 276);

    auto tk_zzzzzz_yyz = pbuffer.data(idx_kin_if + 277);

    auto tk_zzzzzz_yzz = pbuffer.data(idx_kin_if + 278);

    auto tk_zzzzzz_zzz = pbuffer.data(idx_kin_if + 279);

#pragma omp simd aligned(pa_z,              \
                             tk_zzzz_xxx,   \
                             tk_zzzz_xxy,   \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xyy,   \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xzz,   \
                             tk_zzzz_yyy,   \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yzz,   \
                             tk_zzzz_zzz,   \
                             tk_zzzzz_xx,   \
                             tk_zzzzz_xxx,  \
                             tk_zzzzz_xxy,  \
                             tk_zzzzz_xxz,  \
                             tk_zzzzz_xy,   \
                             tk_zzzzz_xyy,  \
                             tk_zzzzz_xyz,  \
                             tk_zzzzz_xz,   \
                             tk_zzzzz_xzz,  \
                             tk_zzzzz_yy,   \
                             tk_zzzzz_yyy,  \
                             tk_zzzzz_yyz,  \
                             tk_zzzzz_yz,   \
                             tk_zzzzz_yzz,  \
                             tk_zzzzz_zz,   \
                             tk_zzzzz_zzz,  \
                             tk_zzzzzz_xxx, \
                             tk_zzzzzz_xxy, \
                             tk_zzzzzz_xxz, \
                             tk_zzzzzz_xyy, \
                             tk_zzzzzz_xyz, \
                             tk_zzzzzz_xzz, \
                             tk_zzzzzz_yyy, \
                             tk_zzzzzz_yyz, \
                             tk_zzzzzz_yzz, \
                             tk_zzzzzz_zzz, \
                             ts_zzzz_xxx,   \
                             ts_zzzz_xxy,   \
                             ts_zzzz_xxz,   \
                             ts_zzzz_xyy,   \
                             ts_zzzz_xyz,   \
                             ts_zzzz_xzz,   \
                             ts_zzzz_yyy,   \
                             ts_zzzz_yyz,   \
                             ts_zzzz_yzz,   \
                             ts_zzzz_zzz,   \
                             ts_zzzzzz_xxx, \
                             ts_zzzzzz_xxy, \
                             ts_zzzzzz_xxz, \
                             ts_zzzzzz_xyy, \
                             ts_zzzzzz_xyz, \
                             ts_zzzzzz_xzz, \
                             ts_zzzzzz_yyy, \
                             ts_zzzzzz_yyz, \
                             ts_zzzzzz_yzz, \
                             ts_zzzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_xxx[i] =
            -10.0 * ts_zzzz_xxx[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxx[i] * fe_0 + tk_zzzzz_xxx[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxx[i] * fz_0;

        tk_zzzzzz_xxy[i] =
            -10.0 * ts_zzzz_xxy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxy[i] * fe_0 + tk_zzzzz_xxy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxy[i] * fz_0;

        tk_zzzzzz_xxz[i] = -10.0 * ts_zzzz_xxz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxz[i] * fe_0 + tk_zzzzz_xx[i] * fe_0 + tk_zzzzz_xxz[i] * pa_z[i] +
                           2.0 * ts_zzzzzz_xxz[i] * fz_0;

        tk_zzzzzz_xyy[i] =
            -10.0 * ts_zzzz_xyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyy[i] * fe_0 + tk_zzzzz_xyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyy[i] * fz_0;

        tk_zzzzzz_xyz[i] = -10.0 * ts_zzzz_xyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyz[i] * fe_0 + tk_zzzzz_xy[i] * fe_0 + tk_zzzzz_xyz[i] * pa_z[i] +
                           2.0 * ts_zzzzzz_xyz[i] * fz_0;

        tk_zzzzzz_xzz[i] = -10.0 * ts_zzzz_xzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xzz[i] * fe_0 + 2.0 * tk_zzzzz_xz[i] * fe_0 +
                           tk_zzzzz_xzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xzz[i] * fz_0;

        tk_zzzzzz_yyy[i] =
            -10.0 * ts_zzzz_yyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyy[i] * fe_0 + tk_zzzzz_yyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyy[i] * fz_0;

        tk_zzzzzz_yyz[i] = -10.0 * ts_zzzz_yyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyz[i] * fe_0 + tk_zzzzz_yy[i] * fe_0 + tk_zzzzz_yyz[i] * pa_z[i] +
                           2.0 * ts_zzzzzz_yyz[i] * fz_0;

        tk_zzzzzz_yzz[i] = -10.0 * ts_zzzz_yzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yzz[i] * fe_0 + 2.0 * tk_zzzzz_yz[i] * fe_0 +
                           tk_zzzzz_yzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yzz[i] * fz_0;

        tk_zzzzzz_zzz[i] = -10.0 * ts_zzzz_zzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_zzz[i] * fe_0 + 3.0 * tk_zzzzz_zz[i] * fe_0 +
                           tk_zzzzz_zzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_zzz[i] * fz_0;
    }
}

}  // namespace kinrec
