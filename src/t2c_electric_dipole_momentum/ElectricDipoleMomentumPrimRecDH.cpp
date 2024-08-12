#include "ElectricDipoleMomentumPrimRecDH.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_dh(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_dh,
                                      const size_t              idx_dip_sh,
                                      const size_t              idx_dip_pg,
                                      const size_t              idx_ovl_ph,
                                      const size_t              idx_dip_ph,
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

    // Set up components of auxiliary buffer : SH

    auto tr_x_0_xxxxx = pbuffer.data(idx_dip_sh);

    auto tr_x_0_xxxxy = pbuffer.data(idx_dip_sh + 1);

    auto tr_x_0_xxxxz = pbuffer.data(idx_dip_sh + 2);

    auto tr_x_0_xxxyy = pbuffer.data(idx_dip_sh + 3);

    auto tr_x_0_xxxyz = pbuffer.data(idx_dip_sh + 4);

    auto tr_x_0_xxxzz = pbuffer.data(idx_dip_sh + 5);

    auto tr_x_0_xxyyy = pbuffer.data(idx_dip_sh + 6);

    auto tr_x_0_xxyyz = pbuffer.data(idx_dip_sh + 7);

    auto tr_x_0_xxyzz = pbuffer.data(idx_dip_sh + 8);

    auto tr_x_0_xxzzz = pbuffer.data(idx_dip_sh + 9);

    auto tr_x_0_xyyyy = pbuffer.data(idx_dip_sh + 10);

    auto tr_x_0_xyyyz = pbuffer.data(idx_dip_sh + 11);

    auto tr_x_0_xyyzz = pbuffer.data(idx_dip_sh + 12);

    auto tr_x_0_xyzzz = pbuffer.data(idx_dip_sh + 13);

    auto tr_x_0_xzzzz = pbuffer.data(idx_dip_sh + 14);

    auto tr_x_0_yyyyy = pbuffer.data(idx_dip_sh + 15);

    auto tr_x_0_yyyyz = pbuffer.data(idx_dip_sh + 16);

    auto tr_x_0_yyyzz = pbuffer.data(idx_dip_sh + 17);

    auto tr_x_0_yyzzz = pbuffer.data(idx_dip_sh + 18);

    auto tr_x_0_yzzzz = pbuffer.data(idx_dip_sh + 19);

    auto tr_x_0_zzzzz = pbuffer.data(idx_dip_sh + 20);

    auto tr_y_0_xxxxx = pbuffer.data(idx_dip_sh + 21);

    auto tr_y_0_xxxxy = pbuffer.data(idx_dip_sh + 22);

    auto tr_y_0_xxxxz = pbuffer.data(idx_dip_sh + 23);

    auto tr_y_0_xxxyy = pbuffer.data(idx_dip_sh + 24);

    auto tr_y_0_xxxyz = pbuffer.data(idx_dip_sh + 25);

    auto tr_y_0_xxxzz = pbuffer.data(idx_dip_sh + 26);

    auto tr_y_0_xxyyy = pbuffer.data(idx_dip_sh + 27);

    auto tr_y_0_xxyyz = pbuffer.data(idx_dip_sh + 28);

    auto tr_y_0_xxyzz = pbuffer.data(idx_dip_sh + 29);

    auto tr_y_0_xxzzz = pbuffer.data(idx_dip_sh + 30);

    auto tr_y_0_xyyyy = pbuffer.data(idx_dip_sh + 31);

    auto tr_y_0_xyyyz = pbuffer.data(idx_dip_sh + 32);

    auto tr_y_0_xyyzz = pbuffer.data(idx_dip_sh + 33);

    auto tr_y_0_xyzzz = pbuffer.data(idx_dip_sh + 34);

    auto tr_y_0_xzzzz = pbuffer.data(idx_dip_sh + 35);

    auto tr_y_0_yyyyy = pbuffer.data(idx_dip_sh + 36);

    auto tr_y_0_yyyyz = pbuffer.data(idx_dip_sh + 37);

    auto tr_y_0_yyyzz = pbuffer.data(idx_dip_sh + 38);

    auto tr_y_0_yyzzz = pbuffer.data(idx_dip_sh + 39);

    auto tr_y_0_yzzzz = pbuffer.data(idx_dip_sh + 40);

    auto tr_y_0_zzzzz = pbuffer.data(idx_dip_sh + 41);

    auto tr_z_0_xxxxx = pbuffer.data(idx_dip_sh + 42);

    auto tr_z_0_xxxxy = pbuffer.data(idx_dip_sh + 43);

    auto tr_z_0_xxxxz = pbuffer.data(idx_dip_sh + 44);

    auto tr_z_0_xxxyy = pbuffer.data(idx_dip_sh + 45);

    auto tr_z_0_xxxyz = pbuffer.data(idx_dip_sh + 46);

    auto tr_z_0_xxxzz = pbuffer.data(idx_dip_sh + 47);

    auto tr_z_0_xxyyy = pbuffer.data(idx_dip_sh + 48);

    auto tr_z_0_xxyyz = pbuffer.data(idx_dip_sh + 49);

    auto tr_z_0_xxyzz = pbuffer.data(idx_dip_sh + 50);

    auto tr_z_0_xxzzz = pbuffer.data(idx_dip_sh + 51);

    auto tr_z_0_xyyyy = pbuffer.data(idx_dip_sh + 52);

    auto tr_z_0_xyyyz = pbuffer.data(idx_dip_sh + 53);

    auto tr_z_0_xyyzz = pbuffer.data(idx_dip_sh + 54);

    auto tr_z_0_xyzzz = pbuffer.data(idx_dip_sh + 55);

    auto tr_z_0_xzzzz = pbuffer.data(idx_dip_sh + 56);

    auto tr_z_0_yyyyy = pbuffer.data(idx_dip_sh + 57);

    auto tr_z_0_yyyyz = pbuffer.data(idx_dip_sh + 58);

    auto tr_z_0_yyyzz = pbuffer.data(idx_dip_sh + 59);

    auto tr_z_0_yyzzz = pbuffer.data(idx_dip_sh + 60);

    auto tr_z_0_yzzzz = pbuffer.data(idx_dip_sh + 61);

    auto tr_z_0_zzzzz = pbuffer.data(idx_dip_sh + 62);

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

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ovl_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ovl_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ovl_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ovl_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ovl_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ovl_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ovl_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ovl_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ovl_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ovl_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ovl_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ovl_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ovl_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ovl_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ovl_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ovl_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ovl_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ovl_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ovl_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ovl_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ovl_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ovl_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ovl_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ovl_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ovl_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ovl_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ovl_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ovl_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ovl_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ovl_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ovl_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ovl_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ovl_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ovl_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ovl_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ovl_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ovl_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ovl_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ovl_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ovl_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ovl_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ovl_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ovl_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ovl_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ovl_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ovl_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ovl_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ovl_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ovl_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ovl_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ovl_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ovl_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ovl_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ovl_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ovl_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ovl_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ovl_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ovl_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ovl_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ovl_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ovl_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ovl_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ovl_ph + 62);

    // Set up components of auxiliary buffer : PH

    auto tr_x_x_xxxxx = pbuffer.data(idx_dip_ph);

    auto tr_x_x_xxxxy = pbuffer.data(idx_dip_ph + 1);

    auto tr_x_x_xxxxz = pbuffer.data(idx_dip_ph + 2);

    auto tr_x_x_xxxyy = pbuffer.data(idx_dip_ph + 3);

    auto tr_x_x_xxxyz = pbuffer.data(idx_dip_ph + 4);

    auto tr_x_x_xxxzz = pbuffer.data(idx_dip_ph + 5);

    auto tr_x_x_xxyyy = pbuffer.data(idx_dip_ph + 6);

    auto tr_x_x_xxyyz = pbuffer.data(idx_dip_ph + 7);

    auto tr_x_x_xxyzz = pbuffer.data(idx_dip_ph + 8);

    auto tr_x_x_xxzzz = pbuffer.data(idx_dip_ph + 9);

    auto tr_x_x_xyyyy = pbuffer.data(idx_dip_ph + 10);

    auto tr_x_x_xyyyz = pbuffer.data(idx_dip_ph + 11);

    auto tr_x_x_xyyzz = pbuffer.data(idx_dip_ph + 12);

    auto tr_x_x_xyzzz = pbuffer.data(idx_dip_ph + 13);

    auto tr_x_x_xzzzz = pbuffer.data(idx_dip_ph + 14);

    auto tr_x_x_yyyyy = pbuffer.data(idx_dip_ph + 15);

    auto tr_x_x_yyyyz = pbuffer.data(idx_dip_ph + 16);

    auto tr_x_x_yyyzz = pbuffer.data(idx_dip_ph + 17);

    auto tr_x_x_yyzzz = pbuffer.data(idx_dip_ph + 18);

    auto tr_x_x_yzzzz = pbuffer.data(idx_dip_ph + 19);

    auto tr_x_x_zzzzz = pbuffer.data(idx_dip_ph + 20);

    auto tr_x_y_xxxxx = pbuffer.data(idx_dip_ph + 21);

    auto tr_x_y_xxxxy = pbuffer.data(idx_dip_ph + 22);

    auto tr_x_y_xxxxz = pbuffer.data(idx_dip_ph + 23);

    auto tr_x_y_xxxyy = pbuffer.data(idx_dip_ph + 24);

    auto tr_x_y_xxxyz = pbuffer.data(idx_dip_ph + 25);

    auto tr_x_y_xxxzz = pbuffer.data(idx_dip_ph + 26);

    auto tr_x_y_xxyyy = pbuffer.data(idx_dip_ph + 27);

    auto tr_x_y_xxyyz = pbuffer.data(idx_dip_ph + 28);

    auto tr_x_y_xxyzz = pbuffer.data(idx_dip_ph + 29);

    auto tr_x_y_xxzzz = pbuffer.data(idx_dip_ph + 30);

    auto tr_x_y_xyyyy = pbuffer.data(idx_dip_ph + 31);

    auto tr_x_y_xyyyz = pbuffer.data(idx_dip_ph + 32);

    auto tr_x_y_xyyzz = pbuffer.data(idx_dip_ph + 33);

    auto tr_x_y_xyzzz = pbuffer.data(idx_dip_ph + 34);

    auto tr_x_y_xzzzz = pbuffer.data(idx_dip_ph + 35);

    auto tr_x_y_yyyyy = pbuffer.data(idx_dip_ph + 36);

    auto tr_x_y_yyyyz = pbuffer.data(idx_dip_ph + 37);

    auto tr_x_y_yyyzz = pbuffer.data(idx_dip_ph + 38);

    auto tr_x_y_yyzzz = pbuffer.data(idx_dip_ph + 39);

    auto tr_x_y_yzzzz = pbuffer.data(idx_dip_ph + 40);

    auto tr_x_y_zzzzz = pbuffer.data(idx_dip_ph + 41);

    auto tr_x_z_xxxxx = pbuffer.data(idx_dip_ph + 42);

    auto tr_x_z_xxxxy = pbuffer.data(idx_dip_ph + 43);

    auto tr_x_z_xxxxz = pbuffer.data(idx_dip_ph + 44);

    auto tr_x_z_xxxyy = pbuffer.data(idx_dip_ph + 45);

    auto tr_x_z_xxxyz = pbuffer.data(idx_dip_ph + 46);

    auto tr_x_z_xxxzz = pbuffer.data(idx_dip_ph + 47);

    auto tr_x_z_xxyyy = pbuffer.data(idx_dip_ph + 48);

    auto tr_x_z_xxyyz = pbuffer.data(idx_dip_ph + 49);

    auto tr_x_z_xxyzz = pbuffer.data(idx_dip_ph + 50);

    auto tr_x_z_xxzzz = pbuffer.data(idx_dip_ph + 51);

    auto tr_x_z_xyyyy = pbuffer.data(idx_dip_ph + 52);

    auto tr_x_z_xyyyz = pbuffer.data(idx_dip_ph + 53);

    auto tr_x_z_xyyzz = pbuffer.data(idx_dip_ph + 54);

    auto tr_x_z_xyzzz = pbuffer.data(idx_dip_ph + 55);

    auto tr_x_z_xzzzz = pbuffer.data(idx_dip_ph + 56);

    auto tr_x_z_yyyyy = pbuffer.data(idx_dip_ph + 57);

    auto tr_x_z_yyyyz = pbuffer.data(idx_dip_ph + 58);

    auto tr_x_z_yyyzz = pbuffer.data(idx_dip_ph + 59);

    auto tr_x_z_yyzzz = pbuffer.data(idx_dip_ph + 60);

    auto tr_x_z_yzzzz = pbuffer.data(idx_dip_ph + 61);

    auto tr_x_z_zzzzz = pbuffer.data(idx_dip_ph + 62);

    auto tr_y_x_xxxxx = pbuffer.data(idx_dip_ph + 63);

    auto tr_y_x_xxxxy = pbuffer.data(idx_dip_ph + 64);

    auto tr_y_x_xxxxz = pbuffer.data(idx_dip_ph + 65);

    auto tr_y_x_xxxyy = pbuffer.data(idx_dip_ph + 66);

    auto tr_y_x_xxxyz = pbuffer.data(idx_dip_ph + 67);

    auto tr_y_x_xxxzz = pbuffer.data(idx_dip_ph + 68);

    auto tr_y_x_xxyyy = pbuffer.data(idx_dip_ph + 69);

    auto tr_y_x_xxyyz = pbuffer.data(idx_dip_ph + 70);

    auto tr_y_x_xxyzz = pbuffer.data(idx_dip_ph + 71);

    auto tr_y_x_xxzzz = pbuffer.data(idx_dip_ph + 72);

    auto tr_y_x_xyyyy = pbuffer.data(idx_dip_ph + 73);

    auto tr_y_x_xyyyz = pbuffer.data(idx_dip_ph + 74);

    auto tr_y_x_xyyzz = pbuffer.data(idx_dip_ph + 75);

    auto tr_y_x_xyzzz = pbuffer.data(idx_dip_ph + 76);

    auto tr_y_x_xzzzz = pbuffer.data(idx_dip_ph + 77);

    auto tr_y_x_yyyyy = pbuffer.data(idx_dip_ph + 78);

    auto tr_y_x_yyyyz = pbuffer.data(idx_dip_ph + 79);

    auto tr_y_x_yyyzz = pbuffer.data(idx_dip_ph + 80);

    auto tr_y_x_yyzzz = pbuffer.data(idx_dip_ph + 81);

    auto tr_y_x_yzzzz = pbuffer.data(idx_dip_ph + 82);

    auto tr_y_x_zzzzz = pbuffer.data(idx_dip_ph + 83);

    auto tr_y_y_xxxxx = pbuffer.data(idx_dip_ph + 84);

    auto tr_y_y_xxxxy = pbuffer.data(idx_dip_ph + 85);

    auto tr_y_y_xxxxz = pbuffer.data(idx_dip_ph + 86);

    auto tr_y_y_xxxyy = pbuffer.data(idx_dip_ph + 87);

    auto tr_y_y_xxxyz = pbuffer.data(idx_dip_ph + 88);

    auto tr_y_y_xxxzz = pbuffer.data(idx_dip_ph + 89);

    auto tr_y_y_xxyyy = pbuffer.data(idx_dip_ph + 90);

    auto tr_y_y_xxyyz = pbuffer.data(idx_dip_ph + 91);

    auto tr_y_y_xxyzz = pbuffer.data(idx_dip_ph + 92);

    auto tr_y_y_xxzzz = pbuffer.data(idx_dip_ph + 93);

    auto tr_y_y_xyyyy = pbuffer.data(idx_dip_ph + 94);

    auto tr_y_y_xyyyz = pbuffer.data(idx_dip_ph + 95);

    auto tr_y_y_xyyzz = pbuffer.data(idx_dip_ph + 96);

    auto tr_y_y_xyzzz = pbuffer.data(idx_dip_ph + 97);

    auto tr_y_y_xzzzz = pbuffer.data(idx_dip_ph + 98);

    auto tr_y_y_yyyyy = pbuffer.data(idx_dip_ph + 99);

    auto tr_y_y_yyyyz = pbuffer.data(idx_dip_ph + 100);

    auto tr_y_y_yyyzz = pbuffer.data(idx_dip_ph + 101);

    auto tr_y_y_yyzzz = pbuffer.data(idx_dip_ph + 102);

    auto tr_y_y_yzzzz = pbuffer.data(idx_dip_ph + 103);

    auto tr_y_y_zzzzz = pbuffer.data(idx_dip_ph + 104);

    auto tr_y_z_xxxxx = pbuffer.data(idx_dip_ph + 105);

    auto tr_y_z_xxxxy = pbuffer.data(idx_dip_ph + 106);

    auto tr_y_z_xxxxz = pbuffer.data(idx_dip_ph + 107);

    auto tr_y_z_xxxyy = pbuffer.data(idx_dip_ph + 108);

    auto tr_y_z_xxxyz = pbuffer.data(idx_dip_ph + 109);

    auto tr_y_z_xxxzz = pbuffer.data(idx_dip_ph + 110);

    auto tr_y_z_xxyyy = pbuffer.data(idx_dip_ph + 111);

    auto tr_y_z_xxyyz = pbuffer.data(idx_dip_ph + 112);

    auto tr_y_z_xxyzz = pbuffer.data(idx_dip_ph + 113);

    auto tr_y_z_xxzzz = pbuffer.data(idx_dip_ph + 114);

    auto tr_y_z_xyyyy = pbuffer.data(idx_dip_ph + 115);

    auto tr_y_z_xyyyz = pbuffer.data(idx_dip_ph + 116);

    auto tr_y_z_xyyzz = pbuffer.data(idx_dip_ph + 117);

    auto tr_y_z_xyzzz = pbuffer.data(idx_dip_ph + 118);

    auto tr_y_z_xzzzz = pbuffer.data(idx_dip_ph + 119);

    auto tr_y_z_yyyyy = pbuffer.data(idx_dip_ph + 120);

    auto tr_y_z_yyyyz = pbuffer.data(idx_dip_ph + 121);

    auto tr_y_z_yyyzz = pbuffer.data(idx_dip_ph + 122);

    auto tr_y_z_yyzzz = pbuffer.data(idx_dip_ph + 123);

    auto tr_y_z_yzzzz = pbuffer.data(idx_dip_ph + 124);

    auto tr_y_z_zzzzz = pbuffer.data(idx_dip_ph + 125);

    auto tr_z_x_xxxxx = pbuffer.data(idx_dip_ph + 126);

    auto tr_z_x_xxxxy = pbuffer.data(idx_dip_ph + 127);

    auto tr_z_x_xxxxz = pbuffer.data(idx_dip_ph + 128);

    auto tr_z_x_xxxyy = pbuffer.data(idx_dip_ph + 129);

    auto tr_z_x_xxxyz = pbuffer.data(idx_dip_ph + 130);

    auto tr_z_x_xxxzz = pbuffer.data(idx_dip_ph + 131);

    auto tr_z_x_xxyyy = pbuffer.data(idx_dip_ph + 132);

    auto tr_z_x_xxyyz = pbuffer.data(idx_dip_ph + 133);

    auto tr_z_x_xxyzz = pbuffer.data(idx_dip_ph + 134);

    auto tr_z_x_xxzzz = pbuffer.data(idx_dip_ph + 135);

    auto tr_z_x_xyyyy = pbuffer.data(idx_dip_ph + 136);

    auto tr_z_x_xyyyz = pbuffer.data(idx_dip_ph + 137);

    auto tr_z_x_xyyzz = pbuffer.data(idx_dip_ph + 138);

    auto tr_z_x_xyzzz = pbuffer.data(idx_dip_ph + 139);

    auto tr_z_x_xzzzz = pbuffer.data(idx_dip_ph + 140);

    auto tr_z_x_yyyyy = pbuffer.data(idx_dip_ph + 141);

    auto tr_z_x_yyyyz = pbuffer.data(idx_dip_ph + 142);

    auto tr_z_x_yyyzz = pbuffer.data(idx_dip_ph + 143);

    auto tr_z_x_yyzzz = pbuffer.data(idx_dip_ph + 144);

    auto tr_z_x_yzzzz = pbuffer.data(idx_dip_ph + 145);

    auto tr_z_x_zzzzz = pbuffer.data(idx_dip_ph + 146);

    auto tr_z_y_xxxxx = pbuffer.data(idx_dip_ph + 147);

    auto tr_z_y_xxxxy = pbuffer.data(idx_dip_ph + 148);

    auto tr_z_y_xxxxz = pbuffer.data(idx_dip_ph + 149);

    auto tr_z_y_xxxyy = pbuffer.data(idx_dip_ph + 150);

    auto tr_z_y_xxxyz = pbuffer.data(idx_dip_ph + 151);

    auto tr_z_y_xxxzz = pbuffer.data(idx_dip_ph + 152);

    auto tr_z_y_xxyyy = pbuffer.data(idx_dip_ph + 153);

    auto tr_z_y_xxyyz = pbuffer.data(idx_dip_ph + 154);

    auto tr_z_y_xxyzz = pbuffer.data(idx_dip_ph + 155);

    auto tr_z_y_xxzzz = pbuffer.data(idx_dip_ph + 156);

    auto tr_z_y_xyyyy = pbuffer.data(idx_dip_ph + 157);

    auto tr_z_y_xyyyz = pbuffer.data(idx_dip_ph + 158);

    auto tr_z_y_xyyzz = pbuffer.data(idx_dip_ph + 159);

    auto tr_z_y_xyzzz = pbuffer.data(idx_dip_ph + 160);

    auto tr_z_y_xzzzz = pbuffer.data(idx_dip_ph + 161);

    auto tr_z_y_yyyyy = pbuffer.data(idx_dip_ph + 162);

    auto tr_z_y_yyyyz = pbuffer.data(idx_dip_ph + 163);

    auto tr_z_y_yyyzz = pbuffer.data(idx_dip_ph + 164);

    auto tr_z_y_yyzzz = pbuffer.data(idx_dip_ph + 165);

    auto tr_z_y_yzzzz = pbuffer.data(idx_dip_ph + 166);

    auto tr_z_y_zzzzz = pbuffer.data(idx_dip_ph + 167);

    auto tr_z_z_xxxxx = pbuffer.data(idx_dip_ph + 168);

    auto tr_z_z_xxxxy = pbuffer.data(idx_dip_ph + 169);

    auto tr_z_z_xxxxz = pbuffer.data(idx_dip_ph + 170);

    auto tr_z_z_xxxyy = pbuffer.data(idx_dip_ph + 171);

    auto tr_z_z_xxxyz = pbuffer.data(idx_dip_ph + 172);

    auto tr_z_z_xxxzz = pbuffer.data(idx_dip_ph + 173);

    auto tr_z_z_xxyyy = pbuffer.data(idx_dip_ph + 174);

    auto tr_z_z_xxyyz = pbuffer.data(idx_dip_ph + 175);

    auto tr_z_z_xxyzz = pbuffer.data(idx_dip_ph + 176);

    auto tr_z_z_xxzzz = pbuffer.data(idx_dip_ph + 177);

    auto tr_z_z_xyyyy = pbuffer.data(idx_dip_ph + 178);

    auto tr_z_z_xyyyz = pbuffer.data(idx_dip_ph + 179);

    auto tr_z_z_xyyzz = pbuffer.data(idx_dip_ph + 180);

    auto tr_z_z_xyzzz = pbuffer.data(idx_dip_ph + 181);

    auto tr_z_z_xzzzz = pbuffer.data(idx_dip_ph + 182);

    auto tr_z_z_yyyyy = pbuffer.data(idx_dip_ph + 183);

    auto tr_z_z_yyyyz = pbuffer.data(idx_dip_ph + 184);

    auto tr_z_z_yyyzz = pbuffer.data(idx_dip_ph + 185);

    auto tr_z_z_yyzzz = pbuffer.data(idx_dip_ph + 186);

    auto tr_z_z_yzzzz = pbuffer.data(idx_dip_ph + 187);

    auto tr_z_z_zzzzz = pbuffer.data(idx_dip_ph + 188);

    // Set up 0-21 components of targeted buffer : DH

    auto tr_x_xx_xxxxx = pbuffer.data(idx_dip_dh);

    auto tr_x_xx_xxxxy = pbuffer.data(idx_dip_dh + 1);

    auto tr_x_xx_xxxxz = pbuffer.data(idx_dip_dh + 2);

    auto tr_x_xx_xxxyy = pbuffer.data(idx_dip_dh + 3);

    auto tr_x_xx_xxxyz = pbuffer.data(idx_dip_dh + 4);

    auto tr_x_xx_xxxzz = pbuffer.data(idx_dip_dh + 5);

    auto tr_x_xx_xxyyy = pbuffer.data(idx_dip_dh + 6);

    auto tr_x_xx_xxyyz = pbuffer.data(idx_dip_dh + 7);

    auto tr_x_xx_xxyzz = pbuffer.data(idx_dip_dh + 8);

    auto tr_x_xx_xxzzz = pbuffer.data(idx_dip_dh + 9);

    auto tr_x_xx_xyyyy = pbuffer.data(idx_dip_dh + 10);

    auto tr_x_xx_xyyyz = pbuffer.data(idx_dip_dh + 11);

    auto tr_x_xx_xyyzz = pbuffer.data(idx_dip_dh + 12);

    auto tr_x_xx_xyzzz = pbuffer.data(idx_dip_dh + 13);

    auto tr_x_xx_xzzzz = pbuffer.data(idx_dip_dh + 14);

    auto tr_x_xx_yyyyy = pbuffer.data(idx_dip_dh + 15);

    auto tr_x_xx_yyyyz = pbuffer.data(idx_dip_dh + 16);

    auto tr_x_xx_yyyzz = pbuffer.data(idx_dip_dh + 17);

    auto tr_x_xx_yyzzz = pbuffer.data(idx_dip_dh + 18);

    auto tr_x_xx_yzzzz = pbuffer.data(idx_dip_dh + 19);

    auto tr_x_xx_zzzzz = pbuffer.data(idx_dip_dh + 20);

#pragma omp simd aligned(pa_x,              \
                             tr_x_0_xxxxx,  \
                             tr_x_0_xxxxy,  \
                             tr_x_0_xxxxz,  \
                             tr_x_0_xxxyy,  \
                             tr_x_0_xxxyz,  \
                             tr_x_0_xxxzz,  \
                             tr_x_0_xxyyy,  \
                             tr_x_0_xxyyz,  \
                             tr_x_0_xxyzz,  \
                             tr_x_0_xxzzz,  \
                             tr_x_0_xyyyy,  \
                             tr_x_0_xyyyz,  \
                             tr_x_0_xyyzz,  \
                             tr_x_0_xyzzz,  \
                             tr_x_0_xzzzz,  \
                             tr_x_0_yyyyy,  \
                             tr_x_0_yyyyz,  \
                             tr_x_0_yyyzz,  \
                             tr_x_0_yyzzz,  \
                             tr_x_0_yzzzz,  \
                             tr_x_0_zzzzz,  \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxxx,  \
                             tr_x_x_xxxxy,  \
                             tr_x_x_xxxxz,  \
                             tr_x_x_xxxy,   \
                             tr_x_x_xxxyy,  \
                             tr_x_x_xxxyz,  \
                             tr_x_x_xxxz,   \
                             tr_x_x_xxxzz,  \
                             tr_x_x_xxyy,   \
                             tr_x_x_xxyyy,  \
                             tr_x_x_xxyyz,  \
                             tr_x_x_xxyz,   \
                             tr_x_x_xxyzz,  \
                             tr_x_x_xxzz,   \
                             tr_x_x_xxzzz,  \
                             tr_x_x_xyyy,   \
                             tr_x_x_xyyyy,  \
                             tr_x_x_xyyyz,  \
                             tr_x_x_xyyz,   \
                             tr_x_x_xyyzz,  \
                             tr_x_x_xyzz,   \
                             tr_x_x_xyzzz,  \
                             tr_x_x_xzzz,   \
                             tr_x_x_xzzzz,  \
                             tr_x_x_yyyy,   \
                             tr_x_x_yyyyy,  \
                             tr_x_x_yyyyz,  \
                             tr_x_x_yyyz,   \
                             tr_x_x_yyyzz,  \
                             tr_x_x_yyzz,   \
                             tr_x_x_yyzzz,  \
                             tr_x_x_yzzz,   \
                             tr_x_x_yzzzz,  \
                             tr_x_x_zzzz,   \
                             tr_x_x_zzzzz,  \
                             tr_x_xx_xxxxx, \
                             tr_x_xx_xxxxy, \
                             tr_x_xx_xxxxz, \
                             tr_x_xx_xxxyy, \
                             tr_x_xx_xxxyz, \
                             tr_x_xx_xxxzz, \
                             tr_x_xx_xxyyy, \
                             tr_x_xx_xxyyz, \
                             tr_x_xx_xxyzz, \
                             tr_x_xx_xxzzz, \
                             tr_x_xx_xyyyy, \
                             tr_x_xx_xyyyz, \
                             tr_x_xx_xyyzz, \
                             tr_x_xx_xyzzz, \
                             tr_x_xx_xzzzz, \
                             tr_x_xx_yyyyy, \
                             tr_x_xx_yyyyz, \
                             tr_x_xx_yyyzz, \
                             tr_x_xx_yyzzz, \
                             tr_x_xx_yzzzz, \
                             tr_x_xx_zzzzz, \
                             ts_x_xxxxx,    \
                             ts_x_xxxxy,    \
                             ts_x_xxxxz,    \
                             ts_x_xxxyy,    \
                             ts_x_xxxyz,    \
                             ts_x_xxxzz,    \
                             ts_x_xxyyy,    \
                             ts_x_xxyyz,    \
                             ts_x_xxyzz,    \
                             ts_x_xxzzz,    \
                             ts_x_xyyyy,    \
                             ts_x_xyyyz,    \
                             ts_x_xyyzz,    \
                             ts_x_xyzzz,    \
                             ts_x_xzzzz,    \
                             ts_x_yyyyy,    \
                             ts_x_yyyyz,    \
                             ts_x_yyyzz,    \
                             ts_x_yyzzz,    \
                             ts_x_yzzzz,    \
                             ts_x_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_xxxxx[i] = tr_x_0_xxxxx[i] * fe_0 + 5.0 * tr_x_x_xxxx[i] * fe_0 + ts_x_xxxxx[i] * fe_0 + tr_x_x_xxxxx[i] * pa_x[i];

        tr_x_xx_xxxxy[i] = tr_x_0_xxxxy[i] * fe_0 + 4.0 * tr_x_x_xxxy[i] * fe_0 + ts_x_xxxxy[i] * fe_0 + tr_x_x_xxxxy[i] * pa_x[i];

        tr_x_xx_xxxxz[i] = tr_x_0_xxxxz[i] * fe_0 + 4.0 * tr_x_x_xxxz[i] * fe_0 + ts_x_xxxxz[i] * fe_0 + tr_x_x_xxxxz[i] * pa_x[i];

        tr_x_xx_xxxyy[i] = tr_x_0_xxxyy[i] * fe_0 + 3.0 * tr_x_x_xxyy[i] * fe_0 + ts_x_xxxyy[i] * fe_0 + tr_x_x_xxxyy[i] * pa_x[i];

        tr_x_xx_xxxyz[i] = tr_x_0_xxxyz[i] * fe_0 + 3.0 * tr_x_x_xxyz[i] * fe_0 + ts_x_xxxyz[i] * fe_0 + tr_x_x_xxxyz[i] * pa_x[i];

        tr_x_xx_xxxzz[i] = tr_x_0_xxxzz[i] * fe_0 + 3.0 * tr_x_x_xxzz[i] * fe_0 + ts_x_xxxzz[i] * fe_0 + tr_x_x_xxxzz[i] * pa_x[i];

        tr_x_xx_xxyyy[i] = tr_x_0_xxyyy[i] * fe_0 + 2.0 * tr_x_x_xyyy[i] * fe_0 + ts_x_xxyyy[i] * fe_0 + tr_x_x_xxyyy[i] * pa_x[i];

        tr_x_xx_xxyyz[i] = tr_x_0_xxyyz[i] * fe_0 + 2.0 * tr_x_x_xyyz[i] * fe_0 + ts_x_xxyyz[i] * fe_0 + tr_x_x_xxyyz[i] * pa_x[i];

        tr_x_xx_xxyzz[i] = tr_x_0_xxyzz[i] * fe_0 + 2.0 * tr_x_x_xyzz[i] * fe_0 + ts_x_xxyzz[i] * fe_0 + tr_x_x_xxyzz[i] * pa_x[i];

        tr_x_xx_xxzzz[i] = tr_x_0_xxzzz[i] * fe_0 + 2.0 * tr_x_x_xzzz[i] * fe_0 + ts_x_xxzzz[i] * fe_0 + tr_x_x_xxzzz[i] * pa_x[i];

        tr_x_xx_xyyyy[i] = tr_x_0_xyyyy[i] * fe_0 + tr_x_x_yyyy[i] * fe_0 + ts_x_xyyyy[i] * fe_0 + tr_x_x_xyyyy[i] * pa_x[i];

        tr_x_xx_xyyyz[i] = tr_x_0_xyyyz[i] * fe_0 + tr_x_x_yyyz[i] * fe_0 + ts_x_xyyyz[i] * fe_0 + tr_x_x_xyyyz[i] * pa_x[i];

        tr_x_xx_xyyzz[i] = tr_x_0_xyyzz[i] * fe_0 + tr_x_x_yyzz[i] * fe_0 + ts_x_xyyzz[i] * fe_0 + tr_x_x_xyyzz[i] * pa_x[i];

        tr_x_xx_xyzzz[i] = tr_x_0_xyzzz[i] * fe_0 + tr_x_x_yzzz[i] * fe_0 + ts_x_xyzzz[i] * fe_0 + tr_x_x_xyzzz[i] * pa_x[i];

        tr_x_xx_xzzzz[i] = tr_x_0_xzzzz[i] * fe_0 + tr_x_x_zzzz[i] * fe_0 + ts_x_xzzzz[i] * fe_0 + tr_x_x_xzzzz[i] * pa_x[i];

        tr_x_xx_yyyyy[i] = tr_x_0_yyyyy[i] * fe_0 + ts_x_yyyyy[i] * fe_0 + tr_x_x_yyyyy[i] * pa_x[i];

        tr_x_xx_yyyyz[i] = tr_x_0_yyyyz[i] * fe_0 + ts_x_yyyyz[i] * fe_0 + tr_x_x_yyyyz[i] * pa_x[i];

        tr_x_xx_yyyzz[i] = tr_x_0_yyyzz[i] * fe_0 + ts_x_yyyzz[i] * fe_0 + tr_x_x_yyyzz[i] * pa_x[i];

        tr_x_xx_yyzzz[i] = tr_x_0_yyzzz[i] * fe_0 + ts_x_yyzzz[i] * fe_0 + tr_x_x_yyzzz[i] * pa_x[i];

        tr_x_xx_yzzzz[i] = tr_x_0_yzzzz[i] * fe_0 + ts_x_yzzzz[i] * fe_0 + tr_x_x_yzzzz[i] * pa_x[i];

        tr_x_xx_zzzzz[i] = tr_x_0_zzzzz[i] * fe_0 + ts_x_zzzzz[i] * fe_0 + tr_x_x_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : DH

    auto tr_x_xy_xxxxx = pbuffer.data(idx_dip_dh + 21);

    auto tr_x_xy_xxxxy = pbuffer.data(idx_dip_dh + 22);

    auto tr_x_xy_xxxxz = pbuffer.data(idx_dip_dh + 23);

    auto tr_x_xy_xxxyy = pbuffer.data(idx_dip_dh + 24);

    auto tr_x_xy_xxxyz = pbuffer.data(idx_dip_dh + 25);

    auto tr_x_xy_xxxzz = pbuffer.data(idx_dip_dh + 26);

    auto tr_x_xy_xxyyy = pbuffer.data(idx_dip_dh + 27);

    auto tr_x_xy_xxyyz = pbuffer.data(idx_dip_dh + 28);

    auto tr_x_xy_xxyzz = pbuffer.data(idx_dip_dh + 29);

    auto tr_x_xy_xxzzz = pbuffer.data(idx_dip_dh + 30);

    auto tr_x_xy_xyyyy = pbuffer.data(idx_dip_dh + 31);

    auto tr_x_xy_xyyyz = pbuffer.data(idx_dip_dh + 32);

    auto tr_x_xy_xyyzz = pbuffer.data(idx_dip_dh + 33);

    auto tr_x_xy_xyzzz = pbuffer.data(idx_dip_dh + 34);

    auto tr_x_xy_xzzzz = pbuffer.data(idx_dip_dh + 35);

    auto tr_x_xy_yyyyy = pbuffer.data(idx_dip_dh + 36);

    auto tr_x_xy_yyyyz = pbuffer.data(idx_dip_dh + 37);

    auto tr_x_xy_yyyzz = pbuffer.data(idx_dip_dh + 38);

    auto tr_x_xy_yyzzz = pbuffer.data(idx_dip_dh + 39);

    auto tr_x_xy_yzzzz = pbuffer.data(idx_dip_dh + 40);

    auto tr_x_xy_zzzzz = pbuffer.data(idx_dip_dh + 41);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxxx,  \
                             tr_x_x_xxxxy,  \
                             tr_x_x_xxxxz,  \
                             tr_x_x_xxxy,   \
                             tr_x_x_xxxyy,  \
                             tr_x_x_xxxyz,  \
                             tr_x_x_xxxz,   \
                             tr_x_x_xxxzz,  \
                             tr_x_x_xxyy,   \
                             tr_x_x_xxyyy,  \
                             tr_x_x_xxyyz,  \
                             tr_x_x_xxyz,   \
                             tr_x_x_xxyzz,  \
                             tr_x_x_xxzz,   \
                             tr_x_x_xxzzz,  \
                             tr_x_x_xyyy,   \
                             tr_x_x_xyyyy,  \
                             tr_x_x_xyyyz,  \
                             tr_x_x_xyyz,   \
                             tr_x_x_xyyzz,  \
                             tr_x_x_xyzz,   \
                             tr_x_x_xyzzz,  \
                             tr_x_x_xzzz,   \
                             tr_x_x_xzzzz,  \
                             tr_x_x_zzzzz,  \
                             tr_x_xy_xxxxx, \
                             tr_x_xy_xxxxy, \
                             tr_x_xy_xxxxz, \
                             tr_x_xy_xxxyy, \
                             tr_x_xy_xxxyz, \
                             tr_x_xy_xxxzz, \
                             tr_x_xy_xxyyy, \
                             tr_x_xy_xxyyz, \
                             tr_x_xy_xxyzz, \
                             tr_x_xy_xxzzz, \
                             tr_x_xy_xyyyy, \
                             tr_x_xy_xyyyz, \
                             tr_x_xy_xyyzz, \
                             tr_x_xy_xyzzz, \
                             tr_x_xy_xzzzz, \
                             tr_x_xy_yyyyy, \
                             tr_x_xy_yyyyz, \
                             tr_x_xy_yyyzz, \
                             tr_x_xy_yyzzz, \
                             tr_x_xy_yzzzz, \
                             tr_x_xy_zzzzz, \
                             tr_x_y_yyyyy,  \
                             tr_x_y_yyyyz,  \
                             tr_x_y_yyyzz,  \
                             tr_x_y_yyzzz,  \
                             tr_x_y_yzzzz,  \
                             ts_y_yyyyy,    \
                             ts_y_yyyyz,    \
                             ts_y_yyyzz,    \
                             ts_y_yyzzz,    \
                             ts_y_yzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_xxxxx[i] = tr_x_x_xxxxx[i] * pa_y[i];

        tr_x_xy_xxxxy[i] = tr_x_x_xxxx[i] * fe_0 + tr_x_x_xxxxy[i] * pa_y[i];

        tr_x_xy_xxxxz[i] = tr_x_x_xxxxz[i] * pa_y[i];

        tr_x_xy_xxxyy[i] = 2.0 * tr_x_x_xxxy[i] * fe_0 + tr_x_x_xxxyy[i] * pa_y[i];

        tr_x_xy_xxxyz[i] = tr_x_x_xxxz[i] * fe_0 + tr_x_x_xxxyz[i] * pa_y[i];

        tr_x_xy_xxxzz[i] = tr_x_x_xxxzz[i] * pa_y[i];

        tr_x_xy_xxyyy[i] = 3.0 * tr_x_x_xxyy[i] * fe_0 + tr_x_x_xxyyy[i] * pa_y[i];

        tr_x_xy_xxyyz[i] = 2.0 * tr_x_x_xxyz[i] * fe_0 + tr_x_x_xxyyz[i] * pa_y[i];

        tr_x_xy_xxyzz[i] = tr_x_x_xxzz[i] * fe_0 + tr_x_x_xxyzz[i] * pa_y[i];

        tr_x_xy_xxzzz[i] = tr_x_x_xxzzz[i] * pa_y[i];

        tr_x_xy_xyyyy[i] = 4.0 * tr_x_x_xyyy[i] * fe_0 + tr_x_x_xyyyy[i] * pa_y[i];

        tr_x_xy_xyyyz[i] = 3.0 * tr_x_x_xyyz[i] * fe_0 + tr_x_x_xyyyz[i] * pa_y[i];

        tr_x_xy_xyyzz[i] = 2.0 * tr_x_x_xyzz[i] * fe_0 + tr_x_x_xyyzz[i] * pa_y[i];

        tr_x_xy_xyzzz[i] = tr_x_x_xzzz[i] * fe_0 + tr_x_x_xyzzz[i] * pa_y[i];

        tr_x_xy_xzzzz[i] = tr_x_x_xzzzz[i] * pa_y[i];

        tr_x_xy_yyyyy[i] = ts_y_yyyyy[i] * fe_0 + tr_x_y_yyyyy[i] * pa_x[i];

        tr_x_xy_yyyyz[i] = ts_y_yyyyz[i] * fe_0 + tr_x_y_yyyyz[i] * pa_x[i];

        tr_x_xy_yyyzz[i] = ts_y_yyyzz[i] * fe_0 + tr_x_y_yyyzz[i] * pa_x[i];

        tr_x_xy_yyzzz[i] = ts_y_yyzzz[i] * fe_0 + tr_x_y_yyzzz[i] * pa_x[i];

        tr_x_xy_yzzzz[i] = ts_y_yzzzz[i] * fe_0 + tr_x_y_yzzzz[i] * pa_x[i];

        tr_x_xy_zzzzz[i] = tr_x_x_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : DH

    auto tr_x_xz_xxxxx = pbuffer.data(idx_dip_dh + 42);

    auto tr_x_xz_xxxxy = pbuffer.data(idx_dip_dh + 43);

    auto tr_x_xz_xxxxz = pbuffer.data(idx_dip_dh + 44);

    auto tr_x_xz_xxxyy = pbuffer.data(idx_dip_dh + 45);

    auto tr_x_xz_xxxyz = pbuffer.data(idx_dip_dh + 46);

    auto tr_x_xz_xxxzz = pbuffer.data(idx_dip_dh + 47);

    auto tr_x_xz_xxyyy = pbuffer.data(idx_dip_dh + 48);

    auto tr_x_xz_xxyyz = pbuffer.data(idx_dip_dh + 49);

    auto tr_x_xz_xxyzz = pbuffer.data(idx_dip_dh + 50);

    auto tr_x_xz_xxzzz = pbuffer.data(idx_dip_dh + 51);

    auto tr_x_xz_xyyyy = pbuffer.data(idx_dip_dh + 52);

    auto tr_x_xz_xyyyz = pbuffer.data(idx_dip_dh + 53);

    auto tr_x_xz_xyyzz = pbuffer.data(idx_dip_dh + 54);

    auto tr_x_xz_xyzzz = pbuffer.data(idx_dip_dh + 55);

    auto tr_x_xz_xzzzz = pbuffer.data(idx_dip_dh + 56);

    auto tr_x_xz_yyyyy = pbuffer.data(idx_dip_dh + 57);

    auto tr_x_xz_yyyyz = pbuffer.data(idx_dip_dh + 58);

    auto tr_x_xz_yyyzz = pbuffer.data(idx_dip_dh + 59);

    auto tr_x_xz_yyzzz = pbuffer.data(idx_dip_dh + 60);

    auto tr_x_xz_yzzzz = pbuffer.data(idx_dip_dh + 61);

    auto tr_x_xz_zzzzz = pbuffer.data(idx_dip_dh + 62);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_x_xxxx,   \
                             tr_x_x_xxxxx,  \
                             tr_x_x_xxxxy,  \
                             tr_x_x_xxxxz,  \
                             tr_x_x_xxxy,   \
                             tr_x_x_xxxyy,  \
                             tr_x_x_xxxyz,  \
                             tr_x_x_xxxz,   \
                             tr_x_x_xxxzz,  \
                             tr_x_x_xxyy,   \
                             tr_x_x_xxyyy,  \
                             tr_x_x_xxyyz,  \
                             tr_x_x_xxyz,   \
                             tr_x_x_xxyzz,  \
                             tr_x_x_xxzz,   \
                             tr_x_x_xxzzz,  \
                             tr_x_x_xyyy,   \
                             tr_x_x_xyyyy,  \
                             tr_x_x_xyyyz,  \
                             tr_x_x_xyyz,   \
                             tr_x_x_xyyzz,  \
                             tr_x_x_xyzz,   \
                             tr_x_x_xyzzz,  \
                             tr_x_x_xzzz,   \
                             tr_x_x_xzzzz,  \
                             tr_x_x_yyyyy,  \
                             tr_x_xz_xxxxx, \
                             tr_x_xz_xxxxy, \
                             tr_x_xz_xxxxz, \
                             tr_x_xz_xxxyy, \
                             tr_x_xz_xxxyz, \
                             tr_x_xz_xxxzz, \
                             tr_x_xz_xxyyy, \
                             tr_x_xz_xxyyz, \
                             tr_x_xz_xxyzz, \
                             tr_x_xz_xxzzz, \
                             tr_x_xz_xyyyy, \
                             tr_x_xz_xyyyz, \
                             tr_x_xz_xyyzz, \
                             tr_x_xz_xyzzz, \
                             tr_x_xz_xzzzz, \
                             tr_x_xz_yyyyy, \
                             tr_x_xz_yyyyz, \
                             tr_x_xz_yyyzz, \
                             tr_x_xz_yyzzz, \
                             tr_x_xz_yzzzz, \
                             tr_x_xz_zzzzz, \
                             tr_x_z_yyyyz,  \
                             tr_x_z_yyyzz,  \
                             tr_x_z_yyzzz,  \
                             tr_x_z_yzzzz,  \
                             tr_x_z_zzzzz,  \
                             ts_z_yyyyz,    \
                             ts_z_yyyzz,    \
                             ts_z_yyzzz,    \
                             ts_z_yzzzz,    \
                             ts_z_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_xxxxx[i] = tr_x_x_xxxxx[i] * pa_z[i];

        tr_x_xz_xxxxy[i] = tr_x_x_xxxxy[i] * pa_z[i];

        tr_x_xz_xxxxz[i] = tr_x_x_xxxx[i] * fe_0 + tr_x_x_xxxxz[i] * pa_z[i];

        tr_x_xz_xxxyy[i] = tr_x_x_xxxyy[i] * pa_z[i];

        tr_x_xz_xxxyz[i] = tr_x_x_xxxy[i] * fe_0 + tr_x_x_xxxyz[i] * pa_z[i];

        tr_x_xz_xxxzz[i] = 2.0 * tr_x_x_xxxz[i] * fe_0 + tr_x_x_xxxzz[i] * pa_z[i];

        tr_x_xz_xxyyy[i] = tr_x_x_xxyyy[i] * pa_z[i];

        tr_x_xz_xxyyz[i] = tr_x_x_xxyy[i] * fe_0 + tr_x_x_xxyyz[i] * pa_z[i];

        tr_x_xz_xxyzz[i] = 2.0 * tr_x_x_xxyz[i] * fe_0 + tr_x_x_xxyzz[i] * pa_z[i];

        tr_x_xz_xxzzz[i] = 3.0 * tr_x_x_xxzz[i] * fe_0 + tr_x_x_xxzzz[i] * pa_z[i];

        tr_x_xz_xyyyy[i] = tr_x_x_xyyyy[i] * pa_z[i];

        tr_x_xz_xyyyz[i] = tr_x_x_xyyy[i] * fe_0 + tr_x_x_xyyyz[i] * pa_z[i];

        tr_x_xz_xyyzz[i] = 2.0 * tr_x_x_xyyz[i] * fe_0 + tr_x_x_xyyzz[i] * pa_z[i];

        tr_x_xz_xyzzz[i] = 3.0 * tr_x_x_xyzz[i] * fe_0 + tr_x_x_xyzzz[i] * pa_z[i];

        tr_x_xz_xzzzz[i] = 4.0 * tr_x_x_xzzz[i] * fe_0 + tr_x_x_xzzzz[i] * pa_z[i];

        tr_x_xz_yyyyy[i] = tr_x_x_yyyyy[i] * pa_z[i];

        tr_x_xz_yyyyz[i] = ts_z_yyyyz[i] * fe_0 + tr_x_z_yyyyz[i] * pa_x[i];

        tr_x_xz_yyyzz[i] = ts_z_yyyzz[i] * fe_0 + tr_x_z_yyyzz[i] * pa_x[i];

        tr_x_xz_yyzzz[i] = ts_z_yyzzz[i] * fe_0 + tr_x_z_yyzzz[i] * pa_x[i];

        tr_x_xz_yzzzz[i] = ts_z_yzzzz[i] * fe_0 + tr_x_z_yzzzz[i] * pa_x[i];

        tr_x_xz_zzzzz[i] = ts_z_zzzzz[i] * fe_0 + tr_x_z_zzzzz[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : DH

    auto tr_x_yy_xxxxx = pbuffer.data(idx_dip_dh + 63);

    auto tr_x_yy_xxxxy = pbuffer.data(idx_dip_dh + 64);

    auto tr_x_yy_xxxxz = pbuffer.data(idx_dip_dh + 65);

    auto tr_x_yy_xxxyy = pbuffer.data(idx_dip_dh + 66);

    auto tr_x_yy_xxxyz = pbuffer.data(idx_dip_dh + 67);

    auto tr_x_yy_xxxzz = pbuffer.data(idx_dip_dh + 68);

    auto tr_x_yy_xxyyy = pbuffer.data(idx_dip_dh + 69);

    auto tr_x_yy_xxyyz = pbuffer.data(idx_dip_dh + 70);

    auto tr_x_yy_xxyzz = pbuffer.data(idx_dip_dh + 71);

    auto tr_x_yy_xxzzz = pbuffer.data(idx_dip_dh + 72);

    auto tr_x_yy_xyyyy = pbuffer.data(idx_dip_dh + 73);

    auto tr_x_yy_xyyyz = pbuffer.data(idx_dip_dh + 74);

    auto tr_x_yy_xyyzz = pbuffer.data(idx_dip_dh + 75);

    auto tr_x_yy_xyzzz = pbuffer.data(idx_dip_dh + 76);

    auto tr_x_yy_xzzzz = pbuffer.data(idx_dip_dh + 77);

    auto tr_x_yy_yyyyy = pbuffer.data(idx_dip_dh + 78);

    auto tr_x_yy_yyyyz = pbuffer.data(idx_dip_dh + 79);

    auto tr_x_yy_yyyzz = pbuffer.data(idx_dip_dh + 80);

    auto tr_x_yy_yyzzz = pbuffer.data(idx_dip_dh + 81);

    auto tr_x_yy_yzzzz = pbuffer.data(idx_dip_dh + 82);

    auto tr_x_yy_zzzzz = pbuffer.data(idx_dip_dh + 83);

#pragma omp simd aligned(pa_y,              \
                             tr_x_0_xxxxx,  \
                             tr_x_0_xxxxy,  \
                             tr_x_0_xxxxz,  \
                             tr_x_0_xxxyy,  \
                             tr_x_0_xxxyz,  \
                             tr_x_0_xxxzz,  \
                             tr_x_0_xxyyy,  \
                             tr_x_0_xxyyz,  \
                             tr_x_0_xxyzz,  \
                             tr_x_0_xxzzz,  \
                             tr_x_0_xyyyy,  \
                             tr_x_0_xyyyz,  \
                             tr_x_0_xyyzz,  \
                             tr_x_0_xyzzz,  \
                             tr_x_0_xzzzz,  \
                             tr_x_0_yyyyy,  \
                             tr_x_0_yyyyz,  \
                             tr_x_0_yyyzz,  \
                             tr_x_0_yyzzz,  \
                             tr_x_0_yzzzz,  \
                             tr_x_0_zzzzz,  \
                             tr_x_y_xxxx,   \
                             tr_x_y_xxxxx,  \
                             tr_x_y_xxxxy,  \
                             tr_x_y_xxxxz,  \
                             tr_x_y_xxxy,   \
                             tr_x_y_xxxyy,  \
                             tr_x_y_xxxyz,  \
                             tr_x_y_xxxz,   \
                             tr_x_y_xxxzz,  \
                             tr_x_y_xxyy,   \
                             tr_x_y_xxyyy,  \
                             tr_x_y_xxyyz,  \
                             tr_x_y_xxyz,   \
                             tr_x_y_xxyzz,  \
                             tr_x_y_xxzz,   \
                             tr_x_y_xxzzz,  \
                             tr_x_y_xyyy,   \
                             tr_x_y_xyyyy,  \
                             tr_x_y_xyyyz,  \
                             tr_x_y_xyyz,   \
                             tr_x_y_xyyzz,  \
                             tr_x_y_xyzz,   \
                             tr_x_y_xyzzz,  \
                             tr_x_y_xzzz,   \
                             tr_x_y_xzzzz,  \
                             tr_x_y_yyyy,   \
                             tr_x_y_yyyyy,  \
                             tr_x_y_yyyyz,  \
                             tr_x_y_yyyz,   \
                             tr_x_y_yyyzz,  \
                             tr_x_y_yyzz,   \
                             tr_x_y_yyzzz,  \
                             tr_x_y_yzzz,   \
                             tr_x_y_yzzzz,  \
                             tr_x_y_zzzz,   \
                             tr_x_y_zzzzz,  \
                             tr_x_yy_xxxxx, \
                             tr_x_yy_xxxxy, \
                             tr_x_yy_xxxxz, \
                             tr_x_yy_xxxyy, \
                             tr_x_yy_xxxyz, \
                             tr_x_yy_xxxzz, \
                             tr_x_yy_xxyyy, \
                             tr_x_yy_xxyyz, \
                             tr_x_yy_xxyzz, \
                             tr_x_yy_xxzzz, \
                             tr_x_yy_xyyyy, \
                             tr_x_yy_xyyyz, \
                             tr_x_yy_xyyzz, \
                             tr_x_yy_xyzzz, \
                             tr_x_yy_xzzzz, \
                             tr_x_yy_yyyyy, \
                             tr_x_yy_yyyyz, \
                             tr_x_yy_yyyzz, \
                             tr_x_yy_yyzzz, \
                             tr_x_yy_yzzzz, \
                             tr_x_yy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_xxxxx[i] = tr_x_0_xxxxx[i] * fe_0 + tr_x_y_xxxxx[i] * pa_y[i];

        tr_x_yy_xxxxy[i] = tr_x_0_xxxxy[i] * fe_0 + tr_x_y_xxxx[i] * fe_0 + tr_x_y_xxxxy[i] * pa_y[i];

        tr_x_yy_xxxxz[i] = tr_x_0_xxxxz[i] * fe_0 + tr_x_y_xxxxz[i] * pa_y[i];

        tr_x_yy_xxxyy[i] = tr_x_0_xxxyy[i] * fe_0 + 2.0 * tr_x_y_xxxy[i] * fe_0 + tr_x_y_xxxyy[i] * pa_y[i];

        tr_x_yy_xxxyz[i] = tr_x_0_xxxyz[i] * fe_0 + tr_x_y_xxxz[i] * fe_0 + tr_x_y_xxxyz[i] * pa_y[i];

        tr_x_yy_xxxzz[i] = tr_x_0_xxxzz[i] * fe_0 + tr_x_y_xxxzz[i] * pa_y[i];

        tr_x_yy_xxyyy[i] = tr_x_0_xxyyy[i] * fe_0 + 3.0 * tr_x_y_xxyy[i] * fe_0 + tr_x_y_xxyyy[i] * pa_y[i];

        tr_x_yy_xxyyz[i] = tr_x_0_xxyyz[i] * fe_0 + 2.0 * tr_x_y_xxyz[i] * fe_0 + tr_x_y_xxyyz[i] * pa_y[i];

        tr_x_yy_xxyzz[i] = tr_x_0_xxyzz[i] * fe_0 + tr_x_y_xxzz[i] * fe_0 + tr_x_y_xxyzz[i] * pa_y[i];

        tr_x_yy_xxzzz[i] = tr_x_0_xxzzz[i] * fe_0 + tr_x_y_xxzzz[i] * pa_y[i];

        tr_x_yy_xyyyy[i] = tr_x_0_xyyyy[i] * fe_0 + 4.0 * tr_x_y_xyyy[i] * fe_0 + tr_x_y_xyyyy[i] * pa_y[i];

        tr_x_yy_xyyyz[i] = tr_x_0_xyyyz[i] * fe_0 + 3.0 * tr_x_y_xyyz[i] * fe_0 + tr_x_y_xyyyz[i] * pa_y[i];

        tr_x_yy_xyyzz[i] = tr_x_0_xyyzz[i] * fe_0 + 2.0 * tr_x_y_xyzz[i] * fe_0 + tr_x_y_xyyzz[i] * pa_y[i];

        tr_x_yy_xyzzz[i] = tr_x_0_xyzzz[i] * fe_0 + tr_x_y_xzzz[i] * fe_0 + tr_x_y_xyzzz[i] * pa_y[i];

        tr_x_yy_xzzzz[i] = tr_x_0_xzzzz[i] * fe_0 + tr_x_y_xzzzz[i] * pa_y[i];

        tr_x_yy_yyyyy[i] = tr_x_0_yyyyy[i] * fe_0 + 5.0 * tr_x_y_yyyy[i] * fe_0 + tr_x_y_yyyyy[i] * pa_y[i];

        tr_x_yy_yyyyz[i] = tr_x_0_yyyyz[i] * fe_0 + 4.0 * tr_x_y_yyyz[i] * fe_0 + tr_x_y_yyyyz[i] * pa_y[i];

        tr_x_yy_yyyzz[i] = tr_x_0_yyyzz[i] * fe_0 + 3.0 * tr_x_y_yyzz[i] * fe_0 + tr_x_y_yyyzz[i] * pa_y[i];

        tr_x_yy_yyzzz[i] = tr_x_0_yyzzz[i] * fe_0 + 2.0 * tr_x_y_yzzz[i] * fe_0 + tr_x_y_yyzzz[i] * pa_y[i];

        tr_x_yy_yzzzz[i] = tr_x_0_yzzzz[i] * fe_0 + tr_x_y_zzzz[i] * fe_0 + tr_x_y_yzzzz[i] * pa_y[i];

        tr_x_yy_zzzzz[i] = tr_x_0_zzzzz[i] * fe_0 + tr_x_y_zzzzz[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : DH

    auto tr_x_yz_xxxxx = pbuffer.data(idx_dip_dh + 84);

    auto tr_x_yz_xxxxy = pbuffer.data(idx_dip_dh + 85);

    auto tr_x_yz_xxxxz = pbuffer.data(idx_dip_dh + 86);

    auto tr_x_yz_xxxyy = pbuffer.data(idx_dip_dh + 87);

    auto tr_x_yz_xxxyz = pbuffer.data(idx_dip_dh + 88);

    auto tr_x_yz_xxxzz = pbuffer.data(idx_dip_dh + 89);

    auto tr_x_yz_xxyyy = pbuffer.data(idx_dip_dh + 90);

    auto tr_x_yz_xxyyz = pbuffer.data(idx_dip_dh + 91);

    auto tr_x_yz_xxyzz = pbuffer.data(idx_dip_dh + 92);

    auto tr_x_yz_xxzzz = pbuffer.data(idx_dip_dh + 93);

    auto tr_x_yz_xyyyy = pbuffer.data(idx_dip_dh + 94);

    auto tr_x_yz_xyyyz = pbuffer.data(idx_dip_dh + 95);

    auto tr_x_yz_xyyzz = pbuffer.data(idx_dip_dh + 96);

    auto tr_x_yz_xyzzz = pbuffer.data(idx_dip_dh + 97);

    auto tr_x_yz_xzzzz = pbuffer.data(idx_dip_dh + 98);

    auto tr_x_yz_yyyyy = pbuffer.data(idx_dip_dh + 99);

    auto tr_x_yz_yyyyz = pbuffer.data(idx_dip_dh + 100);

    auto tr_x_yz_yyyzz = pbuffer.data(idx_dip_dh + 101);

    auto tr_x_yz_yyzzz = pbuffer.data(idx_dip_dh + 102);

    auto tr_x_yz_yzzzz = pbuffer.data(idx_dip_dh + 103);

    auto tr_x_yz_zzzzz = pbuffer.data(idx_dip_dh + 104);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_y_xxxxy,  \
                             tr_x_y_xxxyy,  \
                             tr_x_y_xxyyy,  \
                             tr_x_y_xyyyy,  \
                             tr_x_y_yyyyy,  \
                             tr_x_yz_xxxxx, \
                             tr_x_yz_xxxxy, \
                             tr_x_yz_xxxxz, \
                             tr_x_yz_xxxyy, \
                             tr_x_yz_xxxyz, \
                             tr_x_yz_xxxzz, \
                             tr_x_yz_xxyyy, \
                             tr_x_yz_xxyyz, \
                             tr_x_yz_xxyzz, \
                             tr_x_yz_xxzzz, \
                             tr_x_yz_xyyyy, \
                             tr_x_yz_xyyyz, \
                             tr_x_yz_xyyzz, \
                             tr_x_yz_xyzzz, \
                             tr_x_yz_xzzzz, \
                             tr_x_yz_yyyyy, \
                             tr_x_yz_yyyyz, \
                             tr_x_yz_yyyzz, \
                             tr_x_yz_yyzzz, \
                             tr_x_yz_yzzzz, \
                             tr_x_yz_zzzzz, \
                             tr_x_z_xxxxx,  \
                             tr_x_z_xxxxz,  \
                             tr_x_z_xxxyz,  \
                             tr_x_z_xxxz,   \
                             tr_x_z_xxxzz,  \
                             tr_x_z_xxyyz,  \
                             tr_x_z_xxyz,   \
                             tr_x_z_xxyzz,  \
                             tr_x_z_xxzz,   \
                             tr_x_z_xxzzz,  \
                             tr_x_z_xyyyz,  \
                             tr_x_z_xyyz,   \
                             tr_x_z_xyyzz,  \
                             tr_x_z_xyzz,   \
                             tr_x_z_xyzzz,  \
                             tr_x_z_xzzz,   \
                             tr_x_z_xzzzz,  \
                             tr_x_z_yyyyz,  \
                             tr_x_z_yyyz,   \
                             tr_x_z_yyyzz,  \
                             tr_x_z_yyzz,   \
                             tr_x_z_yyzzz,  \
                             tr_x_z_yzzz,   \
                             tr_x_z_yzzzz,  \
                             tr_x_z_zzzz,   \
                             tr_x_z_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yz_xxxxx[i] = tr_x_z_xxxxx[i] * pa_y[i];

        tr_x_yz_xxxxy[i] = tr_x_y_xxxxy[i] * pa_z[i];

        tr_x_yz_xxxxz[i] = tr_x_z_xxxxz[i] * pa_y[i];

        tr_x_yz_xxxyy[i] = tr_x_y_xxxyy[i] * pa_z[i];

        tr_x_yz_xxxyz[i] = tr_x_z_xxxz[i] * fe_0 + tr_x_z_xxxyz[i] * pa_y[i];

        tr_x_yz_xxxzz[i] = tr_x_z_xxxzz[i] * pa_y[i];

        tr_x_yz_xxyyy[i] = tr_x_y_xxyyy[i] * pa_z[i];

        tr_x_yz_xxyyz[i] = 2.0 * tr_x_z_xxyz[i] * fe_0 + tr_x_z_xxyyz[i] * pa_y[i];

        tr_x_yz_xxyzz[i] = tr_x_z_xxzz[i] * fe_0 + tr_x_z_xxyzz[i] * pa_y[i];

        tr_x_yz_xxzzz[i] = tr_x_z_xxzzz[i] * pa_y[i];

        tr_x_yz_xyyyy[i] = tr_x_y_xyyyy[i] * pa_z[i];

        tr_x_yz_xyyyz[i] = 3.0 * tr_x_z_xyyz[i] * fe_0 + tr_x_z_xyyyz[i] * pa_y[i];

        tr_x_yz_xyyzz[i] = 2.0 * tr_x_z_xyzz[i] * fe_0 + tr_x_z_xyyzz[i] * pa_y[i];

        tr_x_yz_xyzzz[i] = tr_x_z_xzzz[i] * fe_0 + tr_x_z_xyzzz[i] * pa_y[i];

        tr_x_yz_xzzzz[i] = tr_x_z_xzzzz[i] * pa_y[i];

        tr_x_yz_yyyyy[i] = tr_x_y_yyyyy[i] * pa_z[i];

        tr_x_yz_yyyyz[i] = 4.0 * tr_x_z_yyyz[i] * fe_0 + tr_x_z_yyyyz[i] * pa_y[i];

        tr_x_yz_yyyzz[i] = 3.0 * tr_x_z_yyzz[i] * fe_0 + tr_x_z_yyyzz[i] * pa_y[i];

        tr_x_yz_yyzzz[i] = 2.0 * tr_x_z_yzzz[i] * fe_0 + tr_x_z_yyzzz[i] * pa_y[i];

        tr_x_yz_yzzzz[i] = tr_x_z_zzzz[i] * fe_0 + tr_x_z_yzzzz[i] * pa_y[i];

        tr_x_yz_zzzzz[i] = tr_x_z_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : DH

    auto tr_x_zz_xxxxx = pbuffer.data(idx_dip_dh + 105);

    auto tr_x_zz_xxxxy = pbuffer.data(idx_dip_dh + 106);

    auto tr_x_zz_xxxxz = pbuffer.data(idx_dip_dh + 107);

    auto tr_x_zz_xxxyy = pbuffer.data(idx_dip_dh + 108);

    auto tr_x_zz_xxxyz = pbuffer.data(idx_dip_dh + 109);

    auto tr_x_zz_xxxzz = pbuffer.data(idx_dip_dh + 110);

    auto tr_x_zz_xxyyy = pbuffer.data(idx_dip_dh + 111);

    auto tr_x_zz_xxyyz = pbuffer.data(idx_dip_dh + 112);

    auto tr_x_zz_xxyzz = pbuffer.data(idx_dip_dh + 113);

    auto tr_x_zz_xxzzz = pbuffer.data(idx_dip_dh + 114);

    auto tr_x_zz_xyyyy = pbuffer.data(idx_dip_dh + 115);

    auto tr_x_zz_xyyyz = pbuffer.data(idx_dip_dh + 116);

    auto tr_x_zz_xyyzz = pbuffer.data(idx_dip_dh + 117);

    auto tr_x_zz_xyzzz = pbuffer.data(idx_dip_dh + 118);

    auto tr_x_zz_xzzzz = pbuffer.data(idx_dip_dh + 119);

    auto tr_x_zz_yyyyy = pbuffer.data(idx_dip_dh + 120);

    auto tr_x_zz_yyyyz = pbuffer.data(idx_dip_dh + 121);

    auto tr_x_zz_yyyzz = pbuffer.data(idx_dip_dh + 122);

    auto tr_x_zz_yyzzz = pbuffer.data(idx_dip_dh + 123);

    auto tr_x_zz_yzzzz = pbuffer.data(idx_dip_dh + 124);

    auto tr_x_zz_zzzzz = pbuffer.data(idx_dip_dh + 125);

#pragma omp simd aligned(pa_z,              \
                             tr_x_0_xxxxx,  \
                             tr_x_0_xxxxy,  \
                             tr_x_0_xxxxz,  \
                             tr_x_0_xxxyy,  \
                             tr_x_0_xxxyz,  \
                             tr_x_0_xxxzz,  \
                             tr_x_0_xxyyy,  \
                             tr_x_0_xxyyz,  \
                             tr_x_0_xxyzz,  \
                             tr_x_0_xxzzz,  \
                             tr_x_0_xyyyy,  \
                             tr_x_0_xyyyz,  \
                             tr_x_0_xyyzz,  \
                             tr_x_0_xyzzz,  \
                             tr_x_0_xzzzz,  \
                             tr_x_0_yyyyy,  \
                             tr_x_0_yyyyz,  \
                             tr_x_0_yyyzz,  \
                             tr_x_0_yyzzz,  \
                             tr_x_0_yzzzz,  \
                             tr_x_0_zzzzz,  \
                             tr_x_z_xxxx,   \
                             tr_x_z_xxxxx,  \
                             tr_x_z_xxxxy,  \
                             tr_x_z_xxxxz,  \
                             tr_x_z_xxxy,   \
                             tr_x_z_xxxyy,  \
                             tr_x_z_xxxyz,  \
                             tr_x_z_xxxz,   \
                             tr_x_z_xxxzz,  \
                             tr_x_z_xxyy,   \
                             tr_x_z_xxyyy,  \
                             tr_x_z_xxyyz,  \
                             tr_x_z_xxyz,   \
                             tr_x_z_xxyzz,  \
                             tr_x_z_xxzz,   \
                             tr_x_z_xxzzz,  \
                             tr_x_z_xyyy,   \
                             tr_x_z_xyyyy,  \
                             tr_x_z_xyyyz,  \
                             tr_x_z_xyyz,   \
                             tr_x_z_xyyzz,  \
                             tr_x_z_xyzz,   \
                             tr_x_z_xyzzz,  \
                             tr_x_z_xzzz,   \
                             tr_x_z_xzzzz,  \
                             tr_x_z_yyyy,   \
                             tr_x_z_yyyyy,  \
                             tr_x_z_yyyyz,  \
                             tr_x_z_yyyz,   \
                             tr_x_z_yyyzz,  \
                             tr_x_z_yyzz,   \
                             tr_x_z_yyzzz,  \
                             tr_x_z_yzzz,   \
                             tr_x_z_yzzzz,  \
                             tr_x_z_zzzz,   \
                             tr_x_z_zzzzz,  \
                             tr_x_zz_xxxxx, \
                             tr_x_zz_xxxxy, \
                             tr_x_zz_xxxxz, \
                             tr_x_zz_xxxyy, \
                             tr_x_zz_xxxyz, \
                             tr_x_zz_xxxzz, \
                             tr_x_zz_xxyyy, \
                             tr_x_zz_xxyyz, \
                             tr_x_zz_xxyzz, \
                             tr_x_zz_xxzzz, \
                             tr_x_zz_xyyyy, \
                             tr_x_zz_xyyyz, \
                             tr_x_zz_xyyzz, \
                             tr_x_zz_xyzzz, \
                             tr_x_zz_xzzzz, \
                             tr_x_zz_yyyyy, \
                             tr_x_zz_yyyyz, \
                             tr_x_zz_yyyzz, \
                             tr_x_zz_yyzzz, \
                             tr_x_zz_yzzzz, \
                             tr_x_zz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_xxxxx[i] = tr_x_0_xxxxx[i] * fe_0 + tr_x_z_xxxxx[i] * pa_z[i];

        tr_x_zz_xxxxy[i] = tr_x_0_xxxxy[i] * fe_0 + tr_x_z_xxxxy[i] * pa_z[i];

        tr_x_zz_xxxxz[i] = tr_x_0_xxxxz[i] * fe_0 + tr_x_z_xxxx[i] * fe_0 + tr_x_z_xxxxz[i] * pa_z[i];

        tr_x_zz_xxxyy[i] = tr_x_0_xxxyy[i] * fe_0 + tr_x_z_xxxyy[i] * pa_z[i];

        tr_x_zz_xxxyz[i] = tr_x_0_xxxyz[i] * fe_0 + tr_x_z_xxxy[i] * fe_0 + tr_x_z_xxxyz[i] * pa_z[i];

        tr_x_zz_xxxzz[i] = tr_x_0_xxxzz[i] * fe_0 + 2.0 * tr_x_z_xxxz[i] * fe_0 + tr_x_z_xxxzz[i] * pa_z[i];

        tr_x_zz_xxyyy[i] = tr_x_0_xxyyy[i] * fe_0 + tr_x_z_xxyyy[i] * pa_z[i];

        tr_x_zz_xxyyz[i] = tr_x_0_xxyyz[i] * fe_0 + tr_x_z_xxyy[i] * fe_0 + tr_x_z_xxyyz[i] * pa_z[i];

        tr_x_zz_xxyzz[i] = tr_x_0_xxyzz[i] * fe_0 + 2.0 * tr_x_z_xxyz[i] * fe_0 + tr_x_z_xxyzz[i] * pa_z[i];

        tr_x_zz_xxzzz[i] = tr_x_0_xxzzz[i] * fe_0 + 3.0 * tr_x_z_xxzz[i] * fe_0 + tr_x_z_xxzzz[i] * pa_z[i];

        tr_x_zz_xyyyy[i] = tr_x_0_xyyyy[i] * fe_0 + tr_x_z_xyyyy[i] * pa_z[i];

        tr_x_zz_xyyyz[i] = tr_x_0_xyyyz[i] * fe_0 + tr_x_z_xyyy[i] * fe_0 + tr_x_z_xyyyz[i] * pa_z[i];

        tr_x_zz_xyyzz[i] = tr_x_0_xyyzz[i] * fe_0 + 2.0 * tr_x_z_xyyz[i] * fe_0 + tr_x_z_xyyzz[i] * pa_z[i];

        tr_x_zz_xyzzz[i] = tr_x_0_xyzzz[i] * fe_0 + 3.0 * tr_x_z_xyzz[i] * fe_0 + tr_x_z_xyzzz[i] * pa_z[i];

        tr_x_zz_xzzzz[i] = tr_x_0_xzzzz[i] * fe_0 + 4.0 * tr_x_z_xzzz[i] * fe_0 + tr_x_z_xzzzz[i] * pa_z[i];

        tr_x_zz_yyyyy[i] = tr_x_0_yyyyy[i] * fe_0 + tr_x_z_yyyyy[i] * pa_z[i];

        tr_x_zz_yyyyz[i] = tr_x_0_yyyyz[i] * fe_0 + tr_x_z_yyyy[i] * fe_0 + tr_x_z_yyyyz[i] * pa_z[i];

        tr_x_zz_yyyzz[i] = tr_x_0_yyyzz[i] * fe_0 + 2.0 * tr_x_z_yyyz[i] * fe_0 + tr_x_z_yyyzz[i] * pa_z[i];

        tr_x_zz_yyzzz[i] = tr_x_0_yyzzz[i] * fe_0 + 3.0 * tr_x_z_yyzz[i] * fe_0 + tr_x_z_yyzzz[i] * pa_z[i];

        tr_x_zz_yzzzz[i] = tr_x_0_yzzzz[i] * fe_0 + 4.0 * tr_x_z_yzzz[i] * fe_0 + tr_x_z_yzzzz[i] * pa_z[i];

        tr_x_zz_zzzzz[i] = tr_x_0_zzzzz[i] * fe_0 + 5.0 * tr_x_z_zzzz[i] * fe_0 + tr_x_z_zzzzz[i] * pa_z[i];
    }

    // Set up 126-147 components of targeted buffer : DH

    auto tr_y_xx_xxxxx = pbuffer.data(idx_dip_dh + 126);

    auto tr_y_xx_xxxxy = pbuffer.data(idx_dip_dh + 127);

    auto tr_y_xx_xxxxz = pbuffer.data(idx_dip_dh + 128);

    auto tr_y_xx_xxxyy = pbuffer.data(idx_dip_dh + 129);

    auto tr_y_xx_xxxyz = pbuffer.data(idx_dip_dh + 130);

    auto tr_y_xx_xxxzz = pbuffer.data(idx_dip_dh + 131);

    auto tr_y_xx_xxyyy = pbuffer.data(idx_dip_dh + 132);

    auto tr_y_xx_xxyyz = pbuffer.data(idx_dip_dh + 133);

    auto tr_y_xx_xxyzz = pbuffer.data(idx_dip_dh + 134);

    auto tr_y_xx_xxzzz = pbuffer.data(idx_dip_dh + 135);

    auto tr_y_xx_xyyyy = pbuffer.data(idx_dip_dh + 136);

    auto tr_y_xx_xyyyz = pbuffer.data(idx_dip_dh + 137);

    auto tr_y_xx_xyyzz = pbuffer.data(idx_dip_dh + 138);

    auto tr_y_xx_xyzzz = pbuffer.data(idx_dip_dh + 139);

    auto tr_y_xx_xzzzz = pbuffer.data(idx_dip_dh + 140);

    auto tr_y_xx_yyyyy = pbuffer.data(idx_dip_dh + 141);

    auto tr_y_xx_yyyyz = pbuffer.data(idx_dip_dh + 142);

    auto tr_y_xx_yyyzz = pbuffer.data(idx_dip_dh + 143);

    auto tr_y_xx_yyzzz = pbuffer.data(idx_dip_dh + 144);

    auto tr_y_xx_yzzzz = pbuffer.data(idx_dip_dh + 145);

    auto tr_y_xx_zzzzz = pbuffer.data(idx_dip_dh + 146);

#pragma omp simd aligned(pa_x,              \
                             tr_y_0_xxxxx,  \
                             tr_y_0_xxxxy,  \
                             tr_y_0_xxxxz,  \
                             tr_y_0_xxxyy,  \
                             tr_y_0_xxxyz,  \
                             tr_y_0_xxxzz,  \
                             tr_y_0_xxyyy,  \
                             tr_y_0_xxyyz,  \
                             tr_y_0_xxyzz,  \
                             tr_y_0_xxzzz,  \
                             tr_y_0_xyyyy,  \
                             tr_y_0_xyyyz,  \
                             tr_y_0_xyyzz,  \
                             tr_y_0_xyzzz,  \
                             tr_y_0_xzzzz,  \
                             tr_y_0_yyyyy,  \
                             tr_y_0_yyyyz,  \
                             tr_y_0_yyyzz,  \
                             tr_y_0_yyzzz,  \
                             tr_y_0_yzzzz,  \
                             tr_y_0_zzzzz,  \
                             tr_y_x_xxxx,   \
                             tr_y_x_xxxxx,  \
                             tr_y_x_xxxxy,  \
                             tr_y_x_xxxxz,  \
                             tr_y_x_xxxy,   \
                             tr_y_x_xxxyy,  \
                             tr_y_x_xxxyz,  \
                             tr_y_x_xxxz,   \
                             tr_y_x_xxxzz,  \
                             tr_y_x_xxyy,   \
                             tr_y_x_xxyyy,  \
                             tr_y_x_xxyyz,  \
                             tr_y_x_xxyz,   \
                             tr_y_x_xxyzz,  \
                             tr_y_x_xxzz,   \
                             tr_y_x_xxzzz,  \
                             tr_y_x_xyyy,   \
                             tr_y_x_xyyyy,  \
                             tr_y_x_xyyyz,  \
                             tr_y_x_xyyz,   \
                             tr_y_x_xyyzz,  \
                             tr_y_x_xyzz,   \
                             tr_y_x_xyzzz,  \
                             tr_y_x_xzzz,   \
                             tr_y_x_xzzzz,  \
                             tr_y_x_yyyy,   \
                             tr_y_x_yyyyy,  \
                             tr_y_x_yyyyz,  \
                             tr_y_x_yyyz,   \
                             tr_y_x_yyyzz,  \
                             tr_y_x_yyzz,   \
                             tr_y_x_yyzzz,  \
                             tr_y_x_yzzz,   \
                             tr_y_x_yzzzz,  \
                             tr_y_x_zzzz,   \
                             tr_y_x_zzzzz,  \
                             tr_y_xx_xxxxx, \
                             tr_y_xx_xxxxy, \
                             tr_y_xx_xxxxz, \
                             tr_y_xx_xxxyy, \
                             tr_y_xx_xxxyz, \
                             tr_y_xx_xxxzz, \
                             tr_y_xx_xxyyy, \
                             tr_y_xx_xxyyz, \
                             tr_y_xx_xxyzz, \
                             tr_y_xx_xxzzz, \
                             tr_y_xx_xyyyy, \
                             tr_y_xx_xyyyz, \
                             tr_y_xx_xyyzz, \
                             tr_y_xx_xyzzz, \
                             tr_y_xx_xzzzz, \
                             tr_y_xx_yyyyy, \
                             tr_y_xx_yyyyz, \
                             tr_y_xx_yyyzz, \
                             tr_y_xx_yyzzz, \
                             tr_y_xx_yzzzz, \
                             tr_y_xx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_xxxxx[i] = tr_y_0_xxxxx[i] * fe_0 + 5.0 * tr_y_x_xxxx[i] * fe_0 + tr_y_x_xxxxx[i] * pa_x[i];

        tr_y_xx_xxxxy[i] = tr_y_0_xxxxy[i] * fe_0 + 4.0 * tr_y_x_xxxy[i] * fe_0 + tr_y_x_xxxxy[i] * pa_x[i];

        tr_y_xx_xxxxz[i] = tr_y_0_xxxxz[i] * fe_0 + 4.0 * tr_y_x_xxxz[i] * fe_0 + tr_y_x_xxxxz[i] * pa_x[i];

        tr_y_xx_xxxyy[i] = tr_y_0_xxxyy[i] * fe_0 + 3.0 * tr_y_x_xxyy[i] * fe_0 + tr_y_x_xxxyy[i] * pa_x[i];

        tr_y_xx_xxxyz[i] = tr_y_0_xxxyz[i] * fe_0 + 3.0 * tr_y_x_xxyz[i] * fe_0 + tr_y_x_xxxyz[i] * pa_x[i];

        tr_y_xx_xxxzz[i] = tr_y_0_xxxzz[i] * fe_0 + 3.0 * tr_y_x_xxzz[i] * fe_0 + tr_y_x_xxxzz[i] * pa_x[i];

        tr_y_xx_xxyyy[i] = tr_y_0_xxyyy[i] * fe_0 + 2.0 * tr_y_x_xyyy[i] * fe_0 + tr_y_x_xxyyy[i] * pa_x[i];

        tr_y_xx_xxyyz[i] = tr_y_0_xxyyz[i] * fe_0 + 2.0 * tr_y_x_xyyz[i] * fe_0 + tr_y_x_xxyyz[i] * pa_x[i];

        tr_y_xx_xxyzz[i] = tr_y_0_xxyzz[i] * fe_0 + 2.0 * tr_y_x_xyzz[i] * fe_0 + tr_y_x_xxyzz[i] * pa_x[i];

        tr_y_xx_xxzzz[i] = tr_y_0_xxzzz[i] * fe_0 + 2.0 * tr_y_x_xzzz[i] * fe_0 + tr_y_x_xxzzz[i] * pa_x[i];

        tr_y_xx_xyyyy[i] = tr_y_0_xyyyy[i] * fe_0 + tr_y_x_yyyy[i] * fe_0 + tr_y_x_xyyyy[i] * pa_x[i];

        tr_y_xx_xyyyz[i] = tr_y_0_xyyyz[i] * fe_0 + tr_y_x_yyyz[i] * fe_0 + tr_y_x_xyyyz[i] * pa_x[i];

        tr_y_xx_xyyzz[i] = tr_y_0_xyyzz[i] * fe_0 + tr_y_x_yyzz[i] * fe_0 + tr_y_x_xyyzz[i] * pa_x[i];

        tr_y_xx_xyzzz[i] = tr_y_0_xyzzz[i] * fe_0 + tr_y_x_yzzz[i] * fe_0 + tr_y_x_xyzzz[i] * pa_x[i];

        tr_y_xx_xzzzz[i] = tr_y_0_xzzzz[i] * fe_0 + tr_y_x_zzzz[i] * fe_0 + tr_y_x_xzzzz[i] * pa_x[i];

        tr_y_xx_yyyyy[i] = tr_y_0_yyyyy[i] * fe_0 + tr_y_x_yyyyy[i] * pa_x[i];

        tr_y_xx_yyyyz[i] = tr_y_0_yyyyz[i] * fe_0 + tr_y_x_yyyyz[i] * pa_x[i];

        tr_y_xx_yyyzz[i] = tr_y_0_yyyzz[i] * fe_0 + tr_y_x_yyyzz[i] * pa_x[i];

        tr_y_xx_yyzzz[i] = tr_y_0_yyzzz[i] * fe_0 + tr_y_x_yyzzz[i] * pa_x[i];

        tr_y_xx_yzzzz[i] = tr_y_0_yzzzz[i] * fe_0 + tr_y_x_yzzzz[i] * pa_x[i];

        tr_y_xx_zzzzz[i] = tr_y_0_zzzzz[i] * fe_0 + tr_y_x_zzzzz[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : DH

    auto tr_y_xy_xxxxx = pbuffer.data(idx_dip_dh + 147);

    auto tr_y_xy_xxxxy = pbuffer.data(idx_dip_dh + 148);

    auto tr_y_xy_xxxxz = pbuffer.data(idx_dip_dh + 149);

    auto tr_y_xy_xxxyy = pbuffer.data(idx_dip_dh + 150);

    auto tr_y_xy_xxxyz = pbuffer.data(idx_dip_dh + 151);

    auto tr_y_xy_xxxzz = pbuffer.data(idx_dip_dh + 152);

    auto tr_y_xy_xxyyy = pbuffer.data(idx_dip_dh + 153);

    auto tr_y_xy_xxyyz = pbuffer.data(idx_dip_dh + 154);

    auto tr_y_xy_xxyzz = pbuffer.data(idx_dip_dh + 155);

    auto tr_y_xy_xxzzz = pbuffer.data(idx_dip_dh + 156);

    auto tr_y_xy_xyyyy = pbuffer.data(idx_dip_dh + 157);

    auto tr_y_xy_xyyyz = pbuffer.data(idx_dip_dh + 158);

    auto tr_y_xy_xyyzz = pbuffer.data(idx_dip_dh + 159);

    auto tr_y_xy_xyzzz = pbuffer.data(idx_dip_dh + 160);

    auto tr_y_xy_xzzzz = pbuffer.data(idx_dip_dh + 161);

    auto tr_y_xy_yyyyy = pbuffer.data(idx_dip_dh + 162);

    auto tr_y_xy_yyyyz = pbuffer.data(idx_dip_dh + 163);

    auto tr_y_xy_yyyzz = pbuffer.data(idx_dip_dh + 164);

    auto tr_y_xy_yyzzz = pbuffer.data(idx_dip_dh + 165);

    auto tr_y_xy_yzzzz = pbuffer.data(idx_dip_dh + 166);

    auto tr_y_xy_zzzzz = pbuffer.data(idx_dip_dh + 167);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xy_xxxxx, \
                             tr_y_xy_xxxxy, \
                             tr_y_xy_xxxxz, \
                             tr_y_xy_xxxyy, \
                             tr_y_xy_xxxyz, \
                             tr_y_xy_xxxzz, \
                             tr_y_xy_xxyyy, \
                             tr_y_xy_xxyyz, \
                             tr_y_xy_xxyzz, \
                             tr_y_xy_xxzzz, \
                             tr_y_xy_xyyyy, \
                             tr_y_xy_xyyyz, \
                             tr_y_xy_xyyzz, \
                             tr_y_xy_xyzzz, \
                             tr_y_xy_xzzzz, \
                             tr_y_xy_yyyyy, \
                             tr_y_xy_yyyyz, \
                             tr_y_xy_yyyzz, \
                             tr_y_xy_yyzzz, \
                             tr_y_xy_yzzzz, \
                             tr_y_xy_zzzzz, \
                             tr_y_y_xxxx,   \
                             tr_y_y_xxxxx,  \
                             tr_y_y_xxxxy,  \
                             tr_y_y_xxxxz,  \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxxyy,  \
                             tr_y_y_xxxyz,  \
                             tr_y_y_xxxz,   \
                             tr_y_y_xxxzz,  \
                             tr_y_y_xxyy,   \
                             tr_y_y_xxyyy,  \
                             tr_y_y_xxyyz,  \
                             tr_y_y_xxyz,   \
                             tr_y_y_xxyzz,  \
                             tr_y_y_xxzz,   \
                             tr_y_y_xxzzz,  \
                             tr_y_y_xyyy,   \
                             tr_y_y_xyyyy,  \
                             tr_y_y_xyyyz,  \
                             tr_y_y_xyyz,   \
                             tr_y_y_xyyzz,  \
                             tr_y_y_xyzz,   \
                             tr_y_y_xyzzz,  \
                             tr_y_y_xzzz,   \
                             tr_y_y_xzzzz,  \
                             tr_y_y_yyyy,   \
                             tr_y_y_yyyyy,  \
                             tr_y_y_yyyyz,  \
                             tr_y_y_yyyz,   \
                             tr_y_y_yyyzz,  \
                             tr_y_y_yyzz,   \
                             tr_y_y_yyzzz,  \
                             tr_y_y_yzzz,   \
                             tr_y_y_yzzzz,  \
                             tr_y_y_zzzz,   \
                             tr_y_y_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_xxxxx[i] = 5.0 * tr_y_y_xxxx[i] * fe_0 + tr_y_y_xxxxx[i] * pa_x[i];

        tr_y_xy_xxxxy[i] = 4.0 * tr_y_y_xxxy[i] * fe_0 + tr_y_y_xxxxy[i] * pa_x[i];

        tr_y_xy_xxxxz[i] = 4.0 * tr_y_y_xxxz[i] * fe_0 + tr_y_y_xxxxz[i] * pa_x[i];

        tr_y_xy_xxxyy[i] = 3.0 * tr_y_y_xxyy[i] * fe_0 + tr_y_y_xxxyy[i] * pa_x[i];

        tr_y_xy_xxxyz[i] = 3.0 * tr_y_y_xxyz[i] * fe_0 + tr_y_y_xxxyz[i] * pa_x[i];

        tr_y_xy_xxxzz[i] = 3.0 * tr_y_y_xxzz[i] * fe_0 + tr_y_y_xxxzz[i] * pa_x[i];

        tr_y_xy_xxyyy[i] = 2.0 * tr_y_y_xyyy[i] * fe_0 + tr_y_y_xxyyy[i] * pa_x[i];

        tr_y_xy_xxyyz[i] = 2.0 * tr_y_y_xyyz[i] * fe_0 + tr_y_y_xxyyz[i] * pa_x[i];

        tr_y_xy_xxyzz[i] = 2.0 * tr_y_y_xyzz[i] * fe_0 + tr_y_y_xxyzz[i] * pa_x[i];

        tr_y_xy_xxzzz[i] = 2.0 * tr_y_y_xzzz[i] * fe_0 + tr_y_y_xxzzz[i] * pa_x[i];

        tr_y_xy_xyyyy[i] = tr_y_y_yyyy[i] * fe_0 + tr_y_y_xyyyy[i] * pa_x[i];

        tr_y_xy_xyyyz[i] = tr_y_y_yyyz[i] * fe_0 + tr_y_y_xyyyz[i] * pa_x[i];

        tr_y_xy_xyyzz[i] = tr_y_y_yyzz[i] * fe_0 + tr_y_y_xyyzz[i] * pa_x[i];

        tr_y_xy_xyzzz[i] = tr_y_y_yzzz[i] * fe_0 + tr_y_y_xyzzz[i] * pa_x[i];

        tr_y_xy_xzzzz[i] = tr_y_y_zzzz[i] * fe_0 + tr_y_y_xzzzz[i] * pa_x[i];

        tr_y_xy_yyyyy[i] = tr_y_y_yyyyy[i] * pa_x[i];

        tr_y_xy_yyyyz[i] = tr_y_y_yyyyz[i] * pa_x[i];

        tr_y_xy_yyyzz[i] = tr_y_y_yyyzz[i] * pa_x[i];

        tr_y_xy_yyzzz[i] = tr_y_y_yyzzz[i] * pa_x[i];

        tr_y_xy_yzzzz[i] = tr_y_y_yzzzz[i] * pa_x[i];

        tr_y_xy_zzzzz[i] = tr_y_y_zzzzz[i] * pa_x[i];
    }

    // Set up 168-189 components of targeted buffer : DH

    auto tr_y_xz_xxxxx = pbuffer.data(idx_dip_dh + 168);

    auto tr_y_xz_xxxxy = pbuffer.data(idx_dip_dh + 169);

    auto tr_y_xz_xxxxz = pbuffer.data(idx_dip_dh + 170);

    auto tr_y_xz_xxxyy = pbuffer.data(idx_dip_dh + 171);

    auto tr_y_xz_xxxyz = pbuffer.data(idx_dip_dh + 172);

    auto tr_y_xz_xxxzz = pbuffer.data(idx_dip_dh + 173);

    auto tr_y_xz_xxyyy = pbuffer.data(idx_dip_dh + 174);

    auto tr_y_xz_xxyyz = pbuffer.data(idx_dip_dh + 175);

    auto tr_y_xz_xxyzz = pbuffer.data(idx_dip_dh + 176);

    auto tr_y_xz_xxzzz = pbuffer.data(idx_dip_dh + 177);

    auto tr_y_xz_xyyyy = pbuffer.data(idx_dip_dh + 178);

    auto tr_y_xz_xyyyz = pbuffer.data(idx_dip_dh + 179);

    auto tr_y_xz_xyyzz = pbuffer.data(idx_dip_dh + 180);

    auto tr_y_xz_xyzzz = pbuffer.data(idx_dip_dh + 181);

    auto tr_y_xz_xzzzz = pbuffer.data(idx_dip_dh + 182);

    auto tr_y_xz_yyyyy = pbuffer.data(idx_dip_dh + 183);

    auto tr_y_xz_yyyyz = pbuffer.data(idx_dip_dh + 184);

    auto tr_y_xz_yyyzz = pbuffer.data(idx_dip_dh + 185);

    auto tr_y_xz_yyzzz = pbuffer.data(idx_dip_dh + 186);

    auto tr_y_xz_yzzzz = pbuffer.data(idx_dip_dh + 187);

    auto tr_y_xz_zzzzz = pbuffer.data(idx_dip_dh + 188);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_x_xxxxx,  \
                             tr_y_x_xxxxy,  \
                             tr_y_x_xxxyy,  \
                             tr_y_x_xxyyy,  \
                             tr_y_x_xyyyy,  \
                             tr_y_xz_xxxxx, \
                             tr_y_xz_xxxxy, \
                             tr_y_xz_xxxxz, \
                             tr_y_xz_xxxyy, \
                             tr_y_xz_xxxyz, \
                             tr_y_xz_xxxzz, \
                             tr_y_xz_xxyyy, \
                             tr_y_xz_xxyyz, \
                             tr_y_xz_xxyzz, \
                             tr_y_xz_xxzzz, \
                             tr_y_xz_xyyyy, \
                             tr_y_xz_xyyyz, \
                             tr_y_xz_xyyzz, \
                             tr_y_xz_xyzzz, \
                             tr_y_xz_xzzzz, \
                             tr_y_xz_yyyyy, \
                             tr_y_xz_yyyyz, \
                             tr_y_xz_yyyzz, \
                             tr_y_xz_yyzzz, \
                             tr_y_xz_yzzzz, \
                             tr_y_xz_zzzzz, \
                             tr_y_z_xxxxz,  \
                             tr_y_z_xxxyz,  \
                             tr_y_z_xxxz,   \
                             tr_y_z_xxxzz,  \
                             tr_y_z_xxyyz,  \
                             tr_y_z_xxyz,   \
                             tr_y_z_xxyzz,  \
                             tr_y_z_xxzz,   \
                             tr_y_z_xxzzz,  \
                             tr_y_z_xyyyz,  \
                             tr_y_z_xyyz,   \
                             tr_y_z_xyyzz,  \
                             tr_y_z_xyzz,   \
                             tr_y_z_xyzzz,  \
                             tr_y_z_xzzz,   \
                             tr_y_z_xzzzz,  \
                             tr_y_z_yyyyy,  \
                             tr_y_z_yyyyz,  \
                             tr_y_z_yyyz,   \
                             tr_y_z_yyyzz,  \
                             tr_y_z_yyzz,   \
                             tr_y_z_yyzzz,  \
                             tr_y_z_yzzz,   \
                             tr_y_z_yzzzz,  \
                             tr_y_z_zzzz,   \
                             tr_y_z_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xz_xxxxx[i] = tr_y_x_xxxxx[i] * pa_z[i];

        tr_y_xz_xxxxy[i] = tr_y_x_xxxxy[i] * pa_z[i];

        tr_y_xz_xxxxz[i] = 4.0 * tr_y_z_xxxz[i] * fe_0 + tr_y_z_xxxxz[i] * pa_x[i];

        tr_y_xz_xxxyy[i] = tr_y_x_xxxyy[i] * pa_z[i];

        tr_y_xz_xxxyz[i] = 3.0 * tr_y_z_xxyz[i] * fe_0 + tr_y_z_xxxyz[i] * pa_x[i];

        tr_y_xz_xxxzz[i] = 3.0 * tr_y_z_xxzz[i] * fe_0 + tr_y_z_xxxzz[i] * pa_x[i];

        tr_y_xz_xxyyy[i] = tr_y_x_xxyyy[i] * pa_z[i];

        tr_y_xz_xxyyz[i] = 2.0 * tr_y_z_xyyz[i] * fe_0 + tr_y_z_xxyyz[i] * pa_x[i];

        tr_y_xz_xxyzz[i] = 2.0 * tr_y_z_xyzz[i] * fe_0 + tr_y_z_xxyzz[i] * pa_x[i];

        tr_y_xz_xxzzz[i] = 2.0 * tr_y_z_xzzz[i] * fe_0 + tr_y_z_xxzzz[i] * pa_x[i];

        tr_y_xz_xyyyy[i] = tr_y_x_xyyyy[i] * pa_z[i];

        tr_y_xz_xyyyz[i] = tr_y_z_yyyz[i] * fe_0 + tr_y_z_xyyyz[i] * pa_x[i];

        tr_y_xz_xyyzz[i] = tr_y_z_yyzz[i] * fe_0 + tr_y_z_xyyzz[i] * pa_x[i];

        tr_y_xz_xyzzz[i] = tr_y_z_yzzz[i] * fe_0 + tr_y_z_xyzzz[i] * pa_x[i];

        tr_y_xz_xzzzz[i] = tr_y_z_zzzz[i] * fe_0 + tr_y_z_xzzzz[i] * pa_x[i];

        tr_y_xz_yyyyy[i] = tr_y_z_yyyyy[i] * pa_x[i];

        tr_y_xz_yyyyz[i] = tr_y_z_yyyyz[i] * pa_x[i];

        tr_y_xz_yyyzz[i] = tr_y_z_yyyzz[i] * pa_x[i];

        tr_y_xz_yyzzz[i] = tr_y_z_yyzzz[i] * pa_x[i];

        tr_y_xz_yzzzz[i] = tr_y_z_yzzzz[i] * pa_x[i];

        tr_y_xz_zzzzz[i] = tr_y_z_zzzzz[i] * pa_x[i];
    }

    // Set up 189-210 components of targeted buffer : DH

    auto tr_y_yy_xxxxx = pbuffer.data(idx_dip_dh + 189);

    auto tr_y_yy_xxxxy = pbuffer.data(idx_dip_dh + 190);

    auto tr_y_yy_xxxxz = pbuffer.data(idx_dip_dh + 191);

    auto tr_y_yy_xxxyy = pbuffer.data(idx_dip_dh + 192);

    auto tr_y_yy_xxxyz = pbuffer.data(idx_dip_dh + 193);

    auto tr_y_yy_xxxzz = pbuffer.data(idx_dip_dh + 194);

    auto tr_y_yy_xxyyy = pbuffer.data(idx_dip_dh + 195);

    auto tr_y_yy_xxyyz = pbuffer.data(idx_dip_dh + 196);

    auto tr_y_yy_xxyzz = pbuffer.data(idx_dip_dh + 197);

    auto tr_y_yy_xxzzz = pbuffer.data(idx_dip_dh + 198);

    auto tr_y_yy_xyyyy = pbuffer.data(idx_dip_dh + 199);

    auto tr_y_yy_xyyyz = pbuffer.data(idx_dip_dh + 200);

    auto tr_y_yy_xyyzz = pbuffer.data(idx_dip_dh + 201);

    auto tr_y_yy_xyzzz = pbuffer.data(idx_dip_dh + 202);

    auto tr_y_yy_xzzzz = pbuffer.data(idx_dip_dh + 203);

    auto tr_y_yy_yyyyy = pbuffer.data(idx_dip_dh + 204);

    auto tr_y_yy_yyyyz = pbuffer.data(idx_dip_dh + 205);

    auto tr_y_yy_yyyzz = pbuffer.data(idx_dip_dh + 206);

    auto tr_y_yy_yyzzz = pbuffer.data(idx_dip_dh + 207);

    auto tr_y_yy_yzzzz = pbuffer.data(idx_dip_dh + 208);

    auto tr_y_yy_zzzzz = pbuffer.data(idx_dip_dh + 209);

#pragma omp simd aligned(pa_y,              \
                             tr_y_0_xxxxx,  \
                             tr_y_0_xxxxy,  \
                             tr_y_0_xxxxz,  \
                             tr_y_0_xxxyy,  \
                             tr_y_0_xxxyz,  \
                             tr_y_0_xxxzz,  \
                             tr_y_0_xxyyy,  \
                             tr_y_0_xxyyz,  \
                             tr_y_0_xxyzz,  \
                             tr_y_0_xxzzz,  \
                             tr_y_0_xyyyy,  \
                             tr_y_0_xyyyz,  \
                             tr_y_0_xyyzz,  \
                             tr_y_0_xyzzz,  \
                             tr_y_0_xzzzz,  \
                             tr_y_0_yyyyy,  \
                             tr_y_0_yyyyz,  \
                             tr_y_0_yyyzz,  \
                             tr_y_0_yyzzz,  \
                             tr_y_0_yzzzz,  \
                             tr_y_0_zzzzz,  \
                             tr_y_y_xxxx,   \
                             tr_y_y_xxxxx,  \
                             tr_y_y_xxxxy,  \
                             tr_y_y_xxxxz,  \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxxyy,  \
                             tr_y_y_xxxyz,  \
                             tr_y_y_xxxz,   \
                             tr_y_y_xxxzz,  \
                             tr_y_y_xxyy,   \
                             tr_y_y_xxyyy,  \
                             tr_y_y_xxyyz,  \
                             tr_y_y_xxyz,   \
                             tr_y_y_xxyzz,  \
                             tr_y_y_xxzz,   \
                             tr_y_y_xxzzz,  \
                             tr_y_y_xyyy,   \
                             tr_y_y_xyyyy,  \
                             tr_y_y_xyyyz,  \
                             tr_y_y_xyyz,   \
                             tr_y_y_xyyzz,  \
                             tr_y_y_xyzz,   \
                             tr_y_y_xyzzz,  \
                             tr_y_y_xzzz,   \
                             tr_y_y_xzzzz,  \
                             tr_y_y_yyyy,   \
                             tr_y_y_yyyyy,  \
                             tr_y_y_yyyyz,  \
                             tr_y_y_yyyz,   \
                             tr_y_y_yyyzz,  \
                             tr_y_y_yyzz,   \
                             tr_y_y_yyzzz,  \
                             tr_y_y_yzzz,   \
                             tr_y_y_yzzzz,  \
                             tr_y_y_zzzz,   \
                             tr_y_y_zzzzz,  \
                             tr_y_yy_xxxxx, \
                             tr_y_yy_xxxxy, \
                             tr_y_yy_xxxxz, \
                             tr_y_yy_xxxyy, \
                             tr_y_yy_xxxyz, \
                             tr_y_yy_xxxzz, \
                             tr_y_yy_xxyyy, \
                             tr_y_yy_xxyyz, \
                             tr_y_yy_xxyzz, \
                             tr_y_yy_xxzzz, \
                             tr_y_yy_xyyyy, \
                             tr_y_yy_xyyyz, \
                             tr_y_yy_xyyzz, \
                             tr_y_yy_xyzzz, \
                             tr_y_yy_xzzzz, \
                             tr_y_yy_yyyyy, \
                             tr_y_yy_yyyyz, \
                             tr_y_yy_yyyzz, \
                             tr_y_yy_yyzzz, \
                             tr_y_yy_yzzzz, \
                             tr_y_yy_zzzzz, \
                             ts_y_xxxxx,    \
                             ts_y_xxxxy,    \
                             ts_y_xxxxz,    \
                             ts_y_xxxyy,    \
                             ts_y_xxxyz,    \
                             ts_y_xxxzz,    \
                             ts_y_xxyyy,    \
                             ts_y_xxyyz,    \
                             ts_y_xxyzz,    \
                             ts_y_xxzzz,    \
                             ts_y_xyyyy,    \
                             ts_y_xyyyz,    \
                             ts_y_xyyzz,    \
                             ts_y_xyzzz,    \
                             ts_y_xzzzz,    \
                             ts_y_yyyyy,    \
                             ts_y_yyyyz,    \
                             ts_y_yyyzz,    \
                             ts_y_yyzzz,    \
                             ts_y_yzzzz,    \
                             ts_y_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_xxxxx[i] = tr_y_0_xxxxx[i] * fe_0 + ts_y_xxxxx[i] * fe_0 + tr_y_y_xxxxx[i] * pa_y[i];

        tr_y_yy_xxxxy[i] = tr_y_0_xxxxy[i] * fe_0 + tr_y_y_xxxx[i] * fe_0 + ts_y_xxxxy[i] * fe_0 + tr_y_y_xxxxy[i] * pa_y[i];

        tr_y_yy_xxxxz[i] = tr_y_0_xxxxz[i] * fe_0 + ts_y_xxxxz[i] * fe_0 + tr_y_y_xxxxz[i] * pa_y[i];

        tr_y_yy_xxxyy[i] = tr_y_0_xxxyy[i] * fe_0 + 2.0 * tr_y_y_xxxy[i] * fe_0 + ts_y_xxxyy[i] * fe_0 + tr_y_y_xxxyy[i] * pa_y[i];

        tr_y_yy_xxxyz[i] = tr_y_0_xxxyz[i] * fe_0 + tr_y_y_xxxz[i] * fe_0 + ts_y_xxxyz[i] * fe_0 + tr_y_y_xxxyz[i] * pa_y[i];

        tr_y_yy_xxxzz[i] = tr_y_0_xxxzz[i] * fe_0 + ts_y_xxxzz[i] * fe_0 + tr_y_y_xxxzz[i] * pa_y[i];

        tr_y_yy_xxyyy[i] = tr_y_0_xxyyy[i] * fe_0 + 3.0 * tr_y_y_xxyy[i] * fe_0 + ts_y_xxyyy[i] * fe_0 + tr_y_y_xxyyy[i] * pa_y[i];

        tr_y_yy_xxyyz[i] = tr_y_0_xxyyz[i] * fe_0 + 2.0 * tr_y_y_xxyz[i] * fe_0 + ts_y_xxyyz[i] * fe_0 + tr_y_y_xxyyz[i] * pa_y[i];

        tr_y_yy_xxyzz[i] = tr_y_0_xxyzz[i] * fe_0 + tr_y_y_xxzz[i] * fe_0 + ts_y_xxyzz[i] * fe_0 + tr_y_y_xxyzz[i] * pa_y[i];

        tr_y_yy_xxzzz[i] = tr_y_0_xxzzz[i] * fe_0 + ts_y_xxzzz[i] * fe_0 + tr_y_y_xxzzz[i] * pa_y[i];

        tr_y_yy_xyyyy[i] = tr_y_0_xyyyy[i] * fe_0 + 4.0 * tr_y_y_xyyy[i] * fe_0 + ts_y_xyyyy[i] * fe_0 + tr_y_y_xyyyy[i] * pa_y[i];

        tr_y_yy_xyyyz[i] = tr_y_0_xyyyz[i] * fe_0 + 3.0 * tr_y_y_xyyz[i] * fe_0 + ts_y_xyyyz[i] * fe_0 + tr_y_y_xyyyz[i] * pa_y[i];

        tr_y_yy_xyyzz[i] = tr_y_0_xyyzz[i] * fe_0 + 2.0 * tr_y_y_xyzz[i] * fe_0 + ts_y_xyyzz[i] * fe_0 + tr_y_y_xyyzz[i] * pa_y[i];

        tr_y_yy_xyzzz[i] = tr_y_0_xyzzz[i] * fe_0 + tr_y_y_xzzz[i] * fe_0 + ts_y_xyzzz[i] * fe_0 + tr_y_y_xyzzz[i] * pa_y[i];

        tr_y_yy_xzzzz[i] = tr_y_0_xzzzz[i] * fe_0 + ts_y_xzzzz[i] * fe_0 + tr_y_y_xzzzz[i] * pa_y[i];

        tr_y_yy_yyyyy[i] = tr_y_0_yyyyy[i] * fe_0 + 5.0 * tr_y_y_yyyy[i] * fe_0 + ts_y_yyyyy[i] * fe_0 + tr_y_y_yyyyy[i] * pa_y[i];

        tr_y_yy_yyyyz[i] = tr_y_0_yyyyz[i] * fe_0 + 4.0 * tr_y_y_yyyz[i] * fe_0 + ts_y_yyyyz[i] * fe_0 + tr_y_y_yyyyz[i] * pa_y[i];

        tr_y_yy_yyyzz[i] = tr_y_0_yyyzz[i] * fe_0 + 3.0 * tr_y_y_yyzz[i] * fe_0 + ts_y_yyyzz[i] * fe_0 + tr_y_y_yyyzz[i] * pa_y[i];

        tr_y_yy_yyzzz[i] = tr_y_0_yyzzz[i] * fe_0 + 2.0 * tr_y_y_yzzz[i] * fe_0 + ts_y_yyzzz[i] * fe_0 + tr_y_y_yyzzz[i] * pa_y[i];

        tr_y_yy_yzzzz[i] = tr_y_0_yzzzz[i] * fe_0 + tr_y_y_zzzz[i] * fe_0 + ts_y_yzzzz[i] * fe_0 + tr_y_y_yzzzz[i] * pa_y[i];

        tr_y_yy_zzzzz[i] = tr_y_0_zzzzz[i] * fe_0 + ts_y_zzzzz[i] * fe_0 + tr_y_y_zzzzz[i] * pa_y[i];
    }

    // Set up 210-231 components of targeted buffer : DH

    auto tr_y_yz_xxxxx = pbuffer.data(idx_dip_dh + 210);

    auto tr_y_yz_xxxxy = pbuffer.data(idx_dip_dh + 211);

    auto tr_y_yz_xxxxz = pbuffer.data(idx_dip_dh + 212);

    auto tr_y_yz_xxxyy = pbuffer.data(idx_dip_dh + 213);

    auto tr_y_yz_xxxyz = pbuffer.data(idx_dip_dh + 214);

    auto tr_y_yz_xxxzz = pbuffer.data(idx_dip_dh + 215);

    auto tr_y_yz_xxyyy = pbuffer.data(idx_dip_dh + 216);

    auto tr_y_yz_xxyyz = pbuffer.data(idx_dip_dh + 217);

    auto tr_y_yz_xxyzz = pbuffer.data(idx_dip_dh + 218);

    auto tr_y_yz_xxzzz = pbuffer.data(idx_dip_dh + 219);

    auto tr_y_yz_xyyyy = pbuffer.data(idx_dip_dh + 220);

    auto tr_y_yz_xyyyz = pbuffer.data(idx_dip_dh + 221);

    auto tr_y_yz_xyyzz = pbuffer.data(idx_dip_dh + 222);

    auto tr_y_yz_xyzzz = pbuffer.data(idx_dip_dh + 223);

    auto tr_y_yz_xzzzz = pbuffer.data(idx_dip_dh + 224);

    auto tr_y_yz_yyyyy = pbuffer.data(idx_dip_dh + 225);

    auto tr_y_yz_yyyyz = pbuffer.data(idx_dip_dh + 226);

    auto tr_y_yz_yyyzz = pbuffer.data(idx_dip_dh + 227);

    auto tr_y_yz_yyzzz = pbuffer.data(idx_dip_dh + 228);

    auto tr_y_yz_yzzzz = pbuffer.data(idx_dip_dh + 229);

    auto tr_y_yz_zzzzz = pbuffer.data(idx_dip_dh + 230);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_y_xxxxx,  \
                             tr_y_y_xxxxy,  \
                             tr_y_y_xxxy,   \
                             tr_y_y_xxxyy,  \
                             tr_y_y_xxxyz,  \
                             tr_y_y_xxyy,   \
                             tr_y_y_xxyyy,  \
                             tr_y_y_xxyyz,  \
                             tr_y_y_xxyz,   \
                             tr_y_y_xxyzz,  \
                             tr_y_y_xyyy,   \
                             tr_y_y_xyyyy,  \
                             tr_y_y_xyyyz,  \
                             tr_y_y_xyyz,   \
                             tr_y_y_xyyzz,  \
                             tr_y_y_xyzz,   \
                             tr_y_y_xyzzz,  \
                             tr_y_y_yyyy,   \
                             tr_y_y_yyyyy,  \
                             tr_y_y_yyyyz,  \
                             tr_y_y_yyyz,   \
                             tr_y_y_yyyzz,  \
                             tr_y_y_yyzz,   \
                             tr_y_y_yyzzz,  \
                             tr_y_y_yzzz,   \
                             tr_y_y_yzzzz,  \
                             tr_y_yz_xxxxx, \
                             tr_y_yz_xxxxy, \
                             tr_y_yz_xxxxz, \
                             tr_y_yz_xxxyy, \
                             tr_y_yz_xxxyz, \
                             tr_y_yz_xxxzz, \
                             tr_y_yz_xxyyy, \
                             tr_y_yz_xxyyz, \
                             tr_y_yz_xxyzz, \
                             tr_y_yz_xxzzz, \
                             tr_y_yz_xyyyy, \
                             tr_y_yz_xyyyz, \
                             tr_y_yz_xyyzz, \
                             tr_y_yz_xyzzz, \
                             tr_y_yz_xzzzz, \
                             tr_y_yz_yyyyy, \
                             tr_y_yz_yyyyz, \
                             tr_y_yz_yyyzz, \
                             tr_y_yz_yyzzz, \
                             tr_y_yz_yzzzz, \
                             tr_y_yz_zzzzz, \
                             tr_y_z_xxxxz,  \
                             tr_y_z_xxxzz,  \
                             tr_y_z_xxzzz,  \
                             tr_y_z_xzzzz,  \
                             tr_y_z_zzzzz,  \
                             ts_z_xxxxz,    \
                             ts_z_xxxzz,    \
                             ts_z_xxzzz,    \
                             ts_z_xzzzz,    \
                             ts_z_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_xxxxx[i] = tr_y_y_xxxxx[i] * pa_z[i];

        tr_y_yz_xxxxy[i] = tr_y_y_xxxxy[i] * pa_z[i];

        tr_y_yz_xxxxz[i] = ts_z_xxxxz[i] * fe_0 + tr_y_z_xxxxz[i] * pa_y[i];

        tr_y_yz_xxxyy[i] = tr_y_y_xxxyy[i] * pa_z[i];

        tr_y_yz_xxxyz[i] = tr_y_y_xxxy[i] * fe_0 + tr_y_y_xxxyz[i] * pa_z[i];

        tr_y_yz_xxxzz[i] = ts_z_xxxzz[i] * fe_0 + tr_y_z_xxxzz[i] * pa_y[i];

        tr_y_yz_xxyyy[i] = tr_y_y_xxyyy[i] * pa_z[i];

        tr_y_yz_xxyyz[i] = tr_y_y_xxyy[i] * fe_0 + tr_y_y_xxyyz[i] * pa_z[i];

        tr_y_yz_xxyzz[i] = 2.0 * tr_y_y_xxyz[i] * fe_0 + tr_y_y_xxyzz[i] * pa_z[i];

        tr_y_yz_xxzzz[i] = ts_z_xxzzz[i] * fe_0 + tr_y_z_xxzzz[i] * pa_y[i];

        tr_y_yz_xyyyy[i] = tr_y_y_xyyyy[i] * pa_z[i];

        tr_y_yz_xyyyz[i] = tr_y_y_xyyy[i] * fe_0 + tr_y_y_xyyyz[i] * pa_z[i];

        tr_y_yz_xyyzz[i] = 2.0 * tr_y_y_xyyz[i] * fe_0 + tr_y_y_xyyzz[i] * pa_z[i];

        tr_y_yz_xyzzz[i] = 3.0 * tr_y_y_xyzz[i] * fe_0 + tr_y_y_xyzzz[i] * pa_z[i];

        tr_y_yz_xzzzz[i] = ts_z_xzzzz[i] * fe_0 + tr_y_z_xzzzz[i] * pa_y[i];

        tr_y_yz_yyyyy[i] = tr_y_y_yyyyy[i] * pa_z[i];

        tr_y_yz_yyyyz[i] = tr_y_y_yyyy[i] * fe_0 + tr_y_y_yyyyz[i] * pa_z[i];

        tr_y_yz_yyyzz[i] = 2.0 * tr_y_y_yyyz[i] * fe_0 + tr_y_y_yyyzz[i] * pa_z[i];

        tr_y_yz_yyzzz[i] = 3.0 * tr_y_y_yyzz[i] * fe_0 + tr_y_y_yyzzz[i] * pa_z[i];

        tr_y_yz_yzzzz[i] = 4.0 * tr_y_y_yzzz[i] * fe_0 + tr_y_y_yzzzz[i] * pa_z[i];

        tr_y_yz_zzzzz[i] = ts_z_zzzzz[i] * fe_0 + tr_y_z_zzzzz[i] * pa_y[i];
    }

    // Set up 231-252 components of targeted buffer : DH

    auto tr_y_zz_xxxxx = pbuffer.data(idx_dip_dh + 231);

    auto tr_y_zz_xxxxy = pbuffer.data(idx_dip_dh + 232);

    auto tr_y_zz_xxxxz = pbuffer.data(idx_dip_dh + 233);

    auto tr_y_zz_xxxyy = pbuffer.data(idx_dip_dh + 234);

    auto tr_y_zz_xxxyz = pbuffer.data(idx_dip_dh + 235);

    auto tr_y_zz_xxxzz = pbuffer.data(idx_dip_dh + 236);

    auto tr_y_zz_xxyyy = pbuffer.data(idx_dip_dh + 237);

    auto tr_y_zz_xxyyz = pbuffer.data(idx_dip_dh + 238);

    auto tr_y_zz_xxyzz = pbuffer.data(idx_dip_dh + 239);

    auto tr_y_zz_xxzzz = pbuffer.data(idx_dip_dh + 240);

    auto tr_y_zz_xyyyy = pbuffer.data(idx_dip_dh + 241);

    auto tr_y_zz_xyyyz = pbuffer.data(idx_dip_dh + 242);

    auto tr_y_zz_xyyzz = pbuffer.data(idx_dip_dh + 243);

    auto tr_y_zz_xyzzz = pbuffer.data(idx_dip_dh + 244);

    auto tr_y_zz_xzzzz = pbuffer.data(idx_dip_dh + 245);

    auto tr_y_zz_yyyyy = pbuffer.data(idx_dip_dh + 246);

    auto tr_y_zz_yyyyz = pbuffer.data(idx_dip_dh + 247);

    auto tr_y_zz_yyyzz = pbuffer.data(idx_dip_dh + 248);

    auto tr_y_zz_yyzzz = pbuffer.data(idx_dip_dh + 249);

    auto tr_y_zz_yzzzz = pbuffer.data(idx_dip_dh + 250);

    auto tr_y_zz_zzzzz = pbuffer.data(idx_dip_dh + 251);

#pragma omp simd aligned(pa_z,              \
                             tr_y_0_xxxxx,  \
                             tr_y_0_xxxxy,  \
                             tr_y_0_xxxxz,  \
                             tr_y_0_xxxyy,  \
                             tr_y_0_xxxyz,  \
                             tr_y_0_xxxzz,  \
                             tr_y_0_xxyyy,  \
                             tr_y_0_xxyyz,  \
                             tr_y_0_xxyzz,  \
                             tr_y_0_xxzzz,  \
                             tr_y_0_xyyyy,  \
                             tr_y_0_xyyyz,  \
                             tr_y_0_xyyzz,  \
                             tr_y_0_xyzzz,  \
                             tr_y_0_xzzzz,  \
                             tr_y_0_yyyyy,  \
                             tr_y_0_yyyyz,  \
                             tr_y_0_yyyzz,  \
                             tr_y_0_yyzzz,  \
                             tr_y_0_yzzzz,  \
                             tr_y_0_zzzzz,  \
                             tr_y_z_xxxx,   \
                             tr_y_z_xxxxx,  \
                             tr_y_z_xxxxy,  \
                             tr_y_z_xxxxz,  \
                             tr_y_z_xxxy,   \
                             tr_y_z_xxxyy,  \
                             tr_y_z_xxxyz,  \
                             tr_y_z_xxxz,   \
                             tr_y_z_xxxzz,  \
                             tr_y_z_xxyy,   \
                             tr_y_z_xxyyy,  \
                             tr_y_z_xxyyz,  \
                             tr_y_z_xxyz,   \
                             tr_y_z_xxyzz,  \
                             tr_y_z_xxzz,   \
                             tr_y_z_xxzzz,  \
                             tr_y_z_xyyy,   \
                             tr_y_z_xyyyy,  \
                             tr_y_z_xyyyz,  \
                             tr_y_z_xyyz,   \
                             tr_y_z_xyyzz,  \
                             tr_y_z_xyzz,   \
                             tr_y_z_xyzzz,  \
                             tr_y_z_xzzz,   \
                             tr_y_z_xzzzz,  \
                             tr_y_z_yyyy,   \
                             tr_y_z_yyyyy,  \
                             tr_y_z_yyyyz,  \
                             tr_y_z_yyyz,   \
                             tr_y_z_yyyzz,  \
                             tr_y_z_yyzz,   \
                             tr_y_z_yyzzz,  \
                             tr_y_z_yzzz,   \
                             tr_y_z_yzzzz,  \
                             tr_y_z_zzzz,   \
                             tr_y_z_zzzzz,  \
                             tr_y_zz_xxxxx, \
                             tr_y_zz_xxxxy, \
                             tr_y_zz_xxxxz, \
                             tr_y_zz_xxxyy, \
                             tr_y_zz_xxxyz, \
                             tr_y_zz_xxxzz, \
                             tr_y_zz_xxyyy, \
                             tr_y_zz_xxyyz, \
                             tr_y_zz_xxyzz, \
                             tr_y_zz_xxzzz, \
                             tr_y_zz_xyyyy, \
                             tr_y_zz_xyyyz, \
                             tr_y_zz_xyyzz, \
                             tr_y_zz_xyzzz, \
                             tr_y_zz_xzzzz, \
                             tr_y_zz_yyyyy, \
                             tr_y_zz_yyyyz, \
                             tr_y_zz_yyyzz, \
                             tr_y_zz_yyzzz, \
                             tr_y_zz_yzzzz, \
                             tr_y_zz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_xxxxx[i] = tr_y_0_xxxxx[i] * fe_0 + tr_y_z_xxxxx[i] * pa_z[i];

        tr_y_zz_xxxxy[i] = tr_y_0_xxxxy[i] * fe_0 + tr_y_z_xxxxy[i] * pa_z[i];

        tr_y_zz_xxxxz[i] = tr_y_0_xxxxz[i] * fe_0 + tr_y_z_xxxx[i] * fe_0 + tr_y_z_xxxxz[i] * pa_z[i];

        tr_y_zz_xxxyy[i] = tr_y_0_xxxyy[i] * fe_0 + tr_y_z_xxxyy[i] * pa_z[i];

        tr_y_zz_xxxyz[i] = tr_y_0_xxxyz[i] * fe_0 + tr_y_z_xxxy[i] * fe_0 + tr_y_z_xxxyz[i] * pa_z[i];

        tr_y_zz_xxxzz[i] = tr_y_0_xxxzz[i] * fe_0 + 2.0 * tr_y_z_xxxz[i] * fe_0 + tr_y_z_xxxzz[i] * pa_z[i];

        tr_y_zz_xxyyy[i] = tr_y_0_xxyyy[i] * fe_0 + tr_y_z_xxyyy[i] * pa_z[i];

        tr_y_zz_xxyyz[i] = tr_y_0_xxyyz[i] * fe_0 + tr_y_z_xxyy[i] * fe_0 + tr_y_z_xxyyz[i] * pa_z[i];

        tr_y_zz_xxyzz[i] = tr_y_0_xxyzz[i] * fe_0 + 2.0 * tr_y_z_xxyz[i] * fe_0 + tr_y_z_xxyzz[i] * pa_z[i];

        tr_y_zz_xxzzz[i] = tr_y_0_xxzzz[i] * fe_0 + 3.0 * tr_y_z_xxzz[i] * fe_0 + tr_y_z_xxzzz[i] * pa_z[i];

        tr_y_zz_xyyyy[i] = tr_y_0_xyyyy[i] * fe_0 + tr_y_z_xyyyy[i] * pa_z[i];

        tr_y_zz_xyyyz[i] = tr_y_0_xyyyz[i] * fe_0 + tr_y_z_xyyy[i] * fe_0 + tr_y_z_xyyyz[i] * pa_z[i];

        tr_y_zz_xyyzz[i] = tr_y_0_xyyzz[i] * fe_0 + 2.0 * tr_y_z_xyyz[i] * fe_0 + tr_y_z_xyyzz[i] * pa_z[i];

        tr_y_zz_xyzzz[i] = tr_y_0_xyzzz[i] * fe_0 + 3.0 * tr_y_z_xyzz[i] * fe_0 + tr_y_z_xyzzz[i] * pa_z[i];

        tr_y_zz_xzzzz[i] = tr_y_0_xzzzz[i] * fe_0 + 4.0 * tr_y_z_xzzz[i] * fe_0 + tr_y_z_xzzzz[i] * pa_z[i];

        tr_y_zz_yyyyy[i] = tr_y_0_yyyyy[i] * fe_0 + tr_y_z_yyyyy[i] * pa_z[i];

        tr_y_zz_yyyyz[i] = tr_y_0_yyyyz[i] * fe_0 + tr_y_z_yyyy[i] * fe_0 + tr_y_z_yyyyz[i] * pa_z[i];

        tr_y_zz_yyyzz[i] = tr_y_0_yyyzz[i] * fe_0 + 2.0 * tr_y_z_yyyz[i] * fe_0 + tr_y_z_yyyzz[i] * pa_z[i];

        tr_y_zz_yyzzz[i] = tr_y_0_yyzzz[i] * fe_0 + 3.0 * tr_y_z_yyzz[i] * fe_0 + tr_y_z_yyzzz[i] * pa_z[i];

        tr_y_zz_yzzzz[i] = tr_y_0_yzzzz[i] * fe_0 + 4.0 * tr_y_z_yzzz[i] * fe_0 + tr_y_z_yzzzz[i] * pa_z[i];

        tr_y_zz_zzzzz[i] = tr_y_0_zzzzz[i] * fe_0 + 5.0 * tr_y_z_zzzz[i] * fe_0 + tr_y_z_zzzzz[i] * pa_z[i];
    }

    // Set up 252-273 components of targeted buffer : DH

    auto tr_z_xx_xxxxx = pbuffer.data(idx_dip_dh + 252);

    auto tr_z_xx_xxxxy = pbuffer.data(idx_dip_dh + 253);

    auto tr_z_xx_xxxxz = pbuffer.data(idx_dip_dh + 254);

    auto tr_z_xx_xxxyy = pbuffer.data(idx_dip_dh + 255);

    auto tr_z_xx_xxxyz = pbuffer.data(idx_dip_dh + 256);

    auto tr_z_xx_xxxzz = pbuffer.data(idx_dip_dh + 257);

    auto tr_z_xx_xxyyy = pbuffer.data(idx_dip_dh + 258);

    auto tr_z_xx_xxyyz = pbuffer.data(idx_dip_dh + 259);

    auto tr_z_xx_xxyzz = pbuffer.data(idx_dip_dh + 260);

    auto tr_z_xx_xxzzz = pbuffer.data(idx_dip_dh + 261);

    auto tr_z_xx_xyyyy = pbuffer.data(idx_dip_dh + 262);

    auto tr_z_xx_xyyyz = pbuffer.data(idx_dip_dh + 263);

    auto tr_z_xx_xyyzz = pbuffer.data(idx_dip_dh + 264);

    auto tr_z_xx_xyzzz = pbuffer.data(idx_dip_dh + 265);

    auto tr_z_xx_xzzzz = pbuffer.data(idx_dip_dh + 266);

    auto tr_z_xx_yyyyy = pbuffer.data(idx_dip_dh + 267);

    auto tr_z_xx_yyyyz = pbuffer.data(idx_dip_dh + 268);

    auto tr_z_xx_yyyzz = pbuffer.data(idx_dip_dh + 269);

    auto tr_z_xx_yyzzz = pbuffer.data(idx_dip_dh + 270);

    auto tr_z_xx_yzzzz = pbuffer.data(idx_dip_dh + 271);

    auto tr_z_xx_zzzzz = pbuffer.data(idx_dip_dh + 272);

#pragma omp simd aligned(pa_x,              \
                             tr_z_0_xxxxx,  \
                             tr_z_0_xxxxy,  \
                             tr_z_0_xxxxz,  \
                             tr_z_0_xxxyy,  \
                             tr_z_0_xxxyz,  \
                             tr_z_0_xxxzz,  \
                             tr_z_0_xxyyy,  \
                             tr_z_0_xxyyz,  \
                             tr_z_0_xxyzz,  \
                             tr_z_0_xxzzz,  \
                             tr_z_0_xyyyy,  \
                             tr_z_0_xyyyz,  \
                             tr_z_0_xyyzz,  \
                             tr_z_0_xyzzz,  \
                             tr_z_0_xzzzz,  \
                             tr_z_0_yyyyy,  \
                             tr_z_0_yyyyz,  \
                             tr_z_0_yyyzz,  \
                             tr_z_0_yyzzz,  \
                             tr_z_0_yzzzz,  \
                             tr_z_0_zzzzz,  \
                             tr_z_x_xxxx,   \
                             tr_z_x_xxxxx,  \
                             tr_z_x_xxxxy,  \
                             tr_z_x_xxxxz,  \
                             tr_z_x_xxxy,   \
                             tr_z_x_xxxyy,  \
                             tr_z_x_xxxyz,  \
                             tr_z_x_xxxz,   \
                             tr_z_x_xxxzz,  \
                             tr_z_x_xxyy,   \
                             tr_z_x_xxyyy,  \
                             tr_z_x_xxyyz,  \
                             tr_z_x_xxyz,   \
                             tr_z_x_xxyzz,  \
                             tr_z_x_xxzz,   \
                             tr_z_x_xxzzz,  \
                             tr_z_x_xyyy,   \
                             tr_z_x_xyyyy,  \
                             tr_z_x_xyyyz,  \
                             tr_z_x_xyyz,   \
                             tr_z_x_xyyzz,  \
                             tr_z_x_xyzz,   \
                             tr_z_x_xyzzz,  \
                             tr_z_x_xzzz,   \
                             tr_z_x_xzzzz,  \
                             tr_z_x_yyyy,   \
                             tr_z_x_yyyyy,  \
                             tr_z_x_yyyyz,  \
                             tr_z_x_yyyz,   \
                             tr_z_x_yyyzz,  \
                             tr_z_x_yyzz,   \
                             tr_z_x_yyzzz,  \
                             tr_z_x_yzzz,   \
                             tr_z_x_yzzzz,  \
                             tr_z_x_zzzz,   \
                             tr_z_x_zzzzz,  \
                             tr_z_xx_xxxxx, \
                             tr_z_xx_xxxxy, \
                             tr_z_xx_xxxxz, \
                             tr_z_xx_xxxyy, \
                             tr_z_xx_xxxyz, \
                             tr_z_xx_xxxzz, \
                             tr_z_xx_xxyyy, \
                             tr_z_xx_xxyyz, \
                             tr_z_xx_xxyzz, \
                             tr_z_xx_xxzzz, \
                             tr_z_xx_xyyyy, \
                             tr_z_xx_xyyyz, \
                             tr_z_xx_xyyzz, \
                             tr_z_xx_xyzzz, \
                             tr_z_xx_xzzzz, \
                             tr_z_xx_yyyyy, \
                             tr_z_xx_yyyyz, \
                             tr_z_xx_yyyzz, \
                             tr_z_xx_yyzzz, \
                             tr_z_xx_yzzzz, \
                             tr_z_xx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_xxxxx[i] = tr_z_0_xxxxx[i] * fe_0 + 5.0 * tr_z_x_xxxx[i] * fe_0 + tr_z_x_xxxxx[i] * pa_x[i];

        tr_z_xx_xxxxy[i] = tr_z_0_xxxxy[i] * fe_0 + 4.0 * tr_z_x_xxxy[i] * fe_0 + tr_z_x_xxxxy[i] * pa_x[i];

        tr_z_xx_xxxxz[i] = tr_z_0_xxxxz[i] * fe_0 + 4.0 * tr_z_x_xxxz[i] * fe_0 + tr_z_x_xxxxz[i] * pa_x[i];

        tr_z_xx_xxxyy[i] = tr_z_0_xxxyy[i] * fe_0 + 3.0 * tr_z_x_xxyy[i] * fe_0 + tr_z_x_xxxyy[i] * pa_x[i];

        tr_z_xx_xxxyz[i] = tr_z_0_xxxyz[i] * fe_0 + 3.0 * tr_z_x_xxyz[i] * fe_0 + tr_z_x_xxxyz[i] * pa_x[i];

        tr_z_xx_xxxzz[i] = tr_z_0_xxxzz[i] * fe_0 + 3.0 * tr_z_x_xxzz[i] * fe_0 + tr_z_x_xxxzz[i] * pa_x[i];

        tr_z_xx_xxyyy[i] = tr_z_0_xxyyy[i] * fe_0 + 2.0 * tr_z_x_xyyy[i] * fe_0 + tr_z_x_xxyyy[i] * pa_x[i];

        tr_z_xx_xxyyz[i] = tr_z_0_xxyyz[i] * fe_0 + 2.0 * tr_z_x_xyyz[i] * fe_0 + tr_z_x_xxyyz[i] * pa_x[i];

        tr_z_xx_xxyzz[i] = tr_z_0_xxyzz[i] * fe_0 + 2.0 * tr_z_x_xyzz[i] * fe_0 + tr_z_x_xxyzz[i] * pa_x[i];

        tr_z_xx_xxzzz[i] = tr_z_0_xxzzz[i] * fe_0 + 2.0 * tr_z_x_xzzz[i] * fe_0 + tr_z_x_xxzzz[i] * pa_x[i];

        tr_z_xx_xyyyy[i] = tr_z_0_xyyyy[i] * fe_0 + tr_z_x_yyyy[i] * fe_0 + tr_z_x_xyyyy[i] * pa_x[i];

        tr_z_xx_xyyyz[i] = tr_z_0_xyyyz[i] * fe_0 + tr_z_x_yyyz[i] * fe_0 + tr_z_x_xyyyz[i] * pa_x[i];

        tr_z_xx_xyyzz[i] = tr_z_0_xyyzz[i] * fe_0 + tr_z_x_yyzz[i] * fe_0 + tr_z_x_xyyzz[i] * pa_x[i];

        tr_z_xx_xyzzz[i] = tr_z_0_xyzzz[i] * fe_0 + tr_z_x_yzzz[i] * fe_0 + tr_z_x_xyzzz[i] * pa_x[i];

        tr_z_xx_xzzzz[i] = tr_z_0_xzzzz[i] * fe_0 + tr_z_x_zzzz[i] * fe_0 + tr_z_x_xzzzz[i] * pa_x[i];

        tr_z_xx_yyyyy[i] = tr_z_0_yyyyy[i] * fe_0 + tr_z_x_yyyyy[i] * pa_x[i];

        tr_z_xx_yyyyz[i] = tr_z_0_yyyyz[i] * fe_0 + tr_z_x_yyyyz[i] * pa_x[i];

        tr_z_xx_yyyzz[i] = tr_z_0_yyyzz[i] * fe_0 + tr_z_x_yyyzz[i] * pa_x[i];

        tr_z_xx_yyzzz[i] = tr_z_0_yyzzz[i] * fe_0 + tr_z_x_yyzzz[i] * pa_x[i];

        tr_z_xx_yzzzz[i] = tr_z_0_yzzzz[i] * fe_0 + tr_z_x_yzzzz[i] * pa_x[i];

        tr_z_xx_zzzzz[i] = tr_z_0_zzzzz[i] * fe_0 + tr_z_x_zzzzz[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : DH

    auto tr_z_xy_xxxxx = pbuffer.data(idx_dip_dh + 273);

    auto tr_z_xy_xxxxy = pbuffer.data(idx_dip_dh + 274);

    auto tr_z_xy_xxxxz = pbuffer.data(idx_dip_dh + 275);

    auto tr_z_xy_xxxyy = pbuffer.data(idx_dip_dh + 276);

    auto tr_z_xy_xxxyz = pbuffer.data(idx_dip_dh + 277);

    auto tr_z_xy_xxxzz = pbuffer.data(idx_dip_dh + 278);

    auto tr_z_xy_xxyyy = pbuffer.data(idx_dip_dh + 279);

    auto tr_z_xy_xxyyz = pbuffer.data(idx_dip_dh + 280);

    auto tr_z_xy_xxyzz = pbuffer.data(idx_dip_dh + 281);

    auto tr_z_xy_xxzzz = pbuffer.data(idx_dip_dh + 282);

    auto tr_z_xy_xyyyy = pbuffer.data(idx_dip_dh + 283);

    auto tr_z_xy_xyyyz = pbuffer.data(idx_dip_dh + 284);

    auto tr_z_xy_xyyzz = pbuffer.data(idx_dip_dh + 285);

    auto tr_z_xy_xyzzz = pbuffer.data(idx_dip_dh + 286);

    auto tr_z_xy_xzzzz = pbuffer.data(idx_dip_dh + 287);

    auto tr_z_xy_yyyyy = pbuffer.data(idx_dip_dh + 288);

    auto tr_z_xy_yyyyz = pbuffer.data(idx_dip_dh + 289);

    auto tr_z_xy_yyyzz = pbuffer.data(idx_dip_dh + 290);

    auto tr_z_xy_yyzzz = pbuffer.data(idx_dip_dh + 291);

    auto tr_z_xy_yzzzz = pbuffer.data(idx_dip_dh + 292);

    auto tr_z_xy_zzzzz = pbuffer.data(idx_dip_dh + 293);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_x_xxxxx,  \
                             tr_z_x_xxxxz,  \
                             tr_z_x_xxxzz,  \
                             tr_z_x_xxzzz,  \
                             tr_z_x_xzzzz,  \
                             tr_z_xy_xxxxx, \
                             tr_z_xy_xxxxy, \
                             tr_z_xy_xxxxz, \
                             tr_z_xy_xxxyy, \
                             tr_z_xy_xxxyz, \
                             tr_z_xy_xxxzz, \
                             tr_z_xy_xxyyy, \
                             tr_z_xy_xxyyz, \
                             tr_z_xy_xxyzz, \
                             tr_z_xy_xxzzz, \
                             tr_z_xy_xyyyy, \
                             tr_z_xy_xyyyz, \
                             tr_z_xy_xyyzz, \
                             tr_z_xy_xyzzz, \
                             tr_z_xy_xzzzz, \
                             tr_z_xy_yyyyy, \
                             tr_z_xy_yyyyz, \
                             tr_z_xy_yyyzz, \
                             tr_z_xy_yyzzz, \
                             tr_z_xy_yzzzz, \
                             tr_z_xy_zzzzz, \
                             tr_z_y_xxxxy,  \
                             tr_z_y_xxxy,   \
                             tr_z_y_xxxyy,  \
                             tr_z_y_xxxyz,  \
                             tr_z_y_xxyy,   \
                             tr_z_y_xxyyy,  \
                             tr_z_y_xxyyz,  \
                             tr_z_y_xxyz,   \
                             tr_z_y_xxyzz,  \
                             tr_z_y_xyyy,   \
                             tr_z_y_xyyyy,  \
                             tr_z_y_xyyyz,  \
                             tr_z_y_xyyz,   \
                             tr_z_y_xyyzz,  \
                             tr_z_y_xyzz,   \
                             tr_z_y_xyzzz,  \
                             tr_z_y_yyyy,   \
                             tr_z_y_yyyyy,  \
                             tr_z_y_yyyyz,  \
                             tr_z_y_yyyz,   \
                             tr_z_y_yyyzz,  \
                             tr_z_y_yyzz,   \
                             tr_z_y_yyzzz,  \
                             tr_z_y_yzzz,   \
                             tr_z_y_yzzzz,  \
                             tr_z_y_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xy_xxxxx[i] = tr_z_x_xxxxx[i] * pa_y[i];

        tr_z_xy_xxxxy[i] = 4.0 * tr_z_y_xxxy[i] * fe_0 + tr_z_y_xxxxy[i] * pa_x[i];

        tr_z_xy_xxxxz[i] = tr_z_x_xxxxz[i] * pa_y[i];

        tr_z_xy_xxxyy[i] = 3.0 * tr_z_y_xxyy[i] * fe_0 + tr_z_y_xxxyy[i] * pa_x[i];

        tr_z_xy_xxxyz[i] = 3.0 * tr_z_y_xxyz[i] * fe_0 + tr_z_y_xxxyz[i] * pa_x[i];

        tr_z_xy_xxxzz[i] = tr_z_x_xxxzz[i] * pa_y[i];

        tr_z_xy_xxyyy[i] = 2.0 * tr_z_y_xyyy[i] * fe_0 + tr_z_y_xxyyy[i] * pa_x[i];

        tr_z_xy_xxyyz[i] = 2.0 * tr_z_y_xyyz[i] * fe_0 + tr_z_y_xxyyz[i] * pa_x[i];

        tr_z_xy_xxyzz[i] = 2.0 * tr_z_y_xyzz[i] * fe_0 + tr_z_y_xxyzz[i] * pa_x[i];

        tr_z_xy_xxzzz[i] = tr_z_x_xxzzz[i] * pa_y[i];

        tr_z_xy_xyyyy[i] = tr_z_y_yyyy[i] * fe_0 + tr_z_y_xyyyy[i] * pa_x[i];

        tr_z_xy_xyyyz[i] = tr_z_y_yyyz[i] * fe_0 + tr_z_y_xyyyz[i] * pa_x[i];

        tr_z_xy_xyyzz[i] = tr_z_y_yyzz[i] * fe_0 + tr_z_y_xyyzz[i] * pa_x[i];

        tr_z_xy_xyzzz[i] = tr_z_y_yzzz[i] * fe_0 + tr_z_y_xyzzz[i] * pa_x[i];

        tr_z_xy_xzzzz[i] = tr_z_x_xzzzz[i] * pa_y[i];

        tr_z_xy_yyyyy[i] = tr_z_y_yyyyy[i] * pa_x[i];

        tr_z_xy_yyyyz[i] = tr_z_y_yyyyz[i] * pa_x[i];

        tr_z_xy_yyyzz[i] = tr_z_y_yyyzz[i] * pa_x[i];

        tr_z_xy_yyzzz[i] = tr_z_y_yyzzz[i] * pa_x[i];

        tr_z_xy_yzzzz[i] = tr_z_y_yzzzz[i] * pa_x[i];

        tr_z_xy_zzzzz[i] = tr_z_y_zzzzz[i] * pa_x[i];
    }

    // Set up 294-315 components of targeted buffer : DH

    auto tr_z_xz_xxxxx = pbuffer.data(idx_dip_dh + 294);

    auto tr_z_xz_xxxxy = pbuffer.data(idx_dip_dh + 295);

    auto tr_z_xz_xxxxz = pbuffer.data(idx_dip_dh + 296);

    auto tr_z_xz_xxxyy = pbuffer.data(idx_dip_dh + 297);

    auto tr_z_xz_xxxyz = pbuffer.data(idx_dip_dh + 298);

    auto tr_z_xz_xxxzz = pbuffer.data(idx_dip_dh + 299);

    auto tr_z_xz_xxyyy = pbuffer.data(idx_dip_dh + 300);

    auto tr_z_xz_xxyyz = pbuffer.data(idx_dip_dh + 301);

    auto tr_z_xz_xxyzz = pbuffer.data(idx_dip_dh + 302);

    auto tr_z_xz_xxzzz = pbuffer.data(idx_dip_dh + 303);

    auto tr_z_xz_xyyyy = pbuffer.data(idx_dip_dh + 304);

    auto tr_z_xz_xyyyz = pbuffer.data(idx_dip_dh + 305);

    auto tr_z_xz_xyyzz = pbuffer.data(idx_dip_dh + 306);

    auto tr_z_xz_xyzzz = pbuffer.data(idx_dip_dh + 307);

    auto tr_z_xz_xzzzz = pbuffer.data(idx_dip_dh + 308);

    auto tr_z_xz_yyyyy = pbuffer.data(idx_dip_dh + 309);

    auto tr_z_xz_yyyyz = pbuffer.data(idx_dip_dh + 310);

    auto tr_z_xz_yyyzz = pbuffer.data(idx_dip_dh + 311);

    auto tr_z_xz_yyzzz = pbuffer.data(idx_dip_dh + 312);

    auto tr_z_xz_yzzzz = pbuffer.data(idx_dip_dh + 313);

    auto tr_z_xz_zzzzz = pbuffer.data(idx_dip_dh + 314);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xz_xxxxx, \
                             tr_z_xz_xxxxy, \
                             tr_z_xz_xxxxz, \
                             tr_z_xz_xxxyy, \
                             tr_z_xz_xxxyz, \
                             tr_z_xz_xxxzz, \
                             tr_z_xz_xxyyy, \
                             tr_z_xz_xxyyz, \
                             tr_z_xz_xxyzz, \
                             tr_z_xz_xxzzz, \
                             tr_z_xz_xyyyy, \
                             tr_z_xz_xyyyz, \
                             tr_z_xz_xyyzz, \
                             tr_z_xz_xyzzz, \
                             tr_z_xz_xzzzz, \
                             tr_z_xz_yyyyy, \
                             tr_z_xz_yyyyz, \
                             tr_z_xz_yyyzz, \
                             tr_z_xz_yyzzz, \
                             tr_z_xz_yzzzz, \
                             tr_z_xz_zzzzz, \
                             tr_z_z_xxxx,   \
                             tr_z_z_xxxxx,  \
                             tr_z_z_xxxxy,  \
                             tr_z_z_xxxxz,  \
                             tr_z_z_xxxy,   \
                             tr_z_z_xxxyy,  \
                             tr_z_z_xxxyz,  \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxxzz,  \
                             tr_z_z_xxyy,   \
                             tr_z_z_xxyyy,  \
                             tr_z_z_xxyyz,  \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxyzz,  \
                             tr_z_z_xxzz,   \
                             tr_z_z_xxzzz,  \
                             tr_z_z_xyyy,   \
                             tr_z_z_xyyyy,  \
                             tr_z_z_xyyyz,  \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyyzz,  \
                             tr_z_z_xyzz,   \
                             tr_z_z_xyzzz,  \
                             tr_z_z_xzzz,   \
                             tr_z_z_xzzzz,  \
                             tr_z_z_yyyy,   \
                             tr_z_z_yyyyy,  \
                             tr_z_z_yyyyz,  \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyyzz,  \
                             tr_z_z_yyzz,   \
                             tr_z_z_yyzzz,  \
                             tr_z_z_yzzz,   \
                             tr_z_z_yzzzz,  \
                             tr_z_z_zzzz,   \
                             tr_z_z_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_xxxxx[i] = 5.0 * tr_z_z_xxxx[i] * fe_0 + tr_z_z_xxxxx[i] * pa_x[i];

        tr_z_xz_xxxxy[i] = 4.0 * tr_z_z_xxxy[i] * fe_0 + tr_z_z_xxxxy[i] * pa_x[i];

        tr_z_xz_xxxxz[i] = 4.0 * tr_z_z_xxxz[i] * fe_0 + tr_z_z_xxxxz[i] * pa_x[i];

        tr_z_xz_xxxyy[i] = 3.0 * tr_z_z_xxyy[i] * fe_0 + tr_z_z_xxxyy[i] * pa_x[i];

        tr_z_xz_xxxyz[i] = 3.0 * tr_z_z_xxyz[i] * fe_0 + tr_z_z_xxxyz[i] * pa_x[i];

        tr_z_xz_xxxzz[i] = 3.0 * tr_z_z_xxzz[i] * fe_0 + tr_z_z_xxxzz[i] * pa_x[i];

        tr_z_xz_xxyyy[i] = 2.0 * tr_z_z_xyyy[i] * fe_0 + tr_z_z_xxyyy[i] * pa_x[i];

        tr_z_xz_xxyyz[i] = 2.0 * tr_z_z_xyyz[i] * fe_0 + tr_z_z_xxyyz[i] * pa_x[i];

        tr_z_xz_xxyzz[i] = 2.0 * tr_z_z_xyzz[i] * fe_0 + tr_z_z_xxyzz[i] * pa_x[i];

        tr_z_xz_xxzzz[i] = 2.0 * tr_z_z_xzzz[i] * fe_0 + tr_z_z_xxzzz[i] * pa_x[i];

        tr_z_xz_xyyyy[i] = tr_z_z_yyyy[i] * fe_0 + tr_z_z_xyyyy[i] * pa_x[i];

        tr_z_xz_xyyyz[i] = tr_z_z_yyyz[i] * fe_0 + tr_z_z_xyyyz[i] * pa_x[i];

        tr_z_xz_xyyzz[i] = tr_z_z_yyzz[i] * fe_0 + tr_z_z_xyyzz[i] * pa_x[i];

        tr_z_xz_xyzzz[i] = tr_z_z_yzzz[i] * fe_0 + tr_z_z_xyzzz[i] * pa_x[i];

        tr_z_xz_xzzzz[i] = tr_z_z_zzzz[i] * fe_0 + tr_z_z_xzzzz[i] * pa_x[i];

        tr_z_xz_yyyyy[i] = tr_z_z_yyyyy[i] * pa_x[i];

        tr_z_xz_yyyyz[i] = tr_z_z_yyyyz[i] * pa_x[i];

        tr_z_xz_yyyzz[i] = tr_z_z_yyyzz[i] * pa_x[i];

        tr_z_xz_yyzzz[i] = tr_z_z_yyzzz[i] * pa_x[i];

        tr_z_xz_yzzzz[i] = tr_z_z_yzzzz[i] * pa_x[i];

        tr_z_xz_zzzzz[i] = tr_z_z_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : DH

    auto tr_z_yy_xxxxx = pbuffer.data(idx_dip_dh + 315);

    auto tr_z_yy_xxxxy = pbuffer.data(idx_dip_dh + 316);

    auto tr_z_yy_xxxxz = pbuffer.data(idx_dip_dh + 317);

    auto tr_z_yy_xxxyy = pbuffer.data(idx_dip_dh + 318);

    auto tr_z_yy_xxxyz = pbuffer.data(idx_dip_dh + 319);

    auto tr_z_yy_xxxzz = pbuffer.data(idx_dip_dh + 320);

    auto tr_z_yy_xxyyy = pbuffer.data(idx_dip_dh + 321);

    auto tr_z_yy_xxyyz = pbuffer.data(idx_dip_dh + 322);

    auto tr_z_yy_xxyzz = pbuffer.data(idx_dip_dh + 323);

    auto tr_z_yy_xxzzz = pbuffer.data(idx_dip_dh + 324);

    auto tr_z_yy_xyyyy = pbuffer.data(idx_dip_dh + 325);

    auto tr_z_yy_xyyyz = pbuffer.data(idx_dip_dh + 326);

    auto tr_z_yy_xyyzz = pbuffer.data(idx_dip_dh + 327);

    auto tr_z_yy_xyzzz = pbuffer.data(idx_dip_dh + 328);

    auto tr_z_yy_xzzzz = pbuffer.data(idx_dip_dh + 329);

    auto tr_z_yy_yyyyy = pbuffer.data(idx_dip_dh + 330);

    auto tr_z_yy_yyyyz = pbuffer.data(idx_dip_dh + 331);

    auto tr_z_yy_yyyzz = pbuffer.data(idx_dip_dh + 332);

    auto tr_z_yy_yyzzz = pbuffer.data(idx_dip_dh + 333);

    auto tr_z_yy_yzzzz = pbuffer.data(idx_dip_dh + 334);

    auto tr_z_yy_zzzzz = pbuffer.data(idx_dip_dh + 335);

#pragma omp simd aligned(pa_y,              \
                             tr_z_0_xxxxx,  \
                             tr_z_0_xxxxy,  \
                             tr_z_0_xxxxz,  \
                             tr_z_0_xxxyy,  \
                             tr_z_0_xxxyz,  \
                             tr_z_0_xxxzz,  \
                             tr_z_0_xxyyy,  \
                             tr_z_0_xxyyz,  \
                             tr_z_0_xxyzz,  \
                             tr_z_0_xxzzz,  \
                             tr_z_0_xyyyy,  \
                             tr_z_0_xyyyz,  \
                             tr_z_0_xyyzz,  \
                             tr_z_0_xyzzz,  \
                             tr_z_0_xzzzz,  \
                             tr_z_0_yyyyy,  \
                             tr_z_0_yyyyz,  \
                             tr_z_0_yyyzz,  \
                             tr_z_0_yyzzz,  \
                             tr_z_0_yzzzz,  \
                             tr_z_0_zzzzz,  \
                             tr_z_y_xxxx,   \
                             tr_z_y_xxxxx,  \
                             tr_z_y_xxxxy,  \
                             tr_z_y_xxxxz,  \
                             tr_z_y_xxxy,   \
                             tr_z_y_xxxyy,  \
                             tr_z_y_xxxyz,  \
                             tr_z_y_xxxz,   \
                             tr_z_y_xxxzz,  \
                             tr_z_y_xxyy,   \
                             tr_z_y_xxyyy,  \
                             tr_z_y_xxyyz,  \
                             tr_z_y_xxyz,   \
                             tr_z_y_xxyzz,  \
                             tr_z_y_xxzz,   \
                             tr_z_y_xxzzz,  \
                             tr_z_y_xyyy,   \
                             tr_z_y_xyyyy,  \
                             tr_z_y_xyyyz,  \
                             tr_z_y_xyyz,   \
                             tr_z_y_xyyzz,  \
                             tr_z_y_xyzz,   \
                             tr_z_y_xyzzz,  \
                             tr_z_y_xzzz,   \
                             tr_z_y_xzzzz,  \
                             tr_z_y_yyyy,   \
                             tr_z_y_yyyyy,  \
                             tr_z_y_yyyyz,  \
                             tr_z_y_yyyz,   \
                             tr_z_y_yyyzz,  \
                             tr_z_y_yyzz,   \
                             tr_z_y_yyzzz,  \
                             tr_z_y_yzzz,   \
                             tr_z_y_yzzzz,  \
                             tr_z_y_zzzz,   \
                             tr_z_y_zzzzz,  \
                             tr_z_yy_xxxxx, \
                             tr_z_yy_xxxxy, \
                             tr_z_yy_xxxxz, \
                             tr_z_yy_xxxyy, \
                             tr_z_yy_xxxyz, \
                             tr_z_yy_xxxzz, \
                             tr_z_yy_xxyyy, \
                             tr_z_yy_xxyyz, \
                             tr_z_yy_xxyzz, \
                             tr_z_yy_xxzzz, \
                             tr_z_yy_xyyyy, \
                             tr_z_yy_xyyyz, \
                             tr_z_yy_xyyzz, \
                             tr_z_yy_xyzzz, \
                             tr_z_yy_xzzzz, \
                             tr_z_yy_yyyyy, \
                             tr_z_yy_yyyyz, \
                             tr_z_yy_yyyzz, \
                             tr_z_yy_yyzzz, \
                             tr_z_yy_yzzzz, \
                             tr_z_yy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_xxxxx[i] = tr_z_0_xxxxx[i] * fe_0 + tr_z_y_xxxxx[i] * pa_y[i];

        tr_z_yy_xxxxy[i] = tr_z_0_xxxxy[i] * fe_0 + tr_z_y_xxxx[i] * fe_0 + tr_z_y_xxxxy[i] * pa_y[i];

        tr_z_yy_xxxxz[i] = tr_z_0_xxxxz[i] * fe_0 + tr_z_y_xxxxz[i] * pa_y[i];

        tr_z_yy_xxxyy[i] = tr_z_0_xxxyy[i] * fe_0 + 2.0 * tr_z_y_xxxy[i] * fe_0 + tr_z_y_xxxyy[i] * pa_y[i];

        tr_z_yy_xxxyz[i] = tr_z_0_xxxyz[i] * fe_0 + tr_z_y_xxxz[i] * fe_0 + tr_z_y_xxxyz[i] * pa_y[i];

        tr_z_yy_xxxzz[i] = tr_z_0_xxxzz[i] * fe_0 + tr_z_y_xxxzz[i] * pa_y[i];

        tr_z_yy_xxyyy[i] = tr_z_0_xxyyy[i] * fe_0 + 3.0 * tr_z_y_xxyy[i] * fe_0 + tr_z_y_xxyyy[i] * pa_y[i];

        tr_z_yy_xxyyz[i] = tr_z_0_xxyyz[i] * fe_0 + 2.0 * tr_z_y_xxyz[i] * fe_0 + tr_z_y_xxyyz[i] * pa_y[i];

        tr_z_yy_xxyzz[i] = tr_z_0_xxyzz[i] * fe_0 + tr_z_y_xxzz[i] * fe_0 + tr_z_y_xxyzz[i] * pa_y[i];

        tr_z_yy_xxzzz[i] = tr_z_0_xxzzz[i] * fe_0 + tr_z_y_xxzzz[i] * pa_y[i];

        tr_z_yy_xyyyy[i] = tr_z_0_xyyyy[i] * fe_0 + 4.0 * tr_z_y_xyyy[i] * fe_0 + tr_z_y_xyyyy[i] * pa_y[i];

        tr_z_yy_xyyyz[i] = tr_z_0_xyyyz[i] * fe_0 + 3.0 * tr_z_y_xyyz[i] * fe_0 + tr_z_y_xyyyz[i] * pa_y[i];

        tr_z_yy_xyyzz[i] = tr_z_0_xyyzz[i] * fe_0 + 2.0 * tr_z_y_xyzz[i] * fe_0 + tr_z_y_xyyzz[i] * pa_y[i];

        tr_z_yy_xyzzz[i] = tr_z_0_xyzzz[i] * fe_0 + tr_z_y_xzzz[i] * fe_0 + tr_z_y_xyzzz[i] * pa_y[i];

        tr_z_yy_xzzzz[i] = tr_z_0_xzzzz[i] * fe_0 + tr_z_y_xzzzz[i] * pa_y[i];

        tr_z_yy_yyyyy[i] = tr_z_0_yyyyy[i] * fe_0 + 5.0 * tr_z_y_yyyy[i] * fe_0 + tr_z_y_yyyyy[i] * pa_y[i];

        tr_z_yy_yyyyz[i] = tr_z_0_yyyyz[i] * fe_0 + 4.0 * tr_z_y_yyyz[i] * fe_0 + tr_z_y_yyyyz[i] * pa_y[i];

        tr_z_yy_yyyzz[i] = tr_z_0_yyyzz[i] * fe_0 + 3.0 * tr_z_y_yyzz[i] * fe_0 + tr_z_y_yyyzz[i] * pa_y[i];

        tr_z_yy_yyzzz[i] = tr_z_0_yyzzz[i] * fe_0 + 2.0 * tr_z_y_yzzz[i] * fe_0 + tr_z_y_yyzzz[i] * pa_y[i];

        tr_z_yy_yzzzz[i] = tr_z_0_yzzzz[i] * fe_0 + tr_z_y_zzzz[i] * fe_0 + tr_z_y_yzzzz[i] * pa_y[i];

        tr_z_yy_zzzzz[i] = tr_z_0_zzzzz[i] * fe_0 + tr_z_y_zzzzz[i] * pa_y[i];
    }

    // Set up 336-357 components of targeted buffer : DH

    auto tr_z_yz_xxxxx = pbuffer.data(idx_dip_dh + 336);

    auto tr_z_yz_xxxxy = pbuffer.data(idx_dip_dh + 337);

    auto tr_z_yz_xxxxz = pbuffer.data(idx_dip_dh + 338);

    auto tr_z_yz_xxxyy = pbuffer.data(idx_dip_dh + 339);

    auto tr_z_yz_xxxyz = pbuffer.data(idx_dip_dh + 340);

    auto tr_z_yz_xxxzz = pbuffer.data(idx_dip_dh + 341);

    auto tr_z_yz_xxyyy = pbuffer.data(idx_dip_dh + 342);

    auto tr_z_yz_xxyyz = pbuffer.data(idx_dip_dh + 343);

    auto tr_z_yz_xxyzz = pbuffer.data(idx_dip_dh + 344);

    auto tr_z_yz_xxzzz = pbuffer.data(idx_dip_dh + 345);

    auto tr_z_yz_xyyyy = pbuffer.data(idx_dip_dh + 346);

    auto tr_z_yz_xyyyz = pbuffer.data(idx_dip_dh + 347);

    auto tr_z_yz_xyyzz = pbuffer.data(idx_dip_dh + 348);

    auto tr_z_yz_xyzzz = pbuffer.data(idx_dip_dh + 349);

    auto tr_z_yz_xzzzz = pbuffer.data(idx_dip_dh + 350);

    auto tr_z_yz_yyyyy = pbuffer.data(idx_dip_dh + 351);

    auto tr_z_yz_yyyyz = pbuffer.data(idx_dip_dh + 352);

    auto tr_z_yz_yyyzz = pbuffer.data(idx_dip_dh + 353);

    auto tr_z_yz_yyzzz = pbuffer.data(idx_dip_dh + 354);

    auto tr_z_yz_yzzzz = pbuffer.data(idx_dip_dh + 355);

    auto tr_z_yz_zzzzz = pbuffer.data(idx_dip_dh + 356);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yz_xxxxx, \
                             tr_z_yz_xxxxy, \
                             tr_z_yz_xxxxz, \
                             tr_z_yz_xxxyy, \
                             tr_z_yz_xxxyz, \
                             tr_z_yz_xxxzz, \
                             tr_z_yz_xxyyy, \
                             tr_z_yz_xxyyz, \
                             tr_z_yz_xxyzz, \
                             tr_z_yz_xxzzz, \
                             tr_z_yz_xyyyy, \
                             tr_z_yz_xyyyz, \
                             tr_z_yz_xyyzz, \
                             tr_z_yz_xyzzz, \
                             tr_z_yz_xzzzz, \
                             tr_z_yz_yyyyy, \
                             tr_z_yz_yyyyz, \
                             tr_z_yz_yyyzz, \
                             tr_z_yz_yyzzz, \
                             tr_z_yz_yzzzz, \
                             tr_z_yz_zzzzz, \
                             tr_z_z_xxxx,   \
                             tr_z_z_xxxxx,  \
                             tr_z_z_xxxxy,  \
                             tr_z_z_xxxxz,  \
                             tr_z_z_xxxy,   \
                             tr_z_z_xxxyy,  \
                             tr_z_z_xxxyz,  \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxxzz,  \
                             tr_z_z_xxyy,   \
                             tr_z_z_xxyyy,  \
                             tr_z_z_xxyyz,  \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxyzz,  \
                             tr_z_z_xxzz,   \
                             tr_z_z_xxzzz,  \
                             tr_z_z_xyyy,   \
                             tr_z_z_xyyyy,  \
                             tr_z_z_xyyyz,  \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyyzz,  \
                             tr_z_z_xyzz,   \
                             tr_z_z_xyzzz,  \
                             tr_z_z_xzzz,   \
                             tr_z_z_xzzzz,  \
                             tr_z_z_yyyy,   \
                             tr_z_z_yyyyy,  \
                             tr_z_z_yyyyz,  \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyyzz,  \
                             tr_z_z_yyzz,   \
                             tr_z_z_yyzzz,  \
                             tr_z_z_yzzz,   \
                             tr_z_z_yzzzz,  \
                             tr_z_z_zzzz,   \
                             tr_z_z_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_xxxxx[i] = tr_z_z_xxxxx[i] * pa_y[i];

        tr_z_yz_xxxxy[i] = tr_z_z_xxxx[i] * fe_0 + tr_z_z_xxxxy[i] * pa_y[i];

        tr_z_yz_xxxxz[i] = tr_z_z_xxxxz[i] * pa_y[i];

        tr_z_yz_xxxyy[i] = 2.0 * tr_z_z_xxxy[i] * fe_0 + tr_z_z_xxxyy[i] * pa_y[i];

        tr_z_yz_xxxyz[i] = tr_z_z_xxxz[i] * fe_0 + tr_z_z_xxxyz[i] * pa_y[i];

        tr_z_yz_xxxzz[i] = tr_z_z_xxxzz[i] * pa_y[i];

        tr_z_yz_xxyyy[i] = 3.0 * tr_z_z_xxyy[i] * fe_0 + tr_z_z_xxyyy[i] * pa_y[i];

        tr_z_yz_xxyyz[i] = 2.0 * tr_z_z_xxyz[i] * fe_0 + tr_z_z_xxyyz[i] * pa_y[i];

        tr_z_yz_xxyzz[i] = tr_z_z_xxzz[i] * fe_0 + tr_z_z_xxyzz[i] * pa_y[i];

        tr_z_yz_xxzzz[i] = tr_z_z_xxzzz[i] * pa_y[i];

        tr_z_yz_xyyyy[i] = 4.0 * tr_z_z_xyyy[i] * fe_0 + tr_z_z_xyyyy[i] * pa_y[i];

        tr_z_yz_xyyyz[i] = 3.0 * tr_z_z_xyyz[i] * fe_0 + tr_z_z_xyyyz[i] * pa_y[i];

        tr_z_yz_xyyzz[i] = 2.0 * tr_z_z_xyzz[i] * fe_0 + tr_z_z_xyyzz[i] * pa_y[i];

        tr_z_yz_xyzzz[i] = tr_z_z_xzzz[i] * fe_0 + tr_z_z_xyzzz[i] * pa_y[i];

        tr_z_yz_xzzzz[i] = tr_z_z_xzzzz[i] * pa_y[i];

        tr_z_yz_yyyyy[i] = 5.0 * tr_z_z_yyyy[i] * fe_0 + tr_z_z_yyyyy[i] * pa_y[i];

        tr_z_yz_yyyyz[i] = 4.0 * tr_z_z_yyyz[i] * fe_0 + tr_z_z_yyyyz[i] * pa_y[i];

        tr_z_yz_yyyzz[i] = 3.0 * tr_z_z_yyzz[i] * fe_0 + tr_z_z_yyyzz[i] * pa_y[i];

        tr_z_yz_yyzzz[i] = 2.0 * tr_z_z_yzzz[i] * fe_0 + tr_z_z_yyzzz[i] * pa_y[i];

        tr_z_yz_yzzzz[i] = tr_z_z_zzzz[i] * fe_0 + tr_z_z_yzzzz[i] * pa_y[i];

        tr_z_yz_zzzzz[i] = tr_z_z_zzzzz[i] * pa_y[i];
    }

    // Set up 357-378 components of targeted buffer : DH

    auto tr_z_zz_xxxxx = pbuffer.data(idx_dip_dh + 357);

    auto tr_z_zz_xxxxy = pbuffer.data(idx_dip_dh + 358);

    auto tr_z_zz_xxxxz = pbuffer.data(idx_dip_dh + 359);

    auto tr_z_zz_xxxyy = pbuffer.data(idx_dip_dh + 360);

    auto tr_z_zz_xxxyz = pbuffer.data(idx_dip_dh + 361);

    auto tr_z_zz_xxxzz = pbuffer.data(idx_dip_dh + 362);

    auto tr_z_zz_xxyyy = pbuffer.data(idx_dip_dh + 363);

    auto tr_z_zz_xxyyz = pbuffer.data(idx_dip_dh + 364);

    auto tr_z_zz_xxyzz = pbuffer.data(idx_dip_dh + 365);

    auto tr_z_zz_xxzzz = pbuffer.data(idx_dip_dh + 366);

    auto tr_z_zz_xyyyy = pbuffer.data(idx_dip_dh + 367);

    auto tr_z_zz_xyyyz = pbuffer.data(idx_dip_dh + 368);

    auto tr_z_zz_xyyzz = pbuffer.data(idx_dip_dh + 369);

    auto tr_z_zz_xyzzz = pbuffer.data(idx_dip_dh + 370);

    auto tr_z_zz_xzzzz = pbuffer.data(idx_dip_dh + 371);

    auto tr_z_zz_yyyyy = pbuffer.data(idx_dip_dh + 372);

    auto tr_z_zz_yyyyz = pbuffer.data(idx_dip_dh + 373);

    auto tr_z_zz_yyyzz = pbuffer.data(idx_dip_dh + 374);

    auto tr_z_zz_yyzzz = pbuffer.data(idx_dip_dh + 375);

    auto tr_z_zz_yzzzz = pbuffer.data(idx_dip_dh + 376);

    auto tr_z_zz_zzzzz = pbuffer.data(idx_dip_dh + 377);

#pragma omp simd aligned(pa_z,              \
                             tr_z_0_xxxxx,  \
                             tr_z_0_xxxxy,  \
                             tr_z_0_xxxxz,  \
                             tr_z_0_xxxyy,  \
                             tr_z_0_xxxyz,  \
                             tr_z_0_xxxzz,  \
                             tr_z_0_xxyyy,  \
                             tr_z_0_xxyyz,  \
                             tr_z_0_xxyzz,  \
                             tr_z_0_xxzzz,  \
                             tr_z_0_xyyyy,  \
                             tr_z_0_xyyyz,  \
                             tr_z_0_xyyzz,  \
                             tr_z_0_xyzzz,  \
                             tr_z_0_xzzzz,  \
                             tr_z_0_yyyyy,  \
                             tr_z_0_yyyyz,  \
                             tr_z_0_yyyzz,  \
                             tr_z_0_yyzzz,  \
                             tr_z_0_yzzzz,  \
                             tr_z_0_zzzzz,  \
                             tr_z_z_xxxx,   \
                             tr_z_z_xxxxx,  \
                             tr_z_z_xxxxy,  \
                             tr_z_z_xxxxz,  \
                             tr_z_z_xxxy,   \
                             tr_z_z_xxxyy,  \
                             tr_z_z_xxxyz,  \
                             tr_z_z_xxxz,   \
                             tr_z_z_xxxzz,  \
                             tr_z_z_xxyy,   \
                             tr_z_z_xxyyy,  \
                             tr_z_z_xxyyz,  \
                             tr_z_z_xxyz,   \
                             tr_z_z_xxyzz,  \
                             tr_z_z_xxzz,   \
                             tr_z_z_xxzzz,  \
                             tr_z_z_xyyy,   \
                             tr_z_z_xyyyy,  \
                             tr_z_z_xyyyz,  \
                             tr_z_z_xyyz,   \
                             tr_z_z_xyyzz,  \
                             tr_z_z_xyzz,   \
                             tr_z_z_xyzzz,  \
                             tr_z_z_xzzz,   \
                             tr_z_z_xzzzz,  \
                             tr_z_z_yyyy,   \
                             tr_z_z_yyyyy,  \
                             tr_z_z_yyyyz,  \
                             tr_z_z_yyyz,   \
                             tr_z_z_yyyzz,  \
                             tr_z_z_yyzz,   \
                             tr_z_z_yyzzz,  \
                             tr_z_z_yzzz,   \
                             tr_z_z_yzzzz,  \
                             tr_z_z_zzzz,   \
                             tr_z_z_zzzzz,  \
                             tr_z_zz_xxxxx, \
                             tr_z_zz_xxxxy, \
                             tr_z_zz_xxxxz, \
                             tr_z_zz_xxxyy, \
                             tr_z_zz_xxxyz, \
                             tr_z_zz_xxxzz, \
                             tr_z_zz_xxyyy, \
                             tr_z_zz_xxyyz, \
                             tr_z_zz_xxyzz, \
                             tr_z_zz_xxzzz, \
                             tr_z_zz_xyyyy, \
                             tr_z_zz_xyyyz, \
                             tr_z_zz_xyyzz, \
                             tr_z_zz_xyzzz, \
                             tr_z_zz_xzzzz, \
                             tr_z_zz_yyyyy, \
                             tr_z_zz_yyyyz, \
                             tr_z_zz_yyyzz, \
                             tr_z_zz_yyzzz, \
                             tr_z_zz_yzzzz, \
                             tr_z_zz_zzzzz, \
                             ts_z_xxxxx,    \
                             ts_z_xxxxy,    \
                             ts_z_xxxxz,    \
                             ts_z_xxxyy,    \
                             ts_z_xxxyz,    \
                             ts_z_xxxzz,    \
                             ts_z_xxyyy,    \
                             ts_z_xxyyz,    \
                             ts_z_xxyzz,    \
                             ts_z_xxzzz,    \
                             ts_z_xyyyy,    \
                             ts_z_xyyyz,    \
                             ts_z_xyyzz,    \
                             ts_z_xyzzz,    \
                             ts_z_xzzzz,    \
                             ts_z_yyyyy,    \
                             ts_z_yyyyz,    \
                             ts_z_yyyzz,    \
                             ts_z_yyzzz,    \
                             ts_z_yzzzz,    \
                             ts_z_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_xxxxx[i] = tr_z_0_xxxxx[i] * fe_0 + ts_z_xxxxx[i] * fe_0 + tr_z_z_xxxxx[i] * pa_z[i];

        tr_z_zz_xxxxy[i] = tr_z_0_xxxxy[i] * fe_0 + ts_z_xxxxy[i] * fe_0 + tr_z_z_xxxxy[i] * pa_z[i];

        tr_z_zz_xxxxz[i] = tr_z_0_xxxxz[i] * fe_0 + tr_z_z_xxxx[i] * fe_0 + ts_z_xxxxz[i] * fe_0 + tr_z_z_xxxxz[i] * pa_z[i];

        tr_z_zz_xxxyy[i] = tr_z_0_xxxyy[i] * fe_0 + ts_z_xxxyy[i] * fe_0 + tr_z_z_xxxyy[i] * pa_z[i];

        tr_z_zz_xxxyz[i] = tr_z_0_xxxyz[i] * fe_0 + tr_z_z_xxxy[i] * fe_0 + ts_z_xxxyz[i] * fe_0 + tr_z_z_xxxyz[i] * pa_z[i];

        tr_z_zz_xxxzz[i] = tr_z_0_xxxzz[i] * fe_0 + 2.0 * tr_z_z_xxxz[i] * fe_0 + ts_z_xxxzz[i] * fe_0 + tr_z_z_xxxzz[i] * pa_z[i];

        tr_z_zz_xxyyy[i] = tr_z_0_xxyyy[i] * fe_0 + ts_z_xxyyy[i] * fe_0 + tr_z_z_xxyyy[i] * pa_z[i];

        tr_z_zz_xxyyz[i] = tr_z_0_xxyyz[i] * fe_0 + tr_z_z_xxyy[i] * fe_0 + ts_z_xxyyz[i] * fe_0 + tr_z_z_xxyyz[i] * pa_z[i];

        tr_z_zz_xxyzz[i] = tr_z_0_xxyzz[i] * fe_0 + 2.0 * tr_z_z_xxyz[i] * fe_0 + ts_z_xxyzz[i] * fe_0 + tr_z_z_xxyzz[i] * pa_z[i];

        tr_z_zz_xxzzz[i] = tr_z_0_xxzzz[i] * fe_0 + 3.0 * tr_z_z_xxzz[i] * fe_0 + ts_z_xxzzz[i] * fe_0 + tr_z_z_xxzzz[i] * pa_z[i];

        tr_z_zz_xyyyy[i] = tr_z_0_xyyyy[i] * fe_0 + ts_z_xyyyy[i] * fe_0 + tr_z_z_xyyyy[i] * pa_z[i];

        tr_z_zz_xyyyz[i] = tr_z_0_xyyyz[i] * fe_0 + tr_z_z_xyyy[i] * fe_0 + ts_z_xyyyz[i] * fe_0 + tr_z_z_xyyyz[i] * pa_z[i];

        tr_z_zz_xyyzz[i] = tr_z_0_xyyzz[i] * fe_0 + 2.0 * tr_z_z_xyyz[i] * fe_0 + ts_z_xyyzz[i] * fe_0 + tr_z_z_xyyzz[i] * pa_z[i];

        tr_z_zz_xyzzz[i] = tr_z_0_xyzzz[i] * fe_0 + 3.0 * tr_z_z_xyzz[i] * fe_0 + ts_z_xyzzz[i] * fe_0 + tr_z_z_xyzzz[i] * pa_z[i];

        tr_z_zz_xzzzz[i] = tr_z_0_xzzzz[i] * fe_0 + 4.0 * tr_z_z_xzzz[i] * fe_0 + ts_z_xzzzz[i] * fe_0 + tr_z_z_xzzzz[i] * pa_z[i];

        tr_z_zz_yyyyy[i] = tr_z_0_yyyyy[i] * fe_0 + ts_z_yyyyy[i] * fe_0 + tr_z_z_yyyyy[i] * pa_z[i];

        tr_z_zz_yyyyz[i] = tr_z_0_yyyyz[i] * fe_0 + tr_z_z_yyyy[i] * fe_0 + ts_z_yyyyz[i] * fe_0 + tr_z_z_yyyyz[i] * pa_z[i];

        tr_z_zz_yyyzz[i] = tr_z_0_yyyzz[i] * fe_0 + 2.0 * tr_z_z_yyyz[i] * fe_0 + ts_z_yyyzz[i] * fe_0 + tr_z_z_yyyzz[i] * pa_z[i];

        tr_z_zz_yyzzz[i] = tr_z_0_yyzzz[i] * fe_0 + 3.0 * tr_z_z_yyzz[i] * fe_0 + ts_z_yyzzz[i] * fe_0 + tr_z_z_yyzzz[i] * pa_z[i];

        tr_z_zz_yzzzz[i] = tr_z_0_yzzzz[i] * fe_0 + 4.0 * tr_z_z_yzzz[i] * fe_0 + ts_z_yzzzz[i] * fe_0 + tr_z_z_yzzzz[i] * pa_z[i];

        tr_z_zz_zzzzz[i] = tr_z_0_zzzzz[i] * fe_0 + 5.0 * tr_z_z_zzzz[i] * fe_0 + ts_z_zzzzz[i] * fe_0 + tr_z_z_zzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
