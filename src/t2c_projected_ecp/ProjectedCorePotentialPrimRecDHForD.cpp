#include "ProjectedCorePotentialPrimRecDHForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_dh_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_dh_d_0_0_0,
                                        const size_t idx_sh_d_0_0_0,
                                        const size_t idx_ph_d_0_0_0,
                                        const size_t idx_pg_p_0_0_1,
                                        const size_t idx_ph_p_0_0_1,
                                        const size_t idx_sh_d_1_0_0,
                                        const size_t idx_ph_d_1_0_0,
                                        const size_t idx_sh_s_1_0_1,
                                        const size_t idx_ph_s_1_0_1,
                                        const int p,
                                        const size_t idx_sh_d_0_0_1,
                                        const size_t idx_ph_d_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up B center coordinates

    auto rb_x = factors.data(idx_b);

    auto rb_y = factors.data(idx_b + 1);

    auto rb_z = factors.data(idx_b + 2);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0);

    auto tg_0_xxxxy_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 1);

    auto tg_0_xxxxz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 2);

    auto tg_0_xxxyy_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 3);

    auto tg_0_xxxyz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 4);

    auto tg_0_xxxzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 5);

    auto tg_0_xxyyy_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 6);

    auto tg_0_xxyyz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 7);

    auto tg_0_xxyzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 8);

    auto tg_0_xxzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 9);

    auto tg_0_xyyyy_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 10);

    auto tg_0_xyyyz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 11);

    auto tg_0_xyyzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 12);

    auto tg_0_xyzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 13);

    auto tg_0_xzzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 14);

    auto tg_0_yyyyy_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 15);

    auto tg_0_yyyyz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 16);

    auto tg_0_yyyzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 17);

    auto tg_0_yyzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 18);

    auto tg_0_yzzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 19);

    auto tg_0_zzzzz_d_0_0_0 = pbuffer.data(idx_sh_d_0_0_0 + 20);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0);

    auto tg_x_xxxxy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 1);

    auto tg_x_xxxxz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 2);

    auto tg_x_xxxyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 3);

    auto tg_x_xxxyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 4);

    auto tg_x_xxxzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 5);

    auto tg_x_xxyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 6);

    auto tg_x_xxyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 7);

    auto tg_x_xxyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 8);

    auto tg_x_xxzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 9);

    auto tg_x_xyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 10);

    auto tg_x_xyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 11);

    auto tg_x_xyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 12);

    auto tg_x_xyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 13);

    auto tg_x_xzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 14);

    auto tg_x_yyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 15);

    auto tg_x_yyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 16);

    auto tg_x_yyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 17);

    auto tg_x_yyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 18);

    auto tg_x_yzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 19);

    auto tg_x_zzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 20);

    auto tg_y_xxxxx_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 21);

    auto tg_y_xxxxy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 22);

    auto tg_y_xxxxz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 23);

    auto tg_y_xxxyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 24);

    auto tg_y_xxxyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 25);

    auto tg_y_xxxzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 26);

    auto tg_y_xxyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 27);

    auto tg_y_xxyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 28);

    auto tg_y_xxyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 29);

    auto tg_y_xxzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 30);

    auto tg_y_xyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 31);

    auto tg_y_xyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 32);

    auto tg_y_xyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 33);

    auto tg_y_xyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 34);

    auto tg_y_xzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 35);

    auto tg_y_yyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 36);

    auto tg_y_yyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 37);

    auto tg_y_yyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 38);

    auto tg_y_yyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 39);

    auto tg_y_yzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 40);

    auto tg_y_zzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 41);

    auto tg_z_xxxxx_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 42);

    auto tg_z_xxxxy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 43);

    auto tg_z_xxxxz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 44);

    auto tg_z_xxxyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 45);

    auto tg_z_xxxyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 46);

    auto tg_z_xxxzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 47);

    auto tg_z_xxyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 48);

    auto tg_z_xxyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 49);

    auto tg_z_xxyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 50);

    auto tg_z_xxzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 51);

    auto tg_z_xyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 52);

    auto tg_z_xyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 53);

    auto tg_z_xyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 54);

    auto tg_z_xyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 55);

    auto tg_z_xzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 56);

    auto tg_z_yyyyy_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 57);

    auto tg_z_yyyyz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 58);

    auto tg_z_yyyzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 59);

    auto tg_z_yyzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 60);

    auto tg_z_yzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 61);

    auto tg_z_zzzzz_d_0_0_0 = pbuffer.data(idx_ph_d_0_0_0 + 62);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1);

    auto tg_x_xxxy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 1);

    auto tg_x_xxxz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 2);

    auto tg_x_xxyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 3);

    auto tg_x_xxyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 4);

    auto tg_x_xxzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 5);

    auto tg_x_xyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 6);

    auto tg_x_xyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 7);

    auto tg_x_xyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 8);

    auto tg_x_xzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 9);

    auto tg_x_yyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 10);

    auto tg_x_yyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 11);

    auto tg_x_yyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 12);

    auto tg_x_yzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 13);

    auto tg_x_zzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 14);

    auto tg_y_xxxx_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 15);

    auto tg_y_xxxy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 16);

    auto tg_y_xxxz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 17);

    auto tg_y_xxyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 18);

    auto tg_y_xxyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 19);

    auto tg_y_xxzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 20);

    auto tg_y_xyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 21);

    auto tg_y_xyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 22);

    auto tg_y_xyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 23);

    auto tg_y_xzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 24);

    auto tg_y_yyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 25);

    auto tg_y_yyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 26);

    auto tg_y_yyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 27);

    auto tg_y_yzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 28);

    auto tg_y_zzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 29);

    auto tg_z_xxxx_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 30);

    auto tg_z_xxxy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 31);

    auto tg_z_xxxz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 32);

    auto tg_z_xxyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 33);

    auto tg_z_xxyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 34);

    auto tg_z_xxzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 35);

    auto tg_z_xyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 36);

    auto tg_z_xyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 37);

    auto tg_z_xyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 38);

    auto tg_z_xzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 39);

    auto tg_z_yyyy_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 40);

    auto tg_z_yyyz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 41);

    auto tg_z_yyzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 42);

    auto tg_z_yzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 43);

    auto tg_z_zzzz_p_0_0_1 = pbuffer.data(idx_pg_p_0_0_1 + 44);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1);

    auto tg_x_xxxxy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 1);

    auto tg_x_xxxxz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 2);

    auto tg_x_xxxyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 3);

    auto tg_x_xxxyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 4);

    auto tg_x_xxxzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 5);

    auto tg_x_xxyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 6);

    auto tg_x_xxyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 7);

    auto tg_x_xxyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 8);

    auto tg_x_xxzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 9);

    auto tg_x_xyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 10);

    auto tg_x_xyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 11);

    auto tg_x_xyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 12);

    auto tg_x_xyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 13);

    auto tg_x_xzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 14);

    auto tg_x_yyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 15);

    auto tg_x_yyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 16);

    auto tg_x_yyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 17);

    auto tg_x_yyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 18);

    auto tg_x_yzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 19);

    auto tg_x_zzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 20);

    auto tg_y_xxxxx_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 21);

    auto tg_y_xxxxy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 22);

    auto tg_y_xxxxz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 23);

    auto tg_y_xxxyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 24);

    auto tg_y_xxxyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 25);

    auto tg_y_xxxzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 26);

    auto tg_y_xxyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 27);

    auto tg_y_xxyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 28);

    auto tg_y_xxyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 29);

    auto tg_y_xxzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 30);

    auto tg_y_xyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 31);

    auto tg_y_xyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 32);

    auto tg_y_xyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 33);

    auto tg_y_xyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 34);

    auto tg_y_xzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 35);

    auto tg_y_yyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 36);

    auto tg_y_yyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 37);

    auto tg_y_yyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 38);

    auto tg_y_yyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 39);

    auto tg_y_yzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 40);

    auto tg_y_zzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 41);

    auto tg_z_xxxxx_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 42);

    auto tg_z_xxxxy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 43);

    auto tg_z_xxxxz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 44);

    auto tg_z_xxxyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 45);

    auto tg_z_xxxyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 46);

    auto tg_z_xxxzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 47);

    auto tg_z_xxyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 48);

    auto tg_z_xxyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 49);

    auto tg_z_xxyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 50);

    auto tg_z_xxzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 51);

    auto tg_z_xyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 52);

    auto tg_z_xyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 53);

    auto tg_z_xyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 54);

    auto tg_z_xyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 55);

    auto tg_z_xzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 56);

    auto tg_z_yyyyy_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 57);

    auto tg_z_yyyyz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 58);

    auto tg_z_yyyzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 59);

    auto tg_z_yyzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 60);

    auto tg_z_yzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 61);

    auto tg_z_zzzzz_p_0_0_1 = pbuffer.data(idx_ph_p_0_0_1 + 62);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0);

    auto tg_0_xxxxy_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 1);

    auto tg_0_xxxxz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 2);

    auto tg_0_xxxyy_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 3);

    auto tg_0_xxxyz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 4);

    auto tg_0_xxxzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 5);

    auto tg_0_xxyyy_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 6);

    auto tg_0_xxyyz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 7);

    auto tg_0_xxyzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 8);

    auto tg_0_xxzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 9);

    auto tg_0_xyyyy_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 10);

    auto tg_0_xyyyz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 11);

    auto tg_0_xyyzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 12);

    auto tg_0_xyzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 13);

    auto tg_0_xzzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 14);

    auto tg_0_yyyyy_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 15);

    auto tg_0_yyyyz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 16);

    auto tg_0_yyyzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 17);

    auto tg_0_yyzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 18);

    auto tg_0_yzzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 19);

    auto tg_0_zzzzz_d_1_0_0 = pbuffer.data(idx_sh_d_1_0_0 + 20);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0);

    auto tg_x_xxxxy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 1);

    auto tg_x_xxxxz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 2);

    auto tg_x_xxxyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 3);

    auto tg_x_xxxyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 4);

    auto tg_x_xxxzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 5);

    auto tg_x_xxyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 6);

    auto tg_x_xxyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 7);

    auto tg_x_xxyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 8);

    auto tg_x_xxzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 9);

    auto tg_x_xyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 10);

    auto tg_x_xyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 11);

    auto tg_x_xyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 12);

    auto tg_x_xyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 13);

    auto tg_x_xzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 14);

    auto tg_x_yyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 15);

    auto tg_x_yyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 16);

    auto tg_x_yyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 17);

    auto tg_x_yyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 18);

    auto tg_x_yzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 19);

    auto tg_x_zzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 20);

    auto tg_y_xxxxx_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 21);

    auto tg_y_xxxxy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 22);

    auto tg_y_xxxxz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 23);

    auto tg_y_xxxyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 24);

    auto tg_y_xxxyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 25);

    auto tg_y_xxxzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 26);

    auto tg_y_xxyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 27);

    auto tg_y_xxyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 28);

    auto tg_y_xxyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 29);

    auto tg_y_xxzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 30);

    auto tg_y_xyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 31);

    auto tg_y_xyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 32);

    auto tg_y_xyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 33);

    auto tg_y_xyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 34);

    auto tg_y_xzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 35);

    auto tg_y_yyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 36);

    auto tg_y_yyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 37);

    auto tg_y_yyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 38);

    auto tg_y_yyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 39);

    auto tg_y_yzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 40);

    auto tg_y_zzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 41);

    auto tg_z_xxxxx_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 42);

    auto tg_z_xxxxy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 43);

    auto tg_z_xxxxz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 44);

    auto tg_z_xxxyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 45);

    auto tg_z_xxxyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 46);

    auto tg_z_xxxzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 47);

    auto tg_z_xxyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 48);

    auto tg_z_xxyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 49);

    auto tg_z_xxyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 50);

    auto tg_z_xxzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 51);

    auto tg_z_xyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 52);

    auto tg_z_xyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 53);

    auto tg_z_xyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 54);

    auto tg_z_xyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 55);

    auto tg_z_xzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 56);

    auto tg_z_yyyyy_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 57);

    auto tg_z_yyyyz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 58);

    auto tg_z_yyyzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 59);

    auto tg_z_yyzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 60);

    auto tg_z_yzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 61);

    auto tg_z_zzzzz_d_1_0_0 = pbuffer.data(idx_ph_d_1_0_0 + 62);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1);

    auto tg_0_xxxxy_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 1);

    auto tg_0_xxxxz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 2);

    auto tg_0_xxxyy_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 3);

    auto tg_0_xxxyz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 4);

    auto tg_0_xxxzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 5);

    auto tg_0_xxyyy_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 6);

    auto tg_0_xxyyz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 7);

    auto tg_0_xxyzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 8);

    auto tg_0_xxzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 9);

    auto tg_0_xyyyy_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 10);

    auto tg_0_xyyyz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 11);

    auto tg_0_xyyzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 12);

    auto tg_0_xyzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 13);

    auto tg_0_xzzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 14);

    auto tg_0_yyyyy_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 15);

    auto tg_0_yyyyz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 16);

    auto tg_0_yyyzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 17);

    auto tg_0_yyzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 18);

    auto tg_0_yzzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 19);

    auto tg_0_zzzzz_s_1_0_1 = pbuffer.data(idx_sh_s_1_0_1 + 20);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1);

    auto tg_x_xxxxy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 1);

    auto tg_x_xxxxz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 2);

    auto tg_x_xxxyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 3);

    auto tg_x_xxxyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 4);

    auto tg_x_xxxzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 5);

    auto tg_x_xxyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 6);

    auto tg_x_xxyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 7);

    auto tg_x_xxyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 8);

    auto tg_x_xxzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 9);

    auto tg_x_xyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 10);

    auto tg_x_xyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 11);

    auto tg_x_xyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 12);

    auto tg_x_xyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 13);

    auto tg_x_xzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 14);

    auto tg_x_yyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 15);

    auto tg_x_yyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 16);

    auto tg_x_yyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 17);

    auto tg_x_yyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 18);

    auto tg_x_yzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 19);

    auto tg_x_zzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 20);

    auto tg_y_xxxxx_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 21);

    auto tg_y_xxxxy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 22);

    auto tg_y_xxxxz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 23);

    auto tg_y_xxxyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 24);

    auto tg_y_xxxyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 25);

    auto tg_y_xxxzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 26);

    auto tg_y_xxyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 27);

    auto tg_y_xxyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 28);

    auto tg_y_xxyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 29);

    auto tg_y_xxzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 30);

    auto tg_y_xyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 31);

    auto tg_y_xyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 32);

    auto tg_y_xyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 33);

    auto tg_y_xyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 34);

    auto tg_y_xzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 35);

    auto tg_y_yyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 36);

    auto tg_y_yyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 37);

    auto tg_y_yyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 38);

    auto tg_y_yyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 39);

    auto tg_y_yzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 40);

    auto tg_y_zzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 41);

    auto tg_z_xxxxx_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 42);

    auto tg_z_xxxxy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 43);

    auto tg_z_xxxxz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 44);

    auto tg_z_xxxyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 45);

    auto tg_z_xxxyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 46);

    auto tg_z_xxxzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 47);

    auto tg_z_xxyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 48);

    auto tg_z_xxyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 49);

    auto tg_z_xxyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 50);

    auto tg_z_xxzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 51);

    auto tg_z_xyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 52);

    auto tg_z_xyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 53);

    auto tg_z_xyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 54);

    auto tg_z_xyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 55);

    auto tg_z_xzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 56);

    auto tg_z_yyyyy_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 57);

    auto tg_z_yyyyz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 58);

    auto tg_z_yyyzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 59);

    auto tg_z_yyzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 60);

    auto tg_z_yzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 61);

    auto tg_z_zzzzz_s_1_0_1 = pbuffer.data(idx_ph_s_1_0_1 + 62);

    // Set up components of targeted buffer : DH

    auto tg_xx_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0);

    auto tg_xx_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 1);

    auto tg_xx_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 2);

    auto tg_xx_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 3);

    auto tg_xx_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 4);

    auto tg_xx_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 5);

    auto tg_xx_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 6);

    auto tg_xx_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 7);

    auto tg_xx_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 8);

    auto tg_xx_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 9);

    auto tg_xx_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 10);

    auto tg_xx_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 11);

    auto tg_xx_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 12);

    auto tg_xx_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 13);

    auto tg_xx_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 14);

    auto tg_xx_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 15);

    auto tg_xx_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 16);

    auto tg_xx_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 17);

    auto tg_xx_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 18);

    auto tg_xx_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 19);

    auto tg_xx_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 20);

    auto tg_xy_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 21);

    auto tg_xy_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 22);

    auto tg_xy_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 23);

    auto tg_xy_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 24);

    auto tg_xy_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 25);

    auto tg_xy_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 26);

    auto tg_xy_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 27);

    auto tg_xy_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 28);

    auto tg_xy_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 29);

    auto tg_xy_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 30);

    auto tg_xy_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 31);

    auto tg_xy_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 32);

    auto tg_xy_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 33);

    auto tg_xy_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 34);

    auto tg_xy_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 35);

    auto tg_xy_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 36);

    auto tg_xy_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 37);

    auto tg_xy_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 38);

    auto tg_xy_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 39);

    auto tg_xy_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 40);

    auto tg_xy_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 41);

    auto tg_xz_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 42);

    auto tg_xz_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 43);

    auto tg_xz_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 44);

    auto tg_xz_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 45);

    auto tg_xz_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 46);

    auto tg_xz_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 47);

    auto tg_xz_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 48);

    auto tg_xz_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 49);

    auto tg_xz_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 50);

    auto tg_xz_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 51);

    auto tg_xz_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 52);

    auto tg_xz_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 53);

    auto tg_xz_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 54);

    auto tg_xz_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 55);

    auto tg_xz_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 56);

    auto tg_xz_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 57);

    auto tg_xz_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 58);

    auto tg_xz_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 59);

    auto tg_xz_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 60);

    auto tg_xz_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 61);

    auto tg_xz_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 62);

    auto tg_yy_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 63);

    auto tg_yy_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 64);

    auto tg_yy_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 65);

    auto tg_yy_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 66);

    auto tg_yy_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 67);

    auto tg_yy_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 68);

    auto tg_yy_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 69);

    auto tg_yy_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 70);

    auto tg_yy_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 71);

    auto tg_yy_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 72);

    auto tg_yy_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 73);

    auto tg_yy_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 74);

    auto tg_yy_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 75);

    auto tg_yy_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 76);

    auto tg_yy_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 77);

    auto tg_yy_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 78);

    auto tg_yy_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 79);

    auto tg_yy_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 80);

    auto tg_yy_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 81);

    auto tg_yy_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 82);

    auto tg_yy_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 83);

    auto tg_yz_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 84);

    auto tg_yz_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 85);

    auto tg_yz_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 86);

    auto tg_yz_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 87);

    auto tg_yz_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 88);

    auto tg_yz_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 89);

    auto tg_yz_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 90);

    auto tg_yz_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 91);

    auto tg_yz_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 92);

    auto tg_yz_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 93);

    auto tg_yz_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 94);

    auto tg_yz_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 95);

    auto tg_yz_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 96);

    auto tg_yz_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 97);

    auto tg_yz_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 98);

    auto tg_yz_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 99);

    auto tg_yz_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 100);

    auto tg_yz_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 101);

    auto tg_yz_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 102);

    auto tg_yz_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 103);

    auto tg_yz_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 104);

    auto tg_zz_xxxxx_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 105);

    auto tg_zz_xxxxy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 106);

    auto tg_zz_xxxxz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 107);

    auto tg_zz_xxxyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 108);

    auto tg_zz_xxxyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 109);

    auto tg_zz_xxxzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 110);

    auto tg_zz_xxyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 111);

    auto tg_zz_xxyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 112);

    auto tg_zz_xxyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 113);

    auto tg_zz_xxzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 114);

    auto tg_zz_xyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 115);

    auto tg_zz_xyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 116);

    auto tg_zz_xyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 117);

    auto tg_zz_xyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 118);

    auto tg_zz_xzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 119);

    auto tg_zz_yyyyy_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 120);

    auto tg_zz_yyyyz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 121);

    auto tg_zz_yyyzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 122);

    auto tg_zz_yyzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 123);

    auto tg_zz_yzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 124);

    auto tg_zz_zzzzz_d_0_0_0 = pbuffer.data(idx_dh_d_0_0_0 + 125);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxxx_d_0_0_0, tg_0_xxxxx_d_1_0_0, tg_0_xxxxx_s_1_0_1, tg_0_xxxxy_d_0_0_0, tg_0_xxxxy_d_1_0_0, tg_0_xxxxy_s_1_0_1, tg_0_xxxxz_d_0_0_0, tg_0_xxxxz_d_1_0_0, tg_0_xxxxz_s_1_0_1, tg_0_xxxyy_d_0_0_0, tg_0_xxxyy_d_1_0_0, tg_0_xxxyy_s_1_0_1, tg_0_xxxyz_d_0_0_0, tg_0_xxxyz_d_1_0_0, tg_0_xxxyz_s_1_0_1, tg_0_xxxzz_d_0_0_0, tg_0_xxxzz_d_1_0_0, tg_0_xxxzz_s_1_0_1, tg_0_xxyyy_d_0_0_0, tg_0_xxyyy_d_1_0_0, tg_0_xxyyy_s_1_0_1, tg_0_xxyyz_d_0_0_0, tg_0_xxyyz_d_1_0_0, tg_0_xxyyz_s_1_0_1, tg_0_xxyzz_d_0_0_0, tg_0_xxyzz_d_1_0_0, tg_0_xxyzz_s_1_0_1, tg_0_xxzzz_d_0_0_0, tg_0_xxzzz_d_1_0_0, tg_0_xxzzz_s_1_0_1, tg_0_xyyyy_d_0_0_0, tg_0_xyyyy_d_1_0_0, tg_0_xyyyy_s_1_0_1, tg_0_xyyyz_d_0_0_0, tg_0_xyyyz_d_1_0_0, tg_0_xyyyz_s_1_0_1, tg_0_xyyzz_d_0_0_0, tg_0_xyyzz_d_1_0_0, tg_0_xyyzz_s_1_0_1, tg_0_xyzzz_d_0_0_0, tg_0_xyzzz_d_1_0_0, tg_0_xyzzz_s_1_0_1, tg_0_xzzzz_d_0_0_0, tg_0_xzzzz_d_1_0_0, tg_0_xzzzz_s_1_0_1, tg_0_yyyyy_d_0_0_0, tg_0_yyyyy_d_1_0_0, tg_0_yyyyy_s_1_0_1, tg_0_yyyyz_d_0_0_0, tg_0_yyyyz_d_1_0_0, tg_0_yyyyz_s_1_0_1, tg_0_yyyzz_d_0_0_0, tg_0_yyyzz_d_1_0_0, tg_0_yyyzz_s_1_0_1, tg_0_yyzzz_d_0_0_0, tg_0_yyzzz_d_1_0_0, tg_0_yyzzz_s_1_0_1, tg_0_yzzzz_d_0_0_0, tg_0_yzzzz_d_1_0_0, tg_0_yzzzz_s_1_0_1, tg_0_zzzzz_d_0_0_0, tg_0_zzzzz_d_1_0_0, tg_0_zzzzz_s_1_0_1, tg_x_xxxx_p_0_0_1, tg_x_xxxxx_d_0_0_0, tg_x_xxxxx_d_1_0_0, tg_x_xxxxx_p_0_0_1, tg_x_xxxxx_s_1_0_1, tg_x_xxxxy_d_0_0_0, tg_x_xxxxy_d_1_0_0, tg_x_xxxxy_p_0_0_1, tg_x_xxxxy_s_1_0_1, tg_x_xxxxz_d_0_0_0, tg_x_xxxxz_d_1_0_0, tg_x_xxxxz_p_0_0_1, tg_x_xxxxz_s_1_0_1, tg_x_xxxy_p_0_0_1, tg_x_xxxyy_d_0_0_0, tg_x_xxxyy_d_1_0_0, tg_x_xxxyy_p_0_0_1, tg_x_xxxyy_s_1_0_1, tg_x_xxxyz_d_0_0_0, tg_x_xxxyz_d_1_0_0, tg_x_xxxyz_p_0_0_1, tg_x_xxxyz_s_1_0_1, tg_x_xxxz_p_0_0_1, tg_x_xxxzz_d_0_0_0, tg_x_xxxzz_d_1_0_0, tg_x_xxxzz_p_0_0_1, tg_x_xxxzz_s_1_0_1, tg_x_xxyy_p_0_0_1, tg_x_xxyyy_d_0_0_0, tg_x_xxyyy_d_1_0_0, tg_x_xxyyy_p_0_0_1, tg_x_xxyyy_s_1_0_1, tg_x_xxyyz_d_0_0_0, tg_x_xxyyz_d_1_0_0, tg_x_xxyyz_p_0_0_1, tg_x_xxyyz_s_1_0_1, tg_x_xxyz_p_0_0_1, tg_x_xxyzz_d_0_0_0, tg_x_xxyzz_d_1_0_0, tg_x_xxyzz_p_0_0_1, tg_x_xxyzz_s_1_0_1, tg_x_xxzz_p_0_0_1, tg_x_xxzzz_d_0_0_0, tg_x_xxzzz_d_1_0_0, tg_x_xxzzz_p_0_0_1, tg_x_xxzzz_s_1_0_1, tg_x_xyyy_p_0_0_1, tg_x_xyyyy_d_0_0_0, tg_x_xyyyy_d_1_0_0, tg_x_xyyyy_p_0_0_1, tg_x_xyyyy_s_1_0_1, tg_x_xyyyz_d_0_0_0, tg_x_xyyyz_d_1_0_0, tg_x_xyyyz_p_0_0_1, tg_x_xyyyz_s_1_0_1, tg_x_xyyz_p_0_0_1, tg_x_xyyzz_d_0_0_0, tg_x_xyyzz_d_1_0_0, tg_x_xyyzz_p_0_0_1, tg_x_xyyzz_s_1_0_1, tg_x_xyzz_p_0_0_1, tg_x_xyzzz_d_0_0_0, tg_x_xyzzz_d_1_0_0, tg_x_xyzzz_p_0_0_1, tg_x_xyzzz_s_1_0_1, tg_x_xzzz_p_0_0_1, tg_x_xzzzz_d_0_0_0, tg_x_xzzzz_d_1_0_0, tg_x_xzzzz_p_0_0_1, tg_x_xzzzz_s_1_0_1, tg_x_yyyy_p_0_0_1, tg_x_yyyyy_d_0_0_0, tg_x_yyyyy_d_1_0_0, tg_x_yyyyy_p_0_0_1, tg_x_yyyyy_s_1_0_1, tg_x_yyyyz_d_0_0_0, tg_x_yyyyz_d_1_0_0, tg_x_yyyyz_p_0_0_1, tg_x_yyyyz_s_1_0_1, tg_x_yyyz_p_0_0_1, tg_x_yyyzz_d_0_0_0, tg_x_yyyzz_d_1_0_0, tg_x_yyyzz_p_0_0_1, tg_x_yyyzz_s_1_0_1, tg_x_yyzz_p_0_0_1, tg_x_yyzzz_d_0_0_0, tg_x_yyzzz_d_1_0_0, tg_x_yyzzz_p_0_0_1, tg_x_yyzzz_s_1_0_1, tg_x_yzzz_p_0_0_1, tg_x_yzzzz_d_0_0_0, tg_x_yzzzz_d_1_0_0, tg_x_yzzzz_p_0_0_1, tg_x_yzzzz_s_1_0_1, tg_x_zzzz_p_0_0_1, tg_x_zzzzz_d_0_0_0, tg_x_zzzzz_d_1_0_0, tg_x_zzzzz_p_0_0_1, tg_x_zzzzz_s_1_0_1, tg_xx_xxxxx_d_0_0_0, tg_xx_xxxxy_d_0_0_0, tg_xx_xxxxz_d_0_0_0, tg_xx_xxxyy_d_0_0_0, tg_xx_xxxyz_d_0_0_0, tg_xx_xxxzz_d_0_0_0, tg_xx_xxyyy_d_0_0_0, tg_xx_xxyyz_d_0_0_0, tg_xx_xxyzz_d_0_0_0, tg_xx_xxzzz_d_0_0_0, tg_xx_xyyyy_d_0_0_0, tg_xx_xyyyz_d_0_0_0, tg_xx_xyyzz_d_0_0_0, tg_xx_xyzzz_d_0_0_0, tg_xx_xzzzz_d_0_0_0, tg_xx_yyyyy_d_0_0_0, tg_xx_yyyyz_d_0_0_0, tg_xx_yyyzz_d_0_0_0, tg_xx_yyzzz_d_0_0_0, tg_xx_yzzzz_d_0_0_0, tg_xx_zzzzz_d_0_0_0, tg_xy_xxxxx_d_0_0_0, tg_xy_xxxxy_d_0_0_0, tg_xy_xxxxz_d_0_0_0, tg_xy_xxxyy_d_0_0_0, tg_xy_xxxyz_d_0_0_0, tg_xy_xxxzz_d_0_0_0, tg_xy_xxyyy_d_0_0_0, tg_xy_xxyyz_d_0_0_0, tg_xy_xxyzz_d_0_0_0, tg_xy_xxzzz_d_0_0_0, tg_xy_xyyyy_d_0_0_0, tg_xy_xyyyz_d_0_0_0, tg_xy_xyyzz_d_0_0_0, tg_xy_xyzzz_d_0_0_0, tg_xy_xzzzz_d_0_0_0, tg_xy_yyyyy_d_0_0_0, tg_xy_yyyyz_d_0_0_0, tg_xy_yyyzz_d_0_0_0, tg_xy_yyzzz_d_0_0_0, tg_xy_yzzzz_d_0_0_0, tg_xy_zzzzz_d_0_0_0, tg_xz_xxxxx_d_0_0_0, tg_xz_xxxxy_d_0_0_0, tg_xz_xxxxz_d_0_0_0, tg_xz_xxxyy_d_0_0_0, tg_xz_xxxyz_d_0_0_0, tg_xz_xxxzz_d_0_0_0, tg_xz_xxyyy_d_0_0_0, tg_xz_xxyyz_d_0_0_0, tg_xz_xxyzz_d_0_0_0, tg_xz_xxzzz_d_0_0_0, tg_xz_xyyyy_d_0_0_0, tg_xz_xyyyz_d_0_0_0, tg_xz_xyyzz_d_0_0_0, tg_xz_xyzzz_d_0_0_0, tg_xz_xzzzz_d_0_0_0, tg_xz_yyyyy_d_0_0_0, tg_xz_yyyyz_d_0_0_0, tg_xz_yyyzz_d_0_0_0, tg_xz_yyzzz_d_0_0_0, tg_xz_yzzzz_d_0_0_0, tg_xz_zzzzz_d_0_0_0, tg_y_xxxx_p_0_0_1, tg_y_xxxxx_d_0_0_0, tg_y_xxxxx_d_1_0_0, tg_y_xxxxx_p_0_0_1, tg_y_xxxxx_s_1_0_1, tg_y_xxxxy_d_0_0_0, tg_y_xxxxy_d_1_0_0, tg_y_xxxxy_p_0_0_1, tg_y_xxxxy_s_1_0_1, tg_y_xxxxz_d_0_0_0, tg_y_xxxxz_d_1_0_0, tg_y_xxxxz_p_0_0_1, tg_y_xxxxz_s_1_0_1, tg_y_xxxy_p_0_0_1, tg_y_xxxyy_d_0_0_0, tg_y_xxxyy_d_1_0_0, tg_y_xxxyy_p_0_0_1, tg_y_xxxyy_s_1_0_1, tg_y_xxxyz_d_0_0_0, tg_y_xxxyz_d_1_0_0, tg_y_xxxyz_p_0_0_1, tg_y_xxxyz_s_1_0_1, tg_y_xxxz_p_0_0_1, tg_y_xxxzz_d_0_0_0, tg_y_xxxzz_d_1_0_0, tg_y_xxxzz_p_0_0_1, tg_y_xxxzz_s_1_0_1, tg_y_xxyy_p_0_0_1, tg_y_xxyyy_d_0_0_0, tg_y_xxyyy_d_1_0_0, tg_y_xxyyy_p_0_0_1, tg_y_xxyyy_s_1_0_1, tg_y_xxyyz_d_0_0_0, tg_y_xxyyz_d_1_0_0, tg_y_xxyyz_p_0_0_1, tg_y_xxyyz_s_1_0_1, tg_y_xxyz_p_0_0_1, tg_y_xxyzz_d_0_0_0, tg_y_xxyzz_d_1_0_0, tg_y_xxyzz_p_0_0_1, tg_y_xxyzz_s_1_0_1, tg_y_xxzz_p_0_0_1, tg_y_xxzzz_d_0_0_0, tg_y_xxzzz_d_1_0_0, tg_y_xxzzz_p_0_0_1, tg_y_xxzzz_s_1_0_1, tg_y_xyyy_p_0_0_1, tg_y_xyyyy_d_0_0_0, tg_y_xyyyy_d_1_0_0, tg_y_xyyyy_p_0_0_1, tg_y_xyyyy_s_1_0_1, tg_y_xyyyz_d_0_0_0, tg_y_xyyyz_d_1_0_0, tg_y_xyyyz_p_0_0_1, tg_y_xyyyz_s_1_0_1, tg_y_xyyz_p_0_0_1, tg_y_xyyzz_d_0_0_0, tg_y_xyyzz_d_1_0_0, tg_y_xyyzz_p_0_0_1, tg_y_xyyzz_s_1_0_1, tg_y_xyzz_p_0_0_1, tg_y_xyzzz_d_0_0_0, tg_y_xyzzz_d_1_0_0, tg_y_xyzzz_p_0_0_1, tg_y_xyzzz_s_1_0_1, tg_y_xzzz_p_0_0_1, tg_y_xzzzz_d_0_0_0, tg_y_xzzzz_d_1_0_0, tg_y_xzzzz_p_0_0_1, tg_y_xzzzz_s_1_0_1, tg_y_yyyy_p_0_0_1, tg_y_yyyyy_d_0_0_0, tg_y_yyyyy_d_1_0_0, tg_y_yyyyy_p_0_0_1, tg_y_yyyyy_s_1_0_1, tg_y_yyyyz_d_0_0_0, tg_y_yyyyz_d_1_0_0, tg_y_yyyyz_p_0_0_1, tg_y_yyyyz_s_1_0_1, tg_y_yyyz_p_0_0_1, tg_y_yyyzz_d_0_0_0, tg_y_yyyzz_d_1_0_0, tg_y_yyyzz_p_0_0_1, tg_y_yyyzz_s_1_0_1, tg_y_yyzz_p_0_0_1, tg_y_yyzzz_d_0_0_0, tg_y_yyzzz_d_1_0_0, tg_y_yyzzz_p_0_0_1, tg_y_yyzzz_s_1_0_1, tg_y_yzzz_p_0_0_1, tg_y_yzzzz_d_0_0_0, tg_y_yzzzz_d_1_0_0, tg_y_yzzzz_p_0_0_1, tg_y_yzzzz_s_1_0_1, tg_y_zzzz_p_0_0_1, tg_y_zzzzz_d_0_0_0, tg_y_zzzzz_d_1_0_0, tg_y_zzzzz_p_0_0_1, tg_y_zzzzz_s_1_0_1, tg_yy_xxxxx_d_0_0_0, tg_yy_xxxxy_d_0_0_0, tg_yy_xxxxz_d_0_0_0, tg_yy_xxxyy_d_0_0_0, tg_yy_xxxyz_d_0_0_0, tg_yy_xxxzz_d_0_0_0, tg_yy_xxyyy_d_0_0_0, tg_yy_xxyyz_d_0_0_0, tg_yy_xxyzz_d_0_0_0, tg_yy_xxzzz_d_0_0_0, tg_yy_xyyyy_d_0_0_0, tg_yy_xyyyz_d_0_0_0, tg_yy_xyyzz_d_0_0_0, tg_yy_xyzzz_d_0_0_0, tg_yy_xzzzz_d_0_0_0, tg_yy_yyyyy_d_0_0_0, tg_yy_yyyyz_d_0_0_0, tg_yy_yyyzz_d_0_0_0, tg_yy_yyzzz_d_0_0_0, tg_yy_yzzzz_d_0_0_0, tg_yy_zzzzz_d_0_0_0, tg_yz_xxxxx_d_0_0_0, tg_yz_xxxxy_d_0_0_0, tg_yz_xxxxz_d_0_0_0, tg_yz_xxxyy_d_0_0_0, tg_yz_xxxyz_d_0_0_0, tg_yz_xxxzz_d_0_0_0, tg_yz_xxyyy_d_0_0_0, tg_yz_xxyyz_d_0_0_0, tg_yz_xxyzz_d_0_0_0, tg_yz_xxzzz_d_0_0_0, tg_yz_xyyyy_d_0_0_0, tg_yz_xyyyz_d_0_0_0, tg_yz_xyyzz_d_0_0_0, tg_yz_xyzzz_d_0_0_0, tg_yz_xzzzz_d_0_0_0, tg_yz_yyyyy_d_0_0_0, tg_yz_yyyyz_d_0_0_0, tg_yz_yyyzz_d_0_0_0, tg_yz_yyzzz_d_0_0_0, tg_yz_yzzzz_d_0_0_0, tg_yz_zzzzz_d_0_0_0, tg_z_xxxx_p_0_0_1, tg_z_xxxxx_d_0_0_0, tg_z_xxxxx_d_1_0_0, tg_z_xxxxx_p_0_0_1, tg_z_xxxxx_s_1_0_1, tg_z_xxxxy_d_0_0_0, tg_z_xxxxy_d_1_0_0, tg_z_xxxxy_p_0_0_1, tg_z_xxxxy_s_1_0_1, tg_z_xxxxz_d_0_0_0, tg_z_xxxxz_d_1_0_0, tg_z_xxxxz_p_0_0_1, tg_z_xxxxz_s_1_0_1, tg_z_xxxy_p_0_0_1, tg_z_xxxyy_d_0_0_0, tg_z_xxxyy_d_1_0_0, tg_z_xxxyy_p_0_0_1, tg_z_xxxyy_s_1_0_1, tg_z_xxxyz_d_0_0_0, tg_z_xxxyz_d_1_0_0, tg_z_xxxyz_p_0_0_1, tg_z_xxxyz_s_1_0_1, tg_z_xxxz_p_0_0_1, tg_z_xxxzz_d_0_0_0, tg_z_xxxzz_d_1_0_0, tg_z_xxxzz_p_0_0_1, tg_z_xxxzz_s_1_0_1, tg_z_xxyy_p_0_0_1, tg_z_xxyyy_d_0_0_0, tg_z_xxyyy_d_1_0_0, tg_z_xxyyy_p_0_0_1, tg_z_xxyyy_s_1_0_1, tg_z_xxyyz_d_0_0_0, tg_z_xxyyz_d_1_0_0, tg_z_xxyyz_p_0_0_1, tg_z_xxyyz_s_1_0_1, tg_z_xxyz_p_0_0_1, tg_z_xxyzz_d_0_0_0, tg_z_xxyzz_d_1_0_0, tg_z_xxyzz_p_0_0_1, tg_z_xxyzz_s_1_0_1, tg_z_xxzz_p_0_0_1, tg_z_xxzzz_d_0_0_0, tg_z_xxzzz_d_1_0_0, tg_z_xxzzz_p_0_0_1, tg_z_xxzzz_s_1_0_1, tg_z_xyyy_p_0_0_1, tg_z_xyyyy_d_0_0_0, tg_z_xyyyy_d_1_0_0, tg_z_xyyyy_p_0_0_1, tg_z_xyyyy_s_1_0_1, tg_z_xyyyz_d_0_0_0, tg_z_xyyyz_d_1_0_0, tg_z_xyyyz_p_0_0_1, tg_z_xyyyz_s_1_0_1, tg_z_xyyz_p_0_0_1, tg_z_xyyzz_d_0_0_0, tg_z_xyyzz_d_1_0_0, tg_z_xyyzz_p_0_0_1, tg_z_xyyzz_s_1_0_1, tg_z_xyzz_p_0_0_1, tg_z_xyzzz_d_0_0_0, tg_z_xyzzz_d_1_0_0, tg_z_xyzzz_p_0_0_1, tg_z_xyzzz_s_1_0_1, tg_z_xzzz_p_0_0_1, tg_z_xzzzz_d_0_0_0, tg_z_xzzzz_d_1_0_0, tg_z_xzzzz_p_0_0_1, tg_z_xzzzz_s_1_0_1, tg_z_yyyy_p_0_0_1, tg_z_yyyyy_d_0_0_0, tg_z_yyyyy_d_1_0_0, tg_z_yyyyy_p_0_0_1, tg_z_yyyyy_s_1_0_1, tg_z_yyyyz_d_0_0_0, tg_z_yyyyz_d_1_0_0, tg_z_yyyyz_p_0_0_1, tg_z_yyyyz_s_1_0_1, tg_z_yyyz_p_0_0_1, tg_z_yyyzz_d_0_0_0, tg_z_yyyzz_d_1_0_0, tg_z_yyyzz_p_0_0_1, tg_z_yyyzz_s_1_0_1, tg_z_yyzz_p_0_0_1, tg_z_yyzzz_d_0_0_0, tg_z_yyzzz_d_1_0_0, tg_z_yyzzz_p_0_0_1, tg_z_yyzzz_s_1_0_1, tg_z_yzzz_p_0_0_1, tg_z_yzzzz_d_0_0_0, tg_z_yzzzz_d_1_0_0, tg_z_yzzzz_p_0_0_1, tg_z_yzzzz_s_1_0_1, tg_z_zzzz_p_0_0_1, tg_z_zzzzz_d_0_0_0, tg_z_zzzzz_d_1_0_0, tg_z_zzzzz_p_0_0_1, tg_z_zzzzz_s_1_0_1, tg_zz_xxxxx_d_0_0_0, tg_zz_xxxxy_d_0_0_0, tg_zz_xxxxz_d_0_0_0, tg_zz_xxxyy_d_0_0_0, tg_zz_xxxyz_d_0_0_0, tg_zz_xxxzz_d_0_0_0, tg_zz_xxyyy_d_0_0_0, tg_zz_xxyyz_d_0_0_0, tg_zz_xxyzz_d_0_0_0, tg_zz_xxzzz_d_0_0_0, tg_zz_xyyyy_d_0_0_0, tg_zz_xyyyz_d_0_0_0, tg_zz_xyyzz_d_0_0_0, tg_zz_xyzzz_d_0_0_0, tg_zz_xzzzz_d_0_0_0, tg_zz_yyyyy_d_0_0_0, tg_zz_yyyyz_d_0_0_0, tg_zz_yyyzz_d_0_0_0, tg_zz_yyzzz_d_0_0_0, tg_zz_yzzzz_d_0_0_0, tg_zz_zzzzz_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xx_xxxxx_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxx_d_0_0_0[i] * fzi_0 + tg_0_xxxxx_d_1_0_0[i] * fbzi_0 * fbzi_0 + 25.0 / 2.0 * tg_x_xxxx_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxxx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxxx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxxx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxx_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxy_d_0_0_0[i] * fzi_0 + tg_0_xxxxy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_x_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxxy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxxy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxxy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxy_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxz_d_0_0_0[i] * fzi_0 + tg_0_xxxxz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_x_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxxz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxxz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxxz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyy_d_0_0_0[i] * fzi_0 + tg_0_xxxyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_x_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyy_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyz_d_0_0_0[i] * fzi_0 + tg_0_xxxyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_x_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxzz_d_0_0_0[i] * fzi_0 + tg_0_xxxzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_x_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxxzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxxzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyy_d_0_0_0[i] * fzi_0 + tg_0_xxyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_x_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyz_d_0_0_0[i] * fzi_0 + tg_0_xxyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_x_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyzz_d_0_0_0[i] * fzi_0 + tg_0_xxyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_x_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xxzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxzzz_d_0_0_0[i] * fzi_0 + tg_0_xxzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_x_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xxzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xxzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyy_d_0_0_0[i] * fzi_0 + tg_0_xyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_x_yyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xyyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xyyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyz_d_0_0_0[i] * fzi_0 + tg_0_xyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_x_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyzz_d_0_0_0[i] * fzi_0 + tg_0_xyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_x_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyzzz_d_0_0_0[i] * fzi_0 + tg_0_xyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_x_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_xzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xzzzz_d_0_0_0[i] * fzi_0 + tg_0_xzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_x_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_x_xzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_xzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyy_d_0_0_0[i] * fzi_0 + tg_0_yyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_yyyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_yyyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyz_d_0_0_0[i] * fzi_0 + tg_0_yyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_yyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_yyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyzz_d_0_0_0[i] * fzi_0 + tg_0_yyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_yyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_yyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_yyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyzzz_d_0_0_0[i] * fzi_0 + tg_0_yyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_yyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_yyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_yzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yzzzz_d_0_0_0[i] * fzi_0 + tg_0_yzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_yzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_yzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xx_zzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_zzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_zzzzz_d_0_0_0[i] * fzi_0 + tg_0_zzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_x_zzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_x_zzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_zzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_zzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxx_d_0_0_0[i] = -5.0 * tg_x_xxxxx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_x_xxxxx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxxxx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxx_d_0_0_0[i] * a_y * faz_0;

        tg_xy_xxxxy_d_0_0_0[i] = 10.0 * tg_y_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxxy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxxxy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxxxy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxy_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxz_d_0_0_0[i] = -5.0 * tg_x_xxxxz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_x_xxxxz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxxxz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxz_d_0_0_0[i] * a_y * faz_0;

        tg_xy_xxxyy_d_0_0_0[i] = 15.0 / 2.0 * tg_y_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxxyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxxyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyy_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxyz_d_0_0_0[i] = 15.0 / 2.0 * tg_y_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxxyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxxyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxzz_d_0_0_0[i] = -5.0 * tg_x_xxxzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_x_xxxzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxxzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxzz_d_0_0_0[i] * a_y * faz_0;

        tg_xy_xxyyy_d_0_0_0[i] = 5.0 * tg_y_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyyz_d_0_0_0[i] = 5.0 * tg_y_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyzz_d_0_0_0[i] = 5.0 * tg_y_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xxyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xxzzz_d_0_0_0[i] = -5.0 * tg_x_xxzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_x_xxzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxzzz_d_0_0_0[i] * a_y * faz_0;

        tg_xy_xyyyy_d_0_0_0[i] = 5.0 / 2.0 * tg_y_yyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xyyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyyz_d_0_0_0[i] = 5.0 / 2.0 * tg_y_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyzz_d_0_0_0[i] = 5.0 / 2.0 * tg_y_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xyzzz_d_0_0_0[i] = 5.0 / 2.0 * tg_y_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_xyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_xzzzz_d_0_0_0[i] = -5.0 * tg_x_xzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_x_xzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_xy_yyyyy_d_0_0_0[i] = -5.0 * tg_y_yyyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_yyyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyyz_d_0_0_0[i] = -5.0 * tg_y_yyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_yyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyzz_d_0_0_0[i] = -5.0 * tg_y_yyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_yyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_yyzzz_d_0_0_0[i] = -5.0 * tg_y_yyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_yyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_yzzzz_d_0_0_0[i] = -5.0 * tg_y_yzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_yzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xy_zzzzz_d_0_0_0[i] = -5.0 * tg_y_zzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_y_zzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_zzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxx_d_0_0_0[i] = -5.0 * tg_x_xxxxx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_x_xxxxx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxxxx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxx_d_0_0_0[i] * a_z * faz_0;

        tg_xz_xxxxy_d_0_0_0[i] = -5.0 * tg_x_xxxxy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_x_xxxxy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxxxy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxy_d_0_0_0[i] * a_z * faz_0;

        tg_xz_xxxxz_d_0_0_0[i] = 10.0 * tg_z_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxxz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxxxz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxxxz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxyy_d_0_0_0[i] = -5.0 * tg_x_xxxyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_x_xxxyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxxyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyy_d_0_0_0[i] * a_z * faz_0;

        tg_xz_xxxyz_d_0_0_0[i] = 15.0 / 2.0 * tg_z_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxxyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxxyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxzz_d_0_0_0[i] = 15.0 / 2.0 * tg_z_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxxzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxxzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyyy_d_0_0_0[i] = -5.0 * tg_x_xxyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_x_xxyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyy_d_0_0_0[i] * a_z * faz_0;

        tg_xz_xxyyz_d_0_0_0[i] = 5.0 * tg_z_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyzz_d_0_0_0[i] = 5.0 * tg_z_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xxzzz_d_0_0_0[i] = 5.0 * tg_z_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xxzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyyy_d_0_0_0[i] = -5.0 * tg_x_xyyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_x_xyyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xyyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyy_d_0_0_0[i] * a_z * faz_0;

        tg_xz_xyyyz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xyzzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_xzzzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_xzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyyy_d_0_0_0[i] = -5.0 * tg_z_yyyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_yyyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyy_d_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyyz_d_0_0_0[i] = -5.0 * tg_z_yyyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_yyyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyzz_d_0_0_0[i] = -5.0 * tg_z_yyyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_yyyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_yyzzz_d_0_0_0[i] = -5.0 * tg_z_yyzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_yyzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_yzzzz_d_0_0_0[i] = -5.0 * tg_z_yzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_yzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_xz_zzzzz_d_0_0_0[i] = -5.0 * tg_z_zzzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_z_zzzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_zzzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzz_d_0_0_0[i] * a_x * faz_0;

        tg_yy_xxxxx_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxx_d_0_0_0[i] * fzi_0 + tg_0_xxxxx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_xxxxx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxxx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxxx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxx_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxy_d_0_0_0[i] * fzi_0 + tg_0_xxxxy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_y_xxxx_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxxy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxxy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxxy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxy_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxz_d_0_0_0[i] * fzi_0 + tg_0_xxxxz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_xxxxz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxxz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxxz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyy_d_0_0_0[i] * fzi_0 + tg_0_xxxyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_y_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyy_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyz_d_0_0_0[i] * fzi_0 + tg_0_xxxyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_y_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxxyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxzz_d_0_0_0[i] * fzi_0 + tg_0_xxxzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_xxxzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxxzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyy_d_0_0_0[i] * fzi_0 + tg_0_xxyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_y_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxyyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxyyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyy_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyz_d_0_0_0[i] * fzi_0 + tg_0_xxyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_y_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyzz_d_0_0_0[i] * fzi_0 + tg_0_xxyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_y_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xxyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xxzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxzzz_d_0_0_0[i] * fzi_0 + tg_0_xxzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_xxzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xxzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyy_d_0_0_0[i] * fzi_0 + tg_0_xyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_y_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xyyyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyyyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyy_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyz_d_0_0_0[i] * fzi_0 + tg_0_xyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_y_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xyyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyzz_d_0_0_0[i] * fzi_0 + tg_0_xyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_y_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xyyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyzzz_d_0_0_0[i] * fzi_0 + tg_0_xyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_y_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_xyzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xyzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_xzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xzzzz_d_0_0_0[i] * fzi_0 + tg_0_xzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_xzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_xzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyy_d_0_0_0[i] * fzi_0 + tg_0_yyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 25.0 / 2.0 * tg_y_yyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_yyyyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_yyyyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyyyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyy_d_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyz_d_0_0_0[i] * fzi_0 + tg_0_yyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_y_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_yyyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_yyyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyzz_d_0_0_0[i] * fzi_0 + tg_0_yyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_y_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_yyyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_yyyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_yyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyzzz_d_0_0_0[i] * fzi_0 + tg_0_yyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_y_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_yyzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_yyzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_yzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yzzzz_d_0_0_0[i] * fzi_0 + tg_0_yzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_y_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_y_yzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_yzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yy_zzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_zzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_zzzzz_d_0_0_0[i] * fzi_0 + tg_0_zzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_y_zzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_y_zzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_zzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxx_d_0_0_0[i] = -5.0 * tg_z_xxxxx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxxxx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxxx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxx_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxy_d_0_0_0[i] = -5.0 * tg_y_xxxxy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_y_xxxxy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xxxxy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxy_d_0_0_0[i] * a_z * faz_0;

        tg_yz_xxxxz_d_0_0_0[i] = -5.0 * tg_z_xxxxz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxxxz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxxz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxyy_d_0_0_0[i] = -5.0 * tg_y_xxxyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_y_xxxyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xxxyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyy_d_0_0_0[i] * a_z * faz_0;

        tg_yz_xxxyz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxxyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxzz_d_0_0_0[i] = -5.0 * tg_z_xxxzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxxzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyyy_d_0_0_0[i] = -5.0 * tg_y_xxyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_y_xxyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xxyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyy_d_0_0_0[i] * a_z * faz_0;

        tg_yz_xxyyz_d_0_0_0[i] = 5.0 * tg_z_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xxzzz_d_0_0_0[i] = -5.0 * tg_z_xxzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xxzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyyy_d_0_0_0[i] = -5.0 * tg_y_xyyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_y_xyyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xyyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyy_d_0_0_0[i] * a_z * faz_0;

        tg_yz_xyyyz_d_0_0_0[i] = 15.0 / 2.0 * tg_z_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xyyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xyyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyzz_d_0_0_0[i] = 5.0 * tg_z_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xyyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xyyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xyzzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xyzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xyzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_xzzzz_d_0_0_0[i] = -5.0 * tg_z_xzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_xzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyyy_d_0_0_0[i] = -5.0 * tg_y_yyyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_y_yyyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_yyyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyy_d_0_0_0[i] * a_z * faz_0;

        tg_yz_yyyyz_d_0_0_0[i] = 10.0 * tg_z_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_yyyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yyyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyzz_d_0_0_0[i] = 15.0 / 2.0 * tg_z_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_yyyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yyyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_yyzzz_d_0_0_0[i] = 5.0 * tg_z_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_yyzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yyzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_yzzzz_d_0_0_0[i] = 5.0 / 2.0 * tg_z_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_yzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_yz_zzzzz_d_0_0_0[i] = -5.0 * tg_z_zzzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_z_zzzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_zzzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzz_d_0_0_0[i] * a_y * faz_0;

        tg_zz_xxxxx_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxx_d_0_0_0[i] * fzi_0 + tg_0_xxxxx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_xxxxx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxxx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxxx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxx_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxy_d_0_0_0[i] * fzi_0 + tg_0_xxxxy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_xxxxy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxxy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxxy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxy_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxxz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxxz_d_0_0_0[i] * fzi_0 + tg_0_xxxxz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_z_xxxx_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxxz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxxz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxxz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyy_d_0_0_0[i] * fzi_0 + tg_0_xxxyy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_xxxyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyy_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxyz_d_0_0_0[i] * fzi_0 + tg_0_xxxyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_z_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxxzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxxzz_d_0_0_0[i] * fzi_0 + tg_0_xxxzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_z_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxxzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxxzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyy_d_0_0_0[i] * fzi_0 + tg_0_xxyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_xxyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyy_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyyz_d_0_0_0[i] * fzi_0 + tg_0_xxyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_z_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxyyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxyyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxyzz_d_0_0_0[i] * fzi_0 + tg_0_xxyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_z_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxyzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxyzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxyzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xxzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xxzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xxzzz_d_0_0_0[i] * fzi_0 + tg_0_xxzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_z_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xxzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xxzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyy_d_0_0_0[i] * fzi_0 + tg_0_xyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_xyyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xyyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyy_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyyz_d_0_0_0[i] * fzi_0 + tg_0_xyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_z_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xyyyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyyyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyyzz_d_0_0_0[i] * fzi_0 + tg_0_xyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_z_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyyzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xyyzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyyzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xyzzz_d_0_0_0[i] * fzi_0 + tg_0_xyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_z_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xyzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xyzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_xzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_xzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_xzzzz_d_0_0_0[i] * fzi_0 + tg_0_xzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_z_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_xzzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_xzzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xzzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyyy_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyy_d_0_0_0[i] * fzi_0 + tg_0_yyyyy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_z_yyyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_yyyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyy_d_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyyz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyyz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyyz_d_0_0_0[i] * fzi_0 + tg_0_yyyyz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_z_yyyy_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyyyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_yyyyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyyyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyyzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyyzz_d_0_0_0[i] * fzi_0 + tg_0_yyyzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_z_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyyzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_yyyzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyyzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_yyzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yyzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yyzzz_d_0_0_0[i] * fzi_0 + tg_0_yyzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 15.0 / 2.0 * tg_z_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yyzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_yyzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_yzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_yzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_yzzzz_d_0_0_0[i] * fzi_0 + tg_0_yzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 10.0 * tg_z_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_yzzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_yzzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yzzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzz_d_0_0_0[i] * a_z * faz_0;

        tg_zz_zzzzz_d_0_0_0[i] = -5.0 / 2.0 * tg_0_zzzzz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_zzzzz_d_0_0_0[i] * fzi_0 + tg_0_zzzzz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 25.0 / 2.0 * tg_z_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_z_zzzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_z_zzzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_zzzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzz_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SH

        auto tg_0_xxxxx_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1);

        auto tg_0_xxxxy_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 1);

        auto tg_0_xxxxz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 2);

        auto tg_0_xxxyy_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 3);

        auto tg_0_xxxyz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 4);

        auto tg_0_xxxzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 5);

        auto tg_0_xxyyy_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 6);

        auto tg_0_xxyyz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 7);

        auto tg_0_xxyzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 8);

        auto tg_0_xxzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 9);

        auto tg_0_xyyyy_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 10);

        auto tg_0_xyyyz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 11);

        auto tg_0_xyyzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 12);

        auto tg_0_xyzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 13);

        auto tg_0_xzzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 14);

        auto tg_0_yyyyy_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 15);

        auto tg_0_yyyyz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 16);

        auto tg_0_yyyzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 17);

        auto tg_0_yyzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 18);

        auto tg_0_yzzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 19);

        auto tg_0_zzzzz_d_0_0_1 = pbuffer.data(idx_sh_d_0_0_1 + 20);

        // Set up components of auxiliary buffer : PH

        auto tg_x_xxxxx_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1);

        auto tg_x_xxxxy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 1);

        auto tg_x_xxxxz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 2);

        auto tg_x_xxxyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 3);

        auto tg_x_xxxyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 4);

        auto tg_x_xxxzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 5);

        auto tg_x_xxyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 6);

        auto tg_x_xxyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 7);

        auto tg_x_xxyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 8);

        auto tg_x_xxzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 9);

        auto tg_x_xyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 10);

        auto tg_x_xyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 11);

        auto tg_x_xyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 12);

        auto tg_x_xyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 13);

        auto tg_x_xzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 14);

        auto tg_x_yyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 15);

        auto tg_x_yyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 16);

        auto tg_x_yyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 17);

        auto tg_x_yyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 18);

        auto tg_x_yzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 19);

        auto tg_x_zzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 20);

        auto tg_y_xxxxx_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 21);

        auto tg_y_xxxxy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 22);

        auto tg_y_xxxxz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 23);

        auto tg_y_xxxyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 24);

        auto tg_y_xxxyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 25);

        auto tg_y_xxxzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 26);

        auto tg_y_xxyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 27);

        auto tg_y_xxyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 28);

        auto tg_y_xxyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 29);

        auto tg_y_xxzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 30);

        auto tg_y_xyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 31);

        auto tg_y_xyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 32);

        auto tg_y_xyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 33);

        auto tg_y_xyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 34);

        auto tg_y_xzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 35);

        auto tg_y_yyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 36);

        auto tg_y_yyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 37);

        auto tg_y_yyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 38);

        auto tg_y_yyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 39);

        auto tg_y_yzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 40);

        auto tg_y_zzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 41);

        auto tg_z_xxxxx_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 42);

        auto tg_z_xxxxy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 43);

        auto tg_z_xxxxz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 44);

        auto tg_z_xxxyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 45);

        auto tg_z_xxxyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 46);

        auto tg_z_xxxzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 47);

        auto tg_z_xxyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 48);

        auto tg_z_xxyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 49);

        auto tg_z_xxyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 50);

        auto tg_z_xxzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 51);

        auto tg_z_xyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 52);

        auto tg_z_xyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 53);

        auto tg_z_xyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 54);

        auto tg_z_xyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 55);

        auto tg_z_xzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 56);

        auto tg_z_yyyyy_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 57);

        auto tg_z_yyyyz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 58);

        auto tg_z_yyyzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 59);

        auto tg_z_yyzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 60);

        auto tg_z_yzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 61);

        auto tg_z_zzzzz_d_0_0_1 = pbuffer.data(idx_ph_d_0_0_1 + 62);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxxx_d_0_0_1, tg_0_xxxxy_d_0_0_1, tg_0_xxxxz_d_0_0_1, tg_0_xxxyy_d_0_0_1, tg_0_xxxyz_d_0_0_1, tg_0_xxxzz_d_0_0_1, tg_0_xxyyy_d_0_0_1, tg_0_xxyyz_d_0_0_1, tg_0_xxyzz_d_0_0_1, tg_0_xxzzz_d_0_0_1, tg_0_xyyyy_d_0_0_1, tg_0_xyyyz_d_0_0_1, tg_0_xyyzz_d_0_0_1, tg_0_xyzzz_d_0_0_1, tg_0_xzzzz_d_0_0_1, tg_0_yyyyy_d_0_0_1, tg_0_yyyyz_d_0_0_1, tg_0_yyyzz_d_0_0_1, tg_0_yyzzz_d_0_0_1, tg_0_yzzzz_d_0_0_1, tg_0_zzzzz_d_0_0_1, tg_x_xxxxx_d_0_0_1, tg_x_xxxxy_d_0_0_1, tg_x_xxxxz_d_0_0_1, tg_x_xxxyy_d_0_0_1, tg_x_xxxyz_d_0_0_1, tg_x_xxxzz_d_0_0_1, tg_x_xxyyy_d_0_0_1, tg_x_xxyyz_d_0_0_1, tg_x_xxyzz_d_0_0_1, tg_x_xxzzz_d_0_0_1, tg_x_xyyyy_d_0_0_1, tg_x_xyyyz_d_0_0_1, tg_x_xyyzz_d_0_0_1, tg_x_xyzzz_d_0_0_1, tg_x_xzzzz_d_0_0_1, tg_x_yyyyy_d_0_0_1, tg_x_yyyyz_d_0_0_1, tg_x_yyyzz_d_0_0_1, tg_x_yyzzz_d_0_0_1, tg_x_yzzzz_d_0_0_1, tg_x_zzzzz_d_0_0_1, tg_xx_xxxxx_d_0_0_0, tg_xx_xxxxy_d_0_0_0, tg_xx_xxxxz_d_0_0_0, tg_xx_xxxyy_d_0_0_0, tg_xx_xxxyz_d_0_0_0, tg_xx_xxxzz_d_0_0_0, tg_xx_xxyyy_d_0_0_0, tg_xx_xxyyz_d_0_0_0, tg_xx_xxyzz_d_0_0_0, tg_xx_xxzzz_d_0_0_0, tg_xx_xyyyy_d_0_0_0, tg_xx_xyyyz_d_0_0_0, tg_xx_xyyzz_d_0_0_0, tg_xx_xyzzz_d_0_0_0, tg_xx_xzzzz_d_0_0_0, tg_xx_yyyyy_d_0_0_0, tg_xx_yyyyz_d_0_0_0, tg_xx_yyyzz_d_0_0_0, tg_xx_yyzzz_d_0_0_0, tg_xx_yzzzz_d_0_0_0, tg_xx_zzzzz_d_0_0_0, tg_xy_xxxxx_d_0_0_0, tg_xy_xxxxy_d_0_0_0, tg_xy_xxxxz_d_0_0_0, tg_xy_xxxyy_d_0_0_0, tg_xy_xxxyz_d_0_0_0, tg_xy_xxxzz_d_0_0_0, tg_xy_xxyyy_d_0_0_0, tg_xy_xxyyz_d_0_0_0, tg_xy_xxyzz_d_0_0_0, tg_xy_xxzzz_d_0_0_0, tg_xy_xyyyy_d_0_0_0, tg_xy_xyyyz_d_0_0_0, tg_xy_xyyzz_d_0_0_0, tg_xy_xyzzz_d_0_0_0, tg_xy_xzzzz_d_0_0_0, tg_xy_yyyyy_d_0_0_0, tg_xy_yyyyz_d_0_0_0, tg_xy_yyyzz_d_0_0_0, tg_xy_yyzzz_d_0_0_0, tg_xy_yzzzz_d_0_0_0, tg_xy_zzzzz_d_0_0_0, tg_xz_xxxxx_d_0_0_0, tg_xz_xxxxy_d_0_0_0, tg_xz_xxxxz_d_0_0_0, tg_xz_xxxyy_d_0_0_0, tg_xz_xxxyz_d_0_0_0, tg_xz_xxxzz_d_0_0_0, tg_xz_xxyyy_d_0_0_0, tg_xz_xxyyz_d_0_0_0, tg_xz_xxyzz_d_0_0_0, tg_xz_xxzzz_d_0_0_0, tg_xz_xyyyy_d_0_0_0, tg_xz_xyyyz_d_0_0_0, tg_xz_xyyzz_d_0_0_0, tg_xz_xyzzz_d_0_0_0, tg_xz_xzzzz_d_0_0_0, tg_xz_yyyyy_d_0_0_0, tg_xz_yyyyz_d_0_0_0, tg_xz_yyyzz_d_0_0_0, tg_xz_yyzzz_d_0_0_0, tg_xz_yzzzz_d_0_0_0, tg_xz_zzzzz_d_0_0_0, tg_y_xxxxx_d_0_0_1, tg_y_xxxxy_d_0_0_1, tg_y_xxxxz_d_0_0_1, tg_y_xxxyy_d_0_0_1, tg_y_xxxyz_d_0_0_1, tg_y_xxxzz_d_0_0_1, tg_y_xxyyy_d_0_0_1, tg_y_xxyyz_d_0_0_1, tg_y_xxyzz_d_0_0_1, tg_y_xxzzz_d_0_0_1, tg_y_xyyyy_d_0_0_1, tg_y_xyyyz_d_0_0_1, tg_y_xyyzz_d_0_0_1, tg_y_xyzzz_d_0_0_1, tg_y_xzzzz_d_0_0_1, tg_y_yyyyy_d_0_0_1, tg_y_yyyyz_d_0_0_1, tg_y_yyyzz_d_0_0_1, tg_y_yyzzz_d_0_0_1, tg_y_yzzzz_d_0_0_1, tg_y_zzzzz_d_0_0_1, tg_yy_xxxxx_d_0_0_0, tg_yy_xxxxy_d_0_0_0, tg_yy_xxxxz_d_0_0_0, tg_yy_xxxyy_d_0_0_0, tg_yy_xxxyz_d_0_0_0, tg_yy_xxxzz_d_0_0_0, tg_yy_xxyyy_d_0_0_0, tg_yy_xxyyz_d_0_0_0, tg_yy_xxyzz_d_0_0_0, tg_yy_xxzzz_d_0_0_0, tg_yy_xyyyy_d_0_0_0, tg_yy_xyyyz_d_0_0_0, tg_yy_xyyzz_d_0_0_0, tg_yy_xyzzz_d_0_0_0, tg_yy_xzzzz_d_0_0_0, tg_yy_yyyyy_d_0_0_0, tg_yy_yyyyz_d_0_0_0, tg_yy_yyyzz_d_0_0_0, tg_yy_yyzzz_d_0_0_0, tg_yy_yzzzz_d_0_0_0, tg_yy_zzzzz_d_0_0_0, tg_yz_xxxxx_d_0_0_0, tg_yz_xxxxy_d_0_0_0, tg_yz_xxxxz_d_0_0_0, tg_yz_xxxyy_d_0_0_0, tg_yz_xxxyz_d_0_0_0, tg_yz_xxxzz_d_0_0_0, tg_yz_xxyyy_d_0_0_0, tg_yz_xxyyz_d_0_0_0, tg_yz_xxyzz_d_0_0_0, tg_yz_xxzzz_d_0_0_0, tg_yz_xyyyy_d_0_0_0, tg_yz_xyyyz_d_0_0_0, tg_yz_xyyzz_d_0_0_0, tg_yz_xyzzz_d_0_0_0, tg_yz_xzzzz_d_0_0_0, tg_yz_yyyyy_d_0_0_0, tg_yz_yyyyz_d_0_0_0, tg_yz_yyyzz_d_0_0_0, tg_yz_yyzzz_d_0_0_0, tg_yz_yzzzz_d_0_0_0, tg_yz_zzzzz_d_0_0_0, tg_z_xxxxx_d_0_0_1, tg_z_xxxxy_d_0_0_1, tg_z_xxxxz_d_0_0_1, tg_z_xxxyy_d_0_0_1, tg_z_xxxyz_d_0_0_1, tg_z_xxxzz_d_0_0_1, tg_z_xxyyy_d_0_0_1, tg_z_xxyyz_d_0_0_1, tg_z_xxyzz_d_0_0_1, tg_z_xxzzz_d_0_0_1, tg_z_xyyyy_d_0_0_1, tg_z_xyyyz_d_0_0_1, tg_z_xyyzz_d_0_0_1, tg_z_xyzzz_d_0_0_1, tg_z_xzzzz_d_0_0_1, tg_z_yyyyy_d_0_0_1, tg_z_yyyyz_d_0_0_1, tg_z_yyyzz_d_0_0_1, tg_z_yyzzz_d_0_0_1, tg_z_yzzzz_d_0_0_1, tg_z_zzzzz_d_0_0_1, tg_zz_xxxxx_d_0_0_0, tg_zz_xxxxy_d_0_0_0, tg_zz_xxxxz_d_0_0_0, tg_zz_xxxyy_d_0_0_0, tg_zz_xxxyz_d_0_0_0, tg_zz_xxxzz_d_0_0_0, tg_zz_xxyyy_d_0_0_0, tg_zz_xxyyz_d_0_0_0, tg_zz_xxyzz_d_0_0_0, tg_zz_xxzzz_d_0_0_0, tg_zz_xyyyy_d_0_0_0, tg_zz_xyyyz_d_0_0_0, tg_zz_xyyzz_d_0_0_0, tg_zz_xyzzz_d_0_0_0, tg_zz_xzzzz_d_0_0_0, tg_zz_yyyyy_d_0_0_0, tg_zz_yyyyz_d_0_0_0, tg_zz_yyyzz_d_0_0_0, tg_zz_yyzzz_d_0_0_0, tg_zz_yzzzz_d_0_0_0, tg_zz_zzzzz_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_xxxxx_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_zzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_zzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxx_d_0_0_0[i] += tg_y_xxxxx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxy_d_0_0_0[i] += tg_y_xxxxy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxz_d_0_0_0[i] += tg_y_xxxxz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxyy_d_0_0_0[i] += tg_y_xxxyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxyz_d_0_0_0[i] += tg_y_xxxyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxzz_d_0_0_0[i] += tg_y_xxxzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyyy_d_0_0_0[i] += tg_y_xxyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyyz_d_0_0_0[i] += tg_y_xxyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyzz_d_0_0_0[i] += tg_y_xxyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxzzz_d_0_0_0[i] += tg_y_xxzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyyy_d_0_0_0[i] += tg_y_xyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyyz_d_0_0_0[i] += tg_y_xyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyzz_d_0_0_0[i] += tg_y_xyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyzzz_d_0_0_0[i] += tg_y_xyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xzzzz_d_0_0_0[i] += tg_y_xzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyyy_d_0_0_0[i] += tg_y_yyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyyz_d_0_0_0[i] += tg_y_yyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyzz_d_0_0_0[i] += tg_y_yyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyzzz_d_0_0_0[i] += tg_y_yyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yzzzz_d_0_0_0[i] += tg_y_yzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_zzzzz_d_0_0_0[i] += tg_y_zzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxx_d_0_0_0[i] += tg_z_xxxxx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxy_d_0_0_0[i] += tg_z_xxxxy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxz_d_0_0_0[i] += tg_z_xxxxz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxyy_d_0_0_0[i] += tg_z_xxxyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxyz_d_0_0_0[i] += tg_z_xxxyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxzz_d_0_0_0[i] += tg_z_xxxzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyyy_d_0_0_0[i] += tg_z_xxyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyyz_d_0_0_0[i] += tg_z_xxyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyzz_d_0_0_0[i] += tg_z_xxyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxzzz_d_0_0_0[i] += tg_z_xxzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyyy_d_0_0_0[i] += tg_z_xyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyyz_d_0_0_0[i] += tg_z_xyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyzz_d_0_0_0[i] += tg_z_xyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyzzz_d_0_0_0[i] += tg_z_xyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xzzzz_d_0_0_0[i] += tg_z_xzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyyy_d_0_0_0[i] += tg_z_yyyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyyz_d_0_0_0[i] += tg_z_yyyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyzz_d_0_0_0[i] += tg_z_yyyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyzzz_d_0_0_0[i] += tg_z_yyzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yzzzz_d_0_0_0[i] += tg_z_yzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_zzzzz_d_0_0_0[i] += tg_z_zzzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_xxxxx_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_zzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_zzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxx_d_0_0_0[i] += tg_z_xxxxx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxy_d_0_0_0[i] += tg_z_xxxxy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxz_d_0_0_0[i] += tg_z_xxxxz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxyy_d_0_0_0[i] += tg_z_xxxyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxyz_d_0_0_0[i] += tg_z_xxxyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxzz_d_0_0_0[i] += tg_z_xxxzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyyy_d_0_0_0[i] += tg_z_xxyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyyz_d_0_0_0[i] += tg_z_xxyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyzz_d_0_0_0[i] += tg_z_xxyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxzzz_d_0_0_0[i] += tg_z_xxzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyyy_d_0_0_0[i] += tg_z_xyyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyyz_d_0_0_0[i] += tg_z_xyyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyzz_d_0_0_0[i] += tg_z_xyyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyzzz_d_0_0_0[i] += tg_z_xyzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xzzzz_d_0_0_0[i] += tg_z_xzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyyy_d_0_0_0[i] += tg_z_yyyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyyz_d_0_0_0[i] += tg_z_yyyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyzz_d_0_0_0[i] += tg_z_yyyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyzzz_d_0_0_0[i] += tg_z_yyzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yzzzz_d_0_0_0[i] += tg_z_yzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_zzzzz_d_0_0_0[i] += tg_z_zzzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_xxxxx_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xzzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyyy_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyyz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yzzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_zzzzz_d_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_zzzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

