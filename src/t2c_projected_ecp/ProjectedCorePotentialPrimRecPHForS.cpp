#include "ProjectedCorePotentialPrimRecPHForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ph_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_ph_s_0_0_0,
                                        const size_t idx_sh_s_0_0_0,
                                        const size_t idx_sh_s_1_0_0,
                                        const int p,
                                        const size_t idx_sh_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0);

    auto tg_0_xxxxy_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 1);

    auto tg_0_xxxxz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 2);

    auto tg_0_xxxyy_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 3);

    auto tg_0_xxxyz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 4);

    auto tg_0_xxxzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 5);

    auto tg_0_xxyyy_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 6);

    auto tg_0_xxyyz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 7);

    auto tg_0_xxyzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 8);

    auto tg_0_xxzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 9);

    auto tg_0_xyyyy_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 10);

    auto tg_0_xyyyz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 11);

    auto tg_0_xyyzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 12);

    auto tg_0_xyzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 13);

    auto tg_0_xzzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 14);

    auto tg_0_yyyyy_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 15);

    auto tg_0_yyyyz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 16);

    auto tg_0_yyyzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 17);

    auto tg_0_yyzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 18);

    auto tg_0_yzzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 19);

    auto tg_0_zzzzz_s_0_0_0 = pbuffer.data(idx_sh_s_0_0_0 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0);

    auto tg_0_xxxxy_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 1);

    auto tg_0_xxxxz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 2);

    auto tg_0_xxxyy_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 3);

    auto tg_0_xxxyz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 4);

    auto tg_0_xxxzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 5);

    auto tg_0_xxyyy_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 6);

    auto tg_0_xxyyz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 7);

    auto tg_0_xxyzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 8);

    auto tg_0_xxzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 9);

    auto tg_0_xyyyy_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 10);

    auto tg_0_xyyyz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 11);

    auto tg_0_xyyzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 12);

    auto tg_0_xyzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 13);

    auto tg_0_xzzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 14);

    auto tg_0_yyyyy_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 15);

    auto tg_0_yyyyz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 16);

    auto tg_0_yyyzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 17);

    auto tg_0_yyzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 18);

    auto tg_0_yzzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 19);

    auto tg_0_zzzzz_s_1_0_0 = pbuffer.data(idx_sh_s_1_0_0 + 20);

    // Set up components of targeted buffer : PH

    auto tg_x_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0);

    auto tg_x_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 1);

    auto tg_x_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 2);

    auto tg_x_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 3);

    auto tg_x_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 4);

    auto tg_x_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 5);

    auto tg_x_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 6);

    auto tg_x_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 7);

    auto tg_x_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 8);

    auto tg_x_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 9);

    auto tg_x_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 10);

    auto tg_x_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 11);

    auto tg_x_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 12);

    auto tg_x_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 13);

    auto tg_x_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 14);

    auto tg_x_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 15);

    auto tg_x_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 16);

    auto tg_x_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 17);

    auto tg_x_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 18);

    auto tg_x_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 19);

    auto tg_x_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 20);

    auto tg_y_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 21);

    auto tg_y_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 22);

    auto tg_y_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 23);

    auto tg_y_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 24);

    auto tg_y_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 25);

    auto tg_y_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 26);

    auto tg_y_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 27);

    auto tg_y_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 28);

    auto tg_y_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 29);

    auto tg_y_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 30);

    auto tg_y_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 31);

    auto tg_y_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 32);

    auto tg_y_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 33);

    auto tg_y_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 34);

    auto tg_y_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 35);

    auto tg_y_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 36);

    auto tg_y_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 37);

    auto tg_y_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 38);

    auto tg_y_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 39);

    auto tg_y_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 40);

    auto tg_y_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 41);

    auto tg_z_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 42);

    auto tg_z_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 43);

    auto tg_z_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 44);

    auto tg_z_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 45);

    auto tg_z_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 46);

    auto tg_z_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 47);

    auto tg_z_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 48);

    auto tg_z_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 49);

    auto tg_z_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 50);

    auto tg_z_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 51);

    auto tg_z_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 52);

    auto tg_z_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 53);

    auto tg_z_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 54);

    auto tg_z_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 55);

    auto tg_z_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 56);

    auto tg_z_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 57);

    auto tg_z_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 58);

    auto tg_z_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 59);

    auto tg_z_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 60);

    auto tg_z_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 61);

    auto tg_z_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 62);

    #pragma omp simd aligned(b_exps, tg_0_xxxxx_s_0_0_0, tg_0_xxxxx_s_1_0_0, tg_0_xxxxy_s_0_0_0, tg_0_xxxxy_s_1_0_0, tg_0_xxxxz_s_0_0_0, tg_0_xxxxz_s_1_0_0, tg_0_xxxyy_s_0_0_0, tg_0_xxxyy_s_1_0_0, tg_0_xxxyz_s_0_0_0, tg_0_xxxyz_s_1_0_0, tg_0_xxxzz_s_0_0_0, tg_0_xxxzz_s_1_0_0, tg_0_xxyyy_s_0_0_0, tg_0_xxyyy_s_1_0_0, tg_0_xxyyz_s_0_0_0, tg_0_xxyyz_s_1_0_0, tg_0_xxyzz_s_0_0_0, tg_0_xxyzz_s_1_0_0, tg_0_xxzzz_s_0_0_0, tg_0_xxzzz_s_1_0_0, tg_0_xyyyy_s_0_0_0, tg_0_xyyyy_s_1_0_0, tg_0_xyyyz_s_0_0_0, tg_0_xyyyz_s_1_0_0, tg_0_xyyzz_s_0_0_0, tg_0_xyyzz_s_1_0_0, tg_0_xyzzz_s_0_0_0, tg_0_xyzzz_s_1_0_0, tg_0_xzzzz_s_0_0_0, tg_0_xzzzz_s_1_0_0, tg_0_yyyyy_s_0_0_0, tg_0_yyyyy_s_1_0_0, tg_0_yyyyz_s_0_0_0, tg_0_yyyyz_s_1_0_0, tg_0_yyyzz_s_0_0_0, tg_0_yyyzz_s_1_0_0, tg_0_yyzzz_s_0_0_0, tg_0_yyzzz_s_1_0_0, tg_0_yzzzz_s_0_0_0, tg_0_yzzzz_s_1_0_0, tg_0_zzzzz_s_0_0_0, tg_0_zzzzz_s_1_0_0, tg_x_xxxxx_s_0_0_0, tg_x_xxxxy_s_0_0_0, tg_x_xxxxz_s_0_0_0, tg_x_xxxyy_s_0_0_0, tg_x_xxxyz_s_0_0_0, tg_x_xxxzz_s_0_0_0, tg_x_xxyyy_s_0_0_0, tg_x_xxyyz_s_0_0_0, tg_x_xxyzz_s_0_0_0, tg_x_xxzzz_s_0_0_0, tg_x_xyyyy_s_0_0_0, tg_x_xyyyz_s_0_0_0, tg_x_xyyzz_s_0_0_0, tg_x_xyzzz_s_0_0_0, tg_x_xzzzz_s_0_0_0, tg_x_yyyyy_s_0_0_0, tg_x_yyyyz_s_0_0_0, tg_x_yyyzz_s_0_0_0, tg_x_yyzzz_s_0_0_0, tg_x_yzzzz_s_0_0_0, tg_x_zzzzz_s_0_0_0, tg_y_xxxxx_s_0_0_0, tg_y_xxxxy_s_0_0_0, tg_y_xxxxz_s_0_0_0, tg_y_xxxyy_s_0_0_0, tg_y_xxxyz_s_0_0_0, tg_y_xxxzz_s_0_0_0, tg_y_xxyyy_s_0_0_0, tg_y_xxyyz_s_0_0_0, tg_y_xxyzz_s_0_0_0, tg_y_xxzzz_s_0_0_0, tg_y_xyyyy_s_0_0_0, tg_y_xyyyz_s_0_0_0, tg_y_xyyzz_s_0_0_0, tg_y_xyzzz_s_0_0_0, tg_y_xzzzz_s_0_0_0, tg_y_yyyyy_s_0_0_0, tg_y_yyyyz_s_0_0_0, tg_y_yyyzz_s_0_0_0, tg_y_yyzzz_s_0_0_0, tg_y_yzzzz_s_0_0_0, tg_y_zzzzz_s_0_0_0, tg_z_xxxxx_s_0_0_0, tg_z_xxxxy_s_0_0_0, tg_z_xxxxz_s_0_0_0, tg_z_xxxyy_s_0_0_0, tg_z_xxxyz_s_0_0_0, tg_z_xxxzz_s_0_0_0, tg_z_xxyyy_s_0_0_0, tg_z_xxyyz_s_0_0_0, tg_z_xxyzz_s_0_0_0, tg_z_xxzzz_s_0_0_0, tg_z_xyyyy_s_0_0_0, tg_z_xyyyz_s_0_0_0, tg_z_xyyzz_s_0_0_0, tg_z_xyzzz_s_0_0_0, tg_z_xzzzz_s_0_0_0, tg_z_yyyyy_s_0_0_0, tg_z_yyyyz_s_0_0_0, tg_z_yyyzz_s_0_0_0, tg_z_yyzzz_s_0_0_0, tg_z_yzzzz_s_0_0_0, tg_z_zzzzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_x_xxxxx_s_0_0_0[i] = 2.0 * tg_0_xxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxy_s_0_0_0[i] = 2.0 * tg_0_xxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxz_s_0_0_0[i] = 2.0 * tg_0_xxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyy_s_0_0_0[i] = 2.0 * tg_0_xxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyz_s_0_0_0[i] = 2.0 * tg_0_xxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxxzz_s_0_0_0[i] = 2.0 * tg_0_xxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyy_s_0_0_0[i] = 2.0 * tg_0_xxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyz_s_0_0_0[i] = 2.0 * tg_0_xxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxyzz_s_0_0_0[i] = 2.0 * tg_0_xxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xxzzz_s_0_0_0[i] = 2.0 * tg_0_xxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyy_s_0_0_0[i] = 2.0 * tg_0_xyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyz_s_0_0_0[i] = 2.0 * tg_0_xyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xyyzz_s_0_0_0[i] = 2.0 * tg_0_xyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xyzzz_s_0_0_0[i] = 2.0 * tg_0_xyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_xzzzz_s_0_0_0[i] = 2.0 * tg_0_xzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyy_s_0_0_0[i] = 2.0 * tg_0_yyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyz_s_0_0_0[i] = 2.0 * tg_0_yyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_x_yyyzz_s_0_0_0[i] = 2.0 * tg_0_yyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_yyzzz_s_0_0_0[i] = 2.0 * tg_0_yyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_yzzzz_s_0_0_0[i] = 2.0 * tg_0_yzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_x_zzzzz_s_0_0_0[i] = 2.0 * tg_0_zzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_y_xxxxx_s_0_0_0[i] = 2.0 * tg_0_xxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxy_s_0_0_0[i] = 2.0 * tg_0_xxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxz_s_0_0_0[i] = 2.0 * tg_0_xxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyy_s_0_0_0[i] = 2.0 * tg_0_xxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyz_s_0_0_0[i] = 2.0 * tg_0_xxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxxzz_s_0_0_0[i] = 2.0 * tg_0_xxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyy_s_0_0_0[i] = 2.0 * tg_0_xxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyz_s_0_0_0[i] = 2.0 * tg_0_xxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxyzz_s_0_0_0[i] = 2.0 * tg_0_xxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xxzzz_s_0_0_0[i] = 2.0 * tg_0_xxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyy_s_0_0_0[i] = 2.0 * tg_0_xyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyz_s_0_0_0[i] = 2.0 * tg_0_xyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xyyzz_s_0_0_0[i] = 2.0 * tg_0_xyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xyzzz_s_0_0_0[i] = 2.0 * tg_0_xyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_xzzzz_s_0_0_0[i] = 2.0 * tg_0_xzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyy_s_0_0_0[i] = 2.0 * tg_0_yyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyz_s_0_0_0[i] = 2.0 * tg_0_yyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_y_yyyzz_s_0_0_0[i] = 2.0 * tg_0_yyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_yyzzz_s_0_0_0[i] = 2.0 * tg_0_yyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_yzzzz_s_0_0_0[i] = 2.0 * tg_0_yzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_y_zzzzz_s_0_0_0[i] = 2.0 * tg_0_zzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_z_xxxxx_s_0_0_0[i] = 2.0 * tg_0_xxxxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxy_s_0_0_0[i] = 2.0 * tg_0_xxxxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxz_s_0_0_0[i] = 2.0 * tg_0_xxxxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyy_s_0_0_0[i] = 2.0 * tg_0_xxxyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyz_s_0_0_0[i] = 2.0 * tg_0_xxxyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxxzz_s_0_0_0[i] = 2.0 * tg_0_xxxzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyy_s_0_0_0[i] = 2.0 * tg_0_xxyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyz_s_0_0_0[i] = 2.0 * tg_0_xxyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxyzz_s_0_0_0[i] = 2.0 * tg_0_xxyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xxzzz_s_0_0_0[i] = 2.0 * tg_0_xxzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyy_s_0_0_0[i] = 2.0 * tg_0_xyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyz_s_0_0_0[i] = 2.0 * tg_0_xyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xyyzz_s_0_0_0[i] = 2.0 * tg_0_xyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xyzzz_s_0_0_0[i] = 2.0 * tg_0_xyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_xzzzz_s_0_0_0[i] = 2.0 * tg_0_xzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyy_s_0_0_0[i] = 2.0 * tg_0_yyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyz_s_0_0_0[i] = 2.0 * tg_0_yyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_z_yyyzz_s_0_0_0[i] = 2.0 * tg_0_yyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_yyzzz_s_0_0_0[i] = 2.0 * tg_0_yyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_yzzzz_s_0_0_0[i] = 2.0 * tg_0_yzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_z_zzzzz_s_0_0_0[i] = 2.0 * tg_0_zzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SH

        auto tg_0_xxxxx_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1);

        auto tg_0_xxxxy_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 1);

        auto tg_0_xxxxz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 2);

        auto tg_0_xxxyy_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 3);

        auto tg_0_xxxyz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 4);

        auto tg_0_xxxzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 5);

        auto tg_0_xxyyy_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 6);

        auto tg_0_xxyyz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 7);

        auto tg_0_xxyzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 8);

        auto tg_0_xxzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 9);

        auto tg_0_xyyyy_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 10);

        auto tg_0_xyyyz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 11);

        auto tg_0_xyyzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 12);

        auto tg_0_xyzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 13);

        auto tg_0_xzzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 14);

        auto tg_0_yyyyy_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 15);

        auto tg_0_yyyyz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 16);

        auto tg_0_yyyzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 17);

        auto tg_0_yyzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 18);

        auto tg_0_yzzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 19);

        auto tg_0_zzzzz_s_0_0_1 = pbuffer.data(idx_sh_s_0_0_1 + 20);

        #pragma omp simd aligned(b_exps, tg_0_xxxxx_s_0_0_1, tg_0_xxxxy_s_0_0_1, tg_0_xxxxz_s_0_0_1, tg_0_xxxyy_s_0_0_1, tg_0_xxxyz_s_0_0_1, tg_0_xxxzz_s_0_0_1, tg_0_xxyyy_s_0_0_1, tg_0_xxyyz_s_0_0_1, tg_0_xxyzz_s_0_0_1, tg_0_xxzzz_s_0_0_1, tg_0_xyyyy_s_0_0_1, tg_0_xyyyz_s_0_0_1, tg_0_xyyzz_s_0_0_1, tg_0_xyzzz_s_0_0_1, tg_0_xzzzz_s_0_0_1, tg_0_yyyyy_s_0_0_1, tg_0_yyyyz_s_0_0_1, tg_0_yyyzz_s_0_0_1, tg_0_yyzzz_s_0_0_1, tg_0_yzzzz_s_0_0_1, tg_0_zzzzz_s_0_0_1, tg_x_xxxxx_s_0_0_0, tg_x_xxxxy_s_0_0_0, tg_x_xxxxz_s_0_0_0, tg_x_xxxyy_s_0_0_0, tg_x_xxxyz_s_0_0_0, tg_x_xxxzz_s_0_0_0, tg_x_xxyyy_s_0_0_0, tg_x_xxyyz_s_0_0_0, tg_x_xxyzz_s_0_0_0, tg_x_xxzzz_s_0_0_0, tg_x_xyyyy_s_0_0_0, tg_x_xyyyz_s_0_0_0, tg_x_xyyzz_s_0_0_0, tg_x_xyzzz_s_0_0_0, tg_x_xzzzz_s_0_0_0, tg_x_yyyyy_s_0_0_0, tg_x_yyyyz_s_0_0_0, tg_x_yyyzz_s_0_0_0, tg_x_yyzzz_s_0_0_0, tg_x_yzzzz_s_0_0_0, tg_x_zzzzz_s_0_0_0, tg_y_xxxxx_s_0_0_0, tg_y_xxxxy_s_0_0_0, tg_y_xxxxz_s_0_0_0, tg_y_xxxyy_s_0_0_0, tg_y_xxxyz_s_0_0_0, tg_y_xxxzz_s_0_0_0, tg_y_xxyyy_s_0_0_0, tg_y_xxyyz_s_0_0_0, tg_y_xxyzz_s_0_0_0, tg_y_xxzzz_s_0_0_0, tg_y_xyyyy_s_0_0_0, tg_y_xyyyz_s_0_0_0, tg_y_xyyzz_s_0_0_0, tg_y_xyzzz_s_0_0_0, tg_y_xzzzz_s_0_0_0, tg_y_yyyyy_s_0_0_0, tg_y_yyyyz_s_0_0_0, tg_y_yyyzz_s_0_0_0, tg_y_yyzzz_s_0_0_0, tg_y_yzzzz_s_0_0_0, tg_y_zzzzz_s_0_0_0, tg_z_xxxxx_s_0_0_0, tg_z_xxxxy_s_0_0_0, tg_z_xxxxz_s_0_0_0, tg_z_xxxyy_s_0_0_0, tg_z_xxxyz_s_0_0_0, tg_z_xxxzz_s_0_0_0, tg_z_xxyyy_s_0_0_0, tg_z_xxyyz_s_0_0_0, tg_z_xxyzz_s_0_0_0, tg_z_xxzzz_s_0_0_0, tg_z_xyyyy_s_0_0_0, tg_z_xyyyz_s_0_0_0, tg_z_xyyzz_s_0_0_0, tg_z_xyzzz_s_0_0_0, tg_z_xzzzz_s_0_0_0, tg_z_yyyyy_s_0_0_0, tg_z_yyyyz_s_0_0_0, tg_z_yyyzz_s_0_0_0, tg_z_yyzzz_s_0_0_0, tg_z_yzzzz_s_0_0_0, tg_z_zzzzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xxxxx_s_0_0_0[i] += tg_0_xxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxy_s_0_0_0[i] += tg_0_xxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxz_s_0_0_0[i] += tg_0_xxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyy_s_0_0_0[i] += tg_0_xxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyz_s_0_0_0[i] += tg_0_xxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxzz_s_0_0_0[i] += tg_0_xxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyy_s_0_0_0[i] += tg_0_xxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyz_s_0_0_0[i] += tg_0_xxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyzz_s_0_0_0[i] += tg_0_xxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxzzz_s_0_0_0[i] += tg_0_xxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyy_s_0_0_0[i] += tg_0_xyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyz_s_0_0_0[i] += tg_0_xyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyzz_s_0_0_0[i] += tg_0_xyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyzzz_s_0_0_0[i] += tg_0_xyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xzzzz_s_0_0_0[i] += tg_0_xzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyy_s_0_0_0[i] += tg_0_yyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyz_s_0_0_0[i] += tg_0_yyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyzz_s_0_0_0[i] += tg_0_yyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyzzz_s_0_0_0[i] += tg_0_yyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yzzzz_s_0_0_0[i] += tg_0_yzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zzzzz_s_0_0_0[i] += tg_0_zzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xxxxx_s_0_0_0[i] += tg_0_xxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxy_s_0_0_0[i] += tg_0_xxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxz_s_0_0_0[i] += tg_0_xxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyy_s_0_0_0[i] += tg_0_xxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyz_s_0_0_0[i] += tg_0_xxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxzz_s_0_0_0[i] += tg_0_xxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyy_s_0_0_0[i] += tg_0_xxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyz_s_0_0_0[i] += tg_0_xxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyzz_s_0_0_0[i] += tg_0_xxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxzzz_s_0_0_0[i] += tg_0_xxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyy_s_0_0_0[i] += tg_0_xyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyz_s_0_0_0[i] += tg_0_xyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyzz_s_0_0_0[i] += tg_0_xyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyzzz_s_0_0_0[i] += tg_0_xyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xzzzz_s_0_0_0[i] += tg_0_xzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyy_s_0_0_0[i] += tg_0_yyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyz_s_0_0_0[i] += tg_0_yyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyzz_s_0_0_0[i] += tg_0_yyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyzzz_s_0_0_0[i] += tg_0_yyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yzzzz_s_0_0_0[i] += tg_0_yzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zzzzz_s_0_0_0[i] += tg_0_zzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xxxxx_s_0_0_0[i] += tg_0_xxxxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxy_s_0_0_0[i] += tg_0_xxxxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxz_s_0_0_0[i] += tg_0_xxxxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyy_s_0_0_0[i] += tg_0_xxxyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyz_s_0_0_0[i] += tg_0_xxxyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxzz_s_0_0_0[i] += tg_0_xxxzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyy_s_0_0_0[i] += tg_0_xxyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyz_s_0_0_0[i] += tg_0_xxyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyzz_s_0_0_0[i] += tg_0_xxyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxzzz_s_0_0_0[i] += tg_0_xxzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyy_s_0_0_0[i] += tg_0_xyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyz_s_0_0_0[i] += tg_0_xyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyzz_s_0_0_0[i] += tg_0_xyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyzzz_s_0_0_0[i] += tg_0_xyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xzzzz_s_0_0_0[i] += tg_0_xzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyy_s_0_0_0[i] += tg_0_yyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyz_s_0_0_0[i] += tg_0_yyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyzz_s_0_0_0[i] += tg_0_yyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyzzz_s_0_0_0[i] += tg_0_yyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yzzzz_s_0_0_0[i] += tg_0_yzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zzzzz_s_0_0_0[i] += tg_0_zzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

