#include "NuclearPotentialPrimRecPH.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_ph(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_ph,
                               const size_t idx_npot_0_sg,
                               const size_t idx_npot_1_sg,
                               const size_t idx_npot_0_sh,
                               const size_t idx_npot_1_sh,
                               const CSimdArray<double>& factors,
                               const size_t idx_rpa,
                               const size_t idx_rpc,
                               const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_0 = pbuffer.data(idx_npot_0_sg);

    auto ta_0_xxxy_0 = pbuffer.data(idx_npot_0_sg + 1);

    auto ta_0_xxxz_0 = pbuffer.data(idx_npot_0_sg + 2);

    auto ta_0_xxyy_0 = pbuffer.data(idx_npot_0_sg + 3);

    auto ta_0_xxyz_0 = pbuffer.data(idx_npot_0_sg + 4);

    auto ta_0_xxzz_0 = pbuffer.data(idx_npot_0_sg + 5);

    auto ta_0_xyyy_0 = pbuffer.data(idx_npot_0_sg + 6);

    auto ta_0_xyyz_0 = pbuffer.data(idx_npot_0_sg + 7);

    auto ta_0_xyzz_0 = pbuffer.data(idx_npot_0_sg + 8);

    auto ta_0_xzzz_0 = pbuffer.data(idx_npot_0_sg + 9);

    auto ta_0_yyyy_0 = pbuffer.data(idx_npot_0_sg + 10);

    auto ta_0_yyyz_0 = pbuffer.data(idx_npot_0_sg + 11);

    auto ta_0_yyzz_0 = pbuffer.data(idx_npot_0_sg + 12);

    auto ta_0_yzzz_0 = pbuffer.data(idx_npot_0_sg + 13);

    auto ta_0_zzzz_0 = pbuffer.data(idx_npot_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxxy_1 = pbuffer.data(idx_npot_1_sg + 1);

    auto ta_0_xxxz_1 = pbuffer.data(idx_npot_1_sg + 2);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxyz_1 = pbuffer.data(idx_npot_1_sg + 4);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_xyyy_1 = pbuffer.data(idx_npot_1_sg + 6);

    auto ta_0_xyyz_1 = pbuffer.data(idx_npot_1_sg + 7);

    auto ta_0_xyzz_1 = pbuffer.data(idx_npot_1_sg + 8);

    auto ta_0_xzzz_1 = pbuffer.data(idx_npot_1_sg + 9);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyyz_1 = pbuffer.data(idx_npot_1_sg + 11);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_yzzz_1 = pbuffer.data(idx_npot_1_sg + 13);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_0 = pbuffer.data(idx_npot_0_sh);

    auto ta_0_xxxxy_0 = pbuffer.data(idx_npot_0_sh + 1);

    auto ta_0_xxxxz_0 = pbuffer.data(idx_npot_0_sh + 2);

    auto ta_0_xxxyy_0 = pbuffer.data(idx_npot_0_sh + 3);

    auto ta_0_xxxyz_0 = pbuffer.data(idx_npot_0_sh + 4);

    auto ta_0_xxxzz_0 = pbuffer.data(idx_npot_0_sh + 5);

    auto ta_0_xxyyy_0 = pbuffer.data(idx_npot_0_sh + 6);

    auto ta_0_xxyyz_0 = pbuffer.data(idx_npot_0_sh + 7);

    auto ta_0_xxyzz_0 = pbuffer.data(idx_npot_0_sh + 8);

    auto ta_0_xxzzz_0 = pbuffer.data(idx_npot_0_sh + 9);

    auto ta_0_xyyyy_0 = pbuffer.data(idx_npot_0_sh + 10);

    auto ta_0_xyyyz_0 = pbuffer.data(idx_npot_0_sh + 11);

    auto ta_0_xyyzz_0 = pbuffer.data(idx_npot_0_sh + 12);

    auto ta_0_xyzzz_0 = pbuffer.data(idx_npot_0_sh + 13);

    auto ta_0_xzzzz_0 = pbuffer.data(idx_npot_0_sh + 14);

    auto ta_0_yyyyy_0 = pbuffer.data(idx_npot_0_sh + 15);

    auto ta_0_yyyyz_0 = pbuffer.data(idx_npot_0_sh + 16);

    auto ta_0_yyyzz_0 = pbuffer.data(idx_npot_0_sh + 17);

    auto ta_0_yyzzz_0 = pbuffer.data(idx_npot_0_sh + 18);

    auto ta_0_yzzzz_0 = pbuffer.data(idx_npot_0_sh + 19);

    auto ta_0_zzzzz_0 = pbuffer.data(idx_npot_0_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_1 = pbuffer.data(idx_npot_1_sh);

    auto ta_0_xxxxy_1 = pbuffer.data(idx_npot_1_sh + 1);

    auto ta_0_xxxxz_1 = pbuffer.data(idx_npot_1_sh + 2);

    auto ta_0_xxxyy_1 = pbuffer.data(idx_npot_1_sh + 3);

    auto ta_0_xxxyz_1 = pbuffer.data(idx_npot_1_sh + 4);

    auto ta_0_xxxzz_1 = pbuffer.data(idx_npot_1_sh + 5);

    auto ta_0_xxyyy_1 = pbuffer.data(idx_npot_1_sh + 6);

    auto ta_0_xxyyz_1 = pbuffer.data(idx_npot_1_sh + 7);

    auto ta_0_xxyzz_1 = pbuffer.data(idx_npot_1_sh + 8);

    auto ta_0_xxzzz_1 = pbuffer.data(idx_npot_1_sh + 9);

    auto ta_0_xyyyy_1 = pbuffer.data(idx_npot_1_sh + 10);

    auto ta_0_xyyyz_1 = pbuffer.data(idx_npot_1_sh + 11);

    auto ta_0_xyyzz_1 = pbuffer.data(idx_npot_1_sh + 12);

    auto ta_0_xyzzz_1 = pbuffer.data(idx_npot_1_sh + 13);

    auto ta_0_xzzzz_1 = pbuffer.data(idx_npot_1_sh + 14);

    auto ta_0_yyyyy_1 = pbuffer.data(idx_npot_1_sh + 15);

    auto ta_0_yyyyz_1 = pbuffer.data(idx_npot_1_sh + 16);

    auto ta_0_yyyzz_1 = pbuffer.data(idx_npot_1_sh + 17);

    auto ta_0_yyzzz_1 = pbuffer.data(idx_npot_1_sh + 18);

    auto ta_0_yzzzz_1 = pbuffer.data(idx_npot_1_sh + 19);

    auto ta_0_zzzzz_1 = pbuffer.data(idx_npot_1_sh + 20);

    // Set up 0-21 components of targeted buffer : PH

    auto ta_x_xxxxx_0 = pbuffer.data(idx_npot_0_ph);

    auto ta_x_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 1);

    auto ta_x_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 2);

    auto ta_x_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 3);

    auto ta_x_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 4);

    auto ta_x_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 5);

    auto ta_x_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 6);

    auto ta_x_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 7);

    auto ta_x_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 8);

    auto ta_x_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 9);

    auto ta_x_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 10);

    auto ta_x_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 11);

    auto ta_x_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 12);

    auto ta_x_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 13);

    auto ta_x_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 14);

    auto ta_x_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 15);

    auto ta_x_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 16);

    auto ta_x_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 17);

    auto ta_x_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 18);

    auto ta_x_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 19);

    auto ta_x_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 20);

    #pragma omp simd aligned(pa_x, pc_x, ta_0_xxxx_0, ta_0_xxxx_1, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxy_0, ta_0_xxxy_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyz_0, ta_0_xxyz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xyyy_0, ta_0_xyyy_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyz_0, ta_0_xyyz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyzz_0, ta_0_xyzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_yyyy_0, ta_0_yyyy_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yzzz_0, ta_0_yzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_zzzz_0, ta_0_zzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_x_xxxxx_0, ta_x_xxxxy_0, ta_x_xxxxz_0, ta_x_xxxyy_0, ta_x_xxxyz_0, ta_x_xxxzz_0, ta_x_xxyyy_0, ta_x_xxyyz_0, ta_x_xxyzz_0, ta_x_xxzzz_0, ta_x_xyyyy_0, ta_x_xyyyz_0, ta_x_xyyzz_0, ta_x_xyzzz_0, ta_x_xzzzz_0, ta_x_yyyyy_0, ta_x_yyyyz_0, ta_x_yyyzz_0, ta_x_yyzzz_0, ta_x_yzzzz_0, ta_x_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_xxxxx_0[i] = 5.0 * ta_0_xxxx_0[i] * fe_0 - 5.0 * ta_0_xxxx_1[i] * fe_0 + ta_0_xxxxx_0[i] * pa_x[i] - ta_0_xxxxx_1[i] * pc_x[i];

        ta_x_xxxxy_0[i] = 4.0 * ta_0_xxxy_0[i] * fe_0 - 4.0 * ta_0_xxxy_1[i] * fe_0 + ta_0_xxxxy_0[i] * pa_x[i] - ta_0_xxxxy_1[i] * pc_x[i];

        ta_x_xxxxz_0[i] = 4.0 * ta_0_xxxz_0[i] * fe_0 - 4.0 * ta_0_xxxz_1[i] * fe_0 + ta_0_xxxxz_0[i] * pa_x[i] - ta_0_xxxxz_1[i] * pc_x[i];

        ta_x_xxxyy_0[i] = 3.0 * ta_0_xxyy_0[i] * fe_0 - 3.0 * ta_0_xxyy_1[i] * fe_0 + ta_0_xxxyy_0[i] * pa_x[i] - ta_0_xxxyy_1[i] * pc_x[i];

        ta_x_xxxyz_0[i] = 3.0 * ta_0_xxyz_0[i] * fe_0 - 3.0 * ta_0_xxyz_1[i] * fe_0 + ta_0_xxxyz_0[i] * pa_x[i] - ta_0_xxxyz_1[i] * pc_x[i];

        ta_x_xxxzz_0[i] = 3.0 * ta_0_xxzz_0[i] * fe_0 - 3.0 * ta_0_xxzz_1[i] * fe_0 + ta_0_xxxzz_0[i] * pa_x[i] - ta_0_xxxzz_1[i] * pc_x[i];

        ta_x_xxyyy_0[i] = 2.0 * ta_0_xyyy_0[i] * fe_0 - 2.0 * ta_0_xyyy_1[i] * fe_0 + ta_0_xxyyy_0[i] * pa_x[i] - ta_0_xxyyy_1[i] * pc_x[i];

        ta_x_xxyyz_0[i] = 2.0 * ta_0_xyyz_0[i] * fe_0 - 2.0 * ta_0_xyyz_1[i] * fe_0 + ta_0_xxyyz_0[i] * pa_x[i] - ta_0_xxyyz_1[i] * pc_x[i];

        ta_x_xxyzz_0[i] = 2.0 * ta_0_xyzz_0[i] * fe_0 - 2.0 * ta_0_xyzz_1[i] * fe_0 + ta_0_xxyzz_0[i] * pa_x[i] - ta_0_xxyzz_1[i] * pc_x[i];

        ta_x_xxzzz_0[i] = 2.0 * ta_0_xzzz_0[i] * fe_0 - 2.0 * ta_0_xzzz_1[i] * fe_0 + ta_0_xxzzz_0[i] * pa_x[i] - ta_0_xxzzz_1[i] * pc_x[i];

        ta_x_xyyyy_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + ta_0_xyyyy_0[i] * pa_x[i] - ta_0_xyyyy_1[i] * pc_x[i];

        ta_x_xyyyz_0[i] = ta_0_yyyz_0[i] * fe_0 - ta_0_yyyz_1[i] * fe_0 + ta_0_xyyyz_0[i] * pa_x[i] - ta_0_xyyyz_1[i] * pc_x[i];

        ta_x_xyyzz_0[i] = ta_0_yyzz_0[i] * fe_0 - ta_0_yyzz_1[i] * fe_0 + ta_0_xyyzz_0[i] * pa_x[i] - ta_0_xyyzz_1[i] * pc_x[i];

        ta_x_xyzzz_0[i] = ta_0_yzzz_0[i] * fe_0 - ta_0_yzzz_1[i] * fe_0 + ta_0_xyzzz_0[i] * pa_x[i] - ta_0_xyzzz_1[i] * pc_x[i];

        ta_x_xzzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_0_xzzzz_0[i] * pa_x[i] - ta_0_xzzzz_1[i] * pc_x[i];

        ta_x_yyyyy_0[i] = ta_0_yyyyy_0[i] * pa_x[i] - ta_0_yyyyy_1[i] * pc_x[i];

        ta_x_yyyyz_0[i] = ta_0_yyyyz_0[i] * pa_x[i] - ta_0_yyyyz_1[i] * pc_x[i];

        ta_x_yyyzz_0[i] = ta_0_yyyzz_0[i] * pa_x[i] - ta_0_yyyzz_1[i] * pc_x[i];

        ta_x_yyzzz_0[i] = ta_0_yyzzz_0[i] * pa_x[i] - ta_0_yyzzz_1[i] * pc_x[i];

        ta_x_yzzzz_0[i] = ta_0_yzzzz_0[i] * pa_x[i] - ta_0_yzzzz_1[i] * pc_x[i];

        ta_x_zzzzz_0[i] = ta_0_zzzzz_0[i] * pa_x[i] - ta_0_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : PH

    auto ta_y_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 21);

    auto ta_y_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 22);

    auto ta_y_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 23);

    auto ta_y_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 24);

    auto ta_y_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 25);

    auto ta_y_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 26);

    auto ta_y_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 27);

    auto ta_y_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 28);

    auto ta_y_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 29);

    auto ta_y_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 30);

    auto ta_y_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 31);

    auto ta_y_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 32);

    auto ta_y_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 33);

    auto ta_y_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 34);

    auto ta_y_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 35);

    auto ta_y_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 36);

    auto ta_y_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 37);

    auto ta_y_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 38);

    auto ta_y_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 39);

    auto ta_y_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 40);

    auto ta_y_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 41);

    #pragma omp simd aligned(pa_y, pc_y, ta_0_xxxx_0, ta_0_xxxx_1, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxy_0, ta_0_xxxy_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyz_0, ta_0_xxyz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xyyy_0, ta_0_xyyy_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyz_0, ta_0_xyyz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyzz_0, ta_0_xyzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_yyyy_0, ta_0_yyyy_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yzzz_0, ta_0_yzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_zzzz_0, ta_0_zzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_y_xxxxx_0, ta_y_xxxxy_0, ta_y_xxxxz_0, ta_y_xxxyy_0, ta_y_xxxyz_0, ta_y_xxxzz_0, ta_y_xxyyy_0, ta_y_xxyyz_0, ta_y_xxyzz_0, ta_y_xxzzz_0, ta_y_xyyyy_0, ta_y_xyyyz_0, ta_y_xyyzz_0, ta_y_xyzzz_0, ta_y_xzzzz_0, ta_y_yyyyy_0, ta_y_yyyyz_0, ta_y_yyyzz_0, ta_y_yyzzz_0, ta_y_yzzzz_0, ta_y_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_xxxxx_0[i] = ta_0_xxxxx_0[i] * pa_y[i] - ta_0_xxxxx_1[i] * pc_y[i];

        ta_y_xxxxy_0[i] = ta_0_xxxx_0[i] * fe_0 - ta_0_xxxx_1[i] * fe_0 + ta_0_xxxxy_0[i] * pa_y[i] - ta_0_xxxxy_1[i] * pc_y[i];

        ta_y_xxxxz_0[i] = ta_0_xxxxz_0[i] * pa_y[i] - ta_0_xxxxz_1[i] * pc_y[i];

        ta_y_xxxyy_0[i] = 2.0 * ta_0_xxxy_0[i] * fe_0 - 2.0 * ta_0_xxxy_1[i] * fe_0 + ta_0_xxxyy_0[i] * pa_y[i] - ta_0_xxxyy_1[i] * pc_y[i];

        ta_y_xxxyz_0[i] = ta_0_xxxz_0[i] * fe_0 - ta_0_xxxz_1[i] * fe_0 + ta_0_xxxyz_0[i] * pa_y[i] - ta_0_xxxyz_1[i] * pc_y[i];

        ta_y_xxxzz_0[i] = ta_0_xxxzz_0[i] * pa_y[i] - ta_0_xxxzz_1[i] * pc_y[i];

        ta_y_xxyyy_0[i] = 3.0 * ta_0_xxyy_0[i] * fe_0 - 3.0 * ta_0_xxyy_1[i] * fe_0 + ta_0_xxyyy_0[i] * pa_y[i] - ta_0_xxyyy_1[i] * pc_y[i];

        ta_y_xxyyz_0[i] = 2.0 * ta_0_xxyz_0[i] * fe_0 - 2.0 * ta_0_xxyz_1[i] * fe_0 + ta_0_xxyyz_0[i] * pa_y[i] - ta_0_xxyyz_1[i] * pc_y[i];

        ta_y_xxyzz_0[i] = ta_0_xxzz_0[i] * fe_0 - ta_0_xxzz_1[i] * fe_0 + ta_0_xxyzz_0[i] * pa_y[i] - ta_0_xxyzz_1[i] * pc_y[i];

        ta_y_xxzzz_0[i] = ta_0_xxzzz_0[i] * pa_y[i] - ta_0_xxzzz_1[i] * pc_y[i];

        ta_y_xyyyy_0[i] = 4.0 * ta_0_xyyy_0[i] * fe_0 - 4.0 * ta_0_xyyy_1[i] * fe_0 + ta_0_xyyyy_0[i] * pa_y[i] - ta_0_xyyyy_1[i] * pc_y[i];

        ta_y_xyyyz_0[i] = 3.0 * ta_0_xyyz_0[i] * fe_0 - 3.0 * ta_0_xyyz_1[i] * fe_0 + ta_0_xyyyz_0[i] * pa_y[i] - ta_0_xyyyz_1[i] * pc_y[i];

        ta_y_xyyzz_0[i] = 2.0 * ta_0_xyzz_0[i] * fe_0 - 2.0 * ta_0_xyzz_1[i] * fe_0 + ta_0_xyyzz_0[i] * pa_y[i] - ta_0_xyyzz_1[i] * pc_y[i];

        ta_y_xyzzz_0[i] = ta_0_xzzz_0[i] * fe_0 - ta_0_xzzz_1[i] * fe_0 + ta_0_xyzzz_0[i] * pa_y[i] - ta_0_xyzzz_1[i] * pc_y[i];

        ta_y_xzzzz_0[i] = ta_0_xzzzz_0[i] * pa_y[i] - ta_0_xzzzz_1[i] * pc_y[i];

        ta_y_yyyyy_0[i] = 5.0 * ta_0_yyyy_0[i] * fe_0 - 5.0 * ta_0_yyyy_1[i] * fe_0 + ta_0_yyyyy_0[i] * pa_y[i] - ta_0_yyyyy_1[i] * pc_y[i];

        ta_y_yyyyz_0[i] = 4.0 * ta_0_yyyz_0[i] * fe_0 - 4.0 * ta_0_yyyz_1[i] * fe_0 + ta_0_yyyyz_0[i] * pa_y[i] - ta_0_yyyyz_1[i] * pc_y[i];

        ta_y_yyyzz_0[i] = 3.0 * ta_0_yyzz_0[i] * fe_0 - 3.0 * ta_0_yyzz_1[i] * fe_0 + ta_0_yyyzz_0[i] * pa_y[i] - ta_0_yyyzz_1[i] * pc_y[i];

        ta_y_yyzzz_0[i] = 2.0 * ta_0_yzzz_0[i] * fe_0 - 2.0 * ta_0_yzzz_1[i] * fe_0 + ta_0_yyzzz_0[i] * pa_y[i] - ta_0_yyzzz_1[i] * pc_y[i];

        ta_y_yzzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_0_yzzzz_0[i] * pa_y[i] - ta_0_yzzzz_1[i] * pc_y[i];

        ta_y_zzzzz_0[i] = ta_0_zzzzz_0[i] * pa_y[i] - ta_0_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : PH

    auto ta_z_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 42);

    auto ta_z_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 43);

    auto ta_z_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 44);

    auto ta_z_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 45);

    auto ta_z_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 46);

    auto ta_z_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 47);

    auto ta_z_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 48);

    auto ta_z_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 49);

    auto ta_z_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 50);

    auto ta_z_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 51);

    auto ta_z_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 52);

    auto ta_z_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 53);

    auto ta_z_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 54);

    auto ta_z_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 55);

    auto ta_z_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 56);

    auto ta_z_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 57);

    auto ta_z_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 58);

    auto ta_z_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 59);

    auto ta_z_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 60);

    auto ta_z_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 61);

    auto ta_z_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 62);

    #pragma omp simd aligned(pa_z, pc_z, ta_0_xxxx_0, ta_0_xxxx_1, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxy_0, ta_0_xxxy_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyz_0, ta_0_xxyz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xyyy_0, ta_0_xyyy_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyz_0, ta_0_xyyz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyzz_0, ta_0_xyzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_yyyy_0, ta_0_yyyy_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yzzz_0, ta_0_yzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_zzzz_0, ta_0_zzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_z_xxxxx_0, ta_z_xxxxy_0, ta_z_xxxxz_0, ta_z_xxxyy_0, ta_z_xxxyz_0, ta_z_xxxzz_0, ta_z_xxyyy_0, ta_z_xxyyz_0, ta_z_xxyzz_0, ta_z_xxzzz_0, ta_z_xyyyy_0, ta_z_xyyyz_0, ta_z_xyyzz_0, ta_z_xyzzz_0, ta_z_xzzzz_0, ta_z_yyyyy_0, ta_z_yyyyz_0, ta_z_yyyzz_0, ta_z_yyzzz_0, ta_z_yzzzz_0, ta_z_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_xxxxx_0[i] = ta_0_xxxxx_0[i] * pa_z[i] - ta_0_xxxxx_1[i] * pc_z[i];

        ta_z_xxxxy_0[i] = ta_0_xxxxy_0[i] * pa_z[i] - ta_0_xxxxy_1[i] * pc_z[i];

        ta_z_xxxxz_0[i] = ta_0_xxxx_0[i] * fe_0 - ta_0_xxxx_1[i] * fe_0 + ta_0_xxxxz_0[i] * pa_z[i] - ta_0_xxxxz_1[i] * pc_z[i];

        ta_z_xxxyy_0[i] = ta_0_xxxyy_0[i] * pa_z[i] - ta_0_xxxyy_1[i] * pc_z[i];

        ta_z_xxxyz_0[i] = ta_0_xxxy_0[i] * fe_0 - ta_0_xxxy_1[i] * fe_0 + ta_0_xxxyz_0[i] * pa_z[i] - ta_0_xxxyz_1[i] * pc_z[i];

        ta_z_xxxzz_0[i] = 2.0 * ta_0_xxxz_0[i] * fe_0 - 2.0 * ta_0_xxxz_1[i] * fe_0 + ta_0_xxxzz_0[i] * pa_z[i] - ta_0_xxxzz_1[i] * pc_z[i];

        ta_z_xxyyy_0[i] = ta_0_xxyyy_0[i] * pa_z[i] - ta_0_xxyyy_1[i] * pc_z[i];

        ta_z_xxyyz_0[i] = ta_0_xxyy_0[i] * fe_0 - ta_0_xxyy_1[i] * fe_0 + ta_0_xxyyz_0[i] * pa_z[i] - ta_0_xxyyz_1[i] * pc_z[i];

        ta_z_xxyzz_0[i] = 2.0 * ta_0_xxyz_0[i] * fe_0 - 2.0 * ta_0_xxyz_1[i] * fe_0 + ta_0_xxyzz_0[i] * pa_z[i] - ta_0_xxyzz_1[i] * pc_z[i];

        ta_z_xxzzz_0[i] = 3.0 * ta_0_xxzz_0[i] * fe_0 - 3.0 * ta_0_xxzz_1[i] * fe_0 + ta_0_xxzzz_0[i] * pa_z[i] - ta_0_xxzzz_1[i] * pc_z[i];

        ta_z_xyyyy_0[i] = ta_0_xyyyy_0[i] * pa_z[i] - ta_0_xyyyy_1[i] * pc_z[i];

        ta_z_xyyyz_0[i] = ta_0_xyyy_0[i] * fe_0 - ta_0_xyyy_1[i] * fe_0 + ta_0_xyyyz_0[i] * pa_z[i] - ta_0_xyyyz_1[i] * pc_z[i];

        ta_z_xyyzz_0[i] = 2.0 * ta_0_xyyz_0[i] * fe_0 - 2.0 * ta_0_xyyz_1[i] * fe_0 + ta_0_xyyzz_0[i] * pa_z[i] - ta_0_xyyzz_1[i] * pc_z[i];

        ta_z_xyzzz_0[i] = 3.0 * ta_0_xyzz_0[i] * fe_0 - 3.0 * ta_0_xyzz_1[i] * fe_0 + ta_0_xyzzz_0[i] * pa_z[i] - ta_0_xyzzz_1[i] * pc_z[i];

        ta_z_xzzzz_0[i] = 4.0 * ta_0_xzzz_0[i] * fe_0 - 4.0 * ta_0_xzzz_1[i] * fe_0 + ta_0_xzzzz_0[i] * pa_z[i] - ta_0_xzzzz_1[i] * pc_z[i];

        ta_z_yyyyy_0[i] = ta_0_yyyyy_0[i] * pa_z[i] - ta_0_yyyyy_1[i] * pc_z[i];

        ta_z_yyyyz_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + ta_0_yyyyz_0[i] * pa_z[i] - ta_0_yyyyz_1[i] * pc_z[i];

        ta_z_yyyzz_0[i] = 2.0 * ta_0_yyyz_0[i] * fe_0 - 2.0 * ta_0_yyyz_1[i] * fe_0 + ta_0_yyyzz_0[i] * pa_z[i] - ta_0_yyyzz_1[i] * pc_z[i];

        ta_z_yyzzz_0[i] = 3.0 * ta_0_yyzz_0[i] * fe_0 - 3.0 * ta_0_yyzz_1[i] * fe_0 + ta_0_yyzzz_0[i] * pa_z[i] - ta_0_yyzzz_1[i] * pc_z[i];

        ta_z_yzzzz_0[i] = 4.0 * ta_0_yzzz_0[i] * fe_0 - 4.0 * ta_0_yzzz_1[i] * fe_0 + ta_0_yzzzz_0[i] * pa_z[i] - ta_0_yzzzz_1[i] * pc_z[i];

        ta_z_zzzzz_0[i] = 5.0 * ta_0_zzzz_0[i] * fe_0 - 5.0 * ta_0_zzzz_1[i] * fe_0 + ta_0_zzzzz_0[i] * pa_z[i] - ta_0_zzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

