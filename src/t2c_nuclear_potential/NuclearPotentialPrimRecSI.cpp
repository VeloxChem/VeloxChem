#include "NuclearPotentialPrimRecSI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_si(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_si,
                               const size_t              idx_npot_0_sg,
                               const size_t              idx_npot_1_sg,
                               const size_t              idx_npot_0_sh,
                               const size_t              idx_npot_1_sh,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpb,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_0 = pbuffer.data(idx_npot_0_sg);

    auto ta_0_xxyy_0 = pbuffer.data(idx_npot_0_sg + 3);

    auto ta_0_xxzz_0 = pbuffer.data(idx_npot_0_sg + 5);

    auto ta_0_xyyy_0 = pbuffer.data(idx_npot_0_sg + 6);

    auto ta_0_xzzz_0 = pbuffer.data(idx_npot_0_sg + 9);

    auto ta_0_yyyy_0 = pbuffer.data(idx_npot_0_sg + 10);

    auto ta_0_yyzz_0 = pbuffer.data(idx_npot_0_sg + 12);

    auto ta_0_yzzz_0 = pbuffer.data(idx_npot_0_sg + 13);

    auto ta_0_zzzz_0 = pbuffer.data(idx_npot_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_xyyy_1 = pbuffer.data(idx_npot_1_sg + 6);

    auto ta_0_xzzz_1 = pbuffer.data(idx_npot_1_sg + 9);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_yzzz_1 = pbuffer.data(idx_npot_1_sg + 13);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_0 = pbuffer.data(idx_npot_0_sh);

    auto ta_0_xxxxz_0 = pbuffer.data(idx_npot_0_sh + 2);

    auto ta_0_xxxyy_0 = pbuffer.data(idx_npot_0_sh + 3);

    auto ta_0_xxxzz_0 = pbuffer.data(idx_npot_0_sh + 5);

    auto ta_0_xxyyy_0 = pbuffer.data(idx_npot_0_sh + 6);

    auto ta_0_xxzzz_0 = pbuffer.data(idx_npot_0_sh + 9);

    auto ta_0_xyyyy_0 = pbuffer.data(idx_npot_0_sh + 10);

    auto ta_0_xyyzz_0 = pbuffer.data(idx_npot_0_sh + 12);

    auto ta_0_xzzzz_0 = pbuffer.data(idx_npot_0_sh + 14);

    auto ta_0_yyyyy_0 = pbuffer.data(idx_npot_0_sh + 15);

    auto ta_0_yyyyz_0 = pbuffer.data(idx_npot_0_sh + 16);

    auto ta_0_yyyzz_0 = pbuffer.data(idx_npot_0_sh + 17);

    auto ta_0_yyzzz_0 = pbuffer.data(idx_npot_0_sh + 18);

    auto ta_0_yzzzz_0 = pbuffer.data(idx_npot_0_sh + 19);

    auto ta_0_zzzzz_0 = pbuffer.data(idx_npot_0_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_1 = pbuffer.data(idx_npot_1_sh);

    auto ta_0_xxxxz_1 = pbuffer.data(idx_npot_1_sh + 2);

    auto ta_0_xxxyy_1 = pbuffer.data(idx_npot_1_sh + 3);

    auto ta_0_xxxzz_1 = pbuffer.data(idx_npot_1_sh + 5);

    auto ta_0_xxyyy_1 = pbuffer.data(idx_npot_1_sh + 6);

    auto ta_0_xxzzz_1 = pbuffer.data(idx_npot_1_sh + 9);

    auto ta_0_xyyyy_1 = pbuffer.data(idx_npot_1_sh + 10);

    auto ta_0_xyyzz_1 = pbuffer.data(idx_npot_1_sh + 12);

    auto ta_0_xzzzz_1 = pbuffer.data(idx_npot_1_sh + 14);

    auto ta_0_yyyyy_1 = pbuffer.data(idx_npot_1_sh + 15);

    auto ta_0_yyyyz_1 = pbuffer.data(idx_npot_1_sh + 16);

    auto ta_0_yyyzz_1 = pbuffer.data(idx_npot_1_sh + 17);

    auto ta_0_yyzzz_1 = pbuffer.data(idx_npot_1_sh + 18);

    auto ta_0_yzzzz_1 = pbuffer.data(idx_npot_1_sh + 19);

    auto ta_0_zzzzz_1 = pbuffer.data(idx_npot_1_sh + 20);

    // Set up components of targeted buffer : SI

    auto ta_0_xxxxxx_0 = pbuffer.data(idx_npot_0_si);

    auto ta_0_xxxxxy_0 = pbuffer.data(idx_npot_0_si + 1);

    auto ta_0_xxxxxz_0 = pbuffer.data(idx_npot_0_si + 2);

    auto ta_0_xxxxyy_0 = pbuffer.data(idx_npot_0_si + 3);

    auto ta_0_xxxxyz_0 = pbuffer.data(idx_npot_0_si + 4);

    auto ta_0_xxxxzz_0 = pbuffer.data(idx_npot_0_si + 5);

    auto ta_0_xxxyyy_0 = pbuffer.data(idx_npot_0_si + 6);

    auto ta_0_xxxyyz_0 = pbuffer.data(idx_npot_0_si + 7);

    auto ta_0_xxxyzz_0 = pbuffer.data(idx_npot_0_si + 8);

    auto ta_0_xxxzzz_0 = pbuffer.data(idx_npot_0_si + 9);

    auto ta_0_xxyyyy_0 = pbuffer.data(idx_npot_0_si + 10);

    auto ta_0_xxyyyz_0 = pbuffer.data(idx_npot_0_si + 11);

    auto ta_0_xxyyzz_0 = pbuffer.data(idx_npot_0_si + 12);

    auto ta_0_xxyzzz_0 = pbuffer.data(idx_npot_0_si + 13);

    auto ta_0_xxzzzz_0 = pbuffer.data(idx_npot_0_si + 14);

    auto ta_0_xyyyyy_0 = pbuffer.data(idx_npot_0_si + 15);

    auto ta_0_xyyyyz_0 = pbuffer.data(idx_npot_0_si + 16);

    auto ta_0_xyyyzz_0 = pbuffer.data(idx_npot_0_si + 17);

    auto ta_0_xyyzzz_0 = pbuffer.data(idx_npot_0_si + 18);

    auto ta_0_xyzzzz_0 = pbuffer.data(idx_npot_0_si + 19);

    auto ta_0_xzzzzz_0 = pbuffer.data(idx_npot_0_si + 20);

    auto ta_0_yyyyyy_0 = pbuffer.data(idx_npot_0_si + 21);

    auto ta_0_yyyyyz_0 = pbuffer.data(idx_npot_0_si + 22);

    auto ta_0_yyyyzz_0 = pbuffer.data(idx_npot_0_si + 23);

    auto ta_0_yyyzzz_0 = pbuffer.data(idx_npot_0_si + 24);

    auto ta_0_yyzzzz_0 = pbuffer.data(idx_npot_0_si + 25);

    auto ta_0_yzzzzz_0 = pbuffer.data(idx_npot_0_si + 26);

    auto ta_0_zzzzzz_0 = pbuffer.data(idx_npot_0_si + 27);

#pragma omp simd aligned(pb_x,              \
                             pb_y,          \
                             pb_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta_0_xxxx_0,   \
                             ta_0_xxxx_1,   \
                             ta_0_xxxxx_0,  \
                             ta_0_xxxxx_1,  \
                             ta_0_xxxxxx_0, \
                             ta_0_xxxxxy_0, \
                             ta_0_xxxxxz_0, \
                             ta_0_xxxxyy_0, \
                             ta_0_xxxxyz_0, \
                             ta_0_xxxxz_0,  \
                             ta_0_xxxxz_1,  \
                             ta_0_xxxxzz_0, \
                             ta_0_xxxyy_0,  \
                             ta_0_xxxyy_1,  \
                             ta_0_xxxyyy_0, \
                             ta_0_xxxyyz_0, \
                             ta_0_xxxyzz_0, \
                             ta_0_xxxzz_0,  \
                             ta_0_xxxzz_1,  \
                             ta_0_xxxzzz_0, \
                             ta_0_xxyy_0,   \
                             ta_0_xxyy_1,   \
                             ta_0_xxyyy_0,  \
                             ta_0_xxyyy_1,  \
                             ta_0_xxyyyy_0, \
                             ta_0_xxyyyz_0, \
                             ta_0_xxyyzz_0, \
                             ta_0_xxyzzz_0, \
                             ta_0_xxzz_0,   \
                             ta_0_xxzz_1,   \
                             ta_0_xxzzz_0,  \
                             ta_0_xxzzz_1,  \
                             ta_0_xxzzzz_0, \
                             ta_0_xyyy_0,   \
                             ta_0_xyyy_1,   \
                             ta_0_xyyyy_0,  \
                             ta_0_xyyyy_1,  \
                             ta_0_xyyyyy_0, \
                             ta_0_xyyyyz_0, \
                             ta_0_xyyyzz_0, \
                             ta_0_xyyzz_0,  \
                             ta_0_xyyzz_1,  \
                             ta_0_xyyzzz_0, \
                             ta_0_xyzzzz_0, \
                             ta_0_xzzz_0,   \
                             ta_0_xzzz_1,   \
                             ta_0_xzzzz_0,  \
                             ta_0_xzzzz_1,  \
                             ta_0_xzzzzz_0, \
                             ta_0_yyyy_0,   \
                             ta_0_yyyy_1,   \
                             ta_0_yyyyy_0,  \
                             ta_0_yyyyy_1,  \
                             ta_0_yyyyyy_0, \
                             ta_0_yyyyyz_0, \
                             ta_0_yyyyz_0,  \
                             ta_0_yyyyz_1,  \
                             ta_0_yyyyzz_0, \
                             ta_0_yyyzz_0,  \
                             ta_0_yyyzz_1,  \
                             ta_0_yyyzzz_0, \
                             ta_0_yyzz_0,   \
                             ta_0_yyzz_1,   \
                             ta_0_yyzzz_0,  \
                             ta_0_yyzzz_1,  \
                             ta_0_yyzzzz_0, \
                             ta_0_yzzz_0,   \
                             ta_0_yzzz_1,   \
                             ta_0_yzzzz_0,  \
                             ta_0_yzzzz_1,  \
                             ta_0_yzzzzz_0, \
                             ta_0_zzzz_0,   \
                             ta_0_zzzz_1,   \
                             ta_0_zzzzz_0,  \
                             ta_0_zzzzz_1,  \
                             ta_0_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_0_xxxxxx_0[i] = 5.0 * ta_0_xxxx_0[i] * fe_0 - 5.0 * ta_0_xxxx_1[i] * fe_0 + ta_0_xxxxx_0[i] * pb_x[i] - ta_0_xxxxx_1[i] * pc_x[i];

        ta_0_xxxxxy_0[i] = ta_0_xxxxx_0[i] * pb_y[i] - ta_0_xxxxx_1[i] * pc_y[i];

        ta_0_xxxxxz_0[i] = ta_0_xxxxx_0[i] * pb_z[i] - ta_0_xxxxx_1[i] * pc_z[i];

        ta_0_xxxxyy_0[i] = 3.0 * ta_0_xxyy_0[i] * fe_0 - 3.0 * ta_0_xxyy_1[i] * fe_0 + ta_0_xxxyy_0[i] * pb_x[i] - ta_0_xxxyy_1[i] * pc_x[i];

        ta_0_xxxxyz_0[i] = ta_0_xxxxz_0[i] * pb_y[i] - ta_0_xxxxz_1[i] * pc_y[i];

        ta_0_xxxxzz_0[i] = 3.0 * ta_0_xxzz_0[i] * fe_0 - 3.0 * ta_0_xxzz_1[i] * fe_0 + ta_0_xxxzz_0[i] * pb_x[i] - ta_0_xxxzz_1[i] * pc_x[i];

        ta_0_xxxyyy_0[i] = 2.0 * ta_0_xyyy_0[i] * fe_0 - 2.0 * ta_0_xyyy_1[i] * fe_0 + ta_0_xxyyy_0[i] * pb_x[i] - ta_0_xxyyy_1[i] * pc_x[i];

        ta_0_xxxyyz_0[i] = ta_0_xxxyy_0[i] * pb_z[i] - ta_0_xxxyy_1[i] * pc_z[i];

        ta_0_xxxyzz_0[i] = ta_0_xxxzz_0[i] * pb_y[i] - ta_0_xxxzz_1[i] * pc_y[i];

        ta_0_xxxzzz_0[i] = 2.0 * ta_0_xzzz_0[i] * fe_0 - 2.0 * ta_0_xzzz_1[i] * fe_0 + ta_0_xxzzz_0[i] * pb_x[i] - ta_0_xxzzz_1[i] * pc_x[i];

        ta_0_xxyyyy_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + ta_0_xyyyy_0[i] * pb_x[i] - ta_0_xyyyy_1[i] * pc_x[i];

        ta_0_xxyyyz_0[i] = ta_0_xxyyy_0[i] * pb_z[i] - ta_0_xxyyy_1[i] * pc_z[i];

        ta_0_xxyyzz_0[i] = ta_0_yyzz_0[i] * fe_0 - ta_0_yyzz_1[i] * fe_0 + ta_0_xyyzz_0[i] * pb_x[i] - ta_0_xyyzz_1[i] * pc_x[i];

        ta_0_xxyzzz_0[i] = ta_0_xxzzz_0[i] * pb_y[i] - ta_0_xxzzz_1[i] * pc_y[i];

        ta_0_xxzzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_0_xzzzz_0[i] * pb_x[i] - ta_0_xzzzz_1[i] * pc_x[i];

        ta_0_xyyyyy_0[i] = ta_0_yyyyy_0[i] * pb_x[i] - ta_0_yyyyy_1[i] * pc_x[i];

        ta_0_xyyyyz_0[i] = ta_0_yyyyz_0[i] * pb_x[i] - ta_0_yyyyz_1[i] * pc_x[i];

        ta_0_xyyyzz_0[i] = ta_0_yyyzz_0[i] * pb_x[i] - ta_0_yyyzz_1[i] * pc_x[i];

        ta_0_xyyzzz_0[i] = ta_0_yyzzz_0[i] * pb_x[i] - ta_0_yyzzz_1[i] * pc_x[i];

        ta_0_xyzzzz_0[i] = ta_0_yzzzz_0[i] * pb_x[i] - ta_0_yzzzz_1[i] * pc_x[i];

        ta_0_xzzzzz_0[i] = ta_0_zzzzz_0[i] * pb_x[i] - ta_0_zzzzz_1[i] * pc_x[i];

        ta_0_yyyyyy_0[i] = 5.0 * ta_0_yyyy_0[i] * fe_0 - 5.0 * ta_0_yyyy_1[i] * fe_0 + ta_0_yyyyy_0[i] * pb_y[i] - ta_0_yyyyy_1[i] * pc_y[i];

        ta_0_yyyyyz_0[i] = ta_0_yyyyy_0[i] * pb_z[i] - ta_0_yyyyy_1[i] * pc_z[i];

        ta_0_yyyyzz_0[i] = 3.0 * ta_0_yyzz_0[i] * fe_0 - 3.0 * ta_0_yyzz_1[i] * fe_0 + ta_0_yyyzz_0[i] * pb_y[i] - ta_0_yyyzz_1[i] * pc_y[i];

        ta_0_yyyzzz_0[i] = 2.0 * ta_0_yzzz_0[i] * fe_0 - 2.0 * ta_0_yzzz_1[i] * fe_0 + ta_0_yyzzz_0[i] * pb_y[i] - ta_0_yyzzz_1[i] * pc_y[i];

        ta_0_yyzzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_0_yzzzz_0[i] * pb_y[i] - ta_0_yzzzz_1[i] * pc_y[i];

        ta_0_yzzzzz_0[i] = ta_0_zzzzz_0[i] * pb_y[i] - ta_0_zzzzz_1[i] * pc_y[i];

        ta_0_zzzzzz_0[i] = 5.0 * ta_0_zzzz_0[i] * fe_0 - 5.0 * ta_0_zzzz_1[i] * fe_0 + ta_0_zzzzz_0[i] * pb_z[i] - ta_0_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
