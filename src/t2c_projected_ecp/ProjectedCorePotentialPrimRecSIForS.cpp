#include "ProjectedCorePotentialPrimRecSIForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_si_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_si_s_0_0_0,
                                        const size_t idx_sg_s_0_0_0,
                                        const size_t idx_sh_s_0_0_0,
                                        const size_t idx_sg_s_0_1_0,
                                        const size_t idx_sh_s_0_1_0,
                                        const int m,
                                        const size_t idx_sg_s_0_0_1,
                                        const size_t idx_sh_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
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

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0);

    auto tg_0_xxxy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 1);

    auto tg_0_xxxz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 2);

    auto tg_0_xxyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 3);

    auto tg_0_xxyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 4);

    auto tg_0_xxzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 5);

    auto tg_0_xyyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 6);

    auto tg_0_xyyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 7);

    auto tg_0_xyzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 8);

    auto tg_0_xzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 9);

    auto tg_0_yyyy_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 10);

    auto tg_0_yyyz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 11);

    auto tg_0_yyzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 12);

    auto tg_0_yzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 13);

    auto tg_0_zzzz_s_0_0_0 = pbuffer.data(idx_sg_s_0_0_0 + 14);

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

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0);

    auto tg_0_xxxy_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 1);

    auto tg_0_xxxz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 2);

    auto tg_0_xxyy_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 3);

    auto tg_0_xxyz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 4);

    auto tg_0_xxzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 5);

    auto tg_0_xyyy_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 6);

    auto tg_0_xyyz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 7);

    auto tg_0_xyzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 8);

    auto tg_0_xzzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 9);

    auto tg_0_yyyy_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 10);

    auto tg_0_yyyz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 11);

    auto tg_0_yyzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 12);

    auto tg_0_yzzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 13);

    auto tg_0_zzzz_s_0_1_0 = pbuffer.data(idx_sg_s_0_1_0 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0);

    auto tg_0_xxxxy_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 1);

    auto tg_0_xxxxz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 2);

    auto tg_0_xxxyy_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 3);

    auto tg_0_xxxyz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 4);

    auto tg_0_xxxzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 5);

    auto tg_0_xxyyy_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 6);

    auto tg_0_xxyyz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 7);

    auto tg_0_xxyzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 8);

    auto tg_0_xxzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 9);

    auto tg_0_xyyyy_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 10);

    auto tg_0_xyyyz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 11);

    auto tg_0_xyyzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 12);

    auto tg_0_xyzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 13);

    auto tg_0_xzzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 14);

    auto tg_0_yyyyy_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 15);

    auto tg_0_yyyyz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 16);

    auto tg_0_yyyzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 17);

    auto tg_0_yyzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 18);

    auto tg_0_yzzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 19);

    auto tg_0_zzzzz_s_0_1_0 = pbuffer.data(idx_sh_s_0_1_0 + 20);

    // Set up components of targeted buffer : SI

    auto tg_0_xxxxxx_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0);

    auto tg_0_xxxxxy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 1);

    auto tg_0_xxxxxz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 2);

    auto tg_0_xxxxyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 3);

    auto tg_0_xxxxyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 4);

    auto tg_0_xxxxzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 5);

    auto tg_0_xxxyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 6);

    auto tg_0_xxxyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 7);

    auto tg_0_xxxyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 8);

    auto tg_0_xxxzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 9);

    auto tg_0_xxyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 10);

    auto tg_0_xxyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 11);

    auto tg_0_xxyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 12);

    auto tg_0_xxyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 13);

    auto tg_0_xxzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 14);

    auto tg_0_xyyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 15);

    auto tg_0_xyyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 16);

    auto tg_0_xyyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 17);

    auto tg_0_xyyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 18);

    auto tg_0_xyzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 19);

    auto tg_0_xzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 20);

    auto tg_0_yyyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 21);

    auto tg_0_yyyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 22);

    auto tg_0_yyyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 23);

    auto tg_0_yyyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 24);

    auto tg_0_yyzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 25);

    auto tg_0_yzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 26);

    auto tg_0_zzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 27);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_s_0_0_0, tg_0_xxxx_s_0_1_0, tg_0_xxxxx_s_0_0_0, tg_0_xxxxx_s_0_1_0, tg_0_xxxxxx_s_0_0_0, tg_0_xxxxxy_s_0_0_0, tg_0_xxxxxz_s_0_0_0, tg_0_xxxxyy_s_0_0_0, tg_0_xxxxyz_s_0_0_0, tg_0_xxxxz_s_0_0_0, tg_0_xxxxz_s_0_1_0, tg_0_xxxxzz_s_0_0_0, tg_0_xxxyy_s_0_0_0, tg_0_xxxyy_s_0_1_0, tg_0_xxxyyy_s_0_0_0, tg_0_xxxyyz_s_0_0_0, tg_0_xxxyzz_s_0_0_0, tg_0_xxxzz_s_0_0_0, tg_0_xxxzz_s_0_1_0, tg_0_xxxzzz_s_0_0_0, tg_0_xxyy_s_0_0_0, tg_0_xxyy_s_0_1_0, tg_0_xxyyy_s_0_0_0, tg_0_xxyyy_s_0_1_0, tg_0_xxyyyy_s_0_0_0, tg_0_xxyyyz_s_0_0_0, tg_0_xxyyzz_s_0_0_0, tg_0_xxyzzz_s_0_0_0, tg_0_xxzz_s_0_0_0, tg_0_xxzz_s_0_1_0, tg_0_xxzzz_s_0_0_0, tg_0_xxzzz_s_0_1_0, tg_0_xxzzzz_s_0_0_0, tg_0_xyyy_s_0_0_0, tg_0_xyyy_s_0_1_0, tg_0_xyyyy_s_0_0_0, tg_0_xyyyy_s_0_1_0, tg_0_xyyyyy_s_0_0_0, tg_0_xyyyyz_s_0_0_0, tg_0_xyyyzz_s_0_0_0, tg_0_xyyzz_s_0_0_0, tg_0_xyyzz_s_0_1_0, tg_0_xyyzzz_s_0_0_0, tg_0_xyzzzz_s_0_0_0, tg_0_xzzz_s_0_0_0, tg_0_xzzz_s_0_1_0, tg_0_xzzzz_s_0_0_0, tg_0_xzzzz_s_0_1_0, tg_0_xzzzzz_s_0_0_0, tg_0_yyyy_s_0_0_0, tg_0_yyyy_s_0_1_0, tg_0_yyyyy_s_0_0_0, tg_0_yyyyy_s_0_1_0, tg_0_yyyyyy_s_0_0_0, tg_0_yyyyyz_s_0_0_0, tg_0_yyyyz_s_0_0_0, tg_0_yyyyz_s_0_1_0, tg_0_yyyyzz_s_0_0_0, tg_0_yyyzz_s_0_0_0, tg_0_yyyzz_s_0_1_0, tg_0_yyyzzz_s_0_0_0, tg_0_yyzz_s_0_0_0, tg_0_yyzz_s_0_1_0, tg_0_yyzzz_s_0_0_0, tg_0_yyzzz_s_0_1_0, tg_0_yyzzzz_s_0_0_0, tg_0_yzzz_s_0_0_0, tg_0_yzzz_s_0_1_0, tg_0_yzzzz_s_0_0_0, tg_0_yzzzz_s_0_1_0, tg_0_yzzzzz_s_0_0_0, tg_0_zzzz_s_0_0_0, tg_0_zzzz_s_0_1_0, tg_0_zzzzz_s_0_0_0, tg_0_zzzzz_s_0_1_0, tg_0_zzzzzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double fbz_0 = -(a_exp + c_exp) * fzi_0;

        const double fazi_0 = a_exp * fzi_0;

        const double fb_0 = b_exps[i];

        tg_0_xxxxxx_s_0_0_0[i] = 5.0 * tg_0_xxxx_s_0_0_0[i] * fzi_0 + 5.0 * tg_0_xxxx_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxxxx_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxxxy_s_0_0_0[i] = 2.0 * tg_0_xxxxx_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxxxz_s_0_0_0[i] = 2.0 * tg_0_xxxxx_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxxxyy_s_0_0_0[i] = 3.0 * tg_0_xxyy_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_xxyy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxxyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxxyz_s_0_0_0[i] = 2.0 * tg_0_xxxxz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxxzz_s_0_0_0[i] = 3.0 * tg_0_xxzz_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_xxzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxxzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxyyy_s_0_0_0[i] = 2.0 * tg_0_xyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_0_xyyy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxyyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxyyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxyyz_s_0_0_0[i] = 2.0 * tg_0_xxxyy_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxyy_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxxyzz_s_0_0_0[i] = 2.0 * tg_0_xxxzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxzzz_s_0_0_0[i] = 2.0 * tg_0_xzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_0_xzzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xxzzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxzzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyyyy_s_0_0_0[i] = tg_0_yyyy_s_0_0_0[i] * fzi_0 + tg_0_yyyy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xyyyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xyyyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyyyz_s_0_0_0[i] = 2.0 * tg_0_xxyyy_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxyyy_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxyyzz_s_0_0_0[i] = tg_0_yyzz_s_0_0_0[i] * fzi_0 + tg_0_yyzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xyyzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xyyzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyzzz_s_0_0_0[i] = 2.0 * tg_0_xxzzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxzzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxzzzz_s_0_0_0[i] = tg_0_zzzz_s_0_0_0[i] * fzi_0 + tg_0_zzzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_xzzzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xzzzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyyy_s_0_0_0[i] = 2.0 * tg_0_yyyyy_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyyz_s_0_0_0[i] = 2.0 * tg_0_yyyyz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyzz_s_0_0_0[i] = 2.0 * tg_0_yyyzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyzzz_s_0_0_0[i] = 2.0 * tg_0_yyzzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyzzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyzzzz_s_0_0_0[i] = 2.0 * tg_0_yzzzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzzzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xzzzzz_s_0_0_0[i] = 2.0 * tg_0_zzzzz_s_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_s_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_yyyyyy_s_0_0_0[i] = 5.0 * tg_0_yyyy_s_0_0_0[i] * fzi_0 + 5.0 * tg_0_yyyy_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yyyyy_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyyyyz_s_0_0_0[i] = 2.0 * tg_0_yyyyy_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_s_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_yyyyzz_s_0_0_0[i] = 3.0 * tg_0_yyzz_s_0_0_0[i] * fzi_0 + 3.0 * tg_0_yyzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yyyzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyyzzz_s_0_0_0[i] = 2.0 * tg_0_yzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_0_yzzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yyzzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyzzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyzzzz_s_0_0_0[i] = tg_0_zzzz_s_0_0_0[i] * fzi_0 + tg_0_zzzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_yzzzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzzzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yzzzzz_s_0_0_0[i] = 2.0 * tg_0_zzzzz_s_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_s_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_zzzzzz_s_0_0_0[i] = 5.0 * tg_0_zzzz_s_0_0_0[i] * fzi_0 + 5.0 * tg_0_zzzz_s_0_1_0[i] * fazi_0 * fazi_0 + 2.0 * tg_0_zzzzz_s_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_s_0_0_0[i] * rb_z[i] * fbz_0;
    }

    if (m > 0)
    {
        const double fm_0 = (double)m;

        // Set up components of auxiliary buffer : SG

        auto tg_0_xxxx_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1);

        auto tg_0_xxxy_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 1);

        auto tg_0_xxxz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 2);

        auto tg_0_xxyy_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 3);

        auto tg_0_xxyz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 4);

        auto tg_0_xxzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 5);

        auto tg_0_xyyy_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 6);

        auto tg_0_xyyz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 7);

        auto tg_0_xyzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 8);

        auto tg_0_xzzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 9);

        auto tg_0_yyyy_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 10);

        auto tg_0_yyyz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 11);

        auto tg_0_yyzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 12);

        auto tg_0_yzzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 13);

        auto tg_0_zzzz_s_0_0_1 = pbuffer.data(idx_sg_s_0_0_1 + 14);

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

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_s_0_0_1, tg_0_xxxxx_s_0_0_1, tg_0_xxxxxx_s_0_0_0, tg_0_xxxxxy_s_0_0_0, tg_0_xxxxxz_s_0_0_0, tg_0_xxxxyy_s_0_0_0, tg_0_xxxxyz_s_0_0_0, tg_0_xxxxz_s_0_0_1, tg_0_xxxxzz_s_0_0_0, tg_0_xxxyy_s_0_0_1, tg_0_xxxyyy_s_0_0_0, tg_0_xxxyyz_s_0_0_0, tg_0_xxxyzz_s_0_0_0, tg_0_xxxzz_s_0_0_1, tg_0_xxxzzz_s_0_0_0, tg_0_xxyy_s_0_0_1, tg_0_xxyyy_s_0_0_1, tg_0_xxyyyy_s_0_0_0, tg_0_xxyyyz_s_0_0_0, tg_0_xxyyzz_s_0_0_0, tg_0_xxyzzz_s_0_0_0, tg_0_xxzz_s_0_0_1, tg_0_xxzzz_s_0_0_1, tg_0_xxzzzz_s_0_0_0, tg_0_xyyy_s_0_0_1, tg_0_xyyyy_s_0_0_1, tg_0_xyyyyy_s_0_0_0, tg_0_xyyyyz_s_0_0_0, tg_0_xyyyzz_s_0_0_0, tg_0_xyyzz_s_0_0_1, tg_0_xyyzzz_s_0_0_0, tg_0_xyzzzz_s_0_0_0, tg_0_xzzz_s_0_0_1, tg_0_xzzzz_s_0_0_1, tg_0_xzzzzz_s_0_0_0, tg_0_yyyy_s_0_0_1, tg_0_yyyyy_s_0_0_1, tg_0_yyyyyy_s_0_0_0, tg_0_yyyyyz_s_0_0_0, tg_0_yyyyz_s_0_0_1, tg_0_yyyyzz_s_0_0_0, tg_0_yyyzz_s_0_0_1, tg_0_yyyzzz_s_0_0_0, tg_0_yyzz_s_0_0_1, tg_0_yyzzz_s_0_0_1, tg_0_yyzzzz_s_0_0_0, tg_0_yzzz_s_0_0_1, tg_0_yzzzz_s_0_0_1, tg_0_yzzzzz_s_0_0_0, tg_0_zzzz_s_0_0_1, tg_0_zzzzz_s_0_0_1, tg_0_zzzzzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fbi_0 = 1.0 / b_exps[i];

            tg_0_xxxxxx_s_0_0_0[i] += 5.0 / 2.0 * tg_0_xxxx_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxxx_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxxxy_s_0_0_0[i] += tg_0_xxxxx_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxxxz_s_0_0_0[i] += tg_0_xxxxx_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxxxyy_s_0_0_0[i] += 3.0 / 2.0 * tg_0_xxyy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxxyz_s_0_0_0[i] += tg_0_xxxxz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxxzz_s_0_0_0[i] += 3.0 / 2.0 * tg_0_xxzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxyyy_s_0_0_0[i] += tg_0_xyyy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxyyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxyyz_s_0_0_0[i] += tg_0_xxxyy_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxxyzz_s_0_0_0[i] += tg_0_xxxzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxzzz_s_0_0_0[i] += tg_0_xzzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxzzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xyyyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyyyz_s_0_0_0[i] += tg_0_xxyyy_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xyyzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyzzz_s_0_0_0[i] += tg_0_xxzzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xzzzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyyy_s_0_0_0[i] += tg_0_yyyyy_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyyz_s_0_0_0[i] += tg_0_yyyyz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyzz_s_0_0_0[i] += tg_0_yyyzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyzzz_s_0_0_0[i] += tg_0_yyzzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyzzzz_s_0_0_0[i] += tg_0_yzzzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xzzzzz_s_0_0_0[i] += tg_0_zzzzz_s_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_yyyyyy_s_0_0_0[i] += 5.0 / 2.0 * tg_0_yyyy_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyyyy_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyyyyz_s_0_0_0[i] += tg_0_yyyyy_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_yyyyzz_s_0_0_0[i] += 3.0 / 2.0 * tg_0_yyzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyyzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyyzzz_s_0_0_0[i] += tg_0_yzzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyzzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yzzzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yzzzzz_s_0_0_0[i] += tg_0_zzzzz_s_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_zzzzzz_s_0_0_0[i] += 5.0 / 2.0 * tg_0_zzzz_s_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_zzzzz_s_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;
        }
    }
}

} // t2pecp namespace

