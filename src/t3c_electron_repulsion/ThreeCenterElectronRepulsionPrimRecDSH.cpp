#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsh,
                                 size_t idx_eri_0_ssh,
                                 size_t idx_eri_1_ssh,
                                 size_t idx_eri_1_psg,
                                 size_t idx_eri_1_psh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_ssh);

    auto g_0_0_xxxxy_0 = pbuffer.data(idx_eri_0_ssh + 1);

    auto g_0_0_xxxxz_0 = pbuffer.data(idx_eri_0_ssh + 2);

    auto g_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_ssh + 3);

    auto g_0_0_xxxyz_0 = pbuffer.data(idx_eri_0_ssh + 4);

    auto g_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_ssh + 5);

    auto g_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_ssh + 6);

    auto g_0_0_xxyyz_0 = pbuffer.data(idx_eri_0_ssh + 7);

    auto g_0_0_xxyzz_0 = pbuffer.data(idx_eri_0_ssh + 8);

    auto g_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_ssh + 9);

    auto g_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_ssh + 10);

    auto g_0_0_xyyyz_0 = pbuffer.data(idx_eri_0_ssh + 11);

    auto g_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_ssh + 12);

    auto g_0_0_xyzzz_0 = pbuffer.data(idx_eri_0_ssh + 13);

    auto g_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_ssh + 14);

    auto g_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_ssh + 15);

    auto g_0_0_yyyyz_0 = pbuffer.data(idx_eri_0_ssh + 16);

    auto g_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_ssh + 17);

    auto g_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_ssh + 18);

    auto g_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_ssh + 19);

    auto g_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_ssh + 20);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_ssh);

    auto g_0_0_xxxxy_1 = pbuffer.data(idx_eri_1_ssh + 1);

    auto g_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_ssh + 2);

    auto g_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_ssh + 3);

    auto g_0_0_xxxyz_1 = pbuffer.data(idx_eri_1_ssh + 4);

    auto g_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_ssh + 5);

    auto g_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_ssh + 6);

    auto g_0_0_xxyyz_1 = pbuffer.data(idx_eri_1_ssh + 7);

    auto g_0_0_xxyzz_1 = pbuffer.data(idx_eri_1_ssh + 8);

    auto g_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_ssh + 9);

    auto g_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_ssh + 10);

    auto g_0_0_xyyyz_1 = pbuffer.data(idx_eri_1_ssh + 11);

    auto g_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_ssh + 12);

    auto g_0_0_xyzzz_1 = pbuffer.data(idx_eri_1_ssh + 13);

    auto g_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_ssh + 14);

    auto g_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_ssh + 15);

    auto g_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_ssh + 16);

    auto g_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_ssh + 17);

    auto g_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_ssh + 18);

    auto g_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_ssh + 19);

    auto g_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_ssh + 20);

    /// Set up components of auxilary buffer : PSG

    auto g_x_0_xxxx_1 = pbuffer.data(idx_eri_1_psg);

    auto g_x_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 1);

    auto g_x_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 2);

    auto g_x_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 3);

    auto g_x_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 4);

    auto g_x_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 5);

    auto g_x_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 6);

    auto g_x_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 7);

    auto g_x_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 8);

    auto g_x_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 9);

    auto g_x_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 10);

    auto g_x_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 11);

    auto g_x_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 12);

    auto g_x_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 13);

    auto g_x_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 14);

    auto g_y_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 15);

    auto g_y_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 16);

    auto g_y_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 17);

    auto g_y_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 18);

    auto g_y_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 19);

    auto g_y_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 20);

    auto g_y_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 21);

    auto g_y_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 22);

    auto g_y_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 23);

    auto g_y_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 24);

    auto g_y_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 25);

    auto g_y_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 26);

    auto g_y_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 27);

    auto g_y_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 28);

    auto g_y_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 29);

    auto g_z_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 30);

    auto g_z_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 31);

    auto g_z_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 32);

    auto g_z_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 33);

    auto g_z_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 34);

    auto g_z_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 35);

    auto g_z_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 36);

    auto g_z_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 37);

    auto g_z_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 38);

    auto g_z_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 39);

    auto g_z_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 40);

    auto g_z_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 41);

    auto g_z_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 42);

    auto g_z_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 43);

    auto g_z_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 44);

    /// Set up components of auxilary buffer : PSH

    auto g_x_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh);

    auto g_x_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 1);

    auto g_x_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 2);

    auto g_x_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 3);

    auto g_x_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 4);

    auto g_x_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 5);

    auto g_x_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 6);

    auto g_x_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 7);

    auto g_x_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 8);

    auto g_x_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 9);

    auto g_x_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 10);

    auto g_x_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 11);

    auto g_x_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 12);

    auto g_x_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 13);

    auto g_x_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 14);

    auto g_x_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 15);

    auto g_x_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 16);

    auto g_x_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 17);

    auto g_x_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 18);

    auto g_x_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 19);

    auto g_x_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 20);

    auto g_y_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 21);

    auto g_y_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 22);

    auto g_y_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 23);

    auto g_y_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 24);

    auto g_y_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 25);

    auto g_y_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 26);

    auto g_y_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 27);

    auto g_y_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 28);

    auto g_y_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 29);

    auto g_y_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 30);

    auto g_y_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 31);

    auto g_y_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 32);

    auto g_y_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 33);

    auto g_y_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 34);

    auto g_y_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 35);

    auto g_y_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 36);

    auto g_y_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 37);

    auto g_y_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 38);

    auto g_y_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 39);

    auto g_y_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 40);

    auto g_y_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 41);

    auto g_z_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 42);

    auto g_z_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 43);

    auto g_z_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 44);

    auto g_z_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 45);

    auto g_z_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 46);

    auto g_z_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 47);

    auto g_z_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 48);

    auto g_z_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 49);

    auto g_z_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 50);

    auto g_z_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 51);

    auto g_z_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 52);

    auto g_z_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 53);

    auto g_z_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 54);

    auto g_z_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 55);

    auto g_z_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 56);

    auto g_z_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 57);

    auto g_z_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 58);

    auto g_z_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 59);

    auto g_z_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 60);

    auto g_z_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 61);

    auto g_z_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 62);

    /// Set up 0-21 components of targeted buffer : DSH

    auto g_xx_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh);

    auto g_xx_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 1);

    auto g_xx_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 2);

    auto g_xx_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 3);

    auto g_xx_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 4);

    auto g_xx_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 5);

    auto g_xx_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 6);

    auto g_xx_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 7);

    auto g_xx_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 8);

    auto g_xx_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 9);

    auto g_xx_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 10);

    auto g_xx_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 11);

    auto g_xx_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 12);

    auto g_xx_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 13);

    auto g_xx_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 14);

    auto g_xx_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 15);

    auto g_xx_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 16);

    auto g_xx_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 17);

    auto g_xx_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 18);

    auto g_xx_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 19);

    auto g_xx_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 20);

    #pragma omp simd aligned(g_0_0_xxxxx_0, g_0_0_xxxxx_1, g_0_0_xxxxy_0, g_0_0_xxxxy_1, g_0_0_xxxxz_0, g_0_0_xxxxz_1, g_0_0_xxxyy_0, g_0_0_xxxyy_1, g_0_0_xxxyz_0, g_0_0_xxxyz_1, g_0_0_xxxzz_0, g_0_0_xxxzz_1, g_0_0_xxyyy_0, g_0_0_xxyyy_1, g_0_0_xxyyz_0, g_0_0_xxyyz_1, g_0_0_xxyzz_0, g_0_0_xxyzz_1, g_0_0_xxzzz_0, g_0_0_xxzzz_1, g_0_0_xyyyy_0, g_0_0_xyyyy_1, g_0_0_xyyyz_0, g_0_0_xyyyz_1, g_0_0_xyyzz_0, g_0_0_xyyzz_1, g_0_0_xyzzz_0, g_0_0_xyzzz_1, g_0_0_xzzzz_0, g_0_0_xzzzz_1, g_0_0_yyyyy_0, g_0_0_yyyyy_1, g_0_0_yyyyz_0, g_0_0_yyyyz_1, g_0_0_yyyzz_0, g_0_0_yyyzz_1, g_0_0_yyzzz_0, g_0_0_yyzzz_1, g_0_0_yzzzz_0, g_0_0_yzzzz_1, g_0_0_zzzzz_0, g_0_0_zzzzz_1, g_x_0_xxxx_1, g_x_0_xxxxx_1, g_x_0_xxxxy_1, g_x_0_xxxxz_1, g_x_0_xxxy_1, g_x_0_xxxyy_1, g_x_0_xxxyz_1, g_x_0_xxxz_1, g_x_0_xxxzz_1, g_x_0_xxyy_1, g_x_0_xxyyy_1, g_x_0_xxyyz_1, g_x_0_xxyz_1, g_x_0_xxyzz_1, g_x_0_xxzz_1, g_x_0_xxzzz_1, g_x_0_xyyy_1, g_x_0_xyyyy_1, g_x_0_xyyyz_1, g_x_0_xyyz_1, g_x_0_xyyzz_1, g_x_0_xyzz_1, g_x_0_xyzzz_1, g_x_0_xzzz_1, g_x_0_xzzzz_1, g_x_0_yyyy_1, g_x_0_yyyyy_1, g_x_0_yyyyz_1, g_x_0_yyyz_1, g_x_0_yyyzz_1, g_x_0_yyzz_1, g_x_0_yyzzz_1, g_x_0_yzzz_1, g_x_0_yzzzz_1, g_x_0_zzzz_1, g_x_0_zzzzz_1, g_xx_0_xxxxx_0, g_xx_0_xxxxy_0, g_xx_0_xxxxz_0, g_xx_0_xxxyy_0, g_xx_0_xxxyz_0, g_xx_0_xxxzz_0, g_xx_0_xxyyy_0, g_xx_0_xxyyz_0, g_xx_0_xxyzz_0, g_xx_0_xxzzz_0, g_xx_0_xyyyy_0, g_xx_0_xyyyz_0, g_xx_0_xyyzz_0, g_xx_0_xyzzz_0, g_xx_0_xzzzz_0, g_xx_0_yyyyy_0, g_xx_0_yyyyz_0, g_xx_0_yyyzz_0, g_xx_0_yyzzz_0, g_xx_0_yzzzz_0, g_xx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxxx_0[i] = g_0_0_xxxxx_0[i] * fbe_0 - g_0_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_x_0_xxxx_1[i] * fi_acd_0 + g_x_0_xxxxx_1[i] * wa_x[i];

        g_xx_0_xxxxy_0[i] = g_0_0_xxxxy_0[i] * fbe_0 - g_0_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_x_0_xxxy_1[i] * fi_acd_0 + g_x_0_xxxxy_1[i] * wa_x[i];

        g_xx_0_xxxxz_0[i] = g_0_0_xxxxz_0[i] * fbe_0 - g_0_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxz_1[i] * fi_acd_0 + g_x_0_xxxxz_1[i] * wa_x[i];

        g_xx_0_xxxyy_0[i] = g_0_0_xxxyy_0[i] * fbe_0 - g_0_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_x_0_xxyy_1[i] * fi_acd_0 + g_x_0_xxxyy_1[i] * wa_x[i];

        g_xx_0_xxxyz_0[i] = g_0_0_xxxyz_0[i] * fbe_0 - g_0_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyz_1[i] * fi_acd_0 + g_x_0_xxxyz_1[i] * wa_x[i];

        g_xx_0_xxxzz_0[i] = g_0_0_xxxzz_0[i] * fbe_0 - g_0_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxzz_1[i] * fi_acd_0 + g_x_0_xxxzz_1[i] * wa_x[i];

        g_xx_0_xxyyy_0[i] = g_0_0_xxyyy_0[i] * fbe_0 - g_0_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyyy_1[i] * fi_acd_0 + g_x_0_xxyyy_1[i] * wa_x[i];

        g_xx_0_xxyyz_0[i] = g_0_0_xxyyz_0[i] * fbe_0 - g_0_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyz_1[i] * fi_acd_0 + g_x_0_xxyyz_1[i] * wa_x[i];

        g_xx_0_xxyzz_0[i] = g_0_0_xxyzz_0[i] * fbe_0 - g_0_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyzz_1[i] * fi_acd_0 + g_x_0_xxyzz_1[i] * wa_x[i];

        g_xx_0_xxzzz_0[i] = g_0_0_xxzzz_0[i] * fbe_0 - g_0_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzzz_1[i] * fi_acd_0 + g_x_0_xxzzz_1[i] * wa_x[i];

        g_xx_0_xyyyy_0[i] = g_0_0_xyyyy_0[i] * fbe_0 - g_0_0_xyyyy_1[i] * fz_be_0 + g_x_0_yyyy_1[i] * fi_acd_0 + g_x_0_xyyyy_1[i] * wa_x[i];

        g_xx_0_xyyyz_0[i] = g_0_0_xyyyz_0[i] * fbe_0 - g_0_0_xyyyz_1[i] * fz_be_0 + g_x_0_yyyz_1[i] * fi_acd_0 + g_x_0_xyyyz_1[i] * wa_x[i];

        g_xx_0_xyyzz_0[i] = g_0_0_xyyzz_0[i] * fbe_0 - g_0_0_xyyzz_1[i] * fz_be_0 + g_x_0_yyzz_1[i] * fi_acd_0 + g_x_0_xyyzz_1[i] * wa_x[i];

        g_xx_0_xyzzz_0[i] = g_0_0_xyzzz_0[i] * fbe_0 - g_0_0_xyzzz_1[i] * fz_be_0 + g_x_0_yzzz_1[i] * fi_acd_0 + g_x_0_xyzzz_1[i] * wa_x[i];

        g_xx_0_xzzzz_0[i] = g_0_0_xzzzz_0[i] * fbe_0 - g_0_0_xzzzz_1[i] * fz_be_0 + g_x_0_zzzz_1[i] * fi_acd_0 + g_x_0_xzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyy_0[i] = g_0_0_yyyyy_0[i] * fbe_0 - g_0_0_yyyyy_1[i] * fz_be_0 + g_x_0_yyyyy_1[i] * wa_x[i];

        g_xx_0_yyyyz_0[i] = g_0_0_yyyyz_0[i] * fbe_0 - g_0_0_yyyyz_1[i] * fz_be_0 + g_x_0_yyyyz_1[i] * wa_x[i];

        g_xx_0_yyyzz_0[i] = g_0_0_yyyzz_0[i] * fbe_0 - g_0_0_yyyzz_1[i] * fz_be_0 + g_x_0_yyyzz_1[i] * wa_x[i];

        g_xx_0_yyzzz_0[i] = g_0_0_yyzzz_0[i] * fbe_0 - g_0_0_yyzzz_1[i] * fz_be_0 + g_x_0_yyzzz_1[i] * wa_x[i];

        g_xx_0_yzzzz_0[i] = g_0_0_yzzzz_0[i] * fbe_0 - g_0_0_yzzzz_1[i] * fz_be_0 + g_x_0_yzzzz_1[i] * wa_x[i];

        g_xx_0_zzzzz_0[i] = g_0_0_zzzzz_0[i] * fbe_0 - g_0_0_zzzzz_1[i] * fz_be_0 + g_x_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : DSH

    auto g_xy_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 21);

    auto g_xy_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 22);

    auto g_xy_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 23);

    auto g_xy_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 24);

    auto g_xy_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 25);

    auto g_xy_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 26);

    auto g_xy_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 27);

    auto g_xy_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 28);

    auto g_xy_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 29);

    auto g_xy_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 30);

    auto g_xy_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 31);

    auto g_xy_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 32);

    auto g_xy_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 33);

    auto g_xy_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 34);

    auto g_xy_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 35);

    auto g_xy_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 36);

    auto g_xy_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 37);

    auto g_xy_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 38);

    auto g_xy_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 39);

    auto g_xy_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 40);

    auto g_xy_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 41);

    #pragma omp simd aligned(g_x_0_xxxxx_1, g_x_0_xxxxz_1, g_x_0_xxxzz_1, g_x_0_xxzzz_1, g_x_0_xzzzz_1, g_xy_0_xxxxx_0, g_xy_0_xxxxy_0, g_xy_0_xxxxz_0, g_xy_0_xxxyy_0, g_xy_0_xxxyz_0, g_xy_0_xxxzz_0, g_xy_0_xxyyy_0, g_xy_0_xxyyz_0, g_xy_0_xxyzz_0, g_xy_0_xxzzz_0, g_xy_0_xyyyy_0, g_xy_0_xyyyz_0, g_xy_0_xyyzz_0, g_xy_0_xyzzz_0, g_xy_0_xzzzz_0, g_xy_0_yyyyy_0, g_xy_0_yyyyz_0, g_xy_0_yyyzz_0, g_xy_0_yyzzz_0, g_xy_0_yzzzz_0, g_xy_0_zzzzz_0, g_y_0_xxxxy_1, g_y_0_xxxy_1, g_y_0_xxxyy_1, g_y_0_xxxyz_1, g_y_0_xxyy_1, g_y_0_xxyyy_1, g_y_0_xxyyz_1, g_y_0_xxyz_1, g_y_0_xxyzz_1, g_y_0_xyyy_1, g_y_0_xyyyy_1, g_y_0_xyyyz_1, g_y_0_xyyz_1, g_y_0_xyyzz_1, g_y_0_xyzz_1, g_y_0_xyzzz_1, g_y_0_yyyy_1, g_y_0_yyyyy_1, g_y_0_yyyyz_1, g_y_0_yyyz_1, g_y_0_yyyzz_1, g_y_0_yyzz_1, g_y_0_yyzzz_1, g_y_0_yzzz_1, g_y_0_yzzzz_1, g_y_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxxx_0[i] = g_x_0_xxxxx_1[i] * wa_y[i];

        g_xy_0_xxxxy_0[i] = 4.0 * g_y_0_xxxy_1[i] * fi_acd_0 + g_y_0_xxxxy_1[i] * wa_x[i];

        g_xy_0_xxxxz_0[i] = g_x_0_xxxxz_1[i] * wa_y[i];

        g_xy_0_xxxyy_0[i] = 3.0 * g_y_0_xxyy_1[i] * fi_acd_0 + g_y_0_xxxyy_1[i] * wa_x[i];

        g_xy_0_xxxyz_0[i] = 3.0 * g_y_0_xxyz_1[i] * fi_acd_0 + g_y_0_xxxyz_1[i] * wa_x[i];

        g_xy_0_xxxzz_0[i] = g_x_0_xxxzz_1[i] * wa_y[i];

        g_xy_0_xxyyy_0[i] = 2.0 * g_y_0_xyyy_1[i] * fi_acd_0 + g_y_0_xxyyy_1[i] * wa_x[i];

        g_xy_0_xxyyz_0[i] = 2.0 * g_y_0_xyyz_1[i] * fi_acd_0 + g_y_0_xxyyz_1[i] * wa_x[i];

        g_xy_0_xxyzz_0[i] = 2.0 * g_y_0_xyzz_1[i] * fi_acd_0 + g_y_0_xxyzz_1[i] * wa_x[i];

        g_xy_0_xxzzz_0[i] = g_x_0_xxzzz_1[i] * wa_y[i];

        g_xy_0_xyyyy_0[i] = g_y_0_yyyy_1[i] * fi_acd_0 + g_y_0_xyyyy_1[i] * wa_x[i];

        g_xy_0_xyyyz_0[i] = g_y_0_yyyz_1[i] * fi_acd_0 + g_y_0_xyyyz_1[i] * wa_x[i];

        g_xy_0_xyyzz_0[i] = g_y_0_yyzz_1[i] * fi_acd_0 + g_y_0_xyyzz_1[i] * wa_x[i];

        g_xy_0_xyzzz_0[i] = g_y_0_yzzz_1[i] * fi_acd_0 + g_y_0_xyzzz_1[i] * wa_x[i];

        g_xy_0_xzzzz_0[i] = g_x_0_xzzzz_1[i] * wa_y[i];

        g_xy_0_yyyyy_0[i] = g_y_0_yyyyy_1[i] * wa_x[i];

        g_xy_0_yyyyz_0[i] = g_y_0_yyyyz_1[i] * wa_x[i];

        g_xy_0_yyyzz_0[i] = g_y_0_yyyzz_1[i] * wa_x[i];

        g_xy_0_yyzzz_0[i] = g_y_0_yyzzz_1[i] * wa_x[i];

        g_xy_0_yzzzz_0[i] = g_y_0_yzzzz_1[i] * wa_x[i];

        g_xy_0_zzzzz_0[i] = g_y_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 42-63 components of targeted buffer : DSH

    auto g_xz_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 42);

    auto g_xz_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 43);

    auto g_xz_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 44);

    auto g_xz_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 45);

    auto g_xz_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 46);

    auto g_xz_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 47);

    auto g_xz_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 48);

    auto g_xz_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 49);

    auto g_xz_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 50);

    auto g_xz_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 51);

    auto g_xz_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 52);

    auto g_xz_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 53);

    auto g_xz_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 54);

    auto g_xz_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 55);

    auto g_xz_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 56);

    auto g_xz_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 57);

    auto g_xz_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 58);

    auto g_xz_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 59);

    auto g_xz_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 60);

    auto g_xz_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 61);

    auto g_xz_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 62);

    #pragma omp simd aligned(g_x_0_xxxxx_1, g_x_0_xxxxy_1, g_x_0_xxxyy_1, g_x_0_xxyyy_1, g_x_0_xyyyy_1, g_xz_0_xxxxx_0, g_xz_0_xxxxy_0, g_xz_0_xxxxz_0, g_xz_0_xxxyy_0, g_xz_0_xxxyz_0, g_xz_0_xxxzz_0, g_xz_0_xxyyy_0, g_xz_0_xxyyz_0, g_xz_0_xxyzz_0, g_xz_0_xxzzz_0, g_xz_0_xyyyy_0, g_xz_0_xyyyz_0, g_xz_0_xyyzz_0, g_xz_0_xyzzz_0, g_xz_0_xzzzz_0, g_xz_0_yyyyy_0, g_xz_0_yyyyz_0, g_xz_0_yyyzz_0, g_xz_0_yyzzz_0, g_xz_0_yzzzz_0, g_xz_0_zzzzz_0, g_z_0_xxxxz_1, g_z_0_xxxyz_1, g_z_0_xxxz_1, g_z_0_xxxzz_1, g_z_0_xxyyz_1, g_z_0_xxyz_1, g_z_0_xxyzz_1, g_z_0_xxzz_1, g_z_0_xxzzz_1, g_z_0_xyyyz_1, g_z_0_xyyz_1, g_z_0_xyyzz_1, g_z_0_xyzz_1, g_z_0_xyzzz_1, g_z_0_xzzz_1, g_z_0_xzzzz_1, g_z_0_yyyyy_1, g_z_0_yyyyz_1, g_z_0_yyyz_1, g_z_0_yyyzz_1, g_z_0_yyzz_1, g_z_0_yyzzz_1, g_z_0_yzzz_1, g_z_0_yzzzz_1, g_z_0_zzzz_1, g_z_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxxx_0[i] = g_x_0_xxxxx_1[i] * wa_z[i];

        g_xz_0_xxxxy_0[i] = g_x_0_xxxxy_1[i] * wa_z[i];

        g_xz_0_xxxxz_0[i] = 4.0 * g_z_0_xxxz_1[i] * fi_acd_0 + g_z_0_xxxxz_1[i] * wa_x[i];

        g_xz_0_xxxyy_0[i] = g_x_0_xxxyy_1[i] * wa_z[i];

        g_xz_0_xxxyz_0[i] = 3.0 * g_z_0_xxyz_1[i] * fi_acd_0 + g_z_0_xxxyz_1[i] * wa_x[i];

        g_xz_0_xxxzz_0[i] = 3.0 * g_z_0_xxzz_1[i] * fi_acd_0 + g_z_0_xxxzz_1[i] * wa_x[i];

        g_xz_0_xxyyy_0[i] = g_x_0_xxyyy_1[i] * wa_z[i];

        g_xz_0_xxyyz_0[i] = 2.0 * g_z_0_xyyz_1[i] * fi_acd_0 + g_z_0_xxyyz_1[i] * wa_x[i];

        g_xz_0_xxyzz_0[i] = 2.0 * g_z_0_xyzz_1[i] * fi_acd_0 + g_z_0_xxyzz_1[i] * wa_x[i];

        g_xz_0_xxzzz_0[i] = 2.0 * g_z_0_xzzz_1[i] * fi_acd_0 + g_z_0_xxzzz_1[i] * wa_x[i];

        g_xz_0_xyyyy_0[i] = g_x_0_xyyyy_1[i] * wa_z[i];

        g_xz_0_xyyyz_0[i] = g_z_0_yyyz_1[i] * fi_acd_0 + g_z_0_xyyyz_1[i] * wa_x[i];

        g_xz_0_xyyzz_0[i] = g_z_0_yyzz_1[i] * fi_acd_0 + g_z_0_xyyzz_1[i] * wa_x[i];

        g_xz_0_xyzzz_0[i] = g_z_0_yzzz_1[i] * fi_acd_0 + g_z_0_xyzzz_1[i] * wa_x[i];

        g_xz_0_xzzzz_0[i] = g_z_0_zzzz_1[i] * fi_acd_0 + g_z_0_xzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyy_0[i] = g_z_0_yyyyy_1[i] * wa_x[i];

        g_xz_0_yyyyz_0[i] = g_z_0_yyyyz_1[i] * wa_x[i];

        g_xz_0_yyyzz_0[i] = g_z_0_yyyzz_1[i] * wa_x[i];

        g_xz_0_yyzzz_0[i] = g_z_0_yyzzz_1[i] * wa_x[i];

        g_xz_0_yzzzz_0[i] = g_z_0_yzzzz_1[i] * wa_x[i];

        g_xz_0_zzzzz_0[i] = g_z_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 63-84 components of targeted buffer : DSH

    auto g_yy_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 63);

    auto g_yy_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 64);

    auto g_yy_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 65);

    auto g_yy_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 66);

    auto g_yy_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 67);

    auto g_yy_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 68);

    auto g_yy_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 69);

    auto g_yy_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 70);

    auto g_yy_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 71);

    auto g_yy_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 72);

    auto g_yy_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 73);

    auto g_yy_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 74);

    auto g_yy_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 75);

    auto g_yy_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 76);

    auto g_yy_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 77);

    auto g_yy_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 78);

    auto g_yy_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 79);

    auto g_yy_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 80);

    auto g_yy_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 81);

    auto g_yy_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 82);

    auto g_yy_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 83);

    #pragma omp simd aligned(g_0_0_xxxxx_0, g_0_0_xxxxx_1, g_0_0_xxxxy_0, g_0_0_xxxxy_1, g_0_0_xxxxz_0, g_0_0_xxxxz_1, g_0_0_xxxyy_0, g_0_0_xxxyy_1, g_0_0_xxxyz_0, g_0_0_xxxyz_1, g_0_0_xxxzz_0, g_0_0_xxxzz_1, g_0_0_xxyyy_0, g_0_0_xxyyy_1, g_0_0_xxyyz_0, g_0_0_xxyyz_1, g_0_0_xxyzz_0, g_0_0_xxyzz_1, g_0_0_xxzzz_0, g_0_0_xxzzz_1, g_0_0_xyyyy_0, g_0_0_xyyyy_1, g_0_0_xyyyz_0, g_0_0_xyyyz_1, g_0_0_xyyzz_0, g_0_0_xyyzz_1, g_0_0_xyzzz_0, g_0_0_xyzzz_1, g_0_0_xzzzz_0, g_0_0_xzzzz_1, g_0_0_yyyyy_0, g_0_0_yyyyy_1, g_0_0_yyyyz_0, g_0_0_yyyyz_1, g_0_0_yyyzz_0, g_0_0_yyyzz_1, g_0_0_yyzzz_0, g_0_0_yyzzz_1, g_0_0_yzzzz_0, g_0_0_yzzzz_1, g_0_0_zzzzz_0, g_0_0_zzzzz_1, g_y_0_xxxx_1, g_y_0_xxxxx_1, g_y_0_xxxxy_1, g_y_0_xxxxz_1, g_y_0_xxxy_1, g_y_0_xxxyy_1, g_y_0_xxxyz_1, g_y_0_xxxz_1, g_y_0_xxxzz_1, g_y_0_xxyy_1, g_y_0_xxyyy_1, g_y_0_xxyyz_1, g_y_0_xxyz_1, g_y_0_xxyzz_1, g_y_0_xxzz_1, g_y_0_xxzzz_1, g_y_0_xyyy_1, g_y_0_xyyyy_1, g_y_0_xyyyz_1, g_y_0_xyyz_1, g_y_0_xyyzz_1, g_y_0_xyzz_1, g_y_0_xyzzz_1, g_y_0_xzzz_1, g_y_0_xzzzz_1, g_y_0_yyyy_1, g_y_0_yyyyy_1, g_y_0_yyyyz_1, g_y_0_yyyz_1, g_y_0_yyyzz_1, g_y_0_yyzz_1, g_y_0_yyzzz_1, g_y_0_yzzz_1, g_y_0_yzzzz_1, g_y_0_zzzz_1, g_y_0_zzzzz_1, g_yy_0_xxxxx_0, g_yy_0_xxxxy_0, g_yy_0_xxxxz_0, g_yy_0_xxxyy_0, g_yy_0_xxxyz_0, g_yy_0_xxxzz_0, g_yy_0_xxyyy_0, g_yy_0_xxyyz_0, g_yy_0_xxyzz_0, g_yy_0_xxzzz_0, g_yy_0_xyyyy_0, g_yy_0_xyyyz_0, g_yy_0_xyyzz_0, g_yy_0_xyzzz_0, g_yy_0_xzzzz_0, g_yy_0_yyyyy_0, g_yy_0_yyyyz_0, g_yy_0_yyyzz_0, g_yy_0_yyzzz_0, g_yy_0_yzzzz_0, g_yy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxxx_0[i] = g_0_0_xxxxx_0[i] * fbe_0 - g_0_0_xxxxx_1[i] * fz_be_0 + g_y_0_xxxxx_1[i] * wa_y[i];

        g_yy_0_xxxxy_0[i] = g_0_0_xxxxy_0[i] * fbe_0 - g_0_0_xxxxy_1[i] * fz_be_0 + g_y_0_xxxx_1[i] * fi_acd_0 + g_y_0_xxxxy_1[i] * wa_y[i];

        g_yy_0_xxxxz_0[i] = g_0_0_xxxxz_0[i] * fbe_0 - g_0_0_xxxxz_1[i] * fz_be_0 + g_y_0_xxxxz_1[i] * wa_y[i];

        g_yy_0_xxxyy_0[i] = g_0_0_xxxyy_0[i] * fbe_0 - g_0_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxxy_1[i] * fi_acd_0 + g_y_0_xxxyy_1[i] * wa_y[i];

        g_yy_0_xxxyz_0[i] = g_0_0_xxxyz_0[i] * fbe_0 - g_0_0_xxxyz_1[i] * fz_be_0 + g_y_0_xxxz_1[i] * fi_acd_0 + g_y_0_xxxyz_1[i] * wa_y[i];

        g_yy_0_xxxzz_0[i] = g_0_0_xxxzz_0[i] * fbe_0 - g_0_0_xxxzz_1[i] * fz_be_0 + g_y_0_xxxzz_1[i] * wa_y[i];

        g_yy_0_xxyyy_0[i] = g_0_0_xxyyy_0[i] * fbe_0 - g_0_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xxyy_1[i] * fi_acd_0 + g_y_0_xxyyy_1[i] * wa_y[i];

        g_yy_0_xxyyz_0[i] = g_0_0_xxyyz_0[i] * fbe_0 - g_0_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xxyz_1[i] * fi_acd_0 + g_y_0_xxyyz_1[i] * wa_y[i];

        g_yy_0_xxyzz_0[i] = g_0_0_xxyzz_0[i] * fbe_0 - g_0_0_xxyzz_1[i] * fz_be_0 + g_y_0_xxzz_1[i] * fi_acd_0 + g_y_0_xxyzz_1[i] * wa_y[i];

        g_yy_0_xxzzz_0[i] = g_0_0_xxzzz_0[i] * fbe_0 - g_0_0_xxzzz_1[i] * fz_be_0 + g_y_0_xxzzz_1[i] * wa_y[i];

        g_yy_0_xyyyy_0[i] = g_0_0_xyyyy_0[i] * fbe_0 - g_0_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_y_0_xyyy_1[i] * fi_acd_0 + g_y_0_xyyyy_1[i] * wa_y[i];

        g_yy_0_xyyyz_0[i] = g_0_0_xyyyz_0[i] * fbe_0 - g_0_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_y_0_xyyz_1[i] * fi_acd_0 + g_y_0_xyyyz_1[i] * wa_y[i];

        g_yy_0_xyyzz_0[i] = g_0_0_xyyzz_0[i] * fbe_0 - g_0_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_y_0_xyzz_1[i] * fi_acd_0 + g_y_0_xyyzz_1[i] * wa_y[i];

        g_yy_0_xyzzz_0[i] = g_0_0_xyzzz_0[i] * fbe_0 - g_0_0_xyzzz_1[i] * fz_be_0 + g_y_0_xzzz_1[i] * fi_acd_0 + g_y_0_xyzzz_1[i] * wa_y[i];

        g_yy_0_xzzzz_0[i] = g_0_0_xzzzz_0[i] * fbe_0 - g_0_0_xzzzz_1[i] * fz_be_0 + g_y_0_xzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyy_0[i] = g_0_0_yyyyy_0[i] * fbe_0 - g_0_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_y_0_yyyy_1[i] * fi_acd_0 + g_y_0_yyyyy_1[i] * wa_y[i];

        g_yy_0_yyyyz_0[i] = g_0_0_yyyyz_0[i] * fbe_0 - g_0_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_y_0_yyyz_1[i] * fi_acd_0 + g_y_0_yyyyz_1[i] * wa_y[i];

        g_yy_0_yyyzz_0[i] = g_0_0_yyyzz_0[i] * fbe_0 - g_0_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_y_0_yyzz_1[i] * fi_acd_0 + g_y_0_yyyzz_1[i] * wa_y[i];

        g_yy_0_yyzzz_0[i] = g_0_0_yyzzz_0[i] * fbe_0 - g_0_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzzz_1[i] * fi_acd_0 + g_y_0_yyzzz_1[i] * wa_y[i];

        g_yy_0_yzzzz_0[i] = g_0_0_yzzzz_0[i] * fbe_0 - g_0_0_yzzzz_1[i] * fz_be_0 + g_y_0_zzzz_1[i] * fi_acd_0 + g_y_0_yzzzz_1[i] * wa_y[i];

        g_yy_0_zzzzz_0[i] = g_0_0_zzzzz_0[i] * fbe_0 - g_0_0_zzzzz_1[i] * fz_be_0 + g_y_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 84-105 components of targeted buffer : DSH

    auto g_yz_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 84);

    auto g_yz_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 85);

    auto g_yz_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 86);

    auto g_yz_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 87);

    auto g_yz_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 88);

    auto g_yz_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 89);

    auto g_yz_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 90);

    auto g_yz_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 91);

    auto g_yz_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 92);

    auto g_yz_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 93);

    auto g_yz_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 94);

    auto g_yz_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 95);

    auto g_yz_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 96);

    auto g_yz_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 97);

    auto g_yz_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 98);

    auto g_yz_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 99);

    auto g_yz_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 100);

    auto g_yz_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 101);

    auto g_yz_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 102);

    auto g_yz_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 103);

    auto g_yz_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 104);

    #pragma omp simd aligned(g_y_0_xxxxy_1, g_y_0_xxxyy_1, g_y_0_xxyyy_1, g_y_0_xyyyy_1, g_y_0_yyyyy_1, g_yz_0_xxxxx_0, g_yz_0_xxxxy_0, g_yz_0_xxxxz_0, g_yz_0_xxxyy_0, g_yz_0_xxxyz_0, g_yz_0_xxxzz_0, g_yz_0_xxyyy_0, g_yz_0_xxyyz_0, g_yz_0_xxyzz_0, g_yz_0_xxzzz_0, g_yz_0_xyyyy_0, g_yz_0_xyyyz_0, g_yz_0_xyyzz_0, g_yz_0_xyzzz_0, g_yz_0_xzzzz_0, g_yz_0_yyyyy_0, g_yz_0_yyyyz_0, g_yz_0_yyyzz_0, g_yz_0_yyzzz_0, g_yz_0_yzzzz_0, g_yz_0_zzzzz_0, g_z_0_xxxxx_1, g_z_0_xxxxz_1, g_z_0_xxxyz_1, g_z_0_xxxz_1, g_z_0_xxxzz_1, g_z_0_xxyyz_1, g_z_0_xxyz_1, g_z_0_xxyzz_1, g_z_0_xxzz_1, g_z_0_xxzzz_1, g_z_0_xyyyz_1, g_z_0_xyyz_1, g_z_0_xyyzz_1, g_z_0_xyzz_1, g_z_0_xyzzz_1, g_z_0_xzzz_1, g_z_0_xzzzz_1, g_z_0_yyyyz_1, g_z_0_yyyz_1, g_z_0_yyyzz_1, g_z_0_yyzz_1, g_z_0_yyzzz_1, g_z_0_yzzz_1, g_z_0_yzzzz_1, g_z_0_zzzz_1, g_z_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxxx_0[i] = g_z_0_xxxxx_1[i] * wa_y[i];

        g_yz_0_xxxxy_0[i] = g_y_0_xxxxy_1[i] * wa_z[i];

        g_yz_0_xxxxz_0[i] = g_z_0_xxxxz_1[i] * wa_y[i];

        g_yz_0_xxxyy_0[i] = g_y_0_xxxyy_1[i] * wa_z[i];

        g_yz_0_xxxyz_0[i] = g_z_0_xxxz_1[i] * fi_acd_0 + g_z_0_xxxyz_1[i] * wa_y[i];

        g_yz_0_xxxzz_0[i] = g_z_0_xxxzz_1[i] * wa_y[i];

        g_yz_0_xxyyy_0[i] = g_y_0_xxyyy_1[i] * wa_z[i];

        g_yz_0_xxyyz_0[i] = 2.0 * g_z_0_xxyz_1[i] * fi_acd_0 + g_z_0_xxyyz_1[i] * wa_y[i];

        g_yz_0_xxyzz_0[i] = g_z_0_xxzz_1[i] * fi_acd_0 + g_z_0_xxyzz_1[i] * wa_y[i];

        g_yz_0_xxzzz_0[i] = g_z_0_xxzzz_1[i] * wa_y[i];

        g_yz_0_xyyyy_0[i] = g_y_0_xyyyy_1[i] * wa_z[i];

        g_yz_0_xyyyz_0[i] = 3.0 * g_z_0_xyyz_1[i] * fi_acd_0 + g_z_0_xyyyz_1[i] * wa_y[i];

        g_yz_0_xyyzz_0[i] = 2.0 * g_z_0_xyzz_1[i] * fi_acd_0 + g_z_0_xyyzz_1[i] * wa_y[i];

        g_yz_0_xyzzz_0[i] = g_z_0_xzzz_1[i] * fi_acd_0 + g_z_0_xyzzz_1[i] * wa_y[i];

        g_yz_0_xzzzz_0[i] = g_z_0_xzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyy_0[i] = g_y_0_yyyyy_1[i] * wa_z[i];

        g_yz_0_yyyyz_0[i] = 4.0 * g_z_0_yyyz_1[i] * fi_acd_0 + g_z_0_yyyyz_1[i] * wa_y[i];

        g_yz_0_yyyzz_0[i] = 3.0 * g_z_0_yyzz_1[i] * fi_acd_0 + g_z_0_yyyzz_1[i] * wa_y[i];

        g_yz_0_yyzzz_0[i] = 2.0 * g_z_0_yzzz_1[i] * fi_acd_0 + g_z_0_yyzzz_1[i] * wa_y[i];

        g_yz_0_yzzzz_0[i] = g_z_0_zzzz_1[i] * fi_acd_0 + g_z_0_yzzzz_1[i] * wa_y[i];

        g_yz_0_zzzzz_0[i] = g_z_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 105-126 components of targeted buffer : DSH

    auto g_zz_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 105);

    auto g_zz_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 106);

    auto g_zz_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 107);

    auto g_zz_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 108);

    auto g_zz_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 109);

    auto g_zz_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 110);

    auto g_zz_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 111);

    auto g_zz_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 112);

    auto g_zz_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 113);

    auto g_zz_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 114);

    auto g_zz_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 115);

    auto g_zz_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 116);

    auto g_zz_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 117);

    auto g_zz_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 118);

    auto g_zz_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 119);

    auto g_zz_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 120);

    auto g_zz_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 121);

    auto g_zz_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 122);

    auto g_zz_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 123);

    auto g_zz_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 124);

    auto g_zz_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 125);

    #pragma omp simd aligned(g_0_0_xxxxx_0, g_0_0_xxxxx_1, g_0_0_xxxxy_0, g_0_0_xxxxy_1, g_0_0_xxxxz_0, g_0_0_xxxxz_1, g_0_0_xxxyy_0, g_0_0_xxxyy_1, g_0_0_xxxyz_0, g_0_0_xxxyz_1, g_0_0_xxxzz_0, g_0_0_xxxzz_1, g_0_0_xxyyy_0, g_0_0_xxyyy_1, g_0_0_xxyyz_0, g_0_0_xxyyz_1, g_0_0_xxyzz_0, g_0_0_xxyzz_1, g_0_0_xxzzz_0, g_0_0_xxzzz_1, g_0_0_xyyyy_0, g_0_0_xyyyy_1, g_0_0_xyyyz_0, g_0_0_xyyyz_1, g_0_0_xyyzz_0, g_0_0_xyyzz_1, g_0_0_xyzzz_0, g_0_0_xyzzz_1, g_0_0_xzzzz_0, g_0_0_xzzzz_1, g_0_0_yyyyy_0, g_0_0_yyyyy_1, g_0_0_yyyyz_0, g_0_0_yyyyz_1, g_0_0_yyyzz_0, g_0_0_yyyzz_1, g_0_0_yyzzz_0, g_0_0_yyzzz_1, g_0_0_yzzzz_0, g_0_0_yzzzz_1, g_0_0_zzzzz_0, g_0_0_zzzzz_1, g_z_0_xxxx_1, g_z_0_xxxxx_1, g_z_0_xxxxy_1, g_z_0_xxxxz_1, g_z_0_xxxy_1, g_z_0_xxxyy_1, g_z_0_xxxyz_1, g_z_0_xxxz_1, g_z_0_xxxzz_1, g_z_0_xxyy_1, g_z_0_xxyyy_1, g_z_0_xxyyz_1, g_z_0_xxyz_1, g_z_0_xxyzz_1, g_z_0_xxzz_1, g_z_0_xxzzz_1, g_z_0_xyyy_1, g_z_0_xyyyy_1, g_z_0_xyyyz_1, g_z_0_xyyz_1, g_z_0_xyyzz_1, g_z_0_xyzz_1, g_z_0_xyzzz_1, g_z_0_xzzz_1, g_z_0_xzzzz_1, g_z_0_yyyy_1, g_z_0_yyyyy_1, g_z_0_yyyyz_1, g_z_0_yyyz_1, g_z_0_yyyzz_1, g_z_0_yyzz_1, g_z_0_yyzzz_1, g_z_0_yzzz_1, g_z_0_yzzzz_1, g_z_0_zzzz_1, g_z_0_zzzzz_1, g_zz_0_xxxxx_0, g_zz_0_xxxxy_0, g_zz_0_xxxxz_0, g_zz_0_xxxyy_0, g_zz_0_xxxyz_0, g_zz_0_xxxzz_0, g_zz_0_xxyyy_0, g_zz_0_xxyyz_0, g_zz_0_xxyzz_0, g_zz_0_xxzzz_0, g_zz_0_xyyyy_0, g_zz_0_xyyyz_0, g_zz_0_xyyzz_0, g_zz_0_xyzzz_0, g_zz_0_xzzzz_0, g_zz_0_yyyyy_0, g_zz_0_yyyyz_0, g_zz_0_yyyzz_0, g_zz_0_yyzzz_0, g_zz_0_yzzzz_0, g_zz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxxx_0[i] = g_0_0_xxxxx_0[i] * fbe_0 - g_0_0_xxxxx_1[i] * fz_be_0 + g_z_0_xxxxx_1[i] * wa_z[i];

        g_zz_0_xxxxy_0[i] = g_0_0_xxxxy_0[i] * fbe_0 - g_0_0_xxxxy_1[i] * fz_be_0 + g_z_0_xxxxy_1[i] * wa_z[i];

        g_zz_0_xxxxz_0[i] = g_0_0_xxxxz_0[i] * fbe_0 - g_0_0_xxxxz_1[i] * fz_be_0 + g_z_0_xxxx_1[i] * fi_acd_0 + g_z_0_xxxxz_1[i] * wa_z[i];

        g_zz_0_xxxyy_0[i] = g_0_0_xxxyy_0[i] * fbe_0 - g_0_0_xxxyy_1[i] * fz_be_0 + g_z_0_xxxyy_1[i] * wa_z[i];

        g_zz_0_xxxyz_0[i] = g_0_0_xxxyz_0[i] * fbe_0 - g_0_0_xxxyz_1[i] * fz_be_0 + g_z_0_xxxy_1[i] * fi_acd_0 + g_z_0_xxxyz_1[i] * wa_z[i];

        g_zz_0_xxxzz_0[i] = g_0_0_xxxzz_0[i] * fbe_0 - g_0_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxz_1[i] * fi_acd_0 + g_z_0_xxxzz_1[i] * wa_z[i];

        g_zz_0_xxyyy_0[i] = g_0_0_xxyyy_0[i] * fbe_0 - g_0_0_xxyyy_1[i] * fz_be_0 + g_z_0_xxyyy_1[i] * wa_z[i];

        g_zz_0_xxyyz_0[i] = g_0_0_xxyyz_0[i] * fbe_0 - g_0_0_xxyyz_1[i] * fz_be_0 + g_z_0_xxyy_1[i] * fi_acd_0 + g_z_0_xxyyz_1[i] * wa_z[i];

        g_zz_0_xxyzz_0[i] = g_0_0_xxyzz_0[i] * fbe_0 - g_0_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxyz_1[i] * fi_acd_0 + g_z_0_xxyzz_1[i] * wa_z[i];

        g_zz_0_xxzzz_0[i] = g_0_0_xxzzz_0[i] * fbe_0 - g_0_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxzz_1[i] * fi_acd_0 + g_z_0_xxzzz_1[i] * wa_z[i];

        g_zz_0_xyyyy_0[i] = g_0_0_xyyyy_0[i] * fbe_0 - g_0_0_xyyyy_1[i] * fz_be_0 + g_z_0_xyyyy_1[i] * wa_z[i];

        g_zz_0_xyyyz_0[i] = g_0_0_xyyyz_0[i] * fbe_0 - g_0_0_xyyyz_1[i] * fz_be_0 + g_z_0_xyyy_1[i] * fi_acd_0 + g_z_0_xyyyz_1[i] * wa_z[i];

        g_zz_0_xyyzz_0[i] = g_0_0_xyyzz_0[i] * fbe_0 - g_0_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyyz_1[i] * fi_acd_0 + g_z_0_xyyzz_1[i] * wa_z[i];

        g_zz_0_xyzzz_0[i] = g_0_0_xyzzz_0[i] * fbe_0 - g_0_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xyzz_1[i] * fi_acd_0 + g_z_0_xyzzz_1[i] * wa_z[i];

        g_zz_0_xzzzz_0[i] = g_0_0_xzzzz_0[i] * fbe_0 - g_0_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xzzz_1[i] * fi_acd_0 + g_z_0_xzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyy_0[i] = g_0_0_yyyyy_0[i] * fbe_0 - g_0_0_yyyyy_1[i] * fz_be_0 + g_z_0_yyyyy_1[i] * wa_z[i];

        g_zz_0_yyyyz_0[i] = g_0_0_yyyyz_0[i] * fbe_0 - g_0_0_yyyyz_1[i] * fz_be_0 + g_z_0_yyyy_1[i] * fi_acd_0 + g_z_0_yyyyz_1[i] * wa_z[i];

        g_zz_0_yyyzz_0[i] = g_0_0_yyyzz_0[i] * fbe_0 - g_0_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyyz_1[i] * fi_acd_0 + g_z_0_yyyzz_1[i] * wa_z[i];

        g_zz_0_yyzzz_0[i] = g_0_0_yyzzz_0[i] * fbe_0 - g_0_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yyzz_1[i] * fi_acd_0 + g_z_0_yyzzz_1[i] * wa_z[i];

        g_zz_0_yzzzz_0[i] = g_0_0_yzzzz_0[i] * fbe_0 - g_0_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_yzzz_1[i] * fi_acd_0 + g_z_0_yzzzz_1[i] * wa_z[i];

        g_zz_0_zzzzz_0[i] = g_0_0_zzzzz_0[i] * fbe_0 - g_0_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_zzzz_1[i] * fi_acd_0 + g_z_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

