#include "ElectronRepulsionPrimRecSPSI.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_spsi(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_spsi,
                                  size_t                idx_eri_1_sssh,
                                  size_t                idx_eri_0_sssi,
                                  size_t                idx_eri_1_sssi,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_sssh);

    auto g_0_0_0_xxxxy_1 = pbuffer.data(idx_eri_1_sssh + 1);

    auto g_0_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_sssh + 2);

    auto g_0_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_sssh + 3);

    auto g_0_0_0_xxxyz_1 = pbuffer.data(idx_eri_1_sssh + 4);

    auto g_0_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_sssh + 5);

    auto g_0_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_sssh + 6);

    auto g_0_0_0_xxyyz_1 = pbuffer.data(idx_eri_1_sssh + 7);

    auto g_0_0_0_xxyzz_1 = pbuffer.data(idx_eri_1_sssh + 8);

    auto g_0_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_sssh + 9);

    auto g_0_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_sssh + 10);

    auto g_0_0_0_xyyyz_1 = pbuffer.data(idx_eri_1_sssh + 11);

    auto g_0_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_sssh + 12);

    auto g_0_0_0_xyzzz_1 = pbuffer.data(idx_eri_1_sssh + 13);

    auto g_0_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_sssh + 14);

    auto g_0_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_sssh + 15);

    auto g_0_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_sssh + 16);

    auto g_0_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_sssh + 17);

    auto g_0_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_sssh + 18);

    auto g_0_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_sssh + 19);

    auto g_0_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_sssh + 20);

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sssi);

    auto g_0_0_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sssi + 1);

    auto g_0_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sssi + 2);

    auto g_0_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sssi + 3);

    auto g_0_0_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sssi + 4);

    auto g_0_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sssi + 5);

    auto g_0_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sssi + 6);

    auto g_0_0_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sssi + 7);

    auto g_0_0_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sssi + 8);

    auto g_0_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sssi + 9);

    auto g_0_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sssi + 10);

    auto g_0_0_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sssi + 11);

    auto g_0_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sssi + 12);

    auto g_0_0_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sssi + 13);

    auto g_0_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sssi + 14);

    auto g_0_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 15);

    auto g_0_0_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 16);

    auto g_0_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 17);

    auto g_0_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 18);

    auto g_0_0_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 19);

    auto g_0_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 20);

    auto g_0_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 21);

    auto g_0_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 22);

    auto g_0_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 23);

    auto g_0_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 24);

    auto g_0_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 25);

    auto g_0_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 26);

    auto g_0_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 27);

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sssi);

    auto g_0_0_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sssi + 1);

    auto g_0_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sssi + 2);

    auto g_0_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sssi + 3);

    auto g_0_0_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sssi + 4);

    auto g_0_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sssi + 5);

    auto g_0_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sssi + 6);

    auto g_0_0_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sssi + 7);

    auto g_0_0_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sssi + 8);

    auto g_0_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sssi + 9);

    auto g_0_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sssi + 10);

    auto g_0_0_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sssi + 11);

    auto g_0_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sssi + 12);

    auto g_0_0_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sssi + 13);

    auto g_0_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sssi + 14);

    auto g_0_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 15);

    auto g_0_0_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sssi + 16);

    auto g_0_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 17);

    auto g_0_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 18);

    auto g_0_0_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sssi + 19);

    auto g_0_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 20);

    auto g_0_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 21);

    auto g_0_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sssi + 22);

    auto g_0_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 23);

    auto g_0_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 24);

    auto g_0_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sssi + 25);

    auto g_0_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 26);

    auto g_0_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 27);

    /// Set up 0-28 components of targeted buffer : SPSI

    auto g_0_x_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi);

    auto g_0_x_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 1);

    auto g_0_x_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 2);

    auto g_0_x_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 3);

    auto g_0_x_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 4);

    auto g_0_x_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 5);

    auto g_0_x_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 6);

    auto g_0_x_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 7);

    auto g_0_x_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 8);

    auto g_0_x_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 9);

    auto g_0_x_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 10);

    auto g_0_x_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 11);

    auto g_0_x_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 12);

    auto g_0_x_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 13);

    auto g_0_x_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 14);

    auto g_0_x_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 15);

    auto g_0_x_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 16);

    auto g_0_x_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 17);

    auto g_0_x_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 18);

    auto g_0_x_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 19);

    auto g_0_x_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 20);

    auto g_0_x_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 21);

    auto g_0_x_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 22);

    auto g_0_x_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 23);

    auto g_0_x_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 24);

    auto g_0_x_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 25);

    auto g_0_x_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 26);

    auto g_0_x_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 27);

#pragma omp simd aligned(g_0_0_0_xxxxx_1,      \
                             g_0_0_0_xxxxxx_0, \
                             g_0_0_0_xxxxxx_1, \
                             g_0_0_0_xxxxxy_0, \
                             g_0_0_0_xxxxxy_1, \
                             g_0_0_0_xxxxxz_0, \
                             g_0_0_0_xxxxxz_1, \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxyy_0, \
                             g_0_0_0_xxxxyy_1, \
                             g_0_0_0_xxxxyz_0, \
                             g_0_0_0_xxxxyz_1, \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxxzz_0, \
                             g_0_0_0_xxxxzz_1, \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyyy_0, \
                             g_0_0_0_xxxyyy_1, \
                             g_0_0_0_xxxyyz_0, \
                             g_0_0_0_xxxyyz_1, \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxyzz_0, \
                             g_0_0_0_xxxyzz_1, \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxxzzz_0, \
                             g_0_0_0_xxxzzz_1, \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyyy_0, \
                             g_0_0_0_xxyyyy_1, \
                             g_0_0_0_xxyyyz_0, \
                             g_0_0_0_xxyyyz_1, \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyyzz_0, \
                             g_0_0_0_xxyyzz_1, \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxyzzz_0, \
                             g_0_0_0_xxyzzz_1, \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xxzzzz_0, \
                             g_0_0_0_xxzzzz_1, \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyyy_0, \
                             g_0_0_0_xyyyyy_1, \
                             g_0_0_0_xyyyyz_0, \
                             g_0_0_0_xyyyyz_1, \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyyzz_0, \
                             g_0_0_0_xyyyzz_1, \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyyzzz_0, \
                             g_0_0_0_xyyzzz_1, \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xyzzzz_0, \
                             g_0_0_0_xyzzzz_1, \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_xzzzzz_0, \
                             g_0_0_0_xzzzzz_1, \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyyy_0, \
                             g_0_0_0_yyyyyy_1, \
                             g_0_0_0_yyyyyz_0, \
                             g_0_0_0_yyyyyz_1, \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyyzz_0, \
                             g_0_0_0_yyyyzz_1, \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyyzzz_0, \
                             g_0_0_0_yyyzzz_1, \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yyzzzz_0, \
                             g_0_0_0_yyzzzz_1, \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_yzzzzz_0, \
                             g_0_0_0_yzzzzz_1, \
                             g_0_0_0_zzzzz_1,  \
                             g_0_0_0_zzzzzz_0, \
                             g_0_0_0_zzzzzz_1, \
                             g_0_x_0_xxxxxx_0, \
                             g_0_x_0_xxxxxy_0, \
                             g_0_x_0_xxxxxz_0, \
                             g_0_x_0_xxxxyy_0, \
                             g_0_x_0_xxxxyz_0, \
                             g_0_x_0_xxxxzz_0, \
                             g_0_x_0_xxxyyy_0, \
                             g_0_x_0_xxxyyz_0, \
                             g_0_x_0_xxxyzz_0, \
                             g_0_x_0_xxxzzz_0, \
                             g_0_x_0_xxyyyy_0, \
                             g_0_x_0_xxyyyz_0, \
                             g_0_x_0_xxyyzz_0, \
                             g_0_x_0_xxyzzz_0, \
                             g_0_x_0_xxzzzz_0, \
                             g_0_x_0_xyyyyy_0, \
                             g_0_x_0_xyyyyz_0, \
                             g_0_x_0_xyyyzz_0, \
                             g_0_x_0_xyyzzz_0, \
                             g_0_x_0_xyzzzz_0, \
                             g_0_x_0_xzzzzz_0, \
                             g_0_x_0_yyyyyy_0, \
                             g_0_x_0_yyyyyz_0, \
                             g_0_x_0_yyyyzz_0, \
                             g_0_x_0_yyyzzz_0, \
                             g_0_x_0_yyzzzz_0, \
                             g_0_x_0_yzzzzz_0, \
                             g_0_x_0_zzzzzz_0, \
                             wp_x : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxxxxx_0[i] = 6.0 * g_0_0_0_xxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxx_0[i] * pb_x + g_0_0_0_xxxxxx_1[i] * wp_x[i];

        g_0_x_0_xxxxxy_0[i] = 5.0 * g_0_0_0_xxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxy_0[i] * pb_x + g_0_0_0_xxxxxy_1[i] * wp_x[i];

        g_0_x_0_xxxxxz_0[i] = 5.0 * g_0_0_0_xxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxz_0[i] * pb_x + g_0_0_0_xxxxxz_1[i] * wp_x[i];

        g_0_x_0_xxxxyy_0[i] = 4.0 * g_0_0_0_xxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyy_0[i] * pb_x + g_0_0_0_xxxxyy_1[i] * wp_x[i];

        g_0_x_0_xxxxyz_0[i] = 4.0 * g_0_0_0_xxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyz_0[i] * pb_x + g_0_0_0_xxxxyz_1[i] * wp_x[i];

        g_0_x_0_xxxxzz_0[i] = 4.0 * g_0_0_0_xxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzz_0[i] * pb_x + g_0_0_0_xxxxzz_1[i] * wp_x[i];

        g_0_x_0_xxxyyy_0[i] = 3.0 * g_0_0_0_xxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyy_0[i] * pb_x + g_0_0_0_xxxyyy_1[i] * wp_x[i];

        g_0_x_0_xxxyyz_0[i] = 3.0 * g_0_0_0_xxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyz_0[i] * pb_x + g_0_0_0_xxxyyz_1[i] * wp_x[i];

        g_0_x_0_xxxyzz_0[i] = 3.0 * g_0_0_0_xxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzz_0[i] * pb_x + g_0_0_0_xxxyzz_1[i] * wp_x[i];

        g_0_x_0_xxxzzz_0[i] = 3.0 * g_0_0_0_xxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzz_0[i] * pb_x + g_0_0_0_xxxzzz_1[i] * wp_x[i];

        g_0_x_0_xxyyyy_0[i] = 2.0 * g_0_0_0_xyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyy_0[i] * pb_x + g_0_0_0_xxyyyy_1[i] * wp_x[i];

        g_0_x_0_xxyyyz_0[i] = 2.0 * g_0_0_0_xyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyz_0[i] * pb_x + g_0_0_0_xxyyyz_1[i] * wp_x[i];

        g_0_x_0_xxyyzz_0[i] = 2.0 * g_0_0_0_xyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzz_0[i] * pb_x + g_0_0_0_xxyyzz_1[i] * wp_x[i];

        g_0_x_0_xxyzzz_0[i] = 2.0 * g_0_0_0_xyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzz_0[i] * pb_x + g_0_0_0_xxyzzz_1[i] * wp_x[i];

        g_0_x_0_xxzzzz_0[i] = 2.0 * g_0_0_0_xzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzz_0[i] * pb_x + g_0_0_0_xxzzzz_1[i] * wp_x[i];

        g_0_x_0_xyyyyy_0[i] = g_0_0_0_yyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyy_0[i] * pb_x + g_0_0_0_xyyyyy_1[i] * wp_x[i];

        g_0_x_0_xyyyyz_0[i] = g_0_0_0_yyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyz_0[i] * pb_x + g_0_0_0_xyyyyz_1[i] * wp_x[i];

        g_0_x_0_xyyyzz_0[i] = g_0_0_0_yyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzz_0[i] * pb_x + g_0_0_0_xyyyzz_1[i] * wp_x[i];

        g_0_x_0_xyyzzz_0[i] = g_0_0_0_yyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzz_0[i] * pb_x + g_0_0_0_xyyzzz_1[i] * wp_x[i];

        g_0_x_0_xyzzzz_0[i] = g_0_0_0_yzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzz_0[i] * pb_x + g_0_0_0_xyzzzz_1[i] * wp_x[i];

        g_0_x_0_xzzzzz_0[i] = g_0_0_0_zzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzz_0[i] * pb_x + g_0_0_0_xzzzzz_1[i] * wp_x[i];

        g_0_x_0_yyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * pb_x + g_0_0_0_yyyyyy_1[i] * wp_x[i];

        g_0_x_0_yyyyyz_0[i] = g_0_0_0_yyyyyz_0[i] * pb_x + g_0_0_0_yyyyyz_1[i] * wp_x[i];

        g_0_x_0_yyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * pb_x + g_0_0_0_yyyyzz_1[i] * wp_x[i];

        g_0_x_0_yyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * pb_x + g_0_0_0_yyyzzz_1[i] * wp_x[i];

        g_0_x_0_yyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * pb_x + g_0_0_0_yyzzzz_1[i] * wp_x[i];

        g_0_x_0_yzzzzz_0[i] = g_0_0_0_yzzzzz_0[i] * pb_x + g_0_0_0_yzzzzz_1[i] * wp_x[i];

        g_0_x_0_zzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * pb_x + g_0_0_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SPSI

    auto g_0_y_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi + 28);

    auto g_0_y_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 29);

    auto g_0_y_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 30);

    auto g_0_y_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 31);

    auto g_0_y_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 32);

    auto g_0_y_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 33);

    auto g_0_y_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 34);

    auto g_0_y_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 35);

    auto g_0_y_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 36);

    auto g_0_y_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 37);

    auto g_0_y_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 38);

    auto g_0_y_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 39);

    auto g_0_y_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 40);

    auto g_0_y_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 41);

    auto g_0_y_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 42);

    auto g_0_y_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 43);

    auto g_0_y_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 44);

    auto g_0_y_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 45);

    auto g_0_y_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 46);

    auto g_0_y_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 47);

    auto g_0_y_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 48);

    auto g_0_y_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 49);

    auto g_0_y_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 50);

    auto g_0_y_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 51);

    auto g_0_y_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 52);

    auto g_0_y_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 53);

    auto g_0_y_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 54);

    auto g_0_y_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 55);

#pragma omp simd aligned(g_0_0_0_xxxxx_1,      \
                             g_0_0_0_xxxxxx_0, \
                             g_0_0_0_xxxxxx_1, \
                             g_0_0_0_xxxxxy_0, \
                             g_0_0_0_xxxxxy_1, \
                             g_0_0_0_xxxxxz_0, \
                             g_0_0_0_xxxxxz_1, \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxyy_0, \
                             g_0_0_0_xxxxyy_1, \
                             g_0_0_0_xxxxyz_0, \
                             g_0_0_0_xxxxyz_1, \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxxzz_0, \
                             g_0_0_0_xxxxzz_1, \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyyy_0, \
                             g_0_0_0_xxxyyy_1, \
                             g_0_0_0_xxxyyz_0, \
                             g_0_0_0_xxxyyz_1, \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxyzz_0, \
                             g_0_0_0_xxxyzz_1, \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxxzzz_0, \
                             g_0_0_0_xxxzzz_1, \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyyy_0, \
                             g_0_0_0_xxyyyy_1, \
                             g_0_0_0_xxyyyz_0, \
                             g_0_0_0_xxyyyz_1, \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyyzz_0, \
                             g_0_0_0_xxyyzz_1, \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxyzzz_0, \
                             g_0_0_0_xxyzzz_1, \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xxzzzz_0, \
                             g_0_0_0_xxzzzz_1, \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyyy_0, \
                             g_0_0_0_xyyyyy_1, \
                             g_0_0_0_xyyyyz_0, \
                             g_0_0_0_xyyyyz_1, \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyyzz_0, \
                             g_0_0_0_xyyyzz_1, \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyyzzz_0, \
                             g_0_0_0_xyyzzz_1, \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xyzzzz_0, \
                             g_0_0_0_xyzzzz_1, \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_xzzzzz_0, \
                             g_0_0_0_xzzzzz_1, \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyyy_0, \
                             g_0_0_0_yyyyyy_1, \
                             g_0_0_0_yyyyyz_0, \
                             g_0_0_0_yyyyyz_1, \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyyzz_0, \
                             g_0_0_0_yyyyzz_1, \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyyzzz_0, \
                             g_0_0_0_yyyzzz_1, \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yyzzzz_0, \
                             g_0_0_0_yyzzzz_1, \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_yzzzzz_0, \
                             g_0_0_0_yzzzzz_1, \
                             g_0_0_0_zzzzz_1,  \
                             g_0_0_0_zzzzzz_0, \
                             g_0_0_0_zzzzzz_1, \
                             g_0_y_0_xxxxxx_0, \
                             g_0_y_0_xxxxxy_0, \
                             g_0_y_0_xxxxxz_0, \
                             g_0_y_0_xxxxyy_0, \
                             g_0_y_0_xxxxyz_0, \
                             g_0_y_0_xxxxzz_0, \
                             g_0_y_0_xxxyyy_0, \
                             g_0_y_0_xxxyyz_0, \
                             g_0_y_0_xxxyzz_0, \
                             g_0_y_0_xxxzzz_0, \
                             g_0_y_0_xxyyyy_0, \
                             g_0_y_0_xxyyyz_0, \
                             g_0_y_0_xxyyzz_0, \
                             g_0_y_0_xxyzzz_0, \
                             g_0_y_0_xxzzzz_0, \
                             g_0_y_0_xyyyyy_0, \
                             g_0_y_0_xyyyyz_0, \
                             g_0_y_0_xyyyzz_0, \
                             g_0_y_0_xyyzzz_0, \
                             g_0_y_0_xyzzzz_0, \
                             g_0_y_0_xzzzzz_0, \
                             g_0_y_0_yyyyyy_0, \
                             g_0_y_0_yyyyyz_0, \
                             g_0_y_0_yyyyzz_0, \
                             g_0_y_0_yyyzzz_0, \
                             g_0_y_0_yyzzzz_0, \
                             g_0_y_0_yzzzzz_0, \
                             g_0_y_0_zzzzzz_0, \
                             wp_y : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxxxxx_0[i] = g_0_0_0_xxxxxx_0[i] * pb_y + g_0_0_0_xxxxxx_1[i] * wp_y[i];

        g_0_y_0_xxxxxy_0[i] = g_0_0_0_xxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxy_0[i] * pb_y + g_0_0_0_xxxxxy_1[i] * wp_y[i];

        g_0_y_0_xxxxxz_0[i] = g_0_0_0_xxxxxz_0[i] * pb_y + g_0_0_0_xxxxxz_1[i] * wp_y[i];

        g_0_y_0_xxxxyy_0[i] = 2.0 * g_0_0_0_xxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyy_0[i] * pb_y + g_0_0_0_xxxxyy_1[i] * wp_y[i];

        g_0_y_0_xxxxyz_0[i] = g_0_0_0_xxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyz_0[i] * pb_y + g_0_0_0_xxxxyz_1[i] * wp_y[i];

        g_0_y_0_xxxxzz_0[i] = g_0_0_0_xxxxzz_0[i] * pb_y + g_0_0_0_xxxxzz_1[i] * wp_y[i];

        g_0_y_0_xxxyyy_0[i] = 3.0 * g_0_0_0_xxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyy_0[i] * pb_y + g_0_0_0_xxxyyy_1[i] * wp_y[i];

        g_0_y_0_xxxyyz_0[i] = 2.0 * g_0_0_0_xxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyz_0[i] * pb_y + g_0_0_0_xxxyyz_1[i] * wp_y[i];

        g_0_y_0_xxxyzz_0[i] = g_0_0_0_xxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzz_0[i] * pb_y + g_0_0_0_xxxyzz_1[i] * wp_y[i];

        g_0_y_0_xxxzzz_0[i] = g_0_0_0_xxxzzz_0[i] * pb_y + g_0_0_0_xxxzzz_1[i] * wp_y[i];

        g_0_y_0_xxyyyy_0[i] = 4.0 * g_0_0_0_xxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyy_0[i] * pb_y + g_0_0_0_xxyyyy_1[i] * wp_y[i];

        g_0_y_0_xxyyyz_0[i] = 3.0 * g_0_0_0_xxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyz_0[i] * pb_y + g_0_0_0_xxyyyz_1[i] * wp_y[i];

        g_0_y_0_xxyyzz_0[i] = 2.0 * g_0_0_0_xxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzz_0[i] * pb_y + g_0_0_0_xxyyzz_1[i] * wp_y[i];

        g_0_y_0_xxyzzz_0[i] = g_0_0_0_xxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzz_0[i] * pb_y + g_0_0_0_xxyzzz_1[i] * wp_y[i];

        g_0_y_0_xxzzzz_0[i] = g_0_0_0_xxzzzz_0[i] * pb_y + g_0_0_0_xxzzzz_1[i] * wp_y[i];

        g_0_y_0_xyyyyy_0[i] = 5.0 * g_0_0_0_xyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyy_0[i] * pb_y + g_0_0_0_xyyyyy_1[i] * wp_y[i];

        g_0_y_0_xyyyyz_0[i] = 4.0 * g_0_0_0_xyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyz_0[i] * pb_y + g_0_0_0_xyyyyz_1[i] * wp_y[i];

        g_0_y_0_xyyyzz_0[i] = 3.0 * g_0_0_0_xyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzz_0[i] * pb_y + g_0_0_0_xyyyzz_1[i] * wp_y[i];

        g_0_y_0_xyyzzz_0[i] = 2.0 * g_0_0_0_xyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzz_0[i] * pb_y + g_0_0_0_xyyzzz_1[i] * wp_y[i];

        g_0_y_0_xyzzzz_0[i] = g_0_0_0_xzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzz_0[i] * pb_y + g_0_0_0_xyzzzz_1[i] * wp_y[i];

        g_0_y_0_xzzzzz_0[i] = g_0_0_0_xzzzzz_0[i] * pb_y + g_0_0_0_xzzzzz_1[i] * wp_y[i];

        g_0_y_0_yyyyyy_0[i] = 6.0 * g_0_0_0_yyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyy_0[i] * pb_y + g_0_0_0_yyyyyy_1[i] * wp_y[i];

        g_0_y_0_yyyyyz_0[i] = 5.0 * g_0_0_0_yyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyz_0[i] * pb_y + g_0_0_0_yyyyyz_1[i] * wp_y[i];

        g_0_y_0_yyyyzz_0[i] = 4.0 * g_0_0_0_yyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzz_0[i] * pb_y + g_0_0_0_yyyyzz_1[i] * wp_y[i];

        g_0_y_0_yyyzzz_0[i] = 3.0 * g_0_0_0_yyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzz_0[i] * pb_y + g_0_0_0_yyyzzz_1[i] * wp_y[i];

        g_0_y_0_yyzzzz_0[i] = 2.0 * g_0_0_0_yzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzz_0[i] * pb_y + g_0_0_0_yyzzzz_1[i] * wp_y[i];

        g_0_y_0_yzzzzz_0[i] = g_0_0_0_zzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzz_0[i] * pb_y + g_0_0_0_yzzzzz_1[i] * wp_y[i];

        g_0_y_0_zzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * pb_y + g_0_0_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SPSI

    auto g_0_z_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi + 56);

    auto g_0_z_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 57);

    auto g_0_z_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 58);

    auto g_0_z_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 59);

    auto g_0_z_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 60);

    auto g_0_z_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 61);

    auto g_0_z_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 62);

    auto g_0_z_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 63);

    auto g_0_z_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 64);

    auto g_0_z_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 65);

    auto g_0_z_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 66);

    auto g_0_z_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 67);

    auto g_0_z_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 68);

    auto g_0_z_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 69);

    auto g_0_z_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 70);

    auto g_0_z_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 71);

    auto g_0_z_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 72);

    auto g_0_z_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 73);

    auto g_0_z_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 74);

    auto g_0_z_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 75);

    auto g_0_z_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 76);

    auto g_0_z_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 77);

    auto g_0_z_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 78);

    auto g_0_z_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 79);

    auto g_0_z_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 80);

    auto g_0_z_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 81);

    auto g_0_z_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 82);

    auto g_0_z_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 83);

#pragma omp simd aligned(g_0_0_0_xxxxx_1,      \
                             g_0_0_0_xxxxxx_0, \
                             g_0_0_0_xxxxxx_1, \
                             g_0_0_0_xxxxxy_0, \
                             g_0_0_0_xxxxxy_1, \
                             g_0_0_0_xxxxxz_0, \
                             g_0_0_0_xxxxxz_1, \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxyy_0, \
                             g_0_0_0_xxxxyy_1, \
                             g_0_0_0_xxxxyz_0, \
                             g_0_0_0_xxxxyz_1, \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxxzz_0, \
                             g_0_0_0_xxxxzz_1, \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyyy_0, \
                             g_0_0_0_xxxyyy_1, \
                             g_0_0_0_xxxyyz_0, \
                             g_0_0_0_xxxyyz_1, \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxyzz_0, \
                             g_0_0_0_xxxyzz_1, \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxxzzz_0, \
                             g_0_0_0_xxxzzz_1, \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyyy_0, \
                             g_0_0_0_xxyyyy_1, \
                             g_0_0_0_xxyyyz_0, \
                             g_0_0_0_xxyyyz_1, \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyyzz_0, \
                             g_0_0_0_xxyyzz_1, \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxyzzz_0, \
                             g_0_0_0_xxyzzz_1, \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xxzzzz_0, \
                             g_0_0_0_xxzzzz_1, \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyyy_0, \
                             g_0_0_0_xyyyyy_1, \
                             g_0_0_0_xyyyyz_0, \
                             g_0_0_0_xyyyyz_1, \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyyzz_0, \
                             g_0_0_0_xyyyzz_1, \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyyzzz_0, \
                             g_0_0_0_xyyzzz_1, \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xyzzzz_0, \
                             g_0_0_0_xyzzzz_1, \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_xzzzzz_0, \
                             g_0_0_0_xzzzzz_1, \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyyy_0, \
                             g_0_0_0_yyyyyy_1, \
                             g_0_0_0_yyyyyz_0, \
                             g_0_0_0_yyyyyz_1, \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyyzz_0, \
                             g_0_0_0_yyyyzz_1, \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyyzzz_0, \
                             g_0_0_0_yyyzzz_1, \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yyzzzz_0, \
                             g_0_0_0_yyzzzz_1, \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_yzzzzz_0, \
                             g_0_0_0_yzzzzz_1, \
                             g_0_0_0_zzzzz_1,  \
                             g_0_0_0_zzzzzz_0, \
                             g_0_0_0_zzzzzz_1, \
                             g_0_z_0_xxxxxx_0, \
                             g_0_z_0_xxxxxy_0, \
                             g_0_z_0_xxxxxz_0, \
                             g_0_z_0_xxxxyy_0, \
                             g_0_z_0_xxxxyz_0, \
                             g_0_z_0_xxxxzz_0, \
                             g_0_z_0_xxxyyy_0, \
                             g_0_z_0_xxxyyz_0, \
                             g_0_z_0_xxxyzz_0, \
                             g_0_z_0_xxxzzz_0, \
                             g_0_z_0_xxyyyy_0, \
                             g_0_z_0_xxyyyz_0, \
                             g_0_z_0_xxyyzz_0, \
                             g_0_z_0_xxyzzz_0, \
                             g_0_z_0_xxzzzz_0, \
                             g_0_z_0_xyyyyy_0, \
                             g_0_z_0_xyyyyz_0, \
                             g_0_z_0_xyyyzz_0, \
                             g_0_z_0_xyyzzz_0, \
                             g_0_z_0_xyzzzz_0, \
                             g_0_z_0_xzzzzz_0, \
                             g_0_z_0_yyyyyy_0, \
                             g_0_z_0_yyyyyz_0, \
                             g_0_z_0_yyyyzz_0, \
                             g_0_z_0_yyyzzz_0, \
                             g_0_z_0_yyzzzz_0, \
                             g_0_z_0_yzzzzz_0, \
                             g_0_z_0_zzzzzz_0, \
                             wp_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxxxxx_0[i] = g_0_0_0_xxxxxx_0[i] * pb_z + g_0_0_0_xxxxxx_1[i] * wp_z[i];

        g_0_z_0_xxxxxy_0[i] = g_0_0_0_xxxxxy_0[i] * pb_z + g_0_0_0_xxxxxy_1[i] * wp_z[i];

        g_0_z_0_xxxxxz_0[i] = g_0_0_0_xxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxz_0[i] * pb_z + g_0_0_0_xxxxxz_1[i] * wp_z[i];

        g_0_z_0_xxxxyy_0[i] = g_0_0_0_xxxxyy_0[i] * pb_z + g_0_0_0_xxxxyy_1[i] * wp_z[i];

        g_0_z_0_xxxxyz_0[i] = g_0_0_0_xxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyz_0[i] * pb_z + g_0_0_0_xxxxyz_1[i] * wp_z[i];

        g_0_z_0_xxxxzz_0[i] = 2.0 * g_0_0_0_xxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzz_0[i] * pb_z + g_0_0_0_xxxxzz_1[i] * wp_z[i];

        g_0_z_0_xxxyyy_0[i] = g_0_0_0_xxxyyy_0[i] * pb_z + g_0_0_0_xxxyyy_1[i] * wp_z[i];

        g_0_z_0_xxxyyz_0[i] = g_0_0_0_xxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyz_0[i] * pb_z + g_0_0_0_xxxyyz_1[i] * wp_z[i];

        g_0_z_0_xxxyzz_0[i] = 2.0 * g_0_0_0_xxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzz_0[i] * pb_z + g_0_0_0_xxxyzz_1[i] * wp_z[i];

        g_0_z_0_xxxzzz_0[i] = 3.0 * g_0_0_0_xxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzz_0[i] * pb_z + g_0_0_0_xxxzzz_1[i] * wp_z[i];

        g_0_z_0_xxyyyy_0[i] = g_0_0_0_xxyyyy_0[i] * pb_z + g_0_0_0_xxyyyy_1[i] * wp_z[i];

        g_0_z_0_xxyyyz_0[i] = g_0_0_0_xxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyz_0[i] * pb_z + g_0_0_0_xxyyyz_1[i] * wp_z[i];

        g_0_z_0_xxyyzz_0[i] = 2.0 * g_0_0_0_xxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzz_0[i] * pb_z + g_0_0_0_xxyyzz_1[i] * wp_z[i];

        g_0_z_0_xxyzzz_0[i] = 3.0 * g_0_0_0_xxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzz_0[i] * pb_z + g_0_0_0_xxyzzz_1[i] * wp_z[i];

        g_0_z_0_xxzzzz_0[i] = 4.0 * g_0_0_0_xxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzz_0[i] * pb_z + g_0_0_0_xxzzzz_1[i] * wp_z[i];

        g_0_z_0_xyyyyy_0[i] = g_0_0_0_xyyyyy_0[i] * pb_z + g_0_0_0_xyyyyy_1[i] * wp_z[i];

        g_0_z_0_xyyyyz_0[i] = g_0_0_0_xyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyz_0[i] * pb_z + g_0_0_0_xyyyyz_1[i] * wp_z[i];

        g_0_z_0_xyyyzz_0[i] = 2.0 * g_0_0_0_xyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzz_0[i] * pb_z + g_0_0_0_xyyyzz_1[i] * wp_z[i];

        g_0_z_0_xyyzzz_0[i] = 3.0 * g_0_0_0_xyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzz_0[i] * pb_z + g_0_0_0_xyyzzz_1[i] * wp_z[i];

        g_0_z_0_xyzzzz_0[i] = 4.0 * g_0_0_0_xyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzz_0[i] * pb_z + g_0_0_0_xyzzzz_1[i] * wp_z[i];

        g_0_z_0_xzzzzz_0[i] = 5.0 * g_0_0_0_xzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzz_0[i] * pb_z + g_0_0_0_xzzzzz_1[i] * wp_z[i];

        g_0_z_0_yyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * pb_z + g_0_0_0_yyyyyy_1[i] * wp_z[i];

        g_0_z_0_yyyyyz_0[i] = g_0_0_0_yyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyz_0[i] * pb_z + g_0_0_0_yyyyyz_1[i] * wp_z[i];

        g_0_z_0_yyyyzz_0[i] = 2.0 * g_0_0_0_yyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzz_0[i] * pb_z + g_0_0_0_yyyyzz_1[i] * wp_z[i];

        g_0_z_0_yyyzzz_0[i] = 3.0 * g_0_0_0_yyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzz_0[i] * pb_z + g_0_0_0_yyyzzz_1[i] * wp_z[i];

        g_0_z_0_yyzzzz_0[i] = 4.0 * g_0_0_0_yyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzz_0[i] * pb_z + g_0_0_0_yyzzzz_1[i] * wp_z[i];

        g_0_z_0_yzzzzz_0[i] = 5.0 * g_0_0_0_yzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzz_0[i] * pb_z + g_0_0_0_yzzzzz_1[i] * wp_z[i];

        g_0_z_0_zzzzzz_0[i] = 6.0 * g_0_0_0_zzzzz_1[i] * fi_abcd_0 + g_0_0_0_zzzzzz_0[i] * pb_z + g_0_0_0_zzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
