#include "ElectronRepulsionContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pixx(CSimdArray<double>& contr_buffer_pixx,
                                     const CSimdArray<double>& contr_buffer_sixx,
                                     const CSimdArray<double>& contr_buffer_skxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_pixx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_sixx

            const auto si_off = i * dcomps + j;

            auto g_0_xxxxxx = contr_buffer_sixx[si_off + 0 * ccomps * dcomps];

            auto g_0_xxxxxy = contr_buffer_sixx[si_off + 1 * ccomps * dcomps];

            auto g_0_xxxxxz = contr_buffer_sixx[si_off + 2 * ccomps * dcomps];

            auto g_0_xxxxyy = contr_buffer_sixx[si_off + 3 * ccomps * dcomps];

            auto g_0_xxxxyz = contr_buffer_sixx[si_off + 4 * ccomps * dcomps];

            auto g_0_xxxxzz = contr_buffer_sixx[si_off + 5 * ccomps * dcomps];

            auto g_0_xxxyyy = contr_buffer_sixx[si_off + 6 * ccomps * dcomps];

            auto g_0_xxxyyz = contr_buffer_sixx[si_off + 7 * ccomps * dcomps];

            auto g_0_xxxyzz = contr_buffer_sixx[si_off + 8 * ccomps * dcomps];

            auto g_0_xxxzzz = contr_buffer_sixx[si_off + 9 * ccomps * dcomps];

            auto g_0_xxyyyy = contr_buffer_sixx[si_off + 10 * ccomps * dcomps];

            auto g_0_xxyyyz = contr_buffer_sixx[si_off + 11 * ccomps * dcomps];

            auto g_0_xxyyzz = contr_buffer_sixx[si_off + 12 * ccomps * dcomps];

            auto g_0_xxyzzz = contr_buffer_sixx[si_off + 13 * ccomps * dcomps];

            auto g_0_xxzzzz = contr_buffer_sixx[si_off + 14 * ccomps * dcomps];

            auto g_0_xyyyyy = contr_buffer_sixx[si_off + 15 * ccomps * dcomps];

            auto g_0_xyyyyz = contr_buffer_sixx[si_off + 16 * ccomps * dcomps];

            auto g_0_xyyyzz = contr_buffer_sixx[si_off + 17 * ccomps * dcomps];

            auto g_0_xyyzzz = contr_buffer_sixx[si_off + 18 * ccomps * dcomps];

            auto g_0_xyzzzz = contr_buffer_sixx[si_off + 19 * ccomps * dcomps];

            auto g_0_xzzzzz = contr_buffer_sixx[si_off + 20 * ccomps * dcomps];

            auto g_0_yyyyyy = contr_buffer_sixx[si_off + 21 * ccomps * dcomps];

            auto g_0_yyyyyz = contr_buffer_sixx[si_off + 22 * ccomps * dcomps];

            auto g_0_yyyyzz = contr_buffer_sixx[si_off + 23 * ccomps * dcomps];

            auto g_0_yyyzzz = contr_buffer_sixx[si_off + 24 * ccomps * dcomps];

            auto g_0_yyzzzz = contr_buffer_sixx[si_off + 25 * ccomps * dcomps];

            auto g_0_yzzzzz = contr_buffer_sixx[si_off + 26 * ccomps * dcomps];

            auto g_0_zzzzzz = contr_buffer_sixx[si_off + 27 * ccomps * dcomps];

            /// Set up components of auxilary buffer : contr_buffer_skxx

            const auto sk_off = i * dcomps + j;

            auto g_0_xxxxxxx = contr_buffer_skxx[sk_off + 0 * ccomps * dcomps];

            auto g_0_xxxxxxy = contr_buffer_skxx[sk_off + 1 * ccomps * dcomps];

            auto g_0_xxxxxxz = contr_buffer_skxx[sk_off + 2 * ccomps * dcomps];

            auto g_0_xxxxxyy = contr_buffer_skxx[sk_off + 3 * ccomps * dcomps];

            auto g_0_xxxxxyz = contr_buffer_skxx[sk_off + 4 * ccomps * dcomps];

            auto g_0_xxxxxzz = contr_buffer_skxx[sk_off + 5 * ccomps * dcomps];

            auto g_0_xxxxyyy = contr_buffer_skxx[sk_off + 6 * ccomps * dcomps];

            auto g_0_xxxxyyz = contr_buffer_skxx[sk_off + 7 * ccomps * dcomps];

            auto g_0_xxxxyzz = contr_buffer_skxx[sk_off + 8 * ccomps * dcomps];

            auto g_0_xxxxzzz = contr_buffer_skxx[sk_off + 9 * ccomps * dcomps];

            auto g_0_xxxyyyy = contr_buffer_skxx[sk_off + 10 * ccomps * dcomps];

            auto g_0_xxxyyyz = contr_buffer_skxx[sk_off + 11 * ccomps * dcomps];

            auto g_0_xxxyyzz = contr_buffer_skxx[sk_off + 12 * ccomps * dcomps];

            auto g_0_xxxyzzz = contr_buffer_skxx[sk_off + 13 * ccomps * dcomps];

            auto g_0_xxxzzzz = contr_buffer_skxx[sk_off + 14 * ccomps * dcomps];

            auto g_0_xxyyyyy = contr_buffer_skxx[sk_off + 15 * ccomps * dcomps];

            auto g_0_xxyyyyz = contr_buffer_skxx[sk_off + 16 * ccomps * dcomps];

            auto g_0_xxyyyzz = contr_buffer_skxx[sk_off + 17 * ccomps * dcomps];

            auto g_0_xxyyzzz = contr_buffer_skxx[sk_off + 18 * ccomps * dcomps];

            auto g_0_xxyzzzz = contr_buffer_skxx[sk_off + 19 * ccomps * dcomps];

            auto g_0_xxzzzzz = contr_buffer_skxx[sk_off + 20 * ccomps * dcomps];

            auto g_0_xyyyyyy = contr_buffer_skxx[sk_off + 21 * ccomps * dcomps];

            auto g_0_xyyyyyz = contr_buffer_skxx[sk_off + 22 * ccomps * dcomps];

            auto g_0_xyyyyzz = contr_buffer_skxx[sk_off + 23 * ccomps * dcomps];

            auto g_0_xyyyzzz = contr_buffer_skxx[sk_off + 24 * ccomps * dcomps];

            auto g_0_xyyzzzz = contr_buffer_skxx[sk_off + 25 * ccomps * dcomps];

            auto g_0_xyzzzzz = contr_buffer_skxx[sk_off + 26 * ccomps * dcomps];

            auto g_0_xzzzzzz = contr_buffer_skxx[sk_off + 27 * ccomps * dcomps];

            auto g_0_yyyyyyy = contr_buffer_skxx[sk_off + 28 * ccomps * dcomps];

            auto g_0_yyyyyyz = contr_buffer_skxx[sk_off + 29 * ccomps * dcomps];

            auto g_0_yyyyyzz = contr_buffer_skxx[sk_off + 30 * ccomps * dcomps];

            auto g_0_yyyyzzz = contr_buffer_skxx[sk_off + 31 * ccomps * dcomps];

            auto g_0_yyyzzzz = contr_buffer_skxx[sk_off + 32 * ccomps * dcomps];

            auto g_0_yyzzzzz = contr_buffer_skxx[sk_off + 33 * ccomps * dcomps];

            auto g_0_yzzzzzz = contr_buffer_skxx[sk_off + 34 * ccomps * dcomps];

            auto g_0_zzzzzzz = contr_buffer_skxx[sk_off + 35 * ccomps * dcomps];

            /// set up bra offset for contr_buffer_pixx

            const auto pi_off = i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : contr_buffer_pixx

            auto g_x_xxxxxx = contr_buffer_pixx[pi_off + 0 * ccomps * dcomps];

            auto g_x_xxxxxy = contr_buffer_pixx[pi_off + 1 * ccomps * dcomps];

            auto g_x_xxxxxz = contr_buffer_pixx[pi_off + 2 * ccomps * dcomps];

            auto g_x_xxxxyy = contr_buffer_pixx[pi_off + 3 * ccomps * dcomps];

            auto g_x_xxxxyz = contr_buffer_pixx[pi_off + 4 * ccomps * dcomps];

            auto g_x_xxxxzz = contr_buffer_pixx[pi_off + 5 * ccomps * dcomps];

            auto g_x_xxxyyy = contr_buffer_pixx[pi_off + 6 * ccomps * dcomps];

            auto g_x_xxxyyz = contr_buffer_pixx[pi_off + 7 * ccomps * dcomps];

            auto g_x_xxxyzz = contr_buffer_pixx[pi_off + 8 * ccomps * dcomps];

            auto g_x_xxxzzz = contr_buffer_pixx[pi_off + 9 * ccomps * dcomps];

            auto g_x_xxyyyy = contr_buffer_pixx[pi_off + 10 * ccomps * dcomps];

            auto g_x_xxyyyz = contr_buffer_pixx[pi_off + 11 * ccomps * dcomps];

            auto g_x_xxyyzz = contr_buffer_pixx[pi_off + 12 * ccomps * dcomps];

            auto g_x_xxyzzz = contr_buffer_pixx[pi_off + 13 * ccomps * dcomps];

            auto g_x_xxzzzz = contr_buffer_pixx[pi_off + 14 * ccomps * dcomps];

            auto g_x_xyyyyy = contr_buffer_pixx[pi_off + 15 * ccomps * dcomps];

            auto g_x_xyyyyz = contr_buffer_pixx[pi_off + 16 * ccomps * dcomps];

            auto g_x_xyyyzz = contr_buffer_pixx[pi_off + 17 * ccomps * dcomps];

            auto g_x_xyyzzz = contr_buffer_pixx[pi_off + 18 * ccomps * dcomps];

            auto g_x_xyzzzz = contr_buffer_pixx[pi_off + 19 * ccomps * dcomps];

            auto g_x_xzzzzz = contr_buffer_pixx[pi_off + 20 * ccomps * dcomps];

            auto g_x_yyyyyy = contr_buffer_pixx[pi_off + 21 * ccomps * dcomps];

            auto g_x_yyyyyz = contr_buffer_pixx[pi_off + 22 * ccomps * dcomps];

            auto g_x_yyyyzz = contr_buffer_pixx[pi_off + 23 * ccomps * dcomps];

            auto g_x_yyyzzz = contr_buffer_pixx[pi_off + 24 * ccomps * dcomps];

            auto g_x_yyzzzz = contr_buffer_pixx[pi_off + 25 * ccomps * dcomps];

            auto g_x_yzzzzz = contr_buffer_pixx[pi_off + 26 * ccomps * dcomps];

            auto g_x_zzzzzz = contr_buffer_pixx[pi_off + 27 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxy, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxxzz, g_0_xxxxyy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxxzzz, g_0_xxxyyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxxzzzz, g_0_xxyyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xxzzzzz, g_0_xyyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_xzzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_x_xxxxxx, g_x_xxxxxy, g_x_xxxxxz, g_x_xxxxyy, g_x_xxxxyz, g_x_xxxxzz, g_x_xxxyyy, g_x_xxxyyz, g_x_xxxyzz, g_x_xxxzzz, g_x_xxyyyy, g_x_xxyyyz, g_x_xxyyzz, g_x_xxyzzz, g_x_xxzzzz, g_x_xyyyyy, g_x_xyyyyz, g_x_xyyyzz, g_x_xyyzzz, g_x_xyzzzz, g_x_xzzzzz, g_x_yyyyyy, g_x_yyyyyz, g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz, g_x_yzzzzz, g_x_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xxxxxx[k] = -g_0_xxxxxx[k] * ab_x + g_0_xxxxxxx[k];

                g_x_xxxxxy[k] = -g_0_xxxxxy[k] * ab_x + g_0_xxxxxxy[k];

                g_x_xxxxxz[k] = -g_0_xxxxxz[k] * ab_x + g_0_xxxxxxz[k];

                g_x_xxxxyy[k] = -g_0_xxxxyy[k] * ab_x + g_0_xxxxxyy[k];

                g_x_xxxxyz[k] = -g_0_xxxxyz[k] * ab_x + g_0_xxxxxyz[k];

                g_x_xxxxzz[k] = -g_0_xxxxzz[k] * ab_x + g_0_xxxxxzz[k];

                g_x_xxxyyy[k] = -g_0_xxxyyy[k] * ab_x + g_0_xxxxyyy[k];

                g_x_xxxyyz[k] = -g_0_xxxyyz[k] * ab_x + g_0_xxxxyyz[k];

                g_x_xxxyzz[k] = -g_0_xxxyzz[k] * ab_x + g_0_xxxxyzz[k];

                g_x_xxxzzz[k] = -g_0_xxxzzz[k] * ab_x + g_0_xxxxzzz[k];

                g_x_xxyyyy[k] = -g_0_xxyyyy[k] * ab_x + g_0_xxxyyyy[k];

                g_x_xxyyyz[k] = -g_0_xxyyyz[k] * ab_x + g_0_xxxyyyz[k];

                g_x_xxyyzz[k] = -g_0_xxyyzz[k] * ab_x + g_0_xxxyyzz[k];

                g_x_xxyzzz[k] = -g_0_xxyzzz[k] * ab_x + g_0_xxxyzzz[k];

                g_x_xxzzzz[k] = -g_0_xxzzzz[k] * ab_x + g_0_xxxzzzz[k];

                g_x_xyyyyy[k] = -g_0_xyyyyy[k] * ab_x + g_0_xxyyyyy[k];

                g_x_xyyyyz[k] = -g_0_xyyyyz[k] * ab_x + g_0_xxyyyyz[k];

                g_x_xyyyzz[k] = -g_0_xyyyzz[k] * ab_x + g_0_xxyyyzz[k];

                g_x_xyyzzz[k] = -g_0_xyyzzz[k] * ab_x + g_0_xxyyzzz[k];

                g_x_xyzzzz[k] = -g_0_xyzzzz[k] * ab_x + g_0_xxyzzzz[k];

                g_x_xzzzzz[k] = -g_0_xzzzzz[k] * ab_x + g_0_xxzzzzz[k];

                g_x_yyyyyy[k] = -g_0_yyyyyy[k] * ab_x + g_0_xyyyyyy[k];

                g_x_yyyyyz[k] = -g_0_yyyyyz[k] * ab_x + g_0_xyyyyyz[k];

                g_x_yyyyzz[k] = -g_0_yyyyzz[k] * ab_x + g_0_xyyyyzz[k];

                g_x_yyyzzz[k] = -g_0_yyyzzz[k] * ab_x + g_0_xyyyzzz[k];

                g_x_yyzzzz[k] = -g_0_yyzzzz[k] * ab_x + g_0_xyyzzzz[k];

                g_x_yzzzzz[k] = -g_0_yzzzzz[k] * ab_x + g_0_xyzzzzz[k];

                g_x_zzzzzz[k] = -g_0_zzzzzz[k] * ab_x + g_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : contr_buffer_pixx

            auto g_y_xxxxxx = contr_buffer_pixx[pi_off + 28 * ccomps * dcomps];

            auto g_y_xxxxxy = contr_buffer_pixx[pi_off + 29 * ccomps * dcomps];

            auto g_y_xxxxxz = contr_buffer_pixx[pi_off + 30 * ccomps * dcomps];

            auto g_y_xxxxyy = contr_buffer_pixx[pi_off + 31 * ccomps * dcomps];

            auto g_y_xxxxyz = contr_buffer_pixx[pi_off + 32 * ccomps * dcomps];

            auto g_y_xxxxzz = contr_buffer_pixx[pi_off + 33 * ccomps * dcomps];

            auto g_y_xxxyyy = contr_buffer_pixx[pi_off + 34 * ccomps * dcomps];

            auto g_y_xxxyyz = contr_buffer_pixx[pi_off + 35 * ccomps * dcomps];

            auto g_y_xxxyzz = contr_buffer_pixx[pi_off + 36 * ccomps * dcomps];

            auto g_y_xxxzzz = contr_buffer_pixx[pi_off + 37 * ccomps * dcomps];

            auto g_y_xxyyyy = contr_buffer_pixx[pi_off + 38 * ccomps * dcomps];

            auto g_y_xxyyyz = contr_buffer_pixx[pi_off + 39 * ccomps * dcomps];

            auto g_y_xxyyzz = contr_buffer_pixx[pi_off + 40 * ccomps * dcomps];

            auto g_y_xxyzzz = contr_buffer_pixx[pi_off + 41 * ccomps * dcomps];

            auto g_y_xxzzzz = contr_buffer_pixx[pi_off + 42 * ccomps * dcomps];

            auto g_y_xyyyyy = contr_buffer_pixx[pi_off + 43 * ccomps * dcomps];

            auto g_y_xyyyyz = contr_buffer_pixx[pi_off + 44 * ccomps * dcomps];

            auto g_y_xyyyzz = contr_buffer_pixx[pi_off + 45 * ccomps * dcomps];

            auto g_y_xyyzzz = contr_buffer_pixx[pi_off + 46 * ccomps * dcomps];

            auto g_y_xyzzzz = contr_buffer_pixx[pi_off + 47 * ccomps * dcomps];

            auto g_y_xzzzzz = contr_buffer_pixx[pi_off + 48 * ccomps * dcomps];

            auto g_y_yyyyyy = contr_buffer_pixx[pi_off + 49 * ccomps * dcomps];

            auto g_y_yyyyyz = contr_buffer_pixx[pi_off + 50 * ccomps * dcomps];

            auto g_y_yyyyzz = contr_buffer_pixx[pi_off + 51 * ccomps * dcomps];

            auto g_y_yyyzzz = contr_buffer_pixx[pi_off + 52 * ccomps * dcomps];

            auto g_y_yyzzzz = contr_buffer_pixx[pi_off + 53 * ccomps * dcomps];

            auto g_y_yzzzzz = contr_buffer_pixx[pi_off + 54 * ccomps * dcomps];

            auto g_y_zzzzzz = contr_buffer_pixx[pi_off + 55 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxxy, g_0_xxxxxy, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyz, g_0_yyyyyzz, g_0_yyyyzz, g_0_yyyyzzz, g_0_yyyzzz, g_0_yyyzzzz, g_0_yyzzzz, g_0_yyzzzzz, g_0_yzzzzz, g_0_yzzzzzz, g_0_zzzzzz, g_y_xxxxxx, g_y_xxxxxy, g_y_xxxxxz, g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxzz, g_y_xxxyyy, g_y_xxxyyz, g_y_xxxyzz, g_y_xxxzzz, g_y_xxyyyy, g_y_xxyyyz, g_y_xxyyzz, g_y_xxyzzz, g_y_xxzzzz, g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyzz, g_y_xyyzzz, g_y_xyzzzz, g_y_xzzzzz, g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyzz, g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz, g_y_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xxxxxx[k] = -g_0_xxxxxx[k] * ab_y + g_0_xxxxxxy[k];

                g_y_xxxxxy[k] = -g_0_xxxxxy[k] * ab_y + g_0_xxxxxyy[k];

                g_y_xxxxxz[k] = -g_0_xxxxxz[k] * ab_y + g_0_xxxxxyz[k];

                g_y_xxxxyy[k] = -g_0_xxxxyy[k] * ab_y + g_0_xxxxyyy[k];

                g_y_xxxxyz[k] = -g_0_xxxxyz[k] * ab_y + g_0_xxxxyyz[k];

                g_y_xxxxzz[k] = -g_0_xxxxzz[k] * ab_y + g_0_xxxxyzz[k];

                g_y_xxxyyy[k] = -g_0_xxxyyy[k] * ab_y + g_0_xxxyyyy[k];

                g_y_xxxyyz[k] = -g_0_xxxyyz[k] * ab_y + g_0_xxxyyyz[k];

                g_y_xxxyzz[k] = -g_0_xxxyzz[k] * ab_y + g_0_xxxyyzz[k];

                g_y_xxxzzz[k] = -g_0_xxxzzz[k] * ab_y + g_0_xxxyzzz[k];

                g_y_xxyyyy[k] = -g_0_xxyyyy[k] * ab_y + g_0_xxyyyyy[k];

                g_y_xxyyyz[k] = -g_0_xxyyyz[k] * ab_y + g_0_xxyyyyz[k];

                g_y_xxyyzz[k] = -g_0_xxyyzz[k] * ab_y + g_0_xxyyyzz[k];

                g_y_xxyzzz[k] = -g_0_xxyzzz[k] * ab_y + g_0_xxyyzzz[k];

                g_y_xxzzzz[k] = -g_0_xxzzzz[k] * ab_y + g_0_xxyzzzz[k];

                g_y_xyyyyy[k] = -g_0_xyyyyy[k] * ab_y + g_0_xyyyyyy[k];

                g_y_xyyyyz[k] = -g_0_xyyyyz[k] * ab_y + g_0_xyyyyyz[k];

                g_y_xyyyzz[k] = -g_0_xyyyzz[k] * ab_y + g_0_xyyyyzz[k];

                g_y_xyyzzz[k] = -g_0_xyyzzz[k] * ab_y + g_0_xyyyzzz[k];

                g_y_xyzzzz[k] = -g_0_xyzzzz[k] * ab_y + g_0_xyyzzzz[k];

                g_y_xzzzzz[k] = -g_0_xzzzzz[k] * ab_y + g_0_xyzzzzz[k];

                g_y_yyyyyy[k] = -g_0_yyyyyy[k] * ab_y + g_0_yyyyyyy[k];

                g_y_yyyyyz[k] = -g_0_yyyyyz[k] * ab_y + g_0_yyyyyyz[k];

                g_y_yyyyzz[k] = -g_0_yyyyzz[k] * ab_y + g_0_yyyyyzz[k];

                g_y_yyyzzz[k] = -g_0_yyyzzz[k] * ab_y + g_0_yyyyzzz[k];

                g_y_yyzzzz[k] = -g_0_yyzzzz[k] * ab_y + g_0_yyyzzzz[k];

                g_y_yzzzzz[k] = -g_0_yzzzzz[k] * ab_y + g_0_yyzzzzz[k];

                g_y_zzzzzz[k] = -g_0_zzzzzz[k] * ab_y + g_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : contr_buffer_pixx

            auto g_z_xxxxxx = contr_buffer_pixx[pi_off + 56 * ccomps * dcomps];

            auto g_z_xxxxxy = contr_buffer_pixx[pi_off + 57 * ccomps * dcomps];

            auto g_z_xxxxxz = contr_buffer_pixx[pi_off + 58 * ccomps * dcomps];

            auto g_z_xxxxyy = contr_buffer_pixx[pi_off + 59 * ccomps * dcomps];

            auto g_z_xxxxyz = contr_buffer_pixx[pi_off + 60 * ccomps * dcomps];

            auto g_z_xxxxzz = contr_buffer_pixx[pi_off + 61 * ccomps * dcomps];

            auto g_z_xxxyyy = contr_buffer_pixx[pi_off + 62 * ccomps * dcomps];

            auto g_z_xxxyyz = contr_buffer_pixx[pi_off + 63 * ccomps * dcomps];

            auto g_z_xxxyzz = contr_buffer_pixx[pi_off + 64 * ccomps * dcomps];

            auto g_z_xxxzzz = contr_buffer_pixx[pi_off + 65 * ccomps * dcomps];

            auto g_z_xxyyyy = contr_buffer_pixx[pi_off + 66 * ccomps * dcomps];

            auto g_z_xxyyyz = contr_buffer_pixx[pi_off + 67 * ccomps * dcomps];

            auto g_z_xxyyzz = contr_buffer_pixx[pi_off + 68 * ccomps * dcomps];

            auto g_z_xxyzzz = contr_buffer_pixx[pi_off + 69 * ccomps * dcomps];

            auto g_z_xxzzzz = contr_buffer_pixx[pi_off + 70 * ccomps * dcomps];

            auto g_z_xyyyyy = contr_buffer_pixx[pi_off + 71 * ccomps * dcomps];

            auto g_z_xyyyyz = contr_buffer_pixx[pi_off + 72 * ccomps * dcomps];

            auto g_z_xyyyzz = contr_buffer_pixx[pi_off + 73 * ccomps * dcomps];

            auto g_z_xyyzzz = contr_buffer_pixx[pi_off + 74 * ccomps * dcomps];

            auto g_z_xyzzzz = contr_buffer_pixx[pi_off + 75 * ccomps * dcomps];

            auto g_z_xzzzzz = contr_buffer_pixx[pi_off + 76 * ccomps * dcomps];

            auto g_z_yyyyyy = contr_buffer_pixx[pi_off + 77 * ccomps * dcomps];

            auto g_z_yyyyyz = contr_buffer_pixx[pi_off + 78 * ccomps * dcomps];

            auto g_z_yyyyzz = contr_buffer_pixx[pi_off + 79 * ccomps * dcomps];

            auto g_z_yyyzzz = contr_buffer_pixx[pi_off + 80 * ccomps * dcomps];

            auto g_z_yyzzzz = contr_buffer_pixx[pi_off + 81 * ccomps * dcomps];

            auto g_z_yzzzzz = contr_buffer_pixx[pi_off + 82 * ccomps * dcomps];

            auto g_z_zzzzzz = contr_buffer_pixx[pi_off + 83 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxxz, g_0_xxxxxy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxxzz, g_0_xxxxyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxxzzz, g_0_xxxyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxxzzzz, g_0_xxyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xxzzzzz, g_0_xyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_xzzzzzz, g_0_yyyyyy, g_0_yyyyyyz, g_0_yyyyyz, g_0_yyyyyzz, g_0_yyyyzz, g_0_yyyyzzz, g_0_yyyzzz, g_0_yyyzzzz, g_0_yyzzzz, g_0_yyzzzzz, g_0_yzzzzz, g_0_yzzzzzz, g_0_zzzzzz, g_0_zzzzzzz, g_z_xxxxxx, g_z_xxxxxy, g_z_xxxxxz, g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxzz, g_z_xxxyyy, g_z_xxxyyz, g_z_xxxyzz, g_z_xxxzzz, g_z_xxyyyy, g_z_xxyyyz, g_z_xxyyzz, g_z_xxyzzz, g_z_xxzzzz, g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyzz, g_z_xyyzzz, g_z_xyzzzz, g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz, g_z_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xxxxxx[k] = -g_0_xxxxxx[k] * ab_z + g_0_xxxxxxz[k];

                g_z_xxxxxy[k] = -g_0_xxxxxy[k] * ab_z + g_0_xxxxxyz[k];

                g_z_xxxxxz[k] = -g_0_xxxxxz[k] * ab_z + g_0_xxxxxzz[k];

                g_z_xxxxyy[k] = -g_0_xxxxyy[k] * ab_z + g_0_xxxxyyz[k];

                g_z_xxxxyz[k] = -g_0_xxxxyz[k] * ab_z + g_0_xxxxyzz[k];

                g_z_xxxxzz[k] = -g_0_xxxxzz[k] * ab_z + g_0_xxxxzzz[k];

                g_z_xxxyyy[k] = -g_0_xxxyyy[k] * ab_z + g_0_xxxyyyz[k];

                g_z_xxxyyz[k] = -g_0_xxxyyz[k] * ab_z + g_0_xxxyyzz[k];

                g_z_xxxyzz[k] = -g_0_xxxyzz[k] * ab_z + g_0_xxxyzzz[k];

                g_z_xxxzzz[k] = -g_0_xxxzzz[k] * ab_z + g_0_xxxzzzz[k];

                g_z_xxyyyy[k] = -g_0_xxyyyy[k] * ab_z + g_0_xxyyyyz[k];

                g_z_xxyyyz[k] = -g_0_xxyyyz[k] * ab_z + g_0_xxyyyzz[k];

                g_z_xxyyzz[k] = -g_0_xxyyzz[k] * ab_z + g_0_xxyyzzz[k];

                g_z_xxyzzz[k] = -g_0_xxyzzz[k] * ab_z + g_0_xxyzzzz[k];

                g_z_xxzzzz[k] = -g_0_xxzzzz[k] * ab_z + g_0_xxzzzzz[k];

                g_z_xyyyyy[k] = -g_0_xyyyyy[k] * ab_z + g_0_xyyyyyz[k];

                g_z_xyyyyz[k] = -g_0_xyyyyz[k] * ab_z + g_0_xyyyyzz[k];

                g_z_xyyyzz[k] = -g_0_xyyyzz[k] * ab_z + g_0_xyyyzzz[k];

                g_z_xyyzzz[k] = -g_0_xyyzzz[k] * ab_z + g_0_xyyzzzz[k];

                g_z_xyzzzz[k] = -g_0_xyzzzz[k] * ab_z + g_0_xyzzzzz[k];

                g_z_xzzzzz[k] = -g_0_xzzzzz[k] * ab_z + g_0_xzzzzzz[k];

                g_z_yyyyyy[k] = -g_0_yyyyyy[k] * ab_z + g_0_yyyyyyz[k];

                g_z_yyyyyz[k] = -g_0_yyyyyz[k] * ab_z + g_0_yyyyyzz[k];

                g_z_yyyyzz[k] = -g_0_yyyyzz[k] * ab_z + g_0_yyyyzzz[k];

                g_z_yyyzzz[k] = -g_0_yyyzzz[k] * ab_z + g_0_yyyzzzz[k];

                g_z_yyzzzz[k] = -g_0_yyzzzz[k] * ab_z + g_0_yyzzzzz[k];

                g_z_yzzzzz[k] = -g_0_yzzzzz[k] * ab_z + g_0_yzzzzzz[k];

                g_z_zzzzzz[k] = -g_0_zzzzzz[k] * ab_z + g_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace
