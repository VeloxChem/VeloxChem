#include "ElectronRepulsionGeom0100ContrRecSKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_skxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_skxx,
                                            const size_t idx_sixx,
                                            const size_t idx_slxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SISS

            const auto si_off = idx_sixx + i * dcomps + j;

            auto g_0_xxxxxx = cbuffer.data(si_off + 0 * ccomps * dcomps);

            auto g_0_xxxxxy = cbuffer.data(si_off + 1 * ccomps * dcomps);

            auto g_0_xxxxxz = cbuffer.data(si_off + 2 * ccomps * dcomps);

            auto g_0_xxxxyy = cbuffer.data(si_off + 3 * ccomps * dcomps);

            auto g_0_xxxxyz = cbuffer.data(si_off + 4 * ccomps * dcomps);

            auto g_0_xxxxzz = cbuffer.data(si_off + 5 * ccomps * dcomps);

            auto g_0_xxxyyy = cbuffer.data(si_off + 6 * ccomps * dcomps);

            auto g_0_xxxyyz = cbuffer.data(si_off + 7 * ccomps * dcomps);

            auto g_0_xxxyzz = cbuffer.data(si_off + 8 * ccomps * dcomps);

            auto g_0_xxxzzz = cbuffer.data(si_off + 9 * ccomps * dcomps);

            auto g_0_xxyyyy = cbuffer.data(si_off + 10 * ccomps * dcomps);

            auto g_0_xxyyyz = cbuffer.data(si_off + 11 * ccomps * dcomps);

            auto g_0_xxyyzz = cbuffer.data(si_off + 12 * ccomps * dcomps);

            auto g_0_xxyzzz = cbuffer.data(si_off + 13 * ccomps * dcomps);

            auto g_0_xxzzzz = cbuffer.data(si_off + 14 * ccomps * dcomps);

            auto g_0_xyyyyy = cbuffer.data(si_off + 15 * ccomps * dcomps);

            auto g_0_xyyyyz = cbuffer.data(si_off + 16 * ccomps * dcomps);

            auto g_0_xyyyzz = cbuffer.data(si_off + 17 * ccomps * dcomps);

            auto g_0_xyyzzz = cbuffer.data(si_off + 18 * ccomps * dcomps);

            auto g_0_xyzzzz = cbuffer.data(si_off + 19 * ccomps * dcomps);

            auto g_0_xzzzzz = cbuffer.data(si_off + 20 * ccomps * dcomps);

            auto g_0_yyyyyy = cbuffer.data(si_off + 21 * ccomps * dcomps);

            auto g_0_yyyyyz = cbuffer.data(si_off + 22 * ccomps * dcomps);

            auto g_0_yyyyzz = cbuffer.data(si_off + 23 * ccomps * dcomps);

            auto g_0_yyyzzz = cbuffer.data(si_off + 24 * ccomps * dcomps);

            auto g_0_yyzzzz = cbuffer.data(si_off + 25 * ccomps * dcomps);

            auto g_0_yzzzzz = cbuffer.data(si_off + 26 * ccomps * dcomps);

            auto g_0_zzzzzz = cbuffer.data(si_off + 27 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SLSS

            const auto sl_off = idx_slxx + i * dcomps + j;

            auto g_0_xxxxxxxx = cbuffer.data(sl_off + 0 * ccomps * dcomps);

            auto g_0_xxxxxxxy = cbuffer.data(sl_off + 1 * ccomps * dcomps);

            auto g_0_xxxxxxxz = cbuffer.data(sl_off + 2 * ccomps * dcomps);

            auto g_0_xxxxxxyy = cbuffer.data(sl_off + 3 * ccomps * dcomps);

            auto g_0_xxxxxxyz = cbuffer.data(sl_off + 4 * ccomps * dcomps);

            auto g_0_xxxxxxzz = cbuffer.data(sl_off + 5 * ccomps * dcomps);

            auto g_0_xxxxxyyy = cbuffer.data(sl_off + 6 * ccomps * dcomps);

            auto g_0_xxxxxyyz = cbuffer.data(sl_off + 7 * ccomps * dcomps);

            auto g_0_xxxxxyzz = cbuffer.data(sl_off + 8 * ccomps * dcomps);

            auto g_0_xxxxxzzz = cbuffer.data(sl_off + 9 * ccomps * dcomps);

            auto g_0_xxxxyyyy = cbuffer.data(sl_off + 10 * ccomps * dcomps);

            auto g_0_xxxxyyyz = cbuffer.data(sl_off + 11 * ccomps * dcomps);

            auto g_0_xxxxyyzz = cbuffer.data(sl_off + 12 * ccomps * dcomps);

            auto g_0_xxxxyzzz = cbuffer.data(sl_off + 13 * ccomps * dcomps);

            auto g_0_xxxxzzzz = cbuffer.data(sl_off + 14 * ccomps * dcomps);

            auto g_0_xxxyyyyy = cbuffer.data(sl_off + 15 * ccomps * dcomps);

            auto g_0_xxxyyyyz = cbuffer.data(sl_off + 16 * ccomps * dcomps);

            auto g_0_xxxyyyzz = cbuffer.data(sl_off + 17 * ccomps * dcomps);

            auto g_0_xxxyyzzz = cbuffer.data(sl_off + 18 * ccomps * dcomps);

            auto g_0_xxxyzzzz = cbuffer.data(sl_off + 19 * ccomps * dcomps);

            auto g_0_xxxzzzzz = cbuffer.data(sl_off + 20 * ccomps * dcomps);

            auto g_0_xxyyyyyy = cbuffer.data(sl_off + 21 * ccomps * dcomps);

            auto g_0_xxyyyyyz = cbuffer.data(sl_off + 22 * ccomps * dcomps);

            auto g_0_xxyyyyzz = cbuffer.data(sl_off + 23 * ccomps * dcomps);

            auto g_0_xxyyyzzz = cbuffer.data(sl_off + 24 * ccomps * dcomps);

            auto g_0_xxyyzzzz = cbuffer.data(sl_off + 25 * ccomps * dcomps);

            auto g_0_xxyzzzzz = cbuffer.data(sl_off + 26 * ccomps * dcomps);

            auto g_0_xxzzzzzz = cbuffer.data(sl_off + 27 * ccomps * dcomps);

            auto g_0_xyyyyyyy = cbuffer.data(sl_off + 28 * ccomps * dcomps);

            auto g_0_xyyyyyyz = cbuffer.data(sl_off + 29 * ccomps * dcomps);

            auto g_0_xyyyyyzz = cbuffer.data(sl_off + 30 * ccomps * dcomps);

            auto g_0_xyyyyzzz = cbuffer.data(sl_off + 31 * ccomps * dcomps);

            auto g_0_xyyyzzzz = cbuffer.data(sl_off + 32 * ccomps * dcomps);

            auto g_0_xyyzzzzz = cbuffer.data(sl_off + 33 * ccomps * dcomps);

            auto g_0_xyzzzzzz = cbuffer.data(sl_off + 34 * ccomps * dcomps);

            auto g_0_xzzzzzzz = cbuffer.data(sl_off + 35 * ccomps * dcomps);

            auto g_0_yyyyyyyy = cbuffer.data(sl_off + 36 * ccomps * dcomps);

            auto g_0_yyyyyyyz = cbuffer.data(sl_off + 37 * ccomps * dcomps);

            auto g_0_yyyyyyzz = cbuffer.data(sl_off + 38 * ccomps * dcomps);

            auto g_0_yyyyyzzz = cbuffer.data(sl_off + 39 * ccomps * dcomps);

            auto g_0_yyyyzzzz = cbuffer.data(sl_off + 40 * ccomps * dcomps);

            auto g_0_yyyzzzzz = cbuffer.data(sl_off + 41 * ccomps * dcomps);

            auto g_0_yyzzzzzz = cbuffer.data(sl_off + 42 * ccomps * dcomps);

            auto g_0_yzzzzzzz = cbuffer.data(sl_off + 43 * ccomps * dcomps);

            auto g_0_zzzzzzzz = cbuffer.data(sl_off + 44 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_skxx

            const auto sk_geom_01_off = idx_geom_01_skxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_xxxxxxx = cbuffer.data(sk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxy = cbuffer.data(sk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxz = cbuffer.data(sk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxxxxyy = cbuffer.data(sk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxxxxyz = cbuffer.data(sk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxxxxzz = cbuffer.data(sk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xxxxyyy = cbuffer.data(sk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xxxxyyz = cbuffer.data(sk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xxxxyzz = cbuffer.data(sk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xxxxzzz = cbuffer.data(sk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_xxxyyyy = cbuffer.data(sk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_xxxyyyz = cbuffer.data(sk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_xxxyyzz = cbuffer.data(sk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_xxxyzzz = cbuffer.data(sk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_xxxzzzz = cbuffer.data(sk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_0_xxyyyyy = cbuffer.data(sk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_0_xxyyyyz = cbuffer.data(sk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_0_xxyyyzz = cbuffer.data(sk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_0_xxyyzzz = cbuffer.data(sk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_0_xxyzzzz = cbuffer.data(sk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_0_xxzzzzz = cbuffer.data(sk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_0_xyyyyyy = cbuffer.data(sk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_0_xyyyyyz = cbuffer.data(sk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_0_xyyyyzz = cbuffer.data(sk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_0_xyyyzzz = cbuffer.data(sk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_0_xyyzzzz = cbuffer.data(sk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_0_xyzzzzz = cbuffer.data(sk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_0_xzzzzzz = cbuffer.data(sk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_0_yyyyyyy = cbuffer.data(sk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_0_yyyyyyz = cbuffer.data(sk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_0_yyyyyzz = cbuffer.data(sk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_0_yyyyzzz = cbuffer.data(sk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_0_yyyzzzz = cbuffer.data(sk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_0_yyzzzzz = cbuffer.data(sk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_0_yzzzzzz = cbuffer.data(sk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_0_zzzzzzz = cbuffer.data(sk_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzzz, g_0_xxxxxx, g_0_xxxxxxxx, g_0_xxxxxxxy, g_0_xxxxxxxz, g_0_xxxxxxyy, g_0_xxxxxxyz, g_0_xxxxxxzz, g_0_xxxxxy, g_0_xxxxxyyy, g_0_xxxxxyyz, g_0_xxxxxyzz, g_0_xxxxxz, g_0_xxxxxzzz, g_0_xxxxyy, g_0_xxxxyyyy, g_0_xxxxyyyz, g_0_xxxxyyzz, g_0_xxxxyz, g_0_xxxxyzzz, g_0_xxxxzz, g_0_xxxxzzzz, g_0_xxxyyy, g_0_xxxyyyyy, g_0_xxxyyyyz, g_0_xxxyyyzz, g_0_xxxyyz, g_0_xxxyyzzz, g_0_xxxyzz, g_0_xxxyzzzz, g_0_xxxzzz, g_0_xxxzzzzz, g_0_xxyyyy, g_0_xxyyyyyy, g_0_xxyyyyyz, g_0_xxyyyyzz, g_0_xxyyyz, g_0_xxyyyzzz, g_0_xxyyzz, g_0_xxyyzzzz, g_0_xxyzzz, g_0_xxyzzzzz, g_0_xxzzzz, g_0_xxzzzzzz, g_0_xyyyyy, g_0_xyyyyyyy, g_0_xyyyyyyz, g_0_xyyyyyzz, g_0_xyyyyz, g_0_xyyyyzzz, g_0_xyyyzz, g_0_xyyyzzzz, g_0_xyyzzz, g_0_xyyzzzzz, g_0_xyzzzz, g_0_xyzzzzzz, g_0_xzzzzz, g_0_xzzzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xxxxxxx[k] = -7.0 * g_0_xxxxxx[k] + g_0_xxxxxxxx[k];

                g_0_x_0_xxxxxxy[k] = -6.0 * g_0_xxxxxy[k] + g_0_xxxxxxxy[k];

                g_0_x_0_xxxxxxz[k] = -6.0 * g_0_xxxxxz[k] + g_0_xxxxxxxz[k];

                g_0_x_0_xxxxxyy[k] = -5.0 * g_0_xxxxyy[k] + g_0_xxxxxxyy[k];

                g_0_x_0_xxxxxyz[k] = -5.0 * g_0_xxxxyz[k] + g_0_xxxxxxyz[k];

                g_0_x_0_xxxxxzz[k] = -5.0 * g_0_xxxxzz[k] + g_0_xxxxxxzz[k];

                g_0_x_0_xxxxyyy[k] = -4.0 * g_0_xxxyyy[k] + g_0_xxxxxyyy[k];

                g_0_x_0_xxxxyyz[k] = -4.0 * g_0_xxxyyz[k] + g_0_xxxxxyyz[k];

                g_0_x_0_xxxxyzz[k] = -4.0 * g_0_xxxyzz[k] + g_0_xxxxxyzz[k];

                g_0_x_0_xxxxzzz[k] = -4.0 * g_0_xxxzzz[k] + g_0_xxxxxzzz[k];

                g_0_x_0_xxxyyyy[k] = -3.0 * g_0_xxyyyy[k] + g_0_xxxxyyyy[k];

                g_0_x_0_xxxyyyz[k] = -3.0 * g_0_xxyyyz[k] + g_0_xxxxyyyz[k];

                g_0_x_0_xxxyyzz[k] = -3.0 * g_0_xxyyzz[k] + g_0_xxxxyyzz[k];

                g_0_x_0_xxxyzzz[k] = -3.0 * g_0_xxyzzz[k] + g_0_xxxxyzzz[k];

                g_0_x_0_xxxzzzz[k] = -3.0 * g_0_xxzzzz[k] + g_0_xxxxzzzz[k];

                g_0_x_0_xxyyyyy[k] = -2.0 * g_0_xyyyyy[k] + g_0_xxxyyyyy[k];

                g_0_x_0_xxyyyyz[k] = -2.0 * g_0_xyyyyz[k] + g_0_xxxyyyyz[k];

                g_0_x_0_xxyyyzz[k] = -2.0 * g_0_xyyyzz[k] + g_0_xxxyyyzz[k];

                g_0_x_0_xxyyzzz[k] = -2.0 * g_0_xyyzzz[k] + g_0_xxxyyzzz[k];

                g_0_x_0_xxyzzzz[k] = -2.0 * g_0_xyzzzz[k] + g_0_xxxyzzzz[k];

                g_0_x_0_xxzzzzz[k] = -2.0 * g_0_xzzzzz[k] + g_0_xxxzzzzz[k];

                g_0_x_0_xyyyyyy[k] = -g_0_yyyyyy[k] + g_0_xxyyyyyy[k];

                g_0_x_0_xyyyyyz[k] = -g_0_yyyyyz[k] + g_0_xxyyyyyz[k];

                g_0_x_0_xyyyyzz[k] = -g_0_yyyyzz[k] + g_0_xxyyyyzz[k];

                g_0_x_0_xyyyzzz[k] = -g_0_yyyzzz[k] + g_0_xxyyyzzz[k];

                g_0_x_0_xyyzzzz[k] = -g_0_yyzzzz[k] + g_0_xxyyzzzz[k];

                g_0_x_0_xyzzzzz[k] = -g_0_yzzzzz[k] + g_0_xxyzzzzz[k];

                g_0_x_0_xzzzzzz[k] = -g_0_zzzzzz[k] + g_0_xxzzzzzz[k];

                g_0_x_0_yyyyyyy[k] = g_0_xyyyyyyy[k];

                g_0_x_0_yyyyyyz[k] = g_0_xyyyyyyz[k];

                g_0_x_0_yyyyyzz[k] = g_0_xyyyyyzz[k];

                g_0_x_0_yyyyzzz[k] = g_0_xyyyyzzz[k];

                g_0_x_0_yyyzzzz[k] = g_0_xyyyzzzz[k];

                g_0_x_0_yyzzzzz[k] = g_0_xyyzzzzz[k];

                g_0_x_0_yzzzzzz[k] = g_0_xyzzzzzz[k];

                g_0_x_0_zzzzzzz[k] = g_0_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_xxxxxxx = cbuffer.data(sk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxy = cbuffer.data(sk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxz = cbuffer.data(sk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_0_xxxxxyy = cbuffer.data(sk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_0_xxxxxyz = cbuffer.data(sk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_0_xxxxxzz = cbuffer.data(sk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_0_xxxxyyy = cbuffer.data(sk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_0_xxxxyyz = cbuffer.data(sk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_0_xxxxyzz = cbuffer.data(sk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_0_xxxxzzz = cbuffer.data(sk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_0_xxxyyyy = cbuffer.data(sk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_0_xxxyyyz = cbuffer.data(sk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_0_xxxyyzz = cbuffer.data(sk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_0_xxxyzzz = cbuffer.data(sk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_0_xxxzzzz = cbuffer.data(sk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_0_xxyyyyy = cbuffer.data(sk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_0_xxyyyyz = cbuffer.data(sk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_0_xxyyyzz = cbuffer.data(sk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_0_xxyyzzz = cbuffer.data(sk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_0_xxyzzzz = cbuffer.data(sk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_0_xxzzzzz = cbuffer.data(sk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_0_xyyyyyy = cbuffer.data(sk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_0_xyyyyyz = cbuffer.data(sk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_0_xyyyyzz = cbuffer.data(sk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_0_xyyyzzz = cbuffer.data(sk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_0_xyyzzzz = cbuffer.data(sk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_0_xyzzzzz = cbuffer.data(sk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_0_xzzzzzz = cbuffer.data(sk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_0_yyyyyyy = cbuffer.data(sk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_0_yyyyyyz = cbuffer.data(sk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_0_yyyyyzz = cbuffer.data(sk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_0_yyyyzzz = cbuffer.data(sk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_0_yyyzzzz = cbuffer.data(sk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_0_yyzzzzz = cbuffer.data(sk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_0_yzzzzzz = cbuffer.data(sk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_0_zzzzzzz = cbuffer.data(sk_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxxxy, g_0_xxxxxxyy, g_0_xxxxxxyz, g_0_xxxxxy, g_0_xxxxxyyy, g_0_xxxxxyyz, g_0_xxxxxyzz, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyyyy, g_0_xxxxyyyz, g_0_xxxxyyzz, g_0_xxxxyz, g_0_xxxxyzzz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyyyy, g_0_xxxyyyyz, g_0_xxxyyyzz, g_0_xxxyyz, g_0_xxxyyzzz, g_0_xxxyzz, g_0_xxxyzzzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyyyy, g_0_xxyyyyyz, g_0_xxyyyyzz, g_0_xxyyyz, g_0_xxyyyzzz, g_0_xxyyzz, g_0_xxyyzzzz, g_0_xxyzzz, g_0_xxyzzzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyyyy, g_0_xyyyyyyz, g_0_xyyyyyzz, g_0_xyyyyz, g_0_xyyyyzzz, g_0_xyyyzz, g_0_xyyyzzzz, g_0_xyyzzz, g_0_xyyzzzzz, g_0_xyzzzz, g_0_xyzzzzzz, g_0_xzzzzz, g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzzz, g_0_yyyyyy, g_0_yyyyyyyy, g_0_yyyyyyyz, g_0_yyyyyyzz, g_0_yyyyyz, g_0_yyyyyzzz, g_0_yyyyzz, g_0_yyyyzzzz, g_0_yyyzzz, g_0_yyyzzzzz, g_0_yyzzzz, g_0_yyzzzzzz, g_0_yzzzzz, g_0_yzzzzzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xxxxxxx[k] = g_0_xxxxxxxy[k];

                g_0_y_0_xxxxxxy[k] = -g_0_xxxxxx[k] + g_0_xxxxxxyy[k];

                g_0_y_0_xxxxxxz[k] = g_0_xxxxxxyz[k];

                g_0_y_0_xxxxxyy[k] = -2.0 * g_0_xxxxxy[k] + g_0_xxxxxyyy[k];

                g_0_y_0_xxxxxyz[k] = -g_0_xxxxxz[k] + g_0_xxxxxyyz[k];

                g_0_y_0_xxxxxzz[k] = g_0_xxxxxyzz[k];

                g_0_y_0_xxxxyyy[k] = -3.0 * g_0_xxxxyy[k] + g_0_xxxxyyyy[k];

                g_0_y_0_xxxxyyz[k] = -2.0 * g_0_xxxxyz[k] + g_0_xxxxyyyz[k];

                g_0_y_0_xxxxyzz[k] = -g_0_xxxxzz[k] + g_0_xxxxyyzz[k];

                g_0_y_0_xxxxzzz[k] = g_0_xxxxyzzz[k];

                g_0_y_0_xxxyyyy[k] = -4.0 * g_0_xxxyyy[k] + g_0_xxxyyyyy[k];

                g_0_y_0_xxxyyyz[k] = -3.0 * g_0_xxxyyz[k] + g_0_xxxyyyyz[k];

                g_0_y_0_xxxyyzz[k] = -2.0 * g_0_xxxyzz[k] + g_0_xxxyyyzz[k];

                g_0_y_0_xxxyzzz[k] = -g_0_xxxzzz[k] + g_0_xxxyyzzz[k];

                g_0_y_0_xxxzzzz[k] = g_0_xxxyzzzz[k];

                g_0_y_0_xxyyyyy[k] = -5.0 * g_0_xxyyyy[k] + g_0_xxyyyyyy[k];

                g_0_y_0_xxyyyyz[k] = -4.0 * g_0_xxyyyz[k] + g_0_xxyyyyyz[k];

                g_0_y_0_xxyyyzz[k] = -3.0 * g_0_xxyyzz[k] + g_0_xxyyyyzz[k];

                g_0_y_0_xxyyzzz[k] = -2.0 * g_0_xxyzzz[k] + g_0_xxyyyzzz[k];

                g_0_y_0_xxyzzzz[k] = -g_0_xxzzzz[k] + g_0_xxyyzzzz[k];

                g_0_y_0_xxzzzzz[k] = g_0_xxyzzzzz[k];

                g_0_y_0_xyyyyyy[k] = -6.0 * g_0_xyyyyy[k] + g_0_xyyyyyyy[k];

                g_0_y_0_xyyyyyz[k] = -5.0 * g_0_xyyyyz[k] + g_0_xyyyyyyz[k];

                g_0_y_0_xyyyyzz[k] = -4.0 * g_0_xyyyzz[k] + g_0_xyyyyyzz[k];

                g_0_y_0_xyyyzzz[k] = -3.0 * g_0_xyyzzz[k] + g_0_xyyyyzzz[k];

                g_0_y_0_xyyzzzz[k] = -2.0 * g_0_xyzzzz[k] + g_0_xyyyzzzz[k];

                g_0_y_0_xyzzzzz[k] = -g_0_xzzzzz[k] + g_0_xyyzzzzz[k];

                g_0_y_0_xzzzzzz[k] = g_0_xyzzzzzz[k];

                g_0_y_0_yyyyyyy[k] = -7.0 * g_0_yyyyyy[k] + g_0_yyyyyyyy[k];

                g_0_y_0_yyyyyyz[k] = -6.0 * g_0_yyyyyz[k] + g_0_yyyyyyyz[k];

                g_0_y_0_yyyyyzz[k] = -5.0 * g_0_yyyyzz[k] + g_0_yyyyyyzz[k];

                g_0_y_0_yyyyzzz[k] = -4.0 * g_0_yyyzzz[k] + g_0_yyyyyzzz[k];

                g_0_y_0_yyyzzzz[k] = -3.0 * g_0_yyzzzz[k] + g_0_yyyyzzzz[k];

                g_0_y_0_yyzzzzz[k] = -2.0 * g_0_yzzzzz[k] + g_0_yyyzzzzz[k];

                g_0_y_0_yzzzzzz[k] = -g_0_zzzzzz[k] + g_0_yyzzzzzz[k];

                g_0_y_0_zzzzzzz[k] = g_0_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_xxxxxxx = cbuffer.data(sk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxy = cbuffer.data(sk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxz = cbuffer.data(sk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_0_xxxxxyy = cbuffer.data(sk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_0_xxxxxyz = cbuffer.data(sk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_0_xxxxxzz = cbuffer.data(sk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_0_xxxxyyy = cbuffer.data(sk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_0_xxxxyyz = cbuffer.data(sk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_0_xxxxyzz = cbuffer.data(sk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_0_xxxxzzz = cbuffer.data(sk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_0_xxxyyyy = cbuffer.data(sk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_0_xxxyyyz = cbuffer.data(sk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_z_0_xxxyyzz = cbuffer.data(sk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_z_0_xxxyzzz = cbuffer.data(sk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_z_0_xxxzzzz = cbuffer.data(sk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_z_0_xxyyyyy = cbuffer.data(sk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_z_0_xxyyyyz = cbuffer.data(sk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_z_0_xxyyyzz = cbuffer.data(sk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_z_0_xxyyzzz = cbuffer.data(sk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_z_0_xxyzzzz = cbuffer.data(sk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_z_0_xxzzzzz = cbuffer.data(sk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_z_0_xyyyyyy = cbuffer.data(sk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_z_0_xyyyyyz = cbuffer.data(sk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_z_0_xyyyyzz = cbuffer.data(sk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_z_0_xyyyzzz = cbuffer.data(sk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_z_0_xyyzzzz = cbuffer.data(sk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_z_0_xyzzzzz = cbuffer.data(sk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_z_0_xzzzzzz = cbuffer.data(sk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_z_0_yyyyyyy = cbuffer.data(sk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_z_0_yyyyyyz = cbuffer.data(sk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_z_0_yyyyyzz = cbuffer.data(sk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_z_0_yyyyzzz = cbuffer.data(sk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_z_0_yyyzzzz = cbuffer.data(sk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_z_0_yyzzzzz = cbuffer.data(sk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_z_0_yzzzzzz = cbuffer.data(sk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_z_0_zzzzzzz = cbuffer.data(sk_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxxxz, g_0_xxxxxxyz, g_0_xxxxxxzz, g_0_xxxxxy, g_0_xxxxxyyz, g_0_xxxxxyzz, g_0_xxxxxz, g_0_xxxxxzzz, g_0_xxxxyy, g_0_xxxxyyyz, g_0_xxxxyyzz, g_0_xxxxyz, g_0_xxxxyzzz, g_0_xxxxzz, g_0_xxxxzzzz, g_0_xxxyyy, g_0_xxxyyyyz, g_0_xxxyyyzz, g_0_xxxyyz, g_0_xxxyyzzz, g_0_xxxyzz, g_0_xxxyzzzz, g_0_xxxzzz, g_0_xxxzzzzz, g_0_xxyyyy, g_0_xxyyyyyz, g_0_xxyyyyzz, g_0_xxyyyz, g_0_xxyyyzzz, g_0_xxyyzz, g_0_xxyyzzzz, g_0_xxyzzz, g_0_xxyzzzzz, g_0_xxzzzz, g_0_xxzzzzzz, g_0_xyyyyy, g_0_xyyyyyyz, g_0_xyyyyyzz, g_0_xyyyyz, g_0_xyyyyzzz, g_0_xyyyzz, g_0_xyyyzzzz, g_0_xyyzzz, g_0_xyyzzzzz, g_0_xyzzzz, g_0_xyzzzzzz, g_0_xzzzzz, g_0_xzzzzzzz, g_0_yyyyyy, g_0_yyyyyyyz, g_0_yyyyyyzz, g_0_yyyyyz, g_0_yyyyyzzz, g_0_yyyyzz, g_0_yyyyzzzz, g_0_yyyzzz, g_0_yyyzzzzz, g_0_yyzzzz, g_0_yyzzzzzz, g_0_yzzzzz, g_0_yzzzzzzz, g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzzz, g_0_zzzzzz, g_0_zzzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xxxxxxx[k] = g_0_xxxxxxxz[k];

                g_0_z_0_xxxxxxy[k] = g_0_xxxxxxyz[k];

                g_0_z_0_xxxxxxz[k] = -g_0_xxxxxx[k] + g_0_xxxxxxzz[k];

                g_0_z_0_xxxxxyy[k] = g_0_xxxxxyyz[k];

                g_0_z_0_xxxxxyz[k] = -g_0_xxxxxy[k] + g_0_xxxxxyzz[k];

                g_0_z_0_xxxxxzz[k] = -2.0 * g_0_xxxxxz[k] + g_0_xxxxxzzz[k];

                g_0_z_0_xxxxyyy[k] = g_0_xxxxyyyz[k];

                g_0_z_0_xxxxyyz[k] = -g_0_xxxxyy[k] + g_0_xxxxyyzz[k];

                g_0_z_0_xxxxyzz[k] = -2.0 * g_0_xxxxyz[k] + g_0_xxxxyzzz[k];

                g_0_z_0_xxxxzzz[k] = -3.0 * g_0_xxxxzz[k] + g_0_xxxxzzzz[k];

                g_0_z_0_xxxyyyy[k] = g_0_xxxyyyyz[k];

                g_0_z_0_xxxyyyz[k] = -g_0_xxxyyy[k] + g_0_xxxyyyzz[k];

                g_0_z_0_xxxyyzz[k] = -2.0 * g_0_xxxyyz[k] + g_0_xxxyyzzz[k];

                g_0_z_0_xxxyzzz[k] = -3.0 * g_0_xxxyzz[k] + g_0_xxxyzzzz[k];

                g_0_z_0_xxxzzzz[k] = -4.0 * g_0_xxxzzz[k] + g_0_xxxzzzzz[k];

                g_0_z_0_xxyyyyy[k] = g_0_xxyyyyyz[k];

                g_0_z_0_xxyyyyz[k] = -g_0_xxyyyy[k] + g_0_xxyyyyzz[k];

                g_0_z_0_xxyyyzz[k] = -2.0 * g_0_xxyyyz[k] + g_0_xxyyyzzz[k];

                g_0_z_0_xxyyzzz[k] = -3.0 * g_0_xxyyzz[k] + g_0_xxyyzzzz[k];

                g_0_z_0_xxyzzzz[k] = -4.0 * g_0_xxyzzz[k] + g_0_xxyzzzzz[k];

                g_0_z_0_xxzzzzz[k] = -5.0 * g_0_xxzzzz[k] + g_0_xxzzzzzz[k];

                g_0_z_0_xyyyyyy[k] = g_0_xyyyyyyz[k];

                g_0_z_0_xyyyyyz[k] = -g_0_xyyyyy[k] + g_0_xyyyyyzz[k];

                g_0_z_0_xyyyyzz[k] = -2.0 * g_0_xyyyyz[k] + g_0_xyyyyzzz[k];

                g_0_z_0_xyyyzzz[k] = -3.0 * g_0_xyyyzz[k] + g_0_xyyyzzzz[k];

                g_0_z_0_xyyzzzz[k] = -4.0 * g_0_xyyzzz[k] + g_0_xyyzzzzz[k];

                g_0_z_0_xyzzzzz[k] = -5.0 * g_0_xyzzzz[k] + g_0_xyzzzzzz[k];

                g_0_z_0_xzzzzzz[k] = -6.0 * g_0_xzzzzz[k] + g_0_xzzzzzzz[k];

                g_0_z_0_yyyyyyy[k] = g_0_yyyyyyyz[k];

                g_0_z_0_yyyyyyz[k] = -g_0_yyyyyy[k] + g_0_yyyyyyzz[k];

                g_0_z_0_yyyyyzz[k] = -2.0 * g_0_yyyyyz[k] + g_0_yyyyyzzz[k];

                g_0_z_0_yyyyzzz[k] = -3.0 * g_0_yyyyzz[k] + g_0_yyyyzzzz[k];

                g_0_z_0_yyyzzzz[k] = -4.0 * g_0_yyyzzz[k] + g_0_yyyzzzzz[k];

                g_0_z_0_yyzzzzz[k] = -5.0 * g_0_yyzzzz[k] + g_0_yyzzzzzz[k];

                g_0_z_0_yzzzzzz[k] = -6.0 * g_0_yzzzzz[k] + g_0_yzzzzzzz[k];

                g_0_z_0_zzzzzzz[k] = -7.0 * g_0_zzzzzz[k] + g_0_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

