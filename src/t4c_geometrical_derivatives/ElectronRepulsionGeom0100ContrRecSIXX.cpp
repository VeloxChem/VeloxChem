#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_sixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sixx,
                                            const size_t idx_shxx,
                                            const size_t idx_skxx,
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
            /// Set up components of auxilary buffer : SHSS

            const auto sh_off = idx_shxx + i * dcomps + j;

            auto g_0_xxxxx = cbuffer.data(sh_off + 0 * ccomps * dcomps);

            auto g_0_xxxxy = cbuffer.data(sh_off + 1 * ccomps * dcomps);

            auto g_0_xxxxz = cbuffer.data(sh_off + 2 * ccomps * dcomps);

            auto g_0_xxxyy = cbuffer.data(sh_off + 3 * ccomps * dcomps);

            auto g_0_xxxyz = cbuffer.data(sh_off + 4 * ccomps * dcomps);

            auto g_0_xxxzz = cbuffer.data(sh_off + 5 * ccomps * dcomps);

            auto g_0_xxyyy = cbuffer.data(sh_off + 6 * ccomps * dcomps);

            auto g_0_xxyyz = cbuffer.data(sh_off + 7 * ccomps * dcomps);

            auto g_0_xxyzz = cbuffer.data(sh_off + 8 * ccomps * dcomps);

            auto g_0_xxzzz = cbuffer.data(sh_off + 9 * ccomps * dcomps);

            auto g_0_xyyyy = cbuffer.data(sh_off + 10 * ccomps * dcomps);

            auto g_0_xyyyz = cbuffer.data(sh_off + 11 * ccomps * dcomps);

            auto g_0_xyyzz = cbuffer.data(sh_off + 12 * ccomps * dcomps);

            auto g_0_xyzzz = cbuffer.data(sh_off + 13 * ccomps * dcomps);

            auto g_0_xzzzz = cbuffer.data(sh_off + 14 * ccomps * dcomps);

            auto g_0_yyyyy = cbuffer.data(sh_off + 15 * ccomps * dcomps);

            auto g_0_yyyyz = cbuffer.data(sh_off + 16 * ccomps * dcomps);

            auto g_0_yyyzz = cbuffer.data(sh_off + 17 * ccomps * dcomps);

            auto g_0_yyzzz = cbuffer.data(sh_off + 18 * ccomps * dcomps);

            auto g_0_yzzzz = cbuffer.data(sh_off + 19 * ccomps * dcomps);

            auto g_0_zzzzz = cbuffer.data(sh_off + 20 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SKSS

            const auto sk_off = idx_skxx + i * dcomps + j;

            auto g_0_xxxxxxx = cbuffer.data(sk_off + 0 * ccomps * dcomps);

            auto g_0_xxxxxxy = cbuffer.data(sk_off + 1 * ccomps * dcomps);

            auto g_0_xxxxxxz = cbuffer.data(sk_off + 2 * ccomps * dcomps);

            auto g_0_xxxxxyy = cbuffer.data(sk_off + 3 * ccomps * dcomps);

            auto g_0_xxxxxyz = cbuffer.data(sk_off + 4 * ccomps * dcomps);

            auto g_0_xxxxxzz = cbuffer.data(sk_off + 5 * ccomps * dcomps);

            auto g_0_xxxxyyy = cbuffer.data(sk_off + 6 * ccomps * dcomps);

            auto g_0_xxxxyyz = cbuffer.data(sk_off + 7 * ccomps * dcomps);

            auto g_0_xxxxyzz = cbuffer.data(sk_off + 8 * ccomps * dcomps);

            auto g_0_xxxxzzz = cbuffer.data(sk_off + 9 * ccomps * dcomps);

            auto g_0_xxxyyyy = cbuffer.data(sk_off + 10 * ccomps * dcomps);

            auto g_0_xxxyyyz = cbuffer.data(sk_off + 11 * ccomps * dcomps);

            auto g_0_xxxyyzz = cbuffer.data(sk_off + 12 * ccomps * dcomps);

            auto g_0_xxxyzzz = cbuffer.data(sk_off + 13 * ccomps * dcomps);

            auto g_0_xxxzzzz = cbuffer.data(sk_off + 14 * ccomps * dcomps);

            auto g_0_xxyyyyy = cbuffer.data(sk_off + 15 * ccomps * dcomps);

            auto g_0_xxyyyyz = cbuffer.data(sk_off + 16 * ccomps * dcomps);

            auto g_0_xxyyyzz = cbuffer.data(sk_off + 17 * ccomps * dcomps);

            auto g_0_xxyyzzz = cbuffer.data(sk_off + 18 * ccomps * dcomps);

            auto g_0_xxyzzzz = cbuffer.data(sk_off + 19 * ccomps * dcomps);

            auto g_0_xxzzzzz = cbuffer.data(sk_off + 20 * ccomps * dcomps);

            auto g_0_xyyyyyy = cbuffer.data(sk_off + 21 * ccomps * dcomps);

            auto g_0_xyyyyyz = cbuffer.data(sk_off + 22 * ccomps * dcomps);

            auto g_0_xyyyyzz = cbuffer.data(sk_off + 23 * ccomps * dcomps);

            auto g_0_xyyyzzz = cbuffer.data(sk_off + 24 * ccomps * dcomps);

            auto g_0_xyyzzzz = cbuffer.data(sk_off + 25 * ccomps * dcomps);

            auto g_0_xyzzzzz = cbuffer.data(sk_off + 26 * ccomps * dcomps);

            auto g_0_xzzzzzz = cbuffer.data(sk_off + 27 * ccomps * dcomps);

            auto g_0_yyyyyyy = cbuffer.data(sk_off + 28 * ccomps * dcomps);

            auto g_0_yyyyyyz = cbuffer.data(sk_off + 29 * ccomps * dcomps);

            auto g_0_yyyyyzz = cbuffer.data(sk_off + 30 * ccomps * dcomps);

            auto g_0_yyyyzzz = cbuffer.data(sk_off + 31 * ccomps * dcomps);

            auto g_0_yyyzzzz = cbuffer.data(sk_off + 32 * ccomps * dcomps);

            auto g_0_yyzzzzz = cbuffer.data(sk_off + 33 * ccomps * dcomps);

            auto g_0_yzzzzzz = cbuffer.data(sk_off + 34 * ccomps * dcomps);

            auto g_0_zzzzzzz = cbuffer.data(sk_off + 35 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_sixx

            const auto si_geom_01_off = idx_geom_01_sixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_xxxxxx = cbuffer.data(si_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxxxy = cbuffer.data(si_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxxxz = cbuffer.data(si_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxxxyy = cbuffer.data(si_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxxxyz = cbuffer.data(si_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxxxzz = cbuffer.data(si_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xxxyyy = cbuffer.data(si_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xxxyyz = cbuffer.data(si_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xxxyzz = cbuffer.data(si_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xxxzzz = cbuffer.data(si_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_xxyyyy = cbuffer.data(si_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_xxyyyz = cbuffer.data(si_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_xxyyzz = cbuffer.data(si_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_xxyzzz = cbuffer.data(si_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_xxzzzz = cbuffer.data(si_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_0_xyyyyy = cbuffer.data(si_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_0_xyyyyz = cbuffer.data(si_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_0_xyyyzz = cbuffer.data(si_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_0_xyyzzz = cbuffer.data(si_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_0_xyzzzz = cbuffer.data(si_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_0_xzzzzz = cbuffer.data(si_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_0_yyyyyy = cbuffer.data(si_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_0_yyyyyz = cbuffer.data(si_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_0_yyyyzz = cbuffer.data(si_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_0_yyyzzz = cbuffer.data(si_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_0_yyzzzz = cbuffer.data(si_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_0_yzzzzz = cbuffer.data(si_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_0_zzzzzz = cbuffer.data(si_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_0_xxxxx, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxz, g_0_xxxxzzz, g_0_xxxyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyz, g_0_xxxyzzz, g_0_xxxzz, g_0_xxxzzzz, g_0_xxyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyz, g_0_xxyyzzz, g_0_xxyzz, g_0_xxyzzzz, g_0_xxzzz, g_0_xxzzzzz, g_0_xyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyz, g_0_xyyyzzz, g_0_xyyzz, g_0_xyyzzzz, g_0_xyzzz, g_0_xyzzzzz, g_0_xzzzz, g_0_xzzzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xxxxxx[k] = -6.0 * g_0_xxxxx[k] + g_0_xxxxxxx[k];

                g_0_x_0_xxxxxy[k] = -5.0 * g_0_xxxxy[k] + g_0_xxxxxxy[k];

                g_0_x_0_xxxxxz[k] = -5.0 * g_0_xxxxz[k] + g_0_xxxxxxz[k];

                g_0_x_0_xxxxyy[k] = -4.0 * g_0_xxxyy[k] + g_0_xxxxxyy[k];

                g_0_x_0_xxxxyz[k] = -4.0 * g_0_xxxyz[k] + g_0_xxxxxyz[k];

                g_0_x_0_xxxxzz[k] = -4.0 * g_0_xxxzz[k] + g_0_xxxxxzz[k];

                g_0_x_0_xxxyyy[k] = -3.0 * g_0_xxyyy[k] + g_0_xxxxyyy[k];

                g_0_x_0_xxxyyz[k] = -3.0 * g_0_xxyyz[k] + g_0_xxxxyyz[k];

                g_0_x_0_xxxyzz[k] = -3.0 * g_0_xxyzz[k] + g_0_xxxxyzz[k];

                g_0_x_0_xxxzzz[k] = -3.0 * g_0_xxzzz[k] + g_0_xxxxzzz[k];

                g_0_x_0_xxyyyy[k] = -2.0 * g_0_xyyyy[k] + g_0_xxxyyyy[k];

                g_0_x_0_xxyyyz[k] = -2.0 * g_0_xyyyz[k] + g_0_xxxyyyz[k];

                g_0_x_0_xxyyzz[k] = -2.0 * g_0_xyyzz[k] + g_0_xxxyyzz[k];

                g_0_x_0_xxyzzz[k] = -2.0 * g_0_xyzzz[k] + g_0_xxxyzzz[k];

                g_0_x_0_xxzzzz[k] = -2.0 * g_0_xzzzz[k] + g_0_xxxzzzz[k];

                g_0_x_0_xyyyyy[k] = -g_0_yyyyy[k] + g_0_xxyyyyy[k];

                g_0_x_0_xyyyyz[k] = -g_0_yyyyz[k] + g_0_xxyyyyz[k];

                g_0_x_0_xyyyzz[k] = -g_0_yyyzz[k] + g_0_xxyyyzz[k];

                g_0_x_0_xyyzzz[k] = -g_0_yyzzz[k] + g_0_xxyyzzz[k];

                g_0_x_0_xyzzzz[k] = -g_0_yzzzz[k] + g_0_xxyzzzz[k];

                g_0_x_0_xzzzzz[k] = -g_0_zzzzz[k] + g_0_xxzzzzz[k];

                g_0_x_0_yyyyyy[k] = g_0_xyyyyyy[k];

                g_0_x_0_yyyyyz[k] = g_0_xyyyyyz[k];

                g_0_x_0_yyyyzz[k] = g_0_xyyyyzz[k];

                g_0_x_0_yyyzzz[k] = g_0_xyyyzzz[k];

                g_0_x_0_yyzzzz[k] = g_0_xyyzzzz[k];

                g_0_x_0_yzzzzz[k] = g_0_xyzzzzz[k];

                g_0_x_0_zzzzzz[k] = g_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_xxxxxx = cbuffer.data(si_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_0_xxxxxy = cbuffer.data(si_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_0_xxxxxz = cbuffer.data(si_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_0_xxxxyy = cbuffer.data(si_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_0_xxxxyz = cbuffer.data(si_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_0_xxxxzz = cbuffer.data(si_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_0_xxxyyy = cbuffer.data(si_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_0_xxxyyz = cbuffer.data(si_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_0_xxxyzz = cbuffer.data(si_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_0_xxxzzz = cbuffer.data(si_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_0_xxyyyy = cbuffer.data(si_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_0_xxyyyz = cbuffer.data(si_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_0_xxyyzz = cbuffer.data(si_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_0_xxyzzz = cbuffer.data(si_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_0_xxzzzz = cbuffer.data(si_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_0_xyyyyy = cbuffer.data(si_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_0_xyyyyz = cbuffer.data(si_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_0_xyyyzz = cbuffer.data(si_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_0_xyyzzz = cbuffer.data(si_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_0_xyzzzz = cbuffer.data(si_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_0_xzzzzz = cbuffer.data(si_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_0_yyyyyy = cbuffer.data(si_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_0_yyyyyz = cbuffer.data(si_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_0_yyyyzz = cbuffer.data(si_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_0_yyyzzz = cbuffer.data(si_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_0_yyzzzz = cbuffer.data(si_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_0_yzzzzz = cbuffer.data(si_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_0_zzzzzz = cbuffer.data(si_geom_01_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxxxy, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxz, g_0_xxxyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyz, g_0_xxxyzzz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyz, g_0_xxyyzzz, g_0_xxyzz, g_0_xxyzzzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyz, g_0_xyyyzzz, g_0_xyyzz, g_0_xyyzzzz, g_0_xyzzz, g_0_xyzzzzz, g_0_xzzzz, g_0_y_0_xxxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_0_yyyyy, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyz, g_0_yyyyzzz, g_0_yyyzz, g_0_yyyzzzz, g_0_yyzzz, g_0_yyzzzzz, g_0_yzzzz, g_0_yzzzzzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xxxxxx[k] = g_0_xxxxxxy[k];

                g_0_y_0_xxxxxy[k] = -g_0_xxxxx[k] + g_0_xxxxxyy[k];

                g_0_y_0_xxxxxz[k] = g_0_xxxxxyz[k];

                g_0_y_0_xxxxyy[k] = -2.0 * g_0_xxxxy[k] + g_0_xxxxyyy[k];

                g_0_y_0_xxxxyz[k] = -g_0_xxxxz[k] + g_0_xxxxyyz[k];

                g_0_y_0_xxxxzz[k] = g_0_xxxxyzz[k];

                g_0_y_0_xxxyyy[k] = -3.0 * g_0_xxxyy[k] + g_0_xxxyyyy[k];

                g_0_y_0_xxxyyz[k] = -2.0 * g_0_xxxyz[k] + g_0_xxxyyyz[k];

                g_0_y_0_xxxyzz[k] = -g_0_xxxzz[k] + g_0_xxxyyzz[k];

                g_0_y_0_xxxzzz[k] = g_0_xxxyzzz[k];

                g_0_y_0_xxyyyy[k] = -4.0 * g_0_xxyyy[k] + g_0_xxyyyyy[k];

                g_0_y_0_xxyyyz[k] = -3.0 * g_0_xxyyz[k] + g_0_xxyyyyz[k];

                g_0_y_0_xxyyzz[k] = -2.0 * g_0_xxyzz[k] + g_0_xxyyyzz[k];

                g_0_y_0_xxyzzz[k] = -g_0_xxzzz[k] + g_0_xxyyzzz[k];

                g_0_y_0_xxzzzz[k] = g_0_xxyzzzz[k];

                g_0_y_0_xyyyyy[k] = -5.0 * g_0_xyyyy[k] + g_0_xyyyyyy[k];

                g_0_y_0_xyyyyz[k] = -4.0 * g_0_xyyyz[k] + g_0_xyyyyyz[k];

                g_0_y_0_xyyyzz[k] = -3.0 * g_0_xyyzz[k] + g_0_xyyyyzz[k];

                g_0_y_0_xyyzzz[k] = -2.0 * g_0_xyzzz[k] + g_0_xyyyzzz[k];

                g_0_y_0_xyzzzz[k] = -g_0_xzzzz[k] + g_0_xyyzzzz[k];

                g_0_y_0_xzzzzz[k] = g_0_xyzzzzz[k];

                g_0_y_0_yyyyyy[k] = -6.0 * g_0_yyyyy[k] + g_0_yyyyyyy[k];

                g_0_y_0_yyyyyz[k] = -5.0 * g_0_yyyyz[k] + g_0_yyyyyyz[k];

                g_0_y_0_yyyyzz[k] = -4.0 * g_0_yyyzz[k] + g_0_yyyyyzz[k];

                g_0_y_0_yyyzzz[k] = -3.0 * g_0_yyzzz[k] + g_0_yyyyzzz[k];

                g_0_y_0_yyzzzz[k] = -2.0 * g_0_yzzzz[k] + g_0_yyyzzzz[k];

                g_0_y_0_yzzzzz[k] = -g_0_zzzzz[k] + g_0_yyzzzzz[k];

                g_0_y_0_zzzzzz[k] = g_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_xxxxxx = cbuffer.data(si_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_z_0_xxxxxy = cbuffer.data(si_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_z_0_xxxxxz = cbuffer.data(si_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_z_0_xxxxyy = cbuffer.data(si_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_0_xxxxyz = cbuffer.data(si_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_0_xxxxzz = cbuffer.data(si_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_0_xxxyyy = cbuffer.data(si_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_0_xxxyyz = cbuffer.data(si_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_0_xxxyzz = cbuffer.data(si_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_0_xxxzzz = cbuffer.data(si_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_0_xxyyyy = cbuffer.data(si_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_0_xxyyyz = cbuffer.data(si_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_0_xxyyzz = cbuffer.data(si_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_0_xxyzzz = cbuffer.data(si_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_0_xxzzzz = cbuffer.data(si_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_0_xyyyyy = cbuffer.data(si_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_0_xyyyyz = cbuffer.data(si_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_0_xyyyzz = cbuffer.data(si_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_0_xyyzzz = cbuffer.data(si_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_0_xyzzzz = cbuffer.data(si_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_0_xzzzzz = cbuffer.data(si_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_0_yyyyyy = cbuffer.data(si_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_0_yyyyyz = cbuffer.data(si_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_0_yyyyzz = cbuffer.data(si_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_0_yyyzzz = cbuffer.data(si_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_0_yyzzzz = cbuffer.data(si_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_0_yzzzzz = cbuffer.data(si_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_0_zzzzzz = cbuffer.data(si_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxxxz, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxz, g_0_xxxxzzz, g_0_xxxyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyz, g_0_xxxyzzz, g_0_xxxzz, g_0_xxxzzzz, g_0_xxyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyz, g_0_xxyyzzz, g_0_xxyzz, g_0_xxyzzzz, g_0_xxzzz, g_0_xxzzzzz, g_0_xyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyz, g_0_xyyyzzz, g_0_xyyzz, g_0_xyyzzzz, g_0_xyzzz, g_0_xyzzzzz, g_0_xzzzz, g_0_xzzzzzz, g_0_yyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyz, g_0_yyyyzzz, g_0_yyyzz, g_0_yyyzzzz, g_0_yyzzz, g_0_yyzzzzz, g_0_yzzzz, g_0_yzzzzzz, g_0_z_0_xxxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_0_zzzzz, g_0_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xxxxxx[k] = g_0_xxxxxxz[k];

                g_0_z_0_xxxxxy[k] = g_0_xxxxxyz[k];

                g_0_z_0_xxxxxz[k] = -g_0_xxxxx[k] + g_0_xxxxxzz[k];

                g_0_z_0_xxxxyy[k] = g_0_xxxxyyz[k];

                g_0_z_0_xxxxyz[k] = -g_0_xxxxy[k] + g_0_xxxxyzz[k];

                g_0_z_0_xxxxzz[k] = -2.0 * g_0_xxxxz[k] + g_0_xxxxzzz[k];

                g_0_z_0_xxxyyy[k] = g_0_xxxyyyz[k];

                g_0_z_0_xxxyyz[k] = -g_0_xxxyy[k] + g_0_xxxyyzz[k];

                g_0_z_0_xxxyzz[k] = -2.0 * g_0_xxxyz[k] + g_0_xxxyzzz[k];

                g_0_z_0_xxxzzz[k] = -3.0 * g_0_xxxzz[k] + g_0_xxxzzzz[k];

                g_0_z_0_xxyyyy[k] = g_0_xxyyyyz[k];

                g_0_z_0_xxyyyz[k] = -g_0_xxyyy[k] + g_0_xxyyyzz[k];

                g_0_z_0_xxyyzz[k] = -2.0 * g_0_xxyyz[k] + g_0_xxyyzzz[k];

                g_0_z_0_xxyzzz[k] = -3.0 * g_0_xxyzz[k] + g_0_xxyzzzz[k];

                g_0_z_0_xxzzzz[k] = -4.0 * g_0_xxzzz[k] + g_0_xxzzzzz[k];

                g_0_z_0_xyyyyy[k] = g_0_xyyyyyz[k];

                g_0_z_0_xyyyyz[k] = -g_0_xyyyy[k] + g_0_xyyyyzz[k];

                g_0_z_0_xyyyzz[k] = -2.0 * g_0_xyyyz[k] + g_0_xyyyzzz[k];

                g_0_z_0_xyyzzz[k] = -3.0 * g_0_xyyzz[k] + g_0_xyyzzzz[k];

                g_0_z_0_xyzzzz[k] = -4.0 * g_0_xyzzz[k] + g_0_xyzzzzz[k];

                g_0_z_0_xzzzzz[k] = -5.0 * g_0_xzzzz[k] + g_0_xzzzzzz[k];

                g_0_z_0_yyyyyy[k] = g_0_yyyyyyz[k];

                g_0_z_0_yyyyyz[k] = -g_0_yyyyy[k] + g_0_yyyyyzz[k];

                g_0_z_0_yyyyzz[k] = -2.0 * g_0_yyyyz[k] + g_0_yyyyzzz[k];

                g_0_z_0_yyyzzz[k] = -3.0 * g_0_yyyzz[k] + g_0_yyyzzzz[k];

                g_0_z_0_yyzzzz[k] = -4.0 * g_0_yyzzz[k] + g_0_yyzzzzz[k];

                g_0_z_0_yzzzzz[k] = -5.0 * g_0_yzzzz[k] + g_0_yzzzzzz[k];

                g_0_z_0_zzzzzz[k] = -6.0 * g_0_zzzzz[k] + g_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

