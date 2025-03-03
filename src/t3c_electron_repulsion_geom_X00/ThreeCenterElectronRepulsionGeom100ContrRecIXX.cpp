#include "ThreeCenterElectronRepulsionGeom100ContrRecIXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_ixx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_ixx,
                                      const size_t idx_hxx,
                                      const size_t idx_kxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : HSS

            const auto h_off = idx_hxx + i * dcomps + j;

            auto g_xxxxx_0 = cbuffer.data(h_off + 0 * ccomps * dcomps);

            auto g_xxxxy_0 = cbuffer.data(h_off + 1 * ccomps * dcomps);

            auto g_xxxxz_0 = cbuffer.data(h_off + 2 * ccomps * dcomps);

            auto g_xxxyy_0 = cbuffer.data(h_off + 3 * ccomps * dcomps);

            auto g_xxxyz_0 = cbuffer.data(h_off + 4 * ccomps * dcomps);

            auto g_xxxzz_0 = cbuffer.data(h_off + 5 * ccomps * dcomps);

            auto g_xxyyy_0 = cbuffer.data(h_off + 6 * ccomps * dcomps);

            auto g_xxyyz_0 = cbuffer.data(h_off + 7 * ccomps * dcomps);

            auto g_xxyzz_0 = cbuffer.data(h_off + 8 * ccomps * dcomps);

            auto g_xxzzz_0 = cbuffer.data(h_off + 9 * ccomps * dcomps);

            auto g_xyyyy_0 = cbuffer.data(h_off + 10 * ccomps * dcomps);

            auto g_xyyyz_0 = cbuffer.data(h_off + 11 * ccomps * dcomps);

            auto g_xyyzz_0 = cbuffer.data(h_off + 12 * ccomps * dcomps);

            auto g_xyzzz_0 = cbuffer.data(h_off + 13 * ccomps * dcomps);

            auto g_xzzzz_0 = cbuffer.data(h_off + 14 * ccomps * dcomps);

            auto g_yyyyy_0 = cbuffer.data(h_off + 15 * ccomps * dcomps);

            auto g_yyyyz_0 = cbuffer.data(h_off + 16 * ccomps * dcomps);

            auto g_yyyzz_0 = cbuffer.data(h_off + 17 * ccomps * dcomps);

            auto g_yyzzz_0 = cbuffer.data(h_off + 18 * ccomps * dcomps);

            auto g_yzzzz_0 = cbuffer.data(h_off + 19 * ccomps * dcomps);

            auto g_zzzzz_0 = cbuffer.data(h_off + 20 * ccomps * dcomps);

            /// Set up components of auxilary buffer : KSS

            const auto k_off = idx_kxx + i * dcomps + j;

            auto g_xxxxxxx_0 = cbuffer.data(k_off + 0 * ccomps * dcomps);

            auto g_xxxxxxy_0 = cbuffer.data(k_off + 1 * ccomps * dcomps);

            auto g_xxxxxxz_0 = cbuffer.data(k_off + 2 * ccomps * dcomps);

            auto g_xxxxxyy_0 = cbuffer.data(k_off + 3 * ccomps * dcomps);

            auto g_xxxxxyz_0 = cbuffer.data(k_off + 4 * ccomps * dcomps);

            auto g_xxxxxzz_0 = cbuffer.data(k_off + 5 * ccomps * dcomps);

            auto g_xxxxyyy_0 = cbuffer.data(k_off + 6 * ccomps * dcomps);

            auto g_xxxxyyz_0 = cbuffer.data(k_off + 7 * ccomps * dcomps);

            auto g_xxxxyzz_0 = cbuffer.data(k_off + 8 * ccomps * dcomps);

            auto g_xxxxzzz_0 = cbuffer.data(k_off + 9 * ccomps * dcomps);

            auto g_xxxyyyy_0 = cbuffer.data(k_off + 10 * ccomps * dcomps);

            auto g_xxxyyyz_0 = cbuffer.data(k_off + 11 * ccomps * dcomps);

            auto g_xxxyyzz_0 = cbuffer.data(k_off + 12 * ccomps * dcomps);

            auto g_xxxyzzz_0 = cbuffer.data(k_off + 13 * ccomps * dcomps);

            auto g_xxxzzzz_0 = cbuffer.data(k_off + 14 * ccomps * dcomps);

            auto g_xxyyyyy_0 = cbuffer.data(k_off + 15 * ccomps * dcomps);

            auto g_xxyyyyz_0 = cbuffer.data(k_off + 16 * ccomps * dcomps);

            auto g_xxyyyzz_0 = cbuffer.data(k_off + 17 * ccomps * dcomps);

            auto g_xxyyzzz_0 = cbuffer.data(k_off + 18 * ccomps * dcomps);

            auto g_xxyzzzz_0 = cbuffer.data(k_off + 19 * ccomps * dcomps);

            auto g_xxzzzzz_0 = cbuffer.data(k_off + 20 * ccomps * dcomps);

            auto g_xyyyyyy_0 = cbuffer.data(k_off + 21 * ccomps * dcomps);

            auto g_xyyyyyz_0 = cbuffer.data(k_off + 22 * ccomps * dcomps);

            auto g_xyyyyzz_0 = cbuffer.data(k_off + 23 * ccomps * dcomps);

            auto g_xyyyzzz_0 = cbuffer.data(k_off + 24 * ccomps * dcomps);

            auto g_xyyzzzz_0 = cbuffer.data(k_off + 25 * ccomps * dcomps);

            auto g_xyzzzzz_0 = cbuffer.data(k_off + 26 * ccomps * dcomps);

            auto g_xzzzzzz_0 = cbuffer.data(k_off + 27 * ccomps * dcomps);

            auto g_yyyyyyy_0 = cbuffer.data(k_off + 28 * ccomps * dcomps);

            auto g_yyyyyyz_0 = cbuffer.data(k_off + 29 * ccomps * dcomps);

            auto g_yyyyyzz_0 = cbuffer.data(k_off + 30 * ccomps * dcomps);

            auto g_yyyyzzz_0 = cbuffer.data(k_off + 31 * ccomps * dcomps);

            auto g_yyyzzzz_0 = cbuffer.data(k_off + 32 * ccomps * dcomps);

            auto g_yyzzzzz_0 = cbuffer.data(k_off + 33 * ccomps * dcomps);

            auto g_yzzzzzz_0 = cbuffer.data(k_off + 34 * ccomps * dcomps);

            auto g_zzzzzzz_0 = cbuffer.data(k_off + 35 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ixx

            const auto i_geom_100_off = idx_geom_100_ixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxxxxx_0 = cbuffer.data(i_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxxy_0 = cbuffer.data(i_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxxz_0 = cbuffer.data(i_geom_100_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxxyy_0 = cbuffer.data(i_geom_100_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxxyz_0 = cbuffer.data(i_geom_100_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxxzz_0 = cbuffer.data(i_geom_100_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxxyyy_0 = cbuffer.data(i_geom_100_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxxyyz_0 = cbuffer.data(i_geom_100_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxxyzz_0 = cbuffer.data(i_geom_100_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxxzzz_0 = cbuffer.data(i_geom_100_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xxyyyy_0 = cbuffer.data(i_geom_100_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xxyyyz_0 = cbuffer.data(i_geom_100_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xxyyzz_0 = cbuffer.data(i_geom_100_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xxyzzz_0 = cbuffer.data(i_geom_100_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xxzzzz_0 = cbuffer.data(i_geom_100_off + 14 * ccomps * dcomps);

            auto g_x_0_0_xyyyyy_0 = cbuffer.data(i_geom_100_off + 15 * ccomps * dcomps);

            auto g_x_0_0_xyyyyz_0 = cbuffer.data(i_geom_100_off + 16 * ccomps * dcomps);

            auto g_x_0_0_xyyyzz_0 = cbuffer.data(i_geom_100_off + 17 * ccomps * dcomps);

            auto g_x_0_0_xyyzzz_0 = cbuffer.data(i_geom_100_off + 18 * ccomps * dcomps);

            auto g_x_0_0_xyzzzz_0 = cbuffer.data(i_geom_100_off + 19 * ccomps * dcomps);

            auto g_x_0_0_xzzzzz_0 = cbuffer.data(i_geom_100_off + 20 * ccomps * dcomps);

            auto g_x_0_0_yyyyyy_0 = cbuffer.data(i_geom_100_off + 21 * ccomps * dcomps);

            auto g_x_0_0_yyyyyz_0 = cbuffer.data(i_geom_100_off + 22 * ccomps * dcomps);

            auto g_x_0_0_yyyyzz_0 = cbuffer.data(i_geom_100_off + 23 * ccomps * dcomps);

            auto g_x_0_0_yyyzzz_0 = cbuffer.data(i_geom_100_off + 24 * ccomps * dcomps);

            auto g_x_0_0_yyzzzz_0 = cbuffer.data(i_geom_100_off + 25 * ccomps * dcomps);

            auto g_x_0_0_yzzzzz_0 = cbuffer.data(i_geom_100_off + 26 * ccomps * dcomps);

            auto g_x_0_0_zzzzzz_0 = cbuffer.data(i_geom_100_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx_0, g_x_0_0_xxxxxy_0, g_x_0_0_xxxxxz_0, g_x_0_0_xxxxyy_0, g_x_0_0_xxxxyz_0, g_x_0_0_xxxxzz_0, g_x_0_0_xxxyyy_0, g_x_0_0_xxxyyz_0, g_x_0_0_xxxyzz_0, g_x_0_0_xxxzzz_0, g_x_0_0_xxyyyy_0, g_x_0_0_xxyyyz_0, g_x_0_0_xxyyzz_0, g_x_0_0_xxyzzz_0, g_x_0_0_xxzzzz_0, g_x_0_0_xyyyyy_0, g_x_0_0_xyyyyz_0, g_x_0_0_xyyyzz_0, g_x_0_0_xyyzzz_0, g_x_0_0_xyzzzz_0, g_x_0_0_xzzzzz_0, g_x_0_0_yyyyyy_0, g_x_0_0_yyyyyz_0, g_x_0_0_yyyyzz_0, g_x_0_0_yyyzzz_0, g_x_0_0_yyzzzz_0, g_x_0_0_yzzzzz_0, g_x_0_0_zzzzzz_0, g_xxxxx_0, g_xxxxxxx_0, g_xxxxxxy_0, g_xxxxxxz_0, g_xxxxxyy_0, g_xxxxxyz_0, g_xxxxxzz_0, g_xxxxy_0, g_xxxxyyy_0, g_xxxxyyz_0, g_xxxxyzz_0, g_xxxxz_0, g_xxxxzzz_0, g_xxxyy_0, g_xxxyyyy_0, g_xxxyyyz_0, g_xxxyyzz_0, g_xxxyz_0, g_xxxyzzz_0, g_xxxzz_0, g_xxxzzzz_0, g_xxyyy_0, g_xxyyyyy_0, g_xxyyyyz_0, g_xxyyyzz_0, g_xxyyz_0, g_xxyyzzz_0, g_xxyzz_0, g_xxyzzzz_0, g_xxzzz_0, g_xxzzzzz_0, g_xyyyy_0, g_xyyyyyy_0, g_xyyyyyz_0, g_xyyyyzz_0, g_xyyyz_0, g_xyyyzzz_0, g_xyyzz_0, g_xyyzzzz_0, g_xyzzz_0, g_xyzzzzz_0, g_xzzzz_0, g_xzzzzzz_0, g_yyyyy_0, g_yyyyz_0, g_yyyzz_0, g_yyzzz_0, g_yzzzz_0, g_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxxxx_0[k] = -6.0 * g_xxxxx_0[k] + g_xxxxxxx_0[k];

                g_x_0_0_xxxxxy_0[k] = -5.0 * g_xxxxy_0[k] + g_xxxxxxy_0[k];

                g_x_0_0_xxxxxz_0[k] = -5.0 * g_xxxxz_0[k] + g_xxxxxxz_0[k];

                g_x_0_0_xxxxyy_0[k] = -4.0 * g_xxxyy_0[k] + g_xxxxxyy_0[k];

                g_x_0_0_xxxxyz_0[k] = -4.0 * g_xxxyz_0[k] + g_xxxxxyz_0[k];

                g_x_0_0_xxxxzz_0[k] = -4.0 * g_xxxzz_0[k] + g_xxxxxzz_0[k];

                g_x_0_0_xxxyyy_0[k] = -3.0 * g_xxyyy_0[k] + g_xxxxyyy_0[k];

                g_x_0_0_xxxyyz_0[k] = -3.0 * g_xxyyz_0[k] + g_xxxxyyz_0[k];

                g_x_0_0_xxxyzz_0[k] = -3.0 * g_xxyzz_0[k] + g_xxxxyzz_0[k];

                g_x_0_0_xxxzzz_0[k] = -3.0 * g_xxzzz_0[k] + g_xxxxzzz_0[k];

                g_x_0_0_xxyyyy_0[k] = -2.0 * g_xyyyy_0[k] + g_xxxyyyy_0[k];

                g_x_0_0_xxyyyz_0[k] = -2.0 * g_xyyyz_0[k] + g_xxxyyyz_0[k];

                g_x_0_0_xxyyzz_0[k] = -2.0 * g_xyyzz_0[k] + g_xxxyyzz_0[k];

                g_x_0_0_xxyzzz_0[k] = -2.0 * g_xyzzz_0[k] + g_xxxyzzz_0[k];

                g_x_0_0_xxzzzz_0[k] = -2.0 * g_xzzzz_0[k] + g_xxxzzzz_0[k];

                g_x_0_0_xyyyyy_0[k] = -g_yyyyy_0[k] + g_xxyyyyy_0[k];

                g_x_0_0_xyyyyz_0[k] = -g_yyyyz_0[k] + g_xxyyyyz_0[k];

                g_x_0_0_xyyyzz_0[k] = -g_yyyzz_0[k] + g_xxyyyzz_0[k];

                g_x_0_0_xyyzzz_0[k] = -g_yyzzz_0[k] + g_xxyyzzz_0[k];

                g_x_0_0_xyzzzz_0[k] = -g_yzzzz_0[k] + g_xxyzzzz_0[k];

                g_x_0_0_xzzzzz_0[k] = -g_zzzzz_0[k] + g_xxzzzzz_0[k];

                g_x_0_0_yyyyyy_0[k] = g_xyyyyyy_0[k];

                g_x_0_0_yyyyyz_0[k] = g_xyyyyyz_0[k];

                g_x_0_0_yyyyzz_0[k] = g_xyyyyzz_0[k];

                g_x_0_0_yyyzzz_0[k] = g_xyyyzzz_0[k];

                g_x_0_0_yyzzzz_0[k] = g_xyyzzzz_0[k];

                g_x_0_0_yzzzzz_0[k] = g_xyzzzzz_0[k];

                g_x_0_0_zzzzzz_0[k] = g_xzzzzzz_0[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxxxxx_0 = cbuffer.data(i_geom_100_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxxxxy_0 = cbuffer.data(i_geom_100_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxxxxz_0 = cbuffer.data(i_geom_100_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xxxxyy_0 = cbuffer.data(i_geom_100_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xxxxyz_0 = cbuffer.data(i_geom_100_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xxxxzz_0 = cbuffer.data(i_geom_100_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xxxyyy_0 = cbuffer.data(i_geom_100_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xxxyyz_0 = cbuffer.data(i_geom_100_off + 35 * ccomps * dcomps);

            auto g_y_0_0_xxxyzz_0 = cbuffer.data(i_geom_100_off + 36 * ccomps * dcomps);

            auto g_y_0_0_xxxzzz_0 = cbuffer.data(i_geom_100_off + 37 * ccomps * dcomps);

            auto g_y_0_0_xxyyyy_0 = cbuffer.data(i_geom_100_off + 38 * ccomps * dcomps);

            auto g_y_0_0_xxyyyz_0 = cbuffer.data(i_geom_100_off + 39 * ccomps * dcomps);

            auto g_y_0_0_xxyyzz_0 = cbuffer.data(i_geom_100_off + 40 * ccomps * dcomps);

            auto g_y_0_0_xxyzzz_0 = cbuffer.data(i_geom_100_off + 41 * ccomps * dcomps);

            auto g_y_0_0_xxzzzz_0 = cbuffer.data(i_geom_100_off + 42 * ccomps * dcomps);

            auto g_y_0_0_xyyyyy_0 = cbuffer.data(i_geom_100_off + 43 * ccomps * dcomps);

            auto g_y_0_0_xyyyyz_0 = cbuffer.data(i_geom_100_off + 44 * ccomps * dcomps);

            auto g_y_0_0_xyyyzz_0 = cbuffer.data(i_geom_100_off + 45 * ccomps * dcomps);

            auto g_y_0_0_xyyzzz_0 = cbuffer.data(i_geom_100_off + 46 * ccomps * dcomps);

            auto g_y_0_0_xyzzzz_0 = cbuffer.data(i_geom_100_off + 47 * ccomps * dcomps);

            auto g_y_0_0_xzzzzz_0 = cbuffer.data(i_geom_100_off + 48 * ccomps * dcomps);

            auto g_y_0_0_yyyyyy_0 = cbuffer.data(i_geom_100_off + 49 * ccomps * dcomps);

            auto g_y_0_0_yyyyyz_0 = cbuffer.data(i_geom_100_off + 50 * ccomps * dcomps);

            auto g_y_0_0_yyyyzz_0 = cbuffer.data(i_geom_100_off + 51 * ccomps * dcomps);

            auto g_y_0_0_yyyzzz_0 = cbuffer.data(i_geom_100_off + 52 * ccomps * dcomps);

            auto g_y_0_0_yyzzzz_0 = cbuffer.data(i_geom_100_off + 53 * ccomps * dcomps);

            auto g_y_0_0_yzzzzz_0 = cbuffer.data(i_geom_100_off + 54 * ccomps * dcomps);

            auto g_y_0_0_zzzzzz_0 = cbuffer.data(i_geom_100_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxxx_0, g_xxxxxxy_0, g_xxxxxyy_0, g_xxxxxyz_0, g_xxxxy_0, g_xxxxyyy_0, g_xxxxyyz_0, g_xxxxyzz_0, g_xxxxz_0, g_xxxyy_0, g_xxxyyyy_0, g_xxxyyyz_0, g_xxxyyzz_0, g_xxxyz_0, g_xxxyzzz_0, g_xxxzz_0, g_xxyyy_0, g_xxyyyyy_0, g_xxyyyyz_0, g_xxyyyzz_0, g_xxyyz_0, g_xxyyzzz_0, g_xxyzz_0, g_xxyzzzz_0, g_xxzzz_0, g_xyyyy_0, g_xyyyyyy_0, g_xyyyyyz_0, g_xyyyyzz_0, g_xyyyz_0, g_xyyyzzz_0, g_xyyzz_0, g_xyyzzzz_0, g_xyzzz_0, g_xyzzzzz_0, g_xzzzz_0, g_y_0_0_xxxxxx_0, g_y_0_0_xxxxxy_0, g_y_0_0_xxxxxz_0, g_y_0_0_xxxxyy_0, g_y_0_0_xxxxyz_0, g_y_0_0_xxxxzz_0, g_y_0_0_xxxyyy_0, g_y_0_0_xxxyyz_0, g_y_0_0_xxxyzz_0, g_y_0_0_xxxzzz_0, g_y_0_0_xxyyyy_0, g_y_0_0_xxyyyz_0, g_y_0_0_xxyyzz_0, g_y_0_0_xxyzzz_0, g_y_0_0_xxzzzz_0, g_y_0_0_xyyyyy_0, g_y_0_0_xyyyyz_0, g_y_0_0_xyyyzz_0, g_y_0_0_xyyzzz_0, g_y_0_0_xyzzzz_0, g_y_0_0_xzzzzz_0, g_y_0_0_yyyyyy_0, g_y_0_0_yyyyyz_0, g_y_0_0_yyyyzz_0, g_y_0_0_yyyzzz_0, g_y_0_0_yyzzzz_0, g_y_0_0_yzzzzz_0, g_y_0_0_zzzzzz_0, g_yyyyy_0, g_yyyyyyy_0, g_yyyyyyz_0, g_yyyyyzz_0, g_yyyyz_0, g_yyyyzzz_0, g_yyyzz_0, g_yyyzzzz_0, g_yyzzz_0, g_yyzzzzz_0, g_yzzzz_0, g_yzzzzzz_0, g_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxxxx_0[k] = g_xxxxxxy_0[k];

                g_y_0_0_xxxxxy_0[k] = -g_xxxxx_0[k] + g_xxxxxyy_0[k];

                g_y_0_0_xxxxxz_0[k] = g_xxxxxyz_0[k];

                g_y_0_0_xxxxyy_0[k] = -2.0 * g_xxxxy_0[k] + g_xxxxyyy_0[k];

                g_y_0_0_xxxxyz_0[k] = -g_xxxxz_0[k] + g_xxxxyyz_0[k];

                g_y_0_0_xxxxzz_0[k] = g_xxxxyzz_0[k];

                g_y_0_0_xxxyyy_0[k] = -3.0 * g_xxxyy_0[k] + g_xxxyyyy_0[k];

                g_y_0_0_xxxyyz_0[k] = -2.0 * g_xxxyz_0[k] + g_xxxyyyz_0[k];

                g_y_0_0_xxxyzz_0[k] = -g_xxxzz_0[k] + g_xxxyyzz_0[k];

                g_y_0_0_xxxzzz_0[k] = g_xxxyzzz_0[k];

                g_y_0_0_xxyyyy_0[k] = -4.0 * g_xxyyy_0[k] + g_xxyyyyy_0[k];

                g_y_0_0_xxyyyz_0[k] = -3.0 * g_xxyyz_0[k] + g_xxyyyyz_0[k];

                g_y_0_0_xxyyzz_0[k] = -2.0 * g_xxyzz_0[k] + g_xxyyyzz_0[k];

                g_y_0_0_xxyzzz_0[k] = -g_xxzzz_0[k] + g_xxyyzzz_0[k];

                g_y_0_0_xxzzzz_0[k] = g_xxyzzzz_0[k];

                g_y_0_0_xyyyyy_0[k] = -5.0 * g_xyyyy_0[k] + g_xyyyyyy_0[k];

                g_y_0_0_xyyyyz_0[k] = -4.0 * g_xyyyz_0[k] + g_xyyyyyz_0[k];

                g_y_0_0_xyyyzz_0[k] = -3.0 * g_xyyzz_0[k] + g_xyyyyzz_0[k];

                g_y_0_0_xyyzzz_0[k] = -2.0 * g_xyzzz_0[k] + g_xyyyzzz_0[k];

                g_y_0_0_xyzzzz_0[k] = -g_xzzzz_0[k] + g_xyyzzzz_0[k];

                g_y_0_0_xzzzzz_0[k] = g_xyzzzzz_0[k];

                g_y_0_0_yyyyyy_0[k] = -6.0 * g_yyyyy_0[k] + g_yyyyyyy_0[k];

                g_y_0_0_yyyyyz_0[k] = -5.0 * g_yyyyz_0[k] + g_yyyyyyz_0[k];

                g_y_0_0_yyyyzz_0[k] = -4.0 * g_yyyzz_0[k] + g_yyyyyzz_0[k];

                g_y_0_0_yyyzzz_0[k] = -3.0 * g_yyzzz_0[k] + g_yyyyzzz_0[k];

                g_y_0_0_yyzzzz_0[k] = -2.0 * g_yzzzz_0[k] + g_yyyzzzz_0[k];

                g_y_0_0_yzzzzz_0[k] = -g_zzzzz_0[k] + g_yyzzzzz_0[k];

                g_y_0_0_zzzzzz_0[k] = g_yzzzzzz_0[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxxxxx_0 = cbuffer.data(i_geom_100_off + 56 * ccomps * dcomps);

            auto g_z_0_0_xxxxxy_0 = cbuffer.data(i_geom_100_off + 57 * ccomps * dcomps);

            auto g_z_0_0_xxxxxz_0 = cbuffer.data(i_geom_100_off + 58 * ccomps * dcomps);

            auto g_z_0_0_xxxxyy_0 = cbuffer.data(i_geom_100_off + 59 * ccomps * dcomps);

            auto g_z_0_0_xxxxyz_0 = cbuffer.data(i_geom_100_off + 60 * ccomps * dcomps);

            auto g_z_0_0_xxxxzz_0 = cbuffer.data(i_geom_100_off + 61 * ccomps * dcomps);

            auto g_z_0_0_xxxyyy_0 = cbuffer.data(i_geom_100_off + 62 * ccomps * dcomps);

            auto g_z_0_0_xxxyyz_0 = cbuffer.data(i_geom_100_off + 63 * ccomps * dcomps);

            auto g_z_0_0_xxxyzz_0 = cbuffer.data(i_geom_100_off + 64 * ccomps * dcomps);

            auto g_z_0_0_xxxzzz_0 = cbuffer.data(i_geom_100_off + 65 * ccomps * dcomps);

            auto g_z_0_0_xxyyyy_0 = cbuffer.data(i_geom_100_off + 66 * ccomps * dcomps);

            auto g_z_0_0_xxyyyz_0 = cbuffer.data(i_geom_100_off + 67 * ccomps * dcomps);

            auto g_z_0_0_xxyyzz_0 = cbuffer.data(i_geom_100_off + 68 * ccomps * dcomps);

            auto g_z_0_0_xxyzzz_0 = cbuffer.data(i_geom_100_off + 69 * ccomps * dcomps);

            auto g_z_0_0_xxzzzz_0 = cbuffer.data(i_geom_100_off + 70 * ccomps * dcomps);

            auto g_z_0_0_xyyyyy_0 = cbuffer.data(i_geom_100_off + 71 * ccomps * dcomps);

            auto g_z_0_0_xyyyyz_0 = cbuffer.data(i_geom_100_off + 72 * ccomps * dcomps);

            auto g_z_0_0_xyyyzz_0 = cbuffer.data(i_geom_100_off + 73 * ccomps * dcomps);

            auto g_z_0_0_xyyzzz_0 = cbuffer.data(i_geom_100_off + 74 * ccomps * dcomps);

            auto g_z_0_0_xyzzzz_0 = cbuffer.data(i_geom_100_off + 75 * ccomps * dcomps);

            auto g_z_0_0_xzzzzz_0 = cbuffer.data(i_geom_100_off + 76 * ccomps * dcomps);

            auto g_z_0_0_yyyyyy_0 = cbuffer.data(i_geom_100_off + 77 * ccomps * dcomps);

            auto g_z_0_0_yyyyyz_0 = cbuffer.data(i_geom_100_off + 78 * ccomps * dcomps);

            auto g_z_0_0_yyyyzz_0 = cbuffer.data(i_geom_100_off + 79 * ccomps * dcomps);

            auto g_z_0_0_yyyzzz_0 = cbuffer.data(i_geom_100_off + 80 * ccomps * dcomps);

            auto g_z_0_0_yyzzzz_0 = cbuffer.data(i_geom_100_off + 81 * ccomps * dcomps);

            auto g_z_0_0_yzzzzz_0 = cbuffer.data(i_geom_100_off + 82 * ccomps * dcomps);

            auto g_z_0_0_zzzzzz_0 = cbuffer.data(i_geom_100_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxxx_0, g_xxxxxxz_0, g_xxxxxyz_0, g_xxxxxzz_0, g_xxxxy_0, g_xxxxyyz_0, g_xxxxyzz_0, g_xxxxz_0, g_xxxxzzz_0, g_xxxyy_0, g_xxxyyyz_0, g_xxxyyzz_0, g_xxxyz_0, g_xxxyzzz_0, g_xxxzz_0, g_xxxzzzz_0, g_xxyyy_0, g_xxyyyyz_0, g_xxyyyzz_0, g_xxyyz_0, g_xxyyzzz_0, g_xxyzz_0, g_xxyzzzz_0, g_xxzzz_0, g_xxzzzzz_0, g_xyyyy_0, g_xyyyyyz_0, g_xyyyyzz_0, g_xyyyz_0, g_xyyyzzz_0, g_xyyzz_0, g_xyyzzzz_0, g_xyzzz_0, g_xyzzzzz_0, g_xzzzz_0, g_xzzzzzz_0, g_yyyyy_0, g_yyyyyyz_0, g_yyyyyzz_0, g_yyyyz_0, g_yyyyzzz_0, g_yyyzz_0, g_yyyzzzz_0, g_yyzzz_0, g_yyzzzzz_0, g_yzzzz_0, g_yzzzzzz_0, g_z_0_0_xxxxxx_0, g_z_0_0_xxxxxy_0, g_z_0_0_xxxxxz_0, g_z_0_0_xxxxyy_0, g_z_0_0_xxxxyz_0, g_z_0_0_xxxxzz_0, g_z_0_0_xxxyyy_0, g_z_0_0_xxxyyz_0, g_z_0_0_xxxyzz_0, g_z_0_0_xxxzzz_0, g_z_0_0_xxyyyy_0, g_z_0_0_xxyyyz_0, g_z_0_0_xxyyzz_0, g_z_0_0_xxyzzz_0, g_z_0_0_xxzzzz_0, g_z_0_0_xyyyyy_0, g_z_0_0_xyyyyz_0, g_z_0_0_xyyyzz_0, g_z_0_0_xyyzzz_0, g_z_0_0_xyzzzz_0, g_z_0_0_xzzzzz_0, g_z_0_0_yyyyyy_0, g_z_0_0_yyyyyz_0, g_z_0_0_yyyyzz_0, g_z_0_0_yyyzzz_0, g_z_0_0_yyzzzz_0, g_z_0_0_yzzzzz_0, g_z_0_0_zzzzzz_0, g_zzzzz_0, g_zzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxxxx_0[k] = g_xxxxxxz_0[k];

                g_z_0_0_xxxxxy_0[k] = g_xxxxxyz_0[k];

                g_z_0_0_xxxxxz_0[k] = -g_xxxxx_0[k] + g_xxxxxzz_0[k];

                g_z_0_0_xxxxyy_0[k] = g_xxxxyyz_0[k];

                g_z_0_0_xxxxyz_0[k] = -g_xxxxy_0[k] + g_xxxxyzz_0[k];

                g_z_0_0_xxxxzz_0[k] = -2.0 * g_xxxxz_0[k] + g_xxxxzzz_0[k];

                g_z_0_0_xxxyyy_0[k] = g_xxxyyyz_0[k];

                g_z_0_0_xxxyyz_0[k] = -g_xxxyy_0[k] + g_xxxyyzz_0[k];

                g_z_0_0_xxxyzz_0[k] = -2.0 * g_xxxyz_0[k] + g_xxxyzzz_0[k];

                g_z_0_0_xxxzzz_0[k] = -3.0 * g_xxxzz_0[k] + g_xxxzzzz_0[k];

                g_z_0_0_xxyyyy_0[k] = g_xxyyyyz_0[k];

                g_z_0_0_xxyyyz_0[k] = -g_xxyyy_0[k] + g_xxyyyzz_0[k];

                g_z_0_0_xxyyzz_0[k] = -2.0 * g_xxyyz_0[k] + g_xxyyzzz_0[k];

                g_z_0_0_xxyzzz_0[k] = -3.0 * g_xxyzz_0[k] + g_xxyzzzz_0[k];

                g_z_0_0_xxzzzz_0[k] = -4.0 * g_xxzzz_0[k] + g_xxzzzzz_0[k];

                g_z_0_0_xyyyyy_0[k] = g_xyyyyyz_0[k];

                g_z_0_0_xyyyyz_0[k] = -g_xyyyy_0[k] + g_xyyyyzz_0[k];

                g_z_0_0_xyyyzz_0[k] = -2.0 * g_xyyyz_0[k] + g_xyyyzzz_0[k];

                g_z_0_0_xyyzzz_0[k] = -3.0 * g_xyyzz_0[k] + g_xyyzzzz_0[k];

                g_z_0_0_xyzzzz_0[k] = -4.0 * g_xyzzz_0[k] + g_xyzzzzz_0[k];

                g_z_0_0_xzzzzz_0[k] = -5.0 * g_xzzzz_0[k] + g_xzzzzzz_0[k];

                g_z_0_0_yyyyyy_0[k] = g_yyyyyyz_0[k];

                g_z_0_0_yyyyyz_0[k] = -g_yyyyy_0[k] + g_yyyyyzz_0[k];

                g_z_0_0_yyyyzz_0[k] = -2.0 * g_yyyyz_0[k] + g_yyyyzzz_0[k];

                g_z_0_0_yyyzzz_0[k] = -3.0 * g_yyyzz_0[k] + g_yyyzzzz_0[k];

                g_z_0_0_yyzzzz_0[k] = -4.0 * g_yyzzz_0[k] + g_yyzzzzz_0[k];

                g_z_0_0_yzzzzz_0[k] = -5.0 * g_yzzzz_0[k] + g_yzzzzzz_0[k];

                g_z_0_0_zzzzzz_0[k] = -6.0 * g_zzzzz_0[k] + g_zzzzzzz_0[k];
            }
        }
    }
}

} // t3ceri namespace

