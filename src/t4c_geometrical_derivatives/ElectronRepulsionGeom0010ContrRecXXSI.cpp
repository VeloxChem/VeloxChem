#include "ElectronRepulsionGeom0010ContrRecXXSI.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsi,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsi,
                                            const size_t idx_xxsk,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSI

            const auto si_off = idx_xxsi + (i * bcomps + j) * 28;

            auto g_0_xxxxxx = pbuffer.data(si_off + 0);

            auto g_0_xxxxxy = pbuffer.data(si_off + 1);

            auto g_0_xxxxxz = pbuffer.data(si_off + 2);

            auto g_0_xxxxyy = pbuffer.data(si_off + 3);

            auto g_0_xxxxyz = pbuffer.data(si_off + 4);

            auto g_0_xxxxzz = pbuffer.data(si_off + 5);

            auto g_0_xxxyyy = pbuffer.data(si_off + 6);

            auto g_0_xxxyyz = pbuffer.data(si_off + 7);

            auto g_0_xxxyzz = pbuffer.data(si_off + 8);

            auto g_0_xxxzzz = pbuffer.data(si_off + 9);

            auto g_0_xxyyyy = pbuffer.data(si_off + 10);

            auto g_0_xxyyyz = pbuffer.data(si_off + 11);

            auto g_0_xxyyzz = pbuffer.data(si_off + 12);

            auto g_0_xxyzzz = pbuffer.data(si_off + 13);

            auto g_0_xxzzzz = pbuffer.data(si_off + 14);

            auto g_0_xyyyyy = pbuffer.data(si_off + 15);

            auto g_0_xyyyyz = pbuffer.data(si_off + 16);

            auto g_0_xyyyzz = pbuffer.data(si_off + 17);

            auto g_0_xyyzzz = pbuffer.data(si_off + 18);

            auto g_0_xyzzzz = pbuffer.data(si_off + 19);

            auto g_0_xzzzzz = pbuffer.data(si_off + 20);

            auto g_0_yyyyyy = pbuffer.data(si_off + 21);

            auto g_0_yyyyyz = pbuffer.data(si_off + 22);

            auto g_0_yyyyzz = pbuffer.data(si_off + 23);

            auto g_0_yyyzzz = pbuffer.data(si_off + 24);

            auto g_0_yyzzzz = pbuffer.data(si_off + 25);

            auto g_0_yzzzzz = pbuffer.data(si_off + 26);

            auto g_0_zzzzzz = pbuffer.data(si_off + 27);

            /// Set up components of auxilary buffer : SSSK

            const auto sk_off = idx_xxsk + (i * bcomps + j) * 36;

            auto g_0_xxxxxxx = pbuffer.data(sk_off + 0);

            auto g_0_xxxxxxy = pbuffer.data(sk_off + 1);

            auto g_0_xxxxxxz = pbuffer.data(sk_off + 2);

            auto g_0_xxxxxyy = pbuffer.data(sk_off + 3);

            auto g_0_xxxxxyz = pbuffer.data(sk_off + 4);

            auto g_0_xxxxxzz = pbuffer.data(sk_off + 5);

            auto g_0_xxxxyyy = pbuffer.data(sk_off + 6);

            auto g_0_xxxxyyz = pbuffer.data(sk_off + 7);

            auto g_0_xxxxyzz = pbuffer.data(sk_off + 8);

            auto g_0_xxxxzzz = pbuffer.data(sk_off + 9);

            auto g_0_xxxyyyy = pbuffer.data(sk_off + 10);

            auto g_0_xxxyyyz = pbuffer.data(sk_off + 11);

            auto g_0_xxxyyzz = pbuffer.data(sk_off + 12);

            auto g_0_xxxyzzz = pbuffer.data(sk_off + 13);

            auto g_0_xxxzzzz = pbuffer.data(sk_off + 14);

            auto g_0_xxyyyyy = pbuffer.data(sk_off + 15);

            auto g_0_xxyyyyz = pbuffer.data(sk_off + 16);

            auto g_0_xxyyyzz = pbuffer.data(sk_off + 17);

            auto g_0_xxyyzzz = pbuffer.data(sk_off + 18);

            auto g_0_xxyzzzz = pbuffer.data(sk_off + 19);

            auto g_0_xxzzzzz = pbuffer.data(sk_off + 20);

            auto g_0_xyyyyyy = pbuffer.data(sk_off + 21);

            auto g_0_xyyyyyz = pbuffer.data(sk_off + 22);

            auto g_0_xyyyyzz = pbuffer.data(sk_off + 23);

            auto g_0_xyyyzzz = pbuffer.data(sk_off + 24);

            auto g_0_xyyzzzz = pbuffer.data(sk_off + 25);

            auto g_0_xyzzzzz = pbuffer.data(sk_off + 26);

            auto g_0_xzzzzzz = pbuffer.data(sk_off + 27);

            auto g_0_yyyyyyy = pbuffer.data(sk_off + 28);

            auto g_0_yyyyyyz = pbuffer.data(sk_off + 29);

            auto g_0_yyyyyzz = pbuffer.data(sk_off + 30);

            auto g_0_yyyyzzz = pbuffer.data(sk_off + 31);

            auto g_0_yyyzzzz = pbuffer.data(sk_off + 32);

            auto g_0_yyzzzzz = pbuffer.data(sk_off + 33);

            auto g_0_yzzzzzz = pbuffer.data(sk_off + 34);

            auto g_0_zzzzzzz = pbuffer.data(sk_off + 35);

            /// set up bra offset for contr_buffer_xxsi

            const auto si_geom_10_off = idx_geom_10_xxsi + (i * bcomps + j) * 28;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_0_xxxxxx, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxy, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxxzz, g_0_xxxxyy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxxzzz, g_0_xxxyyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxxzzzz, g_0_xxyyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xxzzzzz, g_0_xyyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_xzzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxxxx[k] = -g_0_xxxxxx[k] * cd_x[k] + g_0_xxxxxxx[k];

                g_x_0_0_xxxxxy[k] = -g_0_xxxxxy[k] * cd_x[k] + g_0_xxxxxxy[k];

                g_x_0_0_xxxxxz[k] = -g_0_xxxxxz[k] * cd_x[k] + g_0_xxxxxxz[k];

                g_x_0_0_xxxxyy[k] = -g_0_xxxxyy[k] * cd_x[k] + g_0_xxxxxyy[k];

                g_x_0_0_xxxxyz[k] = -g_0_xxxxyz[k] * cd_x[k] + g_0_xxxxxyz[k];

                g_x_0_0_xxxxzz[k] = -g_0_xxxxzz[k] * cd_x[k] + g_0_xxxxxzz[k];

                g_x_0_0_xxxyyy[k] = -g_0_xxxyyy[k] * cd_x[k] + g_0_xxxxyyy[k];

                g_x_0_0_xxxyyz[k] = -g_0_xxxyyz[k] * cd_x[k] + g_0_xxxxyyz[k];

                g_x_0_0_xxxyzz[k] = -g_0_xxxyzz[k] * cd_x[k] + g_0_xxxxyzz[k];

                g_x_0_0_xxxzzz[k] = -g_0_xxxzzz[k] * cd_x[k] + g_0_xxxxzzz[k];

                g_x_0_0_xxyyyy[k] = -g_0_xxyyyy[k] * cd_x[k] + g_0_xxxyyyy[k];

                g_x_0_0_xxyyyz[k] = -g_0_xxyyyz[k] * cd_x[k] + g_0_xxxyyyz[k];

                g_x_0_0_xxyyzz[k] = -g_0_xxyyzz[k] * cd_x[k] + g_0_xxxyyzz[k];

                g_x_0_0_xxyzzz[k] = -g_0_xxyzzz[k] * cd_x[k] + g_0_xxxyzzz[k];

                g_x_0_0_xxzzzz[k] = -g_0_xxzzzz[k] * cd_x[k] + g_0_xxxzzzz[k];

                g_x_0_0_xyyyyy[k] = -g_0_xyyyyy[k] * cd_x[k] + g_0_xxyyyyy[k];

                g_x_0_0_xyyyyz[k] = -g_0_xyyyyz[k] * cd_x[k] + g_0_xxyyyyz[k];

                g_x_0_0_xyyyzz[k] = -g_0_xyyyzz[k] * cd_x[k] + g_0_xxyyyzz[k];

                g_x_0_0_xyyzzz[k] = -g_0_xyyzzz[k] * cd_x[k] + g_0_xxyyzzz[k];

                g_x_0_0_xyzzzz[k] = -g_0_xyzzzz[k] * cd_x[k] + g_0_xxyzzzz[k];

                g_x_0_0_xzzzzz[k] = -g_0_xzzzzz[k] * cd_x[k] + g_0_xxzzzzz[k];

                g_x_0_0_yyyyyy[k] = -g_0_yyyyyy[k] * cd_x[k] + g_0_xyyyyyy[k];

                g_x_0_0_yyyyyz[k] = -g_0_yyyyyz[k] * cd_x[k] + g_0_xyyyyyz[k];

                g_x_0_0_yyyyzz[k] = -g_0_yyyyzz[k] * cd_x[k] + g_0_xyyyyzz[k];

                g_x_0_0_yyyzzz[k] = -g_0_yyyzzz[k] * cd_x[k] + g_0_xyyyzzz[k];

                g_x_0_0_yyzzzz[k] = -g_0_yyzzzz[k] * cd_x[k] + g_0_xyyzzzz[k];

                g_x_0_0_yzzzzz[k] = -g_0_yzzzzz[k] * cd_x[k] + g_0_xyzzzzz[k];

                g_x_0_0_zzzzzz[k] = -g_0_zzzzzz[k] * cd_x[k] + g_0_xzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 0);

            auto g_y_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 1);

            auto g_y_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 2);

            auto g_y_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 3);

            auto g_y_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 4);

            auto g_y_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 5);

            auto g_y_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 6);

            auto g_y_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 7);

            auto g_y_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 8);

            auto g_y_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 9);

            auto g_y_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 10);

            auto g_y_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 11);

            auto g_y_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 12);

            auto g_y_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 13);

            auto g_y_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 14);

            auto g_y_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 15);

            auto g_y_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 16);

            auto g_y_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 17);

            auto g_y_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 18);

            auto g_y_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 19);

            auto g_y_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 20);

            auto g_y_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 21);

            auto g_y_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 22);

            auto g_y_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 23);

            auto g_y_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 24);

            auto g_y_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 25);

            auto g_y_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 26);

            auto g_y_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_y, g_0_xxxxxx, g_0_xxxxxxy, g_0_xxxxxy, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyz, g_0_yyyyyzz, g_0_yyyyzz, g_0_yyyyzzz, g_0_yyyzzz, g_0_yyyzzzz, g_0_yyzzzz, g_0_yyzzzzz, g_0_yzzzzz, g_0_yzzzzzz, g_0_zzzzzz, g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxxxx[k] = -g_0_xxxxxx[k] * cd_y[k] + g_0_xxxxxxy[k];

                g_y_0_0_xxxxxy[k] = -g_0_xxxxxy[k] * cd_y[k] + g_0_xxxxxyy[k];

                g_y_0_0_xxxxxz[k] = -g_0_xxxxxz[k] * cd_y[k] + g_0_xxxxxyz[k];

                g_y_0_0_xxxxyy[k] = -g_0_xxxxyy[k] * cd_y[k] + g_0_xxxxyyy[k];

                g_y_0_0_xxxxyz[k] = -g_0_xxxxyz[k] * cd_y[k] + g_0_xxxxyyz[k];

                g_y_0_0_xxxxzz[k] = -g_0_xxxxzz[k] * cd_y[k] + g_0_xxxxyzz[k];

                g_y_0_0_xxxyyy[k] = -g_0_xxxyyy[k] * cd_y[k] + g_0_xxxyyyy[k];

                g_y_0_0_xxxyyz[k] = -g_0_xxxyyz[k] * cd_y[k] + g_0_xxxyyyz[k];

                g_y_0_0_xxxyzz[k] = -g_0_xxxyzz[k] * cd_y[k] + g_0_xxxyyzz[k];

                g_y_0_0_xxxzzz[k] = -g_0_xxxzzz[k] * cd_y[k] + g_0_xxxyzzz[k];

                g_y_0_0_xxyyyy[k] = -g_0_xxyyyy[k] * cd_y[k] + g_0_xxyyyyy[k];

                g_y_0_0_xxyyyz[k] = -g_0_xxyyyz[k] * cd_y[k] + g_0_xxyyyyz[k];

                g_y_0_0_xxyyzz[k] = -g_0_xxyyzz[k] * cd_y[k] + g_0_xxyyyzz[k];

                g_y_0_0_xxyzzz[k] = -g_0_xxyzzz[k] * cd_y[k] + g_0_xxyyzzz[k];

                g_y_0_0_xxzzzz[k] = -g_0_xxzzzz[k] * cd_y[k] + g_0_xxyzzzz[k];

                g_y_0_0_xyyyyy[k] = -g_0_xyyyyy[k] * cd_y[k] + g_0_xyyyyyy[k];

                g_y_0_0_xyyyyz[k] = -g_0_xyyyyz[k] * cd_y[k] + g_0_xyyyyyz[k];

                g_y_0_0_xyyyzz[k] = -g_0_xyyyzz[k] * cd_y[k] + g_0_xyyyyzz[k];

                g_y_0_0_xyyzzz[k] = -g_0_xyyzzz[k] * cd_y[k] + g_0_xyyyzzz[k];

                g_y_0_0_xyzzzz[k] = -g_0_xyzzzz[k] * cd_y[k] + g_0_xyyzzzz[k];

                g_y_0_0_xzzzzz[k] = -g_0_xzzzzz[k] * cd_y[k] + g_0_xyzzzzz[k];

                g_y_0_0_yyyyyy[k] = -g_0_yyyyyy[k] * cd_y[k] + g_0_yyyyyyy[k];

                g_y_0_0_yyyyyz[k] = -g_0_yyyyyz[k] * cd_y[k] + g_0_yyyyyyz[k];

                g_y_0_0_yyyyzz[k] = -g_0_yyyyzz[k] * cd_y[k] + g_0_yyyyyzz[k];

                g_y_0_0_yyyzzz[k] = -g_0_yyyzzz[k] * cd_y[k] + g_0_yyyyzzz[k];

                g_y_0_0_yyzzzz[k] = -g_0_yyzzzz[k] * cd_y[k] + g_0_yyyzzzz[k];

                g_y_0_0_yzzzzz[k] = -g_0_yzzzzz[k] * cd_y[k] + g_0_yyzzzzz[k];

                g_y_0_0_zzzzzz[k] = -g_0_zzzzzz[k] * cd_y[k] + g_0_yzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 0);

            auto g_z_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 1);

            auto g_z_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 2);

            auto g_z_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 3);

            auto g_z_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 4);

            auto g_z_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 5);

            auto g_z_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 6);

            auto g_z_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 7);

            auto g_z_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 8);

            auto g_z_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 9);

            auto g_z_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 10);

            auto g_z_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 11);

            auto g_z_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 12);

            auto g_z_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 13);

            auto g_z_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 14);

            auto g_z_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 15);

            auto g_z_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 16);

            auto g_z_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 17);

            auto g_z_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 18);

            auto g_z_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 19);

            auto g_z_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 20);

            auto g_z_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 21);

            auto g_z_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 22);

            auto g_z_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 23);

            auto g_z_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 24);

            auto g_z_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 25);

            auto g_z_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 26);

            auto g_z_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_z, g_0_xxxxxx, g_0_xxxxxxz, g_0_xxxxxy, g_0_xxxxxyz, g_0_xxxxxz, g_0_xxxxxzz, g_0_xxxxyy, g_0_xxxxyyz, g_0_xxxxyz, g_0_xxxxyzz, g_0_xxxxzz, g_0_xxxxzzz, g_0_xxxyyy, g_0_xxxyyyz, g_0_xxxyyz, g_0_xxxyyzz, g_0_xxxyzz, g_0_xxxyzzz, g_0_xxxzzz, g_0_xxxzzzz, g_0_xxyyyy, g_0_xxyyyyz, g_0_xxyyyz, g_0_xxyyyzz, g_0_xxyyzz, g_0_xxyyzzz, g_0_xxyzzz, g_0_xxyzzzz, g_0_xxzzzz, g_0_xxzzzzz, g_0_xyyyyy, g_0_xyyyyyz, g_0_xyyyyz, g_0_xyyyyzz, g_0_xyyyzz, g_0_xyyyzzz, g_0_xyyzzz, g_0_xyyzzzz, g_0_xyzzzz, g_0_xyzzzzz, g_0_xzzzzz, g_0_xzzzzzz, g_0_yyyyyy, g_0_yyyyyyz, g_0_yyyyyz, g_0_yyyyyzz, g_0_yyyyzz, g_0_yyyyzzz, g_0_yyyzzz, g_0_yyyzzzz, g_0_yyzzzz, g_0_yyzzzzz, g_0_yzzzzz, g_0_yzzzzzz, g_0_zzzzzz, g_0_zzzzzzz, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxxxx[k] = -g_0_xxxxxx[k] * cd_z[k] + g_0_xxxxxxz[k];

                g_z_0_0_xxxxxy[k] = -g_0_xxxxxy[k] * cd_z[k] + g_0_xxxxxyz[k];

                g_z_0_0_xxxxxz[k] = -g_0_xxxxxz[k] * cd_z[k] + g_0_xxxxxzz[k];

                g_z_0_0_xxxxyy[k] = -g_0_xxxxyy[k] * cd_z[k] + g_0_xxxxyyz[k];

                g_z_0_0_xxxxyz[k] = -g_0_xxxxyz[k] * cd_z[k] + g_0_xxxxyzz[k];

                g_z_0_0_xxxxzz[k] = -g_0_xxxxzz[k] * cd_z[k] + g_0_xxxxzzz[k];

                g_z_0_0_xxxyyy[k] = -g_0_xxxyyy[k] * cd_z[k] + g_0_xxxyyyz[k];

                g_z_0_0_xxxyyz[k] = -g_0_xxxyyz[k] * cd_z[k] + g_0_xxxyyzz[k];

                g_z_0_0_xxxyzz[k] = -g_0_xxxyzz[k] * cd_z[k] + g_0_xxxyzzz[k];

                g_z_0_0_xxxzzz[k] = -g_0_xxxzzz[k] * cd_z[k] + g_0_xxxzzzz[k];

                g_z_0_0_xxyyyy[k] = -g_0_xxyyyy[k] * cd_z[k] + g_0_xxyyyyz[k];

                g_z_0_0_xxyyyz[k] = -g_0_xxyyyz[k] * cd_z[k] + g_0_xxyyyzz[k];

                g_z_0_0_xxyyzz[k] = -g_0_xxyyzz[k] * cd_z[k] + g_0_xxyyzzz[k];

                g_z_0_0_xxyzzz[k] = -g_0_xxyzzz[k] * cd_z[k] + g_0_xxyzzzz[k];

                g_z_0_0_xxzzzz[k] = -g_0_xxzzzz[k] * cd_z[k] + g_0_xxzzzzz[k];

                g_z_0_0_xyyyyy[k] = -g_0_xyyyyy[k] * cd_z[k] + g_0_xyyyyyz[k];

                g_z_0_0_xyyyyz[k] = -g_0_xyyyyz[k] * cd_z[k] + g_0_xyyyyzz[k];

                g_z_0_0_xyyyzz[k] = -g_0_xyyyzz[k] * cd_z[k] + g_0_xyyyzzz[k];

                g_z_0_0_xyyzzz[k] = -g_0_xyyzzz[k] * cd_z[k] + g_0_xyyzzzz[k];

                g_z_0_0_xyzzzz[k] = -g_0_xyzzzz[k] * cd_z[k] + g_0_xyzzzzz[k];

                g_z_0_0_xzzzzz[k] = -g_0_xzzzzz[k] * cd_z[k] + g_0_xzzzzzz[k];

                g_z_0_0_yyyyyy[k] = -g_0_yyyyyy[k] * cd_z[k] + g_0_yyyyyyz[k];

                g_z_0_0_yyyyyz[k] = -g_0_yyyyyz[k] * cd_z[k] + g_0_yyyyyzz[k];

                g_z_0_0_yyyyzz[k] = -g_0_yyyyzz[k] * cd_z[k] + g_0_yyyyzzz[k];

                g_z_0_0_yyyzzz[k] = -g_0_yyyzzz[k] * cd_z[k] + g_0_yyyzzzz[k];

                g_z_0_0_yyzzzz[k] = -g_0_yyzzzz[k] * cd_z[k] + g_0_yyzzzzz[k];

                g_z_0_0_yzzzzz[k] = -g_0_yzzzzz[k] * cd_z[k] + g_0_yzzzzzz[k];

                g_z_0_0_zzzzzz[k] = -g_0_zzzzzz[k] * cd_z[k] + g_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

