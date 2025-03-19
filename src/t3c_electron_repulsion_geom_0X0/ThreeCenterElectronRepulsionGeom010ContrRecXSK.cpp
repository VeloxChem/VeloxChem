#include "ThreeCenterElectronRepulsionGeom010ContrRecXSK.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xsk(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xsk,
                                        const size_t idx_xsk,
                                        const size_t idx_xsl,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SSK

        const auto sk_off = idx_xsk + i * 36;

        auto g_0_xxxxxxx = cbuffer.data(sk_off + 0);

        auto g_0_xxxxxxy = cbuffer.data(sk_off + 1);

        auto g_0_xxxxxxz = cbuffer.data(sk_off + 2);

        auto g_0_xxxxxyy = cbuffer.data(sk_off + 3);

        auto g_0_xxxxxyz = cbuffer.data(sk_off + 4);

        auto g_0_xxxxxzz = cbuffer.data(sk_off + 5);

        auto g_0_xxxxyyy = cbuffer.data(sk_off + 6);

        auto g_0_xxxxyyz = cbuffer.data(sk_off + 7);

        auto g_0_xxxxyzz = cbuffer.data(sk_off + 8);

        auto g_0_xxxxzzz = cbuffer.data(sk_off + 9);

        auto g_0_xxxyyyy = cbuffer.data(sk_off + 10);

        auto g_0_xxxyyyz = cbuffer.data(sk_off + 11);

        auto g_0_xxxyyzz = cbuffer.data(sk_off + 12);

        auto g_0_xxxyzzz = cbuffer.data(sk_off + 13);

        auto g_0_xxxzzzz = cbuffer.data(sk_off + 14);

        auto g_0_xxyyyyy = cbuffer.data(sk_off + 15);

        auto g_0_xxyyyyz = cbuffer.data(sk_off + 16);

        auto g_0_xxyyyzz = cbuffer.data(sk_off + 17);

        auto g_0_xxyyzzz = cbuffer.data(sk_off + 18);

        auto g_0_xxyzzzz = cbuffer.data(sk_off + 19);

        auto g_0_xxzzzzz = cbuffer.data(sk_off + 20);

        auto g_0_xyyyyyy = cbuffer.data(sk_off + 21);

        auto g_0_xyyyyyz = cbuffer.data(sk_off + 22);

        auto g_0_xyyyyzz = cbuffer.data(sk_off + 23);

        auto g_0_xyyyzzz = cbuffer.data(sk_off + 24);

        auto g_0_xyyzzzz = cbuffer.data(sk_off + 25);

        auto g_0_xyzzzzz = cbuffer.data(sk_off + 26);

        auto g_0_xzzzzzz = cbuffer.data(sk_off + 27);

        auto g_0_yyyyyyy = cbuffer.data(sk_off + 28);

        auto g_0_yyyyyyz = cbuffer.data(sk_off + 29);

        auto g_0_yyyyyzz = cbuffer.data(sk_off + 30);

        auto g_0_yyyyzzz = cbuffer.data(sk_off + 31);

        auto g_0_yyyzzzz = cbuffer.data(sk_off + 32);

        auto g_0_yyzzzzz = cbuffer.data(sk_off + 33);

        auto g_0_yzzzzzz = cbuffer.data(sk_off + 34);

        auto g_0_zzzzzzz = cbuffer.data(sk_off + 35);

        /// Set up components of auxilary buffer : SSL

        const auto sl_off = idx_xsl + i * 45;

        auto g_0_xxxxxxxx = cbuffer.data(sl_off + 0);

        auto g_0_xxxxxxxy = cbuffer.data(sl_off + 1);

        auto g_0_xxxxxxxz = cbuffer.data(sl_off + 2);

        auto g_0_xxxxxxyy = cbuffer.data(sl_off + 3);

        auto g_0_xxxxxxyz = cbuffer.data(sl_off + 4);

        auto g_0_xxxxxxzz = cbuffer.data(sl_off + 5);

        auto g_0_xxxxxyyy = cbuffer.data(sl_off + 6);

        auto g_0_xxxxxyyz = cbuffer.data(sl_off + 7);

        auto g_0_xxxxxyzz = cbuffer.data(sl_off + 8);

        auto g_0_xxxxxzzz = cbuffer.data(sl_off + 9);

        auto g_0_xxxxyyyy = cbuffer.data(sl_off + 10);

        auto g_0_xxxxyyyz = cbuffer.data(sl_off + 11);

        auto g_0_xxxxyyzz = cbuffer.data(sl_off + 12);

        auto g_0_xxxxyzzz = cbuffer.data(sl_off + 13);

        auto g_0_xxxxzzzz = cbuffer.data(sl_off + 14);

        auto g_0_xxxyyyyy = cbuffer.data(sl_off + 15);

        auto g_0_xxxyyyyz = cbuffer.data(sl_off + 16);

        auto g_0_xxxyyyzz = cbuffer.data(sl_off + 17);

        auto g_0_xxxyyzzz = cbuffer.data(sl_off + 18);

        auto g_0_xxxyzzzz = cbuffer.data(sl_off + 19);

        auto g_0_xxxzzzzz = cbuffer.data(sl_off + 20);

        auto g_0_xxyyyyyy = cbuffer.data(sl_off + 21);

        auto g_0_xxyyyyyz = cbuffer.data(sl_off + 22);

        auto g_0_xxyyyyzz = cbuffer.data(sl_off + 23);

        auto g_0_xxyyyzzz = cbuffer.data(sl_off + 24);

        auto g_0_xxyyzzzz = cbuffer.data(sl_off + 25);

        auto g_0_xxyzzzzz = cbuffer.data(sl_off + 26);

        auto g_0_xxzzzzzz = cbuffer.data(sl_off + 27);

        auto g_0_xyyyyyyy = cbuffer.data(sl_off + 28);

        auto g_0_xyyyyyyz = cbuffer.data(sl_off + 29);

        auto g_0_xyyyyyzz = cbuffer.data(sl_off + 30);

        auto g_0_xyyyyzzz = cbuffer.data(sl_off + 31);

        auto g_0_xyyyzzzz = cbuffer.data(sl_off + 32);

        auto g_0_xyyzzzzz = cbuffer.data(sl_off + 33);

        auto g_0_xyzzzzzz = cbuffer.data(sl_off + 34);

        auto g_0_xzzzzzzz = cbuffer.data(sl_off + 35);

        auto g_0_yyyyyyyy = cbuffer.data(sl_off + 36);

        auto g_0_yyyyyyyz = cbuffer.data(sl_off + 37);

        auto g_0_yyyyyyzz = cbuffer.data(sl_off + 38);

        auto g_0_yyyyyzzz = cbuffer.data(sl_off + 39);

        auto g_0_yyyyzzzz = cbuffer.data(sl_off + 40);

        auto g_0_yyyzzzzz = cbuffer.data(sl_off + 41);

        auto g_0_yyzzzzzz = cbuffer.data(sl_off + 42);

        auto g_0_yzzzzzzz = cbuffer.data(sl_off + 43);

        auto g_0_zzzzzzzz = cbuffer.data(sl_off + 44);

        /// set up bra offset for contr_buffer_xxsk

        const auto sk_geom_10_off = idx_geom_10_xsk + i * 36;

        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_x_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_0_xxxxxxx, g_0_xxxxxxxx, g_0_xxxxxxxy, g_0_xxxxxxxz, g_0_xxxxxxy, g_0_xxxxxxyy, g_0_xxxxxxyz, g_0_xxxxxxz, g_0_xxxxxxzz, g_0_xxxxxyy, g_0_xxxxxyyy, g_0_xxxxxyyz, g_0_xxxxxyz, g_0_xxxxxyzz, g_0_xxxxxzz, g_0_xxxxxzzz, g_0_xxxxyyy, g_0_xxxxyyyy, g_0_xxxxyyyz, g_0_xxxxyyz, g_0_xxxxyyzz, g_0_xxxxyzz, g_0_xxxxyzzz, g_0_xxxxzzz, g_0_xxxxzzzz, g_0_xxxyyyy, g_0_xxxyyyyy, g_0_xxxyyyyz, g_0_xxxyyyz, g_0_xxxyyyzz, g_0_xxxyyzz, g_0_xxxyyzzz, g_0_xxxyzzz, g_0_xxxyzzzz, g_0_xxxzzzz, g_0_xxxzzzzz, g_0_xxyyyyy, g_0_xxyyyyyy, g_0_xxyyyyyz, g_0_xxyyyyz, g_0_xxyyyyzz, g_0_xxyyyzz, g_0_xxyyyzzz, g_0_xxyyzzz, g_0_xxyyzzzz, g_0_xxyzzzz, g_0_xxyzzzzz, g_0_xxzzzzz, g_0_xxzzzzzz, g_0_xyyyyyy, g_0_xyyyyyyy, g_0_xyyyyyyz, g_0_xyyyyyz, g_0_xyyyyyzz, g_0_xyyyyzz, g_0_xyyyyzzz, g_0_xyyyzzz, g_0_xyyyzzzz, g_0_xyyzzzz, g_0_xyyzzzzz, g_0_xyzzzzz, g_0_xyzzzzzz, g_0_xzzzzzz, g_0_xzzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz, g_x_0_0_xxxxxxx, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_yyyyyyy, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_0_xxxxxxx[k] = -g_0_xxxxxxx[k] * cd_x[k] + g_0_xxxxxxxx[k];

            g_x_0_0_xxxxxxy[k] = -g_0_xxxxxxy[k] * cd_x[k] + g_0_xxxxxxxy[k];

            g_x_0_0_xxxxxxz[k] = -g_0_xxxxxxz[k] * cd_x[k] + g_0_xxxxxxxz[k];

            g_x_0_0_xxxxxyy[k] = -g_0_xxxxxyy[k] * cd_x[k] + g_0_xxxxxxyy[k];

            g_x_0_0_xxxxxyz[k] = -g_0_xxxxxyz[k] * cd_x[k] + g_0_xxxxxxyz[k];

            g_x_0_0_xxxxxzz[k] = -g_0_xxxxxzz[k] * cd_x[k] + g_0_xxxxxxzz[k];

            g_x_0_0_xxxxyyy[k] = -g_0_xxxxyyy[k] * cd_x[k] + g_0_xxxxxyyy[k];

            g_x_0_0_xxxxyyz[k] = -g_0_xxxxyyz[k] * cd_x[k] + g_0_xxxxxyyz[k];

            g_x_0_0_xxxxyzz[k] = -g_0_xxxxyzz[k] * cd_x[k] + g_0_xxxxxyzz[k];

            g_x_0_0_xxxxzzz[k] = -g_0_xxxxzzz[k] * cd_x[k] + g_0_xxxxxzzz[k];

            g_x_0_0_xxxyyyy[k] = -g_0_xxxyyyy[k] * cd_x[k] + g_0_xxxxyyyy[k];

            g_x_0_0_xxxyyyz[k] = -g_0_xxxyyyz[k] * cd_x[k] + g_0_xxxxyyyz[k];

            g_x_0_0_xxxyyzz[k] = -g_0_xxxyyzz[k] * cd_x[k] + g_0_xxxxyyzz[k];

            g_x_0_0_xxxyzzz[k] = -g_0_xxxyzzz[k] * cd_x[k] + g_0_xxxxyzzz[k];

            g_x_0_0_xxxzzzz[k] = -g_0_xxxzzzz[k] * cd_x[k] + g_0_xxxxzzzz[k];

            g_x_0_0_xxyyyyy[k] = -g_0_xxyyyyy[k] * cd_x[k] + g_0_xxxyyyyy[k];

            g_x_0_0_xxyyyyz[k] = -g_0_xxyyyyz[k] * cd_x[k] + g_0_xxxyyyyz[k];

            g_x_0_0_xxyyyzz[k] = -g_0_xxyyyzz[k] * cd_x[k] + g_0_xxxyyyzz[k];

            g_x_0_0_xxyyzzz[k] = -g_0_xxyyzzz[k] * cd_x[k] + g_0_xxxyyzzz[k];

            g_x_0_0_xxyzzzz[k] = -g_0_xxyzzzz[k] * cd_x[k] + g_0_xxxyzzzz[k];

            g_x_0_0_xxzzzzz[k] = -g_0_xxzzzzz[k] * cd_x[k] + g_0_xxxzzzzz[k];

            g_x_0_0_xyyyyyy[k] = -g_0_xyyyyyy[k] * cd_x[k] + g_0_xxyyyyyy[k];

            g_x_0_0_xyyyyyz[k] = -g_0_xyyyyyz[k] * cd_x[k] + g_0_xxyyyyyz[k];

            g_x_0_0_xyyyyzz[k] = -g_0_xyyyyzz[k] * cd_x[k] + g_0_xxyyyyzz[k];

            g_x_0_0_xyyyzzz[k] = -g_0_xyyyzzz[k] * cd_x[k] + g_0_xxyyyzzz[k];

            g_x_0_0_xyyzzzz[k] = -g_0_xyyzzzz[k] * cd_x[k] + g_0_xxyyzzzz[k];

            g_x_0_0_xyzzzzz[k] = -g_0_xyzzzzz[k] * cd_x[k] + g_0_xxyzzzzz[k];

            g_x_0_0_xzzzzzz[k] = -g_0_xzzzzzz[k] * cd_x[k] + g_0_xxzzzzzz[k];

            g_x_0_0_yyyyyyy[k] = -g_0_yyyyyyy[k] * cd_x[k] + g_0_xyyyyyyy[k];

            g_x_0_0_yyyyyyz[k] = -g_0_yyyyyyz[k] * cd_x[k] + g_0_xyyyyyyz[k];

            g_x_0_0_yyyyyzz[k] = -g_0_yyyyyzz[k] * cd_x[k] + g_0_xyyyyyzz[k];

            g_x_0_0_yyyyzzz[k] = -g_0_yyyyzzz[k] * cd_x[k] + g_0_xyyyyzzz[k];

            g_x_0_0_yyyzzzz[k] = -g_0_yyyzzzz[k] * cd_x[k] + g_0_xyyyzzzz[k];

            g_x_0_0_yyzzzzz[k] = -g_0_yyzzzzz[k] * cd_x[k] + g_0_xyyzzzzz[k];

            g_x_0_0_yzzzzzz[k] = -g_0_yzzzzzz[k] * cd_x[k] + g_0_xyzzzzzz[k];

            g_x_0_0_zzzzzzz[k] = -g_0_zzzzzzz[k] * cd_x[k] + g_0_xzzzzzzz[k];
        }
        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_y_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 36 * acomps  + 0);

        auto g_y_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 1);

        auto g_y_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 2);

        auto g_y_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 3);

        auto g_y_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 4);

        auto g_y_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 5);

        auto g_y_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 6);

        auto g_y_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 7);

        auto g_y_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 8);

        auto g_y_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 9);

        auto g_y_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 10);

        auto g_y_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 11);

        auto g_y_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 12);

        auto g_y_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 13);

        auto g_y_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 14);

        auto g_y_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 15);

        auto g_y_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 16);

        auto g_y_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 17);

        auto g_y_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 18);

        auto g_y_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 19);

        auto g_y_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 20);

        auto g_y_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 21);

        auto g_y_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 22);

        auto g_y_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 23);

        auto g_y_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 24);

        auto g_y_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 25);

        auto g_y_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 26);

        auto g_y_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 27);

        auto g_y_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps  + 28);

        auto g_y_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 29);

        auto g_y_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 30);

        auto g_y_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 31);

        auto g_y_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 32);

        auto g_y_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 33);

        auto g_y_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 34);

        auto g_y_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps  + 35);

        #pragma omp simd aligned(cd_y, g_0_xxxxxxx, g_0_xxxxxxxy, g_0_xxxxxxy, g_0_xxxxxxyy, g_0_xxxxxxyz, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyyy, g_0_xxxxxyyz, g_0_xxxxxyz, g_0_xxxxxyzz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyyy, g_0_xxxxyyyz, g_0_xxxxyyz, g_0_xxxxyyzz, g_0_xxxxyzz, g_0_xxxxyzzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyyy, g_0_xxxyyyyz, g_0_xxxyyyz, g_0_xxxyyyzz, g_0_xxxyyzz, g_0_xxxyyzzz, g_0_xxxyzzz, g_0_xxxyzzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyyy, g_0_xxyyyyyz, g_0_xxyyyyz, g_0_xxyyyyzz, g_0_xxyyyzz, g_0_xxyyyzzz, g_0_xxyyzzz, g_0_xxyyzzzz, g_0_xxyzzzz, g_0_xxyzzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyyy, g_0_xyyyyyyz, g_0_xyyyyyz, g_0_xyyyyyzz, g_0_xyyyyzz, g_0_xyyyyzzz, g_0_xyyyzzz, g_0_xyyyzzzz, g_0_xyyzzzz, g_0_xyyzzzzz, g_0_xyzzzzz, g_0_xyzzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyyy, g_0_yyyyyyyz, g_0_yyyyyyz, g_0_yyyyyyzz, g_0_yyyyyzz, g_0_yyyyyzzz, g_0_yyyyzzz, g_0_yyyyzzzz, g_0_yyyzzzz, g_0_yyyzzzzz, g_0_yyzzzzz, g_0_yyzzzzzz, g_0_yzzzzzz, g_0_yzzzzzzz, g_0_zzzzzzz, g_y_0_0_xxxxxxx, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_yyyyyyy, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_0_xxxxxxx[k] = -g_0_xxxxxxx[k] * cd_y[k] + g_0_xxxxxxxy[k];

            g_y_0_0_xxxxxxy[k] = -g_0_xxxxxxy[k] * cd_y[k] + g_0_xxxxxxyy[k];

            g_y_0_0_xxxxxxz[k] = -g_0_xxxxxxz[k] * cd_y[k] + g_0_xxxxxxyz[k];

            g_y_0_0_xxxxxyy[k] = -g_0_xxxxxyy[k] * cd_y[k] + g_0_xxxxxyyy[k];

            g_y_0_0_xxxxxyz[k] = -g_0_xxxxxyz[k] * cd_y[k] + g_0_xxxxxyyz[k];

            g_y_0_0_xxxxxzz[k] = -g_0_xxxxxzz[k] * cd_y[k] + g_0_xxxxxyzz[k];

            g_y_0_0_xxxxyyy[k] = -g_0_xxxxyyy[k] * cd_y[k] + g_0_xxxxyyyy[k];

            g_y_0_0_xxxxyyz[k] = -g_0_xxxxyyz[k] * cd_y[k] + g_0_xxxxyyyz[k];

            g_y_0_0_xxxxyzz[k] = -g_0_xxxxyzz[k] * cd_y[k] + g_0_xxxxyyzz[k];

            g_y_0_0_xxxxzzz[k] = -g_0_xxxxzzz[k] * cd_y[k] + g_0_xxxxyzzz[k];

            g_y_0_0_xxxyyyy[k] = -g_0_xxxyyyy[k] * cd_y[k] + g_0_xxxyyyyy[k];

            g_y_0_0_xxxyyyz[k] = -g_0_xxxyyyz[k] * cd_y[k] + g_0_xxxyyyyz[k];

            g_y_0_0_xxxyyzz[k] = -g_0_xxxyyzz[k] * cd_y[k] + g_0_xxxyyyzz[k];

            g_y_0_0_xxxyzzz[k] = -g_0_xxxyzzz[k] * cd_y[k] + g_0_xxxyyzzz[k];

            g_y_0_0_xxxzzzz[k] = -g_0_xxxzzzz[k] * cd_y[k] + g_0_xxxyzzzz[k];

            g_y_0_0_xxyyyyy[k] = -g_0_xxyyyyy[k] * cd_y[k] + g_0_xxyyyyyy[k];

            g_y_0_0_xxyyyyz[k] = -g_0_xxyyyyz[k] * cd_y[k] + g_0_xxyyyyyz[k];

            g_y_0_0_xxyyyzz[k] = -g_0_xxyyyzz[k] * cd_y[k] + g_0_xxyyyyzz[k];

            g_y_0_0_xxyyzzz[k] = -g_0_xxyyzzz[k] * cd_y[k] + g_0_xxyyyzzz[k];

            g_y_0_0_xxyzzzz[k] = -g_0_xxyzzzz[k] * cd_y[k] + g_0_xxyyzzzz[k];

            g_y_0_0_xxzzzzz[k] = -g_0_xxzzzzz[k] * cd_y[k] + g_0_xxyzzzzz[k];

            g_y_0_0_xyyyyyy[k] = -g_0_xyyyyyy[k] * cd_y[k] + g_0_xyyyyyyy[k];

            g_y_0_0_xyyyyyz[k] = -g_0_xyyyyyz[k] * cd_y[k] + g_0_xyyyyyyz[k];

            g_y_0_0_xyyyyzz[k] = -g_0_xyyyyzz[k] * cd_y[k] + g_0_xyyyyyzz[k];

            g_y_0_0_xyyyzzz[k] = -g_0_xyyyzzz[k] * cd_y[k] + g_0_xyyyyzzz[k];

            g_y_0_0_xyyzzzz[k] = -g_0_xyyzzzz[k] * cd_y[k] + g_0_xyyyzzzz[k];

            g_y_0_0_xyzzzzz[k] = -g_0_xyzzzzz[k] * cd_y[k] + g_0_xyyzzzzz[k];

            g_y_0_0_xzzzzzz[k] = -g_0_xzzzzzz[k] * cd_y[k] + g_0_xyzzzzzz[k];

            g_y_0_0_yyyyyyy[k] = -g_0_yyyyyyy[k] * cd_y[k] + g_0_yyyyyyyy[k];

            g_y_0_0_yyyyyyz[k] = -g_0_yyyyyyz[k] * cd_y[k] + g_0_yyyyyyyz[k];

            g_y_0_0_yyyyyzz[k] = -g_0_yyyyyzz[k] * cd_y[k] + g_0_yyyyyyzz[k];

            g_y_0_0_yyyyzzz[k] = -g_0_yyyyzzz[k] * cd_y[k] + g_0_yyyyyzzz[k];

            g_y_0_0_yyyzzzz[k] = -g_0_yyyzzzz[k] * cd_y[k] + g_0_yyyyzzzz[k];

            g_y_0_0_yyzzzzz[k] = -g_0_yyzzzzz[k] * cd_y[k] + g_0_yyyzzzzz[k];

            g_y_0_0_yzzzzzz[k] = -g_0_yzzzzzz[k] * cd_y[k] + g_0_yyzzzzzz[k];

            g_y_0_0_zzzzzzz[k] = -g_0_zzzzzzz[k] * cd_y[k] + g_0_yzzzzzzz[k];
        }
        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_z_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 72 * acomps  + 0);

        auto g_z_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 1);

        auto g_z_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 2);

        auto g_z_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 3);

        auto g_z_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 4);

        auto g_z_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 5);

        auto g_z_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 6);

        auto g_z_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 7);

        auto g_z_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 8);

        auto g_z_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 9);

        auto g_z_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 10);

        auto g_z_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 11);

        auto g_z_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 12);

        auto g_z_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 13);

        auto g_z_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 14);

        auto g_z_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 15);

        auto g_z_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 16);

        auto g_z_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 17);

        auto g_z_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 18);

        auto g_z_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 19);

        auto g_z_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 20);

        auto g_z_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 21);

        auto g_z_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 22);

        auto g_z_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 23);

        auto g_z_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 24);

        auto g_z_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 25);

        auto g_z_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 26);

        auto g_z_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 27);

        auto g_z_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps  + 28);

        auto g_z_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 29);

        auto g_z_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 30);

        auto g_z_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 31);

        auto g_z_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 32);

        auto g_z_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 33);

        auto g_z_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 34);

        auto g_z_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_0_xxxxxxx, g_0_xxxxxxxz, g_0_xxxxxxy, g_0_xxxxxxyz, g_0_xxxxxxz, g_0_xxxxxxzz, g_0_xxxxxyy, g_0_xxxxxyyz, g_0_xxxxxyz, g_0_xxxxxyzz, g_0_xxxxxzz, g_0_xxxxxzzz, g_0_xxxxyyy, g_0_xxxxyyyz, g_0_xxxxyyz, g_0_xxxxyyzz, g_0_xxxxyzz, g_0_xxxxyzzz, g_0_xxxxzzz, g_0_xxxxzzzz, g_0_xxxyyyy, g_0_xxxyyyyz, g_0_xxxyyyz, g_0_xxxyyyzz, g_0_xxxyyzz, g_0_xxxyyzzz, g_0_xxxyzzz, g_0_xxxyzzzz, g_0_xxxzzzz, g_0_xxxzzzzz, g_0_xxyyyyy, g_0_xxyyyyyz, g_0_xxyyyyz, g_0_xxyyyyzz, g_0_xxyyyzz, g_0_xxyyyzzz, g_0_xxyyzzz, g_0_xxyyzzzz, g_0_xxyzzzz, g_0_xxyzzzzz, g_0_xxzzzzz, g_0_xxzzzzzz, g_0_xyyyyyy, g_0_xyyyyyyz, g_0_xyyyyyz, g_0_xyyyyyzz, g_0_xyyyyzz, g_0_xyyyyzzz, g_0_xyyyzzz, g_0_xyyyzzzz, g_0_xyyzzzz, g_0_xyyzzzzz, g_0_xyzzzzz, g_0_xyzzzzzz, g_0_xzzzzzz, g_0_xzzzzzzz, g_0_yyyyyyy, g_0_yyyyyyyz, g_0_yyyyyyz, g_0_yyyyyyzz, g_0_yyyyyzz, g_0_yyyyyzzz, g_0_yyyyzzz, g_0_yyyyzzzz, g_0_yyyzzzz, g_0_yyyzzzzz, g_0_yyzzzzz, g_0_yyzzzzzz, g_0_yzzzzzz, g_0_yzzzzzzz, g_0_zzzzzzz, g_0_zzzzzzzz, g_z_0_0_xxxxxxx, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_yyyyyyy, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_0_xxxxxxx[k] = -g_0_xxxxxxx[k] * cd_z[k] + g_0_xxxxxxxz[k];

            g_z_0_0_xxxxxxy[k] = -g_0_xxxxxxy[k] * cd_z[k] + g_0_xxxxxxyz[k];

            g_z_0_0_xxxxxxz[k] = -g_0_xxxxxxz[k] * cd_z[k] + g_0_xxxxxxzz[k];

            g_z_0_0_xxxxxyy[k] = -g_0_xxxxxyy[k] * cd_z[k] + g_0_xxxxxyyz[k];

            g_z_0_0_xxxxxyz[k] = -g_0_xxxxxyz[k] * cd_z[k] + g_0_xxxxxyzz[k];

            g_z_0_0_xxxxxzz[k] = -g_0_xxxxxzz[k] * cd_z[k] + g_0_xxxxxzzz[k];

            g_z_0_0_xxxxyyy[k] = -g_0_xxxxyyy[k] * cd_z[k] + g_0_xxxxyyyz[k];

            g_z_0_0_xxxxyyz[k] = -g_0_xxxxyyz[k] * cd_z[k] + g_0_xxxxyyzz[k];

            g_z_0_0_xxxxyzz[k] = -g_0_xxxxyzz[k] * cd_z[k] + g_0_xxxxyzzz[k];

            g_z_0_0_xxxxzzz[k] = -g_0_xxxxzzz[k] * cd_z[k] + g_0_xxxxzzzz[k];

            g_z_0_0_xxxyyyy[k] = -g_0_xxxyyyy[k] * cd_z[k] + g_0_xxxyyyyz[k];

            g_z_0_0_xxxyyyz[k] = -g_0_xxxyyyz[k] * cd_z[k] + g_0_xxxyyyzz[k];

            g_z_0_0_xxxyyzz[k] = -g_0_xxxyyzz[k] * cd_z[k] + g_0_xxxyyzzz[k];

            g_z_0_0_xxxyzzz[k] = -g_0_xxxyzzz[k] * cd_z[k] + g_0_xxxyzzzz[k];

            g_z_0_0_xxxzzzz[k] = -g_0_xxxzzzz[k] * cd_z[k] + g_0_xxxzzzzz[k];

            g_z_0_0_xxyyyyy[k] = -g_0_xxyyyyy[k] * cd_z[k] + g_0_xxyyyyyz[k];

            g_z_0_0_xxyyyyz[k] = -g_0_xxyyyyz[k] * cd_z[k] + g_0_xxyyyyzz[k];

            g_z_0_0_xxyyyzz[k] = -g_0_xxyyyzz[k] * cd_z[k] + g_0_xxyyyzzz[k];

            g_z_0_0_xxyyzzz[k] = -g_0_xxyyzzz[k] * cd_z[k] + g_0_xxyyzzzz[k];

            g_z_0_0_xxyzzzz[k] = -g_0_xxyzzzz[k] * cd_z[k] + g_0_xxyzzzzz[k];

            g_z_0_0_xxzzzzz[k] = -g_0_xxzzzzz[k] * cd_z[k] + g_0_xxzzzzzz[k];

            g_z_0_0_xyyyyyy[k] = -g_0_xyyyyyy[k] * cd_z[k] + g_0_xyyyyyyz[k];

            g_z_0_0_xyyyyyz[k] = -g_0_xyyyyyz[k] * cd_z[k] + g_0_xyyyyyzz[k];

            g_z_0_0_xyyyyzz[k] = -g_0_xyyyyzz[k] * cd_z[k] + g_0_xyyyyzzz[k];

            g_z_0_0_xyyyzzz[k] = -g_0_xyyyzzz[k] * cd_z[k] + g_0_xyyyzzzz[k];

            g_z_0_0_xyyzzzz[k] = -g_0_xyyzzzz[k] * cd_z[k] + g_0_xyyzzzzz[k];

            g_z_0_0_xyzzzzz[k] = -g_0_xyzzzzz[k] * cd_z[k] + g_0_xyzzzzzz[k];

            g_z_0_0_xzzzzzz[k] = -g_0_xzzzzzz[k] * cd_z[k] + g_0_xzzzzzzz[k];

            g_z_0_0_yyyyyyy[k] = -g_0_yyyyyyy[k] * cd_z[k] + g_0_yyyyyyyz[k];

            g_z_0_0_yyyyyyz[k] = -g_0_yyyyyyz[k] * cd_z[k] + g_0_yyyyyyzz[k];

            g_z_0_0_yyyyyzz[k] = -g_0_yyyyyzz[k] * cd_z[k] + g_0_yyyyyzzz[k];

            g_z_0_0_yyyyzzz[k] = -g_0_yyyyzzz[k] * cd_z[k] + g_0_yyyyzzzz[k];

            g_z_0_0_yyyzzzz[k] = -g_0_yyyzzzz[k] * cd_z[k] + g_0_yyyzzzzz[k];

            g_z_0_0_yyzzzzz[k] = -g_0_yyzzzzz[k] * cd_z[k] + g_0_yyzzzzzz[k];

            g_z_0_0_yzzzzzz[k] = -g_0_yzzzzzz[k] * cd_z[k] + g_0_yzzzzzzz[k];

            g_z_0_0_zzzzzzz[k] = -g_0_zzzzzzz[k] * cd_z[k] + g_0_zzzzzzzz[k];
        }
    }
}

} // t3ceri namespace

