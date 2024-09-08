#include "ElectronRepulsionContrRecPKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pkxx(CSimdArray<double>&   cbuffer,
                                     const size_t          idx_pkxx,
                                     const size_t          idx_skxx,
                                     const size_t          idx_slxx,
                                     const TPoint<double>& r_ab,
                                     const int             c_angmom,
                                     const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
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

            /// set up bra offset for contr_buffer_pkxx

            const auto pk_off = idx_pkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

            auto g_x_xxxxxxx = cbuffer.data(pk_off + 0 * ccomps * dcomps);

            auto g_x_xxxxxxy = cbuffer.data(pk_off + 1 * ccomps * dcomps);

            auto g_x_xxxxxxz = cbuffer.data(pk_off + 2 * ccomps * dcomps);

            auto g_x_xxxxxyy = cbuffer.data(pk_off + 3 * ccomps * dcomps);

            auto g_x_xxxxxyz = cbuffer.data(pk_off + 4 * ccomps * dcomps);

            auto g_x_xxxxxzz = cbuffer.data(pk_off + 5 * ccomps * dcomps);

            auto g_x_xxxxyyy = cbuffer.data(pk_off + 6 * ccomps * dcomps);

            auto g_x_xxxxyyz = cbuffer.data(pk_off + 7 * ccomps * dcomps);

            auto g_x_xxxxyzz = cbuffer.data(pk_off + 8 * ccomps * dcomps);

            auto g_x_xxxxzzz = cbuffer.data(pk_off + 9 * ccomps * dcomps);

            auto g_x_xxxyyyy = cbuffer.data(pk_off + 10 * ccomps * dcomps);

            auto g_x_xxxyyyz = cbuffer.data(pk_off + 11 * ccomps * dcomps);

            auto g_x_xxxyyzz = cbuffer.data(pk_off + 12 * ccomps * dcomps);

            auto g_x_xxxyzzz = cbuffer.data(pk_off + 13 * ccomps * dcomps);

            auto g_x_xxxzzzz = cbuffer.data(pk_off + 14 * ccomps * dcomps);

            auto g_x_xxyyyyy = cbuffer.data(pk_off + 15 * ccomps * dcomps);

            auto g_x_xxyyyyz = cbuffer.data(pk_off + 16 * ccomps * dcomps);

            auto g_x_xxyyyzz = cbuffer.data(pk_off + 17 * ccomps * dcomps);

            auto g_x_xxyyzzz = cbuffer.data(pk_off + 18 * ccomps * dcomps);

            auto g_x_xxyzzzz = cbuffer.data(pk_off + 19 * ccomps * dcomps);

            auto g_x_xxzzzzz = cbuffer.data(pk_off + 20 * ccomps * dcomps);

            auto g_x_xyyyyyy = cbuffer.data(pk_off + 21 * ccomps * dcomps);

            auto g_x_xyyyyyz = cbuffer.data(pk_off + 22 * ccomps * dcomps);

            auto g_x_xyyyyzz = cbuffer.data(pk_off + 23 * ccomps * dcomps);

            auto g_x_xyyyzzz = cbuffer.data(pk_off + 24 * ccomps * dcomps);

            auto g_x_xyyzzzz = cbuffer.data(pk_off + 25 * ccomps * dcomps);

            auto g_x_xyzzzzz = cbuffer.data(pk_off + 26 * ccomps * dcomps);

            auto g_x_xzzzzzz = cbuffer.data(pk_off + 27 * ccomps * dcomps);

            auto g_x_yyyyyyy = cbuffer.data(pk_off + 28 * ccomps * dcomps);

            auto g_x_yyyyyyz = cbuffer.data(pk_off + 29 * ccomps * dcomps);

            auto g_x_yyyyyzz = cbuffer.data(pk_off + 30 * ccomps * dcomps);

            auto g_x_yyyyzzz = cbuffer.data(pk_off + 31 * ccomps * dcomps);

            auto g_x_yyyzzzz = cbuffer.data(pk_off + 32 * ccomps * dcomps);

            auto g_x_yyzzzzz = cbuffer.data(pk_off + 33 * ccomps * dcomps);

            auto g_x_yzzzzzz = cbuffer.data(pk_off + 34 * ccomps * dcomps);

            auto g_x_zzzzzzz = cbuffer.data(pk_off + 35 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxxx,      \
                             g_0_xxxxxxxx, \
                             g_0_xxxxxxxy, \
                             g_0_xxxxxxxz, \
                             g_0_xxxxxxy,  \
                             g_0_xxxxxxyy, \
                             g_0_xxxxxxyz, \
                             g_0_xxxxxxz,  \
                             g_0_xxxxxxzz, \
                             g_0_xxxxxyy,  \
                             g_0_xxxxxyyy, \
                             g_0_xxxxxyyz, \
                             g_0_xxxxxyz,  \
                             g_0_xxxxxyzz, \
                             g_0_xxxxxzz,  \
                             g_0_xxxxxzzz, \
                             g_0_xxxxyyy,  \
                             g_0_xxxxyyyy, \
                             g_0_xxxxyyyz, \
                             g_0_xxxxyyz,  \
                             g_0_xxxxyyzz, \
                             g_0_xxxxyzz,  \
                             g_0_xxxxyzzz, \
                             g_0_xxxxzzz,  \
                             g_0_xxxxzzzz, \
                             g_0_xxxyyyy,  \
                             g_0_xxxyyyyy, \
                             g_0_xxxyyyyz, \
                             g_0_xxxyyyz,  \
                             g_0_xxxyyyzz, \
                             g_0_xxxyyzz,  \
                             g_0_xxxyyzzz, \
                             g_0_xxxyzzz,  \
                             g_0_xxxyzzzz, \
                             g_0_xxxzzzz,  \
                             g_0_xxxzzzzz, \
                             g_0_xxyyyyy,  \
                             g_0_xxyyyyyy, \
                             g_0_xxyyyyyz, \
                             g_0_xxyyyyz,  \
                             g_0_xxyyyyzz, \
                             g_0_xxyyyzz,  \
                             g_0_xxyyyzzz, \
                             g_0_xxyyzzz,  \
                             g_0_xxyyzzzz, \
                             g_0_xxyzzzz,  \
                             g_0_xxyzzzzz, \
                             g_0_xxzzzzz,  \
                             g_0_xxzzzzzz, \
                             g_0_xyyyyyy,  \
                             g_0_xyyyyyyy, \
                             g_0_xyyyyyyz, \
                             g_0_xyyyyyz,  \
                             g_0_xyyyyyzz, \
                             g_0_xyyyyzz,  \
                             g_0_xyyyyzzz, \
                             g_0_xyyyzzz,  \
                             g_0_xyyyzzzz, \
                             g_0_xyyzzzz,  \
                             g_0_xyyzzzzz, \
                             g_0_xyzzzzz,  \
                             g_0_xyzzzzzz, \
                             g_0_xzzzzzz,  \
                             g_0_xzzzzzzz, \
                             g_0_yyyyyyy,  \
                             g_0_yyyyyyz,  \
                             g_0_yyyyyzz,  \
                             g_0_yyyyzzz,  \
                             g_0_yyyzzzz,  \
                             g_0_yyzzzzz,  \
                             g_0_yzzzzzz,  \
                             g_0_zzzzzzz,  \
                             g_x_xxxxxxx,  \
                             g_x_xxxxxxy,  \
                             g_x_xxxxxxz,  \
                             g_x_xxxxxyy,  \
                             g_x_xxxxxyz,  \
                             g_x_xxxxxzz,  \
                             g_x_xxxxyyy,  \
                             g_x_xxxxyyz,  \
                             g_x_xxxxyzz,  \
                             g_x_xxxxzzz,  \
                             g_x_xxxyyyy,  \
                             g_x_xxxyyyz,  \
                             g_x_xxxyyzz,  \
                             g_x_xxxyzzz,  \
                             g_x_xxxzzzz,  \
                             g_x_xxyyyyy,  \
                             g_x_xxyyyyz,  \
                             g_x_xxyyyzz,  \
                             g_x_xxyyzzz,  \
                             g_x_xxyzzzz,  \
                             g_x_xxzzzzz,  \
                             g_x_xyyyyyy,  \
                             g_x_xyyyyyz,  \
                             g_x_xyyyyzz,  \
                             g_x_xyyyzzz,  \
                             g_x_xyyzzzz,  \
                             g_x_xyzzzzz,  \
                             g_x_xzzzzzz,  \
                             g_x_yyyyyyy,  \
                             g_x_yyyyyyz,  \
                             g_x_yyyyyzz,  \
                             g_x_yyyyzzz,  \
                             g_x_yyyzzzz,  \
                             g_x_yyzzzzz,  \
                             g_x_yzzzzzz,  \
                             g_x_zzzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xxxxxxx[k] = -g_0_xxxxxxx[k] * ab_x + g_0_xxxxxxxx[k];

                g_x_xxxxxxy[k] = -g_0_xxxxxxy[k] * ab_x + g_0_xxxxxxxy[k];

                g_x_xxxxxxz[k] = -g_0_xxxxxxz[k] * ab_x + g_0_xxxxxxxz[k];

                g_x_xxxxxyy[k] = -g_0_xxxxxyy[k] * ab_x + g_0_xxxxxxyy[k];

                g_x_xxxxxyz[k] = -g_0_xxxxxyz[k] * ab_x + g_0_xxxxxxyz[k];

                g_x_xxxxxzz[k] = -g_0_xxxxxzz[k] * ab_x + g_0_xxxxxxzz[k];

                g_x_xxxxyyy[k] = -g_0_xxxxyyy[k] * ab_x + g_0_xxxxxyyy[k];

                g_x_xxxxyyz[k] = -g_0_xxxxyyz[k] * ab_x + g_0_xxxxxyyz[k];

                g_x_xxxxyzz[k] = -g_0_xxxxyzz[k] * ab_x + g_0_xxxxxyzz[k];

                g_x_xxxxzzz[k] = -g_0_xxxxzzz[k] * ab_x + g_0_xxxxxzzz[k];

                g_x_xxxyyyy[k] = -g_0_xxxyyyy[k] * ab_x + g_0_xxxxyyyy[k];

                g_x_xxxyyyz[k] = -g_0_xxxyyyz[k] * ab_x + g_0_xxxxyyyz[k];

                g_x_xxxyyzz[k] = -g_0_xxxyyzz[k] * ab_x + g_0_xxxxyyzz[k];

                g_x_xxxyzzz[k] = -g_0_xxxyzzz[k] * ab_x + g_0_xxxxyzzz[k];

                g_x_xxxzzzz[k] = -g_0_xxxzzzz[k] * ab_x + g_0_xxxxzzzz[k];

                g_x_xxyyyyy[k] = -g_0_xxyyyyy[k] * ab_x + g_0_xxxyyyyy[k];

                g_x_xxyyyyz[k] = -g_0_xxyyyyz[k] * ab_x + g_0_xxxyyyyz[k];

                g_x_xxyyyzz[k] = -g_0_xxyyyzz[k] * ab_x + g_0_xxxyyyzz[k];

                g_x_xxyyzzz[k] = -g_0_xxyyzzz[k] * ab_x + g_0_xxxyyzzz[k];

                g_x_xxyzzzz[k] = -g_0_xxyzzzz[k] * ab_x + g_0_xxxyzzzz[k];

                g_x_xxzzzzz[k] = -g_0_xxzzzzz[k] * ab_x + g_0_xxxzzzzz[k];

                g_x_xyyyyyy[k] = -g_0_xyyyyyy[k] * ab_x + g_0_xxyyyyyy[k];

                g_x_xyyyyyz[k] = -g_0_xyyyyyz[k] * ab_x + g_0_xxyyyyyz[k];

                g_x_xyyyyzz[k] = -g_0_xyyyyzz[k] * ab_x + g_0_xxyyyyzz[k];

                g_x_xyyyzzz[k] = -g_0_xyyyzzz[k] * ab_x + g_0_xxyyyzzz[k];

                g_x_xyyzzzz[k] = -g_0_xyyzzzz[k] * ab_x + g_0_xxyyzzzz[k];

                g_x_xyzzzzz[k] = -g_0_xyzzzzz[k] * ab_x + g_0_xxyzzzzz[k];

                g_x_xzzzzzz[k] = -g_0_xzzzzzz[k] * ab_x + g_0_xxzzzzzz[k];

                g_x_yyyyyyy[k] = -g_0_yyyyyyy[k] * ab_x + g_0_xyyyyyyy[k];

                g_x_yyyyyyz[k] = -g_0_yyyyyyz[k] * ab_x + g_0_xyyyyyyz[k];

                g_x_yyyyyzz[k] = -g_0_yyyyyzz[k] * ab_x + g_0_xyyyyyzz[k];

                g_x_yyyyzzz[k] = -g_0_yyyyzzz[k] * ab_x + g_0_xyyyyzzz[k];

                g_x_yyyzzzz[k] = -g_0_yyyzzzz[k] * ab_x + g_0_xyyyzzzz[k];

                g_x_yyzzzzz[k] = -g_0_yyzzzzz[k] * ab_x + g_0_xyyzzzzz[k];

                g_x_yzzzzzz[k] = -g_0_yzzzzzz[k] * ab_x + g_0_xyzzzzzz[k];

                g_x_zzzzzzz[k] = -g_0_zzzzzzz[k] * ab_x + g_0_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

            auto g_y_xxxxxxx = cbuffer.data(pk_off + 36 * ccomps * dcomps);

            auto g_y_xxxxxxy = cbuffer.data(pk_off + 37 * ccomps * dcomps);

            auto g_y_xxxxxxz = cbuffer.data(pk_off + 38 * ccomps * dcomps);

            auto g_y_xxxxxyy = cbuffer.data(pk_off + 39 * ccomps * dcomps);

            auto g_y_xxxxxyz = cbuffer.data(pk_off + 40 * ccomps * dcomps);

            auto g_y_xxxxxzz = cbuffer.data(pk_off + 41 * ccomps * dcomps);

            auto g_y_xxxxyyy = cbuffer.data(pk_off + 42 * ccomps * dcomps);

            auto g_y_xxxxyyz = cbuffer.data(pk_off + 43 * ccomps * dcomps);

            auto g_y_xxxxyzz = cbuffer.data(pk_off + 44 * ccomps * dcomps);

            auto g_y_xxxxzzz = cbuffer.data(pk_off + 45 * ccomps * dcomps);

            auto g_y_xxxyyyy = cbuffer.data(pk_off + 46 * ccomps * dcomps);

            auto g_y_xxxyyyz = cbuffer.data(pk_off + 47 * ccomps * dcomps);

            auto g_y_xxxyyzz = cbuffer.data(pk_off + 48 * ccomps * dcomps);

            auto g_y_xxxyzzz = cbuffer.data(pk_off + 49 * ccomps * dcomps);

            auto g_y_xxxzzzz = cbuffer.data(pk_off + 50 * ccomps * dcomps);

            auto g_y_xxyyyyy = cbuffer.data(pk_off + 51 * ccomps * dcomps);

            auto g_y_xxyyyyz = cbuffer.data(pk_off + 52 * ccomps * dcomps);

            auto g_y_xxyyyzz = cbuffer.data(pk_off + 53 * ccomps * dcomps);

            auto g_y_xxyyzzz = cbuffer.data(pk_off + 54 * ccomps * dcomps);

            auto g_y_xxyzzzz = cbuffer.data(pk_off + 55 * ccomps * dcomps);

            auto g_y_xxzzzzz = cbuffer.data(pk_off + 56 * ccomps * dcomps);

            auto g_y_xyyyyyy = cbuffer.data(pk_off + 57 * ccomps * dcomps);

            auto g_y_xyyyyyz = cbuffer.data(pk_off + 58 * ccomps * dcomps);

            auto g_y_xyyyyzz = cbuffer.data(pk_off + 59 * ccomps * dcomps);

            auto g_y_xyyyzzz = cbuffer.data(pk_off + 60 * ccomps * dcomps);

            auto g_y_xyyzzzz = cbuffer.data(pk_off + 61 * ccomps * dcomps);

            auto g_y_xyzzzzz = cbuffer.data(pk_off + 62 * ccomps * dcomps);

            auto g_y_xzzzzzz = cbuffer.data(pk_off + 63 * ccomps * dcomps);

            auto g_y_yyyyyyy = cbuffer.data(pk_off + 64 * ccomps * dcomps);

            auto g_y_yyyyyyz = cbuffer.data(pk_off + 65 * ccomps * dcomps);

            auto g_y_yyyyyzz = cbuffer.data(pk_off + 66 * ccomps * dcomps);

            auto g_y_yyyyzzz = cbuffer.data(pk_off + 67 * ccomps * dcomps);

            auto g_y_yyyzzzz = cbuffer.data(pk_off + 68 * ccomps * dcomps);

            auto g_y_yyzzzzz = cbuffer.data(pk_off + 69 * ccomps * dcomps);

            auto g_y_yzzzzzz = cbuffer.data(pk_off + 70 * ccomps * dcomps);

            auto g_y_zzzzzzz = cbuffer.data(pk_off + 71 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxxx,      \
                             g_0_xxxxxxxy, \
                             g_0_xxxxxxy,  \
                             g_0_xxxxxxyy, \
                             g_0_xxxxxxyz, \
                             g_0_xxxxxxz,  \
                             g_0_xxxxxyy,  \
                             g_0_xxxxxyyy, \
                             g_0_xxxxxyyz, \
                             g_0_xxxxxyz,  \
                             g_0_xxxxxyzz, \
                             g_0_xxxxxzz,  \
                             g_0_xxxxyyy,  \
                             g_0_xxxxyyyy, \
                             g_0_xxxxyyyz, \
                             g_0_xxxxyyz,  \
                             g_0_xxxxyyzz, \
                             g_0_xxxxyzz,  \
                             g_0_xxxxyzzz, \
                             g_0_xxxxzzz,  \
                             g_0_xxxyyyy,  \
                             g_0_xxxyyyyy, \
                             g_0_xxxyyyyz, \
                             g_0_xxxyyyz,  \
                             g_0_xxxyyyzz, \
                             g_0_xxxyyzz,  \
                             g_0_xxxyyzzz, \
                             g_0_xxxyzzz,  \
                             g_0_xxxyzzzz, \
                             g_0_xxxzzzz,  \
                             g_0_xxyyyyy,  \
                             g_0_xxyyyyyy, \
                             g_0_xxyyyyyz, \
                             g_0_xxyyyyz,  \
                             g_0_xxyyyyzz, \
                             g_0_xxyyyzz,  \
                             g_0_xxyyyzzz, \
                             g_0_xxyyzzz,  \
                             g_0_xxyyzzzz, \
                             g_0_xxyzzzz,  \
                             g_0_xxyzzzzz, \
                             g_0_xxzzzzz,  \
                             g_0_xyyyyyy,  \
                             g_0_xyyyyyyy, \
                             g_0_xyyyyyyz, \
                             g_0_xyyyyyz,  \
                             g_0_xyyyyyzz, \
                             g_0_xyyyyzz,  \
                             g_0_xyyyyzzz, \
                             g_0_xyyyzzz,  \
                             g_0_xyyyzzzz, \
                             g_0_xyyzzzz,  \
                             g_0_xyyzzzzz, \
                             g_0_xyzzzzz,  \
                             g_0_xyzzzzzz, \
                             g_0_xzzzzzz,  \
                             g_0_yyyyyyy,  \
                             g_0_yyyyyyyy, \
                             g_0_yyyyyyyz, \
                             g_0_yyyyyyz,  \
                             g_0_yyyyyyzz, \
                             g_0_yyyyyzz,  \
                             g_0_yyyyyzzz, \
                             g_0_yyyyzzz,  \
                             g_0_yyyyzzzz, \
                             g_0_yyyzzzz,  \
                             g_0_yyyzzzzz, \
                             g_0_yyzzzzz,  \
                             g_0_yyzzzzzz, \
                             g_0_yzzzzzz,  \
                             g_0_yzzzzzzz, \
                             g_0_zzzzzzz,  \
                             g_y_xxxxxxx,  \
                             g_y_xxxxxxy,  \
                             g_y_xxxxxxz,  \
                             g_y_xxxxxyy,  \
                             g_y_xxxxxyz,  \
                             g_y_xxxxxzz,  \
                             g_y_xxxxyyy,  \
                             g_y_xxxxyyz,  \
                             g_y_xxxxyzz,  \
                             g_y_xxxxzzz,  \
                             g_y_xxxyyyy,  \
                             g_y_xxxyyyz,  \
                             g_y_xxxyyzz,  \
                             g_y_xxxyzzz,  \
                             g_y_xxxzzzz,  \
                             g_y_xxyyyyy,  \
                             g_y_xxyyyyz,  \
                             g_y_xxyyyzz,  \
                             g_y_xxyyzzz,  \
                             g_y_xxyzzzz,  \
                             g_y_xxzzzzz,  \
                             g_y_xyyyyyy,  \
                             g_y_xyyyyyz,  \
                             g_y_xyyyyzz,  \
                             g_y_xyyyzzz,  \
                             g_y_xyyzzzz,  \
                             g_y_xyzzzzz,  \
                             g_y_xzzzzzz,  \
                             g_y_yyyyyyy,  \
                             g_y_yyyyyyz,  \
                             g_y_yyyyyzz,  \
                             g_y_yyyyzzz,  \
                             g_y_yyyzzzz,  \
                             g_y_yyzzzzz,  \
                             g_y_yzzzzzz,  \
                             g_y_zzzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xxxxxxx[k] = -g_0_xxxxxxx[k] * ab_y + g_0_xxxxxxxy[k];

                g_y_xxxxxxy[k] = -g_0_xxxxxxy[k] * ab_y + g_0_xxxxxxyy[k];

                g_y_xxxxxxz[k] = -g_0_xxxxxxz[k] * ab_y + g_0_xxxxxxyz[k];

                g_y_xxxxxyy[k] = -g_0_xxxxxyy[k] * ab_y + g_0_xxxxxyyy[k];

                g_y_xxxxxyz[k] = -g_0_xxxxxyz[k] * ab_y + g_0_xxxxxyyz[k];

                g_y_xxxxxzz[k] = -g_0_xxxxxzz[k] * ab_y + g_0_xxxxxyzz[k];

                g_y_xxxxyyy[k] = -g_0_xxxxyyy[k] * ab_y + g_0_xxxxyyyy[k];

                g_y_xxxxyyz[k] = -g_0_xxxxyyz[k] * ab_y + g_0_xxxxyyyz[k];

                g_y_xxxxyzz[k] = -g_0_xxxxyzz[k] * ab_y + g_0_xxxxyyzz[k];

                g_y_xxxxzzz[k] = -g_0_xxxxzzz[k] * ab_y + g_0_xxxxyzzz[k];

                g_y_xxxyyyy[k] = -g_0_xxxyyyy[k] * ab_y + g_0_xxxyyyyy[k];

                g_y_xxxyyyz[k] = -g_0_xxxyyyz[k] * ab_y + g_0_xxxyyyyz[k];

                g_y_xxxyyzz[k] = -g_0_xxxyyzz[k] * ab_y + g_0_xxxyyyzz[k];

                g_y_xxxyzzz[k] = -g_0_xxxyzzz[k] * ab_y + g_0_xxxyyzzz[k];

                g_y_xxxzzzz[k] = -g_0_xxxzzzz[k] * ab_y + g_0_xxxyzzzz[k];

                g_y_xxyyyyy[k] = -g_0_xxyyyyy[k] * ab_y + g_0_xxyyyyyy[k];

                g_y_xxyyyyz[k] = -g_0_xxyyyyz[k] * ab_y + g_0_xxyyyyyz[k];

                g_y_xxyyyzz[k] = -g_0_xxyyyzz[k] * ab_y + g_0_xxyyyyzz[k];

                g_y_xxyyzzz[k] = -g_0_xxyyzzz[k] * ab_y + g_0_xxyyyzzz[k];

                g_y_xxyzzzz[k] = -g_0_xxyzzzz[k] * ab_y + g_0_xxyyzzzz[k];

                g_y_xxzzzzz[k] = -g_0_xxzzzzz[k] * ab_y + g_0_xxyzzzzz[k];

                g_y_xyyyyyy[k] = -g_0_xyyyyyy[k] * ab_y + g_0_xyyyyyyy[k];

                g_y_xyyyyyz[k] = -g_0_xyyyyyz[k] * ab_y + g_0_xyyyyyyz[k];

                g_y_xyyyyzz[k] = -g_0_xyyyyzz[k] * ab_y + g_0_xyyyyyzz[k];

                g_y_xyyyzzz[k] = -g_0_xyyyzzz[k] * ab_y + g_0_xyyyyzzz[k];

                g_y_xyyzzzz[k] = -g_0_xyyzzzz[k] * ab_y + g_0_xyyyzzzz[k];

                g_y_xyzzzzz[k] = -g_0_xyzzzzz[k] * ab_y + g_0_xyyzzzzz[k];

                g_y_xzzzzzz[k] = -g_0_xzzzzzz[k] * ab_y + g_0_xyzzzzzz[k];

                g_y_yyyyyyy[k] = -g_0_yyyyyyy[k] * ab_y + g_0_yyyyyyyy[k];

                g_y_yyyyyyz[k] = -g_0_yyyyyyz[k] * ab_y + g_0_yyyyyyyz[k];

                g_y_yyyyyzz[k] = -g_0_yyyyyzz[k] * ab_y + g_0_yyyyyyzz[k];

                g_y_yyyyzzz[k] = -g_0_yyyyzzz[k] * ab_y + g_0_yyyyyzzz[k];

                g_y_yyyzzzz[k] = -g_0_yyyzzzz[k] * ab_y + g_0_yyyyzzzz[k];

                g_y_yyzzzzz[k] = -g_0_yyzzzzz[k] * ab_y + g_0_yyyzzzzz[k];

                g_y_yzzzzzz[k] = -g_0_yzzzzzz[k] * ab_y + g_0_yyzzzzzz[k];

                g_y_zzzzzzz[k] = -g_0_zzzzzzz[k] * ab_y + g_0_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

            auto g_z_xxxxxxx = cbuffer.data(pk_off + 72 * ccomps * dcomps);

            auto g_z_xxxxxxy = cbuffer.data(pk_off + 73 * ccomps * dcomps);

            auto g_z_xxxxxxz = cbuffer.data(pk_off + 74 * ccomps * dcomps);

            auto g_z_xxxxxyy = cbuffer.data(pk_off + 75 * ccomps * dcomps);

            auto g_z_xxxxxyz = cbuffer.data(pk_off + 76 * ccomps * dcomps);

            auto g_z_xxxxxzz = cbuffer.data(pk_off + 77 * ccomps * dcomps);

            auto g_z_xxxxyyy = cbuffer.data(pk_off + 78 * ccomps * dcomps);

            auto g_z_xxxxyyz = cbuffer.data(pk_off + 79 * ccomps * dcomps);

            auto g_z_xxxxyzz = cbuffer.data(pk_off + 80 * ccomps * dcomps);

            auto g_z_xxxxzzz = cbuffer.data(pk_off + 81 * ccomps * dcomps);

            auto g_z_xxxyyyy = cbuffer.data(pk_off + 82 * ccomps * dcomps);

            auto g_z_xxxyyyz = cbuffer.data(pk_off + 83 * ccomps * dcomps);

            auto g_z_xxxyyzz = cbuffer.data(pk_off + 84 * ccomps * dcomps);

            auto g_z_xxxyzzz = cbuffer.data(pk_off + 85 * ccomps * dcomps);

            auto g_z_xxxzzzz = cbuffer.data(pk_off + 86 * ccomps * dcomps);

            auto g_z_xxyyyyy = cbuffer.data(pk_off + 87 * ccomps * dcomps);

            auto g_z_xxyyyyz = cbuffer.data(pk_off + 88 * ccomps * dcomps);

            auto g_z_xxyyyzz = cbuffer.data(pk_off + 89 * ccomps * dcomps);

            auto g_z_xxyyzzz = cbuffer.data(pk_off + 90 * ccomps * dcomps);

            auto g_z_xxyzzzz = cbuffer.data(pk_off + 91 * ccomps * dcomps);

            auto g_z_xxzzzzz = cbuffer.data(pk_off + 92 * ccomps * dcomps);

            auto g_z_xyyyyyy = cbuffer.data(pk_off + 93 * ccomps * dcomps);

            auto g_z_xyyyyyz = cbuffer.data(pk_off + 94 * ccomps * dcomps);

            auto g_z_xyyyyzz = cbuffer.data(pk_off + 95 * ccomps * dcomps);

            auto g_z_xyyyzzz = cbuffer.data(pk_off + 96 * ccomps * dcomps);

            auto g_z_xyyzzzz = cbuffer.data(pk_off + 97 * ccomps * dcomps);

            auto g_z_xyzzzzz = cbuffer.data(pk_off + 98 * ccomps * dcomps);

            auto g_z_xzzzzzz = cbuffer.data(pk_off + 99 * ccomps * dcomps);

            auto g_z_yyyyyyy = cbuffer.data(pk_off + 100 * ccomps * dcomps);

            auto g_z_yyyyyyz = cbuffer.data(pk_off + 101 * ccomps * dcomps);

            auto g_z_yyyyyzz = cbuffer.data(pk_off + 102 * ccomps * dcomps);

            auto g_z_yyyyzzz = cbuffer.data(pk_off + 103 * ccomps * dcomps);

            auto g_z_yyyzzzz = cbuffer.data(pk_off + 104 * ccomps * dcomps);

            auto g_z_yyzzzzz = cbuffer.data(pk_off + 105 * ccomps * dcomps);

            auto g_z_yzzzzzz = cbuffer.data(pk_off + 106 * ccomps * dcomps);

            auto g_z_zzzzzzz = cbuffer.data(pk_off + 107 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxxx,      \
                             g_0_xxxxxxxz, \
                             g_0_xxxxxxy,  \
                             g_0_xxxxxxyz, \
                             g_0_xxxxxxz,  \
                             g_0_xxxxxxzz, \
                             g_0_xxxxxyy,  \
                             g_0_xxxxxyyz, \
                             g_0_xxxxxyz,  \
                             g_0_xxxxxyzz, \
                             g_0_xxxxxzz,  \
                             g_0_xxxxxzzz, \
                             g_0_xxxxyyy,  \
                             g_0_xxxxyyyz, \
                             g_0_xxxxyyz,  \
                             g_0_xxxxyyzz, \
                             g_0_xxxxyzz,  \
                             g_0_xxxxyzzz, \
                             g_0_xxxxzzz,  \
                             g_0_xxxxzzzz, \
                             g_0_xxxyyyy,  \
                             g_0_xxxyyyyz, \
                             g_0_xxxyyyz,  \
                             g_0_xxxyyyzz, \
                             g_0_xxxyyzz,  \
                             g_0_xxxyyzzz, \
                             g_0_xxxyzzz,  \
                             g_0_xxxyzzzz, \
                             g_0_xxxzzzz,  \
                             g_0_xxxzzzzz, \
                             g_0_xxyyyyy,  \
                             g_0_xxyyyyyz, \
                             g_0_xxyyyyz,  \
                             g_0_xxyyyyzz, \
                             g_0_xxyyyzz,  \
                             g_0_xxyyyzzz, \
                             g_0_xxyyzzz,  \
                             g_0_xxyyzzzz, \
                             g_0_xxyzzzz,  \
                             g_0_xxyzzzzz, \
                             g_0_xxzzzzz,  \
                             g_0_xxzzzzzz, \
                             g_0_xyyyyyy,  \
                             g_0_xyyyyyyz, \
                             g_0_xyyyyyz,  \
                             g_0_xyyyyyzz, \
                             g_0_xyyyyzz,  \
                             g_0_xyyyyzzz, \
                             g_0_xyyyzzz,  \
                             g_0_xyyyzzzz, \
                             g_0_xyyzzzz,  \
                             g_0_xyyzzzzz, \
                             g_0_xyzzzzz,  \
                             g_0_xyzzzzzz, \
                             g_0_xzzzzzz,  \
                             g_0_xzzzzzzz, \
                             g_0_yyyyyyy,  \
                             g_0_yyyyyyyz, \
                             g_0_yyyyyyz,  \
                             g_0_yyyyyyzz, \
                             g_0_yyyyyzz,  \
                             g_0_yyyyyzzz, \
                             g_0_yyyyzzz,  \
                             g_0_yyyyzzzz, \
                             g_0_yyyzzzz,  \
                             g_0_yyyzzzzz, \
                             g_0_yyzzzzz,  \
                             g_0_yyzzzzzz, \
                             g_0_yzzzzzz,  \
                             g_0_yzzzzzzz, \
                             g_0_zzzzzzz,  \
                             g_0_zzzzzzzz, \
                             g_z_xxxxxxx,  \
                             g_z_xxxxxxy,  \
                             g_z_xxxxxxz,  \
                             g_z_xxxxxyy,  \
                             g_z_xxxxxyz,  \
                             g_z_xxxxxzz,  \
                             g_z_xxxxyyy,  \
                             g_z_xxxxyyz,  \
                             g_z_xxxxyzz,  \
                             g_z_xxxxzzz,  \
                             g_z_xxxyyyy,  \
                             g_z_xxxyyyz,  \
                             g_z_xxxyyzz,  \
                             g_z_xxxyzzz,  \
                             g_z_xxxzzzz,  \
                             g_z_xxyyyyy,  \
                             g_z_xxyyyyz,  \
                             g_z_xxyyyzz,  \
                             g_z_xxyyzzz,  \
                             g_z_xxyzzzz,  \
                             g_z_xxzzzzz,  \
                             g_z_xyyyyyy,  \
                             g_z_xyyyyyz,  \
                             g_z_xyyyyzz,  \
                             g_z_xyyyzzz,  \
                             g_z_xyyzzzz,  \
                             g_z_xyzzzzz,  \
                             g_z_xzzzzzz,  \
                             g_z_yyyyyyy,  \
                             g_z_yyyyyyz,  \
                             g_z_yyyyyzz,  \
                             g_z_yyyyzzz,  \
                             g_z_yyyzzzz,  \
                             g_z_yyzzzzz,  \
                             g_z_yzzzzzz,  \
                             g_z_zzzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xxxxxxx[k] = -g_0_xxxxxxx[k] * ab_z + g_0_xxxxxxxz[k];

                g_z_xxxxxxy[k] = -g_0_xxxxxxy[k] * ab_z + g_0_xxxxxxyz[k];

                g_z_xxxxxxz[k] = -g_0_xxxxxxz[k] * ab_z + g_0_xxxxxxzz[k];

                g_z_xxxxxyy[k] = -g_0_xxxxxyy[k] * ab_z + g_0_xxxxxyyz[k];

                g_z_xxxxxyz[k] = -g_0_xxxxxyz[k] * ab_z + g_0_xxxxxyzz[k];

                g_z_xxxxxzz[k] = -g_0_xxxxxzz[k] * ab_z + g_0_xxxxxzzz[k];

                g_z_xxxxyyy[k] = -g_0_xxxxyyy[k] * ab_z + g_0_xxxxyyyz[k];

                g_z_xxxxyyz[k] = -g_0_xxxxyyz[k] * ab_z + g_0_xxxxyyzz[k];

                g_z_xxxxyzz[k] = -g_0_xxxxyzz[k] * ab_z + g_0_xxxxyzzz[k];

                g_z_xxxxzzz[k] = -g_0_xxxxzzz[k] * ab_z + g_0_xxxxzzzz[k];

                g_z_xxxyyyy[k] = -g_0_xxxyyyy[k] * ab_z + g_0_xxxyyyyz[k];

                g_z_xxxyyyz[k] = -g_0_xxxyyyz[k] * ab_z + g_0_xxxyyyzz[k];

                g_z_xxxyyzz[k] = -g_0_xxxyyzz[k] * ab_z + g_0_xxxyyzzz[k];

                g_z_xxxyzzz[k] = -g_0_xxxyzzz[k] * ab_z + g_0_xxxyzzzz[k];

                g_z_xxxzzzz[k] = -g_0_xxxzzzz[k] * ab_z + g_0_xxxzzzzz[k];

                g_z_xxyyyyy[k] = -g_0_xxyyyyy[k] * ab_z + g_0_xxyyyyyz[k];

                g_z_xxyyyyz[k] = -g_0_xxyyyyz[k] * ab_z + g_0_xxyyyyzz[k];

                g_z_xxyyyzz[k] = -g_0_xxyyyzz[k] * ab_z + g_0_xxyyyzzz[k];

                g_z_xxyyzzz[k] = -g_0_xxyyzzz[k] * ab_z + g_0_xxyyzzzz[k];

                g_z_xxyzzzz[k] = -g_0_xxyzzzz[k] * ab_z + g_0_xxyzzzzz[k];

                g_z_xxzzzzz[k] = -g_0_xxzzzzz[k] * ab_z + g_0_xxzzzzzz[k];

                g_z_xyyyyyy[k] = -g_0_xyyyyyy[k] * ab_z + g_0_xyyyyyyz[k];

                g_z_xyyyyyz[k] = -g_0_xyyyyyz[k] * ab_z + g_0_xyyyyyzz[k];

                g_z_xyyyyzz[k] = -g_0_xyyyyzz[k] * ab_z + g_0_xyyyyzzz[k];

                g_z_xyyyzzz[k] = -g_0_xyyyzzz[k] * ab_z + g_0_xyyyzzzz[k];

                g_z_xyyzzzz[k] = -g_0_xyyzzzz[k] * ab_z + g_0_xyyzzzzz[k];

                g_z_xyzzzzz[k] = -g_0_xyzzzzz[k] * ab_z + g_0_xyzzzzzz[k];

                g_z_xzzzzzz[k] = -g_0_xzzzzzz[k] * ab_z + g_0_xzzzzzzz[k];

                g_z_yyyyyyy[k] = -g_0_yyyyyyy[k] * ab_z + g_0_yyyyyyyz[k];

                g_z_yyyyyyz[k] = -g_0_yyyyyyz[k] * ab_z + g_0_yyyyyyzz[k];

                g_z_yyyyyzz[k] = -g_0_yyyyyzz[k] * ab_z + g_0_yyyyyzzz[k];

                g_z_yyyyzzz[k] = -g_0_yyyyzzz[k] * ab_z + g_0_yyyyzzzz[k];

                g_z_yyyzzzz[k] = -g_0_yyyzzzz[k] * ab_z + g_0_yyyzzzzz[k];

                g_z_yyzzzzz[k] = -g_0_yyzzzzz[k] * ab_z + g_0_yyzzzzzz[k];

                g_z_yzzzzzz[k] = -g_0_yzzzzzz[k] * ab_z + g_0_yzzzzzzz[k];

                g_z_zzzzzzz[k] = -g_0_zzzzzzz[k] * ab_z + g_0_zzzzzzzz[k];
            }
        }
    }
}

}  // namespace erirec
