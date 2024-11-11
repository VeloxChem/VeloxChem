#include "ElectronRepulsionContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pixx(CSimdArray<double>&   cbuffer,
                                     const size_t          idx_pixx,
                                     const size_t          idx_sixx,
                                     const size_t          idx_skxx,
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

            /// set up bra offset for contr_buffer_pixx

            const auto pi_off = idx_pixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_xxxxxx = cbuffer.data(pi_off + 0 * ccomps * dcomps);

            auto g_x_xxxxxy = cbuffer.data(pi_off + 1 * ccomps * dcomps);

            auto g_x_xxxxxz = cbuffer.data(pi_off + 2 * ccomps * dcomps);

            auto g_x_xxxxyy = cbuffer.data(pi_off + 3 * ccomps * dcomps);

            auto g_x_xxxxyz = cbuffer.data(pi_off + 4 * ccomps * dcomps);

            auto g_x_xxxxzz = cbuffer.data(pi_off + 5 * ccomps * dcomps);

            auto g_x_xxxyyy = cbuffer.data(pi_off + 6 * ccomps * dcomps);

            auto g_x_xxxyyz = cbuffer.data(pi_off + 7 * ccomps * dcomps);

            auto g_x_xxxyzz = cbuffer.data(pi_off + 8 * ccomps * dcomps);

            auto g_x_xxxzzz = cbuffer.data(pi_off + 9 * ccomps * dcomps);

            auto g_x_xxyyyy = cbuffer.data(pi_off + 10 * ccomps * dcomps);

            auto g_x_xxyyyz = cbuffer.data(pi_off + 11 * ccomps * dcomps);

            auto g_x_xxyyzz = cbuffer.data(pi_off + 12 * ccomps * dcomps);

            auto g_x_xxyzzz = cbuffer.data(pi_off + 13 * ccomps * dcomps);

            auto g_x_xxzzzz = cbuffer.data(pi_off + 14 * ccomps * dcomps);

            auto g_x_xyyyyy = cbuffer.data(pi_off + 15 * ccomps * dcomps);

            auto g_x_xyyyyz = cbuffer.data(pi_off + 16 * ccomps * dcomps);

            auto g_x_xyyyzz = cbuffer.data(pi_off + 17 * ccomps * dcomps);

            auto g_x_xyyzzz = cbuffer.data(pi_off + 18 * ccomps * dcomps);

            auto g_x_xyzzzz = cbuffer.data(pi_off + 19 * ccomps * dcomps);

            auto g_x_xzzzzz = cbuffer.data(pi_off + 20 * ccomps * dcomps);

            auto g_x_yyyyyy = cbuffer.data(pi_off + 21 * ccomps * dcomps);

            auto g_x_yyyyyz = cbuffer.data(pi_off + 22 * ccomps * dcomps);

            auto g_x_yyyyzz = cbuffer.data(pi_off + 23 * ccomps * dcomps);

            auto g_x_yyyzzz = cbuffer.data(pi_off + 24 * ccomps * dcomps);

            auto g_x_yyzzzz = cbuffer.data(pi_off + 25 * ccomps * dcomps);

            auto g_x_yzzzzz = cbuffer.data(pi_off + 26 * ccomps * dcomps);

            auto g_x_zzzzzz = cbuffer.data(pi_off + 27 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxx,      \
                             g_0_xxxxxxx, \
                             g_0_xxxxxxy, \
                             g_0_xxxxxxz, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyy, \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxxzz, \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyy, \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxxzzz, \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyy, \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxxzzzz, \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyy, \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xxzzzzz, \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyy, \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_xzzzzzz, \
                             g_0_yyyyyy,  \
                             g_0_yyyyyz,  \
                             g_0_yyyyzz,  \
                             g_0_yyyzzz,  \
                             g_0_yyzzzz,  \
                             g_0_yzzzzz,  \
                             g_0_zzzzzz,  \
                             g_x_xxxxxx,  \
                             g_x_xxxxxy,  \
                             g_x_xxxxxz,  \
                             g_x_xxxxyy,  \
                             g_x_xxxxyz,  \
                             g_x_xxxxzz,  \
                             g_x_xxxyyy,  \
                             g_x_xxxyyz,  \
                             g_x_xxxyzz,  \
                             g_x_xxxzzz,  \
                             g_x_xxyyyy,  \
                             g_x_xxyyyz,  \
                             g_x_xxyyzz,  \
                             g_x_xxyzzz,  \
                             g_x_xxzzzz,  \
                             g_x_xyyyyy,  \
                             g_x_xyyyyz,  \
                             g_x_xyyyzz,  \
                             g_x_xyyzzz,  \
                             g_x_xyzzzz,  \
                             g_x_xzzzzz,  \
                             g_x_yyyyyy,  \
                             g_x_yyyyyz,  \
                             g_x_yyyyzz,  \
                             g_x_yyyzzz,  \
                             g_x_yyzzzz,  \
                             g_x_yzzzzz,  \
                             g_x_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_y_xxxxxx = cbuffer.data(pi_off + 28 * ccomps * dcomps);

            auto g_y_xxxxxy = cbuffer.data(pi_off + 29 * ccomps * dcomps);

            auto g_y_xxxxxz = cbuffer.data(pi_off + 30 * ccomps * dcomps);

            auto g_y_xxxxyy = cbuffer.data(pi_off + 31 * ccomps * dcomps);

            auto g_y_xxxxyz = cbuffer.data(pi_off + 32 * ccomps * dcomps);

            auto g_y_xxxxzz = cbuffer.data(pi_off + 33 * ccomps * dcomps);

            auto g_y_xxxyyy = cbuffer.data(pi_off + 34 * ccomps * dcomps);

            auto g_y_xxxyyz = cbuffer.data(pi_off + 35 * ccomps * dcomps);

            auto g_y_xxxyzz = cbuffer.data(pi_off + 36 * ccomps * dcomps);

            auto g_y_xxxzzz = cbuffer.data(pi_off + 37 * ccomps * dcomps);

            auto g_y_xxyyyy = cbuffer.data(pi_off + 38 * ccomps * dcomps);

            auto g_y_xxyyyz = cbuffer.data(pi_off + 39 * ccomps * dcomps);

            auto g_y_xxyyzz = cbuffer.data(pi_off + 40 * ccomps * dcomps);

            auto g_y_xxyzzz = cbuffer.data(pi_off + 41 * ccomps * dcomps);

            auto g_y_xxzzzz = cbuffer.data(pi_off + 42 * ccomps * dcomps);

            auto g_y_xyyyyy = cbuffer.data(pi_off + 43 * ccomps * dcomps);

            auto g_y_xyyyyz = cbuffer.data(pi_off + 44 * ccomps * dcomps);

            auto g_y_xyyyzz = cbuffer.data(pi_off + 45 * ccomps * dcomps);

            auto g_y_xyyzzz = cbuffer.data(pi_off + 46 * ccomps * dcomps);

            auto g_y_xyzzzz = cbuffer.data(pi_off + 47 * ccomps * dcomps);

            auto g_y_xzzzzz = cbuffer.data(pi_off + 48 * ccomps * dcomps);

            auto g_y_yyyyyy = cbuffer.data(pi_off + 49 * ccomps * dcomps);

            auto g_y_yyyyyz = cbuffer.data(pi_off + 50 * ccomps * dcomps);

            auto g_y_yyyyzz = cbuffer.data(pi_off + 51 * ccomps * dcomps);

            auto g_y_yyyzzz = cbuffer.data(pi_off + 52 * ccomps * dcomps);

            auto g_y_yyzzzz = cbuffer.data(pi_off + 53 * ccomps * dcomps);

            auto g_y_yzzzzz = cbuffer.data(pi_off + 54 * ccomps * dcomps);

            auto g_y_zzzzzz = cbuffer.data(pi_off + 55 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxx,      \
                             g_0_xxxxxxy, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyy, \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyy, \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyy, \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyy, \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyy, \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_yyyyyy,  \
                             g_0_yyyyyyy, \
                             g_0_yyyyyyz, \
                             g_0_yyyyyz,  \
                             g_0_yyyyyzz, \
                             g_0_yyyyzz,  \
                             g_0_yyyyzzz, \
                             g_0_yyyzzz,  \
                             g_0_yyyzzzz, \
                             g_0_yyzzzz,  \
                             g_0_yyzzzzz, \
                             g_0_yzzzzz,  \
                             g_0_yzzzzzz, \
                             g_0_zzzzzz,  \
                             g_y_xxxxxx,  \
                             g_y_xxxxxy,  \
                             g_y_xxxxxz,  \
                             g_y_xxxxyy,  \
                             g_y_xxxxyz,  \
                             g_y_xxxxzz,  \
                             g_y_xxxyyy,  \
                             g_y_xxxyyz,  \
                             g_y_xxxyzz,  \
                             g_y_xxxzzz,  \
                             g_y_xxyyyy,  \
                             g_y_xxyyyz,  \
                             g_y_xxyyzz,  \
                             g_y_xxyzzz,  \
                             g_y_xxzzzz,  \
                             g_y_xyyyyy,  \
                             g_y_xyyyyz,  \
                             g_y_xyyyzz,  \
                             g_y_xyyzzz,  \
                             g_y_xyzzzz,  \
                             g_y_xzzzzz,  \
                             g_y_yyyyyy,  \
                             g_y_yyyyyz,  \
                             g_y_yyyyzz,  \
                             g_y_yyyzzz,  \
                             g_y_yyzzzz,  \
                             g_y_yzzzzz,  \
                             g_y_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_z_xxxxxx = cbuffer.data(pi_off + 56 * ccomps * dcomps);

            auto g_z_xxxxxy = cbuffer.data(pi_off + 57 * ccomps * dcomps);

            auto g_z_xxxxxz = cbuffer.data(pi_off + 58 * ccomps * dcomps);

            auto g_z_xxxxyy = cbuffer.data(pi_off + 59 * ccomps * dcomps);

            auto g_z_xxxxyz = cbuffer.data(pi_off + 60 * ccomps * dcomps);

            auto g_z_xxxxzz = cbuffer.data(pi_off + 61 * ccomps * dcomps);

            auto g_z_xxxyyy = cbuffer.data(pi_off + 62 * ccomps * dcomps);

            auto g_z_xxxyyz = cbuffer.data(pi_off + 63 * ccomps * dcomps);

            auto g_z_xxxyzz = cbuffer.data(pi_off + 64 * ccomps * dcomps);

            auto g_z_xxxzzz = cbuffer.data(pi_off + 65 * ccomps * dcomps);

            auto g_z_xxyyyy = cbuffer.data(pi_off + 66 * ccomps * dcomps);

            auto g_z_xxyyyz = cbuffer.data(pi_off + 67 * ccomps * dcomps);

            auto g_z_xxyyzz = cbuffer.data(pi_off + 68 * ccomps * dcomps);

            auto g_z_xxyzzz = cbuffer.data(pi_off + 69 * ccomps * dcomps);

            auto g_z_xxzzzz = cbuffer.data(pi_off + 70 * ccomps * dcomps);

            auto g_z_xyyyyy = cbuffer.data(pi_off + 71 * ccomps * dcomps);

            auto g_z_xyyyyz = cbuffer.data(pi_off + 72 * ccomps * dcomps);

            auto g_z_xyyyzz = cbuffer.data(pi_off + 73 * ccomps * dcomps);

            auto g_z_xyyzzz = cbuffer.data(pi_off + 74 * ccomps * dcomps);

            auto g_z_xyzzzz = cbuffer.data(pi_off + 75 * ccomps * dcomps);

            auto g_z_xzzzzz = cbuffer.data(pi_off + 76 * ccomps * dcomps);

            auto g_z_yyyyyy = cbuffer.data(pi_off + 77 * ccomps * dcomps);

            auto g_z_yyyyyz = cbuffer.data(pi_off + 78 * ccomps * dcomps);

            auto g_z_yyyyzz = cbuffer.data(pi_off + 79 * ccomps * dcomps);

            auto g_z_yyyzzz = cbuffer.data(pi_off + 80 * ccomps * dcomps);

            auto g_z_yyzzzz = cbuffer.data(pi_off + 81 * ccomps * dcomps);

            auto g_z_yzzzzz = cbuffer.data(pi_off + 82 * ccomps * dcomps);

            auto g_z_zzzzzz = cbuffer.data(pi_off + 83 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxxxx,      \
                             g_0_xxxxxxz, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxxzz, \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxxzzz, \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxxzzzz, \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xxzzzzz, \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_xzzzzzz, \
                             g_0_yyyyyy,  \
                             g_0_yyyyyyz, \
                             g_0_yyyyyz,  \
                             g_0_yyyyyzz, \
                             g_0_yyyyzz,  \
                             g_0_yyyyzzz, \
                             g_0_yyyzzz,  \
                             g_0_yyyzzzz, \
                             g_0_yyzzzz,  \
                             g_0_yyzzzzz, \
                             g_0_yzzzzz,  \
                             g_0_yzzzzzz, \
                             g_0_zzzzzz,  \
                             g_0_zzzzzzz, \
                             g_z_xxxxxx,  \
                             g_z_xxxxxy,  \
                             g_z_xxxxxz,  \
                             g_z_xxxxyy,  \
                             g_z_xxxxyz,  \
                             g_z_xxxxzz,  \
                             g_z_xxxyyy,  \
                             g_z_xxxyyz,  \
                             g_z_xxxyzz,  \
                             g_z_xxxzzz,  \
                             g_z_xxyyyy,  \
                             g_z_xxyyyz,  \
                             g_z_xxyyzz,  \
                             g_z_xxyzzz,  \
                             g_z_xxzzzz,  \
                             g_z_xyyyyy,  \
                             g_z_xyyyyz,  \
                             g_z_xyyyzz,  \
                             g_z_xyyzzz,  \
                             g_z_xyzzzz,  \
                             g_z_xzzzzz,  \
                             g_z_yyyyyy,  \
                             g_z_yyyyyz,  \
                             g_z_yyyyzz,  \
                             g_z_yyyzzz,  \
                             g_z_yyzzzz,  \
                             g_z_yzzzzz,  \
                             g_z_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
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

}  // namespace erirec
